#define CHECK_RESULTS

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <vector>

// Probabilities are expressed in fixed point, with kProbBits bits of
// resolution. No need to go overboard with this.
static int const kProbBits = 15;
static uint32_t const kProbMax = 1u << kProbBits;

// Type used for buffers.
typedef std::vector<uint8_t> ByteVec;

// Binary arithmetic encoder (Ilya Muravyov's variant)
// Encodes/decodes a string of binary (0/1) events with
// probabilities that are not 1/2.
//
// This code is written for clarity, not performance.
class BinArithEncoder
{
    uint32_t lo, hi;
    ByteVec &bytes;

    // noncopyable
    BinArithEncoder(BinArithEncoder const &);
    BinArithEncoder &operator =(BinArithEncoder const &);

public:
    // Initialize
    explicit BinArithEncoder(ByteVec &target) : lo(0), hi(~0u), bytes(target) { }

    // Finish encoding - flushes remaining codeword
    ~BinArithEncoder()
    {
        // Find shortest encoding that still decodes to the right symbols.
        // The decoder implicitly zero-pads w
        uint32_t round_up = 0xffffffu;
        while (round_up)
        {
            if ((lo | round_up) != ~0u)
            {
                uint32_t rounded = (lo + round_up) & ~round_up;
                if (rounded <= hi) // inside interval, we're good!
                {
                    lo = rounded;
                    break;
                }
            }

            round_up >>= 8;
        }

        while (lo)
        {
            bytes.push_back(lo >> 24);
            lo <<= 8;
        }
    }

    // Encode a binary symbol "bit" with the probability of a 1 being "prob".
    // Note that prob=0 (or prob=1<<kProbBits) really mean that a 1 (or 0,
    // respectively) cannot occur!
    void encode(int bit, uint32_t prob)
    {
        // Midpoint of active probability interval subdivided via prob
        uint32_t x = lo + ((uint64_t(hi - lo) * prob) >> kProbBits);

        if (bit)
            hi = x;
        else
            lo = x + 1;

        // Renormalize: when top byte of lo/hi is same, shift it out.
        while ((lo ^ hi) < (1u << 24))
        {
            bytes.push_back(lo >> 24);
            lo <<= 8;
            hi = (hi << 8) | 0xff;
        }
    }
};

// Corresponding decoder.
class BinArithDecoder
{
    uint32_t code, lo, hi;
    ByteVec const &bytes;
    size_t read_pos, size;

    // noncopyable
    BinArithDecoder(BinArithDecoder const &);
    BinArithDecoder &operator =(BinArithDecoder const &);

    uint8_t getb()
    {
        if (read_pos < size)
            return bytes[read_pos++];
        else
            return 0;
    }

public:
    // Start decoding
    explicit BinArithDecoder(ByteVec const &source)
        : lo(0), hi(~0u), bytes(source), read_pos(0)
    {
        code = 0;
        size = source.size();
        for (int i = 0; i < 4; ++i)
            code = (code << 8) | getb();
    }

    // Decode a binary symbol with the probability of a 1 being "prob".
    int decode(uint32_t prob)
    {
        int bit;

        // Midpoint of active probability interval subdivided via prob
        uint32_t x = lo + ((uint64_t(hi - lo) * prob) >> kProbBits);

        if (code <= x)
        {
            hi = x;
            bit = 1;
        }
        else
        {
            lo = x + 1;
            bit = 0;
        }

        // Renormalize
        while ((lo ^ hi) < (1u << 24))
        {
            code = (code << 8) | getb();
            lo <<= 8;
            hi = (hi << 8) | 0xff;
        }

        return bit;
    }
};

// ---- A few basic models

// NOTE: Again, this is written for clarity and ease of tinkering.
// In practice, you will write more direct code for these once you've
// figured out your coding structure.

// Adaptive binary model. These are pretty good!
// Lower Inertia = faster.
//
// You typically build more sophisticated models out of these
// by having lots of them and choosing the active model based on
// context.
template<int Inertia>
struct BinShiftModel
{
    uint16_t prob;

    BinShiftModel() : prob(kProbMax / 2) {}

    void encode(BinArithEncoder &enc, int bit)
    {
        enc.encode(bit, prob);
        adapt(bit);
    }

    int decode(BinArithDecoder &dec)
    {
        int bit = dec.decode(prob);
        adapt(bit);
        return bit;
    }

    void adapt(int bit)
    {
        // Note prob never his 0 or kProbMax with this update rule!
        if (bit)
            prob += (kProbMax - prob) >> Inertia;
        else
            prob -= prob >> Inertia;
    }
};

// BitTree model. A tree-shaped cascade of BinShiftModels.
// This is the de-facto standard way to build a multi-symbol coder
// (values with NumBits bits) out of binary models.
//
// LZMA (as in 7zip/xz) uses this type of model (backed by a BinShiftModel
// as above) for its literals.
template<typename BitModel, int NumBits>
struct BitTreeModel
{
    static size_t const kNumSyms = 1 << NumBits;
    static size_t const kMSB = kNumSyms / 2;

    BitModel model[kNumSyms - 1];

    void encode(BinArithEncoder &enc, size_t value)
    {
        assert(value < kNumSyms);

        // The first bit sent is the MSB of the value and coded without context
        // Second bit is the bit below the MSB, using the value of the MSB as context
        // and so forth.
        //
        // 1 + 2 + 4 + ... = 2^NumBits - 1 contexts.
        // Numbering the MSB context 1 and then shifting in the coded bits from the
        // bottom is a convenient way to index them. (So ctx is 1-based)
        size_t ctx = 1;
        while (ctx < kNumSyms)
        {
            int bit = (value & kMSB) != 0;
            value += value; // shift value by 1 for next iter
            model[ctx - 1].encode(enc, bit);
            ctx += ctx + bit; // shift in "bit" into context
        }
    }

    size_t decode(BinArithDecoder &dec)
    {
        // Corresponding decoder is nice and easy:
        size_t ctx = 1;
        while (ctx < kNumSyms)
            ctx += ctx + model[ctx - 1].decode(dec);

        return ctx - kNumSyms;
    }
};

// Unsigned exponential Golomb-style model.
template<typename MagModel>
struct UExpGolombModel
{
    BitTreeModel<MagModel, 5> mag;

    void encode(BinArithEncoder &enc, uint32_t value)
    {
        ++value; // we code non-negative values

        // determine magnitude (position of highest 1 bit)
        // and send it in unary.
        // bitscan is the better way to do this.
        uint32_t m = 0;
        while (value >= (2u << m))
            ++m;
        mag.encode(enc, m);

        // send remaining bits flat, MSB->LSB
        uint32_t mask = m ? 1u << (m - 1) : 0;
        while (mask)
        {
            uint32_t bit = (value & mask) != 0;
            enc.encode(bit, kProbMax / 2);
            mask >>= 1;
        }
    }

    uint32_t decode(BinArithDecoder &dec)
    {
        // decode magnitude code
        uint32_t m = (uint32_t) mag.decode(dec);

        // decode value bits
        uint32_t v = 1;
        for (uint32_t i = 0; i < m; ++i)
            v += v + dec.decode(kProbMax / 2);

        return v - 1;
    }
};

// Signed exponential Golomb-style model.
template<typename MagModel>
struct SExpGolombModel
{
    UExpGolombModel<MagModel> abs_coder;

    void encode(BinArithEncoder &enc, int32_t value)
    {
        uint32_t absv = (value < 0) ? -value : value;
        abs_coder.encode(enc, absv);
        if (absv)
            enc.encode(value < 0, kProbMax / 2);
    }

    int32_t decode(BinArithDecoder &dec)
    {
        int32_t v = abs_coder.decode(dec);
        if (v)
        {
            if (dec.decode(kProbMax / 2))
                v = -v;
        }

        return v;
    }
};

// ---- Data format

static const int kNumCubes = 901;

static const int kRefDist = 6; // distance to reference frame
static const int kFrameRate = 60; // just used to calc kbps

struct CubeState
{
    int orientation_largest;
    int orientation_a;
    int orientation_b;
    int orientation_c;
    int position_x;
    int position_y;
    int position_z;
    int interacting;
};

// Prediction state. Not sent in the stream; inferred from the data
// sent to aid coding.
struct PredState
{
    int changing;
    int vel_x;
    int vel_y;
    int vel_z;
};

struct ModelSet
{
    typedef BinShiftModel<5> DefaultBit;
    typedef SExpGolombModel<DefaultBit> SExpGolomb;

    DefaultBit orientation_different[2]; // [refp.changing]
    BitTreeModel<DefaultBit, 2> orientation_largest[4]; // [ref.orientation_largest]
    SExpGolomb orientation_delta;
    SExpGolomb orientation_val;

    DefaultBit pos_different[2]; // [orientation_differs]
    SExpGolomb pos_xy;
    SExpGolomb pos_z;

    DefaultBit interacting[2]; // [ref.interacting]
};

struct Frame
{
    Frame();

    CubeState cubes[kNumCubes];
    ModelSet models; // coding state
    PredState pred[kNumCubes]; // prediction state
};

Frame::Frame()
{
    memset(cubes, 0, sizeof(cubes));
    memset(pred, 0, sizeof(pred));
}

static void encode_frame(ByteVec &dest, Frame *cur, Frame const *ref)
{
    BinArithEncoder coder(dest);
    ModelSet &m = cur->models;

    // Start with ref frame models
    m = ref->models;

    for (int i = 0; i < kNumCubes; ++i)
    {
        CubeState *cube = &cur->cubes[i];
        PredState *pred = &cur->pred[i];
        CubeState const *refc = &ref->cubes[i];
        PredState const *refp = &ref->pred[i];

        bool diff_orient = false;

        if (cube->orientation_largest != refc->orientation_largest ||
            cube->orientation_a != refc->orientation_a ||
            cube->orientation_b != refc->orientation_b ||
            cube->orientation_c != refc->orientation_c)
        {
            diff_orient = true;
            m.orientation_different[refp->changing].encode(coder, 1);
            m.orientation_largest[refc->orientation_largest].encode(coder, cube->orientation_largest);
            if (cube->orientation_largest == refc->orientation_largest)
            {
                int da = cube->orientation_a - refc->orientation_a;
                int db = cube->orientation_b - refc->orientation_b;
                int dc = cube->orientation_c - refc->orientation_c;

                m.orientation_delta.encode(coder, da);
                m.orientation_delta.encode(coder, db);
                m.orientation_delta.encode(coder, dc);
            }
            else
            {
                m.orientation_val.encode(coder, cube->orientation_a - 256);
                m.orientation_val.encode(coder, cube->orientation_b - 256);
                m.orientation_val.encode(coder, cube->orientation_c - 256);
            }
        }
        else
            m.orientation_different[refp->changing].encode(coder, 0);

        int dx = cube->position_x - refc->position_x;
        int dy = cube->position_y - refc->position_y;
        int dz = cube->position_z - refc->position_z;
        if (dx || dy || dz)
        {
            m.pos_different[diff_orient].encode(coder, 1);
            m.pos_xy.encode(coder, dx - refp->vel_x);
            m.pos_xy.encode(coder, dy - refp->vel_y);
            m.pos_z.encode(coder, dz - refp->vel_z);
        }
        else
            m.pos_different[diff_orient].encode(coder, 0);

        m.interacting[refc->interacting].encode(coder, cube->interacting);

        // NOTE: in general, we would need to account for variable frame
        // spacing here. But in this testbed we always predict from 6 frames
        // ago.
        pred->vel_x = dx;
        pred->vel_y = dy;
        pred->vel_z = dz;
        pred->changing = (int(diff_orient) | dx | dy | dz) != 0;
    }
}

static void decode_frame(ByteVec const &src, Frame *cur, Frame const *ref)
{
    BinArithDecoder coder(src);
    ModelSet &m = cur->models;

    // Start with ref frame models
    m = ref->models;

    for (int i = 0; i < kNumCubes; ++i)
    {
        CubeState *cube = &cur->cubes[i];
        PredState *pred = &cur->pred[i];
        CubeState const *refc = &ref->cubes[i];
        PredState const *refp = &ref->pred[i];
        bool diff_orient = false;

        if (m.orientation_different[refp->changing].decode(coder))
        {
            diff_orient = true;
            cube->orientation_largest = (int) m.orientation_largest[refc->orientation_largest].decode(coder);
            if (cube->orientation_largest == refc->orientation_largest)
            {
                cube->orientation_a = refc->orientation_a + m.orientation_delta.decode(coder);
                cube->orientation_b = refc->orientation_b + m.orientation_delta.decode(coder);
                cube->orientation_c = refc->orientation_c + m.orientation_delta.decode(coder);
            }
            else
            {
                cube->orientation_a = m.orientation_val.decode(coder) + 256;
                cube->orientation_b = m.orientation_val.decode(coder) + 256;
                cube->orientation_c = m.orientation_val.decode(coder) + 256;
            }
        }
        else
        {
            cube->orientation_largest = refc->orientation_largest;
            cube->orientation_a = refc->orientation_a;
            cube->orientation_b = refc->orientation_b;
            cube->orientation_c = refc->orientation_c;
        }

        int dx = 0, dy = 0, dz = 0;
        if (m.pos_different[diff_orient].decode(coder))
        {
            dx = refp->vel_x + m.pos_xy.decode(coder);
            dy = refp->vel_y + m.pos_xy.decode(coder);
            dz = refp->vel_z + m.pos_z.decode(coder);
        }

        cube->position_x = refc->position_x + dx;
        cube->position_y = refc->position_y + dy;
        cube->position_z = refc->position_z + dz;
        cube->interacting = m.interacting[refc->interacting].decode(coder);

        pred->vel_x = dx;
        pred->vel_y = dy;
        pred->vel_z = dz;
        pred->changing = (int(diff_orient) | dx | dy | dz) != 0;
    }
}

// ---- I/O and main

static Frame *read_data(char const *filename, int &num_frames, Frame *initial)
{
    FILE *f = fopen(filename, "rb");
    if (!f)
    {
        printf("data missing!\n");
        exit(1);
    }

    fseek(f, 0, SEEK_END);
    num_frames = ftell(f) / (kNumCubes * sizeof(CubeState)) - 1;
    fseek(f, 0, SEEK_SET);

    Frame *frames = new Frame[num_frames];

    // read initial frame cubes
    if (fread(initial->cubes, sizeof(CubeState), kNumCubes, f) != kNumCubes)
    {
        printf("error reading initial frame!\n");
        exit(1);
    }

    for (int frame = 0; frame < num_frames; ++frame)
    {
        if (fread(frames[frame].cubes, sizeof(CubeState), kNumCubes, f) != kNumCubes)
        {
            printf("error reading frame %d!\n", frame);
            exit(1);
        }
    }

    fclose(f);
    return frames;
}

static void write_data(char const *filename, Frame *frames, int num_frames, Frame const *initial)
{
    FILE *f = fopen(filename, "wb");
    if (!f)
    {
        printf("error writing output!\n");
        exit(1);
    }

    fwrite(initial->cubes, sizeof(CubeState), kNumCubes, f);

    for (int frame = 0; frame < num_frames; ++frame)
        fwrite(frames[frame].cubes, sizeof(CubeState), kNumCubes, f);

    fclose(f);
}

int main()
{
    Frame initial_frame;

    // Read the data
    printf("reading...\n");
    int num_frames;
    Frame *frames = read_data("delta_data_realnew.bin", num_frames, &initial_frame);
    printf("done.\n");

    // Coding loop
    ByteVec packet_buf;

    size_t packet_size_sum = 0;
    Frame out;

    for (int frame = 0; frame < num_frames; ++frame)
    {
        Frame *cur = &frames[frame];
        Frame *ref = (frame >= kRefDist) ? &frames[frame - kRefDist] : &initial_frame;

        packet_buf.clear();
        encode_frame(packet_buf, cur, ref);
        decode_frame(packet_buf, &out, ref);

#ifdef CHECK_RESULTS
        if (memcmp(out.cubes, cur->cubes, sizeof(out.cubes)) != 0)
        {
            printf("decode mismatch on frame %d\n", frame);
            return 1;
        }
#endif

        packet_size_sum += packet_buf.size();
    }

    printf("total packed size %d\n", (int)packet_size_sum);

    double bytes_per_frame = (double)packet_size_sum / (double)num_frames;
    double kbps = bytes_per_frame * kFrameRate * 8.0 / 1000.0;

    printf("%.2f bytes/frame\n", bytes_per_frame);
    printf("%.2f kbps\n", kbps);

    // Write output
    write_data("output.bin", frames, num_frames, &initial_frame);

    // Clean up
    delete[] frames;
    return 0;
}

// vim:et:sts=4:sw=4
