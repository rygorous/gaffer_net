#define CHECK_RESULTS

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
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

template<int Inertia0, int Inertia1>
struct TwoBinShiftModel
{
    uint16_t p0, p1;

    TwoBinShiftModel() : p0(kProbMax / 4), p1(kProbMax / 4) {}

    void encode(BinArithEncoder &enc, int bit)
    {
        enc.encode(bit, p0 + p1);
        adapt(bit);
    }

    int decode(BinArithDecoder &dec)
    {
        int bit = dec.decode(p0 + p1);
        adapt(bit);
        return bit;
    }

    void adapt(int bit)
    {
        // Note prob never his 0 or kProbMax with this update rule!
        if (bit)
        {
            p0 += (kProbMax/2 - p0) >> Inertia0;
            p1 += (kProbMax/2 - p1) >> Inertia1;
        }
        else
        {
            p0 -= p0 >> Inertia0;
            p1 -= p1 >> Inertia1;
        }
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
    int orientation[3];
    int position[3];
    int interacting;
};

// Prediction state. Not sent in the stream; inferred from the data
// sent to aid coding.
struct PredState
{
    int changing;
    int orient_delta[3];
    int vel[3];
};

struct ModelSet
{
    static const int kNumPosCtx = 10;

    typedef TwoBinShiftModel<3, 7> DefaultBit;
    typedef SExpGolombModel<DefaultBit> SExpGolomb;

    DefaultBit orientation_different[2]; // [refp.changing]
    BitTreeModel<DefaultBit, 2> orientation_largest[4*4]; // [orient_context]
    SExpGolomb orientation_delta;
    DefaultBit orientation_signflip[2]; // [second_largest_sign]
    SExpGolomb orientation_val;

    DefaultBit pos_different[2]; // [orientation_differs]
    SExpGolomb pos_delta[kNumPosCtx]; // [pos_ctx]

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

static int xyzw_from_abc(int abc_ind, int largest)
{
    return abc_ind + (abc_ind >= largest);
}

static int abc_from_xyzw(int xyzw_ind, int largest)
{
    assert(xyzw_ind != largest);
    return xyzw_ind - (xyzw_ind >= largest);
}

static int orient_context(CubeState const *cube)
{
    // Largest axis is elided. Find index and magnitude of second-largest.
    int v[3];
    for (int i = 0; i < 3; ++i)
        v[i] = abs(cube->orientation[i] - 256);

    int abc_ind;
    if (v[0] >= v[1])
        abc_ind = (v[0] >= v[2]) ? 0 : 2;
    else
        abc_ind = (v[1] >= v[2]) ? 1 : 2;

    int ctx = cube->orientation_largest;
    if (v[abc_ind] >= 128) // second-largest axis is getting closer to cross-over
        ctx += 4 * (abc_ind + 1);

    return ctx;
}

static int pos_context(int dv)
{
    int v = abs(dv);
    int ctx = 0;
    while (v > 1 && ctx < ModelSet::kNumPosCtx - 1)
    {
        ++ctx;
        v /= 2;
    }
    return ctx;
}

static void unpack_quat_prediction(int dest[4], int const src[3], int largest)
{
    for (int i = 0; i < 3; ++i)
        dest[xyzw_from_abc(i, largest)] = src[i];
    dest[largest] = 450;
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

        int diff_orient = (cube->orientation_largest != refc->orientation_largest), diff_pos = 0;
        for (int i = 0; i < 3; ++i)
        {
            pred->orient_delta[i] = cube->orientation[i] - refc->orientation[i];
            pred->vel[i] = cube->position[i] - refc->position[i];
            diff_orient |= pred->orient_delta[i];
            diff_pos |= pred->vel[i];
        }

        m.orientation_different[refp->changing].encode(coder, diff_orient);
        if (diff_orient)
        {
            int orient_ctx = orient_context(refc);
            m.orientation_largest[orient_ctx].encode(coder, cube->orientation_largest);
            if (cube->orientation_largest == refc->orientation_largest)
            {
                for (int i = 0; i < 3; ++i)
                    m.orientation_delta.encode(coder, pred->orient_delta[i]);
            }
            else
            {
                int old_largest = refc->orientation_largest;
                int new_largest = cube->orientation_largest;
                int old[4];
                unpack_quat_prediction(old, refc->orientation, old_largest);

                int sign_context = old[new_largest] < 256;
                if (cube->orientation[abc_from_xyzw(old_largest, new_largest)] < 256)
                {
                    m.orientation_signflip[sign_context].encode(coder, 1);
                    for (int i = 0; i < 4; ++i)
                        old[i] = 512 - old[i];
                }
                else
                    m.orientation_signflip[sign_context].encode(coder, 0);

                for (int i = 0; i < 3; ++i)
                    m.orientation_delta.encode(coder, cube->orientation[i] - old[xyzw_from_abc(i, new_largest)]);
            }
        }

        m.pos_different[diff_orient != 0].encode(coder, diff_pos != 0);
        if (diff_pos)
        {
            for (int i = 0; i < 3; ++i)
            {
                int ctx = pos_context(refp->vel[i]);
                m.pos_delta[ctx].encode(coder, pred->vel[i] - refp->vel[i]);
            }
        }

        m.interacting[refc->interacting].encode(coder, cube->interacting);

        // NOTE: in general, we would need to account for variable frame
        // spacing here. But in this testbed we always predict from 6 frames
        // ago, so no problem.
        pred->changing = (diff_orient | diff_pos) != 0;
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
            int orient_ctx = orient_context(refc);
            cube->orientation_largest = (int) m.orientation_largest[orient_ctx].decode(coder);
            if (cube->orientation_largest == refc->orientation_largest)
            {
                for (int i = 0; i < 3; ++i)
                    cube->orientation[i] = refc->orientation[i] + m.orientation_delta.decode(coder);
            }
            else
            {
                int old_largest = refc->orientation_largest;
                int new_largest = cube->orientation_largest;
                int old[4];
                unpack_quat_prediction(old, refc->orientation, old_largest);

                if (m.orientation_signflip[old[new_largest] < 256].decode(coder))
                {
                    for (int i = 0; i < 4; ++i)
                        old[i] = 512 - old[i];
                }

                for (int i = 0; i < 3; ++i)
                    cube->orientation[i] = m.orientation_delta.decode(coder) + old[xyzw_from_abc(i, new_largest)];
            }
        }
        else
        {
            cube->orientation_largest = refc->orientation_largest;
            for (int i = 0; i < 3; ++i)
                cube->orientation[i] = refc->orientation[i];
        }

        pred->vel[0] = pred->vel[1] = pred->vel[2] = 0;

        if (m.pos_different[diff_orient].decode(coder))
        {
            for (int i = 0; i < 3; ++i)
            {
                int ctx = pos_context(refp->vel[i]);
                pred->vel[i] = refp->vel[i] + m.pos_delta[ctx].decode(coder);
            }
        }

        for (int i = 0; i < 3; ++i)
            cube->position[i] = refc->position[i] + pred->vel[i];

        cube->interacting = m.interacting[refc->interacting].decode(coder);
        pred->changing = (int(diff_orient) | pred->vel[0] | pred->vel[1] | pred->vel[2]) != 0;
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
    size_t packet_count = 0;
    Frame out;

    clock_t enc_start = clock();

    // Gaffer says skip the first 6 frames. Okay.
    for (int frame = 6; frame < num_frames; ++frame)
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
        ++packet_count;
    }

    double enc_time = double(clock() - enc_start) / CLOCKS_PER_SEC;
    printf("processing took %.2fs (%.2fus/frame)\n", enc_time, 1e6*enc_time / (double)packet_count);
    printf("total packed size %d\n", (int)packet_size_sum);

    double bytes_per_frame = (double)packet_size_sum / (double)packet_count;
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
