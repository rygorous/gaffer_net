// Simple byte-aligned binary arithmetic coder (Ilya Muravyov's variant) - public domain - Fabian 'ryg' Giesen 2015
//
// Written for clarity not speed!

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
static int const kProbBits = 12;
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
        for (int i = 0; i < 4; ++i)
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
    size_t read_pos;

    // noncopyable
    BinArithDecoder(BinArithDecoder const &);
    BinArithDecoder &operator =(BinArithDecoder const &);

public:
    // Start decoding
    explicit BinArithDecoder(ByteVec const &source)
        : lo(0), hi(~0u), bytes(source), read_pos(0)
    {
        code = 0;
        for (int i = 0; i < 4; ++i)
            code = (code << 8) | bytes[read_pos++];
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
            code = (code << 8) | bytes[read_pos++];
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

// ---- Data format

static const int kNumCubes = 901;
static const int kNumFrames = 2593;

static const int kRefDist = 6; // distance to reference frame

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

struct Frame
{
    CubeState cubes[kNumCubes];
};

struct ModelSet
{
    typedef BinShiftModel<5> DefaultBit;
    typedef BitTreeModel<DefaultBit, 8> DefaultByte;

    BitTreeModel<DefaultBit, 2> orientation_largest;
    BitTreeModel<DefaultBit, 9> orientation_val;
    DefaultByte pos_xy[3]; // first, second, third byte
    DefaultByte pos_z[2]; // first, second byte
    DefaultBit interacting[2]; // [ref.interacting]
};

typedef BitTreeModel<BinShiftModel<5>, 8> ByteModel;

static void encode_frame(ByteVec &dest, Frame *cur, Frame const *ref)
{
    BinArithEncoder coder(dest);
    ModelSet m;

    for (int i = 0; i < kNumCubes; ++i)
    {
        CubeState *cube = &cur->cubes[i];
        CubeState const *refc = &ref->cubes[i];

        m.orientation_largest.encode(coder, cube->orientation_largest);
        for (int j = 0; j < 3; ++j)
            m.orientation_val.encode(coder, (&cube->orientation_a)[j]);
        for (int j = 0; j < 2; ++j)
        {
            int p = (&cube->position_x)[j];
            m.pos_xy[0].encode(coder, (p >>  0) & 0xff);
            m.pos_xy[1].encode(coder, (p >>  8) & 0xff);
            m.pos_xy[2].encode(coder, (p >> 16) & 0xff);
        }
        {
            int p = cube->position_z;
            m.pos_z[0].encode(coder, (p >> 0) & 0xff);
            m.pos_z[1].encode(coder, (p >> 8) & 0xff);
        }

        m.interacting[refc->interacting].encode(coder, cube->interacting);
    }
}

static void decode_frame(ByteVec const &src, Frame *cur, Frame const *ref)
{
    BinArithDecoder coder(src);
    ModelSet m;

    for (int i = 0; i < kNumCubes; ++i)
    {
        CubeState *cube = &cur->cubes[i];
        CubeState const *refc = &ref->cubes[i];

        cube->orientation_largest = m.orientation_largest.decode(coder);
        for (int j = 0; j < 3; ++j)
            (&cube->orientation_a)[j] = m.orientation_val.decode(coder);

        for (int j = 0; j < 2; ++j)
        {
            int p = 0;
            p |= m.pos_xy[0].decode(coder) <<  0;
            p |= m.pos_xy[1].decode(coder) <<  8;
            p |= m.pos_xy[2].decode(coder) << 16;
            // sign extend
            (&cube->position_x)[j] = (int) (((int32_t) (p << 8)) >> 8);
        }
        {
            int p = 0;
            p |= m.pos_z[0].decode(coder) << 0;
            p |= m.pos_z[1].decode(coder) << 8;
            cube->position_z = p;
        }

        cube->interacting = m.interacting[refc->interacting].decode(coder);
    }
}

// ---- I/O and main

static void read_data(char const *filename, Frame *frames, int num_frames)
{
    CubeState dummy;

    FILE *f = fopen("delta_data.bin", "rb");
    if (!f)
    {
        printf("data missing!\n");
        exit(1);
    }

    for (int frame = 0; frame < num_frames; ++frame)
    {
        for (int cube = 0; cube < kNumCubes; ++cube)
        {
            if (fread(&frames[frame].cubes[cube], sizeof(CubeState), 1, f) != 1)
            {
                printf("error reading frame %d cube %d!\n", frame, cube);
                exit(1);
            }
            fread(&dummy, sizeof(CubeState), 1, f); // skip ref CubeState
        }
    }
    fclose(f);
}

static void write_data(char const *filename, Frame *frames, int num_frames)
{
    CubeState null_state = {};

    FILE *f = fopen("output.bin", "wb");
    if (!f)
    {
        printf("error writing output!\n");
        exit(1);
    }

    for (int frame = 0; frame < num_frames; ++frame)
    {
        Frame *cur_frame = &frames[frame];
        Frame *ref_frame = (frame >= kRefDist) ? &frames[frame - kRefDist] : 0;

        for (int cube = 0; cube < kNumCubes; ++cube)
        {
            // current data
            fwrite(&cur_frame->cubes[cube], sizeof(CubeState), 1, f);
            // ref data from ref frame (0 when no ref frame)
            fwrite(ref_frame ? &ref_frame->cubes[cube] : &null_state, sizeof(CubeState), 1, f);
        }
    }

    fclose(f);
}

int main()
{
    Frame null_frame;
    memset(&null_frame, 0, sizeof(null_frame));

    // Read the data
    Frame *frames = new Frame[kNumFrames];
    printf("reading...\n");
    read_data("delta_data.bin", frames, kNumFrames);
    printf("done.\n");

    // Coding loop
    ByteVec packet_buf;

    size_t packet_size_sum = 0;
    Frame out;

    for (int frame = 0; frame < kNumFrames; ++frame)
    {
        Frame *cur = &frames[frame];
        Frame *ref = (frame >= kRefDist) ? &frames[frame - kRefDist] : &null_frame;

        packet_buf.clear();
        encode_frame(packet_buf, cur, ref);
        decode_frame(packet_buf, &out, ref);

        if (memcmp(out.cubes, cur->cubes, sizeof(out.cubes)) != 0)
        {
            printf("decode mismatch on frame %d\n", frame);
            return 1;
        }

        packet_size_sum += packet_buf.size();
    }

    printf("total packed size %d\n", (int)packet_size_sum);
    printf("%.2f bytes/frame\n", (double)packet_size_sum / (double)kNumFrames);

    // Write output
    write_data("output.bin", frames, kNumFrames);

    // Clean up
    delete[] frames;
    return 0;
}

// vim:et:sts=4:sw=4
