#include <iostream>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <random>
#include <iomanip>
#include <omp.h>
#include <ctime>
#include <cstdlib>
#include <limits>
using namespace std;

// ============================================================
//  SKINNY-128 / single-tweakey, byte-oriented state.
//
//  State convention follows the Midori-style code:
//      state[row][col] = s_{4 * col + row}
//
//  The mask arrays below are therefore column-major:
//      pos = 4 * col + row
//
//  This convention is important for the supplied trail.
// ============================================================

uint8_t dot_table[256][256];

inline int cell_row(int pos) { return pos & 3; }
inline int cell_col(int pos) { return pos >> 2; }

const uint8_t ZERO_TK[16] = { 0 };
const bool USE_ROUND_CONSTANTS = false;
const bool USE_ROUND_TWEAKEY = true;

// SKINNY-128 8-bit S-box.
const uint8_t skinny_sbox8[256] = {
    0x65, 0x4c, 0x6a, 0x42 ,0x4b ,0x63 ,0x43 ,0x6b ,0x55 ,0x75 ,0x5a ,0x7a ,0x53 ,0x73 ,0x5b ,0x7b ,
    0x35, 0x8c, 0x3a, 0x81 ,0x89 ,0x33 ,0x80 ,0x3b ,0x95 ,0x25 ,0x98 ,0x2a ,0x90 ,0x23 ,0x99 ,0x2b ,
    0xe5, 0xcc, 0xe8, 0xc1 ,0xc9 ,0xe0 ,0xc0 ,0xe9 ,0xd5 ,0xf5 ,0xd8 ,0xf8 ,0xd0 ,0xf0 ,0xd9 ,0xf9 ,
    0xa5, 0x1c, 0xa8, 0x12 ,0x1b ,0xa0 ,0x13 ,0xa9 ,0x05 ,0xb5 ,0x0a ,0xb8 ,0x03 ,0xb0 ,0x0b ,0xb9 ,
    0x32, 0x88, 0x3c, 0x85 ,0x8d ,0x34 ,0x84 ,0x3d ,0x91 ,0x22 ,0x9c ,0x2c ,0x94 ,0x24 ,0x9d ,0x2d ,
    0x62, 0x4a, 0x6c, 0x45 ,0x4d ,0x64 ,0x44 ,0x6d ,0x52 ,0x72 ,0x5c ,0x7c ,0x54 ,0x74 ,0x5d ,0x7d ,
    0xa1, 0x1a, 0xac, 0x15 ,0x1d ,0xa4 ,0x14 ,0xad ,0x02 ,0xb1 ,0x0c ,0xbc ,0x04 ,0xb4 ,0x0d ,0xbd ,
    0xe1, 0xc8, 0xec, 0xc5 ,0xcd ,0xe4 ,0xc4 ,0xed ,0xd1 ,0xf1 ,0xdc ,0xfc ,0xd4 ,0xf4 ,0xdd ,0xfd ,
    0x36, 0x8e, 0x38, 0x82 ,0x8b ,0x30 ,0x83 ,0x39 ,0x96 ,0x26 ,0x9a ,0x28 ,0x93 ,0x20 ,0x9b ,0x29 ,
    0x66, 0x4e, 0x68, 0x41 ,0x49 ,0x60 ,0x40 ,0x69 ,0x56 ,0x76 ,0x58 ,0x78 ,0x50 ,0x70 ,0x59 ,0x79 ,
    0xa6, 0x1e, 0xaa, 0x11 ,0x19 ,0xa3 ,0x10 ,0xab ,0x06 ,0xb6 ,0x08 ,0xba ,0x00 ,0xb3 ,0x09 ,0xbb ,
    0xe6, 0xce, 0xea, 0xc2 ,0xcb ,0xe3 ,0xc3 ,0xeb ,0xd6 ,0xf6 ,0xda ,0xfa ,0xd3 ,0xf3 ,0xdb ,0xfb ,
    0x31, 0x8a, 0x3e, 0x86 ,0x8f ,0x37 ,0x87 ,0x3f ,0x92 ,0x21 ,0x9e ,0x2e ,0x97 ,0x27 ,0x9f ,0x2f ,
    0x61, 0x48, 0x6e, 0x46 ,0x4f ,0x67 ,0x47 ,0x6f ,0x51 ,0x71 ,0x5e ,0x7e ,0x57 ,0x77 ,0x5f ,0x7f ,
    0xa2, 0x18, 0xae, 0x16 ,0x1f ,0xa7 ,0x17 ,0xaf ,0x01 ,0xb2 ,0x0e ,0xbe ,0x07 ,0xb7 ,0x0f ,0xbf ,
    0xe2, 0xca, 0xee, 0xc6 ,0xcf ,0xe7 ,0xc7 ,0xef ,0xd2 ,0xf2 ,0xde ,0xfe ,0xd7 ,0xf7 ,0xdf ,0xff
};




// First 62 round constants.
const uint8_t skinny_rc[62] = {
    0x01,0x03,0x07,0x0F,0x1F,0x3E,0x3D,0x3B,0x37,0x2F,0x1E,0x3C,0x39,0x33,0x27,0x0E,
    0x1D,0x3A,0x35,0x2B,0x16,0x2C,0x18,0x30,0x21,0x02,0x05,0x0B,0x17,0x2E,0x1C,0x38,
    0x31,0x23,0x06,0x0D,0x1B,0x36,0x2D,0x1A,0x34,0x29,0x12,0x24,0x08,0x11,0x22,0x04,
    0x09,0x13,0x26,0x0C,0x19,0x32,0x25,0x0A,0x15,0x2A,0x14,0x28,0x10,0x20
};

void SubCells(uint8_t state[4][4])
{
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        state[r][c] = skinny_sbox8[state[r][c]];
    }
}

void AddConstants(uint8_t state[4][4], int round)
{
    if (!USE_ROUND_CONSTANTS) return;

    uint8_t rc = skinny_rc[round];

    // Constants are added to cells s0, s4, s8, i.e. row 0/1/2, column 0.
    state[0][0] ^= (rc & 0x0F);
    state[1][0] ^= ((rc >> 4) & 0x03);
    state[2][0] ^= 0x02;
}

void AddRoundTweakey(uint8_t state[4][4], const uint8_t tk[16])
{
    if (!USE_ROUND_TWEAKEY) return;

    // TK is added to the first two rows.
    // Here tk uses the same column-major convention as the state.
    for (int c = 0; c < 4; c++) {
        for (int r = 0; r < 2; r++) {
            int pos = 4 * c + r;
            state[r][c] ^= tk[pos];
        }
    }
}

void UpdateTweakey(uint8_t tk[16])
{
    // TK1 permutation. With ZERO_TK this has no visible effect.
    const int perm[16] = {
        9,15,8,13,
        10,14,12,11,
        0,1,2,3,
        4,5,6,7
    };

    uint8_t old[16];
    memcpy(old, tk, 16);

    for (int i = 0; i < 16; i++)
        tk[i] = old[perm[i]];
}

void ShiftRows(uint8_t state[4][4])
{
    uint8_t old[4][4];
    memcpy(old, state, 16);

    // Row r is rotated right by r positions.
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            state[r][c] = old[r][(c - r + 4) & 3];
        }
    }
}

void MixColumns(uint8_t state[4][4])
{
    // SKINNY MixColumns matrix:
    // [1 0 1 1]
    // [1 0 0 0]
    // [0 1 1 0]
    // [1 0 1 0]
    for (int c = 0; c < 4; c++) {
        uint8_t x0 = state[0][c];
        uint8_t x1 = state[1][c];
        uint8_t x2 = state[2][c];
        uint8_t x3 = state[3][c];

        state[0][c] = x0 ^ x2 ^ x3;
        state[1][c] = x0;
        state[2][c] = x1 ^ x2;
        state[3][c] = x0 ^ x2;
    }
}

void SkinnyRound(uint8_t state[4][4], int round, uint8_t tk[16])
{
    SubCells(state);
    AddConstants(state, round);
    AddRoundTweakey(state, tk);
    ShiftRows(state);
    MixColumns(state);
    UpdateTweakey(tk);
}

uint8_t dot_state(const uint8_t mask[16], const uint8_t state[4][4])
{
    uint8_t v = 0;

    for (int pos = 0; pos < 16; pos++)
        v ^= dot_table[mask[pos]][state[cell_row(pos)][cell_col(pos)]];

    return v & 1;
}



// ============================================================
//  Given 6-round SKINNY conditional linear trail.
// ============================================================

/*
The S - box output mask for round 0:
0x90, 0xD0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x90, 0x90, 0x00, 0x90, 0x00, 0x00, 0x00, 0xD0, 0x00
The S - box input mask for rounds 1
0x90, 0x00, 0x00, 0x00, 0x00, 0x00, 0xD0, 0x00, 0x00, 0x90, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
The S - box output mask for rounds 1
0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
The S - box input mask for rounds 2
0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00
The S - box output mask for rounds 2
0x00, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00
The S - box input mask for rounds 3
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
The S - box output mask for rounds 3
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0B, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
The S - box input mask for rounds 4
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0B, 0x00, 0x0B
The S - box output mask for rounds 4
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00, 0x28
The S - box input mask for rounds 5
0x00, 0x28, 0x28, 0x28, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00, 0x00, 0x28, 0x00, 0x00, 0x00, 0x00
*/


const uint8_t r1_sbox_output_mask[16] = {
0x90, 0xD0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x90, 0x90, 0x00, 0x90, 0x00, 0x00, 0x00, 0xD0, 0x00
};


const uint8_t r6_sbox_input_mask[16] = {
0x00, 0x28, 0x28, 0x28, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00, 0x00, 0x28, 0x00, 0x00, 0x00, 0x00
};




// Boundary masks and conditions are derived automatically from one-condition CLAT.
uint8_t in_mask[16];
uint8_t in_cond_mask[16];
uint8_t in_cond_value[16];

uint8_t out_mask[16];
uint8_t out_cond_mask[16];
uint8_t out_cond_value[16];

struct CLATChoice {
    uint8_t lin_mask = 0;
    uint8_t cond_mask = 0;
    uint8_t cond_value = 0;
    int numerator = 0;
    double weight = numeric_limits<double>::infinity();
};

// Input-side CLAT:
// fixed output mask lambda_o;
// search input mask lambda_i and condition <T, x> = t.
CLATChoice search_input_side_choice(uint8_t fixed_output_mask)
{
    CLATChoice best;

    if (fixed_output_mask == 0)
    {
        best.weight = 0.0;
        return best;
    }

    for (int T = 1; T < 256; T++) {
        for (int tv = 0; tv <= 1; tv++) {
            for (int li = 0; li < 256; li++) {
                int num = 0;
                int den = 0;

                for (int x = 0; x < 256; x++)
                {
                    if (dot_table[T][x] != tv)
                        continue;

                    den++;

                    uint8_t e = dot_table[li][x] ^ dot_table[fixed_output_mask][skinny_sbox8[x]];

                    num += (e == 0) ? 1 : -1;
                }

                if (num == 0)
                    continue;

                double corr = fabs((double)num / (double)den);
                double w = -2.0 * log2(corr) + 1.0;

                if (w + 1e-12 < best.weight) {
                    best.lin_mask = (uint8_t)li;
                    best.cond_mask = (uint8_t)T;
                    best.cond_value = (uint8_t)tv;
                    best.numerator = num;
                    best.weight = w;

                    if (fabs(w - 1.0) < 1e-12)
                        return best;
                }
            }
        }
    }

    return best;
}

// Output-side CLAT:
// fixed input mask lambda_i;
// search output mask lambda_o and condition <T, S(x)> = t.
CLATChoice search_output_side_choice(uint8_t fixed_input_mask)
{
    CLATChoice best;

    if (fixed_input_mask == 0) {
        best.weight = 0.0;
        return best;
    }

    for (int T = 1; T < 256; T++) {
        for (int tv = 0; tv <= 1; tv++) {
            for (int lo = 0; lo < 256; lo++) {
                int num = 0;
                int den = 0;

                for (int x = 0; x < 256; x++) {
                    uint8_t y = skinny_sbox8[x];

                    if (dot_table[T][y] != tv)
                        continue;

                    den++;

                    uint8_t e =
                        dot_table[fixed_input_mask][x]
                        ^ dot_table[lo][y];

                    num += (e == 0) ? 1 : -1;
                }

                if (num == 0)
                    continue;

                double corr = fabs((double)num / (double)den);
                double w = -2.0 * log2(corr) + 1.0;

                if (w + 1e-12 < best.weight) {
                    best.lin_mask = (uint8_t)lo;
                    best.cond_mask = (uint8_t)T;
                    best.cond_value = (uint8_t)tv;
                    best.numerator = num;
                    best.weight = w;

                    if (fabs(w - 1.0) < 1e-12)
                        return best;
                }
            }
        }
    }

    return best;
}

void init_boundary_choices(bool verbose)
{
    memset(in_mask, 0, sizeof(in_mask));
    memset(in_cond_mask, 0, sizeof(in_cond_mask));
    memset(in_cond_value, 0, sizeof(in_cond_value));

    memset(out_mask, 0, sizeof(out_mask));
    memset(out_cond_mask, 0, sizeof(out_cond_mask));
    memset(out_cond_value, 0, sizeof(out_cond_value));

    double total_boundary_weight = 0.0;
    int condition_count = 0;

    if (verbose)
        cout << "Input-side CLAT choices for round 1:\n";

    for (int pos = 0; pos < 16; pos++) {
        CLATChoice ch = search_input_side_choice(r1_sbox_output_mask[pos]);

        in_mask[pos] = ch.lin_mask;
        in_cond_mask[pos] = ch.cond_mask;
        in_cond_value[pos] = ch.cond_value;

        if (r1_sbox_output_mask[pos] != 0) {
            condition_count++;
            total_boundary_weight += ch.weight;


            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << ": fixed_out=0x" << hex << setw(2) << setfill('0')
                    << (int)r1_sbox_output_mask[pos]
                    << "  in=0x" << setw(2) << (int)ch.lin_mask
                    << "  cond=<0x" << setw(2) << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << setfill(' ') << "\n";
            }
        }
    }

    if (verbose)
        cout << "Output-side CLAT choices for round 6:\n";

    for (int pos = 0; pos < 16; pos++) {
        CLATChoice ch = search_output_side_choice(r6_sbox_input_mask[pos]);

        out_mask[pos] = ch.lin_mask;
        out_cond_mask[pos] = ch.cond_mask;
        out_cond_value[pos] = ch.cond_value;

        if (r6_sbox_input_mask[pos] != 0) {
            condition_count++;
            total_boundary_weight += ch.weight;


            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << ": fixed_in=0x" << hex << setw(2) << setfill('0')
                    << (int)r6_sbox_input_mask[pos]
                    << "  out=0x" << setw(2) << (int)ch.lin_mask
                    << "  cond=<0x" << setw(2) << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << setfill(' ') << "\n";
            }
        }
    }
}

bool check_conditions(const uint8_t state[4][4],
    const uint8_t cond_mask[16],
    const uint8_t cond_value[16])
{
    for (int pos = 0; pos < 16; pos++) {
        if (cond_mask[pos] == 0)
            continue;

        uint8_t x = state[cell_row(pos)][cell_col(pos)];
        uint8_t lhs = dot_table[cond_mask[pos]][x];
        uint8_t rhs = cond_value[pos] & 1;

        if (lhs != rhs)
            return false;
    }

    return true;
}

void random_state(uint8_t state[4][4],
    std::uniform_int_distribution<int>& dist,
    std::mt19937_64& gen)
{
    for (int pos = 0; pos < 16; pos++)
        state[cell_row(pos)][cell_col(pos)] = static_cast<uint8_t>(dist(gen) & 0xFF);
}



void test_Skinny6_conditional()
{
    // Total generated plaintexts = TEST_NUM1 * TEST_NUM2.
    // There are 11 one-bit conditions, so surviving samples are about N / 2^11.
    const int TEST_NUM1 = 2048;
    const int TEST_NUM2 = 1 << 20;

    cout << "log2(TEST_NUM) = "
        << log2((double)TEST_NUM1) + log2((double)TEST_NUM2)
        << endl;

    long long count = 0;
    long long total = 0;

    uint8_t Tweakey[16];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    cout << "Tweakey: ";
    for (int pos = 0; pos < 16; pos++)
    {
        Tweakey[pos] = static_cast<uint8_t>(dist(gen));

        std::cout << "0x" << std::hex << std::uppercase << std::setw(2) << std::setfill('0') << static_cast<int>(Tweakey[pos]) << " ";
    }

    cout << endl << std::dec;


#pragma omp parallel for
    for (int test1 = 0; test1 < TEST_NUM1; test1++) {
        uint64_t seed =
            std::random_device{}()
            ^ ((uint64_t)time(nullptr) << 1)
            ^ ((uint64_t)omp_get_thread_num() << 32)
            ^ (uint64_t)test1;

        std::mt19937_64 local_rng(seed);
        std::uniform_int_distribution<int> local_dist(0, 255);

        long long count_temp = 0;
        long long total_temp = 0;

        uint8_t state[4][4];




        for (int test2 = 0; test2 < TEST_NUM2; test2++) {

            random_state(state, local_dist, local_rng);

            if (!check_conditions(state, in_cond_mask, in_cond_value))
                continue;



            uint8_t result = dot_state(in_mask, state);

            uint8_t tk[16];

            for (int pos = 0; pos < 16; pos++)
                tk[pos] = Tweakey[pos];

            // Rounds 1 to 5 are full rounds.
            for (int r = 0; r < 5; r++)
                SkinnyRound(state, r, tk);

            // Stop after SubCells of round 6.
            // Output-side conditions are imposed on round-6 S-box outputs.
            SubCells(state);

            if (!check_conditions(state, out_cond_mask, out_cond_value))
                continue;

            total_temp++;

            result ^= dot_state(out_mask, state);

            if (result == 0)
                count_temp++;
            else
                count_temp--;
        }

#pragma omp critical
        {
            count += count_temp;
            total += total_temp;
        }
    }

    cout << "log2(Actual_NUM) = " << log2(total) << endl;

    if (total == 0 || count == 0) {
        cout << "Correlation cannot be computed: total="
            << total << ", count=" << count << endl;
        return;
    }

    double correlation_log2 =
        log2((double)llabs(count)) - log2((double)total);

    cout << "Correlation: 2^("
        << fixed << setprecision(6) << correlation_log2 << ")"
        << "    count: " << count << endl;
}



int main()
{
    for (int m = 0; m < 256; m++) {
        for (int v = 0; v < 256; v++) {
            uint8_t temp = (uint8_t)(m & v);
            uint8_t t = 0;

            for (int i = 0; i < 8; i++)
                t ^= (temp >> i) & 1;
            dot_table[m][v] = t;
        }
    }
    init_boundary_choices(true);

    for (int i = 0; i < 128; i++)
    {
        time_t t1, t2;
        time(&t1);
        test_Skinny6_conditional();
        time(&t2);

        cout << "Time measured: " << t2 - t1 << " seconds." << endl;
        cout << "--------------------------------------" << endl;
    }

    return 0;
}
