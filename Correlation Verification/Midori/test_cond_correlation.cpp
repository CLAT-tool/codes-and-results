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

thread_local std::mt19937_64 rng(std::random_device{}());

// ============================================================
//  Midori-128 basic components, byte-oriented state.
//
//  State convention used in this file follows the Midori paper:
//      (s0, s1, ..., s15) is stored column by column, i.e.,
//      state[row][col] = s_{4 * col + row}, row,col in {0,1,2,3}.
//
//  Therefore, mask/condition arrays indexed by t=0..15 refer to Midori cells s_t.
//  If your search program outputs masks in row-major order, transpose the masks
//  before filling the arrays below, or change cell_row/cell_col accordingly.
// ============================================================


uint8_t dot_table[256][256];


inline int cell_row(int pos) { return pos & 3; }
inline int cell_col(int pos) { return pos >> 2; }

// Midori-128 uses four 8-bit S-boxes SSb0..SSb3 generated from Sb1.
const uint8_t midori_sb1[16] = {
    0x1, 0x0, 0x5, 0x3,
    0xE, 0x2, 0xF, 0x7,
    0xD, 0xA, 0x9, 0xB,
    0xC, 0x8, 0x4, 0x6
};

inline void sb1_bits(uint8_t a, uint8_t out[4])
{
    uint8_t y = midori_sb1[a & 0xF];
    out[0] = (y >> 3) & 1;
    out[1] = (y >> 2) & 1;
    out[2] = (y >> 1) & 1;
    out[3] = y & 1;
}

uint8_t SSbi(uint8_t x, int i)
{
    // x0 is the most significant bit and x7 is the least significant bit.
    uint8_t x0 = (x >> 7) & 1;
    uint8_t x1 = (x >> 6) & 1;
    uint8_t x2 = (x >> 5) & 1;
    uint8_t x3 = (x >> 4) & 1;
    uint8_t x4 = (x >> 3) & 1;
    uint8_t x5 = (x >> 2) & 1;
    uint8_t x6 = (x >> 1) & 1;
    uint8_t x7 = x & 1;

    uint8_t n0[4], n1[4], b[8];

    if (i == 0) {
        sb1_bits((x4 << 3) | (x1 << 2) | (x6 << 1) | x3, n0);
        sb1_bits((x0 << 3) | (x5 << 2) | (x2 << 1) | x7, n1);
        uint8_t tmp[8] = { n1[0], n0[1], n1[2], n0[3], n0[0], n1[1], n0[2], n1[3] };
        memcpy(b, tmp, 8);
    }
    else if (i == 1) {
        sb1_bits((x1 << 3) | (x6 << 2) | (x7 << 1) | x0, n0);
        sb1_bits((x5 << 3) | (x2 << 2) | (x3 << 1) | x4, n1);
        uint8_t tmp[8] = { n0[3], n0[0], n1[1], n1[2], n1[3], n1[0], n0[1], n0[2] };
        memcpy(b, tmp, 8);
    }
    else if (i == 2) {
        sb1_bits((x2 << 3) | (x3 << 2) | (x4 << 1) | x1, n0);
        sb1_bits((x6 << 3) | (x7 << 2) | (x0 << 1) | x5, n1);
        uint8_t tmp[8] = { n1[2], n0[3], n0[0], n0[1], n0[2], n1[3], n1[0], n1[1] };
        memcpy(b, tmp, 8);
    }
    else {
        sb1_bits((x7 << 3) | (x4 << 2) | (x1 << 1) | x2, n0);
        sb1_bits((x3 << 3) | (x0 << 2) | (x5 << 1) | x6, n1);
        uint8_t tmp[8] = { n1[1], n0[2], n0[3], n1[0], n0[1], n1[2], n1[3], n0[0] };
        memcpy(b, tmp, 8);
    }

    uint8_t y = 0;
    for (int j = 0; j < 8; j++)
        y = (uint8_t)((y << 1) | (b[j] & 1));
    return y;
}

void SubCell(uint8_t state[4][4])
{
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        state[r][c] = SSbi(state[r][c], pos & 3);
    }
}

void ShuffleCell(uint8_t state[4][4])
{
    // New (s0,...,s15) <= old (s0,s10,s5,s15,s14,s4,s11,s1,s9,s3,s12,s6,s7,s13,s2,s8).
    const int perm[16] = { 0, 10, 5, 15, 14, 4, 11, 1, 9, 3, 12, 6, 7, 13, 2, 8 };
    uint8_t old[16];
    for (int pos = 0; pos < 16; pos++)
        old[pos] = state[cell_row(pos)][cell_col(pos)];
    for (int pos = 0; pos < 16; pos++)
        state[cell_row(pos)][cell_col(pos)] = old[perm[pos]];
}

void MixColumn(uint8_t state[4][4])
{
    // Binary involutory matrix with zero diagonal and one elsewhere.
    // For each column: y_i = XOR of the other three bytes.
    for (int c = 0; c < 4; c++) {
        uint8_t a0 = state[0][c];
        uint8_t a1 = state[1][c];
        uint8_t a2 = state[2][c];
        uint8_t a3 = state[3][c];
        state[0][c] = a1 ^ a2 ^ a3;
        state[1][c] = a0 ^ a2 ^ a3;
        state[2][c] = a0 ^ a1 ^ a3;
        state[3][c] = a0 ^ a1 ^ a2;
    }
}

// First 19 Midori-128 round constants alpha_i. Each entry is XORed to the LSB of one byte.
const uint8_t midori_alpha[19][16] = {
    {0,0,0,1, 0,1,0,1, 1,0,1,1, 0,0,1,1},
    {0,1,1,1, 1,0,0,0, 1,1,0,0, 0,0,0,0},
    {1,0,1,0, 0,1,0,0, 0,0,1,1, 0,1,0,1},
    {0,1,1,0, 0,0,1,0, 0,0,0,1, 0,0,1,1},
    {0,0,0,1, 0,0,0,0, 0,1,0,0, 1,1,1,1},
    {1,1,0,1, 0,0,1,1, 0,0,0,1, 0,1,0,0},
    {0,0,0,0, 0,0,1,0, 0,1,1,0, 0,0,1,0},
    {0,0,0,0, 1,0,1,1, 1,1,0,0, 1,0,0,0},
    {1,0,0,1, 0,1,0,0, 1,0,0,0, 0,0,0,1},
    {0,1,0,0, 0,0,0,0, 1,0,1,1, 1,0,0,0},
    {0,1,1,1, 0,0,0,1, 1,0,0,1, 0,1,1,1},
    {0,0,1,0, 0,0,0,0, 1,1,1,0, 1,0,1,0},
    {0,1,0,1, 0,0,0,1, 0,0,1,1, 0,0,0,0},
    {1,1,1,1, 1,0,0,0, 1,1,1,1, 1,0,0,0},
    {1,1,0,1, 1,1,1,1, 1,0,0,1, 0,0,0,0},
    {0,1,1,1, 1,1,0,0, 1,0,0,0, 0,0,1,1},
    {0,0,0,1, 1,1,0,0, 0,0,1,0, 0,1,0,0},
    {0,0,1,0, 0,0,0,1, 1,0,1,1, 0,1,0,0},
    {0,1,1,0, 0,0,0,0, 1,0,0,0, 1,0,1,0}
};

const bool USE_ROUND_CONSTANTS = true;

// Midori-128 key schedule used here:
//   whitening key WK = K
//   round key RK_i = K xor alpha_i, where alpha_i is injected into the LSB
//   of each byte according to midori_alpha[i][pos].
void KeyAdd(uint8_t state[4][4], const uint8_t key[16])
{
    for (int pos = 0; pos < 16; pos++)
        state[cell_row(pos)][cell_col(pos)] ^= key[pos];
}

void MakeRoundKey(const uint8_t master_key[16], int round, uint8_t round_key[16])
{
    if (round < 0 || round >= 19) {
        cerr << "[key-schedule error] Midori round index out of range: "
            << round << endl;
        std::abort();
    }

    for (int pos = 0; pos < 16; pos++) {
        round_key[pos] = master_key[pos];

        if (USE_ROUND_CONSTANTS)
            round_key[pos] ^= (midori_alpha[round][pos] & 1);
    }
}

void AddRoundKeyAndConstant(uint8_t state[4][4], const uint8_t master_key[16], int round)
{
    uint8_t round_key[16];
    MakeRoundKey(master_key, round, round_key);
    KeyAdd(state, round_key);
}

void MidoriRound(uint8_t state[4][4], int round, const uint8_t master_key[16])
{
    SubCell(state);
    ShuffleCell(state);
    MixColumn(state);
    AddRoundKeyAndConstant(state, master_key, round);
}

void random_midori_key(uint8_t key[16], std::mt19937_64& gen)
{
    std::uniform_int_distribution<int> dist(0, 255);

    for (int pos = 0; pos < 16; pos++)
        key[pos] = static_cast<uint8_t>(dist(gen) & 0xFF);
}

void print_midori_key(const uint8_t key[16])
{
    cout << "Midori-128 master key: ";

    for (int pos = 0; pos < 16; pos++) {
        cout << "0x" << hex << uppercase << setw(2) << setfill('0')
            << (int)key[pos];

        if (pos + 1 != 16)
            cout << " ";
    }

    cout << setfill(' ') << nouppercase << dec << endl;
}

bool check_round_key_add_once(const uint8_t master_key[16], int round)
{
    uint8_t base[4][4];
    uint8_t via_function[4][4];
    uint8_t via_expanded_key[4][4];
    uint8_t round_key[16];

    for (int pos = 0; pos < 16; pos++)
    {
        uint8_t v = static_cast<uint8_t>((0xA5 + 17 * pos + 13 * round) & 0xFF);
        base[cell_row(pos)][cell_col(pos)] = v;
        via_function[cell_row(pos)][cell_col(pos)] = v;
        via_expanded_key[cell_row(pos)][cell_col(pos)] = v;
    }

    MakeRoundKey(master_key, round, round_key);
    AddRoundKeyAndConstant(via_function, master_key, round);
    KeyAdd(via_expanded_key, round_key);

    for (int pos = 0; pos < 16; pos++) {
        uint8_t a = via_function[cell_row(pos)][cell_col(pos)];
        uint8_t b = via_expanded_key[cell_row(pos)][cell_col(pos)];

        if (a != b) {
            cerr << "[round-key-add error] round=" << round
                << " pos=" << pos
                << " AddRoundKeyAndConstant=0x" << hex << (int)a
                << " KeyAdd(expanded)=0x" << (int)b
                << dec << endl;
            return false;
        }
    }

    AddRoundKeyAndConstant(via_function, master_key, round);

    for (int pos = 0; pos < 16; pos++) {
        uint8_t a = via_function[cell_row(pos)][cell_col(pos)];
        uint8_t b = base[cell_row(pos)][cell_col(pos)];

        if (a != b) {
            cerr << "[round-key involution error] round=" << round
                << " pos=" << pos
                << " after double add=0x" << hex << (int)a
                << " original=0x" << (int)b
                << dec << endl;
            return false;
        }
    }

    return true;
}

bool check_midori_key_schedule(const uint8_t master_key[16], bool verbose)
{
    for (int r = 0; r < 19; r++) {
        for (int pos = 0; pos < 16; pos++) {
            if ((midori_alpha[r][pos] & ~1) != 0) {
                cerr << "[key-schedule error] midori_alpha is not binary at round "
                    << r << ", pos " << pos << endl;
                return false;
            }
        }
    }

    for (int r = 0; r < 19; r++) {
        if (!check_round_key_add_once(master_key, r))
            return false;
    }

    for (int r = 0; r < 19; r++) {
        uint8_t round_key[16];
        MakeRoundKey(master_key, r, round_key);

        for (int pos = 0; pos < 16; pos++) {
            uint8_t expected = master_key[pos];

            if (USE_ROUND_CONSTANTS)
                expected ^= (midori_alpha[r][pos] & 1);

            if (round_key[pos] != expected) {
                cerr << "[key-schedule formula error] round=" << r
                    << " pos=" << pos
                    << " got=0x" << hex << (int)round_key[pos]
                    << " expected=0x" << (int)expected
                    << dec << endl;
                return false;
            }
        }
    }

    return true;
}


uint8_t dot_state(const uint8_t mask[16], const uint8_t state[4][4])
{
    uint8_t v = 0;
    for (int pos = 0; pos < 16; pos++)
        v ^= dot_table[mask[pos]][state[cell_row(pos)][cell_col(pos)]];
    return v & 1;
}

// ============================================================
//  Given masks from your 5-round Midori conditional linear trail.
// ============================================================

const uint8_t r1_sbox_output_mask[16] = {
    0x10, 0x00, 0x40, 0x00, 0x00, 0x10, 0x08, 0x00, 0x40, 0x08, 0x10, 0x00, 0x08, 0x40, 0x00, 0x00
};

const uint8_t r2_to_r5_sbox_input_mask[4][16] = {
    { // round 2 S-box input mask
        0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00
    },
    { // round 3 S-box input mask
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00
    },
    { // round 4 S-box input mask
        0x00, 0x00, 0x00, 0x00, 0x10, 0x10, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    },
    { // round 5 S-box input mask
        0x08, 0x08, 0x00, 0x08, 0x40, 0x00, 0x40, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x10, 0x10
    }
};

// ============================================================
//  Boundary masks and conditions are now derived automatically.
//
//  Input side:
//      fixed output mask = r1_sbox_output_mask[pos]
//      search input mask lambda_i and condition <T, x> = t.
//
//  Output side:
//      fixed input mask = r2_to_r5_sbox_input_mask[3][pos]
//      search output mask lambda_o and condition <T, S_i(x)> = t.
//
//  Midori-128 uses four position-dependent 8-bit S-boxes:
//      S_i = SSbi(x, pos & 3)
// ============================================================

uint8_t in_mask[16];
uint8_t in_cond_mask[16];
uint8_t in_cond_value[16];

uint8_t out_mask[16];
uint8_t out_cond_mask[16];
uint8_t out_cond_value[16];

int input_condition_count = 0;
int output_condition_count = 0;

struct CLATChoice {
    uint8_t lin_mask = 0;
    uint8_t cond_mask = 0;
    uint8_t cond_value = 0;
    int numerator = 0;
    double weight = std::numeric_limits<double>::infinity();
};

// Input-side CLAT for one Midori byte:
// fixed output mask lambda_o;
// search input mask lambda_i and one input condition <T, x> = t.
CLATChoice search_input_side_choice(int sbox_index, uint8_t fixed_output_mask)
{
    CLATChoice best;

    if (fixed_output_mask == 0) {
        best.weight = 0.0;
        return best;
    }

    for (int T = 1; T < 256; T++) {
        for (int tv = 0; tv <= 1; tv++) {
            for (int li = 0; li < 256; li++) {
                int num = 0;
                int den = 0;

                for (int x = 0; x < 256; x++) {
                    if (dot_table[T][x] != tv)
                        continue;

                    den++;

                    uint8_t y = SSbi((uint8_t)x, sbox_index);
                    uint8_t e =
                        dot_table[li][x]
                        ^ dot_table[fixed_output_mask][y];

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

                    // For one-bit conditions, weight 1 means |corr| = 1.
                    // It cannot be improved, so we can stop early.
                    if (fabs(w - 1.0) < 1e-12)
                        return best;
                }
            }
        }
    }

    return best;
}

// Output-side CLAT for one Midori byte:
// fixed input mask lambda_i;
// search output mask lambda_o and one output condition <T, S_i(x)> = t.
CLATChoice search_output_side_choice(int sbox_index, uint8_t fixed_input_mask)
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
                    uint8_t y = SSbi((uint8_t)x, sbox_index);

                    if (dot_table[T][y] != tv)
                        continue;

                    den++;

                    uint8_t e = dot_table[fixed_input_mask][x] ^ dot_table[lo][y];

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

int count_active_cells(const uint8_t mask[16])
{
    int cnt = 0;

    for (int i = 0; i < 16; i++) {
        if (mask[i] != 0)
            cnt++;
    }

    return cnt;
}

void init_boundary_choices(bool verbose = true)
{
    memset(in_mask, 0, sizeof(in_mask));
    memset(in_cond_mask, 0, sizeof(in_cond_mask));
    memset(in_cond_value, 0, sizeof(in_cond_value));

    memset(out_mask, 0, sizeof(out_mask));
    memset(out_cond_mask, 0, sizeof(out_cond_mask));
    memset(out_cond_value, 0, sizeof(out_cond_value));

    input_condition_count = 0;
    output_condition_count = 0;

    double boundary_weight = 0.0;

    if (verbose)
        cout << "Input-side CLAT choices for round 1:" << endl;

    for (int pos = 0; pos < 16; pos++) {
        int sbox_index = pos & 3;
        CLATChoice ch =
            search_input_side_choice(sbox_index, r1_sbox_output_mask[pos]);

        in_mask[pos] = ch.lin_mask;
        in_cond_mask[pos] = ch.cond_mask;
        in_cond_value[pos] = ch.cond_value;

        if (r1_sbox_output_mask[pos] != 0) {
            input_condition_count++;
            boundary_weight += ch.weight;

            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << "  S" << sbox_index
                    << ": fixed_out=0x" << hex << uppercase << setw(2) << setfill('0')
                    << (int)r1_sbox_output_mask[pos]
                    << "  in=0x" << setw(2) << (int)ch.lin_mask
                    << "  cond=<0x" << setw(2) << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << setfill(' ') << nouppercase << endl;;
            }
        }
    }

    if (verbose)
        cout << "Output-side CLAT choices for round 5:\n";

    for (int pos = 0; pos < 16; pos++) {
        int sbox_index = pos & 3;
        CLATChoice ch =
            search_output_side_choice(sbox_index, r2_to_r5_sbox_input_mask[3][pos]);

        out_mask[pos] = ch.lin_mask;
        out_cond_mask[pos] = ch.cond_mask;
        out_cond_value[pos] = ch.cond_value;

        if (r2_to_r5_sbox_input_mask[3][pos] != 0) {
            output_condition_count++;
            boundary_weight += ch.weight;

            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << "  S" << sbox_index
                    << ": fixed_in=0x" << hex << uppercase << setw(2) << setfill('0')
                    << (int)r2_to_r5_sbox_input_mask[3][pos]
                    << "  out=0x" << setw(2) << (int)ch.lin_mask
                    << "  cond=<0x" << setw(2) << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << setfill(' ') << nouppercase << endl;
            }
        }
    }
}

bool check_conditions(const uint8_t state[4][4], const uint8_t cond_mask[16], const uint8_t cond_value[16])
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

bool all_zero_mask(const uint8_t mask[16])
{
    for (int i = 0; i < 16; i++)
        if ((mask[i] & 0xFF) != 0)
            return false;
    return true;
}

void random_state(uint8_t state[4][4], std::uniform_int_distribution<int>& dist,
    std::mt19937_64& gen)
{
    for (int pos = 0; pos < 16; pos++)
        state[cell_row(pos)][cell_col(pos)] = static_cast<uint8_t>(dist(gen) & 0xFF);
}

void test_Midori5_conditional()
{
    // Adjust these according to the expected number of surviving samples.
    const int TEST_NUM1 = 4096 * 2;
    const int TEST_NUM2 = 1 << 24;

    cout << "log2(TEST_NUM) = " << log2((double)TEST_NUM1) + log2((double)TEST_NUM2) << endl;

    uint8_t MasterKey[16];
    uint64_t key_seed = std::random_device{}() ^ ((uint64_t)time(nullptr) << 17);
    std::mt19937_64 key_rng(key_seed);
    random_midori_key(MasterKey, key_rng);
    print_midori_key(MasterKey);

    if (!check_midori_key_schedule(MasterKey, true))
    {
        cerr << "Key schedule / round-key-add check failed. Abort this experiment." << endl;
        return;
    }

    long long count = 0;
    long long total = 0;

#pragma omp parallel for
    for (int test1 = 0; test1 < TEST_NUM1; test1++) 
    {
        std::mt19937_64 local_rng(std::random_device{}() ^
            (uint64_t)time(nullptr) ^
            ((uint64_t)omp_get_thread_num() << 32) ^
            (uint64_t)test1);
        std::uniform_int_distribution<int> local_dist(0, 255);

        long long count_temp = 0;
        long long total_temp = 0;
        uint8_t state[4][4];

        for (int test2 = 0; test2 < TEST_NUM2; test2++)
        {
            random_state(state, local_dist, local_rng);

            // Round-1 S-box input equals the state after whitening key.
            KeyAdd(state, MasterKey);

            if (!check_conditions(state, in_cond_mask, in_cond_value))
                continue;

            uint8_t result = dot_state(in_mask, state);

            // Finish rounds 1-4 and stop after SubCell of round 5.
            for (int r = 0; r < 4; r++)
                MidoriRound(state, r, MasterKey);
            SubCell(state);

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

    cout << "Actual_NUM = " << log2(total) << endl;
    if (total == 0 || count == 0) {
        cout << "Correlation cannot be computed: total=" << total
            << ", count=" << count << endl;
        return;
    }

    double correlation_log2 = log2((double)std::llabs(count)) - log2((double)total);
    cout << "Correlation: 2^(" << fixed << setprecision(6) << correlation_log2 << ")"
        << "    count: " << count << endl;
}


int main()
{
    for (int m = 0; m < 256; m++)
    {
        for (int v = 0; v < 256; v++)
        {
            uint8_t temp = m & v;
            uint8_t t = 0;
            for (int i = 0; i < 8; i++)
            {
                t ^= (temp >> i) & 1;
            }
            dot_table[m][v] = t;
        }
    }
    init_boundary_choices(true);
    for (int i = 0; i < 64; i++)
    {
        time_t t1, t2;

        time(&t1);
        test_Midori5_conditional();
        time(&t2);

        cout << "Time measured: " << t2 - t1 << " seconds." << endl;
        cout << "--------------------------------------" << endl;
    }
    return 0;
}



