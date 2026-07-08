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
//  LED 4-round conditional linear correlation experiment.
//
//  Default state convention:
//      row-major: state[row][col] = s_{4 * row + col}
//
//  If your trail/masks are generated in column-major convention:
//      state[row][col] = s_{4 * col + row}
//  then set LED_COLUMN_MAJOR = true.
// ============================================================

const bool LED_COLUMN_MAJOR = true;

// For estimating absolute correlation, constants and keys only affect signs
// or condition values. This implementation checks input conditions after
// AddConstants, so USE_ROUND_CONSTANTS can be enabled if desired.
const bool USE_ROUND_CONSTANTS = true;

// LED injects the master key before round 1 and after each group of 4 rounds.
// Set LED_KEY_NIBBLES to 16 for LED-64, or to 32 for LED-128.
const bool USE_ROUND_KEY = true;
const int LED_KEY_NIBBLES = 16;


uint8_t dot_table[16][16];

inline int cell_row(int pos)
{
    return LED_COLUMN_MAJOR ? (pos & 3) : (pos >> 2);
}

inline int cell_col(int pos)
{
    return LED_COLUMN_MAJOR ? (pos >> 2) : (pos & 3);
}

// LED uses the PRESENT 4-bit S-box.
const uint8_t led_sbox[16] = {
    0xC, 0x5, 0x6, 0xB,
    0x9, 0x0, 0xA, 0xD,
    0x3, 0xE, 0xF, 0x8,
    0x4, 0x7, 0x1, 0x2
};

// LED MixColumnsSerial matrix over GF(2^4).
const uint8_t led_mc_matrix[4][4] = {
    {0x4, 0x1, 0x2, 0x2},
    {0x8, 0x6, 0x5, 0x6},
    {0xB, 0xE, 0xA, 0x9},
    {0x2, 0x2, 0xF, 0xB}
};

// First several LED round constants.
// Only needed when USE_ROUND_CONSTANTS = true.
const uint8_t led_rc[48] = {
    0x01,0x03,0x07,0x0F,0x1F,0x3E,0x3D,0x3B,
    0x37,0x2F,0x1E,0x3C,0x39,0x33,0x27,0x0E,
    0x1D,0x3A,0x35,0x2B,0x16,0x2C,0x18,0x30,
    0x21,0x02,0x05,0x0B,0x17,0x2E,0x1C,0x38,
    0x31,0x23,0x06,0x0D,0x1B,0x36,0x2D,0x1A,
    0x34,0x29,0x12,0x24,0x08,0x11,0x22,0x04
};

uint8_t gf16_mul(uint8_t a, uint8_t b)
{
    uint8_t res = 0;

    for (int i = 0; i < 4; i++) {
        if (b & 1)
            res ^= a;

        uint8_t carry = a & 0x8;
        a = (a << 1) & 0xF;

        // x^4 = x + 1 for polynomial x^4 + x + 1.
        if (carry)
            a ^= 0x3;

        b >>= 1;
    }

    return res & 0xF;
}

void AddConstants(uint8_t state[4][4], int round)
{
    if (!USE_ROUND_CONSTANTS)
        return;

    uint8_t rc = led_rc[round] & 0x3F;
    uint8_t rc_high = (rc >> 3) & 0x7;
    uint8_t rc_low = rc & 0x7;

    // A common LED AddConstants form. If your LED reference code uses a
    // different convention, replace only this function.
    state[0][0] ^= rc_high;
    state[1][0] ^= rc_high ^ 0x1;
    state[2][0] ^= rc_low ^ 0x2;
    state[3][0] ^= rc_low ^ 0x3;

    for (int r = 0; r < 4; r++)
        state[r][0] &= 0xF;
}

void SubCells(uint8_t state[4][4])
{
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        state[r][c] = led_sbox[state[r][c] & 0xF];
    }
}

void ShiftRows(uint8_t state[4][4])
{
    uint8_t old[4][4];
    memcpy(old, state, 16);

    // LED uses AES-like ShiftRows: row r is rotated left by r nibbles.
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++)
            state[r][c] = old[r][(c + r) & 3];
    }
}

void MixColumnsSerial(uint8_t state[4][4])
{
    uint8_t old[4][4];
    memcpy(old, state, 16);

    for (int c = 0; c < 4; c++) {
        for (int r = 0; r < 4; r++) {
            uint8_t v = 0;

            for (int k = 0; k < 4; k++)
                v ^= gf16_mul(led_mc_matrix[r][k], old[k][c]);

            state[r][c] = v & 0xF;
        }
    }
}

void AddRoundKey(uint8_t state[4][4],
    const uint8_t master_key[32],
    int key_nibbles,
    int step)
{
    if (!USE_ROUND_KEY)
        return;

    // LED adds a 64-bit subkey every four rounds.
    // For LED-64, key_nibbles = 16, so the same 16 nibbles are reused.
    // For LED-128, key_nibbles = 32, so the two 64-bit halves alternate.
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        int key_pos = (16 * step + pos) % key_nibbles;
        state[r][c] ^= master_key[key_pos] & 0xF;
        state[r][c] &= 0xF;
    }
}

void random_led_key(uint8_t master_key[32],
    int key_nibbles,
    std::mt19937_64& gen)
{
    std::uniform_int_distribution<int> dist(0, 15);

    memset(master_key, 0, 32);

    for (int i = 0; i < key_nibbles; i++)
        master_key[i] = static_cast<uint8_t>(dist(gen) & 0xF);
}

void print_led_key(const uint8_t master_key[32], int key_nibbles)
{
    cout << "LED master key: ";

    for (int i = 0; i < key_nibbles; i++) {
        cout << "0x" << hex << uppercase << (int)(master_key[i] & 0xF);
        if (i + 1 != key_nibbles)
            cout << " ";
    }

    cout << nouppercase << dec << endl;
}

void LEDRound(uint8_t state[4][4], int round)
{
    AddConstants(state, round);
    SubCells(state);
    ShiftRows(state);
    MixColumnsSerial(state);
}

void LEDEncryptFullRounds(uint8_t state[4][4],
    int total_rounds,
    const uint8_t master_key[32],
    int key_nibbles)
{
    // Helper for full-round experiments. It follows the LED schedule:
    // AddRoundKey, then every 4 rounds add the next 64-bit subkey.
    AddRoundKey(state, master_key, key_nibbles, 0);

    for (int r = 0; r < total_rounds; r++) {
        LEDRound(state, r);

        if (((r + 1) % 4) == 0)
            AddRoundKey(state, master_key, key_nibbles, (r + 1) / 4);
    }
}

uint8_t dot_state(const uint8_t mask[16], const uint8_t state[4][4])
{
    uint8_t v = 0;

    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        v ^= dot_table[mask[pos] & 0xF][state[r][c] & 0xF];
    }

    return v & 1;
}

const uint8_t r1_sbox_output_mask[16] = {
0x01,0x09,0x00,0x00,0x00,0x01,0x04,0x00,0x00,0x00,0x0E,0x0B,0x01,0x00,0x00,0x06
};

const uint8_t r2_to_r4_sbox_input_mask[3][16] = {
    { // round 2 S-box input mask
0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x00
    },
    { // round 3 S-box input mask
0x00, 0x00, 0x00, 0x00, 0x00, 0x0F, 0x01, 0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    },
    { // round 4 S-box input mask
0x09, 0x02, 0x04, 0x05, 0x00, 0x00, 0x00, 0x00, 0x0B, 0x05, 0x05, 0x0B, 0x05, 0x02, 0x03, 0x0A
    }
};





// Boundary masks and conditions derived from one-condition CLAT.
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
    double weight = numeric_limits<double>::infinity();
};

// Input-side CLAT:
// fixed output mask lambda_o;
// search input mask lambda_i and condition <T, x> = t.
CLATChoice search_input_side_choice(uint8_t fixed_output_mask)
{
    CLATChoice best;

    if (fixed_output_mask == 0) {
        best.weight = 0.0;
        return best;
    }

    for (int T = 1; T < 16; T++) {
        for (int tv = 0; tv <= 1; tv++) {
            for (int li = 0; li < 16; li++) {
                int num = 0;
                int den = 0;

                for (int x = 0; x < 16; x++) {
                    if (dot_table[T][x] != tv)
                        continue;

                    den++;

                    uint8_t e = dot_table[li][x] ^ dot_table[fixed_output_mask][led_sbox[x]];

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

    for (int T = 1; T < 16; T++) {
        for (int tv = 0; tv <= 1; tv++) {
            for (int lo = 0; lo < 16; lo++) {
                int num = 0;
                int den = 0;

                for (int x = 0; x < 16; x++) {
                    uint8_t y = led_sbox[x];

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
        cout << "Input-side CLAT choices for round 1:\n";

    for (int pos = 0; pos < 16; pos++) {
        CLATChoice ch = search_input_side_choice(r1_sbox_output_mask[pos]);

        in_mask[pos] = ch.lin_mask;
        in_cond_mask[pos] = ch.cond_mask;
        in_cond_value[pos] = ch.cond_value;

        if (r1_sbox_output_mask[pos] != 0) {
            input_condition_count++;
            boundary_weight += ch.weight;

            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << ": fixed_out=0x" << hex << uppercase
                    << (int)r1_sbox_output_mask[pos]
                    << "  in=0x" << (int)ch.lin_mask
                    << "  cond=<0x" << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << nouppercase << "\n";
            }
        }
    }

    if (verbose)
        cout << "Output-side CLAT choices for round 4:\n";

    for (int pos = 0; pos < 16; pos++) {
        CLATChoice ch = search_output_side_choice(r2_to_r4_sbox_input_mask[2][pos]);

        out_mask[pos] = ch.lin_mask;
        out_cond_mask[pos] = ch.cond_mask;
        out_cond_value[pos] = ch.cond_value;

        if (r2_to_r4_sbox_input_mask[2][pos] != 0) {
            output_condition_count++;
            boundary_weight += ch.weight;

            if (verbose) {
                cout << "  cell " << setw(2) << pos
                    << ": fixed_in=0x" << hex << uppercase
                    << (int)r2_to_r4_sbox_input_mask[2][pos]
                    << "  out=0x" << (int)ch.lin_mask
                    << "  cond=<0x" << (int)ch.cond_mask
                    << "," << dec << (int)ch.cond_value << ">"
                    << "  numerator=" << ch.numerator
                    << "  weight=" << fixed << setprecision(6) << ch.weight
                    << nouppercase << "\n";
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

        uint8_t x = state[cell_row(pos)][cell_col(pos)] & 0xF;
        uint8_t lhs = dot_table[cond_mask[pos] & 0xF][x];
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
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);
        state[r][c] = static_cast<uint8_t>(dist(gen) & 0xF);
    }
}

uint8_t random_nibble_satisfying(uint8_t cond_mask,
    uint8_t cond_value,
    std::uniform_int_distribution<int>& dist,
    std::mt19937_64& gen)
{
    while (true) {
        uint8_t x = static_cast<uint8_t>(dist(gen) & 0xF);

        if (cond_mask == 0)
            return x;

        if (dot_table[cond_mask & 0xF][x] == (cond_value & 1))
            return x;
    }
}

void random_state_satisfying_input_conditions(
    uint8_t state[4][4],
    std::uniform_int_distribution<int>& dist,
    std::mt19937_64& gen)
{
    for (int pos = 0; pos < 16; pos++) {
        int r = cell_row(pos);
        int c = cell_col(pos);

        state[r][c] =
            random_nibble_satisfying(in_cond_mask[pos],
                in_cond_value[pos],
                dist,
                gen);
    }
}

// Test the full conditional approximation:
// input side: round-1 S-box input;
// output side: round-4 S-box output.
void test_LED4_conditional()
{
    // If DIRECT_INPUT_CONDITION_SAMPLING = true, samples are generated
    // uniformly from the input-side conditioned subspace.
    //
    const int TEST_NUM1 = 2048;
    const int TEST_NUM2 = 1 << 23;

    double log2_generated =
        log2((double)TEST_NUM1) + log2((double)TEST_NUM2);

    cout << "log2(TEST_NUM) = " << log2_generated << endl;


    long long count = 0;
    long long total = 0;

    uint8_t MasterKey[32];
    uint64_t key_seed =
        std::random_device{}()
        ^ ((uint64_t)time(nullptr) << 17);
    std::mt19937_64 key_rng(key_seed);
    random_led_key(MasterKey, LED_KEY_NIBBLES, key_rng);
    print_led_key(MasterKey, LED_KEY_NIBBLES);

#pragma omp parallel for
    for (int test1 = 0; test1 < TEST_NUM1; test1++) {
        uint64_t seed =
            std::random_device{}()
            ^ ((uint64_t)time(nullptr) << 1)
            ^ ((uint64_t)omp_get_thread_num() << 32)
            ^ (uint64_t)test1;

        std::mt19937_64 local_rng(seed);
        std::uniform_int_distribution<int> local_dist(0, 15);

        long long count_temp = 0;
        long long total_temp = 0;

        uint8_t state[4][4];

        for (int test2 = 0; test2 < TEST_NUM2; test2++)
        {


            random_state(state, local_dist, local_rng);
            AddRoundKey(state, MasterKey, LED_KEY_NIBBLES, 0);
            AddConstants(state, 0);

            if (!check_conditions(state, in_cond_mask, in_cond_value))
                continue;


            uint8_t result = dot_state(in_mask, state);

            // Round 1: the state is already at the S-box input.
            SubCells(state);
            ShiftRows(state);
            MixColumnsSerial(state);

            // Rounds 2 and 3 are full LED rounds.
            LEDRound(state, 1);
            LEDRound(state, 2);

            // Round 4: stop after SubCells.
            AddConstants(state, 3);
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

    if (total == 0) {
        cout << "Correlation cannot be computed: total=0" << endl;
        return;
    }

    if (count == 0) {
        cout << "Correlation: 0"
            << "    count: 0"
            << "    total: " << total << endl;
        return;
    }

    double correlation_log2 =
        log2((double)llabs(count)) - log2((double)total);

    cout << "Correlation: 2^("
        << fixed << setprecision(6) << correlation_log2 << ")"
        << "    count: " << count
        << "    total: " << total
        << endl;
}


int main()
{
    for (int m = 0; m < 16; m++) {
        for (int v = 0; v < 16; v++) {
            uint8_t temp = (uint8_t)(m & v);
            uint8_t t = 0;

            for (int i = 0; i < 4; i++)
                t ^= (temp >> i) & 1;

            dot_table[m][v] = t & 1;
        }
    }

    init_boundary_choices(true);

    for (int i = 0; i < 64; i++)
    {
        time_t t1, t2;
        time(&t1);

        test_LED4_conditional();

        time(&t2);

        cout << "Time measured: " << t2 - t1 << " seconds." << endl;
        cout << "--------------------------------------" << endl;
    }

    return 0;
}
