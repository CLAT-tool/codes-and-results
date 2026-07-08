#include <iostream>
#include <cstdint>
#include <vector>
#include <cstring>
#include <cmath>
#include <random>
#include <iomanip>
#include <omp.h> 
#include <time.h>
#include <limits>
#include <cstdlib>
using namespace std;
thread_local std::mt19937_64 rng(std::random_device{}());

uint8_t dot_table[256][256];

const uint8_t sbox[256] = {
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};

const uint8_t mul2[256] = {
    0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e,0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e,
    0x20,0x22,0x24,0x26,0x28,0x2a,0x2c,0x2e,0x30,0x32,0x34,0x36,0x38,0x3a,0x3c,0x3e,
    0x40,0x42,0x44,0x46,0x48,0x4a,0x4c,0x4e,0x50,0x52,0x54,0x56,0x58,0x5a,0x5c,0x5e,
    0x60,0x62,0x64,0x66,0x68,0x6a,0x6c,0x6e,0x70,0x72,0x74,0x76,0x78,0x7a,0x7c,0x7e,
    0x80,0x82,0x84,0x86,0x88,0x8a,0x8c,0x8e,0x90,0x92,0x94,0x96,0x98,0x9a,0x9c,0x9e,
    0xa0,0xa2,0xa4,0xa6,0xa8,0xaa,0xac,0xae,0xb0,0xb2,0xb4,0xb6,0xb8,0xba,0xbc,0xbe,
    0xc0,0xc2,0xc4,0xc6,0xc8,0xca,0xcc,0xce,0xd0,0xd2,0xd4,0xd6,0xd8,0xda,0xdc,0xde,
    0xe0,0xe2,0xe4,0xe6,0xe8,0xea,0xec,0xee,0xf0,0xf2,0xf4,0xf6,0xf8,0xfa,0xfc,0xfe,
    0x1b,0x19,0x1f,0x1d,0x13,0x11,0x17,0x15,0x0b,0x09,0x0f,0x0d,0x03,0x01,0x07,0x05,
    0x3b,0x39,0x3f,0x3d,0x33,0x31,0x37,0x35,0x2b,0x29,0x2f,0x2d,0x23,0x21,0x27,0x25,
    0x5b,0x59,0x5f,0x5d,0x53,0x51,0x57,0x55,0x4b,0x49,0x4f,0x4d,0x43,0x41,0x47,0x45,
    0x7b,0x79,0x7f,0x7d,0x73,0x71,0x77,0x75,0x6b,0x69,0x6f,0x6d,0x63,0x61,0x67,0x65,
    0x9b,0x99,0x9f,0x9d,0x93,0x91,0x97,0x95,0x8b,0x89,0x8f,0x8d,0x83,0x81,0x87,0x85,
    0xbb,0xb9,0xbf,0xbd,0xb3,0xb1,0xb7,0xb5,0xab,0xa9,0xaf,0xad,0xa3,0xa1,0xa7,0xa5,
    0xdb,0xd9,0xdf,0xdd,0xd3,0xd1,0xd7,0xd5,0xcb,0xc9,0xcf,0xcd,0xc3,0xc1,0xc7,0xc5,
    0xfb,0xf9,0xff,0xfd,0xf3,0xf1,0xf7,0xf5,0xeb,0xe9,0xef,0xed,0xe3,0xe1,0xe7,0xe5
};

const uint8_t mul3[256] = {
    0x00,0x03,0x06,0x05,0x0c,0x0f,0x0a,0x09,0x18,0x1b,0x1e,0x1d,0x14,0x17,0x12,0x11,
    0x30,0x33,0x36,0x35,0x3c,0x3f,0x3a,0x39,0x28,0x2b,0x2e,0x2d,0x24,0x27,0x22,0x21,
    0x60,0x63,0x66,0x65,0x6c,0x6f,0x6a,0x69,0x78,0x7b,0x7e,0x7d,0x74,0x77,0x72,0x71,
    0x50,0x53,0x56,0x55,0x5c,0x5f,0x5a,0x59,0x48,0x4b,0x4e,0x4d,0x44,0x47,0x42,0x41,
    0xc0,0xc3,0xc6,0xc5,0xcc,0xcf,0xca,0xc9,0xd8,0xdb,0xde,0xdd,0xd4,0xd7,0xd2,0xd1,
    0xf0,0xf3,0xf6,0xf5,0xfc,0xff,0xfa,0xf9,0xe8,0xeb,0xee,0xed,0xe4,0xe7,0xe2,0xe1,
    0xa0,0xa3,0xa6,0xa5,0xac,0xaf,0xaa,0xa9,0xb8,0xbb,0xbe,0xbd,0xb4,0xb7,0xb2,0xb1,
    0x90,0x93,0x96,0x95,0x9c,0x9f,0x9a,0x99,0x88,0x8b,0x8e,0x8d,0x84,0x87,0x82,0x81,
    0x9b,0x98,0x9d,0x9e,0x97,0x94,0x91,0x92,0x83,0x80,0x85,0x86,0x8f,0x8c,0x89,0x8a,
    0xab,0xa8,0xad,0xae,0xa7,0xa4,0xa1,0xa2,0xb3,0xb0,0xb5,0xb6,0xbf,0xbc,0xb9,0xba,
    0xfb,0xf8,0xfd,0xfe,0xf7,0xf4,0xf1,0xf2,0xe3,0xe0,0xe5,0xe6,0xef,0xec,0xe9,0xea,
    0xcb,0xc8,0xcd,0xce,0xc7,0xc4,0xc1,0xc2,0xd3,0xd0,0xd5,0xd6,0xdf,0xdc,0xd9,0xda,
    0x5b,0x58,0x5d,0x5e,0x57,0x54,0x51,0x52,0x43,0x40,0x45,0x46,0x4f,0x4c,0x49,0x4a,
    0x6b,0x68,0x6d,0x6e,0x67,0x64,0x61,0x62,0x73,0x70,0x75,0x76,0x7f,0x7c,0x79,0x7a,
    0x3b,0x38,0x3d,0x3e,0x37,0x34,0x31,0x32,0x23,0x20,0x25,0x26,0x2f,0x2c,0x29,0x2a,
    0x0b,0x08,0x0d,0x0e,0x07,0x04,0x01,0x02,0x13,0x10,0x15,0x16,0x1f,0x1c,0x19,0x1a
};


uint8_t Rcon[10] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36 };


void SubBytes(uint8_t state[4][4])
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            state[i][j] = sbox[state[i][j]];
}


void ShiftRows(uint8_t state[4][4])
{
    uint8_t temp = state[1][0];
    state[1][0] = state[1][1];
    state[1][1] = state[1][2];
    state[1][2] = state[1][3];
    state[1][3] = temp;

    swap(state[2][0], state[2][2]);
    swap(state[2][1], state[2][3]);

    temp = state[3][0];
    state[3][0] = state[3][3];
    state[3][3] = state[3][2];
    state[3][2] = state[3][1];
    state[3][1] = temp;
}

void MixColumns(uint8_t state[4][4])
{
    for (int i = 0; i < 4; ++i)
    {
        uint8_t a = state[0][i];
        uint8_t b = state[1][i];
        uint8_t c = state[2][i];
        uint8_t d = state[3][i];

        state[0][i] = mul2[a] ^ mul3[b] ^ c ^ d;
        state[1][i] = a ^ mul2[b] ^ mul3[c] ^ d;
        state[2][i] = a ^ b ^ mul2[c] ^ mul3[d];
        state[3][i] = mul3[a] ^ b ^ c ^ mul2[d];
    }
}


void KeyExpansion(const uint8_t* key, uint8_t* roundKeys)
{

    memcpy(roundKeys, key, 16);

    for (int i = 4; i < 20; i++)
    {
        uint8_t temp[4];
        memcpy(temp, &roundKeys[(i - 1) * 4], 4);

        if (i % 4 == 0)
        {
            uint8_t tmp = temp[0];
            temp[0] = temp[1];
            temp[1] = temp[2];
            temp[2] = temp[3];
            temp[3] = tmp;

            for (int j = 0; j < 4; j++)
                temp[j] = sbox[temp[j]];

            temp[0] ^= Rcon[i / 4 - 1];
        }

        for (int j = 0; j < 4; j++)
            roundKeys[i * 4 + j] = roundKeys[(i - 4) * 4 + j] ^ temp[j];
    }
}


void AddRoundKey(uint8_t state[4][4], const uint8_t* key)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            state[j][i] ^= key[i * 4 + j];
}


void aesround(int r, uint8_t state[4][4], uint8_t roundKeys[80])
{


    //AddRoundKey(state, roundKeys + 0 * 16);

    for (int i = 0; i < r - 1; i++)
    {
        SubBytes(state);
        ShiftRows(state);
        MixColumns(state);
        AddRoundKey(state, roundKeys + 16 + i * 16);
    }
    SubBytes(state);
    //AddRoundKey(state, roundKeys + 16 + r * 16);

}








// ============================================================
//  Boundary CLAT choices are derived automatically.
//
//  State convention used by this AES code:
//      state[row][col], and AddRoundKey uses AES column-major round-key bytes,
//      i.e. key[4 * col + row] is XORed to state[row][col].
//
//  The 3-round experiment starts at the round-1 S-box input and stops after
//  the round-3 SubBytes layer. The fixed masks below are the middle-trail
//  boundary masks:
//
//      round-1 S-box output masks on state[0][0], state[1][1],
//      state[2][2], state[3][3];
//
//      round-3 S-box input masks on state[0][0], state[1][0],
//      state[2][0], state[3][0].
// ============================================================

const int AES_BOUNDARY_BYTES = 4;
const int AES_CONDITION_NUM = 3;

const uint8_t r1_sbox_output_mask[AES_BOUNDARY_BYTES] = {
    0x80, 0x81, 0x01, 0x01
};

const uint8_t r3_sbox_input_mask[AES_BOUNDARY_BYTES] = {
    0x17, 0xFD, 0x5B, 0x39
};

uint8_t in_mask[AES_BOUNDARY_BYTES];
uint8_t in_cond[AES_BOUNDARY_BYTES][AES_CONDITION_NUM];
uint8_t in_conV[AES_BOUNDARY_BYTES][AES_CONDITION_NUM];

uint8_t out_mask[AES_BOUNDARY_BYTES];
uint8_t out_cond[AES_BOUNDARY_BYTES][AES_CONDITION_NUM];
uint8_t out_conV[AES_BOUNDARY_BYTES][AES_CONDITION_NUM];

struct CLATChoice {
    uint8_t lin_mask = 0;
    uint8_t cond_mask[AES_CONDITION_NUM] = { 0, 0, 0 };
    uint8_t cond_value[AES_CONDITION_NUM] = { 0, 0, 0 };
    int numerator = 0;
    int denominator = 0;
    double weight = std::numeric_limits<double>::infinity();
};

bool independent_from_chosen(const uint8_t chosen[AES_CONDITION_NUM],
    int chosen_count,
    uint8_t candidate)
{
    if (candidate == 0)
        return false;

    int span_size = 1 << chosen_count;

    for (int subset = 0; subset < span_size; subset++) {
        uint8_t v = 0;

        for (int i = 0; i < chosen_count; i++) {
            if ((subset >> i) & 1)
                v ^= chosen[i];
        }

        if (v == candidate)
            return false;
    }

    return true;
}

void score_input_side(uint8_t fixed_output_mask,
    uint8_t lin_mask,
    const uint8_t cond_mask[AES_CONDITION_NUM],
    const uint8_t cond_value[AES_CONDITION_NUM],
    int cond_count,
    int& numerator,
    int& denominator)
{
    numerator = 0;
    denominator = 0;

    for (int x = 0; x < 256; x++) {
        bool ok = true;

        for (int k = 0; k < cond_count; k++) {
            if (dot_table[cond_mask[k]][x] != (cond_value[k] & 1)) {
                ok = false;
                break;
            }
        }

        if (!ok)
            continue;

        denominator++;

        uint8_t e =
            dot_table[lin_mask][x]
            ^ dot_table[fixed_output_mask][sbox[x]];

        numerator += (e == 0) ? 1 : -1;
    }
}

void score_output_side(uint8_t fixed_input_mask,
    uint8_t lin_mask,
    const uint8_t cond_mask[AES_CONDITION_NUM],
    const uint8_t cond_value[AES_CONDITION_NUM],
    int cond_count,
    int& numerator,
    int& denominator)
{
    numerator = 0;
    denominator = 0;

    for (int x = 0; x < 256; x++) {
        uint8_t y = sbox[x];
        bool ok = true;

        for (int k = 0; k < cond_count; k++) {
            if (dot_table[cond_mask[k]][y] != (cond_value[k] & 1)) {
                ok = false;
                break;
            }
        }

        if (!ok)
            continue;

        denominator++;

        uint8_t e =
            dot_table[fixed_input_mask][x]
            ^ dot_table[lin_mask][y];

        numerator += (e == 0) ? 1 : -1;
    }
}

// Greedy three-condition search.
// This mirrors the Midori/LED/SKINNY format, but AES uses three independent
// one-bit conditions per active boundary byte.
CLATChoice search_input_side_choice(uint8_t fixed_output_mask)
{
    CLATChoice cur;

    if (fixed_output_mask == 0) {
        cur.weight = 0.0;
        return cur;
    }

    for (int stage = 0; stage < AES_CONDITION_NUM; stage++) {
        int best_abs_num = -1;
        CLATChoice best = cur;

        for (int li = 0; li < 256; li++) {
            for (int T = 1; T < 256; T++) {
                if (!independent_from_chosen(cur.cond_mask, stage, (uint8_t)T))
                    continue;

                for (int tv = 0; tv <= 1; tv++) {
                    uint8_t trial_mask[AES_CONDITION_NUM] = {
                        cur.cond_mask[0], cur.cond_mask[1], cur.cond_mask[2]
                    };
                    uint8_t trial_value[AES_CONDITION_NUM] = {
                        cur.cond_value[0], cur.cond_value[1], cur.cond_value[2]
                    };

                    trial_mask[stage] = (uint8_t)T;
                    trial_value[stage] = (uint8_t)tv;

                    int num = 0;
                    int den = 0;

                    score_input_side(fixed_output_mask,
                        (uint8_t)li,
                        trial_mask,
                        trial_value,
                        stage + 1,
                        num,
                        den);

                    int abs_num = std::abs(num);

                    if (abs_num > best_abs_num) {
                        best_abs_num = abs_num;
                        best.lin_mask = (uint8_t)li;
                        memcpy(best.cond_mask, trial_mask, AES_CONDITION_NUM);
                        memcpy(best.cond_value, trial_value, AES_CONDITION_NUM);
                        best.numerator = num;
                        best.denominator = den;
                    }
                }
            }
        }

        cur = best;
    }

    double corr = fabs((double)cur.numerator / (double)cur.denominator);
    cur.weight = -2.0 * log2(corr) + AES_CONDITION_NUM;

    return cur;
}

CLATChoice search_output_side_choice(uint8_t fixed_input_mask)
{
    CLATChoice cur;

    if (fixed_input_mask == 0) {
        cur.weight = 0.0;
        return cur;
    }

    for (int stage = 0; stage < AES_CONDITION_NUM; stage++) {
        int best_abs_num = -1;
        CLATChoice best = cur;

        for (int lo = 0; lo < 256; lo++) {
            for (int T = 1; T < 256; T++) {
                if (!independent_from_chosen(cur.cond_mask, stage, (uint8_t)T))
                    continue;

                for (int tv = 0; tv <= 1; tv++) {
                    uint8_t trial_mask[AES_CONDITION_NUM] = {
                        cur.cond_mask[0], cur.cond_mask[1], cur.cond_mask[2]
                    };
                    uint8_t trial_value[AES_CONDITION_NUM] = {
                        cur.cond_value[0], cur.cond_value[1], cur.cond_value[2]
                    };

                    trial_mask[stage] = (uint8_t)T;
                    trial_value[stage] = (uint8_t)tv;

                    int num = 0;
                    int den = 0;

                    score_output_side(fixed_input_mask,
                        (uint8_t)lo,
                        trial_mask,
                        trial_value,
                        stage + 1,
                        num,
                        den);

                    int abs_num = std::abs(num);

                    if (abs_num > best_abs_num) {
                        best_abs_num = abs_num;
                        best.lin_mask = (uint8_t)lo;
                        memcpy(best.cond_mask, trial_mask, AES_CONDITION_NUM);
                        memcpy(best.cond_value, trial_value, AES_CONDITION_NUM);
                        best.numerator = num;
                        best.denominator = den;
                    }
                }
            }
        }

        cur = best;
    }

    double corr = fabs((double)cur.numerator / (double)cur.denominator);
    cur.weight = -2.0 * log2(corr) + AES_CONDITION_NUM;

    return cur;
}

void init_boundary_choices(bool verbose = true)
{
    memset(in_mask, 0, sizeof(in_mask));
    memset(in_cond, 0, sizeof(in_cond));
    memset(in_conV, 0, sizeof(in_conV));

    memset(out_mask, 0, sizeof(out_mask));
    memset(out_cond, 0, sizeof(out_cond));
    memset(out_conV, 0, sizeof(out_conV));

    double boundary_weight = 0.0;

    if (verbose)
        cout << "Input-side CLAT choices for round 1:" << endl;

    for (int i = 0; i < AES_BOUNDARY_BYTES; i++) {
        CLATChoice ch = search_input_side_choice(r1_sbox_output_mask[i]);

        in_mask[i] = ch.lin_mask;

        for (int k = 0; k < AES_CONDITION_NUM; k++) {
            in_cond[i][k] = ch.cond_mask[k];
            in_conV[i][k] = ch.cond_value[k];
        }

        boundary_weight += ch.weight;

        if (verbose) {
            cout << "  byte " << i
                << ": fixed_out=0x" << hex << uppercase << setw(2) << setfill('0')
                << (int)r1_sbox_output_mask[i]
                << "  in=0x" << setw(2) << (int)ch.lin_mask
                << "  conds=";

            for (int k = 0; k < AES_CONDITION_NUM; k++) {
                cout << "<0x" << setw(2) << (int)ch.cond_mask[k]
                    << "," << dec << (int)ch.cond_value[k] << ">";
                if (k + 1 != AES_CONDITION_NUM)
                    cout << ",";
                cout << hex << uppercase << setfill('0');
            }

            cout << dec << setfill(' ')
                << "  numerator=" << ch.numerator
                << "  denominator=" << ch.denominator
                << "  weight=" << fixed << setprecision(6) << ch.weight
                << nouppercase << endl;
        }
    }

    if (verbose)
        cout << "Output-side CLAT choices for round 3:" << endl;

    for (int i = 0; i < AES_BOUNDARY_BYTES; i++) {
        CLATChoice ch = search_output_side_choice(r3_sbox_input_mask[i]);

        out_mask[i] = ch.lin_mask;

        for (int k = 0; k < AES_CONDITION_NUM; k++) {
            out_cond[i][k] = ch.cond_mask[k];
            out_conV[i][k] = ch.cond_value[k];
        }

        boundary_weight += ch.weight;

        if (verbose) {
            cout << "  byte " << i
                << ": fixed_in=0x" << hex << uppercase << setw(2) << setfill('0')
                << (int)r3_sbox_input_mask[i]
                << "  out=0x" << setw(2) << (int)ch.lin_mask
                << "  conds=";

            for (int k = 0; k < AES_CONDITION_NUM; k++) {
                cout << "<0x" << setw(2) << (int)ch.cond_mask[k]
                    << "," << dec << (int)ch.cond_value[k] << ">";
                if (k + 1 != AES_CONDITION_NUM)
                    cout << ",";
                cout << hex << uppercase << setfill('0');
            }

            cout << dec << setfill(' ')
                << "  numerator=" << ch.numerator
                << "  denominator=" << ch.denominator
                << "  weight=" << fixed << setprecision(6) << ch.weight
                << nouppercase << endl;
        }
    }

}


void print_key(const uint8_t key[16])
{
    cout << "AES-128 master key: ";

    for (int i = 0; i < 16; i++) {
        cout << "0x" << hex << uppercase << setw(2) << setfill('0')
            << (int)key[i];
        if (i + 1 != 16)
            cout << " ";
    }

    cout << dec << setfill(' ') << nouppercase << endl;
}


void test_3AES()
{
    const int TEST_NUM1 = 2048;
    const int TEST_NUM2 = 1 << 29;

    cout << "log2(TEST_NUM) =  "
        << log2((double)TEST_NUM1) + log2((double)TEST_NUM2)
        << endl;

    long long count = 0;
    long long total = 0;

    uint64_t key_seed =
        std::random_device{}()
        ^ ((uint64_t)time(nullptr) << 17);

    std::mt19937_64 key_rng(key_seed);
    std::uniform_int_distribution<int> key_dist(0, 255);

    uint8_t key[16];
    for (int k = 0; k < 16; k++)
        key[k] = (uint8_t)key_dist(key_rng);

    print_key(key);

    uint8_t roundKeys[80];
    KeyExpansion(key, roundKeys);

#pragma omp parallel for
    for (int test1 = 0; test1 < TEST_NUM1; test1++)
    {
        uint64_t seed =
            std::random_device{}()
            ^ ((uint64_t)time(nullptr) << 1)
            ^ ((uint64_t)omp_get_thread_num() << 32)
            ^ (uint64_t)test1;

        std::mt19937_64 local_rng(seed);
        std::uniform_int_distribution<int> local_dist(0, 255);

        uint8_t state[4][4];
        long long count_temp = 0;
        long long total_temp = 0;

        for (int test2 = 0; test2 < TEST_NUM2; test2++)
        {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++)
                    state[i][j] = (uint8_t)local_dist(local_rng);
            }

            // Convert plaintext to the round-1 S-box input.
            // The previous version effectively sampled this internal state
            // directly. XORing a random key preserves uniformity but also tests
            // the AddRoundKey implementation in the experiment path.
            AddRoundKey(state, roundKeys + 0 * 16);

            int flag = 0;

            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < AES_CONDITION_NUM; k++) {
                    flag += (dot_table[in_cond[j][k]][state[j][j]]
                        ^ in_conV[j][k]);
                }

                if (flag > 0)
                    break;
            }

            if (flag == 0)
            {
                uint8_t result =
                    dot_table[in_mask[0]][state[0][0]]
                    ^ dot_table[in_mask[1]][state[1][1]]
                    ^ dot_table[in_mask[2]][state[2][2]]
                    ^ dot_table[in_mask[3]][state[3][3]];

                aesround(3, state, roundKeys);

                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < AES_CONDITION_NUM; k++) {
                        flag += (dot_table[out_cond[j][k]][state[j][0]]
                            ^ out_conV[j][k]);
                    }

                    if (flag > 0)
                        break;
                }

                if (flag == 0)
                {
                    total_temp++;

                    result =
                        result
                        ^ dot_table[out_mask[0]][state[0][0]]
                        ^ dot_table[out_mask[1]][state[1][0]]
                        ^ dot_table[out_mask[2]][state[2][0]]
                        ^ dot_table[out_mask[3]][state[3][0]];

                    if (result == 0)
                        count_temp++;
                    else
                        count_temp--;
                }
            }
        }

#pragma omp critical
        {
            count += count_temp;
            total += total_temp;
        }
    }

    if (total == 0) {
        cout << "Correlation cannot be computed: total=0" << endl;
        return;
    }

    cout << "log2(Actual_NUM) =  " << log2((double)total) << endl;

    if (count == 0) {
        cout << "Correlation: 0"
            << "    count: 0"
            << "    total: " << total
            << endl;
        return;
    }

    double correlation = log2((double)std::llabs(count)) - log2((double)total);

    cout << "Correlation: 2^("
        << fixed << setprecision(6) << correlation << ")"
        << "    count: " << count
        << "    total: " << total
        << endl;
}




int main()
{

    for (int m = 0; m < 256; m++) {
        for (int v = 0; v < 256; v++) {
            uint8_t temp = (uint8_t)(m & v);
            uint8_t t = 0;

            for (int i = 0; i < 8; i++)
                t ^= (temp >> i) & 1;

            dot_table[m][v] = t & 1;
        }
    }

    init_boundary_choices(true);

    
    for (int i = 0; i < 64; i++)
    {
        time_t t1, t2;
        time(&t1);

        test_3AES();

        time(&t2);

        cout << "Time measured: " << t2 - t1 << " seconds." << endl;
        cout << "--------------------------------------" << endl;
    }

    return 0;
}
