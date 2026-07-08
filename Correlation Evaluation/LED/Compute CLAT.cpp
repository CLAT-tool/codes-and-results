#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <bitset>
#include <iomanip>  
#include <sstream>
#include <algorithm>
#include <tuple>
using namespace std;





const uint8_t Sbox[16] = { 0xC,0x5,0x6,0xB,0x9,0x0,0xA,0xD,0x3,0xE,0xF,0x8,0x4,0x7,0x1,0x2 };


int main()
{

    uint8_t dot_table[16][16];
    for (int m = 0; m < 16; m++)
    {
        for (int v = 0; v < 16; v++)
        {
            uint8_t temp = m & v;
            uint8_t t = 0;
            for (int i = 0; i < 4; i++)
            {
                t ^= (temp >> i) & 1;
            }
            dot_table[m][v] = t;
        }
    }

    double results[16] = { 0 };

    for (int i = 1; i < 16; i++)
    {
        results[i] = 20;
    }


#pragma omp parallel for schedule(dynamic)
    for (int lambda_o = 1; lambda_o < 16; lambda_o++)
    {
        for (int lambda_i = 1; lambda_i < 16; lambda_i++)
        {
            int temp = 0;

            for (int x = 0; x < 16; x++)
            {
                temp += (dot_table[lambda_i][x] == dot_table[lambda_o][Sbox[x]]) ? 1 : -1;
            }


            if (temp != 0)
            {
                double log_val = -2.0 * log2(abs(temp / 16.0));
#pragma omp critical
                {
                    if (log_val < results[lambda_o])
                    {
                        results[lambda_o] = log_val;
                    }

                }

            }

        }
    }



    for (int i = 0; i < 16; i++)
    {
        cout << i << "       " << results[i] << endl;
    }


#pragma omp parallel for schedule(dynamic)
    for (int lambda_o = 1; lambda_o < 16; lambda_o++)
    {
        for (int lambda_i = 1; lambda_i < 16; lambda_i++)
        {
            for (int con0 = 1; con0 < 16; con0++)
            {

                int temp[2] = { 0 };

                for (int x = 0; x < 16; x++)
                {
                    int t0 = dot_table[con0][x];
                    temp[t0] += (dot_table[lambda_i][x] == dot_table[lambda_o][Sbox[x]]) ? 1 : -1;
                }


                for (int i = 0; i < 2; i++)
                {
                    if (temp[i] != 0)
                    {
                        double log_val = -2.0 * log2(abs(temp[i] / 8.0));
#pragma omp critical
                        {
                            if (log_val + 1 < results[lambda_o])
                            {
                                results[lambda_o] = log_val + 1;
                            }
                        }
                    }
                }
            }
        }
    }



    for (int i = 0; i < 16; i++)
    {
        cout << i << "       " << results[i] << endl;
    }


    for (int i = 0; i < 16; i++)
    {
        cout << results[i] << ",";
    }
    cout << endl;



    return 0;
}


