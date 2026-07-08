
#include <iostream>
#include<fstream>
#include<sstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <map>
#include <iomanip>  
#include <mutex>
#include <time.h>
#include <string>
#include"gurobi_c++.h"
using namespace std;

//rijndael SubBytes
void SC_Sbox(GRBModel& model, GRBVar* s, GRBVar* t)
{
    for (int cell = 0; cell < 32; cell++)
    {
        model.addConstr(s[cell] == t[cell]);
    }
}

//rijndael ShiftRows
void SR_Sbox(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    model.addConstr(t2[0] == t1[0]);
    model.addConstr(t2[4] == t1[4]);
    model.addConstr(t2[8] == t1[8]);
    model.addConstr(t2[12] == t1[12]);
    model.addConstr(t2[16] == t1[16]);
    model.addConstr(t2[20] == t1[20]);
    model.addConstr(t2[24] == t1[24]);
    model.addConstr(t2[28] == t1[28]);

    model.addConstr(t2[1] == t1[5]);
    model.addConstr(t2[5] == t1[9]);
    model.addConstr(t2[9] == t1[13]);
    model.addConstr(t2[13] == t1[17]);
    model.addConstr(t2[17] == t1[21]);
    model.addConstr(t2[21] == t1[25]);
    model.addConstr(t2[25] == t1[29]);
    model.addConstr(t2[29] == t1[1]);

    model.addConstr(t2[2] == t1[14]);
    model.addConstr(t2[14] == t1[26]);
    model.addConstr(t2[26] == t1[6]);
    model.addConstr(t2[6] == t1[18]);
    model.addConstr(t2[18] == t1[30]);
    model.addConstr(t2[30] == t1[10]);
    model.addConstr(t2[10] == t1[22]);
    model.addConstr(t2[22] == t1[2]);

    model.addConstr(t2[3] == t1[19]);
    model.addConstr(t2[19] == t1[3]);

    model.addConstr(t2[7] == t1[23]);
    model.addConstr(t2[23] == t1[7]);

    model.addConstr(t2[11] == t1[27]);
    model.addConstr(t2[27] == t1[11]);

    model.addConstr(t2[15] == t1[31]);
    model.addConstr(t2[31] == t1[15]);
}

//rijndael MixColumns
void MC_Sbox(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 8; col++)
    {
        GRBVar d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(t[4 * col + 0] + t[4 * col + 1] + t[4 * col + 2] + t[4 * col + 3] + s[4 * col + 0] + s[4 * col + 1] + s[4 * col + 2] + s[4 * col + 3] >= 5 * d);
        for (int i = 0; i < 4; i++)
        {
            model.addConstr(d >= t[4 * col + i]);
            model.addConstr(d >= s[4 * col + i]);
        }
    }
}

vector<bitset<300>> Search_Rijndael_Sbox(int round, int Sbox_solution_number)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    //env.set(GRB_DoubleParam_TimeLimit, 100);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);

    GRBVar** s = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t1 = new GRBVar * [round];        //State after the S-box
    GRBVar** t2 = new GRBVar * [round];        //State after the SR
    for (i = 0; i < round; i++)
    {
        s[i] = model.addVars(32, GRB_BINARY);
        t1[i] = model.addVars(32, GRB_BINARY);
        t2[i] = model.addVars(32, GRB_BINARY);
    }
    s[round] = model.addVars(32, GRB_BINARY);


    for (int loc = 0; loc < round; loc++)
    {
        SC_Sbox(model, s[loc], t1[loc]);
        SR_Sbox(model, t1[loc], t2[loc]);
        MC_Sbox(model, t2[loc], s[loc + 1]);
    }


    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 32; i++)
        init_mask_sum += t1[0][i];
    model.addConstr(init_mask_sum >= 1);



    GRBLinExpr cor = 0;
    for (i = 1; i < round - 1; i++)
    {
        for (j = 0; j < 32; j++)
        {
            //6
            cor += 6 * s[i][j];
        }
    }
    for (j = 0; j < 32; j++)
    {
        //3.83 
        cor += 3.83 * s[0][j];

        //3.83 
        cor += 3.83 * s[round - 1][j];
    }

    //model.addConstr(cor <= 50);
    //model.addConstr(cor >= 30);

    model.setObjective(cor, GRB_MINIMIZE);

    vector<bitset<300>> foundSolutions;

    while (foundSolutions.size() < Sbox_solution_number)
    {
        model.update();
        model.optimize();

        bitset<300> temp(0);
        GRBLinExpr cut;
        for (int i = 0; i < round; i++)
        {
            cout << "{";
            for (int n1 = 0; n1 < 32; n1++)
            {
                if (s[i][n1].get(GRB_DoubleAttr_X) > 0.5)
                {
                    temp[32 * i + n1] = 1;
                    cout << "1, ";
                    cut += (1 - s[i][n1]);
                }
                else
                {
                    temp[32 * i + n1] = 0;
                    cout << "0, ";
                    cut += s[i][n1];
                }
            }
            cout << "}," << endl;
        }
        std::cout << "\nOptimal weight: " << cor.getValue() << std::endl;
        model.addConstr(cut >= 1);
        foundSolutions.push_back(temp);

        model.addConstr(cor >= cor.getValue());
    }
    return foundSolutions;

}

void XOR5(GRBModel& model, GRBVar& in1, GRBVar& in2, GRBVar& in3, GRBVar& in4, GRBVar& in5, GRBVar& out)
{
    GRBVar t = model.addVar(0.0, 2.0, 0.0, GRB_INTEGER);
    model.addConstr(in1 + in2 + in3 + in4 + in5 == 2 * t + out);
}
void XOR11(GRBModel& model, GRBVar& in1, GRBVar& in2, GRBVar& in3, GRBVar& in4, GRBVar& in5,
    GRBVar& in6, GRBVar& in7, GRBVar& in8, GRBVar& in9, GRBVar& in10, GRBVar& in11, GRBVar& out)
{
    GRBVar t = model.addVar(0.0, 5.0, 0.0, GRB_INTEGER);
    model.addConstr(in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 == 2 * t + out);
}

// AES Sbox 
// Number of constraints : 
// Input  : a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output : b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb

void SBox8_p_null(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(1 - p >= a0);
    model.addConstr(1 - p >= a1);
    model.addConstr(1 - p >= a2);
    model.addConstr(1 - p >= a3);
    model.addConstr(1 - p >= a4);
    model.addConstr(1 - p >= a5);
    model.addConstr(1 - p >= a6);
    model.addConstr(1 - p >= a7);

    model.addConstr(1 - p >= b0);
    model.addConstr(1 - p >= b1);
    model.addConstr(1 - p >= b2);
    model.addConstr(1 - p >= b3);
    model.addConstr(1 - p >= b4);
    model.addConstr(1 - p >= b5);
    model.addConstr(1 - p >= b6);
    model.addConstr(1 - p >= b7);
}

//weight: 1/64   2^-6
//Number of constraints : 1889
void SBox_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)

    
    std::ifstream infile("milp_read_AES_Sbox_linear_1-64.txt");
    std::string line;
    int constraint_count = 0;
    while (std::getline(infile, line)) {
        
        if (line.empty() || line[0] == '#') continue;

        
        std::istringstream iss(line);
        std::string coeffs_str;
        if (!std::getline(iss, coeffs_str, '|')) continue;

        double rhs_constant;
        if (!(iss >> rhs_constant)) continue;

       
        char separator;
        if (!(iss >> separator) || separator != '|') continue;

        int has_p_flag;
        if (!(iss >> has_p_flag)) continue;

        
        std::istringstream coeff_stream(coeffs_str);
        std::vector<double> coeffs(NUM_VARS, 0.0);
        for (int i = 0; i < NUM_VARS; ++i) {
            if (!(coeff_stream >> coeffs[i])) {
                
                break;
            }
        }

        
        GRBLinExpr lhs = 0;
        for (int i = 0; i < NUM_VARS; ++i) {
            if (coeffs[i] != 0.0) {
                lhs += coeffs[i] * (*vars[i]);
            }
        }

        
        
        
        
        if (has_p_flag == 1) {
            lhs += 1 - p; 
        }

       
        model.addConstr(lhs >= rhs_constant,
            "SBox1_p0_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

//weight: 49/4096   2^-6.385
//Number of constraints : 2498
void SBox_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)

    
    std::ifstream infile("milp_read_AES_Sbox_linear_49-4096.txt");
    std::string line;
    int constraint_count = 0;
    while (std::getline(infile, line)) {
        
        if (line.empty() || line[0] == '#') continue;

        
        std::istringstream iss(line);
        std::string coeffs_str;
        if (!std::getline(iss, coeffs_str, '|')) continue;

        double rhs_constant;
        if (!(iss >> rhs_constant)) continue;

       
        char separator;
        if (!(iss >> separator) || separator != '|') continue;

        int has_p_flag;
        if (!(iss >> has_p_flag)) continue;

        
        std::istringstream coeff_stream(coeffs_str);
        std::vector<double> coeffs(NUM_VARS, 0.0);
        for (int i = 0; i < NUM_VARS; ++i) {
            if (!(coeff_stream >> coeffs[i])) {
                
                break;
            }
        }

        
        GRBLinExpr lhs = 0;
        for (int i = 0; i < NUM_VARS; ++i) {
            if (coeffs[i] != 0.0) {
                lhs += coeffs[i] * (*vars[i]);
            }
        }

        
        
        
        
        if (has_p_flag == 1) {
            lhs += 1 - p; 
        }

       
        model.addConstr(lhs >= rhs_constant,
            "SBox1_p1_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

void SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    SBox8_p_null(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    model.addConstr(p0 + p1 + p2 == 1);
}

void SB(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)
{
    for (int i = 0; i < 32; i++)
    {
        SBox(model, s[8 * i], s[8 * i + 1], s[8 * i + 2], s[8 * i + 3], s[8 * i + 4], s[8 * i + 5], s[8 * i + 6], s[8 * i + 7],
            t[8 * i], t[8 * i + 1], t[8 * i + 2], t[8 * i + 3], t[8 * i + 4], t[8 * i + 5], t[8 * i + 6], t[8 * i + 7],
            p[i][0], p[i][1], p[i][2]);
    }
}

void MC(GRBModel& model, GRBVar* t, GRBVar* s)
{
    int M_02[8][8] = {
        {0,1,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0},
        {0,0,0,1,0,0,0,0},
        {1,0,0,0,1,0,0,0},
        {1,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,1,0},
        {1,0,0,0,0,0,0,1},
        {1,0,0,0,0,0,0,0}
    };
    int I8[8][8] = {
        {1,0,0,0,0,0,0,0},
        {0,1,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0},
        {0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,1}
    };
    int M_mix[32][32] = { 0 };
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            M_mix[i][j] = M_02[i][j];
        }
        for (int j = 8; j < 16; j++)
        {
            M_mix[i][j] = M_02[i][j - 8] ^ I8[i][j - 8];
        }
        for (int j = 16; j < 24; j++)
        {
            M_mix[i][j] = I8[i][j - 16];
        }
        for (int j = 24; j < 32; j++)
        {
            M_mix[i][j] = I8[i][j - 24];
        }
    }
    for (int i = 8; i < 32; i++)
    {
        for (int j = 0; j < 32; j++)
        {
            M_mix[i][j] = M_mix[i - 8][(j + 24) % 32];
        }
    }

    int index[11];
    for (int i = 0; i < 32; i++)
    {
        int k = 0;
        for (int j = 0; j < 32; j++)
        {
            if (M_mix[j][i] == 1)
            {
                index[k] = j;
                k++;
            }
        }
        if (k == 5)
        {
            for (int j = 0; j < 8; j++) 
            {
                XOR5(model, s[index[0] + 32 * j], s[index[1] + 32 * j], s[index[2] + 32 * j], s[index[3] + 32 * j], s[index[4] + 32 * j], t[i + 32 * j]);
            }
        }
        else 
        {
            for (int j = 0; j < 8; j++)
            {
                XOR11(model, s[index[0] + 32 * j], s[index[1] + 32 * j], s[index[2] + 32 * j], s[index[3] + 32 * j], s[index[4] + 32 * j], s[index[5] + 32 * j],
                    s[index[6] + 32 * j], s[index[7] + 32 * j], s[index[8] + 32 * j], s[index[9] + 32 * j], s[index[10] + 32 * j], t[i + 32 * j]);
            }
        }

    }
}

void SR(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    for (int i = 0; i < 8; i++)
    {
        model.addConstr(t2[i] == t1[i]);
        model.addConstr(t2[32 + i] == t1[32 + i]);
        model.addConstr(t2[32 * 2 + i] == t1[32 * 2 + i]);
        model.addConstr(t2[32 * 3 + i] == t1[32 * 3 + i]);
        model.addConstr(t2[32 * 4 + i] == t1[32 * 4 + i]);
        model.addConstr(t2[32 * 5 + i] == t1[32 * 5 + i]);
        model.addConstr(t2[32 * 6 + i] == t1[32 * 6 + i]);
        model.addConstr(t2[32 * 7 + i] == t1[32 * 7 + i]);

        model.addConstr(t2[8 + i] == t1[8 + 32 + i]);
        model.addConstr(t2[8 + 32 + i] == t1[8 + 32 * 2 + i]);
        model.addConstr(t2[8 + 32 * 2 + i] == t1[8 + 32 * 3 + i]);
        model.addConstr(t2[8 + 32 * 3 + i] == t1[8 + 32 * 4 + i]);
        model.addConstr(t2[8 + 32 * 4 + i] == t1[8 + 32 * 5 + i]);
        model.addConstr(t2[8 + 32 * 5 + i] == t1[8 + 32 * 6 + i]);
        model.addConstr(t2[8 + 32 * 6 + i] == t1[8 + 32 * 7 + i]);
        model.addConstr(t2[8 + 32 * 7 + i] == t1[8 + i]);

        model.addConstr(t2[16 + i] == t1[112 + i]);
        model.addConstr(t2[112 + i] == t1[208 + i]);
        model.addConstr(t2[208 + i] == t1[48 + i]);
        model.addConstr(t2[48 + i] == t1[144 + i]);
        model.addConstr(t2[144 + i] == t1[240 + i]);
        model.addConstr(t2[240 + i] == t1[80 + i]);
        model.addConstr(t2[80 + i] == t1[176 + i]);
        model.addConstr(t2[176 + i] == t1[16 + i]);

        model.addConstr(t2[24 + i] == t1[152 + i]);
        model.addConstr(t2[152 + i] == t1[24 + i]);

        model.addConstr(t2[56 + i] == t1[184 + i]);
        model.addConstr(t2[184 + i] == t1[56 + i]);

        model.addConstr(t2[88 + i] == t1[216 + i]);
        model.addConstr(t2[216 + i] == t1[88 + i]);

        model.addConstr(t2[120 + i] == t1[248 + i]);
        model.addConstr(t2[248 + i] == t1[120 + i]);
    }
}

double Search_rijndael_cond_trail(int round, bitset<300> Template, double min)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    //env.set(GRB_DoubleParam_TimeLimit, 3600);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);
    GRBVar** s = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t1 = new GRBVar * [round + 1];        //State after the S-box
    GRBVar** t2 = new GRBVar * [round + 1];        //State after the SR
    for (i = 0; i < round + 1; i++)
    {
        s[i] = model.addVars(256, GRB_BINARY);
        t1[i] = model.addVars(256, GRB_BINARY);
        t2[i] = model.addVars(256, GRB_BINARY);
    }

    GRBVar*** p;
    p = new GRBVar * *[round + 2];
    for (i = 0; i < round + 2; i++)
    {
        p[i] = new GRBVar * [32];
    }
    for (i = 1; i < round + 1; i++)
    {
        for (j = 0; j < 32; j++)
        {
            p[i][j] = model.addVars(3, GRB_BINARY);
        }
    }
    for (j = 0; j < 32; j++)
    {
        p[0][j] = model.addVars(1, GRB_BINARY);
        p[round + 1][j] = model.addVars(1, GRB_BINARY);
        GRBVar or_vars1[] = { t1[0][8 * j], t1[0][8 * j + 1], t1[0][8 * j + 2], t1[0][8 * j + 3], t1[0][8 * j + 4], t1[0][8 * j + 5], t1[0][8 * j + 6], t1[0][8 * j + 7] };
        model.addGenConstrOr(p[0][j][0], or_vars1, 8);
        GRBVar or_vars2[] = { s[round][8 * j], s[round][8 * j + 1], s[round][8 * j + 2], s[round][8 * j + 3], s[round][8 * j + 4], s[round][8 * j + 5], s[round][8 * j + 6], s[round][8 * j + 7] };
        model.addGenConstrOr(p[round + 1][j][0], or_vars2, 8);
    }

    SR(model, t1[0], t2[0]);
    MC(model, t2[0], s[0]);
    for (int loc = 0; loc < round; loc++)
    {
        SB(model, s[loc], t1[loc + 1], p[loc + 1]);
        SR(model, t1[loc + 1], t2[loc + 1]);
        MC(model, t2[loc + 1], s[loc + 1]);
    }


    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 256; i++)
        init_mask_sum += s[0][i];
    model.addConstr(init_mask_sum >= 1);

    int sbox_index[9][32] = { 0 };

    for (i = 0; i < round + 2; i++)
    {
        for (j = 0; j < 32; j++)
        {
            sbox_index[i][j] = Template[32 * i + j];
        }
    }

    for (int cell = 0; cell < 32; cell++)
    {
        if (sbox_index[0][cell] == 0)
        {
            model.addConstr(p[0][cell][0] == 0);
        }
        else
        {
            model.addConstr(p[0][cell][0] == 1);
        }
    }
    for (int i = 1; i < round + 1; i++)
    {
        for (int cell = 0; cell < 32; cell++)
        {
            if (sbox_index[i][cell] == 0)
            {
                model.addConstr(p[i][cell][0] == 1);
            }
            else
            {
                GRBLinExpr sum = 0;
                for (int j = 1; j < 3; j++)
                    sum += p[i][cell][j];
                model.addConstr(sum == 1);
            }
        }
    }
    for (int cell = 0; cell < 32; cell++)
    {
        if (sbox_index[round + 1][cell] == 0)
        {
            model.addConstr(p[round + 1][cell][0] == 0);
        }
        else
        {
            model.addConstr(p[round + 1][cell][0] == 1);
        }
    }



    GRBLinExpr cor = 0;
    for (i = 1; i < round + 1; i++)
    {
        for (j = 0; j < 32; j++)
        {
            //6, 6.385, 6.830, 7.356, 8, 8.830, 10, 12
            cor += 6 * p[i][j][1] + 6.385 * p[i][j][2];
        }
    }
    for (j = 0; j < 32; j++)
    {
        cor += 3.83 * p[0][j][0] + 3.83 * p[round + 1][j][0];
    }

    model.addConstr(cor <= min);

    model.setObjective(cor, GRB_MINIMIZE);

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        string sss = "results.txt";
        ofstream f(sss, ios::app);
        f << "===== Rijndael " << round + 2 << "-round result of the linear distinguisher =====" << std::endl;
        f << "The S-box output mask for round 1: " << std::endl;
        for (int n1 = 0; n1 < 32; n1++)
        {
            int result = 0;
            for (int n2 = 0; n2 < 8; n2++)
            {
                if (t1[0][8 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
                    result = result * 2 + 1;
                else
                    result = result * 2;
            }
            f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
        }
        f << endl;
        f << "The S-box input mask for rounds 2 to " << round + 2 << ": " << std::endl;
        for (int i = 0; i < round + 1; i++)
        {
            for (int n1 = 0; n1 < 32; n1++)
            {
                int result = 0;
                for (int n2 = 0; n2 < 8; n2++)
                {
                    if (s[i][8 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
                        result = result * 2 + 1;
                    else
                        result = result * 2;
                }
                f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
            }
            f << endl;
        }
        f << "\nConditional linear weight: " << cor.getValue() << std::endl;
        f.close();
        return cor.getValue();
    }

    return min;
}

int main()
{
    cout << "====================================================" << endl;
    cout << "Rijndael-256 Conditional Linear Distinguisher Search" << endl;
    cout << "====================================================" << endl;
    int round;
    cout << " ROUND: ";
    cin >> round;
    cin.ignore();

    string input;
    int Sbox_solution_number = 20;
    cout << " Number of Sbox solutions attempted (default is 20): ";
    if (std::getline(std::cin, input)) {
        if (!input.empty()) {

            std::stringstream ss(input);
            if (!(ss >> Sbox_solution_number)) {

                std::cout << "Input is invalid. Default values will be used." << std::endl;
                Sbox_solution_number = 20;
            }
        }
    }


    vector<bitset<300>> foundSolutions = Search_Rijndael_Sbox(round, Sbox_solution_number);
    double min = 200;
    for (int i = 0; i < foundSolutions.size(); i++)
    {
        cout << "-------------------------- Try " << dec << i << " --------------------------" << endl;
        min = Search_rijndael_cond_trail(round - 2, foundSolutions[i], min);
    }

    cout << endl << "===================================================" << endl;
    cout << "The final Conditional linear weight:" << min << endl;
    cout << "===================================================" << endl;

    return 0;
}