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

void XOR2(GRBModel& model, GRBVar& a, GRBVar& b, GRBVar& o)
{
    GRBVar t = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    model.addConstr(t >= a);
    model.addConstr(t >= b);
    model.addConstr(t >= o);
    model.addConstr(a + b + o >= 2 * t);
}

void XOR3(GRBModel& model, GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& o)
{
    GRBVar t = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    model.addConstr(a + b + c == 2 * t + o);
}

//Midori SubCell_Sbox
void SC_Sbox(GRBModel& model, GRBVar* s, GRBVar* t)
{
    for (int i = 0; i < 16; i++)
    {
        model.addConstr(s[i] == t[i]);
    }
}

//Midori ShuffleCell_Sbox
void SR_Sbox(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    model.addConstr(t2[0] == t1[0]);
    model.addConstr(t2[1] == t1[10]);
    model.addConstr(t2[2] == t1[5]);
    model.addConstr(t2[3] == t1[15]);

    model.addConstr(t2[4] == t1[14]);
    model.addConstr(t2[5] == t1[4]);
    model.addConstr(t2[6] == t1[11]);
    model.addConstr(t2[7] == t1[1]);

    model.addConstr(t2[8] == t1[9]);
    model.addConstr(t2[9] == t1[3]);
    model.addConstr(t2[10] == t1[12]);
    model.addConstr(t2[11] == t1[6]);

    model.addConstr(t2[12] == t1[7]);
    model.addConstr(t2[13] == t1[13]);
    model.addConstr(t2[14] == t1[2]);
    model.addConstr(t2[15] == t1[8]);


}

void MC_constraint(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7)
{
    model.addConstr(x1 + x2 + x3 - x4 >= 0);
    model.addConstr(x0 + x2 + x3 - x5 >= 0);
    model.addConstr(x0 + x1 + x3 - x6 >= 0);
    model.addConstr(x1 + x2 - x3 + x4 >= 0);
    model.addConstr(x0 + x2 - x3 + x5 >= 0);
    model.addConstr(x0 + x1 - x3 + x6 >= 0);
    model.addConstr(x0 + x1 + x2 - x7 >= 0);
    model.addConstr(x1 - x2 + x3 + x4 >= 0);
    model.addConstr(x0 - x2 + x3 + x5 >= 0);
    model.addConstr(x0 + x1 - x2 + x4 - x5 >= -1);
    model.addConstr(x0 + x1 - x2 - x4 + x5 >= -1);
    model.addConstr(x0 + x1 - x2 + x7 >= 0);
    model.addConstr(-x1 + x2 + x3 + x4 >= 0);
    model.addConstr(x0 - x1 + x3 + x6 >= 0);
    model.addConstr(x0 - x1 + x2 + x4 - x6 >= -1);
    model.addConstr(x0 - x1 + x2 - x4 + x6 >= -1);
    model.addConstr(x0 - x1 + x2 + x7 >= 0);
    model.addConstr(x0 - x1 + x3 - x4 + x7 >= -1);
    model.addConstr(x0 - x1 + x5 - x6 + x7 >= -1);
    model.addConstr(x0 - x1 - x5 + x6 + x7 >= -1);
    model.addConstr(x0 - x1 - x3 + x4 + x7 >= -1);
    model.addConstr(x0 - x1 + x3 + x4 - x7 >= -1);
    model.addConstr(x0 - x1 - x2 + x4 + x5 - x7 >= -2);
    model.addConstr(x0 - x1 - x2 + x4 + x6 - x7 >= -2);
    model.addConstr(x0 - x1 - x2 + x5 + x6 - x7 >= -2);
    model.addConstr(-x0 + x2 + x3 + x5 >= 0);
    model.addConstr(-x0 + x1 + x3 + x6 >= 0);
    model.addConstr(-x0 + x1 + x2 + x5 - x6 >= -1);
    model.addConstr(-x0 + x1 + x2 - x5 + x6 >= -1);
    model.addConstr(-x0 + x1 + x2 + x7 >= 0);
    model.addConstr(-x0 + x1 + x3 - x5 + x7 >= -1);
    model.addConstr(-x0 + x1 + x4 - x6 + x7 >= -1);
    model.addConstr(-x0 + x1 - x4 + x6 + x7 >= -1);
    model.addConstr(-x0 + x1 - x3 + x5 + x7 >= -1);
    model.addConstr(-x0 + x1 + x3 + x5 - x7 >= -1);
    model.addConstr(-x0 + x1 - x2 + x4 + x5 - x7 >= -2);
    model.addConstr(-x0 + x1 - x2 + x4 + x6 - x7 >= -2);
    model.addConstr(-x0 + x1 - x2 + x5 + x6 - x7 >= -2);
    model.addConstr(-x0 + x2 + x3 - x6 + x7 >= -1);
    model.addConstr(-x0 + x2 + x4 - x5 + x7 >= -1);
    model.addConstr(-x0 + x2 - x4 + x5 + x7 >= -1);
    model.addConstr(-x0 + x2 - x3 + x6 + x7 >= -1);
    model.addConstr(-x0 + x2 + x3 + x6 - x7 >= -1);
    model.addConstr(-x0 - x1 + x2 + x4 + x5 - x7 >= -2);
    model.addConstr(-x0 - x1 + x2 + x4 + x6 - x7 >= -2);
    model.addConstr(-x0 - x1 + x2 + x5 + x6 - x7 >= -2);
    model.addConstr(-x0 - x1 + x3 + x4 + x5 - x6 >= -2);
    model.addConstr(-x0 + x3 + x4 - x5 + x6 >= -1);
    model.addConstr(-x0 - x1 + x3 + x4 + x7 >= -1);
    model.addConstr(-x0 + x3 - x4 + x5 + x6 >= -1);
    model.addConstr(-x0 - x1 + x3 + x5 + x7 >= -1);
    model.addConstr(-x0 - x2 + x3 + x6 + x7 >= -1);
    model.addConstr(-x0 - x1 - x2 - x3 + x4 + x5 + x6 >= -3);
    model.addConstr(-x0 - x1 - x2 + x4 + x5 + x7 >= -2);
    model.addConstr(-x0 - x1 - x2 + x4 + x6 + x7 >= -2);
    model.addConstr(-x0 - x1 - x2 + x5 + x6 + x7 >= -2);

}

//Midori MixColumn_Sbox
void MC_Sbox(GRBModel& model, GRBVar* t, GRBVar* s)
{
    //GRBLinExpr branch_number[4];
    GRBLinExpr branch_number[4];

    for (int col = 0; col < 4; col++)
    {
        GRBVar d0 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        XOR2(model, s[4 * col + 2], s[4 * col + 3], d0);        //d = s2 + s3
        XOR2(model, s[4 * col + 1], d0, t[4 * col + 0]);        //t0 = s1 + s2 + s3 = s0 + d
        XOR2(model, s[4 * col + 0], d0, t[4 * col + 1]);        //t1 = s0 + s2 + s3 = s0 + d

        GRBVar d1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        XOR2(model, s[4 * col + 0], s[4 * col + 1], d1);        //d1 = s0 + s1
        XOR2(model, s[4 * col + 3], d1, t[4 * col + 2]);        //t2 = s0 + s1 + s3 = s3 + d
        XOR2(model, s[4 * col + 2], d1, t[4 * col + 3]);        //t3 = s0 + s1 + s2 = s2 + d

        //GRBVar d2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 1], s[4 * col + 2], d2);      //d2 = s1 + s2
        //XOR2(model, s[4 * col + 3], d2, t[4 * col + 0]);      //t0 = s1 + s2 + s3 = s3 + d
        //XOR2(model, s[4 * col + 0], d2, t[4 * col + 3]);      //t3 = s0 + s1 + s2 = s0 + d

        //GRBVar d3 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 0], s[4 * col + 2], d3);      //d3 = s0 + s2
        //XOR2(model, s[4 * col + 3], d3, t[4 * col + 1]);      //t1 = s0 + s2 + s3 = s3 + d
        //XOR2(model, s[4 * col + 1], d3, t[4 * col + 3]);      //t3 = s0 + s1 + s2 = s1 + d

        //GRBVar d4 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 0], s[4 * col + 3], d4);      //d4 = s0 + s3
        //XOR2(model, s[4 * col + 2], d4, t[4 * col + 1]);      //t1 = s0 + s2 + s3 = s2 + d
        //XOR2(model, s[4 * col + 1], d4, t[4 * col + 2]);      //t2 = s0 + s1 + s3 = s1 + d

        //GRBVar d5 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 1], s[4 * col + 3], d5);      //d4 = s1 + s3
        //XOR2(model, s[4 * col + 1], d5, t[4 * col + 0]);      //t0 = s1 + s2 + s3 = s1 + d
        //XOR2(model, s[4 * col + 0], d5, t[4 * col + 2]);      //t2 = s0 + s1 + s3 = s0 + d
        //
        GRBVar d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        branch_number[col] = t[4 * col] + t[4 * col + 1] + t[4 * col + 2] + t[4 * col + 3] + s[4 * col] + s[4 * col + 1] + s[4 * col + 2] + s[4 * col + 3];
        branch_number[col] -= 4 * d;
        model.addConstr(branch_number[col] >= 0);
        model.addConstr(t[4 * col] + t[4 * col + 1] + t[4 * col + 2] + t[4 * col + 3] - d >= 0);
        model.addConstr(s[4 * col] + s[4 * col + 1] + s[4 * col + 2] + s[4 * col + 3] - d >= 0);
        for (int j = 0; j < 4; j++)
        {
            model.addConstr(d >= t[4 * col + j]);
            model.addConstr(d >= s[4 * col + j]);
        }
    }
    for (int col = 0; col < 4; col++)
    {
        MC_constraint(model, t[4 * col], t[4 * col + 1], t[4 * col + 2], t[4 * col + 3], s[4 * col], s[4 * col + 1], s[4 * col + 2], s[4 * col + 3]);

        //GRBVar** d = new GRBVar * [8]; 
        //for (int cell = 0; cell < 4; cell++)
        //{
        //    d[cell] = model.addVars(4, GRB_BINARY);
        //    OR5(model, d[cell][0], d[cell][1], d[cell][2], d[cell][3], t[4 * col + cell]);
        //}
        //for (int cell = 4; cell < 8; cell++)
        //{
        //    d[cell] = model.addVars(4, GRB_BINARY);
        //    OR5(model, d[cell][0], d[cell][1], d[cell][2], d[cell][3], s[4 * col + cell - 4]);
        //}

        //for (int bit = 0; bit < 4; bit++)
        //{
        //    XOR3(model, d[1][bit], d[2][bit], d[3][bit], d[4][bit]);
        //    XOR3(model, d[0][bit], d[2][bit], d[3][bit], d[5][bit]);
        //    XOR3(model, d[0][bit], d[1][bit], d[3][bit], d[6][bit]);
        //    XOR3(model, d[0][bit], d[1][bit], d[2][bit], d[7][bit]);
        //}
        //
        //branch_number[i] = t[4 * i] + t[4 * i + 1] + t[4 * i + 2] + t[4 * i + 3] + s[4 * i] + s[4 * i + 1] + s[4 * i + 2] + s[4 * i + 3];
        //branch_number[i] -= 4 * d[i];
        //model.addConstr(branch_number[i] >= 0);
        //model.addConstr(t[4 * i] + t[4 * i + 1] + t[4 * i + 2] + t[4 * i + 3] - d[i] >= 0);
        //model.addConstr(s[4 * i] + s[4 * i + 1] + s[4 * i + 2] + s[4 * i + 3] - d[i] >= 0);
    }

}

vector<bitset<300>> Search_Midori_Sbox(int round)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    //env.set(GRB_DoubleParam_TimeLimit, 100);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);
    GRBVar** s = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t1 = new GRBVar * [round];        //State after the S-box
    GRBVar** t2 = new GRBVar * [round];        //State after the SR
    for (i = 0; i < round; i++)
    {
        s[i] = model.addVars(16, GRB_BINARY);
        for (j = 0; j < 16; j++)
        {
            s[i][j].set(GRB_IntAttr_PoolIgnore, 1);
        }
        t1[i] = model.addVars(16, GRB_BINARY);
        t2[i] = model.addVars(16, GRB_BINARY);
    }
    s[round] = model.addVars(16, GRB_BINARY);





    for (int loc = 0; loc < round; loc++)
    {
        SC_Sbox(model, s[loc], t1[loc]);
        SR_Sbox(model, t1[loc], t2[loc]);
        MC_Sbox(model, t2[loc], s[loc + 1]);
    }

    /*
    GRBLinExpr constr = 0;
    for (i = 0; i < 16; i++)
        constr += p0[round / 2][i] + p0[round / 2][i];
    model.addConstr(constr == 1);
    */

    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 16; i++)
        init_mask_sum += s[0][i];
    model.addConstr(init_mask_sum >= 1);

    /*GRBLinExpr mask_sum = 0;
    for (i = 0; i < 16; i++)
        mask_sum += s[int(round/2)][i];
    model.addConstr(mask_sum >= 1);*/

    GRBLinExpr condition_cor = 0;
    for (i = 1; i < round - 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            condition_cor += 2 * s[i][j];
        }
    }
    for (j = 0; j < 16; j++)
    {
        condition_cor += s[0][j];
        condition_cor += s[round - 1][j];
    }


    model.setObjective(condition_cor, GRB_MINIMIZE);

    vector<bitset<300>> foundSolutions;  

    while (foundSolutions.size() < 10)
    {
        model.update();
        model.optimize();

        bitset<300> temp(0);
        GRBLinExpr cut;
        for (int i = 0; i < round; i++)
        {
            cout << "{";
            for (int n1 = 0; n1 < 16; n1++)
            {
                if (s[i][n1].get(GRB_DoubleAttr_X) > 0.5)
                {
                    temp[16 * i + n1] = 1;
                    cout << "1, ";
                    cut += (1 - s[i][n1]);
                }
                else
                {
                    temp[16 * i + n1] = 0;
                    cout << "0, ";
                    cut += s[i][n1];
                }
            }
            cout << "}," << endl;
        }
        std::cout << "\nOptimal weight: " << condition_cor.getValue() << std::endl;
        model.addConstr(cut >= 1);
        foundSolutions.push_back(temp);

        model.addConstr(condition_cor >= condition_cor.getValue());
    }


    return foundSolutions;



}

void SBox_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
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

//weight: 1/4   2^-2
//Number of constraints : 73
void SBox0_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox0_linear_1-4.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox0_p1 constraint_num " << constraint_num << " (should be 73)" << endl;
    infile.close();

}
//weight: 1/16   2^-4
//Number of constraints : 268
void SBox0_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox0_linear_1-16.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox0_p2 constraint_num " << constraint_num << " (should be 268)" << endl;
    infile.close();

}
//weight: 1/64   2^-6
//Number of constraints : 808
void SBox0_p3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox0_linear_1-64.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox0_p3 constraint_num " << constraint_num << " (should be 808)" << endl;
    infile.close();

}
//weight: 1/256   2^-8
//Number of constraints : 82
void SBox0_p4(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox0_linear_1-256.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox0_p4 constraint_num " << constraint_num << " (should be 82)" << endl;
    infile.close();

}

//weight: 1/4   2^-2
//Number of constraints : 73
void SBox1_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox1_linear_1-4.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox1_p1 constraint_num " << constraint_num << " (should be 73)" << endl;
    infile.close();

}
//weight: 1/16   2^-4
//Number of constraints : 268
void SBox1_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox1_linear_1-16.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox1_p2 constraint_num " << constraint_num << " (should be 268)" << endl;
    infile.close();

}
//weight: 1/64   2^-6
//Number of constraints : 823
void SBox1_p3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox1_linear_1-64.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox1_p3 constraint_num " << constraint_num << " (should be 823)" << endl;
    infile.close();

}
//weight: 1/256   2^-8
//Number of constraints : 82
void SBox1_p4(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox1_linear_1-256.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox1_p4 constraint_num " << constraint_num << " (should be 82)" << endl;
    infile.close();

}

//weight: 1/4   2^-2
//Number of constraints : 71
void SBox2_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox2_linear_1-4.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox2_p1 constraint_num " << constraint_num << " (should be 71)" << endl;
    infile.close();

}
//weight: 1/16   2^-4
//Number of constraints : 269
void SBox2_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox2_linear_1-16.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox2_p2 constraint_num " << constraint_num << " (should be 269)" << endl;
    infile.close();

}
//weight: 1/64   2^-6
//Number of constraints : 818
void SBox2_p3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox2_linear_1-64.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox2_p3 constraint_num " << constraint_num << " (should be 818)" << endl;
    infile.close();

}
//weight: 1/256   2^-8
//Number of constraints : 82
void SBox2_p4(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox2_linear_1-256.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox2_p4 constraint_num " << constraint_num << " (should be 82)" << endl;
    infile.close();

}

//weight: 1/4   2^-2
//Number of constraints : 71
void SBox3_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox3_linear_1-4.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox3_p1 constraint_num " << constraint_num << " (should be 71)" << endl;
    infile.close();

}
//weight: 1/16   2^-4
//Number of constraints : 270
void SBox3_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox3_linear_1-16.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox3_p2 constraint_num " << constraint_num << " (should be 270)" << endl;
    infile.close();

}
//weight: 1/64   2^-6
//Number of constraints : 809
void SBox3_p3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox3_linear_1-64.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox3_p3 constraint_num " << constraint_num << " (should be 809)" << endl;
    infile.close();

}
//weight: 1/256   2^-8
//Number of constraints : 82
void SBox3_p4(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)
    int constraint_num = 0;
    
    std::ifstream infile("milp_read_Midori_Sbox3_linear_1-256.txt");
    std::string line;
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

        
        model.addConstr(lhs >= rhs_constant);
        constraint_num++;
    }
    //cout << "SBox3_p4 constraint_num " << constraint_num << " (should be 82)" << endl;
    infile.close();

}


//iCLAT

void cond_SBox_p0(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x7 >= 0 - (1 - p));
    model.addConstr(-x6 >= 0 - (1 - p));
    model.addConstr(-x5 >= 0 - (1 - p));
    model.addConstr(-x4 >= 0 - (1 - p));
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(-x1 >= 0 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 1.000
void cond_SBox0_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x6 - x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x1 - x7 >= -1 - (1 - p));
    model.addConstr(-x5 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x6 >= -1 - (1 - p));
    model.addConstr(-x0 - x6 >= -1 - (1 - p));
    model.addConstr(-x4 - x5 >= -1 - (1 - p));
    model.addConstr(-x3 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 >= -1 - (1 - p));
    model.addConstr(-x2 - x4 >= -1 - (1 - p));
    model.addConstr(-x0 - x4 >= -1 - (1 - p));
    model.addConstr(-x2 - x3 >= -1 - (1 - p));
    model.addConstr(-x0 - x3 >= -1 - (1 - p));
    model.addConstr(-x1 - x2 >= -1 - (1 - p));
    model.addConstr(-x0 - x1 >= -1 - (1 - p));

}
//weight: 2.000
void cond_SBox0_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x2 + x5 + x7 >= 1 - (1 - p));
    model.addConstr(x1 + x3 + x4 + x6 >= 1 - (1 - p));
}
//weight: 1.000
void cond_SBox1_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x5 - x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x2 - x7 >= -1 - (1 - p));
    model.addConstr(-x5 - x6 >= -1 - (1 - p));
    model.addConstr(-x4 - x6 >= -1 - (1 - p));
    model.addConstr(-x3 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x6 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 >= -1 - (1 - p));
    model.addConstr(-x0 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x4 >= -1 - (1 - p));
    model.addConstr(-x0 - x4 >= -1 - (1 - p));
    model.addConstr(-x1 - x3 >= -1 - (1 - p));
    model.addConstr(-x0 - x3 >= -1 - (1 - p));
    model.addConstr(-x1 - x2 >= -1 - (1 - p));
    model.addConstr(-x0 - x2 >= -1 - (1 - p));

}
//weight: 2.000
void cond_SBox1_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x1 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x2 + x3 + x4 + x5 >= 1 - (1 - p));
}
//weight: 1.000
void cond_SBox2_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x4 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x2 - x7 >= -1 - (1 - p));
    model.addConstr(-x1 - x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x6 >= -1 - (1 - p));
    model.addConstr(-x3 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x6 >= -1 - (1 - p));
    model.addConstr(-x1 - x6 >= -1 - (1 - p));
    model.addConstr(-x4 - x5 >= -1 - (1 - p));
    model.addConstr(-x3 - x5 >= -1 - (1 - p));
    model.addConstr(-x2 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 >= -1 - (1 - p));
    model.addConstr(-x0 - x4 >= -1 - (1 - p));
    model.addConstr(-x0 - x3 >= -1 - (1 - p));
    model.addConstr(-x0 - x2 >= -1 - (1 - p));
    model.addConstr(-x0 - x1 >= -1 - (1 - p));
}
//weight: 2.000
void cond_SBox2_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x1 + x2 + x3 + x4 >= 1 - (1 - p));
}
//weight: 1.000
void cond_SBox3_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x6 - x7 >= -1 - (1 - p));
    model.addConstr(-x5 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x0 - x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x6 >= -1 - (1 - p));
    model.addConstr(-x1 - x6 >= -1 - (1 - p));
    model.addConstr(-x4 - x5 >= -1 - (1 - p));
    model.addConstr(-x2 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 >= -1 - (1 - p));
    model.addConstr(-x3 - x4 >= -1 - (1 - p));
    model.addConstr(-x0 - x4 >= -1 - (1 - p));
    model.addConstr(-x2 - x3 >= -1 - (1 - p));
    model.addConstr(-x1 - x3 >= -1 - (1 - p));
    model.addConstr(-x0 - x2 >= -1 - (1 - p));
    model.addConstr(-x0 - x1 >= -1 - (1 - p));
}
//weight: 2.000
void cond_SBox3_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x1 + x2 + x4 + x7 >= 1 - (1 - p));
    model.addConstr(x0 + x3 + x5 + x6 >= 1 - (1 - p));
}

void SBox0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4)
{
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox0_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox0_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    SBox0_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p3);
    SBox0_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p4);

    model.addConstr(p0 + p1 + p2 + p3 + p4 == 1);
}
void SBox1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4)
{
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox1_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox1_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    SBox1_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p3);
    SBox1_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p4);

    model.addConstr(p0 + p1 + p2 + p3 + p4 == 1);
}
void SBox2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4)
{
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox2_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox2_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    SBox2_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p3);
    SBox2_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p4);

    model.addConstr(p0 + p1 + p2 + p3 + p4 == 1);
}
void SBox3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4)
{
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox3_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox3_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    SBox3_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p3);
    SBox3_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p4);

    model.addConstr(p0 + p1 + p2 + p3 + p4 == 1);
}

void cond_SBox0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_SBox0_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_SBox0_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);

    model.addConstr(p0 + p1 + p2 == 1);
}
void cond_SBox1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_SBox1_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_SBox1_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);

    model.addConstr(p0 + p1 + p2 == 1);
}
void cond_SBox2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_SBox2_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_SBox2_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);

    model.addConstr(p0 + p1 + p2 == 1);
}
void cond_SBox3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_SBox3_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_SBox3_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);

    model.addConstr(p0 + p1 + p2 == 1);
}


//Midori SubCell
void SC(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        if (cell % 4 == 0)
        {
            SBox0(model, s[8 * cell + 0], s[8 * cell + 1], s[8 * cell + 2], s[8 * cell + 3], s[8 * cell + 4], s[8 * cell + 5], s[8 * cell + 6], s[8 * cell + 7],
                t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4]);
        }
        else if (cell % 4 == 1)
        {
            SBox1(model, s[8 * cell + 0], s[8 * cell + 1], s[8 * cell + 2], s[8 * cell + 3], s[8 * cell + 4], s[8 * cell + 5], s[8 * cell + 6], s[8 * cell + 7],
                t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4]);
        }
        else if (cell % 4 == 2)
        {
            SBox2(model, s[8 * cell + 0], s[8 * cell + 1], s[8 * cell + 2], s[8 * cell + 3], s[8 * cell + 4], s[8 * cell + 5], s[8 * cell + 6], s[8 * cell + 7],
                t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4]);
        }
        else
        {
            SBox3(model, s[8 * cell + 0], s[8 * cell + 1], s[8 * cell + 2], s[8 * cell + 3], s[8 * cell + 4], s[8 * cell + 5], s[8 * cell + 6], s[8 * cell + 7],
                t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4]);
        }
    }
}
void cond_SC(GRBModel& model, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        if (cell % 4 == 0)
        {
            cond_SBox0(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2]);
        }
        else if (cell % 4 == 1)
        {
            cond_SBox1(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2]);
        }
        else if (cell % 4 == 2)
        {
            cond_SBox2(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2]);
        }
        else
        {
            cond_SBox3(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
                p[cell][0], p[cell][1], p[cell][2]);
        }
    }
}

//Midori ShuffleCell
void SR(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    for (int bit = 0; bit < 8; bit++)
    {
        model.addConstr(t2[0 * 8 + bit] == t1[0 * 8 + bit]);
        model.addConstr(t2[1 * 8 + bit] == t1[10 * 8 + bit]);
        model.addConstr(t2[2 * 8 + bit] == t1[5 * 8 + bit]);
        model.addConstr(t2[3 * 8 + bit] == t1[15 * 8 + bit]);

        model.addConstr(t2[4 * 8 + bit] == t1[14 * 8 + bit]);
        model.addConstr(t2[5 * 8 + bit] == t1[4 * 8 + bit]);
        model.addConstr(t2[6 * 8 + bit] == t1[11 * 8 + bit]);
        model.addConstr(t2[7 * 8 + bit] == t1[1 * 8 + bit]);

        model.addConstr(t2[8 * 8 + bit] == t1[9 * 8 + bit]);
        model.addConstr(t2[9 * 8 + bit] == t1[3 * 8 + bit]);
        model.addConstr(t2[10 * 8 + bit] == t1[12 * 8 + bit]);
        model.addConstr(t2[11 * 8 + bit] == t1[6 * 8 + bit]);

        model.addConstr(t2[12 * 8 + bit] == t1[7 * 8 + bit]);
        model.addConstr(t2[13 * 8 + bit] == t1[13 * 8 + bit]);
        model.addConstr(t2[14 * 8 + bit] == t1[2 * 8 + bit]);
        model.addConstr(t2[15 * 8 + bit] == t1[8 * 8 + bit]);

    }

}


//Skinny MixColumn
void MC(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 4; col++)
    {
        for (int bit = 0; bit < 8; bit++)
        {
            XOR3(model, s[(4 * col + 1) * 8 + bit], s[(4 * col + 2) * 8 + bit], s[(4 * col + 3) * 8 + bit], t[(4 * col + 0) * 8 + bit]);
            XOR3(model, s[(4 * col + 0) * 8 + bit], s[(4 * col + 2) * 8 + bit], s[(4 * col + 3) * 8 + bit], t[(4 * col + 1) * 8 + bit]);
            XOR3(model, s[(4 * col + 0) * 8 + bit], s[(4 * col + 1) * 8 + bit], s[(4 * col + 3) * 8 + bit], t[(4 * col + 2) * 8 + bit]);
            XOR3(model, s[(4 * col + 0) * 8 + bit], s[(4 * col + 1) * 8 + bit], s[(4 * col + 2) * 8 + bit], t[(4 * col + 3) * 8 + bit]);
        }
    }

}

double Search_Midori_trail(int round, bitset<300> Template, double min)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    env.set(GRB_DoubleParam_TimeLimit, 4000);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);
    GRBVar** s = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t1 = new GRBVar * [round + 1];        //State after the S-box
    GRBVar** t2 = new GRBVar * [round + 1];        //State after the SR
    for (i = 0; i < round + 1; i++)
    {
        s[i] = model.addVars(16 * 8, GRB_BINARY);
        t1[i] = model.addVars(16 * 8, GRB_BINARY);
        t2[i] = model.addVars(16 * 8, GRB_BINARY);
    }
    GRBVar*** p;
    p = new GRBVar * *[round + 2];        
    for (i = 1; i < round + 1; i++)
    {
        p[i] = new GRBVar * [16];
        for (j = 0; j < 16; j++)
        {
            p[i][j] = model.addVars(5, GRB_BINARY);
        }
    }
    p[0] = new GRBVar * [16];
    p[round + 1] = new GRBVar * [16];
    for (j = 0; j < 16; j++)
    {
        p[0][j] = model.addVars(3, GRB_BINARY);
        p[round + 1][j] = model.addVars(3, GRB_BINARY);
    }

    cond_SC(model, t1[0], p[0]);        //iCLAT
    SR(model, t1[0], t2[0]);
    MC(model, t2[0], s[0]);
    for (int loc = 0; loc < round; loc++)
    {
        SC(model, s[loc], t1[loc + 1], p[loc + 1]);
        SR(model, t1[loc + 1], t2[loc + 1]);
        MC(model, t2[loc + 1], s[loc + 1]);
    }
    cond_SC(model, s[round], p[round + 1]);     //oCLAT

    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 16 * 8; i++)
        init_mask_sum += t1[0][i];
    model.addConstr(init_mask_sum >= 1);

    int sbox_index[18][16] = { 0 };

    for (i = 0; i < round + 2; i++)
    {
        for (j = 0; j < 16; j++)
        {
            sbox_index[i][j] = Template[16 * i + j];
        }

    }

    for (int cell = 0; cell < 16; cell++)
    {
        if (sbox_index[0][cell] == 0)
        {
            model.addConstr(p[0][cell][0] == 1);
        }
        else
        {
            GRBLinExpr sum = 0;
            for (int j = 1; j < 3; j++)
                sum += p[0][cell][j];
            model.addConstr(sum == 1);
        }
    }
    for (int i = 1; i < round + 1; i++)
    {
        for (int cell = 0; cell < 16; cell++)
        {
            if (sbox_index[i][cell] == 0)
            {
                model.addConstr(p[i][cell][0] == 1);
            }
            else
            {
                GRBLinExpr sum = 0;
                for (int j = 1; j < 5; j++)
                    sum += p[i][cell][j];
                model.addConstr(sum == 1);
            }
        }
    }
    for (int cell = 0; cell < 16; cell++)
    {
        if (sbox_index[round + 1][cell] == 0)
        {
            model.addConstr(p[round + 1][cell][0] == 1);
        }
        else
        {
            GRBLinExpr sum = 0;
            for (int j = 1; j < 3; j++)
                sum += p[round + 1][cell][j];
            model.addConstr(sum == 1);
        }
    }


    GRBLinExpr cor = 0;
    for (i = 1; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            //2  4  6  8
            cor += 2 * p[i][j][1] + 4 * p[i][j][2] + 6 * p[i][j][3] + 8 * p[i][j][4];
        }
    }
    for (j = 0; j < 16; j++)
    {
        //1  2
        cor += p[0][j][1] + 2 * p[0][j][2];

        //1  2
        cor += p[round + 1][j][1] + 2 * p[round + 1][j][2];
    }

    model.addConstr(cor <= min - 0.01);

    model.setObjective(cor, GRB_MINIMIZE);


    model.update();
    model.optimize();


    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        std::cout << "===== Midori " << round + 2 << "-round result of the linear distinguisher =====" << std::endl;
        std::cout << "The S-box output mask for round 1: " << std::endl;
        for (int n1 = 0; n1 < 16; n1++)
        {
            int result = 0;
            for (int n2 = 0; n2 < 8; n2++)
            {
                if (t1[0][8 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
                    result = result * 2 + 1;
                else
                    result = result * 2;
            }
            std::cout << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
        }
        cout << endl;
        std::cout << "The S-box input mask for rounds 2 to " << round + 2 << ": " << std::endl;
        for (int i = 0; i < round + 1; i++)
        {
            for (int n1 = 0; n1 < 16; n1++)
            {
                int result = 0;
                for (int n2 = 0; n2 < 8; n2++)
                {
                    if (s[i][8 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
                        result = result * 2 + 1;
                    else
                        result = result * 2;
                }
                std::cout << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
            }
            cout << endl;
        }
        std::cout << "\nConditional linear weight: " << cor.getValue() << std::endl;
        return cor.getValue();
    }


    return min;
}






int main()
{
    int round;
    cout << " ROUND: ";
    cin >> round;




    vector<bitset<300>> foundSolutions = Search_Midori_Sbox(round);
    double min = 200;
    for (int i = 0; i < foundSolutions.size(); i++)
    {
        cout << "-------------------------- Try " << i << " --------------------------" << endl;
        min = Search_Midori_trail(round - 2, foundSolutions[i], min);
    }

    return 0;
}
