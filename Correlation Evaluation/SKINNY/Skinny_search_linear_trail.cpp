#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <map>
#include <bitset>
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
void OR5(GRBModel& model, GRBVar& in0, GRBVar& in1, GRBVar& in2, GRBVar& in3, GRBVar& out)
{
    model.addConstr(2 * out >= in0 + in1);
    model.addConstr(2 * out >= in2 + in3);
    model.addConstr(out <= in0 + in1 + in2 + in3);
}

//Skinny SubCell_Sbox
void SC_Sbox(GRBModel& model, GRBVar* s, GRBVar* t)
{
    for (int i = 0; i < 16; i++)
    {
        model.addConstr(s[i] == t[i]);
    }
}

//Skinny ShiftRow_Sbox
void SR_Sbox(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    model.addConstr(t2[0] == t1[0]);
    model.addConstr(t2[1] == t1[13]);
    model.addConstr(t2[2] == t1[10]);
    model.addConstr(t2[3] == t1[7]);

    model.addConstr(t2[4] == t1[4]);
    model.addConstr(t2[5] == t1[1]);
    model.addConstr(t2[6] == t1[14]);
    model.addConstr(t2[7] == t1[11]);

    model.addConstr(t2[8] == t1[8]);
    model.addConstr(t2[9] == t1[5]);
    model.addConstr(t2[10] == t1[2]);
    model.addConstr(t2[11] == t1[15]);

    model.addConstr(t2[12] == t1[12]);
    model.addConstr(t2[13] == t1[9]);
    model.addConstr(t2[14] == t1[6]);
    model.addConstr(t2[15] == t1[3]);


}


//Skinny MixColumn_Sbox
void MC_Sbox(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 4; col++)
    {
        GRBVar d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        XOR2(model, s[4 * col + 0], s[4 * col + 3], d);
        XOR2(model, s[4 * col + 1], d, t[4 * col + 0]);
        model.addConstr(s[4 * col + 2] == t[4 * col + 1]);
        XOR2(model, s[4 * col + 2], d, t[4 * col + 2]);
        model.addConstr(s[4 * col + 0] == t[4 * col + 3]);
    }

}

bitset<300> Search_Skinny_Sbox(int round, vector<bitset<300>> foundSolutions)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    //env.set(GRB_DoubleParam_TimeLimit, 100);


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


    for (int loc = 0; loc < foundSolutions.size(); loc++)
    {
        GRBLinExpr cut;
        for (int i = 0; i < round; i++)
        {
            for (int n1 = 0; n1 < 16; n1++)
            {
                if (foundSolutions[loc][16 * i + n1])
                {
                    cut += (1 - s[i][n1]);
                }
                else
                {
                    cut += s[i][n1];
                }
            }
        }
        model.addConstr(cut >= 1);
    }




    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 16; i++)
        init_mask_sum += s[0][i];
    model.addConstr(init_mask_sum >= 1);

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
        condition_cor += 1 * s[0][j];
        condition_cor += 1 * s[round - 1][j];
    }
    model.setObjective(condition_cor, GRB_MINIMIZE);



    model.update();
    model.optimize();
    bitset<300> temp(0);
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL)
    {
        for (int i = 0; i < round; i++)
        {
            for (int n1 = 0; n1 < 16; n1++)
            {
                if (s[i][n1].get(GRB_DoubleAttr_X) > 0.5)
                {
                    temp[16 * i + n1] = 1;
                }
                else
                {
                    temp[16 * i + n1] = 0;
                }
            }
        }
       
    }
    else 
    {
        cout << "No solution!!!!!" << endl;
    }

    return temp;
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
//Number of constraints : 46
void SBox_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a4 - b3 >= 0 - (1 - p));
    model.addConstr(a5 - a7 + b2 >= 0 - (1 - p));
    model.addConstr(-a2 + b0 + b2 >= 0 - (1 - p));
    model.addConstr(-a0 - a2 >= -1 - (1 - p));
    model.addConstr(-a3 - a4 >= -1 - (1 - p));
    model.addConstr(-a0 - a5 >= -1 - (1 - p));
    model.addConstr(-a2 - a5 >= -1 - (1 - p));
    model.addConstr(-a0 - a6 >= -1 - (1 - p));
    model.addConstr(-a2 - a7 >= -1 - (1 - p));
    model.addConstr(-a3 - b0 >= -1 - (1 - p));
    model.addConstr(-a5 - b1 >= -1 - (1 - p));
    model.addConstr(-a3 + b1 >= 0 - (1 - p));
    model.addConstr(-a0 - b2 >= -1 - (1 - p));
    model.addConstr(-a6 - b3 >= -1 - (1 - p));
    model.addConstr(-b4 - b5 >= -1 - (1 - p));
    model.addConstr(a5 - b7 >= 0 - (1 - p));
    model.addConstr(-a7 - b7 >= -1 - (1 - p));
    model.addConstr(-b5 - b7 >= -1 - (1 - p));
    model.addConstr(-a4 - a5 + a7 >= -1 - (1 - p));
    model.addConstr(a4 - a6 - b1 >= -1 - (1 - p));
    model.addConstr(-a6 + b0 - b1 >= -1 - (1 - p));
    model.addConstr(a6 - a7 + b2 >= 0 - (1 - p));
    model.addConstr(a2 - b2 - b3 >= -1 - (1 - p));
    model.addConstr(-a2 + b0 + b3 >= 0 - (1 - p));
    model.addConstr(-a5 - b2 - b4 >= -2 - (1 - p));
    model.addConstr(-a0 - a1 - b5 >= -2 - (1 - p));
    model.addConstr(-b1 + b3 - b5 >= -1 - (1 - p));
    model.addConstr(-a1 + a3 + b5 >= 0 - (1 - p));
    model.addConstr(-a6 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(-b1 + b5 - b6 >= -1 - (1 - p));
    model.addConstr(-a0 + a3 + b6 >= 0 - (1 - p));
    model.addConstr(a1 - b5 + b6 >= 0 - (1 - p));
    model.addConstr(-b0 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(a6 - b4 + b7 >= 0 - (1 - p));
    model.addConstr(a0 - b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a2 - a4 - a6 + b1 >= -2 - (1 - p));
    model.addConstr(a0 + a2 - b0 + b3 >= 0 - (1 - p));
    model.addConstr(b0 - b1 - b2 + b3 >= -1 - (1 - p));
    model.addConstr(a4 - a6 - a7 + b4 >= -1 - (1 - p));
    model.addConstr(-a4 - b0 - b2 + b4 >= -2 - (1 - p));
    model.addConstr(-a4 + b2 + b3 + b4 >= 0 - (1 - p));
    model.addConstr(-a4 + a6 + b1 - b5 >= -1 - (1 - p));
    model.addConstr(-a5 - a6 + a7 + b5 >= -1 - (1 - p));
    model.addConstr(-a4 + a6 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(a7 + b0 - b2 + b3 + b4 >= 0 - (1 - p));
    model.addConstr(a1 + a3 + a7 + b0 + b3 + b4 + b6 + b7 >= 1 - (1 - p));

}

//weight: 225/1024   2^-2.186
//Number of constraints : 16
void SBox_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a0 >= 1 - (1 - p));
    model.addConstr(-a2 >= 0 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));
    model.addConstr(-a4 >= 0 - (1 - p));
    model.addConstr(a5 >= 1 - (1 - p));
    model.addConstr(-a6 >= 0 - (1 - p));
    model.addConstr(-a7 >= 0 - (1 - p));
    model.addConstr(-b0 >= 0 - (1 - p));
    model.addConstr(b1 >= 1 - (1 - p));
    model.addConstr(-b2 >= 0 - (1 - p));
    model.addConstr(-b3 >= 0 - (1 - p));
    model.addConstr(-b5 >= 0 - (1 - p));
    model.addConstr(b6 >= 1 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(a1 - b4 >= 0 - (1 - p));
    model.addConstr(-a1 + b4 >= 0 - (1 - p));

}

//weight: 169/1024   2^-2.599
//Number of constraints : 16
void SBox_p3(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a0 >= 1 - (1 - p));
    model.addConstr(-a2 >= 0 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));
    model.addConstr(-a4 >= 0 - (1 - p));
    model.addConstr(a5 >= 1 - (1 - p));
    model.addConstr(-a6 >= 0 - (1 - p));
    model.addConstr(-a7 >= 0 - (1 - p));
    model.addConstr(-b0 >= 0 - (1 - p));
    model.addConstr(b1 >= 1 - (1 - p));
    model.addConstr(-b2 >= 0 - (1 - p));
    model.addConstr(-b3 >= 0 - (1 - p));
    model.addConstr(-b5 >= 0 - (1 - p));
    model.addConstr(-b6 >= 0 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(a1 - b4 >= 0 - (1 - p));
    model.addConstr(-a1 + b4 >= 0 - (1 - p));

}

//weight: 9/64   2^-2.83
//Number of constraints : 21
void SBox_p4(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(-a2 >= 0 - (1 - p));
    model.addConstr(-b3 >= 0 - (1 - p));
    model.addConstr(-a4 + a6 >= 0 - (1 - p));
    model.addConstr(-a0 - b0 >= -1 - (1 - p));
    model.addConstr(a7 - b2 >= 0 - (1 - p));
    model.addConstr(-a7 + b2 >= 0 - (1 - p));
    model.addConstr(-a5 - b6 >= -1 - (1 - p));
    model.addConstr(a5 + b6 >= 1 - (1 - p));
    model.addConstr(a3 - b1 >= 0 - (1 - p));
    model.addConstr(-a0 + b1 >= 0 - (1 - p));
    model.addConstr(-b0 + b1 >= 0 - (1 - p));
    model.addConstr(-a4 + b2 >= 0 - (1 - p));
    model.addConstr(-b1 - b4 >= -1 - (1 - p));
    model.addConstr(-a3 - b7 >= -1 - (1 - p));
    model.addConstr(-a6 + b7 >= 0 - (1 - p));
    model.addConstr(-b2 + b7 >= 0 - (1 - p));
    model.addConstr(b6 + b7 >= 1 - (1 - p));
    model.addConstr(a0 + a1 + b0 >= 1 - (1 - p));
    model.addConstr(a4 - a6 - b2 >= -1 - (1 - p));
    model.addConstr(a0 + b0 + b5 >= 1 - (1 - p));
    model.addConstr(a0 + b0 + b7 >= 1 - (1 - p));

}

//weight: 121/1024   2^-3.081
//Number of constraints : 17
void SBox_p5(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(-a2 >= 0 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));
    model.addConstr(-a4 >= 0 - (1 - p));
    model.addConstr(a5 >= 1 - (1 - p));
    model.addConstr(-a6 >= 0 - (1 - p));
    model.addConstr(-a7 >= 0 - (1 - p));
    model.addConstr(-b0 >= 0 - (1 - p));
    model.addConstr(b1 >= 1 - (1 - p));
    model.addConstr(-b2 >= 0 - (1 - p));
    model.addConstr(-b3 >= 0 - (1 - p));
    model.addConstr(-b6 >= 0 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(a0 - b5 >= 0 - (1 - p));
    model.addConstr(-a0 - a1 - b4 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 + b4 >= 0 - (1 - p));
    model.addConstr(a1 - b4 + b5 >= 0 - (1 - p));
    model.addConstr(-a1 + b4 + b5 >= 0 - (1 - p));

}

//weight: 81/1024   2^-3.66
//Number of constraints : 28
void SBox_p6(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(-a2 >= 0 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));
    model.addConstr(-b0 >= 0 - (1 - p));
    model.addConstr(b1 >= 1 - (1 - p));
    model.addConstr(-b3 >= 0 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(-a1 - a6 - b4 >= -2 - (1 - p));
    model.addConstr(a1 - a6 + b4 >= 0 - (1 - p));
    model.addConstr(-a4 + a6 >= 0 - (1 - p));
    model.addConstr(a7 - b2 >= 0 - (1 - p));
    model.addConstr(-a4 + b2 >= 0 - (1 - p));
    model.addConstr(a5 + b2 >= 1 - (1 - p));
    model.addConstr(a0 + a4 - a7 >= 0 - (1 - p));
    model.addConstr(-a0 - a6 - b5 >= -2 - (1 - p));
    model.addConstr(a0 - b5 - b6 >= -1 - (1 - p));
    model.addConstr(a5 + b5 - b6 >= 0 - (1 - p));
    model.addConstr(-a5 + b5 + b6 >= 0 - (1 - p));
    model.addConstr(-a0 - a1 - b4 - b5 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 + b4 - b5 >= -1 - (1 - p));
    model.addConstr(-a0 + a6 + b2 + b5 >= 0 - (1 - p));
    model.addConstr(a0 + a1 - b4 + b5 >= 0 - (1 - p));
    model.addConstr(a0 - a1 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(a1 + a5 - b4 + b6 >= 0 - (1 - p));
    model.addConstr(-a1 + a5 + b4 + b6 >= 0 - (1 - p));
    model.addConstr(a1 + a6 + b2 - b4 + b6 >= 0 - (1 - p));
    model.addConstr(-a1 + a6 + b2 + b4 + b6 >= 0 - (1 - p));
    model.addConstr(a1 + a4 - a5 - a7 - b4 - b6 >= -3 - (1 - p));
    model.addConstr(-a1 + a4 - a5 - a7 + b4 - b6 >= -3 - (1 - p));

}

//weight: 1/16   2^-4
//Number of constraints : 101
void SBox_p7(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(-a3 - a6 - b3 >= -2 - (1 - p));
    model.addConstr(-a3 - a7 - b3 >= -2 - (1 - p));
    model.addConstr(-a3 - b4 - b5 >= -2 - (1 - p));
    model.addConstr(-a0 - a7 - b7 >= -2 - (1 - p));
    model.addConstr(-a6 - b3 - b7 >= -2 - (1 - p));
    model.addConstr(-a7 - b3 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a5 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 + b0 + b2 >= -1 - (1 - p));
    model.addConstr(-a1 + a4 + a6 - b3 >= -1 - (1 - p));
    model.addConstr(a0 - a2 + b0 + b3 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 + b0 + b3 >= -1 - (1 - p));
    model.addConstr(a2 + a3 - b1 + b3 >= 0 - (1 - p));
    model.addConstr(a1 - a3 + a6 - b4 >= -1 - (1 - p));
    model.addConstr(a4 + a6 - b3 - b4 >= -1 - (1 - p));
    model.addConstr(a4 + a7 - b3 + b4 >= 0 - (1 - p));
    model.addConstr(a0 + a1 + a3 - b5 >= 0 - (1 - p));
    model.addConstr(-a1 - b3 - b4 - b5 >= -3 - (1 - p));
    model.addConstr(a0 - a1 + a3 + b5 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 + a5 - b6 >= -2 - (1 - p));
    model.addConstr(-a3 - a4 + b3 - b6 >= -2 - (1 - p));
    model.addConstr(-a3 + b1 + b3 - b6 >= -1 - (1 - p));
    model.addConstr(-a2 - a5 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(-a2 - a3 - a5 + b6 >= -2 - (1 - p));
    model.addConstr(-a3 - a7 - b0 + b6 >= -2 - (1 - p));
    model.addConstr(-a0 - a2 + a3 - b7 >= -2 - (1 - p));
    model.addConstr(a1 + a5 + a7 - b7 >= 0 - (1 - p));
    model.addConstr(-a1 + a2 + b1 - b7 >= -1 - (1 - p));
    model.addConstr(a2 - a3 + b1 - b7 >= -1 - (1 - p));
    model.addConstr(-a0 - a6 + b1 - b7 >= -2 - (1 - p));
    model.addConstr(-a7 - b0 + b1 - b7 >= -2 - (1 - p));
    model.addConstr(a2 - b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a7 - b0 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a1 + a5 + b6 - b7 >= -1 - (1 - p));
    model.addConstr(a6 - b4 + b5 + b7 >= 0 - (1 - p));
    model.addConstr(-a1 - a2 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a0 + a3 - b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a3 - a6 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a2 - a7 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a6 - b4 - b6 + b7 >= -1 - (1 - p));
    model.addConstr(-a0 + a3 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a0 + a2 - a3 - b0 + b3 >= -2 - (1 - p));
    model.addConstr(a0 + a2 - b0 + b1 + b3 >= 0 - (1 - p));
    model.addConstr(-a2 - a5 + b0 + b2 - b4 >= -2 - (1 - p));
    model.addConstr(-a5 - b0 - b2 - b3 - b4 >= -4 - (1 - p));
    model.addConstr(-a2 - a7 + b0 - b2 + b4 >= -2 - (1 - p));
    model.addConstr(-a4 - a6 - a7 - b3 + b4 >= -3 - (1 - p));
    model.addConstr(-a7 - b0 + b2 - b3 + b4 >= -2 - (1 - p));
    model.addConstr(-a3 + a7 - b0 + b3 - b5 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - b1 - b6 >= -4 - (1 - p));
    model.addConstr(-a0 - a1 - a6 - b4 - b6 >= -4 - (1 - p));
    model.addConstr(-a0 + a1 - a6 + b4 - b6 >= -2 - (1 - p));
    model.addConstr(-a3 - b1 + b2 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(-a1 - b1 - b3 - b5 - b6 >= -4 - (1 - p));
    model.addConstr(-a3 - a4 + a7 - b5 + b6 >= -2 - (1 - p));
    model.addConstr(a3 + a4 - a6 - b1 - b7 >= -2 - (1 - p));
    model.addConstr(-a1 + a3 + a7 - b1 - b7 >= -2 - (1 - p));
    model.addConstr(a2 + a6 + a7 - b2 - b7 >= -1 - (1 - p));
    model.addConstr(-a4 + a6 + a7 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a1 - a4 + a7 + b4 - b7 >= -2 - (1 - p));
    model.addConstr(a2 + a7 - b2 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(a2 - a7 + b2 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(-a4 + a7 + b3 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(a1 - a5 + a7 - b6 + b7 >= -1 - (1 - p));
    model.addConstr(a0 + b0 + b3 - b6 + b7 >= 0 - (1 - p));
    model.addConstr(a2 - b0 + b3 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a5 + a7 + b5 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a6 + b4 + b5 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(a1 + a2 + a3 + a5 + a6 + b2 >= 1 - (1 - p));
    model.addConstr(-a2 + a3 - a4 + a7 + b1 + b3 >= -1 - (1 - p));
    model.addConstr(a4 + a6 - a7 - b0 + b1 - b4 >= -2 - (1 - p));
    model.addConstr(a0 - a2 + a7 + b0 + b2 + b4 >= 0 - (1 - p));
    model.addConstr(-a1 - a4 + a6 + b1 + b3 + b4 >= -1 - (1 - p));
    model.addConstr(-a4 + a6 - a7 - b1 + b3 - b5 >= -3 - (1 - p));
    model.addConstr(a0 + a2 + a6 + b2 + b3 + b5 >= 1 - (1 - p));
    model.addConstr(a2 - a3 + b2 + b3 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(-a2 + a4 + b1 - b3 + b5 - b6 >= -2 - (1 - p));
    model.addConstr(a7 - b0 - b2 - b3 + b4 + b6 >= -2 - (1 - p));
    model.addConstr(a1 + a2 + a5 - a7 + b2 + b7 >= 0 - (1 - p));
    model.addConstr(a1 + a2 - a5 - b2 - b4 + b7 >= -2 - (1 - p));
    model.addConstr(a4 - a6 - a7 + b3 + b4 + b7 >= -1 - (1 - p));
    model.addConstr(a2 + a5 - a7 + b2 + b5 + b7 >= 0 - (1 - p));
    model.addConstr(a2 - a5 - b2 - b4 + b5 + b7 >= -2 - (1 - p));
    model.addConstr(-a1 + b1 - b3 + b5 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a2 + b1 + b3 + b4 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(a2 + b1 + b3 + b5 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(a3 + a7 + b4 + b5 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(-a0 + a1 + a2 + a3 + b1 - b3 - b5 >= -2 - (1 - p));
    model.addConstr(-a4 - a5 + a6 + a7 - b0 - b2 + b5 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 + a2 + a3 - b1 - b3 + b5 >= -2 - (1 - p));
    model.addConstr(-a1 - a3 - a5 + a6 + a7 + b4 - b6 >= -3 - (1 - p));
    model.addConstr(-a4 + a5 - b1 - b3 - b4 + b5 - b6 >= -4 - (1 - p));
    model.addConstr(a0 + a2 + a5 - a6 + b1 - b2 + b3 + b5 >= -1 - (1 - p));
    model.addConstr(a0 + a2 - a5 - a7 + b1 + b2 + b3 + b5 >= -1 - (1 - p));
    model.addConstr(-a2 + a3 + a4 + a7 - b1 + b3 - b4 + b7 >= -2 - (1 - p));
    model.addConstr(-a2 + a3 + a4 - a6 - b1 - b4 - b5 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 + a3 - a4 - a6 + b1 - b4 - b5 + b7 >= -4 - (1 - p));
    model.addConstr(a2 + b1 + b2 - b3 + b4 + b5 - b6 + b7 >= -1 - (1 - p));
    model.addConstr(-a4 + a7 - b0 - b1 + b3 + b5 + b6 + b7 >= -2 - (1 - p));
    model.addConstr(a4 + a7 + b1 + b3 - b4 + b5 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(a5 - a7 - b0 + b2 - b3 >= -2 - (1 - p));
    model.addConstr(-a2 + a5 - a7 + b0 - b2 >= -2 - (1 - p));

}

//weight: 49/1024   2^-4.385
//Number of constraints : 59
void SBox_p8(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(-a7 + b1 >= 0 - (1 - p));
    model.addConstr(a4 - b3 >= 0 - (1 - p));
    model.addConstr(-a6 - b3 >= -1 - (1 - p));
    model.addConstr(-a7 - b3 >= -1 - (1 - p));
    model.addConstr(a4 - a6 - b0 >= -1 - (1 - p));
    model.addConstr(-b0 - b1 - b2 >= -2 - (1 - p));
    model.addConstr(a0 - a7 + b2 >= 0 - (1 - p));
    model.addConstr(a0 + b1 + b3 >= 1 - (1 - p));
    model.addConstr(a2 + b1 + b3 >= 1 - (1 - p));
    model.addConstr(b0 + b1 + b3 >= 1 - (1 - p));
    model.addConstr(-a7 + b2 - b5 >= -1 - (1 - p));
    model.addConstr(a5 + a7 + b5 >= 1 - (1 - p));
    model.addConstr(a2 - a4 - b0 - b1 >= -2 - (1 - p));
    model.addConstr(-a2 - a4 + a6 + b1 >= -1 - (1 - p));
    model.addConstr(a4 - a6 - b2 - b5 >= -2 - (1 - p));
    model.addConstr(a0 + a4 - b2 + b5 >= 0 - (1 - p));
    model.addConstr(-a1 - a6 - b4 + b5 >= -2 - (1 - p));
    model.addConstr(a1 - a6 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(-a1 - a6 - b4 - b6 >= -3 - (1 - p));
    model.addConstr(a1 - a6 + b4 - b6 >= -1 - (1 - p));
    model.addConstr(a0 - a6 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(a5 - a6 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(a5 + b2 - b5 + b6 >= 0 - (1 - p));
    model.addConstr(-a1 - a4 + a7 + b3 - b4 >= -2 - (1 - p));
    model.addConstr(a1 - a4 + a7 + b3 + b4 >= 0 - (1 - p));
    model.addConstr(a0 - a5 + b0 - b5 - b6 >= -2 - (1 - p));
    model.addConstr(a2 + a6 + b2 + b3 + b6 >= 1 - (1 - p));
    model.addConstr(-a1 + a4 + a5 - b4 + b6 >= -1 - (1 - p));
    model.addConstr(a0 + a1 + a6 - b4 + b6 >= 0 - (1 - p));
    model.addConstr(a1 + a4 + a5 + b4 + b6 >= 1 - (1 - p));
    model.addConstr(a0 - a1 + a6 + b4 + b6 >= 0 - (1 - p));
    model.addConstr(a0 - a4 + a7 + b5 + b6 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 - b1 + b3 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a5 + a6 + a7 + b3 - b5 >= -2 - (1 - p));
    model.addConstr(-a0 + a4 - a5 - b1 + b2 + b5 >= -2 - (1 - p));
    model.addConstr(a0 - a1 - a4 + b2 - b4 + b5 >= -2 - (1 - p));
    model.addConstr(a0 + a1 - a4 + b2 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(a1 + a4 + a5 - b2 - b4 - b6 >= -2 - (1 - p));
    model.addConstr(-a1 + a4 + a5 - b2 + b4 - b6 >= -2 - (1 - p));
    model.addConstr(a4 - a6 + a7 + b2 + b5 - b6 >= -1 - (1 - p));
    model.addConstr(a1 - a5 + b3 - b4 - b5 + b6 >= -2 - (1 - p));
    model.addConstr(-a1 - a5 + b3 + b4 - b5 + b6 >= -2 - (1 - p));
    model.addConstr(a2 - a3 + a6 + a7 - b2 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 - a4 + a6 - b1 - b5 - b7 >= -4 - (1 - p));
    model.addConstr(-a3 - a4 + a7 - b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(a4 - a5 + a6 - b1 - b2 + b5 - b6 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a3 + a6 - b4 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 + a1 - a3 + a6 + b4 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 - a3 + a6 - b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a3 + a6 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a3 - a4 - b1 + b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a1 - a3 + a6 + b0 - b1 + b2 - b4 - b7 >= -4 - (1 - p));
    model.addConstr(a1 - a3 + a6 + b0 - b1 + b2 + b4 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 + a2 - a3 - a4 + a6 + b5 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a2 - a3 + a5 - a6 - a7 + b0 - b1 - b2 + b3 + b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 + a2 - a3 - a4 - a5 - a6 - a7 + b0 - b1 + b3 + b5 - b7 >= -7 - (1 - p));
    model.addConstr(a5 - b0 >= 0 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));

}

//weight: 9/256   2^-4.83
//Number of constraints : 63
void SBox_p9(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a3 + b7 >= 1 - (1 - p));
    model.addConstr(b6 + b7 >= 1 - (1 - p));
    model.addConstr(-a2 - a3 - a7 >= -2 - (1 - p));
    model.addConstr(a2 - a3 + b1 >= 0 - (1 - p));
    model.addConstr(a4 + a7 - b3 >= 0 - (1 - p));
    model.addConstr(a0 + a3 + b3 >= 1 - (1 - p));
    model.addConstr(a2 - a3 + b4 >= 0 - (1 - p));
    model.addConstr(a0 + a3 + b5 >= 1 - (1 - p));
    model.addConstr(-a2 - a5 - b6 >= -2 - (1 - p));
    model.addConstr(-a2 + a5 + b6 >= 0 - (1 - p));
    model.addConstr(-a5 + a7 + b7 >= 0 - (1 - p));
    model.addConstr(a0 + a2 - a3 + b0 >= 0 - (1 - p));
    model.addConstr(a2 + a5 - a7 + b2 >= 0 - (1 - p));
    model.addConstr(-a4 - a6 - a7 - b3 >= -3 - (1 - p));
    model.addConstr(-a0 - a2 + b0 + b3 >= -1 - (1 - p));
    model.addConstr(a3 + b1 - b3 + b5 >= 0 - (1 - p));
    model.addConstr(-a3 + a6 - b4 - b6 >= -2 - (1 - p));
    model.addConstr(-a3 - b3 - b4 - b6 >= -3 - (1 - p));
    model.addConstr(-a3 - a6 + b4 - b6 >= -2 - (1 - p));
    model.addConstr(a2 + a7 - b2 - b7 >= -1 - (1 - p));
    model.addConstr(a2 - a7 + b2 - b7 >= -1 - (1 - p));
    model.addConstr(-a7 - b0 + b2 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 - b1 - b3 - b7 >= -3 - (1 - p));
    model.addConstr(a2 - b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a3 + b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a4 + b2 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(a0 - a5 - b6 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 + a5 - b6 - b7 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a6 - b1 >= -4 - (1 - p));
    model.addConstr(-a0 + a2 - b0 - b1 + b3 >= -2 - (1 - p));
    model.addConstr(-a2 + a4 + a7 + b1 - b4 >= -1 - (1 - p));
    model.addConstr(-a0 - a3 + b0 + b2 + b4 >= -1 - (1 - p));
    model.addConstr(a0 - a2 - b0 + b3 - b6 >= -2 - (1 - p));
    model.addConstr(-a5 - b1 - b2 + b3 - b6 >= -3 - (1 - p));
    model.addConstr(-a4 + a7 - b0 - b2 - b7 >= -3 - (1 - p));
    model.addConstr(-a2 + a7 + b0 + b2 - b7 >= -1 - (1 - p));
    model.addConstr(a7 + b0 - b2 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a0 + a2 - b1 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a3 - a4 - a5 - a6 + a7 >= -4 - (1 - p));
    model.addConstr(-a0 - a4 + a7 - b0 - b2 + b4 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a6 - b1 + b4 - b5 >= -4 - (1 - p));
    model.addConstr(a0 - a3 - b0 + b2 + b4 - b6 >= -2 - (1 - p));
    model.addConstr(a2 + a4 - a6 - b2 + b3 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a1 - a3 + a6 + b4 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a3 + a6 - b3 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 + a6 + b4 - b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a4 - a6 + a7 + b1 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a1 - a3 + a6 + b4 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a1 + a3 - a5 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a1 - a3 - a5 - a6 + a7 - b4 >= -4 - (1 - p));
    model.addConstr(-a0 + a1 - a5 - a6 - b0 - b4 - b5 >= -5 - (1 - p));
    model.addConstr(a0 - a3 - a4 + b0 - b2 + b4 - b6 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a5 + a6 - b4 - b5 + b6 >= -4 - (1 - p));
    model.addConstr(-a1 - a6 - b1 + b3 - b4 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a1 - a6 + b4 >= 0 - (1 - p));
    model.addConstr(a0 + a1 + a6 - b4 >= 0 - (1 - p));
    model.addConstr(a1 - a6 + b4 + b5 >= 0 - (1 - p));
    model.addConstr(a1 + a6 - b4 + b5 >= 0 - (1 - p));
    model.addConstr(-a4 + a6 + b3 >= 0 - (1 - p));
    model.addConstr(-a2 - a7 + b0 - b2 >= -2 - (1 - p));
    model.addConstr(a4 + a6 - b3 >= 0 - (1 - p));
    model.addConstr(a1 + a3 >= 1 - (1 - p));
    model.addConstr(a0 + a5 + b6 >= 1 - (1 - p));

}

//weight: 25/1024   2^-5.356
//Number of constraints : 128
void SBox_p10(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a5 + a7 - b0 >= 0 - (1 - p));
    model.addConstr(a4 + a6 - b3 >= 0 - (1 - p));
    model.addConstr(a4 + a7 - b3 >= 0 - (1 - p));
    model.addConstr(a5 + a7 + b5 >= 1 - (1 - p));
    model.addConstr(a2 - a5 - a7 + b2 >= -1 - (1 - p));
    model.addConstr(a2 + a6 - a7 + b2 >= 0 - (1 - p));
    model.addConstr(a2 + a5 - b0 + b2 >= 0 - (1 - p));
    model.addConstr(-a2 + a4 - a6 + b3 >= -1 - (1 - p));
    model.addConstr(-a4 + a6 - b0 + b3 >= -1 - (1 - p));
    model.addConstr(a5 + b1 - b5 + b6 >= 0 - (1 - p));
    model.addConstr(a5 - b3 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(-a2 - a4 - a6 - a7 - b1 >= -4 - (1 - p));
    model.addConstr(-a0 + a2 + a4 + a5 - b2 >= -1 - (1 - p));
    model.addConstr(-a4 - a6 + b1 - b2 - b3 >= -3 - (1 - p));
    model.addConstr(-a4 - a7 - b0 + b2 - b3 >= -3 - (1 - p));
    model.addConstr(-a4 + a7 - b0 - b1 + b3 >= -2 - (1 - p));
    model.addConstr(a4 - a6 - b0 - b2 + b3 >= -2 - (1 - p));
    model.addConstr(a0 + a1 - a2 + a6 - b4 >= -1 - (1 - p));
    model.addConstr(-a1 - a2 - a6 - b1 - b4 >= -4 - (1 - p));
    model.addConstr(a0 - a1 - a2 + a6 + b4 >= -1 - (1 - p));
    model.addConstr(a1 - a2 - a6 - b1 + b4 >= -2 - (1 - p));
    model.addConstr(a2 + a5 - a6 - b2 + b5 >= -1 - (1 - p));
    model.addConstr(a5 - b0 - b1 + b3 - b6 >= -2 - (1 - p));
    model.addConstr(-a5 - a7 + b1 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(-a5 - a7 - b3 - b5 - b6 >= -4 - (1 - p));
    model.addConstr(-a6 + b1 - b3 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(a2 + a6 + b2 + b5 - b6 >= 0 - (1 - p));
    model.addConstr(a2 - a3 - b1 - b3 - b7 >= -3 - (1 - p));
    model.addConstr(-a2 - a3 + b1 - b3 - b7 >= -3 - (1 - p));
    model.addConstr(-a2 - a3 + b0 + b3 - b7 >= -2 - (1 - p));
    model.addConstr(a2 - a3 + b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(a0 - a1 - a7 + b0 + b2 - b4 >= -2 - (1 - p));
    model.addConstr(a0 - a1 - a4 + a7 + b3 - b4 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 - a4 + a7 + b3 - b4 >= -2 - (1 - p));
    model.addConstr(a0 + a1 - a7 + b0 + b2 + b4 >= 0 - (1 - p));
    model.addConstr(-a0 - a1 - a4 + a7 + b3 + b4 >= -2 - (1 - p));
    model.addConstr(a0 + a1 - a4 + a7 + b3 + b4 >= 0 - (1 - p));
    model.addConstr(a0 - a7 - b1 + b2 + b3 - b5 >= -2 - (1 - p));
    model.addConstr(-a1 - a2 - a7 - b1 - b4 - b5 >= -5 - (1 - p));
    model.addConstr(a0 + a1 + a6 - b3 - b4 - b5 >= -2 - (1 - p));
    model.addConstr(a1 - a2 - a7 - b1 + b4 - b5 >= -3 - (1 - p));
    model.addConstr(a0 - a1 + a6 - b3 + b4 - b5 >= -2 - (1 - p));
    model.addConstr(a6 + a7 - b0 - b1 - b2 + b5 >= -2 - (1 - p));
    model.addConstr(-a0 + a7 - b0 - b2 - b3 - b6 >= -4 - (1 - p));
    model.addConstr(a0 - a5 + b0 - b2 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(a1 + a5 - b2 - b4 + b5 - b6 >= -2 - (1 - p));
    model.addConstr(-a1 + a5 - b2 + b4 + b5 - b6 >= -2 - (1 - p));
    model.addConstr(a6 + a7 - b0 - b1 - b2 + b6 >= -2 - (1 - p));
    model.addConstr(-a1 + a5 - a6 - b2 - b4 + b6 >= -3 - (1 - p));
    model.addConstr(a1 + a5 - a6 - b2 + b4 + b6 >= -1 - (1 - p));
    model.addConstr(a4 + a5 + b0 - b2 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(a2 - a4 + a6 + b2 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(-a0 + a4 - a5 - b2 + b5 + b6 >= -2 - (1 - p));
    model.addConstr(a2 - a3 + a6 + a7 - b2 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a3 - b0 - b1 + b3 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a7 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a4 + a6 + a7 + b0 - b1 + b2 >= -2 - (1 - p));
    model.addConstr(a0 - a1 + a4 - a6 - b2 + b3 - b4 >= -3 - (1 - p));
    model.addConstr(a0 + a1 + a4 - a6 - b2 + b3 + b4 >= -1 - (1 - p));
    model.addConstr(-a1 + a4 + a6 + b0 - b2 - b4 - b5 >= -3 - (1 - p));
    model.addConstr(a1 + a4 + a6 + b0 - b2 + b4 - b5 >= -1 - (1 - p));
    model.addConstr(-a4 + a6 + a7 + b0 - b1 + b2 + b5 >= -1 - (1 - p));
    model.addConstr(a0 - a1 - a4 - a7 - b3 - b4 + b5 >= -4 - (1 - p));
    model.addConstr(a0 + a1 - a6 + b1 - b3 - b4 + b5 >= -2 - (1 - p));
    model.addConstr(a0 + a1 - a4 - a7 - b3 + b4 + b5 >= -2 - (1 - p));
    model.addConstr(a0 - a1 - a6 + b1 - b3 + b4 + b5 >= -2 - (1 - p));
    model.addConstr(a0 - a2 - a6 - b0 - b1 + b2 - b6 >= -4 - (1 - p));
    model.addConstr(a0 - a6 + b0 - b1 - b2 - b3 - b6 >= -4 - (1 - p));
    model.addConstr(-a0 + a7 + b0 - b1 + b2 - b3 - b6 >= -3 - (1 - p));
    model.addConstr(a0 + a2 - a4 + a6 - b1 - b5 - b6 >= -3 - (1 - p));
    model.addConstr(-a0 - a7 + b0 - b1 - b2 + b5 - b6 >= -4 - (1 - p));
    model.addConstr(a0 - a4 - a5 - a6 - a7 - b1 + b6 >= -4 - (1 - p));
    model.addConstr(-a2 + a4 + a5 - b0 - b1 - b2 + b6 >= -3 - (1 - p));
    model.addConstr(a0 - a5 + a6 - b0 - b1 - b2 + b6 >= -3 - (1 - p));
    model.addConstr(a0 - a4 + b0 - b1 + b2 - b3 + b6 >= -2 - (1 - p));
    model.addConstr(a0 + a4 - a5 - a6 + b2 + b3 + b6 >= -1 - (1 - p));
    model.addConstr(-a1 + a4 + a5 + b0 - b4 + b5 + b6 >= -1 - (1 - p));
    model.addConstr(a1 - a5 + a6 - b2 - b4 + b5 + b6 >= -2 - (1 - p));
    model.addConstr(a1 + a4 + a5 + b0 + b4 + b5 + b6 >= 1 - (1 - p));
    model.addConstr(-a1 - a5 + a6 - b2 + b4 + b5 + b6 >= -2 - (1 - p));
    model.addConstr(-a0 - a1 - a3 - a6 - b4 + b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 + a1 - a3 + a6 - b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 - a3 - a6 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a3 + a6 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a3 + a5 - a7 + b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a3 - a5 - a7 - b3 + b5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 + a1 + a4 - a5 + b0 - b1 - b2 - b4 >= -4 - (1 - p));
    model.addConstr(-a0 - a1 + a4 - a5 + b0 - b1 - b2 + b4 >= -4 - (1 - p));
    model.addConstr(-a1 - a4 - a5 - a6 - a7 + b0 - b4 - b5 >= -6 - (1 - p));
    model.addConstr(a1 - a4 - a5 - a6 - a7 + b0 + b4 - b5 >= -4 - (1 - p));
    model.addConstr(-a0 + a6 - a7 + b0 - b1 - b2 - b5 + b6 >= -4 - (1 - p));
    model.addConstr(-a0 - a4 + b0 - b1 - b2 - b3 + b5 + b6 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 + a4 - a5 - b0 - b1 + b2 - b7 >= -5 - (1 - p));
    model.addConstr(a0 + a1 - a3 - a4 + a6 - b1 - b4 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a1 - a3 - a4 + a6 - b1 + b4 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a2 - a3 + a7 - b0 + b3 - b5 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 + a1 - a2 - a3 - a6 - b4 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(a0 - a1 - a3 - a6 + b1 - b4 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a3 + a6 + b1 - b4 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a2 - a3 - a6 + b4 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(a0 + a1 - a3 - a6 + b1 + b4 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 - a3 + a6 + b1 + b4 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a3 - a4 - b0 - b1 - b2 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a2 - a3 - a5 + a6 + a7 - b4 - b7 >= -6 - (1 - p));
    model.addConstr(-a0 + a1 - a2 - a3 - a5 + a6 + a7 + b4 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a3 - a4 - a5 - a6 + a7 - b1 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(a0 - a3 + a4 - a6 + b0 - b1 - b2 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 + a1 - a3 - a5 - a6 + b1 - b4 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a3 - a5 - a6 + b1 + b4 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a2 - a3 - a5 - b0 + b2 + b5 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(a0 - a1 + a2 - a3 - a4 - b4 + b5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a1 + a2 - a3 - a4 + b4 + b5 + b6 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 - a3 + a5 - a6 - a7 + b0 + b2 - b4 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a3 + a5 - a6 - a7 + b0 + b2 + b4 - b7 >= -5 - (1 - p));
    model.addConstr(a1 - a3 - a4 - a7 + b0 - b1 - b2 - b4 - b6 - b7 >= -7 - (1 - p));
    model.addConstr(-a1 - a3 - a4 - a7 + b0 - b1 - b2 + b4 - b6 - b7 >= -7 - (1 - p));
    model.addConstr(-a3 - a5 + a6 + a7 + b0 - b1 + b2 - b5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a3 - a5 + a6 + a7 - b0 - b1 + b2 - b5 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(a1 + a2 - a3 + a4 + a7 + b0 - b1 + b2 + b3 - b4 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a1 + a2 - a3 + a4 + a7 + b0 - b1 + b2 + b3 + b4 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a1 + a2 - a3 + a4 + a7 + b0 - b1 + b2 + b3 - b4 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(a1 + a2 - a3 + a4 + a7 + b0 - b1 + b2 + b3 + b4 + b6 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a1 + a2 - a3 + a4 - a5 + a6 + b0 - b1 + b3 - b4 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(a0 + a1 + a2 - a3 + a4 - a5 + a6 + b0 - b1 + b3 + b4 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(a2 - a3 + a4 - a5 - a6 + a7 + b0 - b1 + b2 + b3 - b5 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(-a4 + a5 + a7 >= 0 - (1 - p));
    model.addConstr(b7 >= 1 - (1 - p));
    model.addConstr(a3 >= 1 - (1 - p));

}

//weight: 1/64   2^-6
//Number of constraints : 223
void SBox_p11(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    model.addConstr(a1 + a3 + a5 + a7 >= 1 - (1 - p));
    model.addConstr(a0 + a1 + a3 - b5 >= 0 - (1 - p));
    model.addConstr(a0 - a1 + a3 + b5 >= 0 - (1 - p));
    model.addConstr(a2 - b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(a1 + a2 + a3 + b7 >= 1 - (1 - p));
    model.addConstr(a2 + a3 + b4 + b7 >= 1 - (1 - p));
    model.addConstr(a0 + a3 - b6 + b7 >= 0 - (1 - p));
    model.addConstr(a3 + b3 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(a3 + b4 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(a3 + b5 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(-a0 - a2 - a3 + b0 + b3 >= -2 - (1 - p));
    model.addConstr(a0 + a3 + a7 + b3 + b5 >= 1 - (1 - p));
    model.addConstr(-a2 - a3 + b0 + b3 + b6 >= -1 - (1 - p));
    model.addConstr(-a0 - a1 + a2 + a3 - b7 >= -2 - (1 - p));
    model.addConstr(a4 + a6 + a7 - b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a2 - a3 - b1 - b3 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a3 - b3 - b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 + a5 + a7 + b5 - b7 >= -1 - (1 - p));
    model.addConstr(a0 + a2 + b3 + b5 - b7 >= 0 - (1 - p));
    model.addConstr(-a3 + a7 + b3 + b6 - b7 >= -1 - (1 - p));
    model.addConstr(a2 - a3 + b1 + b3 + b7 >= 0 - (1 - p));
    model.addConstr(a1 + a3 + a6 - b4 + b7 >= 0 - (1 - p));
    model.addConstr(a2 + a7 - b2 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(a4 + a7 - b3 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(-a4 + a7 + b3 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(-a3 + a6 - b4 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a2 - b0 + b3 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(a2 + b3 + b4 + b6 + b7 >= 1 - (1 - p));
    model.addConstr(a6 - b4 + b5 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a6 + b4 + b5 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 - a4 + a7 - b1 + b3 >= -3 - (1 - p));
    model.addConstr(-a1 + a3 + a4 + a6 - b3 + b4 >= -1 - (1 - p));
    model.addConstr(-a3 + a4 + a6 - b1 - b3 - b5 >= -3 - (1 - p));
    model.addConstr(a1 + a2 + a3 + b1 - b3 - b5 >= -1 - (1 - p));
    model.addConstr(-a3 + a6 + a7 - b1 - b3 + b5 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 + a7 + b1 + b3 - b6 >= -2 - (1 - p));
    model.addConstr(-a0 + a2 + a3 + a6 + a7 - b7 >= -1 - (1 - p));
    model.addConstr(a1 + a2 - a5 - a7 + b2 - b7 >= -2 - (1 - p));
    model.addConstr(a1 + a2 + a6 - a7 + b2 - b7 >= -1 - (1 - p));
    model.addConstr(a0 + a3 + a6 + a7 - b3 - b7 >= -1 - (1 - p));
    model.addConstr(a2 - a5 - a7 + b2 - b3 - b7 >= -3 - (1 - p));
    model.addConstr(a2 + a6 - a7 + b2 - b3 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a2 + a3 + b0 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(a0 - a2 - a4 + b0 + b3 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a2 - a7 + b1 + b3 - b7 >= -3 - (1 - p));
    model.addConstr(a0 + a7 - b0 + b1 + b3 - b7 >= -1 - (1 - p));
    model.addConstr(-a1 + a4 + a7 - b3 - b4 - b7 >= -3 - (1 - p));
    model.addConstr(a0 + a3 - a7 - b1 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 + a2 - b0 + b1 - b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 + a4 + a6 - b3 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a3 + a4 + a7 - b3 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(a3 - a7 + b1 - b3 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a3 + a6 + a7 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 + a2 + a3 - b1 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(a2 + a5 - a6 - b2 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a2 + a3 + a5 - b6 - b7 >= -3 - (1 - p));
    model.addConstr(-a3 + a5 + a7 + b1 + b6 - b7 >= -1 - (1 - p));
    model.addConstr(-a0 - a3 - a4 + b3 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a2 - a3 - a7 + b2 + b4 + b7 >= -1 - (1 - p));
    model.addConstr(-a4 - a6 - a7 - b3 + b4 + b7 >= -3 - (1 - p));
    model.addConstr(a2 + a4 - b2 + b3 + b4 + b7 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 + a7 + b3 + b5 + b7 >= -1 - (1 - p));
    model.addConstr(a2 - b1 + b3 - b4 + b5 + b7 >= -1 - (1 - p));
    model.addConstr(a1 - a2 - a5 + a7 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 + a7 + b4 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a2 - b1 - b4 - b5 - b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a0 - b3 - b4 - b5 + b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a3 - b3 - b4 - b5 + b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a3 + a6 - a7 - b0 + b2 - b3 >= -4 - (1 - p));
    model.addConstr(-a0 - a2 + a3 - a4 + a7 + b1 + b3 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 - a2 + a3 - a6 - a7 + b4 >= -3 - (1 - p));
    model.addConstr(a2 - a3 + a7 - b1 - b2 - b3 - b5 >= -4 - (1 - p));
    model.addConstr(a2 - a3 - a7 - b1 + b2 - b3 - b5 >= -4 - (1 - p));
    model.addConstr(-a3 + a4 - a6 - a7 + b1 + b3 - b5 >= -3 - (1 - p));
    model.addConstr(-a2 - a3 + a4 + a7 + b1 - b4 - b5 >= -3 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a7 - b0 + b3 - b6 >= -4 - (1 - p));
    model.addConstr(a4 + a6 + b1 - b3 + b4 - b5 + b6 >= -1 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a5 + a7 + b5 + b6 >= -3 - (1 - p));
    model.addConstr(-a3 - a4 + a7 + b1 + b4 + b5 + b6 >= -1 - (1 - p));
    model.addConstr(-a2 - a3 + a5 + a7 + b0 + b2 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 + a5 + a7 - b0 - b2 - b3 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a2 - a3 + a4 - a7 + b3 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a3 - a4 + a6 - b1 + b3 - b7 >= -3 - (1 - p));
    model.addConstr(-a1 - a3 - a4 - a6 - b1 - b4 - b7 >= -6 - (1 - p));
    model.addConstr(a1 - a3 - a4 - a6 - b1 + b4 - b7 >= -4 - (1 - p));
    model.addConstr(-a3 - a4 - a6 - a7 - b3 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(-a1 - a3 - b1 - b3 - b4 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(a1 - a3 - b1 - b3 + b4 - b5 - b7 >= -4 - (1 - p));
    model.addConstr(-a1 + a3 - a7 - b1 - b3 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(-a1 + a3 + a7 + b1 - b3 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a1 + a7 - b1 + b3 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a1 - a3 - a6 - a7 - b4 + b5 - b7 >= -5 - (1 - p));
    model.addConstr(a1 - a3 + a6 - a7 - b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a1 - a3 - a6 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a1 - a3 + a6 - a7 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a4 + a5 + a7 + b1 + b2 - b6 - b7 >= -1 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a7 + b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 - a5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a3 + a6 + a7 + b1 + b6 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a1 - a3 - a6 - b4 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a1 - a3 - a6 + b4 + b6 - b7 >= -2 - (1 - p));
    model.addConstr(a0 - a3 - a4 + b1 + b5 + b6 - b7 >= -2 - (1 - p));
    model.addConstr(-a3 - a5 - a6 + a7 - b3 - b4 + b7 >= -4 - (1 - p));
    model.addConstr(-a3 - a6 - b1 - b3 + b4 + b5 + b7 >= -3 - (1 - p));
    model.addConstr(a4 + a7 - b1 + b3 - b4 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a3 - a6 + b1 + b4 - b5 - b6 + b7 >= -3 - (1 - p));
    model.addConstr(a2 + b1 + b2 + b4 - b5 - b6 + b7 >= -1 - (1 - p));
    model.addConstr(a2 + b1 - b3 - b4 + b5 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a2 - a7 + b0 - b2 + b4 + b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a2 + a7 + b0 + b2 + b4 + b6 + b7 >= 0 - (1 - p));
    model.addConstr(a7 - b0 - b2 - b3 + b4 + b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a0 - a2 - a3 + a6 - a7 + b0 - b2 - b3 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 + a4 + a6 - a7 - b1 - b3 - b5 >= -5 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a6 - a7 + b1 - b3 - b5 >= -6 - (1 - p));
    model.addConstr(-a1 - a3 + a4 - a6 - a7 + b3 + b4 - b5 >= -4 - (1 - p));
    model.addConstr(-a1 + a3 + a4 + a6 - a7 + b1 - b3 + b5 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 - a2 - a6 - a7 - b3 + b4 - b6 >= -5 - (1 - p));
    model.addConstr(-a3 - a5 + a6 - a7 - b0 + b2 - b3 + b6 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 - a5 - a7 + b0 - b2 - b5 + b6 >= -5 - (1 - p));
    model.addConstr(-a3 - a5 - a7 - b0 + b2 - b3 - b5 + b6 >= -5 - (1 - p));
    model.addConstr(a1 + a2 + a3 + a5 - a6 - a7 - b2 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a3 - a5 - a7 + b0 + b1 - b2 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a3 - a5 - a7 - b0 + b2 - b3 - b7 >= -6 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 + a7 - b1 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(a1 + a3 - a4 - a6 + b1 - b3 - b5 - b7 >= -4 - (1 - p));
    model.addConstr(a1 - a3 + a4 + a6 - a7 + b4 - b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a4 + a7 - b1 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a2 - a5 - a7 + b0 - b2 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a2 + a6 - a7 + b0 - b2 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a1 - a2 + a5 - a6 + b0 + b2 + b5 - b7 >= -2 - (1 - p));
    model.addConstr(-a2 - a3 + a5 - a6 + b0 + b2 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a1 + a5 - a6 - b0 - b2 - b3 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(-a3 + a5 - a6 - b0 - b2 - b3 + b5 - b7 >= -5 - (1 - p));
    model.addConstr(a0 - a5 - a7 - b0 + b2 - b3 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a6 - a7 - b0 + b2 - b3 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 + a4 - a6 + a7 - b1 + b3 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a1 - a4 - a6 + b1 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a1 - a3 + a4 - a6 - b3 + b4 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 + a1 - a2 - a3 - a7 - b4 - b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a3 - a6 - b3 - b4 - b6 - b7 >= -7 - (1 - p));
    model.addConstr(a0 - a1 - a2 - a3 - a7 + b4 - b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a1 - a3 + a7 - b3 - b4 - b5 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(a1 - a3 + a7 - b3 + b4 - b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(a0 + a1 - a3 - a4 + b1 - b4 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a1 - a3 - a4 + b1 + b4 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a2 - a3 + a5 - a7 + b0 - b2 + b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a3 + a5 - a7 - b0 + b2 - b3 + b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a2 - a3 - a5 + b0 + b2 - b4 + b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a3 - a5 - b0 - b2 - b3 - b4 + b7 >= -6 - (1 - p));
    model.addConstr(-a3 + a4 - a7 - b1 - b3 - b4 - b5 + b7 >= -5 - (1 - p));
    model.addConstr(-a2 - a3 - a4 + a6 + b3 + b4 - b5 + b7 >= -3 - (1 - p));
    model.addConstr(a0 - a2 + a5 - a7 - b0 - b2 - b6 + b7 >= -4 - (1 - p));
    model.addConstr(a0 + a5 - a7 + b0 + b2 - b3 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a0 + a5 + a6 - b0 + b2 + b3 - b6 + b7 >= -1 - (1 - p));
    model.addConstr(a0 - a2 - a5 - b0 + b2 - b4 - b6 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 + a4 - a7 - b3 + b4 - b6 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a7 + b3 + b4 - b6 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 + a4 + a7 - b1 - b4 + b5 - b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a4 + a7 + b1 - b3 - b4 + b5 - b6 + b7 >= -3 - (1 - p));
    model.addConstr(a2 + a4 + a6 + b1 + b4 + b5 - b6 + b7 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 + a5 - a7 + b0 - b2 + b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a3 + a5 - a7 - b0 + b2 - b3 + b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a2 - a5 - a6 + b0 + b2 + b5 + b6 + b7 >= -2 - (1 - p));
    model.addConstr(-a5 - a6 - b0 - b2 - b3 + b5 + b6 + b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 + a4 - a6 - a7 - b1 + b3 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 - a4 + a6 - a7 - b1 + b3 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 + a4 + a7 + b1 - b3 - b5 >= -4 - (1 - p));
    model.addConstr(-a1 - a3 + a4 + a6 - a7 - b0 - b1 - b4 - b5 >= -6 - (1 - p));
    model.addConstr(-a1 - a2 - a4 + a6 - a7 - b1 + b3 + b4 - b5 >= -5 - (1 - p));
    model.addConstr(-a1 - a2 - a4 + a6 - a7 - b1 - b3 - b4 + b5 >= -6 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a5 - a7 - b0 - b2 + b4 - b6 >= -6 - (1 - p));
    model.addConstr(a0 - a3 - a5 - a7 + b0 + b2 - b3 + b4 - b6 >= -4 - (1 - p));
    model.addConstr(-a0 - a1 - a2 - a4 + a6 - a7 + b3 + b4 - b6 >= -5 - (1 - p));
    model.addConstr(a0 - a2 - a3 + a4 - a5 - b0 - b1 - b5 - b6 >= -6 - (1 - p));
    model.addConstr(a0 - a3 - a5 + a7 + b0 + b1 - b2 - b5 - b6 >= -4 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a5 + a7 - b0 + b2 - b5 - b6 >= -5 - (1 - p));
    model.addConstr(-a0 - a2 + a3 + a4 - a6 + a7 - b1 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(a0 - a2 - a3 + a4 + a5 + b0 - b2 - b5 - b7 >= -4 - (1 - p));
    model.addConstr(a1 - a2 - a4 + a6 - a7 - b1 - b4 - b5 - b7 >= -6 - (1 - p));
    model.addConstr(a1 - a3 + a4 - a6 - a7 + b3 - b4 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 + a1 - a4 - a7 + b1 - b3 + b4 - b5 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a2 + a3 + a4 + a7 + b1 - b3 + b5 - b7 >= -3 - (1 - p));
    model.addConstr(a0 + a1 - a3 - a4 - a6 + b1 - b4 + b5 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a3 - a5 - a7 - b0 + b1 - b2 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(-a0 - a3 + a7 + b0 + b1 + b2 - b5 - b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a3 + a7 - b0 - b2 - b3 - b5 - b6 - b7 >= -7 - (1 - p));
    model.addConstr(a5 - a7 + b0 + b1 - b2 - b3 - b5 - b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a3 + a5 - a7 - b0 + b2 - b3 - b5 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(-a0 - a3 + a5 - a6 + b0 + b1 + b2 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(-a1 + a4 - a5 + a7 + b1 + b2 - b5 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a3 - a6 - b0 + b1 + b2 - b5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 + a1 - a3 - a6 + b1 - b4 - b5 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a1 - a3 - a6 + b1 + b4 - b5 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a0 - a2 - a3 + a4 - a7 - b3 + b5 + b6 - b7 >= -5 - (1 - p));
    model.addConstr(a2 - a3 - a5 - a6 - a7 - b2 - b3 - b4 + b7 >= -6 - (1 - p));
    model.addConstr(a2 - a3 + a5 - a6 - a7 + b2 - b3 - b4 + b7 >= -4 - (1 - p));
    model.addConstr(-a1 - a2 + a3 - a4 - a6 + b1 + b3 - b4 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a6 - b1 + b3 - b4 - b5 + b7 >= -6 - (1 - p));
    model.addConstr(-a0 - a3 + a4 + a5 + a6 + b0 + b2 - b6 + b7 >= -2 - (1 - p));
    model.addConstr(a0 - a2 - a5 - a6 + b0 - b2 - b3 - b6 + b7 >= -5 - (1 - p));
    model.addConstr(-a2 + a4 + a6 - a7 + b1 + b3 - b4 - b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a2 - a4 + a7 - b1 - b3 - b4 - b5 - b6 + b7 >= -6 - (1 - p));
    model.addConstr(-a2 - a3 + a4 - a6 - a7 + b1 + b5 - b6 + b7 >= -4 - (1 - p));
    model.addConstr(-a2 - a3 - a4 - a7 - b1 - b4 + b5 - b6 + b7 >= -6 - (1 - p));
    model.addConstr(a4 + a6 - a7 - b0 - b1 + b3 - b4 + b6 + b7 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 - a4 - a6 - a7 - b1 - b3 - b5 >= -8 - (1 - p));
    model.addConstr(-a0 - a1 - a2 - a4 + a6 - a7 + b1 - b3 - b4 - b5 >= -7 - (1 - p));
    model.addConstr(-a0 + a1 - a2 - a4 - a5 - a6 + a7 - b1 - b3 + b5 >= -6 - (1 - p));
    model.addConstr(-a0 - a1 - a2 + a3 - a4 - a6 - a7 + b1 - b3 + b5 >= -6 - (1 - p));
    model.addConstr(-a1 - a2 + a3 + a4 - a6 + b1 - b3 - b4 - b5 - b6 >= -6 - (1 - p));
    model.addConstr(-a2 - a3 - a4 + a6 - a7 - b0 + b1 + b3 - b5 + b6 >= -5 - (1 - p));
    model.addConstr(-a1 + a3 - a4 + a5 - a6 - a7 + b0 - b2 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(-a1 + a4 + a5 + a6 - a7 + b0 - b2 - b5 - b6 - b7 >= -5 - (1 - p));
    model.addConstr(-a1 + a4 - a5 + a6 - a7 + b0 + b1 - b2 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(-a0 - a2 - a3 + a5 - a6 - b0 - b2 - b3 + b6 - b7 >= -7 - (1 - p));
    model.addConstr(a0 - a4 - a5 - a6 - a7 + b0 - b2 - b5 + b6 - b7 >= -6 - (1 - p));
    model.addConstr(a0 - a3 + a5 + a6 - b0 + b1 - b2 - b5 + b6 - b7 >= -4 - (1 - p));
    model.addConstr(a0 - a4 + a5 + b0 + b1 + b2 - b3 - b5 + b6 - b7 >= -3 - (1 - p));
    model.addConstr(a0 - a2 - a3 - a5 - a7 + b0 + b1 + b2 - b3 - b6 - b7 >= -6 - (1 - p));
    model.addConstr(-a1 - a2 + a3 + a4 - a6 - b1 - b4 + b5 >= -4 - (1 - p));
    model.addConstr(a4 - a6 - a7 + b3 + b4 + b7 >= -1 - (1 - p));
    model.addConstr(a2 - a3 + b1 - b7 >= -1 - (1 - p));
    model.addConstr(-a0 + a1 - a2 + a3 - a7 - b7 >= -3 - (1 - p));

}


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
void cond_SBox_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x2 - x3 >= -1 - (1 - p));
    model.addConstr(-x3 - x4 >= -1 - (1 - p));
    model.addConstr(-x3 - x5 >= -1 - (1 - p));
    model.addConstr(x1 + x2 + x3 + x4 + x5 >= 1 - (1 - p));
    model.addConstr(x0 - x3 >= 0 - (1 - p));
    model.addConstr(-x0 + x3 >= 0 - (1 - p));
    model.addConstr(-x4 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x4 >= -1 - (1 - p));
    model.addConstr(-x1 - x2 >= -1 - (1 - p));
    model.addConstr(-x7 >= 0 - (1 - p));
    model.addConstr(-x6 >= 0 - (1 - p));
}
//weight: 1.830
void cond_SBox_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 - x4 + x7 >= 0 - (1 - p));
    model.addConstr(-x2 - x6 + x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x6 + x7 >= -1 - (1 - p));
    model.addConstr(-x5 + x6 + x7 >= 0 - (1 - p));
    model.addConstr(x0 - x2 + x3 + x5 >= 0 - (1 - p));
    model.addConstr(x0 + x3 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x1 - x3 + x5 - x6 >= -2 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(x1 - x3 - x5 >= -1 - (1 - p));
    model.addConstr(-x1 - x5 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x3 + x6 >= -1 - (1 - p));
    model.addConstr(-x0 - x7 >= -1 - (1 - p));
}
//weight: 2.000
void cond_SBox_p3(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x1 + x2 + x3 + x5 >= 1 - (1 - p));
    model.addConstr(x1 + x3 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x0 + x1 + x3 + x4 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x0 + x3 + x5 + x6 >= 0 - (1 - p));
    model.addConstr(x0 - x3 + x4 + x5 >= 0 - (1 - p));
    model.addConstr(x0 - x1 - x4 - x5 >= -2 - (1 - p));
    model.addConstr(-x1 + x6 - x7 >= -1 - (1 - p));
    model.addConstr(-x4 - x6 + x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x6 >= -1 - (1 - p));
    model.addConstr(x2 + x4 + x5 >= 1 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x4 - x5 >= -2 - (1 - p));
    model.addConstr(-x5 - x7 >= -1 - (1 - p));
    model.addConstr(-x0 - x7 >= -1 - (1 - p));
    model.addConstr(x2 - x6 >= 0 - (1 - p));
}
//weight: 2.186
void cond_SBox_p4(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(-x6 >= 0 - (1 - p));
    model.addConstr(-x5 >= 0 - (1 - p));
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(x2 >= 1 - (1 - p));
    model.addConstr(x1 >= 1 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.356
void cond_SBox_p5(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(-x6 >= 0 - (1 - p));
    model.addConstr(x5 >= 1 - (1 - p));
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(x1 >= 1 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.386
void cond_SBox_p6(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x2 - x3 + x4 >= -1 - (1 - p));
    model.addConstr(x1 + x3 + x4 >= 1 - (1 - p));
    model.addConstr(x2 + x3 + x6 >= 1 - (1 - p));
    model.addConstr(x0 - x3 - x7 >= -1 - (1 - p));
    model.addConstr(x1 + x3 - x7 >= 0 - (1 - p));
    model.addConstr(-x0 - x5 - x7 >= -2 - (1 - p));
    model.addConstr(-x2 + x5 - x7 >= -1 - (1 - p));
    model.addConstr(x3 + x5 - x7 >= 0 - (1 - p));
    model.addConstr(-x1 + x4 - x5 + x7 >= -1 - (1 - p));
    model.addConstr(x1 + x4 + x5 + x7 >= 1 - (1 - p));
    model.addConstr(x1 - x3 - x4 - x5 >= -2 - (1 - p));
    model.addConstr(-x1 - x4 + x5 + x7 >= -1 - (1 - p));
    model.addConstr(-x1 + x3 + x7 >= 0 - (1 - p));
    model.addConstr(x6 + x7 >= 1 - (1 - p));
}
//weight: 2.600
void cond_SBox_p7(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(x6 >= 1 - (1 - p));
    model.addConstr(-x5 >= 0 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(x1 >= 1 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.712
void cond_SBox_p8(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(x6 >= 1 - (1 - p));
    model.addConstr(-x5 >= 0 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
    model.addConstr(x2 >= 1 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.830
void cond_SBox_p9(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x3 - x7 >= 0 - (1 - p));
    model.addConstr(x0 - x5 - x7 >= -1 - (1 - p));
    model.addConstr(x1 - x4 + x7 >= 0 - (1 - p));
    model.addConstr(x2 + x4 + x7 >= 1 - (1 - p));
    model.addConstr(-x0 - x3 + x5 - x7 >= -2 - (1 - p));
    model.addConstr(-x1 - x3 + x5 - x6 - x7 >= -3 - (1 - p));
    model.addConstr(-x2 - x3 + x5 - x6 - x7 >= -3 - (1 - p));
    model.addConstr(x6 + x7 >= 1 - (1 - p));
    model.addConstr(-x3 - x4 + x7 >= -1 - (1 - p));
    model.addConstr(x3 + x4 + x7 >= 1 - (1 - p));
}
//weight: 2.952
void cond_SBox_p10(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(-x6 >= 0 - (1 - p));
    model.addConstr(x5 >= 1 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(-x1 >= 0 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 3.000
void cond_SBox_p11(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x3 - x6 >= 0 - (1 - p));
    model.addConstr(x2 + x5 - x7 >= 0 - (1 - p));
    model.addConstr(-x1 - x5 - x6 + x7 >= -2 - (1 - p));
    model.addConstr(x1 + x2 + x6 - x7 >= 0 - (1 - p));
    model.addConstr(x4 + x7 >= 1 - (1 - p));
    model.addConstr(-x0 - x5 - x7 >= -2 - (1 - p));
    model.addConstr(x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x1 + x5 + x7 >= 1 - (1 - p));
    model.addConstr(x0 + x5 - x7 >= 0 - (1 - p));
    model.addConstr(x1 + x3 >= 1 - (1 - p));
    model.addConstr(x3 - x7 >= 0 - (1 - p));
    model.addConstr(-x0 + x3 >= 0 - (1 - p));
}

//oCLAT
//weight: 1.000
void cond_inv_SBox_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x0 + x2 + x4 + x5 + x6 >= 1 - (1 - p));
    model.addConstr(-x0 - x6 >= -1 - (1 - p));
    model.addConstr(-x4 - x5 >= -1 - (1 - p));
    model.addConstr(-x0 - x2 >= -1 - (1 - p));
    model.addConstr(-x5 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x5 >= -1 - (1 - p));
    model.addConstr(-x0 - x5 >= -1 - (1 - p));
    model.addConstr(-x7 >= 0 - (1 - p));
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(-x1 >= 0 - (1 - p));
}
//weight: 1.386
void cond_inv_SBox_p2(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x4 - x6 - x7 >= -1 - (1 - p));
    model.addConstr(x0 - x3 - x5 >= -1 - (1 - p));
    model.addConstr(x1 + x3 >= 1 - (1 - p));
    model.addConstr(-x0 + x3 >= 0 - (1 - p));
    model.addConstr(-x4 + x7 >= 0 - (1 - p));
    model.addConstr(-x4 + x6 >= 0 - (1 - p));
    model.addConstr(-x3 - x7 >= -1 - (1 - p));
    model.addConstr(-x3 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
}
//weight: 1.830
void cond_inv_SBox_p3(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(-x1 >= 0 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.000
void cond_inv_SBox_p4(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x0 - x1 + x3 >= -1 - (1 - p));
    model.addConstr(-x1 + x2 + x3 >= 0 - (1 - p));
    model.addConstr(-x1 + x3 - x4 + x6 >= -1 - (1 - p));
    model.addConstr(-x1 + x3 - x4 + x7 >= -1 - (1 - p));
    model.addConstr(x0 + x1 + x3 + x5 + x7 >= 1 - (1 - p));
    model.addConstr(x2 + x4 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x2 + x3 + x5 + x6 + x7 >= 1 - (1 - p));
    model.addConstr(x0 + x2 + x3 + x4 + x6 >= 1 - (1 - p));
    model.addConstr(-x1 + x4 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x5 - x6 + x7 >= -2 - (1 - p));
    model.addConstr(-x3 - x4 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x3 + x4 + x7 >= -1 - (1 - p));
    model.addConstr(-x0 - x2 - x5 >= -2 - (1 - p));
    model.addConstr(x0 + x2 - x7 >= 0 - (1 - p));
    model.addConstr(-x0 - x2 - x7 >= -2 - (1 - p));
    model.addConstr(-x3 - x6 >= -1 - (1 - p));
    model.addConstr(-x2 - x3 >= -1 - (1 - p));
}
//weight: 2.386
void cond_inv_SBox_p5(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x2 - x3 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x3 + x6 - x7 >= -2 - (1 - p));
    model.addConstr(-x3 - x4 + x6 - x7 >= -2 - (1 - p));
    model.addConstr(x2 - x3 + x6 + x7 >= 0 - (1 - p));
    model.addConstr(x0 - x3 + x4 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(x0 - x3 + x5 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(x0 + x3 - x4 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x2 + x3 >= -1 - (1 - p));
    model.addConstr(x0 + x3 + x4 + x7 >= 1 - (1 - p));
    model.addConstr(x0 - x3 - x5 + x7 >= -1 - (1 - p));
    model.addConstr(x0 + x3 + x4 + x6 >= 1 - (1 - p));
    model.addConstr(x1 + x3 >= 1 - (1 - p));
}
//weight: 2.600
void cond_inv_SBox_p6(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x7 >= 0 - (1 - p));
    model.addConstr(x6 >= 1 - (1 - p));
    model.addConstr(x5 >= 1 - (1 - p));
    model.addConstr(-x4 >= 0 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 2.830
void cond_inv_SBox_p7(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x4 + x5 >= 1 - (1 - p));
    model.addConstr(x4 - x7 >= 0 - (1 - p));
    model.addConstr(-x0 + x7 >= 0 - (1 - p));
    model.addConstr(x2 + x7 >= 1 - (1 - p));
    model.addConstr(-x2 - x3 - x4 >= -2 - (1 - p));
    model.addConstr(-x5 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(-x0 - x6 >= -1 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
}
//weight: 3.000
void cond_inv_SBox_p8(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(-x0 - x1 + x2 >= -1 - (1 - p));
    model.addConstr(-x3 + x4 + x7 >= 0 - (1 - p));
    model.addConstr(-x3 + x5 + x7 >= 0 - (1 - p));
    model.addConstr(-x0 + x2 + x4 - x6 - x7 >= -2 - (1 - p));
    model.addConstr(x0 - x2 + x4 - x6 >= -1 - (1 - p));
    model.addConstr(x1 + x5 + x7 >= 1 - (1 - p));
    model.addConstr(x2 - x4 - x7 >= -1 - (1 - p));
    model.addConstr(-x0 - x3 + x7 >= -1 - (1 - p));
    model.addConstr(x0 + x3 >= 1 - (1 - p));
    model.addConstr(x2 + x6 >= 1 - (1 - p));
}
//weight: 3.386
void cond_inv_SBox_p9(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& x4, GRBVar& x5, GRBVar& x6, GRBVar& x7, GRBVar& p)
{
    model.addConstr(x7 >= 1 - (1 - p));
    model.addConstr(x6 >= 1 - (1 - p));
    model.addConstr(-x4 >= 0 - (1 - p));
    model.addConstr(x3 >= 1 - (1 - p));
    model.addConstr(x2 >= 1 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}

void SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4, GRBVar& p5, GRBVar& p6, GRBVar& p7, GRBVar& p8, GRBVar& p9, GRBVar& p10, GRBVar& p11)
{
    SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    SBox_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p2);
    SBox_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p3);
    SBox_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p4);
    SBox_p5(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p5);
    SBox_p6(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p6);
    SBox_p7(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p7);
    SBox_p8(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p8);
    SBox_p9(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p9);
    SBox_p10(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p10);
    SBox_p11(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p11);
    model.addConstr(p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 == 1);
}

void cond_SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4, GRBVar& p5, GRBVar& p6, GRBVar& p7, GRBVar& p8, GRBVar& p9, GRBVar& p10, GRBVar& p11)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_SBox_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_SBox_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);
    cond_SBox_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, p3);
    cond_SBox_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, p4);
    cond_SBox_p5(model, a0, a1, a2, a3, a4, a5, a6, a7, p5);
    cond_SBox_p6(model, a0, a1, a2, a3, a4, a5, a6, a7, p6);
    cond_SBox_p7(model, a0, a1, a2, a3, a4, a5, a6, a7, p7);
    cond_SBox_p8(model, a0, a1, a2, a3, a4, a5, a6, a7, p8);
    cond_SBox_p9(model, a0, a1, a2, a3, a4, a5, a6, a7, p9);
    cond_SBox_p10(model, a0, a1, a2, a3, a4, a5, a6, a7, p10);
    cond_SBox_p11(model, a0, a1, a2, a3, a4, a5, a6, a7, p11);
    model.addConstr(p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 == 1);
}

void cond_inv_SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3, GRBVar& p4, GRBVar& p5, GRBVar& p6, GRBVar& p7, GRBVar& p8, GRBVar& p9)
{
    cond_SBox_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, p0);
    cond_inv_SBox_p1(model, a0, a1, a2, a3, a4, a5, a6, a7, p1);
    cond_inv_SBox_p2(model, a0, a1, a2, a3, a4, a5, a6, a7, p2);
    cond_inv_SBox_p3(model, a0, a1, a2, a3, a4, a5, a6, a7, p3);
    cond_inv_SBox_p4(model, a0, a1, a2, a3, a4, a5, a6, a7, p4);
    cond_inv_SBox_p5(model, a0, a1, a2, a3, a4, a5, a6, a7, p5);
    cond_inv_SBox_p6(model, a0, a1, a2, a3, a4, a5, a6, a7, p6);
    cond_inv_SBox_p7(model, a0, a1, a2, a3, a4, a5, a6, a7, p7);
    cond_inv_SBox_p8(model, a0, a1, a2, a3, a4, a5, a6, a7, p8);
    cond_inv_SBox_p9(model, a0, a1, a2, a3, a4, a5, a6, a7, p9);
    model.addConstr(p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 == 1);
}


//Skinny SubCell_Sbox
void SC(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        SBox(model, s[8 * cell + 0], s[8 * cell + 1], s[8 * cell + 2], s[8 * cell + 3], s[8 * cell + 4], s[8 * cell + 5], s[8 * cell + 6], s[8 * cell + 7],
            t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
            p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4], p[cell][5], p[cell][6], p[cell][7], p[cell][8], p[cell][9], p[cell][10], p[cell][11]);
    }
}
void cond_SC(GRBModel& model, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        cond_SBox(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
            p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4], p[cell][5], p[cell][6], p[cell][7], p[cell][8], p[cell][9], p[cell][10], p[cell][11]);
    }
}
void cond_inv_SC(GRBModel& model, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        cond_inv_SBox(model, t[8 * cell + 0], t[8 * cell + 1], t[8 * cell + 2], t[8 * cell + 3], t[8 * cell + 4], t[8 * cell + 5], t[8 * cell + 6], t[8 * cell + 7],
            p[cell][0], p[cell][1], p[cell][2], p[cell][3], p[cell][4], p[cell][5], p[cell][6], p[cell][7], p[cell][8], p[cell][9]);
    }
}

//Skinny ShiftRow
void SR(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    for (int bit = 0; bit < 8; bit++)
    {
        model.addConstr(t2[0 * 8 + bit] == t1[0 * 8 + bit]);
        model.addConstr(t2[1 * 8 + bit] == t1[13 * 8 + bit]);
        model.addConstr(t2[2 * 8 + bit] == t1[10 * 8 + bit]);
        model.addConstr(t2[3 * 8 + bit] == t1[7 * 8 + bit]);

        model.addConstr(t2[4 * 8 + bit] == t1[4 * 8 + bit]);
        model.addConstr(t2[5 * 8 + bit] == t1[1 * 8 + bit]);
        model.addConstr(t2[6 * 8 + bit] == t1[14 * 8 + bit]);
        model.addConstr(t2[7 * 8 + bit] == t1[11 * 8 + bit]);

        model.addConstr(t2[8 * 8 + bit] == t1[8 * 8 + bit]);
        model.addConstr(t2[9 * 8 + bit] == t1[5 * 8 + bit]);
        model.addConstr(t2[10 * 8 + bit] == t1[2 * 8 + bit]);
        model.addConstr(t2[11 * 8 + bit] == t1[15 * 8 + bit]);

        model.addConstr(t2[12 * 8 + bit] == t1[12 * 8 + bit]);
        model.addConstr(t2[13 * 8 + bit] == t1[9 * 8 + bit]);
        model.addConstr(t2[14 * 8 + bit] == t1[6 * 8 + bit]);
        model.addConstr(t2[15 * 8 + bit] == t1[3 * 8 + bit]);
    }

}


//Skinny MixColumn
void MC(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 4; col++)
    {
        for (int bit = 0; bit < 8; bit++)
        {
            XOR3(model, s[(4 * col + 0) * 8 + bit], s[(4 * col + 1) * 8 + bit], s[(4 * col + 3) * 8 + bit], t[(4 * col + 0) * 8 + bit]);
            model.addConstr(s[(4 * col + 2) * 8 + bit] == t[(4 * col + 1) * 8 + bit]);
            XOR3(model, s[(4 * col + 0) * 8 + bit], s[(4 * col + 2) * 8 + bit], s[(4 * col + 3) * 8 + bit], t[(4 * col + 2) * 8 + bit]);
            model.addConstr(s[(4 * col + 0) * 8 + bit] == t[(4 * col + 3) * 8 + bit]);
        }
    }

}



double Search_Skinny_trail(int round, bitset<300> Template, double min)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    //env.set(GRB_DoubleParam_TimeLimit, 500);
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
            p[i][j] = model.addVars(12, GRB_BINARY);
        }
    }
    p[0] = new GRBVar * [16];
    p[round + 1] = new GRBVar * [16];
    for (j = 0; j < 16; j++)
    {
        p[0][j] = model.addVars(12, GRB_BINARY);
        p[round + 1][j] = model.addVars(10, GRB_BINARY);
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
    cond_inv_SC(model, s[round], p[round + 1]);     //oCLAT

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
            for (int j = 1; j < 12; j++)
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
                for (int j = 1; j < 12; j++)
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
            for (int j = 1; j < 10; j++)
                sum += p[round + 1][cell][j];
            model.addConstr(sum == 1);
        }
    }


    GRBLinExpr cor = 0;
    for (i = 1; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            //2  2.186  2.599  2.83  3.081  3.66  4  4.385  4.83  5.356  6
            cor += 2 * p[i][j][1] + 2.186 * p[i][j][2] + 2.599 * p[i][j][3] + 2.83 * p[i][j][4] + 3.081 * p[i][j][5] + 3.66 * p[i][j][6]
                + 4 * p[i][j][7] + 4.385 * p[i][j][8] + 4.83 * p[i][j][9] + 5.356 * p[i][j][10] + 6 * p[i][j][11];
        }
    }
    for (j = 0; j < 16; j++)
    {
        //1.000  1.830  2.000  2.186  2.356  2.386  2.600  2.712  2.830  2.952  3.000
        cor += p[0][j][1] + 1.830 * p[0][j][2] + 2.000 * p[0][j][3] + 2.186 * p[0][j][4] + 2.356 * p[0][j][5] + 2.386 * p[0][j][6]
            + 2.600 * p[0][j][7] + 2.712 * p[0][j][8] + 2.830 * p[0][j][9] + 2.952 * p[0][j][10] + 3.000 * p[0][j][11];

        //1.000  1.386  1.830  2.000  2.386  2.600  2.830  3.000  3.386
        cor += p[round + 1][j][1] + 1.386 * p[round + 1][j][2] + 1.830 * p[round + 1][j][3] + 2 * p[round + 1][j][4] + 2.386 * p[round + 1][j][5]
            + 2.6 * p[round + 1][j][6] + 2.83 * p[round + 1][j][7] + 3 * p[round + 1][j][8] + 3.386 * p[round + 1][j][9];
    }


    model.addConstr(cor <= 119);

    model.setObjective(cor, GRB_MINIMIZE);


    model.update();
    model.optimize();



    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {

        string sss = "results.txt";
        ofstream f(sss, ios::app);
        f << "===== Skinny " << round + 2 << "-round result of the linear distinguisher =====" << std::endl;
        f << "The S-box output mask for round 1: " << std::endl;
        for (int n1 = 0; n1 < 16; n1++)
        {
            int result = 0;
            for (int n2 = 0; n2 < 8; n2++)
            {
                result = result * 2 + (int)t1[0][8 * n1 + n2].get(GRB_DoubleAttr_X);
            }
            f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
        }
        f << endl;
        f << "The S-box input mask for rounds 2 to " << round + 2 << std::endl;
        for (int i = 0; i < round + 1; i++)
        {
            for (int n1 = 0; n1 < 16; n1++)
            {
                int result = 0;
                for (int n2 = 0; n2 < 8; n2++)
                {
                    result = result * 2 + (int)s[i][8 * n1 + n2].get(GRB_DoubleAttr_X);
                }
                f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
            }
            f << endl;
        }

        int ccccc = 0;
        for (j = 0; j < 16; j++)
        {
            ccccc = ccccc + 1 - (int)p[0][j][0].get(GRB_DoubleAttr_X) + 1 - (int)p[round + 1][j][0].get(GRB_DoubleAttr_X);
        }

        f << "\nConditional linear weight: " << cor.getValue() << std::endl;
        f << "Active Conditional S-box: " << dec << ccccc << std::endl;

        f.close();

        return cor.getValue();
    }

    return min;
}






int main()
{


    int round;
    cout << " ROUND: ";
    cin >> round;

    vector<bitset<300>> foundSolutions;
    
    double min = 119;
    for (int i = 0; i < 50; i++)
    {
        cout << "-------------------------- Try " << i << " --------------------------" << endl;
        bitset<300> temp = Search_Skinny_Sbox(round, foundSolutions);
        foundSolutions.push_back(temp);
        min = Search_Skinny_trail(round - 2, temp, min);
    }
    cout << "Minimum Conditional linear weight: " << min << endl;
    return 0;
}

