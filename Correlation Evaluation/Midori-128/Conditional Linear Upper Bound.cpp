#include <iostream>
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
void OR5(GRBModel& model, GRBVar& in0, GRBVar& in1, GRBVar& in2, GRBVar& in3, GRBVar& out)
{
    model.addConstr(2 * out >= in0 + in1);
    model.addConstr(2 * out >= in2 + in3);
    model.addConstr(out <= in0 + in1 + in2 + in3);
}

//Midori SubCell_Sbox
void SC(GRBModel& model, GRBVar* s, GRBVar* t)
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
        XOR2(model, s[4 * col + 2], s[4 * col + 3], d0);		//d = s2 + s3
        XOR2(model, s[4 * col + 1], d0, t[4 * col + 0]);		//t0 = s1 + s2 + s3 = s0 + d
        XOR2(model, s[4 * col + 0], d0, t[4 * col + 1]);		//t1 = s0 + s2 + s3 = s0 + d

        GRBVar d1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        XOR2(model, s[4 * col + 0], s[4 * col + 1], d1);		//d1 = s0 + s1
        XOR2(model, s[4 * col + 3], d1, t[4 * col + 2]);		//t2 = s0 + s1 + s3 = s3 + d
        XOR2(model, s[4 * col + 2], d1, t[4 * col + 3]);		//t3 = s0 + s1 + s2 = s2 + d

        //GRBVar d2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 1], s[4 * col + 2], d2);		//d2 = s1 + s2
        //XOR2(model, s[4 * col + 3], d2, t[4 * col + 0]);		//t0 = s1 + s2 + s3 = s3 + d
        //XOR2(model, s[4 * col + 0], d2, t[4 * col + 3]);		//t3 = s0 + s1 + s2 = s0 + d

        //GRBVar d3 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 0], s[4 * col + 2], d3);		//d3 = s0 + s2
        //XOR2(model, s[4 * col + 3], d3, t[4 * col + 1]);		//t1 = s0 + s2 + s3 = s3 + d
        //XOR2(model, s[4 * col + 1], d3, t[4 * col + 3]);		//t3 = s0 + s1 + s2 = s1 + d

        //GRBVar d4 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 0], s[4 * col + 3], d4);		//d4 = s0 + s3
        //XOR2(model, s[4 * col + 2], d4, t[4 * col + 1]);		//t1 = s0 + s2 + s3 = s2 + d
        //XOR2(model, s[4 * col + 1], d4, t[4 * col + 2]);		//t2 = s0 + s1 + s3 = s1 + d

        //GRBVar d5 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //XOR2(model, s[4 * col + 1], s[4 * col + 3], d5);		//d4 = s1 + s3
        //XOR2(model, s[4 * col + 1], d5, t[4 * col + 0]);		//t0 = s1 + s2 + s3 = s1 + d
        //XOR2(model, s[4 * col + 0], d5, t[4 * col + 2]);		//t2 = s0 + s1 + s3 = s0 + d
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


void Search_Midori_Sbox(int round)
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
        t1[i] = model.addVars(16, GRB_BINARY);
        t2[i] = model.addVars(16, GRB_BINARY);
    }
    s[round] = model.addVars(16, GRB_BINARY);





    for (int loc = 0; loc < round; loc++)
    {
        SC(model, s[loc], t1[loc]);
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


    model.update();
    model.optimize();



    std::cout << "===== Midori " << round << "-round minimum S-box result of the linear distinguisher =====" << std::endl;
    int sbox_number = 0;

    for (int i = 0; i < round; i++)
    {
        for (int n1 = 0; n1 < 16; n1++)
        {
            sbox_number += (int)s[i][n1].get(GRB_DoubleAttr_X);
            std::cout << setw(1) << setfill('0') << (int)s[i][n1].get(GRB_DoubleAttr_X) << " ";
        }
        cout << endl;
        for (int n1 = 0; n1 < 16; n1++)
        {
            std::cout << setw(1) << setfill('0') << (int)t1[i][n1].get(GRB_DoubleAttr_X) << " ";
        }
        cout << endl;
        for (int n1 = 0; n1 < 16; n1++)
        {
            std::cout << setw(1) << setfill('0') << (int)t2[i][n1].get(GRB_DoubleAttr_X) << " ";
        }
        cout << endl;
        cout << endl;
    }
    std::cout << "\nNumber of active sboxes: " << sbox_number << std::endl;
    std::cout << "\nOptimal weight: " << condition_cor.getValue() << std::endl;
}






int main()
{
    int round = 10;
    Search_Midori_Sbox(round);

    return 0;
}

