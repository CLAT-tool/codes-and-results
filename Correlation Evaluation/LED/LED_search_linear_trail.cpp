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
    model.addConstr(a + b + o == 2 * t);
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



//LED SubCells
void SC_Sbox(GRBModel& model, GRBVar* s, GRBVar* t)
{
    for (int cell = 0; cell < 16; cell++)
    {
        model.addConstr(s[cell] == t[cell]);
    }
}

//LED ShiftRows
void SR_Sbox(GRBModel& model, GRBVar* t1, GRBVar* t2)
{

    model.addConstr(t2[0] == t1[0]);
    model.addConstr(t2[1] == t1[5]);
    model.addConstr(t2[2] == t1[10]);
    model.addConstr(t2[3] == t1[15]);

    model.addConstr(t2[4] == t1[4]);
    model.addConstr(t2[5] == t1[9]);
    model.addConstr(t2[6] == t1[14]);
    model.addConstr(t2[7] == t1[3]);

    model.addConstr(t2[8] == t1[8]);
    model.addConstr(t2[9] == t1[13]);
    model.addConstr(t2[10] == t1[2]);
    model.addConstr(t2[11] == t1[7]);

    model.addConstr(t2[12] == t1[12]);
    model.addConstr(t2[13] == t1[1]);
    model.addConstr(t2[14] == t1[6]);
    model.addConstr(t2[15] == t1[11]);


}

//LED MixColumns
void MC_Sbox(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 4; col++)
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



bitset<300> Search_LED_Sbox(int round, vector<bitset<300>> foundSolutions)
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
        init_mask_sum += t1[0][i];
    model.addConstr(init_mask_sum >= 1);



    GRBLinExpr cor = 0;
    for (i = 1; i < round - 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            //2  4
            cor += 2 * s[i][j];
        }
    }
    for (j = 0; j < 16; j++)
    {
        //1 
        cor += s[0][j];

        //1 
        cor += s[round - 1][j];
    }

 
    model.setObjective(cor, GRB_MINIMIZE);



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
        cout << "\nOptimal weight: " << cor.getValue() << std::endl;
    }
    else
    {
        cout << "No solution!!!!!" << endl;
    }

    return temp;
}






void SBox_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& p)
{
    model.addConstr(1 - p >= a0);
    model.addConstr(1 - p >= a1);
    model.addConstr(1 - p >= a2);
    model.addConstr(1 - p >= a3);

    model.addConstr(1 - p >= b0);
    model.addConstr(1 - p >= b1);
    model.addConstr(1 - p >= b2);
    model.addConstr(1 - p >= b3);

}

//weight: 1/4   2^-2
//Number of constraints : 
void SBox_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& p)
{
    model.addConstr(-a0 - a3 - b1 >= -2 - (1 - p));
    model.addConstr(a1 + a2 + a3 + b1 >= 1 - (1 - p));
    model.addConstr(a0 + a1 + b0 + b1 >= 1 - (1 - p));
    model.addConstr(a0 + a2 + b0 + b1 >= 1 - (1 - p));
    model.addConstr(-a0 + b0 - b1 - b2 >= -2 - (1 - p));
    model.addConstr(-a0 + a1 + a3 + b2 >= 0 - (1 - p));
    model.addConstr(a0 + a1 + b1 + b2 >= 1 - (1 - p));
    model.addConstr(-a1 - b0 - b1 - b3 >= -3 - (1 - p));
    model.addConstr(-a0 + b0 - b2 - b3 >= -2 - (1 - p));
    model.addConstr(a2 + a3 + b2 - b3 >= 0 - (1 - p));
    model.addConstr(a0 + a1 + a2 + b3 >= 1 - (1 - p));
    model.addConstr(a2 + b0 - b1 + b3 >= 0 - (1 - p));
    model.addConstr(a0 + a3 + b1 + b3 >= 1 - (1 - p));
    model.addConstr(a1 - b1 + b2 + b3 >= 0 - (1 - p));
    model.addConstr(-a0 - a1 + a3 + b0 + b1 >= -1 - (1 - p));
    model.addConstr(-a0 + a1 - a2 - a3 - b2 >= -3 - (1 - p));
    model.addConstr(a0 + a1 + a3 - b1 - b2 >= -1 - (1 - p));
    model.addConstr(a1 - a2 - b0 - b1 - b2 >= -3 - (1 - p));
    model.addConstr(-a0 - a1 - a2 - b0 + b2 >= -3 - (1 - p));
    model.addConstr(a0 - a2 + b0 - b1 - b3 >= -2 - (1 - p));
    model.addConstr(a1 - a3 - b0 + b1 - b3 >= -2 - (1 - p));
    model.addConstr(a2 - a3 - b0 + b1 - b3 >= -2 - (1 - p));
    model.addConstr(-a1 - a2 - b0 - b2 - b3 >= -4 - (1 - p));
    model.addConstr(-a1 + a2 + a3 - b1 + b3 >= -1 - (1 - p));
    model.addConstr(a0 - a1 - a2 + b1 + b3 >= -1 - (1 - p));
    model.addConstr(a1 - a3 + b0 + b1 + b3 >= 0 - (1 - p));
    model.addConstr(-a2 - a3 + b0 + b1 + b3 >= -1 - (1 - p));
    model.addConstr(-a2 + a3 - b0 - b2 + b3 >= -2 - (1 - p));
    model.addConstr(a0 + a2 - b1 - b2 + b3 >= -1 - (1 - p));
    model.addConstr(a2 - a3 + b1 + b2 + b3 >= 0 - (1 - p));
    model.addConstr(a0 - a1 - a3 + b0 - b1 - b2 >= -3 - (1 - p));
    model.addConstr(a0 - a2 - a3 - b0 - b1 + b2 >= -3 - (1 - p));
    model.addConstr(a0 - a1 - a3 + b0 + b2 - b3 >= -2 - (1 - p));
    model.addConstr(-a0 - a1 + a2 - b0 + b1 - b2 + b3 >= -3 - (1 - p));
    model.addConstr(a0 + a3 + b0 + b2 >= 1 - (1 - p));
}
//weight: 1/16   2^-4
//Number of constraints : 
void SBox_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& p)
{
    model.addConstr(-a0 - a1 + a2 + b1 >= -1 - (1 - p));
    model.addConstr(-a0 + a1 - a2 + b1 >= -1 - (1 - p));
    model.addConstr(a0 - b0 - b1 + b2 >= -1 - (1 - p));
    model.addConstr(a1 + a2 + b0 + b2 >= 1 - (1 - p));
    model.addConstr(-a1 - a2 + b0 + b2 >= -1 - (1 - p));
    model.addConstr(a0 + b0 - b1 - b2 >= -1 - (1 - p));
    model.addConstr(a1 + a2 - b0 - b2 >= -1 - (1 - p));
    model.addConstr(-a1 - a2 - b0 - b2 >= -3 - (1 - p));
    model.addConstr(b0 + b1 + b2 >= 1 - (1 - p));
    model.addConstr(-b0 + b1 - b2 >= -1 - (1 - p));
    model.addConstr(a0 + a1 + a2 >= 1 - (1 - p));
    model.addConstr(a0 - a1 - a2 >= -1 - (1 - p));

}


void cond_SBox_p0(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& p)
{
    model.addConstr(-x3 >= 0 - (1 - p));
    model.addConstr(-x2 >= 0 - (1 - p));
    model.addConstr(-x1 >= 0 - (1 - p));
    model.addConstr(-x0 >= 0 - (1 - p));
}
//weight: 1.000
void cond_SBox_p1(GRBModel& model, GRBVar& x0, GRBVar& x1, GRBVar& x2, GRBVar& x3, GRBVar& p)
{
    model.addConstr(x0 + x1 + x2 + x3 >= 1 - (1 - p));
}

void SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& p0, GRBVar& p1, GRBVar& p2)
{
    SBox_p0(model, a0, a1, a2, a3, b0, b1, b2, b3, p0);
    SBox_p1(model, a0, a1, a2, a3, b0, b1, b2, b3, p1);
    SBox_p2(model, a0, a1, a2, a3, b0, b1, b2, b3, p2);

    model.addConstr(p0 + p1 + p2 == 1);
}

void cond_SBox(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& p0, GRBVar& p1)
{
    cond_SBox_p0(model, a0, a1, a2, a3, p0);
    cond_SBox_p1(model, a0, a1, a2, a3, p1);

    model.addConstr(p0 + p1 == 1);
}


//LED SubCells
void SC(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        SBox(model, s[4 * cell + 0], s[4 * cell + 1], s[4 * cell + 2], s[4 * cell + 3],
            t[4 * cell + 0], t[4 * cell + 1], t[4 * cell + 2], t[4 * cell + 3],
            p[cell][0], p[cell][1], p[cell][2]);

    }
}
void cond_SC(GRBModel& model, GRBVar* t, GRBVar** p)
{
    for (int cell = 0; cell < 16; cell++)
    {
        cond_SBox(model, t[4 * cell + 0], t[4 * cell + 1], t[4 * cell + 2], t[4 * cell + 3],
            p[cell][0], p[cell][1]);
    }
}

//LED ShiftRows
void SR(GRBModel& model, GRBVar* t1, GRBVar* t2)
{
    for (int bit = 0; bit < 4; bit++)
    {
        model.addConstr(t2[0 * 4 + bit] == t1[0 * 4 + bit]);
        model.addConstr(t2[1 * 4 + bit] == t1[5 * 4 + bit]);
        model.addConstr(t2[2 * 4 + bit] == t1[10 * 4 + bit]);
        model.addConstr(t2[3 * 4 + bit] == t1[15 * 4 + bit]);

        model.addConstr(t2[4 * 4 + bit] == t1[4 * 4 + bit]);
        model.addConstr(t2[5 * 4 + bit] == t1[9 * 4 + bit]);
        model.addConstr(t2[6 * 4 + bit] == t1[14 * 4 + bit]);
        model.addConstr(t2[7 * 4 + bit] == t1[3 * 4 + bit]);

        model.addConstr(t2[8 * 4 + bit] == t1[8 * 4 + bit]);
        model.addConstr(t2[9 * 4 + bit] == t1[13 * 4 + bit]);
        model.addConstr(t2[10 * 4 + bit] == t1[2 * 4 + bit]);
        model.addConstr(t2[11 * 4 + bit] == t1[7 * 4 + bit]);

        model.addConstr(t2[12 * 4 + bit] == t1[12 * 4 + bit]);
        model.addConstr(t2[13 * 4 + bit] == t1[1 * 4 + bit]);
        model.addConstr(t2[14 * 4 + bit] == t1[6 * 4 + bit]);
        model.addConstr(t2[15 * 4 + bit] == t1[11 * 4 + bit]);

    }

}

const int MC_matrix[16][16] = {
  {0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,1},
  {0,0,1,1,0,1,1,0,1,0,1,0,1,0,0,0},
  {1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0},
  {0,1,0,0,1,0,0,0,1,0,1,1,0,0,1,0},
  {1,0,0,0,0,1,0,1,1,0,0,1,0,0,1,1},
  {0,1,0,0,1,0,1,1,1,1,0,1,1,0,0,0},
  {0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0},
  {0,0,0,1,0,1,1,0,1,1,1,0,0,0,1,0},
  {0,0,1,1,1,1,1,0,1,1,1,1,0,0,0,1},
  {1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,1},
  {0,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1},
  {0,0,1,0,0,1,0,1,1,0,1,0,1,1,1,1},
  {0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1},
  {1,0,0,0,1,0,1,1,0,0,1,0,1,0,1,0},
  {0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,1},
  {0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1}
};

//LED MixColumns
void MC(GRBModel& model, GRBVar* t, GRBVar* s)
{
    for (int col = 0; col < 4; col++)
    {
        for (int bit_out = 0; bit_out < 16; bit_out++)
        {

            std::vector<int> ones_positions;
            for (int j = 0; j < 16; j++)
            {
                if (MC_matrix[bit_out][j] == 1)
                {
                    ones_positions.push_back(j);
                }
            }

            if (ones_positions.empty())
            {
                model.addConstr(t[16 * col + bit_out] == 0, "MC_zero_" +
                    std::to_string(col) + "_" + std::to_string(bit_out));
            }
            else if (ones_positions.size() == 1)
            {

                int pos = ones_positions[0];
                model.addConstr(t[16 * col + bit_out] == s[16 * col + pos], "MC_single_" +
                    std::to_string(col) + "_" + std::to_string(bit_out));
            }
            else
            {

                GRBVar current = s[16 * col + ones_positions[0]];

                for (size_t i = 1; i < ones_positions.size(); i++)
                {
                    if (i == ones_positions.size() - 1)
                    {

                        XOR2(model, current, s[16 * col + ones_positions[i]], t[16 * col + bit_out]);
                    }
                    else
                    {

                        GRBVar temp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                        XOR2(model, current, s[16 * col + ones_positions[i]], temp);
                        current = temp;
                    }
                }
            }
        }
    }

}




double Search_LED_cond_trail(int round, bitset<300> Template, double min)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    env.set(GRB_DoubleParam_TimeLimit, 1000);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);
    GRBVar** s = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t1 = new GRBVar * [round + 1];        //State after the S-box
    GRBVar** t2 = new GRBVar * [round + 1];        //State after the SR
    for (i = 0; i < round + 1; i++)
    {
        s[i] = model.addVars(16 * 4, GRB_BINARY);
        t1[i] = model.addVars(16 * 4, GRB_BINARY);
        t2[i] = model.addVars(16 * 4, GRB_BINARY);
    }
    GRBVar*** p;
    p = new GRBVar * *[round + 2];
    for (i = 1; i < round + 1; i++)
    {
        p[i] = new GRBVar * [16];
        for (j = 0; j < 16; j++)
        {
            p[i][j] = model.addVars(3, GRB_BINARY);
        }
    }
    p[0] = new GRBVar * [16];
    p[round + 1] = new GRBVar * [16];
    for (j = 0; j < 16; j++)
    {
        p[0][j] = model.addVars(2, GRB_BINARY);
        p[round + 1][j] = model.addVars(2, GRB_BINARY);
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
    for (i = 0; i < 16 * 4; i++)
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
            for (int j = 1; j < 2; j++)
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
                for (int j = 1; j < 3; j++)
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
            for (int j = 1; j < 2; j++)
                sum += p[round + 1][cell][j];
            model.addConstr(sum == 1);
        }
    }


    GRBLinExpr cor = 0;
    for (i = 1; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            //2  4
            cor += 2 * p[i][j][1] + 4 * p[i][j][2];
        }
    }
    for (j = 0; j < 16; j++)
    {
        //1 
        cor += p[0][j][1];

        //1 
        cor += p[round + 1][j][1];
    }

    model.addConstr(cor <= min);

    model.setObjective(cor, GRB_MINIMIZE);


    model.update();
    model.optimize();


    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        string sss = "results.txt";
        ofstream f(sss, ios::app);
        f << "===== LED " << round + 2 << "-round result of the linear distinguisher =====" << std::endl;
        f << "The S-box output mask for round 1: " << std::endl;
        for (int n1 = 0; n1 < 16; n1++)
        {
            int result = 0;
            for (int n2 = 0; n2 < 4; n2++)
            {
                if (t1[0][4 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
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
            for (int n1 = 0; n1 < 16; n1++)
            {
                int result = 0;
                for (int n2 = 0; n2 < 4; n2++)
                {
                    if (s[i][4 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
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


    int round;
    cout << " ROUND: ";
    cin >> round;

    vector<bitset<300>> foundSolutions;

    double min = 200;
    for (int i = 0; i < 200; i++)
    {
        cout << "-------------------------- Try " << dec << i << " --------------------------" << endl;
        bitset<300> temp = Search_LED_Sbox(round, foundSolutions);
        foundSolutions.push_back(temp);
        min = Search_LED_cond_trail(round - 2, temp, min);
    }
    cout << endl << "===================================================" << endl;
    cout << "The final Conditional linear weight:" << min << endl;
    cout << "===================================================" << endl;
    return 0;
}

