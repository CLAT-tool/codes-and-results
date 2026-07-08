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
    model.addConstr(a + b + o == 2 * t);
}
void XOR7(GRBModel& model, GRBVar& in0, GRBVar& in1, GRBVar& in2, GRBVar& in3, GRBVar& in4, GRBVar& in5, GRBVar& in6, GRBVar& out)
{
    GRBVar t = model.addVar(0.0, 3.0, 0.0, GRB_INTEGER);
    model.addConstr(in0 + in1 + in2 + in3 + in4 + in5 + in6 == 2 * t + out);
}
void SC_Sbox(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar* active)
{
    for (int cell = 0; cell < 16; cell++)
    {
        GRBLinExpr sum0 = 0;
        GRBLinExpr sum1 = 0;
        for (int bit = 0; bit < 8; bit++)
        {
            model.addConstr(active[cell] - s[8 * cell + bit] >= 0);
            sum0 += s[8 * cell + bit];
            sum1 += t[8 * cell + bit];
        }
        model.addConstr(sum0 - active[cell] >= 0);
        model.addConstr(8 * sum0 - sum1 >= 0);
        model.addConstr(8 * sum1 - sum0 >= 0);
    }
}

void Diffusion_Sbox(GRBModel& model, GRBVar* t, GRBVar* s)
{
    //s: Previous state of the S-box
    //t: State after the S-box
    int MC_matrix[16][16] = {
    {0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0},
    {0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1},
    {0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1},
    {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0},
    {1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1},
    {0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1},
    {1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0},
    {0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0},
    {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1},
    {1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0},
    {0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1},
    {0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0},
    {0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0},
    {1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0},
    {1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0},
    {0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1}
    };


    for (int cell = 0; cell < 16; cell++)  
    {

        std::vector<int> ones_positions;
        for (int j = 0; j < 16; j++)
        {
            if (MC_matrix[cell][j] == 1)
            {
                ones_positions.push_back(j);
            }
        }
        // XOR£ºa XOR b XOR c = ((a XOR b) XOR c)
        for (int bit = 0; bit < 8; bit++)
        {
            XOR7(model, s[ones_positions[0] * 8 + bit], s[ones_positions[1] * 8 + bit], s[ones_positions[2] * 8 + bit],
                s[ones_positions[3] * 8 + bit], s[ones_positions[4] * 8 + bit], s[ones_positions[5] * 8 + bit], s[ones_positions[6] * 8 + bit], t[cell * 8 + bit]);

            //// The first variable
            //GRBVar current = s[8 * ones_positions[0] + bit];

            //// Chain-based calculation XOR
            //for (size_t i = 1; i < ones_positions.size(); i++)
            //{
            //    if (i == ones_positions.size() - 1)
            //    {
            //        // The last XOR operation, and the result is assigned to the output
            //        XOR2(model, current, s[8 * ones_positions[i] + bit], t[8 * cell + bit]);
            //    }
            //    else
            //    {
            //        // Intermediate result requires the creation of a temporary variable
            //        GRBVar temp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            //        XOR2(model, current, s[8 * ones_positions[i] + bit], temp);
            //        current = temp;
            //    }
            //}
        }
    }

    //int index[7];
    //for (int i = 0; i < 16; i++)
    //{
    //    for (int j = 0; j < 16; j++)
    //    {
    //        int k = 0;
    //        if (diffusion_matrix[j][i] == 1)
    //        {
    //            index[k] = j;
    //            k++;
    //        }
    //    }
    //    for (int j = 0; j < 8; j++)
    //    {
    //        XOR8(model, s[index[0] * 8 + j], s[index[1] * 8 + j], s[index[2] * 8 + j], s[index[3] * 8 + j], s[index[4] * 8 + j], s[index[5] * 8 + j], s[index[6] * 8 + j], t[i * 8 + j]);
    //    }
    //}

}


void Search_ARIA_Sbox(int round)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);
    //env.set(GRB_DoubleParam_TimeLimit, 100);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, 2000000000);


    GRBModel model = GRBModel(env);

    GRBVar** sbox = new GRBVar * [round];       //Active indication
    GRBVar** s_state = new GRBVar * [round + 1];    //Previous state of the S-box
    GRBVar** t_state = new GRBVar * [round];        //State after the S-box
    for (i = 0; i < round; i++)
    {
        sbox[i] = model.addVars(16, GRB_BINARY);
        s_state[i] = model.addVars(16 * 8, GRB_BINARY);
        t_state[i] = model.addVars(16 * 8, GRB_BINARY);
    }
    s_state[round] = model.addVars(16 * 8, GRB_BINARY);


    for (int loc = 0; loc < round; loc++)
    {
        SC_Sbox(model, s_state[loc], t_state[loc], sbox[loc]);
        Diffusion_Sbox(model, t_state[loc], s_state[loc + 1]);
    }

    /*
    GRBLinExpr constr = 0;
    for (i = 0; i < 16; i++)
        constr += p0[round / 2][i] + p0[round / 2][i];
    model.addConstr(constr == 1);
    */

    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 16 * 8; i++)
        init_mask_sum += s_state[0][i];
    model.addConstr(init_mask_sum >= 1);

    GRBLinExpr sbox_sum = 0;
    for (i = 0; i < round; i++)
        for (j = 0; j < 16; j++)
            sbox_sum += sbox[i][j];
    if (round % 2 == 1)
        model.addConstr(sbox_sum == 8 * int(round / 2) + 1);
    else
        model.addConstr(sbox_sum == 8 * int(round / 2));

    for (i = 0; i < round - 1; i++)
    {
        GRBLinExpr sbox_sum0 = 0;
        for (int cell = 0; cell < 16; cell++)
        {
            sbox_sum0 += sbox[i][cell] + sbox[i + 1][cell];
        }
        model.addConstr(sbox_sum0 == 8);
    }


    GRBLinExpr condition_cor = 0;
    for (i = 1; i < round - 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            condition_cor += 6 * sbox[i][j];
        }
    }
    for (j = 0; j < 16; j++)
    {
        condition_cor += 3.83 * sbox[0][j];
        condition_cor += 3.83 * sbox[round - 1][j];
    }
    model.setObjective(condition_cor, GRB_MINIMIZE);


    model.update();
    model.optimize();



    std::cout << "===== ARIA " << round << "-round minimum S-box result of the linear distinguisher =====" << std::endl;
    int sbox_number = 0;

    for (int i = 0; i < round; i++)
    {
        for (int n1 = 0; n1 < 16; n1++)
        {
            sbox_number += (int)sbox[i][n1].get(GRB_DoubleAttr_X);
            std::cout << setw(1) << setfill('0') << (int)sbox[i][n1].get(GRB_DoubleAttr_X) << " ";
        }
        cout << endl;
        cout << endl;
    }
    std::cout << "\nNumber of active sboxes: " << sbox_number << std::endl;
    std::cout << "\nOptimal correlation: " << condition_cor.getValue() << std::endl;
}

int main()
{
    int round = 7;

    Search_ARIA_Sbox(round);
    return 0;
}