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
    model.addConstr(a + b + o == 2 * t);
    model.addConstr(t >= a);
    model.addConstr(t >= b);
    model.addConstr(t >= o);
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


void Search_LED_Sbox(int round)
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
        SC_Sbox(model, s[loc], t1[loc]);
        SR_Sbox(model, t1[loc], t2[loc]);
        MC_Sbox(model, t2[loc], s[loc + 1]);
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

    //model.addConstr(cor <= 50);
    //model.addConstr(cor >= 30);

    model.setObjective(cor, GRB_MINIMIZE);


    model.update();
    model.optimize();


    int sbox_number = 0;
    std::cout << "===== LED " << round << "-round Sbox result of the linear distinguisher =====" << std::endl;
    std::cout << "The S-box output mask for round 1: " << std::endl;
    for (int i = 0; i < round; i++)
    {
        for (int n1 = 0; n1 < 16; n1++)
        {
            sbox_number += (int)s[i][n1].get(GRB_DoubleAttr_X);
            std::cout << setw(1) << setfill('0') << (int)s[i][n1].get(GRB_DoubleAttr_X) << " ";
        }
        cout << endl;
        cout << endl;
    }
    std::cout << "\nNumber of active sboxes: " << sbox_number << std::endl;
    std::cout << "\nOptimal weight: " << cor.getValue() << std::endl;
}



int main()
{
    int round = 6;

    Search_LED_Sbox(round);  

    return 0;
}
