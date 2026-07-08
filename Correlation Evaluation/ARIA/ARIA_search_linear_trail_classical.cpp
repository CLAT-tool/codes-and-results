#include <iostream>
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
    model.addConstr(a + b == 2 * t + o);
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

void Diffusion(GRBModel& model, GRBVar* t, GRBVar* s)
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

        for (int bit = 0; bit < 8; bit++)
        {
            XOR7(model, s[ones_positions[0] * 8 + bit], s[ones_positions[1] * 8 + bit], s[ones_positions[2] * 8 + bit],
                s[ones_positions[3] * 8 + bit], s[ones_positions[4] * 8 + bit], s[ones_positions[5] * 8 + bit], s[ones_positions[6] * 8 + bit], t[cell * 8 + bit]);
        }
    }

}


bitset<300> Search_ARIA_Sbox(int round, vector<bitset<300>> foundSolutions)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);



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
        Diffusion(model, t_state[loc], s_state[loc + 1]);
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
                    cut += (1 - sbox[i][n1]);
                }
                else
                {
                    cut += sbox[i][n1];
                }
            }
        }
        model.addConstr(cut >= 1);
    }

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
    for (i = 0; i < round; i++)
    {
        for (j = 0; j < 16; j++)
        {
            condition_cor += 6 * sbox[i][j];
        }
    }
    model.setObjective(condition_cor, GRB_MINIMIZE);



    model.update();
    model.optimize();

    bitset<300> temp(0);

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        for (int i = 0; i < round; i++)
        {
            for (int n1 = 0; n1 < 16; n1++)
            {
                if (sbox[i][n1].get(GRB_DoubleAttr_X) > 0.5)
                {
                    temp[16 * i + n1] = 1;
                }
                else
                {
                    temp[16 * i + n1] = 0;
                }
            }
        }
        cout << "\nOptimal weight: " << condition_cor.getValue() << std::endl;
    }
    else
    {
        cout << "No solution£¡£¡£¡£¡£¡£¡£¡£¡£¡" << endl;
    }
    return temp;
}


//Inactive sbox
void SBox_p_null(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
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

// ARIA Sbox 1 
// Number of constraints : 
// Input  : a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output : b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb

//weight: 1/64   2^-6
//Number of constraints : 1889
void SBox1_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)

    // Read the preprocessed constraint file
    std::ifstream infile("milp_read_AES_Sbox_linear_1-64.txt");
    std::string line;
    int constraint_count = 0;
    while (std::getline(infile, line)) {
        // Skip over blank lines and comment lines
        if (line.empty() || line[0] == '#') continue;

        // Split line: Coefficient part | Right-hand constant | p flag
        std::istringstream iss(line);
        std::string coeffs_str;
        if (!std::getline(iss, coeffs_str, '|')) continue;

        double rhs_constant;
        if (!(iss >> rhs_constant)) continue;

        // Skip the '|' delimiter
        char separator;
        if (!(iss >> separator) || separator != '|') continue;

        int has_p_flag;
        if (!(iss >> has_p_flag)) continue;

        // Analyze 17 coefficients
        std::istringstream coeff_stream(coeffs_str);
        std::vector<double> coeffs(NUM_VARS, 0.0);
        for (int i = 0; i < NUM_VARS; ++i) {
            if (!(coeff_stream >> coeffs[i])) {
                // Error in processing: Insufficient number of coefficients
                break;
            }
        }

        // Construct the linear expression lhs
        GRBLinExpr lhs = 0;
        for (int i = 0; i < NUM_VARS; ++i) {
            if (coeffs[i] != 0.0) {
                lhs += coeffs[i] * (*vars[i]);
            }
        }

        if (has_p_flag == 1) {
            lhs += 1 - p;
        }

        // Add constraints to the model
        model.addConstr(lhs >= rhs_constant,
            "SBox1_p0_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 49/4096   2^-6.385
//Number of constraints : 2498
void SBox1_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
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
//weight: 9/1024   2^-6.830
//Number of constraints : 4164
void SBox1_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_AES_Sbox_linear_9-1024.txt");
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
            "SBox1_p2_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

// ARIA Sbox2 
// Number of constraints : 
// Input  : a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output : b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb

//weight: 1/64   2^-6
//Number of constraints : 2028-> 1059
void SBox2_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_Sbox2_linear_1-64.txt");
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
            "SBox2_p0_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 49/4096   2^-6.385
//Number of constraints : 2508
void SBox2_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_Sbox2_linear_49-4096.txt");
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
            "SBox2_p1_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 9/1024   2^-6.830
//Number of constraints : 4182
void SBox2_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_Sbox2_linear_9-1024.txt");
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
            "SBox2_p2_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

// ARIA inv_Sbox 1 
// Number of constraints : 
// Input  : a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output : b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// 
//weight: 1/64   2^-6
//Number of constraints : 1054
void invSBox1_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_AES_invSbox_linear_1-64.txt");
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
            "invSBox1_p0_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 49/4096   2^-6.385
//Number of constraints : 2494
void invSBox1_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_AES_invSbox_linear_49-4096.txt");
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
            "invSBox1_p1_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 9/1024   2^-6.830
//Number of constraints : 4173
void invSBox1_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_AES_invSbox_linear_9-1024.txt");
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
            "invSBox1_p2_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

// ARIA Sbox2inverse 
// Number of constraints : 
// Input  : a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output : b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// 
//weight: 1/64   2^-6
//Number of constraints : 1063
void invSBox2_p0(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_invSbox2_linear_1-64.txt");
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
            "invSBox2_p0_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 49/4096   2^-6.385
//Number of constraints : 2499
void invSBox2_p1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_invSbox2_linear_49-4096.txt");
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
            "invSBox2_p1_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}
//weight: 9/1024   2^-6.830
//Number of constraints : 4201
void invSBox2_p2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7, GRBVar& p)
{
    std::vector<GRBVar*> vars = { &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                             &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7,
                             &p };
    const int NUM_VARS = 17; // a0-a7 (8) + b0-b7 (8) + p (1)


    std::ifstream infile("milp_read_ARIA_invSbox2_linear_9-1024.txt");
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
            "invSBox2_p2_cons_" + std::to_string(constraint_count));

        constraint_count++;
    }

    infile.close();
}

// Input:   a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output:  b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// P: 6, 6.385, 6.830, 7.356, 8, 8.830, 10, 12
void SBox1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3)
{
    SBox_p_null(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox1_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);

    model.addConstr(p0 + p1  == 1);
}

// Input:   a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output:  b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// P: 6, 6.385, 6.830, 7.356, 8, 8.830, 10, 12
void SBox2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3)
{
    SBox_p_null(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    SBox2_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);

    model.addConstr(p0 + p1  == 1);
}

// Input:   a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output:  b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// P: 6, 6.385, 6.830, 7.356, 8, 8.830, 10, 12
void invSBox1(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3)
{
    SBox_p_null(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    invSBox1_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);

    model.addConstr(p0 + p1  == 1);
}

// Input:   a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7; a0: msb
// Output:  b0 || b1 || b2 || b3 || b4 || b5 || b6 || b7; b0: msb
// P: 6, 6.385, 6.830, 7.356, 8, 8.830, 10, 12
void invSBox2(GRBModel& model, GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
    GRBVar& b0, GRBVar& b1, GRBVar& b2, GRBVar& b3, GRBVar& b4, GRBVar& b5, GRBVar& b6, GRBVar& b7,
    GRBVar& p0, GRBVar& p1, GRBVar& p2, GRBVar& p3)
{
    SBox_p_null(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p0);
    invSBox2_p0(model, a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7, p1);
    model.addConstr(p0 + p1  == 1);
}

//ARIA Substitution Layer
void Substitution1(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)   //Type 1
{
    for (int i = 0; i < 4; i++)
    {
        SBox1(model, s[32 * i], s[32 * i + 1], s[32 * i + 2], s[32 * i + 3], s[32 * i + 4], s[32 * i + 5], s[32 * i + 6], s[32 * i + 7],
            t[32 * i], t[32 * i + 1], t[32 * i + 2], t[32 * i + 3], t[32 * i + 4], t[32 * i + 5], t[32 * i + 6], t[32 * i + 7],
            p[4 * i][0], p[4 * i][1]);
        SBox2(model, s[32 * i + 8], s[32 * i + 1 + 8], s[32 * i + 2 + 8], s[32 * i + 3 + 8], s[32 * i + 4 + 8], s[32 * i + 5 + 8], s[32 * i + 6 + 8], s[32 * i + 7 + 8],
            t[32 * i + 8], t[32 * i + 1 + 8], t[32 * i + 2 + 8], t[32 * i + 3 + 8], t[32 * i + 4 + 8], t[32 * i + 5 + 8], t[32 * i + 6 + 8], t[32 * i + 7 + 8],
            p[4 * i + 1][0], p[4 * i + 1][1]);
        invSBox1(model, s[32 * i + 16], s[32 * i + 1 + 16], s[32 * i + 2 + 16], s[32 * i + 3 + 16], s[32 * i + 4 + 16], s[32 * i + 5 + 16], s[32 * i + 6 + 16], s[32 * i + 7 + 16],
            t[32 * i + 16], t[32 * i + 1 + 16], t[32 * i + 2 + 16], t[32 * i + 3 + 16], t[32 * i + 4 + 16], t[32 * i + 5 + 16], t[32 * i + 6 + 16], t[32 * i + 7 + 16],
            p[4 * i + 2][0], p[4 * i + 2][1]);
        invSBox2(model, s[32 * i + 24], s[32 * i + 1 + 24], s[32 * i + 2 + 24], s[32 * i + 3 + 24], s[32 * i + 4 + 24], s[32 * i + 5 + 24], s[32 * i + 6 + 24], s[32 * i + 7 + 24],
            t[32 * i + 24], t[32 * i + 1 + 24], t[32 * i + 2 + 24], t[32 * i + 3 + 24], t[32 * i + 4 + 24], t[32 * i + 5 + 24], t[32 * i + 6 + 24], t[32 * i + 7 + 24],
            p[4 * i + 3][0], p[4 * i + 3][1]);
    }
}
void Substitution2(GRBModel& model, GRBVar* s, GRBVar* t, GRBVar** p)   //Type 2
{
    for (int i = 0; i < 4; i++)
    {
        invSBox1(model, s[32 * i], s[32 * i + 1], s[32 * i + 2], s[32 * i + 3], s[32 * i + 4], s[32 * i + 5], s[32 * i + 6], s[32 * i + 7],
            t[32 * i], t[32 * i + 1], t[32 * i + 2], t[32 * i + 3], t[32 * i + 4], t[32 * i + 5], t[32 * i + 6], t[32 * i + 7],
            p[4 * i][0], p[4 * i][1]);
        invSBox2(model, s[32 * i + 8], s[32 * i + 1 + 8], s[32 * i + 2 + 8], s[32 * i + 3 + 8], s[32 * i + 4 + 8], s[32 * i + 5 + 8], s[32 * i + 6 + 8], s[32 * i + 7 + 8],
            t[32 * i + 8], t[32 * i + 1 + 8], t[32 * i + 2 + 8], t[32 * i + 3 + 8], t[32 * i + 4 + 8], t[32 * i + 5 + 8], t[32 * i + 6 + 8], t[32 * i + 7 + 8],
            p[4 * i + 1][0], p[4 * i + 1][1]);
        SBox1(model, s[32 * i + 16], s[32 * i + 1 + 16], s[32 * i + 2 + 16], s[32 * i + 3 + 16], s[32 * i + 4 + 16], s[32 * i + 5 + 16], s[32 * i + 6 + 16], s[32 * i + 7 + 16],
            t[32 * i + 16], t[32 * i + 1 + 16], t[32 * i + 2 + 16], t[32 * i + 3 + 16], t[32 * i + 4 + 16], t[32 * i + 5 + 16], t[32 * i + 6 + 16], t[32 * i + 7 + 16],
            p[4 * i + 2][0], p[4 * i + 2][1]);
        SBox2(model, s[32 * i + 24], s[32 * i + 1 + 24], s[32 * i + 2 + 24], s[32 * i + 3 + 24], s[32 * i + 4 + 24], s[32 * i + 5 + 24], s[32 * i + 6 + 24], s[32 * i + 7 + 24],
            t[32 * i + 24], t[32 * i + 1 + 24], t[32 * i + 2 + 24], t[32 * i + 3 + 24], t[32 * i + 4 + 24], t[32 * i + 5 + 24], t[32 * i + 6 + 24], t[32 * i + 7 + 24],
            p[4 * i + 3][0], p[4 * i + 3][1]);
    }
}


double Search_ARIA_cond_trail(int round, bitset<300> Template, double min)
{
    int i, j;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);

    GRBModel model = GRBModel(env);

    //Add an additional layer of diffusion on top.
    GRBVar** s = new GRBVar * [round + 2];    //Previous state of the S-box
    GRBVar** t = new GRBVar * [round + 1];        //State after the S-box
    for (i = 0; i < round + 1; i++)
    {
        s[i] = model.addVars(128, GRB_BINARY);
        t[i] = model.addVars(128, GRB_BINARY);
    }
    s[round+1] = model.addVars(128, GRB_BINARY);

    GRBVar*** p;
    p = new GRBVar * *[round + 1];        //Two layers of conditional linear probabilities are extended up and down
    for (i = 0; i < round + 1; i++)
    {
        p[i] = new GRBVar * [16];
    }
    for (i = 0; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            p[i][j] = model.addVars(2, GRB_BINARY);
        }
    }
    

    for (int loc = 0; loc < round+1; loc++)
    {
        if (loc % 2 == 1)
        {
            Substitution2(model, s[loc], t[loc], p[loc]);
        }
        else
        {
            Substitution1(model, s[loc], t[loc], p[loc]);
        }
        Diffusion(model, t[loc], s[loc + 1]);
    }

    GRBLinExpr init_mask_sum = 0;
    for (i = 0; i < 128; i++)
        init_mask_sum += s[0][i];
    model.addConstr(init_mask_sum >= 1);

    int sbox_index[18][16] = { 0 };

    for (i = 0; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            sbox_index[i][j] = Template[16 * i + j];
        }

    }


    for (int i = 0; i < round + 1; i++)
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
                for (int j = 1; j < 2; j++)
                    sum += p[i][cell][j];
                model.addConstr(sum == 1);
            }
        }
    }


    GRBLinExpr cor = 0;
    for (i = 0; i < round + 1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            cor += 6 * p[i][j][1];
        }
    }


    model.addConstr(cor <= 127);

    model.setObjective(cor, GRB_MINIMIZE);

    model.update();
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        string sss = "results.txt";
        ofstream f(sss, ios::app);
        f << "===== Classical ARIA " << round + 1 << "-round result of the linear distinguisher =====" << std::endl;
        f << "The S-box input and output masks for rounds 1 to " << round + 2 << ": " << std::endl;
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
                f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
            }
            f << endl;
            for (int n1 = 0; n1 < 16; n1++)
            {
                int result = 0;
                for (int n2 = 0; n2 < 8; n2++)
                {
                    if (t[i][8 * n1 + n2].get(GRB_DoubleAttr_X) >= 0.5)
                        result = result * 2 + 1;
                    else
                        result = result * 2;
                }
                f << " 0x" << hex << setw(2) << setfill('0') << uppercase << result << " ";
            }
            f << endl;
            f << endl;
        }

        int ccccc = 0;
        for (j = 0; j < 16; j++)
        {
            ccccc = ccccc + 1 - (int)p[0][j][0].get(GRB_DoubleAttr_X) + 1 - (int)p[round + 1][j][0].get(GRB_DoubleAttr_X);
        }

        f << "\nLinear weight: " << cor.getValue() << std::endl;
        f << "Active Conditional S-box: " << ccccc << std::endl;


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
    for (int i = 0; i < 1000; i++)
    {
        cout << "-------------------------- Try " << dec << i << " --------------------------" << endl;
        bitset<300> temp = Search_ARIA_Sbox(round, foundSolutions);
        foundSolutions.push_back(temp);
        min = Search_ARIA_cond_trail(round - 1, temp, min);
    }
    cout << endl << "===================================================" << endl;
    cout << "The final linear weight:" << min << endl;
    cout << "===================================================" << endl;
    return 0;
}

