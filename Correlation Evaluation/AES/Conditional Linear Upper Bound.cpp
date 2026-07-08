#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <map>
#include <mutex>
#include <time.h>
#include <string>
#include"gurobi_c++.h"
using namespace std;

int main()
{

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_LogToConsole, 1);


	GRBModel model = GRBModel(env);
	GRBVar* x = model.addVars(250, GRB_BINARY);
	GRBVar* d = model.addVars(50, GRB_BINARY);


	model.addConstr(x[0] + x[5] + x[10] + x[15] + x[16] + x[17] + x[18] + x[19] - 5 * d[0] >= 0);
	model.addConstr(d[0] - x[0] >= 0);
	model.addConstr(d[0] - x[5] >= 0);
	model.addConstr(d[0] - x[10] >= 0);
	model.addConstr(d[0] - x[15] >= 0);
	model.addConstr(d[0] - x[16] >= 0);
	model.addConstr(d[0] - x[17] >= 0);
	model.addConstr(d[0] - x[18] >= 0);
	model.addConstr(d[0] - x[19] >= 0);
	model.addConstr(x[1] + x[6] + x[11] + x[12] + x[20] + x[21] + x[22] + x[23] - 5 * d[1] >= 0);
	model.addConstr(d[1] - x[1] >= 0);
	model.addConstr(d[1] - x[6] >= 0);
	model.addConstr(d[1] - x[11] >= 0);
	model.addConstr(d[1] - x[12] >= 0);
	model.addConstr(d[1] - x[20] >= 0);
	model.addConstr(d[1] - x[21] >= 0);
	model.addConstr(d[1] - x[22] >= 0);
	model.addConstr(d[1] - x[23] >= 0);
	model.addConstr(x[2] + x[7] + x[8] + x[13] + x[24] + x[25] + x[26] + x[27] - 5 * d[2] >= 0);
	model.addConstr(d[2] - x[2] >= 0);
	model.addConstr(d[2] - x[7] >= 0);
	model.addConstr(d[2] - x[8] >= 0);
	model.addConstr(d[2] - x[13] >= 0);
	model.addConstr(d[2] - x[24] >= 0);
	model.addConstr(d[2] - x[25] >= 0);
	model.addConstr(d[2] - x[26] >= 0);
	model.addConstr(d[2] - x[27] >= 0);
	model.addConstr(x[3] + x[4] + x[9] + x[14] + x[28] + x[29] + x[30] + x[31] - 5 * d[3] >= 0);
	model.addConstr(d[3] - x[3] >= 0);
	model.addConstr(d[3] - x[4] >= 0);
	model.addConstr(d[3] - x[9] >= 0);
	model.addConstr(d[3] - x[14] >= 0);
	model.addConstr(d[3] - x[28] >= 0);
	model.addConstr(d[3] - x[29] >= 0);
	model.addConstr(d[3] - x[30] >= 0);
	model.addConstr(d[3] - x[31] >= 0);
	model.addConstr(x[16] + x[21] + x[26] + x[31] + x[32] + x[33] + x[34] + x[35] - 5 * d[4] >= 0);
	model.addConstr(d[4] - x[16] >= 0);
	model.addConstr(d[4] - x[21] >= 0);
	model.addConstr(d[4] - x[26] >= 0);
	model.addConstr(d[4] - x[31] >= 0);
	model.addConstr(d[4] - x[32] >= 0);
	model.addConstr(d[4] - x[33] >= 0);
	model.addConstr(d[4] - x[34] >= 0);
	model.addConstr(d[4] - x[35] >= 0);
	model.addConstr(x[20] + x[25] + x[30] + x[19] + x[36] + x[37] + x[38] + x[39] - 5 * d[5] >= 0);
	model.addConstr(d[5] - x[20] >= 0);
	model.addConstr(d[5] - x[25] >= 0);
	model.addConstr(d[5] - x[30] >= 0);
	model.addConstr(d[5] - x[19] >= 0);
	model.addConstr(d[5] - x[36] >= 0);
	model.addConstr(d[5] - x[37] >= 0);
	model.addConstr(d[5] - x[38] >= 0);
	model.addConstr(d[5] - x[39] >= 0);
	model.addConstr(x[24] + x[29] + x[18] + x[23] + x[40] + x[41] + x[42] + x[43] - 5 * d[6] >= 0);
	model.addConstr(d[6] - x[24] >= 0);
	model.addConstr(d[6] - x[29] >= 0);
	model.addConstr(d[6] - x[18] >= 0);
	model.addConstr(d[6] - x[23] >= 0);
	model.addConstr(d[6] - x[40] >= 0);
	model.addConstr(d[6] - x[41] >= 0);
	model.addConstr(d[6] - x[42] >= 0);
	model.addConstr(d[6] - x[43] >= 0);
	model.addConstr(x[28] + x[17] + x[22] + x[27] + x[44] + x[45] + x[46] + x[47] - 5 * d[7] >= 0);
	model.addConstr(d[7] - x[28] >= 0);
	model.addConstr(d[7] - x[17] >= 0);
	model.addConstr(d[7] - x[22] >= 0);
	model.addConstr(d[7] - x[27] >= 0);
	model.addConstr(d[7] - x[44] >= 0);
	model.addConstr(d[7] - x[45] >= 0);
	model.addConstr(d[7] - x[46] >= 0);
	model.addConstr(d[7] - x[47] >= 0);
	model.addConstr(x[32] + x[37] + x[42] + x[47] + x[48] + x[49] + x[50] + x[51] - 5 * d[8] >= 0);
	model.addConstr(d[8] - x[32] >= 0);
	model.addConstr(d[8] - x[37] >= 0);
	model.addConstr(d[8] - x[42] >= 0);
	model.addConstr(d[8] - x[47] >= 0);
	model.addConstr(d[8] - x[48] >= 0);
	model.addConstr(d[8] - x[49] >= 0);
	model.addConstr(d[8] - x[50] >= 0);
	model.addConstr(d[8] - x[51] >= 0);
	model.addConstr(x[36] + x[41] + x[46] + x[35] + x[52] + x[53] + x[54] + x[55] - 5 * d[9] >= 0);
	model.addConstr(d[9] - x[36] >= 0);
	model.addConstr(d[9] - x[41] >= 0);
	model.addConstr(d[9] - x[46] >= 0);
	model.addConstr(d[9] - x[35] >= 0);
	model.addConstr(d[9] - x[52] >= 0);
	model.addConstr(d[9] - x[53] >= 0);
	model.addConstr(d[9] - x[54] >= 0);
	model.addConstr(d[9] - x[55] >= 0);
	model.addConstr(x[40] + x[45] + x[34] + x[39] + x[56] + x[57] + x[58] + x[59] - 5 * d[10] >= 0);
	model.addConstr(d[10] - x[40] >= 0);
	model.addConstr(d[10] - x[45] >= 0);
	model.addConstr(d[10] - x[34] >= 0);
	model.addConstr(d[10] - x[39] >= 0);
	model.addConstr(d[10] - x[56] >= 0);
	model.addConstr(d[10] - x[57] >= 0);
	model.addConstr(d[10] - x[58] >= 0);
	model.addConstr(d[10] - x[59] >= 0);
	model.addConstr(x[44] + x[33] + x[38] + x[43] + x[60] + x[61] + x[62] + x[63] - 5 * d[11] >= 0);
	model.addConstr(d[11] - x[44] >= 0);
	model.addConstr(d[11] - x[33] >= 0);
	model.addConstr(d[11] - x[38] >= 0);
	model.addConstr(d[11] - x[43] >= 0);
	model.addConstr(d[11] - x[60] >= 0);
	model.addConstr(d[11] - x[61] >= 0);
	model.addConstr(d[11] - x[62] >= 0);
	model.addConstr(d[11] - x[63] >= 0);
	model.addConstr(x[48] + x[53] + x[58] + x[63] + x[64] + x[65] + x[66] + x[67] - 5 * d[12] >= 0);
	model.addConstr(d[12] - x[48] >= 0);
	model.addConstr(d[12] - x[53] >= 0);
	model.addConstr(d[12] - x[58] >= 0);
	model.addConstr(d[12] - x[63] >= 0);
	model.addConstr(d[12] - x[64] >= 0);
	model.addConstr(d[12] - x[65] >= 0);
	model.addConstr(d[12] - x[66] >= 0);
	model.addConstr(d[12] - x[67] >= 0);
	model.addConstr(x[52] + x[57] + x[62] + x[51] + x[68] + x[69] + x[70] + x[71] - 5 * d[13] >= 0);
	model.addConstr(d[13] - x[52] >= 0);
	model.addConstr(d[13] - x[57] >= 0);
	model.addConstr(d[13] - x[62] >= 0);
	model.addConstr(d[13] - x[51] >= 0);
	model.addConstr(d[13] - x[68] >= 0);
	model.addConstr(d[13] - x[69] >= 0);
	model.addConstr(d[13] - x[70] >= 0);
	model.addConstr(d[13] - x[71] >= 0);
	model.addConstr(x[56] + x[61] + x[50] + x[55] + x[72] + x[73] + x[74] + x[75] - 5 * d[14] >= 0);
	model.addConstr(d[14] - x[56] >= 0);
	model.addConstr(d[14] - x[61] >= 0);
	model.addConstr(d[14] - x[50] >= 0);
	model.addConstr(d[14] - x[55] >= 0);
	model.addConstr(d[14] - x[72] >= 0);
	model.addConstr(d[14] - x[73] >= 0);
	model.addConstr(d[14] - x[74] >= 0);
	model.addConstr(d[14] - x[75] >= 0);
	model.addConstr(x[60] + x[49] + x[54] + x[59] + x[76] + x[77] + x[78] + x[79] - 5 * d[15] >= 0);
	model.addConstr(d[15] - x[60] >= 0);
	model.addConstr(d[15] - x[49] >= 0);
	model.addConstr(d[15] - x[54] >= 0);
	model.addConstr(d[15] - x[59] >= 0);
	model.addConstr(d[15] - x[76] >= 0);
	model.addConstr(d[15] - x[77] >= 0);
	model.addConstr(d[15] - x[78] >= 0);
	model.addConstr(d[15] - x[79] >= 0);
	model.addConstr(x[64] + x[69] + x[74] + x[79] + x[80] + x[81] + x[82] + x[83] - 5 * d[16] >= 0);
	model.addConstr(d[16] - x[64] >= 0);
	model.addConstr(d[16] - x[69] >= 0);
	model.addConstr(d[16] - x[74] >= 0);
	model.addConstr(d[16] - x[79] >= 0);
	model.addConstr(d[16] - x[80] >= 0);
	model.addConstr(d[16] - x[81] >= 0);
	model.addConstr(d[16] - x[82] >= 0);
	model.addConstr(d[16] - x[83] >= 0);
	model.addConstr(x[68] + x[73] + x[78] + x[67] + x[84] + x[85] + x[86] + x[87] - 5 * d[17] >= 0);
	model.addConstr(d[17] - x[68] >= 0);
	model.addConstr(d[17] - x[73] >= 0);
	model.addConstr(d[17] - x[78] >= 0);
	model.addConstr(d[17] - x[67] >= 0);
	model.addConstr(d[17] - x[84] >= 0);
	model.addConstr(d[17] - x[85] >= 0);
	model.addConstr(d[17] - x[86] >= 0);
	model.addConstr(d[17] - x[87] >= 0);
	model.addConstr(x[72] + x[77] + x[66] + x[71] + x[88] + x[89] + x[90] + x[91] - 5 * d[18] >= 0);
	model.addConstr(d[18] - x[72] >= 0);
	model.addConstr(d[18] - x[77] >= 0);
	model.addConstr(d[18] - x[66] >= 0);
	model.addConstr(d[18] - x[71] >= 0);
	model.addConstr(d[18] - x[88] >= 0);
	model.addConstr(d[18] - x[89] >= 0);
	model.addConstr(d[18] - x[90] >= 0);
	model.addConstr(d[18] - x[91] >= 0);
	model.addConstr(x[76] + x[65] + x[70] + x[75] + x[92] + x[93] + x[94] + x[95] - 5 * d[19] >= 0);
	model.addConstr(d[19] - x[76] >= 0);
	model.addConstr(d[19] - x[65] >= 0);
	model.addConstr(d[19] - x[70] >= 0);
	model.addConstr(d[19] - x[75] >= 0);
	model.addConstr(d[19] - x[92] >= 0);
	model.addConstr(d[19] - x[93] >= 0);
	model.addConstr(d[19] - x[94] >= 0);
	model.addConstr(d[19] - x[95] >= 0);
	model.addConstr(x[80] + x[85] + x[90] + x[95] + x[96] + x[97] + x[98] + x[99] - 5 * d[20] >= 0);
	model.addConstr(d[20] - x[80] >= 0);
	model.addConstr(d[20] - x[85] >= 0);
	model.addConstr(d[20] - x[90] >= 0);
	model.addConstr(d[20] - x[95] >= 0);
	model.addConstr(d[20] - x[96] >= 0);
	model.addConstr(d[20] - x[97] >= 0);
	model.addConstr(d[20] - x[98] >= 0);
	model.addConstr(d[20] - x[99] >= 0);
	model.addConstr(x[84] + x[89] + x[94] + x[83] + x[100] + x[101] + x[102] + x[103] - 5 * d[21] >= 0);
	model.addConstr(d[21] - x[84] >= 0);
	model.addConstr(d[21] - x[89] >= 0);
	model.addConstr(d[21] - x[94] >= 0);
	model.addConstr(d[21] - x[83] >= 0);
	model.addConstr(d[21] - x[100] >= 0);
	model.addConstr(d[21] - x[101] >= 0);
	model.addConstr(d[21] - x[102] >= 0);
	model.addConstr(d[21] - x[103] >= 0);
	model.addConstr(x[88] + x[93] + x[82] + x[87] + x[104] + x[105] + x[106] + x[107] - 5 * d[22] >= 0);
	model.addConstr(d[22] - x[88] >= 0);
	model.addConstr(d[22] - x[93] >= 0);
	model.addConstr(d[22] - x[82] >= 0);
	model.addConstr(d[22] - x[87] >= 0);
	model.addConstr(d[22] - x[104] >= 0);
	model.addConstr(d[22] - x[105] >= 0);
	model.addConstr(d[22] - x[106] >= 0);
	model.addConstr(d[22] - x[107] >= 0);
	model.addConstr(x[92] + x[81] + x[86] + x[91] + x[108] + x[109] + x[110] + x[111] - 5 * d[23] >= 0);
	model.addConstr(d[23] - x[92] >= 0);
	model.addConstr(d[23] - x[81] >= 0);
	model.addConstr(d[23] - x[86] >= 0);
	model.addConstr(d[23] - x[91] >= 0);
	model.addConstr(d[23] - x[108] >= 0);
	model.addConstr(d[23] - x[109] >= 0);
	model.addConstr(d[23] - x[110] >= 0);
	model.addConstr(d[23] - x[111] >= 0);
	model.addConstr(x[96] + x[101] + x[106] + x[111] + x[112] + x[113] + x[114] + x[115] - 5 * d[24] >= 0);
	model.addConstr(d[24] - x[96] >= 0);
	model.addConstr(d[24] - x[101] >= 0);
	model.addConstr(d[24] - x[106] >= 0);
	model.addConstr(d[24] - x[111] >= 0);
	model.addConstr(d[24] - x[112] >= 0);
	model.addConstr(d[24] - x[113] >= 0);
	model.addConstr(d[24] - x[114] >= 0);
	model.addConstr(d[24] - x[115] >= 0);
	model.addConstr(x[100] + x[105] + x[110] + x[99] + x[116] + x[117] + x[118] + x[119] - 5 * d[25] >= 0);
	model.addConstr(d[25] - x[100] >= 0);
	model.addConstr(d[25] - x[105] >= 0);
	model.addConstr(d[25] - x[110] >= 0);
	model.addConstr(d[25] - x[99] >= 0);
	model.addConstr(d[25] - x[116] >= 0);
	model.addConstr(d[25] - x[117] >= 0);
	model.addConstr(d[25] - x[118] >= 0);
	model.addConstr(d[25] - x[119] >= 0);
	model.addConstr(x[104] + x[109] + x[98] + x[103] + x[120] + x[121] + x[122] + x[123] - 5 * d[26] >= 0);
	model.addConstr(d[26] - x[104] >= 0);
	model.addConstr(d[26] - x[109] >= 0);
	model.addConstr(d[26] - x[98] >= 0);
	model.addConstr(d[26] - x[103] >= 0);
	model.addConstr(d[26] - x[120] >= 0);
	model.addConstr(d[26] - x[121] >= 0);
	model.addConstr(d[26] - x[122] >= 0);
	model.addConstr(d[26] - x[123] >= 0);
	model.addConstr(x[108] + x[97] + x[102] + x[107] + x[124] + x[125] + x[126] + x[127] - 5 * d[27] >= 0);
	model.addConstr(d[27] - x[108] >= 0);
	model.addConstr(d[27] - x[97] >= 0);
	model.addConstr(d[27] - x[102] >= 0);
	model.addConstr(d[27] - x[107] >= 0);
	model.addConstr(d[27] - x[124] >= 0);
	model.addConstr(d[27] - x[125] >= 0);
	model.addConstr(d[27] - x[126] >= 0);
	model.addConstr(d[27] - x[127] >= 0);
	model.addConstr(x[112] + x[117] + x[122] + x[127] + x[128] + x[129] + x[130] + x[131] - 5 * d[28] >= 0);
	model.addConstr(d[28] - x[112] >= 0);
	model.addConstr(d[28] - x[117] >= 0);
	model.addConstr(d[28] - x[122] >= 0);
	model.addConstr(d[28] - x[127] >= 0);
	model.addConstr(d[28] - x[128] >= 0);
	model.addConstr(d[28] - x[129] >= 0);
	model.addConstr(d[28] - x[130] >= 0);
	model.addConstr(d[28] - x[131] >= 0);
	model.addConstr(x[116] + x[121] + x[126] + x[115] + x[132] + x[133] + x[134] + x[135] - 5 * d[29] >= 0);
	model.addConstr(d[29] - x[116] >= 0);
	model.addConstr(d[29] - x[121] >= 0);
	model.addConstr(d[29] - x[126] >= 0);
	model.addConstr(d[29] - x[115] >= 0);
	model.addConstr(d[29] - x[132] >= 0);
	model.addConstr(d[29] - x[133] >= 0);
	model.addConstr(d[29] - x[134] >= 0);
	model.addConstr(d[29] - x[135] >= 0);
	model.addConstr(x[120] + x[125] + x[114] + x[119] + x[136] + x[137] + x[138] + x[139] - 5 * d[30] >= 0);
	model.addConstr(d[30] - x[120] >= 0);
	model.addConstr(d[30] - x[125] >= 0);
	model.addConstr(d[30] - x[114] >= 0);
	model.addConstr(d[30] - x[119] >= 0);
	model.addConstr(d[30] - x[136] >= 0);
	model.addConstr(d[30] - x[137] >= 0);
	model.addConstr(d[30] - x[138] >= 0);
	model.addConstr(d[30] - x[139] >= 0);
	model.addConstr(x[124] + x[113] + x[118] + x[123] + x[140] + x[141] + x[142] + x[143] - 5 * d[31] >= 0);
	model.addConstr(d[31] - x[124] >= 0);
	model.addConstr(d[31] - x[113] >= 0);
	model.addConstr(d[31] - x[118] >= 0);
	model.addConstr(d[31] - x[123] >= 0);
	model.addConstr(d[31] - x[140] >= 0);
	model.addConstr(d[31] - x[141] >= 0);
	model.addConstr(d[31] - x[142] >= 0);
	model.addConstr(d[31] - x[143] >= 0);
	model.addConstr(x[128] + x[133] + x[138] + x[143] + x[144] + x[145] + x[146] + x[147] - 5 * d[32] >= 0);
	model.addConstr(d[32] - x[128] >= 0);
	model.addConstr(d[32] - x[133] >= 0);
	model.addConstr(d[32] - x[138] >= 0);
	model.addConstr(d[32] - x[143] >= 0);
	model.addConstr(d[32] - x[144] >= 0);
	model.addConstr(d[32] - x[145] >= 0);
	model.addConstr(d[32] - x[146] >= 0);
	model.addConstr(d[32] - x[147] >= 0);
	model.addConstr(x[132] + x[137] + x[142] + x[131] + x[148] + x[149] + x[150] + x[151] - 5 * d[33] >= 0);
	model.addConstr(d[33] - x[132] >= 0);
	model.addConstr(d[33] - x[137] >= 0);
	model.addConstr(d[33] - x[142] >= 0);
	model.addConstr(d[33] - x[131] >= 0);
	model.addConstr(d[33] - x[148] >= 0);
	model.addConstr(d[33] - x[149] >= 0);
	model.addConstr(d[33] - x[150] >= 0);
	model.addConstr(d[33] - x[151] >= 0);
	model.addConstr(x[136] + x[141] + x[130] + x[135] + x[152] + x[153] + x[154] + x[155] - 5 * d[34] >= 0);
	model.addConstr(d[34] - x[136] >= 0);
	model.addConstr(d[34] - x[141] >= 0);
	model.addConstr(d[34] - x[130] >= 0);
	model.addConstr(d[34] - x[135] >= 0);
	model.addConstr(d[34] - x[152] >= 0);
	model.addConstr(d[34] - x[153] >= 0);
	model.addConstr(d[34] - x[154] >= 0);
	model.addConstr(d[34] - x[155] >= 0);
	model.addConstr(x[140] + x[129] + x[134] + x[139] + x[156] + x[157] + x[158] + x[159] - 5 * d[35] >= 0);
	model.addConstr(d[35] - x[140] >= 0);
	model.addConstr(d[35] - x[129] >= 0);
	model.addConstr(d[35] - x[134] >= 0);
	model.addConstr(d[35] - x[139] >= 0);
	model.addConstr(d[35] - x[156] >= 0);
	model.addConstr(d[35] - x[157] >= 0);
	model.addConstr(d[35] - x[158] >= 0);
	model.addConstr(d[35] - x[159] >= 0);
	model.addConstr(x[144] + x[149] + x[154] + x[159] + x[160] + x[161] + x[162] + x[163] - 5 * d[36] >= 0);
	model.addConstr(d[36] - x[144] >= 0);
	model.addConstr(d[36] - x[149] >= 0);
	model.addConstr(d[36] - x[154] >= 0);
	model.addConstr(d[36] - x[159] >= 0);
	model.addConstr(d[36] - x[160] >= 0);
	model.addConstr(d[36] - x[161] >= 0);
	model.addConstr(d[36] - x[162] >= 0);
	model.addConstr(d[36] - x[163] >= 0);
	model.addConstr(x[148] + x[153] + x[158] + x[147] + x[164] + x[165] + x[166] + x[167] - 5 * d[37] >= 0);
	model.addConstr(d[37] - x[148] >= 0);
	model.addConstr(d[37] - x[153] >= 0);
	model.addConstr(d[37] - x[158] >= 0);
	model.addConstr(d[37] - x[147] >= 0);
	model.addConstr(d[37] - x[164] >= 0);
	model.addConstr(d[37] - x[165] >= 0);
	model.addConstr(d[37] - x[166] >= 0);
	model.addConstr(d[37] - x[167] >= 0);
	model.addConstr(x[152] + x[157] + x[146] + x[151] + x[168] + x[169] + x[170] + x[171] - 5 * d[38] >= 0);
	model.addConstr(d[38] - x[152] >= 0);
	model.addConstr(d[38] - x[157] >= 0);
	model.addConstr(d[38] - x[146] >= 0);
	model.addConstr(d[38] - x[151] >= 0);
	model.addConstr(d[38] - x[168] >= 0);
	model.addConstr(d[38] - x[169] >= 0);
	model.addConstr(d[38] - x[170] >= 0);
	model.addConstr(d[38] - x[171] >= 0);
	model.addConstr(x[156] + x[145] + x[150] + x[155] + x[172] + x[173] + x[174] + x[175] - 5 * d[39] >= 0);
	model.addConstr(d[39] - x[156] >= 0);
	model.addConstr(d[39] - x[145] >= 0);
	model.addConstr(d[39] - x[150] >= 0);
	model.addConstr(d[39] - x[155] >= 0);
	model.addConstr(d[39] - x[172] >= 0);
	model.addConstr(d[39] - x[173] >= 0);
	model.addConstr(d[39] - x[174] >= 0);
	model.addConstr(d[39] - x[175] >= 0);
	model.addConstr(x[160] + x[165] + x[170] + x[175] + x[176] + x[177] + x[178] + x[179] - 5 * d[40] >= 0);
	model.addConstr(d[40] - x[160] >= 0);
	model.addConstr(d[40] - x[165] >= 0);
	model.addConstr(d[40] - x[170] >= 0);
	model.addConstr(d[40] - x[175] >= 0);
	model.addConstr(d[40] - x[176] >= 0);
	model.addConstr(d[40] - x[177] >= 0);
	model.addConstr(d[40] - x[178] >= 0);
	model.addConstr(d[40] - x[179] >= 0);
	model.addConstr(x[164] + x[169] + x[174] + x[163] + x[180] + x[181] + x[182] + x[183] - 5 * d[41] >= 0);
	model.addConstr(d[41] - x[164] >= 0);
	model.addConstr(d[41] - x[169] >= 0);
	model.addConstr(d[41] - x[174] >= 0);
	model.addConstr(d[41] - x[163] >= 0);
	model.addConstr(d[41] - x[180] >= 0);
	model.addConstr(d[41] - x[181] >= 0);
	model.addConstr(d[41] - x[182] >= 0);
	model.addConstr(d[41] - x[183] >= 0);
	model.addConstr(x[168] + x[173] + x[162] + x[167] + x[184] + x[185] + x[186] + x[187] - 5 * d[42] >= 0);
	model.addConstr(d[42] - x[168] >= 0);
	model.addConstr(d[42] - x[173] >= 0);
	model.addConstr(d[42] - x[162] >= 0);
	model.addConstr(d[42] - x[167] >= 0);
	model.addConstr(d[42] - x[184] >= 0);
	model.addConstr(d[42] - x[185] >= 0);
	model.addConstr(d[42] - x[186] >= 0);
	model.addConstr(d[42] - x[187] >= 0);
	model.addConstr(x[172] + x[161] + x[166] + x[171] + x[188] + x[189] + x[190] + x[191] - 5 * d[43] >= 0);
	model.addConstr(d[43] - x[172] >= 0);
	model.addConstr(d[43] - x[161] >= 0);
	model.addConstr(d[43] - x[166] >= 0);
	model.addConstr(d[43] - x[171] >= 0);
	model.addConstr(d[43] - x[188] >= 0);
	model.addConstr(d[43] - x[189] >= 0);
	model.addConstr(d[43] - x[190] >= 0);
	model.addConstr(d[43] - x[191] >= 0);
	model.addConstr(x[176] + x[181] + x[186] + x[191] + x[192] + x[193] + x[194] + x[195] - 5 * d[44] >= 0);
	model.addConstr(d[44] - x[176] >= 0);
	model.addConstr(d[44] - x[181] >= 0);
	model.addConstr(d[44] - x[186] >= 0);
	model.addConstr(d[44] - x[191] >= 0);
	model.addConstr(d[44] - x[192] >= 0);
	model.addConstr(d[44] - x[193] >= 0);
	model.addConstr(d[44] - x[194] >= 0);
	model.addConstr(d[44] - x[195] >= 0);
	model.addConstr(x[180] + x[185] + x[190] + x[179] + x[196] + x[197] + x[198] + x[199] - 5 * d[45] >= 0);
	model.addConstr(d[45] - x[180] >= 0);
	model.addConstr(d[45] - x[185] >= 0);
	model.addConstr(d[45] - x[190] >= 0);
	model.addConstr(d[45] - x[179] >= 0);
	model.addConstr(d[45] - x[196] >= 0);
	model.addConstr(d[45] - x[197] >= 0);
	model.addConstr(d[45] - x[198] >= 0);
	model.addConstr(d[45] - x[199] >= 0);
	model.addConstr(x[184] + x[189] + x[178] + x[183] + x[200] + x[201] + x[202] + x[203] - 5 * d[46] >= 0);
	model.addConstr(d[46] - x[184] >= 0);
	model.addConstr(d[46] - x[189] >= 0);
	model.addConstr(d[46] - x[178] >= 0);
	model.addConstr(d[46] - x[183] >= 0);
	model.addConstr(d[46] - x[200] >= 0);
	model.addConstr(d[46] - x[201] >= 0);
	model.addConstr(d[46] - x[202] >= 0);
	model.addConstr(d[46] - x[203] >= 0);
	model.addConstr(x[188] + x[177] + x[182] + x[187] + x[204] + x[205] + x[206] + x[207] - 5 * d[47] >= 0);
	model.addConstr(d[47] - x[188] >= 0);
	model.addConstr(d[47] - x[177] >= 0);
	model.addConstr(d[47] - x[182] >= 0);
	model.addConstr(d[47] - x[187] >= 0);
	model.addConstr(d[47] - x[204] >= 0);
	model.addConstr(d[47] - x[205] >= 0);
	model.addConstr(d[47] - x[206] >= 0);
	model.addConstr(d[47] - x[207] >= 0);
	model.addConstr(x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] + x[11] + x[12] + x[13] + x[14] + x[15] + x[16] + x[17] + x[18] + x[19] + x[20] + x[21] + x[22] + x[23] + x[24] + x[25] + x[26] + x[27] + x[28] + x[29] + x[30] + x[31] + x[32] + x[33] + x[34] + x[35] + x[36] + x[37] + x[38] + x[39] + x[40] + x[41] + x[42] + x[43] + x[44] + x[45] + x[46] + x[47] + x[48] + x[49] + x[50] + x[51] + x[52] + x[53] + x[54] + x[55] + x[56] + x[57] + x[58] + x[59] + x[60] + x[61] + x[62] + x[63] + x[64] + x[65] + x[66] + x[67] + x[68] + x[69] + x[70] + x[71] + x[72] + x[73] + x[74] + x[75] + x[76] + x[77] + x[78] + x[79] + x[80] + x[81] + x[82] + x[83] + x[84] + x[85] + x[86] + x[87] + x[88] + x[89] + x[90] + x[91] + x[92] + x[93] + x[94] + x[95] + x[96] + x[97] + x[98] + x[99] + x[100] + x[101] + x[102] + x[103] + x[104] + x[105] + x[106] + x[107] + x[108] + x[109] + x[110] + x[111] + x[112] + x[113] + x[114] + x[115] + x[116] + x[117] + x[118] + x[119] + x[120] + x[121] + x[122] + x[123] + x[124] + x[125] + x[126] + x[127] + x[128] + x[129] + x[130] + x[131] + x[132] + x[133] + x[134] + x[135] + x[136] + x[137] + x[138] + x[139] + x[140] + x[141] + x[142] + x[143] + x[144] + x[145] + x[146] + x[147] + x[148] + x[149] + x[150] + x[151] + x[152] + x[153] + x[154] + x[155] + x[156] + x[157] + x[158] + x[159] + x[160] + x[161] + x[162] + x[163] + x[164] + x[165] + x[166] + x[167] + x[168] + x[169] + x[170] + x[171] + x[172] + x[173] + x[174] + x[175] + x[176] + x[177] + x[178] + x[179] + x[180] + x[181] + x[182] + x[183] + x[184] + x[185] + x[186] + x[187] + x[188] + x[189] + x[190] + x[191] >= 1);



	GRBLinExpr val = 0;

	int round = 10;

	for (int i = 16; i < 16 * round - 16; i++)
	{
		val += 6 * x[i];
	}


	for (int i = 0; i < 16; i++)
	{
		val += 3.83 * x[i];
	}
	for (int i = 16 * round - 16; i < 16 * round; i++)
	{
		val += 3.83 * x[i];
	}

	model.setObjective(val, GRB_MINIMIZE);
	model.update();
	model.optimize();

	cout << "GRB_MINIMIZE: " << val.getValue() << std::endl;

	for (int r = 0; r < round; r++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << fabs(x[r * 16 + i + 4 * j].get(GRB_DoubleAttr_Xn));
			}
			cout << endl;
		}
		cout << endl << endl;
	}

	cout << endl;


	return 0;
}




