const int nSTARcentPoints =9;
const int nSTARpTPoints =8;
const int nSTARetaPoints =20;
const int nSTARsystems = 2;

// http://drupal.star.bnl.gov/STAR/files/starpublications/120/data.html

// from Fig. 3 http://prl.aps.org/abstract/PRL/v101/i25/e252301
// "http://drupal.star.bnl.gov/STAR/files/starpublications/120/fig2.html"
float v1_AuAu_200GeV_pT_5_40[nSTARetaPoints][5]=
{
0.226, 0.074, 0.076, -0.096, 0.011,
0.394, 0.106, 0.094, -0.136, 0.010,
0.592, 0.108, 0.92, -0.149, 0.012,
0.832, 0.168, 0.132, -0.152, 0.013,
1.168, 0.232, 0.168, -0.096, 0.018,
1.628, 0.372, 0.228, -0.007, 0.026,
2.278, 0.522, 0.278, 0.237, 0.052,
3.170, 0.830, 0.370, 0.311, 0.127
};

float v1_AuAu_200GeV_pT_40_80[nSTARetaPoints][5]=
{
0.230, 0.070, 0.080, -0.114, 0.016,
0.393, 0.107, 0.093, -0.193, 0.015,
0.591, 0.109, 0.91, -0.244, 0.020,
0.830, 0.170, 0.130, -0.246, 0.022,
1.166, 0.234, 0.166, -0.212, 0.030,
1.626, 0.374, 0.226, -0.133, 0.046,
2.281, 0.519, 0.281, -0.096, 0.091,
3.180, 0.820, 0.380, -0.154, 0.211
};

// from Fig. 3 http://prl.aps.org/abstract/PRL/v101/i25/e252301
// "http://drupal.star.bnl.gov/STAR/files/starpublications/120/fig3.html"
float v1_AuAu_200GeV_eta_30_60[nSTARetaPoints][5]=
{
-3.82, 0.12, 0.18, 2.898, 0.048,
-3.54, 0.14, 0.16, 2.465, 0.033,
-3.30, 0.10, 0.10, 2.076, 0.036,
-3.10, 0.10, 0.10, 1.894, 0.034,
-2.87, 0.17, 0.13, 1.644, 0.030,
-2.64, 0.14, 0.06, 1.417, 0.102,
-1.03, 0.13, 0.27, 0.385, 0.022,
-0.75, 0.15, 0.15, 0.216, 0.018,
-0.45, 0.15, 0.15, 0.101, 0.018,
-0.15, 0.15, 0.15, 0.076, 0.018,
0.15, 0.15, 0.15, -0.010, 0.018,
0.45, 0.15, 0.15, -0.124, 0.018,
0.75, 0.15, 0.15, -0.184, 0.018,
1.03, 0.27, 0.13, -0.330, 0.022,
2.64, 0.06, 0.14, -1.669, 0.103,
2.87, 0.13, 0.17, -1.679, 0.032,
3.10, 0.10, 0.10, -1.837, 0.036,
3.30, 0.10, 0.10, -2.128, 0.038,
3.54, 0.16, 0.14, -2.445, 0.035,
3.82, 0.18, 0.12, -2.611, 0.051
};

float v1_CuCu_200GeV_eta_30_60[12][5]=
{
-3.76, 0.16, 0.24, 2.869, 0.099,
-3.39, 0.19, 0.21, 2.220, 0.071,
-3.00, 0.20, 0.20, 1.663, 0.064,
-2.71, 0.21, 0.09, 1.177, 0.128,
-0.85, 0.25, 0.45, 0.233, 0.040,
-0.30, 0.30, 0.30, 0.167, 0.036,
0.30, 0.30, 0.30, -0.036, 0.036,
0.85, 0.45, 0.25, -0.224, 0.040,
2.71, 0.09, 0.21, -1.381, 0.139,
3.00, 0.20, 0.20, -1.587, 0.067,
3.39, 0.21, 0.19, -2.072, 0.073,
3.76, 0.24, 0.16, -2.672, 0.100
};

float v1_AuAu_62GeV_eta_30_60[nSTARetaPoints][5]=
{
-3.82, 0.12, 0.18, 5.656, 0.244,
-3.54, 0.14, 0.16, 5.572, 0.166,
-3.30, 0.10, 0.10, 5.329, 0.176,
-3.10, 0.10, 0.10, 4.699, 0.162,
-2.87, 0.17, 0.13, 3.824, 0.142,
-2.64, 0.14, 0.06, 2.781, 0.473,
-1.02, 0.12, 0.28, 0.925, 0.087,
-0.75, 0.15, 0.15, 0.669, 0.069,
-0.45, 0.15, 0.15, 0.317, 0.068,
-0.15, 0.15, 0.15, 0.122, 0.067,
0.15, 0.15, 0.15, -0.118, 0.067,
0.45, 0.15, 0.15, -0.357, 0.067,
0.75, 0.15, 0.15, -0.634, 0.068,
1.02, 0.28, 0.12, -0.882, 0.086,
2.64, 0.06, 0.14, -3.453, 0.328,
2.87, 0.13, 0.17, -3.986, 0.122,
3.10, 0.10, 0.10, -4.602, 0.149,
3.30, 0.10, 0.10, -5.176, 0.161,
3.54, 0.16, 0.14, -5.888, 0.156,
3.82, 0.18, 0.12, -6.434, 0.233
};

float v1_CuCu_62GeV_eta_30_60[16][5]=
{
-3.76, 0.16, 0.24, 8.914, 0.853,
-3.44, 0.14, 0.16, 5.966, 0.667,
-3.15, 0.15, 0.15, 6.035, 0.619,
-2.82, 0.32, 0.18, 4.218, 0.493,
-1.04, 0.14, 0.26, 1.027, 0.359,
-0.75, 0.15, 0.15, 0.500, 0.307,
-0.45, 0.15, 0.15, 0.363, 0.306,
-0.15, 0.15, 0.15, 0.282, 0.307,
0.15, 0.15, 0.15, -0.174, 0.306,
0.45, 0.15, 0.15, -0.431, 0.316,
0.75, 0.15, 0.15, -0.466, 0.314,
1.04, 0.26, 0.14, -1.205, 0.381,
2.82, 0.18, 0.32, -3.454, 0.545,
3.15, 0.15, 0.15, -5.165, 0.642,
3.44, 0.16, 0.14, -6.084, 0.742,
3.76, 0.24, 0.16, -6.599, 0.872
};
// from Fig. 4 http://prl.aps.org/abstract/PRL/v101/i25/e252301
// "http://drupal.star.bnl.gov/STAR/files/starpublications/120/Fig4.html"

float star_data_AuAu200_cent[nSTARcentPoints][5]=
{
  //Au+Au@200
  2.5, -0.038, 0.026, -0.642, 0.059,
  7.5, -0.084, 0.015, -0.794, 0.033,
  15, -0.130, 0.009, -1.104, 0.018,
  25, -0.123, 0.008, -1.511, 0.016,
  35, -0.153, 0.010, -1.894, 0.017,
  45, -0.160, 0.012, -2.284, 0.020,
  55, -0.192, 0.015, -2.610, 0.025,
  65, -0.191, 0.034, -2.904, 0.036,
  75, -0.206, 0.038, -3.061, 0.059
};
float star_data_AuAu62_cent[nSTARcentPoints][5]=
{
  //Au+Au@62
  2.5, -999,0,0,0,
  7.5, -999,0,0,0,
  15, -0.280, 0.051, -2.173, 0.117,
  25, -0.352, 0.038, -3.010, 0.082,
  35, -0.463, 0.035, -3.978, 0.071,
  45, -0.467, 0.038, -5.089, 0.074,
  55, -0.468, 0.044, -6.303, 0.083,
  65, -0.553, 0.058, -7.558, 0.105,
  75, -0.531, 0.086, -8.554, 0.152  
};
TGraphErrors *v1_star_AuAu62_cent;
TGraphErrors *v1_star_AuAu200_cent;
TGraphErrors *v1_star_AuAu62_30_60_eta;
TGraphErrors *v1_star_AuAu200_30_60_eta;
TGraphErrors *v1_star_AuAu200_5_40_pT;
TGraphErrors *v1_star_AuAu200_40_80_pT;

TGraphErrors *star_dataGraph_cent[nSTARsystems];
TGraphErrors *star_dataGraph_eta[nSTARsystems];
TGraphErrors *star_dataGraph_pT[nSTARsystems];

void directedFlow_2008_STARdataPt(float nSTARscale1 =0.01,float nSTARscale2 =0.01)
{
  float star_dataX[nSTARsystems][nSTARpTPoints];
  float star_dataY[nSTARsystems][nSTARpTPoints];
  float star_dataXerr[nSTARsystems][nSTARpTPoints];
  float star_dataYerr[nSTARsystems][nSTARpTPoints];
  for (int i=0; i < nSTARpTPoints; ++i)
  {
    star_dataX[0][i] = v1_AuAu_200GeV_pT_5_40[i][0];
    star_dataXerr[0][i] = 0;
    star_dataY[0][i] = nSTARscale1*v1_AuAu_200GeV_pT_5_40[i][3];
    star_dataYerr[0][i] = nSTARscale1*v1_AuAu_200GeV_pT_5_40[i][4];
    
    star_dataX[1][i] = v1_AuAu_200GeV_pT_40_80[i][0];
    star_dataXerr[1][i] = 0;
    star_dataY[1][i] = nSTARscale2*v1_AuAu_200GeV_pT_40_80[i][3];
    star_dataYerr[1][i] = nSTARscale2*v1_AuAu_200GeV_pT_40_80[i][4];
  }
  for (int j=0; j < nSTARsystems; ++j)
  {
    star_dataGraph_pT[j]= new TGraphErrors(nSTARpTPoints,star_dataX[j],star_dataY[j],star_dataXerr[j],star_dataYerr[j]);
    star_dataGraph_pT[j]->SetLineColor(j+1);
    star_dataGraph_pT[j]->SetLineWidth(3);
    star_dataGraph_pT[j]->SetLineStyle(3);
    star_dataGraph_pT[j]->SetMarkerColor(j+1);
    star_dataGraph_pT[j]->SetMarkerStyle(24+j);
    star_dataGraph_pT[j]->SetMarkerSize(3);
  v1_star_AuAu200_5_40_pT = star_dataGraph_pT[0];
  v1_star_AuAu200_40_80_pT = star_dataGraph_pT[1];
  }
}

void directedFlow_2008_STARdataCentrality(float nSTARscale1 =0.01,float nSTARscale2 =0.01)
{
  float star_dataX[nSTARsystems][nSTARcentPoints];
  float star_dataY[nSTARsystems][nSTARcentPoints];
  float star_dataXerr[nSTARsystems][nSTARcentPoints];
  float star_dataYerr[nSTARsystems][nSTARcentPoints];
  for (int i=0; i < nSTARcentPoints; ++i)
  {
    star_dataX[0][i] = star_data_AuAu200_cent[i][0];
    star_dataXerr[0][i] = 0;
    star_dataY[0][i] = nSTARscale1*star_data_AuAu200_cent[i][1];
    star_dataYerr[0][i] = nSTARscale1*star_data_AuAu200_cent[i][2];
    
    star_dataX[1][i] = star_data_AuAu62_cent[i][0];
    star_dataXerr[1][i] = 0;
    star_dataY[1][i] = nSTARscale2*star_data_AuAu62_cent[i][1];
    star_dataYerr[1][i] = nSTARscale2*star_data_AuAu62_cent[i][2];
  }
  for (int j=0; j < nSTARsystems; ++j)
  {
    star_dataGraph_cent[j]= new TGraphErrors(nSTARcentPoints,star_dataX[j],star_dataY[j],star_dataXerr[j],star_dataYerr[j]);
    star_dataGraph_cent[j]->SetLineColor(j+1);
    star_dataGraph_cent[j]->SetLineWidth(3);
    star_dataGraph_cent[j]->SetLineStyle(3);
    star_dataGraph_cent[j]->SetMarkerColor(j+1);
    star_dataGraph_cent[j]->SetMarkerStyle(24+j);
    star_dataGraph_cent[j]->SetMarkerSize(3);
  }
  v1_star_AuAu200_cent = star_dataGraph_cent[0];
  v1_star_AuAu62_cent = star_dataGraph_cent[1];
}

void directedFlow_2008_STARdataEta(float nSTARscale1 =0.01,float nSTARscale2 =0.01)
{
  float star_dataX[nSTARsystems][nSTARetaPoints];
  float star_dataY[nSTARsystems][nSTARetaPoints];
  float star_dataXerr[nSTARsystems][nSTARetaPoints];
  float star_dataYerr[nSTARsystems][nSTARetaPoints];
  for (int i=0; i < nSTARetaPoints; ++i)
  {
    star_dataX[0][i] = v1_AuAu_200GeV_eta_30_60[i][0];
    star_dataXerr[0][i] = 0;
    star_dataY[0][i] = nSTARscale1*v1_AuAu_200GeV_eta_30_60[i][3];
    star_dataYerr[0][i] = nSTARscale1*v1_AuAu_200GeV_eta_30_60[i][4];
    
    star_dataX[1][i] = v1_AuAu_62GeV_eta_30_60[i][0];
    star_dataXerr[1][i] = 0;
    star_dataY[1][i] = nSTARscale2*v1_AuAu_62GeV_eta_30_60[i][3];
    star_dataYerr[1][i] = nSTARscale2*v1_AuAu_62GeV_eta_30_60[i][4];
  }
  for (int j=0; j < nSTARsystems; ++j)
  {
    star_dataGraph_eta[j]= new TGraphErrors(nSTARetaPoints,star_dataX[j],star_dataY[j],star_dataXerr[j],star_dataYerr[j]);
    star_dataGraph_eta[j]->SetLineColor(j+1);
    star_dataGraph_eta[j]->SetLineWidth(3);
    star_dataGraph_eta[j]->SetLineStyle(3);
    star_dataGraph_eta[j]->SetMarkerColor(j+1);
    star_dataGraph_eta[j]->SetMarkerStyle(24+j);
    star_dataGraph_eta[j]->SetMarkerSize(3);
  v1_star_AuAu200_30_60_eta = star_dataGraph_eta[0];
  v1_star_AuAu62_30_60_eta = star_dataGraph_eta[1];
  }
}


