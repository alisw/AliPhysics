// Merged original figures 2a and 2b into one figure. 

void figure2()
{
 // Set style:
 SetFlowStyle();


  Int_t  kDarkGreen = TColor::GetColor(34,139,34);
  Int_t  kDarkCyan = TColor::GetColor(0,180,180);



  const int nbins = 17;

  //-----70-80% ---------------------------------------------------------------------------

  double pt_7080[nbins] = {
    0.29644,0.475364,0.700708,0.895582,1.07612,
    1.30166,1.4961,1.67731,1.90243,2.18302,2.56019,2.96531,
    3.39006,3.93562,5.10687,7.32286, 10.2651
  };
  double errPtLow_7080[nbins] = {
    0.0964398,0.075364,0.100708,0.0955816,
    0.076121,0.101662,0.0960967,0.0773092,0.102431,
    0.183019,0.160193,0.165308,0.19006,0.335623,0.606867,0.822858,0.765051
  };
  double errPtHigh_7080[nbins] = {
    0.10356,0.124636,0.0992923,0.104418,
    0.123879,0.0983377,0.103903,0.122691,0.097569,0.216981,
    0.239806,0.234692,0.209939,0.564377,1.39313,2.17714,1.73495
  };
  double v2_7080[nbins] ={
    6.18034,7.78103,9.58698,16.708,16.1315,14.3096,
    11.8754,18.648,5.1453,6.53926,-6.88641,-13.9267,-52.5485,
    -191.067,-415.917,291.021,45.6615
  };
  double errV2_7080[nbins] = {
    0.566557,0.582667,0.830618,1.38276,1.49568,
    2.06431,3.40244,3.51732,4.83874,5.55657,8.97981,
    13.9386,22.9194,31.4821,76.8651,-9999,-9999
  };

  //-----60-70% ----------------------------------------------------------------------------------

  double pt_6070[nbins] = {
    0.296676,0.47597,0.701133,0.895742,1.07645,
    1.30194,1.49616,1.67738,1.90256,2.18364,2.56013,2.96414,
    3.38945,3.93604,5.09071,7.3686,10.683
  };		    
  double errPtLow_6070[nbins] = {
    0.0966761,0.0759702,0.101133,0.0957418,
    0.0764486,0.101935,0.0961564,0.0773785,0.102563,0.183642,
    0.160132,0.164139,0.189453,0.33604,0.59071,0.868602,1.18298
  };
  double errPtHigh_6070[nbins] = {
    0.103324,0.12403,0.0988673,0.104258,
    0.123551,0.0980649,0.103844,0.122621,0.0974373,
    0.216358,0.239868,0.235861,0.210547,0.563961,1.40929,2.1314,1.31702
  };
  double v2_6070[nbins] = {
    3.49974,5.88458,8.12386,10.0714,12.5145,11.6577,
    11.5439,13.5126,11.4789,10.7059,11.6547,9.49895,
    1.37392,5.47184,-1.3349,-45.6053,156.746
  };
  double errV2_6070[nbins] = {
    0.177159,0.179905,0.250419,0.406175,0.431074,
    0.595281,0.940248,0.987307,1.32458,1.50529,2.40035,
    3.48782,6.0306,7.22813,13.5038,56.712,-9999
  };

 
 //-----50-60% -----------------------------------------------------------------------------------

  double pt_5060[nbins] = {
    0.296932,0.47652,0.701555,0.895881,1.07678,1.30221,
    1.49621,1.67759,1.90265,2.18367,2.55994,2.96436,
    3.38993,3.93277,5.08075,7.40325,10.5265
  };		    
  double errPtLow_5060[nbins] = {
    0.0969318,0.0765199,0.101555,0.0958806,
    0.0767766,0.102209,0.0962076,0.0775871,0.102653,0.183672,
    0.159938,0.164359,0.189935,0.332773,0.580747,0.903252,1.02646
  };
  double errPtHigh_5060[nbins] = {
    0.103068,0.12348,0.0984446,0.104119,
    0.123223,0.097791,0.103792,0.122413,0.0973468,0.216328,
    0.240062,0.235642,0.210065,0.567227,1.41925,2.09675,1.47354
  };
  double v2_5060[nbins] = {
    3.44798,5.90898,8.4642,10.5149,11.9192,13.3181,
    14.0252,14.3745,15.0161,14.8689,13.5232,13.4977,11.6071,
    9.47076,7.71149,17.5352,-5.54787
  };
  double errV2_5060[nbins] = {
    0.0606734,0.0606551,0.0829413,0.132536,
    0.139612,0.190289,0.299361,0.317186,0.414861,0.475373,
    0.769568,1.05829,1.91803,2.09282,3.91973,13.8477,40.6673
  };
  
//-----40-50% ----------------------------------------------------------------------------------

  double pt_4050[nbins] = {
    0.297194,0.47704,0.701961,0.896023,1.07714,1.30245,
    1.4963,1.6777,1.90274,2.18372,2.55965,2.96361,3.38975,
    3.93175,5.08329,7.35969,10.4992
  };
  double errPtLow_4050[nbins] = {
    0.0971937,0.0770399,0.101961,0.0960231,0.0771427,
    0.102452,0.0963024,0.077695,0.10274,0.183724,
    0.159645,0.16361,0.189751,0.331754,0.583295,0.859692,0.999224
  };
  double errPtHigh_4050[nbins] = {
    0.102806,0.12296,0.0980395,0.103977,0.122857,0.0975482,
    0.103698,0.122305,0.0972599,0.216276,0.240355,
    0.23639,0.210248,0.568246,1.41671,2.14031,1.50078
  };
  double v2_4050[nbins] = {
    3.64534,5.97555,8.71339,10.8948,12.3118,14.0208,15.1841,
    16.2186,17.0478,17.4933,17.837,17.5645,16.0996,14.2722,
    11.2284,14.8254,39.8013
  };
  double errV2_4050[nbins] = {
    0.0280054,0.0275449,0.0370108,0.0584979,0.0610596,
    0.0821168,0.129148,0.133787,0.179269,0.201004,0.319255,
    0.465269,0.805499,0.917085,1.69463,5.79271,-9999
  };

//-----30-40% -----------------------------------------------------------------------------------

  double pt_3040[nbins] = {
    0.297459,0.477509,0.702306,0.896152,1.07742,1.30267,
    1.49633,1.67775,1.90278,2.18368,2.55957,2.96322,3.38921,
    3.92869,5.07593,7.41037,10.408
  };
  double errPtLow_3040[nbins] = {
    0.0974586,0.0775094,0.102306,0.0961525,
    0.0774162,0.102673,0.096333,0.077749,0.10278,0.183675,
    0.159572,0.163219,0.189209,0.32869,0.575929,0.910369,0.907955
  };
  double errPtHigh_3040[nbins] = {
    0.102541,0.122491,0.0976941,0.103848,
    0.122584,0.0973271,0.103667,0.122251,0.0972203,0.216325,
    0.240428,0.236781,0.210791,0.57131,1.42407,2.08963,1.59205
  };
  double v2_3040[nbins] = {
    3.32144,5.60995,8.20984,10.1864,11.7526,13.5713,
    14.8011,15.8389,16.7527,17.5895,18.1343,18.2392,18.4816,16.8442,
    14.3587,7.07376,-26.4391
  };
  double errV2_3040[nbins] = {
    0.0171702,0.0166421,0.0220506,0.034594,0.0358166,
    0.0478107,0.0746667,0.0773614,0.103065,0.115688,
    0.184685,0.270488,0.46637,0.53962,1.01023,3.51725,12.6064
  };


//-----20-30% -----------------------------------------------------------------------------------

  double pt_2030[nbins] = {
    0.297722,0.477941,0.702611,0.896268,1.07768,1.30284,
    1.49638,1.67786,1.9028,2.18362,2.55924,2.96217,3.38881,
    3.92699,5.07304,7.42509,10.5134
  };
  double errPtLow_2030[nbins] = {
    0.0977221,0.0779407,0.102611,0.0962682,
    0.0776817,0.102839,0.0963758,0.0778584,0.102797,
    0.183618,0.159242,0.162174,0.188808,0.326986,0.573036,0.925089,1.01339
  };
  double errPtHigh_2030[nbins] = {
    0.102278,0.122059,0.0973894,0.103732,
    0.122318,0.0971605,0.103624,0.122142,0.0972035,0.216382,
    0.240758,0.237826,0.211192,0.573015,1.42696,2.07491,1.48661
  };
  double v2_2030[nbins] = {
    2.95888,4.90692,7.19851,8.89996,10.2714,11.926,13.1472,
    14.1925,15.2307,16.2584,16.8146,16.351,16.1069,14.8237,13.7406,
    6.26249,-2.67304
  };
  double errV2_2030[nbins] = {
    0.0134584,0.0128658,0.016865,0.0262622,
    0.0270303,0.03587,0.0557706,0.0576958,0.0765681,0.0859164,
    0.13768,0.203798,0.357643,0.419151,0.826191,2.70724,8.88004
  };
 
 //-----10-20% -----------------------------------------------------------------------------------

  double pt_1020[nbins] = {
    0.298001,0.478337,0.702859,0.896349,1.07787,1.30299,
    1.49642,1.6779,1.90278,2.18349,2.55915,2.96176,3.38798,
    3.92357,5.07047,7.40607,10.4205
  };
  double errPtLow_1020[nbins] = {
    0.0980014,0.0783374,0.102859,0.0963492,
    0.0778747,0.102987,0.0964151,0.0779025,0.102782,0.183491,
    0.159151,0.161757,0.187976,0.323571,0.57047,0.906065,0.920542
  };
  double errPtHigh_1020[nbins] = {
    0.101999,0.121663,0.0971407,0.103651,
    0.122125,0.097013,0.103585,0.122097,0.0972182,0.216509,
    0.240849,0.238243,0.212024,0.576429,1.42953,2.09393,1.57946
  };
  double v2_1020[nbins] = {
    2.17225,3.66955,5.32934,6.65091,7.67923,8.95035,9.92419,
    10.6285,11.4726,12.4663,13.2892,13.1702,13.383,12.3723,8.93792,
    12.5296,15.0085
  };
  double errV2_1020[nbins] = {
    0.0162515,0.0153472,0.0199456,0.0308611,
    0.031683,0.0419599,0.0651578,0.0672946,0.0898034,0.100566,
    0.16173,0.243491,0.434311,0.50903,1.02401,3.4382,9.82681
  };

 
//-----5-10% -----------------------------------------------------------------------------------

  double pt_510[nbins] = {
    0.298249,0.478648,0.703017,0.896405,1.078,1.30306,1.49641,
    1.6779,1.90274,2.18332,2.55882,2.96126,3.38753,3.9217,5.065,7.43055,10.4663
  };
  double errPtLow_510[nbins] = {
    0.0982494,0.0786479,0.103017,0.0964049,
    0.0779983,0.10306,0.0964067,0.0778999,0.102742,0.183315,
    0.158818,0.161262,0.187531,0.321702,0.565,0.930555,0.966266
  };
  double errPtHigh_510[nbins] = {
    0.101751,0.121352,0.0969831,0.103595,
    0.122002,0.0969397,0.103593,0.1221,0.0972583,0.216685,
    0.241182,0.238738,0.212469,0.578298,1.435,2.06945,1.53373
  };
  double v2_510[nbins] = {
    1.336,2.31457,3.34265,4.2049,4.76444,5.67546,5.92831,
    6.72256,6.95237,7.25684,8.04808,8.10241,5.13828,9.00966,
    1.0595,-7.48113,34.1494
  };
  double errV2_510[nbins] = {
    0.0490133,0.0459249,0.0593962,0.0917575,
    0.0941489,0.124114,0.193714,0.198987,0.266057,0.298391,
    0.482814,0.729011,1.32001,1.56362,3.12691,11.0502,56.8878
  };



  double errPt[nbins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  double v2_7080_ecc[nbins], errV2_7080_ecc[nbins];
  double v2_6070_ecc[nbins], errV2_6070_ecc[nbins];
  double v2_5060_ecc[nbins], errV2_5060_ecc[nbins];
  double v2_4050_ecc[nbins], errV2_4050_ecc[nbins];
  double v2_3040_ecc[nbins], errV2_3040_ecc[nbins];
  double v2_2030_ecc[nbins], errV2_2030_ecc[nbins];
  double v2_1020_ecc[nbins], errV2_1020_ecc[nbins];
  double v2_510_ecc[nbins], errV2_510_ecc[nbins];
  
 
 for (int i= 0; i< nbins; i++) {
  v2_7080_ecc[i] = 0.01*v2_7080[i];
  errV2_7080_ecc[i] = 0.01*errV2_7080[i];

  v2_6070_ecc[i] = 0.01*v2_6070[i];
  errV2_6070_ecc[i] = 0.01*errV2_6070[i];

  v2_5060_ecc[i] = 0.01*v2_5060[i];
  errV2_5060_ecc[i] = 0.01*errV2_5060[i];

  v2_4050_ecc[i] = 0.01*v2_4050[i];
  errV2_4050_ecc[i] = 0.01*errV2_4050[i];

  v2_3040_ecc[i] = 0.01*v2_3040[i];
  errV2_3040_ecc[i] = 0.01*errV2_3040[i];

  v2_2030_ecc[i] = 0.01*v2_2030[i];
  errV2_2030_ecc[i] = 0.01*errV2_2030[i];

  v2_1020_ecc[i] = 0.01*v2_1020[i];
  errV2_1020_ecc[i] = 0.01*errV2_1020[i];

  v2_510_ecc[i] = 0.01*v2_510[i];
  errV2_510_ecc[i] = 0.01*errV2_510[i];
 }
   
 //===================================================================================================================
 // centrality 10-20% 
 Double_t xCumulant4th1020ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
5.500000,7.000000,9.000000};
 Double_t yCumulant4th1020ALICE[] = {0.000000,0.000000,0.024075,0.031505,0.040413,0.044981,0.055358,0.060563,0.063378,0.070030,
0.082692,0.091611,0.099641,0.107223,0.122376,0.131240,0.137425,0.146050,0.131365,0.124708,
0.000000,0.000000,0.000000};
 Double_t xErrCumulant4th1020ALICE[23] = {0.};
 Double_t yErrCumulant4th1020ALICE[] = {0.000000,0.000000,0.002413,0.002931,0.003444,0.003950,0.004338,0.004835,0.005059,0.005586,
0.004521,0.005278,0.005999,0.007072,0.008260,0.007279,0.011897,0.017409,0.023995,0.025701,
0.000000,0.000000,0.000000};
 Int_t nPointsCumulant4th1020ALICE = sizeof(xCumulant4th1020ALICE)/sizeof(Double_t);                                      
 TGraphErrors *Cumulant4th1020ALICE = new TGraphErrors(nPointsCumulant4th1020ALICE,xCumulant4th1020ALICE,yCumulant4th1020ALICE,
                                                       xErrCumulant4th1020ALICE,yErrCumulant4th1020ALICE);
 Cumulant4th1020ALICE->SetMarkerStyle(kFullCircle);
 Cumulant4th1020ALICE->SetMarkerColor(kBlue);
 Cumulant4th1020ALICE->SetMarkerSize(1.2);
 Cumulant4th1020ALICE->SetFillStyle(1001);
 Cumulant4th1020ALICE->SetFillColor(kBlue-10);
 //===================================================================================================================
   
 //===================================================================================================================
 // centrality 20-30% 
 Double_t xCumulant4th2030ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
5.500000,7.000000,9.000000};
 Double_t yCumulant4th2030ALICE[] = {0.000000,0.000000,0.030926,0.041076,0.052063,0.059429,0.070500,0.084461,0.086745,0.099254,
0.109691,0.116398,0.130831,0.141959,0.158932,0.169680,0.171387,0.178858,0.171475,0.140358,
0.000000,0.000000,0.000000};
 Double_t xErrCumulant4th2030ALICE[23] = {0.};
 Double_t yErrCumulant4th2030ALICE[] = {0.000000,0.000000,0.002857,0.003451,0.003567,0.003859,0.004609,0.004976,0.005412,0.006277,
0.004748,0.005808,0.006896,0.007987,0.008683,0.008080,0.013278,0.018413,0.024873,0.026057,
0.000000,0.000000,0.000000};
 Int_t nPointsCumulant4th2030ALICE = sizeof(xCumulant4th2030ALICE)/sizeof(Double_t);                                      
 TGraphErrors *Cumulant4th2030ALICE = new TGraphErrors(nPointsCumulant4th2030ALICE,xCumulant4th2030ALICE,yCumulant4th2030ALICE,
                                                       xErrCumulant4th2030ALICE,yErrCumulant4th2030ALICE);
 Cumulant4th2030ALICE->SetMarkerStyle(kFullSquare);
 Cumulant4th2030ALICE->SetMarkerColor(kRed);
 Cumulant4th2030ALICE->SetMarkerSize(1.2);
 Cumulant4th2030ALICE->SetFillStyle(1001);
 Cumulant4th2030ALICE->SetFillColor(kRed-10);
 //ShiftAlongXaxis(Cumulant4th2030ALICE,0.03); 
 //===================================================================================================================
 
 //===================================================================================================================
 // centrality 30-40% 
 Double_t xCumulant4th3040ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
5.500000,7.000000,9.000000};
 Double_t yCumulant4th3040ALICE[] = {0.000000,0.000000,0.037071,0.048566,0.061083,0.070910,0.078831,0.091396,0.102026,0.109691,
0.124449,0.139819,0.155561,0.165701,0.173678,0.191149,0.202015,0.204540,0.212560,0.195885,
0.000000,0.000000,0.000000};
 Double_t xErrCumulant4th3040ALICE[23] = {0.};
 Double_t yErrCumulant4th3040ALICE[] = {0.000000,0.000000,0.002992,0.003364,0.003669,0.003931,0.004698,0.005261,0.005446,0.006151,
0.004980,0.005741,0.007198,0.008576,0.010868,0.009926,0.015269,0.020691,0.027601,0.031834,
0.000000,0.000000,0.000000};
 Int_t nPointsCumulant4th3040ALICE = sizeof(xCumulant4th3040ALICE)/sizeof(Double_t);                                      
 TGraphErrors *Cumulant4th3040ALICE = new TGraphErrors(nPointsCumulant4th3040ALICE,xCumulant4th3040ALICE,yCumulant4th3040ALICE,
                                                       xErrCumulant4th3040ALICE,yErrCumulant4th3040ALICE);
 Cumulant4th3040ALICE->SetMarkerStyle(kFullTriangleUp);
 Cumulant4th3040ALICE->SetMarkerColor(kGreen+2);
 Cumulant4th3040ALICE->SetMarkerSize(1.2);
 Cumulant4th3040ALICE->SetFillStyle(1001);
 Cumulant4th3040ALICE->SetFillColor(kGreen+2);
 //ShiftAlongXaxis(Cumulant4th3040ALICE,0.0); 
 //===================================================================================================================
     
 //===================================================================================================================
 // centrality 40-50% 
 // v2{2}
 Double_t xCumulant2nd4050ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
5.500000,7.000000,9.000000};
 Double_t yCumulant2nd4050ALICE[] = {0.000000,0.000000,0.043400,0.059911,0.073516,0.089756,0.105486,0.117391,0.128199,0.138013,
0.158271,0.177726,0.196383,0.208277,0.216648,0.242954,0.249961,0.240131,0.269006,0.207796,
0.000000,0.000000,0.000000};
 Double_t xErrCumulant2nd4050ALICE[23] = {0.};
 Double_t yErrCumulant2nd4050ALICE[] = {0.000000,0.000000,0.001496,0.001611,0.001700,0.001895,0.002081,0.002335,0.002598,0.002839,
0.002385,0.002908,0.003639,0.004457,0.005427,0.004813,0.007465,0.011109,0.016086,0.018067,
0.000000,0.000000,0.000000};
 Int_t nPointsCumulant2nd4050ALICE = sizeof(xCumulant2nd4050ALICE)/sizeof(Double_t);                                      
 TGraphErrors *Cumulant2nd4050ALICE = new TGraphErrors(nPointsCumulant2nd4050ALICE,xCumulant2nd4050ALICE,yCumulant2nd4050ALICE,
                                                       xErrCumulant2nd4050ALICE,yErrCumulant2nd4050ALICE);
 Cumulant2nd4050ALICE->SetMarkerStyle(kStar);
 Cumulant2nd4050ALICE->SetMarkerColor(kBlue);
 Cumulant2nd4050ALICE->SetMarkerSize(1.2); 
 Cumulant2nd4050ALICE->SetFillStyle(1001);
 Cumulant2nd4050ALICE->SetFillColor(kBlue-10);
 //ShiftAlongXaxis(Cumulant2nd4050ALICE,0.044); 

 // v2{4}
 Double_t xCumulant4th4050ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
5.500000,7.000000,9.000000};
 Double_t yCumulant4th4050ALICE[] = {0.000000,0.000000,0.038646,0.049824,0.066662,0.075856,0.081583,0.099778,0.104674,0.118545,
0.131874,0.152959,0.155348,0.169751,0.179052,0.178532,0.198851,0.185737,0.239901,0.186098,
0.000000,0.000000,0.000000};
 Double_t xErrCumulant4th4050ALICE[23] = {0.};
 Double_t yErrCumulant4th4050ALICE[] = {0.000000,0.000000,0.004853,0.006008,0.006304,0.006887,0.007797,0.008258,0.009678,0.010500,
0.008706,0.010322,0.012939,0.015934,0.019074,0.016231,0.023799,0.033236,0.046821,0.044138,
0.000000,0.000000,0.000000};
 Int_t nPointsCumulant4th4050ALICE = sizeof(xCumulant4th4050ALICE)/sizeof(Double_t);                                      
 TGraphErrors *Cumulant4th4050ALICE = new TGraphErrors(nPointsCumulant4th4050ALICE,xCumulant4th4050ALICE,yCumulant4th4050ALICE,
                                                       xErrCumulant4th4050ALICE,yErrCumulant4th4050ALICE);
 Cumulant4th4050ALICE->SetMarkerStyle(kFullTriangleDown);
 Cumulant4th4050ALICE->SetMarkerColor(kRed);
 Cumulant4th4050ALICE->SetMarkerSize(1.2); 
 Cumulant4th4050ALICE->SetFillStyle(1001);
 Cumulant4th4050ALICE->SetFillColor(kRed-10);
 //ShiftAlongXaxis(Cumulant4th4050ALICE,0.044); 
 //===================================================================================================================
 
 v2_1020_gr = new TGraphErrors(nbins,pt_1020, v2_1020_ecc, errPt, errV2_1020_ecc);
 v2_1020_gr->SetMarkerStyle(22);
 v2_1020_gr->SetMarkerSize(1.2);
 v2_1020_gr->SetLineColor(kBlack);
 v2_1020_gr->SetMarkerColor(4);
 v2_1020_gr->SetLineWidth(2);
 v2_1020_gr->SetFillStyle(3244);
 v2_1020_gr->SetFillColor(kBlack);
	
 v2_2030_gr = new TGraphErrors(nbins,pt_2030, v2_2030_ecc, errPt, errV2_2030_ecc);
 v2_2030_gr->SetMarkerStyle(22);
 v2_2030_gr->SetMarkerSize(1.2);
 v2_2030_gr->SetLineColor(kBlack);
 v2_2030_gr->SetMarkerColor(4);
 v2_2030_gr->SetLineWidth(2);
 v2_2030_gr->SetFillStyle(3144);
 v2_2030_gr->SetFillColor(kBlack);
	
 v2_3040_gr = new TGraphErrors(nbins,pt_3040, v2_3040_ecc, errPt, errV2_3040_ecc);
 v2_3040_gr->SetMarkerStyle(22);
 v2_3040_gr->SetMarkerSize(1.2);
 v2_3040_gr->SetLineColor(kBlack);
 v2_3040_gr->SetMarkerColor(4);
 v2_3040_gr->SetLineWidth(2);
 v2_3040_gr->SetFillStyle(1001);
 v2_3040_gr->SetFillColor(kGray+1);

 v2_4050_gr = new TGraphErrors(nbins,pt_4050, v2_4050_ecc, errPt, errV2_4050_ecc);
 v2_4050_gr->SetMarkerStyle(22);
 v2_4050_gr->SetMarkerSize(1.2);
 v2_4050_gr->SetLineColor(kBlack);
 v2_4050_gr->SetMarkerColor(4);
 v2_4050_gr->SetLineWidth(2);
 v2_4050_gr->SetFillStyle(1001);
 v2_4050_gr->SetFillColor(kGray+1);
	
 // Legend in uppper panel:
 TLegend *legendUp = new TLegend(0.20,0.57,0.37,0.94);
 legendUp->SetFillStyle(0); // white legendUp background
 legendUp->AddEntry(Cumulant2nd4050ALICE,"v_{2}{2}","p"); 
 legendUp->AddEntry(Cumulant4th4050ALICE,"v_{2}{4}","p"); 
 legendUp->AddEntry(v2_4050_gr,"v_{2}{4} (STAR)","f"); 
 
 legendUp->SetTextSize(0.050);
 legendUp->SetTextFont(22); // 22 = Times New Roman (bold)
 
 // Text in uppper panel:
 TPaveText *paveTextUp = new TPaveText(2.49,0.04,2.74,0.07,"br");
 paveTextUp->AddText("centrality 40-50%");
 paveTextUp->SetTextFont(22); // 22 = Times New Roman (bold)
 paveTextUp->SetTextSize(0.10);
 paveTextUp->SetFillColor(0);  // white pave text background
 
 // Final drawing u=in upper panel (order is important!):
 TCanvas *v2pt = new TCanvas("v2pt","v2pt",10,10,750,700);
 v2pt->Divide(1,2,0,0);
 v2pt->cd(1);
 StyleHistogramVsPtUp()->Draw(); 
 v2_4050_gr->Draw("same3"); // mesh
 v2_4050_gr->Draw("lsameX"); // line
 //gPad->SetTopMargin(0);  
 gPad->SetBottomMargin(0.0);  
 Cumulant2nd4050ALICE->Draw("psame");
 Cumulant4th4050ALICE->Draw("psame");
 legendUp->Draw("same");
 paveTextUp->Draw("same");
 v2pt->cd(2);
 StyleHistogramVsPtUp()->Draw(); 
 v2_4050_gr->Draw("same3"); // mesh
 v2_4050_gr->Draw("lsameX"); // line
 gPad->SetTopMargin(0);  
 Cumulant2nd4050ALICE->Draw("psame");
 Cumulant4th4050ALICE->Draw("psame");

 // Legend in lower panel:
 v2pt->cd(2);
 StyleHistogramVsPtDown()->Draw(); 
 TLegend *legendDown = new TLegend(0.20,0.59,0.39,0.96);
 legendDown->SetFillStyle(0); // white legend background
 legendDown->AddEntry(Cumulant4th1020ALICE,"10-20%","p"); 
 legendDown->AddEntry(Cumulant4th2030ALICE,"20-30%","p"); 
 legendDown->AddEntry(Cumulant4th3040ALICE,"30-40%","p"); 
 legendDown->AddEntry(v2_1020_gr,"10-20% (STAR)","f"); 
 legendDown->AddEntry(v2_2030_gr,"20-30% (STAR)","f"); 
 legendDown->AddEntry(v2_3040_gr,"30-40% (STAR)","f"); 
 v2_3040_gr->Draw("same3"); // mash  
 v2_3040_gr->Draw("lsameX"); // line  
 v2_2030_gr->Draw("same3"); // mash
 v2_2030_gr->Draw("lsameX"); // line 
 v2_1020_gr->Draw("same3"); // mash
 v2_1020_gr->Draw("lsameX"); // line 
 Cumulant4th1020ALICE->Draw("psame");
 Cumulant4th2030ALICE->Draw("psame");
 Cumulant4th3040ALICE->Draw("psame");
 legendDown->Draw("same");
 
 legendDown->SetTextSize(0.038);
 legendDown->SetTextFont(22); // 22 = Times New Roman (bold) 
 
 //v2pt->Print("fig2.pdf");    
 
} // end of void figure2b()

//==================================================================================

TH1D* StyleHistogramVsPtUp()
{
 // Style histogram:
 
 Int_t nPtBins = 100; // to be improved - hardwired 100
 Double_t ptMin = 0.; // in [GeV/c] // to be improved - hardwired 0
 Double_t ptMax = 5.044; // in [GeV/c] // to be improved - hardwired 5.75
 TString xTitle  = "p_{t} (GeV/#font[72]{c})";
 TString yTitle  = "v_{2}";
   
 TH1D* hist = new TH1D("","",nPtBins,ptMin,ptMax);
 hist->SetXTitle(xTitle.Data());
 hist->SetYTitle(yTitle.Data());
 hist->GetYaxis()->SetLabelSize(0.06);
 hist->GetYaxis()->SetTitleSize(0.07);
 hist->GetYaxis()->SetTitleOffset(1.10);
 hist->SetMinimum(0.00001); // minimum y-value - to be improved
 hist->SetMaximum(0.34499); // maximum y-value - to be improved
 //hist->SetLineStyle(7);
 return hist;
 
} // end of TH1D* StyleHistogramVsPtUp()
  
// =====================================================================================

TH1D* StyleHistogramVsPtDown()
{
 // Style histogram:
 
 Int_t nPtBins = 100; // to be improved - hardwired 100
 Double_t ptMin = 0.; // in [GeV/c] // to be improved - hardwired 0
 Double_t ptMax = 5.044; // in [GeV/c] // to be improved - hardwired 5.75
 TString xTitle  = "p_{t} (GeV/#font[72]{c})";
 TString yTitle  = "v_{2}{4}";
   
 TH1D* hist = new TH1D("","",nPtBins,ptMin,ptMax);
 hist->SetXTitle(xTitle.Data());
 hist->SetYTitle(yTitle.Data());
 hist->SetMinimum(0.00001); // minimum y-value - to be improved
 hist->SetMaximum(0.29999); // maximum y-value - to be improved
 //hist->SetLineStyle(7);
 return hist;
 
} // end of TH1D* StyleHistogramVsPtDown()
  
// =====================================================================================


void RemoveZeroPoints(TGraphErrors *ge)
{
 // Remove zero points from TGraphErrors.

 if(!ge){return;}

 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 {
  ge->GetPoint(p,x,y);
  if(TMath::Abs(y)<1.e-10)
  {
   ge->RemovePoint(p);
   cout<<Form(" WARNING (%s): point %d is < 1.e-10 and it was removed from the plot !!!!",ge->GetName(),p+1)<<endl;
  }
 } // end of for(Int_t p=0;p<nPoints;p++)
 
 cout<<endl;
 return;
 
} // end of void RemoveZeroPoints(TGraphErrors *ge)

// =====================================================================================

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
 // Shift original TGraphErrors along x-axis by amount determined by 'shift'.
 
 if(!ge)
 {
  printf("\n WARNING: ge is NULL in ShiftAlongXaxis() !!!! \n\n");
  return;
 }
  
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 { 
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 } // end of for(Int_t p=0;p<nPoints;p++)

} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

// =====================================================================================

void SetFlowStyle()
{
 // Set style which will affect all plots.
 
 gStyle->Reset();
 // gStyle->SetOptitle(0);
 // gStyle->SetOptStat(0);
 //gStyle->SetOptDate(1);
 // gStyle->SetPalette(8,0);  // (1,0)
 gStyle->SetPalette(1);  // (1,0)
 gStyle->SetDrawBorder(0);
 // gStyle->SetFillColor(0);  // kills palete ???
 gStyle->SetCanvasColor(0);
 gStyle->SetPadColor(0);
 // gStyle->SetFillColor(0); // otherwize it affects Fill colors later
 gStyle->SetFrameFillColor(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetFrameLineWidth(2);
 // gStyle->SetFrameFillStyle(4000);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetPadBottomMargin(0.15);
 gStyle->SetPadLeftMargin(0.15);
 gStyle->SetHistLineWidth(2);
 gStyle->SetFuncWidth(2);
 gStyle->SetLineWidth(2);
 gStyle->SetLabelSize(0.05,"xyz");
 gStyle->SetLabelOffset(0.01,"y");
 gStyle->SetLabelColor(kBlack,"xyz");
 gStyle->SetTitleSize(0.06,"xyz");
 gStyle->SetTitleOffset(1.3,"y");
 gStyle->SetTitleFillColor(0);
 gStyle->SetLineWidth(2);  
 gStyle->SetHistLineColor(1);
 gStyle->SetTextColor(1);
 gStyle->SetTitleTextColor(1);
 TGaxis::SetMaxDigits(4);
 gStyle->SetOptStat(0); // removes stat. box from all histos
 gROOT->ForceStyle();

} // end of void SetFlowStyle()

