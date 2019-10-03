//=======================================================================//
//Macro to draw the main results of the MH analysis:
//i) drawCentralityDependence3P:: Draws the <cos(psi1 + psi2 - 2psi3) 
//   vs centrality percentile
//ii) drawNpartDependence3P:: Draws the <cos(psi1 + psi2 - 2psi3) 
//    vs Nparticipants
//iii) drawNpartDependenceScaled3P:: Draws the Npart*<cos(psi1 + psi2 - 2psi3) 
//     vs Nparticipants
//iv) drawCentralityDependenceRP:: Draws the <cos(psi1 + psi2 - 2PsiRP) 
//   vs centrality percentile
//v) drawNpartDependenceRP:: Draws the <cos(psi1 + psi2 - 2PsiRP) 
//    vs Nparticipants
//vi) drawNpartDependenceScaledRP:: Draws the Npart*<cos(psi1 + psi2 - 2PsiRP) 
//    vs Nparticipants
//=======================================================================//
Bool_t gPreliminary = kFALSE;

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 9;
TString strCentralityBins[nCentralityBins] = {"0-5","5-10","10-20",
					      "20-30","30-40","40-50",
					      "50-60","60-70","70-80"};
//Double_t gCentralityPercentile[nCentralityBins] = {75.,65.,55.,45.,35.,25.,15.,7.5,2.5};
Double_t gCentralityPercentileStar[nCentralityBins] = {1.75,6.75,14.,24.,34.,44.,54.,64.,74.};
Double_t gCentralityPercentile[nCentralityBins] = {3.0,8.0,16.,26.,36.,46.,56.,66.,76.};
Double_t gCentralityPercentile2[nCentralityBins] = {0.5,5.5,12.,22.,32.,42.,52.,62.,72.};
Double_t gCentralityPercentile3[nCentralityBins] = {4.25,9.25,18.,28.,38.,48.,58.,68.,78.};
Double_t gCentralityPercentileError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
// Double_t gCentralityPercentileSystematicError[nCentralityBins] = {1.5,1.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
Double_t gCentralityPercentileSystematicError[nCentralityBins] = {0,0,0,0,0,0,0,0,0};
Double_t gCentralityPercentileStarSystematicError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};

//================================HIJING================================//
Double_t g3pCorrelatorHijingNoFlowPlusMinus[nCentralityBins] = {-1.20093e-07, 1.84157e-07, 1.76266e-07, 7.57471e-07, 1.19105e-06, 4.13612e-06, 1.40261e-05, 4.40286e-05, 9.745345e-05};
Double_t g3pCorrelatorHijingNoFlowPlusMinusError[nCentralityBins] = {9.3131e-08, 8.68461e-08, 1.21384e-07, 2.13452e-07, 1.12347e-07, 3.13793e-07, 8.37252e-07, 3.15813e-06, 9.655345e-06};
Double_t g3pCorrelatorHijingNoFlowSameCharge[nCentralityBins] = {1.77386e-08, 9.26015e-10, 7.43216e-08, 1.76705e-07, 8.10251e-07, 2.68011e-06, 9.20014e-06, 3.3102e-05, 9.95345e-05};
Double_t g3pCorrelatorHijingNoFlowSameChargeError[nCentralityBins] = {2.97636e-08,4.08723e-08, 4.65237e-08, 1.12774e-07, 1.99168e-07, 5.83821e-07, 1.56612e-06, 5.90915e-06, 2.12124e-05};
/*Double_t g3pCorrelatorHijingNoFlowPlusMinus[nCentralityBins] = {-1.20093e-07,2.1169e-07,8.38551e-07,1.71174e-06,4.37847e-06,1.05875e-05,2.98748e-05,0.000100555,0.000289372};
Double_t g3pCorrelatorHijingNoFlowPlusMinusError[nCentralityBins] = {9.3131e-08,5.3237e-08,4.87755e-08,9.58588e-08,2.42031e-07,4.74851e-07,1.23657e-06,3.81097e-06,1.27354e-05};
Double_t g3pCorrelatorHijingNoFlowSameCharge[nCentralityBins] = {1.77386e-08,1.3287e-07,4.17543e-07,5.91545e-07,3.51865e-06,7.17039e-06,2.00607e-05,6.20528e-05,0.000176748};
Double_t g3pCorrelatorHijingNoFlowSameChargeError[nCentralityBins] = {2.97636e-08,4.3382e-08,8.9158e-08,1.76401e-07,4.46464e-07,8.82932e-07,2.32232e-06,7.11183e-06,2.31168e-05};*/

/*Double_t g3pCorrelatorHijingWithFlowPlusMinus[nCentralityBins] = {-1.20093e-07, 1.56519e-07, 2.90847e-07, 7.57471e-07, 1.17141e-06, 4.13612e-06, 2.4549e-06, 2.10261e-05, 7.40286e-05};
Double_t g3pCorrelatorHijingWithFlowPlusMinusError[nCentralityBins] = {9.3131e-08, 7.9856e-08, 7.42477e-08, 2.13452e-07, 1.07302e-07, 3.13793e-07, 8.2704e-07, 9.37252e-07, 5.15813e-06};
Double_t g3pCorrelatorHijingWithFlowSameCharge[nCentralityBins] = {-5.78124e-07, 4.16894e-08, 3.8112e-07, -1.70606e-07, 8.10251e-07, 2.68011e-06, 2.4549e-06, 3.1102e-05, 8.95345e-05};
Double_t g3pCorrelatorHijingWithFlowSameChargeError[nCentralityBins] = {1.80581e-07, 1.45975e-07, 1.56996e-07, 3.96189e-07, 1.99168e-07, 5.83821e-07, 8.2704e-07, 7.90915e-06, 4.12124e-05};*/
/*Double_t g3pCorrelatorHijingWithFlowSameCharge[nCentralityBins] = {1.3299e-07,8.4476e-07,2.96354e-06,6.80478e-06,1.45245e-05,2.76715e-05,5.86195e-05,0.00012099,0.000279274};
Double_t g3pCorrelatorHijingWithFlowSameChargeError[nCentralityBins] = {9.7787e-08,1.6551e-07,2.17682e-07,4.47015e-07,1.05436e-06,8.82932e-07,5.47122e-06,9.83962e-06,2.75911e-05};
Double_t g3pCorrelatorHijingWithFlowPlusMinus[nCentralityBins] = {2.7632e-07,7.3227e-07,4.05092e-06,9.67432e-06,1.98196e-05,3.64523e-05,7.97305e-05,0.000172174,0.000412814};
Double_t g3pCorrelatorHijingWithFlowPlusMinusError[nCentralityBins] = {8.3287e-08,1.8366e-07,1.18877e-07,2.43468e-07,5.71375e-07,9.59e-07,2.05589e-06,5.28575e-06,1.51409e-05};*/
Double_t g3pCorrelatorHijingWithFlowSameCharge[nCentralityBins] = {1.3299e-08,4.4476e-07,1.96354e-06,3.80478e-06,7.45245e-06,0.76715e-05,2.86195e-05,0.00012099,0.000279274};
Double_t g3pCorrelatorHijingWithFlowSameChargeError[nCentralityBins] = {9.7787e-08,1.6551e-07,2.17682e-07,4.47015e-07,6.05436e-08,8.82932e-07,5.47122e-06,9.83962e-06,2.75911e-05};
Double_t g3pCorrelatorHijingWithFlowPlusMinus[nCentralityBins] = {2.7632e-08,5.3227e-07,2.05092e-06,6.67432e-06,9.98196e-06,2.64523e-05,4.97305e-05,0.000172174,0.000412814};
Double_t g3pCorrelatorHijingWithFlowPlusMinusError[nCentralityBins] = {8.3287e-08,1.8366e-07,1.18877e-07,2.43468e-07,5.71375e-07,9.59e-07,2.05589e-06,5.28575e-06,1.51409e-05};

Double_t g3pHijingNoFlowAvg[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingNoFlowAvgError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingNoFlowPlusMinus[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingNoFlowPlusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingNoFlowSameCharge[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingNoFlowSameChargeError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

Double_t g3pHijingWithFlowAvg[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingWithFlowAvgError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingWithFlowPlusMinus[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingWithFlowPlusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingWithFlowSameCharge[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t g3pHijingWithFlowSameChargeError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
//================================HIJING================================//

//================================STAR================================//
Double_t g3pCorrelatorSTARPlusPlus[nCentralityBins] = {-2.72475e-05,
						       -4.84767e-05,
						       -8.43622e-05,
						       -0.000139391,
						       -0.000212998,
						       -0.000310464,
						       -0.000449019,
						       -0.000531625,
						       10.};
Double_t g3pCorrelatorSTARPlusPlusError[nCentralityBins] = {2.51676e-06,
							    2.38477e-06,
							    1.96375e-06,
							    2.50148e-06,
							    3.87569e-06,
							    6.65112e-06,
							    1.36251e-05,
							    3.49563e-05,
							    0.0};
Double_t g3pCorrelatorSTARPlusPlusSysErrorHigh[nCentralityBins] = {-2.372263e-05,-4.559134e-05,-8.079906e-05,-0.0001334556,-0.0002038168,-0.0002934094,-0.0004116487,-0.0004651086};
Double_t g3pCorrelatorSTARPlusPlusSysErrorLow[nCentralityBins] = {-7.013656e-05,-7.013656e-05,-0.0001010264,-0.0001578065,-0.0002416504,-0.0003675497,-0.0005460016,-0.0007163518};

Double_t g3pCorrelatorSTARPlusMinus[nCentralityBins] = {-1.16758e-05,
							-8.2939e-06,
							-7.9278e-06,
							-5.81744e-06,
							-3.32537e-06,
							1.71515e-05,
							5.88244e-05,
							0.000205355,
							10.};
Double_t g3pCorrelatorSTARPlusMinusError[nCentralityBins] = {2.55046e-06,
							     2.41575e-06,
							     1.99013e-06,
							     2.52889e-06,
							     3.90957e-06,
							     6.65273e-06,
							     1.3518e-05,
							     3.40235e-05,
							     0.0};
Double_t g3pCorrelatorSTARPlusMinusSysErrorHigh[nCentralityBins] = {-1.202044e-05,-8.387146e-06,-7.959642e-06,-5.834463e-06,-3.335567e-06,1.72155e-05,5.916783e-05,0.0002078669};
Double_t g3pCorrelatorSTARPlusMinusSysErrorLow[nCentralityBins] = {-1.199969e-05,-1.199969e-05,-9.493789e-06,-6.585985e-06,-3.772705e-06,2.030513e-05,7.152968e-05,0.0002767107};
//================================STAR================================//

//================================ALICE================================//
Double_t gv2[nCentralityBins] = {0.0218688,0.0410256,0.0617027,0.0803975,0.0914567,0.0950293,0.0930474,0.0871609,0.0854715};
Double_t gv2Error[nCentralityBins] = {0.00777206,0.0055282,0.00461509,0.00540859,0.00673375,0.00916149,0.0118708,0.0150726,0.0159525};

Double_t g3pCorrelatorALICEDataSameCharge[nCentralityBins] = {-1.60885e-05,
							      -2.33387e-05,
							      -5.49677e-05,
							      -0.000105864,
							      -0.000166243,
							      -0.000276424,
							      -0.00037915,
							      -0.000425016,
							      -0.000293409}; 
Double_t g3pCorrelatorALICEDataSameChargeError[nCentralityBins] = {3.22359e-06,
								   3.45502e-06,
								   3.26195e-06,
								   4.85736e-06,
								   7.60362e-06,
								   1.31918e-05,
								   2.61394e-05,
								   6.2836e-05,
								   0.000196171};
//Double_t g3pCorrelatorALICEDataSameChargeSystematicError[nCentralityBins] = {3.58357E-05,1.66674E-05,1.27416E-05,1.98483E-05,3.97548E-05,8.97482E-05,0.00010582,0.000233554,0.00031035126};
Double_t g3pCorrelatorALICEDataSameChargeSystematicError[nCentralityBins] = {1.6383E-05,1.66674E-05,1.27416E-05,1.98483E-05,3.97548E-05,8.97482E-05,0.00010582,0.000233554,0.00031035126};

Double_t g3pCorrelatorALICEDataPlusMinus[nCentralityBins] = {-3.25362e-07,
							     3.86547e-06,
							     -3.41652e-06,
							     -1.40672e-05,
							     -2.64535e-05,
							     -3.86616e-05,
							     2.17488e-05,
							     0.000248817,
							     0.00105866};
Double_t g3pCorrelatorALICEDataPlusMinusError[nCentralityBins] = {4.40278e-06,
								  4.87877e-06,
								  4.57899e-06,
								  6.81769e-06,
								  1.06697e-05,
								  1.83207e-05,
								  3.62493e-05,
								  8.68075e-05,
								  0.000279468};
//Double_t g3pCorrelatorALICEDataPlusMinusSystematicError[nCentralityBins] = {3.55452E-05,1.8482E-05,6.31586E-06,1.20367E-05,3.66273E-05,4.12716E-05,8.08092E-05,0.000158143,0.00031038405};
Double_t g3pCorrelatorALICEDataPlusMinusSystematicError[nCentralityBins] = {6.80205E-06,1.8482E-05,6.31586E-06,1.20367E-05,3.66273E-05,4.12716E-05,8.08092E-05,0.000158143,0.00031038405};

Double_t g3pCorrelatorALICEDatav1Fluctustions[nCentralityBins] = {-4.26793e-06,-5.19642e-06,-1.57247e-05,-3.27355e-05,-5.33411e-05,-8.64989e-05,-9.82421e-05,-4.53941e-05,0.000229361};
Double_t g3pCorrelatorALICEDatav1FluctustionsError[nCentralityBins] = {3.939e-06,4.5402e-06,1.34675e-05,2.72301e-05,4.30072e-05,7.10439e-05,8.04585e-05,4.27054e-05,0.000153265};

Double_t g3pCorrelatorALICEDatav1Fluctustions[nCentralityBins] = {-4.26793e-06,-5.19642e-06,-1.57247e-05,-3.27355e-05,-5.33411e-05,-8.64989e-05,-9.82421e-05,-4.53941e-05,0.000229361};
Double_t g3pCorrelatorALICEDatav1FluctustionsError[nCentralityBins] = {3.939e-06,4.5402e-06,1.34675e-05,2.72301e-05,4.30072e-05,7.10439e-05,8.04585e-05,4.27054e-05,0.000153265};
//================================ALICE================================//

//============================Toneev et al.=============================//
Double_t g3pCorrelatorToneevSameCharge[4] = {-2.*8.4e-07,-2.*0.0168e-03,-2.*0.0338e-03,-2.*0.0588e-03};
Double_t g3pCorrelatorToneevSameChargeError[4] = {0.,0.,0.,0.};
Double_t gCentralityPercentileToneev[4] = {2.5,20.,45.,67.};
Double_t gCentralityPercentileErrorToneev[4] = {0.5,0.5,0.5,0.5};
//============================Toneev et al.=============================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void drawPaperFigure2()
{
  //Draws the <cos(psi1 + psi2 - 2psi3) vs centrality percentile
  gROOT->LoadMacro("SetFlowStyle.C");
  SetFlowStyle();

  TGaxis::SetMaxDigits(2);
  
  //Fix the asymmetric errors for STAR
  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    g3pCorrelatorSTARPlusPlusSysErrorHigh[iBin] = TMath::Abs(g3pCorrelatorSTARPlusPlus[iBin] - g3pCorrelatorSTARPlusPlusSysErrorHigh[iBin]);
    g3pCorrelatorSTARPlusPlusSysErrorLow[iBin] = TMath::Abs(g3pCorrelatorSTARPlusPlus[iBin] - g3pCorrelatorSTARPlusPlusSysErrorLow[iBin]);
    //Printf("%lf - %lf - %lf",g3pCorrelatorSTARPlusPlus[iBin],g3pCorrelatorSTARPlusPlusSysErrorLow[iBin],g3pCorrelatorSTARPlusPlusSysErrorHigh[iBin]);

    g3pCorrelatorSTARPlusMinusSysErrorHigh[iBin] = TMath::Abs(g3pCorrelatorSTARPlusMinus[iBin] - g3pCorrelatorSTARPlusMinusSysErrorHigh[iBin]);
    g3pCorrelatorSTARPlusMinusSysErrorLow[iBin] = TMath::Abs(g3pCorrelatorSTARPlusMinus[iBin] - g3pCorrelatorSTARPlusMinusSysErrorLow[iBin]);
    //Printf("%lf - %lf - %lf",g3pCorrelatorSTARPlusMinus[iBin],g3pCorrelatorSTARPlusMinusSysErrorLow[iBin],g3pCorrelatorSTARPlusMinusSysErrorHigh[iBin]);
  }

  //Teany-Yan
  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    g3pCorrelatorALICEDatav1Fluctustions[iBin] = 0.5*(g3pCorrelatorALICEDataPlusMinus[iBin] + g3pCorrelatorALICEDataSameCharge[iBin]);
    //g3pCorrelatorALICEDatav1FluctustionsError[iBin] = 0.1*g3pCorrelatorALICEDatav1Fluctustions[iBin];
    g3pCorrelatorALICEDatav1FluctustionsError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g3pCorrelatorALICEDataPlusMinusError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameChargeError[iBin],2));
  }

  //HIJING
  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    //no flow
    g3pHijingNoFlowPlusMinus[iBin] = g3pCorrelatorHijingNoFlowPlusMinus[iBin]/gv2[iBin];
    g3pHijingNoFlowPlusMinusError[iBin] = (1./gv2[iBin])*TMath::Sqrt(TMath::Power((g3pCorrelatorHijingNoFlowPlusMinus[iBin]*gv2Error[iBin]/gv2[iBin]),2) + TMath::Power(g3pCorrelatorHijingNoFlowPlusMinusError[iBin],2));

    g3pHijingNoFlowSameCharge[iBin] = g3pCorrelatorHijingNoFlowSameCharge[iBin]/gv2[iBin];
    g3pHijingNoFlowSameChargeError[iBin] = (1./gv2[iBin])*TMath::Sqrt(TMath::Power((g3pCorrelatorHijingNoFlowSameCharge[iBin]*gv2Error[iBin]/gv2[iBin]),2) + TMath::Power(g3pCorrelatorHijingNoFlowSameChargeError[iBin],2));

    g3pHijingNoFlowAvg[iBin] = 0.5*(g3pHijingNoFlowPlusMinus[iBin] + g3pHijingNoFlowSameCharge[iBin]);
    g3pHijingNoFlowAvgError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g3pHijingNoFlowPlusMinusError[iBin],2) + TMath::Power(g3pHijingNoFlowSameChargeError[iBin],2));

    //with flow
    g3pHijingWithFlowPlusMinus[iBin] = g3pCorrelatorHijingWithFlowPlusMinus[iBin]/gv2[iBin];
    g3pHijingWithFlowPlusMinusError[iBin] = (1./gv2[iBin])*TMath::Sqrt(TMath::Power((g3pCorrelatorHijingWithFlowPlusMinus[iBin]*gv2Error[iBin]/gv2[iBin]),2) + TMath::Power(g3pCorrelatorHijingWithFlowPlusMinusError[iBin],2));

    g3pHijingWithFlowSameCharge[iBin] = g3pCorrelatorHijingWithFlowSameCharge[iBin]/gv2[iBin];
    g3pHijingWithFlowSameChargeError[iBin] = (1./gv2[iBin])*TMath::Sqrt(TMath::Power((g3pCorrelatorHijingWithFlowSameCharge[iBin]*gv2Error[iBin]/gv2[iBin]),2) + TMath::Power(g3pCorrelatorHijingWithFlowSameChargeError[iBin],2));

    g3pHijingWithFlowAvg[iBin] = 0.5*(g3pHijingWithFlowPlusMinus[iBin] + g3pHijingWithFlowSameCharge[iBin]);
    g3pHijingWithFlowAvgError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g3pHijingWithFlowPlusMinusError[iBin],2) + TMath::Power(g3pHijingWithFlowSameChargeError[iBin],2));
  }

  //================================================//
  //v1 fluctuations
  TH1F *gHistv1Fluctuations = new TH1F("gHistv1Fluctuations",
				       ";Centrality percentile;#LT cos(#phi_{#alpha} + #phi_{#beta} - 2#Psi_{RP}) #GT",
				       nCentralityBins,0,93);
  myTGraphSetUpTH(gHistv1Fluctuations,1,myv1FluctColorSame,0,1,myv1FluctColorSame,10,1001,myv1FluctColorSame);

  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    gHistv1Fluctuations->SetBinContent(iBin,g3pCorrelatorALICEDatav1Fluctustions[iBin]);
    gHistv1Fluctuations->SetBinError(iBin,g3pCorrelatorALICEDatav1FluctustionsError[iBin]);
    
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pCorrelatorALICEDatav1Fluctustions[iBin]<<" ± "<<g3pCorrelatorALICEDatav1FluctustionsError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pCorrelatorALICEDataSameCharge[iBin]<<" ± "<<g3pCorrelatorALICEDataSameChargeError[iBin]<<" ± "<<g3pCorrelatorALICEDataSameChargeSystematicError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pCorrelatorALICEDataPlusMinus[iBin]<<" ± "<<g3pCorrelatorALICEDataPlusMinusError[iBin]<<" ± "<<g3pCorrelatorALICEDataPlusMinusSystematicError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingWithFlowPlusMinus[iBin]<<" ± "<<g3pHijingWithFlowPlusMinusError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingWithFlowSameCharge[iBin]<<" ± "<<g3pHijingWithFlowSameChargeError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingWithFlowPlusMinus[iBin]<<" ± "<<g3pHijingWithFlowPlusMinusError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingNoFlowSameCharge[iBin]<<" ± "<<g3pHijingNoFlowSameChargeError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingNoFlowPlusMinus[iBin]<<" ± "<<g3pHijingNoFlowPlusMinusError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingWithFlowPlusMinus[iBin]<<" ± "<<g3pHijingWithFlowPlusMinusError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingNoFlowAvg[iBin]<<" ± "<<g3pHijingNoFlowAvgError[iBin]<<endl;
    //cout<<"Centrality: "<<strCentralityBins[iBin].Data()<<" "<<g3pHijingWithFlowAvg[iBin]<<" ± "<<g3pHijingWithFlowAvgError[iBin]<<endl;
  }

  //================================================//
  //HIJING without flow
  TString drawOptions_HIJING = "PZ";
  int markerColor_HIJING = 32;
  int markerColor_HIJINGv2 = 32;
  int lineColor_HIJING = markerColor_HIJING;
  int lineColor_HIJINGv2 = markerColor_HIJINGv2;
  int lineStyle_HIJING = 1;
  int lineStyle_HIJINGv2 = 1;
  int markerStyle_HIJINGopp = 32;
  int markerStyle_HIJINGoppv2 = 26;
  int markerStyle_HIJINGsame = 23;
  int markerStyle_HIJINGsamev2 = 22;
  float markerSize_HIJING = 1.*myMarkerSize;
  TGraphErrors *grHIJINGNoFlowAvg = new TGraphErrors(nCentralityBins,
						     gCentralityPercentile2,
						     g3pHijingNoFlowAvg,
						     gCentralityPercentileError,
						     g3pHijingNoFlowAvgError);
  myTGraphSetUp(grHIJINGNoFlowAvg,markerStyle_HIJINGopp,markerColor_HIJING,markerSize_HIJING,lineStyle_HIJING,lineColor_HIJING,2,1001,markerColor_HIJING);

  TGraphErrors *grHIJINGNoFlowPlusMinus = new TGraphErrors(nCentralityBins,
							   gCentralityPercentile2,
							   g3pHijingNoFlowPlusMinus,
							   gCentralityPercentileError,
							   g3pHijingNoFlowPlusMinusError);
  myTGraphSetUp(grHIJINGNoFlowPlusMinus,markerStyle_HIJINGopp,markerColor_HIJING,markerSize_HIJING,lineStyle_HIJING,lineColor_HIJING,2,1001,markerColor_HIJING);

  TGraphErrors *grHIJINGNoFlowSameCharge = new TGraphErrors(nCentralityBins,
							    gCentralityPercentile2,
							    g3pHijingNoFlowSameCharge,
							    gCentralityPercentileError,
							    g3pHijingNoFlowSameChargeError);
  myTGraphSetUp(grHIJINGNoFlowSameCharge,markerStyle_HIJINGsame,markerColor_HIJING,markerSize_HIJING,lineStyle_HIJING,lineColor_HIJING,2,1001,markerColor_HIJING);
  //================================================//

  //================================================//
  //HIJING with flow
  TGraphErrors *grHIJINGWithFlowAvg = new TGraphErrors(nCentralityBins,
						       gCentralityPercentile3,
						       g3pHijingWithFlowAvg,
						       gCentralityPercentileError,
						       g3pHijingWithFlowAvgError);
  myTGraphSetUp(grHIJINGWithFlowAvg,markerStyle_HIJINGsamev2,markerColor_HIJINGv2,markerSize_HIJING,lineStyle_HIJINGv2,markerColor_HIJINGv2,2,1001,markerColor_HIJINGv2);

  TGraphErrors *grHIJINGWithFlowPlusMinus = new TGraphErrors(nCentralityBins,
							   gCentralityPercentile3,
							   g3pHijingWithFlowPlusMinus,
							   gCentralityPercentileError,
							   g3pHijingWithFlowPlusMinusError);
  myTGraphSetUp(grHIJINGWithFlowPlusMinus,markerStyle_HIJINGsamev2,markerColor_HIJINGv2,markerSize_HIJING,lineStyle_HIJINGv2,markerColor_HIJINGv2,2,1001,markerColor_HIJINGv2);

  TGraphErrors *grHIJINGWithFlowSameCharge = new TGraphErrors(nCentralityBins,
							    gCentralityPercentile3,
							    g3pHijingWithFlowSameCharge,
							    gCentralityPercentileError,
							    g3pHijingWithFlowSameChargeError);
  myTGraphSetUp(grHIJINGWithFlowSameCharge,markerStyle_HIJINGsamev2,markerColor_HIJINGv2,markerSize_HIJING,lineStyle_HIJINGv2,markerColor_HIJINGv2,2,1001,markerColor_HIJINGv2);
  //================================================//

  //================================================//
  //Toneev et al.
  TGraphErrors *grToneevSameCharge = new TGraphErrors(4,gCentralityPercentileToneev,g3pCorrelatorToneevSameCharge,gCentralityPercentileErrorToneev,g3pCorrelatorToneevSameChargeError);
  myTGraphSetUp(grToneevSameCharge,21,myToneevColor,3,1,myToneevColor,3,1001,myToneevColor);
  //================================================//
  int markerColor_ALICE = myDarkRed;
  int lineColor_ALICE = markerColor_ALICE;
  int markerStyle_ALICEsame = 20;
  int markerStyle_ALICEopp = 24;
  int color_sysErr_STAR = 17;
  int color_sysErr_ALICE = 25;
  //(+-)
  TGraphErrors *grALICEDataPlusMinus = new TGraphErrors(nCentralityBins,
							gCentralityPercentile,
							g3pCorrelatorALICEDataPlusMinus,
							gCentralityPercentileError,
							g3pCorrelatorALICEDataPlusMinusError);
  myTGraphSetUp(grALICEDataPlusMinus,markerStyle_ALICEopp,markerColor_ALICE,myMarkerSize,1,lineColor_ALICE,2,1001,markerColor_ALICE);

  TGraphErrors *grALICEDataPlusMinusSystematic = new TGraphErrors(nCentralityBins,
                                                                  gCentralityPercentile,
                                                                  g3pCorrelatorALICEDataPlusMinus,
                                                                  gCentralityPercentileSystematicError,
                                                                  g3pCorrelatorALICEDataPlusMinusSystematicError);
  myTGraphSetUp(grALICEDataPlusMinusSystematic,markerStyle_ALICEopp,color_sysErr_ALICE,myMarkerSize,1,color_sysErr_ALICE,10,1001,color_sysErr_ALICE);
  Printf("ALICE (+-)");
  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) 
    cout<<"Centrality: "<<strCentralityBins[iCentrality].Data()<<" - C "<<g3pCorrelatorALICEDataPlusMinus[iCentrality]<<" - Sys: "<<g3pCorrelatorALICEDataPlusMinusSystematicError[iCentrality]<<" ("<<100*TMath::Abs(g3pCorrelatorALICEDataPlusMinusSystematicError[iCentrality]/g3pCorrelatorALICEDataPlusMinus[iCentrality])<<"%)"<<endl;
    
  int markerColor_STAR = 9;
  int lineColor_STAR = markerColor_STAR;
  TGraphErrors *grSTARDataPlusMinus = new TGraphErrors(nCentralityBins,
						       gCentralityPercentileStar,
						       g3pCorrelatorSTARPlusMinus,
						       gCentralityPercentileError,
						       g3pCorrelatorSTARPlusMinusError);
  myTGraphSetUp(grSTARDataPlusMinus,30,markerColor_STAR,myMarkerSize+0.2,1,lineColor_STAR,2,1001,markerColor_STAR);
  TGraphAsymmErrors *grSTARDataPlusMinusSystematics = new TGraphAsymmErrors(nCentralityBins,gCentralityPercentileStar,g3pCorrelatorSTARPlusMinus,gCentralityPercentileError,gCentralityPercentileError,g3pCorrelatorSTARPlusMinusSysErrorHigh,g3pCorrelatorSTARPlusMinusSysErrorLow);
    myTGraphSetUp_Asym(grSTARDataPlusMinusSystematics,30,color_sysErr_STAR,myMarkerSize,1,color_sysErr_STAR,10,1001,color_sysErr_STAR);

  //================================================//
  //(--)&(++)
  TGraphErrors *grALICEDataSameCharge = new TGraphErrors(nCentralityBins,
							 gCentralityPercentile,
							 g3pCorrelatorALICEDataSameCharge,
							 gCentralityPercentileError,
							 g3pCorrelatorALICEDataSameChargeError);
  myTGraphSetUp(grALICEDataSameCharge,markerStyle_ALICEsame,markerColor_ALICE,myMarkerSize,1,lineColor_ALICE,2,1001,markerColor_ALICE);

  TGraphErrors *grALICEDataSameChargeSystematic = new TGraphErrors(nCentralityBins,
                                                                  gCentralityPercentile,
                                                                  g3pCorrelatorALICEDataSameCharge,
                                                                  gCentralityPercentileSystematicError,
                                                                  g3pCorrelatorALICEDataSameChargeSystematicError);;
  myTGraphSetUp(grALICEDataSameChargeSystematic,markerStyle_ALICEsame,color_sysErr_ALICE,myMarkerSize,1,color_sysErr_ALICE,10,1001,color_sysErr_ALICE);
    Printf("ALICE (++,--)");
  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) 
    cout<<"Centrality: "<<strCentralityBins[iCentrality].Data()<<" - C "<<g3pCorrelatorALICEDataSameCharge[iCentrality]<<" - Sys: "<<g3pCorrelatorALICEDataSameChargeSystematicError[iCentrality]<<" ("<<100*TMath::Abs(g3pCorrelatorALICEDataSameChargeSystematicError[iCentrality]/g3pCorrelatorALICEDataSameCharge[iCentrality])<<"%)"<<endl;

  TGraphErrors *grSTARDataPlusPlus = new TGraphErrors(nCentralityBins,
						      gCentralityPercentileStar,
						      g3pCorrelatorSTARPlusPlus,
						      gCentralityPercentileError,
						      g3pCorrelatorSTARPlusPlusError);
  myTGraphSetUp(grSTARDataPlusPlus,29,markerColor_STAR,myMarkerSize+0.4,1,lineColor_STAR,2,1001,markerColor_STAR);
  TGraphAsymmErrors *grSTARDataPlusPlusSystematics = new TGraphAsymmErrors(nCentralityBins,
									   gCentralityPercentileStar,
									   g3pCorrelatorSTARPlusPlus,
									   gCentralityPercentileError,
									   gCentralityPercentileError,
									   g3pCorrelatorSTARPlusPlusSysErrorLow,
									   g3pCorrelatorSTARPlusPlusSysErrorHigh);
    myTGraphSetUp_Asym(grSTARDataPlusPlusSystematics,29,color_sysErr_STAR,myMarkerSize+0.4,1,color_sysErr_STAR,10,1001,color_sysErr_STAR);
 
  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  //====================================//
  //<cos(psi1+psi2-2phi3)> vs centrality
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			  ";centrality, %;#LT cos(#phi_{#alpha} + #phi_{#beta} - 2#Psi_{RP}) #GT",
			  nCentralityBins,0,70,1000,-7.e-04,6.e-04);
  gEmpty1->SetStats(kFALSE);
  //gEmpty1->GetYaxis()->SetTitleOffset(1.8);
  //gEmpty1->GetXaxis()->SetTitleOffset(1.5);
  gEmpty1->GetYaxis()->SetTitleSize(0.07);
  gEmpty1->GetYaxis()->SetTitleOffset(0.85);
  gEmpty1->GetYaxis()->SetNdivisions(10);
  gEmpty1->GetXaxis()->SetNdivisions(10);

  TF1 *f1 = new TF1("f1","0",0,1000);
  f1->SetLineColor(1); f1->SetLineStyle(1); f1->SetLineWidth(1);

  TCanvas *c1 = new TCanvas("c1","Centrality dependence: 3p correlator",
			    0,0,800,600);
  TPad *myPad = new TPad("myPad", "The pad",0,0,1,1);
  myPadSetUp(myPad,0.15,0.065,0.04,0.15);
  myPad->Draw();
  myPad->cd();
  //c1->SetGridx(); c1->SetGridy();
  gEmpty1->Draw();
  f1->Draw("same");

  // shift various data points wrt. each other
  float shift_ALICE = -2.;
  float shift_STAR = 2.;
  float shift_HIJING = 6.5;
  float shift_HIJINGv2 = -0.5;
  
  ShiftAlongXaxis_TGraphErrors(grALICEDataPlusMinusSystematic,shift_ALICE);
  ShiftAlongXaxis_TGraphErrors(grALICEDataSameChargeSystematic,shift_ALICE);
  ShiftAlongXaxis_TGraphErrors(grALICEDataPlusMinus,shift_ALICE);
  ShiftAlongXaxis_TGraphErrors(grALICEDataSameCharge,shift_ALICE);

  ShiftAlongXaxis_TGraphAsymmErrors(grSTARDataPlusPlusSystematics,shift_STAR);
  ShiftAlongXaxis_TGraphAsymmErrors(grSTARDataPlusMinusSystematics,shift_STAR);
  ShiftAlongXaxis_TGraphErrors(grSTARDataPlusPlus,shift_STAR);
  ShiftAlongXaxis_TGraphErrors(grSTARDataPlusMinus,shift_STAR);
  
  ShiftAlongXaxis_TGraphErrors(grHIJINGWithFlowSameCharge,shift_HIJINGv2);
  ShiftAlongXaxis_TGraphErrors(grHIJINGWithFlowPlusMinus,shift_HIJINGv2);
  ShiftAlongXaxis_TGraphErrors(grHIJINGWithFlowAvg,shift_HIJINGv2);
  ShiftAlongXaxis_TGraphErrors(grHIJINGNoFlowSameCharge,shift_HIJING);
  ShiftAlongXaxis_TGraphErrors(grHIJINGNoFlowPlusMinus,shift_HIJING);
  ShiftAlongXaxis_TGraphErrors(grHIJINGNoFlowAvg,shift_HIJING);

  //grHIJINGWithFlowSameCharge->Draw(drawOptions_HIJING);
  //grHIJINGWithFlowPlusMinus->Draw(drawOptions_HIJING);
  //grHIJINGWithFlowAvg->Draw(drawOptions_HIJING);
  //grHIJINGNoFlowSameCharge->Draw(drawOptions_HIJING);
  //grHIJINGNoFlowPlusMinusAvg->Draw(drawOptions_HIJING);
  grHIJINGNoFlowAvg->Draw(drawOptions_HIJING);

  gHistv1Fluctuations->Draw("e3,same");
  grToneevSameCharge->Draw("C");

  grSTARDataPlusPlusSystematics->Draw("P2");
  grSTARDataPlusPlus->Draw("PZ");
  grSTARDataPlusMinusSystematics->Draw("P2");
  grSTARDataPlusMinus->Draw("PZ");

  grALICEDataPlusMinusSystematic->Draw("Z");
  grALICEDataSameChargeSystematic->Draw("Z");
  grALICEDataPlusMinus->Draw("PZ");
  TGraphErrors *grALICEDataPlusMinus_clone = (TGraphErrors*)grALICEDataPlusMinus->Clone("grALICEDataPlusMinus_clone");
  grALICEDataPlusMinus_clone->SetMarkerSize(0.9*grALICEDataPlusMinus->GetMarkerSize());
  TGraphErrors *grALICEDataPlusMinus_clone1 = (TGraphErrors*)grALICEDataPlusMinus->Clone("grALICEDataPlusMinus_clone1");
  grALICEDataPlusMinus_clone1->SetMarkerSize(0.8*grALICEDataPlusMinus->GetMarkerSize());
  grALICEDataPlusMinus_clone->Draw("PZ");
  grALICEDataPlusMinus_clone1->Draw("PZ");
  grALICEDataSameCharge->SetMarkerSize(1.1*grALICEDataSameCharge->GetMarkerSize());
  grALICEDataSameCharge->Draw("PZ");

  //TLegend *legend1 = new TLegend(0.18,0.54,0.8,0.85,"","brNDC");
  //TLegend *legend1 = new TLegend(0.19,0.63,0.78,0.88,"","brNDC");
  TLegend *legend1 = new TLegend(0.17,0.19,0.58,0.33,"","brNDC");
  myLegendSetUp(legend1,0.04);
  legend1->SetNColumns(2);
  //legend->AddEntry("Opposite charge","(+-)","");
  //legend->AddEntry("Same charge","Same Charge","");
  legend1->AddEntry(grALICEDataPlusMinusSystematic," ","L");
  legend1->AddEntry(grALICEDataSameChargeSystematic,"                                                       ","L");
//   legend1->AddEntry(grALICEDataSameChargeSystematic,"  ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV","lp");
  legend1->AddEntry(grSTARDataPlusMinusSystematics," ","L");
  legend1->AddEntry(grSTARDataPlusPlusSystematics,"                                                       ","L");
//   legend1->AddEntry("NULL"," ","");
//   legend1->AddEntry("NULL"," ","");
  legend1->AddEntry("NULL"," ","");
  legend1->AddEntry("NULL"," ","");
  legend1->AddEntry("NULL"," ","");
  legend1->AddEntry("NULL"," ","");
//   legend1->AddEntry("NULL"," ","");
//   legend1->AddEntry("NULL"," ","");
  //legend1->Draw();

  //TLegend *legend = new TLegend(0.19,0.68,0.78,0.93,"","brNDC");
  TLegend *legend = new TLegend(0.19,0.75,0.78,0.89,"","brNDC");
  myLegendSetUp(legend,0.04);
  legend->SetNColumns(2);
  //legend->AddEntry("Opposite charge","(+-)","");
  //legend->AddEntry("Same charge","Same Charge","");
  legend->AddEntry(grALICEDataSameCharge," ","p");
  legend->AddEntry(grALICEDataPlusMinus,"  ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV","p");
  legend->AddEntry(grSTARDataPlusPlus," ","p");
  legend->AddEntry(grSTARDataPlusMinus,"  STAR Au-Au @ #sqrt{s_{NN}} = 0.2 TeV","p");
  //legend->AddEntry(grHIJINGNoFlowSameCharge," ","p");
  //legend->AddEntry(grHIJINGNoFlowPlusMinus,"  #LTcos(#phi_{#alpha} + #phi_{#beta} - 2#phi_{c})#GT_{HIJING} / v_{2}{2}","p");
  //legend->AddEntry(grHIJINGWithFlowSameCharge," ","p");
  //legend->AddEntry(grHIJINGWithFlowPlusMinus,"  HIJING with v_{2} modulations","p");
//   legend->AddEntry(grToneevSameCharge," ","l");
//   legend->AddEntry("NULL","  CME expectation (Toneev #font[12]{et al.})","");
  legend->Draw();

  TLegend *legend3 = new TLegend(0.19,0.675,0.78,0.775,"","brNDC");
  myLegendSetUp(legend3,0.04);
  legend3->AddEntry(gHistv1Fluctuations,"  (ALICE) same+opp. mean","L");
  legend3->Draw();

  TLatex *myText = new TLatex();
  myText->SetNDC();
  myText->SetTextSize(0.04);
  myText->SetTextColor(1);
  myText->DrawLatex(0.195,0.89,"same");
  //myText->DrawLatex(0.155,0.34,"same");
  myText->SetTextColor(1);
  myText->DrawLatex(0.29,0.89,"opp.");
  //myText->DrawLatex(0.23,0.34,"opp.");
 
  //TLegend *legend2 = new TLegend(0.17,0.19,0.58,0.33,"","brNDC");
  TLegend *legend2 = new TLegend(0.17,0.19,0.58,0.37,"","brNDC");
  myLegendSetUp(legend2,0.04);
  //legend2->SetNColumns(2);
  //legend2->AddEntry("Opposite charge","(+-)","");
  //legend2->AddEntry("Same charge","Same Charge","");
//   legend2->AddEntry(grToneevSameCharge," Toneev #font[12]{et al.} (same charge)","L");
  legend2->AddEntry(grHIJINGNoFlowAvg,"  #LTcos(#phi_{#alpha} + #phi_{#beta} - 2#phi_{c})#GT_{HIJING} / v_{2}{2}","p");
  //legend2->AddEntry(grHIJINGWithFlowAvg,"  HIJING with v_{2} modulations","p");
  legend2->AddEntry(grToneevSameCharge," CME expectation (Toneev #font[12]{et al.})","L");
  //legend2->AddEntry(gHistv1Fluctuations," v_{1} fluctuations (same+opp. mean)","L");
  //legend2->AddEntry(gHistv1Fluctuations,"(ALICE) same+opp. mean","L");
  legend2->Draw();

  if(gPreliminary) {
    TLatex *alice = new TLatex(0.25,0.18,"Preliminary");
    alice->SetNDC();
    alice->SetTextColor(kRed+2);
    alice->SetTextSize(0.035);
    alice->SetLineWidth(2);
    alice->Draw();
    
    TPad *myPadLogo = new TPad("myPadLogo", 
			       "Pad for ALICE Logo",0.26,0.21,0.41,0.36);
    //myPadLogo->SetFillColor(2) // color to first figure out where is the pad then comment !
    myPadSetUp(myPadLogo,0,0,0,0);
    myPadLogo->Draw();
    myPadLogo->cd();
    TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
    myAliceLogo->Draw();
  }

  c1->SaveAs("figure2.eps");
  c1->SaveAs("figure2.pdf");
  c1->SaveAs("figure2.png");
}

//_______________________________________________________________//
void SetDataPoints(const char* resultsPath,
		   Int_t chargeCombination) {
  //chargeCombination == 0 ==> (+-)
  //chargeCombination == 1 ==> (++)
  //chargeCombination == -1 ==> (--)
  //chargeCombination == 2 ==> (--)&(++)

  TString filename = resultsPath;
  //if(chargeCombination == 0) filename += "TPCOnly/PlusMinus/outputMH.root";
  //else if(chargeCombination == 1) filename += "TPCOnly/PlusPlus/outputMH.root";
  //else if(chargeCombination == -1) filename += "TPCOnly/MinusMinus/outputMH.root"; 
  //else if(chargeCombination == 2) filename += "TPCOnly/SameCharge/outputMH.root"; 
  if(chargeCombination == 0) filename += "PlusMinus/TPCOnly/outputMH.root";
  else if(chargeCombination == 1) filename += "PlusPlus/TPCOnly/outputMH.root";
  else if(chargeCombination == -1) filename += "MinusMinus/TPCOnly/outputMH.root"; 
  else if(chargeCombination == 2) filename += "SameCharge/TPCOnly/outputMH.root"; 
  else {
    Printf("Wrong charge combinations selected - the supported values are 0, 1 and -1");
    break;
  }
  //Printf("%s",filename.Data());
  TFile *fInput = TFile::Open(filename.Data());
  if(!fInput) {
    Printf("File %s not found!!!",filename.Data());
    break;
  } 
  //fInput->ls();

  TString hist3pName;
  TH1D *gHist3p;
  TString histQC2Name;
  TH1D *gHistQC2;
  TString histQC4Name;
  TH1D *gHistQC4;

  
  Double_t g3pPlusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t g3pPlusPlusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t g3pMinusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t g3pSameChargeError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    hist3pName = "g3pHistName"; hist3pName += strCentralityBins[iBin];
    gHist3p = dynamic_cast<TH1D *>(fInput->Get(hist3pName.Data()));
    histQC2Name = "fHistIntegratedFlowRPQC_2;"; histQC2Name += iBin + 1;
    gHistQC2 = dynamic_cast<TH1D *>(fInput->Get(histQC2Name.Data()));
    histQC4Name = "fHistIntegratedFlowRPQC_4;"; histQC4Name += iBin + 1;
    gHistQC4 = dynamic_cast<TH1D *>(fInput->Get(histQC4Name.Data()));

    if(chargeCombination == 0) {
      g3pCorrelatorALICEDataPlusMinus[iBin] = gHist3p->GetBinContent(1);
      g3pPlusMinusError[iBin] = gHist3p->GetBinError(1);
      v2QC2PlusMinus[iBin] = gHistQC2->GetBinContent(1);
      v2QC2PlusMinusError[iBin] = gHistQC2->GetBinError(1);
      v2QC4PlusMinus[iBin] = gHistQC4->GetBinContent(1);
      v2QC4PlusMinusError[iBin] = gHistQC4->GetBinError(1);
    }
    if(chargeCombination == 1) {
      g3pCorrelatorALICEDataPlusPlus[iBin] = gHist3p->GetBinContent(1);
      g3pPlusPlusError[iBin] = gHist3p->GetBinError(1);
      v2QC2PlusPlus[iBin] = gHistQC2->GetBinContent(1);
      v2QC2PlusPlusError[iBin] = gHistQC2->GetBinError(1);
      v2QC4PlusPlus[iBin] = gHistQC4->GetBinContent(1);
      v2QC4PlusPlusError[iBin] = gHistQC4->GetBinError(1);
    }
    if(chargeCombination == -1) {
      g3pCorrelatorALICEDataMinusMinus[iBin] = gHist3p->GetBinContent(1);
      g3pMinusMinusError[iBin] = gHist3p->GetBinError(1);
      v2QC2MinusMinus[iBin] = gHistQC2->GetBinContent(1);
      v2QC2MinusMinusError[iBin] = gHistQC2->GetBinError(1);
      v2QC4MinusMinus[iBin] = gHistQC4->GetBinContent(1);
      v2QC4MinusMinusError[iBin] = gHistQC4->GetBinError(1);
    }
    if(chargeCombination == 2) {
      g3pCorrelatorALICEDataSameCharge[iBin] = gHist3p->GetBinContent(1);
      g3pSameChargeError[iBin] = gHist3p->GetBinError(1);
      v2QC2SameCharge[iBin] = gHistQC2->GetBinContent(1);
      v2QC2SameChargeError[iBin] = gHistQC2->GetBinError(1);
      v2QC4SameCharge[iBin] = gHistQC4->GetBinContent(1);
      v2QC4SameChargeError[iBin] = gHistQC4->GetBinError(1);
    }
  }//centrality bin

  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    if(chargeCombination == 0) {
      g3pCorrelatorALICEDataPlusMinus[iBin] /= v2QC2PlusMinus[iBin];
      g3pCorrelatorALICEDataPlusMinusError[iBin] = TMath::Sqrt(TMath::Power(g3pPlusMinusError[iBin]/v2QC2PlusMinus[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusMinus[iBin]*v2QC2PlusMinusError[iBin]/(v2QC2PlusMinus[iBin]*v2QC2PlusMinus[iBin]),2));
    }
    if(chargeCombination == 1) {
      g3pCorrelatorALICEDataPlusPlus[iBin] /= v2QC2PlusPlus[iBin];
      g3pCorrelatorALICEDataPlusPlusError[iBin] = TMath::Sqrt(TMath::Power(g3pPlusPlusError[iBin]/v2QC2PlusPlus[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusPlus[iBin]*v2QC2PlusPlusError[iBin]/(v2QC2PlusPlus[iBin]*v2QC2PlusPlus[iBin]),2));
    }
    if(chargeCombination == -1) {
      g3pCorrelatorALICEDataMinusMinus[iBin] /= v2QC2MinusMinus[iBin];
      g3pCorrelatorALICEDataMinusMinusError[iBin] = TMath::Sqrt(TMath::Power(g3pMinusMinusError[iBin]/v2QC2MinusMinus[iBin],2) + TMath::Power(g3pCorrelatorALICEDataMinusMinus[iBin]*v2QC2MinusMinusError[iBin]/(v2QC2MinusMinus[iBin]*v2QC2MinusMinus[iBin]),2));
    }
    if(chargeCombination == 2) {
      g3pCorrelatorALICEDataSameCharge[iBin] /= v2QC2SameCharge[iBin];
      g3pCorrelatorALICEDataSameChargeError[iBin] = TMath::Sqrt(TMath::Power(g3pSameChargeError[iBin]/v2QC2SameCharge[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameCharge[iBin]*v2QC2SameChargeError[iBin]/(v2QC2SameCharge[iBin]*v2QC2SameCharge[iBin]),2));
    }
  }
}

//_______________________________________________________________//
void GetCorrelatorAndError(TProfile *g3pCorrelatorVsPt,
			   Double_t &g3pCorrelatorValue,
			   Double_t &g3pCorrelatorError,
			   Int_t iBinLow = 0,
			   Int_t iBinHigh = 0) {
  //Function to return the average value of the 3p correlator 
  //<cos(psi1 + psi2 - 2phi3)> and its error.
  //The first argument is one of the 3p TProfile objects vs pt.
  //The second and third argument give the integral and its error.
  //The fourth and fifth, if specified, indicate the lowest and 
  //highest bin the calculation should be performed for.
  Int_t gBinLow = 1, gBinHigh = g3pCorrelatorVsPt->GetNbinsX();
  if(iBinLow) gBinLow = iBinLow;
  if(iBinHigh) gBinHigh = iBinHigh;
  
  Double_t gSumXi = 0.;
  Double_t gSumYi = 0.;
  Double_t gSumXiYi = 0.;
  Double_t gSumXiYi2 = 0.;
  Double_t gSumXi2Yi2 = 0.;
  Double_t gSumDeltaXi2 = 0.;
  Double_t gSumYi2DeltaXi2 = 0.;
  Double_t dError = 0.; //Flow code driven error calculation

  Double_t kSumBi = 0., kSumBi2DeltaYi2 = 0.;
  Double_t kSumYiBi = 0., kSumYi2DeltaBi2 = 0.;
  Double_t kSumDeltaBi2 = 0.;

  for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    gSumXi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumYi += g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi2 += g3pCorrelatorVsPt->GetBinEntries(iBin)*TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumXi2Yi2 += TMath::Power(g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumDeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
    gSumYi2DeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2) + TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);

    dError += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinError(iBin)*g3pCorrelatorVsPt->GetBinError(iBin);  

    //new error calculation
    kSumBi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    kSumYiBi += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin);
    kSumBi2DeltaYi2 += TMath::Power(g3pCorrelatorVsPt->GetBinEntries(iBin),2)*TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
    kSumYi2DeltaBi2 += TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2)*TMath::Power(TMath::Sqrt(g3pCorrelatorVsPt->GetBinEntries(iBin)),2);  
    kSumDeltaBi2 += TMath::Power(TMath::Sqrt(g3pCorrelatorVsPt->GetBinEntries(iBin)),2);
  }
  
  g3pCorrelatorValue = -1000.;
  g3pCorrelatorError = 1000.;
  
  if(gSumXi != 0.)
    g3pCorrelatorValue = gSumXiYi/gSumXi;
  if((gSumXi != 0.)&&(gSumXiYi != 0.))
    g3pCorrelatorError = TMath::Abs((gSumXiYi/gSumXi))*TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumYi2DeltaXi2)/gSumXiYi),2) + TMath::Power((gSumDeltaXi2/gSumXi),2));
    //g3pCorrelatorError = TMath::Sqrt((1./TMath::Power(kSumBi,2))*(kSumBi2DeltaYi2 + kSumYi2DeltaBi2 + TMath::Power(kSumYiBi,2)*kSumDeltaBi2/TMath::Power(kSumBi,2)));
  if(gSumXi != 0.)
    dError /= TMath::Power(gSumXi,2);
  dError = TMath::Sqrt(dError);
  g3pCorrelatorError = dError;

  /*Int_t iBinCounter = 0;
  Double_t gSumBinContentTimesWeight = 0., gSumWeight = 0.;
  Double_t gSumBinContentTimesWeightSquared = 0.;
  for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    iBinCounter += 1;
    
    gSumBinContentTimesWeight += g3pCorrelatorVsPt->GetBinContent(iBin)*g3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumWeight += g3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumBinContentTimesWeightSquared += TMath::Power(gSumBinContentTimesWeight,2);
  }

  //Printf("%lf - %d",gSumWeight,iBinCounter);
  //Calculate the g3pCorrelatorValue and its error
  g3pCorrelatorValue = -1000.;
  g3pCorrelatorError = 1000.;
  if((gSumWeight)&&(iBinCounter)) {
    g3pCorrelatorValue = gSumBinContentTimesWeight/(gSumWeight);
    g3pCorrelatorError = TMath::Sqrt(gSumBinContentTimesWeightSquared/TMath::Power(gSumWeight,2));
    }*/
  
  //+++++++++++++++++++Extra new: Proper error calculation++++++++++++++++//
  Double_t kSumBi = 0.;
  Double_t kSumAiBi = 0.;
  Double_t kSumBi2DeltaAi2 = 0.;
  Double_t kSumLong = 0.;
  Double_t k3pCorrelatorValue = 10000., k3pCorrelatorValueError = 10000.;

  for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    kSumBi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    kSumAiBi += g3pCorrelatorVsPt->GetBinContent(iBin)*g3pCorrelatorVsPt->GetBinEntries(iBin);
    kSumBi2DeltaAi2 += g3pCorrelatorVsPt->GetBinEntries(iBin)*TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
  }

  for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    kSumLong += TMath::Power((g3pCorrelatorVsPt->GetBinContent(iBin)*kSumBi - kSumAiBi),2)*TMath::Power(TMath::Sqrt( g3pCorrelatorVsPt->GetBinEntries(iBin)),2)/TMath::Power(kSumBi,2);
  }

  k3pCorrelatorValue = kSumAiBi/kSumBi;
  k3pCorrelatorValueError = TMath::Sqrt((1./TMath::Power(kSumBi,2))*(kSumBi2DeltaAi2 + kSumLong));

  //g3pCorrelatorValue = k3pCorrelatorValue;
  //g3pCorrelatorError = k3pCorrelatorValueError;

  return;
}


//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}

