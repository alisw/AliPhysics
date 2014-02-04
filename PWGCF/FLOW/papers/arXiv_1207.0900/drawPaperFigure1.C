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

// #include "SetFlowStyle.C"

Bool_t gPreliminary = kFALSE;

float myMarkerSize = 2.0;
static  int      myDarkRed  = TColor::GetColor(128,0,0);
static  int      myLightRed  = TColor::GetColor(128,0,0);
static  int      myBlue     = 9;
static  int      myGreen     = kGreen+3;

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 9;
TString strCentralityBins[nCentralityBins] = {"0-5","5-10","10-20",
					      "20-30","30-40","40-50",
					      "50-60","60-70","70-80"};
//Double_t gCentralityPercentile[nCentralityBins] = {75.,65.,55.,45.,35.,25.,15.,7.5,2.5};
//Double_t gCentralityPercentile[nCentralityBins] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
//Double_t gCentralityPercentileError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
Double_t gCentralityPercentileSystematicError[nCentralityBins] = {1.5,1.5,3.0,3.0,3.0,3.0,3.0,3.0,3.0};

Double_t gCentralityPercentileError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
Double_t gCentralityPercentileVZERO[nCentralityBins] = {3.,8.,16.,26.,36.,46.,56.,66.,76.};
Double_t gCentralityPercentileZDC[nCentralityBins] = {4.25,9.25,18.,28.,38.,48.,58.,68.,78.};
Double_t gCentralityPercentile[nCentralityBins] = {0.5,5.5,12.,22.,32.,42.,52.,62.,74.};
Double_t gCentralityPercentileTPC[nCentralityBins] = {1.75,6.75,14.,24.,34.,44.,54.,64.,74.};

//================================ALICE================================//
//================================Panos-VZERO================================//
Double_t g3pCorrelatorPanosCorrectedPlusMinus[nCentralityBins] = {1.09401e-07,-1.09401e-06,-1.41212e-05,-2.49265e-05,-6.4138e-05,-2.57942e-05,0.000107605,0.000314927,0.000989089};
Double_t g3pCorrelatorPanosCorrectedPlusMinusError[nCentralityBins] = {4.05503e-06,6.16342e-06,9.18652e-06,1.42343e-05,2.32727e-05,4.07319e-05,8.00524e-05,0.00012171819,0.0002171819};
Double_t g3pCorrelatorPanosCorrectedSameCharge[nCentralityBins] = {-1.09401e-07,-9.82253e-06,-4.64251e-05,-9.80566e-05,-0.000147225,-0.000248039,-0.000366027,-0.000370332,-0.000393409};
Double_t g3pCorrelatorPanosCorrectedSameChargeError[nCentralityBins] = {4.10413e-06,6.23425e-06,9.33561e-06,1.4481e-05,2.3662e-05,4.13592e-05,8.14661e-05,0.000176888,0.000376888};
//================================Panos-VZERO================================//

//================================Ilya-ZDC================================//
//Double_t g3pCorrelatorIlyaCorrectedPlusMinus[nCentralityBins] = {4.297763e-05,-2.408162e-05,-1.637881e-05,-4.876254e-06,-2.316564e-05,-7.809004e-05,6.207905e-06,-0.0004647401,-0.0001632693};  
Double_t g3pCorrelatorIlyaCorrectedPlusMinus[nCentralityBins] = {4.297763e+15,-2.408162e-05,-1.637881e-05,-4.876254e-06,-2.316564e-05,-7.809004e-05,6.207905e-06,-0.0004647401,-0.0001632693};  
Double_t g3pCorrelatorIlyaCorrectedPlusMinusError[nCentralityBins] = {3.5165e-05,2.527027e-05,1.044865e-05,1.15894e-05,1.778335e-05,3.373117e-05,7.60854e-05,0.000225265,0.0008951721};
//Double_t g3pCorrelatorIlyaCorrectedSameCharge[nCentralityBins] = {-5.446542e-05,-3.563536e-05,-5.561554e-05,-9.352835e-05,-0.0001713429,-0.0003199126,-0.0003791019,-0.0006610145,-0.00121334};
Double_t g3pCorrelatorIlyaCorrectedSameCharge[nCentralityBins] = {-5.446542e+15,-3.563536e-05,-5.561554e-05,-9.352835e-05,-0.0001713429,-0.0003199126,-0.0003791019,-0.0006610145,-0.00121334};
Double_t g3pCorrelatorIlyaCorrectedSameChargeError[nCentralityBins] = {3.556981e-05,2.425196e-05,1.01993e-05,1.134225e-05,1.746077e-05,3.325377e-05,7.538573e-05,0.0002243383,0.0008956909};
//================================Ilya-ZDC================================//

//==============================Alexandru-TPC==============================//
Double_t g3pCorrelatorAlexandruCorrectedPlusMinus[nCentralityBins] = {-6.41377e-06,1.87768e-06,-5.07504e-06,-1.13659e-05,-1.03051e-05,-3.0122e-05,-9.96488e-06,0.000246629,0.00108092};
Double_t g3pCorrelatorAlexandruCorrectedPlusMinusError[nCentralityBins] = {3.47328e-06,3.67446e-06,3.31881e-06,4.77789e-06,7.23694e-06,1.17143e-05,2.16551e-05,4.5914e-05,0.000123041};
Double_t g3pCorrelatorAlexandruCorrectedSameCharge[nCentralityBins] = {-1.03186e-05,-2.34825e-05,-5.40217e-05,-9.6981e-05,-0.000159935,-0.000268521,-0.000371931,-0.000445415,-0.000214926};
Double_t g3pCorrelatorAlexandruCorrectedSameChargeError[nCentralityBins] = {3.47395e-06,3.67532e-06,3.31981e-06,4.77997e-06,7.24162e-06,1.17259e-05,2.16903e-05,4.60457e-05,0.000123736};
//==============================Alexandru-TPC==============================//


//================================ALICE================================//
//<cos*cos>
Double_t gCosCosALICEDataSameCharge[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gCosCosALICEDataSameChargeError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gCosCosALICEDataSameChargeSystematicError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gCosCosALICEDataPlusMinus[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gCosCosALICEDataPlusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gCosCosALICEDataPlusMinusSystematicError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

//<sin*sin>
Double_t gSinSinALICEDataSameCharge[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gSinSinALICEDataSameChargeError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gSinSinALICEDataSameChargeSystematicError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gSinSinALICEDataPlusMinus[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gSinSinALICEDataPlusMinusError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gSinSinALICEDataPlusMinusSystematicError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
//================================ALICE================================//

//================================ALICE================================//
Double_t g2pCorrelatorALICEDataPlusMinus[nCentralityBins] = {0.000456523, 0.000560778, 0.000731902, 0.00102103, 0.00143198, 0.00212054, 0.00330624, 0.00544651, 0.00919686};
Double_t g2pCorrelatorALICEDataPlusMinusError[nCentralityBins] = {2.30045e-06,
								  2.81061e-06,
								  2.64844e-06,
								  3.9621e-06,
								  6.07271e-06,
								  9.91051e-06,
								  1.75048e-05,
								  3.37665e-05,
								  7.36336e-05};
Double_t g2pCorrelatorALICEDataPlusMinusSystematicError[nCentralityBins] = {1.35E-07,
                                                                            2.16E-07,
                                                                            1.55E-07,
                                                                            3.98E-07,
                                                                            1.54E-06,
                                                                            1.49E-06,
                                                                            2.90E-06,
                                                                            2.60E-06,
                                                                            2.27E-05};

Double_t g2pCorrelatorALICEDataSameCharge[nCentralityBins] = {0.000139404,
							      0.000177643,
							      0.000241073,
							      0.000356669,
							      0.000501588,
							      0.000756712,
							      0.00116784,
							      0.00209485,
							      0.00392973};
Double_t g2pCorrelatorALICEDataSameChargeError[nCentralityBins] = {1.63072e-06,
								   1.99034e-06,
								   1.87724e-06,
								   2.80842e-06,
								   4.29531e-06,
								   7.01318e-06,
								   1.24019e-05,
								   2.39888e-05,
								   5.25678e-05};
Double_t g2pCorrelatorALICEDataSameChargeSystematicError[nCentralityBins] = {1.38541E-07,
                                                                             2.82213E-07,
                                                                             6.48121E-07,
                                                                             1.35494E-06,
                                                                             2.78063E-06,
                                                                             6.51713E-06,
                                                                             7.92961E-06,
                                                                             5.92337E-06,
                                                                             6.47308E-05};
//================================ALICE================================//
Double_t g3pCorrelatorALICEDataSameCharge[nCentralityBins] = {-1.60885e-05,-2.33387e-05,-5.49677e-05,-0.000105864,-0.000166243,-0.000276424,-0.00037915,-0.000425016,-0.000293409}; 
Double_t g3pCorrelatorALICEDataSameChargeError[nCentralityBins] = {3.22359e-06,
3.45502e-06,3.26195e-06,4.85736e-06,7.60362e-06,1.31918e-05,2.61394e-05,6.2836e-05,0.000196171};
Double_t g3pCorrelatorALICEDataSameChargeSystematicError[nCentralityBins] = {3.58357E-05,1.66674E-05,1.27416E-05,1.98483E-05,3.97548E-05,8.97482E-05,0.00010582,
0.000233554,0.00031035126};

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
Double_t g3pCorrelatorALICEDataPlusMinusSystematicError[nCentralityBins] = {3.55452E-05, 1.8482E-05,
									    6.31586E-06,
									    1.20367E-05,
									    3.66273E-05,
									    4.12716E-05,
									    8.08092E-05,
									    0.000158143,
									    0.00031038405};

//================================STAR================================//
Double_t g2pCorrelatorSTARPlusMinus[nCentralityBins] = {0.000345933,
							    0.000426757,
							    0.000546398,
							    0.000743832,
							    0.00102597,
							    0.00140117,
							    0.00195394,
							    0.0027375,
							    0.00363606};
Double_t g2pCorrelatorSTARPlusMinusError[nCentralityBins] = {2.06688e-06,
							     2.5451e-06,
							     2.43789e-06,
							     3.24192e-06,
							     4.82798e-06,
							     7.43894e-06,
							     1.22702e-05,
							     2.25468e-05,
							     4.74887e-05};

Double_t g2pCorrelatorSTARPlusPlus[nCentralityBins] = {-0.000162705,
						       -0.000184491,
						       -0.000201566,
						       -0.00022444,
						       -0.000262823,
						       -0.000328381,
						       -0.000364104,
						       -0.000444511,
						       -0.000558773};
Double_t g2pCorrelatorSTARPlusPlusError[nCentralityBins] = {2.85059e-06,
							    3.50367e-06,
							    3.35781e-06,
							    4.46332e-06,
							    6.67312e-06,
							    1.03454e-05,
							    1.72965e-05,
							    3.21068e-05,
							    7.06495e-05};
//================================STAR================================//

//================================HIJING================================//
Double_t g2pCorrelatorHIJINGPlusMinus[nCentralityBins] = {-10.,0.00041372,0.000590745,0.000879697,0.00133003,0.00214314,0.00357934,0.00626945,0.0111598};
Double_t g2pCorrelatorHIJINGPlusMinusError[nCentralityBins] = {0.,1.0236e-06,8.68066e-07,1.33362e-06,2.56541e-06,3.58639e-06,6.34634e-06,1.22583e-05,2.49233e-05};

Double_t g2pCorrelatorHIJINGSameCharge[nCentralityBins] = {-10.,0.00038675,0.000515641,0.000779283,0.001169,0.00190631,0.00321236,0.00564697,0.0104565};
Double_t g2pCorrelatorHIJINGSameChargeError[nCentralityBins] = {0.,1.73833e-06,1.87754e-06,2.89076e-06,5.52952e-06,7.75963e-06,1.37514e-05,2.6596e-05,5.41739e-05};
//================================HIJING================================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void drawPaperFigure1() {
  //Draws the <cos*cos> and <sin*sin> vs centrality percentile
  //gROOT->LoadMacro("SetFlowStyle.C");
  //SetFlowStyle();
  gROOT->LoadMacro("SetPlotStyle.C");
  SetPlotStyle();

  TGaxis::SetMaxDigits(5);

  Double_t arrOpp[4] = {g3pCorrelatorALICEDataPlusMinus[0],g3pCorrelatorPanosCorrectedPlusMinus[0],g3pCorrelatorIlyaCorrectedPlusMinus[0],g3pCorrelatorAlexandruCorrectedPlusMinus[0]};
  Double_t arrSame[4] = {g3pCorrelatorALICEDataSameCharge[0],g3pCorrelatorPanosCorrectedSameCharge[0],g3pCorrelatorIlyaCorrectedSameCharge[0],g3pCorrelatorAlexandruCorrectedSameCharge[0]};
  Printf("(Same charge) Systematic: %3.10lf",TMath::Abs(arrSame[TMath::LocMax(4,arrSame)]-arrSame[TMath::LocMin(3,arrSame)])/2.);
  Printf("(Opp charge) Systematic: %3.10lf",TMath::Abs(arrOpp[TMath::LocMax(4,arrOpp)]-arrOpp[TMath::LocMin(3,arrOpp)])/2.);
  
  //Calculate the coscos and sinsin terms and their errors
  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    //coscos (same charge)
    gCosCosALICEDataSameCharge[iBin] = 0.5*(g2pCorrelatorALICEDataSameCharge[iBin] + g3pCorrelatorALICEDataSameCharge[iBin]);
    gCosCosALICEDataSameChargeError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataSameChargeError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameChargeError[iBin],2));
    gCosCosALICEDataSameChargeSystematicError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataSameChargeSystematicError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameChargeSystematicError[iBin],2));

    //coscos (opposite charge)
    gCosCosALICEDataPlusMinus[iBin] = 0.5*(g2pCorrelatorALICEDataPlusMinus[iBin] + g3pCorrelatorALICEDataPlusMinus[iBin]);
    gCosCosALICEDataPlusMinusError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataPlusMinusError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusMinusError[iBin],2));
    gCosCosALICEDataPlusMinusSystematicError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataPlusMinusSystematicError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusMinusSystematicError[iBin],2));

    //sinsin (same charge)
    gSinSinALICEDataSameCharge[iBin] = 0.5*(g2pCorrelatorALICEDataSameCharge[iBin] - g3pCorrelatorALICEDataSameCharge[iBin]);
    gSinSinALICEDataSameChargeError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataSameChargeError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameChargeError[iBin],2));
    gSinSinALICEDataSameChargeSystematicError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataSameChargeSystematicError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataSameChargeSystematicError[iBin],2));

    //sinsin (opposite charge)
    gSinSinALICEDataPlusMinus[iBin] = 0.5*(g2pCorrelatorALICEDataPlusMinus[iBin] - g3pCorrelatorALICEDataPlusMinus[iBin]);
    gSinSinALICEDataPlusMinusError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataPlusMinusError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusMinusError[iBin],2));
    gSinSinALICEDataPlusMinusSystematicError[iBin] = 0.5*TMath::Sqrt(TMath::Power(g2pCorrelatorALICEDataPlusMinusSystematicError[iBin],2) + TMath::Power(g3pCorrelatorALICEDataPlusMinusSystematicError[iBin],2));
    
  }

  //================================================//
  //(+-)
  TGraphErrors *grALICEDataPlusMinus = new TGraphErrors(nCentralityBins,
							gCentralityPercentile,
							g3pCorrelatorALICEDataPlusMinus,
							gCentralityPercentileError,
							g3pCorrelatorALICEDataPlusMinusError);
  myTGraphSetUp(grALICEDataPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);

  TGraphErrors *grAlexandruCorrectedPlusMinus = new TGraphErrors(nCentralityBins,
								 gCentralityPercentileTPC,
								 g3pCorrelatorAlexandruCorrectedPlusMinus,
								 gCentralityPercentileError,
								 g3pCorrelatorAlexandruCorrectedPlusMinusError);
  myTGraphSetUp(grAlexandruCorrectedPlusMinus,25,myGreen,myMarkerSize,1,myGreen,2,1001,myGreen);

  TGraphErrors *grPanosCorrectedPlusMinus = new TGraphErrors(nCentralityBins,
							       gCentralityPercentileVZERO,
							       g3pCorrelatorPanosCorrectedPlusMinus,
							       gCentralityPercentileError,
							       g3pCorrelatorPanosCorrectedPlusMinusError);
  myTGraphSetUp(grPanosCorrectedPlusMinus,26,myBlue,myMarkerSize,1,myBlue,2,1001,myBlue);

  TGraphErrors *grIlyaCorrectedPlusMinus = new TGraphErrors(nCentralityBins,
							    gCentralityPercentileZDC,
							    g3pCorrelatorIlyaCorrectedPlusMinus,
							    gCentralityPercentileError,
							    g3pCorrelatorIlyaCorrectedPlusMinusError);
  myTGraphSetUp(grIlyaCorrectedPlusMinus,30,1,myMarkerSize,1,1,2,1001,1);

  //================================================//
  //(--)&(++)
  TGraphErrors *grALICEDataSameCharge = new TGraphErrors(nCentralityBins,
							 gCentralityPercentile,
							 g3pCorrelatorALICEDataSameCharge,
							 gCentralityPercentileError,
							 g3pCorrelatorALICEDataSameChargeError);
  myTGraphSetUp(grALICEDataSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  
  TGraphErrors *grAlexandruCorrectedSameCharge = new TGraphErrors(nCentralityBins,
								  gCentralityPercentileTPC,
								  g3pCorrelatorAlexandruCorrectedSameCharge,
								  gCentralityPercentileError,
								  g3pCorrelatorAlexandruCorrectedSameChargeError);
  myTGraphSetUp(grAlexandruCorrectedSameCharge,21,myGreen,myMarkerSize,1,myGreen,2,1001,myGreen);
  
  TGraphErrors *grPanosCorrectedSameCharge = new TGraphErrors(nCentralityBins,
								gCentralityPercentileVZERO,
								g3pCorrelatorPanosCorrectedSameCharge,
								gCentralityPercentileError,
								g3pCorrelatorPanosCorrectedSameChargeError);
  myTGraphSetUp(grPanosCorrectedSameCharge,22,myBlue,myMarkerSize,1,myBlue,2,1001,myBlue);

  TGraphErrors *grIlyaCorrectedSameCharge = new TGraphErrors(nCentralityBins,
							     gCentralityPercentileZDC,
							     g3pCorrelatorIlyaCorrectedSameCharge,
							     gCentralityPercentileError,
							     g3pCorrelatorIlyaCorrectedSameChargeError);
  myTGraphSetUp(grIlyaCorrectedSameCharge,29,1,myMarkerSize,1,1,2,1001,1);

  //================================================//
  //STAR
  TGraphErrors *gr2pSTARDataPlusMinus = new TGraphErrors(nCentralityBins,
							 gCentralityPercentileVZERO,
							 g2pCorrelatorSTARPlusMinus,
							 gCentralityPercentileError,
							 g2pCorrelatorSTARPlusMinusError);
  myTGraphSetUp(gr2pSTARDataPlusMinus,30,myGreen,myMarkerSize+0.5,1,myGreen,2,1001,myBlue);
  
  TGraphErrors *gr2pSTARDataPlusPlus = new TGraphErrors(nCentralityBins,
							gCentralityPercentileVZERO,
							g2pCorrelatorSTARPlusPlus,
							gCentralityPercentileError,
							g2pCorrelatorSTARPlusPlusError);
  myTGraphSetUp(gr2pSTARDataPlusPlus,29,myGreen,myMarkerSize+0.5,1,myDarkRed,2,1001,myGreen);
  //================================================//

  //================================================//
  //(+-)
  TGraphErrors *gr2pALICEDataPlusMinus = new TGraphErrors(nCentralityBins-1,
							  gCentralityPercentileVZERO,
							  g2pCorrelatorALICEDataPlusMinus,
							  gCentralityPercentileError,
							  g2pCorrelatorALICEDataPlusMinusError);
  myTGraphSetUp(gr2pALICEDataPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,5,1001,myBlue);
  
  TGraphErrors *gr2pALICEDataPlusMinusSystematic = new TGraphErrors(nCentralityBins-1,
								    gCentralityPercentileVZERO,
								    g2pCorrelatorALICEDataPlusMinus,
								    gCentralityPercentileSystematicError,
								    g2pCorrelatorALICEDataPlusMinusSystematicError);
  myTGraphSetUp(gr2pALICEDataPlusMinusSystematic,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myBlue);

  //(++)&(--)
  TGraphErrors *gr2pALICEDataSameCharge = new TGraphErrors(nCentralityBins-1,
							   gCentralityPercentileVZERO,
							   g2pCorrelatorALICEDataSameCharge,
							   gCentralityPercentileError,
							   g2pCorrelatorALICEDataSameChargeError);
  myTGraphSetUp(gr2pALICEDataSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,5,1001,myDarkRed);
  
  TGraphErrors *gr2pALICEDataSameChargeSystematic = new TGraphErrors(nCentralityBins-1,
								     gCentralityPercentileVZERO,
								     g2pCorrelatorALICEDataSameCharge,
								     gCentralityPercentileSystematicError,
								     g2pCorrelatorALICEDataSameChargeSystematicError);
  myTGraphSetUp(gr2pALICEDataSameChargeSystematic,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  //================================================//

  //================================================//
  //(+-) 
  TGraphErrors *grCosCosALICEDataPlusMinus = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gCosCosALICEDataPlusMinus,gCentralityPercentileError,gCosCosALICEDataPlusMinusError);
  myTGraphSetUp(grCosCosALICEDataPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,5,1001,myBlue);

  TGraphErrors *grCosCosALICEDataPlusMinusSystematic = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gCosCosALICEDataPlusMinus,gCentralityPercentileError,gCosCosALICEDataPlusMinusSystematicError);
  myTGraphSetUp(grCosCosALICEDataPlusMinusSystematic,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myBlue);

  TGraphErrors *grSinSinALICEDataPlusMinus = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gSinSinALICEDataPlusMinus,gCentralityPercentileError,gSinSinALICEDataPlusMinusError);
  myTGraphSetUp(grSinSinALICEDataPlusMinus,25,myBlue,myMarkerSize,1,myBlue,5,1001,myBlue);

  TGraphErrors *grSinSinALICEDataPlusMinusSystematic = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gSinSinALICEDataPlusMinus,gCentralityPercentileError,gSinSinALICEDataPlusMinusSystematicError);
  myTGraphSetUp(grSinSinALICEDataPlusMinusSystematic,25,myBlue,myMarkerSize,1,myBlue,2,1001,myBlue);


  //================================================//
  //(++)&(--)
  TGraphErrors *grCosCosALICEDataSameCharge = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gCosCosALICEDataSameCharge,gCentralityPercentileError,gCosCosALICEDataSameChargeError);
  myTGraphSetUp(grCosCosALICEDataSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,5,1001,myBlue);

  TGraphErrors *grCosCosALICEDataSameChargeSystematic = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gCosCosALICEDataSameCharge,gCentralityPercentileError,gCosCosALICEDataSameChargeSystematicError);
  myTGraphSetUp(grCosCosALICEDataSameChargeSystematic,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myBlue);

  TGraphErrors *grSinSinALICEDataSameCharge = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gSinSinALICEDataSameCharge,gCentralityPercentileError,gSinSinALICEDataSameChargeError);
  myTGraphSetUp(grSinSinALICEDataSameCharge,21,myBlue,myMarkerSize,1,myBlue,5,1001,myBlue);

  TGraphErrors *grSinSinALICEDataSameChargeSystematic = new TGraphErrors(nCentralityBins,gCentralityPercentileVZERO,gSinSinALICEDataSameCharge,gCentralityPercentileError,gSinSinALICEDataSameChargeSystematicError);
  myTGraphSetUp(grSinSinALICEDataSameChargeSystematic,21,myBlue,myMarkerSize,1,myBlue,2,1001,myBlue);
  
  //================================================//
  //HIJING (+-)
  TGraphErrors *gr2pHIJINGPlusMinus = new TGraphErrors(nCentralityBins-1,
						       gCentralityPercentileVZERO,
						       g2pCorrelatorHIJINGPlusMinus,
						       gCentralityPercentileError,
						       g2pCorrelatorHIJINGPlusMinusError);
  myTGraphSetUp(gr2pHIJINGPlusMinus,26,myBlue,myMarkerSize,1,myBlue,5,1001,myBlue);

  //(++)&(--)
  TGraphErrors *gr2pHIJINGSameCharge = new TGraphErrors(nCentralityBins-1,
							gCentralityPercentileVZERO,
							g2pCorrelatorHIJINGSameCharge,
							gCentralityPercentileError,
							g2pCorrelatorHIJINGSameChargeError);
  myTGraphSetUp(gr2pHIJINGSameCharge,22,myBlue,myMarkerSize,1,myBlue,5,1001,myBlue);
  //================================================//

  //_____________________________________________________//
  //Draw the results

  //====================================//
  //<cos(psi1+psi2-2phi3)> vs centrality
  //TH2F *gEmpty1 = new TH2F("gEmpty1",";centrality, %;#LT cos(#Delta #phi_{#alpha}) cos(#Delta #phi_{#beta}) #GT", nCentralityBins,0,80,1000,-1.5e-03,0.01);
  TH2F *gEmpty1 = new TH2F("gEmpty1",";centrality, %;", nCentralityBins,0,72,1000,-1.5e-03,0.015);
  gEmpty1->SetStats(kFALSE);
  gEmpty1->GetYaxis()->SetTitleSize(0.07);
  gEmpty1->GetYaxis()->SetTitleOffset(0.95);
  gEmpty1->GetYaxis()->SetNdivisions(10);
  gEmpty1->GetXaxis()->SetNdivisions(10);

  TF1 *f1 = new TF1("f1","0",0,1000);
  f1->SetLineColor(1); f1->SetLineStyle(1); f1->SetLineWidth(1);

  //================================================//
  TCanvas *c1 = new TCanvas("c1","Centrality dependence",0,0,700,1200);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->Divide(1,3,0.99,0.0,10);
  //TPad *myPad1 = new TPad("myPad1","myPad1",0,0.667,1,1);
  //c1->cd();
  //myPadSetUp(myPad1,0.13,0.07,0.04,0.0,10);
  //myPad1->Draw();
  c1->cd(1)->SetLeftMargin(0.19);
  c1->cd(1)->SetTopMargin(0.083);
  c1->cd(1)->SetRightMargin(0.01);
  //myPad1->cd();
  gEmpty1->GetYaxis()->SetLabelSize(0.075);
  gEmpty1->GetXaxis()->SetLabelSize(0.075);
  gEmpty1->GetYaxis()->SetTitleSize(0.095);
  gEmpty1->GetYaxis()->SetTitleOffset(1.01);
  gEmpty1->GetXaxis()->SetTitleSize(0.075);
  gEmpty1->GetYaxis()->SetNdivisions(4);
  gEmpty1->GetXaxis()->SetNdivisions(0);
  gEmpty1->GetYaxis()->SetRangeUser(-9.e-04,5.e-04);
  gEmpty1->GetYaxis()->SetTitle("#LT cos(#phi_{#alpha} + #phi_{#beta} - 2#Psi_{RP}) #GT");
  //gEmpty1->GetYaxis()->SetNdivisions(3);
  gEmpty1->DrawCopy();
  f1->Draw("same");
  grAlexandruCorrectedPlusMinus->Draw("P,Z");
  grAlexandruCorrectedSameCharge->Draw("P,Z");
  grIlyaCorrectedPlusMinus->Draw("P,Z");
  grIlyaCorrectedSameCharge->Draw("eZ,P");
  grPanosCorrectedPlusMinus->Draw("eZ,P");
  grPanosCorrectedSameCharge->Draw("eZ,P");
  grALICEDataPlusMinus->Draw("eZ,P");
  grALICEDataSameCharge->Draw("eZ,P");

  TLegend *legend1 = new TLegend(0.23,0.03,0.85,0.35,"","brNDC");
  myLegendSetUp(legend1,0.06);
  legend1->SetNColumns(2);
  legend1->AddEntry(grALICEDataSameCharge,"      ","lp");
  legend1->AddEntry(grALICEDataPlusMinus,"    TPC (cumulants)                                    ","lp");
  legend1->AddEntry(grAlexandruCorrectedSameCharge,"      ","lp");
  legend1->AddEntry(grAlexandruCorrectedPlusMinus,"    TPC","lp");
  legend1->AddEntry(grPanosCorrectedSameCharge,"      ","lp");
  legend1->AddEntry(grPanosCorrectedPlusMinus,"    VZERO","lp");
  legend1->AddEntry(grIlyaCorrectedSameCharge,"      ","lp");
  legend1->AddEntry(grIlyaCorrectedPlusMinus,"    ZDC","lp");
  legend1->Draw();
  TLatex *myText1 = new TLatex();
  myText1->SetNDC();
  myText1->SetTextSize(0.07);
  myText1->SetTextColor(1);
  myText1->DrawLatex(0.23,0.37,"same");
  myText1->SetTextColor(1);
  myText1->DrawLatex(0.35,0.37,"opp.");
  myText1->SetTextColor(1);
  myText1->DrawLatex(0.37,0.85,"ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV");
  myText1->DrawLatex(0.285,0.85,"(a)");

  //================================================//
  //TPad *myPad2 = new TPad("myPad2","myPad2",0,0.333,1,0.667);
  //c1->cd();
  //myPadSetUp(myPad2,0.13,0.00,0.04,0.00,10);
  //myPad2->Draw();
  //c1->cd();
  //myPad2->cd();
  c1->cd(2)->SetLeftMargin(0.19); c1->cd(2)->SetRightMargin(0.01);
  gEmpty1->GetYaxis()->SetLabelSize(0.075);
  gEmpty1->GetXaxis()->SetLabelSize(0.075);
  gEmpty1->GetYaxis()->SetTitleSize(0.095);
  gEmpty1->GetXaxis()->SetTitleSize(0.075);
  gEmpty1->GetXaxis()->SetNdivisions(0);
  gEmpty1->GetYaxis()->SetNdivisions(5);
  gEmpty1->GetYaxis()->SetTitle("#LT cos(#phi_{#alpha}-#phi_{#beta}) #GT");
  gEmpty1->GetYaxis()->CenterTitle();
  gEmpty1->GetYaxis()->SetRangeUser(-1.5e-03,0.007);
  gEmpty1->DrawCopy();
  f1->Draw("same");
  gr2pSTARDataPlusMinus->Draw("P,eZ");
  gr2pSTARDataPlusPlus->Draw("P,eZ");
  gr2pALICEDataPlusMinusSystematic->Draw("E3");
  gr2pALICEDataPlusMinusSystematic->Draw("L,X0,same");
  gr2pALICEDataPlusMinus->Draw("P,eZ");
  gr2pALICEDataSameChargeSystematic->Draw("E3");
  gr2pALICEDataSameChargeSystematic->Draw("L,X0,same");
  gr2pALICEDataSameCharge->Draw("P,eZ");

  //gr2pHIJINGPlusMinus->Draw("P,eZ");
  //gr2pHIJINGSameCharge->Draw("P,eZ");

  TLatex *myText2 = new TLatex();
  myText2->SetNDC();
  myText2->SetTextSize(0.07);
  myText2->SetTextColor(1);
  myText2->DrawLatex(0.22,0.86,"same");
  myText2->SetTextColor(1);
  myText2->DrawLatex(0.325,0.86,"opp.");
  myText2->DrawLatex(0.285,0.93,"(b)");

  TLegend *legend2 = new TLegend(0.22,0.63,0.89,0.83,"","brNDC");
  myLegendSetUp(legend2,0.06);
  legend2->SetNColumns(2);
  legend2->AddEntry(gr2pALICEDataSameCharge," ","PL");
  legend2->AddEntry(gr2pALICEDataPlusMinus,"  ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV","PL");
  legend2->AddEntry(gr2pSTARDataPlusPlus," ","P");
  legend2->AddEntry(gr2pSTARDataPlusMinus,"  STAR Au-Au @ #sqrt{s_{NN}} = 0.2 TeV","P");
  legend2->Draw();

  //================================================//
  //cos*cos & sin*sin
  //TPad *myPad3 = new TPad("myPad3","myPad3",0,0,1,0.333);
  //c1->cd();
  //myPadSetUp(myPad3,0.13,0.00,0.04,0.15,10);
  //myPad3->Draw();
  //c1->cd();
  //myPad3->cd();
  //gEmpty1->GetYaxis()->SetTitle("#LT cos(#phi_{#alpha}-#phi_{#beta}) #GT #pm #LT cos(#phi_{#alpha} + #phi_{#beta} - 2#Psi_{RP}) #GT");
  c1->cd(3)->SetLeftMargin(0.19); c1->cd(3)->SetRightMargin(0.01);
  gEmpty1->GetYaxis()->SetLabelSize(0.065);
  gEmpty1->GetXaxis()->SetLabelSize(0.065);
  gEmpty1->GetYaxis()->SetTitleSize(0.065);
  gEmpty1->GetXaxis()->SetTitleSize(0.075);
  gEmpty1->GetXaxis()->SetNdivisions(10);
  gEmpty1->GetYaxis()->SetTitle("");
  gEmpty1->GetYaxis()->SetRangeUser(0.0,0.0035);
  gEmpty1->GetYaxis()->SetNdivisions(4);
  gEmpty1->DrawCopy();
  f1->Draw("same");
  gStyle->SetErrorX(0);
  
  grCosCosALICEDataSameCharge->Draw("P,eZ");
  grCosCosALICEDataSameChargeSystematic->Draw("Z");
  grCosCosALICEDataPlusMinus->Draw("P,eZ");
  grCosCosALICEDataPlusMinusSystematic->Draw("Z");
  grSinSinALICEDataSameCharge->Draw("P,eZ");
  grSinSinALICEDataSameChargeSystematic->Draw("Z");
  grSinSinALICEDataPlusMinus->Draw("P,eZ");
  grSinSinALICEDataPlusMinusSystematic->Draw("Z");

  TLatex *myText3 = new TLatex();
  myText3->SetNDC();
  myText3->SetTextSize(0.055);
  myText3->SetTextColor(1);
  myText3->DrawLatex(0.22,0.86,"same");
  myText3->SetTextColor(1);
  myText3->DrawLatex(0.325,0.86,"opp.");
  myText3->DrawLatex(0.285,0.93,"(c)");
  myText3->DrawLatex(0.37,0.93,"ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV");

  TLegend *legend3 = new TLegend(0.22,0.63,0.89,0.83,"","brNDC");
  myLegendSetUp(legend3,0.06);
  legend3->SetNColumns(2);
  legend3->AddEntry(grCosCosALICEDataSameCharge," ","P");
  legend3->AddEntry(grCosCosALICEDataPlusMinus,"  #LT cos(#Delta #phi_{#alpha}) cos(#Delta #phi_{#beta}) #GT","P");
  legend3->AddEntry(grSinSinALICEDataSameCharge," ","P");
  legend3->AddEntry(grSinSinALICEDataPlusMinus,"  #LT sin(#Delta #phi_{#alpha}) sin(#Delta #phi_{#beta}) #GT","P");
  legend3->Draw();

 if(gPreliminary) {
    TLatex *alice = new TLatex(0.23,0.38,"Preliminary");
    alice->SetNDC();
    alice->SetTextColor(myDarkRed);
    alice->SetTextSize(0.035);
    alice->SetLineWidth(2);
    alice->Draw();
    
    TPad *myPadLogo = new TPad("myPadLogo", 
			       "Pad for ALICE Logo",0.26,0.41,0.36,0.51);
    //myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
    myPadSetUp(myPadLogo,0,0,0,0);
    myPadLogo->Draw();
    myPadLogo->cd();
    TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
    myAliceLogo->Draw();
  }

 c1->SaveAs("figure1.eps");
 c1->SaveAs("figure1.pdf");
 c1->SaveAs("figure1.png");
}

void SetPadSetUp(TPad *currentPad, 
		 float currentLeft=0.11, float currentTop=0.04, 
		 float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  currentPad->SetFillColor(10);

  return;
}

void myTGraphSetUp(TGraphErrors *currentGraph=0,
		   int myMarkerStyle=8,
		   int myMarkerColor=1,
		   float myMarkerSize=1,
		   int myLineStyle=1,
		   int myLineColor=1,
		   float myLineWidth=1,
		   int myFillStyle =1001,
		   int myFillColor =1) {
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
