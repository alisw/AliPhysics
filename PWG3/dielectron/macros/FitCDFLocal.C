TH1F *psproperMCsecJPSI();
void LoadLib();

void FitCDFLocal(){

  ///////////////////////////////////////////////////////////////////
  //
  // Example macro to read local N-Tuples of JPSI 
  // and bkg candidates and perform log-likelihood 
  // minimization using the minimization handler class
  //
  // Origin: C. Di Giglio
  //
  ///////////////////////////////////////////////////////////////////
  LoadLib();
  AliDielectronBtoJPSItoEleCDFfitFCN *likely_obj = 0x0;
  
  // Those initial parameters are evaluated from the fits done separately on
  // invariant mass and pseudoproper distributions. For the moment are setted by
  // hand 

  Double_t *resParam = new Double_t[6];
  resParam[0] = 2586.32; resParam[1]= -0.163616; resParam[2] = 31.026770; resParam[3] = 433.85; 
  resParam[4] = 4.67; resParam[5] = 73.31;

  Double_t *bkgParam  = new Double_t[7];
  bkgParam[0] = 7063.8; bkgParam[1] = 23214.45; bkgParam[2] = 24742.65; bkgParam[3] = 15188.57; 
  bkgParam[4] = 0.009579; bkgParam[5] = 0.010422; bkgParam[6] = 0.000423;

  Double_t *invMassParam = new Double_t[9];
  invMassParam[0] = 3.092; invMassParam[1] =  0.0306; invMassParam[2] =  0.4939; 
  invMassParam[3] = 1.36; invMassParam[4] = 41.97; invMassParam[5] = 80.43; invMassParam[6] = 0.975;
  invMassParam[7] = 1.61; invMassParam[8] = 2.129;
  
  Double_t Fsig = 0.454;

  TH1F* hCsiMCPithya = new TH1F();

  hCsiMCPithya = psproperMCsecJPSI();  
  Double_t integral = 0;
  for(int i=1;i<hCsiMCPithya->GetNbinsX()+1; i++) integral += (hCsiMCPithya->GetBinContent(i)*hCsiMCPithya->GetBinWidth(i));
  hCsiMCPithya->Scale(1./integral);

  Double_t* x=0x0; Double_t* m=0x0; Int_t n=0;
  AliDielectronBtoJPSItoEle* aBtoJPSItoEle =new AliDielectronBtoJPSItoEle();
  Double_t paramInputValues [20] =   
                                   /*paramInputValues[0]   ----> fWeightRes;
                                     paramInputValues[1]   ----> fPos;
                                     paramInputValues[2]   ----> fNeg;
                                     paramInputValues[3]   ----> fSym;
                                     paramInputValues[4]   ----> fOneOvLamPlus;
                                     paramInputValues[5]   ----> fOneOvLamMinus;
                                     paramInputValues[6]   ----> fOneOvLamSym;
                                     paramInputValues[7]   ----> fB;
                                     paramInputValues[8]   ----> fFsig;
                                     paramInputValues[9]   ----> fMmean;
                                     paramInputValues[10]  ----> fNexp;
                                     paramInputValues[11]  ----> fSigma;
                                     paramInputValues[12]  ----> fAlpha;
                                     paramInputValues[13]  ----> fNorm;
                                     paramInputValues[14]  ----> fBkgNorm;
                                     paramInputValues[15]  ----> fBkgMean;
                                     paramInputValues[16]  ----> fBkgSlope;
                                     paramInputValues[17]  ----> fBkgConst;
                                     paramInputValues[18]  ----> fGaus1Norm; 
                                     paramInputValues[19]  ----> fGaus2Norm;*/                                     
                                   { bkgParam[9],  // fWeightRes [0]
                                     bkgParam[10], // Fpos [1]
                                     bkgParam[11], // FNeg [2]
                                     bkgParam[12], // FSym [3]
                                     bkgParam[6], // LamdaPos [4]
                                     bkgParam[7], // LambdaNeg [5]
                                     bkgParam[8], // LamdaSym [6]
                                     0.10, // fB [7]
                                     Fsig, // FSig [8]
                                     invMassParam[0], // FMean [9]
                                     invMassParam[4], // fNexp [10]
                                     invMassParam[1], // fNsigma [11]
                                     invMassParam[2], // fAlpha [12]
                                     invMassParam[3], //fNorm [13]
                                     invMassParam[5], // fBkgNorm [14]
                                     invMassParam[6], // fBkgMean [15]
                                     invMassParam[7], // fBkgSlope [16]
                                     invMassParam[8], // fBkgConst [17]
                                     resParam[0], // Gaus1Norm [18]
                                     resParam[3] // Gaus2Norm [19]
                                      }; // Starting values for parameters 

  // retrieve the TNtupla from the file result.root
  TFile f("result.root");
  TList * list = (TList*)f.Get("resultAOD");
  TNtuple *nt=(TNtuple*)list->FindObject("fNtupleJPSI");
  nt->ls();
  aBtoJPSItoEle->ReadCandidates(nt,x,m,n); // read N-Tuples
  printf("+++\n+++ Number of total candidates (prim J/psi + secondary J/psi + bkg) ---> %d candidates \n+++\n",n);

  aBtoJPSItoEle->SetFitHandler(x,m,n); // Set the fit handler with given values of x, m, # of candidates 

  aBtoJPSItoEle->CloneMCtemplate(hCsiMCPithya);    // clone MC template and copy internally
                                                   // in this way any model can be setted from outside

  aBtoJPSItoEle->SetCsiMC(); // Pass the MC template to the CDF fit function

  AliDielectronBtoJPSItoEleCDFfitHandler* fithandler = aBtoJPSItoEle->GetCDFFitHandler(); // Get the fit handler

  //
  // Set some fit options through the handler class
  //
  fithandler->SetErrorDef(0.5); // tells Minuit that the error interval is the one in which
                                // the function differs from the minimum for less than setted value

  fithandler->SetPrintStatus(kTRUE);

  fithandler->SetParamStartValues(paramInputValues);

  Double_t *resConst = new Double_t[4];
  
  //resolution constants: fixed from MC
  resConst[0]= resParam[1]; // mean1
  resConst[1]= resParam[2]; // sigma1
  resConst[2]= resParam[4]; // mean2
  resConst[3]= resParam[5]; // sigma2

  
  fithandler->SetResolutionConstants(resConst); 
  fithandler->SetCrystalBallFunction(kTRUE);
  
  // fix parameter
  fithandler->FixParam(0,kTRUE);
  fithandler->FixParam(1,kTRUE); // Fix Bkg  weights  
  fithandler->FixParam(2,kTRUE);
  fithandler->FixParam(3,kTRUE);

  fithandler->FixParam(4,kTRUE);
  fithandler->FixParam(5,kTRUE); // Fix bkg exponential 
  fithandler->FixParam(6,kTRUE);
  //fithandler->FixParam(8,kTRUE); // Fsig fix
  
  fithandler->FixParam(9,kTRUE);
  fithandler->FixParam(10,kTRUE); // Fix CristalBall param
  fithandler->FixParam(11,kTRUE);
  fithandler->FixParam(12,kTRUE);  
  fithandler->FixParam(13,kTRUE); 

  fithandler->FixParam(14,kTRUE);
  fithandler->FixParam(15,kTRUE); // Fix Bkg invMass param
  fithandler->FixParam(16,kTRUE);
  fithandler->FixParam(17,kTRUE);
  
  fithandler->FixParam(18,kTRUE); // Fix resolution weights
  fithandler->FixParam(19,kTRUE);
  
  //invariant mass function normalized between invMassMin, invMassMax  -->
  Double_t invMassMin = 2.0; // GeV
  Double_t invMassMax = 4.0; // GeV 
  Double_t massLow = (TDatabasePDG::Instance()->GetParticle(443)->Mass()) - invMassMin;
  Double_t massHigh = invMassMax - (TDatabasePDG::Instance()->GetParticle(443)->Mass());
  fithandler->SetMassWndLow(massLow);
  fithandler->SetMassWndHigh(massHigh);

  likely_obj = fithandler->LikelihoodPointer();
  likely_obj->SetAllParameters(paramInputValues); 
  likely_obj->ComputeMassIntegral();
  likely_obj->PrintStatus(); 
 
  aBtoJPSItoEle->DoMinimization();
  f.Close();
  return;
}


TH1F *psproperMCsecJPSI(){
   TH1F *hDecayTimeMCjpsifromB = new TH1F("hDecayTimeMCjpsifromB","B --> J/#Psi MC pseudo proper decay length",150,-3000,3000
); 
   hDecayTimeMCjpsifromB->SetBinContent(67,2);
   hDecayTimeMCjpsifromB->SetBinContent(68,2);
   hDecayTimeMCjpsifromB->SetBinContent(69,2);
   hDecayTimeMCjpsifromB->SetBinContent(70,4);
   hDecayTimeMCjpsifromB->SetBinContent(71,5);
   hDecayTimeMCjpsifromB->SetBinContent(72,17);
   hDecayTimeMCjpsifromB->SetBinContent(73,35);
   hDecayTimeMCjpsifromB->SetBinContent(74,96);
   hDecayTimeMCjpsifromB->SetBinContent(75,439);
   hDecayTimeMCjpsifromB->SetBinContent(76,11428);
   hDecayTimeMCjpsifromB->SetBinContent(77,9094);
   hDecayTimeMCjpsifromB->SetBinContent(78,7769);
   hDecayTimeMCjpsifromB->SetBinContent(79,6475);
   hDecayTimeMCjpsifromB->SetBinContent(80,5695);
   hDecayTimeMCjpsifromB->SetBinContent(81,4884);
   hDecayTimeMCjpsifromB->SetBinContent(82,4294);
   hDecayTimeMCjpsifromB->SetBinContent(83,3780);
   hDecayTimeMCjpsifromB->SetBinContent(84,3321);
   hDecayTimeMCjpsifromB->SetBinContent(85,2827);
   hDecayTimeMCjpsifromB->SetBinContent(86,2531);
   hDecayTimeMCjpsifromB->SetBinContent(87,2256);
   hDecayTimeMCjpsifromB->SetBinContent(88,2099);
   hDecayTimeMCjpsifromB->SetBinContent(89,1782);
   hDecayTimeMCjpsifromB->SetBinContent(90,1592);
   hDecayTimeMCjpsifromB->SetBinContent(91,1478);
   hDecayTimeMCjpsifromB->SetBinContent(92,1286);
   hDecayTimeMCjpsifromB->SetBinContent(93,1145);
   hDecayTimeMCjpsifromB->SetBinContent(94,980);
   hDecayTimeMCjpsifromB->SetBinContent(95,933);
   hDecayTimeMCjpsifromB->SetBinContent(96,865);
   hDecayTimeMCjpsifromB->SetBinContent(97,761);
   hDecayTimeMCjpsifromB->SetBinContent(98,654);
   hDecayTimeMCjpsifromB->SetBinContent(99,622);
   hDecayTimeMCjpsifromB->SetBinContent(100,587);
   hDecayTimeMCjpsifromB->SetBinContent(101,572);
   hDecayTimeMCjpsifromB->SetBinContent(102,449);
   hDecayTimeMCjpsifromB->SetBinContent(103,416);
   hDecayTimeMCjpsifromB->SetBinContent(104,346);
   hDecayTimeMCjpsifromB->SetBinContent(105,361);
   hDecayTimeMCjpsifromB->SetBinContent(105,361);
   hDecayTimeMCjpsifromB->SetBinContent(106,315);
   hDecayTimeMCjpsifromB->SetBinContent(107,265);
   hDecayTimeMCjpsifromB->SetBinContent(108,239);
   hDecayTimeMCjpsifromB->SetBinContent(109,247);
   hDecayTimeMCjpsifromB->SetBinContent(110,189);
   hDecayTimeMCjpsifromB->SetBinContent(111,223);
   hDecayTimeMCjpsifromB->SetBinContent(112,174);
   hDecayTimeMCjpsifromB->SetBinContent(113,184);
   hDecayTimeMCjpsifromB->SetBinContent(114,171);
   hDecayTimeMCjpsifromB->SetBinContent(115,153);
   hDecayTimeMCjpsifromB->SetBinContent(116,151);
   hDecayTimeMCjpsifromB->SetBinContent(117,126);
   hDecayTimeMCjpsifromB->SetBinContent(118,105);
   hDecayTimeMCjpsifromB->SetBinContent(119,98);
   hDecayTimeMCjpsifromB->SetBinContent(120,88);
   hDecayTimeMCjpsifromB->SetBinContent(121,92);
   hDecayTimeMCjpsifromB->SetBinContent(122,86);
   hDecayTimeMCjpsifromB->SetBinContent(123,73);
   hDecayTimeMCjpsifromB->SetBinContent(124,78);
   hDecayTimeMCjpsifromB->SetBinContent(125,63);
   hDecayTimeMCjpsifromB->SetBinContent(126,69);
   hDecayTimeMCjpsifromB->SetBinContent(127,55);
   hDecayTimeMCjpsifromB->SetBinContent(128,43);
   hDecayTimeMCjpsifromB->SetBinContent(129,46);
   hDecayTimeMCjpsifromB->SetBinContent(130,38);
   hDecayTimeMCjpsifromB->SetBinContent(131,34);
   hDecayTimeMCjpsifromB->SetBinContent(132,33);
   hDecayTimeMCjpsifromB->SetBinContent(133,51);
   hDecayTimeMCjpsifromB->SetBinContent(134,25);
   hDecayTimeMCjpsifromB->SetBinContent(135,44);
   hDecayTimeMCjpsifromB->SetBinContent(136,24);
   hDecayTimeMCjpsifromB->SetBinContent(137,25);
   hDecayTimeMCjpsifromB->SetBinContent(138,21);
   hDecayTimeMCjpsifromB->SetBinContent(139,23);
   hDecayTimeMCjpsifromB->SetBinContent(140,16);
   hDecayTimeMCjpsifromB->SetBinContent(141,24);
   hDecayTimeMCjpsifromB->SetBinContent(142,14);
   hDecayTimeMCjpsifromB->SetBinContent(143,21);
   hDecayTimeMCjpsifromB->SetBinContent(144,22);
   hDecayTimeMCjpsifromB->SetBinContent(145,22);
   hDecayTimeMCjpsifromB->SetBinContent(146,13);
   hDecayTimeMCjpsifromB->SetBinContent(147,10);
   hDecayTimeMCjpsifromB->SetBinContent(148,15);
   hDecayTimeMCjpsifromB->SetBinContent(149,12);
   hDecayTimeMCjpsifromB->SetBinContent(150,12);
   hDecayTimeMCjpsifromB->SetBinContent(151,231);
   return hDecayTimeMCjpsifromB;
}

void LoadLib(){
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3dielectron.so");
 }
