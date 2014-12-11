TH1F *psproperMCsecJPSI();
void LoadLib();
char * inputDistr = "data/result_data_70any.root";

Double_t *resParamFF = new Double_t[9];
Double_t *resParamFS = new Double_t[9];
Double_t *resParamSS = new Double_t[9];
Double_t *invMassParam = new Double_t[9];
Double_t *bkgParam = new Double_t[10];
Double_t Fb = 0.1;
Double_t Fsig = 0.291701; // 2.4 - 4

Bool_t range = 0; // 1 signal (2.92 - 3.16) 
                  // 0 all (2.4 - 4) 

//select types of candidates used for likelihood fit
TString resType = "FF;FS;SS";

Double_t weightType[] = {0.,0.,0.};

void FitCDFLocal();
AliDielectronBtoJPSItoEleCDFfitFCN *likely_obj = 0x0;

void  FitCDFLocal(){

  ///////////////////////////////////////////////////////////////////
  //
  // Example macro to read local N-Tuples of JPSI 
  // and bkg candidates and perform log-likelihood 
  // minimization using the minimization handler class
  //
  // Origin: C. Di Giglio
  //
  ///////////////////////////////////////////////////////////////////

  TH1F* hCsiMCPithya = new TH1F();
  TH1F* hCsiMCevtgen = new TH1F();
  // background parameters
  bkgParam[0] = 3.07667e-04;
  bkgParam[1] = 3.44581e-04;
  bkgParam[2] = 9.31570e-03;
  bkgParam[3] = 1.76852e+03;
  bkgParam[4] = 6.57758e+03;
  bkgParam[5] = 6.10483e+03;
  bkgParam[6] = 1.63737e+04;
  bkgParam[7] = 1.33739e-01;
  bkgParam[8] = 1.02785e+00;
  bkgParam[9] = 1.07043e+00;

  // resolution parameters FF
  resParamFF[0] = 1.30781e+03;
  resParamFF[1] = 0.;
  resParamFF[2] = 2.65464e+01;
  resParamFF[3] = 2.09853e+03;
  resParamFF[4] = 0.;
  resParamFF[5] = 5.75582e+01;
  resParamFF[6] = 1.77041e+02;
  resParamFF[7] = 3.49328e+00;
  resParamFF[8] = 3.74067e+02;

  // resolution parameters FS
  resParamFS[0] = 2.88183e+03;
  resParamFS[1] = 1.65460e+00;
  resParamFS[2] = 1.11174e+02;
  resParamFS[3] = 2.10971e+03;
  resParamFS[4] = 1.11906e+00;
  resParamFS[5] = 4.87024e+01;
  resParamFS[6] = 3.56243e+02;
  resParamFS[7] = 3.22939e+00;
  resParamFS[8] = 5.94102e+02;
 

  // resolution parameters SS
  resParamSS[0] = 5.77864e+03;
  resParamSS[1] = 1.38674e+00;
  resParamSS[2] = 2.80921e+02;
  resParamSS[3] = 1.17984e+04;
  resParamSS[4] = -5.67447e-01;
  resParamSS[5] = 1.04591e+02;
  resParamSS[6] = 1.14000e+03;
  resParamSS[7] = 1.25609e+00;
  resParamSS[8] = 1.43915e+03;

  // invariant mass parameters
  invMassParam[0] = 3.08547e+00;
  invMassParam[1] = 3.21099e-02;
  invMassParam[2] = 5.17599e-01;
  invMassParam[3] = 9.62091e+01;
  invMassParam[4] = 1.36800e+00;
  invMassParam[5] = 1.24918e+02;
  invMassParam[6] = 1.12386e+00;
  invMassParam[7] = 1.63711e+00;
  invMassParam[8] = -1.20962e+01; 

  hCsiMCPithya = psproperMCsecJPSI();  
  Double_t integral = 0;
  for(int i=1;i<hCsiMCPithya->GetNbinsX()+1; i++) integral += (hCsiMCPithya->GetBinContent(i)*hCsiMCPithya->GetBinWidth(i));
  hCsiMCPithya->Scale(1./integral);

  Double_t* x=0x0; Double_t* m=0x0; Int_t* type=0; Int_t n=0; 
  AliDielectronBtoJPSItoEle* aBtoJPSItoEle =new AliDielectronBtoJPSItoEle();
  Double_t paramInputValues[45] =  { bkgParam[3], // fWeightRes [0]
                                     bkgParam[4], // Fpos [1]
                                     bkgParam[5], // FNeg [2]
                                     bkgParam[6], // FSym [3]
                                     bkgParam[0], // LamdaPos [4]
                                     bkgParam[1], // LambdaNeg [5]
                                     bkgParam[2], // LamdaSym [6]
                                     Fb, // fB [7]
                                     Fsig, // FSig [8]
                                     invMassParam[0], // FMean [9]
                                     invMassParam[4], // fNexp [10]
                                     invMassParam[1], // fNsigma [11]
                                     invMassParam[2], // fAlpha [12]
                                     invMassParam[3],//fNorm [13]
                                     invMassParam[5], // fBkgNorm [14]
                                     invMassParam[6], // fBkgMean [15]
                                     invMassParam[7], // fBkgSlope [16]
                                     invMassParam[8], // fBkgConst [17]
                                     resParamFF[0], // Gaus1Norm [18]
                                     resParamFF[3], // Gaus2Norm [19]
                                     resParamFF[1], // Mean1Res[20]
                                     resParamFF[2], // Sigma1Res[21]
                                     resParamFF[4], // Mean2Res[22]
                                     resParamFF[5],  // Sigma2Res[23]
                                     resParamFF[6], // alfa res [24]
                                     resParamFF[7], // lambda res [25]
                                     resParamFF[8], // norm res [26]
                                     resParamFS[0], // Gaus1Norm [27]
                                     resParamFS[3], // Gaus2Norm [28]
                                     resParamFS[1], // Mean1Res[29]
                                     resParamFS[2], // Sigma1Res[30]
                                     resParamFS[4], // Mean2Res[31]
                                     resParamFS[5],  // Sigma2Res[32]
                                     resParamFS[6], // alfa res [33]
                                     resParamFS[7], // lambda res [34]
                                     resParamFS[8], // norm res [35]
                                     resParamSS[0], // Gaus1Norm [36]
                                     resParamSS[3], // Gaus2Norm [37]
                                     resParamSS[1], // Mean1Res[38]
                                     resParamSS[2], // Sigma1Res[39]
                                     resParamSS[4], // Mean2Res[40]
                                     resParamSS[5], // Sigma2Res[41]
                                     resParamSS[6], // alfa res [42]
                                     resParamSS[7], // lambda res [43]
                                     resParamSS[8]  // norm res [44]
                                     }; // Starting values for parameters 
 
  TFile f(inputDistr);
  TList * list = (TList*)f.Get("resultAOD");
  TNtuple *nt;
  if(range == 0) nt=(TNtuple*)list->FindObject("fNtupleJPSItype");
  if(range == 1) nt=(TNtuple*)list->FindObject("fNtupleJPSItype_signal");
  nt->ls();

  aBtoJPSItoEle->SetResTypeAnalysis(resType);
  aBtoJPSItoEle->ReadCandidates(nt,x,m,type,n); // read N-Tuples
  printf("+++\n+++ Number of total candidates (prim J/psi + secondary J/psi + bkg) ---> %d candidates \n+++\n",n);

  aBtoJPSItoEle->SetFitHandler(x,m,type,n); // Set the fit handler with given values of x, m, # of candidates 

  aBtoJPSItoEle->CloneMCtemplate(hCsiMCPithya);    // clone MC template and copy internally
                                                   // in this way any model can be setted from outside
  //aBtoJPSItoEle->CloneMCtemplate(hCsiMCEvtGen);  // clone MC template and copy internally
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
  fithandler->SetCrystalBallFunction(kTRUE);
  
  fithandler->FixParam(0,kTRUE); // bkg weights
  fithandler->FixParam(1,kTRUE); 
  fithandler->FixParam(2,kTRUE);
  fithandler->FixParam(3,kTRUE);
  fithandler->FixParam(4,kTRUE); // lPos
  fithandler->FixParam(5,kTRUE); // lNeg  Fix bkg exponential   
  fithandler->FixParam(6,kTRUE); // lSym
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
  fithandler->FixParam(18,kTRUE);
  fithandler->FixParam(19,kTRUE); // norm2
  fithandler->FixParam(20,kTRUE);
  fithandler->FixParam(21,kTRUE);
  fithandler->FixParam(22,kTRUE); // mean2
  fithandler->FixParam(23,kTRUE); // sigma2
  fithandler->FixParam(24,kTRUE); // alfaRes
  fithandler->FixParam(25,kTRUE); // lambdaRes
  fithandler->FixParam(26,kTRUE); // ConstRes
  // resolution func for FS and SS types
  fithandler->FixParam(27,kTRUE);
  fithandler->FixParam(28,kTRUE); // norm2
  fithandler->FixParam(29,kTRUE);
  fithandler->FixParam(30,kTRUE);
  fithandler->FixParam(31,kTRUE); // mean2
  fithandler->FixParam(32,kTRUE); // sigma2
  fithandler->FixParam(33,kTRUE); // alfaRes
  fithandler->FixParam(34,kTRUE); // lambdaRes
  fithandler->FixParam(35,kTRUE); // ConstRes
  fithandler->FixParam(36,kTRUE);
  fithandler->FixParam(37,kTRUE); // norm2
  fithandler->FixParam(38,kTRUE);
  fithandler->FixParam(39,kTRUE);
  fithandler->FixParam(40,kTRUE); // mean2
  fithandler->FixParam(41,kTRUE); // sigma2
  fithandler->FixParam(42,kTRUE); // alfaRes
  fithandler->FixParam(43,kTRUE); // lambdaRes
  fithandler->FixParam(44,kTRUE); // ConstRes

  //invMass funcions normalized between bandLow - bandUp -->
  Double_t bandLow;
  Double_t bandUp;

  if(range==1){
   bandLow = 2.92;
   bandUp = 3.16;
  }
  
  if(range==0){
  bandLow = 2.4;
  bandUp = 4.;
  }
  
  Double_t massLow = (TDatabasePDG::Instance()->GetParticle(443)->Mass()) - bandLow;
  Double_t massHigh = bandUp - (TDatabasePDG::Instance()->GetParticle(443)->Mass());
  fithandler->SetMassWndLow(massLow);
  fithandler->SetMassWndHigh(massHigh);

  likely_obj = fithandler->LikelihoodPointer();
  likely_obj->SetAllParameters(paramInputValues); 
  likely_obj->ComputeMassIntegral();
  likely_obj->PrintStatus(); 

  aBtoJPSItoEle->DoMinimization();
 
  Double_t Fsig_fromFit = fithandler->GetParameter(8);
  Double_t Fb_fromFit = fithandler->GetParameter(7);

  Double_t Fsig_err = fithandler->GetParameterError(8);
  Double_t Fb_err = fithandler->GetParameterError(7);  

  //fill histos from Ntupla
  TH1F *fInvMass = new TH1F("fInvMass","Invariant Mass; InvMass[GeV]; Entries/40MeV",40,2.4,2.4+40*.04); // step 40MeV
  TH1F *fpsproperSignal = new TH1F("psproper_decay_length",Form("psproper_decay_length_distrib(%f < M < %f);X [#mu m];Entries/40#mu m",bandLow,bandUp),150,-3000.,3000.); 
  
  Float_t mass , psproper, typeCand;
  TString arrType[]={"SS","FS","FF"}; 
  nt->SetBranchAddress("Xdecaytime",&psproper);
  nt->SetBranchAddress("Mass",&mass);
  nt->SetBranchAddress("Type",&typeCand);
  Int_t fNcurrent=0; Double_t nCandSel = 0;
  Int_t nb = (Int_t)nt->GetEvent(fNcurrent);

  for (Int_t iev=0; iev<(nt->GetEntries()); iev++){
   if(resType.Contains(arrType[(Int_t)typeCand])){
   nCandSel += 1;
   weightType[(Int_t)typeCand] += 1.;
   fInvMass->Fill(mass);
   fpsproperSignal->Fill(psproper);}
   fNcurrent++;
   nb = (Int_t) nt->GetEvent(fNcurrent);
   }
   likely_obj->SetWeightType(weightType[2]/nCandSel,weightType[1]/nCandSel,weightType[0]/nCandSel);
  
  // draw psproper total
  TCanvas *d6 = new TCanvas();
  d6->SetLogy();
  Double_t norm1 = ((Double_t)fpsproperSignal->GetEntries())*fpsproperSignal->GetBinWidth(1);
  TF1 *psproperTot = likely_obj->GetEvaluateCDFDecayTimeTotalDistrAllTypes(-1.e+04, 1.e+04,norm1);
  TFitResultPtr rPsproper = fpsproperSignal->Fit(psproperTot->GetName(),"S");
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.035);
  fpsproperSignal->DrawCopy("E");
  lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.3f ",(rPsproper->Chi2()/(Double_t)rPsproper->Ndf())));
  if(range == 0) lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, Fsig_fromFit, Fsig_err));
  else lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f](fixed from LS) = %4.3f", bandLow, bandUp, Fsig_fromFit));
  lat->DrawLatex(0.53, 0.62, Form("F_{B} = %4.3f #pm %4.3f", Fb_fromFit, Fb_err)); 

  //prompt jpsi
  Double_t normPrompt = (psproperTot->GetParameter(0))*Fsig_fromFit*(1 - Fb_fromFit);
  TF1 *prompt = likely_obj->GetResolutionFuncAllTypes(-1.e+04,1.e+04,normPrompt);
  prompt->SetLineColor(2);
  prompt->Draw("same");
 
  //secondary Jpsi 
  Double_t normSec = (psproperTot->GetParameter(0))*Fsig_fromFit*Fb_fromFit;
  TF1 *templateMC = likely_obj->GetFunBAllTypes(-1.e+04,1.e+04,normSec);
  templateMC->SetLineColor(6);
  templateMC->SetFillColor(6);
  templateMC->SetFillStyle(3005);
  templateMC->Draw("same");
 
  //bkg
  Double_t normBkg = (psproperTot->GetParameter(0))*(1 - Fsig_fromFit);
  TF1 *psProperBack = likely_obj->GetEvaluateCDFDecayTimeBkgDistrAllTypes(-1.e+04,1.e+04,normBkg);
  psProperBack->SetLineColor(3);
  psProperBack->Draw("same");

  fpsproperSignal->Sumw2();
  fpsproperSignal->DrawCopy("same"); 
  
  //legend
  TLegend *leg=new TLegend(0.17,0.72,0.42,0.88);
  leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextFont(42);
  leg->SetFillStyle(0); leg->SetMargin(0.25); //separation symbol-text
  leg->SetEntrySeparation(0.15);
  leg->AddEntry(psproperTot, "all","l");
  leg->AddEntry(prompt, "prompt J/#psi","l");
  leg->AddEntry(templateMC, "secondary J/#psi","l");
  leg->AddEntry(psProperBack, "bkg","l");
  leg->Draw("same");

  // draw invariant mass 
  Double_t norm =((Double_t)fInvMass->GetEntries())*fInvMass->GetBinWidth(1);
  TF1 *invMassFunc = likely_obj->GetEvaluateCDFInvMassTotalDistr(bandLow,bandUp,norm);
  TCanvas *d5 = new TCanvas();
  Double_t intTot = invMassFunc->Integral(bandLow,bandUp); 
  printf("intTot (%f-%f)= %f \n",bandLow,bandUp,intTot);
  TFitResultPtr rMass = fInvMass->Fit(invMassFunc->GetName(),"S");
  fInvMass->DrawCopy("E"); 
  lat->DrawLatex(0.53, 0.82, Form("#chi^{2}/dof = %4.2f ",(rMass->Chi2()/(Double_t)rMass->Ndf())));
  if(range == 0) lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f] = %4.3f #pm %4.3f", bandLow, bandUp, Fsig_fromFit, Fsig_err));
  else lat->DrawLatex(0.53, 0.72, Form("F_{Sig}[%1.2f - %1.2f](fixed from LS) = %4.3f", bandLow, bandUp, Fsig_fromFit));

  TF1 *invMassSig = likely_obj->GetEvaluateCDFInvMassSigDistr(2.4,4,norm*Fsig_fromFit);
  TF1 *invMassBkg = likely_obj->GetEvaluateCDFInvMassBkgDistr(2.4,4,norm*(1.-Fsig_fromFit));

  invMassSig->SetLineColor(2);
  invMassBkg->SetLineColor(3);

  invMassSig->Draw("same");
  invMassBkg->Draw("same");
  
  if(range == 0){ 
  Double_t integSig = invMassSig->Integral(2.92,3.16);
  Double_t integBkg = invMassBkg->Integral(2.92,3.16);
  Double_t Fsig_new = integSig/(integSig+integBkg); 
  printf(" Fsig(rescaled) = %f \n",Fsig_new);
  }
  
  f.Close();
  return;
}

void LoadLib(){
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG3dielectron");
 }


TH1F *psproperMCsecJPSI(){
   TH1F *hDecayTimeMCjpsifromB = new TH1F("hDecayTimeMCjpsifromB","B --> J/#Psi MC pseudo proper decay length",150,-3000,3000);
   hDecayTimeMCjpsifromB->SetBinContent(67,1);
   hDecayTimeMCjpsifromB->SetBinContent(68,2);
   hDecayTimeMCjpsifromB->SetBinContent(69,1);
   hDecayTimeMCjpsifromB->SetBinContent(70,1);
   hDecayTimeMCjpsifromB->SetBinContent(71,1);
   hDecayTimeMCjpsifromB->SetBinContent(72,9);
   hDecayTimeMCjpsifromB->SetBinContent(73,16);
   hDecayTimeMCjpsifromB->SetBinContent(74,51);
   hDecayTimeMCjpsifromB->SetBinContent(75,204);
   hDecayTimeMCjpsifromB->SetBinContent(76,6494);
   hDecayTimeMCjpsifromB->SetBinContent(77,5249);
   hDecayTimeMCjpsifromB->SetBinContent(78,4530);
   hDecayTimeMCjpsifromB->SetBinContent(79,3786);
   hDecayTimeMCjpsifromB->SetBinContent(80,3376);
   hDecayTimeMCjpsifromB->SetBinContent(81,2892);
   hDecayTimeMCjpsifromB->SetBinContent(82,2585);
   hDecayTimeMCjpsifromB->SetBinContent(83,2261);
   hDecayTimeMCjpsifromB->SetBinContent(84,1997);
   hDecayTimeMCjpsifromB->SetBinContent(85,1712);
   hDecayTimeMCjpsifromB->SetBinContent(86,1521);
   hDecayTimeMCjpsifromB->SetBinContent(87,1318);
   hDecayTimeMCjpsifromB->SetBinContent(88,1293);
   hDecayTimeMCjpsifromB->SetBinContent(89,1085);
   hDecayTimeMCjpsifromB->SetBinContent(90,936);
   hDecayTimeMCjpsifromB->SetBinContent(91,872);
   hDecayTimeMCjpsifromB->SetBinContent(92,768);
   hDecayTimeMCjpsifromB->SetBinContent(93,686);
   hDecayTimeMCjpsifromB->SetBinContent(94,594);
   hDecayTimeMCjpsifromB->SetBinContent(95,545);
   hDecayTimeMCjpsifromB->SetBinContent(96,545);
   hDecayTimeMCjpsifromB->SetBinContent(97,462);
   hDecayTimeMCjpsifromB->SetBinContent(98,405);
   hDecayTimeMCjpsifromB->SetBinContent(99,376);
   hDecayTimeMCjpsifromB->SetBinContent(100,367);
   hDecayTimeMCjpsifromB->SetBinContent(101,355);
   hDecayTimeMCjpsifromB->SetBinContent(102,252);
   hDecayTimeMCjpsifromB->SetBinContent(103,250);
   hDecayTimeMCjpsifromB->SetBinContent(104,207);
   hDecayTimeMCjpsifromB->SetBinContent(105,204);
   hDecayTimeMCjpsifromB->SetBinContent(106,182);
   hDecayTimeMCjpsifromB->SetBinContent(107,168);
   hDecayTimeMCjpsifromB->SetBinContent(108,125);
   hDecayTimeMCjpsifromB->SetBinContent(109,142);
   hDecayTimeMCjpsifromB->SetBinContent(110,116);
   hDecayTimeMCjpsifromB->SetBinContent(111,132);
   hDecayTimeMCjpsifromB->SetBinContent(112,100);
   hDecayTimeMCjpsifromB->SetBinContent(113,115);
   hDecayTimeMCjpsifromB->SetBinContent(114,93);
   hDecayTimeMCjpsifromB->SetBinContent(115,85);
   hDecayTimeMCjpsifromB->SetBinContent(116,96);
   hDecayTimeMCjpsifromB->SetBinContent(117,73);
   hDecayTimeMCjpsifromB->SetBinContent(118,76);
   hDecayTimeMCjpsifromB->SetBinContent(119,56);
   hDecayTimeMCjpsifromB->SetBinContent(120,53);
   hDecayTimeMCjpsifromB->SetBinContent(121,46);
   hDecayTimeMCjpsifromB->SetBinContent(122,50);
   hDecayTimeMCjpsifromB->SetBinContent(123,42);
   hDecayTimeMCjpsifromB->SetBinContent(124,40);
   hDecayTimeMCjpsifromB->SetBinContent(125,36);
   hDecayTimeMCjpsifromB->SetBinContent(126,43);
   hDecayTimeMCjpsifromB->SetBinContent(127,31);
   hDecayTimeMCjpsifromB->SetBinContent(128,20);
   hDecayTimeMCjpsifromB->SetBinContent(129,18);
   hDecayTimeMCjpsifromB->SetBinContent(130,23);
   hDecayTimeMCjpsifromB->SetBinContent(131,19);
   hDecayTimeMCjpsifromB->SetBinContent(132,18);
   hDecayTimeMCjpsifromB->SetBinContent(133,29);
   hDecayTimeMCjpsifromB->SetBinContent(134,20);
   hDecayTimeMCjpsifromB->SetBinContent(135,26);
   hDecayTimeMCjpsifromB->SetBinContent(136,13);
   hDecayTimeMCjpsifromB->SetBinContent(137,12);
   hDecayTimeMCjpsifromB->SetBinContent(138,12);
   hDecayTimeMCjpsifromB->SetBinContent(139,14);
   hDecayTimeMCjpsifromB->SetBinContent(140,8);
   hDecayTimeMCjpsifromB->SetBinContent(141,14);
   hDecayTimeMCjpsifromB->SetBinContent(142,9);
   hDecayTimeMCjpsifromB->SetBinContent(143,14);
   hDecayTimeMCjpsifromB->SetBinContent(144,13);
   hDecayTimeMCjpsifromB->SetBinContent(145,13);
   hDecayTimeMCjpsifromB->SetBinContent(146,11);
   hDecayTimeMCjpsifromB->SetBinContent(147,5);
   hDecayTimeMCjpsifromB->SetBinContent(148,11);
   hDecayTimeMCjpsifromB->SetBinContent(149,7);
   hDecayTimeMCjpsifromB->SetBinContent(150,6);
   hDecayTimeMCjpsifromB->SetBinContent(151,136);
   hDecayTimeMCjpsifromB->SetBinError(67,1);
   hDecayTimeMCjpsifromB->SetBinError(68,1.414214);
   hDecayTimeMCjpsifromB->SetBinError(69,1);
   hDecayTimeMCjpsifromB->SetBinError(70,1);
   hDecayTimeMCjpsifromB->SetBinError(71,1);
   hDecayTimeMCjpsifromB->SetBinError(72,3);
   hDecayTimeMCjpsifromB->SetBinError(73,4);
   hDecayTimeMCjpsifromB->SetBinError(74,7.141428);
   hDecayTimeMCjpsifromB->SetBinError(75,14.28286);
   hDecayTimeMCjpsifromB->SetBinError(76,80.58536);
   hDecayTimeMCjpsifromB->SetBinError(77,72.44998);
   hDecayTimeMCjpsifromB->SetBinError(78,67.30527);
   hDecayTimeMCjpsifromB->SetBinError(79,61.53048);
   hDecayTimeMCjpsifromB->SetBinError(80,58.10336);
   hDecayTimeMCjpsifromB->SetBinError(81,53.77732);
   hDecayTimeMCjpsifromB->SetBinError(82,50.8429);
   hDecayTimeMCjpsifromB->SetBinError(83,47.54997);
   hDecayTimeMCjpsifromB->SetBinError(84,44.68781);
   hDecayTimeMCjpsifromB->SetBinError(85,41.37632);
   hDecayTimeMCjpsifromB->SetBinError(86,39);
   hDecayTimeMCjpsifromB->SetBinError(87,36.30427);
   hDecayTimeMCjpsifromB->SetBinError(88,35.95831);
   hDecayTimeMCjpsifromB->SetBinError(89,32.93934);
   hDecayTimeMCjpsifromB->SetBinError(90,30.59412);
   hDecayTimeMCjpsifromB->SetBinError(91,29.52965);
   hDecayTimeMCjpsifromB->SetBinError(92,27.71281);
   hDecayTimeMCjpsifromB->SetBinError(93,26.1916);
   hDecayTimeMCjpsifromB->SetBinError(94,24.37212);
   hDecayTimeMCjpsifromB->SetBinError(95,23.34524);
   hDecayTimeMCjpsifromB->SetBinError(96,23.34524);
   hDecayTimeMCjpsifromB->SetBinError(97,21.49419);
   hDecayTimeMCjpsifromB->SetBinError(98,20.12461);
   hDecayTimeMCjpsifromB->SetBinError(99,19.39072);
   hDecayTimeMCjpsifromB->SetBinError(100,19.15724);
   hDecayTimeMCjpsifromB->SetBinError(101,18.84144);
   hDecayTimeMCjpsifromB->SetBinError(102,15.87451);
   hDecayTimeMCjpsifromB->SetBinError(103,15.81139);
   hDecayTimeMCjpsifromB->SetBinError(104,14.38749);
   hDecayTimeMCjpsifromB->SetBinError(105,14.28286);
   hDecayTimeMCjpsifromB->SetBinError(106,13.49074);
   hDecayTimeMCjpsifromB->SetBinError(107,12.96148);
   hDecayTimeMCjpsifromB->SetBinError(108,11.18034);
   hDecayTimeMCjpsifromB->SetBinError(109,11.91638);
   hDecayTimeMCjpsifromB->SetBinError(110,10.77033);
   hDecayTimeMCjpsifromB->SetBinError(111,11.48913);
   hDecayTimeMCjpsifromB->SetBinError(112,10);
   hDecayTimeMCjpsifromB->SetBinError(113,10.72381);
   hDecayTimeMCjpsifromB->SetBinError(114,9.643651);
   hDecayTimeMCjpsifromB->SetBinError(115,9.219544);
   hDecayTimeMCjpsifromB->SetBinError(116,9.797959);
   hDecayTimeMCjpsifromB->SetBinError(117,8.544004);
   hDecayTimeMCjpsifromB->SetBinError(118,8.717798);
   hDecayTimeMCjpsifromB->SetBinError(119,7.483315);
   hDecayTimeMCjpsifromB->SetBinError(120,7.28011);
   hDecayTimeMCjpsifromB->SetBinError(121,6.78233);
   hDecayTimeMCjpsifromB->SetBinError(122,7.071068);
   hDecayTimeMCjpsifromB->SetBinError(123,6.480741);
   hDecayTimeMCjpsifromB->SetBinError(124,6.324555);
   hDecayTimeMCjpsifromB->SetBinError(125,6);
   hDecayTimeMCjpsifromB->SetBinError(126,6.557439);
   hDecayTimeMCjpsifromB->SetBinError(127,5.567764);
   hDecayTimeMCjpsifromB->SetBinError(128,4.472136);
   hDecayTimeMCjpsifromB->SetBinError(129,4.242641);
   hDecayTimeMCjpsifromB->SetBinError(130,4.795832);
   hDecayTimeMCjpsifromB->SetBinError(131,4.358899);
   hDecayTimeMCjpsifromB->SetBinError(132,4.242641);
   hDecayTimeMCjpsifromB->SetBinError(133,5.385165);
   hDecayTimeMCjpsifromB->SetBinError(134,4.472136);
   hDecayTimeMCjpsifromB->SetBinError(135,5.09902);
   hDecayTimeMCjpsifromB->SetBinError(136,3.605551);
   hDecayTimeMCjpsifromB->SetBinError(137,3.464102);
   hDecayTimeMCjpsifromB->SetBinError(138,3.464102);
   hDecayTimeMCjpsifromB->SetBinError(139,3.741657);
   hDecayTimeMCjpsifromB->SetBinError(140,2.828427);
   hDecayTimeMCjpsifromB->SetBinError(141,3.741657);
   hDecayTimeMCjpsifromB->SetBinError(142,3);
   hDecayTimeMCjpsifromB->SetBinError(143,3.741657);
   hDecayTimeMCjpsifromB->SetBinError(144,3.605551);
   hDecayTimeMCjpsifromB->SetBinError(145,3.605551);
   hDecayTimeMCjpsifromB->SetBinError(146,3.316625);
   hDecayTimeMCjpsifromB->SetBinError(147,2.236068);
   hDecayTimeMCjpsifromB->SetBinError(148,3.316625);
   hDecayTimeMCjpsifromB->SetBinError(149,2.645751);
   hDecayTimeMCjpsifromB->SetBinError(150,2.44949);
   hDecayTimeMCjpsifromB->SetBinError(151,11.6619);
   return hDecayTimeMCjpsifromB;
}


