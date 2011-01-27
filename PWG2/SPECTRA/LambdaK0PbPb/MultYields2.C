#include "FitControl.h"

void MultYields2(TH3F *PtMassMult, Int_t particleMode, Int_t MultBin, Char_t* label){

  TString *tLabel = new TString(label);
    // Do .L on command line first
    //gROOT->LoadMacro("macros/PtMassAna2.C");
/* Modifications to produce a single uncorrected spectrum in a particular mult. bin for K0 or Lambda
Old ratio code preserved in MultYieldsRatio.C
Dec 2003 */
  Char_t* part; //for name of particle
  Float_t BR; //branching ratio - check them
  if (particleMode==0){
    part = "K0";
      BR=0.686;
  } else if (particleMode==1){
    part = "Lambda";
    BR=0.639;
  } else if (particleMode==2){
      part = "AntiLambda";
      BR=0.639;
  } else if (particleMode==3){
    part = "Lambda+antiLambda";
    BR=0.639;
  } else if (particleMode==4){
    part = "Xi";
    BR=1.0; // Should be Lam ->p pi Br
  } else if (particleMode ==6){
    part = "Omega";
    BR=1.0; // Should be Lam->p pi  * Om -> K Lam
  }

  
//   //Setup. Multiplicity bin boundaries, N_binary scaling numbers
//   //These are latest numbers from M. Miller Aug03
//   Float_t NBin[4] = {7.5,3.99,10.2,15.04};
//   Float_t Veff[4] = {0.925,0.875,0.999,1.0};
// // These are dAuFilt numbers
// //  Float_t Nev1=5.03505e+06, Nev2=2.4501e+06, Nev3=2.23543e+06, NevMB=9.722229e+06;    
//   Float_t MultLo[4]={0,0,10,17};
//   Float_t MultHi[4]={120,9,16,120};
// //These are dAuFilt2 numbers (discrepancy with above, changed GoodEvent BUT that shouldn't decrease Nevent)
// //New N events figures from the 25cm cut file
//  if (tLabel->Contains("25cm")) {
//    Float_t Nev[4]={4.586559e+06,2.37952e+06,1.11382e+06,1.09322e+06};
//  }
//  else {
//    Float_t Nev[4]={9.590239e+06,5.0683e+06,2.31685e+06,2.20508e+06};
//  }

// Redo the possible multiplicity bin stuff for Cu+Cu SVT test
// As a start just use one bin
  Float_t NBin[1] = {63.5}; // Not really that useful for this study
                            // used 0-80% bin from Johan's page
  Float_t Veff[1] = {0.91}; // TBD from data?
  Float_t Nev[1] = {34712502.0}; // For MB
//  Float_t Nev[1] = {2727356.0}; // For HM
  Float_t MultLo[1] = {0};
  Float_t MultHi[1] = {300};
 
TString title[1]={"MinimumBias"};
  //Make 2D projections from the 3D histogram
  //Minbias (i.e. everything)
  PtMassMult->GetZaxis()->SetRange(MultLo[MultBin],MultHi[MultBin]);
  TH2F* hPtMass = (TH2F*)PtMassMult->Project3D("XY");// FIX:MF
  hPtMass->SetTitle("PtMass");

  hPtMass->Draw();

  TObjArray *controllerArray = new TObjArray(40,1); // 2nd arg means can count from 1!

  //Here probably need switch-case, depending on mult bin and particle
  // LoPt, HiPt, polynomial order, rebinning factor

  /////  LAMBDA and LAMBDA+ANTILAMBDA (Combination)

  if(particleMode == 1 || particleMode == 3){ // Lambda or Lambda+Anti-Lambda
  //  controllerArray->AddLast(new FitControl(0.1,0.2, 1,2));
  //  controllerArray->AddLast(new FitControl(0.2,0.3, 1,2));
  //  controllerArray->AddLast(new FitControl(0.3,0.4, 1,1, 1.096,1.16))
  //  controllerArray->AddLast(new FitControl(0.4,0.5, 2,1, 1.095,1.15));
    if(tLabel->Contains("SVT")){
      controllerArray->AddLast(new FitControl(0.4,0.5, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.5,0.6, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.6,0.7, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.7,0.8, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.8,1.0, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(1.0,1.2, 2,1));
      controllerArray->AddLast(new FitControl(1.2,1.4, 2,1));
      controllerArray->AddLast(new FitControl(1.4,1.6, 2,1));
      controllerArray->AddLast(new FitControl(1.6,1.8, 2,1));
      controllerArray->AddLast(new FitControl(1.8,2.0, 2,2));
      controllerArray->AddLast(new FitControl(2.0,2.2, 2,2, 1.095,1.17));
      controllerArray->AddLast(new FitControl(2.2,2.4, 2,2, 1.095,1.16));
//       controllerArray->AddLast(new FitControl(2.4,2.6, 2,1, 1.095,1.17));
//       controllerArray->AddLast(new FitControl(2.6,2.8, 1,2, 1.095,1.17));
//       controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 1.095,1.17));
//       controllerArray->AddLast(new FitControl(3.0,3.5, 1,2, 1.098,1.16));
//       controllerArray->AddLast(new FitControl(3.5,4.0, 1,2, 1.095,1.17));
    }else{
      controllerArray->AddLast(new FitControl(0.5,0.6, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.6,0.7, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.7,0.8, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(0.8,1.0, 2,1, 1.095,1.15));
      controllerArray->AddLast(new FitControl(1.0,1.2, 2,1));
      controllerArray->AddLast(new FitControl(1.2,1.4, 2,1));
      controllerArray->AddLast(new FitControl(1.4,1.6, 2,1));
      controllerArray->AddLast(new FitControl(1.6,1.8, 2,1));
      controllerArray->AddLast(new FitControl(1.8,2.0, 2,1));
      controllerArray->AddLast(new FitControl(2.0,2.2, 2,1, 1.095,1.17));
      controllerArray->AddLast(new FitControl(2.2,2.4, 2,1, 1.095,1.16));
      controllerArray->AddLast(new FitControl(2.4,2.6, 2,1, 1.095,1.17));
      controllerArray->AddLast(new FitControl(2.6,2.8, 1,2, 1.095,1.17));
      controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 1.095,1.17));
      controllerArray->AddLast(new FitControl(3.0,3.5, 1,2, 1.098,1.16));
      controllerArray->AddLast(new FitControl(3.5,4.0, 1,2, 1.095,1.17));
    }
//   if (MultBin==2){controllerArray->AddLast(new FitControl(4.0,4.5, 1,2, 1.090,1.17));}
//   else if (MultBin==3){controllerArray->AddLast(new FitControl(4.0,4.5, 0,3, 1.090,1.17));}
//   else {controllerArray->AddLast(new FitControl(4.0,4.5, 1,2, 1.095,1.17));}
//   controllerArray->AddLast(new FitControl(4.5,5.0, 0,4, 1.095,1.17));
//   if (MultBin==1){  controllerArray->AddLast(new FitControl(5.0,5.5, 0,5, 1.095,1.17));}
//   else if (MultBin==2){controllerArray->AddLast(new FitControl(5.0,5.5, 0,6, 1.085,1.16));} 
//   else {controllerArray->AddLast(new FitControl(5.0,5.5, 0,4, 1.095,1.17));}
//   if(MultBin==3){
//     controllerArray->AddLast(new FitControl(5.5,6.0, 0,5, 1.085,1.165));
//   } else if(MultBin==2){
//     controllerArray->AddLast(new FitControl(5.5,6.0, 0,5, 1.085,1.17));
//   } else{
//     controllerArray->AddLast(new FitControl(5.5,6.0, 0,6, 1.085,1.17));
//   }  
  //  controllerArray->AddLast(new FitControl(6.0,6.5, 0,6));
  //  controllerArray->AddLast(new FitControl(6.5,7.0, 0,10));
  // The above, bins from 500 MeV to 6 GeV works for CENTRAL Lambdas
  // i.e all chisq/dof for backgrounds ~<2, no wild fits, perhaps sig/bkgrd too
  // high for 2.8-3.0 GeV bin
  }
  /// ANTI LAMBDA ---->
  else if (particleMode == 2){ // Anti-Lambdas
  controllerArray->AddLast(new FitControl(0.3,0.4, 2,1, 1.095,1.15));
  if(MultBin==3) { controllerArray->AddLast(new FitControl(0.4,0.5, 2,1, 1.095,1.16));}
  else {  controllerArray->AddLast(new FitControl(0.4,0.5, 2,1, 1.095,1.15));}
  controllerArray->AddLast(new FitControl(0.5,0.6, 2,1, 1.095,1.15));
  controllerArray->AddLast(new FitControl(0.6,0.7, 2,1, 1.095,1.15));
  controllerArray->AddLast(new FitControl(0.7,0.8, 2,1, 1.090,1.15));
  controllerArray->AddLast(new FitControl(0.8,1.0, 2,1, 1.095,1.15));
  controllerArray->AddLast(new FitControl(1.0,1.2, 2,1, 1.090,1.17));
  
  if (MultBin==1){  controllerArray->AddLast(new FitControl(1.2,1.4, 2,1, 1.090,1.16));}
  else{controllerArray->AddLast(new FitControl(1.2,1.4, 2,1, 1.090,1.17));}
  controllerArray->AddLast(new FitControl(1.4,1.6, 2,1, 1.090,1.17));
  controllerArray->AddLast(new FitControl(1.6,1.8, 2,1, 1.090,1.17));
  controllerArray->AddLast(new FitControl(1.8,2.0, 2,1, 1.090,1.17));
  controllerArray->AddLast(new FitControl(2.0,2.2, 2,1, 1.095,1.17));
  controllerArray->AddLast(new FitControl(2.2,2.4, 2,1, 1.095,1.16));
  if(MultBin==1){  controllerArray->AddLast(new FitControl(2.4,2.6, 1,1, 1.095,1.17));}
  else if(MultBin==3){  controllerArray->AddLast(new FitControl(2.4,2.6, 1,1, 1.095,1.17));}
  else {  controllerArray->AddLast(new FitControl(2.4,2.6, 2,1, 1.095,1.17));}
  controllerArray->AddLast(new FitControl(2.6,2.8, 1,2, 1.095,1.17));
  if(MultBin==2){  controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 1.090,1.17));}
  else if(MultBin==3){  controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 1.090,1.17));}
  else{  controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 1.095,1.17));}
  controllerArray->AddLast(new FitControl(3.0,3.5, 1,2, 1.095,1.17));
  controllerArray->AddLast(new FitControl(3.5,4.0, 1,2, 1.095,1.17));
  controllerArray->AddLast(new FitControl(4.0,4.5, 1,2, 1.090,1.17));
  if (MultBin==1){
    //  controllerArray->AddLast(new FitControl(4.5,5.0, 0,5, 1.095,1.17));
   controllerArray->AddLast(new FitControl(4.5,5.0, 0,3, 1.085,1.17));
  controllerArray->AddLast(new FitControl(5.0,5.5, 0,3, 1.090,1.17));
  //controllerArray->AddLast(new FitControl(5.5,6.0, 0,6, 1.085,1.17));
  } else if (MultBin==2){
  controllerArray->AddLast(new FitControl(4.5,5.0, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.0,5.5, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.5,6.0, 0,4, 1.085,1.17));
  }  else if (MultBin==3){
  controllerArray->AddLast(new FitControl(4.5,5.0, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.0,5.5, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.5,6.0, 0,4, 1.085,1.18));
  }  else{
  controllerArray->AddLast(new FitControl(4.5,5.0, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.0,5.5, 0,4, 1.095,1.17));
  controllerArray->AddLast(new FitControl(5.5,6.0, 0,6, 1.085,1.17));
  }
  } // end if anti-Lambda
  else if (particleMode == 0){ // K0s case
    if(tLabel->Contains("SVT")){
      controllerArray->AddLast(new FitControl(0.4,0.5, 1,1, 0.441,0.558));
      controllerArray->AddLast(new FitControl(0.5,0.6, 1,1, 0.442,0.559));
      controllerArray->AddLast(new FitControl(0.6,0.7, 1,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(0.7,0.8, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(0.8,1.0, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.0,1.2, 2,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(1.2,1.4, 2,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.4,1.6, 2,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(1.6,1.8, 2,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.8,2.0, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.0,2.2, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.2,2.4, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.4,2.6, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.6,2.8, 1,4, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.8,3.0, 1,4, 0.44,0.56));
    }else{
      //controllerArray->AddLast(new FitControl(0.1,0.2, 1,2, 0.44,0.56));
      //controllerArray->AddLast(new FitControl(0.2,0.3, 1,2, 0.44,0.56));
      //controllerArray->AddLast(new FitControl(0.3,0.4, 1,1, 0.442,0.559));
      controllerArray->AddLast(new FitControl(0.4,0.5, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(0.5,0.6, 1,1, 0.442,0.559));
      controllerArray->AddLast(new FitControl(0.6,0.7, 1,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(0.7,0.8, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(0.8,1.0, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.0,1.2, 2,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(1.2,1.4, 2,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.4,1.6, 2,1, 0.442,0.558));
      controllerArray->AddLast(new FitControl(1.6,1.8, 2,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(1.8,2.0, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.0,2.2, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.2,2.4, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.4,2.6, 1,1, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.6,2.8, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(2.8,3.0, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(3.0,3.5, 1,2, 0.44,0.56));
      controllerArray->AddLast(new FitControl(3.5,4.0, 1,2, 0.44,0.56));
      //   controllerArray->AddLast(new FitControl(4.0,4.5, 1,2, 0.418,0.578));
      //   controllerArray->AddLast(new FitControl(4.5,5.0, 0,2, 0.418,0.578));
      //   controllerArray->AddLast(new FitControl(5.0,5.5, 0,4, 0.418,0.578));
      //   controllerArray->AddLast(new FitControl(5.5,6.0, 0,4, 0.418,0.578));
      //   controllerArray->AddLast(new FitControl(5.5,6.0, 0,4, 0.418,0.568));
      //   controllerArray->AddLast(new FitControl(6.0,6.5, 0,4, 0.418,0.578));

    }
  // Above values all good for PERIPHERAL K0
  // For MINBIAS K0 there is a problem for the 1.4-1.6 GeV bin - this bin has possibly the most asymmetric background
  // For MIDCENTRAL K0, bins 1.4-1.6 and 1.6-1.8 GeV have problems also bin 4.5-5 GeV should probably be rebinned (peak looks split)
  // For CENTRAL K0s bin 1.4-1.6 fails. The two bins before that have S/N ratios which deviate from the trend.
  //Reduced fit range for 1.4-1.6 bin means working OK for CENTRAL now.

  //  controllerArray->AddLast(new FitControl(6.5,7.0, 0,10));
  }  else if (particleMode == 4) { //Xi case
    //controllerArray->AddLast(new FitControl(0.5,0.7, 1,1, 1.28,1.45)); //signal not visible with 
    controllerArray->AddLast(new FitControl(0.7,0.9, 1,1, 1.285,1.45));
    controllerArray->AddLast(new FitControl(0.9,1.1, 1,1, 1.285,1.45));
    controllerArray->AddLast(new FitControl(1.1,1.3, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(1.3,1.5, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(1.5,1.7, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(1.7,2.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(2.0,2.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new FitControl(2.5,3.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(3.0,3.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new FitControl(3.5,4.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new FitControl(4.0,4.5, 2,1, 1.27,1.44));
    controllerArray->AddLast(new FitControl(4.5,5.0, 2,1, 1.27,1.45));
  }
  TString *tLabel = new TString(label);
  // Special setup for doing the comparison to dEdx kaons (different binning) 
  if(tLabel->Contains("KdEdxComp")){
    controllerArray->Delete(); // removes all that went before!
    controllerArray->AddLast(new FitControl(0.25,0.35,1,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.35,0.45,1,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.45,0.55,1,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.55,0.65,1,2,0.418,0.568));
    controllerArray->AddLast(new FitControl(0.65,0.75,2,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.75,0.85,2,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.85,0.95,2,2,0.418,0.568));
    }
  // Special setup for doing comparison to Hai's K0s (different binning)
  if(tLabel->Contains("HJBins")){
    controllerArray->Delete(); // removes all that went before!
    controllerArray->AddLast(new FitControl(0.4,0.6,1,1,0.418,0.578));
    controllerArray->AddLast(new FitControl(0.6,0.8,1,1,0.418,0.568));
    controllerArray->AddLast(new FitControl(0.8,1.0,1,1,0.418,0.578));
    controllerArray->AddLast(new FitControl(1.0,1.2,2,1,0.418,0.578));
    controllerArray->AddLast(new FitControl(1.2,1.4,2,1,0.418,0.578));
    controllerArray->AddLast(new FitControl(1.4,1.6,2,1,0.418,0.568));
    controllerArray->AddLast(new FitControl(1.6,1.8,2,1,0.418,0.578));
    controllerArray->AddLast(new FitControl(1.8,2.0,2,1,0.418,0.568));
    controllerArray->AddLast(new FitControl(2.0,2.5,2,2,0.418,0.558));
    controllerArray->AddLast(new FitControl(2.5,3.0,1,2,0.418,0.568));
    controllerArray->AddLast(new FitControl(3.0,3.5,1,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(3.5,4.5,1,2,0.418,0.578));
    controllerArray->AddLast(new FitControl(4.5,6.0,0,4,0.418,0.578));
  }
  controllerArray->Sort();

  // Make the proper label to pass in for use in labelling the diagnostics
  // saved canvas files and histos
  Char_t fulllabel[80];
  sprintf(fulllabel,"%s%s",title[MultBin].Data(),label);

  //Slice up projection into mass histograms to extract yield
TH1F* hYield = (TH1F*)PtMassAna2(hPtMass,particleMode,controllerArray->GetEntries(),controllerArray,fulllabel);

// CORRECTIONS - only do number of events for now
hYield->Scale(1.0/Nev[MultBin]); //Divide by the number of events
//hYield->Scale(Veff[MultBin]); //Multiply by the vertex efficiency effectively increaing number of events 
                              //(since Veff<1) therefore decreases yield
//hYield->Scale(1.0/BR);  //Divide by branching ratio (again increases yield since BR<1)
 //hYield->Scale(1.0/(2*TMath::Pi())); // Always plot 1/2pi ...

Char_t yieldTitle[80];
sprintf(yieldTitle,"Uncorrected %s yield: %s",part,title[MultBin].Data());
hYield->SetTitle(yieldTitle);
hYield->SetXTitle("p_{t} / [GeV/c]");
hYield->SetYTitle("1/Nev.dN/dp_{t}");

// Create plots

Char_t fileNameBase[80];
sprintf(fileNameBase,"Masses%s_%s%s",part,title[MultBin].Data(),label);
Char_t fileNamePng[80];
sprintf(fileNamePng,"%s.png",fileNameBase);
Char_t fileNameEps[80];
sprintf(fileNameEps,"%s.eps",fileNameBase);
Char_t fileNamePdf[80];
sprintf(fileNamePdf,"%s.pdf",fileNameBase);

c1->SaveAs(fileNamePng);
c1->SaveAs(fileNameEps);
c1->SaveAs(fileNamePdf);

//c1->Clear();

//c1->SetLogy();
TCanvas *cYield = new TCanvas("Yield","Corrected Yield",600,400);
cYield->cd();
//cYield->SetLogy();
hYield->SetStats(kFALSE);
hYield->Draw();
 cYield->SetLogy();
 cYield->Update();
//  hRC_MB->SetMarkerStyle(20);
//  hRC_MB->SetMarkerColor(4);
//  hRC_MB->Scale(NBinMB/NBin3);
Char_t fnametext[80];
sprintf(fnametext,"Yield%s_%s%s",part,title[MultBin].Data(),label);
Char_t fnamePng[80];
sprintf(fnamePng,"%s.png",fnametext);
c1->SaveAs(fnamePng);
Char_t fnameEps[80];
sprintf(fnameEps,"%s.eps",fnametext);
c1->SaveAs(fnameEps);
  
TH1F* hScYield = hYield->Clone("ScYield");
Char_t scYieldTitle[80];
sprintf(scYieldTitle,"<N_{bin}> Scaled %s",hYield->GetTitle());
hScYield->SetTitle(scYieldTitle);
//SCALING for scaled yield only  divide by mean Nbin (scaled yield is therefore smaller)
hScYield->Scale(1/NBin[MultBin]);

Char_t fnameRoot[80];
sprintf(fnameRoot,"%s.root",fnametext);
TFile *YieldFile = new TFile(fnameRoot,"RECREATE");
hYield->Write();
hScYield->Write();
YieldFile->Close();

 controllerArray->Delete();
}
