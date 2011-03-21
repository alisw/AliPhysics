#include "AliMassFitControl.h"

// Two functions now.
// MultYields3 takes a 3D histogram and calls the 2D function MultYields2
// MultYields2 could be called directly if we only have a 2D histogram

void MultYields3(TH3F *PtMassMult, Int_t particleMode, Int_t MultBin, Char_t* label){

  // Leave open possibility to choose different values depending on MultBin - LSB
  Float_t MultLo[1] = {0};
  Float_t MultHi[1] = {300};

  //Make 2D projections from the 3D histogram
  PtMassMult->GetZaxis()->SetRange(MultLo[MultBin],MultHi[MultBin]);
  TH2F* hPtMass = (TH2F*)PtMassMult->Project3D("XY");// FIX:MF
  hPtMass->SetTitle("PtMass");
  
  //hPtMass->Draw(); // Drawing invokes default c1 canvas making it inaccessible in MultYields2 - TO FIX
  MultYields2(hPtMass, particleMode,MultBin,label);
  
}

void MultYields2(TH2F *hPtMass, Int_t particleMode, Int_t MultBin, Char_t* label){

  hPtMass->Draw();
  
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

  
 
 
TString title[1]={"MinimumBias"};
  //Make 2D projections from the 3D histogram
  //Minbias (i.e. everything)
  
  TObjArray *controllerArray = new TObjArray(40,1); // 2nd arg means can count from 1!

  //Here probably need switch-case, depending on mult bin and particle
  // LoPt, HiPt, polynomial order, rebinning factor

  /////  LAMBDA and LAMBDA+ANTILAMBDA (Combination)

  if(particleMode == 1 || particleMode == 3){ // Lambda or Lambda+Anti-Lambda
  //  controllerArray->AddLast(new AliMassFitControl(0.1,0.2, 1,2));
  //  controllerArray->AddLast(new AliMassFitControl(0.2,0.3, 1,2));
  //  controllerArray->AddLast(new AliMassFitControl(0.3,0.4, 1,1, 1.096,1.16))
  //  controllerArray->AddLast(new AliMassFitControl(0.4,0.5, 2,1, 1.095,1.15));
  //  controllerArray->AddLast(new AliMassFitControl(0.5,0.6, 2,1, 1.095,1.15));
  //  controllerArray->AddLast(new AliMassFitControl(0.6,0.7, 2,1, 1.095,1.15));
  //  controllerArray->AddLast(new AliMassFitControl(0.7,0.8, 2,1, 1.095,1.15));
      controllerArray->AddLast(new AliMassFitControl(0.8,1.0, 2,1, 1.095,1.15));
      controllerArray->AddLast(new AliMassFitControl(1.0,1.2, 2,1, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(1.2,1.4, 2,1));
      controllerArray->AddLast(new AliMassFitControl(1.4,1.6, 2,1));
      controllerArray->AddLast(new AliMassFitControl(1.6,1.8, 2,1));
      controllerArray->AddLast(new AliMassFitControl(1.8,2.0, 2,1));
      controllerArray->AddLast(new AliMassFitControl(2.0,2.2, 2,1, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(2.2,2.4, 2,1, 1.095,1.16));
      controllerArray->AddLast(new AliMassFitControl(2.4,2.6, 2,1, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(2.6,2.8, 1,2, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(2.8,3.0, 1,2, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(3.0,3.5, 1,2, 1.098,1.16));
      controllerArray->AddLast(new AliMassFitControl(3.5,4.0, 1,2, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 1,2, 1.095,1.17));
      controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 1,2, 1.095,1.165));
      controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 1,2, 1.095,1.17));
 //     controllerArray->AddLast(new AliMassFitControl(5.5,6.0, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(6.0,6.5, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(6.5,7.0, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(7.0,7.5, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(7.5,8.0, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(8.0,9.0, 1,2, 1.095,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(9.0,10.0, 1,2, 1.090,1.17));
 //   controllerArray->AddLast(new AliMassFitControl(10.0,12.0, 1,2, 1.082,1.165));
 }
  /// ANTI LAMBDA ---->
  else if (particleMode == 2){ // Anti-Lambdas
  controllerArray->AddLast(new AliMassFitControl(0.3,0.4, 2,1, 1.095,1.15));
  controllerArray->AddLast(new AliMassFitControl(0.4,0.5, 2,1, 1.095,1.15));
  controllerArray->AddLast(new AliMassFitControl(0.5,0.6, 2,1, 1.095,1.15));
  controllerArray->AddLast(new AliMassFitControl(0.6,0.7, 2,1, 1.095,1.15));
  controllerArray->AddLast(new AliMassFitControl(0.7,0.8, 2,1, 1.090,1.15));
  controllerArray->AddLast(new AliMassFitControl(0.8,1.0, 2,1, 1.095,1.15));
  controllerArray->AddLast(new AliMassFitControl(1.0,1.2, 2,1, 1.090,1.17));
  controllerArray->AddLast(new AliMassFitControl(1.2,1.4, 2,1, 1.090,1.16));
  controllerArray->AddLast(new AliMassFitControl(1.4,1.6, 2,1, 1.090,1.17));
  controllerArray->AddLast(new AliMassFitControl(1.6,1.8, 2,1, 1.090,1.17));
  controllerArray->AddLast(new AliMassFitControl(1.8,2.0, 2,1, 1.090,1.17));
  controllerArray->AddLast(new AliMassFitControl(2.0,2.2, 2,1, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(2.2,2.4, 2,1, 1.095,1.16));
  controllerArray->AddLast(new AliMassFitControl(2.4,2.6, 1,1, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(2.6,2.8, 1,2, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(2.8,3.0, 1,2, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(3.0,3.5, 1,2, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(3.5,4.0, 1,2, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 1,2, 1.090,1.17));
  controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 0,4, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 0,4, 1.095,1.17));
  controllerArray->AddLast(new AliMassFitControl(5.5,6.0, 0,6, 1.085,1.17));
  } // end if anti-Lambda
  else if (particleMode == 0){ // K0s case
      //controllerArray->AddLast(new AliMassFitControl(0.1,0.2, 1,2, 0.44,0.56));
      //controllerArray->AddLast(new AliMassFitControl(0.2,0.3, 1,2, 0.44,0.56));
      //controllerArray->AddLast(new AliMassFitControl(0.3,0.4, 1,1, 0.442,0.559));
      //controllerArray->AddLast(new AliMassFitControl(0.4,0.5, 1,1, 0.44,0.56));
      //controllerArray->AddLast(new AliMassFitControl(0.5,0.6, 1,1, 0.442,0.559));
      //controllerArray->AddLast(new AliMassFitControl(0.6,0.7, 1,1, 0.442,0.558));
      //controllerArray->AddLast(new AliMassFitControl(0.7,0.8, 1,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(0.8,1.0, 2,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(1.0,1.2, 2,1, 0.45,0.56));
      controllerArray->AddLast(new AliMassFitControl(1.2,1.4, 2,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(1.4,1.6, 2,1, 0.45,0.56));
      controllerArray->AddLast(new AliMassFitControl(1.6,1.8, 2,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(1.8,2.0, 1,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(2.0,2.2, 1,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(2.2,2.4, 1,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(2.4,2.6, 1,1, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(2.6,2.8, 1,2, 0.44,0.54));
      controllerArray->AddLast(new AliMassFitControl(2.8,3.0, 1,2, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(3.0,3.5, 1,2, 0.44,0.56));
      controllerArray->AddLast(new AliMassFitControl(3.5,4.0, 1,2, 0.46,0.56));
      controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 1,2, 0.46,0.56));
      controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 1,2, 0.418,0.52));
      //controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 0,4, 0.418,0.53));
//      controllerArray->AddLast(new AliMassFitControl(5.5,6.0, 1,4, 0.418,0.54));
//      controllerArray->AddLast(new AliMassFitControl(6.0,6.5, 1,4, 0.418,0.56));
      //controllerArray->AddLast(new AliMassFitControl(6.5,7.0, 0,4, 0.418,0.578));
  }  else if (particleMode == 4) { //Xi case
    //controllerArray->AddLast(new AliMassFitControl(0.5,0.7, 1,1, 1.28,1.45)); //signal not visible with 
    controllerArray->AddLast(new AliMassFitControl(0.7,0.9, 1,1, 1.285,1.45));
    controllerArray->AddLast(new AliMassFitControl(0.9,1.1, 1,1, 1.285,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.1,1.3, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.3,1.5, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.5,1.7, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.7,2.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(2.0,2.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(2.5,3.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(3.0,3.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(3.5,4.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 2,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 2,1, 1.27,1.45));
  }
  controllerArray->Sort();

  // Make the proper label to pass in for use in labelling the diagnostics
  // saved canvas files and histos
  Char_t fulllabel[80];
  sprintf(fulllabel,"%s%s",title[MultBin].Data(),label);

  //Slice up projection into mass histograms to extract yield
TH1F* hYield = (TH1F*)PtMassAna2(hPtMass,particleMode,controllerArray->GetEntries(),controllerArray,fulllabel);

// CORRECTIONS
  Int_t Nev = 2800;  // Obviously we need a better way to access the proper number FIX ME
 hYield->Scale(1.0/Nev); //Divide by the number of events
// Better not to hide these other corrections deep in here
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

// This section for yield scaled by number of binary collisions.
// Could add array of values and do scaling according to 'MultBin' index
TH1F* hScYield = hYield->Clone("ScYield");
Char_t scYieldTitle[80];
sprintf(scYieldTitle,"<N_{bin}> Scaled %s",hYield->GetTitle());
hScYield->SetTitle(scYieldTitle);
//SCALING for scaled yield only  divide by mean Nbin (scaled yield is therefore smaller)
//hScYield->Scale(1/NBin[MultBin]);

Char_t fnameRoot[80];
sprintf(fnameRoot,"%s.root",fnametext);
TFile *YieldFile = new TFile(fnameRoot,"RECREATE");
hYield->Write();
hScYield->Write();
YieldFile->Close();

 controllerArray->Delete();
}
