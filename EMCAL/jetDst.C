//*-- Author: Aleksei Pavlinov (WSU)
//*-- 
#include "jetDst.h"
#if !defined(__CINT__) || defined(__MAKECINT__) 
// These files are needed for compilation only
// (see directory /auto/u/pavlinov/macros on PDSF)
// Using with caution if you are not PAI
#include "macroIncludePai.h"
#include "macroIncludeAlice.h"
extern "C++" {void loadlibs();}
#endif

// For convinience only
class AliEMCALJetFinder; AliEMCALJetFinder* jetFinder=0;
TH2F *hLego, *hLegoW, *hTrig;
TH1F *hETrigPatch, *hJetPartPt;
TCanvas *c1=0;
// directory with file and number of files 
TString dirIn("/auto/alice/sahal/PPR3/11/"), nFileIn;
Int_t file1=1, file2=20, nevMax=0;
// change dirOut if you want.
TString dirOut("RES/FILE/"), nFileOut;
TFile *file=0;
// Default value for jet fineer - 11-feb-2202
Int_t   debugJet=1;
Float_t coneRadius = 0.3;
Float_t etSeed     = 4.;
Float_t minJetEt   = 40.;
Float_t minCellEt  = 1.; // be carefull
Float_t minPtCut   = 2.; // for including charge particles to file 
// For BG subtraction
Int_t modeSubtraction = 1;
Float_t minMove = 0.05;
Float_t maxMove = 0.15;
Float_t precBg  = 0.035;
// key 13-feb-2002
Bool_t hadronCorrection = kTRUE;
//
class AliEMCALJetMicroDst; AliEMCALJetMicroDst *dst=0;
Int_t nevSum = 0, nevF;
TArrayI nevInFiles(100);


void jetDst(Int_t mode, Int_t var, Int_t nmax) 
{ 
//
// mode - variable for selection input mode (see function defineInput());
// var  - variable for selection jetfinder mode (see function defineJetFinder()).
//
// Dynamically link some shared libs                    
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("../macros/loadlibs.C");
        loadlibs();
    }

    if(!defineInput(mode)) return;

    dst = new AliEMCALJetMicroDst;
    dst->SetDebug(1);
    if(!nFileOut.Contains(".root")){
      nFileOut += var;
      nFileOut += ".root";
    }
    dst->Create(nFileOut.Data()); 
    gROOT->GetListOfBrowsables()->Add(dst);

    for(Int_t j=file1; j<=file2; j++){
      nevF = 0;
      if(file2>1) {// multiple file
        nFileIn = dirIn; 
        if       (mode>=11 && mode<=11) {
          nFileIn += j;
          nFileIn += ("/galice.root");
        } else if(mode>=21 && mode<=51){
          nFileIn += "galice"; 
          if(j>0) nFileIn += j;
          nFileIn += ".root"; 
        }
      }
// Connect the Root Galice file containing Geometry, Kine and Hits
      file = new TFile(nFileIn.Data(), "READ");
      printf("\n FILE %i -> %s\n", j, nFileIn.Data());
    
      if (file == 0) continue; 
    
// Get AliRun object from file or create it if not on file
      gAlice = 0;
      gAlice = (AliRun*)(file->Get("gAlice"));
      if(!gAlice)  return; 
      printf("AliRun object found on file\n");

      defineJetFinder(mode,var);
//
//   Loop over events in file 
//

//      if(!c1) c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
      Int_t nhit=0;
      TString sw;
      if(nmax>0) nevMax = nmax;
      for (Int_t nev = 0; nev<= nevMax; nev++) {
	jetFinder->ResetJets();	
	Int_t nparticles = gAlice->GetEvent(nev);
	if (nparticles <= 0) break; // last event
        nevSum++;
        nevF++;
// Partons jet
/*
        jetFinder->FillFromPartons();
	jetFinder->Find();
        for(Int_t ij=0; ij<jetFinder->Njets(); ij++){
          hJetPartPt->Fill(jetFinder->JetEnergy(ij));
        }
	jetFinder->ResetJets();	
*/
// ECAL information
	jetFinder->FillFromHits();

	//        gROOT->cd(); be carefull with this staff
        hLegoW = new TH2F((*jetFinder->GetLego())); // for saving
        sw = hLegoW->GetName(); 
        sw += nev; 
        hLegoW->SetName(sw.Data());
        sw = "Energy in ECAL:Event"; 
        sw += nev; 
        hLegoW->SetTitle(sw.Data());

        if(hTrig==0) {
           hLego = jetFinder->GetLego();
           hTrig = bookTrigHist(hLego);
        }
        fillTriggerPatch(hLegoW, hTrig, hETrigPatch);

//	jetFinder->FillFromDigits();
// TPC  information
	jetFinder->FillFromTracks(1, 0);
//                                ^ preserves info from hit

// TPC  information from Hits associated to the EMCAL
//      jetFinder->FillFromHitFlaggedTracks(1);
//  
 	jetFinder->Find();
        dst->Fill(gAlice, jetFinder);
//
// Look at results
	printf("\n Events : %d  Jets : %d \n", nev, jetFinder->Njets());
	Int_t njet = jetFinder->Njets();
	for (Int_t nj=0; nj<njet; nj++)
	{
	    printf("\n Jet Energy %5d %8.2f \n", 
		   nj, jetFinder->JetEnergy(nj));
	    printf("Jet Phi    %5d %8.2f %8.2f \n", 
		   nj, jetFinder->JetPhiL(nj), jetFinder->JetPhiW(nj));
	    printf("Jet Eta    %5d %8.2f %8.2f \n", 
		   nj, jetFinder->JetEtaL(nj), jetFinder->JetEtaW(nj));
	}

        Int_t in;
        if (c1 && nev <-10 && njet==1 && jetFinder->JetEnergy(0)>150) {
	   c1->cd();
	   cout << "Event "<< nev <<endl;
           hLego->Draw("lego");
           c1->Update();
           cin  >> in;
           if(in<=0) return;
        }
      } // event loop
      nevInFiles[j-1] = nevF;
      if (gAlice) {
         delete gAlice;  // don't delete after closing file
         gAlice = 0;
      }
      file->Close();
      if(c1) hLego->Draw("lego");
    } // file loop
    if(nmax==0) dst->Close();
    nevInFiles.Set(file2);
    printf("\n Summ.events %6i -> %6i\n", nevSum, (Int_t) nevInFiles.GetSum());

    printf("Events    ");  
    for(Int_t j=file1; j<=file2; j++){
      printf(" %i4 ",nevInFiles[j-1]);
      if((j-file1)%10==0 && j!=file1) printf("\n          ");
    }
}

// ###################################################################
TString* defineInput(Int_t mode)
{// 7-feb-2002
  Char_t ctmp[200];
  switch(mode){
// case -5 needing for testing after merging
  case -5: 
  case -4:
    hadronCorrection = kFALSE; // for checking Tom's idea
  case -3: 
  case -2: 
  case -1:
    if     (mode==-2) minPtCut = 1.; 
    else if(mode<=-3 || mode>=-5) {
      minCellEt = 0.0; // 13-feb-2002
      minPtCut  = 0.0;
    }
 
    coneRadius = 0.5;
    minJetEt   = 20.0;
    file1=1; 
    file2=1;
    modeSubtraction = 0;
    nevMax = 1;
    nevMax = 999;
    dirIn    = ".";
    nFileIn  = "galice.root";
    sprintf(ctmp,"PythiaMode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);
    nFileOut = dirOut + ctmp;
    if(!hadronCorrection) nFileOut.ReplaceAll(".root","_NoHadrCorr.root"); 
    break; 
  case 1: 
    nFileOut = dirOut + "dstTest.root"; 
    break;

  case 11: 
    // jj100 41 events per file - see Log file 
    dirIn = "/auto/alice/sahal/PPR4/1/";
    file1 = 1;
    file2 = 30; // max 75
    nevMax = 41;
    sprintf(ctmp,"HijingMode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;
    //=============================================================
    // central events with jet trigger - 15 feb 2002
    //=============================================================
  case 21: 
    dirIn  = "/auto/alice/sahal/kHijing_jj50/";
    file1  = 0;
    file2  = 9;
    nevMax = 41;
    sprintf(ctmp,"Hijing_jj50Mode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;
  case 22: 
    dirIn  = "/auto/alice/sahal/kHijing_jj75/";
    file1  = 0;
    file2  = 9;
    nevMax = 41;
    sprintf(ctmp,"Hijing_jj75Mode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;
  case 23: 
    dirIn  = "/auto/alice/sahal/kHijing_jj100/";
    file1  = 0;
    file2  = 4;
    nevMax = 41;
    sprintf(ctmp,"Hijing_jj100Mode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;
  case 24: 
    dirIn  = "/auto/alice/sahal/kHijing_jj125/";
    file1  = 0;
    file2  = 4;
    nevMax = 41;
    sprintf(ctmp,"Hijing_jj125Mode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;

  case 51: 
    // BG - HIJING events
    dirIn = "/auto/alice/sahal/background/";
    file1 = 1;
    file2 = 9;
    nevMax = 100; // ??
    sprintf(ctmp,"HijingBGMode%iR%3.1fSeed%3.1fJEt%4.1fCEt%3.1fPtCut%3.1f.root",
	    mode, coneRadius, etSeed, minJetEt, minCellEt, minPtCut);  
    nFileOut = dirOut + ctmp; 
    break;
  default: 
     printf("\n<I> Wrong input mode => %i\ncout", mode);
     return 0; 
  }
  printf("<I> dir %s file1 %3i file2 %3i\n", dirIn.Data(), file1, file2);
  return &nFileOut; 
}

void  defineJetFinder(Int_t mode, Int_t var)
{// 11-feb-2002
//  Create and configure jet finder
//  gAlice must be defined already !!! 
  if(!jetFinder) {
    jetFinder =  new AliEMCALJetFinder("UA1 Jet Finder", "Pre Production");
    jetFinder->Init();

    jetFinder->SetDebug(debugJet);
    jetFinder->SetConeRadius(coneRadius);
    jetFinder->SetEtSeed(etSeed);
    jetFinder->SetMinJetEt(minJetEt);
    jetFinder->SetMinCellEt(minCellEt);
    jetFinder->SetPtCut(minPtCut);

    if(var>=2) jetFinder->SetMomentumSmearing(1);
    if(var>=3) jetFinder->SetEfficiencySim(1);

    jetFinder->SetHadronCorrection(0);
    if(var>=4 && hadronCorrection) {
      jetFinder->SetHadronCorrection(1);
      jetFinder->SetHadronCorrector(AliEMCALHadronCorrectionv0::Instance());
    }
    if(var>=5) { // 11-feb-2002
      jetFinder->SetParametersForBgSubtraction
      (modeSubtraction,minMove,maxMove,precBg); 
    }

    if(mode<0) jetFinder->SetSamplingFraction(1.);//must use only for PYTHIA file
    else       jetFinder->SetSamplingFraction(12.9);
    gROOT->GetListOfBrowsables()->Add(jetFinder);
  }
}

void   testDST(char* fname)
{// 8-feb-2002
  dst = new AliEMCALJetMicroDst;
  dst->Open(fname);
  dst->Test();
  gROOT->GetListOfBrowsables()->Add(dst);
}

//----------------- 4-feb-2002 --------------
TH2F* bookTrigHist(TH2F *hid)
{
  if(hid==0) return 0;
  TAxis *xax  = hid->GetXaxis();
  Double_t xmin = xax->GetXmin();
  Double_t xmax = xax->GetXmax();
  Int_t nx = 6;
  TAxis *yax  = hid->GetYaxis();
  Double_t ymin = yax->GetXmin();
  Double_t ymax = yax->GetXmax();
  Int_t ny = 6;
  TH2F *hout = new TH2F("hTrig","Energy grid for trigger patch",
  nx,xmin,xmax, ny,ymin,ymax);
  hout->AddDirectory(0);
  
  // time solution
  hETrigPatch = new TH1F("hETrigPatch","Energy in trigger patch",300, 0.0, 300.);
  hJetPartPt  = new TH1F("hJetPartPt", "E_{T} for partons jet", 300, 0.0, 300.);
  hETrigPatch->AddDirectory(0);
  return hout;
}

void fillTriggerPatch(TH2F *hid, TH2F *htrig, TH1F *hE)
{
  TAxis *xax  = hid->GetXaxis();
  TAxis *yax  = hid->GetYaxis();
  Int_t nx = xax->GetNbins();
  Int_t ny = yax->GetNbins();

  htrig->Reset();
  Double_t x, y, e;
  for(Int_t ix=1; ix<=nx; ix++){
     x = xax->GetBinCenter(ix);
     for(Int_t iy=1; iy<=ny; iy++){
        y = yax->GetBinCenter(iy);
        e = hid->GetBinContent(ix, iy);
        htrig->Fill(x,y,e);
     } 
  }

  xax  = htrig->GetXaxis();
  yax  = htrig->GetYaxis();
  nx = xax->GetNbins();
  ny = yax->GetNbins();
  for(Int_t ix=1; ix<=nx; ix++){
     for(Int_t iy=1; iy<=ny; iy++){
       hE->Fill(htrig->GetBinContent(ix, iy)); 
     }
  }
}
