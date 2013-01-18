#define AliV0CutVariations_cxx
// The class definition in AliV0CutVariations.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AliV0CutVariations.C")
// Root > T->Process("AliV0CutVariations.C","some options")
// Root > T->Process("AliV0CutVariations.C+")
//

#include "AliV0CutVariations.h"

#include "Riostream.h"
#include "TCanvas.h"
#include "TPDGCode.h"
#include <TDatabasePDG.h>
#include <TH2.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

AliV0CutVariations::AliV0CutVariations(TTree * /*tree*/ ) : 
fChain(0), 
fIsMC(kFALSE),
fSelectNonInjected(kFALSE),
fCMin(0.),
fCMax(90.),
fCPA(0.9975),
fDCA(1.0),
fTPCcr(70.),
fTPCcrfd(0.8),
fDCApv(0.1),
fMult(0),
fdEdx(0),
fdEdxPid(0),
fCosPA(0),
fDtrDCA(0),
fTPCrows(0),
fTPCratio(0),
fPrimDCA(0),

fK0sM(0),
fK0sSi(0),
fK0sMC(0),
fK0sAs(0),

fLambdaM(0),
fLambdaSi(0),
fLambdaMC(0),
fLambdaAs(0),

fLambdaEff(0),
fLambdaPt(0),

fLambdaBarM(0),
fLambdaBarSi(0),
fLambdaBarMC(0),
fLambdaBarAs(0),

fLambdaBarEff(0),
fLambdaBarPt(0),

fLambdaFromXi(0),
fXiM(0),
fXiSiP(0),

fLambdaBarFromXiBar(0),
fXiBarM(0),
fXiBarSiP(0)
{ }

void AliV0CutVariations::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void AliV0CutVariations::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  Init(tree);

  //cout<<tree->GetEntries()<<endl;

  TString option = GetOption();

  Int_t    lbins=100;  // number of bins in lt
  Int_t    kbins=33;  // number of bins in pt of Xi
  Double_t ltMax=100.;
  Double_t ptMax=12.;

  const Double_t xBins[]={
    0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
    1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
    2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
    4.5,5.0,5.5,6.5,8.0,10.0,12.0
  };
  const Int_t nBins=sizeof(xBins)/sizeof(Double_t) - 1;

  Double_t yBins[lbins+1]; 
  for (Int_t i=0; i<=(lbins+1); i++) yBins[i]=i*ltMax/lbins;
  Double_t zBins[kbins+1]; 
  for (Int_t i=0; i<=(kbins+1); i++) zBins[i]=i*(ptMax+2)/kbins; 


  //fMult=new TH1F("fMult","Multiplicity",1100,0.,3300);
  //fMult->GetXaxis()->SetTitle("N tracks"); 
  fMult=new TH1F("fMult","Centrality",100,0.,100);
  fMult->GetXaxis()->SetTitle("Centrality (%)"); 
  fOutput->Add(fMult);


  fCosPA=new TH1F("fCosPA","Cos(PA) distribution",50,0.9975,1.0005);
  fCosPA->GetXaxis()->SetTitle("Cos(PA)"); 
  fOutput->Add(fCosPA);

  fDtrDCA=new TH1F("fDtrDCA","DCA between V0 daughters",50,0.0,1.5);
  fDtrDCA->GetXaxis()->SetTitle("DCA (rel. u.)"); 
  fOutput->Add(fDtrDCA);

  fTPCrows=new TH1F("fTPCrows","TPC crossed pad rows",180,0.,180.);
  fTPCrows->GetXaxis()->SetTitle("TPC crossed pad rows"); 
  fOutput->Add(fTPCrows);

  fTPCratio=new TH1F("fTPCratio","TPC crossed/findable pad rows",50,0.0,1.5);
  fTPCratio->GetXaxis()->SetTitle("TPC crossed/findable pad rows"); 
  fOutput->Add(fTPCratio);

  fPrimDCA=new TH1F("fPrimDCA","DCA wrt the primary vertex",50,0.0,1.5);
  fPrimDCA->GetXaxis()->SetTitle("DCA wrt the PV (cm)"); 
  fOutput->Add(fPrimDCA);


  fdEdx=new TH2F("fdEdx","dE/dx",50,0.2,3,50,0.,6.);
  fOutput->Add(fdEdx);

  fdEdxPid=new TH2F("fdEdxPid","dE/dx with PID",50,0.2,3,50,0.,6.);
  fOutput->Add(fdEdxPid);

  fK0sM = 
  new TH2F("fK0sM", "Mass for K^{0}_{s}", 50,0.448,0.548,nBins,xBins);
  fK0sM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fK0sM);

  fK0sSi = 
  new TH2F("fK0sSi","L_{T} vs p_{T} for K^{0}_{s}, side-band subtracted",
  nBins,xBins,lbins,0.,ltMax);
  fK0sSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sSi);

  fK0sMC = 
  new TH2F("fK0sMC","L_{T} vs p_{T} for K^{0}_{s}, from MC stack", 
  nBins,xBins,lbins,0.,ltMax);
  fK0sMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sMC);

  fK0sAs = 
  new TH2F("fK0sAs", "L_{T} vs p_{T} for K^{0}_{s}, associated", 
  nBins,xBins,lbins,0.,ltMax);
  fK0sAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sAs);

  //----------------------

  fLambdaM = 
  new TH2F("fLambdaM","Mass for \\Lambda", 100, 1.065, 1.165,nBins,xBins);
  fLambdaM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fLambdaM);

 fLambdaSi = 
  new TH2F("fLambdaSi","L_{T} vs p_{T} for \\Lambda, side-band subtructed",
  nBins,xBins,lbins,0.,ltMax);
  fLambdaSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaSi);

  fLambdaMC = 
  new TH2F("fLambdaMC","L_{T} vs p_{T} for \\Lambda, from MC stack", 
  nBins,xBins,lbins,0.,ltMax);
  fLambdaMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaMC);

  fLambdaAs = 
  new TH2F("fLambdaAs","L_{T} vs p_{T} for \\Lambda, associated",
  nBins,xBins,lbins,0.,ltMax);
  fLambdaAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaAs);

  //----------------------

  fLambdaEff=fLambdaAs->ProjectionX();
  fLambdaEff->SetName("fLambdaEff");
  fLambdaEff->SetTitle("Efficiency for #Lambda");
  fOutput->Add(fLambdaEff);

  fLambdaPt=fLambdaAs->ProjectionX();
  fLambdaPt->SetName("fLambdaPt");
  fLambdaPt->SetTitle("Raw #Lambda pT spectrum");
  fOutput->Add(fLambdaPt);
  //----------------------

  fLambdaBarM = 
  new TH2F("fLambdaBarM","Mass for anti-\\Lambda",100,1.065,1.165,nBins,xBins);
  fLambdaBarM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fLambdaBarM);

  fLambdaBarSi = 
  new TH2F("fLambdaBarSi","L_{T} vs p_{T} for anti-\\Lambda, side-band subtructed",
  nBins,xBins,lbins,0.,ltMax);
  fLambdaBarSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarSi);

  fLambdaBarMC = 
  new TH2F("fLambdaBarMC","L_{T} vs p_{T} for anti-\\Lambda, from MC stack", 
  nBins,xBins,lbins,0.,ltMax);
  fLambdaBarMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarMC);

  fLambdaBarAs = 
  new TH2F("fLambdaBarAs","L_{T} vs p_{T} for anti-\\Lambda, associated",
  nBins,xBins,lbins,0.,ltMax);
  fLambdaBarAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarAs);


  //----------------------

  fLambdaBarEff=fLambdaBarAs->ProjectionX();
  fLambdaBarEff->SetName("fLambdaBarEff");
  fLambdaBarEff->SetTitle("Efficiency for anti-#Lambda");
  fOutput->Add(fLambdaBarEff);

  fLambdaBarPt=fLambdaBarAs->ProjectionX();
  fLambdaBarPt->SetName("fLambdaBarPt");
  fLambdaBarPt->SetTitle("Raw anti-#Lambda pT spectrum");
  fOutput->Add(fLambdaBarPt);

  //----------------------
  fLambdaFromXi=new TH3F("fLambdaFromXi","L_{T} vs p_{T} vs p_{T} of \\Xi for \\Lambda from \\Xi",
  nBins,xBins,lbins,yBins,kbins,zBins);
  fOutput->Add(fLambdaFromXi);

  fXiM  = 
  new TH2F("fXiM", "\\Xi mass distribution", 50, 1.271, 1.371,33,0.,ptMax+2);
  fOutput->Add(fXiM);

  fXiSiP  = new TH1F("fXiSiP", "Pt for \\Xi, side-band subracted",
  33,0.,ptMax+2);
  fOutput->Add(fXiSiP);


  fLambdaBarFromXiBar=new TH3F("fLambdaBarFromXiBar","L_{T} vs p_{T} vs p_{T} of anti-\\Xi for anti-\\Lambda from anti-\\Xi",
  nBins,xBins,lbins,yBins,kbins,zBins);
  fOutput->Add(fLambdaBarFromXiBar);

  fXiBarM  = 
  new TH2F("fXiBarM", "anti-\\Xi mass distribution", 50, 1.271, 1.371,33,0.,ptMax+2);
  fOutput->Add(fXiBarM);

  fXiBarSiP  = new TH1F("fXiBarSiP", "Pt for anti-\\Xi, side-band subracted",
  33,0.,ptMax+2);
  fOutput->Add(fXiBarSiP);

}

Bool_t AliV0CutVariations::AcceptTracks() {
  //if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  //if (t->GetKinkIndex(0)>0) return kFALSE;

  /* 
  if (fTreeVariableLeastNbrCrossedRows < fTPCcr) return kFALSE;
  if (fTreeVariableLeastRatioCrossedRowsOverFindable < fTPCcrfd) return kFALSE;
  */
  if (TMath::Abs(fTreeVariableNegEta) > 0.8) return kFALSE;
  if (TMath::Abs(fTreeVariablePosEta) > 0.8) return kFALSE;

  fTPCrows->Fill(fTreeVariableLeastNbrCrossedRows);
  fTPCratio->Fill(fTreeVariableLeastRatioCrossedRowsOverFindable);

  return kTRUE;   
}

Bool_t AliV0CutVariations::AcceptV0()
{
  //if (v0->GetOnFlyStatus()) return kFALSE;

  if (fTreeVariableV0CosineOfPointingAngle < fCPA) return kFALSE;
  if (fTreeVariableDcaV0Daughters > fDCA) return kFALSE;

  if (!AcceptTracks()) return kFALSE;

  if (TMath::Abs(fTreeVariableDcaNegToPrimVertex) < fDCApv) return kFALSE;
  if (TMath::Abs(fTreeVariableDcaPosToPrimVertex) < fDCApv) return kFALSE;

  if (fTreeVariableV0Radius < 0.9) return kFALSE;
  if (fTreeVariableV0Radius > 100) return kFALSE;

  fCosPA->Fill(fTreeVariableV0CosineOfPointingAngle);
  fDtrDCA->Fill(fTreeVariableDcaV0Daughters);
  fPrimDCA->Fill(fTreeVariableDcaNegToPrimVertex); 
  fPrimDCA->Fill(fTreeVariableDcaPosToPrimVertex);

  return kTRUE;
}

Bool_t AliV0CutVariations::AcceptPID(Int_t code) 
{
  /*
  if (code > 0) {
     if (fTreeVariablePosTransvMomentum > 1.) return kTRUE;
  } else {
     if (fTreeVariableNegTransvMomentum > 1.) return kTRUE;
  }
  */
  if (fTreeVariablePt > 1.2) return kTRUE;

  if (fIsMC) {
    // MC PID

    if (fSelectNonInjected && (!fTreeVariableIsNonInjected)) return kFALSE;

    if (code > 0) {
      if (fTreeVariablePIDPositive == code) return kTRUE;
    } else {
      if (fTreeVariablePIDNegative == code) return kTRUE;
    }
  } else {
    // Real PID
    if (code > 0) {
       if (TMath::Abs(fTreeVariableNSigmasPosProton) < 3.) return kTRUE;
    } else {
       if (TMath::Abs(fTreeVariableNSigmasNegProton) < 3.) return kTRUE;
    }
  }
  
  return kFALSE; 
}


Bool_t AliV0CutVariations::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either AliV0CutVariations::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   const Double_t yMax=0.5;
   const Double_t pMin=0.0;
   //const Double_t lMax=0.001;

   fChain->GetTree()->GetEntry(entry);

   //cout<<entry<<'\r';

   if (fTreeVariableMultiplicity<fCMin) return kFALSE;
   if (fTreeVariableMultiplicity>fCMax) return kFALSE;

   fMult->Fill(fTreeVariableMultiplicity);

   Double_t pt=0;
   Double_t lt=0;
   //+++++++ MC
   if (fIsMC) {
      if (fSelectNonInjected && (!fTreeVariableIsNonInjected)) goto real;

      Int_t code=fTreeVariablePID;

      if (code != kK0Short)
	 if (code != kLambda0)
	    if (code != kLambda0Bar) goto real;

       pt=fTreeVariablePtMC;
       if (pt < pMin) goto real;

       if (TMath::Abs(fTreeVariableRapMC) > yMax) goto real;

       //if (TMath::Abs(fTreeVariableV0CreationRadius) > lMax) goto real;
       if (fTreeVariablePrimaryStatus != 1) goto real;

       lt=fTreeVariableV0Radius;
       switch (code) {
       case kK0Short:
          fK0sMC->Fill(pt,lt);
          break;
       case kLambda0:
          fLambdaMC->Fill(pt,lt);
          break;
       case kLambda0Bar:
          fLambdaBarMC->Fill(pt,lt);
          break;
       default: break;
       }

  }
   //+++++++

 real:

   pt=fTreeVariablePt;
   if (pt < pMin) return kFALSE;

   if (!AcceptV0()) return kFALSE;

   lt=fTreeVariableV0Radius;
   if (lt/pt > 3*7.89/1.1157) return kFALSE;  

   //--- V0 switches
   Bool_t isK0s=kTRUE;
   Bool_t isLambda=kTRUE;
   Bool_t isLambdaBar=kTRUE;

   if (0.4977*lt/pt > 3*2.68) isK0s=kFALSE;
   if (1.1157*lt/pt > 3*7.89) isLambdaBar=isLambda=kFALSE;

   if (fTreeVariablePtArmV0<0.2*TMath::Abs(fTreeVariableAlphaV0)) isK0s=kFALSE;

   if (!AcceptPID(kProton)) isLambda=kFALSE;
   if (!AcceptPID(kProtonBar)) isLambdaBar=kFALSE;

   Double_t yK0s=TMath::Abs(fTreeVariableRapK0Short);
   Double_t yLam=TMath::Abs(fTreeVariableRapLambda);
   if (yK0s > yMax) isK0s=kFALSE;
   if (yLam > yMax) isLambda=isLambdaBar=kFALSE;
   //---

   Double_t mass=0., m=0., s=0.;
   Double_t peakWidth=0;
   if (isK0s) {
      mass=fTreeVariableInvMassK0s;
      fK0sM->Fill(mass,pt);

      m=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
      //s=0.0029 + 0.00077*pt;
      s=0.0041 + 0.00056*pt;
      peakWidth=5*s;
      if (TMath::Abs(m-mass) < peakWidth) {
         fK0sSi->Fill(pt,lt);
      } else {
	 isK0s=kFALSE;
      }
      if (TMath::Abs(m-mass + 1.5*peakWidth) < 0.5*peakWidth) {
         fK0sSi->Fill(pt,lt,-1);
      }
      if (TMath::Abs(m-mass - 1.5*peakWidth) < 0.5*peakWidth) {
         fK0sSi->Fill(pt,lt,-1);
      }
   }

   if (isLambda) {
      mass=fTreeVariableInvMassLambda;
      fLambdaM->Fill(mass,pt);

      m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
      //s=0.0021 + 0.00022*pt;
      s=0.0018 + 0.00030*pt;
      peakWidth = 3*s + 0.002;
      if (TMath::Abs(m-mass) < peakWidth) {
         fLambdaSi->Fill(pt,lt);
      } else {
	 isLambda=kFALSE;
      } 
      if (TMath::Abs(m-mass + 1.5*peakWidth) < 0.5*peakWidth) {
         fLambdaSi->Fill(pt,lt,-1);
      }
      if (TMath::Abs(m-mass - 1.5*peakWidth) < 0.5*peakWidth) {
         fLambdaSi->Fill(pt,lt,-1);
      }
   }

   if (isLambdaBar) {
      mass=fTreeVariableInvMassAntiLambda;
      fLambdaBarM->Fill(mass,pt);

      m=TDatabasePDG::Instance()->GetParticle(kLambda0Bar)->Mass();
      //s=0.0021 + 0.00022*pt;
      s=0.0018 + 0.00035*pt;
      peakWidth = 3*s + 0.002;
      if (TMath::Abs(m-mass) < peakWidth) {
         fLambdaBarSi->Fill(pt,lt);
      } else {
	 isLambdaBar=kFALSE;
      }
      if (TMath::Abs(m-mass + 1.5*peakWidth) < 0.5*peakWidth) {
         fLambdaBarSi->Fill(pt,lt,-1);
      }
      if (TMath::Abs(m-mass - 1.5*peakWidth) < 0.5*peakWidth) {
         fLambdaBarSi->Fill(pt,lt,-1);
      }
   }

   if (!fIsMC) return kFALSE;
   if (fSelectNonInjected && (!fTreeVariableIsNonInjected)) return kFALSE;


   //++++++ MC 
   if (!isK0s)
      if (!isLambda)
         if (!isLambdaBar) return kFALSE;//check MC only for the accepted V0s 

   Int_t code=fTreeVariablePID;
   if (code != kK0Short)
      if (code != kLambda0)
         if (code != kLambda0Bar) return kFALSE;

   Double_t ptAs=fTreeVariablePtMC;
   if (ptAs < pMin) return kFALSE;

   Double_t yAs=fTreeVariableRapMC;
   if (TMath::Abs(yAs) > yMax) return kFALSE;

   Double_t ltAs=fTreeVariableV0Radius;
      if (fTreeVariablePrimaryStatus == 1) {
	 switch (code) {
         case kK0Short:
            if (isK0s)       fK0sAs->Fill(ptAs,ltAs);
            break;
	 case kLambda0:
            if (isLambda)    fLambdaAs->Fill(ptAs,ltAs);
            break;
	 case kLambda0Bar:
            if (isLambdaBar) fLambdaBarAs->Fill(ptAs,ltAs);
            break;
         default: break;
	 }
      } else {
	 if (code==kK0Short) return kTRUE;

         Int_t xcode=fTreeVariablePIDMother;
         Double_t xpt=fTreeVariablePtXiMother;         

	 switch (code) {
	 case kLambda0:
            if (isLambda)
	    if ((xcode==kXiMinus) || (xcode==3322)) 
               fLambdaFromXi->Fill(ptAs,ltAs,xpt);
            break;
	 case kLambda0Bar:
            if (isLambdaBar)
	    if ((xcode==kXiPlusBar)||(xcode==-3322)) 
               fLambdaBarFromXiBar->Fill(ptAs,ltAs,xpt);
            break;
         default: break;
	 }
   }
 
   return kTRUE;
}

void AliV0CutVariations::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void AliV0CutVariations::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  cout<<endl;

  if (!fOutput) {
     Printf("ERROR: fOutput not available");
     return;
  }
 
  fMult = dynamic_cast<TH1F*>(fOutput->FindObject("fMult")) ; 
  if (!fMult) {
     Printf("ERROR: fMult not available");
     return;
  }

  fdEdx = dynamic_cast<TH2F*>(fOutput->FindObject("fdEdx")) ; 
  if (!fdEdx) {
     Printf("ERROR: fdEdx not available");
     return;
  }

  fdEdxPid = dynamic_cast<TH2F*>(fOutput->FindObject("fdEdxPid")) ; 
  if (!fdEdxPid) {
     Printf("ERROR: fdEdxPid not available");
     return;
  }


  fK0sMC = dynamic_cast<TH2F*>(fOutput->FindObject("fK0sMC")) ; 
  if (!fK0sMC) {
     Printf("ERROR: fK0sMC not available");
     return;
  }
  TH1D *k0sMcPx=fK0sMC->ProjectionX(); k0sMcPx->Sumw2();
  fK0sAs = dynamic_cast<TH2F*>(fOutput->FindObject("fK0sAs")) ; 
  if (!fK0sAs) {
     Printf("ERROR: fK0sAs not available");
     return;
  }
  TH1D *k0sAsPx=fK0sAs->ProjectionX(); 
  k0sAsPx->Sumw2(); //k0sAsPx->Scale(0.69);



  // Lambda histograms 
  fLambdaFromXi = dynamic_cast<TH3F*>(fOutput->FindObject("fLambdaFromXi")) ; 
  if (!fLambdaFromXi) {
     Printf("ERROR: fLambdaFromXi not available");
     return;
  }
  TH1D *lambdaFromXiPx=fLambdaFromXi->ProjectionX(); lambdaFromXiPx->Sumw2();

  fLambdaMC = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaMC")) ; 
  if (!fLambdaMC) {
     Printf("ERROR: fLambdaMC not available");
     return;
  }
  TH1D *lambdaMcPx=fLambdaMC->ProjectionX(); lambdaMcPx->Sumw2();

  fLambdaAs = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaAs")) ; 
  if (!fLambdaAs) {
     Printf("ERROR: fLambdaAs not available");
     return;
  }
  TH1D *lambdaAsPx=fLambdaAs->ProjectionX(); 
  lambdaAsPx->Sumw2(); //lambdaAsPx->Scale(0.64);

  fLambdaSi = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaSi")) ; 
  if (!fLambdaSi) {
     Printf("ERROR: fLambdaSi not available");
     return;
  }
  TH1D *lambdaSiPx=fLambdaSi->ProjectionX(); 
  lambdaSiPx->SetName("fLambdaPt");
  lambdaSiPx->Sumw2();

  fLambdaEff = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaEff")) ; 
  if (!fLambdaEff) {
     Printf("ERROR: fLambdaEff not available");
     return;
  }
  fLambdaPt = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaPt")) ; 
  if (!fLambdaPt) {
     Printf("ERROR: fLambdaPt not available");
     return;
  }


  // anti-Lambda histograms 
  fLambdaBarFromXiBar = 
  dynamic_cast<TH3F*>(fOutput->FindObject("fLambdaBarFromXiBar")) ; 
  if (!fLambdaBarFromXiBar) {
     Printf("ERROR: fLambdaBarFromXiBar not available");
     return;
  }
  TH1D *lambdaBarFromXiBarPx=
  fLambdaBarFromXiBar->ProjectionX("Bar"); lambdaBarFromXiBarPx->Sumw2();

  fLambdaBarMC = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarMC")) ; 
  if (!fLambdaBarMC) {
     Printf("ERROR: fLambdaBarMC not available");
     return;
  }
  TH1D *lambdaBarMcPx=fLambdaBarMC->ProjectionX(); lambdaBarMcPx->Sumw2();

  fLambdaBarAs = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarAs")) ; 
  if (!fLambdaBarAs) {
     Printf("ERROR: fLambdaBarAs not available");
     return;
  }
  TH1D *lambdaBarAsPx=fLambdaBarAs->ProjectionX(); 
  lambdaBarAsPx->Sumw2(); //lambdaBarAsPx->Scale(0.64);

  fLambdaBarSi = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarSi")) ; 
  if (!fLambdaBarSi) {
     Printf("ERROR: fLambdaBarSi not available");
     return;
  }
  TH1D *lambdaBarSiPx=fLambdaBarSi->ProjectionX(); 
  lambdaBarSiPx->SetName("fLambdaBarPt");
  lambdaBarSiPx->Sumw2();

  fLambdaBarEff = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaBarEff")) ; 
  if (!fLambdaBarEff) {
     Printf("ERROR: fLambdaBarEff not available");
     return;
  }
  fLambdaBarPt = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaBarPt")) ; 
  if (!fLambdaBarPt) {
     Printf("ERROR: fLambdaBarPt not available");
     return;
  }


  if (!gROOT->IsBatch()) {

    TCanvas *c1 = new TCanvas("c1","Mulitplicity");
    c1->SetLogy();
    fMult->DrawCopy() ;

    new TCanvas("c2","dE/dx");
    fdEdx->DrawCopy() ;

    new TCanvas("c3","dE/dx with PID");
    fdEdxPid->DrawCopy() ;

    if (fIsMC) {
       /*
       TH1D effK(*k0sAsPx); effK.SetTitle("Efficiency for K0s");
       effK.Divide(k0sAsPx,k0sMcPx,1,1,"b");
       new TCanvas("c4","Efficiency for K0s");
       effK.DrawCopy("E") ;
       */

       //+++ Lambdas
       fLambdaEff->Divide(lambdaAsPx,lambdaMcPx,1,1,"b");
       new TCanvas("c5","Efficiency for #Lambda");
       fLambdaEff->DrawCopy("E") ;

       lambdaSiPx->Add(lambdaFromXiPx,-1);
       lambdaSiPx->Divide(fLambdaEff);

       new TCanvas("c6","Corrected #Lambda pt");
       lambdaSiPx->SetTitle("Corrected #Lambda pt");
      *fLambdaPt = *lambdaSiPx; 
       fLambdaPt->SetLineColor(2);
       fLambdaPt->DrawCopy("E");
    
       lambdaMcPx->DrawCopy("same");
 

       //+++ anti-Lambdas
       fLambdaBarEff->Divide(lambdaBarAsPx,lambdaBarMcPx,1,1,"b");
       new TCanvas("c7","Efficiency for anti-#Lambda");
       fLambdaBarEff->DrawCopy("E") ;

       lambdaBarSiPx->Add(lambdaBarFromXiBarPx,-1);
       lambdaBarSiPx->Divide(fLambdaBarEff);

       new TCanvas("c8","Corrected anti-#Lambda pt");
       lambdaBarSiPx->SetTitle("Corrected anti-#Lambda pt");
      *fLambdaBarPt = *lambdaBarSiPx; 
       fLambdaBarPt->SetLineColor(2);
       fLambdaBarPt->DrawCopy("E");
    
       lambdaBarMcPx->DrawCopy("same");
    } else {
       new TCanvas("c6","Raw #Lambda pt");
       lambdaSiPx->SetTitle("Raw #Lambda pt");
      *fLambdaPt = *lambdaSiPx; 
       fLambdaPt->SetLineColor(2);
       fLambdaPt->DrawCopy("E");


       new TCanvas("c7","Raw anti-#Lambda pt");
       lambdaBarSiPx->SetTitle("Raw anti-#Lambda pt");
      *fLambdaBarPt = *lambdaBarSiPx; 
       fLambdaBarPt->SetLineColor(2);
       fLambdaBarPt->DrawCopy("E");
    }
  }

  TFile *file = 0;
  if (fIsMC)
    if (fSelectNonInjected) 
       file=TFile::Open("AliV0CutVariationsMC_nonInj.root", "RECREATE");
    else  
       file=TFile::Open("AliV0CutVariationsMC.root", "RECREATE");
  else
     file=TFile::Open("AliV0CutVariations.root", "RECREATE");

  fOutput->Write() ; 
  file->Close() ;
  delete file ;

}
