// $Id$
//
// Jet spectrum task.
//
// Author: R.Reed, M.Connors

#include "AliAnalysisTaskEmcalJetSpectra.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"
#include "AliRhoParameter.h"
#include "AliEmcalParticle.h"
#include "AliLocalRhoParameter.h"
#include "AliAnalysisTaskLocalRho.h"

ClassImp(AliAnalysisTaskEmcalJetSpectra)

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra() : 
  AliAnalysisTaskEmcalJet("spectra",kFALSE), 
  fLocalRhoVal(0),
  fHistRhovsCent(0),
  fHistNjetvsCent(0),
  fHistGLvsLOCrho(0),
  fHistRhovsdEPLOC(0), fHistRhovsdEPGL(0),
  fHistJetPtvsdEPLOC(0), fHistJetPtvsdEPGL(0),
  fHistRhovsEPLOC(0), fHistRhovsEPGL(0),
  fHistJetPtvsEPLOC(0), fHistJetPtvsEPGL(0),
  fHistCorJetPt(0), fHistCorJetPtGL(0)
{
  // Default constructor.
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistJetPtvsEP[i]           = 0;
    fHistJetPtvsEPBias[i]       = 0;
    fHistRhovsEP[i]             = 0;
    fHistCorJetPtfromLocalRho[i]= 0;
    fHistCorJetPtfromGlobalRho[i] = 0;

    fHistCorJetPtfromLocalRhoIN[i]  = 0;
    fHistCorJetPtfromLocalRhoOUT[i] = 0;
    fHistCorJetPtfromGlobalRhoIN[i] = 0;
    fHistCorJetPtfromGlobalRhoOUT[i]= 0;
    fHistRhodEPcentLOC[i]       = 0;
    fHistRhodEPcentGL[i]        = 0;
    fHistCorJetPtdEPcentLOC[i]  = 0;
    fHistCorJetPtdEPcentGL[i]   = 0;
    fHistRhoEPcentLOC[i]        = 0;
    fHistRhoEPcentGL[i]         = 0;
    fHistCorJetPtEPcentLOC[i]   = 0;
    fHistCorJetPtEPcentGL[i]    = 0;

  }
  
  SetMakeGeneralHistograms(kTRUE);
    
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fLocalRhoVal(0),
  fHistRhovsCent(0),
  fHistNjetvsCent(0),
  fHistGLvsLOCrho(0),
  fHistRhovsdEPLOC(0), fHistRhovsdEPGL(0),
  fHistJetPtvsdEPLOC(0), fHistJetPtvsdEPGL(0),
  fHistRhovsEPLOC(0), fHistRhovsEPGL(0),
  fHistJetPtvsEPLOC(0), fHistJetPtvsEPGL(0),
  fHistCorJetPt(0), fHistCorJetPtGL(0)
 { 
   for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistJetPtvsEP[i]           = 0;
    fHistJetPtvsEPBias[i]       = 0;
    fHistRhovsEP[i]             = 0;
    fHistCorJetPtfromLocalRho[i]= 0;
    fHistCorJetPtfromGlobalRho[i] = 0;

    fHistCorJetPtfromLocalRhoIN[i]  = 0;
    fHistCorJetPtfromLocalRhoOUT[i] = 0;
    fHistCorJetPtfromGlobalRhoIN[i] = 0;
    fHistCorJetPtfromGlobalRhoOUT[i]= 0;
    fHistRhodEPcentLOC[i]       = 0;
    fHistRhodEPcentGL[i]        = 0;
    fHistCorJetPtdEPcentLOC[i]  = 0;
    fHistCorJetPtdEPcentGL[i]   = 0;
    fHistRhoEPcentLOC[i]        = 0;
    fHistRhoEPcentGL[i]         = 0;
    fHistCorJetPtEPcentLOC[i]   = 0;
    fHistCorJetPtEPcentGL[i]    = 0;

   }
 
   SetMakeGeneralHistograms(kTRUE);
 }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSpectra::UserCreateOutputObjects()
{

  if (! fCreateHisto) return;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRhovsCent             = new TH2F("RhovsCent",              "RhovsCent",             100, 0.0, 100.0, 500, 0, 500);
  fHistNjetvsCent            = new TH2F("NjetvsCent",             "NjetvsCent",            100, 0.0, 100.0, 100, 0, 100);

  fHistGLvsLOCrho            = new TH2F("GLvsLOCrho", "GLvsLOCrho", 400, 0.0, 400.0, 400, 0.0, 400.0);

  fHistRhovsdEPLOC           = new TH2F("RhovsdEPLOC", "RhovsdEPLOC",400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistRhovsdEPGL            = new TH2F("RhovsdEPGL", "RhovsdEPGL",400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistJetPtvsdEPLOC         = new TH2F("JetPtvsdEPLOC", "JetPtvsdEPLOC", 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistJetPtvsdEPGL          = new TH2F("JetPtvsdEPGL", "JetPtvsdEPGL", 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi()); 
  fHistRhovsEPLOC            = new TH2F("RhovsEPLOC", "RhovsEPLOC",400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistRhovsEPGL             = new TH2F("RhovsEPGL", "RhovsEPGL",400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistJetPtvsEPLOC          = new TH2F("JetPtvsEPLOC", "JetPtvsEPLOC", 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
  fHistJetPtvsEPGL           = new TH2F("JetPtvsEPGL", "JetPtvsEPGL", 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());

  fHistCorJetPtGL            = new TH1F("NjetvsCorrJetPtGL", "NjetvsCorrJetPtGL", 500, -250, 250);
  fHistCorJetPt              = new TH1F("NjetvsCorrJetPt", "NjetvsCorrJetPt", 500, -250, 250);

  TString name;
  TString title;
  for (Int_t i = 0;i<6;++i){
    name = TString(Form("JetPtvsTrackPt_%i",i));
    title = TString(Form("Jet pT vs Leading Track pT cent bin %i",i));
    fHistJetPtvsTrackPt[i] = new TH2F(name,title,1000,-500,500,100,0,100);
    fOutput->Add(fHistJetPtvsTrackPt[i]);
    name = TString(Form("RawJetPtvsTrackPt_%i",i));
    title = TString(Form("Raw Jet pT vs Leading Track pT cent bin %i",i));
    fHistRawJetPtvsTrackPt[i] = new TH2F(name,title,1000,-500,500,100,0,100);
    fOutput->Add(fHistRawJetPtvsTrackPt[i]);
    name = TString(Form("TrackPt_%i",i));
    title = TString(Form("Track pT cent bin %i",i));
    fHistTrackPt[i] = new TH1F(name,title,1000,0,200);
    fOutput->Add(fHistTrackPt[i]);
   
    name = TString(Form("EP0_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0[i]);
    name = TString(Form("EP0A_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0A[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0A[i]);
    name = TString(Form("EP0C_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0C[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0C[i]);
    name = TString(Form("EPAvsC_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEPAvsC[i] = new TH2F(name,title,100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEPAvsC[i]);
    name = TString(Form("JetPtvsdEP_%i",i));
    title = TString(Form("Jet pt vs dEP cent bin %i",i));
    fHistJetPtvsdEP[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEP[i]);
    name = TString(Form("JetPtvsdEPBias_%i",i));
    title = TString(Form("Bias Jet pt vs dEP cent bin %i",i));
    fHistJetPtvsdEPBias[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEPBias[i]);
    name = TString(Form("JetPtvsEP_%i",i));
    title = TString(Form("Jet pt vs EP cent bin %i",i));
    fHistJetPtvsEP[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsEP[i]);
    name = TString(Form("JetPtvsEPBias_%i",i));
    title = TString(Form("Bias Jet pt vs EP cent bin %i",i));
    fHistJetPtvsEPBias[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsEPBias[i]);
    name = TString(Form("RhovsEP_%i",i));
    title = TString(Form("Rho vs EP cent bin %i",i));
    fHistRhovsEP[i] = new TH2F(name,title,500,0,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistRhovsEP[i]);
      
      name = TString(Form("NjetvsCorrJetPtfromLocalRho_%i",i));
      title = TString(Form("Njets vs Corrected jet pT from Local Rho cent bin %i",i));
      fHistCorJetPtfromLocalRho[i] = new TH1F(name,title, 500, -250,250);
      fOutput->Add(fHistCorJetPtfromLocalRho[i]);
  
      name = TString(Form("NjetvsCorrJetPtfromGlobalRho_%i",i));
      title = TString(Form("Njets vs Corrected jet pT from Global Rho cent bin %i",i));
      fHistCorJetPtfromGlobalRho[i] = new TH1F(name,title, 500, -250,250);
      fOutput->Add(fHistCorJetPtfromGlobalRho[i]);

      name = TString(Form("NjetvsCorrJetPtfromGlobalRhoIN_%i",i));
      title = TString(Form("Njets vs Corrected jet pT from Global Rho IN PLANE cent bin %i",i));
      fHistCorJetPtfromGlobalRhoIN[i] = new TH1F(name,title, 500, -250,250);
      fOutput->Add(fHistCorJetPtfromGlobalRhoIN[i]);

      name = TString(Form("NjetvsCorrJetPtfromGlobalRhoOUT_%i",i));
      title = TString(Form("Njets vs Corrected jet pT from Global Rho OUT PLANE cent bin %i",i));
      fHistCorJetPtfromGlobalRhoOUT[i] = new TH1F(name,title, 500, -250, 250);
      fOutput->Add(fHistCorJetPtfromGlobalRhoOUT[i]);

      name = TString(Form("RhovsdEPcentGLOBAL_%i",i));
      title = TString(Form("Rho vs delta Event Plane angle for Global Rho cent bin %i",i));
      fHistRhodEPcentGL[i] = new TH2F(name,title,400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
      fOutput->Add(fHistRhodEPcentGL[i]);

      name = TString(Form("JePtvsdEPcentGLOBAL_%i",i));
      title = TString(Form("Jet Pt vs delta Event Plane angle for Global Rho cent bin %i",i));
      fHistCorJetPtdEPcentGL[i] = new TH2F(name,title, 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
      fOutput->Add(fHistCorJetPtdEPcentGL[i]);

    name = TString(Form("NjetvsCorrJetPtfromLocalRhoIN_%i",i));
    title = TString(Form("Njets vs Corrected jet pT from Local Rho IN PLANE cent bin %i",i));
    fHistCorJetPtfromLocalRhoIN[i] = new TH1F(name,title, 500, -250, 250);
    fOutput->Add(fHistCorJetPtfromLocalRhoIN[i]);

    name = TString(Form("NjetvsCorrJetPtfromLocalRhoOUT_%i",i));
    title = TString(Form("Njets vs Corrected jet pT from Local Rho OUT PLANE cent bin %i",i));
    fHistCorJetPtfromLocalRhoOUT[i] = new TH1F(name,title, 500, -250, 250);
    fOutput->Add(fHistCorJetPtfromLocalRhoOUT[i]);

    name = TString(Form("RhovsdEPcentLOCAL_%i",i));
    title = TString(Form("Rho vs delta Event Plane angle for Local Rho cent bin %i",i));
    fHistRhodEPcentLOC[i] = new TH2F(name,title,400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
    fOutput->Add(fHistRhodEPcentLOC[i]);

    name = TString(Form("JetPtvsdEPcentLOCAL_%i",i));
    title = TString(Form("Jet Pt vs delta Event Plane angle for Local Rho cent bin %i",i));
    fHistCorJetPtdEPcentLOC[i] = new TH2F(name,title, 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
    fOutput->Add(fHistCorJetPtdEPcentLOC[i]);


      name = TString(Form("RhovsEPcentGLOBAL_%i",i));
      title = TString(Form("Rho vs Event Plane angle for Global Rho cent bin %i",i));
      fHistRhoEPcentGL[i] = new TH2F(name,title,400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
      fOutput->Add(fHistRhoEPcentGL[i]);

    name = TString(Form("RhovsEPcentLOCAL_%i",i));
    title = TString(Form("Rho vs Event Plane angle for Local Rho cent bin %i",i));
    fHistRhoEPcentLOC[i] = new TH2F(name,title,400,0,400, 144,-1*TMath::Pi(),1*TMath::Pi());
    fOutput->Add(fHistRhoEPcentLOC[i]);

    name = TString(Form("JetPtvsEPcentLOCAL_%i",i));
    title = TString(Form("Jet Pt vs Event Plane angle for Local Rho cent bin %i",i));
    fHistCorJetPtEPcentLOC[i] = new TH2F(name,title, 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
    fOutput->Add(fHistCorJetPtEPcentLOC[i]);


      name = TString(Form("JetPtvsEPcentGLOBAL_%i",i));
      title = TString(Form("Jet Pt vs Event Plane angle for Global Rho cent bin %i",i));
      fHistCorJetPtEPcentGL[i] = new TH2F(name,title, 500, -250, 250, 144,-1*TMath::Pi(),1*TMath::Pi());
      fOutput->Add(fHistCorJetPtEPcentGL[i]);

  }
 
  fOutput->Add(fHistRhovsCent);
  fOutput->Add(fHistNjetvsCent);

  fOutput->Add(fHistGLvsLOCrho);

  fOutput->Add(fHistRhovsdEPLOC);
  fOutput->Add(fHistRhovsdEPGL);
  fOutput->Add(fHistJetPtvsdEPLOC);
  fOutput->Add(fHistJetPtvsdEPGL);

  fOutput->Add(fHistRhovsEPLOC);
  fOutput->Add(fHistJetPtvsEPLOC);

  fOutput->Add(fHistCorJetPtGL);
  fOutput->Add(fHistCorJetPt);
  fOutput->Add(fHistRhovsEPGL);
  fOutput->Add(fHistJetPtvsEPGL);

   PostData(1, fOutput);
}

//________________________________________________________________________

Int_t AliAnalysisTaskEmcalJetSpectra::GetCentBin(Double_t cent) const 
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0;
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}

//________________________________________________________________________

Float_t AliAnalysisTaskEmcalJetSpectra:: RelativePhi(Double_t mphi,Double_t vphi) const
{
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;//dphi in [-Pi, Pi]                                                                                                    
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSpectra::Run()
{
  Int_t centbin = GetCentBin(fCent);
  //for pp analyses we will just use the first centrality bin
  if (centbin == -1)
    centbin = 0;

  if (!fTracks)
    return kTRUE;
  
  const Int_t nTrack = fTracks->GetEntriesFast();
  for (int i = 0;i<nTrack;i++){
    AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
    if (! track)
      continue;
    fHistTrackPt[centbin]->Fill(track->Pt());
  }
    
    if(!fLocalRho) {
        AliWarning(Form("%s: No LocalRho object found, attempting to get it from Event based on name!",GetName()));
        fLocalRho = GetLocalRhoFromEvent(fLocalRhoName);
    }

  fHistEP0[centbin]->Fill(fEPV0);
  fHistEP0A[centbin]->Fill(fEPV0A);
  fHistEP0C[centbin]->Fill(fEPV0C);
  fHistEPAvsC[centbin]->Fill(fEPV0A,fEPV0C);
  fRho = GetRhoFromEvent(fRhoName);
  fRhoVal = fRho->GetVal();
  fHistRhovsCent->Fill(fCent,fRhoVal);
  fHistRhovsEP[centbin]->Fill(fRhoVal,fEPV0);
  const Int_t Njets = fJets->GetEntriesFast();

  Int_t NjetAcc = 0;
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
     AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
     if (!jet)
       continue; 
     if (jet->Area()==0)
       continue;
     if (jet->Pt()<0.1)
       continue;
     if (jet->MaxTrackPt()>100)
       continue;
     if (! AcceptJet(jet))
       continue;
     //  jets.push_back(jet);
     NjetAcc++;
     Double_t jetPt = -500;
     jetPt = jet->Pt()-jet->Area()*fRhoVal;    
     fHistJetPtvsTrackPt[centbin]->Fill(jetPt,jet->MaxTrackPt());
     fHistRawJetPtvsTrackPt[centbin]->Fill(jet->Pt(),jet->MaxTrackPt());
     fHistJetPtvsdEP[centbin]->Fill(jetPt,RelativePhi((fEPV0+TMath::Pi()),jet->Phi()));
     fHistJetPtvsEP[centbin]->Fill(jetPt,fEPV0);

     // get local rho value
     fLocalRhoVal = fLocalRho->GetLocalVal(jet->Phi(), 0.2);
     Double_t jetPtLocal = jet->Pt() - jet->Area()*fLocalRhoVal;
     Double_t jetPtGlobal = jet->Pt() - jet->Area()*fRhoVal;

     // calculate relative angle between jet and event plane       
     Float_t dEP = -500;              // initialize angle between jet and event plane
     dEP = RelativeEPJET(jet->Phi(),fEPV0);

     fHistCorJetPtfromLocalRho[centbin]->Fill(jetPtLocal);
     fHistCorJetPtfromGlobalRho[centbin]->Fill(jetPtGlobal);

     fHistRhodEPcentLOC[centbin]->Fill(fLocalRhoVal,dEP);
     fHistRhodEPcentGL[centbin]->Fill(fRhoVal,dEP);
     fHistRhovsdEPLOC->Fill(fLocalRhoVal,dEP);
     fHistRhovsdEPGL->Fill(fRhoVal,dEP);

     fHistRhoEPcentLOC[centbin]->Fill(fLocalRhoVal,fEPV0);
     fHistRhoEPcentGL[centbin]->Fill(fRhoVal,fEPV0);
     fHistRhovsEPLOC->Fill(fLocalRhoVal,fEPV0);
     fHistRhovsEPGL->Fill(fRhoVal,fEPV0);

     fHistCorJetPtdEPcentLOC[centbin]->Fill(jetPtLocal,dEP);
     fHistCorJetPtdEPcentGL[centbin]->Fill(jetPt,dEP);
     fHistJetPtvsdEPLOC->Fill(jetPtLocal,dEP);
     fHistJetPtvsdEPGL->Fill(jetPt,dEP);

     fHistCorJetPtEPcentLOC[centbin]->Fill(jetPtLocal,fEPV0);
     fHistCorJetPtEPcentGL[centbin]->Fill(jetPt,fEPV0);
     fHistJetPtvsEPLOC->Fill(jetPtLocal,fEPV0);
     fHistJetPtvsEPGL->Fill(jetPt,fEPV0);

     fHistCorJetPt->Fill(jetPtLocal);
     fHistCorJetPtGL->Fill(jetPt);

     fHistGLvsLOCrho->Fill(fRhoVal,fLocalRhoVal);

     // apply max track bias
     if (jet->MaxTrackPt()>5.0){
       fHistJetPtvsdEPBias[centbin]->Fill(jetPt,RelativePhi((fEPV0+TMath::Pi()),jet->Phi()));
       fHistJetPtvsEPBias[centbin]->Fill(jetPt,fEPV0);
     }
  
     // in plane and out of plane histo's
     if( dEP>0 && dEP<=(TMath::Pi()/6) ){
        // we are IN plane, lets fill some histo's
        fHistCorJetPtfromLocalRhoIN[centbin]->Fill(jetPtLocal);
        fHistCorJetPtfromGlobalRhoIN[centbin]->Fill(jetPt);
     }else if( dEP>(TMath::Pi()/3) && dEP<=(TMath::Pi()/2) ){
        // we are OUT of PLANE, lets fill some histo's
        fHistCorJetPtfromLocalRhoOUT[centbin]->Fill(jetPtLocal);
        fHistCorJetPtfromGlobalRhoOUT[centbin]->Fill(jetPt);
     }

  }

  fHistNjetvsCent->Fill(fCent,NjetAcc);
  return kTRUE;
}      

//_________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetSpectra:: RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{ // function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi = (EPAng - jetAng);
    
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  }

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  }

  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 )
    AliWarning(Form("%s: dPHI is outside of restricted range: [0,Pi/2]",GetName()));

  return dphi;   // dphi in [0, Pi/2]
}




