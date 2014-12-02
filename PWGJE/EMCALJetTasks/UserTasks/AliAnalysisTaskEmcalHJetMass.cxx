//
// Jet mass analysis task for jets recoiling from high pT hadron
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TRandom3.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalHJetMass.h"

ClassImp(EmcalHJetMassAnalysis::AliAnalysisTaskEmcalHJetMass)

namespace EmcalHJetMassAnalysis {

  //________________________________________________________________________
  AliAnalysisTaskEmcalHJetMass::AliAnalysisTaskEmcalHJetMass() : 
    AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHJetMass", kTRUE),
    fContainerBase(0),
    fContainerUnsub(1),
    fMinFractionShared(0),
    fUseUnsubJet(0),
    fJetMassType(kRaw),
    fDPhiHJetMax(0.6),
    fTriggerTrackType(kInclusive),
    fPtTTMin(0),
    fPtTTMax(0),
    fRandom(0),
    fh1PtHadron(0),
    fh3PtHPtJDPhi(0),
    fh3PtJet1VsMassVsHPtAllSel(0),
    fh3PtJet1VsMassVsHPtAllSelMatch(0),
    fh3PtJet1VsMassVsHPtTagged(0),
    fh3PtJet1VsMassVsHPtTaggedMatch(0),
    fh3PtJet1VsRatVsHPtAllSel(0),
    fh3PtJet1VsRatVsHPtAllSelMatch(0),
    fh3PtJet1VsRatVsHPtTagged(0),
    fh3PtJet1VsRatVsHPtTaggedMatch(0)
  {
    // Default constructor.

    fh1PtHadron                       = new TH1F*[fNcentBins];
    fh3PtHPtJDPhi                     = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtAllSel        = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtAllSelMatch   = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtTagged        = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtTaggedMatch   = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtAllSel         = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtAllSelMatch    = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtTagged         = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtTaggedMatch    = new TH3F*[fNcentBins];

    for (Int_t i = 0; i < fNcentBins; i++) {
      fh1PtHadron[i]                       = 0;
      fh3PtHPtJDPhi[i]                     = 0;
      fh3PtJet1VsMassVsHPtAllSel[i]        = 0;
      fh3PtJet1VsMassVsHPtAllSelMatch[i]   = 0;
      fh3PtJet1VsMassVsHPtTagged[i]        = 0;
      fh3PtJet1VsMassVsHPtTaggedMatch[i]   = 0;
      fh3PtJet1VsRatVsHPtAllSel[i]         = 0;
      fh3PtJet1VsRatVsHPtAllSelMatch[i]    = 0;
      fh3PtJet1VsRatVsHPtTagged[i]         = 0;
      fh3PtJet1VsRatVsHPtTaggedMatch[i]    = 0;
    }

    fPtTTMin = new TArrayF();
    fPtTTMax = new TArrayF();

    SetMakeGeneralHistograms(kTRUE);
  }

  //________________________________________________________________________
  AliAnalysisTaskEmcalHJetMass::AliAnalysisTaskEmcalHJetMass(const char *name) : 
    AliAnalysisTaskEmcalJet(name, kTRUE),  
    fContainerBase(0),
    fContainerUnsub(1),
    fMinFractionShared(0),
    fUseUnsubJet(0),
    fJetMassType(kRaw),
    fDPhiHJetMax(0.6),
    fTriggerTrackType(kInclusive),
    fPtTTMin(0),
    fPtTTMax(0),
    fRandom(0),
    fh1PtHadron(0),
    fh3PtHPtJDPhi(0),
    fh3PtJet1VsMassVsHPtAllSel(0),
    fh3PtJet1VsMassVsHPtAllSelMatch(0),
    fh3PtJet1VsMassVsHPtTagged(0),
    fh3PtJet1VsMassVsHPtTaggedMatch(0),
    fh3PtJet1VsRatVsHPtAllSel(0),
    fh3PtJet1VsRatVsHPtAllSelMatch(0),
    fh3PtJet1VsRatVsHPtTagged(0),
    fh3PtJet1VsRatVsHPtTaggedMatch(0)
  {
    // Standard constructor.

    fh1PtHadron                       = new TH1F*[fNcentBins];
    fh3PtHPtJDPhi                     = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtAllSel        = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtAllSelMatch   = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtTagged        = new TH3F*[fNcentBins];
    fh3PtJet1VsMassVsHPtTaggedMatch   = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtAllSel         = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtAllSelMatch    = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtTagged         = new TH3F*[fNcentBins];
    fh3PtJet1VsRatVsHPtTaggedMatch    = new TH3F*[fNcentBins];
 
    for (Int_t i = 0; i < fNcentBins; i++) {
      fh1PtHadron[i]                       = 0;
      fh3PtHPtJDPhi[i]                     = 0;
      fh3PtJet1VsMassVsHPtAllSel[i]        = 0;
      fh3PtJet1VsMassVsHPtAllSelMatch[i]   = 0;
      fh3PtJet1VsMassVsHPtTagged[i]        = 0;
      fh3PtJet1VsMassVsHPtTaggedMatch[i]   = 0;
      fh3PtJet1VsRatVsHPtAllSel[i]         = 0;
      fh3PtJet1VsRatVsHPtAllSelMatch[i]    = 0;
      fh3PtJet1VsRatVsHPtTagged[i]         = 0;
      fh3PtJet1VsRatVsHPtTaggedMatch[i]    = 0;
    }

    fPtTTMin = new TArrayF();
    fPtTTMax = new TArrayF();

    SetMakeGeneralHistograms(kTRUE);
  }

  //________________________________________________________________________
  AliAnalysisTaskEmcalHJetMass::~AliAnalysisTaskEmcalHJetMass()
  {
    // Destructor.

    if(fRandom)      delete fRandom;
    if(fPtTTMin)     delete fPtTTMin;
    if(fPtTTMax)     delete fPtTTMax;
  }

  //________________________________________________________________________
  void AliAnalysisTaskEmcalHJetMass::UserCreateOutputObjects()
  {
    // Create user output.

    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    const Int_t nBinsPt  = 200;
    const Double_t minPt = -50.;
    const Double_t maxPt = 150.;

    const Int_t nBinsM  = 100;
    const Double_t minM = -20.;
    const Double_t maxM = 80.;

    const Int_t nBinsR  = 100;
    const Double_t minR = -0.2;
    const Double_t maxR = 0.8;

    const Int_t nBinsPtH  = 100;
    const Double_t minPtH = 0.;
    const Double_t maxPtH = 100.;

    const Int_t nBinsPhi  = 18*4;
    const Double_t minPhi = -0.5*TMath::Pi();
    const Double_t maxPhi = 1.5*TMath::Pi();

    TString histName = "";
    TString histTitle = "";
    for (Int_t i = 0; i < fNcentBins; i++) {
      histName = TString::Format("fh1PtHadron_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,h}",histName.Data());
      fh1PtHadron[i] = new TH1F(histName.Data(),histTitle.Data(),200.,0.,200.);
      fOutput->Add(fh1PtHadron[i]);

      histName = TString::Format("fh3PtHPtJDPhi_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,h};#it{p}_{T,jet};#Delta#varphi_{h,jet}",histName.Data());
      fh3PtHPtJDPhi[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPtH,minPtH,maxPtH,nBinsPt,minPt,maxPt,nBinsPhi,minPhi,maxPhi);
      fOutput->Add(fh3PtHPtJDPhi[i]);

      histName = TString::Format("fh3PtJet1VsMassVsHPtAllSel_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsMassVsHPtAllSel[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsMassVsHPtAllSel[i]);

      histName = TString::Format("fh3PtJet1VsMassVsHPtAllSelMatch_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsMassVsHPtAllSelMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsMassVsHPtAllSelMatch[i]);

      histName = TString::Format("fh3PtJet1VsMassVsHPtTagged_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsMassVsHPtTagged[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsMassVsHPtTagged[i]);

      histName = TString::Format("fh3PtJet1VsMassVsHPtTaggedMatch_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsMassVsHPtTaggedMatch[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsMassVsHPtTaggedMatch[i]);

      //
      histName = TString::Format("fh3PtJet1VsRatVsHPtAllSel_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsRatVsHPtAllSel[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsRatVsHPtAllSel[i]);

      histName = TString::Format("fh3PtJet1VsRatVsHPtAllSelMatch_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsRatVsHPtAllSelMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsRatVsHPtAllSelMatch[i]);

      histName = TString::Format("fh3PtJet1VsRatVsHPtTagged_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsRatVsHPtTagged[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsRatVsHPtTagged[i]);

      histName = TString::Format("fh3PtJet1VsRatVsHPtTaggedMatch_%d",i);
      histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,h}",histName.Data());
      fh3PtJet1VsRatVsHPtTaggedMatch[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,nBinsPtH,minPtH,maxPtH);
      fOutput->Add(fh3PtJet1VsRatVsHPtTaggedMatch[i]);
    }

    TH1::AddDirectory(oldStatus);

    fRandom = new TRandom3(0);

    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  }

  //________________________________________________________________________
  Bool_t AliAnalysisTaskEmcalHJetMass::Run()
  {
    // Run analysis code here, if needed. It will be executed before FillHistograms().

    AliParticleContainer *pCont = GetParticleContainer(0);
    AliJetContainer      *jCont = GetJetContainer(fContainerBase);
    if(!pCont || !jCont) return kFALSE;

    AliVParticle *vp = NULL;

    if(fTriggerTrackType==kInclusive) {
      pCont->ResetCurrentID();
      while((vp = pCont->GetNextAcceptParticle())) {
        fh1PtHadron[fCentBin]->Fill(vp->Pt()); //all hadrons
        AliEmcalJet* jet = NULL;
        if(jCont) {
          jCont->ResetCurrentID();
          while((jet = jCont->GetNextAcceptJet())) {
            Double_t dphi = GetDeltaPhi(vp,jet)-TMath::Pi();
            fh3PtHPtJDPhi[fCentBin]->Fill(vp->Pt(),jet->Pt() - GetRhoVal(fContainerBase)*jet->Area(),dphi);
            if(dphi>fDPhiHJetMax) continue;
            FillHJetHistograms(vp->Pt(),jet);
          }
        }
      }
    }
    else if(fTriggerTrackType==kSingleInclusive) {
     for(Int_t it = 0; it<fPtTTMin->GetSize(); it++) {
       vp = GetSingleInclusiveTT(pCont,fPtTTMin->At(it),fPtTTMax->At(it));
       if(!vp) continue;
       fh1PtHadron[fCentBin]->Fill(vp->Pt()); //all trigger tracks
       AliEmcalJet* jet = NULL;
       jCont->ResetCurrentID();
       while((jet = jCont->GetNextAcceptJet())) {
         Double_t dphi = GetDeltaPhi(vp,jet)-TMath::Pi();
         fh3PtHPtJDPhi[fCentBin]->Fill(vp->Pt(),jet->Pt() - GetRhoVal(fContainerBase)*jet->Area(),dphi);
         if(dphi>fDPhiHJetMax) continue;
         FillHJetHistograms(vp->Pt(),jet);
       }
     }//trigger track types
    }
    return kTRUE;
  }

  //________________________________________________________________________
  AliVParticle* AliAnalysisTaskEmcalHJetMass::GetSingleInclusiveTT(AliParticleContainer *pCont, Double_t ptmin, Double_t ptmax) const {
    AliVParticle *vp;
    TArrayI arr; arr.Set(pCont->GetNParticles());
    arr.Reset();
    Int_t counter = -1;
    pCont->ResetCurrentID();
    while((vp = pCont->GetNextAcceptParticle())) {
      if(vp->Pt()>=ptmin && vp->Pt()<ptmax ) {
        counter++;
        arr.SetAt(pCont->GetCurrentID(),counter);
      }
    }
    if(counter<0) return NULL;
    //select trigger track randomly
    fRandom->SetSeed(arr.At(0)); //random selection reproducible
    Double_t rnd = fRandom->Uniform() * counter;
    Int_t trigID = arr.At(TMath::FloorNint(rnd));
    vp = pCont->GetParticle(trigID);
    return vp;
  }

  //________________________________________________________________________
  Bool_t AliAnalysisTaskEmcalHJetMass::FillHJetHistograms(Double_t pt, const AliEmcalJet *jet)
  {
    // Fill hadron-jet histograms.
    Double_t ptJet = jet->Pt() - GetRhoVal(fContainerBase)*jet->Area();
    Double_t mJet = GetJetMass(jet);
    Double_t rat = -1.;
    if(ptJet<0. || ptJet>0.) rat = mJet/ptJet;

    fh3PtJet1VsMassVsHPtAllSel[fCentBin]->Fill(ptJet,mJet,pt);
    fh3PtJet1VsRatVsHPtAllSel[fCentBin]->Fill(ptJet,rat,pt);

    Double_t fraction = -1.;
    if(fUseUnsubJet) {
      AliEmcalJet *jetUS = NULL;
      AliJetContainer *jetContUS = GetJetContainer(fContainerUnsub);
      Int_t ifound = 0;
      Int_t ilab = -1;
      for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
        jetUS = jetContUS->GetJet(i);
        if(jetUS->GetLabel()==jet->GetLabel()) {
          ifound++;
          if(ifound==1) ilab = i;
        }
      }
      if(ifound>1) AliDebug(2,Form("Found %d partners",ifound));
      if(ifound==0) jetUS = 0x0;
      else          jetUS = jetContUS->GetJet(ilab);
      fraction = jetContUS->GetFractionSharedPt(jetUS);
    } else {
      AliJetContainer *jetCont = GetJetContainer(fContainerBase);
      fraction = jetCont->GetFractionSharedPt(jet);
    }
  
    if(fMinFractionShared>0. && fraction>fMinFractionShared) {
      fh3PtJet1VsMassVsHPtAllSelMatch[fCentBin]->Fill(ptJet,mJet,pt);
      fh3PtJet1VsRatVsHPtAllSelMatch[fCentBin]->Fill(ptJet,rat,pt);
    }

    if(jet->GetTagStatus()<1 || !jet->GetTaggedJet())
      return kFALSE;

    fh3PtJet1VsMassVsHPtTagged[fCentBin]->Fill(ptJet,mJet,pt);
    fh3PtJet1VsRatVsHPtTagged[fCentBin]->Fill(ptJet,rat,pt);

    if(fMinFractionShared>0. && fraction>fMinFractionShared) {
      fh3PtJet1VsMassVsHPtTaggedMatch[fCentBin]->Fill(ptJet,mJet,pt);
      fh3PtJet1VsRatVsHPtTaggedMatch[fCentBin]->Fill(ptJet,rat,pt);
    }
    return kTRUE;
  }

  //________________________________________________________________________
  Double_t AliAnalysisTaskEmcalHJetMass::GetJetMass(const AliEmcalJet *jet) const {
    //calc subtracted jet mass
    if(fJetMassType==kRaw)
      return jet->M();
    else if(fJetMassType==kDeriv)
      return jet->GetSecondOrderSubtracted();
  
    return 0;
  }

  //________________________________________________________________________
  Double_t AliAnalysisTaskEmcalHJetMass::GetDeltaPhi(const AliVParticle *vp, const AliEmcalJet* jet) const {
    // Calculate azimuthal angle between particle and jet. range:[-0.5\pi,1.5\pi]
    return GetDeltaPhi(vp->Phi(),jet->Phi());
  }

  //________________________________________________________________________
  Double_t AliAnalysisTaskEmcalHJetMass::GetDeltaPhi(Double_t phi1, Double_t phi2) const {
    // Calculate azimuthal angle between phi1 and phi2. range:[-0.5\pi,1.5\pi]
    Double_t dPhi = phi1-phi2;
    if(dPhi <-0.5*TMath::Pi())  dPhi += TMath::TwoPi();
    if(dPhi > 1.5*TMath::Pi())  dPhi -= TMath::TwoPi();
    return dPhi;
  }


  //________________________________________________________________________
  Bool_t AliAnalysisTaskEmcalHJetMass::RetrieveEventObjects() {
    //
    // retrieve event objects
    //
    if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
      return kFALSE;

    return kTRUE;
  }

  //________________________________________________________________________
  void AliAnalysisTaskEmcalHJetMass::AddTriggerTrackPtCuts(Float_t min, Float_t max) {
    if(!fPtTTMin) fPtTTMin = new TArrayF();
    if(!fPtTTMax) fPtTTMax = new TArrayF();
    Int_t newSize = fPtTTMin->GetSize()+1;
    fPtTTMin->Set(newSize);
    fPtTTMax->Set(newSize);
    fPtTTMin->AddAt(min,newSize-1);
    fPtTTMax->AddAt(max,newSize-1);
  }

  //_______________________________________________________________________
  void AliAnalysisTaskEmcalHJetMass::Terminate(Option_t *) 
  {
    // Called once at the end of the analysis.
  }
}

