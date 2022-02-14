/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

/// \file AliCSPairAnalysis.cxx
/// \brief Implementation of the AliCSPairAnalysis class

#include <TList.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TVector3.h>
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliCSPairAnalysis.h"

#define CSPAIRANALYSISNOOFTRACKS 3000 ///< the initial number of tracks

const Int_t AliCSPairAnalysis::kgHistosDimension = 4;

/// Default constructor for object serialization
AliCSPairAnalysis::AliCSPairAnalysis() :
    TNamed(),
    fOutput(NULL),

    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT1 bins */
    fNBins_pt(18),
    fMin_pt(0.2),
    fMax_pt(2.0),
    fWidth_pt(0.1),
    /* phi1 bins */
    fNBins_phi(72),
    fMin_phi(0.0),
    fMax_phi(TMath::Pi()*2.0),
    fWidth_phi(TMath::Pi()*2.0/72.0),
    /* eta1 bins */
    fNBins_eta(20),
    fMin_eta(-1.0),
    fMax_eta(1.0),
    fWidth_eta(0.1),

    /* the arrays with positive and negative tracks references */
    fPlusTracksArray(NULL),
    fMinusTracksArray(NULL),
    /* the singles efficiency structNures */
    fPlusEfficiency(NULL),
    fMinusEfficiency(NULL),
    fEventPlusEfficiency(NULL),
    fEventMinusEfficiency(NULL),
    /* the pairs efficiencies */
    fPairEfficiency_PP(NULL),
    fPairEfficiency_PM(NULL),
    fPairEfficiency_MM(NULL),
    fPairEfficiency_MP(NULL),
    /* histograms */
    fhNplus(NULL),
    fhNminus(NULL),
    fhPPDeltaEtaDeltaPhi(NULL),
    fhPMDeltaEtaDeltaPhi(NULL),
    fhMPDeltaEtaDeltaPhi(NULL),
    fhMMDeltaEtaDeltaPhi(NULL)
{

}

/// Normal constructor
/// \param name the name for the object instance
AliCSPairAnalysis::AliCSPairAnalysis(const char *name) :
    TNamed(name,name),
    fOutput(NULL),

    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT1 bins */
    fNBins_pt(18),
    fMin_pt(0.2),
    fMax_pt(2.0),
    fWidth_pt(0.1),
    /* phi1 bins */
    fNBins_phi(72),
    fMin_phi(0.0),
    fMax_phi(TMath::Pi()*2.0),
    fWidth_phi(TMath::Pi()*2.0/72.0),
    /* eta1 bins */
    fNBins_eta(20),
    fMin_eta(-1.0),
    fMax_eta(1.0),
    fWidth_eta(0.1),

    /* the arrays with positive and negative tracks references */
    fPlusTracksArray(NULL),
    fMinusTracksArray(NULL),
    /* the singles efficiency structures */
    fPlusEfficiency(NULL),
    fMinusEfficiency(NULL),
    fEventPlusEfficiency(NULL),
    fEventMinusEfficiency(NULL),
    /* the pairs efficiencies */
    fPairEfficiency_PP(NULL),
    fPairEfficiency_PM(NULL),
    fPairEfficiency_MM(NULL),
    fPairEfficiency_MP(NULL),
    /* histograms */
    fhNplus(NULL),
    fhNminus(NULL),
    fhPPDeltaEtaDeltaPhi(NULL),
    fhPMDeltaEtaDeltaPhi(NULL),
    fhMPDeltaEtaDeltaPhi(NULL),
    fhMMDeltaEtaDeltaPhi(NULL)
{

}

/// \brief Default destructor
/// Deallocates the allocated memory
AliCSPairAnalysis::~AliCSPairAnalysis() {
  if (fPlusTracksArray != NULL) delete fPlusTracksArray;
  if (fMinusTracksArray != NULL) delete fMinusTracksArray;
  if (fEventPlusEfficiency != NULL) {
    for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
      for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
        delete [] fEventPlusEfficiency[ieta][iphi];
      }
      delete [] fEventPlusEfficiency[ieta];
    }
    delete [] fEventPlusEfficiency;
  }
  if (fEventMinusEfficiency != NULL) {
    for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
      for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
        delete [] fEventMinusEfficiency[ieta][iphi];
      }
      delete [] fEventMinusEfficiency[ieta];
    }
    delete [] fEventMinusEfficiency;
  }
}


/// \brief Establishes the binning configuration
/// \param confstring string containing the binning configuration parameters
Bool_t AliCSPairAnalysis::ConfigureBinning(const char *confstring) {

  Double_t min_pt, max_pt, width_pt;
  Double_t min_eta, max_eta, width_eta;
  Int_t    nBins_phi = 0;
  Char_t buffer[24];

  /* we are a bit lazy here because the format should have been checked somewhere else */
  sscanf(confstring, "halfsymm:%3s;phishift:%lf;%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
      buffer,
      &fNBinsPhiShift,
      &fMin_vertexZ, &fMax_vertexZ, &fWidth_vertexZ,
      &min_pt, &max_pt, &width_pt,
      &min_eta, &max_eta, &width_eta, &nBins_phi);

  fMin_pt = min_pt;
  fMax_pt = max_pt;
  fWidth_pt = width_pt;
  fMin_eta = min_eta;
  fMax_eta = max_eta;
  fWidth_eta = width_eta;
  fNBins_phi = nBins_phi;
  AliInfo("=====================================================");
  AliInfo(Form("Configured binning: %s", GetBinningConfigurationString().Data()));
  AliInfo("=====================================================");
  return kTRUE;
}

/// \brief Build the configuration string
/// \return the configuration string corresponding to the current configuration
TString AliCSPairAnalysis::GetBinningConfigurationString() const {

  return TString(Form("Binning:phishift:%.1f;%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
      fNBinsPhiShift,
      fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
      fMin_pt, fMax_pt, fWidth_pt,
      fMin_eta, fMax_eta, fWidth_eta,fNBins_phi));
}



/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void AliCSPairAnalysis::Initialize()
{
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fPlusTracksArray = new TObjArray(CSPAIRANALYSISNOOFTRACKS);
  fPlusTracksArray->SetOwner(kFALSE);
  fMinusTracksArray = new TObjArray(CSPAIRANALYSISNOOFTRACKS);
  fMinusTracksArray->SetOwner(kFALSE);

  fEventPlusEfficiency = new Float_t**[fNBins_eta];
  fEventMinusEfficiency = new Float_t**[fNBins_eta];
  for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
    fEventPlusEfficiency[ieta] = new Float_t*[fNBins_phi];
    fEventMinusEfficiency[ieta] = new Float_t*[fNBins_phi];
    for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
      fEventPlusEfficiency[ieta][iphi] = new Float_t[fNBins_pt];
      fEventMinusEfficiency[ieta][iphi] = new Float_t[fNBins_pt];
    }
  }

  fOutput = new TList();
  fOutput->SetName(GetName());
  fOutput->SetOwner();

  fNBins_vertexZ     = Int_t(0.5 + (fMax_vertexZ - fMin_vertexZ) / fWidth_vertexZ);
  fNBins_pt          = Int_t(0.5 + (fMax_pt - fMin_pt ) / fWidth_pt);
  fNBins_eta         = Int_t(0.5 + (fMax_eta - fMin_eta) / fWidth_eta);
  fWidth_phi         = (fMax_phi  - fMin_phi) / fNBins_phi;

  fhNplus = new TH1F("Nplus", "(+) multiplicity; tracks", 400, 0, CSPAIRANALYSISNOOFTRACKS);
  fhNminus = new TH1F("Nminus", "(-) multiplicity; tracks", 400, 0, CSPAIRANALYSISNOOFTRACKS);

  Int_t nbinsaxes[kgHistosDimension] = {2*fNBins_eta-1,fNBins_phi,fNBins_pt,fNBins_pt};
  Double_t axeslow[kgHistosDimension] = {fMin_eta-fMax_eta,fMin_phi,fMin_pt,fMin_pt};
  Double_t axesup[kgHistosDimension] = {fMax_eta-fMin_eta,fMax_phi,fMax_pt,fMax_pt};

  /* we adapt the values of the phi limits after taking the delta phi values */
  fMax_phi           = fMax_phi - fWidth_phi * fNBinsPhiShift;
  fMin_phi           = fMin_phi - fWidth_phi * fNBinsPhiShift;

  fhPPDeltaEtaDeltaPhi = new THnF("PPDeltaEtaDeltaPhi","(+ +) #Delta #eta #Delta #varphi separation; #Delta #eta_{1};#Delta #varphi_{2};P_{T_{1}} (GeV/c)",
      kgHistosDimension, nbinsaxes,axeslow,axesup);
  fhPPDeltaEtaDeltaPhi->GetAxis(3)->SetTitle("P_{T_{2}} (GeV/c)");
  fhPMDeltaEtaDeltaPhi = new THnF("PMDeltaEtaDeltaPhi","(+ -) #Delta #eta #Delta #varphi separation; #Delta #eta_{1};#Delta #varphi_{2};P_{T_{1}} (GeV/c)",
      kgHistosDimension, nbinsaxes,axeslow,axesup);
  fhPMDeltaEtaDeltaPhi->GetAxis(3)->SetTitle("P_{T_{2}} (GeV/c)");
  fhMPDeltaEtaDeltaPhi = new THnF("MPDeltaEtaDeltaPhi","(- +) #Delta #eta #Delta #varphi separation; #Delta #eta_{1};#Delta #varphi_{2};P_{T_{1}} (GeV/c)",
      kgHistosDimension, nbinsaxes,axeslow,axesup);
  fhMPDeltaEtaDeltaPhi->GetAxis(3)->SetTitle("P_{T_{2}} (GeV/c)");
  fhMMDeltaEtaDeltaPhi = new THnF("MMDeltaEtaDeltaPhi","(- -) #Delta #eta #Delta #varphi separation; #Delta #eta_{1};#Delta #varphi_{2};P_{T_{1}} (GeV/c)",
      kgHistosDimension, nbinsaxes,axeslow,axesup);
  fhMMDeltaEtaDeltaPhi->GetAxis(3)->SetTitle("P_{T_{2}} (GeV/c)");

  fOutput->Add(fhNplus);
  fOutput->Add(fhNminus);
  fOutput->Add(fhPPDeltaEtaDeltaPhi);
  fOutput->Add(fhPMDeltaEtaDeltaPhi);
  fOutput->Add(fhMPDeltaEtaDeltaPhi);
  fOutput->Add(fhMMDeltaEtaDeltaPhi);

  TH1::AddDirectory(oldstatus);
}

/// \brief initializes the object instance to start a new event
Bool_t AliCSPairAnalysis::StartEvent(Float_t zVertex) {

  /* reset the tracks arrays */
  fPlusTracksArray->Clear();
  fMinusTracksArray->Clear();

  if (fPlusEfficiency != NULL) {
    /* initialize event efficiency from the efficiency histogram */
    Int_t zvtxbin = fPlusEfficiency->GetXaxis()->FindBin(zVertex);
    for (Int_t ipt = 0; ipt < fNBins_pt; ipt++) {
      for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
        Int_t etaphibinbase = ieta*fNBins_phi;
        for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
          fEventPlusEfficiency[ieta][iphi][ipt] = fPlusEfficiency->GetBinContent(zvtxbin,etaphibinbase+iphi+1,ipt+1);
        }
      }
    }
  }
  else {
    /* initialize the efficiency to one */
    for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
      for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
        for (Int_t ipt = 0; ipt < fNBins_pt; ipt++) {
          fEventPlusEfficiency[ieta][iphi][ipt] = 1.0;
        }
      }
    }
  }

  if (fMinusEfficiency != NULL) {
    /* initialize event efficiency from the efficiency histogram */
    Int_t zvtxbin = fMinusEfficiency->GetXaxis()->FindBin(zVertex);
    for (Int_t ipt = 0; ipt < fNBins_pt; ipt++) {
      for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
        Int_t etaphibinbase = ieta*fNBins_phi;
        for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
          fEventMinusEfficiency[ieta][iphi][ipt] = fMinusEfficiency->GetBinContent(zvtxbin,etaphibinbase+iphi+1,ipt+1);
        }
      }
    }
  }
  else {
    /* initialize the efficiency to one */
    for (Int_t ieta = 0; ieta < fNBins_eta; ieta++) {
      for (Int_t iphi = 0; iphi < fNBins_phi; iphi++) {
        for (Int_t ipt = 0; ipt < fNBins_pt; ipt++) {
          fEventMinusEfficiency[ieta][iphi][ipt] = 1.0;
        }
      }
    }
  }

  return kTRUE;
}

/// \brief Stores track reference in the corresponding array according to its charge
/// \param trk the intended track
/// \return kTRUE if the track is properly stored kFALSE otherwise
Bool_t AliCSPairAnalysis::ProcessTrack(Int_t, AliVTrack *trk) {

  if (trk->Charge() == 0) return kFALSE;

  /* add the track to the corresponding array */
  if (trk->Pt() < fMin_pt || fMax_pt < trk->Pt()) return kFALSE;
  if (trk->Eta() < fMin_eta || fMax_eta < trk->Eta()) return kFALSE;
  if (trk->Charge() > 0) fPlusTracksArray->Add(trk);
  else fMinusTracksArray->Add(trk);

  return kTRUE;
}

/// \brief Stores particle reference in the corresponding array according to its charge
/// \param par the intended particle
/// \return kTRUE if the track is properly stored kFALSE otherwise
Bool_t AliCSPairAnalysis::ProcessTrack(Int_t, AliVParticle *par) {

  if (par->Charge() == 0) return kFALSE;

  /* add the track to the corresponding array */
  if (par->Pt() < fMin_pt || fMax_pt < par->Pt()) return kFALSE;
  if (par->Eta() < fMin_eta || fMax_eta < par->Eta()) return kFALSE;
  if (par->Charge() > 0) fPlusTracksArray->Add(par);
  else fMinusTracksArray->Add(par);

  return kTRUE;
}

/// \brief Perform the pairs processing for the whole event

void AliCSPairAnalysis::ProcessEventData() {

  fhNplus->Fill(fPlusTracksArray->GetEntriesFast());
  fhNminus->Fill(fMinusTracksArray->GetEntriesFast());

  Bool_t trackdata = kTRUE;
  if (fPlusTracksArray->GetEntriesFast() != 0 || fMinusTracksArray->GetEntriesFast() != 0) {
    /* track data or particle data */
    if (fPlusTracksArray->GetEntriesFast() != 0)
      trackdata = fPlusTracksArray->At(0)->IsA()->InheritsFrom("AliVTrack");
    if (fMinusTracksArray->GetEntriesFast() != 0)
      trackdata = fMinusTracksArray->At(0)->IsA()->InheritsFrom("AliVTrack");
  }
  else {
    /* no tracks stored */
    return;
  }

  if (trackdata) {
    for (Int_t ixplus1 = 0; ixplus1 < fPlusTracksArray->GetEntriesFast(); ixplus1++) {
      Int_t fillbins[kgHistosDimension];
      Int_t fillbinssw[kgHistosDimension];

      /* the plus-plus pairs analysis */
      AliVTrack *trk1 = (AliVTrack *) fPlusTracksArray->At(ixplus1);
      Int_t ixeta1 = Int_t ((trk1->Eta() - fMin_eta) / fWidth_eta);
      Int_t ixphi1 = Int_t ((((trk1->Phi() < fMax_phi) ? trk1->Phi() : trk1->Phi() - 2*TMath::Pi())- fMin_phi) / fWidth_phi);
      Int_t ixpt1 = fhPPDeltaEtaDeltaPhi->GetAxis(2)->FindBin(trk1->Pt()) - 1;

      /* start correcting for the efficiency on the first track, positive */
      Float_t fillweight1 = fEventPlusEfficiency[ixeta1][ixphi1][ixpt1];

      fillbins[2] = fillbinssw[3] = ixpt1+1;

      for (Int_t ixplus2 = ixplus1+1; ixplus2 < fPlusTracksArray->GetEntriesFast(); ixplus2++) {
        AliVTrack *trk2 = (AliVTrack *) fPlusTracksArray->At(ixplus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhPPDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        /* correct for the efficiency on the second track, positive */
        Float_t fillweight = fillweight1 * fEventPlusEfficiency[ixeta2][ixphi2][ixpt2];

        fillbins[3] = fillbinssw[2] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdetasw = ixeta2 - ixeta1 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;
        Int_t ixdphisw = ixphi2 - ixphi1; if (ixdphisw < 0) ixdphisw += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;
        fillbinssw[0] = ixdetasw+1;
        fillbinssw[1] = ixdphisw+1;

        /* correct for the pair efficiency if applicable */
        if (fPairEfficiency_PP != 0) {
          fillweight /= fPairEfficiency_PP->GetBinContent(fillbins);
        }

        fhPPDeltaEtaDeltaPhi->AddBinContent(fillbins,fillweight);
        fhPPDeltaEtaDeltaPhi->AddBinContent(fillbinssw,fillweight);
      }
      /* the plus-minus pairs analysis */
      for (Int_t ixminus2 = 0; ixminus2 < fMinusTracksArray->GetEntriesFast(); ixminus2++) {
        AliVTrack *trk2 = (AliVTrack *) fMinusTracksArray->At(ixminus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhPMDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        /* correct for the efficiency on the second track, negative */
        Float_t fillweight = fillweight1 * fEventMinusEfficiency[ixeta2][ixphi2][ixpt2];

        fillbins[3] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;

        /* correct for the pair efficiency if applicable */
        if (fPairEfficiency_PM != 0) {
          fillweight /= fPairEfficiency_PM->GetBinContent(fillbins);
        }

        fhPMDeltaEtaDeltaPhi->AddBinContent(fillbins,fillweight);
      }
    }


    for (Int_t ixminus1 = 0; ixminus1 < fMinusTracksArray->GetEntriesFast(); ixminus1++) {
      Int_t fillbins[kgHistosDimension];
      Int_t fillbinssw[kgHistosDimension];

      /* the minus-minus pairs analysis */
      AliVTrack *trk1 = (AliVTrack *) fMinusTracksArray->At(ixminus1);
      Int_t ixeta1 = Int_t ((trk1->Eta() - fMin_eta) / fWidth_eta);
      Int_t ixphi1 = Int_t ((((trk1->Phi() < fMax_phi) ? trk1->Phi() : trk1->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
      Int_t ixpt1 = fhMMDeltaEtaDeltaPhi->GetAxis(2)->FindBin(trk1->Pt()) - 1;

      /* start correcting for the efficiency on the first track, negative */
      Float_t fillweight1 = fEventMinusEfficiency[ixeta1][ixphi1][ixpt1];

      fillbins[2] = fillbinssw[3] = ixpt1+1;

      for (Int_t ixminus2 = ixminus1+1; ixminus2 < fMinusTracksArray->GetEntriesFast(); ixminus2++) {
        AliVTrack *trk2 = (AliVTrack *) fMinusTracksArray->At(ixminus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhMMDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        /* correct for the efficiency on the second track, negative */
        Float_t fillweight = fillweight1 * fEventMinusEfficiency[ixeta2][ixphi2][ixpt2];

        fillbins[3] = fillbinssw[2] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdetasw = ixeta2 - ixeta1 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;
        Int_t ixdphisw = ixphi2 - ixphi1; if (ixdphisw < 0) ixdphisw += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;
        fillbinssw[0] = ixdetasw+1;
        fillbinssw[1] = ixdphisw+1;

        /* correct for the pair efficiency if applicable */
        if (fPairEfficiency_MM != 0) {
          fillweight /= fPairEfficiency_MM->GetBinContent(fillbins);
        }

        fhMMDeltaEtaDeltaPhi->AddBinContent(fillbins,fillweight);
        fhMMDeltaEtaDeltaPhi->AddBinContent(fillbinssw,fillweight);
      }

      /* the minus-plus pairs analysis */
      for (Int_t ixplus2 = 0; ixplus2 < fPlusTracksArray->GetEntriesFast(); ixplus2++) {
        AliVTrack *trk2 = (AliVTrack *) fPlusTracksArray->At(ixplus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhMPDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        /* correct for the efficiency on the second track, positive */
        Float_t fillweight = fillweight1 * fEventPlusEfficiency[ixeta2][ixphi2][ixpt2];

        fillbins[3] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;

        /* correct for the pair efficiency if applicable */
        if (fPairEfficiency_MP != 0) {
          fillweight /= fPairEfficiency_MP->GetBinContent(fillbins);
        }

        fhMPDeltaEtaDeltaPhi->AddBinContent(fillbins,fillweight);
      }
    }
  }
  else {
    /* for ESD MC processing */
    for (Int_t ixplus1 = 0; ixplus1 < fPlusTracksArray->GetEntriesFast(); ixplus1++) {
      Int_t fillbins[kgHistosDimension];
      Int_t fillbinssw[kgHistosDimension];

      /* the plus-plus pairs analysis */
      AliVParticle *trk1 = (AliVParticle *) fPlusTracksArray->At(ixplus1);
      Int_t ixeta1 = Int_t ((trk1->Eta() - fMin_eta) / fWidth_eta);
      Int_t ixphi1 = Int_t ((((trk1->Phi() < fMax_phi) ? trk1->Phi() : trk1->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
      Int_t ixpt1 = fhPPDeltaEtaDeltaPhi->GetAxis(2)->FindBin(trk1->Pt()) - 1;

      fillbins[2] = fillbinssw[3] = ixpt1+1;

      for (Int_t ixplus2 = ixplus1+1; ixplus2 < fPlusTracksArray->GetEntriesFast(); ixplus2++) {
        AliVParticle *trk2 = (AliVParticle *) fPlusTracksArray->At(ixplus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhPPDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        fillbins[3] = fillbinssw[2] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdetasw = ixeta2 - ixeta1 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;
        Int_t ixdphisw = ixphi2 - ixphi1; if (ixdphisw < 0) ixdphisw += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;
        fillbinssw[0] = ixdetasw+1;
        fillbinssw[1] = ixdphisw+1;

        fhPPDeltaEtaDeltaPhi->AddBinContent(fillbins,1.0);
        fhPPDeltaEtaDeltaPhi->AddBinContent(fillbinssw,1.0);
      }
      /* the plus-minus pairs analysis */
      for (Int_t ixminus2 = 0; ixminus2 < fMinusTracksArray->GetEntriesFast(); ixminus2++) {
        AliVParticle *trk2 = (AliVParticle *) fMinusTracksArray->At(ixminus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhPMDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        fillbins[3] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;

        fhPMDeltaEtaDeltaPhi->AddBinContent(fillbins,1.0);
      }
    }

    for (Int_t ixminus1 = 0; ixminus1 < fMinusTracksArray->GetEntriesFast(); ixminus1++) {
      Int_t fillbins[kgHistosDimension];
      Int_t fillbinssw[kgHistosDimension];

      /* the minus-minus pairs analysis */
      AliVParticle *trk1 = (AliVParticle *) fMinusTracksArray->At(ixminus1);
      Int_t ixeta1 = Int_t ((trk1->Eta() - fMin_eta) / fWidth_eta);
      Int_t ixphi1 = Int_t ((((trk1->Phi() < fMax_phi) ? trk1->Phi() : trk1->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
      Int_t ixpt1 = fhMMDeltaEtaDeltaPhi->GetAxis(2)->FindBin(trk1->Pt()) - 1;

      fillbins[2] = fillbinssw[3] = ixpt1+1;

      for (Int_t ixminus2 = ixminus1+1; ixminus2 < fMinusTracksArray->GetEntriesFast(); ixminus2++) {
        AliVParticle *trk2 = (AliVParticle *) fMinusTracksArray->At(ixminus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhMMDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        fillbins[3] = fillbinssw[2] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdetasw = ixeta2 - ixeta1 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;
        Int_t ixdphisw = ixphi2 - ixphi1; if (ixdphisw < 0) ixdphisw += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;
        fillbinssw[0] = ixdetasw+1;
        fillbinssw[1] = ixdphisw+1;

        fhMMDeltaEtaDeltaPhi->AddBinContent(fillbins,1.0);
        fhMMDeltaEtaDeltaPhi->AddBinContent(fillbinssw,1.0);
      }
      /* the minus-plus pairs analysis */
      for (Int_t ixplus2 = 0; ixplus2 < fPlusTracksArray->GetEntriesFast(); ixplus2++) {
        AliVParticle *trk2 = (AliVParticle *) fPlusTracksArray->At(ixplus2);
        Int_t ixeta2 = Int_t ((trk2->Eta() - fMin_eta) / fWidth_eta);
        Int_t ixphi2 = Int_t ((((trk2->Phi() < fMax_phi) ? trk2->Phi() : trk2->Phi() - 2*TMath::Pi() )- fMin_phi) / fWidth_phi);
        Int_t ixpt2 = fhMPDeltaEtaDeltaPhi->GetAxis(3)->FindBin(trk2->Pt()) - 1;

        fillbins[3] = ixpt2+1;

        Int_t ixdeta = ixeta1 - ixeta2 + fNBins_eta - 1;
        Int_t ixdphi = ixphi1 - ixphi2; if (ixdphi < 0) ixdphi += fNBins_phi;

        fillbins[0] = ixdeta+1;
        fillbins[1] = ixdphi+1;

        fhMPDeltaEtaDeltaPhi->AddBinContent(fillbins,1.0);
      }
    }
  }
}

/// \brief Do the final process for the end of analysis
void  AliCSPairAnalysis::FinalizeProcess()
{
}


/// \cond CLASSIMP
ClassImp(AliCSPairAnalysis);
/// \endcond
