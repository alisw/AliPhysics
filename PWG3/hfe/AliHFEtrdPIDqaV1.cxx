/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */

//
// Class AliHFEtrdPIDqaV1
// Monitoring TRD PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   TRD electron likelihood
//   TRD dE/dx (Absolute values)
//   TPC dE/dx (Number of sigmas, control histogram)
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//    Markus Fasel <M.Fasel@gsi.de>
//

#include <TAxis.h>
#include <TBrowser.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"

#include "AliHFEcollection.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTRD.h"
#include "AliHFEtrdPIDqaV1.h"

ClassImp(AliHFEtrdPIDqaV1)

//____________________________________________________________
AliHFEtrdPIDqaV1::AliHFEtrdPIDqaV1():
    AliHFEdetPIDqa(),
    fHistos(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliHFEtrdPIDqaV1::AliHFEtrdPIDqaV1(const Char_t *name):
    AliHFEdetPIDqa(name, "QA for TRD"),
    fHistos(NULL)
{
  //
  // Default constructor
  //
}

//____________________________________________________________
AliHFEtrdPIDqaV1::AliHFEtrdPIDqaV1(const AliHFEtrdPIDqaV1 &o):
    AliHFEdetPIDqa(o),
    fHistos(NULL)
{
  //
  // Copy constructor
  //
}

//____________________________________________________________
AliHFEtrdPIDqaV1 &AliHFEtrdPIDqaV1::operator=(const AliHFEtrdPIDqaV1 &o){
  //
  // Make assignment
  //
  AliHFEdetPIDqa::operator=(o);
  fHistos = o.fHistos;
  
  return *this;
}

//_________________________________________________________
Long64_t AliHFEtrdPIDqaV1::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEtrdPIDqaV1 *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEtrdPIDqaV1 *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEtrdPIDqaV1::Browse(TBrowser *b){
  //
  // Browse the PID QA
  //
  if(b){
    if(fHistos){
      b->Add(fHistos, fHistos->GetName());

      // Make Projections of the dE/dx Spectra and add them to a new Folder
      TString specnames[4] = {"All", "Electrons", "Pions", "Protons"};
      Int_t specind[4] = {-1, AliPID::kElectron, AliPID::kPion, AliPID::kProton};
      TList *listTM = new TList;
      listTM->SetOwner();
      TList *listLike = new TList;
      listLike->SetOwner();
      TList *listCharge = new TList;
      listCharge->SetOwner();
      TList *listTPCnsigma = new TList;
      listTPCnsigma->SetOwner();

      TH2 *hptr = NULL; 
      for(Int_t ispec = 0; ispec < 4; ispec++){
        for(Int_t istep = 0; istep < 2; istep++){
          hptr = MakeTRDspectrumTM(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec]);
          if(hptr){
            hptr->SetName(Form("hTRDtm%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After"));
            listTM->Add(hptr);
          }
          hptr = MakeTRDlikelihoodDistribution(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec]);
          hptr->SetName(Form("hTRDlike%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After"));
          listLike->Add(hptr);
          hptr = MakeTRDchargeDistribution(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec]);
          hptr->SetName(Form("hTRDcharge%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After"));
          listCharge->Add(hptr);
          hptr = MakeTPCspectrumNsigma(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec]);
          hptr->SetName(Form("hTPCspectrum%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After"));
          listTPCnsigma->Add(hptr);
        }
      }
      
      b->Add(listTM, "Projections Truncated Mean");
      b->Add(listLike, "Projections Likelihood distribution");
      b->Add(listCharge, "Projections Tracklet Charge");
      b->Add(listTPCnsigma, "Projections TPC spectra");
    }
  }
}

//____________________________________________________________
void AliHFEtrdPIDqaV1::Initialize(){
  //
  // Initialize QA histos for TRD PID
  //

  AliDebug(1, "Initializing PID QA for TRD");
  // Make common binning
  const Int_t kPIDbins = AliPID::kSPECIES + 1;
  const Int_t kPbins = 100;
  const Int_t kSteps = 2;
  const Double_t kMinPID = -1;
  const Double_t kMinP = 0.;
  const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
  const Double_t kMaxP = 20.;

  fHistos = new AliHFEcollection("trdqahistos", "Collection of TRD QA histograms");
  
  // Create Control Histogram monitoring the TPC sigma between the selection steps
  const Int_t kTPCSigmaBins = 140;
  Int_t nBinsTPCSigma[4] = {kPIDbins, kPbins, kTPCSigmaBins, kSteps};
  Double_t minTPCSigma[4] = {kMinPID, kMinP, -12., 0};
  Double_t maxTPCSigma[4] = {kMaxPID, kMaxP, 12., 2.};
  fHistos->CreateTHnSparse("hTPCsigma", "TPC sigma; species p [GeV/c]; TPC dEdx - <dE/dx>|_{el} [#sigma]; selection step", 4, nBinsTPCSigma, minTPCSigma, maxTPCSigma);
  fHistos->Sumw2("hTPCsigma");
  // Create Monitoring histogram for the Likelihood distribution
  const Int_t kTRDLikelihoodBins = 100;
  Int_t nBinsTRDlike[4] = {kPIDbins, kPbins, kTRDLikelihoodBins, kSteps};
  Double_t minTRDlike[4] = {kMinPID, kMinP, 0., 0.};
  Double_t maxTRDlike[4] = {kMaxPID, kMaxP, 1., 2.};
  fHistos->CreateTHnSparse("hTRDlikelihood", "TRD Likelihood Distribution; species; p [GeV/c]; TRD electron Likelihood; selection step", 4, nBinsTRDlike, minTRDlike, maxTRDlike);
  fHistos->Sumw2("hTRDlikelihood");
  // Create Monitoring histogram for the TRD total charge
  const Int_t kTRDchargeBins = 1000;
  Int_t nBinsTRDcharge[4] = {kPIDbins, kPbins, kTRDchargeBins, kSteps};
  Double_t minTRDcharge[4] = {kMinPID, kMinP, 0., 0.};
  Double_t maxTRDcharge[4] = {kMaxPID, kMaxP, 100000., 2.};
  fHistos->CreateTHnSparse("hTRDcharge", "Total TRD charge; species; p [GeV/c]; TRD charge [a.u.]; selection step", 4, nBinsTRDcharge, minTRDcharge, maxTRDcharge);
  fHistos->Sumw2("hTRDcharge");
  // Monitoring of the TRD truncated mean according to version 1
  const Int_t kTRDtmBins = 1000;
  Int_t nBinsTRDtm[4] = {kPIDbins, kPbins, kTRDtmBins, kSteps};
  Double_t minTRDtm[4] = {kMinPID, kMinP, 0., 0.};
  Double_t maxTRDtm[4] = {kMaxPID, kMaxP, 20000., 2.};
  fHistos->CreateTHnSparse("hTRDtruncatedMean", "TRD truncated Mean; species; p [GeV/c]; TRD signal [a.u.]; selection step", 4, nBinsTRDtm, minTRDtm, maxTRDtm);
  fHistos->Sumw2("hTRDtruncatedMean");
}

//____________________________________________________________
void AliHFEtrdPIDqaV1::ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step){
  //
  // Process the track, fill the containers 
  //
  AliDebug(1, Form("QA started for TRD PID for step %d", (Int_t)step));
  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;

  AliHFEpidTRD *trdpid = dynamic_cast<AliHFEpidTRD *>(fQAmanager->GetDetectorPID(AliHFEpid::kTRDpid));
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid));
 
  Double_t container[4];
  container[0] = species;
  container[1] = trdpid ? trdpid->GetP(track->GetRecTrack(), anatype) : 0.;
  container[2] = tpcpid ? tpcpid->NumberOfSigmas(track->GetRecTrack(), AliPID::kElectron, anatype) : 0.;
  container[3] = step;
  fHistos->Fill("hTPCsigma", container);

  container[2] = trdpid->GetElectronLikelihood(track->GetRecTrack(), anatype);
  fHistos->Fill("hTRDlikelihood", container);

  if(track->IsESDanalysis()){
    container[2] = trdpid->GetTRDSignalV1(dynamic_cast<const AliESDtrack *>(track->GetRecTrack()));
    fHistos->Fill("hTRDtruncatedMean", container);
  }
  for(UInt_t ily = 0; ily < 6; ily++){
    container[2] = trdpid->GetChargeLayer(track->GetRecTrack(), ily, anatype);
    fHistos->Fill("hTRDcharge", container);
  }
}

//_________________________________________________________
TH2 *AliHFEtrdPIDqaV1::MakeTPCspectrumNsigma(AliHFEdetPIDqa::EStep_t step, Int_t species){
  //
  // Get the TPC control histogram for the TRD selection step (either before or after PID)
  //
  THnSparseF *histo = dynamic_cast<THnSparseF *>(fHistos->Get("hTPCsigma"));
  if(!histo){
    AliError("QA histogram monitoring TPC nSigma not available");
    return NULL;
  }
  if(species > -1 && species < AliPID::kSPECIES){
    // cut on species (if available)
    histo->GetAxis(0)->SetRange(species + 2, species + 2); // undef + underflow
  }
  histo->GetAxis(3)->SetRange(step + 1, step + 1); 

  TH2 *hSpec = histo->Projection(2, 1);
  // construct title and name
  TString stepname = step == AliHFEdetPIDqa::kBeforePID ? "before" : "after";
  TString speciesname = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "all Particles";
  TString specID = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "unid";
  TString histname = Form("hSigmaTPC%s%s", specID.Data(), stepname.Data());
  TString histtitle = Form("TPC Sigma for %s %s PID", speciesname.Data(), stepname.Data());
  hSpec->SetName(histname.Data());
  hSpec->SetTitle(histtitle.Data());

  // Unset range on the original histogram
  histo->GetAxis(0)->SetRange(0, histo->GetAxis(0)->GetNbins());
  histo->GetAxis(2)->SetRange(0, histo->GetAxis(2)->GetNbins());
  return hSpec; 
}

//_________________________________________________________
TH2 *AliHFEtrdPIDqaV1::MakeTRDspectrumTM(AliHFEdetPIDqa::EStep_t step, Int_t species){
  //
  // Get the TPC control histogram for the TRD selection step (either before or after PID)
  //
  THnSparseF *histo = dynamic_cast<THnSparseF *>(fHistos->Get("hTRDtruncatedMean"));
  if(!histo){
    AliError("QA histogram monitoring TPC nSigma not available");
    return NULL;
  }
  if(species > -1 && species < AliPID::kSPECIES){
    // cut on species (if available)
    histo->GetAxis(0)->SetRange(species + 2, species + 2); // undef + underflow
  }
  histo->GetAxis(3)->SetRange(step + 1, step + 1); 

  TH2 *hSpec = histo->Projection(2, 1);
  // construct title and name
  TString stepname = step == AliHFEdetPIDqa::kBeforePID ? "before" : "after";
  TString speciesname = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "all Particles";
  TString specID = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "unid";
  TString histname = Form("hTMTRD%s%s", specID.Data(), stepname.Data());
  TString histtitle = Form("TRD Truncated Mean for %s %s PID", speciesname.Data(), stepname.Data());
  hSpec->SetName(histname.Data());
  hSpec->SetTitle(histtitle.Data());

  // Unset range on the original histogram
  histo->GetAxis(0)->SetRange(0, histo->GetAxis(0)->GetNbins());
  histo->GetAxis(2)->SetRange(0, histo->GetAxis(2)->GetNbins());
  return hSpec; 
}

//_________________________________________________________
TH2 *AliHFEtrdPIDqaV1::MakeTRDlikelihoodDistribution(AliHFEdetPIDqa::EStep_t step, Int_t species){
  //
  // Make Histogram for TRD Likelihood distribution
  //
  THnSparseF *histo = dynamic_cast<THnSparseF *>(fHistos->Get("hTRDlikelihood"));
  if(!histo){
    AliError("QA histogram monitoring TRD Electron Likelihood not available");
    return NULL;
  }
  if(species > -1 && species < AliPID::kSPECIES){
    // cut on species (if available)
    histo->GetAxis(0)->SetRange(species + 2, species + 2); // undef + underflow
  }
  histo->GetAxis(3)->SetRangeUser(step + 1, step + 1);

  TH2 *hSpec = histo->Projection(2, 1);
  // construct title and name
  TString stepname = step == AliHFEdetPIDqa::kBeforePID ? "before" : "after";
  TString speciesname = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "all Particles";
  TString specID = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "unid";
  TString histname = Form("hLikeElTRD%s%s", specID.Data(), stepname.Data());
  TString histtitle = Form("TRD electron Likelihood for %s %s PID", speciesname.Data(), stepname.Data());
  hSpec->SetName(histname.Data());
  hSpec->SetTitle(histtitle.Data());

  // Unset range on the original histogram
  histo->GetAxis(0)->SetRange(0, histo->GetAxis(0)->GetNbins());
  histo->GetAxis(2)->SetRange(0, histo->GetAxis(2)->GetNbins());
  return hSpec; 
}

//_________________________________________________________
TH2 *AliHFEtrdPIDqaV1::MakeTRDchargeDistribution(AliHFEdetPIDqa::EStep_t step, Int_t species){
  //
  // Make Histogram for TRD Likelihood distribution
  //
  THnSparseF *histo = dynamic_cast<THnSparseF *>(fHistos->Get("hTRDcharge"));
  if(!histo){
    AliError("QA histogram monitoring TRD total charge not available");
    return NULL;
  }
  if(species > -1 && species < AliPID::kSPECIES){
    // cut on species (if available)
    histo->GetAxis(0)->SetRange(species + 2, species + 2); // undef + underflow
  }
  histo->GetAxis(3)->SetRange(step + 1, step + 1);

  TH2 *hSpec = histo->Projection(2, 1);
  // construct title and name
  TString stepname = step == AliHFEdetPIDqa::kBeforePID ? "before" : "after";
  TString speciesname = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "all Particles";
  TString specID = species > -1 && species < AliPID::kSPECIES ? AliPID::ParticleName(species) : "unid";
  TString histname = Form("hChargeTRD%s%s", specID.Data(), stepname.Data());
  TString histtitle = Form("TRD total charge for %s %s PID", speciesname.Data(), stepname.Data());
  hSpec->SetName(histname.Data());
  hSpec->SetTitle(histtitle.Data());

  // Unset range on the original histogram
  histo->GetAxis(0)->SetRange(0, histo->GetAxis(0)->GetNbins());
  histo->GetAxis(2)->SetRange(0, histo->GetAxis(2)->GetNbins());
  return hSpec; 
}

