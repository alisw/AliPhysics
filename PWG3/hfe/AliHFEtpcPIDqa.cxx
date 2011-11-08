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
//
// Class AliHFEtpcPIDqa
// Monitoring TPC PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   TPC dE/dx (Number of sigmas)
//   TPC dE/dx (Absolute values)
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//    Markus Fasel <M.Fasel@gsi.de>
//
#include <TBrowser.h>
#include <TClass.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliHFEcollection.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTPC.h"
#include "AliHFEtools.h"
#include "AliHFEtpcPIDqa.h"

ClassImp(AliHFEtpcPIDqa)

//_________________________________________________________
AliHFEtpcPIDqa::AliHFEtpcPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
  , fBrowseCentrality(-1)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEtpcPIDqa::AliHFEtpcPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for TPC")
  , fHistos(NULL)
  , fBrowseCentrality(-1)
{
  //
  // Default constructor
  //
}

//_________________________________________________________
AliHFEtpcPIDqa::AliHFEtpcPIDqa(const AliHFEtpcPIDqa &o):
    AliHFEdetPIDqa(o)
  , fHistos(NULL)
  , fBrowseCentrality(o.fBrowseCentrality)
{
  //
  // Copy constructor
  //
  o.Copy(*this);
}

//_________________________________________________________
AliHFEtpcPIDqa &AliHFEtpcPIDqa::operator=(const AliHFEtpcPIDqa &o){
  //
  // Do assignment
  //
  AliHFEdetPIDqa::operator=(o);
  if(&o != this) o.Copy(*this);
  return *this;
}

//_________________________________________________________
AliHFEtpcPIDqa::~AliHFEtpcPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}

//_________________________________________________________
void AliHFEtpcPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEtpcPIDqa &target = dynamic_cast<AliHFEtpcPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);
  target.fBrowseCentrality = fBrowseCentrality;
}

//_________________________________________________________
Long64_t AliHFEtpcPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEtpcPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEtpcPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEtpcPIDqa::Browse(TBrowser *b){
  //
  // Browse the PID QA
  //
  if(b){
    if(fHistos){
      b->Add(fHistos, fHistos->GetName());

      // Make Projections of the dE/dx Spectra and add them to a new Folder
      TString specnames[AliPID::kSPECIES+1] = {"All", "Electrons", "Muon", "Pions", "Kaon", "Protons"};
      Int_t specind[AliPID::kSPECIES+1] = {-1, AliPID::kElectron, AliPID::kMuon, AliPID::kPion, AliPID::kKaon, AliPID::kProton};
      TList *listdEdx = new TList;
      listdEdx->SetOwner();
      TList *listNsigma = new TList;
      listNsigma->SetOwner();

      TH2 *hptr = NULL; 
      for(Int_t ispec = 0; ispec < AliPID::kSPECIES+1; ispec++){
        for(Int_t istep = 0; istep < 2; istep++){
          hptr = MakeSpectrumdEdx(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec], fBrowseCentrality);
          hptr->SetName(Form("hTPCdEdx%s%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After" , fBrowseCentrality == -1 ? "MinBias" : Form("Cent%d", fBrowseCentrality)));
          listdEdx->Add(hptr);
          hptr = MakeSpectrumNSigma(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec], fBrowseCentrality);
          hptr->SetName(Form("hTPCnsigma%s%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After", fBrowseCentrality == -1 ? "MinBias" : Form("Cent%d", fBrowseCentrality)));
          listNsigma->Add(hptr);
        }
      }
      
      b->Add(listdEdx, "Projections dE/dx");
      b->Add(listNsigma, "Projections NSigma");
    }
  }
}

//_________________________________________________________
void AliHFEtpcPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("tpcqahistos", "Collection of TPC QA histograms");

  // Make common binning
  const Int_t kNdim = 5;
  const Int_t kPIDbins = AliPID::kSPECIES + 1;
  const Int_t kSteps = 2;
  const Int_t kCentralityBins = 11;
  const Double_t kMinPID = -1;
  const Double_t kMinP = 0.;
  const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
  const Double_t kMaxP = 20.;
  // Quantities where one can switch between low and high resolution
  Int_t kPbins = fQAmanager->HasHighResolutionHistos() ? 1000 : 100;
  Int_t kDedxbins = fQAmanager->HasHighResolutionHistos() ? 400 : 200;
  Int_t kSigmaBins = fQAmanager->HasHighResolutionHistos() ? 1400 : 240;
 
  // 1st histogram: TPC dEdx: (species, p, dEdx, step)
  Int_t nBinsdEdx[kNdim] = {kPIDbins, kPbins, kDedxbins, kSteps, kCentralityBins};
  Double_t mindEdx[kNdim] =  {kMinPID, kMinP, 0., 0., 0.};
  Double_t maxdEdx[kNdim] =  {kMaxPID, kMaxP, 200, 2., 11.}; 
  fHistos->CreateTHnSparse("tpcDedx", "TPC signal; species; p [GeV/c]; TPC signal [a.u.]; Selection Step; Centrality", kNdim, nBinsdEdx, mindEdx, maxdEdx);
  // 2nd histogram: TPC sigmas: (species, p nsigma, step)
  Int_t nBinsSigma[kNdim] = {kPIDbins, kPbins, kSigmaBins, kSteps, kCentralityBins};
  Double_t minSigma[kNdim] = {kMinPID, kMinP, -12., 0., 0.};
  Double_t maxSigma[kNdim] = {kMaxPID, kMaxP, 12., 2., 11.};
  fHistos->CreateTHnSparse("tpcnSigma", "TPC signal; species; p [GeV/c]; TPC signal [a.u.]; Selection Step; Centrality", kNdim, nBinsSigma, minSigma, maxSigma);

  // General TPC QA
}

//_________________________________________________________
void AliHFEtpcPIDqa::ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step){
  //
  // Fill TPC histograms
  //
  AliDebug(1, Form("QA started for TPC PID for step %d", (Int_t)step));
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t centrality = track->GetCentrality();
  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;

  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid));
  const AliPIDResponse *pidResponse = tpcpid ? tpcpid->GetPIDResponse() : NULL;
  Double_t contentSignal[5];
  contentSignal[0] = species;
  contentSignal[1] = tpcpid ? tpcpid->GetP(track->GetRecTrack(), anatype) : 0.;
  contentSignal[2] = GetTPCsignal(track->GetRecTrack(), anatype);
  contentSignal[3] = step;
  contentSignal[4] = centrality;
  fHistos->Fill("tpcDedx", contentSignal);

  contentSignal[2] = pidResponse ? pidResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron) : 0.; 
  fHistos->Fill("tpcnSigma", contentSignal);
}

//_________________________________________________________
Double_t AliHFEtpcPIDqa::GetTPCsignal(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype){
  //
  // Get TPC signal for ESD and AOD track
  //
  Double_t tpcSignal = 0.;
  if(anatype == AliHFEpidObject::kESDanalysis){
    const AliESDtrack *esdtrack = static_cast<const AliESDtrack *>(track);
    tpcSignal = esdtrack->GetTPCsignal();
  } else {
    const AliAODTrack *aodtrack = static_cast<const AliAODTrack *>(track);
    tpcSignal = aodtrack->GetDetPid() ? aodtrack->GetDetPid()->GetTPCsignal() : 0.;
  }
  return tpcSignal;
}

//_________________________________________________________
TH2 *AliHFEtpcPIDqa::MakeSpectrumdEdx(AliHFEdetPIDqa::EStep_t istep, Int_t species, Int_t centralityClass){
  //
  // Plot the Spectrum
  //
  THnSparseF *hSignal = dynamic_cast<THnSparseF *>(fHistos->Get("tpcDedx"));
  if(!hSignal) return NULL;
  hSignal->GetAxis(3)->SetRange(istep + 1, istep + 1);
  if(species >= 0 && species < AliPID::kSPECIES)
    hSignal->GetAxis(0)->SetRange(2 + species, 2 + species);
  TString hname = Form("hTPCsignal%s", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after"), 
          htitle = Form("TPC dE/dx Spectrum %s selection", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after");
  if(species > -1){
    hname += AliPID::ParticleName(species);
    htitle += Form(" for %ss", AliPID::ParticleName(species));
  }
  TString centname, centtitle;
  Bool_t hasCentralityInfo = kTRUE;
  if(centralityClass > -1){
    if(hSignal->GetNdimensions() < 5){
      AliError("Centrality Information not available");
      centname = centtitle = "MinBias";
      hasCentralityInfo= kFALSE;
    } else {
      // Project centrality classes
      // -1 is Min. Bias
      hSignal->GetAxis(4)->SetRange(centralityClass+1, centralityClass+1);
      centname = Form("Cent%d", centralityClass);
      centtitle = Form(" Centrality %d", centralityClass);
    }
  } else {
    centname = centtitle = "MinBias";
    hasCentralityInfo= kFALSE;
  }
  hname += centtitle;
  htitle += centtitle;

  TH2 *hTmp = hSignal->Projection(2,1);
  hTmp->SetName(hname.Data());
  hTmp->SetTitle(htitle.Data());
  hTmp->SetStats(kFALSE);
  hTmp->GetXaxis()->SetTitle("p [GeV/c]");
  hTmp->GetYaxis()->SetTitle("TPC signal [a.u.]");
  hSignal->GetAxis(3)->SetRange(0, hSignal->GetAxis(3)->GetNbins());
  hSignal->GetAxis(0)->SetRange(0, hSignal->GetAxis(0)->GetNbins());
  if(hasCentralityInfo) hSignal->GetAxis(4)->SetRange(0, hSignal->GetAxis(4)->GetNbins());
  return hTmp;
}

//_________________________________________________________
TH2 *AliHFEtpcPIDqa::MakeSpectrumNSigma(AliHFEdetPIDqa::EStep_t istep, Int_t species, Int_t centralityClass){
  //
  // Plot the Spectrum
  //
  THnSparseF *hSignal = dynamic_cast<THnSparseF *>(fHistos->Get("tpcnSigma"));
  if(!hSignal) return NULL;
  hSignal->GetAxis(3)->SetRange(istep + 1, istep + 1);
  if(species >= 0 && species < AliPID::kSPECIES)
    hSignal->GetAxis(0)->SetRange(2 + species, 2 + species);
  TString hname = Form("hTPCsigma%s", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after"), 
          htitle = Form("TPC dE/dx Spectrum[#sigma] %s selection", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after");
  if(species > -1){
    hname += AliPID::ParticleName(species);
    htitle += Form(" for %ss", AliPID::ParticleName(species));
  }
  TString centname, centtitle;
  Bool_t hasCentralityInfo = kTRUE;
  if(centralityClass > -1){
    if(hSignal->GetNdimensions() < 5){
      AliError("Centrality Information not available");
      centname = centtitle = "MinBias";
      hasCentralityInfo= kFALSE;
    } else {
      // Project centrality classes
      // -1 is Min. Bias
      hSignal->GetAxis(4)->SetRange(centralityClass+1, centralityClass+1);
      centname = Form("Cent%d", centralityClass);
      centtitle = Form(" Centrality %d", centralityClass);
    }
  } else {
    centname = centtitle = "MinBias";
    hasCentralityInfo= kFALSE;
  }
  hname += centtitle;
  htitle += centtitle;

  TH2 *hTmp = hSignal->Projection(2,1);
  hTmp->SetName(hname.Data());
  hTmp->SetTitle(htitle.Data());
  hTmp->SetStats(kFALSE);
  hTmp->GetXaxis()->SetTitle("p [GeV/c]");
  hTmp->GetYaxis()->SetTitle("TPC dE/dx - <dE/dx>|_{el} [#sigma]");
  hSignal->GetAxis(3)->SetRange(0, hSignal->GetAxis(3)->GetNbins());
  hSignal->GetAxis(0)->SetRange(0, hSignal->GetAxis(0)->GetNbins());
  if(hasCentralityInfo) hSignal->GetAxis(4)->SetRange(0, hSignal->GetAxis(4)->GetNbins());
  return hTmp;
}

