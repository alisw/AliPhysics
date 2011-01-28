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
// Class AliHFEtofPIDqa
//
// Monitoring TOF PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   TOF time distance to electron hypothesis (Number of sigmas)
//   TOF time distance to the pion hypothesis (Absolute value)
//   TPC dE/dx (Number of sigmas, as control histogram)
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//
//   Markus Fasel <M.Fasel@gsi.de>
#include <TClass.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"

#include "AliHFEcollection.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTOF.h"
#include "AliHFEtools.h"
#include "AliHFEtofPIDqa.h"

ClassImp(AliHFEpidTOF)

//_________________________________________________________
AliHFEtofPIDqa::AliHFEtofPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEtofPIDqa::AliHFEtofPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for TOF")
  , fHistos(NULL)
{
  //
  // Default constructor
  //
}

//_________________________________________________________
AliHFEtofPIDqa::AliHFEtofPIDqa(const AliHFEtofPIDqa &o):
    AliHFEdetPIDqa(o)
  , fHistos()
{
  //
  // Copy constructor
  //
  o.Copy(*this);
}

//_________________________________________________________
AliHFEtofPIDqa &AliHFEtofPIDqa::operator=(const AliHFEtofPIDqa &o){
  //
  // Do assignment
  //
  AliHFEdetPIDqa::operator=(o);
  if(&o != this) o.Copy(*this);
  return *this;
}

//_________________________________________________________
AliHFEtofPIDqa::~AliHFEtofPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}

//_________________________________________________________
void AliHFEtofPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEtofPIDqa &target = dynamic_cast<AliHFEtofPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);
}

//_________________________________________________________
Long64_t AliHFEtofPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEtofPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEtofPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEtofPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("tofqahistos", "Collection of TOF QA histograms");

  // Make common binning
  const Int_t kPIDbins = AliPID::kSPECIES + 1;
  const Int_t kSteps = 2;
  const Double_t kMinPID = -1;
  const Double_t kMinP = 0.;
  const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
  const Double_t kMaxP = 20.;
  // Quantities where one can switch between low and high resolution
  Int_t kPbins = fQAmanager->HasHighResolutionHistos() ?  1000 : 100;
  Int_t kSigmaBins = fQAmanager->HasHighResolutionHistos() ? 1400 : 240;

  // 1st histogram: TOF sigmas: (species, p nsigma, step)
  Int_t nBinsSigma[4] = {kPIDbins, kPbins, kSigmaBins, kSteps};
  Double_t minSigma[4] = {kMinPID, kMinP, -12., 0};
  Double_t maxSigma[4] = {kMaxPID, kMaxP, 12., 2.};
  fHistos->CreateTHnSparse("tofnSigma", "TOF signal; species; p [GeV/c]; TOF signal [a.u.]; Selection Step", 4, nBinsSigma, minSigma, maxSigma);

  // 2nd histogram: TOF time - pion hypothesis (TOF Time Resolution Monitor)
  fHistos->CreateTH2F("tofTimeRes", "Difference between measured and expected time for Pions; p [GeV/c]; #Deltat [ps]", 100, 0.1, 10, 100, -200, 200); 
  
  // 3rd histogram: TPC sigmas to the electron line: (species, p nsigma, step - only filled if apriori PID information available)
  fHistos->CreateTHnSparse("tofMonitorTPC", "TPC signal; species; p [GeV/c]; TPC signal [a.u.]; Selection Step", 4, nBinsSigma, minSigma, maxSigma);
}

//_________________________________________________________
void AliHFEtofPIDqa::ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step){
  //
  // Fill TPC histograms
  //
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;

  AliDebug(1, Form("Monitoring particle of type %d for step %d", species, step));
  AliHFEpidTOF *tofpid= dynamic_cast<AliHFEpidTOF *>(fQAmanager->GetDetectorPID(AliHFEpid::kTOFpid));
  AliHFEpidTPC *tpcpid= dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid));
  
  Double_t contentSignal[4];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = tofpid ? tofpid->NumberOfSigmas(track->GetRecTrack(), AliPID::kElectron, anatype): 0.;
  contentSignal[3] = step;
  fHistos->Fill("tofnSigma", contentSignal);
  if(tofpid){
    Double_t timeTof = tofpid->GetTOFsignal(track->GetRecTrack(), anatype);
    Double_t time0 = tofpid->GetTime0(anatype);
    Double_t tof = timeTof - time0;
    Double_t times[AliPID::kSPECIES]; tofpid->GetIntegratedTimes(track->GetRecTrack(), times, anatype);
    fHistos->Fill("tofTimeRes",contentSignal[1], tof - times[AliPID::kPion]);
  }
  if(species > -1){
    contentSignal[2] = tpcpid ? tpcpid->NumberOfSigmas(track->GetRecTrack(), AliPID::kElectron, anatype) : 0.;
    fHistos->Fill("tofMonitorTPC", contentSignal);
  }
}

//_________________________________________________________
TH2 *AliHFEtofPIDqa::MakeSpectrumNSigma(AliHFEdetPIDqa::EStep_t istep, Int_t species){
  //
  // Plot the Spectrum
  //
  THnSparseF *hSignal = dynamic_cast<THnSparseF *>(fHistos->Get("tofnSigma"));
  if(!hSignal) return NULL;
  hSignal->GetAxis(3)->SetRange(istep + 1, istep + 1);
  if(species > 0 && species < AliPID::kSPECIES)
    hSignal->GetAxis(0)->SetRange(2 + species, 2 + species);
  TH2 *hTmp = hSignal->Projection(2,1);
  TString hname = Form("hTPCsigma%s", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after"), 
          htitle = Form("TPC dE/dx Spectrum[#sigma] %s selection", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after");
  if(species > -1){
    hname += AliPID::ParticleName(species);
    htitle += Form(" for %ss", AliPID::ParticleName(species));
  }
  hTmp->SetName(hname.Data());
  hTmp->SetTitle(htitle.Data());
  hTmp->SetStats(kFALSE);
  hTmp->GetXaxis()->SetTitle("p [GeV/c]");
  hTmp->GetYaxis()->SetTitle("TOF time|_{el} - expected time|_{el} [#sigma]");
  hSignal->GetAxis(3)->SetRange(0, hSignal->GetAxis(3)->GetNbins());
  hSignal->GetAxis(0)->SetRange(0, hSignal->GetAxis(0)->GetNbins());
  return hTmp;
}

//_________________________________________________________
TH1 *AliHFEtofPIDqa::GetHistogram(const char *name){
  // 
  // Get the histogram
  //
  if(!fHistos) return NULL;
  TObject *histo = fHistos->Get(name);
  if(!histo->InheritsFrom("TH1")) return NULL;
  return dynamic_cast<TH1 *>(histo);
}

