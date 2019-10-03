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
// Class AliHFEitsPIDqa
// Monitoring ITS PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   ITS dE/dx (Number of sigmas)
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//    Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
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
#include "AliHFEpidITS.h"
#include "AliHFEtools.h"
#include "AliHFEitsPIDqa.h"

ClassImp(AliHFEitsPIDqa)

//_________________________________________________________
AliHFEitsPIDqa::AliHFEitsPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
  , fBrowseCentrality(-1)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEitsPIDqa::AliHFEitsPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for ITS")
  , fHistos(NULL)
  , fBrowseCentrality(-1)
{
  //
  // Default constructor
  //
}

//_________________________________________________________
AliHFEitsPIDqa::AliHFEitsPIDqa(const AliHFEitsPIDqa &o):
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
AliHFEitsPIDqa &AliHFEitsPIDqa::operator=(const AliHFEitsPIDqa &o){
  //
  // Do assignment
  //
  AliHFEdetPIDqa::operator=(o);
  if(&o != this) o.Copy(*this);
  return *this;
}

//_________________________________________________________
AliHFEitsPIDqa::~AliHFEitsPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}

//_________________________________________________________
void AliHFEitsPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEitsPIDqa &target = dynamic_cast<AliHFEitsPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);
  target.fBrowseCentrality = fBrowseCentrality;
}

//_________________________________________________________
Long64_t AliHFEitsPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEitsPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEitsPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEitsPIDqa::Browse(TBrowser *b){
  //
  // Browse the PID QA
  //
  if(b){
    if(fHistos){
      b->Add(fHistos, fHistos->GetName());

      // Make Projections of the nsigma Spectra and add them to a new Folder
      TString specnames[AliPID::kSPECIES+1] = {"All", "Electrons", "Muon", "Pions", "Kaon", "Protons"};
      Int_t specind[AliPID::kSPECIES+1] = {-1, AliPID::kElectron, AliPID::kMuon, AliPID::kPion, AliPID::kKaon, AliPID::kProton};
      TList *listNsigma = new TList;
      listNsigma->SetOwner();

      TH2 *hptr = NULL; 
      for(Int_t ispec = 0; ispec < AliPID::kSPECIES+1; ispec++){
        for(Int_t istep = 0; istep < 2; istep++){
          hptr = MakeSpectrumNSigma(static_cast<AliHFEdetPIDqa::EStep_t>(istep), specind[ispec], fBrowseCentrality);
          hptr->SetName(Form("hITSnsigma%s%s%s", specnames[ispec].Data(), istep == 0 ? "Before" : "After", fBrowseCentrality == -1 ? "MinBias" : Form("Cent%d", fBrowseCentrality)));
          listNsigma->Add(hptr);
        }
      }
      
      b->Add(listNsigma, "Projections NSigma");
    }
  }
}

//_________________________________________________________
void AliHFEitsPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("itsqahistos", "Collection of ITS QA histograms");

  // Make common binning
  const Int_t kNdim = 5;
  const Int_t kPIDbins = AliPID::kSPECIES + 1;
  const Int_t kSteps = 2;
  const Int_t kCentralityBins = 11;
  const Double_t kMinPID = -1;
  const Double_t kMinP = 0.;
  const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
  const Double_t kMaxP = 2.;

  // Quantities where one can switch between low and high resolution
  Int_t kPbins = fQAmanager->HasHighResolutionHistos() ? 100 : 10;
  Int_t kSigmaBins = fQAmanager->HasHighResolutionHistos() ? 1400 : 240;
 
  // 1st histogram: ITS sigmas: (species, p, nsigma, step, centrality)
  Int_t nBinsSigma[kNdim] = {kPIDbins, kPbins, kSigmaBins, kSteps, kCentralityBins};
  Double_t minSigma[kNdim] = {kMinPID, kMinP, -12., 0., 0.};
  Double_t maxSigma[kNdim] = {kMaxPID, kMaxP, 12., 2., 11.};
  fHistos->CreateTHnSparse("itsnSigma", "ITS signal; species; p [GeV/c]; ITS signal [a.u.]; Selection Step; Centrality; Eta", kNdim, nBinsSigma, minSigma, maxSigma);

  // General ITS QA
}

//_________________________________________________________
void AliHFEitsPIDqa::ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step){
  //
  // Fill ITS histograms
  //
  AliDebug(1, Form("QA started for ITS PID for step %d", (Int_t)step));
  //AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t centrality = track->GetCentrality();
  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;

  AliHFEpidITS *itspid = dynamic_cast<AliHFEpidITS *>(fQAmanager->GetDetectorPID(AliHFEpid::kITSpid));
  Double_t contentSignal[5];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = itspid ? itspid->GetITSNsigmaCorrected(track->GetRecTrack()) : -100.; 
  contentSignal[3] = step;
  contentSignal[4] = centrality;
 
  fHistos->Fill("itsnSigma", contentSignal);
}

//_________________________________________________________
TH2 *AliHFEitsPIDqa::MakeSpectrumNSigma(AliHFEdetPIDqa::EStep_t istep, Int_t species, Int_t centralityClass){
  //
  // Plot the Spectrum
  //
  THnSparseF *hSignal = dynamic_cast<THnSparseF *>(fHistos->Get("itsnSigma"));
  if(!hSignal) return NULL;
  hSignal->GetAxis(3)->SetRange(istep + 1, istep + 1);
  if(species >= 0 && species < AliPID::kSPECIES)
    hSignal->GetAxis(0)->SetRange(2 + species, 2 + species);
  TString hname = Form("hITSsigma%s", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after"), 
          htitle = Form("ITS dE/dx Spectrum[#sigma] %s selection", istep == AliHFEdetPIDqa::kBeforePID ? "before" : "after");
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
  hTmp->GetYaxis()->SetTitle("ITS dE/dx - <dE/dx>|_{el} [#sigma]");
  hSignal->GetAxis(3)->SetRange(0, hSignal->GetAxis(3)->GetNbins());
  hSignal->GetAxis(0)->SetRange(0, hSignal->GetAxis(0)->GetNbins());
  if(hasCentralityInfo) hSignal->GetAxis(4)->SetRange(0, hSignal->GetAxis(4)->GetNbins());
  return hTmp;
}

