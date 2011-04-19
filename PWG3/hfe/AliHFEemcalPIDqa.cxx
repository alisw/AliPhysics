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
// Class AliHFEemcalPIDqa
//
// Monitoring EMCAL PID in the HFE PID montioring framework. The following
// quantities are monitored:
//   TPC Signal distribution 
// (Always as function of momentum, particle species and centrality 
// before and after cut)
// More information about the PID monitoring framework can be found in
// AliHFEpidQAmanager.cxx and AliHFEdetPIDqa.cxx
//
// Author:
//
//   Shingo Sakai <ssakai@lbl.gov>
#include <TClass.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliHFEcollection.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidEMCAL.h"
#include "AliHFEtools.h"
#include "AliHFEemcalPIDqa.h"

ClassImp(AliHFEpidEMCAL)

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for EMCAL")
  , fHistos(NULL)
{
  //
  // Default constructor
  //
}

//_________________________________________________________
AliHFEemcalPIDqa::AliHFEemcalPIDqa(const AliHFEemcalPIDqa &o):
    AliHFEdetPIDqa(o)
  , fHistos()
{
  //
  // Copy constructor
  //
  o.Copy(*this);
}

//_________________________________________________________
AliHFEemcalPIDqa &AliHFEemcalPIDqa::operator=(const AliHFEemcalPIDqa &o){
  //
  // Do assignment
  //
  AliHFEdetPIDqa::operator=(o);
  if(&o != this) o.Copy(*this);
  return *this;
}

//_________________________________________________________
AliHFEemcalPIDqa::~AliHFEemcalPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}

//_________________________________________________________
void AliHFEemcalPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEemcalPIDqa &target = dynamic_cast<AliHFEemcalPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);
}

//_________________________________________________________
Long64_t AliHFEemcalPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEemcalPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEemcalPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}

//_________________________________________________________
void AliHFEemcalPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("emcalqahistos", "Collection of EMCAL QA histograms");

  // Make common binning
  const Int_t kCentralityBins = 11;
  const Double_t kMinP = 0.;
  const Double_t kMaxP = 20.;
  const Double_t kTPCSigMim = 40.;
  const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, Centrality)
  Int_t nBins[4] = {20, 20, 400, kCentralityBins};
  Double_t min[4] = {kMinP, kMinP, kTPCSigMim,  0};
  Double_t max[4] = {kMaxP, kMaxP, kTPCSigMax,  11.};
  fHistos->CreateTHnSparse("EMCAL_TPCdedx", "EMCAL signal; species; p [GeV/c]; pT [GeV/c] ; TPC signal [a.u.]; Centrality", 4, nBins, min, max);
  
}

//_________________________________________________________
void AliHFEemcalPIDqa::ProcessTrack(const AliHFEpidObject *track,AliHFEdetPIDqa::EStep_t /*step*/){
  //
  // Fill TPC histograms
  //
  //AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t centrality = track->GetCentrality();

  const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(vtrack);

  Double_t contentSignal[4];
  contentSignal[0] = track->GetRecTrack()->P();
  contentSignal[1] = track->GetRecTrack()->Pt();
  contentSignal[2] = esdtrack->GetTPCsignal(); //?
  contentSignal[3] = centrality;
  //printf("ProcessTrack ; Print Content %g; %g; %g; %g \n",contentSignal[0],contentSignal[1],contentSignal[2],contentSignal[3]); 
  fHistos->Fill("EMCAL_TPCdedx", contentSignal);
}

//_________________________________________________________
TH1 *AliHFEemcalPIDqa::GetHistogram(const char *name){
  // 
  // Get the histogram
  //
  if(!fHistos) return NULL;
  TObject *histo = fHistos->Get(name);
  if(!histo->InheritsFrom("TH1")) return NULL;
  return dynamic_cast<TH1 *>(histo);
}

