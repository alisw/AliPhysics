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
// QA class for Bayes PID
// Plot likelihoods vs various observables
// combined likelihood as well as indiv. detector response
// also plot mass spectrum using TOF for prior determination
//
// Author:
//   Yvonne Pachmayer <Y.Pachmayer@gsi.de>
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
#include "AliHFEpidTOF.h"
#include "AliHFEpidTRD.h"
#include "AliHFEpidBayes.h"
#include "AliHFEtools.h"
#include "AliHFEbayesPIDqa.h"
#include "AliHFEdetPIDqa.h"

ClassImp(AliHFEbayesPIDqa)

//_________________________________________________________
AliHFEbayesPIDqa::AliHFEbayesPIDqa():
    AliHFEdetPIDqa()
  , fHistos(NULL)
{
  //
  // Dummy constructor
  //
}

//_________________________________________________________
AliHFEbayesPIDqa::AliHFEbayesPIDqa(const char* name):
    AliHFEdetPIDqa(name, "QA for Bayes")
  , fHistos(NULL)
{
  //
  // Default constructor
  //
}

//____________________________________________________________
AliHFEbayesPIDqa::AliHFEbayesPIDqa(const AliHFEbayesPIDqa &o):
    AliHFEdetPIDqa(o),
    fHistos(NULL)
{
  //
  // Copy constructor
  //
}

//____________________________________________________________
AliHFEbayesPIDqa &AliHFEbayesPIDqa::operator=(const AliHFEbayesPIDqa &o){
  //
  // Make assignment
  //
  if(this != &o){
    AliHFEdetPIDqa::operator=(o);
    fHistos = o.fHistos;
  }
  return *this;
}


//_________________________________________________________
AliHFEbayesPIDqa::~AliHFEbayesPIDqa(){
  //
  // Destructor
  //
  if(fHistos) delete fHistos;
}
//_________________________________________________________
void AliHFEbayesPIDqa::Copy(TObject &o) const {
  //
  // Make copy
  //
  AliHFEbayesPIDqa &target = dynamic_cast<AliHFEbayesPIDqa &>(o);
  if(target.fHistos){
    delete target.fHistos;
    target.fHistos = NULL;
  }
  if(fHistos) target.fHistos = new AliHFEcollection(*fHistos);

}

//_________________________________________________________
Long64_t AliHFEbayesPIDqa::Merge(TCollection *coll){
  //
  // Merge with other objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;

  TIter it(coll);
  AliHFEbayesPIDqa *refQA = NULL;
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEbayesPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count + 1;
}
//_________________________________________________________
void AliHFEbayesPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("BAYESqahistos", "Collection of Bayes QA histograms");

  CreateDetectorSignalHistograms();

//  fBAYESpid = new AliHFEpidBayes("QAbayesPID");

}


void AliHFEbayesPIDqa::CreateDetectorSignalHistograms(){
  //
  // Create Histogram for Probability Studies
    //

    Int_t kPbins = 1000;  //fQAmanager->HasHighResolutionHistos() ? 1000 : 100;
    Int_t kSigmaBins = 300; //fQAmanager->HasHighResolutionHistos() ? 600 : 150;
    Int_t trdLikelihoodBins = 100; // fQAmanager->HasHighResolutionHistos() ? 200 : 100;
    const Int_t kPIDbins = AliPID::kSPECIES + 1;
    const Int_t kProbbins = 100;
    const Int_t kSteps = 2;
    const Int_t kCentralityBins = 11;
    const Double_t kMinPID = -1;
    const Double_t kMinP = 0.;
    const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
    const Double_t kMaxP = 10.;

    Int_t nBinsContr[8] =  {kPIDbins, kPbins, kProbbins, kSigmaBins, kSigmaBins, trdLikelihoodBins, kSteps, kCentralityBins};
    Double_t minContr[8] = {kMinPID, kMinP, 0,  -10., -7., 0., 0., 0.};
    Double_t maxContr[8] = {kMaxPID, kMaxP, 1,    5.,   8., 1., 2., 11.};
    fHistos->CreateTHnSparse("control", "Control; species; p [GeV/c]; Comb prob; TPC sigma; TOF sigma; TRD electron Likelihood; selection step; centrality", 8, nBinsContr, minContr,maxContr);



    Int_t nBinsTOFmass[5] = {kPIDbins, kPbins, 150, kSteps, kCentralityBins};
    Double_t minTOFmass[5] = {kMinPID, kMinP, 0., 0., 0.};
    Double_t maxTOFmass[5] = {kMaxPID, kMaxP, 1.5, 2., 11.};
    fHistos->CreateTHnSparse("tofmass", "TOF mass; species; p [GeV/c]; TOF mass; selection step", 5, nBinsTOFmass, minTOFmass, maxTOFmass);


}
//_________________________________________________________
void AliHFEbayesPIDqa::ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step){
  //
  // Fill TPC histograms
  //
  AliDebug(1, Form("QA started for BAYES PID for step %d", (Int_t)step));
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t centrality = track->GetCentrality();
  Int_t species = track->GetAbInitioPID();
  if(species >= AliPID::kSPECIES) species = -1;

  Double_t probComb[AliPID::kSPECIES]={0.};
  AliHFEpidBayes *bayespid = dynamic_cast<AliHFEpidBayes *>(fQAmanager->GetDetectorPID(AliHFEpid::kBAYESpid));
  const AliPIDResponse *pidResponseBayes = bayespid ? bayespid->GetPIDResponse() : NULL;
  if(!pidResponseBayes){
    AliError("No PID Response available");
    return;
  }
  bayespid->CalcCombProb(track,pidResponseBayes, probComb);

  Double_t contentSignal[8];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = probComb[AliPID::kElectron];
//  contentSignal[2] = contentSignal[2] = pidResponse ? pidResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron) : 0.;

  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fQAmanager->GetDetectorPID(AliHFEpid::kTOFpid));
  const AliPIDResponse *pidResponseTOF = tofpid ? tofpid->GetPIDResponse() : NULL;
  if(!pidResponseTOF){
    AliError("No PID response available");
    return;
  }
  contentSignal[4] = pidResponseTOF->NumberOfSigmasTOF(track->GetRecTrack(), AliPID::kElectron);

  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid));
  const AliPIDResponse *pidResponse = tpcpid ? tpcpid->GetPIDResponse() : NULL;
  if(!pidResponse){
    AliError("No PID response available");
    return;
  }
  contentSignal[3] = pidResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron);

  AliHFEpidTRD *trdpid = dynamic_cast<AliHFEpidTRD *>(fQAmanager->GetDetectorPID(AliHFEpid::kTRDpid));
  contentSignal[5] = trdpid ? trdpid->GetElectronLikelihood(static_cast<const AliVTrack*>(track->GetRecTrack()), anatype) : -10;

  contentSignal[6] = step;
  contentSignal[7] = centrality;
  fHistos->Fill("control", contentSignal);

  Double_t contentFill[5];
  contentFill[0]=contentSignal[0];
  contentFill[1]=contentSignal[1];
  contentFill[3]=contentSignal[6];
  contentFill[4]=contentSignal[7];

  Double_t masscalcfromtof=CalcTOFMass(track);
  contentFill[2]=masscalcfromtof;
  fHistos->Fill("tofmass", contentFill);
}

Double_t AliHFEbayesPIDqa::CalcTOFMass(const AliHFEpidObject *track){
    //
    // Calc TOF mass
    //

   Double_t mass=-99; //GeV
   Double_t length=((AliESDtrack*)track->GetRecTrack())->GetIntegratedLength();
   if (length<=0) return mass;
   Double_t tofTime=((AliESDtrack*)track->GetRecTrack())->GetTOFsignal();//in ps
   Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector
   if (tof<=0) return mass;
   Double_t c=TMath::C()*1.E-9;// m/ns
   length =length*0.01; // in meters
   tof=tof*c;
   Double_t fact= (tof/length)*(tof/length) -1.;
   Double_t mom=track->GetRecTrack()->P();
   if(mom==0) return mass;

   if(fact<=0) {
       mass = -mom*TMath::Sqrt(-fact);
   }else{
       mass = mom*TMath::Sqrt(fact);
   }

   return mass;
}

