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

//----------------------------------------------------------------
// Definition of common binning
const Int_t kNdim = 5;
const Int_t kPIDbins = AliPID::kSPECIES + 1;
const Int_t kProbbins = 100;
const Int_t kSteps = 2;
const Int_t kCentralityBins = 11;
const Double_t kMinPID = -1;
const Double_t kMinP = 0.;
const Double_t kMaxPID = (Double_t)AliPID::kSPECIES;
const Double_t kMaxP = 20.;

const Int_t AliHFEbayesPIDqa::fgkNBinsProb[kNdim] = {
    kPIDbins,                     // species
    1000,                       // p-bins
    kProbbins,                    // probability
    kSteps,                       // step
    kCentralityBins               // centrality
};
const Double_t AliHFEbayesPIDqa::fgkMinBinsProb[kNdim] = {
    -1.,                          // species
    0.,                           // p-bins
    0.,                           // probability
    0.,                           // step
    0.                            // centrality
};
const Double_t AliHFEbayesPIDqa::fgkMaxBinsProb[kNdim] = {
    (Double_t)AliPID::kSPECIES,   // species
    20.,                          // p-bins
    1.,                           // probability
    2.,                           // step
    11.                           // centrality
};
//----------------------------------------------------------------

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
  AliHFEdetPIDqa::operator=(o);
  fHistos = o.fHistos;
  
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
void AliHFEbayesPIDqa::Initialize(){
  //
  // Define Histograms
  //

  fHistos = new AliHFEcollection("BAYESqahistos", "Collection of Bayes QA histograms");

  CreateProbabilityHistograms();
  CreateDetectorSignalHistograms();

//  fBAYESpid = new AliHFEpidBayes("QAbayesPID");

}

//__________________________________________________________________
void AliHFEbayesPIDqa::CreateProbabilityHistograms(){
  //
  // Create Histogram for Probability Studies
  //

  fHistos->CreateTHnSparse("combprob", "Comb prob; species; p [GeV/c]; Comb prob; Selection Step; Centrality", kNdim, fgkNBinsProb, fgkMinBinsProb, fgkMaxBinsProb);


}

void AliHFEbayesPIDqa::CreateDetectorSignalHistograms(){
  //
  // Create Histogram for Probability Studies
  //
    Int_t kPbins = fQAmanager->HasHighResolutionHistos() ? 1000 : 100;
    Int_t kSigmaBins = fQAmanager->HasHighResolutionHistos() ? 1400 : 240;
    Int_t nBinsSigma[kNdim] = {kPIDbins, kPbins, kSigmaBins, kSteps, kCentralityBins};
    Double_t minSigma[kNdim] = {kMinPID, kMinP, -12., 0., 0.};
    Double_t maxSigma[kNdim] = {kMaxPID, kMaxP, 12., 2., 11.};
    fHistos->CreateTHnSparse("tpcsigma", "TPC sigma; species; p [GeV/c]; TPC sigma; Selection Step; Centrality", kNdim, nBinsSigma, minSigma,maxSigma);
    fHistos->CreateTHnSparse("tofsigma", "TOF sigma; species; p [GeV/c]; TOF sigma; Selection Step; Centrality", kNdim, nBinsSigma, minSigma,maxSigma);

    Int_t trdLikelihoodBins = fQAmanager->HasHighResolutionHistos() ? 200 : 100;
    Int_t nBinsTRDlike[5] = {kPIDbins, kPbins, trdLikelihoodBins, kSteps, kCentralityBins};
    Double_t minTRDlike[5] = {kMinPID, kMinP, 0., 0., 0.};
    Double_t maxTRDlike[5] = {kMaxPID, kMaxP, 1., 2., 11.};
    fHistos->CreateTHnSparse("trdlikeli", "TRD Likelihood Distribution; species; p [GeV/c]; TRD electron Likelihood; selection step", 5, nBinsTRDlike, minTRDlike, maxTRDlike);

    Int_t nBinsContr[6] = {kPIDbins, kPbins, kProbbins, kSigmaBins, kSigmaBins, trdLikelihoodBins};
    Double_t minContr[6] = {kMinPID, kMinP, 0,  -12., -12., 0.};
    Double_t maxContr[6] = {kMaxPID, kMaxP, 1,   12.,  12., 1.};
    fHistos->CreateTHnSparse("control", "Control; species; p [GeV/c]; Comb prob; TPC sigma; TOF sigma; TRD electron Likelihood", 6, nBinsContr, minContr,maxContr);

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
  bayespid->CalcCombProb(track,pidResponseBayes, probComb);

  Double_t contentSignal[5];
  contentSignal[0] = species;
  contentSignal[1] = track->GetRecTrack()->P();
  contentSignal[2] = probComb[AliPID::kElectron];
//  contentSignal[2] = contentSignal[2] = pidResponse ? pidResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron) : 0.;
  contentSignal[3] = step;
  contentSignal[4] = centrality;
  fHistos->Fill("combprob", contentSignal);

  Double_t contentContr[6];
  contentContr[0]=contentSignal[0];
  contentContr[1]=contentSignal[1];
  contentContr[2]=contentSignal[2];

  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fQAmanager->GetDetectorPID(AliHFEpid::kTOFpid));
  const AliPIDResponse *pidResponseTOF = tofpid ? tofpid->GetPIDResponse() : NULL;
  contentSignal[2] = pidResponseTOF ? pidResponseTOF->NumberOfSigmasTOF(track->GetRecTrack(), AliPID::kElectron): -10.;
  fHistos->Fill("tofsigma", contentSignal);
  contentContr[3]=contentSignal[2];

  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fQAmanager->GetDetectorPID(AliHFEpid::kTPCpid));
  const AliPIDResponse *pidResponse = tpcpid ? tpcpid->GetPIDResponse() : NULL;
  contentSignal[2] = pidResponse ? pidResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron) : -10.;
  fHistos->Fill("tpcsigma", contentSignal);
  contentContr[4]=contentSignal[2];


  AliHFEpidTRD *trdpid = dynamic_cast<AliHFEpidTRD *>(fQAmanager->GetDetectorPID(AliHFEpid::kTRDpid));
  contentSignal[2] = trdpid ? trdpid->GetElectronLikelihood(static_cast<const AliVTrack*>(track->GetRecTrack()), anatype) : -10;
  fHistos->Fill("trdlikeli", contentSignal);

  contentContr[5]=contentSignal[2];
  fHistos->Fill("control", contentContr);

  Double_t masscalcfromtof=CalcTOFMass(track);
  contentSignal[2]=masscalcfromtof;
  fHistos->Fill("tofmass", contentSignal);
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

