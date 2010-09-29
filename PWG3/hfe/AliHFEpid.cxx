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
// PID Steering Class 
// Interface to the user task
// Several strategies for Electron identification implemented.
// In addition users can combine different detectors or use
// single detector PID
// 
// Authors:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TAxis.h>
#include <TClass.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TIterator.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"

#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidITS.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTRD.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidMC.h"

ClassImp(AliHFEpid)

//____________________________________________________________
AliHFEpid::AliHFEpid():
  TNamed(),
  fEnabledDetectors(0),
  fPIDstrategy(kUndefined),
  fQAlist(0x0),
  fDebugLevel(0),
  fCommonObjects(NULL)
{
  //
  // Default constructor
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
}

//____________________________________________________________
AliHFEpid::AliHFEpid(const Char_t *name):
  TNamed(name, ""),
  fEnabledDetectors(0),
  fPIDstrategy(kUndefined),
  fQAlist(NULL),
  fDebugLevel(0),
  fCommonObjects(NULL)
{
  //
  // Default constructor
  // Create PID objects for all detectors
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
  fDetectorPID[kMCpid] = new AliHFEpidMC("MCPID");
  fDetectorPID[kTPCpid] = new AliHFEpidTPC("TRDPID");
  fDetectorPID[kTRDpid] = new AliHFEpidTRD("TRDPID");
  fDetectorPID[kTOFpid] = new AliHFEpidTOF("TOFPID");

}

//____________________________________________________________
AliHFEpid::AliHFEpid(const AliHFEpid &c):
  TNamed(c),
  fEnabledDetectors(c.fEnabledDetectors),
  fPIDstrategy(kUndefined),
  fQAlist(NULL),
  fDebugLevel(c.fDebugLevel),
  fCommonObjects(NULL)
{
  //
  // Copy Constructor
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
  c.Copy(*this);
}

//____________________________________________________________
AliHFEpid& AliHFEpid::operator=(const AliHFEpid &c){
  //
  // Assignment operator
  //
  if(&c != this) c.Copy(*this);
  return *this;
}

//____________________________________________________________
AliHFEpid::~AliHFEpid(){
  //
  // Destructor
  //
  for(Int_t idet = 0; idet < kNdetectorPID; idet++)
    if(fDetectorPID[idet]) delete fDetectorPID[idet];
  if(fQAlist) delete fQAlist; fQAlist = NULL;  // Each detector has to care about its Histograms
  ClearCommonObjects();
}

//____________________________________________________________
void AliHFEpid::Copy(TObject &o) const{
  //
  // Make copy
  //
  
  TNamed::Copy(o);
  AliHFEpid &target = dynamic_cast<AliHFEpid &>(o);
  target.ClearCommonObjects();

  target.fEnabledDetectors = fEnabledDetectors;
  target.fPIDstrategy = fPIDstrategy;
  if(target.fQAlist){
    delete target.fQAlist; target.fQAlist = NULL;
  }
  if(fQAlist) target.fQAlist = new TList;
  target.fQAlist = 0x0;
  target.fDebugLevel = fDebugLevel;
 
  // Copy detector PIDs
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    //Cleanup pointers in case of assignment
    if(target.fDetectorPID[idet])  
      delete target.fDetectorPID[idet];     
    if(fDetectorPID[idet]) 
      target.fDetectorPID[idet] = dynamic_cast<AliHFEpidBase *>(fDetectorPID[idet]->Clone());
  }
}

//____________________________________________________________
void AliHFEpid::AddCommonObject(TObject * const o){
  //
  // Add common object to the garbage collection
  //
  if(!fCommonObjects) fCommonObjects = new TObjArray;
  fCommonObjects->Add(o);
}

//____________________________________________________________
void AliHFEpid::ClearCommonObjects(){
  //
  // empty garbage collection
  //
  if(fCommonObjects){
    fCommonObjects->Delete();
    delete fCommonObjects;
    fCommonObjects = NULL;
  }
}

//____________________________________________________________
Bool_t AliHFEpid::InitializePID(TString arg){
  //
  // Initializes PID Object:
  // + Defines which detectors to use
  // + Initializes Detector PID objects
  // + Handles QA
  //
  
  Bool_t initFail = kFALSE;
  if(arg.BeginsWith("Strategy")){
    // Initialize detector PIDs according to PID Strategies
    arg.ReplaceAll("Strategy", "");
    fPIDstrategy = arg.Atoi();
    AliDebug(1, Form("%s - PID Strategy %d enabled", GetName(), fPIDstrategy));
    switch(fPIDstrategy){
      case 0: SwitchOnDetector(kMCpid); break;    // Pure MC PID - only valid in MC mode
      case 1: InitStrategy1(); break;
      case 2: InitStrategy2(); break;
      case 3: InitStrategy3(); break;
      case 4: InitStrategy4(); break;
      case 5: InitStrategy5(); break;
      case 6: InitStrategy6(); break;
      case 7: InitStrategy7(); break;
      case 8: InitStrategy8(); break;
      default: initFail = kFALSE;
    }
  } else {
    // No Strategy defined, Initialize according to detectors specified
    AliDebug(1, Form("%s - Doing InitializePID for Detectors %s end", GetName(), arg.Data()));
  
    TObjArray *detsEnabled = arg.Tokenize(":");
    TIterator *detIterator = detsEnabled->MakeIterator();
    TObjString *det = NULL;
    Int_t detector = -1;
    TString detectors[kNdetectorPID] = {"MC", "ESD", "ITS", "TPC", "TRD", "TOF"};
    Int_t nDetectors = 0;
    while((det = dynamic_cast<TObjString *>(detIterator->Next()))){
      TString &detstring = det->String();
      detector = -1;
      for(Int_t idet = 0; idet < kNdetectorPID; idet++){
        if(!detstring.CompareTo(detectors[idet])){
          detector = idet;
          break;
        }
      }
      if(detector > -1){
        SwitchOnDetector(detector);
        nDetectors++;
      } else AliError(Form("Detector %s not implemented (yet)", detstring.Data()));
    }
    if(!nDetectors) initFail = kTRUE;
  }
  if(initFail){
    AliError("Initializaion of the PID Failed");
    return kFALSE;
  }

  // Initialize PID Objects
  Bool_t status = kTRUE;
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(!IsDetectorOn(idet)) continue;
    if(fDetectorPID[idet]){ 
      status &= fDetectorPID[idet]->InitializePID();
      if(IsQAOn() && status) fDetectorPID[idet]->SetQAOn(fQAlist);
      if(HasMCData() && status) fDetectorPID[idet]->SetHasMCData();
    }
  }
  PrintStatus();
  return status;
  AliDebug(1, Form("%s - Done", GetName()));
}

//____________________________________________________________
Bool_t AliHFEpid::IsSelected(AliHFEpidObject *track){
  //
  // Steers PID decision for single detectors respectively combined
  // PID decision
  //
  if(!track->fRecTrack){
    // MC Event
    return (TMath::Abs(fDetectorPID[kMCpid]->IsSelected(track)) == 11);
  }
  if(fPIDstrategy < 9){
    AliDebug(1, Form("%s - PID Strategy %d", GetName(), fPIDstrategy));
    Int_t pid = 0;
    switch(fPIDstrategy){
      case 0: pid = IdentifyStrategy0(track); break;
      case 1: pid = IdentifyStrategy1(track); break;
      case 2: pid = IdentifyStrategy2(track); break;
      case 3: pid = IdentifyStrategy3(track); break;
      case 4: pid = IdentifyStrategy4(track); break;
      case 5: pid = IdentifyStrategy5(track); break;
      case 6: pid = IdentifyStrategy6(track); break;
      case 7: pid = IdentifyStrategy7(track); break;
      case 8: pid = IdentifyStrategy8(track); break;
      default: break;
    }
    return pid;
  }
  if(TESTBIT(fEnabledDetectors, kTPCpid)){
    if(IsQAOn() && fDebugLevel > 1){ 
      AliInfo("Filling QA plots");
      MakePlotsItsTpc(track);  // First fill the QA histograms
    }
    if(TESTBIT(fEnabledDetectors, kTOFpid)){
      // case TPC-TOF
      return MakePidTpcTof(track);
    } else if(TESTBIT(fEnabledDetectors, kTRDpid)){
      // case TPC-TRD with low level detector Signals
      return MakePidTpcTrd(track);
    } else
      return (TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) ==11);
  } else if(TESTBIT(fEnabledDetectors, kTRDpid)){
    return (TMath::Abs(fDetectorPID[kTRDpid]->IsSelected(track)) ==11);
  } else if(TESTBIT(fEnabledDetectors, kTOFpid)){
    return (TMath::Abs(fDetectorPID[kTOFpid]->IsSelected(track)) ==11);
  }
  
  return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEpid::MakePidTpcTof(AliHFEpidObject *track){
  //
  // Combines TPC and TOF PID decision
  //
  if(track->fAnalysisType != AliHFEpidObject::kESDanalysis) return kFALSE;

  AliHFEpidTOF *tofPID = dynamic_cast<AliHFEpidTOF*>(fDetectorPID[kTOFpid]);
  if(!tofPID){
    AliWarning("TOF pid object is NULL");
    return kFALSE;
  }

  AliHFEpidTPC *tpcPID = dynamic_cast<AliHFEpidTPC*>(fDetectorPID[kTPCpid]);
  if(!tpcPID){
    AliWarning("TPC pid object is NULL");
    return kFALSE;
  }
  
  // Use TOF PID to select particles with a sigma to the electron line > defined in initialize PID
  if(TMath::Abs(tofPID->IsSelected(track)) != 11) return kFALSE;
  // request TOF PID information, reject if no TOF information is available or if TOF identified as Proton or Kaon
  // apply cut only up to certain upper momentum cut
  /*Int_t pidTOF = tofPID->IsSelected(track);
  Bool_t isRejected = kFALSE;
  switch(TMath::Abs(pidTOF)){
    case 321:   if(track->fRecTrack->P() < 1.5) isRejected = kTRUE; break;
    case 2212:  if(track->fRecTrack->P() < 3) isRejected = kTRUE; break;
    case 0:     if(track->fRecTrack->P() < 3) isRejected = kTRUE; break;  // No TOF information available
    default: break;
  };
  if(isRejected) return kFALSE;*/

  // Particle passed TOF, let TPC decide, no line crossings defined anymore
 
  // request the tpc PID information
  return (TMath::Abs(tpcPID->IsSelected(track)) == 11);
}

//____________________________________________________________
Bool_t AliHFEpid::MakePidTpcTrd(AliHFEpidObject *track){
  //
  // Combination of TPC and TRD PID
  // Fills Histograms TPC Signal vs. TRD signal for different
  // momentum slices
  //
  if(track->fAnalysisType != AliHFEpidObject::kESDanalysis) return kFALSE; //AOD based detector PID combination not yet implemented
  AliDebug(1, "Analysis Type OK, do PID");
  AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(track->fRecTrack);
  AliHFEpidTRD *trdPid = dynamic_cast<AliHFEpidTRD *>(fDetectorPID[kTRDpid]);
  Int_t pdg = TMath::Abs(fDetectorPID[kMCpid]->IsSelected(track));
  Int_t pid = -1;
  switch(pdg){
    case 11:    pid = AliPID::kElectron; break;
    case 13:    pid = AliPID::kMuon; break;
    case 211:   pid = AliPID::kPion; break;
    case 321:   pid = AliPID::kKaon; break;
    case 2212:  pid = AliPID::kProton; break;
    default:    pid = -1;
  };
  Double_t content[10];
  content[0] = pid;
  content[1] = esdTrack->P();
  content[2] = esdTrack->GetTPCsignal();
  content[3] = trdPid->GetTRDSignalV1(esdTrack, pid);
  if(IsQAOn() && fDebugLevel > 1){
    content[4] = trdPid->GetTRDSignalV2(esdTrack, pid);
    AliDebug(1, Form("Momentum: %f, TRD Signal: Method 1[%f], Method 2[%f]", content[1], content[3], content[4]));
    (dynamic_cast<THnSparseF *>(fQAlist->At(kTRDSignal)))->Fill(content);
  }
  if(content[1] > 2){ // perform combined
    AliDebug(1, "Momentum bigger 2 GeV/c, doing combined PID"); 
    if(content[2] > 65 && content[3] > 500) return kTRUE;
    else return kFALSE;
  }
  else {
    AliDebug(1, "Momentum smaller 2GeV/c, doing TPC alone PID");
    return fDetectorPID[kTPCpid]->IsSelected(track) == 11;
  }
}

//____________________________________________________________
void AliHFEpid::MakePlotsItsTpc(AliHFEpidObject *track){
  //
  // Make a plot ITS signal - TPC signal for several momentum bins
  //
  if(track->fAnalysisType != AliHFEpidObject::kESDanalysis) return; //AOD based detector PID combination not yet implemented
  AliESDtrack * esdTrack = dynamic_cast<AliESDtrack *>(track->fRecTrack);
   // Fill My Histograms for MC PID
  Int_t pdg = TMath::Abs(fDetectorPID[kMCpid]->IsSelected(track));
  Int_t pid = -1;
  switch(pdg){
    case 11:    pid = AliPID::kElectron; break;
    case 13:    pid = AliPID::kMuon; break;
    case 211:   pid = AliPID::kPion; break;
    case 321:   pid = AliPID::kKaon; break;
    case 2212:  pid = AliPID::kProton; break;
    default:    pid = -1;
  };
  if(IsQAOn() && fDebugLevel > 0){
    Double_t content[10];
    content[0] = pid;
    content[1] = esdTrack->GetTPCInnerParam() ? esdTrack->GetTPCInnerParam()->P() : esdTrack->P();
    content[2] = (dynamic_cast<AliHFEpidITS *>(fDetectorPID[kITSpid]))->GetITSSignalV1(track->fRecTrack, pid);
    content[3] = esdTrack->GetTPCsignal();
    AliDebug(1, Form("Momentum %f, TPC Signal %f, ITS Signal %f", content[1], content[2], content[3]));
    (dynamic_cast<THnSparseF *>(fQAlist->At(kITSSignal)))->Fill(content);
  }
}

//____________________________________________________________
void AliHFEpid::SetESDpid(AliESDpid *pid){
  //
  // Set ESD PID to the Detector PID objects
  //
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]) fDetectorPID[idet]->SetESDpid(pid);
  }
}

//____________________________________________________________
void AliHFEpid::SetQAOn(){
  //
  // Switch on QA
  //
  SetBit(kIsQAOn, kTRUE);
  AliInfo("QA switched on");
  if(fQAlist) return;
  fQAlist = new TList;
  fQAlist->SetName("PIDqa");
  THnSparseF *histo = NULL;

  // Prepare axis for QA histograms
  const Int_t kMomentumBins = 41;
  const Double_t kPtMin = 0.1;
  const Double_t kPtMax = 10.;
  Double_t momentumBins[kMomentumBins];
  for(Int_t ibin = 0; ibin < kMomentumBins; ibin++)
    momentumBins[ibin] = static_cast<Double_t>(TMath::Power(10,TMath::Log10(kPtMin) + (TMath::Log10(kPtMax)-TMath::Log10(kPtMin))/(kMomentumBins-1)*static_cast<Double_t>(ibin)));

  // Add Histogram for combined TPC-TRD PID
  if(fDebugLevel > 1){
    AliDebug(1, "Adding histogram for ITS-TPC investigation");
    const Int_t kDimensionsTRDsig = 5;
    Int_t kNbinsTRDsig[kDimensionsTRDsig] = {AliPID::kSPECIES + 1, kMomentumBins - 1, 200, 3000, 3000};
    Double_t binMinTRDsig[kDimensionsTRDsig] = {-1., 0.1, 0, 0, 0};
    Double_t binMaxTRDsig[kDimensionsTRDsig] = {AliPID::kSPECIES, 10., 200., 3000., 3000.};
    fQAlist->AddAt((histo = new THnSparseF("fCombTPCTRDpid", "Combined TPC-TRD PID", kDimensionsTRDsig, kNbinsTRDsig, binMinTRDsig, binMaxTRDsig)), kTRDSignal);
    histo->GetAxis(1)->Set(kMomentumBins - 1, momentumBins);
    histo->GetAxis(0)->SetTitle("Particle Species");
    histo->GetAxis(1)->SetTitle("p / GeV/c");
    histo->GetAxis(2)->SetTitle("TPC Signal / a.u.");
    histo->GetAxis(3)->SetTitle("TRD Signal / a.u.");
    histo->GetAxis(4)->SetTitle("TRD Signal / a.u.");
  }

  // Add Histogram for combined TPC-ITS PID
  if(fDebugLevel > 0){
    AliDebug(1, "Adding histogram for TPC-TRD investigation");
    const Int_t kDimensionsITSsig = 4;
    Int_t kNbinsITSsig[kDimensionsITSsig] = {AliPID::kSPECIES + 1, kMomentumBins - 1, 300, 3000};
    Double_t binMinITSsig[kDimensionsITSsig] = {-1., 0.1, 0., 0.};
    Double_t binMaxITSsig[kDimensionsITSsig] = {AliPID::kSPECIES, 10., 300., 300.};
    fQAlist->AddAt((histo = new THnSparseF("fCombTPCITSpid", "Combined TPC-ITS PID", kDimensionsITSsig, kNbinsITSsig, binMinITSsig, binMaxITSsig)), kITSSignal);
    histo->GetAxis(1)->Set(kMomentumBins - 1, momentumBins);
    histo->GetAxis(0)->SetTitle("Particle Species");
    histo->GetAxis(1)->SetTitle("p / GeV/c");
    histo->GetAxis(2)->SetTitle("ITS Signal / a.u.");
    histo->GetAxis(3)->SetTitle("TPC Signal / a.u.");
  }
}

//____________________________________________________________
void AliHFEpid::InitStrategy1(){
  //
  // TPC alone, 3-sigma cut
  //
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  pid->SetTPCnSigma(1);
  SwitchOnDetector(kTPCpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy2(){
  //
  // TPC alone, symmetric 3 sigma cut and asymmetric sigma cut in the momentum region between 2GeV/c and 10 GeV/c and sigma between -1 and 100
  //
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  pid->SetTPCnSigma(3);
  pid->SetAsymmetricTPCsigmaCut(2., 10., 0., 4.);
  SwitchOnDetector(kTPCpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy3(){
  //
  // TPC alone, symmetric 3 sigma cut and 2 - -100 sigma pion rejection
  //   
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  pid->SetTPCnSigma(3);
  pid->SetRejectParticle(AliPID::kPion, 0., -100., 10., 1.);
  SwitchOnDetector(kTPCpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy4(){
  //
  // TPC and TRD combined, TPC 3 sigma cut and TRD NN 90% el efficiency level above 2 GeV/c
  //
  InitStrategy1();
  AliHFEpidTRD *trdpid = dynamic_cast<AliHFEpidTRD *>(fDetectorPID[kTRDpid]);
  trdpid->SetPIDMethod(AliHFEpidTRD::kLQ);
  trdpid->SetElectronEfficiency(0.71);
  SwitchOnDetector(kTRDpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy5(){
  //
  // TPC and TRD combined, TPC 3 sigma cut and TRD NN 90% el efficiency level above 2 GeV/c
  //
  InitStrategy1();
  SwitchOnDetector(kTRDpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy6(){
  //
  // Combined TPC-TOF PID, combination is discribed in the funtion MakePidTpcTof
  //
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fDetectorPID[kTOFpid]);
  tpcpid->SetTPCnSigma(2);
  tofpid->SetTOFnSigma(3);
  //TF1 *upperCut = new TF1("upperCut", "[0] * TMath::Exp([1]*x)", 0, 20);
  TF1 *upperCut = new TF1("upperCut", "[0]", 0, 20); // Use constant upper cut
  TF1 *lowerCut = new TF1("lowerCut", "[0] * TMath::Exp([1]*x) + [2]", 0, 20);
  upperCut->SetParameter(0, 5.);
  //upperCut->SetParameter(0, 2.7);
  //upperCut->SetParameter(1, -0.4357);
  lowerCut->SetParameter(0, -2.65);
  lowerCut->SetParameter(1, -0.8757);
//  lowerCut->SetParameter(2, -1);
  if(HasMCData()) lowerCut->SetParameter(2, -0.997);
  else lowerCut->SetParameter(2, -0.9);
  tpcpid->SetUpperSigmaCut(upperCut);
  tpcpid->SetLowerSigmaCut(lowerCut);
  AddCommonObject(upperCut);
  AddCommonObject(lowerCut);
  SwitchOnDetector(kTPCpid);
  SwitchOnDetector(kTOFpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy7(){
  //
  // TPC alone, symmetric 3 sigma cut and 2 - -100 sigma pion rejection
  //   
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  pid->SetTPCnSigma(2);
  pid->SetRejectParticle(AliPID::kProton, 0., -3., 10., 3.);
  pid->SetRejectParticle(AliPID::kKaon, 0., -3., 10., 3.);
  SwitchOnDetector(kTPCpid);
}

//____________________________________________________________
void AliHFEpid::InitStrategy8(){
  //
  // TOF, TRD and TPC together
  // 
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fDetectorPID[kTOFpid]);
  AliHFEpidTRD *trdpid = dynamic_cast<AliHFEpidTRD *>(fDetectorPID[kTRDpid]);

  tpcpid->SetTPCnSigma(3);
  tofpid->SetTOFnSigma(3);
  trdpid->SetPIDMethod(AliHFEpidTRD::kLQ);
  trdpid->SetElectronEfficiency(0.71);
  SwitchOnDetector(kTPCpid);
  SwitchOnDetector(kTOFpid);
  SwitchOnDetector(kTRDpid);
}


//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy0(AliHFEpidObject *track){
  return TMath::Abs(fDetectorPID[kMCpid]->IsSelected(track)) == 11;
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy1(AliHFEpidObject *track){
  return TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11;
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy2(AliHFEpidObject *track){
  return TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11;
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy3(AliHFEpidObject *track){
  return TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11;
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy4(AliHFEpidObject *track){
  Int_t trdpid = TMath::Abs(fDetectorPID[kTRDpid]->IsSelected(track));
  return (trdpid == 0 || trdpid == 11) && (TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11);
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy5(AliHFEpidObject *track){
  return MakePidTpcTrd(track);
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy6(AliHFEpidObject *track){
  return MakePidTpcTof(track);
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy7(AliHFEpidObject *track){
  //
  // Do PID in strategy 7: Proton and Kaon rejection and 
  // lower cut on TPC Signal
  //
  if(!(TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11))
    return kFALSE;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliHFEpid::IdentifyStrategy8(AliHFEpidObject *track){
  // 
  // Identify TPC, TRD, TOF
  //
  if(TMath::Abs(fDetectorPID[kTOFpid]->IsSelected(track)) != 11) return kFALSE;
  Int_t trdpid = TMath::Abs(fDetectorPID[kTRDpid]->IsSelected(track));
  return (trdpid == 0 || trdpid == 11) && (TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) == 11);
}

//____________________________________________________________
void AliHFEpid::PrintStatus() const {
  //
  // Print the PID configuration
  //
  printf("\n%s: Printing configuration\n", GetName());
  printf("===============================================\n");
  printf("PID Strategy: %d\n", fPIDstrategy);
  printf("PID Detectors: \n");
  Int_t npid = 0;
  TString detectors[kNdetectorPID] = {"MC", "ESD", "ITS", "TPC", "TRD", "TOF"};
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(IsDetectorOn(idet)){
      printf("\t%s\n", detectors[idet].Data());
      npid++;
    }
  }
  if(!npid) printf("\tNone\n");
  printf("\n");
}
