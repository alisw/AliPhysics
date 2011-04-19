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
#include <TIterator.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliAODpidUtil.h"
#include "AliESDpid.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"

#include "AliHFEcontainer.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidITS.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTRD.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidEMCAL.h"
#include "AliHFEpidMC.h"
#include "AliHFEvarManager.h"

ClassImp(AliHFEpid)

const Char_t* AliHFEpid::fgkDetectorName[AliHFEpid::kNdetectorPID + 1] = {
  "MCPID",
  "ESDPID",
  "ITSPID",
  "TPCPID",
  "TRDPID",
  "TOFPID",
  "EMCALPID",
  "UndefinedPID"
};

//____________________________________________________________
AliHFEpid::AliHFEpid():
  TNamed(),
  fEnabledDetectors(0),
  fNPIDdetectors(0),
  fVarManager(NULL),
  fCommonObjects(NULL)
{
  //
  // Default constructor
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
  memset(fDetectorOrder, kUndefined, sizeof(UInt_t) * kNdetectorPID);
  memset(fDetectorOrder, 0, sizeof(UInt_t) * kNdetectorPID);
}

//____________________________________________________________
AliHFEpid::AliHFEpid(const Char_t *name):
  TNamed(name, ""),
  fEnabledDetectors(0),
  fNPIDdetectors(0),
  fVarManager(NULL),
  fCommonObjects(NULL)
{
  //
  // Default constructor
  // Create PID objects for all detectors
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
  memset(fDetectorOrder, kUndefined, sizeof(UInt_t) * kNdetectorPID);
  memset(fSortedOrder, 0, sizeof(UInt_t) * kNdetectorPID);

  fDetectorPID[kMCpid] = new AliHFEpidMC("MCPID");
  fDetectorPID[kTPCpid] = new AliHFEpidTPC("TPCPID");
  fDetectorPID[kTRDpid] = new AliHFEpidTRD("TRDPID");
  fDetectorPID[kTOFpid] = new AliHFEpidTOF("TOFPID");
  fDetectorPID[kEMCALpid] = new AliHFEpidEMCAL("EMCALPID");

}

//____________________________________________________________
AliHFEpid::AliHFEpid(const AliHFEpid &c):
  TNamed(c),
  fEnabledDetectors(c.fEnabledDetectors),
  fNPIDdetectors(c.fNPIDdetectors),
  fVarManager(c.fVarManager),
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
  target.fNPIDdetectors = fNPIDdetectors;
  target.fVarManager = fVarManager;
 
  // Copy detector PIDs
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    //Cleanup pointers in case of assignment
    if(target.fDetectorPID[idet])  
      delete target.fDetectorPID[idet];     
    if(fDetectorPID[idet]) 
      target.fDetectorPID[idet] = dynamic_cast<AliHFEpidBase *>(fDetectorPID[idet]->Clone());
  }
  memcpy(target.fDetectorOrder, fDetectorOrder, sizeof(UInt_t) * kNdetectorPID);
  memcpy(target.fSortedOrder, fSortedOrder, sizeof(UInt_t) * kNdetectorPID);
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
void AliHFEpid::AddDetector(TString detector, UInt_t position){
  //
  // Add Detector in position 
  //
  UInt_t detectorID = kUndefined;
  detector.ToUpper();
  if(!detector.CompareTo("MC")) detectorID = kMCpid;
  else if(!detector.CompareTo("TPC")) detectorID = kTPCpid;
  else if(!detector.CompareTo("TRD")) detectorID = kTRDpid;
  else if(!detector.CompareTo("TOF")) detectorID = kTOFpid;
  else if(!detector.CompareTo("EMCAL")) detectorID = kEMCALpid;
  else AliError("Detector not available");

  if(detectorID == kUndefined) return;
  if(IsDetectorOn(detectorID)) return;
  SwitchOnDetector(detectorID);
  fDetectorOrder[detectorID] = position;
  fNPIDdetectors++;
}

//____________________________________________________________
Bool_t AliHFEpid::InitializePID(){
  //
  // Initializes the PID object
  //

  TMath::Sort(static_cast<UInt_t>(kNdetectorPID), fDetectorOrder, fSortedOrder, kFALSE);
  // Initialize PID Objects
  Bool_t status = kTRUE;
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(!IsDetectorOn(idet)) continue;
    if(fDetectorPID[idet]){ 
      status &= fDetectorPID[idet]->InitializePID();
      if(HasMCData() && status) fDetectorPID[idet]->SetHasMCData();
    }
  }
  PrintStatus();
  return status;
}

//____________________________________________________________
Bool_t AliHFEpid::IsSelected(AliHFEpidObject *track, AliHFEcontainer *cont, const Char_t *contname, AliHFEpidQAmanager *pidqa){
  //
  // Select Tracks
  //
  Bool_t isSelected = kTRUE;
  AliDebug(1, Form("Particle used for PID, QA available: %s", pidqa ? "Yes" : "No"));
  for(UInt_t idet = 0; idet < fNPIDdetectors; idet++){
    AliDebug(2, Form("Using Detector %s\n", SortedDetectorName(idet)));
    if(TMath::Abs(fDetectorPID[fSortedOrder[idet]]->IsSelected(track, pidqa)) != 11){
      isSelected = kFALSE;
      break;
    }
    AliDebug(2, "Particlae selected by detector");
    if(fVarManager && cont){
      TString reccontname = contname; reccontname += "Reco";
      AliDebug(2, Form("Filling container %s", reccontname.Data()));
      if(fVarManager->IsSignalTrack())
        fVarManager->FillContainerStepname(cont, reccontname.Data(), SortedDetectorName(idet));
      if(HasMCData()){
        TString mccontname = contname; mccontname += "MC";
        AliDebug(2, Form("MC Information available, Filling container %s", mccontname.Data()));
        if(fVarManager->IsSignalTrack()) {
          fVarManager->FillContainerStepname(cont, mccontname.Data(), SortedDetectorName(idet), kTRUE);
	        if(cont->GetCorrelationMatrix("correlationstepafterTOF")){
	          TString tstept("TOFPID"); 
	          if(!tstept.CompareTo(SortedDetectorName(idet))) {
	            fVarManager->FillCorrelationMatrix(cont->GetCorrelationMatrix("correlationstepafterTOF"));
	            //printf("Step %s\n",(const char*) SortedDetectorName(idet));
	          }
	        }
	      }
      }
      // The PID will NOT fill the double counting information
    }
  }
  return isSelected;
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
void AliHFEpid::SetAODpid(AliAODpidUtil *pid){
  //
  // Set ESD PID to the Detector PID objects
  //
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]) fDetectorPID[idet]->SetAODpid(pid);
  }
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCasymmetric(Double_t pmin, Double_t pmax, Double_t sigmamin, Double_t sigmamax){
  //
  // TPC alone, symmetric 3 sigma cut and asymmetric sigma cut in the momentum region between 2GeV/c and 10 GeV/c and sigma between -1 and 100
  //
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  if(pid){
    pid->SetTPCnSigma(3);
    pid->SetAsymmetricTPCsigmaCut(pmin, pmax, sigmamin, sigmamax);
  }
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCrejectionSimple(){
  //
  // TPC alone, symmetric 3 sigma cut and 2 - -100 sigma pion rejection
  //   
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  if(pid){
    pid->SetTPCnSigma(3);
    pid->SetRejectParticle(AliPID::kPion, 0., -100., 10., 1.);
  }
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCrejection(const char *lowerCutParam, Double_t * const params, Float_t upperTPCCut, Float_t TOFCut){
  //
  // Combined TPC-TOF PID
  // if no function parameterizaion is given, then the default one (exponential) is chosen
  //
  if(HasMCData()) AliInfo("Configuring TPC for MC\n");
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fDetectorPID[kTOFpid]);
  if(tofpid) tofpid->SetTOFnSigma(TOFCut);

  //TF1 *upperCut = new TF1("upperCut", "[0] * TMath::Exp([1]*x)", 0, 20);
  TF1 *upperCut = new TF1("upperCut", "[0]", 0, 20); // Use constant upper cut
  TF1 *lowerCut = new TF1("lowerCut", lowerCutParam == NULL ? "[0] * TMath::Exp([1]*x) + [2]": lowerCutParam, 0, 20);

  upperCut->SetParameter(0, upperTPCCut); // pp

//  upperCut->SetParameter(0, 3.5); //PbPb
//  printf("upper %f \n",upperCut->GetParameter(0));
  //upperCut->SetParameter(0, 2.7);
  //upperCut->SetParameter(1, -0.4357);
  if(params){
    for(Int_t ipar = 0; ipar < lowerCut->GetNpar(); ipar++) lowerCut->SetParameter(ipar, params[ipar]);
  } else {
    // Set default parameterization
      if(HasMCData()) lowerCut->SetParameter(0, -2.5);
      else lowerCut->SetParameter(0, -4.03);  //pp
//      else lowerCut->SetParameter(0, -3.83);  //PbPb
  
     
    //else lowerCut->SetParameter(0, -3.71769);
    //else lowerCut->SetParameter(0, -3.7);

    lowerCut->SetParameter(1, -0.22); // pp
//     lowerCut->SetParameter(1,-0.36 ); // PbPb
    //lowerCut->SetParameter(1, -0.40263);
    //lowerCut->SetParameter(1, -0.8);
    
    if(HasMCData()) lowerCut->SetParameter(2, -2.2);
    else lowerCut->SetParameter(2, 0.92); //pp
//     else lowerCut->SetParameter(2, 0.27); //PbPb
   // else lowerCut->SetParameter(2, 0.92); //pp
//     printf("lower0 %f \n",lowerCut->GetParameter(0));
//     printf("lower1 %f \n",lowerCut->GetParameter(1));
//     printf("lower2 %f \n",lowerCut->GetParameter(2));
    //else lowerCut->SetParameter(2, 0.267857);
    //else lowerCut->SetParameter(2, -0.35);
  }


  if(tpcpid){
    tpcpid->SetTPCnSigma(2);
    tpcpid->SetUpperSigmaCut(upperCut);
    tpcpid->SetLowerSigmaCut(lowerCut);
  }
  AddCommonObject(upperCut);
  AddCommonObject(lowerCut);
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCstrategyParis(){
  //
  // TPC alone, symmetric 3 sigma cut and 2 - -100 sigma pion rejection
  //   
  AliHFEpidTPC *pid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  if(pid){
    pid->SetTPCnSigma(2);
    pid->SetRejectParticle(AliPID::kProton, 0., -3., 10., 3.);
    pid->SetRejectParticle(AliPID::kKaon, 0., -3., 10., 3.);
  }
}

//____________________________________________________________
void AliHFEpid::PrintStatus() const {
  //
  // Print the PID configuration
  //
  printf("\n%s: Printing configuration\n", GetName());
  printf("===============================================\n");
  printf("PID Detectors: \n");
  Int_t npid = 0;
  TString detectors[kNdetectorPID] = {"MC", "ESD", "ITS", "TPC", "TRD", "TOF", "EMCAL"};
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(IsDetectorOn(idet)){
      printf("\t%s\n", detectors[idet].Data());
      npid++;
    }
  }
  if(!npid) printf("\tNone\n");
  printf("\n");
}
