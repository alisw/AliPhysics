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
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

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
#include "AliHFEpidBayes.h"
#include "AliHFEvarManager.h"

ClassImp(AliHFEpid)

const Char_t* AliHFEpid::fgkDetectorName[AliHFEpid::kNdetectorPID + 1] = {
  "MCPID",
  "BAYESPID",
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
  memset(fSortedOrder, 0, sizeof(UInt_t) * kNdetectorPID);
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
  fDetectorPID[kBAYESpid] = new AliHFEpidBayes("BAYESPID");
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
  else if(!detector.CompareTo("BAYES")) detectorID = kBAYESpid;
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
Bool_t AliHFEpid::InitializePID(Int_t run){
  //
  // Initializes the PID object
  //

  if(!TestBit(kDetectorsSorted)) SortDetectors();
  Bool_t status = kTRUE;
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(!IsDetectorOn(idet)) continue;
    if(fDetectorPID[idet]){ 
      status &= fDetectorPID[idet]->InitializePID(run);
      if(HasMCData() && status) fDetectorPID[idet]->SetHasMCData();
    }
  }
  SetBit(kIsInit);
  PrintStatus();
  return status;
}

//____________________________________________________________
Bool_t AliHFEpid::IsSelected(const AliHFEpidObject * const track, AliHFEcontainer *cont, const Char_t *contname, AliHFEpidQAmanager *pidqa){
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
void AliHFEpid::SortDetectors(){
  //
  // Make sorted list of detectors
  //
  if(TestBit(kDetectorsSorted)) return; // Don't sort detectors when they are already sorted
  TMath::Sort(static_cast<UInt_t>(kNdetectorPID), fDetectorOrder, fSortedOrder, kFALSE);
  SetBit(kDetectorsSorted);
}

//____________________________________________________________
void AliHFEpid::SetPIDResponse(const AliPIDResponse * const pid){
  //
  // Set ESD PID to the Detector PID objects
  //
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]) fDetectorPID[idet]->SetPIDResponse(pid);
  }
}

//____________________________________________________________
const AliPIDResponse *AliHFEpid::GetPIDResponse() const {
  //
  // Return PID response function
  //
  const AliPIDResponse *response = NULL;
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]){
      response = fDetectorPID[idet]->GetPIDResponse();
      break;
    }
  } 
  return response;
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
void AliHFEpid::ConfigureTOF(Float_t TOFCut){
  //
  // Set Number of sigmas for TOF PID
  //
  AliHFEpidTOF *tofpid = dynamic_cast<AliHFEpidTOF *>(fDetectorPID[kTOFpid]);
  if(tofpid) tofpid->SetTOFnSigma(TOFCut);
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCcentralityCut(Int_t centralityBin, const char *lowerCutParam, const Double_t * const params, Float_t upperTPCCut){
  //
  // Cofigure centrality dependent cut function for TPC PID 
  //
  ConfigureTPCcut(centralityBin, lowerCutParam, params, upperTPCCut);
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCdefaultCut(const char *lowerCutParam, const Double_t * const params, Float_t upperTPCCut){
  //
  // Cofigure default cut function for TPC PID 
  //
  ConfigureTPCcut(-1, lowerCutParam, params, upperTPCCut);
}

//____________________________________________________________
void AliHFEpid::ConfigureTPCcut(Int_t centralityBin, const char *lowerCutParam, const Double_t * const params, Float_t upperTPCCut){
  //
  // Cofigure cut function for TPC PID 
  // if no function parameterizaion is given, then the default one (exponential) is chosen
  //
  
  if(HasMCData()) AliInfo("Configuring TPC for MC\n");
  AliHFEpidTPC *tpcpid = dynamic_cast<AliHFEpidTPC *>(fDetectorPID[kTPCpid]);
  //TF1 *upperCut = new TF1("upperCut", "[0] * TMath::Exp([1]*x)", 0, 20);
  TF1 *upperCut = new TF1(Form("upperCut%s", centralityBin <  0 ? "Default" : Form("Bin%d", centralityBin)), "[0]", 0, 20); // Use constant upper cut
  TF1 *lowerCut = new TF1(Form("lowerCut%s", centralityBin <  0 ? "Default" : Form("Bin%d", centralityBin)), lowerCutParam == NULL ? "[0] * TMath::Exp([1]*x) + [2]": lowerCutParam, 0, 20);

  upperCut->SetParameter(0, upperTPCCut); // pp

  if(params){
      for(Int_t ipar = 0; ipar < lowerCut->GetNpar(); ipar++)
      {
	  lowerCut->SetParameter(ipar, params[ipar]);
        //  printf("printout %i %s %f \n", centralityBin, lowerCutParam, params[ipar]);
      }
  } else {
    // Set default parameterization
    if(HasMCData()) lowerCut->SetParameter(0, -2.5);
    else lowerCut->SetParameter(0, -4.03);  //pp
    lowerCut->SetParameter(1, -0.22); // pp
    
    if(HasMCData()) lowerCut->SetParameter(2, -2.2);
    else lowerCut->SetParameter(2, 0.92); //pp
  }


  if(tpcpid){
    tpcpid->SetTPCnSigma(2);
    if(centralityBin < 0){
      tpcpid->SetUpperSigmaCutDefault(upperCut);
      tpcpid->SetLowerSigmaCutDefault(lowerCut);
    } else {
      tpcpid->SetUpperSigmaCutCentrality(upperCut, centralityBin);
      tpcpid->SetLowerSigmaCutCentrality(lowerCut, centralityBin);
    }
  }
  AddCommonObject(upperCut);
  AddCommonObject(lowerCut);
}

//____________________________________________________________
void AliHFEpid::ConfigureBayesDetectorMask(Int_t detmask){
  //
  // Configure detector mask for Bayes PID
  // if no detector mask is set the default mask is chosen
  //
  
  if(HasMCData()) AliInfo("Configuring Bayes for MC\n");
  AliHFEpidBayes *bayespid = dynamic_cast<AliHFEpidBayes *>(fDetectorPID[kBAYESpid]);

  if(bayespid)
  {
      bayespid->SetBayesDetectorMask(detmask);
      printf("detector mask in pid class %i \n",detmask);
  }

}

//____________________________________________________________
void AliHFEpid::ConfigureBayesPIDThreshold(Float_t pidthres){
  //
  // Configure pid threshold for Bayes PID
  // if no threshold is set the default threshold is chosen
  //
  
  if(HasMCData()) AliInfo("Configuring Bayes for MC\n");
  AliHFEpidBayes *bayespid = dynamic_cast<AliHFEpidBayes *>(fDetectorPID[kBAYESpid]);

  if(bayespid)
  {
      bayespid->SetBayesPIDThreshold(pidthres);
      printf("combined pidthreshold %f \n",pidthres);
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
  TString detectors[kNdetectorPID] = {"MC", "BAYES", "ITS", "TPC", "TRD", "TOF", "EMCAL"};
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(IsDetectorOn(idet)){
      printf("\t%s\n", detectors[idet].Data());
      npid++;
    }
  }
  if(!npid) printf("\tNone\n");
  printf("\n");
}
