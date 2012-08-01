/***************************************************************************
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

/*
 $Id: AliITSDetTypeSimUpg.cxx 56993 2012-06-08 08:24:17Z fca $
*/


/////////////////////////////////////////////////////////////////////
// Base simulation functions for ITS upgrade                       //
//                                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////          
#include "TBranch.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TTree.h"

#include "AliRun.h"

#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliITSdigit.h"
#include "AliITSdigitPixUpg.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeSimUpg.h"
#include "AliITSpListItem.h"
#include "AliITSCalibration.h"
#include "AliITSsimulation.h"
#include "AliITSTriggerConditions.h"
#include "AliITSsegmentation.h"
#include "AliITSsegmentationPixUpg.h"
#include "AliBaseLoader.h"

using std::endl;
using std::cout;
ClassImp(AliITSDetTypeSimUpg)

const char* AliITSDetTypeSimUpg::fgkDetTypeName[AliITSDetTypeSimUpg::kNDetTypes] = {"PixUpg"};


//----------------------------------------------------------------------
AliITSDetTypeSimUpg::AliITSDetTypeSimUpg():
TObject(),
fSimulation(),   // [NDet]
fSegmentation(), // [NDet]
fCalibration(),     // [NMod]
fNSDigits(0),    //! number of SDigits
fSDigits("AliITSpListItem",1000),   
fNDigits(0),     //! number of Digits
fRunNumber(0),   //! Run number (to access DB)
fDigits(),       //! [NMod][NDigits]
fSimuPar(0),
fkDigClassName(), // String with digit class name.
fLoader(0),      // local pointer to loader
fFirstcall(kTRUE),
fTriggerConditions(NULL)
{ 
  // Default Constructor
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    A properly zero-ed AliITSDetTypeSimUpg class.
  for (int i=kNDetTypes;i--;) fDetNModules[i] = 0;
  fSimulation   = new TObjArray(kNDetTypes);
  fSegmentation = new TObjArray(kNDetTypes);
  fSegmentation->SetOwner(kTRUE);
  fDigits = new TObjArray(kNDetTypes);
  fNDigits = new Int_t[kNDetTypes];
  fSimuPar= new AliITSSimuParam();
  SetRunNumber();
}

//----------------------------------------------------------------------
AliITSDetTypeSimUpg::~AliITSDetTypeSimUpg()
{
  // Destructor
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    Nothing.
  
  if(fSimulation){
    fSimulation->Delete();
    delete fSimulation;
  }
  //
  if(fSegmentation){
    fSegmentation->Delete();
    delete fSegmentation;
  }
  //
  if(fCalibration && fRunNumber<0){
    fCalibration->Delete();
    delete fCalibration;
  }
  //
  delete fSimuPar;
  delete[] fNDigits;
  //
  if (fLoader)fLoader->GetModulesFolder()->Remove(this); // Not deleting it.
  fSDigits.Delete();
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }
  //
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::SetITSgeom(AliITSgeom *geom)
{
  // Sets/replaces the existing AliITSgeom object kept in AliITSLoader
    // 
    // Inputs:
    //   AliITSgoem   *geom  The AliITSgeom object to be used.
    // Output:
    //   none.
    // Return:
    //   none.
  if(!fLoader){
    Error("SetITSgeom","No pointer to loader - nothing done");
    return;
  }
  else {
    fLoader->SetITSgeom(geom);  // protections in AliITSLoader::SetITSgeom
  }
 
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::SetLoader(AliITSLoader *loader)
{
  // Sets the local copy of the AliITSLoader, and passes on the
  // AliITSgeom object as needed.
  // Inputs
  //   AliITSLoader  *loader pointer to AliITSLoader for local use
  // Outputs:
  //   none.
  // Return:
  //  none.
  
  if(fLoader==loader) return; // Same do nothing
  if(fLoader){ // alread have an existing loader
    Error("SetLoader","Already have an exisiting loader ptr=%p Nothing done",fLoader);
  } // end if
  fLoader = loader;
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::SetSimulationModel(Int_t dettype,AliITSsimulation *sim)
{
  //Set simulation model for detector type
  if(fSimulation==0) fSimulation = new TObjArray(kNDetTypes);
  fSimulation->AddAt(sim,dettype);
}

//______________________________________________________________________
AliITSsimulation* AliITSDetTypeSimUpg::GetSimulationModel(Int_t dettype) const 
{ 
  //Get simulation model for detector type
  if(fSimulation==0)  {
    Warning("GetSimulationModel","fSimulation is 0!");
    return 0;     
  }
  return (AliITSsimulation*)(fSimulation->At(dettype));
}

//______________________________________________________________________
AliITSsimulation* AliITSDetTypeSimUpg::GetSimulationModelByModule(Int_t module) const 
{
  //Get simulation model by module number
  if(GetITSgeom()==0) {
    Warning("GetSimulationModelByModule","GetITSgeom() is 0!");
    return 0;
  }
  
  return GetSimulationModel(GetITSgeom()->GetModuleType(module));
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaultSegmentation(Int_t idet)
{
  // Set default segmentation model objects
  AliITSsegmentation *seg;
  //  
  if(fSegmentation==0x0){
    fSegmentation = new TObjArray(kNDetTypes);
    fSegmentation->SetOwner(kTRUE);
  }
  if (GetSegmentationModel(idet)) delete (AliITSsegmentation*)fSegmentation->RemoveAt(idet);
  //
  switch (idet) 
    {
    case kDetPixUpg: seg = new AliITSsegmentationPixUpg(); break;
    default        : AliFatal(Form("Uknown detector type %d",idet));
  };
  SetSegmentationModel(idet,seg);
  //
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::SetSegmentationModel(Int_t dettype, AliITSsegmentation *seg)
{
  // Set segmentation model for detector type
  if(fSegmentation==0x0){
    fSegmentation = new TObjArray(kNDetTypes);
    fSegmentation->SetOwner(kTRUE);
  }
  fSegmentation->AddAt(seg,dettype);
}

//______________________________________________________________________
AliITSsegmentation* AliITSDetTypeSimUpg::GetSegmentationModel(Int_t dettype) const
{
  //Get segmentation model for detector type
  
  if(fSegmentation==0) {
    Warning("GetSegmentationModel","fSegmentation is 0!");
    return 0; 
  } 
  return (AliITSsegmentation*)(fSegmentation->At(dettype));
}

//_______________________________________________________________________
AliITSsegmentation* AliITSDetTypeSimUpg::GetSegmentationModelByModule(Int_t module) const
{
  //Get segmentation model by module number
  if(GetITSgeom()==0) {
    Warning("GetSegmentationModelByModule","GetITSgeom() is 0!");
    return 0;
  }     
  return GetSegmentationModel(GetITSgeom()->GetModuleType(module));
  //
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::CreateCalibrationArray() 
{
  //Create the container of calibration functions with correct size
  if (fCalibration) {
    Warning("CreateCalibration","pointer to calibration object exists\n");
    fCalibration->Delete();
    delete fCalibration;
  }
  //
  Int_t nModTot = GetITSgeom()->GetIndexMax();
  fCalibration = new TObjArray(nModTot);
  fCalibration->SetOwner(kTRUE);
  fCalibration->Clear();
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetCalibrationModel(Int_t iMod, AliITSCalibration *resp)
{
  //Set response model for modules
  
  if (fCalibration==0) CreateCalibrationArray();
  
  if (fCalibration->At(iMod)!=0) delete (AliITSCalibration*) fCalibration->At(iMod);
  fCalibration->AddAt(resp, iMod);
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::ResetCalibrationArray() 
{
  //resets response array
  fCalibration->Clear();
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::ResetSegmentation()
{
  //Resets segmentation array
  if(fSegmentation) fSegmentation->Clear();
}

//_______________________________________________________________________
AliITSCalibration* AliITSDetTypeSimUpg::GetCalibrationModel(Int_t iMod) const 
{
  //Get response model for module number iMod 
  //
  if(fCalibration==0) {
    AliError("fCalibration is 0!");
    return 0; 
  }
  return (AliITSCalibration*)fCalibration->At(iMod);
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaults()
{
  //Set defaults for segmentation and response
  
  if(GetITSgeom()==0){
    Warning("SetDefaults","GetITSgeom() is 0!");
    return;
  } // end if
  if (fCalibration==0) {
    CreateCalibrationArray();
  } // end if
  
  ResetSegmentation();
  if(!GetCalibration()){AliFatal("Exit"); exit(0);}
  
  SetDigitClassName(0,"AliITSdigitPixUpg");
  //  
  for(Int_t idet=0;idet<kNDetTypes;idet++){
    if(!GetSegmentationModel(idet)) SetDefaultSegmentation(idet);
  }
  //
}

//______________________________________________________________________
Bool_t AliITSDetTypeSimUpg::GetCalibration() {
  // Get Default calibration if a storage is not defined.

  if(!fFirstcall){
    AliITSCalibration* cal = GetCalibrationModel(0);
    if(cal)return kTRUE;
  }
  else {
    fFirstcall = kFALSE;
  }
  
  SetRunNumber((Int_t)AliCDBManager::Instance()->GetRun());
  Int_t run=GetRunNumber();
  
  Bool_t isCacheActive = kTRUE;
  if (run) {
    isCacheActive=kFALSE;
    fCalibration->SetOwner(kTRUE);
  }
  else{
    fCalibration->SetOwner(kFALSE);
  }
  
  AliCDBManager::Instance()->SetCacheFlag(isCacheActive);
  //
  // TO DO, RS
  return kTRUE;
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaultSimulation()
{
  //Set default simulation for detector type
  //
  if (GetITSgeom()==0) {
    Warning("SetDefaultSimulation","GetITSgeom() is 0!");
    return;
  }
  if (fCalibration==0) {
    Warning("SetDefaultSimulation","fCalibration is 0!");
    return;
  }
  if (fSegmentation==0) {
    Warning("SetDefaultSimulation","fSegmentation is 0!");
    for (Int_t i=0;i<kNDetTypes;i++) SetDefaultSegmentation(i);
  } else for(Int_t i=0;i<kNDetTypes;i++) if(!GetSegmentationModel(i)){
	Warning("SetDefaultSimulation",
		"Segmentation not defined for det %d - Default taken\n!",i);
	SetDefaultSegmentation(i);
      }
  //
  AliITSsimulation* sim;
  //
  for (Int_t idet=0;idet<kNDetTypes;idet++) {
    //OixUpg
    if(idet==0){
      sim = GetSimulationModel(idet); 
      if(!sim){
	sim = new AliITSsimulationPixUpg(this);
	SetSimulationModel(idet,sim);
      }
    }
  }
}

//___________________________________________________________________
void AliITSDetTypeSimUpg::SetTreeAddressS(TTree* treeS, const Char_t* name)
{
  // Set branch address for the ITS summable digits Trees.  

  if(!treeS) return;
  TBranch *branch;
  branch = treeS->GetBranch(name);
  TClonesArray *sdigi = &fSDigits;
  if (branch) branch->SetAddress(&sdigi);

}

//___________________________________________________________________
void AliITSDetTypeSimUpg::SetTreeAddressD(TTree* treeD, const Char_t* name)
{
  // Set branch address for the digit Trees.
  TBranch *branch; 
  TString branchname;
  //
  if(!treeD) return;
  //
  if(!fDigits) fDigits = new TObjArray(kNDetTypes); 
  //
  for(Int_t i=0;i<kNDetTypes;i++){
    const Char_t* digclass = GetDigitClassName(i);
    if(digclass==0x0){
      if(i==0) SetDigitClassName(i,Form("AliITSdigit%s",fgkDetTypeName[i]));
      digclass = GetDigitClassName(i);
    }
    TString classn = digclass;
    if(!(fDigits->At(i))){
      fDigits->AddAt(new TClonesArray(classn.Data(),1000),i);
    }else{
      ResetDigits(i);
    }
    //
    if(kNDetTypes==3) branchname.Form("%sDigits%s",name,fgkDetTypeName[i]);
    branch = treeD->GetBranch(branchname.Data());
    if(branch) branch->SetAddress(&((*fDigits)[i]));    
  }
}

//___________________________________________________________________
void AliITSDetTypeSimUpg::ResetDigits()
{
  // Reset number of digits and the digits array for the ITS detector.  
  //
  if(!fDigits){
    Error("ResetDigits","fDigits is null!");
    return;
  }
  for(Int_t i=0;i<kNDetTypes;i++) ResetDigits(i);
  //
}

//___________________________________________________________________
void AliITSDetTypeSimUpg::ResetDigits(Int_t branch)
{
  // Reset number of digits and the digits array for this branch.
  //
  if (fDigits->At(branch)) ((TClonesArray*)fDigits->At(branch))->Clear();
  if(fNDigits) fNDigits[branch]=0;
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SDigitsToDigits(Option_t* opt, Char_t* name)
{
  // Standard Summable digits to Digits function.
  if(!GetITSgeom()){
    Warning("SDigitsToDigits","GetITSgeom() is null!!");
    return;
  }
  //
  static Bool_t setDef = kTRUE;
  if(setDef) SetDefaultSimulation();
  setDef = kFALSE;
  //
  AliITSsimulation *sim =0;
  TTree* trees = fLoader->TreeS();
  if( !(trees && GetSDigits()) ){
    Error("SDigits2Digits","Error: No trees or SDigits. Returning.");
    return;
  } 
  //
  TBranch* brchSDigits = trees->GetBranch(name);
  //
  Int_t id;
  for(Int_t module=0;module<GetITSgeom()->GetIndexMax();module++) {
    //
    id = GetITSgeom()->GetModuleType(module);
    sim = (AliITSsimulation*)GetSimulationModel(id);
    if(!sim){
      AliFatal(Form("The simulation class was not instanciated for module %d type %s!",
		    module,GetITSgeom()->GetModuleTypeName(module)));
    }
    sim->InitSimulationModule(module,gAlice->GetEvNumber());
    //
    fSDigits.Clear();
    brchSDigits->GetEvent(module);
    sim->AddSDigitsToModule(&fSDigits,0);
    sim->FinishSDigitiseModule();
    fLoader->TreeD()->Fill();
    ResetDigits();
  }
  //
  //  WriteFOSignals(); 
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  fLoader->TreeD()->Reset();
}

//_________________________________________________________
void AliITSDetTypeSimUpg::AddSumDigit(AliITSpListItem &sdig)
{  
  //Adds the module full of summable digits to the summable digits tree.
  new(fSDigits[fNSDigits++]) AliITSpListItem(sdig);
}

//__________________________________________________________
void AliITSDetTypeSimUpg::AddSimDigit(Int_t branch, const AliITSdigit* d)
{  
  // Add a simulated digit.
  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  switch(branch){
  case kDetPixUpg:
    new(ldigits[fNDigits[branch]++]) AliITSdigitPixUpg(*((AliITSdigitPixUpg*)d));
    break;
  default:
    AliFatal(Form("Digit for unknown detector type %d",branch));
  }  
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::AddSimDigit(Int_t branch,Float_t phys,Int_t *digits,
				      Int_t *tracks,Int_t *hits,Float_t *charges, 
				      Int_t sigexpanded)
{
  //   Add a simulated digit to the list.

  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  switch(branch){
  case kDetPixUpg:
    new(ldigits[fNDigits[branch]++]) AliITSdigitPixUpg(digits,tracks,hits);
    break;
  default:
    AliFatal(Form("Digit for unknown detector type %d",branch));
  } 
}

//_______________________________________________________________________
AliITSTriggerConditions* AliITSDetTypeSimUpg::GetTriggerConditions() 
{
  // Get Pixel Trigger Conditions (separate method since it is used only when simulating trigger)
  if (fTriggerConditions==NULL) { // read from db
    fRunNumber = ((Int_t)AliCDBManager::Instance()->GetRun());
    Bool_t origCacheStatus = AliCDBManager::Instance()->GetCacheFlag();
    Bool_t isCacheActive;
    if (fRunNumber<0) isCacheActive=kFALSE;
    else              isCacheActive=kTRUE;
    AliCDBManager::Instance()->SetCacheFlag(isCacheActive);
    AliCDBEntry *pitCond = AliCDBManager::Instance()->Get("TRIGGER/SPD/PITConditions", fRunNumber);
    if (!pitCond) {
      AliError("Trigger conditions retrieval failed! ");
      return NULL;
    }
    fTriggerConditions = (AliITSTriggerConditions*) pitCond->GetObject();
    if (!isCacheActive) pitCond->SetObject(NULL);
    pitCond->SetOwner(kTRUE);
    if (!isCacheActive) {
      delete pitCond;
    }
    AliCDBManager::Instance()->SetCacheFlag(origCacheStatus);
    if (fTriggerConditions==NULL) {
      AliWarning("fTriggerConditions is NULL!");
    }
  }
  return fTriggerConditions;
}

