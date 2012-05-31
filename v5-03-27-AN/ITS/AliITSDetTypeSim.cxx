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
 $Id$
*/


/////////////////////////////////////////////////////////////////////
// Base simulation functions for ITS                               //
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
#include "AliITSdigitSPD.h"
#include "AliITSdigitSDD.h"
#include "AliITSdigitSSD.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeSim.h"
#include "AliITSpListItem.h"
#include "AliITSCalibration.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSMapSDD.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSGainSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSNoiseSSD.h"
#include "AliITSGainSSD.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulation.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSresponseSDD.h"
#include "AliITSDDLModuleMapSDD.h"
#include "AliITSTriggerConditions.h"
#include "AliBaseLoader.h"

const Int_t AliITSDetTypeSim::fgkNdettypes = 3;
const Int_t AliITSDetTypeSim::fgkDefaultNModulesSPD =  240;
const Int_t AliITSDetTypeSim::fgkDefaultNModulesSDD =  260;
const Int_t AliITSDetTypeSim::fgkDefaultNModulesSSD = 1698;

ClassImp(AliITSDetTypeSim)

//----------------------------------------------------------------------
AliITSDetTypeSim::AliITSDetTypeSim():
TObject(),
fSimulation(),   // [NDet]
fSegmentation(), // [NDet]
fCalibration(),     // [NMod]
fSSDCalibration(0),
fSPDNoisy(0),
fSPDSparseDead(0),
fNSDigits(0),    //! number of SDigits
fSDigits("AliITSpListItem",1000),   
fNDigits(0),     //! number of Digits
fRunNumber(0),   //! Run number (to access DB)
fDigits(),       //! [NMod][NDigits]
fSimuPar(0),
fDDLMapSDD(0),
fRespSDD(0),
fAveGainSDD(0),
fkDigClassName(), // String with digit class name.
fLoader(0),      // local pointer to loader
fFirstcall(kTRUE),
fFOGenerator(),
fTriggerConditions(NULL)
{ 
    // Default Constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A properly zero-ed AliITSDetTypeSim class.

  fSimulation = new TObjArray(fgkNdettypes);
  fSegmentation = new TObjArray(fgkNdettypes);
  fSegmentation->SetOwner(kTRUE);
  fDigits = new TObjArray(fgkNdettypes);
  fNDigits = new Int_t[fgkNdettypes];
  fDDLMapSDD=new AliITSDDLModuleMapSDD();
  fSimuPar= new AliITSSimuParam();
  fSSDCalibration=new AliITSCalibrationSSD();
  SetRunNumber();
}
//----------------------------------------------------------------------
AliITSDetTypeSim::~AliITSDetTypeSim(){
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
    fSimulation = 0;
    if(fSegmentation){
	fSegmentation->Delete();
	delete fSegmentation;
    }
    fSegmentation = 0;
    if(fCalibration && fRunNumber<0){

	fCalibration->Delete();
	delete fCalibration;
    }
    fCalibration = 0;
    if(fSSDCalibration) {
      if(!(AliCDBManager::Instance()->GetCacheFlag())) {
	delete fSSDCalibration;
	fSSDCalibration = NULL;
      }
    }
    if(fSPDNoisy){
    if(!(AliCDBManager::Instance()->GetCacheFlag())) {
       fSPDNoisy->Delete();
       delete fSPDNoisy;
       fSPDNoisy = 0;
      }
    }
    if(fSPDSparseDead){
    if(!(AliCDBManager::Instance()->GetCacheFlag())) {
       fSPDSparseDead->Delete();
       delete fSPDSparseDead;
       fSPDSparseDead = 0;
      }
    }
    if(fSimuPar) delete fSimuPar;
    if(fRespSDD){
      if(!(AliCDBManager::Instance()->GetCacheFlag())){
	delete fRespSDD;
	fRespSDD=0;
      }
    }
    if(fDDLMapSDD) delete fDDLMapSDD;
    if(fNDigits) delete [] fNDigits;
    fNDigits = 0;
    if (fLoader)fLoader->GetModulesFolder()->Remove(this);
    fLoader = 0; // Not deleting it.
    fSDigits.Delete();
    if (fDigits) {
      fDigits->Delete();
      delete fDigits;
    }
    fDigits=0;
}
//----------------------------------------------------------------------
// AliITSDetTypeSim::AliITSDetTypeSim(const AliITSDetTypeSim &source) : TObject(source),
// fSimulation(source.fSimulation),   // [NDet]
// fSegmentation(source.fSegmentation), // [NDet]
// fCalibration(source.fCalibration),     // [NMod]
// fSSDCalibration(source.fSSDCalibration),
// fSPDNoisy(source.fSPDNoisy),
// fSPDSparseDead(source.fSPDSparseDead),
// fNSDigits(source.fNSDigits),    //! number of SDigits
// fSDigits(*((TClonesArray*)source.fSDigits.Clone())),
// fNDigits(source.fNDigits),     //! number of Digits
// fRunNumber(source.fRunNumber),   //! Run number (to access DB)
// fDigits(source.fDigits),       //! [NMod][NDigits]
// fSimuPar(source.fSimuPar),
// fDDLMapSDD(source.fDDLMapSDD),
// fRespSDD(source.fRespSDD),
// fAveGainSDD(source.fAveGainSDD),
// fkDigClassName(), // String with digit class name.
// fLoader(source.fLoader),      // local pointer to loader
// fFirstcall(source.fFirstcall),
// fFOGenerator(source.fFOGenerator),
// fTriggerConditions(source.fTriggerConditions) 
// {
//     // Copy Constructor for object AliITSDetTypeSim not allowed
//   for(Int_t i=0;i<fgkNdettypes;i++){
//     fkDigClassName[i] = source.fkDigClassName[i];
//   }
// }
// //----------------------------------------------------------------------
// AliITSDetTypeSim& AliITSDetTypeSim::operator=(const AliITSDetTypeSim &source){
//     // The = operator for object AliITSDetTypeSim
 
//   this->~AliITSDetTypeSim();
//   new(this) AliITSDetTypeSim(source);
//   return *this;
// }

//______________________________________________________________________
void AliITSDetTypeSim::SetITSgeom(AliITSgeom *geom){
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
void AliITSDetTypeSim::SetLoader(AliITSLoader *loader){
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
	Error("SetLoader",
		"Already have an exisiting loader ptr=%p Nothing done",
		fLoader);
    } // end if
    fLoader = loader;
}
//______________________________________________________________________
void AliITSDetTypeSim::SetSimulationModel(Int_t dettype,AliITSsimulation *sim){

  //Set simulation model for detector type

  if(fSimulation==0) fSimulation = new TObjArray(fgkNdettypes);
  fSimulation->AddAt(sim,dettype);
}
//______________________________________________________________________
AliITSsimulation* AliITSDetTypeSim::GetSimulationModel(Int_t dettype) const { 

  //Get simulation model for detector type
  if(fSimulation==0)  {
    Warning("GetSimulationModel","fSimulation is 0!");
    return 0;     
  }
  return (AliITSsimulation*)(fSimulation->At(dettype));
}
//______________________________________________________________________
AliITSsimulation* AliITSDetTypeSim::GetSimulationModelByModule(Int_t module) const {

  //Get simulation model by module number
  if(GetITSgeom()==0) {
    Warning("GetSimulationModelByModule","GetITSgeom() is 0!");
    return 0;
  }
  
  return GetSimulationModel(GetITSgeom()->GetModuleType(module));
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetDefaultSegmentation(Int_t idet){
    // Set default segmentation model objects
    AliITSsegmentation *seg;

    if(fSegmentation==0x0){
	fSegmentation = new TObjArray(fgkNdettypes);
	fSegmentation->SetOwner(kTRUE);
    }
    if(GetSegmentationModel(idet))
	delete (AliITSsegmentation*)fSegmentation->At(idet);
    if(idet==0){
	seg = new AliITSsegmentationSPD();
    }else if(idet==1){
      seg = new AliITSsegmentationSDD();
      AliITSCalibrationSDD* cal=(AliITSCalibrationSDD*)GetCalibrationModel(fgkDefaultNModulesSPD+1);
      if(cal->IsAMAt20MHz()){ 
	seg->SetPadSize(seg->Dpz(0),20.);
	seg->SetNPads(seg->Npz()/2,128);
      }
    }else {
	seg = new AliITSsegmentationSSD();
    }
    SetSegmentationModel(idet,seg);
}
//______________________________________________________________________
void AliITSDetTypeSim::SetSegmentationModel(Int_t dettype,
					    AliITSsegmentation *seg){
   
  //Set segmentation model for detector type
  if(fSegmentation==0x0){
    fSegmentation = new TObjArray(fgkNdettypes);
    fSegmentation->SetOwner(kTRUE);
  }
  fSegmentation->AddAt(seg,dettype);
}
//______________________________________________________________________
AliITSsegmentation* AliITSDetTypeSim::GetSegmentationModel(Int_t dettype) const{
  //Get segmentation model for detector type
   
   if(fSegmentation==0) {
       Warning("GetSegmentationModel","fSegmentation is 0!");
       return 0; 
   } 
   return (AliITSsegmentation*)(fSegmentation->At(dettype));
}
//_______________________________________________________________________
AliITSsegmentation* AliITSDetTypeSim::GetSegmentationModelByModule(Int_t module) const{
    //Get segmentation model by module number
    if(GetITSgeom()==0){
	Warning("GetSegmentationModelByModule","GetITSgeom() is 0!");
	return 0;
    }     
    return GetSegmentationModel(GetITSgeom()->GetModuleType(module));
}
//_______________________________________________________________________
void AliITSDetTypeSim::CreateCalibrationArray() {
    //Create the container of calibration functions with correct size
    if (fCalibration) {
	Warning("CreateCalibration","pointer to calibration object exists\n");
	fCalibration->Delete();
	delete fCalibration;
    }

    Int_t nModTot = GetITSgeom()->GetIndexMax();
    fCalibration = new TObjArray(nModTot);
    fCalibration->SetOwner(kTRUE);
    fCalibration->Clear();
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetCalibrationModel(Int_t iMod, AliITSCalibration *resp){
    //Set response model for modules

    if (fCalibration==0) CreateCalibrationArray();
 
    if (fCalibration->At(iMod)!=0)
	delete (AliITSCalibration*) fCalibration->At(iMod);
    fCalibration->AddAt(resp, iMod);
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetSPDNoisyModel(Int_t iMod, AliITSCalibration *cal){
  //Set noisy pixel info for the SPD module iMod
  if (fSPDNoisy==0) {
    fSPDNoisy = new TObjArray(fgkDefaultNModulesSPD);
    fSPDNoisy->SetOwner(kTRUE);
    fSPDNoisy->Clear();
  }

  if (fSPDNoisy->At(iMod) != 0)
    delete (AliITSCalibration*) fSPDNoisy->At(iMod);
  fSPDNoisy->AddAt(cal,iMod);
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetSPDSparseDeadModel(Int_t iMod, AliITSCalibration *cal){
  //Set sparse dead pixel info for the SPD module iMod
  if (fSPDSparseDead==0) {
    fSPDSparseDead = new TObjArray(fgkDefaultNModulesSPD);
    fSPDSparseDead->SetOwner(kTRUE);
    fSPDSparseDead->Clear();
  }

  if (fSPDNoisy->At(iMod) != 0)
    delete (AliITSCalibration*) fSPDNoisy->At(iMod);
  fSPDNoisy->AddAt(cal,iMod);
}
//______________________________________________________________________
void AliITSDetTypeSim::ResetCalibrationArray(){
    //resets response array
    if(fCalibration && fRunNumber<0){  // if fRunNumber<0 fCalibration is owner
      /*
	AliITSresponse* rspd = ((AliITSCalibration*)fCalibration->At(
                                GetITSgeom()->GetStartSPD()))->GetResponse();
	AliITSresponse* rssd = ((AliITSCalibration*)fCalibration->At(
                                GetITSgeom()->GetStartSSD()))->GetResponse();
	if(rspd) delete rspd;
	if(rssd) delete rssd;
      */
	fCalibration->Clear();
    }else if (fCalibration && fRunNumber>=0){
	fCalibration->Clear();
    }
}
//______________________________________________________________________
void AliITSDetTypeSim::ResetSegmentation(){
    //Resets segmentation array
    if(fSegmentation) fSegmentation->Clear();
}
//_______________________________________________________________________
AliITSCalibration* AliITSDetTypeSim::GetCalibrationModel(Int_t iMod) const {
    //Get response model for module number iMod 
 
    if(fCalibration==0) {
	AliError("fCalibration is 0!");
	return 0; 
    }
  if(iMod<fgkDefaultNModulesSPD+fgkDefaultNModulesSDD){
    return (AliITSCalibration*)fCalibration->At(iMod);
  }else{
    Int_t i=iMod-(fgkDefaultNModulesSPD+fgkDefaultNModulesSDD);
    fSSDCalibration->SetModule(i);
    return (AliITSCalibration*)fSSDCalibration;
  }

}
//_______________________________________________________________________
AliITSCalibration* AliITSDetTypeSim::GetSPDNoisyModel(Int_t iMod) const {
  //Get SPD noisy calib for module iMod 
  if(fSPDNoisy==0) {
    AliWarning("fSPDNoisy is 0!");
    return 0; 
  }
  return (AliITSCalibration*)fSPDNoisy->At(iMod);
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetDefaults(){
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

    SetDigitClassName(0,"AliITSdigitSPD");
    SetDigitClassName(1,"AliITSdigitSDD");
    SetDigitClassName(2,"AliITSdigitSSD");

    for(Int_t idet=0;idet<fgkNdettypes;idet++){
      if(!GetSegmentationModel(idet)) SetDefaultSegmentation(idet);
    }
}
//______________________________________________________________________
Bool_t AliITSDetTypeSim::GetCalibration() {
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

  Bool_t origCacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  Bool_t isCacheActive = kTRUE;
  if(GetRunNumber()<0){
    isCacheActive=kFALSE;
    fCalibration->SetOwner(kTRUE);
  }
  else{
    fCalibration->SetOwner(kFALSE);
  }

  AliCDBManager::Instance()->SetCacheFlag(isCacheActive);

  AliCDBEntry *deadSPD = AliCDBManager::Instance()->Get("ITS/Calib/SPDDead", run);
  AliCDBEntry *noisySPD = AliCDBManager::Instance()->Get("ITS/Calib/SPDNoisy", run);
  AliCDBEntry *foEffSPD = AliCDBManager::Instance()->Get("ITS/Calib/SPDFOEfficiency", run);
  AliCDBEntry *foNoiSPD = AliCDBManager::Instance()->Get("ITS/Calib/SPDFONoise", run);
  AliCDBEntry *entrySDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD", run);
  AliCDBEntry *entry2SDD = AliCDBManager::Instance()->Get("ITS/Calib/RespSDD");
  AliCDBEntry *drSpSDD = AliCDBManager::Instance()->Get("ITS/Calib/DriftSpeedSDD",run);
  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD",run);
  //AliCDBEntry *mapASDD = AliCDBManager::Instance()->Get("ITS/Calib/MapsAnodeSDD",run);
  AliCDBEntry *mapTSDD = AliCDBManager::Instance()->Get("ITS/Calib/MapsTimeSDD",run);
  // AliCDBEntry *entrySSD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSSD", run);
  AliCDBEntry *entryNoiseSSD = AliCDBManager::Instance()->Get("ITS/Calib/NoiseSSD");
  AliCDBEntry *entryGainSSD = AliCDBManager::Instance()->Get("ITS/Calib/GainSSD");
  AliCDBEntry *entryBadChannelsSSD = AliCDBManager::Instance()->Get("ITS/Calib/BadChannelsSSD");

if(!deadSPD || !noisySPD || !foEffSPD || !foNoiSPD 
     || !entrySDD  || !entry2SDD || !entryNoiseSSD || !entryGainSSD || !entryBadChannelsSSD 
     || !drSpSDD || !ddlMapSDD || !mapTSDD){
    AliFatal("Calibration object retrieval failed! ");
    return kFALSE;
  }  	  	
	

  TObjArray *calDeadSPD = (TObjArray*) deadSPD->GetObject();
  if (!isCacheActive) deadSPD->SetObject(NULL);
  deadSPD->SetOwner(kTRUE);

  TObjArray *calNoisySPD = (TObjArray*) noisySPD->GetObject();
  if (!isCacheActive) noisySPD->SetObject(NULL);
  noisySPD->SetOwner(kTRUE);

  AliITSFOEfficiencySPD *calFoEffSPD = (AliITSFOEfficiencySPD*) foEffSPD->GetObject();
  if (!isCacheActive) foEffSPD->SetObject(NULL);
  foEffSPD->SetOwner(kTRUE);

  AliITSFONoiseSPD *calFoNoiSPD = (AliITSFONoiseSPD*) foNoiSPD->GetObject();
  if (!isCacheActive) foNoiSPD->SetObject(NULL);
  foNoiSPD->SetOwner(kTRUE);
   
  TObjArray *calSDD = (TObjArray *)entrySDD->GetObject();
  if(!isCacheActive)entrySDD->SetObject(NULL);
  entrySDD->SetOwner(kTRUE);

  AliITSresponseSDD *pSDD = (AliITSresponseSDD*)entry2SDD->GetObject();
  if(!isCacheActive)entry2SDD->SetObject(NULL);
  entry2SDD->SetOwner(kTRUE);

  TObjArray *drSp = (TObjArray *)drSpSDD->GetObject();
  if(!isCacheActive)drSpSDD->SetObject(NULL);
  drSpSDD->SetOwner(kTRUE);

  AliITSDDLModuleMapSDD *ddlsdd=(AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  if(!isCacheActive)ddlMapSDD->SetObject(NULL);
  ddlMapSDD->SetOwner(kTRUE);

//   TObjArray *mapAn = (TObjArray *)mapASDD->GetObject();
//   if(!isCacheActive)mapASDD->SetObject(NULL);
//   mapASDD->SetOwner(kTRUE);

  TObjArray *mapT = (TObjArray *)mapTSDD->GetObject();
  if(!isCacheActive)mapTSDD->SetObject(NULL);
  mapTSDD->SetOwner(kTRUE);

  /*
  TObjArray *calSSD = (TObjArray *)entrySSD->GetObject();
  if(!isCacheActive)entrySSD->SetObject(NULL);
  entrySSD->SetOwner(kTRUE);
  */

  TObject *emptyssd = 0; TString ssdobjectname;
  AliITSNoiseSSDv2 *noiseSSD = NULL; 
  emptyssd = (TObject *)entryNoiseSSD->GetObject();
  ssdobjectname = emptyssd->GetName();
  if(ssdobjectname=="TObjArray") {
    noiseSSD = new AliITSNoiseSSDv2(); 
    TObjArray *noiseSSDOld = (TObjArray *)entryNoiseSSD->GetObject();
    ReadOldSSDNoise(noiseSSDOld, noiseSSD);
  }
  else if(ssdobjectname=="AliITSNoiseSSDv2")
    noiseSSD = (AliITSNoiseSSDv2 *)entryNoiseSSD->GetObject();
  if(!isCacheActive)entryNoiseSSD->SetObject(NULL);
  entryNoiseSSD->SetOwner(kTRUE);

  AliITSGainSSDv2 *gainSSD = NULL;
  emptyssd = (TObject *)entryGainSSD->GetObject();
  ssdobjectname = emptyssd->GetName();
  if(ssdobjectname=="Gain") {
    TObjArray *gainSSDOld = (TObjArray *)entryGainSSD->GetObject();
    gainSSD = new AliITSGainSSDv2();
    ReadOldSSDGain(gainSSDOld, gainSSD);
  }
  else if(ssdobjectname=="AliITSGainSSDv2")
    gainSSD = (AliITSGainSSDv2 *)entryGainSSD->GetObject();
  if(!isCacheActive)entryGainSSD->SetObject(NULL);
  entryGainSSD->SetOwner(kTRUE);

  AliITSBadChannelsSSDv2 *badChannelsSSD = NULL;
  emptyssd = (TObject *)entryBadChannelsSSD->GetObject();
  ssdobjectname = emptyssd->GetName();
  if(ssdobjectname=="TObjArray") {
    TObjArray *badChannelsSSDOld = (TObjArray *)entryBadChannelsSSD->GetObject();
    badChannelsSSD = new AliITSBadChannelsSSDv2();
    ReadOldSSDBadChannels(badChannelsSSDOld, badChannelsSSD);
  }
  else if(ssdobjectname=="AliITSBadChannelsSSDv2")
    badChannelsSSD = (AliITSBadChannelsSSDv2*)entryBadChannelsSSD->GetObject();
  if(!isCacheActive)entryBadChannelsSSD->SetObject(NULL);
  entryBadChannelsSSD->SetOwner(kTRUE);

  /*AliITSNoiseSSDv2 *noiseSSD = (AliITSNoiseSSDv2 *)entryNoiseSSD->GetObject();
  if(!isCacheActive)entryNoiseSSD->SetObject(NULL);
  entryNoiseSSD->SetOwner(kTRUE);

  AliITSGainSSDv2 *gainSSD = (AliITSGainSSDv2 *)entryGainSSD->GetObject();
  if(!isCacheActive)entryGainSSD->SetObject(NULL);
  entryGainSSD->SetOwner(kTRUE);

  AliITSBadChannelsSSDv2 *badchannelsSSD = 
    (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
  if(!isCacheActive)entryBadChannelsSSD->SetObject(NULL);
  entryBadChannelsSSD->SetOwner(kTRUE);*/

  // DB entries are deleted. In this way metadeta objects are deleted as well
  if(!isCacheActive){
    delete deadSPD;
    delete noisySPD;
    delete foEffSPD;
    delete foNoiSPD;
    delete entrySDD;
    delete entry2SDD;
    delete entryNoiseSSD;
    delete entryGainSSD;
    delete entryBadChannelsSSD;
//    delete mapASDD;   
    delete mapTSDD;
    delete drSpSDD;
    delete ddlMapSDD;
  }
  
  AliCDBManager::Instance()->SetCacheFlag(origCacheStatus);

 if ((!calDeadSPD) || (!calNoisySPD) || (!calFoEffSPD) || (!calFoNoiSPD) 
      || (!calSDD) || (!pSDD)|| (!drSp) || (!ddlsdd)  
      || (!mapT) || (!noiseSSD)|| (!gainSSD)|| (!badChannelsSSD)) {
    AliWarning("Can not get calibration from calibration database !");
    return kFALSE;
  }


  Int_t nmods0 = calDeadSPD->GetEntries();
  Int_t nmods1 = calSDD->GetEntries();

  AliDebug(1,Form("%i SPD, %i SDD and %i SSD in calibration database",
		  nmods0,nmods1,fgkDefaultNModulesSSD));
  AliITSCalibration* cal;
  for (Int_t i=0; i<nmods0; i++) {
    cal = (AliITSCalibration*) calDeadSPD->At(i);
    SetCalibrationModel(i, cal);
    cal = (AliITSCalibration*) calNoisySPD->At(i);
    SetSPDNoisyModel(i, cal);
  }

  fFOGenerator.SetEfficiency(calFoEffSPD); // this cal object is used only by the generator
  fFOGenerator.SetNoise(calFoNoiSPD); // this cal object is used only by the generator
  
  fDDLMapSDD->SetDDLMap(ddlsdd);
  fRespSDD=pSDD;
  Float_t avegain=0.;
  Float_t nGdAnodes=0;
  Bool_t oldMapFormat=kFALSE;
  TObject* objmap=(TObject*)mapT->At(0);
  TString cname(objmap->ClassName());
  if(cname.CompareTo("AliITSMapSDD")==0){ 
    oldMapFormat=kTRUE;
    AliInfo("SDD Maps converted to new format");
  }

  for (Int_t i=0; i<fgkDefaultNModulesSDD; i++) {
    Int_t iddl,icarlos;
    Int_t iMod = i + fgkDefaultNModulesSPD;
    fDDLMapSDD->FindInDDLMap(iMod,iddl,icarlos);
    if(iddl<0){ 
      AliITSCalibrationSDD* calsdddead=new AliITSCalibrationSDD();
      calsdddead->SetBad();      
      AliITSDriftSpeedSDD* driftspdef = new AliITSDriftSpeedSDD();
      AliITSDriftSpeedArraySDD* arrdrsp=new AliITSDriftSpeedArraySDD(1);
      arrdrsp->AddDriftSpeed(driftspdef);
      calsdddead->SetDriftSpeed(0,arrdrsp);
      calsdddead->SetDriftSpeed(1,arrdrsp);
      SetCalibrationModel(iMod, calsdddead);
      AliWarning(Form("SDD module %d not present in DDL map: set it as dead",iMod));
    }else{
      cal = (AliITSCalibration*) calSDD->At(i);
      for(Int_t iAnode=0;iAnode< ((AliITSCalibrationSDD*)cal)->NOfAnodes(); iAnode++){
	if(((AliITSCalibrationSDD*)cal)->IsBadChannel(iAnode)) continue;
	avegain+= ((AliITSCalibrationSDD*)cal)->GetChannelGain(iAnode);
	nGdAnodes++;
      }
      Int_t i0=2*i;
      Int_t i1=1+2*i;
      AliITSDriftSpeedArraySDD* arr0 = (AliITSDriftSpeedArraySDD*) drSp->At(i0);
      AliITSDriftSpeedArraySDD* arr1 = (AliITSDriftSpeedArraySDD*) drSp->At(i1);

      AliITSCorrMapSDD* mt0 = 0;
      AliITSCorrMapSDD* mt1 = 0;
      if(oldMapFormat){ 
	AliITSMapSDD* oldmap0=(AliITSMapSDD*)mapT->At(i0);
	AliITSMapSDD* oldmap1=(AliITSMapSDD*)mapT->At(i1);
	mt0=oldmap0->ConvertToNewFormat();
	mt1=oldmap1->ConvertToNewFormat();
      }else{
	mt0=(AliITSCorrMapSDD*)mapT->At(i0);
	mt1=(AliITSCorrMapSDD*)mapT->At(i1);
      }
      cal->SetDriftSpeed(0,arr0);
      cal->SetDriftSpeed(1,arr1);
      cal->SetMapT(0,mt0);
      cal->SetMapT(1,mt1);
      SetCalibrationModel(iMod, cal);
    }
  }
  if(nGdAnodes) fAveGainSDD=avegain/nGdAnodes;
  AliDebug(3,Form("SDD average gain=%f\n",fAveGainSDD));
  fSSDCalibration->SetNoise(noiseSSD);
  fSSDCalibration->SetGain(gainSSD);
  fSSDCalibration->SetBadChannels(badChannelsSSD);
  //fSSDCalibration->FillBadChipMap();



  return kTRUE;
}
//_______________________________________________________________________
void AliITSDetTypeSim::SetDefaultSimulation(){
  //Set default simulation for detector type

  if(GetITSgeom()==0){
    Warning("SetDefaultSimulation","GetITSgeom() is 0!");
    return;
  }
  if(fCalibration==0){
    Warning("SetDefaultSimulation","fCalibration is 0!");
    return;
  }
  if(fSegmentation==0){
    Warning("SetDefaultSimulation","fSegmentation is 0!");
    for(Int_t i=0;i<fgkNdettypes;i++) SetDefaultSegmentation(i);
  }else for(Int_t i=0;i<fgkNdettypes;i++) if(!GetSegmentationModel(i)){
      Warning("SetDefaultSimulation",
	      "Segmentation not defined for det %d - Default taken\n!",i);
      SetDefaultSegmentation(i);
  }
  AliITSsimulation* sim;

  for(Int_t idet=0;idet<fgkNdettypes;idet++){
   //SPD
    if(idet==0){
      sim = GetSimulationModel(idet); 
      if(!sim){
	sim = new AliITSsimulationSPD(this);
	SetSimulationModel(idet,sim);
      }
    }
    //SDD
    if(idet==1){
      sim = GetSimulationModel(idet);
      if(!sim){
	sim = new AliITSsimulationSDD(this);
	SetSimulationModel(idet,sim);
      }      
    }
    //SSD
    if(idet==2){
      sim = GetSimulationModel(idet);
      if(!sim){
	sim = new AliITSsimulationSSD(this);
	SetSimulationModel(idet,sim);
      }
    }

  }
}
//___________________________________________________________________
void AliITSDetTypeSim::SetTreeAddressS(TTree* treeS, const Char_t* name){
  // Set branch address for the ITS summable digits Trees.  

  if(!treeS){
    return;
  }
  TBranch *branch;
  branch = treeS->GetBranch(name);
  TClonesArray *sdigi = &fSDigits;
  if (branch) branch->SetAddress(&sdigi);

}
//___________________________________________________________________
void AliITSDetTypeSim::SetTreeAddressD(TTree* treeD, const Char_t* name){
  // Set branch address for the digit Trees.
  
  const char *det[3] = {"SPD","SDD","SSD"};
  TBranch *branch;
  
  TString branchname;
  
  if(!treeD){
    return;
  }
  if(!fDigits){
    fDigits = new TObjArray(fgkNdettypes); 
  }
  for(Int_t i=0;i<fgkNdettypes;i++){
    const Char_t* digclass = GetDigitClassName(i);
    if(digclass==0x0){
      if(i==0) SetDigitClassName(i,"AliITSdigitSPD");
      if(i==1) SetDigitClassName(i,"AliITSdigitSDD");
      if(i==2) SetDigitClassName(i,"AliITSdigitSSD");
      digclass = GetDigitClassName(i);
    }
    TString classn = digclass;
    if(!(fDigits->At(i))){
      fDigits->AddAt(new TClonesArray(classn.Data(),1000),i);
    }else{
      ResetDigits(i);
    }
    
    if(fgkNdettypes==3) branchname.Form("%sDigits%s",name,det[i]);
    else branchname.Form("%sDigits%d",name,i+1);
    branch = treeD->GetBranch(branchname.Data());
    if(branch) branch->SetAddress(&((*fDigits)[i]));    
  }
}
//___________________________________________________________________
void AliITSDetTypeSim::ResetDigits(){
  // Reset number of digits and the digits array for the ITS detector.  

  if(!fDigits){
    Error("ResetDigits","fDigits is null!");
    return;
  }
  for(Int_t i=0;i<fgkNdettypes;i++){
    ResetDigits(i);
  }
}
//___________________________________________________________________
void AliITSDetTypeSim::ResetDigits(Int_t branch){
  // Reset number of digits and the digits array for this branch.

  if(fDigits->At(branch)){
    ((TClonesArray*)fDigits->At(branch))->Clear();
  }
  if(fNDigits) fNDigits[branch]=0;

}



//_______________________________________________________________________
void AliITSDetTypeSim::SDigitsToDigits(Option_t* opt, Char_t* name){
  // Standard Summable digits to Digits function.
  if(!GetITSgeom()){
    Warning("SDigitsToDigits","GetITSgeom() is null!!");
    return;
  }
  
  const char *all = strstr(opt,"All");
  const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
			strstr(opt,"SSD")};
  if(!all && !det[0] && !det[1] && !det[2] ) all = "All";

  static Bool_t setDef = kTRUE;
  if(setDef) SetDefaultSimulation();
  setDef = kFALSE;
  
  AliITSsimulation *sim =0;
  TTree* trees = fLoader->TreeS();
  if( !(trees && GetSDigits()) ){
    Error("SDigits2Digits","Error: No trees or SDigits. Returning.");
    return;
  } 

  TBranch* brchSDigits = trees->GetBranch(name);
  
  Int_t id;
  for(Int_t module=0;module<GetITSgeom()->GetIndexMax();module++){
     id = GetITSgeom()->GetModuleType(module);
    if (!all && !det[id]) continue;
    sim = (AliITSsimulation*)GetSimulationModel(id);
    if(!sim){
      Error("SDigit2Digits","The simulation class was not "
	    "instanciated for module %d type %s!",module,
	    GetITSgeom()->GetModuleTypeName(module));
      exit(1);
    }
    sim->InitSimulationModule(module,gAlice->GetEvNumber());
    
    fSDigits.Clear();
    brchSDigits->GetEvent(module);
    sim->AddSDigitsToModule(&fSDigits,0);
    sim->FinishSDigitiseModule();
    fLoader->TreeD()->Fill();
    ResetDigits();
  }

  WriteFOSignals(); 
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  fLoader->TreeD()->Reset();
}
//_________________________________________________________
void AliITSDetTypeSim::AddSumDigit(AliITSpListItem &sdig){  
  //Adds the module full of summable digits to the summable digits tree.

  new(fSDigits[fNSDigits++]) AliITSpListItem(sdig);
}
//__________________________________________________________
void AliITSDetTypeSim::AddSimDigit(Int_t branch, const AliITSdigit* d){  
  //    Add a simulated digit.

  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  switch(branch){
  case 0:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSPD(*((AliITSdigitSPD*)d));
    break;
  case 1:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSDD(*((AliITSdigitSDD*)d));
    break;
  case 2:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSSD(*((AliITSdigitSSD*)d));
    break;
  }  
}
//______________________________________________________________________
void AliITSDetTypeSim::AddSimDigit(Int_t branch,Float_t phys,Int_t *digits,
				   Int_t *tracks,Int_t *hits,Float_t *charges, 
				   Int_t sigexpanded){
  //   Add a simulated digit to the list.

  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  switch(branch){
  case 0:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSPD(digits,tracks,hits);
    break;
  case 1:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSDD(phys,digits,tracks,
						   hits,charges,sigexpanded);
    break;
  case 2:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSSD(digits,tracks,hits);
    break;
  } 
}
//______________________________________________________________________
void AliITSDetTypeSim::ReadOldSSDNoise(const TObjArray *array, 
				       AliITSNoiseSSDv2 *noiseSSD) {
  //Reads the old SSD calibration object and converts it to the new format
  const Int_t fgkSSDSTRIPSPERMODULE = 1536;
  const Int_t fgkSSDPSIDESTRIPSPERMODULE = 768;

  Int_t gNMod = array->GetEntries();
  cout<<"Converting old calibration object for noise..."<<endl;

  //NOISE
  Double_t noise = 0.0;
  for (Int_t iModule = 0; iModule < gNMod; iModule++) {
    AliITSNoiseSSD *noiseModule = (AliITSNoiseSSD*) (array->At(iModule));
    for(Int_t iStrip = 0; iStrip < fgkSSDSTRIPSPERMODULE; iStrip++) {
      noise = (iStrip < fgkSSDPSIDESTRIPSPERMODULE) ? noiseModule->GetNoiseP(iStrip) : noiseModule->GetNoiseN(1535 - iStrip);
      if(iStrip < fgkSSDPSIDESTRIPSPERMODULE)
	noiseSSD->AddNoiseP(iModule,iStrip,noise);
      if(iStrip >= fgkSSDPSIDESTRIPSPERMODULE)
	noiseSSD->AddNoiseN(iModule,1535 - iStrip,noise);
    }//loop over strips
  }//loop over modules      
}

//______________________________________________________________________
void AliITSDetTypeSim::ReadOldSSDBadChannels(const TObjArray *array, 
					     AliITSBadChannelsSSDv2 *badChannelsSSD) {
  //Reads the old SSD calibration object and converts it to the new format
  Int_t nMod = array->GetEntries();
  cout<<"Converting old calibration object for bad channels..."<<endl;
  for (Int_t iModule = 0; iModule < nMod; iModule++) {
    //for (Int_t iModule = 0; iModule < 1; iModule++) {
    AliITSBadChannelsSSD *bad = (AliITSBadChannelsSSD*) (array->At(iModule));
    TArrayI arrayPSide = bad->GetBadPChannelsList();
    for(Int_t iPCounter = 0; iPCounter < arrayPSide.GetSize(); iPCounter++) 
      badChannelsSSD->AddBadChannelP(iModule,
				     iPCounter,
				     (Char_t)arrayPSide.At(iPCounter));
        
    TArrayI arrayNSide = bad->GetBadNChannelsList();
    for(Int_t iNCounter = 0; iNCounter < arrayNSide.GetSize(); iNCounter++) 
      badChannelsSSD->AddBadChannelN(iModule,
				     iNCounter,
				     (Char_t)arrayNSide.At(iNCounter));
    
  }//loop over modules      
}

//______________________________________________________________________
void AliITSDetTypeSim::ReadOldSSDGain(const TObjArray *array, 
				      AliITSGainSSDv2 *gainSSD) {
  //Reads the old SSD calibration object and converts it to the new format

  Int_t nMod = array->GetEntries();
  cout<<"Converting old calibration object for gain..."<<endl;

  //GAIN
  for (Int_t iModule = 0; iModule < nMod; iModule++) {
    AliITSGainSSD *gainModule = (AliITSGainSSD*) (array->At(iModule));
    TArrayF arrayPSide = gainModule->GetGainP();
    for(Int_t iPCounter = 0; iPCounter < arrayPSide.GetSize(); iPCounter++)
      gainSSD->AddGainP(iModule,
			iPCounter,
			arrayPSide.At(iPCounter));
    TArrayF arrayNSide = gainModule->GetGainN();
    for(Int_t iNCounter = 0; iNCounter < arrayNSide.GetSize(); iNCounter++)
      gainSSD->AddGainN(iModule,
			iNCounter,
			arrayNSide.At(iNCounter));
  }//loop over modules 
}
//______________________________________________________________________
void AliITSDetTypeSim::ProcessSPDDigitForFastOr(UInt_t module, UInt_t colM, UInt_t rowM) {
  // Processes wether a single fired pixel will give rise to a fast-or signal
  fFOGenerator.ProcessPixelHitM(module,colM,rowM);
}
//_______________________________________________________________________
AliITSTriggerConditions* AliITSDetTypeSim::GetTriggerConditions() {
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
//_______________________________________________________________________
void AliITSDetTypeSim::WriteFOSignals() {
  // write fo signals to event

  if (!fLoader) {
    AliError("ITS loader is NULL.");
    return;
  }

  if(!fLoader->TreeD()){
   AliError("No TreeD available");
   return;
  }
  TTree *tree = fLoader->TreeD();
  AliITSFOSignalsSPD *foSignals = new AliITSFOSignalsSPD(*GetFOSignals()); 
  tree->GetUserInfo()->Add(foSignals);
  fFOGenerator.ResetSignals();
}

