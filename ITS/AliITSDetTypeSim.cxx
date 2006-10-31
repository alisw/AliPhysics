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
#include "AliITSDetTypeSim.h"
#include "AliITSpListItem.h"
#include "AliITSresponseSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulation.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"


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
fPreProcess(),   // [] e.g. Fill fHitModule with hits
fPostProcess(),  // [] e.g. Wright Raw data
fNSDigits(0),    //! number of SDigits
fSDigits(),      //! [NMod][NSDigits]
fNDigits(0),     //! number of Digits
fRunNumber(0),   //! Run number (to access DB)
fDigits(),       //! [NMod][NDigits]
fHitClassName(), // String with Hit class name.
fSDigClassName(),// String with SDigit class name.
fDigClassName(), // String with digit class name.
fLoader(0),      // local pointer to loader
fFirstcall(kTRUE){ // flag
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
  fSDigits = new TClonesArray("AliITSpListItem",1000);
  fDigits = new TObjArray(fgkNdettypes);
  fNDigits = new Int_t[fgkNdettypes];
  fNMod[0] = fgkDefaultNModulesSPD;
  fNMod[1] = fgkDefaultNModulesSDD;
  fNMod[2] = fgkDefaultNModulesSSD;
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
	AliITSresponse* rspd = ((AliITSCalibration*)fCalibration->At(
                               GetITSgeom()->GetStartSPD()))->GetResponse();
	AliITSresponse* rsdd = ((AliITSCalibration*)fCalibration->At(
                               GetITSgeom()->GetStartSDD()))->GetResponse();
	AliITSresponse* rssd = ((AliITSCalibration*)fCalibration->At(
                               GetITSgeom()->GetStartSSD()))->GetResponse();
	if(rspd) delete rspd;
	if(rsdd) delete rsdd;
	if(rssd) delete rssd;
	fCalibration->Delete();
	delete fCalibration;
    }
    fCalibration = 0;
    if(fPreProcess){
	fPreProcess->Delete();
	delete fPreProcess;
    }
    fPreProcess = 0;
    if(fPostProcess){
	fPostProcess->Delete();
	delete fPostProcess;
    }
    fPostProcess = 0;
    if(fNDigits) delete [] fNDigits;
    fNDigits = 0;
    if (fLoader)fLoader->GetModulesFolder()->Remove(this);
    fLoader = 0; // Not deleting it.
    if (fSDigits) {
	fSDigits->Delete();
	delete fSDigits;
    }
    fSDigits=0;
    if (fDigits) {
      fDigits->Delete();
      delete fDigits;
    }
    fDigits=0;
}
//----------------------------------------------------------------------
AliITSDetTypeSim::AliITSDetTypeSim(const AliITSDetTypeSim &source) : TObject(source),
fSimulation(source.fSimulation),   // [NDet]
fSegmentation(source.fSegmentation), // [NDet]
fCalibration(source.fCalibration),     // [NMod]
fPreProcess(source.fPreProcess),   // [] e.g. Fill fHitModule with hits
fPostProcess(source.fPostProcess),  // [] e.g. Wright Raw data
fNSDigits(source.fNSDigits),    //! number of SDigits
fSDigits(source.fSDigits),      //! [NMod][NSDigits]
fNDigits(source.fNDigits),     //! number of Digits
fRunNumber(source.fRunNumber),   //! Run number (to access DB)
fDigits(source.fDigits),       //! [NMod][NDigits]
fHitClassName(source.fHitClassName), // String with Hit class name.
fSDigClassName(source.fSDigClassName),// String with SDigit class name.
fDigClassName(), // String with digit class name.
fLoader(source.fLoader),      // local pointer to loader
fFirstcall(source.fFirstcall)
{
    // Copy Constructor for object AliITSDetTypeSim not allowed
  for(Int_t i=0;i<fgkNdettypes;i++){
    fDigClassName[i] = source.fDigClassName[i];
  }
}
//----------------------------------------------------------------------
AliITSDetTypeSim& AliITSDetTypeSim::operator=(const AliITSDetTypeSim &source){
    // The = operator for object AliITSDetTypeSim
 
  this->~AliITSDetTypeSim();
  new(this) AliITSDetTypeSim(source);
  return *this;
}

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
AliITSsimulation* AliITSDetTypeSim::GetSimulationModel(Int_t dettype){

  //Get simulation model for detector type
  if(fSimulation==0)  {
    Warning("GetSimulationModel","fSimulation is 0!");
    return 0;     
  }
  return (AliITSsimulation*)(fSimulation->At(dettype));
}
//______________________________________________________________________
AliITSsimulation* AliITSDetTypeSim::GetSimulationModelByModule(Int_t module){

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
	seg = new AliITSsegmentationSPD(GetITSgeom());
    }else if(idet==1){
	AliITSCalibration* res=GetCalibrationModel(
	    GetITSgeom()->GetStartSDD());
	seg = new AliITSsegmentationSDD(GetITSgeom(),res);
    }else {
	seg = new AliITSsegmentationSSD(GetITSgeom());
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
AliITSsegmentation* AliITSDetTypeSim::GetSegmentationModel(Int_t dettype){
  //Get segmentation model for detector type
   
   if(fSegmentation==0) {
       Warning("GetSegmentationModel","fSegmentation is 0!");
       return 0; 
   } 
   return (AliITSsegmentation*)(fSegmentation->At(dettype));
}
//_______________________________________________________________________
AliITSsegmentation* AliITSDetTypeSim::GetSegmentationModelByModule(Int_t module){
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
//______________________________________________________________________
void AliITSDetTypeSim::ResetCalibrationArray(){
    //resets response array
    if(fCalibration && fRunNumber<0){  // if fRunNumber<0 fCalibration is owner
	AliITSresponse* rspd = ((AliITSCalibration*)fCalibration->At(
                                GetITSgeom()->GetStartSPD()))->GetResponse();
	AliITSresponse* rsdd = ((AliITSCalibration*)fCalibration->At(
                                GetITSgeom()->GetStartSDD()))->GetResponse();
	AliITSresponse* rssd = ((AliITSCalibration*)fCalibration->At(
                                GetITSgeom()->GetStartSSD()))->GetResponse();
	if(rspd) delete rspd;
	if(rsdd) delete rsdd;
	if(rssd) delete rssd;
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
AliITSCalibration* AliITSDetTypeSim::GetCalibrationModel(Int_t iMod){
    //Get response model for module number iMod 
 
    if(fCalibration==0) {
	AliError("fCalibration is 0!");
	return 0; 
    }
  return (AliITSCalibration*)(fCalibration->At(iMod));
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

    for(Int_t idet=0;idet<fgkNdettypes;idet++){
	//SPD
	if(idet==0){
	    if(!GetSegmentationModel(idet)) SetDefaultSegmentation(idet);
	    const char *kData0=(GetCalibrationModel(GetITSgeom()->GetStartSPD()))->DataType();
	    if (strstr(kData0,"real")) {
		SetDigitClassName(idet,"AliITSdigit");
	    }else {
		SetDigitClassName(idet,"AliITSdigitSPD");
	    } // end if
	} // end if idet==0
	//SDD
	if(idet==1){
	    if(!GetSegmentationModel(idet)) SetDefaultSegmentation(idet);
	    AliITSCalibrationSDD* rsp = 
		(AliITSCalibrationSDD*)GetCalibrationModel(
		    GetITSgeom()->GetStartSDD());
	    const char *kopt = ((AliITSresponseSDD*)rsp->GetResponse())->
		ZeroSuppOption();
	    if((!strstr(kopt,"2D"))&&(!strstr(kopt,"1D"))) {
		SetDigitClassName(idet,"AliITSdigit");
	    }else {
		SetDigitClassName(idet,"AliITSdigitSDD");
	    } // end if
	} // end if idet==1
	//SSD
	if(idet==2){
	    if(!GetSegmentationModel(idet))SetDefaultSegmentation(idet);
	    const char *kData2 = (GetCalibrationModel(
                                  GetITSgeom()->GetStartSSD())->DataType());
	    if (strstr(kData2,"real")) {
		SetDigitClassName(idet,"AliITSdigit");
	    }else {
		SetDigitClassName(idet,"AliITSdigitSSD");
	    } // end if
	} // end if idet==2
    }// end for idet
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
  if(run<0)run=0;   // if the run number is not yet set, use fake run # 0

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

  AliCDBEntry *entrySPD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSPD", run);
  AliCDBEntry *entrySDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD", run);
  AliCDBEntry *entrySSD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSSD", run);
  AliCDBEntry *entry2SPD = AliCDBManager::Instance()->Get("ITS/Calib/RespSPD", run);
  AliCDBEntry *entry2SDD = AliCDBManager::Instance()->Get("ITS/Calib/RespSDD", run);
  AliCDBEntry *entry2SSD = AliCDBManager::Instance()->Get("ITS/Calib/RespSSD", run);

  if(!entrySPD || !entrySDD || !entrySSD || !entry2SPD || !entry2SDD || !entry2SSD){
  	AliWarning("Calibration object retrieval failed! Dummy calibration will be used.");
	AliCDBStorage *localStor = 
		AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT");
	
  	entrySPD = localStor->Get("ITS/Calib/CalibSPD", run);
  	entrySDD = localStor->Get("ITS/Calib/CalibSDD", run);
  	entrySSD = localStor->Get("ITS/Calib/CalibSSD", run);
 	entry2SPD = localStor->Get("ITS/Calib/RespSPD", run);
  	entry2SDD = localStor->Get("ITS/Calib/RespSDD", run);
  	entry2SSD = localStor->Get("ITS/Calib/RespSSD", run);
  }

  if(!entrySPD || !entrySDD || !entrySSD || !entry2SPD || !entry2SDD || !entry2SSD){
    AliError("Calibration data was not found in $ALICE_ROOT!");
    return kFALSE;
  }

  TObjArray *calSPD = (TObjArray *)entrySPD->GetObject();
  if(!isCacheActive)entrySPD->SetObject(NULL);
  entrySPD->SetOwner(kTRUE);

  AliITSresponseSPD *pSPD = (AliITSresponseSPD*)entry2SPD->GetObject();
  if(!isCacheActive)entry2SPD->SetObject(NULL);
  entry2SPD->SetOwner(kTRUE);
   
  TObjArray *calSDD = (TObjArray *)entrySDD->GetObject();
  if(!isCacheActive)entrySDD->SetObject(NULL);
  entrySDD->SetOwner(kTRUE);

  AliITSresponseSDD *pSDD = (AliITSresponseSDD*)entry2SDD->GetObject();
  if(!isCacheActive)entry2SDD->SetObject(NULL);
  entry2SDD->SetOwner(kTRUE);

  TObjArray *calSSD = (TObjArray *)entrySSD->GetObject();
  if(!isCacheActive)entrySSD->SetObject(NULL);
  entrySSD->SetOwner(kTRUE);

  AliITSresponseSSD *pSSD = (AliITSresponseSSD*)entry2SSD->GetObject();
  if(!isCacheActive)entry2SSD->SetObject(NULL);
  entry2SSD->SetOwner(kTRUE);
  
  // DB entries are deleted. In this way metadeta objects are deleted as well
  if(!isCacheActive){
    delete entrySPD;
    delete entrySDD;
    delete entrySSD;
    delete entry2SPD;
    delete entry2SDD;
    delete entry2SSD;
  }
  
  AliCDBManager::Instance()->SetCacheFlag(origCacheStatus);

  if ((!pSPD)||(!pSDD)||(!pSSD)) {
    AliWarning("Can not get calibration from calibration database !");
    return kFALSE;
  }
  if ((! calSPD)||(! calSDD)||(! calSSD)) {
    AliWarning("Can not get calibration from calibration database !");
    return kFALSE;
  }
  fNMod[0] = calSPD->GetEntries();
  fNMod[1] = calSDD->GetEntries();
  fNMod[2] = calSSD->GetEntries();
  AliInfo(Form("%i SPD, %i SDD and %i SSD in calibration database",
	       fNMod[0], fNMod[1], fNMod[2]));
  AliITSCalibration* cal;
  for (Int_t i=0; i<fNMod[0]; i++) {
    cal = (AliITSCalibration*) calSPD->At(i);
    cal->SetResponse(pSPD);
    SetCalibrationModel(i, cal);
 }
  for (Int_t i=0; i<fNMod[1]; i++) {
    cal = (AliITSCalibration*) calSDD->At(i);
    cal->SetResponse(pSDD);
    Int_t iMod = i + fNMod[0];
    SetCalibrationModel(iMod, cal);
 }
  for (Int_t i=0; i<fNMod[2]; i++) {
    cal = (AliITSCalibration*) calSSD->At(i);
    cal->SetResponse(pSSD);
    Int_t iMod = i + fNMod[0] + fNMod[1];
    SetCalibrationModel(iMod, cal);
 }
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
void AliITSDetTypeSim::SetTreeAddressS(TTree* treeS, Char_t* name){
  // Set branch address for the ITS summable digits Trees.  
  char branchname[30];

  if(!treeS){
    return;
  }
  if (fSDigits ==  0x0){
    fSDigits = new TClonesArray("AliITSpListItem",1000);
  }
  TBranch *branch;
  sprintf(branchname,"%s",name);
  branch = treeS->GetBranch(branchname);
  if (branch) branch->SetAddress(&fSDigits);

}
//___________________________________________________________________
void AliITSDetTypeSim::SetTreeAddressD(TTree* treeD, Char_t* name){
  // Set branch address for the digit Trees.
  
  const char *det[3] = {"SPD","SDD","SSD"};
  TBranch *branch;
  
  char branchname[30];
  
  if(!treeD){
    return;
  }
  if(!fDigits){
    fDigits = new TObjArray(fgkNdettypes); 
  }
  for(Int_t i=0;i<fgkNdettypes;i++){
    Char_t* digclass = GetDigitClassName(i);
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
    
    if(fgkNdettypes==3) sprintf(branchname,"%sDigits%s",name,det[i]);
    else sprintf(branchname,"%sDigits%d",name,i+1);
    if(fDigits){
      branch = treeD->GetBranch(branchname);
      if(branch) branch->SetAddress(&((*fDigits)[i]));
    }
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
  if( !det[0] && !det[1] && !det[2] ) all = "All";
  else all = 0;
  static Bool_t setDef = kTRUE;
  if(setDef) SetDefaultSimulation();
  setDef = kFALSE;
  
  AliITSsimulation *sim =0;
  TTree* trees = fLoader->TreeS();
  if( !(trees && GetSDigits()) ){
    Error("SDigits2Digits","Error: No trees or SDigits. Returning.");
    return;
  } 
  sprintf(name,"%s",name);
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
    
    fSDigits->Clear();
    brchSDigits->GetEvent(module);
    sim->AddSDigitsToModule(fSDigits,0);
    sim->FinishSDigitiseModule();
    fLoader->TreeD()->Fill();
    ResetDigits();
  }
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  fLoader->TreeD()->Reset();
}
//_________________________________________________________
void AliITSDetTypeSim::AddSumDigit(AliITSpListItem &sdig){  
  //Adds the module full of summable digits to the summable digits tree.

  TClonesArray &lsdig = *fSDigits;
  new(lsdig[fNSDigits++]) AliITSpListItem(sdig);
}
//__________________________________________________________
void AliITSDetTypeSim::AddRealDigit(Int_t branch, Int_t *digits){
  //   Add a real digit - as coming from data.

  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  new(ldigits[fNDigits[branch]++]) AliITSdigit(digits); 
}
//__________________________________________________________
void AliITSDetTypeSim::AddSimDigit(Int_t branch, AliITSdigit* d){  
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
				   Int_t *tracks,Int_t *hits,Float_t *charges){
  //   Add a simulated digit to the list.

  TClonesArray &ldigits = *((TClonesArray*)fDigits->At(branch));
  AliITSCalibrationSDD *resp = 0;
  switch(branch){
  case 0:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSPD(digits,tracks,hits);
    break;
  case 1:
    resp = (AliITSCalibrationSDD*)GetCalibrationModel(
	GetITSgeom()->GetStartSDD());
    new(ldigits[fNDigits[branch]++]) AliITSdigitSDD(phys,digits,tracks,
						   hits,charges,resp);
    break;
  case 2:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSSD(digits,tracks,hits);
    break;
  } 
}
//______________________________________________________________________
void AliITSDetTypeSim::StoreCalibration(Int_t firstRun, Int_t lastRun,
					AliCDBMetaData &md) {
  // Store calibration in the calibration database
  // The database must be created in an external piece of code (i.e. 
  // a configuration macro )

  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliWarning("No storage set! Will use dummy one");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  }

  if (!fCalibration) {
    AliError("AliITSCalibration classes are not defined - nothing done");
    return;
  }
  AliCDBId idRespSPD("ITS/Calib/CalibSPD",firstRun, lastRun);
  AliCDBId idRespSDD("ITS/Calib/CalibSDD",firstRun, lastRun);
  AliCDBId idRespSSD("ITS/Calib/CalibSSD",firstRun, lastRun);

  TObjArray respSPD(fNMod[0]);
  TObjArray respSDD(fNMod[1]-fNMod[0]);
  TObjArray respSSD(fNMod[2]-fNMod[1]);
  respSPD.SetOwner(kFALSE);
  respSSD.SetOwner(kFALSE);
  respSSD.SetOwner(kFALSE);

  Int_t index[fgkNdettypes];
  for (Int_t i = 0; i<fgkNdettypes; i++ ) {
    index[i] = 0;
    for (Int_t j = 0; j<=i; j++ )
      index[i]+=fNMod[j];
  }

  for (Int_t i = 0; i<index[0]; i++ )
    respSPD.Add(fCalibration->At(i));

  for (Int_t i = index[0]; i<index[1]; i++ )
    respSDD.Add(fCalibration->At(i));

  for (Int_t i = index[1]; i<index[2]; i++ )
    respSSD.Add(fCalibration->At(i));

  AliCDBManager::Instance()->Put(&respSPD, idRespSPD, &md);
  AliCDBManager::Instance()->Put(&respSDD, idRespSDD, &md);
  AliCDBManager::Instance()->Put(&respSSD, idRespSSD, &md);
}


