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
#include "AliITSgeom.h"
#include "AliITSpListItem.h"
#include "AliITSresponseSSD.h"
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
fGeom(),         //
fSimulation(),   // [NDet]
fSegmentation(), // [NDet]
fResponse(),     // [NMod]
fPreProcess(),   // [] e.g. Fill fHitModule with hits
fPostProcess(),  // [] e.g. Wright Raw data
fNSDigits(0),    //! number of SDigits
fSDigits(),      //! [NMod][NSDigits]
fNDigits(0),     //! number of Digits
fDigits(),       //! [NMod][NDigits]
fHitClassName(), // String with Hit class name.
fSDigClassName(),// String with SDigit class name.
fDigClassName(){ // String with digit class name.
    // Default Constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A properly zero-ed AliITSDetTypeSim class.
  fGeom = 0;
  fSimulation = new TObjArray(fgkNdettypes);
  fSegmentation = new TObjArray(fgkNdettypes);
  fResponse = 0;
  fPreProcess = 0;
  fPostProcess = 0;
  fNSDigits = 0;
  fSDigits = new TClonesArray("AliITSpListItem",1000);
  fDigits = new TObjArray(fgkNdettypes);
  fNDigits = new Int_t[fgkNdettypes];
  fLoader = 0;
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
  
    if(fGeom) delete fGeom;
    if(fSimulation){
      fSimulation->Delete();
      delete fSimulation;
      fSimulation = 0;
    }
    
    if(fSegmentation){
      fSegmentation->Delete();
      delete fSegmentation;
      fSegmentation = 0;
    }
    
    if(fResponse){

      fResponse->Delete();
      delete fResponse;
      fResponse = 0;
    }
    
    if(fPreProcess){
      fPreProcess->Delete();
      delete fPreProcess;
      fPreProcess = 0;
    }
    
    if(fPostProcess){
      fPostProcess->Delete();
      delete fPostProcess;
      fPostProcess = 0;
    }
    
    if (fLoader)
      {
	fLoader->GetModulesFolder()->Remove(this);
      }
    
    if (fSDigits) {
      fSDigits->Delete();
      delete fSDigits;
      fSDigits=0;
    }
    if (fDigits) {
      fDigits->Delete();
      delete fDigits;
      fDigits=0;
    }
  
}
//----------------------------------------------------------------------
AliITSDetTypeSim::AliITSDetTypeSim(const AliITSDetTypeSim &source) : TObject(source){
    // Copy Constructor for object AliITSDetTypeSim not allowed
  if(this==&source) return;
  Error("Copy constructor",
	"You are not allowed to make a copy of the AliITSDetTypeSim");
  exit(1);

        
}
//----------------------------------------------------------------------
AliITSDetTypeSim& AliITSDetTypeSim::operator=(const AliITSDetTypeSim &source){
    // The = operator for object AliITSDetTypeSim
 
    if(&source==this) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSDetTypeSIm");
    exit(1);
    return *this;
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
  if(fGeom==0) {
    Warning("GetSimulationModelByModule","fGeom is 0!");
    return 0;
  }
  
  return GetSimulationModel(fGeom->GetModuleType(module));
}
//______________________________________________________________________
void AliITSDetTypeSim::SetSegmentationModel(Int_t dettype,AliITSsegmentation *seg){
   
  //Set segmentation model for detector type
  if(fSegmentation==0x0) fSegmentation = new TObjArray(fgkNdettypes);
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
   if(fGeom==0){
     Warning("GetSegmentationModelByModule","fGeom is 0!");
     return 0;
   }     
   return GetSegmentationModel(fGeom->GetModuleType(module));

}
//_______________________________________________________________________
void AliITSDetTypeSim::CreateResponses() {

  //Create the container of response functions with correct size
  if (fResponse) {
    fResponse->Delete();
    delete fResponse;
  }

  Int_t nModTot = 0;
  for (Int_t i=0; i<fgkNdettypes; i++) nModTot += fNMod[i];
  fResponse = new TObjArray(nModTot);
  fResponse->SetOwner(kTRUE);
  fResponse->Clear();

}
//_______________________________________________________________________
void AliITSDetTypeSim::SetResponseModel(Int_t iMod, AliITSresponse *resp){

  //Set response model for modules

  if (fResponse==0) CreateResponses();
 
  if (fResponse->At(iMod)!=0)
    delete (AliITSresponse*) fResponse->At(iMod);
  fResponse->AddAt(resp, iMod);

}
/*
//_______________________________________________________________________
void AliITSDetTypeSim::SetResponse(Int_t dettype, Int_t iMod, AliITSresponse *resp){

  //Set response for the module iMod of type dettype
  if (fResponse==0) CreateResponses();

  Int_t nModBefore = 0;
  for (Int_t i=0; i<dettype; i++) nModBefore += fNMod[i];

  if (fResponse->At(nModBefore+iMod) != 0)
    delete (AliITSresponse*) fResponse->At(nModBefore+iMod);
  fResponse->AddAt(resp, nModBefore+iMod);

}
*/
//______________________________________________________________________
void AliITSDetTypeSim::ResetResponse(){

  //resets response array
  if(fResponse)
    fResponse->Clear();

}
//______________________________________________________________________
void AliITSDetTypeSim::ResetSegmentation(){
 
 //Resets segmentation array
  if(fSegmentation){
    for(Int_t i=0;i<fgkNdettypes;i++){
      if(fSegmentation->At(i))
	delete (AliITSsegmentation*)fSegmentation->At(i);
    }
  }
}

//_______________________________________________________________________
AliITSresponse* AliITSDetTypeSim::GetResponseModel(Int_t iMod){
   //Get response model for module number iMod 
 
  if(fResponse==0) {
    AliError("fResponse is 0!");
    return 0; 
  }

  return (AliITSresponse*)(fResponse->At(iMod));

}
//_______________________________________________________________________
void AliITSDetTypeSim::SetDefaults(){

  //Set defaults for segmentation and response


  if(fGeom==0){
    Warning("SetDefaults","fGeom is 0!");
    return;
  }

  if (fResponse==0) CreateResponses();

  AliITSsegmentation* seg;
  AliITSresponse* res;
  ResetResponse();
  ResetSegmentation();
   
  if(!GetCalibration()){AliFatal("Exit"); exit(0);}

  for(Int_t idet=0;idet<fgkNdettypes;idet++){
    //SPD
    if(idet==0){
      if(!GetSegmentationModel(idet)){
	seg = new AliITSsegmentationSPD(fGeom);
	SetSegmentationModel(idet,seg);
      }
      const char *kData0=(GetResponseModel(fGeom->GetStartSPD()))->DataType();
      if (strstr(kData0,"real")) {
	SetDigitClassName(idet,"AliITSdigit");
      }
      else {
	SetDigitClassName(idet,"AliITSdigitSPD");
      }
    }
    //SDD
    if(idet==1){
      if(!GetSegmentationModel(idet)){
	res = GetResponseModel(fGeom->GetStartSDD());   
	seg = new AliITSsegmentationSDD(fGeom,res);  
	SetSegmentationModel(idet,seg);
      }
      const char *kopt = GetResponseModel(fGeom->GetStartSDD())->ZeroSuppOption();
      if((!strstr(kopt,"2D"))&&(!strstr(kopt,"1D"))) {
	SetDigitClassName(idet,"AliITSdigit");
      }
      else {
	SetDigitClassName(idet,"AliITSdigitSDD");
      }
    
    }
    //SSD
    if(idet==2){
      if(!GetSegmentationModel(idet)){
	seg = new AliITSsegmentationSSD(fGeom);
	SetSegmentationModel(idet,seg);
      }
      const char *kData2 = (GetResponseModel(fGeom->GetStartSSD())->DataType());
      if (strstr(kData2,"real")) {
	SetDigitClassName(idet,"AliITSdigit");
      }
      else {
	SetDigitClassName(idet,"AliITSdigitSSD");
      }
    }
  }
 
}

//______________________________________________________________________
Bool_t AliITSDetTypeSim::GetCalibration() {
  // Get Default calibration if a storage is not defined.

  Bool_t deleteManager = kFALSE;
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    deleteManager = kTRUE;
  }
  AliCDBStorage *storage = AliCDBManager::Instance()->GetDefaultStorage();

  AliCDBEntry *entrySPD = storage->Get("ITS/Calib/RespSPD", fRunNumber);
  TObjArray *respSPD = (TObjArray *)entrySPD->GetObject();
  entrySPD->SetObject(NULL);
  entrySPD->SetOwner(kTRUE);
  AliCDBEntry *entrySDD = storage->Get("ITS/Calib/RespSDD", fRunNumber);
  TObjArray *respSDD = (TObjArray *)entrySDD->GetObject();
  entrySDD->SetObject(NULL);
  entrySDD->SetOwner(kTRUE);
  AliCDBEntry *entrySSD = storage->Get("ITS/Calib/RespSSD", fRunNumber);
  TObjArray *respSSD = (TObjArray *)entrySSD->GetObject();
  entrySSD->SetObject(NULL);
  entrySSD->SetOwner(kTRUE);
  // DB entries are dleted. In this waymetadeta objects are deleted as well
  delete entrySPD;
  delete entrySDD;
  delete entrySSD;
  if(deleteManager){
    AliCDBManager::Instance()->Destroy();
    AliCDBManager::Instance()->RemoveDefaultStorage();
    storage = 0;   // the storage is killed by AliCDBManager::Instance()->Destroy()
  }

  if ((! respSPD)||(! respSDD)||(! respSSD)) {
    AliWarning("Can not get calibration from calibration database !");
    return kFALSE;
  }
  fNMod[0] = respSPD->GetEntries();
  fNMod[1] = respSDD->GetEntries();
  fNMod[2] = respSSD->GetEntries();
  AliInfo(Form("%i SPD, %i SDD and %i SSD in calibration database",
	       fNMod[0], fNMod[1], fNMod[2]));
  AliITSresponse* res;
  for (Int_t i=0; i<fNMod[0]; i++) {
    res = (AliITSresponse*) respSPD->At(i);
    SetResponseModel(i, res);
 }
  for (Int_t i=0; i<fNMod[1]; i++) {
    res = (AliITSresponse*) respSDD->At(i);
    Int_t iMod = i + fNMod[0];
    SetResponseModel(iMod, res);
 }
  for (Int_t i=0; i<fNMod[2]; i++) {
    res = (AliITSresponse*) respSSD->At(i);
    Int_t iMod = i + fNMod[0] + fNMod[1];
    SetResponseModel(iMod, res);
 }

  return kTRUE;
}



//_______________________________________________________________________
void AliITSDetTypeSim::SetDefaultSimulation(){

  //Set default simulation for detector type

  if(fGeom==0){
    Warning("SetDefaultSimulation","fGeom is 0!");
    return;
  }
  if(fResponse==0){
    Warning("SetDefaultSimulation","fResponse is 0!");
    return;
  }

  AliITSsimulation* sim;

  for(Int_t idet=0;idet<fgkNdettypes;idet++){
   //SPD
    if(idet==0){
      sim = GetSimulationModel(idet);
 
      if(!sim){
	sim = new AliITSsimulationSPD(this);
	SetSimulationModel(idet,sim);
      } else{
	// loop over all SPD modules
	Int_t nSPD = fGeom->GetStartSDD()-fGeom->GetStartSPD();
	for(Int_t nsp=fGeom->GetStartSPD();nsp<nSPD;nsp++){
	  sim->SetResponseModel(nsp,GetResponseModel(nsp));
	}
	sim->SetSegmentationModel(0,(AliITSsegmentationSPD*)GetSegmentationModel(idet));
	sim->Init();
      }
    }
    //SDD
    if(idet==1){
      sim = GetSimulationModel(idet);
      if(!sim){
	sim = new AliITSsimulationSDD(this);
	SetSimulationModel(idet,sim);
      } else {
	Int_t nSDD = fGeom->GetStartSSD()-fGeom->GetStartSDD();
	for(Int_t nsd=fGeom->GetStartSDD();nsd<nSDD;nsd++){
	  sim->SetResponseModel(nsd,GetResponseModel(nsd));
	}

	sim->SetSegmentationModel(1,(AliITSsegmentationSDD*)GetSegmentationModel(idet));
	sim->Init();
      }
      
    }
    //SSD
    if(idet==2){
      sim = GetSimulationModel(idet);
      if(!sim){
	sim = new AliITSsimulationSSD(this);
	SetSimulationModel(idet,sim);

      } else{

	Int_t nSSD = fGeom->GetLastSSD()-fGeom->GetStartSSD()+1;
	for(Int_t nss=fGeom->GetStartSSD();nss<nSSD;nss++){
	  sim->SetResponseModel(nss,GetResponseModel(nss));
	}

	sim->SetResponseModel(fGeom->GetStartSSD(),(AliITSresponseSSD*)GetResponseModel(fGeom->GetStartSSD()));
	sim->SetSegmentationModel(2,(AliITSsegmentationSSD*)GetSegmentationModel(idet));

    
	sim->Init();
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
  if (fSDigits == 0x0){
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
  if(!fGeom){
    Warning("SDigitsToDigits","fGeom is null!!");
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
  for(Int_t module=0;module<fGeom->GetIndexMax();module++){
     id = fGeom->GetModuleType(module);
    if (!all && !det[id]) continue;
    sim = (AliITSsimulation*)GetSimulationModel(id);
    printf("module=%d name=%s\n",module,sim->ClassName());
    if(!sim){
      Error("SDigit2Digits","The simulation class was not "
	    "instanciated for module %d type %s!",module,
	    fGeom->GetModuleTypeName(module));
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
  AliITSresponseSDD *resp = 0;
  switch(branch){
  case 0:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSPD(digits,tracks,hits);
    break;
  case 1:
    resp = (AliITSresponseSDD*)GetResponseModel(fGeom->GetStartSDD());
    new(ldigits[fNDigits[branch]++]) AliITSdigitSDD(phys,digits,tracks,
						   hits,charges,resp);
    break;
  case 2:
    new(ldigits[fNDigits[branch]++]) AliITSdigitSSD(digits,tracks,hits);
    break;
  } 
}



//______________________________________________________________________
void AliITSDetTypeSim::StoreCalibration(Int_t firstRun, Int_t lastRun, AliCDBMetaData &md) {

  // Store calibration in the calibration database

  // The database must be created in an external piece of code (i.e. 
  // a configuration macro )

  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    //AliError("No storage set!");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    //return;
  }

  if (!fResponse) {
    AliError("AliITSresponse classes are not defined - nothing done");
    return;
  }
  AliCDBId idRespSPD("ITS/Calib/RespSPD",firstRun, lastRun);
  AliCDBId idRespSDD("ITS/Calib/RespSDD",firstRun, lastRun);
  AliCDBId idRespSSD("ITS/Calib/RespSSD",firstRun, lastRun);

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
    respSPD.Add(fResponse->At(i));

  for (Int_t i = index[0]; i<index[1]; i++ )
    respSDD.Add(fResponse->At(i));

  for (Int_t i = index[1]; i<index[2]; i++ )
    respSSD.Add(fResponse->At(i));

  AliCDBManager::Instance()->GetDefaultStorage()->Put(&respSPD, idRespSPD, &md);

  AliCDBManager::Instance()->GetDefaultStorage()->Put(&respSDD, idRespSDD, &md);
  
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&respSSD, idRespSSD, &md);
}


