/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Conributors are mentioned in the code where appropriate.              *
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

////////////////////////////////////////////////////////////////////////
// This class defines the "Standard" reconstruction for the ITS       // 
// detector.                                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////
#include "TObjArray.h"
#include "TTree.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliITSClusterFinder.h"
#include "AliITSClusterFinderV2.h"
#include "AliITSClusterFinderV2SPD.h"
#include "AliITSClusterFinderV2SDD.h"
#include "AliITSClusterFinderV2SSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeom.h"
#include "AliITSRawCluster.h"
#include "AliITSRawClusterSPD.h"
#include "AliITSRawClusterSDD.h"
#include "AliITSRawClusterSSD.h"
#include "AliITSRecPoint.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliLog.h"


const Int_t AliITSDetTypeRec::fgkNdettypes = 3;
const Int_t AliITSDetTypeRec::fgkDefaultNModulesSPD =  240;
const Int_t AliITSDetTypeRec::fgkDefaultNModulesSDD =  260;
const Int_t AliITSDetTypeRec::fgkDefaultNModulesSSD = 1698;

ClassImp(AliITSDetTypeRec)

//________________________________________________________________
AliITSDetTypeRec::AliITSDetTypeRec(): TObject(),
fNMod(0),
fITSgeom(0),
fReconstruction(0),
fSegmentation(0),
fCalibration(0),
fPreProcess(0),
fPostProcess(0),
fDigits(0),
fNdtype(0),
fCtype(0),
fNctype(0),
fRecPoints(0),
fNRecPoints(0),
fSelectedVertexer(),
fFirstcall(kTRUE){
    // Standard Constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //   

  fReconstruction = new TObjArray(fgkNdettypes);
  fDigits = new TObjArray(fgkNdettypes);
  for(Int_t i=0; i<3; i++){
    fClusterClassName[i]=0;
    fDigClassName[i]=0;
    fRecPointClassName[i]=0;
  }
  fNdtype = new Int_t[fgkNdettypes];
  fCtype = new TObjArray(fgkNdettypes);
  fNctype = new Int_t[fgkNdettypes];
  fNMod = new Int_t [fgkNdettypes];
  fNMod[0] = fgkDefaultNModulesSPD;
  fNMod[1] = fgkDefaultNModulesSDD;
  fNMod[2] = fgkDefaultNModulesSSD;
  fRecPoints = new TClonesArray("AliITSRecPoint",3000);
  fNRecPoints = 0;
  
  for(Int_t i=0;i<fgkNdettypes;i++){
    fNdtype[i]=0;
    fNctype[i]=0;
  }
  
  SelectVertexer(" "); 
}

//______________________________________________________________________
AliITSDetTypeRec::AliITSDetTypeRec(const AliITSDetTypeRec & rec):TObject(rec),
fNMod(rec.fNMod),
fITSgeom(rec.fITSgeom),
fReconstruction(rec.fReconstruction),
fSegmentation(rec.fSegmentation),
fCalibration(rec.fCalibration),
fPreProcess(rec.fPreProcess),
fPostProcess(rec.fPostProcess),
fDigits(rec.fDigits),
fNdtype(rec.fNdtype),
fCtype(rec.fCtype),
fNctype(rec.fNctype),
fRecPoints(rec.fRecPoints),
fNRecPoints(rec.fNRecPoints),
fSelectedVertexer(rec.fSelectedVertexer),
fFirstcall(rec.fFirstcall)
{

  // Copy constructor. 

}
//______________________________________________________________________
AliITSDetTypeRec& AliITSDetTypeRec::operator=(const AliITSDetTypeRec& source){
    // Assignment operator. 
    this->~AliITSDetTypeRec();
    new(this) AliITSDetTypeRec(source);
    return *this;

}

//_____________________________________________________________________
AliITSDetTypeRec::~AliITSDetTypeRec(){
  //Destructor

  if(fReconstruction){
    fReconstruction->Delete();
    delete fReconstruction;
    fReconstruction = 0;
  }
  if(fSegmentation){
    fSegmentation->Delete();
    delete fSegmentation;
    fSegmentation = 0;
  }
  if(fCalibration){
    AliITSresponse* rspd = ((AliITSCalibration*)fCalibration->At(GetITSgeom()->GetStartSPD()))->GetResponse();    
    AliITSresponse* rsdd = ((AliITSCalibration*)fCalibration->At(GetITSgeom()->GetStartSDD()))->GetResponse();
    AliITSresponse* rssd = ((AliITSCalibration*)fCalibration->At(GetITSgeom()->GetStartSSD()))->GetResponse();
    if(rspd) delete rspd;
    if(rsdd) delete rsdd;
    if(rssd) delete rssd;
    fCalibration->Delete();
    delete fCalibration;
    fCalibration = 0;
  }
  if(fPreProcess) delete fPreProcess;
  if(fPostProcess) delete fPostProcess;

  if(fDigits){
    fDigits->Delete();
    delete fDigits;
    fDigits=0;
  }
  if(fRecPoints){
    fRecPoints->Delete();
    delete fRecPoints;
    fRecPoints=0;
  }
  if(fCtype) {
    fCtype->Delete();
    delete fCtype;
    fCtype = 0;
  }
  delete [] fNctype;
  delete [] fNdtype;
  delete [] fNMod;
  
  if (fITSgeom) delete fITSgeom;
 
}

//___________________________________________________________________
void AliITSDetTypeRec::SetReconstructionModel(Int_t dettype,AliITSClusterFinder *clf){

  //Set reconstruction model for detector type

  if(fReconstruction==0) fReconstruction = new TObjArray(fgkNdettypes);
  if(fReconstruction->At(dettype)!=0) delete fReconstruction->At(dettype);
  fReconstruction->AddAt(clf,dettype);
}
//______________________________________________________________________
AliITSClusterFinder* AliITSDetTypeRec::GetReconstructionModel(Int_t dettype){

  //Get reconstruction model for detector type
  if(fReconstruction==0)  {
    Warning("GetReconstructionModel","fReconstruction is 0!");
    return 0;     
  }
  return (AliITSClusterFinder*)fReconstruction->At(dettype);
}

//______________________________________________________________________
void AliITSDetTypeRec::SetSegmentationModel(Int_t dettype,AliITSsegmentation *seg){
   
  //Set segmentation model for detector type
  
  if(fSegmentation==0) fSegmentation = new TObjArray(fgkNdettypes);
  if(fSegmentation->At(dettype)!=0) delete fSegmentation->At(dettype);
  fSegmentation->AddAt(seg,dettype);

}
//______________________________________________________________________
AliITSsegmentation* AliITSDetTypeRec::GetSegmentationModel(Int_t dettype){

  //Get segmentation model for detector type
   
   if(fSegmentation==0) {
     Warning("GetSegmentationModel","fSegmentation is 0!");
     return 0; 
   } 
   return (AliITSsegmentation*)fSegmentation->At(dettype);

}
//_______________________________________________________________________
void AliITSDetTypeRec::SetCalibrationModel(Int_t iMod, AliITSCalibration *cal){

  //Set calibration (response) for the module iMod of type dettype
  if (fCalibration==0) {
    fCalibration = new TObjArray(GetITSgeom()->GetIndexMax());
    fCalibration->SetOwner(kTRUE);
    fCalibration->Clear();
  }

  if (fCalibration->At(iMod) != 0)
    delete (AliITSCalibration*) fCalibration->At(iMod);
  fCalibration->AddAt(cal,iMod);

}
//_______________________________________________________________________
AliITSCalibration* AliITSDetTypeRec::GetCalibrationModel(Int_t iMod){
  
  //Get calibration model for module type
  
  if(fCalibration==0) {
    Warning("GetalibrationModel","fCalibration is 0!");
    return 0; 
  }  

  return (AliITSCalibration*)fCalibration->At(iMod);
}

//______________________________________________________________________
void AliITSDetTypeRec::SetTreeAddressD(TTree *treeD){
    // Set branch address for the tree of digits.

    const char *det[4] = {"SPD","SDD","SSD","ITS"};
    TBranch *branch;
    Char_t* digclass;
    Int_t i;
    char branchname[30];

    if(!treeD) return;
    if (fDigits == 0x0) fDigits = new TObjArray(fgkNdettypes);
    for (i=0; i<fgkNdettypes; i++) {
        digclass = GetDigitClassName(i);
	if(!(fDigits->At(i))) {
            fDigits->AddAt(new TClonesArray(digclass,1000),i);
        }else{
            ResetDigits(i);
        } 
        if (fgkNdettypes==3) sprintf(branchname,"%sDigits%s",det[3],det[i]);
        else  sprintf(branchname,"%sDigits%d",det[3],i+1);
        if (fDigits) {
            branch = treeD->GetBranch(branchname);
            if (branch) branch->SetAddress(&((*fDigits)[i]));
        } 
    } 
}

//_______________________________________________________________________
TBranch* AliITSDetTypeRec::MakeBranchInTree(TTree *tree, const char* name, 
                                       const char *classname, 
                                       void* address,Int_t size, 
                                       Int_t splitlevel)
{ 
//
// Makes branch in given tree and diverts them to a separate file
// 
//
//
    
  if (tree == 0x0) {
    Error("MakeBranchInTree","Making Branch %s Tree is NULL",name);
    return 0x0;
  }
  TBranch *branch = tree->GetBranch(name);
  if (branch) {  
    return branch;
  }
  if (classname){
    branch = tree->Branch(name,classname,address,size,splitlevel);
  }
  else {
    branch = tree->Bronch(name, "TClonesArray", address, size, splitlevel);
  }
  
  return branch;
}

//____________________________________________________________________
void AliITSDetTypeRec::SetDefaults(){
  
  //Set defaults for segmentation and response

  if(!GetITSgeom()){
    Warning("SetDefaults","null pointer to AliITSgeomGeom !");
    return;
  }

  AliITSsegmentation* seg;
  if(!GetCalibration()) {AliFatal("Exit");exit(0);}  

  for(Int_t dettype=0;dettype<fgkNdettypes;dettype++){
    if(dettype==0){
      seg = new AliITSsegmentationSPD(GetITSgeom());
      SetSegmentationModel(dettype,seg);
      SetDigitClassName(dettype,"AliITSdigitSPD");
      SetClusterClassName(dettype,"AliITSRawClusterSPD");

    }
    if(dettype==1){
      AliITSCalibrationSDD* res=(AliITSCalibrationSDD*) GetCalibrationModel(GetITSgeom()->GetStartSDD()); 
      seg = new AliITSsegmentationSDD(GetITSgeom(),res);
      SetSegmentationModel(dettype,seg);
      const char *kopt = ((AliITSresponseSDD*)res->GetResponse())->ZeroSuppOption();
      if((!strstr(kopt,"2D"))&&(!strstr(kopt,"1D"))) SetDigitClassName(dettype,"AliITSdigit");
      else SetDigitClassName(dettype,"AliITSdigitSDD");
      SetClusterClassName(dettype,"AliITSRawClusterSDD");

    }
    if(dettype==2){
      AliITSsegmentationSSD* seg2 = new AliITSsegmentationSSD(GetITSgeom());
      seg2->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
      seg2->SetAnglesLay5(0.0075,0.0275); // strip angels rad P and N side.
      seg2->SetAnglesLay6(0.0275,0.0075); // strip angels rad P and N side.
      SetSegmentationModel(dettype,seg2);
      SetDigitClassName(dettype,"AliITSdigitSSD");
      SetClusterClassName(dettype,"AliITSRawClusterSSD");
    }
  }
  
}
//______________________________________________________________________
Bool_t AliITSDetTypeRec::GetCalibration() {
  // Get Default calibration if a storage is not defined.

  if(!fFirstcall){
    AliITSCalibration* cal = GetCalibrationModel(0);
    if(cal)return kTRUE;
  }
  else {
    fFirstcall = kFALSE;
  }

  //  SetRunNumber((Int_t)AliCDBManager::Instance()->GetRun());
  //  Int_t run=GetRunNumber();

  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  if (fCalibration==0) {
    fCalibration = new TObjArray(GetITSgeom()->GetIndexMax());
    fCalibration->SetOwner(!cacheStatus);
    fCalibration->Clear();
  }

  AliCDBEntry *entrySPD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSPD");
  AliCDBEntry *entrySDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD");
  AliCDBEntry *entrySSD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSSD");
  AliCDBEntry *entry2SPD = AliCDBManager::Instance()->Get("ITS/Calib/RespSPD");
  AliCDBEntry *entry2SDD = AliCDBManager::Instance()->Get("ITS/Calib/RespSDD");
  AliCDBEntry *entry2SSD = AliCDBManager::Instance()->Get("ITS/Calib/RespSSD");

  if(!entrySPD || !entrySDD || !entrySSD || !entry2SPD || !entry2SDD || !entry2SSD){
    AliFatal("Calibration object retrieval failed! ");
    return kFALSE;
  }  	

  TObjArray *calSPD = (TObjArray *)entrySPD->GetObject();
  if(!cacheStatus)entrySPD->SetObject(NULL);
  entrySPD->SetOwner(kTRUE);
 
  AliITSresponseSPD *pSPD = (AliITSresponseSPD*)entry2SPD->GetObject();
  if(!cacheStatus)entry2SPD->SetObject(NULL);
  entry2SPD->SetOwner(kTRUE);
    
  TObjArray *calSDD = (TObjArray *)entrySDD->GetObject();
  if(!cacheStatus)entrySDD->SetObject(NULL);
  entrySDD->SetOwner(kTRUE);
 
  AliITSresponseSDD *pSDD = (AliITSresponseSDD*)entry2SDD->GetObject();
  if(!cacheStatus)entry2SDD->SetObject(NULL);
  entry2SDD->SetOwner(kTRUE);

  TObjArray *calSSD = (TObjArray *)entrySSD->GetObject();
  if(!cacheStatus)entrySSD->SetObject(NULL);
  entrySSD->SetOwner(kTRUE);

  AliITSresponseSSD *pSSD = (AliITSresponseSSD*)entry2SSD->GetObject();
  if(!cacheStatus)entry2SSD->SetObject(NULL);
  entry2SSD->SetOwner(kTRUE);

  // DB entries are deleted. In this way metadeta objects are deleted as well
  if(!cacheStatus){
    delete entrySPD;
    delete entrySDD;
    delete entrySSD;
    delete entry2SPD;
    delete entry2SDD;
    delete entry2SSD;
  }

  if ((!pSPD)||(!pSDD)||(!pSSD) || (!calSPD) || (!calSDD) || (!calSSD)) {
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
    cal->SetResponse((AliITSresponse*)pSPD);
    SetCalibrationModel(i, cal);
 }
  for (Int_t i=0; i<fNMod[1]; i++) {
    cal = (AliITSCalibration*) calSDD->At(i);
    cal->SetResponse((AliITSresponse*)pSDD);
    Int_t iMod = i + fNMod[0];
    SetCalibrationModel(iMod, cal);
 }
  for (Int_t i=0; i<fNMod[2]; i++) {
    cal = (AliITSCalibration*) calSSD->At(i);
    cal->SetResponse((AliITSresponse*)pSSD);
    Int_t iMod = i + fNMod[0] + fNMod[1];
    SetCalibrationModel(iMod, cal);
 }

  return kTRUE;
}


//________________________________________________________________
void AliITSDetTypeRec::SetDefaultClusterFinders(){
  
  //set defaults for standard cluster finder

  if(!GetITSgeom()){
    Warning("SetDefaults","null pointer to AliITSgeom!");
    return;
  }

  AliITSClusterFinder *clf; 

  for(Int_t dettype=0;dettype<fgkNdettypes;dettype++){
    //SPD
    if(dettype==0){
      if(!GetReconstructionModel(dettype)){
	TClonesArray *dig0 = DigitsAddress(0);
	TClonesArray *rec0 = ClustersAddress(0);
	clf = new AliITSClusterFinderSPD(this,dig0,rec0);
	SetReconstructionModel(dettype,clf);

      }
    }
   
    //SDD
    if(dettype==1){
      if(!GetReconstructionModel(dettype)){
	TClonesArray *dig1 = DigitsAddress(1);
	TClonesArray *rec1 = ClustersAddress(1);
	clf = new AliITSClusterFinderSDD(this,dig1,rec1);
	SetReconstructionModel(dettype,clf);
      }

    }
    //SSD
    if(dettype==2){
      if(!GetReconstructionModel(dettype)){
	TClonesArray* dig2 = DigitsAddress(2);
	clf = new AliITSClusterFinderSSD(this,dig2);
	SetReconstructionModel(dettype,clf);
      }
    }

 }
 
  
}

//________________________________________________________________
void AliITSDetTypeRec::SetDefaultClusterFindersV2(Bool_t rawdata){

  //Set defaults for cluster finder V2

  if(!GetITSgeom()){
    Warning("SetDefaults","Null pointer to AliITSgeom !");
    return;
  }

  AliITSClusterFinder *clf; 

  for(Int_t dettype=0;dettype<fgkNdettypes;dettype++){
    //SPD
    if(dettype==0){
      if(!GetReconstructionModel(dettype)){
	clf = new AliITSClusterFinderV2SPD(this);
	clf->InitGeometry();
	if(!rawdata) clf->SetDigits(DigitsAddress(0));
	SetReconstructionModel(dettype,clf);

      }
    }
    //SDD
    if(dettype==1){
      if(!GetReconstructionModel(dettype)){
	clf = new AliITSClusterFinderV2SDD(this);
	clf->InitGeometry();
	if(!rawdata) clf->SetDigits(DigitsAddress(1));
	SetReconstructionModel(dettype,clf);
      }

    }

    //SSD
    if(dettype==2){
      if(!GetReconstructionModel(dettype)){
	clf = new AliITSClusterFinderV2SSD(this);
	clf->InitGeometry();
	if(!rawdata) clf->SetDigits(DigitsAddress(2));
	SetReconstructionModel(dettype,clf);
      }
    }

 }
   
}
//______________________________________________________________________
void AliITSDetTypeRec::MakeBranch(TTree* tree, Option_t* option){

  //Creates branches for clusters and recpoints
  Bool_t cR = (strstr(option,"R")!=0);
  Bool_t cRF = (strstr(option,"RF")!=0);
  
  if(cRF)cR = kFALSE;

  if(cR) MakeBranchR(tree);
  if(cRF) MakeBranchRF(tree);

}

//___________________________________________________________________
void AliITSDetTypeRec::AddCluster(Int_t id, AliITSRawCluster *c){

  // Adds a raw cluster to the list
  TClonesArray &lc = *((TClonesArray*)fCtype->At(id));  
  switch(id){
  case 0:
    new(lc[fNctype[id]++]) AliITSRawClusterSPD(*((AliITSRawClusterSPD*)c));
    break;
  case 1:
    new(lc[fNctype[id]++]) AliITSRawClusterSDD(*((AliITSRawClusterSDD*)c));
    break;
  case 2:
    new(lc[fNctype[id]++]) AliITSRawClusterSSD(*((AliITSRawClusterSSD*)c));
    break;
  } 
}
//___________________________________________________________________
void AliITSDetTypeRec::ResetDigits(){
  // Reset number of digits and the digits array for the ITS detector.
  
  if(!fDigits) return;
  for(Int_t i=0;i<fgkNdettypes;i++){
    ResetDigits(i);
  }
}
//___________________________________________________________________
void AliITSDetTypeRec::ResetDigits(Int_t branch){
  // Reset number of digits and the digits array for this branch.
  
  if(fDigits->At(branch)) ((TClonesArray*)fDigits->At(branch))->Clear();
  if(fNdtype) fNdtype[branch]=0;

}

//__________________________________________________________________
void AliITSDetTypeRec::ResetClusters(){

  //Resets number of clusters and the cluster array 
  for(Int_t i=0;i<fgkNdettypes;i++){
    ResetClusters(i);
  }
}

//__________________________________________________________________
void AliITSDetTypeRec::ResetClusters(Int_t i){

  //Resets number of clusters and the cluster array for this branch

  if (fCtype->At(i))    ((TClonesArray*)fCtype->At(i))->Clear();
  if (fNctype)  fNctype[i]=0;
}
//__________________________________________________________________
void AliITSDetTypeRec::MakeBranchR(TTree *treeR, Option_t *opt){

  //Creates tree branches for recpoints
  // Inputs:
  //      cont char *file  File name where RecPoints branch is to be written
  //                       to. If blank it write the SDigits to the same
  //                       file in which the Hits were found.

  Int_t buffsz = 4000;
  char branchname[30];

  // only one branch for rec points for all detector types
  Bool_t oFast= (strstr(opt,"Fast")!=0);
  
  Char_t detname[10] = "ITS";
 
  
  if(oFast){
    sprintf(branchname,"%sRecPointsF",detname);
  } else {
    sprintf(branchname,"%sRecPoints",detname);
  }
  
  if(!fRecPoints)fRecPoints = new TClonesArray("AliITSRecPoint",1000);
  if (treeR)
    MakeBranchInTree(treeR,branchname,0,&fRecPoints,buffsz,99);
}
//______________________________________________________________________
void AliITSDetTypeRec::SetTreeAddressR(TTree *treeR){
    // Set branch address for the Reconstructed points Trees.
    // Inputs:
    //      TTree *treeR   Tree containing the RecPoints.
    // Outputs:
    //      none.
    // Return:

    char branchname[30];
    Char_t namedet[10]="ITS";

    if(!treeR) return;
    if(fRecPoints==0x0) fRecPoints = new TClonesArray("AliITSRecPoint",1000);
    TBranch *branch;
    sprintf(branchname,"%sRecPoints",namedet);
    branch = treeR->GetBranch(branchname);
    if (branch) {
      branch->SetAddress(&fRecPoints);
    }else {
      sprintf(branchname,"%sRecPointsF",namedet);
      branch = treeR->GetBranch(branchname);
      if (branch) {
	branch->SetAddress(&fRecPoints);
      }

    }
}
//____________________________________________________________________
void AliITSDetTypeRec::AddRecPoint(const AliITSRecPoint &r){
    // Add a reconstructed space point to the list
    // Inputs:
    //      const AliITSRecPoint &r RecPoint class to be added to the tree
    //                              of reconstructed points TreeR.
    // Outputs:
    //      none.
    // Return:
    //      none.

    TClonesArray &lrecp = *fRecPoints;
    new(lrecp[fNRecPoints++]) AliITSRecPoint(r);
}

//______________________________________________________________________
void AliITSDetTypeRec::DigitsToRecPoints(TTree *treeD,TTree *treeR,Int_t lastentry,Option_t *opt, Bool_t v2){
  // cluster finding and reconstruction of space points
  // the condition below will disappear when the geom class will be
  // initialized for all versions - for the moment it is only for v5 !
  // 7 is the SDD beam test version
  // Inputs:
  //      TTree *treeD     Digits tree
  //      TTree *treeR     Clusters tree
  //      Int_t lastentry  Offset for module when not all of the modules
  //                       are processed.
  //      Option_t *opt    String indicating which ITS sub-detectors should
  //                       be processed. If ="All" then all of the ITS
  //                       sub detectors are processed.

  const char *all = strstr(opt,"All");
  const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
                        strstr(opt,"SSD")};
  if(!v2) {
    SetDefaultClusterFinders();
    AliInfo("Original cluster finder has been selected\n");
  }
  else   { 
    SetDefaultClusterFindersV2();
    AliInfo("V2 cluster finder has been selected \n");
  }

  AliITSClusterFinder *rec     = 0;
  Int_t id,module,first=0;
  for(module=0;module<GetITSgeom()->GetIndexMax();module++){
      id       = GetITSgeom()->GetModuleType(module);
      if (!all && !det[id]) continue;
      if(det[id]) first = GetITSgeom()->GetStartDet(id);
      rec = (AliITSClusterFinder*)GetReconstructionModel(id);
      TClonesArray *itsDigits  = DigitsAddress(id);
      if (!rec)
          AliFatal("The reconstruction class was not instanciated!");
      ResetDigits();  // MvL: Not sure we neeed this when rereading anyways
      if (all) {
          treeD->GetEvent(lastentry+module);
      }else {
          treeD->GetEvent(lastentry+(module-first));
      }
      Int_t ndigits = itsDigits->GetEntriesFast();
      if(ndigits>0){
	rec->SetDetTypeRec(this);
	rec->SetDigits(DigitsAddress(id));
	//	rec->SetClusters(ClustersAddress(id));
	rec->FindRawClusters(module);
      } // end if
      treeR->Fill();
      ResetRecPoints();
      ResetClusters();
  } 
}
//______________________________________________________________________
void AliITSDetTypeRec::DigitsToRecPoints(AliRawReader* rawReader,TTree *treeR,Option_t *opt){
  // cluster finding and reconstruction of space points
  // the condition below will disappear when the geom class will be
  // initialized for all versions - for the moment it is only for v5 !
  // 7 is the SDD beam test version
  // Inputs:
  //      AliRawReader *rawReader  Pointer to the raw-data reader
  //      TTree *treeR             Clusters tree
  // Outputs:
  //      none.
  // Return:
  //      none.
  const char *all = strstr(opt,"All");
  const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
                        strstr(opt,"SSD")};
  AliITSClusterFinderV2 *rec     = 0;
  Int_t id=0;

  TClonesArray *array=new TClonesArray("AliITSRecPoint",1000);
  TBranch *branch = treeR->Branch("ITSRecPoints",&array);
  delete array;
 
  TClonesArray** clusters = new TClonesArray*[GetITSgeom()->GetIndexMax()]; 
  for (Int_t iModule = 0; iModule < GetITSgeom()->GetIndexMax(); iModule++) {
    clusters[iModule] = NULL;
  }
  for(id=0;id<3;id++){
    if (!all && !det[id]) continue;
    rec = (AliITSClusterFinderV2*)GetReconstructionModel(id);
    if (!rec)
      AliFatal("The reconstruction class was not instanciated");
    rec->SetDetTypeRec(this);
    rec->RawdataToClusters(rawReader,clusters);    
  } 
  Int_t nClusters =0;
  TClonesArray *emptyArray=new TClonesArray("AliITSRecPoint");
  for(Int_t iModule=0;iModule<GetITSgeom()->GetIndexMax();iModule++){
    id = GetITSgeom()->GetModuleType(iModule);
    if (!all && !det[id]) continue;
    array = clusters[iModule];
    if(!array){
      AliDebug(1,Form("data for module %d missing!",iModule));
      array = emptyArray;
    }
    branch->SetAddress(&array);
    treeR->Fill();
    nClusters+=array->GetEntriesFast();

    if (array != emptyArray) {
      array->Delete();
      delete array;
    }
  }
  delete emptyArray;

  delete[] clusters;
  Info("DigitsToRecPoints", "total number of found recpoints in ITS: %d\n", 
       nClusters);
  
}


