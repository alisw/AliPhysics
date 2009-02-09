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

/*
$Log$
Revision 1.1  2007/10/01 14:12:45  pchrist
Class creating the aligmnent object fro the surveyor measurements.

*/ 

//
//  Creates the SSD align object
//
//
#include "Riostream.h"
#include "TFile.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoPhysicalNode.h"

#include "AliITSSurveyToAlignSSD.h"
#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"
#include "AliAlignObjParams.h"

#include "AliLog.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"

ClassImp(AliITSSurveyToAlignSSD)

//________________________________________________________________________
AliITSSurveyToAlignSSD::AliITSSurveyToAlignSSD() :
  TObject(),
  fRun(0),
  fFileLocal(0x0),
  fFileGrid(0x0),
  fSurveyPoints(0),
  fSSDAlignObj(new TClonesArray("AliAlignObjParams",100)),
  fSSDAlignObjParam(0),
  fDebug(0) {
  //
  //  default constructor
  //
}   

//________________________________________________________________________
AliITSSurveyToAlignSSD::AliITSSurveyToAlignSSD(Int_t reportId, Int_t runNumber) :
  TObject(),
  fRun(0),
  fFileLocal(0x0),
  fFileGrid(0x0),
  fSurveyPoints(0),
  fSSDAlignObj(new TClonesArray("AliAlignObjParams",100)),
  fSSDAlignObjParam(0),
  fDebug(0) {
  //
  // constructor - defines data files
  //
  fRun = runNumber;
  fFileLocal = new Char_t[80];
  fFileGrid = new Char_t[80];
  Char_t path[50];

  sprintf(path,gSystem->Getenv("ALICE_ROOT")); 
  sprintf(fFileLocal,"%s/ITS/Survey_SSD_%d.txt",path,reportId);  
}

//_________________________________________________________________________
AliITSSurveyToAlignSSD::AliITSSurveyToAlignSSD(const AliITSSurveyToAlignSSD &align) :
  TObject(),
  fRun(0),
  fFileLocal(0x0),
  fFileGrid(0x0),
  fSurveyPoints(0),
  fSSDAlignObj(new TClonesArray("AliAlignObjParams",100)),
  fSSDAlignObjParam(0),
  fDebug(0) {
  //
  //  copy constructor - dummy
  
  fDebug = align.fDebug;
}

//__________________________________________________________________________
AliITSSurveyToAlignSSD & AliITSSurveyToAlignSSD::operator =(const AliITSSurveyToAlignSSD & align) {
  //
  // assignment operator - dummy
  //
  fDebug = align.fDebug;
  return (*this);
}

//__________________________________________________________________________
AliITSSurveyToAlignSSD::~AliITSSurveyToAlignSSD() {
  //
  // destructor
  //
  if(fSurveyPoints) delete fSurveyPoints;
  if(fSSDAlignObj) delete fSSDAlignObj;
  if(fSSDAlignObjParam) delete fSSDAlignObjParam;
}

//______________________________________________________________________
void AliITSSurveyToAlignSSD::Run() { 
  //
  // runs the full chain
  //
  SetDebug(0);
  Bool_t flag = LoadSurveyData();
  if(!flag) {
    cout<<"Missing points"<<endl;
    return;
  }
  CreateAlignObj();
  StoreAlignObj();
}

//__________________________________________________________________________
Bool_t AliITSSurveyToAlignSSD::LoadSurveyData() {
  //
  // for a time being it loads from the local file the surveyed point
  // and has the ideal points hardwired. I am waiting until Ricardo
  // completes his job
  // 
  
  //Load survey data from the local file
  AliSurveyObj * s1 = new AliSurveyObj();
  if(s1->FillFromLocalFile(fFileLocal))
    fSurveyPoints = s1->GetData();
  else 
    return kFALSE;
  
  cout<<"Number of SSD survey points: "<<fSurveyPoints->GetEntries()<<endl;  
  
  return kTRUE;
}

//_______________________________________________________________________
void AliITSSurveyToAlignSSD::CreateAlignObj(){
  //
  // This method creates AliAlignObj and fills it with Euler angles
  // and shifts. The units are degrees and cm.
  // 
  TClonesArray &alobj = *fSSDAlignObj;

  //load ideal geometry from the OCDB
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(fRun);
  AliCDBEntry* entry = cdb->Get("GRP/Geometry/Data");
  AliGeomManager::SetGeometry(gGeoManager);

  AliSurveyPoint* pt = 0;
  TGeoHMatrix* hm = 0;
  TGeoPNEntry* pne = 0;
  Double_t* tr;
  Int_t uid;
  const char* symname;
  Double_t sx, sy, sz;
  Int_t ilayer, imodule;

  for(Int_t imod = 0; imod < fSurveyPoints->GetEntries(); imod++) {
    pt = (AliSurveyPoint*) fSurveyPoints->At(imod);
    if(!pt) continue;
    sx = pt->GetX();
    sy = pt->GetY();
    sz = pt->GetZ();

    ilayer = (imod < 748) ? AliGeomManager::kSSD1 : AliGeomManager::kSSD2;
    imodule = (imod < 748) ? imod : imod - 748;

    uid = AliGeomManager::LayerToVolUID(ilayer,imodule);
    symname = AliGeomManager::SymName(uid);
    pne = gGeoManager->GetAlignableEntryByUID(uid);
    hm = pne->GetGlobalOrig();
    tr = hm->GetTranslation();

    //Printf("symname %s", symname);
    //Printf("x,y,z from survey: %f %f %f", sx, sy, sz);
    //Printf("x,y,z from ideal : %f %f %f", tr[0], tr[1], tr[2]);
    
   new(alobj[imod]) AliAlignObjParams(symname, uid, sx-tr[0], sy-tr[1], sz-tr[2], 0., 0., 0., kTRUE);
  }//module loop

  delete entry;
}

//_________________________________________________________________________
void AliITSSurveyToAlignSSD::StoreAlignObj(){
  // Stores the TClonesArray of AliAlignObj in
  // $ALICE_ROOT/ITS/Align/Data/SSDfromSurvey.root
  // In a later version these objects will be merged 
  // with the SDD ones and will be put in the OCDB
  const char* filename = "$ALICE_ROOT/ITS/Align/Data/SSDfromSurvey.root";
  TFile *f = TFile::Open(filename,"RECREATE");
  if(!f){
    Error(filename,"cannot open file for output\n");
    return;
  }
  cout<<"Saving alignment objects to the file "<<filename<<endl;
  f->cd();
  f->WriteObject(fSSDAlignObj,"ITSAlignObjs","kSingleKey");
  f->Close();
}
