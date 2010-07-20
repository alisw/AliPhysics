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

#include "Riostream.h"
#include "TFile.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include "AliSurveyToAlignObjs.h"
#include "AliSurveyPoint.h"
#include "AliAlignObjParams.h"

#include "AliLog.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"

ClassImp(AliSurveyToAlignObjs)

//________________________________________________________________________
AliSurveyToAlignObjs::AliSurveyToAlignObjs() :
  TObject(),
  fSurveyObj(NULL),
  fSurveyPoints(NULL),
  fAlignObjArray(NULL),
  fAlignObj(0){
  //
  //  default constructor
  fSurveyObj = new AliSurveyObj();
  fAlignObjArray = new TClonesArray("AliAlignObjParams",10);
}   

//_________________________________________________________________________
AliSurveyToAlignObjs::AliSurveyToAlignObjs(const AliSurveyToAlignObjs &s2aObj) :
  TObject(s2aObj),
  fSurveyObj(s2aObj.fSurveyObj),
  fSurveyPoints(s2aObj.fSurveyPoints),
  fAlignObjArray(s2aObj.fAlignObjArray),
  fAlignObj(s2aObj.fAlignObj)
{
  // copy constructor
}

//__________________________________________________________________________
AliSurveyToAlignObjs & AliSurveyToAlignObjs::operator= (const AliSurveyToAlignObjs &s2aObj) {
  //
  // assignment operator
  fSurveyObj = s2aObj.fSurveyObj;
  fSurveyPoints = s2aObj.fSurveyPoints;
  fAlignObjArray = s2aObj.fAlignObjArray;
  fAlignObj = s2aObj.fAlignObj;
  return (*this);
}

//__________________________________________________________________________
AliSurveyToAlignObjs::~AliSurveyToAlignObjs() {
  //
  // destructor
  //
  if(fSurveyObj) delete fSurveyObj;
  if(fSurveyPoints) delete fSurveyPoints;
  if(fAlignObjArray) delete fAlignObjArray;
  if(fAlignObj) delete fAlignObj;
}

//__________________________________________________________________________
Bool_t AliSurveyToAlignObjs::LoadSurveyFromLocalFile(const char* filename) {
  // Load survey data from a formatted text file
  // residing locally
  //
  
  //Load survey data from the local file
  if(fSurveyObj->FillFromLocalFile(filename))
    fSurveyPoints = fSurveyObj->GetData();
  else 
    return kFALSE;
  
  AliInfo(Form("%d survey points read",fSurveyPoints->GetEntries()));  
  
  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliSurveyToAlignObjs::LoadSurveyFromAlienFile(const char* det, Int_t repNum, Int_t repVersion) {
  // Load survey data from the formatted text file
  // residing in the default location in alien
  //
  
  const char* alienUser=gSystem->Getenv("alien_API_USER");
  if(fSurveyObj->Fill(det, repNum, repVersion, alienUser))
  {
    fSurveyPoints = fSurveyObj->GetData();
  }else{
    AliError("Error reading survey file from alien!");
    return kFALSE;
  }
  
  AliInfo(Form("%d survey points read",fSurveyPoints->GetEntries()));  
  
  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliSurveyToAlignObjs::StoreAlignObjToFile(const char* filename, const char* det){
  // Stores the TClonesArray of alignment objects into the
  // file specified as argument
  //
  TFile *f = TFile::Open(filename,"RECREATE");
  if(!f){
    AliError(Form("cannot open file %s\n",filename));
    return kFALSE;
  }
  AliInfo(Form("Saving alignment objects into the file %s",filename));
  TString arrayname(det);
  arrayname+="AlignObjs";
      
  f->cd();
  f->WriteObject(fAlignObjArray,arrayname,"kSingleKey");
  f->Close();

  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliSurveyToAlignObjs::StoreAlignObjToCDB(const char* cdbFolder, const char* det){
  // Stores the TClonesArray of alignment objects into a
  // CDB entry in the CDB folder specified by the argument
  //

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbFolder);
  cdb->SetRun(0);

  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetComment(Form("Misalignment for subdetector %s from survey",det));
  TString path(det);
  path+="/Align/Data";
  AliCDBId id(path.Data(),0,AliCDBRunRange::Infinity());
  cdb->Put(fAlignObjArray,id,md);

  return kTRUE;
}


