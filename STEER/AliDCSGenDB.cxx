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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class to generate DCS data base entries                                  //
//  Author: Haavard Helstrup                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////





// TTimeStamp startTime(2006,10,18,0,0,0,0,kFALSE)
// TTimeStamp endTime(2006,10,19,0,0,0,0,kFALSE)
// Int_t run=2546
// AliDCSGenDB db
// db->Init(run,"TPC/Config/Temperature","TPC/*/*")
// db->MakeCalib("PressureSensor.txt","DCSMap.root",startTime,endTime,firstRun,lastRun,"TPC/Calib/Temperature")


#include "AliDCSGenDB.h"

const char *kDefaultStorage="local:///afs/cern.ch/alice/tpctest/AliRoot/HEAD";
const char *kSpecificStorage="local:///afs/cern.ch/alice/tpctest/Calib/";
const char *kSensorClass = "AliDCSSensorArray";
const Int_t kBeamPeriod=2;

ClassImp(AliDCSGenDB)

//______________________________________________________________________________________________

AliDCSGenDB::AliDCSGenDB():
   fFirstRun(0),
   fLastRun(0),
   fSpecificStorage(kSpecificStorage),
   fDefaultStorage(kDefaultStorage),
   fSensor(0),
   fStorLoc(0),
   fCalib(0),
   fMetaData(0),
   fConfTree(0)
//
//  standard constructor
//
{}

//______________________________________________________________________________________________

AliDCSGenDB::AliDCSGenDB(const AliDCSGenDB& org):
  TObject(org),
  fFirstRun(org.fFirstRun),
  fLastRun(org.fLastRun),
  fSpecificStorage(org.fSpecificStorage),
  fDefaultStorage(org.fDefaultStorage),
  fSensor(0),
  fStorLoc(0),
  fCalib(0),
  fMetaData(0),
  fConfTree(0)
{
//
//  Copy constructor
//

 ((AliDCSGenDB &) org).Copy(*this);
}

//______________________________________________________________________________________________
AliDCSGenDB::~AliDCSGenDB(){
//
// destructor
//
   fCalib->Terminate();
   delete fSensor;
   delete fMetaData;
   delete fConfTree;
}

//______________________________________________________________________________________________
AliDCSGenDB& AliDCSGenDB::operator= (const AliDCSGenDB& org )
{
 //
 // assignment operator
 //
 if (&org == this) return *this;

 new (this) AliDCSGenDB(org);
 return *this;
}

//______________________________________________________________________________________________
void AliDCSGenDB::Copy(TObject &c) const
{
  //
  // Copy function
  //

  TObject::Copy(c);
}

//______________________________________________________________________________________________
void AliDCSGenDB::MakeCalib(const char *list, const char *mapDCS,
                             const TTimeStamp& startTime,
			     const TTimeStamp& endTime,
			     Int_t firstRun, Int_t lastRun, const char *calibDir )
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TClonesArray *arr = ReadList(list);
   AliDCSSensorArray *fSensor = new AliDCSSensorArray(arr);
   fSensor->SetStartTime(startTime);
   fSensor->SetEndTime(endTime);
   TMap* map = SetGraphFile(mapDCS);
   if (map) {
     fSensor->MakeSplineFit(map);
   }
   delete map;
   map=0;
   mapDCS=0;

   SetFirstRun(firstRun);
   SetLastRun(lastRun);

   StoreObject(calibDir, fSensor, fMetaData);
}

//______________________________________________________________________________________________
void AliDCSGenDB::MakeConfig(const char *file, Int_t firstRun, Int_t lastRun, const char *confDir )
{
   //
   // Store Configuration file to OCDB
   //

   TTree *tree = ReadListTree(file);
   SetConfTree(tree);
   SetFirstRun(firstRun);
   SetLastRun(lastRun);

   StoreObject(confDir, fConfTree, fMetaData);
}




//______________________________________________________________________________________________
AliCDBMetaData* AliDCSGenDB::CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData();
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Haavard Helstrup");
  md1->SetBeamPeriod(kBeamPeriod);
  md1->SetAliRootVersion(gSystem->Getenv("ARVERSION"));

  return md1;
}

//______________________________________________________________________________________________
void AliDCSGenDB::StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{

  AliCDBId id1(cdbPath, fFirstRun, fLastRun);
  if (fStorLoc) fStorLoc->Put(object, id1, metaData);
}

//______________________________________________________________________________________________
void AliDCSGenDB::Init(Int_t run, const char *configDir, const char *specificDir)
{

   fMetaData = CreateMetaObject(kSensorClass);
   AliCDBManager *man = AliCDBManager::Instance();
   man->SetDefaultStorage(fDefaultStorage);
   man->SetRun(run);
   man->SetSpecificStorage(specificDir,fSpecificStorage);
   AliCDBEntry *config = man->Get(configDir);
   if (config) fConfTree = (TTree*)config->GetObject();
   fStorLoc = man->GetStorage(fSpecificStorage);
   if (!fStorLoc)    return;

   fCalib = AliTPCcalibDB::Instance();

}

//______________________________________________________________________________________________


//_____________________________________________________________________________
TMap* AliDCSGenDB::SetGraphFile(const char *fname)
{
  //
  // Read DCS maps from file given by fname
  //
  TFile file(fname);
  TMap * map = (TMap*)file.Get("DCSMap");
  return map;
}

//______________________________________________________________________________________________

TClonesArray * AliDCSGenDB::ReadList(const char *fname, const char *title) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree(title,title);
  tree->ReadFile(fname,"");
  TClonesArray *arr = AliDCSSensor::ReadTree(tree);
  delete tree;
  return arr;
}

//______________________________________________________________________________________________

TTree * AliDCSGenDB::ReadListTree(const char *fname, const char *title) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree(title,title);
  tree->ReadFile(fname,"");
  TClonesArray *arr = AliDCSSensor::ReadTree(tree);
  arr->Delete();
  delete arr;
  return tree;
}



