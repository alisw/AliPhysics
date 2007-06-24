/**
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCDBTemp.C+
TTimeStamp startTime(2006,10,18,0,0,0,0,kFALSE)
TTimeStamp endTime(2006,10,19,0,0,0,0,kFALSE)
Int_t run=2546
AliTPCDBTemp db
db->Init(run)
db->MakeCalib("TempSensor.txt","DCSMap.root",startTime,endTime,run)


**/
#include "AliTPCDBTemp.h"

ClassImp(AliTPCDBTemp)

const Int_t kValCut = 100;         // discard temperatures > 100 degrees
const Int_t kDiffCut = 5;	   // discard temperature differences > 5 degrees

//______________________________________________________________________________________________

AliTPCDBTemp::AliTPCDBTemp(): 
   fFirstRun(0),
   fLastRun(0),
   fTemperature(0),
   fStorLoc(0),
   fCalib(0),
   fMetaData(0),
   fConfTree(0)
{}
//______________________________________________________________________________________________

AliTPCDBTemp::AliTPCDBTemp(const AliTPCDBTemp& org):
  TObject(org),
  fFirstRun(org.fFirstRun),
  fLastRun(org.fLastRun),
  fTemperature(0),
  fStorLoc(0),
  fCalib(0),
  fMetaData(0),
  fConfTree(0)
{
//
//  Copy constructor
//

 ((AliTPCDBTemp &) org).Copy(*this);
}

//______________________________________________________________________________________________
AliTPCDBTemp::~AliTPCDBTemp(){
//
// destructor
//
   fCalib->Terminate();
   delete fTemperature;
   delete fMetaData;
   delete fConfTree;
}

//______________________________________________________________________________________________
AliTPCDBTemp& AliTPCDBTemp::operator= (const AliTPCDBTemp& org )
{
 //
 // assignment operator
 //
 if (&org == this) return *this;

 new (this) AliTPCDBTemp(org);
 return *this;
} 

//______________________________________________________________________________________________
void AliTPCDBTemp::Copy(TObject &c) const
{
  //
  // Copy function
  //

  TObject::Copy(c);
}


//______________________________________________________________________________________________

void AliTPCDBTemp::MakeCalib(const char *fList, const char *fMap,
                             const TTimeStamp& startTime, 
			     const TTimeStamp& endTime,
			     Int_t run )
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   AliTPCSensorTempArray *fTemperature = new AliTPCSensorTempArray(fList);
   fTemperature->SetStartTime(startTime);
   fTemperature->SetEndTime(endTime);
   fTemperature->SetValCut(kValCut);
   fTemperature->SetDiffCut(kDiffCut);
   TMap* map = SetGraphFile(fMap);
   if (map) {
     fTemperature->MakeSplineFit(map);
   }
   delete map;
   map=0;
   fMap=0;

   SetFirstRun(run);
   SetLastRun(run);   		    
   StoreObject("TPC/Calib/Temperature",fTemperature, fMetaData);
}

//______________________________________________________________________________________________
void AliTPCDBTemp::MakeConfig(const char *file, Int_t firstRun, Int_t lastRun )
{
   //
   // Store Configuration file to OCDB
   //

   TTree *tree = ReadListTree(file);
   SetConfTree(tree);
   SetFirstRun(firstRun);
   SetLastRun(lastRun);   		    
   
   AliCDBMetaData* metaConf=CreateMetaObject("TTree");      
   StoreObject("TPC/Config/Temperature",fConfTree, metaConf);
}


//______________________________________________________________________________________________

AliCDBMetaData* AliTPCDBTemp::CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Haavard Helstrup");
  md1->SetBeamPeriod(2);
  md1->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md1->SetComment("Temperature");
  
  return md1;
}
//______________________________________________________________________________________________

void AliTPCDBTemp::StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{

  AliCDBId id1(cdbPath, fFirstRun, fLastRun); 
  if (fStorLoc) fStorLoc->Put(object, id1, metaData); 
}
//______________________________________________________________________________________________

void AliTPCDBTemp::Init(Int_t run){

//   Int_t kLastRun=4000;
   Long64_t longRun;
   
   SetFirstRun(run);
   SetLastRun(run); 
       
   InitDB(run);
   fCalib = AliTPCcalibDB::Instance();    
   longRun=run;
   fCalib->SetRun(longRun);
   fTemperature = fCalib->GetTemperature();
     
}
//______________________________________________________________________________________________

void AliTPCDBTemp::InitDB(Int_t run)
{ 
   //   Data base generation
   
   char   *CDBpath="local:///afs/cern.ch/alice/tpctest/Calib/";

   fMetaData = CreateMetaObject("AliTPCSensorTempArray");
   AliCDBManager *man = AliCDBManager::Instance();
   man->SetDefaultStorage("local:///afs/cern.ch/alice/tpctest/AliRoot/HEAD"); 
   man->SetRun(run);
   man->SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/alice/tpctest/Calib");
   AliCDBEntry *config = man->Get("TPC/Config/Temperature");
   if (config) fConfTree = (TTree*)config->GetObject();
   fStorLoc = man->GetStorage(CDBpath);
   if (!fStorLoc)    return;
}
//_____________________________________________________________________________
TMap* AliTPCDBTemp::SetGraphFile(const char *fname)
{
  // 
  // Read DCS maps from file given by fname 
  //
  TFile file(fname);
  TMap * map = (TMap*)file.Get("DCSMap");
  return map;
}
//______________________________________________________________________________________________

TClonesArray * AliTPCDBTemp::ReadList(const char *fname) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree("tempConf","tempConf");
  tree->ReadFile(fname,"");
  TClonesArray *arr = AliTPCSensorTemp::ReadTree(tree);
  return arr;
}

//______________________________________________________________________________________________

TTree * AliTPCDBTemp::ReadListTree(const char *fname) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree("tempConf","tempConf");
  tree->ReadFile(fname,"");
  TClonesArray *arr = AliTPCSensorTemp::ReadTree(tree);
  arr->Delete();
  delete arr;
  return tree;
}
