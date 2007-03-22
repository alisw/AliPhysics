/**
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliSplineFit.cxx+
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliDCSSensor.cxx+
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliDCSSensorArray.cxx+
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCSensorTemp.cxx+
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCSensorTempArray.cxx+
.L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCDBTemp.C+
TTimeStamp startTime(2006,10,18,0,0,0,0,kFALSE)
TTimeStamp endTime(2006,10,19,0,0,0,0,kFALSE)
Int_t run=2546
AliTPCDBTemp db
db->Init(run)
db->MakeCalib("TempSensor.txt","DCSMap.root",startTime,endTime,run)


**/
#include "AliTPCDBTemp.h"

AliTPCDBTemp::AliTPCDBTemp(): 
   fFirstRun(0),
   fLastRun(0),
   fTemperature(0),
   fStorLoc(0),
   fCalib(0),
   fMetaData(0)
{}



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
   TMap* map = SetGraphFile(fMap);
   if (map) {
     fTemperature->MakeSplineFit(map);
   }
   delete map;

   SetFirstRun(run);
   SetLastRun(run);   		    
   StoreObject("TPC/Calib/Temperature",fTemperature, fMetaData);
}


AliCDBMetaData* AliTPCDBTemp::CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Haavard Helstrup");
  md1->SetBeamPeriod(2);
  md1->SetAliRootVersion("05-13-04"); //root version
  md1->SetComment("Temperature values");
  
  return md1;
}

void AliTPCDBTemp::StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{

  AliCDBId id1(cdbPath, fFirstRun, fLastRun); 
  if (fStorLoc) fStorLoc->Put(object, id1, metaData); 
}

void AliTPCDBTemp::Init(Int_t run){

//   Int_t kLastRun=4000;
//   Long64_t longRun;
   
   SetFirstRun(run);
   SetLastRun(run); 
       
   InitDB(run);
//   fCalib = AliTPCcalibDB::Instance();    
//   longRun=run;
//   fCalib->SetRun(longRun);
//   fTemperature = fCalib->GetTemperature();
     
}

void AliTPCDBTemp::InitDB(Int_t run)
{ 
   //   Data base generation
   
//   printf ("Data base creation started.. \n");
   char   *CDBpath="local:///afs/cern.ch/alice/tpctest/Calib/";

   fMetaData = CreateMetaObject("AliTPCSensorTempArray");
   AliCDBManager *man = AliCDBManager::Instance();
   man->SetDefaultStorage("local:///afs/cern.ch/alice/tpctest/AliRoot/HEAD"); 
   man->SetRun(run);
   man->SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/alice/tpctest/Calib");
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
