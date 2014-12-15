
// TTimeStamp startTime(2006,10,18,0,0,0,0,kFALSE)
// TTimeStamp endTime(2006,10,19,0,0,0,0,kFALSE)
// Int_t run=2546
// AliTPCGenDBTemp db
// db->Init(run,"TPC/Config/Temperature","TPC/*/*")
// db->MakeCalib("TempSensor.txt","DCSMap.root",startTime,endTime,run)

//  Data base entry generation:
  
//  AliTPCGenDBTemp db
//  db->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//  db->SetSpecificStorage("local:///afs/cern.ch/alice/tpctest/Calib/");
//  db->Init(0,"TPC/Config/Temperature","TPC/*/*")
//  db->MakeConfig("TempSensor.txt",0,999999999,"TPC/Config/Temperature")

  




#include "AliTPCGenDBTemp.h"
#include "AliLog.h"

ClassImp(AliTPCGenDBTemp)

const Int_t kValCut = 100;         // discard temperatures > 100 degrees
const Int_t kDiffCut = 5;	   // discard temperature differences > 5 degrees

//______________________________________________________________________________________________

AliTPCGenDBTemp::AliTPCGenDBTemp():
   AliDCSGenDB()
{
}

//______________________________________________________________________________________________

AliTPCGenDBTemp::AliTPCGenDBTemp(const char *defaultStorage, const char *specificStorage) :
   AliDCSGenDB(defaultStorage,specificStorage)
{
}

//______________________________________________________________________________________________

AliTPCGenDBTemp::AliTPCGenDBTemp(const AliTPCGenDBTemp& ) : AliDCSGenDB()
{

//
//  Copy constructor
//
 AliError("copy constructor not implemented");

}

//______________________________________________________________________________________________
AliTPCGenDBTemp::~AliTPCGenDBTemp(){
//
// destructor
//

}
//______________________________________________________________________________________________
AliTPCGenDBTemp& AliTPCGenDBTemp::operator= (const AliTPCGenDBTemp&  )
{
 //
 // assignment operator
 //
 AliError("assignment operator not implemented");
 return *this;
}


//______________________________________________________________________________________________

void AliTPCGenDBTemp::MakeCalib(const char *fList, const char *fMap,
                             const TTimeStamp& startTime,
			     const TTimeStamp& endTime,
			     Int_t run, const TString& amandaString )
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   AliTPCSensorTempArray *temperature=0;
   if ( amandaString.Length()== 0 ) {
    temperature = new AliTPCSensorTempArray(fList);
   } else {
    temperature = new AliTPCSensorTempArray(fList,amandaString);
   }
   temperature->SetStartTime(startTime);
   temperature->SetEndTime(endTime);
   temperature->SetValCut(kValCut);
   temperature->SetDiffCut(kDiffCut);
   TMap* map = SetGraphFile(fMap);
   if (map) {
     temperature->MakeSplineFit(map);
   }
   delete map;
   map=0;
   fMap=0;

   SetFirstRun(run);
   SetLastRun(run);
   SetSensorArray(temperature);
   StoreObject("TPC/Calib/Temperature",temperature, fMetaData);
}

//______________________________________________________________________________________________

TClonesArray * AliTPCGenDBTemp::ReadList(const char *fname, const char *title,
                       const TString& amandaString) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree(title,title);
  tree->ReadFile(fname,"");
  TClonesArray *arr;
  if ( amandaString.Length()== 0 ) {
    arr = AliTPCSensorTemp::ReadTree(tree);
  } else {
    arr = AliTPCSensorTemp::ReadTree(tree,amandaString);
  }
  delete tree;
  return arr;
}

//______________________________________________________________________________________________

TTree * AliTPCGenDBTemp::ReadListTree(const char *fname, const char *title) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree(title,title);
  tree->ReadFile(fname,"");
  return tree;
}

//______________________________________________________________________________________________
void AliTPCGenDBTemp::MakeConfig(const char *file, Int_t firstRun, Int_t lastRun, 
                             const char *confDir)
{
   //
   // Store Configuration file to OCDB
   //

   TTree *tree = ReadListTree(file,"tempConf");
   SetConfTree(tree);
   SetFirstRun(firstRun);
   SetLastRun(lastRun);

   StoreObject(confDir, fConfTree, fMetaData);
}


