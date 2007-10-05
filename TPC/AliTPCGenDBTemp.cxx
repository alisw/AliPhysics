
// .L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCGenDBTemp.C+
// TTimeStamp startTime(2006,10,18,0,0,0,0,kFALSE)
// TTimeStamp endTime(2006,10,19,0,0,0,0,kFALSE)
// Int_t run=2546
// AliTPCGenDBTemp db
// db->Init(run,"TPC/Config/Temperature","TPC/*/*")
// db->MakeCalib("TempSensor.txt","DCSMap.root",startTime,endTime,run)


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

AliTPCGenDBTemp::AliTPCGenDBTemp(const AliTPCGenDBTemp& org) : AliDCSGenDB(org)
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
AliTPCGenDBTemp& AliTPCGenDBTemp::operator= (const AliTPCGenDBTemp& org )
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
			     Int_t run )
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   AliTPCSensorTempArray *temperature = new AliTPCSensorTempArray(fList);
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

TClonesArray * AliTPCGenDBTemp::ReadList(const char *fname) {
  //
  // read values from ascii file
  //
  TTree* tree = new TTree("tempConf","tempConf");
  tree->ReadFile(fname,"");
  TClonesArray *arr = AliTPCSensorTemp::ReadTree(tree);
  delete tree;
  return arr;
}

//______________________________________________________________________________________________

TTree * AliTPCGenDBTemp::ReadListTree(const char *fname) {
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
