
// .L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCGenDBConf.C+
// Int_t run=2546
// Int_t firstRun =0
// Int_t lastRun = 9999999
// AliTPCGenDBConf db
// db->Init(run,"TPC/Config/Preprocessor","TPC/*/*")
// db->MakeConfig("Preprocessor.txt",firstRun,lastRun,"TPC/Config/Preprocessor")


#include "AliTPCGenDBConf.h"
#include "AliLog.h"

ClassImp(AliTPCGenDBConf)


//______________________________________________________________________________________________

AliTPCGenDBConf::AliTPCGenDBConf():
   AliDCSGenDB()
{
}

//______________________________________________________________________________________________

AliTPCGenDBConf::AliTPCGenDBConf(const char *defaultStorage, const char *specificStorage):
   AliDCSGenDB(defaultStorage, specificStorage)
{
}

//______________________________________________________________________________________________

AliTPCGenDBConf::AliTPCGenDBConf(const AliTPCGenDBConf& org) : AliDCSGenDB(org)
{

//
//  Copy constructor
//
 AliError("copy constructor not implemented");
}

//______________________________________________________________________________________________
AliTPCGenDBConf::~AliTPCGenDBConf(){
//
// destructor
//

}
//______________________________________________________________________________________________
AliTPCGenDBConf& AliTPCGenDBConf::operator= (const AliTPCGenDBConf& /*org*/ )
{
 //
 // assignment operator
 //
 AliError("assignment operator not implemented");
 return *this;
}


//______________________________________________________________________________________________

void AliTPCGenDBConf::MakeConfig(const char *file, Int_t firstRun, Int_t lastRun, const char *confDir )
{
   //
   // Store Configuration file to OCDB
   //

   TEnv *confEnv = new TEnv(file);
   SetFirstRun(firstRun);
   SetLastRun(lastRun);

   StoreObject(confDir, confEnv, fMetaData);
}

//______________________________________________________________________________________________

