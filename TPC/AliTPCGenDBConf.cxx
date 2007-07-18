
// .L /afs/cern.ch/user/h/haavard/alice/tpc/temperature/AliTPCGenDBConf.C+
// Int_t run=2546
// Int_t firstRun =0
// Int_t lastRun = 9999999
// AliTPCGenDBConf db
// db->Init(run,"TPC/Config/Preprocessor","TPC/*/*")
// db->MakeConfig("Preprocessor.txt",firstRun,lastRun,"TPC/Config/Preprocessor")


#include "AliTPCGenDBConf.h"

ClassImp(AliTPCGenDBConf)


//______________________________________________________________________________________________

AliTPCGenDBConf::AliTPCGenDBConf():
   AliDCSGenDB()
{
}

//______________________________________________________________________________________________

AliTPCGenDBConf::AliTPCGenDBConf(const AliTPCGenDBConf& org):
  AliDCSGenDB(org)
{

//
//  Copy constructor
//

 ((AliTPCGenDBConf &) org).Copy(*this);
}

//______________________________________________________________________________________________
AliTPCGenDBConf::~AliTPCGenDBConf(){
//
// destructor
//

}
//______________________________________________________________________________________________
AliTPCGenDBConf& AliTPCGenDBConf::operator= (const AliTPCGenDBConf& org )
{
 //
 // assignment operator
 //
 if (&org == this) return *this;

 new (this) AliTPCGenDBConf(org);
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

