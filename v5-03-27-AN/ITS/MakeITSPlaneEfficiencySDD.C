#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRandom3.h>
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBRunRange.h"
#include "AliCDBId.h"
#include "AliITSPlaneEffSDD.h"
#endif
void MakeITSPlaneEfficiencySDD(Int_t firstRun=0,Int_t lastRun=AliCDBRunRange::Infinity()){
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSPlaneEff");
  md1->SetResponsible("Giuseppe Bruno");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 16/01/08"); //root version

  AliCDBId idplaneeffSDD("ITS/PlaneEff/PlaneEffSDD",firstRun, lastRun);
  
  AliITSPlaneEffSDD* planeeffSDD = new AliITSPlaneEffSDD();
  TRandom3 *gran = new TRandom3();

//  planeeffSDD->SetOwner(kFALSE);
  

  // loop over SDD basic block
  Bool_t BFound=kFALSE;
  //for(Int_t key=0;key<2080;key++){
  for(UInt_t mod=0;mod<260;mod++){
  for(UInt_t chip=0;chip<4;chip++){
  for(UInt_t wing=0;wing<2;wing++){
  for(UInt_t subw=0;subw<1;subw++){
  // suppose to have 1000 tracks in each block and an average efficiency of 99%
    for(Int_t j=0; j<1000; j++) {
      BFound=kFALSE;
      //if (gRandom->Uniform(0,1000)>10) BFound=kTRUE;
      if (1000*gran->Uniform()>10) BFound=kTRUE;
      planeeffSDD->UpDatePlaneEff(BFound,mod,chip,wing,subw);
    }
  }}}}
  if(AliCDBManager::Instance()->GetDefaultStorage()->Put(planeeffSDD, idplaneeffSDD, md1))
  printf("Local CDB file with random SDD plane efficiencies written \n");
  delete gran;
  delete planeeffSDD;
  delete md1;
}
