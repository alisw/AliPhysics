#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRandom3.h>
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBRunRange.h"
#include "AliCDBId.h"
#include "AliITSPlaneEffSPD.h"
#endif

void MakeITSPlaneEfficiencySPD(Int_t firstRun=0, Int_t lastRun=AliCDBRunRange::Infinity(), 
 Double_t eff=0.99, Int_t nTried=1000){
  
  if(eff<0 || eff > 1) {
   printf("Efficiency must be in the range [0,1]: nothing done");
  }
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSPlaneEff");
  md1->SetResponsible("Giuseppe Bruno");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 16/02/08"); //root version

  AliCDBId idplaneeffSPD("ITS/PlaneEff/PlaneEffSPD",firstRun, lastRun);
  
  AliITSPlaneEffSPD* planeeffSPD = new AliITSPlaneEffSPD();
  TRandom3 *gran = new TRandom3();

//  planeeffSPD->SetOwner(kFALSE);
  
//  Int_t nTried=1000;
  Double_t limit=nTried;
  limit*=(1-eff);
  printf("limit = %f",limit);
  // loop over SPD chip
  Bool_t BFound=kFALSE;
  for(UInt_t key=0;key<planeeffSPD->Nblock();key++){
  //for(UInt_t mod=0;mod<240;mod++){
  //for(UInt_t chip=0;chip<5;chip++){
  // suppose to have 1000 tracks in each chip and an average efficiency of 99%
    for(Int_t j=0; j<nTried; j++) {
      BFound=kFALSE;
      //if (gRandom->Uniform(0,1000)>10) BFound=kTRUE;
      if (nTried*gran->Uniform()>limit) BFound=kTRUE;
      //planeeffSPD->UpDatePlaneEff(BFound,mod,chip);
      planeeffSPD->UpDatePlaneEff(BFound,key);
    }
  //}}
  }
  if(AliCDBManager::Instance()->GetDefaultStorage()->Put(planeeffSPD, idplaneeffSPD, md1))
  printf("Local CDB file with random SPD plane efficiencies written \n");
 delete gran;
 delete planeeffSPD;
 delete md1;
}
