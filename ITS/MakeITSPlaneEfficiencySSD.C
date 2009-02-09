#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRandom3.h>
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBRunRange.h"
#include "AliCDBId.h"
#include "AliITSPlaneEffSSD.h"
#endif

void MakeITSPlaneEfficiencySSD(Int_t firstRun=0,Int_t lastRun=AliCDBRunRange::Infinity()){
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSPlaneEff");
  md1->SetResponsible("Giuseppe Bruno");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 02/01/08"); //root version

  AliCDBId idplaneeffSSD("ITS/PlaneEff/PlaneEffSSD",firstRun, lastRun);
  
  AliITSPlaneEffSSD* planeeffSSD = new AliITSPlaneEffSSD();
  TRandom3 *gran = new TRandom3();

//  planeeffSSD->SetOwner(kFALSE);

  // loop over SSD modules
  Bool_t BFound=kFALSE;
  for(UInt_t mod=0;mod<1698;mod++){
  // suppose to have 1000 tracks in each module and an average efficiency of 99%
    for(Int_t j=0; j<1000; j++) {
      BFound=kFALSE;
      //if (gRandom->Uniform(0,1000)>10) BFound=kTRUE;
      if (1000*gran->Uniform()>10) BFound=kTRUE;
      planeeffSSD->UpDatePlaneEff(BFound,mod);
    }
  }
  if(AliCDBManager::Instance()->GetDefaultStorage()->Put(planeeffSSD, idplaneeffSSD, md1))
  printf("Local CDB file with random SSD plane efficiencies written \n");
 delete gran;
 delete planeeffSSD;
 delete md1;
}
