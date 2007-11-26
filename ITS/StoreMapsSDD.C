#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSMapSDD.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include <TObjArray.h>
#include <TRandom3.h>
#endif

void StoreMapsSDD(Int_t firstRun=0,Int_t lastRun=9999999 ){
  ///////////////////////////////////////////////////////////////////////
  // Macro to generate and store the residual maps for SDD             //
  // Generates:                                                        //
  //  1 file with 520 AliITSMapSDD anode maps (MapsAnodeSDD)           //
  //  1 file with 520 AliITSMapSDD drift coordinate maps (MapsTimeSDD) //
  ///////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
  }
  

  AliCDBMetaData *md3 = new AliCDBMetaData();
  md3->SetObjectClassName("AliITSMapSDD");
  md3->SetResponsible("Elisabetta Crescio, Francesco Prino");
  md3->SetBeamPeriod(0);
  md3->SetAliRootVersion("Head 24 sept. 2007");
  md3->SetComment("This is a test");

  AliCDBMetaData *md4 = new AliCDBMetaData();
  md4->SetObjectClassName("AliITSMapSDD");
  md4->SetResponsible("Elisabetta Crescio, Francesco Prino");
  md4->SetBeamPeriod(0);
  md4->SetAliRootVersion("Head 24 sept. 2007");
  md4->SetComment("This is a test");


  AliCDBId mapA("ITS/Calib/MapsAnodeSDD",firstRun,lastRun);
  TObjArray anmap(520);
  anmap.SetOwner(kFALSE);

  AliCDBId mapT("ITS/Calib/MapsTimeSDD",firstRun,lastRun);
  TObjArray tmap(520);
  tmap.SetOwner(kFALSE);

  TRandom3 *gran = new TRandom3();
  
  for(Int_t mod=0;mod<260;mod++){
    // maps
    Char_t name[20];
    sprintf(name,"AnodeMap_%d_%d\n",mod,0);
    AliITSMapSDD* mapAnodes0 = new AliITSMapSDD(name);
    sprintf(name,"DriftTimeMap_%d_%d\n",mod,0);
    AliITSMapSDD* mapTime0 = new AliITSMapSDD(name);
    sprintf(name,"AnodeMap_%d_%d\n",mod,1);
    AliITSMapSDD* mapAnodes1 = new AliITSMapSDD(name);
    sprintf(name,"DriftTimeMap_%d_%d\n",mod,1);
    AliITSMapSDD* mapTime1 = new AliITSMapSDD(name);

    for(Int_t nan = 0;nan<256;nan++){
      for(Int_t nt = 0;nt<36*2;nt++){
	mapAnodes0->SetCellContent(nan,nt,gran->Gaus(0,20));
	mapTime0->SetCellContent(nan,nt,gran->Gaus(0,20));
	mapAnodes1->SetCellContent(nan,nt,gran->Gaus(0,20));
	mapTime1->SetCellContent(nan,nt,gran->Gaus(0,20));		     
      }
    }
    anmap.Add(mapAnodes0);
    tmap.Add(mapTime0);
    anmap.Add(mapAnodes1);
    tmap.Add(mapTime1); 
    printf("Added module %d\n",mod);
  }
    
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&anmap, mapA, md3);
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&tmap, mapT, md4);

}
