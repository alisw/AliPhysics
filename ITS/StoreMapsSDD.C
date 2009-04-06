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
  

  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetObjectClassName("TObjArray");
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("Simulated data");



  AliCDBId mapT("ITS/Calib/MapsTimeSDD",firstRun,lastRun);
  TObjArray tmap(520);
  tmap.SetOwner(kFALSE);

  TRandom3 *gran = new TRandom3();
  
  AliITSMapSDD* mapTime0;
  AliITSMapSDD* mapTime1;
  for(Int_t mod=0;mod<260;mod++){
    // maps
    Char_t name[20];
    sprintf(name,"DriftTimeMap_%d_%d\n",mod,0);
    Int_t nbinsan=1;
    if(mod==10 || mod==240){
      nbinsan=256;
      sprintf(name,"DriftTimeMap_%d_%d\n",mod,0);
      mapTime0 = new AliITSMap2DSDD(name,nbinsan,72);
      sprintf(name,"DriftTimeMap_%d_%d\n",mod,1);
      mapTime1 = new AliITSMap2DSDD(name,nbinsan,72);
    }else{
      sprintf(name,"DriftTimeMap_%d_%d\n",mod,0);
      mapTime0 = new AliITSMap1DSDD(name,72);
      sprintf(name,"DriftTimeMap_%d_%d\n",mod,1);
      mapTime1 = new AliITSMap1DSDD(name,72);
    }
    for(Int_t nan = 0;nan< nbinsan;nan++){
      for(Int_t nt = 0;nt<36*2;nt++){
	mapTime0->SetCellContent(nan,nt,gran->Gaus(0,20));
	mapTime1->SetCellContent(nan,nt,gran->Gaus(0,20));		     
      }
    }
    tmap.Add(mapTime0);
    tmap.Add(mapTime1); 
    printf("Added module %d\n",mod);
  }
    
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&tmap, mapT, md);

}
