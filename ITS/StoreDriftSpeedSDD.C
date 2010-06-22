#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include <TObjArray.h>
#include <TRandom3.h>
#endif

void StoreDriftSpeedSDD(Int_t firstRun=0,Int_t lastRun=AliCDBRunRange::Infinity()){
  ///////////////////////////////////////////////////////////////////////
  // Macro to generate and store the drift speed files for SDD         //
  // Generates:                                                        //
  //  1 file with 520 AliITSDriftSpeedArraySDD objects with            //
  ///////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
  }
  

  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("TObjArray");
  md1->SetResponsible("Francesco Prino");
  md1->SetBeamPeriod(0);
  md1->SetComment("Simulated data");

  AliCDBId drSpeed("ITS/Calib/DriftSpeedSDD",firstRun, lastRun);
  TObjArray vdrift(520);
  vdrift.SetOwner(kFALSE);



  Double_t drVelParam[4]={6.533,0.001285,-0.000005,0};
  Double_t edrVelParam[4]={0.1,0,0,0};
  Double_t drVel[4];
  TRandom3 *gran = new TRandom3();
  
  for(Int_t mod=0;mod<260;mod++){
    AliITSDriftSpeedArraySDD *arr0 = new AliITSDriftSpeedArraySDD(5);
    AliITSDriftSpeedArraySDD *arr1 = new AliITSDriftSpeedArraySDD(5);
    for(Int_t iev=0; iev<5;iev++){
      for(Int_t ic=0;ic<4;ic++) drVel[ic]=gran->Gaus(drVelParam[ic],edrVelParam[ic]);
      AliITSDriftSpeedSDD *v0=new AliITSDriftSpeedSDD(iev*20,iev+1000,3,drVel);
      arr0->AddDriftSpeed(v0);
      for(Int_t ic=0;ic<4;ic++) drVel[ic]=gran->Gaus(drVelParam[ic],edrVelParam[ic]);
      AliITSDriftSpeedSDD *v1=new AliITSDriftSpeedSDD(iev*20,iev+1000,3,drVel);
      arr1->AddDriftSpeed(v1);
    }
    vdrift.Add(arr0);
    vdrift.Add(arr1);
    printf("Added module %d\n",mod);
  }
    
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&vdrift, drSpeed, md1);   
}
