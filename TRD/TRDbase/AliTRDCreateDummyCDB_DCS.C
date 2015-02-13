#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>

#include "AliTRDgeometry.h"

#include "AliTRDCalDCS.h"
#include "AliTRDCalDCSFEE.h"

#endif

// ===================================================
// Modified version of the macro AliTRDCreateDummyCDB 
// to create the default object for the DCS CDB entry
// Modifications done by Frederick Kramer 2010-05-11
// ===================================================



// Run numbers for the dummy file
const Int_t    gkDummyRunBeg = 0;
const Int_t    gkDummyRunEnd = 999999999;
AliCDBStorage *gStorLoc      = 0;



//_____________________________________________________________________________
TObjArray *CreateDCSObject()
{
  const Int_t nROC = 540;
  const Int_t nROB = 8;
  const Int_t nMCM = 18;


  AliTRDCalDCS* fCalDCSObjSOR = new AliTRDCalDCS();
  AliTRDCalDCS* fCalDCSObjEOR = new AliTRDCalDCS();
  TObjArray*    fFEEArrSOR    = new TObjArray(nROC);
  TObjArray*    fFEEArrEOR    = new TObjArray(nROC);
  fFEEArrSOR->SetOwner();
  fFEEArrEOR->SetOwner();
  fCalDCSObjSOR->SetObjectStat(1);
  fCalDCSObjEOR->SetObjectStat(1);
  fFEEArrSOR->SetObjectStat(1);
  fFEEArrEOR->SetObjectStat(1);


  for (Int_t iROC=0; iROC<nROC; iROC++) {
    AliTRDCalDCSFEE* fDCSFEEObjSOR = new AliTRDCalDCSFEE();
    fDCSFEEObjSOR->SetObjectStat(1);
    AliTRDgeometry aliGeo;
    Int_t sm    = aliGeo.GetSector(iROC);
    Int_t stack = aliGeo.GetStack(iROC);
    Int_t layer = aliGeo.GetLayer(iROC);

    fDCSFEEObjSOR->SetConfigName("cf_p_zs-s16-deh_tb30_csmtrk_ptrg");
    fDCSFEEObjSOR->SetConfigTag(1681);
    fDCSFEEObjSOR->SetConfigVersion("r3303");
    fDCSFEEObjSOR->SetNumberOfTimeBins(30);
    fDCSFEEObjSOR->SetSM(sm);
    fDCSFEEObjSOR->SetStack(stack);
    fDCSFEEObjSOR->SetLayer(layer);
    fDCSFEEObjSOR->SetFilterType("p");
    fDCSFEEObjSOR->SetReadoutParam("zs");
    fDCSFEEObjSOR->SetTrackletMode("csmtrk");
    fDCSFEEObjSOR->SetTriggerSetup("ptrg");

    AliTRDCalDCSFEE* fDCSFEEObjEOR = (AliTRDCalDCSFEE*)fDCSFEEObjSOR->Clone();

    for (Int_t iROB=0; iROB<nROB; iROB++) {
      for (Int_t iMCM=0; iMCM<nMCM; iMCM++) {
	// SOR
	fDCSFEEObjSOR->SetMCMGlobalState(iROB, iMCM, 3);
	fDCSFEEObjSOR->SetMCMStateNI(iROB, iMCM, 0);
	fDCSFEEObjSOR->SetMCMEventCnt(iROB, iMCM, 0);
	fDCSFEEObjSOR->SetMCMPtCnt(iROB, iMCM, 0);
	// EOR
	fDCSFEEObjEOR->SetMCMGlobalState(iROB, iMCM, 3);
	fDCSFEEObjEOR->SetMCMStateNI(iROB, iMCM, 0);
	fDCSFEEObjEOR->SetMCMEventCnt(iROB, iMCM, 1000);
	fDCSFEEObjEOR->SetMCMPtCnt(iROB, iMCM, 1000);
      } //iMCM
    } // iROB

    fFEEArrSOR->AddAt(fDCSFEEObjSOR, iROC);
    fFEEArrEOR->AddAt(fDCSFEEObjEOR, iROC);
  } // iROC


  fCalDCSObjSOR->SetFEEArr(fFEEArrSOR);
  fCalDCSObjSOR->EvaluateGlobalParameters();
  fCalDCSObjEOR->SetFEEArr(fFEEArrEOR);
  fCalDCSObjEOR->EvaluateGlobalParameters();

  TObjArray* fCalObjArray = new TObjArray(2);
  fCalObjArray->SetOwner();
  fCalObjArray->AddAt(fCalDCSObjSOR,0);
  fCalObjArray->AddAt(fCalDCSObjEOR,1);


  return fCalObjArray;
}


//_____________________________________________________________________________
AliCDBMetaData *CreateMetaObject(const char *objectClassName)
{

  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Frederick Kramer");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-26-00b"); //root version
  md1->SetComment("Ideal DCS configuration data in two AliTRDCalDCS objects in one TObjArray (0:SOR, 1:EOR).");
  
  return md1;

}

//_____________________________________________________________________________
void StoreObject(const char *cdbPath, TObjArray *object, AliCDBMetaData *metaData)
{

  AliCDBId id1(cdbPath,gkDummyRunBeg,gkDummyRunEnd); 
  gStorLoc->Put(object,id1,metaData); 

}
    

//_____________________________________________________________________________
void AliTRDCreateDummyCDB_DCS()
{

  cout << endl 
       << "TRD :: Creating dummy CDB for the runs " 
       << gkDummyRunBeg
       << " -- " 
       << gkDummyRunEnd
       << endl;
  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) {
    return;
  }

  TObjArray      *obj      = 0;
  AliCDBMetaData *metaData = 0;

  //
  // Status objects
  //

  metaData = CreateMetaObject("DCS");
  obj = CreateDCSObject();
  StoreObject("TRD/Calib/DCS", obj, metaData);
  
}

