#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>

#include "AliTRDgeometry.h"

#include "AliTRDCalROC.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalFEE.h"

#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalPadStatus.h"
#include "AliTRDCalSingleChamberStatus.h"

#include "AliTRDCalPIDLQ.h"
#include "AliTRDCalMonitoring.h"

#endif

// Run numbers for the dummy file
const Int_t    gkDummyRunBeg = 0;
const Int_t    gkDummyRunEnd = 999999999;
AliCDBStorage *gStorLoc      = 0;

//_____________________________________________________________________________
TObject *CreatePadObject(const char *shortName, const char *description
                       , Float_t value)
{

  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det = 0; det < AliTRDgeometry::kNdet; ++det) {
    AliTRDCalROC *calROC = calPad->GetCalROC(det);
    for (Int_t channel = 0; channel < calROC->GetNchannels(); ++channel) {
      calROC->SetValue(channel, value);
    }
  }
  return calPad;

}

//_____________________________________________________________________________
TObject *CreateDetObject(const char *shortName, const char *description
                       , Float_t value)
{

  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det = 0; det < AliTRDgeometry::kNdet; ++det) {
    object->SetValue(det, value);
  }
  return object;

}

//_____________________________________________________________________________
TObject *CreateFEEObject() 
{

  AliTRDCalFEE *object = new AliTRDCalFEE("FEE"
                                         ,"TRD FEE calibration parameters");

  object->SetNumberOfTimeBins(24);

  object->SetTailCancelationTau1(0);
  object->SetTailCancelationTau2(0);
  object->SetTailCancelationAmp(0);

  object->SetPedestal(0);

  return object;

}

//_____________________________________________________________________________
TObject *CreatePRFWidthObject()
{

  AliTRDCalPad *calPad = new AliTRDCalPad("PRFWidth","PRFWidth");
  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::kNlayer; ++iLayer) {

    Float_t value = 0.0;
    switch (iLayer) {
      case 0: 
        value = 0.515; 
        break;
      case 1: 
        value = 0.502; 
        break;
      case 2: 
        value = 0.491; 
        break;
      case 3: 
        value = 0.481; 
        break;
      case 4: 
        value = 0.471; 
        break;
      case 5: 
        value = 0.463; 
        break;
      default: 
        cout << "CreatePRFWidthObject: UNEXPECTED" << endl; 
        return 0;
    }

    for (Int_t iStack = 0; iStack < AliTRDgeometry::kNstack; ++iStack) {
      for (Int_t iSector = 0; iSector < AliTRDgeometry::kNsector; ++iSector) {
        AliTRDCalROC *calROC = calPad->GetCalROC(iLayer,iStack,iSector);
        for (Int_t iChannel = 0; iChannel < calROC->GetNchannels(); ++iChannel) {
          calROC->SetValue(iChannel, value);
	}
      }
    }

  }
      
  return calPad;

}

//_____________________________________________________________________________
AliTRDCalChamberStatus *CreateChamberStatusObject()
{

  AliTRDCalChamberStatus *obj = new AliTRDCalChamberStatus("chamberstatus"
                                                          ,"chamberstatus");

  for (Int_t i = 0; i < AliTRDgeometry::kNdet; ++i) {
    obj->SetStatus(i, AliTRDCalChamberStatus::kInstalled);
  }

  return obj;

}

//_____________________________________________________________________________
AliTRDCalPadStatus *CreatePadStatusObject()
{
  AliTRDCalPadStatus     *obj = new AliTRDCalPadStatus("padstatus"
                                                      ,"padstatus");

  return obj;

}

//_____________________________________________________________________________
AliTRDCalMonitoring *CreateMonitoringObject()
{
  AliTRDCalMonitoring    *obj = new AliTRDCalMonitoring();

  return obj;

}

//_____________________________________________________________________________
TClonesArray* CreateRecoParamObject()
{
  TClonesArray *recos = new TClonesArray("AliTRDrecoParam", 3);
  recos->SetOwner(kTRUE);

  AliTRDrecoParam *rec = 0x0;
  AliTRDrecoParam *reco = new((*recos)[0]) AliTRDrecoParam(*(rec = AliTRDrecoParam::GetLowFluxParam()));
  delete rec;
  // further settings for low flux reco param
  // reco->SetThisAndThat()

  reco = new((*recos)[1]) AliTRDrecoParam(*(rec = AliTRDrecoParam::GetHighFluxParam()));
  delete rec;

  reco = new((*recos)[2]) AliTRDrecoParam(*(rec = AliTRDrecoParam::GetCosmicTestParam()));
  delete rec;

  return recos;

}

//_____________________________________________________________________________
AliCDBMetaData *CreateMetaObject(const char *objectClassName)
{

  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Christoph Blume");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-16-00"); //root version
  md1->SetComment("Ideal calibration values");
  
  return md1;

}

//_____________________________________________________________________________
void StoreObject(const char *cdbPath, TObject *object, AliCDBMetaData *metaData)
{

  AliCDBId id1(cdbPath,gkDummyRunBeg,gkDummyRunEnd); 
  gStorLoc->Put(object,id1,metaData); 

}
    

//_____________________________________________________________________________
void AliTRDCreateDummyCDB()
{

  cout << endl 
       << "TRD :: Creating dummy CDB for the runs " 
       << gkDummyRunBeg
       << " -- " 
       << gkDummyRunEnd
       << endl;
  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage("local://$ALICE_ROOT");
  if (!gStorLoc) {
    return;
  }

  TObject        *obj      = 0;
  AliCDBMetaData *metaData = 0;

  //
  // Reco Param Object
  //

  metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TClonesArray");
  metaData->SetResponsible("Alexandru Bercuci");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-19-04"); //root version
  metaData->SetComment("Ideal reconstruction parameters for low, high and cosmic runs");
  StoreObject("TRD/Calib/RecoParam",CreateRecoParamObject(),metaData);
  return;

  //
  // Pad objects
  //

  metaData = CreateMetaObject("AliTRDCalPad");  

  obj = CreatePadObject("LocalVdrift"       ,"TRD drift velocities (local variations)",1);
  StoreObject("TRD/Calib/LocalVdrift"       ,obj,metaData);
  
  obj = CreatePadObject("LocalT0"           ,"T0 (local variations)",0);
  StoreObject("TRD/Calib/LocalT0"           ,obj,metaData);

  obj = CreatePadObject("GainFactor"        ,"GainFactor (local variations)",1);
  StoreObject("TRD/Calib/LocalGainFactor"   ,obj,metaData);

  obj = CreatePRFWidthObject();
  StoreObject("TRD/Calib/PRFWidth"          ,obj,metaData);

  //
  // Detector objects
  //

  metaData = CreateMetaObject("AliTRDCalDet");
  
  obj = CreateDetObject("ChamberVdrift"     ,"TRD drift velocities (detector value)", 1.5);
  StoreObject("TRD/Calib/ChamberVdrift"     ,obj,metaData);
  
  obj = CreateDetObject("ChamberT0"         ,"T0 (detector value)",0);
  StoreObject("TRD/Calib/ChamberT0"         ,obj,metaData);
  
  obj = CreateDetObject("ChamberGainFactor" ,"GainFactor (detector value)", 1);
  StoreObject("TRD/Calib/ChamberGainFactor" ,obj,metaData);
  
  //
  // FEE object
  //

  metaData = CreateMetaObject("AliTRDCalFEE");
  obj = CreateFEEObject();
  StoreObject("TRD/Calib/FEE"               ,obj,metaData);
  
  //
  // Status objects
  //

  metaData = CreateMetaObject("AliTRDCalChamberStatus");
  obj = CreateChamberStatusObject();
  StoreObject("TRD/Calib/ChamberStatus"     ,obj,metaData);

  metaData = CreateMetaObject("AliTRDCalPadStatus");
  obj = CreatePadStatusObject();
  StoreObject("TRD/Calib/PadStatus"         ,obj,metaData);

  //
  // Monitoring object
  //

  metaData = CreateMetaObject("AliTRDCalMonitoring");
  obj = CreateMonitoringObject();
  StoreObject("TRD/Calib/MonitoringData"    ,obj,metaData);

}
