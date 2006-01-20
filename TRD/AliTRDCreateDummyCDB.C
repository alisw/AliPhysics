#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>

#include "AliTRDgeometry.h"

#include "AliTRDCalROC.h"
#include "AliTRDCalChamber.h"
#include "AliTRDCalStack.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalGlobals.h"

#endif

// run number for the dummy file
const Int_t gkDummyRun = 0;
AliCDBStorage* gStorLoc = 0;

TObject* CreatePadObject(const char* shortName, const char* description, Float_t value)
{
  AliTRDCalPad *calPad = new AliTRDCalPad(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
  {
    AliTRDCalROC *calROC = calPad->GetCalROC(det);
    for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel)
      calROC->SetValue(channel, value);
  }
  return calPad;
}

TObject* CreateDetObject(const char* shortName, const char* description, Float_t value)
{
  AliTRDCalDet *object = new AliTRDCalDet(shortName, description);
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
    object->SetValue(det, value);
  return object;
}

TObject* CreateGlobalsObject()
{
  AliTRDCalGlobals *object = new AliTRDCalGlobals("Globals", "Global TRD calibration parameters");
  
  object->SetSamplingFrequency(10.0);
  object->SetNumberOfTimeBins(22);
  
  return object;
}

TObject* CreateChamberObject()
{
  AliTRDCalChamber *object = new AliTRDCalChamber("Chamber", "TRD chamber positions");
  
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
  {
    object->SetPos(det, 0, 0, 0);
    object->SetRot(det, 0, 0, 0);
  }
  
  return object;
}

TObject* CreateStackObject()
{
  AliTRDCalStack *object = new AliTRDCalStack("Stack", "TRD stack positions");
  
  for (Int_t sect=0; sect<AliTRDgeometry::kNsect; ++sect)
  {
    for (Int_t chamber=0; chamber<AliTRDgeometry::kNcham; ++chamber)
    {
      object->SetPos(chamber, sect, 0, 0, 0);
      object->SetRot(chamber, sect, 0, 0, 0);
    }
  }
    
  return object;
}

TObject* CreatePRFWidthObject()
{
  AliTRDCalPad *calPad = new AliTRDCalPad("PRFWidth","PRFWidth");
  for (Int_t plane=0; plane<AliTRDgeometry::kNplan; ++plane)
  {
    Float_t value = 0;
    switch (plane)
    {
      case 0: value = 0.515; break;
      case 1: value = 0.502; break;
      case 2: value = 0.491; break;
      case 3: value = 0.481; break;
      case 4: value = 0.471; break;
      case 5: value = 0.463; break;
      default: cout << "CreatePRFWidthObject: UNEXPECTED" << endl; return 0;
    }
    for (Int_t chamber=0; chamber<AliTRDgeometry::kNcham; ++chamber)
    {
      for (Int_t sector=0; sector<AliTRDgeometry::kNsect; ++sector)
      {
        AliTRDCalROC *calROC = calPad->GetCalROC(plane, chamber, sector);
        for (Int_t channel=0; channel<calROC->GetNchannels(); ++channel)
          calROC->SetValue(channel, value);
      }
    }
  }
      
  return calPad;
}

AliTRDCalPIDLQ* CreatePIDLQObject()
{
  AliTRDCalPIDLQ* pid = new AliTRDCalPIDLQ();
  pid->ReadData("$ALICE_ROOT/TRD/TRDdEdxHistogramsV1.root");
  pid->SetMeanChargeRatio(1.0); // The factor is the ratio of Mean of pi charge dist.
                    // for the New TRD code divided by the Mean of pi charge
                    // dist. given in AliTRDCalPIDLQ object
  
  return pid;
}

AliCDBMetaData* CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Jan Fiete Grosse-Oetringhaus");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-06-00"); //root version
  md1->SetComment("The dummy values in this calibration file are for testing only");
  
  return md1;
}

void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{
  AliCDBId id1(cdbPath, gkDummyRun, gkDummyRun); 
  gStorLoc->Put(object, id1, metaData); 
}
    

void AliTRDCreateDummyCDB()
{
  cout << endl << "TRD :: Creating dummy CDB with event number " << gkDummyRun << endl;
  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage("local://$ALICE_ROOT");
  if (!gStorLoc)
    return;

  TObject* obj = 0;
  AliCDBMetaData* metaData = 0;
  
  metaData = CreateMetaObject("AliTRDCalPad");
  
  obj = CreatePadObject("LocalVdrift","TRD drift velocities (local variations)", 1);
  StoreObject("TRD/Calib/LocalVdrift", obj, metaData);
  
  obj = CreatePadObject("LocalT0","T0 (local variations)", 1);
  StoreObject("TRD/Calib/LocalT0", obj, metaData);
  
  obj = CreatePadObject("GainFactor","GainFactor", 1);
  StoreObject("TRD/Calib/GainFactor", obj, metaData);

  obj = CreatePRFWidthObject();
  StoreObject("TRD/Calib/PRFWidth", obj, metaData);

  metaData = CreateMetaObject("AliTRDCalDet");
  
  obj = CreateDetObject("ChamberVdrift","TRD drift velocities (detector value)", 1.5);
  StoreObject("TRD/Calib/ChamberVdrift", obj, metaData);
  
  obj = CreateDetObject("ChamberT0","T0 (detector value)", 0);
  StoreObject("TRD/Calib/ChamberT0", obj, metaData);
  
  metaData = CreateMetaObject("AliTRDCalGlobals");
  obj = CreateGlobalsObject();
  StoreObject("TRD/Calib/Globals", obj, metaData);
  
  metaData = CreateMetaObject("AliTRDCalChamber");
  obj = CreateChamberObject();
  StoreObject("TRD/Calib/Chamber", obj, metaData);
  
  metaData = CreateMetaObject("AliTRDCalStack");
  obj = CreateStackObject();
  StoreObject("TRD/Calib/Stack", obj, metaData);
  
  metaData = CreateMetaObject("AliTRDCalPIDLQ");
  obj = CreatePIDLQObject();
  StoreObject("TRD/Calib/PIDLQ", obj, metaData);
}
