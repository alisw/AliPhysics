#include <AliCDBEntry.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliDCSValue.h>
#include <AliEMCALPreprocessor.h>
#include <AliEMCALSensorTempArray.h>
#include <AliEMCALTriggerDCSConfig.h>
#include <AliEMCALTriggerSTUDCSConfig.h>
#include <AliEMCALTriggerTRUDCSConfig.h>
#include <AliLog.h>
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TMap.h>
#include <TNtuple.h>
#include <stdio.h>
#include "TInfo.h"

TInfo *readOCDB_Temperature(Int_t runNb  = 286350, Bool_t debug=1)
{
  TGrid::Connect("alien://");

  AliCDBManager*  cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);

  AliCDBEntry *en=cdb->Get("EMCAL/Calib/Temperature/");
  if (!en) {
    printf("no entry found!\n");
    return 0;
  }

  TInfo *info = new TInfo(runNb);

  AliEMCALSensorTempArray *arr = dynamic_cast< AliEMCALSensorTempArray *> (en->GetObject());
  //arr->Print();
  cout << " NumSensors " << arr->NumSensors()
       << " GetFirstIdDCS() " << arr->GetFirstIdDCS()
       << " GetLastIdDCS() " << arr->GetLastIdDCS()
       << endl;

  // info for each sensor
  const int kNumSens = 160;
  int np[kNumSens] = {0};
  double min[kNumSens] = {0};
  double max[kNumSens] = {0};

  Double_t avTime = 0;
  UInt_t fTime = -1;
  UInt_t lTime = 0;

  for (int isensor=0; isensor<kNumSens; isensor++) {
    AliEMCALSensorTemp *o = arr->GetSensor(isensor);
    if (!o) 
      continue;

    UInt_t startt = o->GetStartTime();
    UInt_t stopt = o->GetEndTime();

    if (debug)
      cout << "Sensor " << isensor 
	   << " " << o->GetStringID() << ":" 
	   << " side " << o->GetSide()
	   << " sector " << o->GetSector()
	   << " num " << o->GetNum()
	   << " startTime " << startt
	   << " endTime " << stopt
	   << endl;

    AliSplineFit *f = o->GetFit();
    np[isensor] = 0;
    min[isensor] = +100;
    max[isensor] = -100;

    if (f) {
      np[isensor] = f->GetKnots();
      if (debug)
	cout << " np " << np[isensor] << endl;
      Double_t *x = f->GetX();
      Double_t *y0 = f->GetY0();
      Double_t *y1 = f->GetY1();
      for (int i=0; i<np[isensor]; i++) {
	if (debug)
	  cout << " i " << i
	       << " x " << x[i]
	       << " y0 " << y0[i]
	       << " y1 " << y1[i]
	       << endl;
	
	if (min[isensor]>y0[i]) min[isensor]=y0[i];
	if (max[isensor]<y0[i]) max[isensor]=y0[i];
      }
    }
    if (np[isensor]>0) {
      if (debug)
	cout << "Min Temp: " << min[isensor] << " Max Temp: " << max[isensor] << endl;
      info->Set(isensor, min[isensor], max[isensor]);
      avTime += startt + stopt;
      if (startt<fTime)
	fTime=startt;
      if (stopt>lTime)
	lTime=stopt;
    }

    if (0)
      cout << "Tree: " << runNb 
	   << " " << isensor
	   << " " << o->GetSide() 
	   << " " << o->GetSector()
	   << " " << o->GetNum()
	   << " " << o->GetStartTime()
	   << " " << o->GetEndTime()
           << " " << min[isensor] 
	   << " " << max[isensor] << endl;
  }
  Int_t nv = info->Nvalid();
  if (nv>0) {
    avTime /= 2*nv;
    info->SetTime(avTime,fTime,lTime);
  }
  return info;
}

void testOCDB_Temperature(Int_t runNb  = 286350)
{
  TInfo *i = readOCDB_Temperature(runNb,1);
  i->Print();
  return;
}

void read_LHC18d() 
{
  Int_t runs[] = {285978,285979,285980,286014,286018,286025,286026,286027,286030,286064,286124,286127,286129,286130,286154,286157,286159,286198,286201,286202,286203,286229,286230,286231,286254,286255,286256,286257,286258,286261,286263,286282,286284,286287,286288,286289,286308,286309,286310,286311,286312,286313,286314,286336,286337,286340,286341,286345,286348,286349,286350};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=0;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << i << " " << rn << endl;
    TInfo *ti = readOCDB_Temperature(rn,0);
    ti->Print();
    arr.Add(ti);
  }

  TFile *outf = TFile::Open("temperatures.root","update");
  arr.Write("temperatures_lhc18d",TObject::kSingleKey);
  outf->ls();
  outf->Close();
}
