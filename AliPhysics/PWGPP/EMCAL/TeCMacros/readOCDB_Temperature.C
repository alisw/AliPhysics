#if !defined( __CINT__) || defined(__MAKECINT__)
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
#include <TGrid.h>
#include <TMap.h>
#include <TNtuple.h>
#include <stdio.h>
#include "TInfo.h"
TInfo *readOCDB_Temperature(Int_t runNb  = 286350, Bool_t debug=1);
void testOCDB_Temperature(Int_t runNb = 286350);
#endif

#ifndef _tempfuncs_
#define _tempfuncs_
TInfo *readOCDB_Temperature(Int_t runNb, Bool_t debug)
{
  if (!gGrid)
    TGrid::Connect("alien://");

  AliCDBManager*  cdb = AliCDBManager::Instance();
  if (cdb->GetDefaultStorage()==0)
    cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);

  AliCDBEntry *en=0;
  try {
    en=cdb->Get("EMCAL/Calib/Temperature/");
  } catch (...) { ; }
  if (!en) {
    printf("no entry found for run %d!\n", runNb);
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
  const Int_t kNumSens = 160;
  Int_t np[kNumSens] = {0};
  Double_t min[kNumSens] = {0};
  Double_t max[kNumSens] = {0};
  Double_t avg[kNumSens] = {0};
  Double_t rms[kNumSens] = {0};

  Double_t avTime = 0;
  UInt_t fTime = -1;
  UInt_t lTime = 0;
  for (Int_t isensor=0; isensor<kNumSens; isensor++) {
    AliEMCALSensorTemp *o = arr->GetSensor(isensor);
    if (!o)
      continue;

    UInt_t startt = o->GetStartTime();
    UInt_t stopt  = o->GetEndTime();

    if (debug)
      cout  << "Sensor " << isensor
            << " " << o->GetStringID() << ":"
            << " side " << o->GetSide()
            << " sector " << o->GetSector()
            << " num " << o->GetNum()
            << " startTime " << startt
            << " endTime " << stopt
            << endl;

    np[isensor]  = 0;
    min[isensor] = +100;
    max[isensor] = -100;
    avg[isensor] = 0;
    rms[isensor] = 0;

    AliSplineFit *f = o->GetFit();
    if (f) {
      np[isensor] = f->GetKnots();
      if (debug)
        cout << " np " << np[isensor] << endl;
      Double_t *x = f->GetX();
      Double_t *y0 = f->GetY0();
      Double_t *y1 = f->GetY1();
      for (Int_t i=0; i<np[isensor]; ++i) {
        if (debug)
          cout  << " i " << i
                << " x " << x[i]
                << " y0 " << y0[i]
                << " y1 " << y1[i]
                << endl;
	avg[isensor] += y0[i];
	rms[isensor] += y0[i]*y0[i];
        if (min[isensor]>y0[i]) min[isensor]=y0[i];
        if (max[isensor]<y0[i]) max[isensor]=y0[i];
      }
    } else {
      TGraph *g = o->GetGraph();
      if (g) {
        np[isensor] = g->GetN();
        if (debug)
          cout << " np " << np[isensor] << endl;
        Double_t *x = g->GetX();
        Double_t *y0 = g->GetY();
        for (Int_t i=0; i<np[isensor]; i++) {
          if (debug)
            cout  << " i " << i
                  << " x " << x[i]
                  << " y0 " << y0[i]
                  << endl;
	  avg[isensor] += y0[i];
	  rms[isensor] += y0[i]*y0[i];
          if (min[isensor]>y0[i]) min[isensor]=y0[i];
          if (max[isensor]<y0[i]) max[isensor]=y0[i];
        }
      }
    }

    if (np[isensor]>0) {
      avg[isensor] /= np[isensor];
      if (np[isensor]>1) {
	rms[isensor] /= np[isensor];
	rms[isensor] = TMath::Sqrt(rms[isensor]-avg[isensor]*avg[isensor])/(np[isensor]-1);
      } else {
	rms[isensor] = 0;
      }
      if (debug)
        cout << "Avg Temp: " << avg[isensor] << " RMS Temp: " << rms[isensor] << " Min Temp: " << min[isensor] << " Max Temp: " << max[isensor] << endl;
      info->Set(isensor, avg[isensor], rms[isensor], min[isensor], max[isensor]);
      avTime += startt + stopt;
      if (startt<fTime)
        fTime=startt;
      if (stopt>lTime)
        lTime=stopt;
    }

    if (0)
      cout  << "Tree: " << runNb
            << " " << isensor
            << " " << o->GetSide()
            << " " << o->GetSector()
            << " " << o->GetNum()
            << " " << o->GetStartTime()
            << " " << o->GetEndTime()
            << " " << avg[isensor]
            << " " << rms[isensor]
            << " " << max[isensor]
            << " " << max[isensor] << endl;
  }
  Int_t nv = info->Nvalid();
  if (nv>0) {
    avTime /= 2*nv;
    info->SetTime(avTime,fTime,lTime);
  }
  return info;
}

void testOCDB_Temperature(Int_t runNb)
{
  TInfo *i = readOCDB_Temperature(runNb,1);
  i->Print();
  return;
}
#endif
