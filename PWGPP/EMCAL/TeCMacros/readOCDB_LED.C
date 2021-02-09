#if !defined( __CINT__) || defined(__MAKECINT__)
#include <AliCDBEntry.h>
#include <AliCDBEntry.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliCaloCalibSignal.h>
#include <AliDCSValue.h>
#include <AliEMCALPreprocessor.h>
#include <AliEMCALSensorTempArray.h>
#include <AliEMCALTriggerDCSConfig.h>
#include <AliEMCALTriggerSTUDCSConfig.h>
#include <AliEMCALTriggerTRUDCSConfig.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMap.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include "LInfo.h"

LInfo *readOCDB_LED(Int_t runNb  = 286350, Bool_t debug=1);
void testOCDB_LED(Int_t runNb  = 286350);
#endif

#ifndef _ledfuncs_
#define _ledfuncs_
LInfo *readOCDB_LED(Int_t runNb, Bool_t debug)
{
  if (!gGrid)
    TGrid::Connect("alien://");

  AliCDBManager*  cdb = AliCDBManager::Instance();
  if (cdb->GetDefaultStorage()==0)
    cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);

  AliCDBEntry *en = 0;
  try {
    en=cdb->Get("EMCAL/Calib/LED/");
  } catch (...) { ; }
  if (!en) {
    printf("no entry found for run %d!\n", runNb);
    return 0;
  }

  AliCaloCalibSignal *emcCalibSignal=dynamic_cast<AliCaloCalibSignal*>(en->GetObject());
  if (!emcCalibSignal) {
    printf("can not cast to calib signal!\n");
    return 0;
  }
  //emcCalibSignal->Print();

  LInfo *ret = new LInfo(runNb);

  TTree *treeLed = emcCalibSignal->GetTreeAvgAmpVsTime();
  Int_t tea = treeLed->GetEntries();
  TTree *treeMon = emcCalibSignal->GetTreeLEDAvgAmpVsTime();
  Int_t tem = treeMon->GetEntries();
  if (debug) {
    cout << " LED tree entries " << tea << endl;
    treeLed->Print();
    cout << " LEDMon entries " << tem << endl;
    treeMon->Print();
  }

  Int_t mod   = 0;
  Int_t col   = 0;
  Int_t row   = 0;
  Int_t gain  = 0;
  Int_t strip = 0;

  // associate variables for LEDMON
  Int_t    fRefMon  = 0;
  Double_t fAmpMon  = 0;
  Double_t fRMSMon  = 0;
  Double_t fHourMon = 0;
  treeMon->SetBranchAddress("fRefNum",&fRefMon);
  treeMon->SetBranchAddress("fAvgAmp",&fAmpMon);
  treeMon->SetBranchAddress("fRMS",&fRMSMon);
  treeMon->SetBranchAddress("fHour",&fHourMon);
  for (Int_t i=0; i<tem; ++i) {
    treeMon->GetEntry(i);
    emcCalibSignal->DecodeRefNum(fRefMon, &mod, &strip, &gain);
    if (debug) 
      Printf("fRefNum=%d, mod= %d, strip=%d, amp=%.3f, rms=%.3f, gain=%d", fRefMon, mod, strip, fAmpMon, fRMSMon, gain);
    ret->FillStrip(mod, gain, strip, fAmpMon, fRMSMon);
  }

  // associate variables for LED
  Int_t    fChanLed  = 0;
  Double_t fAmpLed  = 0;
  Double_t fRMSLed  = 0;
  Double_t fHourLed = 0;
  treeLed->SetBranchAddress("fChannelNum",&fChanLed);
  treeLed->SetBranchAddress("fHour",&fHourLed);
  treeLed->SetBranchAddress("fAvgAmp",&fAmpLed);
  treeLed->SetBranchAddress("fRMS",&fRMSLed);
  for (Int_t i=0; i<tea; ++i) {
    treeLed->GetEntry(i);
    emcCalibSignal->DecodeChannelNum(fChanLed, &mod, &col, &row, &gain);
    if (debug)
      Printf("fChannelNum=%d, mod=%d, col=%d, row=%d, amp=%.3f, rms=%.3f, gain=%d", fChanLed, mod, col, row, fAmpLed, fRMSLed, gain);
    ret->FillLed(mod,gain,col, row, fAmpLed, fRMSLed);
  }

  return ret;
}

void testOCDB_LED(Int_t runNb)
{
  LInfo *i = readOCDB_LED(runNb,1);
  i->Print();
  return;
}
#endif
