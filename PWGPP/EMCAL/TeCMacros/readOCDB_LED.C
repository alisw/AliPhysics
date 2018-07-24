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

  TTree *treeAmp = emcCalibSignal->GetTreeAvgAmpVsTime();
  TTree *treeLEDAmp = emcCalibSignal->GetTreeLEDAvgAmpVsTime();
  Int_t tea = treeAmp->GetEntries();
  Int_t tem = treeLEDAmp->GetEntries();
  if (debug) {
    cout << " Tree entries " << tea << endl;
    cout << " TreeLEDMon entries " << tem << endl;
  }

  Int_t    fNum     = 0;
  Double_t fAvgAmp  = 0;
  Double_t fHour    = 0;
  Double_t fRMS     = 0;
  treeAmp->SetBranchAddress("fChannelNum",&fNum);
  treeAmp->SetBranchAddress("fHour",&fHour);
  treeAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  treeAmp->SetBranchAddress("fRMS",&fRMS);

  // associate variables  also for LED
  treeLEDAmp->SetBranchAddress("fRefNum",&fNum);
  treeLEDAmp->SetBranchAddress("fHour",&fHour);
  treeLEDAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  treeLEDAmp->SetBranchAddress("fRMS",&fRMS);

  Int_t mod   = 0;
  Int_t col   = 0;
  Int_t row   = 0;
  Int_t gain  = 0;
  Int_t strip = 0;

  for (Int_t i=0; i<tem; ++i) {
    treeLEDAmp->GetEntry(i);
    emcCalibSignal->DecodeRefNum(fNum, &mod, &strip, &gain);
    if (debug) 
      Printf("fRefNum = %d, mod = %d, strip = %d, amp = %f, gain = %d", fNum, mod, strip, fAvgAmp, gain);
    ret->FillStrip(mod,gain,strip, fAvgAmp, fRMS);
  }

  for (Int_t i=0; i<tea; ++i) {
    treeAmp->GetEntry(i);
    emcCalibSignal->DecodeChannelNum(fNum, &mod, &col, &row, &gain);
    if (debug)
      Printf("fChannelNum = %d, mod = %d, col = %d, amp = %f, gain = %d", fNum, mod, col, fAvgAmp, gain);
    ret->FillLed(mod,gain,col, row, fAvgAmp, fRMS);
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
