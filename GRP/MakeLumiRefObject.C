#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TString.h>
#include "AliLumiRef.h"
#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliLog.h"
#endif

void AddLumiRef(TObjArray* arr, int run, float sig,float eff,const TString trig, const TString comment);
TObjArray* refArr = 0;

void MakeLumiRefObject(TString cdbPath="local://./")
{
  // Macro to create OCDB entry with lumi.trigger references
  refArr = new TObjArray();

  AddLumiRef(refArr,      0,   62.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 118502,   47.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_0.90: 47mb=52 mb *0.91=sigma(INEL)*R(INT1/INEL) (arxiv: 1208.4968, fig.10 + table 3)");
  AddLumiRef(refArr, 118903,   62.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 121039,   47.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_0.90: 47mb=52 mb *0.91=sigma(INEL)*R(INT1/INEL) (arxiv: 1208.4968, fig.10 + table 3)");
  AddLumiRef(refArr, 121041,   62.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 126438,   62.0, 1.00, "CINT1-B-NOPF-ALLNOTRD",  "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 127719,   62.0, 1.00, "CINT1B-ABCE-NOPF-ALL",   "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 127731,   62.0, 1.00, "CINT1-B-NOPF-ALLNOTRD",  "pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 136849, 5970.0, 0.78, "C0SMH-B-NOPF-ALL",       "PbPb_2.76: (Oyama,2011-05-20,RunCond), sigma_hardronic = 7.64 b");
  AddLumiRef(refArr, 139328, 5970.0, 0.78, "C0SMH-B-NOPF-ALLNOTRD",  "PbPb_2.76: (Oyama,2011-05-20,RunCond), sigma_hardronic = 7.64 b");
  AddLumiRef(refArr, 145289,   57.0, 1.00, "CINT1-B-NOPF-ALLNOTRD",  "pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 146808,   57.0, 1.00, "CINT1-B-NOPF-ALL",       "pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 145815,   57.0, 1.00, "CINT1-B-NOPF-ALLNOTRD",  "pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 146857,   57.0, 1.00, "CINT1-B-NOPF-ALL",       "pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 146858,   57.0, 1.00, "CINT1-B-NOPF-ALLNOTRD",  "pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 148370,   54.0, 1.00, "CVBAND-B-NOPF-ALLNOTRD", "pp_7.00: 54.3mb (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 157079,   24.0, 0.44, "C0TVX-B-NOPF-ALLNOTRD",  "pp_7.00: 24mb=54.3mb*0.44=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 166477, 4100.0, 0.54, "CVLN-B-NOPF-ALLNOTRD",   "PbPb_2.76: (Martino,2013-03-15,RunCond)");
  AddLumiRef(refArr, 176658,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177146,   25.0, 0.33, "C0TVX-B-NOPF-CENTNOTRD", "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177148,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177150,   25.0, 0.33, "C0TVX-B-NOPF-CENTNOTRD", "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177166,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177167,   25.0, 0.33, "C0TVX-B-NOPF-CENTNOTRD", "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177168,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177169,   25.0, 0.33, "C0TVX-B-NOPF-CENTNOTRD", "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177173,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177174,   25.0, 0.33, "C0TVX-B-NOPF-CENTNOTRD", "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 177507,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 178018,   67.0, 0.90, "CINT1-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), CINT1/C0TVX=2.7 from 178052");
  AddLumiRef(refArr, 178030,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 178055,   25.0, 0.33, "C0TVX-B-NOPF-ALL",       "pp_8.00: (Artem, 2013-10-04,RunCond), vdM");
  AddLumiRef(refArr, 178062,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 179444,   25.0, 0.33, "C0TVX-S-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 180716,   56.0, 0.75, "CINT7-S-NOPF-ALLNOTRD",  "no C0TVX in these runs, taking VBAND cross section");
  AddLumiRef(refArr, 180721,   25.0, 0.33, "C0TVX-S-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 184845,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 184991,   25.0, 0.33, "C0TVX-S-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 188230, 1590.0, 0.76, "C0TVX-B-NOPF-ALLNOTRD",  "pPb_5.02: pilot. arxiv:1405.1849");
  AddLumiRef(refArr, 188367,   25.0, 0.33, "C0TVX-S-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 193693,   25.0, 0.33, "C0TVX-B-NOPF-ALLNOTRD",  "pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb");
  AddLumiRef(refArr, 195344, 1590.0, 0.76, "C0TVX-B-NOPF-ALLNOTRD",  "pPb_5.02: arxiv:1405.1849");
  AddLumiRef(refArr, 197470,   18.0, 0.39, "C0TVX-B-NOPF-ALLNOTRD",  "pp_2.76: 18mb=47.7mb*0.39=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)");
  AddLumiRef(refArr, 221835,   16.8, 0.32, "CADAND-B-NOPF-ALLNOTRD", "estimates from Martino");
  AddLumiRef(refArr, 221670,   30.0, 0.40, "C0TVX-B-NOPF-ALLNOTRD",  "estimates from Martino and MC");
  AddLumiRef(refArr, 223984,   50.0, 0.66, "CADAND-B-NOPF-ALLNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 223985,   30.0, 0.40, "C0TVX-B-NOPF-ALLNOTRD",  "estimates from Martino and MC");
  AddLumiRef(refArr, 226111,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 226116,   30.0, 0.40, "C0TVX-B-NOPF-ALLNOTRD",  "estimates from Martino and MC");
  AddLumiRef(refArr, 228910,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 229386,   30.0, 0.40, "C0TVX-B-NOPF-MUON",      "estimates from Martino and MC");
  AddLumiRef(refArr, 229409,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 229416,   30.0, 0.40, "C0TVX-B-NOPF-MUON",      "estimates from Martino and MC");
  AddLumiRef(refArr, 229894,   30.0, 0.40, "C0TVX-B-NOPF-ALLNOTRD",  "estimates from Martino and MC");
  AddLumiRef(refArr, 229942,   30.0, 0.40, "C0TVX-B-NOPF-MUON",      "estimates from Martino and MC");
  AddLumiRef(refArr, 232914,   30.0, 0.40, "C0TVX-B-NOPF-CENT",      "estimates from Martino and MC");
  AddLumiRef(refArr, 233910,   30.0, 0.40, "C0TVX-B-NOPF-ALLNOTRD",  "estimates from Martino and MC");
  AddLumiRef(refArr, 234051,   30.0, 0.40, "C0TVX-B-NOPF-CENT",      "estimates from Martino and MC");
  AddLumiRef(refArr, 238670,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 240151,   30.0, 0.40, "C0TVX-B-NOPF-MUON",      "estimates from Martino and MC");
  AddLumiRef(refArr, 240152,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 243374,   21.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 243399, 6700.0, 0.90, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 243985,   21.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 244913, 4600.0, 0.60, "C0V0M-B-NOPF-CENTNOTRD", "estimates from Cvetan and Alberica");
  AddLumiRef(refArr, 246995,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 256148,   30.0, 0.40, "C0TVX-B-NOPF-CENT",      "estimates from Martino and MC");
  AddLumiRef(refArr, 256158,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 265304, 1590.0, 0.76, "C0TVX-B-NOPF-CENTNOTRD", "taken from old p-Pb");
  AddLumiRef(refArr, 265587, 1715.0, 0.82, "C0TVX-B-NOPF-CENTNOTRD", "T0 efficiency estimate from MC, refSigma assuming cs(INEL) = 2092");
  AddLumiRef(refArr, 267132, 1590.0, 0.76, "C0TVX-B-NOPF-CENTNOTRD", "taken from old p-Pb");
  AddLumiRef(refArr, 267167,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "back to pp at 13TeV");
  AddLumiRef(refArr, 275924,   63.0, 0.84, "CINT7-B-NOPF-CENTNOTRD", "no T0, estimate from Martino");
  AddLumiRef(refArr, 276097,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "T0 recovered");
  AddLumiRef(refArr, 280179, 4000.0, 0.71, "C0V0M-B-NOPF-CENTNOTRD", "Xe-Xe, Martino: obtained as 5.6 b * L0b_0V0M / (L2a_INT7ZAC / live_time )");
  AddLumiRef(refArr, 280244,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "estimates from Martino and MC");
  AddLumiRef(refArr, 281962,   21.6, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "switch to pp at 5TeV");
  AddLumiRef(refArr, 282445,   30.0, 0.40, "C0TVX-B-NOPF-CENTNOTRD", "back to pp @ 13 TeV");
  
  // add new entries above
  //=======================================================================================
  refArr->Sort(); // Very important
  
  if (cdbPath.IsNull()) cdbPath = "local://./";

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbPath.Data());
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Ruben Shahoyan");
  metaData->SetComment("Ref. triggers and x-sections for lumi estimate");
  AliCDBId idLumiRef("GRP/CTP/LumiRef",0,AliCDBRunRange::Infinity());
  man->Put(refArr,idLumiRef, metaData);
  //
}

void AddLumiRef(TObjArray* arr, int run, float sig,float eff,const TString trig,const TString comment)
{
  if (sig<=0. || eff<=0. || trig.IsNull()) {
    AliFatalGeneralF("AddLumiRef","Wrong data: %d %f %f %s %s",run,sig,eff,trig.Data(),comment.Data());
  }
  AliLumiRef* ref = new AliLumiRef(trig,comment,run,sig,eff);
  arr->Add(ref);
}
