// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the T0 detector                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliEveT0Module.h"

#include <EveBase/AliEveEventManager.h>

#include <AliT0digit.h>
#include <AliRawReader.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>
#include <AliT0RawReader.h>
#include <AliCDBManager.h>
#include <AliCDBStorage.h>

#include <TArrayI.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom.h>
#include <TEveManager.h>


ClassImp(AliEveT0Module)

/******************************************************************************/
AliEveT0Module::AliEveT0Module(const Text_t* n, Int_t sigType, AliT0digit *digits, AliT0RawReader *start)
  : TEveQuadSet(n), fSigType(sigType), fDigits(digits), fStart(start)
{
  //
  // Default constructor
  //

}

/******************************************************************************/
AliEveT0Module::~AliEveT0Module()
{

}

void AliEveT0Module::LoadRaw(TString fileName, Int_t ievt)
{
  AliT0digit *digits = 0;
  AliRawReader *reader = new AliRawReaderRoot(fileName,ievt);
  reader->LoadEquipmentIdsMap("T0map.txt");
  reader->RequireHeader(kTRUE);
  AliT0RawReader *start = new AliT0RawReader(reader);
  Int_t allData[110][5];
  TRandom r(0);
  //  cout<<ievt<<endl;
  TEveRGBAPalette* rawPalette  = new TEveRGBAPalette(0, 3000);
  rawPalette->SetLimits(1, 3000); // Set proper raw time range.
  TEveQuadSet* raw_a = new AliEveT0Module("T0_RAW_A", 2,digits, start); raw_a->SetPalette(rawPalette);
  raw_a->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  TEveQuadSet* raw_c = new AliEveT0Module("T0_RAW_C", 3,digits, start); raw_c->SetPalette(rawPalette);
  raw_c->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  Float_t angle  = 2 * TMath::Pi() / 12;
  start->Next();
  for (Int_t i=0; i<110; i++)
  {
    for (Int_t iHit=0; iHit<5; iHit++)
    {
      allData[i][iHit]= start->GetData(i,iHit);
      if (allData[i][iHit] != 0)  cout<<"event"<<ievt<<" i "<< i<<" "<<allData[i][iHit] - allData[0][0]<<endl;
    }
  }
  for (Int_t i=0; i<12; i++)
  {
    Float_t x = 6.5 * TMath::Sin(i * angle);
    Float_t y = 6.5 * TMath::Cos(i * angle);
    raw_a->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    raw_a->QuadValue(start->GetData(i+37,0)-start->GetData(0,0));
    raw_c->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    raw_c->QuadValue(start->GetData(i+25,0)-start->GetData(0,0));
  }


  raw_a->RefitPlex();
  raw_c->RefitPlex();

  TEveTrans& ta_a = raw_a->RefHMTrans();
  ta_a.SetPos(0, 0, 373);
  TEveTrans& tc_c = raw_c->RefHMTrans();
  tc_c.SetPos(0, 0, -69.7);

  gEve->AddElement(raw_a);
  gEve->AddElement(raw_c);
  gEve->Redraw3D();
}

/******************************************************************************/
void AliEveT0Module::MakeModules(AliT0digit *digits)
{
  TRandom r(0);
  TArrayI ADC(24);
  TArrayI TDC(24);

  digits->GetQT1(ADC);
  digits->GetTimeCFD(TDC);
  //    printf("%3d\n",besttimeright);
  for (Int_t i=0;i<24; i++){
    printf("%3d %3d\n  ",ADC[i], TDC[i]);
  }

  TEveRGBAPalette* adcPalette  = new TEveRGBAPalette(5, 1024);
  adcPalette->SetLimits(1, 1024); // Set proper ADC range.
  TEveRGBAPalette* tdcPalette  = new TEveRGBAPalette(0, 9999);
  tdcPalette->SetLimits(1, 9999); // Set proper TDC range.

  TEveQuadSet* qa = new AliEveT0Module("T0A_ADC", 0, digits); qa->SetPalette(adcPalette);
  TEveQuadSet* qc = new AliEveT0Module("T0C_ADC", 0, digits); qc->SetPalette(adcPalette);
  TEveQuadSet* qat = new AliEveT0Module("T0A_TDC", 1, digits); qat->SetPalette(tdcPalette);
  TEveQuadSet* qct = new AliEveT0Module("T0C_TDC", 1, digits); qct->SetPalette(tdcPalette);

  Float_t angle  = 2 * TMath::Pi() / 12;

  qa->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  qc->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  qat->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  qct->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);



  for (Int_t i=0; i<12; i++) {
    Float_t x = 6.5 * TMath::Sin(i * angle);
    Float_t y = 6.5 * TMath::Cos(i * angle);

    qa->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qa->QuadValue(ADC[i+12]);
    //    qa->QuadId(new TNamed(Form("PMT   with idx=%d", i), "PMT's aplitude in side A."));

    qat->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qat->QuadValue(TDC[i+12]);
    //    qat->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's time in side A."));

    qc->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qc->QuadValue(ADC[i]);
    // qc->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's amplitude in side C."));

    qct->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qct->QuadValue(TDC[i]);
    // qct->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's time in side C."));
  }

  qa->RefitPlex();
  qc->RefitPlex();
  qat->RefitPlex();
  qct->RefitPlex();

  TEveTrans& ta = qa->RefHMTrans();
  ta.SetPos(0, 0, 373);
  TEveTrans& tc = qc->RefHMTrans();
  tc.SetPos(0, 0, -69.7);

  TEveTrans& tat = qat->RefHMTrans();
  tat.SetPos(0, 0, 373);
  TEveTrans& tct = qct->RefHMTrans();
  tct.SetPos(0, 0, -69.7);

  gEve->AddElement(qa);
  gEve->AddElement(qc);
  gEve->AddElement(qat);
  gEve->AddElement(qct);

  gEve->Redraw3D();
}

/******************************************************************************/

void AliEveT0Module::DigitSelected(Int_t idx)
{
  // Override control-click from TEveQuadSet

  DigitBase_t* qb   = GetDigit(idx);
  if (fSigType == 0) { //ADC
    printf("adc====================\n");
    Int_t   besttimeright = fDigits->BestTimeA();
    Int_t   besttimeleft = fDigits->BestTimeC();
    Int_t   meantime = fDigits->MeanTime();
    Int_t   timediff = fDigits->TimeDiff();
    Int_t   summult = fDigits->SumMult();

    printf("besttimeA=%3d\n",besttimeright);
    printf("besttimeC=%3d\n",besttimeleft);
    printf("meantime=%3d\n",meantime);
    printf("timediff=%3d\n",timediff);
    printf("summult=%3d\n",summult);

    printf("  idx=%d, amplitude=%d\n",  idx, qb->fValue);

  }
  if (fSigType == 1) {
    printf("tdc====================\n");

    Int_t   besttimeright = fDigits->BestTimeA();
    Int_t   besttimeleft = fDigits->BestTimeC();
    Int_t   meantime = fDigits->MeanTime();
    Int_t   timediff = fDigits->TimeDiff();
    Int_t   summult = fDigits->SumMult();

    printf("besttimeA=%3d\n",besttimeright);
    printf("besttimeC=%3d\n",besttimeleft);
    printf("meantime=%3d\n",meantime);
    printf("timediff=%3d\n",timediff);
    printf("summult=%3d\n",summult);

    printf("  idx=%d, amplitude=%d\n",  idx, qb->fValue);
  }
  if (fSigType == 2) {
    printf("raw====================\n");
    printf("besttimeA=%3d\n",fStart->GetData(51,0)-fStart->GetData(0,0));
    printf("besttimeC=%3d\n",fStart->GetData(52,0)-fStart->GetData(0,0));
    printf("meantime=%3d\n",fStart->GetData(49,0)-fStart->GetData(0,0));
    printf("amplitude= %3d\n",fStart->GetData(idx+1,0));

    printf("  idx=%d, time %d\n",  idx, qb->fValue);
  }
  if (fSigType == 3) {
    printf("raw====================\n");
    printf("besttimeA=%3d\n",fStart->GetData(51,0)-fStart->GetData(0,0));
    printf("besttimeC=%3d\n",fStart->GetData(52,0)-fStart->GetData(0,0));
    printf("meantime=%3d\n",fStart->GetData(49,0)-fStart->GetData(0,0));
    printf("amplitude= %3d\n",fStart->GetData(idx+13,0));

    printf("  idx=%d, time %d\n",  idx, qb->fValue);
  }

}
