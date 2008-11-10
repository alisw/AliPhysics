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
#include <AliRawReaderRoot.h>
#include <AliT0RawReader.h>
#include <AliCDBStorage.h>

#include <TEveTrans.h>
#include <TEveManager.h>
#include <TTree.h>
#include <TArrayI.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom.h>


ClassImp(AliEveT0Module)

/******************************************************************************/
AliEveT0Module::AliEveT0Module(const Text_t* n, Int_t sigType, AliT0digit *digits, AliT0RawReader *start)
  : TEveQuadSet(n), fSigType(sigType), fDigits(digits), fStart(start)
{
  //
  // Default constructor
  //
}

void AliEveT0Module::PrintEventInfo()
{
  printf("Blabla\n");
}

/******************************************************************************/

void AliEveT0Module::LoadRaw(AliRawReader* reader)
{
  // Load raw-data from file.

  AliT0digit *digits = 0;

  // ??? How / when is this called during reco from raw?
  reader->LoadEquipmentIdsMap("T0map.txt");
  reader->RequireHeader(kTRUE);

  AliT0RawReader *start = new AliT0RawReader(reader);
  Int_t allData[110][5];
  TRandom r(0);

  TEveRGBAPalette *rawPalette    = new TEveRGBAPalette(0, 3000);
  TEveRGBAPalette *vertexPalette = new TEveRGBAPalette(-100, 100);

  TEveFrameBox *box = new TEveFrameBox();
  {
    Float_t  frame[3*36];
    Float_t *p = frame;
    for (Int_t i = 0; i < 36; ++i, p += 3) {
      p[0] = 8.0f * TMath::Cos(TMath::TwoPi()*i/36);
      p[1] = 8.0f * TMath::Sin(TMath::TwoPi()*i/36);
      p[2] = 0;
    }
    box->SetQuadByPoints(frame, 36);
  }
  box->SetFrameColor(kGray);


  rawPalette->SetLimits(1, 3000); // Set proper raw time range.
  TEveQuadSet* rawA = new AliEveT0Module("T0_RAW_A", 2, digits, start);
  rawA->SetPalette(rawPalette);
  rawA->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
  rawA->SetFrame(box);
  TEveQuadSet* rawC = new AliEveT0Module("T0_RAW_C", 3, digits, start);
  rawC->SetPalette(rawPalette);
  rawC->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32); 
  rawC->SetFrame(box);

  TEveQuadSet* vertexT0 = new AliEveT0Module("T0_Vertex", 5, digits, start);
  vertexT0->SetPalette(vertexPalette);
  vertexT0->Reset(TEveQuadSet::kQT_HexagonXY, kFALSE, 32);
 
  Float_t angle  = 2 * TMath::Pi() / 12;
  start->Next();
  for (Int_t i=0; i<110; i++)
  {
    for (Int_t iHit=0; iHit<5; iHit++)
    {
      allData[i][iHit]= start->GetData(i,iHit);
      if (allData[i][iHit] != 0) {
        using namespace std;
        cout <<"event "<< AliEveEventManager::GetCurrent()->GetEventId()
	     <<" i "<< i <<" "<< allData[i][iHit] - allData[0][0] <<endl;
      }
    }
  }
  Float_t zvertex= (allData[51][0] - allData[52][0])/2*25*2.99752/100;
  using namespace std;
  cout<<"zvertex= "<< zvertex <<endl;
  for (Int_t i=0; i<12; i++)
  {
    Float_t x = 6.5 * TMath::Sin(i * angle);
    Float_t y = 6.5 * TMath::Cos(i * angle);
    rawA->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    rawA->QuadValue(start->GetData(i+37,0)-start->GetData(0,0));
    rawC->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    rawC->QuadValue(start->GetData(i+25,0)-start->GetData(0,0));
  }
    vertexT0->AddHexagon(0, 0, 0, 1);
    vertexT0->QuadValue(TMath::Nint(zvertex));

  rawA->RefitPlex();
  rawC->RefitPlex();
  vertexT0->RefitPlex();

  TEveTrans& taA = rawA->RefMainTrans();
  taA.SetPos(0, 0, 373);
  TEveTrans& tcC = rawC->RefMainTrans();
  tcC.SetPos(0, 0, -69.7);

  TEveTrans& tver = vertexT0->RefMainTrans();
  tver.SetPos(0, 0, zvertex);

  gEve->AddElement(rawA);
  gEve->AddElement(rawC);
  gEve->AddElement(vertexT0);
  gEve->Redraw3D();
}

/******************************************************************************/

void AliEveT0Module::MakeModules(AliT0digit *digits)
{
  // Make modules for digits.

  TRandom r(0);
  TArrayI adc(24);
  TArrayI tdc(24);

  digits->GetQT1(adc);
  digits->GetTimeCFD(tdc);
  //    printf("%3d\n",besttimeright);
  for (Int_t i=0;i<24; i++){
    printf("%3d %3d\n  ",adc[i], tdc[i]);
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
    qa->QuadValue(adc[i+12]);
    //    qa->QuadId(new TNamed(Form("PMT   with idx=%d", i), "PMT's aplitude in side A."));

    qat->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qat->QuadValue(tdc[i+12]);
    //    qat->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's time in side A."));

    qc->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qc->QuadValue(adc[i]);
    // qc->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's amplitude in side C."));

    qct->AddHexagon(x, y, r.Uniform(-0.1, 0.1), 1.0);
    qct->QuadValue(tdc[i]);
    // qct->QuadId(new TNamed(Form("Quad with idx=%d", i), "PMT's time in side C."));
  }

  qa->RefitPlex();
  qc->RefitPlex();
  qat->RefitPlex();
  qct->RefitPlex();

  TEveTrans& ta = qa->RefMainTrans();
  ta.SetPos(0, 0, 373);
  TEveTrans& tc = qc->RefMainTrans();
  tc.SetPos(0, 0, -69.7);

  TEveTrans& tat = qat->RefMainTrans();
  tat.SetPos(0, 0, 373);
  TEveTrans& tct = qct->RefMainTrans();
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
  if (fSigType == 5) {

    printf("vertex====================\n");
    printf("  idx=%d, zvertex pozition %d\n",  idx, qb->fValue);

  }
}
