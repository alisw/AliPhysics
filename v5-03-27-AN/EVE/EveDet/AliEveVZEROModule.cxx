/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The drawing module for the VZERO detector                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliEveVZEROModule.h"

#include <TH1.h>
#include <TEveManager.h>
#include <EveBase/AliEveEventManager.h>
#include <TEveTrans.h>

#include <AliRawReader.h>
#include <AliVZERORawStream.h>
#include <AliESDEvent.h>

static const Float_t RadC[] = { 4.6, 7.1, 11.5, 18.9, 31.4 };
static const Float_t RadA[] = { 4.2, 7.6, 13.8, 22.7, 41.4 };
static const Float_t RadEps = 0.4;
static const Float_t PhiEps = 0.025;
static const Float_t PhiStp = TMath::TwoPi()/8.0;

ClassImp(AliEveVZEROModule)

/******************************************************************************/
AliEveVZEROModule::AliEveVZEROModule(const Text_t* n, Bool_t side)
  : TEveQuadSet(n),
    fStream(NULL),
    fSampleIndex(AliVZERORawStream::kNEvOfInt/2),
    fIsASide(side)
{
  //
  // Default constructor
  //
  TEveRGBAPalette* rawPalette  = new TEveRGBAPalette(0, 1023);
  rawPalette->SetLimits(0, 1023);
  SetPalette(rawPalette);
  Reset(TEveQuadSet::kQT_FreeQuad, kFALSE, 32);
}

/******************************************************************************/
AliEveVZEROModule::~AliEveVZEROModule()
{
  //
  // Destructor
  //
  delete fStream;
}

/******************************************************************************/
void AliEveVZEROModule::LoadRaw(AliRawReader *rawReader)
{
  //
  // Load VZERO raw-data
  //
  if (fStream) delete fStream;
  fStream = new AliVZERORawStream(rawReader);
  if (!fStream->Next()) {
    delete fStream;
    fStream = NULL;
    return;
  }

  for (Int_t iChannel=0; iChannel < AliVZERORawStream::kNChannels; ++iChannel) {
    Int_t offChannel = fStream->GetOfflineChannel(iChannel);
    Float_t minR,maxR,minP,maxP;
    if (fIsASide) {
      if (offChannel < 32) continue;

      Int_t ri = (offChannel-32) / 8;
      Int_t pi = (offChannel-32) % 8;
      minR = RadA[ri]  + RadEps, maxR = RadA[ri+1] - RadEps;
      minP = pi*PhiStp + PhiEps, maxP = (pi+1)*PhiStp - PhiEps;
    }
    else {
      if (offChannel >= 32) continue;

      Int_t ri = offChannel / 8;
      Int_t pi = offChannel % 8;
      minR = RadC[ri]  + RadEps, maxR = RadC[ri+1] - RadEps;
      minP = pi*PhiStp + PhiEps, maxP = (pi+1)*PhiStp - PhiEps;
    }
    Float_t v[12];
    v[ 0] = minR*TMath::Cos(minP); v[ 1] = minR*TMath::Sin(minP); v[ 2] = 0;
    v[ 3] = maxR*TMath::Cos(minP); v[ 4] = maxR*TMath::Sin(minP); v[ 5] = 0;
    v[ 6] = maxR*TMath::Cos(maxP); v[ 7] = maxR*TMath::Sin(maxP); v[ 8] = 0;
    v[ 9] = minR*TMath::Cos(maxP); v[10] = minR*TMath::Sin(maxP); v[11] = 0;

    AddQuad(v);
    QuadValue(fStream->GetPedestal(iChannel,fSampleIndex));
  }

  if (fIsASide) 
    RefMainTrans().SetPos(0, 0, 324);
  else
    RefMainTrans().SetPos(0, 0, -84);

  gEve->AddElement(this);
  gEve->Redraw3D();
}

/******************************************************************************/
void AliEveVZEROModule::DigitSelected(Int_t idx)
{
  //
  // Override control-click from TEveQuadSet
  //
  if (!fStream) return;

  Int_t iPMT = idx;
  if (!fIsASide) iPMT += 32;
  printf("PMT = %2d ADC = ",iPMT);
  TH1S *hADC = new TH1S(Form("VZERO_ADC_%d",iPMT),
			Form("ADC samples for PMT %d, Time/Width = %3.1f/%3.1f ns",
			     iPMT,
			     (Float_t)fStream->GetTime(iPMT)/10.,(Float_t)fStream->GetWidth(iPMT)/10.),
			AliVZERORawStream::kNEvOfInt,-0.5,(Float_t)AliVZERORawStream::kNEvOfInt-0.5);
  hADC->SetXTitle("Sample index");
  hADC->SetYTitle("ADC value");
  hADC->SetStats(kFALSE);
  for (Int_t iEv = 0; iEv < AliVZERORawStream::kNEvOfInt; ++iEv) {
    printf("%4d ",fStream->GetPedestal(iPMT,iEv));
    hADC->SetBinContent(iEv+1,fStream->GetPedestal(iPMT,iEv));
  }
  printf("\nTime = %3.1f ns  Width = %3.1f ns\n",(Float_t)fStream->GetTime(iPMT)/10.,(Float_t)fStream->GetWidth(iPMT)/10.);
  hADC->Draw();
  gPad->Modified();
  gPad->Update();
}

/******************************************************************************/
void AliEveVZEROModule::SetSampleIndex(Int_t index)
{
  fSampleIndex = index;
  if (!fStream) return;

  for(Int_t idx = 0; idx < (AliVZERORawStream::kNChannels/2); ++idx) {
    Int_t iPMT = idx;
    if (!fIsASide) iPMT += 32;
    DigitBase_t *qb  = GetDigit(idx);
    qb->fValue = fStream->GetPedestal(iPMT,fSampleIndex);
  }
}
