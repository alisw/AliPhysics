/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoMatrix.h>
#include <TStyle.h>
#include <TEveElement.h>
#include <TEveFrameBox.h>
#include <TEveManager.h>
#include <TEveRGBAPalette.h>
#include <TEveTrans.h>
#include <TEveQuadSet.h>

#include <AliHMPIDDigit.h>
#include <AliHMPIDv3.h>
#include <AliHMPIDRawStream.h>
#include <AliRawReader.h>
#include <AliEveEventManager.h>
#endif

void hmpid_raw()
{
    printf("*** RAW HMPID ***");
    
  const Char_t *name[] = { "HMPID0", "HMPID1", "HMPID2", "HMPID3",
			   "HMPID4", "HMPID5", "HMPID6" };

  AliRawReader *rawReader = AliEveEventManager::AssertRawReader();
  AliHMPIDRawStream stream(rawReader);    

  TEveElementList* list = new TEveElementList("HMPID Raw");
  gEve->AddElement(list);

  gStyle->SetPalette(1, 0);

  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 3000);
  pal->SetMax(1000);
  TEveFrameBox    *box = new TEveFrameBox();
  box->SetAAQuadXY(0, 0, 0, 144, 121);
  box->SetFrameColor(kGray);

  TEveQuadSet* ms[7];
  for (Int_t iCh = 0; iCh < 7; ++iCh)
  {
    ms[iCh] = new TEveQuadSet(Form("Chamber %d", iCh));

    TEveQuadSet* q = ms[iCh];
    q->SetOwnIds(kTRUE);
    q->SetPalette(pal);
    q->SetFrame(box);
    q->SetAntiFlick(kTRUE);
    q->SetPickable(kTRUE);

    q->Reset(TEveQuadSet::kQT_RectangleXYFixedDimZ, kFALSE, 64);
    q->SetDefCoord(0);
    q->SetDefHeight(0.84f);
    q->SetDefWidth(0.8f);
  }

  while (stream.Next())
  {
    Int_t ch = AliHMPIDParam::DDL2C(stream.GetDDLNumber());
    TEveQuadSet* q = ms[ch];

    for (Int_t iPad = 0; iPad < stream.GetNPads(); ++iPad)
    {
      AliHMPIDDigit dig(stream.GetPadArray()[iPad],stream.GetChargeArray()[iPad]);

      q->AddQuad(dig.PadChX()*0.8f,  dig.PadChY()*0.84f);
      q->QuadValue(TMath::Nint(dig.Q()));
      q->QuadId(new AliHMPIDDigit(dig));
    }
  }

  for (Int_t iCh = 0; iCh < 7; ++iCh)
  {
    TEveQuadSet* q = ms[iCh];

    q->RefitPlex();

    TGeoHMatrix mat;
    AliHMPIDv3::IdealPosition(iCh, &mat);
    q->RefMainTrans().SetFrom(mat);
    q->RefMainTrans().Move3LF(-0.5*144, -0.5*121, 0);

    list->AddElement(q);
  }

  gEve->Redraw3D();
}
