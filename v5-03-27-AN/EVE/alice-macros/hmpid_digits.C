/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TBranch.h>
#include <TGeoMatrix.h>
#include <TTree.h>
#include <TStyle.h>
#include <TEveElement.h>
#include <TEveFrameBox.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TEveTrans.h>
#include <TEveQuadSet.h>

#include <AliHMPIDDigit.h>
#include <AliHMPIDv3.h>
#include <AliCluster3D.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

void hmpid_digits()
{
  const Char_t *name[]={ "HMPID0", "HMPID1", "HMPID2", "HMPID3",
			 "HMPID4", "HMPID5", "HMPID6" };

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("HMPID");

  TTree *dTree = rl->GetTreeD("HMPID", kFALSE);
  if (!dTree) return;

  TEveElementList* list = new TEveElementList("HMPID Digits");
  gEve->AddElement(list);

  gStyle->SetPalette(1, 0);

  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 3000);
  pal->SetMax(1000);
  TEveFrameBox    *box = new TEveFrameBox();
  box->SetAAQuadXY(0, 0, 0, 144, 121);
  box->SetFrameColor(kGray);

  TClonesArray* digits = new TClonesArray("AliHMPIDDigit");
  for (Int_t iCh = 0; iCh < 7; ++iCh)
  {
    TBranch *br = dTree->GetBranch(name[iCh]);
    br->SetAddress(&digits);
    br->GetEntry(0);

    TEveQuadSet* q = new TEveQuadSet(Form("Chamber %d", iCh));
    q->SetOwnIds(kTRUE);
    q->SetPalette(pal);
    q->SetFrame(box);
    q->SetAntiFlick(kTRUE);
    q->SetPickable(kTRUE);

    q->Reset(TEveQuadSet::kQT_RectangleXYFixedDimZ, kFALSE, 64);
    q->SetDefCoord(0);
    q->SetDefHeight(0.84f);
    q->SetDefWidth(0.8f);

    for(Int_t iDig = 0; iDig < digits->GetEntriesFast(); ++iDig)
    {
      AliHMPIDDigit *pDig = (AliHMPIDDigit*) digits->At(iDig);

      q->AddQuad(pDig->PadChX()*0.8f,  pDig->PadChY()*0.84f);
      q->QuadValue(TMath::Nint(pDig->Q()));
      q->QuadId(new AliHMPIDDigit(*pDig));
    }

    q->RefitPlex();

    TGeoHMatrix mat;
    AliHMPIDv3::IdealPosition(iCh, &mat);
    q->RefMainTrans().SetFrom(mat);
    q->RefMainTrans().Move3LF(-0.5f*144, -0.5f*121, 0);

    list->AddElement(q);
  }

  delete digits;
  rl->UnloadDigits("HMPID");

  gEve->Redraw3D();
}
