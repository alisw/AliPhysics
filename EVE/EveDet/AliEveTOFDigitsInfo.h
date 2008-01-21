// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_TOFDigitsInfo_H
#define ALIEVE_TOFDigitsInfo_H

#include <TEveVSD.h>

//#include <map>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliTOF.h>
#include <AliTOFGeometry.h>
#include <AliTOFDigitMap.h>


class AliEveTOFDigitsInfo : public TObject, public TEveRefCnt
  {
    AliEveTOFDigitsInfo(const AliEveTOFDigitsInfo&);            // Not implemented
    AliEveTOFDigitsInfo& operator=(const AliEveTOFDigitsInfo&); // Not implemented

  private:

  protected:

    void        SetTOFSegmentation();

  public:
    TTree*                fTree;
    TTree*                fNewTree;
    AliTOFGeometry*       fGeom;
    AliTOFDigitMap*       fTOFdigitMap;

    AliEveTOFDigitsInfo();
    virtual ~AliEveTOFDigitsInfo();

    void SetTree(TTree* tree);
    void LoadDigits();

    //TClonesArray* GetDigits(Int_t nSector,
    void GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip,
		   Int_t nPadZ, Int_t nPadX,
		   Int_t indexDigit[3]);
    TClonesArray* GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip);
    TClonesArray* GetDigits(Int_t nSector);
    void GetDigits();

    ClassDef(AliEveTOFDigitsInfo, 1);
  }; // endclass AliEveTOFDigitsInfo

#endif
