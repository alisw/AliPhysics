// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveTOFDigitsInfo_H
#define AliEveTOFDigitsInfo_H

#include <TEveVSD.h>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliRawReader.h>

#include <AliTOF.h>
#include <AliTOFGeometry.h>
#include <AliTOFDigitMap.h>

  
class AliEveTOFDigitsInfo : public TObject, public TEveRefCnt
{
public:
  TTree*                fTree;
  TTree*                fNewTree;
  AliTOFGeometry*       fGeom;
  AliTOFDigitMap*       fTOFdigitMap;

  AliEveTOFDigitsInfo();
  virtual ~AliEveTOFDigitsInfo();
    
  void SetTree(TTree* tree);
  void ReadRaw(AliRawReader* rawReader, Int_t nEvent, Bool_t newDecoder=kTRUE);
  void LoadDigits();

  void GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip,
                 Int_t nPadZ, Int_t nPadX,
                 Int_t indexDigit[3]);
  TClonesArray* GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip);
  TClonesArray* GetDigits(Int_t nSector);

  Int_t IsStripFilled(Int_t iSector, Int_t iPlate, Int_t iStrip);
  Int_t GetTOFInfos() const;
  //void GetDigits();

  TTree* GetTree() { return fTree; }

protected:
  void SetTOFSegmentation();
 
private:
  AliEveTOFDigitsInfo(const AliEveTOFDigitsInfo&);            // Not implemented
  AliEveTOFDigitsInfo& operator=(const AliEveTOFDigitsInfo&); // Not implemented

  ClassDef(AliEveTOFDigitsInfo, 0);
}; // endclass AliEveTOFDigitsInfo

#endif
