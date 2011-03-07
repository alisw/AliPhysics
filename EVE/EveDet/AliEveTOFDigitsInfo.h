#ifndef ALIEVETOFDIGITSINFO_H
#define ALIEVETOFDIGITSINFO_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

//
// Class to map TOF digit/raw data information
//

#include <TObject.h>

#include <TEveUtil.h>

class TClonesArray;
class TTree;

class AliRawReader;

class AliTOFGeometry;
class AliTOFDigitMap;

class AliEveTOFDigitsInfo : public TObject, public TEveRefCnt
  {

  public:

    AliEveTOFDigitsInfo();
    virtual ~AliEveTOFDigitsInfo();
    
    void SetTree(TTree * const tree);
    void ReadRaw(AliRawReader* rawReader, Int_t newDecoder=2);
    void LoadDigits();

    void GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip,
		   Int_t nPadZ, Int_t nPadX,
		   Int_t indexDigit[3]);
    TClonesArray* GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip);
    TClonesArray* GetDigits(Int_t nSector);

    Int_t IsStripFilled(Int_t iSector, Int_t iPlate, Int_t iStrip);
    Int_t GetTOFInfos() const;
    AliTOFGeometry * GetTOFgeometry() const {return fGeom;};
    //void GetDigits();

    TTree* GetTree() {return fTree;};

    AliTOFDigitMap* GetTOFdigitMap() const { return fTOFdigitMap;};

  protected:

    AliEveTOFDigitsInfo(const AliEveTOFDigitsInfo&);            // Not implemented
    AliEveTOFDigitsInfo& operator=(const AliEveTOFDigitsInfo&); // Not implemented

  private:

    TTree*           fTree;        // pointer to TOF digit tree
    TTree*           fNewTree;     // pointer to TOF digit tree
    AliTOFGeometry*  fGeom;        // pointer to AliTOFGeometry class
    AliTOFDigitMap*  fTOFdigitMap; // pointer to AliTOFDIgitMap class

    ClassDef(AliEveTOFDigitsInfo, 1);
}; // endclass AliEveTOFDigitsInfo

#endif
