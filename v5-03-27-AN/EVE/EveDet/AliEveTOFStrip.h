#ifndef ALIEVETOFSTRIP_H
#define ALIEVETOFSTRIP_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

//
// Class to visualize the TOF digit information
// in TOF strip frame
//

#include <TEveQuadSet.h>

class TClonesArray;

class AliTOFGeometry;

class AliEveTOFStrip : public TEveQuadSet
{
public:

  virtual void InitModule();
  virtual void SetTrans();

  AliEveTOFStrip(const Text_t* n="AliEveTOFStrip", const Text_t* t=0);
  AliEveTOFStrip(TGeoManager *localGeoManager,
		 Int_t nSector, Int_t nPlate, Int_t nStrip);

  AliEveTOFStrip(TGeoManager *localGeoManager,
		 Int_t nSector, Int_t nPlate, Int_t nStrip,
		 TClonesArray *tofArray);
  virtual ~AliEveTOFStrip();

  void SetThreshold(Short_t t);
  void SetMaxVal(Int_t mv);
  Short_t GetThreshold() const { return fThreshold; }
  Int_t   GetMaxVal()    const { return fMaxVal; }
  virtual void DigitSelected(Int_t idx);

protected:
  static Bool_t    fgStaticInitDone; // Has initialization of static variables been done.
  static void      InitStatics();    // Initialize static variables.

  static TEveFrameBox    *fgTOFstripFrameBox; // Shared box-frame for all strips.
  static TEveRGBAPalette *fgTOFstripPalette;  // Shared palette.

private:
  void LoadQuads();

  AliEveTOFStrip(const AliEveTOFStrip&);            // Not implemented
  AliEveTOFStrip& operator=(const AliEveTOFStrip&); // Not implemented


  AliTOFGeometry *fTOFgeometry; // pointer to TOF geometry

  TClonesArray   *fTOFarray;    // pointer to TOF digits array

  Short_t fThreshold; // threshold to cut on visualization
  Int_t   fMaxVal;    // max value to cut on visualization
  Int_t   fSector;    // TOF sector index
  Int_t   fPlate;     // TOF module index
  Int_t   fStrip;     // TOF strip index

  Float_t fDx;    // x position of TOF digit (in TOF strip RF)
  Float_t fDz;    // z position of TOF digit (in TOF strip RF)

  TGeoManager *fGeoManager; // pointer to the ALICE geometry

  ClassDef(AliEveTOFStrip, 0); // Representation of a TOF strip.
};
#endif
