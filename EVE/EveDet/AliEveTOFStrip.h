// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveTOFStrip_H
#define AliEveTOFStrip_H

#include <TEveQuadSet.h>
#include <TEveElement.h>

#include <TEveRGBAPalette.h>
#include <TEveFrameBox.h>

#include <TGeoManager.h>
#include <TClonesArray.h>

#include <AliTOFGeometry.h>


class AliEveTOFStrip : public TEveQuadSet
{
protected:

  AliTOFGeometry *fTOFgeometry;

  TClonesArray   *fTOFarray;

  Short_t fThreshold;
  Int_t   fMaxVal;
  Int_t   fSector;
  Int_t   fPlate;
  Int_t   fStrip;

  Float_t fDx;
  Float_t fDz;

  TGeoManager *fGeoManager;

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

  ClassDef(AliEveTOFStrip, 0); // Representation of a TOF strip.
};
#endif
