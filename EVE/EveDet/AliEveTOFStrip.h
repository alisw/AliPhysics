// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_TOFStrip_H
#define ALIEVE_TOFStrip_H

#include <TEveQuadSet.h>
#include <TEveElement.h>

#include <TEveRGBAPalette.h>
#include <TEveFrameBox.h>

#include <TGeoManager.h>
#include <TClonesArray.h>

#include <AliTOFGeometry.h>


class AliEveTOFStrip : public TEveQuadSet
{
  AliEveTOFStrip(const AliEveTOFStrip&);            // Not implemented
  AliEveTOFStrip& operator=(const AliEveTOFStrip&); // Not implemented

private:
  void LoadQuads();

protected:
  virtual void InitModule();
  virtual void SetTrans();

  AliTOFGeometry *fTOFgeometry;

  TClonesArray   *fTOFarray;

  Int_t fSector;
  Int_t fPlate;
  Int_t fStrip;

  Float_t  fDx;
  Float_t  fDz;

  TGeoManager *fGeoManager;

public:
  AliEveTOFStrip(const Text_t* n="AliEveTOFStrip", const Text_t* t=0);
  AliEveTOFStrip(TGeoManager *localGeoManager,
		 Int_t nSector, Int_t nPlate, Int_t nStrip);

  AliEveTOFStrip(TGeoManager *localGeoManager,
		 Int_t nSector, Int_t nPlate, Int_t nStrip,
		 TClonesArray *tofArray);
  virtual ~AliEveTOFStrip();

  static Bool_t    fgStaticInitDone;
  static void      InitStatics();

  static TEveFrameBox* fgTOFstripFrameBox;

  static TEveRGBAPalette* fgTOFstripPalette;

  ClassDef(AliEveTOFStrip, 1);
};
#endif
