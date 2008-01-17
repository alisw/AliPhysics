#ifndef ALIEVE_TOFStrip_H
#define ALIEVE_TOFStrip_H

#include <TEveQuadSet.h>
#include <TEveElement.h>

#include <TEveRGBAPalette.h>
#include <TEveFrameBox.h>

#include <TGeoManager.h>
#include <TClonesArray.h>

#include <AliTOFGeometry.h>

namespace Alieve {

class TOFStrip : public TEveQuadSet
{
  TOFStrip(const TOFStrip&);            // Not implemented
  TOFStrip& operator=(const TOFStrip&); // Not implemented

private:
  void LoadQuads();
  
protected:
  virtual void InitModule();
  virtual void SetTrans();

  AliTOFGeometry* fTOFgeometry;

  TClonesArray *fTOFarray;

  Int_t fSector;
  Int_t fPlate;
  Int_t fStrip;

  Float_t  fDx;
  Float_t  fDz;

  TGeoManager *fGeoManager;

public:
  TOFStrip(const Text_t* n="TOFStrip", const Text_t* t=0);
  TOFStrip(TGeoManager *localGeoManager,
	    Int_t nSector, Int_t nPlate, Int_t nStrip);

  TOFStrip(TGeoManager *localGeoManager,
	    Int_t nSector, Int_t nPlate, Int_t nStrip,
	    TClonesArray *tofArray);
  virtual ~TOFStrip();

  static Bool_t    fgStaticInitDone;
  static void      InitStatics();

  static TEveFrameBox* fgTOFstripFrameBox;

  static TEveRGBAPalette* fgTOFstripPalette;

  ClassDef(TOFStrip, 1);
}; 
}
#endif
