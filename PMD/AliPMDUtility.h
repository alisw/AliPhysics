#ifndef PMDUtility_H
#define PMDUtility_H
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Utility class for PMD                              //
//                                                     //
//-----------------------------------------------------//

#include <math.h>
#include "Riostream.h"
#include "TMath.h"
#include "Rtypes.h"

class AliPMDUtility
{
  
 protected:
  Float_t fPx, fPy, fPz;
  Float_t fTheta, fEta, fPhi;

 public:
  AliPMDUtility();
  AliPMDUtility(Float_t /* Px */, Float_t /* Py */, Float_t /* Pz */);
  virtual ~AliPMDUtility();

  void HexGeomCellPos(Int_t /* ism */, Int_t /* xpad */, Int_t /* ypad */,
		Float_t & /* xpos */, Float_t & /* ypos */);
  void RectGeomCellPos(Int_t /* ism */, Int_t /* ium */, 
		       Int_t /* xpad */, Int_t /* ypad */,
		       Float_t & /* xpos */, Float_t & /* ypos */);
  void SetPxPyPz(Float_t /* Px */, Float_t /* Py */, Float_t /* Pz */);
  void SetXYZ(Float_t /* xPos */, Float_t /* yPos */, Float_t /* zPos */);
  void CalculateEta();
  void CalculatePhi();
  void CalculateEtaPhi();
  Float_t GetTheta() const;
  Float_t GetEta() const;
  Float_t GetPhi() const;
  
  ClassDef(AliPMDUtility,1)
};

#endif
