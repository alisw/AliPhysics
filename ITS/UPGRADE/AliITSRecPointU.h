#ifndef ALIITSRECPOINTU_H
#define ALIITSRECPOINTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////////////////////
//  Class to set the local coordinates in ITS Upgrade recpoint //
/////////////////////////////////////////////////////////////////

#include <AliITSRecPoint.h>


class AliITSRecPointU : public AliITSRecPoint {

public :
  AliITSRecPointU();
  virtual ~AliITSRecPointU() {}; // distructor
  AliITSRecPointU(const AliITSRecPointU& pt);
  AliITSRecPointU& operator=(const AliITSRecPointU &source);

  void SetLocalCoord(Float_t x, Float_t z) {fXloc=x; fZloc=z;}
  void SetModule(Int_t i){fModule=i;} 
  Int_t GetModule(){return fModule;}

 protected:
  Int_t fModule;

  ClassDef(AliITSRecPointU,1)  // AliITSRecPointU class

};
#endif
