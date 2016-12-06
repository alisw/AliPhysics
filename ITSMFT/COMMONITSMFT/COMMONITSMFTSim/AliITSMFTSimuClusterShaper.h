#ifndef AliITSMFTSIMUCLUSTERSHAPER_H
#define AliITSMFTSIMUCLUSTERSHAPER_H
/* Copyright(c) 2007-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////
//                                                               //
// Class to generate the cluster shape in the ITSU simulation    //
// Author: Davide Pagano                                         //
///////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <sstream>
#include "AliITSMFTClusterShape.h"

class AliITSMFTSimuClusterShaper : public TObject {

 public:
  AliITSMFTSimuClusterShaper();
  AliITSMFTSimuClusterShaper(const UInt_t &cs);
  virtual ~AliITSMFTSimuClusterShaper();
  void FillClusterRandomly();
  void AddNoisePixel();

  inline UInt_t  GetNRows() {return fCShape->GetNRows();}
  inline UInt_t  GetNCols() {return fCShape->GetNCols();}
  inline UInt_t* GetShape() {return fCShape->GetShape();}

  inline std::string ShapeSting(UInt_t cs, UInt_t *cshape) const {
    std::stringstream out;
    for (Int_t i = 0; i < cs; ++i) {
      out << cshape[i];
      if (i < cs-1) out << " ";
    }
    return out.str();
  }

 private:
  UInt_t fNpixOn;
  AliITSMFTClusterShape *fCShape;

  ClassDef(AliITSMFTSimuClusterShaper,1)
};
#endif
