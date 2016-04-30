#ifndef AliITSMFTSIMUCLUSTERSHAPER_H
#define AliITSMFTSIMUCLUSTERSHAPER_H
/* Copyright(c) 2007-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////
//                                                               //
// Class to generate the cluster shape in the ITSU simulation    //
// Author: Davide Pagano                                         //
///////////////////////////////////////////////////////////////////
#include <TRandom.h>
#include <TObject.h>
#include <sstream>

class AliITSMFTSimuClusterShaper : public TObject {

 public:
  AliITSMFTSimuClusterShaper();
  AliITSMFTSimuClusterShaper(const Int_t &cs);
  virtual ~AliITSMFTSimuClusterShaper();
  void FillClusterRandomly(Int_t *clusterConf);

  void AddNoisePixel();

  Int_t GetNRows() const {return fNrows;}
  Int_t GetNCols() const {return fNcols;}
  
  // just for debug...to be removed    
  inline std::string ShapeSting(Int_t cs, Int_t *cshape) const {
    std::stringstream out;
    for (Int_t i = 0; i < cs; ++i) {
      out << cshape[i];
      if (i < cs-1) out << " ";
    }
    return out.str();
  }    
    
 private:
  Int_t fNrows;
  Int_t fNcols;
  Int_t fNpixOn;
  
  ClassDef(AliITSMFTSimuClusterShaper,1)
};
#endif
