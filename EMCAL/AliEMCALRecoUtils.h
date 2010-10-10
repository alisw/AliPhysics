#ifndef ALIEMCALRECOUTILS_H
#define ALIEMCALRECOUTILS_H

/* $Id: AliEMCALRecoUtils.h 33808 2009-07-15 09:48:08Z gconesab $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class AliEMCALRecoUtils
// Some utilities to recalculate the cluster position or energy linearity
//
//
// Author:  Gustavo Conesa (LPSC- Grenoble) 
///////////////////////////////////////////////////////////////////////////////

//Root includes
#include "TNamed.h"

//AliRoot includes
class AliVCluster;
class AliVCaloCells;
#include "AliLog.h"
class AliEMCALGeoUtils;

class AliEMCALRecoUtils : public TNamed {
  
public:
  
  AliEMCALRecoUtils();
  AliEMCALRecoUtils(const AliEMCALRecoUtils&); 
  AliEMCALRecoUtils& operator=(const AliEMCALRecoUtils&); 
  virtual ~AliEMCALRecoUtils() {;}
  
  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3};
  
  //Position recalculation
  void     RecalculateClusterPosition(AliEMCALGeoUtils *geom, AliVCaloCells* cells, AliVCluster* clu, const Int_t iParticle); 
  void     GetMaxEnergyCell(AliEMCALGeoUtils *geom, AliVCaloCells* cells, AliVCluster* clu, 
                            Int_t & absId,  Int_t& iSupMod, Int_t& ieta, Int_t& iphi);
  
  Float_t  GetMisalShift(const Int_t i) const {
    if(i < 15 ){return fMisalShift[i]; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i)); return 0.;}
  }
  Float_t  *GetMisalShiftArray() {return fMisalShift; }

  void     SetMisalShift(const Int_t i, const Float_t shift) {
    if(i < 15 ){fMisalShift[i] = shift; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i));}
  }
  void     SetMisalShiftArray(Float_t * misal) 
  { for(Int_t i = 0; i < 15; i++)fMisalShift[i] = misal[i]; }

  //Non Linearity
  
  Float_t CorrectClusterEnergyLinearity(AliVCluster* clu);
  
  Float_t  GetNonLinearityParam(const Int_t i) const {
    if(i < 6 ){return fNonLinearityParams[i]; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i)); return 0.;}
  }
  void     SetNonLinearityParam(const Int_t i, const Float_t param) {
    if(i < 6 ){fNonLinearityParams[i] = param; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i));}
  }
  
  Int_t GetNonLinearityFunction() const {return fNonLinearityFunction;}
  void  SetNonLinearityFunction(Int_t fun) {fNonLinearityFunction = fun ;}
  
  void Print(const Option_t*) const;
  
private:
  
  Float_t fMisalShift[15];        // Shift parameters
  Int_t   fNonLinearityFunction;  // Non linearity function choice
  Float_t fNonLinearityParams[6]; // Parameters for the non linearity function

  ClassDef(AliEMCALRecoUtils, 1)
  
};

#endif // ALIEMCALRECOUTILS_H


