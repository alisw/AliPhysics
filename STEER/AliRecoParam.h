#ifndef ALIRECOPARAM_H
#define ALIRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Steering Class for reconstruction parameters                              //
// Revision: cvetan.cheshkov@cern.ch 12/06/2008                              //
// Its structure has been revised and it is interfaced to AliRunInfo and     //
// AliEventInfo.                                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class AliDetectorRecoParam;
class AliRunInfo;
class AliEventInfo;

class AliRecoParam : public TObject
{

 public: 
  AliRecoParam();
  virtual ~AliRecoParam();  
  //
  enum {
    kNSpecies = 5,   // number of event species
    kNDetectors = 16 // number of detectors (last one is GRP!)
  };
  enum EventSpecie_t {kDefault = 1,
		      kLowMult = 2,
		      kHighMult = 4,
		      kCosmic = 8,
		      kCalib = 16};

  virtual void                  Print(Option_t *option="") const;
  const TObjArray              *GetDetRecoParamArray(Int_t iDet) const { return fDetRecoParams[iDet]; }
  void                          SetEventSpecie(const AliRunInfo*/*runInfo*/, const AliEventInfo &/*evInfo*/);
  EventSpecie_t                 GetEventSpecie() const { return fEventSpecie; }
  const AliDetectorRecoParam   *GetDetRecoParam(Int_t iDet) const;
  void                          AddDetRecoParam(Int_t iDet, AliDetectorRecoParam* param);
  Bool_t                        AddDetRecoParamArray(Int_t iDet, TObjArray* parArray);

  AliRecoParam(const AliRecoParam&);
  AliRecoParam& operator=(const AliRecoParam&);

private:

  Int_t      fDetRecoParamsIndex[kNSpecies][kNDetectors]; // index to fDetRecoParams arrays
  TObjArray *fDetRecoParams[kNDetectors];   // array with reconstruction-parameter objects for all detectors
  EventSpecie_t fEventSpecie;               // current event specie

  ClassDef(AliRecoParam, 5)
};


#endif
