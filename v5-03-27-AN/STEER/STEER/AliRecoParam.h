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
class THashTable;

class AliRecoParam : public TObject
{

 public: 
  AliRecoParam();
  virtual ~AliRecoParam();  
  //
  enum {
    kNSpecies = 5,   // number of event species
// #ifdef MFT_UPGRADE  
//     kNDetectors = 17 // number of detectors (last one is MFT!)
// #else
//     kNDetectors = 16 // number of detectors (last one is GRP!)
// #endif
    kNDetectors = 17 // number of detectors (last one is MFT!)   // AU
  };
  enum EventSpecie_t {kDefault = 1,
		      kLowMult = 2,
		      kHighMult = 4,
		      kCosmic = 8,
		      kCalib = 16};

  static Int_t                  AConvert(EventSpecie_t es) ; 
  static EventSpecie_t          Convert(Int_t ies) ; 
  static EventSpecie_t          ConvertIndex(Int_t index) ;

  virtual void                  Print(Option_t *option="") const;
  const TObjArray              *GetDetRecoParamArray(Int_t iDet) const { return fDetRecoParams[iDet]; }
  void                          SetEventSpecie(const AliRunInfo*runInfo, const AliEventInfo &evInfo,
					       const THashTable*cosmicTriggersList);
  EventSpecie_t                 GetEventSpecie() const { return fEventSpecie; }
  static const char*            GetEventSpecieName(EventSpecie_t es);
  static const char*            GetEventSpecieName(Int_t esIndex);
  const char*                   PrintEventSpecie() const;
  const AliDetectorRecoParam   *GetDetRecoParam(Int_t iDet) const;
  void                          AddDetRecoParam(Int_t iDet, AliDetectorRecoParam* param);
  Bool_t                        AddDetRecoParamArray(Int_t iDet, TObjArray* parArray);

  AliRecoParam(const AliRecoParam&);
  AliRecoParam& operator=(const AliRecoParam&);


private:

  Int_t      fDetRecoParamsIndex[kNSpecies][kNDetectors]; // index to fDetRecoParams arrays
  TObjArray *fDetRecoParams[kNDetectors];   // array with reconstruction-parameter objects for all detectors
  EventSpecie_t fEventSpecie;               // current event specie
  static TString fkgEventSpecieName[] ; // the names of the event species
  ClassDef(AliRecoParam, 6)
};


#endif
