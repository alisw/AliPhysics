// -*- mode: c++ -*-
#ifndef ALIEMCALRAWUTILS_H
#define ALIEMCALRAWUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliEMCALRawUtils
/// \brief Handling of raw data.
///
///  Utility Class for handling Raw data
///  Does all transitions from Digits to Raw and vice versa,
///  for simu and reconstruction
///
///  Only one raw signal per digit is generated;
///  either high-gain or low-gain
///  No pedestal is added to the raw signal.
///
/// \author Marco van Leeuwen <Marco.Van.Leeuwen@cern.ch>, LBL. First implementation.
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. Major refactoring.
/// \author David Silvermyr <David.Silvermyr@cern.ch>, Oak Ridge. Trimming and real data adjustments.
//_________________________________________________________________________

#include "TObject.h" // for ROOT types
#include <TString.h>

class AliAnalysisManager;
class AliCaloRawStreamV3;
class AliAltroMapping;
class TGraph;
class AliRawReader;
class AliEMCALGeometry;
class AliCaloCalibPedestal;
class AliCaloRawAnalyzer;
class AliEMCALTriggerRawDigitMaker;
class AliEMCALTriggerData;
class TClonesArray;

#include "AliCaloConstants.h"

class AliEMCALRawUtils : public TObject {
 
public:
  
  AliEMCALRawUtils(Algo::fitAlgorithm fitAlgo = Algo::kStandard);
  
  virtual ~AliEMCALRawUtils();
  
  void     Digits2Raw();
  
  void     Raw2Digits(AliRawReader *reader, TClonesArray *digitsArr, const AliCaloCalibPedestal* pedbadmap,
				      TClonesArray *digitsTRG=0x0, TClonesArray *trgData=0x0);
  
  void     AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Float_t amp, Float_t time, Float_t chi2, Int_t ndf);
  
  void     TrimDigits(TClonesArray *digitsArr);
  
  Int_t    GetNoiseThreshold()          const { return fNoiseThreshold    ; }
    
  Int_t    GetNPedSamples()             const { return fNPedSamples       ; }
    
  Bool_t   GetRemoveBadChannels()       const { return fRemoveBadChannels ; }
    
  Int_t    GetFittingAlgorithm()        const { return fFittingAlgorithm  ; }
    
  Float_t  GetTimeMax()                 const { return fTimeMax           ; }
    
  Float_t  GetTimeMin()                 const { return fTimeMin           ; }
    
  Bool_t   UseFALTRO()                  const { return fUseFALTRO         ; }
  
  void     SetNoiseThreshold(Int_t val)       { fNoiseThreshold=val       ; }
    
  void     SetNPedSamples(Int_t val)          { fNPedSamples=val          ; }
    
  void     SetRemoveBadChannels(Bool_t val)   { fRemoveBadChannels=val    ; }
    
  void     SetFittingAlgorithm(Int_t val) ;
    
  void     SetTimeMin(Float_t t)              { fTimeMin   = t            ; }
    
  void     SetTimeMax(Float_t t)              { fTimeMax   = t            ; }
    
  void     SetFALTROUsage(Bool_t val)         { fUseFALTRO = val          ; }
  
  void     SetL1PhaseUsage(Bool_t val)        { fUseL1Phase = val         ; }
  
  AliCaloRawAnalyzer *GetRawAnalyzer()  const { return fRawAnalyzer       ; }
  
  virtual Option_t* GetOption()         const { return fOption.Data()     ; }
    
  void    SetOption(const Option_t* opt)      { fOption = opt             ; }
  
private:
    
  AliEMCALRawUtils            (const AliEMCALRawUtils& rawUtils);
    
  AliEMCALRawUtils& operator =(const AliEMCALRawUtils& rawUtils);
  
  Int_t   fNoiseThreshold;              ///< Threshold to consider signal or noise.
    
  Int_t   fNPedSamples;                 ///< Number of samples to use in pedestal calculation.
    
  AliEMCALGeometry* fGeom;              ///< Geometry.
    
  AliAltroMapping*  fMapping[4];        ///< What is the array size?
    
  TString fOption;                      ///< Option passed from Reconstructor.
    
  Bool_t  fRemoveBadChannels;           ///< Select if bad channels are removed before fitting.
    
  Int_t   fFittingAlgorithm;            ///< Select the fitting algorithm.
    
  Float_t fTimeMin;                     ///< Minimum threshold for the time of the signal.
    
  Float_t fTimeMax;                     ///< Maximum threshold for the time of the signal.
    
  Bool_t  fUseFALTRO;                   ///< Use FALTRO and pass it to the digits.
  
  Bool_t  fUseL1Phase;                  ///< Use L1Phase time shift.
    
  AliCaloRawAnalyzer *fRawAnalyzer;     ///< e.g. for sample selection for fits.
    
  AliEMCALTriggerRawDigitMaker* fTriggerRawDigitMaker;	///< Trigger raw digit info.
 
  /// \cond CLASSIMP
  ClassDef(AliEMCALRawUtils,8) ;
  /// \endcond

};

#endif
