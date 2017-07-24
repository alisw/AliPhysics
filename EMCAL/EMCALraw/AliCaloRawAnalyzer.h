// -*- mode: c++ -*-
#ifndef ALICALORAWANALYZER_H
#define ALICALORAWANALYZER_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzer
/// \ingroup EMCALraw
/// \brief  Base class for extraction of signal amplitude and peak position
///
/// Base class for extraction 
/// of signal amplitude and peak position
/// from CALO Calorimeter RAW data.
/// It contains some utilities for preparing / selecting
/// signals suitable for signal extraction
/// by derived classes
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________

// standard library
#include "Rtypes.h"
#include <vector>

// ROOT system
#include "TObject.h"
#include "TObjArray.h"

// Calo system
#include "AliCaloFitResults.h"
#include "AliCaloConstants.h"
using namespace ALTRO;
using namespace CALO;

class AliCaloBunchInfo;

class  AliCaloRawAnalyzer : public TObject
{
  
public:
  
  AliCaloRawAnalyzer(const char *name="AliCaloRawAnalyzer", const char *nameshort="RawAna");
  virtual ~AliCaloRawAnalyzer() { ; }
  
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &/*bunchvector*/, 
                                     UInt_t /*altrocfg1*/,  UInt_t /*altrocfg2*/ )  = 0;
  
  static void PrintBunches( const std::vector<AliCaloBunchInfo> &bunchvector );
  static void PrintBunch  ( const AliCaloBunchInfo &bunch );
  
  int      PreFitEvaluateSamples( const std::vector<AliCaloBunchInfo>  &bunchvector, 
                                 UInt_t altrocfg1, UInt_t altrocfg2, Int_t & index,
                                 Float_t & maxf, short & maxamp, short & maxampindex,
                                 Float_t & ped, int & first, int & last, int acut);
  
  void     SetTimeConstraint  (int min, int max );
  void     SetVerbose         (bool verbose = true){ fVerbose = verbose      ; }
  void     SetIsZeroSuppressed(bool iszs = true)   { fIsZerosupressed = iszs ; }
  
  void     SetAmpCut     (Float_t cut)    { fAmpCut      = cut  ; }
  void     SetFitArrayCut(Int_t cut)      { fFitArrayCut = cut  ; }
  void     SetNsampleCut (Int_t cut)      { fNsampleCut  = cut  ; }
  void     SetOverflowCut(Int_t cut)      { fOverflowCut = cut  ; }
  void     SetNsamplePed (Int_t i)        { fNsamplePed  = i    ; }
  void     SetL1Phase    (Double_t phase) { fL1Phase     = phase ; }
  
  bool     GetIsZeroSuppressed()    const { return fIsZerosupressed ; }
  Float_t  GetAmpCut()              const { return fAmpCut      ; }
  Int_t    GetFitArrayCut()         const { return fFitArrayCut ; }
  Int_t    GetNsampleCut()          const { return fNsampleCut  ; }
  Int_t    GetOverflowCut()         const { return fOverflowCut ; }
  Int_t    GetNsamplePed()          const { return fNsamplePed  ; }
  
  // access to array info
  Double_t GetReversed(const int i) const { return fReversed[i] ; }
  const char * GetAlgoName()        const { return fName        ; }
  const char * GetAlgoAbbr()        const { return fNameShort   ; }
  Algo::fitAlgorithm GetAlgo()      const { return fAlgo        ; }
  
  // Used in AliCaloRawAnalyzerFitter
  //
  Float_t  GetTau()                 const { return fTau         ; }
  void     SetTau   (Float_t tau)         { fTau = tau          ; }  
  Bool_t   GetFixTau()      const         { return fFixTau      ; }
  void     SetFixTau(Bool_t b)            { fFixTau = b         ; }
  //
  
  Double_t CalculateChi2(const Double_t amp, const Double_t time,
                         const Int_t first, const Int_t last,
                         const Double_t adcErr=1,
                         const Double_t tau=2.35) const;
  
  void     CalculateMeanAndRMS(const Int_t first, const Int_t last,
                               Double_t & mean, Double_t & rms);
  
  short    Max( const AliCaloBunchInfo *const bunch, int * maxindex ) const;
  
  UShort_t Max( const UShort_t *data, const int length ) const;
  
  bool     CheckBunchEdgesForMax( const AliCaloBunchInfo *const bunch ) const;
  
  bool     IsInTimeRange( const int maxindex, const int maxtime, const int mintime ) const;
  
  Float_t  ReverseAndSubtractPed( const AliCaloBunchInfo *bunch,
                                 UInt_t altrocfg1,  UInt_t altrocfg2,
                                 double *  outarray ) const;
  
  int      SelectBunch( const std::vector<AliCaloBunchInfo> &bunchvector,
                       short * maxampbin, short * maxamplitude );
  
  void     SelectSubarray( const Double_t *date, int length, short maxindex,
                          int * first, int * last, int cut ) const;
  
  Float_t  EvaluatePedestal( const UShort_t * const data, const int length ) const;
  
protected:
  
  Double_t fReversed[ALTROMAXSAMPLES]; ///< Reversed sequence of samples (pedestalsubtracted)
  
  int      fMinTimeIndex;              ///< The timebin of the max signal value must be between fMinTimeIndex and fMaxTimeIndex
  int      fMaxTimeIndex;              ///< The timebin of the max signal value must be between fMinTimeIndex and fMaxTimeIndex
  
  int      fFitArrayCut;               ///< Cut on ADC value (after ped. subtraction) for signals used for fit
  
  Float_t  fAmpCut;                    ///< Max ADC - pedestal must be higher than this befor attemting to extract the amplitude 
  
  int      fNsampleCut;                ///< Minimum number of sample require before attemting to extract signal parameters 
  int      fOverflowCut;               ///< Value when ADC starts to saturate
  
  int      fNsamplePed;                ///< Number of samples used for pedestal calculation (first in bunch) 
  
  bool     fIsZerosupressed;           ///< Wether or not the data is zeros supressed, by default its assumed that the baseline is also subtracted if set to true
  
  bool     fVerbose;                   ///< Print debug information to std out if set to true
  
  char     fName[256];                 ///< Name of the algorithm
  char     fNameShort[256];            ///< Abbrevation for the name
  
  Algo::fitAlgorithm fAlgo;            ///< Which algorithm to use
  
  Double_t fL1Phase;                   ///< Phase of the ADC sampling clock relative to the LHC clock
  
  Double_t fAmp;                       ///< The amplitude in entities of ADC counts
  
  Double_t fTof;                       ///< The amplitude in entities of ADC counts
  
  Float_t  fTau;                       ///< Rise time of the signal (peak position = t0 +tau), by defauly it is 235 ns
  
  Bool_t   fFixTau;                    ///< Fixed fit parameter or not, used in AliCaloRawAnalyzerFitter
  
  /// \cond CLASSIMP
  ClassDef(AliCaloRawAnalyzer, 3) ;
  /// \endcond
  
};

#endif //ALICALORAWANALYZER_H
