#ifndef ALICALORAWANALYZERCOMPARISON_H
#define ALICALORAWANALYZERCOMPARISON_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerComparison
/// \ingroup EMCAL_TestData
/// \brief  Raw data fitting comparison
///
/// Compare the different fitters
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________


#define NANALYZERS 5

#include <vector>
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"

#define NZCOLSSMOD   48     
#define NXROWSSMOD   24    

class AliCaloRawAnalyzer;
class TH2D;
class TH1D;

class  AliCaloRawAnalyzerComparison
{
 public:
  
  AliCaloRawAnalyzerComparison();
  virtual ~AliCaloRawAnalyzerComparison() {;}
  
  void Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
                 const UInt_t altrocfg1,  const UInt_t altrocfg2, const int event, const int col, const int row );
  
  void EventChanged();

  void WriteHistograms();

 private:
  
  AliCaloRawAnalyzerComparison                ( const AliCaloRawAnalyzerComparison  & );
  AliCaloRawAnalyzerComparison   & operator = ( const AliCaloRawAnalyzerComparison  & );

  void InitHistograms( std::vector <AliCaloRawAnalyzer*> analyzers, AliCaloRawAnalyzer* ref );

  TH1D *fAmpHistograms[NANALYZERS][NZCOLSSMOD][NXROWSSMOD]; ///< amplitude histos

  TH2D *fAmplitudeVsEvent  [NANALYZERS]; ///< Amplitude vs envent number
  TH2D *fTofVsEvent        [NANALYZERS]; ///< Tof vs event number
  TH2D *fRefAmpVsAnalyzers [NANALYZERS]; ///< Amplidue from give analyzer vs reference
  TH2D *fRefTofVsAnalyzers [NANALYZERS]; ///< Amplidue from give analyzer vs reference
  TH1D *fAmpDiff           [NANALYZERS]; ///< Difference in amplitude between reference
  TH1D *fTofDiff           [NANALYZERS]; ///< Difference in tof between reference
  TH1D *fTofResDifferential[NANALYZERS]; ///< Differntial tof resolution 
  TH1D *fTofResAbsolute    [NANALYZERS]; ///< Differntial tof resolution

  std::vector <AliCaloRawAnalyzer*> fRawAnalyzers; ///< Raw analyzers
  AliCaloRawAnalyzer *fReferenceAnalyzer;          ///< Reference analyzer
  
  int fMod;     ///< SuperModule index
  int fMonCol1; ///< column index, for tower 1
  int fMonRow1; ///< row index, for tower 1
  int fMonCol2; ///< column index, for tower 2
  int fMonRow2; ///< row index, for tower 1

  AliCaloFitResults fMon1[NANALYZERS]; ///< results for tower 1
  AliCaloFitResults fMon2[NANALYZERS]; ///< results for tower 2

};

#endif
