#ifndef ALICALORAWANALYZERCOMPARISON_H
#define ALICALORAWANALYZERCOMPARISON_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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
  virtual ~AliCaloRawAnalyzerComparison();
  void Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
		 const UInt_t altrocfg1,  const UInt_t altrocfg2, const int event, const int col, const int row ); 
  
  void EventChanged();

  void WriteHistograms();

 private:
  AliCaloRawAnalyzerComparison( const AliCaloRawAnalyzerComparison   & );
  AliCaloRawAnalyzerComparison   & operator = ( const  AliCaloRawAnalyzerComparison  & );

  void IntiHistograms( std::vector <AliCaloRawAnalyzer*> analyzers, AliCaloRawAnalyzer* ref );

  TH1D *fAmpHistograms[NANALYZERS][NZCOLSSMOD][NXROWSSMOD];

  TH2D *fAmplitudeVsEvent[NANALYZERS];  // Amplitude vs envent number
  TH2D *fTofVsEvent[NANALYZERS];        // Tof vs event number
  TH2D *fRefAmpVsAnalyzers[NANALYZERS]; // Amplidue from give analyzer vs reference
  TH2D *fRefTofVsAnalyzers[NANALYZERS]; // Amplidue from give analyzer vs reference
  TH1D *fAmpDiff[NANALYZERS]; // Difference in amplitude between reference
  TH1D *fTofDiff[NANALYZERS]; // Difference in tof between reference
  TH1D *fTofResDifferential[NANALYZERS]; //differntial tof resolution 
  TH1D *fTofResAbsolute[NANALYZERS]; //differntial tof resolution 
  

  std::vector <AliCaloRawAnalyzer*> fRawAnalyzers; // raw analyzers
  AliCaloRawAnalyzer *fReferenceAnalyzer; // reference analyzer
  
  int fMod; // SuperModule index
  int fMonCol1; // column index, for tower 1
  int fMonRow1; // row index, for tower 1
  int fMonCol2; // column index, for tower 2
  int fMonRow2; // row index, for tower 1

  AliCaloFitResults fMon1[NANALYZERS]; // results for tower 1
  AliCaloFitResults fMon2[NANALYZERS]; // results for tower 2

};

#endif
