#ifndef ALICALORAWANALYZERNN_H
#define ALICALORAWANALYZERNN_H

/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <perthomas.hille@yale.edu>            *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Evaluation of peak position
// and amplitude using Neural Networks (NN)
// ------------------

#include "AliCaloRawAnalyzer.h"

class AliCaloBunchInfo;
class AliCaloFitResults;
class AliCaloNeuralFit;


class  AliCaloRawAnalyzerNN : public AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory; // self explanatory
 public:
  // AliCaloRawAnalyzerNN();
  virtual ~AliCaloRawAnalyzerNN();
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
				       const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  //  virtual void SelectSubarray( const Double_t *fData, const int length, const short maxindex, int *const  first, int *const last ) const;

 private:
  AliCaloRawAnalyzerNN();
  AliCaloRawAnalyzerNN( const AliCaloRawAnalyzerNN   & );
  AliCaloRawAnalyzerNN   & operator = ( const  AliCaloRawAnalyzerNN  & );
  AliCaloNeuralFit *fNeuralNet; // pointer to the class whick actually implements the Neural Network for EMCAL
  Double_t fNNInput[5]; // The 5 input Neurons to the network ( mix bin + to samples on each side )
  ClassDef( AliCaloRawAnalyzerNN, 1 )

};

#endif
