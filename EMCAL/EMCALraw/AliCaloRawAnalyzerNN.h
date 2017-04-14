#ifndef ALICALORAWANALYZERNN_H
#define ALICALORAWANALYZERNN_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerNN
/// \ingroup EMCALraw
/// \brief Raw data fitting: Neural network
///
/// Evaluation of peak position
/// and amplitude using Neural Networks (NN)
///
/// \author Paola La Rocca (Catania) 
//_________________________________________________________________________

#include "AliCaloRawAnalyzer.h"

class AliCaloBunchInfo;
class AliCaloFitResults;
class AliCaloNeuralFit;

class  AliCaloRawAnalyzerNN : public AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory; // self explanatory
  
public:
  
  virtual ~AliCaloRawAnalyzerNN();
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
                                     UInt_t altrocfg1, UInt_t altrocfg2 );
  
private:
  
  AliCaloRawAnalyzerNN();
  AliCaloRawAnalyzerNN(                 const AliCaloRawAnalyzerNN & );
  AliCaloRawAnalyzerNN   & operator = ( const AliCaloRawAnalyzerNN & );
  
  AliCaloNeuralFit * fNeuralNet;  ///< Pointer to the class whick actually implements the Neural Network for EMCAL
  Double_t           fNNInput[5]; ///< The 5 input Neurons to the network ( mix bin + to samples on each side )
  
  /// \cond CLASSIMP
  ClassDef( AliCaloRawAnalyzerNN, 1 ) ;
  /// \endcond
  
};

#endif //ALICALORAWANALYZERNN_H
