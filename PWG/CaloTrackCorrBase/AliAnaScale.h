#ifndef ALIANASCALE_H
#define ALIANASCALE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaScale
/// \ingroup CaloTrackCorrelationsBase 
/// \brief Scale input histograms
///
/// Take histograms output of other tasks and scale them.
///
/// \author Yves Schutz <Yves.Schutz@cern.ch>, CNRS-CERN
//_________________________________________________________________________

#include "AliAnalysisTask.h"  
class TH1F ; 

class AliAnaScale : public AliAnalysisTask {

public:
  
  AliAnaScale() ;
  
  AliAnaScale(const char *name) ;
  
  /// Destructor not implemented.
  virtual ~AliAnaScale() { ; }
   
  virtual void ConnectInputData(Option_t * = "");
  
  virtual void CreateOutputObjects(); 
  
  virtual void Init() ; 	
  
  virtual void LocalInit()                { Init()         ; }
  
  virtual void Exec(Option_t * opt = "") ;
  
  void         Set(Double_t val)          { fScale = val   ; }
  
  void         SetDebugLevel(Int_t level) { fDebug = level ; }

  void         MakeSumw2(Bool_t sum)      { fSumw2 = sum   ; }

private:
  
  /// Copy constructor not implemented.
  AliAnaScale(           const AliAnaScale&); 
  
  /// Assignment operator not implemented.
  AliAnaScale& operator=(const AliAnaScale&); 

  Int_t     fDebug ;      ///< Debug flag.
  
  Float_t   fScale ;      ///< Scaling factor. 

  // Histograms

  TList   * fInputList  ; //!<! Input data list.
  
  TList   * fOutputList ; //!<! Output data list.
  
  Bool_t    fSumw2 ;      ///<  Compute sum of squares of weights for bin content error calculation.
  
  TH1F *    fhCount;      //!<! Counter histogram for file merging.

  /// \cond CLASSIMP
  ClassDef(AliAnaScale, 2) ; 
  /// \endcond

};

#endif // ALIANASCALE_H
