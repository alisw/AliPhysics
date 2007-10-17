#ifndef ALIANASCALE_H
#define ALIANASCALE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An post event loop analysis task that scales the input histograms 
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class TH1D ; 
class TH1I ; 

class AliAnaScale : public AliAnalysisTask {

public:
  AliAnaScale() ;
  AliAnaScale(const char *name) ;
  virtual ~AliAnaScale() ;
   
  virtual void ConnectInputData(Option_t * = "");
  virtual void CreateOutputObjects(); 
  virtual void Init() ; 	
  virtual void LocalInit() { Init() ; }
  virtual void Exec(Option_t * opt = "") ;
  void Set(const Double_t val) { fScale = val ; }
  void SetDebugLevel(Int_t level) { fDebug = level ; }
  virtual void Terminate(Option_t * opt = "") ;

private:
  // input and output
  Int_t          fDebug ;         // Debug flag
  // task parameters
  Float_t   fScale ;  // Scaling factor 

  // Histograms
  TList   * fInputList ;  //! input data list
  TList   * fOutputList ; //! output data list
   
  ClassDef(AliAnaScale, 1); // a post event loop scaling 
};
#endif // ALIANASCALE_H
