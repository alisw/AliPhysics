#ifndef ALIPHOSCLUSTERIZER_H
#define ALIPHOSCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// --- ROOT system ---

#include "TTask.h" 

// --- Standard library ---

// --- AliRoot header files ---

//#include "AliPHOSDigit.h"

class AliPHOSClusterizer : public TTask {

public:

  AliPHOSClusterizer() ;          // ctor
  AliPHOSClusterizer(const char * headerFile,const char * digitsBrancheTitle=0);
  virtual ~AliPHOSClusterizer() ; // dtor

  virtual Float_t GetEmcClusteringThreshold()const = 0 ; 
  virtual Float_t GetEmcLocalMaxCut()const = 0 ; 
  virtual Float_t GetEmcLogWeight()const = 0 ; 
  virtual Float_t GetCpvClusteringThreshold()const = 0 ; 
  virtual Float_t GetCpvLocalMaxCut()const = 0 ; 
  virtual Float_t GetCpvLogWeight()const = 0 ; 
  virtual Float_t GetPpsdClusteringThreshold()const = 0 ; 
  virtual char *  GetRecPointsBranch() const = 0 ;
  virtual char *  GetDigitsBranch() const = 0 ;

  virtual void MakeClusters() = 0 ; 
  virtual void Print(Option_t * option)const = 0;

  virtual void SetEmcClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetEmcLocalMaxCut(Float_t cut) = 0 ; 
  virtual void SetEmcLogWeight(Float_t w) = 0 ; 
  virtual void SetCpvClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetCpvLocalMaxCut(Float_t cut) = 0 ; 
  virtual void SetCpvLogWeight(Float_t w) = 0 ; 
  virtual void SetPpsdClusteringThreshold(Float_t cluth) = 0 ; 

  virtual void SetDigitsBranch(const char * title) = 0 ;
  virtual void SetRecPointsBranch(const char *title) = 0 ;

  virtual void SetUnfolding(Bool_t toUnfold ) = 0 ;

 
  ClassDef(AliPHOSClusterizer,1)  // Clusterization algorithm class 

} ;

#endif // AliPHOSCLUSTERIZER_H
