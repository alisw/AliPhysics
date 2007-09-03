#ifndef ALIEMCALCLUSTERIZER_H
#define ALIEMCALCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
// --- ROOT system ---

#include "TObject.h" 
class TTree;

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALClusterizer : public TObject {

public:

  AliEMCALClusterizer() ;        // default ctor
  virtual ~AliEMCALClusterizer() ; // dtorEM

  virtual void    Digits2Clusters(Option_t *option) = 0;

  virtual Float_t GetTimeCut() const = 0;

  virtual void SetECAClusteringThreshold(Float_t) = 0;
  virtual void SetECALocalMaxCut(Float_t)         = 0;
  virtual void SetECALogWeight(Float_t)           = 0;
  virtual void SetTimeCut(Float_t)                = 0;
  virtual void SetUnfolding(Bool_t)               = 0;

  virtual const char * Version() const {Warning("Version", "Not Defined") ; return 0 ; } 

  virtual void SetInput(TTree *digitsTree);
  virtual void SetOutput(TTree *clustersTree);

protected:

  virtual void MakeClusters(char* opt) = 0;

  TClonesArray *fDigitsArr; // Array with EMCAL digits
  TTree *fTreeR;            // Tree with output clusters
  TObjArray    *fRecPoints; // Array with EMCAL clusters

private:
  AliEMCALClusterizer(const AliEMCALClusterizer &); //copy ctor
  AliEMCALClusterizer & operator = (const AliEMCALClusterizer &);

  ClassDef(AliEMCALClusterizer,1)  // Clusterization algorithm class 
} ;

#endif // AliEMCALCLUSTERIZER_H
