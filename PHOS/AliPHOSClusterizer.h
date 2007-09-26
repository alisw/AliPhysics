#ifndef ALIPHOSCLUSTERIZER_H
#define ALIPHOSCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.42  2007/08/28 12:55:07  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.41  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.40  2006/08/25 16:56:30  kharlov
 * Compliance with Effective C++
 *
 * Revision 1.39  2006/03/30 13:04:56  hristov
 * AliRawReader is not persistent
 *
 * Revision 1.38  2005/07/25 15:53:09  kharlov
 * Set raw data reader
 *
 * Revision 1.37  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Base class for the clusterization algorithm 
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// --- ROOT system ---

#include <TObject.h>

class TTree;

class AliPHOSGeometry;

class AliPHOSClusterizer : public TObject {

public:

  AliPHOSClusterizer() ;        // default ctor
  AliPHOSClusterizer(AliPHOSGeometry *geom);
  virtual ~AliPHOSClusterizer() ; // dtor

  virtual void    Digits2Clusters(Option_t *option) = 0;
  virtual Float_t GetEmcClusteringThreshold()const = 0;
  virtual Float_t GetEmcLocalMaxCut()const = 0;
  virtual Float_t GetEmcLogWeight()const = 0;
  virtual Float_t GetEmcTimeGate() const = 0;
  virtual Float_t GetCpvClusteringThreshold()const = 0;
  virtual Float_t GetCpvLocalMaxCut()const = 0;
  virtual Float_t GetCpvLogWeight()const = 0;

  virtual void Print(const Option_t * = "") const = 0;

  virtual void SetEmcClusteringThreshold(Float_t) = 0;
  virtual void SetEmcLocalMaxCut(Float_t )        = 0;
    
  virtual void SetEmcLogWeight(Float_t)           = 0;
  virtual void SetEmcTimeGate(Float_t)            = 0;
  virtual void SetCpvClusteringThreshold(Float_t) = 0;
  virtual void SetCpvLocalMaxCut(Float_t)         = 0;
  virtual void SetCpvLogWeight(Float_t)           = 0;
  virtual void SetUnfolding(Bool_t)               = 0;

  virtual const char * Version() const = 0;

  virtual void SetInput(TTree *digitsTree);
  virtual void SetOutput(TTree *clustersTree);

protected:

  AliPHOSGeometry *fGeom; // Pointer to PHOS geometry
  TClonesArray *fDigitsArr; // Array with input digits
  TTree *fTreeR; // Tree with output clusters
  TObjArray *fEMCRecPoints; // Array with EMC clusters
  TObjArray *fCPVRecPoints; // Array with CPV clusters

private:
  AliPHOSClusterizer(const AliPHOSClusterizer & clusterizer);
  AliPHOSClusterizer & operator = (const AliPHOSClusterizer &clusterer);

  ClassDef(AliPHOSClusterizer,5)  // Clusterization algorithm class 

} ;

#endif // AliPHOSCLUSTERIZER_H
