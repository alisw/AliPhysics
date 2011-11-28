//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVETRD_H
#define ALIHLTEVETRD_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveTRD.h
/// @author Svein Lindal
/// @brief  TRD Instance of Eve display processor

#include "AliHLTEveBase.h"
class TEvePointSetArray;
class TEvePointSet;
class TH1F;
class TClonesArray;

class AliHLTEveTRD : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveTRD();

  /** Destructor **/
 ~AliHLTEveTRD();
  
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:

  /** copy constructor prohibited */
  AliHLTEveTRD(const AliHLTEveTRD&);
  /** assignment operator prohibited */
  AliHLTEveTRD& operator = (const AliHLTEveTRD& );

  /** Create clusters pointset */
  TEvePointSet * CreatePointSet();
  /** Create point set array for colour coded clusters */
  TEvePointSetArray * CreatePointSetArray();

  /** Proces clusters block */
  Int_t ProcessClusters( AliHLTHOMERBlockDesc * block, TEvePointSetArray * contCol );

  /** Inherited from AliHLTEveBase */
  void AddHistogramsToCanvas(AliHLTHOMERBlockDesc* block, TCanvas * canvas, Int_t &cdCount );

  TEvePointSet * fEveClusters;         //clusters pointset
  TEvePointSetArray * fEveColClusters; //Color coded clusters pointset
  const Int_t fNColorBins;    //Number of colorbins for the colored clusters
  TClonesArray* fClusterArray;

  ClassDef(AliHLTEveTRD, 0);
};

#endif
