//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVETPC_H
#define ALIHLTEVETPC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveTPC.h
/// @author Svein Lindal
/// @brief  TPC Instance of Eve display processor

#include "AliHLTEveBase.h"
class TEvePointSetArray;
class TEvePointSet;
class TH1F;

class AliHLTEveTPC : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveTPC();

  /** Destructor **/
 ~AliHLTEveTPC();
  
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:

  /** copy constructor prohibited */
  AliHLTEveTPC(const AliHLTEveTPC&);
  /** assignment operator prohibited */
  AliHLTEveTPC& operator = (const AliHLTEveTPC& );

  /** Create point set for clusters */
  TEvePointSet * CreatePointSet();
  /** Create point set array for colour coded clusters */
  TEvePointSetArray * CreatePointSetArray();

  /** Process clusters block */
  Int_t ProcessClusters( AliHLTHOMERBlockDesc * block, TEvePointSet * cont, TEvePointSetArray * contCol );

  /** Draw the TPC histograms */
  void DrawHistograms();

  TEvePointSet * fEveClusters;          //Clusters pointset
  TEvePointSetArray * fEveColClusters;  //Color coded clusters pointset
  const Int_t fNColorBins;             //Number of colorbins for the colored clusters

  TH1F * fHistCharge;                  //Histo
  TH1F * fHistQMax;                    //Histo
  TH1F * fHistQMaxOverCharge;          //Histo

  ClassDef(AliHLTEveTPC, 0);

};

#endif
