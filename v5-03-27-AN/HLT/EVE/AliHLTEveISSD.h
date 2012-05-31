//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEISSD_H
#define ALIHLTEVEISSD_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveISSD.h
/// @author Svein Lindal
/// @brief  SDD Instance of Eve display processor

#include "AliHLTEveITS.h"
class TEvePointSet;
class TCanvas;

class AliHLTEveISSD : public AliHLTEveITS {

public:
  
  /** Constructor  **/
  AliHLTEveISSD();

  /** Destructor **/
 ~AliHLTEveISSD();

  /** Inherited from AliHLTEveBase */
  virtual void UpdateElements();
  /** Inherited from AliHLTEveBase */
  virtual void ResetElements();


private:

  /** copy constructor prohibited */
  AliHLTEveISSD(const AliHLTEveISSD&);
  /** assignment operator prohibited */
  AliHLTEveISSD& operator = (const AliHLTEveISSD& );
  
  /** Inherited from AliHLTEveITS */
  void SetUpPointSet(TEvePointSet* ps);
  
  /** Inherited from AliHLTEveBase */
  void AddHistogramsToCanvas(AliHLTHOMERBlockDesc* block, TCanvas *canvas, Int_t &cdCount );
  
  TCanvas * f2DCanvas;  //Canvas containing 2D QA histograms

  Int_t f2DHistoCount;  //Counter tracking where to draw histogram

  ClassDef(AliHLTEveISSD, 0);

};

#endif
