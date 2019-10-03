#ifndef ALIANALYSISNONPRIMARYVERTICES_H
#define ALIANALYSISNONPRIMARYVERTICES_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

#include "AliAnalysisCuts.h"

/// \brief Class to filter out vertices which are not primary ones (can be pileup ones though)
/// \author L. Aphecetche (Subatech)
///

class AliAnalysisNonPrimaryVertices : public AliAnalysisCuts
{
public:
  AliAnalysisNonPrimaryVertices();
  virtual ~AliAnalysisNonPrimaryVertices() {}
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  
  ClassDef(AliAnalysisNonPrimaryVertices,1); // Select primary (and pileup) vertices
};

#endif
