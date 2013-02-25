// -*- mode: C++ -*- 
#ifndef ALIESDUTILS_H
#define ALIESDUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//   AliESDUtils - This is a namespace that temporary provides general 
//                 purpose ESD utilities to avoid unnecessary dependencies
//                 between PWG libraries. To be removed/replaced by AODB
//                 framework.
//      
//-------------------------------------------------------------------------
// Author: Andrei Gheata

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif
class AliVEvent;
class AliESDEvent;
class AliVertexerTracks;

namespace AliESDUtils {

  /********************************/
  /* Centrality-related corrections */
  /********************************/

  Float_t GetCorrV0(const AliVEvent* esd, Float_t &v0CorrResc, Float_t *v0multChCorr = NULL, Float_t *v0multChCorrResc = NULL);
  Float_t GetCorrSPD2(Float_t spd2raw,Float_t zv);
  Bool_t  RefitESDVertexTracks(AliESDEvent* esdEv, Int_t algo=6, const Double_t* cuts=0);
  Float_t GetCorrV0A(Float_t v0araw,Float_t zv);
  Float_t GetCorrV0C(Float_t v0craw,Float_t zv);
}  

#endif
