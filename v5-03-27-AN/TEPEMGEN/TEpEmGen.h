#ifndef ROOT_TEpEmGen
#define ROOT_TEpEmGen
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//------------------------------------------------------------------------
// TEpEmGen is an interface class to fortran event generator of
// single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 GeV/c
//%
// Yuri.Kharlov@cern.ch
// 9 October 2002
//------------------------------------------------------------------------

#include "TGenerator.h"

// c++ interface to the f77 program - event generator of
// e+e- pair production in ultraperipheral ion collisions
// Author: Yuri Kharlov, 20 September 2002

class TEpEmGen : public TGenerator {

public:
  TEpEmGen();
  virtual ~TEpEmGen();

  void Initialize    (Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmax);
#if ROOT_VERSION_CODE >= 330500 // Root 5.11/04
  void GenerateEvent() {TGenerator::GenerateEvent();};
#endif
  void GenerateEvent (Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmax,
	       	      Double_t &yElectron, Double_t &yPositron,
		      Double_t &xElectron, Double_t &xPositron,
		      Double_t &phi12,     Double_t &weight);
  Double_t GetXsection();
  Double_t GetDsection();

  ClassDef(TEpEmGen,1)  //Interface to EpEmGen Event Generator
};

#endif
