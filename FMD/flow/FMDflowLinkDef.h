// -*- mode: c++ -*- 
#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/** @file    FMDbaseLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:46 2006
    @brief   Link specifications for base library 
*/
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliFMDFlowAxis+;
#pragma link C++ class AliFMDFlowBin+;
#pragma link C++ class AliFMDFlowBin;+;
#pragma link C++ class AliFMDFlowBinned1D+;
#pragma link C++ class AliFMDFlowBin;+;
#pragma link C++ class AliFMDFlowBinned2D+;
#pragma link C++ class AliFMDFlowEventPlane+;
#pragma link C++ class AliFMDFlowHarmonic+;
#pragma link C++ class AliFMDFlowResolution+;
#pragma link C++ class AliFMDFlowResolutionStar+;
#pragma link C++ class AliFMDFlowResolutionTDR+;
#pragma link C++ class AliFMDFlowStat+;
#pragma link C++ class AliFMDFlowTrueBin+;
#pragma link C++ class AliFMDFlowTrue1D+;
#pragma link C++ function  NormalizeAngle(Double_t);
  

#pragma link C++ namespace AliFMDFlowBessel;
#pragma link C++ function  AliFMDFlowBessel::I01(Double_t,Double_t&,Double_t&,Double_t&,Double_t&)
#pragma link C++ function  AliFMDFlowBessel::Ihalf(Int_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::Iwhole(Int_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::Inu(Double_t,Double_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::I(Double_t,Double_t);
#pragma link C++ function  AliFMDFlowBessel::DiffI(Double_t,Double_t);
  

#else
# error Not for compilation 
#endif
//
// EOF
//
