// -*- mode: c++ -*- 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    FMDsimLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:19:24 2006
    @brief   Link specifications fo simulation library
*/
/* $Id$ */
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// #pragma link C++ class  std::pair<Float_t,UShort_t>;
// #pragma link C++ class  AliFMDMap<std::pair<Float_t,UShort_t> >;
// #pragma link C++ typedef  AliFMDEdepMap;
#pragma link C++ class  AliFMDEdepMap+;
#pragma link C++ class  AliFMDEdepHitPair+;
#pragma link C++ class  AliFMDHit+;
#pragma link C++ class  AliFMD+;
#pragma link C++ class  AliFMDv0+;
#pragma link C++ class  AliFMDv1+;
#pragma link C++ class  AliFMDBaseDigitizer+;
#pragma link C++ class  AliFMDDigitizer+;
#pragma link C++ class  AliFMDSDigitizer+;
#pragma link C++ class  AliFMDSSDigitizer+;
#pragma link C++ class  AliFMDRawWriter+;
#pragma link C++ class  AliFMDQADataMakerSim+;

#else
# error Not for compilation 

#endif
//
// EOF
//
