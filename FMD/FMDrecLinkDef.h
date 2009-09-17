// -*- mode: c++ -*- 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice
 */
/* $Id$ */
/** @file    FMDrecLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:19:08 2006
    @brief   Link specifications for reconstruction library
    
*/
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
// #pragma link C++ class  AliFMDMap<UShort_t>;
// #pragma link C++ typedef AliFMDAdcMap;
#pragma link C++ class  AliFMDReconstructor+;
#pragma link C++ class  AliFMDRecoParam+;
#pragma link C++ class  AliFMDRecPoint+;
#pragma link C++ class  AliFMDRawStream+;
#pragma link C++ class  AliFMDRawReader+;
#pragma link C++ class  AliFMDQADataMakerRec+;
#pragma link C++ class  AliFMDOfflineTrigger+;

#else
# error Not for compilation 
#endif
//
// EOF
//
