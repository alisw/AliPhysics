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
 
#pragma link C++ class  AliFMDIndex+;
#pragma link C++ class  AliFMDObjIndex+;
#pragma link C++ class  AliFMDBaseDigit+;
#pragma link C++ class  AliFMDDigit+;
#pragma link C++ class  AliFMDSDigit+;
#pragma link C++ class  AliFMDBoolMap+;
#pragma link C++ class  AliFMDUShortMap+;
#pragma link C++ class  AliFMD1+;
#pragma link C++ class  AliFMD2+;
#pragma link C++ class  AliFMD3+;
#pragma link C++ class  AliFMDRing+;
#pragma link C++ class  AliFMDDetector+;
#pragma link C++ class  AliFMDGeometry+;
#pragma link C++ class  AliFMDParameters+;
#pragma link C++ class  AliFMDCalibPedestal+;
#pragma link C++ class  AliFMDCalibGain+;
#pragma link C++ class  AliFMDCalibSampleRate+;
#pragma link C++ class  AliFMDCalibStripRange+;
#pragma link C++ class  AliFMDAltroMapping+;
#pragma link C++ class  AliFMDPreprocessor+;
#pragma link C++ class  AliFMDQAChecker+;
#pragma link C++ class  AliFMDGeometryBuilder+;
#pragma link C++ class  AliFMDSurveyToAlignObjs+;

// #pragma link C++ class  AliFMDAltroIO+;
// #pragma link C++ class  AliFMDAltroReader+;
// #pragma link C++ class  AliFMDAltroWriter+;

#else
# error Not for compilation 
#endif
//
// EOF
//
