#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: ITSbaseLinkDef.h 36856 2009-11-16 16:17:07Z masera $ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// ITS upgrade classes 
// v0 part >>>
#pragma link C++ class  AliITSDigitUpgrade+;
#pragma link C++ class  AliITSsegmentationUpgrade+;
#pragma link C++ class  AliITSRecPointU+;
// v0 part <<<
//
#pragma link C++ class  AliITSULoader+;
#pragma link C++ class  AliITSUGeomTGeo+;
#pragma link C++ class  AliITSUCalibrationPix+;
#pragma link C++ class  AliITSUSegmentationPix+;
#pragma link C++ class  AliITSUSensMap+;
#pragma link C++ class  AliITSUSDigit+;
#pragma link C++ class  AliITSUParamList+;
#pragma link C++ namespace  AliITSUAux;
#endif
