/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id $

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// segmentation (change in progress)
#pragma link C++ class AliMUONSegFactory+;
#pragma link C++ class AliMUONSt12QuadrantSegmentation+; 
#pragma link C++ class AliMUONSt345SlatSegmentation+;
#pragma link C++ class AliMUONTriggerSegmentation+;

// geometry 
#pragma link C++ class AliMUONMathieson+; 

// info classes 
#pragma link C++ class AliMUONConstants+; 
#pragma link C++ class AliMUONLogger+;

// containers
#pragma link C++ class AliMUONData+; 
#pragma link C++ class AliMUONDataIterator+; 
#pragma link C++ class AliMUONDataDigitIterator+; 
#pragma link C++ class AliMUONRawCluster+;
#pragma link C++ class AliMUONDigit+; 
#pragma link C++ class AliMUONHitMapA1+; 

// raw data
#pragma link C++ class AliMUONDigitMaker+;
#pragma link C++ class AliMUONRawWriter+;

// calibration access
#pragma link C++ class AliMUONCalibrationData+;
#pragma link C++ class AliMUONCDB+;
#pragma link C++ class AliMUONHVNamer+;

#endif


