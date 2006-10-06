#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUON+; 
#pragma link C++ class AliMUONv1+; 

// segmentation (change in progress)
#pragma link C++ class AliMUONSegFactory+;
#pragma link C++ class AliMUONSt12QuadrantSegmentation+; 
#pragma link C++ class AliMUONSt345SlatSegmentation+;
#pragma link C++ class AliMUONSt345SlatSegmentationV2+;
#pragma link C++ class AliMUONTriggerSegmentation+;
#pragma link C++ class AliMUONTriggerSegmentationV2+;

// geometry 
#pragma link C++ class AliMUONMathieson+; 

// info classes 
#pragma link C++ class AliMUONConstants+; 
#pragma link C++ class AliMUONTriggerConstants+; 
#pragma link C++ class AliMUONDataInterface+; 
#pragma link C++ class AliMUONLoader+; 
#pragma link C++ class AliMUONChamber+; 
#pragma link C++ class AliMUONChamberTrigger+; 
#pragma link C++ class AliMUONTriggerCircuitNew+; 
#pragma link C++ class AliMUONTriggerCrateStore+; 

// containers
#pragma link C++ class AliMUONData+; 
#pragma link C++ class AliMUONDataIterator+; 
#pragma link C++ class AliMUONVDataIterator+; 
#pragma link C++ class AliMUONDataDigitIterator+; 
#pragma link C++ class AliMUONPoints+; 
#pragma link C++ class AliMUONHit+; 
#pragma link C++ class AliMUONRawCluster+;
#pragma link C++ class AliMUONDigit+; 
#pragma link C++ class AliMUONTransientDigit+;
#pragma link C++ class AliMUONGlobalTrigger+; 
#pragma link C++ class AliMUONRegionalTrigger+; 
#pragma link C++ class AliMUONLocalTrigger+; 
#pragma link C++ class AliMUONTriggerLut+; 

// calibration
#pragma link C++ class AliMUONV2DStore+;
#pragma link C++ class AliMUONV1DStore+;
#pragma link C++ class AliMUON2DMap+;
#pragma link C++ class AliMUON2DMapIterator+;
#pragma link C++ class AliMUON1DArray+;
#pragma link C++ class AliMUONVCalibParam+;
#pragma link C++ class AliMUONCalibParam1I+;
#pragma link C++ class AliMUONCalibParam2F+;
#pragma link C++ class AliMUONCalibrationData+;
#pragma link C++ class AliMUONTriggerEfficiencyCells+;

// raw data
#pragma link C++ class AliMUONDigitMaker+;
#pragma link C++ class AliMUONRawStreamTracker+;
#pragma link C++ class AliMUONRawStreamTrigger+;

// debug
#pragma link C++ class AliMUONCheck+;

#endif


