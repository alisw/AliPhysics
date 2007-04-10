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
#pragma link C++ class AliMUONTriggerSegmentation+;

// geometry 
#pragma link C++ class AliMUONMathieson+; 

// info classes 
#pragma link C++ class AliMUONConstants+; 
#pragma link C++ class AliMUONDataInterface+; 
#pragma link C++ class AliMUONLoader+; 
#pragma link C++ class AliMUONChamber+; 
#pragma link C++ class AliMUONChamberTrigger+; 
#pragma link C++ class AliMUONTriggerCircuit+; 
#pragma link C++ class AliMUONTriggerCrateStore+; 
#pragma link C++ class AliMUONLogger+;

// containers
#pragma link C++ class AliMUONData+; 
#pragma link C++ class AliMUONDataIterator+; 
#pragma link C++ class AliMUONDataDigitIterator+; 
#pragma link C++ class AliMUONPoints+; 
#pragma link C++ class AliMUONHit+; 
#pragma link C++ class AliMUONRawCluster+;
#pragma link C++ class AliMUONDigit+; 
#pragma link C++ class AliMUONGlobalTrigger+; 
#pragma link C++ class AliMUONRegionalTrigger+; 
#pragma link C++ class AliMUONLocalTrigger+; 
#pragma link C++ class AliMUONTriggerChamberEff+;

// raw data
#pragma link C++ class AliMUONDigitMaker+;

// calibration access
#pragma link C++ class AliMUONCalibrationData+;
#pragma link C++ class AliMUONCDB+;

#endif


