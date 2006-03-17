#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUON+; 
#pragma link C++ class AliMUONv1+; 

// mapping & segmentation (change in progress)
#pragma link C++ class AliMUONSt12QuadrantSegmentation+; 
#pragma link C++ class AliMUONSt345SlatSegmentation+;
#pragma link C++ class AliMUONSt345SlatSegmentationV2+;
#pragma link C++ class AliMUONTriggerSegmentation+;
#pragma link C++ class AliMUONTriggerSegmentationV2+;

// geometry 
#pragma link C++ class AliMUONMathieson+; 
#pragma link C++ class AliMUONCommonGeometryBuilder+;
#pragma link C++ class AliMUONSt1GeometryBuilder+; 
#pragma link C++ class AliMUONSt1GeometryBuilderV2+; 
#pragma link C++ class AliMUONSt2GeometryBuilder+; 
#pragma link C++ class AliMUONSt2GeometryBuilderV2+; 
#pragma link C++ class AliMUONSlatGeometryBuilder+; 
#pragma link C++ class AliMUONTriggerGeometryBuilder+; 

// info classes 
#pragma link C++ class AliMUONConstants+; 
#pragma link C++ class AliMUONTriggerConstants+; 
#pragma link C++ class AliMUONDataInterface+; 
#pragma link C++ class AliMUONLoader+; 
#pragma link C++ class AliMUONChamber+; 
#pragma link C++ class AliMUONChamberTrigger+; 
#pragma link C++ class AliMUONTriggerCircuit+; 

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
#pragma link C++ class AliMUONLocalTrigger+; 
#pragma link C++ class AliMUONTriggerLut+; 

// calibration
#pragma link C++ class AliMUONV2DStore+;
#pragma link C++ class AliMUONV1DStore+;
#pragma link C++ class AliMUON2DMap+;
#pragma link C++ class AliMUON1DArray+;
#pragma link C++ class AliMUONVCalibParam+;
#pragma link C++ class AliMUONCalibParam1I+;
#pragma link C++ class AliMUONCalibParam2F+;
#pragma link C++ class AliMUONCalibrationData+;

// display
#pragma link C++ class AliMUONDisplay+; 
#pragma link C++ class AliMUONRecoCheck+; 

// segmentation
#pragma link C++ class AliMUONSegFactory+;

// raw data
#pragma link C++ class AliMUONRawReader+;
#pragma link C++ class AliMUONRawStream+;
#pragma link C++ class AliMUONSubEventTracker+;
#pragma link C++ class AliMUONDDLTrigger+;
#pragma link C++ class AliMUONDDLTracker+;
#pragma link C++ class AliMUONSubEventTrigger+;
#pragma link C++ class AliMUONScalerEventTrigger+;

// debug
#pragma link C++ class AliMUONCheck+;

#endif


