#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUON+; 
#pragma link C++ class AliMUONv1+; 

// mapping & segmentation (change in progress)
#pragma link C++ class AliMUONSt12QuadrantSegmentation-; 
#pragma link C++ class AliMUONSt345SlatSegmentation+;
#pragma link C++ class AliMUONSt345SlatSegmentationV2-;
#pragma link C++ class AliMUONTriggerSegmentation+;

// geometry 
#pragma link C++ class AliMUONMathieson+; 
#pragma link C++ class AliMUONGeometryDEIndexing+;
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
#pragma link C++ class AliMUONPoints+; 
#pragma link C++ class AliMUONHit+; 
#pragma link C++ class AliMUONRawCluster+;
#pragma link C++ class AliMUONDigit+; 
#pragma link C++ class AliMUONTransientDigit+;
#pragma link C++ class AliMUONGlobalTrigger+; 
#pragma link C++ class AliMUONLocalTrigger+; 
#pragma link C++ class AliMUONTriggerLut+; 

// display
#pragma link C++ class AliMUONDisplay+; 
#pragma link C++ class AliMUONRecoCheck+; 

#pragma link C++ class AliMUONSegmentationManager+;
#endif


