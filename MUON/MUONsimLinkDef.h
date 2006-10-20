#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// geometry 
#pragma link C++ class AliMUONCommonGeometryBuilder+;
#pragma link C++ class AliMUONSt1SpecialMotif+; 
#pragma link C++ class AliMUONSt1GeometryBuilder+; 
#pragma link C++ class AliMUONSt1GeometryBuilderV2+; 
#pragma link C++ class AliMUONSt2GeometryBuilder+; 
#pragma link C++ class AliMUONSt2GeometryBuilderV2+; 
#pragma link C++ class AliMUONSlatGeometryBuilder+; 
#pragma link C++ class AliMUONTriggerGeometryBuilder+; 

// builder
#pragma link C++ class AliMUONResponseFactory+;  

// digitizer
#pragma link C++ class AliMUONHitMapA1+; 
#pragma link C++ class AliMUONDigitizerV3+; 
#pragma link C++ class AliMUONSDigitizerV2+;  
#pragma link C++ class AliMUONTest+; 

// response
#pragma link C++ class AliMUONResponse+; 
#pragma link C++ class AliMUONResponseV0+;
#pragma link C++ class AliMUONResponseTrigger+; 
#pragma link C++ class AliMUONResponseTriggerV1+;

// trigger

#pragma link C++ class AliMUONTrigger+;
#pragma link C++ class AliMUONLocalTriggerBoard+;
#pragma link C++ class AliMUONRegionalTriggerBoard+;
#pragma link C++ class AliMUONTriggerBoard+;
#pragma link C++ class AliMUONTriggerCrate+;
#pragma link C++ class AliMUONTriggerElectronics+;
#pragma link C++ class AliMUONGlobalTriggerBoard+;

// raw

#pragma link C++ class AliMUONRawWriter+;

#endif


