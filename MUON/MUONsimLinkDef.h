#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// builder
#pragma link C++ class AliMUONResponseFactory+;  

// digitizer
#pragma link C++ class AliMUONHitMapA1+; 
#pragma link C++ class AliMUONDigitizer+; 
#pragma link C++ class AliMUONDigitizerv2+; 
#pragma link C++ class AliMUONDigitizerV3+; 
#pragma link C++ class AliMUONSDigitizerv1+;  
#pragma link C++ class AliMUONSDigitizerV2+;  
#pragma link C++ class AliMUONTriggerDecision+; 
#pragma link C++ class AliMUONTriggerDecisionV1+; 


#pragma link C++ class AliMUONTest+; 

// response
#pragma link C++ class AliMUONSt1Response+;
#pragma link C++ class AliMUONSt1ElectronicElement+; 
#pragma link C++ class AliMUONSt1SpecialMotif+; 
#pragma link C++ class AliMUONSt1ResponseParameter+; 
#pragma link C++ class AliMUONSt1ResponseRule+; 
#pragma link C++ class AliMUONSt1IniReader+; 
#pragma link C++ namespace decoder; 
#pragma link C++ class AliMUONResponse+; 
#pragma link C++ class AliMUONResponseV0+;
#pragma link C++ class AliMUONResponseTrigger+; 
#pragma link C++ class AliMUONResponseTriggerV1+;

// trigger

#pragma link C++ class AliMUONTrigger+;
#pragma link C++ class AliMUONTriggerEfficiencyCells+;
#pragma link C++ class AliMUONLocalTriggerBoard+;
#pragma link C++ class AliMUONRegionalTriggerBoard+;
#pragma link C++ class AliMUONTriggerBoard+;
#pragma link C++ class AliMUONTriggerCrate+;
#pragma link C++ class AliMUONTriggerElectronics+;
#pragma link C++ class AliMUONGlobalTriggerBoard+;

// raw

#pragma link C++ class AliMUONRawWriter+;

#endif


