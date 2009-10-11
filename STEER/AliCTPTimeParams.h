#ifndef ALICTPTIMEPARAMS_H
#define ALICTPTIMEPARAMS_H


#include<TNamed.h>
#include<TObjArray.h>

class TNamed;

class AliCTPTimeParams : public TNamed {

public:
			AliCTPTimeParams();
              virtual   ~AliCTPTimeParams();
	
                        AliCTPTimeParams(const AliCTPTimeParams &timeparams);
      AliCTPTimeParams& operator=(const AliCTPTimeParams& timeparams);		

             
	static AliCTPTimeParams* LoadCTPTimeParams(TString filename);	
	static AliCTPTimeParams* LoadCTPTimeParamsFromString(const char* timeparams);
		Bool_t ProcessCTPTimeParamsLine(const char* line);
		Bool_t WriteCTPTimeParamsOCDB();
		Bool_t GetCTPTimeParamsDAQLog();
		  void AddInput( TString& inputName, UInt_t& inputLevel, UInt_t inputDelay, TString inputEdge );
		  void AddDelayL0L1L2(UInt_t delayL1L0, UInt_t delayL2L0);
	  virtual void Print(const Option_t* opt="") const;
		//Setters

		//Getters
		UInt_t GetDelayL1L0()   const { return fDelayL1L0; }
	        UInt_t GetDelayL2L0()   const { return fDelayL2L0; }
/*		UInt_t* GetDelayInputs() { return fDelayInputs; }
		Bool_t* GetInputIndex()  { return fInputIndex; }		
		UInt_t* GetInputLevel()  { return fInputLevel; }*/

                enum {kNMaxInputs = 60}; //CTP can manage up to 60 trigger detector inputs
private:
			UInt_t fDelayL1L0;
			UInt_t fDelayL2L0;
			TObjArray fCTPInputTimeParams;
/*
			UInt_t fInputDelay[kNMaxInputs];
			Bool_t fInputFlag[kNMaxInputs];
			UInt_t fInputLevel[kNMaxInputs];
		       TString fInputName[kNMaxInputs];
                        Char_t fInputEdge[kNMaxInputs];
*/
          static const TString fgkCTPTimeParamsFileName;   //Name of file containing the CTPTimeParams

  ClassDef( AliCTPTimeParams, 1 ) 
};

#endif
