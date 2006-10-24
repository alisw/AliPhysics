#ifndef ALI_START_PREPROCESSOR_H
#define ALI_START_PREPRECESSOR_H

#include "AliPreprocessor.h"

//
//// Example of a Shuttle Preprocessor
////
//
class AliSTARTPreprocessor: public AliPreprocessor 
{
  public:
        AliSTARTPreprocessor(const char* detector, AliShuttleInterface* shuttle);
	virtual ~AliSTARTPreprocessor();
  
  protected:
//        virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
	virtual UInt_t Process(TMap* dcsAliasMap);

  private:
	//AliSTARTCalc *fData;
	
	ClassDef(AliSTARTPreprocessor, 1);

};

#endif
