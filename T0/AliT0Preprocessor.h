#ifndef ALI_T0_PREPROCESSOR_H
#define ALI_T0_PREPRECESSOR_H

#include "AliPreprocessor.h"

//
//// Example of a Shuttle Preprocessor
////
//
class AliT0Preprocessor: public AliPreprocessor 
{
  public:
        AliT0Preprocessor(const char* detector, AliShuttleInterface* shuttle);
	virtual ~AliT0Preprocessor();
  
  protected:
//        virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
	virtual UInt_t Process(TMap* dcsAliasMap);

  private:
	//AliT0Calc *fData;
	
	ClassDef(AliT0Preprocessor, 1)
};

typedef AliT0Preprocessor AliSTARTPreprocessor; // for backward compatibility

#endif
