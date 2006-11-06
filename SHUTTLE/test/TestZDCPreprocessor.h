#ifndef TEST_ZDC_PRE_PROCESSOR_SPD_H
#define TEST_ZDC_PRE_PRECESSOR_SPD_H

#include "AliPreprocessor.h"

//
// Example of a Shuttle Preprocessor
//

class TestZDCPreprocessor: public AliPreprocessor {
public:
	TestZDCPreprocessor();
	TestZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime, 
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestZDCPreprocessor, 0);
};

#endif
