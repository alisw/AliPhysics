#ifndef TEST_ITS_PRE_PROCESSOR_H
#define TEST_ITS_PRE_PRECESSOR_H

#include "AliPreprocessor.h"

//
// Example of a Shuttle Preprocessor
//

class TestITSPreprocessor: public AliPreprocessor {
public:
	TestITSPreprocessor();
	TestITSPreprocessor(const char* detector, AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime, 
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestITSPreprocessor, 0);
};

#endif
