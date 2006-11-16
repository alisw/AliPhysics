#ifndef TEST_ITS_PRE_PROCESSOR_SPD_H
#define TEST_ITS_PRE_PRECESSOR_SPD_H

#include "AliPreprocessor.h"

//
// Example of a Shuttle Preprocessor
//

class TestITSPreprocessorSPD: public AliPreprocessor {
public:
	TestITSPreprocessorSPD();
	TestITSPreprocessorSPD(AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime, 
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestITSPreprocessorSPD, 0);
};

#endif
