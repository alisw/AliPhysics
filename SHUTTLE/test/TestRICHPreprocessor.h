#ifndef TEST_HMPID_PRE_PROCESSOR_H
#define TEST_HMPID_PRE_PRECESSOR_H

#include "AliPreprocessor.h"

//
// Prototype of HMPID Preprocessor
//

class TestHMPIDPreprocessor: public AliPreprocessor {
public:
	TestHMPIDPreprocessor();
	TestHMPIDPreprocessor(AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime,
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestHMPIDPreprocessor, 0);
};

#endif
