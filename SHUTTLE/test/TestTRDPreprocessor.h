#ifndef TEST_TRD_PRE_PROCESSOR_H
#define TEST_TRD_PRE_PRECESSOR_H

#include "AliPreprocessor.h"

//
// Prototype of TRD Preprocessor
//

class TestTRDPreprocessor: public AliPreprocessor {
public:
	TestTRDPreprocessor();
	TestTRDPreprocessor(AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime,
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestTRDPreprocessor, 0);
};

#endif
