#ifndef TEST_RICH_PRE_PROCESSOR_H
#define TEST_RICH_PRE_PRECESSOR_H

#include "AliPreprocessor.h"

//
// Prototype of RICH Preprocessor
//

class TestRICHPreprocessor: public AliPreprocessor {
public:
	TestRICHPreprocessor();
	TestRICHPreprocessor(AliShuttleInterface* shuttle);

protected:

        virtual void Initialize(Int_t run, UInt_t startTime,
                        UInt_t endTime);

        virtual UInt_t Process(TMap* valueSet);

	ClassDef(TestRICHPreprocessor, 0);
};

#endif
