#ifndef TEST_ITS_PRE_PROCESSOR_H
#define TEST_ITS_PRE_PRECESSOR_H

#include "AliCDBPreProcessor.h"

class TestITSPreProcessor: public AliCDBPreProcessor {
public:
	TestITSPreProcessor();

protected:

        virtual void Initialize(Int_t run, UInt_t startTime, 
                        UInt_t endTime);

        virtual void Finalize();

        virtual void Process(const char* alias, TObjArray& valueSet, 
                        Bool_t hasError);

	ClassDef(TestITSPreProcessor, 0);
};

#endif
