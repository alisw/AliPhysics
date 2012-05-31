#ifndef TEST_TPC_PRE_PROCESSOR_H
#define TEST_TPC_PRE_PRECESSOR_H

#include "AliPreprocessor.h"
#include "AliTPCDataDCS.h"

//
// Example of a Shuttle Preprocessor
//

class TestTPCPreprocessor: public AliPreprocessor {
public:

	TestTPCPreprocessor();
	TestTPCPreprocessor(AliShuttleInterface* shuttle);
	~TestTPCPreprocessor();

protected:

        virtual void Initialize(Int_t run, UInt_t startTime,
                        UInt_t endTime);

        virtual UInt_t Process(TMap* aliasMap);

private:
	Int_t fRun; 		// Run number
	Int_t fStartTime;	// Start time
	Int_t fEndTime;		// End time

	AliTPCDataDCS *fData; 	// pointer to the data container and analyzer class

	ClassDef(TestTPCPreprocessor, 0);
};

#endif
