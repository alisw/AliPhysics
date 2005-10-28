#include "TestITSPreProcessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TList.h>
#include <TTimeStamp.h>

ClassImp(TestITSPreProcessor)

TestITSPreProcessor::TestITSPreProcessor():
	AliCDBPreProcessor("ITS")
{

}

void TestITSPreProcessor::Initialize(Int_t run, UInt_t startTime, 
	UInt_t endTime) 
{

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, 
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));
}

void TestITSPreProcessor::Finalize() {
	AliInfo("Finalizing...");
}

void TestITSPreProcessor::Process(const char* alias, TList& valueSet, 
	Bool_t hasError)
{
	AliInfo(Form("Alias %s, hasError: %d", alias, hasError));	

	TString output;
	
	TIter iter(&valueSet);
	AliDCSValue* aValue;
	while ((aValue = (AliDCSValue*) iter.Next())) {
		output += aValue->ToString();
		output += '\n';
	}	
	output += '\n';

	AliInfo(output);

	AliCDBMetaData metaData;
	metaData.SetComment("This is a test!");
	
	Store(alias, &valueSet, &metaData);	
}

