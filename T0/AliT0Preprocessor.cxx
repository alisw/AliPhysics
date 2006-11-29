#include "AliT0Preprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliT0Calc.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TNamed.h>
#include "AliT0Dqclass.h"

ClassImp(AliT0Preprocessor)

//____________________________________________________
AliT0Preprocessor::AliT0Preprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle)
{

}

AliT0Preprocessor::~AliT0Preprocessor()
{

}

UInt_t AliT0Preprocessor::Process(TMap* dcsAliasMap )
{

	if(!dcsAliasMap) return 0;

        TObjArray *aliasArr;
       // AliDCSValue *aValue;
        Float_t hv[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for(int j=0; j<24; j++){
		TString aliasName =Form("T0HV%d",j);
                // printf("aliasname: %s\n",aliasName.Data());
                aliasArr = dynamic_cast<TObjArray*> (dcsAliasMap->GetValue(aliasName.Data()));
                if(!aliasArr){
                        AliError(Form("Alias %s not found!", aliasName.Data()));
                        continue;
                }
                AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(0));
                // printf("I'm here! %f %x\n", aValue->GetFloat(), aValue->GetTimeStamp());
               hv[j]= aValue->GetFloat()*100;
	       //Float_t timestamp= (Float_t) (aValue->GetTimeStamp());
		// printf("hello! hv = %f timestamp = %f\n" ,hv[j], timestamp);

	}
	Float_t numbers[24];

	AliT0Calc *calibdata = new AliT0Calc();	
	
	const char* TimefileName = GetFile(kDAQ, "TIME", "LDC0");
	
	
	if(TimefileName){
		TFile *file = TFile::Open(TimefileName);
		if(!file || !file->IsOpen()) 
		{
		  printf("File from DAQ does not exist.");
		  return 0;
		} 
		AliT0Dqclass *tempdata = dynamic_cast<AliT0Dqclass*> (file->Get("Time"));
		for(Int_t i=0;i<24;i++){
			numbers[i] = tempdata->GetTime(i);
		//	printf("\nnumbers: %f\n",numbers[i]);
		}
        	file->Close();
		delete tempdata;
	}
	else {return 0;}
	calibdata->SetTime(numbers, hv);
	calibdata->Print();        

	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Tomek&Michal");
	metaData.SetComment("This preprocessor returns time to be used for reconstruction.");
	
	UInt_t result = Store("Calib","Data", calibdata, &metaData);
	delete calibdata;
	return result;

}

	
