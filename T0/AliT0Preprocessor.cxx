#include "AliT0Preprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliT0Calc.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TObjString.h>
#include <TNamed.h>
#include "AliT0Dqclass.h"

ClassImp(AliT0Preprocessor)

//____________________________________________________
AliT0Preprocessor::AliT0Preprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("T00", shuttle)
{

}

AliT0Preprocessor::~AliT0Preprocessor()
{

}

UInt_t AliT0Preprocessor::Process(TMap* dcsAliasMap )
{

	if(!dcsAliasMap) return 1;

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
	
        TList* list = GetFileSources(kDAQ, "TIME");
        if (list)
        {
		TIter iter(list);
      		TObjString *source;
      		while ((source = dynamic_cast<TObjString *> (iter.Next()))) 
		{
			const char* TimefileName = GetFile(kDAQ, "TIME", source->GetName());
	
			if (TimefileName)
			{
				Log(Form("File with Id TIME found in source %s!", source->GetName()));
				TFile *file = TFile::Open(TimefileName);
				if(!file || !file->IsOpen()) 
				{
		  			Log(Form("Error opening file with Id TIME from source %s!", source->GetName()));
		 			return 1;
				} 
				AliT0Dqclass *tempdata = dynamic_cast<AliT0Dqclass*> (file->Get("Time"));
				if (!tempdata) 
				{
					Log("Could not find key \"Time\" in DAQ file!");
					return 1;
				}
				for(Int_t i=0;i<24;i++){
					numbers[i] = tempdata->GetTime(i);
					//	printf("\nnumbers: %f\n",numbers[i]);
				}
        			file->Close();
				delete tempdata;
			} else {
		  		Log(Form("Could not find file with Id TIME in source %s!", source->GetName()));
				return 1;
			}
			calibdata->SetTime(numbers, hv);
			calibdata->Print();
		} 
	} else {
		Log("No sources for Id TIME found!");
	}        

	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Tomek&Michal");
	metaData.SetComment("This preprocessor returns time to be used for reconstruction.");

	Bool_t result = Store("Calib","Data", calibdata, &metaData);
	delete calibdata;
	if(result == kTRUE) 
	  {return 0;}
	else {return 1;}
}

	
