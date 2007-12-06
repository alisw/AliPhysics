#include "AliT0DataDCS.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TCanvas.h>

ClassImp(AliT0DataDCS)

//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fIsProcessed(kFALSE)
{
}

//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fIsProcessed(kFALSE)
{
	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

	Init();

}

//---------------------------------------------------------------
AliT0DataDCS::~AliT0DataDCS() {

}

//---------------------------------------------------------------
Bool_t AliT0DataDCS::ProcessData(TMap& aliasMap)
{
 	        Float_t t0_scaler[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		TObjArray *aliasArr;
		AliDCSValue* aValue;
        	for(int j=0; j<kNAliases; j++)
        	{
                  aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
                  if(!aliasArr)
                  {
                        AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
                        continue;
                  }
                  // Introduce(j, aliasArr);

                  if(aliasArr->GetEntries()<2)
                  {
                        AliError(Form("Alias %s has just %d entries!",
                                        fAliasNames[j].Data(),aliasArr->GetEntries()));
                        continue;
                  }

		  for(int j=0; j<32; j++)
		  {
		    TString aliasName =Form("t00_ac_scaler_%d",j);
                    // printf("aliasname: %s\n",aliasName.Data());
                    //aliasArr = dynamic_cast<TObjArray*> (aliasMap->GetValue(aliasName.Data()));
                    if(!aliasArr)
		    {
                      AliError(Form("Alias %s not found!", aliasName.Data()));
                      continue;
                    }
                  AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(j));
                  // printf("I'm here! %f %x\n", aValue->GetFloat(), aValue->GetTimeStamp());
		  t0_scaler[j]= aValue->GetFloat();
	          }
	        }
	CalcScalerMean(t0_scaler);
	fIsProcessed=kTRUE;
}

//---------------------------------------------------------------
void AliT0DataDCS::Init()
{

	TH1::AddDirectory(kFALSE);

	for(int i=0;i<kNAliases;i++)
	{
		fAliasNames[i] = "DCSAlias";
		fAliasNames[i] += i;
	}

}

//---------------------------------------------------------------
void AliT0DataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr){

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}

//---------------------------------------------------------------

void AliT0DataDCS::CalcScalerMean(Float_t *t0_scal)
{
  for (Int_t i=0; i<23; i++)
    {
        SetScalerMean(i,t0_scal[i]);
    }
}


