/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *  
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log: AliT0DataDCS.cxx,v $
 Revision   2008/01/30 
Fetching data points from DCS, calculating mean and storing data to Reference DB 
 
 Version 1.1  2006/10
Preliminary test version (T.Malkiewicz)
*/

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

// AliT0DataDCS class
// declaring DCS aliases for T0
// fetching T0 data points from DCS, 
// calculating mean values for the entire run
// and storing the result to Reference DB

ClassImp(AliT0DataDCS)

//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS():
	TObject(),
	fRun(0),	
	fStartTime(0),
	fEndTime(0),
	fStartTimeDCSQuery(0),
	fEndTimeDCSQuery(0),
	fAtten(0.),
	fIsProcessed(kFALSE)
{
  for(Int_t i=0;i<kScalers;i++) 
  {
    fScalerMean[i]=0.;
    fScalerSecMean[i]=0.;
  }
  for(Int_t i=0;i<kHV;i++)
  {
    fHViA[i]=0.;
    fHVvA[i]=0.;
    fHViC[i]=0.;
    fHVvC[i]=0.;
  }
  for(Int_t i=0;i<kLV;i++)
  {
    fLViA[i]=0.;
    fLVvA[i]=0.;
    fLViC[i]=0.;
    fLVvC[i]=0.;
  }
  for(Int_t i=0;i<kTRM;i++)
  {
    fTRM[i]=0.;
  }
  for(Int_t i=0;i<kDRM;i++)
  {
    fDRM[i]=0.;
  }

}
//---------------------------------------------------------------
AliT0DataDCS::AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fStartTimeDCSQuery(startTimeDCSQuery),
	fEndTimeDCSQuery(endTimeDCSQuery),
	fAtten(0.),
	fIsProcessed(kFALSE)
{
  for(Int_t i=0;i<kScalers;i++)
  {
    fScalerMean[i]=0.;
    fScalerSecMean[i]=0.;
  }
  for(Int_t i=0;i<kHV;i++)
  {
    fHViA[i]=0.;
    fHVvA[i]=0.;
    fHViC[i]=0.;
    fHVvC[i]=0.;
  }
  for(Int_t i=0;i<kLV;i++)
  {
    fLViA[i]=0.;
    fLVvA[i]=0.;
    fLViC[i]=0.;
    fLVvC[i]=0.;
  }
  for(Int_t i=0;i<kTRM;i++)
  {
    fTRM[i]=0.;
  }
  for(Int_t i=0;i<kDRM;i++)
  {
    fDRM[i]=0.;
  }

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s \n\tStartTime DCS Query %s \n\tEndTime DCS Query %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString(),
	TTimeStamp(startTimeDCSQuery).AsString(), 
        TTimeStamp(endTimeDCSQuery).AsString()));

	Init();

}
//---------------------------------------------------------------

AliT0DataDCS::AliT0DataDCS(const AliT0DataDCS & data):
  TObject(), 
  fRun(0),
  fStartTime(0),
  fEndTime(0),
  fStartTimeDCSQuery(0),
  fEndTimeDCSQuery(0),
  fAtten(0.),
  fIsProcessed(kFALSE)
{

// copy constructor

  fRun=data.fRun;
  fStartTime=data.fStartTime;
  fEndTime=data.fEndTime;
  fStartTimeDCSQuery=data.fStartTimeDCSQuery;
  fEndTimeDCSQuery=data.fEndTimeDCSQuery;
  fIsProcessed=data.fIsProcessed;
  fAtten=data.fAtten;

  for(int i=0;i<kNAliases;i++) 
  {
    fAliasNames[i]=data.fAliasNames[i];
  }
  
  for(Int_t i=0;i<kScalers;i++)
  {
    fScalerMean[i]=data.fScalerMean[i];
    fScalerSecMean[i]=data.fScalerSecMean[i];
  }
  for(Int_t i=0;i<kHV;i++)
  {
    fHViA[i]=data.fHViA[i];
    fHVvA[i]=data.fHVvA[i];
    fHViC[i]=data.fHViC[i];
    fHVvC[i]=data.fHVvC[i];
  }
  for(Int_t i=0;i<kLV;i++)
  {
    fLViA[i]=data.fLViA[i];
    fLVvA[i]=data.fLVvA[i];
    fLViC[i]=data.fLViC[i];
    fLVvC[i]=data.fLVvC[i];
  }
  for(Int_t i=0;i<kTRM;i++)
  {
    fTRM[i]=data.fTRM[i];
  }
  for(Int_t i=0;i<kDRM;i++)
  {
    fDRM[i]=data.fDRM[i];
  }
}
//---------------------------------------------------------------

AliT0DataDCS& AliT0DataDCS:: operator=(const AliT0DataDCS & data) { 

// assignment operator

  this->fRun=data.fRun;
  this->fStartTime=data.fStartTime;
  this->fEndTime=data.fEndTime;
  this->fStartTimeDCSQuery=data.fStartTimeDCSQuery;
  this->fEndTimeDCSQuery=data.fEndTimeDCSQuery;

  for(int i=0;i<kNAliases;i++)
  {
    this->fAliasNames[i]=data.fAliasNames[i];
  }

  return *this;
}

//---------------------------------------------------------------
AliT0DataDCS::~AliT0DataDCS() 
{
  delete fScalerMean;
  delete fScalerSecMean;
}

//---------------------------------------------------------------
Bool_t AliT0DataDCS::ProcessData(TMap& aliasMap)
{
 	        Int_t t0_scaler[kScalers];
		Int_t t0_scaler_sec[kScalers];
		Int_t aliasEntr[kNAliases];
		Float_t t0_a_hv_imon[kHV];
		Float_t t0_a_hv_vmon[kHV];
		Float_t t0_a_lv_imon[kLV];
                Float_t t0_a_lv_vmon[kLV];
                Float_t t0_c_hv_imon[kHV];
                Float_t t0_c_hv_vmon[kHV];
                Float_t t0_c_lv_imon[kLV];
                Float_t t0_c_lv_vmon[kLV];
		Float_t t0_a_cfd_thre[kCFD];
                Float_t t0_a_cfd_walk[kCFD];
                Float_t t0_c_cfd_thre[kCFD];
                Float_t t0_c_cfd_walk[kCFD];
                Float_t t0_ac_trm[kTRM];
                Float_t t0_ac_drm[kDRM];
		Float_t t0_atten=0.;

		TObjArray *aliasArr;
        	for(Int_t k=0; k<kScalers; k++)
		{
		   t0_scaler[k]=0;
                   t0_scaler_sec[k]=0;
		}

		for(Int_t k=0; k<kHV; k++)
                {
                   t0_a_hv_imon[k]=0.;
                   t0_a_hv_vmon[k]=0.;
		   t0_c_hv_imon[k]=0.;
                   t0_c_hv_vmon[k]=0.;
                }
		for(Int_t k=0; k<kHV; k++)
                {
                   t0_a_lv_imon[k]=0.;
                   t0_a_lv_vmon[k]=0.;
                   t0_c_lv_imon[k]=0.;
                   t0_c_lv_vmon[k]=0.;
                }
		for(Int_t k=0; k<kCFD; k++)
                {
                   t0_a_cfd_thre[k]=0.;
                   t0_a_cfd_walk[k]=0.;
 		   t0_c_cfd_thre[k]=0.;
                   t0_c_cfd_walk[k]=0.;
                }
		for(Int_t k=0; k<kTRM; k++)
                {
                   t0_ac_trm[k]=0.;
		}
		for(Int_t k=0; k<kDRM; k++)
                {
                   t0_ac_drm[k]=0.;
                }

// here starts the main loop
		for(Int_t j=0; j<kNAliases; j++)
        	{
		  aliasEntr[j]=0;
		  for (Int_t k=0;k<32;k++) 
		  {
		    t0_scaler[k]=0;
		    t0_scaler_sec[k]=0;	
		    
      		  }
		  aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
                  if(!aliasArr)
                  {
                        AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
                        continue;
                  }
		  Introduce(j, aliasArr);
                  if(aliasArr->GetEntries()<2)
                  {
                        AliError(Form("Alias %s has just %d entries!",
                                        fAliasNames[j].Data(),aliasArr->GetEntries()));
                        continue;
                  }
		  if (j < kScalers)
		  { 
		    aliasEntr[j] = aliasArr->GetEntries();
		    for(Int_t l=0; l<aliasEntr[j]; l++)
		    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
		      t0_scaler[j]+= aValue->GetUInt();
	            }
		    fScalerMean[j] = t0_scaler[j] / aliasEntr[j];
		  }
		  else if (j < 2*kScalers)
		  {
		    aliasEntr[j] = aliasArr->GetEntries();		
   	            for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_scaler_sec[j]+= aValue->GetUInt();
                    }
		    fScalerSecMean[j] = t0_scaler_sec[j] / aliasEntr[j];
		  }
		  else if (j < 2*kScalers+kHV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_hv_imon[j]+= aValue->GetFloat();
                    }
                    fHViA[j] = t0_a_hv_imon[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_hv_vmon[j]+= aValue->GetFloat();
                    }
                    fHVvA[j] = t0_a_hv_vmon[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV+kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_lv_imon[j]+= aValue->GetFloat();
                    }
                    fLViA[j] = t0_a_lv_imon[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_lv_vmon[j]+= aValue->GetFloat();
                    }
                    fLVvA[j] = t0_a_lv_vmon[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+3*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_hv_imon[j]+= aValue->GetFloat();
                    }
                    fHViC[j] = t0_c_hv_imon[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_hv_vmon[j]+= aValue->GetFloat();
                    }
                    fHVvC[j] = t0_c_hv_vmon[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+3*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_lv_imon[j]+= aValue->GetFloat();
                    }
                    fLViC[j] = t0_c_lv_imon[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_lv_vmon[j]+= aValue->GetFloat();
                    }
                    fLVvC[j] = t0_c_lv_vmon[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_cfd_thre[j]+= aValue->GetFloat();
                    }
                    fCFDtA[j] = t0_a_cfd_thre[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+2*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_cfd_walk[j]+= aValue->GetFloat();
                    }
                    fCFDwA[j] = t0_a_cfd_walk[j] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+3*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_cfd_thre[j]+= aValue->GetFloat();
                    }
                    fCFDtC[j] = t0_c_cfd_thre[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_cfd_walk[j]+= aValue->GetFloat();
                    }
                    fCFDwC[j] = t0_c_cfd_walk[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_ac_trm[j]+= aValue->GetFloat();
                    }
                    fTRM[j] = t0_ac_trm[j] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_ac_drm[j]+= aValue->GetFloat();
                    }
                    fDRM[j] = t0_ac_drm[j] / aliasEntr[j];
                  }
                  else
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
		    for(Int_t l=0; l<aliasEntr[j]; l++)
		    {		
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_atten += aValue->GetInt();
                    }
                    fTRM[j] = t0_atten / aliasEntr[j];
                  }
			

	        }
	fIsProcessed=kTRUE;
	return kTRUE;
}

//---------------------------------------------------------------
void AliT0DataDCS::Init()
{
	TString sindex;
	for(int i=0;i<kNAliases;i++)
	{	
		if (i<kScalers)
		{
		  fAliasNames[i] = "t00_ac_scaler_";
		  sindex.Form("%02d",i);
	          fAliasNames[i] += sindex;
		}
		else if (i < 2*kScalers)
		{
		  fAliasNames[i] = "t00_ac_scaler_sec_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+kHV)
                {
                  fAliasNames[i] = "t00_a_hv_imon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV)
                {
                  fAliasNames[i] = "t00_a_hv_vmon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV+kLV)
                {
                  fAliasNames[i] = "t00_a_lv_imon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_a_lv_vmon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+3*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_c_hv_imon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_c_hv_vmon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+3*kLV)
                {
                  fAliasNames[i] = "t00_c_lv_imon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+4*kLV)
                {
                  fAliasNames[i] = "t00_c_lv_vmon_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+kCFD)
                {
                  fAliasNames[i] = "t00_a_cfd_thre_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+2*kCFD)
                {
                  fAliasNames[i] = "t00_a_cfd_walk_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		  else if (i < 2*kScalers+4*kHV+4*kLV+3*kCFD)
                {
                  fAliasNames[i] = "t00_c_cfd_thre_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD)
                {
                  fAliasNames[i] = "t00_c_cfd_walk_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)
                {
                  fAliasNames[i] = "t00_ac_trm_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		 else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM)
                {
                  fAliasNames[i] = "t00_ac_drm_";
                  sindex.Form("%02d",i);
                  fAliasNames[i] += sindex;
                }
		else
		{
                  fAliasNames[i] = "t00_atten";
                }

	}

}

//---------------------------------------------------------------
void AliT0DataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)const
{

	int entries=aliasArr->GetEntries();
	AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
	AliInfo(Form("    	%d DP values collected",entries));

}



