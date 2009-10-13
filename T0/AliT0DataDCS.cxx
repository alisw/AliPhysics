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
#include <iostream>

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
	fMPDcentA(0),
        fMPDcentC(0),
        fMPDsemiCentA(0),
        fMPDsemiCentC(0),
        fTVDCtop(0),
        fTVDCbottom(0),
	fMPDmode(0),
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
	fMPDcentA(0),
	fMPDcentC(0),
	fMPDsemiCentA(0),
	fMPDsemiCentC(0),
	fTVDCtop(0),
	fTVDCbottom(0),
	fMPDmode(0),
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
  fMPDcentA(0),
  fMPDcentC(0),
  fMPDsemiCentA(0),
  fMPDsemiCentC(0),
  fTVDCtop(0),
  fTVDCbottom(0),
  fMPDmode(0),
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
  fMPDcentA=data.fMPDcentA;
  fMPDcentC=data.fMPDcentC;
  fMPDsemiCentA=data.fMPDsemiCentA;
  fMPDsemiCentC=data.fMPDsemiCentC;
  fTVDCtop=data.fTVDCtop;
  fTVDCbottom=data.fTVDCbottom;
  fMPDmode=data.fMPDmode;
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
		Int_t t0_MPDcentA=0;
	        Int_t t0_MPDcentC=0;
	        Int_t t0_MPDsemiCentA=0;
	        Int_t t0_MPDsemiCentC=0;
        	Int_t t0_TVDCtop=0;
	        Int_t t0_TVDCbottom=0;  		
		Int_t t0_MPDmode=0;
	
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
		for(Int_t k=0; k<kLV; k++)
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
                      t0_scaler_sec[j-kScalers]+= aValue->GetUInt();
                    }
		    fScalerSecMean[j-kScalers] = t0_scaler_sec[j-kScalers] / aliasEntr[j];
		  }
		  else if (j < 2*kScalers+kHV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_hv_imon[j-2*kScalers]+= aValue->GetFloat();
                    }
                    fHViA[j-2*kScalers] = t0_a_hv_imon[j-2*kScalers] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_hv_vmon[j-(2*kScalers+kHV)]+= aValue->GetFloat();
                    }
                    fHVvA[j-(2*kScalers+kHV)] = t0_a_hv_vmon[j-(2*kScalers+kHV)] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV+kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_lv_imon[j-(2*kScalers+2*kHV)]+= aValue->GetFloat();
                    }
                    fLViA[j-(2*kScalers+2*kHV)] = t0_a_lv_imon[j-(2*kScalers+2*kHV)] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+2*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_lv_vmon[j-(2*kScalers+2*kHV+kLV)]+= aValue->GetFloat();
                    }
                    fLVvA[j-(2*kScalers+2*kHV+kLV)] = t0_a_lv_vmon[j-(2*kScalers+2*kHV+kLV)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+3*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_hv_imon[j-(2*kScalers+2*kHV+2*kLV)]+= aValue->GetFloat();
                    }
                    fHViC[j-(2*kScalers+2*kHV+2*kLV)] = t0_c_hv_imon[j-(2*kScalers+2*kHV+2*kLV)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+2*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_hv_vmon[j-(2*kScalers+3*kHV+2*kLV)]+= aValue->GetFloat();
                    }
                    fHVvC[j-(2*kScalers+3*kHV+2*kLV)] = t0_c_hv_vmon[j-(2*kScalers+3*kHV+2*kLV)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+3*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_lv_imon[j-(2*kScalers+4*kHV+2*kLV)]+= aValue->GetFloat();
                    }
                    fLViC[j-(2*kScalers+4*kHV+2*kLV)] = t0_c_lv_imon[j-(2*kScalers+4*kHV+2*kLV)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_lv_vmon[j-(2*kScalers+4*kHV+3*kLV)]+= aValue->GetFloat();
                    }
                    fLVvC[j-(2*kScalers+4*kHV+3*kLV)] = t0_c_lv_vmon[j-(2*kScalers+4*kHV+3*kLV)] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_cfd_thre[j-(2*kScalers+4*kHV+4*kLV)]+= aValue->GetFloat();
                    }
                    fCFDtA[j-(2*kScalers+4*kHV+4*kLV)] = t0_a_cfd_thre[j-(2*kScalers+4*kHV+4*kLV)] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+2*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_a_cfd_walk[j-(2*kScalers+4*kHV+4*kLV+kCFD)]+= aValue->GetFloat();
                    }
                    fCFDwA[j-(2*kScalers+4*kHV+4*kLV+kCFD)] = t0_a_cfd_walk[j-(2*kScalers+4*kHV+4*kLV+kCFD)] / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+3*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_cfd_thre[j-(2*kScalers+4*kHV+4*kLV+2*kCFD)]+= aValue->GetFloat();
                    }
                    fCFDtC[j-(2*kScalers+4*kHV+4*kLV+2*kCFD)] = t0_c_cfd_thre[j-(2*kScalers+4*kHV+4*kLV+2*kCFD)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_c_cfd_walk[j-(2*kScalers+4*kHV+4*kLV+3*kCFD)]+= aValue->GetFloat();
                    }
                    fCFDwC[j-(2*kScalers+4*kHV+4*kLV+3*kCFD)] = t0_c_cfd_walk[j-(2*kScalers+4*kHV+4*kLV+3*kCFD)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_ac_trm[j-(2*kScalers+4*kHV+4*kLV+4*kCFD)]+= aValue->GetFloat();
                    }
                    fTRM[j-(2*kScalers+4*kHV+4*kLV+4*kCFD)] = t0_ac_trm[j-(2*kScalers+4*kHV+4*kLV+4*kCFD)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_ac_drm[j-(2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)]+= aValue->GetFloat();
                    }
                    fDRM[j-(2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)] = t0_ac_drm[j-(2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)] / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
		    for(Int_t l=0; l<aliasEntr[j]; l++)
		    {		
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_atten += aValue->GetInt();
                    }
                    fAtten = t0_atten / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+2*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_MPDcentA += aValue->GetInt();
                    }
                    fMPDcentA = t0_MPDcentA / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+3*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_MPDcentC += aValue->GetInt();
                    }
                    fMPDcentC = t0_MPDcentC / aliasEntr[j];
                  }
                  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+4*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_MPDsemiCentA += aValue->GetInt();
                    }
                    fMPDsemiCentA = t0_MPDsemiCentA / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+5*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_MPDsemiCentC += aValue->GetInt();
                    }
                    fMPDsemiCentC = t0_MPDsemiCentC / aliasEntr[j];
                  }
		  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+6*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_TVDCtop += aValue->GetInt();
                    }
                    fTVDCtop = t0_TVDCtop / aliasEntr[j];
                  }
  		  else if (j < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+7*kAtten)
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_TVDCbottom += aValue->GetInt();
                    }
                    fTVDCbottom = t0_TVDCbottom / aliasEntr[j];
                  }
		  else
                  {
                    aliasEntr[j] = aliasArr->GetEntries();
                    for(Int_t l=0; l<aliasEntr[j]; l++)
                    {
                      AliDCSValue *aValue=dynamic_cast<AliDCSValue*> (aliasArr->At(l));
                      t0_MPDmode += aValue->GetInt();
                    }
                    fMPDmode = t0_MPDmode / aliasEntr[j];
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
                  sindex.Form("%02d",i-kScalers);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+kHV)
                {
                  fAliasNames[i] = "t00_a_hv_imon_";
                  sindex.Form("%02d",i-2*kScalers);
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV)
                {
                  fAliasNames[i] = "t00_a_hv_vmon_";
                  sindex.Form("%02d",i-(2*kScalers+kHV));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV+kLV)
                {
                  fAliasNames[i] = "t00_a_lv_imon_";
                  sindex.Form("%01d",i-(2*kScalers+2*kHV));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+2*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_a_lv_vmon_";
                  sindex.Form("%01d",i-(2*kScalers+2*kHV+kLV));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+3*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_c_hv_imon_";
                  sindex.Form("%02d",i-(2*kScalers+2*kHV+2*kLV));
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+2*kLV)
                {
                  fAliasNames[i] = "t00_c_hv_vmon_";
                  sindex.Form("%02d",i-(2*kScalers+3*kHV+2*kLV));
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+3*kLV)
                {
                  fAliasNames[i] = "t00_c_lv_imon_";
                  sindex.Form("%01d",i-(2*kScalers+4*kHV+2*kLV));
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+4*kLV)
                {
                  fAliasNames[i] = "t00_c_lv_vmon_";
                  sindex.Form("%01d",i-(2*kScalers+4*kHV+3*kLV));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+kCFD)
                {
                  fAliasNames[i] = "t00_a_cfd_thre_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+2*kCFD)
                {
                  fAliasNames[i] = "t00_a_cfd_walk_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV+kCFD));
                  fAliasNames[i] += sindex;
                }
		  else if (i < 2*kScalers+4*kHV+4*kLV+3*kCFD)
                {
                  fAliasNames[i] = "t00_c_cfd_thre_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV+2*kCFD));
                  fAliasNames[i] += sindex;
                }
                else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD)
                {
                  fAliasNames[i] = "t00_c_cfd_walk_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV+3*kCFD));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM)
                {
                  fAliasNames[i] = "t00_ac_trm_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV+4*kCFD));
                  fAliasNames[i] += sindex;
                }
		 else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM)
                {
                  fAliasNames[i] = "t00_ac_drm_";
                  sindex.Form("%02d",i-(2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM));
                  fAliasNames[i] += sindex;
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+kAtten)
		{
                  fAliasNames[i] = "t00_ac_atten";
                }
                else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+2*kAtten)
                {
                  fAliasNames[i] = "t00_a_mpd_cent";
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+3*kAtten)
                {
                  fAliasNames[i] = "t00_c_mpd_cent";
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+4*kAtten)
                {
                  fAliasNames[i] = "t00_a_mpd_scent";
                }
                else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+5*kAtten)
                {
                  fAliasNames[i] = "t00_c_mpd_scent";
                }
		else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+6*kAtten)
                {
                  fAliasNames[i] = "t00_ac_tvdc_top";
                }
                else if (i < 2*kScalers+4*kHV+4*kLV+4*kCFD+kTRM+kDRM+7*kAtten)
                {
                  fAliasNames[i] = "t00_ac_tvdc_bottom";
                }
		else
		{
                  fAliasNames[i] = "t00_ac_mpd_mode";
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



