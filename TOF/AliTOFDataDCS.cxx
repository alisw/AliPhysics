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
$Log$
Revision 1.2  2006/12/18 18:17:38  arcelli
Updated Aliases for DCS TOF datapoints (C.Zampolli)

Revision 1.1  2006/10/26 09:10:52  arcelli
Class for handling the TOF DCS data in the Shuttle (C.Zampolli)

*/  

#include "AliTOFDataDCS.h"

#include "AliDCSValue.h"
#include "AliLog.h"

#include "TString.h"
#include "AliTOFFormatDCS.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTimeStamp.h"
#include "TMap.h"

class TH2;
class AliCDBMetaData;
class TDatime;

// AliTOFDataDCS class
// main aim to introduce the aliases for the TOF DCS
// data points to be then
// stored in the OCDB, and to process them. 
// Process() method called by TOFPrepr

ClassImp(AliTOFDataDCS)

//---------------------------------------------------------------
AliTOFDataDCS::AliTOFDataDCS():
	TObject(),
	fRun(0),
	fStartTime(0),
	fEndTime(0),
	fIsProcessed(kFALSE)
{

  // main constructor 

  for(int i=0;i<kNHV;i++) {
    fHVvp[i]=0x0;
    fHVvn[i]=0x0;
    fHVip[i]=0x0;
    fHVin[i]=0x0;
  }
  
  for(int i=0;i<kNLV;i++) {
    fLVv[i]=0x0;
    fLVi[i]=0x0;
  }
  
  for(int i=0;i<kNLV33;i++) {
    fLVv33[i]=0x0;
    fLVi33[i]=0x0;
  }
  
  for(int i=0;i<kNLV50;i++) {
    fLVv50[i]=0x0;
    fLVi50[i]=0x0;
  }
  
  for(int i=0;i<kNLV48;i++) {
    fLVv48[i]=0x0;
    fLVi48[i]=0x0;
  }
  
  for(int i=0;i<kNFEEthr;i++) {
    fFEEthr[i]=0x0;
  }
  
  for(int i=0;i<kNFEEtfeac;i++) {
    fFEEtfeac[i]=0x0;
  }
  
  for(int i=0;i<kNFEEttrm;i++) {
    fFEEttrm[i]=0x0;
  }
  
  for(int i=0;i<3;i++) {
    fT[i]=0;
    fP[i]=0;
  }
  
}

//---------------------------------------------------------------
AliTOFDataDCS::AliTOFDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fIsProcessed(kFALSE)
{

  // constructor with arguments

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString()));

	Init();

}

//---------------------------------------------------------------

AliTOFDataDCS::AliTOFDataDCS(const AliTOFDataDCS & data):
  TObject(), 
  fRun(0),
  fStartTime(0),
  fEndTime(0),
  fIsProcessed(kFALSE)

{

// copy constructor

  fRun=data.fRun;
  fStartTime=data.fStartTime;
  fEndTime=data.fEndTime;
  fIsProcessed=data.fIsProcessed;

  for(int i=0;i<kNAliases;i++) {
    fAliasNames[i]=data.fAliasNames[i];
  }
 
  for(int i=0;i<kNHV;i++) {
    fHVvp[i]=data.fHVvp[i];
    fHVvn[i]=data.fHVvn[i];
    fHVip[i]=data.fHVip[i];
    fHVin[i]=data.fHVin[i];
  }
  
  for(int i=0;i<kNLV;i++) {
    fLVv[i]=data.fLVv[i];
    fLVi[i]=data.fLVi[i];
  }

  for(int i=0;i<kNLV33;i++) {
    fLVv33[i]=data.fLVv33[i];
    fLVi33[i]=data.fLVi33[i];
  }

  for(int i=0;i<kNLV50;i++) {
    fLVv50[i]=data.fLVv50[i];
    fLVi50[i]=data.fLVi50[i];
  }

  for(int i=0;i<kNLV48;i++) {
    fLVv48[i]=data.fLVv48[i];
    fLVi48[i]=data.fLVi48[i];
  }

  for(int i=0;i<kNFEEthr;i++) {
    fFEEthr[i]=data.fFEEthr[i];
  }

  for(int i=0;i<kNFEEtfeac;i++) {
    fFEEtfeac[i]=data.fFEEtfeac[i];
  }

  for(int i=0;i<kNFEEttrm;i++) {
    fFEEttrm[i]=data.fFEEttrm[i];
  }
  
  for(int i=0;i<3;i++) {
    fT[i]=data.fT[i];
    fP[i]=data.fP[i];
  }
  
}
//---------------------------------------------------------------

AliTOFDataDCS& AliTOFDataDCS:: operator=(const AliTOFDataDCS & data) { 

// assignment operator

  this->fRun=data.GetRun();
  this->fStartTime=data.GetStartTime();
  this->fEndTime=data.GetEndTime();

  for(int i=0;i<kNAliases;i++) {
    this->fAliasNames[i]=data.GetAliasName(i);
  }

  for(int i=0;i<3;i++) {
    this->fT[i]=data.GetT(i);
    this->fP[i]=data.GetP(i);
  }


  for(int i=0;i<kNHV;i++) {
    this->fHVvp[i]=data.GetHVvp(i);
    this->fHVvn[i]=data.GetHVvn(i);
    this->fHVip[i]=data.GetHVip(i);
    this->fHVin[i]=data.GetHVin(i);
  }

  for(int i=0;i<kNLV;i++) {
    this->fLVv[i]=data.GetLVv(i);
    this->fLVi[i]=data.GetLVi(i);
  }

  for(int i=0;i<kNLV33;i++) {
    this->fLVv33[i]=data.GetLVv33(i);
    this->fLVi33[i]=data.GetLVi33(i);
  }

  for(int i=0;i<kNLV50;i++) {
    this->fLVv50[i]=data.GetLVv50(i);
    this->fLVi50[i]=data.GetLVi50(i);
  }

  for(int i=0;i<kNLV48;i++) {
    this->fLVv48[i]=data.GetLVv48(i);
    this->fLVi48[i]=data.GetLVi48(i);
  }

  for(int i=0;i<kNFEEthr;i++) {
    this->fFEEthr[i]=data.GetFEEthr(i);
  }

  for(int i=0;i<kNFEEtfeac;i++) {
    this->fFEEtfeac[i]=data.GetFEEtfeac(i);
  }

  for(int i=0;i<kNFEEttrm;i++) {
    this->fFEEttrm[i]=data.GetFEEttrm(i);
  }

  this->fIsProcessed=data.fIsProcessed;

  return *this;
}
//---------------------------------------------------------------
AliTOFDataDCS::~AliTOFDataDCS() {

  // destructor

  for(int i=0;i<kNHV;i++) {
    delete fHVvp[i];
    fHVvp[i]=0;
    delete fHVvn[i];
    fHVvn[i]=0;
    delete fHVip[i];
    fHVip[i]=0;
    delete fHVin[i];
    fHVin[i]=0;
  }
  
  for(int i=0;i<kNLV;i++) {
    delete fLVv[i];
    fLVv[i]=0;
    delete fLVi[i];
    fLVi[i]=0;
  }
  
  for(int i=0;i<kNLV33;i++) {
    delete fLVv33[i];
    fLVv33[i]=0;
    delete fLVi33[i];
    fLVi33[i]=0;
  }
  
  for(int i=0;i<kNLV50;i++) {
    delete fLVv50[i];
    fLVv50[i]=0;
    delete fLVi50[i];
    fLVi50[i]=0;
  }
  
  for(int i=0;i<kNLV48;i++) {
    delete fLVv48[i];
    fLVv48[i]=0;
    delete fLVi48[i];
    fLVi48[i]=0;
  }
  
  for(int i=0;i<kNFEEthr;i++) {
    delete fFEEthr[i];
    fFEEthr[i]=0;
  }
  
  for(int i=0;i<kNFEEtfeac;i++) {
    delete fFEEtfeac[i];
    fFEEtfeac[i]=0;
  }
  
  for(int i=0;i<kNFEEttrm;i++) {
    delete fFEEttrm[i];
    fFEEttrm[i]=0;
  }
}

//-----------------------------------------------------------------------------
Float_t* AliTOFDataDCS::GetT()const {

  // method to retrieve environment temperature info

  Float_t* t=0;
  for (Int_t i=0;i<3;i++){
    t[i]=this->fT[i];
  }
  return t;
}
//-----------------------------------------------------------------------------
Float_t* AliTOFDataDCS::GetP() const{

  // method to retrieve environment pressure info

  Float_t* p=0;
  for (Int_t i=0;i<3;i++){
    p[i]=this->fP[i];
  }
  return p;
}

//---------------------------------------------------------------
Bool_t AliTOFDataDCS::ProcessData(TMap& aliasMap){

  if(!(fAliasNames[0])) Init();

  Float_t timeMin = (Float_t)fStartTime;
  Float_t timeMax = (Float_t)fEndTime;
  Float_t val=0;
  Float_t val1=0;
  Float_t time=0; 
  Float_t delta[2];
  Float_t timedelta[2];

  TObjArray *aliasArr;
  AliDCSValue* aValue;
  AliDCSValue* aValue1;
  TH1F * histoT=0x0;
  TH1F * histoP=0x0;

  // starting loop on aliases
  for(int j=0; j<kNAliases; j++){
    for (Int_t k=0;k<2;k++) {
      delta[k]=0;
      timedelta[k]=0;
    }
    //AliInfo(Form("j = %i, with alias = %s",j,fAliasNames[j].Data()));
    aliasArr = (TObjArray*) aliasMap.GetValue(fAliasNames[j].Data());
    if(!aliasArr){
      AliError(Form("Alias %s not found!", fAliasNames[j].Data()));
      return kFALSE;
    }

    Introduce(j, aliasArr);
    
    if(aliasArr->GetEntries()<3){
      AliError(Form("Alias %s has just %d entries!",
		    fAliasNames[j].Data(),aliasArr->GetEntries()));
      continue;
    }
    
    TIter iterarray(aliasArr);
    
    Int_t nentries = aliasArr->GetEntries();
    Int_t deltaTimeStamp = (Int_t) nentries/3;
    Int_t deltaTimeStamp1 = (Int_t) nentries/2;
    AliDCSValue *lastDCSvalue = (AliDCSValue*) aliasArr->At(nentries-1);
    Float_t maxTimeStamp = (Float_t) (lastDCSvalue->GetTimeStamp());
    Float_t minTimeStamp = 0;

    // filling aliases with 10 floats+1 Usign
    if (j < kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm){
      Int_t index = 0;
      for (Int_t k=0;k<3;k++){
	index = deltaTimeStamp*k;
	if (k==0) {
	  index=0;
	}
	else if (k==1) {
	  index=deltaTimeStamp1;
	} 
	else if (k==2) {
	  index=nentries-1; 
	}
	aValue = (AliDCSValue*) aliasArr->At(index);
	val = aValue->GetFloat();
	time = (Float_t) (aValue->GetTimeStamp());
	if (j<kNHV){
	  fHVvp[j]->SetFloat(k,val);
	  fHVvp[j]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*2){
	  fHVvn[j-kNHV]->SetFloat(k,val);
	  fHVvn[j-kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*3){
	  fHVip[j-2*kNHV]->SetFloat(k,val);
	  fHVip[j-2*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4){
	  fHVin[j-3*kNHV]->SetFloat(k,val);
	  fHVin[j-3*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV){
	  fLVv[j-4*kNHV]->SetFloat(k,val);
	  fLVv[j-4*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2){
	  fLVi[j-4*kNHV-kNLV]->SetFloat(k,val);
	  fLVi[j-4*kNHV-kNLV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33){
	  fLVv33[j-4*kNHV-2*kNLV]->SetFloat(k,val);
	  fLVv33[j-4*kNHV-2*kNLV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2){
	  fLVi33[j-4*kNHV-2*kNLV-kNLV33]->SetFloat(k,val);
	  fLVi33[j-4*kNHV-2*kNLV-kNLV33]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50){
	  fLVv50[j-4*kNHV-2*kNLV-2*kNLV33]->SetFloat(k,val);
	  fLVv50[j-4*kNHV-2*kNLV-2*kNLV33]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2){
	  fLVi50[j-4*kNHV-2*kNLV-2*kNLV33-kNLV50]->SetFloat(k,val);
	  fLVi50[j-4*kNHV-2*kNLV-2*kNLV33-kNLV50]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48){
	  fLVv48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50]->SetFloat(k,val);
	  fLVv48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2){
	  fLVi48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48]->SetFloat(k,val);
	  fLVi48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr){
	  fFEEthr[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48]->SetFloat(k,val);
	  fFEEthr[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac){
	  fFEEtfeac[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr]->SetFloat(k,val);
	  fFEEtfeac[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr]->SetTimeStampFloat(k,time);
	}
	else {
	  fFEEttrm[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac]->SetFloat(k,val);
	  fFEEttrm[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac]->SetTimeStampFloat(k,time);
	}
      }
    }
  
    //filling Temperature and Pressure aliases
    
    else {
      Int_t entriesT=0;
      Int_t entriesP=0;
      if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+1){
	histoT=new TH1F("histoT","histoT",nentries,minTimeStamp,maxTimeStamp);
      }
      else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+2) {
	histoP=new TH1F("histoP","histoP",nentries,minTimeStamp,maxTimeStamp);
      }
      while ((aValue = (AliDCSValue*) iterarray.Next())) {
	val = aValue->GetFloat();
	time = (Float_t) (aValue->GetTimeStamp());
	if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+1){
	  histoT->Fill(time,val);
	  entriesT = (Int_t)(histoT->GetEntries());
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+2){
	  histoP->Fill(time,val);
	  entriesP = (Int_t)(histoP->GetEntries());
	}
      }
    
      if (j==kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+1){
	entriesT = (Int_t)(histoT->GetEntries());
	histoT->Fit("pol1","Q","");
	histoP->Fit("pol1","Q","");
      
	TF1 *tempFunc = histoT->GetFunction("pol1");
	TF1 *pressFunc = histoP->GetFunction("pol1");
      
	SetInterceptT((Float_t)tempFunc->GetParameter(0));
	SetSlopeT((Float_t)tempFunc->GetParameter(1));
	SetMaxT((Float_t)histoT->GetMaximum());
	SetInterceptP((Float_t)pressFunc->GetParameter(0));
	SetSlopeP((Float_t)pressFunc->GetParameter(1));
	SetMaxP((Float_t)histoP->GetMaximum());
      
	TCanvas *chT;
	TString canvasHistoNameT="HistosT";
	chT=new TCanvas(canvasHistoNameT,canvasHistoNameT,20,20,600,600);
	chT->cd();
	histoT->Draw();
	TCanvas *chP;
	TString canvasHistoNameP="HistosP";
	chP=new TCanvas(canvasHistoNameP,canvasHistoNameP,20,20,600,600);
	chP->cd();
	histoP->Draw();
      }
    }
 
    //computing the most significant variations

    Int_t deltamin = (Int_t)(60/(timeMax-timeMin)*nentries);
    Int_t klast = nentries-deltamin;
    
    for (Int_t k=0;k<klast;k++){
      aValue = (AliDCSValue*) aliasArr->At(k);
      aValue1 = (AliDCSValue*) aliasArr->At(k+deltamin);
      val = aValue->GetFloat();
      val1 = aValue1->GetFloat();
      if (delta[0]<=TMath::Abs(val1-val)) {
	delta[0]=TMath::Abs(val1-val);
	timedelta[0] = (Float_t)k;
      }
      if (delta[1]<=delta[0]) {
	Float_t temp = delta[1];
	Float_t timetemp = timedelta[1];
	delta[1]=delta[0];
	delta[0]=temp;
	timedelta[1]=timedelta[0];
	timedelta[0]=timetemp;
      }
    }
    
    for (Int_t kk=0;kk<2;kk++){
      if (j < kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac){
	if (j<kNHV){
	  fHVvp[j]->SetDelta(kk,delta[kk]);
	  fHVvp[j]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*2){
	  fHVvn[j-kNHV]->SetDelta(kk,delta[kk]);
	  fHVvn[j-kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*3){
	  fHVip[j-2*kNHV]->SetDelta(kk,delta[kk]);
	  fHVip[j-2*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4){
	  fHVin[j-3*kNHV]->SetDelta(kk,delta[kk]);
	  fHVin[j-3*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV){
	  fLVv[j-4*kNHV]->SetDelta(kk,delta[kk]);
	  fLVv[j-4*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2){
	  fLVi[j-4*kNHV-kNLV]->SetDelta(kk,delta[kk]);
	  fLVi[j-4*kNHV-kNLV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33){
	  fLVv33[j-4*kNHV-2*kNLV]->SetDelta(kk,delta[kk]);
	  fLVv33[j-4*kNHV-2*kNLV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2){
	  fLVi33[j-4*kNHV-2*kNLV-kNLV33]->SetDelta(kk,delta[kk]);
	  fLVi33[j-4*kNHV-2*kNLV-kNLV33]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50){
	  fLVv50[j-4*kNHV-2*kNLV-2*kNLV33]->SetDelta(kk,delta[kk]);
	  fLVv50[j-4*kNHV-2*kNLV-2*kNLV33]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2){
	  fLVi50[j-4*kNHV-2*kNLV-2*kNLV33-kNLV50]->SetDelta(kk,delta[kk]);
	  fLVi50[j-4*kNHV-2*kNLV-2*kNLV33-kNLV50]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48){
	  fLVv48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50]->SetDelta(kk,delta[kk]);
	  fLVv48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2){
	  fLVi48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48]->SetDelta(kk,delta[kk]);
	  fLVi48[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr){
	  fFEEthr[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48]->SetDelta(kk,delta[kk]);
	  fFEEthr[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac){
	  fFEEtfeac[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr]->SetDelta(kk,delta[kk]);
	  fFEEtfeac[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm){
	  fFEEttrm[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac]->SetDelta(kk,delta[kk]);
	  fFEEttrm[j-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
      }
       
      
    //filling for temperature and pressure
    
      else if (j==kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm){
	fT[2]=delta[1];
      }
      else if (j==kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+1){
	fP[2]=delta[1];
      }
      
    }
  }
    
  fIsProcessed=kTRUE;

  return kTRUE;
}

//---------------------------------------------------------------
void AliTOFDataDCS::Init(){

  // initialization of aliases and DCS data

  TString sindex;
  for(int i=0;i<kNAliases;i++){
    //HV, v
    if (i<kNHV){
	fAliasNames[i] = "tof_hv_vp_";
	sindex.Form("%02i",i);
	fAliasNames[i] += sindex;
	fHVvp[i] = new AliTOFFormatDCS();
    }
    else if (i<kNHV*2){
	fAliasNames[i] = "tof_hv_vn_";
	sindex.Form("%02i",i-kNHV);
	fAliasNames[i] += sindex;
	fHVvn[i-kNHV] = new AliTOFFormatDCS();
    }
    //HV, i
    else if (i<kNHV*3){
	fAliasNames[i] = "tof_hv_ip_";
	sindex.Form("%02i",i-2*kNHV);
	fAliasNames[i] += sindex;
	fHVip[i-2*kNHV] = new AliTOFFormatDCS();
    }
    else if (i<kNHV*4){
	fAliasNames[i] = "tof_hv_in_";
	sindex.Form("%02i",i-3*kNHV);
	fAliasNames[i] += sindex;
	fHVin[i-3*kNHV] = new AliTOFFormatDCS();
    }
    //LV, v
    else if (i<(kNHV*4+kNLV)){
	fAliasNames[i] = "tof_lv_vfea_";
	sindex.Form("%03i",i-4*kNHV);
	fAliasNames[i] += sindex;
	fLVv[i-4*kNHV] = new AliTOFFormatDCS();
    }
    //LV, i
    else if (i<(kNHV*4+kNLV*2)){
	fAliasNames[i] = "tof_lv_ifea_";
	sindex.Form("%03i",i-4*kNHV-kNLV);
	fAliasNames[i] += sindex;
	fLVi[i-4*kNHV-kNLV] = new AliTOFFormatDCS();
    }
    //LV 3.3, v
    else if (i<(kNHV*4+kNLV*2+kNLV33)){
	fAliasNames[i] = "tof_lv_v33_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV);
	fAliasNames[i] += sindex;
	fLVv33[i-4*kNHV-2*kNLV] = new AliTOFFormatDCS();
    }
    //LV 3.3, i
    else if (i<(kNHV*4+kNLV*2+kNLV33*2)){
	fAliasNames[i] = "tof_lv_i33_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV-kNLV33);
	fAliasNames[i] += sindex;
	fLVi33[i-4*kNHV-2*kNLV-kNLV33] = new AliTOFFormatDCS();
    }
    //LV 5.0, v
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50)){
	fAliasNames[i] = "tof_lv_v50_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV-2*kNLV33);
	fAliasNames[i] += sindex;
	fLVv50[i-4*kNHV-2*kNLV-2*kNLV33] = new AliTOFFormatDCS();
    }
    //LV 5.0, i
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2)){
	fAliasNames[i] = "tof_lv_i50_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV-2*kNLV33-kNLV50);
	fAliasNames[i] += sindex;
	fLVi50[i-4*kNHV-2*kNLV-2*kNLV33-kNLV50] = new AliTOFFormatDCS();
    }
    //LV 48, v
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48)){
	fAliasNames[i] = "tof_lv_v48_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50);
	fAliasNames[i] += sindex;
	fLVv48[i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50] = new AliTOFFormatDCS();
    }
    //LV 48, i
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2)){
	fAliasNames[i] = "tof_lv_i48_";
	sindex.Form("%02i",i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48);
	fAliasNames[i] += sindex;
	fLVi48[i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-kNLV48] = new AliTOFFormatDCS();
    }
    //FEE thresholds
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr)){
	fAliasNames[i] = "tof_fee_th_";
	sindex.Form("%04i",i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48);
	fAliasNames[i] += sindex;
 	fFEEthr[i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48] = new AliTOFFormatDCS();
    }
    //FEE FEAC temperatures
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac)){
	fAliasNames[i] = "tof_fee_tfeac_";
	sindex.Form("%03i",i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr);
	fAliasNames[i] += sindex;
 	fFEEtfeac[i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr] = new AliTOFFormatDCS();
    }
    //FEE trms temperatures
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm)){
      AliInfo(Form("**before temperature "));
	fAliasNames[i] = "tof_fee_ttrm_";
	sindex.Form("%04i",i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac);
	fAliasNames[i] += sindex;
 	fFEEttrm[i-4*kNHV-2*kNLV-2*kNLV33-2*kNLV50-2*kNLV48-kNFEEthr-kNFEEtfeac] = new AliTOFFormatDCS();
    }
    //environment temperature
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+1)){
      AliInfo(Form("**temperature "));

	fAliasNames[i] = "temperature";
    }
    //environment pressure
    else if (i<(kNHV*4+kNLV*2+kNLV33*2+kNLV50*2+kNLV48*2+kNFEEthr+kNFEEtfeac+kNFEEttrm+2)){
	fAliasNames[i] = "pressure";
    }
  }
}


//---------------------------------------------------------------
void AliTOFDataDCS::Introduce(UInt_t numAlias, const TObjArray* aliasArr)const
{

  // method to introduce new aliases

  int entries=0;
  entries = aliasArr->GetEntries();
  int nal=0;
  nal=numAlias;
  AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
  AliInfo(Form("    	%d DP values collected",entries));

}

//---------------------------------------------------------------
void AliTOFDataDCS::Draw(const Option_t* /*option*/)
{
// Draw all histos and graphs

  if(!fIsProcessed) return;

  TCanvas *ch;
  TString canvasHistoName="Histos";
  ch=new TCanvas(canvasHistoName,canvasHistoName,20,20,600,600);
  ch->cd();

  // to be implemented

}

