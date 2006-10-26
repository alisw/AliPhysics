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
	  fHVvpos[i]=0x0;
	  fHVvneg[i]=0x0;
	  fHVcpos[i]=0x0;
	  fHVcneg[i]=0x0;
	}

	for(int i=0;i<kNLV;i++) {
	  fLVv[i]=0x0;
	  fLVc[i]=0x0;
	}

	for(int i=0;i<kNFEEthr;i++) {
	  fFEEthr[i]=0x0;
	}

	for(int i=0;i<kNFEEt;i++) {
	  fFEEt[i]=0x0;
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
    fHVvpos[i]=data.fHVvpos[i];
    fHVvneg[i]=data.fHVvneg[i];
    fHVcpos[i]=data.fHVcpos[i];
    fHVcneg[i]=data.fHVcneg[i];
  }
  
  for(int i=0;i<kNLV;i++) {
    fLVv[i]=data.fLVv[i];
    fLVc[i]=data.fLVc[i];
  }

  for(int i=0;i<kNFEEthr;i++) {
    fFEEthr[i]=data.fFEEthr[i];
  }

  for(int i=0;i<kNFEEt;i++) {
    fFEEt[i]=data.fFEEt[i];
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

  //this->fCal=data.GetCal();

  for(int i=0;i<kNHV;i++) {
    this->fHVvpos[i]=data.GetHVvpos(i);
    this->fHVvneg[i]=data.GetHVvneg(i);
    this->fHVcpos[i]=data.GetHVcpos(i);
    this->fHVcneg[i]=data.GetHVcneg(i);
  }

  for(int i=0;i<kNLV;i++) {
    this->fLVv[i]=data.GetLVv(i);
    this->fLVc[i]=data.GetLVc(i);
  }

  for(int i=0;i<kNFEEthr;i++) {
    this->fFEEthr[i]=data.GetFEEthr(i);
  }

  for(int i=0;i<kNFEEt;i++) {
    this->fFEEt[i]=data.GetFEEt(i);
  }

  this->fIsProcessed=data.fIsProcessed;

  return *this;
}
//---------------------------------------------------------------
AliTOFDataDCS::~AliTOFDataDCS() {

  // destructor

	for(int i=0;i<kNHV;i++) {
	  delete fHVvpos[i];
	  fHVvpos[i]=0;
	  delete fHVvneg[i];
	  fHVvneg[i]=0;
	  delete fHVcpos[i];
	  fHVcpos[i]=0;
	  delete fHVcneg[i];
	  fHVcneg[i]=0;
	}

	for(int i=0;i<kNLV;i++) {
	  delete fLVv[i];
	  fLVv[i]=0;
	  delete fLVc[i];
	  fLVc[i]=0;
	}

	for(int i=0;i<kNFEEthr;i++) {
	  delete fFEEthr[i];
	  fFEEthr[i]=0;
	}

	for(int i=0;i<kNFEEt;i++) {
	  delete fFEEt[i];
	  fFEEt[i]=0;
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
void AliTOFDataDCS::ProcessData(TMap& aliasMap){

  if(!(fAliasNames[0])) Init();

  Float_t timeMin = (Float_t)fStartTime;
  Float_t timeMax = (Float_t)fEndTime;
  Int_t nminutes = (Int_t)((timeMax-timeMin)/60);
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
      continue;
    }

    Introduce(j, aliasArr);
    
    if(aliasArr->GetEntries()<3){
      AliError(Form("Alias %s has just %d entries!",
		    fAliasNames[j].Data(),aliasArr->GetEntries()));
      continue;
    }
    
    TIter iterarray(aliasArr);
    
    Int_t nentries = aliasArr->GetEntries();
    //AliInfo(Form("entries = %i",nentries));
    Int_t deltaTimeStamp = (Int_t) nentries/3;
    Int_t deltaTimeStamp1 = (Int_t) nentries/2;
    AliDCSValue *lastDCSvalue = (AliDCSValue*) aliasArr->At(nentries-1);
    Float_t maxTimeStamp = (Float_t) (lastDCSvalue->GetTimeStamp());
    Float_t minTimeStamp = 0;

    // filling aliases with 10 floats+1 Usign
    if (j < kNHV*4+kNLV*2+kNFEEthr+kNFEEt){
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
	  fHVvpos[j]->SetFloat(k,val);
	  fHVvpos[j]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*2){
	  fHVvneg[j-kNHV]->SetFloat(k,val);
	  fHVvneg[j-kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*3){
	  fHVcpos[j-2*kNHV]->SetFloat(k,val);
	  fHVcpos[j-2*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4){
	  fHVcneg[j-3*kNHV]->SetFloat(k,val);
	  fHVcneg[j-3*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV){
	  fLVv[j-4*kNHV]->SetFloat(k,val);
	  fLVv[j-4*kNHV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2){
	  fLVc[j-4*kNHV-kNLV]->SetFloat(k,val);
	  fLVc[j-4*kNHV-kNLV]->SetTimeStampFloat(k,time);
	}
	else if (j<kNHV*4+kNLV*2+kNFEEthr){
	  fFEEthr[j-4*kNHV-2*kNLV]->SetFloat(k,val);
	  fFEEthr[j-4*kNHV-2*kNLV]->SetTimeStampFloat(k,time);
	}
	else {
	  fFEEt[j-4*kNHV-2*kNLV-kNFEEthr]->SetFloat(k,val);
	  fFEEt[j-4*kNHV-2*kNLV-kNFEEthr]->SetTimeStampFloat(k,time);
	}
      }
    }
  
    //filling Temperature and Pressure aliases
    
    else {
      Int_t entriesT=0;
      Int_t entriesP=0;
      if (j<kNHV*4+kNLV*2+kNFEEthr+kNFEEt+1){
	histoT=new TH1F("histoT","histoT",nentries,minTimeStamp,maxTimeStamp);
      }
      else if (j<kNHV*4+kNLV*2+kNFEEthr+kNFEEt+2) {
	histoP=new TH1F("histoP","histoP",nentries,minTimeStamp,maxTimeStamp);
      }
      while ((aValue = (AliDCSValue*) iterarray.Next())) {
	val = aValue->GetFloat();
	time = (Float_t) (aValue->GetTimeStamp());
	if (j<kNHV*4+kNLV*2+kNFEEthr+kNFEEt+1){
	  histoT->Fill(time,val);
	  entriesT = (Int_t)(histoT->GetEntries());
	}
	else if (j<kNHV*4+kNLV*2+kNFEEthr+kNFEEt+2){
	  histoP->Fill(time,val);
	  entriesP = (Int_t)(histoP->GetEntries());
	}
      }
    
      if (j==kNHV*4+kNLV*2+kNFEEthr+kNFEEt+1){
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
      if (j < kNHV*4+kNLV*2+kNFEEthr){
	if (j<kNHV){
	  fHVvpos[j]->SetDelta(kk,delta[kk]);
	  fHVvpos[j]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*2){
	  fHVvneg[j-kNHV]->SetDelta(kk,delta[kk]);
	  fHVvneg[j-kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*3){
	  fHVcpos[j-2*kNHV]->SetDelta(kk,delta[kk]);
	  fHVcpos[j-2*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4){
	  fHVcneg[j-3*kNHV]->SetDelta(kk,delta[kk]);
	  fHVcneg[j-3*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV){
	  fLVv[j-4*kNHV]->SetDelta(kk,delta[kk]);
	  fLVv[j-4*kNHV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2){
	  fLVc[j-4*kNHV-kNLV]->SetDelta(kk,delta[kk]);
	  fLVc[j-4*kNHV-kNLV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNFEEthr){
	  fFEEthr[j-4*kNHV-2*kNLV]->SetDelta(kk,delta[kk]);
	  fFEEthr[j-4*kNHV-2*kNLV]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (j<kNHV*4+kNLV*2+kNFEEthr+kNFEEt){
	  fFEEt[j-4*kNHV-2*kNLV+kNFEEthr]->SetDelta(kk,delta[kk]);
	  fFEEt[j-4*kNHV-2*kNLV+kNFEEthr]->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
      }
       
      
    //filling for temperature and pressure
    
      else if (j==kNHV*4+kNLV*2+kNFEEthr+kNFEEt){
	fT[2]=delta[1];
      }
      else if (j==kNHV*4+kNLV*2+kNFEEthr+kNFEEt+1){
	fP[2]=delta[1];
      }
      
        }
  }
    
  fIsProcessed=kTRUE;

}

//---------------------------------------------------------------
void AliTOFDataDCS::Init(){

  // initialization of aliases and DCS data

  for(int i=0;i<kNAliases;i++){
    if (i<kNHV){
	fAliasNames[i] = "HVvpos";
	fAliasNames[i] += i;
	//AliInfo(Form("i = %i, alias name = %s ", i, fAliasNames[i].Data())); 
	fHVvpos[i] = new AliTOFFormatDCS();
    }
    else if (i<kNHV*2){
	fAliasNames[i] = "HVvneg";
	fAliasNames[i] += i-kNHV;
	fHVvneg[i-kNHV] = new AliTOFFormatDCS();
    }
    else if (i<kNHV*3){
	fAliasNames[i] = "HVcpos";
	fAliasNames[i] += i-2*kNHV;
	fHVcpos[i-2*kNHV] = new AliTOFFormatDCS();
    }
    else if (i<kNHV*4){
	fAliasNames[i] = "HVcneg";
	fAliasNames[i] += i-3*kNHV;
	fHVcneg[i-3*kNHV] = new AliTOFFormatDCS();
    }
    else if (i<(kNHV*4+kNLV)){
	fAliasNames[i] = "LVv";
	fAliasNames[i] += i-4*kNHV;
	fLVv[i-4*kNHV] = new AliTOFFormatDCS();
    }
    else if (i<(kNHV*4+kNLV*2)){
	fAliasNames[i] = "LVc";
	fAliasNames[i] += i-4*kNHV-kNLV;
	fLVc[i-4*kNHV-kNLV] = new AliTOFFormatDCS();
    }
    else if (i<(kNHV*4+kNLV*2+kNFEEthr)){
	fAliasNames[i] = "FEEthr";
	fAliasNames[i] += i-4*kNHV-2*kNLV;
 	fFEEthr[i-4*kNHV-2*kNLV] = new AliTOFFormatDCS();
    }
    else if (i<(kNHV*4+kNLV*2+kNFEEthr+kNFEEt)){
	fAliasNames[i] = "FEEt";
	fAliasNames[i] += i-4*kNHV-2*kNLV-kNFEEthr;
 	fFEEt[i-4*kNHV-2*kNLV-kNFEEthr] = new AliTOFFormatDCS();
    }
    else if (i<(kNHV*4+kNLV*2+kNFEEthr+kNFEEt+1)){
	fAliasNames[i] = "Temperature";
    }
    else if (i<(kNHV*4+kNLV*2+kNFEEthr+kNFEEt+2)){
	fAliasNames[i] = "Pressure";
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
  //AliInfo(Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
  //AliInfo(Form("    	%d DP values collected",entries));

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

