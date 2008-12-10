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
$Log: AliTOFDataDCS.cxx,v $
Revision 1.9  2007/05/04 14:02:45  decaro
AliTOFDataDCS::Draw(Option_t *) method declared const: compiling warning suppression

Revision 1.8  2007/05/03 09:45:09  decaro
Coding convention: RN11 violation -> suppression

Revision 1.7  2007/05/02 14:09:39  arcelli
Retrieval of Env. Temperature removed (will get it from the GRP)

Revision 1.6  2007/04/04 17:19:19  arcelli
Moved some printout to debug level

Revision 1.5  2007/02/20 15:57:00  decaro
Raw data update: to read the TOF raw data defined in UNPACKED mode

Revision 1.4  2007/02/19 15:41:55  decaro
Coding convention: few corrections

Revision 1.3  2007/01/24 11:19:58  arcelli
Modify ProcessData to return a logical (CZ)

Revision 1.2  2006/12/18 18:17:38  arcelli
Updated Aliases for DCS TOF datapoints (C.Zampolli)

Revision 1.1  2006/10/26 09:10:52  arcelli
Class for handling the TOF DCS data in the Shuttle (C.Zampolli)

*/  

#include "TString.h"
//#include "TF1.h"
//#include "TH1F.h"
#include "TTimeStamp.h"
#include "TMap.h"
#include "TCanvas.h"

#include "AliDCSValue.h"
#include "AliLog.h"

#include "AliTOFDataDCS.h"
#include "AliTOFFormatDCS.h"

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
	fStartTimeDCSQuery(0),
	fEndTimeDCSQuery(0),
	fIsProcessed(kFALSE),
	fFDR(kFALSE)
{

  // main constructor 

  for(int i=0;i<kNHV;i++) {
    fHVvp[i]=0x0;
    fHVvn[i]=0x0;
    fHVip[i]=0x0;
    fHVin[i]=0x0;
  }
    
  
}

//---------------------------------------------------------------
AliTOFDataDCS::AliTOFDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery):
	TObject(),
	fRun(nRun),
	fStartTime(startTime),
	fEndTime(endTime),
	fStartTimeDCSQuery(startTimeDCSQuery),
	fEndTimeDCSQuery(endTimeDCSQuery),
	fIsProcessed(kFALSE),
	fFDR(kFALSE)
{

  // constructor with arguments

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s \n\tStartTime DCS Query %s \n\tEndTime DCS Query %s", nRun,
	TTimeStamp(startTime).AsString(),
	TTimeStamp(endTime).AsString(), 
	TTimeStamp(startTimeDCSQuery).AsString(), 
        TTimeStamp(endTimeDCSQuery).AsString()));

	Init();

}

//---------------------------------------------------------------

AliTOFDataDCS::AliTOFDataDCS(const AliTOFDataDCS & data):
  TObject(data), 
  fRun(data.fRun),
  fStartTime(data.fStartTime),
  fEndTime(data.fEndTime),
  fStartTimeDCSQuery(data.fStartTimeDCSQuery),
  fEndTimeDCSQuery(data.fEndTimeDCSQuery),
  fIsProcessed(data.fIsProcessed),
  fFDR(data.fFDR)

{

// copy constructor

  for(int i=0;i<kNAliases;i++) {
    fAliasNames[i]=data.fAliasNames[i];
  }
 
  for(int i=0;i<kNHV;i++) {
    fHVvp[i]=data.fHVvp[i];
    fHVvn[i]=data.fHVvn[i];
    fHVip[i]=data.fHVip[i];
    fHVin[i]=data.fHVin[i];
  }
  
    
}
//---------------------------------------------------------------

AliTOFDataDCS& AliTOFDataDCS:: operator=(const AliTOFDataDCS & data) { 

// assignment operator

  if (this == &data)
    return *this;

  TObject::operator=(data);
  fRun=data.GetRun();
  fStartTime=data.GetStartTime();
  fEndTime=data.GetEndTime();
  fStartTimeDCSQuery=data.GetStartTimeDCSQuery();
  fEndTimeDCSQuery=data.GetEndTimeDCSQuery();

  for(int i=0;i<kNAliases;i++) {
    fAliasNames[i]=data.GetAliasName(i);
  }

  for(int i=0;i<kNHV;i++) {
    fHVvp[i]=data.GetHVvp(i);
    fHVvn[i]=data.GetHVvn(i);
    fHVip[i]=data.GetHVip(i);
    fHVin[i]=data.GetHVin(i);
  }


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
  
}

//---------------------------------------------------------------
Bool_t AliTOFDataDCS::ProcessData(TMap& aliasMap){

  // method to process the data

  if(!(fAliasNames[0])) Init();

  Float_t val=0;
  Float_t val1=0;
  Float_t time=0; 
  Float_t delta[2];
  Float_t timedelta[2];

  AliInfo(Form(" Start Time = %i",fStartTime));
  AliInfo(Form(" End Time = %i",fEndTime));
  AliInfo(Form(" Start Time DCS Query= %i",fStartTimeDCSQuery));
  AliInfo(Form(" End Time DCS Query= %i",fEndTimeDCSQuery));

  if (fEndTime==fStartTime){
    AliError(Form(" Run with null time length: start time = %i = end time = %i",fStartTime,fEndTime));
    return kFALSE;
  }

  TObjArray *aliasArr;
  AliDCSValue* aValue;
  AliDCSValue* aValue1;

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
      if (!fFDR){
	return kFALSE;    // returning only in case we are not in a FDR run
      }
      else {
	continue;
      }
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

    // filling aliases with 10 floats+1 Usign
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
    }
  
    //computing the most significant variations

    //Float_t timeDiff = (Float_t)(fEndTime-fStartTime);
    Float_t timeDiff = (Float_t)(fEndTimeDCSQuery-fStartTimeDCSQuery);
    Int_t deltamin = (Int_t)(60/timeDiff*nentries); //sampling every minute
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
  AliDebug(2,Form("************ Alias: %s **********",fAliasNames[numAlias].Data()));
  AliDebug(2,Form("    	%d DP values collected",entries));

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

