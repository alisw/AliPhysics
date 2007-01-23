/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

#include <TH1.h>
#include <TObjArray.h>
#include <Riostream.h>

#include "AliHBTasCorrFctn.h"

///////////////////////////////////////////////////////
//                                                   //
// AliHBTasCorrFctn.h                             //
//                                                   //
// Class for calculating 3D as correlation        //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTasCorrFctn)

     AliHBTasCorrFctn::AliHBTasCorrFctn(const char* name, const char* title):
 AliHBTOnePairFctn1D(name,title),

     fNum(0x0),
     fDen(0x0),
     fRat(0x0)

{
//ctor
}
     
/******************************************************************/
AliHBTasCorrFctn::AliHBTasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTOnePairFctn1D(name,title,nbins,maxXval,minXval),


     fNum(new TObjArray()),
     fDen(new TObjArray()),
     fRat(new TObjArray())
{
     SetParams(nbins,maxXval, minXval);
}

/******************************************************************/
AliHBTasCorrFctn::AliHBTasCorrFctn(const AliHBTasCorrFctn& in):
 AliHBTOnePairFctn1D(in),



     fNum((in.fNum)?(TObjArray*)in.fNum->Clone():0x0),
     fDen((in.fDen)?(TObjArray*)in.fDen->Clone():0x0),
     fRat((in.fRat)?(TObjArray*)in.fRat->Clone():0x0)
 {
//ctor
}

/******************************************************************/

AliHBTasCorrFctn::~AliHBTasCorrFctn()
{
 //dtor

     delete fNum;
     delete fDen;
     delete fRat;
     
}

/******************************************************************/
Int_t AliHBTasCorrFctn::Write(const char*,Int_t, Int_t)
{
//out    
     Int_t i;
//     Int_t n=GetNumberOfIntervals();
     Double_t scale;

     for(i=0;i<fNumberOfIntervals;i++){
	  TH1D *num = ((TH1D*)fNum->At(i));
	  TH1D *den = ((TH1D*)fDen->At(i));
	  TH1D &rat = *((TH1D*)fRat->At(i));
	  scale = Scale(num,den);
	  Info("Write():","Scale in interval %d = %lf",i,scale);
	  rat.Divide(num,den,scale);
	  
	  
	  num->Write();
	  den->Write();
	  rat.Write();
     }
     return 0;
}

//-------------------------------------
void AliHBTasCorrFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
    //Fills the numerator using pair from the same event
    pair = CheckPair(pair);
    if(pair == 0x0) return;
  //   int n = GetNumberOfIntervals();
     int n = fNumberOfIntervals;
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(pair->Particle1()->Phi()+pair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
     phi=phi*360/(2*TMath::Pi());
     Double_t q=GetValue(pair);
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));
     
     TH1D *num = ((TH1D*)fNum->At(ntv));

     num->Fill(q);
  
}

/****************************************************************/
void AliHBTasCorrFctn::Init()
{
     BuildHistos();
}
/****************************************************************/


void AliHBTasCorrFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//fills denumerator with data from different events 
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(pair->Particle1()->Phi()+pair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
     phi=phi*360/(2*TMath::Pi());
     Double_t qout=GetValue(pair);

//     int n=GetNumberOfIntervals();   
     int n = fNumberOfIntervals;
     
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));

     TH1D &den = *((TH1D*)fDen->At(ntv));

     
     den.Fill(qout);
}


/******************************************************************/


void AliHBTasCorrFctn::SetParams(Int_t nbins, Float_t maxXval, Float_t minXval)
{
//sets parameters for histograms
   	fnbins=nbins;
     fmaxXval= maxXval;
     fminXval=minXval;
}
TH1* AliHBTasCorrFctn::GetResult()
{
//does nothing, it can as well have nothing inside, its just for class to work normally
     TH1D *den = ((TH1D*)fDen->UncheckedAt(1));
     return den;
 }



ClassImp(AliHBTQOutasCorrFctn)
     
     AliHBTQOutasCorrFctn::AliHBTQOutasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQOutasCorrFctn::BuildHistos()
{
//builds histograms
     Int_t i;
     int n=GetNumberOfIntervals();

     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *num;
     TH1D *den;
     TH1D *rat;
     
     TString nameNum = "NumOut";
     TString nameDen = "DenOut";
     TString nameRat = "RatOut";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  num->Sumw2();
	  den->Sumw2();
	  rat->Sumw2();
	  
	  num->Reset();
	  den->Reset();
	  rat->Reset();
	  
	  fNum->Add(num);
	  fDen->Add(den);
	  fRat->Add(rat);
	  
	  nameNum = TString("NumOut");
	  nameDen = TString("DenOut");
	  nameRat = TString("RatOut");
	  
     }

     
 }

ClassImp(AliHBTQSideasCorrFctn)
     
     AliHBTQSideasCorrFctn::AliHBTQSideasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQSideasCorrFctn::BuildHistos()
{
//builds histograms
     Int_t i;
     int n=GetNumberOfIntervals();
     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *num;
     TH1D *den;
     TH1D *rat;
     
     TString nameNum = "NumSide";
     TString nameDen = "DenSide";
     TString nameRat = "RatSide";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  num->Sumw2();
	  den->Sumw2();
	  rat->Sumw2();
	  
	  num->Reset();
	  den->Reset();
	  rat->Reset();
	  
	  fNum->Add(num);
	  fDen->Add(den);
	  fRat->Add(rat);
	  
	  nameNum = TString("NumSide");
	  nameDen = TString("DenSide");
	  nameRat = TString("RatSide");
	  
     }

     
 }

ClassImp(AliHBTQLongasCorrFctn)
     
     AliHBTQLongasCorrFctn::AliHBTQLongasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQLongasCorrFctn::BuildHistos()
{
//builds histograms
     Int_t i;
     int n=GetNumberOfIntervals();
     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *num;
     TH1D *den;
     TH1D *rat;
     
     TString nameNum = "NumLong";
     TString nameDen = "DenLong";
     TString nameRat = "RatLong";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  num->Sumw2();
	  den->Sumw2();
	  rat->Sumw2();
	  
	  num->Reset();
	  den->Reset();
	  rat->Reset();
	  
	  fNum->Add(num);
	  fDen->Add(den);
	  fRat->Add(rat);
	  
	  nameNum = TString("NumLong");
	  nameDen = TString("DenLong");
	  nameRat = TString("RatLong");
	  
     }

     
 }
