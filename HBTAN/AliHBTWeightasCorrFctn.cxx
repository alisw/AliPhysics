#include "AliHBTWeightasCorrFctn.h"
#include <TH1.h>
#include <TObjArray.h>

///////////////////////////////////////////////////////
//                                                   //
// AliHBTWeightasCorrFctn.h                             //
//                                                   //
// Class for calculating 3D Weightas correlation        //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTWeightasCorrFctn)

     AliHBTWeightasCorrFctn::AliHBTWeightasCorrFctn(const char* name, const char* title):
 AliHBTTwoPairFctn1D(name,title),

     fNum(0x0),
     fDen(0x0),
     fRat(0x0)

{
//ctor
}
     
/******************************************************************/
AliHBTWeightasCorrFctn::AliHBTWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTTwoPairFctn1D(name,title,nbins,maxXval,minXval),


     fNum(new TObjArray()),
     fDen(new TObjArray()),
     fRat(new TObjArray())
{
     SetParams(nbins,maxXval, minXval);
}

/******************************************************************/
AliHBTWeightasCorrFctn::AliHBTWeightasCorrFctn(const AliHBTWeightasCorrFctn& in):
 AliHBTTwoPairFctn1D(in),



     fNum((in.fNum)?(TObjArray*)in.fNum->Clone():0x0),
     fDen((in.fDen)?(TObjArray*)in.fDen->Clone():0x0),
     fRat((in.fRat)?(TObjArray*)in.fRat->Clone():0x0)
 {
//ctor
}

/******************************************************************/

AliHBTWeightasCorrFctn::~AliHBTWeightasCorrFctn()
{
 //dtor

     delete fNum;
     delete fDen;
     delete fRat;
     
}

/******************************************************************/
void AliHBTWeightasCorrFctn::Write()
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

}

//-------------------------------------
void AliHBTWeightasCorrFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
    //Fills the numerator using pair from the same event
     trackpair = CheckPair(trackpair);
     if(partpair == 0x0) return;
     if(trackpair == 0x0) return;
     
     Double_t weight = 1.0;

     int n = fNumberOfIntervals;
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka

     phi=phi*360/(2*TMath::Pi());
     Double_t q=GetValue(trackpair, partpair);
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));
     
     TH1D *num = ((TH1D*)fNum->At(ntv));
     
        
     if ( trackpair && partpair)
     {
	  if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
	       ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
	  {
	       weight=partpair->GetWeight();
//	       Info("ProcessSameEvent","weight=%lf",weight);
	  }
	  num->Fill(q,weight);
     }
}

/****************************************************************/
void AliHBTWeightasCorrFctn::Init()
{
     BuildHistos();
}
/****************************************************************/


void AliHBTWeightasCorrFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
     
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
     phi=phi*360/(2*TMath::Pi());
     Double_t qout=GetValue(trackpair, partpair);

     int n = fNumberOfIntervals;
     
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));

     TH1D &den = *((TH1D*)fDen->At(ntv));
     if ( trackpair && partpair)
     {
          den.Fill(qout);
     }
}


/******************************************************************/


void AliHBTWeightasCorrFctn::SetParams(Int_t nbins, Float_t maxXval, Float_t minXval){
     fnbins=nbins;
     fmaxXval= maxXval;
     fminXval=minXval;
}
TH1* AliHBTWeightasCorrFctn::GetResult()
{
       
     TH1D *den = ((TH1D*)fDen->UncheckedAt(1));
     return den;
 }




ClassImp(AliHBTQOutWeightasCorrFctn)
     
     AliHBTQOutWeightasCorrFctn::AliHBTQOutWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTWeightasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQOutWeightasCorrFctn::BuildHistos()
{
    
     Int_t i;
     int n=GetNumberOfIntervals();

     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *Num;
     TH1D *Den;
     TH1D *Rat;
     
     TString nameNum = "NumOut";
     TString nameDen = "DenOut";
     TString nameRat = "RatOut";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  Num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  Den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  Rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  Num->Sumw2();
	  Den->Sumw2();
	  Rat->Sumw2();
	  
	  Num->Reset();
	  Den->Reset();
	  Rat->Reset();
	  
	  fNum->Add(Num);
	  fDen->Add(Den);
	  fRat->Add(Rat);
	  
	  nameNum = TString("NumOut");
	  nameDen = TString("DenOut");
	  nameRat = TString("RatOut");
	  
     }

     
 }

ClassImp(AliHBTQSideWeightasCorrFctn)
     
     AliHBTQSideWeightasCorrFctn::AliHBTQSideWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTWeightasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQSideWeightasCorrFctn::BuildHistos()
{
    
     Int_t i;
     int n=GetNumberOfIntervals();
     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *Num;
     TH1D *Den;
     TH1D *Rat;
     
     TString nameNum = "NumSide";
     TString nameDen = "DenSide";
     TString nameRat = "RatSide";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  Num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  Den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  Rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  Num->Sumw2();
	  Den->Sumw2();
	  Rat->Sumw2();
	  
	  Num->Reset();
	  Den->Reset();
	  Rat->Reset();
	  
	  fNum->Add(Num);
	  fDen->Add(Den);
	  fRat->Add(Rat);
	  
	  nameNum = TString("NumSide");
	  nameDen = TString("DenSide");
	  nameRat = TString("RatSide");
	  
     }

     
 }

ClassImp(AliHBTQLongWeightasCorrFctn)
     
     AliHBTQLongWeightasCorrFctn::AliHBTQLongWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
AliHBTWeightasCorrFctn(name,title,nbins,maxXval,minXval)

{
//ct0r
}

void AliHBTQLongWeightasCorrFctn::BuildHistos()
{
    
     Int_t i;
     int n=GetNumberOfIntervals();
     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *Num;
     TH1D *Den;
     TH1D *Rat;
     
     TString nameNum = "NumLong";
     TString nameDen = "DenLong";
     TString nameRat = "RatLong";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  Num = new TH1D(nameNum.Data(),nameNum.Data(),nbins,min,max);
	  Den = new TH1D(nameDen.Data(),nameDen.Data(),nbins,min,max);
	  Rat = new TH1D(nameRat.Data(),nameRat.Data(),nbins,min,max);
	  
	  Num->Sumw2();
	  Den->Sumw2();
	  Rat->Sumw2();
	  
	  Num->Reset();
	  Den->Reset();
	  Rat->Reset();
	  
	  fNum->Add(Num);
	  fDen->Add(Den);
	  fRat->Add(Rat);
	  
	  nameNum = TString("NumLong");
	  nameDen = TString("DenLong");
	  nameRat = TString("RatLong");
	  
     }

     
 }
