#include "AliHBTWeightasCorrFctn.h"
#include <TH1.h>
#include <Riostream.h>

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
     //process particles from different events 
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
     phi=phi*360/(2*TMath::Pi());
     trackpair = CheckPair(trackpair);
     

     int n = fNumberOfIntervals;
     
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));

     TH1D &den = *((TH1D*)fDen->At(ntv));
     if ( trackpair && partpair)
     {
     Double_t qout=GetValue(trackpair, partpair);	
	     den.Fill(qout);	
     }
}


/******************************************************************/


void AliHBTWeightasCorrFctn::SetParams(Int_t nbins, Float_t maxXval, Float_t minXval)
{
//set histogram parameters
     fnbins=nbins;
     fmaxXval= maxXval;
     fminXval=minXval;
}
TH1* AliHBTWeightasCorrFctn::GetResult()
{
//just for compatibility       
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
//buils histograms (allocates memory etc)
     Int_t i;
     int n=GetNumberOfIntervals();

     int nbins=Getnbins();

     double max = GetmaxXval();
     double min = GetminXval();
     char buff[10];
     

     
     TH1D *num;
     TH1D *den;
     TH1D *rat;
     
     TString nameNum = "WeightNumOut";
     TString nameDen = "WeightDenOut";
     TString nameRat = "WeightRatOut";
     
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
	  
	  nameNum = TString("WeightNumOut");
	  nameDen = TString("WeightDenOut");
	  nameRat = TString("WeightRatOut");
	  
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
     
     
     TString nameNum = "WeightNumSide";
     TString nameDen = "WeightDenSide";
     TString nameRat = "WeightRatSide";
     
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
	  
	  nameNum = TString("WeightNumSide");
	  nameDen = TString("WeightDenSide");
	  nameRat = TString("WeightRatSide");
	  
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
     
     TString nameNum = "WeightNumLong";
     TString nameDen = "WeightDenLong";
     TString nameRat = "WeightRatLong";
     
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
	  
	  nameNum = TString("WeightNumLong");
	  nameDen = TString("WeightDenLong");
	  nameRat = TString("WeightRatLong");
	  
     }

     
 }

/********************************************************************/
ClassImp(AliHBTasWeightQOSLCorrFctn)
     AliHBTasWeightQOSLCorrFctn::AliHBTasWeightQOSLCorrFctn(const char* name, const char* title):
 AliHBTTwoPairFctn3D(name,title),

    fNum(0x0),
     fDen(0x0),
     fRat(0x0)
//	fNumberOfIntervals(1)
    

{
//ctor
}
     
/******************************************************************/
AliHBTasWeightQOSLCorrFctn::AliHBTasWeightQOSLCorrFctn(const char* name, const char* title, Int_t nXbins, Double_t maxXval, Double_t minXval, Int_t nYbins, Double_t maxYval, Double_t minYval,Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(name,title,nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval),


     fNum(new TObjArray()),
     fDen(new TObjArray()),
     fRat(new TObjArray())
//	fNumberOfIntervals(1)
     
{
     SetParams(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval);
}

/******************************************************************/
AliHBTasWeightQOSLCorrFctn::AliHBTasWeightQOSLCorrFctn(const AliHBTasWeightQOSLCorrFctn& in):
 AliHBTTwoPairFctn3D(in),



     fNum((in.fNum)?(TObjArray*)in.fNum->Clone():0x0),
     fDen((in.fDen)?(TObjArray*)in.fDen->Clone():0x0),
     fRat((in.fRat)?(TObjArray*)in.fRat->Clone():0x0)
 {
//ctor
}

/******************************************************************/

AliHBTasWeightQOSLCorrFctn::~AliHBTasWeightQOSLCorrFctn()
{
 //dtor

     delete fNum;
     delete fDen;
     delete fRat;
     
}

/******************************************************************/
void AliHBTasWeightQOSLCorrFctn::Write()
{
//out    
     Int_t i;
//     Int_t n=GetNumberOfIntervals();
    // Double_t scale;

     for(i=0;i<fNumberOfIntervals;i++){
Info("Write()","%d",i);
	     TH3D *num = ((TH3D*)fNum->UncheckedAt(i));
	  TH3D *den = ((TH3D*)fDen->UncheckedAt(i));
	  TH3D &rat = *((TH3D*)fRat->At(i));
	 // scale = Scale(num,den);
//	  Info("Write():","Scale in interval %d = %lf",i,scale);
	  rat.Divide(num,den);//,scale);
		  

	  num->Write();
	  den->Write();
	  rat.Write();
     }

}

//-------------------------------------
void AliHBTasWeightQOSLCorrFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
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
     Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
     Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
     Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
	       
     
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));
     
     TH3D *num = ((TH3D*)fNum->At(ntv));
     
        
     if ( trackpair && partpair)
     {
	  if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
	       ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
	  {
	       weight=partpair->GetWeight();
//	       Info("ProcessSameEvent","weight=%lf",weight);
	  }
	  num->Fill(out,side,lon,weight);
     }
}

/****************************************************************/
void AliHBTasWeightQOSLCorrFctn::Init()
{
     BuildHistos();
}
/****************************************************************/


void AliHBTasWeightQOSLCorrFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//fill denominator     
     Double_t rplane=0.;   //reaction plane angle - 2 B determined
     Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
     phi=phi*360/(2*TMath::Pi());

     Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
     Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
     Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());

     int n = fNumberOfIntervals;
     
     Int_t ntv;
     ntv =  (int)(phi*n/(360.));

     TH3D &den = *((TH3D*)fDen->At(ntv));
     if ( trackpair && partpair)
     {
          den.Fill(out,side,lon);
     }
}


/******************************************************************/


void AliHBTasWeightQOSLCorrFctn::SetParams(Int_t nXbins, Float_t maxXval, Float_t minXval,Int_t nYbins, Float_t maxYval, Float_t minYval,Int_t nZbins, Float_t maxZval, Float_t minZval)
{
//set histogram parametera
     fnXbins=nXbins;
     fmaxXval= maxXval;
     fminXval=minXval;
     fnYbins=nYbins;
     fmaxYval= maxYval;
     fminYval=minYval;

     fnZbins=nZbins;
     fmaxZval= maxZval;
     fminZval=minZval;


}
TH1* AliHBTasWeightQOSLCorrFctn::GetResult()
{
 //just for compatibility      
     TH1D *den = ((TH1D*)fDen->UncheckedAt(1));
     return den;
 }

void AliHBTasWeightQOSLCorrFctn::BuildHistos()
{
    //builds histograms
     Int_t i;
     int n=GetNumberOfIntervals();
     int nXbins=GetnXbins();
     double maxX = GetmaxXval();
     double minX = GetminXval();
     
     int nYbins=GetnYbins();
     double maxY = GetmaxXval();
     double minY = GetminXval();
     int nZbins=GetnZbins();
     double maxZ = GetmaxZval();
     double minZ = GetminZval();
     char buff[10];
        
     TH3D *num;
     TH3D *den;
     TH3D *rat;
     
     TString nameNum = "OSLWeightNum";
     TString nameDen = "OSLWeightDen";
     TString nameRat = "OSLWeightRat";
     
     for(i=0;i<n;i++){
	  
	  sprintf(buff,"%d",i);

	  nameNum +=TString(buff);

	  nameDen +=TString(buff);
	  nameRat +=TString(buff);
	  
	  
	  num = new TH3D(nameNum.Data(),nameNum.Data(),nXbins,minX,maxX,nYbins,minY,maxY,nZbins,minZ,maxZ);
	  den = new TH3D(nameDen.Data(),nameDen.Data(),nXbins,minX,maxX,nYbins,minY,maxY,nZbins,minZ,maxZ);
	  rat = new TH3D(nameRat.Data(),nameRat.Data(),nXbins,minX,maxX,nYbins,minY,maxY,nZbins,minZ,maxZ);
	  
	  num->Sumw2();
	  den->Sumw2();
	  rat->Sumw2();
	  
	  num->Reset();
	  den->Reset();
	  rat->Reset();
	  
	  fNum->Add(num);
	  fDen->Add(den);
	  fRat->Add(rat);
	  
	  nameNum = TString("OSLWeightNumLong");
	  nameDen = TString("OSLWeightDenLong");
	  nameRat = TString("OSLWeightRatLong");
	  
     }
}
Double_t AliHBTasWeightQOSLCorrFctn::Scale(TH3D *num, TH3D *den)
{
//scale histograms method, the default one didn't like me so I've done this one :)
if (AliVAODParticle::GetDebug()>0) Info("Scale","Enetered Scale()");
  if(!num)
     {
         Error("Scale","No numerator");
         return 0.0;
     }
  if(!den)
     {
         Error("Scale","No denominator");
         return 0.0;
     }
                                       
   if( (fNBinsToScaleX < 1) || (fNBinsToScaleY < 1) || (fNBinsToScaleZ < 1))
    {
         Error("Scale","Number of bins for scaling is smaller thnan 1");
         return 0.0;
					    }
    UInt_t nbinsX = num->GetNbinsX();
  if (fNBinsToScaleX > nbinsX)
    {
        Error("Scale","Number of X bins for scaling is bigger thnan number of bins in histograms");
     return 0.0;
    }
		     
    UInt_t nbinsY = num->GetNbinsX();
   if (fNBinsToScaleY > nbinsY)
        {
         Error("Scale","Number of Y bins for scaling is bigger thnan number of bins in histograms");
        return 0.0;
        }

   UInt_t nbinsZ = num->GetNbinsZ();
  if (fNBinsToScaleZ > nbinsZ)
    {
      Error("Scale","Number of Z bins for scaling is bigger thnan number of bins in histograms");
   return 0.0;
     }
		    
  if (AliVAODParticle::GetDebug()>0) Info("Scale","No errors detected");
    Int_t offsetX = nbinsX - fNBinsToScaleX - 1; //bin that we start loop over bins in axis X
    Int_t offsetY = nbinsY - fNBinsToScaleY - 1; //bin that we start loop over bins in axis Y
    Int_t offsetZ = nbinsZ - fNBinsToScaleZ - 1; //bin that we start loop over bins in axis Z
    Double_t densum = 0.0;
    Double_t numsum = 0.0;

   for (UInt_t k = offsetZ; k<nbinsZ; k++)
     for (UInt_t j = offsetY; j<nbinsY; j++)
       for (UInt_t i = offsetX; i<nbinsX; i++)
         {
            if ( num->GetBinContent(i,j,k) > 0.0 )
              {
              densum += den->GetBinContent(i,j,k);
              numsum += num->GetBinContent(i,j,k);
               }
         }
     if(AliVAODParticle::GetDebug() > 0)
     Info("Scale","numsum=%f densum=%f fNBinsToScaleX=%d fNBinsToScaleY=%d fNBinsToScaleZ=%d",
              numsum,densum,fNBinsToScaleX,fNBinsToScaleY,fNBinsToScaleZ);
   if (numsum == 0) return 0.0;
     Double_t ret = densum/numsum;
  if(AliVAODParticle::GetDebug() > 0) Info("Scale","returning %f",ret);
  return ret;
}

