#ifndef ALIHBTWEIGHTASCORRFCTN_H
#define ALIHBTWEIGHTASCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTWeightasCorrFctn.h                           //
//                                                   //
// Class for calculating 3D Weightas correlation       //
// functions                                         //
// author: Grzegorz.Galazka@cern.ch                  //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"

 
class AliHBTWeightasCorrFctn: public AliHBTTwoPairFctn1D
{
public:
     AliHBTWeightasCorrFctn(const char* name = "asejdzbitiCF", 
		       const char* title= "as Correlation Function");

     AliHBTWeightasCorrFctn(const char* name, const char* title,
		       Int_t nbins, Float_t maxXval, Float_t minXval);
     AliHBTWeightasCorrFctn(const AliHBTWeightasCorrFctn& in);
     
     virtual ~AliHBTWeightasCorrFctn();
     
     void Init();
     void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
     void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
    
     void Write();
     
     void SetNumberOfIntervals(Int_t N){fNumberOfIntervals = N;}
     
     Int_t GetNumberOfIntervals() const {return fNumberOfIntervals;}
    
     
     TH1*     GetResult();
   
protected:
     
     virtual Double_t GetValue(AliHBTPair* tpair, AliHBTPair* ppair) const = 0;
     virtual void BuildHistos() = 0;
     
     void SetParams(Int_t nbins,Float_t maxXval,Float_t minXval);
     int Getnbins() const { return fnbins;}
     double GetmaxXval() const {return fmaxXval;}
     double GetminXval() const {return fminXval;}
     
     TObjArray* fNum; // numerators array
     TObjArray* fDen; // denominators array
     TObjArray* fRat;// correl. fnctns array
     
     
private:
     int fnbins; //number of bins in histograms
     Int_t fNumberOfIntervals;   //number of intervals
     double fmaxXval; //max histogram's X value 
     double fminXval; //min histogram's X value
     
     ClassDef(AliHBTWeightasCorrFctn,1)
};

class AliHBTQOutWeightasCorrFctn: public AliHBTWeightasCorrFctn{
public:
     AliHBTQOutWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);
     
     virtual  ~AliHBTQOutWeightasCorrFctn(){};
    
     
protected:
     Double_t GetValue(AliHBTPair* pair, AliHBTPair* ppair) const {ppair=0; return TMath::Abs(pair->GetQOutLCMS());}
     void BuildHistos();
private:   
     ClassDef(AliHBTQOutWeightasCorrFctn,1)
};

class AliHBTQSideWeightasCorrFctn: public AliHBTWeightasCorrFctn{
public:
     AliHBTQSideWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);
     virtual  ~AliHBTQSideWeightasCorrFctn(){};
         
protected:
     Double_t GetValue(AliHBTPair* pair, AliHBTPair* ppair) const {ppair=0;return TMath::Abs(pair->GetQSideLCMS());}
     void BuildHistos();
private:   
     ClassDef(AliHBTQSideWeightasCorrFctn,1)
};


class AliHBTQLongWeightasCorrFctn: public AliHBTWeightasCorrFctn{
public:
     AliHBTQLongWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);

     virtual  ~AliHBTQLongWeightasCorrFctn(){};
     
protected:
     Double_t GetValue(AliHBTPair* pair,AliHBTPair* ppair) const {ppair=0;return TMath::Abs(pair->GetQLongLCMS());}
     void BuildHistos();
private:   
     ClassDef(AliHBTQLongWeightasCorrFctn,1)
};

/********************************************************************************************************/
class AliHBTasWeightQOSLCorrFctn: public AliHBTTwoPairFctn3D
{
	public:
     AliHBTasWeightQOSLCorrFctn(const char* name = "asejdzbiti3dCF", 
		       const char* title= "as 3d Correlation Function");

     AliHBTasWeightQOSLCorrFctn(const char* name, const char* title,Int_t nXbins , Double_t maxXval , Double_t minXval ,
		     Int_t nYbins, Double_t maxYval, Double_t minYval ,
		     Int_t nZbins, Double_t maxZval, Double_t minZval);
     
     AliHBTasWeightQOSLCorrFctn(const AliHBTasWeightQOSLCorrFctn& in);
     
     virtual ~AliHBTasWeightQOSLCorrFctn();
     
     void Init();
     void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
     void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
    
     void Write();
	     
    void SetNumberOfIntervals(Int_t N){fNumberOfIntervals = N;}
     
     Int_t GetNumberOfIntervals(){return fNumberOfIntervals;}
     
     
     TH1*     GetResult();
   
protected:
     
     Double_t GetValue(AliHBTPair* tpair, AliHBTPair* ppair);
     void BuildHistos();
     void GetValues(AliHBTPair*, AliHBTPair*, Double_t&,Double_t&, Double_t&) const {};
     void SetParams(Int_t nXbins,Float_t maxXval,Float_t minXval,Int_t nYbins,Float_t maxYval,Float_t minYval,Int_t nZbins,Float_t maxZval,Float_t minZval);
     
     int GetnXbins(){ return fnXbins;}
     double GetmaxXval(){return fmaxXval;}
     double GetminXval(){return fminXval;}
     
     int GetnYbins(){ return fnYbins;}
     double GetmaxYval(){return fmaxYval;}
     double GetminYval(){return fminYval;}
     
     int GetnZbins(){ return fnZbins;}
     double GetmaxZval(){return fmaxZval;}
     double GetminZval(){return fminZval;}
     
Double_t Scale(TH3D* num, TH3D *den);	

     TObjArray* fNum; // numerators array
     TObjArray* fDen; // denominators array
     TObjArray* fRat;// correl. fnctns array
     
     
private:
     int fnXbins; //number of bins in histograms
     Int_t fNumberOfIntervals;   //number of intervals
     double fmaxXval; //max histogram's X value 
     double fminXval; //min histogram's X value
     int fnYbins;
     int fnZbins;
      double fmaxYval; //max histogram's X value 
     double fminYval; //min histogram's X value
  double fmaxZval; //max histogram's X value 
     double fminZval; //min histogram's X value
 
     
     ClassDef(AliHBTasWeightQOSLCorrFctn,1)
};


#endif
