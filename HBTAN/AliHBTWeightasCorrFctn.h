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
     
     Int_t GetNumberOfIntervals(){return fNumberOfIntervals;}
    
     
     TH1*     GetResult();
   
protected:
     
     virtual Double_t GetValue(AliHBTPair* tpair, AliHBTPair* ppair) const = 0;
     virtual void BuildHistos() = 0;
     
     void SetParams(Int_t nbins,Float_t maxXval,Float_t minXval);
     int Getnbins(){ return fnbins;}
     double GetmaxXval(){return fmaxXval;}
     double GetminXval(){return fminXval;}
     
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
     Double_t GetValue(AliHBTPair* pair, AliHBTPair* ppair) const {ppair=0; return pair->GetQOutLCMS();}
     void BuildHistos();
private:   
     ClassDef(AliHBTQOutWeightasCorrFctn,1)
};

class AliHBTQSideWeightasCorrFctn: public AliHBTWeightasCorrFctn{
public:
     AliHBTQSideWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);
     virtual  ~AliHBTQSideWeightasCorrFctn(){};
         
protected:
     Double_t GetValue(AliHBTPair* pair, AliHBTPair* ppair) const {ppair=0;return pair->GetQSideLCMS();}
     void BuildHistos();
private:   
     ClassDef(AliHBTQSideWeightasCorrFctn,1)
};


class AliHBTQLongWeightasCorrFctn: public AliHBTWeightasCorrFctn{
public:
     AliHBTQLongWeightasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);

     virtual  ~AliHBTQLongWeightasCorrFctn(){};
     
protected:
     Double_t GetValue(AliHBTPair* pair,AliHBTPair* ppair) const {ppair=0;return pair->GetQLongLCMS();}
     void BuildHistos();
private:   
     ClassDef(AliHBTQLongWeightasCorrFctn,1)
};


#endif
