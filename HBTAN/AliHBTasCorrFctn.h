#ifndef ALIHBTASCORRFCTN_H
#define ALIHBTASCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTasCorrFctn.h                                //
//                                                   //
// Class for calculating 3D as correlation           //
// functions                                         //
//author: Grzegorz.Galazka@cern.ch                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"

 
class AliHBTasCorrFctn: public AliHBTOnePairFctn1D
{
public:
     AliHBTasCorrFctn(const char* name = "asejdzbitiCF", 
		       const char* title= "as Correlation Function");

     AliHBTasCorrFctn(const char* name, const char* title,
		       Int_t nbins, Float_t maxXval, Float_t minXval);
     AliHBTasCorrFctn(const AliHBTasCorrFctn& in);
     
     virtual ~AliHBTasCorrFctn();
     
     void Init();
     void ProcessSameEventParticles(AliHBTPair* pair);
     void ProcessDiffEventParticles(AliHBTPair* pair);
     void Write();
     
     void SetNumberOfIntervals(Int_t N){fNumberOfIntervals = N;}
     
     Int_t GetNumberOfIntervals(){return fNumberOfIntervals;}
    
     
     TH1*     GetResult();
   
protected:
     
     virtual Double_t GetValue(AliHBTPair* pair) const = 0; 
     virtual void BuildHistos() = 0; 
     int Getnbins(){ return fnbins;}         // this are workarounds for my lame coding
     double GetmaxXval(){return fmaxXval;}   // these methods are uset to build histograms 
     double GetminXval(){return fminXval;}   // with set by user number of bins etc. 
     void SetParams(Int_t nbins,Float_t maxXval,Float_t minXval); //
     
     TObjArray* fNum; // numerators array
     TObjArray* fDen; // denominators array
     TObjArray* fRat;// correl. fnctns array
     
     
private:
     int fnbins; //number of bins in histograms
     Int_t fNumberOfIntervals;   //number of intervals
     double fmaxXval; //max histogram's X value 
     double fminXval; //min histogram's X value
     
     ClassDef(AliHBTasCorrFctn,1)
};

class AliHBTQOutasCorrFctn: public AliHBTasCorrFctn{
public:
     AliHBTQOutasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);
     
     virtual  ~AliHBTQOutasCorrFctn(){};
    
     
protected:
     Double_t GetValue(AliHBTPair* pair) const {return pair->GetQOutLCMS();}
     void BuildHistos();
private:   
     ClassDef(AliHBTQOutasCorrFctn,1)
};

class AliHBTQSideasCorrFctn: public AliHBTasCorrFctn{
public:
     AliHBTQSideasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);
//     AliHBTSideasCorrFctn(const AliHBTasCorrFctn& in);
     
     virtual  ~AliHBTQSideasCorrFctn(){};
         
protected:
     Double_t GetValue(AliHBTPair* pair) const {return pair->GetQSideLCMS();} 
     void BuildHistos();
private:   
     ClassDef(AliHBTQSideasCorrFctn,1)
};


class AliHBTQLongasCorrFctn: public AliHBTasCorrFctn{
public:
     AliHBTQLongasCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval);

     virtual  ~AliHBTQLongasCorrFctn(){};
     
protected:
     Double_t GetValue(AliHBTPair* pair) const {return pair->GetQLongLCMS();}
     void BuildHistos();
private:   
     ClassDef(AliHBTQLongasCorrFctn,1)
};


#endif
