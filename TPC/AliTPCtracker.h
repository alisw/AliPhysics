#ifndef ALITPCTRACKER_H
#define ALITPCTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                       TPC tracker
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "AliTPCtrack.h"
#include "AliTPCParam.h"

#define kMAXCLUSTER 2500

class TFile;

class AliTPCtracker {
public:
   static Int_t Clusters2Tracks(const AliTPCParam *par, TFile *of);

//**************** Internal tracker class ********************** 
   class AliTPCRow {
   public:
     AliTPCRow() {fN=0;}
     void InsertCluster(const AliTPCcluster *c, UInt_t index);
     operator int() const {return fN;}
     const AliTPCcluster* operator[](Int_t i) const {return fClusters[i];}
     UInt_t GetIndex(Int_t i) const {return fIndex[i];}
     Int_t Find(Double_t y) const; 

   private:
     unsigned fN;                                 //number of clusters 
     const AliTPCcluster *fClusters[kMAXCLUSTER]; //pointers to clusters
     UInt_t fIndex[kMAXCLUSTER];                  //indeces of clusters

   private:
     AliTPCRow(const AliTPCRow& r);            //dummy copy constructor
     AliTPCRow &operator=(const AliTPCRow& r); //dummy assignment operator
   };

//**************** Internal tracker class ********************** 
   class AliTPCSector {
   public:
     AliTPCSector() { fN=0; fRow = 0; }
     virtual ~AliTPCSector() { delete[] fRow; }
     static void SetParam(const AliTPCParam *p) { fgParam=p; }
     AliTPCRow& operator[](Int_t i) const { return *(fRow+i); }
     Int_t GetNRows() const { return fN; }
     virtual Double_t GetX(Int_t l) const = 0;
     virtual Double_t GetMaxY(Int_t l) const = 0;
     virtual Double_t GetAlpha() const = 0;
     virtual Double_t GetAlphaShift() const = 0;
     virtual Int_t GetRowNumber(Double_t x) const = 0;
     virtual Double_t GetPadPitchWidth() const = 0;

   protected:
     unsigned fN;                        //number of pad rows 
     AliTPCRow *fRow;                    //array of pad rows
     static const AliTPCParam *fgParam;  //TPC parameters

   private:
     AliTPCSector(const AliTPCSector &s);           //dummy copy contructor
     AliTPCSector& operator=(const AliTPCSector &s);//dummy assignment operator
   };

//**************** Internal tracker class ********************** 
   class AliTPCSSector : public AliTPCSector {
   public:
      AliTPCSSector();
      Double_t GetX(Int_t l) const { return fgParam->GetPadRowRadiiLow(l); }
      Double_t GetMaxY(Int_t l) const { 
         return GetX(l)*TMath::Tan(0.5*GetAlpha()); 
      }
      Double_t GetAlpha() const {return fgParam->GetInnerAngle();}
      Double_t GetAlphaShift() const {return fgParam->GetInnerAngleShift();}
      Double_t GetPadPitchWidth() const {
         return fgParam->GetInnerPadPitchWidth();
      }
      Int_t GetRowNumber(Double_t x) const;
   };

//**************** Internal tracker class ********************** 
   class AliTPCLSector : public AliTPCSector {
   public:
      AliTPCLSector();
      Double_t GetX(Int_t l) const { return fgParam->GetPadRowRadiiUp(l); }
      Double_t GetMaxY(Int_t l) const { 
         return GetX(l)*TMath::Tan(0.5*GetAlpha()); 
      }
      Double_t GetAlpha() const {return fgParam->GetOuterAngle();}
      Double_t GetAlphaShift() const {return fgParam->GetOuterAngleShift();}
      Double_t GetPadPitchWidth() const {
         return fgParam->GetOuterPadPitchWidth();
      }
      Int_t GetRowNumber(Double_t x) const;
   };

//**************** Internal tracker class ********************** 
   class AliTPCseed : public AliTPCtrack {
   public:
     AliTPCseed():AliTPCtrack(){}
     AliTPCseed(UInt_t index, const Double_t xx[5], 
                const Double_t cc[15], Double_t xr, Double_t alpha): 
                AliTPCtrack(index, xx, cc, xr, alpha) {}
     void SetSampledEdx(Float_t q, Int_t i) {
        Double_t c=GetC(), e=GetEta(), t=GetTgl(), x=GetX();
        q *= TMath::Sqrt((1-(c*x-e)*(c*x-e))/(1+t*t));
        fdEdx[i]=q;
     }
     void UseClusters(AliTPCClustersArray *ca, Int_t n=0);
     void CookdEdx(Double_t low=0.05, Double_t up=0.70);

   private:
     Float_t fdEdx[200]; //array of dE/dx samples 
   };

private:
   static Int_t 
   FindProlongation(AliTPCseed& t,const AliTPCSector *sec,Int_t s,Int_t rf=0);
   static void
   MakeSeeds(TObjArray &sd,const AliTPCSector *s,Int_t max,Int_t i1,Int_t i2);
};

inline AliTPCtracker::AliTPCSSector::AliTPCSSector() {
  //default constructor
  if (!fgParam) {
     fprintf(stderr,"AliTPCSSector: parameters are not set !\n");
     return;
  }
  fN=fgParam->GetNRowLow();
  fRow=new AliTPCRow[fN];
}

inline 
Int_t AliTPCtracker::AliTPCSSector::GetRowNumber(Double_t x) const {
  //return pad row number for this x
  Double_t r=fgParam->GetPadRowRadiiLow(fgParam->GetNRowLow()-1);
  if (x > r) return fgParam->GetNRowLow();
  r=fgParam->GetPadRowRadiiLow(0);
  if (x < r) return -1;
  return Int_t((x-r)/fgParam->GetInnerPadPitchLength() + 0.5);
}

inline AliTPCtracker::AliTPCLSector::AliTPCLSector(){
  //default constructor
  if (!fgParam) {
     fprintf(stderr,"AliTPCLSector: parameters are not set !\n");
     return;
  }
  fN=fgParam->GetNRowUp();
  fRow=new AliTPCRow[fN];
}

inline 
Int_t AliTPCtracker::AliTPCLSector::GetRowNumber(Double_t x) const {
  //return pad row number for this x
  Double_t r=fgParam->GetPadRowRadiiUp(fgParam->GetNRowUp()-1);
  if (x > r) return fgParam->GetNRowUp();
  r=fgParam->GetPadRowRadiiUp(0);
  if (x < r) return -1;
  return Int_t((x-r)/fgParam->GetOuterPadPitchLength() + 0.5);
}

//-----------------------------------------------------------------

#endif


