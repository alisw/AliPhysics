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
#include "AliTracker.h"
#include "AliTPCtrack.h"
#include "AliTPCreco.h"
#include "AliTPCcluster.h"

class TFile;
class TTree;
class TBranch;
class AliTPCParam;
class TObjArray;
class TClonesArray;
class AliESDEvent;

class AliTPCtracker : public AliTracker {
public:
   AliTPCtracker();
   AliTPCtracker(const AliTPCParam *par);
  ~AliTPCtracker();

   Int_t ReadSeeds(const TFile *in);

   Int_t LoadClusters(TTree *c);
   void UnloadClusters();

   AliCluster *GetCluster(Int_t index) const;
   Int_t Clusters2Tracks(AliESDEvent *event);
   Int_t PropagateBack(AliESDEvent *event);
   Int_t RefitInward(AliESDEvent *event);

   virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 

public:
//**************** Internal tracker class ********************** 
   class AliTPCRow {
   public:
     AliTPCRow():
         fN(0), fSize(kMaxClusterPerRow/8),
         fClusterArray(new AliTPCcluster[fSize]), fX(0) {}
     ~AliTPCRow() {delete[] fClusterArray;}
     void InsertCluster(const AliTPCcluster *c, Int_t sec, Int_t row);
     void ResetClusters() {fN=0; delete[] fClusterArray; fClusterArray=0;}
     operator int() const {return fN;}
     const AliTPCcluster *operator[](Int_t i) const {return fClusters[i];}
     const AliTPCcluster *GetUnsortedCluster(Int_t i) const {
       if ((i < 0) || (i >= fN)) return NULL;
       return fClusterArray+i;
     }
     UInt_t GetIndex(Int_t i) const {return fIndex[i];}
     Int_t Find(Double_t y) const; 
     void SetX(Double_t x) {fX=x;}
     Double_t GetX() const {return fX;}

   private:
     AliTPCRow(const AliTPCRow& r);            //dummy copy constructor
     AliTPCRow &operator=(const AliTPCRow& r); //dummy assignment operator
     Int_t fN;                                          //number of clusters 
     const AliTPCcluster *fClusters[kMaxClusterPerRow]; //pointers to clusters
     Int_t fSize;                                 //size of array of clusters
     AliTPCcluster *fClusterArray;                      //array of clusters
     UInt_t fIndex[kMaxClusterPerRow];                  //indeces of clusters
     Double_t fX;                                 //X-coordinate of this row
   };

//**************** Internal tracker class ********************** 
   class AliTPCSector {
   public:
     AliTPCSector(): 
        fN(0),fRow(0),fAlpha(0),fAlphaShift(0),
        fPadPitchWidth(0),f1PadPitchLength(0),f2PadPitchLength(0){}
    ~AliTPCSector() { delete[] fRow; }
     AliTPCRow& operator[](Int_t i) const { return *(fRow+i); }
     Int_t GetNRows() const { return fN; }
     void Setup(const AliTPCParam *par, Int_t flag);
     Double_t GetX(Int_t l) const {return fRow[l].GetX();}
     Double_t GetMaxY(Int_t l) const {
         return GetX(l)*TMath::Tan(0.5*GetAlpha());
     } 
     Double_t GetAlpha() const {return fAlpha;}
     Double_t GetAlphaShift() const {return fAlphaShift;}
     Int_t GetRowNumber(Double_t x) const {
        //return pad row number for this x
       Double_t r;
       if (fN < 64){
	 r=fRow[fN-1].GetX();
	 if (x > r) return fN;
	 r=fRow[0].GetX();
	 if (x < r) return -1;
	 return Int_t((x-r)/f1PadPitchLength + 0.5);}
       else{	
	   r=fRow[fN-1].GetX();
	   if (x > r) return fN;
	   r=fRow[0].GetX();
	   if (x < r) return -1;
	  Double_t r1=fRow[64].GetX();
	  if(x<r1){	  
	    return Int_t((x-r)/f1PadPitchLength + 0.5);}
	  else{
	    return (Int_t((x-r1)/f2PadPitchLength + 0.5)+64);} 
       }
     }
     Double_t GetPadPitchWidth()  const {return fPadPitchWidth;}

   private:
     Int_t fN;                        //number of pad rows 
     AliTPCRow *fRow;                    //array of pad rows
     Double_t fAlpha;                    //opening angle
     Double_t fAlphaShift;               //shift angle;
     Double_t fPadPitchWidth;            //pad pitch width
     Double_t f1PadPitchLength;          //pad pitch length
     Double_t f2PadPitchLength;          //pad pitch length
   private:
     AliTPCSector(const AliTPCSector &s);           //dummy copy contructor
     AliTPCSector& operator=(const AliTPCSector &s);//dummy assignment operator
   };

//**************** Internal tracker class **********************
   class AliTPCseed : public AliTPCtrack {
   public:
     AliTPCseed():AliTPCtrack(){Init();}
     AliTPCseed(const AliTPCtrack &t):AliTPCtrack(t){Init();}
     AliTPCseed(Double_t xr, Double_t alpha, const Double_t xx[5], 
                const Double_t cc[15], Int_t index): 
                AliTPCtrack(xr, alpha, xx, cc, index) {Init();}
     void SetSampledEdx(Float_t q, Int_t i) {
        Double_t s=GetSnp(), t=GetTgl();
        q *= TMath::Sqrt((1-s*s)/(1+t*t));
        fdEdxSample[i]=q;
     }
     void CookdEdx(Double_t low=0.05, Double_t up=0.70);

   private:
     void Init() {for (Int_t i=0; i<200; i++) fdEdxSample[i]=0.;}
     Float_t fdEdxSample[200]; //array of dE/dx samples 
   };
private:

   void MakeSeeds(Int_t i1, Int_t i2);
   Int_t FollowProlongation(AliTPCseed& t, Int_t rf=0);
   Int_t FollowBackProlongation(AliTPCseed &s, const AliTPCtrack &t);
   Int_t FollowRefitInward(AliTPCseed *seed, AliTPCtrack *track);

   AliTPCtracker(const AliTPCtracker& r);           //dummy copy constructor
   AliTPCtracker &operator=(const AliTPCtracker& r);//dummy assignment operator

   const Int_t fkNIS;        //number of inner sectors
   AliTPCSector *fInnerSec;  //array of inner sectors
   const Int_t fkNOS;        //number of outer sectors
   AliTPCSector *fOuterSec;  //array of outer sectors

   Int_t fN;               //number of "active" sectors
   AliTPCSector *fSectors; //pointer to "active" sectors;
   
   AliTPCParam *fParam;      //! TPC parameters for outer reference plane [SR, GSI, 18.02.2003]
   TObjArray *fSeeds;        //array of track seeds

  ClassDef(AliTPCtracker,2)   // Time Projection Chamber tracker
};

#endif


