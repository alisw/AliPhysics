#ifndef TPC_H
#define TPC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC                     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h" 
#include "AliTPCParam.h"
#include <TMatrix.h>
#include <TTree.h>
#include <TClonesArray.h>

class AliTPCcluster;
class AliTPCtrack;
class AliTPCParam;

class AliTPCDigitsArray;
class AliTPCClustersArray;

class AliTPC : public AliDetector {
protected:
  Int_t          fSens;             // ISENS
  Int_t          fSecAL;            // Upper sector selector
  Int_t          fSecAU;            // Lower sector selector
  Int_t          fSecLows[6];       // List of lower sectors selected
  Int_t          fSecUps[12];       // List of upper sectors selected
  Int_t          fNsectors;         // Number of sectors in TPC
  Int_t          fNclusters;        // Number of clusters in TPC
  Int_t          fNtracks;          // Number of tracks in TPC
  TClonesArray   *fClusters;        // List of clusters for all sectors
  TClonesArray   *fTracks;          // List of reconstructed tracks
  //MI changes
  AliTPCDigitsArray * fDigitsArray;              //detector digit object  
  AliTPCClustersArray * fClustersArray; //detector cluster object
  AliTPCParam *fTPCParam;

  //MK changes

  Float_t        fSide;  // selects left(-1), right(+1), or both(0) sides of the TPC
  Int_t          fNoComp; // number of a drift gas components
  Int_t          fMixtComp[3]; // drift gas components
  Float_t        fMixtProp[3]; // mixture proportions

public:
  AliTPC();
  AliTPC(const char *name, const char *title);
  virtual      ~AliTPC();
  virtual void  AddCluster(Float_t*, Int_t*);
  /*virtual*/void  AddCluster(const AliTPCcluster&);
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  AddTrack(Float_t*);
  /*virtual*/void  AddTrack(const AliTPCtrack&);
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials();
  virtual void  Hits2Clusters();
  virtual void  Hits2ExactClustersSector(Int_t isec); // MI change calculate "exact" cluster position

  virtual void  Hits2Digits();   //MI change
  virtual void Hits2DigitsSector(Int_t isec);  //MI change
  virtual void  Init();
  virtual Int_t IsVersion() const =0;
  virtual void  Digits2Clusters();
  virtual void  Clusters2Tracks();
  TClonesArray  *Clusters() {return fClusters;}
  TClonesArray  *Tracks()   {return fTracks;}

  Int_t         GetNsectors()       {return fNsectors;}
  virtual void  MakeBranch(Option_t *opt=" ");
  virtual void  ResetDigits();
  virtual void  SetSecAL(Int_t sec);
  virtual void  SetSecAU(Int_t sec);
  virtual void  SetSecLows(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6);
  virtual void  SetSecUps (Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6,
			   Int_t s7,Int_t s8,Int_t s9,Int_t s10, Int_t s11, Int_t s12);
  virtual void  SetSens(Int_t sens);

  //MK changes

  //MK changes

  virtual void  SetSide(Float_t side);
  virtual void  SetGasMixt(Int_t nc,Int_t c1,Int_t c2,Int_t c3,Float_t p1,
                           Float_t p2,Float_t p3); 

  virtual void  StepManager()=0;
  virtual void  DrawDetector() {}
  AliTPCDigitsArray*  GetDigitsArray() {return fDigitsArray;} //MI change
  AliTPCClustersArray* GetClustersArray(){return fClustersArray;} //MI change
  AliTPCParam *GetParam(){return fTPCParam;} // M.K, M.I changes
  void SetParam(AliTPCParam *param){fTPCParam=param;} // M.K, M.I changes
  void SetDigitsArray(AliTPCDigitsArray* param) {fDigitsArray=param;}  //MI change
  void SetClustersArray(AliTPCClustersArray *clusters) {fClustersArray = clusters;} //MI change
private:
  //
  void DigitizeRow(Int_t irow,Int_t isec,TObjArray **rowTriplet);
  Float_t GetSignal(TObjArray *p1, Int_t ntr, TMatrix *m1, TMatrix *m2,
                    Int_t *IndexRange);
  void GetList (Float_t label,Int_t np,TMatrix *m,Int_t *IndexRange,
                Float_t **pList);
  void MakeSector(Int_t isec,Int_t nrows,TTree *TH,Stat_t ntracks,TObjArray **row);
  void TransportElectron(Float_t *xyz, Int_t *index);
  Int_t fCurrentIndex[4];// index[0] indicates coordinate system, 
                         // index[1] sector number, 
                         // index[2] pad row number  
                         // index[3] pad row number for which signal is calculated
  
  ClassDef(AliTPC,2)  // Time Projection Chamber class
};

//_____________________________________________________________________________

class AliTPCcluster : public TObject {
public:
  Int_t     fTracks[3];//labels of overlapped tracks
  Int_t     fSector;   //sector number
  Int_t     fPadRow;   //PadRow number
  Float_t   fY ;       //Y of cluster
  Float_t   fZ ;       //Z of cluster
  Float_t   fQ ;       //Q of cluster (in ADC counts)
  Float_t   fdEdX;     //dE/dX inside this cluster
  Float_t   fSigmaY2;  //Sigma Y square of cluster
  Float_t   fSigmaZ2;  //Sigma Z square of cluster
  
public:
  AliTPCcluster() {
    fTracks[0]=fTracks[1]=fTracks[2]=0; 
    fSector=fPadRow=0;
    fY=fZ=fQ=fdEdX=fSigmaY2=fSigmaZ2=0.;
  }
  AliTPCcluster(Float_t *hits, Int_t*);
  virtual ~AliTPCcluster() {;}
  void Use() {fQ=-fQ;} //if fQ<0 cluster is already associated with a track
  Int_t IsUsed() const {return (fQ<0) ? 1 : 0;}
  void GetXYZ(Float_t *x, const AliTPCParam *) const; //Get global x,y,z
  Bool_t IsSortable() const;
  Int_t Compare(TObject *o) ;
  ClassDef(AliTPCcluster,1)  // Time Projection Chamber clusters
};


//_____________________________________________________________________________

class AliTPCdigit : public AliDigit {
public:
   Int_t     fSector;     //array of volumes
   Int_t     fPadRow;     //Row number
   Int_t     fPad ;       //Pad number
   Int_t     fTime;       //Time bucket
   Int_t     fSignal;     //Signal amplitude
 
public:
   AliTPCdigit() {}
   AliTPCdigit(Int_t *tracks, Int_t *digits);
   virtual ~AliTPCdigit() {}
 
   ClassDef(AliTPCdigit,1)  // Time Projection Chamber digits
};
 
 
//_____________________________________________________________________________
 
class AliTPChit : public AliHit {
public:
   Int_t     fSector;     //sector number
   Int_t     fPadRow;     //Pad Row number
   Float_t   fQ ;         //charge
 
public:
   AliTPChit() {}
   AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
   virtual ~AliTPChit() {}
 
   ClassDef(AliTPChit,1)  // Time Projection Chamber hits
};
 
 
//_____________________________________________________________________________
class AliTPCtrack : public TObject {
//-----------------------------------------------------------------
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
   Double_t fAlpha;          // rotation angle
   Double_t fX;              // X-coordinate of this track (reference plane)
   TVector x;                // vector of track parameters
   TMatrix C;                // covariance matrix of track parameters
   TObjArray fClusters;      // clusters belonging to this track
   Double_t fChi2;           // total chi2 value for this track
public:
   AliTPCtrack(): x(5), C(5,5), fClusters(200) {fAlpha=fX=fChi2=0.;}
   AliTPCtrack(Float_t *hits);
   AliTPCtrack(const AliTPCcluster *c, const TVector& xx, const TMatrix& CC,
               Double_t xr, Double_t alpha); 
   AliTPCtrack(const AliTPCtrack& t);
   Int_t Compare(TObject *o);
   Int_t PropagateTo(Double_t xr,
                     Double_t x0=28.94,Double_t rho=0.9e-3,Double_t pm=0.139);
   void PropagateToVertex(
                   Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
   void Update(const AliTPCcluster* c, Double_t chi2);
   Int_t Rotate(Double_t angle);

   Bool_t IsSortable() const {return kTRUE;}
   void UseClusters() const ;
   Double_t GetPredictedChi2(const AliTPCcluster*) const ;
   Int_t GetLabel(Int_t nrows) const ;
   void GetPxPyPz(Double_t&, Double_t&, Double_t&) const ;
   Double_t GetdEdX(Double_t low, Double_t up) const ;

   Double_t GetX() const {return fX;}
   Double_t GetY() const {return x(0);}
   Double_t GetZ() const {return x(1);}
   Double_t GetC() const {return x(2);}
   Double_t GetEta() const {return x(3);}
   Double_t GetTgl() const {return x(4);}
   Double_t GetPt() const {return 0.3*0.2/GetC()/100;}
   Double_t GetP() const {
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetSigmaY2() const {return C(0,0);}
   Double_t GetSigmaZ2() const {return C(1,1);}
   Double_t GetSigmaC2() const {return C(2,2);}
   Double_t GetSigmaTgl2() const {return C(4,4);}
   Double_t GetAlpha() const {return fAlpha;}
   Double_t GetChi2() const {return fChi2;}
   operator Int_t() const {return fClusters.GetEntriesFast();}
   const AliTPCcluster *GetCluster(Int_t i) const {
      return (const AliTPCcluster *)fClusters.UncheckedAt(i);
   }
 
   ClassDef(AliTPCtrack,1)  // Time Projection Chamber reconstructed tracks
};



//-----------------------------------------------------------------
// Classes for internal tracking use.
//-----------------------------------------------------------------
const unsigned MAX_CLUSTER_PER_ROW=3500;

class AliTPCRow {
//-----------------------------------------------------------------
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
   unsigned num_of_clusters;
   const AliTPCcluster *clusters[MAX_CLUSTER_PER_ROW];
public:
   AliTPCRow() {num_of_clusters=0;}
   void InsertCluster(const AliTPCcluster*);

   operator Int_t() const {return num_of_clusters;}
   const AliTPCcluster* operator[](Int_t i) const {return clusters[i];}
   Int_t Find(Double_t y) const; 
};

class AliTPCSector {
//-----------------------------------------------------------------
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
protected:
   unsigned num_of_rows; 
   AliTPCRow *row;
   static AliTPCParam *param; 
public:
   AliTPCSector() { row = 0; num_of_rows=0; }
   virtual ~AliTPCSector() { delete[] row; }
   static void SetParam(AliTPCParam *p) { param=p; }
   AliTPCRow& operator[](Int_t i) const { return *(row+i); }
   Int_t GetNRows() const { return num_of_rows; }
   virtual Double_t GetX(Int_t l) const = 0;
   virtual Double_t GetMaxY(Int_t l) const = 0;
   virtual Double_t GetAlpha() const = 0;
   virtual Double_t GetAlphaShift() const = 0;
   virtual Int_t GetRowNumber(Double_t x) const = 0;
   virtual Double_t GetPadPitchWidth() const = 0;
};

class AliTPCSSector : public AliTPCSector {
//-----------------------------------------------------------------
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
public:
   AliTPCSSector(){
     if (!param) {
        fprintf(stderr,"AliTPCSSector: parameters are not set !\n");
        return;
     }
     num_of_rows=param->GetNRowLow();
     row=new AliTPCRow[num_of_rows];
   }
   virtual ~AliTPCSSector() {}
   Double_t GetX(Int_t l) const { return param->GetPadRowRadiiLow(l); }
   Double_t GetMaxY(Int_t l) const { return GetX(l)*tan(0.5*GetAlpha()); }
   Double_t GetAlpha() const {return param->GetInnerAngle();}
   Double_t GetAlphaShift() const {return param->GetInnerAngleShift();}
   Double_t GetPadPitchWidth() const {
      return param->GetInnerPadPitchWidth();
   }
   Int_t GetRowNumber(Double_t x) const {
      Double_t r=param->GetPadRowRadiiLow(param->GetNRowLow()-1);
      if (x > r) return param->GetNRowLow();
      r=param->GetPadRowRadiiLow(0);
      if (x < r) return -1;
      return Int_t((x-r)/param->GetInnerPadPitchLength() + 0.5);
   }
};

class AliTPCLSector : public AliTPCSector {
//-----------------------------------------------------------------
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
public:
   AliTPCLSector(){
     if (!param) {
        fprintf(stderr,"AliTPCLSector: parameters are not set !\n");
        return;
     }
     num_of_rows=param->GetNRowUp();
     row=new AliTPCRow[num_of_rows];
   }
   virtual ~AliTPCLSector() {}
   Double_t GetX(Int_t l) const { return param->GetPadRowRadiiUp(l); }
   Double_t GetMaxY(Int_t l) const { return GetX(l)*tan(0.5*GetAlpha()); }
   Double_t GetAlpha() const {return param->GetOuterAngle();}
   Double_t GetAlphaShift() const {return param->GetOuterAngleShift();}
   Double_t GetPadPitchWidth() const {
      return param->GetOuterPadPitchWidth();
   }
   Int_t GetRowNumber(Double_t x) const {
      Double_t r=param->GetPadRowRadiiUp(param->GetNRowUp()-1);
      if (x > r) return param->GetNRowUp();
      r=param->GetPadRowRadiiUp(0);
      if (x < r) return -1;
      return Int_t((x-r)/param->GetOuterPadPitchLength() + 0.5);
   }
};

#endif

