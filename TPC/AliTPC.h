#ifndef TPC_H
#define TPC_H
////////////////////////////////////////////////
//  Manager class for TPC                     //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h" 
#include "AliTPCSecGeo.h" 
#include "AliTPCParam.h"
#include <TMatrix.h>
#include <TTree.h>
#include <TClonesArray.h>

#define MAXTPCTBK 500

class AliTPCcluster;
class AliTPCtrack;
class AliTPCParam;
//MI changes
class AliTPCD;

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
  Int_t          *fClustersIndex;   // Index for each sector in fClusters
  Int_t          *fDigitsIndex;     // Index for each sector in fDigits

  //MK changes

  Float_t        fSide;  // selects left(-1), right(+1), or both(0) sides of the TPC
  Int_t          fNoComp; // number of a drift gas components
  Int_t          fMixtComp[3]; // drift gas components
  Float_t        fMixtProp[3]; // mixture proportions

  //   

  TClonesArray   *fClusters;        // List of clusters for all sectors
  TClonesArray   *fTracks;          // List of reconstructed tracks
  //MI changes
  AliTPCD * fDigParam;              //detector parameters  
public:
  AliTPC();
  AliTPC(const char *name, const char *title);
  virtual      ~AliTPC();
  virtual void  AddCluster(Float_t*, Int_t*);
  /*virtual*/void  AddCluster(const AliTPCcluster&);
  virtual void  AddDigit(Int_t*, Int_t*);
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  AddTrack(Float_t*);
  /*virtual*/void  AddTrack(const AliTPCtrack&);
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials();
  virtual void  Hits2Clusters();

  virtual void  Hits2Digits();   //MI change
  virtual void Hits2DigitsSector(Int_t isec);  //MI change
  virtual void  Init();
  virtual Int_t IsVersion() const =0;
  virtual void  Digits2Clusters();
  virtual void  Clusters2Tracks();
  TClonesArray  *Clusters() {return fClusters;}
  TClonesArray  *Tracks()   {return fTracks;}
  Int_t         *GetClustersIndex() {return fClustersIndex;}
  Int_t         *GetDigitsIndex()   {return fDigitsIndex;}
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

  virtual void  SetSide(Float_t side);
  virtual void  SetGasMixt(Int_t nc,Int_t c1,Int_t c2,Int_t c3,Float_t p1,
                           Float_t p2,Float_t p3); 

  virtual void  StepManager()=0;
  virtual void  DrawDetector() {}
  AliTPCD*  GetDigParam() {return fDigParam;} //MI change8
  void SetDigParam(AliTPCD* param) {fDigParam=param;}  //MI must think about it
private:
  //
  void ElDiff(Float_t *xyz);
  void MakeTriplet(Int_t row,TObjArray **rowTriplet, 
                   TObjArray **prow);

  void ExB(Float_t *xyz);
  void DigitizeRow(Int_t irow,Int_t isec,TObjArray **rowTriplet);
  Float_t GetSignal(TObjArray *p1, Int_t ntr, Int_t np, TMatrix *m1, TMatrix *m2,
                    Int_t *IndexRange);
  void GetList (Float_t label,Int_t np,TMatrix *m,Int_t *IndexRange,
                Float_t **pList);
  void MakeSector(Int_t isec,Int_t nrows,TTree *TH,Stat_t ntracks,TObjArray **row);
  void GetCrossTalk (Int_t iFlag,TObjArray *p,Int_t ntracks,Int_t *npads,
                     TMatrix *m);

  
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
  Float_t   fSigmaY2;  //Sigma Y square of cluster
  Float_t   fSigmaZ2;  //Sigma Z square of cluster
  
public:
  AliTPCcluster() {
    fTracks[0]=fTracks[1]=fTracks[2]=0; 
    fSector=fPadRow=0;
    fY=fZ=fQ=fSigmaY2=fSigmaZ2=0.;
  }
  AliTPCcluster(Float_t *hits, Int_t*);
  virtual ~AliTPCcluster() {;}
  void Use() {fTracks[0]=-fTracks[0];}
  int IsUsed() const {return (fTracks[0]<0) ? 1 : 0;}
  void GetXYZ(Float_t *x, const AliTPCParam *) const; //Get global x,y,z
  Bool_t    IsSortable() const;
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
 
const unsigned MAX_CLUSTER_PER_ROW=1500;
const Double_t FIELD=0.2;

class AliTPCtrack : public TObject {
   Double_t fAlpha;          // rotation angle
   Double_t ref;             // track reference plane (X-coordinate)
   TVector x;                // vector of track parameters
   TMatrix C;                // covariance matrix of track parameters
   TObjArray clusters;       // pointers to clusters belonging to this track
   Double_t chi2;            // total chi2 value for this track
public:
   AliTPCtrack(Float_t *hits);
   AliTPCtrack(const AliTPCcluster& c, const TVector& xx, const TMatrix& CC,
               const AliTPCParam *); 
   AliTPCtrack(const AliTPCtrack& t);
   int PropagateTo(Double_t x,
                    Double_t x0=28.94,Double_t rho=0.9e-3,Double_t pm=0.139);
   void PropagateToVertex(
                    Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
   void Update(const AliTPCcluster* c, Double_t chi2);
   int Rotate(Double_t angle);

   void UseClusters() const ;
   Double_t GetPredictedChi2(const AliTPCcluster*) const ;
   Double_t GetX() const {return ref;}
   Double_t GetY() const {return x(0);}
   Double_t GetC() const {return x(2);}
   Double_t GetY(Double_t x) const;
   Double_t GetZ() const {return x(1);}
   Double_t GetTgl() const {return x(4);}
   Double_t GetPt() const {return 0.3*FIELD/x(2)/100;}
   int GetLab() const ;
   Double_t GetSigmaY2() const {return C(0,0);}
   Double_t GetSigmaZ2() const {return C(1,1);}
   Double_t GetSigmaC2() const {return C(2,2);}
   Double_t GetSigmaTgl2() const {return C(4,4);}
   Double_t GetAlpha() const {return fAlpha;}
   Double_t GetChi2() const {return chi2;}
   operator int() const {return clusters.GetEntriesFast();}
   AliTPCcluster& operator[](int i) {
      return *((AliTPCcluster*)clusters.UncheckedAt(i));
   } 
   void GetPxPyPz(Double_t&, Double_t&, Double_t&) const ;
   void GetXYZ(Double_t& X,Double_t& Y,Double_t& Z) const {X=ref;Y=x(0);Z=x(1);}
 
   ClassDef(AliTPCtrack,1)  // Time Projection Chamber reconstructed tracks
};

//_____Classes for internal tracking use ______________________________________

class TrackSeed : public AliTPCtrack {
public:
   TrackSeed(const AliTPCcluster& c, const TVector& x, const TMatrix& C,
   const AliTPCParam *p) : AliTPCtrack(c,x,C,p) {}
   Bool_t IsSortable() const {return kTRUE;}
   Int_t Compare(TObject *o) {
      AliTPCtrack *t=(AliTPCtrack*)o;
      Double_t c =GetSigmaY2();
      Double_t co=t->GetSigmaY2();
      if (c>co) return 1;
      else if (c<co) return -1;
      return 0;
   }
};

class AliTPCRow {
   unsigned num_of_clusters;
   const AliTPCcluster *clusters[MAX_CLUSTER_PER_ROW];
public:
   AliTPCRow() {num_of_clusters=0;}
   void InsertCluster(const AliTPCcluster*);

   operator int() const {return num_of_clusters;}
   const AliTPCcluster* operator[](int i) const {return clusters[i];}
   int Find(Double_t y) const; 
};


class AliTPCSector {
protected:
   unsigned num_of_rows; 
   AliTPCRow *row;
   const AliTPCParam *param; 
public:
   AliTPCSector() { 
      row = 0; num_of_rows=0;
      param=0;
   }
   virtual void SetUp(AliTPCParam *p) {param=p;}
   virtual ~AliTPCSector() { delete[] row; }
   AliTPCRow& operator[](int i) const { return *(row+i); }
   virtual Double_t GetX(int l) const = 0;
   virtual Double_t GetMaxY(int l) const = 0;
   virtual Double_t GetAlpha() const = 0;
};

class AliTPCSSector : public AliTPCSector {
public:
   AliTPCSSector(){}
   virtual ~AliTPCSSector() {}
   virtual void SetUp(AliTPCParam *p) {
      param=p;
      num_of_rows=p->GetNRowLow();
      row=new AliTPCRow[num_of_rows];
   }
   Double_t GetX(int l) const { return param->GetPadRowRadiiLow(l); }
   Double_t GetAlpha() const {return alpha_low;}
   Double_t GetMaxY(int l) const { return GetX(l)*tan(0.5*GetAlpha()); }
};

class AliTPCLSector : public AliTPCSector {
public:
   AliTPCLSector(){}
   virtual ~AliTPCLSector() {}
   virtual void SetUp(AliTPCParam *p) {
      param=p;
      num_of_rows=p->GetNRowUp();
      row=new AliTPCRow[num_of_rows];
   }
   Double_t GetX(int l) const { return param->GetPadRowRadiiUp(l); }
   Double_t GetAlpha() const {return alpha_up;}
   Double_t GetMaxY(int l) const { return GetX(l)*tan(0.5*GetAlpha()); }
};




#endif

