#ifndef MUON_H
#define MUON_H
////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliHit.h"
#include "AliMUONConst.h"
#include "AliDigit.h"
#include "AliMUONchamber.h"
#include "AliMUONSegRes.h"
#include <TVector.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TFile.h>
#include <TTree.h>
typedef enum {simple, medium, big} Cluster_t;

static const int NCH=14;

class AliMUONcluster;
class AliMUONRecCluster;
class AliMUONRawCluster;
class AliMUONClusterFinder;
class AliMUONcorrelation;


//----------------------------------------------

 
class AliMUONcluster : public TObject {
public:

   Int_t     fHitNumber;    // Hit number
   Int_t     fCathode;      // Cathode number
   Int_t     fQ  ;          // Total charge      
   Int_t     fPadX  ;       // Pad number along X
   Int_t     fPadY  ;       // Pad number along Y
   Int_t     fQpad  ;       // Charge per pad
   Int_t     fRSec  ;       // R -sector of pad
 
public:
   AliMUONcluster() {
      fHitNumber=fQ=fPadX=fPadY=fQpad=fRSec=0;   
}
   AliMUONcluster(Int_t *clhits);
   virtual ~AliMUONcluster() {;}
 
   ClassDef(AliMUONcluster,1)  //Cluster object for set:MUON
};

 
class AliMUONreccluster : public TObject {
public:

   Int_t     fTracks[3];      //labels of overlapped tracks

   Int_t       fQ  ;          // Q of cluster (in ADC counts)     
   Float_t     fX  ;          // X of cluster
   Float_t     fY  ;          // Y of cluster

public:
   AliMUONreccluster() {
       fTracks[0]=fTracks[1]=fTracks[2]=0; 
       fQ=0; fX=fY=0; 
   }
   virtual ~AliMUONreccluster() {;}
 
   ClassDef(AliMUONreccluster,1)  //Cluster object for set:MUON
};

//_____________________________________________________________________________

class AliMUONdigit : public TObject {
 public:
    Int_t     fPadX;        // Pad number along x
    Int_t     fPadY ;       // Pad number along y
    Int_t     fSignal;      // Signal amplitude
    Int_t     fTcharges[10];   // charge per track making this digit (up to 10)
    Int_t     fTracks[10];     // primary tracks making this digit (up to 10)
    Int_t     fPhysics;        // physics contribution to signal 
    Int_t     fHit;            // hit number - temporary solution


 
 public:
    AliMUONdigit() {}
    AliMUONdigit(Int_t *digits);
    AliMUONdigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliMUONdigit();
    
    ClassDef(AliMUONdigit,1)  //Digits for set:MUON
};
//_____________________________________________________________________________

class AliMUONlist : public AliMUONdigit {
 public:
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList; 
 public:
    AliMUONlist() {fTrackList=0;}
    AliMUONlist(Int_t rpad, Int_t *digits);
    virtual ~AliMUONlist() {delete fTrackList;}
    TObjArray  *TrackList()   {return fTrackList;}
    ClassDef(AliMUONlist,1)  //Digits for set:MUON
};
//___________________________________________


//___________________________________________
 
class AliMUONhit : public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas   
    Int_t     fPHfirst;       // first padhit
    Int_t     fPHlast;        // last padhit

// modifs perso
    Float_t   fPTot;          // hit momentum P
    Float_t   fCxHit;            // Px/P
    Float_t   fCyHit;            // Py/P
    Float_t   fCzHit;            // Pz/P

 public:
    AliMUONhit() {}
    AliMUONhit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliMUONhit() {}
    
    ClassDef(AliMUONhit,1)  //Hits object for set:MUON
};

class AliMUON : public  AliDetector {
 public:
    AliMUON();
    AliMUON(const char *name, const char *title);
    virtual       ~AliMUON();
    virtual void   AddHit(Int_t, Int_t*, Float_t*);
    virtual void   AddCluster(Int_t*);
    virtual void   AddDigits(Int_t, Int_t*, Int_t*, Int_t*);
    virtual void   AddRawCluster(Int_t, const AliMUONRawCluster&);
    virtual void   AddRecCluster(Int_t iCh, Int_t iCat,
				 AliMUONRecCluster* Cluster);
    virtual void   AddCathCorrel(Int_t, Int_t*, Float_t*, Float_t*);
    virtual void   BuildGeometry();
    virtual void   CreateGeometry() {}
    virtual void   CreateMaterials() {}
    virtual void   StepManager();
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const =0;
//
    TClonesArray  *Clusters() {return fClusters;}
    virtual  void  MakeTreeC(Option_t *option="C");
    void           GetTreeC(Int_t);
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetRawClusters();
    virtual void   ResetRecClusters();
    virtual void   ResetCorrelation();
    virtual void   FindClusters(Int_t,Int_t);
    virtual void   Digitise(Int_t,Int_t,Option_t *opt1=" ",Option_t *opt2=" ",Text_t *name=" ");
    virtual void   CathodeCorrelation(Int_t);
    virtual void   SortTracks(Int_t *,Int_t *,Int_t);
//
// modifs perso

    void     Init(Double_t &, Double_t &, Double_t &);
    void     Reconst(Int_t &,Int_t &,Int_t,Int_t &,Int_t&,Int_t&, Option_t *option,Text_t *filename);
    void     FinishEvent();
    void     Close();
    void     SetCutPxz(Double_t p) {fSPxzCut=p;}
    void     SetSigmaCut(Double_t p) {fSSigmaCut=p;}
    void     SetXPrec(Double_t p) {fSXPrec=p;}
    void     SetYPrec(Double_t p) {fSYPrec=p;}
    Double_t GetCutPxz() {return fSPxzCut;}
    Double_t GetSigmaCut() {return fSSigmaCut;}
    Double_t GetXPrec() {return fSXPrec;}
    Double_t GetYPrec() {return fSYPrec;}
// fin modifs perso 
    
// Configuration Methods (per station id)
//
// Set Chamber Segmentation Parameters
// id refers to the station and isec to the cathode plane   
    virtual void   SetPADSIZ(Int_t id, Int_t isec, Float_t p1, Float_t p2);

// Set Signal Generation Parameters
    virtual void   SetSigmaIntegration(Int_t id, Float_t p1);
    virtual void   SetChargeSlope(Int_t id, Float_t p1);
    virtual void   SetChargeSpread(Int_t id, Float_t p1, Float_t p2);
    virtual void   SetMaxAdc(Int_t id, Float_t p1);
// Set Segmentation and Response Model
    virtual void   SetSegmentationModel(Int_t id, Int_t isec, AliMUONsegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, AliMUONresponse *response);
    virtual void   SetNsec(Int_t id, Int_t nsec);
// Set Reconstruction Model
    virtual void   SetReconstructionModel(Int_t id, AliMUONClusterFinder *reconstruction);
// Set Stepping Parameters
    virtual void   SetMaxStepGas(Float_t p1);
    virtual void   SetMaxStepAlu(Float_t p1);
    virtual void   SetMaxDestepGas(Float_t p1);
    virtual void   SetMaxDestepAlu(Float_t p1);
    virtual void   SetMuonAcc(Bool_t acc=0, Float_t angmin=2, Float_t angmax=9);
// Response Simulation
    virtual void   MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id);
// Return reference to Chamber #id
    virtual AliMUONchamber& Chamber(Int_t id) {return *((AliMUONchamber *) (*fChambers)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliMUONcluster* FirstPad(AliMUONhit *, TClonesArray *);
    virtual AliMUONcluster* NextPad(TClonesArray *);
// Return pointers to digits 
    TObjArray            *Dchambers() {return fDchambers;}
    Int_t                *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id) {return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to reconstructed clusters
  //    virtual TObjArray *RecClusters(Int_t iCh, Int_t iCat) 
  //	{return ( (TObjArray*) (*fRecClusters)[iCh+iCat*10]);}

    TObjArray            *RawClusters() {return fRawClusters;}
    Int_t                *Nrawch() {return fNrawch;}
    virtual TClonesArray *RawClustAddress(Int_t id) {return ((TClonesArray *) (*fRawClusters)[id]);}

// modifs perso
    AliMUONRawCluster *RawCluster(Int_t ichamber, Int_t icathod, Int_t icluster);
    
    
    // Return pointers to list of correlated clusters
    TObjArray            *CathCorrel() {return fCathCorrel;}
    Int_t                *Ncorch() {return fNcorch;}
    virtual TClonesArray *CathCorrelAddress(Int_t id)
	{return ((TClonesArray *) (*fCathCorrel)[id]);}

// modifs perso
//  virtual TClonesArray *CathCorrelAddress2(Int_t id)
//	{return ((TClonesArray *) (*fCathCorrel2)[id]);}    
    
// Return pointer to TreeC
    TTree      *TreeC() {return fTreeC;} 
 protected:
    TObjArray            *fChambers;           // List of Tracking Chambers
    Int_t                fNclusters;           // Number of clusters
    TClonesArray         *fClusters;           // List of clusters
    TObjArray            *fDchambers;          // List of digits
    Int_t                *fNdch;               // Number of digits

  //    TObjArray            *fRecClusters;        // List of clusters

    TObjArray            *fRawClusters;            // List of raw clusters
    Int_t                *fNrawch;                 // Number of raw clusters
    TObjArray            *fCathCorrel;             // List of correlated clusters
    Int_t                *fNcorch;                 // Number of correl clusters
    TTree                *fTreeC;                  // Cathode correl index tree

//
    Bool_t   fAccCut;          //Transport acceptance cut
    Float_t  fAccMin;          //Minimum acceptance cut used during transport
    Float_t  fAccMax;          //Minimum acceptance cut used during transport
//  

//  Stepping Parameters
   Float_t fMaxStepGas;      // Maximum step size inside the chamber gas
   Float_t fMaxStepAlu;      // Maximum step size inside the chamber aluminum
   Float_t fMaxDestepGas;    // Maximum relative energy loss in gas
   Float_t fMaxDestepAlu;    // Maximum relative energy loss in aluminum
//
// modifs perso
//  Parameters for reconstruction program
   Double_t fSPxzCut;        // Pxz cut  (GeV/c) to begin the track finding
   Double_t fSSigmaCut;      // Number of sig. delimiting the searching areas
   Double_t fSXPrec;         // Chamber precision in X (cm) 
   Double_t fSYPrec;         // Chamber precision in Y (cm)

   Text_t *fFileName;
   
 protected:

   ClassDef(AliMUON,1)  //Hits manager for set:MUON
};
//___________________________________________
class AliMUONRecCluster : public TObject {
public:
   AliMUONRecCluster() ;
   AliMUONRecCluster(Int_t FirstDigit,Int_t Ichamber, Int_t Icathod) ;
   virtual ~AliMUONRecCluster();
   virtual void  AddDigit(Int_t Digit);
   virtual Int_t FirstDigitIndex();
   virtual Int_t NextDigitIndex();
   virtual Int_t InvalidDigitIndex() {return -1;}

   virtual Int_t NDigits();
   virtual void  Finish();    // Nothing yet ...
   virtual Int_t GetCathod()  {return fCathod;}
   virtual Int_t GetChamber() {return fChamber;}

public:
   Float_t fX; // reconstructed x
   Float_t fY; // reconstructed y

protected:
   TArrayI *fDigits;    // List of digits indexes for that cluster
   Int_t fNdigit;       // Number of digits indexes stored;
   Int_t fCathod;       // Number of the cathod to be used;
   Int_t fChamber;      // Number of the chamber to be used;
   Int_t fCurrentDigit; // Current Digit inside an iteration

   ClassDef(AliMUONRecCluster,1)  //Cluster object for set:MUON
};
//___________________________________________

class AliMUONRawCluster : public TObject {
public:

   Int_t     fTracks[3];      //labels of overlapped tracks
   Int_t       fQ  ;          // Q of cluster (in ADC counts)     
   Float_t     fX  ;          // X of cluster
   Float_t     fY  ;          // Y of cluster
   Int_t       fPeakSignal;
   Int_t       fIndexMap[50];  //indeces of digits
   Int_t       fOffsetMap[50]; //Emmanuel special
   Float_t     fContMap[50];   //Contribution from digit
   Int_t       fPhysicsMap[50];
   Int_t       fMultiplicity;  //cluster multiplicity
   Int_t       fNcluster[2];
   Int_t       fClusterType;   
 public:
   AliMUONRawCluster() {
       fTracks[0]=fTracks[1]=fTracks[2]=-1; 
       fQ=0; fX=fY=0; fMultiplicity=0;
       for (int k=0;k<50;k++) {
           fIndexMap[k]=-1;
           fOffsetMap[k]=0;
	   fContMap[k]=0;
	   fPhysicsMap[k]=-1;
       }
       fNcluster[0]=fNcluster[1]=-1;
   }
   virtual ~AliMUONRawCluster() {}

   Float_t GetRadius() {return TMath::Sqrt(fX*fX+fY*fY);}

   Bool_t IsSortable() const {return kTRUE;}
   Int_t  Compare(TObject *obj);
   Int_t PhysicsContribution();
   static Int_t BinarySearch(Float_t r, TArrayF, Int_t from, Int_t upto);
   static void  SortMin(Int_t *,Float_t *,Float_t *,Float_t *,Float_t *,Int_t);
 
   ClassDef(AliMUONRawCluster,1)  //Cluster object for set:MUON
};

//___________________________________________
class AliMUONcorrelation : public TObject {
public:

  // correlation starts from the 1-st cathode  
  // last number in arrays corresponds to cluster on 1-st cathode

   Int_t       fCorrelIndex[4];  // entry number in TreeR for the associated 
                                 // cluster candidates on the 2-nd cathode
   Float_t     fX[4]  ;          // X of clusters on the 2-nd cathode  
   Float_t     fY[4]  ;          // Y of clusters

public:
   AliMUONcorrelation() {
       fCorrelIndex[0]=fCorrelIndex[1]=fCorrelIndex[2]=fCorrelIndex[3]=0;
       fX[0]=fX[1]=fX[2]=fX[3]=0; fY[0]=fY[1]=fY[2]=fY[3]=0; 
   }
   AliMUONcorrelation(Int_t *idx, Float_t *x, Float_t *y);
   virtual ~AliMUONcorrelation() {}
   ClassDef(AliMUONcorrelation,1)  //Cathode correlation object for set:MUON
};

#endif















