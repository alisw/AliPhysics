#ifndef MUON_H
#define MUON_H
////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliHit.h"
#include "AliMUONConst.h"
#include "AliDigit.h" 
#include <TVector.h>
#include <TObjArray.h>


static const int NCH=14;

class AliMUONcluster;
class AliMUONchamber;
class AliMUONRecCluster;
//----------------------------------------------
class AliMUONgeometry 
{
 public:
    AliMUONgeometry(){}
    virtual ~AliMUONgeometry(){}
    void    InitGeo(Float_t z);
    Float_t fdGas; // half gaz gap
    Float_t fdAlu; // half Alu width
    Float_t frMin; // innermost sensitive radius
    Float_t frMax; // outermost sensitive radius
    ClassDef(AliMUONgeometry,1)
};
//----------------------------------------------
//
// Chamber segmentation virtual base class
//
class AliMUONsegmentation :
public TObject {
    
 public:

    // Set Chamber Segmentation Parameters
    virtual void    SetPADSIZ(Float_t p1, Float_t p2)  =0;
    virtual void    SetDAnod(Float_t D)                =0;
    // Transform from pad (wire) to real coordinates and vice versa
    virtual Float_t GetAnod(Float_t xhit)              =0;
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy)=0;
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y )=0;
    //
    // Initialisation
    virtual void Init(AliMUONchamber*)                 =0;
    //
    // Get member data
    virtual Float_t Dpx()                              =0;
    virtual Float_t Dpy()                              =0;
    virtual Int_t Npx()                                =0;
    virtual Int_t Npy()                                =0;
    //
    // Iterate over pads
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy) =0;
    virtual void  NextPad()=0;
    virtual Int_t MorePads()                           =0;
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])     =0;
    // Provisory RecCluster coordinates reconstructor
    virtual void FitXY(AliMUONRecCluster* Cluster,TClonesArray* MUONdigits)    =0;
    //
    // Current pad cursor during disintegration 
    virtual Int_t  Ix()                                =0;
    virtual Int_t  Iy()                                =0;
    virtual Int_t  ISector()                           =0;
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z) = 0;
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z) = 0;
    virtual void  IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2)  = 0;
    //
    // Identification
    virtual char* YourName()                           =0;
    ClassDef(AliMUONsegmentation,1)
};
//----------------------------------------------
//
// Chamber response virtual base class
//
class AliMUONresponse :
public TObject {
 public:
    //
    // Configuration methods
    virtual void   SetRSIGM(Float_t p1)                =0;
    virtual void   SetMUCHSP(Float_t p1)               =0;
    virtual void   SetMUSIGM(Float_t p1, Float_t p2)   =0;
    virtual void   SetMAXADC(Float_t p1)               =0;
    //
    // Get member data
    virtual Float_t Chslope()                          =0;
    virtual Float_t ChwX()                             =0;
    virtual Float_t ChwY()                             =0;
    virtual Float_t Nsigma()                           =0;
    virtual Float_t adc_satm()                         =0;
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss)               =0;
    // Charge disintegration
    virtual Float_t IntXY(AliMUONsegmentation *)       =0;
    //
    // Identification
    virtual char* YourName()                           =0;
    ClassDef(AliMUONresponse,1)
};

//----------------------------------------------
class AliMUONchamber :
public TObject,
public AliMUONgeometry{
 public:
    AliMUONchamber();
    ~AliMUONchamber(){}
//
// Set and get GEANT id  
  Int_t   GetGid()         {return fGid;}
  void    SetGid(Int_t id) {fGid=id;}
//  
// Initialisation and z-Position
  void    Init();
  void    SetZPOS(Float_t p1) {fzPos=p1;}
  Float_t ZPosition()         {return fzPos;}
//  
// Configure response model
  void    ResponseModel(AliMUONresponse* thisResponse) {fResponse=thisResponse;}
//  
// Configure segmentation model
  void    SegmentationModel(Int_t i, AliMUONsegmentation* thisSegmentation) {
      (*fSegmentation)[i-1] = thisSegmentation;
  }
//  
//  Get reference to response model
  AliMUONresponse*     &GetResponseModel(){return fResponse;}
//  
//  Get reference to segmentation model
  AliMUONsegmentation*  GetSegmentationModel(Int_t isec) {
      return (AliMUONsegmentation *) (*fSegmentation)[isec-1];
  }
  Int_t Nsec()              {return fnsec;}
  void  SetNsec(Int_t nsec) {fnsec=nsec;}
//
// Member function forwarding to the segmentation and response models
//
// Calculate pulse height from energy loss  
  Float_t IntPH(Float_t eloss) {return fResponse->IntPH(eloss);}
//  
// Ask segmentation if signal should be generated  
  Int_t   SigGenCond(Float_t x, Float_t y, Float_t z)
      {
	  if (fnsec==1) {
	      return ((AliMUONsegmentation*) (*fSegmentation)[0])
		  ->SigGenCond(x, y, z) ;
	  } else {
	      return (((AliMUONsegmentation*) (*fSegmentation)[0])
		      ->SigGenCond(x, y, z)) ||
		  (((AliMUONsegmentation*) (*fSegmentation)[1])
		   ->SigGenCond(x, y, z)) ;
	  }
  }
//
// Initialisation of segmentation for hit  
  void    SigGenInit(Float_t x, Float_t y, Float_t z)
      {
	  
	  if (fnsec==1) {
	      ((AliMUONsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	  } else {
	      ((AliMUONsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	      ((AliMUONsegmentation*) (*fSegmentation)[1])->SigGenInit(x, y, z) ;
	  }
      }
  
// Configuration forwarding
//
  void   SetRSIGM(Float_t p1)              {fResponse->SetRSIGM(p1);}
  void   SetMUCHSP(Float_t p1)             {fResponse->SetMUCHSP(p1);}
  void   SetMUSIGM(Float_t p1, Float_t p2) {fResponse->SetMUSIGM(p1,p2);}
  void   SetMAXADC(Float_t p1)             {fResponse->SetMAXADC(p1);}

  void   SetPADSIZ(Int_t isec, Float_t p1, Float_t p2) {
      ((AliMUONsegmentation*) (*fSegmentation)[isec-1])->SetPADSIZ(p1,p2);
  }
//  
// Cluster formation method
  void   DisIntegration(Float_t, Float_t, Float_t, Int_t&x, Float_t newclust[6][500]);
    ClassDef(AliMUONchamber,1)

 private:
// GEANT volume if for sensitive volume of this chamber
  Int_t   fGid;
// z-position of this chamber
  Float_t fzPos;
// The segmentation models for the cathode planes
// fnsec=1: one plane segmented, fnsec=2: both planes are segmented.
  Int_t   fnsec;
  TObjArray           *fSegmentation;
  AliMUONresponse     *fResponse;

};


 
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
    

    Int_t     fTcharges[10];  // charge per track making this digit (up to 10)
    Int_t     fTracks[10];    // tracks making this digit (up to 10)


 
 public:
    AliMUONdigit() {}
    AliMUONdigit(Int_t *digits);
    AliMUONdigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliMUONdigit() {}

 
    ClassDef(AliMUONdigit,1)  //Digits for set:MUON
};
//_____________________________________________________________________________

class AliMUONlist : public AliMUONdigit {
 public:
    
    Int_t          fRpad;       // r_pos of pad
    
    TObjArray     *fTrackList; 

 
 public:
    AliMUONlist() {fTrackList=0;}
    AliMUONlist(Int_t rpad, Int_t *digits);
    virtual ~AliMUONlist() {}

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
    virtual void   AddRecCluster(Int_t iCh, Int_t iCat,
				 AliMUONRecCluster* Cluster);
    virtual void   BuildGeometry();
    virtual void   CreateGeometry() {}
    virtual void   CreateMaterials() {}
    virtual void   StepManager();
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const =0;
//
    TClonesArray  *Clusters() {return fClusters;}
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetRecClusters();
    virtual void   ReconstructClusters();
// 
// Configuration Methods (per station id)
//
// Set Chamber Segmentation Parameters
// id refers to the station and isec to the cathode plane   
    virtual void   SetPADSIZ(Int_t id, Int_t isec, Float_t p1, Float_t p2);

// Set Signal Generation Parameters
    virtual void   SetRSIGM(Int_t id, Float_t p1);
    virtual void   SetMUCHSP(Int_t id, Float_t p1);
    virtual void   SetMUSIGM(Int_t id, Float_t p1, Float_t p2);
    virtual void   SetMAXADC(Int_t id, Float_t p1);
// Set Segmentation and Response Model
    virtual void   SetSegmentationModel(Int_t id, Int_t isec, AliMUONsegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, AliMUONresponse *response);
    virtual void   SetNsec(Int_t id, Int_t nsec);
// Set Stepping Parameters
    virtual void   SetSMAXAR(Float_t p1);
    virtual void   SetSMAXAL(Float_t p1);
    virtual void   SetDMAXAR(Float_t p1);
    virtual void   SetDMAXAL(Float_t p1);
    virtual void   SetMUONACC(Bool_t acc=0, Float_t angmin=2, Float_t angmax=9);
// Response Simulation
    virtual void   MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id);
// Return reference to Chamber #id
    virtual AliMUONchamber& Chamber(Int_t id) {return *((AliMUONchamber *) (*fChambers)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliMUONcluster* FirstPad(AliMUONhit *);
    virtual AliMUONcluster* NextPad();
// Return pointers to digits 
    TObjArray            *Dchambers() {return fDchambers;}
    Int_t                *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id) {return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to reconstructed clusters
    virtual TObjArray *RecClusters(Int_t iCh, Int_t iCat) 
	{return ( (TObjArray*) (*fRecClusters)[iCh+iCat*10]);}

   
 protected:
    TObjArray            *fChambers;           // List of Tracking Chambers
    Int_t                fNclusters;           // Number of clusters
    TClonesArray         *fClusters;           // List of clusters
    TObjArray            *fDchambers;          // List of digits
    TObjArray            *fRecClusters;        // List of clusters
    Int_t                *fNdch;               // Number of digits
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
#endif















