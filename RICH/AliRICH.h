#ifndef RICH_H
#define RICH_H
////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliHit.h"
#include "AliRICHConst.h"
#include "AliDigit.h" 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TRotMatrix.h>

static const int NCH=7;
typedef enum {mip, cerenkov} Response_t;

class AliRICHcluster;
class AliRICHchamber;
class AliRICHRecCluster;
class AliRICHCerenkov;

//----------------------------------------------
//----------------------------------------------
//
// Chamber segmentation virtual base class
//
class AliRICHsegmentation :
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
    virtual void Init(AliRICHchamber*)                 =0;
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
    virtual void FitXY(AliRICHRecCluster* Cluster,TClonesArray* RICHdigits)    =0;
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
    ClassDef(AliRICHsegmentation,1)
	};
//----------------------------------------------
//
// Chamber response virtual base class
//
class AliRICHresponse :
public TObject {
 public:
    //
    // Configuration methods
    virtual void   SetRSIGM(Float_t p1)                          =0;
    virtual void   SetMUCHSP(Float_t p1)                         =0;
    virtual void   SetMUSIGM(Float_t p1, Float_t p2)             =0;
    virtual void   SetMAXADC(Float_t p1)                         =0;
    //
    // Get member data
    virtual Float_t Chslope()                                    =0;
    virtual Float_t ChwX()                                       =0;
    virtual Float_t ChwY()                                       =0;
    virtual Float_t Nsigma()                                     =0;
    virtual Float_t adc_satm()                                   =0;
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss=0)                       =0;
  //    virtual Float_t IntPH()                                      =0;
    virtual Int_t FeedBackPhotons(Float_t *source, Float_t qtot) =0;
    // Charge disintegration
    virtual Float_t IntXY(AliRICHsegmentation *)                 =0;
    //
    // Identification
    virtual char* YourName()                                     =0;
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1)                        =0;
    virtual void   SetKx2(Float_t p1)                            =0;
    virtual void   SetKx4(Float_t p1)                            =0;
    virtual void   SetSqrtKy3(Float_t p1)                        =0;
    virtual void   SetKy2(Float_t p1)                            =0;
    virtual void   SetKy4(Float_t p1)                            =0;
    virtual void   SetPitch(Float_t p1)                          =0;
    ClassDef(AliRICHresponse,1)
	};
	
//----------------------------------------------
class AliRICHchamber :
public TObject
{
    
 public:
    
//Rotation matrices for each chamber
    
    TRotMatrix *fChamberMatrix;
    Float_t fChamberTrans[3];
    
 public:
    AliRICHchamber();
    ~AliRICHchamber(){}
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
    
// Set inner radius of sensitive volume 
    void SetRInner(Float_t rmin) {frMin=rmin;}
// Set outer radius of sensitive volum  
    void SetROuter(Float_t rmax) {frMax=rmax;}  
    
// Return inner radius of sensitive volume 
    Float_t RInner()            {return frMin;}
// Return outer radius of sensitive volum  
    Float_t ROuter()            {return frMax;}
    
//Transformation from Global to local coordinates, chamber-dependant
    void LocaltoGlobal(Float_t pos[3],Float_t Localpos[3]);
    
//Setting chamber specific rotation matrices
    
    void SetChamberTransform(Float_t Trans1,Float_t Trans2,Float_t Trans3,TRotMatrix *Matrix)
	
	{
	    fChamberMatrix=Matrix;
	    fChamberTrans[0]=Trans1;
	    fChamberTrans[1]=Trans2;
	    fChamberTrans[2]=Trans3;
	}
    
// Configure response model
    void    ResponseModel(Response_t res, AliRICHresponse* thisResponse);
    
    //  
// Configure segmentation model
    void    SegmentationModel(Int_t i, AliRICHsegmentation* thisSegmentation) {
	(*fSegmentation)[i-1] = thisSegmentation;
    }
//  
//  Get reference to response model
    AliRICHresponse*     GetResponseModel(Response_t res);
//  
//  Get reference to segmentation model
    AliRICHsegmentation*  GetSegmentationModel(Int_t isec) {
	return (AliRICHsegmentation *) (*fSegmentation)[isec-1];
    }
    Int_t Nsec()              {return fnsec;}
    void  SetNsec(Int_t nsec) {fnsec=nsec;}
//
// Member function forwarding to the segmentation and response models
//
// Calculate pulse height from energy loss  
    Float_t IntPH(Float_t eloss) {return ((AliRICHresponse*) (*fResponse)[0])->IntPH(eloss);}
    Float_t IntPH()              {return ((AliRICHresponse*) (*fResponse)[1])->IntPH(); }
//  
// Ask segmentation if signal should be generated  
    Int_t   SigGenCond(Float_t x, Float_t y, Float_t z)
	{
	    if (fnsec==1) {
		return ((AliRICHsegmentation*) (*fSegmentation)[0])
		    ->SigGenCond(x, y, z) ;
	    } else {
		return (((AliRICHsegmentation*) (*fSegmentation)[0])
			->SigGenCond(x, y, z)) ||
		    (((AliRICHsegmentation*) (*fSegmentation)[1])
		     ->SigGenCond(x, y, z)) ;
	  }
	}
//
// Initialisation of segmentation for hit  
    void    SigGenInit(Float_t x, Float_t y, Float_t z)
	{
	  
	    if (fnsec==1) {
		((AliRICHsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	    } else {
		((AliRICHsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
		((AliRICHsegmentation*) (*fSegmentation)[1])->SigGenInit(x, y, z) ;
	    }
	}
    
// Configuration forwarding
//
    void   SetRSIGM(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetRSIGM(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetRSIGM(p);
	}
    void   SetMUCHSP(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetMUCHSP(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetMUCHSP(p);
	}
    void   SetMUSIGM(Float_t p1, Float_t p2)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetMUSIGM(p1,p2);
	    ((AliRICHresponse*) (*fResponse)[1])->SetMUSIGM(p1,p2);
	}
    void   SetMAXADC(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetMAXADC(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetMAXADC(p);
	}
    void   SetSqrtKx3(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetSqrtKx3(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetSqrtKx3(p);
	}
    void   SetKx2(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetKx2(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetKx2(p);
	}
    void   SetKx4(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetKx4(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetKx4(p);
	}
    void   SetSqrtKy3(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetSqrtKy3(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetSqrtKy3(p);
	}
    void   SetKy2(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetKy2(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetKy2(p);
	}
    void   SetKy4(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetKy4(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetKy4(p);
	}
    
    void   SetPitch(Float_t p)
	{
	    ((AliRICHresponse*) (*fResponse)[0])->SetPitch(p);
	    ((AliRICHresponse*) (*fResponse)[1])->SetPitch(p);
	}
    
    void   SetPADSIZ(Int_t isec, Float_t p1, Float_t p2) {
	((AliRICHsegmentation*) (*fSegmentation)[isec-1])->SetPADSIZ(p1,p2);
    }
//  
// Cluster formation method
    void   DisIntegration(Float_t, Float_t, Float_t, Int_t&x, Float_t newclust[6][500], Response_t res);
    ClassDef(AliRICHchamber,1)
	
	private:
    
// Maximum and Minimum Chamber size
    Float_t frMin;    
    Float_t frMax;
// GEANT volume if for sensitive volume of this chamber
    Int_t   fGid;
// z-position of this chamber
    Float_t fzPos;
// The segmentation models for the cathode planes
// fnsec=1: one plane segmented, fnsec=2: both planes are segmented.
    Int_t   fnsec;
    TObjArray           *fSegmentation;
    TObjArray           *fResponse;
    
};


 
class AliRICHcluster : public TObject {
public:
    
    Int_t     fHitNumber;    // Hit number
    Int_t     fCathode;      // Cathode number
    Int_t     fQ  ;          // Total charge      
    Int_t     fPadX  ;       // Pad number along X
    Int_t     fPadY  ;       // Pad number along Y
    Int_t     fQpad  ;       // Charge per pad
    Int_t     fRSec  ;       // R -sector of pad
    
public:
    AliRICHcluster() {
	fHitNumber=fQ=fPadX=fPadY=fQpad=fRSec=0;   
    }
    AliRICHcluster(Int_t *clhits);
    virtual ~AliRICHcluster() {;}
    
    ClassDef(AliRICHcluster,1)  //Cluster object for set:RICH
	};
	
	
	class AliRICHreccluster : public TObject {
public:
	    
	    Int_t     fTracks[3];      //labels of overlapped tracks
	    
	    Int_t       fQ  ;          // Q of cluster (in ADC counts)     
	    Float_t     fX  ;          // X of cluster
	    Float_t     fY  ;          // Y of cluster
	    
	public:
	    AliRICHreccluster() {
		fTracks[0]=fTracks[1]=fTracks[2]=0; 
		fQ=0; fX=fY=0;   
	    }
	    virtual ~AliRICHreccluster() {;}
	    
	    ClassDef(AliRICHreccluster,1)  //Cluster object for set:RICH
		};
		
//_____________________________________________________________________________

class AliRICHdigit : public TObject {
 public:
    Int_t     fPadX;        // Pad number along x
    Int_t     fPadY ;       // Pad number along y
    Int_t     fSignal;      // Signal amplitude
    
    
    Int_t     fTcharges[10];  // charge per track making this digit (up to 10)
    Int_t     fTracks[10];    // tracks making this digit (up to 10)
    
    
    
 public:
    AliRICHdigit() {}
    AliRICHdigit(Int_t *digits);
    AliRICHdigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliRICHdigit() {}
    
    
    ClassDef(AliRICHdigit,1)  //Digits for set:RICH
	};
//_____________________________________________________________________________

class AliRICHlist : public AliRICHdigit {
 public:
    
    Int_t          fRpad;       // r_pos of pad
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList; 

    
 public:
    AliRICHlist() {fTrackList=0;}
    AliRICHlist(Int_t ich, Int_t *digits);
    virtual ~AliRICHlist() {}

    TObjArray  *TrackList()   {return fTrackList;}
    
    ClassDef(AliRICHlist,1)  //Digits for set:RICH
	};
//___________________________________________


//___________________________________________

class AliRICHhit : public AliHit {
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
    AliRICHhit() {}
    AliRICHhit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliRICHhit() {}
    
    ClassDef(AliRICHhit,1)  //Hits object for set:RICH
	};
	
//------------------------------------------------
// Cerenkov photon  object
//------------------------------------------------

class AliRICHCerenkov: public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Int_t     fPHfirst;       // first padhit
    Int_t     fPHlast;        // last padhit
 public:
    AliRICHCerenkov() {}
    AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
    virtual ~AliRICHCerenkov() {}
    
    ClassDef(AliRICHCerenkov,1)  //Cerenkovs object for set:RICH
	};
	
//--------------------------------------------------

class AliRICH : public  AliDetector {
 public:
    AliRICH();
    AliRICH(const char *name, const char *title);
    virtual       ~AliRICH();
    virtual void   AddHit(Int_t, Int_t*, Float_t*);
    virtual void   AddCerenkov(Int_t, Int_t*, Float_t*);
    virtual void   AddCluster(Int_t*);
    virtual void   AddDigits(Int_t, Int_t*, Int_t*, Int_t*);
    virtual void   AddRecCluster(Int_t iCh, Int_t iCat,
				 AliRICHRecCluster* Cluster);
    virtual void   BuildGeometry();
    virtual void   CreateGeometry() {}
    virtual void   CreateMaterials() {}
    virtual void   StepManager();
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const =0;
//
    TClonesArray  *Clusters() {return fClusters;}
    TClonesArray  *Cerenkovs() {return fCerenkovs;}    
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetRecClusters();
    virtual void   ReconstructClusters();
    virtual void   Digitise(Int_t,Option_t *opt=" ",Text_t *name=" ");
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
    virtual void   SetSegmentationModel(Int_t id, Int_t isec, AliRICHsegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, Response_t res, AliRICHresponse *response);
    virtual void   SetNsec(Int_t id, Int_t nsec);
// Set Stepping Parameters
    virtual void   SetSMAXAR(Float_t p1);
    virtual void   SetSMAXAL(Float_t p1);
    virtual void   SetDMAXAR(Float_t p1);
    virtual void   SetDMAXAL(Float_t p1);
    virtual void   SetRICHACC(Bool_t acc=0, Float_t angmin=2, Float_t angmax=9);
// Response Simulation
    virtual void   MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id,Response_t res);
// Return reference to Chamber #id
    virtual AliRICHchamber& Chamber(Int_t id) {return *((AliRICHchamber *) (*fChambers)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliRICHcluster* FirstPad(AliRICHhit *, TClonesArray *);
    virtual AliRICHcluster* NextPad(TClonesArray *);
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
    Int_t                fNcerenkovs;          // Number of cerenkovs
    TClonesArray         *fClusters;           // List of clusters
    TObjArray            *fDchambers;          // List of digits
    TObjArray            *fRecClusters;        // List of clusters
    TClonesArray         *fCerenkovs;          // List of cerenkovs
    Int_t                *fNdch;               // Number of digits
    Text_t               *fFileName;           // Filename for event mixing
    

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
    
 protected:
    
    ClassDef(AliRICH,1)  //Hits manager for set:RICH
	};
//___________________________________________
class AliRICHRecCluster : public TObject {
 public:
    AliRICHRecCluster() ;
    AliRICHRecCluster(Int_t FirstDigit,Int_t Ichamber, Int_t Icathod) ;
    virtual ~AliRICHRecCluster();
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
    
    ClassDef(AliRICHRecCluster,1)  //Cluster object for set:RICH
	};
//___________________________________________
#endif















