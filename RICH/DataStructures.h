#ifndef RICHDataStructures_H
#define RICHDataStructures_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliHit.h"
#include "AliDigit.h"
#include <TObjArray.h>
#include <TArrayF.h>
#include <TMath.h>


class AliRICHHit : public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas   
    Float_t   fPHfirst;       // first padhit
    Float_t   fPHlast;        // last padhit
    Float_t   fLoss;          // did it hit the freon?
    Float_t   fMomX;            // Local Momentum
    Float_t   fMomY;            // Local Momentum
    Float_t   fMomZ;            // Local Momentum
 public:
    AliRICHHit() {}
    AliRICHHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliRICHHit() {}
    
    ClassDef(AliRICHHit,1)  //Hits object for set:RICH
};
	
//------------------------------------------------
// Cerenkov photon  object
//------------------------------------------------

class AliRICHCerenkov: public AliHit {
 public:
    Int_t     fChamber;         // Chamber number
    Float_t   fTheta ;          // Incident theta angle in degrees      
    Float_t   fPhi   ;          // Incident phi angle in degrees
    Float_t   fTlength;         // Track length inside the chamber
    Float_t   fEloss;           // ionisation energy loss in gas
    Int_t     fPHfirst;         // first padhit
    Int_t     fPHlast;          // last padhit
    Int_t     fCMother;         // index of mother particle
    Float_t   fLoss;            // nature of particle loss
    Float_t   fIndex;           // Index of photon
    Float_t   fProduction;      // Point of production
    Float_t   fMomX;            // Local Momentum
    Float_t   fMomY;            // Local Momentum
    Float_t   fMomZ;            // Local Momentum
 public:
    AliRICHCerenkov() {}
    AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
    virtual ~AliRICHCerenkov() {}
    
    ClassDef(AliRICHCerenkov,1)  //Cerenkovs object for set:RICH
};

// 
class AliRICHPadHit : public TObject {
 public:
    
    Int_t     fHitNumber;    // Hit number
    Int_t     fCathode;      // Cathode number
    Int_t     fQ  ;          // Total charge      
    Int_t     fPadX  ;       // Pad number along X
    Int_t     fPadY  ;       // Pad number along Y
    Int_t     fQpad  ;       // Charge per pad
    Int_t     fRSec  ;       // R -sector of pad
    
 public:
    AliRICHPadHit() {
	fHitNumber=fQ=fPadX=fPadY=fQpad=fRSec=0;   
    }
    AliRICHPadHit(Int_t *clhits);
    virtual ~AliRICHPadHit() {;}
    
    ClassDef(AliRICHPadHit,1)  //Cluster object for set:RICH
};

class AliRICHDigit : public TObject {
 public:
    Int_t     fPadX;        // Pad number along x
    Int_t     fPadY ;       // Pad number along y
    Int_t     fSignal;      // Signal amplitude
    
    
    Int_t     fTcharges[100];  // charge per track making this digit (up to 10)
    Int_t     fTracks[100];    // tracks making this digit (up to 10)
    Int_t     fPhysics;        // physics contribution to signal 
    Int_t     fHit;            // hit number - temporary solution
    
 public:
    AliRICHDigit() {}
    AliRICHDigit(Int_t *digits);
    AliRICHDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliRICHDigit() {}
    ClassDef(AliRICHDigit,1)  //Digits for set:RICH
};
//_____________________________________________________________________________

class AliRICHTransientDigit : public AliRICHDigit {
 public:
    Int_t          fRpad;       // r_pos of pad
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList; 
 public:
    AliRICHTransientDigit() {fTrackList=0;}
    AliRICHTransientDigit(Int_t ich, Int_t *digits);
    virtual ~AliRICHTransientDigit() {}
    
    TObjArray  *TrackList()   {return fTrackList;}
    
    ClassDef(AliRICHTransientDigit,1)  //Digits for set:RICH
};
//___________________________________________
class AliRICHRawCluster : public TObject {
public:
    Int_t       fTracks[3];      //labels of overlapped tracks
    Int_t       fQ  ;            // Q of cluster (in ADC counts)     
    Float_t     fX  ;            // X of cluster
    Float_t     fY  ;            // Y of cluster
    Int_t       fPeakSignal;
    Int_t       fIndexMap[50];   //indeces of digits
    Int_t       fOffsetMap[50];  
    Float_t     fContMap[50];    //Contribution from digit
    Int_t       fPhysicsMap[50];
    Int_t       fMultiplicity;   //cluster multiplicity
    Int_t       fNcluster[2];
    Int_t       fClusterType;
    Int_t       fCtype;          //CL0, CL1, etc...
 public:
    AliRICHRawCluster() {
	fTracks[0]=fTracks[1]=fTracks[2]=-1; 
	fQ=0; fX=fY=0; fMultiplicity=0;
	for (int k=0;k<50;k++) {
	    fIndexMap[k]=-1;
	    fOffsetMap[k]=0;
	    fContMap[k]=0;
	    fPhysicsMap[k]=-1;
	    fCtype=-1;
	}
	fNcluster[0]=fNcluster[1]=-1;
    }
    virtual ~AliRICHRawCluster() {}
    
    Float_t GetRadius() {return TMath::Sqrt(fX*fX+fY*fY);}
    
    Bool_t IsSortable() const {return kTRUE;}
    Int_t  Compare(TObject *obj);
    Int_t PhysicsContribution();
    static Int_t BinarySearch(Float_t r, TArrayF, Int_t from, Int_t upto);
    static void  SortMin(Int_t *,Float_t *,Float_t *,Float_t *,Float_t *,Int_t);
    
   ClassDef(AliRICHRawCluster,1)  //Cluster object for set:RICH
};
//___________________________________________
class AliRICHRecHit : public TObject {
public:
  Float_t     Theta  ;            //Incidence Angle theta
  Float_t     Phi  ;              //Incidence Angle phi
  Float_t     Omega;              //Cherenkov angle omega
  Float_t     fX;                 //Impact coordinate x
  Float_t     fY;                 //Impact coordinate y
 public:
    AliRICHRecHit() {
      Theta=Phi=Omega=0;
    }
    AliRICHRecHit(Int_t id, Float_t* rechit);
    virtual ~AliRICHRecHit() {}
    ClassDef(AliRICHRecHit,1)  //Reconstructed hit object for set:RICH
};

#endif
