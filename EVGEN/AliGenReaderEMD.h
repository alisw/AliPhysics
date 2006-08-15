#ifndef ALIGENREADEREMD_H
#define ALIGENREADEREMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class to read events from external (TNtupla) file
// Events -> neutron removal by EM dissociation of Pb nuclei
// Data from RELDIS code (by I. Pshenichov)
//
#include "AliGenReader.h"


// --------------------------------------------------------------------------
class AliGenReaderEMD : public AliGenReader
{
 public:
    AliGenReaderEMD();
    
    AliGenReaderEMD(const AliGenReaderEMD &reader);
    virtual ~AliGenReaderEMD();
    AliGenReaderEMD & operator=(const AliGenReaderEMD & rhs);
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent(){;}
    // Setter
    void SetIPSide(Int_t iside) {fIPSide = iside;}
    void SetPcToTrack(Int_t ipart) {fPcToTrack = ipart;}
    void SetStartEvent(Int_t nev) {fStartEvent = nev;}
    
 protected:
    Int_t           fStartEvent;      	// points to the first event to read
    Int_t           fNcurrent;      	// points to the current event to read
    Int_t           fNparticle;     	// number of particles
    TTree           *fTreeNtuple;   	// pointer to the TTree
    //
    Int_t 	    fIPSide;		// ZDC arm relative to IP
    Int_t 	    fPcToTrack;		// flag for particles to be tracked
    //
    // --- Declaration of leaves types
    // **** neutrons
    Int_t           fNnLeft;		// No. of neutrons emitted on left side            
    Float_t         fEnLeft; 		// Forward energy Left side
    Float_t	    fPxnLeft[70];     	// momentum x component - left side
    Float_t	    fPynLeft[70];     	// momentum y component - left side     
    Float_t	    fPznLeft[70];     	// momentum z component - left side   
    //  
    Int_t           fNnRight;		// No. of neutrons emitted on right side            
    Float_t         fEnRight; 		// Forward energy Right side
    Float_t	    fPxnRight[70];     	// momentum x component - right side
    Float_t	    fPynRight[70];     	// momentum y component - right side     
    Float_t	    fPznRight[70];     	// momentum z component - right side     
    //
    // **** protons
    Int_t           fNpLeft;		// No. of protons emitted on left side            
    Float_t         fEtapLeft; 		// Forward energy Left side
    Float_t	    fPxpLeft[70];     	// momentum x component - left side
    Float_t	    fPypLeft[70];     	// momentum y component - left side     
    Float_t	    fPzpLeft[70];     	// momentum z component - left side   
    //  
    Int_t           fNpRight;		// No. of protons emitted on right side            
    Float_t         fEtapRight; 	// Forward energy Right side
    Float_t	    fPxpRight[70];     	// momentum x component - right side
    Float_t	    fPypRight[70];     	// momentum y component - right side     
    Float_t	    fPzpRight[70];     	// momentum z component - right side
    
 private:
    void Copy(TObject&) const;
    ClassDef(AliGenReaderEMD,1) // Class to read EMD data
};
#endif






