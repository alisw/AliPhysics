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
    Int_t           fNnAside;		// No. of neutrons emitted on left side            
    Float_t         fEnAside; 		// Forward energy Aside side
    Float_t	    fPxnAside[70];     	// momentum x component - left side
    Float_t	    fPynAside[70];     	// momentum y component - left side     
    Float_t	    fPznAside[70];     	// momentum z component - left side   
    //  
    Int_t           fNnCside;		// No. of neutrons emitted on right side            
    Float_t         fEnCside; 		// Forward energy Cside side
    Float_t	    fPxnCside[70];     	// momentum x component - right side
    Float_t	    fPynCside[70];     	// momentum y component - right side     
    Float_t	    fPznCside[70];     	// momentum z component - right side     
    //
    // **** protons
    Int_t           fNpAside;		// No. of protons emitted on left side            
    Float_t         fEtapAside; 		// Forward energy Aside side
    Float_t	    fPxpAside[70];     	// momentum x component - left side
    Float_t	    fPypAside[70];     	// momentum y component - left side     
    Float_t	    fPzpAside[70];     	// momentum z component - left side   
    //  
    Int_t           fNpCside;		// No. of protons emitted on right side            
    Float_t         fEtapCside; 	// Forward energy Cside side
    Float_t	    fPxpCside[70];     	// momentum x component - right side
    Float_t	    fPypCside[70];     	// momentum y component - right side     
    Float_t	    fPzpCside[70];     	// momentum z component - right side
    
 private:
    void Copy(TObject&) const;
    ClassDef(AliGenReaderEMD,1) // Class to read EMD data
};
#endif






