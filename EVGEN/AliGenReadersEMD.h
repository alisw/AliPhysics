#ifndef ALIGENREADERSEMD_H
#define ALIGENREADERSEMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class to read events from external (TNtupla) file
// Events -> single EM dissociation of Pb nuclei
// Data from RELDIS code (by I. Pshenichov)
//
#include "AliGenReader.h"


// --------------------------------------------------------------------------
class AliGenReadersEMD : public AliGenReader
{
 public:
    enum TrackedPc{kAll=0, kNucleons=1, kOnlyNeutrons=2};
    
    AliGenReadersEMD();
    
    AliGenReadersEMD(const AliGenReadersEMD &reader);
    virtual ~AliGenReadersEMD();
    AliGenReadersEMD & operator=(const AliGenReadersEMD & rhs);
    // Initialise 
    virtual void Init();
    // Reader
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent();
    // Setters
    void TrackNucleons()  {fPcToTrack = kNucleons;}
    void TrackOnlyNeutrons() {fPcToTrack = kOnlyNeutrons;}
    void TrackAll() {fPcToTrack = kAll;}
    void SetStartEvent(Int_t nev) {fStartEvent = nev;}
    
 protected:
    Int_t           fStartEvent;      	// points to the first event to read
    Int_t           fNcurrent;      	// points to the current event to read
    Int_t           fNparticle;     	// number of particles
    TTree           *fTreeNtuple;   	// pointer to the TTree
    //
    Int_t 	    fPcToTrack;		// flag for particles to be tracked
    Int_t           fOffset;		// Needed to correctly read next particle
    //
    // --- Declaration of leaves types
    // **** neutrons
    Int_t           fNneu;		// No. of neutrons emitte on left side            
    Float_t         fEneu; 		// Energy
    Float_t	    fPxneu[70];	   	// momentum x component neutrons
    Float_t	    fPyneu[70];	   	// momentum y component neutrons
    Float_t	    fPzneu[70];	   	// momentum z component neutrons
    //  
    Float_t	    fPxfrag;	    	// momentum x component fragments
    Float_t	    fPyfrag;	    	// momentum y component fragments
    Float_t	    fPzfrag;	    	// momentum z component fragments
    Float_t	    fAfrag;	    	// A fragments
    Float_t	    fZfrag;	    	// Z fragments
    //
    // **** protons
    Int_t           fNpro;		// No. of protons emitted on left side            
    Float_t         fEpro; 		// Forward energy A side
    Float_t	    fPxpro[50];	   	// momentum x component A side
    Float_t	    fPypro[50];	   	// momentum y component A side    
    Float_t	    fPzpro[50];	   	// momentum z component A side  
   
 private:
    void Copy(TObject&) const;
    ClassDef(AliGenReadersEMD, 1) // Class to read EMD data
};
#endif






