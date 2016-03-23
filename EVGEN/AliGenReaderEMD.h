#ifndef ALIGENREADEREMD_H
#define ALIGENREADEREMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class to read events from external (TNtupla) file
// Events -> EM dissociation of Pb nuclei
// Data from RELDIS code (by I. Pshenichov)
//
#include "AliGenReader.h"


// --------------------------------------------------------------------------
class AliGenReaderEMD : public AliGenReader
{
 public:
    enum TrackedPc{kAll=0, kOnlyNucleons=1, kNotNucleons=2};
    
    AliGenReaderEMD();
    
    AliGenReaderEMD(const AliGenReaderEMD &reader);
    virtual ~AliGenReaderEMD();
    AliGenReaderEMD & operator=(const AliGenReaderEMD & rhs);
    // Initialise 
    virtual void Init();
    // Reader
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent();
    // Setters
    void TrackNotNucleons()  {fPcToTrack = kNotNucleons;}
    void TrackOnlyNucleons() {fPcToTrack = kOnlyNucleons;}
    void TrackAllParticles() {fPcToTrack = kAll;}
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
    Int_t           fNnAside;		// No. of neutrons emitted on left side            
    Float_t         fEnAside; 		// Forward energy A side
    Int_t           fnPDGCode;		// PDG code
    Float_t	    fPxnAside[70];     	// momentum x component A side
    Float_t	    fPynAside[70];     	// momentum y component A side	  
    Float_t	    fPznAside[70];     	// momentum z component A side	
    //  
    Int_t           fNnCside;		// No. of neutrons emitted on right side            
    Float_t         fEnCside; 		// Forward energy C side
    Float_t	    fPxnCside[70];     	// momentum x component C side
    Float_t	    fPynCside[70];     	// momentum y component C side	  
    Float_t	    fPznCside[70];     	// momentum z component C side	  
    //
    // **** protons
    Int_t           fNpAside;		// No. of protons emitted on left side            
    Float_t         fEtapAside; 	// Forward energy A side
    Int_t           fpPDGCode;		// PDG code
    Float_t	    fPxpAside[50];     	// momentum x component A side
    Float_t	    fPypAside[50];     	// momentum y component A side	  
    Float_t	    fPzpAside[50];     	// momentum z component A side	
    //  
    Int_t           fNpCside;		// No. of protons emitted on right side            
    Float_t         fEtapCside; 	// Forward energy C side
    Float_t	    fPxpCside[50];     	// momentum x component C side
    Float_t	    fPypCside[50];     	// momentum y component C side	  
    Float_t	    fPzpCside[50];     	// momentum z component C side
    //
    // **** pi +
    Int_t           fNppAside;		// No. of pi+ emitted pi+ on A side            
    Float_t         fEtappAside; 	// Forward energy pi+ A side
    Int_t           fppPDGCode;	// PDG code
    Float_t	    fPxppAside[30];     // momentum x component pi+ A side
    Float_t	    fPyppAside[30];     // momentum y component pi+ A side	  
    Float_t	    fPzppAside[30];     // momentum z component pi+ A side	
    //  
    Int_t           fNppCside;		// No. of pi+ emitted on C side            
    Float_t         fEtappCside; 	// Forward energy pi+ C side
    Float_t	    fPxppCside[30];     // momentum x component pi+ C side
    Float_t	    fPyppCside[30];     // momentum y component pi+ C side	  
    Float_t	    fPzppCside[30];     // momentum z component pi+ C side
    //
    // **** pi -
    Int_t           fNpmAside;		// No. of pi- emitted on A side            
    Float_t         fEtapmAside; 	// Forward energy pi- A side
    Int_t           fpmPDGCode;	// PDG code
    Float_t	    fPxpmAside[30];     // momentum x component pi- A side
    Float_t	    fPypmAside[30];     // momentum y component pi- A side	  
    Float_t	    fPzpmAside[30];     // momentum z component pi- A side	
    //  
    Int_t           fNpmCside;		// No. of pi- emitted on C side            
    Float_t         fEtapmCside; 	// Forward energy pi- C side
    Float_t	    fPxpmCside[30];     // momentum x component pi- C side
    Float_t	    fPypmCside[30];     // momentum y component pi- C side	  
    Float_t	    fPzpmCside[30];     // momentum z component pi- C side
    //
    // **** pi0
    Int_t           fNp0Aside;		// No. of pi0 emitted on A side            
    Float_t         fEtap0Aside; 	// Forward energy pi0 A side
    Int_t           fp0PDGCode;	// PDG code
    Float_t	    fPxp0Aside[30];     // momentum x component pi0 A side
    Float_t	    fPyp0Aside[30];     // momentum y component pi0 A side	  
    Float_t	    fPzp0Aside[30];     // momentum z component pi0 A side	
    //  
    Int_t           fNp0Cside;		// No. of pi0 emitted on C side            
    Float_t         fEtap0Cside; 	// Forward energy pi0 C side
    Float_t	    fPxp0Cside[30];     // momentum x component pi0 C side
    Float_t	    fPyp0Cside[30];     // momentum y component pi0 C side	  
    Float_t	    fPzp0Cside[30];     // momentum z component pi0 C side
    //
    // **** eta
    Int_t           fNetaAside;		// No. of eta emitted on A side            
    Float_t         fEtaetaAside; 	// Forward energy eta A side
    Int_t           fetaPDGCode;	// PDG code
    Float_t	    fPxetaAside[15];    // momentum x component eta A side
    Float_t	    fPyetaAside[15];    // momentum y component eta A side	 
    Float_t	    fPzetaAside[15];    // momentum z component eta A side     
    //  
    Int_t           fNetaCside;		// No. of eta emitted on C side            
    Float_t         fEtaetaCside; 	// Forward energy eta C side
    Float_t	    fPxetaCside[15];    // momentum x component eta C side
    Float_t	    fPyetaCside[15];    // momentum y component eta C side	 
    Float_t	    fPzetaCside[15];    // momentum z component eta C side
    //
    // **** omega
    Int_t           fNomegaAside;	// No. of omega emitted on A side            
    Float_t         fEtaomegaAside; 	// Forward energy omega A side
    Int_t           fomegaPDGCode;	// PDG code
    Float_t	    fPxomegaAside[15];  // momentum x component omega A side
    Float_t	    fPyomegaAside[15];  // momentum y component omega A side   
    Float_t	    fPzomegaAside[15];  // momentum z component omega A side	 
    //  
    Int_t           fNomegaCside;	// No. of omega emitted on C side            
    Float_t         fEtaomegaCside; 	// Forward energy omega C side
    Float_t	    fPxomegaCside[15];  // momentum x component omega C side
    Float_t	    fPyomegaCside[15];  // momentum y component omega C side	 
    Float_t	    fPzomegaCside[15];  // momentum z component omega C side
   
 private:
    void Copy(TObject&) const;
    ClassDef(AliGenReaderEMD, 2) // Class to read EMD data
};
#endif






