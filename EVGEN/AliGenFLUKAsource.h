#ifndef ALIGENFLUKASOURCE_H
#define ALIGENFLUKASOURCE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Read background particles from a FLUKA boundary source file
// This is a very special generator that works for background studies for the muon-spectrometer 
// Ask: andreas.morsch@cern.ch

#include "AliGenerator.h"
class TChain;
class TTree;
class AliGenFLUKAsource : public AliGenerator

{
public:
    enum constants {kAll = 6, kGammas = 7, kNeutrons = 8, kCharged = 9, kNoNeutron = 10};

    AliGenFLUKAsource();
    AliGenFLUKAsource(Int_t npart);
     virtual ~AliGenFLUKAsource();
    // Initialise 
    virtual void Init() {}
    // Initialise fluka data 
    virtual void FlukaInit();
    // choose particle type
    virtual void SetPartFlag(Int_t ikine) {fIkine=ikine;}
    // set time cut 
    virtual void SetAgeMax(Float_t agemax) {fAgeMax=agemax;}
    // use additional weight on neutrals
    virtual void SetAddWeight(Float_t addwgt) {fAddWeight=addwgt;}
    // z-shift of vertex
    virtual void SetZshift(Float_t zshift) {fZshift=zshift;}
    // set file name of data file
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    // set source
    virtual void SetSourceId(Int_t id=-1){fSourceId=id;}
    // add a new source file	  
    virtual void AddFile(const Text_t *filname) ;  
    // read only fraction of data  
    virtual void SetFraction(Float_t frac=1.){fFrac=frac;}
    // generate event
    virtual void Generate();

 protected:

    Int_t       fIkine;         // Flag to choose type of particles to be read
    Float_t     fAgeMax;        // Maximum age of particle
    Float_t     fAddWeight;     // Add weight for neutrons 
    Float_t     fZshift;        // Shift the Z of impact point by this quantity
    Float_t     fFrac;          // Fraction of file that corresponds to one event
    Int_t       fSourceId;      // Source identifier (-1: all sources)
  
  
    const Text_t    *fFileName;          //!Choose the file
    TChain          *fTreeChain;         //file chaining
    TTree           *fTreeFluka;         //pointer to the TTree
//Declaration of variables read from the file -- TTree type
    Float_t         fIp;     // Particle type
    Float_t         fIpp;    // Primary particle type
    Float_t         fXi;     // x-Impact 
    Float_t         fYi;     // y-Impact
    Float_t         fZi;     // z-Impact
    Float_t         fPx;     // Direction cosine x
    Float_t         fPy;     // Direction cosine y
    Float_t         fPz;     // Direction cosine z
    Float_t         fEkin;   // Kinetic energy
    Float_t         fZv;     // z-Position of particle vertex
    Float_t         fRv;     // r-Position of particle vertex
    Float_t         fItra;   // Primary track number
    Float_t         fIgas;   // Volume identifier
    Float_t         fWgt;    // Particle weight
    Float_t         fEtag;   // Pseudorapidity of primary particle
    Float_t         fPtg;    // Pt of primary particle
    Float_t         fAge;    // Time of flight

 private:
    AliGenFLUKAsource(const AliGenFLUKAsource &FLUKAsource);
    AliGenFLUKAsource & operator=(const AliGenFLUKAsource & rhs);

    ClassDef(AliGenFLUKAsource,1) //Boundary source
};
#endif






