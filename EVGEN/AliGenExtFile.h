#ifndef ALIGENEXTFILE_H
#define ALIGENEXTFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Event generator that can read the old ALICE event format based on CW-ntuples
// http://consult.cern.ch/alice/Internal_Notes/1995/32/abstract
// Author: andreas.morsch@cern.ch

#include "AliGenerator.h"
class TTree;

class AliGenExtFile : public AliGenerator
{
 public:
    AliGenExtFile();
    AliGenExtFile(Int_t npart);
    AliGenExtFile(const AliGenExtFile &cocktail);
    
    virtual ~AliGenExtFile();
    // Initialise 
    virtual void Init() {}
    // Initialise fluka data 
    virtual void NtupleInit();
    // set file name of data file
    virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
    // generate event
    virtual void Generate();
    AliGenExtFile & operator=(const AliGenExtFile & rhs);
    enum Code_t {kPDG, kGEANT3};
    void SetParticleCode(Code_t code) {fCode = code;}
    
protected:
    const Text_t     *fFileName;      //! Choose the file
    Int_t             fNcurrent;      // points to the next entry
    TTree            *fTreeNtuple;    // pointer to the TTree
    //Declaration of leaves types
    Code_t          fCode;            // Particle code type
    Int_t           fNihead;          // Number of entries in integer header  
    Int_t           fIhead[12];       // Integer header
    Int_t           fNrhead;          // Number of entries in float header
    Float_t         fRhead[6];        // Float header
    UInt_t          fIdpart;          // Particle type
    Float_t         fTheta;           // Theta 
    Float_t         fPhi;             // Phi
    Float_t         fP;               // Total momentum
    Float_t         fE;               // Total energy
    
  ClassDef(AliGenExtFile,1) //Generate particles from external file
};
#endif






