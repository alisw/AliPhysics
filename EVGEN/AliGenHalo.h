#ifndef ALIGENHALO_H
#define ALIGENHALO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"
#include <TString.h>

// Read beam halo background particles from a boundary source
// Boundary source is in the LHCb format
// http://www.hep.manchester.ac.uk/u/robert/LHC_backgrounds/Note-MIBStudies.pdf
// and has been provided by Robert Appleby
// Author: andreas.morsch@cern.ch

class AliGenHalo : public AliGenerator
{
public:
    enum constants{kY1Day0, kY1Day70, kY2D0, kY2D10, kY3D90, kLHCPR674Startup, kLHCPR674Conditioned};
    AliGenHalo();
    AliGenHalo(Int_t npart);
    virtual ~AliGenHalo();
    virtual void Init();
    virtual void SetFileName(TString filename) {fFileName=TString(filename);}
    virtual void Generate();
    virtual Float_t GasPressureWeight(Float_t zPrimary);
    virtual void SetSide(Int_t flag = 1) {fSide = flag;}
    virtual void SetNskip(Int_t nskip) {fNskip = nskip;}
    virtual void SetRunPeriod(Int_t t = kY3D90) {fRunPeriod = t;}
    virtual void SetTimePerEvent(Float_t t = 1.e-4) {fTimePerEvent = t;}
    virtual void Draw(Option_t * opt="");
    virtual void  CountEvents();
 private:
    virtual void  SkipEvents();
    virtual Int_t ReadNextParticle();
 protected:
    FILE*    fFile;                       // ! Pointer to file
    TString  fFileName;                   //   Choose the file
    Int_t    fSide;                       //   Muon arm side (1) / Castor side (-1)
    Int_t    fRunPeriod;                  //   LHC Running Period
    Float_t  fTimePerEvent;               //   Time corresponding to one event [s]
    Int_t    fNskip;                      //   Number of entries to skip
    Float_t* fZ1;                         // ! z-positions for gas pressure tables
    Float_t* fZ2;                         // ! z-positions for gas pressure tables 
    Float_t* fG1;                         // ! gas pressures
    Float_t* fG2;                         // ! gas pressures
    Int_t    fGPASize;                    // ! Size of arrays
    Int_t    fLossID;                     // ! unique loss ID
    Int_t    fLossA;                      // ! atomic number of scatterer
    Int_t    fPdg;                        // ! pdg code 
    Float_t  fLossT0;                     // ! relative time
    Float_t  fLossZ;                      // ! scaterring position 
    Float_t  fLossW;                      // ! weight of proton loss
    Float_t  fXS;                         // ! x-position on scoring plane 
    Float_t  fYS;                         // ! y-position on scoring plane
    Float_t  fZS;                         // ! z-position on scoring plane
    Float_t  fDX;                         // ! direction cosine x
    Float_t  fDY;                         // ! direction cosine y
    Float_t  fEkin;                       // ! kinetic energy
    Float_t  fTS;                         // ! relative arrival time
    Float_t  fWS;                         // ! weight
 private:
    AliGenHalo(const AliGenHalo &Halo);
    AliGenHalo & operator=(const AliGenHalo & rhs);
    
    ClassDef(AliGenHalo,1)        // LHC beam halo boundary source
};
#endif






