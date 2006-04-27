#ifndef ALILHCPROCESSIBS_H
#define ALILHCPROCESSIBS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Realisation of AliLhcProcess for the fast simulation of the
// Intra Beam Scattering process
// in transverse and longitudinal direction.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include "AliLhcProcess.h"
class TList;
class AliLhcBeam;

class AliLhcProcessIBS : public AliLhcProcess
{
 public:
    AliLhcProcessIBS(AliLHC* lhc, const char* name, const char* title);
    virtual ~AliLhcProcessIBS();
    virtual void SetCrossSection(Float_t sig) {fCrossSection = sig*1.e-24;}
    virtual void Init();
    virtual void Evolve(Float_t dt);
    virtual void  SetMonitor(Int_t n);
    virtual void  Record();
    virtual void  DrawPlots();
    AliLhcProcessIBS & operator=(const AliLhcProcessIBS & rhs);
 protected:
    Float_t fCrossSection; // Interaction cross section 
    TList * fIRegions;     // Interaction Regions
    AliLhcBeam* fBeam[2];  // Beams
    Float_t fR[2];           // elem. ion radius
    Float_t fE[2];           // ion radius
    Float_t fTaux;           // transverse time constant 
    Float_t fTaue;           // longitudinal time constant       
    Float_t* fTauxArray;     // [fNmax]
    Float_t* fTaueArray;     // [fNmax]
    
//
    ClassDef(AliLhcProcessIBS,1) // LHC Process: Intra Beam Scattering
};

#endif
