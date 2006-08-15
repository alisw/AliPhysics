#ifndef ALILHCPROCESSBT_H
#define ALILHCPROCESSBT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "AliLhcProcess.h"
class TList;
class AliLhcBeam;

class AliLhcProcessBT : public AliLhcProcess
{
 public:
    AliLhcProcessBT(AliLHC* lhc, const char* name, const char* title);
    AliLhcProcessBT(const AliLhcProcessBT& bt);
    virtual ~AliLhcProcessBT();
    virtual void SetCrossSection(Float_t sig) {fCrossSection = sig*1.e-24;}
    virtual void Init();
    virtual void Evolve(Float_t dt);
    virtual void SetBetaMin(Float_t b) {fBetaMin = b;}
    
    AliLhcProcessBT & operator=(const AliLhcProcessBT & rhs);
 protected:
    Float_t fCrossSection; // Interaction cross section 
    TList * fIRegions;     // Interaction Regions
    Float_t fBetaMin;      // Minimal allowed beta*
//
    ClassDef(AliLhcProcessBT,1) // LHC Process: Beta Tuning 
};

#endif
