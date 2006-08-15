#ifndef ALILHCPROCESSBB_H
#define ALILHCPROCESSBB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "AliLhcProcess.h"
class TList;
class AliLhcBeam;

class AliLhcProcessBB : public AliLhcProcess
{
 public:
    AliLhcProcessBB(AliLHC* lhc, const char* name, const char* title);
    AliLhcProcessBB(const AliLhcProcessBB& bb);
    virtual ~AliLhcProcessBB();
    virtual void SetCrossSection(Float_t sig) {fCrossSection = sig*1.e-24;}
    virtual void Init();
    virtual void Evolve(Float_t dt);
    
    AliLhcProcessBB & operator=(const AliLhcProcessBB & rhs);
 protected:
    Float_t fCrossSection; // Interaction cross section 
    TList * fIRegions;     // Interaction Regions
    AliLhcBeam* fBeam1;    // Beam1
    AliLhcBeam* fBeam2;    // Beam2    

//
    ClassDef(AliLhcProcessBB,1) // LHC Process: Beam-Beam 
};

#endif
