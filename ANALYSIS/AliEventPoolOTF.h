#ifndef ALIEVENTPOOLOTF_H
#define ALIEVENTPOOLOTF_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Realisation of an AliVEventPool via
// on the flight (OTF) generation of the bin using AliTagAnalysis.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include <AliVEventPool.h>
class AliRunTagCuts;
class AliLHCTagCuts;
class AliDetectorTagCuts;
class AliEventTagCuts;
class AliTagAnalysis;

typedef enum {kMultiplicity, kZVertex, kEventPlane, kLeadingParticleEta, kUser1, kUser2}  EventPoolAxis_t;

class AliEventPoolOTF : public AliVEventPool
{
 public:
    AliEventPoolOTF();
    AliEventPoolOTF(const char* name, const char* title = "AOD");

    virtual ~AliEventPoolOTF() {;}
    // Interface
    virtual TChain* GetNextChain();
    virtual void  GetCurrentBin(Float_t* /*bin*/);
    virtual Int_t GetDimension();
    virtual void  Init();
    virtual void  SetMultiplicityBin(Int_t min, Int_t max, Int_t step)
	{fValueMin[kMultiplicity] = Float_t(min); fValueMax[kMultiplicity] = Float_t(max); fValueStep[kMultiplicity] = Float_t(step);}

    virtual void  SetMultiplicityBinning(Float_t min, Float_t max, Float_t step)
	{fValueMin[kMultiplicity] = min; fValueMax[kMultiplicity] = max; fValueStep[kMultiplicity] = step;}
    virtual void  SetZVertexBinning(Float_t min, Float_t max, Float_t step)
	{fValueMin[kZVertex] = min; fValueMax[kZVertex] = max; fValueStep[kZVertex] = step;}
    virtual void  SetEventPlaneBinning(Float_t min, Float_t max, Float_t step)
	{fValueMin[kEventPlane] = min; fValueMax[kEventPlane] = max; fValueStep[kEventPlane] = step;}
    virtual void  SetLeadingParticleEtaBinning(Float_t min, Float_t max, Float_t step)
	{fValueMin[kLeadingParticleEta] = min; fValueMax[kLeadingParticleEta] = max; fValueStep[kLeadingParticleEta] = step;}

    //
    void SetTagDirectory(const char* dirname) {fTagDirectory = dirname;};
    virtual Int_t BinNumber() const {return fBinNumber;}
	    
 private:
    AliEventPoolOTF(const AliEventPoolOTF& obj);
    AliEventPoolOTF& operator=(const AliEventPoolOTF& other);
    void InitArrays();
    
 protected:
    AliTagAnalysis*      fTagAnalysis;  // Pointer to tag analysis
    AliRunTagCuts*       fRunCuts;      // Run      cuts
    AliLHCTagCuts*       fLHCCuts;      // LHC      cuts
    AliDetectorTagCuts*  fDetectorCuts; // Detector cuts
    AliEventTagCuts*     fEventCuts;    // Event    cuts
    const char*          fTagDirectory; // Directory with local tag files
    // Common pool cuts
    // Multiplicity
    Float_t              fValueMin[6];  // Minimum value
    Float_t              fValueMax[6];  // Maximum value
    Float_t              fValueStep[6]; // Step size
    Float_t              fValue[6];     // Current value
    //
    Int_t                fBinNumber;    // Current bin number
    Bool_t               fNoMore;       // No more bins 
    
    ClassDef(AliEventPoolOTF, 0); 
};
 
#endif
