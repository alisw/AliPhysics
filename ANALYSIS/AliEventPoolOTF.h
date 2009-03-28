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
	{fMultMin = min; fMultMax = max; fMultStep = step;}
    //
    void SetTagDirectory(const char* dirname) {fTagDirectory = dirname;};
    virtual Int_t BinNumber() const {return fBinNumber;}
	    
 private:
    AliEventPoolOTF(const AliEventPoolOTF& obj);
    AliEventPoolOTF& operator=(const AliEventPoolOTF& other);
 protected:
    AliTagAnalysis*      fTagAnalysis;  // Pointer to tag analysis
    AliRunTagCuts*       fRunCuts;      // Run      cuts
    AliLHCTagCuts*       fLHCCuts;      // LHC      cuts
    AliDetectorTagCuts*  fDetectorCuts; // Detector cuts
    AliEventTagCuts*     fEventCuts;    // Event    cuts
    const char*                fTagDirectory; // Directory with local tag files
    Int_t                fMultMin;      // Minimum multiplicity
    Int_t                fMultMax;      // Maximum multiplicity
    Int_t                fMultStep;     // Multiplicity step-size 
    Int_t                fMultiplicity; // Current multiplicity
    Int_t                fBinNumber;    // Current bin number
    
    ClassDef(AliEventPoolOTF, 0); 
};
 
#endif
