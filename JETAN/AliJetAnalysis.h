#ifndef ALIJETANALYSIS_H
#define ALIJETANALYSIS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// JetAnalysis class 
// Perform Jet Analysis
// Author: amdreas.morsch@cern.ch
//---------------------------------------------------------------------
#include <TObject.h> 
class AliJetAnalysis : public TObject
{
 public:
 
    AliJetAnalysis();
    ~AliJetAnalysis(){;}

    void Analyse();
    // Setter
    void SetDirectory(char* directory) {fDirectory = directory;}
    void SetEventRange(Int_t imin, Int_t imax) {fEventMin = imin; fEventMax = imax;}
    void SetRunRange(Int_t imin, Int_t imax) {fRunMin = imin; fRunMax = imax;}
 private:
    char*  fDirectory;   // Directory
    Int_t  fEventMin;    // Minimum event number
    Int_t  fEventMax;    // Maximum event number
    Int_t  fRunMin;      // Minimum run number 
    Int_t  fRunMax;      // Maximum run number
    
	
    ClassDef(AliJetAnalysis,1)
};
 
#endif
