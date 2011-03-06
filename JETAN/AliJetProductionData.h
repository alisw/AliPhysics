#ifndef ALIJETPRODUCTIONDATA_H
#define ALIJETPRODUCTIONDATA_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Service class for jet production data 
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//---------------------------------------------------------------------

#include <TObject.h>
#include <TString.h>
 
class AliJetProductionData : public TObject
{
 public:
    AliJetProductionData();
    ~AliJetProductionData();
    Int_t   GetNumberOfPtHardBins() const {return fNbins;}
    void    GetPtHardLimits(Int_t bin, Float_t& ptmin, Float_t& ptmax);
    TString GetRunTitle(Int_t bin);
    Float_t GetWeight(Int_t bin);
 protected:
    Int_t     fNbins;         // Number of pt_hard bins used in the production
    Float_t*  fPtHardLimits;  //[fNbins+1]
    Float_t*  fWeights;       //[fNbins]
    TString*  fRunTitles;     //[fNbins]
    
 private: 
    AliJetProductionData(const AliJetProductionData& rJetPD);
    AliJetProductionData& operator = (const AliJetProductionData& rjpd);

    ClassDef(AliJetProductionData, 1)
};
 
#endif
