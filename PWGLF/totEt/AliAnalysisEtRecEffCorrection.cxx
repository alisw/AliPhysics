//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selection class for EMCAL
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#include "AliAnalysisEtRecEffCorrection.h"

ClassImp(AliAnalysisEtRecEffCorrection);

AliAnalysisEtRecEffCorrection::AliAnalysisEtRecEffCorrection() : TNamed("RecEff","RecEff")
							       ,fEnergyCorrection()//"ReCorr","1",0.01)
    ,fMaxEnergy(0)
{
  fEnergyCorrection.SetName("ReCorr");
}

AliAnalysisEtRecEffCorrection::AliAnalysisEtRecEffCorrection(TString name, const TF1 &correction, const Double_t maxEnergy) : TNamed(name, name)
    ,fEnergyCorrection(correction)
    ,fMaxEnergy(maxEnergy)
{}

//! Copy constructor
AliAnalysisEtRecEffCorrection::AliAnalysisEtRecEffCorrection(const AliAnalysisEtRecEffCorrection &obj) : TNamed(obj)
    ,fEnergyCorrection(obj.fEnergyCorrection)
    ,fMaxEnergy(obj.fMaxEnergy)
{}

//! Destructor
AliAnalysisEtRecEffCorrection::~AliAnalysisEtRecEffCorrection()
{
}

//! Assignment operator
AliAnalysisEtRecEffCorrection& AliAnalysisEtRecEffCorrection::operator=(const AliAnalysisEtRecEffCorrection &other)
{
    if (this != &other)
    {
        fEnergyCorrection = other.fEnergyCorrection;
        fMaxEnergy = other.fMaxEnergy;
    }
    return *this;
}

//! Equality operator
bool AliAnalysisEtRecEffCorrection::operator==(const AliAnalysisEtRecEffCorrection &other) const
{
    if (this == &other) return true;
    return false;
    //return (fMaxEnergy == other.fMaxEnergy && fEnergyCorrection == other.fEnergyCorrection);
}

Double_t AliAnalysisEtRecEffCorrection::CorrectedEnergy(Double_t energy)
{
 
  return fEnergyCorrection.Eval(energy);
  
}
