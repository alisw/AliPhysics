#ifndef ALIANALYSISETRECONSTRUCTEDPHOS_H
#define ALIANALYSISETRECONSTRUCTEDPHOS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for PHOS
//  - reconstruction output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtReconstructed.h"


class AliAnalysisEtReconstructedPhos : public AliAnalysisEtReconstructed
{

public:

    AliAnalysisEtReconstructedPhos();
    virtual ~AliAnalysisEtReconstructedPhos();

    virtual void Init();

    void CreateHistograms();

    void SetChargedContributionParameters(Double_t par[3])
    {
        fChargedContributionCorrectionParameters[0] = par[0];
        fChargedContributionCorrectionParameters[1] = par[1];
	fChargedContributionCorrectionParameters[2] = par[2];
    }
    void SetNeutralContributionParameters(Double_t par[3])
    {
        fNeutralContributionCorrectionParameters[0] = par[0];
        fNeutralContributionCorrectionParameters[1] = par[1];
	fNeutralContributionCorrectionParameters[2] = par[2];
    }
    void SetRemovedGammaContributionParameters(Double_t par[3])
    {
        fRemovedGammaContributionCorrectionParameters[0] = par[0];
        fRemovedGammaContributionCorrectionParameters[1] = par[1];
	fRemovedGammaContributionCorrectionParameters[2] = par[2];
    }
    void SetSecondaryContributionParameters(Double_t par[3])
    {
        fSecondaryContributionCorrectionParameters[0] = par[0];
        fSecondaryContributionCorrectionParameters[1] = par[1];
	fSecondaryContributionCorrectionParameters[2] = par[2];
    }
    

protected:

    virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);


    virtual Double_t GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr);//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low
    
private:

    Double_t fChargedContributionCorrectionParameters[3]; // Parametrization of the charged contribution as function of cluster multiplicity
    Double_t fNeutralContributionCorrectionParameters[3]; // Parametrization of the neutral contribution as function of cluster multiplicity
    Double_t fRemovedGammaContributionCorrectionParameters[3]; // Parametrization of the negative contribution from removed gammas as function of cluster multiplicity

    Double_t fSecondaryContributionCorrectionParameters[3]; // Parametrization of the positive contribution of secondary particles


    ClassDef(AliAnalysisEtReconstructedPhos, 1);
};

#endif // ALIANALYSISETRECONSTRUCTEDPHOS_H
