

#include "AliAnalysisEtTrackMatchCorrections.h"

ClassImp(AliAnalysisEtTrackMatchCorrections);

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections() : TNamed("TMCorr","TMCorr")
    ,fChargedContr()
    ,fNeutralContr()
    ,fGammaContr()
    ,fSecondaryContr()
    ,fMeanCharged(0)
    ,fMeanNeutral(0)
    ,fMeanGamma(0)
    ,fMeanSecondary(0)
{}

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections(TString name, TF1 chargedContr, TF1 neutralContr, TF1 gammaContr, TF1 secondaryContr,
								       Double_t meanCharged, Double_t meanNeutral, Double_t meanGammas, Double_t meanSecondary) : TNamed(name,name)
    ,fChargedContr(chargedContr)
    ,fNeutralContr(neutralContr)
    ,fGammaContr(gammaContr)
    ,fSecondaryContr(secondaryContr)
    ,fMeanCharged(meanCharged)
    ,fMeanNeutral(meanNeutral)
    ,fMeanGamma(meanGammas)
    ,fMeanSecondary(meanSecondary)
{}

//! Copy constructor
AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections(const AliAnalysisEtTrackMatchCorrections &obj) : TNamed(obj)
    ,fChargedContr(obj.fChargedContr)
    ,fNeutralContr(obj.fNeutralContr)
    ,fGammaContr(obj.fGammaContr)
    ,fSecondaryContr(obj.fSecondaryContr)
    ,fMeanCharged(obj.fMeanCharged)
    ,fMeanNeutral(obj.fMeanNeutral)
    ,fMeanGamma(obj.fMeanGamma)
    ,fMeanSecondary(obj.fMeanSecondary)
    
{}

//! Destructor
AliAnalysisEtTrackMatchCorrections::~AliAnalysisEtTrackMatchCorrections()
{}

//! Assignment operator
AliAnalysisEtTrackMatchCorrections& AliAnalysisEtTrackMatchCorrections::operator=(const AliAnalysisEtTrackMatchCorrections &other)
{
    if (this != &other)
    {
        fChargedContr = other.fChargedContr;
        fNeutralContr = other.fNeutralContr;
        fGammaContr = other.fGammaContr;
        fSecondaryContr = other.fSecondaryContr;
	fMeanCharged = other.fMeanCharged;
	fMeanNeutral = other.fMeanNeutral;
	fMeanGamma = other.fMeanGamma;
	fMeanSecondary = other.fMeanSecondary;
	
    }
    return *this;
}

