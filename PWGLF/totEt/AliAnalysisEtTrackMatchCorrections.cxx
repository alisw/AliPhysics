

#include "AliAnalysisEtTrackMatchCorrections.h"

ClassImp(AliAnalysisEtTrackMatchCorrections);

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections() : TNamed("TMCorr","TMCorr")
    ,fChargedContr(new TF1)
    ,fNeutralContr(new TF1)
    ,fGammaContr(new TF1)
    ,fSecondaryContr(new TF1)
    ,fMeanCharged(0)
    ,fMeanNeutral(0)
    ,fMeanGamma(0)
    ,fMeanSecondary(0)
{}

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections(const TString name, const TF1 &chargedContr, const TF1 &neutralContr, const TF1 &gammaContr, const TF1 &secondaryContr,
								       Double_t meanCharged, Double_t meanNeutral, Double_t meanGammas, Double_t meanSecondary) : TNamed(name,name)
    ,fChargedContr(new TF1(chargedContr))
    ,fNeutralContr(new TF1(neutralContr))
    ,fGammaContr(new TF1(gammaContr))
    ,fSecondaryContr(new TF1(secondaryContr))
    ,fMeanCharged(meanCharged)
    ,fMeanNeutral(meanNeutral)
    ,fMeanGamma(meanGammas)
    ,fMeanSecondary(meanSecondary)
{}

//! Copy constructor
AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections(const AliAnalysisEtTrackMatchCorrections &obj) : TNamed(obj)
    ,fChargedContr(new TF1(*(obj.fChargedContr)))
    ,fNeutralContr(new TF1(*(obj.fNeutralContr)))
    ,fGammaContr(new TF1(*(obj.fGammaContr)))
    ,fSecondaryContr(new TF1(*(obj.fSecondaryContr)))
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
        *fChargedContr = *(other.fChargedContr);
        *fNeutralContr = *(other.fNeutralContr);
        *fGammaContr = *(other.fGammaContr);
        *fSecondaryContr = *(other.fSecondaryContr);
	fMeanCharged = other.fMeanCharged;
	fMeanNeutral = other.fMeanNeutral;
	fMeanGamma = other.fMeanGamma;
	fMeanSecondary = other.fMeanSecondary;
	
    }
    return *this;
}

