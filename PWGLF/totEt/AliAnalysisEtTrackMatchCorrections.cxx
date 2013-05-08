

#include "AliAnalysisEtTrackMatchCorrections.h"

ClassImp(AliAnalysisEtTrackMatchCorrections);

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections() : TNamed("TMCorr","TMCorr")
    ,fChargedContr(new TF1)
    ,fNeutralContr(new TF1)
    ,fGammaContr(new TF1)
    ,fSecondaryContr(new TF1)
									 ,fRecoEff(0)
    ,fMeanCharged(0)
    ,fMeanNeutral(0)
    ,fMeanGamma(0)
    ,fMeanSecondary(0)
{}

AliAnalysisEtTrackMatchCorrections::AliAnalysisEtTrackMatchCorrections(const TString name, const TF1 &chargedContr, const TF1 &neutralContr, const TF1 &gammaContr, const TF1 &secondaryContr, const TH2F &recoEff,
								       Double_t meanCharged, Double_t meanNeutral, Double_t meanGammas, Double_t meanSecondary) : TNamed(name,name)
																				,fChargedContr(new TF1(chargedContr))
																				,fNeutralContr(new TF1(neutralContr))
																				,fGammaContr(new TF1(gammaContr))
																				,fSecondaryContr(new TF1(secondaryContr))										    	,fRecoEff(new TH2F(recoEff))
								      //,fRecoEff(0)
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
    ,fRecoEff(new TH2F(*(obj.fRecoEff)))
    ,fMeanCharged(obj.fMeanCharged)
    ,fMeanNeutral(obj.fMeanNeutral)
    ,fMeanGamma(obj.fMeanGamma)
    ,fMeanSecondary(obj.fMeanSecondary)
    
{}

//! Destructor
AliAnalysisEtTrackMatchCorrections::~AliAnalysisEtTrackMatchCorrections()
{
//   delete fChargedContr;
//   delete fNeutralContr;
//   delete fGammaContr;
//   delete fSecondaryContr;
//   delete fRecoEff;
}

//! Assignment operator
AliAnalysisEtTrackMatchCorrections& AliAnalysisEtTrackMatchCorrections::operator=(const AliAnalysisEtTrackMatchCorrections &other)
{
    if (this != &other)
    {
        *fChargedContr = *(other.fChargedContr);
        *fNeutralContr = *(other.fNeutralContr);
        *fGammaContr = *(other.fGammaContr);
	;        *fSecondaryContr = *(other.fSecondaryContr);
	fMeanCharged = other.fMeanCharged;
	fMeanNeutral = other.fMeanNeutral;
	fMeanGamma = other.fMeanGamma;
	fMeanSecondary = other.fMeanSecondary;
	
    }
    return *this;
}

Double_t AliAnalysisEtTrackMatchCorrections::TrackMatchingEfficiency(float pT, int mult) const{
  Double_t eff = 1.0;
  if(fRecoEff) eff =  fRecoEff->GetBinContent(fRecoEff->GetXaxis()->FindBin(pT),fRecoEff->GetXaxis()->FindBin(mult));
  //cout <<"eff "<<eff<<endl;
  return eff;
}
