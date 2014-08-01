#ifndef ALIANALYSISETTRACKMATCHCORRECTIONS_H
#define ALIANALYSISETTRACKMATCHCORRECTIONS_H


#include "TNamed.h"
#include "TF1.h"
#include "TH2F.h"

class AliAnalysisEtTrackMatchCorrections : public TNamed
{

public:

//! Default constructor
    AliAnalysisEtTrackMatchCorrections();

//! Constructor
    AliAnalysisEtTrackMatchCorrections(const TString name, const TF1 &chargedContr, const TF1 &neutralContr, const TF1 &gammaContr, const TF1 &secondaryContr, const TH2F &recoEff, Double_t meanCharged, Double_t meanNeutral, Double_t meanGammas, Double_t meanSecondary);

//! Copy constructor
    AliAnalysisEtTrackMatchCorrections(const AliAnalysisEtTrackMatchCorrections &obj);

//! Destructor
    virtual ~AliAnalysisEtTrackMatchCorrections();

//! Assignment operator
    AliAnalysisEtTrackMatchCorrections& operator=(const AliAnalysisEtTrackMatchCorrections& other);

// Getters

    TF1 ChargedContr() const {
        return *fChargedContr;
    }

    TF1 NeutralContr() const {
        return *fNeutralContr;
    }

    TF1 GammaContr() const {
        return *fGammaContr;
    }

    TF1 SecondaryContr() const {
        return *fSecondaryContr;
    }

    TH2F *TrackMatchingEfficiency() const {
      return fRecoEff;
    }
    
    Double_t ChargedContr(int mult) const {
        return fChargedContr->Eval(mult)*fMeanCharged;
    }

    Double_t NeutralContr(int mult) const {
        return fNeutralContr->Eval(mult)*fMeanNeutral;
    }

    Double_t GammaContr(int mult) const {
        return -fGammaContr->Eval(mult)*fMeanGamma;
    }

    Double_t SecondaryContr(int mult) const {
        return fSecondaryContr->Eval(mult)*fMeanSecondary;
    }

    Double_t TrackMatchingEfficiency(Float_t pT, Int_t cent) const;


// Setters

    void SetChargedcontr(const TF1 &chargedContr) {
        *fChargedContr = chargedContr;
    }

    void SetNeutralcontr(const TF1 &neutralContr) {
        *fNeutralContr = neutralContr;
    }

    void SetGammacontr(const TF1 &gammaContr) {
        *fGammaContr = gammaContr;
    }

    void SetSecondarycontr(const TF1 &secondaryContr) {
        *fSecondaryContr = secondaryContr;
    }

    void SetReconstructionEfficiency(const TH2F &recoEff) {
        *fRecoEff = recoEff;
    }

    void SetMinEtCorrection(Int_t cb, Double_t val){fMinEtCorrection[cb] = val;}
    void SetNeutronCorrection(Int_t cb, Double_t val){fNeutronCorrection[cb] = val;}
    void SetHadronCorrection(Int_t cb, Double_t val){fHadronCorrection[cb] = val;}
    void SetKaonCorrection(Int_t cb, Double_t val){fKaonCorrection[cb] = val;}
    void SetSecondaryCorrection(Int_t cb, Double_t val){fSecondaryCorrection[cb] = val;}

    Double_t GetMinEtCorrection(Int_t cb){if(cb>=0 && cb <20){return fMinEtCorrection[cb];}else{return 0.0;}}
    Double_t GetNeutronCorrection(Int_t cb){if(cb>=0 && cb <20){return fNeutronCorrection[cb];}else{return 0.0;}}
    Double_t GetHadronCorrection(Int_t cb){if(cb>=0 && cb <20){return fHadronCorrection[cb];}else{return 0.0;}}
    Double_t GetKaonCorrection(Int_t cb){if(cb>=0 && cb <20){return fKaonCorrection[cb];}else{return 0.0;}}
    Double_t GetSecondaryCorrection(Int_t cb){if(cb>=0 && cb <20){return fSecondaryCorrection[cb];}else{return 0.0;}}

private:

    // ChargedContr
    TF1 *fChargedContr;
    // NeutralContr
    TF1 *fNeutralContr;
    // GammaContr
    TF1 *fGammaContr;	
    // SecondaryContr
    TF1 *fSecondaryContr;

    TH2F *fRecoEff;//Reconstruction efficiency, x axis = pT, y axis = multiplicity, z = efficiency
    
    // Mean deposited energy from charged particles
    Double_t fMeanCharged;
    // Mean deposited energy from neutral particles
    Double_t fMeanNeutral;
    // Mean deposited energy from gammas 
    Double_t fMeanGamma;
    // Mean deposited energy from secondaries 
    Double_t fMeanSecondary;

    Double_t fMinEtCorrection[20];//NeutronCorrection
    Double_t fNeutronCorrection[20];//NeutronCorrection
    Double_t fHadronCorrection[20];//NeutronCorrection
    Double_t fKaonCorrection[20];//NeutronCorrection
    Double_t fSecondaryCorrection[20];//NeutronCorrection

    // Prohibited
//! Equality operator
    bool operator==(const AliAnalysisEtTrackMatchCorrections &other) const;
    ClassDef(AliAnalysisEtTrackMatchCorrections, 2);
};

#endif //ALIANALYSISETTRACKMATCHCORRECTIONS_H

