#ifndef ALIANALYSISETTRACKMATCHCORRECTIONS_H
#define ALIANALYSISETTRACKMATCHCORRECTIONS_H


#include "TNamed.h"
#include "TF1.h"

class AliAnalysisEtTrackMatchCorrections : public TNamed
{

public:

//! Default constructor
    AliAnalysisEtTrackMatchCorrections();

//! Constructor
    AliAnalysisEtTrackMatchCorrections(TString name, TF1 chargedContr, TF1 neutralContr, TF1 gammaContr, TF1 secondaryContr, Double_t meanCharged, Double_t meanNeutral, Double_t meanGammas, Double_t meanSecondary);

//! Copy constructor
    AliAnalysisEtTrackMatchCorrections(const AliAnalysisEtTrackMatchCorrections &obj);

//! Destructor
    virtual ~AliAnalysisEtTrackMatchCorrections();

//! Assignment operator
    AliAnalysisEtTrackMatchCorrections& operator=(const AliAnalysisEtTrackMatchCorrections& other);

// Getters

    TF1 ChargedContr() const {
        return fChargedContr;
    }

    TF1 NeutralContr() const {
        return fNeutralContr;
    }

    TF1 GammaContr() const {
        return fGammaContr;
    }

    TF1 SecondaryContr() const {
        return fSecondaryContr;
    }
    
    Double_t ChargedContr(int mult) const {
        return fChargedContr.Eval(mult)*fMeanCharged;
    }

    Double_t NeutralContr(int mult) const {
        return fNeutralContr.Eval(mult)*fMeanNeutral;
    }

    Double_t GammaContr(int mult) const {
        return -fGammaContr.Eval(mult)*fMeanGamma;
    }

    Double_t SecondaryContr(int mult) const {
        return fSecondaryContr.Eval(mult)*fMeanSecondary;
    }


// Setters

    void SetChargedcontr(TF1 chargedContr) {
        fChargedContr = chargedContr;
    }

    void SetNeutralcontr(TF1 neutralContr) {
        fNeutralContr = neutralContr;
    }

    void SetGammacontr(TF1 gammaContr) {
        fGammaContr = gammaContr;
    }

    void SetSecondarycontr(TF1 secondaryContr) {
        fSecondaryContr = secondaryContr;
    }


private:

    // ChargedContr
    TF1 fChargedContr;
    // NeutralContr
    TF1 fNeutralContr;
    // GammaContr
    TF1 fGammaContr;	
    // SecondaryContr
    TF1 fSecondaryContr;
    
    // Mean deposited energy from charged particles
    Double_t fMeanCharged;
    // Mean deposited energy from neutral particles
    Double_t fMeanNeutral;
    // Mean deposited energy from gammas 
    Double_t fMeanGamma;
    // Mean deposited energy from secondaries 
    Double_t fMeanSecondary;

    // Prohibited
//! Equality operator
    bool operator==(const AliAnalysisEtTrackMatchCorrections &other) const;
    ClassDef(AliAnalysisEtTrackMatchCorrections, 1);
};

#endif //ALIANALYSISETTRACKMATCHCORRECTIONS_H

