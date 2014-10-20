//_________________________________________________________________________
//  Utility Class for transverse energy studies
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#ifndef ALIANALYSISETRECEFFCORRECTION_H
#define ALIANALYSISETRECEFFCORRECTION_H

#include "Rtypes.h"
#include "TArrayD.h"
#include "TNamed.h"
#include "TF1.h"
#include "TH2F.h"

class AliAnalysisEtRecEffCorrection : public TNamed
{

public:

//! Default constructor
    AliAnalysisEtRecEffCorrection();

//! Constructor
    AliAnalysisEtRecEffCorrection(TString name, const TF1& correction, const TH2F &recoEff, Double_t maxEnergy);

//! Copy constructor
    AliAnalysisEtRecEffCorrection(const AliAnalysisEtRecEffCorrection &obj);

//! Destructor
    virtual ~AliAnalysisEtRecEffCorrection();

//! Assignment operator
    AliAnalysisEtRecEffCorrection& operator=(const AliAnalysisEtRecEffCorrection& other);

//! Equality operator
    bool operator==(const AliAnalysisEtRecEffCorrection &other) const;

// Getters

    TF1 EnergyCorrection() const {
        return *fEnergyCorrection;
    }


    TH2F *ReconstructionEfficiency() const {
      return fRecoEff;
    }
    Double_t ReconstructionEfficiency(float energy, int cent) const;

    Double_t MaxEnergy() const {
        return fMaxEnergy;
    }

// Setters

    void SetCorrections(const TF1 &corrections) {
        *fEnergyCorrection = corrections;
    }

    void SetMaxenergy(Double_t maxEnergy) {
        fMaxEnergy = maxEnergy;
    }

    void SetReconstructionEfficiency(const TH2F &recoEff) {
        *fRecoEff = recoEff;
    }


    Double_t CorrectedEnergy(Double_t energy); // Calculate corrected cluster E_T 
    Double_t CorrectedEnergy(Double_t energy, int cent); // Calculate corrected cluster E_T 
    
private:

    // Energy correction function
    TF1 *fEnergyCorrection;//
    TH2F *fRecoEff;//Reconstruction efficiency, x axis = pT, y axis = multiplicity, z = efficiency
    
    
    
// MaxEnergy
    Double_t fMaxEnergy; // MaxEnergy
    
    ClassDef(AliAnalysisEtRecEffCorrection, 1);
};

#endif //ALIANALYSISETRECEFFCORRECTION_H

