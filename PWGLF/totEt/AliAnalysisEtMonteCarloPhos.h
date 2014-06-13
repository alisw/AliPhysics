#ifndef ALIANALYSISETMONTECARLOPHOS_H
#define ALIANALYSISETMONTECARLOPHOS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for PHOS
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtMonteCarlo.h"
class AliPHOSGeoUtils;
class TH2I;

class AliAnalysisEtMonteCarloPhos : public AliAnalysisEtMonteCarlo
{

public:
   
    AliAnalysisEtMonteCarloPhos();
    virtual ~AliAnalysisEtMonteCarloPhos();

    virtual void Init();
    void CreateHistograms();
protected:
      virtual Bool_t TooCloseToBadChannel(const AliESDCaloCluster &cluster) const;

 private:
   
    TH2I *fBadMapM2; //! Bad map
    TH2I *fBadMapM3; //! Bad map
    TH2I *fBadMapM4; //! Bad map

    AliPHOSGeoUtils *fGeoUtils;// // Geo utils
    
    // Prohibited
    AliAnalysisEtMonteCarloPhos & operator = (const AliAnalysisEtMonteCarloPhos&) ;//cpy assignment
    AliAnalysisEtMonteCarloPhos(const AliAnalysisEtMonteCarloPhos&) ; // cpy ctor
   ClassDef(AliAnalysisEtMonteCarloPhos, 1); 
};

#endif // ALIANALYSISETRECONSTRUCTEDPHOS_H
