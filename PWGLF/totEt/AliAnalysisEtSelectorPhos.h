#ifndef ALIANALYSISETSELECTORPHOS_H
#define ALIANALYSISETSELECTORPHOS_H

#include "AliAnalysisEtSelector.h"

class TH2I;
class AliPHOSGeometry;

class AliAnalysisEtSelectorPhos : public AliAnalysisEtSelector
{

public:

    AliAnalysisEtSelectorPhos(AliAnalysisEtCuts *cuts);
    virtual ~AliAnalysisEtSelectorPhos();
    
    virtual TRefArray* GetClusters();
    virtual Bool_t CutMinEnergy(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutDistanceToBadChannel(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutTrackMatching(const AliESDCaloCluster& cluster, Double_t &r) const;
    
    virtual Int_t Init(int runNumber);
  
private:


    int LoadGeometry();
    int LoadBadMaps();
    
    AliPHOSGeometry *fGeoUtils;
    
    TH2I *fBadMapM2; // Bad map
    TH2I *fBadMapM3; // Bad map
    TH2I *fBadMapM4; // Bad map
    
    AliAnalysisEtSelectorPhos();
    AliAnalysisEtSelectorPhos(const AliAnalysisEtSelectorPhos& other); // Prohibited
    AliAnalysisEtSelectorPhos& operator=(const AliAnalysisEtSelectorPhos& other); // Prohibited
    bool operator==(const AliAnalysisEtSelectorPhos& other) const; // Prohibited
    
    ClassDef(AliAnalysisEtSelectorPhos, 1);
};

#endif // ALIANALYSISETSELECTORPHOS_H
