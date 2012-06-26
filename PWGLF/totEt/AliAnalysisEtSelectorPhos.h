#ifndef ALIANALYSISETSELECTORPHOS_H
#define ALIANALYSISETSELECTORPHOS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selector Base class for PHOS
//  - 
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________
#include "AliAnalysisEtSelector.h"

class TH2I;
class TParticle;
class AliPHOSGeometry;
class AliESDEvent;

class AliAnalysisEtSelectorPhos : public AliAnalysisEtSelector
{

public:

    AliAnalysisEtSelectorPhos(AliAnalysisEtCuts *cuts);
    virtual ~AliAnalysisEtSelectorPhos();
    
    virtual TRefArray* GetClusters();
    virtual Bool_t CutMinEnergy(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutMinEnergy(const TParticle& part) const;
    virtual Bool_t CutDistanceToBadChannel(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutTrackMatching(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutGeometricalAcceptance(const TParticle& part) const;    
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack& part) const;    
    virtual void Init() {}
    virtual Int_t Init(const AliESDEvent *ev);

    virtual void SetEvent(const AliESDEvent* event);
    
private:


    int LoadGeometry(); // load geometry
    int LoadBadMaps(); // load bad maps
    
    AliPHOSGeometry *fGeoUtils; // geo utils
    
    TH2I *fBadMapM2; // Bad map
    TH2I *fBadMapM3; // Bad map
    TH2I *fBadMapM4; // Bad map

    Bool_t fInitialized; // matrix initialized
    Bool_t fMatrixInitialized; // matrix initialized
    
    AliAnalysisEtSelectorPhos();
    AliAnalysisEtSelectorPhos(const AliAnalysisEtSelectorPhos& other); // Prohibited
    AliAnalysisEtSelectorPhos& operator=(const AliAnalysisEtSelectorPhos& other); // Prohibited
    bool operator==(const AliAnalysisEtSelectorPhos& other) const; // Prohibited
    
    ClassDef(AliAnalysisEtSelectorPhos, 1);
};

#endif // ALIANALYSISETSELECTORPHOS_H
