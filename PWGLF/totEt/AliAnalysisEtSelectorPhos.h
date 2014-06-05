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
class AliStack;

class AliAnalysisEtSelectorPhos : public AliAnalysisEtSelector
{

public:

    AliAnalysisEtSelectorPhos(AliAnalysisEtCuts *cuts);
    AliAnalysisEtSelectorPhos();
    virtual ~AliAnalysisEtSelectorPhos();
    
    virtual TRefArray* GetClusters();
    virtual Bool_t PassMinEnergyCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t PassMinEnergyCut(const TParticle& part) const;
    virtual Bool_t PassMinEnergyCut(Double_t e) const;
    virtual Bool_t PassDistanceToBadChannelCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t PassTrackMatchingCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutGeometricalAcceptance(const TParticle& part);    
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack& part);    
    virtual Bool_t CutGeometricalAcceptance(const AliESDCaloCluster& cluster);    
    virtual void Init() {}
    virtual Int_t Init(const AliESDEvent *ev);

    virtual Bool_t IsDetectorCluster(const AliESDCaloCluster& cluster) const {return cluster.IsPHOS();}

     virtual UInt_t GetLabel(const AliESDCaloCluster *cluster, AliStack *stack);
    
private:


    int LoadGeometry(); // load geometry
    int LoadBadMaps(); // load bad maps
    
    AliPHOSGeometry *fGeoUtils; //! geo utils
    
    TH2I *fBadMapM2; //! Bad map
    TH2I *fBadMapM3; //! Bad map
    TH2I *fBadMapM4; //! Bad map

    Bool_t fMatrixInitialized; // matrix initialized
    
    //AliAnalysisEtSelectorPhos();
    AliAnalysisEtSelectorPhos(const AliAnalysisEtSelectorPhos& other); // Prohibited
    AliAnalysisEtSelectorPhos& operator=(const AliAnalysisEtSelectorPhos& other); // Prohibited
    bool operator==(const AliAnalysisEtSelectorPhos& other) const; // Prohibited
    
    ClassDef(AliAnalysisEtSelectorPhos, 1);
};

#endif // ALIANALYSISETSELECTORPHOS_H
