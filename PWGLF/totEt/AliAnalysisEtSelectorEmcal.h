//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selection class for EMCAL
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#ifndef ALIANALYSISETSELECTOREMCAL_H
#define ALIANALYSISETSELECTOREMCAL_H

#include "AliAnalysisEtSelector.h"


class AliAnalysisEtSelectorEmcal : public AliAnalysisEtSelector
{

public:

    AliAnalysisEtSelectorEmcal(AliAnalysisEtCuts* cuts);
    
    virtual ~AliAnalysisEtSelectorEmcal();
    
    virtual TRefArray* GetClusters();
    virtual Bool_t CutMinEnergy(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutMinEnergy(const TParticle& part) const;
    virtual Bool_t CutDistanceToBadChannel(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutTrackMatching(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutGeometricalAcceptance(const TParticle& part) const;    
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack& part) const;    
    virtual void Init();
    virtual Int_t Init(const AliESDEvent *ev);
    
    virtual void SetEvent(const AliESDEvent* event);

private:
  
    Bool_t fInitialized; // matrix initialized

    AliAnalysisEtSelectorEmcal(); // Prohibited
    AliAnalysisEtSelectorEmcal(const AliAnalysisEtSelectorEmcal& other); // Prohibited
    AliAnalysisEtSelectorEmcal& operator=(const AliAnalysisEtSelectorEmcal& other); // Prohibited
    bool operator==(const AliAnalysisEtSelectorEmcal& other) const; // Prohibited
    
    ClassDef(AliAnalysisEtSelectorEmcal, 1);
};

#endif // ALIANALYSISETSELECTOREMCAL_H
