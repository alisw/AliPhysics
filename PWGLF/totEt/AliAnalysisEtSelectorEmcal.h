//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selection class for EMCAL
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________


#ifndef ALIANALYSISETSELECTOREMCAL_H
#define ALIANALYSISETSELECTOREMCAL_H

#include "AliAnalysisEtSelector.h"

class AliStack;

class AliAnalysisEtSelectorEmcal : public AliAnalysisEtSelector
{

public:

    AliAnalysisEtSelectorEmcal(AliAnalysisEtCuts* cuts);
    AliAnalysisEtSelectorEmcal();
    
    virtual ~AliAnalysisEtSelectorEmcal();
    
    virtual TRefArray* GetClusters();
    virtual Bool_t PassMinEnergyCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t PassMinEnergyCut(const TParticle& part) const;
    virtual Bool_t PassMinEnergyCut(Double_t e) const;
    virtual Bool_t PassDistanceToBadChannelCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t PassTrackMatchingCut(const AliESDCaloCluster& cluster) const;
    virtual Bool_t CutGeometricalAcceptance(const TParticle& part);    
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack& part);        
    virtual Bool_t CutGeometricalAcceptance(const AliESDCaloCluster& cluster);    
    virtual void Init();
    virtual Int_t Init(const AliESDEvent *ev);
    
    virtual Bool_t IsDetectorCluster(const AliESDCaloCluster& cluster) const { return cluster.IsEMCAL(); }

    virtual UInt_t GetLabel(const AliESDCaloCluster *cluster, AliStack *stack);
private:
  
    
    Double_t CalcTrackClusterDistance(const Float_t clsPos[3], Int_t* trkMatchId) const;
    
    //AliAnalysisEtSelectorEmcal(); // Prohibited
    AliAnalysisEtSelectorEmcal(const AliAnalysisEtSelectorEmcal& other); // Prohibited
    AliAnalysisEtSelectorEmcal& operator=(const AliAnalysisEtSelectorEmcal& other); // Prohibited
    bool operator==(const AliAnalysisEtSelectorEmcal& other) const; // Prohibited
    
    ClassDef(AliAnalysisEtSelectorEmcal, 1);
};

#endif // ALIANALYSISETSELECTOREMCAL_H
