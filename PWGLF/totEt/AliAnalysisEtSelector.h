#ifndef ALIANALYSISETSELECTOR_H
#define ALIANALYSISETSELECTOR_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Selector Base class
//  -  
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________
#include <Rtypes.h>
#include "AliAnalysisEtCommon.h"
#include "AliESDEvent.h"

class AliESDCaloCluster;
class AliVTrack;
class TRefArray;
class AliAnalysisEtCuts;
class TParticle;
class TParticlePDG;
class AliStack;

class AliAnalysisEtSelector : public AliAnalysisEtCommon
{

public:
  
    // Constructor takes cuts object
    AliAnalysisEtSelector(AliAnalysisEtCuts *cuts);
    AliAnalysisEtSelector();
    
    // Destructor 
    virtual ~AliAnalysisEtSelector();
    
    // Set the current event
    virtual void SetEvent(const AliESDEvent *event);
    
    // Init
    virtual void Init() {} 
    
    // Init with event
    virtual Int_t Init(const AliESDEvent *event) { fRunNumber = event->GetRunNumber(); return 0; }
    
    // Return CaloClusters for the detector
    virtual TRefArray* GetClusters() { return 0; }

    virtual Float_t ShiftAngle(Float_t phi);//Always returns an angle in radians between -pi<phi<pi
    
    // Return true if cluster has energy > cut
    virtual Bool_t PassMinEnergyCut(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Return true if cluster has energy > cut
    virtual Bool_t PassMinEnergyCut(const TParticle &/*part*/) const { return true; }

    virtual Bool_t PassMinEnergyCut(Double_t e) const;
    
    // Cut on distance to bad channel
    virtual Bool_t PassDistanceToBadChannelCut(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Cut on track matching
    virtual Bool_t PassTrackMatchingCut(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Cut on neutral monte carlo particle
    virtual Bool_t IsNeutralMcParticle(Int_t pIdx, AliStack& s, const TParticlePDG& pdg) const;
    
    // Is it an EM E_T particle
    virtual Bool_t IsEmEtParticle(const Int_t pdgCode) const;
    
    // Does the particle come from an EM E_T primary ?
    virtual Bool_t PrimaryIsEmEtParticle(const Int_t pIdx, AliStack &stack) const;

    // Get the index of primary particle for the particle
    Int_t GetPrimary(const Int_t partIdx, AliStack &stack) const;
    
    // Cut on geometrical acceptance 
    virtual Bool_t CutGeometricalAcceptance(const TParticle &/*part*/) { return true; }
    
    // Cut on geometrical acceptance 
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack &/*part*/) { return true; }
    
    // Cut on geometrical acceptance 
    virtual Bool_t CutGeometricalAcceptance(const AliESDCaloCluster &/*part*/) { return true; }
    
    // From secondary vertex?
    //virtual Bool_t FromSecondaryInteraction(const TParticle& part, AliStack& stack) const;
    virtual Bool_t FromSecondaryInteraction(Int_t partID, AliStack& stack) const;

    Int_t GetMother(Int_t partID, AliStack& stack) const;
    Bool_t IsFromDetectorCover(Int_t partID, AliStack& stack) const;
    Int_t GetFirstMotherNotFromDetectorCover(Int_t partID, AliStack& stack) const;
    
    // Cluster is in correct detector
    virtual Bool_t IsDetectorCluster(const AliESDCaloCluster &cluster) const = 0;

    // Get correct cluster label - PHOS needs different method
    virtual UInt_t GetLabel(const AliESDCaloCluster *cluster, AliStack *stack){if(!stack){return 0;}else{return TMath::Abs(cluster->GetLabel());}}
    
    AliAnalysisEtCuts * GetCuts() const { return fCuts; }
protected:
  
    const AliVEvent *fEvent; //! Pointer to current event

    TRefArray *fClusterArray; //! Array of clusters

    AliAnalysisEtCuts *fCuts; //! Pointer to the cuts object; DS: also in base class?
    
    Bool_t SuspiciousDecayInChain(const UInt_t suspectMotherPdg, const UInt_t suspectDaughterPdg, const TParticle& part, AliStack& stack) const;
    
    Int_t fRunNumber;

    Bool_t fInitialized; // matrix initialized

private:

    //AliAnalysisEtSelector(); // Prohibited
    AliAnalysisEtSelector(const AliAnalysisEtSelector& other);// Prohibited
    AliAnalysisEtSelector& operator=(const AliAnalysisEtSelector& other);// Prohibited
    bool operator==(const AliAnalysisEtSelector& other) const;// Prohibited
    
    ClassDef(AliAnalysisEtSelector, 1);
};

#endif // ALIANALYSISETSELECTOR_H
