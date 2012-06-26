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
    
    // Destructor 
    virtual ~AliAnalysisEtSelector();
    
    // Set the current event
    virtual void SetEvent(const AliESDEvent *event) = 0;// { fEvent = event; }
    
    // Init
    virtual void Init() {} 
    
    // Init with event
    virtual Int_t Init(const AliESDEvent *event) { fRunNumber = event->GetRunNumber(); return 0; }
    
    // Return CaloClusters for the detector
    virtual TRefArray* GetClusters() { return 0; }
    
    // Return true if cluster has energy > cut
    virtual Bool_t CutMinEnergy(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Return true if cluster has energy > cut
    virtual Bool_t CutMinEnergy(const TParticle &/*part*/) const { return true; }
    
    // Cut on distance to bad channel
    virtual Bool_t CutDistanceToBadChannel(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Cut on track matching
    virtual Bool_t CutTrackMatching(const AliESDCaloCluster &/*cluster*/) const { return true; }
    
    // Cut on neutral monte carlo particle
    virtual Bool_t CutNeutralMcParticle(Int_t pIdx, AliStack& s, const TParticlePDG& pdg) const;
    
    // Is it an EM E_T particle
    virtual Bool_t IsEmEtParticle(const Int_t pdgCode) const;
    
    // Does the particle come from an EM E_T primary ?
    virtual Bool_t PrimaryIsEmEtParticle(const Int_t pIdx, AliStack &stack) const;

    // Get the index of primary particle for the particle
    Int_t GetPrimary(const Int_t partIdx, AliStack &stack) const;
    
    // Cut on geometrical acceptance 
    virtual Bool_t CutGeometricalAcceptance(const TParticle &/*part*/) const { return true; }
    
    // Cut on geometrical acceptance 
    virtual Bool_t CutGeometricalAcceptance(const AliVTrack &/*part*/) const { return true; }
    
    
protected:
  
    const AliVEvent *fEvent; // Pointer to current event

    TRefArray *fClusterArray; // Array of clusters

    AliAnalysisEtCuts *fCuts; // Pointer to the cuts object; DS: also in base class?

    Int_t fRunNumber; // run number
    

    
private:

    AliAnalysisEtSelector(); // Prohibited
    AliAnalysisEtSelector(const AliAnalysisEtSelector& other);// Prohibited
    AliAnalysisEtSelector& operator=(const AliAnalysisEtSelector& other);// Prohibited
    bool operator==(const AliAnalysisEtSelector& other) const;// Prohibited
    
    ClassDef(AliAnalysisEtSelector, 1);
};

#endif // ALIANALYSISETSELECTOR_H
