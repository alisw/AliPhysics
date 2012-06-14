#ifndef ALIANALYSISETSELECTOR_H
#define ALIANALYSISETSELECTOR_H
#include <Rtypes.h>

class AliVEvent;
class AliESDCaloCluster;
class TRefArray;
class AliAnalysisEtCuts;
class AliAnalysisEtSelector
{

public:
  
    // Constructor takes cuts object
    AliAnalysisEtSelector(AliAnalysisEtCuts *cuts);
    
    // Destructor 
    virtual ~AliAnalysisEtSelector();
    
    // Set the current event
    void SetEvent(const AliVEvent *event) { fEvent = event; }
    
    // Init 
    virtual Int_t Init(Int_t runNumber) { fRunNumber = runNumber; return 0; }
    
    // Return CaloClusters for the detector
    virtual TRefArray* GetClusters() = 0;
    
    // Return true if cluster has energy > cut
    virtual Bool_t CutMinEnergy(const AliESDCaloCluster & /*cluster*/) const { return true; }
    
    // Cut on distance to bad channel
    virtual Bool_t CutDistanceToBadChannel(const AliESDCaloCluster & /*cluster*/) const { return true; }
    
    // Cut on track matching
    virtual Bool_t CutTrackMatching(const AliESDCaloCluster& /*cluster*/, Double_t &/*r*/) const { return true; }
  
protected:
  
    const AliVEvent *fEvent; // Pointer to current event

    AliAnalysisEtCuts *fCuts; // Pointer to the cuts object
    
    Int_t fRunNumber;
    
private:

    AliAnalysisEtSelector(); // Prohibited
    AliAnalysisEtSelector(const AliAnalysisEtSelector& other);// Prohibited
    AliAnalysisEtSelector& operator=(const AliAnalysisEtSelector& other);// Prohibited
    bool operator==(const AliAnalysisEtSelector& other) const;// Prohibited
    
    ClassDef(AliAnalysisEtSelector, 1);
};

#endif // ALIANALYSISETSELECTOR_H
