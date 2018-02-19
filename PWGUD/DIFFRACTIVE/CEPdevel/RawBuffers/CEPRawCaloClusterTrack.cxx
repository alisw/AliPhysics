/// ////////////////////////////////////////////////////////////////////////////
///
/// CEP Calorimeter buffer
///
/// structure to hold information on the calorimeter
///
// ____________________________________________________________________________
#include "CEPTrackBuffer.h"
#include "CEPRawCaloClusterTrack.h"

ClassImp(CEPRawCaloClusterTrack)

// ____________________________________________________________________________
CEPRawCaloClusterTrack::CEPRawCaloClusterTrack()
  : TObject()
  , fE(0.0)
  , fShapeDisp(0.0)
  , fChi2(CEPTrackBuffer::kdumval)
  , fCaloCpvDist(CEPTrackBuffer::kdumval)
  , fIsPHOS(kFALSE)  
  , fIsEMCAL(kFALSE)  
{

}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::Reset()
{
    fE           = 0.0;
    fShapeDisp   = 0.0;
    fChi2        = CEPTrackBuffer::kdumval;
    fCaloCpvDist = CEPTrackBuffer::kdumval;
    fIsPHOS      = kFALSE;
    fIsEMCAL     = kFALSE; 
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::SetCaloClusterVariables(AliESDCaloCluster* ClusterObj)
{
    // Global member variable setter
    this->SetCaloClusterE(ClusterObj->E());
    this->SetCaloClusterShapeDisp(ClusterObj->GetDispersion());
    this->SetCaloClusterChi2(ClusterObj->Chi2());
    this->SetCaloClusterCPVDist(ClusterObj->GetEmcCpvDistance());
    this->SetCaloClusterIsEMCAL(ClusterObj->IsEMCAL());
    this->SetCaloClusterIsPHOS(ClusterObj->IsPHOS());
}

// ____________________________________________________________________________
