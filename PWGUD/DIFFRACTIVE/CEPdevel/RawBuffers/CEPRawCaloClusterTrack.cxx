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
  , fIsEMCAL(kFALSE)  
  , fIsPHOS(kFALSE)  
  , fM20(CEPTrackBuffer::kdumval)
  , fM02(CEPTrackBuffer::kdumval)
  , fTime(CEPTrackBuffer::kdumval)
{
    this->Reset(); 
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::Reset()
{
    fE            = 0.0;
    fShapeDisp    = 0.0;
    fTime         = -999.;
    fChi2         = -999.;
    fCaloCpvDist  = -999.;
    fIsPHOS       = kFALSE;
    fIsEMCAL      = kFALSE; 
    fGlobalPos[0] = -999.;
    fGlobalPos[1] = -999.;
    fGlobalPos[2] = -999.;
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::SetCaloClusterGlobalPosition(Float_t *x)
{
    if (x){
        fGlobalPos[0] = x[0];
        fGlobalPos[1] = x[1];
        fGlobalPos[2] = x[2];
    } else {
        fGlobalPos[0] = -999.;
        fGlobalPos[1] = -999.;
        fGlobalPos[2] = -999.;
    }
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::GetCaloClusterGlobalPosition(Float_t &x, Float_t &y, Float_t &z)
{
    x = fGlobalPos[0];
    y = fGlobalPos[1];
    z = fGlobalPos[2];
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::SetCaloClusterVariables(AliESDCaloCluster* ClusterObj, AliESDCaloCells* CaloCells)
{
    // Global member variable setter
    this->SetCaloClusterE(ClusterObj->E());
    this->SetCaloClusterShapeDisp(ClusterObj->GetDispersion());
    this->SetCaloClusterChi2(ClusterObj->Chi2());
    this->SetCaloClusterCPVDist(ClusterObj->GetEmcCpvDistance());
    this->SetCaloClusterIsEMCAL(ClusterObj->IsEMCAL());
    this->SetCaloClusterIsPHOS(ClusterObj->IsPHOS());
    this->SetCaloClusterM20(ClusterObj->GetM20());
    this->SetCaloClusterM02(ClusterObj->GetM02());
    Float_t x[3];
    ClusterObj->GetPosition(x);
    this->SetCaloClusterGlobalPosition(x);

    if (ClusterObj->IsEMCAL() && ClusterObj->GetNCells()>0){
        // variable preperation for easier readability
        Double_t cellAmpl_max, cellAmpl;
        Short_t cellNb, cellNb_max_ampl;
        // minimum cell ampl = 0
        cellAmpl_max = 0.;
        // Looping thru all cells that fired 
        for (UInt_t kk(0); kk<ClusterObj->GetNCells(); kk++) {
            cellNb = ClusterObj->GetCellAbsId(kk);
            cellAmpl = CaloCells->GetCellAmplitude(cellNb);
            if (cellAmpl>=cellAmpl_max) { cellAmpl_max = cellAmpl; cellNb_max_ampl = cellNb; }
        }
        fTime = CaloCells->GetCellTime(cellNb_max_ampl);
    } else fTime = -999.;
}

// ____________________________________________________________________________
