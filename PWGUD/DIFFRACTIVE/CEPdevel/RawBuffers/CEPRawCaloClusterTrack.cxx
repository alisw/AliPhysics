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
  , fPhiEtaDistToNearestTrack(999.)
  , fHasTrackToMatch(kFALSE)
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
    fPhiEtaDistToNearestTrack = 999.;
    fHasTrackToMatch = kFALSE;
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
Bool_t CEPRawCaloClusterTrack::IsClusterFromBG(Double_t cutdPhiEta) const
{
    if (fHasTrackToMatch && fPhiEtaDistToNearestTrack<cutdPhiEta) return kFALSE;
    else return kTRUE;
}

// ____________________________________________________________________________
void CEPRawCaloClusterTrack::SetCaloClusterVariables(AliESDCaloCluster* ClusterObj, AliESDCaloCells* CaloCells, AliESDEvent* ESDobj, TArrayI* TTindices)
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
    // get the time in the most energetic cell in the cluster
    if (ClusterObj->IsEMCAL() && ClusterObj->GetNCells()>0) {
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
    
        // /////////////////////////////////////////////////////////////////////////////
        // -------------- get the distance to the nearest track ------------------------
        // contruct the track-TObjArray here
        //TObjArray* track_arr = new TObjArray();
        //track_arr->SetOwner(kFALSE);
        //AliESDtrack* trk = 0x0;
        //for (Int_t kk(0); kk<ESDobj->GetNumberOfTracks(); kk++){
        //    trk = (AliESDtrack*) ESDobj->GetTrack(kk); 
        //    if (!trk) continue;
        //    track_arr->Add(trk);
        //}
        AliESDtrack *tmptrk = NULL;
        Float_t x[3];
        ClusterObj->GetPosition(x);
        
        // v3 phi is in the range [-pi,pi) -> map it to [0, 2pi)
        TVector3 v3(x[0], x[1], x[2]);
        Double_t cluster_phi = (v3.Phi()>0.) ? v3.Phi() : v3.Phi() + 2.*TMath::Pi();
        Double_t cluster_eta = v3.Eta();
        
        Double_t dPhiEtaMin = 999.;
        for (Int_t kk(0); kk<TTindices->GetSize(); kk++) {
            Int_t trkIndex = TTindices->At(kk);
            tmptrk = (AliESDtrack*) ESDobj->GetTrack(trkIndex); 
            if (!tmptrk) continue;
            
            // track position on emcal
            Double_t trkPhiOnEmc = tmptrk->GetTrackPhiOnEMCal();
            
            // Map phi to [0,2pi)
            trkPhiOnEmc = (trkPhiOnEmc>0.) ? trkPhiOnEmc : trkPhiOnEmc+2.*TMath::Pi();
            if (trkPhiOnEmc<0.) trkPhiOnEmc=-999.;
            Double_t trkEtaOnEmc = tmptrk->GetTrackEtaOnEMCal();
            
            // no matching if at least one has value -999.
            if (trkPhiOnEmc==-999. || trkEtaOnEmc==-999.) continue;
            
            // in case of track on emcal: calculate distance in phi and eta
            Double_t dEta = trkEtaOnEmc - cluster_eta;
            Double_t dPhi = trkPhiOnEmc - cluster_phi;
            Double_t dPhiEta = TMath::Sqrt( dEta*dEta + dPhi*dPhi );

            dPhiEtaMin = (dPhiEtaMin<dPhiEta) ? dPhiEtaMin : dPhiEta;
        }
        fPhiEtaDistToNearestTrack = dPhiEtaMin;
        fHasTrackToMatch = (dPhiEtaMin==999.) ? kFALSE : kTRUE;
        // clear tracks

    } 
}

// ____________________________________________________________________________
