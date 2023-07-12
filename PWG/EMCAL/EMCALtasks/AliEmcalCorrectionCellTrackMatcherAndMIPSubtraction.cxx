// AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction
//

#include "AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction.h"

#include <iostream>

#include "AliClusterContainer.h"
#include "AliParticleContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVParticle.h"
#include "AliEmcalParticle.h"
#include "AliEMCALGeometry.h"
#include "AliMCEvent.h"



/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction> AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::reg("AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction");

/**
 * Default constructor
 */
AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction"),
  fEmipData(0.2356),
  fEmipMC(0.2824),
  fDoPropagation(0),
  fDoMatchPrimaryOnly(0),
  fMaxDistToVtxPrim(0),
  fSubtractMultiple(0),
  fCellTrackMatchdEtadPhi(0),
  fCellTrackMatchdEtaDiff(0),
  fCellTrackMatchdPhiDiff(0),
  fCellNTrackMatch(0),
  fCellTrackMatchEbefore(0),
  fCellTrackMatchEafter(0),
  fHistoMatchingEfficiency(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::~AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  GetProperty("EmipData", fEmipData);
  GetProperty("EmipMC", fEmipMC);
  GetProperty("DoPropagation", fDoPropagation);
  GetProperty("MatchOnlyPrimaries", fDoMatchPrimaryOnly); // match only primaries in order to be able to fully reject conversions
  GetProperty("MaxDistToVtxPrim", fMaxDistToVtxPrim);
  GetProperty("SubtractMultiple", fSubtractMultiple);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellTrackMatchdEtadPhi = new TH2F("hCellTrackMatchdEtadPhi","hCellTrackMatchdEtadPhi",100,-0.04,0.04,100,-0.04,0.04);
    fOutput->Add(fCellTrackMatchdEtadPhi);

    if(fEsdMode || fDoPropagation){
      fCellTrackMatchdEtaDiff = new TH1D("hCellTrackMatchdEtaDiff","hCellTrackMatchdEtaDiff",100,-0.04,0.04);
      fOutput->Add(fCellTrackMatchdEtaDiff);
      fCellTrackMatchdPhiDiff = new TH1D("hCellTrackMatchdPhiDiff","hCellTrackMatchdPhiDiff",100,-0.04,0.04);
      fOutput->Add(fCellTrackMatchdPhiDiff);
    }

    fCellNTrackMatch = new TH2F("hCellNTrackMatch","hCellNTrackMatch",5,-0.5,4.5,500,0,20);
    fOutput->Add(fCellNTrackMatch);
    fCellTrackMatchEbefore = new TH1F("hCellTrackMatchEbefore","hCellTrackMatchEbefore",500,0,20);
    fOutput->Add(fCellTrackMatchEbefore);
    fCellTrackMatchEafter = new TH1F("hCellTrackMatchEafter","hCellTrackMatchEafter",500,0,20);
    fOutput->Add(fCellTrackMatchEafter);

    fHistoMatchingEfficiency = new TH2F("HistoMatchingEfficiency","HistoMatchingEfficiency",2,-0.5,1.5,500, 0, 50);
    fHistoMatchingEfficiency->GetXaxis()->SetBinLabel(1, "matched");
    fHistoMatchingEfficiency->GetXaxis()->SetBinLabel(2, "not matched");
    fHistoMatchingEfficiency->GetYaxis()->SetTitle("p_{T, track} (GeV/c)");
    fOutput->Add(fHistoMatchingEfficiency);

  }


}

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::ExecOnce()
{
  //fParticleContainerIndexMap.CopyMappingFrom(AliParticleContainer::GetEmcalContainerIndexMap(), fParticleCollArray);
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::Run()
{
  AliEmcalCorrectionComponent::Run();

  AliVTrack* track = 0;


  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    auto partItCont = partCont->accepted_momentum();
    for (AliParticleIterableMomentumContainer::iterator partIterator = partItCont.begin(); partIterator != partItCont.end(); ++partIterator) {
      track = static_cast<AliVTrack *>(partIterator->second);

      if(fDoMatchPrimaryOnly){
        Float_t bPos[2];
        Float_t bCovPos[3];
        track->GetImpactParameters(bPos,bCovPos);
        Float_t dcaToVertexXYPos = bPos[0];
        // Float_t dcaToVertexZPos  = bPos[1];
        // radius has to be larger than 5 cm
        if(dcaToVertexXYPos > fMaxDistToVtxPrim) continue;
      }

      // only propagate tracks for ESDs. Take existing propagation for AODs or force propagation
      if(fEsdMode || fDoPropagation){

        // Check if the track is within the EMCal acceptance
        if(abs(track->Eta()) > 0.8) continue;
        if(!IsTrackInEmcalAcceptance(track)) continue;

        float TrackEtaOnEMCBefore = track->GetTrackEtaOnEMCal();
        float TrackPhiOnEMCBefore = track->GetTrackPhiOnEMCal();

        // Propagate the track
        AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track, 440, 0.1396, 20, 0.35, kFALSE, kTRUE, kFALSE);

        if (fCreateHisto){
          fCellTrackMatchdEtaDiff->Fill(TrackEtaOnEMCBefore - track->GetTrackEtaOnEMCal());
          fCellTrackMatchdPhiDiff->Fill(TrackPhiOnEMCBefore - track->GetTrackPhiOnEMCal());
        }
      }
      // Here check if there is any match
      Bool_t matched = kFALSE;
      Double_t matchWindow = 0.0132/2.; //half a cell diagonal
      Int_t matchedCellID = -1;
      Int_t CellID = -1;
      Double_t cellEta, cellPhi;

      if(track->GetTrackEtaOnEMCal() == -999 || track->GetTrackPhiOnEMCal() == -999) continue;

      fGeom->GetAbsCellIdFromEtaPhi(track->GetTrackEtaOnEMCal(),track->GetTrackPhiOnEMCal(),CellID);
      fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
      if( abs(cellEta-track->GetTrackEtaOnEMCal())<matchWindow && abs(cellPhi-track->GetTrackPhiOnEMCal())<matchWindow ){
        matchedCellID = CellID;
        matched = kTRUE;
      }
      //the preceding doesn't always match, because it returns an approximate location
      //so here we look around the track to double check the match
      Double_t shiftEta[3] = {-0.0132/2., 0., +0.0132/2.};
      Double_t shiftPhi[3] = {-0.0132/2., 0., +0.0132/2.};
      if(!matched){
        for(Int_t iEta = 0; iEta < 3; iEta++){
          for(Int_t iPhi = 0; iPhi < 3; iPhi++){
            fGeom->GetAbsCellIdFromEtaPhi(track->GetTrackEtaOnEMCal()+shiftEta[iEta],track->GetTrackPhiOnEMCal()+shiftPhi[iPhi],CellID);
            fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
            if( abs(cellEta-track->GetTrackEtaOnEMCal())<matchWindow && abs(cellPhi-track->GetTrackPhiOnEMCal())<matchWindow ){
              matchedCellID = CellID;
              matched = kTRUE;
              break;
            }
          }
          if(matched) break;
        }
      }
      if(matched){
        //here we subtract the energy from the matched cell
        Short_t  iCell = fCaloCells->GetCellPosition(matchedCellID);
        Short_t  absId  =-1;
        Double_t ecell = 0;
        Double_t tcell = 0;
        Double_t efrac = 0;
        Int_t  mclabel = -1;
        Bool_t cellHighGain = fCaloCells->GetHighGain(iCell);
        fCaloCells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);
        fGeom->EtaPhiFromIndex(absId, cellEta, cellPhi);
        if (fCreateHisto) {
          fCellNTrackMatch->Fill(0.,ecell);

          // calculate matching efficiency
          // this is done on cell level. However, if track-hit is part of another cluster and the energy of another particle in this cell is high, it would count as not matched.
          // Information of all particles in the cell is not available (only for the highest). ToDo: Is there a better way?
          if(fMCEvent){
            Int_t MCLabelCell = fCaloCells->GetMCLabel(iCell);
            Int_t MCLabelTrack = track->GetLabel();
            if(MCLabelCell == MCLabelTrack){
              fHistoMatchingEfficiency->Fill(0., track->P());
            } else {
              fHistoMatchingEfficiency->Fill(1., track->P());
            }
          }
        }
        if(ecell>0){
          if (fCreateHisto) {
            fCellTrackMatchdEtadPhi->Fill(cellEta-track->GetTrackEtaOnEMCal(),cellPhi-track->GetTrackPhiOnEMCal());
            fCellTrackMatchEbefore->Fill(ecell);
          }
          //std::cout << "track eta = " << track->Eta() << "\t track phi = " << track->Phi() << std::endl;
          //std::cout << "matchedCellID = " << matchedCellID << std::endl;
          //std::cout << "cell eta = " << cellEta << "\t cell phi = " << cellPhi << std::endl;
          if (!fMCEvent && ecell > fEmipData) { //event is data, subtract fEmipData
            fCaloCells->SetCell(iCell, absId, ecell-fEmipData, tcell, mclabel, efrac, cellHighGain);
            if (fCreateHisto){
              fCellNTrackMatch->Fill(1.,ecell);
              fCellTrackMatchEafter->Fill(ecell-fEmipData);
            }
          } else if (fMCEvent && ecell > fEmipMC) { //event is MC, subtract fEmipMC
            fCaloCells->SetCell(iCell, absId, ecell-fEmipMC, tcell, mclabel, efrac, cellHighGain);
            if (fCreateHisto){
              fCellNTrackMatch->Fill(1.,ecell);
              fCellTrackMatchEafter->Fill(ecell-fEmipMC);
            }
          } else {
            if(fSubtractMultiple){
              double restEnergy = 0.;
              if(fMCEvent){
                restEnergy = fEmipMC - ecell;
              }
              else{
                restEnergy = fEmipData - ecell;
              }
              fCaloCells->SetCell(iCell, absId, 0., tcell, mclabel, efrac, cellHighGain);
              if (fCreateHisto){
                fCellNTrackMatch->Fill(2.,ecell);
                fCellTrackMatchEafter->Fill(0);
              }
              FindAndSubtractNextCell(absId, track, restEnergy);

            }
            else{
              fCaloCells->SetCell(iCell, absId, 0, tcell, mclabel, efrac, cellHighGain);
              if (fCreateHisto){
                fCellNTrackMatch->Fill(2.,ecell);
                fCellTrackMatchEafter->Fill(0);
              }
            }
          }
        }
      }

    }
  }
  return kTRUE;
}

/**
 * Determines if a track is inside the EMCal acceptance, using \f$\eta\f$/\f$\phi\f$ at the vertex (no propagation).
 * Includes +/- edges. Useful to determine whether track propagation should be attempted.
 * @param[in] part Particle to check
 * @param[in] edges Size of the edges in \f$\phi\f$ excluded from the EMCAL acceptance
 * @return True if a particle is inside the EMCAL acceptance, false otherwise
 */
Bool_t AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::IsTrackInEmcalAcceptance(AliVTrack* part, Double_t edges) const
{

  if (!fGeom) {
    AliWarning(Form("%s - AliAnalysisTaskEmcal::IsTrackInEmcalAcceptance - Geometry is not available!", GetName()));
    return kFALSE;
  }

  Double_t minPhi = fGeom->GetArm1PhiMin()*2.*TMath::Pi()/360. - edges;
  Double_t maxPhi = fGeom->GetArm1PhiMax()*2.*TMath::Pi()/360. + edges;

  if (part->Phi() > minPhi && part->Phi() < maxPhi) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}



Bool_t AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction::FindAndSubtractNextCell(Int_t absID, AliVTrack* track, double restEnergy) const
{
  // get the phi direction of the track
  Float_t dPhi = track->GetTrackPhiOnEMCal()-track->Phi();

  // get the original cell that did not have enough energy
  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  fGeom->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

  // Get close cells (in phi only) index, energy and time, not in corners
  Int_t absID1 = -1;
  Int_t absID2 = -1;

  if ( iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if ( iphi > 0 )                                absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  Float_t ecell1 = 0, ecell2 = 0;

  ecell1 = fCaloCells->GetCellAmplitude(absID1);
  ecell2 = fCaloCells->GetCellAmplitude(absID2);

  // if both next cells don't have energy do nothing, because we can not subtract energy from either
  if( ecell1 <= 0 && ecell2 <= 0){
    return false;
  }

  // Check phi direction between original cell and next cells
  Float_t cellEta, cellPhi, cellEta1, cellEta2, cellPhi1, cellPhi2;
  fGeom->EtaPhiFromIndex(absID, cellEta, cellPhi);
  fGeom->EtaPhiFromIndex(absID1, cellEta1, cellPhi1);
  fGeom->EtaPhiFromIndex(absID2, cellEta2, cellPhi2);

  float dPhi1 = cellPhi1 - cellPhi;

  Short_t  iCell = -1;
  Short_t  absId  =-1;
  Double_t ecell = 0;
  Double_t tcell = 0;
  Double_t efrac = 0;
  Int_t  mclabel = -1;
  Bool_t cellHighGain = kFALSE;


  // try to go into the phi direction that the track points to.
  // if cell in that direction has no energy, try in other direction.
  if(signbit(dPhi) == signbit(dPhi1)){
    if(ecell1 <= 0.){
      if(ecell2 > 0.){
        iCell = fCaloCells->GetCellPosition(absID2);
      }
      else{
        return false;
      }
    }
    else{
      iCell = fCaloCells->GetCellPosition(absID1);
    } 
  }
  else{
    if(ecell2 <= 0.){
      if(ecell2 > 0.){
        iCell = fCaloCells->GetCellPosition(absID1);
      }
      else{
        return false;
      }
    }
    else{
      iCell = fCaloCells->GetCellPosition(absID2);
    }
  }

  cellHighGain = fCaloCells->GetHighGain(iCell);
  fCaloCells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);

  if(restEnergy > ecell){
    fCaloCells->SetCell(iCell, absId, 0., tcell, mclabel, efrac, cellHighGain);
    if (fCreateHisto){
      fCellNTrackMatch->Fill(4.,ecell);
      // fCellTrackMatchEafter->Fill(0);
    }
  }
  else{
    fCaloCells->SetCell(iCell, absId, ecell - restEnergy, tcell, mclabel, efrac, cellHighGain);
    if (fCreateHisto){
      fCellNTrackMatch->Fill(3.,ecell);
      // fCellTrackMatchEafter->Fill(ecell - restEnergy);
    }
  }

  return true;
}
