// AliEmcalCorrectionCellTrackMatcher
//

#include "AliEmcalCorrectionCellTrackMatcher.h"

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
ClassImp(AliEmcalCorrectionCellTrackMatcher);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellTrackMatcher> AliEmcalCorrectionCellTrackMatcher::reg("AliEmcalCorrectionCellTrackMatcher");

/**
 * Default constructor
 */
AliEmcalCorrectionCellTrackMatcher::AliEmcalCorrectionCellTrackMatcher() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellTrackMatcher"),
  fMipE(290.),
  fEmcalTracks(0),
  fCellTrackMatchdEtadPhi(0),
  fCellNTrackMatch(0),
  fCellTrackMatchEbefore(0),
  fCellTrackMatchEafter(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellTrackMatcher::~AliEmcalCorrectionCellTrackMatcher()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellTrackMatcher::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  GetProperty("mipE", fMipE);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellTrackMatcher::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  // Determine all particle container array names for naming the AliEmcalParticle array
  std::string particleContainerNames = "";
  bool firstLoop = true;
  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    if (firstLoop != true) {
      particleContainerNames += "_";
    }
    else {
      firstLoop = false;
    }
    particleContainerNames += partCont->GetArrayName();
  }

  fEmcalTracks = new TClonesArray("AliEmcalParticle");
  fEmcalTracks->SetName(Form("EmcalTracks_%s", particleContainerNames.c_str()));

  if (fCreateHisto){
    fCellTrackMatchdEtadPhi = new TH2F("hCellTrackMatchdEtadPhi","hCellTrackMatchdEtadPhi",100,-0.04,0.04,100,-0.04,0.04);
    fOutput->Add(fCellTrackMatchdEtadPhi);
    fCellNTrackMatch = new TH2F("hCellNTrackMatch","hCellNTrackMatch",3,-0.5,2.5,500,0,20);
    fOutput->Add(fCellNTrackMatch);
    fCellTrackMatchEbefore = new TH1F("hCellTrackMatchEbefore","hCellTrackMatchEbefore",500,0,20);
    fOutput->Add(fCellTrackMatchEbefore);
    fCellTrackMatchEafter = new TH1F("hCellTrackMatchEafter","hCellTrackMatchEafter",500,0,20);
    fOutput->Add(fCellTrackMatchEafter);
  }


}

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionCellTrackMatcher::ExecOnce()
{
  fParticleContainerIndexMap.CopyMappingFrom(AliParticleContainer::GetEmcalContainerIndexMap(), fParticleCollArray);
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellTrackMatcher::Run()
{
  AliEmcalCorrectionComponent::Run();

  fEmcalTracks->Delete();
  AliVTrack* track = 0;

  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    auto partItCont = partCont->accepted_momentum();
    for (AliParticleIterableMomentumContainer::iterator partIterator = partItCont.begin(); partIterator != partItCont.end(); ++partIterator) {
      track = static_cast<AliVTrack *>(partIterator->second);
      // Propagate the track
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track, 440, 0.1396, 20, 0.35, kFALSE, kTRUE, kFALSE);
      //std::cout << "found a track with pt = " << track->Pt() << "\t eta = " << track->Eta() << "\t phi = " << track->Phi() << std::endl;
      // Check if the track is within the EMCal acceptance
      if(abs(track->Eta()) > 0.8) continue;
      if(!IsTrackInEmcalAcceptance(track)) continue;

      //std::cout << "track eta = " << track->Eta() << "\t track phi = " << track->Phi() << std::endl;

      // Here check if there is any match
      Bool_t matched = kFALSE;
      Double_t matchWindow = 0.0132/2.; //half a cell diagonal
      Int_t matchedCellID = -1;
      Int_t CellID = -1;
      Double_t cellEta, cellPhi;

      fGeom->GetAbsCellIdFromEtaPhi(track->Eta(),track->Phi(),CellID);
      fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
      if( abs(cellEta-track->Eta())<matchWindow && abs(cellPhi-track->Phi())<matchWindow ){
        matchedCellID = CellID;
        matched = kTRUE;
      }
      //this doesn't always already match the correct cell, because it returns an approximate location
      //so here we look around the track for the real match
      Double_t shiftEta;
      Double_t shiftPhi;
      if(!matched){
        shiftEta = -0.0132/2.;
        shiftPhi = -0.0132/2.;
        for(Int_t i=0; i<3 ; i++){
          fGeom->GetAbsCellIdFromEtaPhi(track->Eta()+shiftEta,track->Phi()+shiftPhi,CellID);
          fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
          if( abs(cellEta-track->Eta())<matchWindow && abs(cellPhi-track->Phi())<matchWindow ){
            matchedCellID = CellID;
            matched = kTRUE;
            break;
          }
          shiftPhi += 0.0132/2.;
        }
      }
      if(!matched){
        shiftEta = 0;
        shiftPhi = -0.0132/2.;
        for(Int_t i=0; i<3 ; i++){
          fGeom->GetAbsCellIdFromEtaPhi(track->Eta()+shiftEta,track->Phi()+shiftPhi,CellID);
          fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
          if( abs(cellEta-track->Eta())<matchWindow && abs(cellPhi-track->Phi())<matchWindow ){
            matchedCellID = CellID;
            matched = kTRUE;
            break;
          }
          shiftPhi += 0.0132/2.;
        }
      }
      if(!matched){
        shiftEta = 0.0132/2.;
        shiftPhi = -0.0132/2.;
        for(Int_t i=0; i<3 ; i++){
          fGeom->GetAbsCellIdFromEtaPhi(track->Eta()+shiftEta,track->Phi()+shiftPhi,CellID);
          fGeom->EtaPhiFromIndex(CellID, cellEta, cellPhi);
          if( abs(cellEta-track->Eta())<matchWindow && abs(cellPhi-track->Phi())<matchWindow ){
            matchedCellID = CellID;
            matched = kTRUE;
            break;
          }
          shiftPhi += 0.0132/2.;
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
        fCellNTrackMatch->Fill(0.,ecell);
        if(ecell>0){
          fCellTrackMatchdEtadPhi->Fill(cellEta-track->Eta(),cellPhi-track->Phi());
          //std::cout << "track eta = " << track->Eta() << "\t track phi = " << track->Phi() << std::endl;
          //std::cout << "matchedCellID = " << matchedCellID << std::endl;
          //std::cout << "cell eta = " << cellEta << "\t cell phi = " << cellPhi << std::endl;
          fCellTrackMatchEbefore->Fill(ecell);
          if (ecell > fMipE) {
            fCaloCells->SetCell(iCell, absId, ecell-fMipE, tcell, mclabel, efrac, cellHighGain);
            fCellNTrackMatch->Fill(1.,ecell);
            fCellTrackMatchEafter->Fill(ecell-fMipE);
          } else {
            fCaloCells->SetCell(iCell, absId, 0, tcell, mclabel, efrac, cellHighGain);
            fCellNTrackMatch->Fill(2.,ecell);
            fCellTrackMatchEafter->Fill(0);
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
Bool_t AliEmcalCorrectionCellTrackMatcher::IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges) const
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
