/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- ROOT system ---
#include <TClass.h>
#include <TMath.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>

//---- ANALYSIS system ----
#include "AliNeutralMesonSelection.h"
#include "AliAnaParticleHadronCorrelation.h"
#include "AliCaloTrackReader.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliFiducialCut.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliMCAnalysisUtils.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEventplane.h"
#include "AliIsolationCut.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

// --- Detectors ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

/// \cond CLASSIMP
ClassImp(AliAnaParticleHadronCorrelation) ;
/// \endcond


//___________________________________________________________________
/// Default Constructor. Initialize parameters.
//___________________________________________________________________
AliAnaParticleHadronCorrelation::AliAnaParticleHadronCorrelation():
AliAnaCaloTrackCorrBaseClass(),
fFillAODWithReferences(0),      fCheckLeadingWithNeutralClusters(0),
fMaxAssocPt(1000.),             fMinAssocPt(0.),
fDeltaPhiMaxCut(0.),            fDeltaPhiMinCut(0.),
fSelectIsolated(0),             fMakeSeveralUE(0),
fUeDeltaPhiMaxCut(0.),          fUeDeltaPhiMinCut(0.),
fPi0AODBranchName(""),          fAODNamepTInConeHisto(""), fNeutralCorr(0),
fPi0Trigger(0),                 fDecayTrigger(0),
fNDecayBits(0),                 fDecayBits(),
fNBkgBin(0),
fMakeAbsoluteLeading(0),        fMakeNearSideLeading(0),
fLeadingTriggerIndex(-1),       fHMPIDCorrelation(0),  fFillBradHisto(0),
fNAssocPtBins(0),               fAssocPtBinLimit(),
fCorrelVzBin(0),
fListMixTrackEvents(),          fListMixCaloEvents(),
fUseMixStoredInReader(0),       fFillNeutralEventMixPool(0),
fM02MaxCut(0),                  fM02MinCut(0),
fSelectLeadingHadronAngle(0),   fFillLeadHadOppositeHisto(0),
fMinLeadHadPhi(0),              fMaxLeadHadPhi(0),
fMinLeadHadPt(0),               fMaxLeadHadPt(0),
fFillEtaGapsHisto(1),           fFillMomImbalancePtAssocBinsHisto(0),
fFillInvMassHisto(0),           fFillBkgBinsHisto(0),
fFillTaggedDecayHistograms(0),  fDecayTagsM02Cut(0),
fMCGenTypeMin(0),               fMCGenTypeMax(0),
fTrackVector(),                 fMomentum(),           fMomentumIM(),
fDecayMom1(),                   fDecayMom2(),
//Histograms
fhPtTriggerInput(0),            fhPtTriggerSSCut(0),
fhPtTriggerIsoCut(0),           fhPtTriggerFidCut(0),
fhPtTrigger(0),                 fhPtTriggerVtxBC0(0),
fhPtTriggerVzBin(0),            fhPtTriggerBin(0),
fhPhiTrigger(0),                fhEtaTrigger(0),
fhPtTriggerMC(),
fhPtDecayTrigger(),             fhPtDecayTriggerMC(),
fhPtTriggerCentrality(0),       fhPtTriggerEventPlane(0),
fhTriggerEventPlaneCentrality(0),
fhPtTriggerMixed(0),            fhPtTriggerMixedVzBin(0), fhPtTriggerMixedBin(0),
fhPhiTriggerMixed(0),           fhEtaTriggerMixed(0),
fhPtLeadingOppositeHadron(0),   fhPtDiffPhiLeadingOppositeHadron(0), fhPtDiffEtaLeadingOppositeHadron(0),
fhPtNoLeadingOppositeHadron(0), fhEtaPhiNoLeadingOppositeHadron(0),
fhDeltaPhiDeltaEtaCharged(0),
fhPhiCharged(0),                fhEtaCharged(0),
fhDeltaPhiCharged(0),           fhDeltaEtaCharged(0),
fhDeltaPhiChargedPt(0),         fhDeltaPhiUeChargedPt(0),
fhUePart(0),
fhXECharged(0),                 fhXECharged_Cone2(0),      fhXEUeCharged(0),
fhXEUeChargedSmallCone(0),      fhXEUeChargedMediumCone(0), fhXEUeChargedLargeCone(0),
fhXEPosCharged(0),              fhXENegCharged(0),
fhPtHbpXECharged(0),            fhPtHbpXECharged_Cone2(0), fhPtHbpXEUeCharged(0),
fhZTCharged(0),                 fhZTUeCharged(0),
fhZTPosCharged(0),              fhZTNegCharged(0),
fhPtHbpZTCharged(0),            fhPtHbpZTUeCharged(0),
fhXEChargedMC(),                fhDeltaPhiChargedMC(),
fhXEUeChargedRightMC(),         fhXEUeChargedLeftMC(),
fhDeltaPhiDeltaEtaChargedPtA3GeV(0),
fhDeltaPhiChargedPtA3GeV(0),    fhDeltaEtaChargedPtA3GeV(0),
// Pile-Up
fhDeltaPhiChargedPileUp(),      fhDeltaEtaChargedPileUp(),
fhDeltaPhiChargedPtA3GeVPileUp(), fhDeltaEtaChargedPtA3GeVPileUp(),
fhXEChargedPileUp(),            fhXEUeChargedPileUp(),
fhZTChargedPileUp(),            fhZTUeChargedPileUp(),
fhPtTrigChargedPileUp(),
fhDeltaPhiChargedOtherBC(),     fhDeltaPhiChargedPtA3GeVOtherBC(),
fhXEChargedOtherBC(),           fhXEUeChargedOtherBC(),
fhZTChargedOtherBC(),           fhZTUeChargedOtherBC(),
fhPtTrigChargedOtherBC(),
fhDeltaPhiChargedBC0(),         fhDeltaPhiChargedPtA3GeVBC0(),
fhXEChargedBC0(),               fhXEUeChargedBC0(),
fhZTChargedBC0(),               fhZTUeChargedBC0(),
fhPtTrigChargedBC0(),
fhDeltaPhiChargedVtxBC0(),      fhDeltaPhiChargedPtA3GeVVtxBC0(),
fhXEChargedVtxBC0(),            fhXEUeChargedVtxBC0(),
fhZTChargedVtxBC0(),            fhZTUeChargedVtxBC0(),
fhPtTrigChargedVtxBC0(),
fhDeltaPhiUeLeftCharged(0),
fhDeltaPhiUeLeftUpCharged(0),   fhDeltaPhiUeRightUpCharged(0),
fhDeltaPhiUeLeftDownCharged(0), fhDeltaPhiUeRightDownCharged(0),
fhXEUeLeftCharged(0),
fhXEUeLeftUpCharged(0),         fhXEUeRightUpCharged(0),
fhXEUeLeftDownCharged(0),       fhXEUeRightDownCharged(0),
fhPtHbpXEUeLeftCharged(0),      fhZTUeLeftCharged(0),
fhPtHbpZTUeLeftCharged(0),
fhPtTrigPout(0),                fhPtTrigCharged(0),
fhDeltaPhiChargedMult(0x0),     fhDeltaEtaChargedMult(0x0),
fhXEMult(0x0),                  fhXEUeMult(0x0),
fhZTMult(0x0),                  fhZTUeMult(0x0),
fhAssocPtBkg(0),                fhDeltaPhiDeltaEtaAssocPtBin(0),
fhDeltaPhiAssocPtBin(0),
fhDeltaPhiAssocPtBinDEta08(0),  fhDeltaPhiAssocPtBinDEta0(0),
fhDeltaPhiAssocPtBinHMPID(0),   fhDeltaPhiAssocPtBinHMPIDAcc(0),
fhDeltaPhiBradAssocPtBin(0),    fhDeltaPhiBrad(0),
fhXEAssocPtBin(0),              fhZTAssocPtBin(0),
fhXEVZ(0),                      fhZTVZ(0),
fhDeltaPhiDeltaEtaNeutral(0),
fhPhiNeutral(0),                fhEtaNeutral(0),
fhDeltaPhiNeutral(0),           fhDeltaEtaNeutral(0),
fhDeltaPhiNeutralPt(0),         fhDeltaPhiUeNeutralPt(0),
fhXENeutral(0),                 fhXEUeNeutral(0),
fhPtHbpXENeutral(0),            fhPtHbpXEUeNeutral(0),
fhZTNeutral(0),                 fhZTUeNeutral(0),
fhPtHbpZTNeutral(0),            fhPtHbpZTUeNeutral(0),
fhDeltaPhiUeLeftNeutral(0),     fhXEUeLeftNeutral(0),
fhPtHbpXEUeLeftNeutral(0),      fhZTUeLeftNeutral(0),
fhPtHbpZTUeLeftNeutral(0),      fhPtPi0DecayRatio(0),
fhDeltaPhiPi0DecayCharged(0),   fhXEPi0DecayCharged(0),        fhZTPi0DecayCharged(0),
fhDeltaPhiPi0DecayNeutral(0),   fhXEPi0DecayNeutral(0),        fhZTPi0DecayNeutral(0),
fhDeltaPhiDecayCharged(),       fhXEDecayCharged(),            fhZTDecayCharged(),
fhDeltaPhiDecayChargedAssocPtBin(),
fhMCPtTrigger(),                fhMCPhiTrigger(),              fhMCEtaTrigger(),
fhMCPtTriggerNotLeading(),      fhMCPhiTriggerNotLeading(),    fhMCEtaTriggerNotLeading(),
fhMCEtaCharged(),               fhMCPhiCharged(),
fhMCDeltaEtaCharged(),          fhMCDeltaPhiCharged(),
fhMCDeltaPhiDeltaEtaCharged(),  fhMCDeltaPhiChargedPt(),
fhMCPtXECharged(),              fhMCPtXEUeCharged(),
fhMCPtXEUeLeftCharged(),
fhMCPtHbpXECharged(),           fhMCPtHbpXEUeCharged(),
fhMCPtHbpXEUeLeftCharged(),
fhMCUePart(),
fhMCPtZTCharged(),              fhMCPtZTUeCharged(),
fhMCPtZTUeLeftCharged(),
fhMCPtHbpZTCharged(),           fhMCPtHbpZTUeCharged(),
fhMCPtHbpZTUeLeftCharged(),
fhMCPtTrigPout(),               fhMCPtAssocDeltaPhi(),
// Mixing
fhNEventsTrigger(0),            fhNtracksMB(0),                 fhNclustersMB(0),
fhMixDeltaPhiCharged(0),        fhMixDeltaPhiDeltaEtaCharged(0),
fhMixXECharged(0),              fhMixXEUeCharged(0),            fhMixHbpXECharged(0),
fhMixDeltaPhiChargedAssocPtBin(),
fhMixDeltaPhiChargedAssocPtBinDEta08(),
fhMixDeltaPhiChargedAssocPtBinDEta0(),
fhMixDeltaPhiDeltaEtaChargedAssocPtBin(),
fhEventBin(0),                  fhEventMixBin(0),               fhEventMBBin(0),
fhMassPtTrigger(0),             fhMCMassPtTrigger(),
fhPtLeadInConeBin(),            fhPtSumInConeBin(),
fhPtLeadConeBinDecay(),         fhSumPtConeBinDecay(),
fhPtLeadConeBinMC(),            fhSumPtConeBinMC(),
fhTrackResolution(0),           fhTrackResolutionUE(0)
{
  InitParameters();
  
  for(Int_t i = 0; i < fgkNmcTypes; i++)
  {
    fhPtTriggerMC[i] = 0;
    fhXEChargedMC[i] = 0;
    fhDeltaPhiChargedMC[i] = 0;
    fhMCMassPtTrigger  [i] = 0 ;
    fhXEUeChargedRightMC[i] = 0;
    fhXEUeChargedLeftMC[i]  = 0;
    
    for(Int_t ib = 0; ib < AliNeutralMesonSelection::fgkMaxNDecayBits; ib++)  fhPtDecayTriggerMC[ib][i] = 0;
  }
  
  for(Int_t ib = 0; ib < AliNeutralMesonSelection::fgkMaxNDecayBits; ib++)  fhPtDecayTrigger[ib] = 0;
  
  for(Int_t i = 0; i < 7; i++)
  {
    fhPtTriggerPileUp             [i] = 0 ;
    fhDeltaPhiChargedPileUp       [i] = 0 ; fhDeltaEtaChargedPileUp       [i] = 0 ;
    fhXEChargedPileUp             [i] = 0 ; fhXEUeChargedPileUp           [i] = 0 ;
    fhZTChargedPileUp             [i] = 0 ; fhZTUeChargedPileUp           [i] = 0 ;
    fhPtTrigChargedPileUp         [i] = 0 ;
    fhDeltaPhiChargedPtA3GeVPileUp[i] = 0 ; fhDeltaEtaChargedPtA3GeVPileUp[i] = 0 ;
  }
}

//_________________________________________________________________
/// Destructor. Remove event containers.
//_________________________________________________________________
AliAnaParticleHadronCorrelation::~AliAnaParticleHadronCorrelation()
{
  if(DoOwnMix())
  {
    if(fListMixTrackEvents)
    {
      for(Int_t iz=0; iz < GetNZvertBin(); iz++)
      {
        for(Int_t ic=0; ic < GetNCentrBin(); ic++)
        {
          for(Int_t irp=0; irp<GetNRPBin(); irp++)
          {
            Int_t bin = GetEventMixBin(ic, iz, irp);
            fListMixTrackEvents[bin]->Delete() ;
            delete fListMixTrackEvents[bin] ;
          }
        }
      }
    }
    
    delete[] fListMixTrackEvents;
    
    if(fListMixCaloEvents)
    {
      for(Int_t iz=0; iz < GetNZvertBin(); iz++)
      {
        for(Int_t ic=0; ic < GetNCentrBin(); ic++)
        {
          for(Int_t irp=0; irp<GetNRPBin(); irp++)
          {
            Int_t bin = GetEventMixBin(ic, iz, irp);
            fListMixCaloEvents[bin]->Delete() ;
            delete fListMixCaloEvents[bin] ;
          }
        }
      }
    }
    
    delete[] fListMixCaloEvents;
  }
}

//____________________________________________________________________________________________________________________________________
/// Fill angular correlation related histograms.
//____________________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedAngularCorrelationHistograms(Float_t ptAssoc,  Float_t ptTrig,      Int_t   bin,
                                                                              Float_t phiAssoc, Float_t phiTrig,     Float_t deltaPhi,
                                                                              Float_t etaAssoc, Float_t etaTrig,
                                                                              Int_t   decayTag, Float_t hmpidSignal, Int_t  outTOF,
                                                                              Int_t   cen,      Int_t   mcTag)
{
  Float_t deltaEta    = etaTrig-etaAssoc;
  Float_t deltaPhiOrg = phiTrig-phiAssoc;
  
  fhEtaCharged       ->Fill(ptAssoc, etaAssoc, GetEventWeight());
  fhPhiCharged       ->Fill(ptAssoc, phiAssoc, GetEventWeight());
  fhDeltaEtaCharged  ->Fill(ptTrig , deltaEta, GetEventWeight());
  fhDeltaPhiCharged  ->Fill(ptTrig , deltaPhi, GetEventWeight());
  fhDeltaPhiChargedPt->Fill(ptAssoc, deltaPhi, GetEventWeight());
  fhDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta, GetEventWeight());
  
  if(ptAssoc > 3 )
  {
    fhDeltaEtaChargedPtA3GeV        ->Fill(ptTrig  , deltaEta, GetEventWeight());
    fhDeltaPhiChargedPtA3GeV        ->Fill(ptTrig  , deltaPhi, GetEventWeight());
    fhDeltaPhiDeltaEtaChargedPtA3GeV->Fill(deltaPhi, deltaEta, GetEventWeight());
  }
  
  // Pile up studies
  
  if(IsPileUpAnalysisOn())
  {
    if     (outTOF==1)
    {
      fhDeltaPhiChargedOtherBC->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVOtherBC->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    else if(outTOF==0)
    {
      fhDeltaPhiChargedBC0->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVBC0->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    
    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhDeltaPhiChargedVtxBC0->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      if(ptAssoc > 3 ) fhDeltaPhiChargedPtA3GeVVtxBC0->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    
    if(GetReader()->IsPileUpFromSPD())
    {
      fhDeltaEtaChargedPileUp[0]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[0]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromEMCal())
    {
      fhDeltaEtaChargedPileUp[1]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[1]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromSPDOrEMCal())
    {
      fhDeltaEtaChargedPileUp[2]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[2]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromSPDAndEMCal())
    {
      fhDeltaEtaChargedPileUp[3]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[3]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())
    {
      fhDeltaEtaChargedPileUp[4]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[4]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())
    {
      fhDeltaEtaChargedPileUp[5]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[5]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
    {
      fhDeltaEtaChargedPileUp[6]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
      fhDeltaPhiChargedPileUp[6]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
    }
    
    if(ptAssoc > 3 )
    {
      if(GetReader()->IsPileUpFromSPD())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[0]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[0]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromEMCal())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[1]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[1]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromSPDOrEMCal())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[2]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[2]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromSPDAndEMCal())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[3]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[3]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[4]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[4]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[5]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[5]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
      {
        fhDeltaEtaChargedPtA3GeVPileUp[6]->Fill(ptTrig, deltaEta, GetEventWeight()) ;
        fhDeltaPhiChargedPtA3GeVPileUp[6]->Fill(ptTrig, deltaPhi, GetEventWeight()) ;
      }
    }
  }
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
    fhDeltaPhiChargedMC[mcIndex]->Fill(ptTrig , deltaPhi, GetEventWeight());
    
    if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
    {
      // check index in GetMCTagIndexHistogram
      if     ( mcIndex == 2 ) fhDeltaPhiChargedMC[8]->Fill(ptTrig, deltaPhi, GetEventWeight()); // pi0 decay
      else if( mcIndex == 4 ) fhDeltaPhiChargedMC[9]->Fill(ptTrig, deltaPhi, GetEventWeight()); // eta decay
    }
  }
  
  if(fDecayTrigger && decayTag > 0)
  {
    for(Int_t ibit = 0; ibit<fNDecayBits; ibit++)
    {
      if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
        fhDeltaPhiDecayCharged[ibit]->Fill(ptTrig, deltaPhi, GetEventWeight());
    }
  }
  
  Double_t  dphiBrad = -100;
  if(fFillBradHisto)
  {
    dphiBrad = atan2(sin(deltaPhiOrg), cos(deltaPhiOrg))/TMath::Pi();//-1 to 1
    if( TMath::Abs(dphiBrad) > 0.325 && TMath::Abs(dphiBrad) < 0.475 )  // Hardcoded values, BAD, FIXME
    {
      fhAssocPtBkg->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    
    if( dphiBrad < -1./3 ) dphiBrad += 2;
    fhDeltaPhiBrad->Fill(ptTrig, dphiBrad, GetEventWeight());
  }
  
  // Fill histograms in bins of associated particle pT
  if ( bin >= 0 )
  {
    fhDeltaPhiDeltaEtaAssocPtBin[bin]->Fill(deltaPhi,deltaEta, GetEventWeight());
    
    fhDeltaPhiAssocPtBin        [bin]->Fill(ptTrig,  deltaPhi, GetEventWeight());
    
    if(fFillEtaGapsHisto)
    {
      if(TMath::Abs(deltaEta)> 0.8)
        fhDeltaPhiAssocPtBinDEta08[bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
      
      if(TMath::Abs(deltaEta)< 0.01)
        fhDeltaPhiAssocPtBinDEta0 [bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
    }
    
    if (fFillBradHisto)
      fhDeltaPhiBradAssocPtBin    [bin]->Fill(ptTrig, dphiBrad, GetEventWeight());
    
    if(fDecayTrigger && decayTag > 0 && fNDecayBits > 0 &&
       GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
    {
      fhDeltaPhiDecayChargedAssocPtBin[bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
    }
    
    if(fHMPIDCorrelation)
    {
      if( hmpidSignal > 0 )
      {
        //printf("Track pt %f with HMPID signal %f \n",pt,hmpidSignal);
        fhDeltaPhiAssocPtBinHMPID[bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
      }
      
      if(phiAssoc > 5*TMath::DegToRad() && phiAssoc < 20*TMath::DegToRad())
      {
        //printf("Track pt %f in HMPID acceptance phi %f \n ",pt,phi*TMath::RadToDeg() );
        fhDeltaPhiAssocPtBinHMPIDAcc[bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
      }
    }
  }
  
  // Fill different multiplicity/centrality histogram
  if(IsHighMultiplicityAnalysisOn() && cen >= 0 && cen < GetNCentrBin())
  {
    fhDeltaPhiChargedMult[cen]->Fill(ptTrig, deltaPhi, GetEventWeight());
    fhDeltaEtaChargedMult[cen]->Fill(ptTrig, deltaEta, GetEventWeight());
  }
}

//___________________________________________________________________________________________________________________________________
/// Fill MC histograms independently of AOD or ESD.
//___________________________________________________________________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::FillChargedMCCorrelationHistograms(Float_t mcAssocPt, Float_t mcAssocPhi, Float_t mcAssocEta,
                                                                           Float_t mcTrigPt,  Float_t mcTrigPhi,  Float_t mcTrigEta,
                                                                           Int_t histoIndex, Bool_t lostDecayPair)
{
  Bool_t lead = kTRUE;
  
  // In case we requested the trigger to be a leading particle,
  // check if this is true at the MC level.
  // Not sure if it is correct to skip or not skip this.
  // Absolute leading?
  if( fMakeAbsoluteLeading && mcAssocPt > mcTrigPt )     lead = kFALSE; // skip event
  
  // Skip this event if near side associated particle pt larger than trigger
  if( mcAssocPhi < 0 ) mcAssocPhi+=TMath::TwoPi();
  
  Float_t mcdeltaPhi= mcTrigPhi-mcAssocPhi;
  if(mcdeltaPhi <= -TMath::PiOver2()) mcdeltaPhi+=TMath::TwoPi();
  if(mcdeltaPhi > 3*TMath::PiOver2()) mcdeltaPhi-=TMath::TwoPi();
  
  if( fMakeNearSideLeading)
  {
    if( mcAssocPt > mcTrigPt && mcdeltaPhi < TMath::PiOver2() ) lead = kFALSE; // skip event
  }
  
  //
  // Select only hadrons in pt range
  if ( mcAssocPt < fMinAssocPt || mcAssocPt > fMaxAssocPt ) return lead ; // exclude but continue
  if ( mcAssocPt < GetReader()->GetCTSPtMin())              return lead ;
  
  
  //
  // Remove trigger itself for correlation when use charged triggers
  if(TMath::Abs(mcAssocPt -mcTrigPt ) < 1e-6 &&
     mcdeltaPhi                       < 1e-6 &&
     TMath::Abs(mcAssocEta-mcTrigEta) < 1e-6)            return lead ; // exclude but continue
  
  Float_t mcxE    =-mcAssocPt/mcTrigPt*TMath::Cos(mcdeltaPhi);// -(mcAssocPx*pxprim+mcAssocPy*pyprim)/(mcTrigPt*mcTrigPt);
  Float_t mchbpXE =-100 ;
  if(mcxE > 0 ) mchbpXE = TMath::Log(1./mcxE);
  
  Float_t mczT    = mcAssocPt/mcTrigPt ;
  Float_t mchbpZT =-100 ;
  if(mczT > 0 ) mchbpZT = TMath::Log(1./mczT);
  
  Double_t mcpout = mcAssocPt*TMath::Sin(mcdeltaPhi) ;
  
  AliDebug(1,Form("Charged hadron: track Pt %f, track Phi %f, phi trigger %f. Cuts:  delta phi  %2.2f < %2.2f < %2.2f",
                  mcAssocPt,mcAssocPhi, mcTrigPhi,fDeltaPhiMinCut, mcdeltaPhi, fDeltaPhiMaxCut));
  
  // Fill Histograms
  fhMCEtaCharged     [histoIndex]->Fill(mcAssocPt, mcAssocEta, GetEventWeight());
  fhMCPhiCharged     [histoIndex]->Fill(mcAssocPt, mcAssocPhi, GetEventWeight());
  fhMCDeltaPhiCharged[histoIndex]->Fill(mcTrigPt , mcdeltaPhi, GetEventWeight());
  fhMCPtAssocDeltaPhi[histoIndex]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
  
  fhMCDeltaEtaCharged        [histoIndex]->Fill(mcTrigPt  , mcTrigEta-mcAssocEta, GetEventWeight());
  fhMCDeltaPhiDeltaEtaCharged[histoIndex]->Fill(mcdeltaPhi, mcTrigEta-mcAssocEta, GetEventWeight());
  
  // Delta phi cut for correlation
  if( (mcdeltaPhi > fDeltaPhiMinCut) && (mcdeltaPhi < fDeltaPhiMaxCut) )
  {
    fhMCDeltaPhiChargedPt[histoIndex]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
    fhMCPtXECharged      [histoIndex]->Fill(mcTrigPt , mcxE      , GetEventWeight());
    fhMCPtHbpXECharged   [histoIndex]->Fill(mcTrigPt , mchbpXE   , GetEventWeight());
    fhMCPtZTCharged      [histoIndex]->Fill(mcTrigPt , mczT      , GetEventWeight());
    fhMCPtHbpZTCharged   [histoIndex]->Fill(mcTrigPt , mchbpZT   , GetEventWeight());
    fhMCPtTrigPout       [histoIndex]->Fill(mcTrigPt , mcpout    , GetEventWeight());
  }
  
  if(lostDecayPair)
  {
    // check index in GetMCTagIndexHistogram
    if(histoIndex == 2 && 8 >= fMCGenTypeMin && 8 <= fMCGenTypeMax )
    {
      // pi0 decay
      fhMCEtaCharged     [8]->Fill(mcAssocPt, mcAssocEta, GetEventWeight());
      fhMCPhiCharged     [8]->Fill(mcAssocPt, mcAssocPhi, GetEventWeight());
      fhMCDeltaPhiCharged[8]->Fill(mcTrigPt , mcdeltaPhi, GetEventWeight());
      fhMCPtAssocDeltaPhi[8]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
      
      fhMCDeltaEtaCharged        [8]->Fill(mcTrigPt  , mcTrigEta-mcAssocEta, GetEventWeight());
      fhMCDeltaPhiDeltaEtaCharged[8]->Fill(mcdeltaPhi, mcTrigEta-mcAssocEta, GetEventWeight());
      
      //delta phi cut for correlation
      if( (mcdeltaPhi > fDeltaPhiMinCut) && (mcdeltaPhi < fDeltaPhiMaxCut) )
      {
        fhMCDeltaPhiChargedPt[8]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
        fhMCPtXECharged      [8]->Fill(mcTrigPt , mcxE      , GetEventWeight());
        fhMCPtHbpXECharged   [8]->Fill(mcTrigPt , mchbpXE   , GetEventWeight());
        fhMCPtZTCharged      [8]->Fill(mcTrigPt , mczT      , GetEventWeight());
        fhMCPtHbpZTCharged   [8]->Fill(mcTrigPt , mchbpZT   , GetEventWeight());
        fhMCPtTrigPout       [8]->Fill(mcTrigPt , mcpout    , GetEventWeight());
      }
    } // pi0 decay lost pair
    
    if(histoIndex == 4 && 9 >= fMCGenTypeMin && 9 <= fMCGenTypeMax )
    {
      // eta decay
      fhMCEtaCharged     [9]->Fill(mcAssocPt, mcAssocEta, GetEventWeight());
      fhMCPhiCharged     [9]->Fill(mcAssocPt, mcAssocPhi, GetEventWeight());
      fhMCDeltaPhiCharged[9]->Fill(mcTrigPt , mcdeltaPhi, GetEventWeight());
      fhMCPtAssocDeltaPhi[9]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
      
      fhMCDeltaEtaCharged        [9]->Fill(mcTrigPt  , mcTrigEta-mcAssocEta, GetEventWeight());
      fhMCDeltaPhiDeltaEtaCharged[9]->Fill(mcdeltaPhi, mcTrigEta-mcAssocEta, GetEventWeight());
      
      // Delta phi cut for correlation
      if( (mcdeltaPhi > fDeltaPhiMinCut) && (mcdeltaPhi < fDeltaPhiMaxCut) )
      {
        fhMCDeltaPhiChargedPt[9]->Fill(mcAssocPt, mcdeltaPhi, GetEventWeight());
        fhMCPtXECharged      [9]->Fill(mcTrigPt , mcxE      , GetEventWeight());
        fhMCPtHbpXECharged   [9]->Fill(mcTrigPt , mchbpXE   , GetEventWeight());
        fhMCPtZTCharged      [9]->Fill(mcTrigPt , mczT      , GetEventWeight());
        fhMCPtHbpZTCharged   [9]->Fill(mcTrigPt , mchbpZT   , GetEventWeight());
        fhMCPtTrigPout       [9]->Fill(mcTrigPt , mcpout    , GetEventWeight()) ;
      }
    } // eta decay lost pair
  }
  // Underlying event
  
  // Right
  if ( (mcdeltaPhi > fUeDeltaPhiMinCut) && (mcdeltaPhi < fUeDeltaPhiMaxCut) )
  {
    //Double_t randomphi = gRandom->Uniform(TMath::Pi()/2,3*TMath::Pi()/2);
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t mcUexE = -(mcAssocPt/mcTrigPt)*TMath::Cos(randomphi);
    Double_t mcUezT =   mcAssocPt/mcTrigPt;
    
    if(mcUexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for right UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      mcUexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhMCPtXEUeCharged[histoIndex]->Fill(mcTrigPt, mcUexE, GetEventWeight());
    if(mcUexE > 0) fhMCPtHbpXEUeCharged[histoIndex]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
    
    fhMCPtZTUeCharged[histoIndex]->Fill(mcTrigPt, mcUezT, GetEventWeight());
    if(mcUezT > 0) fhMCPtHbpZTUeCharged[histoIndex]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
    
    fhMCUePart[histoIndex]->Fill(mcTrigPt, GetEventWeight());
    
    if(lostDecayPair)
    {
      // check index in GetMCTagIndexHistogram
      if(histoIndex == 2 && 8 >= fMCGenTypeMin && 8 <= fMCGenTypeMax )
      {
        // pi0 decay
        fhMCPtXEUeCharged[8]->Fill(mcTrigPt, mcUexE, GetEventWeight());
        if(mcUexE > 0) fhMCPtHbpXEUeCharged[8]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
        
        fhMCPtZTUeCharged[8]->Fill(mcTrigPt, mcUezT, GetEventWeight());
        if(mcUezT > 0) fhMCPtHbpZTUeCharged[8]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
        
        fhMCUePart[8]->Fill(mcTrigPt, GetEventWeight());
      }
      else if(histoIndex == 4 && 9 >= fMCGenTypeMin && 9 <= fMCGenTypeMax )
      {
        // eta decay
        fhMCPtXEUeCharged[9]->Fill(mcTrigPt, mcUexE, GetEventWeight());
        if(mcUexE > 0) fhMCPtHbpXEUeCharged[9]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
        
        fhMCPtZTUeCharged[9]->Fill(mcTrigPt, mcUezT, GetEventWeight());
        if(mcUezT > 0) fhMCPtHbpZTUeCharged[9]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
        
        fhMCUePart[9]->Fill(mcTrigPt, GetEventWeight());
      }
    }
  }
  
  if(fMakeSeveralUE)
  {
    // Left
    if((mcdeltaPhi<-fUeDeltaPhiMinCut) || (mcdeltaPhi >2*fUeDeltaPhiMaxCut))
    {
      Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
      Double_t mcUexE = -(mcAssocPt/mcTrigPt)*TMath::Cos(randomphi);
      Double_t mcUezT =   mcAssocPt/mcTrigPt;
      
      if(mcUexE < 0.)
        AliWarning(Form("Careful!!, negative xE %2.2f for left UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                        mcUexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
      
      fhMCPtXEUeLeftCharged[histoIndex]->Fill(mcTrigPt, mcUexE, GetEventWeight());
      if(mcUexE > 0) fhMCPtHbpXEUeLeftCharged[histoIndex]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
      
      fhMCPtZTUeLeftCharged[histoIndex]->Fill(mcTrigPt, mcUezT, GetEventWeight());
      if(mcUezT > 0) fhMCPtHbpZTUeLeftCharged[histoIndex]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
      
      if(lostDecayPair)
      {
        // Check index in GetMCTagIndexHistogram
        if(histoIndex == 2  && 8 >= fMCGenTypeMin && 8 <= fMCGenTypeMax )
        {
          // pi0 decay
          fhMCPtXEUeLeftCharged[8]->Fill(mcTrigPt,mcUexE, GetEventWeight());
          if(mcUexE > 0) fhMCPtHbpXEUeLeftCharged[8]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
          
          fhMCPtZTUeLeftCharged[8]->Fill(mcTrigPt,mcUezT, GetEventWeight());
          if(mcUezT > 0) fhMCPtHbpZTUeLeftCharged[8]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
        }
        else if(histoIndex == 4  && 9 >= fMCGenTypeMin && 9 <= fMCGenTypeMax )
        {
          // eta decay
          fhMCPtXEUeLeftCharged[9]->Fill(mcTrigPt, mcUexE, GetEventWeight());
          if(mcUexE > 0) fhMCPtHbpXEUeLeftCharged[9]->Fill(mcTrigPt, TMath::Log(1/mcUexE), GetEventWeight());
          
          fhMCPtZTUeLeftCharged[9]->Fill(mcTrigPt, mcUezT, GetEventWeight());
          if(mcUezT > 0) fhMCPtHbpZTUeLeftCharged[9]->Fill(mcTrigPt, TMath::Log(1/mcUezT), GetEventWeight());
        }
      }
    }
  }
  
  return lead;
}

//______________________________________________________________________________________________________________
/// Fill mostly momentum imbalance related histograms.
//______________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedMomentumImbalanceHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                                             Float_t deltaPhi,
                                                                             Int_t   cen,      Int_t   charge,
                                                                             Int_t   bin,      Int_t   decayTag,
                                                                             Int_t   outTOF,   Int_t   mcTag)

{
  Float_t zT =   ptAssoc/ptTrig ;
  Float_t xE   =-ptAssoc/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  Float_t pout = ptAssoc*TMath::Sin(deltaPhi) ;
  
  if(xE < 0.)
    AliWarning(Form("Careful!!, negative xE %2.2f for right UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                    xE,deltaPhi*TMath::RadToDeg(),TMath::Cos(deltaPhi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
  
  Float_t hbpXE = -100;
  Float_t hbpZT = -100;
  
  if(xE > 0 ) hbpXE = TMath::Log(1./xE);
  if(zT > 0 ) hbpZT = TMath::Log(1./zT);
  
  fhXECharged         ->Fill(ptTrig, xE     , GetEventWeight());
  fhPtHbpXECharged    ->Fill(ptTrig, hbpXE  , GetEventWeight());
  fhZTCharged         ->Fill(ptTrig, zT     , GetEventWeight());
  fhPtHbpZTCharged    ->Fill(ptTrig, hbpZT  , GetEventWeight());
  fhPtTrigPout        ->Fill(ptTrig, pout   , GetEventWeight());
  fhPtTrigCharged     ->Fill(ptTrig, ptAssoc, GetEventWeight());
  
  if((deltaPhi > 5*TMath::Pi()/6.)   && (deltaPhi < 7*TMath::Pi()/6.))
  {
    fhXECharged_Cone2         ->Fill(ptTrig, xE   , GetEventWeight());
    fhPtHbpXECharged_Cone2    ->Fill(ptTrig, hbpXE, GetEventWeight());
  }
  
  // MC
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
    fhXEChargedMC[mcIndex]->Fill(ptTrig, xE, GetEventWeight());
    
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
    {
      // check index in GetMCTagIndexHistogram
      if      ( mcIndex == 2 ) fhXEChargedMC[8]->Fill(ptTrig, xE, GetEventWeight()); // pi0 decay
      else if ( mcIndex == 4 ) fhXEChargedMC[9]->Fill(ptTrig, xE, GetEventWeight()); // eta decay
    }
  }
  
  // Pile-up studies
  if(IsPileUpAnalysisOn())
  {
    if     (outTOF==1)
    {
      fhXEChargedOtherBC    ->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedOtherBC    ->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedOtherBC->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    else if(outTOF==0)
    {
      fhXEChargedBC0    ->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedBC0    ->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedBC0->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    
    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhXEChargedVtxBC0    ->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedVtxBC0    ->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedVtxBC0->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    
    if(GetReader()->IsPileUpFromSPD())
    {
      fhXEChargedPileUp    [0]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [0]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[0]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromEMCal())
    {
      fhXEChargedPileUp    [1]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [1]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[1]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDOrEMCal())
    {
      fhXEChargedPileUp    [2]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [2]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[2]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDAndEMCal())
    {
      fhXEChargedPileUp    [3]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [3]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[3]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())
    {
      fhXEChargedPileUp    [4]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [4]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[4]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())
    {
      fhXEChargedPileUp    [5]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [5]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[5]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
    {
      fhXEChargedPileUp    [6]->Fill(ptTrig, xE     , GetEventWeight());
      fhZTChargedPileUp    [6]->Fill(ptTrig, zT     , GetEventWeight());
      fhPtTrigChargedPileUp[6]->Fill(ptTrig, ptAssoc, GetEventWeight());
    }
  }
  
  if(fDecayTrigger && decayTag > 0)
  {
    for(Int_t ibit = 0; ibit<fNDecayBits; ibit++)
    {
      if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
      {
        fhXEDecayCharged[ibit]->Fill(ptTrig, xE, GetEventWeight());
        fhZTDecayCharged[ibit]->Fill(ptTrig, zT, GetEventWeight());
      }
    }
  } // photon decay pi0/eta trigger
  
  if(bin >= 0 && fFillMomImbalancePtAssocBinsHisto)//away side
  {
    fhXEAssocPtBin[bin]->Fill(ptTrig, xE, GetEventWeight()) ;
    fhZTAssocPtBin[bin]->Fill(ptTrig, zT, GetEventWeight()) ;
  }
  
  if(fCorrelVzBin)
  {
    Int_t vz = GetEventVzBin();
    fhXEVZ[vz]->Fill(ptTrig, xE, GetEventWeight()) ;
    fhZTVZ[vz]->Fill(ptTrig, zT, GetEventWeight()) ;
  }
  
  if(charge > 0)
  {
    fhXEPosCharged->Fill(ptTrig, xE, GetEventWeight()) ;
    fhZTPosCharged->Fill(ptTrig, zT, GetEventWeight()) ;
  }
  else
  {
    fhXENegCharged->Fill(ptTrig, xE, GetEventWeight()) ;
    fhZTNegCharged->Fill(ptTrig, zT, GetEventWeight()) ;
  }
  
  // Fill different multiplicity/centrality histogram
  if(IsHighMultiplicityAnalysisOn() && cen >= 0 && cen < GetNCentrBin())
  {
    fhXEMult[cen]->Fill(ptTrig, xE, GetEventWeight());
    fhZTMult[cen]->Fill(ptTrig, zT, GetEventWeight());
  } // multiplicity/centrality events selection
}

//___________________________________________________________________________________________________________________
/// Fill underlying event histograms.
//___________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                                           Float_t deltaPhi, Int_t cen, Int_t outTOF, Int_t   mcTag)
{
  fhUePart->Fill(ptTrig, GetEventWeight());
  
  fhDeltaPhiUeChargedPt->Fill(ptAssoc, deltaPhi, GetEventWeight());
  
  Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
  Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
  Double_t uezT =   ptAssoc/ptTrig;
  
  if(uexE < 0.)
    AliWarning(Form("Careful!!, negative xE %2.2f for right UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                    uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
  
  fhXEUeCharged->Fill(ptTrig, uexE, GetEventWeight());
  if(deltaPhi > TMath::DegToRad()*80 && deltaPhi < TMath::DegToRad()*100)fhXEUeChargedSmallCone->Fill(ptTrig, uexE, GetEventWeight());
  if(deltaPhi > TMath::DegToRad()*70 && deltaPhi < TMath::DegToRad()*110)fhXEUeChargedMediumCone->Fill(ptTrig, uexE, GetEventWeight());
  if(deltaPhi > TMath::DegToRad()*60 && deltaPhi < TMath::DegToRad()*120)fhXEUeChargedLargeCone->Fill(ptTrig, uexE, GetEventWeight());
  
  if(uexE > 0) fhPtHbpXEUeCharged->Fill(ptTrig, TMath::Log(1/uexE), GetEventWeight());
  
  fhZTUeCharged->Fill(ptTrig, uezT, GetEventWeight());
  if(uezT > 0) fhPtHbpZTUeCharged->Fill(ptTrig, TMath::Log(1/uezT), GetEventWeight());
  
  if(IsDataMC())
  {
    Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
    fhXEUeChargedRightMC[mcIndex]->Fill(ptTrig, uexE, GetEventWeight());
    
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
    {
      // check index in GetMCTagIndexHistogram
      if      ( mcIndex == 2 ) fhXEUeChargedRightMC[8]->Fill(ptTrig, uexE, GetEventWeight()); // pi0 decay
      else if ( mcIndex == 4 ) fhXEUeChargedRightMC[9]->Fill(ptTrig, uexE, GetEventWeight()); // eta decay
    }
  }
  
  // Pile-up studies
  
  if(IsPileUpAnalysisOn())
  {
    if     (outTOF==1)
    {
      fhXEUeChargedOtherBC->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedOtherBC->Fill(ptTrig, uezT, GetEventWeight());
    }
    else if(outTOF==0)
    {
      fhXEUeChargedBC0->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedBC0->Fill(ptTrig, uezT, GetEventWeight());
    }
    
    Int_t vtxBC = GetReader()->GetVertexBC();
    if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)
    {
      fhXEUeChargedVtxBC0->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedVtxBC0->Fill(ptTrig, uezT, GetEventWeight());
    }
    
    if(GetReader()->IsPileUpFromSPD())
    {
      fhXEUeChargedPileUp[0]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[0]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromEMCal())
    {
      fhXEUeChargedPileUp[1]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[1]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDOrEMCal())
    {
      fhXEUeChargedPileUp[2]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[2]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDAndEMCal())
    {
      fhXEUeChargedPileUp[3]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[3]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())
    {
      fhXEUeChargedPileUp[4]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[4]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())
    {
      fhXEUeChargedPileUp[5]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[5]->Fill(ptTrig, uezT, GetEventWeight());
    }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
    {
      fhXEUeChargedPileUp[6]->Fill(ptTrig, uexE, GetEventWeight());
      fhZTUeChargedPileUp[6]->Fill(ptTrig, uezT, GetEventWeight());
    }
  }
  
  // Fill different multiplicity/centrality histogram
  if(IsHighMultiplicityAnalysisOn() && cen >= 0 && cen < GetNCentrBin())
  {
    fhXEUeMult[cen]->Fill(ptTrig, uexE, GetEventWeight());
    fhZTUeMult[cen]->Fill(ptTrig, uezT, GetEventWeight());
  }
}

//_______________________________________________________________________________________________
/// Fill underlying event histograms to the left and right of trigger
/// Right cone is the default UE.
//_______________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedUnderlyingEventSidesHistograms(Float_t ptTrig,
                                                                                Float_t ptAssoc,
                                                                                Float_t deltaPhi,
                                                                                Int_t   mcTag)
{
  if((deltaPhi<-fUeDeltaPhiMinCut) || (deltaPhi >2*fUeDeltaPhiMaxCut))
  {
    fhDeltaPhiUeLeftCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    Double_t uezT =   ptAssoc/ptTrig;
    
    if(uexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for left UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhXEUeLeftCharged->Fill(ptTrig, uexE, GetEventWeight());
    if(uexE > 0) fhPtHbpXEUeLeftCharged->Fill(ptTrig, TMath::Log(1/uexE), GetEventWeight());
    
    fhZTUeLeftCharged->Fill(ptTrig, uezT, GetEventWeight());
    if(uezT > 0) fhPtHbpZTUeLeftCharged->Fill(ptTrig, TMath::Log(1/uezT), GetEventWeight());
    
    fhDeltaPhiUeLeftCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    if(IsDataMC())
    {
      Int_t mcIndex = GetMCTagHistogramIndex(mcTag);
      fhXEUeChargedLeftMC[mcIndex]->Fill(ptTrig, uexE, GetEventWeight());
      
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
      {
        // check index in GetMCTagIndexHistogram
        if      ( mcIndex == 2 ) fhXEUeChargedLeftMC[8]->Fill(ptTrig, uexE, GetEventWeight()); // pi0 decay
        else if ( mcIndex == 4 ) fhXEUeChargedLeftMC[9]->Fill(ptTrig, uexE, GetEventWeight()); // eta decay
      }
    }
  }
  
  if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-TMath::Pi()/2))
  {
    fhDeltaPhiUeLeftDownCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for left-down UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhXEUeLeftDownCharged->Fill(ptTrig, uexE, GetEventWeight());
  }
  
  if((deltaPhi>2*fUeDeltaPhiMaxCut) && (deltaPhi <3*TMath::Pi()/2))
  {
    fhDeltaPhiUeLeftUpCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for left-up UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhXEUeLeftUpCharged->Fill(ptTrig, uexE, GetEventWeight());
  }
  
  if((deltaPhi > TMath::Pi()/2) && (deltaPhi < fUeDeltaPhiMaxCut))
  {
    fhDeltaPhiUeRightUpCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for right-up UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhXEUeRightUpCharged->Fill(ptTrig, uexE, GetEventWeight());
  }
  
  if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < TMath::Pi()/2))
  {
    fhDeltaPhiUeRightDownCharged->Fill(ptAssoc, deltaPhi, GetEventWeight());
    
    Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
    Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
    
    if(uexE < 0.)
      AliWarning(Form("Careful!!, negative xE %2.2f for right-down UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                      uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
    
    fhXEUeRightDownCharged->Fill(ptTrig, uexE, GetEventWeight());
  }
}

//_____________________________________________________________________________________________________________________________________
/// Do correlation with decay photons of triggered pi0 or eta.
//_____________________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillDecayPhotonCorrelationHistograms(Float_t ptAssoc, Float_t phiAssoc, Bool_t bChargedOrNeutral)
{
  // Calculate the correlation parameters
  Float_t ptDecay1 = fDecayMom1.Pt();
  Float_t ptDecay2 = fDecayMom2.Pt();
  
  Float_t zTDecay1 = -100, zTDecay2 = -100;
  if(ptDecay1 > 0) zTDecay1 = ptAssoc/ptDecay1 ;
  if(ptDecay2 > 0) zTDecay2 = ptAssoc/ptDecay2 ;
  
  Float_t deltaPhiDecay1 = fDecayMom1.Phi()-phiAssoc;
  if(deltaPhiDecay1< -TMath::PiOver2()) deltaPhiDecay1+=TMath::TwoPi();
  if(deltaPhiDecay1>3*TMath::PiOver2()) deltaPhiDecay1-=TMath::TwoPi();
  
  Float_t deltaPhiDecay2 = fDecayMom2.Phi()-phiAssoc;
  if(deltaPhiDecay2< -TMath::PiOver2()) deltaPhiDecay2+=TMath::TwoPi();
  if(deltaPhiDecay2>3*TMath::PiOver2()) deltaPhiDecay2-=TMath::TwoPi();
  
  Float_t xEDecay1   =-zTDecay1*TMath::Cos(deltaPhiDecay1); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  Float_t xEDecay2   =-zTDecay2*TMath::Cos(deltaPhiDecay2); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  
  if(bChargedOrNeutral) // correlate with charges
  {
    fhDeltaPhiPi0DecayCharged->Fill(ptDecay1, deltaPhiDecay1, GetEventWeight());
    fhDeltaPhiPi0DecayCharged->Fill(ptDecay2, deltaPhiDecay2, GetEventWeight());
    
    AliDebug(2,Form("deltaPhoton1 = %f, deltaPhoton2 = %f", deltaPhiDecay1, deltaPhiDecay2));
    
    if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
    {
      fhZTPi0DecayCharged->Fill(ptDecay1, zTDecay1, GetEventWeight());
      fhXEPi0DecayCharged->Fill(ptDecay1, xEDecay1, GetEventWeight());
    }
    if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
    {
      fhZTPi0DecayCharged->Fill(ptDecay2, zTDecay2, GetEventWeight());
      fhXEPi0DecayCharged->Fill(ptDecay2, xEDecay2, GetEventWeight());
    }
  }
  else // correlate with neutrals
  {
    fhDeltaPhiPi0DecayNeutral->Fill(ptDecay1, deltaPhiDecay1, GetEventWeight());
    fhDeltaPhiPi0DecayNeutral->Fill(ptDecay2, deltaPhiDecay2, GetEventWeight());
    
    AliDebug(2,Form("deltaPhoton1 = %f, deltaPhoton2 = %f", deltaPhiDecay1, deltaPhiDecay2));
    
    if( (deltaPhiDecay1 > fDeltaPhiMinCut) && ( deltaPhiDecay1 < fDeltaPhiMaxCut) )
    {
      fhZTPi0DecayNeutral->Fill(ptDecay1, zTDecay1, GetEventWeight());
      fhXEPi0DecayNeutral->Fill(ptDecay1, xEDecay1, GetEventWeight());
    }
    if( (deltaPhiDecay2 > fDeltaPhiMinCut) && ( deltaPhiDecay2 < fDeltaPhiMaxCut) )
    {
      fhZTPi0DecayNeutral->Fill(ptDecay2, zTDecay2, GetEventWeight());
      fhXEPi0DecayNeutral->Fill(ptDecay2, xEDecay2, GetEventWeight());
    }
  }
}

//_____________________________________________________________________________________________________________________________
/// Fill underlying event histograms to the left of trigger.
/// Right is the default case.
//_____________________________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::FillNeutralUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                                                Float_t zT,       Float_t hbpZT,
                                                                                Float_t deltaPhi)
{
  Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
  
  Float_t xE  =-ptAssoc/ptTrig*TMath::Cos(randomphi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
  Float_t hbpXE = -100;
  if(xE > 0 ) hbpXE = TMath::Log(1./xE);
  
  if((deltaPhi<-fUeDeltaPhiMinCut) && (deltaPhi >-fUeDeltaPhiMaxCut))
  {
    fhDeltaPhiUeLeftNeutral->Fill(ptAssoc, deltaPhi, GetEventWeight());
    fhXEUeLeftNeutral      ->Fill(ptTrig , xE      , GetEventWeight());
    fhPtHbpXEUeLeftNeutral ->Fill(ptTrig , hbpXE   , GetEventWeight());
    fhZTUeLeftNeutral      ->Fill(ptTrig , zT      , GetEventWeight());
    fhPtHbpZTUeLeftNeutral ->Fill(ptTrig , hbpZT   , GetEventWeight());
  }
}

//______________________________________________________
/// Fill the pool with tracks or clusters if requested.
//______________________________________________________
void AliAnaParticleHadronCorrelation::FillEventMixPool()
{
  if ( !DoOwnMix() ) return;
  
  FillChargedEventMixPool();
  
  // Do the cluster pool filling only if requested
  // or in case of isolation cut using clusters in the cone.
  Bool_t isoCase = OnlyIsolated() && (GetIsolationCut()->GetParticleTypeInCone() != AliIsolationCut::kOnlyCharged);
  
  if( !fFillNeutralEventMixPool && !isoCase) return;
  
  FillNeutralEventMixPool();
}

//_____________________________________________________________
/// Mixed event pool filling for tracks.
//_____________________________________________________________
void AliAnaParticleHadronCorrelation::FillChargedEventMixPool()
{
  if(fUseMixStoredInReader && GetReader()->GetLastTracksMixedEvent() == GetEventNumber())
  {
    //printf("%s : Pool already filled for this event !!!\n",GetInputAODName().Data());
    return ; // pool filled previously for another trigger
  }
  
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return ;
  
  // Do mixing only with MB event (or the chosen mask), if not skip
  if( !(inputHandler->IsEventSelected( ) & GetReader()->GetMixEventTriggerMask()) ) return ;
  
  Int_t eventBin = GetEventMixBin();
  
  //Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  fhEventMBBin->Fill(eventBin, GetEventWeight());
  
  TObjArray * mixEventTracks = new TObjArray;
  
  if(fUseMixStoredInReader)
  {
    fListMixTrackEvents[eventBin] = GetReader()->GetListWithMixedEventsForTracks(eventBin);
  }
  
  if(!fListMixTrackEvents[eventBin]) fListMixTrackEvents[eventBin] = new TList();
  
  //printf("%s ***** Pool Event bin : %d - nTracks %d\n",GetInputAODName().Data(),eventBin, GetCTSTracks()->GetEntriesFast());
  
  TList * pool = fListMixTrackEvents[eventBin];
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    Float_t pt   = fTrackVector.Pt();
    
    // Select only hadrons in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    AliAODPWG4Particle * mixedTrack = new AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),0);
    mixedTrack->SetDetectorTag(kCTS);
    mixedTrack->SetChargedBit(track->Charge()>0);
    mixEventTracks->Add(mixedTrack);
  }
  
  fhNtracksMB->Fill(mixEventTracks->GetEntriesFast(), eventBin, GetEventWeight());
  
  // Set the event number where the last event was added, to avoid double pool filling
  GetReader()->SetLastTracksMixedEvent(GetEventNumber());
  
  //printf("Add event to pool with %d tracks \n ",mixEventTracks->GetEntries());
  pool->AddFirst(mixEventTracks);
  mixEventTracks = 0;
  
  //printf("Pool size %d, max %d\n",pool->GetSize(), GetNMaxEvMix());
  
  if(pool->GetSize() > GetNMaxEvMix())
  {
    // Remove last event
    TClonesArray * tmp = static_cast<TClonesArray*>(pool->Last()) ;
    pool->RemoveLast() ;
    delete tmp ;
  }
}

//_____________________________________________________________
/// Mixed event pool filling for neutral clusters.
/// Right now only for EMCAL and in isolation case.
//_____________________________________________________________
void AliAnaParticleHadronCorrelation::FillNeutralEventMixPool()
{
  //printf("FillNeutralEventMixPool for %s\n",GetInputAODName().Data());
  
  if(fUseMixStoredInReader && GetReader()->GetLastCaloMixedEvent() == GetEventNumber())
  {
    //printf("%s : Pool already filled for this event !!!\n",GetInputAODName().Data());
    return ; // pool filled previously for another trigger
  }
  
  TObjArray * pl = GetEMCALClusters();
  //if (GetAODObjArrayName.Contains("PHOS") )pl    = GetPHOSClusters();
  //else                                     pl    = GetEMCALClusters();
  
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return ;
  
  // Do mixing only with MB event (or the chosen mask), if not skip
  if( !(inputHandler->IsEventSelected( ) & GetReader()->GetMixEventTriggerMask()) ) return ;
  
  Int_t eventBin = GetEventMixBin();
  
  // Check that the bin exists, if not (bad determination of RP, centrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  TObjArray * mixEventCalo = new TObjArray;
  
  if(fUseMixStoredInReader)
  {
    fListMixCaloEvents[eventBin] = GetReader()->GetListWithMixedEventsForCalo(eventBin);
  }
  
  if(!fListMixCaloEvents[eventBin]) fListMixCaloEvents[eventBin] = new TList();
  
  TList * poolCalo = fListMixCaloEvents[eventBin];
  
  for(Int_t ipr = 0;ipr <  pl->GetEntriesFast() ; ipr ++ )
  {
    AliVCluster * calo = (AliVCluster *) (pl->At(ipr)) ;
    
    // Remove matched clusters
    if( IsTrackMatched( calo, GetReader()->GetInputEvent() ) ) continue ;
    
    // Cluster momentum calculation
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      calo->GetMomentum(fMomentum,GetVertex(0)) ;
    }// Assume that come from vertex in straight line
    else
    {
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(fMomentum,vertex) ;
    }
    
    Float_t pt = fMomentum.Pt();
    
    // Select only clusters in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    AliAODPWG4Particle * mixedCalo = new AliAODPWG4Particle(fMomentum);
    mixedCalo->SetDetectorTag(kEMCAL);
    mixEventCalo->Add(mixedCalo);
  }
  
  fhNclustersMB->Fill(mixEventCalo->GetEntriesFast(), eventBin, GetEventWeight());
  
  // Set the event number where the last event was added, to avoid double pool filling
  GetReader()->SetLastCaloMixedEvent(GetEventNumber());
  
  //printf("Add event to pool with %d clusters \n ",mixEventCalo->GetEntries());
  poolCalo->AddFirst(mixEventCalo);
  mixEventCalo = 0;
  
  //printf("Pool size %d, max %d\n",poolCalo->GetSize(), GetNMaxEvMix());
  
  if(poolCalo->GetSize() > GetNMaxEvMix())
  {
    // Remove last event
    TClonesArray * tmp = static_cast<TClonesArray*>(poolCalo->Last()) ;
    poolCalo->RemoveLast() ;
    delete tmp ;
  }
}

//_________________________________________________________________________________________________________________
/// Select events where the leading charged particle in the opposite hemisphere
/// to the trigger particle is in a window centered at 180 from the trigger.
//_________________________________________________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::FindLeadingOppositeHadronInWindow(AliAODPWG4ParticleCorrelation * particle)
{
  Float_t etaTrig    = particle->Eta();
  Float_t ptTrig     = particle->Pt();
  Float_t phiTrig    = particle->Phi();
  if(phiTrig < 0 ) phiTrig+= TMath::TwoPi();
  
  Float_t ptLeadHad  = 0 ;
  Float_t dphiLeadHad= -100 ;
  Float_t phiLeadHad = -100 ;
  Float_t etaLeadHad = -100 ;
  Int_t   nTrack     = 0;
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    
    Float_t pt   = fTrackVector.Pt();
    Float_t phi  = fTrackVector.Phi() ;
    if(phi < 0 ) phi+= TMath::TwoPi();
    
    Float_t deltaPhi = phiTrig-phi;
    //
    // Calculate deltaPhi shift so that for the particles on the opposite side
    // it is defined between 90 and 270 degrees
    // Shift [-360,-90]  to [0, 270]
    // and [270,360] to [-90,0]
    if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
    
    if(pt > ptLeadHad && deltaPhi > TMath::PiOver2()) // in opposite hemisphere
    {
      ptLeadHad  = pt ;
      phiLeadHad = phi;
      dphiLeadHad= deltaPhi;
      etaLeadHad = fTrackVector.Eta();
      nTrack++;
    }
  }// track loop
  
  if(fFillLeadHadOppositeHisto)
  {
    if(nTrack == 0)
    {
      fhPtNoLeadingOppositeHadron    ->Fill(ptTrig ,          GetEventWeight());
      fhEtaPhiNoLeadingOppositeHadron->Fill(etaTrig, phiTrig, GetEventWeight());
    }
    else
    {
      fhPtLeadingOppositeHadron       ->Fill(ptTrig,   ptLeadHad,         GetEventWeight());
      fhPtDiffPhiLeadingOppositeHadron->Fill(ptTrig, dphiLeadHad,         GetEventWeight());
      fhPtDiffEtaLeadingOppositeHadron->Fill(ptTrig,  etaLeadHad-etaTrig, GetEventWeight());
    }
  }
  
  
  AliDebug(1,Form("pT %2.2f, phi %2.2f, eta %2.2f, nTracks away %d, total tracks %d",
                  ptLeadHad,phiLeadHad*TMath::RadToDeg(),etaLeadHad,nTrack, GetTrackMultiplicity()));
  AliDebug(1,Form("\t pT trig %2.2f, Dphi (trigger-hadron) %2.2f, Deta (trigger-hadron) %2.2f",
                  ptTrig, dphiLeadHad*TMath::RadToDeg(), etaLeadHad-etaTrig));
  AliDebug(1,Form("\t cuts pT: min %2.2f, max %2.2f; DPhi: min %2.2f, max %2.2f",
                  fMinLeadHadPt,fMaxLeadHadPt,fMinLeadHadPhi*TMath::RadToDeg(),fMaxLeadHadPhi*TMath::RadToDeg()));
  
  
  // reject the trigger if the leading hadron is not in the requested pt or phi window and
  
  if( nTrack == 0 ) return kFALSE; // No track found in opposite hemisphere
  
  if( ptLeadHad < fMinLeadHadPt || ptLeadHad > fMaxLeadHadPt ) return kFALSE;
  
  //printf("Accept leading hadron pT \n");
  
  if( dphiLeadHad < fMinLeadHadPhi || dphiLeadHad > fMaxLeadHadPhi ) return kFALSE;
  
  //printf("Accept leading hadron phi \n");
  
  return kTRUE ;
}

//____________________________________________________________
/// Save parameters used for analysis.
//____________________________________________________________
TObjString* AliAnaParticleHadronCorrelation::GetAnalysisCuts()
{
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 560;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPaticleHadronCorrelation ---:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Pt associated < %3.2f; ", fMinAssocPt,   fMaxAssocPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Phi trigger particle-Hadron < %3.2f; ",    fDeltaPhiMinCut,   fDeltaPhiMaxCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," %3.2f < Phi trigger particle-UeHadron <  %3.2f; ", fUeDeltaPhiMinCut, fUeDeltaPhiMaxCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Isolated Trigger?  %d;",    fSelectIsolated) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Several UE?  %d;",          fMakeSeveralUE) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Name of AOD Pi0 Branch %s;", fPi0AODBranchName.Data());
  parList+=onePar ;
  snprintf(onePar,buffersize,"Do Decay-hadron correlation ?  pi0 %d, decay %d;", fPi0Trigger, fDecayTrigger) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Select absolute leading for cluster triggers ? %d or Near Side Leading %d;",
           fMakeAbsoluteLeading, fMakeNearSideLeading) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Associated particle pt bins  %d: ", fNAssocPtBins) ;
  parList+=onePar ;
  for (Int_t ibin = 0; ibin<fNAssocPtBins; ibin++) {
    snprintf(onePar,buffersize,"bin %d = [%2.1f,%2.1f];", ibin, fAssocPtBinLimit[ibin], fAssocPtBinLimit[ibin+1]) ;
  }
  parList+=onePar ;
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList()
  
  return new TObjString(parList) ;
}

//________________________________________________________________
/// Create histograms to be saved in output file and
/// store them in fOutputContainer.
//________________________________________________________________
TList *  AliAnaParticleHadronCorrelation::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName("CorrelationHistos") ;
  
  Int_t   nptbins = GetHistogramRanges()->GetHistoPtBins(); Int_t  nphibins = GetHistogramRanges()->GetHistoPhiBins(); Int_t   netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax   = GetHistogramRanges()->GetHistoPtMax();  Float_t phimax  = GetHistogramRanges()->GetHistoPhiMax();  Float_t etamax   = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin   = GetHistogramRanges()->GetHistoPtMin();  Float_t phimin  = GetHistogramRanges()->GetHistoPhiMin();  Float_t etamin   = GetHistogramRanges()->GetHistoEtaMin();
  
  Int_t  ndeltaphibins = GetHistogramRanges()->GetHistoDeltaPhiBins(); Int_t   ndeltaetabins = GetHistogramRanges()->GetHistoDeltaEtaBins();
  Float_t deltaphimax  = GetHistogramRanges()->GetHistoDeltaPhiMax();  Float_t deltaetamax   = GetHistogramRanges()->GetHistoDeltaEtaMax();
  Float_t deltaphimin  = GetHistogramRanges()->GetHistoDeltaPhiMin();  Float_t deltaetamin   = GetHistogramRanges()->GetHistoDeltaEtaMin();
  
  Int_t ntrbins = GetHistogramRanges()->GetHistoTrackMultiplicityBins(); Int_t nclbins = GetHistogramRanges()->GetHistoNClustersBins();
  Int_t trmax   = GetHistogramRanges()->GetHistoTrackMultiplicityMax();  Int_t clmax   = GetHistogramRanges()->GetHistoNClustersMax();
  Int_t trmin   = GetHistogramRanges()->GetHistoTrackMultiplicityMin();  Int_t clmin   = GetHistogramRanges()->GetHistoNClustersMin();
  
  Int_t  nxeztbins = GetHistogramRanges()->GetHistoRatioBins();  Int_t  nhbpbins = GetHistogramRanges()->GetHistoHBPBins();
  Float_t xeztmax  = GetHistogramRanges()->GetHistoRatioMax();   Float_t hbpmax  = GetHistogramRanges()->GetHistoHBPMax();
  Float_t xeztmin  = GetHistogramRanges()->GetHistoRatioMin();   Float_t hbpmin  = GetHistogramRanges()->GetHistoHBPMin();
  
  Int_t nMixBins = GetNCentrBin()*GetNZvertBin()*GetNRPBin();
  
  Int_t   nmassbins = GetHistogramRanges()->GetHistoMassBins();
  Float_t massmin   = GetHistogramRanges()->GetHistoMassMin();
  Float_t massmax   = GetHistogramRanges()->GetHistoMassMax();
  
  TString nameMC[]     = {"Photon","Pi0","Pi0Decay","Eta","EtaDecay","OtherDecay","Electron","Hadron","Pi0DecayLostPair","EtaDecayLostPair"};
  TString mcPartType[] = {"#gamma", "#pi^{0} (merged #gamma)", "#gamma_{#pi decay}","#eta (merged #gamma)"   , "#gamma_{#eta decay}", "#gamma_{other decay}", 
    "e^{#pm}" , "hadrons?" , "#gamma_{#pi decay} lost companion", "#gamma_{#eta decay} lost companion"} ;
  TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
  
  TString parTitle = Form("#it{R} = %2.2f",GetIsolationCut()->GetConeSize());
  
  // For vz dependent histograms, if option ON
  Int_t   nz  = 1  ;
  if(fCorrelVzBin) nz = GetNZvertBin();
  TString sz  = "" ;
  TString tz  = "" ;
  
  // Fill histograms for neutral clusters in mixing?
  Bool_t isoCase = OnlyIsolated() && (GetIsolationCut()->GetParticleTypeInCone() != AliIsolationCut::kOnlyCharged);
  Bool_t neutralMix = fFillNeutralEventMixPool || isoCase ;
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kESD )
  {
    fhTrackResolution  = new TH2F ("hTrackResolution","Track resolution: #sigma_{#it{p}_{T}} vs #it{p}_{T}, away side, ESDs", 
                                   nptbins,ptmin,ptmax,600,0,0.3);
    fhTrackResolution->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTrackResolution->SetYTitle("#sigma_{#it{p}_{T}} / #it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhTrackResolution);

    fhTrackResolutionUE  = new TH2F ("hTrackResolutionUE","Track resolution: #sigma_{#it{p}_{T}} vs #it{p}_{T}, UE, ESDs", 
                                   nptbins,ptmin,ptmax,600,0,0.3);
    fhTrackResolutionUE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTrackResolutionUE->SetYTitle("#sigma_{#it{p}_{T}} / #it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhTrackResolutionUE);
  }
  
  fhPtTriggerInput  = new TH1F("hPtTriggerInput","Input trigger #it{p}_{T}", nptbins,ptmin,ptmax);
  fhPtTriggerInput->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
  outputContainer->Add(fhPtTriggerInput);
  
  if( fM02MaxCut > 0 && fM02MinCut > 0 )
  {
    fhPtTriggerSSCut  = new TH1F("hPtTriggerSSCut","Trigger #it{p}_{T} after #lambda^{2}_{0} cut", nptbins,ptmin,ptmax);
    fhPtTriggerSSCut->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    outputContainer->Add(fhPtTriggerSSCut);
  }
  
  if( OnlyIsolated() )
  {
    fhPtTriggerIsoCut  = new TH1F("hPtTriggerIsoCut","Trigger #it{p}_{T} after isolation (and #lambda^{2}_{0} cut)", nptbins,ptmin,ptmax);
    fhPtTriggerIsoCut->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    outputContainer->Add(fhPtTriggerIsoCut);
  }
  
  fhPtTriggerFidCut  = new TH1F("hPtTriggerFidCut","Trigger #it{p}_{T} after fiducial (isolation and #lambda^{2}_{0}) cut", nptbins,ptmin,ptmax);
  fhPtTriggerFidCut->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
  outputContainer->Add(fhPtTriggerFidCut);
  
  fhPtTrigger  = new TH1F("hPtTrigger","#it{p}_{T} distribution of trigger particles (after opposite hadron leading cut and rest)", nptbins,ptmin,ptmax);
  fhPtTrigger->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
  outputContainer->Add(fhPtTrigger);
  
  if(fFillInvMassHisto)
  {
    fhMassPtTrigger  = new TH2F("hMassPtTrigger","2 photons invariant mass vs p_{T}^{trig}",
                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
    fhMassPtTrigger->SetYTitle("M_{#gamma#gamma} (GeV/#it{c}^{2})");
    fhMassPtTrigger->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhMassPtTrigger);
  }
  
  if(fFillBkgBinsHisto)
  {
    
    fhPtLeadInConeBin    = new TH1F*[fNBkgBin] ;
    
    fhPtSumInConeBin     = new TH1F*[fNBkgBin] ;
    
    fhPtLeadConeBinDecay = new TH1F*[fNBkgBin*fNDecayBits] ;
    
    fhSumPtConeBinDecay  = new TH1F*[fNBkgBin*fNDecayBits] ;
    
    fhPtLeadConeBinMC    = new TH1F*[fNBkgBin*fgkNmcTypes] ;
    
    fhSumPtConeBinMC     = new TH1F*[fNBkgBin*fgkNmcTypes] ;
    
    for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
    {
      fhPtLeadInConeBin[ibin]  = new TH1F
      (Form("hPtLeadCone_Bin%d",ibin),
       Form("cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
            fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
      fhPtLeadInConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
      fhPtLeadInConeBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLeadInConeBin[ibin]) ;
      
      fhPtSumInConeBin[ibin]  = new TH1F
      (Form("hSumPtCone_Bin%d",ibin),
       Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
            fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
      fhPtSumInConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
      fhPtSumInConeBin[ibin]->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtSumInConeBin[ibin]) ;
      
      if(fFillTaggedDecayHistograms)
      {
        for(Int_t idecay = 0; idecay < fNDecayBits; idecay++)
        {
          Int_t bindecay = ibin+idecay*fNBkgBin;
          
          fhPtLeadConeBinDecay[bindecay]  = new TH1F
          (Form("hPtLeadCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
           Form("Decay bit %d, cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
                fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
          fhPtLeadConeBinDecay[bindecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhPtLeadConeBinDecay[bindecay]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLeadConeBinDecay[bindecay]) ;
          
          fhSumPtConeBinDecay[bindecay]  = new TH1F
          (Form("hSumPtCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
           Form("Decay bit %d, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c},  %s",
                fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
          fhSumPtConeBinDecay[bindecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhSumPtConeBinDecay[bindecay]->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSumPtConeBinDecay[bindecay]) ;
        }
      }
      
      if(IsDataMC())
      {
        for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
        {
          Int_t binmc = ibin+imc*fNBkgBin;
          fhPtLeadConeBinMC[binmc]  = new TH1F
          (Form("hPtLeadCone_Bin%d_MC%s",ibin, nameMC[imc].Data()),
           Form("in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, MC %s, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax);
          fhPtLeadConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhPtLeadConeBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLeadConeBinMC[binmc]) ;
          
          fhSumPtConeBinMC[binmc]  = new TH1F
          (Form("hSumPtCone_Bin%d_MC%s",ibin,nameMC[imc].Data()),
           Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax);
          fhSumPtConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhSumPtConeBinMC[binmc]->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSumPtConeBinMC[binmc]) ;
        } // MC particle loop
      }
    }
  }
  
  if(IsDataMC())
  {
    for(Int_t i=0; i < fgkNmcTypes; i++)
    {
      fhPtTriggerMC[i]  = new TH1F(Form("hPtTrigger_MC%s",nameMC[i].Data()),
                                   Form("#it{p}_{T} distribution of trigger particles, trigger origin is %s",nameMC[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhPtTriggerMC[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      outputContainer->Add(fhPtTriggerMC[i]);
      
      if(fFillInvMassHisto)
      {
        fhMCMassPtTrigger[i]  = new TH2F(Form("hMassPtTrigger_MC%s",nameMC[i].Data()),
                                         Form("2 photons invariant mass, trigger origin is %s",nameMC[i].Data()),
                                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax);
        fhMCMassPtTrigger[i]->SetYTitle("m_{#gamma#gamma}");
        fhMCMassPtTrigger[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        outputContainer->Add(fhMCMassPtTrigger[i]) ;
      }
    }
  }
  
  if(fDecayTrigger)
  {
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      fhPtDecayTrigger[ibit]  = new TH1F(Form("hPtDecayTrigger_bit%d",fDecayBits[ibit]),
                                         Form("#it{p}_{T} distribution of trigger particles, decay Bit %d",fDecayBits[ibit]),
                                         nptbins,ptmin,ptmax);
      fhPtDecayTrigger[ibit]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      outputContainer->Add(fhPtDecayTrigger[ibit]);
      
      if(IsDataMC())
      {
        for(Int_t i=0; i < fgkNmcTypes; i++)
        {
          fhPtDecayTriggerMC[ibit][i]  = new TH1F(Form("hPtDecayTrigger_bit%d_MC%s",fDecayBits[ibit], nameMC[i].Data()),
                                                  Form("#it{p}_{T} distribution of trigger particles, decay Bit %d, trigger origin is %s",fDecayBits[ibit], nameMC[i].Data()),
                                                  nptbins,ptmin,ptmax);
          fhPtDecayTriggerMC[ibit][i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          outputContainer->Add(fhPtDecayTriggerMC[ibit][i]);
        }
      }
    }
  }
  
  if(fCorrelVzBin)
  {
    fhPtTriggerVzBin  = new TH2F("hPtTriggerVzBin","#it{p}_{T} distribution of trigger particles vs vz bin", nptbins,ptmin,ptmax,GetNZvertBin(),0,GetNZvertBin());
    fhPtTriggerVzBin->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtTriggerVzBin->SetYTitle("#it{v}_{#it{z}} bin");
    outputContainer->Add(fhPtTriggerVzBin);
  }
  
  fhPtTriggerBin  = new TH2F ("hPtTriggerBin","#it{p}_{T} distribution of trigger particles", nptbins,ptmin,ptmax,nMixBins,0,nMixBins);
  fhPtTriggerBin->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
  fhPtTriggerBin->SetYTitle("Bin");
  outputContainer->Add(fhPtTriggerBin);
  
  fhPhiTrigger  = new TH2F ("hPhiTrigger","#varphi distribution of trigger Particles",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
  fhPhiTrigger->SetYTitle("#varphi (rad)");
  outputContainer->Add(fhPhiTrigger);
  
  fhEtaTrigger  = new TH2F ("hEtaTrigger","#eta distribution of trigger",nptbins,ptmin,ptmax, netabins,etamin,etamax);
  fhEtaTrigger->SetYTitle("#eta ");
  outputContainer->Add(fhEtaTrigger);
  
  if(IsHighMultiplicityAnalysisOn())
  {
    fhPtTriggerCentrality   = new TH2F("hPtTriggerCentrality","Trigger particle #it{p}_{T} vs centrality",nptbins,ptmin,ptmax,100,0.,100) ;
    fhPtTriggerCentrality->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtTriggerCentrality->SetYTitle("Centrality (%)");
    outputContainer->Add(fhPtTriggerCentrality) ;
    
    fhPtTriggerEventPlane  = new TH2F("hPtTriggerEventPlane","Trigger particle #it{p}_{T} vs event plane angle",nptbins,ptmin,ptmax, 100,0.,TMath::Pi()) ;
    fhPtTriggerEventPlane->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtTriggerEventPlane->SetXTitle("EP angle (rad)");
    outputContainer->Add(fhPtTriggerEventPlane) ;
    
    fhTriggerEventPlaneCentrality  = new TH2F("hTriggerEventPlaneCentrality","Trigger particle centrality vs event plane angle",100,0.,100,100,0.,TMath::Pi()) ;
    fhTriggerEventPlaneCentrality->SetXTitle("Centrality (%)");
    fhTriggerEventPlaneCentrality->SetYTitle("EP angle (rad)");
    outputContainer->Add(fhTriggerEventPlaneCentrality) ;
  }
  
  // Leading hadron in oposite side
  if(fFillLeadHadOppositeHisto)
  {
    fhPtLeadingOppositeHadron  = new TH2F("hPtTriggerPtLeadingOppositeHadron","Leading hadron opposite to trigger vs trigger #it{p}_{T}",
                                          nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPtLeadingOppositeHadron->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtLeadingOppositeHadron->SetYTitle("#it{p}_{T}^{lead hadron} (GeV/#it{c})");
    outputContainer->Add(fhPtLeadingOppositeHadron);
    
    fhPtNoLeadingOppositeHadron  = new TH1F("hPtTriggerNoLeadingOppositeHadron","No Leading hadron opposite to trigger #it{p}_{T}",
                                            nptbins,ptmin,ptmax);
    fhPtNoLeadingOppositeHadron->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    outputContainer->Add(fhPtNoLeadingOppositeHadron);
    
    fhEtaPhiNoLeadingOppositeHadron  = new TH2F("hEtaPhiTriggerNoLeadingOppositeHadron","No Leading hadron opposite to trigger #eta:#varphi",
                                                netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiNoLeadingOppositeHadron->SetXTitle("#eta");
    fhEtaPhiNoLeadingOppositeHadron->SetYTitle("#varphi");
    outputContainer->Add(fhEtaPhiNoLeadingOppositeHadron);
    
    
    fhPtDiffPhiLeadingOppositeHadron  = new TH2F("hPtTriggerDiffPhiTriggerLeadingOppositeHadron","#varphi_{trigger}-#varphi_{leading opposite hadron} vs #it{p}_{T}^{trig}",
                                                 nptbins,ptmin,ptmax,ndeltaphibins,deltaphimin,deltaphimax);
    fhPtDiffPhiLeadingOppositeHadron->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtDiffPhiLeadingOppositeHadron->SetYTitle("#varphi_{trigger}-#varphi_{leading opposite hadron} (rad)");
    outputContainer->Add(fhPtDiffPhiLeadingOppositeHadron);
    
    fhPtDiffEtaLeadingOppositeHadron  = new TH2F("hPtTriggerDiffEtaTriggerPhiLeadingOppositeHadron","#eta_{trigger}-#eta_{leading opposite hadron} vs #it{p}_{T}^{trig}",
                                                 nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);
    fhPtDiffEtaLeadingOppositeHadron->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtDiffEtaLeadingOppositeHadron->SetYTitle("#eta_{trigger}-#eta_{leading opposite hadron}");
    outputContainer->Add(fhPtDiffEtaLeadingOppositeHadron);
  }
  
  // Correlation with charged hadrons
  
  fhDeltaPhiDeltaEtaCharged  = new TH2F
  ("hDeltaPhiDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs #varphi_{trigger} - #varphi_{h^{#pm}}",
   ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax);
  fhDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #varphi (rad)");
  fhDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");
  
  fhDeltaPhiDeltaEtaChargedPtA3GeV  = new TH2F
  ("hDeltaPhiDeltaEtaChargedPtA3GeV","#eta_{trigger} - #eta_{h^{#pm}} vs #varphi_{trigger} - #varphi_{h^{#pm}, #it{p}_{TA}>3 GeV/#it{c}}",
   ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax);
  fhDeltaPhiDeltaEtaChargedPtA3GeV->SetXTitle("#Delta #varphi (rad)");
  fhDeltaPhiDeltaEtaChargedPtA3GeV->SetYTitle("#Delta #eta");
  
  fhPhiCharged  = new TH2F
  ("hPhiCharged","#varphi_{h^{#pm}}  vs #it{p}_{T #pm}",
   nptbins,ptmin,ptmax,180,0,TMath::TwoPi());
  fhPhiCharged->SetYTitle("#varphi_{h^{#pm}} (rad)");
  fhPhiCharged->SetXTitle("#it{p}_{T #pm} (GeV/#it{c})");
  
  fhEtaCharged  = new TH2F
  ("hEtaCharged","#eta_{h^{#pm}}  vs #it{p}_{T #pm}",
   nptbins,ptmin,ptmax,100,-1.,1.);
  fhEtaCharged->SetYTitle("#eta_{h^{#pm}} (rad)");
  fhEtaCharged->SetXTitle("#it{p}_{T #pm} (GeV/#it{c})");
  
  fhDeltaPhiCharged  = new TH2F
  ("hDeltaPhiCharged","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}",
   nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
  fhDeltaPhiCharged->SetYTitle("#Delta #varphi (rad)");
  fhDeltaPhiCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhDeltaPhiChargedPtA3GeV  = new TH2F
  ("hDeltaPhiChargedPtA3GeV","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}",
   nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
  fhDeltaPhiChargedPtA3GeV->SetYTitle("#Delta #varphi (rad)");
  fhDeltaPhiChargedPtA3GeV->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  
  fhDeltaPhiChargedPt  = new TH2F
  ("hDeltaPhiChargedPt","#varphi_{trigger} - #varphi_{#h^{#pm}} vs #it{p}_{T h^{#pm}}",
   nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
  fhDeltaPhiChargedPt->SetYTitle("#Delta #varphi (rad)");
  fhDeltaPhiChargedPt->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
  
  fhDeltaEtaCharged  = new TH2F
  ("hDeltaEtaCharged","#eta_{trigger} - #eta_{h^{#pm}} vs #it{p}_{T trigger}",
   nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);
  fhDeltaEtaCharged->SetYTitle("#Delta #eta");
  fhDeltaEtaCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhDeltaEtaChargedPtA3GeV  = new TH2F
  ("hDeltaEtaChargedPtA3GeV","#eta_{trigger} - #eta_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}",
   nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);
  fhDeltaEtaChargedPtA3GeV->SetYTitle("#Delta #eta");
  fhDeltaEtaChargedPtA3GeV->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXECharged  =
  new TH2F("hXECharged","#it{x}_{#it{E}} for charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXECharged->SetYTitle("#it{x}_{#it{E}}");
  fhXECharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXECharged_Cone2  =
  new TH2F("hXECharged_Cone2","#it{x}_{#it{E}} for charged tracks in cone 2 (5#pi/6-7#pi/6)",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXECharged_Cone2->SetYTitle("#it{x}_{#it{E}}");
  fhXECharged_Cone2->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXEPosCharged  =
  new TH2F("hXEPositiveCharged","#it{x}_{#it{E}} for positive charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXEPosCharged->SetYTitle("#it{x}_{#it{E}}");
  fhXEPosCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXENegCharged  =
  new TH2F("hXENegativeCharged","#it{x}_{#it{E}} for negative charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXENegCharged->SetYTitle("#it{x}_{#it{E}}");
  fhXENegCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtHbpXECharged  =
  new TH2F("hHbpXECharged","#xi = ln(1/#it{x}_{#it{E}}) with charged hadrons",
           nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
  fhPtHbpXECharged->SetYTitle("ln(1/#it{x}_{#it{E}})");
  fhPtHbpXECharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtHbpXECharged_Cone2  =
  new TH2F("hHbpXECharged_Cone2","#xi = ln(1/#it{x}_{#it{E}}) with charged hadrons in cone 2 (5#pi/6-7#pi/6)",
           nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
  fhPtHbpXECharged_Cone2->SetYTitle("ln(1/#it{x}_{#it{E}})");
  fhPtHbpXECharged_Cone2->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhZTCharged  =
  new TH2F("hZTCharged","#it{z}_{T} for charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhZTCharged->SetYTitle("#it{z}_{T}");
  fhZTCharged->SetXTitle("#it{p}_{T trigger}");
  
  fhZTPosCharged  =
  new TH2F("hZTPositiveCharged","#it{z}_{T} for positive charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhZTPosCharged->SetYTitle("#it{z}_{T}");
  fhZTPosCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhZTNegCharged  =
  new TH2F("hZTNegativeCharged","#it{z}_{T} for negative charged tracks",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhZTNegCharged->SetYTitle("#it{z}_{T}");
  fhZTNegCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtHbpZTCharged  =
  new TH2F("hHbpZTCharged","#xi = ln(1/#it{z}_{T}) with charged hadrons",
           nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
  fhPtHbpZTCharged->SetYTitle("ln(1/#it{z}_{T})");
  fhPtHbpZTCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtTrigPout  =
  new TH2F("hPtTrigPout","Pout with triggers",
           nptbins,ptmin,ptmax,nptbins,-1.*ptmax/2.,ptmax/2.);
  fhPtTrigPout->SetYTitle("#it{p}_{out} (GeV/#it{c})");
  fhPtTrigPout->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtTrigCharged  =
  new TH2F("hPtTrigCharged","trigger and charged tracks pt distribution",
           nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fhPtTrigCharged->SetYTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
  fhPtTrigCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  outputContainer->Add(fhDeltaPhiDeltaEtaCharged);
  outputContainer->Add(fhDeltaPhiDeltaEtaChargedPtA3GeV);
  outputContainer->Add(fhPhiCharged) ;
  outputContainer->Add(fhEtaCharged) ;
  outputContainer->Add(fhDeltaPhiCharged) ;
  outputContainer->Add(fhDeltaPhiChargedPtA3GeV) ;
  outputContainer->Add(fhDeltaEtaCharged) ;
  outputContainer->Add(fhDeltaEtaChargedPtA3GeV) ;
  outputContainer->Add(fhDeltaPhiChargedPt) ;
  
  outputContainer->Add(fhXECharged) ;
  outputContainer->Add(fhXECharged_Cone2) ;
  
  if(IsDataMC())
  {
    for(Int_t i=0; i < fgkNmcTypes; i++)
    {
      
      fhDeltaPhiChargedMC[i]  = new TH2F(Form("hDeltaPhiCharged_MC%s",nameMC[i].Data()),
                                         Form("#Delta #varphi for charged tracks, trigger origin is %s",nameMC[i].Data()),
                                         nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedMC[i]->SetYTitle("#it{x}_{#it{E}}");
      fhDeltaPhiChargedMC[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhDeltaPhiChargedMC[i]) ;
      
      fhXEChargedMC[i]  = new TH2F(Form("hXECharged_MC%s",nameMC[i].Data()),
                                   Form("#it{x}_{#it{E}} for charged tracks, trigger origin is %s",nameMC[i].Data()),
                                   nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEChargedMC[i]->SetYTitle("#it{x}_{#it{E}}");
      fhXEChargedMC[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEChargedMC[i]) ;
      
      fhXEUeChargedRightMC[i]  = new TH2F(Form("hXEUeChargedRight_MC%s",nameMC[i].Data()),
                                          Form("#it{x}_{#it{E}} for charged tracks in right UE cone, trigger origin is %s",nameMC[i].Data()),
                                          nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEUeChargedRightMC[i]->SetYTitle("#it{x}_{#it{E}}");
      fhXEUeChargedRightMC[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEUeChargedRightMC[i]) ;
      
      fhXEUeChargedLeftMC[i]  = new TH2F(Form("hXEUeChargedLeft_MC%s",nameMC[i].Data()),
                                         Form("#it{x}_{#it{E}} for charged tracks in left UE cone, trigger origin is %s",nameMC[i].Data()),
                                         nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEUeChargedLeftMC[i]->SetYTitle("#it{x}_{#it{E}}");
      fhXEUeChargedLeftMC[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEUeChargedLeftMC[i]) ;
    }
  }
  
  outputContainer->Add(fhXEPosCharged) ;
  outputContainer->Add(fhXENegCharged) ;
  outputContainer->Add(fhPtHbpXECharged) ;
  outputContainer->Add(fhPtHbpXECharged_Cone2) ;
  
  outputContainer->Add(fhZTCharged) ;
  outputContainer->Add(fhZTPosCharged) ;
  outputContainer->Add(fhZTNegCharged) ;
  outputContainer->Add(fhPtHbpZTCharged) ;
  
  outputContainer->Add(fhPtTrigPout) ;
  outputContainer->Add(fhPtTrigCharged) ;
  
  TString right = "";
  if(fMakeSeveralUE) right = "Right";
  
  fhUePart  =  new TH1F("hUePart","UE particles distribution vs pt trig",
                        nptbins,ptmin,ptmax);
  fhUePart->SetYTitle("dNch");
  fhUePart->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhDeltaPhiUeChargedPt  = new TH2F
  (Form("hDeltaPhiUe%sChargedPt",right.Data()),"#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}}",
   nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
  fhDeltaPhiUeChargedPt->SetYTitle("#Delta #varphi (rad)");
  fhDeltaPhiUeChargedPt->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
  
  fhXEUeCharged  =
  new TH2F(Form("hXEUeCharged%s",right.Data()),"#it{x}_{#it{E}} for Underlying Event",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXEUeCharged->SetYTitle("#it{x}_{#it{E}}");
  fhXEUeCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXEUeChargedSmallCone  =
  new TH2F(Form("hXEUeChargedSmallCone%s",right.Data()),"#it{x}_{#it{E}} for Underlying Event in a cone [80,100] deg",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXEUeChargedSmallCone->SetYTitle("#it{x}_{#it{E}}");
  fhXEUeChargedSmallCone->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXEUeChargedMediumCone  =
  new TH2F(Form("hXEUeChargedMediumCone%s",right.Data()),"#it{x}_{#it{E}} for Underlying Event in a cone [70,110] deg",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXEUeChargedMediumCone->SetYTitle("#it{x}_{#it{E}}");
  fhXEUeChargedMediumCone->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhXEUeChargedLargeCone  =
  new TH2F(Form("hXEUeChargedLargeCone%s",right.Data()),"#it{x}_{#it{E}} for Underlying Event in a cone [60,120] deg",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhXEUeChargedLargeCone->SetYTitle("#it{x}_{#it{E}}");
  fhXEUeChargedLargeCone->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtHbpXEUeCharged  =
  new TH2F(Form("hHbpXEUeCharged%s",right.Data()),"#xi = ln(1/#it{x}_{#it{E}}) for Underlying Event",
           nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
  fhPtHbpXEUeCharged->SetYTitle("ln(1/#it{x}_{#it{E}})");
  fhPtHbpXEUeCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhZTUeCharged  =
  new TH2F(Form("hZTUeCharged%s",right.Data()),"#it{z}_{T} for Underlying Event",
           nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
  fhZTUeCharged->SetYTitle("#it{z}_{T}");
  fhZTUeCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  fhPtHbpZTUeCharged  =
  new TH2F(Form("hHbpZTUeCharged%s",right.Data()),"#xi = ln(1/#it{z}_{T}) for Underlying Event",
           nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
  fhPtHbpZTUeCharged->SetYTitle("ln(1/#it{x}_{#it{E}})");
  fhPtHbpZTUeCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
  
  outputContainer->Add(fhUePart);
  outputContainer->Add(fhDeltaPhiUeChargedPt) ;
  outputContainer->Add(fhXEUeCharged) ;
  outputContainer->Add(fhXEUeChargedSmallCone) ;
  outputContainer->Add(fhXEUeChargedMediumCone) ;
  outputContainer->Add(fhXEUeChargedLargeCone) ;
  outputContainer->Add(fhPtHbpXEUeCharged) ;
  outputContainer->Add(fhZTUeCharged) ;
  outputContainer->Add(fhPtHbpZTUeCharged) ;
  
  if(fMakeSeveralUE)
  {
    fhDeltaPhiUeLeftCharged  = new TH2F
    ("hDeltaPhiUeLeftChargedPt","#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}} with UE left side range of trigger particles",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeLeftCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeLeftCharged->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    outputContainer->Add(fhDeltaPhiUeLeftCharged) ;
    
    fhDeltaPhiUeLeftUpCharged  = new TH2F
    ("hDeltaPhiUeLeftUpChargedPt","#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}} with UE left Up side range of trigger particles",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeLeftUpCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeLeftUpCharged->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    outputContainer->Add(fhDeltaPhiUeLeftUpCharged) ;
    
    fhDeltaPhiUeRightUpCharged  = new TH2F
    ("hDeltaPhiUeRightUpChargedPt","#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}} with UE right Up side range of trigger particles",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeRightUpCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeRightUpCharged->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    outputContainer->Add(fhDeltaPhiUeRightUpCharged) ;
    
    fhDeltaPhiUeLeftDownCharged  = new TH2F
    ("hDeltaPhiUeLeftDownChargedPt","#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}} with UE left Down side range of trigger particles",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeLeftDownCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeLeftDownCharged->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    outputContainer->Add(fhDeltaPhiUeLeftDownCharged) ;
    
    fhDeltaPhiUeRightDownCharged  = new TH2F
    ("hDeltaPhiUeRightDownChargedPt","#varphi_{trigger} - #varphi_{#Ueh^{#pm}} vs #it{p}_{T Ueh^{#pm}} with UE right Down side range of trigger particles",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeRightDownCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeRightDownCharged->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    outputContainer->Add(fhDeltaPhiUeRightDownCharged) ;
    
    fhXEUeLeftCharged  =
    new TH2F("hXEUeChargedLeft","#it{x}_{#it{E}} with UE left side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeLeftCharged->SetYTitle("#it{x}_{#it{E} Ueh^{#pm}}");
    fhXEUeLeftCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhXEUeLeftCharged) ;
    
    fhXEUeLeftUpCharged  =
    new TH2F("hXEUeChargedLeftUp","#it{x}_{#it{E}} with UE left Up side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeLeftUpCharged->SetYTitle("#it{x}_{#it{E} Ueh^{#pm}}");
    fhXEUeLeftUpCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhXEUeLeftUpCharged) ;
    
    fhXEUeRightUpCharged  =
    new TH2F("hXEUeChargedRightUp","#it{x}_{#it{E} h^{#pm}} with UE right Up side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeRightUpCharged->SetYTitle("#it{z}_{trigger Ueh^{#pm}}");
    fhXEUeRightUpCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhXEUeRightUpCharged) ;
    
    fhXEUeLeftDownCharged  =
    new TH2F("hXEUeChargedLeftDown","#it{x}_{#it{E}} with UE left Down side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeLeftDownCharged->SetYTitle("#it{x}_{#it{E} Ueh^{#pm}}");
    fhXEUeLeftDownCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhXEUeLeftDownCharged) ;
    
    fhXEUeRightDownCharged  =
    new TH2F("hXEUeChargedRightDown","#it{x}_{#it{E} h^{#pm}} with UE right Down side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeRightDownCharged->SetYTitle("#it{z}_{trigger Ueh^{#pm}}");
    fhXEUeRightDownCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhXEUeRightDownCharged) ;
    
    fhPtHbpXEUeLeftCharged  =
    new TH2F("hHbpXEUeChargedLeft","#xi = ln(1/#it{x}_{#it{E}}) with charged UE left side of trigger",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpXEUeLeftCharged->SetYTitle("ln(1/#it{x}_{#it{E}})");
    fhPtHbpXEUeLeftCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhPtHbpXEUeLeftCharged) ;
    
    fhZTUeLeftCharged  =
    new TH2F("hZTUeChargedLeft","#it{z}_{trigger h^{#pm}} = #it{p}_{T Ueh^{#pm}} / #it{p}_{T trigger} with UE left side of trigger",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTUeLeftCharged->SetYTitle("#it{z}_{trigger Ueh^{#pm}}");
    fhZTUeLeftCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhZTUeLeftCharged) ;
    
    fhPtHbpZTUeLeftCharged  =
    new TH2F("hHbpZTUeChargedLeft","#xi = ln(1/#it{z}_{T}) with charged UE left side of trigger",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpZTUeLeftCharged->SetYTitle("ln(1/#it{z}_{T})");
    fhPtHbpZTUeLeftCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhPtHbpZTUeLeftCharged) ;
  }
  
  if(IsPileUpAnalysisOn())
  {
    fhDeltaPhiChargedOtherBC  = new TH2F
    ("hDeltaPhiChargedOtherBC","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, track BC!=0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedOtherBC->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhDeltaPhiChargedPtA3GeVOtherBC  = new TH2F
    ("hDeltaPhiChargedPtA3GeVOtherBC","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}, track BC!=0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPtA3GeVOtherBC->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedPtA3GeVOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtTrigChargedOtherBC  =
    new TH2F("hPtTrigChargedOtherBC","trigger and charged tracks pt distribution, track BC!=0",
             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPtTrigChargedOtherBC->SetYTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    fhPtTrigChargedOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEChargedOtherBC  =
    new TH2F("hXEChargedOtherBC","#it{x}_{#it{E}} for charged tracks, track BC!=0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEChargedOtherBC->SetYTitle("#it{x}_{#it{E}}");
    fhXEChargedOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEUeChargedOtherBC  =
    new TH2F("hXEUeChargedOtherBC","#it{x}_{#it{E}} for Underlying Event, track BC!=0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeChargedOtherBC->SetYTitle("#it{x}_{#it{E}}");
    fhXEUeChargedOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhZTChargedOtherBC  =
    new TH2F("hZTChargedOtherBC","#it{z}_{T} for charged tracks, track BC!=0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTChargedOtherBC->SetYTitle("#it{z}_{T}");
    fhZTChargedOtherBC->SetXTitle("#it{p}_{T trigger}");
    
    fhZTUeChargedOtherBC  =
    new TH2F("hZTUeChargedOtherBC","#it{z}_{T} for Underlying Event, track BC!=0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTUeChargedOtherBC->SetYTitle("#it{z}_{T}");
    fhZTUeChargedOtherBC->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    outputContainer->Add(fhDeltaPhiChargedOtherBC) ;
    outputContainer->Add(fhDeltaPhiChargedPtA3GeVOtherBC) ;
    outputContainer->Add(fhXEChargedOtherBC) ;
    outputContainer->Add(fhXEUeChargedOtherBC) ;
    outputContainer->Add(fhZTChargedOtherBC) ;
    outputContainer->Add(fhZTUeChargedOtherBC) ;
    outputContainer->Add(fhPtTrigChargedOtherBC) ;
    
    fhDeltaPhiChargedBC0  = new TH2F
    ("hDeltaPhiChargedBC0","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, track BC==0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedBC0->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhDeltaPhiChargedPtA3GeVBC0  = new TH2F
    ("hDeltaPhiChargedPtA3GeVBC0","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}, track BC==0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPtA3GeVBC0->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedPtA3GeVBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtTrigChargedBC0  =
    new TH2F("hPtTrigChargedBC0","trigger and charged tracks pt distribution, track BC==0",
             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPtTrigChargedBC0->SetYTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    fhPtTrigChargedBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEChargedBC0  =
    new TH2F("hXEChargedBC0","#it{x}_{#it{E}} for charged tracks, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEChargedBC0->SetYTitle("#it{x}_{#it{E}}");
    fhXEChargedBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEUeChargedBC0  =
    new TH2F("hXEUeChargedBC0","#it{x}_{#it{E}} for Underlying Event, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeChargedBC0->SetYTitle("#it{x}_{#it{E}}");
    fhXEUeChargedBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhZTChargedBC0  =
    new TH2F("hZTChargedBC0","#it{z}_{T} for charged tracks, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTChargedBC0->SetYTitle("#it{z}_{T}");
    fhZTChargedBC0->SetXTitle("#it{p}_{T trigger}");
    
    fhZTUeChargedBC0  =
    new TH2F("hZTUeChargedBC0","#it{z}_{T} for Underlying Event, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTUeChargedBC0->SetYTitle("#it{z}_{T}");
    fhZTUeChargedBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    outputContainer->Add(fhDeltaPhiChargedBC0) ;
    outputContainer->Add(fhDeltaPhiChargedPtA3GeVBC0) ;
    outputContainer->Add(fhXEChargedBC0) ;
    outputContainer->Add(fhXEUeChargedBC0) ;
    outputContainer->Add(fhZTChargedBC0) ;
    outputContainer->Add(fhZTUeChargedBC0) ;
    outputContainer->Add(fhPtTrigChargedBC0) ;
    
    fhPtTriggerVtxBC0  = new TH1F("hPtTriggerVtxBC0","#it{p}_{T} distribution of trigger particles", nptbins,ptmin,ptmax);
    fhPtTriggerVtxBC0->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    
    fhDeltaPhiChargedVtxBC0  = new TH2F
    ("hDeltaPhiChargedVtxBC0","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, track BC==0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedVtxBC0->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhDeltaPhiChargedPtA3GeVVtxBC0  = new TH2F
    ("hDeltaPhiChargedPtA3GeVVtxBC0","#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}, track BC==0",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiChargedPtA3GeVVtxBC0->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiChargedPtA3GeVVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtTrigChargedVtxBC0  =
    new TH2F("hPtTrigChargedVtxBC0","trigger and charged tracks pt distribution, track BC==0",
             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPtTrigChargedVtxBC0->SetYTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
    fhPtTrigChargedVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEChargedVtxBC0  =
    new TH2F("hXEChargedVtxBC0","#it{x}_{#it{E}} for charged tracks, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEChargedVtxBC0->SetYTitle("#it{x}_{#it{E}}");
    fhXEChargedVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXEUeChargedVtxBC0  =
    new TH2F("hXEUeChargedVtxBC0","#it{x}_{#it{E}} for Underlying Event, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeChargedVtxBC0->SetYTitle("#it{x}_{#it{E}}");
    fhXEUeChargedVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhZTChargedVtxBC0  =
    new TH2F("hZTChargedVtxBC0","#it{z}_{T} for charged tracks, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTChargedVtxBC0->SetYTitle("#it{z}_{T}");
    fhZTChargedVtxBC0->SetXTitle("#it{p}_{T trigger}");
    
    fhZTUeChargedVtxBC0  =
    new TH2F("hZTUeChargedVtxBC0","#it{z}_{T} for Underlying Event, track BC==0",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTUeChargedVtxBC0->SetYTitle("#it{z}_{T}");
    fhZTUeChargedVtxBC0->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    outputContainer->Add(fhPtTriggerVtxBC0);
    outputContainer->Add(fhDeltaPhiChargedVtxBC0) ;
    outputContainer->Add(fhDeltaPhiChargedPtA3GeVVtxBC0) ;
    outputContainer->Add(fhXEChargedVtxBC0) ;
    outputContainer->Add(fhXEUeChargedVtxBC0) ;
    outputContainer->Add(fhZTChargedVtxBC0) ;
    outputContainer->Add(fhZTUeChargedVtxBC0) ;
    outputContainer->Add(fhPtTrigChargedVtxBC0) ;
    
    for(Int_t i = 0 ; i < 7 ; i++)
    {
      fhPtTriggerPileUp[i]  = new TH1F(Form("hPtTriggerPileUp%s",pileUpName[i].Data()),
                                       Form("#it{p}_{T} distribution of trigger particles, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptmax);
      fhPtTriggerPileUp[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      outputContainer->Add(fhPtTriggerPileUp[i]);
      
      fhDeltaPhiChargedPileUp[i]  = new TH2F(Form("hDeltaPhiChargedPileUp%s",pileUpName[i].Data()),
                                             Form("#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, %s Pile-Up event",pileUpName[i].Data()),
                                             nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedPileUp[i]->SetYTitle("#Delta #varphi (rad)");
      fhDeltaPhiChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhDeltaPhiChargedPileUp[i]) ;
      
      fhDeltaPhiChargedPtA3GeVPileUp[i]  = new TH2F(Form("hDeltaPhiChargedPtA3GeVPileUp%s",pileUpName[i].Data()),
                                                    Form("#varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}, %s Pile-Up event",pileUpName[i].Data()),
                                                    nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedPtA3GeVPileUp[i]->SetYTitle("#Delta #varphi (rad)");
      fhDeltaPhiChargedPtA3GeVPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhDeltaPhiChargedPtA3GeVPileUp[i]) ;
      
      fhDeltaEtaChargedPileUp[i]  = new TH2F(Form("hDeltaEtaChargedPileUp%s",pileUpName[i].Data()),
                                             Form("#eta_{trigger} - #eta_{h^{#pm}} vs #it{p}_{T trigger}, %s Pile-Up event",pileUpName[i].Data()),
                                             nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);
      fhDeltaEtaChargedPileUp[i]->SetYTitle("#Delta #eta");
      fhDeltaEtaChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhDeltaEtaChargedPileUp[i]) ;
      
      fhDeltaEtaChargedPtA3GeVPileUp[i]  = new TH2F(Form("hDeltaEtaChargedPtA3GeVPileUp%s",pileUpName[i].Data()),
                                                    Form("#eta_{trigger} - #eta_{h^{#pm}} vs #it{p}_{T trigger}, #it{p}_{TA}>3 GeV/#it{c}, %s Pile-Up event",pileUpName[i].Data()),
                                                    nptbins,ptmin,ptmax,ndeltaetabins,deltaetamin,deltaetamax);
      fhDeltaEtaChargedPtA3GeVPileUp[i]->SetYTitle("#Delta #eta");
      fhDeltaEtaChargedPtA3GeVPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhDeltaEtaChargedPtA3GeVPileUp[i]) ;
      
      fhXEChargedPileUp[i]  = new TH2F(Form("hXEChargedPileUp%s",pileUpName[i].Data()),
                                       Form("#it{x}_{#it{E}} for charged tracks, %s Pile-Up event",pileUpName[i].Data()),
                                       nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEChargedPileUp[i]->SetYTitle("#it{x}_{#it{E}}");
      fhXEChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEChargedPileUp[i]) ;
      
      fhXEUeChargedPileUp[i]  = new TH2F(Form("hXEUeChargedPileUp%s",pileUpName[i].Data()),
                                         Form("#it{x}_{#it{E}} for Underlying Event, %s Pile-Up event",pileUpName[i].Data()),
                                         nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEUeChargedPileUp[i]->SetYTitle("#it{x}_{#it{E}}");
      fhXEUeChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEUeChargedPileUp[i]) ;
      
      fhZTChargedPileUp[i]  = new TH2F(Form("hZTChargedPileUp%s",pileUpName[i].Data()),
                                       Form("#it{z}_{T} for charged tracks, %s Pile-Up event",pileUpName[i].Data()),
                                       nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTChargedPileUp[i]->SetYTitle("#it{z}_{T}");
      fhZTChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhZTChargedPileUp[i]) ;
      
      fhZTUeChargedPileUp[i]  = new TH2F(Form("hZTUeChargedPileUp%s",pileUpName[i].Data()),
                                         Form("#it{z}_{T} for Underlying Event, %s Pile-Up event",pileUpName[i].Data()),
                                         nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTUeChargedPileUp[i]->SetYTitle("#it{z}_{T}");
      fhZTUeChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhZTUeChargedPileUp[i]) ;
      
      fhPtTrigChargedPileUp[i]  = new TH2F(Form("hPtTrigChargedPileUp%s",pileUpName[i].Data()),
                                           Form("trigger and charged tracks pt distribution, %s Pile-Up event",pileUpName[i].Data()),
                                           nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhPtTrigChargedPileUp[i]->SetYTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
      fhPtTrigChargedPileUp[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhPtTrigChargedPileUp[i]) ;
    }
  }
  
  if(IsHighMultiplicityAnalysisOn())
  {
    Int_t nMultiBins = GetNCentrBin();
    fhDeltaPhiChargedMult = new TH2F*[nMultiBins] ;
    fhDeltaEtaChargedMult = new TH2F*[nMultiBins] ;
    fhXEMult              = new TH2F*[nMultiBins] ;
    fhXEUeMult            = new TH2F*[nMultiBins] ;
    fhZTMult              = new TH2F*[nMultiBins] ;
    fhZTUeMult            = new TH2F*[nMultiBins] ;
    
    for(Int_t im=0; im<nMultiBins; im++)
    {
      fhDeltaPhiChargedMult[im]  = new TH2F
      (Form("hDeltaPhiCharged_Mult%d",im),Form("#Delta #varphi charged Mult bin %d",im), nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiChargedMult[im]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhDeltaPhiChargedMult[im]->SetYTitle("#Delta #varphi (rad)");
      
      fhDeltaEtaChargedMult[im]  = new TH2F
      (Form("hDeltaEtaCharged_Mult%d",im),Form("#Delta #eta charged Mult bin %d",im), nptbins,ptmin,ptmax, ndeltaetabins ,deltaetamin,deltaetamax);
      fhDeltaEtaChargedMult[im]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhDeltaEtaChargedMult[im]->SetYTitle("#Delta #eta");
      
      fhXEMult[im]  = new TH2F
      (Form("hXECharged_Mult%d",im),Form("#it{x}_{E} charged Mult bin %d",im), nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEMult[im]->SetYTitle("#it{x}_{E}");
      fhXEMult[im]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhXEUeMult[im]  = new TH2F
      (Form("hXEUeCharged_Mult%d",im),Form("#it{x}_{E} UE charged Mult bin %d",im), nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEUeMult[im]->SetYTitle("#it{x}_{E}");
      fhXEUeMult[im]->SetXTitle("#it{p}_{T trigger}(GeV/#it{c})");
      
      fhZTMult[im]  = new TH2F
      (Form("hZTCharged_Mult%d",im),Form("#it{z}_{T} charged  Mult bin %d",im), nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTMult[im]->SetYTitle("#it{z}_{T}");
      fhZTMult[im]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhZTUeMult[im]  = new TH2F
      (Form("hZTUeCharged_Mult%d",im),Form("#it{z}_{T} UE charged  Mult bin %d",im), nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTUeMult[im]->SetYTitle("#it{z}_{T}");
      fhZTUeMult[im]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      outputContainer->Add(fhDeltaPhiChargedMult[im]) ;
      outputContainer->Add(fhDeltaEtaChargedMult[im]) ;
      outputContainer->Add(fhXEMult  [im]);
      outputContainer->Add(fhXEUeMult[im]);
      outputContainer->Add(fhZTMult  [im]);
      outputContainer->Add(fhZTUeMult[im]);
    }
  }
  
  if(fFillBradHisto)
  {
    fhAssocPtBkg        = new TH2F("hAssocPtBkg", " Trigger #it{p}_{T} vs associated hadron #it{p}_{T} from background",
                                   nptbins, ptmin, ptmax,nptbins,ptmin,ptmax);
    fhAssocPtBkg->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    fhAssocPtBkg->SetYTitle("#it{p}_{T associated} (GeV/#it{c})");
    outputContainer->Add(fhAssocPtBkg) ;
    
    fhDeltaPhiBrad = new TH2F("hDeltaPhiBrad","atan2(sin(#Delta #varphi), cos(#Delta #varphi))/#pi vs #it{p}_{T trigger} ",
                              nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
    fhDeltaPhiBrad->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    fhDeltaPhiBrad->SetYTitle("atan2(sin(#Delta #varphi), cos(#Delta #varphi))/#pi");
    outputContainer->Add(fhDeltaPhiBrad) ;
  }
  
  fhDeltaPhiDeltaEtaAssocPtBin = new TH2F*[fNAssocPtBins*nz];
  fhDeltaPhiAssocPtBin       = new TH2F*[fNAssocPtBins*nz];
  fhDeltaPhiAssocPtBinDEta08 = new TH2F*[fNAssocPtBins*nz];
  fhDeltaPhiAssocPtBinDEta0  = new TH2F*[fNAssocPtBins*nz];
  if(fFillMomImbalancePtAssocBinsHisto)
  {
    fhXEAssocPtBin           = new TH2F*[fNAssocPtBins*nz];
    fhZTAssocPtBin           = new TH2F*[fNAssocPtBins*nz];
  }
  
  if(fCorrelVzBin)
  {
    fhXEVZ = new TH2F*[nz];
    fhZTVZ = new TH2F*[nz];
  }
  
  if(fFillBradHisto)
    fhDeltaPhiBradAssocPtBin = new TH2F*[fNAssocPtBins*nz];
  
  
  fhDeltaPhiAssocPtBin       = new TH2F*[fNAssocPtBins*nz];
  if(fFillEtaGapsHisto)fhDeltaPhiAssocPtBinDEta08       = new TH2F*[fNAssocPtBins*nz];
  if(fDecayTrigger)    fhDeltaPhiDecayChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
  
  
  if(fHMPIDCorrelation)
  {
    fhDeltaPhiAssocPtBinHMPID   = new TH2F*[fNAssocPtBins*nz];
    fhDeltaPhiAssocPtBinHMPIDAcc= new TH2F*[fNAssocPtBins*nz];
  }
  
  for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
  {
    for(Int_t z = 0 ; z < nz ; z++)
    {
      Int_t bin = i*nz+z;
      
      if(fCorrelVzBin)
      {
        sz = Form("_vz%d",z);
        tz = Form(", #it{v}_{#it{z}} bin %d",z);
      }
      
      //printf("iAssoc %d, Vz %d, bin %d - sz %s, tz %s	\n",i,z,bin,sz.Data(),tz.Data());
      
      fhDeltaPhiDeltaEtaAssocPtBin[bin]  = new TH2F(Form("hDeltaPhiDeltaEtaPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                    Form("#Delta #varphi vs #Delta #eta for associated #it{p}_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                    ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins,deltaetamin,deltaetamax);
      fhDeltaPhiDeltaEtaAssocPtBin[bin]->SetXTitle("#Delta #varphi (rad)");
      fhDeltaPhiDeltaEtaAssocPtBin[bin]->SetYTitle("#Delta #eta");
      
      fhDeltaPhiAssocPtBin[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                           Form("#Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                           nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiAssocPtBin[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhDeltaPhiAssocPtBin[bin]->SetYTitle("#Delta #varphi (rad)");
      
      outputContainer->Add(fhDeltaPhiDeltaEtaAssocPtBin[bin]) ;
      outputContainer->Add(fhDeltaPhiAssocPtBin[bin]) ;
      
      if(fFillEtaGapsHisto)
      {
        fhDeltaPhiAssocPtBinDEta08[bin] = new TH2F(Form("hDeltaPhiDeltaEta0.8PtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                   Form("#Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta > 0.8", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                   nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinDEta08[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhDeltaPhiAssocPtBinDEta08[bin]->SetYTitle("#Delta #varphi (rad)");
        
        fhDeltaPhiAssocPtBinDEta0[bin] = new TH2F(Form("hDeltaPhiDeltaEta0PtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                  Form("#Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta = 0.", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                  nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinDEta0[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhDeltaPhiAssocPtBinDEta0[bin]->SetYTitle("#Delta #varphi (rad)");
        
        outputContainer->Add(fhDeltaPhiAssocPtBinDEta08[bin]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBinDEta0[bin]) ;
      }
      
      if(fDecayTrigger)
      {
        fhDeltaPhiDecayChargedAssocPtBin[bin] = new TH2F(Form("hDeltaPhiPtDecayChargedAssocPt%2.1f_%2.1f%s_bit%d", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data(),fDecayBits[0]),
                                                         Form("#Delta #varphi vs #it{p}_{T trigger} tagged as decay for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, Bit %d", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data(),fDecayBits[0]),
                                                         nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiDecayChargedAssocPtBin[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhDeltaPhiDecayChargedAssocPtBin[bin]->SetYTitle("#Delta #varphi (rad)");
        
        outputContainer->Add(fhDeltaPhiDecayChargedAssocPtBin[bin]) ;
      }
      
      if(fFillBradHisto)
      {
        fhDeltaPhiBradAssocPtBin[bin] = new TH2F(Form("hDeltaPhiBradPtAssocPt%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                 Form("atan2(sin(#Delta #varphi), cos(#Delta #varphi))/#pi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                 nptbins, ptmin, ptmax,288, -1.0/3.0, 5.0/3.0);
        fhDeltaPhiBradAssocPtBin[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhDeltaPhiBradAssocPtBin[bin]->SetYTitle("atan2(sin(#Delta #varphi), cos(#Delta #varphi))/#pi");
        outputContainer->Add(fhDeltaPhiBradAssocPtBin[bin]) ;
      }
      
      if(fHMPIDCorrelation)
      {
        fhDeltaPhiAssocPtBinHMPID[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%sHMPID", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                  Form("#Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, with track having HMPID signal", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                  nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinHMPID[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})" );
        fhDeltaPhiAssocPtBinHMPID[bin]->SetYTitle("#Delta #varphi (rad)");
        
        fhDeltaPhiAssocPtBinHMPIDAcc[bin] = new TH2F(Form("hDeltaPhiPtAssocPt%2.1f_%2.1f%sHMPIDAcc", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                     Form("#Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, with track within 5<phi<20 deg", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                     nptbins, ptmin, ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
        fhDeltaPhiAssocPtBinHMPIDAcc[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhDeltaPhiAssocPtBinHMPIDAcc[bin]->SetYTitle("#Delta #varphi (rad)");
        
        outputContainer->Add(fhDeltaPhiAssocPtBinHMPID   [bin]) ;
        outputContainer->Add(fhDeltaPhiAssocPtBinHMPIDAcc[bin]) ;
      }
    }
  }
  
  if(fFillMomImbalancePtAssocBinsHisto)
  {
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      fhXEAssocPtBin[i]       = new TH2F(Form("hXEAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]),
                                         Form("#it{x}_{#it{E}} vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]),
                                         nptbins, ptmin, ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEAssocPtBin[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhXEAssocPtBin[i]->SetYTitle("#it{x}_{#it{E}}");
      
      fhZTAssocPtBin[i]       = new TH2F(Form("hZTAssocPtBin%1.f_%1.f", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]),
                                         Form("#it{z}_{T} vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1]),
                                         nptbins, ptmin, ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTAssocPtBin[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhZTAssocPtBin[i]->SetYTitle("#it{z}_{T}");
      
      
      outputContainer->Add(fhXEAssocPtBin[i]);
      outputContainer->Add(fhZTAssocPtBin[i]);
    }
  }
  
  if(fCorrelVzBin)
  {
    for(Int_t z = 0 ; z < nz ; z++)
    {
      sz = Form("_vz%d",z);
      tz = Form(", #it{v}_{#it{z}} bin %d",z);
      
      fhXEVZ[z]       = new TH2F(Form("hXE%s", sz.Data()),
                                 Form("#it{x}_{#it{E}} vs #it{p}_{T trigger}%s", tz.Data()),
                                 nptbins, ptmin, ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEVZ[z]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhXEVZ[z]->SetYTitle("#it{x}_{#it{E}}");
      
      fhZTVZ[z]       = new TH2F(Form("hZT%s",sz.Data()),
                                 Form("#it{z}_{T} vs #it{p}_{T trigger}%s", tz.Data()),
                                 nptbins, ptmin, ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTVZ[z]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      fhZTVZ[z]->SetYTitle("#it{z}_{T}");
      
      outputContainer->Add(fhXEVZ[z]);
      outputContainer->Add(fhZTVZ[z]);
    }
  }
  
  if(fPi0Trigger)
  {
    fhPtPi0DecayRatio  = new TH2F
    ("hPtPi0DecayRatio","#it{p}_{T} of #pi^{0} and the ratio of pt for two decay",
     nptbins,ptmin,ptmax, 100,0.,2.);
    fhPtPi0DecayRatio->SetXTitle("#it{p}_{T}^{#pi^{0}} (GeV/#it{c})");
    fhPtPi0DecayRatio->SetYTitle("#it{p}_{T}^{Decay}/#it{p}_{T}^{#pi^{0}}");
    outputContainer->Add(fhPtPi0DecayRatio) ;
    
    fhDeltaPhiPi0DecayCharged  = new TH2F
    ("hDeltaPhiPi0DecayCharged","#varphi_{Decay} - #varphi_{h^{#pm}} vs #it{p}_{T Decay}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiPi0DecayCharged->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiPi0DecayCharged->SetXTitle("#it{p}_{T Decay} (GeV/#it{c})");
    
    fhXEPi0DecayCharged  =
    new TH2F("hXEPi0DecayCharged","#it{x}_{#it{E}}  Decay",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEPi0DecayCharged->SetYTitle("#it{x}_{#it{E}}");
    fhXEPi0DecayCharged->SetXTitle("#it{p}_{T decay} (GeV/#it{c})");
    
    fhZTPi0DecayCharged  =
    new TH2F("hZTPi0DecayCharged","#it{z}_{trigger h^{#pm}} = #it{p}_{T h^{#pm}} / #it{p}_{T Decay}",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTPi0DecayCharged->SetYTitle("#it{z}_{decay h^{#pm}}");
    fhZTPi0DecayCharged->SetXTitle("#it{p}_{T decay} (GeV/#it{c})");
    
    outputContainer->Add(fhDeltaPhiPi0DecayCharged) ;
    outputContainer->Add(fhXEPi0DecayCharged) ;
    outputContainer->Add(fhZTPi0DecayCharged) ;
  }
  
  if(fDecayTrigger)
  {
    for(Int_t ibit = 0; ibit< fNDecayBits; ibit++)
    {
      fhDeltaPhiDecayCharged[ibit]  = new TH2F
      (Form("hDeltaPhiDecayCharged_bit%d",fDecayBits[ibit]),
       Form("#varphi_{Decay} - #varphi_{h^{#pm}} vs #it{p}_{T Decay}, Bit %d",fDecayBits[ibit]),
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiDecayCharged[ibit]->SetYTitle("#Delta #varphi (rad)");
      fhDeltaPhiDecayCharged[ibit]->SetXTitle("#it{p}_{T Decay} (GeV/#it{c})");
      
      fhXEDecayCharged[ibit]  =
      new TH2F(Form("hXEDecayCharged_bit%d",fDecayBits[ibit]),
               Form("#it{x}_{#it{E}}  Decay, Bit %d",fDecayBits[ibit]),
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEDecayCharged[ibit]->SetYTitle("#it{x}_{#it{E}}");
      fhXEDecayCharged[ibit]->SetXTitle("#it{p}_{T decay} (GeV/#it{c})");
      
      fhZTDecayCharged[ibit]  =
      new TH2F(Form("hZTDecayCharged_bit%d",fDecayBits[ibit]),
               Form("#it{z}_{trigger h^{#pm}} = #it{p}_{T h^{#pm}} / #it{p}_{T Decay}, Bit %d",fDecayBits[ibit]),
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTDecayCharged[ibit]->SetYTitle("#it{z}_{decay h^{#pm}}");
      fhZTDecayCharged[ibit]->SetXTitle("#it{p}_{T decay} (GeV/#it{c})");
      
      outputContainer->Add(fhDeltaPhiDecayCharged[ibit]) ;
      outputContainer->Add(fhXEDecayCharged[ibit]) ;
      outputContainer->Add(fhZTDecayCharged[ibit]) ;
    }
  }
  
  //Correlation with neutral hadrons
  if(fNeutralCorr)
  {
    fhDeltaPhiDeltaEtaNeutral  = new TH2F
    ("hDeltaPhiDeltaEtaNeutral","#varphi_{trigger} - #varphi_{h^{0}} vs #eta_{trigger} - #eta_{h^{0}}",
     ndeltaphibins ,deltaphimin,deltaphimax, ndeltaetabins ,deltaetamin,deltaetamax);
    fhDeltaPhiDeltaEtaNeutral->SetXTitle("#Delta #varphi (rad)");
    fhDeltaPhiDeltaEtaNeutral->SetYTitle("#Delta #eta");
    
    fhPhiNeutral  = new TH2F
    ("hPhiNeutral","#varphi_{#pi^{0}}  vs #it{p}_{T #pi^{0}}",
     nptbins,ptmin,ptmax,180,0,TMath::TwoPi());
    fhPhiNeutral->SetYTitle("#varphi_{#pi^{0}} (rad)");
    fhPhiNeutral->SetXTitle("#it{p}_{T #pi^{0}} (GeV/#it{c})");
    
    fhEtaNeutral  = new TH2F
    ("hEtaNeutral","#eta_{#pi^{0}}  vs #it{p}_{T #pi^{0}}",
     nptbins,ptmin,ptmax,200,-1.,1.);
    fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
    fhEtaNeutral->SetXTitle("#it{p}_{T #pi^{0}} (GeV/#it{c})");
    
    fhDeltaPhiNeutral  = new TH2F
    ("hDeltaPhiNeutral","#varphi_{trigger} - #varphi_{#pi^{0}} vs #it{p}_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiNeutral->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhDeltaPhiNeutralPt  = new TH2F
    ("hDeltaPhiNeutralPt","#varphi_{trigger} - #varphi_{#pi^{0}} vs #it{p}_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiNeutralPt->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiNeutralPt->SetXTitle("#it{p}_{T h^{0}} (GeV/#it{c})");
    
    fhDeltaEtaNeutral  = new TH2F
    ("hDeltaEtaNeutral","#eta_{trigger} - #eta_{#pi^{0}} vs #it{p}_{T trigger}",
     nptbins,ptmin,ptmax, ndeltaetabins ,deltaetamin,deltaetamax);
    fhDeltaEtaNeutral->SetYTitle("#Delta #eta");
    fhDeltaEtaNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhXENeutral  =
    new TH2F("hXENeutral","#it{x}_{#it{E}} for #pi^{0} associated",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXENeutral->SetYTitle("#it{x}_{#it{E}}");
    fhXENeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtHbpXENeutral  =
    new TH2F("hHbpXENeutral","#xi = ln(1/#it{x}_{#it{E}})for #pi^{0} associated",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpXENeutral->SetYTitle("ln(1/#it{x}_{#it{E}})");
    fhPtHbpXENeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhZTNeutral  =
    new TH2F("hZTNeutral","#it{z}_{trigger #pi} = #it{p}_{T #pi^{0}} / #it{p}_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTNeutral->SetYTitle("#it{z}_{trigger #pi^{0}}");
    fhZTNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtHbpZTNeutral  =
    new TH2F("hHbpZTNeutral","#xi = ln(1/#it{x}_{#it{E}}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpZTNeutral->SetYTitle("ln(1/#it{z}_{T})");
    fhPtHbpZTNeutral->SetXTitle("#it{p}_{T trigger}");
    
    outputContainer->Add(fhDeltaPhiDeltaEtaNeutral);
    outputContainer->Add(fhPhiNeutral) ;
    outputContainer->Add(fhEtaNeutral) ;
    outputContainer->Add(fhDeltaPhiNeutral) ;
    outputContainer->Add(fhDeltaPhiNeutralPt) ;
    outputContainer->Add(fhDeltaEtaNeutral) ;
    outputContainer->Add(fhXENeutral) ;
    outputContainer->Add(fhPtHbpXENeutral) ;
    outputContainer->Add(fhZTNeutral) ;
    outputContainer->Add(fhPtHbpZTNeutral) ;
    
    fhDeltaPhiUeNeutralPt  = new TH2F
    (Form("hDeltaPhiUe%sNeutralPt",right.Data()),"#varphi_{trigger} - #varphi_{#pi^{0}} vs #it{p}_{T #pi^{0}}}",
     nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
    fhDeltaPhiUeNeutralPt->SetYTitle("#Delta #varphi (rad)");
    fhDeltaPhiUeNeutralPt->SetXTitle("#it{p}_{T h^{0}} (GeV/#it{c})");
    
    fhXEUeNeutral  =
    new TH2F(Form("hXEUeNeutral%s",right.Data()),"#it{x}_{#it{E}} for #pi^{0} associated",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhXEUeNeutral->SetYTitle("#it{x}_{#it{E}}");
    fhXEUeNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtHbpXEUeNeutral  =
    new TH2F(Form("hHbpXEUeNeutral%s",right.Data()),"#xi = ln(1/#it{x}_{#it{E}}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/#it{x}_{#it{E}})");
    fhPtHbpXEUeNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhZTUeNeutral  =
    new TH2F(Form("hZTUeNeutral%s",right.Data()),"#it{z}_{trigger #pi} = #it{p}_{T #pi^{0}} / #it{p}_{T trigger} for #pi^{0} associated",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhZTUeNeutral->SetYTitle("#it{z}_{trigger #pi^{0}}");
    fhZTUeNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    fhPtHbpZTUeNeutral  =
    new TH2F(Form("hHbpZTUeNeutral%s",right.Data()),"#xi = ln(1/#it{x}_{#it{E}}) for #pi^{0} associated",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhPtHbpXEUeNeutral->SetYTitle("ln(1/#it{z}_{T})");
    fhPtHbpXEUeNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    
    outputContainer->Add(fhDeltaPhiUeNeutralPt) ;
    outputContainer->Add(fhXEUeNeutral) ;
    outputContainer->Add(fhPtHbpXEUeNeutral) ;
    outputContainer->Add(fhZTUeNeutral) ;
    outputContainer->Add(fhPtHbpZTUeNeutral) ;
    
    if(fMakeSeveralUE)
    {
      fhDeltaPhiUeLeftNeutral  = new TH2F
      ("hDeltaPhiUeLeftNeutralPt","#varphi_{trigger} - #varphi_{#Ueh^{0}} vs #it{p}_{T h^{0}} with neutral UE left side range of trigger particles",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiUeLeftNeutral->SetYTitle("#Delta #varphi (rad)");
      fhDeltaPhiUeLeftNeutral->SetXTitle("#it{p}_{T h^{0}} (GeV/#it{c})");
      outputContainer->Add(fhDeltaPhiUeLeftNeutral) ;
      
      fhXEUeLeftNeutral  =
      new TH2F("hXEUeNeutralLeft","#it{x}_{#it{E}} = #it{p}_{T Ueh^{0}} / #it{p}_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEUeLeftNeutral->SetYTitle("#it{z}_{trigger Ueh^{0}}");
      fhXEUeLeftNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhXEUeLeftNeutral) ;
      
      fhPtHbpXEUeLeftNeutral  =
      new TH2F("hHbpXEUeNeutralLeft","#xi = ln(1/#it{x}_{#it{E}}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhPtHbpXEUeLeftNeutral->SetYTitle("ln(1/#it{x}_{#it{E}})");
      fhPtHbpXEUeLeftNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhPtHbpXEUeLeftNeutral) ;
      
      fhZTUeLeftNeutral  =
      new TH2F("hZTUeNeutralLeft","#it{z}_{trigger h^{0}} = #it{p}_{T Ueh^{0}} / #it{p}_{T trigger} with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTUeLeftNeutral->SetYTitle("#it{z}_{trigger Ueh^{0}}");
      fhZTUeLeftNeutral->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      outputContainer->Add(fhZTUeLeftNeutral) ;
      
      fhPtHbpZTUeLeftNeutral  =
      new TH2F("hHbpZTUeNeutralLeft","#xi = ln(1/#it{z}_{T}) with neutral UE left side of trigger",
               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhPtHbpZTUeLeftNeutral->SetYTitle("ln(1/#it{z}_{T})");
      fhPtHbpZTUeLeftNeutral->SetXTitle("#it{p}_{T trigger}");
      outputContainer->Add(fhPtHbpZTUeLeftNeutral) ;
    }
    
    if(fPi0Trigger)
    {
      fhDeltaPhiPi0DecayNeutral  = new TH2F
      ("hDeltaPhiPi0DecayNeutral","#varphi_{Decay} - #varphi_{h^{0}} vs #it{p}_{T Decay}",
       nptbins,ptmin,ptmax, ndeltaphibins ,deltaphimin,deltaphimax);
      fhDeltaPhiPi0DecayNeutral->SetYTitle("#Delta #varphi (rad)");
      fhDeltaPhiPi0DecayNeutral->SetXTitle("#it{p}_{T Decay} (GeV/#it{c})");
      
      fhXEPi0DecayNeutral  =
      new TH2F("hXEPi0DecayNeutral","#it{x}_{#it{E}} for decay trigger",
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhXEPi0DecayNeutral->SetYTitle("#it{x}_{#it{E}}");
      fhXEPi0DecayNeutral->SetXTitle("#it{p}_{T decay}");
      
      fhZTPi0DecayNeutral  =
      new TH2F("hZTPi0DecayNeutral","#it{z}_{trigger h^{0}} = #it{p}_{T h^{0}} / #it{p}_{T Decay}",
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhZTPi0DecayNeutral->SetYTitle("#it{z}_{h^{0}}");
      fhZTPi0DecayNeutral->SetXTitle("#it{p}_{T decay}");
      
      outputContainer->Add(fhDeltaPhiPi0DecayNeutral) ;
      outputContainer->Add(fhXEPi0DecayNeutral) ;
      outputContainer->Add(fhZTPi0DecayNeutral) ;
    }
  }//Correlation with neutral hadrons
  
  // If data is MC, fill more histograms, depending on origin
  if(IsDataMC())
  {
    for(Int_t i= fMCGenTypeMin; i <= fMCGenTypeMax; i++)
    {
      fhMCPtTrigger[i]  = new TH1F (Form("hMCPtTrigger_%s",nameMC[i].Data()),
                                    Form("MC %s: generated trigger #it{p}_{T}",nameMC[i].Data()),
                                    nptbins,ptmin,ptmax);
      fhMCPtTrigger[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      
      fhMCPhiTrigger[i]  = new TH2F (Form("hMCPhiTrigger_%s",nameMC[i].Data()),
                                     Form("MC %s: generated trigger #varphi",nameMC[i].Data()),
                                     nptbins,ptmin,ptmax, nphibins,phimin,phimax);
      fhMCPhiTrigger[i]->SetYTitle("#varphi (rad)");
      fhMCPhiTrigger[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      
      fhMCEtaTrigger[i]  = new TH2F (Form("hMCEtaTrigger_%s",nameMC[i].Data()),
                                     Form("MC %s: generated trigger #eta",nameMC[i].Data()),
                                     nptbins,ptmin,ptmax, netabins,etamin,etamax);
      fhMCEtaTrigger[i]->SetYTitle("#eta");
      fhMCEtaTrigger[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      
      if(fMakeAbsoluteLeading || fMakeNearSideLeading)
      {
        fhMCPtTriggerNotLeading[i]  = new TH1F (Form("hMCPtTriggerNotLeading_%s",nameMC[i].Data()),
                                                Form("MC %s: generated trigger #it{p}_{T}, when not leading of primaries",nameMC[i].Data()),
                                                nptbins,ptmin,ptmax);
        fhMCPtTriggerNotLeading[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        
        fhMCPhiTriggerNotLeading[i]  = new TH2F (Form("hMCPhiTriggerNotLeading_%s",nameMC[i].Data()),
                                                 Form("MC %s: generated trigger #varphi, when not leading of primaries",nameMC[i].Data()),
                                                 nptbins,ptmin,ptmax, nphibins,phimin,phimax);
        fhMCPhiTriggerNotLeading[i]->SetYTitle("#varphi (rad)");
        fhMCPhiTriggerNotLeading[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        
        
        fhMCEtaTriggerNotLeading[i]  = new TH2F (Form("hMCEtaTriggerNotLeading_%s",nameMC[i].Data()),
                                                 Form("MC %s: generated triogger #eta, when not leading of primaries",nameMC[i].Data()),
                                                 nptbins,ptmin,ptmax, netabins,etamin,etamax);
        fhMCEtaTriggerNotLeading[i]->SetYTitle("#eta ");
        fhMCEtaTriggerNotLeading[i]->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      }
      
      fhMCEtaCharged[i]  = new TH2F (Form("hMCEtaCharged_%s",nameMC[i].Data()),
                                     Form("MC %s: #eta_{h^{#pm}}  vs #it{p}_{T #pm}",nameMC[i].Data()),
                                     nptbins,ptmin,ptmax,100,-1.,1.);
      fhMCEtaCharged[i]->SetYTitle("#eta_{h^{#pm}} (rad)");
      fhMCEtaCharged[i]->SetXTitle("#it{p}_{T #pm} (GeV/#it{c})");
      
      fhMCPhiCharged[i]  = new TH2F(Form("hMCPhiCharged_%s",nameMC[i].Data()),
                                    Form("MC %s: phi_{h^{#pm}}  vs #it{p}_{T #pm}",nameMC[i].Data()),
                                    nptbins,ptmin,ptmax,180,0,TMath::TwoPi());
      fhMCPhiCharged[i]->SetYTitle("MC #varphi_{h^{#pm}} (rad)");
      fhMCPhiCharged[i]->SetXTitle("#it{p}_{T #pm} (GeV/#it{c})");
      
      fhMCDeltaPhiDeltaEtaCharged[i]  = new TH2F (Form("hMCDeltaPhiDeltaEtaCharged_%s",nameMC[i].Data()),
                                                  Form("MC %s: phi_{trigger} - #varphi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",nameMC[i].Data()),
                                                  ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax);
      fhMCDeltaPhiDeltaEtaCharged[i]->SetXTitle("#Delta #varphi (rad)");
      fhMCDeltaPhiDeltaEtaCharged[i]->SetYTitle("#Delta #eta");
      
      fhMCDeltaEtaCharged[i]  = new TH2F (Form("hMCDeltaEtaCharged_%s",nameMC[i].Data()),
                                          Form("MC %s: #eta_{trigger} - #eta_{h^{#pm}} vs #it{p}_{T trigger} and #it{p}_{T assoc}",nameMC[i].Data()),
                                          nptbins,ptmin,ptmax,ndeltaetabins ,deltaetamin,deltaetamax);
      fhMCDeltaEtaCharged[i]->SetYTitle("#Delta #eta");
      fhMCDeltaEtaCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCDeltaPhiCharged[i]  = new TH2F  (Form("hMCDeltaPhiCharged_%s",nameMC[i].Data()),
                                           Form("MC %s: #varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}",nameMC[i].Data()),
                                           nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
      fhMCDeltaPhiCharged[i]->SetYTitle("#Delta #varphi (rad)");
      fhMCDeltaPhiCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCDeltaPhiChargedPt[i]  = new TH2F (Form("hMCDeltaPhiChargedPt_%s",nameMC[i].Data()),
                                            Form("MC %s: #varphi_{trigger} - #varphi_{#h^{#pm}} vs #it{p}_{T h^{#pm}}",nameMC[i].Data()),
                                            nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
      fhMCDeltaPhiChargedPt[i]->SetYTitle("#Delta #varphi (rad)");
      fhMCDeltaPhiChargedPt[i]->SetXTitle("#it{p}_{T h^{#pm}} (GeV/#it{c})");
      
      fhMCPtXECharged[i]  = new TH2F (Form("hMCPtXECharged_%s",nameMC[i].Data()),
                                      Form("MC %s: #it{x}_{#it{E}} with charged hadrons",nameMC[i].Data()),
                                      nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhMCPtXECharged[i]->SetYTitle("#it{x}_{#it{E}}");
      fhMCPtXECharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtHbpXECharged[i]  = new TH2F(Form("hMCHbpXECharged_%s",nameMC[i].Data()),
                                        Form("MC %s: #xi = ln(1/#it{x}_{#it{E}}) with charged hadrons",nameMC[i].Data()),
                                        nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhMCPtHbpXECharged[i]->SetYTitle("ln(1/#it{x}_{#it{E}})");
      fhMCPtHbpXECharged[i]->SetXTitle("#it{p}_{T trigger}");
      
      fhMCPtZTCharged[i]  = new TH2F(Form("hMCPtZTCharged_%s",nameMC[i].Data()),
                                     Form("MC %s: #it{z}_{T} with charged hadrons",nameMC[i].Data()),
                                     nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhMCPtZTCharged[i]->SetYTitle("#it{z}_{T}");
      fhMCPtZTCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtHbpZTCharged[i]  = new TH2F(Form("hMCHbpZTCharged_%s",nameMC[i].Data()),
                                        Form("MC %s: #xi = ln(1/#it{z}_{T}) with charged hadrons",nameMC[i].Data()),
                                        nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhMCPtHbpZTCharged[i]->SetYTitle("ln(1/#it{z}_{T})");
      fhMCPtHbpZTCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtTrigPout[i]  = new TH2F(Form("hMCPtTrigPout_%s",nameMC[i].Data()),
                                    Form("MC %s: #it{p}_{out} with triggers",nameMC[i].Data()),
                                    nptbins,ptmin,ptmax,nptbins,-1.*ptmax/2.,ptmax/2.);
      fhMCPtTrigPout[i]->SetYTitle("#it{p}_{out} (GeV/#it{c})");
      fhMCPtTrigPout[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtAssocDeltaPhi[i]  = new TH2F(Form("hMCPtAssocDeltaPhi_%s",nameMC[i].Data()),
                                         Form("MC %s: #Delta #varphi with associated charged hadrons",nameMC[i].Data()),
                                         nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
      fhMCPtAssocDeltaPhi[i]->SetYTitle("#Delta #varphi (rad)");
      fhMCPtAssocDeltaPhi[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      outputContainer->Add(fhMCPtTrigger[i]);
      outputContainer->Add(fhMCPhiTrigger[i]);
      outputContainer->Add(fhMCEtaTrigger[i]);
      
      if(fMakeAbsoluteLeading || fMakeNearSideLeading)
      {
        outputContainer->Add(fhMCPtTriggerNotLeading[i]);
        outputContainer->Add(fhMCPhiTriggerNotLeading[i]);
        outputContainer->Add(fhMCEtaTriggerNotLeading[i]);
      }
      
      outputContainer->Add(fhMCDeltaPhiDeltaEtaCharged[i]);
      outputContainer->Add(fhMCPhiCharged[i]) ;
      outputContainer->Add(fhMCEtaCharged[i]) ;
      outputContainer->Add(fhMCDeltaEtaCharged[i]) ;
      outputContainer->Add(fhMCDeltaPhiCharged[i]) ;
      
      outputContainer->Add(fhMCDeltaPhiChargedPt[i]) ;
      outputContainer->Add(fhMCPtXECharged[i]) ;
      outputContainer->Add(fhMCPtZTCharged[i]) ;
      outputContainer->Add(fhMCPtHbpXECharged[i]) ;
      outputContainer->Add(fhMCPtHbpZTCharged[i]) ;
      outputContainer->Add(fhMCPtTrigPout[i]) ;
      outputContainer->Add(fhMCPtAssocDeltaPhi[i]) ;
      
      // Underlying event
      
      fhMCUePart[i]  =
      new TH1F(Form("hMCUePart_%s",nameMC[i].Data()),
               Form("MC %s: UE particles distribution vs #it{p}_{T trigger}",nameMC[i].Data()),
               nptbins,ptmin,ptmax);
      fhMCUePart[i]->SetYTitle("#it{dN}^{ch}");
      fhMCUePart[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtXEUeCharged[i]  =
      new TH2F(Form("hMCPtXEUeCharged%s_%s",right.Data(),nameMC[i].Data()),
               Form("MC %s: #it{x}_{#it{E}} with charged hadrons, Underlying Event %s",nameMC[i].Data(),right.Data()),
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhMCPtXEUeCharged[i]->SetYTitle("#it{x}_{#it{E}}");
      fhMCPtXEUeCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtHbpXEUeCharged[i] =
      new TH2F(Form("hMCPtHbpXEUeCharged%s_%s",right.Data(),nameMC[i].Data()),
               Form("MC %s: #xi = ln(1/#it{x}_{#it{E}}) with charged hadrons, Underlying Event %s",nameMC[i].Data(),right.Data()),
               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhMCPtHbpXEUeCharged[i]->SetYTitle("ln(1/#it{x}_{#it{E}})");
      fhMCPtHbpXEUeCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtZTUeCharged[i] =
      new TH2F(Form("hMCPtZTUeCharged%s_%s",right.Data(),nameMC[i].Data()),
               Form("MC %s: #it{z}_{T} with charged hadrons, Underlying Event %s",nameMC[i].Data(),right.Data()),
               nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
      fhMCPtZTUeCharged[i]->SetYTitle("#it{z}_{T}");
      fhMCPtZTUeCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      fhMCPtHbpZTUeCharged[i] =
      new TH2F(Form("hMCPtHbpZTUeCharged%s_%s",right.Data(),nameMC[i].Data()),
               Form("MC %s: #xi = ln(1/#it{z}_{T}) with charged hadrons, Underlying Event %s",nameMC[i].Data(),right.Data()),
               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
      fhMCPtHbpZTUeCharged[i]->SetYTitle("ln(1/#it{z}_{T})");
      fhMCPtHbpZTUeCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
      
      outputContainer->Add(fhMCUePart[i]);
      outputContainer->Add(fhMCPtXEUeCharged[i]) ;
      outputContainer->Add(fhMCPtZTUeCharged[i]) ;
      outputContainer->Add(fhMCPtHbpZTUeCharged[i]);
      outputContainer->Add(fhMCPtHbpXEUeCharged[i]);
      
      if(fMakeSeveralUE)
      {
        fhMCPtXEUeLeftCharged[i]  = new TH2F(Form("hMCPtXEUeChargedLeft_%s",nameMC[i].Data()),
                                             Form("MC %s: #it{x}_{#it{E}} with charged hadrons, with UE left side range of trigger particles",nameMC[i].Data()),
                                             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
        fhMCPtXEUeLeftCharged[i]->SetYTitle("#it{x}_{#it{E}}");
        fhMCPtXEUeLeftCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        
        fhMCPtHbpXEUeLeftCharged[i] = new TH2F(Form("hMCPtHbpXEUeChargedLeft_%s",nameMC[i].Data()),
                                               Form("MC %s: #xi = ln(1/#it{x}_{#it{E}}) with charged hadrons, with UE left side range of trigger particles",nameMC[i].Data()),
                                               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
        fhMCPtHbpXEUeLeftCharged[i]->SetYTitle("ln(1/#it{x}_{#it{E}})");
        fhMCPtHbpXEUeLeftCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        
        fhMCPtZTUeLeftCharged[i]  = new TH2F(Form("hMCPtZTUeChargedLeft_%s",nameMC[i].Data()),
                                             Form("MC %s: #it{z}_{T} with charged hadrons, with UE left side range of trigger particles",nameMC[i].Data()),
                                             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
        fhMCPtZTUeLeftCharged[i]->SetYTitle("#it{z}_{T}");
        fhMCPtZTUeLeftCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        
        fhMCPtHbpZTUeLeftCharged[i] = new TH2F(Form("hMCPtHbpZTUeChargedLeft_%s",nameMC[i].Data()),
                                               Form("MC %s: #xi = ln(1/#it{z}_{T}) with charged hadrons, with UE left side range of trigger particles",nameMC[i].Data()),
                                               nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
        fhMCPtHbpZTUeLeftCharged[i]->SetYTitle("ln(1/#it{z}_{T})");
        fhMCPtHbpZTUeLeftCharged[i]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        
        outputContainer->Add(fhMCPtXEUeLeftCharged[i]) ;
        outputContainer->Add(fhMCPtZTUeLeftCharged[i]) ;
        outputContainer->Add(fhMCPtHbpXEUeLeftCharged[i]);
        outputContainer->Add(fhMCPtHbpZTUeLeftCharged[i]) ;
        
      }
    }
  } // For MC histogram
  
  if(DoOwnMix())
  {
    // Create event containers
    
    if(!fUseMixStoredInReader || (fUseMixStoredInReader && !GetReader()->ListWithMixedEventsForTracksExists()))
    {
      Int_t nvz = GetNZvertBin();
      Int_t nrp = GetNRPBin();
      Int_t nce = GetNCentrBin();
      
      fListMixTrackEvents= new TList*[nvz*nrp*nce] ;
      
      for( Int_t ice = 0 ; ice < nce ; ice++ )
      {
        for( Int_t ivz = 0 ; ivz < nvz ; ivz++ )
        {
          for( Int_t irp = 0 ; irp < nrp ; irp++ )
          {
            Int_t bin = GetEventMixBin(ice,ivz,irp); //ic*nvz*nrp+iz*nrp+irp;
            
            //printf("GetCreateOutputObjects - Bins : cent %d, vz %d, RP %d, event %d/%d\n",
            //       ic,iz, irp, bin);
            
            fListMixTrackEvents[bin] = new TList() ;
            fListMixTrackEvents[bin]->SetOwner(kFALSE);
          }
        }
      }
    }
    
    fhPtTriggerMixed  = new TH1F ("hPtTriggerMixed","#it{p}_{T} distribution of trigger particles, used for mixing", nptbins,ptmin,ptmax);
    fhPtTriggerMixed->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    
    if(fCorrelVzBin)
    {
      fhPtTriggerMixedVzBin  = new TH2F ("hPtTriggerMixedVzBin","#it{p}_{T} distribution of trigger particles, used for mixing", nptbins,ptmin,ptmax,GetNZvertBin(),0,GetNZvertBin());
      fhPtTriggerMixedVzBin->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
      fhPtTriggerMixedVzBin->SetYTitle("#it{v}_{#it{z}} bin");
      outputContainer->Add(fhPtTriggerMixedVzBin);
    }
    
    fhPtTriggerMixedBin  = new TH2F ("hPtTriggerMixedBin","#it{p}_{T} distribution of trigger particles vs mixing bin", nptbins,ptmin,ptmax,nMixBins,0,nMixBins);
    fhPtTriggerMixedBin->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPtTriggerMixedBin->SetYTitle("Bin");
    
    fhPhiTriggerMixed  = new TH2F ("hPhiTriggerMixed","#varphi distribution of trigger Particles, used for mixing",nptbins,ptmin,ptmax, nphibins,phimin,phimax);
    fhPhiTriggerMixed->SetYTitle("#varphi (rad)");
    
    fhEtaTriggerMixed  = new TH2F ("hEtaTriggerMixed","#eta distribution of trigger, used for mixing",nptbins,ptmin,ptmax, netabins,etamin,etamax);
    fhEtaTriggerMixed->SetYTitle("#eta ");
    
    outputContainer->Add(fhPtTriggerMixed);
    outputContainer->Add(fhPtTriggerMixedBin);
    outputContainer->Add(fhPhiTriggerMixed);
    outputContainer->Add(fhEtaTriggerMixed);
    
    // Fill the cluster pool only in isolation analysis or if requested
    if( neutralMix && (!fUseMixStoredInReader || (fUseMixStoredInReader && !GetReader()->ListWithMixedEventsForCaloExists())))
    {
      Int_t nvz = GetNZvertBin();
      Int_t nrp = GetNRPBin();
      Int_t nce = GetNCentrBin();
      
      fListMixCaloEvents= new TList*[nvz*nrp*nce] ;
      
      for( Int_t ice = 0 ; ice < nce ; ice++ )
      {
        for( Int_t ivz = 0 ; ivz < nvz ; ivz++ )
        {
          for( Int_t irp = 0 ; irp < nrp ; irp++ )
          {
            Int_t bin = GetEventMixBin(ice,ivz,irp); //ic*nvz*nrp+iz*nrp+irp;
            
            //printf("GetCreateOutputObjects - Bins : cent %d, vz %d, RP %d, event %d/%d\n",
            //       ic,iz, irp, bin);
            
            fListMixCaloEvents[bin] = new TList() ;
            fListMixCaloEvents[bin]->SetOwner(kFALSE);
          }
        }
      }
    }
    
    // Init the list in the reader if not done previously
    if(fUseMixStoredInReader)
    {
      if( !GetReader()->ListWithMixedEventsForTracksExists() )
        GetReader()->SetListWithMixedEventsForTracks(fListMixTrackEvents);
      
      if( !GetReader()->ListWithMixedEventsForCaloExists()   )
        GetReader()->SetListWithMixedEventsForCalo  (fListMixCaloEvents );
    }
    
    fhEventBin=new TH1I("hEventBin","Number of triggers per bin(cen,vz,rp)",
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                        GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventBin->SetXTitle("event bin");
    outputContainer->Add(fhEventBin) ;
    
    fhEventMixBin=new TH1I("hEventMixBin","Number of triggers mixed per event bin(cen,vz,rp)",
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                           GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventMixBin->SetXTitle("event bin");
    outputContainer->Add(fhEventMixBin) ;
    
    fhEventMBBin=new TH1I("hEventMBBin","Number of min bias events per bin(cen,vz,rp)",
                          GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                          GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhEventMBBin->SetXTitle("event bin");
    outputContainer->Add(fhEventMBBin) ;
    
    fhNtracksMB=new TH2F("hNtracksMBEvent","Number of filtered tracks in MB event per event bin",ntrbins,trmin,trmax,
                         GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                         GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
    fhNtracksMB->SetYTitle("event bin");
    fhNtracksMB->SetXTitle("#it{N}_{track}");
    outputContainer->Add(fhNtracksMB);
    
    if( neutralMix )
    {
      fhNclustersMB=new TH2F("hNclustersMBEvent","Number of filtered clusters in MB events per event bin",nclbins,clmin,clmax,
                             GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1,0,
                             GetNCentrBin()*GetNZvertBin()*GetNRPBin()+1) ;
      fhNclustersMB->SetYTitle("event bin");
      fhNclustersMB->SetXTitle("#it{N}_{cluster}");
      outputContainer->Add(fhNclustersMB);
    }
    
    fhMixDeltaPhiCharged  = new TH2F
    ("hMixDeltaPhiCharged","Mixed event : #varphi_{trigger} - #varphi_{h^{#pm}} vs #it{p}_{T trigger}",
     nptbins,ptmin,ptmax,ndeltaphibins ,deltaphimin,deltaphimax);
    fhMixDeltaPhiCharged->SetYTitle("#Delta #varphi (rad)");
    fhMixDeltaPhiCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhMixDeltaPhiCharged);
    
    fhMixDeltaPhiDeltaEtaCharged  = new TH2F
    ("hMixDeltaPhiDeltaEtaCharged","Mixed event : #varphi_{trigger} - #varphi_{h^{#pm}} vs #eta_{trigger} - #eta_{h^{#pm}}",
     ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax);
    fhMixDeltaPhiDeltaEtaCharged->SetXTitle("#Delta #varphi (rad)");
    fhMixDeltaPhiDeltaEtaCharged->SetYTitle("#Delta #eta");
    outputContainer->Add(fhMixDeltaPhiDeltaEtaCharged);
    
    fhMixXECharged  =
    new TH2F("hMixXECharged","Mixed event : #it{x}_{#it{E}} for charged tracks",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhMixXECharged->SetYTitle("#it{x}_{#it{E}}");
    fhMixXECharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhMixXECharged);
    
    fhMixXEUeCharged  =
    new TH2F("hMixXEUeCharged","Mixed event : #it{x}_{#it{E}} for charged tracks in Ue region",
             nptbins,ptmin,ptmax,nxeztbins,xeztmin,xeztmax);
    fhMixXEUeCharged->SetYTitle("#it{x}_{#it{E}}");
    fhMixXEUeCharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhMixXEUeCharged);
    
    fhMixHbpXECharged  =
    new TH2F("hMixHbpXECharged","mixed event : #xi = ln(1/#it{x}_{#it{E}}) with charged hadrons",
             nptbins,ptmin,ptmax,nhbpbins,hbpmin,hbpmax);
    fhMixHbpXECharged->SetYTitle("ln(1/#it{x}_{#it{E}})");
    fhMixHbpXECharged->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
    outputContainer->Add(fhMixHbpXECharged);
    
    fhMixDeltaPhiChargedAssocPtBin         = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiChargedAssocPtBinDEta08   = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiChargedAssocPtBinDEta0    = new TH2F*[fNAssocPtBins*nz];
    fhMixDeltaPhiDeltaEtaChargedAssocPtBin = new TH2F*[fNAssocPtBins*nz];
    
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      for(Int_t z = 0 ; z < nz ; z++)
      {
        Int_t bin = i*nz+z;
        
        if(fCorrelVzBin)
        {
          sz = Form("_vz%d",z);
          tz = Form(", #it{v}_{#it{z}} bin %d",z);
        }
        
        //printf("MIX : iAssoc %d, Vz %d, bin %d - sz %s, tz %s	\n",i,z,bin,sz.Data(),tz.Data());
        
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin] = new TH2F(Form("hMixDeltaPhiDeltaEtaChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                               Form("Mixed event #Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                               ndeltaphibins ,deltaphimin,deltaphimax,ndeltaetabins ,deltaetamin,deltaetamax);
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->SetXTitle("#Delta #varphi (rad)");
        fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->SetYTitle("#Delta #eta");
        
        outputContainer->Add(fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]);
        
        fhMixDeltaPhiChargedAssocPtBin[bin] = new TH2F(Form("hMixDeltaPhiChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                       Form("Mixed event #Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                       nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
        fhMixDeltaPhiChargedAssocPtBin[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
        fhMixDeltaPhiChargedAssocPtBin[bin]->SetYTitle("#Delta #varphi (rad)");
        
        outputContainer->Add(fhMixDeltaPhiChargedAssocPtBin[bin]);
        
        if(fFillEtaGapsHisto)
        {
          fhMixDeltaPhiChargedAssocPtBinDEta08[bin] = new TH2F(Form("hMixDeltaPhiDeltaEta0.8ChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                               Form("Mixed event #Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta > 0.8", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                               nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
          fhMixDeltaPhiChargedAssocPtBinDEta08[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
          fhMixDeltaPhiChargedAssocPtBinDEta08[bin]->SetYTitle("#Delta #varphi (rad)");
          
          fhMixDeltaPhiChargedAssocPtBinDEta0[bin] = new TH2F(Form("hMixDeltaPhiDeltaEta0ChargedAssocPtBin%2.1f_%2.1f%s", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],sz.Data()),
                                                              Form("Mixed event #Delta #varphi vs #it{p}_{T trigger} for associated #it{p}_{T} bin [%2.1f,%2.1f]%s, for #Delta #eta = 0", fAssocPtBinLimit[i], fAssocPtBinLimit[i+1],tz.Data()),
                                                              nptbins, ptmin, ptmax,  ndeltaphibins ,deltaphimin,deltaphimax);
          fhMixDeltaPhiChargedAssocPtBinDEta0[bin]->SetXTitle("#it{p}_{T trigger} (GeV/#it{c})");
          fhMixDeltaPhiChargedAssocPtBinDEta0[bin]->SetYTitle("#Delta #varphi (rad)");
          
          outputContainer->Add(fhMixDeltaPhiChargedAssocPtBinDEta08[bin]);
          outputContainer->Add(fhMixDeltaPhiChargedAssocPtBinDEta0[bin]);
        }
      }
    }
  }
  
  return outputContainer;
}

//_____________________________________________________________________________________________________________________
/// Get the momentum of the pi0/eta assigned decay photons.
/// In case of pi0/eta trigger, we may want to check their decay correlation,
/// get their decay children.
/// \return kTRUE if decay photon indeces are found.
/// \param indexPhoton1: ID of photon cluster.
/// \param indexPhoton2: ID of photon cluster.
/// \param idetector: ID of detector.
//_____________________________________________________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::GetDecayPhotonMomentum(Int_t indexPhoton1, Int_t indexPhoton2, Int_t idetector)
{
  if(indexPhoton1!=-1 || indexPhoton2!=-1) return kFALSE;
  
  AliDebug(1,Form("indexPhoton1 = %d, indexPhoton2 = %d", indexPhoton1, indexPhoton2));
  
  TObjArray * clusters  = 0x0 ;
  if(idetector==kEMCAL) clusters = GetEMCALClusters() ;
  else                  clusters = GetPHOSClusters()  ;
  
  for(Int_t iclus = 0; iclus < clusters->GetEntriesFast(); iclus++)
  {
    AliVCluster * photon =  (AliVCluster*) (clusters->At(iclus));
    
    if(photon->GetID()==indexPhoton1) photon->GetMomentum(fDecayMom1,GetVertex(0)) ;
    if(photon->GetID()==indexPhoton2) photon->GetMomentum(fDecayMom2,GetVertex(0)) ;
    
    AliDebug(1,Form("Photon1 = %f, Photon2 = %f", fDecayMom1.Pt(), fDecayMom2.Pt()));
  } //cluster loop
  
  return kTRUE;
}

//________________________________________________________________________
/// \return Index of MC histograms depending on MC trigger particle origin.
/// \param mcTag: MC tag (AliCaloPID) of cluster/track triggering
//________________________________________________________________________
Int_t AliAnaParticleHadronCorrelation::GetMCTagHistogramIndex(Int_t mcTag)
{
  if     ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt) ||
          GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation) ||
          GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCISR)          ) return 0;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)          ) return 1;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay)     ) return 2;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)          ) return 3;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay)     ) return 4;
  else if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)       ) return 5; // other decays
  else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)     ) return 6;
  else                                                                                    return 7;
}

//_________________________________________
/// Init. If tracks are not loaded, abort.
//_________________________________________
void AliAnaParticleHadronCorrelation::Init()
{
  if(!GetReader()->IsCTSSwitchedOn())
    AliFatal("STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!");
}

//____________________________________________________
/// Initialize the default parameters of the analysis.
//____________________________________________________
void AliAnaParticleHadronCorrelation::InitParameters()
{
  SetInputAODName("Particle");
  SetAODObjArrayName("Hadrons");
  AddToHistogramsName("AnaHadronCorr_");
  
  SetPtCutRange(0.,300);
  fDeltaPhiMinCut       = TMath::DegToRad()*120.;
  fDeltaPhiMaxCut       = TMath::DegToRad()*240. ;
  fSelectIsolated       = kFALSE;
  fMakeSeveralUE        = kFALSE;
  fUeDeltaPhiMinCut     = TMath::DegToRad()*60.;
  fUeDeltaPhiMaxCut     = TMath::DegToRad()*120 ;
  
  fNeutralCorr          = kFALSE ;
  fPi0Trigger           = kFALSE ;
  fDecayTrigger         = kFALSE ;
  fHMPIDCorrelation     = kFALSE ;
  
  fMakeAbsoluteLeading  = kTRUE;
  fMakeNearSideLeading  = kFALSE;
  
  fNAssocPtBins         = 9   ;
  fAssocPtBinLimit[0]   = 0.2 ;
  fAssocPtBinLimit[1]   = 0.5 ;
  fAssocPtBinLimit[2]   = 1.0 ;
  fAssocPtBinLimit[3]   = 2.0 ;
  fAssocPtBinLimit[4]   = 3.0 ;
  fAssocPtBinLimit[5]   = 4.0 ;
  fAssocPtBinLimit[6]   = 5.0 ;
  fAssocPtBinLimit[7]   = 6.0 ;
  fAssocPtBinLimit[8]   = 7.0 ;
  fAssocPtBinLimit[9]   = 8.0 ;
  fAssocPtBinLimit[10]  = 9.0 ;
  fAssocPtBinLimit[11]  = 10.0 ;
  fAssocPtBinLimit[12]  = 12.0 ;
  fAssocPtBinLimit[13]  = 14.0 ;
  fAssocPtBinLimit[14]  = 16.0 ;
  fAssocPtBinLimit[15]  = 20.0 ;
  fAssocPtBinLimit[16]  = 30.0 ;
  fAssocPtBinLimit[17]  = 40.0 ;
  fAssocPtBinLimit[18]  = 50.0 ;
  fAssocPtBinLimit[19]  = 200.0 ;
  
  fUseMixStoredInReader = kTRUE;
  
  fM02MinCut   = -1 ;
  fM02MaxCut   = -1 ;
  
  fDecayTagsM02Cut = 0.27;
  
  fSelectLeadingHadronAngle = kFALSE;
  fFillLeadHadOppositeHisto = kFALSE;
  fMinLeadHadPhi = 150*TMath::DegToRad();
  fMaxLeadHadPhi = 210*TMath::DegToRad();
  
  fMinLeadHadPt  = 1;
  fMaxLeadHadPt  = 100;
  
  fMCGenTypeMin =  0;
  fMCGenTypeMax = fgkNmcTypes-1;
  
  fNDecayBits = 1;
  fDecayBits[0] = AliNeutralMesonSelection::kPi0;
  fDecayBits[1] = AliNeutralMesonSelection::kEta;
  fDecayBits[2] = AliNeutralMesonSelection::kPi0RightSide;
  fDecayBits[3] = AliNeutralMesonSelection::kEtaRightSide;
  fDecayBits[4] = AliNeutralMesonSelection::kEtaLeftSide ;
  fDecayBits[5] = AliNeutralMesonSelection::kEtaBothSides;
  fDecayBits[6] = AliNeutralMesonSelection::kPi0LeftSide ; // Leave it last since likely not used
  fDecayBits[7] = AliNeutralMesonSelection::kPi0BothSides; // Leave it last since likely not used
  
  fNBkgBin = 11;
  fBkgBinLimit[ 0] = 00.0; fBkgBinLimit[ 1] = 00.2; fBkgBinLimit[ 2] = 00.3; fBkgBinLimit[ 3] = 00.4; fBkgBinLimit[ 4] = 00.5;
  fBkgBinLimit[ 5] = 01.0; fBkgBinLimit[ 6] = 01.5; fBkgBinLimit[ 7] = 02.0; fBkgBinLimit[ 8] = 03.0; fBkgBinLimit[ 9] = 05.0;
  fBkgBinLimit[10] = 10.0; fBkgBinLimit[11] = 100.;
  for(Int_t ibin = fNBkgBin+1; ibin < 20; ibin++)
  {
    fBkgBinLimit[ibin] = 00.0;
  }
  
}

//________________________________________________________________________________________________________
/// Do the invariant mass of the trigger particle with other clusters.
/// Activate for photon analysis.
/// \param trigger: Trigger cluster/track kinematics and more, AliAODPWG4ParticleCorrelation object.
/// \param mcIndex: Histograms MC index of trigger particle
//________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::InvMassHisto(AliAODPWG4ParticleCorrelation * trigger, Int_t mcIndex)
{
  Int_t    idTrig = trigger->GetCaloLabel(0);
  Float_t  ptTrig = trigger->Pt();
  Float_t tofTrig = trigger->GetTime();
  
  fMomentum = *(trigger->Momentum());
  
  // Loop on second cluster to be associated to trigger particle to have the invariant mass
  Int_t nphoton = GetInputAODBranch()->GetEntriesFast();
  
  for(Int_t iphoton = 0; iphoton < nphoton; iphoton++)
  {
    AliAODPWG4ParticleCorrelation * photon1 =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iphoton));
    
    if(idTrig == photon1->GetCaloLabel(0)) continue;        //Do not the invariant mass for the same particle
    
    fMomentumIM = *(photon1->Momentum());
    
    // Select secondary photon with proper invariant mass
    if(fM02MaxCut > 0 && fM02MinCut > 0) //clID1 > 0 && clID2 < 0 &&
    {
      Float_t m02 = photon1->GetM02();
      
      if(m02 > fM02MaxCut || m02 < fM02MinCut) continue ;
    }
    
    // Select clusters with good time window difference
    Double_t tdiff = tofTrig - photon1->GetTime();
    if(TMath::Abs(tdiff) > GetPairTimeCut()) continue;
    
    // Add momenta of trigger and possible secondary decay
    fMomentumIM += fMomentum;
    
    Double_t mass  = fMomentumIM.M();
    //    Double_t epair = fMomentumIM.E();
    //    Float_t ptpair = fMomentumIM.Pt();
    
    // Mass of all pairs
    fhMassPtTrigger->Fill(ptTrig, mass, GetEventWeight());
    if ( IsDataMC() ) fhMCMassPtTrigger[mcIndex]->Fill(ptTrig, mass, GetEventWeight());
  }
}

//_________________________________________________________________________
/// Check if the what of the selected triggers is leading particle comparing
/// with all the triggers, all the tracks or all the clusters (if requested for the clusters).
//_________________________________________________________________________
Bool_t AliAnaParticleHadronCorrelation::IsTriggerTheEventLeadingParticle()
{
  Double_t ptTrig      = GetMinPt();
  Double_t phiTrig     = 0 ;
  fLeadingTriggerIndex =-1 ;
  Int_t index          =-1 ;
  AliAODPWG4ParticleCorrelation* pLeading = 0;
  
  // Loop on stored AOD particles, find leading trigger on the selected list, with at least min pT.
  
  for(Int_t iaod = 0; iaod < GetInputAODBranch()->GetEntriesFast() ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    particle->SetLeadingParticle(kFALSE); // set it later
    
    // Vertex cut in case of mixing
    Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
    if(check ==  0) continue;
    if(check == -1) return kFALSE; // not sure if it is correct.
    
    // Find the leading particles with highest momentum
    if (particle->Pt() > ptTrig)
    {
      ptTrig   = particle->Pt() ;
      phiTrig  = particle->Phi();
      index    = iaod     ;
      pLeading = particle ;
    }
  }// Finish search of leading trigger particle on the AOD branch.
  
  if(index < 0) return kFALSE;
  
  //printf("AOD leading pT %2.2f, ID %d\n", pLeading->Pt(),pLeading->GetCaloLabel(0));
  
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  // Compare if it is the leading of all tracks
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
    // do not count the candidate or the daughters of the candidate
    // in the isolation conte
    if ( pLeading->GetDetectorTag() == kCTS ) // make sure conversions are tagged as kCTS!!!
    {
      Int_t  trackID   = GetReader()->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
      Bool_t contained = kFALSE;
      
      for(Int_t i = 0; i < 4; i++) 
      {
        if( trackID == pLeading->GetTrackLabel(i) ) contained = kTRUE;
      }
      
      if ( contained ) continue ;
    }
    
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    Float_t pt   = fTrackVector.Pt();
    Float_t phi  = fTrackVector.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    
    //jump out this event if near side associated particle pt larger than trigger
    if (fMakeNearSideLeading)
    {
      Float_t deltaPhi = phiTrig-phi;
      if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
      if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
      
      if(pt > ptTrig && deltaPhi < TMath::PiOver2())  return kFALSE;
    }
    //jump out this event if there is any other particle with pt larger than trigger
    else
    {
      if(pt > ptTrig)  return kFALSE ;
    }
  }// track loop
  
  // Compare if it is leading of all calorimeter clusters
  
  if(fCheckLeadingWithNeutralClusters)
  {
    // Select the calorimeter cluster list
    TObjArray * nePl = 0x0;
    if      (pLeading->GetDetectorTag() == kPHOS )
      nePl = GetPHOSClusters();
    else
      nePl = GetEMCALClusters();
    
    if(!nePl) return kTRUE; // Do the selection just with the tracks if no calorimeter is available.
    
    for(Int_t ipr = 0;ipr < nePl->GetEntriesFast() ; ipr ++ )
    {
      AliVCluster * cluster = (AliVCluster *) (nePl->At(ipr)) ;
      
      if(cluster->GetID() == pLeading->GetCaloLabel(0) || cluster->GetID() == pLeading->GetCaloLabel(1) ) continue ;
      
      cluster->GetMomentum(fMomentum,GetVertex(0));
      
      Float_t pt   = fMomentum.Pt();
      Float_t phi  = fMomentum.Phi() ;
      if(phi < 0) phi+=TMath::TwoPi();
      
      if(IsTrackMatched(cluster,GetReader()->GetInputEvent())) continue ; // avoid charged clusters, already covered by tracks, or cluster merging with track.
      
      //jump out this event if near side associated particle pt larger than trigger
      // not really needed for calorimeter, unless DCal is included
      if (fMakeNearSideLeading)
      {
        Float_t deltaPhi = phiTrig-phi;
        if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
        if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
        
        if(pt > ptTrig && deltaPhi < TMath::PiOver2()) return kFALSE ;
      }
      //jump out this event if there is any other particle with pt larger than trigger
      else
      {
        if(pt > ptTrig)  return kFALSE ;
      }
    }// cluster loop
  } // check neutral clusters
  
  fLeadingTriggerIndex = index ;
  pLeading->SetLeadingParticle(kTRUE);
  
  AliDebug(1,Form("\t particle AOD with index %d is leading with pT %2.2f", fLeadingTriggerIndex, pLeading->Pt()));
  
  return kTRUE;
}

//_________________________________________________________________
/// Main method, loop on list of trigger particles and fill histograms
/// correlating with tracks. Steps:
/// * If selected, check if the trigger is leading particle (in near hemisphere or all)
/// * Loop on trigger
///    * Select trigger within a shower shape window when requiered.
///    * Selecte isolated triggers (tagged in other task) if requested
///    * Selectr triggers in a fiducial region if requested
///    * Find leading particle in opposite hemisphere and do correlation if found, only if requested.
///    * Make the correlation with charged hadrons
///    * Make correlation with charged at generated level
///    * Make mixed event correlation
///    * Make the correlation with neutral, if requested
///    * Fill some histograms releated to the trigger particle.
//_________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms()
{
  if(!GetInputAODBranch())
  {
    AliFatal(Form("No input particles in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
    return ; // coverity
  }
  
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if( naod == 0 )
  {
    AliDebug(1,"No particle AOD found!");
    return ; // no trigger particles found.
  }
  
  AliDebug(1,Form("Begin hadron correlation analysis, fill histograms"));
  AliDebug(1,Form("n particle branch aod entries %d", naod));
  AliDebug(1,Form("In CTS aod entries %d",GetCTSTracks()->GetEntriesFast()));
  
  //------------------------------------------------------
  // Find leading trigger if analysis request only leading,
  // if there is no leading trigger, then skip the event
  
  Int_t iaod = 0 ;
  if( fMakeAbsoluteLeading || fMakeNearSideLeading )
  {
    Bool_t leading = IsTriggerTheEventLeadingParticle();
    
    AliDebug(1,Form("AOD Leading trigger? %d, with index %d",leading,fLeadingTriggerIndex));
    
    if(!leading)
    {
      AliDebug(1,"Leading was requested and not found");
      return ;
    }
    else
    {
      // Select only the leading in the trigger AOD loop
      naod = fLeadingTriggerIndex+1 ;
      iaod = fLeadingTriggerIndex   ;
    }
  }
  
  //------------------------------------------------------
  // Get event multiplicity and bins
  
  Float_t cen         = GetEventCentrality();
  Float_t ep          = GetEventPlaneAngle();
  if(IsHighMultiplicityAnalysisOn()) fhTriggerEventPlaneCentrality->Fill(cen, ep, GetEventWeight());
  
  Int_t   mixEventBin = GetEventMixBin();
  Int_t   vzbin       = GetEventVzBin();
  
  //------------------------------------------------------
  // Loop on trigger AOD
  
  for( iaod = 0; iaod < naod; iaod++ )
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //
    // Trigger particle selection criteria:
    //
    Float_t pt = particle->Pt();
    
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ;
    
    fhPtTriggerInput->Fill(pt, GetEventWeight());
    
    //
    // check if it was a calorimeter cluster
    // and if the shower shape cut was requested apply it.
    // Not needed if already done at the particle identification level,
    // but for isolation studies, it is preferred not to remove so we do it here
    //
    
    AliDebug(1,Form("%s Trigger : min %f, max %f, det %d",
                    GetInputAODName().Data(),fM02MinCut,fM02MaxCut,particle->GetDetectorTag()));
    
    if(fM02MaxCut > 0 && fM02MinCut > 0) //clID1 > 0 && clID2 < 0 &&
    {
      //      Int_t iclus = -1;
      //      TObjArray* clusters = 0x0;
      //      if     (particle->GetDetectorTag() == kEMCAL) clusters = GetEMCALClusters();
      //      else if(particle->GetDetectorTag() == kPHOS ) clusters = GetPHOSClusters();
      //
      //      if(clusters)
      //      {
      //        AliVCluster *cluster = FindCluster(clusters,clID1,iclus);
      //        Float_t m02 = cluster->GetM02();
      //        if(m02 > fM02MaxCut || m02 < fM02MinCut) continue ;
      //      }
      
      Float_t m02 = particle->GetM02();
      
      if(m02 > fM02MaxCut || m02 < fM02MinCut) continue ;
      
      fhPtTriggerSSCut->Fill(pt, GetEventWeight());
      
      AliDebug(1,"Pass the shower shape cut");
    }
    
    //
    // Check if the particle is isolated or if we want to take the isolation into account
    // This bool is set in AliAnaParticleIsolation
    //
    if(OnlyIsolated())
    {
      if( !particle->IsIsolated() ) continue;
      
      fhPtTriggerIsoCut->Fill(pt, GetEventWeight());
      
      AliDebug(1,"Pass the isolation cut");
    }
    
    //
    // Check if trigger is in fiducial region
    //
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(particle->Eta(),particle->Phi(),particle->GetDetectorTag()) ;
      
      if(! in ) continue ;
      
      AliDebug(1,"Pass the fiducial cut");
    }
    
    fhPtTriggerFidCut->Fill(pt, GetEventWeight());
    
    //---------------------------------------
    // Make correlation
    
    // Find the leading hadron in the opposite hemisphere to the trigger
    // and accept the trigger if leading is in defined window.
    Bool_t okLeadHad = kTRUE;
    if(fSelectLeadingHadronAngle || fFillLeadHadOppositeHisto)
    {
      okLeadHad = FindLeadingOppositeHadronInWindow(particle);
      if(!okLeadHad && fSelectLeadingHadronAngle) continue;
    }
    
    //
    // Charged particles correlation
    //
    MakeChargedCorrelation(particle);
    
    // MC
    Int_t mcIndex = -1;
    Int_t mcTag   = particle->GetTag();
    Bool_t lostDecayPair = kFALSE;
    if(IsDataMC())
    {
      mcIndex = GetMCTagHistogramIndex(mcTag);
      lostDecayPair = GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost);
      MakeMCChargedCorrelation(particle->GetLabel(), mcIndex,lostDecayPair);
    }
    
    // Do own mixed event with charged,
    // add event and remove previous or fill the mixed histograms
    if(DoOwnMix())
      MakeChargedMixCorrelation(particle);
    
    //
    // Neutral particles correlation
    //
    if(fNeutralCorr)
      MakeNeutralCorrelation(particle);
    
    //----------------------------------------------------------------
    // Fill trigger pT related histograms if not absolute leading
    
    //
    // pT of the trigger, vs trigger origin if MC
    //
    fhPtTrigger->Fill(pt, GetEventWeight());
    
    
    if(IsDataMC() && mcIndex >=0 && mcIndex < fgkNmcTypes)
    {
      fhPtTriggerMC[mcIndex]->Fill(pt, GetEventWeight());
      if( lostDecayPair )
      {
        // check index in GetMCTagIndexHistogram
        if     ( mcIndex == 2 ) fhPtTriggerMC[8]->Fill(pt, GetEventWeight()); // pi0 decay
        else if( mcIndex == 4 ) fhPtTriggerMC[9]->Fill(pt, GetEventWeight()); // eta decay
      }
    }
    
    //
    // pT lead cone, pT sum cone background bins distribution
    //
    if(fFillBkgBinsHisto)
    {
      Float_t m02 = particle->GetM02();
      Float_t pTLeadTrackInCone    = 0;
      Float_t pTSumTrackInCone     = 0;
      Float_t pTLeadClusterInCone  = 0;
      Float_t pTSumClusterInCone   = 0;
      
      pTLeadTrackInCone   = particle->GetChargedLeadPtInCone();
      pTLeadClusterInCone = particle->GetNeutralLeadPtInCone();
      
      pTSumTrackInCone    = particle->GetChargedPtSumInCone();
      pTSumClusterInCone  = particle->GetNeutralPtSumInCone();
      
      Float_t pTLeadInCone = pTLeadTrackInCone;
      if(pTLeadClusterInCone > pTLeadInCone) pTLeadInCone = pTLeadClusterInCone;
      Float_t pTSumInCone = pTSumTrackInCone + pTSumClusterInCone;
      
      // Get the background bin for this cone and trigger
      Int_t pTSumBin  = -1;
      Int_t pTLeadBin = -1;
      
      for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
      {
        if( pTSumInCone  >= fBkgBinLimit[ibin] && pTSumInCone  < fBkgBinLimit[ibin+1]) pTSumBin  = ibin;
        if( pTLeadInCone >= fBkgBinLimit[ibin] && pTLeadInCone < fBkgBinLimit[ibin+1]) pTLeadBin = ibin;
      }
      if(pTSumBin > 0)
      {
        fhPtSumInConeBin[pTSumBin]->Fill(pt, GetEventWeight());
      }
      
      if(pTLeadBin > 0)
      {
        fhPtLeadInConeBin[pTLeadBin]->Fill(pt, GetEventWeight());
      }
      
      // Check if it was a decay
      if( fFillTaggedDecayHistograms && m02 < fDecayTagsM02Cut )
      {
        Int_t decayTag = particle->DecayTag();
        if(decayTag < 0) decayTag = 0;
        
        for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
          {
            Int_t pTLeadBinDecay = pTLeadBin+ibit*fNBkgBin;
            Int_t  pTSumBinDecay =  pTSumBin+ibit*fNBkgBin;
            if( pTLeadBin >=0 ) fhPtLeadConeBinDecay[pTLeadBinDecay]->Fill(pt, GetEventWeight());
            if( pTSumBin  >=0 ) fhSumPtConeBinDecay [pTSumBinDecay] ->Fill(pt, GetEventWeight());
          }
        }
      }
      
      if(IsDataMC())
      {
        Int_t pTLeadBinMC = pTLeadBin+mcIndex*fNBkgBin;
        Int_t pTSumBinMC  =  pTSumBin+mcIndex*fNBkgBin;
        
        if( pTLeadBin >=0 )
        {
          fhPtLeadConeBinMC[pTLeadBinMC]->Fill(pt, GetEventWeight());
        }
        
        if( pTSumBin  >=0 )
        {
          fhSumPtConeBinMC [pTSumBinMC]->Fill(pt, GetEventWeight());
        }
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
        {
          pTLeadBinMC = pTLeadBin+kmcPhoton*fNBkgBin;
          pTSumBinMC  =  pTSumBin+kmcPhoton*fNBkgBin;
          if( pTLeadBin >=0 )
          {
            fhPtLeadConeBinMC[pTLeadBinMC]->Fill(pt, GetEventWeight());
          }
          
          if( pTSumBin  >=0 )
          {
            fhSumPtConeBinMC [ pTSumBinMC]->Fill(pt, GetEventWeight());
          }
        }
        
        // Check if decay and if pair is lost
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay )
          {
            pTLeadBinMC = pTLeadBin+kmcPi0DecayLostPair*fNBkgBin;
            pTSumBinMC  =  pTSumBin+kmcPi0DecayLostPair*fNBkgBin;
          }
          else if(mcIndex == kmcEtaDecay)
          {
            pTLeadBinMC = pTLeadBin+kmcEtaDecayLostPair*fNBkgBin;
            pTSumBinMC  =  pTSumBin+kmcEtaDecayLostPair*fNBkgBin;
          }
          else
            AliFatal(Form("Lost decay Bit assigned to bad case, mcIndex %d",mcIndex));
          
          if( pTLeadBin >=0 )
          {
            fhPtLeadConeBinMC[pTLeadBinMC]->Fill(pt, GetEventWeight());
          }
          
          if( pTSumBin  >=0 )
          {
            fhSumPtConeBinMC [ pTSumBinMC]->Fill(pt, GetEventWeight());
          }
          
        }
        
      }
    }
    
    // Invariant mass of trigger particle
    if(fFillInvMassHisto) InvMassHisto(particle,mcIndex);
    
    if(fDecayTrigger)
    {
      Int_t decayTag = particle->DecayTag();
      if(decayTag < 0) decayTag = 0;
      
      for(Int_t ibit = 0; ibit<fNDecayBits; ibit++)
      {
        if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
        {
          fhPtDecayTrigger[ibit]->Fill(pt, GetEventWeight());
          
          if(IsDataMC() && mcIndex >=0 && mcIndex < fgkNmcTypes)
          {
            fhPtDecayTriggerMC[ibit][mcIndex]->Fill(pt, GetEventWeight());
            if( lostDecayPair )
            {
              // check index in GetMCTagIndexHistogram
              if( mcIndex == 2 )  fhPtDecayTriggerMC[ibit][8]->Fill(pt, GetEventWeight()); // pi0 decay
              if( mcIndex == 4 )  fhPtDecayTriggerMC[ibit][9]->Fill(pt, GetEventWeight()); // eta decay
            }
          }
        }// check bit
      }// bit loop
    }
    
    //
    // Acceptance of the trigger
    //
    Float_t phi = particle->Phi();
    if( phi < 0 ) phi+=TMath::TwoPi();
    fhPhiTrigger->Fill(pt, phi, GetEventWeight());
    
    fhEtaTrigger->Fill(pt, particle->Eta(), GetEventWeight());
    //printf("AliAnaParticleHadronCorrelation::MakeAnalysisFillHistograms() - Trigger particle : pt %f, eta %f, phi %f\n",particle->Pt(),particle->Eta(),phi);
    
    //----------------------------------
    // Trigger particle pT vs event bins
    
    fhPtTriggerBin->Fill(pt, mixEventBin, GetEventWeight());
    if(fCorrelVzBin)
      fhPtTriggerVzBin->Fill(pt, vzbin, GetEventWeight());
    
    if(IsHighMultiplicityAnalysisOn())
    {
      fhPtTriggerCentrality->Fill(pt, cen, GetEventWeight());
      fhPtTriggerEventPlane->Fill(pt, ep , GetEventWeight());
    }
    
    //----------------------------------
    // Trigger particle pT vs pile-up
    
    if(IsPileUpAnalysisOn())
    {
      Int_t vtxBC = GetReader()->GetVertexBC();
      if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA)     fhPtTriggerVtxBC0->Fill(pt, GetEventWeight());
      
      if(GetReader()->IsPileUpFromSPD())               fhPtTriggerPileUp[0]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCal())             fhPtTriggerPileUp[1]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtTriggerPileUp[2]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtTriggerPileUp[3]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtTriggerPileUp[4]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtTriggerPileUp[5]->Fill(pt, GetEventWeight());
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtTriggerPileUp[6]->Fill(pt, GetEventWeight());
    }
  } // AOD trigger loop
  
  // Re-init for next event
  fLeadingTriggerIndex = -1;
  
  AliDebug(1,"End fill histograms");
}

//_______________________________________________________________________________________________________
/// Loop on tracks. Fill histograms correlating them with the trigger particle:
///  * Azimuthal correlations
///  * Opposite hemisphere xE, zT, pT ... correlation
///  * Underlying event xE, zT etc on perperndicular cones
///  * Put the correlated charged tracks in a reference array in the trigger particle object if requested.
//_______________________________________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle)
{
  AliDebug(1,"Make trigger particle - charged hadron correlation");
  
  Float_t phiTrig = aodParticle->Phi();
  Float_t etaTrig = aodParticle->Eta();
  Float_t ptTrig  = aodParticle->Pt();
  Int_t    mcTag  = aodParticle->GetTag();
  Double_t bz     = GetReader()->GetInputEvent()->GetMagneticField();
  
  
  Int_t   decayTag = 0;
  if(fDecayTrigger)
  {
    //decay = aodParticle->IsTagged();
    decayTag = aodParticle->DecayTag();
    if(decayTag < 0) decayTag = 0;
    //    printf("Correlation: pT %2.2f, BTag %d, Tagged %d\n",ptTrig, decayTag, aodParticle->IsTagged());
    //    printf("\t check bit Pi0 %d, Eta %d,  Pi0Side %d, EtaSide %d\n",
    //           GetNeutralMesonSelection()->CheckDecayBit(decayTag,AliNeutralMesonSelection::kPi0),
    //           GetNeutralMesonSelection()->CheckDecayBit(decayTag,AliNeutralMesonSelection::kEta),
    //           GetNeutralMesonSelection()->CheckDecayBit(decayTag,AliNeutralMesonSelection::kPi0Side),
    //           GetNeutralMesonSelection()->CheckDecayBit(decayTag,AliNeutralMesonSelection::kEtaSide));
  }
  
  Float_t pt       = -100. ;
  Float_t phi      = -100. ;
  Float_t eta      = -100. ;
  Float_t deltaPhi = -100. ;
  
  TObjArray * reftracks = 0x0;
  Int_t nrefs           = 0;
  
  // Mixed event settings
  Int_t evtIndex11   = -1 ; // cluster trigger or pi0 trigger
  Int_t evtIndex12   = -1 ; // pi0 trigger
  Int_t evtIndex13   = -1 ; // charged trigger
  
  if (GetMixedEvent())
  {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;
    evtIndex13 = GetMixedEvent()->EventIndex(aodParticle->GetTrackLabel(0)) ;
  }
  
  // Track multiplicity or cent bin
  Int_t cenbin = 0;
  if(IsHighMultiplicityAnalysisOn()) cenbin = GetEventCentralityBin();
  
  //
  // In case of pi0/eta trigger, we may want to check their decay correlation,
  // get their decay children
  //
  
  Bool_t decayFound = kFALSE;
  if( fPi0Trigger )
  {
    decayFound = GetDecayPhotonMomentum(aodParticle->GetCaloLabel(0),aodParticle->GetCaloLabel(1),aodParticle->GetDetectorTag());
    if(decayFound)
    {
      fhPtPi0DecayRatio->Fill(ptTrig, fDecayMom1.Pt()/ptTrig, GetEventWeight());
      fhPtPi0DecayRatio->Fill(ptTrig, fDecayMom2.Pt()/ptTrig, GetEventWeight());
    }
  }
  
  //-----------------------------------------------------------------------
  // Track loop, select tracks with good pt, phi and fill AODs or histograms
  //-----------------------------------------------------------------------
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    pt   = fTrackVector.Pt();
    eta  = fTrackVector.Eta();
    phi  = fTrackVector.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    
    //Select only hadrons in pt range
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
    // do not count the candidate or the daughters of the candidate
    // in the isolation conte
    if ( aodParticle->GetDetectorTag() == kCTS ) // make sure conversions are tagged as kCTS!!!
    {
      Int_t  trackID   = GetReader()->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
      Bool_t contained = kFALSE;
      
      for(Int_t i = 0; i < 4; i++) 
      {
        if( trackID == aodParticle->GetTrackLabel(i) ) contained = kTRUE;
      }
      
      if ( contained ) continue ;
    }
    
    //Only for mixed event frame
    Int_t evtIndex2 = 0 ;
    if (GetMixedEvent())
    {
      evtIndex2 = GetMixedEvent()->EventIndex(GetReader()->GetTrackID(track)) ;
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 || evtIndex13 == evtIndex2 ) // photon and track from different events
        continue ;
      //vertex cut
      if (TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut())
        continue;
    }
    
    AliDebug(2,Form("Selected charge for momentum imbalance: pt %2.2f, phi %2.2f, eta %2.2f",pt,phi,eta));
    
    // ------------------------------
    // Track type bin or bits setting
    //
    
    //
    // * Set the pt associated bin for the defined bins *
    //
    Int_t assocBin   = -1;
    for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
    {
      if(pt > fAssocPtBinLimit[i] && pt < fAssocPtBinLimit[i+1]) assocBin= i;
    }
    
    //
    // * Assign to the histogram array a bin corresponding
    // to a combination of pTa and vz bins *
    //
    Int_t nz = 1;
    Int_t vz = 0;
    
    if(fCorrelVzBin)
    {
      nz = GetNZvertBin();
      vz = GetEventVzBin();
    }
    
    Int_t bin = assocBin*nz+vz;
    
    //printf("assoc Bin = %d, vZ bin  = %d, bin = %d \n", assocBin,GetEventVzBin(),bin);
    
    //
    // * Get the status of the TOF bit *
    //
    ULong_t status = track->GetStatus();
    Bool_t okTOF = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
    //Double32_t tof = track->GetTOFsignal()*1e-3;
    Int_t trackBC = track->GetTOFBunchCrossing(bz);
    
    Int_t outTOF = -1;
    if     (okTOF && trackBC!=0) outTOF = 1;
    else if(okTOF && trackBC==0) outTOF = 0;
    
    //----------------
    // Fill Histograms
    
    //
    // Azimuthal Angle histograms
    //
    
    deltaPhi = phiTrig-phi;
    
    //
    // Calculate deltaPhi shift so that for the particles on the opposite side
    // it is defined between 90 and 270 degrees
    // Shift [-360,-90]  to [0, 270]
    // and [270,360] to [-90,0]
    if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
    
    FillChargedAngularCorrelationHistograms(pt,  ptTrig,  bin, phi, phiTrig,  deltaPhi,
                                            eta, etaTrig, decayTag, track->GetHMPIDsignal(),
                                            outTOF, cenbin, mcTag);
    
    //
    // Imbalance zT/xE/pOut histograms
    //
    
    AliESDtrack * esdTrack = dynamic_cast<AliESDtrack*>(track);
    
    //
    // Delta phi cut for momentum imbalance correlation
    //
    if  ( (deltaPhi > fDeltaPhiMinCut)   && (deltaPhi < fDeltaPhiMaxCut)   )
    {
      // Check track resolution
      if ( esdTrack )
      {
        fhTrackResolution->Fill(pt, TMath::Sqrt(esdTrack->GetCovariance()[14])*pt, 
                                GetEventWeight());
      }

      FillChargedMomentumImbalanceHistograms(ptTrig, pt, deltaPhi, cenbin, track->Charge(),
                                             assocBin, decayTag, outTOF, mcTag);
    }
    
    //
    // Underlying event, right side, default case
    //
    if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )
    {
      // Check track resolution
      if ( esdTrack )
      {
        fhTrackResolutionUE->Fill(pt, TMath::Sqrt(esdTrack->GetCovariance()[14])*pt, 
                                  GetEventWeight());
      }
      
      FillChargedUnderlyingEventHistograms(ptTrig, pt, deltaPhi, cenbin, outTOF,mcTag);
    }
    
    //
    // Several UE calculation,  in different perpendicular regions, up to 6:
    // left, right, upper-left, lower left, upper-right, lower-right
    //
    if(fMakeSeveralUE)
      FillChargedUnderlyingEventSidesHistograms(ptTrig,pt,deltaPhi,mcTag);
    
    //
    if(fPi0Trigger && decayFound)
      FillDecayPhotonCorrelationHistograms(pt, phi, kTRUE) ;
    
    //
    // Add track reference to array
    //
    if(fFillAODWithReferences)
    {
      nrefs++;
      if(nrefs==1)
      {
        reftracks = new TObjArray(0);
        TString trackname = Form("%sTracks", GetAODObjArrayName().Data());
        reftracks->SetName(trackname.Data());
        reftracks->SetOwner(kFALSE);
      }
      
      reftracks->Add(track);
    }// reference track to AOD
  }// track loop
  
  //Fill AOD with reference tracks, if not filling histograms
  if(fFillAODWithReferences && reftracks)
  {
    aodParticle->AddObjArray(reftracks);
  }
}

//_________________________________________________________________________________________________________
/// Mix current trigger with tracks in another Minimum Bias event.
//_________________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation(AliAODPWG4ParticleCorrelation *aodParticle)
{
  AliDebug(1,Form("Make trigger particle - charged hadron mixed event correlation"));
  
  if(GetMixedEvent()) return;  // This is not the mixed event from general mixing frame
  
  // Get the event with similar caracteristics
  //printf("MakeChargedMixCorrelation for %s\n",GetInputAODName().Data());
  
  AliAnalysisManager   * manager      = AliAnalysisManager::GetAnalysisManager();
  
  AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  if(!inputHandler) return;
  
  if(!(inputHandler->IsEventSelected( ) & GetReader()->GetEventTriggerMask())) return;
  
  // Get the pool, check if it exits
  Int_t eventBin = GetEventMixBin();
  
  //Check that the bin exists, if not (bad determination of RP, ntrality or vz bin) do nothing
  if(eventBin < 0) return;
  
  fhEventBin->Fill(eventBin, GetEventWeight());
  
  // Get neutral clusters pool?
  Bool_t isoCase = OnlyIsolated() && (GetIsolationCut()->GetParticleTypeInCone() != AliIsolationCut::kOnlyCharged);
  Bool_t neutralMix = fFillNeutralEventMixPool || isoCase ;
  
  TList * pool     = 0;
  TList * poolCalo = 0;
  if(fUseMixStoredInReader)
  {
    pool     = GetReader()->GetListWithMixedEventsForTracks(eventBin);
    if(neutralMix) poolCalo = GetReader()->GetListWithMixedEventsForCalo  (eventBin);
  }
  else
  {
    pool     = fListMixTrackEvents[eventBin];
    if(neutralMix) poolCalo = fListMixCaloEvents [eventBin];
  }
  
  if(!pool) return ;
  
  if( neutralMix && !poolCalo )
    AliWarning("Careful, cluster pool not available");
  
  Double_t ptTrig  = aodParticle->Pt();
  Double_t etaTrig = aodParticle->Eta();
  Double_t phiTrig = aodParticle->Phi();
  if(phiTrig < 0.) phiTrig+=TMath::TwoPi();
  
  AliDebug(1,Form("Pool bin %d size %d, trigger trigger pt=%f, phi=%f, eta=%f",
                  eventBin,pool->GetSize(), ptTrig,phiTrig,etaTrig));
  
  Double_t ptAssoc  = -999.;
  Double_t phiAssoc = -999.;
  Double_t etaAssoc = -999.;
  Double_t deltaPhi = -999.;
  Double_t deltaEta = -999.;
  Double_t xE       = -999.;
  
  // Start from first event in pool except if in this same event the pool was filled
  Int_t ev0 = 0;
  if(GetReader()->GetLastTracksMixedEvent() == GetEventNumber()) ev0 = 1;
  
  for(Int_t ev=ev0; ev < pool->GetSize(); ev++)
  {
    //
    // Recover the lists of tracks or clusters
    //
    TObjArray* bgTracks = static_cast<TObjArray*>(pool->At(ev));
    TObjArray* bgCalo   = 0;
    
    // Recover the clusters list if requested
    if( neutralMix && poolCalo )
    {
      if(pool->GetSize()!=poolCalo->GetSize())
        AliWarning("Different size of calo and track pools");
      
      bgCalo = static_cast<TObjArray*>(poolCalo->At(ev));
      
      if(!bgCalo) AliDebug(1,Form("Event %d in calo pool not available?",ev));
    }
    
    //
    // Isolate the trigger in the mixed event with mixed tracks and clusters
    //
    if( OnlyIsolated() )
    {
      Int_t   n=0, nfrac = 0;
      Bool_t  isolated = kFALSE;
      Float_t coneptsum = 0, coneptlead = 0;
      GetIsolationCut()->MakeIsolationCut(bgTracks,bgCalo,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, aodParticle, "",
                                          n,nfrac,coneptsum,coneptlead,isolated);
      
      //printf("AliAnaParticleHadronCorrelation::MakeChargedMixCorrelation() - Isolated? %d - cone %f, ptthres %f",
      //       isolated,GetIsolationCut()->GetConeSize(),GetIsolationCut()->GetPtThreshold());
      //if(bgTracks)printf(" - n track %d", bgTracks->GetEntriesFast());
      //printf("\n");
      
      if(!isolated) continue ;
    }
    
    //
    // Check if the trigger is leading of mixed event
    //
    Int_t nTracks=bgTracks->GetEntriesFast();
    
    if(fMakeNearSideLeading || fMakeAbsoluteLeading)
    {
      Bool_t leading = kTRUE;
      for(Int_t jlead = 0;jlead < nTracks; jlead++ )
      {
        AliAODPWG4Particle *track = (AliAODPWG4Particle*) bgTracks->At(jlead) ;
        
        ptAssoc  = track->Pt();
        phiAssoc = track->Phi() ;
        if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
        
        if (fMakeNearSideLeading)
        {
          deltaPhi = phiTrig-phiAssoc;
          if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
          if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
          
          if(ptAssoc > ptTrig && deltaPhi < TMath::PiOver2())
          {
            leading = kFALSE;
            break;
          }
        }
        //jump out this event if there is any other particle with pt larger than trigger
        else if(fMakeAbsoluteLeading)
        {
          if(ptAssoc > ptTrig)
          {
            leading = kFALSE;
            break;
          }
        }
      }
      
      if( !neutralMix && fCheckLeadingWithNeutralClusters )
        AliWarning("Leading of clusters requested but no clusters in mixed event");
      
      if(neutralMix && fCheckLeadingWithNeutralClusters && bgCalo)
      {
        Int_t nClusters=bgCalo->GetEntriesFast();
        for(Int_t jlead = 0;jlead <nClusters; jlead++ )
        {
          AliAODPWG4Particle *cluster= (AliAODPWG4Particle*) bgCalo->At(jlead) ;
          
          ptAssoc  = cluster->Pt();
          phiAssoc = cluster->Phi() ;
          if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
          
          if (fMakeNearSideLeading)
          {
            deltaPhi = phiTrig-phiAssoc;
            if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
            if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
            
            if(ptAssoc > ptTrig && deltaPhi < TMath::PiOver2())
            {
              leading = kFALSE;
              break;
            }
          }
          //jump out this event if there is any other particle with pt larger than trigger
          else if(fMakeAbsoluteLeading)
          {
            if(ptAssoc > ptTrig)
            {
              leading = kFALSE;
              break;
            }
          }
        }
      }
      
      if(!leading) continue; // not leading, check the next event in pool
    }
    
    //
    // Fill histograms for selected triggers
    //
    
    fhEventMixBin->Fill(eventBin, GetEventWeight());
    
    //printf("\t Read Pool event %d, nTracks %d\n",ev,nTracks);
    
    fhPtTriggerMixed   ->Fill(ptTrig,           GetEventWeight());
    fhPhiTriggerMixed  ->Fill(ptTrig,  phiTrig, GetEventWeight());
    fhEtaTriggerMixed  ->Fill(ptTrig,  etaTrig, GetEventWeight());
    fhPtTriggerMixedBin->Fill(ptTrig, eventBin, GetEventWeight());
    
    if ( fCorrelVzBin ) fhPtTriggerMixedVzBin->Fill(ptTrig, GetEventVzBin(), GetEventWeight());
    
    //
    // Correlation histograms
    //
    for(Int_t j1 = 0;j1 <nTracks; j1++ )
    {
      AliAODPWG4Particle *track = (AliAODPWG4Particle*) bgTracks->At(j1) ;
      
      if(!track) continue;
      
      ptAssoc  = track->Pt();
      etaAssoc = track->Eta();
      phiAssoc = track->Phi() ;
      if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();
      
      deltaPhi = phiTrig-phiAssoc;
      if(deltaPhi < -TMath::PiOver2())  deltaPhi+=TMath::TwoPi();
      if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
      deltaEta = etaTrig-etaAssoc;
      
      AliDebug(1,Form("deltaPhi= %f, deltaEta=%f",deltaPhi, deltaEta));
      
      // Angular correlation
      fhMixDeltaPhiCharged        ->Fill(ptTrig  , deltaPhi, GetEventWeight());
      fhMixDeltaPhiDeltaEtaCharged->Fill(deltaPhi, deltaEta, GetEventWeight());
      
      //
      // Momentum imbalance
      //
      if ( (deltaPhi > fDeltaPhiMinCut)   && (deltaPhi < fDeltaPhiMaxCut)   )
      {
        xE = -ptAssoc/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
        
        if( xE < 0.)
          AliWarning(Form("Careful!!, negative xE %2.2f for right UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                          xE,deltaPhi*TMath::RadToDeg(),TMath::Cos(deltaPhi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
        
        fhMixXECharged->Fill(ptTrig, xE, GetEventWeight());
        if ( xE > 0 ) fhMixHbpXECharged->Fill(ptTrig, TMath::Log(1./xE), GetEventWeight());
      }
      
      //
      // Underlying event momentum imbalance
      //
      if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )
      {
        //Underlying event region
        Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
        Double_t uexE = -(ptAssoc/ptTrig)*TMath::Cos(randomphi);
        
        if(uexE < 0.)
          AliWarning(Form("Careful!!, negative xE %2.2f for left UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                          uexE,randomphi*TMath::RadToDeg(),TMath::Cos(randomphi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
        
        fhMixXEUeCharged->Fill(ptTrig, uexE, GetEventWeight());
      }
      
      // Set the pt associated bin for the defined bins
      Int_t assocBin   = -1;
      for(Int_t i = 0 ; i < fNAssocPtBins ; i++)
      {
        if(ptAssoc > fAssocPtBinLimit[i] && ptAssoc < fAssocPtBinLimit[i+1]) assocBin= i;
      }
      
      //
      // Assign to the histogram array a bin corresponding to a combination of pTa and vz bins
      //
      Int_t nz = 1;
      Int_t vz = 0;
      
      if(fCorrelVzBin)
      {
        nz = GetNZvertBin();
        vz = GetEventVzBin();
      }
      
      Int_t bin = assocBin*nz+vz;
      
      if(bin < 0) continue ; // this pt bin was not considered
      
      fhMixDeltaPhiChargedAssocPtBin        [bin]->Fill(ptTrig,   deltaPhi, GetEventWeight());
      fhMixDeltaPhiDeltaEtaChargedAssocPtBin[bin]->Fill(deltaPhi, deltaEta, GetEventWeight());
      
      if(fFillEtaGapsHisto)
      {
        if(TMath::Abs(deltaEta) > 0.8)
          fhMixDeltaPhiChargedAssocPtBinDEta08  [bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
        if(TMath::Abs(deltaEta) < 0.01)
          fhMixDeltaPhiChargedAssocPtBinDEta0   [bin]->Fill(ptTrig, deltaPhi, GetEventWeight());
      }
      
    } // track loop
  } // mixed event loop
}

//_______________________________________________________________________________________________________
/// Correlate the trigger with pi0 selected in another task (AliAnaPi0EbE).
/// Similar procedure as for charged tracks.
//_______________________________________________________________________________________________________
void AliAnaParticleHadronCorrelation::MakeNeutralCorrelation(AliAODPWG4ParticleCorrelation * aodParticle)
{
  TObjArray * pi0list = (TObjArray*) GetAODBranch(fPi0AODBranchName); // For the future, foresee more possible pi0 lists
  if(!pi0list) return ;
  
  Int_t npi0 = pi0list->GetEntriesFast();
  if(npi0 == 0) return ;
  
  AliDebug(1,Form("Particle - pi0 correlation, %d pi0's",npi0));
  
  Int_t evtIndex11 = 0 ;
  Int_t evtIndex12 = 0 ;
  if (GetMixedEvent())
  {
    evtIndex11 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(0)) ;
    evtIndex12 = GetMixedEvent()->EventIndexForCaloCluster(aodParticle->GetCaloLabel(1)) ;
  }
  
  Float_t pt   = -100. ;
  Float_t zT   = -100. ;
  Float_t phi  = -100. ;
  Float_t eta  = -100. ;
  Float_t xE   = -100. ;
  Float_t hbpXE= -100. ;
  Float_t hbpZT= -100. ;
  
  Float_t ptTrig  = aodParticle->Pt();
  Float_t phiTrig = aodParticle->Phi();
  Float_t etaTrig = aodParticle->Eta();
  Float_t deltaPhi= -100. ;
  Float_t deltaEta= -100. ;
  
  // In case of pi0/eta trigger, we may want to check their decay correlation,
  // get their decay children
  
  Bool_t decayFound = kFALSE;
  if(fPi0Trigger) decayFound = GetDecayPhotonMomentum(aodParticle->GetCaloLabel(0),aodParticle->GetCaloLabel(1),aodParticle->GetDetectorTag());
  
  TObjArray * refpi0 = 0x0;
  Int_t nrefs        = 0;
  
  // Loop on stored AOD pi0
  
  for(Int_t iaod = 0; iaod < npi0 ; iaod++)
  {
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (pi0list->At(iaod));
    
    Int_t evtIndex2 = 0 ;
    Int_t evtIndex3 = 0 ;
    if (GetMixedEvent())
    {
      evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(0)) ;
      evtIndex3 = GetMixedEvent()->EventIndexForCaloCluster(pi0->GetCaloLabel(1)) ;
      
      if (evtIndex11 == evtIndex2 || evtIndex12 == evtIndex2 ||
          evtIndex11 == evtIndex3 || evtIndex12 == evtIndex3) // trigger and pi0 are not from different events
        continue ;
    }
    
    pt  = pi0->Pt();
    
    if(pt < fMinAssocPt || pt > fMaxAssocPt) continue ;
    
    // Remove trigger itself for correlation when use charged triggers
    if(aodParticle->GetCaloLabel(0) >= 0 &&
       (pi0->GetCaloLabel(0) == aodParticle->GetCaloLabel(0) || pi0->GetCaloLabel(1) == aodParticle->GetCaloLabel(0))) continue ;
    
    if( aodParticle->GetCaloLabel(1) >= 0 &&
       (pi0->GetCaloLabel(0) == aodParticle->GetCaloLabel(1) || pi0->GetCaloLabel(1) == aodParticle->GetCaloLabel(1))) continue ;
    
    //
    // Angular correlations
    //
    phi      = pi0->Phi() ;
    eta      = pi0->Eta() ;
    deltaEta = etaTrig-eta;
    deltaPhi = phiTrig-phi;
    if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
    
    fhEtaNeutral     ->Fill(pt    , eta     , GetEventWeight());
    fhPhiNeutral     ->Fill(pt    , phi     , GetEventWeight());
    fhDeltaEtaNeutral->Fill(ptTrig, deltaEta, GetEventWeight());
    fhDeltaPhiNeutral->Fill(ptTrig, deltaPhi, GetEventWeight());
    
    if ( pt > 2 ) fhDeltaPhiDeltaEtaNeutral->Fill(deltaPhi, deltaEta, GetEventWeight());
    
    //
    // Momentum imbalance
    //
    zT  = pt/ptTrig ;
    
    hbpZT = -100;
    hbpXE = -100;
    
    if(zT > 0 ) hbpZT = TMath::Log(1./zT);
    
    // Delta phi cut for correlation
    if( (deltaPhi > fDeltaPhiMinCut) && ( deltaPhi < fDeltaPhiMaxCut) )
    {
      xE  =-pt/ptTrig*TMath::Cos(deltaPhi); // -(px*pxTrig+py*pyTrig)/(ptTrig*ptTrig);
      
      if( xE < 0. )
        AliWarning(Form("Careful!!, negative xE %2.2f for right UE cos(dPhi %2.2f) = %2.2f, check correlation dPhi limits %f to %f",
                        xE,deltaPhi*TMath::RadToDeg(),TMath::Cos(deltaPhi),fDeltaPhiMinCut*TMath::RadToDeg(),fDeltaPhiMaxCut*TMath::RadToDeg()));
      
      if( xE > 0 ) hbpXE = TMath::Log(1./xE);
      
      fhDeltaPhiNeutralPt->Fill(pt    , deltaPhi, GetEventWeight());
      fhXENeutral        ->Fill(ptTrig, xE      , GetEventWeight());
      fhPtHbpXENeutral   ->Fill(ptTrig, hbpXE   , GetEventWeight());
      fhZTNeutral        ->Fill(ptTrig, zT      , GetEventWeight());
      fhPtHbpZTNeutral   ->Fill(ptTrig, hbpZT   , GetEventWeight());
    }
    else if ( (deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi < fUeDeltaPhiMaxCut) )
    {
      // Randomize angle for xE calculation
      Double_t randomphi = gRandom->Uniform(fDeltaPhiMinCut,fDeltaPhiMaxCut);
      
      xE = -(pt/ptTrig)*TMath::Cos(randomphi);
      if(xE > 0 ) hbpXE = TMath::Log(1./xE);
      
      fhDeltaPhiUeNeutralPt->Fill(pt    , deltaPhi, GetEventWeight());
      fhZTUeNeutral        ->Fill(ptTrig, zT      , GetEventWeight());
      fhPtHbpZTUeNeutral   ->Fill(ptTrig, hbpZT   , GetEventWeight());
      fhXEUeNeutral        ->Fill(ptTrig, xE      , GetEventWeight());
      fhPtHbpXEUeNeutral   ->Fill(ptTrig, hbpXE   , GetEventWeight());
    }
    
    // Several UE calculation, not sure it is useful
    // with partical calorimter acceptance
    if(fMakeSeveralUE) FillNeutralUnderlyingEventSidesHistograms(ptTrig,pt,zT,hbpZT,deltaPhi);
    
    //
    // Decay photon correlations
    //
    if(fPi0Trigger && decayFound)
      FillDecayPhotonCorrelationHistograms(pt, phi, kFALSE) ;
    
    if(fFillAODWithReferences)
    {
      nrefs++;
      if(nrefs==1)
      {
        refpi0 = new TObjArray(0);
        refpi0->SetName(GetAODObjArrayName()+"Pi0s");
        refpi0->SetOwner(kFALSE);
      }
      refpi0->Add(pi0);
    } // put references in trigger AOD
    
    AliDebug(1,Form("Selected pi0: pt %2.2f, phi %2.2f, eta %2.2f",pt,phi,eta));
    
  } // loop
  
  // Fill AOD with reference tracks, if not filling histograms
  if(fFillAODWithReferences && refpi0)
  {
    aodParticle->AddObjArray(refpi0);
  }
}

//__________________________________________________________________________________________________________________
/// Make the trigger-charged particles correlation at the generator level.
//__________________________________________________________________________________________________________________
void  AliAnaParticleHadronCorrelation::MakeMCChargedCorrelation(Int_t label, Int_t histoIndex, Bool_t lostDecayPair)
{
  AliDebug(1,"Make trigger particle - charged hadron correlation in AOD MC level");
  
  if( label < 0 )
  {
    AliDebug(1,Form(" *** bad label ***:  label %d", label));
    return;
  }
  
  if ( !GetMC() )
  {
    AliFatal("Stack not available, is the MC handler called? STOP");
    return;
  }
  
  //Int_t nTracks = GetMC()->GetNumberOfTracks() ;
  Int_t nTracks = GetMC()->GetNumberOfPrimaries();
  
  if( label >= nTracks )
  {
    if(GetDebug() > 2)
      AliInfo(Form(" *** large label ***:  label %d, n prim tracks %d", label,nTracks));
    
    return;
  }

  // Do MC correlation for a given particle type range.
  // Types defined in GetMCTagHistogramIndex:
  // 0 direct gamma; 1 pi0; 2 pi0 decay; 3 eta decay; 4 other decay; 5 electron; 6 other (hadron)
  if(histoIndex < fMCGenTypeMin || histoIndex > fMCGenTypeMax) return ;
    
  Double_t eprim   = 0 ;
  Double_t ptprim  = 0 ;
  Double_t phiprim = 0 ;
  Double_t etaprim = 0 ;
  Int_t iParticle  = 0 ;
    
  // Get the particle
  AliVParticle * primary = GetMC()->GetTrack(label);
  if( !primary )
  {
    AliInfo(Form(" *** no primary ***:  label %d", label));
    return;
  }
  
  eprim   = primary->E();
  ptprim  = primary->Pt();
  etaprim = primary->Eta();
  phiprim = primary->Phi();
  if(phiprim < 0) phiprim+=TMath::TwoPi();
  
  if(ptprim < 0.01 || eprim < 0.01) return ;
  
  Bool_t leadTrig = kTRUE; 

  for (iParticle = 0; iParticle < nTracks; iParticle++)
  {
    if ( !GetReader()->AcceptParticleMCLabel( iParticle ) ) continue ;
    
    AliVParticle * particle = GetMC()->GetTrack(iParticle);
    
    if ( !particle->IsPhysicalPrimary() ) continue; // same as particle->MCStatusCode() !=1
    
    if (  particle->Charge() == 0 ) continue;
    
    particle->Momentum(fMomentum);
    //fMomentum.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->E());
    
    // Particles in CTS acceptance, make sure to use the same selection as in the reader
    Bool_t inCTS = GetReader()->GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kCTS);
    //printf("Accepted? %d; pt %2.2f, eta %2.2f, phi %2.2f\n",inCTS,fMomentum.Pt(),fMomentum.Eta(),fMomentum.Phi()*TMath::RadToDeg());
    if( !inCTS ) continue;
    
    // Remove conversions
    Int_t indexmother = particle->GetMother();
    if ( indexmother > -1 )
    {
      Int_t pdg = particle->PdgCode();
      Int_t mPdg = (GetMC()->GetTrack(indexmother))->PdgCode();
      if (TMath::Abs(pdg) == 11 && mPdg == 22) continue;
    }
    
    if ( label == iParticle ) continue; // avoid trigger particle
    
    Float_t phi = particle->Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    
    Bool_t lead = FillChargedMCCorrelationHistograms(particle->Pt(),phi,particle->Eta(),ptprim,phiprim,etaprim, histoIndex,lostDecayPair);
    
    if ( !lead ) leadTrig = kFALSE;
    //    if ( !lead && (fMakeAbsoluteLeading || fMakeNearSideLeading)) return;
  }  //MC particle loop
  
  // Trigger MC particle histograms
  //if (!lead  && (fMakeAbsoluteLeading || fMakeNearSideLeading)) return;
  
  fhMCPtTrigger [histoIndex]->Fill(ptprim,          GetEventWeight());
  fhMCPhiTrigger[histoIndex]->Fill(ptprim, phiprim, GetEventWeight());
  fhMCEtaTrigger[histoIndex]->Fill(ptprim, etaprim, GetEventWeight());
  
  if(lostDecayPair)
  {
    // check index in GetMCTagIndexHistogram
    if     (histoIndex == 2 && 8 >= fMCGenTypeMin && 8 <= fMCGenTypeMax )
    {
      // pi0 decay
      fhMCPtTrigger [8]->Fill(ptprim,          GetEventWeight());
      fhMCPhiTrigger[8]->Fill(ptprim, phiprim, GetEventWeight());
      fhMCEtaTrigger[8]->Fill(ptprim, etaprim, GetEventWeight());
    }
    else if(histoIndex == 4 && 9 >= fMCGenTypeMin && 9 <= fMCGenTypeMax )
    {
      // eta decay
      fhMCPtTrigger [9]->Fill(ptprim,          GetEventWeight());
      fhMCPhiTrigger[9]->Fill(ptprim, phiprim, GetEventWeight());
      fhMCEtaTrigger[9]->Fill(ptprim, etaprim, GetEventWeight());
    }
  }
  
  if(!leadTrig && (fMakeAbsoluteLeading || fMakeNearSideLeading) )
  {
    AliDebug(1,Form("Not leading primary trigger: pT %2.2f, phi %2.2f, eta %2.2f",
                    ptprim,phiprim*TMath::RadToDeg(),etaprim));
    
    fhMCPtTriggerNotLeading [histoIndex]->Fill(ptprim,          GetEventWeight());
    fhMCPhiTriggerNotLeading[histoIndex]->Fill(ptprim, phiprim, GetEventWeight());
    fhMCEtaTriggerNotLeading[histoIndex]->Fill(ptprim, etaprim, GetEventWeight());
    
    if(lostDecayPair)
    {
      // check index in GetMCTagIndexHistogram
      if      (histoIndex == 2  && 8 >= fMCGenTypeMin && 8 <= fMCGenTypeMax )
      {
        // pi0 decay
        fhMCPtTriggerNotLeading [8]->Fill(ptprim,          GetEventWeight());
        fhMCPhiTriggerNotLeading[8]->Fill(ptprim, phiprim, GetEventWeight());
        fhMCEtaTriggerNotLeading[8]->Fill(ptprim, etaprim, GetEventWeight());
      }
      else  if(histoIndex == 4  && 9 >= fMCGenTypeMin && 9 <= fMCGenTypeMax )
      {
        // eta decay
        fhMCPtTriggerNotLeading [9]->Fill(ptprim,          GetEventWeight());
        fhMCPhiTriggerNotLeading[9]->Fill(ptprim, phiprim, GetEventWeight());
        fhMCEtaTriggerNotLeading[9]->Fill(ptprim, etaprim, GetEventWeight());
      }
    }
  }
}

//_____________________________________________________________________
/// Print some relevant parameters set for the analysis.
//_____________________________________________________________________
void AliAnaParticleHadronCorrelation::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  printf("Pt trigger > %2.2f; < %2.2f\n", GetMinPt() , GetMaxPt()) ;
  printf("Pt associa > %2.2f; < %2.2f\n", fMinAssocPt, fMaxAssocPt) ;
  printf("Phi trigger particle-Hadron   <  %3.2f\n", fDeltaPhiMaxCut) ;
  printf("Phi trigger particle-Hadron   >  %3.2f\n", fDeltaPhiMinCut) ;
  printf("Phi trigger particle-UeHadron <  %3.2f\n", fUeDeltaPhiMaxCut) ;
  printf("Phi trigger particle-UeHadron >  %3.2f\n", fUeDeltaPhiMinCut) ;
  printf("Isolated Trigger?  %d\n"     , fSelectIsolated) ;
  printf("Several UE?  %d\n"           , fMakeSeveralUE) ;
  printf("Name of AOD Pi0 Branch %s \n", fPi0AODBranchName.Data());
  printf("Do Decay-hadron correlation ?  %d\n", fPi0Trigger) ;
  printf("Select absolute leading for cluster triggers ?  %d or Near Side %d\n",
         fMakeAbsoluteLeading, fMakeNearSideLeading) ;
  printf("Trigger pt bins  %d\n", fNAssocPtBins) ;
  for (Int_t ibin = 0; ibin<fNAssocPtBins; ibin++) {
    printf("\t bin %d = [%2.1f,%2.1f]\n", ibin, fAssocPtBinLimit[ibin], fAssocPtBinLimit[ibin+1]) ;
  }
  
}

//____________________________________________________________
/// Set number of associated charged (neutral) hadrons pT bins.
//____________________________________________________________
void AliAnaParticleHadronCorrelation::SetNAssocPtBins(Int_t n)
{
  fNAssocPtBins  = n ;
  
  if(n < 20 && n > 0)
  {
    fNAssocPtBins  = n ;
  }
  else
  {
    AliWarning("n = larger than 19 or too small, set to 19");
    fNAssocPtBins = 19;
  }
}

//______________________________________________________________________________
/// Set the list of pT limits for the of associated charged (neutral) hadrons.
//______________________________________________________________________________
void AliAnaParticleHadronCorrelation::SetAssocPtBinLimit(Int_t ibin, Float_t pt)
{
  if(ibin <= fNAssocPtBins || ibin >= 0)
  {
    fAssocPtBinLimit[ibin] = pt  ;
  }
  else
  {
    AliWarning(Form("Bin  number too large %d > %d or small, nothing done", ibin, fNAssocPtBins)) ;
  }
}

