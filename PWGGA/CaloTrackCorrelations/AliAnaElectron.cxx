 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>
#include "TDatabasePDG.h"
#include "AliVTrack.h"

// --- Analysis system --- 
#include "AliAnaElectron.h" 
#include "AliCaloTrackReader.h"
#include "AliMCEvent.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
#include "AliPIDResponse.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"

/// \cond CLASSIMP
ClassImp(AliAnaElectron) ;
/// \endcond

//________________________________
/// Default constructor. Initialize parameters.
//________________________________
AliAnaElectron::AliAnaElectron() :
AliAnaCaloTrackCorrBaseClass(),
fMinDist(0.),                        
fTimeCutMin(-1),                      fTimeCutMax(999999),         
fNCellsCut(0),                        fNLMCutMin(-1),                        fNLMCutMax(10),
fFillSSHistograms(kFALSE),            fFillOnlySimpleSSHisto(1),
fFillWeightHistograms(kFALSE),        fNOriginHistograms(8), 
fdEdxMin(0.),                         fdEdxMax(400.), 
fEOverPMin(0),                        fEOverPMax(2),
fNSigmaMin(-100),                     fNSigmaMax(100),
fM20Min(-1),                          fM20Max(100),
fM02Min(0),                           fM02Max(100),
fdEdxMinHad(0.),                      fdEdxMaxHad(400.), 
fEOverPMinHad(0),                     fEOverPMaxHad(100),
fNSigmaMinHad(-100),                  fNSigmaMaxHad(100),
fAODParticle(0),
fMomentum(),                          fMomentumMC(),                         fProdVertex(),

// Histograms
fhdEdxvsE(0),                         fhdEdxvsP(0),                 
fhdEdxvsECutM02(0),                   fhdEdxvsPCutM02(0),
fhdEdxvsECutM02AndM20(0),             fhdEdxvsPCutM02AndM20(0),
fhdEdxvsECutEOverP(0),                fhdEdxvsPCutEOverP(0),
fhdEdxvsECutNSigma(0),                fhdEdxvsPCutNSigma(0),
fhdEdxvsECutM02CutNSigma(0),          fhdEdxvsPCutM02CutNSigma(0),
fhdEdxvsECutM02AndM20CutNSigma(0),    fhdEdxvsPCutM02AndM20CutNSigma(0),

fhEOverPvsE(0),                       fhEOverPvsP(0),
fhEOverPvsECutM02(0),                 fhEOverPvsPCutM02(0),
fhEOverPvsECutM02AndM20(0),           fhEOverPvsPCutM02AndM20(0),
fhEOverPvsECutdEdx(0),                fhEOverPvsPCutdEdx(0),
fhEOverPvsECutM02CutdEdx(0),          fhEOverPvsPCutM02CutdEdx(0),
fhEOverPvsECutM02AndM20CutdEdx(0),    fhEOverPvsPCutM02AndM20CutdEdx(0),
fhEOverPvsECutNSigma(0),              fhEOverPvsPCutNSigma(0),
fhEOverPvsECutM02CutNSigma(0),        fhEOverPvsPCutM02CutNSigma(0),
fhEOverPvsECutM02AndM20CutNSigma(0),  fhEOverPvsPCutM02AndM20CutNSigma(0),

fhNSigmavsE(0),                       fhNSigmavsP(0),
fhNSigmavsECutM02(0),                 fhNSigmavsPCutM02(0),
fhNSigmavsECutM02AndM20(0),           fhNSigmavsPCutM02AndM20(0),
fhNSigmavsECutdEdx(0),                fhNSigmavsPCutdEdx(0),
fhNSigmavsECutM02CutdEdx(0),          fhNSigmavsPCutM02CutdEdx(0),
fhNSigmavsECutM02AndM20CutdEdx(0),    fhNSigmavsPCutM02AndM20CutdEdx(0),
fhNSigmavsECutEOverP(0),              fhNSigmavsPCutEOverP(0),

// Weight studies
fhECellClusterRatio(0),               fhECellClusterLogRatio(0),                 
fhEMaxCellClusterRatio(0),            fhEMaxCellClusterLogRatio(0), 

// MC histograms
// Electron SS MC histograms
fhMCElectronELambda0NoOverlap(0),    
fhMCElectronELambda0TwoOverlap(0),    fhMCElectronELambda0NOverlap(0),

//Embedding
fhEmbeddedSignalFractionEnergy(0),     
fhEmbedElectronELambda0FullSignal(0), fhEmbedElectronELambda0MostlySignal(0),  
fhEmbedElectronELambda0MostlyBkg(0),  fhEmbedElectronELambda0FullBkg(0)        
{
  for(Int_t index = 0; index < 2; index++)
  {
    fhNCellsE [index] = 0;
    fhNLME    [index] = 0;
    fhTimeE   [index] = 0;
    fhMaxCellDiffClusterE[index] = 0;
    fhE       [index] = 0;    
    fhPt      [index] = 0;                        
    fhPhi     [index] = 0;                      
    fhEta     [index] = 0; 
    fhEtaPhi  [index] = 0;                   
    fhEtaPhi05[index] = 0;
    
    // Shower shape histograms
    fhDispE   [index] = 0;                    
    fhLam0E   [index] = 0;                    
    fhLam1E   [index] = 0; 
    fhDispETRD[index] = 0;                 
    fhLam0ETRD[index] = 0;                 
    fhLam1ETRD[index] = 0;
    fhNCellsLam0LowE [index] = 0;           
    fhNCellsLam0HighE[index] = 0;       
    fhEtaLam0LowE    [index] = 0;              
    fhPhiLam0LowE    [index] = 0; 
    fhEtaLam0HighE   [index] = 0;             
    fhPhiLam0HighE   [index] = 0; 
    
    fhDispEtaE       [index] = 0;                
    fhDispPhiE       [index] = 0;
    fhSumEtaE        [index] = 0;                
    fhSumPhiE        [index] = 0;                
    fhSumEtaPhiE     [index] = 0;
    fhDispEtaPhiDiffE[index] = 0;         
    fhSphericityE    [index] = 0;
    
    for(Int_t i = 0; i < 10; i++)
    {
      fhMCPt     [index][i] = 0;
      fhMCE      [index][i] = 0;
      fhMCPhi    [index][i] = 0;
      fhMCEta    [index][i] = 0;
      fhMCDeltaE [index][i] = 0;                
      fhMC2E     [index][i] = 0;
      fhMCdEdxvsE       [i] = 0;
      fhMCdEdxvsP       [i] = 0;     
      fhMCNSigmavsE     [i] = 0;
      fhMCNSigmavsP     [i] = 0;
      fhMCEOverPvsE     [i] = 0;
      fhMCEOverPvsP     [i] = 0;    
      fhMCEOverPvsEAfterCuts[i][index] = 0;
      fhMCEOverPvsPAfterCuts[i][index] = 0;
    }
    
    for(Int_t i = 0; i < 6; i++)
    {
      fhMCELambda0       [index][i] = 0;
      fhMCELambda1       [index][i] = 0;
      fhMCEDispEta       [index][i] = 0;
      fhMCEDispPhi       [index][i] = 0;
      fhMCESumEtaPhi     [index][i] = 0;
      fhMCEDispEtaPhiDiff[index][i] = 0;
      fhMCESphericity    [index][i] = 0;
    }
    
    for(Int_t i = 0; i < 5; i++)
    {
      fhDispEtaDispPhiEBin[index][i] = 0 ;
    }
  }
  
  // Mathching Residuals
  for(Int_t indexPID = 0; indexPID < 3; indexPID++)
  {
    for(Int_t ich = 0; ich < 2; ich++)
    { 
      fhDEtavsE [indexPID][ich] = 0 ;      
      fhDPhivsE [indexPID][ich] = 0 ;     
      fhDEtavsP [indexPID][ich] = 0 ;       
      fhDPhivsP [indexPID][ich] = 0 ;       
      fhDEtaDPhi[indexPID][ich] = 0 ;
    }
  }
  
  // Weight studies
  for(Int_t i =0; i < 14; i++)
  {
    fhLambda0ForW0[i] = 0;
    //fhLambda1ForW0[i] = 0;
  }
  
  // Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
/// Select calorimeter clusters if they pass different cuts:
///  * Energy (if stricter cut than in AliCaloTrackReader)
///  * Time (but usually it is already done in AliCaloTrackReader)
///  * Number of cells in cluster
///  * Number of local maxima in cluster
///  * Fiducial cut, eta-phi acceptance cut via AliFiducialCut
///  * Charged clusters are rejected (if requested)
///  * Reject clusters close to a bad channel
///
/// Fill for each of the cuts a 1 dimensional histogram with either the energy
/// or the transverse momentum of the cluster. Also track-matching control histograms
/// can be filled with residuals of the matching.
/// \return kTRUE of cluster is accepted
/// \param calo: cluster pointer.
/// \param nMaxima: number of local maxima.
//____________________________________________________________________________
Bool_t  AliAnaElectron::ClusterSelected(AliVCluster* calo, Int_t nMaxima)
{
  AliDebug(1,Form("Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f",
                  GetReader()->GetEventNumber(),fMomentum.E(),fMomentum.Pt(),
                  calo->E(),GetPhi(fMomentum.Phi())*TMath::RadToDeg(),fMomentum.Eta()));
  
  //.......................................
  // If too small or big energy, skip it
  if ( fMomentum.E() < GetMinEnergy() || 
       fMomentum.E() > GetMaxEnergy()   ) return kFALSE ; 
  AliDebug(2,Form("\t Cluster %d Pass E Cut",calo->GetID()));
  
  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if ( tof < fTimeCutMin || tof > fTimeCutMax ) return kFALSE;
  AliDebug(2,Form("\t Cluster %d Pass Time Cut",calo->GetID()));
  
  //.......................................
  if ( calo->GetNCells() <= fNCellsCut && 
       GetReader()->GetDataType() != AliCaloTrackReader::kMC ) return kFALSE;
  AliDebug(2,Form("\t Cluster %d Pass NCell Cut",calo->GetID()));
  
  //.......................................
  // Check acceptance selection
  if ( IsFiducialCutOn() )
  {
    Bool_t in = GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),GetCalorimeter()) ;
    if ( !in ) return kFALSE ;
  }
  AliDebug(2,"\t Fiducial cut passed");
  
  //.......................................
  // Skip not matched clusters with tracks
  if ( !IsTrackMatched(calo, GetReader()->GetInputEvent()) )
  {
      AliDebug(1,"\t Reject non track-matched clusters");
      return kFALSE ;
  }
  else AliDebug(2,"\t Track-matching cut passed");
  
  //...........................................
  // Skip clusters with too many maxima
  if ( nMaxima < fNLMCutMin || nMaxima > fNLMCutMax ) return kFALSE ;
  AliDebug(2,Form("\t Cluster %d pass NLM %d of out of range",calo->GetID(), nMaxima));
  
  //.......................................
  // Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if ( distBad < 0.       ) distBad=9999. ; //workout strange convension dist = -1. ;
  if ( distBad < fMinDist ) //In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  else AliDebug(2,Form("\t Bad channel cut passed %4.2f > %2.2f",distBad, fMinDist));
  //printf("Cluster %d Pass Bad Dist Cut \n",icalo);

  AliDebug(1,Form("Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f",
                  GetReader()->GetEventNumber(), 
                  fMomentum.E(), fMomentum.Pt(),fMomentum.E(),
                  GetPhi(fMomentum.Phi())*TMath::RadToDeg(),fMomentum.Eta()));
  
  // All checks passed, cluster selected
  return kTRUE;
    
}

//______________________________________________________________________________________________
/// Fill cluster Shower Shape histograms.
//______________________________________________________________________________________________
void  AliAnaElectron::FillShowerShapeHistograms(AliVCluster* cluster, Int_t mcTag, Int_t pidTag)
{
  if ( !fFillSSHistograms || GetMixedEvent() ) return;
  
  Int_t pidIndex = 0;// Electron
  if      ( pidTag == AliCaloPID::kElectron      ) pidIndex = 0;
  else if ( pidTag == AliCaloPID::kChargedHadron ) pidIndex = 1;
  else return;

  Float_t energy  = cluster->E();
  Int_t   ncells  = cluster->GetNCells();
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  Float_t disp    = cluster->GetDispersion()*cluster->GetDispersion();
  
  Float_t l0   = 0., l1   = 0.;
  Float_t dispp= 0., dEta = 0., dPhi    = 0.; 
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  
  Float_t eta = fMomentum.Eta();
  Float_t phi = GetPhi(fMomentum.Phi());
  
  fhLam0E[pidIndex] ->Fill(energy, lambda0, GetEventWeight());
  fhLam1E[pidIndex] ->Fill(energy, lambda1, GetEventWeight());
   if ( !fFillOnlySimpleSSHisto ) fhDispE[pidIndex] ->Fill(energy, disp   , GetEventWeight());
  
  if ( GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 &&
       GetModuleNumber(cluster) >= GetFirstSMCoveredByTRD() )
  {
    fhLam0ETRD[pidIndex]->Fill(energy, lambda0, GetEventWeight());
    fhLam1ETRD[pidIndex]->Fill(energy, lambda1, GetEventWeight());
     if ( !fFillOnlySimpleSSHisto ) fhDispETRD[pidIndex]->Fill(energy, disp   , GetEventWeight());
  }
  
  if ( !fFillOnlySimpleSSHisto )
  {
    if ( energy < 2 )
    {
      fhNCellsLam0LowE[pidIndex] ->Fill(ncells, lambda0, GetEventWeight());
      fhEtaLam0LowE[pidIndex]    ->Fill(eta,    lambda0, GetEventWeight());
      fhPhiLam0LowE[pidIndex]    ->Fill(phi,    lambda0, GetEventWeight());
    }
    else 
    {
      fhNCellsLam0HighE[pidIndex]->Fill(ncells, lambda0, GetEventWeight());
      fhEtaLam0HighE[pidIndex]   ->Fill(eta,    lambda0, GetEventWeight());
      fhPhiLam0HighE[pidIndex]   ->Fill(phi,    lambda0, GetEventWeight());
    }
    
    if ( GetCalorimeter() == kEMCAL )
    {
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                                   l0, l1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);
      fhDispEtaE        [pidIndex]-> Fill(energy, dEta     , GetEventWeight());
      fhDispPhiE        [pidIndex]-> Fill(energy, dPhi     , GetEventWeight());
      fhSumEtaE         [pidIndex]-> Fill(energy, sEta     , GetEventWeight());
      fhSumPhiE         [pidIndex]-> Fill(energy, sPhi     , GetEventWeight());
      fhSumEtaPhiE      [pidIndex]-> Fill(energy, sEtaPhi  , GetEventWeight());
      fhDispEtaPhiDiffE [pidIndex]-> Fill(energy, dPhi-dEta, GetEventWeight());
      if ( dEta+dPhi > 0 )
          fhSphericityE [pidIndex]-> Fill(energy, (dPhi-dEta)/(dEta+dPhi), GetEventWeight());
      
      if      (energy < 2 ) fhDispEtaDispPhiEBin[pidIndex][0]->Fill(dEta, dPhi, GetEventWeight());
      else if (energy < 4 ) fhDispEtaDispPhiEBin[pidIndex][1]->Fill(dEta, dPhi, GetEventWeight());
      else if (energy < 6 ) fhDispEtaDispPhiEBin[pidIndex][2]->Fill(dEta, dPhi, GetEventWeight());
      else if (energy < 10) fhDispEtaDispPhiEBin[pidIndex][3]->Fill(dEta, dPhi, GetEventWeight());
      else                  fhDispEtaDispPhiEBin[pidIndex][4]->Fill(dEta, dPhi, GetEventWeight());
    }
  }
  
  if ( IsDataMC() )
  {
    AliVCaloCells* cells = 0;
    if ( GetCalorimeter() == kEMCAL ) cells = GetEMCALCells();
    else                              cells = GetPHOSCells();
    
    // Fill histograms to check shape of embedded clusters
    Float_t fraction = 0;
    if (  IsEmbedingAnalysisOn() )
    {
      //Only working for EMCAL
      Float_t clusterE = 0; // recalculate in case corrections applied.
      Float_t cellE    = 0;
      if ( !GetReader()->IsEmbeddedMCEventUsed() )
      {
        for(Int_t icell = 0; icell < cluster->GetNCells(); icell++)
        {
          cellE    = cells->GetCellAmplitude(cluster->GetCellAbsId(icell));
          clusterE+=cellE;  
          fraction+=cellE*cluster->GetCellAmplitudeFraction(icell);
        }
        
        // Fraction of total energy due to the embedded signal
        fraction/=clusterE;
      }
      else if ( !GetReader()->IsEmbeddedInputEventUsed() )
      {
        Float_t sigCellE    = 0; 
        Float_t sigClusterE = 0; 
        for(Int_t icell  = 0; icell < cluster->GetNCells(); icell++)
        {
          Int_t id = cluster->GetCellAbsId(icell);
          
          cellE    = cells->GetCellAmplitude(id);
          sigCellE = cells->GetCellEFraction(id); // MC signal
          clusterE   += cellE   ;
          sigClusterE+= sigCellE;
        }
        
        // Fraction of total energy due to the embedded signal
        fraction = sigClusterE / clusterE;
      }
      
      AliDebug(1,Form("Energy fraction of embedded signal %2.3f, Energy %2.3f",fraction, clusterE));
      
      fhEmbeddedSignalFractionEnergy->Fill(clusterE, fraction, GetEventWeight());
    }  // embedded fraction    
    
    // Check the origin and fill histograms
    Int_t index = -1;

    if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) &&
        !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0) &&
        !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta) )
    {
      index = kmcssPhoton;
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron) && 
              !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion))
    {
      index = kmcssElectron;
       
      if ( !IsEmbedingAnalysisOn() )
      {
        //Check particle overlaps in cluster
        
        //Compare the primary depositing more energy with the rest, if no photon/electron as comon ancestor (conversions), count as other particle
        Int_t ancPDG = 0, ancStatus = -1;
        Int_t ancLabel = 0;
        Int_t noverlaps = 1;      
        for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++)
        {
          ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(cluster->GetLabels()[0],cluster->GetLabels()[ilab], GetMC(),
                                                               ancPDG,ancStatus,fMomentumMC,fProdVertex);
          if ( ancPDG!=22 && TMath::Abs(ancPDG)!=11 ) noverlaps++;
        }
        
        if      ( noverlaps == 1 ) {
          fhMCElectronELambda0NoOverlap  ->Fill(energy, lambda0, GetEventWeight());
        }
        else if ( noverlaps == 2 ) {        
          fhMCElectronELambda0TwoOverlap ->Fill(energy, lambda0, GetEventWeight());
        }
        else if ( noverlaps >  2 ) {          
          fhMCElectronELambda0NOverlap   ->Fill(energy, lambda0, GetEventWeight());
        }
        else
        {
          AliWarning(Form("N overlaps = %d for ancestor %d!!", noverlaps, ancLabel));
        }
      }//No embedding
      
      //Fill histograms to check shape of embedded clusters
      if ( IsEmbedingAnalysisOn() )
      {
        if      ( fraction > 0.9 ) 
        {
          fhEmbedElectronELambda0FullSignal   ->Fill(energy, lambda0, GetEventWeight());
        }
        else if ( fraction > 0.5 )
        {
          fhEmbedElectronELambda0MostlySignal ->Fill(energy, lambda0, GetEventWeight());
        }
        else if ( fraction > 0.1 )
        { 
          fhEmbedElectronELambda0MostlyBkg    ->Fill(energy, lambda0, GetEventWeight());
        }
        else
        {
          fhEmbedElectronELambda0FullBkg      ->Fill(energy, lambda0, GetEventWeight());
        }
      } // embedded      
    } // electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron) && 
               GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
    {
      index = kmcssConversion;
    } // conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  )
    {
      index = kmcssPi0;
    } // pi0
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)  )
    {
      index = kmcssEta;      
    } // eta
    else 
    {
      index = kmcssOther;
    } // other particles
    
    fhMCELambda0[pidIndex][index]    ->Fill(energy, lambda0, GetEventWeight());
    fhMCELambda1[pidIndex][index]    ->Fill(energy, lambda1, GetEventWeight());
    
    if ( GetCalorimeter() == kEMCAL && !fFillOnlySimpleSSHisto )
    {
      fhMCEDispEta        [pidIndex][index]-> Fill(energy, dEta     , GetEventWeight());
      fhMCEDispPhi        [pidIndex][index]-> Fill(energy, dPhi     , GetEventWeight());
      fhMCESumEtaPhi      [pidIndex][index]-> Fill(energy, sEtaPhi  , GetEventWeight());
      fhMCEDispEtaPhiDiff [pidIndex][index]-> Fill(energy, dPhi-dEta, GetEventWeight());
      if ( dEta+dPhi > 0 )
          fhMCESphericity [pidIndex][index]-> Fill(energy, (dPhi-dEta)/(dEta+dPhi), GetEventWeight());
    }
  } // MC data
}

//_____________________________________________
/// Save parameters used for analysis.
//_____________________________________________
TObjString *  AliAnaElectron::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaElectron ---: ") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," %2.2f < dEdx < %2.2f;",fdEdxMin,fdEdxMax) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize," %2.2f <  E/P < %2.2f;",fEOverPMin, fEOverPMax) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize," %2.2f <  M20 < %2.2f;",fM20Min, fM20Max) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize," %2.2f < dEdx < %2.2f;",fdEdxMinHad,fdEdxMaxHad) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize," %2.2f <  E/P < %2.2f;",fEOverPMinHad, fEOverPMaxHad) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"fMinDist =%2.2f;",fMinDist) ;
  parList+=onePar ;
 
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  //parList += GetCaloPID()->GetPIDParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList()
  
  return new TObjString(parList) ;
}

//_______________________________________________
// Create histograms to be saved in output file and
// store them in outputContainer.
//_______________________________________________
TList *  AliAnaElectron::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 
	
  Int_t nptbins     = GetHistogramRanges()->GetHistoPtBins();           Float_t ptmax     = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin     = GetHistogramRanges()->GetHistoPtMin(); 
  Int_t nphibins    = GetHistogramRanges()->GetHistoPhiBins();          Float_t phimax    = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin    = GetHistogramRanges()->GetHistoPhiMin(); 
  Int_t netabins    = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax    = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin    = GetHistogramRanges()->GetHistoEtaMin();	
  Int_t ssbins      = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax     = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin     = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t nbins       = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nmax      = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nmin      = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  Int_t ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();         Float_t dedxmax   = GetHistogramRanges()->GetHistodEdxMax();         Float_t dedxmin   = GetHistogramRanges()->GetHistodEdxMin();
  Int_t nPoverEbins = GetHistogramRanges()->GetHistoEOverPBins();       Float_t pOverEmax = GetHistogramRanges()->GetHistoEOverPMax();       Float_t pOverEmin = GetHistogramRanges()->GetHistoEOverPMin();
  Int_t tbins       = GetHistogramRanges()->GetHistoTimeBins() ;        Float_t tmax      = GetHistogramRanges()->GetHistoTimeMax();         Float_t tmin      = GetHistogramRanges()->GetHistoTimeMin();
  Int_t nNSigmabins = GetHistogramRanges()->GetHistoNSigmaBins();       Float_t nSigmamax = GetHistogramRanges()->GetHistoNSigmaMax();       Float_t nSigmamin = GetHistogramRanges()->GetHistoNSigmaMin();
  
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  // MC labels, titles, for originator particles
  TString ptypess[] = { "#gamma","hadron?","#pi^{0}","#eta","#gamma->e^{#pm}","e^{#pm}"} ;
  TString pnamess[] = { "Photon","Hadron" ,"Pi0"    ,"Eta" ,"Conversion"     ,"Electron"} ;
  TString ptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}", "#pi^{0}","#eta",
    "e^{#pm}","#gamma->e^{#pm}","hadron?","Anti-N","Anti-P"                    } ;
  
  TString pname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Pi0","Eta","Electron",
    "Conversion", "Hadron", "AntiNeutron","AntiProton"                        } ;
  
  TString pidParticle[] = {"Electron","ChargedHadron","NoPID"} ;
  
  //
  // E/p histograms
  //
  fhEOverPvsE  = new TH2F 
  ("hEOverP_clusE","matched track #it{E/p} vs cluster #it{E}", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsE->SetXTitle("#it{E} (GeV)");
  fhEOverPvsE->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsE);  
  
  fhEOverPvsP  = new TH2F 
  ("hEOverP_TraP","matched track #it{E/p} vs track #it{p}", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsP->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsP->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsP);  
  
  fhEOverPvsECutM02  = new TH2F 
  ("hEOverP_clusE_CutM02",
   "matched track #it{E/p} vs cluster #it{E}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02);
  
  fhEOverPvsPCutM02  = new TH2F 
  ("hEOverP_TraP_CutM02",
   "matched track #it{E/p} vs track #it{p}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02);
  
  fhEOverPvsECutdEdx  = new TH2F 
  ("hEOverP_clusE_CutdEdx","matched track #it{E/p} vs cluster #it{E}, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsECutdEdx->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutdEdx);  
  
  fhEOverPvsPCutdEdx  = new TH2F 
  ("hEOverP_TraP_CutdEdx","matched track #it{E/p} vs track #it{p}, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsPCutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutdEdx);  
  
  fhEOverPvsECutM02AndM20  = new TH2F 
  ("hEOverP_clusE_CutM02AndM20",
   "matched track #it{E/p} vs cluster #it{E}, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02AndM20->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02AndM20->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02AndM20);
  
  fhEOverPvsPCutM02AndM20  = new TH2F 
  ("hEOverP_TraP_CutM02AndM20",
   "matched track #it{E/p} vs track #it{p}, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02AndM20->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02AndM20->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02AndM20);
  
  fhEOverPvsECutM02CutdEdx  = new TH2F 
  ("hEOverP_clusE_CutM02CutdEdx",
   "matched track #it{E/p} vs cluster #it{E}, <d#it{E}/d#it{x}> cut, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02CutdEdx->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02CutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02CutdEdx);
  
  fhEOverPvsPCutM02CutdEdx  = new TH2F 
  ("hEOverP_TraP_CutM02CutdEdx",
   "matched track #it{E/p} vs track #it{p}, <d#it{E}/d#it{x}> cut, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02CutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02CutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02CutdEdx);
  
  fhEOverPvsECutM02AndM20CutdEdx  = new TH2F 
  ("hEOverP_clusE_CutM02AndM20CutdEdx",
   "matched track #it{E/p} vs cluster #it{E}, <d#it{E}/d#it{x}> cut, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut",
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02AndM20CutdEdx->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02AndM20CutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02AndM20CutdEdx);
  
  fhEOverPvsPCutM02AndM20CutdEdx  = new TH2F 
  ("hEOverP_TraP_CutM02AndM20CutdEdx",
   "matched track #it{E/p} vs track #it{p}, <d#it{E}/d#it{x}> cut, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut",
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02AndM20CutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02AndM20CutdEdx->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02AndM20CutdEdx);
  
  fhEOverPvsECutNSigma  = new TH2F 
  ("hEOverP_clusE_CutNSigma","matched track #it{E/p} vs cluster #it{E}, n#sigma cut cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsECutNSigma->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutNSigma);  
  
  fhEOverPvsPCutNSigma  = new TH2F 
  ("hEOverP_TraP_CutNSigma","matched track #it{E/p} vs track #it{p}, n#sigma cut cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsPCutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutNSigma);  
  
  fhEOverPvsECutM02CutNSigma  = new TH2F 
  ("hEOverP_clusE_CutM02CutNSigma",
   "matched track #it{E/p} vs cluster #it{E}, n#sigma cut, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02CutNSigma->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02CutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02CutNSigma);
  
  fhEOverPvsPCutM02CutNSigma  = new TH2F 
  ("hEOverP_TraP_CutM02CutNSigma",
   "matched track #it{E/p} vs track #it{p}, n#sigma cut, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02CutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02CutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02CutNSigma);
  
  fhEOverPvsECutM02AndM20CutNSigma  = new TH2F 
  ("hEOverP_clusE_CutM02AndM20CutNSigma",
   "matched track #it{E/p} vs cluster #it{E}, n#sigma cut, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut",
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02AndM20CutNSigma->SetXTitle("#it{E} (GeV)");
  fhEOverPvsECutM02AndM20CutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsECutM02AndM20CutNSigma);
  
  fhEOverPvsPCutM02AndM20CutNSigma  = new TH2F 
  ("hEOverP_TraP_CutM02AndM20CutNSigma",
   "matched track #it{E/p} vs track #it{p}, n#sigma cut, #sigma_{long}^{2} cut, #sigma_{short}^{2} cut",
   nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02AndM20CutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhEOverPvsPCutM02AndM20CutNSigma->SetYTitle("#it{E/p}");
  outputContainer->Add(fhEOverPvsPCutM02AndM20CutNSigma);
  
   // dE/dX histograms
  
  fhdEdxvsE  = new TH2F 
  ("hdEdx_clusE","matched track <d#it{E}/d#it{x}> vs cluster #it{E} ", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsE->SetXTitle("#it{E} (GeV)");
  fhdEdxvsE->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsE);  
  
  fhdEdxvsP  = new TH2F 
  ("hdEdx_TraP","matched track <d#it{E}/d#it{x}> vs track #it{p} ", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
  fhdEdxvsP->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsP->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsP);  
  
  fhdEdxvsECutM02  = new TH2F 
  ("hdEdx_clusE_CutM02",
   "matched track <d#it{E}/d#it{x}> vs cluster #it{E}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutM02->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutM02->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutM02);
  
  fhdEdxvsPCutM02  = new TH2F 
  ("hdEdx_TraP_CutM02",
   "matched track <d#it{E}/d#it{x}> vs track #it{p}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutM02->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutM02->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutM02);
  
  fhdEdxvsECutM02AndM20  = new TH2F 
  ("hdEdx_clusE_CutM02AndM20",
   "matched track <d#it{E}/d#it{x}> vs cluster #it{E}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutM02AndM20->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutM02AndM20->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutM02AndM20);
  
  fhdEdxvsPCutM02AndM20  = new TH2F 
  ("hdEdx_TraP_CutM02AndM20",
   "matched track <d#it{E}/d#it{x}> vs track #it{p}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutM02AndM20->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutM02AndM20->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutM02AndM20);
  
  fhdEdxvsECutEOverP  = new TH2F 
  ("hdEdx_clusE_CutEOverP","matched track <d#it{E}/d#it{x}> vs cluster #it{E}, cut on #it{E/p}", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutEOverP->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutEOverP->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutEOverP);
  
  fhdEdxvsPCutEOverP  = new TH2F 
  ("hdEdx_TraP_CutEOverP","matched track <d#it{E}/d#it{x}> vs track #it{p}, cut on #it{E/p}", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutEOverP->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutEOverP->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutEOverP);
  
  fhdEdxvsECutNSigma  = new TH2F 
  ("hdEdx_clusE_CutNSigma","matched track <d#it{E}/d#it{x}> vs cluster #it{E}, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutNSigma->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutNSigma);  
  
  fhdEdxvsPCutNSigma  = new TH2F 
  ("hdEdx_TraP_CutNSigma","matched track <d#it{E}/d#it{x}> vs track #it{p}, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
  fhdEdxvsPCutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutNSigma);  
  
  fhdEdxvsECutM02CutNSigma  = new TH2F 
  ("hdEdx_clusE_CutM02CutNSigma",
   "matched track <d#it{E}/d#it{x}> vs cluster #it{E}, #sigma_{long}^{2} cut, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutM02CutNSigma->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutM02CutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutM02CutNSigma);
  
  fhdEdxvsPCutM02CutNSigma  = new TH2F 
  ("hdEdx_TraP_CutM02CutNSigma",
   "matched track <d#it{E}/d#it{x}> vs track #it{p}, #sigma_{long}^{2} cut, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutM02CutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutM02CutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutM02CutNSigma);
  
  fhdEdxvsECutM02AndM20CutNSigma  = new TH2F 
  ("hdEdx_clusE_CutM02AndM20CutNSigma",
   "matched track <d#it{E}/d#it{x}> vs cluster #it{E}, #sigma_{long}^{2} cut, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutM02AndM20CutNSigma->SetXTitle("#it{E} (GeV)");
  fhdEdxvsECutM02AndM20CutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsECutM02AndM20CutNSigma);
  
  fhdEdxvsPCutM02AndM20CutNSigma  = new TH2F 
  ("hdEdx_TraP_CutM02AndM20CutNSigma",
   "matched track <d#it{E}/d#it{x}> vs track #it{p}, #sigma_{long}^{2} cut, n#sigma cut", 
   nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutM02AndM20CutNSigma->SetXTitle("#it{p} (GeV/#it{c})");
  fhdEdxvsPCutM02AndM20CutNSigma->SetYTitle("<d#it{E}/d#it{x}>");
  outputContainer->Add(fhdEdxvsPCutM02AndM20CutNSigma);
  
   // TPC nSigma histograms
  
  fhNSigmavsE  = new TH2F 
  ("hNSigma_clusE","matched track <d#it{E}/d#it{x}> vs cluster #it{E} ", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsE->SetXTitle("#it{E} (GeV)");
  fhNSigmavsE->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsE);  
  
  fhNSigmavsP  = new TH2F 
  ("hNSigma_TraP","matched track n#sigma vs track #it{p} ", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax); 
  fhNSigmavsP->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsP->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsP);  
  
  fhNSigmavsECutM02  = new TH2F 
  ("hNSigma_clusE_CutM02",
   "matched track n#sigma vs cluster #it{E}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutM02->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutM02->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutM02);
  
  fhNSigmavsPCutM02  = new TH2F 
  ("hNSigma_TraP_CutM02",
   "matched track n#sigma vs track #it{p}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsPCutM02->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutM02->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutM02);
  
  fhNSigmavsECutM02AndM20  = new TH2F 
  ("hNSigma_clusE_CutM02AndM20",
   "matched track n#sigma vs cluster #it{E}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutM02AndM20->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutM02AndM20->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutM02AndM20);
  
  fhNSigmavsPCutM02AndM20  = new TH2F 
  ("hNSigma_TraP_CutM02AndM20",
   "matched track n#sigma vs track #it{p}, #sigma_{long}^{2} cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsPCutM02AndM20->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutM02AndM20->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutM02AndM20);
  
  fhNSigmavsECutEOverP  = new TH2F 
  ("hNSigma_clusE_CutEOverP","matched track n#sigma vs cluster #it{E}, cut on #it{E/p}", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutEOverP->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutEOverP->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutEOverP);
  
  fhNSigmavsPCutEOverP  = new TH2F 
  ("hNSigma_TraP_CutEOverP","matched track n#sigma vs track #it{p}, cut on #it{E/p}", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsPCutEOverP->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutEOverP->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutEOverP);
  
  fhNSigmavsECutdEdx  = new TH2F 
  ("hNSigma_clusE_CutdEdx","matched track n#sigma vs cluster #it{E}, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutdEdx->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutdEdx);  
  
  fhNSigmavsPCutdEdx  = new TH2F 
  ("hNSigma_TraP_CutdEdx","matched track n#sigma vs track #it{p}, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax); 
  fhNSigmavsPCutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutdEdx);  
  
  fhNSigmavsECutM02CutdEdx  = new TH2F 
  ("hNSigma_clusE_CutM02CutdEdx",
   "matched track n#sigma vs cluster #it{E}, #sigma_{long}^{2} cut, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutM02CutdEdx->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutM02CutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutM02CutdEdx);
  
  fhNSigmavsPCutM02CutdEdx  = new TH2F 
  ("hNSigma_TraP_CutM02CutdEdx",
   "matched track n#sigma vs track #it{p}, #sigma_{long}^{2} cut, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsPCutM02CutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutM02CutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutM02CutdEdx);
  
  fhNSigmavsECutM02AndM20CutdEdx  = new TH2F 
  ("hNSigma_clusE_CutM02AndM20CutdEdx",
   "matched track n#sigma vs cluster #it{E}, #sigma_{long}^{2} cut, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsECutM02AndM20CutdEdx->SetXTitle("#it{E} (GeV)");
  fhNSigmavsECutM02AndM20CutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsECutM02AndM20CutdEdx);
  
  fhNSigmavsPCutM02AndM20CutdEdx  = new TH2F 
  ("hNSigma_TraP_CutM02AndM20CutdEdx",
   "matched track n#sigma vs track #it{p}, #sigma_{long}^{2} cut, <d#it{E}/d#it{x}> cut", 
   nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
  fhNSigmavsPCutM02AndM20CutdEdx->SetXTitle("#it{p} (GeV/#it{c})");
  fhNSigmavsPCutM02AndM20CutdEdx->SetYTitle("n#sigma");
  outputContainer->Add(fhNSigmavsPCutM02AndM20CutdEdx);
  
  
  if ( IsDataMC() )
  {
    for(Int_t i = 0; i < fNOriginHistograms; i++)
    {
      fhMCdEdxvsE[i]  = new TH2F
      (Form("hdEdx_clusE_MC%s",pname[i].Data()),
       Form("matched track <d#it{E}/d#it{x}> vs cluster #it{E} from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhMCdEdxvsE[i]->SetXTitle("#it{E} (GeV)");
      fhMCdEdxvsE[i]->SetYTitle("<d#it{E}/d#it{x}>");
      outputContainer->Add(fhMCdEdxvsE[i]) ;
      
      fhMCdEdxvsP[i]  = new TH2F
      (Form("hdEdx_TraP_MC%s",pname[i].Data()),
       Form("matched track <d#it{E}/d#it{x}> vs track P from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhMCdEdxvsP[i]->SetXTitle("#it{E} (GeV)");
      fhMCdEdxvsP[i]->SetYTitle("<d#it{E}/d#it{x}>");
      outputContainer->Add(fhMCdEdxvsP[i]) ;
 
      fhMCNSigmavsE[i]  = new TH2F
      (Form("hNSigma_clusE_MC%s",pname[i].Data()),
       Form("matched track n#sigma vs cluster #it{E} from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
      fhMCNSigmavsE[i]->SetXTitle("#it{E} (GeV)");
      fhMCNSigmavsE[i]->SetYTitle("n#sigma");
      outputContainer->Add(fhMCNSigmavsE[i]) ;
      
      fhMCNSigmavsP[i]  = new TH2F
      (Form("hNSigma_TraP_MC%s",pname[i].Data()),
       Form("matched track n#sigma vs track P from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax, nNSigmabins, nSigmamin, nSigmamax);
      fhMCNSigmavsP[i]->SetXTitle("#it{E} (GeV)");
      fhMCNSigmavsP[i]->SetYTitle("n#sigma");
      outputContainer->Add(fhMCNSigmavsP[i]) ;
      
      fhMCEOverPvsE[i]  = new TH2F
      (Form("hEOverP_clusE_MC%s",pname[i].Data()),
       Form("matched track #it{E/p} vs cluster #it{E} from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhMCEOverPvsE[i]->SetXTitle("#it{E} (GeV)");
      fhMCEOverPvsE[i]->SetYTitle("<d#it{E}/d#it{x}>");
      outputContainer->Add(fhMCEOverPvsE[i]) ;
      
      fhMCEOverPvsP[i]  = new TH2F
      (Form("hEOverP_TraP_MC%s",pname[i].Data()),
       Form("matched track #it{E/p} vs track P from %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhMCEOverPvsP[i]->SetXTitle("#it{E} (GeV)");
      fhMCEOverPvsP[i]->SetYTitle("<d#it{E}/d#it{x}>");
      outputContainer->Add(fhMCEOverPvsP[i]) ;
      
      for(Int_t j = 0; j < 2; j++)
      {
        fhMCEOverPvsEAfterCuts[i][j]  = new TH2F
        (Form("hEOverP_clusE_%s_MC%s",pidParticle[j].Data(), pname[i].Data()),
         Form("matched track #it{E/p} vs cluster #it{E}, id %s from MC %s",pidParticle[j].Data(),ptype[i].Data()),
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
        fhMCEOverPvsEAfterCuts[i][j]->SetXTitle("#it{E} (GeV)");
        fhMCEOverPvsEAfterCuts[i][j]->SetYTitle("<d#it{E}/d#it{x}>");
        outputContainer->Add(fhMCEOverPvsEAfterCuts[i][j]) ;
        
        fhMCEOverPvsPAfterCuts[i][j]  = new TH2F
        (Form("hEOverP_TraP_%s_MC%s",pidParticle[j].Data(),pname[i].Data()),
         Form("matched track #it{E/p} vs track P, id %s, from MC %s",pidParticle[j].Data(),ptype[i].Data()),
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
        fhMCEOverPvsPAfterCuts[i][j]->SetXTitle("#it{E} (GeV)");
        fhMCEOverPvsPAfterCuts[i][j]->SetYTitle("<d#it{E}/d#it{x}>");
        outputContainer->Add(fhMCEOverPvsPAfterCuts[i][j]) ;
      }
      
    } // MC particle loop
  } // Is MC
  
  
  // Matching residuals
  TString sCharge[] = {"Positive","Negative"};
  for(Int_t indexPID = 0; indexPID < 3; indexPID++)
  {
    for(Int_t ich = 0; ich < 2; ich++)
    { 
      fhDEtavsE[indexPID][ich]  = new TH2F
      (Form("hDEta_clusE_%s_%s",sCharge[ich].Data(),pidParticle[indexPID].Data()),
       Form("#Delta #eta of cluster - %s track vs #it{E}_{cluster}, ID: %s",
            sCharge[ich].Data(),pidParticle[indexPID].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhDEtavsE[indexPID][ich]->SetYTitle("#Delta #eta");
      fhDEtavsE[indexPID][ich]->SetXTitle("#it{E}_{cluster} (GeV)");
      outputContainer->Add(fhDEtavsE[indexPID][ich]) ;

      fhDPhivsE[indexPID][ich]  = new TH2F
      (Form("hDPhi_clusE_%s_%s",sCharge[ich].Data(),pidParticle[indexPID].Data()),
       Form("#Delta #varphi of cluster - %s track vs #it{E}_{cluster}, ID: %s",
            sCharge[ich].Data(),pidParticle[indexPID].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhDPhivsE[indexPID][ich]->SetYTitle("#Delta #varphi (rad)");
      fhDPhivsE[indexPID][ich]->SetXTitle("#it{E}_{cluster} (GeV)");
      outputContainer->Add(fhDPhivsE[indexPID][ich]) ;

      fhDEtavsP[indexPID][ich]  = new TH2F
      (Form("hDEta_TraP_%s_%s",sCharge[ich].Data(),pidParticle[indexPID].Data()),
       Form("#Delta #eta of cluster - %s track vs #it{p}_{track}, ID: %s",
            sCharge[ich].Data(),pidParticle[indexPID].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhDEtavsP[indexPID][ich]->SetYTitle("#Delta #eta");
      fhDEtavsP[indexPID][ich]->SetXTitle("#it{p}_{track} (GeV/#it{c})");
      outputContainer->Add(fhDEtavsP[indexPID][ich]) ;
      
      fhDPhivsP[indexPID][ich]  = new TH2F
      (Form("hDPhi_TraP_%s_%s",sCharge[ich].Data(),pidParticle[indexPID].Data()),
       Form("#Delta #varphi of cluster - %s track vs #it{p}_{track}, ID: %s",
            sCharge[ich].Data(),pidParticle[indexPID].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhDPhivsP[indexPID][ich]->SetYTitle("#Delta #varphi (rad)");
      fhDPhivsP[indexPID][ich]->SetXTitle("#it{p}_{track} (GeV/#it{c})");
      outputContainer->Add(fhDPhivsP[indexPID][ich]) ;
      
      fhDEtaDPhi[indexPID][ich]  = new TH2F
      (Form("hDEtaDPhi_%s_%s",sCharge[ich].Data(),pidParticle[indexPID].Data()),
       Form("#Delta #eta vs #Delta #varphi of cluster - %s track, ID: %s",
            sCharge[ich].Data(),pidParticle[indexPID].Data()),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhDEtaDPhi[indexPID][ich]->SetYTitle("#Delta #varphi (rad)");
      fhDEtaDPhi[indexPID][ich]->SetXTitle("#Delta #eta");
      outputContainer->Add(fhDEtaDPhi[indexPID][ich]) ;
    }
  }
  
      
      
  if ( fFillWeightHistograms )
  {
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy, for selected electrons",
                                     nptbins,ptmin,ptmax, 100,0,1.); 
    fhECellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV)");
    fhECellClusterRatio->SetYTitle("#it{E}_{cell i}/#it{E}_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy, for selected electrons",
                                        nptbins,ptmin,ptmax, 100,-10,0); 
    fhECellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV)");
    fhECellClusterLogRatio->SetYTitle("Log (#it{E}_{max cell}/#it{E}_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy, for selected electrons",
                                        nptbins,ptmin,ptmax, 100,0,1.); 
    fhEMaxCellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV)");
    fhEMaxCellClusterRatio->SetYTitle("#it{E}_{max cell}/#it{E}_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy, for selected electrons",
                                           nptbins,ptmin,ptmax, 100,-10,0); 
    fhEMaxCellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV)");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (#it{E}_{max cell}/#it{E}_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #sigma_{long}_{2} vs E, w0 = %1.1f, for selected electrons",1+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLambda0ForW0[iw]->SetXTitle("#it{E}_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#sigma_{long}_{2}");
      outputContainer->Add(fhLambda0ForW0[iw]); 
      
      //        fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #sigma_{short}_{2} vs E, w0 = %1.1f, for selected electrons",1+0.5*iw),
      //                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      //        fhLambda1ForW0[iw]->SetXTitle("#it{E}_{cluster}");
      //        fhLambda1ForW0[iw]->SetYTitle("#sigma_{short}_{2}");
      //        outputContainer->Add(fhLambda1ForW0[iw]);
    }
  }  // weight
  
  for(Int_t pidIndex = 0; pidIndex < 2; pidIndex++)
  {
    // Shower shape
    if ( fFillSSHistograms )
    {
      fhLam0E[pidIndex]  = new TH2F (Form("h%sLam0E",pidParticle[pidIndex].Data()),
                                     Form("%s: #sigma_{long}^{2} vs E",pidParticle[pidIndex].Data()), 
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam0E[pidIndex]->SetYTitle("#sigma_{long}^{2}");
      fhLam0E[pidIndex]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam0E[pidIndex]);  
      
      fhLam1E[pidIndex]  = new TH2F (Form("h%sLam1E",pidParticle[pidIndex].Data()),
                                     Form("%s: #sigma_{short}^{2} vs E",pidParticle[pidIndex].Data()), 
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam1E[pidIndex]->SetYTitle("#sigma_{short}^{2}");
      fhLam1E[pidIndex]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam1E[pidIndex]);  
      
      if ( !fFillOnlySimpleSSHisto )
      {
        fhDispE[pidIndex]  = new TH2F (Form("h%sDispE",pidParticle[pidIndex].Data()),
                                       Form("%s: dispersion^{2} vs E",pidParticle[pidIndex].Data()), 
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhDispE[pidIndex]->SetYTitle("#it{D}^{2}");
        fhDispE[pidIndex]->SetXTitle("#it{E} (GeV) ");
        outputContainer->Add(fhDispE[pidIndex]);
      }
      
      if ( GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >=0 )
      {
        fhLam0ETRD[pidIndex]  = new TH2F (Form("h%sLam0ETRD",pidParticle[pidIndex].Data()),
                                          Form("%s: #sigma_{long}^{2} vs E, EMCAL SM covered by TRD",pidParticle[pidIndex].Data()), 
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam0ETRD[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        fhLam0ETRD[pidIndex]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhLam0ETRD[pidIndex]);  
        
        fhLam1ETRD[pidIndex]  = new TH2F (Form("h%sLam1ETRD",pidParticle[pidIndex].Data()),
                                          Form("%s: #sigma_{short}^{2} vs E, EMCAL SM covered by TRD",pidParticle[pidIndex].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam1ETRD[pidIndex]->SetYTitle("#sigma_{short}^{2}");
        fhLam1ETRD[pidIndex]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhLam1ETRD[pidIndex]);  
        
        if ( !fFillOnlySimpleSSHisto )
        {
          fhDispETRD[pidIndex]  = new TH2F (Form("h%sDispETRD",pidParticle[pidIndex].Data()),
                                            Form("%s: dispersion^{2} vs #it{E}, EMCal SM covered by TRD",pidParticle[pidIndex].Data()), 
                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhDispETRD[pidIndex]->SetYTitle("Dispersion^{2}");
          fhDispETRD[pidIndex]->SetXTitle("#it{E} (GeV) ");
          outputContainer->Add(fhDispETRD[pidIndex]);  
        }
      } 
      
      if ( !fFillOnlySimpleSSHisto )
      {
        fhNCellsLam0LowE[pidIndex]  = new TH2F (Form("h%sNCellsLam0LowE",pidParticle[pidIndex].Data()),
                                                Form("%s: N_{cells} in cluster vs #sigma_{long}^{2}, #it{E} < 2 GeV",pidParticle[pidIndex].Data()),
                                                nbins,nmin, nmax, ssbins,ssmin,ssmax); 
        fhNCellsLam0LowE[pidIndex]->SetXTitle("N_{cells}");
        fhNCellsLam0LowE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        outputContainer->Add(fhNCellsLam0LowE[pidIndex]);  
        
        fhNCellsLam0HighE[pidIndex]  = new TH2F (Form("h%sNCellsLam0HighE",pidParticle[pidIndex].Data()),
                                                 Form("%s: N_{cells} in cluster vs #sigma_{long}^{2}, #it{E} > 2 GeV",pidParticle[pidIndex].Data()), 
                                                 nbins,nmin, nmax, ssbins,ssmin,ssmax); 
        fhNCellsLam0HighE[pidIndex]->SetXTitle("#it{N}_{cells}");
        fhNCellsLam0HighE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        outputContainer->Add(fhNCellsLam0HighE[pidIndex]);  
        
        
        fhEtaLam0LowE[pidIndex]  = new TH2F (Form("h%sEtaLam0LowE",pidParticle[pidIndex].Data()),
                                             Form("%s: #eta vs #sigma_{long}^{2}, #it{E} < 2 GeV",pidParticle[pidIndex].Data()), 
                                             netabins,etamin,etamax, ssbins,ssmin,ssmax); 
        fhEtaLam0LowE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        fhEtaLam0LowE[pidIndex]->SetXTitle("#eta");
        outputContainer->Add(fhEtaLam0LowE[pidIndex]);  
        
        fhPhiLam0LowE[pidIndex]  = new TH2F (Form("h%sPhiLam0LowE",pidParticle[pidIndex].Data()),
                                             Form("%s: #varphi vs #sigma_{long}^{2}, #it{E} < 2 GeV",pidParticle[pidIndex].Data()), 
                                             nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
        fhPhiLam0LowE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        fhPhiLam0LowE[pidIndex]->SetXTitle("#varphi (rad)");
        outputContainer->Add(fhPhiLam0LowE[pidIndex]);  
        
        fhEtaLam0HighE[pidIndex]  = new TH2F (Form("h%sEtaLam0HighE",pidParticle[pidIndex].Data()),
                                              Form("%s: #eta vs #sigma_{long}^{2}, #it{E} > 2 GeV",pidParticle[pidIndex].Data()),
                                              netabins,etamin,etamax, ssbins,ssmin,ssmax); 
        fhEtaLam0HighE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        fhEtaLam0HighE[pidIndex]->SetXTitle("#eta");
        outputContainer->Add(fhEtaLam0HighE[pidIndex]);  
        
        fhPhiLam0HighE[pidIndex]  = new TH2F (Form("h%sPhiLam0HighE",pidParticle[pidIndex].Data()),
                                              Form("%s: #varphi vs #sigma_{long}^{2}, #it{E} > 2 GeV",pidParticle[pidIndex].Data()), 
                                              nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
        fhPhiLam0HighE[pidIndex]->SetYTitle("#sigma_{long}^{2}");
        fhPhiLam0HighE[pidIndex]->SetXTitle("#varphi (rad)");
        outputContainer->Add(fhPhiLam0HighE[pidIndex]);  
        
        if ( GetCalorimeter() == kEMCAL )
        {
          fhDispEtaE[pidIndex]  = new TH2F (Form("h%sDispEtaE",pidParticle[pidIndex].Data()),
                                            Form("%s: #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",pidParticle[pidIndex].Data()),  
                                            nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhDispEtaE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhDispEtaE[pidIndex]->SetYTitle("#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhDispEtaE[pidIndex]);     
          
          fhDispPhiE[pidIndex]  = new TH2F (Form("h%sDispPhiE",pidParticle[pidIndex].Data()),
                                            Form("%s: #sigma^{2}_{#varphi #varphi} = #Sigma w_{i}(#varphi_{i} - <#varphi>)^{2} / #Sigma w_{i} vs E",pidParticle[pidIndex].Data()),  
                                            nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhDispPhiE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhDispPhiE[pidIndex]->SetYTitle("#sigma^{2}_{#varphi #varphi}");
          outputContainer->Add(fhDispPhiE[pidIndex]);  
          
          fhSumEtaE[pidIndex]  = new TH2F (Form("h%sSumEtaE",pidParticle[pidIndex].Data()),
                                           Form("%s: #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i})^{2} / #Sigma w_{i} - <#eta>^{2} vs E",pidParticle[pidIndex].Data()),  
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhSumEtaE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhSumEtaE[pidIndex]->SetYTitle("#delta^{2}_{#eta #eta}");
          outputContainer->Add(fhSumEtaE[pidIndex]);     
          
          fhSumPhiE[pidIndex]  = new TH2F (Form("h%sSumPhiE",pidParticle[pidIndex].Data()),
                                           Form("%s: #sigma^{2}_{#varphi #varphi} = #Sigma w_{i}(#varphi_{i})^{2}/ #Sigma w_{i} - <#varphi>^{2} vs E",pidParticle[pidIndex].Data()),  
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhSumPhiE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhSumPhiE[pidIndex]->SetYTitle("#delta^{2}_{#varphi #varphi}");
          outputContainer->Add(fhSumPhiE[pidIndex]);  
          
          fhSumEtaPhiE[pidIndex]  = new TH2F (Form("h%sSumEtaPhiE",pidParticle[pidIndex].Data()),
                                              Form("%s: #delta^{2}_{#eta #varphi} = #Sigma w_{i}(#varphi_{i} #eta_{i} ) / #Sigma w_{i} - <#varphi><#eta> vs E",pidParticle[pidIndex].Data()),  
                                              nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax); 
          fhSumEtaPhiE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhSumEtaPhiE[pidIndex]->SetYTitle("#delta^{2}_{#eta #varphi}");
          outputContainer->Add(fhSumEtaPhiE[pidIndex]);
          
          fhDispEtaPhiDiffE[pidIndex]  = new TH2F (Form("h%sDispEtaPhiDiffE",pidParticle[pidIndex].Data()),
                                                   Form("%s: #sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta} vs E",pidParticle[pidIndex].Data()), 
                                                   nptbins,ptmin,ptmax,200, -10,10); 
          fhDispEtaPhiDiffE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhDispEtaPhiDiffE[pidIndex]->SetYTitle("#sigma^{2}_{#varphi #varphi}-#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhDispEtaPhiDiffE[pidIndex]);    
          
          fhSphericityE[pidIndex]  = new TH2F (Form("h%sSphericityE",pidParticle[pidIndex].Data()),
                                               Form("%s: (#sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#varphi #varphi}) vs E",pidParticle[pidIndex].Data()),  
                                               nptbins,ptmin,ptmax, 200, -1,1); 
          fhSphericityE[pidIndex]->SetXTitle("#it{E} (GeV)");
          fhSphericityE[pidIndex]->SetYTitle("s = (#sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#varphi #varphi})");
          outputContainer->Add(fhSphericityE[pidIndex]);
          
          Int_t bin[] = {0,2,4,6,10,1000};
          for(Int_t i = 0; i < 5; i++)
          {
            fhDispEtaDispPhiEBin[pidIndex][i] = new TH2F (Form("h%sDispEtaDispPhi_EBin%d",pidParticle[pidIndex].Data(),i),
                                                          Form("%s: #sigma^{2}_{#varphi #varphi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pidParticle[pidIndex].Data(),bin[i],bin[i+1]), 
                                                          ssbins,ssmin,ssmax , ssbins,ssmin,ssmax); 
            fhDispEtaDispPhiEBin[pidIndex][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
            fhDispEtaDispPhiEBin[pidIndex][i]->SetYTitle("#sigma^{2}_{#varphi #varphi}");
            outputContainer->Add(fhDispEtaDispPhiEBin[pidIndex][i]); 
          }
        }
      }
    } // Shower shape
        
    if ( IsDataMC() )
    {
      if ( fFillSSHistograms )
      {
        for(Int_t i = 0; i < 6; i++)
        { 
          fhMCELambda0[pidIndex][i]  = new TH2F(Form("h%sELambda0_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                Form("%s like cluster from %s : #it{E} vs #sigma_{long}^{2}",pidParticle[pidIndex].Data(),ptypess[i].Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhMCELambda0[pidIndex][i]->SetYTitle("#sigma_{long}^{2}");
          fhMCELambda0[pidIndex][i]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhMCELambda0[pidIndex][i]) ; 
  
          fhMCELambda1[pidIndex][i]  = new TH2F(Form("h%sELambda1_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                Form("%s like cluster from %s : #it{E} vs #sigma_{short}^{2}",pidParticle[pidIndex].Data(),ptypess[i].Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhMCELambda1[pidIndex][i]->SetYTitle("#sigma_{short}^{2}");
          fhMCELambda1[pidIndex][i]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhMCELambda1[pidIndex][i]) ; 
          
          if ( GetCalorimeter()==kEMCAL && !fFillOnlySimpleSSHisto )
          {
            fhMCEDispEta[pidIndex][i]  = new TH2F (Form("h%sEDispEtaE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                   Form("cluster from %s : %s like, #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs #it{E}",
                                                        ptypess[i].Data(),pidParticle[pidIndex].Data()),
                                                   nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
            fhMCEDispEta[pidIndex][i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispEta[pidIndex][i]->SetYTitle("#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEta[pidIndex][i]);     
            
            fhMCEDispPhi[pidIndex][i]  = new TH2F (Form("h%sEDispPhiE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                   Form("cluster from %s : %s like, #sigma^{2}_{#varphi #varphi} = #Sigma w_{i}(#varphi_{i} - <#varphi>)^{2} / #Sigma w_{i} vs #it{E}",
                                                        ptypess[i].Data(),pidParticle[pidIndex].Data()),
                                                   nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
            fhMCEDispPhi[pidIndex][i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispPhi[pidIndex][i]->SetYTitle("#sigma^{2}_{#varphi #varphi}");
            outputContainer->Add(fhMCEDispPhi[pidIndex][i]);  
            
            fhMCESumEtaPhi[pidIndex][i]  = new TH2F (Form("h%sESumEtaPhiE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                     Form("cluster from %s : %s like, #delta^{2}_{#eta #varphi} = #Sigma w_{i}(#varphi_{i} #eta_{i} ) / #Sigma w_{i} - <#varphi><#eta> vs E",
                                                          ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                     nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax); 
            fhMCESumEtaPhi[pidIndex][i]->SetXTitle("#it{E} (GeV)");
            fhMCESumEtaPhi[pidIndex][i]->SetYTitle("#delta^{2}_{#eta #varphi}");
            outputContainer->Add(fhMCESumEtaPhi[pidIndex][i]);
            
            fhMCEDispEtaPhiDiff[pidIndex][i]  = new TH2F (Form("h%sEDispEtaPhiDiffE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                          Form("cluster from %s : %s like, #sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta} vs #it{E}",
                                                               ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                          nptbins,ptmin,ptmax,200,-10,10); 
            fhMCEDispEtaPhiDiff[pidIndex][i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispEtaPhiDiff[pidIndex][i]->SetYTitle("#sigma^{2}_{#varphi #varphi}-#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEtaPhiDiff[pidIndex][i]);    
            
            fhMCESphericity[pidIndex][i]  = new TH2F (Form("h%sESphericity_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                      Form("cluster from %s : %s like, (#sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#varphi #varphi}) vs #it{E}",
                                                           ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                      nptbins,ptmin,ptmax, 200,-1,1); 
            fhMCESphericity[pidIndex][i]->SetXTitle("#it{E} (GeV)");
            fhMCESphericity[pidIndex][i]->SetYTitle("#it{s} = (#sigma^{2}_{#varphi #varphi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#varphi #varphi})");
            outputContainer->Add(fhMCESphericity[pidIndex][i]);
          }
          
        }// loop    
      }   
    }
    
    //if ( IsCaloPIDOn() && pidIndex > 0 ) continue;
    
    fhNCellsE[pidIndex]  = new TH2F (Form("h%sNCellsE",pidParticle[pidIndex].Data()),
                                     Form("N cells in %s cluster vs #it{E}",pidParticle[pidIndex].Data()),
                                     nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhNCellsE[pidIndex]->SetXTitle("#it{E} (GeV)");
    fhNCellsE[pidIndex]->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhNCellsE[pidIndex]);  
    
    fhNLME[pidIndex]  = new TH2F (Form("h%sNLME",pidParticle[pidIndex].Data()),
                                     Form("NLM in %s cluster vs #it{E}",pidParticle[pidIndex].Data()),
                                     nptbins,ptmin,ptmax, 10,0,10);
    fhNLME[pidIndex]->SetXTitle("#it{E} (GeV)");
    fhNLME[pidIndex]->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhNLME[pidIndex]);
    
    fhTimeE[pidIndex] = new TH2F(Form("h%sTimeE",pidParticle[pidIndex].Data()),
                                 Form("Time in %s cluster vs #it{E}",pidParticle[pidIndex].Data())
                                 ,nptbins,ptmin,ptmax, tbins,tmin,tmax);
    fhTimeE[pidIndex]->SetXTitle("#it{E} (GeV)");
    fhTimeE[pidIndex]->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhTimeE[pidIndex]);  
    
    fhMaxCellDiffClusterE[pidIndex]  = new TH2F (Form("h%sMaxCellDiffClusterE",pidParticle[pidIndex].Data()),
                                                 Form("%s: energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",pidParticle[pidIndex].Data()),
                                                 nptbins,ptmin,ptmax, 500,0,1.); 
    fhMaxCellDiffClusterE[pidIndex]->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhMaxCellDiffClusterE[pidIndex]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
    outputContainer->Add(fhMaxCellDiffClusterE[pidIndex]);  
    
    fhE[pidIndex]  = new TH1F(Form("h%sE",pidParticle[pidIndex].Data()),
                              Form("Number of %s over calorimeter vs energy",pidParticle[pidIndex].Data()),
                              nptbins,ptmin,ptmax); 
    fhE[pidIndex]->SetYTitle("#it{N}");
    fhE[pidIndex]->SetXTitle("#it{E}(GeV)");
    outputContainer->Add(fhE[pidIndex]) ;   
    
    fhPt[pidIndex]  = new TH1F(Form("h%sPt",pidParticle[pidIndex].Data()),
                               Form("Number of %s over calorimeter vs #it{p}_{T}",pidParticle[pidIndex].Data()),
                               nptbins,ptmin,ptmax); 
    fhPt[pidIndex]->SetYTitle("#it{N}");
    fhPt[pidIndex]->SetXTitle("#it{p} (GeV/#it{c})");
    outputContainer->Add(fhPt[pidIndex]) ; 
    
    fhPhi[pidIndex]  = new TH2F(Form("h%sPhi",pidParticle[pidIndex].Data()),
                                Form("%s: #varphi vs #it{p}_{T}",pidParticle[pidIndex].Data()),
                                nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhi[pidIndex]->SetYTitle("#varphi (rad)");
    fhPhi[pidIndex]->SetXTitle("#it{p} (GeV/#it{c})");
    outputContainer->Add(fhPhi[pidIndex]) ; 
    
    fhEta[pidIndex]  = new TH2F(Form("h%sEta",pidParticle[pidIndex].Data()),
                                Form("%s: #eta vs #it{p}_{T}",pidParticle[pidIndex].Data()),
                                nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEta[pidIndex]->SetYTitle("#eta");
    fhEta[pidIndex]->SetXTitle("#it{p} (GeV/#it{c})");
    outputContainer->Add(fhEta[pidIndex]) ;
    
    fhEtaPhi[pidIndex]  = new TH2F(Form("h%sEtaPhi",pidParticle[pidIndex].Data()),
                                   Form("%s: #eta vs #varphi",pidParticle[pidIndex].Data()),
                                   netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi[pidIndex]->SetYTitle("#varphi (rad)");
    fhEtaPhi[pidIndex]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi[pidIndex]) ;
    
    if ( GetMinPt() < 0.5 )
    {
      fhEtaPhi05[pidIndex]  = new TH2F(Form("h%sEtaPhi05",pidParticle[pidIndex].Data()),
                                       Form("%s: #eta vs #varphi, #it{E} > 0.5 GeV",pidParticle[pidIndex].Data()),
                                       netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05[pidIndex]->SetYTitle("#varphi (rad)");
      fhEtaPhi05[pidIndex]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05[pidIndex]) ;
    }
    
    
    if ( IsDataMC() )
    {      
      for(Int_t i = 0; i < fNOriginHistograms; i++)
      { 
        fhMCE[pidIndex][i]  = new TH1F(Form("h%sE_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                       Form("%s like cluster from %s : E ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                       nptbins,ptmin,ptmax); 
        fhMCE[pidIndex][i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCE[pidIndex][i]) ; 
        
        fhMCPt[pidIndex][i]  = new TH1F(Form("h%sPt_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                        Form("%s like cluster from %s : #it{p}_{T} ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                        nptbins,ptmin,ptmax); 
        fhMCPt[pidIndex][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCPt[pidIndex][i]) ;
        
        fhMCEta[pidIndex][i]  = new TH2F(Form("h%sEta_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like cluster from %s : #eta ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                         nptbins,ptmin,ptmax,netabins,etamin,etamax); 
        fhMCEta[pidIndex][i]->SetYTitle("#eta");
        fhMCEta[pidIndex][i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCEta[pidIndex][i]) ;
        
        fhMCPhi[pidIndex][i]  = new TH2F(Form("h%sPhi_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like cluster from %s : #varphi ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                         nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
        fhMCPhi[pidIndex][i]->SetYTitle("#varphi (rad)");
        fhMCPhi[pidIndex][i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCPhi[pidIndex][i]) ;
        
        
        fhMCDeltaE[pidIndex][i]  = new TH2F (Form("h%sDeltaE_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                             Form("%s like MC - Reco E from %s",pidParticle[pidIndex].Data(),pname[i].Data()), 
                                             nptbins,ptmin,ptmax, 200,-50,50); 
        fhMCDeltaE[pidIndex][i]->SetXTitle("#Delta #it{E} (GeV)");
        outputContainer->Add(fhMCDeltaE[pidIndex][i]);
        
        fhMC2E[pidIndex][i]  = new TH2F (Form("h%s2E_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like E distribution, reconstructed vs generated from %s",pidParticle[pidIndex].Data(),pname[i].Data()), 
                                         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMC2E[pidIndex][i]->SetXTitle("#it{E}_{rec} (GeV)");
        fhMC2E[pidIndex][i]->SetYTitle("#it{E}_{gen} (GeV)");
        outputContainer->Add(fhMC2E[pidIndex][i]);          
        
      }
    } // MC
  }// pid Index
  
  
  if ( fFillSSHistograms )
  {
    if ( IsDataMC() )
    {
      if ( !IsEmbedingAnalysisOn() )
      {
        fhMCElectronELambda0NoOverlap  = new TH2F("hELambda0_MCElectron_NoOverlap",
                                                  "cluster from Electron : E vs #sigma_{long}^{2}",
                                                  nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0NoOverlap->SetYTitle("#sigma_{long}^{2}");
        fhMCElectronELambda0NoOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCElectronELambda0NoOverlap) ; 
        
        fhMCElectronELambda0TwoOverlap  = new TH2F("hELambda0_MCElectron_TwoOverlap",
                                                   "cluster from Electron : E vs #sigma_{long}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0TwoOverlap->SetYTitle("#sigma_{long}^{2}");
        fhMCElectronELambda0TwoOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCElectronELambda0TwoOverlap) ; 
        
        fhMCElectronELambda0NOverlap  = new TH2F("hELambda0_MCElectron_NOverlap",
                                                 "cluster from Electron : E vs #sigma_{long}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0NOverlap->SetYTitle("#sigma_{long}^{2}");
        fhMCElectronELambda0NOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCElectronELambda0NOverlap) ;
      } // No embedding
      
      // Fill histograms to check shape of embedded clusters
      if ( IsEmbedingAnalysisOn() )
      {
        fhEmbeddedSignalFractionEnergy  = new TH2F("hEmbeddedSignal_FractionEnergy",
                                                   "Energy Fraction of embedded signal versus cluster energy",
                                                   nptbins,ptmin,ptmax,100,0.,1.); 
        fhEmbeddedSignalFractionEnergy->SetYTitle("Fraction");
        fhEmbeddedSignalFractionEnergy->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbeddedSignalFractionEnergy) ; 
        
        fhEmbedElectronELambda0FullSignal  = new TH2F("hELambda0_EmbedElectron_FullSignal",
                                                      "cluster from Electron embedded with more than 90% energy in cluster : E vs #sigma_{long}^{2}",
                                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0FullSignal->SetYTitle("#sigma_{long}^{2}");
        fhEmbedElectronELambda0FullSignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0FullSignal) ; 
        
        fhEmbedElectronELambda0MostlySignal  = new TH2F("hELambda0_EmbedElectron_MostlySignal",
                                                        "cluster from Electron embedded with 50% to 90% energy in cluster : E vs #sigma_{long}^{2}",
                                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0MostlySignal->SetYTitle("#sigma_{long}^{2}");
        fhEmbedElectronELambda0MostlySignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0MostlySignal) ; 
        
        fhEmbedElectronELambda0MostlyBkg  = new TH2F("hELambda0_EmbedElectron_MostlyBkg",
                                                     "cluster from Electron embedded with 10% to 50% energy in cluster : E vs #sigma_{long}^{2}",
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0MostlyBkg->SetYTitle("#sigma_{long}^{2}");
        fhEmbedElectronELambda0MostlyBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0MostlyBkg) ; 
        
        fhEmbedElectronELambda0FullBkg  = new TH2F("hELambda0_EmbedElectron_FullBkg",
                                                   "cluster from Electronm embedded with 0% to 10% energy in cluster : E vs #sigma_{long}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0FullBkg->SetYTitle("#sigma_{long}^{2}");
        fhEmbedElectronELambda0FullBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0FullBkg) ;
      } // Embedded histograms
    } // Histos with MC
  } // Fill SS MC histograms
    
  return outputContainer ;
}

//_________________________
/// Init. Check if requested calorimeter is on, if not, abort
//_________________________
void AliAnaElectron::Init()
{
  if       ( GetCalorimeter() == kPHOS  && !GetReader()->IsPHOSSwitchedOn()  && NewOutputAOD() )
    AliFatal("STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!");
  else  if ( GetCalorimeter() == kEMCAL && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD() )
    AliFatal("STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!");
}

//___________________________________
/// Initialize the parameters of the analysis with default values.
//___________________________________
void AliAnaElectron::InitParameters()
{
  AddToHistogramsName("AnaElectron_");
  
  fMinDist     = 2.;
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  
  fNCellsCut   = 1;
  
//  fM20Min = 0.00;
//  fM20Max = 0.30;
  
  fM02Min = 0.05;
  fM02Max = 0.35;
  
//fdEdxMin     = 75;//76. for LHC11a, but for LHC11c pass1 56.                
//fdEdxMax     = 95;//85. for LHC11a, but for LHC11c pass1 64.   

  fEOverPMin   = 0.9; // 0.8 for LHC11a, but for LHC11c pass1 0.9                  
  fEOverPMax   = 1.2; // for LHC11a and LHC11c pass1
  
  fNSigmaMax   = 3;
  fNSigmaMin   =-1;
  
  fNSigmaMaxHad=-4;
  fNSigmaMinHad=-10;
  
}

//_________________________________________
/// Do photon analysis selecting electron clusters (or charged non electron) and fill aods.
//_________________________________________
void  AliAnaElectron::MakeAnalysisFillAOD()
{
  // Get the vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  // Select the Calorimeter of the photon
  TObjArray * pl = 0x0; 
  if      ( GetCalorimeter() == kPHOS  ) pl = GetPHOSClusters ();
  else if ( GetCalorimeter() == kEMCAL ) pl = GetEMCALClusters();
  
  if ( !pl )
  {
    AliWarning(Form("TObjArray with %s clusters is NULL!",GetCalorimeterString().Data()));
    return;
  }
  Int_t dataType = GetReader()->GetDataType() ;
  
  AliPIDResponse *pidResponse = 0;
  
  if      ( dataType == AliCaloTrackReader::kESD ) 
    pidResponse = (dynamic_cast<AliESDInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetPIDResponse();
  else if ( dataType == AliCaloTrackReader::kAOD ) 
    pidResponse = (dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetPIDResponse();
  
  if ( !pidResponse )
  {
    AliFatal("AliPIDResponse not available, did you initialize the task?");
    return; // not needed, coverity ...
  }
  
  //Init arrays, variables, get number of clusters
  Int_t nCaloClusters = pl->GetEntriesFast();
  //List to be used in conversion analysis, to tag the cluster as candidate for conversion
  
  AliDebug(1,Form("Input %s cluster entries %d", GetCalorimeterString().Data(), nCaloClusters));
  
  //----------------------------------------------------
  // Fill AOD with PHOS/EMCAL AliCaloTrackParticle objects
  //----------------------------------------------------
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));	
    //printf("calo %d, %f\n",icalo,cluE);
    
    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ; 
    if ( GetMixedEvent() )
    {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      //Get the vertex and check it is not too large in z
      if ( TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut() ) continue;
    }
    
    //Cluster selection, not charged, with photon id and in fiducial cut	  
    if ( GetReader()->GetDataType() != AliCaloTrackReader::kMC )
    {
      calo->GetMomentum(fMomentum,GetVertex(evtIndex)) ;
    }//Assume that come from vertex in straight line
    else
    {
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(fMomentum,vertex) ;
    }
    
    //--------------------------------------
    // Cluster selection
    //--------------------------------------
    AliVCaloCells* cells = 0;
    if ( GetCalorimeter() == kEMCAL ) cells = GetEMCALCells();
    else                              cells = GetPHOSCells();
    
    Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells); // NLM
    if ( !ClusterSelected(calo,nMaxima) ) continue;
    
    // Select only clusters with MC signal and data background
    //
    if ( SelectEmbededSignal() && IsDataMC() )
    {
      if ( calo->GetNLabels() == 0 || calo->GetLabel() < 0 ) continue;
    }
    
    //-------------------------------------
    // PID selection
    //-------------------------------------
    
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(calo, GetReader()->GetInputEvent());

    if ( !track )
    {
      AliWarning("Null track");
      continue;
    }
    
    Float_t cluE = calo ->E();
    Float_t traP = track->P();
    
    if ( traP < 0.1 )
    {
      AliWarning(Form("Too Low P track %f GeV/#it{c}, continue\n",traP));
      continue;
    }
    
    //
    // Get Identification cuts
    // 
    Float_t dEdx   = track->GetTPCsignal();
    Float_t eOverp = cluE/traP;
    
    Float_t m02    = calo->GetM02();
    Float_t m20    = calo->GetM20();
    
    Double_t nSigma = pidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    
    AliDebug(1,Form("Cluster E %2.2f, Track P %2.2f, E/P %2.2f, dEdx %2.2f, n sigma %2.2f",
                    cluE,traP,eOverp, dEdx,nSigma));
    
    //printf("Cluster E %2.2f, Track P %2.2f, E/P %2.2f, dEdx %2.2f, n sigma %2.2f\n",
    //       cluE,traP,eOverp, dEdx,nSigma);
    
    //
    // Plot different track PID distributions without and with different cut combinations
    //
    fhdEdxvsE  ->Fill(cluE,   dEdx, GetEventWeight());
    fhdEdxvsP  ->Fill(traP,   dEdx, GetEventWeight());
    
    fhEOverPvsE->Fill(cluE, eOverp, GetEventWeight());
    fhEOverPvsP->Fill(traP, eOverp, GetEventWeight()); 
 
    fhNSigmavsE->Fill(cluE, nSigma, GetEventWeight());
    fhNSigmavsP->Fill(traP, nSigma, GetEventWeight()); 
    
    if ( eOverp < fEOverPMax && eOverp > fEOverPMin )
    {
      fhdEdxvsECutEOverP  ->Fill(cluE, dEdx  , GetEventWeight());
      fhdEdxvsPCutEOverP  ->Fill(traP, dEdx  , GetEventWeight());
      
      fhNSigmavsECutEOverP->Fill(cluE, nSigma, GetEventWeight());
      fhNSigmavsPCutEOverP->Fill(traP, nSigma, GetEventWeight());
    }
    
    if ( m02 > fM02Min && m02 < fM02Max )
    {
      fhdEdxvsECutM02  ->Fill(cluE, dEdx  , GetEventWeight());
      fhdEdxvsPCutM02  ->Fill(traP, dEdx  , GetEventWeight());
      fhEOverPvsECutM02->Fill(cluE, eOverp, GetEventWeight());
      fhEOverPvsPCutM02->Fill(traP, eOverp, GetEventWeight());
      fhNSigmavsECutM02->Fill(cluE, nSigma, GetEventWeight());
      fhNSigmavsPCutM02->Fill(traP, nSigma, GetEventWeight());
      
      if ( m20 > fM20Min && m20 < fM20Max )
      {
        fhdEdxvsECutM02AndM20  ->Fill(cluE, dEdx  , GetEventWeight());
        fhdEdxvsPCutM02AndM20  ->Fill(traP, dEdx  , GetEventWeight());
        fhEOverPvsECutM02AndM20->Fill(cluE, eOverp, GetEventWeight());
        fhEOverPvsPCutM02AndM20->Fill(traP, eOverp, GetEventWeight());
        fhNSigmavsECutM02AndM20->Fill(cluE, nSigma, GetEventWeight());
        fhNSigmavsPCutM02AndM20->Fill(traP, nSigma, GetEventWeight());
      }
    } // shower shape cut
    
    if ( dEdx < fdEdxMax && dEdx > fdEdxMin )
    {
      fhEOverPvsECutdEdx->Fill(cluE, eOverp, GetEventWeight());
      fhEOverPvsPCutdEdx->Fill(traP, eOverp, GetEventWeight());
      fhNSigmavsECutdEdx->Fill(cluE, nSigma, GetEventWeight());
      fhNSigmavsPCutdEdx->Fill(traP, nSigma, GetEventWeight());
      
      if ( m02 > fM02Min && m02 < fM02Max )
      {
        fhEOverPvsECutM02CutdEdx->Fill(cluE, eOverp, GetEventWeight());
        fhEOverPvsPCutM02CutdEdx->Fill(traP, eOverp, GetEventWeight());
        fhNSigmavsECutM02CutdEdx->Fill(cluE, nSigma, GetEventWeight());
        fhNSigmavsPCutM02CutdEdx->Fill(traP, nSigma, GetEventWeight());
        
        if ( m20 > fM20Min && m20 < fM20Max )
        {
          fhEOverPvsECutM02AndM20CutdEdx->Fill(cluE, eOverp, GetEventWeight());
          fhEOverPvsPCutM02AndM20CutdEdx->Fill(traP, eOverp, GetEventWeight());        
          fhNSigmavsECutM02AndM20CutdEdx->Fill(cluE, eOverp, GetEventWeight());
          fhNSigmavsPCutM02AndM20CutdEdx->Fill(traP, eOverp, GetEventWeight());
        }
      } // shower shape cut
    }// dE/dx
    
    if ( nSigma < fNSigmaMax && nSigma > fNSigmaMin )
    {
      fhdEdxvsECutNSigma  ->Fill(cluE, dEdx, GetEventWeight());
      fhdEdxvsPCutNSigma  ->Fill(traP, dEdx, GetEventWeight());
      fhEOverPvsECutNSigma->Fill(cluE, eOverp, GetEventWeight());
      fhEOverPvsPCutNSigma->Fill(traP, eOverp, GetEventWeight());     
      
      if ( m02 > fM02Min && m02 < fM02Max )
      {
        fhdEdxvsECutM02CutNSigma  ->Fill(cluE, dEdx, GetEventWeight());
        fhdEdxvsPCutM02CutNSigma  ->Fill(traP, dEdx, GetEventWeight());
        fhEOverPvsECutM02CutNSigma->Fill(cluE, eOverp, GetEventWeight());
        fhEOverPvsPCutM02CutNSigma->Fill(traP, eOverp, GetEventWeight());   
        
        if ( m20 > fM20Min && m20 < fM20Max )
        {
          fhdEdxvsECutM02AndM20CutNSigma  ->Fill(cluE, dEdx, GetEventWeight());
          fhdEdxvsPCutM02AndM20CutNSigma  ->Fill(traP, dEdx, GetEventWeight());
          fhEOverPvsECutM02AndM20CutNSigma->Fill(cluE, eOverp, GetEventWeight());
          fhEOverPvsPCutM02AndM20CutNSigma->Fill(traP, eOverp, GetEventWeight());   
        }
      } // shower shape cut
    }// n sigma
    
    //--------------------------------------------------------------------------
    // Select as hadron or electron or none
    //--------------------------------------------------------------------------
    Int_t pid      = -1;
    Int_t pidNoSS  = -1;
    Int_t pidIndex = -1;

    //
    // Electron?
    //
    
    if ( dEdx   < fdEdxMax   && dEdx   > fdEdxMin   &&
         eOverp < fEOverPMax && eOverp > fEOverPMin && 
         nSigma < fNSigmaMax && nSigma > fNSigmaMin   )
    {
      pidNoSS = AliCaloPID::kElectron;
    }
    
    if ( dEdx   < fdEdxMax   && dEdx   > fdEdxMin   &&
         eOverp < fEOverPMax && eOverp > fEOverPMin && 
         nSigma < fNSigmaMax && nSigma > fNSigmaMin &&
         m02    < fM02Max    && m02    > fM02Min    &&
         m20    < fM20Max    && m20    > fM20Min        )
    {
      pid      = AliCaloPID::kElectron;
      pidIndex = 0;
    } 
  
    //
    // Hadron?
    //
    else if ( dEdx   < fdEdxMaxHad   && dEdx   > fdEdxMinHad    && 
              eOverp < fEOverPMaxHad && eOverp > fEOverPMinHad  &&
              nSigma < fNSigmaMaxHad && nSigma > fNSigmaMinHad    )
    {
        pid      = AliCaloPID::kChargedHadron;
        pidNoSS  = AliCaloPID::kChargedHadron;
        pidIndex = 1;
    }
    
    //--------------------------------------------------------------------------
    // Histograms with track matching residuals, with and without PID selection
    //--------------------------------------------------------------------------
    
    Float_t dEta = calo->GetTrackDz();
    Float_t dPhi = calo->GetTrackDx();
    Bool_t positive = kFALSE;
    if(track) positive = (track->Charge()>0);
    
    AliDebug(1,Form("\tt charge %d, dEta %2.2f, dPhi %2.2f",
                    positive,dEta,dPhi));
    
    //printf("\t charge %d, dEta %2.2f, dPhi %2.2f\n",positive,dEta,dPhi);
    
    // Any accepted cluster
    fhDEtavsE[2][positive]->Fill(cluE, dEta, GetEventWeight());
    fhDPhivsE[2][positive]->Fill(cluE, dPhi, GetEventWeight());
    fhDEtavsP[2][positive]->Fill(traP, dEta, GetEventWeight());
    fhDPhivsP[2][positive]->Fill(traP, dPhi, GetEventWeight());
    
    if ( cluE > 1 && track->P() > 1 ) 
      fhDEtaDPhi[2][positive]->Fill(dEta,dPhi, GetEventWeight());
    
    // Identified clusters
    if ( pidIndex != -1 )
    {
      fhDEtavsE[pidIndex][positive]->Fill(cluE, dEta, GetEventWeight());
      fhDPhivsE[pidIndex][positive]->Fill(cluE, dPhi, GetEventWeight());
      fhDEtavsP[pidIndex][positive]->Fill(traP, dEta, GetEventWeight());
      fhDPhivsP[pidIndex][positive]->Fill(traP, dPhi, GetEventWeight());
      
      if ( cluE > 1 && track->P() > 1 ) 
        fhDEtaDPhi[pidIndex][positive]->Fill(dEta,dPhi, GetEventWeight());
    }
    
    //--------------------------------------------------------------------------
    // Play with the MC stack if available
    //--------------------------------------------------------------------------
    
    // Check origin of the candidates
    Int_t tag = -1 ;
    if ( IsDataMC() )
    {
      tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(), 
                                              calo->GetClusterMCEdepFraction(),
                                              calo->GetNLabels(), GetMC(),
                                              GetReader()->GetNameOfMCEventHederGeneratorToAccept(),
                                              cluE);
      
      AliDebug(1,Form("Origin of candidate, bit map %d",tag));
         
      Int_t mcFlag = -1;
      if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
      {
        fhMCdEdxvsE  [kmcPhoton]->Fill(cluE, dEdx  , GetEventWeight());
        fhMCdEdxvsP  [kmcPhoton]->Fill(traP, dEdx  , GetEventWeight());
        fhMCNSigmavsE[kmcPhoton]->Fill(cluE, nSigma, GetEventWeight());
        fhMCNSigmavsP[kmcPhoton]->Fill(traP, nSigma, GetEventWeight());
        fhMCEOverPvsE[kmcPhoton]->Fill(cluE, eOverp, GetEventWeight());
        fhMCEOverPvsP[kmcPhoton]->Fill(traP, eOverp, GetEventWeight());

        if(pidIndex!=-1)
        {
          fhMCEOverPvsEAfterCuts[kmcPhoton][pidIndex]->Fill(cluE, eOverp, GetEventWeight());
          fhMCEOverPvsPAfterCuts[kmcPhoton][pidIndex]->Fill(traP, eOverp, GetEventWeight());
        }
        
        if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) )
        {
          mcFlag = kmcConversion;
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) &&
                 !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)        )
        {
          mcFlag = kmcPi0Decay;
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) )
        {
          mcFlag = kmcPi0;
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) )
        {
          mcFlag = kmcEta;
        }
        else if ( fhMCE[pidIndex][kmcOtherDecay] )
        {
          mcFlag = kmcOtherDecay;
        }
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron)  )
      {
        mcFlag = kmcAntiNeutron;
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) )
      {
        mcFlag = kmcAntiProton;
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) )
      {
        mcFlag = kmcElectron;
      }
      else if ( fhMCE[pidIndex][kmcOther] )
      {
        mcFlag = kmcOther;
      }
      
      if ( mcFlag > 0 && fhMCdEdxvsE[mcFlag])
      {
        fhMCdEdxvsE  [mcFlag]->Fill(cluE, dEdx  , GetEventWeight());
        fhMCdEdxvsP  [mcFlag]->Fill(traP, dEdx  , GetEventWeight());
        fhMCNSigmavsE[mcFlag]->Fill(cluE, nSigma, GetEventWeight());
        fhMCNSigmavsP[mcFlag]->Fill(traP, nSigma, GetEventWeight());
        fhMCEOverPvsE[mcFlag]->Fill(cluE, eOverp, GetEventWeight());
        fhMCEOverPvsP[mcFlag]->Fill(traP, eOverp, GetEventWeight());
        
        if(pidIndex!=-1)
        {
          fhMCEOverPvsEAfterCuts[mcFlag][pidIndex]->Fill(cluE, eOverp, GetEventWeight());
          fhMCEOverPvsPAfterCuts[mcFlag][pidIndex]->Fill(traP, eOverp, GetEventWeight());       
        }
      }
      
    }// set MC tag and fill Histograms with MC
    
    //---------------------------------
    // Fill some shower shape histograms
    //---------------------------------

    FillShowerShapeHistograms(calo,tag,pidNoSS);
  
    if ( pidNoSS == AliCaloPID::kElectron )
      WeightHistograms(calo);
    
// Do not rely for the moment on AliCaloPID method, all electron/hadron 
// selection implemented here in the lines above
//
//   //-----------------------------------------
//   // PID Shower Shape selection or bit setting
//   //-----------------------------------------
//   
//    // Data, PID check on
//    if ( IsCaloPIDOn() )
//    {
//      // Get most probable PID, 2 options check bayesian PID weights or redo PID
//      // By default, redo PID
//    
//      if ( GetCaloPID()->GetIdentifiedParticleType(calo)!=AliCaloPID::kPhoton )
//      {
//        if ( fAODParticle == AliCaloPID::kElectron )
//          continue;
//        
//        if ( fAODParticle == 0 )
//          pid = AliCaloPID::kChargedHadron ;
//      }
//      
//      AliDebug(1,Form("PDG of identified particle %d",pid));
//    }
        
    
    //--------------------------------------------------------------------------
    // From no on select either identified electron or hadron clusters
    //--------------------------------------------------------------------------
    if ( pidIndex < 0 )  continue;
    
    AliDebug(1,Form("Selection cuts passed: pT %3.2f, type %d",fMomentum.Pt(),pidIndex));
    
    Float_t maxCellFraction = 0;
    Int_t absID = GetCaloUtils()->GetMaxEnergyCell(cells, calo,maxCellFraction);
    if ( absID >= 0 )
      fhMaxCellDiffClusterE[pidIndex]->Fill(fMomentum.E(), maxCellFraction, GetEventWeight());
    
    fhNCellsE[pidIndex] ->Fill(fMomentum.E(), calo->GetNCells()  , GetEventWeight());
    fhNLME   [pidIndex] ->Fill(fMomentum.E(), nMaxima            , GetEventWeight());
    fhTimeE  [pidIndex] ->Fill(fMomentum.E(), calo->GetTOF()*1.e9, GetEventWeight());
    
    //----------------------------
    // Create AOD for analysis
    //----------------------------

    // Add AOD with electron/hadron object to aod branch
    if ( pid == fAODParticle || fAODParticle == 0 ) 
    {
      AliCaloTrackParticle aodpart = AliCaloTrackParticle(fMomentum);
      
      //...............................................
      //Set the indeces of the original caloclusters (MC, ID), and calorimeter
      Int_t label = calo->GetLabel();
      aodpart.SetLabel(label);
      aodpart.SetCaloLabel (calo->GetID(),-1);
      aodpart.SetTrackLabel(GetReader()->GetTrackID(track),-1); // needed instead of track->GetID() since AOD needs some manipulations

      aodpart.SetDetectorTag(GetCalorimeter());
      //printf("Index %d, Id %d, iaod %d\n",icalo, calo->GetID(),GetOutputAODBranch()->GetEntriesFast());
      
      aodpart.SetM02(calo->GetM02());
      aodpart.SetM20(calo->GetM20());
      aodpart.SetNLM(nMaxima);
      aodpart.SetTime(calo->GetTOF()*1e9);
      aodpart.SetNCells(calo->GetNCells());
      Int_t nSM = GetModuleNumber(calo);
      aodpart.SetSModNumber(nSM);
      
      // MC tag
      aodpart.SetTag(tag);
      
      // PID tag
      aodpart.SetIdentifiedParticleType(pid);
      
      AddAODParticle(aodpart);
    }

  }//loop
  
  AliDebug(1,Form("End fill AODs, with %d entries",GetOutputAODBranch()->GetEntriesFast()));
  
}

//________________________________________________
/// Fill histograms for selected clusters.
//________________________________________________
void  AliAnaElectron::MakeAnalysisFillHistograms()
{
  // Access MC information in stack if requested, check that it exists.

  AliVParticle * primary     = 0x0;   
  
  if ( IsDataMC() && !GetMC() )
  {
    AliFatal("MCEvent not available! STOP");
    return;
  } 
  
  // Get vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  //fhVertex->Fill(v[0], v[1], v[2], GetEventWeight());
  if ( TMath::Abs(v[2]) > GetZvertexCut() ) return ; // done elsewhere for Single Event analysis, but there for mixed event
  
  //----------------------------------
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  AliDebug(1,Form("AOD branch entries %d", naod));
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliCaloTrackParticle* ph =  (AliCaloTrackParticle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ph->GetIdentifiedParticleType();

    Int_t pidIndex = 0;// Electron
    if      ( pdg == AliCaloPID::kElectron      ) pidIndex = 0;
    else if ( pdg == AliCaloPID::kChargedHadron ) pidIndex = 1;
    else                                       continue    ;
          
    if ( ((Int_t) ph->GetDetectorTag()) != GetCalorimeter() ) continue;
    
    AliDebug(1,Form("ID Electron: pt %f, phi %f, eta %f", ph->Pt(),ph->Phi(),ph->Eta())) ;
    
    //................................
    //Fill photon histograms 
    Float_t ptcluster  = ph->Pt();
    Float_t phicluster = ph->Phi();
    Float_t etacluster = ph->Eta();
    Float_t ecluster   = ph->E();
    
    fhE [pidIndex]  ->Fill(ecluster,  GetEventWeight());
    fhPt[pidIndex]  ->Fill(ptcluster, GetEventWeight());
      
    fhPhi[pidIndex] ->Fill(ptcluster, phicluster, GetEventWeight());
    fhEta[pidIndex] ->Fill(ptcluster, etacluster, GetEventWeight());
      
    if      ( ecluster   > 0.5 ) fhEtaPhi  [pidIndex]->Fill(etacluster, phicluster, GetEventWeight());
    else if ( GetMinPt() < 0.5 ) fhEtaPhi05[pidIndex]->Fill(etacluster, phicluster, GetEventWeight());
  
    //.......................................
    //Play with the MC data if available
    if ( IsDataMC() )
    {
      //....................................................................
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      if ( label < 0 )
      {
        AliDebug(1,Form("*** bad label ***:  label %d", label));
        continue;
      }
      
      Int_t nprim = GetMC()->GetNumberOfTracks();
      if ( label >=  nprim )
      {
        AliDebug(1,Form("*** large label ***:  label %d, n tracks %d", label, nprim));
        continue ;
      }
      
      Float_t eprim   = 0;
      //Float_t ptprim  = 0;
      primary = GetMC()->GetTrack(label);
      
      if ( !primary )
      {
        AliWarning(Form("*** no primary ***:  label %d", label));
        continue;
      }
      
      eprim   = primary->E();
      //ptprim  = aodprimary->Pt();
      
      Int_t tag =ph->GetTag();
      
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[pidIndex][kmcPhoton] )
      {
        fhMCE  [pidIndex][kmcPhoton] ->Fill(ecluster , GetEventWeight());
        fhMCPt [pidIndex][kmcPhoton] ->Fill(ptcluster, GetEventWeight());
          
        fhMCPhi[pidIndex][kmcPhoton] ->Fill(ecluster, phicluster, GetEventWeight());
        fhMCEta[pidIndex][kmcPhoton] ->Fill(ecluster, etacluster, GetEventWeight());
        
        fhMC2E    [pidIndex][kmcPhoton] ->Fill(ecluster, eprim         , GetEventWeight());
        fhMCDeltaE[pidIndex][kmcPhoton] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        
        if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[pidIndex][kmcConversion] )
        {
          fhMCE  [pidIndex][kmcConversion] ->Fill(ecluster , GetEventWeight());
          fhMCPt [pidIndex][kmcConversion] ->Fill(ptcluster, GetEventWeight());
            
          fhMCPhi[pidIndex][kmcConversion] ->Fill(ecluster, phicluster, GetEventWeight());
          fhMCEta[pidIndex][kmcConversion] ->Fill(ecluster, etacluster, GetEventWeight());
          
          fhMC2E[pidIndex][kmcConversion]     ->Fill(ecluster, eprim         , GetEventWeight());
          fhMCDeltaE[pidIndex][kmcConversion] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) && 
                 !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[pidIndex][kmcPi0Decay] )
        {
          fhMCE  [pidIndex][kmcPi0Decay] ->Fill(ecluster , GetEventWeight());
          fhMCPt [pidIndex][kmcPi0Decay] ->Fill(ptcluster, GetEventWeight());
            
          fhMCPhi[pidIndex][kmcPi0Decay] ->Fill(ecluster, phicluster, GetEventWeight());
          fhMCEta[pidIndex][kmcPi0Decay] ->Fill(ecluster, etacluster, GetEventWeight());
          
          fhMC2E[pidIndex][kmcPi0Decay]     ->Fill(ecluster, eprim         , GetEventWeight());
          fhMCDeltaE[pidIndex][kmcPi0Decay] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        }
        else if ( ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                    GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[pidIndex][kmcOtherDecay] )
        {
          fhMCE  [pidIndex][kmcOtherDecay] ->Fill(ecluster , GetEventWeight());
          fhMCPt [pidIndex][kmcOtherDecay] ->Fill(ptcluster, GetEventWeight());
            
          fhMCPhi[pidIndex][kmcOtherDecay] ->Fill(ecluster, phicluster, GetEventWeight());
          fhMCEta[pidIndex][kmcOtherDecay] ->Fill(ecluster, etacluster, GetEventWeight());
          
          fhMC2E    [pidIndex][kmcOtherDecay] ->Fill(ecluster, eprim         , GetEventWeight());
          fhMCDeltaE[pidIndex][kmcOtherDecay] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [pidIndex][kmcPi0] )
        {
          fhMCE  [pidIndex][kmcPi0] ->Fill(ecluster , GetEventWeight());
          fhMCPt [pidIndex][kmcPi0] ->Fill(ptcluster, GetEventWeight());
            
          fhMCPhi[pidIndex][kmcPi0] ->Fill(ecluster, phicluster, GetEventWeight());
          fhMCEta[pidIndex][kmcPi0] ->Fill(ecluster, etacluster, GetEventWeight());
          
          fhMC2E[pidIndex][kmcPi0]     ->Fill(ecluster, eprim         , GetEventWeight());
          fhMCDeltaE[pidIndex][kmcPi0] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        }
        else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[pidIndex][kmcEta] )
        {
          fhMCE  [pidIndex][kmcEta] ->Fill(ecluster , GetEventWeight());
          fhMCPt [pidIndex][kmcEta] ->Fill(ptcluster, GetEventWeight());
            
          fhMCPhi[pidIndex][kmcEta] ->Fill(ecluster, phicluster, GetEventWeight());
          fhMCEta[pidIndex][kmcEta] ->Fill(ecluster, etacluster, GetEventWeight());
          
          fhMC2E    [pidIndex][kmcEta] ->Fill(ecluster, eprim         , GetEventWeight());
          fhMCDeltaE[pidIndex][kmcEta] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
        }
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[pidIndex][kmcAntiNeutron] )
      {
        fhMCE  [pidIndex][kmcAntiNeutron] ->Fill(ecluster , GetEventWeight());
        fhMCPt [pidIndex][kmcAntiNeutron] ->Fill(ptcluster, GetEventWeight());
          
        fhMCPhi[pidIndex][kmcAntiNeutron] ->Fill(ecluster, phicluster, GetEventWeight());
        fhMCEta[pidIndex][kmcAntiNeutron] ->Fill(ecluster, etacluster, GetEventWeight());
        
        fhMC2E[pidIndex][kmcAntiNeutron]     ->Fill(ecluster, eprim         , GetEventWeight());
        fhMCDeltaE[pidIndex][kmcAntiNeutron] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[pidIndex][kmcAntiProton] )
      {
        fhMCE  [pidIndex][kmcAntiProton] ->Fill(ecluster , GetEventWeight());
        fhMCPt [pidIndex][kmcAntiProton] ->Fill(ptcluster, GetEventWeight());
          
        fhMCPhi[pidIndex][kmcAntiProton] ->Fill(ecluster, phicluster, GetEventWeight());
        fhMCEta[pidIndex][kmcAntiProton] ->Fill(ecluster, etacluster, GetEventWeight());

        fhMC2E    [pidIndex][kmcAntiProton] ->Fill(ecluster, eprim         , GetEventWeight());
        fhMCDeltaE[pidIndex][kmcAntiProton] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[pidIndex][kmcElectron] )
      {
        fhMCE  [pidIndex][kmcElectron] ->Fill(ecluster , GetEventWeight());
        fhMCPt [pidIndex][kmcElectron] ->Fill(ptcluster, GetEventWeight());
          
        fhMCPhi[pidIndex][kmcElectron] ->Fill(ecluster, phicluster, GetEventWeight());
        fhMCEta[pidIndex][kmcElectron] ->Fill(ecluster, etacluster, GetEventWeight());
        
        fhMC2E[pidIndex][kmcElectron]     ->Fill(ecluster, eprim         , GetEventWeight());
        fhMCDeltaE[pidIndex][kmcElectron] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
      }
      else if ( fhMCE[pidIndex][kmcOther] )
      {
        fhMCE  [pidIndex][kmcOther] ->Fill(ecluster , GetEventWeight());
        fhMCPt [pidIndex][kmcOther] ->Fill(ptcluster, GetEventWeight());
          
        fhMCPhi[pidIndex][kmcOther] ->Fill(ecluster, phicluster, GetEventWeight());
        fhMCEta[pidIndex][kmcOther] ->Fill(ecluster, etacluster, GetEventWeight());
        
        fhMC2E    [pidIndex][kmcOther] ->Fill(ecluster, eprim         , GetEventWeight());
        fhMCDeltaE[pidIndex][kmcOther] ->Fill(ecluster, eprim-ecluster, GetEventWeight());
      }
    } // Histograms with MC
  } // aod loop
}

//____________________________________________________
/// Print some relevant parameters set for the analysis.
//____________________________________________________
void AliAnaElectron::Print(const Option_t * opt) const
{
  if ( !opt )
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  //AliAnaCaloTrackCorrBaseClass::Print(" ");

  printf("Calorimeter = %s\n", GetCalorimeterString().Data()) ;
  printf("Select particle type (0-both): %d\n",fAODParticle);
  
  printf("Basic cuts:");
  printf("\t Dist. to bad channel > %2.1f \n",fMinDist);
  printf("\t %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("\t %d < NLM < %d  \n",fNLMCutMin,fNLMCutMax) ;
  printf("\t N cells > %d \n", fNCellsCut);
  
  printf("Electron cuts:\n");
  printf(" \t %2.2f < dEdx < %2.2f  \n",fdEdxMin,fdEdxMax) ;
  printf(" \t %2.2f <  E/P < %2.2f  \n",fEOverPMin,fEOverPMax) ;
  printf(" \t %2.2f < nSig < %2.2f  \n",fNSigmaMin,fNSigmaMax) ;
  printf(" \t %2.2f <  M02 < %2.2f  \n",fM02Min,fM02Max) ;
  printf(" \t %2.2f <  M20 < %2.2f  \n",fM20Min,fM20Max) ;
  
  printf("Hadron cuts:\n");
  printf(" \t %2.2f < dEdx < %2.2f  \n",fdEdxMinHad,fdEdxMaxHad) ;
  printf(" \t %2.2f <  E/P < %2.2f  \n",fEOverPMinHad,fEOverPMaxHad) ;
  printf(" \t %2.2f < nSig < %2.2f  \n",fNSigmaMinHad,fNSigmaMaxHad) ;  
} 

//______________________________________________________
/// Calculate weights and fill histograms.
//______________________________________________________
void AliAnaElectron::WeightHistograms(AliVCluster *clus)
{
  if ( !fFillWeightHistograms || GetMixedEvent() ) return;
  
  AliVCaloCells* cells = 0;
  if ( GetCalorimeter() == kEMCAL ) cells = GetEMCALCells();
  else                              cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) 
  {
    Int_t id = clus->GetCellsAbsId()[ipos];
    
    // Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,GetCalorimeter(), id);
    
    energy    += amp;
    
    if ( amp> ampMax ) 
      ampMax = amp;
  } // energy loop
  
  if ( energy <= 0 )
  {
    AliWarning(Form("Wrong calculated energy %f",energy));
    return;
  }

  //printf("AliAnaElectron::WeightHistograms() - energy %f, ampmax %f, rat %f, lograt %f\n",energy,ampMax,ampMax/energy,TMath::Log(ampMax/energy));
  fhEMaxCellClusterRatio   ->Fill(energy, ampMax/energy            , GetEventWeight());
  fhEMaxCellClusterLogRatio->Fill(energy, TMath::Log(ampMax/energy), GetEventWeight());
  
  // Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp, GetCalorimeter(), id);

    //printf("energy %f, amp %f, rat %f, lograt %f\n",energy,amp,amp/energy,TMath::Log(amp/energy));
    fhECellClusterRatio   ->Fill(energy, amp/energy            , GetEventWeight());
    fhECellClusterLogRatio->Fill(energy, TMath::Log(amp/energy), GetEventWeight());
  }        
  
  // Recalculate shower shape for different W0
  if ( GetCalorimeter()==kEMCAL )
  {
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(1+iw*0.5); 
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      fhLambda0ForW0[iw]->Fill(energy, clus->GetM02(), GetEventWeight());
//    fhLambda1ForW0[iw]->Fill(energy, clus->GetM20(), GetEventWeight());
      
      //printf("\t w %1.1f, l0 %f, l1 %f,\n",3+iw*0.5,clus->GetM02(),clus->GetM20());
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
  } // EMCAL
}
  

