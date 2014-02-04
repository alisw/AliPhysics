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

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
// Copy of AliAnaPhoton just add electron id.
//
// -- Author: Gustavo Conesa (LPSC-IN2P3-CRNS) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "AliVTrack.h"

// --- Analysis system --- 
#include "AliAnaElectron.h" 
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"


ClassImp(AliAnaElectron)
  
//________________________________
AliAnaElectron::AliAnaElectron() : 
    AliAnaCaloTrackCorrBaseClass(),       fCalorimeter(""), 
    fMinDist(0.),                         fMinDist2(0.),                         fMinDist3(0.), 
    fTimeCutMin(-1),                      fTimeCutMax(999999),         
    fNCellsCut(0),                        fNLMCutMin(-1),                        fNLMCutMax(10),
    fFillSSHistograms(kFALSE),             fFillOnlySimpleSSHisto(1),
    fFillWeightHistograms(kFALSE),        fNOriginHistograms(8), 
    fdEdxMin(0.),                         fdEdxMax (200.), 
    fEOverPMin(0),                        fEOverPMax (2),
    fAODParticle(0),
    // Histograms
    fhdEdxvsE(0),                         fhdEdxvsP(0),                 
    fhEOverPvsE(0),                       fhEOverPvsP(0),
    fhdEdxvsECutM02(0),                   fhdEdxvsPCutM02(0),
    fhEOverPvsECutM02(0),                 fhEOverPvsPCutM02(0),
    fhdEdxvsECutEOverP(0),                fhdEdxvsPCutEOverP(0),
    fhEOverPvsECutM02CutdEdx(0),          fhEOverPvsPCutM02CutdEdx(0),
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
  //default ctor
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
      fhMCEOverPvsE     [i] = 0;
      fhMCEOverPvsP     [i] = 0;
    }
    
    for(Int_t i = 0; i < 6; i++)
    {
      fhMCELambda0       [index][i] = 0;
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
  
  //Weight studies
  for(Int_t i =0; i < 14; i++)
  {
    fhLambda0ForW0[i] = 0;
    //fhLambda1ForW0[i] = 0;
  }
  
  //Initialize parameters
  InitParameters();
  
}

//____________________________________________________________________________
Bool_t  AliAnaElectron::ClusterSelected(AliVCluster* calo, TLorentzVector mom, Int_t nMaxima)
{
  //Select clusters if they pass different cuts
  if(GetDebug() > 2) 
    printf("AliAnaElectron::ClusterSelected() Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //.......................................
  //If too small or big energy, skip it
  if(mom.E() < GetMinEnergy() || mom.E() > GetMaxEnergy() ) return kFALSE ; 
  if(GetDebug() > 2) printf("\t Cluster %d Pass E Cut \n",calo->GetID());
  
  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if(tof < fTimeCutMin || tof > fTimeCutMax) return kFALSE;
  if(GetDebug() > 2)  printf("\t Cluster %d Pass Time Cut \n",calo->GetID());
  
  //.......................................
  if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) return kFALSE;
  if(GetDebug() > 2) printf("\t Cluster %d Pass NCell Cut \n",calo->GetID());
  
  //.......................................
  //Check acceptance selection
  if(IsFiducialCutOn()){
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(! in ) return kFALSE ;
  }
  if(GetDebug() > 2) printf("Fiducial cut passed \n");
  
  //.......................................
  //Skip not matched clusters with tracks
  if(!IsTrackMatched(calo, GetReader()->GetInputEvent())) {
      if(GetDebug() > 2) printf("\t Reject non track-matched clusters\n");
      return kFALSE ;
  }
  else if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  
  //...........................................
  // skip clusters with too many maxima
  if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax) return kFALSE ;
  if(GetDebug() > 2) printf(" \t Cluster %d pass NLM %d of out of range \n",calo->GetID(), nMaxima);
  
  //.......................................
  //Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  if(distBad < fMinDist) {//In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  }
  else if(GetDebug() > 2) printf("\t Bad channel cut passed %4.2f > %2.2f \n",distBad, fMinDist);
  //printf("Cluster %d Pass Bad Dist Cut \n",icalo);

  if(GetDebug() > 0) 
    printf("AliAnaElectron::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(), 
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //All checks passed, cluster selected
  return kTRUE;
    
}

//______________________________________________________________________________________________
void  AliAnaElectron::FillShowerShapeHistograms(AliVCluster* cluster, Int_t mcTag, Int_t pidTag)
{
  
  //Fill cluster Shower Shape histograms
  
  if(!fFillSSHistograms || GetMixedEvent()) return;
  
  Int_t pidIndex = 0;// Electron
  if     (pidTag == AliCaloPID::kElectron)      pidIndex = 0;
  else if(pidTag == AliCaloPID::kChargedHadron) pidIndex = 1;
  else return;

  Float_t energy  = cluster->E();
  Int_t   ncells  = cluster->GetNCells();
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  Float_t disp    = cluster->GetDispersion()*cluster->GetDispersion();
  
  Float_t l0   = 0., l1   = 0.;
  Float_t dispp= 0., dEta = 0., dPhi    = 0.; 
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;    
  
  TLorentzVector mom;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
  {
    cluster->GetMomentum(mom,GetVertex(0)) ;}//Assume that come from vertex in straight line
  else{
    Double_t vertex[]={0,0,0};
    cluster->GetMomentum(mom,vertex) ;
  }
  
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  
  fhLam0E[pidIndex] ->Fill(energy,lambda0);
  fhLam1E[pidIndex] ->Fill(energy,lambda1);
  fhDispE[pidIndex] ->Fill(energy,disp);
  
  if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5)
  {
    fhLam0ETRD[pidIndex]->Fill(energy,lambda0);
    fhLam1ETRD[pidIndex]->Fill(energy,lambda1);
    fhDispETRD[pidIndex]->Fill(energy,disp);
  }
  
  if(!fFillOnlySimpleSSHisto)
  {
    if(energy < 2)
    {
      fhNCellsLam0LowE[pidIndex] ->Fill(ncells,lambda0);
      fhEtaLam0LowE[pidIndex]    ->Fill(eta,   lambda0);
      fhPhiLam0LowE[pidIndex]    ->Fill(phi,   lambda0);
    }
    else 
    {
      fhNCellsLam0HighE[pidIndex]->Fill(ncells,lambda0);
      fhEtaLam0HighE[pidIndex]   ->Fill(eta,   lambda0);
      fhPhiLam0HighE[pidIndex]   ->Fill(phi,   lambda0);
    }
    
    if(fCalorimeter == "EMCAL")
    {
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                                   l0, l1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);
      fhDispEtaE        [pidIndex]-> Fill(energy,dEta);
      fhDispPhiE        [pidIndex]-> Fill(energy,dPhi);
      fhSumEtaE         [pidIndex]-> Fill(energy,sEta);
      fhSumPhiE         [pidIndex]-> Fill(energy,sPhi);
      fhSumEtaPhiE      [pidIndex]-> Fill(energy,sEtaPhi);
      fhDispEtaPhiDiffE [pidIndex]-> Fill(energy,dPhi-dEta);
      if(dEta+dPhi>0)fhSphericityE     [pidIndex]-> Fill(energy,(dPhi-dEta)/(dEta+dPhi));
      
      if      (energy < 2 ) fhDispEtaDispPhiEBin[pidIndex][0]->Fill(dEta,dPhi);
      else if (energy < 4 ) fhDispEtaDispPhiEBin[pidIndex][1]->Fill(dEta,dPhi);
      else if (energy < 6 ) fhDispEtaDispPhiEBin[pidIndex][2]->Fill(dEta,dPhi);
      else if (energy < 10) fhDispEtaDispPhiEBin[pidIndex][3]->Fill(dEta,dPhi);
      else                  fhDispEtaDispPhiEBin[pidIndex][4]->Fill(dEta,dPhi);
      
    }
  }
  
  if(IsDataMC())
  {
    AliVCaloCells* cells = 0;
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    //Fill histograms to check shape of embedded clusters
    Float_t fraction = 0;
    if(GetReader()->IsEmbeddedClusterSelectionOn()){//Only working for EMCAL

      Float_t clusterE = 0; // recalculate in case corrections applied.
      Float_t cellE    = 0;
      for(Int_t icell = 0; icell < cluster->GetNCells(); icell++){
        cellE    = cells->GetCellAmplitude(cluster->GetCellAbsId(icell));
        clusterE+=cellE;  
        fraction+=cellE*cluster->GetCellAmplitudeFraction(icell);
      }
      
      //Fraction of total energy due to the embedded signal
      fraction/=clusterE;
      
      if(GetDebug() > 1 ) printf("AliAnaElectron::FillShowerShapeHistogram() - Energy fraction of embedded signal %2.3f, Energy %2.3f\n",fraction, clusterE);
      
      fhEmbeddedSignalFractionEnergy->Fill(clusterE,fraction);
      
    }  // embedded fraction    
    
    // Check the origin and fill histograms
    Int_t index = -1;

    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta))
    {
      index = kmcssPhoton;
            
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron && 
              !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion)))
    {
      index = kmcssElectron;
       
      if(!GetReader()->IsEmbeddedClusterSelectionOn())
      {
        //Check particle overlaps in cluster
        
        //Compare the primary depositing more energy with the rest, if no photon/electron as comon ancestor (conversions), count as other particle
        Int_t ancPDG = 0, ancStatus = -1;
        TLorentzVector momentum; TVector3 prodVertex;
        Int_t ancLabel = 0;
        Int_t noverlaps = 1;      
        for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ ) {
          ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(cluster->GetLabels()[0],cluster->GetLabels()[ilab], GetReader(),ancPDG,ancStatus,momentum,prodVertex);
          if(ancPDG!=22 && TMath::Abs(ancPDG)!=11) noverlaps++;
        }
        
        if(noverlaps == 1){
          fhMCElectronELambda0NoOverlap  ->Fill(energy, lambda0);
        }
        else if(noverlaps == 2){        
          fhMCElectronELambda0TwoOverlap ->Fill(energy, lambda0);
        }
        else if(noverlaps > 2){          
          fhMCElectronELambda0NOverlap   ->Fill(energy, lambda0);
        }
        else {
          printf("AliAnaElectron::FillShowerShapeHistogram() - n overlaps = %d for ancestor %d!!", noverlaps, ancLabel);
        }
      }//No embedding
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        if     (fraction > 0.9) 
        {
          fhEmbedElectronELambda0FullSignal   ->Fill(energy, lambda0);
        }
        else if(fraction > 0.5)
        {
          fhEmbedElectronELambda0MostlySignal ->Fill(energy, lambda0);
        }
        else if(fraction > 0.1)
        { 
          fhEmbedElectronELambda0MostlyBkg    ->Fill(energy, lambda0);
        }
        else
        {
          fhEmbedElectronELambda0FullBkg      ->Fill(energy, lambda0);
        }
      } // embedded      
    }//electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron) && 
               GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
    {
      index = kmcssConversion;
    }//conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  )
    {
      index = kmcssPi0;
    }//pi0
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)  )
    {
      index = kmcssEta;      
    }//eta    
    else 
    {
      index = kmcssOther;
    }//other particles 
    
    fhMCELambda0[pidIndex][index]    ->Fill(energy, lambda0);
    
    if(fCalorimeter == "EMCAL" && !fFillOnlySimpleSSHisto)
    {
      fhMCEDispEta        [pidIndex][index]-> Fill(energy,dEta);
      fhMCEDispPhi        [pidIndex][index]-> Fill(energy,dPhi);
      fhMCESumEtaPhi      [pidIndex][index]-> Fill(energy,sEtaPhi);
      fhMCEDispEtaPhiDiff [pidIndex][index]-> Fill(energy,dPhi-dEta);
      if(dEta+dPhi>0)fhMCESphericity     [pidIndex][index]-> Fill(energy,(dPhi-dEta)/(dEta+dPhi));
    }
    
  }//MC data
  
}

//_____________________________________________
TObjString *  AliAnaElectron::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaElectron ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," %2.2f < dEdx < %2.2f  \n",fdEdxMin,fdEdxMax) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize," %2.2f <  E/P < %2.2f  \n",fEOverPMin, fEOverPMax) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
  parList+=onePar ;  
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  parList += GetCaloPID()->GetPIDParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList()
  
  return new TObjString(parList) ;
}

//_______________________________________________
TList *  AliAnaElectron::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 
	
  Int_t nptbins     = GetHistogramRanges()->GetHistoPtBins();           Float_t ptmax     = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin     = GetHistogramRanges()->GetHistoPtMin(); 
  Int_t nphibins    = GetHistogramRanges()->GetHistoPhiBins();          Float_t phimax    = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin    = GetHistogramRanges()->GetHistoPhiMin(); 
  Int_t netabins    = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax    = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin    = GetHistogramRanges()->GetHistoEtaMin();	
  Int_t ssbins      = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax     = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin     = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t nbins       = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nmax      = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nmin      = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  Int_t ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();         Float_t dedxmax   = GetHistogramRanges()->GetHistodEdxMax();         Float_t dedxmin   = GetHistogramRanges()->GetHistodEdxMin();
  Int_t nPoverEbins = GetHistogramRanges()->GetHistoPOverEBins();       Float_t pOverEmax = GetHistogramRanges()->GetHistoPOverEMax();       Float_t pOverEmin = GetHistogramRanges()->GetHistoPOverEMin();
  Int_t tbins       = GetHistogramRanges()->GetHistoTimeBins() ;        Float_t tmax      = GetHistogramRanges()->GetHistoTimeMax();         Float_t tmin      = GetHistogramRanges()->GetHistoTimeMin();

  
  // MC labels, titles, for originator particles
  TString ptypess[] = { "#gamma","hadron?","#pi^{0}","#eta","#gamma->e^{#pm}","e^{#pm}"} ;
  TString pnamess[] = { "Photon","Hadron" ,"Pi0"    ,"Eta" ,"Conversion"     ,"Electron"} ;
  TString ptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}", "#pi^{0}","#eta",
    "e^{#pm}","#gamma->e^{#pm}","hadron?","Anti-N","Anti-P"                    } ;
  
  TString pname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Pi0","Eta","Electron",
    "Conversion", "Hadron", "AntiNeutron","AntiProton"                        } ;

  
  fhdEdxvsE  = new TH2F ("hdEdxvsE","matched track <dE/dx> vs cluster E ", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsE->SetXTitle("E (GeV)");
  fhdEdxvsE->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsE);  
  
  fhdEdxvsP  = new TH2F ("hdEdxvsP","matched track <dE/dx> vs track P ", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
  fhdEdxvsP->SetXTitle("P (GeV/c)");
  fhdEdxvsP->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsP);  
  
  fhEOverPvsE  = new TH2F ("hEOverPvsE","matched track E/p vs cluster E ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsE->SetXTitle("E (GeV)");
  fhEOverPvsE->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsE);  
  
  fhEOverPvsP  = new TH2F ("hEOverPvsP","matched track E/p vs track P ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
  fhEOverPvsP->SetXTitle("P (GeV/c)");
  fhEOverPvsP->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsP);  
  
  
  fhdEdxvsECutM02  = new TH2F ("hdEdxvsECutM02","matched track <dE/dx> vs cluster E, mild #lambda_{0}^{2} cut", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutM02->SetXTitle("E (GeV)");
  fhdEdxvsECutM02->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsECutM02);
  
  fhdEdxvsPCutM02  = new TH2F ("hdEdxvsPCutM02","matched track <dE/dx> vs track P, mild #lambda_{0}^{2} cut", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutM02->SetXTitle("P (GeV/c)");
  fhdEdxvsPCutM02->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsPCutM02);
  
  fhEOverPvsECutM02  = new TH2F ("hEOverPvsECutM02","matched track E/p vs cluster E, mild #lambda_{0}^{2} cut", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02->SetXTitle("E (GeV)");
  fhEOverPvsECutM02->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsECutM02);
  
  fhEOverPvsPCutM02  = new TH2F ("hEOverPvsPCutM02","matched track E/p vs track P, mild #lambda_{0}^{2} cut", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02->SetXTitle("P (GeV/c)");
  fhEOverPvsPCutM02->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsPCutM02);

  
  fhdEdxvsECutEOverP  = new TH2F ("hdEdxvsECutEOverP","matched track <dE/dx> vs cluster E, cut on E/p", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsECutEOverP->SetXTitle("E (GeV)");
  fhdEdxvsECutEOverP->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsECutEOverP);
  
  fhdEdxvsPCutEOverP  = new TH2F ("hdEdxvsPCutEOverP","matched track <dE/dx> vs track P, cut on E/p", nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
  fhdEdxvsPCutEOverP->SetXTitle("P (GeV/c)");
  fhdEdxvsPCutEOverP->SetYTitle("<dE/dx>");
  outputContainer->Add(fhdEdxvsPCutEOverP);
  
  fhEOverPvsECutM02CutdEdx  = new TH2F ("hEOverPvsECutM02CutdEdx","matched track E/p vs cluster E, dEdx cut, mild #lambda_{0}^{2} cut", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsECutM02CutdEdx->SetXTitle("E (GeV)");
  fhEOverPvsECutM02CutdEdx->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsECutM02CutdEdx);
  
  fhEOverPvsPCutM02CutdEdx  = new TH2F ("hEOverPvsPCutM02CutdEdx","matched track E/p vs track P, dEdx cut, mild #lambda_{0}^{2} cut ", nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
  fhEOverPvsPCutM02CutdEdx->SetXTitle("P (GeV/c)");
  fhEOverPvsPCutM02CutdEdx->SetYTitle("E/p");
  outputContainer->Add(fhEOverPvsPCutM02CutdEdx);

  if(IsDataMC())
  {
    for(Int_t i = 0; i < fNOriginHistograms; i++)
    {
      fhMCdEdxvsE[i]  = new TH2F(Form("hdEdxvsE_MC%s",pname[i].Data()),
                                     Form("matched track <dE/dx> vs cluster E from %s : E ",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhMCdEdxvsE[i]->SetXTitle("E (GeV)");
      fhMCdEdxvsE[i]->SetYTitle("<dE/dx>");
      outputContainer->Add(fhMCdEdxvsE[i]) ;
      
      fhMCdEdxvsP[i]  = new TH2F(Form("hdEdxvsP_MC%s",pname[i].Data()),
                                 Form("matched track <dE/dx> vs track P from %s : E ",ptype[i].Data()),
                                 nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhMCdEdxvsP[i]->SetXTitle("E (GeV)");
      fhMCdEdxvsP[i]->SetYTitle("<dE/dx>");
      outputContainer->Add(fhMCdEdxvsP[i]) ;

      
      fhMCEOverPvsE[i]  = new TH2F(Form("hEOverPvsE_MC%s",pname[i].Data()),
                                 Form("matched track E/p vs cluster E from %s : E ",ptype[i].Data()),
                                 nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhMCEOverPvsE[i]->SetXTitle("E (GeV)");
      fhMCEOverPvsE[i]->SetYTitle("<dE/dx>");
      outputContainer->Add(fhMCEOverPvsE[i]) ;
      
      fhMCEOverPvsP[i]  = new TH2F(Form("hEOverPvsP_MC%s",pname[i].Data()),
                                 Form("matched track E/pvs track P from %s : E ",ptype[i].Data()),
                                 nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhMCEOverPvsP[i]->SetXTitle("E (GeV)");
      fhMCEOverPvsP[i]->SetYTitle("<dE/dx>");
      outputContainer->Add(fhMCEOverPvsP[i]) ;
      
    }
  }
  
  TString pidParticle[] = {"Electron","ChargedHadron"} ;
  
  if(fFillWeightHistograms)
  {
    
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy, for selected electrons",
                                     nptbins,ptmin,ptmax, 100,0,1.); 
    fhECellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("E_{cell i}/E_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy, for selected electrons",
                                        nptbins,ptmin,ptmax, 100,-10,0); 
    fhECellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy, for selected electrons",
                                        nptbins,ptmin,ptmax, 100,0,1.); 
    fhEMaxCellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("E_{max cell}/E_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy, for selected electrons",
                                           nptbins,ptmin,ptmax, 100,-10,0); 
    fhEMaxCellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 14; iw++)
    {
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for selected electrons",1+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLambda0ForW0[iw]->SetXTitle("E_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#lambda^{2}_{0}");
      outputContainer->Add(fhLambda0ForW0[iw]); 
      
      //        fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for selected electrons",1+0.5*iw),
      //                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      //        fhLambda1ForW0[iw]->SetXTitle("E_{cluster}");
      //        fhLambda1ForW0[iw]->SetYTitle("#lambda^{2}_{1}");
      //        outputContainer->Add(fhLambda1ForW0[iw]); 
      
    }
  }  
  
  for(Int_t pidIndex = 0; pidIndex < 2; pidIndex++)
  {
    //Shower shape
    if(fFillSSHistograms)
    {
      fhLam0E[pidIndex]  = new TH2F (Form("h%sLam0E",pidParticle[pidIndex].Data()),
                                     Form("%s: #lambda_{0}^{2} vs E",pidParticle[pidIndex].Data()), 
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam0E[pidIndex]->SetYTitle("#lambda_{0}^{2}");
      fhLam0E[pidIndex]->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam0E[pidIndex]);  
      
      fhLam1E[pidIndex]  = new TH2F (Form("h%sLam1E",pidParticle[pidIndex].Data()),
                                     Form("%s: #lambda_{1}^{2} vs E",pidParticle[pidIndex].Data()), 
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam1E[pidIndex]->SetYTitle("#lambda_{1}^{2}");
      fhLam1E[pidIndex]->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam1E[pidIndex]);  
      
      fhDispE[pidIndex]  = new TH2F (Form("h%sDispE",pidParticle[pidIndex].Data()),
                                     Form("%s: dispersion^{2} vs E",pidParticle[pidIndex].Data()), 
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhDispE[pidIndex]->SetYTitle("D^{2}");
      fhDispE[pidIndex]->SetXTitle("E (GeV) ");
      outputContainer->Add(fhDispE[pidIndex]);
      
      if(fCalorimeter == "EMCAL")
      {
        fhLam0ETRD[pidIndex]  = new TH2F (Form("h%sLam0ETRD",pidParticle[pidIndex].Data()),
                                          Form("%s: #lambda_{0}^{2} vs E, EMCAL SM covered by TRD",pidParticle[pidIndex].Data()), 
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam0ETRD[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        fhLam0ETRD[pidIndex]->SetXTitle("E (GeV)");
        outputContainer->Add(fhLam0ETRD[pidIndex]);  
        
        fhLam1ETRD[pidIndex]  = new TH2F (Form("h%sLam1ETRD",pidParticle[pidIndex].Data()),
                                          Form("%s: #lambda_{1}^{2} vs E, EMCAL SM covered by TRD",pidParticle[pidIndex].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam1ETRD[pidIndex]->SetYTitle("#lambda_{1}^{2}");
        fhLam1ETRD[pidIndex]->SetXTitle("E (GeV)");
        outputContainer->Add(fhLam1ETRD[pidIndex]);  
        
        fhDispETRD[pidIndex]  = new TH2F (Form("h%sDispETRD",pidParticle[pidIndex].Data()),
                                          Form("%s: dispersion^{2} vs E, EMCAL SM covered by TRD",pidParticle[pidIndex].Data()), 
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhDispETRD[pidIndex]->SetYTitle("Dispersion^{2}");
        fhDispETRD[pidIndex]->SetXTitle("E (GeV) ");
        outputContainer->Add(fhDispETRD[pidIndex]);   
      } 
      
      if(!fFillOnlySimpleSSHisto)
      {
        
        fhNCellsLam0LowE[pidIndex]  = new TH2F (Form("h%sNCellsLam0LowE",pidParticle[pidIndex].Data()),
                                                Form("%s: N_{cells} in cluster vs #lambda_{0}^{2}, E < 2 GeV",pidParticle[pidIndex].Data()),
                                                nbins,nmin, nmax, ssbins,ssmin,ssmax); 
        fhNCellsLam0LowE[pidIndex]->SetXTitle("N_{cells}");
        fhNCellsLam0LowE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhNCellsLam0LowE[pidIndex]);  
        
        fhNCellsLam0HighE[pidIndex]  = new TH2F (Form("h%sNCellsLam0HighE",pidParticle[pidIndex].Data()),
                                                 Form("%s: N_{cells} in cluster vs #lambda_{0}^{2}, E > 2 GeV",pidParticle[pidIndex].Data()), 
                                                 nbins,nmin, nmax, ssbins,ssmin,ssmax); 
        fhNCellsLam0HighE[pidIndex]->SetXTitle("N_{cells}");
        fhNCellsLam0HighE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhNCellsLam0HighE[pidIndex]);  
        
        
        fhEtaLam0LowE[pidIndex]  = new TH2F (Form("h%sEtaLam0LowE",pidParticle[pidIndex].Data()),
                                             Form("%s: #eta vs #lambda_{0}^{2}, E < 2 GeV",pidParticle[pidIndex].Data()), 
                                             netabins,etamin,etamax, ssbins,ssmin,ssmax); 
        fhEtaLam0LowE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        fhEtaLam0LowE[pidIndex]->SetXTitle("#eta");
        outputContainer->Add(fhEtaLam0LowE[pidIndex]);  
        
        fhPhiLam0LowE[pidIndex]  = new TH2F (Form("h%sPhiLam0LowE",pidParticle[pidIndex].Data()),
                                             Form("%s: #phi vs #lambda_{0}^{2}, E < 2 GeV",pidParticle[pidIndex].Data()), 
                                             nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
        fhPhiLam0LowE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        fhPhiLam0LowE[pidIndex]->SetXTitle("#phi");
        outputContainer->Add(fhPhiLam0LowE[pidIndex]);  
        
        fhEtaLam0HighE[pidIndex]  = new TH2F (Form("h%sEtaLam0HighE",pidParticle[pidIndex].Data()),
                                              Form("%s: #eta vs #lambda_{0}^{2}, E > 2 GeV",pidParticle[pidIndex].Data()),
                                              netabins,etamin,etamax, ssbins,ssmin,ssmax); 
        fhEtaLam0HighE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        fhEtaLam0HighE[pidIndex]->SetXTitle("#eta");
        outputContainer->Add(fhEtaLam0HighE[pidIndex]);  
        
        fhPhiLam0HighE[pidIndex]  = new TH2F (Form("h%sPhiLam0HighE",pidParticle[pidIndex].Data()),
                                              Form("%s: #phi vs #lambda_{0}^{2}, E > 2 GeV",pidParticle[pidIndex].Data()), 
                                              nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
        fhPhiLam0HighE[pidIndex]->SetYTitle("#lambda_{0}^{2}");
        fhPhiLam0HighE[pidIndex]->SetXTitle("#phi");
        outputContainer->Add(fhPhiLam0HighE[pidIndex]);  
        
        if(fCalorimeter == "EMCAL")
        {
          fhDispEtaE[pidIndex]  = new TH2F (Form("h%sDispEtaE",pidParticle[pidIndex].Data()),
                                            Form("%s: #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",pidParticle[pidIndex].Data()),  
                                            nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhDispEtaE[pidIndex]->SetXTitle("E (GeV)");
          fhDispEtaE[pidIndex]->SetYTitle("#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhDispEtaE[pidIndex]);     
          
          fhDispPhiE[pidIndex]  = new TH2F (Form("h%sDispPhiE",pidParticle[pidIndex].Data()),
                                            Form("%s: #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",pidParticle[pidIndex].Data()),  
                                            nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhDispPhiE[pidIndex]->SetXTitle("E (GeV)");
          fhDispPhiE[pidIndex]->SetYTitle("#sigma^{2}_{#phi #phi}");
          outputContainer->Add(fhDispPhiE[pidIndex]);  
          
          fhSumEtaE[pidIndex]  = new TH2F (Form("h%sSumEtaE",pidParticle[pidIndex].Data()),
                                           Form("%s: #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i})^{2} / #Sigma w_{i} - <#eta>^{2} vs E",pidParticle[pidIndex].Data()),  
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhSumEtaE[pidIndex]->SetXTitle("E (GeV)");
          fhSumEtaE[pidIndex]->SetYTitle("#delta^{2}_{#eta #eta}");
          outputContainer->Add(fhSumEtaE[pidIndex]);     
          
          fhSumPhiE[pidIndex]  = new TH2F (Form("h%sSumPhiE",pidParticle[pidIndex].Data()),
                                           Form("%s: #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i})^{2}/ #Sigma w_{i} - <#phi>^{2} vs E",pidParticle[pidIndex].Data()),  
                                           nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
          fhSumPhiE[pidIndex]->SetXTitle("E (GeV)");
          fhSumPhiE[pidIndex]->SetYTitle("#delta^{2}_{#phi #phi}");
          outputContainer->Add(fhSumPhiE[pidIndex]);  
          
          fhSumEtaPhiE[pidIndex]  = new TH2F (Form("h%sSumEtaPhiE",pidParticle[pidIndex].Data()),
                                              Form("%s: #delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",pidParticle[pidIndex].Data()),  
                                              nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax); 
          fhSumEtaPhiE[pidIndex]->SetXTitle("E (GeV)");
          fhSumEtaPhiE[pidIndex]->SetYTitle("#delta^{2}_{#eta #phi}");
          outputContainer->Add(fhSumEtaPhiE[pidIndex]);
          
          fhDispEtaPhiDiffE[pidIndex]  = new TH2F (Form("h%sDispEtaPhiDiffE",pidParticle[pidIndex].Data()),
                                                   Form("%s: #sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",pidParticle[pidIndex].Data()), 
                                                   nptbins,ptmin,ptmax,200, -10,10); 
          fhDispEtaPhiDiffE[pidIndex]->SetXTitle("E (GeV)");
          fhDispEtaPhiDiffE[pidIndex]->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhDispEtaPhiDiffE[pidIndex]);    
          
          fhSphericityE[pidIndex]  = new TH2F (Form("h%sSphericityE",pidParticle[pidIndex].Data()),
                                               Form("%s: (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",pidParticle[pidIndex].Data()),  
                                               nptbins,ptmin,ptmax, 200, -1,1); 
          fhSphericityE[pidIndex]->SetXTitle("E (GeV)");
          fhSphericityE[pidIndex]->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
          outputContainer->Add(fhSphericityE[pidIndex]);
          
          Int_t bin[] = {0,2,4,6,10,1000};
          for(Int_t i = 0; i < 5; i++)
          {
            fhDispEtaDispPhiEBin[pidIndex][i] = new TH2F (Form("h%sDispEtaDispPhi_EBin%d",pidParticle[pidIndex].Data(),i),
                                                          Form("%s: #sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pidParticle[pidIndex].Data(),bin[i],bin[i+1]), 
                                                          ssbins,ssmin,ssmax , ssbins,ssmin,ssmax); 
            fhDispEtaDispPhiEBin[pidIndex][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
            fhDispEtaDispPhiEBin[pidIndex][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhDispEtaDispPhiEBin[pidIndex][i]); 
          }
        }
      }
    } // Shower shape
        
    if(IsDataMC())
    {
      if(fFillSSHistograms)
      {
        for(Int_t i = 0; i < 6; i++)
        { 
          fhMCELambda0[pidIndex][i]  = new TH2F(Form("h%sELambda0_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                Form("%s like cluster from %s : E vs #lambda_{0}^{2}",pidParticle[pidIndex].Data(),ptypess[i].Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhMCELambda0[pidIndex][i]->SetYTitle("#lambda_{0}^{2}");
          fhMCELambda0[pidIndex][i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCELambda0[pidIndex][i]) ; 
          
          if(fCalorimeter=="EMCAL" && !fFillOnlySimpleSSHisto)
          {
            fhMCEDispEta[pidIndex][i]  = new TH2F (Form("h%sEDispEtaE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                   Form("cluster from %s : %s like, #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",ptypess[i].Data(),pidParticle[pidIndex].Data()),
                                                   nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
            fhMCEDispEta[pidIndex][i]->SetXTitle("E (GeV)");
            fhMCEDispEta[pidIndex][i]->SetYTitle("#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEta[pidIndex][i]);     
            
            fhMCEDispPhi[pidIndex][i]  = new TH2F (Form("h%sEDispPhiE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                   Form("cluster from %s : %s like, #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",ptypess[i].Data(),pidParticle[pidIndex].Data()),
                                                   nptbins,ptmin,ptmax, ssbins,ssmin,ssmax); 
            fhMCEDispPhi[pidIndex][i]->SetXTitle("E (GeV)");
            fhMCEDispPhi[pidIndex][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhMCEDispPhi[pidIndex][i]);  
            
            fhMCESumEtaPhi[pidIndex][i]  = new TH2F (Form("h%sESumEtaPhiE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                     Form("cluster from %s : %s like, #delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                     nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax); 
            fhMCESumEtaPhi[pidIndex][i]->SetXTitle("E (GeV)");
            fhMCESumEtaPhi[pidIndex][i]->SetYTitle("#delta^{2}_{#eta #phi}");
            outputContainer->Add(fhMCESumEtaPhi[pidIndex][i]);
            
            fhMCEDispEtaPhiDiff[pidIndex][i]  = new TH2F (Form("h%sEDispEtaPhiDiffE_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                          Form("cluster from %s : %s like, #sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                          nptbins,ptmin,ptmax,200,-10,10); 
            fhMCEDispEtaPhiDiff[pidIndex][i]->SetXTitle("E (GeV)");
            fhMCEDispEtaPhiDiff[pidIndex][i]->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEtaPhiDiff[pidIndex][i]);    
            
            fhMCESphericity[pidIndex][i]  = new TH2F (Form("h%sESphericity_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                      Form("cluster from %s : %s like, (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",ptypess[i].Data(),pidParticle[pidIndex].Data()),  
                                                      nptbins,ptmin,ptmax, 200,-1,1); 
            fhMCESphericity[pidIndex][i]->SetXTitle("E (GeV)");
            fhMCESphericity[pidIndex][i]->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
            outputContainer->Add(fhMCESphericity[pidIndex][i]);
          }
          
        }// loop    
      }   
    }
    
    //if(IsCaloPIDOn() && pidIndex > 0) continue;
    
    fhNCellsE[pidIndex]  = new TH2F (Form("h%sNCellsE",pidParticle[pidIndex].Data()),
                                     Form("N cells in %s cluster vs E ",pidParticle[pidIndex].Data()),
                                     nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhNCellsE[pidIndex]->SetXTitle("E (GeV)");
    fhNCellsE[pidIndex]->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhNCellsE[pidIndex]);  
    
    fhNLME[pidIndex]  = new TH2F (Form("h%sNLME",pidParticle[pidIndex].Data()),
                                     Form("NLM in %s cluster vs E ",pidParticle[pidIndex].Data()),
                                     nptbins,ptmin,ptmax, 10,0,10);
    fhNLME[pidIndex]->SetXTitle("E (GeV)");
    fhNLME[pidIndex]->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhNLME[pidIndex]);
    
    fhTimeE[pidIndex] = new TH2F(Form("h%sTimeE",pidParticle[pidIndex].Data()),
                                 Form("Time in %s cluster vs E ",pidParticle[pidIndex].Data())
                                 ,nptbins,ptmin,ptmax, tbins,tmin,tmax);
    fhTimeE[pidIndex]->SetXTitle("E (GeV)");
    fhTimeE[pidIndex]->SetYTitle(" t (ns)");
    outputContainer->Add(fhTimeE[pidIndex]);  
    
    fhMaxCellDiffClusterE[pidIndex]  = new TH2F (Form("h%sMaxCellDiffClusterE",pidParticle[pidIndex].Data()),
                                                 Form("%s: energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",pidParticle[pidIndex].Data()),
                                                 nptbins,ptmin,ptmax, 500,0,1.); 
    fhMaxCellDiffClusterE[pidIndex]->SetXTitle("E_{cluster} (GeV) ");
    fhMaxCellDiffClusterE[pidIndex]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
    outputContainer->Add(fhMaxCellDiffClusterE[pidIndex]);  
    
    fhE[pidIndex]  = new TH1F(Form("h%sE",pidParticle[pidIndex].Data()),
                              Form("Number of %s over calorimeter vs energy",pidParticle[pidIndex].Data()),
                              nptbins,ptmin,ptmax); 
    fhE[pidIndex]->SetYTitle("N");
    fhE[pidIndex]->SetXTitle("E_{#gamma}(GeV)");
    outputContainer->Add(fhE[pidIndex]) ;   
    
    fhPt[pidIndex]  = new TH1F(Form("h%sPtElectron",pidParticle[pidIndex].Data()),
                               Form("Number of %s over calorimeter vs p_{T}",pidParticle[pidIndex].Data()),
                               nptbins,ptmin,ptmax); 
    fhPt[pidIndex]->SetYTitle("N");
    fhPt[pidIndex]->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPt[pidIndex]) ; 
    
    fhPhi[pidIndex]  = new TH2F(Form("h%sPhiElectron",pidParticle[pidIndex].Data()),
                                Form("%s: #phi vs p_{T}",pidParticle[pidIndex].Data()),
                                nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhi[pidIndex]->SetYTitle("#phi (rad)");
    fhPhi[pidIndex]->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhi[pidIndex]) ; 
    
    fhEta[pidIndex]  = new TH2F(Form("h%sEta",pidParticle[pidIndex].Data()),
                                Form("%s: #eta vs p_{T}",pidParticle[pidIndex].Data()),
                                nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEta[pidIndex]->SetYTitle("#eta");
    fhEta[pidIndex]->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEta[pidIndex]) ;
    
    fhEtaPhi[pidIndex]  = new TH2F(Form("h%sEtaPhi",pidParticle[pidIndex].Data()),
                                   Form("%s: #eta vs #phi",pidParticle[pidIndex].Data()),
                                   netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi[pidIndex]->SetYTitle("#phi (rad)");
    fhEtaPhi[pidIndex]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi[pidIndex]) ;
    if(GetMinPt() < 0.5)
    {
      fhEtaPhi05[pidIndex]  = new TH2F(Form("h%sEtaPhi05",pidParticle[pidIndex].Data()),
                                       Form("%s: #eta vs #phi, E > 0.5",pidParticle[pidIndex].Data()),
                                       netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05[pidIndex]->SetYTitle("#phi (rad)");
      fhEtaPhi05[pidIndex]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05[pidIndex]) ;
    }
    
    
    if(IsDataMC())
    {      
      for(Int_t i = 0; i < fNOriginHistograms; i++)
      { 
        fhMCE[pidIndex][i]  = new TH1F(Form("h%sE_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                       Form("%s like cluster from %s : E ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                       nptbins,ptmin,ptmax); 
        fhMCE[pidIndex][i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCE[pidIndex][i]) ; 
        
        fhMCPt[pidIndex][i]  = new TH1F(Form("h%sPt_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                        Form("%s like cluster from %s : p_{T} ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                        nptbins,ptmin,ptmax); 
        fhMCPt[pidIndex][i]->SetXTitle("p_{T} (GeV/c)");
        outputContainer->Add(fhMCPt[pidIndex][i]) ;
        
        fhMCEta[pidIndex][i]  = new TH2F(Form("h%sEta_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like cluster from %s : #eta ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                         nptbins,ptmin,ptmax,netabins,etamin,etamax); 
        fhMCEta[pidIndex][i]->SetYTitle("#eta");
        fhMCEta[pidIndex][i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEta[pidIndex][i]) ;
        
        fhMCPhi[pidIndex][i]  = new TH2F(Form("h%sPhi_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like cluster from %s : #phi ",pidParticle[pidIndex].Data(),ptype[i].Data()),
                                         nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
        fhMCPhi[pidIndex][i]->SetYTitle("#phi (rad)");
        fhMCPhi[pidIndex][i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCPhi[pidIndex][i]) ;
        
        
        fhMCDeltaE[pidIndex][i]  = new TH2F (Form("h%sDeltaE_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                             Form("%s like MC - Reco E from %s",pidParticle[pidIndex].Data(),pname[i].Data()), 
                                             nptbins,ptmin,ptmax, 200,-50,50); 
        fhMCDeltaE[pidIndex][i]->SetXTitle("#Delta E (GeV)");
        outputContainer->Add(fhMCDeltaE[pidIndex][i]);
        
        fhMC2E[pidIndex][i]  = new TH2F (Form("h%s2E_MC%s",pidParticle[pidIndex].Data(),pname[i].Data()),
                                         Form("%s like E distribution, reconstructed vs generated from %s",pidParticle[pidIndex].Data(),pname[i].Data()), 
                                         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMC2E[pidIndex][i]->SetXTitle("E_{rec} (GeV)");
        fhMC2E[pidIndex][i]->SetYTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMC2E[pidIndex][i]);          
        
      }
    } // MC
    
  }// pid Index
  
  
  if(fFillSSHistograms)
  {
    if(IsDataMC())
    {
      if(!GetReader()->IsEmbeddedClusterSelectionOn())
      {
        fhMCElectronELambda0NoOverlap  = new TH2F("hELambda0_MCElectron_NoOverlap",
                                                  "cluster from Electron : E vs #lambda_{0}^{2}",
                                                  nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0NoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCElectronELambda0NoOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCElectronELambda0NoOverlap) ; 
        
        fhMCElectronELambda0TwoOverlap  = new TH2F("hELambda0_MCElectron_TwoOverlap",
                                                   "cluster from Electron : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0TwoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCElectronELambda0TwoOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCElectronELambda0TwoOverlap) ; 
        
        fhMCElectronELambda0NOverlap  = new TH2F("hELambda0_MCElectron_NOverlap",
                                                 "cluster from Electron : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCElectronELambda0NOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCElectronELambda0NOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCElectronELambda0NOverlap) ; 
        
      } //No embedding     
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        
        fhEmbeddedSignalFractionEnergy  = new TH2F("hEmbeddedSignal_FractionEnergy",
                                                   "Energy Fraction of embedded signal versus cluster energy",
                                                   nptbins,ptmin,ptmax,100,0.,1.); 
        fhEmbeddedSignalFractionEnergy->SetYTitle("Fraction");
        fhEmbeddedSignalFractionEnergy->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbeddedSignalFractionEnergy) ; 
        
        fhEmbedElectronELambda0FullSignal  = new TH2F("hELambda0_EmbedElectron_FullSignal",
                                                      "cluster from Electron embedded with more than 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0FullSignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedElectronELambda0FullSignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0FullSignal) ; 
        
        fhEmbedElectronELambda0MostlySignal  = new TH2F("hELambda0_EmbedElectron_MostlySignal",
                                                        "cluster from Electron embedded with 50% to 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0MostlySignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedElectronELambda0MostlySignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0MostlySignal) ; 
        
        fhEmbedElectronELambda0MostlyBkg  = new TH2F("hELambda0_EmbedElectron_MostlyBkg",
                                                     "cluster from Electron embedded with 10% to 50% energy in cluster : E vs #lambda_{0}^{2}",
                                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0MostlyBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedElectronELambda0MostlyBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0MostlyBkg) ; 
        
        fhEmbedElectronELambda0FullBkg  = new TH2F("hELambda0_EmbedElectron_FullBkg",
                                                   "cluster from Electronm embedded with 0% to 10% energy in cluster : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedElectronELambda0FullBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedElectronELambda0FullBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedElectronELambda0FullBkg) ; 
        
        
      }// embedded histograms
      
    }//Histos with MC
    
  }// Fill SS MC histograms
  
  return outputContainer ;
  
}

//_________________________
void AliAnaElectron::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaElectron::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaElectron::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//___________________________________
void AliAnaElectron::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaElectron_");
  
  fCalorimeter = "EMCAL" ;
  fMinDist     = 2.;
  fMinDist2    = 4.;
  fMinDist3    = 5.;
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCut   = 0;
  
  fdEdxMin     = 76.; // for LHC11a, but for LHC11c pass1 56.                
  fdEdxMax     = 85.; // for LHC11a, but for LHC11c pass1 64.   

  fEOverPMin   = 0.8; // for LHC11a, but for LHC11c pass1 0.9                  
  fEOverPMax   = 1.2; // for LHC11a and LHC11c pass1  
  
}

//_________________________________________
void  AliAnaElectron::MakeAnalysisFillAOD() 
{
  //Do photon analysis and fill aods
  
  //Get the vertex 
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  //Select the Calorimeter of the photon
  TObjArray * pl = 0x0; 
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  if(!pl) {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
  
  //Init arrays, variables, get number of clusters
  TLorentzVector mom, mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  //List to be used in conversion analysis, to tag the cluster as candidate for conversion
  
  if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  
  //----------------------------------------------------
  // Fill AOD with PHOS/EMCAL AliAODPWG4Particle objects
  //----------------------------------------------------
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++){    
	  
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));	
    //printf("calo %d, %f\n",icalo,calo->E());
    
    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      //Get the vertex and check it is not too large in z
      if(TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut()) continue;
    }
    
    //Cluster selection, not charged, with photon id and in fiducial cut	  
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;}//Assume that come from vertex in straight line
    else{
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
    
    //--------------------------------------
    // Cluster selection
    //--------------------------------------
    AliVCaloCells* cells    = 0;
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells); // NLM
    if(!ClusterSelected(calo,mom,nMaxima)) continue;
    
    //----------------------------
    //Create AOD for analysis
    //----------------------------
    AliAODPWG4Particle aodpart = AliAODPWG4Particle(mom);
    
    //...............................................
    //Set the indeces of the original caloclusters (MC, ID), and calorimeter  
    Int_t label = calo->GetLabel();
    aodpart.SetLabel(label);
    aodpart.SetCaloLabel(calo->GetID(),-1);
    aodpart.SetDetector(fCalorimeter);
    //printf("Index %d, Id %d, iaod %d\n",icalo, calo->GetID(),GetOutputAODBranch()->GetEntriesFast());
    
    //...............................................
    //Set bad channel distance bit
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if     (distBad > fMinDist3) aodpart.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodpart.SetDistToBad(1) ; 
    else                         aodpart.SetDistToBad(0) ;
    //printf("DistBad %f Bit %d\n",distBad, aodpart.DistToBad());
    
    //-------------------------------------
    //PID selection via dEdx
    //-------------------------------------
    
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(calo, GetReader()->GetInputEvent());

    if(!track)
    {
      printf("AliAnaElectron::MakeAnalysisFillAOD() - Null track");
      continue;
    }
    
    //printf("track dedx %f, p %f, cluster E %f\n",track->GetTPCsignal(),track->P(),calo->E());
    Float_t dEdx = track->GetTPCsignal();
    Float_t eOverp = calo->E()/track->P();
    
    fhdEdxvsE->Fill(calo->E(), dEdx);
    fhdEdxvsP->Fill(track->P(),dEdx);
    
    if( eOverp < fEOverPMax && eOverp > fEOverPMin)
    {
      fhdEdxvsECutEOverP  ->Fill(calo->E(), dEdx);
      fhdEdxvsPCutEOverP  ->Fill(track->P(),dEdx);
    }
    
    //Apply a mild cut on the cluster SS and check the value of dEdX and EOverP
    Float_t m02 = calo->GetM02();
    if(m02 > 0.1 && m02 < 0.4)
    {
      fhdEdxvsECutM02  ->Fill(calo->E(), dEdx);
      fhdEdxvsPCutM02  ->Fill(track->P(),dEdx);
      fhEOverPvsECutM02->Fill(calo->E(),  eOverp);
      fhEOverPvsPCutM02->Fill(track->P(), eOverp);
    }
    
    Int_t pid  = AliCaloPID::kChargedHadron;
    
    if( dEdx < fdEdxMax && dEdx > fdEdxMin)
    {
      fhEOverPvsE->Fill(calo->E(),  eOverp);
      fhEOverPvsP->Fill(track->P(), eOverp);
      
      if(m02 > 0.1 && m02 < 0.4)
      {
        fhEOverPvsECutM02CutdEdx->Fill(calo->E(),  eOverp);
        fhEOverPvsPCutM02CutdEdx->Fill(track->P(), eOverp);
      }
      
      if( eOverp < fEOverPMax && eOverp > fEOverPMin)
      {
        pid  = AliCaloPID::kElectron;
      } //E/p
      
    }// dEdx
    
    aodpart.SetIdentifiedParticleType(pid);    
    
    Int_t pidIndex = 0;// Electron
    if(pid == AliCaloPID::kChargedHadron) pidIndex = 1;
  
    //--------------------------------------------------------------------------------------
    //Play with the MC stack if available
    //--------------------------------------------------------------------------------------
    
    //Check origin of the candidates
    if(IsDataMC())
    {
      Int_t tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader());
      aodpart.SetTag(tag);
      
      if(GetDebug() > 0)
        printf("AliAnaElectron::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodpart.GetTag());
         
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[pidIndex][kmcPhoton])
      {
        fhMCdEdxvsE  [kmcPhoton]->Fill(calo ->E(), dEdx);
        fhMCdEdxvsP  [kmcPhoton]->Fill(track->P(), dEdx);
        fhMCEOverPvsE[kmcPhoton]->Fill(calo ->E(), eOverp);
        fhMCEOverPvsP[kmcPhoton]->Fill(track->P(), eOverp);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[pidIndex][kmcConversion])
        {
          fhMCdEdxvsE  [kmcConversion]->Fill(calo ->E(), dEdx);
          fhMCdEdxvsP  [kmcConversion]->Fill(track->P(), dEdx);
          fhMCEOverPvsE[kmcConversion]->Fill(calo ->E(), eOverp);
          fhMCEOverPvsP[kmcConversion]->Fill(track->P(), eOverp);
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) &&
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[pidIndex][kmcPi0Decay])
        {
          fhMCdEdxvsE  [kmcPi0Decay]->Fill(calo ->E(), dEdx);
          fhMCdEdxvsP  [kmcPi0Decay]->Fill(track->P(), dEdx);
          fhMCEOverPvsE[kmcPi0Decay]->Fill(calo ->E(), eOverp);
          fhMCEOverPvsP[kmcPi0Decay]->Fill(track->P(), eOverp);
        }
        else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) ||
                  GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[pidIndex][kmcOtherDecay])
        {
          fhMCdEdxvsE  [kmcOtherDecay]->Fill(calo ->E(), dEdx);
          fhMCdEdxvsP  [kmcOtherDecay]->Fill(track->P(), dEdx);
          fhMCEOverPvsE[kmcOtherDecay]->Fill(calo ->E(), eOverp);
          fhMCEOverPvsP[kmcOtherDecay]->Fill(track->P(), eOverp);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [pidIndex][kmcPi0])
        {
          fhMCdEdxvsE  [kmcPi0]->Fill(calo ->E(), dEdx);
          fhMCdEdxvsP  [kmcPi0]->Fill(track->P(), dEdx);
          fhMCEOverPvsE[kmcPi0]->Fill(calo ->E(), eOverp);
          fhMCEOverPvsP[kmcPi0]->Fill(track->P(), eOverp);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[pidIndex][kmcEta])
        {
          fhMCdEdxvsE  [kmcEta]->Fill(calo ->E(), dEdx);
          fhMCdEdxvsP  [kmcEta]->Fill(track->P(), dEdx);
          fhMCEOverPvsE[kmcEta]->Fill(calo ->E(), eOverp);
          fhMCEOverPvsP[kmcEta]->Fill(track->P(), eOverp);
        }
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[pidIndex][kmcAntiNeutron])
      {
        fhMCdEdxvsE  [kmcAntiNeutron]->Fill(calo ->E(), dEdx);
        fhMCdEdxvsP  [kmcAntiNeutron]->Fill(track->P(), dEdx);
        fhMCEOverPvsE[kmcAntiNeutron]->Fill(calo ->E(), eOverp);
        fhMCEOverPvsP[kmcAntiNeutron]->Fill(track->P(), eOverp);
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[pidIndex][kmcAntiProton])
      {
        fhMCdEdxvsE  [kmcAntiProton]->Fill(calo ->E(), dEdx);
        fhMCdEdxvsP  [kmcAntiProton]->Fill(track->P(), dEdx);
        fhMCEOverPvsE[kmcAntiProton]->Fill(calo ->E(), eOverp);
        fhMCEOverPvsP[kmcAntiProton]->Fill(track->P(), eOverp);
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[pidIndex][kmcElectron])
      {
        fhMCdEdxvsE  [kmcElectron]->Fill(calo ->E(), dEdx);
        fhMCdEdxvsP  [kmcElectron]->Fill(track->P(), dEdx);
        fhMCEOverPvsE[kmcElectron]->Fill(calo ->E(), eOverp);
        fhMCEOverPvsP[kmcElectron]->Fill(track->P(), eOverp);
      }
      else if( fhMCE[pidIndex][kmcOther])
      {
        fhMCdEdxvsE  [kmcOther]->Fill(calo ->E(), dEdx);
        fhMCdEdxvsP  [kmcOther]->Fill(track->P(), dEdx);
        fhMCEOverPvsE[kmcOther]->Fill(calo ->E(), eOverp);
        fhMCEOverPvsP[kmcOther]->Fill(track->P(), eOverp);
      }
    }// set MC tag and fill Histograms with MC
    
    //---------------------------------
    //Fill some shower shape histograms
    //---------------------------------

    FillShowerShapeHistograms(calo,aodpart.GetTag(),pid);
  
    if(pid == AliCaloPID::kElectron)
      WeightHistograms(calo);
    
    //-----------------------------------------
    //PID Shower Shape selection or bit setting
    //-----------------------------------------
    
    // Data, PID check on
    if(IsCaloPIDOn())
    {
      // Get most probable PID, 2 options check bayesian PID weights or redo PID
      // By default, redo PID
    
      if(GetCaloPID()->GetIdentifiedParticleType(calo)!=AliCaloPID::kPhoton)
      {
        continue;
      }
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodpart.GetIdentifiedParticleType());
      
    }
        
    if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",
                              aodpart.Pt(), aodpart.GetIdentifiedParticleType());
    
    //FIXME, this to MakeAnalysisFillHistograms ...
    Float_t maxCellFraction = 0;
    
    Int_t absID = GetCaloUtils()->GetMaxEnergyCell(cells, calo,maxCellFraction);
    if(absID>=0)fhMaxCellDiffClusterE[pidIndex]->Fill(aodpart.E(),maxCellFraction);
    fhNCellsE[pidIndex]            ->Fill(aodpart.E(),calo->GetNCells());
    fhNLME[pidIndex]               ->Fill(aodpart.E(),nMaxima);
    fhTimeE[pidIndex]              ->Fill(aodpart.E(),calo->GetTOF()*1.e9);

    //Add AOD with electron/hadron object to aod branch
    if ( pid == fAODParticle || fAODParticle == 0 ) 
    {
      AddAODParticle(aodpart);
    }

  }//loop
  
  if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//________________________________________________
void  AliAnaElectron::MakeAnalysisFillHistograms() 
{
  //Fill histograms
  
  //-------------------------------------------------------------------
  // Access MC information in stack if requested, check that it exists.	
  AliStack         * stack       = 0x0;
  TParticle        * primary     = 0x0;   
  TClonesArray     * mcparticles = 0x0;
  AliAODMCParticle * aodprimary  = 0x0; 
  
  if(IsDataMC()){
    
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;
      if(!stack) {
        printf("AliAnaElectron::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
        abort();
      }
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      //Get the list of MC particles
      mcparticles = GetReader()->GetAODMCParticles();
      if(!mcparticles && GetDebug() > 0) 	{
        printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
      }	
    }
  }// is data and MC
  
  
  // Get vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  //fhVertex->Fill(v[0],v[1],v[2]);  
  if(TMath::Abs(v[2]) > GetZvertexCut()) return ; // done elsewhere for Single Event analysis, but there for mixed event
  
  //----------------------------------
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ph->GetIdentifiedParticleType();

    Int_t pidIndex = 0;// Electron
    if     (pdg == AliCaloPID::kElectron)      pidIndex = 0;
    else if(pdg == AliCaloPID::kChargedHadron) pidIndex = 1;
    else                                       continue    ;
          
    if(ph->GetDetector() != fCalorimeter) continue;
    
    if(GetDebug() > 2) 
      printf("AliAnaElectron::MakeAnalysisFillHistograms() - ID Electron: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
    
    //................................
    //Fill photon histograms 
    Float_t ptcluster  = ph->Pt();
    Float_t phicluster = ph->Phi();
    Float_t etacluster = ph->Eta();
    Float_t ecluster   = ph->E();
    
    fhE[pidIndex]   ->Fill(ecluster);
    fhPt[pidIndex]  ->Fill(ptcluster);
    fhPhi[pidIndex] ->Fill(ptcluster,phicluster);
    fhEta[pidIndex] ->Fill(ptcluster,etacluster);    
    if     (ecluster > 0.5)   fhEtaPhi[pidIndex]  ->Fill(etacluster, phicluster);
    else if(GetMinPt() < 0.5) fhEtaPhi05[pidIndex]->Fill(etacluster, phicluster);
  
    //.......................................
    //Play with the MC data if available
    if(IsDataMC()){
      
      //....................................................................
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      if(label < 0) {
        if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
        continue;
      }
      
      Float_t eprim   = 0;
      //Float_t ptprim  = 0;
      if(GetReader()->ReadStack()){
        
        if(label >=  stack->GetNtrack()) {
          if(GetDebug() > 2)  printf("AliAnaElectron::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
          continue ;
        }
        
        primary = stack->Particle(label);
        if(!primary){
          printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        
        eprim   = primary->Energy();
        //ptprim  = primary->Pt();
        
      }
      else if(GetReader()->ReadAODMCParticles()){
        //Check which is the input
        if(ph->GetInputFileIndex() == 0){
          if(!mcparticles) continue;
          if(label >=  mcparticles->GetEntriesFast()) {
            if(GetDebug() > 2)  printf("AliAnaElectron::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                       label, mcparticles->GetEntriesFast());
            continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles->At(label);
          
        }
        
        if(!aodprimary){
          printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        
        eprim   = aodprimary->E();
        //ptprim  = aodprimary->Pt();
        
      }
      
      Int_t tag =ph->GetTag();
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[pidIndex][kmcPhoton])
      {
        fhMCE  [pidIndex][kmcPhoton] ->Fill(ecluster);
        fhMCPt [pidIndex][kmcPhoton] ->Fill(ptcluster);
        fhMCPhi[pidIndex][kmcPhoton] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][kmcPhoton] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][kmcPhoton]     ->Fill(ecluster, eprim); 
        fhMCDeltaE[pidIndex][kmcPhoton] ->Fill(ecluster,eprim-ecluster);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[pidIndex][kmcConversion])
        {
          fhMCE  [pidIndex][kmcConversion] ->Fill(ecluster);
          fhMCPt [pidIndex][kmcConversion] ->Fill(ptcluster);
          fhMCPhi[pidIndex][kmcConversion] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][kmcConversion] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][kmcConversion]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][kmcConversion] ->Fill(ecluster,eprim-ecluster);
          
        }			        
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[pidIndex][kmcPi0Decay])
        {
          fhMCE  [pidIndex][kmcPi0Decay] ->Fill(ecluster);
          fhMCPt [pidIndex][kmcPi0Decay] ->Fill(ptcluster);
          fhMCPhi[pidIndex][kmcPi0Decay] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][kmcPi0Decay] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][kmcPi0Decay]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][kmcPi0Decay] ->Fill(ecluster,eprim-ecluster);
        }
        else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                  GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[pidIndex][kmcOtherDecay])
        {
          fhMCE  [pidIndex][kmcOtherDecay] ->Fill(ecluster);
          fhMCPt [pidIndex][kmcOtherDecay] ->Fill(ptcluster);
          fhMCPhi[pidIndex][kmcOtherDecay] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][kmcOtherDecay] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][kmcOtherDecay]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][kmcOtherDecay] ->Fill(ecluster,eprim-ecluster);          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [pidIndex][kmcPi0])
        {
          fhMCE  [pidIndex][kmcPi0] ->Fill(ecluster);
          fhMCPt [pidIndex][kmcPi0] ->Fill(ptcluster);
          fhMCPhi[pidIndex][kmcPi0] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][kmcPi0] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][kmcPi0]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][kmcPi0] ->Fill(ecluster,eprim-ecluster);
          
        } 
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[pidIndex][kmcEta])
        {
          fhMCE  [pidIndex][kmcEta] ->Fill(ecluster);
          fhMCPt [pidIndex][kmcEta] ->Fill(ptcluster);
          fhMCPhi[pidIndex][kmcEta] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][kmcEta] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][kmcEta]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][kmcEta] ->Fill(ecluster,eprim-ecluster);
          
        }      
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[pidIndex][kmcAntiNeutron])
      {
        fhMCE  [pidIndex][kmcAntiNeutron] ->Fill(ecluster);
        fhMCPt [pidIndex][kmcAntiNeutron] ->Fill(ptcluster);
        fhMCPhi[pidIndex][kmcAntiNeutron] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][kmcAntiNeutron] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][kmcAntiNeutron]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][kmcAntiNeutron] ->Fill(ecluster,eprim-ecluster);
        
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[pidIndex][kmcAntiProton])
      {
        fhMCE  [pidIndex][kmcAntiProton] ->Fill(ecluster);
        fhMCPt [pidIndex][kmcAntiProton] ->Fill(ptcluster);
        fhMCPhi[pidIndex][kmcAntiProton] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][kmcAntiProton] ->Fill(ecluster,etacluster);

        fhMC2E[pidIndex][kmcAntiProton]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][kmcAntiProton] ->Fill(ecluster,eprim-ecluster);
        
      } 
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[pidIndex][kmcElectron])
      {
        fhMCE  [pidIndex][kmcElectron] ->Fill(ecluster);
        fhMCPt [pidIndex][kmcElectron] ->Fill(ptcluster);
        fhMCPhi[pidIndex][kmcElectron] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][kmcElectron] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][kmcElectron]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][kmcElectron] ->Fill(ecluster,eprim-ecluster);
        
      }     
      else if( fhMCE[pidIndex][kmcOther]){
        fhMCE  [pidIndex][kmcOther] ->Fill(ecluster);
        fhMCPt [pidIndex][kmcOther] ->Fill(ptcluster);
        fhMCPhi[pidIndex][kmcOther] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][kmcOther] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][kmcOther]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][kmcOther] ->Fill(ecluster,eprim-ecluster);
        
      }
      
    }//Histograms with MC
    
  }// aod loop
  
}

//____________________________________________________
void AliAnaElectron::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");

  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf(" %2.2f < dEdx < %2.2f  \n",fdEdxMin,fdEdxMax) ;
  printf(" %2.2f <  E/P < %2.2f  \n",fEOverPMin,fEOverPMax) ;
  printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
  printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
  printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
} 

//______________________________________________________
void AliAnaElectron::WeightHistograms(AliVCluster *clus)
{
  // Calculate weights and fill histograms
  
  if(!fFillWeightHistograms || GetMixedEvent()) return;
  
  AliVCaloCells* cells = 0;
  if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
  else                        cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    energy    += amp;
    
    if(amp> ampMax) 
      ampMax = amp;
    
  } // energy loop       
  
  if(energy <=0 ) {
    printf("AliAnaCalorimeterQA::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  //printf("AliAnaElectron::WeightHistograms() - energy %f, ampmax %f, rat %f, lograt %f\n",energy,ampMax,ampMax/energy,TMath::Log(ampMax/energy));
  fhEMaxCellClusterRatio   ->Fill(energy,ampMax/energy);
  fhEMaxCellClusterLogRatio->Fill(energy,TMath::Log(ampMax/energy));
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp, fCalorimeter, id);

    //printf("energy %f, amp %f, rat %f, lograt %f\n",energy,amp,amp/energy,TMath::Log(amp/energy));
    fhECellClusterRatio   ->Fill(energy,amp/energy);
    fhECellClusterLogRatio->Fill(energy,TMath::Log(amp/energy));
  }        
  
  //Recalculate shower shape for different W0
  if(fCalorimeter=="EMCAL"){
    
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    
    for(Int_t iw = 0; iw < 14; iw++){
      
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(1+iw*0.5); 
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      fhLambda0ForW0[iw]->Fill(energy,clus->GetM02());
      //fhLambda1ForW0[iw]->Fill(energy,clus->GetM20());
      
      //printf("\t w %1.1f, l0 %f, l1 %f,\n",3+iw*0.5,clus->GetM02(),clus->GetM20());
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
  }// EMCAL
}
  

