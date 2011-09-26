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
/* $Id: AliAnaElectron.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

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
  
//____________________________________________________________________________
AliAnaElectron::AliAnaElectron() : 
    AliAnaPartCorrBaseClass(),    fCalorimeter(""), 
    fMinDist(0.),                 fMinDist2(0.),                fMinDist3(0.), 
    fTimeCutMin(-1),              fTimeCutMax(999999),         
    fNCellsCut(0),                fFillSSHistograms(kFALSE),    fNOriginHistograms(8), 
    fdEdxMin(0.),                 fdEdxMax (200.), 
    fEOverPMin(0),                fEOverPMax (2),
    // Histograms
    fhdEdxvsE(0),                 fhdEdxvsP(0),                 
    fhEOverPvsE(0),               fhEOverPvsP(0),
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
  for(Int_t index = 0; index < 2; index++){
    
    fhNCellsE [index] = 0;                  
    fhMaxCellDiffClusterE[index] = 0;
    fhE       [index] = 0;                        
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
    
    for(Int_t i = 0; i < 10; i++){
      fhMCPt     [index][i] = 0;
      fhMCE      [index][i] = 0;
      fhMCPhi    [index][i] = 0;
      fhMCEta    [index][i] = 0;
      fhMCDeltaE [index][i] = 0;                
      fhMC2E     [index][i] = 0;              
    }
    
    for(Int_t i = 0; i < 6; i++){
      fhMCELambda0[index][i]= 0;
    }
    
  }
  //Initialize parameters
  InitParameters();
  
}

//__________________________________________________________________
Bool_t  AliAnaElectron::ClusterSelected(AliVCluster* calo, TLorentzVector mom) 
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
  if(!IsTrackMatched(calo)) {
      if(GetDebug() > 2) printf("\t Reject non track-matched clusters\n");
      return kFALSE ;
  }
  else if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  
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

//__________________________________________________________________
void  AliAnaElectron::FillShowerShapeHistograms(AliVCluster* cluster, const Int_t mcTag, const Int_t pidTag){
  
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
  
  TLorentzVector mom;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
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
  
  if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5){
    fhLam0ETRD[pidIndex]->Fill(energy,lambda0);
    fhLam1ETRD[pidIndex]->Fill(energy,lambda1);
    fhDispETRD[pidIndex]->Fill(energy,disp);
  }
  
  if(energy < 2){
    fhNCellsLam0LowE[pidIndex] ->Fill(ncells,lambda0);
    fhEtaLam0LowE[pidIndex]    ->Fill(eta,   lambda0);
    fhPhiLam0LowE[pidIndex]    ->Fill(phi,   lambda0);
  }
  else {
    fhNCellsLam0HighE[pidIndex]->Fill(ncells,lambda0);
    fhEtaLam0HighE[pidIndex]   ->Fill(eta,   lambda0);
    fhPhiLam0HighE[pidIndex]   ->Fill(phi,   lambda0);
  }
  
  if(IsDataMC()){
    
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
    
    // Get the fraction of the cluster energy that carries the cell with highest energy
    Int_t absID             =-1 ;
    Float_t maxCellFraction = 0.;
    
    absID = GetCaloUtils()->GetMaxEnergyCell(cells, cluster,maxCellFraction);
    
    // Check the origin and fill histograms
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) && 
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)){
      fhMCELambda0[pidIndex][mcssPhoton]    ->Fill(energy, lambda0);

            
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron && 
              !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion))){
      fhMCELambda0[pidIndex][mcssElectron]    ->Fill(energy, lambda0);
       
      if(!GetReader()->IsEmbeddedClusterSelectionOn()){
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
          printf("AliAnaElectron::FillShowerShapeHistogram() - n overlaps = %d!!", noverlaps);
        }
      }//No embedding
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn()){
        
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
               GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) ){
      fhMCELambda0[pidIndex][mcssConversion]    ->Fill(energy, lambda0);      
    }//conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  ){
      fhMCELambda0[pidIndex][mcssPi0]    ->Fill(energy, lambda0);      
    }//pi0
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)  ){
      fhMCELambda0[pidIndex][mcssEta]    ->Fill(energy, lambda0);
      
    }//eta    
    else {
      fhMCELambda0[pidIndex][mcssOther]    ->Fill(energy, lambda0);      
    }//other particles 
    
  }//MC data
  
}

//________________________________________________________________________
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

//________________________________________________________________________
TList *  AliAnaElectron::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 
	
  Int_t nptbins     = GetHistoPtBins();           Float_t ptmax     = GetHistoPtMax();           Float_t ptmin     = GetHistoPtMin(); 
  Int_t nphibins    = GetHistoPhiBins();          Float_t phimax    = GetHistoPhiMax();          Float_t phimin    = GetHistoPhiMin(); 
  Int_t netabins    = GetHistoEtaBins();          Float_t etamax    = GetHistoEtaMax();          Float_t etamin    = GetHistoEtaMin();	
  Int_t ssbins      = GetHistoShowerShapeBins();  Float_t ssmax     = GetHistoShowerShapeMax();  Float_t ssmin     = GetHistoShowerShapeMin();
  Int_t nbins       = GetHistoNClusterCellBins(); Int_t   nmax      = GetHistoNClusterCellMax(); Int_t   nmin      = GetHistoNClusterCellMin(); 
  Int_t ndedxbins   = GetHistodEdxBins();         Float_t dedxmax   = GetHistodEdxMax();         Float_t dedxmin   = GetHistodEdxMin();
  Int_t nPoverEbins = GetHistoPOverEBins();       Float_t pOverEmax = GetHistoPOverEMax();       Float_t pOverEmin = GetHistoPOverEMin();

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
  
  
  TString pidParticle[] = {"Electron","ChargedHadron"} ;
  
  for(Int_t pidIndex = 0; pidIndex < 2; pidIndex++){
    
    //Shower shape
    if(fFillSSHistograms){
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
      
      if(fCalorimeter == "EMCAL"){
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
      
    } // Shower shape
    
    
    if(IsDataMC()){
      
      if(fFillSSHistograms){
        
        TString ptypess[] = { "#gamma","hadron?","#pi^{0}","#eta","#gamma->e^{#pm}","e^{#pm}"} ; 
        
        TString pnamess[] = { "Photon","Hadron","Pi0","Eta","Conversion","Electron"} ;
        
        for(Int_t i = 0; i < 6; i++){ 
          
          fhMCELambda0[pidIndex][i]  = new TH2F(Form("h%sELambda0_MC%s",pidParticle[pidIndex].Data(),pnamess[i].Data()),
                                                Form("%s like cluster from %s : E vs #lambda_{0}^{2}",pidParticle[pidIndex].Data(),ptypess[i].Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhMCELambda0[pidIndex][i]->SetYTitle("#lambda_{0}^{2}");
          fhMCELambda0[pidIndex][i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhMCELambda0[pidIndex][i]) ; 
          
        }// loop    
      }   
    }
    
    if(IsCaloPIDOn() && pidIndex > 0) continue;
    
    fhNCellsE[pidIndex]  = new TH2F (Form("h%sNCellsE",pidParticle[pidIndex].Data()),
                                     Form("N cells in %s cluster vs E ",pidParticle[pidIndex].Data()),
                                     nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhNCellsE[pidIndex]->SetXTitle("E (GeV)");
    fhNCellsE[pidIndex]->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhNCellsE[pidIndex]);  
    
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
    if(GetMinPt() < 0.5){
      fhEtaPhi05[pidIndex]  = new TH2F(Form("h%sEtaPhi05",pidParticle[pidIndex].Data()),
                                       Form("%s: #eta vs #phi, E > 0.5",pidParticle[pidIndex].Data()),
                                       netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05[pidIndex]->SetYTitle("#phi (rad)");
      fhEtaPhi05[pidIndex]->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05[pidIndex]) ;
    }
    
    
    if(IsDataMC()){
      
      
      TString ptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}", "#pi^{0}","#eta",
        "e^{#pm}","#gamma->e^{#pm}","hadron?","Anti-N","Anti-P"                    } ; 
      
      TString pname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Pi0","Eta","Electron",
        "Conversion", "Hadron", "AntiNeutron","AntiProton"                        } ;
      
      for(Int_t i = 0; i < fNOriginHistograms; i++){ 
                
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
  
  
  if(fFillSSHistograms){
    
    if(IsDataMC()){
      
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
  
  
  //Store calo PID histograms
  if(IsCaloPIDOn()){
    TList * caloPIDHistos = GetCaloPID()->GetCreateOutputObjects() ;
    for(Int_t i = 0; i < caloPIDHistos->GetEntries(); i++) {
      outputContainer->Add(caloPIDHistos->At(i)) ;
    }
    delete caloPIDHistos;
  }
  return outputContainer ;
  
}

//____________________________________________________________________________
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

//____________________________________________________________________________
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

//__________________________________________________________________
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
    if(!ClusterSelected(calo,mom)) continue;
    
    //----------------------------
    //Create AOD for analysis
    //----------------------------
    AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
    
    //...............................................
    //Set the indeces of the original caloclusters (MC, ID), and calorimeter  
    Int_t label = calo->GetLabel();
    aodph.SetLabel(label);
    aodph.SetCaloLabel(calo->GetID(),-1);
    aodph.SetDetector(fCalorimeter);
    //printf("Index %d, Id %d, iaod %d\n",icalo, calo->GetID(),GetOutputAODBranch()->GetEntriesFast());
    
    //...............................................
    //Set bad channel distance bit
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if     (distBad > fMinDist3) aodph.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodph.SetDistToBad(1) ; 
    else                         aodph.SetDistToBad(0) ;
    //printf("DistBad %f Bit %d\n",distBad, aodph.DistToBad());
    
    //--------------------------------------------------------------------------------------
    //Play with the MC stack if available
    //--------------------------------------------------------------------------------------
    
    //Check origin of the candidates
    if(IsDataMC()){
      aodph.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex()));
      
      if(GetDebug() > 0)
        printf("AliAnaElectron::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
    
    
    //-------------------------------------
    //PID selection via dedx
    //-------------------------------------
    Int_t pid  = AliCaloPID::kChargedHadron;
    AliVTrack *track = 0;
    if(!strcmp("AliESDCaloCluster",Form("%s",calo->ClassName()))){
      Int_t iESDtrack = calo->GetTrackMatchedIndex();
      if(iESDtrack<0) printf("AliAnaElectron::MakeAnalysisFillAOD() - Wrong track index\n");
      AliVEvent * event = GetReader()->GetInputEvent();
      track = dynamic_cast<AliVTrack*> (event->GetTrack(iESDtrack));
    }
    else {
      track = dynamic_cast<AliVTrack*>(calo->GetTrackMatched(0));
    }
    
    if(!track) {
      printf("AliAnaElectron::MakeAnalysisFillAOD() - Null track");
      continue;
    }
    
    Float_t dEdx = track->GetTPCsignal();
    fhdEdxvsE->Fill(calo->E(), dEdx);
    fhdEdxvsP->Fill(track->P(),dEdx);
    
    if( dEdx < fdEdxMax && dEdx > fdEdxMin) {
      
      Float_t eOverp = calo->E()/track->P();
      fhEOverPvsE->Fill(calo->E(),  eOverp);
      fhEOverPvsP->Fill(track->P(), eOverp);
      
      if( eOverp < fEOverPMax && eOverp > fEOverPMin) {
        
        pid  = AliCaloPID::kElectron;
        
      } //E/p
      
    }// dEdx
    
    aodph.SetIdentifiedParticleType(pid);    
    
    Int_t pidIndex = 0;// Electron
    if     (pid == AliCaloPID::kElectron)      pidIndex = 0;
    else if(pid == AliCaloPID::kChargedHadron) pidIndex = 1;
    else                                       continue    ;
    
    //--------------------------------------------------------------------------------------
    //Fill some shower shape histograms
    //--------------------------------------------------------------------------------------
    
    FillShowerShapeHistograms(calo,aodph.GetTag(),pid);
    
    //-------------------------------------
    //PID Shower Shape selection or bit setting
    //-------------------------------------
    // Data, PID check on
    if(IsCaloPIDOn()){
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn()){
        if(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E())!=AliCaloPID::kPhoton) continue;
      }
      else{
        if(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E())!=AliCaloPID::kPhoton) continue;
      }
    
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());
      
    }
    
    //    //...............................................
    //    // Data, PID check off
    //    else{
    //      //Set PID bits for later selection (AliAnaPi0 for example)
    //      //GetPDG already called in SetPIDBits.
    //      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph, GetCaloUtils());
    //      if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - PID Bits set \n");		
    //    }
    
    if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetIdentifiedParticleType());
    
    
    //FIXME, this to MakeAnalysisFillHistograms ...
    Int_t absID             = 0; 
    Float_t maxCellFraction = 0;
    AliVCaloCells* cells    = 0;
    
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    absID = GetCaloUtils()->GetMaxEnergyCell(cells, calo,maxCellFraction);
    
    fhMaxCellDiffClusterE[pidIndex]->Fill(aodph.E(),maxCellFraction);
    fhNCellsE[pidIndex]            ->Fill(aodph.E(),calo->GetNCells());
    
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
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
      mcparticles = GetReader()->GetAODMCParticles(0);
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
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
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
      Float_t ptprim  = 0;
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
        ptprim  = primary->Pt();		
        
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
        ptprim  = aodprimary->Pt();
        
      }
      
      Int_t tag =ph->GetTag();
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[pidIndex][mcPhoton])
      {
        fhMCE  [pidIndex][mcPhoton] ->Fill(ecluster);
        fhMCPt [pidIndex][mcPhoton] ->Fill(ptcluster);
        fhMCPhi[pidIndex][mcPhoton] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][mcPhoton] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][mcPhoton]     ->Fill(ecluster, eprim); 
        fhMCDeltaE[pidIndex][mcPhoton] ->Fill(ecluster,eprim-ecluster);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[pidIndex][mcConversion])
        {
          fhMCE  [pidIndex][mcConversion] ->Fill(ecluster);
          fhMCPt [pidIndex][mcConversion] ->Fill(ptcluster);
          fhMCPhi[pidIndex][mcConversion] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][mcConversion] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][mcConversion]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][mcConversion] ->Fill(ecluster,eprim-ecluster);
          
        }			        
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[pidIndex][mcPi0Decay])
        {
          fhMCE  [pidIndex][mcPi0Decay] ->Fill(ecluster);
          fhMCPt [pidIndex][mcPi0Decay] ->Fill(ptcluster);
          fhMCPhi[pidIndex][mcPi0Decay] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][mcPi0Decay] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][mcPi0Decay]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][mcPi0Decay] ->Fill(ecluster,eprim-ecluster);
        }
        else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                  GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[pidIndex][mcOtherDecay])
        {
          fhMCE  [pidIndex][mcOtherDecay] ->Fill(ecluster);
          fhMCPt [pidIndex][mcOtherDecay] ->Fill(ptcluster);
          fhMCPhi[pidIndex][mcOtherDecay] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][mcOtherDecay] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][mcOtherDecay]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][mcOtherDecay] ->Fill(ecluster,eprim-ecluster);          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [pidIndex][mcPi0])
        {
          fhMCE  [pidIndex][mcPi0] ->Fill(ecluster);
          fhMCPt [pidIndex][mcPi0] ->Fill(ptcluster);
          fhMCPhi[pidIndex][mcPi0] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][mcPi0] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][mcPi0]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][mcPi0] ->Fill(ecluster,eprim-ecluster);
          
        } 
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[pidIndex][mcEta])
        {
          fhMCE  [pidIndex][mcEta] ->Fill(ecluster);
          fhMCPt [pidIndex][mcEta] ->Fill(ptcluster);
          fhMCPhi[pidIndex][mcEta] ->Fill(ecluster,phicluster);
          fhMCEta[pidIndex][mcEta] ->Fill(ecluster,etacluster);
          
          fhMC2E[pidIndex][mcEta]     ->Fill(ecluster, eprim);
          fhMCDeltaE[pidIndex][mcEta] ->Fill(ecluster,eprim-ecluster);
          
        }      
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[pidIndex][mcAntiNeutron])
      {
        fhMCE  [pidIndex][mcAntiNeutron] ->Fill(ecluster);
        fhMCPt [pidIndex][mcAntiNeutron] ->Fill(ptcluster);
        fhMCPhi[pidIndex][mcAntiNeutron] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][mcAntiNeutron] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][mcAntiNeutron]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][mcAntiNeutron] ->Fill(ecluster,eprim-ecluster);
        
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[pidIndex][mcAntiProton])
      {
        fhMCE  [pidIndex][mcAntiProton] ->Fill(ecluster);
        fhMCPt [pidIndex][mcAntiProton] ->Fill(ptcluster);
        fhMCPhi[pidIndex][mcAntiProton] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][mcAntiProton] ->Fill(ecluster,etacluster);

        fhMC2E[pidIndex][mcAntiProton]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][mcAntiProton] ->Fill(ecluster,eprim-ecluster);
        
      } 
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[pidIndex][mcElectron])
      {
        fhMCE  [pidIndex][mcElectron] ->Fill(ecluster);
        fhMCPt [pidIndex][mcElectron] ->Fill(ptcluster);
        fhMCPhi[pidIndex][mcElectron] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][mcElectron] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][mcElectron]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][mcElectron] ->Fill(ecluster,eprim-ecluster);
        
      }     
      else if( fhMCE[pidIndex][mcOther]){
        fhMCE  [pidIndex][mcOther] ->Fill(ecluster);
        fhMCPt [pidIndex][mcOther] ->Fill(ptcluster);
        fhMCPhi[pidIndex][mcOther] ->Fill(ecluster,phicluster);
        fhMCEta[pidIndex][mcOther] ->Fill(ecluster,etacluster);
        
        fhMC2E[pidIndex][mcOther]     ->Fill(ecluster, eprim);
        fhMCDeltaE[pidIndex][mcOther] ->Fill(ecluster,eprim-ecluster);
        
      }
      
    }//Histograms with MC
    
  }// aod loop
  
}


//__________________________________________________________________
void AliAnaElectron::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

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
