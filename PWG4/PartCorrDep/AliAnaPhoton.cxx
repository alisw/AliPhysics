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
/* $Id: AliAnaPhoton.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
// Produces input for other analysis classes like AliAnaPi0, 
// AliAnaParticleHadronCorrelation ... 
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>
#include "TParticle.h"
#include "TDatabasePDG.h"

// --- Analysis system --- 
#include "AliAnaPhoton.h" 
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"


ClassImp(AliAnaPhoton)
  
//____________________________________________________________________________
AliAnaPhoton::AliAnaPhoton() : 
    AliAnaPartCorrBaseClass(),    fCalorimeter(""), 
    fMinDist(0.),                 fMinDist2(0.),                fMinDist3(0.), 
    fRejectTrackMatch(0),         fTimeCutMin(-1),              fTimeCutMax(999999),         
    fNCellsCut(0),                fFillSSHistograms(kFALSE), 
    fCheckConversion(kFALSE),     fRemoveConvertedPair(kFALSE), 
    fAddConvertedPairsToAOD(kFALSE), 
    fMassCut(0),                  fConvAsymCut(1.),             fConvDEtaCut(2.),
    fConvDPhiMinCut(-1.),         fConvDPhiMaxCut(7.), 

    // Histograms
    fhNCellsE(0),                 fhEPhoton(0),                 fhPtPhoton(0),  
    fhPhiPhoton(0),               fhEtaPhoton(0), 
    fhEtaPhiPhoton(0),            fhEtaPhi05Photon(0),

    // Conversion histograms
    fhPtPhotonConv(0),            fhEtaPhiPhotonConv(0),        fhEtaPhi05PhotonConv(0),
    fhConvDeltaEta(0),            fhConvDeltaPhi(0),            fhConvDeltaEtaPhi(0), 
    fhConvAsym(0),                fhConvPt(0),
    fhConvDistEta(0),             fhConvDistEn(0),              fhConvDistMass(0),     
    fhConvDistEtaCutEta(0),       fhConvDistEnCutEta(0),        fhConvDistMassCutEta(0),
    fhConvDistEtaCutMass(0),      fhConvDistEnCutMass(0), 
    fhConvDistEtaCutAsy(0),       fhConvDistEnCutAsy(0),

    //Shower shape histograms
    fhDispE(0),                   fhLam0E(0),                   fhLam1E(0), 
    fhdDispE(0),                  fhdLam0E(0),                  fhdLam1E(0), 
    fhDispETRD(0),                fhLam0ETRD(0),                fhLam1ETRD(0),
    fhdDispETRD(0),               fhdLam0ETRD(0),               fhdLam1ETRD(0),

    fhNCellsLam0LowE(0),          fhNCellsLam1LowE(0),          fhNCellsDispLowE(0),  
    fhNCellsLam0HighE(0),         fhNCellsLam1HighE(0),         fhNCellsDispHighE(0),
    fhNCellsdLam0LowE(0),         fhNCellsdLam1LowE(0),         fhNCellsdDispLowE(0),  
    fhNCellsdLam0HighE(0),        fhNCellsdLam1HighE(0),        fhNCellsdDispHighE(0),

    fhEtaLam0LowE(0),             fhPhiLam0LowE(0), 
    fhEtaLam0HighE(0),            fhPhiLam0HighE(0), 
    fhLam0DispLowE(0),            fhLam0DispHighE(0), 
    fhLam1Lam0LowE(0),            fhLam1Lam0HighE(0), 
    fhDispLam1LowE(0),            fhDispLam1HighE(0),
    fhEtadLam0LowE(0),            fhPhidLam0LowE(0), 
    fhEtadLam0HighE(0),           fhPhidLam0HighE(0), 
    fhdLam0dDispLowE(0),          fhdLam0dDispHighE(0), 
    fhdLam1dLam0LowE(0),          fhdLam1dLam0HighE(0), 
    fhdDispdLam1LowE(0),          fhdDispdLam1HighE(0),

    //MC histograms
    fhDeltaE(0),                  fhDeltaPt(0),
    fhRatioE(0),                  fhRatioPt(0),
    fh2E(0),                      fh2Pt(0),

    // Conversion MC histograms
    fhPtConversionTagged(0),           fhPtAntiNeutronTagged(0),       
    fhPtAntiProtonTagged(0),           fhPtUnknownTagged(0),
    fhEtaPhiConversion(0),             fhEtaPhi05Conversion(0),

    fhConvDeltaEtaMCConversion(0),     fhConvDeltaPhiMCConversion(0),  fhConvDeltaEtaPhiMCConversion(0),
    fhConvAsymMCConversion(0),         fhConvPtMCConversion(0),           
    fhConvDispersionMCConversion(0),   fhConvM02MCConversion(0),

    fhConvDeltaEtaMCAntiNeutron(0),    fhConvDeltaPhiMCAntiNeutron(0), fhConvDeltaEtaPhiMCAntiNeutron(0), 
    fhConvAsymMCAntiNeutron(0),        fhConvPtMCAntiNeutron(0), 
    fhConvDispersionMCAntiNeutron(0),  fhConvM02MCAntiNeutron(0),
    fhConvDeltaEtaMCAntiProton(0),     fhConvDeltaPhiMCAntiProton(0),  fhConvDeltaEtaPhiMCAntiProton(0),  
    fhConvAsymMCAntiProton(0),         fhConvPtMCAntiProton(0),  
    fhConvDispersionMCAntiProton(0),   fhConvM02MCAntiProton(0),
    fhConvDeltaEtaMCString(0),         fhConvDeltaPhiMCString(0),      fhConvDeltaEtaPhiMCString(0),      
    fhConvAsymMCString(0),             fhConvPtMCString(0),      
    fhConvDispersionMCString(0),       fhConvM02MCString(0),
    fhConvDistMCConversion(0),         fhConvDistMCConversionCuts(0)
{
  //default ctor
  
  for(Int_t i = 0; i < 12; i++){
    fhPtMC[i]  = 0;
    fhEMC[i]   = 0;
    fhPhiMC[i] = 0;
    fhEtaMC[i] = 0;
  }
  
  for(Int_t i = 0; i < 5; i++){
    fhEMCLambda0[i]     = 0;
    fhEMCLambda1[i]     = 0;
    fhEMCDispersion[i]  = 0;
    fhEMCdLambda0[i]    = 0;
    fhEMCdLambda1[i]    = 0;
    fhEMCdDispersion[i] = 0;
  }
  
  //Initialize parameters
  InitParameters();

}//____________________________________________________________________________
AliAnaPhoton::~AliAnaPhoton() 
{
  //dtor

}

//__________________________________________________________________
Bool_t  AliAnaPhoton::ClusterSelected(AliVCluster* calo, TLorentzVector mom) 
{
  //Select clusters if they pass different cuts
  if(GetDebug() > 2) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
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
  //Skip matched clusters with tracks
  if(fRejectTrackMatch){
    if(IsTrackMatched(calo)) {
      if(GetDebug() > 2) printf("\t Reject track-matched clusters\n");
      return kFALSE ;
    }
    else  
      if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  }// reject matched clusters
  
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
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(), 
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //All checks passed, cluster selected
  return kTRUE;
    
}


//__________________________________________________________________
void  AliAnaPhoton::FillShowerShapeHistograms(AliVCluster* cluster, const Int_t mcTag){
  
  //Fill cluster Shower Shape histograms
  
  if(!fFillSSHistograms || GetMixedEvent()) return;
  
  Float_t energy  = cluster->E();
  Int_t   ncells  = cluster->GetNCells();
  Int_t   ncells2 = ncells*ncells;
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
  
  fhLam0E ->Fill(energy,lambda0);
  fhLam1E ->Fill(energy,lambda1);
  fhDispE ->Fill(energy,disp);
  fhdLam0E->Fill(energy,lambda0/ncells2);
  fhdLam1E->Fill(energy,lambda1/ncells2);
  fhdDispE->Fill(energy,disp/ncells2);
  
  if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5){
    fhLam0ETRD->Fill(energy,lambda0);
    fhLam1ETRD->Fill(energy,lambda1);
    fhDispETRD->Fill(energy,disp);
    fhdLam0ETRD->Fill(energy,lambda0/ncells2);
    fhdLam1ETRD->Fill(energy,lambda1/ncells2);
    fhdDispETRD->Fill(energy,disp/ncells2);
  }
  
  if(energy < 2){
    fhNCellsLam0LowE ->Fill(ncells,lambda0);
    fhNCellsLam1LowE ->Fill(ncells,lambda1);
    fhNCellsDispLowE ->Fill(ncells,disp);
    fhNCellsdLam0LowE->Fill(ncells,lambda0/ncells2);
    fhNCellsdLam1LowE->Fill(ncells,lambda1/ncells2);
    fhNCellsdDispLowE->Fill(ncells,disp/ncells2);
    
    fhLam1Lam0LowE  ->Fill(lambda1,lambda0);
    fhLam0DispLowE  ->Fill(lambda0,disp);
    fhDispLam1LowE  ->Fill(disp,lambda1);
    fhEtaLam0LowE   ->Fill(eta,lambda0);
    fhPhiLam0LowE   ->Fill(phi,lambda0);  
    
    fhdLam1dLam0LowE->Fill(lambda1/ncells2,lambda0/ncells2);
    fhdLam0dDispLowE->Fill(lambda0/ncells2,disp/ncells2);
    fhdDispdLam1LowE->Fill(disp/ncells2,lambda1/ncells2);
    fhEtadLam0LowE  ->Fill(eta,lambda0/ncells2);
    fhPhidLam0LowE  ->Fill(phi,lambda0/ncells2);  
  }
  else {
    fhNCellsLam0HighE ->Fill(ncells,lambda0);
    fhNCellsLam1HighE ->Fill(ncells,lambda1);
    fhNCellsDispHighE ->Fill(ncells,disp);
    fhNCellsdLam0HighE->Fill(ncells,lambda0/ncells2);
    fhNCellsdLam1HighE->Fill(ncells,lambda1/ncells2);
    fhNCellsdDispHighE->Fill(ncells,disp/ncells2);
    
    fhLam1Lam0HighE  ->Fill(lambda1,lambda0);
    fhLam0DispHighE  ->Fill(lambda0,disp);
    fhDispLam1HighE  ->Fill(disp,lambda1);
    fhEtaLam0HighE   ->Fill(eta, lambda0);
    fhPhiLam0HighE   ->Fill(phi, lambda0);
    
    fhdLam1dLam0HighE->Fill(lambda1/ncells2,lambda0/ncells2);
    fhdLam0dDispHighE->Fill(lambda0/ncells2,disp/ncells2);
    fhdDispdLam1HighE->Fill(disp/ncells2,lambda1/ncells2);
    fhEtadLam0HighE  ->Fill(eta, lambda0/ncells2);
    fhPhidLam0HighE  ->Fill(phi, lambda0/ncells2);
  }
  
  if(IsDataMC()){
  
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) && 
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) ){
      fhEMCLambda0[0]    ->Fill(energy, lambda0);
      fhEMCdLambda0[0]   ->Fill(energy, lambda0/ncells2);
      fhEMCLambda1[0]    ->Fill(energy, lambda1);
      fhEMCdLambda1[0]   ->Fill(energy, lambda1/ncells2);
      fhEMCDispersion[0] ->Fill(energy, disp);
      fhEMCdDispersion[0]->Fill(energy, disp/ncells2);                
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)){
      fhEMCLambda0[3]    ->Fill(energy, lambda0);
      fhEMCdLambda0[3]   ->Fill(energy, lambda0/ncells2);
      fhEMCLambda1[3]    ->Fill(energy, lambda1);
      fhEMCdLambda1[3]   ->Fill(energy, lambda1/ncells2);
      fhEMCDispersion[3] ->Fill(energy, disp);
      fhEMCdDispersion[3]->Fill(energy, disp/ncells2);
    }//electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) && 
              GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) ){
      fhEMCLambda0[1]    ->Fill(energy, lambda0);
      fhEMCdLambda0[1]   ->Fill(energy, lambda0/ncells2);
      fhEMCLambda1[1]    ->Fill(energy, lambda1);
      fhEMCdLambda1[1]   ->Fill(energy, lambda1/ncells2);
      fhEMCDispersion[1] ->Fill(energy, disp);
      fhEMCdDispersion[1]->Fill(energy, disp/ncells2);           
    }//conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  ){
      fhEMCLambda0[2]    ->Fill(energy, lambda0);
      fhEMCdLambda0[2]   ->Fill(energy, lambda0/ncells2);
      fhEMCLambda1[2]    ->Fill(energy, lambda1);
      fhEMCdLambda1[2]   ->Fill(energy, lambda1/ncells2);
      fhEMCDispersion[2] ->Fill(energy, disp);
      fhEMCdDispersion[2]->Fill(energy, disp/ncells2);                
    }//pi0
    else {
      fhEMCLambda0[4]    ->Fill(energy, lambda0);
      fhEMCdLambda0[4]   ->Fill(energy, lambda0/ncells2);
      fhEMCLambda1[4]    ->Fill(energy, lambda1);
      fhEMCdLambda1[4]   ->Fill(energy, lambda1/ncells2);
      fhEMCDispersion[4] ->Fill(energy, disp);
      fhEMCdDispersion[4]->Fill(energy, disp/ncells2);    
    }//other particles 
    
  }//MC data
  
}

//________________________________________________________________________
TObjString *  AliAnaPhoton::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPhoton ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRejectTrackMatch: %d\n",fRejectTrackMatch) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"Conversion Selection: fConvAsymCut %1.2f, fConvDEtaCut %1.2f fConvDPhiCut (%1.2f,%1.2f)\n",
           fConvAsymCut, fConvDEtaCut, fConvDPhiMinCut, fConvDPhiMaxCut) ;
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
TList *  AliAnaPhoton::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("PhotonHistos") ; 
	
  Int_t nptbins  = GetHistoPtBins();  Float_t ptmax  = GetHistoPtMax();  Float_t ptmin  = GetHistoPtMin(); 
  Int_t nphibins = GetHistoPhiBins(); Float_t phimax = GetHistoPhiMax(); Float_t phimin = GetHistoPhiMin(); 
  Int_t netabins = GetHistoEtaBins(); Float_t etamax = GetHistoEtaMax(); Float_t etamin = GetHistoEtaMin();	
  Int_t ssbins   = GetHistoShowerShapeBins();  Float_t ssmax = GetHistoShowerShapeMax();  Float_t ssmin = GetHistoShowerShapeMin();

  fhNCellsE  = new TH2F ("hNCellsE","# of cells in cluster vs E of clusters", 100,0, 20, 20,0,20); 
  fhNCellsE->SetXTitle("E (GeV)");
  fhNCellsE->SetYTitle("# of cells in cluster");
  outputContainer->Add(fhNCellsE);  
  
  fhEPhoton  = new TH1F("hEPhoton","Number of #gamma over calorimeter vs energy",nptbins,ptmin,ptmax); 
  fhEPhoton->SetYTitle("N");
  fhEPhoton->SetXTitle("E_{#gamma}(GeV)");
  outputContainer->Add(fhEPhoton) ;   
  
  fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter vs p_{T}",nptbins,ptmin,ptmax); 
  fhPtPhoton->SetYTitle("N");
  fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtPhoton) ; 
  
  fhPhiPhoton  = new TH2F
    ("hPhiPhoton","#phi_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiPhoton->SetYTitle("#phi (rad)");
  fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiPhoton) ; 
  
  fhEtaPhoton  = new TH2F
    ("hEtaPhoton","#eta_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaPhoton->SetYTitle("#eta");
  fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaPhoton) ;
  
  fhEtaPhiPhoton  = new TH2F
  ("hEtaPhiPhoton","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiPhoton->SetYTitle("#phi (rad)");
  fhEtaPhiPhoton->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiPhoton) ;
  if(GetMinPt() < 0.5){
    fhEtaPhi05Photon  = new TH2F
    ("hEtaPhi05Photon","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi05Photon->SetYTitle("#phi (rad)");
    fhEtaPhi05Photon->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05Photon) ;
  }
  
  //Conversion
  if(fCheckConversion){
    
    fhEtaPhiConversion  = new TH2F
    ("hEtaPhiConversion","#eta vs #phi for converted clusters",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiConversion->SetYTitle("#phi (rad)");
    fhEtaPhiConversion->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiConversion) ;
    if(GetMinPt() < 0.5){
      fhEtaPhi05Conversion  = new TH2F
      ("hEtaPhi05Conversion","#eta vs #phi, E > 0.5, for converted clusters",netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05Conversion->SetYTitle("#phi (rad)");
      fhEtaPhi05Conversion->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05Conversion) ;
    }    
    
    fhPtPhotonConv  = new TH1F("hPtPhotonConv","Number of #gamma over calorimeter, conversion",nptbins,ptmin,ptmax); 
    fhPtPhotonConv->SetYTitle("N");
    fhPtPhotonConv->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPhotonConv) ; 
    
    fhEtaPhiPhotonConv  = new TH2F
    ("hEtaPhiPhotonConv","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiPhotonConv->SetYTitle("#phi (rad)");
    fhEtaPhiPhotonConv->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiPhotonConv) ;
    if(GetMinPt() < 0.5){
      fhEtaPhi05PhotonConv  = new TH2F
      ("hEtaPhi05PhotonConv","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05PhotonConv->SetYTitle("#phi (rad)");
      fhEtaPhi05PhotonConv->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05PhotonConv) ;
    }
    
    fhConvDeltaEta  = new TH2F
    ("hConvDeltaEta","#Delta #eta of selected conversion pairs",100,0,fMassCut,netabins*2,-0.5,0.5); 
    fhConvDeltaEta->SetYTitle("#Delta #eta");
    fhConvDeltaEta->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEta) ;
    
    fhConvDeltaPhi  = new TH2F
    ("hConvDeltaPhi","#Delta #phi of selected conversion pairs",100,0,fMassCut,nphibins*2,-0.5,0.5); 
    fhConvDeltaPhi->SetYTitle("#Delta #phi");
    fhConvDeltaPhi->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhi) ;
    
    fhConvDeltaEtaPhi  = new TH2F
    ("hConvDeltaEtaPhi","#Delta #eta vs #Delta #phi of selected conversion pairs",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhi->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhi->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhi) ;
    
    fhConvAsym  = new TH2F
    ("hConvAsym","Asymmetry of selected conversion pairs",100,0,fMassCut,100,0,1); 
    fhConvAsym->SetYTitle("Asymmetry");
    fhConvAsym->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsym) ;  
    
    fhConvPt  = new TH2F
    ("hConvPt","p_{T} of selected conversion pairs",100,0,fMassCut,100,0.,10.); 
    fhConvPt->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPt->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPt) ;
    
    fhConvDistEta  = new TH2F
    ("hConvDistEta","distance to conversion vertex",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEta->SetXTitle("#eta");
    fhConvDistEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEta) ;
    
    fhConvDistEn  = new TH2F
    ("hConvDistEn","distance to conversion vertex",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEn->SetXTitle("E (GeV)");
    fhConvDistEn->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEn) ;

    fhConvDistMass  = new TH2F
    ("hConvDistMass","distance to conversion vertex",100,0,fMassCut,100,0.,5.); 
    fhConvDistMass->SetXTitle("m (GeV/c^2)");
    fhConvDistMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistMass) ;
    
    fhConvDistEtaCutEta  = new TH2F
    ("hConvDistEtaCutEta","distance to conversion vertex, dEta < 0.05",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutEta->SetXTitle("#eta");
    fhConvDistEtaCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutEta) ;
    
    fhConvDistEnCutEta  = new TH2F
    ("hConvDistEnCutEta","distance to conversion vertex, dEta < 0.05",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutEta->SetXTitle("E (GeV)");
    fhConvDistEnCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutEta) ;
    
    fhConvDistMassCutEta  = new TH2F
    ("hConvDistMassCutEta","distance to conversion vertex, dEta < 0.05",100,0,fMassCut,100,0.,5.); 
    fhConvDistMassCutEta->SetXTitle("m (GeV/c^2)");
    fhConvDistMassCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistMassCutEta) ;
    
    fhConvDistEtaCutMass  = new TH2F
    ("hConvDistEtaCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutMass->SetXTitle("#eta");
    fhConvDistEtaCutMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutMass) ;
    
    fhConvDistEnCutMass  = new TH2F
    ("hConvDistEnCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutMass->SetXTitle("E (GeV)");
    fhConvDistEnCutMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutMass) ;

    fhConvDistEtaCutAsy  = new TH2F
    ("hConvDistEtaCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutAsy->SetXTitle("#eta");
    fhConvDistEtaCutAsy->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutAsy) ;
    
    fhConvDistEnCutAsy  = new TH2F
    ("hConvDistEnCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutAsy->SetXTitle("E (GeV)");
    fhConvDistEnCutAsy->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutAsy) ;
    
  } // check conversion
  
  
  //Shower shape
  if(fFillSSHistograms){
    
    fhLam0E  = new TH2F ("hLam0E","#lambda_{0}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLam0E->SetYTitle("#lambda_{0}^{2}");
    fhLam0E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam0E);  
    
    fhLam1E  = new TH2F ("hLam1E","#lambda_{1}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLam1E->SetYTitle("#lambda_{1}^{2}");
    fhLam1E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam1E);  
    
    fhDispE  = new TH2F ("hDispE"," dispersion^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhDispE->SetYTitle("D^{2}");
    fhDispE->SetXTitle("E (GeV) ");
    outputContainer->Add(fhDispE);
    
    fhdLam0E  = new TH2F ("hdLam0E","#lambda_{0}^{2}/N_{cells}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhdLam0E->SetYTitle("d#lambda_{0}^{2}");
    fhdLam0E->SetXTitle("E (GeV)");
    outputContainer->Add(fhdLam0E);  
    
    fhdLam1E  = new TH2F ("hdLam1E","#lambda_{1}^{2}/N_{cells}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhdLam1E->SetYTitle("d#lambda_{1}^{2}");
    fhdLam1E->SetXTitle("E (GeV)");
    outputContainer->Add(fhdLam1E);  
    
    fhdDispE  = new TH2F ("hdDispE"," dispersion^{2}/N_{cells}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhdDispE->SetYTitle("dD^{2}");
    fhdDispE->SetXTitle("E (GeV) ");
    outputContainer->Add(fhdDispE);
    
    
    if(fCalorimeter == "EMCAL"){
      fhLam0ETRD  = new TH2F ("hLam0ETRD","#lambda_{0}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam0ETRD->SetYTitle("#lambda_{0}^{2}");
      fhLam0ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam0ETRD);  
      
      fhLam1ETRD  = new TH2F ("hLam1ETRD","#lambda_{1}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam1ETRD->SetYTitle("#lambda_{1}^{2}");
      fhLam1ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam1ETRD);  
      
      fhDispETRD  = new TH2F ("hDispETRD"," dispersion^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhDispETRD->SetYTitle("Dispersion^{2}");
      fhDispETRD->SetXTitle("E (GeV) ");
      outputContainer->Add(fhDispETRD);   
      
      fhdLam0ETRD  = new TH2F ("hdLam0ETRD","#lambda_{0}^{2}/N_{cells}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
      fhdLam0ETRD->SetYTitle("d#lambda_{0}^{2}");
      fhdLam0ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhdLam0ETRD);  
      
      fhdLam1ETRD  = new TH2F ("hdLam1ETRD","#lambda_{1}^{2}/N_{cells}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
      fhdLam1ETRD->SetYTitle("d#lambda_{1}^{2}");
      fhdLam1ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhdLam1ETRD);  
      
      fhdDispETRD  = new TH2F ("hdDispETRD"," dispersion^{2}/N_{cells}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
      fhdDispETRD->SetYTitle("dD^{2}");
      fhdDispETRD->SetXTitle("E (GeV) ");
      outputContainer->Add(fhdDispETRD);
      
    } 
    
    fhNCellsLam0LowE  = new TH2F ("hNCellsLam0LowE","N_{cells} in cluster vs #lambda_{0}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsLam0LowE->SetXTitle("N_{cells}");
    fhNCellsLam0LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam0LowE);  
    
    fhNCellsLam0HighE  = new TH2F ("hNCellsLam0HighE","N_{cells} in cluster vs #lambda_{0}^{2}, E > 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsLam0HighE->SetXTitle("N_{cells}");
    fhNCellsLam0HighE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam0HighE);  
    
    fhNCellsLam1LowE  = new TH2F ("hNCellsLam1LowE","N_{cells} in cluster vs #lambda_{1}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsLam1LowE->SetXTitle("N_{cells}");
    fhNCellsLam1LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam1LowE);  
    
    fhNCellsLam1HighE  = new TH2F ("hNCellsLam1HighE","N_{cells} in cluster vs #lambda_{1}^{2}, E > 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsLam1HighE->SetXTitle("N_{cells}");
    fhNCellsLam1HighE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam1HighE);  
    
    fhNCellsDispLowE  = new TH2F ("hNCellsDispLowE","N_{cells} in cluster vs dispersion^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsDispLowE->SetXTitle("N_{cells}");
    fhNCellsDispLowE->SetYTitle("D^{2}");
    outputContainer->Add(fhNCellsDispLowE);  
    
    fhNCellsDispHighE  = new TH2F ("hNCellsDispHighE","N_{cells} in cluster vs dispersion^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax); 
    fhNCellsDispHighE->SetXTitle("N_{cells}");
    fhNCellsDispHighE->SetYTitle("D^{2}");
    outputContainer->Add(fhNCellsDispHighE);  
    
    
    fhNCellsdLam0LowE  = new TH2F ("hNCellsdLam0LowE","N_{cells} in cluster vs #lambda_{0}^{2}/N_{cells}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdLam0LowE->SetXTitle("N_{cells}");
    fhNCellsdLam0LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsdLam0LowE);  
    
    fhNCellsdLam0HighE  = new TH2F ("hNCellsdLam0HighE","N_{cells} in cluster vs #lambda_{0}^{2}/N_{cells}^{2}, E > 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdLam0HighE->SetXTitle("N_{cells}");
    fhNCellsdLam0HighE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsdLam0HighE);  
    
    fhNCellsdLam1LowE  = new TH2F ("hNCellsdLam1LowE","N_{cells} in cluster vs #lambda_{1}^{2}/N_{cells}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdLam1LowE->SetXTitle("N_{cells}");
    fhNCellsdLam1LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsdLam1LowE);  
    
    fhNCellsdLam1HighE  = new TH2F ("hNCellsdLam1HighE","N_{cells} in cluster vs #lambda_{1}^{2}/N_{cells}^{2}, E > 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdLam1HighE->SetXTitle("N_{cells}");
    fhNCellsdLam1HighE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsdLam1HighE);  
    
    fhNCellsdDispLowE  = new TH2F ("hNCellsdDispLowE","N_{cells} in cluster vs dispersion^{2}/N_{cells}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdDispLowE->SetXTitle("N_{cells}");
    fhNCellsdDispLowE->SetYTitle("D^{2}");
    outputContainer->Add(fhNCellsdDispLowE);  
    
    fhNCellsdDispHighE  = new TH2F ("hNCellsdDispHighE","N_{cells} in cluster vs dispersion^{2}/N_{cells}^{2}, E < 2 GeV", 20,0, 20, ssbins,ssmin,ssmax/50); 
    fhNCellsdDispHighE->SetXTitle("N_{cells}");
    fhNCellsdDispHighE->SetYTitle("D^{2}");
    outputContainer->Add(fhNCellsdDispHighE);      
    
    
    fhEtaLam0LowE  = new TH2F ("hEtaLam0LowE","#eta vs #lambda_{0}^{2}, E < 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax); 
    fhEtaLam0LowE->SetYTitle("#lambda_{0}^{2}");
    fhEtaLam0LowE->SetXTitle("#eta");
    outputContainer->Add(fhEtaLam0LowE);  
    
    fhPhiLam0LowE  = new TH2F ("hPhiLam0LowE","#phi vs #lambda_{0}^{2}, E < 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
    fhPhiLam0LowE->SetYTitle("#lambda_{0}^{2}");
    fhPhiLam0LowE->SetXTitle("#phi");
    outputContainer->Add(fhPhiLam0LowE);  
    
    fhEtaLam0HighE  = new TH2F ("hEtaLam0HighE","#eta vs #lambda_{0}^{2}, E > 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax); 
    fhEtaLam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhEtaLam0HighE->SetXTitle("#eta");
    outputContainer->Add(fhEtaLam0HighE);  
    
    fhPhiLam0HighE  = new TH2F ("hPhiLam0HighE","#phi vs #lambda_{0}^{2}, E > 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
    fhPhiLam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhPhiLam0HighE->SetXTitle("#phi");
    outputContainer->Add(fhPhiLam0HighE);  
    
    fhLam1Lam0LowE  = new TH2F ("hLam1Lam0LowE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam1Lam0LowE->SetYTitle("#lambda_{0}^{2}");
    fhLam1Lam0LowE->SetXTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhLam1Lam0LowE);  
    
    fhLam1Lam0HighE  = new TH2F ("hLam1Lam0HighE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam1Lam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhLam1Lam0HighE->SetXTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhLam1Lam0HighE);  
    
    fhLam0DispLowE  = new TH2F ("hLam0DispLowE","#lambda_{0}^{2} vs dispersion^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam0DispLowE->SetXTitle("#lambda_{0}^{2}");
    fhLam0DispLowE->SetYTitle("D^{2}");
    outputContainer->Add(fhLam0DispLowE);  
    
    fhLam0DispHighE  = new TH2F ("hLam0DispHighE","#lambda_{0}^{2} vs dispersion^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam0DispHighE->SetXTitle("#lambda_{0}^{2}");
    fhLam0DispHighE->SetYTitle("D^{2}");
    outputContainer->Add(fhLam0DispHighE);  
    
    fhDispLam1LowE  = new TH2F ("hDispLam1LowE","Dispersion^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhDispLam1LowE->SetXTitle("D^{2}");
    fhDispLam1LowE->SetYTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhDispLam1LowE);  
    
    fhDispLam1HighE  = new TH2F ("hDispLam1HighE","Dispersion^{2} vs #lambda_{1^{2}} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhDispLam1HighE->SetXTitle("D^{2}");
    fhDispLam1HighE->SetYTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhDispLam1HighE);  
    
    
    
    fhEtadLam0LowE  = new TH2F ("hEtadLam0LowE","#eta vs #lambda_{0}^{2}/N_{cells}^{2}, E < 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax/50); 
    fhEtadLam0LowE->SetYTitle("d#lambda_{0}^{2}");
    fhEtadLam0LowE->SetXTitle("#eta");
    outputContainer->Add(fhEtadLam0LowE);  
    
    fhPhidLam0LowE  = new TH2F ("hPhidLam0LowE","#phi vs #lambda_{0}^{2}/N_{cells}^{2}, E < 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax/50); 
    fhPhidLam0LowE->SetYTitle("d#lambda_{0}^{2}");
    fhPhidLam0LowE->SetXTitle("#phi");
    outputContainer->Add(fhPhidLam0LowE);  
    
    fhEtadLam0HighE  = new TH2F ("hEtadLam0HighE","#eta vs #lambda_{0}^{2}/N_{cells}^{2}", netabins,etamin,etamax, ssbins,ssmin,ssmax/50); 
    fhEtadLam0HighE->SetYTitle("d#lambda_{0}^{2}");
    fhEtadLam0HighE->SetXTitle("#eta");
    outputContainer->Add(fhEtadLam0HighE);  
    
    fhPhidLam0HighE  = new TH2F ("hPhidLam0HighE","#phi vs #lambda_{0}^{2}/N_{cells}^{2}, E > 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax/50); 
    fhPhidLam0HighE->SetYTitle("d#lambda_{0}^{2}");
    fhPhidLam0HighE->SetXTitle("#phi");
    outputContainer->Add(fhPhidLam0HighE);  
    
    fhdLam1dLam0LowE  = new TH2F ("hdLam1dLam0LowE","#lambda_{0}^{2}/N_{cells}^{2} vs #lambda_{1}^{2}/N_{cells}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdLam1dLam0LowE->SetYTitle("d#lambda_{0}^{2}");
    fhdLam1dLam0LowE->SetXTitle("d#lambda_{1}^{2}");
    outputContainer->Add(fhdLam1dLam0LowE);  
    
    fhdLam1dLam0HighE  = new TH2F ("hdLam1dLam0HighE","#lambda_{0}^{2}/N_{cells}^{2}vs #lambda_{1}^{2}/N_{cells}^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdLam1dLam0HighE->SetYTitle("d#lambda_{0}^{2}");
    fhdLam1dLam0HighE->SetXTitle("d#lambda_{1}^{2}");
    outputContainer->Add(fhdLam1dLam0HighE);  
    
    fhdLam0dDispLowE  = new TH2F ("hdLam0dDispLowE","#lambda_{0}^{2}/N_{cells}^{2} vs dispersion^{2}/N_{cells}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdLam0dDispLowE->SetXTitle("d#lambda_{0}^{2}");
    fhdLam0dDispLowE->SetYTitle("dD^{2}");
    outputContainer->Add(fhdLam0dDispLowE);  
    
    fhdLam0dDispHighE  = new TH2F ("hdLam0dDispHighE","#lambda_{0}^{2}/N_{cells}^{2} vs dispersion^{2}/N_{cells}^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdLam0dDispHighE->SetXTitle("d#lambda_{0}^{2}");
    fhdLam0dDispHighE->SetYTitle("dD^{2}");
    outputContainer->Add(fhdLam0dDispHighE);  
    
    fhdDispdLam1LowE  = new TH2F ("hdDispddLam1LowE","Dispersion^{2}/N_{cells}^{2} vs #lambda_{1}^{2}/N_{cells}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdDispdLam1LowE->SetXTitle("dD^{2}");
    fhdDispdLam1LowE->SetYTitle("d#lambda_{1}^{2}");
    outputContainer->Add(fhdDispdLam1LowE);  
    
    fhdDispdLam1HighE  = new TH2F ("hdDispdLam1HighE","Dispersion^{2}/N_{cells}^{2} vs #lambda_{1^{2}}/N_{cells}^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax/50, ssbins,ssmin,ssmax/50); 
    fhdDispdLam1HighE->SetXTitle("dD^{2}");
    fhdDispdLam1HighE->SetYTitle("d#lambda_{1}^{2}");
    outputContainer->Add(fhdDispdLam1HighE);  
    
    
  } // Shower shape
  
  
  if(IsDataMC()){
    fhDeltaE  = new TH1F ("hDeltaE","MC - Reco E ", 200,-50,50); 
    fhDeltaE->SetXTitle("#Delta E (GeV)");
    outputContainer->Add(fhDeltaE);
                
    fhDeltaPt  = new TH1F ("hDeltaPt","MC - Reco p_{T} ", 200,-50,50); 
    fhDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
    outputContainer->Add(fhDeltaPt);

    fhRatioE  = new TH1F ("hRatioE","Reco/MC E ", 200,0,2); 
    fhRatioE->SetXTitle("E_{reco}/E_{gen}");
    outputContainer->Add(fhRatioE);
    
    fhRatioPt  = new TH1F ("hRatioPt","Reco/MC p_{T} ", 200,0,2); 
    fhRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
    outputContainer->Add(fhRatioPt);    

    fh2E  = new TH2F ("h2E","E distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2E->SetXTitle("E_{rec} (GeV)");
    fh2E->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fh2E);          
    
    fh2Pt  = new TH2F ("h2Pt","p_T distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2Pt->SetXTitle("p_{T,rec} (GeV/c)");
    fh2Pt->SetYTitle("p_{T,gen} (GeV/c)");
    outputContainer->Add(fh2Pt);
   
    TString ptype[] ={"#gamma","#gamma->e^{#pm}","#pi^{0}","e^{#pm}", "hadron",
      "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}","#gamma_{#pi decay}",
      "Anti-N","Anti-P","String"}; 
    TString pname[] ={"Photon","Conversion","Pi0",    "Electron","Hadron", 
      "PhotonPrompt","PhotonFragmentation","PhotonISR","PhotonPi0Decay",
      "AntiNeutron","AntiProton","String"};
    
    for(Int_t i = 0; i < 12; i++){ 
      
      fhEMC[i]  = new TH1F(Form("hE_MC%s",pname[i].Data()),
                                Form("cluster from %s : E ",ptype[i].Data()),
                                nptbins,ptmin,ptmax); 
      fhEMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEMC[i]) ; 
      
      fhPtMC[i]  = new TH1F(Form("hPt_MC%s",pname[i].Data()),
                           Form("cluster from %s : p_{T} ",ptype[i].Data()),
                           nptbins,ptmin,ptmax); 
      fhPtMC[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtMC[i]) ;
      
      fhEtaMC[i]  = new TH2F(Form("hEta_MC%s",pname[i].Data()),
                           Form("cluster from %s : #eta ",ptype[i].Data()),
                           nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaMC[i]->SetYTitle("#eta");
      fhEtaMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEtaMC[i]) ;
      
      fhPhiMC[i]  = new TH2F(Form("hPhi_MC%s",pname[i].Data()),
                           Form("cluster from %s : #phi ",ptype[i].Data()),
                           nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiMC[i]->SetYTitle("#phi (rad)");
      fhPhiMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhPhiMC[i]) ;
      
    }
    	
    if(fCheckConversion){  
      fhPtConversionTagged  = new TH1F("hPtMCConversionTagged","Number of converted #gamma over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtConversionTagged->SetYTitle("N");
      fhPtConversionTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtConversionTagged) ; 
      
      
      fhPtAntiNeutronTagged  = new TH1F("hPtMCAntiNeutronTagged","Number of AntiNeutron id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtAntiNeutronTagged->SetYTitle("N");
      fhPtAntiNeutronTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtAntiNeutronTagged) ; 
      
      fhPtAntiProtonTagged  = new TH1F("hPtMCAntiProtonTagged","Number of AntiProton id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtAntiProtonTagged->SetYTitle("N");
      fhPtAntiProtonTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtAntiProtonTagged) ; 
      
      fhPtUnknownTagged  = new TH1F("hPtMCUnknownTagged","Number of Unknown id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtUnknownTagged->SetYTitle("N");
      fhPtUnknownTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtUnknownTagged) ;     
      
      fhConvDeltaEtaMCConversion  = new TH2F
      ("hConvDeltaEtaMCConversion","#Delta #eta of selected conversion pairs from real conversions",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCConversion->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCConversion) ;
      
      fhConvDeltaPhiMCConversion  = new TH2F
      ("hConvDeltaPhiMCConversion","#Delta #phi of selected conversion pairs from real conversions",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCConversion->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCConversion) ;
      
      fhConvDeltaEtaPhiMCConversion  = new TH2F
      ("hConvDeltaEtaPhiMCConversion","#Delta #eta vs #Delta #phi of selected conversion pairs, from real conversions",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCConversion->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCConversion->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCConversion) ;
      
      fhConvAsymMCConversion  = new TH2F
      ("hConvAsymMCConversion","Asymmetry of selected conversion pairs from real conversions",100,0,fMassCut,100,0,1); 
      fhConvAsymMCConversion->SetYTitle("Asymmetry");
      fhConvAsymMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCConversion) ;
      
      fhConvPtMCConversion  = new TH2F
      ("hConvPtMCConversion","p_{T} of selected conversion pairs from real conversions",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCConversion->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCConversion) ;    
      
      fhConvDispersionMCConversion  = new TH2F
      ("hConvDispersionMCConversion","p_{T} of selected conversion pairs from real conversions",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCConversion->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCConversion->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCConversion) ;   
      
      fhConvM02MCConversion  = new TH2F
      ("hConvM02MCConversion","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCConversion->SetYTitle("M02 cluster 1");
      fhConvM02MCConversion->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCConversion) ;           
      
      fhConvDeltaEtaMCAntiNeutron  = new TH2F
      ("hConvDeltaEtaMCAntiNeutron","#Delta #eta of selected conversion pairs from anti-neutrons",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCAntiNeutron->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCAntiNeutron) ;
      
      fhConvDeltaPhiMCAntiNeutron  = new TH2F
      ("hConvDeltaPhiMCAntiNeutron","#Delta #phi of selected conversion pairs from anti-neutrons",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCAntiNeutron) ;
      
      fhConvDeltaEtaPhiMCAntiNeutron  = new TH2F
      ("hConvDeltaEtaPhiMCAntiNeutron","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-neutrons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCAntiNeutron->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCAntiNeutron) ;    
      
      fhConvAsymMCAntiNeutron  = new TH2F
      ("hConvAsymMCAntiNeutron","Asymmetry of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0,1); 
      fhConvAsymMCAntiNeutron->SetYTitle("Asymmetry");
      fhConvAsymMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCAntiNeutron) ;
      
      fhConvPtMCAntiNeutron  = new TH2F
      ("hConvPtMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCAntiNeutron->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCAntiNeutron) ;    
      
      fhConvDispersionMCAntiNeutron  = new TH2F
      ("hConvDispersionMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCAntiNeutron->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCAntiNeutron->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCAntiNeutron) ;       
      
      fhConvM02MCAntiNeutron  = new TH2F
      ("hConvM02MCAntiNeutron","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCAntiNeutron->SetYTitle("M02 cluster 1");
      fhConvM02MCAntiNeutron->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCAntiNeutron) ;  
      
      fhConvDeltaEtaMCAntiProton  = new TH2F
      ("hConvDeltaEtaMCAntiProton","#Delta #eta of selected conversion pairs from anti-protons",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCAntiProton->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCAntiProton) ;
      
      fhConvDeltaPhiMCAntiProton  = new TH2F
      ("hConvDeltaPhiMCAntiProton","#Delta #phi of selected conversion pairs from anti-protons",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCAntiProton->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCAntiProton) ;
      
      fhConvDeltaEtaPhiMCAntiProton  = new TH2F
      ("hConvDeltaEtaPhiMCAntiProton","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-protons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCAntiProton->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCAntiProton->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCAntiProton) ;    
      
      fhConvAsymMCAntiProton  = new TH2F
      ("hConvAsymMCAntiProton","Asymmetry of selected conversion pairs from anti-protons",100,0,fMassCut,100,0,1); 
      fhConvAsymMCAntiProton->SetYTitle("Asymmetry");
      fhConvAsymMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCAntiProton) ;
      
      fhConvPtMCAntiProton  = new TH2F
      ("hConvPtMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCAntiProton->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCAntiProton) ;
      
      fhConvDispersionMCAntiProton  = new TH2F
      ("hConvDispersionMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCAntiProton->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCAntiProton->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCAntiProton) ;       
      
      fhConvM02MCAntiProton  = new TH2F
      ("hConvM02MCAntiProton","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCAntiProton->SetYTitle("M02 cluster 1");
      fhConvM02MCAntiProton->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCAntiProton) ;       
      
      fhConvDeltaEtaMCString  = new TH2F
      ("hConvDeltaEtaMCString","#Delta #eta of selected conversion pairs from string",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCString->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCString) ;
      
      fhConvDeltaPhiMCString  = new TH2F
      ("hConvDeltaPhiMCString","#Delta #phi of selected conversion pairs from string",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCString->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCString) ;
      
      fhConvDeltaEtaPhiMCString  = new TH2F
      ("hConvDeltaEtaPhiMCString","#Delta #eta vs #Delta #phi of selected conversion pairs from string",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCString->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCString->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCString) ;    
      
      fhConvAsymMCString  = new TH2F
      ("hConvAsymMCString","Asymmetry of selected conversion pairs from string",100,0,fMassCut,100,0,1); 
      fhConvAsymMCString->SetYTitle("Asymmetry");
      fhConvAsymMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCString) ;
      
      fhConvPtMCString  = new TH2F
      ("hConvPtMCString","p_{T} of selected conversion pairs from string",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCString->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCString) ;
      
      fhConvDispersionMCString  = new TH2F
      ("hConvDispersionMCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCString->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCString->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCString) ;       
      
      fhConvM02MCString  = new TH2F
      ("hConvM02MCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCString->SetYTitle("M02 cluster 1");
      fhConvM02MCString->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCString) ; 
      
      fhConvDistMCConversion  = new TH2F
      ("hConvDistMCConversion","calculated conversion distance vs real vertes for MC conversion",100,0.,5.,100,0.,5.); 
      fhConvDistMCConversion->SetYTitle("distance");
      fhConvDistMCConversion->SetXTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversion) ; 
      
      fhConvDistMCConversionCuts  = new TH2F
      ("hConvDistMCConversionCuts","calculated conversion distance vs real vertes for MC conversion, deta < 0.05, m < 10 MeV, asym < 0.1",100,0.,5.,100,0.,5.); 
      fhConvDistMCConversionCuts->SetYTitle("distance");
      fhConvDistMCConversionCuts->SetXTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversionCuts) ; 
    }
    
    if(fFillSSHistograms){
      
      for(Int_t i = 0; i < 5; i++){ 
        
        fhEMCLambda0[i]  = new TH2F(Form("hELambda0_MC%s",pname[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{0}^{2}",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEMCLambda0[i]->SetYTitle("#lambda_{0}^{2}");
        fhEMCLambda0[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCLambda0[i]) ; 
        
        fhEMCdLambda0[i]  = new TH2F(Form("hEdLambda0_MC%s",pname[i].Data()),
                                     Form("cluster from %s : E vs #lambda_{0}^{2}/N_{cells}^{2}",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
        fhEMCdLambda0[i]->SetYTitle("d#lambda_{0}^{2}");
        fhEMCdLambda0[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCdLambda0[i]) ; 
        
        fhEMCLambda1[i]  = new TH2F(Form("hELambda1_MC%s",pname[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{1}^{2}",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEMCLambda1[i]->SetYTitle("#lambda_{1}^{2}");
        fhEMCLambda1[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCLambda1[i]) ; 
        
        fhEMCdLambda1[i]  = new TH2F(Form("hEdLambda1_MC%s",pname[i].Data()),
                                     Form("cluster from %s : E vs d#lambda_{1}^{2}/N_{cells}^{2}",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
        fhEMCdLambda1[i]->SetYTitle("d#lambda_{1}^{2}");
        fhEMCdLambda1[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCdLambda1[i]) ; 
        
        fhEMCDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pname[i].Data()),
                                       Form("cluster from %s : E vs dispersion^{2}",ptype[i].Data()),
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEMCDispersion[i]->SetYTitle("D^{2}");
        fhEMCDispersion[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCDispersion[i]) ; 
        
        fhEMCdDispersion[i]  = new TH2F(Form("hEdDispersion_MC%s",pname[i].Data()),
                                        Form("cluster from %s : E vs dispersion^{2}/N_{cells}^{2}",ptype[i].Data()),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
        fhEMCdDispersion[i]->SetYTitle("dD^{2}");
        fhEMCdDispersion[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhEMCdDispersion[i]) ;   
        
      }// loop       
      
    }// Fill SS MC histograms
    
  }//Histos with MC
    
  //Store calo PID histograms
  if(fRejectTrackMatch){
    TList * caloPIDHistos = GetCaloPID()->GetCreateOutputObjects() ;
    for(Int_t i = 0; i < caloPIDHistos->GetEntries(); i++) {
      outputContainer->Add(caloPIDHistos->At(i)) ;
    }
    delete caloPIDHistos;
  }
  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaPhoton::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPhoton::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPhoton::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//____________________________________________________________________________
void AliAnaPhoton::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaPhoton_");
  
  fCalorimeter = "EMCAL" ;
  fMinDist     = 2.;
  fMinDist2    = 4.;
  fMinDist3    = 5.;
  fMassCut     = 0.03; //30 MeV
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCut   = 0;
	
  fRejectTrackMatch       = kTRUE ;
  fCheckConversion        = kFALSE;
  fRemoveConvertedPair    = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
	
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillAOD() 
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
  Bool_t * indexConverted = 0x0;
  if(fCheckConversion){
    indexConverted = new Bool_t[nCaloClusters];
    for (Int_t i = 0; i < nCaloClusters; i++) 
      indexConverted[i] = kFALSE;
	}
  
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  
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
        printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
    
    //--------------------------------------------------------------------------------------
    //Fill some shower shape histograms before PID is applied
    //--------------------------------------------------------------------------------------
    
    FillShowerShapeHistograms(calo,aodph.GetTag());
    
    //-------------------------------------
    //PID selection or bit setting
    //-------------------------------------
    // MC
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());	 
      //If primary is not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;
    }	
    //...............................................
    // Data, PID check on
    else if(IsCaloPIDOn()){
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn())
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      else
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;			
      
    }
    //...............................................
    // Data, PID check off
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph, GetCaloUtils());
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetIdentifiedParticleType());
    
    
    //--------------------------------------------------------------------------------------
    // Conversions pairs analysis
    // Check if cluster comes from a conversion in the material in front of the calorimeter
    // Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
    //--------------------------------------------------------------------------------------
    
    // Do analysis only if there are more than one cluster
    if( nCaloClusters > 1 && fCheckConversion){
      Bool_t bConverted = kFALSE;
      Int_t id2 = -1;
		  
      //Check if set previously as converted couple, if so skip its use.
      if (indexConverted[icalo]) continue;
		  
      // Second cluster loop
      for(Int_t jcalo = icalo + 1 ; jcalo < nCaloClusters ; jcalo++) {
        //Check if set previously as converted couple, if so skip its use.
        if (indexConverted[jcalo]) continue;
        //printf("Check Conversion indeces %d and %d\n",icalo,jcalo);
        AliVCluster * calo2 =  (AliVCluster*) (pl->At(jcalo));  //Get cluster kinematics
        
        
        //Mixed event, get index of event
        Int_t evtIndex2 = 0 ; 
        if (GetMixedEvent()) {
          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo2->GetID()) ; 
          
        }      
        
        //Get kinematics of second cluster
        if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
          calo2->GetMomentum(mom2,GetVertex(evtIndex2)) ;}//Assume that come from vertex in straight line
        else{
          Double_t vertex[]={0,0,0};
          calo2->GetMomentum(mom2,vertex) ;
        }
        
        //--------------------------------------
        // Cluster selection
        //--------------------------------------
        
        if(!ClusterSelected(calo2,mom2)) continue;  
        
        //................................................
        // Get TOF of each cluster in pair, calculate difference if small, 
        // take this pair. Way to reject clusters from hadrons (or pileup?)
        Double_t t12diff = calo2->GetTOF()-calo->GetTOF()*1e9;
        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
        
        //................................................
        //Get mass of pair, if small, take this pair.
        Float_t pairM     = (mom+mom2).M();
        //printf("\t both in calo, mass %f, cut %f\n",pairM,fMassCut);
        if(pairM < fMassCut){  
          aodph.SetTagged(kFALSE);
          id2 = calo2->GetID();
          indexConverted[icalo]=kTRUE;
          indexConverted[jcalo]=kTRUE; 
          Float_t asymmetry = TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E());
          Float_t dPhi      = mom.Phi()-mom2.Phi();
          Float_t dEta      = mom.Eta()-mom2.Eta();  
          
          //...............................................
          //Fill few histograms with kinematics of the pair
          //FIXME, move all this to MakeAnalysisFillHistograms ...
          
          fhConvDeltaEta   ->Fill( pairM, dPhi      );
          fhConvDeltaPhi   ->Fill( pairM, dEta      );
          fhConvAsym       ->Fill( pairM, asymmetry );
          fhConvDeltaEtaPhi->Fill( dEta , dPhi      );
          fhConvPt         ->Fill( pairM, (mom+mom2).Pt());          
          
          //Estimate conversion distance, T. Awes, M. Ivanov
          //Under the assumption that the pair has zero mass, and that each electron 
          //of the pair has the same momentum, they will each have the same bend radius 
          //given by R=p/(qB) = p / (300 B) with p in [MeV/c], B in [Tesla] and R in [m]. 
          //With nominal ALICE magnet current of 30kA B=0.5T, and so with E_cluster=p,  
          //R = E/1.5 [cm].  Under these assumptions, the distance from the conversion 
          //point to the EMCal can be related to the separation distance, L=2y, on the EMCal 
          //as d = sqrt(R^2 -(R-y)^2) = sqrt(2Ry - y^2). And since R>>y we can write as 
          //d = sqrt(E*L/1.5) where E is the cluster energy and L is the distance in cm between 
          //the clusters.
          Float_t pos1[3];
          calo->GetPosition(pos1); 
          Float_t pos2[3];
          calo2->GetPosition(pos2); 
          Float_t clustDist = TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+
                                          (pos1[1]-pos2[1])*(pos1[1]-pos2[1])+
                                          (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
          
          Float_t convDist  = TMath::Sqrt(mom.E() *clustDist*0.01/0.15);
          Float_t convDist2 = TMath::Sqrt(mom2.E()*clustDist*0.01/0.15);
          //printf("l = %f, e1 = %f, d1=%f, e2 = %f, d2=%f\n",clustDist,mom.E(),convDist,mom2.E(),convDist2);
          if(GetDebug() > 2)
            printf("AliAnaPhoton::MakeAnalysisFillAOD(): Pair with mass %2.3f < %2.3f, %1.2f < dPhi %2.2f < %2.2f, dEta %f < %2.2f, asymmetry %2.2f< %2.2f; \n    cluster1 id %d, e %2.3f  SM %d, eta %2.3f, phi %2.3f ; \n    cluster2 id %d, e %2.3f, SM %d,eta %2.3f, phi %2.3f\n",
                   pairM,fMassCut,fConvDPhiMinCut, dPhi, fConvDPhiMaxCut, dEta, fConvDEtaCut, asymmetry, fConvAsymCut,
                   calo->GetID(),calo->E(),GetCaloUtils()->GetModuleNumber(calo), mom.Eta(), mom.Phi(),
                   id2, calo2->E(), GetCaloUtils()->GetModuleNumber(calo2),mom2.Eta(), mom2.Phi());
          
          fhConvDistEta ->Fill(mom .Eta(), convDist );
          fhConvDistEta ->Fill(mom2.Eta(), convDist2);
          fhConvDistEn  ->Fill(mom .E(),   convDist );
          fhConvDistEn  ->Fill(mom2.E(),   convDist2);        
          fhConvDistMass->Fill((mom+mom2).M(), convDist );
          //dEta cut
          if(dEta<0.05){
            fhConvDistEtaCutEta ->Fill(mom .Eta(), convDist );
            fhConvDistEtaCutEta ->Fill(mom2.Eta(), convDist2);
            fhConvDistEnCutEta  ->Fill(mom .E(),   convDist );
            fhConvDistEnCutEta  ->Fill(mom2.E(),   convDist2);        
            fhConvDistMassCutEta->Fill((mom+mom2).M(), convDist );
            //mass cut
            if(pairM<0.01){//10 MeV
              fhConvDistEtaCutMass ->Fill(mom .Eta(), convDist );
              fhConvDistEtaCutMass ->Fill(mom2.Eta(), convDist2);
              fhConvDistEnCutMass  ->Fill(mom .E(),   convDist );
              fhConvDistEnCutMass  ->Fill(mom2.E(),   convDist2);        
              // asymmetry cut
              if(asymmetry<0.1){
                fhConvDistEtaCutAsy ->Fill(mom .Eta(), convDist );
                fhConvDistEtaCutAsy ->Fill(mom2.Eta(), convDist2);
                fhConvDistEnCutAsy  ->Fill(mom .E(),   convDist );
                fhConvDistEnCutAsy  ->Fill(mom2.E(),   convDist2); 
              }//asymmetry cut
            }//mass cut            
          }//dEta cut
          
          //...............................................
          //Select pairs in a eta-phi window
          if(TMath::Abs(dEta) < fConvDEtaCut    && 
             TMath::Abs(dPhi) < fConvDPhiMaxCut &&
             TMath::Abs(dPhi) > fConvDPhiMinCut && 
             asymmetry        < fConvAsymCut       ){
            bConverted = kTRUE;          
          }
          //printf("Accepted? %d\n",bConverted);
          //...........................................
          //Fill more histograms, simulated data
          //FIXME, move all this to MakeAnalysisFillHistograms ...
          if(IsDataMC()){
            
            //Check the origin of the pair, look for conversion, antinucleons or jet correlations (strings)
            Int_t ancPDG    = 0;
            Int_t ancStatus = 0;
            TLorentzVector momentum;
            TVector3 prodVertex;
            Int_t ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(calo->GetLabel(), calo2->GetLabel(), 
                                                                        GetReader(), ancPDG, ancStatus, momentum, prodVertex);
            
            // printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
            //                          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
            
            Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(calo2->GetLabels(),calo2->GetNLabels(),GetReader(), 0);
            if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCConversion)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && (ancPDG==22 || TMath::Abs(ancPDG)==11) && ancLabel > -1){
                fhConvDeltaEtaMCConversion   ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCConversion   ->Fill( pairM, dPhi      );
                fhConvAsymMCConversion       ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCConversion->Fill( dEta , dPhi      );
                fhConvPtMCConversion         ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCConversion ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCConversion        ->Fill( calo->GetM02(), calo2->GetM02());
                fhConvDistMCConversion       ->Fill( convDist , prodVertex.Mag() );
                fhConvDistMCConversion       ->Fill( convDist2, prodVertex.Mag() );
                
                if(dEta<0.05 && pairM<0.01 && asymmetry<0.1){
                  fhConvDistMCConversionCuts->Fill( convDist , prodVertex.Mag() );
                  fhConvDistMCConversionCuts->Fill( convDist2, prodVertex.Mag() );
                }
                
              }              
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiNeutron)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiNeutron) && ancPDG==-2112 && ancLabel > -1){
                fhConvDeltaEtaMCAntiNeutron    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiNeutron    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiNeutron        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiNeutron ->Fill( dEta , dPhi      );
                fhConvPtMCAntiNeutron          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiNeutron  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiNeutron         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiProton)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiProton) && ancPDG==-2212 && ancLabel > -1){
                fhConvDeltaEtaMCAntiProton    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiProton    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiProton        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiProton ->Fill( dEta , dPhi      );
                fhConvPtMCAntiProton          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiProton  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiProton         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            
            //Pairs coming from fragmenting pairs.
            if(ancPDG < 22 && ancLabel > 7 && (ancStatus == 11 || ancStatus == 12) ){
              fhConvDeltaEtaMCString    ->Fill( pairM, dPhi);
              fhConvDeltaPhiMCString    ->Fill( pairM, dPhi);
              fhConvAsymMCString        ->Fill( pairM, TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E()) );
              fhConvDeltaEtaPhiMCString ->Fill( dPhi, dPhi );
              fhConvPtMCString          ->Fill( pairM, (mom+mom2).Pt());
              fhConvDispersionMCString  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
              fhConvM02MCString         ->Fill( calo->GetM02(), calo2->GetM02());
            }
            
          }// Data MC
          
          break;
        }
			  
      }//Mass loop
		  
      //..........................................................................................................
      //Pair selected as converted, remove both clusters or recombine them into a photon and put them in the AOD
      if(bConverted){ 
        //Add to AOD
        if(fAddConvertedPairsToAOD){
          //Create AOD of pair analysis
          TLorentzVector mpair = mom+mom2;
          AliAODPWG4Particle aodpair = AliAODPWG4Particle(mpair);
          aodpair.SetLabel(aodph.GetLabel());
          //aodpair.SetInputFileIndex(input);
          
          //printf("Index %d, Id %d\n",icalo, calo->GetID());
          //Set the indeces of the original caloclusters  
          aodpair.SetCaloLabel(calo->GetID(),id2);
          aodpair.SetDetector(fCalorimeter);
          aodpair.SetIdentifiedParticleType(aodph.GetIdentifiedParticleType());
          aodpair.SetTag(aodph.GetTag());
          aodpair.SetTagged(kTRUE);
          //Add AOD with pair object to aod branch
          AddAODParticle(aodpair);
          //printf("\t \t both added pair\n");
        }
        
        //Do not add the current calocluster
        if(fRemoveConvertedPair) continue;
        else {
          //printf("TAGGED\n");
          //Tag this cluster as likely conversion
          aodph.SetTagged(kTRUE);
        }
      }//converted pair
    }//check conversion
    //printf("\t \t added single cluster %d\n",icalo);
	  
    //FIXME, this to MakeAnalysisFillHistograms ...
    fhNCellsE->Fill(aodph.E(),calo->GetNCells());
    
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  delete [] indexConverted;
	
  if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms() 
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
	printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
	abort();
      }
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      //Get the list of MC particles
      mcparticles = GetReader()->GetAODMCParticles(0);
      if(!mcparticles && GetDebug() > 0) 	{
	printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
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
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ph->GetIdentifiedParticleType();
    
    if(GetDebug() > 3) 
      printf("AliAnaPhoton::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ph->GetIdentifiedParticleType(),ph->GetTag(), (ph->GetDetector()).Data()) ;
    
    //If PID used, fill histos with photons in Calorimeter fCalorimeter
    if(IsCaloPIDOn() && pdg != AliCaloPID::kPhoton) continue; 
    if(ph->GetDetector() != fCalorimeter) continue;
    
    if(GetDebug() > 2) 
      printf("AliAnaPhoton::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
    
    //................................
    //Fill photon histograms 
    Float_t ptcluster  = ph->Pt();
    Float_t phicluster = ph->Phi();
    Float_t etacluster = ph->Eta();
    Float_t ecluster   = ph->E();
    
    fhEPhoton   ->Fill(ecluster);
    fhPtPhoton  ->Fill(ptcluster);
    fhPhiPhoton ->Fill(ptcluster,phicluster);
    fhEtaPhoton ->Fill(ptcluster,etacluster);    
    if     (ecluster > 0.5)   fhEtaPhiPhoton  ->Fill(etacluster, phicluster);
    else if(GetMinPt() < 0.5) fhEtaPhi05Photon->Fill(etacluster, phicluster);
    
    if(fCheckConversion &&ph->IsTagged()){
      fhPtPhotonConv->Fill(ptcluster);
      if(ecluster > 0.5)        fhEtaPhiPhotonConv  ->Fill(etacluster, phicluster);
      else if(GetMinPt() < 0.5) fhEtaPhi05PhotonConv->Fill(etacluster, phicluster);
    }
    
    //.......................................
    //Play with the MC data if available
    if(IsDataMC()){
      
      Int_t tag =ph->GetTag();
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
	{
	  fhEMC  [mcPhoton] ->Fill(ecluster);
	  fhPtMC [mcPhoton] ->Fill(ptcluster);
	  fhPhiMC[mcPhoton] ->Fill(ecluster,phicluster);
	  fhEtaMC[mcPhoton] ->Fill(ecluster,etacluster);
	  
	  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
	    {
	      fhEMC  [mcConversion] ->Fill(ecluster);
	      fhPtMC [mcConversion] ->Fill(ptcluster);
	      fhPhiMC[mcConversion] ->Fill(ecluster,phicluster);
	      fhEtaMC[mcConversion] ->Fill(ecluster,etacluster);
	      
	      if(fCheckConversion){
		if(ph->IsTagged()) fhPtConversionTagged ->Fill(ptcluster);
		if(ptcluster > 0.5)fhEtaPhiConversion   ->Fill(etacluster,phicluster);
		else               fhEtaPhi05Conversion ->Fill(etacluster,phicluster);
	      }
	    }			
	  
	  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)){
	    fhEMC  [mcPrompt] ->Fill(ecluster);
	    fhPtMC [mcPrompt] ->Fill(ptcluster);
	    fhPhiMC[mcPrompt] ->Fill(ecluster,phicluster);
	    fhEtaMC[mcPrompt] ->Fill(ecluster,etacluster);          
	  }
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
	    {
	      fhEMC  [mcFragmentation] ->Fill(ecluster);
	      fhPtMC [mcFragmentation] ->Fill(ptcluster);
	      fhPhiMC[mcFragmentation] ->Fill(ecluster,phicluster);
	      fhEtaMC[mcFragmentation] ->Fill(ecluster,etacluster);
	      
	    }
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
	    {
	      fhEMC  [mcISR] ->Fill(ecluster);
	      fhPtMC [mcISR] ->Fill(ptcluster);
	      fhPhiMC[mcISR] ->Fill(ecluster,phicluster);
	      fhEtaMC[mcISR] ->Fill(ecluster,etacluster);          
	    }
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
	    {
	      fhEMC  [mcPi0Decay] ->Fill(ecluster);
	      fhPtMC [mcPi0Decay] ->Fill(ptcluster);
	      fhPhiMC[mcPi0Decay] ->Fill(ecluster,phicluster);
	      fhEtaMC[mcPi0Decay] ->Fill(ecluster,etacluster);
	    }
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
	    {
	      fhEMC  [mcOtherDecay] ->Fill(ecluster);
	      fhPtMC [mcOtherDecay] ->Fill(ptcluster);
	      fhPhiMC[mcOtherDecay] ->Fill(ecluster,phicluster);
	      fhEtaMC[mcOtherDecay] ->Fill(ecluster,etacluster);
	    }
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron))
	{
	  fhEMC  [mcAntiNeutron] ->Fill(ecluster);
	  fhPtMC [mcAntiNeutron] ->Fill(ptcluster);
	  fhPhiMC[mcAntiNeutron] ->Fill(ecluster,phicluster);
	  fhEtaMC[mcAntiNeutron] ->Fill(ecluster,etacluster);
	  if(ph->IsTagged() && fCheckConversion) fhPtAntiNeutronTagged ->Fill(ptcluster);
	  
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))
	{
	  fhEMC  [mcAntiProton] ->Fill(ecluster);
	  fhPtMC [mcAntiProton] ->Fill(ptcluster);
	  fhPhiMC[mcAntiProton] ->Fill(ecluster,phicluster);
	  fhEtaMC[mcAntiProton] ->Fill(ecluster,etacluster);
	  if(ph->IsTagged() && fCheckConversion) fhPtAntiProtonTagged ->Fill(ptcluster);
	} 
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))
	{
	  fhEMC  [mcElectron] ->Fill(ecluster);
	  fhPtMC [mcElectron] ->Fill(ptcluster);
	  fhPhiMC[mcElectron] ->Fill(ecluster,phicluster);
	  fhEtaMC[mcElectron] ->Fill(ecluster,etacluster);
	}     
      else{
        fhEMC  [mcOther] ->Fill(ecluster);
        fhPtMC [mcOther] ->Fill(ptcluster);
        fhPhiMC[mcOther] ->Fill(ecluster,phicluster);
        fhEtaMC[mcOther] ->Fill(ecluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtUnknownTagged ->Fill(ptcluster);
        
	
        //		 printf(" AliAnaPhoton::MakeAnalysisFillHistograms() - Label %d, pT %2.3f Unknown, bits set: ",
        //					ph->GetLabel(),ph->Pt());
        //		  for(Int_t i = 0; i < 20; i++) {
        //			  if(GetMCAnalysisUtils()->CheckTagBit(tag,i)) printf(" %d, ",i);
        //		  }
        //		  printf("\n");
        
      }
      
      //....................................................................
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      if(label < 0) {
	printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
	continue;
      }
      
      Float_t eprim   = 0;
      Float_t ptprim  = 0;
      if(GetReader()->ReadStack()){
	
	if(label >=  stack->GetNtrack()) {
          if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
          continue ;
	}
	
	primary = stack->Particle(label);
	if(!primary){
	  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
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
	    if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                       label, mcparticles->GetEntriesFast());
	    continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles->At(label);
          
       }
	
       if(!aodprimary){
	 printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
	 continue;
       }
       
       eprim   = aodprimary->E();
       ptprim  = aodprimary->Pt();
       
      }
      
      fh2E     ->Fill(ecluster, eprim);
      fh2Pt    ->Fill(ptcluster, ptprim);     
      fhDeltaE ->Fill(eprim-ecluster);
      fhDeltaPt->Fill(ptprim-ptcluster);     
      if(eprim > 0)  fhRatioE  ->Fill(ecluster/eprim);
      if(ptprim > 0) fhRatioPt ->Fill(ptcluster/ptprim); 		
      
    }//Histograms with MC
    
  }// aod loop
  
}


//__________________________________________________________________
void AliAnaPhoton::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
  printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
  printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  printf("Reject clusters with a track matched = %d\n",fRejectTrackMatch);
  printf("Check Pair Conversion                = %d\n",fCheckConversion);
  printf("Add conversion pair to AOD           = %d\n",fAddConvertedPairsToAOD);
  printf("Conversion pair mass cut             = %f\n",fMassCut);
  printf("Conversion selection cut : A < %1.2f; %1.3f < Dphi < %1.3f; Deta < %1.3f\n",
         fConvAsymCut,fConvDPhiMinCut, fConvDPhiMaxCut, fConvDEtaCut);
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
} 
