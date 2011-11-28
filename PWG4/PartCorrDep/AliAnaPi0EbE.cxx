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
/* $Id: AliAnaPi0EbE.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
// Class for the analysis of high pT pi0 event by event
// Pi0 identified by one of the following:
//  -Invariant mass of 2 cluster in calorimeter
//  -Shower shape analysis in calorimeter
//  -Invariant mass of one cluster in calorimeter and one photon reconstructed in CTS
//
// -- Author: Gustavo Conesa (LNF-INFN) &  Raphaelle Ichou (SUBATECH)
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TList.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TH3F.h>
//#include "Riostream.h"

// --- Analysis system --- 
#include "AliAnaPi0EbE.h" 
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnaPi0EbE)
  
//____________________________
AliAnaPi0EbE::AliAnaPi0EbE() : 
    AliAnaPartCorrBaseClass(),    fAnaType(kIMCalo),            fCalorimeter(""),
    fMinDist(0.),fMinDist2(0.),   fMinDist3(0.),	              fFillWeightHistograms(kFALSE),
    fInputAODGammaConvName(""),
    //Histograms
    fhPtPi0(0),                   fhEPi0(0),                    fhEEtaPhiPi0(0),
    //Shower shape histos
    fhEDispersion(0),             fhELambda0(0),                fhELambda1(0), 
    fhELambda0NoTRD(0),           fhELambda0FracMaxCellCut(0),  
    fhEFracMaxCell(0),            fhEFracMaxCellNoTRD(0),            
    fhENCells(0),                 fhETime(0),                   fhEPairDiffTime(0),    
    //MC histos
    fhPtMCNoPi0(0),               fhPhiMCNoPi0(0),              fhEtaMCNoPi0(0), 
    fhPtMCPi0(0),                 fhPhiMCPi0(0),                fhEtaMCPi0(0),
    // Weight studies
    fhECellClusterRatio(0),       fhECellClusterLogRatio(0),                 
    fhEMaxCellClusterRatio(0),    fhEMaxCellClusterLogRatio(0)
{
  //default ctor
  
  for(Int_t i = 0; i < 6; i++){
    fhEMCLambda0[i]     = 0;
    fhEMCLambda0NoTRD[i]= 0;
    fhEMCLambda0FracMaxCellCut[i]= 0;
    fhEMCFracMaxCell[i] = 0;
    fhEMCLambda1[i]     = 0;
    fhEMCDispersion[i]  = 0;
  }
  
  //Weight studies
  for(Int_t i =0; i < 14; i++){
    fhLambda0ForW0[i] = 0;
    //fhLambda1ForW0[i] = 0;
  }
  
  //Initialize parameters
  InitParameters();
  
}

//_____________________________________________________________________________________
void AliAnaPi0EbE::FillSelectedClusterHistograms(AliVCluster* cluster, const Int_t tag){
  
  // Fill shower shape, timing and other histograms for selected clusters from decay
  
  Float_t e    = cluster->E();
  Float_t disp = cluster->GetDispersion()*cluster->GetDispersion();
  Float_t l0   = cluster->GetM02();
  Float_t l1   = cluster->GetM20(); 
  Int_t   nSM  = GetModuleNumber(cluster);
  
  AliVCaloCells * cell = 0x0; 
  if(fCalorimeter == "PHOS") 
    cell = GetPHOSCells();
  else		              
    cell = GetEMCALCells();
  
  Float_t maxCellFraction = 0;
  GetCaloUtils()->GetMaxEnergyCell(cell, cluster, maxCellFraction);
  fhEFracMaxCell->Fill(e,maxCellFraction);  
  
  FillWeightHistograms(cluster);
  
  fhEDispersion->Fill(e, disp);   
  fhELambda0   ->Fill(e, l0  );  
  fhELambda1   ->Fill(e, l1  );  
  
  if(fCalorimeter=="EMCAL" && nSM < 6) {
    fhELambda0NoTRD->Fill(e, l0  );
    fhEFracMaxCellNoTRD->Fill(e,maxCellFraction);  
  }
  
  if(maxCellFraction < 0.5) 
    fhELambda0FracMaxCellCut->Fill(e, l0  );  
  
  fhETime  ->Fill(e, cluster->GetTOF()*1.e9);
  fhENCells->Fill(e, cluster->GetNCells());
  
  if(IsDataMC()) {
    //Photon1
    if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  ){
      fhEMCLambda0[mcPi0]    ->Fill(e, l0);
      fhEMCLambda1[mcPi0]    ->Fill(e, l1);
      fhEMCDispersion[mcPi0] ->Fill(e, disp);
      
      fhEMCFracMaxCell[mcPi0]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcPi0]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcPi0]->Fill(e, l0  );  
      
    }//pi0
    else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  ){
      fhEMCLambda0[mcEta]    ->Fill(e, l0);
      fhEMCLambda1[mcEta]    ->Fill(e, l1);
      fhEMCDispersion[mcEta] ->Fill(e, disp);
      fhEMCFracMaxCell[mcEta]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcEta]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcEta]->Fill(e, l0  );  
    }//eta          
    else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) &&
              GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) ){
      fhEMCLambda0[mcConversion]    ->Fill(e, l0);
      fhEMCLambda1[mcConversion]    ->Fill(e, l1);
      fhEMCDispersion[mcConversion] ->Fill(e, disp);
      fhEMCFracMaxCell[mcConversion]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcConversion]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcConversion]->Fill(e, l0  );  
    }//conversion photon
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) ){
      fhEMCLambda0[mcPhoton]    ->Fill(e, l0);
      fhEMCLambda1[mcPhoton]    ->Fill(e, l1);
      fhEMCDispersion[mcPhoton] ->Fill(e, disp);
      fhEMCFracMaxCell[mcPhoton]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcPhoton]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcPhoton]->Fill(e, l0  );  
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)){
      fhEMCLambda0[mcElectron]    ->Fill(e, l0);
      fhEMCLambda1[mcElectron]    ->Fill(e, l1);
      fhEMCDispersion[mcElectron] ->Fill(e, disp);
      fhEMCFracMaxCell[mcElectron]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcElectron]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcElectron]->Fill(e, l0  );  
    }//electron
    else {
      fhEMCLambda0[mcHadron]    ->Fill(e, l0);
      fhEMCLambda1[mcHadron]    ->Fill(e, l1);
      fhEMCDispersion[mcHadron] ->Fill(e, disp);
      fhEMCFracMaxCell[mcHadron]->Fill(e,maxCellFraction);  
      if(fCalorimeter=="EMCAL" && nSM < 6) 
        fhEMCLambda0NoTRD[mcHadron]->Fill(e, l0  );
      if(maxCellFraction < 0.5) 
        fhEMCLambda0FracMaxCellCut[mcHadron]->Fill(e, l0  );  
    }//other particles 
  }//MC
}

//________________________________________________________
void AliAnaPi0EbE::FillWeightHistograms(AliVCluster *clus)
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
    RecalibrateCellAmplitude(amp,id);
    
    energy    += amp;
    
    if(amp> ampMax) 
      ampMax = amp;
    
  } // energy loop       
  
  if(energy <=0 ) {
    printf("AliAnaPi0EbE::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  fhEMaxCellClusterRatio   ->Fill(energy,ampMax/energy);
  fhEMaxCellClusterLogRatio->Fill(energy,TMath::Log(ampMax/energy));
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    RecalibrateCellAmplitude(amp,id);
    
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
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
  }// EMCAL
}

//___________________________________________
TObjString *  AliAnaPi0EbE::GetAnalysisCuts()
{	
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPi0EbE ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fAnaType=%d (Pi0 selection type) \n",fAnaType) ;
  parList+=onePar ;
  
  if(fAnaType == kSSCalo){
    snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
    parList+=onePar ;
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  if(fAnaType == kSSCalo) parList += GetCaloPID()->GetPIDParametersList() ;
  
  return new TObjString(parList) ;
}

//_____________________________________________
TList *  AliAnaPi0EbE::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("Pi0EbEHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins();           Float_t ptmax  = GetHistoPtMax();           Float_t ptmin  = GetHistoPtMin();
  Int_t nphibins = GetHistoPhiBins();          Float_t phimax = GetHistoPhiMax();          Float_t phimin = GetHistoPhiMin();
  Int_t netabins = GetHistoEtaBins();          Float_t etamax = GetHistoEtaMax();          Float_t etamin = GetHistoEtaMin();
  Int_t ssbins   = GetHistoShowerShapeBins();  Float_t ssmax  = GetHistoShowerShapeMax();  Float_t ssmin  = GetHistoShowerShapeMin();
  Int_t tdbins   = GetHistoDiffTimeBins() ;    Float_t tdmax  = GetHistoDiffTimeMax();     Float_t tdmin  = GetHistoDiffTimeMin();
  Int_t tbins    = GetHistoTimeBins() ;        Float_t tmax   = GetHistoTimeMax();         Float_t tmin   = GetHistoTimeMin();
  Int_t nbins    = GetHistoNClusterCellBins(); Int_t   nmax   = GetHistoNClusterCellMax(); Int_t   nmin   = GetHistoNClusterCellMin(); 

  fhPtPi0  = new TH1F("hPtPi0","Number of identified  #pi^{0} decay",nptbins,ptmin,ptmax); 
  fhPtPi0->SetYTitle("N");
  fhPtPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
  outputContainer->Add(fhPtPi0) ; 
  
  fhEPi0  = new TH1F("hEPi0","Number of identified  #pi^{0} decay",nptbins,ptmin,ptmax); 
  fhEPi0->SetYTitle("N");
  fhEPi0->SetXTitle("E  #pi^{0}(GeV)");
  outputContainer->Add(fhEPi0) ; 
  
  fhEEtaPhiPi0  = new TH3F
  ("hEEtaPhiPi0","Selected #pi^{0} pairs: E vs #eta vs #phi",nptbins,ptmin,ptmax,netabins,etamin,etamax, nphibins,phimin,phimax); 
  fhEEtaPhiPi0->SetZTitle("#phi");
  fhEEtaPhiPi0->SetYTitle("#eta");
  fhEEtaPhiPi0->SetXTitle("E (GeV)");
  outputContainer->Add(fhEEtaPhiPi0) ; 
  
  ////////
  
  if(fAnaType == kIMCalo || fAnaType == kIMCaloTracks ){
    
    fhEDispersion  = new TH2F
    ("hEDispersion","Selected #pi^{0} pairs: E vs dispersion",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhEDispersion->SetYTitle("D^{2}");
    fhEDispersion->SetXTitle("E (GeV)");
    outputContainer->Add(fhEDispersion) ; 
    
    fhELambda0  = new TH2F
    ("hELambda0","Selected #pi^{0} pairs: E vs #lambda_{0}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhELambda0->SetYTitle("#lambda_{0}^{2}");
    fhELambda0->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda0) ; 

    fhELambda1  = new TH2F
    ("hELambda1","Selected #pi^{0} pairs: E vs #lambda_{1}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhELambda1->SetYTitle("#lambda_{1}^{2}");
    fhELambda1->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda1) ; 
        
    fhELambda0FracMaxCellCut  = new TH2F
    ("hELambda0FracMaxCellCut","Selected #pi^{0} pairs: E vs #lambda_{0}, Max cell fraction of energy < 0.5",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhELambda0FracMaxCellCut->SetYTitle("#lambda_{0}^{2}");
    fhELambda0FracMaxCellCut->SetXTitle("E (GeV)");
    outputContainer->Add(fhELambda0FracMaxCellCut) ; 

    fhEFracMaxCell  = new TH2F
    ("hEFracMaxCell","Selected #pi^{0} pairs: E vs #lambda_{0}, Max cell fraction of energy",nptbins,ptmin,ptmax,100,0,1); 
    fhEFracMaxCell->SetYTitle("Fraction");
    fhEFracMaxCell->SetXTitle("E (GeV)");
    outputContainer->Add(fhEFracMaxCell) ; 

    if(fCalorimeter=="EMCAL"){
      fhELambda0NoTRD  = new TH2F
      ("hELambda0NoTRD","Selected #pi^{0} pairs: E vs #lambda_{0}, not behind TRD",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhELambda0NoTRD->SetYTitle("#lambda_{0}^{2}");
      fhELambda0NoTRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhELambda0NoTRD) ; 
      
      fhEFracMaxCellNoTRD  = new TH2F
      ("hEFracMaxCellNoTRD","Selected #pi^{0} pairs: E vs #lambda_{0}, Max cell fraction of energy, not behind TRD",nptbins,ptmin,ptmax,100,0,1); 
      fhEFracMaxCellNoTRD->SetYTitle("Fraction");
      fhEFracMaxCellNoTRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhEFracMaxCellNoTRD) ; 
    }
    
    fhENCells  = new TH2F ("hENCells","N cells in cluster vs E ", nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhENCells->SetXTitle("E (GeV)");
    fhENCells->SetYTitle("# of cells in cluster");
    outputContainer->Add(fhENCells);  
    
    fhETime = new TH2F("hETime","cluster time vs pair E",nptbins,ptmin,ptmax, tbins,tmin,tmax);
    fhETime->SetXTitle("E (GeV)");
    fhETime->SetYTitle(" t (ns)");
    outputContainer->Add(fhETime);    
    
    fhEPairDiffTime = new TH2F("hEPairDiffTime","cluster pair time difference vs E",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhEPairDiffTime->SetXTitle("E_{pair} (GeV)");
    fhEPairDiffTime->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhEPairDiffTime);
    
    
  }// Invariant mass analysis in calorimeters only
  
  if(fFillWeightHistograms){
    
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy, for selected decay photons from neutral meson",
                                     nptbins,ptmin,ptmax, 100,0,1.); 
    fhECellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("E_{cell i}/E_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,-10,0); 
    fhECellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy, for selected decay photons from neutral meson",
                                        nptbins,ptmin,ptmax, 100,0,1.); 
    fhEMaxCellClusterRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("E_{max cell}/E_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy, for selected decay photons from neutral meson",
                                           nptbins,ptmin,ptmax, 100,-10,0); 
    fhEMaxCellClusterLogRatio->SetXTitle("E_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (E_{max cell}/E_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    for(Int_t iw = 0; iw < 14; iw++){
      fhLambda0ForW0[iw]  = new TH2F (Form("hLambda0ForW0%d",iw),Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for selected decay photons from neutral meson",1+0.5*iw),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLambda0ForW0[iw]->SetXTitle("E_{cluster}");
      fhLambda0ForW0[iw]->SetYTitle("#lambda^{2}_{0}");
      outputContainer->Add(fhLambda0ForW0[iw]); 
      
//      fhLambda1ForW0[iw]  = new TH2F (Form("hLambda1ForW0%d",iw),Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for selected decay photons from neutral meson",0.5+0.5*iw),
//                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
//      fhLambda1ForW0[iw]->SetXTitle("E_{cluster}");
//      fhLambda1ForW0[iw]->SetYTitle("#lambda^{2}_{1}");
//      outputContainer->Add(fhLambda1ForW0[iw]); 
      
    }
  }  
  
  if(IsDataMC()) {
    if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
       GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      
      fhPtMCPi0  = new TH1F("hPtMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax); 
      fhPtMCPi0->SetYTitle("N");
      fhPtMCPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
      outputContainer->Add(fhPtMCPi0) ; 
      
      fhPhiMCPi0  = new TH2F
      ("hPhiMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiMCPi0->SetYTitle("#phi");
      fhPhiMCPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhPhiMCPi0) ; 
      
      fhEtaMCPi0  = new TH2F
      ("hEtaMCPi0","Identified pi0 from pi0",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaMCPi0->SetYTitle("#eta");
      fhEtaMCPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhEtaMCPi0) ;
      
      fhPtMCNoPi0  = new TH1F("hPtMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax); 
      fhPtMCNoPi0->SetYTitle("N");
      fhPtMCNoPi0->SetXTitle("p_{T #pi^{0}}(GeV/c)");
      outputContainer->Add(fhPtMCNoPi0) ; 
      
      fhPhiMCNoPi0  = new TH2F
      ("hPhiMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiMCNoPi0->SetYTitle("#phi");
      fhPhiMCNoPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhPhiMCNoPi0) ; 
      
      fhEtaMCNoPi0  = new TH2F
      ("hEtaMCNoPi0","Identified pi0 not from pi0",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaMCNoPi0->SetYTitle("#eta");
      fhEtaMCNoPi0->SetXTitle("p_{T #pi^{0}} (GeV/c)");
      outputContainer->Add(fhEtaMCNoPi0) ;
      
      if(fAnaType == kIMCalo){
        TString ptype[] ={"#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron"}; 
        TString pname[] ={"Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron"};
        for(Int_t i = 0; i < 6; i++){ 
          
          fhEMCLambda0[i]  = new TH2F(Form("hELambda0_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCLambda0[i]->SetYTitle("#lambda_{0}^{2}");
          fhEMCLambda0[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda0[i]) ; 
          
          if(fCalorimeter=="EMCAL"){
            fhEMCLambda0NoTRD[i]  = new TH2F(Form("hELambda0NoTRD_MC%s",pname[i].Data()),
                                             Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}, NoTRD",ptype[i].Data()),
                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhEMCLambda0NoTRD[i]->SetYTitle("#lambda_{0}^{2}");
            fhEMCLambda0NoTRD[i]->SetXTitle("E (GeV)");
            outputContainer->Add(fhEMCLambda0NoTRD[i]) ; 
          }
          
          fhEMCLambda0FracMaxCellCut[i]  = new TH2F(Form("hELambda0FracMaxCellCut_MC%s",pname[i].Data()),
                                                    Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}, Max cell fraction of energy < 0.5 ",ptype[i].Data()),
                                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCLambda0FracMaxCellCut[i]->SetYTitle("#lambda_{0}^{2}");
          fhEMCLambda0FracMaxCellCut[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda0FracMaxCellCut[i]) ; 
          
          fhEMCFracMaxCell[i]  = new TH2F(Form("hEFracMaxCell_MC%s",pname[i].Data()),
                                          Form("Selected pair, cluster from %s : E vs Max cell fraction of energy",ptype[i].Data()),
                                          nptbins,ptmin,ptmax,100,0,1); 
          fhEMCFracMaxCell[i]->SetYTitle("Fraction");
          fhEMCFracMaxCell[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCFracMaxCell[i]) ;           

          fhEMCLambda1[i]  = new TH2F(Form("hELambda1_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{1}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCLambda1[i]->SetYTitle("#lambda_{1}^{2}");
          fhEMCLambda1[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda1[i]) ; 
                    
          fhEMCDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pname[i].Data()),
                                         Form("Selected pair, cluster from %s : E vs dispersion^{2}",ptype[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCDispersion[i]->SetYTitle("D^{2}");
          fhEMCDispersion[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCDispersion[i]) ; 
      
        }//
        
      }//kIMCalo
    } //Not MC reader
  }//Histos with MC
  
  
  //Keep neutral meson selection histograms if requiered
  //Setting done in AliNeutralMesonSelection
  
  if(fAnaType!=kSSCalo && GetNeutralMesonSelection()){
    
    TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
    delete nmsHistos;
	  
  }
  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaPi0EbE::Init()
{ 
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPi0EbE::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//____________________________________________________________________________
void AliAnaPi0EbE::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0EbE_");
  
  fInputAODGammaConvName = "PhotonsCTS" ;
  fAnaType = kIMCalo ;
  fCalorimeter = "EMCAL" ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  
  switch(fAnaType) 
  {
    case kIMCalo:
      MakeInvMassInCalorimeter();
      break;
      
    case kSSCalo:
      MakeShowerShapeIdentification();
      break;
      
    case kIMCaloTracks:
      MakeInvMassInCalorimeterAndCTS();
      break;
      
  }
}

//____________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeter() 
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  Int_t tag1 = 0;
  Int_t tag2 = 0;
  Int_t tag  = 0;
  
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - No input calo photons in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  //Get shower shape information of clusters
  TObjArray *clusters = 0;
  if     (fCalorimeter=="EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter=="PHOS")  clusters = GetPHOSClusters() ;
  
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast()-1; iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    
    //Vertex cut in case of mixed events
    Int_t evtIndex1 = 0 ; 
    if(GetMixedEvent())
      evtIndex1 = GetMixedEvent()->EventIndexForCaloCluster(photon1->GetCaloLabel(0)) ;
    if(TMath::Abs(GetVertex(evtIndex1)[2]) > GetZvertexCut()) continue ;  //vertex cut
    mom1 = *(photon1->Momentum());
    
    //Get original cluster, to recover some information
    Int_t iclus = -1;
    AliVCluster *cluster1 = FindCluster(clusters,photon1->GetCaloLabel(0),iclus); 
    
    if(!cluster1){
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - First cluster not found\n");
      return;
    }
    
    for(Int_t jphoton = iphoton+1; jphoton < GetInputAODBranch()->GetEntriesFast(); jphoton++){
      
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(jphoton));
      Int_t evtIndex2 = 0 ; 
      if(GetMixedEvent())
        evtIndex2 = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
      if(GetMixedEvent() && (evtIndex1 == evtIndex2))
        continue ; 
      if(TMath::Abs(GetVertex(evtIndex2)[2]) > GetZvertexCut()) continue ;  //vertex cut
      mom2 = *(photon2->Momentum());
      
      //Get original cluster, to recover some information
      Int_t iclus2;
      AliVCluster *cluster2 = FindCluster(clusters,photon2->GetCaloLabel(0),iclus2,iclus+1); 
      
      if(!cluster2){
        printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Second cluster not found\n");
        return;
      }
      
      Float_t e1    = photon1->E();      
      Float_t e2    = photon2->E();
      
      //Select clusters with good time window difference
      Float_t tof1  = cluster1->GetTOF()*1e9;;
      Float_t tof2  = cluster2->GetTOF()*1e9;;
      Double_t t12diff = tof1-tof2;
      fhEPairDiffTime->Fill(e1+e2,    t12diff);
      if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2,fCalorimeter))
      {
        if(GetDebug()>1) 
          printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Selected gamma pair: pt %f, phi %f, eta%f \n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
        
        //Play with the MC stack if available
        if(IsDataMC()){
          //Check origin of the candidates
          Int_t  label1 = photon1->GetLabel();
          Int_t  label2 = photon2->GetLabel();
          if(label1>=0)tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader(), photon1->GetInputFileIndex());
          if(label2>=0)tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(), photon2->GetInputFileIndex());
          
          if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
          if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) && 
             GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)){
            
            //Check if pi0 mother is the same
            if(GetReader()->ReadStack()){ 
              if(label1>=0){
                TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
                label1 = mother1->GetFirstMother();
                //mother1 = GetMCStack()->Particle(label1);//pi0
              }
              if(label2>=0){
                TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
                label2 = mother2->GetFirstMother();
                //mother2 = GetMCStack()->Particle(label2);//pi0
              }
            }
            else if(GetReader()->ReadAODMCParticles()){//&& (input > -1)){
              if(label1>=0){
                AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon1->GetInputFileIndex()))->At(label1);//photon in kine tree
                label1 = mother1->GetMother();
                //mother1 = GetMCStack()->Particle(label1);//pi0
              }
              if(label2>=0){
                AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon2->GetInputFileIndex()))->At(label2);//photon in kine tree
                label2 = mother2->GetMother();
                //mother2 = GetMCStack()->Particle(label2);//pi0
              }
            }
            
            //printf("mother1 %d, mother2 %d\n",label1,label2);
            if(label1 == label2 && label1>=0)
              GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
          }
        }//Work with stack also   
        
        
        //Fill some histograms about shower shape
        if(clusters && GetReader()->GetDataType()!=AliCaloTrackReader::kMC){
          FillSelectedClusterHistograms(cluster1, tag1);
          FillSelectedClusterHistograms(cluster2, tag2);
        }
        
        // Tag both photons as decay
        photon1->SetTagged(kTRUE);
        photon2->SetTagged(kTRUE);

        //Create AOD for analysis
        mom = mom1+mom2;
        AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
        //pi0.SetLabel(calo->GetLabel());
        pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
        pi0.SetDetector(photon1->GetDetector());
        pi0.SetTag(tag);  
        //Set the indeces of the original caloclusters  
        pi0.SetCaloLabel(photon1->GetCaloLabel(0), photon2->GetCaloLabel(0));
        //pi0.SetInputFileIndex(input);
        AddAODParticle(pi0);
      }//pi0
      
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeter() - End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() 
{
  //Do analysis and fill aods
  //Search for the photon decay in calorimeters
  //Read photon list from AOD, produced in class AliAnaPhoton and AliGammaConversion
  //Check if 2 photons have the mass of the pi0.
  
  TLorentzVector mom1;
  TLorentzVector mom2;
  TLorentzVector mom ;
  Int_t tag1 = 0;
  Int_t tag2 = 0;
  Int_t tag  = 0;
  Int_t evtIndex = 0;
  
  // Check calorimeter input
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input calo photons in AOD branch with name < %s > , STOP\n",GetInputAODName().Data());
    abort();
  }
  
  // Get the array with conversion photons
  TClonesArray * inputAODGammaConv = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaConvName);
  if(!inputAODGammaConv) {
    
    inputAODGammaConv = (TClonesArray *) GetReader()->GetInputEvent()->FindListObject(fInputAODGammaConvName);
    
    if(!inputAODGammaConv) {
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input gamma conversions in AOD branch with name < %s >\n",fInputAODGammaConvName.Data());
      
      return;
    }
  }  
  
  //Get shower shape information of clusters
  TObjArray *clusters = 0;
  if     (fCalorimeter=="EMCAL") clusters = GetEMCALClusters();
  else if(fCalorimeter=="PHOS")  clusters = GetPHOSClusters() ;  
  
  Int_t nCTS  = inputAODGammaConv->GetEntriesFast();
  Int_t nCalo = GetInputAODBranch()->GetEntriesFast();
  if(nCTS<=0 || nCalo <=0) {
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - nCalo %d, nCTS %d, cannot loop\n",nCalo,nCTS);
    return;
  }
  
  if(GetDebug() > 1)
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Number of conversion photons %d\n",nCTS);
  
  // Do the loop, first calo, second CTS
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    mom1 = *(photon1->Momentum());
    
    //Get original cluster, to recover some information
    Int_t iclus = -1;
    AliVCluster *cluster = FindCluster(clusters,photon1->GetCaloLabel(0),iclus);     
    
    for(Int_t jphoton = 0; jphoton < nCTS; jphoton++){
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (inputAODGammaConv->At(jphoton));
      if(GetMixedEvent())
        evtIndex = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
      if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
      
      mom2 = *(photon2->Momentum());
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2,fCalorimeter)){
        if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Selected gamma pair: pt %f, phi %f, eta%f\n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
        
        if(IsDataMC()){
          Int_t	label1 = photon1->GetLabel();
          Int_t	label2 = photon2->GetLabel();
          if(label1>=0)tag1 = GetMCAnalysisUtils()->CheckOrigin(label1, GetReader(), photon1->GetInputFileIndex());
          if(label2>=0)tag2 = GetMCAnalysisUtils()->CheckOrigin(label2, GetReader(), photon2->GetInputFileIndex());
          if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - Origin of: photon1 %d; photon2 %d \n",tag1, tag2);
          if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0Decay) && 
             GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0Decay)){
            //Check if pi0 mother is the same
            
            if(GetReader()->ReadStack()){ 
              if(label1>=0){
                TParticle * mother1 = GetMCStack()->Particle(label1);//photon in kine tree
                label1 = mother1->GetFirstMother();
                //mother1 = GetMCStack()->Particle(label1);//pi0
              }
              if(label2>=0){
                TParticle * mother2 = GetMCStack()->Particle(label2);//photon in kine tree
                label2 = mother2->GetFirstMother();
                //mother2 = GetMCStack()->Particle(label2);//pi0
              }
            }
            else if(GetReader()->ReadAODMCParticles()&& label1>=0 && label2>=0){ //&& (input > -1)){
              if(label1>=0){
                AliAODMCParticle * mother1 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon1->GetInputFileIndex()))->At(label1);//photon in kine tree
                label1 = mother1->GetMother();
                //mother1 = GetMCStack()->Particle(label1);//pi0
              }
              if(label2>=0){
                AliAODMCParticle * mother2 = (AliAODMCParticle *) (GetReader()->GetAODMCParticles(photon2->GetInputFileIndex()))->At(label2);//photon in kine tree
                label2 = mother2->GetMother();
                //mother2 = GetMCStack()->Particle(label2);//pi0
              }
            }
            
            //printf("mother1 %d, mother2 %d\n",label1,label2);
            if(label1 == label2 && label1>=0)
              GetMCAnalysisUtils()->SetTagBit(tag,AliMCAnalysisUtils::kMCPi0);
          }
        }//Work with stack also   
        
        //Fill some histograms about shower shape
        if(cluster && GetReader()->GetDataType()!=AliCaloTrackReader::kMC){
          FillSelectedClusterHistograms(cluster, tag1);
        }        
        
        // Tag both photons as decay
        photon1->SetTagged(kTRUE);
        photon2->SetTagged(kTRUE);        
        
        //Create AOD for analysis
        mom = mom1+mom2;
        AliAODPWG4Particle pi0 = AliAODPWG4Particle(mom);
        //pi0.SetLabel(calo->GetLabel());
        pi0.SetIdentifiedParticleType(AliCaloPID::kPi0);
        pi0.SetDetector(photon1->GetDetector());
        pi0.SetTag(tag);
        //Set the indeces of the original tracks or caloclusters  
        pi0.SetCaloLabel(photon1->GetCaloLabel(0), -1);
        pi0.SetTrackLabel(photon2->GetTrackLabel(0), photon2->GetTrackLabel(1));
        //pi0.SetInputFileIndex(input);
        AddAODParticle(pi0);
      }//pi0
    }//2n photon loop
    
  }//1st photon loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - End fill AODs \n");  
  
}


//__________________________________________________________________
void  AliAnaPi0EbE::MakeShowerShapeIdentification() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl = 0x0; 
  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  if(!pl) {
    Info("MakeShowerShapeIdentification","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }  
	
  TLorentzVector mom ;
  for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++){
    AliVCluster * calo = (AliVCluster*) (pl->At(icalo));	
    
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
    }
    if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
    
    //Get Momentum vector, 
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;}//Assume that come from vertex in straight line
    else{
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
	  
    //If too small or big pt, skip it
    if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 
    //Check acceptance selection
    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue ;
    }
    
    //Create AOD for analysis
    AliAODPWG4Particle aodpi0 = AliAODPWG4Particle(mom);
    aodpi0.SetLabel(calo->GetLabel());
    //Set the indeces of the original caloclusters  
    aodpi0.SetCaloLabel(calo->GetID(),-1);
    aodpi0.SetDetector(fCalorimeter);
    if(GetDebug() > 1) 
      printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",aodpi0.Pt(),aodpi0.Phi(),aodpi0.Eta());	
    
    //Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist) //In bad channel (PHOS cristal size 2.2x2.2 cm)
      continue ;
    
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: Bad channel cut passed %4.2f\n",distBad);
    
    if(distBad > fMinDist3) aodpi0.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodpi0.SetDistToBad(1) ; 
    else aodpi0.SetDistToBad(0) ;
    
    //Check PID
    //PID selection or bit setting
    if(IsCaloPIDOn()){
      //Skip matched clusters with tracks
      if(IsTrackMatched(calo, GetReader()->GetInputEvent())) continue ;
      
      // Get most probable PID, 2 options check bayesian PID weights or redo PID
      // By default, redo PID
     
      aodpi0.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PDG of identified particle %d\n",aodpi0.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not pi0, skip it.
      if(aodpi0.GetIdentifiedParticleType() != AliCaloPID::kPi0) continue ;			
      
    }
    else{
      //Set PID bits for later selection 
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodpi0, GetCaloUtils(), GetReader()->GetInputEvent());
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - Pi0 selection cuts passed: pT %3.2f, pdg %d\n",aodpi0.Pt(), aodpi0.GetIdentifiedParticleType());
    
    //Play with the MC stack if available
    //Check origin of the candidates
    if(IsDataMC()){
      if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
         GetReader()->GetDataType() != AliCaloTrackReader::kMC){
        //aodpi0.SetInputFileIndex(input);
        Int_t tag	=0;
        tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabel(),GetReader(), aodpi0.GetInputFileIndex());
        //GetMCAnalysisUtils()->CheckMultipleOrigin(calo->GetLabels(),calo->GetNLabels(), GetReader(), aodpi0.GetInputFileIndex(), tag);
        aodpi0.SetTag(tag);
        if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - Origin of candidate %d\n",aodpi0.GetTag());
      }
    }//Work with stack also   
    
    //Add AOD with pi0 object to aod branch
    AddAODParticle(aodpi0);
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - End fill AODs \n");  
  
}
//__________________________________________________________________
void  AliAnaPi0EbE::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  if(!GetOutputAODBranch()){
    printf("AliAnaPi0EbE::MakeAnalysisFillHistograms()  - No output pi0 in AOD branch with name < %s >,STOP \n",GetOutputAODName().Data());
    abort();
  }
  //Loop on stored AOD pi0
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPi0EbE::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    
    AliAODPWG4Particle* pi0 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = pi0->GetIdentifiedParticleType();
	  
    if(IsCaloPIDOn() && pdg != AliCaloPID::kPi0) continue;              
    
    //Fill pi0 histograms 
    Float_t ener  = pi0->E();
    Float_t pt    = pi0->Pt();
    Float_t phi   = pi0->Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    Float_t eta = pi0->Eta();
    
    fhPtPi0      ->Fill(pt);
    fhEPi0       ->Fill(ener);
    fhEEtaPhiPi0 ->Fill(ener,eta,phi);
    
    if(IsDataMC()){
      if((GetReader()->GetDataType() == AliCaloTrackReader::kMC && fAnaType!=kSSCalo) || 
         GetReader()->GetDataType() != AliCaloTrackReader::kMC){
        if(GetMCAnalysisUtils()->CheckTagBit(pi0->GetTag(), AliMCAnalysisUtils::kMCPi0)){
          fhPtMCPi0  ->Fill(pt);
          fhPhiMCPi0 ->Fill(pt,phi);
          fhEtaMCPi0 ->Fill(pt,eta);
        }
        else{
          fhPtMCNoPi0  ->Fill(pt);
          fhPhiMCNoPi0 ->Fill(pt,phi);
          fhEtaMCNoPi0 ->Fill(pt,eta);
        }
      }
    }//Histograms with MC
    
  }// aod loop
  
}

//__________________________________________________________________
void AliAnaPi0EbE::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print("");
  printf("Analysis Type = %d \n",  fAnaType) ;
  if(fAnaType == kSSCalo){     
    printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
    printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
    printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
    printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3); 
  } 
  printf("    \n") ;
  
} 

//___________________________________________________________________________________
void AliAnaPi0EbE::RecalibrateCellAmplitude(Float_t & amp, const Int_t id)
{
  //Recaculate cell energy if recalibration factor
  
  Int_t icol     = -1; Int_t irow     = -1; Int_t iRCU     = -1;
  Int_t nModule  = GetModuleNumberCellIndexes(id,fCalorimeter, icol, irow, iRCU);
  
  if (GetCaloUtils()->IsRecalibrationOn()) {
    if(fCalorimeter == "PHOS") {
      amp *= GetCaloUtils()->GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
    }
    else		                   {
      amp *= GetCaloUtils()->GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
    }
  }
}


