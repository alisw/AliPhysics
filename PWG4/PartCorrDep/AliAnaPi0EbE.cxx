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
  
//____________________________________________________________________________
AliAnaPi0EbE::AliAnaPi0EbE() : 
    AliAnaPartCorrBaseClass(),    fAnaType(kIMCalo),            fCalorimeter(""),
    fMinDist(0.),fMinDist2(0.),   fMinDist3(0.),	
    fInputAODGammaConv(0x0),      fInputAODGammaConvName(""),
    //Histograms
    fhPtPi0(0),                   fhEPi0(0),                    fhEEtaPhiPi0(0),
    //Shower shape histos
    fhEDispersion(0),             fhELambda0(0),                fhELambda1(0),             
    fhEdDispersion(0),            fhEdLambda0(0),               fhEdLambda1(0),            
    //Time histograms
    fhClusterPairDiffTimeE(0),    fhClusterPairDiffTimeAsy(0),
    //MC histos
    fhPtMCNoPi0(0),               fhPhiMCNoPi0(0),              fhEtaMCNoPi0(0), 
    fhPtMCPi0(0),                 fhPhiMCPi0(0),                fhEtaMCPi0(0)
{
  //default ctor
  
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
  
}

//____________________________________________________________________________
AliAnaPi0EbE::~AliAnaPi0EbE() 
{
  //dtor
  if(fInputAODGammaConv){
    fInputAODGammaConv->Clear() ; 
    delete fInputAODGammaConv ;
  }
}

//________________________________________________________________________
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

//________________________________________________________________________
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
  
  if(fAnaType == kIMCalo){
    
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
    
    
    fhEdDispersion  = new TH2F
    ("hEdDispersion","Selected #pi^{0} pairs: E vs dispersion^{2}/N_{cells}^{2}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhEdDispersion->SetYTitle("D^{2}");
    fhEdDispersion->SetXTitle("E (GeV)");
    outputContainer->Add(fhEdDispersion) ; 
    
    fhEdLambda0  = new TH2F
    ("hEdLambda0","Selected #pi^{0} pairs: E vs #lambda_{0}^{2}/N_{cells}^{2}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhEdLambda0->SetYTitle("#lambda_{0}^{2}");
    fhEdLambda0->SetXTitle("E (GeV)");
    outputContainer->Add(fhEdLambda0) ; 
    
    fhEdLambda1  = new TH2F
    ("hEdLambda1","Selected #pi^{0} pairs: E vs #lambda_{1}^{2}/N_{cells}^{2}",nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
    fhEdLambda1->SetYTitle("#lambda_{1}^{2}");
    fhEdLambda1->SetXTitle("E (GeV)");
    outputContainer->Add(fhEdLambda1) ; 
    
    
    fhClusterPairDiffTimeE = new TH2F("hClusterPairDiffTimeE","cluster pair time difference vs E",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhClusterPairDiffTimeE->SetXTitle("E_{pair} (GeV)");
    fhClusterPairDiffTimeE->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhClusterPairDiffTimeE);
    
    fhClusterPairDiffTimeAsy = new TH2F("hClusterPairDiffTime","cluster pair time difference vs pair asymmetry",100,0,1, tdbins,tdmin,tdmax);
    fhClusterPairDiffTimeAsy->SetXTitle("Asymmetry");
    fhClusterPairDiffTimeAsy->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhClusterPairDiffTimeAsy);    
    
  }// Invariant mass analysis in calorimeters only
  
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
        TString ptype[] ={"#gamma","#gamma->e^{#pm}","#pi^{0}","e^{#pm}", "hadron"}; 
        TString pname[] ={"Photon","Conversion",     "Pi0",    "Electron","Hadron"};
        for(Int_t i = 0; i < 5; i++){ 
          
          fhEMCLambda0[i]  = new TH2F(Form("hELambda0_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCLambda0[i]->SetYTitle("#lambda_{0}^{2}");
          fhEMCLambda0[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda0[i]) ; 
          
          fhEMCdLambda0[i]  = new TH2F(Form("hEdLambda0_MC%s",pname[i].Data()),
                                       Form("Selected pair, cluster from %s : E vs #lambda_{0}^{2}/N_{cells}^{2}",ptype[i].Data()),
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
          fhEMCdLambda0[i]->SetYTitle("d#lambda_{0}^{2}");
          fhEMCdLambda0[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCdLambda0[i]) ; 
          
          fhEMCLambda1[i]  = new TH2F(Form("hELambda1_MC%s",pname[i].Data()),
                                      Form("Selected pair, cluster from %s : E vs #lambda_{1}^{2}",ptype[i].Data()),
                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCLambda1[i]->SetYTitle("#lambda_{1}^{2}");
          fhEMCLambda1[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCLambda1[i]) ; 
          
          fhEMCdLambda1[i]  = new TH2F(Form("hEdLambda1_MC%s",pname[i].Data()),
                                       Form("Selected pair, cluster from %s : E vs d#lambda_{1}^{2}/N_{cells}^{2}",ptype[i].Data()),
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
          fhEMCdLambda1[i]->SetYTitle("d#lambda_{1}^{2}");
          fhEMCdLambda1[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCdLambda1[i]) ; 
          
          fhEMCDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pname[i].Data()),
                                         Form("Selected pair, cluster from %s : E vs dispersion^{2}",ptype[i].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhEMCDispersion[i]->SetYTitle("D^{2}");
          fhEMCDispersion[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCDispersion[i]) ; 
          
          fhEMCdDispersion[i]  = new TH2F(Form("hEdDispersion_MC%s",pname[i].Data()),
                                          Form("Selected pair, cluster from %s : E vs dispersion^{2}/N_{cells}^{2}",ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/50); 
          fhEMCdDispersion[i]->SetYTitle("dD^{2}");
          fhEMCdDispersion[i]->SetXTitle("E (GeV)");
          outputContainer->Add(fhEMCdDispersion[i]) ;   
          
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
  
  fInputAODGammaConvName = "gammaconv" ;
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

//__________________________________________________________________
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
  
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast()-1; iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    
    //Vertex cut in case of mixed events
    Int_t evtIndex1 = 0 ; 
    if(GetMixedEvent())
      evtIndex1 = GetMixedEvent()->EventIndexForCaloCluster(photon1->GetCaloLabel(0)) ;
    if(TMath::Abs(GetVertex(evtIndex1)[2]) > GetZvertexCut()) continue ;  //vertex cut
    mom1 = *(photon1->Momentum());
    
    //Get shower shape information of clusters
    TObjArray *clusters = 0;
    if     (fCalorimeter="EMCAL") clusters = GetEMCALClusters();
    else if(fCalorimeter="PHOS" ) clusters = GetPHOSClusters() ;
    
    Bool_t bFound1        = kFALSE;
    Int_t  caloLabel1     = photon1->GetCaloLabel(0);
    Bool_t iclus1         = -1;    
    AliVCluster *cluster1 = 0; 
    if(clusters) {
      //Get original cluster, to recover some information
      for(Int_t iclus = 0; iclus < clusters->GetEntriesFast(); iclus++){
        AliVCluster *cluster= dynamic_cast<AliVCluster*> (clusters->At(iclus));
        if(cluster){
          if     (cluster->GetID()==caloLabel1) {
            bFound1  = kTRUE  ;
            cluster1 = cluster;
            iclus1   = iclus;
          }
        }      
        if(bFound1) break;
      }// calorimeter clusters loop
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
      
      //Photon1
      Float_t e1    = photon1->E();
      Float_t tof1  = -1;
      Float_t disp1 = -1;
      Float_t l01   = -1;
      Float_t l11   = -1; 
      Int_t ncells1 = 0;
      
      //Photon2
      Float_t e2    = photon2->E();
      Float_t disp2 = -1;
      Float_t tof2  = -1;
      Float_t l02   = -1;
      Float_t l12   = -1; 
      Int_t ncells2 = 0;
      Bool_t bFound2        = kFALSE;
      Int_t  caloLabel2     = photon2->GetCaloLabel(0);
      AliVCluster *cluster2 = 0; 
      if(clusters){      
        for(Int_t iclus = iclus1+1; iclus < clusters->GetEntriesFast(); iclus++){
          AliVCluster *cluster= dynamic_cast<AliVCluster*> (clusters->At(iclus));
          if(cluster){
            if(cluster->GetID()==caloLabel2) {
              bFound2  = kTRUE  ;
              cluster2 = cluster;
            }          
          }      
          if(bFound2) break;
        }// calorimeter clusters loop
        
        //Photon/Cluster 1
        if(cluster1 && bFound1){
          disp1   = cluster1->GetDispersion()*cluster1->GetDispersion();
          l11     = cluster1->GetM20();
          l01     = cluster1->GetM02();
          tof1    = cluster1->GetTOF()*1e9;
          ncells1 = cluster1->GetNCells()*cluster1->GetNCells();
        }
        //        else printf("cluster1 not available: calo label %d / %d, cluster ID %d\n",
        //                     photon2->GetCaloLabel(0),(GetReader()->GetInputEvent())->GetNumberOfCaloClusters()-1,cluster1->GetID());
        
        //Photon/Cluster 2
        if(cluster2 && bFound2){
          disp2   = cluster2->GetDispersion()*cluster2->GetDispersion();
          l12     = cluster2->GetM20();
          l02     = cluster2->GetM02();
          tof2    = cluster2->GetTOF()*1e9;
          ncells2 = cluster2->GetNCells()*cluster2->GetNCells();
        }
        //        else printf("cluster2 not available: calo label %d / %d, cluster ID %d\n",
        //                     photon2->GetCaloLabel(0),(GetReader()->GetInputEvent())->GetNumberOfCaloClusters()-1,cluster2->GetID());    
        
        //Select clusters with good time window difference
        Double_t t12diff = tof1-tof2;
        Float_t asymmetry = TMath::Abs(e1-e2)/(e1+e2);
        fhClusterPairDiffTimeE  ->Fill(e1+e2,    t12diff);
        fhClusterPairDiffTimeAsy->Fill(asymmetry,t12diff);
        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
      }
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2))
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
          //Photon1 
          
          //printf("Signal Cl1: e %f, pt %f, disp %f, l1 %f, l0 %f, eta %f, phi %f \n",
          //       e,pt,disp,l1,l0,photon2->Eta(),photon2->Phi());
          
          fhEDispersion   ->Fill(e1, disp1);   
          fhELambda0->Fill(e1, l01  );  
          fhELambda1->Fill(e1, l11  );  
          
          fhEdDispersion   ->Fill(e1, disp1/ncells1);   
          fhEdLambda0->Fill(e1, l01/ncells1  );  
          fhEdLambda1->Fill(e1, l11/ncells1  );  
          
          //Photon2
          //printf("Signal Cl2: e %f, pt %f, disp %f, l1 %f, l0 %f, eta %f, phi %f \n",e
          //     ,pt,disp,l1,l0,photon2->Eta(),photon2->Phi());
          
          fhEDispersion   ->Fill(e2, disp2);   
          fhELambda0->Fill(e2, l02  ); 
          fhELambda1->Fill(e2, l12  ); 
          
          fhEdDispersion   ->Fill(e2, disp2/ncells2);   
          fhEdLambda0->Fill(e2, l02/ncells2  ); 
          fhEdLambda1->Fill(e2, l12/ncells2  ); 
          
          if(IsDataMC()) {
            //Photon1
            if( GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCConversion) ){
              fhEMCLambda0[0]    ->Fill(e1, l01);
              fhEMCdLambda0[0]   ->Fill(e1, l01/ncells1);
              fhEMCLambda1[0]    ->Fill(e1, l11);
              fhEMCdLambda1[0]   ->Fill(e1, l11/ncells1);
              fhEMCDispersion[0] ->Fill(e1, disp1);
              fhEMCdDispersion[0]->Fill(e1, disp1/ncells1);                
            }//photon   no conversion
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCElectron)){
              fhEMCLambda0[3]    ->Fill(e1, l01);
              fhEMCdLambda0[3]   ->Fill(e1, l01/ncells1);
              fhEMCLambda1[3]    ->Fill(e1, l11);
              fhEMCdLambda1[3]   ->Fill(e1, l11/ncells1);
              fhEMCDispersion[3] ->Fill(e1, disp1);
              fhEMCdDispersion[3]->Fill(e1, disp1/ncells1);
            }//electron
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCConversion) && 
                      GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCConversion) ){
              fhEMCLambda0[1]    ->Fill(e1, l01);
              fhEMCdLambda0[1]   ->Fill(e1, l01/ncells1);
              fhEMCLambda1[1]    ->Fill(e1, l11);
              fhEMCdLambda1[1]   ->Fill(e1, l11/ncells1);
              fhEMCDispersion[1] ->Fill(e1, disp1);
              fhEMCdDispersion[1]->Fill(e1, disp1/ncells1);           
            }//conversion photon
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCPi0)  ){
              fhEMCLambda0[2]    ->Fill(e1, l01);
              fhEMCdLambda0[2]   ->Fill(e1, l01/ncells1);
              fhEMCLambda1[2]    ->Fill(e1, l11);
              fhEMCdLambda1[2]   ->Fill(e1, l11/ncells1);
              fhEMCDispersion[2] ->Fill(e1, disp1);
              fhEMCdDispersion[2]->Fill(e1, disp1/ncells1);     
            }//pi0
            else {
              fhEMCLambda0[4]    ->Fill(e1, l01);
              fhEMCdLambda0[4]   ->Fill(e1, l01/ncells1);
              fhEMCLambda1[4]    ->Fill(e1, l11);
              fhEMCdLambda1[4]   ->Fill(e1, l11/ncells1);
              fhEMCDispersion[4] ->Fill(e1, disp1);
              fhEMCdDispersion[4]->Fill(e1, disp1/ncells1);    
            }//other particles 
            
            //Photon 2
            if( GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) ){
              fhEMCLambda0[0]    ->Fill(e2, l02);
              fhEMCdLambda0[0]   ->Fill(e2, l02/ncells2);
              fhEMCLambda1[0]    ->Fill(e2, l12);
              fhEMCdLambda1[0]   ->Fill(e2, l12/ncells2);
              fhEMCDispersion[0] ->Fill(e2, disp2);
              fhEMCdDispersion[0]->Fill(e2, disp2/ncells2);            
            }//photon   no conversion
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCElectron)){
              fhEMCLambda0[3]    ->Fill(e2, l02);
              fhEMCdLambda0[3]   ->Fill(e2, l02/ncells2);
              fhEMCLambda1[3]    ->Fill(e2, l12);
              fhEMCdLambda1[3]   ->Fill(e2, l12/ncells2);
              fhEMCDispersion[3] ->Fill(e2, disp2);
              fhEMCdDispersion[3]->Fill(e2, disp2/ncells2);    
            }//electron
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && 
                      GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) ){
              fhEMCLambda0[1]    ->Fill(e2, l02);
              fhEMCdLambda0[1]   ->Fill(e2, l02/ncells2);
              fhEMCLambda1[1]    ->Fill(e2, l12);
              fhEMCdLambda1[1]   ->Fill(e2, l12/ncells2);
              fhEMCDispersion[1] ->Fill(e2, disp2);
              fhEMCdDispersion[1]->Fill(e2, disp2/ncells2);                
            }//conversion photon
            else if  ( GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCPi0)  ){
              fhEMCLambda0[2]    ->Fill(e2, l02);
              fhEMCdLambda0[2]   ->Fill(e2, l02/ncells2);
              fhEMCLambda1[2]    ->Fill(e2, l12);
              fhEMCdLambda1[2]   ->Fill(e2, l12/ncells2);
              fhEMCDispersion[2] ->Fill(e2, disp2);
              fhEMCdDispersion[2]->Fill(e2, disp2/ncells2);                
            }//pi0
            else {
              fhEMCLambda0[4]    ->Fill(e2, l02);
              fhEMCdLambda0[4]   ->Fill(e2, l02/ncells2);
              fhEMCLambda1[4]    ->Fill(e2, l12);
              fhEMCdLambda1[4]   ->Fill(e2, l12/ncells2);
              fhEMCDispersion[4] ->Fill(e2, disp2);
              fhEMCdDispersion[4]->Fill(e2, disp2/ncells2);                
            }//other particles 
          }//is datamc
        }//MC histograms
        
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
  if(!GetInputAODBranch()){
    printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input calo photons in AOD branch with name < %s > , STOP\n",GetInputAODName().Data());
    abort();
  }
  
  for(Int_t iphoton = 0; iphoton < GetInputAODBranch()->GetEntriesFast(); iphoton++){
    AliAODPWG4Particle * photon1 =  (AliAODPWG4Particle*) (GetInputAODBranch()->At(iphoton));
    mom1 = *(photon1->Momentum());
    
    //Play with the MC stack if available
    fInputAODGammaConv = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaConvName);
    if(!fInputAODGammaConv) {
      printf("AliAnaPi0EbE::MakeInvMassInCalorimeterAndCTS() - No input gamma conversions in AOD branch with name < %s >, STOP\n",fInputAODGammaConvName.Data());
      abort();	
    }
    for(Int_t jphoton = 0; jphoton < fInputAODGammaConv->GetEntriesFast(); jphoton++){
      AliAODPWG4Particle * photon2 =  (AliAODPWG4Particle*) (fInputAODGammaConv->At(jphoton));
      if(GetMixedEvent())
        evtIndex = GetMixedEvent()->EventIndexForCaloCluster(photon2->GetCaloLabel(0)) ;
      if(TMath::Abs(GetVertex(evtIndex)[2]) > GetZvertexCut()) continue ;  //vertex cut
      
      mom2 = *(photon2->Momentum());
      
      //Int_t input = -1;	//if -1 photons come from different files, not a pi0
      //if(photon1->GetInputFileIndex() == photon2->GetInputFileIndex()) input = photon1->GetInputFileIndex();
      
      //Select good pair (good phi, pt cuts, aperture and invariant mass)
      if(GetNeutralMesonSelection()->SelectPair(mom1, mom2)){
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
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodpi0.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      if(GetDebug() > 1) 
        printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - FillAOD: PDG of identified particle %d\n",aodpi0.GetIdentifiedParticleType());
      //If primary is not pi0, skip it.
      if(aodpi0.GetIdentifiedParticleType() != AliCaloPID::kPi0) continue ;
    }					
    else if(IsCaloPIDOn()){
      //Skip matched clusters with tracks
      if(IsTrackMatched(calo)) continue ;
      
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn())
        aodpi0.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      else
        aodpi0.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPi0EbE::MakeShowerShapeIdentification() - PDG of identified particle %d\n",aodpi0.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not pi0, skip it.
      if(aodpi0.GetIdentifiedParticleType() != AliCaloPID::kPi0) continue ;			
      
    }
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodpi0, GetCaloUtils());
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
