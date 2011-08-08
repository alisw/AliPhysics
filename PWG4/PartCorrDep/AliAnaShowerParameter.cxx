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
// Class cloned from AliAnaPhoton, main aim is shower shape studies
// 
// 
//
//-- Author: Jocelyn Mlynarz (WSU) and Gustavo Conesa (LPSC-CNRS)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system --- 
#include <TH3F.h>
#include <TH2F.h>
#include <TClonesArray.h>
//#include <TObjString.h>
#include <Riostream.h>
#include "TParticle.h"
//#include <fstream>

// --- Analysis system --- 
#include "AliAnaShowerParameter.h" 
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliAODCaloCluster.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskParticleCorrelation.h"
#include "AliAODEvent.h"


ClassImp(AliAnaShowerParameter)

//____________________________________________________________________________
AliAnaShowerParameter::AliAnaShowerParameter() : 
AliAnaPartCorrBaseClass(), fCalorimeter(""), fNCellsCutMin(0),
fNCellsCutMax(0), fLambdaCut(0), fTimeCutMin(-1), fTimeCutMax(9999999),
fhNClusters(0), fhNCellCluster(0), fhEtaPhiPtCluster(0), 
fhLambdaPtCluster(0), 

//MC

fhLambdaPtPhoton(0), fhLambdaPtPi0(0), fhLambdaPtPion(0), fhPtTruthPi0(0)

{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
}

//____________________________________________________________________________
AliAnaShowerParameter::~AliAnaShowerParameter() 
{
  //dtor
  
}

//________________________________________________________________________
TObjString *  AliAnaShowerParameter::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaShowerParameter ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNCellsCutMin: (Cut in the minimum number of cells) %e\n",fNCellsCutMin) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"fNCellsCutMax: (Cut in the maximum number of cells) %e\n",fNCellsCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fLambdaCut: (Cut in the minimum lambda) %e\n",fLambdaCut) ;
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
TList *  AliAnaShowerParameter::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("PhotonHistos") ; 
	
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  //General non-MC Cluster histograms
  fhNClusters = new TH1F ("hNClusters","NClusters",21,-0.5,20.5); 
  fhNClusters->SetXTitle("N_{Clusters}");
  outputContainer->Add(fhNClusters) ;
  
  fhNCellCluster = new TH2F ("hNCellCluster","Number of cell per cluster",nptbins,ptmin,ptmax,21,-0.5,20.5); 
  fhNCellCluster->SetYTitle("N_{Cells}");  
  fhNCellCluster->SetXTitle("p_{T}");
  outputContainer->Add(fhNCellCluster) ;
  
  fhEtaPhiPtCluster = new TH3F
  ("hEtaPhiPtCluster","#phi_{Cluster} and #eta_{Cluster}",nptbins,ptmin,ptmax,nphibins,phimin,phimax,netabins,etamin,etamax); 
  fhEtaPhiPtCluster->SetZTitle("#eta");  
  fhEtaPhiPtCluster->SetYTitle("#phi");
  fhEtaPhiPtCluster->SetXTitle("pT_{Cluster} (GeV/c)");
  outputContainer->Add(fhEtaPhiPtCluster) ; 
  
  fhLambdaPtCluster  = new TH2F
  ("hLambdaCluster","#lambda_{Cluster}",nptbins,ptmin,ptmax,300,0,3); 
  fhLambdaPtCluster->SetYTitle("#lambda_{0}^{2}");
  fhLambdaPtCluster->SetXTitle("p_{T Cluster} (GeV/c)");
  outputContainer->Add(fhLambdaPtCluster) ;
  
  if(IsDataMC()){
    
    fhLambdaPtPhoton  = new TH2F
    ("hLambdaPtPhoton","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2); 
    fhLambdaPtPhoton->SetYTitle("#lambda_{0}^{2}");
    fhLambdaPtPhoton->SetXTitle("pT_{#gamma, Reco} (GeV)");
    outputContainer->Add(fhLambdaPtPhoton) ;
    
    fhLambdaPtPi0  = new TH2F
    ("hLambdaPtPi0","#lambda_{#pi^{0}}",nptbins,ptmin,ptmax,200,0,2); 
    fhLambdaPtPi0->SetYTitle("#lambda_{0}^{2}");
    fhLambdaPtPi0->SetXTitle("pT_{#pi^{0}, Reco} (GeV)");
    outputContainer->Add(fhLambdaPtPi0) ;
    
    fhLambdaPtPion  = new TH2F
    ("hLambdaPtPion","#lambda_{#pi^{+}}",nptbins,ptmin,ptmax,200,0,2); 
    fhLambdaPtPion->SetYTitle("#lambda_{0}^{2}");
    fhLambdaPtPion->SetXTitle("pT_{#pi^{+}, Reco} (GeV)");
    outputContainer->Add(fhLambdaPtPion) ;
   
    fhPtTruthPi0  = new TH1D("hPtTruthPi0","#pi^{0} MC truth pT",nptbins,ptmin,ptmax) ;
    fhPtTruthPi0->SetXTitle("pT_{#pi^{0}, Reco} (GeV)");
    outputContainer->Add(fhPtTruthPi0) ;
    
  }//Histos with MC
  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaShowerParameter::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaShowerParameter::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaShowerParameter::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//____________________________________________________________________________
void AliAnaShowerParameter::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaLambda_");
  
  fCalorimeter = "EMCAL" ;
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCutMin = 0 ;
  fNCellsCutMax = 0 ;
  fLambdaCut = 0.01 ;
	
}

//__________________________________________________________________
void  AliAnaShowerParameter::MakeAnalysisFillHistograms() 
{
  
  //Do analysis and fill histograms
  if ( GetDebug() > 0) printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Starting analysis.\n");
  
  // Access MC information in stack if requested, check that it exists.	
  AliStack * stack = 0x0;
  TClonesArray * mcparticles0 = 0x0;
  Int_t NClusters = 0 ;
  TLorentzVector momCluster ;
  
  //Check if the stack is available when analysing MC data.
  if(IsDataMC()){
    
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;
      if(!stack) {
        printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
				abort();
      }
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      //Get the list of MC particles
      mcparticles0 = GetReader()->GetAODMCParticles(0);
      if(!mcparticles0 && GetDebug() > 0) 	{
        printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() -  Standard MCParticles not available !\n"); 
      }
    }// is data and MC
  }	
  //Loop on stored AOD photons
  
  TClonesArray *  clustArray = 0x0 ; 
  clustArray = GetAODCaloClusters() ;
  NClusters = clustArray->GetEntriesFast() ;    
  fhNClusters->Fill(NClusters) ;
  
  if(GetDebug() > 0) printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - aod branch entries %d\n", NClusters);
  
  for(Int_t iCluster = 0 ; iCluster < NClusters ; iCluster++){    
    Int_t input = 0 ;
    
    AliAODCaloCluster * caloCluster =  (AliAODCaloCluster*) (clustArray->At(iCluster)) ;
    caloCluster->GetMomentum(momCluster,GetVertex(0)) ;
    AliAODPWG4Particle * aodCluster = new AliAODPWG4Particle(momCluster) ;
    Int_t labelCluster = caloCluster->GetLabel() ;
    aodCluster->SetLabel(labelCluster) ;;
    aodCluster->SetInputFileIndex(input) ;
    aodCluster->SetCaloLabel(caloCluster->GetID(),-1) ;
    aodCluster->SetDetector(fCalorimeter) ;
    aodCluster->SetTag(GetMCAnalysisUtils()->CheckOrigin(caloCluster->GetLabels(),caloCluster->GetNLabels(),GetReader(), aodCluster->GetInputFileIndex())) ;
    
    if (GetDebug() > 3) printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Cluster particle label = %d\n",labelCluster) ;
    
    //Fill Cluster histograms 
    Float_t ptcluster  = aodCluster->Pt() ;
    Float_t phicluster = aodCluster->Phi() ;
    Float_t etacluster = aodCluster->Eta() ;
    Float_t lambdaMainCluster = caloCluster->GetM02() ;
    Int_t iNumCell = caloCluster->GetNCells() ;	
    
    if (iNumCell>=fNCellsCutMin&&iNumCell<=fNCellsCutMax&&lambdaMainCluster>=fLambdaCut){
      
      //Fill the basic non-MC information about the cluster.
      fhNCellCluster->Fill(ptcluster,iNumCell) ;
      fhEtaPhiPtCluster->Fill(ptcluster,phicluster,etacluster) ;
      fhLambdaPtCluster->Fill(ptcluster,lambdaMainCluster) ;
      
      //Play with the MC data if available
      if(IsDataMC()){
        
        //Get the tag from AliMCAnalysisUtils for PID
        Int_t tag = aodCluster->GetTag();
        if (GetDebug() > 3) printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Cluster particle tag= %d\n",tag) ;
        if ( aodCluster->GetLabel() < 0){
          if(GetDebug() > 0) 
            printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Label is -1; problem in the MC ESD? ");
          continue;
        }
        
        //Check if the tag matches one of the different particle types and fill the corresponding histograms
        if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)&&!(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))) {
          fhLambdaPtPhoton->Fill(ptcluster,lambdaMainCluster) ;
        }//kMCPhoton
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)) {     
          fhLambdaPtPi0->Fill(ptcluster,lambdaMainCluster) ;        
        }
        if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPion)) {
          fhLambdaPtPion->Fill(ptcluster,lambdaMainCluster) ;
        }
        
        // Access MC information in stack if requested, check that it exists.
        Int_t label = aodCluster->GetLabel();
        if(label < 0) {
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
          continue;
        }
      }
    }
  }
  
  Float_t fPtGener = 0 ;
  if(IsDataMC() && stack){
    TParticle* pGener = stack->Particle(0) ;
    fPtGener = pGener->Pt() ;
  }
  
  
}


//__________________________________________________________________
void AliAnaShowerParameter::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("Min number of cells in cluster is        > %f \n", fNCellsCutMin);
  printf("Max number of cells in cluster is        > %f \n", fNCellsCutMax);
  printf("Min lambda of cluster is        > %f \n", fLambdaCut);
  printf("    \n") ;
  
} 
