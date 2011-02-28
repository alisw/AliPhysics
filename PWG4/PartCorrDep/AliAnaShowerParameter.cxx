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
#include "AliAnaShowerParameter.h" 
#include "AliEMCALGeoUtils.h"
#include "AliAODCaloCells.h"
#include "AliAODEvent.h"


ClassImp(AliAnaShowerParameter)

//____________________________________________________________________________
AliAnaShowerParameter::AliAnaShowerParameter() : 
AliAnaPartCorrBaseClass(), fCalorimeter(""), 
fMinDist(0.),fMinDist2(0.),fMinDist3(0.),fRejectTrackMatch(0),
fCheckConversion(kFALSE),fAddConvertedPairsToAOD(kFALSE), fMassCut(0), fNCellsCut(0),
fTimeCutMin(-1), fTimeCutMax(9999999),fEMCALGeo(0),fNumClusters(0),
fhNClusters(0),fhNCellCluster(0),fhPtCluster(0),fhPhiCluster(0),fhEtaCluster(0),fhDeltaPhiClusters(0),fhDeltaEtaClusters(0),
fhLambdaCluster(0),fhDispersionCluster(0),fhELambdaCluster(0),fhELambdaCellCluster(0),
fhNCellPhoton(0),fhLambdaPhoton(0),fhDispersionPhoton(0),fhELambdaPhoton(0),fhELambdaCellPhoton(0),
fhNCellPi0(0),fhLambdaPi0(0),fhDispersionPi0(0),fhELambdaPi0(0),fhELambdaCellPi0(0),
fhNCellChargedHadron(0),fhLambdaChargedHadron(0),fhDispersionChargedHadron(0),fhELambdaChargedHadron(0),fhELambdaCellChargedHadron(0),  
//MC
fhDeltaE(0), fhDeltaPt(0),fhRatioE(0), fhRatioPt(0),fh2E(0),fh2Pt(0),fhMCPdg(0),
fhEMCPhoton(0),fhPhiMCPhoton(0),fhEtaMCPhoton(0),fhLambdaMCPhoton(0),fhPhotTrueE(0),fhRatioEPhoton(0),
fhEMCPi0(0),fhPhiMCPi0(0),fhEtaMCPi0(0),fhLambdaMCPi0(0),fhPi0TrueE(0),
fhProductionDistance(0),fhProductionRadius(0),fhDecayAngle(0),fhRatioEPi0(0),
fhEMCPion(0),fhPhiMCPion(0),fhEtaMCPion(0),fhLambdaMCPion(0),fhPionTrueE(0),fhRatioEPion(0), 
fhEMCProton(0),fhPhiMCProton(0),fhEtaMCProton(0),fhLambdaMCProton(0),fhProtonTrueE(0),fhRatioEProton(0),
fhEMCAntiProton(0),fhPhiMCAntiProton(0),fhEtaMCAntiProton(0),fhLambdaMCAntiProton(0),fhAntiProtonTrueE(0),fhRatioEAntiProton(0),
fhEMCNeutron(0),fhPhiMCNeutron(0),fhEtaMCNeutron(0),fhLambdaMCNeutron(0),fhNeutronTrueE(0),fhRatioENeutron(0),
fhEMCEta(0),fhPhiMCEta(0),fhEtaMCEta(0),fhLambdaMCEta(0),
fhPrimPt(0x0)
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
  snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRejectTrackMatch: %d\n",fRejectTrackMatch) ;
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
  fhNCellCluster->SetXTitle("p_{T}");
  fhNCellCluster->SetYTitle("N_{Cells}");  
  outputContainer->Add(fhNCellCluster) ;
  
  fhPtCluster  = new TH1F("hPtCluster","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
  fhPtCluster->SetXTitle("p_{Cluster}(GeV/c)");
  outputContainer->Add(fhPtCluster) ; 
  
  fhPhiCluster  = new TH2F
  ("hPhiCluster","#phi_{Cluster}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiCluster->SetYTitle("#phi");
  fhPhiCluster->SetXTitle("E_{Cluster} (GeV/c)");
  outputContainer->Add(fhPhiCluster) ; 
  
  fhEtaCluster  = new TH2F
  ("hEtaCluster","#phi_{Cluster}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaCluster->SetYTitle("#eta");  
  fhEtaCluster->SetXTitle("E_{Cluster} (GeV/c)");
  outputContainer->Add(fhEtaCluster) ;

  fhDeltaEtaClusters  = new TH1F("hDeltaEtaClusters","#Delta #eta of clusters over calorimeter",netabins,etamin,etamax); 
  fhDeltaEtaClusters->SetXTitle("#Delta #eta_{Clusters}");
  outputContainer->Add(fhDeltaEtaClusters) ; 

  fhDeltaPhiClusters  = new TH1F("hDeltaPhiClusters","#Delta #phi of clusters over calorimeter",nphibins,phimin,phimax); 
  fhDeltaPhiClusters->SetXTitle("#Delta #phi_{Clusters}");
  outputContainer->Add(fhDeltaPhiClusters) ; 
  
  fhLambdaCluster  = new TH3F
  ("hLambdaCluster","#lambda_{Cluster}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhLambdaCluster->SetZTitle("#lambda_{1}^{2}"); 
  fhLambdaCluster->SetYTitle("#lambda_{0}^{2}");
  fhLambdaCluster->SetXTitle("p_{T Cluster} (GeV/c)");
  outputContainer->Add(fhLambdaCluster) ;
  
  fhELambdaCluster  = new TH3F
  ("hELambdaCluster","#lambda_{Cluster}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhELambdaCluster->SetZTitle("#lambda_{1}^{2}"); 
  fhELambdaCluster->SetYTitle("#lambda_{0}^{2}");
  fhELambdaCluster->SetXTitle("p_{T Cluster} (GeV)");
  outputContainer->Add(fhELambdaCluster) ;
  
  fhDispersionCluster  = new TH2F
  ("hDispersionCluster","Dispersion_{Cluster}",nptbins,ptmin,ptmax,200,0,2); 
  fhDispersionCluster->SetYTitle("Dispersion");
  fhDispersionCluster->SetXTitle("p_{T Cluster} (GeV/c)");
  outputContainer->Add(fhDispersionCluster) ;
  
  fhELambdaCellCluster  = new TH3F
  ("hELambdaCellCluster","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
  fhELambdaCellCluster->SetZTitle("N Cells"); 
  fhELambdaCellCluster->SetYTitle("#lambda_{0}^{2}");
  fhELambdaCellCluster->SetXTitle("E_{Cluster} (GeV)");
  outputContainer->Add(fhELambdaCellCluster) ;
  
  //Identified photon histograms
  fhNCellPhoton = new TH2F ("hNCellPhoton","Number of cell per cluster",nptbins,ptmin,ptmax,21,-0.5,20.5); 
  fhNCellPhoton->SetXTitle("E");
  fhNCellPhoton->SetYTitle("N_{Cells}");  
  outputContainer->Add(fhNCellPhoton) ;
  
  fhLambdaPhoton  = new TH3F
  ("hLambdaPhoton","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhLambdaPhoton->SetZTitle("#lambda_{1}^{2}"); 
  fhLambdaPhoton->SetYTitle("#lambda_{0}^{2}");
  fhLambdaPhoton->SetXTitle("p_{T Photon} (GeV/c)");
  outputContainer->Add(fhLambdaPhoton) ;
  
  fhELambdaPhoton  = new TH3F
  ("hELambdaPhoton","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhELambdaPhoton->SetZTitle("#lambda_{1}^{2}"); 
  fhELambdaPhoton->SetYTitle("#lambda_{0}^{2}");
  fhELambdaPhoton->SetXTitle("E_{Photon} (GeV)");
  outputContainer->Add(fhELambdaPhoton) ;
  
  fhDispersionPhoton  = new TH2F
  ("hLambdaDispersion","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2); 
  fhDispersionPhoton->SetYTitle("Dispersion");
  fhDispersionPhoton->SetXTitle("p_{T Photon} (GeV/c)");
  outputContainer->Add(fhDispersionPhoton) ;
  
  fhELambdaCellPhoton  = new TH3F
  ("hELambdaCellPhoton","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
  fhELambdaCellPhoton->SetZTitle("N Cells"); 
  fhELambdaCellPhoton->SetYTitle("#lambda_{0}^{2}");
  fhELambdaCellPhoton->SetXTitle("E_{Photon} (GeV)");
  outputContainer->Add(fhELambdaCellPhoton) ;
  
  //Identified Pi0 histograms
  fhNCellPi0 = new TH2F ("hNCellPi0","Number of cell per cluster",nptbins,ptmin,ptmax,21,-0.5,20.5); 
  fhNCellPi0->SetXTitle("E");
  fhNCellPi0->SetYTitle("N_{Cells}");  
  outputContainer->Add(fhNCellPi0) ;
  
  fhLambdaPi0  = new TH3F
  ("hLambdaPi0","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhLambdaPi0->SetZTitle("#lambda_{1}^{2}"); 
  fhLambdaPi0->SetYTitle("#lambda_{0}^{2}");
  fhLambdaPi0->SetXTitle("p_{T Pi0} (GeV/c)");
  outputContainer->Add(fhLambdaPi0) ;
  
  fhELambdaPi0  = new TH3F
  ("hELambdaPi0","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhELambdaPi0->SetZTitle("#lambda_{1}^{2}"); 
  fhELambdaPi0->SetYTitle("#lambda_{0}^{2}");
  fhELambdaPi0->SetXTitle("E_{Pi0} (GeV)");
  outputContainer->Add(fhELambdaPi0) ;
  
  fhDispersionPi0  = new TH2F
  ("hLambdaDispersion","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2); 
  fhDispersionPi0->SetYTitle("Dispersion");
  fhDispersionPi0->SetXTitle("p_{T Pi0} (GeV/c)");
  outputContainer->Add(fhDispersionPi0) ;
  
  fhELambdaCellPi0  = new TH3F
  ("hELambdaCellPi0","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
  fhELambdaCellPi0->SetZTitle("N Cells"); 
  fhELambdaCellPi0->SetYTitle("#lambda_{0}^{2}");
  fhELambdaCellPi0->SetXTitle("E_{Pi0} (GeV)");
  outputContainer->Add(fhELambdaCellPi0) ;
  
  //Identified charged hadron histograms
  fhNCellChargedHadron = new TH2F ("hNCellChargedHadron","Number of cell per cluster",nptbins,ptmin,ptmax,21,-0.5,20.5); 
  fhNCellChargedHadron->SetXTitle("p_{T}");
  fhNCellChargedHadron->SetYTitle("N_{Cells}");  
  outputContainer->Add(fhNCellChargedHadron) ;
  
  fhLambdaChargedHadron  = new TH3F
  ("hLambdaChargedHadron","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhLambdaChargedHadron->SetZTitle("#lambda_{1}^{2}"); 
  fhLambdaChargedHadron->SetYTitle("#lambda_{0}^{2}");
  fhLambdaChargedHadron->SetXTitle("p_{T ChargedHadron} (GeV/c)");
  outputContainer->Add(fhLambdaChargedHadron) ;
  
  fhELambdaChargedHadron  = new TH3F
  ("hELambdaChargedHadron","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,200,0,2); 
  fhELambdaChargedHadron->SetZTitle("#lambda_{1}^{2}"); 
  fhELambdaChargedHadron->SetYTitle("#lambda_{0}^{2}");
  fhELambdaChargedHadron->SetXTitle("E_{ChargedHadron} (GeV)");
  outputContainer->Add(fhELambdaChargedHadron) ;
  
  fhDispersionChargedHadron  = new TH2F
  ("hLambdaDispersion","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2); 
  fhDispersionChargedHadron->SetYTitle("Dispersion");
  fhDispersionChargedHadron->SetXTitle("p_{T ChargedHadron} (GeV/c)");
  outputContainer->Add(fhDispersionChargedHadron) ;
  
  fhELambdaCellChargedHadron  = new TH3F
  ("hELambdaCellChargedHadron","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
  fhELambdaCellChargedHadron->SetZTitle("N Cells"); 
  fhELambdaCellChargedHadron->SetYTitle("#lambda_{0}^{2}");
  fhELambdaCellChargedHadron->SetXTitle("E_{ChargedHadron} (GeV)");
  outputContainer->Add(fhELambdaCellChargedHadron) ;
  
  if(IsDataMC()){
    fhDeltaE  = new TH1F ("hDeltaE","MC - Reco E ", 200,-50,50); 
    fhDeltaE->SetXTitle("#Delta E (GeV)");
    outputContainer->Add(fhDeltaE);
    
    fhDeltaPt  = new TH1F ("hDeltaPt","MC - Reco p_{T} ", 200,-50,50); 
    fhDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
    outputContainer->Add(fhDeltaPt);
    
    fhRatioE  = new TH1F ("hRatioE","Reco/MC E ", 1000,0,2); 
    fhRatioE->SetXTitle("E_{reco}/E_{gen}");
    outputContainer->Add(fhRatioE);
    
    fhRatioPt  = new TH1F ("hRatioPt","Reco/MC p_{T} ", 1000,0,2); 
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

    fhMCPdg  = new TH2F ("hMCPdg","PDG Code distribution",nptbins,ptmin,ptmax,6001,-3000.5,3000.5); 
    fhMCPdg->SetXTitle("E_{Reco Cluster} (GeV)");
    fhMCPdg->SetYTitle("PDG Code");
    outputContainer->Add(fhMCPdg);
    
    fhEMCPhoton  = new TH1F("hEMCPhoton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    //fhEMCPhoton->SetYTitle("N");
    fhEMCPhoton->SetXTitle("E_{#gamma}(GeV)");
    outputContainer->Add(fhEMCPhoton) ; 
    
    fhPhiMCPhoton  = new TH2F
    ("hPhiMCPhoton","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCPhoton->SetYTitle("#phi");
    fhPhiMCPhoton->SetXTitle("E_{ #gamma} (GeV)");
    outputContainer->Add(fhPhiMCPhoton) ; 
    
    fhEtaMCPhoton  = new TH2F
    ("hEtaMCPhoton","#phi_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCPhoton->SetYTitle("#eta");
    fhEtaMCPhoton->SetXTitle("E_{ #gamma} (GeV)");
    outputContainer->Add(fhEtaMCPhoton) ;
    
    fhLambdaMCPhoton  = new TH3F
    ("hLambdaMCPhoton","#lambda_{#gamma}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCPhoton->SetZTitle("N_{Cells}");
    fhLambdaMCPhoton->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCPhoton->SetXTitle("E_{#gamma} (GeV)");
    outputContainer->Add(fhLambdaMCPhoton) ;
    
    fhPhotTrueE  = new TH3F
    ("hPhotTrueE","#lambda_{#gamma}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhPhotTrueE->SetZTitle("#lambda_{0}^{2}");
    fhPhotTrueE->SetYTitle("Recons. E_{#gamma}");
    fhPhotTrueE->SetXTitle("MC Truth E_{#gamma} (GeV)");
    outputContainer->Add(fhPhotTrueE) ;
    
    fhRatioEPhoton  = new TH1F ("hRatioEPhoton","E_{Reco #gamma}/E_{MC truth #gamma}", 1000,0,2); 
    fhRatioEPhoton->SetXTitle("E_{reco #gamma}/E_{gen #gamma}");
    outputContainer->Add(fhRatioEPhoton);
    
    fhEMCPi0  = new TH1F("hEMCPi0","Number of #pi^0 over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCPi0->SetYTitle("N");
    fhEMCPi0->SetXTitle("E_{ #pi^0}(GeV)");
    outputContainer->Add(fhEMCPi0) ; 
    
    fhPhiMCPi0  = new TH2F
    ("hPhiMCPi0","#phi_{#pi^0}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCPi0->SetYTitle("#phi");
    fhPhiMCPi0->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhPhiMCPi0) ; 
    
    fhEtaMCPi0  = new TH2F
    ("hEtaMCPi0","#phi_{#pi^0}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCPi0->SetYTitle("#eta");
    fhEtaMCPi0->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhEtaMCPi0) ;
    
    fhLambdaMCPi0  = new TH3F
    ("hLambdaMCPi0","#lambda_{#pi^0}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCPi0->SetZTitle("N_{Cells}");
    fhLambdaMCPi0->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCPi0->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhLambdaMCPi0) ;
    
    fhProductionDistance  = new TH2F
    ("hProductionDistance","#lambda_{#pi^0}",nptbins,ptmin,ptmax,160,0,800); 
    fhProductionDistance->SetYTitle("Distance of production of Pi0 from vertex (cm)");
    fhProductionDistance->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhProductionDistance) ;
    
    fhProductionRadius  = new TH2F
    ("hProductionRadius","#lambda_{#pi^0}",nptbins,ptmin,ptmax,160,0,800); 
    fhProductionRadius->SetYTitle("Radius of production of Pi0 from vertex (cm)");
    fhProductionRadius->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhProductionRadius) ;
    
    fhDecayAngle = new TH2F
    ("hDecayAngle","#lambda_{#pi^0}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhDecayAngle->SetYTitle("Opening angle of the #pi^{0} decay (radians)");
    fhDecayAngle->SetXTitle("E_{ #pi^0} (GeV)");
    outputContainer->Add(fhDecayAngle) ;
    
    fhPi0TrueE  = new TH3F
    ("hPi0TrueE","#lambda_{#pi^0}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhPi0TrueE->SetZTitle("#lambda_{0}^{2}");
    fhPi0TrueE->SetYTitle("Recons. E_{#pi^0}");
    fhPi0TrueE->SetXTitle("MC Truth E_{#pi^0} (GeV)");
    outputContainer->Add(fhPi0TrueE) ;
    
    fhRatioEPi0  = new TH1F ("hRatioEPi0","E_{Reco #pi^{0}}/E_{MC truth #pi^{0}}", 1000,0,2); 
    fhRatioEPi0->SetXTitle("E_{reco #pi^{0}}/E_{gen #pi^{0}}");
    outputContainer->Add(fhRatioEPi0);
    
    fhEMCPion  = new TH1F("hEMCPion","Number of #pi over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCPion->SetYTitle("N");
    fhEMCPion->SetXTitle("E_{ #pi}(GeV)");
    outputContainer->Add(fhEMCPion) ; 
    
    fhPhiMCPion  = new TH2F
    ("hPhiMCPion","#phi_{#pi}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCPion->SetYTitle("#phi");
    fhPhiMCPion->SetXTitle("E_{ #pi} (GeV)");
    outputContainer->Add(fhPhiMCPion) ; 
    
    fhEtaMCPion  = new TH2F
    ("hEtaMCPion","#phi_{#pi}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCPion->SetYTitle("#eta");
    fhEtaMCPion->SetXTitle("E_{ #pi} (GeV)");
    outputContainer->Add(fhEtaMCPion) ;
    
    fhLambdaMCPion  = new TH3F
    ("hLambdaMCPion","#lambda_{#pi}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCPion->SetZTitle("N_{Cells}");
    fhLambdaMCPion->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCPion->SetXTitle("E_{ #pi} (GeV)");
    outputContainer->Add(fhLambdaMCPion) ;
    
    fhPionTrueE  = new TH3F
    ("hPionTrueE","#lambda_{#pi}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhPionTrueE->SetZTitle("#lambda_{0}^{2}");
    fhPionTrueE->SetYTitle("Recons. E_{#pi}");
    fhPionTrueE->SetXTitle("MC Truth E_{#pi} (GeV)");
    outputContainer->Add(fhPionTrueE) ;
    
    fhRatioEPion  = new TH1F ("hRatioEPion","E_{Reco #pi^{#pm}}/E_{MC truth #pi^{#pm}}", 1000,0,2); 
    fhRatioEPion->SetXTitle("E_{reco #pi^{#pm}}/E_{gen #pi^{#pm}}");
    outputContainer->Add(fhRatioEPion);
    
    fhEMCProton  = new TH1F("hEMCProton","Number of p over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCProton->SetYTitle("N");
    fhEMCProton->SetXTitle("E_{ p}(GeV)");
    outputContainer->Add(fhEMCProton) ; 
    
    fhPhiMCProton  = new TH2F
    ("hPhiMCProton","#phi_{p}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCProton->SetYTitle("#phi");
    fhPhiMCProton->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhPhiMCProton) ; 
    
    fhEtaMCProton  = new TH2F
    ("hEtaMCProton","#phi_{p}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCProton->SetYTitle("#eta");
    fhEtaMCProton->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhEtaMCProton) ;
    
    fhLambdaMCProton  = new TH3F
    ("hLambdaMCProton","#lambda_{p}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCProton->SetZTitle("N_{Cells}");
    fhLambdaMCProton->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCProton->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhLambdaMCProton) ;
    
    fhProtonTrueE  = new TH3F
    ("hProtonTrueE","#lambda_{p}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhProtonTrueE->SetZTitle("#lambda_{0}^{2}");
    fhProtonTrueE->SetYTitle("Recons. E_{p}");
    fhProtonTrueE->SetXTitle("MC Truth E_{p} (GeV)");
    outputContainer->Add(fhProtonTrueE) ;
    
    fhRatioEProton  = new TH1F ("hRatioEProton","E_{Reco p}/E_{MC truth p}", 1000,0,2); 
    fhRatioEProton->SetXTitle("E_{reco p}/E_{gen p}");
    outputContainer->Add(fhRatioEProton);
    
    fhEMCAntiProton  = new TH1F("hEMCAntiProton","Number of #bar{p} over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCAntiProton->SetYTitle("N");
    fhEMCAntiProton->SetXTitle("E_{#bar{p}}(GeV)");
    outputContainer->Add(fhEMCAntiProton) ; 
    
    fhPhiMCAntiProton  = new TH2F
    ("hPhiMCAntiProton","#phi_{#bar{p}}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCAntiProton->SetYTitle("#phi");
    fhPhiMCAntiProton->SetXTitle("E_{#bar{p}} (GeV)");
    outputContainer->Add(fhPhiMCAntiProton) ; 
    
    fhEtaMCAntiProton  = new TH2F
    ("hEtaMCAntiProton","#phi_{#bar{p}}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCAntiProton->SetYTitle("#eta");
    fhEtaMCAntiProton->SetXTitle("E_{#bar{p}} (GeV)");
    outputContainer->Add(fhEtaMCAntiProton) ;
    
    fhLambdaMCAntiProton  = new TH3F
    ("hLambdaMCAntiProton","#lambda_{#bar{p}}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCAntiProton->SetZTitle("N_{Cells}");
    fhLambdaMCAntiProton->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCAntiProton->SetXTitle("E_{#bar{p}} (GeV)");
    outputContainer->Add(fhLambdaMCAntiProton) ;
    
    fhAntiProtonTrueE  = new TH3F
    ("hAntiProtonTrueE","#lambda_{#bar{p}}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhAntiProtonTrueE->SetZTitle("#lambda_{0}^{2}");
    fhAntiProtonTrueE->SetYTitle("Recons. E_{#bar{p}}");
    fhAntiProtonTrueE->SetXTitle("MC Truth E_{#bar{p}} (GeV)");
    outputContainer->Add(fhAntiProtonTrueE) ;
    
    fhRatioEAntiProton  = new TH1F ("hRatioEAntiProton","E_{Reco #bar{p}}/E_{MC truth #bar{p}}", 1000,0,2); 
    fhRatioEAntiProton->SetXTitle("E_{reco #bar{p}}/E_{gen #bar{p}}");
    outputContainer->Add(fhRatioEAntiProton);
    
    fhEMCNeutron  = new TH1F("hEMCNeutron","Number of p over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCNeutron->SetYTitle("N");
    fhEMCNeutron->SetXTitle("E_{ p}(GeV)");
    outputContainer->Add(fhEMCNeutron) ; 
    
    fhPhiMCNeutron  = new TH2F
    ("hPhiMCNeutron","#phi_{p}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCNeutron->SetYTitle("#phi");
    fhPhiMCNeutron->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhPhiMCNeutron) ; 
    
    fhEtaMCNeutron  = new TH2F
    ("hEtaMCNeutron","#phi_{p}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCNeutron->SetYTitle("#eta");
    fhEtaMCNeutron->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhEtaMCNeutron) ;
    
    fhLambdaMCNeutron  = new TH3F
    ("hLambdaMCNeutron","#lambda_{p}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCNeutron->SetZTitle("N_{Cells}");
    fhLambdaMCNeutron->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCNeutron->SetXTitle("E_{ p} (GeV)");
    outputContainer->Add(fhLambdaMCNeutron) ;
    
    fhNeutronTrueE  = new TH3F
    ("hNeutronTrueE","#lambda_{n}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax,200,0,2); 
    fhNeutronTrueE->SetZTitle("#lambda_{0}^{2}");
    fhNeutronTrueE->SetYTitle("Recons. E_{n}");
    fhNeutronTrueE->SetXTitle("MC Truth E_{n} (GeV)");
    outputContainer->Add(fhNeutronTrueE) ;
    
    fhRatioENeutron  = new TH1F ("hRatioENeutron","E_{Reco n}/E_{MC truth n}", 1000,0,2); 
    fhRatioENeutron->SetXTitle("E_{reco n}/E_{gen n}");
    outputContainer->Add(fhRatioENeutron);
    
    fhEMCEta  = new TH1F("hEMCEta","Number of #eta over calorimeter",nptbins,ptmin,ptmax); 
    fhEMCEta->SetYTitle("N");
    fhEMCEta->SetXTitle("E_{ #eta}(GeV)");
    outputContainer->Add(fhEMCEta) ; 
    
    fhPhiMCEta  = new TH2F
    ("hPhiMCEta","#phi_{#eta}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCEta->SetYTitle("#phi");
    fhPhiMCEta->SetXTitle("E_{ #eta} (GeV)");
    outputContainer->Add(fhPhiMCEta) ; 
    
    fhEtaMCEta  = new TH2F
    ("hEtaMCEta","#phi_{#eta}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCEta->SetYTitle("#eta");
    fhEtaMCEta->SetXTitle("E_{ #eta} (GeV)");
    outputContainer->Add(fhEtaMCEta) ;
    
    fhLambdaMCEta  = new TH3F
    ("hLambdaMCEta","#lambda_{#eta}",nptbins,ptmin,ptmax,200,0,2,21,-0.5,20.5); 
    fhLambdaMCEta->SetZTitle("N_{Cells}");
    fhLambdaMCEta->SetYTitle("#lambda_{0}^{2}");
    fhLambdaMCEta->SetXTitle("E_{ #eta} (GeV)");
    outputContainer->Add(fhLambdaMCEta) ;
    
    fhPrimPt     = new TH1D("hPrimPt","Primary pi0 pt",nptbins,ptmin,ptmax) ;
    outputContainer->Add(fhPrimPt) ;
    
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
  AddToHistogramsName("AnaPhoton_");
  
  fCalorimeter = "PHOS" ;
  fMinDist  = 2.;
  fMinDist2 = 4.;
  fMinDist3 = 5.;
  fMassCut  = 0.03; //30 MeV
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCut = 0;
	
  fRejectTrackMatch       = kTRUE ;
  fCheckConversion        = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
  
  fNumClusters = -1; // By Default, select all events.
	
}

//__________________________________________________________________
void  AliAnaShowerParameter::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  //Search for photons in fCalorimeter 
  
  //Get vertex for photon momentum calculation
  
  for (Int_t iev = 0; iev < GetNMixedEvent(); iev++) {
    if (!GetMixedEvent()) 
      GetReader()->GetVertex(GetVertex(iev));
    else 
      GetMixedEvent()->GetVertexOfEvent(iev)->GetXYZ(GetVertex(iev)); 
  } 
  
  //Select the Calorimeter of the photon
  TObjArray * pl = 0x0; 
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  if(!pl){
    printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - Careful cluster array NULL!!\n");
    return;
  }
  
  //Fill AODCaloClusters and AODParticle with PHOS/EMCAL aods
  TLorentzVector mom, mom2 ; 
  Int_t nCaloClusters = pl->GetEntriesFast();   
  //Cut on the number of clusters in the event.
  if ((fNumClusters !=-1) && (nCaloClusters != fNumClusters)) return;
  Bool_t * indexConverted = new Bool_t[nCaloClusters];
  for (Int_t i = 0; i < nCaloClusters; i++) 
    indexConverted[i] = kFALSE;
  
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++){    
    
    AliAODCaloCluster * calo =  (AliAODCaloCluster*) (pl->At(icalo));	
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
    }
    //Cluster selection, not charged, with photon id and in fiducial cut
	  
    //Input from second AOD?
    Int_t input = 0;
    
    //Get Momentum vector, 
    if (input == 0) 
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;//Assume that come from vertex in straight line
    
    //Skip the cluster if it doesn't fit inside the cuts.
    if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 
    Double_t tof = calo->GetTOF()*1e9;    
    if(tof < fTimeCutMin || tof > fTimeCutMax) continue;	  
    if(calo->GetNCells() <= fNCellsCut) continue;
    
    //Check acceptance selection
    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue ;
    }
    
    //Create AOD for analysis
    AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
    Int_t label = calo->GetLabel();
    aodph.SetLabel(label);
    aodph.SetInputFileIndex(input);
    
    //Set the indices of the original caloclusters  
    aodph.SetCaloLabel(calo->GetID(),-1);
    aodph.SetDetector(fCalorimeter);
    if(GetDebug() > 1) 
      printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",aodph.Pt(),aodph.Phi(),aodph.Eta());	
    
    //Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist) //In bad channel (PHOS cristal size 2.2x2.2 cm)
      continue ;
    
    if(GetDebug() > 1) printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - Bad channel cut passed %4.2f\n",distBad);
    
    if     (distBad > fMinDist3) aodph.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodph.SetDistToBad(1) ; 
    else                         aodph.SetDistToBad(0) ;
    
    //Skip matched clusters with tracks
    if(fRejectTrackMatch && IsTrackMatched(calo)) continue ;
    if(GetDebug() > 1) printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - TrackMatching cut passed \n");
    
    //Set PID bits for later selection (AliAnaPi0 for example)
    //GetPDG already called in SetPIDBits.
    GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph, GetCaloUtils());
    if(GetDebug() > 1) printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - PID Bits set \n");		    
    
    //Play with the MC stack if available
    //Check origin of the candidates
    if(IsDataMC()){
      aodph.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex()));
      if(GetDebug() > 0) printf("AliAnaShowerParameter::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }
    
    // Check if cluster comes from a conversion in the material in front of the calorimeter
    // Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
    
    if(fCheckConversion && nCaloClusters > 1){
      Bool_t bConverted = kFALSE;
      Int_t id2 = -1;
		  
      //Check if set previously as converted couple, if so skip its use.
      if (indexConverted[icalo]) continue;
		  
      for(Int_t jcalo = icalo + 1 ; jcalo < nCaloClusters ; jcalo++) {
        //Check if set previously as converted couple, if so skip its use.
        if (indexConverted[jcalo]) continue;
        //printf("Check Conversion indeces %d and %d\n",icalo,jcalo);
        AliAODCaloCluster * calo2 =  (AliAODCaloCluster*) (pl->At(jcalo));              //Get cluster kinematics
        Int_t evtIndex2 = 0 ; 
        if (GetMixedEvent()) {
          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo2->GetID()) ; 
        }        
        calo2->GetMomentum(mom2,GetVertex(evtIndex2));
        //Check only certain regions
        Bool_t in2 = kTRUE;
        if(IsFiducialCutOn()) in2 =  GetFiducialCut()->IsInFiducialCut(mom2,fCalorimeter) ;
        if(!in2) continue;      
        
        //Get mass of pair, if small, take this pair.
        //printf("\t both in calo, mass %f, cut %f\n",(mom+mom2).M(),fMassCut);
        if((mom+mom2).M() < fMassCut){  
          bConverted = kTRUE;
          id2 = calo2->GetID();
          indexConverted[jcalo]=kTRUE;
          break;
        }
			  
      }//Mass loop
		  
      if(bConverted){ 
        if(fAddConvertedPairsToAOD){
          //Create AOD of pair analysis
          TLorentzVector mpair = mom+mom2;
          AliAODPWG4Particle aodpair = AliAODPWG4Particle(mpair);
          aodpair.SetLabel(aodph.GetLabel());
          aodpair.SetInputFileIndex(input);
          
          //printf("Index %d, Id %d\n",icalo, calo->GetID());
          //Set the indeces of the original caloclusters  
          aodpair.SetCaloLabel(calo->GetID(),id2);
          aodpair.SetDetector(fCalorimeter);
          aodpair.SetPdg(aodph.GetPdg());
          aodpair.SetTag(aodph.GetTag());
          
          //Add AOD with pair object to aod branch
          AddAODParticle(aodpair);
          //printf("\t \t both added pair\n");
        }
        
        //Do not add the current calocluster
        continue;
      }//converted pair
    }//check conversion
	  
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop;
  delete [] indexConverted;
	
  if(GetDebug() > 1) printf("AliAnaShowerParameter::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
void  AliAnaShowerParameter::MakeAnalysisFillHistograms() 
{
  
  //Do analysis and fill histograms
  
  // Access MC information in stack if requested, check that it exists.	
  AliStack * stack = 0x0;
  TParticle * primary = 0x0;   
  TClonesArray * mcparticles0 = 0x0;
  //TClonesArray * mcparticles1 = 0x0;
  AliAODMCParticle * aodprimary = 0x0; 
  TObjArray * pl = 0x0;
  Int_t iNumCell=0;
  
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
        printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n"); 
      }
    }// is data and MC
  }	
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  fhNClusters->Fill(naod);
  if(GetDebug() > 0) printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    
    if(ph->GetDetector() != fCalorimeter) continue;
    
    if(GetDebug() > 2) 
      printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
    
    //Fill Cluster histograms 
    Float_t ptcluster  = ph->Pt();
    Float_t phicluster = ph->Phi();
    Float_t etacluster = ph->Eta();
    Float_t ecluster   = ph->E();
    Float_t lambdaMainCluster = -1; //Can't set their values here, but need a default -1 value in case something is wrong
    Float_t lambdaSecondCluster = -1;
    Float_t dispcluster = -1;
    
    //Get the list of clusters from the AOD and loop over them to find the onces corresponding to the reconstructed 'photons'.
    if(fCalorimeter == "PHOS")
      pl = GetPHOSClusters();
    else if (fCalorimeter == "EMCAL")
      pl = GetEMCALClusters();
    if(pl){
      //Some values are stored in AliAODCaloCluster objects only; we need to fetch them.
      for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++){
        AliAODCaloCluster * calo =  (AliAODCaloCluster*) (pl->At(icalo));
        if (calo->GetLabel()==ph->GetLabel()) {  //The Cluster is the right one for this particle
          lambdaMainCluster = calo->GetM02();   //lambda_0
          lambdaSecondCluster = calo->GetM20();     //lambda_1
          dispcluster = calo->GetDispersion();    //Dispersion
          iNumCell = calo->GetNCells();	
          if(GetDebug() > 2) 
            printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Cluster Lambda0: %3.2f, Lambda1: %3.2f, Dispersion: %3.2f, NCells: %d \n",lambdaMainCluster,lambdaSecondCluster,dispcluster,iNumCell) ;
        }
      }
    }
    
    //Fill the basic non-MC information about the cluster.
    fhPtCluster  ->Fill(ecluster);
    fhPhiCluster ->Fill(ecluster,phicluster);
    fhEtaCluster ->Fill(ecluster,etacluster);
    fhDispersionCluster->Fill(ptcluster,dispcluster);
    fhLambdaCluster ->Fill(ptcluster,lambdaMainCluster,lambdaSecondCluster);
    fhELambdaCluster ->Fill(ecluster,lambdaMainCluster,lambdaSecondCluster);
    fhELambdaCellCluster ->Fill(ecluster,lambdaMainCluster,iNumCell);
    fhNCellCluster->Fill(ecluster,iNumCell);
    
    //In the case of an event with several clusters, calculate DeltaEta and DeltaPhi for the cluster pairs.
    for (Int_t jaod=iaod+1;jaod<naod;jaod++){
      AliAODPWG4Particle* phSecond = (AliAODPWG4Particle*) (GetOutputAODBranch()->At(jaod));	
      if(phSecond->GetDetector() != fCalorimeter) continue;
      fhDeltaPhiClusters->Fill(phicluster-(phSecond->Phi()));
      fhDeltaEtaClusters->Fill(etacluster-(phSecond->Eta()));  
    }
    
    //Fill the lambda distribution for identified particles
    if(ph->GetPdg() ==  AliCaloPID::kPhoton){
      fhDispersionPhoton->Fill(ptcluster,dispcluster);
      fhLambdaPhoton ->Fill(ptcluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaPhoton ->Fill(ecluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaCellPhoton ->Fill(ecluster,lambdaMainCluster,iNumCell);
      fhNCellPhoton->Fill(ecluster,iNumCell);
    }
    if(ph->GetPdg() ==  AliCaloPID::kPi0){
      fhDispersionPi0->Fill(ptcluster,dispcluster);
      fhLambdaPi0 ->Fill(ptcluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaPi0 ->Fill(ecluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaCellPi0 ->Fill(ecluster,lambdaMainCluster,iNumCell);
      fhNCellPi0->Fill(ecluster,iNumCell);
    }
    if(ph->GetPdg() ==  AliCaloPID::kChargedHadron){
      fhDispersionChargedHadron->Fill(ptcluster,dispcluster);
      fhLambdaChargedHadron ->Fill(ptcluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaChargedHadron ->Fill(ecluster,lambdaMainCluster,lambdaSecondCluster);
      fhELambdaCellChargedHadron ->Fill(ecluster,lambdaMainCluster,iNumCell);
      fhNCellChargedHadron->Fill(ecluster,iNumCell);
    }
    
    //Play with the MC data if available
    if(IsDataMC()){
      if(GetReader()->ReadStack() && !stack) return;
      if(GetReader()->ReadAODMCParticles() && !mcparticles0) return;

      //Get the tag from AliMCAnalysisUtils for PID
      Int_t tag =ph->GetTag();
      if ( ph->GetLabel() < 0){
        if(GetDebug() > 0) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - Label is -1; problem in the MC ESD? ");
        continue;
      }
      
      
      //Check if the tag matches on of the different particle types and fill the corresponding histograms
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)&&!(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))) //kMCPhoton; making sure that this is not a Pi0 cluster (it should not happen).
      {
        fhEMCPhoton  ->Fill(ecluster);
        fhPhiMCPhoton ->Fill(ecluster,phicluster);
        fhEtaMCPhoton ->Fill(ecluster,etacluster);
        fhLambdaMCPhoton ->Fill(ecluster,lambdaMainCluster,iNumCell);
        Double_t ePhot;
        //Get the true energy of the photon
        Int_t iCurrent = ph->GetLabel();
        TParticle* pCurrent = stack->Particle(iCurrent);
        ePhot = pCurrent->Energy(); 	  
        fhPhotTrueE->Fill(ePhot,ecluster,lambdaMainCluster);
        fhRatioEPhoton->Fill(ecluster/ePhot);
        fhMCPdg->Fill(ecluster,22);
        if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCPhoton: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",ePhot,ecluster,lambdaMainCluster);
      }//kMCPhoton
      
      if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)) //kMCPi0 : single cluster created by 2 photons from the same Pi0 decay
      {     
        
        //Fill the basic information about the Pi0 cluster
        fhEMCPi0  ->Fill(ecluster);
        fhPhiMCPi0 ->Fill(ecluster,phicluster);
        fhEtaMCPi0 ->Fill(ecluster,etacluster);
        fhLambdaMCPi0 ->Fill(ecluster,lambdaMainCluster,iNumCell);
        
        //Check the position of the Pi0: does it come from the vertex or was it created by hadron-material collision?
        Int_t iCurrent = ph->GetLabel();
        Double_t ePi0 = 0;
        TParticle* pCurrent = stack->Particle(iCurrent);
        while (pCurrent->GetPdgCode()!=111)
        {
          iCurrent = pCurrent->GetFirstMother();
          pCurrent = stack->Particle(iCurrent);
        }
        Float_t fDistance = pCurrent->Rho();
        Float_t fRadius = pCurrent->R();
        ePi0 = pCurrent->Energy();
        fhProductionDistance->Fill(ecluster,fDistance);
        fhProductionRadius->Fill(ecluster,fRadius);
        //Compare the cluster energy and true energy
        fhPi0TrueE->Fill(ePi0,ecluster,lambdaMainCluster);
        fhRatioEPi0->Fill(ecluster/ePi0);
        fhMCPdg->Fill(ecluster,111);
        if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCPi0: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",ePi0,ecluster,lambdaMainCluster);
      }
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPion))//kMCPion
      {
        fhEMCPion  ->Fill(ecluster);
        fhPhiMCPion ->Fill(ecluster,phicluster);
        fhEtaMCPion ->Fill(ecluster,etacluster);
        fhLambdaMCPion ->Fill(ecluster,lambdaMainCluster,iNumCell);
        
        Int_t iCurrent = ph->GetLabel();
        Double_t ePion = 0;
        TParticle* pCurrent = stack->Particle(iCurrent);
        while ((TMath::Abs(pCurrent->GetPdgCode())!=211)&&(iCurrent>=0))
        {
          if (iCurrent!=pCurrent->GetFirstMother()) iCurrent = pCurrent->GetFirstMother(); //Avoiding infinite loops
          if (iCurrent>=0) pCurrent = stack->Particle(iCurrent);
        } 
        if ((TMath::Abs(pCurrent->GetPdgCode())==211)&&(iCurrent>=0))
        {
          ePion = pCurrent->Energy();
          fhPionTrueE->Fill(ePion,ecluster,lambdaMainCluster);
          fhRatioEPion->Fill(ecluster/ePion);
        }
        fhMCPdg->Fill(ecluster,211);
        if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCPion: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",ePion,ecluster,lambdaMainCluster);
      }	
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCProton))//kMCProton
      {
        fhEMCProton  ->Fill(ecluster);
        fhPhiMCProton ->Fill(ecluster,phicluster);
        fhEtaMCProton ->Fill(ecluster,etacluster);
        fhLambdaMCProton ->Fill(ecluster,lambdaMainCluster,iNumCell);
        
        Int_t iCurrent = ph->GetLabel();
        Double_t eProton = 0;
        TParticle* pCurrent = stack->Particle(iCurrent);
        while ((TMath::Abs(pCurrent->GetPdgCode())!=2212)&&(iCurrent>=0))
        {
          if (iCurrent!=pCurrent->GetFirstMother()) iCurrent = pCurrent->GetFirstMother(); //Avoiding infinite loops
          if (iCurrent>=0) pCurrent = stack->Particle(iCurrent);
        } 
        if ((TMath::Abs(pCurrent->GetPdgCode())==2212)&&(iCurrent>=0))
        {
          eProton = pCurrent->Energy();
          fhProtonTrueE->Fill(eProton,ecluster,lambdaMainCluster);
          fhRatioEProton->Fill(ecluster/eProton);
        }
        fhMCPdg->Fill(ecluster,2212);
        if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCProton: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",eProton,ecluster,lambdaMainCluster);
      } 
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))//kMCAntiProton
      {
        fhEMCAntiProton  ->Fill(ecluster);
        fhPhiMCAntiProton ->Fill(ecluster,phicluster);
        fhEtaMCAntiProton ->Fill(ecluster,etacluster);
        fhLambdaMCAntiProton ->Fill(ecluster,lambdaMainCluster,iNumCell);
        
        Int_t iCurrent = ph->GetLabel();
        Double_t eAntiProton = 0;
        TParticle* pCurrent = stack->Particle(iCurrent);
        while ((TMath::Abs(pCurrent->GetPdgCode())!=-2212)&&(iCurrent>=0))
        {
          if (iCurrent!=pCurrent->GetFirstMother()) iCurrent = pCurrent->GetFirstMother(); //Avoiding infinite loops
          if (iCurrent>=0) pCurrent = stack->Particle(iCurrent);
        } 
        if ((TMath::Abs(pCurrent->GetPdgCode())==-2212)&&(iCurrent>=0))
        {
          eAntiProton = pCurrent->Energy();
          fhAntiProtonTrueE->Fill(eAntiProton,ecluster,lambdaMainCluster);
          fhRatioEAntiProton->Fill(ecluster/eAntiProton);
        }	
        fhMCPdg->Fill(ecluster,-2212);
        if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCAntiProton: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",eAntiProton,ecluster,lambdaMainCluster);
      }
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCNeutron))//kMCNeutron
      {
        fhEMCNeutron  ->Fill(ecluster);
        fhPhiMCNeutron ->Fill(ecluster,phicluster);
        fhEtaMCNeutron ->Fill(ecluster,etacluster);
        fhLambdaMCNeutron ->Fill(ecluster,lambdaMainCluster,iNumCell);
        
        Int_t iCurrent = ph->GetLabel();
        Double_t eNeutron = 0;
        TParticle* pCurrent = stack->Particle(iCurrent);
        while ((TMath::Abs(pCurrent->GetPdgCode())!=2112)&&(iCurrent>=0))
        {
          if (iCurrent!=pCurrent->GetFirstMother()) iCurrent = pCurrent->GetFirstMother(); //Avoiding infinite loops
          if (iCurrent>=0) pCurrent = stack->Particle(iCurrent);
        } 
        if ((TMath::Abs(pCurrent->GetPdgCode())==2112)&&(iCurrent>=0))
        {
          eNeutron = pCurrent->Energy();
          fhNeutronTrueE->Fill(eNeutron,ecluster,lambdaMainCluster);
          fhRatioENeutron->Fill(ecluster/eNeutron);
        }
        fhMCPdg->Fill(ecluster,2112);
	  	  if(GetDebug() > 3) 
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() - MC Cluster identified as kMCNeutron: True E: %3.2f, Reco E: %3.2f, Lambda0:  %3.2f\n",eNeutron,ecluster,lambdaMainCluster);
      }
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta))//kMCEta
      {
        fhEMCEta  ->Fill(ecluster);
        fhPhiMCEta ->Fill(ecluster,phicluster);
        fhEtaMCEta ->Fill(ecluster,etacluster);
        fhLambdaMCEta ->Fill(ecluster,lambdaMainCluster,iNumCell);
      }      
      
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      if(label < 0) {
        printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
        continue;
      }
      
      //Calculate Pi0 decay angles
      if(stack && (IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC)) ){
        for(Int_t i=0 ; i<stack->GetNprimary(); i++){
          TParticle * prim = stack->Particle(i) ;
          if(prim->GetPdgCode() == 111){
            TLorentzVector mom1 ;
            (stack->Particle((prim->GetFirstDaughter())))->Momentum(mom1);
            TLorentzVector mom2 ;
            (stack->Particle((prim->GetLastDaughter())))->Momentum(mom2);
            Float_t fDecayAngle = mom1.Angle(mom2.Vect());
            fhDecayAngle->Fill(ecluster,fDecayAngle);
          }
        }
      }
      
      
      Float_t eprim   = 0;
      Float_t ptprim  = 0;
      if(GetReader()->ReadStack()){
        
        if(label >=  stack->GetNtrack()) {
          if(GetDebug() > 2)  printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
          continue ;
        }
        
        primary = stack->Particle(label);
        if(!primary){
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        eprim   = primary->Energy();
        ptprim  = primary->Pt();
      } 
      else if(GetReader()->ReadAODMCParticles()){
        //Check which is the input
        if(ph->GetInputFileIndex() == 0){
          if(!mcparticles0) continue;
          if(label >=  mcparticles0->GetEntriesFast()) {
            if(GetDebug() > 2)  printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n",label, mcparticles0->GetEntriesFast());
            continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles0->At(label);
          
        }
        //        else {//Second input
        //          if(!mcparticles1) continue;
        //          if(label >=  mcparticles1->GetEntriesFast()) {
        //            if(GetDebug() > 2)  printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n",label, mcparticles1->GetEntriesFast());
        //            continue ;
        //          }
        //          //Get the particle
        //          aodprimary = (AliAODMCParticle*) mcparticles1->At(label); 
        //        }//second input
        
        if(!aodprimary){
          printf("AliAnaShowerParameter::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        
        eprim   = aodprimary->E();
        ptprim  = aodprimary->Pt();
      }
      //Write down MC truth information concerning all clusters
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
void AliAnaShowerParameter::Print(const Option_t * opt) const
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
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %f \n", fNCellsCut);
  printf("Number of clusters in envent is        > %d \n", fNumClusters);
  printf("    \n") ;
  
} 
