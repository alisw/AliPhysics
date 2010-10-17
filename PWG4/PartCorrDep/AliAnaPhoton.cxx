 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, providGetMixedEvent()ed that the above copyright notice appears in all   *
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
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>
#include "TParticle.h"

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
    AliAnaPartCorrBaseClass(), fCalorimeter(""), 
    fMinDist(0.),fMinDist2(0.),fMinDist3(0.),fRejectTrackMatch(0),
    fCheckConversion(kFALSE),fAddConvertedPairsToAOD(kFALSE), fMassCut(0),
    fTimeCutMin(-1), fTimeCutMax(9999999), fNCellsCut(0),
    fhPtPhoton(0),fhPhiPhoton(0),fhEtaPhoton(0),
    //MC
    fhDeltaE(0), fhDeltaPt(0),fhRatioE(0), fhRatioPt(0),fh2E(0),fh2Pt(0),
    fhPtMCPhoton(0),fhPhiMCPhoton(0),fhEtaMCPhoton(0), 
    fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0), 
    fhPtFragmentation(0),fhPhiFragmentation(0),fhEtaFragmentation(0), 
    fhPtISR(0),fhPhiISR(0),fhEtaISR(0), 
    fhPtPi0Decay(0),fhPhiPi0Decay(0),fhEtaPi0Decay(0), 
    fhPtOtherDecay(0),fhPhiOtherDecay(0),fhEtaOtherDecay(0), 
    fhPtConversion(0),fhPhiConversion(0),fhEtaConversion(0), 
    fhPtUnknown(0),fhPhiUnknown(0),fhEtaUnknown(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}//____________________________________________________________________________
AliAnaPhoton::~AliAnaPhoton() 
{
  //dtor

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
	
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  //Histograms of highest Photon identified in Event
  fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
  fhPtPhoton->SetYTitle("N");
  fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtPhoton) ; 
  
  fhPhiPhoton  = new TH2F
    ("hPhiPhoton","#phi_{#gamma}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiPhoton->SetYTitle("#phi");
  fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiPhoton) ; 
  
  fhEtaPhoton  = new TH2F
    ("hEtaPhoton","#eta_{#gamma}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaPhoton->SetYTitle("#eta");
  fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaPhoton) ;
  
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
   
    fhPtMCPhoton  = new TH1F("hPtMCPhoton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtMCPhoton->SetYTitle("N");
    fhPtMCPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtMCPhoton) ; 
    
    fhPhiMCPhoton  = new TH2F
      ("hPhiMCPhoton","#phi_{#gamma}, #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCPhoton->SetYTitle("#phi");
    fhPhiMCPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiMCPhoton) ; 
    
    fhEtaMCPhoton  = new TH2F
      ("hEtaMCPhoton","#eta_{#gamma}, #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCPhoton->SetYTitle("#eta");
    fhEtaMCPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaMCPhoton) ;
    
    fhPtPrompt  = new TH1F("hPtMCPrompt","Number of prompt #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtPrompt->SetYTitle("N");
    fhPtPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPrompt) ; 
    
    fhPhiPrompt  = new TH2F
      ("hPhiMCPrompt","#phi_{#gamma}, prompt #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiPrompt->SetYTitle("#phi");
    fhPhiPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPrompt) ; 
    
    fhEtaPrompt  = new TH2F
      ("hEtaMCPrompt","#eta_{#gamma}, prompt #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaPrompt->SetYTitle("#eta");
    fhEtaPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPrompt) ;
    
    fhPtFragmentation  = new TH1F("hPtMCFragmentation","Number of fragmentation #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtFragmentation->SetYTitle("N");
    fhPtFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtFragmentation) ; 
    
    fhPhiFragmentation  = new TH2F
      ("hPhiMCFragmentation","#phi_{#gamma}, fragmentation #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiFragmentation->SetYTitle("#phi");
    fhPhiFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiFragmentation) ; 
    
    fhEtaFragmentation  = new TH2F
      ("hEtaMCFragmentation","#eta_{#gamma}, fragmentation #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaFragmentation->SetYTitle("#eta");
    fhEtaFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaFragmentation) ;
    
    fhPtISR  = new TH1F("hPtMCISR","Number of initial state radiation #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtISR->SetYTitle("N");
    fhPtISR->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtISR) ; 
    
    fhPhiISR  = new TH2F
      ("hPhiMCISR","#phi_{#gamma} initial state radiation",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiISR->SetYTitle("#phi");
    fhPhiISR->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiISR) ; 
    
    fhEtaISR  = new TH2F
      ("hEtaMCISR","#eta_{#gamma} initial state radiation",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaISR->SetYTitle("#eta");
    fhEtaISR->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaISR) ;
    
    fhPtPi0Decay  = new TH1F("hPtMCPi0Decay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtPi0Decay->SetYTitle("N");
    fhPtPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPi0Decay) ; 
    
    fhPhiPi0Decay  = new TH2F
      ("hPhiMCPi0Decay","#phi_{#gamma}, #pi^{0} decay #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiPi0Decay->SetYTitle("#phi");
    fhPhiPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPi0Decay) ; 
    
    fhEtaPi0Decay  = new TH2F
      ("hEtaMCPi0Decay","#eta_{#gamma}, #pi^{0} #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaPi0Decay->SetYTitle("#eta");
    fhEtaPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPi0Decay) ;
    
    fhPtOtherDecay  = new TH1F("hPtMCOtherDecay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtOtherDecay->SetYTitle("N");
    fhPtOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtOtherDecay) ; 
    
    fhPhiOtherDecay  = new TH2F
      ("hPhiMCOtherDecay","#phi_{#gamma}, other decay #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiOtherDecay->SetYTitle("#phi");
    fhPhiOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiOtherDecay) ; 
    
    fhEtaOtherDecay  = new TH2F
      ("hEtaMCOtherDecay","#eta_{#gamma}, other decay #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaOtherDecay->SetYTitle("#eta");
    fhEtaOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaOtherDecay) ;
    
    fhPtConversion  = new TH1F("hPtMCConversion","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtConversion->SetYTitle("N");
    fhPtConversion->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtConversion) ; 
    
    fhPhiConversion  = new TH2F
      ("hPhiMCConversion","#phi_{#gamma}, conversion #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiConversion->SetYTitle("#phi");
    fhPhiConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiConversion) ; 
    
    fhEtaConversion  = new TH2F
      ("hEtaMCConversion","#eta_{#gamma}, conversion #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaConversion->SetYTitle("#eta");
    fhEtaConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaConversion) ;
    
    fhPtUnknown  = new TH1F("hPtMCUnknown","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtUnknown->SetYTitle("N");
    fhPtUnknown->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtUnknown) ; 
    
    fhPhiUnknown  = new TH2F
      ("hPhiMCUnknown","#phi_{#gamma}, unknown origin",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiUnknown->SetYTitle("#phi");
    fhPhiUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiUnknown) ; 
    
    fhEtaUnknown  = new TH2F
      ("hEtaMCUnknown","#eta_{#gamma}, unknown origin",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaUnknown->SetYTitle("#eta");
    fhEtaUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaUnknown) ;
	
  }//Histos with MC
    
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
	
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillAOD() 
{
  //Do photon analysis and fill aods
  
  //  Double_t vertex2[] = {0,0,0} ; //vertex from second input aod
  //Get the vertex and check it is not too large in z, cut for SE
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if(!GetMixedEvent() && TMath::Abs(v[2]) > GetZvertexCut()) return ;  
  
  //Select the Calorimeter of the photon
  TObjArray * pl = 0x0; 
  if(fCalorimeter == "PHOS")
    pl = GetAODPHOS();
  else if (fCalorimeter == "EMCAL")
    pl = GetAODEMCAL();
  
  if(!pl) {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }

  //Fill AODParticle with PHOS/EMCAL aods
  TLorentzVector mom, mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  Bool_t * indexConverted = new Bool_t[nCaloClusters];
  for (Int_t i = 0; i < nCaloClusters; i++) 
    indexConverted[i] = kFALSE;
	
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++){    
	  
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));	
    //printf("calo %d, %f\n",icalo,calo->E());
    
    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      //Get the vertex and check it is not too large in z
      if(TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut()) continue ;  
    }

    //Cluster selection, not charged, with photon id and in fiducial cut
	  
    //Input from second AOD?
    //Int_t input = 0;
    //    if (fCalorimeter == "EMCAL" && GetReader()->GetAODEMCALNormalInputEntries() <= icalo) 
    //      input = 1 ;
    //    else if(fCalorimeter == "PHOS"  && GetReader()->GetAODPHOSNormalInputEntries()  <= icalo) 
    //      input = 1;
	  
    //Get Momentum vector, 
    //if (input == 0) 
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;}//Assume that come from vertex in straight line
    else{
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
    //printf("AliAnaPhoton::MakeAnalysisFillAOD(): Vertex : %f,%f,%f\n",GetVertex(evtIndex)[0] ,GetVertex(evtIndex)[1],GetVertex(evtIndex)[2]);

    //    else if(input == 1) 
    //      calo->GetMomentum(mom,vertex2);//Assume that come from vertex in straight line  
    
    //If too small or big pt, skip it
    if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 
    Double_t tof = calo->GetTOF()*1e9;
    
    if(tof < fTimeCutMin || tof > fTimeCutMax) continue;
	  
    if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) continue;
	  
    //printf("AliAnaPhoton::Current Event %d; Current File Name : %s, E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",GetReader()->GetEventNumber(),(GetReader()->GetCurrentFileName()).Data(), 
      //     mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
    
    //Check acceptance selection
    if(IsFiducialCutOn()){
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue ;
    }
    //printf("Fiducial cut passed \n");
    
    //Create AOD for analysis
    AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
    Int_t label = calo->GetLabel();
    aodph.SetLabel(label);
    //aodph.SetInputFileIndex(input);
    
    //printf("Index %d, Id %d\n",icalo, calo->GetID());
    //Set the indeces of the original caloclusters  
    aodph.SetCaloLabel(calo->GetID(),-1);
    aodph.SetDetector(fCalorimeter);
    if(GetDebug() > 1) 
      printf("AliAnaPhoton::MakeAnalysisFillAOD() - Min pt cut and fiducial cut passed: pt %3.2f, phi %2.2f, eta %1.2f\n",aodph.Pt(),aodph.Phi(),aodph.Eta());	
    
    //Check Distance to Bad channel, set bit.
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad < fMinDist) //In bad channel (PHOS cristal size 2.2x2.2 cm)
      continue ;
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Bad channel cut passed %4.2f\n",distBad);
    
    if     (distBad > fMinDist3) aodph.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodph.SetDistToBad(1) ; 
    else                         aodph.SetDistToBad(0) ;
    
    //Skip matched clusters with tracks
    if(fRejectTrackMatch && IsTrackMatched(calo)) continue ;
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - TrackMatching cut passed \n");
    
    //Check PID
    //PID selection or bit setting
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetPdg());	 
      //If primary is not photon, skip it.
      if(aodph.GetPdg() != AliCaloPID::kPhoton) continue ;
    }					
    else if(IsCaloPIDOn()){
      
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn())
        aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      else
        aodph.SetPdg(GetCaloPID()->GetPdg(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetPdg());
      
      //If cluster does not pass pid, not photon, skip it.
      if(aodph.GetPdg() != AliCaloPID::kPhoton) continue ;			
      
    }
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph);
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetPdg());
    
    //Play with the MC stack if available
    //Check origin of the candidates
    if(IsDataMC()){
      
      aodph.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex()));
      if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
    
    
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
        AliVCluster * calo2 =  (AliVCluster*) (pl->At(jcalo));              //Get cluster kinematics
        Int_t evtIndex2 = 0 ; 
        if (GetMixedEvent()) {
          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo2->GetID()) ; 
          
        }        

        if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
          calo->GetMomentum(mom2,GetVertex(evtIndex2)) ;}//Assume that come from vertex in straight line
        else{
          Double_t vertex[]={0,0,0};
          calo->GetMomentum(mom2,vertex) ;
        }
        
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
          //aodpair.SetInputFileIndex(input);
          
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
    //printf("\t \t added single cluster %d\n",icalo);
	  
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  delete [] indexConverted;
	
  if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
	// Access MC information in stack if requested, check that it exists.	
	AliStack * stack = 0x0;
	TParticle * primary = 0x0;   
	TClonesArray * mcparticles0 = 0x0;
	//TClonesArray * mcparticles1 = 0x0;
	AliAODMCParticle * aodprimary = 0x0; 
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
			mcparticles0 = GetReader()->GetAODMCParticles(0);
			if(!mcparticles0 && GetDebug() > 0) 	{
				printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
			}	
      //			if(GetReader()->GetSecondInputAODTree()){
      //				mcparticles1 = GetReader()->GetAODMCParticles(1);
      //				if(!mcparticles1 && GetDebug() > 0) 	{
      //					printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
      //				}
      //			}		
			
		}
	}// is data and MC
	
	//Loop on stored AOD photons
	Int_t naod = GetOutputAODBranch()->GetEntriesFast();
	if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
	
	for(Int_t iaod = 0; iaod < naod ; iaod++){
	  AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
	  Int_t pdg = ph->GetPdg();
	  
	  if(GetDebug() > 3) 
	    printf("AliAnaPhoton::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ph->GetPdg(),ph->GetTag(), (ph->GetDetector()).Data()) ;
	  
	  //If PID used, fill histos with photons in Calorimeter fCalorimeter
	  if(IsCaloPIDOn() && pdg != AliCaloPID::kPhoton) continue; 
	  if(ph->GetDetector() != fCalorimeter) continue;
	  
	  if(GetDebug() > 2) 
	    printf("AliAnaPhoton::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
	  
	  //Fill photon histograms 
	  Float_t ptcluster  = ph->Pt();
	  Float_t phicluster = ph->Phi();
	  Float_t etacluster = ph->Eta();
	  Float_t ecluster   = ph->E();
	  
	  fhPtPhoton  ->Fill(ptcluster);
	  fhPhiPhoton ->Fill(ptcluster,phicluster);
	  fhEtaPhoton ->Fill(ptcluster,etacluster);
	  
	  //Play with the MC data if available
	  if(IsDataMC()){
	    
	    Int_t tag =ph->GetTag();
	    
	    if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
      {
        fhPtMCPhoton  ->Fill(ptcluster);
        fhPhiMCPhoton ->Fill(ptcluster,phicluster);
        fhEtaMCPhoton ->Fill(ptcluster,etacluster);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
        {
          fhPtConversion  ->Fill(ptcluster);
          fhPhiConversion ->Fill(ptcluster,phicluster);
          fhEtaConversion ->Fill(ptcluster,etacluster);
        }			
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)){
          fhPtPrompt  ->Fill(ptcluster);
          fhPhiPrompt ->Fill(ptcluster,phicluster);
          fhEtaPrompt ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
        {
          fhPtFragmentation  ->Fill(ptcluster);
          fhPhiFragmentation ->Fill(ptcluster,phicluster);
          fhEtaFragmentation ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
        {
          fhPtISR  ->Fill(ptcluster);
          fhPhiISR ->Fill(ptcluster,phicluster);
          fhEtaISR ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhPtPi0Decay  ->Fill(ptcluster);
          fhPhiPi0Decay ->Fill(ptcluster,phicluster);
          fhEtaPi0Decay ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhPtOtherDecay  ->Fill(ptcluster);
          fhPhiOtherDecay ->Fill(ptcluster,phicluster);
          fhEtaOtherDecay ->Fill(ptcluster,etacluster);
        }
      }
	    else{
	      fhPtUnknown  ->Fill(ptcluster);
	      fhPhiUnknown ->Fill(ptcluster,phicluster);
	      fhEtaUnknown ->Fill(ptcluster,etacluster);
	      
        //		 printf(" AliAnaPhoton::MakeAnalysisFillHistograms() - Label %d, pT %2.3f Unknown, bits set: ",
        //					ph->GetLabel(),ph->Pt());
        //		  for(Int_t i = 0; i < 20; i++) {
        //			  if(GetMCAnalysisUtils()->CheckTagBit(tag,i)) printf(" %d, ",i);
        //		  }
        //		  printf("\n");
        
	    }
	    
	    
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
          if(!mcparticles0) continue;
          if(label >=  mcparticles0->GetEntriesFast()) {
            if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                       label, mcparticles0->GetEntriesFast());
            continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles0->At(label);
          
	      }
//	      else {//Second input
//          if(!mcparticles1) continue;
//          if(label >=  mcparticles1->GetEntriesFast()) {
//            if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
//                                       label, mcparticles1->GetEntriesFast());
//            continue ;
//          }
//          //Get the particle
//          aodprimary = (AliAODMCParticle*) mcparticles1->At(label);
//          
//	      }//second input
	      
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
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
} 
