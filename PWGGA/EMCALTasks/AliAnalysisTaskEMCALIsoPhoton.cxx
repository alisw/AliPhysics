// $Id$

#include "AliAnalysisTaskEMCALIsoPhoton.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliV0vertexer.h"
#include "AliVCluster.h"

ClassImp(AliAnalysisTaskEMCALIsoPhoton)

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton() : 
  AliAnalysisTaskSE(), 
  fCaloClusters(0),
  fSelPrimTracks(0),
  fEMCalCells(0),
  fPrTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fESD(0),
  fOutputList(0),
  fEvtSel(0),
  fPVtxZ(0),                 
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),        
  fM02Et(0),                 
  fM02EtTM(0),               
  fM02EtCeIso1(0),           
  fM02EtCeIso2(0),           
  fM02EtCeIso5(0),           
  fM02EtTrIso1(0),           
  fM02EtTrIso2(0),           
  fM02EtTrIso5(0),           
  fM02EtAllIso1(0),          
  fM02EtAllIso2(0),          
  fM02EtAllIso5(0),          
  fCeIsoVsEtPho(0),          
  fTrIsoVsEtPho(0),          
  fAllIsoVsEtPho(0),         
  fCeIsoVsEtPi0(0),          
  fTrIsoVsEtPi0(0),          
  fAllIsoVsEtPi0(0),         
  fM02EtCeIso1TM(0),         
  fM02EtCeIso2TM(0),         
  fM02EtCeIso5TM(0),         
  fM02EtTrIso1TM(0),         
  fM02EtTrIso2TM(0),         
  fM02EtTrIso5TM(0),         
  fM02EtAllIso1TM(0),        
  fM02EtAllIso2TM(0),        
  fM02EtAllIso5TM(0),        
  fCeIsoVsEtPhoTM(0),        
  fTrIsoVsEtPhoTM(0),        
  fAllIsoVsEtPhoTM(0),    
  fCeIsoVsEtPi0TM(0),        
  fTrIsoVsEtPi0TM(0),        
  fAllIsoVsEtPi0TM(0)    
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton(const char *name) : 
  AliAnalysisTaskSE(name), 
  fCaloClusters(0),
  fSelPrimTracks(0),
  fEMCalCells(0),
  fPrTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fESD(0),
  fOutputList(0),
  fEvtSel(0),
  fPVtxZ(0),            
  fCellAbsIdVsAmpl(0),  
  fNClusHighClusE(0),   
  fM02Et(0),            
  fM02EtTM(0),          
  fM02EtCeIso1(0),      
  fM02EtCeIso2(0),      
  fM02EtCeIso5(0),      
  fM02EtTrIso1(0),      
  fM02EtTrIso2(0),      
  fM02EtTrIso5(0),      
  fM02EtAllIso1(0),     
  fM02EtAllIso2(0),     
  fM02EtAllIso5(0),     
  fCeIsoVsEtPho(0),     
  fTrIsoVsEtPho(0),     
  fAllIsoVsEtPho(0),    
  fCeIsoVsEtPi0(0),          
  fTrIsoVsEtPi0(0),          
  fAllIsoVsEtPi0(0),         
  fM02EtCeIso1TM(0),
  fM02EtCeIso2TM(0),        
  fM02EtCeIso5TM(0),        
  fM02EtTrIso1TM(0),        
  fM02EtTrIso2TM(0),        
  fM02EtTrIso5TM(0),        
  fM02EtAllIso1TM(0),       
  fM02EtAllIso2TM(0),       
  fM02EtAllIso5TM(0),       
  fCeIsoVsEtPhoTM(0),       
  fTrIsoVsEtPhoTM(0),       
  fAllIsoVsEtPhoTM(0),       
  fCeIsoVsEtPi0TM(0),        
  fTrIsoVsEtPi0TM(0),        
  fAllIsoVsEtPi0TM(0)    
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::UserCreateOutputObjects()
{
  // Create histograms, called once.
    
  fCaloClusters = new TRefArray();
  fSelPrimTracks = new TObjArray();

  
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
    
  
  fPVtxZ = new TH1F("hPVtxZ","primary vertex Z before cut;prim-vz(cm) ;",200,-100,100);
  fOutputList->Add(fPVtxZ);
  
  fCellAbsIdVsAmpl = new TH2F("hCellAbsIdVsAmpl","cell abs id vs cell amplitude (energy);E (GeV);ID",200,0,100,24*48*10,-0.5,24*48*10-0.5);
  fOutputList->Add(fPVtxZ);

  fNClusHighClusE = new TH2F("hNClusHighClusE","total number of clusters vs. highest clus energy in the event;E (GeV);NClus",200,0,100,301,-0.5,300.5);
  fOutputList->Add(fNClusHighClusE);

  fM02Et = new TH2F("fM02Et","M02 vs Et for all clusters;E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02Et->Sumw2();
  fOutputList->Add(fM02Et);

  fM02EtTM = new TH2F("fM02EtTM","M02 vs Et for all track-matched clusters;E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTM->Sumw2();
  fOutputList->Add(fM02EtTM);

  fM02EtCeIso1 = new TH2F("fM02EtCeIso1","M02 vs Et for all clusters (ISO_{EMC}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso1->Sumw2();
  fOutputList->Add(fM02EtCeIso1);

  fM02EtCeIso2 = new TH2F("fM02EtCeIso2","M02 vs Et for all clusters (ISO_{EMC}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso2->Sumw2();
  fOutputList->Add(fM02EtCeIso2);

  fM02EtCeIso5 = new TH2F("fM02EtCeIso5","M02 vs Et for all clusters (ISO_{EMC}<5GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso5->Sumw2();
  fOutputList->Add(fM02EtCeIso5);

  fM02EtTrIso1 = new TH2F("fM02EtTrIso1","M02 vs Et for all clusters (ISO_{TRK}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso1->Sumw2();
  fOutputList->Add(fM02EtTrIso1);

  fM02EtTrIso2 = new TH2F("fM02EtTrIso2","M02 vs Et for all clusters (ISO_{TRK}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso2->Sumw2();
  fOutputList->Add(fM02EtTrIso2);

  fM02EtTrIso5 = new TH2F("fM02EtTrIso5","M02 vs Et for all clusters (ISO_{TRK}<5GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso5->Sumw2();
  fOutputList->Add(fM02EtTrIso5);

  fM02EtAllIso1 = new TH2F("fM02EtAllIso1","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso1->Sumw2();
  fOutputList->Add(fM02EtAllIso1);

  fM02EtAllIso2 = new TH2F("fM02EtAllIso2","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso2->Sumw2();
  fOutputList->Add(fM02EtAllIso2);

  fM02EtAllIso5 = new TH2F("fM02EtAllIso5","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso5->Sumw2();
  fOutputList->Add(fM02EtAllIso5);

  fCeIsoVsEtPho = new TH2F("fCeIsoVsEtPho","ISO_{EMC} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{EMC}",1000,0,100,1000,-10,190);
  fCeIsoVsEtPho->Sumw2();
  fOutputList->Add(fCeIsoVsEtPho);

  fTrIsoVsEtPho = new TH2F("fTrIsoVsEtPho","ISO_{TRK} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{TRK}",1000,0,100,1000,-10,190);
  fTrIsoVsEtPho->Sumw2();
  fOutputList->Add(fTrIsoVsEtPho);

  fAllIsoVsEtPho = new TH2F("fAllIsoVsEtPho","ISO_{EMC+TRK} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{EMC+TRK}",1000,0,100,1000,-10,190);
  fAllIsoVsEtPho->Sumw2();
  fOutputList->Add(fAllIsoVsEtPho);

  fCeIsoVsEtPi0 = new TH2F("fCeIsoVsEtPi0","ISO_{EMC} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{EMC}",1000,0,100,1000,-10,190);
  fCeIsoVsEtPi0->Sumw2();
  fOutputList->Add(fCeIsoVsEtPi0);

  fTrIsoVsEtPi0 = new TH2F("fTrIsoVsEtPi0","ISO_{TRK} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{TRK}",1000,0,100,1000,-10,190);
  fTrIsoVsEtPi0->Sumw2();
  fOutputList->Add(fTrIsoVsEtPi0);

  fAllIsoVsEtPi0 = new TH2F("fAllIsoVsEtPi0","ISO_{EMC+TRK} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{EMC+TRK}",1000,0,100,1000,-10,190);
  fAllIsoVsEtPi0->Sumw2();
  fOutputList->Add(fAllIsoVsEtPi0);

  fM02EtCeIso1TM = new TH2F("fM02EtCeIso1TM","M02 vs Et for all clusters (ISO_{EMC}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso1TM->Sumw2();
  fOutputList->Add(fM02EtCeIso1TM);

  fM02EtCeIso2TM = new TH2F("fM02EtCeIso2TM","M02 vs Et for all clusters (ISO_{EMC}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso2TM->Sumw2();
  fOutputList->Add(fM02EtCeIso2TM);

  fM02EtCeIso5TM = new TH2F("fM02EtCeIso5TM","M02 vs Et for all clusters (ISO_{EMC}<5GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtCeIso5TM->Sumw2();
  fOutputList->Add(fM02EtCeIso5TM);

  fM02EtTrIso1TM = new TH2F("fM02EtTrIso1TM","M02 vs Et for all clusters (ISO_{TRK}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso1TM->Sumw2();
  fOutputList->Add(fM02EtTrIso1TM);

  fM02EtTrIso2TM = new TH2F("fM02EtTrIso2TM","M02 vs Et for all clusters (ISO_{TRK}<1GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso2TM->Sumw2();
  fOutputList->Add(fM02EtTrIso2TM);

  fM02EtTrIso5TM = new TH2F("fM02EtTrIso5TM","M02 vs Et for all clusters (ISO_{TRK}<5GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtTrIso5TM->Sumw2();
  fOutputList->Add(fM02EtTrIso5TM);

  fM02EtAllIso1TM = new TH2F("fM02EtAllIso1TM","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso1TM->Sumw2();
  fOutputList->Add(fM02EtAllIso1TM);

  fM02EtAllIso2TM = new TH2F("fM02EtAllIso2TM","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso2TM->Sumw2();
  fOutputList->Add(fM02EtAllIso2TM);

  fM02EtAllIso5TM = new TH2F("fM02EtAllIso5TM","M02 vs Et for all clusters (ISO_{EMC+TRK}<2GeV);E_{T} (GeV);M02",1000,0,100,400,0,4);
  fM02EtAllIso5TM->Sumw2();
  fOutputList->Add(fM02EtAllIso5TM);

  fCeIsoVsEtPhoTM = new TH2F("fCeIsoVsEtPhoTM","ISO_{EMC} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{EMC}",1000,0,100,1000,-10,190);
  fCeIsoVsEtPhoTM->Sumw2();
  fOutputList->Add(fCeIsoVsEtPhoTM);

  fTrIsoVsEtPhoTM = new TH2F("fTrIsoVsEtPhoTM","ISO_{TRK} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{TRK}",1000,0,100,1000,-10,190);
  fTrIsoVsEtPhoTM->Sumw2();
  fOutputList->Add(fTrIsoVsEtPhoTM);

  fAllIsoVsEtPhoTM = new TH2F("fAllIsoVsEtPhoTM","ISO_{EMC+TRK} vs. E_{T}^{clus} (0.1<#lambda_{0}^{2}<0.3);E_{T} (GeV);ISO_{EMC+TRK}",1000,0,100,1000,-10,190);
  fAllIsoVsEtPhoTM->Sumw2();
  fOutputList->Add(fAllIsoVsEtPhoTM);

  fCeIsoVsEtPi0TM = new TH2F("fCeIsoVsEtPi0TM","ISO_{EMC} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{EMC}",1000,0,100,1000,-10,190);
  fCeIsoVsEtPi0TM->Sumw2();
  fOutputList->Add(fCeIsoVsEtPi0TM);

  fTrIsoVsEtPi0TM = new TH2F("fTrIsoVsEtPi0TM","ISO_{TRK} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{TRK}",1000,0,100,1000,-10,190);
  fTrIsoVsEtPi0TM->Sumw2();
  fOutputList->Add(fTrIsoVsEtPi0TM);

  fAllIsoVsEtPi0TM = new TH2F("fAllIsoVsEtPi0TM","ISO_{EMC+TRK} vs. E_{T}^{clus} (#pi^{0} selection for BG);E_{T} (GeV);ISO_{EMC+TRK}",1000,0,100,1000,-10,190);
  fAllIsoVsEtPi0TM->Sumw2();
  fOutputList->Add(fAllIsoVsEtPi0TM);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // event trigger selection
  Bool_t isSelected = 0;
  if(fPeriod.Contains("11a")){
    if(fTrigBit.Contains("kEMC"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
    else
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  }
  else{
    if(fTrigBit.Contains("kEMC"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7);
    else
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  }
  if(!isSelected )
        return; 

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  fEvtSel->Fill(0);
  
  AliESDVertex *pv = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(!pv) 
    return;
  fPVtxZ->Fill(pv->GetZ());
  if(TMath::Abs(pv->GetZ())>15)
    return;

  fEvtSel->Fill(1);

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(iTracks);
    if (!track)
      continue;
    if (fPrTrCuts && fPrTrCuts->IsSelected(track)){
      fSelPrimTracks->Add(track);
      //printf("pt,eta,phi:%1.1f,%1.1f,%1.1f \n",track->Pt(),track->Eta(), track->Phi());
    }
  }

  if(!fIsTrain){
    for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
      if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
        break;
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
    }
  }
  fESD->GetEMCALClusters(fCaloClusters);
  fEMCalCells = fESD->GetEMCALCells();
  
    
  FillClusHists(); 

  fCaloClusters->Clear();
  fSelPrimTracks->Clear();
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FillClusHists()
{
  // Fill cluster histograms.

  if(!fCaloClusters)
    return;
  const Int_t nclus = fCaloClusters->GetEntries();
  if(nclus==0)
    return;
  Double_t maxE = 0;
  for(Int_t ic=0;ic<nclus;ic++){
    maxE=0;
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(fCaloClusters->At(ic));
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    Short_t id;
    Double_t Emax = GetMaxCellEnergy( c, id);
    Double_t Ecross = GetCrossEnergy( c, id);
    if((1-Ecross/Emax)>fExoticCut)
      continue;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
    Float_t ceiso, cephiband, cecore;
    Float_t triso, trphiband, trcore;
    Float_t alliso, allphiband, allcore;
    Float_t phibandArea = (1.4 - 2*fIsoConeR)*2*fIsoConeR;
    Float_t netConeArea = TMath::Pi()*(fIsoConeR*fIsoConeR - 0.04*0.04);
    GetCeIso(clsVec, ceiso, cephiband, cecore);
    GetTrIso(clsVec, triso, trphiband, trcore);
    alliso = ceiso + triso;
    allphiband = cephiband + trphiband;
    allcore = cecore + trcore;
    Float_t ceisoue =  cephiband/phibandArea*netConeArea;
    Float_t trisoue =  trphiband/phibandArea*netConeArea;
    Float_t allisoue =  allphiband/phibandArea*netConeArea;
    Double_t dR = TMath::Sqrt(pow(c->GetTrackDx(),2)+pow(c->GetTrackDz(),2));
    Double_t M02u;
    if(Et<12)
      M02u = 0.02486*Et*Et - 0.7289*Et + 6.266;
    else
      M02u = 14.32/Et - 0.09863;
    if(M02u<0.65)
      M02u = 0.65;
    Double_t M02l = 12.88/Et - 0.3194;
    if(M02l<0.4)
      M02l = 0.4;
    if(dR>0.028){
      fM02Et->Fill(Et, c->GetM02());
      if(ceiso-cecore<1)
	fM02EtCeIso1->Fill(Et, c->GetM02());
      if(ceiso-cecore<2)
	fM02EtCeIso2->Fill(Et, c->GetM02());
      if(ceiso-cecore<5)
	fM02EtCeIso5->Fill(Et, c->GetM02());
      if(triso-trcore<1)
	fM02EtTrIso1->Fill(Et, c->GetM02());
      if(triso-trcore<2)
	fM02EtTrIso2->Fill(Et, c->GetM02());
      if(triso-trcore<5)
	fM02EtTrIso5->Fill(Et, c->GetM02());
      if(alliso-allcore<1)
	fM02EtAllIso1->Fill(Et, c->GetM02());
      if(alliso-allcore<2)
	fM02EtAllIso2->Fill(Et, c->GetM02());
      if(alliso-allcore<5)
	fM02EtAllIso5->Fill(Et, c->GetM02());
      if(c->GetM02()>0.1 && c->GetM02()<0.3){
	fCeIsoVsEtPho->Fill(Et, ceiso - cecore - ceisoue);
	fTrIsoVsEtPho->Fill(Et, triso - trcore - trisoue);
	fAllIsoVsEtPho->Fill(Et, alliso - allcore - allisoue);
      }
      if(c->GetM02()>M02l && c->GetM02()<M02u){
	fCeIsoVsEtPi0->Fill(Et, ceiso - cecore - ceisoue);
	fTrIsoVsEtPi0->Fill(Et, triso - trcore - trisoue);
	fAllIsoVsEtPi0->Fill(Et, alliso - allcore - allisoue);
      }
    }
    else
      {
      fM02EtTM->Fill(Et, c->GetM02());
      if(ceiso-cecore<1)
	fM02EtCeIso1TM->Fill(Et, c->GetM02());
      if(ceiso-cecore<2)
	fM02EtCeIso2TM->Fill(Et, c->GetM02());
      if(ceiso-cecore<5)
	fM02EtCeIso5TM->Fill(Et, c->GetM02());
      if(triso-trcore<1)
	fM02EtTrIso1TM->Fill(Et, c->GetM02());
      if(triso-trcore<2)
	fM02EtTrIso2TM->Fill(Et, c->GetM02());
      if(triso-trcore<5)
	fM02EtTrIso5TM->Fill(Et, c->GetM02());
      if(alliso-allcore<1)
	fM02EtAllIso1TM->Fill(Et, c->GetM02());
      if(alliso-allcore<2)
	fM02EtAllIso2TM->Fill(Et, c->GetM02());
      if(alliso-allcore<5)
	fM02EtAllIso5TM->Fill(Et, c->GetM02());
      if(c->GetM02()>0.1 && c->GetM02()<0.3){
	fCeIsoVsEtPhoTM->Fill(Et, ceiso - cecore - ceisoue);
	fTrIsoVsEtPhoTM->Fill(Et, triso - trcore - trisoue);
	fAllIsoVsEtPhoTM->Fill(Et, alliso - allcore - allisoue);
      }
      if(c->GetM02()>M02l && c->GetM02()<M02u){
	fCeIsoVsEtPi0TM->Fill(Et, ceiso - cecore - ceisoue);
	fTrIsoVsEtPi0TM->Fill(Et, triso - trcore - trisoue);
	fAllIsoVsEtPi0TM->Fill(Et, alliso - allcore - allisoue);
      }
    }
    if(c->E()>maxE)
      maxE = c->E();
  }
  fNClusHighClusE->Fill(maxE,nclus);
} 

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetCeIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core)
{
  // Get cell isolation.

  if(!fEMCalCells)
    return;
  const Int_t ncells = fEMCalCells->GetNumberOfCells();
  Float_t totiso=0;
  Float_t totband=0;
  Float_t totcore=0;
  Float_t etacl = vec.Eta();
  Float_t phicl = vec.Phi();
  Float_t thetacl = vec.Theta();
  if(phicl<0)
    phicl+=TMath::TwoPi();
  Int_t absid = -1;
  Float_t eta=-1, phi=-1;  
  for(int icell=0;icell<ncells;icell++){
    absid = TMath::Abs(fEMCalCells->GetCellNumber(icell));
    if(!fGeom)
      return;
    fGeom->EtaPhiFromIndex(absid,eta,phi);
    Float_t dphi = TMath::Abs(phi-phicl);
    Float_t deta = TMath::Abs(eta-etacl);
    Float_t R = TMath::Sqrt(deta*deta + dphi*dphi);
    Float_t etcell = fEMCalCells->GetCellAmplitude(absid)*TMath::Sin(thetacl);
    if(R<fIsoConeR){
      totiso += etcell;
      if(R<0.04)
	totcore += etcell;
    }
    else{
      if(dphi>fIsoConeR)
	continue;
      if(deta<fIsoConeR)
	continue;
      totband += fEMCalCells->GetCellAmplitude(absid)*TMath::Sin(thetacl);
    }
  }
  iso = totiso;
  phiband = totband;
  core = totcore;
} 

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetTrIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core)
{
  // Get track isolation.

  if(!fSelPrimTracks)
    return;
  const Int_t ntracks = fSelPrimTracks->GetEntries();
  Double_t totiso=0;
  Double_t totband=0;
  Double_t totcore=0;
  Double_t etacl = vec.Eta();
  Double_t phicl = vec.Phi();
  if(phicl<0)
    phicl+=TMath::TwoPi();
  for(int itrack=0;itrack<ntracks;itrack++){
    AliESDtrack *track = static_cast<AliESDtrack*> (fSelPrimTracks->At(itrack));
    if(!track)
      continue;
    Double_t dphi = TMath::Abs(track->Phi()-phicl);
    Double_t deta = TMath::Abs(track->Eta()-etacl);
    Double_t R = TMath::Sqrt(deta*deta + dphi*dphi);
    Double_t pt = track->Pt();
    if(R<fIsoConeR){
      totiso += track->Pt();
      if(R<0.04)
	totcore += pt;
    }
    else{
      if(dphi>fIsoConeR)
	continue;
      if(deta<fIsoConeR)
	continue;
      totband += track->Pt();
    }
  }
  iso = totiso;
  phiband = totband;
  core = totcore;
} 

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = 0;
  cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  if (!fGeom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
	(aphidiff==0 && aetadiff==1) ) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = 0;
  cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
