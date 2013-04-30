// $Id$

#include "AliAnalysisTaskEMCALIsoPhoton.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TList.h>

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
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliV0vertexer.h"
#include "AliVCluster.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALIsoPhoton)

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton() : 
  AliAnalysisTaskSE(), 
  fCaloClusters(0),
  fSelPrimTracks(0),
  fTracks(0),
  fEMCalCells(0),
  fPrTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fIsMc(0),
  fDebug(0),
  fPathStrOpt("/"),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
  fHigherPtCone(0),
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fRecoPV(0),
  fPVtxZ(0),
  fTrMultDist(0),
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),
  fHigherPtConeM02(0),
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fMcPtInConeBG(0),
  fMcPtInConeSBG(0),
  fMcPtInConeBGnoUE(0),
  fMcPtInConeSBGnoUE(0),
  fAllIsoEtMcGamma(0),
  fAllIsoNoUeEtMcGamma(0),
  fMCDirPhotonPtEtaPhiNoClus(0),
  fHnOutput(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton(const char *name) : 
  AliAnalysisTaskSE(name), 
  fCaloClusters(0),
  fSelPrimTracks(0),
  fTracks(0),
  fEMCalCells(0),
  fPrTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fIsMc(0),
  fDebug(0),
  fPathStrOpt("/"),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
  fHigherPtCone(0),
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fRecoPV(0),
  fPVtxZ(0),            
  fTrMultDist(0),
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),   
  fHigherPtConeM02(0),
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fMcPtInConeBG(0),
  fMcPtInConeSBG(0),
  fMcPtInConeBGnoUE(0),
  fMcPtInConeSBGnoUE(0),
  fAllIsoEtMcGamma(0),
  fAllIsoNoUeEtMcGamma(0),
  fMCDirPhotonPtEtaPhiNoClus(0),
  fHnOutput(0)
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
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName.Data());
  
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
  
  fNClusEt10 = new TH1F("hNClusEt10","# of cluster with E_{T}>10 per event;E;",101,-0.5,100.5);
  fOutputList->Add(fNClusEt10);
  
  fRecoPV = new TH1F("hRecoPV","Prim. vert. reconstruction (yes or no);reco (0=no, 1=yes) ;",2,-0.5,1.5);
  fOutputList->Add(fRecoPV);

  fPVtxZ = new TH1F("hPVtxZ","primary vertex Z before cut;prim-vz(cm) ;",200,-100,100);
  fOutputList->Add(fPVtxZ);

  fTrMultDist = new TH1F("fTrMultDist","track multiplicity;tracks/event;#events",200,0.5,200.5);
  fOutputList->Add(fTrMultDist);

  fMCDirPhotonPtEtaPhi = new TH3F("hMCDirPhotonPtEtaPhi","photon (gq->#gammaq) p_{T}, #eta, #phi;GeV/c;#eta;#phi",1000,0,100,154,-0.77,0.77,130,1.38,3.20);
  fMCDirPhotonPtEtaPhi->Sumw2();
  fOutputList->Add(fMCDirPhotonPtEtaPhi);

  fDecayPhotonPtMC = new TH1F("hDecayPhotonPtMC","decay photon p_{T};GeV/c;dN/dp_{T} (c GeV^{-1})",1000,0,100);
  fDecayPhotonPtMC->Sumw2();
  fOutputList->Add(fDecayPhotonPtMC);

  fCellAbsIdVsAmpl = new TH2F("hCellAbsIdVsAmpl","cell abs id vs cell amplitude (energy);E (GeV);ID",200,0,100,24*48*10,-0.5,24*48*10-0.5);
  fOutputList->Add(fCellAbsIdVsAmpl);

  fNClusHighClusE = new TH2F("hNClusHighClusE","total number of clusters vs. highest clus energy in the event;E (GeV);NClus",200,0,100,301,-0.5,300.5);
  fOutputList->Add(fNClusHighClusE);

  fHigherPtConeM02 = new TH2F("hHigherPtConeM02","#lambda_{0}^{2} vs. in-cone-p_{T}^{max};p_{T}^{max} (GeV/c, in the cone);#lambda_{0}^{2}",1000,0,100,400,0,4);
  fOutputList->Add(fHigherPtConeM02);

  fClusEtMcPt = new TH2F("hClusEtMcPt","E_{T}^{clus} vs. p_{T}^{mc}; p_{T}^{mc};E_{T}^{clus}",500,0,100,500,0,100);
  fOutputList->Add(fClusEtMcPt);

  fClusMcDetaDphi = new TH2F("hClusMcDetaDphi","#Delta#phi vs. #Delta#eta (reco-mc);#Delta#eta;#Delta#phi",100,-.7,.7,100,-.7,.7);
  fOutputList->Add(fClusMcDetaDphi);

  fNClusPerPho = new TH2F("hNClusPerPho","Number of clusters per prompt photon;p_{T}^{MC};N_{clus}",500,0,100,11,-0.5,10.5);
  fOutputList->Add(fNClusPerPho);

  fMcPtInConeBG = new TH2F("hMcPtInConeBG","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (BG template);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeBG);

  fMcPtInConeSBG  = new TH2F("hMcPtInConeSBG","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (SBG range);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeSBG);

  fMcPtInConeBGnoUE = new TH2F("hMcPtInConeBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (BG template);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeBGnoUE);

  fMcPtInConeSBGnoUE  = new TH2F("hMcPtInConeSBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (SBG range);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeSBGnoUE);

  fAllIsoEtMcGamma  = new TH2F("hAllIsoEtMcGammaE","ISO^{TRK+EMC} vs. E_{T}^{clus} for clusters comming from a MC prompt #gamma; E_{T}^{clus} (GeV);ISO^{TRK+EMC} (GeV);",1000,0,100,600,-10,50);
  fOutputList->Add(fAllIsoEtMcGamma);

  fAllIsoNoUeEtMcGamma  = new TH2F("hAllIsoNoUeEtMcGammaE","ISO^{TRK+EMC}_{noue} vs. E_{T}^{clus} for clusters comming from a MC prompt #gamma; E_{T}^{clus} (GeV);ISO^{TRK+EMC}_{noue} (GeV);",1000,0,100,600,-10,50);
  fOutputList->Add(fAllIsoNoUeEtMcGamma);


  fMCDirPhotonPtEtaPhiNoClus = new TH3F("hMCDirPhotonPhiEtaNoClus","p_{T}, #eta and  #phi of prompt photons with no reco clusters;p_{T};#eta;#phi",1000,0,100,154,-0.77,0.77,130,1.38,3.20);
  fOutputList->Add(fMCDirPhotonPtEtaPhiNoClus);

  Int_t nEt=1000, nM02=400, nCeIso=1000, nTrIso=1000,  nAllIso=1000,  nCeIsoNoUE=1000,  nAllIsoNoUE=1000, nTrClDphi=200, nTrClDeta=100, nClEta=140, nClPhi=128, nTime=60, nMult=100, nPhoMcPt=101;
  Int_t bins[] = {nEt, nM02, nCeIso, nTrIso, nAllIso, nCeIsoNoUE, nAllIsoNoUE, nTrClDphi, nTrClDeta,nClEta,nClPhi,nTime,nMult,nPhoMcPt};
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  Double_t xmin[] = { 0.,   0.,  -10.,   -10., -10., -10., -10., -0.1,-0.05, -0.7, 1.4,-0.15e-06,-0.5,-1.5};
  Double_t xmax[] = { 100., 4., 190., 190., 190.,  190., 190., 0.1, 0.05, 0.7, 3.192, 0.15e-06,99.5,99.5};
  fHnOutput =  new THnSparseF("fHnOutput","Output matrix: E_{T},M02,CeIso,TrIso,AllIso, CeIsoNoUESub, AllIsoNoUESub, d#phi_{trk},d#eta_{trk},#eta_{clus},#phi_{clus},T_{max},mult,mc-p_{T}^{#gamma}", ndims, bins, xmin, xmax);
  fOutputList->Add(fHnOutput);



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
      if(fTrigBit.Contains("kMB"))
	isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      else
	isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  }
  if(fIsMc)
    isSelected = kTRUE;
  if(fDebug)
    printf("isSelected = %d, fIsMC=%d\n", isSelected, fIsMc);
  if(!isSelected )
        return; 
  if(fIsMc){
    TTree *tree = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetTree();
    TFile *file = (TFile*)tree->GetCurrentFile();
    TString filename = file->GetName();
    if(!filename.Contains(fPathStrOpt.Data()))
      return;
  }
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  fEvtSel->Fill(0);
  if(fDebug)
    printf("fESD is ok\n");
  
  AliESDVertex *pv = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(!pv)
    return;
  if(!pv->GetStatus())
    fRecoPV->Fill(0);
  else
    fRecoPV->Fill(1);
  fPVtxZ->Fill(pv->GetZ());
  if(TMath::Abs(pv->GetZ())>15)
    return;
  if(fDebug)
    printf("passed vertex cut\n");

  fEvtSel->Fill(1);

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
  if(!fTracks){
    AliError("Track array in event is NULL!");
    return;
  }
  // Track loop to fill a pT spectrum
  const Int_t Ntracks = fTracks->GetEntriesFast();
  for (Int_t iTracks = 0;  iTracks < Ntracks; ++iTracks) {
    //  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    //AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(iTracks);
    AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(iTracks));
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
  AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  fTrackMult = fTrackCuts->GetReferenceMultiplicity(fESD);//kTrackletsITSTPC ,0.5); 
  fTrMultDist->Fill(fTrackMult);

  fESD->GetEMCALClusters(fCaloClusters);
  fEMCalCells = fESD->GetEMCALCells();
  
    

  fMCEvent = MCEvent();
  if(fMCEvent){
    if(fDebug)
      std::cout<<"MCevent exists\n";
    fStack = (AliStack*)fMCEvent->Stack();
  }
  else{
    if(fDebug)
      std::cout<<"ERROR: NO MC EVENT!!!!!!\n";
  }
  FollowGamma();
  if(fDebug)
    printf("passed calling of FollowGamma\n");
  FillClusHists(); 
  if(fDebug)
    printf("passed calling of FillClusHists\n");
  FillMcHists();
  if(fDebug)
    printf("passed calling of FillMcHists\n");

  fCaloClusters->Clear();
  fSelPrimTracks->Clear();
  fNClusForDirPho = 0;

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
  if(fDebug)
    printf("Inside FillClusHists and there are %d clusters\n",nclus);
  Double_t maxE = 0;
  Int_t nclus10 = 0;
  Double_t ptmc=-1;
  for(Int_t ic=0;ic<nclus;ic++){
    maxE=0;
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(fCaloClusters->At(ic));
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    if(c->E()<fECut)
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
    if(fDebug)
      printf("\tcluster eta=%1.1f,phi=%1.1f,E=%1.1f\n",clsVec.Eta(),clsVec.Phi(),c->E());
    if(Et>10)
      nclus10++;
    Float_t ceiso, cephiband, cecore;
    Float_t triso, trphiband, trcore;
    Float_t alliso, allphiband, allcore;
    Float_t phibandArea = (1.4 - 2*fIsoConeR)*2*fIsoConeR;
    Float_t netConeArea = TMath::Pi()*(fIsoConeR*fIsoConeR - 0.04*0.04);
    GetCeIso(clsVec, id, ceiso, cephiband, cecore);
    GetTrIso(clsVec, triso, trphiband, trcore);
    Double_t dr = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx() + c->GetTrackDz()*c->GetTrackDz());
    if(Et>10 && Et<15 && dr>0.025){
      fHigherPtConeM02->Fill(fHigherPtCone,c->GetM02());
      if(fDebug)
	printf("\t\tHigher track pt inside the cone: %1.1f GeV/c; M02=%1.2f\n",fHigherPtCone,c->GetM02());
    }
    alliso = ceiso + triso;
    allphiband = cephiband + trphiband;
    allcore = cecore + trcore;
    Float_t ceisoue =  cephiband/phibandArea*netConeArea;
    Float_t trisoue =  trphiband/phibandArea*netConeArea;
    Float_t allisoue =  allphiband/phibandArea*netConeArea;
    Float_t mcptsum = GetMcPtSumInCone(clsVec.Eta(), clsVec.Phi(),fIsoConeR); 
    if(fDebug && Et>10)
      printf("\t alliso=%1.1f, Et=%1.1f=-=-=-=-=\n",alliso,Et);
    if(c->GetM02()>0.5 && c->GetM02()<2.0){
      fMcPtInConeBG->Fill(alliso-Et-allisoue, mcptsum);
      fMcPtInConeBGnoUE->Fill(alliso-Et, mcptsum);
    }
    Double_t r = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx()+c->GetTrackDz()*c->GetTrackDz());
    if(c->GetM02()>0.1 && c->GetM02()<0.3 && r>0.03){
      fMcPtInConeSBG->Fill(alliso-Et-allisoue, mcptsum);
      fMcPtInConeSBGnoUE->Fill(alliso-Et, mcptsum);
      if(fMcIdFamily.Contains((Form("%d",c->GetLabel())))){
	fAllIsoEtMcGamma->Fill(Et, alliso-/*Et*/cecore-allisoue);
	fAllIsoNoUeEtMcGamma->Fill(Et, alliso-/*Et*/cecore);
      }
    }
    const Int_t ndims =   fNDimensions;
    Double_t outputValues[ndims];
    ptmc = GetClusSource(c);
    outputValues[0] = Et;
    outputValues[1] = c->GetM02();
    outputValues[2] = ceiso-Et/*cecore*/-ceisoue;
    outputValues[3] = triso-trisoue;
    outputValues[4] = alliso-Et/*cecore*/-allisoue;
    outputValues[5] = ceiso-Et;
    outputValues[6] = alliso-Et;
    outputValues[7] = c->GetTrackDx();
    outputValues[8] = c->GetTrackDz();
    outputValues[9] = clsVec.Eta();
    outputValues[10] = clsVec.Phi();
    outputValues[11] = fEMCalCells->GetCellTime(id);
    outputValues[12] = fTrackMult;
    outputValues[13] = ptmc;
    fHnOutput->Fill(outputValues);
    if(c->E()>maxE)
      maxE = c->E();
  }
  fNClusHighClusE->Fill(maxE,nclus);
  fNClusEt10->Fill(nclus10);
  fNClusPerPho->Fill(fDirPhoPt,fNClusForDirPho);
} 

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetCeIso(TVector3 vec, Int_t maxid, Float_t &iso, Float_t &phiband, Float_t &core)
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
  Float_t maxtcl = fEMCalCells->GetCellTime(maxid);
  if(phicl<0)
    phicl+=TMath::TwoPi();
  Int_t absid = -1;
  Float_t eta=-1, phi=-1;  
  for(int icell=0;icell<ncells;icell++){
    absid = TMath::Abs(fEMCalCells->GetCellNumber(icell));
    Float_t celltime = fEMCalCells->GetCellTime(absid);
    //if(TMath::Abs(celltime)>2e-8 && (!fIsMc))
    if(TMath::Abs(celltime-maxtcl)>2e-8 )
      continue;
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
  fHigherPtCone = 0;
  const Int_t ntracks = fSelPrimTracks->GetEntries();
  Double_t totiso=0;
  Double_t totband=0;
  Double_t totcore=0;
  Double_t etacl = vec.Eta();
  Double_t phicl = vec.Phi();
  if(phicl<0)
    phicl+=TMath::TwoPi();
  for(int itrack=0;itrack<ntracks;itrack++){
    AliVTrack *track = static_cast<AliVTrack*> (fSelPrimTracks->At(itrack));
    if(!track)
      continue;
    Double_t dphi = TMath::Abs(track->Phi()-phicl);
    Double_t deta = TMath::Abs(track->Eta()-etacl);
    Double_t R = TMath::Sqrt(deta*deta + dphi*dphi);
    Double_t pt = track->Pt();
    if(pt>fHigherPtCone)
      fHigherPtCone = pt;
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
void AliAnalysisTaskEMCALIsoPhoton ::FillMcHists()
{
  if(!fStack)
    return;
  Int_t nTracks = fStack->GetNtrack();
  if(fDebug)
    printf("Inside FillMcHists and there are %d mcparts\n",nTracks);
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++){
    TParticle *mcp = static_cast<TParticle*>(fStack->Particle(iTrack));  
    if(!mcp)
      continue;  
    Int_t pdg = mcp->GetPdgCode();
    if(pdg!=22)
      continue;
    if(TMath::Abs(mcp->Eta())>0.7 ||mcp->Phi()<1.4 || mcp->Phi()>3.2)
      continue;
    Int_t imom = mcp->GetMother(0);
    if(imom<0 || imom>nTracks)
      continue;
    TParticle *mcmom = static_cast<TParticle*>(fStack->Particle(imom));  
    if(!mcmom)
      continue;
    Int_t pdgMom = mcmom->GetPdgCode();
    if((imom==6 || imom==7) && pdgMom==22) {
      fMCDirPhotonPtEtaPhi->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
      if(fNClusForDirPho==0)
	fMCDirPhotonPtEtaPhiNoClus->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
      if(fDebug){
	printf("Found \"photonic\" parton at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d,\n",imom,mcmom->Pt(), mcmom->Eta(), mcmom->Phi(), mcmom->GetStatusCode());
	printf("with a final photon at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d\n",iTrack,mcp->Pt(), mcp->Eta(), mcp->Phi(),mcp->GetStatusCode());
      }
    }
    else{
      if(TMath::Abs(pdgMom)>100 && TMath::Abs(pdgMom)<1000)
	fDecayPhotonPtMC->Fill(mcp->Pt());
    }
  }
}
//________________________________________________________________________
Float_t AliAnalysisTaskEMCALIsoPhoton::GetClusSource(const AliVCluster *c)
{
  if(!c)
    return -1;
  if(!fStack)
    return -1;
  Int_t nmcp = fStack->GetNtrack();
  Int_t clabel = c->GetLabel();
  if(fDebug && fMcIdFamily.Contains(Form("%d",clabel)))
    printf("\n\t ==== Label %d is a descendent of the prompt photon ====\n\n",clabel);
  if(!fMcIdFamily.Contains(Form("%d",clabel)))
    return -1;
  fNClusForDirPho++;
  TString partonposstr = (TSubString)fMcIdFamily.operator()(0,1);
  Int_t partonpos = partonposstr.Atoi();
  if(fDebug)
    printf("\t^^^^ parton position = %d ^^^^\n",partonpos);
  if(clabel<0 || clabel>nmcp)
    return 0;
  Float_t clsPos[3] = {0,0,0};
  c->GetPosition(clsPos);
  TVector3 clsVec(clsPos);
  Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
  TParticle *mcp = static_cast<TParticle*>(fStack->Particle(partonpos));
  if(!mcp)
    return -1;
  if(fDebug){
    printf("\tclus mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcp->Eta(),mcp->Phi(),mcp->Energy(),mcp->GetPdgCode(),clabel);
  }
  Int_t lab1 =  mcp->GetFirstDaughter();
  if(lab1<0 || lab1>nmcp)
    return -1;
  TParticle *mcd = static_cast<TParticle*>(fStack->Particle(lab1));
  if(!mcd)
    return -1;
  if(fDebug)
    printf("\t\tmom mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcd->Eta(),mcd->Phi(),mcd->Energy(),mcd->GetPdgCode(),lab1);
  if(mcd->GetPdgCode()==22){
    fClusEtMcPt->Fill(mcd->Pt(), Et);
    fClusMcDetaDphi->Fill(clsVec.Eta() - mcd->Eta(), clsVec.Phi() - mcd->Phi());
  }
  else{
    printf("Warning: daughter of photon parton is not a photon\n");
    return -1;
  }
  return fDirPhoPt;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FollowGamma()
{
  if(!fStack)
    return;
  Int_t selfid = 6;
  TParticle *mcp = static_cast<TParticle*>(fStack->Particle(selfid));
  if(!mcp)
    return;
  if(mcp->GetPdgCode()!=22){
    mcp = static_cast<TParticle*>(fStack->Particle(++selfid));
    if(!mcp)
      return;
  }  
  Int_t daug0f =  mcp->GetFirstDaughter();
  Int_t daug0l =  mcp->GetLastDaughter();
  Int_t nd0 = daug0l - daug0f;
  if(fDebug)
    printf("\n\tGenerated gamma (%d) eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, n-daug=%d\n",selfid,mcp->Eta(),mcp->Phi(),mcp->Energy(),mcp->GetPdgCode(),nd0+1);
  fMcIdFamily = Form("%d,",selfid);
  GetDaughtersInfo(daug0f,daug0l, selfid,"");
  if(fDebug){
    printf("\t---- end of the prompt  gamma's family tree ----\n\n");
    printf("Family id string = %s,\n\n",fMcIdFamily.Data());
  }
  TParticle *md = static_cast<TParticle*>(fStack->Particle(daug0f));
  if(!md)
    return;
  fDirPhoPt = md->Pt();
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetDaughtersInfo(int firstd, int lastd, int selfid, const char *inputind)
{
  int nmcp = fStack->GetNtrack();
  if(firstd<0 || firstd>nmcp)
    return;
  if(lastd<0 || lastd>nmcp)
    return;
  if(firstd>lastd){
    printf("WARNING: First daughter > last (%d,%d)\n",firstd,lastd);
    return;
  }
  TString indenter = Form("\t%s",inputind);
  TParticle *dp = 0x0;
  if(fDebug)
    printf("\t%s--- Daughters of particle %d ---\n", indenter.Data(), selfid);
  for(int id=firstd; id<lastd+1; id++){
    dp =  static_cast<TParticle*>(fStack->Particle(id));
    if(!dp)
      continue;
    Int_t dfd = dp->GetFirstDaughter(); 
    Int_t dld = dp->GetLastDaughter();
    Int_t dnd =  dld - dfd + 1;
    Float_t vr = TMath::Sqrt(dp->Vx()*dp->Vx()+dp->Vy()*dp->Vy());
    if(fDebug)
      printf("\t%sParticle daughter(%d) eta=%1.1f,phi=%1.1f,E=%1.1f, vR=%1.1f, pdgcode=%d, n-daug=%d(%d,%d)\n", indenter.Data(),id, dp->Eta(), dp->Phi(), dp->Energy(), vr, dp->GetPdgCode(), dnd, dfd, dld);
    fMcIdFamily += Form("%d,",id);
    GetDaughtersInfo(dfd,dld,id,indenter.Data());
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskEMCALIsoPhoton::GetMcPtSumInCone(Float_t etaclus, Float_t phiclus, Float_t R)
{
  if(!fStack)
    return 0;
  if(fDebug)
    printf("Inside GetMcPtSumInCone!!\n");
  Int_t nTracks = fStack->GetNtrack();
  Float_t ptsum = 0;
  TString addpartlabels = "";
  for(Int_t iTrack=0;iTrack<nTracks;iTrack++){
    TParticle *mcp = static_cast<TParticle*>(fStack->Particle(iTrack));  
    if(!mcp)
      continue;  
    Int_t firstd = mcp->GetFirstDaughter();
    Int_t lastd = mcp->GetLastDaughter();
    if((iTrack==6 || iTrack==7) && mcp->GetPdgCode()==22){
     for(Int_t id=firstd;id<=lastd;id++)
      addpartlabels += Form("%d,",id);
     continue;
    }
   if(mcp->Rho()>2.5)
      continue;
    else {
      if(fDebug)
	printf("      >>>> mcp Rho, Vx, Vy, Vz = %1.1f,%1.1f,%1.1f,%1.1f.......\n",mcp->Rho(),mcp->Vx(), mcp->Vy(),mcp->Vz());
    }
    Float_t dphi = mcp->Phi() - phiclus;
    Float_t deta = mcp->Eta() - etaclus;
    if(fDebug)
      printf("      >>>> mcphi = %1.1f, mceta = %1.1f\n>>>> dphi = %1.1f, deta = %1.1f\n", mcp->Phi(), mcp->Eta(),dphi,deta);
    if(deta>10)
      continue;
    Float_t dR = TMath::Sqrt(dphi*dphi +  deta*deta);
    if(dR>R)
      continue;
    if(addpartlabels.Contains(Form("%d",iTrack))){
      if(fDebug)
	printf("      >>>> list of particles included (and their daughters) %s\n",addpartlabels.Data());
      continue;
    }
    addpartlabels += Form("%d,",iTrack);
    for(Int_t id=firstd;id<=lastd;id++)
      addpartlabels += Form("%d,",id);
    ptsum += mcp->Pt();
  }
  return ptsum;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
