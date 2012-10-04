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
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fPVtxZ(0),  
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),
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
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fPVtxZ(0),            
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),   
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
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
  
  fNClusEt10 = new TH1F("hNClusEt10","# of cluster with E_{T}>10 per event;E;",101,-0.5,100.5);
  fOutputList->Add(fNClusEt10);
  
  fPVtxZ = new TH1F("hPVtxZ","primary vertex Z before cut;prim-vz(cm) ;",200,-100,100);
  fOutputList->Add(fPVtxZ);

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

  Int_t nEt=1000, nM02=400, nCeIso=1000, nTrIso=1000,  nAllIso=1000,  nCeIsoNoUE=1000,  nAllIsoNoUE=1000, nTrClDphi=200, nTrClDeta=100, nClEta=140, nClPhi=128, nTime=60, nMult=20;
  Int_t bins[] = {nEt, nM02, nCeIso, nTrIso, nAllIso, nCeIsoNoUE, nAllIsoNoUE, nTrClDphi, nTrClDeta,nClEta,nClPhi,nTime,nMult};
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  Double_t xmin[] = { 0.,   0.,  -10.,   -10., -10., -10., -10., -0.1,-0.05, -0.7, 1.4,-0.15e-06,-0.5};
  Double_t xmax[] = { 100., 4., 190., 190., 190.,  190., 190., 0.1, 0.05, 0.7, 3.192, 0.15e-06,99.5};
  fHnOutput =  new THnSparseF("fHnOutput","Output matrix: E_{T},M02,CeIso,TrIso,AllIso, CeIsoNoUESub, AllIsoNoUESub, d#phi_{trk},d#eta_{trk},#eta_{clus},#phi_{clus},T_{max},mult", ndims, bins, xmin, xmax);
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
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  }
  if(fIsMc)
    isSelected = kTRUE;
  if(fDebug)
    printf("isSelected = %d, fIsMC=%d\n", isSelected, fIsMc);
  if(!isSelected )
        return; 

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
  fPVtxZ->Fill(pv->GetZ());
  if(TMath::Abs(pv->GetZ())>15)
    return;
  if(fDebug)
    printf("passed vertex cut\n");

  fEvtSel->Fill(1);

  if (!fTracks)  
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
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

  fESD->GetEMCALClusters(fCaloClusters);
  fEMCalCells = fESD->GetEMCALCells();
  
    
  FillClusHists(); 
  if(fDebug)
    printf("passed calling of FillClusHists\n");

  fCaloClusters->Clear();
  fSelPrimTracks->Clear();

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
  FillMcHists();
  if(fDebug)
    printf("passed calling of FillMcHists\n");


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
    if(Et>10)
      nclus10++;
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
    const Int_t ndims =   fNDimensions;
    Double_t outputValues[ndims];
    outputValues[0] = Et;
    outputValues[1] = c->GetM02();
    outputValues[2] = ceiso-cecore-ceisoue;
    outputValues[3] = triso-trisoue;
    outputValues[4] = alliso-cecore-allisoue;
    outputValues[5] = ceiso-Et;
    outputValues[6] = alliso-Et;
    outputValues[7] = c->GetTrackDx();
    outputValues[8] = c->GetTrackDz();
    outputValues[9] = clsVec.Eta();
    outputValues[10] = clsVec.Phi();
    outputValues[11] = fEMCalCells->GetCellTime(id);
    outputValues[12] = fTrackMult;
    fHnOutput->Fill(outputValues);
    if(c->E()>maxE)
      maxE = c->E();
  }
  fNClusHighClusE->Fill(maxE,nclus);
  fNClusEt10->Fill(nclus10);
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
    AliVTrack *track = static_cast<AliVTrack*> (fSelPrimTracks->At(itrack));
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
void AliAnalysisTaskEMCALIsoPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
