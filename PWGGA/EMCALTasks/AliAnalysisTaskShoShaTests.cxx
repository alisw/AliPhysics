// $Id$
//
//task for EMCAL cluster shower shape investigations
//
// Authors: M.Cosentino

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDUtils.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "TLorentzVector.h"
#include "AliVCluster.h"
#include "AliAnalysisTaskShoShaTests.h"
#include "TFile.h"
#include "AliOADBContainer.h"


ClassImp(AliAnalysisTaskShoShaTests)

AliAnalysisTaskShoShaTests::AliAnalysisTaskShoShaTests() : 
  AliAnalysisTaskSE(), 
  fCaloClusters(0),
  fTracks(0),
  fEMCalCells(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fPeriod("LHC11c"),
  fIsTrain(0),
  fIsMC(0),
  fTrigThresh(4.8),
  fExoticCut(0.97),
  fEClusCut(0.5),
  fLowPi0MCut(0.12),
  fHighPi0MCut(0.16),
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fPIDResponse(0),
  fOutputList(0),
  fPvPos(0x0),
  fEvtSel(0),
  fPVZ(0),
  fClusEt(0),
  fClusEtTM(0),
  fClusEtLead(0),
  fClusEtSubLead(0),
  fClusEtLeadTM(0),
  fClusEtSubLeadTM(0),
  fClusEtExotic(0), 
  fClusEtExoticTM(0),
  fClusEtSingleExotic(0),
  fCellEnergy(0),
  fInvMassEMCNN(0),
  fM02Et(0),
  fM02EtPi0MassClCl(0),
  fM02EtPi0MassClClTruPi0(0),
  fM02EtPi0MassClClTruPiC(0),
  fM02EtPi0MassClClTruEta(0),
  fM02EtPi0MassClClTruK_0(0),
  fM02EtPi0MassClClTruK_C(0),
  fM02EtPi0MassClClTruPro(0),
  fM02EtPi0MassClClTruNeu(0),
  fM02EtTM(0),
  fM02EtExot(0),
  fM02EtExotTM(0),
  fElecNSigmaVsP(0),
  fPionNSigmaVsP(0),
  fKaonNSigmaVsP(0),
  fProtNSigmaVsP(0),
  fM02EtPion(0),
  fM02EtKaon(0),
  fM02EtProt(0),
  fM02EtElec(0)
{
  // Default constructor.
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
}

//________________________________________________________________________
AliAnalysisTaskShoShaTests::AliAnalysisTaskShoShaTests(const char *name) : 
  AliAnalysisTaskSE(name), 
  fCaloClusters(0),
  fTracks(0),
  fEMCalCells(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fPeriod("LHC11c"),
  fIsTrain(0),
  fIsMC(0),
  fTrigThresh(4.8),
  fExoticCut(0.97),
  fEClusCut(0.5),
  fLowPi0MCut(0.12),
  fHighPi0MCut(0.16),
  fESD(0),
  fMCEvent(0),
  fStack(0),
  fPIDResponse(0),
  fOutputList(0),
  fPvPos(0x0),
  fEvtSel(0),
  fPVZ(0),
  fClusEt(0),
  fClusEtTM(0),
  fClusEtLead(0),
  fClusEtSubLead(0),
  fClusEtLeadTM(0),
  fClusEtSubLeadTM(0),
  fClusEtExotic(0),
  fClusEtExoticTM(0), 
  fClusEtSingleExotic(0),
  fCellEnergy(0),
  fInvMassEMCNN(0),
  fM02Et(0),
  fM02EtPi0MassClCl(0),
  fM02EtPi0MassClClTruPi0(0),
  fM02EtPi0MassClClTruPiC(0),
  fM02EtPi0MassClClTruEta(0),
  fM02EtPi0MassClClTruK_0(0),
  fM02EtPi0MassClClTruK_C(0),
  fM02EtPi0MassClClTruPro(0),
  fM02EtPi0MassClClTruNeu(0),
  fM02EtTM(0),
  fM02EtExot(0),
  fM02EtExotTM(0),
  fElecNSigmaVsP(0),
  fPionNSigmaVsP(0),
  fKaonNSigmaVsP(0),
  fProtNSigmaVsP(0),
  fM02EtPion(0),
  fM02EtKaon(0),
  fM02EtProt(0),
  fM02EtElec(0)
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
void AliAnalysisTaskShoShaTests::UserCreateOutputObjects()
{
  // Create histograms, called once.
    
  fCaloClusters = new TRefArray();
  
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName.Data());
  fOADBContainer = new AliOADBContainer("AliEMCALgeo");
  fOADBContainer->InitFromFile(Form("$ALICE_PHYSICS/OADB/EMCAL/EMCALlocal2master.root"),"AliEMCALgeo");
  
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);

  fPVZ = new TH1F("hPVZ","distribution of pv-z;z (cm);counts",200,-20,20);
  fOutputList->Add(fPVZ);
  
  fClusEt = new TH1F("hClusEt","Clusters E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEt);
  
  fClusEtTM = new TH1F("hClusEtTM","Clusters (track-matched) E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtTM);
  
  fClusEtLead = new TH1F("hClusEtLead","Clusters (leading-trig) E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtLead);
  
  fClusEtSubLead = new TH1F("hClusEtSubLead","Clusters (subleading-trig) E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtSubLead);
  
  fClusEtLeadTM = new TH1F("hClusEtLeadTM","Clusters (leading-trig, TM) E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtLeadTM);
  
  fClusEtSubLeadTM = new TH1F("hClusEtSubLeadTM","Clusters (subleading-trig, TM) E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtSubLeadTM);
  
  fClusEtExotic = new TH1F("hClusEtExotic","Exotic trigger clusters  E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtExotic);
  
  fClusEtExoticTM = new TH1F("hClusEtExoticTM","Exotic trigger clusters (TM)  E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtExoticTM);
  
  fClusEtSingleExotic = new TH1F("hClusEtSingleExotic","Exotic trigger clusters (only this above thrs.)  E_{T} ;E_{T} ;dN/dE_{T}",400,0,200);
  fOutputList->Add(fClusEtSingleExotic);

  fCellEnergy = new TH1F("hCellE","cell energy spectrum;E_{cell} (GeV);entries",200,0,20);
  fOutputList->Add(fCellEnergy);

  fInvMassEMCNN = new TH2F("hInvMassEMCNN","inv mass of neutral EMC clusters pairs;m_{cc};n-entries",400,0,200,200,0,1);
  fOutputList->Add(fInvMassEMCNN);
  
  fM02Et = new TH2F("hM02Et","#lambda_{0}^{2} vs. E_{T} for trigger clusters ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02Et);
  
  fM02EtPi0MassClCl = new TH2F("fM02EtPi0MassClCl","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClCl);

  fM02EtPi0MassClClTruPi0  = new TH2F("fM02EtPi0MassClClTruPi0","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (pi0 mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruPi0);

  fM02EtPi0MassClClTruPiC  = new TH2F("fM02EtPi0MassClClTruPiC","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (pi+- mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruPiC);

  fM02EtPi0MassClClTruEta  = new TH2F("fM02EtPi0MassClClTruEta","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (eta mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruEta);

  fM02EtPi0MassClClTruK_0  = new TH2F("fM02EtPi0MassClClTruK_0","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (k0 mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruK_0);

  fM02EtPi0MassClClTruK_C  = new TH2F("fM02EtPi0MassClClTruK_C","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (k+- mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruK_C);

  fM02EtPi0MassClClTruPro  = new TH2F("fM02EtPi0MassClClTruPro","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (proton mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruPro);

  fM02EtPi0MassClClTruNeu  = new TH2F("fM02EtPi0MassClClTruNeu","#lambda_{0}^{2} vs. E_{T} for #pi^0 tagged clusters (neutron mc truth);E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPi0MassClClTruNeu);

  fM02EtTM = new TH2F("hM02EtTM","#lambda_{0}^{2} vs. E_{T} for trigger clusters(TM) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtTM);
  
  fM02EtExot = new TH2F("hM02EtExot","#lambda_{0}^{2} vs. E_{T} for trigger clusters(Exotic) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtExot);

  fM02EtExotTM = new TH2F("hM02EtExotTM","#lambda_{0}^{2} vs. E_{T} for trigger clusters(TM+Exotic) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtExotTM);

  fElecNSigmaVsP = new TH2F("hElecNSigmaVsP","n-#sigma electron vs momentum;track P (GeV/c);n-#sigma (electron)",200,0,20,200,-10,10);
  fOutputList->Add(fElecNSigmaVsP);

  fPionNSigmaVsP = new TH2F("hPionNSigmaVsP","n-#sigma pion vs momentum;track P (GeV/c);n-#sigma (electron)",200,0,20,200,-10,10);
  fOutputList->Add(fPionNSigmaVsP);

  fKaonNSigmaVsP = new TH2F("hKaonNSigmaVsP","n-#sigma kaon vs momentum;track P (GeV/c);n-#sigma (electron)",200,0,20,200,-10,10);
  fOutputList->Add(fKaonNSigmaVsP);

  fProtNSigmaVsP = new TH2F("hProtNSigmaVsP","n-#sigma proton vs momentum;track P (GeV/c);n-#sigma (electron)",200,0,20,200,-10,10);
  fOutputList->Add(fProtNSigmaVsP);
  
  fM02EtPion  = new TH2F("hM02EtPion","#lambda_{0}^{2} vs. E_{T} for #pi^#pm clusters (n-#sigma based) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtPion);

  fM02EtKaon  = new TH2F("hM02EtKaon","#lambda_{0}^{2} vs. E_{T} for kaon clusters(n-#sigma based) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtKaon);

  fM02EtProt  = new TH2F("hM02EtProt","#lambda_{0}^{2} vs. E_{T} for proton clusters(n-#sigma based) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtProt);

  fM02EtElec  = new TH2F("hM02EtElec","#lambda_{0}^{2} vs. E_{T} for electron clusters(n-#sigma based) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtElec);
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskShoShaTests::UserExec(Option_t *) 
{
  // User exec. Called once per event.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr == NULL) return;

  AliESDInputHandler *esdHandler = (AliESDInputHandler*)mgr->GetInputEventHandler();
  if (esdHandler == NULL) return;
  fPIDResponse = (AliPIDResponse*)esdHandler->GetPIDResponse();

  Bool_t isSelected = 0;
  /*if(fPeriod.Contains("11a"))
    isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
  else
    isSelected =  ((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral) ||
		   (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral));
  if(!isSelected )
    return; */

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  Int_t   runnumber = InputEvent()->GetRunNumber() ;
  fEvtSel->Fill(0);
  AliESDVertex *pv = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(!pv) 
    return;
  fPVZ->Fill(pv->GetZ());
  Double_t pvxyz[3] = {pv->GetX(), pv->GetY(), pv->GetZ()};
  fPvPos = (Double_t*)pvxyz;
  if(TMath::Abs(pv->GetZ())>15)
    return;
  fEvtSel->Fill(1);
  TObjArray *matEMCAL=(TObjArray*)fOADBContainer->GetObject(runnumber,"EmcalMatrices");
  for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
    if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
      break;
    fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod);
    fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
  }
  fESD->GetEMCALClusters(fCaloClusters);
  fEMCalCells = fESD->GetEMCALCells();
  for(int i=0;i<fEMCalCells->GetNumberOfCells();i++){
    Double_t e = fEMCalCells->GetCellAmplitude(TMath::Abs(fEMCalCells->GetAmplitude(i)));
    fCellEnergy->Fill(e);
  }

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
  fMCEvent = MCEvent(); 
  if(fMCEvent){
    fStack = (AliStack*)fMCEvent->Stack(); 
    fIsMC = kTRUE;
  }
  FillClusHists(); 
  
  fCaloClusters->Clear();
  PostData(1, fOutputList);
}
 
//________________________________________________________________________
void AliAnalysisTaskShoShaTests::FillClusHists()
{
  // Fill cluster histograms.

  if(!fCaloClusters)
    return;
  const Int_t nclus = fCaloClusters->GetEntries();
  if(nclus==0)
    return;
  Double_t EtArray;
  Bool_t isTM;
  Bool_t isEx;
  Int_t index;
  //printf("+++++++++\nstarting cluster loop!\n++++++++++++\n");
  for(int ic = 0; ic<nclus; ic++){
    AliESDCaloCluster *c1 = static_cast<AliESDCaloCluster*>(fCaloClusters->At(ic));
    if(!c1)
      continue;
    if(!c1->IsEMCAL())
      continue;
    if(c1->E()<fEClusCut)
      continue;
    Short_t id;
    Double_t Emax = GetMaxCellEnergy( c1, id);
    Double_t Ecross = GetCrossEnergy( c1, id);
    //printf("Ecross/Emax:%1.3f =================\n",Ecross/Emax);
    if((1.0-Ecross/Emax)>fExoticCut)
      isEx = 1;
    Float_t clsPos[3] = {0,0,0};
    c1->GetPosition(clsPos);
    TVector3 clsVecToZero(clsPos);
    TVector3 pvVec(fPvPos[0],fPvPos[1], fPvPos[2]);
    TVector3 clsVecToPv = clsVecToZero - pvVec;
    Double_t Et1 = c1->E()*TMath::Sin(clsVecToPv.Theta());
    fClusEt->Fill(Et1);
    fM02Et->Fill(Et1, c1->GetM02());
    if(isEx){
      fClusEtExotic->Fill(Et1);
      fM02EtExot->Fill(Et1,c1->GetM02()); 
    }
    Double_t dR = TMath::Sqrt(pow(c1->GetTrackDx(),2)+pow(c1->GetTrackDz(),2));
    if(dR<0.025){
      isTM=kTRUE;
      fClusEtTM->Fill(Et1);
      fM02EtTM->Fill(Et1, c1->GetM02());
      if(c1->GetNTracksMatched()<2){
	//TArrayI *fTracksMatched = (TArrayI*)c1->GetTracksMatched();
	FillChargedClustersShoSha(c1->GetTrackMatchedIndex(0), c1, Et1);
      }
      if(isEx){
	fClusEtExoticTM->Fill(Et1);
	fM02EtExotTM->Fill(Et1,c1->GetM02());
      }
    }
    Double_t clclmass=0;
    //cluster pairs both neutral and !exotic
    //if( (ic<(nclus-1)) && !isTM )// && !isEx ) 
      clclmass = NeutClusPairInvMass(c1, ic);
  }
}
 

//________________________________________________________________________
Double_t AliAnalysisTaskShoShaTests::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
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
Double_t AliAnalysisTaskShoShaTests ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
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
Double_t AliAnalysisTaskShoShaTests::NeutClusPairInvMass(const AliVCluster *cl1, Int_t ic)
{
  //printf("\nA neutral good cluster, looking for a pair.... \n\n");
  Double_t mass = -1;
  if(!cl1)
    return mass;
  if(!fCaloClusters)
    return mass;
  const Int_t nclus = fCaloClusters->GetEntries();
  if(ic==(nclus-1))
    return mass;
  if(nclus==0)
    return mass;
  for(int jc = ic+1; jc<nclus; jc++){
    AliESDCaloCluster *cl2 = static_cast<AliESDCaloCluster*>(fCaloClusters->At(jc));
    if(!cl2)
      continue;
    if(!cl2->IsEMCAL())
      continue;
    if(cl2->E()<fEClusCut)
      continue;
    Short_t id;
    Double_t Emax = GetMaxCellEnergy( cl2, id);
    Double_t Ecross = GetCrossEnergy( cl2, id);
    if((1-Ecross/Emax)>fExoticCut)
      continue;
    Double_t dR = TMath::Sqrt(pow(cl2->GetTrackDx(),2)+pow(cl2->GetTrackDz(),2));
    if(dR<0.03)
      continue;
    Float_t clsPos1[3] = {0,0,0};
    Float_t clsPos2[3] = {0,0,0};
    TVector3 pvVec(fPvPos[0],fPvPos[1], fPvPos[2]);
    cl1->GetPosition(clsPos1);
    TVector3 clsVec1(clsPos1);
    TVector3 clsVec1ToPv = clsVec1 - pvVec;
    Double_t Et1 = cl1->E()*TMath::Sin(clsVec1ToPv.Theta());
    cl2->GetPosition(clsPos2);
    TVector3 clsVec2(clsPos2);
    TVector3 clsVec2ToPv = clsVec2 - pvVec;
    Double_t Et2 = cl2->E()*TMath::Sin(clsVec2ToPv.Theta());
    TLorentzVector lv1, lv2, lvm;
    lv1.SetPtEtaPhiM(Et1,clsVec1ToPv.Eta(),clsVec1ToPv.Phi(),0.0);
    lv2.SetPtEtaPhiM(Et2,clsVec2ToPv.Eta(),clsVec2ToPv.Phi(),0.0);
    lvm = lv1 + lv2;
    mass = lvm.M();
    //printf("the pair mass is %1.2f ++++++++++++++++\n",mass);
    fInvMassEMCNN->Fill(lvm.Pt(),lvm.M());
    if(mass>fLowPi0MCut && mass<fHighPi0MCut){
      fM02EtPi0MassClCl->Fill(Et1,cl1->GetM02());
      fM02EtPi0MassClCl->Fill(Et2,cl2->GetM02());  
      if(fIsMC){
	Int_t ancesPdg = GetAncestorPdg(cl1->GetLabel());
	if(ancesPdg == 111 || ancesPdg == 11 || ancesPdg == 22)fM02EtPi0MassClClTruPi0->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 211)fM02EtPi0MassClClTruPiC->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 221)fM02EtPi0MassClClTruEta->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 311)fM02EtPi0MassClClTruK_0->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 321)fM02EtPi0MassClClTruK_C->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 2212)fM02EtPi0MassClClTruPro->Fill(Et1,cl1->GetM02());
	if(ancesPdg == 2112)fM02EtPi0MassClClTruNeu->Fill(Et1,cl1->GetM02());
	ancesPdg = GetAncestorPdg(cl2->GetLabel()); 
	if(ancesPdg == 111 || ancesPdg == 11 || ancesPdg == 2)fM02EtPi0MassClClTruPi0->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 211)fM02EtPi0MassClClTruPiC->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 221)fM02EtPi0MassClClTruEta->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 311)fM02EtPi0MassClClTruK_0->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 321)fM02EtPi0MassClClTruK_C->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 2212)fM02EtPi0MassClClTruPro->Fill(Et1,cl2->GetM02());
	if(ancesPdg == 2112)fM02EtPi0MassClClTruNeu->Fill(Et1,cl2->GetM02());
      }
    }
  }
  return mass;
}
//________________________________________________________________________
Int_t AliAnalysisTaskShoShaTests::GetAncestorPdg(const Int_t label)
{
  Int_t abspdg = 0, fpdg = 0, imom=label;
  //possible codes: pi0, pi+-, eta, k0, K+-, proton, neutron
  TString codes = "11 22 111 211 221 311 321 2212 2112";
  TString thecode ;
  if(!fStack)
    return fpdg;
  Int_t nTracks = fStack->GetNtrack();
  if(label<0 || label > nTracks)
    return fpdg;
  TParticle *mcp = 0x0;
  Int_t nbacksteps = 0;
  while(nbacksteps<3){
    mcp = static_cast<TParticle*>(fStack->Particle(imom));
    if(!mcp){
      nbacksteps++;
      continue;
    }
    abspdg = TMath::Abs(mcp->GetPdgCode());
    thecode = Form(" %d",abspdg);
    if(codes.Contains(thecode)){
      return abspdg;
    }
    imom = mcp->GetMother(0);
    nbacksteps++;
  }
  return fpdg;
}
//________________________________________________________________________
void AliAnalysisTaskShoShaTests::FillChargedClustersShoSha(Int_t trackLabel, AliESDCaloCluster *c, Double_t Et)
{
  if(!fTracks)
    return;
  AliESDtrack *track = static_cast<AliESDtrack*>(fTracks->At(trackLabel));
  if(!track)
    return;
  if(!fPIDResponse){
    return;
  }
  Float_t nselec=fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 
  Float_t nspion=fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion); 
  Float_t nskaon=fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon); 
  Float_t nsprot=fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton); 
  fElecNSigmaVsP->Fill(track->P(),nselec);
  fPionNSigmaVsP->Fill(track->P(),nspion);
  fKaonNSigmaVsP->Fill(track->P(),nskaon);
  fProtNSigmaVsP->Fill(track->P(),nsprot);
  if(TMath::Abs(nselec)>3 && TMath::Abs(nspion)<2 && TMath::Abs(nskaon)>3 && TMath::Abs(nsprot)>3)
    fM02EtPion->Fill(Et,c->GetM02());
  if(TMath::Abs(nselec)>3 && TMath::Abs(nspion)>3 && TMath::Abs(nskaon)<2 && TMath::Abs(nsprot)>3)
    fM02EtKaon->Fill(Et,c->GetM02());
  if(TMath::Abs(nselec)>3 && TMath::Abs(nspion)>3 && TMath::Abs(nskaon)>3 && TMath::Abs(nsprot)<2)
    fM02EtProt->Fill(Et,c->GetM02());
  if(TMath::Abs(nselec)<2 && TMath::Abs(nspion)>3 && TMath::Abs(nskaon)>3 && TMath::Abs(nsprot)>3)
    fM02EtElec->Fill(Et,c->GetM02());
  return;
}
//________________________________________________________________________
void AliAnalysisTaskShoShaTests::Terminate(Option_t *) 
{
  // Called once at the end of the query
}
