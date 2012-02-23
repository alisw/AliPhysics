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


#include "AliAnalysisTaskTrgContam.h"
#include "TFile.h"


ClassImp(AliAnalysisTaskTrgContam)

//________________________________________________________________________
AliAnalysisTaskTrgContam::AliAnalysisTaskTrgContam(const char *name) 
  : AliAnalysisTaskSE(name), 
  fCaloClusters(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fIsTrain(0),
  fTrigThresh(4.8),
  fExoticCut(0.97),

  
  fESD(0),
  
  fOutputList(0),
  
  fEvtSel(0),

  fClusEt(0),
  fClusEtTM(0),
  fClusEtLead(0),
  fClusEtSubLead(0),
  fClusEtLeadTM(0),
  fClusEtSubLeadTM(0),
  fClusEtExotic(0),
  fClusEtExoticTM(0), 
  fClusEtSingleExotic(0),
  fM02Et(0),
  fM02EtTM(0),
  fM02EtExot(0),
  fM02EtExotTM(0)


  
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
void AliAnalysisTaskTrgContam::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    
  fCaloClusters = new TRefArray();
  
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
    
  
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
  
  fM02Et = new TH2F("hM02Et","#lambda_{0}^{2} vs. E_{T} for trigger clusters ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02Et);
  
  fM02EtTM = new TH2F("hM02EtTM","#lambda_{0}^{2} vs. E_{T} for trigger clusters(TM) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtTM);
  
  fM02EtExot = new TH2F("hM02EtExot","#lambda_{0}^{2} vs. E_{T} for trigger clusters(Exotic) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtExot);

  fM02EtExotTM = new TH2F("hM02EtExotTM","#lambda_{0}^{2} vs. E_{T} for trigger clusters(TM+Exotic) ;E_{T} ;#lambda_{0}^{2}",400,0,200, 400,0,4);
  fOutputList->Add(fM02EtExotTM);
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskTrgContam::UserExec(Option_t *) 
{
  //event trigger selection
  Bool_t isSelected = 0;
  if(fPeriod.Contains("11a"))
    isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
  else
    isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7);
  if(!isSelected )
        return; 
  // Main loop
  // Called for each event

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  fEvtSel->Fill(0);
  
  AliESDVertex *pv = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(!pv) 
    return;
  if(TMath::Abs(pv->GetZ())>15)
    return;

  fEvtSel->Fill(1);


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
  PostData(1, fOutputList);
}      
//________________________________________________________________________
void AliAnalysisTaskTrgContam::FillClusHists()
{
  if(!fCaloClusters)
    return;
  const Int_t nclus = fCaloClusters->GetEntries();
  if(nclus==0)
    return;
  Double_t EtArray[nclus];
  Bool_t isTM[nclus];
  Bool_t isEx[nclus];
  Int_t index[nclus];
  Int_t nthresholds = 0;
  for(Int_t ic=0;ic<nclus;ic++){
    EtArray[ic]=0;
    isTM[ic] = 0;
    isEx[ic] = 0;
    index[ic]=0;
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(fCaloClusters->At(ic));
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    if(c->E()<fTrigThresh)
      continue;
    nthresholds++;
    Short_t id;
    Double_t Emax = GetMaxCellEnergy( c, id);
    Double_t Ecross = GetCrossEnergy( c, id);
    if((1-Ecross/Emax)>fExoticCut)
      isEx[ic] = 1;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
    EtArray[ic] = Et;
    fClusEt->Fill(Et);
    fM02Et->Fill(Et, c->GetM02());
    if(isEx[ic]){
      fClusEtExotic->Fill(Et);
      fM02EtExot->Fill(Et,c->GetM02()); 
     }
    Double_t dR = TMath::Sqrt(pow(c->GetTrackDx(),2)+pow(c->GetTrackDz(),2));
    if(dR<0.025){
      isTM[ic]=1;
      fClusEtTM->Fill(Et);
      fM02EtTM->Fill(Et, c->GetM02());
      if(isEx[ic]){
	fClusEtExoticTM->Fill(Et);
	fM02EtExotTM->Fill(Et,c->GetM02());
      }
    }
  }
  TMath::Sort(nclus,EtArray,index, kTRUE);
  if(EtArray[index[0]]>0){
    fClusEtLead->Fill(EtArray[index[0]]);
    if(nthresholds==1 && isEx[index[0]])
       fClusEtSingleExotic->Fill(EtArray[index[0]]);
  }
  if(nclus>1)if(EtArray[index[1]]>0)
    fClusEtSubLead->Fill(EtArray[index[1]]);
  if(isTM[index[0]] && EtArray[index[0]]>0)
    fClusEtLeadTM->Fill(EtArray[index[0]]);
  if(nclus>1)if(isTM[index[1]] && EtArray[index[1]]>0)
    fClusEtSubLeadTM->Fill(EtArray[index[1]]);
} 
//________________________________________________________________________
Double_t AliAnalysisTaskTrgContam::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
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
Double_t AliAnalysisTaskTrgContam ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
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
void AliAnalysisTaskTrgContam::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

}
