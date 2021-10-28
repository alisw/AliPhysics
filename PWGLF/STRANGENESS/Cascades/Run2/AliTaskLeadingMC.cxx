#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliTrackReference.h"
#include "AliTaskLeadingMC.h"
#include "AliGenPythiaEventHeader.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliPWG0Helper.h"
#include "AliAnalysisFilter.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliVEvent.h"


#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

ClassImp(AliTaskLeadingMC)

//________________________________________________________________________
AliTaskLeadingMC::AliTaskLeadingMC(const char *name) :
AliAnalysisTask(name, ""),
fTrackFilter(0x0),
fSpherocity(-10),
fNTracksSpherocity(0)
{
  // Constructor

  SetZDCPGeo();
  SetZDCNGeo();
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  DefineOutput(0, TTree::Class());// crea una slot di output "0" per un oggetto della classe TH1D
  DefineOutput(1, TH3F::Class());// crea una slot di output "0" per un oggetto della classe TH1D
  DefineOutput(2, TH3F::Class());// crea una slot di output "0" per un oggetto della classe TH1D
  DefineOutput(3, TH2F::Class());// crea una slot di output "0" per un oggetto della classe TH1D
  DefineOutput(4, TH2F::Class());// crea una slot di output "0" per un oggetto della classe TH1D
}
//________________________________________________________________________
void AliTaskLeadingMC::SetZDCPGeo(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax){
  fZDCP_Xmin = xmin;
  fZDCP_Xmax = xmax;
  fZDCP_Ymin = ymin;
  fZDCP_Ymax = ymax;
  fZDCP_Zmin = zmin;
  fZDCP_Zmax = zmax;
}
//________________________________________________________________________
void AliTaskLeadingMC::SetZDCNGeo(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax){
  fZDCN_Xmin = xmin;
  fZDCN_Xmax = xmax;
  fZDCN_Ymin = ymin;
  fZDCN_Ymax = ymax;
  fZDCN_Zmin = zmin;
  fZDCN_Zmax = zmax;
}
//________________________________________________________________________
void AliTaskLeadingMC::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
//     tree->SetBranchStatus("*", kTRUE);
//     tree->SetBranchStatus("Tracks.*", kTRUE);
//     tree->SetBranchStatus("*SPDV*", kTRUE);
//     tree->SetBranchStatus("AliESDTZERO*", kTRUE);
//     tree->SetBranchStatus("AliESDRun*", kTRUE);
//     tree->SetBranchStatus("AliESDHeader*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = (AliESDEvent *) esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliTaskLeadingMC::CreateOutputObjects()
{
  // Called once
  fDB = TDatabasePDG::Instance();

  //initialize quality trackcuts for spherocity
  if(!fTrackFilter){	
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    SetTrackCuts(fTrackFilter);
  }
  
  //Create output tree
  fTree=new TTree("fTree", "ZDC");

  // check for availability of track refs
  fTree->Branch("isTR", &fIsTrackRef, "isTR/I");

  // global observables MC generator
  fTree->Branch("v0mPerc", &fV0Perc, "v0mPerc/F");
  fTree->Branch("zdcPerc", &fZdcPerc, "zdcPerc/F");
  fTree->Branch("zdcPercFired", &fZdcPercFired, "zdcPercFired/F");
  fTree->Branch("multRef5", &fMultRef5, "multRef5/F");
  fTree->Branch("multRef8", &fMultRef8, "multRef8/F");
  fTree->Branch("multSPDcl", &fMultSPDcl, "multSPDcl/F");
  fTree->Branch("multSPDtr", &fMultSPDtr, "multSPDtr/F");
  fTree->Branch("inelGT0",&fInelGT0,"inelGT0/I");
  fTree->Branch("SPDtracklets", &fSPDtracklets, "SPDtracklets/I");
  fTree->Branch("SPDtrackletsA", &fSPDtrackletsA, "SPDtrackletsA/I");
  fTree->Branch("SPDtrackletsC", &fSPDtrackletsC, "SPDtrackletsC/I");
  fTree->Branch("TOFclusters", &fTOFclusters, "TOFclusters/I");
  fTree->Branch("TOFclustersTrg", &fTOFclustersTrg, "TOFclustersTrg/I");
  fTree->Branch("nch", &fNch, "nch/I");
  fTree->Branch("nchEta", &fNchEta, "nchEta/I");
  fTree->Branch("nchEtaA", &fNchEtaA, "nchEtaA/I");
  fTree->Branch("nchEtaC", &fNchEtaC, "nchEtaC/I");
  fTree->Branch("energyhEta", &fEnergyEta, "energyEta/F");  
  fTree->Branch("nmpi", &fNMPI, "nmpi/I");
  fTree->Branch("nLambdaEta", &fNLambdaEta, "nLambdaEta/I");
  fTree->Branch("nXiEta", &fNXiEta, "nXiEta/I");
  fTree->Branch("ptXiEta", fPtXiEta, "ptXiEta[nXiEta]/F");
  fTree->Branch("nAntiXiEta", &fNAntiXiEta, "nAntiXiEta/I");
  fTree->Branch("ptAntiXiEta", fPtAntiXiEta, "ptAntiXiEta[nAntiXiEta]/F");
  fTree->Branch("nXiEtaFrag", &fNXiEtaFrag, "nXiEtaFrag/I");
  fTree->Branch("ptXiEtaFrag", fPtXiEtaFrag, "ptEtaFrag[nXiEtaFrag]/F");
  fTree->Branch("nXiEtaUp", &fNXiEtaUp, "nXiEtaUp/I");
  fTree->Branch("ptXiEtaUp", fPtXiEtaUp, "ptEtaUp[nXiEtaUp]/F");
  fTree->Branch("nXiEtaDown", &fNXiEtaDown, "nXiEtaDown/I");
  fTree->Branch("ptXiEtaDown", fPtXiEtaDown, "ptEtaDown[nXiEtaDown]/F");
  fTree->Branch("nOmegaEta", &fNOmegaEta, "nOmegaEta/I");
  fTree->Branch("nPiEta", &fNPiEta, "nPiEta/I");
  fTree->Branch("nPi0Eta", &fNPi0Eta, "nPi0Eta/I");
  fTree->Branch("nKchEta", &fNKchEta, "nKchEta/I");
  fTree->Branch("nK0Eta", &fNK0Eta, "nK0Eta/I");
  fTree->Branch("sumLambdaXi", &fSumPtLambdaEta, "sumPtLambda/F");
  fTree->Branch("sumPtXi", &fSumPtXiEta, "sumPtXi/F");
  fTree->Branch("sumPtOmega", &fSumPtOmegaEta, "sumPtOmega/F");
  fTree->Branch("sumPtPi", &fSumPtPiEta, "sumPtPi/F");
  fTree->Branch("maxChargePt", &fMaxChargePt, "maxChargePt/F");
  fTree->Branch("effEnergy",&fEffEnergy, "effEnergy/F");
  fTree->Branch("spherocity",&fSpherocity, "spherocity/F");
  fTree->Branch("ntracksspherocity",&fNTracksSpherocity, "ntracksspherocity/I");


  // ZDC recon infos
  fTree->Branch("adcZDCN1",fAdcZDCN1,"adcZDCN1[5]/F");
  fTree->Branch("adcZDCP1",fAdcZDCP1,"adcZDCP1[5]/F");
  fTree->Branch("adcZDCN2",fAdcZDCN2,"adcZDCN2[5]/F");
  fTree->Branch("adcZDCP2",fAdcZDCP2,"adcZDCP2[5]/F");

  // particles with hits
  fTree->Branch("fP1hits", &fP1hits, "fP1hits/I");
  fTree->Branch("fP2hits", &fP2hits, "fP2hits/I");
  fTree->Branch("fN1hits", &fN1hits, "fN1hits/I");
  fTree->Branch("fN2hits", &fN2hits, "fN2hits/I");
  fTree->Branch("E_p_leadC", fE_p_leadC, "E_p_leadC[fP1hits]/F");
  fTree->Branch("E_p_leadA", fE_p_leadA, "E_p_leadA[fP2hits]/F");
  fTree->Branch("E_n_leadC", fE_n_leadC, "E_n_leadC[fN1hits]/F");
  fTree->Branch("E_n_leadA", fE_n_leadA, "E_n_leadA[fN2hits]/F");
  fTree->Branch("pdg_leadingP1", fPdg_leadingP1, "pdg_leadingP1[fP1hits]/I");
  fTree->Branch("pdg_leadingP2", fPdg_leadingP2, "pdg_leadingP2[fP2hits]/I");
  fTree->Branch("pdg_leadingN1", fPdg_leadingN1, "pdg_leadingN1[fN1hits]/I");
  fTree->Branch("pdg_leadingN2", fPdg_leadingN2, "pdg_leadingN2[fN2hits]/I");
  fTree->Branch("labelHitP1", fLabel_leadP1, "labelHitP1[fP1hits]/I");
  fTree->Branch("labelHitP2", fLabel_leadP2, "labelHitP2[fP2hits]/I");
  fTree->Branch("labelHitN1", fLabel_leadN1, "labelHitN1[fN1hits]/I");
  fTree->Branch("labelHitN2", fLabel_leadN2, "labelHitN2[fN2hits]/I");	
  fTree->Branch("labelMP1", fLabelMP1, "labelMP1[fP1hits]/I");
  fTree->Branch("labelMP2", fLabelMP2, "labelMP2[fP2hits]/I");
  fTree->Branch("labelMN1", fLabelMN1, "labelMN1[fN1hits]/I");
  fTree->Branch("labelMN2", fLabelMN2, "labelMN2[fN2hits]/I");
  fTree->Branch("ZP1", fZP1, "ZP1[fP1hits]/F");
  fTree->Branch("ZP2", fZP2, "ZP2[fP2hits]/F");
  fTree->Branch("ZN1", fZN1, "ZN1[fN1hits]/F");
  fTree->Branch("ZN2", fZN2, "ZN2[fN2hits]/F");
  fTree->Branch("XP1", fXP1, "XP1[fP1hits]/F");
  fTree->Branch("XP2", fXP2, "XP2[fP2hits]/F");
  fTree->Branch("XN1", fXN1, "XN1[fN1hits]/F");
  fTree->Branch("XN2", fXN2, "XN2[fN2hits]/F");	
  fTree->Branch("YP1", fYP1, "YP1[fP1hits]/F");
  fTree->Branch("YP2", fYP2, "YP2[fP2hits]/F");
  fTree->Branch("YN1", fYN1, "YN1[fN1hits]/F");
  fTree->Branch("YN2", fYN2, "YN2[fN2hits]/F");
  fTree->Branch("momHitN1x", fMomN1x, "momHitN1x[fN1hits]/F");
  fTree->Branch("momHitN2x", fMomN2x, "momHitN2x[fN2hits]/F");
  fTree->Branch("momHitP1x", fMomP1x, "momHitP1x[fP1hits]/F");
  fTree->Branch("momHitP2x", fMomP2x, "momHitP2x[fP2hits]/F");
  fTree->Branch("momHitN1y", fMomN1y, "momHitN1y[fN1hits]/F");
  fTree->Branch("momHitN2y", fMomN2y, "momHitN2y[fN2hits]/F");
  fTree->Branch("momHitP1y", fMomP1y, "momHitP1y[fP1hits]/F");
  fTree->Branch("momHitP2y", fMomP2y, "momHitP2y[fP2hits]/F");
  fTree->Branch("momHitN1z", fMomN1z, "momHitN1z[fN1hits]/F");
  fTree->Branch("momHitN2z", fMomN2z, "momHitN2z[fN2hits]/F");
  fTree->Branch("momHitP1z", fMomP1z, "momHitP1z[fP1hits]/F");
  fTree->Branch("momHitP2z", fMomP2z, "momHitP2z[fP2hits]/F");
  
  // all candidate with high pseudorapidities
  fTree->Branch("p_cand_leadC", &fP_cand_leadC, "p_cand_leadC/I");
  fTree->Branch("p_cand_leadA", &fP_cand_leadA, "p_cand_leadA/I");
  fTree->Branch("n_cand_leadC", &fN_cand_leadC, "n_cand_leadC/I");
  fTree->Branch("n_cand_leadA", &fN_cand_leadA, "n_cand_leadA/I");
  fTree->Branch("E_p_cand_leadC", fE_p_cand_leadC, "E_p_cand_leadC[p_cand_leadC]/F");
  fTree->Branch("E_p_cand_leadA", fE_p_cand_leadA, "E_p_cand_leadA[p_cand_leadA]/F"); //energia dei protoni candidati leading nella regione anticlockwise
  fTree->Branch("E_n_cand_leadC", fE_n_cand_leadC, "E_n_cand_leadC[n_cand_leadC]/F");
  fTree->Branch("E_n_cand_leadA", fE_n_cand_leadA, "E_n_cand_leadA[n_cand_leadA]/F");
  fTree->Branch("pdg_cand_leadP1", fPdg_cand_leadP1, "pdg_cand_leadP1[p_cand_leadC]/I");
  fTree->Branch("pdg_cand_leadP2", fPdg_cand_leadP2, "pdg_cand_leadP2[p_cand_leadA]/I");
  fTree->Branch("pdg_cand_leadN1", fPdg_cand_leadN1, "pdg_cand_leadN1[n_cand_leadC]/I");
  fTree->Branch("pdg_cand_leadN2", fPdg_cand_leadN2, "pdg_cand_leadN2[n_cand_leadA]/I");
  fTree->Branch("iLabel_cand_P1", fLabel_cand_P1, "iLabel_cand_P1[p_cand_leadC]/I");
  fTree->Branch("iLabel_cand_P2", fLabel_cand_P2, "iLabel_cand_P2[p_cand_leadA]/I");
  fTree->Branch("iLabel_cand_N1", fLabel_cand_N1, "iLabel_cand_N1[n_cand_leadC]/I");
  fTree->Branch("iLabel_cand_N2", fLabel_cand_N2, "iLabel_cand_N2[n_cand_leadA]/I");
  fTree->Branch("mom_cand_N1x", fMomCN1x, "mom_cand_N1x[n_cand_leadC]/F");
  fTree->Branch("mom_cand_N1y", fMomCN1y, "mom_cand_N1y[n_cand_leadC]/F");
  fTree->Branch("mom_cand_N1z", fMomCN1z, "mom_cand_N1z[n_cand_leadC]/F");
  fTree->Branch("mom_cand_N2x", fMomCN2x, "mom_cand_N2x[n_cand_leadA]/F");
  fTree->Branch("mom_cand_N2y", fMomCN2y, "mom_cand_N2y[n_cand_leadA]/F");
  fTree->Branch("mom_cand_N2z", fMomCN2z, "mom_cand_N2z[n_cand_leadA]/F");
  fTree->Branch("mom_cand_P1x", fMomCP1x, "mom_cand_P1x[p_cand_leadC]/F");
  fTree->Branch("mom_cand_P1y", fMomCP1y, "mom_cand_P1y[p_cand_leadC]/F");
  fTree->Branch("mom_cand_P1z", fMomCP1z, "mom_cand_P1z[p_cand_leadC]/F");
  fTree->Branch("mom_cand_P2x", fMomCP2x, "mom_cand_P2x[p_cand_leadA]/F");
  fTree->Branch("mom_cand_P2y", fMomCP2y, "mom_cand_P2y[p_cand_leadA]/F");
  fTree->Branch("mom_cand_P2z", fMomCP2z, "mom_cand_P2z[p_cand_leadA]/F");
  fTree->Branch("iLabelM_cand_P1", fLabelM_cand_P1, "iLabelM_cand_P1[p_cand_leadC]/I");
  fTree->Branch("iLabelM_cand_P2", fLabelM_cand_P2, "iLabelM_cand_P2[p_cand_leadA]/I");
  fTree->Branch("iLabelM_cand_N1", fLabelM_cand_N1, "iLabelM_cand_N1[n_cand_leadC]/I");
  fTree->Branch("iLabelM_cand_N2", fLabelM_cand_N2, "iLabelM_cand_N2[n_cand_leadA]/I");
  fTree->Branch("pdgM_cand_leadP1", fPdgM_cand_leadP1, "pdgM_cand_leadP1[p_cand_leadC]/I");
  fTree->Branch("pdgM_cand_leadP2", fPdgM_cand_leadP2, "pdgM_cand_leadP2[p_cand_leadA]/I");
  fTree->Branch("pdgM_cand_leadN1", fPdgM_cand_leadN1, "pdgM_cand_leadN1[n_cand_leadC]/I");
  fTree->Branch("pdgM_cand_leadN2", fPdgM_cand_leadN2, "pdgM_cand_leadN2[n_cand_leadA]/I");
  
  PostData(0, fTree);	

  fH_SPD_VZERO = new TH3F("hSPD_VZERO",";#eta;SPD perc (%); VOM perc (%)",200,-10,10,100,0,100,100,0,100);
  PostData(1, fH_SPD_VZERO);
  fH_SPD_ZDC = new TH3F("hSPD_ZDC",";#eta;SPD perc (%); ZDC perc (%)",200,-10,10,100,0,100,100,0,100);
  PostData(2, fH_SPD_ZDC);
  fH_SPD_VZERO_ev = new TH2F("hSPD_VZERO_ev",";SPD perc (%); VOM perc (%)",100,0,100,100,0,100);
  PostData(3, fH_SPD_VZERO_ev);
  fH_SPD_ZDC_ev = new TH2F("hSPD_ZDC_ev",";SPD perc (%); ZDC perc (%)",100,0,100,100,0,100);
  PostData(4, fH_SPD_ZDC_ev);
}
//________________________________________________________________________
void AliTaskLeadingMC::Exec(Option_t *) 
{	
  AliMCEvent* mcEvent=NULL;
  AliMCEventHandler* eventHandler=NULL;
  TTree* treeTR=NULL;
  fIsTrackRef = fAskTrackRef;
  
  // recupera le info MC
  eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    Printf("ERROR: Could not retrieve MC event handler");
    return;
  }
  
  mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  // Get TrackRefs
  if(fIsTrackRef){
    treeTR = eventHandler->TreeTR();
    if(!treeTR){
      fIsTrackRef = 0;
    }
    else{
      Printf("TrackReference found!");
    }
  }

  // Check ESD event
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }  

  fTOFclusters = fESD->GetTOFHeader()->GetNumberOfTOFclusters();
  fTOFclustersTrg = fESD->GetTOFHeader()->GetNumberOfTOFtrgPads();

  fInelGT0 = 0;

  // percentile selections ($ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C required)
  fV0Perc = 500;
  fZdcPerc = 500;
  AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
  if( !MultSelection) {
     //If you get this warning please check that the AliMultSelectionTask actually ran (before your task)
     AliWarning("AliMultSelection object not found!");
  } else {
     //V0M Multiplicity 
     fV0Perc = MultSelection->GetMultiplicityPercentile("V0M");
     fZdcPercFired = MultSelection->GetMultiplicityPercentile("ZPNACpp");
     fZdcPerc = MultSelection->GetMultiplicityPercentile("ZPNACTowerpp");
     fMultRef5 = MultSelection->GetMultiplicityPercentile("RefMult05");
     fMultRef8 = MultSelection->GetMultiplicityPercentile("RefMult08");
     fMultSPDcl = MultSelection->GetMultiplicityPercentile("SPDClusters");
     fMultSPDtr = MultSelection->GetMultiplicityPercentile("SPDTracklets");
  }

  //SPDTracklets to get a midrapidity Nch estimator
  AliMultiplicity* multiplicity =  fESD->GetMultiplicity();
  fSPDtracklets = 0;
  fSPDtrackletsA = 0;
  fSPDtrackletsC = 0;
  for (auto it = 0; it<multiplicity->GetNumberOfTracklets(); it++) {
     Double_t eta = multiplicity->GetEta(it);
     if ( abs(eta) < 0.5 ){
       fSPDtracklets++;
       if(eta < 0) fSPDtrackletsA++;
       else fSPDtrackletsC++;
     }
  }

  fSpherocity = ComputeSpherocity(); //Calculates event spherocity

  // fill signal reconstructed in ZDCs
  fillZDCreco();

  // MC generator general infos
  // to be check if working also for other generators
  AliGenPythiaEventHeader *mcheader = (AliGenPythiaEventHeader *) mcEvent->GenEventHeader();
  fNMPI = mcheader->GetNMPI(); // N. multi-parton interactions

  loopMC(mcEvent);

  //  loopTrack(mcEvent);
    
  // look into track refs for ZDC hits
  if(fIsTrackRef) loopTrackRef(treeTR, mcEvent);
  else{
    // reset variables
    fP2hits=0;
    fP1hits=0;
    fN2hits=0;
    fN1hits=0;
  }
  
  fTree->Fill();
  PostData(0, fTree);
  PostData(1, fH_SPD_VZERO);
  PostData(2, fH_SPD_ZDC);
  PostData(3, fH_SPD_VZERO_ev);
  PostData(4, fH_SPD_ZDC_ev);
}  
//________________________________________________________________________
void AliTaskLeadingMC::loopMC(AliMCEvent *mcEvent){
  // # of reco tracks and primaries from MC
  Int_t ntra = mcEvent->GetNumberOfTracks();

  fMaxChargePt = 0;
  fEffEnergy = 0;
  
  // reset variables
  fP_cand_leadA=0;
  fN_cand_leadA=0;
  fP_cand_leadC=0;
  fN_cand_leadC=0;

  fNch=0,fNchEta=0,fNchEtaA=0,fNchEtaC=0,fNLambdaEta=0,fNXiEta=0,fNAntiXiEta=0,fNXiEtaFrag=0,fNXiEtaUp=0,fNXiEtaDown=0,fNOmegaEta=0,fNPiEta=0,fNPi0Eta=0,fNKchEta=0,fNK0Eta=0;
  fSumPtLambdaEta=fSumPtXiEta=fSumPtOmegaEta=fSumPtPiEta=0;
  fEnergyEta=0;

  if(fSPDtracklets){
    fH_SPD_VZERO_ev->Fill(fMultSPDcl,fV0Perc);
    fH_SPD_ZDC_ev->Fill(fMultSPDcl,fZdcPerc);
  }

  // loop MC particles
  Int_t nPrim = mcEvent->Stack()->GetNprimary();
  for (Int_t i = 0; i < nPrim; i++){
    TParticle *part = mcEvent->Stack()->Particle(i);
    if (!AliPWG0Helper::IsPrimaryCharged(part, nPrim)) continue;
    Double_t eta = part->Eta();
    if (fabs(eta) < 1.0) fInelGT0 = true;
  }

  Float_t px,py,pz,pt;
  for(Int_t i = 0; i < ntra ;i++){
    // get particle from stack
    AliMCParticle *MCpart = (AliMCParticle *) mcEvent->GetTrack(i);
    TParticle *part = MCpart->Particle(); 
    
    Int_t status = (part->GetStatusCode() == 1);
    
    // get momentum
    px = part->Px();
    py = part->Py();
    pz = part->Pz();
    pt = sqrt(px*px + py*py);
    
    Int_t charge = 0;

    if(fDB->GetParticle(part->GetPdgCode())) charge = Int_t(fDB->GetParticle(part->GetPdgCode())->Charge());
    else continue; // skip particles with undefined pdg

    if(charge) {
      status = AliPWG0Helper::IsPrimaryCharged(part, nPrim); // official definition of charged primary
    }
    else{
      if (part->GetFirstDaughter() != -1 && part->GetFirstDaughter() < nPrim) status = 0;
    }

    if(status && charge){ // charged particles
      fNch++;

      if(fSPDtracklets){
        if(part->Eta() < -9.999){
          fH_SPD_VZERO->Fill(-9.999,fMultSPDcl,fV0Perc);
          fH_SPD_ZDC->Fill(-9.999,fMultSPDcl,fZdcPerc);
        }
        else if(part->Eta() > 9.999){
          fH_SPD_VZERO->Fill(9.999,fMultSPDcl,fV0Perc);
          fH_SPD_ZDC->Fill(9.999,fMultSPDcl,fZdcPerc);
        }
        else{
          fH_SPD_VZERO->Fill(part->Eta(),fMultSPDcl,fV0Perc);
          fH_SPD_ZDC->Fill(part->Eta(),fMultSPDcl,fZdcPerc);
        }
      }

      if(TMath::Abs(part->Eta())<fEtaBarrel){
	fEnergyEta+=part->Energy();
	fNchEta++;
        if(part->Eta() < 0) fNchEtaA++;
        else fNchEtaC++;

	if(pt > fMaxChargePt) fMaxChargePt = pt;
	
	if(TMath::Abs(part->GetPdgCode()) == 211){ // charged pions
	  fNPiEta++;
	  fSumPtPiEta += pt;
	}
        else if(TMath::Abs(part->GetPdgCode()) == 321){ // neutral kaons
          fNKchEta++;
        }
      }
    }
    else if(! charge){
      if(TMath::Abs(part->Eta())<fEtaBarrel){
        if(TMath::Abs(part->GetPdgCode()) == 111){ // neutral pions
          fNPi0Eta++;
        }
        else if(TMath::Abs(part->GetPdgCode()) == 311){ // K0s neutral kaons
          fNK0Eta++;
        }
      }
    }

    if(TMath::Abs(part->Eta())<fEtaBarrel){
      if(TMath::Abs(part->GetPdgCode()) == 3122){ // lambda
	fNLambdaEta++;
	fSumPtLambdaEta += pt;
      }
      if(TMath::Abs(part->GetPdgCode()) == 3312){ // Xi
        if(part->GetPdgCode() > 0){
          fPtXiEta[fNXiEta] = pt;
          fNXiEta++;
        }
        else{
          fPtXiEta[fNAntiXiEta] = pt;
          fNAntiXiEta++;
        }
	fSumPtXiEta += pt;
        // look if it comes from fragmenetation
        Int_t imoth = part->GetFirstMother();
        AliMCParticle *partMCM;
        TParticle *partM;
        while(imoth >= 0 && fNXiEtaFrag < 100){
          partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);  
          partM = partMCM->Particle();
          if(TMath::Abs(partM->GetPdgCode()) == 3){
            fPtXiEtaFrag[fNXiEtaFrag] = pt;
            fNXiEtaFrag++;
            imoth = -1;
          }
          else if(TMath::Abs(partM->GetPdgCode()) == 1){
            fPtXiEtaUp[fNXiEtaUp] = pt;
            fNXiEtaUp++;
            imoth = -1;
          }
          else if(TMath::Abs(partM->GetPdgCode()) == 2){
            fPtXiEtaDown[fNXiEtaDown] = pt;
            fNXiEtaDown++;
            imoth = -1;
          }
          else{
            imoth = partM->GetFirstMother();
          }
       }
      }
      if(TMath::Abs(part->GetPdgCode()) == 3334){ // Omega
	fNOmegaEta++;
	fSumPtOmegaEta += pt;
     }
    }
    
    // keep only stable particles at forward pseudorapidity
    if(status!=1 || TMath::Abs(part->Eta())<fEtaThreshold || part->Energy() < fEnergyThreshold) continue;
    
    Int_t imoth = part->GetFirstMother();
    
    Int_t statusM = 0;
    Int_t pdgM = 0;
    
    AliMCParticle *partMCM;
    if(imoth>7){ // first 8 usually refered to initial incident protons
      partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);  
      TParticle *partM = partMCM->Particle();
      
      statusM = partM->GetStatusCode();
      
      if(statusM==11 &&  TMath::Abs(partM->GetPdgCode())>99) pdgM = partM->GetPdgCode();
      else statusM = 0;
    }
    
    if(statusM != 11) imoth = -1; // only decay accepted to define mother

    fEffEnergy -= part->Energy();
    
    if((part->Eta())>0.0) { // positive Z (A side)
      if(charge && fP_cand_leadA < fgkDim){ // charged candidate
	fE_p_cand_leadA[fP_cand_leadA]=part->Energy();
	fPdg_cand_leadP2[fP_cand_leadA]=part->GetPdgCode();
	fPdgM_cand_leadP2[fP_cand_leadA]=pdgM;
	fLabel_cand_P2[fP_cand_leadA]=i;
	fLabelM_cand_P2[fP_cand_leadA]=imoth;
	fMomCP2x[fP_cand_leadA] = px;
	fMomCP2y[fP_cand_leadA] = py;
	fMomCP2z[fP_cand_leadA] = pz;
	fP_cand_leadA++;
      } 
      else if(!charge && fN_cand_leadA < fgkDim){ // neutral candidate
	fE_n_cand_leadA[fN_cand_leadA]=part->Energy();
	fPdg_cand_leadN2[fN_cand_leadA]=part->GetPdgCode(); 
	fPdgM_cand_leadN2[fN_cand_leadA]=pdgM;
	fLabel_cand_N2[fN_cand_leadA]=i;
	fLabelM_cand_N2[fN_cand_leadA]=imoth;
	fMomCN2x[fN_cand_leadA] = px;
	fMomCN2y[fN_cand_leadA] = py;
	fMomCN2z[fN_cand_leadA] = pz;
	fN_cand_leadA++;
      }
    }
    else{ // negative Z (C side)
      if(charge && fP_cand_leadC < fgkDim){ // charged candidate
	fE_p_cand_leadC[fP_cand_leadC]=part->Energy(); 
	fPdg_cand_leadP1[fP_cand_leadC]=part->GetPdgCode(); 
	fPdgM_cand_leadP1[fP_cand_leadC]=pdgM;
	fLabel_cand_P1[fP_cand_leadC]=i; 
	fLabelM_cand_P1[fP_cand_leadC]=imoth;
	fMomCP1x[fP_cand_leadC] = px;
	fMomCP1y[fP_cand_leadC] = py;
	fMomCP1z[fP_cand_leadC] = pz;
	fP_cand_leadC++;
      }
      else if(!charge && fN_cand_leadC < fgkDim){ // neutral candidate
	fE_n_cand_leadC[fN_cand_leadC]=part->Energy(); 
	fPdg_cand_leadN1[fN_cand_leadC]=part->GetPdgCode(); 
	fPdgM_cand_leadN1[fN_cand_leadC]=pdgM;
	fLabel_cand_N1[fN_cand_leadC]=i; 
	fLabelM_cand_N1[fN_cand_leadC]=imoth;
	fMomCN1x[fN_cand_leadC] = px;
	fMomCN1y[fN_cand_leadC] = py;
	fMomCN1z[fN_cand_leadC] = pz;
	fN_cand_leadC++;
      }
    }
  } 
}
//________________________________________________________________________
void AliTaskLeadingMC::loopTrackRef(TTree *treeTR, AliMCEvent *mcEvent){
  // reset variables
  fP2hits=0;
  fP1hits=0;
  fN2hits=0;
  fN1hits=0;

  Int_t nevTR = treeTR->GetEntries();
  Int_t imoth;
  Float_t px,py,pz,pt;
  
  for(Int_t j=0; j<nevTR; j++){
    treeTR->GetEvent(j);
    // z momentum
    fTRpz = treeTR->GetLeaf("TrackReferences.fPz")->GetValue();
    // label in MC stack
    Int_t itr = Int_t(treeTR->GetLeaf("TrackReferences.fTrack")->GetValue());
    // hit position
    fX = treeTR->GetLeaf("TrackReferences.fX")->GetValue();
    fY = treeTR->GetLeaf("TrackReferences.fY")->GetValue();
    fZ = treeTR->GetLeaf("TrackReferences.fZ")->GetValue();
    
    // get particles from stack
    AliMCParticle *partMC = (AliMCParticle *) mcEvent->GetTrack(itr);  
    TParticle *part = partMC->Particle();
    
    px = part->Px();
    py = part->Py();
    pz = part->Pz();
    imoth = part->GetFirstMother();

    /*
    // recursive loop over mothers (to be checked)
    AliMCParticle *partMCM;
    if(imoth>7){
      partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);  
      TParticle *partM = partMCM->Particle();
      
      Int_t status = partM->GetStatusCode();
      
      if(status==11 &&  TMath::Abs(partM->GetPdgCode())>99) status = 1;
      if(status==1&&partM->GetFirstMother()>7){
	AliMCParticle *MCpartMM = (AliMCParticle *) mcEvent->GetTrack(partM->GetFirstMother());
	TParticle *partMM = MCpartMM->Particle(); 
	if(partMM->GetStatusCode()==11 &&  TMath::Abs(partMM->GetPdgCode())>99) status = 0;
      }
      
      while(partM && status != 1 && partM->GetFirstMother() != -1 && partM->GetFirstMother() > 7){
	imoth = partM->GetFirstMother();
	partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);
	partM = partMCM->Particle();
	
	status = partM->GetStatusCode();
	if(status==11 &&  TMath::Abs(partM->GetPdgCode())>99) status = 1;
	if(status==1&&partM->GetFirstMother()>7){
	  AliMCParticle *MCpartMM = (AliMCParticle *) mcEvent->GetTrack(partM->GetFirstMother());
	  TParticle *partMM = MCpartMM->Particle(); 
	  if(partMM->GetStatusCode()==11 &&  TMath::Abs(partMM->GetPdgCode())>99) status = 0;
	}
      }
    }
    else imoth = -1;
    */
    
    // check if hit in ZDC acceptance
    if(fY> fZDCP_Ymin  && fY < fZDCP_Ymax && fX > fZDCP_Xmin && fX < fZDCP_Xmax && TMath::Abs(fZ) > fZDCP_Zmin && TMath::Abs(fZ) <  fZDCP_Zmax){ // proton caloriemters
      if(fZ > 0.0 && fP2hits < fgkDim) { // A side
	fLabel_leadP2[fP2hits] = itr;
	fE_p_leadA[fP2hits]=part->Energy();
	fPdg_leadingP2[fP2hits]=part->GetPdgCode();
	fZP2[fP2hits]=fZ;
	fXP2[fP2hits]=fX;
	fYP2[fP2hits]=fY;
	fMomP2x[fP2hits] = px;
	fMomP2y[fP2hits] = py;
	fMomP2z[fP2hits] = pz;
	fLabelMP2[fP2hits] = imoth;
	fP2hits++;
      }
      else if(fZ < 0.0 && fP1hits < fgkDim){ // C side
	fLabel_leadP1[fP1hits] = itr;
	fE_p_leadC[fP1hits]=part->Energy(); 
	fPdg_leadingP1[fP1hits]=part->GetPdgCode();
	fZP1[fP1hits]=fZ;
	fXP1[fP1hits]=fX;
	fYP1[fP1hits]=fY;
	fMomP1x[fP1hits] = px;
	fMomP1y[fP1hits] = py;
	fMomP1z[fP1hits] = pz;
	fLabelMP1[fP1hits] = imoth;
	fP1hits++;
      }
    }
    else if(fY>fZDCN_Ymin && fY < fZDCN_Ymax && fX > fZDCN_Xmin && fX < fZDCN_Xmax && TMath::Abs(fZ) > fZDCN_Zmin && TMath::Abs(fZ) < fZDCN_Zmax){ // neutron calorimeters
      if(fZ > 0.0 && fN2hits < fgkDim) { 
	fLabel_leadN2[fN2hits] = itr;	
	fE_n_leadA[fN2hits]=part->Energy();
	fPdg_leadingN2[fN2hits]=part->GetPdgCode();
	fZN2[fN2hits]=fZ;
	fXN2[fN2hits]=fX;
	fYN2[fN2hits]=fY;
	fMomN2x[fN2hits] = px;
	fMomN2y[fN2hits] = py;
	fMomN2z[fN2hits] = pz;
	fLabelMN2[fN2hits] = imoth;
	fN2hits++;
      }
      else if(fZ < 0.0 && fN1hits < fgkDim) { 
	fLabel_leadN1[fN1hits] = itr;
	fE_n_leadC[fN1hits]=part->Energy();
	fPdg_leadingN1[fN1hits]=part->GetPdgCode();
	fZN1[fN1hits]=fZ;
	fXN1[fN1hits]=fX;
	fYN1[fN1hits]=fY;
	fMomN1x[fN1hits] = px;
	fMomN1y[fN1hits] = py;
	fMomN1z[fN1hits] = pz;
	fLabelMN1[fN1hits] = imoth;
	fN1hits++;
      }
    }    
  }
}
//________________________________________________________________________
void AliTaskLeadingMC::loopTrack(AliMCEvent *mcEvent){
  // do nothing
}
//________________________________________________________________________
void AliTaskLeadingMC::fillZDCreco(){
  // Get reco ZDC info from ZDC
  const Double_t *aZDCP2 = fESD->GetESDZDC()->GetZP2TowerEnergy();
  fAdcZDCP2[0] = aZDCP2[0];
  fAdcZDCP2[1] = aZDCP2[1];
  fAdcZDCP2[2] = aZDCP2[2];
  fAdcZDCP2[3] = aZDCP2[3];
  fAdcZDCP2[4] = aZDCP2[4];
  
  const Double_t *aZDCP1 = fESD->GetESDZDC()->GetZP1TowerEnergy();
  fAdcZDCP1[0] = aZDCP1[0];
  fAdcZDCP1[1] = aZDCP1[1];
  fAdcZDCP1[2] = aZDCP1[2];
  fAdcZDCP1[3] = aZDCP1[3];
  fAdcZDCP1[4] = aZDCP1[4];
  
  const Double_t *aZDCN2 = fESD->GetESDZDC()->GetZN2TowerEnergy();
  fAdcZDCN2[0] = aZDCN2[0];
  fAdcZDCN2[1] = aZDCN2[1];
  fAdcZDCN2[2] = aZDCN2[2];
  fAdcZDCN2[3] = aZDCN2[3];
  fAdcZDCN2[4] = aZDCN2[4];
  
  const Double_t *aZDCN1 = fESD->GetESDZDC()->GetZN1TowerEnergy();
  fAdcZDCN1[0] = aZDCN1[0];
  fAdcZDCN1[1] = aZDCN1[1];
  fAdcZDCN1[2] = aZDCN1[2];
  fAdcZDCN1[3] = aZDCN1[3];
  fAdcZDCN1[4] = aZDCN1[4];
}
//________________________________________________________________________
void AliTaskLeadingMC::Terminate(Option_t *) 
{
  fTree = dynamic_cast<TTree*> (GetOutputData(0));
  if (!fTree) {
    Printf("ERROR: fTree not available");
    return;
  }

  system("touch ok.job");
}

//__________________________________________________________________________________________________
/// Computes event spherocity. From: AliPhysics/PWGLF/RESONANCES/AliRsnMiniAnalysisTask.cxx 
/// 
/// For AODs, the filter bit 5 is used to select tracks for the computation. 
/// For tESDs, the track filter set externally is used to select tracks. 
///
/// \return Spherocity value if at least 10 good tracks are used for the computation. If not, returns -10.
/// 
Double_t AliTaskLeadingMC::ComputeSpherocity()
{
  AliVEvent * evTypeS = fESD;
  Int_t ntracksLoop = evTypeS->GetNumberOfTracks();
  fNTracksSpherocity = 0;
  Float_t spherocity = -10.0;
  Float_t pFull = 0;
  Float_t Spherocity = 2;
  Float_t pt[10000],phi[1000];
  
  //computing total pt
  Float_t sumapt = 0;
  for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
    AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
    AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
    AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
    if (aodt) if (!aodt->TestFilterBit(5)) continue;
    if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
    if (track->Pt() < 0.15) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    //pt[i1] = track->Pt();
    pt[i1] = 1.0;
    sumapt += pt[i1];
    fNTracksSpherocity++;
  }
  if (fNTracksSpherocity < 3) return -10.0;
  //Getting thrust
  for(Int_t i = 0; i < 360/0.1; ++i){
	Float_t numerador = 0;
	Float_t phiparam  = 0;
	Float_t nx = 0;
	Float_t ny = 0;
	phiparam=( (TMath::Pi()) * i * 0.1 ) / 180; // parametrization of the angle
	nx = TMath::Cos(phiparam);            // x component of an unitary vector n
	ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
	  AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
	  AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
	  AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
	  if (aodt) if (!aodt->TestFilterBit(5)) continue;
	  if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
	  if (track->Pt() < 0.15) continue;
	  if(TMath::Abs(track->Eta()) > 0.8) continue;
	  //pt[i1] = track->Pt();
	  pt[i1] = 1.0;
	  phi[i1] = track->Phi();
	  Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
	  Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );
	  numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
	}
	pFull=TMath::Power( (numerador / sumapt),2 );
	if(pFull < Spherocity)//maximization of pFull
	  {
	    Spherocity = pFull;
	  }
  }
  spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
  if (fNTracksSpherocity > 2) return spherocity;
  else return -10.0;
}

//----------------------------------------------------------------------------------
/// Defines track cuts for Spherocity computation. From: AliPhysics/PWGLF/RESONANCES/AliRsnMiniAnalysisTask.cxx 
///
/// Uses a track filter to set the quality cuts for the spherocity calculation.
/// Implementation only used for ESDs. AOD filtering is done via the filter bit. 
///
/// \param fTrackFilter Pointer to the AliAnalysisFilter 
///
void AliTaskLeadingMC::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
	//TPC Only
	esdTrackCuts->SetMinNClustersTPC(50);
	esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetMaxDCAToVertexZ(3.2);
	esdTrackCuts->SetMaxDCAToVertexXY(2.4);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
	
	esdTrackCuts->SetRequireTPCRefit(kTRUE);// TPC Refit
	esdTrackCuts->SetRequireITSRefit(kTRUE);// ITS Refit
	fTrackFilter->AddCuts(esdTrackCuts);
   return;
}

