//
// Jet trigger QA analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"
#include "AliEmcalParticle.h"
#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliVCaloCells.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskEmcalJetTriggerQA.h"

ClassImp(AliAnalysisTaskEmcalJetTriggerQA)

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::AliAnalysisTaskEmcalJetTriggerQA() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetTriggerQA", kTRUE),
  fDebug(kFALSE),
  fUseAnaUtils(kTRUE),
  fAnalysisUtils(0),
  fTriggerClass(""),
  fBitJ1((1<<8)),
  fBitJ2((1<<11)),
  fContainerFull(0),
  fContainerCharged(1),
  fMaxPatchEnergy(0),
  fTriggerType(-1),
  fNFastOR(16),
  fhNEvents(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fh2NJetsPtFull(0),
  fh2NJetsPtCharged(0),
  fh3PtEtaAreaJetFull(0),
  fh3PtEtaAreaJetCharged(0),
  fh2PtNConstituentsCharged(0),
  fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0),
  fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0),
  fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtzCharged(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0),
  fh3PtLeadJet2VsPatchEnergy(0),
  fh3PatchEnergyEtaPhiCenter(0),
  fh2CellEnergyVsTime(0)
{
  // Default constructor.

  
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::AliAnalysisTaskEmcalJetTriggerQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(kFALSE),
  fUseAnaUtils(kTRUE),
  fAnalysisUtils(0),
  fTriggerClass(""),
  fBitJ1((1<<8)),
  fBitJ2((1<<11)),
  fContainerFull(0),
  fContainerCharged(1),
  fMaxPatchEnergy(0),
  fTriggerType(-1),
  fNFastOR(16),
  fhNEvents(0),
  fh3PtEtaPhiJetFull(0),
  fh3PtEtaPhiJetCharged(0),
  fh2NJetsPtFull(0),
  fh2NJetsPtCharged(0),
  fh3PtEtaAreaJetFull(0),
  fh3PtEtaAreaJetCharged(0),
  fh2PtNConstituentsCharged(0),
  fh2PtNConstituents(0),
  fh2PtMeanPtConstituentsCharged(0),
  fh2PtMeanPtConstituentsNeutral(0),
  fh2PtNEF(0),
  fh3NEFEtaPhi(0),
  fh2NEFNConstituentsCharged(0),
  fh2NEFNConstituentsNeutral(0),
  fh2Ptz(0),
  fh2PtzCharged(0),
  fh2PtLeadJet1VsLeadJet2(0),
  fh3EEtaPhiCluster(0),
  fh3PtLeadJet1VsPatchEnergy(0),
  fh3PtLeadJet2VsPatchEnergy(0),
  fh3PatchEnergyEtaPhiCenter(0),
  fh2CellEnergyVsTime(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetTriggerQA::~AliAnalysisTaskEmcalJetTriggerQA()
{
  // Destructor.
}

//----------------------------------------------------------------------------------------------
void AliAnalysisTaskEmcalJetTriggerQA::InitOnce() {
  //
  // only initialize once
  //

  // Initialize analysis util class for vertex selection
  if(fUseAnaUtils) {
    fAnalysisUtils = new AliAnalysisUtils();
    fAnalysisUtils->SetMinVtxContr(2);
    fAnalysisUtils->SetMaxVtxZ(10.);
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::SelectEvent() {
  //
  // Decide if event should be selected for analysis
  //

  fhNEvents->Fill(0.5);
  
  if(fAnalysisUtils) {
    if(!fAnalysisUtils->IsVertexSelected2013pA(InputEvent()))
      return kFALSE;

    fhNEvents->Fill(2.5);

    if(fAnalysisUtils->IsPileUpEvent(InputEvent()))
      return kFALSE;
  }
  else{
    if(fUseAnaUtils)
      AliError("fAnalysisUtils not initialized. Call AliAnalysisTaskEmcalJetTriggerQA::InitOnce()");
  }

  fhNEvents->Fill(3.5);

  if(!fTriggerClass.IsNull()) {
    //Check if requested trigger was fired
    TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();
    if(!firedTrigClass.Contains(fTriggerClass))
      return kFALSE;
  }

  fhNEvents->Fill(1.5);

  return kTRUE;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::FindTriggerPatch() {

  //Code to get position of trigger

  if(!fGeom)
    fGeom = AliEMCALGeometry::GetInstance();


  TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();

  AliAODCaloTrigger *trg = dynamic_cast<AliAODCaloTrigger*>(InputEvent()->GetCaloTrigger("EMCAL"));
  trg->Reset();

  int col, row; //FASTOR position

  Int_t nPatchNotEmpty = 0; //counter number of patches which are not empty
  fMaxPatchEnergy = 0.;
  while (trg->Next())
    {

      trg->GetPosition(col, row); //col (0 to 63), row (0 to 47)
      
      if (col > -1 && row > -1)
	{
	  Int_t id = -1; //FASTOR index
	  Int_t cellIndex[4] = {-1};
	  fGeom->GetAbsFastORIndexFromPositionInEMCAL(col,row,id); //phi is column, eta is row
	  if(id<0) continue;

	  fGeom->GetCellIndexFromFastORIndex(id,cellIndex);

	  Int_t ts = 0;
	  trg->GetL1TimeSum(ts);

	  Float_t ampTrg = 0.;
	  trg->GetAmplitude(ampTrg);

	  //L1 analysis
	  Int_t bit = 0;
	  trg->GetTriggerBits(bit);

	  Bool_t bTrigJ = TestFilterBit(bit,fBitJ1|fBitJ2);
	  if(!bTrigJ)
	    continue;

	  Bool_t bTrigJ1 = TestFilterBit(bit,fBitJ1);	   
	  Bool_t bTrigJ2 = TestFilterBit(bit,fBitJ2);	   

	  if(bTrigJ1) fTriggerType = 0;
	  if(bTrigJ2) fTriggerType = 1;
	  if(!bTrigJ1 && !bTrigJ2) fTriggerType = -1;

	  Double_t minPhiPatch =  10.;
	  Double_t maxPhiPatch = -10.;
	  Double_t minEtaPatch =  10.;
	  Double_t maxEtaPatch = -10.;

	  //Get energy in trigger patch 8x8 FASTOR
	  Double_t patchEnergy = 0.;
	  Double_t sumAmp = 0.;
	  //	  const Int_t nFastOR = 8;//16;
	  for(Int_t fastrow = 0; fastrow<fNFastOR; fastrow++) {
	    for(Int_t fastcol = 0; fastcol<fNFastOR; fastcol++) {
	      Int_t nrow = row+fastrow;
	      Int_t ncol = col+fastcol;
	      fGeom->GetAbsFastORIndexFromPositionInEMCAL(ncol,nrow,id);

	      if(id<0) {
		AliWarning(Form("%s: id smaller than 0 %d",GetName(),id));
		continue;
	      }
	      
	      fGeom->GetCellIndexFromFastORIndex(id,cellIndex);
	      for(Int_t icell=0; icell<4; icell++) {

		if(!fGeom->CheckAbsCellId(cellIndex[icell])) continue;

		Double_t amp =0., time = 0., efrac = 0;
		Int_t mclabel = -1;
		Short_t absId = -1;
		Int_t nSupMod = -1, nModule = -1, nIphi = -1, nIeta = -1;
 
		fGeom->GetCellIndex(cellIndex[icell], nSupMod, nModule, nIphi, nIeta );
		
		fCaloCells->GetCell(cellIndex[icell], absId, amp, time,mclabel,efrac);

		Double_t eta, phi;
		fGeom->EtaPhiFromIndex(cellIndex[icell], eta, phi);

		if(phi<minPhiPatch) minPhiPatch = phi;
		if(phi>maxPhiPatch) maxPhiPatch = phi;
		if(eta<minEtaPatch) minEtaPatch = eta;
		if(eta>maxEtaPatch) maxEtaPatch = eta;

		sumAmp+=fCaloCells->GetAmplitude(cellIndex[icell]);
		patchEnergy+=amp;

	      }//cells in fastor loop
	    }//fastor col loop
	  }//fastor row loop

	  Double_t etaCent = minEtaPatch + 0.5*(maxEtaPatch-minEtaPatch);
	  Double_t phiCent = minPhiPatch + 0.5*(maxPhiPatch-minPhiPatch);
	  fh3PatchEnergyEtaPhiCenter->Fill(patchEnergy,etaCent,phiCent);

	  if(patchEnergy>0)
	    if(fDebug>2) AliInfo(Form("%s: patch edges eta: %f - %f   phi: %f - %f ampTrg: %f patchEnergy: %f sumAmp: %f",GetName(),minEtaPatch,maxEtaPatch,minPhiPatch,maxPhiPatch,ampTrg,patchEnergy,sumAmp));

	  if(patchEnergy>0) nPatchNotEmpty++;
	  
	  if(patchEnergy>fMaxPatchEnergy) fMaxPatchEnergy=patchEnergy;
 
	}
    }

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::LoadExtraBranches() {
  //
  // load extra brances
  //

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::UserCreateOutputObjects()
{
  // Create user output.

  InitOnce();

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhNEvents = new TH1F("fhNEvents","fhNEvents;selection;N_{evt}",5,0,5);
  fOutput->Add(fhNEvents);

  Int_t nBinsPt = 120;
  Double_t minPt = -20.;
  Double_t maxPt = 100.;
  Int_t nBinsEta = 100;
  Double_t minEta = -1.;
  Double_t maxEta = 1.;
  Int_t nBinsPhi = 18*8;
  Double_t minPhi = 0.;
  Double_t maxPhi = TMath::TwoPi();
  Int_t nBinsArea = 100;
  Double_t minArea = 0.;
  Double_t maxArea = 1.;

  fh3PtEtaPhiJetFull = new TH3F("fh3PtEtaPhiJetFull","fh3PtEtaPhiJetFull;#it{p}_{T}^{jet};#eta;#varphi",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3PtEtaPhiJetFull);

  fh3PtEtaPhiJetCharged = new TH3F("fh3PtEtaPhiJetCharged","fh3PtEtaPhiJetCharged;#it{p}_{T}^{jet};#eta;#varphi",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3PtEtaPhiJetCharged);

  fh2NJetsPtFull = new TH2F("fh2NJetsPtFull","fh2NJetsPtFull;N_{jets};#it{p}_{T}^{jet}",20,-0.5,19.5,nBinsPt,minPt,maxPt);
  fOutput->Add(fh2NJetsPtFull);

  fh2NJetsPtCharged = new TH2F("fh2NJetsPtCharged","fh2NJetsPtCharged;N_{jets};#it{p}_{T}^{jet}",20,-0.5,19.5,nBinsPt,minPt,maxPt);
  fOutput->Add(fh2NJetsPtCharged);

  fh3PtEtaAreaJetFull = new TH3F("fh3PtEtaAreaJetFull","fh3PtEtaAreaJetFull;#it{p}_{T}^{jet};#eta;A",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsArea,minArea,maxArea);
  fOutput->Add(fh3PtEtaAreaJetFull);

  fh3PtEtaAreaJetCharged = new TH3F("fh3PtEtaAreaJetCharged","fh3PtEtaAreaJetCharged;#it{p}_{T}^{jet};#eta;A",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsArea,minArea,maxArea);
  fOutput->Add(fh3PtEtaAreaJetCharged);

  Int_t nBinsConst =100;
  Double_t minConst = 0.;
  Double_t maxConst = 100.;

  Int_t nBinsMeanPt = 200;
  Double_t minMeanPt = 0.;
  Double_t maxMeanPt = 20.;

  Int_t nBinsNEF = 101;
  Double_t minNEF = 0.;
  Double_t maxNEF = 1.01;

  Int_t nBinsz = 101;
  Double_t minz = 0.;
  Double_t maxz = 1.01;

  Int_t nBinsECluster =100;
  Double_t minECluster = 0.;
  Double_t maxECluster = 100.;


  fh2PtNConstituentsCharged = new TH2F("fh2PtNConstituentsCharged","fh2PtNConstituentsCharged;#it{p}_{T}^{jet};N_{charged constituents}",nBinsPt,minPt,maxPt,nBinsConst,minConst,maxConst);
  fOutput->Add(fh2PtNConstituentsCharged);

  fh2PtNConstituents = new TH2F("fh2PtNConstituents","fh2PtNConstituents;#it{p}_{T}^{jet};N_{constituents}",nBinsPt,minPt,maxPt,nBinsConst,minConst,maxConst);
  fOutput->Add(fh2PtNConstituents);

  fh2PtMeanPtConstituentsCharged = new TH2F("fh2PtMeanPtConstituentsCharged","fh2PtMeanPtConstituentsCharged;#it{p}_{T}^{jet};charged #langle #it{p}_{T} #rangle",nBinsPt,minPt,maxPt,nBinsMeanPt,minMeanPt,maxMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsCharged);

  fh2PtMeanPtConstituentsNeutral = new TH2F("fh2PtMeanPtConstituentsNeutral","fh2PtMeanPtConstituentsNeutral;#it{p}_{T}^{jet};neutral langle #it{p}_{T} #rangle",nBinsPt,minPt,maxPt,nBinsMeanPt,minMeanPt,maxMeanPt);
  fOutput->Add(fh2PtMeanPtConstituentsNeutral);

  fh2PtNEF = new TH2F("fh2PtNEF","fh2PtNEF;#it{p}_{T}^{jet};NEF",nBinsPt,minPt,maxPt,nBinsNEF,minNEF,maxNEF);
  fOutput->Add(fh2PtNEF);

  fh3NEFEtaPhi = new TH3F("fh3NEFEtaPhi","fh3NEFEtaPhi;NEF;#eta;#varphi",nBinsNEF,minNEF,maxNEF,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3NEFEtaPhi);

  fh2NEFNConstituentsCharged = new TH2F("fh2NEFNConstituentsCharged","fh2NEFNConstituentsCharged;NEF;N_{charged constituents}",nBinsNEF,minNEF,maxNEF,nBinsConst,minConst,maxConst);
  fOutput->Add(fh2NEFNConstituentsCharged);

  fh2NEFNConstituentsNeutral = new TH2F("fh2NEFNConstituentsNeutral","fh2NEFNConstituentsNeutral;NEF;N_{clusters}",nBinsNEF,minNEF,maxNEF,nBinsConst,minConst,maxConst);
  fOutput->Add(fh2NEFNConstituentsNeutral);

  fh2Ptz = new TH2F("fh2Ptz","fh2Ptz;#it{p}_{T}^{jet};z=p_{t,trk}^{proj}/p_{jet}",nBinsPt,minPt,maxPt,nBinsz,minz,maxz);
  fOutput->Add(fh2Ptz);

  fh2PtzCharged = new TH2F("fh2PtzCharged","fh2Ptz;#it{p}_{T}^{ch jet};z=p_{t,trk}^{proj}/p_{ch jet}",nBinsPt,minPt,maxPt,nBinsz,minz,maxz);
  fOutput->Add(fh2PtzCharged);

  fh2PtLeadJet1VsLeadJet2 = new TH2F("fh2PtLeadJet1VsLeadJet2","fh2PtLeadJet1VsLeadJet2;#it{p}_{T}^{jet 1};#it{p}_{T}^{jet 2}",nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
  fOutput->Add(fh2PtLeadJet1VsLeadJet2);

  fh3EEtaPhiCluster = new TH3F("fh3EEtaPhiCluster","fh3EEtaPhiCluster;E_{clus};#eta;#phi",nBinsECluster,minECluster,maxECluster,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3EEtaPhiCluster);

  fh3PtLeadJet1VsPatchEnergy = new TH3F("fh3PtLeadJet1VsPatchEnergy","fh3PtLeadJet1VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,2,-0.5,1.5);
  fOutput->Add(fh3PtLeadJet1VsPatchEnergy);
  fh3PtLeadJet2VsPatchEnergy = new TH3F("fh3PtLeadJet2VsPatchEnergy","fh3PtLeadJet2VsPatchEnergy;#it{p}_{T}^{jet 1};Amplitude_{patch};trig type",nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,2,-0.5,1.5);
  fOutput->Add(fh3PtLeadJet2VsPatchEnergy);

  fh3PatchEnergyEtaPhiCenter = new TH3F("fh3PatchEnergyEtaPhiCenter","fh3PatchEnergyEtaPhiCenter;E_{patch};#eta;#phi",nBinsPt,minPt,maxPt,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi);
  fOutput->Add(fh3PatchEnergyEtaPhiCenter);

  fh2CellEnergyVsTime = new TH2F("fh2CellEnergyVsTime","fh2CellEnergyVsTime;E_{cell};time",100,0.,100.,700,-400,1000);
  fOutput->Add(fh2CellEnergyVsTime);


  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
    if (h2){
      h2->Sumw2();
      continue;
    }
    TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
    if (h3){
      h3->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::FillHistograms()
{
  // Fill histograms.
  
  AliClusterContainer  *clusCont = GetClusterContainer(0);

  if (clusCont) {
    Int_t nclusters = clusCont->GetNClusters();
    for (Int_t ic = 0; ic < nclusters; ic++) {
      AliVCluster *cluster = static_cast<AliVCluster*>(clusCont->GetCluster(ic));
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", ic));
	continue;
      }
      if (!cluster->IsEMCAL()) {
	AliDebug(11,Form("%s: Cluster is not emcal",GetName()));
	continue;
      }

      TLorentzVector lp;
      cluster->GetMomentum(lp, const_cast<Double_t*>(fVertex));

      //Fill eta,phi,E of clusters here
      fh3EEtaPhiCluster->Fill(lp.E(),lp.Eta(),lp.Phi());
    }
  }

  if(fCaloCells) {
    const Short_t nCells   = fCaloCells->GetNumberOfCells();

    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = fCaloCells->GetCellNumber(iCell);
      Double_t cellE = fCaloCells->GetCellAmplitude(cellId);
      Double_t cellT = fCaloCells->GetCellTime(cellId);

      AliDebug(2,Form("cell energy = %f  time = %f",cellE,cellT*1e9));
      fh2CellEnergyVsTime->Fill(cellE,cellT*1e9);
    
    }

  }

  Double_t ptLeadJet1 = 0.;
  Double_t ptLeadJet2 = 0.;

  TArrayI *nJetsArr = new TArrayI(fh2NJetsPtFull->GetNbinsY()+1);
  nJetsArr->Reset(0);
  nJetsArr->Set(fh2NJetsPtFull->GetNbinsY()+1);

  if (GetJetContainer(fContainerFull)) {
    const Int_t njets = GetNJets(fContainerFull);
    for (Int_t ij = 0; ij < njets; ij++) {

      AliEmcalJet* jet = GetAcceptJetFromArray(ij,fContainerFull);
      if (!jet)
	continue; //jet not selected
      
      Double_t jetPt = jet->Pt();
      if(jetPt>ptLeadJet1) ptLeadJet1=jetPt;
      fh3PtEtaPhiJetFull->Fill(jetPt,jet->Eta(),jet->Phi());
      fh3PtEtaAreaJetFull->Fill(jetPt,jet->Eta(),jet->Area());
      
      //count jets above certain pT threshold
      Int_t ptbin = fh2NJetsPtFull->GetYaxis()->FindBin(jetPt);
      for(Int_t iptbin = ptbin; iptbin<=fh2NJetsPtFull->GetNbinsY(); iptbin++)
	nJetsArr->AddAt(nJetsArr->At(iptbin)+1,iptbin);
      
      fh2PtNConstituentsCharged->Fill(jetPt,jet->GetNumberOfTracks());
      fh2PtNConstituents->Fill(jetPt,jet->GetNumberOfConstituents());

      fh2PtNEF->Fill(jetPt,jet->NEF());
      fh3NEFEtaPhi->Fill(jet->NEF(),jet->Eta(),jet->Phi());
      fh2NEFNConstituentsCharged->Fill(jet->NEF(),jet->GetNumberOfTracks());
      fh2NEFNConstituentsNeutral->Fill(jet->NEF(),jet->GetNumberOfClusters());

      AliVParticle *vp;
      Double_t sumPtCh = 0.;
      for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
	vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
	if(!vp) continue;
	fh2Ptz->Fill(jetPt,GetZ(vp,jet));
	sumPtCh+=vp->Pt();
      }
      
      if(jet->GetNumberOfTracks()>0)
	fh2PtMeanPtConstituentsCharged->Fill(jetPt,sumPtCh/(double)(jet->GetNumberOfTracks()) );


      AliVCluster *vc = 0x0;
      Double_t sumPtNe = 0.;
      if (clusCont) {
	for(Int_t icc=0; icc<jet->GetNumberOfClusters(); icc++) {
	  vc = static_cast<AliVCluster*>(clusCont->GetCluster(icc));
	  if(!vc) continue;

	  TLorentzVector lp;
	  vc->GetMomentum(lp, const_cast<Double_t*>(fVertex));
	  sumPtNe+=lp.Pt();
	  
	}

	if(jet->GetNumberOfClusters()>0)
	  fh2PtMeanPtConstituentsNeutral->Fill(jetPt,sumPtNe/(double)(jet->GetNumberOfClusters()) );
      }
    } //full jet loop

    for(Int_t i=1; i<=fh2NJetsPtFull->GetNbinsY(); i++) {
      Int_t nJetsInEvent = nJetsArr->At(i);
      fh2NJetsPtFull->Fill(nJetsInEvent,fh2NJetsPtFull->GetYaxis()->GetBinCenter(i));
    }

  }

  //Reset array to zero to also count charged jets
  nJetsArr->Reset(0);
  
  //Loop over charged jets
  if (GetJetContainer(fContainerCharged)) {
    const Int_t njets = GetNJets(fContainerCharged);
    for (Int_t ij = 0; ij < njets; ij++) {

      AliEmcalJet* jet = GetAcceptJetFromArray(ij,fContainerCharged);
      if (!jet)
	continue; //jet not selected

      Double_t jetPt = jet->Pt();
      if(jetPt>ptLeadJet2) ptLeadJet2=jetPt;
      fh3PtEtaPhiJetCharged->Fill(jetPt,jet->Eta(),jet->Phi());
      fh3PtEtaAreaJetCharged->Fill(jetPt,jet->Eta(),jet->Area());

      AliVParticle *vp;
      for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
	vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
	if(!vp) continue;
	fh2PtzCharged->Fill(jetPt,GetZ(vp,jet));
      }
      
      //count jets above certain pT threshold
      Int_t ptbin = fh2NJetsPtCharged->GetYaxis()->FindBin(jetPt);
      for(Int_t iptbin = ptbin; iptbin<=fh2NJetsPtCharged->GetNbinsY(); iptbin++)
	nJetsArr->AddAt(nJetsArr->At(iptbin)+1,iptbin);
      
    }//ch jet loop
    for(Int_t i=1; i<=fh2NJetsPtCharged->GetNbinsY(); i++) {
      Int_t nJetsInEvent = nJetsArr->At(i);
      fh2NJetsPtCharged->Fill(nJetsInEvent,fh2NJetsPtCharged->GetYaxis()->GetBinCenter(i));
    }
  }

  if(GetJetContainer(fContainerFull) && GetJetContainer(fContainerCharged)) {
    fh2PtLeadJet1VsLeadJet2->Fill(ptLeadJet1,ptLeadJet2);
  }

  fh3PtLeadJet1VsPatchEnergy->Fill(ptLeadJet1,fMaxPatchEnergy,fTriggerType);
  fh3PtLeadJet2VsPatchEnergy->Fill(ptLeadJet2,fMaxPatchEnergy,fTriggerType);

  if(nJetsArr) delete nJetsArr;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetTriggerQA::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  //Check if event is selected (vertex & pile-up)
  if(!SelectEvent())
    return kFALSE;
  
  LoadExtraBranches();

  if(!fTriggerClass.IsNull())
    FindTriggerPatch();

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetTriggerQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetZ(const AliVParticle *trk, const AliEmcalJet *jet)          const
{  
  // Get Z of constituent trk

  return GetZ(trk->Px(),trk->Py(),trk->Pz(),jet->Px(),jet->Py(),jet->Pz());
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetTriggerQA::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const
{
  // 
  // Get the z of a constituent inside of a jet
  //

  Double_t pJetSq = jetPx*jetPx+jetPy*jetPy+jetPz*jetPz;

  if(pJetSq>0.)
    return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/pJetSq;
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f",GetName(), pJetSq)); 
    return 0;
  }
}
