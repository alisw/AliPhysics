//

//_________________________________________________________________________
//  Task for transverse energy studies
//
//*-- Author: Marcelo G. Munhoz (USP)
//_________________________________________________________________________

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"

#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliCentrality.h"
#include "AliTrackerBase.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliKFVertex.h"
#include "AliKFParticle.h"

#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTrack.h"

#include "AliAnalysisTaskCaloElectron.h"

// Electron analysis: M. G. Munhoz
// based on example of an analysis task creating a p_t spectrum


ClassImp(AliAnalysisTaskCaloElectron)

//________________________________________________________________________
AliAnalysisTaskCaloElectron::AliAnalysisTaskCaloElectron(const char *name, Bool_t isMc) : AliAnalysisTaskSE(name), 
fIsMc(isMc) ,fSelectTrigger(0), fTrigger(0), fCurrentRunNum(-1), fESD(0), fVertex(0), fStack(0), fPID(0), fGeoUt(0), fGeom(0), fRecoUtil(0), fCaloClusters(0), fCells(0), 
fVtxStat(0), fEta(0), fPt(0), fP(0), fDedx(0), fNSigma(0), fDcaXY(0), fDcaZ(0), fE(0), fEtaClu(0), fPhiClu(0), fRes(0), fInvMass(0),
fEtaMCIn(0), fPtMCIn(0), fPMCIn(0), fEMCIn(0), fvXMCIn(0), fvYMCIn(0), fvZMCIn(0), fIsPhyPrimIn(0), fMomIn(0),
fEtaMC(0), fPtMC(0), fPMC(0), fDedxMC(0), fNSigmaMC(0), fDcaXYMC(0), fDcaZMC(0), fEMC(0), fEtaCluMC(0), fPhiCluMC(0), fResMC(0), fEMCMatch(0), fPhiCluMCMatch(0), fEtaCluMCMatch(0), fResMCMatch(0), fInvMassMC(0), fIsPhyPrim(0), fMom(0), 
fOutputList(0), fTree(0), fTreeMCIn(0), fTreeMC(0),
fHistVertexZ(0), fHistVertexZAll(0), 
fHistEtaPt(0), fHistEtaPteTPC(0), fHistEtaPteTPCMatch(0), fHistEtaPteEMCAL(0), 
fHistEtaPtStackeMC(0), fHistEtaPteMC(0), fHistEtaPteTPCeMC(0), fHistEtaPteTPCMatcheMC(0), fHistEtaPteEMCALeMC(0), fHistEtaPteMCMatch(0), fHistEtaPteTPCeMCMatch(0),
fHistEtaPtStackePhMC(0), fHistEtaPtePhMC(0), fHistEtaPteTPCePhMC(0), fHistEtaPteTPCMatchePhMC(0), fHistEtaPteEMCALePhMC(0), fHistEtaPtePhMCMatch(0), fHistEtaPteTPCePhMCMatch(0),
fHistCluE(0), fHistCluEeTPCMatch(0), fHistCluEeEMCAL(0),
fHistCluEeMC(0), fHistCluEeTPCMatcheMC(0), fHistCluEeEMCALeMC(0),
fHistCluEePhMC(0), fHistCluEeTPCMatchePhMC(0), fHistCluEeEMCALePhMC(0),
fHistCluEtaPhi(0), fHistCluEtaPhieTPCMatch(0), fHistCluEtaPhieEMCAL(0),
fHistCluEtaPhieMC(0), fHistCluEtaPhieTPCMatcheMC(0), fHistCluEtaPhieEMCALeMC(0),
fHistCluEtaPhiePhMC(0), fHistCluEtaPhieTPCMatchePhMC(0), fHistCluEtaPhieEMCALePhMC(0),
fHistVtxXYeMC(0),
fHistdEdx(0), fHistdEdxeTPC(0), fHistdEdxeTPCMatch(0), fHistdEdxeEMCAL(0),
fHistdEdxeMC(0), fHistdEdxeTPCeMC(0), fHistdEdxeTPCMatcheMC(0), fHistdEdxeEMCALeMC(0), fHistdEdxeMCMatch(0), fHistdEdxeTPCeMCMatch(0),
fHistdEdxePhMC(0), fHistdEdxeTPCePhMC(0), fHistdEdxeTPCMatchePhMC(0), fHistdEdxeEMCALePhMC(0), fHistdEdxePhMCMatch(0), fHistdEdxeTPCePhMCMatch(0),
fHistDeltaRZeTPCMatch(0), 
fHistDeltaRZeTPCMatcheMC(0), fHistDeltaRZeMCMatch(0), fHistDeltaRZeTPCeMCMatch(0),
fHistDeltaRZeTPCMatchePhMC(0), fHistDeltaRZePhMCMatch(0), fHistDeltaRZeTPCePhMCMatch(0),
fHistResPteTPCMatch(0), 
fHistResPteTPCMatcheMC(0), fHistResPteMCMatch(0), fHistResPteTPCeMCMatch(0),
fHistResPteTPCMatchePhMC(0), fHistResPtePhMCMatch(0), fHistResPteTPCePhMCMatch(0),
fHistEoverPPteTPCMatch(0), fHistEoverPPteEMCAL(0), fHistEoverPPeTPCMatch(0), fHistEoverPPeEMCAL(0), fHistEoverPPeEMCALtrigger(0), fHistEoverPEeTPCMatch(0), fHistEoverPEeEMCAL(0), 
fHistEoverPPteTPCMatcheMC(0), fHistEoverPPteEMCALeMC(0), 
fHistEoverPPteMCMatch(0), fHistEoverPPteEMCALeMCMatch(0), fHistEoverPPteTPCeMCMatch(0), fHistEoverPPteTPCeEMCALeMCMatch(0),
fHistEoverPPteTPCMatchePhMC(0), fHistEoverPPteEMCALePhMC(0), 
fHistEoverPPtePhMCMatch(0), fHistEoverPPteEMCALePhMCMatch(0), fHistEoverPPteTPCePhMCMatch(0), fHistEoverPPteTPCeEMCALePhMCMatch(0),
fHistMassLikePteTPC(0), fHistMassUnLikePteTPC(0),
fHistMassLikePteEMCAL(0), fHistMassUnLikePteEMCAL(0), fHistMassLikePeEMCAL(0), fHistMassUnLikePeEMCAL(0), fHistMassLikeEeEMCAL(0), fHistMassUnLikeEeEMCAL(0),
fHistMassLikePteTPCePhMC(0), fHistMassUnLikePteTPCePhMC(0),
fHistMassLikePteEMCALePhMC(0), fHistMassUnLikePteEMCALePhMC(0), fHistMassLikePeEMCALePhMC(0), fHistMassUnLikePeEMCALePhMC(0), fHistMassLikeEeEMCALePhMC(0), fHistMassUnLikeEeEMCALePhMC(0),
fHistMassLikePteTPCeNonPhMC(0), fHistMassUnLikePteTPCeNonPhMC(0),
fHistMassLikePteEMCALeNonPhMC(0), fHistMassUnLikePteEMCALeNonPhMC(0), fHistMassLikePeEMCALeNonPhMC(0), fHistMassUnLikePeEMCALeNonPhMC(0), fHistMassLikeEeEMCALeNonPhMC(0), fHistMassUnLikeEeEMCALeNonPhMC(0)
{
	// Constructor
	fgMassElectron = 0.000510998910;
	fCutRMatch = 0.02;
	fCutEoverPMin = 0.8;
	fCutEoverPMax = 1.2;
	//fCutdEdxMin = 75; // real data
	fCutdEdxMin = 68; // MC
	fCutdEdxMax = 95;

	fEsdTrackCuts = GetTrackCuts();
	fEsdTrackPhCuts = GetTrackPhCuts();
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 id reserved by the base class for AOD
	// Output slot #1 writes into a TH1 container
	DefineOutput(1, TList::Class());
	
	TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
	//TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
	TGeoManager::Import("geometry.root");
}

AliAnalysisTaskCaloElectron::~AliAnalysisTaskCaloElectron() 
{
	//Destructor

	delete fESD;
	delete fGeoUt;
	delete fRecoUtil;
    delete fEsdTrackCuts;
    delete fEsdTrackPhCuts;
	
	fOutputList->Clear();
}

//________________________________________________________________________
void AliAnalysisTaskCaloElectron::UserCreateOutputObjects()
{
	// set EMCal geometry
	fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE1");
	InitCaloUtil();
	
	// Create histograms
	// Called once

	Int_t nBinsPt = 100;
	Float_t minPt = 0.;
	Float_t maxPt = 10.;
	
	Int_t nBinsDelta = 200;
	Float_t minDelta = -0.1;
	Float_t maxDelta = 0.1;	
		
	Int_t nBinsRes = 200;
	Float_t minRes = 0.;
	Float_t maxRes = 0.1;	
	
	Int_t nBinsEoverP = 500;
	Float_t minEoverP = 0.;
	Float_t maxEoverP = 10.;
	
	Int_t nBinsEta = 200;
	Float_t minEta = -1.;
	Float_t maxEta = 1.;	
	
	Int_t nBinsPhi = 240;
	Float_t minPhi = 70.;
	Float_t maxPhi = 190.;	
	
	Int_t nBinsMass = 1000;
	Float_t minMass = 0.;
	Float_t maxMass = 1.0;

	fOutputList = new TList();
	
	fTree = new TTree("fTree","electrons");
    fTree->Branch("vtxStatus",&fVtxStat,"fVtxStat/F");
    fTree->Branch("eta",&fEta,"fEta/F");
    fTree->Branch("pt",&fPt,"fPt/F");
    fTree->Branch("p",&fP,"fP/F");
    fTree->Branch("dedx",&fDedx,"fDedx/F");
    fTree->Branch("nSigma",&fNSigma,"fNSigma/F");
    fTree->Branch("dcaXY",&fDcaXY,"fDcaXY/F");
    fTree->Branch("dcaZ",&fDcaZ,"fDcaZ/F");
    fTree->Branch("E",&fE,"fE/F");
    fTree->Branch("etaClu",&fEtaClu,"fEtaClu/F");
    fTree->Branch("phiClu",&fPhiClu,"fPhiClu/F");
    fTree->Branch("res",&fRes,"fRes/F");
	fTree->Branch("invMass",&fInvMass,"fInvMass/F");

	fHistVertexZ = new TH2F("fHistVertexZ","Z vertex position",2,-0.5,1.5,100,-50.,50.);
	fHistVertexZAll = new TH2F("fHistVertexZAll","Z vertex position, all events",2,-0.5,1.5,100,-50.,50.);
	
	
	fHistEtaPt = new TH2F("fHistEtaPt","eta vs pT distibution",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
	fHistEtaPteTPC = new TH2F("fHistEtaPteTPC","electron eta vs pT distibution, TPC ID",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
	fHistEtaPteTPCMatch = new TH2F("fHistEtaPteTPCMatch","electron eta vs pT distibution, TPC ID",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
	fHistEtaPteEMCAL = new TH2F("fHistEtaPteEMCAL","eta vs pT distibution, EMCal ID",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);

	fHistdEdx = new TH2F("fHistdEdx","TPC dEdx vs P",nBinsPt,minPt,maxPt,100,0.,200.);
	fHistdEdxeTPC = new TH2F("fHistdEdxeTPC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
	fHistdEdxeTPCMatch = new TH2F("fHistdEdxeTPCMatch","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
	fHistdEdxeEMCAL = new TH2F("fHistdEdxeEMCAL","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);

	
	fHistCluE = new TH1F("fHistCluE","cluster E",nBinsPt,minPt,maxPt);
	fHistCluEeTPCMatch = new TH1F("fHistCluEeTPCMatch","cluster E, TPC id",nBinsPt,minPt,maxPt);
	fHistCluEeEMCAL = new TH1F("fHistCluEeEMCAL","cluster E, EMCAL id",nBinsPt,minPt,maxPt);


	fHistCluEtaPhi = new TH2F("fHistCluEtaPhi","cluster phi vs eta",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
	fHistCluEtaPhieTPCMatch = new TH2F("fHistCluEtaPhieTPCMatch","cluster phi vs eta, TPC id",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
	fHistCluEtaPhieEMCAL = new TH2F("fHistCluEtaPhieEMCAL","cluster phi vs eta, EMCAL id",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);

	
	fHistDeltaRZeTPCMatch = new TH2F("fHistDeltaRZeTPCMatch","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);

	
	fHistResPteTPCMatch = new TH2F("fHistResPteTPCMatch","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);	
	
	
	fHistEoverPPteTPCMatch = new TH2F("fHistEoverPPteTPCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	fHistEoverPPteEMCAL = new TH2F("fHistEoverPPteEMCAL","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	
	fHistEoverPPeTPCMatch = new TH2F("fHistEoverPPeTPCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	fHistEoverPPeEMCAL = new TH2F("fHistEoverPPeEMCAL","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	//fHistEoverPPeEMCALtrigger = new TH2F("fHistEoverPPeEMCALtrigger","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	
	fHistEoverPEeTPCMatch = new TH2F("fHistEoverPEeTPCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	fHistEoverPEeEMCAL = new TH2F("fHistEoverPEeEMCAL","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	
	
	fHistMassLikePteTPC = new TH2F("fHistMassLikePteTPC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassUnLikePteTPC = new TH2F("fHistMassUnLikePteTPC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassLikePteEMCAL = new TH2F("fHistMassLikePteEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassUnLikePteEMCAL = new TH2F("fHistMassUnLikePteEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassLikePeEMCAL = new TH2F("fHistMassLikePeEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassUnLikePeEMCAL = new TH2F("fHistMassUnLikePeEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassLikeEeEMCAL = new TH2F("fHistMassLikeEeEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	fHistMassUnLikeEeEMCAL = new TH2F("fHistMassUnLikeEeEMCAL","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	
	if (fIsMc)
	{
		fTreeMCIn = new TTree("fTreeMCIn","electrons MC input");
		fTreeMCIn->Branch("etaMCIn",&fEtaMCIn,"fEtaMCIn/F");
		fTreeMCIn->Branch("ptMCIn",&fPtMCIn,"fPtMCIn/F");
		fTreeMCIn->Branch("pMCIn",&fPMCIn,"fPMCIn/F");
		fTreeMCIn->Branch("EMCIn",&fEMCIn,"fEMCIn/F");
		fTreeMCIn->Branch("vXMCIn",&fvXMCIn,"fvXMCIn/F");
		fTreeMCIn->Branch("vYMCIn",&fvYMCIn,"fvYMCIn/F");
		fTreeMCIn->Branch("vZMCIn",&fvZMCIn,"fvZMCIn/F");
		fTreeMCIn->Branch("IsPhyPrimIn",&fIsPhyPrimIn,"fIsPhyPrimIn/F");		
		fTreeMCIn->Branch("motherIn",&fMomIn,"fMomIn/F");
		
		fTreeMC = new TTree("fTreeMC","electrons MC");
		fTreeMC->Branch("etaMC",&fEtaMC,"fEtaMC/F");
		fTreeMC->Branch("ptMC",&fPtMC,"fPtMC/F");
		fTreeMC->Branch("pMC",&fPMC,"fPMC/F");
		fTreeMC->Branch("dedxMC",&fDedxMC,"fDedxMC/F");
		fTreeMC->Branch("nSigmaMC",&fNSigmaMC,"fNSigmaMC/F");
		fTreeMC->Branch("dcaXYMC",&fDcaXYMC,"fDcaXYMC/F");
		fTreeMC->Branch("dcaZMC",&fDcaZMC,"fDcaZMC/F");		
		fTreeMC->Branch("EMC",&fEMC,"fEMC/F");
		fTreeMC->Branch("etaCluMC",&fEtaCluMC,"fEtaCluMC/F");
		fTreeMC->Branch("phiCluMC",&fPhiCluMC,"fPhiCluMC/F");
		fTreeMC->Branch("resMC",&fResMC,"fResMC/F");
		fTreeMC->Branch("EMCMatch",&fEMCMatch,"fEMCMatch/F");
		fTreeMC->Branch("etaCluMCMatch",&fEtaCluMCMatch,"fEtaCluMCMatch/F");
		fTreeMC->Branch("phiCluMCMatch",&fPhiCluMCMatch,"fPhiCluMCMatch/F");
		fTreeMC->Branch("resMCMatch",&fResMCMatch,"fResMCMatch/F");
		fTreeMC->Branch("invMassMC",&fInvMassMC,"fInvMassMC/F");
		fTreeMC->Branch("IsPhyPrim",&fIsPhyPrim,"fIsPhyPrim/F");		
		fTreeMC->Branch("mother",&fMom,"fMom/F");

		fHistEtaPtStackeMC = new TH2F("fHistEtaPtStackeMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteMC = new TH2F("fHistEtaPteMC","eta vs pT distibution, MC electron, good tracks",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCeMC = new TH2F("fHistEtaPteTPCeMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCMatcheMC = new TH2F("fHistEtaPteTPCMatcheMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteEMCALeMC = new TH2F("fHistEtaPteEMCALeMC","eta vs pT distibution, EMCal ID",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteMCMatch = new TH2F("fHistEtaPteMCMatch","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCeMCMatch = new TH2F("fHistEtaPteTPCeMCMatch","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);

		fHistEtaPtStackePhMC = new TH2F("fHistEtaPtStackePhMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPtePhMC = new TH2F("fHistEtaPtePhMC","eta vs pT distibution, MC electron, good tracks",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCePhMC = new TH2F("fHistEtaPteTPCePhMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCMatchePhMC = new TH2F("fHistEtaPteTPCMatchePhMC","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteEMCALePhMC = new TH2F("fHistEtaPteEMCALePhMC","eta vs pT distibution, EMCal ID",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPtePhMCMatch = new TH2F("fHistEtaPtePhMCMatch","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		fHistEtaPteTPCePhMCMatch = new TH2F("fHistEtaPteTPCePhMCMatch","eta vs pT distibution, MC electron",nBinsPt,minPt,maxPt,nBinsEta, minEta, maxEta);
		
		fHistCluEeMC = new TH1F("fHistCluEeMC","cluster E, MC electrons",nBinsPt,minPt,maxPt);
		fHistCluEeTPCMatcheMC = new TH1F("fHistCluEeTPCMatcheMC","cluster E, TPC id, MC electrons",nBinsPt,minPt,maxPt);
		fHistCluEeEMCALeMC = new TH1F("fHistCluEeEMCALeMC","cluster E, EMCAL id, MC electrons",nBinsPt,minPt,maxPt);
		
		fHistCluEePhMC = new TH1F("fHistCluEePhMC","cluster E, photonic MC electrons",nBinsPt,minPt,maxPt);
		fHistCluEeTPCMatchePhMC = new TH1F("fHistCluEeTPCMatchePhMC","cluster E, TPC id, photonic MC electrons",nBinsPt,minPt,maxPt);
		fHistCluEeEMCALePhMC = new TH1F("fHistCluEeEMCALePhMC","cluster E, EMCAL id, photonic MC electrons",nBinsPt,minPt,maxPt);
		
		fHistCluEtaPhieMC = new TH2F("fHistCluEtaPhieMC","cluster phi vs eta, MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
		fHistCluEtaPhieTPCMatcheMC = new TH2F("fHistCluEtaPhieTPCMatcheMC","cluster phi vs eta, TPC id, MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
		fHistCluEtaPhieEMCALeMC = new TH2F("fHistCluEtaPhieEMCALeMC","cluster phi vs eta, EMCAL id, MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
				
		fHistCluEtaPhiePhMC = new TH2F("fHistCluEtaPhiePhMC","cluster phi vs eta, photonic MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
		fHistCluEtaPhieTPCMatchePhMC = new TH2F("fHistCluEtaPhieTPCMatchePhMC","cluster phi vs eta, TPC id, photonic MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
		fHistCluEtaPhieEMCALePhMC = new TH2F("fHistCluEtaPhieEMCALePhMC","cluster phi vs eta, EMCAL id, photonic MC electrons",nBinsPhi,minPhi,maxPhi,nBinsEta,minEta,maxEta);
		
		//fHistVtxXYeMC = new TH2F("fHistVtxXYeMC","Electron Mother X,Y vertex position",1201,-600.5,600.5,1201,-600.5,600.5);
		
		fHistdEdxeMC = new TH2F("fHistdEdxeMC","TPC dEdx vs P, MC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCeMC = new TH2F("fHistdEdxeTPCeMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCMatcheMC = new TH2F("fHistdEdxeTPCMatcheMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeEMCALeMC = new TH2F("fHistdEdxeEMCALeMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeMCMatch = new TH2F("fHistdEdxeMCMatch","TPC dEdx vs P, MC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCeMCMatch = new TH2F("fHistdEdxeTPCeMCMatch","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		
		fHistdEdxePhMC = new TH2F("fHistdEdxePhMC","TPC dEdx vs P, MC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCePhMC = new TH2F("fHistdEdxeTPCePhMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCMatchePhMC = new TH2F("fHistdEdxeTPCMatchePhMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeEMCALePhMC = new TH2F("fHistdEdxeEMCALePhMC","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxePhMCMatch = new TH2F("fHistdEdxePhMCMatch","TPC dEdx vs P, MC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		fHistdEdxeTPCePhMCMatch = new TH2F("fHistdEdxeTPCePhMCMatch","TPC dEdx vs P, TPC electrons",nBinsPt,minPt,maxPt,100,0.,200.);
		
		fHistDeltaRZeTPCMatcheMC = new TH2F("fHistDeltaRZeTPCMatcheMC","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);		
		fHistDeltaRZeMCMatch = new TH2F("fHistDeltaRZeMCMatch","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);
		fHistDeltaRZeTPCeMCMatch = new TH2F("fHistDeltaRZeTPCeMCMatch","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);
		
		fHistDeltaRZeTPCMatchePhMC = new TH2F("fHistDeltaRZeTPCMatchePhMC","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);		
		fHistDeltaRZePhMCMatch = new TH2F("fHistDeltaRZePhMCMatch","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);
		fHistDeltaRZeTPCePhMCMatch = new TH2F("fHistDeltaRZeTPCePhMCMatch","#Delta#phi vs #Delta#eta (track projection - cluster position), electrons",nBinsDelta,minDelta,maxDelta,nBinsDelta,minDelta,maxDelta);
		
		fHistResPteTPCMatcheMC = new TH2F("fHistResPteTPCMatcheMC","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);		
		fHistResPteMCMatch = new TH2F("fHistResPteMCMatch","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);	
		fHistResPteTPCeMCMatch = new TH2F("fHistResPteTPCeMCMatch","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);	
		
		fHistResPteTPCMatchePhMC = new TH2F("fHistResPteTPCMatchePhMC","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);		
		fHistResPtePhMCMatch = new TH2F("fHistResPtePhMCMatch","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);	
		fHistResPteTPCePhMCMatch = new TH2F("fHistResPteTPCePhMCMatch","match radius vs pt, electrons",nBinsPt,minPt,maxPt,nBinsRes,minRes,maxRes);	
		
		fHistEoverPPteTPCMatcheMC = new TH2F("fHistEoverPPteTPCMatcheMC","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteEMCALeMC = new TH2F("fHistEoverPPteEMCALeMC","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);		
		fHistEoverPPteMCMatch = new TH2F("fHistEoverPPteMCMatch","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteEMCALeMCMatch = new TH2F("fHistEoverPPteEMCALeMCMatch","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteTPCeMCMatch = new TH2F("fHistEoverPPteTPCeMCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteTPCeEMCALeMCMatch = new TH2F("fHistEoverPPteTPCeEMCALeMCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
	
		fHistEoverPPteTPCMatchePhMC = new TH2F("fHistEoverPPteTPCMatchePhMC","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteEMCALePhMC = new TH2F("fHistEoverPPteEMCALePhMC","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);		
		fHistEoverPPtePhMCMatch = new TH2F("fHistEoverPPtePhMCMatch","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteEMCALePhMCMatch = new TH2F("fHistEoverPPteEMCALePhMCMatch","E over p, primary electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteTPCePhMCMatch = new TH2F("fHistEoverPPteTPCePhMCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		fHistEoverPPteTPCeEMCALePhMCMatch = new TH2F("fHistEoverPPteTPCeEMCALePhMCMatch","E over p, electrons",nBinsPt,minPt,maxPt,nBinsEoverP,minEoverP,maxEoverP);
		
		fHistMassLikePteTPCePhMC = new TH2F("fHistMassLikePteTPCePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePteTPCePhMC = new TH2F("fHistMassUnLikePteTPCePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikePteEMCALePhMC = new TH2F("fHistMassLikePteEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePteEMCALePhMC = new TH2F("fHistMassUnLikePteEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikePeEMCALePhMC = new TH2F("fHistMassLikePeEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePeEMCALePhMC = new TH2F("fHistMassUnLikePeEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikeEeEMCALePhMC = new TH2F("fHistMassLikeEeEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikeEeEMCALePhMC = new TH2F("fHistMassUnLikeEeEMCALePhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		
		fHistMassLikePteTPCeNonPhMC = new TH2F("fHistMassLikePteTPCeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePteTPCeNonPhMC = new TH2F("fHistMassUnLikePteTPCeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikePteEMCALeNonPhMC = new TH2F("fHistMassLikePteEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePteEMCALeNonPhMC = new TH2F("fHistMassUnLikePteEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikePeEMCALeNonPhMC = new TH2F("fHistMassLikePeEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikePeEMCALeNonPhMC = new TH2F("fHistMassUnLikePeEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassLikeEeEMCALeNonPhMC = new TH2F("fHistMassLikeEeEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
		fHistMassUnLikeEeEMCALeNonPhMC = new TH2F("fHistMassUnLikeEeEMCALeNonPhMC","Invariant mass of e+/e-",nBinsPt,minPt,maxPt,nBinsMass,minMass,maxMass);
	}

	fOutputList->Add(fTree);

	fOutputList->Add(fHistVertexZAll);
	fOutputList->Add(fHistVertexZ);
	
	fOutputList->Add(fHistEtaPt);
	fOutputList->Add(fHistEtaPteTPC);
	fOutputList->Add(fHistEtaPteTPCMatch);
	fOutputList->Add(fHistEtaPteEMCAL);

	fOutputList->Add(fHistdEdx);
	fOutputList->Add(fHistdEdxeTPC);
	fOutputList->Add(fHistdEdxeTPCMatch);
	fOutputList->Add(fHistdEdxeEMCAL);

	fOutputList->Add(fHistCluE);
	fOutputList->Add(fHistCluEeTPCMatch);
	fOutputList->Add(fHistCluEeEMCAL);

	fOutputList->Add(fHistCluEtaPhi);
	fOutputList->Add(fHistCluEtaPhieTPCMatch);
	fOutputList->Add(fHistCluEtaPhieEMCAL);
	
	fOutputList->Add(fHistDeltaRZeTPCMatch);
	
	fOutputList->Add(fHistResPteTPCMatch);

	fOutputList->Add(fHistEoverPPteTPCMatch);
	fOutputList->Add(fHistEoverPPteEMCAL);

	fOutputList->Add(fHistEoverPPeTPCMatch);
	fOutputList->Add(fHistEoverPPeEMCAL);
	//fOutputList->Add(fHistEoverPPeEMCALtrigger);
	
	fOutputList->Add(fHistEoverPEeTPCMatch);
	fOutputList->Add(fHistEoverPEeEMCAL);
	
	fOutputList->Add(fHistMassLikePteTPC);
	fOutputList->Add(fHistMassUnLikePteTPC);
	fOutputList->Add(fHistMassLikePteEMCAL);
	fOutputList->Add(fHistMassUnLikePteEMCAL);
	fOutputList->Add(fHistMassLikePeEMCAL);
	fOutputList->Add(fHistMassUnLikePeEMCAL);
	fOutputList->Add(fHistMassLikeEeEMCAL);
	fOutputList->Add(fHistMassUnLikeEeEMCAL);
	
	if (fIsMc)
	{
		fOutputList->Add(fTreeMCIn);
		fOutputList->Add(fTreeMC);

		fOutputList->Add(fHistEtaPtStackeMC);
		fOutputList->Add(fHistEtaPteMC);
		fOutputList->Add(fHistEtaPteTPCeMC);
		fOutputList->Add(fHistEtaPteTPCMatcheMC);
		fOutputList->Add(fHistEtaPteEMCALeMC);
		fOutputList->Add(fHistEtaPteMCMatch);
		fOutputList->Add(fHistEtaPteTPCeMCMatch);
	
		fOutputList->Add(fHistEtaPtStackePhMC);
		fOutputList->Add(fHistEtaPtePhMC);
		fOutputList->Add(fHistEtaPteTPCePhMC);
		fOutputList->Add(fHistEtaPteTPCMatchePhMC);
		fOutputList->Add(fHistEtaPteEMCALePhMC);
		fOutputList->Add(fHistEtaPtePhMCMatch);
		fOutputList->Add(fHistEtaPteTPCePhMCMatch);
		
		//fOutputList->Add(fHistVtxXYeMC);

		fOutputList->Add(fHistdEdxeMC);
		fOutputList->Add(fHistdEdxeTPCeMC);
		fOutputList->Add(fHistdEdxeTPCMatcheMC);
		fOutputList->Add(fHistdEdxeEMCALeMC);
		fOutputList->Add(fHistdEdxeMCMatch);
		fOutputList->Add(fHistdEdxeTPCeMCMatch);

		fOutputList->Add(fHistdEdxePhMC);
		fOutputList->Add(fHistdEdxeTPCePhMC);
		fOutputList->Add(fHistdEdxeTPCMatchePhMC);
		fOutputList->Add(fHistdEdxeEMCALePhMC);
		fOutputList->Add(fHistdEdxePhMCMatch);
		fOutputList->Add(fHistdEdxeTPCePhMCMatch);

		fOutputList->Add(fHistCluEeMC);
		fOutputList->Add(fHistCluEeTPCMatcheMC);
		fOutputList->Add(fHistCluEeEMCALeMC);

		fOutputList->Add(fHistCluEePhMC);
		fOutputList->Add(fHistCluEeTPCMatchePhMC);
		fOutputList->Add(fHistCluEeEMCALePhMC);
		
		fOutputList->Add(fHistCluEtaPhieMC);
		fOutputList->Add(fHistCluEtaPhieTPCMatcheMC);
		fOutputList->Add(fHistCluEtaPhieEMCALeMC);
		
		fOutputList->Add(fHistCluEtaPhiePhMC);
		fOutputList->Add(fHistCluEtaPhieTPCMatchePhMC);
		fOutputList->Add(fHistCluEtaPhieEMCALePhMC);
		
		fOutputList->Add(fHistDeltaRZeTPCMatcheMC);
		fOutputList->Add(fHistDeltaRZeMCMatch);
		fOutputList->Add(fHistDeltaRZeTPCeMCMatch);

		fOutputList->Add(fHistDeltaRZeTPCMatchePhMC);
		fOutputList->Add(fHistDeltaRZePhMCMatch);
		fOutputList->Add(fHistDeltaRZeTPCePhMCMatch);
		
		fOutputList->Add(fHistResPteTPCMatcheMC);
		fOutputList->Add(fHistResPteMCMatch);
		fOutputList->Add(fHistResPteTPCeMCMatch);
		
		fOutputList->Add(fHistResPteTPCMatchePhMC);
		fOutputList->Add(fHistResPtePhMCMatch);
		fOutputList->Add(fHistResPteTPCePhMCMatch);
		
		fOutputList->Add(fHistEoverPPteTPCMatcheMC);
		fOutputList->Add(fHistEoverPPteEMCALeMC);
		fOutputList->Add(fHistEoverPPteMCMatch);
		fOutputList->Add(fHistEoverPPteEMCALeMCMatch);
		fOutputList->Add(fHistEoverPPteTPCeMCMatch);
		fOutputList->Add(fHistEoverPPteTPCeEMCALeMCMatch);
		
		fOutputList->Add(fHistEoverPPteTPCMatchePhMC);
		fOutputList->Add(fHistEoverPPteEMCALePhMC);
		fOutputList->Add(fHistEoverPPtePhMCMatch);
		fOutputList->Add(fHistEoverPPteEMCALePhMCMatch);
		fOutputList->Add(fHistEoverPPteTPCePhMCMatch);
		fOutputList->Add(fHistEoverPPteTPCeEMCALePhMCMatch);
		
		fOutputList->Add(fHistMassLikePteTPCePhMC);
		fOutputList->Add(fHistMassUnLikePteTPCePhMC);
		fOutputList->Add(fHistMassLikePteEMCALePhMC);
		fOutputList->Add(fHistMassUnLikePteEMCALePhMC);
		fOutputList->Add(fHistMassLikePeEMCALePhMC);
		fOutputList->Add(fHistMassUnLikePeEMCALePhMC);
		fOutputList->Add(fHistMassLikeEeEMCALePhMC);
		fOutputList->Add(fHistMassUnLikeEeEMCALePhMC);

		fOutputList->Add(fHistMassLikePteTPCeNonPhMC);
		fOutputList->Add(fHistMassUnLikePteTPCeNonPhMC);
		fOutputList->Add(fHistMassLikePteEMCALeNonPhMC);
		fOutputList->Add(fHistMassUnLikePteEMCALeNonPhMC);
		fOutputList->Add(fHistMassLikePeEMCALeNonPhMC);
		fOutputList->Add(fHistMassUnLikePeEMCALeNonPhMC);
		fOutputList->Add(fHistMassLikeEeEMCALeNonPhMC);
		fOutputList->Add(fHistMassUnLikeEeEMCALeNonPhMC);		
	}		
}


//________________________________________________________________________
void AliAnalysisTaskCaloElectron::UserExec(Option_t *) 
{
	// Main loop
	// Called for each event
	
	// Post output data.
	fESD = dynamic_cast<AliESDEvent*>(InputEvent());
	if (!fESD) {
		printf("ERROR: ESD not available\n");
		return;
	}
	
	if (!IsPhysicsSelected())
		return;
		
	fVertex = (AliESDVertex*)fESD->GetPrimaryVertex();
	if (!fVertex) {
		printf("ERROR: vertex not available\n");
		return;
	}
	
	Double_t zVertex = fVertex->GetZ();
	fVtxStat = fVertex->GetStatus();

	fHistVertexZAll->Fill(fVtxStat,zVertex);

	if (!IsGoodEvent())
		return;
		
	//for PID
	fPID = new AliESDpid();
	if (!fPID) {
		printf("ERROR: pID not available\n");
		return;
	}
	fPID->MakePID(fESD);
	
	// for MC
	if (fIsMc)
	{
		AliMCEvent* mcEvent = MCEvent();
		if (!mcEvent) {
			printf("ERROR: MC Event not available\n");
			return;
		}
	
		fStack = mcEvent->Stack();
	}
	
	// calorimeter
	// define geometry utils
	fGeoUt = new AliEMCALGeometry("EMCAL_COMPLETE1","EMCAL");
	fGeoUt->SetMisalMatrix(fESD->GetEMCALMatrix(0),0);

	// retrieve clusters
	fCaloClusters = new TRefArray();
	fESD->GetEMCALClusters( fCaloClusters );	
	fCells = fESD->GetEMCALCells();	
	if (!fCells) {
		printf("ERROR: cells not available\n");
		return;
	}
	
	//Int_t res = CheckPhysicsSelection(fESD->GetRunNumber());

	//AliCentrality *cent = fESD->GetCentrality();

	//if (res == 0 && cent)
	//{
		//if (IsPhysicsSelected())
		//{
			fHistVertexZ->Fill(fVtxStat,zVertex);
			
			Run();
		//}
	//}
	
	delete fPID;
	delete fCaloClusters;
	delete fGeoUt;
	
	PostData(1, fOutputList);	
}

//________________________________________________________________________
void AliAnalysisTaskCaloElectron::Terminate(Option_t *) 
{
	// Draw result to the screen
	// Called once at the end of the query
	
	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList) {
		printf("ERROR: Output list not available\n");
		return;
	}
}

//________________________________________________________________________
void AliAnalysisTaskCaloElectron::Run() 
{
	// for MC
	Int_t nStackTracks = 0;
	Int_t iPart = 0;
	Int_t pdgMomCode = 0;
	AliESDCaloCluster* caloClusterMCMatch = 0;
	TVector3 caloPosMCMatch(0,0,0);
	Double_t resMCMatch = 0;
	
	// cluster variables
	AliESDCaloCluster* caloCluster = 0;
	TVector3 caloPos(0,0,0);
	Float_t pos[3] = {0};
	Double_t raddeg = TMath::RadToDeg();
	
	// tracking
	AliESDtrack* track = 0;
	AliEMCALTrack *emcTrack = 0;	
	TVector3 trackPos(0,0,0);
	Double_t res = 0;
	Float_t dEdx = 0;
	Double_t invMassTPC = 0, invMassEMCAL = 0;
	
	// loop over the stack to fill MC in TTree
	if (fIsMc)
	{
		if (fStack)
			nStackTracks = fStack->GetNtrack();
		
		for (Int_t iiPart = 0; iiPart < nStackTracks; iiPart++)
		{
			TParticle *part = fStack->Particle(iiPart);
			
			if (!part)
			{
				Printf("ERROR: Could not get particle %d", iiPart);
				continue;
			}
			
			TParticlePDG *pdg = part->GetPDG(0);
			
			if (!pdg)
			{
				Printf("ERROR: Could not get particle PDG %d", iiPart);
				continue;
			}		
			
			pdgMomCode = GetMCMother(part);	
			
			if ( ((pdg->PdgCode() == -11) || (pdg->PdgCode() == 11)) && (TMath::Abs(part->Eta())<1.) )
			{
				if (fStack->IsPhysicalPrimary(iiPart))
					fHistEtaPtStackeMC->Fill(part->Pt(),part->Eta());
				
				if (pdgMomCode == 22)
					fHistEtaPtStackePhMC->Fill(part->Pt(),part->Eta());
				
				if (part->Pt()>1.)
				{
					fEtaMCIn = part->Eta();
					fPtMCIn = part->Pt();
					fPMCIn = part->P();
					fEMCIn = part->Energy();
					fvXMCIn = part->Vx();
					fvYMCIn = part->Vy();
					fvZMCIn = part->Vz();
					fIsPhyPrimIn = fStack->IsPhysicalPrimary(iiPart);
					fMomIn = pdgMomCode;
					
					fTreeMCIn->Fill();
				}
				//fHistVtxXYeMC->Fill(part->Vx(),part->Vy());			
			} // end of if electron
		} // end of loop over stack		
	} // end of fill MC in TTree
		
	// loop over clusters
	Int_t nCluster = fCaloClusters->GetEntries();
	for (Int_t iCluster = 0; iCluster < nCluster; iCluster++ ) 
	{
		caloCluster = ( AliESDCaloCluster* )fCaloClusters->At( iCluster );
		
		if (!caloCluster) continue;
		
		if(!caloCluster->IsEMCAL()) continue;
		if(!fRecoUtil->CheckCellFiducialRegion(fGeom, caloCluster, fCells)) continue;
		if(fRecoUtil->ClusterContainsBadChannel(fGeom, caloCluster->GetCellsAbsId(), caloCluster->GetNCells() )) continue; 	
		
		caloCluster->GetPosition(pos);		
		caloPos.SetXYZ(pos[0],pos[1],pos[2]);
	
		fHistCluE->Fill(caloCluster->E());
		fHistCluEtaPhi->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
		
		if (fIsMc)
		{		
			iPart = (Int_t)TMath::Abs(caloCluster->GetLabel());
			
			if (IsElectronMC(iPart))
			{
				fHistCluEeMC->Fill(caloCluster->E());
				fHistCluEtaPhieMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
			}
				
			if (IsElectronMC(iPart, "ePhMC"))
			{
				fHistCluEePhMC->Fill(caloCluster->E());
				fHistCluEtaPhiePhMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
			}		
		}
	}			
	
	// do matching
	//AliEMCALRecoUtils::FindMatches(fESD, fCaloClusters, dataType)
	
	// loop over tracks
	TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
	//TObjArray* list = fEsdTrackPhCuts->GetAcceptedTracks(fESD);
	Int_t nGoodTracks = list->GetEntries();
	
	for (Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) 
	{
		track = dynamic_cast<AliESDtrack*> (list->At(iTracks));
		emcTrack = new AliEMCALTrack(*track);
		
		fHistEtaPt->Fill(track->Pt(),track->Eta());
		
		dEdx = track->GetTPCsignal();
		fHistdEdx->Fill(track->P(),dEdx);
		
		if (fIsMc)
		{
			iPart = (Int_t)TMath::Abs(track->GetLabel());
			pdgMomCode = GetMCMother(iPart);	
			
			if (IsElectronMC(iPart))
			{
				fHistEtaPteMC->Fill(track->Pt(),track->Eta());
				fHistdEdxeMC->Fill(track->P(),dEdx);
			}
			if (IsElectronMC(iPart, "ePhMC"))
			{
				fHistEtaPtePhMC->Fill(track->Pt(),track->Eta());
				fHistdEdxePhMC->Fill(track->P(),dEdx);
			}			
		}
				
		if (IsElectronTPC(track))
		{
			fHistEtaPteTPC->Fill(track->Pt(),track->Eta());
			fHistdEdxeTPC->Fill(track->P(),dEdx);

			if (fIsMc)
			{
				if (IsElectronMC(iPart))
				{
					fHistEtaPteTPCeMC ->Fill(track->Pt(),track->Eta());
					fHistdEdxeTPCeMC->Fill(track->P(),dEdx);
				}
				if (IsElectronMC(iPart, "ePhMC"))
				{
					fHistEtaPteTPCePhMC->Fill(track->Pt(),track->Eta());
					fHistdEdxeTPCePhMC->Fill(track->P(),dEdx);
				}			
			}			
			
			invMassTPC = RunPhotonic(track, "TPC");

			caloCluster = FindMatch(emcTrack,trackPos);

			if (caloCluster)
			{
				fHistEtaPteTPCMatch->Fill(track->Pt(),track->Eta());
				fHistdEdxeTPCMatch->Fill(track->P(),dEdx);

				caloCluster->GetPosition(pos);		
				caloPos.SetXYZ(pos[0],pos[1],pos[2]);
				
				fHistCluEeTPCMatch->Fill(caloCluster->E());
				fHistCluEtaPhieTPCMatch->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
				
				res = sqrt(pow(trackPos.Phi()-caloPos.Phi(),2)+pow(trackPos.Eta()-caloPos.Eta(),2));
				
				fHistDeltaRZeTPCMatch->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
				fHistResPteTPCMatch->Fill(track->Pt(),res);
				fHistEoverPPteTPCMatch->Fill(track->Pt(),caloCluster->E()/track->P());
				fHistEoverPPeTPCMatch->Fill(track->P(),caloCluster->E()/track->P());
				fHistEoverPEeTPCMatch->Fill(caloCluster->E(),caloCluster->E()/track->P());

				if (fIsMc)
				{
					if (IsElectronMC(iPart))
					{
						fHistEtaPteTPCMatcheMC->Fill(track->Pt(),track->Eta());
						fHistdEdxeTPCMatcheMC->Fill(track->P(),dEdx);

						fHistCluEeTPCMatcheMC->Fill(caloCluster->E());
						fHistCluEtaPhieTPCMatcheMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());

						fHistDeltaRZeTPCMatcheMC->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
						fHistResPteTPCMatcheMC->Fill(track->Pt(),res);
						fHistEoverPPteTPCMatcheMC->Fill(track->Pt(),caloCluster->E()/track->P());
					}
					if (IsElectronMC(iPart, "ePhMC"))
					{
						fHistEtaPteTPCMatchePhMC->Fill(track->Pt(),track->Eta());
						fHistdEdxeTPCMatchePhMC->Fill(track->P(),dEdx);
						
						fHistCluEeTPCMatchePhMC->Fill(caloCluster->E());
						fHistCluEtaPhieTPCMatchePhMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());

						fHistDeltaRZeTPCMatchePhMC->Fill(trackPos.Phi()-caloPos.Phi(),trackPos.Eta()-caloPos.Eta());
						fHistResPteTPCMatchePhMC->Fill(track->Pt(),res);
						fHistEoverPPteTPCMatchePhMC->Fill(track->Pt(),caloCluster->E()/track->P());
					}
				}
				
				// check for matching
				//if ((Res < fCutRMatch) && IsElectronEMCAL(caloCluster,track))
				if (res < fCutRMatch)
				{
					fHistEoverPPteEMCAL->Fill(track->Pt(),caloCluster->E()/track->P());
					fHistEoverPPeEMCAL->Fill(track->P(),caloCluster->E()/track->P());
					fHistEoverPEeEMCAL->Fill(caloCluster->E(),caloCluster->E()/track->P());

					//if (caloCluster->E()>3.5)
					//	fHistEoverPPeEMCALtrigger->Fill(track->P(),caloCluster->E()/track->P());

					if (fIsMc)
					{
						if (IsElectronMC(iPart))
						{
							fHistEoverPPteEMCALeMC->Fill(track->Pt(),caloCluster->E()/track->P());
						}
						if (IsElectronMC(iPart, "ePhMC"))
						{
							fHistEoverPPteEMCALePhMC->Fill(track->Pt(),caloCluster->E()/track->P());	
						}
					}
					
					if (IsElectronEMCAL(caloCluster,track))
					{
						fHistEtaPteEMCAL->Fill(track->Pt(),track->Eta());
						fHistdEdxeEMCAL->Fill(track->P(),dEdx);
						
						fHistCluEeEMCAL->Fill(caloCluster->E());
						fHistCluEtaPhieEMCAL->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
						
						invMassEMCAL = RunPhotonic(track, "EMCAL", caloCluster);
						
						if (fIsMc)
						{
							if (IsElectronMC(iPart))
							{
								fHistEtaPteEMCALeMC->Fill(track->Pt(),track->Eta());
								fHistdEdxeEMCALeMC->Fill(track->P(),dEdx);
								
								fHistCluEeEMCALeMC->Fill(caloCluster->E());
								fHistCluEtaPhieEMCALeMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
							}
							if (IsElectronMC(iPart, "ePhMC"))
							{
								fHistEtaPteEMCALePhMC->Fill(track->Pt(),track->Eta());
								fHistdEdxeEMCALePhMC->Fill(track->P(),dEdx);
								
								fHistCluEeEMCALePhMC->Fill(caloCluster->E());
								fHistCluEtaPhieEMCALePhMC->Fill(raddeg*caloPos.Phi(),caloPos.Eta());
							}
						}
					}
				}
			}
		} // end of is electron TPC
		
		//Fill TTree variables
		if ((dEdx>60) && (dEdx<110) && (track->Pt()>1.))
		{
			fEta = track->Eta();
			fPt = track->Pt(); 
			fP=track->P(); 
			fDedx = dEdx; 
			fNSigma = fPID->NumberOfSigmasTPC(track,AliPID::kElectron);
			track->GetImpactParameters(fDcaXY,fDcaZ);
			
			if (caloCluster)
			{
				fE=caloCluster->E(); 
				fPhiClu = raddeg*caloPos.Phi();
				fEtaClu = caloPos.Eta();
				fRes=res;
			}
			else
			{
				fE = -1;
				fPhiClu = -1;
				fEtaClu = -1;
				fRes = -1;
			}
					
			fInvMass = RunPhotonic(track, "NONE");
			
			fTree->Fill();
		}
		
		// find cluster that matches the same (MC) track
		if (fIsMc)
		{
			caloClusterMCMatch = FindMCMatch(iPart);
			
			if (caloClusterMCMatch)
			{
				if (IsElectronMC(iPart))
				{
					fHistEtaPteMCMatch->Fill(track->Pt(),track->Eta());
					fHistdEdxeMCMatch->Fill(track->P(),dEdx);
				}
				if (IsElectronMC(iPart, "ePhMC"))
				{
					fHistEtaPtePhMCMatch->Fill(track->Pt(),track->Eta());
					fHistdEdxePhMCMatch->Fill(track->P(),dEdx);
				}
				
				caloClusterMCMatch->GetPosition(pos);		
				caloPosMCMatch.SetXYZ(pos[0],pos[1],pos[2]);
				
				if (GetTrackProjection(emcTrack,trackPos,caloPosMCMatch))
				{
					resMCMatch = sqrt(pow(trackPos.Phi()-caloPosMCMatch.Phi(),2)+pow(trackPos.Eta()-caloPosMCMatch.Eta(),2));
					
					if (IsElectronMC(iPart))
					{
						fHistDeltaRZeMCMatch->Fill(trackPos.Phi()-caloPosMCMatch.Phi(),trackPos.Eta()-caloPosMCMatch.Eta());
						fHistResPteMCMatch->Fill(track->Pt(),resMCMatch);	
						fHistEoverPPteMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
					}
					if (IsElectronMC(iPart, "ePhMC"))
					{
						fHistDeltaRZePhMCMatch->Fill(trackPos.Phi()-caloPosMCMatch.Phi(),trackPos.Eta()-caloPosMCMatch.Eta());
						fHistResPtePhMCMatch->Fill(track->Pt(),resMCMatch);	
						fHistEoverPPtePhMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
					}
					
					//if ((ResMCMatch < fCutRMatch) && (IsElectronEMCAL(caloClusterMCMatch,track)))
					if (resMCMatch < fCutRMatch)
					{
						if (IsElectronMC(iPart))
						{
							fHistEoverPPteEMCALeMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
						}
						if (IsElectronMC(iPart, "ePhMC"))
						{
							fHistEoverPPteEMCALePhMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
						}					
					}
					
					if (IsElectronTPC(track))
					{
						if (IsElectronMC(iPart))
						{
							fHistEtaPteTPCeMCMatch->Fill(track->Pt(),track->Eta());
							fHistdEdxeTPCeMCMatch->Fill(track->P(),dEdx);
							
							fHistDeltaRZeTPCeMCMatch->Fill(trackPos.Phi()-caloPosMCMatch.Phi(),trackPos.Eta()-caloPosMCMatch.Eta());
							fHistResPteTPCeMCMatch->Fill(track->Pt(),resMCMatch);	
							fHistEoverPPteTPCeMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
						}
						if (IsElectronMC(iPart, "ePhMC"))
						{
							fHistEtaPteTPCePhMCMatch->Fill(track->Pt(),track->Eta());
							fHistdEdxeTPCePhMCMatch->Fill(track->P(),dEdx);
							
							fHistDeltaRZeTPCePhMCMatch->Fill(trackPos.Phi()-caloPosMCMatch.Phi(),trackPos.Eta()-caloPosMCMatch.Eta());
							fHistResPteTPCePhMCMatch->Fill(track->Pt(),resMCMatch);	
							fHistEoverPPteTPCePhMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());
						}
						
						//if ((ResMCMatch < fCutRMatch) && (IsElectronEMCAL(caloClusterMCMatch,track)))
						if (resMCMatch < fCutRMatch)
						{
							if (IsElectronMC(iPart))
							{
								fHistEoverPPteTPCeEMCALeMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());					
							}
							if (IsElectronMC(iPart, "ePhMC"))
							{
								fHistEoverPPteTPCeEMCALePhMCMatch->Fill(track->Pt(),caloClusterMCMatch->E()/track->P());					
							}
						}
					}
				}
			} // end of caloClusterMCMatch			
			
			//Fill MC TTree variables
			if ( (IsElectronMC(iPart)) && (track->Pt()>1.) )
			{
				fEtaMC = track->Eta();
				fPtMC = track->Pt(); 
				fPMC = track->P(); 
				fDedxMC = dEdx; 
				fNSigmaMC = fPID->NumberOfSigmasTPC(track,AliPID::kElectron);
				track->GetImpactParameters(fDcaXYMC,fDcaZMC);
				
				if (caloCluster)
				{
					fEMC = caloCluster->E(); 
					fPhiCluMC = raddeg*caloPos.Phi();
					fEtaCluMC = caloPos.Eta();
					fResMC = res;
				}
				else
				{
					fEMC = -1;
					fPhiCluMC = -1;
					fEtaCluMC = -1;					
					fResMC = -1;
				}
				
				if (caloClusterMCMatch)
				{
					fEMCMatch = caloClusterMCMatch->E(); 
					fPhiCluMCMatch = raddeg*caloPosMCMatch.Phi();
					fEtaCluMCMatch = caloPosMCMatch.Eta();
					fResMCMatch = resMCMatch;
				}
				else
				{
					fEMCMatch = -1;
					fPhiCluMCMatch = -1;
					fEtaCluMCMatch = -1;
					fResMCMatch = -1;
				}
								
				fInvMassMC = RunPhotonic(track, "ALL");
				
				fIsPhyPrim = fStack->IsPhysicalPrimary(iPart);
				fMom = pdgMomCode;
				
				fTreeMC->Fill();
			} // end of Fill MC TTree variables
		} // end of IsMc
		
	} //end of track loop 
}      

//________________________________________________________________________
Double_t AliAnalysisTaskCaloElectron::RunPhotonic(AliESDtrack * trackIN, TString type, AliESDCaloCluster *caloCluster) 
{//Marcelo add comment
	AliESDtrack* track = 0;
	
	Double_t mass = 0, massMin = 999, sMass = 0;

	// Define first KF particle
	Int_t pid = 0;
	if (trackIN->Charge() <0)
		pid = 11;
	else if (trackIN->Charge() > 0)
		pid = -11;

	AliKFParticle trackIN_KF(*trackIN, pid);	
	
	// MC part
	Int_t iPartIN  = (Int_t)TMath::Abs(trackIN->GetLabel());
	
	// do matching
	//AliEMCALRecoUtils::FindMatches(fESD, fCaloClusters, dataType)
	
	// loop over tracks
	TObjArray* list = fEsdTrackPhCuts->GetAcceptedTracks(fESD);
	Int_t nGoodTracks = list->GetEntries();
	
	for (Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) 
	{
		track = fESD->GetTrack(iTracks);
		
		if (IsElectronTPC(track))
		{
			// define second KF particle and photon decay
			if (track->Charge() <0)
				pid = 11;
			else if (track->Charge() > 0)
				pid = -11;
			
			AliKFParticle track_KF(*track, pid);	

			AliKFParticle photon(trackIN_KF, track_KF);
			photon.SetProductionVertex(AliKFVertex(*fVertex));
			//photon.SetVtxGuess(0,0,0);
			photon.SetMassConstraint(0.0);		
			
			//mass = CalcMass(trackIN, track);
			if (photon.GetMass(mass,sMass))
				continue;
			
			if ( (track->Charge() == trackIN->Charge()) && (track->GetLabel() != trackIN->GetLabel()) )
			{
				if (type == "TPC")
				{
					fHistMassLikePteTPC->Fill(trackIN->Pt(),mass);

					if (fIsMc)
					{
						if (IsElectronMC(iPartIN,"ePhMC"))
						{
							fHistMassLikePteTPCePhMC->Fill(trackIN->Pt(),mass);
						}
						else if (IsElectronMC(iPartIN))
						{
							fHistMassLikePteTPCeNonPhMC->Fill(trackIN->Pt(),mass);							
						}
					}
				} // end of TPC
				else if (type == "EMCAL")
				{
					fHistMassLikePteEMCAL->Fill(trackIN->Pt(),mass);					
					fHistMassLikePeEMCAL->Fill(trackIN->P(),mass);					
					fHistMassLikeEeEMCAL->Fill(caloCluster->E(),mass);					
				
					if (fIsMc)
					{
						if (IsElectronMC(iPartIN,"ePhMC"))
						{
							fHistMassLikePteEMCALePhMC->Fill(trackIN->Pt(),mass);					
							fHistMassLikePeEMCALePhMC->Fill(trackIN->P(),mass);					
							fHistMassLikeEeEMCALePhMC->Fill(caloCluster->E(),mass);											
						}
						else if (IsElectronMC(iPartIN))
						{
							fHistMassLikePteEMCALeNonPhMC->Fill(trackIN->Pt(),mass);					
							fHistMassLikePeEMCALeNonPhMC->Fill(trackIN->P(),mass);					
							fHistMassLikeEeEMCALeNonPhMC->Fill(caloCluster->E(),mass);											
						}
					}
				} // end of EMCAL
				
			} // end of like charge
			else if ( (track->Charge() != trackIN->Charge()) )
			{
				if (mass < massMin)
					massMin = mass;
				
				if (type == "TPC")
				{
					fHistMassUnLikePteTPC->Fill(trackIN->Pt(),mass);

					if (fIsMc)
					{
						if (IsElectronMC(iPartIN,"ePhMC"))
						{
							fHistMassUnLikePteTPCePhMC->Fill(trackIN->Pt(),mass);
						}
						else if (IsElectronMC(iPartIN))
						{
							fHistMassUnLikePteTPCeNonPhMC->Fill(trackIN->Pt(),mass);							
						}
					}
				} // end of TPC
				else if (type == "EMCAL")
				{
					fHistMassUnLikePteEMCAL->Fill(trackIN->Pt(),mass);					
					fHistMassUnLikePeEMCAL->Fill(trackIN->P(),mass);					
					fHistMassUnLikeEeEMCAL->Fill(caloCluster->E(),mass);					

					if (fIsMc)
					{
						if (IsElectronMC(iPartIN,"ePhMC"))
						{
							fHistMassUnLikePteEMCALePhMC->Fill(trackIN->Pt(),mass);					
							fHistMassUnLikePeEMCALePhMC->Fill(trackIN->P(),mass);					
							fHistMassUnLikeEeEMCALePhMC->Fill(caloCluster->E(),mass);											
						}
						else if (IsElectronMC(iPartIN))
						{
							fHistMassUnLikePteEMCALeNonPhMC->Fill(trackIN->Pt(),mass);					
							fHistMassUnLikePeEMCALeNonPhMC->Fill(trackIN->P(),mass);					
							fHistMassUnLikeEeEMCALeNonPhMC->Fill(caloCluster->E(),mass);											
						}
					}
				} // end of EMCAL			
			} // end of un-like charge
		} // end of - is Electron TPC	
	}
	
	return massMin;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsGoodEvent()
{//Is this a good event?
	if (fESD)
	{
		//if (fVertex->GetStatus())
		if ((fESD->GetNumberOfTracks() > 1) && (fVtxStat>0))
		{
			if (TMath::Abs(fVertex->GetZ()) < 10)		
				return kTRUE;
		}
	}
	
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsGoodTrack(AliESDtrack *track)
{//Is this a good track?
	Float_t dcaXY=0, dcaZ=0;
	Float_t chi2=0;
	
	if (track)
	{
		track->GetImpactParameters(dcaXY,dcaZ);
		if (track->GetTPCNcls() > 0)
			chi2 = track->GetTPCchi2()/track->GetTPCNcls();
		else
			chi2 = 999;
		
		if ((track->GetTPCNcls() > 80) && (track->GetTPCsignalN() > 80) && (dcaXY < 0.05) && (dcaZ < 0.05) && (chi2 < 3.5)) 
			return kTRUE;
	}
	
	return kFALSE;
}

//____________________________________________________________________
AliESDtrackCuts* AliAnalysisTaskCaloElectron::GetTrackCuts() const
{//Get the track cuts
	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
	
	// TPC  
	esdTrackCuts->SetMinNClustersTPC(80);
	esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	// ITS
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
										   AliESDtrackCuts::kAny);
	esdTrackCuts->SetMaxDCAToVertexXY(0.05);
	esdTrackCuts->SetMaxDCAToVertexZ(0.05);
	//esdTrackCuts->SetDCAToVertex2D(kFALSE);
	esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
	esdTrackCuts->SetEtaRange(-0.9,0.9);
	
	return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliAnalysisTaskCaloElectron::GetTrackPhCuts() const
{//Get the track Ph Cuts
	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
	
	// TPC  
	esdTrackCuts->SetMinNClustersTPC(50);
	esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	// ITS
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	//esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	//									   AliESDtrackCuts::kAny);
	esdTrackCuts->SetEtaRange(-0.9,0.9);
	
	return esdTrackCuts;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsElectronTPC(AliESDtrack *track)
{//Is the track an electron using the TPC dEdx
	Float_t nSigmaElectron;
	Float_t dEdx = track->GetTPCsignal();
	
	nSigmaElectron = TMath::Abs(fPID->NumberOfSigmasTPC(track,AliPID::kElectron));
	
	//if (nSigmaElectron<3.0)
	if ((dEdx>fCutdEdxMin) && (dEdx<fCutdEdxMax))
		return kTRUE;
	else
		return kFALSE;	
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsElectronEMCAL(AliESDCaloCluster *caloCluster, AliESDtrack *track)
{//Is the electron an electron from the EMCAL
	Double_t eOverP = 0;
	Double_t caloE = caloCluster->E();
	Double_t p = track->P();
	
	if (p>0)
	{
		eOverP = caloE/p;
		
		if ((eOverP > fCutEoverPMin) && (eOverP < fCutEoverPMax))
			return kTRUE;
	}

	return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsElectronMC(Int_t iPart, TString type)
{	//Is the electron a MC electron
	if (fStack)
	{
		Int_t nStackTracks = fStack->GetNtrack();
		
		if ((iPart>=0) && (iPart < nStackTracks))
		{
			TParticle *part = fStack->Particle(iPart);
			TParticlePDG *pdg = 0;
			
			if (part)
			{
				pdg = part->GetPDG(0);
				
				if (pdg)
				{
					if ((pdg->PdgCode() == -11) || (pdg->PdgCode() == 11))
						//if ((pdg->PdgCode() == -211) || (pdg->PdgCode() == 211))
					{
						if (type == "ePhMC")
						{
							if (GetMCMother(part) == 22)
							{
								return kTRUE;
							}
						}
						else
							return kTRUE;
					}
				}
			}
		}
	}
	return kFALSE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskCaloElectron::GetMCMother(TParticle *part)
{//What is the MC mother track
	TParticlePDG *pdg;
	Int_t pdgMomCode = -999;	
	Int_t nStackTracks = fStack->GetNtrack();
	
	if (part)
	{
		pdg = part->GetPDG(0);
		
		if (pdg)
		{
			TParticle *partMom = 0;
			
			Int_t iPartMom = part->GetMother(0);
			
			TParticlePDG *pdgMom = 0;
			
			if ((iPartMom>=0) && (iPartMom < nStackTracks))
			{
				partMom = fStack->Particle(iPartMom);
				if (partMom)
					pdgMom = partMom->GetPDG(0);
			}
			
			if (pdgMom)
				pdgMomCode = pdgMom->PdgCode();
		}
	}
	
	return pdgMomCode;
}

//________________________________________________________________________
Int_t AliAnalysisTaskCaloElectron::GetMCMother(Int_t iPart)
{//Get the MC mother track
	Int_t pdgMomCode = -999;	

	if (fStack)
	{
		Int_t nStackTracks = fStack->GetNtrack();
		
		if ((iPart>=0) && (iPart < nStackTracks))
		{
			TParticle *part = fStack->Particle(iPart);
			TParticlePDG *pdg = 0;
			
			if (part)
			{
				pdg = part->GetPDG(0);
				
				if (pdg)
				{
					TParticle *partMom = 0;
					
					Int_t iPartMom = part->GetMother(0);
					
					TParticlePDG *pdgMom = 0;
					
					if ((iPartMom>=0) && (iPartMom < nStackTracks))
					{
						partMom = fStack->Particle(iPartMom);
						if (partMom)
							pdgMom = partMom->GetPDG(0);
					}
					
					if (pdgMom)
						pdgMomCode = pdgMom->PdgCode();
				}
			}
		}
	}
	return pdgMomCode;
}

//________________________________________________________________________
void AliAnalysisTaskCaloElectron::InitCaloUtil()
{//Initiate calo utils
	fRecoUtil = new AliEMCALRecoUtils();

	//	fRecoUtil->InitEMCALGeometry();
	fRecoUtil->SetNumberOfCellsFromEMCALBorder(1);
	//    cu->SetNumberOfCellsFromPHOSBorder(2);
    fRecoUtil->SwitchOnNoFiducialBorderInEMCALEta0();
    fRecoUtil->SwitchOnBadChannelsRemoval(); 
	
	const Int_t nTowers=80;
	Int_t hotChannels[nTowers]={103, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 368, 369, 370, 371, 372, 373, 374,375, 376, 377, 378, 379, 380, 381, 382, 383, 1275, 1288, 1519, 1860, 1967, 2022, 2026, 2047, 2117, 2298, 2776, 6095, 6111, 6592, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7371, 7425, 7430, 7457, 7491, 7709, 8352, 8353, 8356, 8357, 8808, 8810, 8812, 8814, 9056, 9815, 9837};
	
	Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
	for(Int_t i=0; i<nTowers; i++)
    {
		fGeom->GetCellIndex(hotChannels[i],nSupMod,nModule,nIphi,nIeta);
		fGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
		fRecoUtil->SetEMCALChannelStatus(nSupMod, ieta, iphi);
    }
	
	Int_t myChannels[9] = {74,152,917,1595,2540,3135,3764,6481,9769};
	for(Int_t i=0; i<9; i++)
    {
		fGeom->GetCellIndex(myChannels[i],nSupMod,nModule,nIphi,nIeta);
		fGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
		fRecoUtil->SetEMCALChannelStatus(nSupMod, ieta, iphi);
    }
	
	fRecoUtil->Print("");
}

//________________________________________________________________________
AliESDCaloCluster* AliAnalysisTaskCaloElectron::FindMatch(AliEMCALTrack *track, TVector3 &trackPos)
{	//find matched tracks
	AliESDCaloCluster* caloCluster = 0;
	AliESDCaloCluster* caloClusterMatch = 0;

	Int_t nCluster = fCaloClusters->GetEntries();

	Double_t res=0, resMin=999;
	
	Float_t pos[3] = {0};
	
	TVector3 caloPos(0,0,0);

	// loop over clusters
	for (Int_t iCluster = 0; iCluster < nCluster; iCluster++ ) 
	{
		caloCluster = ( AliESDCaloCluster* )fCaloClusters->At( iCluster );
		
		if (!caloCluster) continue;
		
		if(!caloCluster->IsEMCAL()) continue;
		if(!fRecoUtil->CheckCellFiducialRegion(fGeom, caloCluster, fCells)) continue;
		if(fRecoUtil->ClusterContainsBadChannel(fGeom, caloCluster->GetCellsAbsId(), caloCluster->GetNCells() )) continue; 	
		
		caloCluster->GetPosition(pos);		
		caloPos.SetXYZ(pos[0],pos[1],pos[2]);
		
		if (GetTrackProjection(track,trackPos,caloPos))
		{
			res = sqrt(pow(trackPos.Phi()-caloPos.Phi(),2)+pow(trackPos.Eta()-caloPos.Eta(),2));

			if (res < resMin)
			{
				resMin = res;
				caloClusterMatch = caloCluster;
			}
		}
	}
	
	return caloClusterMatch;
}

//________________________________________________________________________
AliESDCaloCluster* AliAnalysisTaskCaloElectron::FindMCMatch(Int_t iPart)
{	//find MC matched tracks
	AliESDCaloCluster* caloCluster = 0;
	Int_t nCluster = fCaloClusters->GetEntries();
	Int_t iPartClu = 0;

	// loop over clusters
	for (Int_t iCluster = 0; iCluster < nCluster; iCluster++ ) 
	{
		caloCluster = ( AliESDCaloCluster* )fCaloClusters->At( iCluster );
		
		if (!caloCluster) continue;

		iPartClu = (Int_t)TMath::Abs(caloCluster->GetLabel());
		
		if (iPartClu == iPart) // found match
			return caloCluster;
	}
	
	return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos)
{//Get the track projection
	Bool_t proj = kFALSE;
	
	if (emcTrack)
	{	
		Double_t trkPos[3] = {0};
		
		emcTrack->PropagateToGlobal(clusPos.X(),clusPos.Y(),clusPos.Z(),0.,0.);
		emcTrack->GetXYZ(trkPos);
		
		trackPos.SetXYZ(trkPos[0],trkPos[1],trkPos[2]);
		
		proj = kTRUE;
	}
	
	return proj;
}

//________________________________________________________________________
Int_t AliAnalysisTaskCaloElectron::CheckPhysicsSelection(Int_t runNumber)
{//Check the physics selection
	// Check if the physics selection is valid, if not load a new one
    if (runNumber == fCurrentRunNum || fIsMc)
    {
        return 0;
    }
    else
    {
        AliPhysicsSelection *selection = 0;
		AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
		if(!mgr){
			AliError("Error: no analysis manager");
			return -1;
		}
        AliPhysicsSelectionTask *physSelTask = dynamic_cast<AliPhysicsSelectionTask*>(mgr->GetTask("AliPhysicsSelectionTask"));
        if (physSelTask)
        {	
			// define selection
			selection = new AliPhysicsSelection();
			selection->AddCollisionTriggerClass("+CINT1-B-NOPF-ALLNOTRD");
			//selection->AddCollisionTriggerClass("+CEMC1-B-NOPF-ALLNOTRD");
        }
        else
        {
            AliError("Could not get physics selection task from manager, undefined or no selection will be used");
            return -1;
        }
        AliInfo("Changing the physics selection");
        AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (mgr->GetInputEventHandler());
        // The physics selection task has a bit weird implementation, setting the the physics selection will not update the handler, so we do it manually
		if(!handler){
            AliError("Analysis manager does not exist!");
            return -1;
		}
        physSelTask->SetPhysicsSelection(selection);
        handler->SetEventSelection(selection);
        fCurrentRunNum = runNumber;		
    }
	
    return 1;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCaloElectron::IsPhysicsSelected() 
{	//Is this event selected?
	// See header file for class documentation
    if(fIsMc) return kTRUE;
	
    //return ( !(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kFastOnly) && 
	//		  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB) );
    return (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);

    //TString trgClasses = fESD->GetFiredTriggerClasses();
    //if ( (!fSelectTrigger) || (fSelectTrigger && trgClasses.Contains(fTrigger.Data())) )
	//	return kTRUE;

	//return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskCaloElectron::SetTriggerSelection(Bool_t triggerFlag, const char *trigger)
{//set the trigger selection
	fSelectTrigger = triggerFlag;
	fTrigger = trigger;
}

//________________________________________________________________________
Double_t AliAnalysisTaskCaloElectron::CalcMass(AliESDtrack *trackIN, AliESDtrack *track)
{//calculate the track mass
	Double_t e1 = TMath::Sqrt( TMath::Power(trackIN->P(),2.0) + TMath::Power(fgMassElectron,2.0) );
	Double_t e2 = TMath::Sqrt( TMath::Power(track->P(),2.0) + TMath::Power(fgMassElectron,2.0) );
	Double_t mass = TMath::Sqrt( TMath::Power(e1 + e2,2.0) - ( TMath::Power(trackIN->Px()+track->Px(),2.0) + TMath::Power(trackIN->Py()+track->Py(),2.0) + TMath::Power(trackIN->Pz()+track->Pz(),2.0) ) );
	
	return mass;
}
