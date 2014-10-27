/* Leading Charged Track+V0 Correlation.(Works for Real,Monte Carlo Data)
 *                            Sandun Jayarathna
 *                          University of Houston
 *                      sandun.pahula.hewage@cern.ch
 *****************************************************************************************/
#include <TROOT.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TRandom.h>
#include <THnSparse.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliVParticle.h"
#include "AliMultiplicity.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliExternalTrackParam.h"
#include "AliAnalyseLeadingTrackUE.h"

#include "AliLeadingV0Correlation.h"

#define CorrBinsX 240
#define CorrBinsY 260


Double_t PI =TMath::Pi();

ClassImp(AliLeadingV0Correlation)
ClassImp(V0Correlationparticle)

//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation()
   : AliAnalysisTaskSE(),
	fAODEvent					(0x0),
	fPoolMgr					(0x0),
	fPIDResponse				(0x0),
	fAnalyseUE					(0x0),
	fPoolMaxNEvents				(0), 
	fPoolMinNTracks				(0), 
	fMinEventsToMix				(0),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fcollidingSys				(""),
	fpvzcut						(0),
	fTrackEtaCut				(0.9),
	fFilterBit					(128),
	fAnalysisMC					(0),
	fCase						(0),
	fRemoveAutoCorr				(0),
	fRapidityCut				(0),
	fV0radius					(0),
	fV0PostoPVz					(0),
	fV0NegtoPVz					(0),
	fDCAV0Daughters				(0),
	fCPAK0						(0),
	fCPALam						(0),
	fRejectLamK0				(0),
	fRejectK0Lam				(0),
    fSigmaPID					(0),
	fCutCTK0					(0),
	fCutCTLa					(0),
	fMassCutK0					(0),
	fMassCutLa					(0),
	fTriglow					(0),
	fTrighigh					(0),
	fTPCClusters				(0),					
	fTPCfindratio				(0),
	fUseChargeHadrons			(kTRUE), 
	fPtMin						(0.15),
	fOutputList					(0),
	fHist_Mult_B4_Trg_Sel		(0),
	fHist_Mult_Af_Trg_Sel		(0),
	fHist_Mult_PVz_Cut			(0),
	fHist_Mult_SPD_PVz			(0),
	fHist_Mult_SPD_PVz_Pileup	(0),
	fHistPVx					(0),
	fHistPVy					(0),
	fHistPVz					(0),
	fHistPVxAnalysis			(0),
	fHistPVyAnalysis			(0),
	fHistPVzAnalysis			(0),
	fHistEventViceGen			(0),
	fHistEventViceReconst		(0),
	fHistMCGenK0				(0),
	fHistMCGenLAM				(0),
	fHistMCGenALAM				(0),
	fHistMCGenLAMXIPLS			(0),
	fHistMCGenLAMXI 			(0),
	fHistReconstK0				(0),
	fHistReconstLA				(0),
	fHistReconstALA				(0),
	fHistMCAssoK0				(0),
	fHistMCAssoLA				(0),
	fHistMCAssoALA				(0),
	fHistMCAssoLAXI				(0),
	fHistMCAssoALAXiPlus		(0),
	fHistReconstSib				(0),
	fHistReconstMix				(0),
	fHistReconstSibGEN			(0),
	fHistReconstMixGEN			(0),
	fHistReconstSibASO			(0),
	fHistReconstMixASO			(0),
	fHistReconstSibFEED			(0),
	fHistReconstMixFEED			(0),
	fHistTriggerSib				(0),
	fHistTriggerMix				(0),
	fHistTriggerSibGEN			(0),
	fHistTriggerMixGEN			(0),
	fHistTriggerSibASO			(0),
	fHistTriggerMixASO			(0)					
{	

  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
}
//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation(const char *name)
   : AliAnalysisTaskSE(name),
	fAODEvent					(0x0),
	fPoolMgr					(0x0),
    fPIDResponse				(0x0),
	fAnalyseUE					(0x0),
	fPoolMaxNEvents				(0), 
	fPoolMinNTracks				(0), 
	fMinEventsToMix				(0),
	fNzVtxBins					(0), 
	fNCentBins					(0),
	fcollidingSys				(""),
	fpvzcut						(0),
    fTrackEtaCut				(0.9),
    fFilterBit					(128),
	fAnalysisMC					(0),
	fCase						(0),
	fRemoveAutoCorr				(0),
    fRapidityCut				(0),
	fV0radius					(0),
	fV0PostoPVz					(0),
	fV0NegtoPVz					(0),
	fDCAV0Daughters				(0),
	fCPAK0						(0),
	fCPALam						(0),
	fRejectLamK0				(0),
	fRejectK0Lam				(0),
    fSigmaPID					(0),
	fCutCTK0					(0),
	fCutCTLa					(0),
	fMassCutK0					(0),
	fMassCutLa					(0),
	fTriglow					(0),
	fTrighigh					(0),
	fTPCClusters				(0),					
	fTPCfindratio				(0),
	fUseChargeHadrons			(kTRUE), 
	fPtMin						(0.15),
	fOutputList					(0),
	fHist_Mult_B4_Trg_Sel		(0),
	fHist_Mult_Af_Trg_Sel		(0),
	fHist_Mult_PVz_Cut			(0),
	fHist_Mult_SPD_PVz			(0),
	fHist_Mult_SPD_PVz_Pileup	(0),
	fHistPVx					(0),
	fHistPVy					(0),
	fHistPVz					(0),
	fHistPVxAnalysis			(0),
	fHistPVyAnalysis			(0),
	fHistPVzAnalysis			(0),
	fHistEventViceGen			(0),
	fHistEventViceReconst		(0),
	fHistMCGenK0				(0),
	fHistMCGenLAM				(0),
	fHistMCGenALAM				(0),
	fHistMCGenLAMXIPLS			(0),
	fHistMCGenLAMXI 			(0),
	fHistReconstK0				(0),
	fHistReconstLA				(0),
	fHistReconstALA				(0),
	fHistMCAssoK0				(0),
	fHistMCAssoLA				(0),
	fHistMCAssoALA				(0),
	fHistMCAssoLAXI				(0),
	fHistMCAssoALAXiPlus		(0),
	fHistReconstSib				(0),
	fHistReconstMix				(0),
	fHistReconstSibGEN			(0),
	fHistReconstMixGEN			(0),
	fHistReconstSibASO			(0),
	fHistReconstMixASO			(0),
	fHistReconstSibFEED			(0),
	fHistReconstMixFEED			(0),
	fHistTriggerSib				(0),
	fHistTriggerMix				(0),
	fHistTriggerSibGEN			(0),
	fHistTriggerMixGEN			(0),
	fHistTriggerSibASO			(0),
	fHistTriggerMixASO			(0)


{	
  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  DefineOutput(1, TList::Class());                                            
}

//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::~AliLeadingV0Correlation()
{
   if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutputList;
   }
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::UserCreateOutputObjects()
{	
	fAnalyseUE =new AliAnalyseLeadingTrackUE();
	if(!fAnalysisMC)
	{
	fAnalyseUE->SetParticleSelectionCriteria(fFilterBit,fUseChargeHadrons,fTrackEtaCut,fPtMin);
	fAnalyseUE->DefineESDCuts(fFilterBit);
	}
	
	fOutputList = new TList();
	fOutputList->SetOwner();
	
	fHist_Mult_B4_Trg_Sel = new TH2F("fHist_Mult_B4_Trg_Sel","Tracks per event;Nbr of Tracks;Events", 1000, 0, 10000, 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_B4_Trg_Sel);
	
	fHist_Mult_Af_Trg_Sel = new TH2F("fHist_Mult_Af_Trg_Sel","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000, 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_Af_Trg_Sel);
	
	fHist_Mult_PVz_Cut = new TH2F("fHist_Mult_PVz_Cut","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000, 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_PVz_Cut);
	
	fHist_Mult_SPD_PVz = new TH2F("fHist_Mult_SPD_PVz","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000, 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_SPD_PVz);
	
	fHist_Mult_SPD_PVz_Pileup = new TH2F("fHist_Mult_SPD_PVz_Pileup","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000, 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_SPD_PVz_Pileup);
	
	fHistPVx = new TH1F("fHistPVx","PV x position;Nbr of Evts;x", 200, -0.5, 0.5); 		
	fOutputList->Add(fHistPVx);
	
	fHistPVy = new TH1F("fHistPVy","PV y position;Nbr of Evts;y",200, -0.5, 0.5); 		
	fOutputList->Add(fHistPVy);
	
	fHistPVz = new TH1F("fHistPVz","PV z position;Nbr of Evts;z",400, -20, 20); 		
	fOutputList->Add(fHistPVz);
	
	fHistPVxAnalysis = new TH1F("fHistPVxAnalysis","PV x position;Nbr of Evts;x", 200, -0.5, 0.5); 		
	fOutputList->Add(fHistPVxAnalysis);
	
	fHistPVyAnalysis = new TH1F("fHistPVyAnalysis","PV y position;Nbr of Evts;y",200, -0.5, 0.5); 		
	fOutputList->Add(fHistPVyAnalysis);
	
	fHistPVzAnalysis = new TH1F("fHistPVzAnalysis","PV z position;Nbr of Evts;z",400, -20, 20); 		
	fOutputList->Add(fHistPVzAnalysis);
	
	//---------------------------------------------- Events histograms -----------------------------------------------------//
	
	fHistEventViceGen= new TH2F("fHistEventViceGen", "fHistEventViceGen", 200, -20, 20, 10,0,1000);
	fOutputList->Add(fHistEventViceGen);
	
	fHistEventViceReconst= new TH2F("fHistEventViceReconst", "fHistEventViceReconst", 200, -20, 20, 10,0,1000);
	fOutputList->Add(fHistEventViceReconst);
	
	fHistMCGenLAM  = new TH2F("fHistMCGenLAM" , "fHistMCGenLAM" ,140,1.06,1.2, 120, 0, fTriglow);
	fOutputList->Add(fHistMCGenLAM);
	
	fHistMCGenALAM = new TH2F("fHistMCGenALAM", "fHistMCGenALAM",140,1.06,1.2, 120, 0, fTriglow);
	fOutputList->Add(fHistMCGenALAM);
	
	fHistMCGenK0   = new TH2F("fHistMCGenK0"  , "fHistMCGenK0"  ,200,0.4,0.6, 120, 0, fTriglow);
	fOutputList->Add(fHistMCGenK0);
	
	fHistMCGenLAMXIPLS = new TH2F("fHistMCGenLAMXIPLS", "fHistMCGenLAMXIPLS",140,1.06,1.2, 120, 0, fTriglow);
	fOutputList->Add(fHistMCGenLAMXIPLS);
	
	fHistMCGenLAMXI   = new TH2F("fHistMCGenLAMXI"  , "fHistMCGenLAMXI"  ,140,1.06,1.2, 120, 0, fTriglow);
	fOutputList->Add(fHistMCGenLAMXI);
	
	//New dimension for feed down corection 
	
	const Int_t ndimsK0 = 4;       
	Int_t    binsK0[ndimsK0] = {200, 120,500,1000};
	Double_t xminK0[ndimsK0] = {0.4,   0,  0,0.99};
	Double_t xmaxK0[ndimsK0] = {0.6,   fTriglow, 10,   1};
	
	const Int_t ndimsLA = 4;       
	Int_t    binsLA[ndimsLA] = { 140, 120,500,1000};
	Double_t xminLA[ndimsLA] = {1.06,   0,  0,0.99};
	Double_t xmaxLA[ndimsLA] = { 1.2,   fTriglow, 10,   1};
	
	fHistReconstK0= new THnSparseD("fHistReconstK0"  , "fHistReconstK0",ndimsK0,binsK0,xminK0,xmaxK0);
	fHistReconstK0->Sumw2();
	fOutputList->Add(fHistReconstK0);
	
	fHistReconstLA= new THnSparseD("fHistReconstLA"  , "fHistReconstLA",ndimsLA,binsLA,xminLA,xmaxLA);
	fHistReconstLA->Sumw2();
	fOutputList->Add(fHistReconstLA);
	
	fHistReconstALA= new THnSparseD("fHistReconstALA", "fHistReconstALA",ndimsLA,binsLA,xminLA,xmaxLA);
	fHistReconstALA->Sumw2();
	fOutputList->Add(fHistReconstALA);
	
	fHistMCAssoK0= new THnSparseD("fHistMCAssoK0"   , "fHistMCAssoK0"   ,ndimsK0,binsK0,xminK0,xmaxK0);
	fHistMCAssoK0->Sumw2();
	fOutputList->Add(fHistMCAssoK0);
	
	fHistMCAssoLA= new THnSparseD("fHistMCAssoLA"   , "fHistMCAssoLA"   ,ndimsLA,binsLA,xminLA,xmaxLA);
	fHistMCAssoLA->Sumw2();
	fOutputList->Add(fHistMCAssoLA);
	
	fHistMCAssoALA= new THnSparseD("fHistMCAssoALA" , "fHistMCAssoALA" , ndimsLA,binsLA,xminLA,xmaxLA);
	fHistMCAssoALA->Sumw2();
	fOutputList->Add(fHistMCAssoALA);
	
	fHistMCAssoLAXI= new THnSparseD("fHistMCAssoLAXI" , "fHistMCAssoLAXI" , ndimsLA,binsLA,xminLA,xmaxLA);
	fHistMCAssoLAXI->Sumw2();
	fOutputList->Add(fHistMCAssoLAXI);
	
	fHistMCAssoALAXiPlus= new THnSparseD("fHistMCAssoALAXiPlus" , "fHistMCAssoALAXiPlus" , ndimsLA,binsLA,xminLA,xmaxLA);
	fHistMCAssoALAXiPlus->Sumw2();
	fOutputList->Add(fHistMCAssoALAXiPlus);
	
	//--------------------------------------------Correlation Histos -----------------------------------------------------//
	
	//0-pTK0,1-PhiK0,2-EtaK0,3-DPhiK0,4-DEtaK0,5-TYPE,6-CutSet
	const Int_t ndimsv0CORR = 8;       
	Int_t    binsv0CORR[ndimsv0CORR] = {120, 200,          200,CorrBinsX,      CorrBinsY,4,500,1000};
	
	Double_t xminv0CORR[ndimsv0CORR] = {  0,   0,-fTrackEtaCut,    -PI/2,-2*fTrackEtaCut,0,  0,0.99};
	
	Double_t xmaxv0CORR[ndimsv0CORR] = {  fTriglow,2*PI, fTrackEtaCut,   3*PI/2, 2*fTrackEtaCut,4, 10,   1};
	
	fHistReconstSib= new THnSparseD("fHistReconstSib", "fHistReconstSib", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstSib->Sumw2();
	fOutputList->Add(fHistReconstSib);
	
	fHistReconstMix= new THnSparseD("fHistReconstMix", "fHistReconstMix", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstMix->Sumw2();
	fOutputList->Add(fHistReconstMix);
	
	fHistReconstSibGEN= new THnSparseD("fHistReconstSibGEN", "fHistReconstSibGEN", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstSibGEN->Sumw2();
	fOutputList->Add(fHistReconstSibGEN);
	
	fHistReconstMixGEN= new THnSparseD("fHistReconstMixGEN", "fHistReconstMixGEN", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstMixGEN->Sumw2();
	fOutputList->Add(fHistReconstMixGEN);
	
	fHistReconstSibASO= new THnSparseD("fHistReconstSibASO", "fHistReconstSibASO", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstSibASO->Sumw2();
	fOutputList->Add(fHistReconstSibASO);
	
	fHistReconstMixASO= new THnSparseD("fHistReconstMixASO", "fHistReconstMixASO", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstMixASO->Sumw2();
	fOutputList->Add(fHistReconstMixASO);
	
	fHistReconstSibFEED= new THnSparseD("fHistReconstSibFEED", "fHistReconstSibFEED", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstSibFEED->Sumw2();
	fOutputList->Add(fHistReconstSibFEED);
	
	fHistReconstMixFEED= new THnSparseD("fHistReconstMixFEED", "fHistReconstMixFEED", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstMixFEED->Sumw2();
	fOutputList->Add(fHistReconstMixFEED);
	
	
	fHistTriggerSib= new TH3F("fHistTriggerSib", "fHistTriggerSib", 100, fTriglow, fTrighigh,200,0,2*PI,200,-fTrackEtaCut,fTrackEtaCut);
	fHistTriggerSib->Sumw2();
	fOutputList->Add(fHistTriggerSib);
	
	fHistTriggerMix= new TH1F("fHistTriggerMix", "fHistTriggerMix", 100, fTriglow, fTrighigh);
	fHistTriggerMix->Sumw2();
	fOutputList->Add(fHistTriggerMix);
	
	fHistTriggerSibGEN= new TH3F("fHistTriggerSibGEN", "fHistTriggerSibGEN", 100, fTriglow, fTrighigh,200,0,2*PI,200,-fTrackEtaCut,fTrackEtaCut);
	fHistTriggerSibGEN->Sumw2();
	fOutputList->Add(fHistTriggerSibGEN);
	
	fHistTriggerMixGEN= new TH1F("fHistTriggerMixGEN", "fHistTriggerMixGEN", 100, fTriglow, fTrighigh);
	fHistTriggerMixGEN->Sumw2();
	fOutputList->Add(fHistTriggerMixGEN);
	
	fHistTriggerSibASO= new TH3F("fHistTriggerSibASO", "fHistTriggerSibASO", 100, fTriglow, fTrighigh,200,0,2*PI,200,-fTrackEtaCut,fTrackEtaCut);
	fHistTriggerSibASO->Sumw2();
	fOutputList->Add(fHistTriggerSibASO);
	
	fHistTriggerMixASO= new TH1F("fHistTriggerMixASO", "fHistTriggerMixASO", 100, fTriglow, fTrighigh);
	fHistTriggerMixASO->Sumw2();
	fOutputList->Add(fHistTriggerMixASO);
	
	//----------------------------------------------Event Pool-----------------------------------------------------//
	fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
	if(!fPoolMgr) return;
	
	PostData(1, fOutputList);
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::UserExec(Option_t *)
{
	
    AliAnalysisManager   *mgr      = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inEvMain = (AliInputEventHandler*)(mgr->GetInputEventHandler());
	if (!inEvMain) return;
	
	// Pointers to PID Response objects.	
	fPIDResponse = inEvMain->GetPIDResponse();
	if(!fPIDResponse) return;
	
    fAODEvent = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
	if(!fAODEvent) return;
	
	Int_t  ltrackMultiplicity        = 0;
	Int_t  lrefMultiplicity          = 0;

	//------------------------------------------------
	// Before Physics Selection
	//------------------------------------------------ 
	ltrackMultiplicity   = (InputEvent())->GetNumberOfTracks();
        AliAODHeader * header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
	lrefMultiplicity     = header->GetRefMultiplicity();
	
	fHist_Mult_B4_Trg_Sel->Fill(ltrackMultiplicity,lrefMultiplicity);
	
	Double_t * CentBins = fCentBins;
	Double_t poolmin    = CentBins[0];
	Double_t poolmax    = CentBins[fNCentBins];
	
	//----------------------------------------------------------
	// Efficency denomenator comes before the physics selection
	//----------------------------------------------------------
	
	Double_t  dimEventviceMC[2];
	if(fAnalysisMC)    //Efficency denomenator comes before the physics selection
	{
		AliAODMCHeader *aodMCheader = (AliAODMCHeader*)fAODEvent->FindListObject(AliAODMCHeader::StdBranchName());
		if(!aodMCheader) return;
		Float_t mcZv = aodMCheader->GetVtxZ();
		
		if (TMath::Abs(mcZv) >= fpvzcut) return;
		
		dimEventviceMC[0]=aodMCheader->GetVtxZ();
		
		TClonesArray *mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray)return;
		
		Int_t nMCTracks = mcArray->GetEntriesFast();
		
		dimEventviceMC[1]=nMCTracks;
		fHistEventViceGen->Fill(dimEventviceMC[0],dimEventviceMC[1]);
		
		TObjArray *selectedTracksLeadingMC=fAnalyseUE->FindLeadingObjects(mcArray);
		if(!selectedTracksLeadingMC) return;
		selectedTracksLeadingMC->SetOwner(kTRUE);
		
		TObjArray * selectedV0sMC =new TObjArray;
		selectedV0sMC->SetOwner(kTRUE);
		
		TObjArray * selectedV0sMCXI =new TObjArray;
		selectedV0sMCXI->SetOwner(kTRUE);
		
		for (Int_t iMC = 0; iMC<nMCTracks; iMC++)
		{
			AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(iMC);
			if (!mcTrack) continue;
			// Charged track Generated level
			Double_t mcTrackPt  = mcTrack->Pt();
			if ((mcTrackPt<fPtMin)||(mcTrackPt>fTriglow)) continue;
			
			Double_t mcTrackEta = mcTrack->Eta();
			Double_t mcTrackPhi = mcTrack->Phi();
			Bool_t TrIsPrime    = mcTrack->IsPhysicalPrimary();
			Bool_t TrPtMin      = mcTrackPt>fPtMin;
			Bool_t TrCharge     = (mcTrack->Charge())!=0;
			
			if (!TrIsPrime  && !TrPtMin && !TrCharge) continue;  //Check Point 1
			
			// V0 Generated level
			Int_t mcPartPdg		  = mcTrack->GetPdgCode();
			
			Double_t mcRapidity   = mcTrack->Y();
			Bool_t V0RapMax       = TMath::Abs(mcRapidity)<fRapidityCut;
			Bool_t V0EtaMax       = TMath::Abs(mcTrackEta)<fTrackEtaCut;
			Double_t mcMass       = mcTrack->M();
			
			Double_t mcK0[3] = {mcMass,mcTrackPt,static_cast<Double_t>(nMCTracks)};
			Double_t mcLa[3] = {mcMass,mcTrackPt,static_cast<Double_t>(nMCTracks)};
			Double_t mcAl[3] = {mcMass,mcTrackPt,static_cast<Double_t>(nMCTracks)};
			
			Int_t myTrackMotherLabel = mcTrack->GetMother();
			
			AliAODMCParticle *mcMother = (AliAODMCParticle*)mcArray->At(myTrackMotherLabel);
			if (!mcMother) continue;
			Int_t MotherPdg            = mcMother->GetPdgCode();
			Bool_t IsK0         = mcPartPdg==310;
			Bool_t IsLambda     = mcPartPdg==3122;
			Bool_t IsAntiLambda = mcPartPdg==-3122;
			
			Bool_t IsXImin  =MotherPdg== 3312;
			Bool_t IsXIPlus =MotherPdg==-3312;
			Bool_t IsXizero =MotherPdg== 3322;
			Bool_t IsOmega  =MotherPdg== 3334;
			
			switch (fCase) {
				case 1:
					
					if (IsK0) 
					{
						fHistMCGenK0->Fill(mcK0[0],mcK0[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
					} 
					
					if (IsLambda) 
					{
						fHistMCGenLAM->Fill(mcLa[0],mcLa[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda) 
					{	
						fHistMCGenALAM->Fill(mcAl[0],mcAl[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
					}
					
					if (IsLambda && (IsXizero || IsXImin)) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
						fHistMCGenLAMXI->Fill(mcLa[0],mcLa[1]);
					}
					
					if (IsLambda && IsOmega) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda && IsXIPlus) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
						fHistMCGenLAMXIPLS->Fill(mcAl[0],mcAl[1]);
					}
					
					break;
					
				case 2:

					if (IsK0 && V0RapMax && TrIsPrime) 
					{
						fHistMCGenK0->Fill(mcK0[0],mcK0[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
					} 
					
					if (IsLambda && V0RapMax && TrIsPrime) 
					{
						fHistMCGenLAM->Fill(mcLa[0],mcLa[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda && V0RapMax && TrIsPrime) 
					{	
						fHistMCGenALAM->Fill(mcAl[0],mcAl[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
					}
					
					if (IsLambda && V0RapMax && (IsXizero || IsXImin)) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
						fHistMCGenLAMXI->Fill(mcLa[0],mcLa[1]);
					}
					
					if (IsLambda && V0RapMax && IsOmega) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda && V0RapMax && IsXIPlus) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
						fHistMCGenLAMXIPLS->Fill(mcAl[0],mcAl[1]);
					}
					
					break;
					
				case 3:

					if (IsK0 && V0EtaMax && TrIsPrime) 
					{
						fHistMCGenK0->Fill(mcK0[0],mcK0[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
					} 
					
					if (IsLambda && V0EtaMax && TrIsPrime) 
					{
						fHistMCGenLAM->Fill(mcLa[0],mcLa[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda && V0EtaMax && TrIsPrime) 
					{	
						fHistMCGenALAM->Fill(mcAl[0],mcAl[1]);
						selectedV0sMC->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
					}
					
					if (IsLambda && V0EtaMax && (IsXizero || IsXImin)) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,1,0,0));
						fHistMCGenLAMXI->Fill(mcLa[0],mcLa[1]);
					}
					
					if (IsLambda && V0EtaMax && IsOmega) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,2,0,0));
					}
					
					if (IsAntiLambda && V0EtaMax && IsXIPlus) 
					{	
						selectedV0sMCXI->Add(new V0Correlationparticle(mcTrackEta,mcTrackPhi,mcTrackPt,3,0,0));
						fHistMCGenLAMXIPLS->Fill(mcAl[0],mcAl[1]);
					}
					break;
					
				default:
					AliInfo(Form("No case selected"));
					break;
			}
		}
		
		FillCorrelationSibling(nMCTracks,selectedTracksLeadingMC,selectedV0sMC,fHistTriggerSibGEN,fHistReconstSibGEN);
		FillCorrelationMixing(nMCTracks,mcZv,poolmax,poolmin,selectedTracksLeadingMC,selectedV0sMC,fHistTriggerMixGEN,fHistReconstMixGEN);
		
		FillCorrelationSibling(nMCTracks,selectedTracksLeadingMC,selectedV0sMCXI,0,fHistReconstSibFEED);
		FillCorrelationMixing(nMCTracks,mcZv,poolmax,poolmin,selectedTracksLeadingMC,selectedV0sMCXI,0,fHistReconstMixFEED);
    }
	
	// End Loop over MC condition
	
	//------------------------------------------------
	// Physics Selection
	//------------------------------------------------ 
	UInt_t maskIsSelected = inEvMain->IsEventSelected();
	Bool_t isSelected = ((maskIsSelected & AliVEvent::kMB)== AliVEvent::kMB);
    if (!isSelected) return;
	
	//------------------------------------------------
	// After Trigger Selection
	//------------------------------------------------
	
	fHist_Mult_Af_Trg_Sel->Fill(ltrackMultiplicity,lrefMultiplicity);
	
	//------------------------------------------------
	// Getting: Primary Vertex + MagField Info
	//------------------------------------------------
	Double_t  dimEventviceReal[3];
	Double_t  lBestPrimaryVtxPos[3];
	Double_t  tPrimaryVtxPosition[3];
	Double_t  lV0Position[3];
	
	AliAODVertex *lPrimaryBestAODVtx = fAODEvent->GetPrimaryVertex();
	if (!lPrimaryBestAODVtx) return;
	// get the best primary vertex available for the event
	// As done in AliCascadeVertexer, we keep the one which is the best one available.
	// between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
	// This one will be used for next calculations (DCA essentially)
	lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);
	
	const AliVVertex *primaryVtx = fAODEvent->GetPrimaryVertex();
	if(!primaryVtx)return;
	tPrimaryVtxPosition[0] = primaryVtx->GetX();
	tPrimaryVtxPosition[1] = primaryVtx->GetY();
	tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	fHistPVx->Fill( tPrimaryVtxPosition[0] );
	fHistPVy->Fill( tPrimaryVtxPosition[1] );
	fHistPVz->Fill( tPrimaryVtxPosition[2] );
	
	//------------------------------------------------
	// Primary Vertex Z position: SKIP
	//------------------------------------------------
	
	Double_t lPVx = lBestPrimaryVtxPos[0];
	Double_t lPVy = lBestPrimaryVtxPos[1];
	Double_t lPVz = lBestPrimaryVtxPos[2];
	
	if ((TMath::Abs(lPVz)) >= fpvzcut) return ;
	if (TMath::Abs(lPVx)<10e-5 && TMath::Abs(lPVy)<10e-5 && TMath::Abs(lPVz)<10e-5) return;
	fHist_Mult_PVz_Cut->Fill(ltrackMultiplicity,lrefMultiplicity);
	
	//------------------------------------------------
	// Only look at events with well-established PV
	//------------------------------------------------
	
	const AliAODVertex *lPrimaryTrackingAODVtxCheck = fAODEvent->GetPrimaryVertex();
	const AliAODVertex *lPrimarySPDVtx = fAODEvent->GetPrimaryVertexSPD();
	if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtxCheck )return;
	
	fHist_Mult_SPD_PVz->Fill(ltrackMultiplicity,lrefMultiplicity);
	//------------------------------------------------
	// Pileup Rejection
	//------------------------------------------------
	
	// FIXME : quality selection regarding pile-up rejection 
	if(fAODEvent->IsPileupFromSPD()) return;
	fHist_Mult_SPD_PVz_Pileup->Fill(ltrackMultiplicity,lrefMultiplicity);
	
	fHistPVxAnalysis->Fill(tPrimaryVtxPosition[0]);
	fHistPVyAnalysis->Fill(tPrimaryVtxPosition[1]);
	fHistPVzAnalysis->Fill(tPrimaryVtxPosition[2]);
	
    dimEventviceReal[0]=tPrimaryVtxPosition[2];
	dimEventviceReal[1]=ltrackMultiplicity;
	
	fHistEventViceReconst->Fill(dimEventviceReal[0],dimEventviceReal[1]);

	//---------------------------------------------------------------------------------------------
	
	Double_t lDcaPosToPrimVertex = 0;Double_t lDcaNegToPrimVertex = 0;Double_t lDcaV0Daughters     = 0;
	Double_t lV0cosPointAngle    = 0;Double_t lV0DecayLength      = 0;Double_t lV0Radius           = 0;
	Double_t lcTauLambda         = 0;Double_t lcTauAntiLambda     = 0;   
	Double_t lcTauK0s            = 0; 
	Double_t lDCAV0PVz           = 0; 
	
	Double_t lInvMassK0   = 0, lInvMassLambda    = 0, lInvMassAntiLambda = 0;
	Double_t lPtV0s       = 0; Double_t lPhiV0s  = 0; Double_t lEtaV0s   = 0;
	Double_t lRapK0       = 0, lRapLambda        = 0, lRapAntiLambda     = 0;
	Double_t lPzV0s       = 0; 
	Double_t lPV0s        = 0;
	
	TObjArray *selectedTracksLeading=0;
	selectedTracksLeading=fAnalyseUE->FindLeadingObjects(fAODEvent);
	if(!selectedTracksLeading) return;
	selectedTracksLeading->SetOwner(kTRUE);
	
	TObjArray * selectedV0s = new TObjArray;
	selectedV0s->SetOwner(kTRUE);
	
	TObjArray * selectedV0sAssoc = new TObjArray;
	selectedV0sAssoc->SetOwner(kTRUE);
	
	Int_t nV0s = fAODEvent->GetNumberOfV0s();
	
	for (Int_t i = 0; i < nV0s; i++) 
	{ // start of V0 slection loop
		AliAODv0* aodV0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(i));
		if (!aodV0) continue;
		
		if (((aodV0->Pt())<fPtMin)||((aodV0->Pt())>fTriglow)) continue;
		
		// get daughters
   	    AliAODTrack *myTrackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        AliAODTrack *myTrackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
		
		if (!myTrackPos || !myTrackNeg) continue;
		
        if (!IsAcseptedV0(aodV0,myTrackPos,myTrackNeg)) continue;
		
		// VO's main characteristics to check the reconstruction cuts
		lDcaV0Daughters    = aodV0->DcaV0Daughters();
		lV0cosPointAngle   = aodV0->CosPointingAngle(lBestPrimaryVtxPos);
		
		aodV0->GetXYZ(lV0Position);
		
		lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
		lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - tPrimaryVtxPosition[0],2) +
									 TMath::Power(lV0Position[1] - tPrimaryVtxPosition[1],2) +
									 TMath::Power(lV0Position[2] - tPrimaryVtxPosition[2],2));
		
		// DCA between daughter and Primary Vertex:
		if (myTrackPos) lDcaPosToPrimVertex = aodV0->DcaPosToPrimVertex();
		if (myTrackNeg) lDcaNegToPrimVertex = aodV0->DcaNegToPrimVertex();   
		lDCAV0PVz   = aodV0->DcaV0ToPrimVertex(); 
		
		// Quality tracks cuts:
		if ( !(IsAcseptedDaughterTrack(myTrackPos)) || !(IsAcseptedDaughterTrack(myTrackNeg)) ) { continue;}
		
		// Invariant mass
		lInvMassK0         = aodV0->MassK0Short();
		lInvMassLambda     = aodV0->MassLambda();
		lInvMassAntiLambda = aodV0->MassAntiLambda();
		
		lPtV0s = aodV0->Pt();
		lPhiV0s= aodV0->Phi();
		lEtaV0s= aodV0->Eta();
		lPzV0s = aodV0->Pz();
		
		// Rapidity:
		lRapK0     = aodV0->RapK0Short();
		lRapLambda = aodV0->RapLambda();
		lRapAntiLambda = aodV0->Y(-3122);
		
		if (lPtV0s==0) {continue;}
		
        Float_t nSigmaPosPion   = 0.;
        Float_t nSigmaNegPion   = 0.;
        Float_t nSigmaPosProton = 0.;
        Float_t nSigmaNegProton = 0.;
		
        const AliAODPid *pPid = myTrackPos->GetDetPid();
        const AliAODPid *nPid = myTrackNeg->GetDetPid();
		
        if (pPid)
        {
            Double_t pdMom = pPid->GetTPCmomentum();
            if (pdMom<1.0 && (fcollidingSys=="PbPb"))
            {
                nSigmaPosPion   = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
			
			if (fcollidingSys=="PP")
            {
                nSigmaPosPion   = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
        }
		
        if (nPid)
        {
            Double_t ndMom = nPid->GetTPCmomentum();
            if (ndMom<1.0 && (fcollidingSys=="PbPb"))
            {
                nSigmaNegPion   = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion);
                nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton);
            }
			
			if (fcollidingSys=="PP")
            {
                nSigmaNegPion   = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion);
                nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton);
            }
        }
		Bool_t bpPion   = TMath::Abs(nSigmaPosPion)   <= fSigmaPID;
        Bool_t bpProton = TMath::Abs(nSigmaPosProton) <= fSigmaPID;
        Bool_t bnPion   = TMath::Abs(nSigmaNegPion)   <= fSigmaPID;
        Bool_t bnProton = TMath::Abs(nSigmaNegProton) <= fSigmaPID;
		
        Bool_t cutK0Pid         = (bpPion   && bnPion)  ;
        Bool_t cutLambdaPid     = (bpProton && bnPion)  ;
        Bool_t cutAntiLambdaPid = (bpPion   && bnProton);
        //--------------------------------------------------
		
		lPV0s = TMath::Sqrt(lPzV0s*lPzV0s + lPtV0s*lPtV0s);
		
		if(lPV0s > 0) lcTauLambda     = (lV0DecayLength*lInvMassLambda)/lPV0s;
		if(lPV0s > 0) lcTauAntiLambda = (lV0DecayLength*lInvMassAntiLambda)/lPV0s; 
		if(lPV0s > 0) lcTauK0s        = (lV0DecayLength*lInvMassK0)/lPV0s;	
		
		Bool_t k0ctcut = (lcTauK0s        < fCutCTK0);
		Bool_t lactcut = (lcTauLambda     < fCutCTLa);
		Bool_t alactcut= (lcTauAntiLambda < fCutCTLa);
		
		Bool_t k0Rapcut = (TMath::Abs(lRapK0)         < fRapidityCut);
		Bool_t laRapcut = (TMath::Abs(lRapLambda)     < fRapidityCut);
		Bool_t alaRapcut= (TMath::Abs(lRapAntiLambda) < fRapidityCut);
		
		Bool_t V0EtaMax= (TMath::Abs(lEtaV0s) < fTrackEtaCut);
		
		Bool_t k0cutset = IsAcseptedK0(lV0Radius,lDcaPosToPrimVertex,lDcaNegToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassLambda,lInvMassAntiLambda);
		Bool_t lacutset = IsAcseptedLA(lV0Radius,lDcaPosToPrimVertex,lDcaNegToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassK0);
		Bool_t alacutset= IsAcseptedLA(lV0Radius,lDcaNegToPrimVertex,lDcaPosToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassK0);
		
		Double_t spK0[4] = {lInvMassK0,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
		Double_t spLa[4] = {lInvMassLambda,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
		Double_t spAl[4] = {lInvMassAntiLambda,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
	
		switch (fCase) {
			case 1:
				fHistReconstK0->Fill(spK0); 
				if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
				
				fHistReconstLA->Fill(spLa); 
				if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
				
				fHistReconstALA->Fill(spAl);
				if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
				
				break;
				
			case 2:
				if(k0ctcut && k0Rapcut && k0cutset && cutK0Pid)
				{
					fHistReconstK0->Fill(spK0); 
					if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
				}
				
				if (lactcut && laRapcut && lacutset && cutLambdaPid)
				{
					fHistReconstLA->Fill(spLa); 
					if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
				}
				
				if (alactcut && alaRapcut && alacutset && cutAntiLambdaPid)
				{
					fHistReconstALA->Fill(spAl);
					if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
				}

				break;
				
			case 3:
				if(k0ctcut && V0EtaMax && k0cutset && cutK0Pid)
				{
					fHistReconstK0->Fill(spK0); 
					if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
				}
				
				if (lactcut && V0EtaMax && lacutset && cutLambdaPid)
				{
					fHistReconstLA->Fill(spLa); 
					if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
				}
				
				if (alactcut && V0EtaMax && alacutset && cutAntiLambdaPid)
				{
					fHistReconstALA->Fill(spAl);
					if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
				}
				break;
				
			default:
				AliInfo(Form("No case selected"));
				break;
		}
		
    	if (fAnalysisMC)
    	{
      		TClonesArray *mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
      		if(!mcArray)return;
      		
			Int_t myTrackPosLabel        = TMath::Abs(myTrackPos->GetLabel());
			Int_t myTrackNegLabel        = TMath::Abs(myTrackNeg->GetLabel());
			
			AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
			if(!mcPosTrack)continue;
			AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
			if(!mcNegTrack)continue;
			
			Int_t PosDaughterPdg = mcPosTrack->GetPdgCode();
			Int_t NegDaughterPdg = mcNegTrack->GetPdgCode();
			
			Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
			Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();
			
			if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) continue;
			if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;
			
			AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
			if(!mcPosMother)continue;
			Int_t MotherPdg  = mcPosMother->GetPdgCode();
			Bool_t IsPrime   = mcPosMother->IsPhysicalPrimary();
			
			Int_t myGrandMotherLabel = mcPosMother->GetMother();
			AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)mcArray->At(myGrandMotherLabel);
			Int_t GrandMotherPdg     = mcGrandMother->GetPdgCode();
			
			Double_t rcK0[4] = {lInvMassK0,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
			Double_t rcLa[4] = {lInvMassLambda,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
			Double_t rcAl[4] = {lInvMassAntiLambda,lPtV0s,lDCAV0PVz,lV0cosPointAngle};
			
			switch (fCase) {
				case 1:
					fHistMCAssoK0->Fill(rcK0);
					if(IsK0InvMass(lInvMassK0))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
					
					fHistMCAssoLA->Fill(rcLa);
					if(IsLambdaInvMass(lInvMassLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
					
					fHistMCAssoALA->Fill(rcAl);
					if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
					
					break;
					
				case 2:
					if ((k0ctcut && k0Rapcut && k0cutset)&&(MotherPdg     ==  310 && 
															PosDaughterPdg==  211 && 
															NegDaughterPdg== -211 &&
															IsPrime))
					{
						fHistMCAssoK0->Fill(rcK0);
						if(IsK0InvMass(lInvMassK0))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((lactcut && laRapcut && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															IsPrime)) 
					{
						fHistMCAssoLA->Fill(rcLa);
						if(IsLambdaInvMass(lInvMassLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((alactcut && alaRapcut && alacutset)&&(MotherPdg     == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   IsPrime))
					{
						fHistMCAssoALA->Fill(rcAl);
						if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((lactcut && laRapcut && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															(GrandMotherPdg==3322 ||GrandMotherPdg==3312))) 
					{
						fHistMCAssoLAXI->Fill(rcLa);
					}
					
					if ((alactcut && alaRapcut && alacutset)&&(MotherPdg      == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   GrandMotherPdg== -3312))
					{
						fHistMCAssoALAXiPlus->Fill(rcAl);
					}
					
					break;
					
				case 3:
					if ((k0ctcut && V0EtaMax && k0cutset)&&(MotherPdg     ==  310 && 
																	   PosDaughterPdg==  211 && 
																	   NegDaughterPdg== -211 &&
																	   IsPrime))
					{
						fHistMCAssoK0->Fill(rcK0); 
						if(IsK0InvMass(lInvMassK0))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((lactcut && V0EtaMax && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															IsPrime)) 
					{
						fHistMCAssoLA->Fill(rcLa);
						if(IsLambdaInvMass(lInvMassLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((alactcut && V0EtaMax && alacutset)&&(MotherPdg     == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   IsPrime))
					{
						fHistMCAssoALA->Fill(rcAl);
						if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0sAssoc->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3,lDCAV0PVz,lV0cosPointAngle));
					}
					
					if ((lactcut && V0EtaMax && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															(GrandMotherPdg==3322 ||GrandMotherPdg==3312))) 
					{
						fHistMCAssoLAXI->Fill(rcLa);
					}
					
					if ((alactcut && V0EtaMax && alacutset)&&(MotherPdg      == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   GrandMotherPdg== -3312))
					{
						fHistMCAssoALAXiPlus->Fill(rcAl);
					}
					break;
					
				default:
					AliInfo(Form("No case selected"));
					break;
			}	
    	}
	} 	
	
	FillCorrelationSibling(ltrackMultiplicity,selectedTracksLeading,selectedV0s,fHistTriggerSib,fHistReconstSib);
	FillCorrelationMixing(ltrackMultiplicity,tPrimaryVtxPosition[2],poolmax,poolmin,selectedTracksLeading,selectedV0s,fHistTriggerMix,fHistReconstMix);
	
	FillCorrelationSibling(ltrackMultiplicity,selectedTracksLeading,selectedV0sAssoc,fHistTriggerSibASO,fHistReconstSibASO);
	FillCorrelationMixing(ltrackMultiplicity,lPVz,poolmax,poolmin,selectedTracksLeading,selectedV0sAssoc,fHistTriggerMixASO,fHistReconstMixASO);
	
	PostData(1,fOutputList);
}	
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::Terminate(Option_t *)
{
	//No need in the grid
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedDaughterTrack(const AliAODTrack *itrack)
{	
	if(fCase==1 || fCase==2)
	if(TMath::Abs(itrack->Eta())>fTrackEtaCut)return kFALSE;
	
	if (!itrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	
	Float_t nCrossedRowsTPC = itrack->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < fTPCClusters) return kFALSE;
	
	Int_t findable=itrack->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	
	if (nCrossedRowsTPC/findable < fTPCfindratio) return kFALSE;
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedV0(const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg)
{
	if (!aodV0) return kFALSE;
	
	// Offline reconstructed V0 only
    if (aodV0->GetOnFlyStatus()) return kFALSE;
	
    // Get daughters and check them
	myTrackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
	myTrackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
	
	if (!myTrackPos||!myTrackNeg) return kFALSE;
	// Unlike signs of daughters
    if (myTrackPos->Charge() == myTrackNeg->Charge()) return kFALSE;
	
	// Track cuts for daughers
    if ( !(IsAcseptedDaughterTrack(myTrackPos)) || !(IsAcseptedDaughterTrack(myTrackNeg)) ) return kFALSE;
	
	// Minimum pt of daughters
    Double_t lPtPos = myTrackPos->Pt();
    Double_t lPtNeg = myTrackNeg->Pt();
	
	if (lPtPos<fPtMin || lPtNeg<fPtMin) return kFALSE;
	
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedK0(Double_t v0rad,
							Double_t dcaptp,
							Double_t dcantp,
							Double_t dcav0d,
							Double_t cpa,
							Double_t massLa,
							Double_t massALa)
{	
			if(v0rad  >=fV0radius		&&
			   dcaptp >=fV0PostoPVz		&&
			   dcantp >=fV0NegtoPVz		&&
			   dcav0d <=fDCAV0Daughters	&&
			   cpa    >=fCPAK0			&&
			   TMath::Abs(massLa  - 1.115683) > fRejectLamK0 &&
			   TMath::Abs(massALa - 1.115683) > fRejectLamK0 )return kTRUE;
	return kFALSE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedLA(Double_t v0rad,
							Double_t dcaptp,
							Double_t dcantp,
							Double_t dcav0d,
							Double_t cpa,
							Double_t massK0)
{
	if(v0rad  >=fV0radius		&&
	   dcaptp >=fV0PostoPVz		&&
	   dcantp >=fV0NegtoPVz		&&
	   dcav0d <=fDCAV0Daughters	&&
	   cpa    >=fCPALam			&&
	   TMath::Abs(massK0  - 0.4976) > fRejectK0Lam &&
	   TMath::Abs(massK0  - 0.4976) > fRejectK0Lam )return kTRUE;
	return kFALSE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsK0InvMass(const Double_t mass) const 
{
	const Float_t massK0            = 0.497; 
	
	return ((massK0-fMassCutK0)<=mass && mass<=(massK0 + fMassCutK0))?1:0;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsLambdaInvMass(const Double_t mass) const 
{
	const Float_t massLambda        = 1.116; 
	
	return ((massLambda-fMassCutLa)<=mass && mass<=(massLambda + fMassCutLa))?1:0;
}
//---------------------------------------------------------------------------------------
Double_t AliLeadingV0Correlation::RangePhi(Double_t DPhi)
{
	if (DPhi < -TMath::Pi()/2)  DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	
	return DPhi;	
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsTrackFromV0(AliAODTrack* track)
{
	Int_t atrID = track->GetID();

	for(int i=0; i<fAODEvent->GetNumberOfV0s(); i++){ // loop over V0s
		AliAODv0* aodV0 = fAODEvent->GetV0(i);
		
		AliAODTrack *trackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        AliAODTrack *trackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
			
		if ( !(IsAcseptedDaughterTrack(trackPos)) || !(IsAcseptedDaughterTrack(trackNeg)) ) continue;
		//----------------------------------
		Int_t negID = trackNeg->GetID();
		Int_t posID = trackPos->GetID();
		
		if ((TMath::Abs(negID)+1)==(TMath::Abs(atrID))){ return kTRUE;}
		if ((TMath::Abs(posID)+1)==(TMath::Abs(atrID))){ return kTRUE;}
		//----------------------------------
	}
	return kFALSE;
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::FillCorrelationSibling(Double_t MultipOrCent,
									  TObjArray*triggerArray,
									  TObjArray*selectedV0Array,
									  TH3F*triggerHist,
									  THnSparse*associateHist)
{
	Double_t  binsv0CORR[8];
	Double_t  binsTrigSib[2];
	Int_t counterSibMCA=0;
	
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++)
	{
		AliAODTrack* trigger = (AliAODTrack*)triggerArray->At(0);
		if(!trigger)continue;
		
		if(fRemoveAutoCorr) 
		if(IsTrackFromV0(trigger))continue;
			
		Double_t triggerPt  = trigger->Pt();
		Double_t triggerPhi = trigger->Phi();
		Double_t triggerEta = trigger->Eta();
		
		if(triggerPt<fTriglow||triggerPt>fTrighigh)continue;
		counterSibMCA++;
		
		if(counterSibMCA==triggerArray->GetEntriesFast()){
			
			binsTrigSib[0]=triggerPt;
			binsTrigSib[1]=MultipOrCent;
			
			if(triggerHist)triggerHist->Fill(binsTrigSib[0],triggerPhi,triggerEta);
			
			for (Int_t j=0; j<selectedV0Array->GetEntriesFast(); j++){
				
				V0Correlationparticle* associate = (V0Correlationparticle*) selectedV0Array->At(j);
				if(!associate)continue;
				
				binsv0CORR[0]= associate->Pt();
				binsv0CORR[1]= associate->Phi();
				binsv0CORR[2]= associate->Eta();
				
				if(binsv0CORR[0]>triggerPt) continue;
				
				binsv0CORR[3]=RangePhi(triggerPhi-binsv0CORR[1]);
				binsv0CORR[4]=triggerEta-binsv0CORR[2];
				binsv0CORR[5]= associate->WhichCandidate();
				binsv0CORR[6]= associate->DCAPostoP();
				binsv0CORR[7]= associate->DCANegtoP();
				
				associateHist->Fill(binsv0CORR);
			}
		  }
		}
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::FillCorrelationMixing(Double_t MultipOrCentMix,
								   Double_t pvxMix,
								   Double_t poolmax,
								   Double_t poolmin,
								   TObjArray*triggerArray,
								   TObjArray*selectedV0Array,
								   TH1F*triggerHist,
								   THnSparse*associateHist)
{
	if(TMath::Abs(pvxMix)>=fpvzcut || MultipOrCentMix>poolmax || MultipOrCentMix < poolmin)
	{
		if(fcollidingSys=="PP")AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",pvxMix,MultipOrCentMix));
		return;
	}
	
	Double_t  binsv0CORRMix[8];
	Double_t  binsTrigMix[2];
	Double_t  counterMix=0;
	
	AliEventPool* pool = fPoolMgr->GetEventPool(MultipOrCentMix, pvxMix);
	if (!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", MultipOrCentMix, pvxMix));
	
    if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks  || pool->GetCurrentNEvents() > fMinEventsToMix)
	{
		Int_t nMix = pool->GetCurrentNEvents();
		for (Int_t jMix=0; jMix<nMix; jMix++){
			
			TObjArray* mixEvents = pool->GetEvent(jMix);
			for (Int_t i=0; i<triggerArray->GetEntriesFast(); i++){
				
				AliAODTrack* trig = (AliAODTrack*)triggerArray->At(0);
				if(!trig)continue;
				
				if(fRemoveAutoCorr) 
				if(IsTrackFromV0(trig))continue;
				
				Double_t trigPhi  = trig->Phi();
				Double_t trigEta  = trig->Eta();
				Double_t trigPt   = trig->Pt();
				
				if(trigPt<fTriglow||trigPt>fTrighigh)continue;
				counterMix++;
				
				if(counterMix==triggerArray->GetEntriesFast()){
					
					binsTrigMix[0]=trigPt;
					binsTrigMix[1]=MultipOrCentMix;
					
					if(triggerHist)triggerHist->Fill(binsTrigMix[0]);
					
					for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
						
						V0Correlationparticle* associate = (V0Correlationparticle*) mixEvents->At(j);
						if(!associate)continue;
						
						binsv0CORRMix[0]= associate->Pt();
						binsv0CORRMix[1]= associate->Phi();
						binsv0CORRMix[2]= associate->Eta();
						
						if(binsv0CORRMix[0]>trigPt) continue;
						
						binsv0CORRMix[3]=RangePhi(trigPhi-binsv0CORRMix[1]);
						binsv0CORRMix[4]=trigEta-binsv0CORRMix[2];
						binsv0CORRMix[5]=associate->WhichCandidate();
						binsv0CORRMix[6]=associate->DCAPostoP();
						binsv0CORRMix[7]=associate->DCANegtoP();
						
						associateHist->Fill(binsv0CORRMix);
						}
					}
				}
			}
		}
	
	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);
	
	for (Int_t i=0; i<selectedV0Array->GetEntriesFast(); i++)
	{
		V0Correlationparticle* particle = (V0Correlationparticle*) selectedV0Array->At(i);
		tracksClone->Add(new V0Correlationparticle(particle->Eta(), 
											  particle->Phi(), 
											  particle->Pt(),
											  particle->WhichCandidate(),
											  particle->DCAPostoP(),
											  particle->DCANegtoP()));
	};
	pool->UpdatePool(tracksClone);
}
//---------------------------------------------------------------------------------------											
								
					
												
												
