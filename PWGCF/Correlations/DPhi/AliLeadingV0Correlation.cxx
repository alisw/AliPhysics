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

#define CorrBinsX 24
#define CorrBinsY 26


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
	fHistReconstK0				(0),
	fHistReconstLA				(0),
	fHistReconstALA				(0),
	fHistMCAssoK0				(0),
	fHistMCAssoLA				(0),
	fHistMCAssoALA				(0),
	fHistReconstSib				(0),
	fHistReconstMix				(0),
	fHistTriggerSib				(0),
	fHistTriggerMix				(0)
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
	fHistReconstK0				(0),
	fHistReconstLA				(0),
	fHistReconstALA				(0),
	fHistMCAssoK0				(0),
	fHistMCAssoLA				(0),
	fHistMCAssoALA				(0),
	fHistReconstSib				(0),
	fHistReconstMix				(0),
	fHistTriggerSib				(0),
	fHistTriggerMix				(0)

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
	
	fHist_Mult_B4_Trg_Sel = new TH1F("fHist_Mult_B4_Trg_Sel","Tracks per event;Nbr of Tracks;Events", 1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_B4_Trg_Sel);
	
	fHist_Mult_Af_Trg_Sel = new TH1F("fHist_Mult_Af_Trg_Sel","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_Af_Trg_Sel);
	
	fHist_Mult_PVz_Cut = new TH1F("fHist_Mult_PVz_Cut","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_PVz_Cut);
	
	fHist_Mult_SPD_PVz = new TH1F("fHist_Mult_SPD_PVz","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000); 		
	fOutputList->Add(fHist_Mult_SPD_PVz);
	
	fHist_Mult_SPD_PVz_Pileup = new TH1F("fHist_Mult_SPD_PVz_Pileup","Tracks per event;Nbr of Tracks;Events",1000, 0, 10000); 		
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
	//0-PVx,1-PVy,2-PVz,3-MULT,4-CENT
	const Int_t ndimsEV = 3;
	Int_t    binsEV[ndimsEV] = { 200,  100,100};
	Double_t xminEV[ndimsEV] = {-20 ,    0,  0};
	Double_t xmaxEV[ndimsEV] = { 20 ,  300,100};
	
	fHistEventViceGen= new THnSparseD("fHistEventViceGen", "fHistEventViceGen", ndimsEV, binsEV, xminEV, xmaxEV);
	fOutputList->Add(fHistEventViceGen);
	
	fHistEventViceReconst= new THnSparseD("fHistEventViceReconst", "fHistEventViceReconst", ndimsEV, binsEV, xminEV, xmaxEV);
	fOutputList->Add(fHistEventViceReconst);
	
	//0-YK0,1-Pt
	const Int_t ndimsGenMC = 3;
	Int_t    binsGenMCLA[ndimsGenMC] = {120, 140,100};
	Double_t xminGenMCLA[ndimsGenMC] = {  0,1.06,  0};
	Double_t xmaxGenMCLA[ndimsGenMC] = {  6, 1.2,100};  
	
	Int_t    binsGenMCK0[ndimsGenMC] = {120, 200,100};
	Double_t xminGenMCK0[ndimsGenMC] = {  0, 0.4,  0};
	Double_t xmaxGenMCK0[ndimsGenMC] = {  6, 0.6,100};
	
	fHistMCGenLAM  = new THnSparseD("fHistMCGenLAM" , "fHistMCGenLAM" , ndimsGenMC, binsGenMCLA, xminGenMCLA, xmaxGenMCLA);
	fOutputList->Add(fHistMCGenLAM);
	
	fHistMCGenALAM = new THnSparseD("fHistMCGenALAM", "fHistMCGenALAM", ndimsGenMC, binsGenMCLA, xminGenMCLA, xmaxGenMCLA);
	fOutputList->Add(fHistMCGenALAM);
	
	fHistMCGenK0   = new THnSparseD("fHistMCGenK0"  , "fHistMCGenK0"  , ndimsGenMC, binsGenMCK0, xminGenMCK0, xmaxGenMCK0);
	fOutputList->Add(fHistMCGenK0);
	
	const Int_t ndims=3;    //MK0  mLA  MALA PT   cent
	Int_t    binsK0[ndims] = {  200, 120  ,100};
	Double_t xminK0[ndims] = {  0.4,   0  ,  0};
	Double_t xmaxK0[ndims] = {  0.6,   6  ,100};
	 
	Int_t    binsLA[ndims] = {  140, 120  ,100};
	Double_t xminLA[ndims] = { 1.06,   0  ,  0};
	Double_t xmaxLA[ndims] = {  1.2,   6  ,100};
	
	
	fHistReconstK0= new THnSparseD("fHistReconstK0"  , "fHistReconstK0",  ndims, binsK0, xminK0, xmaxK0);
	fHistReconstK0->Sumw2();
	fOutputList->Add(fHistReconstK0);
	
	fHistReconstLA= new THnSparseD("fHistReconstLA"  , "fHistReconstLA",  ndims, binsLA, xminLA, xmaxLA);
	fHistReconstLA->Sumw2();
	fOutputList->Add(fHistReconstLA);
	
	fHistReconstALA= new THnSparseD("fHistReconstALA", "fHistReconstALA", ndims, binsLA, xminLA, xmaxLA);
	fHistReconstALA->Sumw2();
	fOutputList->Add(fHistReconstALA);
	
	fHistMCAssoK0= new THnSparseD("fHistMCAssoK0"   , "fHistMCAssoK0"   , ndims, binsK0, xminK0, xmaxK0);
	fHistMCAssoK0->Sumw2();
	fOutputList->Add(fHistMCAssoK0);
	
	fHistMCAssoLA= new THnSparseD("fHistMCAssoLA"   , "fHistMCAssoLA"   , ndims, binsLA, xminLA, xmaxLA);
	fHistMCAssoLA->Sumw2();
	fOutputList->Add(fHistMCAssoLA);
	
	fHistMCAssoALA= new THnSparseD("fHistMCAssoALA" , "fHistMCAssoALA" ,  ndims, binsLA, xminLA, xmaxLA);
	fHistMCAssoALA->Sumw2();
	fOutputList->Add(fHistMCAssoALA);
	
	//--------------------------------------------Correlation Histos -----------------------------------------------------//
	
	//0-pTK0,1-PhiK0,2-EtaK0,3-DPhiK0,4-DEtaK0,5-TYPE,6-CutSet
	const Int_t ndimsv0CORR = 7;       
	Int_t    binsv0CORR[ndimsv0CORR] = {120, 200,          200,CorrBinsX,      CorrBinsY,4,100};
	
	Double_t xminv0CORR[ndimsv0CORR] = {  0,   0,-fTrackEtaCut,    -PI/2,-2*fTrackEtaCut,0,0};
	
	Double_t xmaxv0CORR[ndimsv0CORR] = {  6,2*PI, fTrackEtaCut,   3*PI/2, 2*fTrackEtaCut,4,100};
	
	fHistReconstSib= new THnSparseD("fHistReconstSib", "fHistReconstSib", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstSib->Sumw2();
	fOutputList->Add(fHistReconstSib);
	
	fHistReconstMix= new THnSparseD("fHistReconstMix", "fHistReconstMix", ndimsv0CORR, binsv0CORR, xminv0CORR, xmaxv0CORR);
	fHistReconstMix->Sumw2();
	fOutputList->Add(fHistReconstMix);
	
	//0-pt,1-PHI,2-Eta
	const Int_t triggerdims       =2;
	Int_t binsTrig[triggerdims]   ={       100, 100};
	Double_t xminTrig[triggerdims]={  fTriglow,   0};
	Double_t xmaxTrig[triggerdims]={ fTrighigh, 100};
	
	fHistTriggerSib= new THnSparseD("fHistTriggerSib", "fHistTriggerSib", triggerdims, binsTrig, xminTrig, xmaxTrig);
	fHistTriggerSib->Sumw2();
	fOutputList->Add(fHistTriggerSib);
	
	fHistTriggerMix= new THnSparseD("fHistTriggerMix", "fHistTriggerMix", triggerdims, binsTrig, xminTrig, xmaxTrig);
	fHistTriggerMix->Sumw2();
	fOutputList->Add(fHistTriggerMix);
	
	
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
	
	Int_t multiplicity    = -1;
	Int_t multiplicityMC  = -1;
	Double_t MultipOrCent = -1; 
	Double_t CentPecentMC = -1;
	Double_t CentPecentAfterPhySel    = -1;
	Int_t    nTrackMultiplicity       = -1;
	Float_t lPrimaryTrackMultiplicity = 0;
	
	nTrackMultiplicity              = (InputEvent())->GetNumberOfTracks();
    for (Int_t itrack = 0; itrack<nTrackMultiplicity; itrack++) {
		AliAODTrack* track = fAODEvent->GetTrack(itrack);
		if(!fAnalysisMC) if (track->TestFilterBit(fFilterBit)) lPrimaryTrackMultiplicity++;
		lPrimaryTrackMultiplicity++;
    }

	fHist_Mult_B4_Trg_Sel->Fill(lPrimaryTrackMultiplicity);
	
	if(fcollidingSys=="PbPb"){   
	AliCentrality *centralityObjMC = fAODEvent->GetHeader()->GetCentralityP();
		CentPecentMC  = centralityObjMC->GetCentralityPercentileUnchecked("V0M");
		if ((CentPecentMC < 0.)||(CentPecentMC > 90)) return;
	}
	
	Double_t * CentBins = fCentBins;
	Double_t poolmin    = CentBins[0];
	Double_t poolmax    = CentBins[fNCentBins];
	
	//----------------------------------------------------------
	// Efficency denomenator comes before the physics selection
	//----------------------------------------------------------
	
	Double_t  dimEventviceMC[3];
	if(fAnalysisMC)    //Efficency denomenator comes before the physics selection
	{
		AliAODMCHeader *aodMCheader = (AliAODMCHeader*)fAODEvent->FindListObject(AliAODMCHeader::StdBranchName());
		Float_t mcZv = aodMCheader->GetVtxZ();
		
		if (TMath::Abs(mcZv) >= fpvzcut) return;
		
		dimEventviceMC[0]=aodMCheader->GetVtxZ();
		
		TClonesArray *mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray)return;
		
		Int_t nMCTracks = mcArray->GetEntriesFast();
		
		if(fcollidingSys=="PbPb") multiplicityMC=CentPecentMC;
		if(fcollidingSys=="PP")   multiplicityMC=nMCTracks;
		
		dimEventviceMC[1]=nMCTracks;
		dimEventviceMC[2]=CentPecentMC;
		fHistEventViceGen->Fill(dimEventviceMC);
		
		for (Int_t iMC = 0; iMC<nMCTracks; iMC++)
		{
			AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(iMC);
			if (!mcTrack) continue;
			// Charged track Generated level
			Double_t mcTrackPt  = mcTrack->Pt();
			if ((mcTrackPt<fPtMin)||(mcTrackPt>6.0)) continue;
			
			Bool_t TrIsPrime    = mcTrack->IsPhysicalPrimary();
			Bool_t TrPtMin      = mcTrackPt>fPtMin;
			Bool_t TrCharge     = (mcTrack->Charge())!=0;
			
			if (!TrIsPrime  && !TrPtMin && !TrCharge) continue;  //Check Point 1
			
			// V0 Generated level
			Int_t mcPartPdg		  = mcTrack->GetPdgCode();
			
			Double_t mcRapidity   = mcTrack->Y();
			Bool_t V0RapMax       = TMath::Abs(mcRapidity)<fRapidityCut;
			Double_t mcMass       = mcTrack->M();
			
			Double_t mcK0[3] = {mcTrackPt,mcMass,multiplicityMC};
			Double_t mcLa[3] = {mcTrackPt,mcMass,multiplicityMC};
			Double_t mcAl[3] = {mcTrackPt,mcMass,multiplicityMC};
			
			
			Bool_t IsK0 = mcPartPdg==310;
			if (IsK0 && V0RapMax && TrIsPrime) 
			{
				fHistMCGenK0->Fill(mcK0);
			} 
			
			Bool_t IsLambda = mcPartPdg==3122;
			if (IsLambda && V0RapMax && TrIsPrime) 
			{
				fHistMCGenLAM->Fill(mcLa);
			}
			
			Bool_t IsAntiLambda = mcPartPdg==-3122;
			if (IsAntiLambda && V0RapMax && TrIsPrime) 
			{	
				fHistMCGenALAM->Fill(mcAl);
			}			
		}
    }
	
	// End Loop over MC condition
	
	//------------------------------------------------
	// Physics Selection
	//------------------------------------------------ 
	UInt_t maskIsSelected = inEvMain->IsEventSelected();
	Bool_t isSelected = ((maskIsSelected & AliVEvent::kMB)== AliVEvent::kMB 
					  || (maskIsSelected & AliVEvent::kCentral)== AliVEvent::kCentral 
					  || (maskIsSelected & AliVEvent::kSemiCentral)== AliVEvent::kSemiCentral);
    if (!isSelected) return;
	
	//------------------------------------------------
	// After Trigger Selection
	//------------------------------------------------
	
	fHist_Mult_Af_Trg_Sel->Fill(lPrimaryTrackMultiplicity);
	
	//------------------------------------------------
	// Getting: Primary Vertex + MagField Info
	//------------------------------------------------
	Double_t  dimEventviceReal[3];
	Double_t  lBestPrimaryVtxPos[3];
	Double_t  tPrimaryVtxPosition[3];
	Double_t  lV0Position[3];
	
	
	if(fcollidingSys=="PbPb"){  //
	AliCentrality *centralityObj = fAODEvent->GetHeader()->GetCentralityP();
		CentPecentAfterPhySel  = centralityObj->GetCentralityPercentileUnchecked("V0M");
		if ((CentPecentAfterPhySel < 0.)||(CentPecentAfterPhySel > 90)) return;
	} //
	
	AliAODVertex *lPrimaryBestAODVtx = fAODEvent->GetPrimaryVertex();
	if (!lPrimaryBestAODVtx) return;
	// get the best primary vertex available for the event
	// As done in AliCascadeVertexer, we keep the one which is the best one available.
	// between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
	// This one will be used for next calculations (DCA essentially)
	lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);
	
	const AliVVertex *primaryVtx = fAODEvent->GetPrimaryVertex();
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
	fHist_Mult_PVz_Cut->Fill(lPrimaryTrackMultiplicity);
	
	//------------------------------------------------
	// Only look at events with well-established PV
	//------------------------------------------------
	
	const AliAODVertex *lPrimaryTrackingAODVtxCheck = fAODEvent->GetPrimaryVertex();
	const AliAODVertex *lPrimarySPDVtx = fAODEvent->GetPrimaryVertexSPD();
	if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtxCheck )return;
	
	fHist_Mult_SPD_PVz->Fill(lPrimaryTrackMultiplicity);

	
	//------------------------------------------------
	// Pileup Rejection
	//------------------------------------------------
	
	// FIXME : quality selection regarding pile-up rejection 
	if(fAODEvent->IsPileupFromSPD()) return;
	fHist_Mult_SPD_PVz_Pileup->Fill(lPrimaryTrackMultiplicity);
	
	fHistPVxAnalysis->Fill(tPrimaryVtxPosition[0]);
	fHistPVyAnalysis->Fill(tPrimaryVtxPosition[1]);
	fHistPVzAnalysis->Fill(tPrimaryVtxPosition[2]);
	
    dimEventviceReal[0]=tPrimaryVtxPosition[2];
	multiplicity       = fAODEvent->GetNTracks();
	
	dimEventviceReal[1]=multiplicity;
	dimEventviceReal[2]=CentPecentAfterPhySel;
	
	fHistEventViceReconst->Fill(dimEventviceReal);
	
	if(fcollidingSys=="PP")MultipOrCent=multiplicity;
	if(fcollidingSys=="PbPb")MultipOrCent=CentPecentAfterPhySel;

	//---------------------------------------------------------------------------------------------
	
	Double_t lDcaPosToPrimVertex = 0;Double_t lDcaNegToPrimVertex = 0;Double_t lDcaV0Daughters     = 0;
	Double_t lV0cosPointAngle    = 0;Double_t lV0DecayLength      = 0;Double_t lV0Radius           = 0;
	Double_t lcTauLambda         = 0;Double_t lcTauAntiLambda     = 0;   
	Double_t lcTauK0s            = 0;   
	
	Double_t lInvMassK0   = 0, lInvMassLambda    = 0, lInvMassAntiLambda = 0;
	Double_t lPtV0s       = 0; Double_t lPhiV0s  = 0; Double_t lEtaV0s   = 0;
	Double_t lRapK0s      = 0, lRapLambda        = 0, lRapAntiLambda     = 0;
	Double_t lPzV0s       = 0; Double_t lAlphaV0 = 0, lPtArmV0           = 0;
	Double_t lPV0s        = 0;
	
	TObjArray *selectedTracksLeading=0;
	selectedTracksLeading=fAnalyseUE->FindLeadingObjects(fAODEvent);
	if(!selectedTracksLeading) return;
	selectedTracksLeading->SetOwner(kTRUE);
	
	TObjArray * selectedV0s = new TObjArray;
	selectedV0s->SetOwner(kTRUE);
	
	Int_t nV0s = fAODEvent->GetNumberOfV0s();
	
	for (Int_t i = 0; i < nV0s; i++) 
	{ // start of V0 slection loop
		AliAODv0* aodV0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(i));
		if (!aodV0) continue;
		
		if (((aodV0->Pt())<fPtMin)||((aodV0->Pt())>6.0)) continue;
		
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
		
		// Quality tracks cuts:
		if ( !(IsAcseptedDaughterTrack(myTrackPos)) || !(IsAcseptedDaughterTrack(myTrackNeg)) ) { continue;}
		
		// Armenteros variables:
		lAlphaV0      =  aodV0->AlphaV0();
		lPtArmV0      =  aodV0->PtArmV0();
		
		// Invariant mass
		lInvMassK0         = aodV0->MassK0Short();
		lInvMassLambda     = aodV0->MassLambda();
		lInvMassAntiLambda = aodV0->MassAntiLambda();
		
		lPtV0s = aodV0->Pt();
		lPhiV0s= aodV0->Phi();
		lEtaV0s= aodV0->Eta();
		lPzV0s = aodV0->Pz();
		
		// Rapidity:
		lRapK0s    = aodV0->RapK0Short();
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

		Bool_t k0APcut = (lPtArmV0>(TMath::Abs(0.2*lAlphaV0)));
		
		Bool_t k0Rapcut = (TMath::Abs(lRapK0s)        < fRapidityCut);
		Bool_t laRapcut = (TMath::Abs(lRapLambda)     < fRapidityCut);
		Bool_t alaRapcut= (TMath::Abs(lRapAntiLambda) < fRapidityCut);
		
		if(fcollidingSys=="PbPb")if(lV0Radius>=100) continue;
		
		Bool_t k0cutset = IsAcseptedK0(lV0Radius,lDcaPosToPrimVertex,lDcaNegToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassLambda,lInvMassAntiLambda);
		Bool_t lacutset = IsAcseptedLA(lV0Radius,lDcaPosToPrimVertex,lDcaNegToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassK0);
		Bool_t alacutset= IsAcseptedLA(lV0Radius,lDcaNegToPrimVertex,lDcaPosToPrimVertex,lDcaV0Daughters,lV0cosPointAngle,lInvMassK0);
		
		Double_t spK0[3] = {lInvMassK0, lPtV0s,MultipOrCent};
		Double_t spLa[3] = {lInvMassLambda,lPtV0s,MultipOrCent};
		Double_t spAl[3] = {lInvMassAntiLambda,lPtV0s,MultipOrCent};
	
		switch (fCase) {
			case 1:
				fHistReconstK0->Fill(spK0); 
				if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1));
				
				fHistReconstLA->Fill(spLa); 
				if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2));
				
				fHistReconstALA->Fill(spAl);
				if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3));
				
				break;
				
			case 2:
				if(k0ctcut && k0Rapcut && k0cutset && cutK0Pid)
				{
					fHistReconstK0->Fill(spK0); 
					if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1));
				}
				
				if (lactcut && laRapcut && lacutset && cutLambdaPid)
				{
					fHistReconstLA->Fill(spLa); 
					if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2));
				}
				
				if (alactcut && alaRapcut && alacutset && cutAntiLambdaPid)
				{
					fHistReconstALA->Fill(spAl);
					if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3));
				}

				break;
				
			case 3:
				if(k0ctcut && k0Rapcut && k0cutset && cutK0Pid && k0APcut)
				{
					fHistReconstK0->Fill(spK0); 
					if(IsK0InvMass(lInvMassK0))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,1));
				}
				
				if (lactcut && laRapcut && lacutset && cutLambdaPid)
				{
					fHistReconstLA->Fill(spLa); 
					if(IsLambdaInvMass(lInvMassLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,2));
				}
				
				if (alactcut && alaRapcut && alacutset && cutAntiLambdaPid)
				{
					fHistReconstALA->Fill(spAl);
					if(IsLambdaInvMass(lInvMassAntiLambda))selectedV0s->Add(new V0Correlationparticle(lEtaV0s,lPhiV0s,lPtV0s,3));
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
			AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
			
			Int_t PosDaughterPdg = mcPosTrack->GetPdgCode();
			Int_t NegDaughterPdg = mcNegTrack->GetPdgCode();
			
			Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
			Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();
			
			if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) continue;
			if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;
			
			AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
			Int_t MotherPdg  = mcPosMother->GetPdgCode();
			Bool_t IsPrime   = mcPosMother->IsPhysicalPrimary();
			
			Double_t rcK0[3] = {lInvMassK0, lPtV0s,MultipOrCent};
			Double_t rcLa[3] = {lInvMassLambda,lPtV0s,MultipOrCent};
			Double_t rcAl[3] = {lInvMassAntiLambda,lPtV0s,MultipOrCent};
			
			switch (fCase) {
				case 1:
					fHistMCAssoK0->Fill(rcK0); 
					fHistMCAssoLA->Fill(rcLa);
					fHistMCAssoALA->Fill(rcAl);
					
					break;
					
				case 2:
					if ((k0ctcut && k0Rapcut && k0cutset)&&(MotherPdg     ==  310 && 
															PosDaughterPdg==  211 && 
															NegDaughterPdg== -211 &&
															IsPrime))
					{
						fHistMCAssoK0->Fill(rcK0); 
					}
					
					if ((lactcut && laRapcut && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															IsPrime)) 
					{
						fHistMCAssoLA->Fill(rcLa);
					}
					
					if ((alactcut && alaRapcut && alacutset)&&(MotherPdg     == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   IsPrime))
					{
						fHistMCAssoALA->Fill(rcAl);
					}
					
					break;
					
				case 3:
					if ((k0ctcut && k0Rapcut && k0cutset && k0APcut)&&(MotherPdg     ==  310 && 
																	   PosDaughterPdg==  211 && 
																	   NegDaughterPdg== -211 &&
																	   IsPrime))
					{
						fHistMCAssoK0->Fill(rcK0); 
					}
					
					if ((lactcut && laRapcut && lacutset)&&(MotherPdg     == 3122 && 
															PosDaughterPdg== 2212 && 
															NegDaughterPdg== -211 &&
															IsPrime)) 
					{
						fHistMCAssoLA->Fill(rcLa);
					}
					
					if ((alactcut && alaRapcut && alacutset)&&(MotherPdg     == -3122 && 
															   PosDaughterPdg==   211 && 
															   NegDaughterPdg== -2212 &&
															   IsPrime))
					{
						fHistMCAssoALA->Fill(rcAl);
					}
					break;
					
				default:
					AliInfo(Form("No case selected"));
					break;
			}	
    	}
	} 	
	
	FillCorrelationSibling(MultipOrCent,selectedTracksLeading,selectedV0s,fHistTriggerSib,fHistReconstSib);
	FillCorrelationMixing(MultipOrCent,tPrimaryVtxPosition[2],poolmax,poolmin,selectedTracksLeading,selectedV0s,fHistTriggerMix,fHistReconstMix);
	
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
	if(TMath::Abs(itrack->Eta())>fTrackEtaCut)return kFALSE;
	
	if (!itrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	
	Float_t nCrossedRowsTPC = itrack->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < 70) return kFALSE;
	
	Int_t findable=itrack->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	
	if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
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
									  THnSparse*triggerHist,
									  THnSparse*associateHist)
{
	Double_t  binsv0CORR[7];
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
			
			if(triggerHist)triggerHist->Fill(binsTrigSib);
			
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
				binsv0CORR[6]= MultipOrCent;
				
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
								   THnSparse*triggerHist,
								   THnSparse*associateHist)
{
	if(TMath::Abs(pvxMix)>=fpvzcut || MultipOrCentMix>poolmax || MultipOrCentMix < poolmin)
	{
		if(fcollidingSys=="PP")AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",pvxMix,MultipOrCentMix));
		if(fcollidingSys=="PbPb") AliInfo(Form("PbPb Event with Zvertex = %.2f cm and centrality = %.1f  out of pool bounds, SKIPPING",pvxMix,MultipOrCentMix));
		return;
	}
	
	Double_t  binsv0CORRMix[7];
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
					
					if(triggerHist)triggerHist->Fill(binsTrigMix);
					
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
						binsv0CORRMix[6]=MultipOrCentMix;
						
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
											  particle->WhichCandidate()));
	};
	pool->UpdatePool(tracksClone);
}
//---------------------------------------------------------------------------------------											
								
					
												
												
