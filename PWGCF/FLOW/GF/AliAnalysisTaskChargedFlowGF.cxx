#include "AliAnalysisTaskChargedFlowGF.h"

// # include <TList.h>
// # include <TH1.h>
// # include <TH2.h>
// # include <TH3.h>
// # include <TProfile.h>
// # include <TComplex.h>
// # include <TBits.h> 
// AliRoot includes
// # include "AliESDEvent.h"
// # include "AliAODEvent.h"
// # include "AliVEvent.h"
// # include "AliVTrack.h"
// # include "AliVVertex.h"
// # include "AliAnalysisFilter.h"
// # include "AliESDtrackCuts.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskChargedFlowGF)
//___________________________________________________________________________
AliAnalysisTaskChargedFlowGF::AliAnalysisTaskChargedFlowGF():
	AliAnalysisTaskSE(),
	fAOD(0),
	fFilterbit(96),
	fEtaCut(0.8),
	fVtxCut(10.0),  
	fMinPt(0.2),
	fMaxPt(5.0),
	fTPCclusters(70),
	fUseDCAzCut(0),
	fDCAz(0),
	fUseDCAxyCut(0),
	fDCAxy(0),
	fSample(1),
	fTrigger(0),
	fNUE(0),
	fNUA(0),
		//..MC
	fIsMC(false),
		//....
	fPeriod(""),

	fListOfObjects(0),
	fListOfObjectsMC(0),

	fMultTOFLowCut(0),
	fMultTOFHighCut(0),
	fMultCentLowCut(0),

	fTrackEfficiency(0),
	hTrackEfficiency(0),
	hTrackEfficiencyRun(0),

	fPhiWeight(0),
	hPhiWeight(0),

	hEventCount(0),
	hMult(0), 
	fVtxAfterCuts(0),
	fCentralityDis(0),  
	fV0CentralityDis(0),  
	hMultV0vsNtrksAfterCuts(0),
	hNtrksVSmultPercentile(0),
	fCentralityV0MCL1(0),
	fCentralityV0MCL0(0),
	fCentralityCL0CL1(0),
  fMultvsCentr(0),
	fMult128vsCentr(0),
	fMultTPCvsTOF(0),
	fMultTPCvsESD(0),

	fPhiDis1D(0),
	fPhiDis(0), 
	fEtaDis(0),
	fEtaBefore(0),
	fPtDis(0),
	fPtBefore(0),
	hDCAxyBefore(0),
	hDCAzBefore(0),
	hITSclustersBefore(0),
	hChi2Before(0),
	hDCAxy(0),
	hDCAz(0),
	hITSclusters(0),
	hChi2(0),

	fChsc4242(0),
	fChsc4242Gap0(0),
	fChsc4242_3sub(0),
	fChsc4224_3sub(0),
	fChsc4242_3sub_perm2(0),
	fChsc4224_3sub_perm2(0),
	fChsc4242_3sub_perm3(0),
	fChsc4224_3sub_perm3(0),
	fChsc3232(0),
	fChsc3232Gap0(0),
	fChsc3232_3sub(0),
	fChsc3223_3sub(0),
	fChsc3232_3sub_perm2(0),
	fChsc3223_3sub_perm2(0),
	fChsc3232_3sub_perm3(0),
	fChsc3223_3sub_perm3(0),

		//..MC
	hMultMC(0),
	fPhiDisTruth(0),
	fEtaDisTruth(0),
	fPtDisTruth(0),

	hReco(0),
	hTruth(0),
	hNtrksRecoNtrksTruth(0),

	hDCAptMC(0),
	hDCAptMC_material(0),
	hDCAptMC_weak(0),

	fChMCsc4242(0),
	fChMCsc4242Gap0(0),
	fChMCsc4242_3sub(0),
	fChMCsc4224_3sub(0),
	fChMCsc3232(0),
	fChMCsc3232Gap0(0),
	fChMCsc3232_3sub(0),
	fChMCsc3223_3sub(0)
{
 
	for(int h=0; h<3; h++)
	{
 
		fChcn2Ntrks1bin[h] 									= 0;
		fChcn2Gap0Ntrks1bin[h] 							= 0;
		fChcn2Gap2Ntrks1bin[h] 							= 0;
		fChcn2Gap4Ntrks1bin[h] 							= 0;
		fChcn2Gap8Ntrks1bin[h] 							= 0;
		fChcn2GapNtrks1bin[h] 							= 0;
		fChcn2Gap14Ntrks1bin[h] 						= 0;

		fChcn2_3subLMNtrks1bin[h] 					= 0;
		fChcn2_3subRMNtrks1bin[h] 					= 0;
		fChcn2_3subLM_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subLR_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subRM_perm3Ntrks1bin[h] 		= 0;
		fChcn2_3subRL_perm3Ntrks1bin[h] 		= 0;

		fChcn4Ntrks1bin[h] 									= 0;
		fChcn4Gap0Ntrks1bin[h] 							= 0;
		fChcn4_3subNtrks1bin[h] 						= 0;
		fChcn4_3sub_perm2Ntrks1bin[h]				= 0;
		fChcn4_3sub_perm3Ntrks1bin[h]				= 0;

		fChcn6Ntrks1bin[h] 									= 0;
		fChcn6Gap0Ntrks1bin[h] 							= 0;

		fChcn8Ntrks1bin[h] 									= 0;
		fChcn8Gap0Ntrks1bin[h] 							= 0;

			//..MC truth
		fChMCcn2Ntrks1bin[h] 								= 0;
		fChMCcn2Gap0Ntrks1bin[h] 						= 0;
		fChMCcn2Gap2Ntrks1bin[h] 						= 0;
		fChMCcn2Gap4Ntrks1bin[h] 						= 0;
		fChMCcn2Gap8Ntrks1bin[h] 						= 0;
		fChMCcn2GapNtrks1bin[h] 						= 0;
		fChMCcn2Gap14Ntrks1bin[h] 					= 0;

		fChMCcn2_3subLMNtrks1bin[h] 				= 0;
		fChMCcn2_3subRMNtrks1bin[h] 				= 0;

		fChMCcn4Ntrks1bin[h] 								= 0;
		fChMCcn4Gap0Ntrks1bin[h] 						= 0;
		fChMCcn4_3subNtrks1bin[h] 					= 0;

		fChMCcn6Ntrks1bin[h] 								= 0;
		fChMCcn6Gap0Ntrks1bin[h] 						= 0;

		fChMCcn8Ntrks1bin[h] 								= 0;
		fChMCcn8Gap0Ntrks1bin[h] 						= 0;

 	}

	for(int k=0; k<fSample+1; k++)
	{			

		fsc4242[k] 													= 0;
		fsc4242Gap0[k] 											= 0;
		fsc4242_3sub[k] 										= 0;
		fsc4224_3sub[k] 										= 0;
		fsc4242_3sub_perm2[k] 						= 0;
		fsc4224_3sub_perm2[k] 						= 0;
		fsc4242_3sub_perm3[k]							= 0;
		fsc4224_3sub_perm3[k]							= 0;

		fsc3232[k] 													= 0;
		fsc3232Gap0[k] 											= 0;
		fsc3232_3sub[k] 										= 0;
		fsc3223_3sub[k] 										= 0;
		fsc3232_3sub_perm2[k] 						= 0;
		fsc3223_3sub_perm2[k] 						= 0;
		fsc3232_3sub_perm3[k] 						= 0;
		fsc3223_3sub_perm3[k] 						= 0;

		fMCsc4242[k] 												= 0;
		fMCsc4242Gap0[k] 										= 0;
		fMCsc4242_3sub[k] 									= 0;
		fMCsc4224_3sub[k] 									= 0;

		fMCsc3232[k] 												= 0;
		fMCsc3232Gap0[k] 										= 0;
		fMCsc3232_3sub[k] 									= 0;
		fMCsc3223_3sub[k] 									= 0;
	
		for(Int_t h = 0; h < 3; h++)	
		{		

			fcn2Ntrks1bin[h][k] 							= 0;
			fcn2GapNtrks1bin[h][k] 						= 0;
			fcn2Gap0Ntrks1bin[h][k] 					= 0;
			fcn2Gap2Ntrks1bin[h][k] 					= 0;
			fcn2Gap4Ntrks1bin[h][k] 					= 0;
			fcn2Gap8Ntrks1bin[h][k] 					= 0;
			fcn2Gap14Ntrks1bin[h][k] 					= 0;

			fcn2_3subLMNtrks1bin[h][k] 				= 0;
			fcn2_3subRMNtrks1bin[h][k] 				= 0;
			fcn2_3subLM_perm2Ntrks1bin[h][k] = 0;
			fcn2_3subLR_perm2Ntrks1bin[h][k] = 0;
			fcn2_3subRM_perm3Ntrks1bin[h][k] = 0;
			fcn2_3subRL_perm3Ntrks1bin[h][k] = 0;

			fcn4Ntrks1bin[h][k] 							= 0;
			fcn4Gap0Ntrks1bin[h][k] 					= 0;
			fcn4_3subNtrks1bin[h][k] 					= 0;	
			fcn4_3sub_perm2Ntrks1bin[h][k]		= 0;
			fcn4_3sub_perm3Ntrks1bin[h][k]		= 0;

			fcn6Ntrks1bin[h][k] 							= 0;
			fcn6Gap0Ntrks1bin[h][k] 					= 0;

			fcn8Ntrks1bin[h][k] 							= 0;
			fcn8Gap0Ntrks1bin[h][k] 					= 0;

			//..MC truth
			fMCcn2Ntrks1bin[h][k] 						= 0;
			fMCcn2GapNtrks1bin[h][k] 					= 0;
			fMCcn2Gap0Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap2Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap4Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap8Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap14Ntrks1bin[h][k] 				= 0;

			fMCcn2_3subLMNtrks1bin[h][k] 			= 0;
			fMCcn2_3subRMNtrks1bin[h][k] 			= 0;
	
			fMCcn4Ntrks1bin[h][k] 						= 0;
			fMCcn4Gap0Ntrks1bin[h][k] 				= 0;
			fMCcn4_3subNtrks1bin[h][k] 				= 0;

			fMCcn6Ntrks1bin[h][k] 						= 0;
			fMCcn6Gap0Ntrks1bin[h][k] 				= 0;
		
			fMCcn8Ntrks1bin[h][k] 						= 0;
			fMCcn8Gap0Ntrks1bin[h][k] 				= 0;

	 }
       
  }
    
}
//______________________________________________________________________________
AliAnalysisTaskChargedFlowGF::AliAnalysisTaskChargedFlowGF(const char *name):
  AliAnalysisTaskSE(name),
	fAOD(0),
	fFilterbit(96),
	fEtaCut(0.8),
	fVtxCut(10.0),  
	fMinPt(0.2),
	fMaxPt(5.0),
	fTPCclusters(70),
	fUseDCAzCut(0),
	fDCAz(0),
	fUseDCAxyCut(0),
	fDCAxy(0),
	fSample(1),
	fTrigger(0),
	fNUE(0),
	fNUA(0),
		//..MC
	fIsMC(false),
		//....
	fPeriod(""),

	fListOfObjects(0),
	fListOfObjectsMC(0),

	fMultTOFLowCut(0),
	fMultTOFHighCut(0),
	fMultCentLowCut(0),

	fTrackEfficiency(0),
	hTrackEfficiency(0),
	hTrackEfficiencyRun(0),

	fPhiWeight(0),
	hPhiWeight(0),

	hEventCount(0),
	hMult(0), 
	fVtxAfterCuts(0),
	fCentralityDis(0),  
	fV0CentralityDis(0),  
	hMultV0vsNtrksAfterCuts(0),
	hNtrksVSmultPercentile(0),
	fCentralityV0MCL1(0),
	fCentralityV0MCL0(0),
	fCentralityCL0CL1(0),
  fMultvsCentr(0),
	fMult128vsCentr(0),
	fMultTPCvsTOF(0),
	fMultTPCvsESD(0),

	fPhiDis1D(0),
	fPhiDis(0), 
	fEtaDis(0),
	fEtaBefore(0),
	fPtDis(0),
	fPtBefore(0),
	hDCAxyBefore(0),
	hDCAzBefore(0),
	hITSclustersBefore(0),
	hChi2Before(0),
	hDCAxy(0),
	hDCAz(0),
	hITSclusters(0),
	hChi2(0),

	fChsc4242(0),
	fChsc4242Gap0(0),
	fChsc4242_3sub(0),
	fChsc4224_3sub(0),
	fChsc4242_3sub_perm2(0),
	fChsc4224_3sub_perm2(0),
	fChsc4242_3sub_perm3(0),
	fChsc4224_3sub_perm3(0),
	fChsc3232(0),
	fChsc3232Gap0(0),
	fChsc3232_3sub(0),
	fChsc3223_3sub(0),
	fChsc3232_3sub_perm2(0),
	fChsc3223_3sub_perm2(0),
	fChsc3232_3sub_perm3(0),
	fChsc3223_3sub_perm3(0),

		//..MC
	hMultMC(0),
	fPhiDisTruth(0),
	fEtaDisTruth(0),
	fPtDisTruth(0),

	hReco(0),
	hTruth(0),
	hNtrksRecoNtrksTruth(0),

	hDCAptMC(0),
	hDCAptMC_material(0),
	hDCAptMC_weak(0),

	fChMCsc4242(0),
	fChMCsc4242Gap0(0),
	fChMCsc4242_3sub(0),
	fChMCsc4224_3sub(0),
	fChMCsc3232(0),
	fChMCsc3232Gap0(0),
	fChMCsc3232_3sub(0),
	fChMCsc3223_3sub(0)
{

	for(int h=0; h<3; h++)
	{
 
		fChcn2Ntrks1bin[h] 									= 0;
		fChcn2Gap0Ntrks1bin[h] 							= 0;
		fChcn2Gap2Ntrks1bin[h] 							= 0;
		fChcn2Gap4Ntrks1bin[h] 							= 0;
		fChcn2Gap8Ntrks1bin[h] 							= 0;
		fChcn2GapNtrks1bin[h] 							= 0;
		fChcn2Gap14Ntrks1bin[h] 						= 0;

		fChcn2_3subLMNtrks1bin[h] 					= 0;
		fChcn2_3subRMNtrks1bin[h] 					= 0;
		fChcn2_3subLM_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subLR_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subRM_perm3Ntrks1bin[h] 		= 0;
		fChcn2_3subRL_perm3Ntrks1bin[h] 		= 0;

		fChcn4Ntrks1bin[h] 									= 0;
		fChcn4Gap0Ntrks1bin[h] 							= 0;
		fChcn4_3subNtrks1bin[h] 						= 0;
		fChcn4_3sub_perm2Ntrks1bin[h]				= 0;
		fChcn4_3sub_perm3Ntrks1bin[h]				= 0;

		fChcn6Ntrks1bin[h] 									= 0;
		fChcn6Gap0Ntrks1bin[h] 							= 0;

		fChcn8Ntrks1bin[h] 									= 0;
		fChcn8Gap0Ntrks1bin[h] 							= 0;

			//..MC truth
		fChMCcn2Ntrks1bin[h] 								= 0;
		fChMCcn2Gap0Ntrks1bin[h] 						= 0;
		fChMCcn2Gap2Ntrks1bin[h] 						= 0;
		fChMCcn2Gap4Ntrks1bin[h] 						= 0;
		fChMCcn2Gap8Ntrks1bin[h] 						= 0;
		fChMCcn2GapNtrks1bin[h] 						= 0;
		fChMCcn2Gap14Ntrks1bin[h] 					= 0;

		fChMCcn2_3subLMNtrks1bin[h] 				= 0;
		fChMCcn2_3subRMNtrks1bin[h] 				= 0;

		fChMCcn4Ntrks1bin[h] 								= 0;
		fChMCcn4Gap0Ntrks1bin[h] 						= 0;
		fChMCcn4_3subNtrks1bin[h] 					= 0;

		fChMCcn6Ntrks1bin[h] 								= 0;
		fChMCcn6Gap0Ntrks1bin[h] 						= 0;

		fChMCcn8Ntrks1bin[h] 								= 0;
		fChMCcn8Gap0Ntrks1bin[h] 						= 0;

 	}

	for(int k=0; k<fSample+1; k++)
	{			

		fsc4242[k] 													= 0;
		fsc4242Gap0[k] 											= 0;
		fsc4242_3sub[k] 										= 0;
		fsc4224_3sub[k] 										= 0;
		fsc4242_3sub_perm2[k] 						= 0;
		fsc4224_3sub_perm2[k] 						= 0;
		fsc4242_3sub_perm3[k]							= 0;
		fsc4224_3sub_perm3[k]							= 0;

		fsc3232[k] 													= 0;
		fsc3232Gap0[k] 											= 0;
		fsc3232_3sub[k] 										= 0;
		fsc3223_3sub[k] 										= 0;
		fsc3232_3sub_perm2[k] 						= 0;
		fsc3223_3sub_perm2[k] 						= 0;
		fsc3232_3sub_perm3[k] 						= 0;
		fsc3223_3sub_perm3[k] 						= 0;

		fMCsc4242[k] 												= 0;
		fMCsc4242Gap0[k] 										= 0;
		fMCsc4242_3sub[k] 									= 0;
		fMCsc4224_3sub[k] 									= 0;

		fMCsc3232[k] 												= 0;
		fMCsc3232Gap0[k] 										= 0;
		fMCsc3232_3sub[k] 									= 0;
		fMCsc3223_3sub[k] 									= 0;
	
		for(Int_t h = 0; h < 3; h++)	
		{		

			fcn2Ntrks1bin[h][k] 							= 0;
			fcn2GapNtrks1bin[h][k] 						= 0;
			fcn2Gap0Ntrks1bin[h][k] 					= 0;
			fcn2Gap2Ntrks1bin[h][k] 					= 0;
			fcn2Gap4Ntrks1bin[h][k] 					= 0;
			fcn2Gap8Ntrks1bin[h][k] 					= 0;
			fcn2Gap14Ntrks1bin[h][k] 					= 0;

			fcn2_3subLMNtrks1bin[h][k] 				= 0;
			fcn2_3subRMNtrks1bin[h][k] 				= 0;
			fcn2_3subLM_perm2Ntrks1bin[h][k] 				= 0;
			fcn2_3subLR_perm2Ntrks1bin[h][k] 				= 0;
			fcn2_3subRM_perm3Ntrks1bin[h][k] 				= 0;
			fcn2_3subRL_perm3Ntrks1bin[h][k] 				= 0;

			fcn4Ntrks1bin[h][k] 							= 0;
			fcn4Gap0Ntrks1bin[h][k] 					= 0;
			fcn4_3subNtrks1bin[h][k] 					= 0;	
			fcn4_3sub_perm2Ntrks1bin[h][k]		= 0;
			fcn4_3sub_perm3Ntrks1bin[h][k]		= 0;

			fcn6Ntrks1bin[h][k] 							= 0;
			fcn6Gap0Ntrks1bin[h][k] 					= 0;

			fcn8Ntrks1bin[h][k]			 					= 0;
			fcn8Gap0Ntrks1bin[h][k] 					= 0;

			//..MC truth
			fMCcn2Ntrks1bin[h][k] 						= 0;
			fMCcn2GapNtrks1bin[h][k] 					= 0;
			fMCcn2Gap0Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap2Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap4Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap8Ntrks1bin[h][k] 				= 0;
			fMCcn2Gap14Ntrks1bin[h][k] 				= 0;

			fMCcn2_3subLMNtrks1bin[h][k] 			= 0;
			fMCcn2_3subRMNtrks1bin[h][k] 			= 0;
	
			fMCcn4Ntrks1bin[h][k] 						= 0;
			fMCcn4Gap0Ntrks1bin[h][k] 				= 0;
			fMCcn4_3subNtrks1bin[h][k] 				= 0;

			fMCcn6Ntrks1bin[h][k] 						= 0;
			fMCcn6Gap0Ntrks1bin[h][k] 				= 0;
		
			fMCcn8Ntrks1bin[h][k] 						= 0;
			fMCcn8Gap0Ntrks1bin[h][k] 				= 0;

	 }
       
  }
    
	// Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
	DefineOutput(2, TList::Class());
	
}

//_____________________________________________________________________________
AliAnalysisTaskChargedFlowGF::~AliAnalysisTaskChargedFlowGF()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects) 
    delete fListOfObjects;
	if(fListOfObjectsMC)
		delete fListOfObjectsMC;

}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlowGF::UserCreateOutputObjects()
{ 
	
  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

	if(fIsMC == true)
	{
  	//OpenFile(1);
  	fListOfObjectsMC = new TList();
  	fListOfObjectsMC->SetOwner();
	}

	//..defining LHC15o cuts
	fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
	fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);

	fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);

	fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
	fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);

	fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
		//........

	//..Settings for AliEventCuts: 
	//..This adds QA plots to the output
	fEventCuts.AddQAplotsToList(fListOfObjects);
	//..kINT7 is set in the class as default, if I want to have kHigHMultV0 in pp, I have to switch to manual mode
	if(fTrigger == 0){ //	kHighMultV0
		fEventCuts.SetManualMode();

  	fEventCuts.fRequireTrackVertex = true;
  	fEventCuts.fMinVtz = -10.f;
  	fEventCuts.fMaxVtz = 10.f;
  	fEventCuts.fMaxDeltaSpdTrackAbsolute = 0.5f;
  	fEventCuts.fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  	fEventCuts.fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  	fEventCuts.fMaxResolutionSPDvertex = 0.25f;

  	fEventCuts.fRejectDAQincomplete = true;

  	fEventCuts.fRequiredSolenoidPolarity = 0;

  	fEventCuts.fUseMultiplicityDependentPileUpCuts = true; // If user specify a value it is not overwritten
  	fEventCuts.fSPDpileupMinZdist = 0.8;
  	fEventCuts.fSPDpileupNsigmaZdist = 3.;
  	fEventCuts.fSPDpileupNsigmaDiamXY = 2.;
  	fEventCuts.fSPDpileupNsigmaDiamZ = 5.;
  	fEventCuts.fTrackletBGcut = true;

  	fEventCuts.fFB128vsTrklLinearCut[0] = 32.077;
  	fEventCuts.fFB128vsTrklLinearCut[1] = 0.932;

  	if(fTrigger == 0) fEventCuts.fTriggerMask = AliVEvent::kHighMultV0;

  	fEventCuts.fUtils.SetMaxPlpChi2MV(5);
  	fEventCuts.fUtils.SetMinWDistMV(15);
  	fEventCuts.fUtils.SetCheckPlpFromDifferentBCMV(kFALSE);
  	fEventCuts.fPileUpCutMV = true;
	}

	// range on Xaxis: 
	//		for pp, pPb and for lowM PbPb and XeXe the range is 200 in unit bins
	//		for PbPb and XeXe the range is 5000 in bin width 50 (low M part which won't be correct due to mult. fluctuations will be done from unit bins in separate running)
	int range = 0;
	int nbins = 0;
	if(fPeriod == "LHC15oLR" || fPeriod == "LHC15oHR" || fPeriod == "LHC17n") {
		range = 3500;
		nbins = range/50;
	}
	else {
		range = 200;
		nbins = 200;
	}
	const int nn = nbins;

	/////////////////////////////////////////////
	if(fPeriod == "LHC16k" || fPeriod == "LHC16l" || fPeriod == "LHC16o")
	{
		if(fIsMC == false)
		{
	
		//..tracking efficiency
			fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pp13TeV_%s.root", fPeriod.Data()));
			if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA", fPeriod.Data()));
			//	if any of the below cases are switched on (for systematics), then the histogram is overwritten
			if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_FB768", fPeriod.Data()));
			if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters100", fPeriod.Data()));

			//..phi weight
			fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_%s.root", fPeriod.Data()));
			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s", fPeriod.Data())); 
			// for systematic checks
			if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters100", fPeriod.Data()));

		}else{
			//	these are used for MC closure test
    	//..tracking efficiency
    	fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiencyMC_%s.root", fPeriod.Data()));
    	if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
    	  printf("file does not exist");
    	  return;
    	}    

    	hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_default", fPeriod.Data()));

			//..phi weights from MonteCarlo
			fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeightMC_%s.root", fPeriod.Data()));
			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_default");

		}

	}

	/////////////////////////////////////////////
	if(fPeriod == "LHC17m" || fPeriod == "LHC17o")
	{
		//..tracking efficiency
			fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2017/TrackingEfficiency_pp13TeV_%s.root", fPeriod.Data()));
			if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA", fPeriod.Data()));
			//	if any of the below cases are switched on (for systematics), then the histogram is overwritten
			if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_FB768", fPeriod.Data()));
			if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters100", fPeriod.Data()));

			//..phi weight
			fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2017/PhiWeight_%s.root", fPeriod.Data()));
			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s", fPeriod.Data())); 
			// for systematic checks
			if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters100", fPeriod.Data()));

	}

	////////////////////////////////////////////////
	if(fPeriod == "LHC16t_CENT_wSDD" || fPeriod == "LHC16t_CENT_woSDD" || fPeriod == "LHC16t_FAST")
 	{
		//	this data set doens't need run-by-run weight
		//	CENT_woSDD and FAST reconstructions are identical -> they have the same weights

		//..phi weight
		fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_LHC16t.root");
		if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
			printf("file does not exist");
			return;
		}

		if(fPeriod == "LHC16t_CENT_wSDD") 
		{
			//	default
			hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_wSDD");
			//	for systematics	
			if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_wSDD_TPCclusters80");
			if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_wSDD_TPCclusters90");
			if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_wSDD_TPCclusters100");
		}
		if(fPeriod == "LHC16t_CENT_woSDD" || fPeriod == "LHC16t_FAST")
		{
			//	default
			hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_woSDD_FAST");
			//	for systematics
			if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_woSDD_FAST_TPCclusters80");
			if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_woSDD_FAST_TPCclusters90");
			if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_CENT_woSDD_FAST_TPCclusters100");
		}

		//..tracking efficiency
		fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_LHC16t.root");
		if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
			printf("file does not exist");
			return;
		}

		if(fPeriod == "LHC16t_CENT_wSDD"){
			//	default
			hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_wSDD");		
			//	for systematics
			if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_wSDD_FB768"); 
			if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_wSDD_TPCclusters80");
			if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_wSDD_TPCclusters90");
			if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_wSDD_TPCclusters100");
		}
		if(fPeriod == "LHC16t_CENT_woSDD" || fPeriod == "LHC16t_FAST")
		{
			//	default
			hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_woSDD_FAST");		
			//	for systematics
			if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_woSDD_FAST_FB768"); 
			if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_woSDD_FAST_TPCclusters80");
			if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_woSDD_FAST_TPCclusters90");
			if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_DPMJET_CENT_woSDD_FAST_TPCclusters100");
		}
	
	}

	////////////////////////////////////////////////
	if(fPeriod == "LHC16q_CENT_wSDD" || fPeriod == "LHC16q_CENT_woSDD" || fPeriod == "LHC16q_FAST")
	{//	this data set has run-by-run weight. Here I just define files, I get the histograms later in GetWeight(), GetPtWeight()

		if(fIsMC == false)
		{
			//..phi weight
			//	default
			fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_%s.root", fPeriod.Data()));
			//	for systematic checks
			if(fTPCclusters == 80) fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_%s_TPCclusters80.root", fPeriod.Data()));
			if(fTPCclusters == 90) fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_%s_TPCclusters90.root", fPeriod.Data()));
			if(fTPCclusters == 100) fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_%s_TPCclusters100.root", fPeriod.Data()));

			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			//..tracking efficiency
			//	default
			fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_%s.root", fPeriod.Data()));
			//	for systematic checks
			if(fFilterbit == 768) 
				fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_%s_FB768.root", fPeriod.Data()));
			if(fTPCclusters == 80) 
				fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_%s_TPCclusters80.root", fPeriod.Data()));
			if(fTPCclusters == 90) 
				fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_%s_TPCclusters90.root", fPeriod.Data()));
			if(fTPCclusters == 100) 
				fTrackEfficiency = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeV_%s_TPCclusters100.root", fPeriod.Data()));

			if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
				printf("file does not exist");
				return;
			}

		} else{
			//	these are used for MC closure test
			//..phi weight
			fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeightMC_LHC16q.root");
			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			if(fPeriod == "LHC16q_CENT_wSDD") hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_LHC16q_DPMJET_CENT_wSDD_GoodRuns");
			if(fPeriod == "LHC16q_CENT_woSDD") hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_LHC16q_DPMJET_CENT_woSDD_GoodRuns");
			if(fPeriod == "LHC16q_FAST") hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_LHC16q_DPMJET_FAST_GoodRuns");
	
			//..tracking efficiency
			//..NOTE: only LHC16q period has Monte Carlo production
			fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/TrackingEfficiency_pPb5TeVLHC16q.root");
			if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
				printf("file does not exist");
				return;
			}

			if(fPeriod == "LHC16q_CENT_wSDD") hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_3D_FB96_LHC16q_eta08_DPMJET_CENT_wSDD_GoodRuns");
			if(fPeriod == "LHC16q_CENT_woSDD") hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_3D_FB96_LHC16q_eta08_DPMJET_CENT_woSDD_GoodRuns");
			if(fPeriod == "LHC16q_FAST") hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_3D_FB96_LHC16q_eta08_DPMJET_FAST_GoodRuns");

		}	

	}

	//////////////////////////////////////////////
	if(fPeriod == "LHC15oHR")
	{

			//..phi weight: it is done run-by-run, the histograms are obtained in GetWeight() function
			//	default
			fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR.root");
			//	for systematics
			if(fFilterbit == 768) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR_FB768.root");
			if(fTPCclusters == 80) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR_TPCclusters80.root");
			if(fTPCclusters == 90) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR_TPCclusters90.root");
			if(fTPCclusters == 100) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR_TPCclusters100.root");

			//..tracking efficiency: also done run-by-run, the histograms are obtained in GetPtWeight() function
			//	AMPT
			//	default
			fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT.root");
			// for systematics
			if(fFilterbit == 768) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT_FB768.root");
			if(fTPCclusters == 80) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT_TPCclusters80.root");
			if(fTPCclusters == 90) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT_TPCclusters90.root");
			if(fTPCclusters == 100) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT_TPCclusters100.root");

	}

	//////////////////////////////////////////////
	if(fPeriod == "LHC15oLR")
	{

			//..phi weight: it is done run-by-run, the histograms are obtained in GetWeight() function
			//	default
			fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_LIR.root");
			//	for systematics
			if(fFilterbit == 768) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_LIR_FB768.root");
			if(fTPCclusters == 80) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_LIR_TPCclusters80.root");
			if(fTPCclusters == 90) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_LIR_TPCclusters90.root");
			if(fTPCclusters == 100) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_LIR_TPCclusters100.root");

			//..tracking efficiency: also done run-by-run, the histograms are obtained in GetPtWeight() function
			//	default
			fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_LIR.root");
			// for systematics
			if(fFilterbit == 768) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_LIR_FB768.root");
			if(fTPCclusters == 80) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_LIR_TPCclusters80.root");
			if(fTPCclusters == 90) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_LIR_TPCclusters90.root");
			if(fTPCclusters == 100) 
				fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_LIR_TPCclusters100.root");
	
	}

	//////////////////////////////////////////////
	if(fPeriod == "LHC17n")
	{
		fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2017/PhiWeight_XeXe.root");
		if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
			printf("file does not exist");
			return;
		}
		hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight");
		if(fFilterbit == 768) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_FB768");
		if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_TPCclusters80");
		if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_TPCclusters90");
		if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get("fPhiWeight_TPCclusters100");
		
		fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2017/TrackingEfficiency_XeXe.root");
		if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
			printf("file does not exist");
			return;
		}

		hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_HIJING");
		if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_HIJING_FB768");
		if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_HIJING_TPCclusters80");
		if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_HIJING_TPCclusters90");
		if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get("eff_HIJING_TPCclusters100");

	}


	//////////////////////////////////////////////
	if(fPeriod == "LHC15i")
	{

			//..tracking efficiency
			fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_pp13TeV_LHC15i.root");
			if((!fTrackEfficiency) || (!fTrackEfficiency->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_default", fPeriod.Data()));
			//	if any of the below cases are switched on (for systematics), then the histogram is overwritten
			if(fFilterbit == 768) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_FB768", fPeriod.Data()));
			if(fTPCclusters == 80) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hTrackEfficiency = (TH3F*)fTrackEfficiency->Get(Form("eff_%s_PYTHIA_TPCclusters100", fPeriod.Data()));

			//..phi weight
			fPhiWeight = TFile::Open(Form("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_%s.root", fPeriod.Data()));
			if((!fPhiWeight) || (!fPhiWeight->IsOpen())){
				printf("file does not exist");
				return;
			}

			//	default
			hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s", fPeriod.Data())); 
			// for systematic checks
			if(fFilterbit == 768) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_FB768", fPeriod.Data()));
			if(fTPCclusters == 80) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters80", fPeriod.Data()));
			if(fTPCclusters == 90) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters90", fPeriod.Data()));
			if(fTPCclusters == 100) hPhiWeight = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%s_TPCclusters100", fPeriod.Data()));

	}


  // Event Histograms
	hEventCount = new TH1D("hEventCount", "; centrality;;", 1, 0, 1);
	fListOfObjects->Add(hEventCount);

	hMult = new TH1F("hMult", ";number of tracks; entries", nn, 0, range);
	hMult->Sumw2();
	fListOfObjects->Add(hMult);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fVtxAfterCuts->Sumw2();
  fListOfObjects->Add(fVtxAfterCuts);

	if(fTrigger == 0)
	{
  	fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 0.5);
  	fListOfObjects->Add(fCentralityDis);
	}else{
 		fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
 		fListOfObjects->Add(fCentralityDis);
	}   
 
  fV0CentralityDis = new TH1F("fV0CentralityDis", "centrality V0/<V0> distribution; centrality; Counts", 200, 0, 10);
  fListOfObjects->Add(fV0CentralityDis);
  
	hMultV0vsNtrksAfterCuts = new TH2F("hMultV0vsNtrksAfterCuts","V0 mult vs. number of tracks; V0 mult; number of tracks", 100, 0, 10, nn, 0, range);
	fListOfObjects->Add(hMultV0vsNtrksAfterCuts);

	hNtrksVSmultPercentile = new TH2F("hNtrksVSmultPercentile", ";Multiplicity percentile;ITSsa tracks", 100, 0, 100, 1000, 0, 2000);
	fListOfObjects->Add(hNtrksVSmultPercentile);

	//..histograms for PbPb selection
	if(fPeriod == "LHC15oHR" || fPeriod == "LHC15oLR" || fPeriod == "LHC17n")
	{
		fCentralityV0MCL1 = new TH2F("fCentralityV0MCL1", "; centrality V0M; centrality CL1", 100, 0, 100, 100, 0, 100);
		fListOfObjects->Add(fCentralityV0MCL1);

		fCentralityV0MCL0 = new TH2F("fCentralityV0MCL0", "; centrality V0M; centrality CL0", 100, 0, 100, 100, 0, 100);
		fListOfObjects->Add(fCentralityV0MCL0);

		fCentralityCL0CL1 = new TH2F("fCentralityCL0CL1", "; centrality CL0; centrality CL1", 100, 0, 100, 100, 0, 100);
		fListOfObjects->Add(fCentralityCL0CL1);

  	fMultvsCentr = new TH2F("fMultvsCentr", "; centrality V0M; TPC multiplicity (FB32)", 100, 0, 100, 1000, 0, 2000);
		fListOfObjects->Add(fMultvsCentr);

  	fMult128vsCentr = new TH2F("fMult128vsCentr", "; centrality V0M; TPC multiplicity (FB128)", 100, 0, 100, 1000, 0, 2000);
		fListOfObjects->Add(fMult128vsCentr);

		fMultTPCvsTOF = new TH2F("fMultTPCvsTOF", "; TPC FB32 multiplicity; TOF multiplicity", 400, 0, 4000, 200, 0, 2000);
		fListOfObjects->Add(fMultTPCvsTOF);

		fMultTPCvsESD = new TH2F("fMultTPCvsESD", "; TPC FB128 multiplicity; ESD multiplicity", 700, 0, 7000, 1000, -1000, 25000);
		fListOfObjects->Add(fMultTPCvsESD);
	}

	// Track histograms  
  fPhiDis1D = new TH1F("fPhiDis1D", "phi distribution; #phi; Counts", 360, 0, 2*TMath::Pi()); //..18 sectors of TPC, 20 bins per each sector
	fPhiDis1D->Sumw2();
  fListOfObjects->Add(fPhiDis1D);
    
  fPhiDis = new TH3F("fPhiDis", "phi distribution; #phi; Counts", 360, 0, 2*TMath::Pi(), 100, -2, 2, 40, -10, 10);
	fPhiDis->Sumw2();
  fListOfObjects->Add(fPhiDis);
    
  fEtaDis = new TH1F("fEtaDis", "eta distribution; #eta; Counts", 200, -2., 2.);
	fEtaDis->Sumw2();
  fListOfObjects->Add(fEtaDis);
   
  fEtaBefore = new TH1F("fEtaBefore", "eta distribution; #eta; Counts", 200, -2., 2.);
	fEtaBefore->Sumw2();
  fListOfObjects->Add(fEtaBefore);
   
	fPtDis = new TH1F("fPtDis", "pt distribution; p_{T}; Counts", 100, 0, 10);
	fPtDis->Sumw2();
	fListOfObjects->Add(fPtDis);
 
	fPtBefore = new TH1F("fPtBefore", "pt distribution; p_{T}; Counts", 100, 0, 10);
	fPtBefore->Sumw2();
	fListOfObjects->Add(fPtBefore);
 
	hDCAxyBefore = new TH1F("hDCAxyBefore", "; DCAxyBefore", 20, 0, 2); 
	hDCAxyBefore->Sumw2();
	fListOfObjects->Add(hDCAxyBefore);

	hDCAzBefore = new TH1F("hDCAzBefore", "; DCAzBefore", 60, -3, 3); 
	hDCAzBefore->Sumw2();
	fListOfObjects->Add(hDCAzBefore);

	hITSclustersBefore = new TH1F("hITSclustersBefore", "; ITSclustersBefore", 8, 0, 8); 
	hITSclustersBefore->Sumw2();
	fListOfObjects->Add(hITSclustersBefore);

	hChi2Before = new TH1F("hChi2Before", "; Chi2Before", 50, 0, 5); 
	hChi2Before->Sumw2();
	fListOfObjects->Add(hChi2Before);

	hDCAxy = new TH1F("hDCAxy", "; DCAxy", 20, 0, 2); 
	hDCAxy->Sumw2();
	fListOfObjects->Add(hDCAxy);

	hDCAz = new TH1F("hDCAz", "; DCAz", 60, -3, 3); 
	hDCAz->Sumw2();
	fListOfObjects->Add(hDCAz);

	hITSclusters = new TH1F("hITSclusters", "; ITSclusters", 8, 0, 8); 
	hITSclusters->Sumw2();
	fListOfObjects->Add(hITSclusters);

	hChi2 = new TH1F("hChi2", "; Chi2", 50, 0, 5); 
	hChi2->Sumw2();
	fListOfObjects->Add(hChi2);

	// Physics profiles
	// SC(n,m)
	fChsc4242 = new TProfile("fChsc4242", "# of tracks", nn, 0, range);
	fChsc4242->Sumw2();
	fListOfObjects->Add(fChsc4242);

	fChsc4242Gap0 = new TProfile("fChsc4242Gap0", "# of tracks", nn, 0, range);
	fChsc4242Gap0->Sumw2();
	fListOfObjects->Add(fChsc4242Gap0);

	fChsc4242_3sub = new TProfile("fChsc4242_3sub", "# of tracks", nn, 0, range);
	fChsc4242_3sub->Sumw2();
	fListOfObjects->Add(fChsc4242_3sub);

	fChsc4224_3sub = new TProfile("fChsc4224_3sub", "# of tracks", nn, 0, range);
	fChsc4224_3sub->Sumw2();
	fListOfObjects->Add(fChsc4224_3sub);

	fChsc4242_3sub_perm2 = new TProfile("fChsc4242_3sub_perm2", "# of tracks", nn, 0, range);
	fChsc4242_3sub_perm2->Sumw2();
	fListOfObjects->Add(fChsc4242_3sub_perm2);

	fChsc4224_3sub_perm2 = new TProfile("fChsc4224_3sub_perm2", "# of tracks", nn, 0, range);
	fChsc4224_3sub_perm2->Sumw2();
	fListOfObjects->Add(fChsc4224_3sub_perm2);

	fChsc4242_3sub_perm3 = new TProfile("fChsc4242_3sub_perm3", "# of tracks", nn, 0, range);
	fChsc4242_3sub_perm3->Sumw2();
	fListOfObjects->Add(fChsc4242_3sub_perm3);

	fChsc4224_3sub_perm3 = new TProfile("fChsc4224_3sub_perm3", "# of tracks", nn, 0, range);
	fChsc4224_3sub_perm3->Sumw2();
	fListOfObjects->Add(fChsc4224_3sub_perm3);


	fChsc3232 = new TProfile("fChsc3232", "# of tracks", nn, 0, range);
	fChsc3232->Sumw2();
	fListOfObjects->Add(fChsc3232);

	fChsc3232Gap0 = new TProfile("fChsc3232Gap0", "# of tracks", nn, 0, range);
	fChsc3232Gap0->Sumw2();
	fListOfObjects->Add(fChsc3232Gap0);

	fChsc3232_3sub = new TProfile("fChsc3232_3sub", "# of tracks", nn, 0, range);
	fChsc3232_3sub->Sumw2();
	fListOfObjects->Add(fChsc3232_3sub);

	fChsc3223_3sub = new TProfile("fChsc3223_3sub", "# of tracks", nn, 0, range);
	fChsc3223_3sub->Sumw2();
	fListOfObjects->Add(fChsc3223_3sub);

	fChsc3232_3sub_perm2 = new TProfile("fChsc3232_3sub_perm2", "# of tracks", nn, 0, range);
	fChsc3232_3sub_perm2->Sumw2();
	fListOfObjects->Add(fChsc3232_3sub_perm2);

	fChsc3223_3sub_perm2 = new TProfile("fChsc3223_3sub_perm2", "# of tracks", nn, 0, range);
	fChsc3223_3sub_perm2->Sumw2();
	fListOfObjects->Add(fChsc3223_3sub_perm2);

	fChsc3232_3sub_perm3 = new TProfile("fChsc3232_3sub_perm3", "# of tracks", nn, 0, range);
	fChsc3232_3sub_perm3->Sumw2();
	fListOfObjects->Add(fChsc3232_3sub_perm3);

	fChsc3223_3sub_perm3 = new TProfile("fChsc3223_3sub_perm3", "# of tracks", nn, 0, range);
	fChsc3223_3sub_perm3->Sumw2();
	fListOfObjects->Add(fChsc3223_3sub_perm3);


	for(int h=0; h<3; h++)
	{
	  fChcn2Ntrks1bin[h] = new TProfile(Form("fChc%d2Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Ntrks1bin[h]);

    fChcn2Gap0Ntrks1bin[h] = new TProfile(Form("fChc%d2Gap0Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Gap0Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Gap0Ntrks1bin[h]);

    fChcn2Gap2Ntrks1bin[h] = new TProfile(Form("fChc%d2Gap2Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Gap2Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Gap2Ntrks1bin[h]);

    fChcn2Gap4Ntrks1bin[h] = new TProfile(Form("fChc%d2Gap4Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Gap4Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Gap4Ntrks1bin[h]);

    fChcn2Gap8Ntrks1bin[h] = new TProfile(Form("fChc%d2Gap8Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Gap8Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Gap8Ntrks1bin[h]);

    fChcn2GapNtrks1bin[h] = new TProfile(Form("fChc%d2GapNtrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2GapNtrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2GapNtrks1bin[h]);

    fChcn2Gap14Ntrks1bin[h] = new TProfile(Form("fChc%d2Gap14Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    fChcn2Gap14Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2Gap14Ntrks1bin[h]);


    fChcn2_3subLMNtrks1bin[h] = new TProfile(Form("fChc%d2_3subLMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, 0, range);
    fChcn2_3subLMNtrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subLMNtrks1bin[h]);

    fChcn2_3subRMNtrks1bin[h] = new TProfile(Form("fChc%d2_3subRMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, 0, range);
    fChcn2_3subRMNtrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subRMNtrks1bin[h]);

    fChcn2_3subLM_perm2Ntrks1bin[h] = new TProfile(Form("fChc%d2_3subLM_perm2Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, 0, range);
    fChcn2_3subLM_perm2Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subLM_perm2Ntrks1bin[h]);

    fChcn2_3subLR_perm2Ntrks1bin[h] = new TProfile(Form("fChc%d2_3subLR_perm2Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, 0, range);
    fChcn2_3subLR_perm2Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subLR_perm2Ntrks1bin[h]);

    fChcn2_3subRL_perm3Ntrks1bin[h] = new TProfile(Form("fChc%d2_3subRL_perm3Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, 0, range);
    fChcn2_3subRL_perm3Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subRL_perm3Ntrks1bin[h]);

    fChcn2_3subRM_perm3Ntrks1bin[h] = new TProfile(Form("fChc%d2_3subRM_perm3Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, 0, range);
    fChcn2_3subRM_perm3Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn2_3subRM_perm3Ntrks1bin[h]);


	  fChcn4Ntrks1bin[h] = new TProfile(Form("fChc%d4Ntrks1bin", h+2), "<<4>> Re; # of tracks", nn, 0, range);
    fChcn4Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn4Ntrks1bin[h]);

    fChcn4Gap0Ntrks1bin[h] = new TProfile(Form("fChc%d4Gap0Ntrks1bin", h+2), "<<4>> Re; # of tracks", nn, 0, range);
    fChcn4Gap0Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn4Gap0Ntrks1bin[h]);

		fChcn4_3subNtrks1bin[h] = new TProfile(Form("fChc%d4_3subNtrks1bin", h+2), "<<4>> 3-subevent method; # of tracks", nn, 0, range);
		fChcn4_3subNtrks1bin[h]->Sumw2();
		fListOfObjects->Add(fChcn4_3subNtrks1bin[h]);
	
		fChcn4_3sub_perm2Ntrks1bin[h] = new TProfile(Form("fChc%d4_3sub_perm2Ntrks1bin", h+2), "<<4>> for 3-subevent method", nn, 0, range);
		fChcn4_3sub_perm2Ntrks1bin[h]->Sumw2();
		fListOfObjects->Add(fChcn4_3sub_perm2Ntrks1bin[h]);

		fChcn4_3sub_perm3Ntrks1bin[h] = new TProfile(Form("fChc%d4_3sub_perm3Ntrks1bin", h+2), "<<4>> for 3-subevent method", nn, 0, range);
		fChcn4_3sub_perm3Ntrks1bin[h]->Sumw2();
		fListOfObjects->Add(fChcn4_3sub_perm3Ntrks1bin[h]);


	  fChcn6Ntrks1bin[h] = new TProfile(Form("fChc%d6Ntrks1bin", h+2), "<<6>> Re; # of tracks", nn, 0, range);
    fChcn6Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn6Ntrks1bin[h]);

	  fChcn6Gap0Ntrks1bin[h] = new TProfile(Form("fChc%d6Gap0Ntrks1bin", h+2), "<<6>> Gap0 Re; # of tracks", nn, 0, range);
    fChcn6Gap0Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn6Gap0Ntrks1bin[h]);


	  fChcn8Ntrks1bin[h] = new TProfile(Form("fChc%d8Ntrks1bin", h+2), "<<8>>  Re; # of tracks", nn, 0, range);
    fChcn8Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn8Ntrks1bin[h]);

	  fChcn8Gap0Ntrks1bin[h] = new TProfile(Form("fChc%d8Gap0Ntrks1bin", h+2), "<<8>> Gap0 Re; # of tracks", nn, 0, range);
    fChcn8Gap0Ntrks1bin[h]->Sumw2();
    fListOfObjects->Add(fChcn8Gap0Ntrks1bin[h]);

	}// harmonics

	for(int j=0; j<fSample+1; j++)
	{

		fsc4242[j] = new TProfile(Form("fsc4242_number%d", j), "; Ntrks", nn, 0, range);
		fsc4242[j]->Sumw2();
		fListOfObjects->Add(fsc4242[j]);

		fsc4242Gap0[j] = new TProfile(Form("fsc4242Gap0_number%d", j), "; Ntrks", nn, 0, range);
		fsc4242Gap0[j]->Sumw2();
		fListOfObjects->Add(fsc4242Gap0[j]);

		fsc4242_3sub[j] = new TProfile(Form("fsc4242_3sub_number%d", j), "; Ntrks", nn, 0, range);
		fsc4242_3sub[j]->Sumw2();
		fListOfObjects->Add(fsc4242_3sub[j]);

		fsc4224_3sub[j] = new TProfile(Form("fsc4224_3sub_number%d", j), "; Ntrks", nn, 0, range);
		fsc4224_3sub[j]->Sumw2();
		fListOfObjects->Add(fsc4224_3sub[j]);

		fsc4242_3sub_perm2[j] = new TProfile(Form("fsc4242_3sub_perm2_number%d", j), "; Ntrks", nn, 0, range);
		fsc4242_3sub_perm2[j]->Sumw2();
		fListOfObjects->Add(fsc4242_3sub_perm2[j]);

		fsc4224_3sub_perm2[j] = new TProfile(Form("fsc4224_3sub_perm2_number%d", j), "; Ntrks", nn, 0, range);
		fsc4224_3sub_perm2[j]->Sumw2();
		fListOfObjects->Add(fsc4224_3sub_perm2[j]);

		fsc4242_3sub_perm3[j] = new TProfile(Form("fsc4242_3sub_perm3_number%d", j), "; Ntrks", nn, 0, range);
		fsc4242_3sub_perm3[j]->Sumw2();
		fListOfObjects->Add(fsc4242_3sub_perm3[j]);

		fsc4224_3sub_perm3[j] = new TProfile(Form("fsc4224_3sub_perm3_number%d", j), "; Ntrks", nn, 0, range);
		fsc4224_3sub_perm3[j]->Sumw2();
		fListOfObjects->Add(fsc4224_3sub_perm3[j]);


		fsc3232[j] = new TProfile(Form("fsc3232_number%d", j), "; Ntrks", nn, 0, range);
		fsc3232[j]->Sumw2();
		fListOfObjects->Add(fsc3232[j]);

		fsc3232Gap0[j] = new TProfile(Form("fsc3232Gap0_number%d", j), "; Ntrks", nn, 0, range);
		fsc3232Gap0[j]->Sumw2();
		fListOfObjects->Add(fsc3232Gap0[j]);

		fsc3232_3sub[j] = new TProfile(Form("fsc3232_3sub_number%d", j), "; Ntrks", nn, 0, range);
		fsc3232_3sub[j]->Sumw2();
		fListOfObjects->Add(fsc3232_3sub[j]);

		fsc3223_3sub[j] = new TProfile(Form("fsc3223_3sub_number%d", j), "; Ntrks", nn, 0, range);
		fsc3223_3sub[j]->Sumw2();
		fListOfObjects->Add(fsc3223_3sub[j]);

		fsc3232_3sub_perm2[j] = new TProfile(Form("fsc3232_3sub_perm2_number%d", j), "; Ntrks", nn, 0, range);
		fsc3232_3sub_perm2[j]->Sumw2();
		fListOfObjects->Add(fsc3232_3sub_perm2[j]);

		fsc3223_3sub_perm2[j] = new TProfile(Form("fsc3223_3sub_perm2_number%d", j), "; Ntrks", nn, 0, range);
		fsc3223_3sub_perm2[j]->Sumw2();
		fListOfObjects->Add(fsc3223_3sub_perm2[j]);

		fsc3232_3sub_perm3[j] = new TProfile(Form("fsc3232_3sub_perm3_number%d", j), "; Ntrks", nn, 0, range);
		fsc3232_3sub_perm3[j]->Sumw2();
		fListOfObjects->Add(fsc3232_3sub_perm3[j]);

		fsc3223_3sub_perm3[j] = new TProfile(Form("fsc3223_3sub_perm3_number%d", j), "; Ntrks", nn, 0, range);
		fsc3223_3sub_perm3[j]->Sumw2();
		fListOfObjects->Add(fsc3223_3sub_perm3[j]);

		for(int h=0; h<3; h++)
		{
	    fcn2Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dNtrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Ntrks1bin[h][j]);

      fcn2Gap0Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGap0Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Gap0Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Gap0Ntrks1bin[h][j]);

      fcn2Gap2Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGap2Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Gap2Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Gap2Ntrks1bin[h][j]);

      fcn2Gap4Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGap4Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Gap4Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Gap4Ntrks1bin[h][j]);

      fcn2Gap8Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGap8Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Gap8Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Gap8Ntrks1bin[h][j]);

      fcn2GapNtrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGapNtrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2GapNtrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2GapNtrks1bin[h][j]);

      fcn2Gap14Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%dGap14Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      fcn2Gap14Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2Gap14Ntrks1bin[h][j]);


      fcn2_3subLMNtrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subLMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks", nn, 0, range);
      fcn2_3subLMNtrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subLMNtrks1bin[h][j]);

      fcn2_3subRMNtrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subRMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks", nn, 0, range);
      fcn2_3subRMNtrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subRMNtrks1bin[h][j]);

      fcn2_3subLM_perm2Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subLM_perm2Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks", nn, 0, range);
      fcn2_3subLM_perm2Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subLM_perm2Ntrks1bin[h][j]);

      fcn2_3subLR_perm2Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subLR_perm2Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks", nn, 0, range);
      fcn2_3subLR_perm2Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subLR_perm2Ntrks1bin[h][j]);

      fcn2_3subRM_perm3Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subRM_perm3Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks", nn, 0, range);
      fcn2_3subRM_perm3Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subRM_perm3Ntrks1bin[h][j]);

      fcn2_3subRL_perm3Ntrks1bin[h][j] = new TProfile(Form("fc%d2_number%d_3subRL_perm3Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks", nn, 0, range);
      fcn2_3subRL_perm3Ntrks1bin[h][j]->Sumw2();
      fListOfObjects->Add(fcn2_3subRL_perm3Ntrks1bin[h][j]);


			fcn4Ntrks1bin[h][j] = new TProfile(Form("fc%d4_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks", nn, 0, range);
			fcn4Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn4Ntrks1bin[h][j]);
	
			fcn4Gap0Ntrks1bin[h][j] = new TProfile(Form("fc%d4Gap0_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks", nn, 0, range);
			fcn4Gap0Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn4Gap0Ntrks1bin[h][j]);
	
			fcn4_3subNtrks1bin[h][j] = new TProfile(Form("fc%d4_3sub_number%dNtrks1bin", h+2, j), "<<4>> 3-subevent method; # of tracks", nn, 0, range);
			fcn4_3subNtrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn4_3subNtrks1bin[h][j]);

			fcn4_3sub_perm2Ntrks1bin[h][j] = new TProfile(Form("fc%d4_number%d_3sub_perm2Ntrks1bin", h+2, j), "<<4>> for 3-subevent method", nn, 0, range);
			fcn4_3sub_perm2Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn4_3sub_perm2Ntrks1bin[h][j]);

			fcn4_3sub_perm3Ntrks1bin[h][j] = new TProfile(Form("fc%d4_number%d_3sub_perm3Ntrks1bin", h+2, j), "<<4>> for 3-subevent method", nn, 0, range);
			fcn4_3sub_perm3Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn4_3sub_perm3Ntrks1bin[h][j]);


			fcn6Ntrks1bin[h][j] = new TProfile(Form("fc%d6_number%dNtrks1bin", h+2, j), "<<6>>; # of tracks", nn, 0, range);
			fcn6Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn6Ntrks1bin[h][j]);
	
			fcn6Gap0Ntrks1bin[h][j] = new TProfile(Form("fc%d6Gap0_number%dNtrks1bin", h+2, j), "<<6>> Gap0; # of tracks", nn, 0, range);
			fcn6Gap0Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn6Gap0Ntrks1bin[h][j]);
	

			fcn8Ntrks1bin[h][j] = new TProfile(Form("fc%d8_number%dNtrks1bin", h+2, j), "<<8>> ; # of tracks", nn, 0, range);
			fcn8Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn8Ntrks1bin[h][j]);

			fcn8Gap0Ntrks1bin[h][j] = new TProfile(Form("fc%d8Gap0_number%dNtrks1bin", h+2, j), "<<8>> Gap0; # of tracks", nn, 0, range);
			fcn8Gap0Ntrks1bin[h][j]->Sumw2();
			fListOfObjects->Add(fcn8Gap0Ntrks1bin[h][j]);

		}// harmonics

	}// samples

	//........................
	//***MC plots
	if(fIsMC == true)
	{
		hMultMC = new TH1F("hMultMC", ";number of tracks; entries", nn, 0, range);
		hMultMC->Sumw2();
		fListOfObjectsMC->Add(hMultMC);

  	fPhiDisTruth = new TH1F("fPhiDisTruth", "phi distribution; #phi; Counts", 100, 0, 6.29);
		fPhiDisTruth->Sumw2();
  	fListOfObjectsMC->Add(fPhiDisTruth);
  	  
  	fEtaDisTruth = new TH1F("fEtaDisTruth", "eta distribution; #eta; Counts", 200, -2., 2.);
		fEtaDisTruth->Sumw2();
  	fListOfObjectsMC->Add(fEtaDisTruth);
  	 
		fPtDisTruth = new TH1F("fPtDisTruth", "pt distribution; p_{T}; Counts", 100, 0, 10);
		fPtDisTruth->Sumw2();
		fListOfObjectsMC->Add(fPtDisTruth);
 
		//..bin in pT
		double binX[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0};//, 4.0, 5.0
		const int nX = 11;
		//..bin in eta
		double binY[] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		const int nY = 20;
		//..bin in Vz
		double binZ[] = {-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
		const int nZ = 20;

		hReco = new TH3F("hReco", "MC reconstructed; p_{T}; #eta; V_{Z}", nX, binX, nY, binY, nZ, binZ);
		hReco->Sumw2();
		fListOfObjectsMC->Add(hReco);
 
		hTruth = new TH3F("hTruth", "MC truth; p_{T}; #eta; V_{Z}", nX, binX, nY, binY, nZ, binZ);
		hTruth->Sumw2();
		fListOfObjectsMC->Add(hTruth);


		// DCA distributions
		hDCAptMC = new TH2F("hDCAptMC", "pt; DCAz", 50, 0, 5, 100, -1.0, 1.0);
		fListOfObjectsMC->Add(hDCAptMC);

		hDCAptMC_material = new TH2F("hDCAptMC_material", "; pt; DCAz", 50, 0, 5, 100, -1.0, 1.0);
		fListOfObjectsMC->Add(hDCAptMC_material);

		hDCAptMC_weak = new TH2F("hDCAptMC_weak", "; pt ;DCAz", 50, 0, 5, 100, -1.0, 1.0);
		fListOfObjectsMC->Add(hDCAptMC_weak);

		// Correlation histogram of Ntrks
		hNtrksRecoNtrksTruth = new TH2F("hNtrksRecoNtrksTruth", "; Ntrks_reco; Ntrks_truth", range, 0, range, range, 0, range);
		fListOfObjectsMC->Add(hNtrksRecoNtrksTruth);

		// Physics profiles
		// SC(m,n) 
		fChMCsc4242 = new TProfile("fChMCsc4242", "# of tracks", nn, 0, range);
		fChMCsc4242->Sumw2();
		fListOfObjectsMC->Add(fChMCsc4242);

		fChMCsc4242Gap0 = new TProfile("fChMCsc4242Gap0", "# of tracks", nn, 0, range);
		fChMCsc4242Gap0->Sumw2();
		fListOfObjectsMC->Add(fChMCsc4242Gap0);

		fChMCsc4242_3sub = new TProfile("fChMCsc4242_3sub", "# of tracks", nn, 0, range);
		fChMCsc4242_3sub->Sumw2();
		fListOfObjectsMC->Add(fChMCsc4242_3sub);

		fChMCsc4224_3sub = new TProfile("fChMCsc4224_3sub", "# of tracks", nn, 0, range);
		fChMCsc4224_3sub->Sumw2();
		fListOfObjectsMC->Add(fChMCsc4224_3sub);

		fChMCsc3232 = new TProfile("fChMCsc3232", "# of tracks", nn, 0, range);
		fChMCsc3232->Sumw2();
		fListOfObjectsMC->Add(fChMCsc3232);

		fChMCsc3232Gap0 = new TProfile("fChMCsc3232Gap0", "# of tracks", nn, 0, range);
		fChMCsc3232Gap0->Sumw2();
		fListOfObjectsMC->Add(fChMCsc3232Gap0);

		fChMCsc3232_3sub = new TProfile("fChMCsc3232_3sub", "# of tracks", nn, 0, range);
		fChMCsc3232_3sub->Sumw2();
		fListOfObjectsMC->Add(fChMCsc3232_3sub);

		fChMCsc3223_3sub = new TProfile("fChMCsc3223_3sub", "# of tracks", nn, 0, range);
		fChMCsc3223_3sub->Sumw2();
		fListOfObjectsMC->Add(fChMCsc3223_3sub);

		for(int h=0; h<3; h++)
		{
	  	fChMCcn2Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Ntrks1bin[h]);

    	fChMCcn2Gap0Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Gap0Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Gap0Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Gap0Ntrks1bin[h]);

    	fChMCcn2Gap2Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Gap2Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Gap2Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Gap2Ntrks1bin[h]);

    	fChMCcn2Gap4Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Gap4Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Gap4Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Gap4Ntrks1bin[h]);

    	fChMCcn2Gap8Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Gap8Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Gap8Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Gap8Ntrks1bin[h]);

    	fChMCcn2GapNtrks1bin[h] = new TProfile(Form("fChMCc%d2GapNtrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2GapNtrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2GapNtrks1bin[h]);

    	fChMCcn2Gap14Ntrks1bin[h] = new TProfile(Form("fChMCc%d2Gap14Ntrks1bin", h+2), "<<2>> Re; # of tracks", nn, 0, range);
    	fChMCcn2Gap14Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2Gap14Ntrks1bin[h]);


    	fChMCcn2_3subLMNtrks1bin[h] = new TProfile(Form("fChMCc%d2_3subLMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, 0, range);
    	fChMCcn2_3subLMNtrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2_3subLMNtrks1bin[h]);

    	fChMCcn2_3subRMNtrks1bin[h] = new TProfile(Form("fChMCc%d2_3subRMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, 0, range);
    	fChMCcn2_3subRMNtrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn2_3subRMNtrks1bin[h]);


	  	fChMCcn4Ntrks1bin[h] = new TProfile(Form("fChMCc%d4Ntrks1bin", h+2), "<<4>> Re; # of tracks", nn, 0, range);
    	fChMCcn4Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn4Ntrks1bin[h]);

    	fChMCcn4Gap0Ntrks1bin[h] = new TProfile(Form("fChMCc%d4Gap0Ntrks1bin", h+2), "<<4>> Re; # of tracks", nn, 0, range);
    	fChMCcn4Gap0Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn4Gap0Ntrks1bin[h]);

			fChMCcn4_3subNtrks1bin[h] = new TProfile(Form("fChMCc%d4_3subNtrks1bin", h+2), "<<4>> 3-subevent method; # of tracks", nn, 0, range);
			fChMCcn4_3subNtrks1bin[h]->Sumw2();
			fListOfObjectsMC->Add(fChMCcn4_3subNtrks1bin[h]);
	

	  	fChMCcn6Ntrks1bin[h] = new TProfile(Form("fChMCc%d6Ntrks1bin", h+2), "<<6>> Re; # of tracks", nn, 0, range);
    	fChMCcn6Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn6Ntrks1bin[h]);

	  	fChMCcn6Gap0Ntrks1bin[h] = new TProfile(Form("fChMCc%d6Gap0Ntrks1bin", h+2), "<<6>> Gap0 Re; # of tracks", nn, 0, range);
    	fChMCcn6Gap0Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn6Gap0Ntrks1bin[h]);


	  	fChMCcn8Ntrks1bin[h] = new TProfile(Form("fChMCc%d8Ntrks1bin", h+2), "<<8>>  Re; # of tracks", nn, 0, range);
    	fChMCcn8Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn8Ntrks1bin[h]);

	  	fChMCcn8Gap0Ntrks1bin[h] = new TProfile(Form("fChMCc%d8Gap0Ntrks1bin", h+2), "<<8>> Gap0 Re; # of tracks", nn, 0, range);
    	fChMCcn8Gap0Ntrks1bin[h]->Sumw2();
    	fListOfObjectsMC->Add(fChMCcn8Gap0Ntrks1bin[h]);

		}// harmonics

	  for (Int_t j=1; j < fSample+1; j++) 
		{

			fMCsc4242[j] = new TProfile(Form("fMCsc4242_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc4242[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc4242[j]);

			fMCsc4242Gap0[j] = new TProfile(Form("fMCsc4242Gap0_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc4242Gap0[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc4242Gap0[j]);

			fMCsc4242_3sub[j] = new TProfile(Form("fMCsc4242_3sub_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc4242_3sub[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc4242_3sub[j]);

			fMCsc4224_3sub[j] = new TProfile(Form("fMCsc4224_3sub_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc4224_3sub[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc4224_3sub[j]);

			fMCsc3232[j] = new TProfile(Form("fMCsc3232_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc3232[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc3232[j]);

			fMCsc3232Gap0[j] = new TProfile(Form("fMCsc3232Gap0_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc3232Gap0[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc3232Gap0[j]);

			fMCsc3232_3sub[j] = new TProfile(Form("fMCsc3232_3sub_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc3232_3sub[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc3232_3sub[j]);

			fMCsc3223_3sub[j] = new TProfile(Form("fMCsc3223_3sub_number%d", j), "; Ntrks", nn, 0, range);
			fMCsc3223_3sub[j]->Sumw2();
			fListOfObjectsMC->Add(fMCsc3223_3sub[j]);

			for(Int_t h=0; h<3; h++)
			{
	    	fMCcn2Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dNtrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Ntrks1bin[h][j]);

      	fMCcn2Gap0Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGap0Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Gap0Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Gap0Ntrks1bin[h][j]);

      	fMCcn2Gap2Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGap2Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Gap2Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Gap2Ntrks1bin[h][j]);

      	fMCcn2Gap4Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGap4Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Gap4Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Gap4Ntrks1bin[h][j]);

      	fMCcn2Gap8Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGap8Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Gap8Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Gap8Ntrks1bin[h][j]);

      	fMCcn2GapNtrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGapNtrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2GapNtrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2GapNtrks1bin[h][j]);

      	fMCcn2Gap14Ntrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%dGap14Ntrks1bin", h+2, j), "<<2>> Re; # of tracks", nn, 0, range);
      	fMCcn2Gap14Ntrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2Gap14Ntrks1bin[h][j]);


      	fMCcn2_3subLMNtrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%d_3subLMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks", nn, 0, range);
      	fMCcn2_3subLMNtrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2_3subLMNtrks1bin[h][j]);

      	fMCcn2_3subRMNtrks1bin[h][j] = new TProfile(Form("fMCc%d2_number%d_3subRMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks", nn, 0, range);
      	fMCcn2_3subRMNtrks1bin[h][j]->Sumw2();
      	fListOfObjectsMC->Add(fMCcn2_3subRMNtrks1bin[h][j]);


				fMCcn4Ntrks1bin[h][j] = new TProfile(Form("fMCc%d4_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks", nn, 0, range);
				fMCcn4Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn4Ntrks1bin[h][j]);
	
				fMCcn4Gap0Ntrks1bin[h][j] = new TProfile(Form("fMCc%d4Gap0_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks", nn, 0, range);
				fMCcn4Gap0Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn4Gap0Ntrks1bin[h][j]);
	
				fMCcn4_3subNtrks1bin[h][j] = new TProfile(Form("fMCc%d4_3sub_number%dNtrks1bin", h+2, j), "<<4>> 3-subevent method; # of tracks", nn, 0, range);
				fMCcn4_3subNtrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn4_3subNtrks1bin[h][j]);


				fMCcn6Ntrks1bin[h][j] = new TProfile(Form("fMCc%d6_number%dNtrks1bin", h+2, j), "<<6>>; # of tracks", nn, 0, range);
				fMCcn6Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn6Ntrks1bin[h][j]);
	
				fMCcn6Gap0Ntrks1bin[h][j] = new TProfile(Form("fMCc%d6Gap0_number%dNtrks1bin", h+2, j), "<<6>> Gap0; # of tracks", nn, 0, range);
				fMCcn6Gap0Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn6Gap0Ntrks1bin[h][j]);
	

				fMCcn8Ntrks1bin[h][j] = new TProfile(Form("fMCc%d8_number%dNtrks1bin", h+2, j), "<<8>> ; # of tracks", nn, 0, range);
				fMCcn8Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn8Ntrks1bin[h][j]);

				fMCcn8Gap0Ntrks1bin[h][j] = new TProfile(Form("fMCc%d8Gap0_number%dNtrks1bin", h+2, j), "<<8>> Gap0; # of tracks", nn, 0, range);
				fMCcn8Gap0Ntrks1bin[h][j]->Sumw2();
				fListOfObjectsMC->Add(fMCcn8Gap0Ntrks1bin[h][j]);

			}// harmonics

		}// random sample

	}// MC

  // Post output data.
  PostData(1, fListOfObjects);
	if(fIsMC == true) PostData(2, fListOfObjectsMC);
}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlowGF::UserExec(Option_t *) 
{

	//..apply physics selection
	UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	Bool_t isTrigselected = false;
	if(fTrigger == 0) isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
	if(fTrigger == 1) isTrigselected = fSelectMask&AliVEvent::kINT7;
	if(isTrigselected == false) return;

	//..check if I have AOD
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  
	hEventCount->Fill("after PS", 1.);

	//..do the event selection
	if(!fEventCuts.AcceptEvent(fAOD))
	{// automatic event selection for Run2
		PostData(1,fListOfObjects);
		if(fIsMC == true) PostData(2,fListOfObjectsMC);
		return;
	}

	//..filling Vz distribution
	AliVVertex *vtx = fAOD->GetPrimaryVertex();
	float fVtxZ = vtx->GetZ();	  
	if(TMath::Abs(fVtxZ) > fVtxCut) return;
	fVtxAfterCuts->Fill(fVtxZ);

	//..standard event plots (cent. percentiles, mult-vs-percentile)
	float fMultV0Meq = 0;
	float fMultMeanV0M = 0;
	float centrV0 = 0;
	float cent = 0;
	float v0Centr = 0;
	float cl1Centr = 0;
	float cl0Centr = 0;

	if(fPeriod != "LHC17m" && fPeriod != "LHC17o")
	{
		AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
		MultSelection->GetEstimator("V0M")->SetUseAnchor(kTRUE); // required to have same turn-on of HM trigger after I cut on 0-0.1% events
		fMultV0Meq = MultSelection->GetEstimator("V0M")->GetValue();
		fMultMeanV0M = MultSelection->GetEstimator("V0M")->GetMean();
		centrV0 = MultSelection->GetMultiplicityPercentile("V0M");
		if(fPeriod == "LHC15oHR" || fPeriod == "LHC15oLR" || fPeriod == "LHC17n")
		{
			v0Centr = centrV0;
			cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
			cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
		}

		cent = fMultV0Meq/fMultMeanV0M;
	}

	if(fPeriod == "LHC15oHR" || fPeriod == "LHC15oLR" || fPeriod == "LHC17n")
	{

		fCentralityV0MCL1->Fill(v0Centr, cl1Centr);
		fCentralityV0MCL0->Fill(v0Centr, cl0Centr);
		fCentralityCL0CL1->Fill(cl0Centr, cl1Centr);

		int multTPC = 0;
		int multTPC32 = 0;
		int multTOF = 0;
		int multESD = 0;
		int multTrk = 0;
		const int nTracks = fAOD->GetNumberOfTracks();
		for(int it=0; it<nTracks; it++)
		{
			AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(it);
			if(!track)
			{
				delete track;
				continue;
			}

			if(track->TestFilterBit(128)) multTPC++;

			if(track->TestFilterBit(32))	{
				multTPC32++;
				if ( TMath::Abs(track->GetTOFsignalDz()) <= 10. && track->GetTOFsignal() >= 12000. && track->GetTOFsignal() <= 25000.) multTOF++;
				if((TMath::Abs(track->Eta())) < fEtaCut && (track->GetTPCNcls() >= fTPCclusters) && (track->Pt() >= fMinPt) && (track->Pt() < fMaxPt)) multTrk++;
			}
		}
		multESD = ((AliAODHeader*)fInputEvent->GetHeader())->GetNumberOfESDTracks();
		//..do the Alex event selection
		float multESDTPCdif = 0;
		if(fPeriod == "LHC17n") multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC + 0.000126397*multTPC*multTPC);
		else multESDTPCdif = multESD - multTPC*3.38;

		if(fPeriod == "LHC17n") {
			if(multESDTPCdif > 1000) return;	
			if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return;
		}
		else {
			if(multESDTPCdif > 15000) return;
			if(float(multTOF) < fMultTOFLowCut->Eval(float(multTPC32))) return;
			if(float(multTOF) > fMultTOFHighCut->Eval(float(multTPC32))) return;
			if(float(multTrk) < fMultCentLowCut->Eval(v0Centr)) return;
			if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return;
		}
		//..after Alex event selection

		//..fill QA plots
		fCentralityV0MCL1->Fill(v0Centr, cl1Centr);
		fCentralityV0MCL0->Fill(v0Centr, cl0Centr);
		fCentralityCL0CL1->Fill(cl0Centr, cl1Centr);

		fMult128vsCentr->Fill(v0Centr, multTPC);
		fMultvsCentr->Fill(v0Centr, multTrk);
		fMultTPCvsTOF->Fill(multTPC32, multTOF);
		fMultTPCvsESD->Fill(multTPC, multESD);

	}

	fCentralityDis->Fill(centrV0);
	fV0CentralityDis->Fill(cent);

	AnalyzeAOD(fInputEvent, centrV0, cent, fVtxZ);

   // Post output data.
		PostData(1, fListOfObjects);
		if(fIsMC == true) PostData(2, fListOfObjectsMC);
}

//________________________________________________________________________
void AliAnalysisTaskChargedFlowGF::AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float fVtxZ)
{  

	const int nAODTracks = aod->GetNumberOfTracks();
	
	double Ntrks = 0; //..count number of tracks after track cuts (and possible primary, pions, etc.)
	double NtrksAfter = 0;
	double NtrksAfterGap0M = 0;
	double NtrksAfterGap0P = 0;
	double NtrksAfterGap2M = 0;
	double NtrksAfterGap2P = 0;
	double NtrksAfterGap4M = 0;
	double NtrksAfterGap4P = 0;
	double NtrksAfterGap8M = 0;
	double NtrksAfterGap8P = 0;
	double NtrksAfterGapM = 0;
	double NtrksAfterGapP = 0;
	double NtrksAfterGap14M = 0;
	double NtrksAfterGap14P = 0;
	double NtrksAfter3subL = 0;
	double NtrksAfter3subM = 0;
	double NtrksAfter3subR = 0;

	//..for DCA	
	double pos[3], vz, vx, vy;
	vz = aod->GetPrimaryVertex()->GetZ();
	vx = aod->GetPrimaryVertex()->GetX();
	vy = aod->GetPrimaryVertex()->GetY();

	//..loop to get number of tracks in an event
	//........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++)
  {
  
    AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);
			
    if (!aodTrk){
      delete aodTrk;
      continue;
   	}

		aodTrk->GetXYZ(pos);
		double dcaX = 100;
		double dcaY = 100;
		double dcaZ = 100;
		dcaX = pos[0] - vx;
		dcaY = pos[1] - vy;
		dcaZ = pos[2] - vz;
		double dcaXY = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);

		int nClustersITS = 0;
		nClustersITS = aodTrk->GetITSNcls();
		float chi2PerClusterITS = -1;
		chi2PerClusterITS = aodTrk->GetITSchi2()/float(nClustersITS);

		//..cut on filter bit
		if(!(aodTrk->TestFilterBit(fFilterbit))) continue;

		if(fUseDCAzCut == true) 
			if(TMath::Abs(dcaZ) > fDCAz) continue;

		if(fUseDCAxyCut == true) 
			if(TMath::Abs(dcaXY) > fDCAxy) continue;

		if(aodTrk->GetTPCNcls() < fTPCclusters) continue;

    if (TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

		if(aodTrk->Pt() < fMinPt) continue;
		if(aodTrk->Pt() > fMaxPt) continue;
		
		Ntrks += 1;

	}

	double Ntruth = 0;
	if(fIsMC == true)
	{
		Ntruth = ProcessMCTruth(aod, fVtxZ);
	
		//..fill NTrks_truth vs. Ntrks_reco for the Xaxis correction when doing MCclosure test (to avoid multiplicity fluctuations)
		hNtrksRecoNtrksTruth->Fill(Ntrks, Ntruth);
	}

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGapM[20][20] = {0};
	double QsinGapM[20][20] = {0};
	double QcosGapP[20][20] = {0};
	double QsinGapP[20][20] = {0};

	double QcosGap0M[20][20] = {0};
	double QsinGap0M[20][20] = {0};
	double QcosGap0P[20][20] = {0};
	double QsinGap0P[20][20] = {0};
  
	double QcosGap2M[20][20] = {0};
	double QsinGap2M[20][20] = {0};
	double QcosGap2P[20][20] = {0};
	double QsinGap2P[20][20] = {0};
  
	double QcosGap4M[20][20] = {0};
	double QsinGap4M[20][20] = {0};
	double QcosGap4P[20][20] = {0};
	double QsinGap4P[20][20] = {0};
	
	double QcosGap8M[20][20] = {0};
	double QsinGap8M[20][20] = {0};
	double QcosGap8P[20][20] = {0};
	double QsinGap8P[20][20] = {0};
    
	double QcosGap14M[20][20] = {0};
	double QsinGap14M[20][20] = {0};
	double QcosGap14P[20][20] = {0};
	double QsinGap14P[20][20] = {0};
    
	double QcosSubLeft[20][20] = {0};
	double QsinSubLeft[20][20] = {0};
	double QcosSubMiddle[20][20] = {0};
	double QsinSubMiddle[20][20] = {0};
	double QcosSubRight[20][20] = {0};
	double QsinSubRight[20][20] = {0};
    
  Int_t fBin = 0;
  TRandom3 rr(0);
  Double_t ranNum = 0.;

		if(fSample == 1) fBin = 1;
    
    if(fSample == 10)
		{
    	ranNum = rr.Rndm();
        
    	if (ranNum <= 0.1 && ranNum > 0. ) fBin = 1;
    	else if (ranNum <= 0.2 && ranNum > 0.1) fBin = 2;
    	else if (ranNum <= 0.3 && ranNum > 0.2) fBin = 3;
    	else if (ranNum <= 0.4 && ranNum > 0.3) fBin = 4;
    	else if (ranNum <= 0.5 && ranNum > 0.4) fBin = 5;
    	else if (ranNum <= 0.6 && ranNum > 0.5) fBin = 6;
    	else if (ranNum <= 0.7 && ranNum > 0.6) fBin = 7;
    	else if (ranNum <= 0.8 && ranNum > 0.7) fBin = 8;
    	else if (ranNum <= 0.9 && ranNum > 0.8) fBin = 9;
    	else if (ranNum > 0.9 && ranNum <= 1.) fBin = 10;
    	else {
         cout << ranNum << endl;
         
         cout << "This is very strange!!!" << endl;
         fBin = 11;
    	}

    }

	double runNumber = fInputEvent->GetRunNumber();

	//..LOOP OVER TRACKS........	
	//........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++)
  {
 
    AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);
			
    if (!aodTrk){
      delete aodTrk;
      continue;
   	}

		aodTrk->GetXYZ(pos);
		double dcaZ = 100;
		double dcaX = 100;
		double dcaY = 100;
		double dcaXY = 100;
		dcaZ = pos[2] - vz;
		dcaX = pos[0] - vx;
		dcaY = pos[1] - vy;
		dcaXY = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);

		int nClustersITS = 0;
		nClustersITS = aodTrk->GetITSNcls();
		float chi2PerClusterITS = -1;
		if(nClustersITS != 0) chi2PerClusterITS = aodTrk->GetITSchi2()/float(nClustersITS);

		hDCAxyBefore->Fill(dcaXY);
		hDCAzBefore->Fill(dcaZ);
		hITSclustersBefore->Fill(nClustersITS);
		hChi2Before->Fill(chi2PerClusterITS);

		//..cut on filter bit
		if(!(aodTrk->TestFilterBit(fFilterbit))) continue;

		fEtaBefore->Fill(aodTrk->Eta());
		fPtBefore->Fill(aodTrk->Pt());

		if(aodTrk->Pt() < fMinPt) continue;
		if(aodTrk->Pt() > fMaxPt) continue;

    if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

		if(fUseDCAzCut == true) 
			if(TMath::Abs(dcaZ) > fDCAz) continue;

		if(fUseDCAxyCut == true) 
			if(TMath::Abs(dcaXY) > fDCAxy) continue;

		if(aodTrk->GetTPCNcls() < fTPCclusters) continue;
	
		hDCAxy->Fill(dcaXY);
		hDCAz->Fill(dcaZ);
		hITSclusters->Fill(nClustersITS);
		hChi2->Fill(chi2PerClusterITS);

		fPhiDis1D->Fill(aodTrk->Phi());
    fPhiDis->Fill(aodTrk->Phi(), aodTrk->Eta(), fVtxZ);
    fEtaDis->Fill(aodTrk->Eta());
		fPtDis->Fill(aodTrk->Pt());     
 
		if(fIsMC == true)
		{
			hReco->Fill(aodTrk->Pt(), aodTrk->Eta(), fVtxZ);
		}

		NtrksAfter += 1;

		//..get phi-weight for NUA correction
		double weight = 1;
		if(fNUA == 1) {
			weight = GetWeight(aodTrk->Phi(), aodTrk->Eta(), fVtxZ, runNumber);
		}
		double weightPt = 1;
		if(fNUE == 1) weightPt = GetPtWeight(aodTrk->Pt(), aodTrk->Eta(), fVtxZ, runNumber);



			//..calculate Q-vectors
			//..no eta gap
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}

      //..Gap > 0.0
      if (aodTrk->Eta() > 0.)
      {
				NtrksAfterGap0P += 1;
         
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
 
      }
      if (aodTrk->Eta() < -0.)
      {
				NtrksAfterGap0M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
		
      //..Gap > 0.2
      if (aodTrk->Eta() > 0.1)
      {
        NtrksAfterGap2P += 1; 
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
 
      }
      if (aodTrk->Eta() < -0.1)
      {
				NtrksAfterGap2M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
		
      //..Gap > 0.4
      if (aodTrk->Eta() > 0.2)
      {
				NtrksAfterGap4P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
      if (aodTrk->Eta() < -0.2)
      {
				NtrksAfterGap4M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
	 	
      //..Gap > 0.8
      if (aodTrk->Eta() > 0.4)
      {
				NtrksAfterGap8P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
      if (aodTrk->Eta() < -0.4)
      {
				NtrksAfterGap8M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}

      }

			//..Gap > 1.0
			if(aodTrk->Eta() < -0.5)
			{
				NtrksAfterGapM += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapM[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGapM[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
			} 
			if(aodTrk->Eta() > 0.5)
			{
				NtrksAfterGapP += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapP[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGapP[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());	
					}
				}
			}

      //..Gap > 1.4
      if (aodTrk->Eta() > 0.7)
      {
				NtrksAfterGap14P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
          
      if (aodTrk->Eta() < -0.7)
      {
				NtrksAfterGap14M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}

      }

			//..3-subevent method
			if(aodTrk->Eta() < -0.4)
			{//..left part
				NtrksAfter3subL += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
			}
			if(aodTrk->Eta() >= -0.4 && aodTrk->Eta() <= 0.4)
			{//..middle part
				NtrksAfter3subM += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
			}
			if(aodTrk->Eta() > 0.4)
			{//..right part
				NtrksAfter3subR += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
			}


  } // end loop of all track

	hMult->Fill(NtrksAfter);

	hMultV0vsNtrksAfterCuts->Fill(cent, NtrksAfter); //..cent = M/<M>
	hNtrksVSmultPercentile->Fill(centrV0, NtrksAfter);

	//............................
	//..GENERIC FRAMEWORK RP
	//............................
 
	//..calculate Q-vector for each harmonics n and power p
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);
			QvectorM[iharm][ipow] = TComplex(QcosGapM[iharm][ipow], QsinGapM[iharm][ipow]);
			QvectorP[iharm][ipow] = TComplex(QcosGapP[iharm][ipow], QsinGapP[iharm][ipow]);
			Qvector0M[iharm][ipow] = TComplex(QcosGap0M[iharm][ipow], QsinGap0M[iharm][ipow]);
			Qvector0P[iharm][ipow] = TComplex(QcosGap0P[iharm][ipow], QsinGap0P[iharm][ipow]);
			Qvector2M[iharm][ipow] = TComplex(QcosGap2M[iharm][ipow], QsinGap2M[iharm][ipow]);
			Qvector2P[iharm][ipow] = TComplex(QcosGap2P[iharm][ipow], QsinGap2P[iharm][ipow]);
			Qvector4M[iharm][ipow] = TComplex(QcosGap4M[iharm][ipow], QsinGap4M[iharm][ipow]);
			Qvector4P[iharm][ipow] = TComplex(QcosGap4P[iharm][ipow], QsinGap4P[iharm][ipow]);
			Qvector8M[iharm][ipow] = TComplex(QcosGap8M[iharm][ipow], QsinGap8M[iharm][ipow]);
			Qvector8P[iharm][ipow] = TComplex(QcosGap8P[iharm][ipow], QsinGap8P[iharm][ipow]);
			Qvector14M[iharm][ipow] = TComplex(QcosGap14M[iharm][ipow], QsinGap14M[iharm][ipow]);
			Qvector14P[iharm][ipow] = TComplex(QcosGap14P[iharm][ipow], QsinGap14P[iharm][ipow]);
			QvectorSubLeft[iharm][ipow] = TComplex(QcosSubLeft[iharm][ipow], QsinSubLeft[iharm][ipow]);
			QvectorSubRight[iharm][ipow] = TComplex(QcosSubRight[iharm][ipow], QsinSubRight[iharm][ipow]);
			QvectorSubMiddle[iharm][ipow] = TComplex(QcosSubMiddle[iharm][ipow], QsinSubMiddle[iharm][ipow]);
		}
	}
	
	//..calculate 2-particle correlations
	//..................................
	double Dn2 = Two(0, 0).Re();
	double Dn2Gap = TwoGap(0, 0).Re();
	double Dn2Gap0 = TwoGap0(0, 0).Re();
	double Dn2Gap2 = TwoGap2(0, 0).Re();
	double Dn2Gap4 = TwoGap4(0, 0).Re();
	double Dn2Gap8 = TwoGap8(0, 0).Re();
	double Dn2Gap14 = TwoGap14(0, 0).Re();
	double Dn2_3subLM = Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = Two_3SubRM(0, 0).Re();
	double Dn2_3subLM_perm2 = Two_3SubLM_perm2(0, 0).Re();
	double Dn2_3subLR_perm2 = Two_3SubLR_perm2(0, 0).Re();
	double Dn2_3subRM_perm3 = Two_3SubRM_perm3(0, 0).Re();
	double Dn2_3subRL_perm3 = Two_3SubRL_perm3(0, 0).Re();

	if(NtrksAfter > 1 && Dn2 != 0)
	{	

		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		fChcn2Ntrks1bin[0]->Fill(NtrksAfter, v22Re, Dn2);
		fcn2Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		fChcn2Ntrks1bin[1]->Fill(NtrksAfter, v32Re, Dn2);
		fcn2Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		fChcn2Ntrks1bin[2]->Fill(NtrksAfter, v42Re, Dn2);
		fcn2Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42Re, Dn2);

	}	

	if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 0 && Dn2Gap0 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap0 = TwoGap0(2, -2);
		double v22ReGap0 = v22Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[0]->Fill(NtrksAfter, v22ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap0, Dn2Gap0);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap0 = TwoGap0(3, -3);
		double v32ReGap0 = v32Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[1]->Fill(NtrksAfter, v32ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap0, Dn2Gap0);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap0 = TwoGap0(4, -4);
		double v42ReGap0 = v42Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[2]->Fill(NtrksAfter, v42ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap0, Dn2Gap0);

	}

	if(NtrksAfterGap2M > 0 && NtrksAfterGap2P > 0 && Dn2Gap2 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap2 = TwoGap2(2, -2);
		double v22ReGap2 = v22Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[0]->Fill(NtrksAfter, v22ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap2, Dn2Gap2);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap2 = TwoGap2(3, -3);
		double v32ReGap2 = v32Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[1]->Fill(NtrksAfter, v32ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap2, Dn2Gap2);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap2 = TwoGap2(4, -4);
		double v42ReGap2 = v42Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[2]->Fill(NtrksAfter, v42ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap2, Dn2Gap2);

	}

	if(NtrksAfterGap4M > 0 && NtrksAfterGap4P > 0 && Dn2Gap4 != 0)
	{
		//..v2{2} with eta Gap > 0.4
		TComplex v22Gap4 = TwoGap4(2, -2);
		double v22ReGap4 = v22Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[0]->Fill(NtrksAfter, v22ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap4, Dn2Gap4);

		//..v3{2} with eta Gap > 0.4
		TComplex v32Gap4 = TwoGap4(3, -3);
		double v32ReGap4 = v32Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[1]->Fill(NtrksAfter, v32ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap4, Dn2Gap4);

		//..v4{2} with eta Gap > 0.4
		TComplex v42Gap4 = TwoGap4(4, -4);
		double v42ReGap4 = v42Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[2]->Fill(NtrksAfter, v42ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap4, Dn2Gap4);

	}

	if(NtrksAfterGap8M > 0 && NtrksAfterGap8P > 0 && Dn2Gap8 != 0)
	{
		//..v2{2} with eta Gap > 0.8
		TComplex v22Gap8 = TwoGap8(2, -2);
		double v22ReGap8 = v22Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[0]->Fill(NtrksAfter, v22ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap8, Dn2Gap8);

		//..v3{2} with eta Gap > 0.8
		TComplex v32Gap8 = TwoGap8(3, -3);
		double v32ReGap8 = v32Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[1]->Fill(NtrksAfter, v32ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap8, Dn2Gap8);

		//..v4{2} with eta Gap > 0.8
		TComplex v42Gap8 = TwoGap8(4, -4);
		double v42ReGap8 = v42Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[2]->Fill(NtrksAfter, v42ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap8, Dn2Gap8);

	}

	if(NtrksAfterGapM > 0 && NtrksAfterGapP > 0 && Dn2Gap != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap = TwoGap(2, -2);
		double v22ReGap = v22Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[0]->Fill(NtrksAfter, v22ReGap, Dn2Gap);
		fcn2GapNtrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap, Dn2Gap);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap = TwoGap(3, -3);
		double v32ReGap = v32Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[1]->Fill(NtrksAfter, v32ReGap, Dn2Gap);
		fcn2GapNtrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap, Dn2Gap);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap = TwoGap(4, -4);
		double v42ReGap = v42Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[2]->Fill(NtrksAfter, v42ReGap, Dn2Gap);
		fcn2GapNtrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap, Dn2Gap);

	}

	if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[0]->Fill(NtrksAfter, v22ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22ReGap14, Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[1]->Fill(NtrksAfter, v32ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32ReGap14, Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[2]->Fill(NtrksAfter, v42ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42ReGap14, Dn2Gap14);

	}


	//..for 3-subevent method, Gap0
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[0]->Fill(NtrksAfter, v22Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subLM, Dn2_3subLM);

		TComplex v32_3subLM = Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[1]->Fill(NtrksAfter, v32Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subLM, Dn2_3subLM);

		TComplex v42_3subLM = Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[2]->Fill(NtrksAfter, v42Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subLM, Dn2_3subLM);

	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[0]->Fill(NtrksAfter, v22Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subRM, Dn2_3subRM);

		TComplex v32_3subRM = Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[1]->Fill(NtrksAfter, v32Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subRM, Dn2_3subRM);

		TComplex v42_3subRM = Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[2]->Fill(NtrksAfter, v42Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subRM, Dn2_3subRM);

	}

	//..for 3-subevent method, perm 2
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM_perm2 != 0)
	{//..left+middle
		TComplex v22_3subLM_perm2 = Two_3SubLM_perm2(2, -2);
		double v22Re_3subLM_perm2 = v22_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[0]->Fill(NtrksAfter, v22Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subLM_perm2, Dn2_3subLM_perm2);

		TComplex v32_3subLM_perm2 = Two_3SubLM_perm2(3, -3);
		double v32Re_3subLM_perm2 = v32_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[1]->Fill(NtrksAfter, v32Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subLM_perm2, Dn2_3subLM_perm2);

		TComplex v42_3subLM_perm2 = Two_3SubLM_perm2(4, -4);
		double v42Re_3subLM_perm2 = v42_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[2]->Fill(NtrksAfter, v42Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subLM_perm2, Dn2_3subLM_perm2);

	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR_perm2 != 0)
	{//..left+right
		TComplex v22_3subLR_perm2 = Two_3SubLR_perm2(2, -2);
		double v22Re_3subLR_perm2 = v22_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[0]->Fill(NtrksAfter, v22Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subLR_perm2, Dn2_3subLR_perm2);

		TComplex v32_3subLR_perm2 = Two_3SubLR_perm2(3, -3);
		double v32Re_3subLR_perm2 = v32_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[1]->Fill(NtrksAfter, v32Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subLR_perm2, Dn2_3subLR_perm2);

		TComplex v42_3subLR_perm2 = Two_3SubLR_perm2(4, -4);
		double v42Re_3subLR_perm2 = v42_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[2]->Fill(NtrksAfter, v42Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subLR_perm2, Dn2_3subLR_perm2);

	}

	//..for 3-subevent method, perm 3
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subRL_perm3 != 0)
	{//..left+right
		TComplex v22_3subRL_perm3 = Two_3SubRL_perm3(2, -2);
		double v22Re_3subRL_perm3 = v22_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[0]->Fill(NtrksAfter, v22Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subRL_perm3, Dn2_3subRL_perm3);

		TComplex v32_3subRL_perm3 = Two_3SubRL_perm3(3, -3);
		double v32Re_3subRL_perm3 = v32_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[1]->Fill(NtrksAfter, v32Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subRL_perm3, Dn2_3subRL_perm3);

		TComplex v42_3subRL_perm3 = Two_3SubRL_perm3(4, -4);
		double v42Re_3subRL_perm3 = v42_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[2]->Fill(NtrksAfter, v42Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subRL_perm3, Dn2_3subRL_perm3);

	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM_perm3 != 0)
	{//..middle+right
		TComplex v22_3subRM_perm3 = Two_3SubRM_perm3(2, -2);
		double v22Re_3subRM_perm3 = v22_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[0]->Fill(NtrksAfter, v22Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, v22Re_3subRM_perm3, Dn2_3subRM_perm3);

		TComplex v32_3subRM_perm3 = Two_3SubRM_perm3(3, -3);
		double v32Re_3subRM_perm3 = v32_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[1]->Fill(NtrksAfter, v32Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, v32Re_3subRM_perm3, Dn2_3subRM_perm3);

		TComplex v42_3subRM_perm3 = Two_3SubRM_perm3(4, -4);
		double v42Re_3subRM_perm3 = v42_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[2]->Fill(NtrksAfter, v42Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, v42Re_3subRM_perm3, Dn2_3subRM_perm3);

	}


	//..calculate 4-particle correlations
	//................................
	double Dn4 = Four(0, 0, 0, 0).Re();
	double Dn4Gap0 = FourGap0(0, 0, 0, 0).Re();
	double Dn4_3sub = Four_3SubEvts(0, 0, 0, 0).Re();
	double Dn4_3sub_perm2 = Four_3SubEvts_perm2(0, 0, 0, 0).Re();
	double Dn4_3sub_perm3 = Four_3SubEvts_perm3(0, 0, 0, 0).Re();

	if(NtrksAfter > 3 && Dn4 != 0)
	{
		
		TComplex v24 = Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		fChcn4Ntrks1bin[0]->Fill(NtrksAfter, v24Re, Dn4);
		fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Re, Dn4);

		TComplex v34 = Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		fChcn4Ntrks1bin[1]->Fill(NtrksAfter, v34Re, Dn4);
		fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Re, Dn4);

		TComplex v44 = Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		fChcn4Ntrks1bin[2]->Fill(NtrksAfter, v44Re, Dn4);
		fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Re, Dn4);

	}

	if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 !=0)
	{

		TComplex v24Gap0 = FourGap0(2, 2, -2, -2);
		double v24Gap0Re = v24Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[0]->Fill(NtrksAfter, v24Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Gap0Re, Dn4Gap0);

		TComplex v34Gap0 = FourGap0(3, 3, -3, -3);
		double v34Gap0Re = v34Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[1]->Fill(NtrksAfter, v34Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Gap0Re, Dn4Gap0);

		TComplex v44Gap0 = FourGap0(4, 4, -4, -4);
		double v44Gap0Re = v44Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[2]->Fill(NtrksAfter, v44Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Gap0Re, Dn4Gap0);

	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 1 && NtrksAfter3subR > 0 && Dn4_3sub != 0)
	{
		TComplex v24_3sub = Four_3SubEvts(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[0]->Fill(NtrksAfter, v24_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[0][fBin]->Fill(NtrksAfter, v24_3subRe, Dn4_3sub);

		TComplex v34_3sub = Four_3SubEvts(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[1]->Fill(NtrksAfter, v34_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[1][fBin]->Fill(NtrksAfter, v34_3subRe, Dn4_3sub);

		TComplex v44_3sub = Four_3SubEvts(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[2]->Fill(NtrksAfter, v44_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[2][fBin]->Fill(NtrksAfter, v44_3subRe, Dn4_3sub);
	}

	//	c24_3sub perm 2
	if(NtrksAfter3subL > 1 && NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn4_3sub_perm2 != 0)
	{
		TComplex v24_3sub_perm2 = Four_3SubEvts_perm2(2, 2, -2, -2);
		double v24_3sub_perm2Re = v24_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[0]->Fill(NtrksAfter, v24_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24_3sub_perm2Re, Dn4_3sub_perm2);

		TComplex v34_3sub_perm2 = Four_3SubEvts_perm2(3, 3, -3, -3);
		double v34_3sub_perm2Re = v34_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[1]->Fill(NtrksAfter, v34_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34_3sub_perm2Re, Dn4_3sub_perm2);

		TComplex v44_3sub_perm2 = Four_3SubEvts_perm2(4, 4, -4, -4);
		double v44_3sub_perm2Re = v44_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[2]->Fill(NtrksAfter, v44_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44_3sub_perm2Re, Dn4_3sub_perm2);
	}

	//	c24_3sub perm 2
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && NtrksAfter3subR > 1 && Dn4_3sub_perm3 != 0)
	{
		TComplex v24_3sub_perm3 = Four_3SubEvts_perm3(2, 2, -2, -2);
		double v24_3sub_perm3Re = v24_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[0]->Fill(NtrksAfter, v24_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24_3sub_perm3Re, Dn4_3sub_perm3);

		TComplex v34_3sub_perm3 = Four_3SubEvts_perm3(3, 3, -3, -3);
		double v34_3sub_perm3Re = v34_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[1]->Fill(NtrksAfter, v34_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34_3sub_perm3Re, Dn4_3sub_perm3);

		TComplex v44_3sub_perm3 = Four_3SubEvts_perm3(4, 4, -4, -4);
		double v44_3sub_perm3Re = v44_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[2]->Fill(NtrksAfter, v44_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44_3sub_perm3Re, Dn4_3sub_perm3);
	}


	//..SC(n,m)
	//................................
	if(NtrksAfter > 3 && Dn4 != 0)
	{
		//..SC(4,2,-4,-2)	
		TComplex sc4242 = Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		fChsc4242->Fill(NtrksAfter, sc4242Re, Dn4);
		fsc4242[fBin]->Fill(NtrksAfter, sc4242Re, Dn4);

		//..SC(3,2,-3,-2)	
		TComplex sc3232 = Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		fChsc3232->Fill(NtrksAfter, sc3232Re, Dn4);
		fsc3232[fBin]->Fill(NtrksAfter, sc3232Re, Dn4);

	}	

	if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 != 0)
	{

		TComplex sc4242Gap0 = FourGap0(4, 2, -4, -2);
		double sc4242Gap0Re = sc4242Gap0.Re()/Dn4Gap0;
		fChsc4242Gap0->Fill(NtrksAfter, sc4242Gap0Re, Dn4Gap0);
		fsc4242Gap0[fBin]->Fill(NtrksAfter, sc4242Gap0Re, Dn4Gap0);

		TComplex sc3232Gap0 = FourGap0(3, 2, -3, -2);
		double sc3232Gap0Re = sc3232Gap0.Re()/Dn4Gap0;
		fChsc3232Gap0->Fill(NtrksAfter, sc3232Gap0Re, Dn4Gap0);
		fsc3232Gap0[fBin]->Fill(NtrksAfter, sc3232Gap0Re, Dn4Gap0);

	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 1 && NtrksAfter3subR > 0 && Dn4_3sub != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub = Four_3SubEvts(4, 2, -4, -2);
		double sc4242_3subRe = sc4242_3sub.Re()/Dn4_3sub;
		fChsc4242_3sub->Fill(NtrksAfter, sc4242_3subRe, Dn4_3sub);
		fsc4242_3sub[fBin]->Fill(NtrksAfter, sc4242_3subRe, Dn4_3sub);

		TComplex sc3232_3sub = Four_3SubEvts(3, 2, -3, -2);
		double sc3232_3subRe = sc3232_3sub.Re()/Dn4_3sub;
		fChsc3232_3sub->Fill(NtrksAfter, sc3232_3subRe, Dn4_3sub);
		fsc3232_3sub[fBin]->Fill(NtrksAfter, sc3232_3subRe, Dn4_3sub);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub = Four_3SubEvts(4, 2, -2, -4);
		double sc4224_3subRe = sc4224_3sub.Re()/Dn4_3sub;
		fChsc4224_3sub->Fill(NtrksAfter, sc4224_3subRe, Dn4_3sub);
		fsc4224_3sub[fBin]->Fill(NtrksAfter, sc4224_3subRe, Dn4_3sub);

		TComplex sc3223_3sub = Four_3SubEvts(3, 2, -2, -3);
		double sc3223_3subRe = sc3223_3sub.Re()/Dn4_3sub;
		fChsc3223_3sub->Fill(NtrksAfter, sc3223_3subRe, Dn4_3sub);
		fsc3223_3sub[fBin]->Fill(NtrksAfter, sc3223_3subRe, Dn4_3sub);

	}

	//	perm 2
	if(NtrksAfter3subL > 1 && NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn4_3sub_perm2 != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub_perm2 = Four_3SubEvts_perm2(4, 2, -4, -2);
		double sc4242_3subRe_perm2 = sc4242_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc4242_3sub_perm2->Fill(NtrksAfter, sc4242_3subRe_perm2, Dn4_3sub_perm2);
		fsc4242_3sub_perm2[fBin]->Fill(NtrksAfter, sc4242_3subRe_perm2, Dn4_3sub_perm2);

		TComplex sc3232_3sub_perm2 = Four_3SubEvts_perm2(3, 2, -3, -2);
		double sc3232_3subRe_perm2 = sc3232_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc3232_3sub_perm2->Fill(NtrksAfter, sc3232_3subRe_perm2, Dn4_3sub_perm2);
		fsc3232_3sub_perm2[fBin]->Fill(NtrksAfter, sc3232_3subRe_perm2, Dn4_3sub_perm2);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub_perm2 = Four_3SubEvts_perm2(4, 2, -2, -4);
		double sc4224_3subRe_perm2 = sc4224_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc4224_3sub_perm2->Fill(NtrksAfter, sc4224_3subRe_perm2, Dn4_3sub_perm2);
		fsc4224_3sub_perm2[fBin]->Fill(NtrksAfter, sc4224_3subRe_perm2, Dn4_3sub_perm2);

		TComplex sc3223_3sub_perm2 = Four_3SubEvts_perm2(3, 2, -2, -3);
		double sc3223_3subRe_perm2 = sc3223_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc3223_3sub_perm2->Fill(NtrksAfter, sc3223_3subRe_perm2, Dn4_3sub_perm2);
		fsc3223_3sub_perm2[fBin]->Fill(NtrksAfter, sc3223_3subRe_perm2, Dn4_3sub_perm2);

	}

	//	perm 3
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && NtrksAfter3subR > 1 && Dn4_3sub_perm3 != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub_perm3 = Four_3SubEvts_perm3(4, 2, -4, -2);
		double sc4242_3subRe_perm3 = sc4242_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc4242_3sub_perm3->Fill(NtrksAfter, sc4242_3subRe_perm3, Dn4_3sub_perm3);
		fsc4242_3sub_perm3[fBin]->Fill(NtrksAfter, sc4242_3subRe_perm3, Dn4_3sub_perm3);

		TComplex sc3232_3sub_perm3 = Four_3SubEvts_perm3(3, 2, -3, -2);
		double sc3232_3subRe_perm3 = sc3232_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc3232_3sub_perm3->Fill(NtrksAfter, sc3232_3subRe_perm3, Dn4_3sub_perm3);
		fsc3232_3sub_perm3[fBin]->Fill(NtrksAfter, sc3232_3subRe_perm3, Dn4_3sub_perm3);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub_perm3 = Four_3SubEvts_perm3(4, 2, -2, -4);
		double sc4224_3subRe_perm3 = sc4224_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc4224_3sub_perm3->Fill(NtrksAfter, sc4224_3subRe_perm3, Dn4_3sub_perm3);
		fsc4224_3sub_perm3[fBin]->Fill(NtrksAfter, sc4224_3subRe_perm3, Dn4_3sub_perm3);

		TComplex sc3223_3sub_perm3 = Four_3SubEvts_perm3(3, 2, -2, -3);
		double sc3223_3subRe_perm3 = sc3223_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc3223_3sub_perm3->Fill(NtrksAfter, sc3223_3subRe_perm3, Dn4_3sub_perm3);
		fsc3223_3sub_perm3[fBin]->Fill(NtrksAfter, sc3223_3subRe_perm3, Dn4_3sub_perm3);

	}

	//..calculate 6-particle correlations	
	//...................................
	double Dn6 = Six(0, 0, 0, 0, 0, 0).Re();
	double Dn6Gap0 = SixGap0(0, 0, 0, 0, 0, 0).Re();

	if(NtrksAfter > 5 && Dn6 != 0)
	{

		TComplex v26 = Six(2, 2, 2, -2, -2, -2);	
		double v26Re = v26.Re()/Dn6;
		fChcn6Ntrks1bin[0]->Fill(NtrksAfter, v26Re, Dn6);
		fcn6Ntrks1bin[0][fBin]->Fill(NtrksAfter, v26Re, Dn6);

		TComplex v36 = Six(3, 3, 3, -3, -3, -3);
		double v36Re = v36.Re()/Dn6;
		fChcn6Ntrks1bin[1]->Fill(NtrksAfter, v36Re, Dn6);
		fcn6Ntrks1bin[1][fBin]->Fill(NtrksAfter, v36Re, Dn6);

		TComplex v46 = Six(4, 4, 4, -4, -4, -4);
		double v46Re = v46.Re()/Dn6;
		fChcn6Ntrks1bin[2]->Fill(NtrksAfter, v46Re, Dn6);
		fcn6Ntrks1bin[2][fBin]->Fill(NtrksAfter, v46Re, Dn6);

	}

	if(NtrksAfterGap0M > 2 && NtrksAfterGap0P > 2 && Dn6Gap0 != 0)
	{

		TComplex v26Gap0 = SixGap0(2, 2, 2, -2, -2, -2);	
		double v26Gap0Re = v26Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[0]->Fill(NtrksAfter, v26Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, v26Gap0Re, Dn6Gap0);

		TComplex v36Gap0 = SixGap0(3, 3, 3, -3, -3, -3);
		double v36Gap0Re = v36Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[1]->Fill(NtrksAfter, v36Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, v36Gap0Re, Dn6Gap0);

		TComplex v46Gap0 = SixGap0(4, 4, 4, -4, -4, -4);
		double v46Gap0Re = v46Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[2]->Fill(NtrksAfter, v46Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, v46Gap0Re, Dn6Gap0);

	}

  //..calculate 8-particle correlations 
  //...................................
	double Dn8 = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
  double Dn8Gap0 = EightGap0(0, 0, 0, 0, 0, 0, 0, 0).Re();

  if(NtrksAfter > 7 && Dn8 != 0)                                                                                                                                                                    
  {

    TComplex v28 = Eight(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Re = v28.Re()/Dn8;
    fChcn8Ntrks1bin[0]->Fill(NtrksAfter, v28Re, Dn8);
    fcn8Ntrks1bin[0][fBin]->Fill(NtrksAfter, v28Re, Dn8);

    TComplex v38 = Eight(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Re = v38.Re()/Dn8;
    fChcn8Ntrks1bin[1]->Fill(NtrksAfter, v38Re, Dn8);
    fcn8Ntrks1bin[1][fBin]->Fill(NtrksAfter, v38Re, Dn8);

    TComplex v48 = Eight(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Re = v48.Re()/Dn8;
    fChcn8Ntrks1bin[2]->Fill(NtrksAfter, v48Re, Dn8);
    fcn8Ntrks1bin[2][fBin]->Fill(NtrksAfter, v48Re, Dn8);

  } 

  if(NtrksAfterGap0M > 3 && NtrksAfterGap0P > 3 && Dn8Gap0 != 0)                                                                                                                                                                    
  {

    TComplex v28Gap0 = EightGap0(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Gap0Re = v28Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[0]->Fill(NtrksAfter, v28Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, v28Gap0Re, Dn8Gap0);

    TComplex v38Gap0 = EightGap0(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Gap0Re = v38Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[1]->Fill(NtrksAfter, v38Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, v38Gap0Re, Dn8Gap0);

    TComplex v48Gap0 = EightGap0(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Gap0Re = v48Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[2]->Fill(NtrksAfter, v48Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, v48Gap0Re, Dn8Gap0);

  } 


}
//____________________________________________________________________
double AliAnalysisTaskChargedFlowGF::ProcessMCTruth(AliVEvent *aod, float fVtxZ)
{

	TClonesArray* farray = (TClonesArray*)aod->FindListObject("mcparticles");
	if (!farray) return 0;

	Int_t ntrks = farray->GetEntries();

	double NtrksAfterMC = 0;
	double NtrksAfterGap0MMC = 0;
	double NtrksAfterGap0PMC = 0;
	double NtrksAfterGap2MMC = 0;
	double NtrksAfterGap2PMC = 0;
	double NtrksAfterGap4MMC = 0;
	double NtrksAfterGap4PMC = 0;
	double NtrksAfterGap8MMC = 0;
	double NtrksAfterGap8PMC = 0;
	double NtrksAfterGapMMC = 0;
	double NtrksAfterGapPMC = 0;
	double NtrksAfterGap14MMC = 0;
	double NtrksAfterGap14PMC = 0;
	double NtrksAfter3subLMC = 0;
	double NtrksAfter3subMMC = 0;
	double NtrksAfter3subRMC = 0;

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGapM[20][20] = {0};
	double QsinGapM[20][20] = {0};
	double QcosGapP[20][20] = {0};
	double QsinGapP[20][20] = {0};

	double QcosGap0M[20][20] = {0};
	double QsinGap0M[20][20] = {0};
	double QcosGap0P[20][20] = {0};
	double QsinGap0P[20][20] = {0};
  
	double QcosGap2M[20][20] = {0};
	double QsinGap2M[20][20] = {0};
	double QcosGap2P[20][20] = {0};
	double QsinGap2P[20][20] = {0};
  
	double QcosGap4M[20][20] = {0};
	double QsinGap4M[20][20] = {0};
	double QcosGap4P[20][20] = {0};
	double QsinGap4P[20][20] = {0};
	
	double QcosGap8M[20][20] = {0};
	double QsinGap8M[20][20] = {0};
	double QcosGap8P[20][20] = {0};
	double QsinGap8P[20][20] = {0};
    
	double QcosGap14M[20][20] = {0};
	double QsinGap14M[20][20] = {0};
	double QcosGap14P[20][20] = {0};
	double QsinGap14P[20][20] = {0};
    
	double QcosSubLeft[20][20] = {0};
	double QsinSubLeft[20][20] = {0};
	double QcosSubMiddle[20][20] = {0};
	double QsinSubMiddle[20][20] = {0};
	double QcosSubRight[20][20] = {0};
	double QsinSubRight[20][20] = {0};
    
    
  Int_t fBin = 0;
  TRandom3 rr(0);
  Double_t ranNum = 0.;

		if(fSample == 1) fBin = 1;
    
    if(fSample == 10)
		{
    	ranNum = rr.Rndm();
        
    	if (ranNum <= 0.1 && ranNum > 0. ) fBin = 1;
    	else if (ranNum <= 0.2 && ranNum > 0.1) fBin = 2;
    	else if (ranNum <= 0.3 && ranNum > 0.2) fBin = 3;
    	else if (ranNum <= 0.4 && ranNum > 0.3) fBin = 4;
    	else if (ranNum <= 0.5 && ranNum > 0.4) fBin = 5;
    	else if (ranNum <= 0.6 && ranNum > 0.5) fBin = 6;
    	else if (ranNum <= 0.7 && ranNum > 0.6) fBin = 7;
    	else if (ranNum <= 0.8 && ranNum > 0.7) fBin = 8;
    	else if (ranNum <= 0.9 && ranNum > 0.8) fBin = 9;
    	else if (ranNum > 0.9 && ranNum <= 1.) fBin = 10;
    	else {
         cout << ranNum << endl;
         
         cout << "This is very strange!!!" << endl;
         fBin = 11;
    	}

    }

	//..loop to get Ntrks	
  for(Int_t nt = 0; nt < ntrks; nt++)
    {
      AliAODMCParticle *trk = (AliAODMCParticle*) farray->At(nt);
			
       if (!(trk->IsPhysicalPrimary())){
         continue;
			}
	
			if(trk->Charge() == 0) continue;

			if(TMath::Abs(trk->Eta()) > fEtaCut) continue;
			if(trk->Pt() < fMinPt) continue;
			if(trk->Pt() > fMaxPt) continue;

			NtrksAfterMC += 1;

		}

	hMultMC->Fill(NtrksAfterMC);

	double countMC = 0;
	double vtxz = fInputEvent->GetPrimaryVertex()->GetZ();
	//..LOOP OVER TRACKS........	
	//........................................
  for(Int_t nt = 0; nt < ntrks; nt++)
    {
      AliAODMCParticle *trk = (AliAODMCParticle*) farray->At(nt);
			
      if (!(trk->IsPhysicalPrimary())){
         continue;
			}
  
			if(trk->Charge() == 0) continue;

			if(TMath::Abs(trk->Eta()) > fEtaCut) continue;
			if(trk->Pt() < fMinPt) continue;
			if(trk->Pt() > fMaxPt) continue;

			double posZ = trk->Zv();
			double dca = posZ - vtxz;
	
			if(trk->IsSecondaryFromMaterial()) hDCAptMC_material->Fill(trk->Pt(), dca);
			if(trk->IsSecondaryFromWeakDecay()) hDCAptMC_weak->Fill(trk->Pt(), dca);

			hDCAptMC->Fill(trk->Pt(), dca);

			countMC += 1;

			//..fill phi distribution
			fPhiDisTruth->Fill(trk->Phi());
			fEtaDisTruth->Fill(trk->Eta());
			fPtDisTruth->Fill(trk->Pt());

			//..for efficiency     
			hTruth->Fill(trk->Pt(), trk->Eta(), fVtxZ);

			double weight = 1;

			//..calculate Q-vectors
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					Qcos[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
					Qsin[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
				}
			}


      //Gap = 0.0
      if (trk->Eta() > 0.)
      {
        NtrksAfterGap0PMC += 1; 
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap0P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
 
      }
      if (trk->Eta() < -0.)
      {
				NtrksAfterGap0MMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap0M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
		
      //Gap = 0.2
      if (trk->Eta() > 0.1)
      {
        NtrksAfterGap2PMC += 1; 
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
 
      }
      if (trk->Eta() < -0.1)
      {
				NtrksAfterGap2MMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
		
      // Gap = 0.4
      if (trk->Eta() > 0.2)
      {
				NtrksAfterGap4PMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap4P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
      if (trk->Eta() < -0.2)
      {
				NtrksAfterGap4MMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap4M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
	 	
      // Gap = 0.8
      if (trk->Eta() > 0.4)
      {
				NtrksAfterGap8PMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
          
      if (trk->Eta() < -0.4)
      {
				NtrksAfterGap8MMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}

      }


			//..if eta < -0.5
			if(trk->Eta() < -0.5)
			{
				NtrksAfterGapMMC += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapM[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGapM[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
			}
	
			//..if eta > 0.5
			if(trk->Eta() > 0.5)
			{
				NtrksAfterGapPMC += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapP[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGapP[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());	
					}
				}
			} 

      // Gap = 1.4
      if (trk->Eta() > 0.7)
      {
				NtrksAfterGap14PMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap14P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
      }
          
      if (trk->Eta() < -0.7)
      {
				NtrksAfterGap14MMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinGap14M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}

      }
       
			//..3-subevent method
			if(trk->Eta() < -0.4)
			{//..left part
				NtrksAfter3subLMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubLeft[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinSubLeft[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
			}
			if(trk->Eta() >= -0.4 && trk->Eta() <= 0.4)
			{//..middle part
				NtrksAfter3subMMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubMiddle[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinSubMiddle[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
			}
			if(trk->Eta() > 0.4)
			{//..right part
				NtrksAfter3subRMC += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubRight[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*trk->Phi());
						QsinSubRight[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*trk->Phi());
					}
				}
			}

  } // end loop of all track

	//............................
	//..GENERIC FRAMEWORK RP
	//............................

	//..calculate Q-vector for each harmonics n and power p
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);
			QvectorM[iharm][ipow] = TComplex(QcosGapM[iharm][ipow], QsinGapM[iharm][ipow]);
			QvectorP[iharm][ipow] = TComplex(QcosGapP[iharm][ipow], QsinGapP[iharm][ipow]);
			Qvector0M[iharm][ipow] = TComplex(QcosGap0M[iharm][ipow], QsinGap0M[iharm][ipow]);
			Qvector0P[iharm][ipow] = TComplex(QcosGap0P[iharm][ipow], QsinGap0P[iharm][ipow]);
			Qvector2M[iharm][ipow] = TComplex(QcosGap2M[iharm][ipow], QsinGap2M[iharm][ipow]);
			Qvector2P[iharm][ipow] = TComplex(QcosGap2P[iharm][ipow], QsinGap2P[iharm][ipow]);
			Qvector4M[iharm][ipow] = TComplex(QcosGap4M[iharm][ipow], QsinGap4M[iharm][ipow]);
			Qvector4P[iharm][ipow] = TComplex(QcosGap4P[iharm][ipow], QsinGap4P[iharm][ipow]);
			Qvector8M[iharm][ipow] = TComplex(QcosGap8M[iharm][ipow], QsinGap8M[iharm][ipow]);
			Qvector8P[iharm][ipow] = TComplex(QcosGap8P[iharm][ipow], QsinGap8P[iharm][ipow]);
			Qvector14M[iharm][ipow] = TComplex(QcosGap14M[iharm][ipow], QsinGap14M[iharm][ipow]);
			Qvector14P[iharm][ipow] = TComplex(QcosGap14P[iharm][ipow], QsinGap14P[iharm][ipow]);
			QvectorSubLeft[iharm][ipow] = TComplex(QcosSubLeft[iharm][ipow], QsinSubLeft[iharm][ipow]);
			QvectorSubRight[iharm][ipow] = TComplex(QcosSubRight[iharm][ipow], QsinSubRight[iharm][ipow]);
			QvectorSubMiddle[iharm][ipow] = TComplex(QcosSubMiddle[iharm][ipow], QsinSubMiddle[iharm][ipow]);
		}
	}

	//..calculate 2-particle correlations
	//..................................
	double Dn2 = Two(0, 0).Re();
	double Dn2Gap = TwoGap(0, 0).Re();
	double Dn2Gap0 = TwoGap0(0, 0).Re();
	double Dn2Gap2 = TwoGap2(0, 0).Re();
	double Dn2Gap4 = TwoGap4(0, 0).Re();
	double Dn2Gap8 = TwoGap8(0, 0).Re();
	double Dn2Gap14 = TwoGap14(0, 0).Re();
	double Dn2_3subLM = Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = Two_3SubRM(0, 0).Re();

	if(NtrksAfterMC > 1 && Dn2 != 0)
	{	

		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		fChMCcn2Ntrks1bin[0]->Fill(NtrksAfterMC, v22Re, Dn2);
		fMCcn2Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		fChMCcn2Ntrks1bin[1]->Fill(NtrksAfterMC, v32Re, Dn2);
		fMCcn2Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		fChMCcn2Ntrks1bin[2]->Fill(NtrksAfterMC, v42Re, Dn2);
		fMCcn2Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42Re, Dn2);

	}	

	if(NtrksAfterGap0MMC > 0 && NtrksAfterGap0PMC > 0 && Dn2Gap0 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap0 = TwoGap0(2, -2);
		double v22ReGap0 = v22Gap0.Re()/Dn2Gap0;
		fChMCcn2Gap0Ntrks1bin[0]->Fill(NtrksAfterMC, v22ReGap0, Dn2Gap0);
		fMCcn2Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap0, Dn2Gap0);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap0 = TwoGap0(3, -3);
		double v32ReGap0 = v32Gap0.Re()/Dn2Gap0;
		fChMCcn2Gap0Ntrks1bin[1]->Fill(NtrksAfterMC, v32ReGap0, Dn2Gap0);
		fMCcn2Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap0, Dn2Gap0);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap0 = TwoGap0(4, -4);
		double v42ReGap0 = v42Gap0.Re()/Dn2Gap0;
		fChMCcn2Gap0Ntrks1bin[2]->Fill(NtrksAfterMC, v42ReGap0, Dn2Gap0);
		fMCcn2Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap0, Dn2Gap0);

	}

	if(NtrksAfterGap2MMC > 0 && NtrksAfterGap2PMC > 0 && Dn2Gap2 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap2 = TwoGap2(2, -2);
		double v22ReGap2 = v22Gap2.Re()/Dn2Gap2;
		fChMCcn2Gap2Ntrks1bin[0]->Fill(NtrksAfterMC, v22ReGap2, Dn2Gap2);
		fMCcn2Gap2Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap2, Dn2Gap2);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap2 = TwoGap2(3, -3);
		double v32ReGap2 = v32Gap2.Re()/Dn2Gap2;
		fChMCcn2Gap2Ntrks1bin[1]->Fill(NtrksAfterMC, v32ReGap2, Dn2Gap2);
		fMCcn2Gap2Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap2, Dn2Gap2);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap2 = TwoGap2(4, -4);
		double v42ReGap2 = v42Gap2.Re()/Dn2Gap2;
		fChMCcn2Gap2Ntrks1bin[2]->Fill(NtrksAfterMC, v42ReGap2, Dn2Gap2);
		fMCcn2Gap2Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap2, Dn2Gap2);

	}

	if(NtrksAfterGap4MMC > 0 && NtrksAfterGap4PMC > 0 && Dn2Gap4 != 0)
	{
		//..v2{2} with eta Gap > 0.4
		TComplex v22Gap4 = TwoGap4(2, -2);
		double v22ReGap4 = v22Gap4.Re()/Dn2Gap4;
		fChMCcn2Gap4Ntrks1bin[0]->Fill(NtrksAfterMC, v22ReGap4, Dn2Gap4);
		fMCcn2Gap4Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap4, Dn2Gap4);

		//..v3{2} with eta Gap > 0.4
		TComplex v32Gap4 = TwoGap4(3, -3);
		double v32ReGap4 = v32Gap4.Re()/Dn2Gap4;
		fChMCcn2Gap4Ntrks1bin[1]->Fill(NtrksAfterMC, v32ReGap4, Dn2Gap4);
		fMCcn2Gap4Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap4, Dn2Gap4);

		//..v4{2} with eta Gap > 0.4
		TComplex v42Gap4 = TwoGap4(4, -4);
		double v42ReGap4 = v42Gap4.Re()/Dn2Gap4;
		fChMCcn2Gap4Ntrks1bin[2]->Fill(NtrksAfterMC, v42ReGap4, Dn2Gap4);
		fMCcn2Gap4Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap4, Dn2Gap4);

	}

	if(NtrksAfterGap8MMC > 0 && NtrksAfterGap0PMC > 0 && Dn2Gap8 != 0)
	{
		//..v2{2} with eta Gap > 0.8
		TComplex v22Gap8 = TwoGap8(2, -2);
		double v22ReGap8 = v22Gap8.Re()/Dn2Gap8;
		fChMCcn2Gap8Ntrks1bin[0]->Fill(NtrksAfterMC, v22ReGap8, Dn2Gap8);
		fMCcn2Gap8Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap8, Dn2Gap8);

		//..v3{2} with eta Gap > 0.8
		TComplex v32Gap8 = TwoGap8(3, -3);
		double v32ReGap8 = v32Gap8.Re()/Dn2Gap8;
		fChMCcn2Gap8Ntrks1bin[1]->Fill(NtrksAfterMC, v32ReGap8, Dn2Gap8);
		fMCcn2Gap8Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap8, Dn2Gap8);

		//..v4{2} with eta Gap > 0.8
		TComplex v42Gap8 = TwoGap8(4, -4);
		double v42ReGap8 = v42Gap8.Re()/Dn2Gap8;
		fChMCcn2Gap8Ntrks1bin[2]->Fill(NtrksAfterMC, v42ReGap8, Dn2Gap8);
		fMCcn2Gap8Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap8, Dn2Gap8);

	}

	if(NtrksAfterGapMMC > 0 && NtrksAfterGapPMC > 0 && Dn2Gap != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap = TwoGap(2, -2);
		double v22ReGap = v22Gap.Re()/Dn2Gap;
		fChMCcn2GapNtrks1bin[0]->Fill(NtrksAfterMC, v22ReGap, Dn2Gap);
		fMCcn2GapNtrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap, Dn2Gap);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap = TwoGap(3, -3);
		double v32ReGap = v32Gap.Re()/Dn2Gap;
		fChMCcn2GapNtrks1bin[1]->Fill(NtrksAfterMC, v32ReGap, Dn2Gap);
		fMCcn2GapNtrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap, Dn2Gap);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap = TwoGap(4, -4);
		double v42ReGap = v42Gap.Re()/Dn2Gap;
		fChMCcn2GapNtrks1bin[2]->Fill(NtrksAfterMC, v42ReGap, Dn2Gap);
		fMCcn2GapNtrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap, Dn2Gap);

	}

	if(NtrksAfterGap14MMC > 0 && NtrksAfterGap14PMC > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		fChMCcn2Gap14Ntrks1bin[0]->Fill(NtrksAfterMC, v22ReGap14, Dn2Gap14);
		fMCcn2Gap14Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v22ReGap14, Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		fChMCcn2Gap14Ntrks1bin[1]->Fill(NtrksAfterMC, v32ReGap14, Dn2Gap14);
		fMCcn2Gap14Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v32ReGap14, Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		fChMCcn2Gap14Ntrks1bin[2]->Fill(NtrksAfterMC, v42ReGap14, Dn2Gap14);
		fMCcn2Gap14Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v42ReGap14, Dn2Gap14);

	}

	//..for 3-subevent method, Gap0
	if(NtrksAfter3subMMC > 1 && NtrksAfter3subLMC > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		fChMCcn2_3subLMNtrks1bin[0]->Fill(NtrksAfterMC, v22Re_3subLM, Dn2_3subLM);
		fMCcn2_3subLMNtrks1bin[0][fBin]->Fill(NtrksAfterMC, v22Re_3subLM, Dn2_3subLM);

		TComplex v32_3subLM = Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		fChMCcn2_3subLMNtrks1bin[1]->Fill(NtrksAfterMC, v32Re_3subLM, Dn2_3subLM);
		fMCcn2_3subLMNtrks1bin[1][fBin]->Fill(NtrksAfterMC, v32Re_3subLM, Dn2_3subLM);

		TComplex v42_3subLM = Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		fChMCcn2_3subLMNtrks1bin[2]->Fill(NtrksAfterMC, v42Re_3subLM, Dn2_3subLM);
		fMCcn2_3subLMNtrks1bin[2][fBin]->Fill(NtrksAfterMC, v42Re_3subLM, Dn2_3subLM);

	}

	if(NtrksAfter3subRMC > 0 && NtrksAfter3subMMC > 1 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		fChMCcn2_3subRMNtrks1bin[0]->Fill(NtrksAfterMC, v22Re_3subRM, Dn2_3subRM);
		fMCcn2_3subRMNtrks1bin[0][fBin]->Fill(NtrksAfterMC, v22Re_3subRM, Dn2_3subRM);

		TComplex v32_3subRM = Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		fChMCcn2_3subRMNtrks1bin[1]->Fill(NtrksAfterMC, v32Re_3subRM, Dn2_3subRM);
		fMCcn2_3subRMNtrks1bin[1][fBin]->Fill(NtrksAfterMC, v32Re_3subRM, Dn2_3subRM);

		TComplex v42_3subRM = Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		fChMCcn2_3subRMNtrks1bin[2]->Fill(NtrksAfterMC, v42Re_3subRM, Dn2_3subRM);
		fMCcn2_3subRMNtrks1bin[2][fBin]->Fill(NtrksAfterMC, v42Re_3subRM, Dn2_3subRM);

	}

	//..calculate 4-particle correlations
	//................................
	double Dn4 = Four(0, 0, 0, 0).Re();
	double Dn4Gap0 = FourGap0(0, 0, 0, 0).Re();
	double Dn4_3sub = Four_3SubEvts(0, 0, 0, 0).Re();

	if(NtrksAfterMC > 3 && Dn4 != 0)
	{
		
		TComplex v24 = Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		fChMCcn4Ntrks1bin[0]->Fill(NtrksAfterMC, v24Re, Dn4);
		fMCcn4Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v24Re, Dn4);

		TComplex v34 = Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		fChMCcn4Ntrks1bin[1]->Fill(NtrksAfterMC, v34Re, Dn4);
		fMCcn4Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v34Re, Dn4);

		TComplex v44 = Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		fChMCcn4Ntrks1bin[2]->Fill(NtrksAfterMC, v44Re, Dn4);
		fMCcn4Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v44Re, Dn4);

	}

	if(NtrksAfterGap0MMC > 1 && NtrksAfterGap0PMC > 1 && Dn4Gap0 !=0)
	{

		TComplex v24Gap0 = FourGap0(2, 2, -2, -2);
		double v24Gap0Re = v24Gap0.Re()/Dn4Gap0;
		fChMCcn4Gap0Ntrks1bin[0]->Fill(NtrksAfterMC, v24Gap0Re, Dn4Gap0);
		fMCcn4Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v24Gap0Re, Dn4Gap0);

		TComplex v34Gap0 = FourGap0(3, 3, -3, -3);
		double v34Gap0Re = v34Gap0.Re()/Dn4Gap0;
		fChMCcn4Gap0Ntrks1bin[1]->Fill(NtrksAfterMC, v34Gap0Re, Dn4Gap0);
		fMCcn4Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v34Gap0Re, Dn4Gap0);

		TComplex v44Gap0 = FourGap0(4, 4, -4, -4);
		double v44Gap0Re = v44Gap0.Re()/Dn4Gap0;
		fChMCcn4Gap0Ntrks1bin[2]->Fill(NtrksAfterMC, v44Gap0Re, Dn4Gap0);
		fMCcn4Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v44Gap0Re, Dn4Gap0);

	}

	//..3-subevent method
	if(NtrksAfter3subLMC > 0 && NtrksAfter3subRMC > 0 && NtrksAfter3subMMC > 1 && Dn4_3sub != 0)
	{
		TComplex v24_3sub = Four_3SubEvts(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3sub;
		fChMCcn4_3subNtrks1bin[0]->Fill(NtrksAfterMC, v24_3subRe, Dn4_3sub);
		fMCcn4_3subNtrks1bin[0][fBin]->Fill(NtrksAfterMC, v24_3subRe, Dn4_3sub);

		TComplex v34_3sub = Four_3SubEvts(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3sub;
		fChMCcn4_3subNtrks1bin[1]->Fill(NtrksAfterMC, v34_3subRe, Dn4_3sub);
		fMCcn4_3subNtrks1bin[1][fBin]->Fill(NtrksAfterMC, v34_3subRe, Dn4_3sub);

		TComplex v44_3sub = Four_3SubEvts(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3sub;
		fChMCcn4_3subNtrks1bin[2]->Fill(NtrksAfterMC, v44_3subRe, Dn4_3sub);
		fMCcn4_3subNtrks1bin[2][fBin]->Fill(NtrksAfterMC, v44_3subRe, Dn4_3sub);
	}


	//..SC(n,m)
	//................................
	if(NtrksAfterMC > 3 && Dn4 != 0)
	{
		//..SC(4,2,-4,-2)	
		TComplex sc4242 = Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		fChMCsc4242->Fill(NtrksAfterMC, sc4242Re, Dn4);
		fMCsc4242[fBin]->Fill(NtrksAfterMC, sc4242Re, Dn4);

		//..SC(3,2,-3,-2)	
		TComplex sc3232 = Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		fChMCsc3232->Fill(NtrksAfterMC, sc3232Re, Dn4);
		fMCsc3232[fBin]->Fill(NtrksAfterMC, sc3232Re, Dn4);

	}	

	if(NtrksAfterGap0MMC > 1 && NtrksAfterGap0PMC > 1 && Dn4Gap0 != 0)
	{

		TComplex sc4242Gap0 = FourGap0(4, 2, -4, -2);
		double sc4242Gap0Re = sc4242Gap0.Re()/Dn4Gap0;
		fChMCsc4242Gap0->Fill(NtrksAfterMC, sc4242Gap0Re, Dn4Gap0);
		fMCsc4242Gap0[fBin]->Fill(NtrksAfterMC, sc4242Gap0Re, Dn4Gap0);

		TComplex sc3232Gap0 = FourGap0(3, 2, -3, -2);
		double sc3232Gap0Re = sc3232Gap0.Re()/Dn4Gap0;
		fChMCsc3232Gap0->Fill(NtrksAfterMC, sc3232Gap0Re, Dn4Gap0);
		fMCsc3232Gap0[fBin]->Fill(NtrksAfterMC, sc3232Gap0Re, Dn4Gap0);

	}


	if(NtrksAfter3subLMC > 0 && NtrksAfter3subRMC > 0 && NtrksAfter3subMMC > 1 && Dn4_3sub != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub = Four_3SubEvts(4, 2, -4, -2);
		double sc4242_3subRe = sc4242_3sub.Re()/Dn4_3sub;
		fChMCsc4242_3sub->Fill(NtrksAfterMC, sc4242_3subRe, Dn4_3sub);
		fMCsc4242_3sub[fBin]->Fill(NtrksAfterMC, sc4242_3subRe, Dn4_3sub);

		TComplex sc3232_3sub = Four_3SubEvts(3, 2, -3, -2);
		double sc3232_3subRe = sc3232_3sub.Re()/Dn4_3sub;
		fChMCsc3232_3sub->Fill(NtrksAfterMC, sc3232_3subRe, Dn4_3sub);
		fMCsc3232_3sub[fBin]->Fill(NtrksAfterMC, sc3232_3subRe, Dn4_3sub);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub = Four_3SubEvts(4, 2, -2, -4);
		double sc4224_3subRe = sc4224_3sub.Re()/Dn4_3sub;
		fChMCsc4224_3sub->Fill(NtrksAfterMC, sc4224_3subRe, Dn4_3sub);
		fMCsc4224_3sub[fBin]->Fill(NtrksAfterMC, sc4224_3subRe, Dn4_3sub);

		TComplex sc3223_3sub = Four_3SubEvts(3, 2, -2, -3);
		double sc3223_3subRe = sc3223_3sub.Re()/Dn4_3sub;
		fChMCsc3223_3sub->Fill(NtrksAfterMC, sc3223_3subRe, Dn4_3sub);
		fMCsc3223_3sub[fBin]->Fill(NtrksAfterMC, sc3223_3subRe, Dn4_3sub);

	}


	//..calculate 6-particle correlations	
	//...................................
	double Dn6 = Six(0, 0, 0, 0, 0, 0).Re();
	double Dn6Gap0 = SixGap0(0, 0, 0, 0, 0, 0).Re();

	if(NtrksAfterMC > 5 && Dn6 != 0)
	{

		TComplex v26 = Six(2, 2, 2, -2, -2, -2);	
		double v26Re = v26.Re()/Dn6;
		fChMCcn6Ntrks1bin[0]->Fill(NtrksAfterMC, v26Re, Dn6);
		fMCcn6Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v26Re, Dn6);

		TComplex v36 = Six(3, 3, 3, -3, -3, -3);
		double v36Re = v36.Re()/Dn6;
		fChMCcn6Ntrks1bin[1]->Fill(NtrksAfterMC, v36Re, Dn6);
		fMCcn6Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v36Re, Dn6);

		TComplex v46 = Six(4, 4, 4, -4, -4, -4);
		double v46Re = v46.Re()/Dn6;
		fChMCcn6Ntrks1bin[2]->Fill(NtrksAfterMC, v46Re, Dn6);
		fMCcn6Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v46Re, Dn6);

	}

	if(NtrksAfterGap0MMC > 2 && NtrksAfterGap0PMC > 2 && Dn6Gap0 != 0)
	{

		TComplex v26Gap0 = SixGap0(2, 2, 2, -2, -2, -2);	
		double v26Gap0Re = v26Gap0.Re()/Dn6Gap0;
		fChMCcn6Gap0Ntrks1bin[0]->Fill(NtrksAfterMC, v26Gap0Re, Dn6Gap0);
		fMCcn6Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v26Gap0Re, Dn6Gap0);

		TComplex v36Gap0 = SixGap0(3, 3, 3, -3, -3, -3);
		double v36Gap0Re = v36Gap0.Re()/Dn6Gap0;
		fChMCcn6Gap0Ntrks1bin[1]->Fill(NtrksAfterMC, v36Gap0Re, Dn6Gap0);
		fMCcn6Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v36Gap0Re, Dn6Gap0);

		TComplex v46Gap0 = SixGap0(4, 4, 4, -4, -4, -4);
		double v46Gap0Re = v46Gap0.Re()/Dn6Gap0;
		fChMCcn6Gap0Ntrks1bin[2]->Fill(NtrksAfterMC, v46Gap0Re, Dn6Gap0);
		fMCcn6Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v46Gap0Re, Dn6Gap0);

	}

  //..calculate 8-particle correlations 
  //...................................
  double Dn8 = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
  double Dn8Gap0 = EightGap0(0, 0, 0, 0, 0, 0, 0, 0).Re();

  if(NtrksAfterMC > 7 && Dn8 != 0)                                                                                                                                                                    
  {

    TComplex v28 = Eight(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Re = v28.Re()/Dn8;
    fChMCcn8Ntrks1bin[0]->Fill(NtrksAfterMC, v28Re, Dn8);
    fMCcn8Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v28Re, Dn8);

    TComplex v38 = Eight(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Re = v38.Re()/Dn8;
    fChMCcn8Ntrks1bin[1]->Fill(NtrksAfterMC, v38Re, Dn8);
    fMCcn8Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v38Re, Dn8);

    TComplex v48 = Eight(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Re = v48.Re()/Dn8;
    fChMCcn8Ntrks1bin[2]->Fill(NtrksAfterMC, v48Re, Dn8);
    fMCcn8Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v48Re, Dn8);

  } 

  if(NtrksAfterGap0MMC > 3 && NtrksAfterGap0PMC > 3 && Dn8Gap0 != 0)                                                                                                                                                                    
  {

    TComplex v28Gap0 = EightGap0(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Gap0Re = v28Gap0.Re()/Dn8Gap0;
    fChMCcn8Gap0Ntrks1bin[0]->Fill(NtrksAfterMC, v28Gap0Re, Dn8Gap0);
    fMCcn8Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfterMC, v28Gap0Re, Dn8Gap0);

    TComplex v38Gap0 = EightGap0(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Gap0Re = v38Gap0.Re()/Dn8Gap0;
    fChMCcn8Gap0Ntrks1bin[1]->Fill(NtrksAfterMC, v38Gap0Re, Dn8Gap0);
    fMCcn8Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfterMC, v38Gap0Re, Dn8Gap0);

    TComplex v48Gap0 = EightGap0(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Gap0Re = v48Gap0.Re()/Dn8Gap0;
    fChMCcn8Gap0Ntrks1bin[2]->Fill(NtrksAfterMC, v48Gap0Re, Dn8Gap0);
    fMCcn8Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfterMC, v48Gap0Re, Dn8Gap0);

  } 


	return NtrksAfterMC;

}
//____________________________________________________________________
double AliAnalysisTaskChargedFlowGF::GetPtWeight(double pt, double eta, float vz, double runNumber)
{

	double weight = 1;

	if(fPeriod == "LHC16q_CENT_wSDD" || fPeriod == "LHC16q_CENT_woSDD" || fPeriod == "LHC16q_FAST")
	{//	periods that need run-by-run efficiency
		hTrackEfficiencyRun = (TH3F*)fTrackEfficiency->Get(Form("eff_%.0lf", runNumber));	

		double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
		double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
		double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
		//..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
		double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
		double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

		if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
		else{
			TRandom3 r(0);
			double efficiency = 0;
			efficiency = r.Gaus(eff, error);
			weight = 1./efficiency; //..taking into account errors
			//weight = 1./eff;
		}
	}else if(fPeriod == "LHC15oHR" || fPeriod == "LHC15oLR")
	{//	periods that need run-by-run efficiency
		hTrackEfficiencyRun = (TH3F*)fTrackEfficiency->Get(Form("eff_LHC15o_AMPT_%.0lf", runNumber));	

		double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
		double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
		double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
		//..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
		double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
		double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

		if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
		else{
			TRandom3 r(0);
			double efficiency = 0;
			efficiency = r.Gaus(eff, error);
			weight = 1./efficiency; //..taking into account errors
			//weight = 1./eff;
		}
	}else{//	periods that doesn't need run-by-run efficiency
		double binPt = hTrackEfficiency->GetXaxis()->FindBin(pt);
		double binEta = hTrackEfficiency->GetYaxis()->FindBin(eta);
		double binVz = hTrackEfficiency->GetZaxis()->FindBin(vz);
		//..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
		double eff = hTrackEfficiency->GetBinContent(binPt, binEta, binVz);
		double error = hTrackEfficiency->GetBinError(binPt, binEta, binVz);

		if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
		else{
			TRandom3 r(0);
			double efficiency = 0;
			efficiency = r.Gaus(eff, error);
			weight = 1./efficiency; //..taking into account errors
			//weight = 1./eff;
		}
	}

	return weight;

}
//____________________________________________________________________
double AliAnalysisTaskChargedFlowGF::GetWeight(double phi, double eta, double vz, double runNumber)
{

	double weight = 1;

	if(fPeriod == "LHC15i")
	{
		weight = hPhiWeight->GetBinContent(hPhiWeight->GetXaxis()->FindBin(phi),
																				hPhiWeight->GetYaxis()->FindBin(eta),
																				hPhiWeight->GetZaxis()->FindBin(vz));
	}

	///////////////////////////////////////////////
	if(fPeriod == "LHC16k" || fPeriod == "LHC16l" || fPeriod == "LHC16p" || 
		fPeriod == "LHC16o" || fPeriod == "LHC16j" || fPeriod == "LHC16i" || 
		fPeriod == "LHC16h" || fPeriod == "LHC17m" || fPeriod == "LHC17o")
	{
		weight = hPhiWeight->GetBinContent(hPhiWeight->GetXaxis()->FindBin(phi),
																				hPhiWeight->GetYaxis()->FindBin(eta),
																				hPhiWeight->GetZaxis()->FindBin(vz));
	}


	///////////////////////////////////////////////
	if(fPeriod == "LHC16q_CENT_wSDD" || fPeriod == "LHC16q_CENT_woSDD" || fPeriod == "LHC16q_FAST"){
		hPhiWeightRun = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%0.lf", runNumber));
		weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
																					hPhiWeightRun->GetYaxis()->FindBin(eta),
																					hPhiWeightRun->GetZaxis()->FindBin(vz));
	}
	
	///////////////////////////////////////////////
	if(fPeriod == "LHC16t_CENT_wSDD" || fPeriod == "LHC16t_CENT_woSDD" || fPeriod == "LHC16t_FAST") {
		weight = hPhiWeight->GetBinContent(hPhiWeight->GetXaxis()->FindBin(phi),
		                                   hPhiWeight->GetYaxis()->FindBin(eta),
	   		                               hPhiWeight->GetZaxis()->FindBin(vz)); 
	}
	
	///////////////////////////////////////////////
	if(fPeriod == "LHC15oHR" || fPeriod == "LHC15oLR") {
		hPhiWeightRun = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%0.lf", runNumber));
		weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
				                                  hPhiWeightRun->GetYaxis()->FindBin(eta),
	      		                              hPhiWeightRun->GetZaxis()->FindBin(vz)); 
	}
	
	///////////////////////////////////////////////
	if(fPeriod == "LHC17n")
	{
		weight = hPhiWeight->GetBinContent(hPhiWeight->GetXaxis()->FindBin(phi),
		                                   hPhiWeight->GetYaxis()->FindBin(eta),
	   		                               hPhiWeight->GetZaxis()->FindBin(vz)); 
	}


	return weight;

}
//_____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGapM(int n, int p)
{

	if(n>=0) return QvectorM[n][p];
  else return TComplex::Conjugate(QvectorM[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGapP(int n, int p)
{

	if(n>=0) return QvectorP[n][p];
  else return TComplex::Conjugate(QvectorP[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap0M(int n, int p)
{

	if(n>=0) return Qvector0M[n][p];
  else return TComplex::Conjugate(Qvector0M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap0P(int n, int p)
{

	if(n>=0) return Qvector0P[n][p];
  else return TComplex::Conjugate(Qvector0P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap2M(int n, int p)
{

	if(n>=0) return Qvector2M[n][p];
  else return TComplex::Conjugate(Qvector2M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap2P(int n, int p)
{

	if(n>=0) return Qvector2P[n][p];
  else return TComplex::Conjugate(Qvector2P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap4M(int n, int p)
{

	if(n>=0) return Qvector4M[n][p];
  else return TComplex::Conjugate(Qvector4M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap4P(int n, int p)
{

	if(n>=0) return Qvector4P[n][p];
  else return TComplex::Conjugate(Qvector4P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap8M(int n, int p)
{

	if(n>=0) return Qvector8M[n][p];
  else return TComplex::Conjugate(Qvector8M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap8P(int n, int p)
{

	if(n>=0) return Qvector8P[n][p];
  else return TComplex::Conjugate(Qvector8P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap14M(int n, int p)
{

	if(n>=0) return Qvector14M[n][p];
  else return TComplex::Conjugate(Qvector14M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QGap14P(int n, int p)
{

	if(n>=0) return Qvector14P[n][p];
  else return TComplex::Conjugate(Qvector14P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QsubLeft(int n, int p)
{

	if(n>=0) return QvectorSubLeft[n][p];
  else return TComplex::Conjugate(QvectorSubLeft[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QsubRight(int n, int p)
{

	if(n>=0) return QvectorSubRight[n][p];
  else return TComplex::Conjugate(QvectorSubRight[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::QsubMiddle(int n, int p)
{

	if(n>=0) return QvectorSubMiddle[n][p];
  else return TComplex::Conjugate(QvectorSubMiddle[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two(int n1, int n2)
{

	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap0(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap2(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap4(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap8(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap(int n1, int n2)
{

	TComplex formula = QGapM(n1,1)*QGapP(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::TwoGap14(int n1, int n2)
{

	TComplex formula = QGap14M(n1,1)*QGap14P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubLM(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubRM(int n1, int n2)
{

	TComplex formula = QsubMiddle(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubLM_perm2(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubLR_perm2(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubRL_perm3(int n1, int n2)
{

	TComplex formula = QsubRight(n1,1)*QsubLeft(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Two_3SubRM_perm3(int n1, int n2)
{

	TComplex formula = QsubRight(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
 		                 - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::ThreeGap0M(int n1, int n2, int n3)                                                                                                                                 
{

  TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)
                     - QGap0M(n1,1)*QGap0M(n2+n3,2)+2.*QGap0M(n1+n2+n3,3);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::ThreeGap0P(int n1, int n2, int n3)
{

  TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)
                     - QGap0P(n1,1)*QGap0P(n2+n3,2)+2.*QGap0P(n1+n2+n3,3);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Four(int n1, int n2, int n3, int n4)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                 		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                 		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::FourGap0(int n1, int n2, int n3, int n4)
{
	
	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0P(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)
                    -QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3+n4,2)+QGap0P(n1+n2,2)*QGap0M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::FourGap0M(int n1, int n2, int n3, int n4)                                                                                                                          
{

  TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)*QGap0M(n4,1)
                    - QGap0M(n1,1)*QGap0M(n2+n3,2)*QGap0M(n4,1)+2.*QGap0M(n1+n2+n3,3)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n1+n4,2)
                    + QGap0M(n2+n3,2)*QGap0M(n1+n4,2)-QGap0M(n1,1)*QGap0M(n3,1)*QGap0M(n2+n4,2)+QGap0M(n1+n3,2)*QGap0M(n2+n4,2)
                    + 2.*QGap0M(n3,1)*QGap0M(n1+n2+n4,3)-QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3+n4,2)+QGap0M(n1+n2,2)*QGap0M(n3+n4,2)
                    + 2.*QGap0M(n2,1)*QGap0M(n1+n3+n4,3)+2.*QGap0M(n1,1)*QGap0M(n2+n3+n4,3)-6.*QGap0M(n1+n2+n3+n4,4);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::FourGap0P(int n1, int n2, int n3, int n4)
{

  TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)*QGap0P(n4,1)
                    - QGap0P(n1,1)*QGap0P(n2+n3,2)*QGap0P(n4,1)+2.*QGap0P(n1+n2+n3,3)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n1+n4,2)
                    + QGap0P(n2+n3,2)*QGap0P(n1+n4,2)-QGap0P(n1,1)*QGap0P(n3,1)*QGap0P(n2+n4,2)+QGap0P(n1+n3,2)*QGap0P(n2+n4,2)
                    + 2.*QGap0P(n3,1)*QGap0P(n1+n2+n4,3)-QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3+n4,2)+QGap0P(n1+n2,2)*QGap0P(n3+n4,2)
                    + 2.*QGap0P(n2,1)*QGap0P(n1+n3+n4,3)+2.*QGap0P(n1,1)*QGap0P(n2+n3+n4,3)-6.*QGap0P(n1+n2+n3+n4,4);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Four_3SubEvts(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubMiddle(n1,1)*QsubMiddle(n2,1)*QsubLeft(n3,1)*QsubRight(n4,1)-QsubMiddle(n1+n2,2)*QsubLeft(n3,1)*QsubRight(n4,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Four_3SubEvts_perm2(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubLeft(n1,1)*QsubLeft(n2,1)*QsubMiddle(n3,1)*QsubRight(n4,1)-QsubLeft(n1+n2,2)*QsubMiddle(n3,1)*QsubRight(n4,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Four_3SubEvts_perm3(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubRight(n1,1)*QsubRight(n2,1)*QsubMiddle(n3,1)*QsubLeft(n4,1)-QsubRight(n1+n2,2)*QsubMiddle(n3,1)*QsubLeft(n4,1);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
                 - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
                 + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
                 + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
                 + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
                 - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
                 + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
                 - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
                 + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
                 + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
                 - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
                 + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
                 - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
                 - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
                 + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
                 + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
                 + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
                 + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
                 - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
                 + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
                 + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
                 + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
                 + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
                 - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
                 - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
                 - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5); 
	return formula;

}
//___________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{

	
 TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::SixGap0(int n1, int n2, int n3, int n4, int n5, int n6)
{                                                                                                                                                                                                       

  TComplex formula = ThreeGap0M(n1, n2, n3)*ThreeGap0P(n4, n5, n6); 
  return formula;

}
//_________________________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
										Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
												Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
														Narray[int(array[3])], Narray[int(array[4])])*
												Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
															Narray[int(array[3])])*
													Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4
		
		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
													Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3
	
		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Two(Narray[int(array[0])], Narray[int(array[1])])*
													Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
														+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation 
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
				Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

} 
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
										Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
												Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
														Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
												Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
															Narray[int(array[3])], Narray[int(array[4])])*
													Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5
		
		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
													Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4
		
		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
													Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3
	
		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Two(Narray[int(array[0])], Narray[int(array[1])])*
													Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
														+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation 
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
				Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlowGF::EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)                                                                                  
{

  TComplex formula = FourGap0M(n1, n2, n3, n4)*FourGap0P(n5, n6, n7, n8);
  return formula;

}
//_____________________________________________________________________________
void AliAnalysisTaskChargedFlowGF::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
