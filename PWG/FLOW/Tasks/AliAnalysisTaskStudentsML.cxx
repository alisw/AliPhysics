/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

/********************************************** 
* template class for student projects         *
* Marcel Lesch (marcel.lesch@cern.ch)         *
**********************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskStudentsML.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "TRandom.h" 
#include "TComplex.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include <TArrayD.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliHeader.h"
#include "TExMap.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStudentsML)

//================================================================================================================

AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 // Control histograms:
 fCentralityHistogramBefore(NULL),
 //SelectionCuts
 bDoAnalysis(kTRUE),
 bUseRecoKineTable(kTRUE),
 bUseWeakSecondaries(kTRUE),
 bSaveAllQA(kTRUE),
 bMultCut(kFALSE),
 fMainFilter(0),
 fSecondFilter(0),
 fSlopeUpperLine(0.),
 fAxisUpperLine(0.),
 fSlopeLowerLine(0.),
 fAxisLowerLine(0.),
 fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
 fCentralityBins(16),
 //Global Quality
 bCutOnVertexX(kTRUE), 
 bCutOnVertexY(kTRUE),
 bCutOnVertexZ(kTRUE), 
 fMinVertexX(-44.),
 fMaxVertexX(-44.),
 fMinVertexY(-44.),
 fMaxVertexY(-44.),
 fMinVertexZ(-10.),
 fMaxVertexZ(10.),
 fCentralityEstimator("CL1"), 
 //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 bNumberTPCCluster(kTRUE),
 bNumberITSCluster(kTRUE),
 bChiSquareTPC(kTRUE),
 bDCAz(kTRUE),
 bDCAxy(kTRUE),
 bChargeCut(kFALSE),		
 bChargePos(kTRUE), 
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fMinTPCCluster(70.),
 fMinITSCluster(2.),
 fChooseChiSquareMethod(1), 
 fMinChiSquareTPC(0.1),
 fMaxChiSquareTPC(4.0),
 fMaxDCAz(3.2),
 fMaxDCAxy(2.4),
 bDoFisherYates(kFALSE),
 fFisherYatesCutOff(1.),
 //Weights
 bUseWeights(kFALSE),
 bUsePtWeights(kFALSE),
 bUsePhiWeights(kFALSE), 
 bUseEtaWeights(kFALSE),
 fNumberRuns(90),
 //Variables for the correlation
 fMaxCorrelator(14),
 fNumber(0),  		//number of correlation first correlator
 fNumberSecond(0), 	//number of correlation second correlator
 fNumberThird(0),	//number of correlation second correlator
 fNumberFourth(0),	//number of correlation fourth correlator
 fNumberFifth(0),	//number of correlation fifth correlator
 fNumberSixth(0),	//number of correlation sixth correlator
 fNumberSeventh(0),	//number of correlation seventh correlator
 fNumberEighth(0),	//number of correlation eighth correlator
 fMinNumberPart(14),
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), //first set of harmonics
 fb1(0), fb2(0), fb3(0), fb4(0), fb5(0), fb6(0), fb7(0), //second set of harmonics
 fd1(0), fd2(0), fd3(0), fd4(0), fd5(0), fd6(0), fd7(0), //third set of harmonics
 fe1(0), fe2(0), fe3(0), fe4(0), fe5(0), fe6(0), fe7(0), //fourth set of harmonics
 ff1(0), ff2(0), ff3(0), ff4(0), ff5(0), ff6(0), ff7(0),  //fifth set of harmonics 
 fg1(0), fg2(0), fg3(0), fg4(0), fg5(0), fg6(0), fg7(0),  // sixth set of harmonics 
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), //seventh set of harmonics
 fi1(0), fi2(0), fi3(0), fi4(0), fi5(0), fi6(0), fi7(0),  //eighth set of harmonics 
 bDoMixed(kFALSE),
 bDifferentCharge(kTRUE),
 bSetSameChargePositiv(kTRUE),			
 fMixedHarmonic(0),
 fCounterHistogram(NULL),
 fProfileEventCuts(NULL),  	
 fProfileTrackCuts(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("outputStudentAnalysis");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reserved              
  // Output slot #1 writes into a TList container

  DefineOutput(1, TList::Class());  

  if(useParticleWeights)
  {
   // not needed for the time being
  }

} // AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights): 

//==========================================================================================================================================================================

AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 // Control histograms:
 fCentralityHistogramBefore(NULL),
 //SelectionCuts
 bDoAnalysis(kTRUE),
 bUseRecoKineTable(kTRUE),
 bUseWeakSecondaries(kTRUE),
 bSaveAllQA(kTRUE),
 bMultCut(kFALSE),
 fMainFilter(0),
 fSecondFilter(0),
 fSlopeUpperLine(0.),
 fAxisUpperLine(0.),
 fSlopeLowerLine(0.),
 fAxisLowerLine(0.),
 fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
 fCentralityBins(16),
 //Global Quality
 bCutOnVertexX(kTRUE), 
 bCutOnVertexY(kTRUE),
 bCutOnVertexZ(kTRUE), 
 fMinVertexX(-44.),
 fMaxVertexX(-44.),
 fMinVertexY(-44.),
 fMaxVertexY(-44.),
 fMinVertexZ(-10.),
 fMaxVertexZ(10.),
 fCentralityEstimator("CL1"),
 //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 bNumberTPCCluster(kTRUE),
 bNumberITSCluster(kTRUE),
 bChiSquareTPC(kTRUE),
 bDCAz(kTRUE),
 bDCAxy(kTRUE),
 bChargeCut(kFALSE),		
 bChargePos(kTRUE), 
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fMinTPCCluster(70.),
 fMinITSCluster(2.),
 fChooseChiSquareMethod(1), 
 fMinChiSquareTPC(0.1),
 fMaxChiSquareTPC(4.0),
 fMaxDCAz(3.2),
 fMaxDCAxy(2.4),
 bDoFisherYates(kFALSE),
 fFisherYatesCutOff(1.),
 //Weights
 bUseWeights(kFALSE),
 bUsePtWeights(kFALSE),
 bUsePhiWeights(kFALSE), 
 bUseEtaWeights(kFALSE),
 fNumberRuns(90),
 //Variables for the correlation
 fMaxCorrelator(14),
 fNumber(0),  		//number of correlation first correlator
 fNumberSecond(0), 	//number of correlation second correlator
 fNumberThird(0),	//number of correlation second correlator
 fNumberFourth(0),	//number of correlation fourth correlator
 fNumberFifth(0),	//number of correlation fifth correlator
 fNumberSixth(0),	//number of correlation sixth correlator
 fNumberSeventh(0),	//number of correlation seventh correlator
 fNumberEighth(0),	//number of correlation eighth correlator
 fMinNumberPart(14),
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), //first set of harmonics
 fb1(0), fb2(0), fb3(0), fb4(0), fb5(0), fb6(0), fb7(0), //second set of harmonics
 fd1(0), fd2(0), fd3(0), fd4(0), fd5(0), fd6(0), fd7(0), //third set of harmonics
 fe1(0), fe2(0), fe3(0), fe4(0), fe5(0), fe6(0), fe7(0), //fourth set of harmonics
 ff1(0), ff2(0), ff3(0), ff4(0), ff5(0), ff6(0), ff7(0),  //fifth set of harmonics 
 fg1(0), fg2(0), fg3(0), fg4(0), fg5(0), fg6(0), fg7(0),  // sixth set of harmonics 
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), //seventh set of harmonics
 fi1(0), fi2(0), fi3(0), fi4(0), fi5(0), fi6(0), fi7(0),  //eighth set of harmonics 
 // Final results:
 bDoMixed(kFALSE),
 bDifferentCharge(kTRUE),
 bSetSameChargePositiv(kTRUE),		
 fMixedHarmonic(0),
 fCounterHistogram(NULL),
 fProfileEventCuts(NULL),  	
 fProfileTrackCuts(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML()");

  this->InitializeArrays();

} // AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML():

//==========================================================================================================================================================================

AliAnalysisTaskStudentsML::~AliAnalysisTaskStudentsML()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskStudentsML::~AliAnalysisTaskStudentsML()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
 // d) 
  
 // a) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 if(bSaveAllQA){this->BookControlHistograms();}
 this->BookFinalResultsHistograms();

 // d) Save cut values 
 if (fCentralityEstimator == "CL1") {fProfileEventCuts->Fill(0.5, 1);} 
 else if ((fCentralityEstimator == "V0M")) {fProfileEventCuts->Fill(0.5, 2);} // 1: CL1, 2: V0M 
 //if (fSndCentralityEstim == "CL1") {fProfileEventCuts->Fill(1.5, 1);} 
 //else {fProfileEventCuts->Fill(1.5, 2);}  
 if (bCutOnVertexX) {fProfileEventCuts->Fill(2.5, fMinVertexX); fProfileEventCuts->Fill(3.5, fMaxVertexX);} 
 if (bCutOnVertexY) {fProfileEventCuts->Fill(4.5, fMinVertexY); fProfileEventCuts->Fill(5.5, fMaxVertexY);} 
 if (bCutOnVertexZ) {fProfileEventCuts->Fill(6.5, fMinVertexZ); fProfileEventCuts->Fill(7.5, fMaxVertexZ);} 
 fProfileEventCuts->Fill(8.5, fMinNumberPart); 
 fProfileEventCuts->Fill(9.5, fMainFilter); 
 fProfileEventCuts->Fill(10.5, fSecondFilter); 
 if (bMultCut) 
 {
   fProfileEventCuts->Fill(11.5, fSlopeLowerLine); 
   fProfileEventCuts->Fill(12.5, fAxisLowerLine); 
   fProfileEventCuts->Fill(13.5, fSlopeUpperLine); 
   fProfileEventCuts->Fill(14.5, fAxisUpperLine); 
 }
 if(bUseRecoKineTable){ fProfileEventCuts->Fill(15.5, 1); } 




 if (bCutOnPt) {fProfileTrackCuts->Fill(0.5, fMinPtCut); fProfileTrackCuts->Fill(1.5, fMaxPtCut);} 
 if (bCutOnEta) {fProfileTrackCuts->Fill(2.5, fMinEtaCut); fProfileTrackCuts->Fill(3.5, fMaxEtaCut);} 
 if (bNumberTPCCluster) {fProfileTrackCuts->Fill(4.5, fMinTPCCluster);} 
 if (bChiSquareTPC) {fProfileTrackCuts->Fill(5.5, fMinChiSquareTPC); fProfileTrackCuts->Fill(6.5, fMaxChiSquareTPC); fProfileTrackCuts->Fill(7.5, fChooseChiSquareMethod);} 
 if (bNumberITSCluster) {fProfileTrackCuts->Fill(8.5, fMinITSCluster);} 
 if (bDCAxy) {fProfileTrackCuts->Fill(9.5, fMaxDCAxy);} 
 if (bDCAz) {fProfileTrackCuts->Fill(10.5, fMaxDCAz);} 
 if (bChargeCut)
 {     
    fProfileTrackCuts->Fill(11.5, 1);
    if(bChargePos) { fProfileTrackCuts->Fill(12.5, 1); } 
    else{ fProfileTrackCuts->Fill(12.5, -1); }

 }
 if(bDoFisherYates){ fProfileTrackCuts->Fill(13.5, 1); fProfileTrackCuts->Fill(14.5, fFisherYatesCutOff); } 
 if (bUsePtWeights) {fProfileTrackCuts->Fill(15.5, 1);} 
 if (bUsePhiWeights) {fProfileTrackCuts->Fill(16.5, 1);} 
 if (bUseEtaWeights) {fProfileTrackCuts->Fill(17.5, 1);} 

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsML::UserCreateOutputObjects() 

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a) Get pointer to AOD event and MCEvent;
 // b) Do Analysis or Weights
 // c) PostData.

 fCounterHistogram->Fill(0.5); // Checks if User Exec is entered proberly

 //a) Get pointer to AOD event and MCEvent;

 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent());
 AliMCEvent *aMCEvent = MCEvent();

 //b) Do Analysis or Weights
 if(bDoAnalysis){PhysicsAnalysis(aAOD);}
 else{GetKineDist(aAOD,aMCEvent);}

 // c) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsML::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsML::PhysicsAnalysis(AliAODEvent *aAODEvent)
{

 // a) Global QA (Multiplicity and Vertex cut) 
 // b.0) Start analysis over AODs;
 // 	b.1) Frist Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut) 
 //	     -> Starting first loop over tracks: Used to fill Histograms + Determination of Multiplicity (after Trackselection)
 // 	b.2) Second Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut) 
 //	     -> Fill the Phi Arrays
 // c) Calculus
 // d) Reset event-by-event objects;
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //a) Global QA (Centrality check, Vertex cut and high multiplicity outlier)
 if(!GlobalQualityAssurance(aAODEvent)){return;}

 AliAODVertex *avtx = (AliAODVertex*)aAODEvent->GetPrimaryVertex();

 //Select Centrality Bin (has to be valid, as it passed the GlobalQualityAssurance)
 Int_t CentralityBin = SelectCentrality(aAODEvent);
	
 if(bSaveAllQA) 
 {
 	fVertexXHistogram[CentralityBin][1]->Fill(avtx->GetX());
 	fVertexYHistogram[CentralityBin][1]->Fill(avtx->GetY());
 	fVertexZHistogram[CentralityBin][1]->Fill(avtx->GetZ());
 } //if(bSaveAllQA)

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 //Access Data

 //b.0) Start analysis over AODs
 Int_t nTracks = 0;  // number of all tracks in current event 
 nTracks = aAODEvent->GetNumberOfTracks(); 

 Double_t* angles_A = NULL; // Azimuthal angles
 Double_t* pt_A = NULL;
 Double_t* eta_A = NULL;
 Double_t* weights_A = NULL;
 Int_t Mult_A = 0.;

 Bool_t *PassedTrackSelection = new Bool_t[nTracks]; //Holds information if a track passed the track selection. kTRUE if yes, kFALSE otherwise

 Double_t* angles_B = NULL; 
 Int_t Mult_B = 0.; 
 Int_t CounterSameCharge = 0.; //used if bDifferentCharge = kTRUE


 if(bSaveAllQA)
 { 
	if(nTracks>0){fMultHistogram[CentralityBin][1]->Fill(nTracks);} //multiplicity distribution before track selection
 } //if(bSaveAllQA)


 AliAODVertex *primaryVertex = (AliAODVertex*)aAODEvent->GetPrimaryVertex(); 

 //b.1) Frist Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut)
 //	Starting first loop over tracks: Used to fill Histograms + Determination of Multiplicity (after Trackselection)
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) 
 {
	AliAODTrack *aAODTrack = NULL;
    	aAODTrack = dynamic_cast<AliAODTrack*>(aAODEvent->GetTrack(iTrack));
    
	if(!aAODTrack){continue;} // protection against NULL pointers
  	if(!aAODTrack->TestFilterBit(fMainFilter)){continue;} //Check if in Filter
    
    	Double_t phi = aAODTrack->Phi(); // azimuthal angle
     	Double_t eta = aAODTrack->Eta(); // pseudorapidity
     	Double_t pt = aAODTrack->Pt(); // Pt (transverse momentum)
     	Double_t charge = aAODTrack->Charge(); // charge
     	Int_t NumberOfTPCClusters = aAODTrack->GetTPCNcls(); //number of TPC clusters of the track
     	Int_t NumberOfITSClusters = aAODTrack->GetITSNcls(); //number of ITS clusters of the track


	//Get the DCA information (cf PWGCF/EBYE/BalanceFunctions/AliAnalysisTaskBFPsi.cxx)
	Float_t ValueDCAxy = 999.;   // DCA in the xy-plane.
	Float_t ValueDCAz = 999.;    // DCA along z.

	if (fMainFilter == 128)  // These methods work only for constrained TPConly tracks.
	{ //These two quantities are the DCA from global tracks but not what we will cut on.
	  ValueDCAxy = aAODTrack->DCA();
	  ValueDCAz = aAODTrack->ZAtDCA();
	}
	else  //For the unconstrained tracks.
	{
	  Double_t v[3];    //Coordinates of the PV?
	  Double_t pos[3];  //Coordinates of the track closest to PV?

          primaryVertex->GetXYZ(v);
	  aAODTrack->GetXYZ(pos);
	  ValueDCAxy = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
	  ValueDCAz = pos[2] - v[2];
	}

	//Get the chi^2 per TPC cluster.
	Float_t ChiSquareInTPC = 999.;
	/// Personal method, should be equal to GetTPCchi2perCluster()
	if (fChooseChiSquareMethod == 1) {ChiSquareInTPC = (aAODTrack->GetTPCchi2())/(aAODTrack->GetNcls(1));}
	else if (fChooseChiSquareMethod == 2) {ChiSquareInTPC = aAODTrack->GetTPCchi2perCluster();}
	else if (fChooseChiSquareMethod == 3) {ChiSquareInTPC = aAODTrack->GetTPCchi2perNDF();}
	else {ChiSquareInTPC = aAODTrack->Chi2perNDF();}

  	//............................................................................................
	//Fill control histograms with the particles before track selection:
     	if(bSaveAllQA) 
     	{
		if(!bDoMixed)
		{
		  fPhiHistogram[CentralityBin][0]->Fill(phi); 
		  fEtaHistogram[CentralityBin][0]->Fill(eta);
		  fPTHistogram[CentralityBin][0]->Fill(pt);
		}//if(!bDoMixed)

		if(bDoMixed)
		{
		  if(charge>0.)
		  {
		    fPhiHistogram[CentralityBin][0]->Fill(phi); 
		    fEtaHistogram[CentralityBin][0]->Fill(eta);
		    fPTHistogram[CentralityBin][0]->Fill(pt);
		  }//if(charge>0.)

		  if(charge<0.)
		  {
		    fPhiHistogram[CentralityBin][2]->Fill(phi); 
		    fEtaHistogram[CentralityBin][2]->Fill(eta);
		    fPTHistogram[CentralityBin][2]->Fill(pt);
		  } //if(charge<0.)
		} //if(bDoMixed)
   	 

     		fTPCClustersHistogram[CentralityBin][0]->Fill(NumberOfTPCClusters);
     		fITSClustersHistogram[CentralityBin][0]->Fill(NumberOfITSClusters);
    		fChiSquareTPCHistogram[CentralityBin][0]->Fill(ChiSquareInTPC);
     		fDCAzHistogram[CentralityBin][0]->Fill(ValueDCAz);
		fDCAxyHistogram[CentralityBin][0]->Fill(ValueDCAxy);
		fChargeHistogram[CentralityBin][0]->Fill(charge); 

        } //if(bSaveAllQA)

	//............................................................................................
	//Check if track passes track selection. If no, continue
      	if(!TrackSelection(primaryVertex, aAODTrack)){ PassedTrackSelection[iTrack]=kFALSE; continue;} 

	PassedTrackSelection[iTrack]=kTRUE; 

	//............................................................................................
     	// Fill control histograms with the particles after track selection:
     	if(bSaveAllQA) 
     	{
		if(!bDoMixed)
		{
		fPhiHistogram[CentralityBin][1]->Fill(phi); 
		fEtaHistogram[CentralityBin][1]->Fill(eta);
		fPTHistogram[CentralityBin][1]->Fill(pt);
		}//if(!bDoMixed)

		if(bDoMixed)
		{
		  if(charge>0.)
		  {
		    fPhiHistogram[CentralityBin][1]->Fill(phi); 
		    fEtaHistogram[CentralityBin][1]->Fill(eta);
		    fPTHistogram[CentralityBin][1]->Fill(pt);
		  }//if(charge>0.)
		  if(charge<0.)
		  {
		    fPhiHistogram[CentralityBin][3]->Fill(phi); 
		    fEtaHistogram[CentralityBin][3]->Fill(eta);
		    fPTHistogram[CentralityBin][3]->Fill(pt);
		  }//if(charge<0.)
		}//if(bDoMixed)

	     	fTPCClustersHistogram[CentralityBin][1]->Fill(NumberOfTPCClusters);
	     	fITSClustersHistogram[CentralityBin][1]->Fill(NumberOfITSClusters);
	     	fChiSquareTPCHistogram[CentralityBin][1]->Fill(ChiSquareInTPC);
	     	fDCAzHistogram[CentralityBin][1]->Fill(ValueDCAz);
	     	fDCAxyHistogram[CentralityBin][1]->Fill(ValueDCAxy);
		fChargeHistogram[CentralityBin][1]->Fill(charge); 

        } //if(bSaveAllQA)

	//............................................................................................

     	if(!bDoMixed) {Mult_A+=1.;}
    	if(bDoMixed)
     	{
		if(bDifferentCharge)
		{
			if(charge>0.){Mult_A+=1.;}
			if(charge<0.){Mult_B+=1.;}
		}//if(bDifferentCharge)

		if(!bDifferentCharge)
		{
			Double_t UsedCharge = 0.;
			if(bSetSameChargePositiv) {UsedCharge = charge;}
			if(!bSetSameChargePositiv) {UsedCharge = -1.*charge;}

			if(UsedCharge>0.)
			{
		  	   if(0 == CounterSameCharge%2) {Mult_A+=1.; CounterSameCharge++; }
		  	   else {Mult_B+=1.; CounterSameCharge++; }
			}

		}//if(!bDifferentCharge)
     	}//if(bDoMixed)


 } // First Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut)



 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Fisher Yales

 Int_t Before_FisherYates_Mult = Mult_A;  //Multiplicity directly after track selection, before cut off 
 Int_t After_FisherYates_Mult=0;
 Int_t *RandomIndexArray = new Int_t[Before_FisherYates_Mult]; //Array to hold array positions

 if(bDoFisherYates)
 { 
	this->FisherYatesRandomizing(Before_FisherYates_Mult, RandomIndexArray); //Random Index Generator 
	Mult_A = (Int_t)(fFisherYatesCutOff*(Float_t)Mult_A); //Cut of the Multiplicity
	After_FisherYates_Mult = Mult_A;
 } 


 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Checking if required minimum particles per event is reached and filling multiplicity Distributions after Track selection

 if(!bDoMixed)
 {
   if(bSaveAllQA)
   {
	if(Mult_A>0){fMultHistogram[CentralityBin][2]->Fill(Mult_A);} //multiplicity distribution after track selection 
   }//if(bSaveAllQA)

   if(Mult_A< fMinNumberPart) { return; } 

 }//if(!bDoMixed)
 
 if(bDoMixed)
 {
   if(bSaveAllQA)
   {
  	if(Mult_A>0){fMultHistogram[CentralityBin][2]->Fill(Mult_A);} //multiplicity distribution after track selection
   	if(Mult_B>0){fMultHistogram[CentralityBin][3]->Fill(Mult_B);}
   }//if(bSaveAllQA)

   if(Mult_A < fMinNumberPart || Mult_B < fMinNumberPart ){ return; }

 }//if(bDoMixed)



 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Declare arrays for Phi, pt, eta and weights

 angles_A = new Double_t[Mult_A]; 
 pt_A = new Double_t[Mult_A];
 eta_A = new Double_t[Mult_A];
 weights_A = new Double_t[Mult_A];

 Mult_A = 0.; 
       
 if(Mult_B>0){angles_B = new Double_t[Mult_B];}
 else{angles_B = new Double_t[1]; angles_B[0]=0.;} //Dummy 
 Mult_B = 0.; //Reset for filling

 CounterSameCharge = 0.; //reset the same charge counter


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//b.2) Second Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut) 
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {  
    AliAODTrack *aAODTrack = NULL;
    aAODTrack = dynamic_cast<AliAODTrack*>(aAODEvent->GetTrack(iTrack));
    

    Double_t phi = 0.; // azimuthal angle
    Double_t pt = 0.; // Pt (transverse momentum)
    Double_t eta = 0.;
    Double_t charge = 0.; // charge


    if(!aAODTrack){continue;} // protection against NULL pointers
    if(!aAODTrack->TestFilterBit(fMainFilter)){continue;} //Check if in Filter

    phi = aAODTrack->Phi(); // azimuthal angle
    pt = aAODTrack->Pt(); // Pt (transverse momentum)
    eta = aAODTrack->Eta();
    charge = aAODTrack->Charge(); // charge
    
     if(!PassedTrackSelection[iTrack]){continue;} //Track did not pass physics selection 
     

     if(!bDoMixed)
     {	
	if(bDoFisherYates)
	{
		Int_t NewIndex = RandomIndexArray[Mult_A];
		if(NewIndex>=After_FisherYates_Mult) continue; //this is outside our Array Size
		
		angles_A[NewIndex] = phi; 
		pt_A[NewIndex] = pt; 
		eta_A[NewIndex] = eta; 
		weights_A[NewIndex] = 1.;
	}
	else
	{
		angles_A[Mult_A] = phi; 
		pt_A[Mult_A] = pt; 
		eta_A[Mult_A] = eta; 
		weights_A[Mult_A] = 1.;
	} 

        Mult_A+=1.;

     }//if(!bDoMixed)

     if(bDoMixed)
     {
	if(bDifferentCharge)
	{	
		if(charge>0.){angles_A[Mult_A] = phi; Mult_A+=1.;}
		if(charge<0.){angles_B[Mult_B] = phi; Mult_B+=1.;}
	}//if(bDifferentCharge)

	if(!bDifferentCharge)
	{
		Double_t UsedCharge = 0.;
		if(bSetSameChargePositiv) {UsedCharge = charge;}
		if(!bSetSameChargePositiv) {UsedCharge = -1.*charge;}

		if(UsedCharge>0.)
		{
		   if(0 == CounterSameCharge%2) {angles_A[Mult_A] = phi; Mult_A+=1.; CounterSameCharge++; }
		   else {angles_B[Mult_B] = phi; Mult_B+=1.; CounterSameCharge++;}
		}

	}//if(!bDifferentCharge)

     }//if(bDoMixed)

 } //Second Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut) 


 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Add Member Function Get Weights
 if(bUseWeights)
 {
  Int_t CallRunNumber = aAODEvent->GetRunNumber();
  CalculateWeight(CentralityBin, CallRunNumber, weights_A, Mult_A, angles_A, pt_A, eta_A);
 }
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //c) Calculus
 
 if(!bDoMixed)
 {
   this->MainTask(CentralityBin, Mult_A, angles_A,weights_A); //Actual Multi-Particle Correlation
 }//if(!bDoMixed)
 
 if(bDoMixed)
 {
   this->MixedParticle(CentralityBin, fMixedHarmonic, Mult_A, angles_A, Mult_B, angles_B);
 }//if(bDoMixed)
 

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // d) Reset event-by-event objects:
 delete [] PassedTrackSelection;
 Mult_A =0.;
 Before_FisherYates_Mult=0.;
 After_FisherYates_Mult=0.;
 delete [] RandomIndexArray;
 delete [] angles_A; 
 delete [] pt_A;
 delete [] eta_A;
 delete [] weights_A;
 Mult_B =0.;
 delete [] angles_B;
 CounterSameCharge = 0.; //reset the same charge counter

} //void AliAnalysisTaskStudentsML::PhysicsAnalysis()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::GetKineDist(AliAODEvent *aAODEve, AliMCEvent *aMCEve)
{ 

 TString sMethod = "void AliAnalysisTaskStudentsML::GetKineDist(AliAODEvent *aAODEve, AliMCEvent *aMCEve)";

 //Used for gaining weights
 //a) Global QA for Reco and Kine
 //b) Run over Reco + Generate Table 
 //c) Run over Kine + Get Distributions

 //a) Global QA for Reco and Kine (for Reco only if we use the RecoKineTable)
 if(bUseRecoKineTable){ if(!GlobalQualityAssurance(aAODEve)){return;} }

 //Get Centrality from AOD event
 if(!aAODEve){return;} //In case we did not apply the global QA on the AOD, protect against NULL pointer
 
 Int_t CentralityBin = SelectCentrality(aAODEve);

 if(!bUseRecoKineTable) //In case we did not apply the global QA on the AOD, do the centrality selection 
 { 
   AliMultSelection *ams = (AliMultSelection*)aAODEve->FindListObject("MultSelection");
   if(!ams){return;}

   Float_t CentralityValue = 0.;

   if(fCentralityEstimator == "V0M" ) { CentralityValue = ams->GetMultiplicityPercentile("V0M"); }
   else if(fCentralityEstimator == "CL1") { CentralityValue = ams->GetMultiplicityPercentile("CL1"); }
   else {Fatal(sMethod.Data(), "FATAL: no valid centrality estimator!");} 

   fCentralityHistogramBefore->Fill(CentralityValue);

   if(CentralityBin<0){return;} //No valid centrality bin

   if(bSaveAllQA){fCentralityHistogram[CentralityBin]->Fill(CentralityValue);}
 }

 AliAODVertex *primaryVertex = (AliAODVertex*)aAODEve->GetPrimaryVertex(); 


 if(!GlobalQualityAssurance(CentralityBin, aMCEve)){return;} 

 AliMCVertex *avtx = (AliMCVertex*)aMCEve->GetPrimaryVertex();

 if(bSaveAllQA)
 {
	 fVertexXHistogram[CentralityBin][1]->Fill(avtx->GetX()); 
	 fVertexYHistogram[CentralityBin][1]->Fill(avtx->GetY()); 
	 fVertexZHistogram[CentralityBin][1]->Fill(avtx->GetZ()); 
 }//if(bSaveAllQA)

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //b)Run over Reco + Generate Table 
 Int_t nTracks = 0;
 nTracks = aAODEve->GetNumberOfTracks();

 TExMap *RecoKineMap = new TExMap();

 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) 
 {
	AliAODTrack *aAODTrack = NULL;
    	aAODTrack = dynamic_cast<AliAODTrack*>(aAODEve->GetTrack(iTrack));

	if(!aAODTrack){continue;} // protection against NULL pointers

  	if(!aAODTrack->TestFilterBit(fMainFilter)){continue;} //Check if in Filter
   
	//........................................................................
	//Track did not pass physics selection 
      	if(!TrackSelection(primaryVertex, aAODTrack))
	{
          // The reco track did not pass the selection: MC-key tagged with "-1" in the table.
          if (aAODTrack->GetLabel() >= 0)
          {
            if ( RecoKineMap->GetValue(aAODTrack->GetLabel()) == 0 ) {RecoKineMap->Add(aAODTrack->GetLabel(), -1);}
          }

	} //if(!TrackSelection(aAODTrack)) 
        else
        {
          // The reco track passed the selection: MC-key tagged with "1" in the table.
          if (aAODTrack->GetLabel() >= 0)
          {
            if ( RecoKineMap->GetValue(aAODTrack->GetLabel()) == 0 ) {RecoKineMap->Add(aAODTrack->GetLabel(), 1);}
          }
        }

 } // Loop over the Reco tracks to Get Reco-Kine Table

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //c) Run over Kine + Get Distributions
 nTracks = 0;
 nTracks = aMCEve->GetNumberOfTracks();

 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) 
 {
	AliAODMCParticle *aMCTrack = NULL;
    	aMCTrack = dynamic_cast<AliAODMCParticle*>(aMCEve->GetTrack(iTrack));

	if(!aMCTrack){continue;} // protection against NULL pointers
   
	Double_t phi = aMCTrack->Phi(); // azimuthal angle
     	Double_t eta = aMCTrack->Eta(); // pseudorapidity
     	Double_t pt = aMCTrack->Pt(); // Pt (transverse momentum)
     	Double_t charge = aMCTrack->Charge(); // charge
	Bool_t isPrimary = aMCTrack->IsPhysicalPrimary(); // kTRUE: the particle is a primary particle.
	Bool_t isWeakSecondary = aMCTrack->IsSecondaryFromWeakDecay(); 
	//Bool_t isMaterialSecondary = aMCTrack->IsSecondaryFromMaterial();

	if(charge == 0){continue;}

	if (bUseWeakSecondaries)   // Primaries and weak secondaries can be selected.
  	{ if (!isPrimary && !isWeakSecondary) { continue;} }
  	else { if (!isPrimary) { continue;} } // Only the primaries can be selected.


 	if(bSaveAllQA)
 	{
		fPhiHistogram[CentralityBin][0]->Fill(phi); 
		fEtaHistogram[CentralityBin][0]->Fill(eta); 
		fPTHistogram[CentralityBin][0]->Fill(pt); 
	}//if(bSaveAllQA)
	//........................................................................
	//Track did not pass physics selection 
      	if(!TrackSelection(aMCTrack)){continue;} 

	if(bUseRecoKineTable)
	{

		if(RecoKineMap->GetValue(iTrack)==-1){continue;} // GetValue = -1: the corresponding reco track has not been selected by its track selection.
		else
		{	// GetValue = 0: the kine track has been lost in ALICE --> Included in the pT distribution.
        		// GetValue = 1: the reco track has been selected.
 			if(bSaveAllQA)
 			{
				fPhiHistogram[CentralityBin][1]->Fill(phi); 
				fEtaHistogram[CentralityBin][1]->Fill(eta); 
				fPTHistogram[CentralityBin][1]->Fill(pt); 
			}//if(bSaveAllQA)
		} //else
	} //if(bUseRecoKineTable)
	else
	{
 		if(bSaveAllQA)
 		{
	 		fPhiHistogram[CentralityBin][1]->Fill(phi); 
			fEtaHistogram[CentralityBin][1]->Fill(eta); 
			fPTHistogram[CentralityBin][1]->Fill(pt); 
		}//if(bSaveAllQA)
	}

 } // Loop over Kine Tracks with Histogram Filling

 delete RecoKineMap;

} // AliAnalysisTaskStudentsML::GetKineDist(AliAODEvent *aAODEve, AliMCEvent *aMCEve)

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::Terminate(Option_t *)
{
 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:
 // ...


} // end of void AliAnalysisTaskStudentsML::Terminate(Option_t *)

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::InitializeArrays()
{
  // Initialize all data members which are arrays in this method.

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Multiparticle-Correlations

   for(Int_t js=0;js<113;js++) 
   {
     for(Int_t j=0;j<15;j++)
     {
      fQvector[js][j] = TComplex(0.,0.); //! 
     } 
   } 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Run-by-Run List and Weight-Histograms
  for (Int_t iRun = 0; iRun < 90; iRun++)
  {
	 fListRuns[iRun] = 0;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Centrality Dependend Objects: Lists, QC-Histograms and Output Histograms

  for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
  {


	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // Run-by-Run List and Weight-Histograms
	  for (Int_t iRun = 0; iRun < 90; iRun++)
	  {
	    fHistoPtWeight[icent][iRun] = NULL;
	    fHistoEtaWeight[icent][iRun] = NULL;
	    fHistoPhiWeight[icent][iRun] = NULL;
	  }

	//Lists
	fCentralityList[icent] = NULL;
	fControlHistogramsList[icent] = NULL;
	fFinalResultsList[icent] = NULL;

	//QC-Histograms
	 for(Int_t i=0; i<2; i++)
	 {
	   fVertexXHistogram[icent][i] = NULL;
	   fVertexYHistogram[icent][i] = NULL;
	   fVertexZHistogram[icent][i] = NULL;
	   fTPCClustersHistogram[icent][i] = NULL;
	   fITSClustersHistogram[icent][i] = NULL;
	   fChiSquareTPCHistogram[icent][i] = NULL;
	   fDCAzHistogram[icent][i] = NULL;
	   fDCAxyHistogram[icent][i] = NULL;
           fChargeHistogram[icent][i] = NULL; 
	   
	 }

	 fCentralityHistogram[icent] = NULL;

	 for(Int_t i=0; i<5; i++)
	 {
	   fPTHistogram[icent][i] = NULL;
	   fPhiHistogram[icent][i] = NULL;
	   fEtaHistogram[icent][i] = NULL;
	 }   

	 for(Int_t i=0; i<4; i++)
	 {
	   fMultHistogram[icent][i] = NULL;
	 }

	//Output Histograms
	fResults[icent] = NULL;
	fCovResults[icent] = NULL; 
	fJoinedCovResults[icent] = NULL; 
  	fMixedParticleHarmonics[icent] = NULL;

  }

} // void AliAnalysisTaskStudentsML::InitializeArrays()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::SetInitializeCentralityArray()
{

 TString sMethodName = "void AliAnalysisTaskStudentsML::BookAndNestAllLists()";

  Float_t ListCentralities[17] = { fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16 };

  for(Int_t i=0; i<17; i++) { fcentralityArray[i] = ListCentralities[i]; }

  //Protections
  if(fcentralityArray[0] < 0 || fcentralityArray[1] < 0) { Fatal(sMethodName.Data(),"First Centrality bin not defined"); } //The need at least one well defined centrality bin

  for(Int_t icent=0; icent<fCentralityBins; icent++)
  {
	//The next bin should be a valid boundery, i.e. it is > 0. but it is also smaller than the previous boundery -> Wrong ordering
	if( fcentralityArray[icent+1] > 0. && fcentralityArray[icent+1] < fcentralityArray[icent] ) { Fatal(sMethodName.Data(),"Wrong ordering of centrality bounderies"); }
  }

} //void AliAnalysisTaskStudentsML::InitializeCentralityArray()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskStudentsML::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
 {
	//Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	fCentralityList[icent] = new TList();
	fCentralityList[icent] ->SetName(Form("MultCut_task=>%.1f-%.1f",fcentralityArray[icent],fcentralityArray[icent+1]));
	fCentralityList[icent]->SetOwner(kTRUE);
	fHistList->Add(fCentralityList[icent]);

	 
	
	// a) Book and nest lists for control histograms:
	fControlHistogramsList[icent] = new TList();
	fControlHistogramsList[icent]->SetName("ControlHistograms");
	fControlHistogramsList[icent]->SetOwner(kTRUE);
	if(bSaveAllQA){	fCentralityList[icent]->Add(fControlHistogramsList[icent]); }
	

	// b) Book and nest lists for final results:
	fFinalResultsList[icent] = new TList();
	fFinalResultsList[icent]->SetName("FinalResults");
	fFinalResultsList[icent]->SetOwner(kTRUE);
	fCentralityList[icent]->Add(fFinalResultsList[icent]);

 }


} // void AliAnalysisTaskStudentsML::BookAndNestAllLists()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra
 // b) Book histogram to hold phi spectra
 // c) Book histogram to hold eta spectra
 // d) Book histogam to hold multiplicty distributions 
 // e) Book histogam for Vertex X 
 // f) Book histogam for Vertex Y 
 // g) Book histogam for Vertex Z 
 // h) Book histogram to debug
 // i) Book histogram for number of TPC clustes 
 // j) Book histogram for number of ITC clusters 
 // k) Book histogram for chi square TPC 
 // l) Book histogram for DCAz 
 // m) Book histogram for DCAxy 
 // n) Book histogram for Charge Cut 
 // o) Book histogram Centrality 

 for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
 {

	//Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	 // a) Book histogram to hold pt spectra:
	 fPTHistogram[icent][0] = new TH1F("fPTHistBeforeTrackSelection","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][0]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][0]);

	 fPTHistogram[icent][1] = new TH1F("fPTHistAfterTrackSelection","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][1]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][1]);

	 fPTHistogram[icent][4] = new TH1F("fPTHistAfterTrackSelectionWeighted","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][4]->Sumw2(); 
	 fPTHistogram[icent][4]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][4]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][4]);

	 fPTHistogram[icent][2] = new TH1F("fPTHistBeforeTrackSelectionSecond","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][2]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][2]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPTHistogram[icent][2]); } 

	 fPTHistogram[icent][3] = new TH1F("fPTHistAfterTrackSelectionSecond","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][3]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][3]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPTHistogram[icent][3]); } 

	 
	 // b) Book histogram to hold phi spectra
	 fPhiHistogram[icent][0] = new TH1F("fPhiHistBeforeTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][0]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][0]);

	 fPhiHistogram[icent][1] = new TH1F("fPhiHistAfterTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][1]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][1]);

	 fPhiHistogram[icent][4] = new TH1F("fPhiHistAfterTrackSelectionWeighted","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][4]->Sumw2(); 
	 fPhiHistogram[icent][4]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][4]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][4]);

	 fPhiHistogram[icent][2] = new TH1F("fPhiHistBeforeTrackSelectionSecond","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][2]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][2]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPhiHistogram[icent][2]); } 

	 fPhiHistogram[icent][3] = new TH1F("fPhiHistAfterTrackSelectionSecond","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][3]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][3]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPhiHistogram[icent][3]); } 

	 // c) Book histogram to hold eta distribution before track selection:
	 fEtaHistogram[icent][0] = new TH1F("fEtaHistBeforeTrackSelection","Eta Distribution",1000,-1.,1.); 
	 fEtaHistogram[icent][0]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][0]);

	 fEtaHistogram[icent][1] = new TH1F("fEtaHistAfterTrackSelection","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][1]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][1]);

	 fEtaHistogram[icent][4] = new TH1F("fEtaHistAfterTrackSelectionWeighted","Eta Distribution",1000,-1.,1.); 
	 fEtaHistogram[icent][4]->Sumw2(); 
	 fEtaHistogram[icent][4]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][4]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][4]);

	 fEtaHistogram[icent][2] = new TH1F("fEtaHistBeforeTrackSelectionSecond","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][2]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][2]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fEtaHistogram[icent][2]); } 

	 fEtaHistogram[icent][3] = new TH1F("fEtaHistAfterTrackSelectionSecond","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][3]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][3]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fEtaHistogram[icent][3]); } 

	 // d) Book histogam to hold multiplicty distributions 
	 fMultHistogram[icent][0] = new TH1F("fMultiHistoBeforeMultCut","Multiplicity",30000,0.,30000.); 
	 fMultHistogram[icent][0]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][0]);
	 
	 fMultHistogram[icent][1] = new TH1F("fMultiHistoBeforeTrackSelection","Multiplicity",30000,0.,30000.); 
	 fMultHistogram[icent][1]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][1]);

	 fMultHistogram[icent][2] = new TH1F("fMultiHistoAfterTrackSelection","Multiplicity",30000,0.,30000.);
	 fMultHistogram[icent][2]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][2]);

	 fMultHistogram[icent][3] = new TH1F("fMultiHistoAfterTrackSelectionSecond","Multiplicity",30000,0.,30000.);
	 fMultHistogram[icent][3]->GetXaxis()->SetTitle("Multiplicity M");
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fMultHistogram[icent][3]); }

	 // e) Book histogam for Vertex X 
	 fVertexXHistogram[icent][0] = new TH1F("fVertexXBefore","VertexXBefore",1000,-20.,20.); 
	 fVertexXHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][0]);

	 fVertexXHistogram[icent][1] = new TH1F("fVertexXAfter","VertexXAfter",1000,-20.,20.); 
	 fVertexXHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][1]);

	 // f) Book histogam for Vertex Y 
	 fVertexYHistogram[icent][0] = new TH1F("fVertexYBefore","VertexYBefore",1000,-20.,20.); 
	 fVertexYHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][0]);

	 fVertexYHistogram[icent][1] = new TH1F("fVertexYAfter","VertexYAfter",1000,-20.,20.); 
	 fVertexYHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][1]);

	 // g) Book histogam for Vertex Z 
	 fVertexZHistogram[icent][0] = new TH1F("fVertexZBefore","VertexZBefore",1000,-20.,20.); 
	 fVertexZHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][0]);

	 fVertexZHistogram[icent][1] = new TH1F("fVertexZAfter","VertexZAfter",1000,-20.,20.); 
	 fVertexZHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][1]);



	 // i) Book histogram for number of TPC clustes 
	 fTPCClustersHistogram[icent][0] = new TH1F("fTPCClustersBeforeCut","TPCClustersBeforeCut",170,0.,170.); 
	 fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][0]);

	 fTPCClustersHistogram[icent][1] = new TH1F("fTPCClustersAfterCut","TPCClustersAfterCut",170,0.,170.); 
	 fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][1]);

	 //j) Book histogram for number of ITC clusters 
	 fITSClustersHistogram[icent][0] = new TH1F("fITSClustersBeforeCut","ITSClustersBeforeCut",10,0.,10.); 
	 fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][0]);

	 fITSClustersHistogram[icent][1] = new TH1F("fITSClustersAfterCut","ITSClustersAfterCut",10,0.,10.); 
	 fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][1]);

	 // k) Book histogram for chi square TPC 
	 fChiSquareTPCHistogram[icent][0] = new TH1F("fChiSquareTPCBeforeCut","ChiSquareTPCBeforeCut",1000,0.,20.); 
	 fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][0]);

	 fChiSquareTPCHistogram[icent][1] = new TH1F("fChiSquareTPCAfterCut","ChiSquareTPCAfterCut",1000,0.,20.); 
	 fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][1]);

	  // l) Book histogram for DCAz
	 fDCAzHistogram[icent][0] = new TH1F("fDCAzBeforeCut","DCAzBeforeCut",1000,-10.,10.);  
	 fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][0]);

	 fDCAzHistogram[icent][1] = new TH1F("fDCAzAfterCut","DCAzAfterCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][1]);
	 
	 // m) Book histogram for DCAxy
	 fDCAxyHistogram[icent][0] = new TH1F("fDCAxyBeforeCut","DCAxyBeforeCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][0]);

	 fDCAxyHistogram[icent][1] = new TH1F("fDCAxyAfterCut","DCAxyAfterCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][1]); 

	 // n) Book histogram for Charge
	 fChargeHistogram[icent][0] = new TH1I("ChargeBeforeCut","ChargeBeforeCut",11,-5.5,5.5); 
	 fControlHistogramsList[icent]->Add(fChargeHistogram[icent][0]);

	 fChargeHistogram[icent][1] = new TH1I("ChargeAfterCut","DCAxyAfterCut",11,-5.5,5.5); 
	 fControlHistogramsList[icent]->Add(fChargeHistogram[icent][1]); 

	 // n) Book histogram Centrality 
	 fCentralityHistogram[icent]= new TH1F("fCentralityHistogramAfter","CentralityHistogramAfter",22,0.,110.);
	 fCentralityHistogram[icent]->GetXaxis()->SetTitle("Centrality");
	 fCentralityHistogram[icent]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fCentralityHistogram[icent]);

  }//for(Int_t icent=0; icent<fCentralityBins; icent++)

} //void AliAnalysisTaskStudentsML::BookControlHistograms()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

  for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
 {
	//Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	 fResults[icent] = new TProfile("fResults","Result Analysis First Set Correlators",16,0.,16.,"s"); //centrality dependent output
	 fResults[icent]->GetXaxis()->SetTitle("");
	 fResults[icent]->GetYaxis()->SetTitle("");
	 fResults[icent]->Sumw2();
	 fFinalResultsList[icent]->Add(fResults[icent]);

	 fCovResults[icent] = new TProfile("fCovResults","Result for Covariance Terms",32,0.,32.,"s"); //centrality dependent output
	 fCovResults[icent]->GetXaxis()->SetTitle("");
	 fCovResults[icent]->GetYaxis()->SetTitle("");
	 fCovResults[icent]->Sumw2();
	 fFinalResultsList[icent]->Add(fCovResults[icent]); 

	 fJoinedCovResults[icent] = new TProfile("fJoinedCovResults","Result joined Covariance term calculated as one correlator <z> not product of two correlators <x*y>",16,0.,16.); //centrality dependent output 
	 fJoinedCovResults[icent]->GetXaxis()->SetTitle("");
	 fJoinedCovResults[icent]->GetYaxis()->SetTitle("");
	 fJoinedCovResults[icent]->Sumw2();
	 fFinalResultsList[icent]->Add(fJoinedCovResults[icent]);

	 fMixedParticleHarmonics[icent] = new TProfile("fMixedParticleHarmonics","fMixedParticleHarmonics",2,0.,2.,"s"); //centrality dependent output
	 fMixedParticleHarmonics[icent]->GetXaxis()->SetTitle("");
	 fMixedParticleHarmonics[icent]->GetYaxis()->SetTitle("");
	 fMixedParticleHarmonics[icent]->Sumw2(); 
	 fFinalResultsList[icent]->Add(fMixedParticleHarmonics[icent]);
 } 
 // Book histogram to debug
 fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks",3,0.,3.);
 fHistList->Add(fCounterHistogram);


  //Centrality distribution before cuts
  fCentralityHistogramBefore = new TH1F("fCentralityHistogramBefore","fCentralityHistogramBefore",22,0.,110.);
  fCentralityHistogramBefore->GetXaxis()->SetTitle("Centrality");
  fCentralityHistogramBefore->SetLineColor(4);
  fHistList->Add(fCentralityHistogramBefore);

 //Profiles to save the current cut values 
 //Profile to save the cut values for event selection
  fProfileEventCuts = new TProfile("", "", 16, 0., 16.);
  fProfileEventCuts->SetName("fProfileEventCuts");
  fProfileEventCuts->SetTitle("Configuration of the event selection");
  fProfileEventCuts->SetStats(kFALSE);
  fProfileEventCuts->GetXaxis()->SetBinLabel(1, "1st centrality");
  fProfileEventCuts->GetXaxis()->SetBinLabel(2, "2nd centrality");
  fProfileEventCuts->GetXaxis()->SetBinLabel(3, "PV_{x} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(4, "PV_{x} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(5, "PV_{y} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(6, "PV_{y} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(7, "PV_{z} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(8, "PV_{z} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(9, "Multiplicity min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(10, "1st filter");
  fProfileEventCuts->GetXaxis()->SetBinLabel(11, "2nd filter");
  fProfileEventCuts->GetXaxis()->SetBinLabel(12, "HMO slope lower line");
  fProfileEventCuts->GetXaxis()->SetBinLabel(13, "HMO slope upper line");
  fProfileEventCuts->GetXaxis()->SetBinLabel(14, "HMO axis lower line");
  fProfileEventCuts->GetXaxis()->SetBinLabel(15, "HMO axis upper line");
  fProfileEventCuts->GetXaxis()->SetBinLabel(16, "Use Reco-Kine-Table");
  fHistList->Add(fProfileEventCuts);

 //Profile to save the cut values for track selection
  fProfileTrackCuts = new TProfile("", "", 18, 0., 18.);
  fProfileTrackCuts->SetName("fProfileTrackCuts");
  fProfileTrackCuts->SetTitle("Configuration of the track selection");
  fProfileTrackCuts->SetStats(kFALSE);
  fProfileTrackCuts->GetXaxis()->SetBinLabel(1, "p_{T} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(2, "p_{T} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(3, "#eta min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(4, "#eta max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(5, "N_{TPC} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(6, "#chi^{2} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(7, "#chi^{2} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(8, "#chi^{2} Method");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(9, "N_{ITS} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(10, "DCA_{xy} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(11, "DCA_{z} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(12, "Charge Cut?");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(13, "Charge"); 
  fProfileTrackCuts->GetXaxis()->SetBinLabel(14, "Fisher Yates?"); 
  fProfileTrackCuts->GetXaxis()->SetBinLabel(15, "Keeping Percentege");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(16, "p_{T} weights?");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(17, "#eta weights?");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(18, "#phi weights?"); 
  fHistList->Add(fProfileTrackCuts);

 Cosmetics();
 
} // void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::Cosmetics()
{
 // Book everything here.
  

} // void Cosmetics()

//==========================================================================================================================================================================

Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliAODEvent *aAODevent){

  //a) Protection against NULL-Pointers
  //b) Check Centrality
  //c) Cuts on AliAODVertex:
  //d) remove high multiplicity outliers

 
  TString sMethod = "Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliAODEvent *aAODevent)";
 
  //a) Protection against NULL-Pointers
  if(!aAODevent){return kFALSE;}
  fCounterHistogram->Fill(1.5); // counter hist 2nd bin

  //b) Check Centrality
  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if(!ams){return kFALSE;}
  fCounterHistogram->Fill(2.5); // counter hist 3rd bin


  //Get Centrality Bin
  Int_t CentBin = SelectCentrality(aAODevent);
  Float_t CentralityValue = 0.;

  if(fCentralityEstimator == "V0M" ) { CentralityValue = ams->GetMultiplicityPercentile("V0M"); }
  else if(fCentralityEstimator == "CL1") { CentralityValue = ams->GetMultiplicityPercentile("CL1"); }
  else {Fatal(sMethod.Data(), "FATAL: no valid centrality estimator!");} 


  fCentralityHistogramBefore->Fill(CentralityValue);

  if (CentBin < 0){ return kFALSE; }
  
  if(bSaveAllQA){fCentralityHistogram[CentBin]->Fill(CentralityValue);}
  
  // c) Cuts on AliAODVertex:
  AliAODVertex *avtx = (AliAODVertex*)aAODevent->GetPrimaryVertex();
 
  if(bSaveAllQA){
	  fVertexXHistogram[CentBin][0]->Fill(avtx->GetX());
	  fVertexYHistogram[CentBin][0]->Fill(avtx->GetY());
	  fVertexZHistogram[CentBin][0]->Fill(avtx->GetZ());
  }//if(bSaveAllQA)

  if(bCutOnVertexX)
  {
   if(avtx->GetX() < fMinVertexX) return kFALSE;
   if(avtx->GetX() > fMaxVertexX) return kFALSE;
  }
  if(bCutOnVertexY)
  {
   if(avtx->GetY() < fMinVertexY) return kFALSE;
   if(avtx->GetY() > fMaxVertexY) return kFALSE;
  }
  if(bCutOnVertexZ)
  {
   if(avtx->GetZ() < fMinVertexZ) return kFALSE;
   if(avtx->GetZ() > fMaxVertexZ) return kFALSE;
  }

  //d) remove high multiplicity outliers

  if(bMultCut)
  {
  	Int_t nTracks = aAODevent->GetNumberOfTracks(); // number of all tracks in current event 
  	Int_t nCounterMainFilter=0; //Counter for MainFilter
  	Int_t nCounterSecondFilter=0; //Counter for SecondFilter

	if(bSaveAllQA){fMultHistogram[CentBin][0]->Fill(nTracks);} //multiplicity distribution before high multiplicity outlier removal
  	for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 	{
  	  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack)); 
   	  if(!aTrack){continue;} 
   	  if(aTrack->TestFilterBit(fMainFilter)){nCounterMainFilter++; } //one more track with main filter
   	  if(aTrack->TestFilterBit(fSecondFilter)){ nCounterSecondFilter++; } //one more track with second filter
  	}//for(Int_t iTrack=0;iTrack<nTracks;iTrack++)
	
	if( (Float_t)nCounterMainFilter > (fSlopeUpperLine*(Float_t)nCounterSecondFilter + fAxisUpperLine) ) return kFALSE;
	if( (Float_t)nCounterMainFilter < (fSlopeLowerLine*(Float_t)nCounterSecondFilter + fAxisLowerLine) ) return kFALSE;

	nCounterMainFilter=0;
  	nCounterSecondFilter=0;

  }//end if(bMultCut)

  return kTRUE;
 }//end  Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliAODEvent *aAODevent)

//==========================================================================================================================================================================

 Bool_t AliAnalysisTaskStudentsML::TrackSelection(AliAODVertex *aPrimaryVertex, AliAODTrack *aTrack) 
 {

        // example variables for each track:
 	/*Double_t px = aTrack->Px(); // x-component of momenta
 	Double_t py = aTrack->Py(); // y-component of momenta
 	Double_t pz = aTrack->Pz(); // z-component of momenta
 	Double_t e = aTrack->E();  // energy
 	Double_t phi = aTrack->Phi(); // azimuthal angle*/

 	Double_t eta = aTrack->Eta(); // pseudorapidity
 	Double_t pt = aTrack->Pt(); // Pt (transverse momentum)
	Int_t NumberOfTPCClusters = aTrack->GetTPCNcls(); //number of TPC clusters of the track
	Int_t NumberOfITSClusters = aTrack->GetITSNcls(); //number of ITS clusters of the track
	Double_t charge = aTrack->Charge(); // charge

	// Get the DCA information (cf PWGCF/EBYE/BalanceFunctions/AliAnalysisTaskBFPsi.cxx)
	Float_t ValueDCAxy = 999.;   // DCA in the xy-plane.
	Float_t ValueDCAz = 999.;    // DCA along z.

	if (fMainFilter == 128)  // These methods work only for constrained TPConly tracks.
	{ // These two quantities are the DCA from global tracks but not what we will cut on.
	  ValueDCAxy = aTrack->DCA();
	  ValueDCAz = aTrack->ZAtDCA();
	}
	else  // For the unconstrained tracks.
	{
	  Double_t v[3];    // Coordinates of the PV?
	  Double_t pos[3];  // Coordinates of the track closest to PV?

          aPrimaryVertex->GetXYZ(v);
	  aTrack->GetXYZ(pos);
	  ValueDCAxy = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
	  ValueDCAz = pos[2] - v[2];
	}

	// Get the chi^2 per TPC cluster.
	Float_t ChiSquareInTPC = 999.;
	/// Personal method, should be equal to GetTPCchi2perCluster()
	if (fChooseChiSquareMethod == 1) {ChiSquareInTPC = (aTrack->GetTPCchi2())/(aTrack->GetNcls(1));}
	else if (fChooseChiSquareMethod == 2) {ChiSquareInTPC = aTrack->GetTPCchi2perCluster();}
	else if (fChooseChiSquareMethod == 3) {ChiSquareInTPC = aTrack->GetTPCchi2perNDF();}
	else {ChiSquareInTPC = aTrack->Chi2perNDF();}


	if(bCutOnEta) 
	{
	  if(eta<fMinEtaCut) return kFALSE;
	  if(eta>fMaxEtaCut) return kFALSE;
	}

	if(bCutOnPt) 
	{
	  if(pt<fMinPtCut) return kFALSE;
	  if(pt>fMaxPtCut) return kFALSE;
	}

	if(bNumberTPCCluster) 
	{
	  if(NumberOfTPCClusters<fMinTPCCluster) return kFALSE;
	  
	}

	if(bNumberITSCluster) 
	{
	  if(NumberOfITSClusters<fMinITSCluster) return kFALSE;
	}

	if(bChiSquareTPC) 
	{
	  if(ChiSquareInTPC<fMinChiSquareTPC) return kFALSE;
	  if(ChiSquareInTPC>fMaxChiSquareTPC) return kFALSE;
	}

	if(bDCAz) 
	{
	  if(TMath::Abs(ValueDCAz)>fMaxDCAz) return kFALSE;
	}

	if(bDCAxy) 
	{
	  if(TMath::Abs(ValueDCAxy)>fMaxDCAxy) return kFALSE;
	}

	if(bChargeCut)
	{
	  if(bChargePos) { if(charge<=0.) return kFALSE; } 
	  else { if(charge>=0.) return kFALSE; } 
	}

    return kTRUE;

 }// end AliAnalysisTaskStudentsML::PhysicsSelection()

//==========================================================================================================================================================================
// Overloading functions for MC Events

 Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(Int_t CentBin, AliMCEvent *aMCKineEvent)
 {

  //a) Reset Histo's
  //b) Protection against NULL-Pointers
  //c) Cuts on AliAODVertex:
  //d) remove high multiplicity outliers

  
  //a) Reset Histo's (are already filled in QA for AOD)
  fCounterHistogram->SetBinContent(2,0);
  fCounterHistogram->SetBinContent(3,0);
  if(bSaveAllQA)
  {
	  fVertexXHistogram[CentBin][0]->Reset(); 
	  fVertexYHistogram[CentBin][0]->Reset();  
	  fVertexZHistogram[CentBin][0]->Reset(); 
	  fMultHistogram[CentBin][0]->Reset(); 
  }//if(bSaveAllQA)

  //b) Protection against NULL-Pointers
  if(!aMCKineEvent){return kFALSE;}
  fCounterHistogram->Fill(1.5); // counter hist 2nd bin

  // c) Cuts on AliAODVertex:
  AliMCVertex *avtx = (AliMCVertex*)aMCKineEvent->GetPrimaryVertex();
 
  if(bSaveAllQA)
  {
	  fVertexXHistogram[CentBin][0]->Fill(avtx->GetX()); 
	  fVertexYHistogram[CentBin][0]->Fill(avtx->GetY()); 
	  fVertexZHistogram[CentBin][0]->Fill(avtx->GetZ()); 
 }//if(bSaveAllQA)

  if(bCutOnVertexX)
  {
   if(avtx->GetX() < fMinVertexX) return kFALSE;
   if(avtx->GetX() > fMaxVertexX) return kFALSE;
  }
  if(bCutOnVertexY)
  {
   if(avtx->GetY() < fMinVertexY) return kFALSE;
   if(avtx->GetY() > fMaxVertexY) return kFALSE;
  }
  if(bCutOnVertexZ)
  {
   if(avtx->GetZ() < fMinVertexZ) return kFALSE;
   if(avtx->GetZ() > fMaxVertexZ) return kFALSE;
  }

  //d) removal of high multiplicity outliers not needed in MC Kine level

  return kTRUE;
 }//end  AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliMCEvent *aMCKineEvent)

//==========================================================================================================================================================================
 Bool_t AliAnalysisTaskStudentsML::TrackSelection(AliAODMCParticle *aMCtrack)
 {
   	// example variables for each track:
	Double_t eta = aMCtrack->Eta(); // pseudorapidity
 	Double_t pt = aMCtrack->Pt(); // Pt (transverse momentum)

	if(bCutOnEta) 
	{
	  if(eta<fMinEtaCut) return kFALSE;
	  if(eta>fMaxEtaCut) return kFALSE;
	}

	if(bCutOnPt) 
	{
	  if(pt<fMinPtCut) return kFALSE;
	  if(pt>fMaxPtCut) return kFALSE;
	}

    return kTRUE;

 }// end AliAnalysisTaskStudentsML::PhysicsSelection(AliAODMCParticle *aMCtrack)

//==========================================================================================================================================================================

Int_t AliAnalysisTaskStudentsML::SelectCentrality(AliAODEvent *aAODevent)
{

  TString sMethod = "Int_t AliAnalysisTaskStudentsML::SelectCentrality(AliAODEvent *aAODevent)";

  //if this functions returns a negative value -> Error. No Centrality could be selected

  if(!aAODevent){return -1;}

  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if(!ams){return -1;}

  //Access the Centrality of the event
  Float_t CentralityValue = 0.;

  if(fCentralityEstimator == "V0M" ) { CentralityValue = ams->GetMultiplicityPercentile("V0M"); }
  else if(fCentralityEstimator == "CL1") { CentralityValue = ams->GetMultiplicityPercentile("CL1"); }
  else {Fatal(sMethod.Data(), "FATAL: no valid centrality estimator!");} 

  //Check for centrality bin
  for(Int_t icent=0; icent<fCentralityBins+1; icent++) //loop over all centrality bins
  {
	if(fcentralityArray[icent]<0){return -1;}
	if(CentralityValue >= fcentralityArray[icent]) { continue; } 
	else { return icent-1; } 
  }	

  //We went through all centrality edges without returning. This means, that the measured value is bigger than the maximum centrality that we want for our analyis
  return -1;

} //Int_t AliAnalysisTaskStudentsML::SelectCentrality(AliAODEvent *aAODevent)

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::FisherYatesRandomizing(Int_t Mult, Int_t *RandomIndex)
{
  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);
        
  if(Mult <= 0) return;
            

  for(Int_t i=0;i<Mult;i++)
  {
        RandomIndex[i] = i;
  }
    
  for(Int_t i=Mult-1;i>=1;i--)
  {
         Int_t j = gRandom->Integer(i+1);
            
         Int_t Temp_Storage = RandomIndex[j];
         RandomIndex[j] = RandomIndex[i];
         RandomIndex[i] = Temp_Storage;
        
  } // end of for(Int_t i=nPrim-1;i>=1;i--)

  return;
}

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::CalculateWeight(Int_t CentBin, Int_t RunNumber, Double_t* weights, Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta)
{
  //a) Get Array-Index of the Run
  Int_t RunIndex = GetRunIndex(RunNumber);

  //b) Loop over all particles in the event and asign weights
  for(Int_t i=0; i<Multi; i++)
  {
    //First, by default, we have unit weights
    Double_t weight_phi = 1.;
    Double_t weight_pt = 1.;
    Double_t weight_eta = 1.;
    Int_t iBin = 0;


    if(bUsePhiWeights)
    { 
      iBin = fHistoPhiWeight[CentBin][RunIndex]->FindBin(angles[i]); 
      weight_phi = fHistoPhiWeight[CentBin][RunIndex]->GetBinContent(iBin); 
      if(bSaveAllQA){ fPhiHistogram[CentBin][4]->Fill(angles[i], weight_phi); }
    }
    if(bUsePtWeights)
    {
      iBin = fHistoPtWeight[CentBin][RunIndex]->FindBin(pt[i]); 
      weight_pt = fHistoPtWeight[CentBin][RunIndex]->GetBinContent(iBin);
      if(bSaveAllQA){ fPTHistogram[CentBin][4]->Fill(pt[i], weight_pt); }
    }
    if(bUseEtaWeights)
    { 
      iBin = fHistoEtaWeight[CentBin][RunIndex]->FindBin(eta[i]); 
      weight_eta = fHistoEtaWeight[CentBin][RunIndex]->GetBinContent(iBin);
      if(bSaveAllQA){ fEtaHistogram[CentBin][4]->Fill(eta[i], weight_eta); } 
    }
   
    //Final overall weight
    weights[i] = weight_phi*weight_pt*weight_eta;

  }//for(Int_t i=0; i<Multi; i++)

} //void AliAnalysisTaskStudentsML::CalculateWeight(Int_t CentBin, Int_t RunNumber, Double_t* weights, Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta)

//==========================================================================================================================================================================

Int_t AliAnalysisTaskStudentsML::GetRunIndex(Int_t runNumber)
{
/* Return for the given run the index in the run-by-run arrays.                              */
  TString sMethod = "Int_t AliAnalysisTaskStudentsML::GetRunIndex()";
  Int_t cRun = -1; // Current index in the loop.

// Find the position of the given run into the list of runs.
  for (Int_t iRun = 0; iRun < fNumberRuns; iRun++)
  {
    if (fListRuns[iRun] == runNumber)
    {
      cRun = iRun;
      break;
    } // End: for (Int_t iRun = 0; iRun < fNumberRuns; iRun++).
  } // End: iRun.

 if(cRun == -1){Fatal(sMethod.Data(), "FATAL: Run Number not in List of Runs!");}  

  return cRun;
} // End: Int_t GetRunIndex(Int_t).

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::SetListOfRuns(TString dataPeriod) //Called in Constructor
{
/* Set the list of runs to use according to the chosen data-taking period................... */
  TString sMethod = "void AliAnalysisTaskStudentsML::SetListOfRuns()";

  if (dataPeriod == "LHC10h")
  {
    fNumberRuns = 90;
    Int_t listRuns[90] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
    for (Int_t i = 0; i < fNumberRuns; i++) {fListRuns[i] = listRuns[i];}
  } // End: if (dataPeriod == "LHC10h").
  else {Fatal(sMethod.Data(), "FATAL: not a valid data period!");}

} // End: void SetListOfRuns(TString).

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::SetInputParticleWeights(TString fileWeight)
{
// Setter to open the external file with the particle weights and import them in the task....
// a.)  Open the external file.                                                                  
// b.)  Parse the runs.                                                                          
// b.1) Open the TDirectoryFile for the current run.                                            
// b.2) Open the list for the current centrality range.                                                                                         
// c.)  Close the external file.                                                               

  TString sMethod = "void AliAnalysisTaskStudentsML::SetInputParticleWeights()"; 

  // a.) Open the external file.
  TFile *weightsFile = TFile::Open(Form("%s", fileWeight.Data()), "READ");
  if (!weightsFile) {Fatal(sMethod.Data(), "ERROR 404: File not found");}


 for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
 {

	//Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	  // b.) Parse the runs.
	  for (Int_t iRun = 0; iRun < fNumberRuns; iRun++)
	  {
	    // b.1) Open the TDirectoryFile for the current run.
	    Int_t runNumber = fListRuns[iRun];
	    //printf("Run number: %d\n", runNumber);
	    TDirectoryFile *runTDF = dynamic_cast<TDirectoryFile*>(weightsFile->Get(Form("%d", runNumber)));
	    if (!runTDF) {Fatal(sMethod.Data(), "ERROR: Directory not found");}

	    // b.2) Open the list for the current centrality range.
	    TList *centralityList = dynamic_cast<TList*>(runTDF->Get(Form("Centrality-%.1f-%.1f", fcentralityArray[icent], fcentralityArray[icent+1])));
	    if (!centralityList) {Fatal(sMethod.Data(), "ERROR: List not found");}

	      fHistoPtWeight[icent][iRun] = dynamic_cast<TH1F*>(centralityList->FindObject("pt-weight"));
	      if (!fHistoPtWeight[icent][iRun]) {Fatal(sMethod.Data(), "ERROR: pt-weight histogram not found");}
	      else {fHistoPtWeight[icent][iRun]->SetDirectory(0);} // Kill the default ownership.

	      fHistoEtaWeight[icent][iRun] = dynamic_cast<TH1F*>(centralityList->FindObject("eta-weight"));
	      if (!fHistoEtaWeight[icent][iRun]) { Fatal(sMethod.Data(), "ERROR: eta-weight histogram not found"); }
	      else {fHistoEtaWeight[icent][iRun]->SetDirectory(0);}  // Kill the default ownership.

	      fHistoPhiWeight[icent][iRun] = dynamic_cast<TH1F*>(centralityList->FindObject("phi-weight"));
	      if (!fHistoPhiWeight[icent][iRun]) {Fatal(sMethod.Data(), "ERROR: phi-weight histogram not found");}
	      else {fHistoPhiWeight[icent][iRun]->SetDirectory(0);}  // Kill the default ownership.
	    

	     delete centralityList;
	     delete runTDF;

	  } // End: iRun.

 }//for(Int_t icent=0; icent<fCentralityBins; icent++)

  // c.) Close the external file.
  weightsFile->Close();
  delete weightsFile;

} // void AliAnalysisTaskStudentsML::SetInputParticleWeights(TString fileWeight)

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights)
{
 // Calculate Q-vectors.
 // a) Make sure all Q-vectors are initially zero;
 // b) Calculate Q-vectors for available angles and weights. 

 // a) Make sure all Q-vectors are initially zero:
 for(Int_t h=0;h<113;h++)
 {
  for(Int_t p=0;p<15;p++)
  {
   fQvector[h][p] = TComplex(0.,0.);
  } //  for(Int_t p=0;p<kMaxPower;p++)
 } // for(Int_t h=0;h<kMaxHarmonic;h++)

 // b) Calculate Q-vectors for available angles and weights: 
 Double_t dPhi2 = 0.; // particle angle
 Double_t wPhi = 1.; // particle weight
 Double_t wPhiToPowerP = 1.; // particle weight raised to power p
 for(Int_t i=0;i<CalculateQvectors_nParticles;i++) // loop over particles
 {
  dPhi2 = CalculateQvectors_angles[i];
  if(bUseWeights){wPhi = CalculateQvectors_weights[i];} //Change some point
  for(Int_t h=0;h<113;h++)
  {
   for(Int_t p=0;p<15;p++)
   {
    if(bUseWeights){wPhiToPowerP = pow(wPhi,p);}
    fQvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi2),wPhiToPowerP*TMath::Sin(h*dPhi2));
   } //  for(Int_t p=0;p<kMaxPower;p++)
  } // for(Int_t h=0;h<kMaxHarmonic;h++)
 } //  for(Int_t i=0;i<fParticles;i++) // loop over particles


} // void AliAnalysisTaskStudentsML::CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights)

//==========================================================================================================================================================================

TComplex AliAnalysisTaskStudentsML::CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B){

    //Reset
    TComplex Q_A = TComplex(0.,0.);
    TComplex Q_2A = TComplex(0.,0.);
    TComplex Q_B = TComplex(0.,0.);
    TComplex Q_2B = TComplex(0.,0.);
    TComplex v_AB = TComplex(0.,0.);
    
    for(Int_t i=0; i<M_A; i++){
        Q_A += TComplex(TMath::Cos(Harm*Ang_A[i]),TMath::Sin(Harm*Ang_A[i]));
        Q_2A += TComplex(TMath::Cos(2.*Harm*Ang_A[i]),TMath::Sin(2.*Harm*Ang_A[i]));
    }
    
    for(Int_t i=0; i<M_B; i++){
        Q_B += TComplex(TMath::Cos(Harm*Ang_B[i]),TMath::Sin(Harm*Ang_B[i]));
        Q_2B += TComplex(TMath::Cos(2.*Harm*Ang_B[i]),TMath::Sin(2.*Harm*Ang_B[i]));
    }
    
  v_AB = Q_A*Q_A*TComplex::Conjugate(Q_B)*TComplex::Conjugate(Q_B) - Q_2A*TComplex::Conjugate(Q_B)*TComplex::Conjugate(Q_B) - Q_A*Q_A*TComplex::Conjugate(Q_2B) + Q_2A*TComplex::Conjugate(Q_2B);
    

    return v_AB;

 } //TComplex AliAnalysisTaskStudentsML::CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B)

//==========================================================================================================================================================================

TComplex AliAnalysisTaskStudentsML::Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][p];} 
 return TComplex::Conjugate(fQvector[-n][p]);
 
} // TComplex AliAnalysisTaskStudentsML::Q(Int_t n, Int_t p)


//==========================================================================================================================================================================


TComplex AliAnalysisTaskStudentsML::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Int_t h13, Int_t h14, Double_t *Correlation_Data)
{
	
       if(h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12+h13+h14!=0.){return;} //protection against anisotropic correlators
	
       // Calculate n-particle correlations from Q-vectors (using recursion):	
         
        if(2==Number)
        {
         Int_t harmonicsTwoNum[2] = {h1,h2};     
         Int_t harmonicsTwoDen[2] = {0,0};      
         Double_t wTwoRecursion = Recursion(2,harmonicsTwoDen).Re(); 
         TComplex twoRecursion = Recursion(2,harmonicsTwoNum)/wTwoRecursion;

	 Correlation_Data[0] = twoRecursion.Re(); //<cos(h1*phi1+h2*phi2)>
	 Correlation_Data[1] = wTwoRecursion;	  // weight
	 Correlation_Data[2] = twoRecursion.Im(); // <sin(h1*phi1+h2*phi2)>
	
         }//  2-p correlation
        
        if(3==Number)
        {
         Int_t harmonicsThreeNum[3] = {h1,h2,h3};       
         Int_t harmonicsThreeDen[3] = {0,0,0}; 
         Double_t wThreeRecursion = Recursion(3,harmonicsThreeDen).Re();      
         TComplex threeRecursion = Recursion(3,harmonicsThreeNum)/wThreeRecursion;

	 Correlation_Data[0] = threeRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3)>
	 Correlation_Data[1] = wThreeRecursion;	    // weight
	 Correlation_Data[2] = threeRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3)>

         } //  3-p correlation
        
        if(4==Number)
        {
         Int_t harmonicsFourNum[4] = {h1,h2,h3,h4};       
         Int_t harmonicsFourDen[4] = {0,0,0,0}; 
         Double_t wFourRecursion = Recursion(4,harmonicsFourDen).Re();      
         TComplex fourRecursion = Recursion(4,harmonicsFourNum)/wFourRecursion;

	 Correlation_Data[0] = fourRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
	 Correlation_Data[1] = wFourRecursion;     // weight
	 Correlation_Data[2] = fourRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>

         }//  4-p correlation
        
        if(5==Number)
        {
         Int_t harmonicsFiveNum[5] = {h1,h2,h3,h4,h5};       
         Int_t harmonicsFiveDen[5] = {0,0,0,0,0};       
         Double_t wFiveRecursion = Recursion(5,harmonicsFiveDen).Re();
         TComplex fiveRecursion = Recursion(5,harmonicsFiveNum)/wFiveRecursion;

	 Correlation_Data[0] = fiveRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
	 Correlation_Data[1] = wFiveRecursion;     // weight
	 Correlation_Data[2] = fiveRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>

        }//  5-p correlation

        if(6==Number)
        {
         Int_t harmonicsSixNum[6] = {h1,h2,h3,h4,h5,h6};       
         Int_t harmonicsSixDen[6] = {0,0,0,0,0,0}; 
         Double_t wSixRecursion = Recursion(6,harmonicsSixDen).Re();      
         TComplex sixRecursion = Recursion(6,harmonicsSixNum)/wSixRecursion;

	 Correlation_Data[0] = sixRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
	 Correlation_Data[1] = wSixRecursion;     // weight
	 Correlation_Data[2] = sixRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>

         }//  6-p correlation
        
        
        if(7==Number)
        {
         Int_t harmonicsSevenNum[7] = {h1,h2,h3,h4,h5,h6,h7};       
         Int_t harmonicsSevenDen[7] = {0,0,0,0,0,0,0};  
         Double_t wSevenRecursion = Recursion(7,harmonicsSevenDen).Re();     
         TComplex sevenRecursion = Recursion(7,harmonicsSevenNum)/wSevenRecursion;

	 Correlation_Data[0] = sevenRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
	 Correlation_Data[1] = wSevenRecursion;     // weight
	 Correlation_Data[2] = sevenRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>

        }//  7-p correlation
        
        
        if(8==Number)
        {
         Int_t harmonicsEightNum[8] = {h1,h2,h3,h4,h5,h6,h7,h8};       
         Int_t harmonicsEightDen[8] = {0,0,0,0,0,0,0,0};       
         Double_t wEightRecursion = Recursion(8,harmonicsEightDen).Re();
         TComplex eightRecursion = Recursion(8,harmonicsEightNum)/wEightRecursion;

	 Correlation_Data[0] = eightRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
	 Correlation_Data[1] = wEightRecursion;     // weight
	 Correlation_Data[2] = eightRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
        
        }//  8-p correlation


        if(9==Number)
        {
         Int_t harmonicsNineNum[9] = {h1,h2,h3,h4,h5,h6,h7,h8,h9};       
         Int_t harmonicsNineDen[9] = {0,0,0,0,0,0,0,0,0};    
         Double_t wNineRecursion = Recursion(9,harmonicsNineDen).Re();   
         TComplex nineRecursion = Recursion(9,harmonicsNineNum)/wNineRecursion;

	 Correlation_Data[0] = nineRecursion.Re(); 
	 Correlation_Data[1] = wNineRecursion;
	 Correlation_Data[2] = nineRecursion.Im();        

        }//  9-p correlation

        if(10==Number)
        {
         Int_t harmonicsTenNum[10] = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10};       
         Int_t harmonicsTenDen[10] = {0,0,0,0,0,0,0,0,0,0};  
         Double_t wTenRecursion = Recursion(10,harmonicsTenDen).Re();     
         TComplex tenRecursion = Recursion(10,harmonicsTenNum)/wTenRecursion;

	 Correlation_Data[0] = tenRecursion.Re(); 
	 Correlation_Data[1] = wTenRecursion;
	 Correlation_Data[2] = tenRecursion.Im();
        
        }//  10-p correlation

        if(12==Number)
        {
         Int_t harmonicsTwelveNum[12] = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12};       
         Int_t harmonicsTwelveDen[12] = {0,0,0,0,0,0,0,0,0,0,0,0};     
         Double_t wTwelveRecursion = Recursion(12,harmonicsTwelveDen).Re();  
         TComplex twelveRecursion = Recursion(12,harmonicsTwelveNum)/wTwelveRecursion;

	 Correlation_Data[0] = twelveRecursion.Re(); 
	 Correlation_Data[1] = wTwelveRecursion;
	 Correlation_Data[2] = twelveRecursion.Im();
        
        }//  12-p correlation

        if(14==Number)
        {
         Int_t harmonicsFourteenNum[14] = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14};       
         Int_t harmonicsFourteenDen[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
         Double_t wFourteenRecursion = Recursion(14,harmonicsFourteenDen).Re();      
         TComplex fourteenRecursion = Recursion(14,harmonicsFourteenNum)/wFourteenRecursion;

         Correlation_Data[0] = fourteenRecursion.Re(); 
	 Correlation_Data[1] = wFourteenRecursion;
	 Correlation_Data[2] = fourteenRecursion.Im();

        }//  14-p correlation

        if(Number!=2 && Number!=3 && Number!=4 && Number!=5 && Number!=6 && Number!=7 && Number!=8 && Number!=9 && Number!=10 && Number!=12 && Number!=14) { return; }
      
 }//void Correlation() 

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array)
{

 if(MainTask_Mult>=fMinNumberPart) //do the correlation only if there are more than 8 particles in the event
 { 

   // Calculate Q-vectors for available angles and weights;
   this->CalculateQvectors(MainTask_Mult, MainTask_Angle_Array, MainTask_Weight_Array);


    //Array for temp. data storing
    Double_t* Data_Correlation = new Double_t[3]; //cos, weight, sin
    
     
    Double_t CorrelationNum[8]={0.};
    Double_t Weight_CorrelationNum[8]={0.};
    Double_t CorrelationDenom[8]={0.};
    Double_t Weight_CorrelationDenom[8]={0.};

    Double_t CorrelationJoinedCov[8]={0.};
    Double_t Weight_CorrelationJoinedCov[8]={0.};

    if(fNumber!=0)
    {
	this->Correlation(fNumber, fa1,fa2,fa3,fa4,fa5,fa6,fa7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

   	CorrelationNum[0] = Data_Correlation[0];
	Weight_CorrelationNum[0] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumber;

	this->Correlation(Number_Denom, fa1,-1.*fa1,fa2,-1.*fa2,fa3,-1.*fa3,fa4,-1.*fa4,fa5,-1.*fa5,fa6,-1.*fa6,fa7,-1.*fa7, Data_Correlation);  

	CorrelationDenom[0] = Data_Correlation[0];
	Weight_CorrelationDenom[0] = Data_Correlation[1];


	if(fNumber <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumber;

		this->Correlation(Number_JoinedCov, fa1,fa2,fa3,fa4,fa1,-1.*fa1,fa2,-1.*fa2,fa3,-1.*fa3,fa4,-1.*fa4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[0] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[0] = Data_Correlation[1];
	}


    } //if(fNumber!=0)



    if(fNumberSecond!=0)
    {
	this->Correlation(fNumberSecond, fb1,fb2,fb3,fb4,fb5,fb6,fb7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

   	CorrelationNum[1] = Data_Correlation[0];
	Weight_CorrelationNum[1] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberSecond;

	this->Correlation(Number_Denom, fb1,-1.*fb1,fb2,-1.*fb2,fb3,-1.*fb3,fb4,-1.*fb4,fb5,-1.*fb5,fb6,-1.*fb6,fb7,-1.*fb7, Data_Correlation);  

	CorrelationDenom[1] = Data_Correlation[0];
	Weight_CorrelationDenom[1] = Data_Correlation[1];

	if(fNumberSecond <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberSecond;

		this->Correlation(Number_JoinedCov, fb1,fb2,fb3,fb4,fb1,-1.*fb1,fb2,-1.*fb2,fb3,-1.*fb3,fb4,-1.*fb4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[1] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[1] = Data_Correlation[1];
	}

    } //if(fNumberSecond!=0)



    if(fNumberThird!=0)
    {
	this->Correlation(fNumberThird, fd1,fd2,fd3,fd4,fd5,fd6,fd7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

   	CorrelationNum[2] = Data_Correlation[0];
	Weight_CorrelationNum[2] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberThird;

	this->Correlation(Number_Denom, fd1,-1.*fd1,fd2,-1.*fd2,fd3,-1.*fd3,fd4,-1.*fd4,fd5,-1.*fd5,fd6,-1.*fd6,fd7,-1.*fd7, Data_Correlation);  

	CorrelationDenom[2] = Data_Correlation[0];
	Weight_CorrelationDenom[2] = Data_Correlation[1];

	if(fNumberThird <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberThird;

		this->Correlation(Number_JoinedCov, fd1,fd2,fd3,fd4,fd1,-1.*fd1,fd2,-1.*fd2,fd3,-1.*fd3,fd4,-1.*fd4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[2] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[2] = Data_Correlation[1];
	}


    } //if(fNumberThird!=0)



    if(fNumberFourth!=0)
    {
	this->Correlation(fNumberFourth, fe1,fe2,fe3,fe4,fe5,fe6,fe7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

   	CorrelationNum[3] = Data_Correlation[0];
	Weight_CorrelationNum[3] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberFourth;

	this->Correlation(Number_Denom, fe1,-1.*fe1,fe2,-1.*fe2,fe3,-1.*fe3,fe4,-1.*fe4,fe5,-1.*fe5,fe6,-1.*fe6,fe7,-1.*fe7, Data_Correlation);  

   	CorrelationDenom[3] = Data_Correlation[0];
	Weight_CorrelationDenom[3] = Data_Correlation[1];

	if(fNumberFourth <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberFourth;

		this->Correlation(Number_JoinedCov, fe1,fe2,fe3,fe4,fe1,-1.*fe1,fe2,-1.*fe2,fe3,-1.*fe3,fe4,-1.*fe4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[3] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[3] = Data_Correlation[1];
	}


    }//if(fNumberFourth!=0)


    if(fNumberFifth!=0)
    {
	this->Correlation(fNumberFifth, ff1,ff2,ff3,ff4,ff5,ff6,ff7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

	CorrelationNum[4] = Data_Correlation[0];
	Weight_CorrelationNum[4] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberFifth;

	this->Correlation(Number_Denom, ff1,-1.*ff1,ff2,-1.*ff2,ff3,-1.*ff3,ff4,-1.*ff4,ff5,-1.*ff5,ff6,-1.*ff6,ff7,-1.*ff7, Data_Correlation);  

	CorrelationDenom[4] = Data_Correlation[0];
	Weight_CorrelationDenom[4] = Data_Correlation[1];

	if(fNumberFifth <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberFifth;

		this->Correlation(Number_JoinedCov, ff1,ff2,ff3,ff4,ff1,-1.*ff1,ff2,-1.*ff2,ff3,-1.*ff3,ff4,-1.*ff4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[4] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[4] = Data_Correlation[1];
	}

    } //if(fNumberFifth!=0)



    if(fNumberSixth!=0)
    {
	this->Correlation(fNumberSixth, fg1,fg2,fg3,fg4,fg5,fg6,fg7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

	CorrelationNum[5] = Data_Correlation[0];
	Weight_CorrelationNum[5] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberSixth;

	this->Correlation(Number_Denom, fg1,-1.*fg1,fg2,-1.*fg2,fg3,-1.*fg3,fg4,-1.*fg4,fg5,-1.*fg5,fg6,-1.*fg6,fg7,-1.*fg7, Data_Correlation);  

	CorrelationDenom[5] = Data_Correlation[0];
	Weight_CorrelationDenom[5] = Data_Correlation[1];

	if(fNumberSixth <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberSixth;

		this->Correlation(Number_JoinedCov, fg1,fg2,fg3,fg4,fg1,-1.*fg1,fg2,-1.*fg2,fg3,-1.*fg3,fg4,-1.*fg4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[5] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[5] = Data_Correlation[1];
	}

    } //if(fNumberSixth!=0)


    if(fNumberSeventh!=0)
    {
	this->Correlation(fNumberSeventh, fh1,fh2,fh3,fh4,fh5,fh6,fh7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

	CorrelationNum[6] = Data_Correlation[0];
	Weight_CorrelationNum[6] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberSeventh;

	this->Correlation(Number_Denom, fh1,-1.*fh1,fh2,-1.*fh2,fh3,-1.*fh3,fh4,-1.*fh4,fh5,-1.*fh5,fh6,-1.*fh6,fh7,-1.*fh7, Data_Correlation);  

	CorrelationDenom[6] = Data_Correlation[0];
	Weight_CorrelationDenom[6] = Data_Correlation[1];

	if(fNumberSeventh <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberSeventh;

		this->Correlation(Number_JoinedCov, fh1,fh2,fh3,fh4,fh1,-1.*fh1,fh2,-1.*fh2,fh3,-1.*fh3,fh4,-1.*fh4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[6] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[6] = Data_Correlation[1];
	}

    } //if(fNumberSeventh!=0)


    if(fNumberEighth!=0)
    {
	this->Correlation(fNumberEighth, fi1,fi2,fi3,fi4,fi5,fi6,fi7,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  

	CorrelationNum[7] = Data_Correlation[0];
	Weight_CorrelationNum[7] = Data_Correlation[1];

	Int_t Number_Denom = 2*fNumberEighth;

	this->Correlation(Number_Denom, fi1,-1.*fi1,fi2,-1.*fi2,fi3,-1.*fi3,fi4,-1.*fi4,fi5,-1.*fi5,fi6,-1.*fi6,fi7,-1.*fi7, Data_Correlation);  

	CorrelationDenom[7] = Data_Correlation[0];
	Weight_CorrelationDenom[7] = Data_Correlation[1];

	if(fNumberEighth <= 4.) //calculated the joined product of numerator and denominator as one single correlator. Only possible for order of correlation in numerator <= 4
	{
		Int_t Number_JoinedCov = 3*fNumberEighth;

		this->Correlation(Number_JoinedCov, fi1,fi2,fi3,fi4,fi1,-1.*fi1,fi2,-1.*fi2,fi3,-1.*fi3,fi4,-1.*fi4,0.,0., Data_Correlation);  

		CorrelationJoinedCov[7] = Data_Correlation[0];
		Weight_CorrelationJoinedCov[7] = Data_Correlation[1];
	}

    } //if(fNumberEighth!=0)


    for(Int_t i=0; i<8;i++)
   {
     fResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+0.5,CorrelationNum[i],Weight_CorrelationNum[i]); //safe output first set of harmonics
     fResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+1.5,CorrelationDenom[i],Weight_CorrelationDenom[i]); //safe output first set of harmonics

     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+0.5,CorrelationNum[i]*CorrelationDenom[i],Weight_CorrelationNum[i]*Weight_CorrelationDenom[i]); //w_D*N*w_D*D
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+1.5,Weight_CorrelationNum[i]*Weight_CorrelationDenom[i],1.); //w_N*w_D
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+2.5,Weight_CorrelationNum[i],1.); //w_N
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+3.5,Weight_CorrelationDenom[i],1.); //w_D 

     fJoinedCovResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+0.5,CorrelationJoinedCov[i],Weight_CorrelationJoinedCov[i]); //Joined Cov Term z //GANESHA
     fJoinedCovResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+1.5,Weight_CorrelationJoinedCov[i],1.); //w_z safe output first set of harmonics
	
   } 

	delete [] Data_Correlation; 

  } //if(fParticles>=fMinNumberPart)

} //void AliAnalysisTaskStudentsML::MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array)

//==========================================================================================================================================================================

  void AliAnalysisTaskStudentsML::MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)
  {

 if(Mixed_Mult_A>=fMinNumberPart && Mixed_Mult_B>=fMinNumberPart) 
 { 
    Double_t* Data_Correlation = new Double_t[3]; //cos, weight, sin

    Double_t FirstCorrelation=0.;
    Double_t Weight_FirstCorrelation=0.;
    Double_t SecondCorrelation=0.;
    Double_t Weight_SecondCorrelation=0.;

    //Dummy weights, may be implemented later
    Double_t* Dummy_Weights_A = new Double_t[Mixed_Mult_A]; 
    Double_t* Dummy_Weights_B = new Double_t[Mixed_Mult_B]; 
    for(Int_t i=0; i<Mixed_Mult_A; i++){Dummy_Weights_A[i]=1.;}
    for(Int_t i=0; i<Mixed_Mult_B; i++){Dummy_Weights_B[i]=1.;}

    // Calculate Q-vectors group A, then do correlation. Then same for B
    this->CalculateQvectors(Mixed_Mult_A, Mixed_Angle_A, Dummy_Weights_A);
    this->Correlation(4.,Harmonicus, -Harmonicus, Harmonicus, -Harmonicus,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., Data_Correlation);  //do the correlation for the first set

    FirstCorrelation=Data_Correlation[0];
    Weight_FirstCorrelation=Data_Correlation[1];

    //~~~~~~~~~~~~~~~~~
    this->CalculateQvectors(Mixed_Mult_B, Mixed_Angle_B, Dummy_Weights_B);
    this->Correlation(4.,Harmonicus, -Harmonicus, Harmonicus, -Harmonicus,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., Data_Correlation); 

    SecondCorrelation=Data_Correlation[0];
    Weight_SecondCorrelation=Data_Correlation[1];
    

    //~~~~~~~~~~~~~~~~~~
     Double_t Special_Weight = (Double_t)Mixed_Mult_A*((Double_t)Mixed_Mult_A-1.)*(Double_t)Mixed_Mult_B*((Double_t)Mixed_Mult_B-1.);
     TComplex Mixed = TComplex(0.,0.);
     Mixed = CalculateMixedQVectors((Double_t)Harmonicus, Mixed_Mult_A, Mixed_Mult_B, Mixed_Angle_A, Mixed_Angle_B);

    //~~~~~~~~~~~~~~~~~~

   fResults[MP_CentBin]->Fill(0.5,(1./Special_Weight)*Mixed.Re(),Special_Weight); //safe output first set of harmonics
   fResults[MP_CentBin]->Fill(1.5,FirstCorrelation*SecondCorrelation,Weight_FirstCorrelation*Weight_SecondCorrelation); //safe output second set of harmonics    
   
   fMixedParticleHarmonics[MP_CentBin]->Fill(0.5,FirstCorrelation,Weight_FirstCorrelation);
   fMixedParticleHarmonics[MP_CentBin]->Fill(1.5,SecondCorrelation,Weight_SecondCorrelation);

   delete [] Dummy_Weights_A; 
   delete [] Dummy_Weights_B;  
   delete [] Data_Correlation; 

 }  
 
} //void AliAnalysisTaskStudentsML::MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)

//==========================================================================================================================================================================




