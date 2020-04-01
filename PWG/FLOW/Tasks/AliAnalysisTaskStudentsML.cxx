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
 fControlHistogramsList(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fPhiHistBeforeTrackSeletion(NULL),
 fEtaHistBeforeTrackSeletion(NULL),
 fPTHistBeforeTrackSeletction(NULL),
 fPhiHistBeforeTrackSeletionSecond(NULL),
 fEtaHistBeforeTrackSeletionSecond(NULL),
 fPTHistBeforeTrackSeletctionSecond(NULL),
 fTotalMultBeforeTrackSeletion(NULL),
 fMultiHistoBeforeTrackSeletion(NULL),
 fPhiHistAfterTrackSeletion(NULL),
 fEtaHistAfterTrackSeletion(NULL),
 fPTHistAfterTrackSeletction(NULL),
 fPhiHistAfterTrackSeletionSecond(NULL),
 fEtaHistAfterTrackSeletionSecond(NULL),
 fPTHistAfterTrackSeletctionSecond(NULL),
 fTotalMultAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion_Second(NULL),
 fMultiHistoBeforeMultCut(NULL),
 fTPCClustersBeforeCut(NULL),
 fTPCClustersAfterCut(NULL),
 fITSClustersBeforeCut(NULL),
 fITSClustersAfterCut(NULL),
 fChiSquareTPCBeforeCut(NULL),
 fChiSquareTPCAfterCut(NULL),
 fDCAzBeforeCut(NULL),
 fDCAzAfterCut(NULL),
 fDCAxyBeforeCut(NULL),
 fDCAxyAfterCut(NULL),
 //SelectionCuts
 bDoAnalysis(kTRUE),
 bUseRecoKineTable(kTRUE),
 bMultCut(kFALSE),
 fMainFilter(0),
 fSecondFilter(0),
 fSlopeUpperLine(0.),
 fAxisUpperLine(0.),
 fSlopeLowerLine(0.),
 fAxisLowerLine(0.),
 fMinCentrality(0.),
 fMaxCentrality(100.),
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
  fVertexXBefore(NULL),
 fVertexXAfter(NULL),
 fVertexYBefore(NULL),
 fVertexYAfter(NULL),
 fVertexZBefore(NULL),
 fVertexZAfter(NULL),
 fCentralityfromVZero(kTRUE),
 //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 bNumberTPCCluster(kTRUE),
 bNumberITSCluster(kTRUE),
 bChiSquareTPC(kTRUE),
 bDCAz(kTRUE),
 bDCAxy(kTRUE),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fMinTPCCluster(70.),
 fMinITSCluster(2.),
 fMinChiSquareTPC(0.1),
 fMaxChiSquareTPC(4.0),
 fMaxDCAz(3.2),
 fMaxDCAxy(2.4),
 //Weights
 bUseWeights(kFALSE),
 bUsePtWeights(kFALSE),
 bUsePhiWeights(kFALSE), 
 bUseEtaWeights(kFALSE),
 bGridWeights(kTRUE),
 fPeriodUsedForWeight(""),
 //Variables for the correlation
 fMaxCorrelator(12),
 fNumber(6),  //number of correlation first correlator
 fNumberSecond(6), //number of correlation second correlator
 fNumberThird(6),
 bDoThirdCorrelation(kFALSE),
 fMinNumberPart(10),
 bUseRatioWeight(kTRUE),
 fDenominatorMinValue(1.0e-16),
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), fh8(0), fh9(0), fh10(0), fh11(0), fh12(0),  //harmonics
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), fa8(0), fa9(0), fa10(0), fa11(0), fa12(0), //second set of harmonics
 fb1(0), fb2(0), fb3(0), fb4(0), fb5(0), fb6(0), fb7(0), fb8(0), fb9(0), fb10(0), fb11(0), fb12(0), //third set of harmonics
 fCentrality(NULL),
 fCentralitySecond(NULL),
 fCentralityThird(NULL),
 fEvCentrality(NULL),
 bDoEbERatio(kFALSE),
 fMixedParticleHarmonics(NULL),
 bDoMixed(kFALSE),
 bDifferentCharge(kTRUE),
 bSetSameChargePositiv(kTRUE),			
 fMixedHarmonic(0),
 fCounterHistogram(NULL),
 // Final results:
 fFinalResultsList(NULL)
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
 fControlHistogramsList(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fPhiHistBeforeTrackSeletion(NULL),
 fEtaHistBeforeTrackSeletion(NULL),
 fPTHistBeforeTrackSeletction(NULL),
 fPhiHistBeforeTrackSeletionSecond(NULL),
 fEtaHistBeforeTrackSeletionSecond(NULL),
 fPTHistBeforeTrackSeletctionSecond(NULL),
 fTotalMultBeforeTrackSeletion(NULL),
 fMultiHistoBeforeTrackSeletion(NULL),
 fPhiHistAfterTrackSeletion(NULL),
 fEtaHistAfterTrackSeletion(NULL),
 fPTHistAfterTrackSeletction(NULL),
 fPhiHistAfterTrackSeletionSecond(NULL),
 fEtaHistAfterTrackSeletionSecond(NULL),
 fPTHistAfterTrackSeletctionSecond(NULL),
 fTotalMultAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion_Second(NULL),
 fMultiHistoBeforeMultCut(NULL),
 fTPCClustersBeforeCut(NULL),
 fTPCClustersAfterCut(NULL),
 fITSClustersBeforeCut(NULL),
 fITSClustersAfterCut(NULL),
 fChiSquareTPCBeforeCut(NULL),
 fChiSquareTPCAfterCut(NULL),
 fDCAzBeforeCut(NULL),
 fDCAzAfterCut(NULL),
 fDCAxyBeforeCut(NULL),
 fDCAxyAfterCut(NULL),
 //SelectionCuts
 bDoAnalysis(kTRUE),
 bUseRecoKineTable(kTRUE),
 bMultCut(kFALSE),
 fMainFilter(0),
 fSecondFilter(0),
 fSlopeUpperLine(0.),
 fAxisUpperLine(0.),
 fSlopeLowerLine(0.),
 fAxisLowerLine(0.),
 fMinCentrality(0.),
 fMaxCentrality(100.),
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
 fVertexXBefore(NULL),
 fVertexXAfter(NULL),
 fVertexYBefore(NULL),
 fVertexYAfter(NULL),
 fVertexZBefore(NULL),
 fVertexZAfter(NULL),
 fCentralityfromVZero(kTRUE),
 //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 bNumberTPCCluster(kTRUE),
 bNumberITSCluster(kTRUE),
 bChiSquareTPC(kTRUE),
 bDCAz(kTRUE),
 bDCAxy(kTRUE),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fMinTPCCluster(70.),
 fMinITSCluster(2.),
 fMinChiSquareTPC(0.1),
 fMaxChiSquareTPC(4.0),
 fMaxDCAz(3.2),
 fMaxDCAxy(2.4),
 //Weights
 bUseWeights(kFALSE),
 bUsePtWeights(kFALSE),
 bUsePhiWeights(kFALSE), 
 bUseEtaWeights(kFALSE),
 bGridWeights(kTRUE),
 fPeriodUsedForWeight(""),
 //Variables for the correlation
 fMaxCorrelator(12),
 fNumber(6),  //number of correlation first correlator
 fNumberSecond(6), //number of correlation second correlator
 fNumberThird(6),
 bDoThirdCorrelation(kFALSE),
 fMinNumberPart(10),
 bUseRatioWeight(kTRUE),
 fDenominatorMinValue(1.0e-16),
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), fh8(0), fh9(0), fh10(0), fh11(0), fh12(0),  //harmonics
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), fa8(0), fa9(0), fa10(0), fa11(0), fa12(0), //second set of harmonics
 fb1(0), fb2(0), fb3(0), fb4(0), fb5(0), fb6(0), fb7(0), fb8(0), fb9(0), fb10(0), fb11(0), fb12(0), //third set of harmonics
 // Final results:
 fCentrality(NULL),
 fCentralitySecond(NULL),
 fCentralityThird(NULL),
 fEvCentrality(NULL),
 bDoEbERatio(kFALSE),
 fMixedParticleHarmonics(NULL),
 bDoMixed(kFALSE),
 bDifferentCharge(kTRUE),
 bSetSameChargePositiv(kTRUE),		
 fMixedHarmonic(0),
 fCounterHistogram(NULL),
 fFinalResultsList(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML()");

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
 this->BookControlHistograms();
 this->BookFinalResultsHistograms();

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
	
 fVertexXAfter->Fill(avtx->GetX());
 fVertexYAfter->Fill(avtx->GetY());
 fVertexZAfter->Fill(avtx->GetZ());
 

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 //Access Data
 
 //b.0) Start analysis over AODs
 Int_t nTracks = 0;  // number of all tracks in current event 
 nTracks = aAODEvent->GetNumberOfTracks();

 Double_t* angles_A = NULL; // Azimuthal angles
 Double_t* pt_A = NULL;
 Double_t* eta_A = NULL;
 Double_t* weights_A = NULL;
 Int_t Multi_Ang_A = 0.;

 Double_t* angles_B = NULL; 
 Int_t Multi_Ang_B = 0.; 

 Int_t CounterSameCharge = 0.; //used if bDifferentCharge = kTRUE

 if(nTracks>0){fMultiHistoBeforeTrackSeletion->Fill(nTracks);} //multiplicity distribution before track selection
 for(Int_t u=0;u<nTracks;u++){fTotalMultBeforeTrackSeletion->Fill(0.5);} //total number of particles in whole centrality class before track sel.

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
     	Double_t ChiSquareInTPC = (aAODTrack->GetTPCchi2())/(aAODTrack->GetNcls(1)); //chi square in the TPC
     	Double_t ValueDCAz = aAODTrack->ZAtDCA();  //z-coordinate of DCA
        Double_t ValueDCAy = aAODTrack->YAtDCA();  //x-coordinate of DCA
        Double_t ValueDCAx = aAODTrack->XAtDCA();  //y-coordinate of DCA
	Double_t ValueDCAxy = TMath::Sqrt(ValueDCAx*ValueDCAx + ValueDCAy*ValueDCAy);
 

  	//............................................................................................
	//Fill control histograms with the particles before track selection:

	if(!bDoMixed)
	{
          fPhiHistBeforeTrackSeletion->Fill(phi); 
          fEtaHistBeforeTrackSeletion->Fill(eta);
          fPTHistBeforeTrackSeletction->Fill(pt);
	}//if(!bDoMixed)

        if(bDoMixed)
	{
          if(charge>0.)
	  {
	    fPhiHistBeforeTrackSeletion->Fill(phi); 
            fEtaHistBeforeTrackSeletion->Fill(eta);
            fPTHistBeforeTrackSeletction->Fill(pt);
          }//if(charge>0.)

	  if(charge<0.)
	  {
	    fPhiHistBeforeTrackSeletionSecond->Fill(phi); 
            fEtaHistBeforeTrackSeletionSecond->Fill(eta);
            fPTHistBeforeTrackSeletctionSecond->Fill(pt);
	  } //if(charge<0.)
	} //if(bDoMixed)
    

     	fTPCClustersBeforeCut->Fill(NumberOfTPCClusters);
     	fITSClustersBeforeCut->Fill(NumberOfITSClusters);
    	fChiSquareTPCBeforeCut->Fill(ChiSquareInTPC);
     	fDCAzBeforeCut->Fill(ValueDCAz);
	fDCAxyBeforeCut->Fill(ValueDCAxy);

	//............................................................................................
	//Track did not pass physics selection 
      	if(!TrackSelection(aAODTrack)){continue;} 

	//............................................................................................
     	// Fill control histograms with the particles after track selection:

	if(!bDoMixed)
	{
        fPhiHistAfterTrackSeletion->Fill(phi); 
        fEtaHistAfterTrackSeletion->Fill(eta);
        fPTHistAfterTrackSeletction->Fill(pt);
	}//if(!bDoMixed)

        if(bDoMixed)
	{
          if(charge>0.)
	  {
	    fPhiHistAfterTrackSeletion->Fill(phi); 
            fEtaHistAfterTrackSeletion->Fill(eta);
            fPTHistAfterTrackSeletction->Fill(pt);
          }//if(charge>0.)
	  if(charge<0.)
	  {
	    fPhiHistAfterTrackSeletionSecond->Fill(phi); 
            fEtaHistAfterTrackSeletionSecond->Fill(eta);
            fPTHistAfterTrackSeletctionSecond->Fill(pt);
	  }//if(charge<0.)
        }//if(bDoMixed)

     	fTPCClustersAfterCut->Fill(NumberOfTPCClusters);
     	fITSClustersAfterCut->Fill(NumberOfITSClusters);
     	fChiSquareTPCAfterCut->Fill(ChiSquareInTPC);
     	fDCAzAfterCut->Fill(ValueDCAz);
     	fDCAxyAfterCut->Fill(ValueDCAxy);

	//............................................................................................

     	if(!bDoMixed) {Multi_Ang_A+=1.;}
    	if(bDoMixed)
     	{
		if(bDifferentCharge)
		{
			if(charge>0.){Multi_Ang_A+=1.;}
			if(charge<0.){Multi_Ang_B+=1.;}
		}//if(bDifferentCharge)

		if(!bDifferentCharge)
		{
			Double_t UsedCharge = 0.;
			if(bSetSameChargePositiv) {UsedCharge = charge;}
			if(!bSetSameChargePositiv) {UsedCharge = -1.*charge;}

			if(UsedCharge>0.)
			{
		  	   if(0 == CounterSameCharge%2) {Multi_Ang_A+=1.; CounterSameCharge++; }
		  	   else {Multi_Ang_B+=1.; CounterSameCharge++; }
			}

		}//if(!bDifferentCharge)
     	}//if(bDoMixed)


 } // Frist Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut)



 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Declare arrays for Phi

 angles_A = new Double_t[Multi_Ang_A]; 
 pt_A = new Double_t[Multi_Ang_A];
 eta_A = new Double_t[Multi_Ang_A];
 weights_A = new Double_t[Multi_Ang_A];

 Multi_Ang_A = 0.; //Reset for filling
       
 if(Multi_Ang_B>0){angles_B = new Double_t[Multi_Ang_B];}
 else{angles_B = new Double_t[1]; angles_B[0]=0.;} //Dummy 
 Multi_Ang_B = 0.; //Reset for filling

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
    
     if(!TrackSelection(aAODTrack)){continue;} //Track did not pass physics selection 
     

     if(!bDoMixed)
     {
	angles_A[Multi_Ang_A] = phi; 
	pt_A[Multi_Ang_A] = pt; 
	eta_A[Multi_Ang_A] = eta; 
	weights_A[Multi_Ang_A] = 1.;
        Multi_Ang_A+=1.;
     }//if(!bDoMixed)

     if(bDoMixed)
     {
	if(bDifferentCharge)
	{	
		if(charge>0.){angles_A[Multi_Ang_A] = phi; Multi_Ang_A+=1.;}
		if(charge<0.){angles_B[Multi_Ang_B] = phi; Multi_Ang_B+=1.;}
	}//if(bDifferentCharge)

	if(!bDifferentCharge)
	{
		Double_t UsedCharge = 0.;
		if(bSetSameChargePositiv) {UsedCharge = charge;}
		if(!bSetSameChargePositiv) {UsedCharge = -1.*charge;}

		if(UsedCharge>0.)
		{
		   if(0 == CounterSameCharge%2) {angles_A[Multi_Ang_A] = phi; Multi_Ang_A+=1.; CounterSameCharge++; }
		   else {angles_B[Multi_Ang_B] = phi; Multi_Ang_B+=1.; CounterSameCharge++;}
		}

	}//if(!bDifferentCharge)

     }//if(bDoMixed)

 } //Second Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut) 



 //Add Member Function Get Weights
 if(bUseWeights)
 {
  Int_t CallRunNumber = aAODEvent->GetRunNumber();
  CalculateWeight(CallRunNumber, weights_A, Multi_Ang_A, angles_A, pt_A, eta_A);
 }


 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //c) Calculus
 
 if(!bDoMixed)
 {
   if(Multi_Ang_A>0){fMultiHistoAfterTrackSeletion->Fill(Multi_Ang_A);} //multiplicity distribution after track selection
   for(Int_t u=0;u<Multi_Ang_A;u++){fTotalMultAfterTrackSeletion->Fill(0.5);} //total number of particles in whole centrality class after track sel.

   this->MainTask(Multi_Ang_A, angles_A,weights_A); //Actual Multi-Particle Correlation

 }//if(!bDoMixed)
 
 if(bDoMixed)
 {
   if(Multi_Ang_A>0){fMultiHistoAfterTrackSeletion->Fill(Multi_Ang_A);} //multiplicity distribution after track selection
   for(Int_t u=0;u<Multi_Ang_A+Multi_Ang_B;u++){fTotalMultAfterTrackSeletion->Fill(0.5);} //total number of particles in whole centrality class after track sel.
   if(Multi_Ang_B>0){fMultiHistoAfterTrackSeletion_Second->Fill(Multi_Ang_B);}

   this->MixedParticle(fMixedHarmonic, Multi_Ang_A, angles_A, Multi_Ang_B, angles_B);

 }//if(bDoMixed)
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // d) Reset event-by-event objects:
 Multi_Ang_A =0.;
 delete [] angles_A; 
 delete [] pt_A;
 delete [] eta_A;
 delete [] weights_A;
 Multi_Ang_B =0.;
 delete [] angles_B;
 CounterSameCharge = 0.; //reset the same charge counter



} //void AliAnalysisTaskStudentsML::PhysicsAnalysis()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::GetKineDist(AliAODEvent *aAODEve, AliMCEvent *aMCEve)
{ 
 //Used for gaining weights
 //a) Global QA for Reco and Kine
 //b) Run over Reco + Generate Table 
 //c) Run over Kine + Get Distributions

 //a) Global QA for Reco and Kine
 if(!GlobalQualityAssurance(aAODEve)){return;}

 if(!GlobalQualityAssurance(aMCEve)){return;} 
 
 AliMCVertex *avtx = (AliMCVertex*)aMCEve->GetPrimaryVertex();
	
 fVertexXAfter->Fill(avtx->GetX());
 fVertexYAfter->Fill(avtx->GetY());
 fVertexZAfter->Fill(avtx->GetZ());

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
      	if(!TrackSelection(aAODTrack))
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
	if(!isPrimary && !isWeakSecondary){continue;}

	fPhiHistBeforeTrackSeletion->Fill(phi); 
        fEtaHistBeforeTrackSeletion->Fill(eta);
        fPTHistBeforeTrackSeletction->Fill(pt);

	//........................................................................
	//Track did not pass physics selection 
      	if(!TrackSelection(aMCTrack)){continue;} 

	if(bUseRecoKineTable)
	{

		if(RecoKineMap->GetValue(iTrack)==-1){continue;} // GetValue = -1: the corresponding reco track has not been selected by its track selection.
		else
		{	// GetValue = 0: the kine track has been lost in ALICE --> Included in the pT distribution.
        		// GetValue = 1: the reco track has been selected.
			fPhiHistAfterTrackSeletion->Fill(phi); 
        		fEtaHistAfterTrackSeletion->Fill(eta);
        		fPTHistAfterTrackSeletction->Fill(pt);
		} //else
	} //if(bUseRecoKineTable)
	else
	{
 		fPhiHistAfterTrackSeletion->Fill(phi); 
        	fEtaHistAfterTrackSeletion->Fill(eta);
        	fPTHistAfterTrackSeletction->Fill(pt);
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

 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskStudentsML::Terminate(Option_t *)

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

   for(Int_t cs=0;cs<2;cs++) 
   {
     for(Int_t c=0;c<fMaxCorrelator;c++)
     {
   
      fRecursion[cs][c] = NULL; //! [cs]: real (0) or imaginary part (1) ....
   
     }  
    }  //for(Int_t cs=0;cs<2;cs++)

   for(Int_t js=0;js<97;js++) 
   {
     for(Int_t j=0;j<13;j++)
     {
   
      fQvector[js][j] = TComplex(0.,0.); //! 
   
     } 
   } 


  for(Int_t i=0; i<3; i++)
  {
	fWeightsHist[i] = NULL;
  }
   
} // void AliAnalysisTaskStudentsML::InitializeArrays()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskStudentsML::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);


} // void AliAnalysisTaskStudentsML::BookAndNestAllLists()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra
 // b) Book histogram to hold phi distribution before track selection
 // c) Book histogram to hold eta distribution before track selection
 // d) Book Mult. Histo before before track selection
 // e) Book histogam to hold multiplicty distribution before track selection
 // f) Book histogram to hold phi distribution after track selection
 // g) Book histogram to hold eta distribution after track selection:
 // h) Book Mult. Histo before after track selection
 // i) Book histogam to hold multiplicty distribution after track selection:
 // j) Book histogam to hold multiplicty distribution before high multiplicity outlier cut:
 // k) Book histogam for Vertex X before Cut
 // l) Book histogam for Vertex X after Cut
 // m) Book histogam for Vertex Y before Cut
 // n) Book histogam for Vertex Y after Cut
 // o) Book histogam for Vertex Z before Cut
 // p) Book histogam for Vertex Z after Cut
 // q) Book histogram to debug
 // r) Book histogram for number of TPC clustes before cut
 // s) Book histogram for number of TPC clustes after cut
 // t) Book histogram for number of ITC clusters before the cut
 // u) Book histogram for number of ITC clusters after the cut
 // v) Book histogram for chi square TPC before cut
 // w) Book histogram for chi square TPC after cut
 // x) Book histogram for DCAz before + after cut
 // y) Book histogram for DCAxy before + after cut



 // a) Book histogram to hold pt spectra:
 fPTHistBeforeTrackSeletction = new TH1F("fPTHistBeforeTrackSeletction","Pt Distribution",1000,0.,10.);
 fPTHistBeforeTrackSeletction->GetXaxis()->SetTitle("P_t");
 fPTHistBeforeTrackSeletction->SetLineColor(4);
 fControlHistogramsList->Add(fPTHistBeforeTrackSeletction);

 fPTHistBeforeTrackSeletctionSecond = new TH1F("fPTHistBeforeTrackSeletctionSecond","Pt Distribution",1000,0.,10.);
 fPTHistBeforeTrackSeletctionSecond->GetXaxis()->SetTitle("P_t");
 fPTHistBeforeTrackSeletctionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fPTHistBeforeTrackSeletctionSecond);

 fPTHistAfterTrackSeletction = new TH1F("fPTHistAfterTrackSeletction","Pt Distribution",1000,0.,10.);
 fPTHistAfterTrackSeletction->GetXaxis()->SetTitle("P_t");
 fPTHistAfterTrackSeletction->SetLineColor(4);
 fControlHistogramsList->Add(fPTHistAfterTrackSeletction);

 fPTHistAfterTrackSeletctionSecond = new TH1F("fPTHistAfterTrackSeletctionSecond","Pt Distribution",1000,0.,10.);
 fPTHistAfterTrackSeletctionSecond->GetXaxis()->SetTitle("P_t");
 fPTHistAfterTrackSeletctionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fPTHistAfterTrackSeletctionSecond);


 // b) Book histogram to hold phi distribution before track selection:
 fPhiHistBeforeTrackSeletion = new TH1F("fPhiHistBeforeTrackSeletion","Phi Distribution",1000,0.,6.3);
 fPhiHistBeforeTrackSeletion->GetXaxis()->SetTitle("Phi");
 fPhiHistBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistBeforeTrackSeletion);

 fPhiHistBeforeTrackSeletionSecond = new TH1F("fPhiHistBeforeTrackSeletionSecond","Phi Distribution",1000,0.,6.3);
 fPhiHistBeforeTrackSeletionSecond->GetXaxis()->SetTitle("Phi");
 fPhiHistBeforeTrackSeletionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistBeforeTrackSeletionSecond);

 // c) Book histogram to hold eta distribution before track selection:
 fEtaHistBeforeTrackSeletion = new TH1F("fEtaHistBeforeTrackSeletion","Eta Distribution",1000,-1.,1.);
 fEtaHistBeforeTrackSeletion->GetXaxis()->SetTitle("Eta");
 fEtaHistBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistBeforeTrackSeletion);

 fEtaHistBeforeTrackSeletionSecond = new TH1F("fEtaHistBeforeTrackSeletionSecond","Eta Distribution",1000,-1.,1.);
 fEtaHistBeforeTrackSeletionSecond->GetXaxis()->SetTitle("Eta");
 fEtaHistBeforeTrackSeletionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistBeforeTrackSeletionSecond);

 // d) Book Mult. Histo before before track selection
 fTotalMultBeforeTrackSeletion = new TH1F("fTotalMultBeforeTrackSeletion","Mult. Counts per Class before brute cut",1,0.,1.);
 fTotalMultBeforeTrackSeletion->GetYaxis()->SetTitle("Counts");
 fTotalMultBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fTotalMultBeforeTrackSeletion);
 
 // e) Book histogam to hold multiplicty distribution before track selection:
 fMultiHistoBeforeTrackSeletion = new TH1F("fMultiHistoBeforeTrackSeletion","Multiplicity",5000,0.,5000.); 
 fMultiHistoBeforeTrackSeletion->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoBeforeTrackSeletion);

 // f) Book histogram to hold phi distribution after track selection:
 fPhiHistAfterTrackSeletion = new TH1F("fPhiHistAfterTrackSeletion","Phi Distribution",1000,0.,6.3);
 fPhiHistAfterTrackSeletion->GetXaxis()->SetTitle("Phi");
 fPhiHistAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistAfterTrackSeletion);

 fPhiHistAfterTrackSeletionSecond = new TH1F("fPhiHistAfterTrackSeletionSecond","Phi Distribution",1000,0.,6.3);
 fPhiHistAfterTrackSeletionSecond->GetXaxis()->SetTitle("Phi");
 fPhiHistAfterTrackSeletionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistAfterTrackSeletionSecond);

 // g) Book histogram to hold eta distribution after track selection:
 fEtaHistAfterTrackSeletion = new TH1F("fEtaHistAfterTrackSeletion","Eta Distribution",1000,-1.,1.);
 fEtaHistAfterTrackSeletion->GetXaxis()->SetTitle("Eta");
 fEtaHistAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistAfterTrackSeletion);

 fEtaHistAfterTrackSeletionSecond = new TH1F("fEtaHistAfterTrackSeletionSecond","Eta Distribution",1000,-1.,1.);
 fEtaHistAfterTrackSeletionSecond->GetXaxis()->SetTitle("Eta");
 fEtaHistAfterTrackSeletionSecond->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistAfterTrackSeletionSecond);

 // h) Book Mult. Histo before after track selection
  fTotalMultAfterTrackSeletion = new TH1F("fTotalMultAfterTrackSeletion","Mult. Counts per Class before brute cut",1,0.,1.);
 fTotalMultAfterTrackSeletion->GetYaxis()->SetTitle("Counts");
 fTotalMultAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fTotalMultAfterTrackSeletion);
 
 // i) Book histogam to hold multiplicty distribution after track selection:
 fMultiHistoAfterTrackSeletion = new TH1F("fMultiHistoAfterTrackSeletion","Multiplicity",5000,0.,5000.);
 fMultiHistoAfterTrackSeletion->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoAfterTrackSeletion);

 fMultiHistoAfterTrackSeletion_Second = new TH1F("fMultiHistoAfterTrackSeletion_Second","Multiplicity",5000,0.,5000.);
 fMultiHistoAfterTrackSeletion_Second->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoAfterTrackSeletion_Second);

 // j) Book histogam to hold multiplicty distribution before high multiplicity outlier cut:
 fMultiHistoBeforeMultCut = new TH1F("fMultiHistoBeforeMultCut","Multiplicity",5000,0.,5000.); 
 fMultiHistoBeforeMultCut->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoBeforeMultCut);

 // k) Book histogam for Vertex X before Cut
 fVertexXBefore = new TH1F("fVertexXBefore","fVertexXBefore",1000,-20.,20.); 
 fVertexXBefore->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexXBefore);

  // l) Book histogam for Vertex X after Cut
 fVertexXAfter = new TH1F("fVertexXAfter","fVertexXAfter",1000,-20.,20.); 
 fVertexXAfter->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexXAfter);

 // m) Book histogam for Vertex Y before Cut
 fVertexYBefore = new TH1F("fVertexYBefore","fVertexYBefore",1000,-20.,20.); 
 fVertexYBefore->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexYBefore);

 // n) Book histogam for Vertex Y after Cut
 fVertexYAfter = new TH1F("fVertexYAfter","fVertexYAfter",1000,-20.,20.); 
 fVertexYAfter->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexYAfter);

  // o) Book histogam for Vertex Y after Cut
 fVertexZBefore = new TH1F("fVertexZBefore","fVertexZBefore",1000,-20.,20.); 
 fVertexZBefore->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexZBefore);

  // p) Book histogam for Vertex Y after Cut
 fVertexZAfter = new TH1F("fVertexZAfter","fVertexZAfter",1000,-20.,20.); 
 fVertexZAfter->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexZAfter);

 // q) Book histogram to debug
 fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks",3,0.,3.); 
 fControlHistogramsList->Add(fCounterHistogram);

 // r) Book histogram for number of TPC clustes before cut
 fTPCClustersBeforeCut = new TH1F("fTPCClustersBeforeCut","fTPCClustersBeforeCut",170,0.,170.); 
 fControlHistogramsList->Add(fTPCClustersBeforeCut);

 // s) Book histogram for number of TPC clustes after cut
 fTPCClustersAfterCut = new TH1F("fTPCClustersAfterCut","fTPCClustersAfterCut",170,0.,170.); 
 fControlHistogramsList->Add(fTPCClustersAfterCut);

 // t) Book histogram for number of ITC clusters before the cut
 fITSClustersBeforeCut = new TH1F("fITSClustersBeforeCut","fITSClustersBeforeCut",10,0.,10.); 
 fControlHistogramsList->Add(fITSClustersBeforeCut);

 // u) Book histogram for number of ITC clusters after the cut
 fITSClustersAfterCut = new TH1F("fITSClustersAfterCut","fITSClustersAfterCut",10,0.,10.); 
 fControlHistogramsList->Add(fITSClustersAfterCut);

 // v) Book histogram for chi square TPC before cut
 fChiSquareTPCBeforeCut = new TH1F("fChiSquareTPCBeforeCut","fChiSquareTPCBeforeCut",1000,0.,20.); 
 fControlHistogramsList->Add(fChiSquareTPCBeforeCut);

 // w) Book histogram for chi square TPC after cut
 fChiSquareTPCAfterCut = new TH1F("fChiSquareTPCAfterCut","fChiSquareTPCAfterCut",1000,0.,20.); 
 fControlHistogramsList->Add(fChiSquareTPCAfterCut);

 // x) Book histogram for DCAz before + after cut
 fDCAzBeforeCut = new TH1F("fDCAzBeforeCut","fDCAzBeforeCut",1000,0.,10.); 
 fControlHistogramsList->Add(fDCAzBeforeCut);

 fDCAzAfterCut = new TH1F("fDCAzAfterCut","fDCAzAfterCut",1000,0.,10.); 
 fControlHistogramsList->Add(fDCAzAfterCut);
 // y) Book histogram for DCAxy before + after cut
 fDCAxyBeforeCut = new TH1F("fDCAxyBeforeCut","fDCAxyBeforeCut",1000,0.,10.); 
 fControlHistogramsList->Add(fDCAxyBeforeCut);

 fDCAxyAfterCut = new TH1F("fDCAxyAfterCut","fDCAxyAfterCut",1000,0.,10.); 
 fControlHistogramsList->Add(fDCAxyAfterCut);

} //void AliAnalysisTaskStudentsML::BookControlHistograms()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.
  
 fCentrality = new TProfile("fCentrality","Result Analysis First Set Correlators",2,0.,2.); //centrality dependet output
 fCentrality->GetXaxis()->SetTitle("");
 fCentrality->GetYaxis()->SetTitle("flow");
 fCentrality->Sumw2();
 fFinalResultsList->Add(fCentrality);

 fCentralitySecond = new TProfile("fCentralitySecond","Result Analysis Second Set Correlators",2,0.,2.); //centrality dependet output
 fCentralitySecond->GetXaxis()->SetTitle("");
 fCentralitySecond->GetYaxis()->SetTitle("flow");
 fCentralitySecond->Sumw2(); 
 fFinalResultsList->Add(fCentralitySecond);

 fCentralityThird = new TProfile("fCentralityThird","Result Analysis Second Third Correlators",2,0.,2.); //centrality dependet output
 fCentralityThird->GetXaxis()->SetTitle("");
 fCentralityThird->GetYaxis()->SetTitle("flow");
 fCentralityThird->Sumw2(); 
 fFinalResultsList->Add(fCentralityThird);

 fEvCentrality = new TProfile("fEvCentrality","Result Analysis EbE Method",1,0.,1.); //centrality dependet output
 fEvCentrality->GetXaxis()->SetTitle("");
 fEvCentrality->Sumw2();  
 fFinalResultsList->Add(fEvCentrality);

 fMixedParticleHarmonics = new TProfile("fMixedParticleHarmonics","fMixedParticleHarmonics",2,0.,2.); //centrality dependet output
 fMixedParticleHarmonics->GetXaxis()->SetTitle("");
 fMixedParticleHarmonics->GetYaxis()->SetTitle("");
 fMixedParticleHarmonics->Sumw2(); 
 fFinalResultsList->Add(fMixedParticleHarmonics);


 Cosmetics();
 
} // void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::Cosmetics()
{
 // Book everything here.
  
 for(Int_t cs=0;cs<2;cs++) 
 {
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   
   fRecursion[cs][c] = new TProfile("","",2,0.,2.); 
   //fRecursion[cs][c]->Sumw2();
 
   //NOTE for fRecursion: 1.) [cs] will say if its the real (0) or imaginary part (1) 
   // 2.) [c] gives gives the kind of correlation. [n] is the (n+2)-particle correlation
   //3.) first bin holds value of correlation (single event). Second bin holds the weight! 
   
  } // end of for(Int_t c=0;c<fMaxCorrelator;c++) 
 } // end of for(Int_t cs=0;cs<2;cs++) 

} // void Cosmetics()

//==========================================================================================================================================================================

Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliAODEvent *aAODevent){

  //a) Protection against NULL-Pointers
  //b) Check Centrality
  //c) Cuts on AliAODVertex:
  //d) remove high multiplicity outliers


  //a) Protection against NULL-Pointers
  if(!aAODevent){return kFALSE;}
  fCounterHistogram->Fill(1.5); // counter hist 2nd bin

  //b) Check Centrality
  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if(!ams){return kFALSE;}
  fCounterHistogram->Fill(2.5); // counter hist 3rd bin
 
  if(fCentralityfromVZero)
  {
  	if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality){ }
  	else{ return kFALSE; } // this event do not belong to the centrality class specified for this particular analysis 
  }

  if(!fCentralityfromVZero)
  {
  	if(ams->GetMultiplicityPercentile("CL1") >= fMinCentrality && ams->GetMultiplicityPercentile("CL1") < fMaxCentrality){ }
  	else{ return kFALSE; } // this event do not belong to the centrality class specified for this particular analysis 
  }
 

  // c) Cuts on AliAODVertex:
  AliAODVertex *avtx = (AliAODVertex*)aAODevent->GetPrimaryVertex();
 
  fVertexXBefore->Fill(avtx->GetX());
  fVertexYBefore->Fill(avtx->GetY());
  fVertexZBefore->Fill(avtx->GetZ());


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

	fMultiHistoBeforeMultCut->Fill(nTracks); //multiplicity distribution before high multiplicity outlier removal
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

 Bool_t AliAnalysisTaskStudentsML::TrackSelection(AliAODTrack *aTrack)
 {

        // example variables for each track:
 	/*Double_t px = aTrack->Px(); // x-component of momenta
 	Double_t py = aTrack->Py(); // y-component of momenta
 	Double_t pz = aTrack->Pz(); // z-component of momenta
 	Double_t e = aTrack->E();  // energy
        Double_t charge = aTrack->Charge(); // charge
 	Double_t phi = aTrack->Phi(); // azimuthal angle*/
 	Double_t eta = aTrack->Eta(); // pseudorapidity
 	Double_t pt = aTrack->Pt(); // Pt (transverse momentum)
	Int_t NumberOfTPCClusters = aTrack->GetTPCNcls(); //number of TPC clusters of the track
	Int_t NumberOfITSClusters = aTrack->GetITSNcls(); //number of ITS clusters of the track
	Double_t ChiSquareInTPC = (aTrack->GetTPCchi2())/(aTrack->GetNcls(1)); //chi square in the TPC
	Double_t ValueDCAz = aTrack->ZAtDCA();  //z-coordinate of DCA
	Double_t ValueDCAy = aTrack->YAtDCA();  //x-coordinate of DCA
        Double_t ValueDCAx = aTrack->XAtDCA();  //y-coordinate of DCA
	Double_t ValueDCAxy = TMath::Sqrt(ValueDCAx*ValueDCAx + ValueDCAy*ValueDCAy);


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
	  if(ValueDCAz>fMaxDCAz) return kFALSE;
	}

	if(bDCAxy) 
	{
	  if(ValueDCAxy>fMaxDCAxy) return kFALSE;
	}

    return kTRUE;

 }// end AliAnalysisTaskStudentsML::PhysicsSelection()

//==========================================================================================================================================================================
// Overloading functions for MC Events

 Bool_t AliAnalysisTaskStudentsML::GlobalQualityAssurance(AliMCEvent *aMCKineEvent)
 {

  //a) Reset Histo's
  //b) Protection against NULL-Pointers
  //c) Cuts on AliAODVertex:
  //d) remove high multiplicity outliers

  
  //a) Reset Histo's (are already filled in QA for AOD)
  fCounterHistogram->SetBinContent(2,0);
  fCounterHistogram->SetBinContent(3,0);
  fVertexXBefore->Reset();
  fVertexYBefore->Reset();
  fVertexZBefore->Reset();
  fMultiHistoBeforeMultCut->Reset();

  //b) Protection against NULL-Pointers
  if(!aMCKineEvent){return kFALSE;}
  fCounterHistogram->Fill(1.5); // counter hist 2nd bin

  // c) Cuts on AliAODVertex:
  AliMCVertex *avtx = (AliMCVertex*)aMCKineEvent->GetPrimaryVertex();
 
  fVertexXBefore->Fill(avtx->GetX());
  fVertexYBefore->Fill(avtx->GetY());
  fVertexZBefore->Fill(avtx->GetZ());

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

void AliAnalysisTaskStudentsML::CalculateWeight(Int_t RunNumber, Double_t* weights, Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta)
{
  //GetHistograms

  if(bUsePhiWeights){fWeightsHist[0] = GetHistogramWithWeights(RunNumber, "phi");}
  if(bUsePtWeights){fWeightsHist[1] = GetHistogramWithWeights(RunNumber, "pt");}
  if(bUseEtaWeights){fWeightsHist[2] =GetHistogramWithWeights(RunNumber, "eta");}

  for(Int_t i=0; i<Multi; i++)
  {
    Double_t weight_phi = 1.;
    Double_t weight_pt = 1.;
    Double_t weight_eta = 1.;
    Int_t iBin = 0;

    if(bUsePhiWeights){ iBin = fWeightsHist[0]->FindBin(angles[i]); weight_phi = fWeightsHist[0]->GetBinContent(iBin); }
    if(bUsePtWeights){ iBin = fWeightsHist[1]->FindBin(pt[i]); weight_pt = fWeightsHist[1]->GetBinContent(iBin); }
    if(bUseEtaWeights){ iBin = fWeightsHist[2]->FindBin(eta[i]); weight_eta = fWeightsHist[2]->GetBinContent(iBin); }

    weights[i] = weight_phi*weight_pt*weight_eta;
  }

} //void AliAnalysisTaskStudentsML::CalculateWeight(Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta)

//==========================================================================================================================================================================

TH1F* AliAnalysisTaskStudentsML::GetHistogramWithWeights(Int_t RunNumber, const char *variable)
{
 // Access from external ROOT file the desired histogram with weights.

 // a) Return value;
 // b) Basic protection for arguments;
 // c) Check if the external ROOT file exists at specified path;
 // d) Access the external ROOT file and fetch the desired histogram with weights;
 // e) Close the external ROOT file.

 TH1F* hist = NULL;
 TString sMethod = "void AliAnalysisTaskStudentsML::GetHistogramWithWeights()";

 // b) Basic protection for some arguments:
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  // bail out with your favorite error message
   Fatal(sMethod.Data(), "ERROR: not the correct variable");
 }

 // d) Access the external ROOT file and fetch the desired histogram with weights:
 TFile *weightsFile = NULL;
 if(bGridWeights){weightsFile = TFile::Open(Form("/alice/cern.ch/user/m/mlesch/Weights/%s/%d/Weights.root",fPeriodUsedForWeight.Data(),RunNumber),"READ");}
 if(!bGridWeights) {weightsFile = TFile::Open(Form("%s/Weights.root",fPeriodUsedForWeight.Data()),"READ");}
 if(!weightsFile){ Fatal(sMethod.Data(), "ERROR 404 File not found"); } 
 
 TDirectoryFile *directoryFile = dynamic_cast<TDirectoryFile*>(weightsFile->Get(Form("%s_Weights", variable)));
 if(!directoryFile){Fatal(sMethod.Data(), "Directory not found");} //checking if pointer is set correctly

 TList *List= dynamic_cast<TList*>(directoryFile->Get(Form("%s_Weight=>%.1f-%.1f",variable, fMinCentrality,fMaxCentrality))); 
 if(!List){Fatal(sMethod.Data(), "ERROR: List not found");} 


 // Finally, access the desired histogram:
 hist = dynamic_cast<TH1F*>(List->FindObject(Form("%s_Weight",variable)));
 if(!hist)
 {
  // bail out with your favorite error message
  Fatal(sMethod.Data(), "ERROR: Hist not found");
 } 
 else 
 { 
  hist->SetDirectory(0); // kill the default ownership
 }

 // e) Close the external ROOT file:
 weightsFile->Close(); delete weightsFile;

 return hist;
} // TH1F* AliAnalysisTaskStudentsML::GetHistogramWithWeights( const char *variable, const char *dataset)

//==========================================================================================================================================================================
void AliAnalysisTaskStudentsML::CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights)
{
 // Calculate Q-vectors.

 // a) Make sure all Q-vectors are initially zero;
 // b) Calculate Q-vectors for available angles and weights. 

 // a) Make sure all Q-vectors are initially zero:
 for(Int_t h=0;h<97;h++)
 {
  for(Int_t p=0;p<13;p++)
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
  for(Int_t h=0;h<97;h++)
  {
   for(Int_t p=0;p<13;p++)
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

 }


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

void AliAnalysisTaskStudentsML::Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Double_t* Correlation_Angle, Int_t Correlation_Mult, Double_t* Correlation_Weight)
{
	
       if(h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12!=0.){return;} //protection against anisotropic correlators
	
       // Calculate n-particle correlations from Q-vectors (using recursion):	

	this->CalculateQvectors(Correlation_Mult, Correlation_Angle,Correlation_Weight);
         
        if(2==Number)
        {
         Int_t harmonicsTwoNum[2] = {h1,h2};     
         Int_t harmonicsTwoDen[2] = {0,0};       
         TComplex twoRecursion = Recursion(2,harmonicsTwoNum)/Recursion(2,harmonicsTwoDen).Re();
         Double_t wTwoRecursion = Recursion(2,harmonicsTwoDen).Re();
         fRecursion[0][0]->Fill(0.5,twoRecursion.Re()); // <<cos(h1*phi1+h2*phi2)>>
         fRecursion[0][0]->Fill(1.5,wTwoRecursion); //weight 
         fRecursion[1][0]->Fill(0.5,twoRecursion.Im()); // <<sin(h1*phi1+h2*phi2)>>
         fRecursion[1][0]->Fill(1.5,wTwoRecursion); //weight
	
         }//  2-p correlation
        
        if(3==Number)
        {
         Int_t harmonicsThreeNum[3] = {h1,h2,h3};       
         Int_t harmonicsThreeDen[3] = {0,0,0};       
         TComplex threeRecursion = Recursion(3,harmonicsThreeNum)/Recursion(3,harmonicsThreeDen).Re();
         Double_t wThreeRecursion = Recursion(3,harmonicsThreeDen).Re();
         fRecursion[0][1]->Fill(0.5,threeRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[0][1]->Fill(1.5,wThreeRecursion); //weight
         fRecursion[1][1]->Fill(0.5,threeRecursion.Im()); // <<sin(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[1][1]->Fill(1.5,wThreeRecursion); //weight

         } //  3-p correlation
        
        if(4==Number)
        {
         Int_t harmonicsFourNum[4] = {h1,h2,h3,h4};       
         Int_t harmonicsFourDen[4] = {0,0,0,0};       
         TComplex fourRecursion = Recursion(4,harmonicsFourNum)/Recursion(4,harmonicsFourDen).Re();
         Double_t wFourRecursion = Recursion(4,harmonicsFourDen).Re();
         fRecursion[0][2]->Fill(0.5,fourRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[0][2]->Fill(1.5,wFourRecursion); //weight
         fRecursion[1][2]->Fill(0.5,fourRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[1][2]->Fill(1.5,wFourRecursion); //weight

         }//  4-p correlation
        
        if(5==Number)
        {
         Int_t harmonicsFiveNum[5] = {h1,h2,h3,h4,h5};       
         Int_t harmonicsFiveDen[5] = {0,0,0,0,0};       
         TComplex fiveRecursion = Recursion(5,harmonicsFiveNum)/Recursion(5,harmonicsFiveDen).Re();
         Double_t wFiveRecursion = Recursion(5,harmonicsFiveDen).Re();
         fRecursion[0][3]->Fill(0.5,fiveRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[0][3]->Fill(1.5,wFiveRecursion);
         fRecursion[1][3]->Fill(0.5,fiveRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[1][3]->Fill(1.5,wFiveRecursion);
        }//  5-p correlation

        if(6==Number)
        {
         Int_t harmonicsSixNum[6] = {h1,h2,h3,h4,h5,h6};       
         Int_t harmonicsSixDen[6] = {0,0,0,0,0,0};       
         TComplex sixRecursion = Recursion(6,harmonicsSixNum)/Recursion(6,harmonicsSixDen).Re();
         Double_t wSixRecursion = Recursion(6,harmonicsSixDen).Re();
         fRecursion[0][4]->Fill(0.5,sixRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[0][4]->Fill(1.5,wSixRecursion);
         fRecursion[1][4]->Fill(0.5,sixRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[1][4]->Fill(1.5,wSixRecursion);

         }//  6-p correlation
        
        
        if(7==Number)
        {
         Int_t harmonicsSevenNum[7] = {h1,h2,h3,h4,h5,h6,h7};       
         Int_t harmonicsSevenDen[7] = {0,0,0,0,0,0,0};       
         TComplex sevenRecursion = Recursion(7,harmonicsSevenNum)/Recursion(7,harmonicsSevenDen).Re();
         Double_t wSevenRecursion = Recursion(7,harmonicsSevenDen).Re();
         fRecursion[0][5]->Fill(0.5,sevenRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[0][5]->Fill(1.5,wSevenRecursion);
         fRecursion[1][5]->Fill(0.5,sevenRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[1][5]->Fill(1.5,wSevenRecursion);
  

        }//  7-p correlation
        
        
        if(8==Number)
        {
         Int_t harmonicsEightNum[8] = {h1,h2,h3,h4,h5,h6,h7,h8};       
         Int_t harmonicsEightDen[8] = {0,0,0,0,0,0,0,0};       
         TComplex eightRecursion = Recursion(8,harmonicsEightNum)/Recursion(8,harmonicsEightDen).Re();
         Double_t wEightRecursion = Recursion(8,harmonicsEightDen).Re();
         fRecursion[0][6]->Fill(0.5,eightRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][6]->Fill(1.5,wEightRecursion);
         fRecursion[1][6]->Fill(0.5,eightRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][6]->Fill(1.5,wEightRecursion);
        
        }//  8-p correlation


        if(9==Number)
        {
         Int_t harmonicsNineNum[9] = {h1,h2,h3,h4,h5,h6,h7,h8,h9};       
         Int_t harmonicsNineDen[9] = {0,0,0,0,0,0,0,0,0};       
         TComplex nineRecursion = Recursion(9,harmonicsNineNum)/Recursion(9,harmonicsNineDen).Re();
         Double_t wnineRecursion = Recursion(9,harmonicsNineDen).Re();
         fRecursion[0][7]->Fill(0.5,nineRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][7]->Fill(1.5,wnineRecursion);
         fRecursion[1][7]->Fill(0.5,nineRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][7]->Fill(1.5,wnineRecursion);
        
        }//  8-p correlation

        if(10==Number)
        {
         Int_t harmonicsTenNum[10] = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10};       
         Int_t harmonicsTenDen[10] = {0,0,0,0,0,0,0,0,0,0};       
         TComplex tenRecursion = Recursion(10,harmonicsTenNum)/Recursion(10,harmonicsTenDen).Re();
         Double_t wtenRecursion = Recursion(10,harmonicsTenDen).Re();
         fRecursion[0][8]->Fill(0.5,tenRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][8]->Fill(1.5,wtenRecursion);
         fRecursion[1][8]->Fill(0.5,tenRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][8]->Fill(1.5,wtenRecursion);
        
        }//  10-p correlation

        if(12==Number)
        {
         Int_t harmonicsTwelveNum[12] = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12};       
         Int_t harmonicsTwelveDen[12] = {0,0,0,0,0,0,0,0,0,0,0,0};       
         TComplex twelveRecursion = Recursion(12,harmonicsTwelveNum)/Recursion(12,harmonicsTwelveDen).Re();
         Double_t wtwelveRecursion = Recursion(12,harmonicsTwelveDen).Re();
         fRecursion[0][10]->Fill(0.5,twelveRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][10]->Fill(1.5,wtwelveRecursion);
         fRecursion[1][10]->Fill(0.5,twelveRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][10]->Fill(1.5,wtwelveRecursion);
        
        }//  12-p correlation

        if(Number!=2 && Number!=3 && Number!=4 && Number!=5 && Number!=6 && Number!=7 && Number!=8 && Number!=9 && Number!=10 && Number!=12) { return; }
      
 }//void Correlation() 

//==========================================================================================================================================================================

void AliAnalysisTaskStudentsML::MainTask(Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array)
{

    if(MainTask_Mult>=fMinNumberPart) //do the correlation only if there are more than 8 particles in the event
    { 

    // Calculate Q-vectors for available angles and weights;
     

    Double_t FirstCorrelation=0.;
    Double_t Weight_FirstCorrelation=0.;
    Double_t SecondCorrelation=0.;
    Double_t Weight_SecondCorrelation=0.;

    Double_t FirstCorrelation_Im=0.;
    Double_t Weight_FirstCorrelation_Im=0.;
    Double_t SecondCorrelation_Im=0.;
    Double_t Weight_SecondCorrelation_Im=0.;

    Double_t ThirdCorrelation=0.;
    Double_t Weight_ThirdCorrelation=0.;
    Double_t ThirdCorrelation_Im=0.;
    Double_t Weight_ThirdCorrelation_Im=0.;
    
    //~~~~~~~~~~~~~~~~~

    this->Correlation(fNumber,fh1,fh2,fh3,fh4,fh5,fh6,fh7,fh8,fh9,fh10,fh11,fh12,MainTask_Angle_Array,MainTask_Mult,MainTask_Weight_Array);  
    //do the correlation for the first set

    FirstCorrelation=fRecursion[0][fNumber-2]->GetBinContent(1);
    Weight_FirstCorrelation=fRecursion[0][fNumber-2]->GetBinContent(2);
    FirstCorrelation_Im=fRecursion[1][fNumber-2]->GetBinContent(1);
    Weight_FirstCorrelation_Im=fRecursion[1][fNumber-2]->GetBinContent(2);

    fRecursion[0][fNumber-2]->Reset(); //Reset
    fRecursion[1][fNumber-2]->Reset(); //Reset

    //~~~~~~~~~~~~~~~~~

    this->Correlation(fNumberSecond,fa1,fa2,fa3,fa4,fa5,fa6,fa7,fa8,fa9,fa10,fa11,fa12,MainTask_Angle_Array,MainTask_Mult,MainTask_Weight_Array);  //do the correlation for the second set

    SecondCorrelation=fRecursion[0][fNumberSecond-2]->GetBinContent(1);
    Weight_SecondCorrelation=fRecursion[0][fNumberSecond-2]->GetBinContent(2);
    SecondCorrelation_Im=fRecursion[1][fNumberSecond-2]->GetBinContent(1);
    Weight_SecondCorrelation_Im=fRecursion[1][fNumberSecond-2]->GetBinContent(2);
    
    fRecursion[0][fNumberSecond-2]->Reset(); //Reset
    fRecursion[1][fNumberSecond-2]->Reset(); //Reset


    if(bDoThirdCorrelation)
    {
	this->Correlation(fNumberThird,fb1,fb2,fb3,fb4,fb5,fb6,fb7,fb8,fb9,fb10,fb11,fb12,MainTask_Angle_Array,MainTask_Mult,MainTask_Weight_Array);  
   	 //do the correlation for the first set

   	 ThirdCorrelation=fRecursion[0][fNumberThird-2]->GetBinContent(1);
   	 Weight_ThirdCorrelation=fRecursion[0][fNumberThird-2]->GetBinContent(2);
   	 ThirdCorrelation_Im=fRecursion[1][fNumberThird-2]->GetBinContent(1);
   	 Weight_ThirdCorrelation_Im=fRecursion[1][fNumberThird-2]->GetBinContent(2);
	
  	 fRecursion[0][fNumberThird-2]->Reset(); //Reset
   	 fRecursion[1][fNumberThird-2]->Reset(); //Reset
    }


    //~~~~~~~~~~~~~~~~~

    if(bDoEbERatio){
    	if(TMath::Abs(SecondCorrelation)>=fDenominatorMinValue)
	{
    	   if(bUseRatioWeight){ fEvCentrality->Fill(0.5,(FirstCorrelation)/(SecondCorrelation),Weight_SecondCorrelation); } 
   	   else { fEvCentrality->Fill(0.5,(FirstCorrelation)/(SecondCorrelation),1.); } 
        } //protection against 0, we will come back to this later
    } //if(bDoEbERatio)


    fCentrality->Fill(0.5,FirstCorrelation,Weight_FirstCorrelation); //safe output first set of harmonics
    fCentralitySecond->Fill(0.5,SecondCorrelation,Weight_SecondCorrelation); //safe output second set of harmonics    

    fCentrality->Fill(1.5,FirstCorrelation_Im,Weight_FirstCorrelation_Im); //safe output first set of harmonics
    fCentralitySecond->Fill(1.5,SecondCorrelation_Im,Weight_SecondCorrelation_Im); //safe output second set of harmonics



    if(bDoThirdCorrelation)
    {
	fCentralityThird->Fill(0.5,ThirdCorrelation,Weight_ThirdCorrelation); //safe output second set of harmonics    
   	fCentralityThird->Fill(1.5,ThirdCorrelation_Im,Weight_ThirdCorrelation_Im); //safe output second set of harmonics
    }


  } //if(fParticles>=fMinNumberPart)



} //void AliAnalysisTaskStudentsML::MainTask(Int_t MainTask_Mult, Double_t* MainTask_Angle_Array)

//==========================================================================================================================================================================

  void AliAnalysisTaskStudentsML::MixedParticle(Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)
  {

   if(Mixed_Mult_A>=fMinNumberPart && Mixed_Mult_B>=fMinNumberPart) 
    { 

    Double_t FirstCorrelation=0.;
    Double_t Weight_FirstCorrelation=0.;
    Double_t SecondCorrelation=0.;
    Double_t Weight_SecondCorrelation=0.;

    //Dummy weights, may be implemented later
    Double_t* Dummy_Weights_A = new Double_t[Mixed_Mult_A]; 
    Double_t* Dummy_Weights_B = new Double_t[Mixed_Mult_B]; 
    for(Int_t i=0; i<Mixed_Mult_A; i++){Dummy_Weights_A[i]=1.;}
    for(Int_t i=0; i<Mixed_Mult_B; i++){Dummy_Weights_B[i]=1.;}


   //Calculas for particle group A and B
    this->Correlation(4.,Harmonicus, -Harmonicus, Harmonicus, -Harmonicus,0.,0.,0.,0.,0.,0.,0.,0.,Mixed_Angle_A,Mixed_Mult_A,Dummy_Weights_A);  //do the correlation for the first set

    FirstCorrelation=fRecursion[0][2]->GetBinContent(1);
    Weight_FirstCorrelation=fRecursion[0][2]->GetBinContent(2);

    fRecursion[0][2]->Reset(); //Reset
    fRecursion[1][2]->Reset(); //Reset

    //~~~~~~~~~~~~~~~~~

   this->Correlation(4.,Harmonicus, -Harmonicus, Harmonicus, -Harmonicus,0.,0.,0.,0.,0.,0.,0.,0.,Mixed_Angle_B,Mixed_Mult_B,Dummy_Weights_B); 

    SecondCorrelation=fRecursion[0][2]->GetBinContent(1);
    Weight_SecondCorrelation=fRecursion[0][2]->GetBinContent(2);
    
    fRecursion[0][2]->Reset(); //Reset
    fRecursion[1][2]->Reset(); //Reset

    //~~~~~~~~~~~~~~~~~~
     Double_t Special_Weight = (Double_t)Mixed_Mult_A*((Double_t)Mixed_Mult_A-1.)*(Double_t)Mixed_Mult_B*((Double_t)Mixed_Mult_B-1.);
     TComplex Mixed = TComplex(0.,0.);
     Mixed = CalculateMixedQVectors((Double_t)Harmonicus, Mixed_Mult_A, Mixed_Mult_B, Mixed_Angle_A, Mixed_Angle_B);

    //~~~~~~~~~~~~~~~~~~

   fCentrality->Fill(0.5,(1./Special_Weight)*Mixed.Re(),Special_Weight); //safe output first set of harmonics
   fCentralitySecond->Fill(0.5,FirstCorrelation*SecondCorrelation,Weight_FirstCorrelation*Weight_SecondCorrelation); //safe output second set of harmonics    
   
   fMixedParticleHarmonics->Fill(0.5,FirstCorrelation,Weight_FirstCorrelation);
   fMixedParticleHarmonics->Fill(1.5,SecondCorrelation,Weight_SecondCorrelation);

   delete [] Dummy_Weights_A; 
   delete [] Dummy_Weights_B;  

   }  
 
} //void AliAnalysisTaskStudentsML::MixedParticle(Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)

//==========================================================================================================================================================================




