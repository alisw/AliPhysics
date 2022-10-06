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
#include "AliAnalysisSPC.h"
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
#include <TExMap.h>
#include "TDirectoryFile.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include <TClonesArray.h>

using std::cout;
using std::endl;

ClassImp(AliAnalysisSPC)

//================================================================================================================

AliAnalysisSPC::AliAnalysisSPC(const char *name, Bool_t useParticleWeights): 
 fHistList(NULL),
 fInputList(0),
 fCentrality(0.),
 fDebugLevel(0),
  //SelectionCuts
 bSaveAllQA(kTRUE),
 fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
 fCentralityBins(16),
 bDoFisherYates(kFALSE),
 fFisherYatesCutOff(1.),
 //Weights
 bUseWeightsNUE(kTRUE),
 bUseWeightsNUA(kFALSE),
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
 bComputeEtaGap(kFALSE),
 fEtaGap(0.8),
 bSetSameChargePositive(kTRUE),
 fMixedHarmonic(0),
 fCounterHistogram(NULL),  	
 fProfileTrackCuts(NULL)
 {
  // Constructor.
  printf("AliAnalysisSPC::AliAnalysisSPC(const char *name, Bool_t useParticleWeights)\n");


  // Base list:
  fHistList = new TList();
  fHistList->SetName("outputStudentAnalysis");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

  //DefineOutput(1, TList::Class());  

} // AliAnalysisSPC::AliAnalysisSPC(const char *name, Bool_t useParticleWeights): 

//==========================================================================================================================================================================

AliAnalysisSPC::AliAnalysisSPC():
 fHistList(NULL),
 fInputList(0),
 fCentrality(0.),
 fDebugLevel(0),
 //SelectionCuts
 bSaveAllQA(kTRUE),
 fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
 fCentralityBins(16),
 bDoFisherYates(kFALSE),
 fFisherYatesCutOff(1.),
 //Weights
 bUseWeightsNUE(kTRUE),
 bUseWeightsNUA(kFALSE),
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
 bComputeEtaGap(kFALSE),
 fEtaGap(0.8),
 bSetSameChargePositive(kTRUE),
 fMixedHarmonic(0),
 fCounterHistogram(NULL),  	
 fProfileTrackCuts(NULL)
{
  // Dummy constructor.
  printf("AliAnalysisSPC::AliAnalysisSPC()\n");

  this->InitializeArrays();

} // AliAnalysisSPC::AliAnalysisSPC():

//==========================================================================================================================================================================

AliAnalysisSPC::~AliAnalysisSPC()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisSPC::~AliAnalysisSPC()

//==========================================================================================================================================================================

void AliAnalysisSPC::UserCreateOutputObjects()
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // d) Save cut values
 // a*) Trick to avoid name clashes, part 2.
    
 // a) Trick to avoid name clashes, part 1:
 //Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 //TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 if(bSaveAllQA){this->BookControlHistograms();}
 this->BookFinalResultsHistograms();

 // d) Save cut values 
 if(bDoFisherYates){ fProfileTrackCuts->Fill(0.5, 1); fProfileTrackCuts->Fill(1.5, fFisherYatesCutOff); } 
 fProfileTrackCuts->Fill(2.5, fMinNumberPart);
 if(bUseWeightsNUE){ fProfileTrackCuts->Fill(3.5, 1); } 
 if(bUseWeightsNUA){ fProfileTrackCuts->Fill(4.5, 1); } 

 // *) Trick to avoid name clashes, part 2:
 //TH1::AddDirectory(oldHistAddStatus);

 //PostData(1,fHistList);

} // void AliAnalysisSPC::UserCreateOutputObjects() 

//==========================================================================================================================================================================

void AliAnalysisSPC::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 fCounterHistogram->Fill(0.5); // Checks if User Exec is entered proberly
    
 //Do Analysis
 PhysicsAnalysis();

} // void AliAnalysisSPC::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisSPC::PhysicsAnalysis()
{

 // a) Get Centrality 
 // b.0) Start analysis over AODs;
 // 	b.1) First Loop needed if tracks should be seperated into 2 groups
 // 	b.2) Second loop: Assign track quantities into arrays
 // c) Calculus
 // d) Reset event-by-event objects;
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //a)Select Centrality Bin (has to be valid, as it passed the GlobalQualityAssurance)
 Int_t CentralityBin = SelectCentrality(fCentrality); 
  
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 //Access Data

 //b.0) Start analysis over AODs
 Int_t nTracks = 0;  // number of all tracks in current event 
 nTracks = fInputList->GetEntriesFast(); 

 Double_t* angles_A = NULL; // Azimuthal angles
 Double_t* pt_A = NULL;
 Double_t* eta_A = NULL;
 Double_t* weights_A = NULL;
 Int_t Mult_A = 0.;

 Double_t* angles_B = NULL; 
 Int_t Mult_B = 0.; 
 Int_t CounterSameCharge = 0.; //used if bDifferentCharge = kTRUE



 //b.1) First Loop needed if tracks should be seperated into 2 groups
 if(bDoMixed){ 
	 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) 
	 {

		AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack); // load track

		Double_t charge = aTrack->GetCharge();

		if(bDifferentCharge)
		{
			if(charge>0.){Mult_A+=1.;}
			if(charge<0.){Mult_B+=1.;}
		}//if(bDifferentCharge)

		if(!bDifferentCharge)
		{
			Double_t UsedCharge = 0.;
			if(bSetSameChargePositive) {UsedCharge = charge;}
			if(!bSetSameChargePositive) {UsedCharge = -1.*charge;}

			if(UsedCharge>0.)
			{
			  if(0 == CounterSameCharge%2) {Mult_A+=1.; CounterSameCharge++; }
			  else {Mult_B+=1.; CounterSameCharge++; }
			}

		}//if(!bDifferentCharge)


	 } // First Loop over the tracks if tracks should be seperated into 2 groups
 }
 else { Mult_A = nTracks; } 


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
	if(Mult_A>0){fMultHistogram[CentralityBin][0]->Fill(Mult_A);} //multiplicity distribution after track selection 
   }//if(bSaveAllQA)

   if(Mult_A< fMinNumberPart) { return; } 

 }//if(!bDoMixed)
 else
 {
   if(bSaveAllQA)
   {
  	if(Mult_A>0){fMultHistogram[CentralityBin][0]->Fill(Mult_A);} //multiplicity distribution after track selection
   	if(Mult_B>0){fMultHistogram[CentralityBin][1]->Fill(Mult_B);}
   }//if(bSaveAllQA)

   if(Mult_A < fMinNumberPart || Mult_B < fMinNumberPart ){ return; }

 }//if(bDoMixed)

fCentralityHistogram[CentralityBin]->Fill(fCentrality);

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
//b.2) Second loop: Assign track quantities into arrays
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {  
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack); // load track
    
    Double_t phi = 0.; // azimuthal angle
    Double_t pt = 0.; // Pt (transverse momentum)
    Double_t eta = 0.;
    Double_t charge = 0.; // charge
    Double_t weight = 1.;

    if(!aTrack){continue;} // protection against NULL pointers

    phi = aTrack->Phi(); // azimuthal angle
    pt = aTrack->Pt(); // Pt (transverse momentum)
    eta = aTrack->Eta();
    charge = aTrack->GetCharge(); // charge

    Double_t iEffCorr = 1.;//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
    Double_t iEffInverse = 1.;
    Double_t phi_module_corr = 1.;// doing it in AliJCatalyst while filling track information.

    if (bUseWeightsNUE || bUseWeightsNUA)
    {
      if(bUseWeightsNUE)
      {
  	iEffCorr = aTrack->GetTrackEff();//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
        iEffInverse = 1.0/iEffCorr;
      }
      if(bUseWeightsNUA)
      {
  	phi_module_corr = aTrack->GetWeight();// doing it in AliJCatalyst while filling track information.
      }

      //printf("iEffCorr: %.6f iPhiModuleCorr: %.6f \n", iEffCorr, phi_module_corr);
      
      weight = iEffInverse/phi_module_corr;

    } // End: if (fUseJEfficiency).

     if(!bDoMixed)
     {	
	if(bDoFisherYates)
	{
		Int_t NewIndex = RandomIndexArray[Mult_A];
		if(NewIndex>=After_FisherYates_Mult) continue; //this is outside our Array Size
		
		angles_A[NewIndex] = phi; 
		pt_A[NewIndex] = pt; 
		eta_A[NewIndex] = eta; 
		weights_A[NewIndex] = weight;
	}
	else
	{
		angles_A[Mult_A] = phi; 
		pt_A[Mult_A] = pt; 
		eta_A[Mult_A] = eta; 
		weights_A[Mult_A] = weight;
	} 

        Mult_A+=1.;

     }//if(!bDoMixed)
     else
     {
	if(bDifferentCharge)
	{	
		if(charge>0.){angles_A[Mult_A] = phi; Mult_A+=1.;}
		if(charge<0.){angles_B[Mult_B] = phi; Mult_B+=1.;}
	}//if(bDifferentCharge)
        else
	{
		Double_t UsedCharge = 0.;
		if(bSetSameChargePositive) {UsedCharge = charge;}
		if(!bSetSameChargePositive) {UsedCharge = -1.*charge;}

		if(UsedCharge>0.)
		{
		   if(0 == CounterSameCharge%2) {angles_A[Mult_A] = phi; Mult_A+=1.; CounterSameCharge++; }
		   else {angles_B[Mult_B] = phi; Mult_B+=1.; CounterSameCharge++;}
		}

	}//if(!bDifferentCharge)

     }//if(bDoMixed)


   //............................................................................................
   // Fill control histograms with the particles after track selection:
     if(bSaveAllQA) 
     {
     	if(!bDoMixed)
	{
		fPhiHistogram[CentralityBin][0]->Fill(phi, (1./phi_module_corr)); 
 		fPhiWeightProfile[CentralityBin]->Fill(phi,(1./phi_module_corr));
		fEtaHistogram[CentralityBin][0]->Fill(eta);
		fPTHistogram[CentralityBin][0]->Fill(pt, (1./iEffCorr));
    		fChargeHistogram[CentralityBin]->Fill(charge); 
		
	}//if(!bDoMixed)

	else
	{
		if(charge>0.)
		{
			fPhiHistogram[CentralityBin][0]->Fill(phi); 
			fEtaHistogram[CentralityBin][0]->Fill(eta);
			fPTHistogram[CentralityBin][0]->Fill(pt);
		}//if(charge>0.)
		if(charge<0.)
		{
			fPhiHistogram[CentralityBin][1]->Fill(phi); 
			fEtaHistogram[CentralityBin][1]->Fill(eta);
			fPTHistogram[CentralityBin][1]->Fill(pt);
		}//if(charge<0.)
	}//else 

    } //if(bSaveAllQA)

		


 } //Second loop: Assign track quantities into arrays
 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //c) Calculus
 
 if(!bDoMixed)
 { 
   //add boolean for two particle correlation for eta gaps + calculation alla cindy 
   this->MainTask(CentralityBin, Mult_A, angles_A,weights_A); //Actual Multi-Particle Correlation
 }//if(!bDoMixed)
 
 if(bDoMixed)
 {
   this->MixedParticle(CentralityBin, fMixedHarmonic, Mult_A, angles_A, Mult_B, angles_B);
 }//if(bDoMixed)
 
 if(bComputeEtaGap){ComputeTPCWithEtaGaps(CentralityBin, Mult_A, angles_A, weights_A, eta_A);} 

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // d) Reset event-by-event objects:
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

} //void AliAnalysisSPC::PhysicsAnalysis()

//==========================================================================================================================================================================

void AliAnalysisSPC::Terminate(Option_t *)
{
 // Accessing the merged output list. 
} // end of void AliAnalysisSPC::Terminate(Option_t *)

//==========================================================================================================================================================================

void AliAnalysisSPC::InitializeArrays()
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
  //Centrality Dependend Objects: Lists, QC-Histograms and Output Histograms

  for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
  {
  	//Lists
  	fCentralityList[icent] = NULL;
  	fControlHistogramsList[icent] = NULL;
  	fFinalResultsList[icent] = NULL;

  	//QC-Histograms
    fChargeHistogram[icent] = NULL; 
  	fCentralityHistogram[icent] = NULL;

  	for(Int_t i=0; i<2; i++)
  	{
  	  fPTHistogram[icent][i] = NULL;
  	  fPhiHistogram[icent][i] = NULL;
  	  fEtaHistogram[icent][i] = NULL;
  	  fMultHistogram[icent][i] = NULL;
  	}   

        fPhiWeightProfile[icent] = NULL;

  	//Output Histograms
  	fResults[icent] = NULL;
        fResultsAlternativeError[icent] = NULL; 
  	fCovResults[icent] = NULL; 
  	fJoinedCovResults[icent] = NULL; 
    fMixedParticleHarmonics[icent] = NULL;

        fProfileTPCEta[icent] = NULL;
        
  }

} // void AliAnalysisSPC::InitializeArrays()

//==========================================================================================================================================================================

void AliAnalysisSPC::SetInitializeCentralityArray()
{

 TString sMethodName = "void AliAnalysisSPC::BookAndNestAllLists()";

  Float_t ListCentralities[17] = { fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16 };

  for(Int_t i=0; i<17; i++) { fcentralityArray[i] = ListCentralities[i]; }

  //Protections
  if(fcentralityArray[0] < 0 || fcentralityArray[1] < 0) { Fatal(sMethodName.Data(),"First Centrality bin not defined"); } //The need at least one well defined centrality bin

  for(Int_t icent=0; icent<fCentralityBins; icent++)
  {
	//The next bin should be a valid boundery, i.e. it is > 0. but it is also smaller than the previous boundery -> Wrong ordering
	if( fcentralityArray[icent+1] > 0. && fcentralityArray[icent+1] < fcentralityArray[icent] ) { Fatal(sMethodName.Data(),"Wrong ordering of centrality bounderies"); }
  }

} //void AliAnalysisSPC::InitializeCentralityArray()

//==========================================================================================================================================================================

void AliAnalysisSPC::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisSPC::BookAndNestAllLists()";
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


} // void AliAnalysisSPC::BookAndNestAllLists()

//==========================================================================================================================================================================

void AliAnalysisSPC::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra
 // b) Book histogram to hold phi spectra
 // c) Book histogram to hold eta spectra
 // d) Book histogam to hold multiplicty distributions
 // e) Book histogram for Charge Cut
 // f) Book histogram Centrality

 for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
 {

	//Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	 // a) Book histogram to hold pt spectra:
	 fPTHistogram[icent][0] = new TH1F("fPTHistAfterTrackSelection","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][0]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][0]); 

	 fPTHistogram[icent][1] = new TH1F("fPTHistAfterTrackSelectionSecond","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][1]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][1]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPTHistogram[icent][1]); } 

	 
	 // b) Book histogram to hold phi spectra
	 fPhiHistogram[icent][0] = new TH1F("fPhiHistAfterTrackSelection","Phi Distribution",1000,-TMath::Pi(),TMath::Pi());
	 fPhiHistogram[icent][0]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][0]);

	 fPhiHistogram[icent][1] = new TH1F("fPhiHistAfterTrackSelectionSecond","Phi Distribution",1000,-TMath::Pi(),TMath::Pi()); 
	 fPhiHistogram[icent][1]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][1]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fPhiHistogram[icent][1]); } 

         fPhiWeightProfile[icent] = new TProfile("fPhiWeightProfile","Phi Weights",100,-TMath::Pi(),TMath::Pi()); //centrality dependent output
	 fPhiWeightProfile[icent]->GetXaxis()->SetTitle("#varphi");
	 fPhiWeightProfile[icent]->GetYaxis()->SetTitle("weight");
	 fControlHistogramsList[icent]->Add(fPhiWeightProfile[icent]);

	 // c) Book histogram to hold eta distribution before track selection:
	 fEtaHistogram[icent][0] = new TH1F("fEtaHistAfterTrackSelection","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][0]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][0]);

	 fEtaHistogram[icent][1] = new TH1F("fEtaHistAfterTrackSelectionSecond","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][1]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][1]->SetLineColor(4);
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fEtaHistogram[icent][1]); } 

	 // d) Book histogam to hold multiplicty distributions
	 fMultHistogram[icent][0] = new TH1F("fMultiHistoAfterTrackSelection","Multiplicity",30000,0.,30000.);
	 fMultHistogram[icent][0]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][0]);

	 fMultHistogram[icent][1] = new TH1F("fMultiHistoAfterTrackSelectionSecond","Multiplicity",30000,0.,30000.);
	 fMultHistogram[icent][1]->GetXaxis()->SetTitle("Multiplicity M");
	 if(bDoMixed) { fControlHistogramsList[icent]->Add(fMultHistogram[icent][1]); }

	 // e) Book histogram for Charge
	 fChargeHistogram[icent] = new TH1I("ChargeAfterCut","ChargeAfterCut",11,-5.5,5.5); 
	 fControlHistogramsList[icent]->Add(fChargeHistogram[icent]); 

	 // f) Book histogram Centrality
	 fCentralityHistogram[icent]= new TH1F("fCentralityHistogramAfter","CentralityHistogramAfter",22,0.,110.);
	 fCentralityHistogram[icent]->GetXaxis()->SetTitle("Centrality");
	 fCentralityHistogram[icent]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fCentralityHistogram[icent]);

  }//for(Int_t icent=0; icent<fCentralityBins; icent++)

} //void AliAnalysisSPC::BookControlHistograms()

//==========================================================================================================================================================================

void AliAnalysisSPC::BookFinalResultsHistograms()
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

	  fResultsAlternativeError[icent] = new TProfile("fResultsAlternativeError","Result Analysis First Set Correlators",16,0.,16.); //centrality dependent output
	  fResultsAlternativeError[icent]->GetXaxis()->SetTitle("");
	  fResultsAlternativeError[icent]->GetYaxis()->SetTitle("");
	  fResultsAlternativeError[icent]->Sumw2();
	  fFinalResultsList[icent]->Add(fResultsAlternativeError[icent]);

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


          fProfileTPCEta[icent] = new TProfile("fProfileTPCEta","fProfileTPCEta",9,0.,9.,"s"); 
	  fProfileTPCEta[icent]->GetXaxis()->SetTitle("");
	  fProfileTPCEta[icent]->GetYaxis()->SetTitle("");
	  fProfileTPCEta[icent]->Sumw2();
	  if(bComputeEtaGap){fFinalResultsList[icent]->Add(fProfileTPCEta[icent]);} 
         
          
          
  }

  // Book histogram to debug
  fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks",3,0.,3.);
  fHistList->Add(fCounterHistogram);


 //Profile to save the cut values for track selection
  fProfileTrackCuts = new TProfile("", "", 5, 0., 5.);
  fProfileTrackCuts->SetName("fProfileTrackCuts");
  fProfileTrackCuts->SetTitle("Configuration of the track selection");
  fProfileTrackCuts->SetStats(kFALSE);
  fProfileTrackCuts->GetXaxis()->SetBinLabel(1, "Fisher Yates?"); 
  fProfileTrackCuts->GetXaxis()->SetBinLabel(2, "Keeping Percentage");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(3, "Multiplicity min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(4, "NUE-Weights");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(5, "NUA-Weights");
  fHistList->Add(fProfileTrackCuts);

 
 
} // void AliAnalysisSPC::BookFinalResultsHistograms()

//==========================================================================================================================================================================

Int_t AliAnalysisSPC::SelectCentrality(Double_t CentralityValue)
{

  //Check for centrality bin
  for(Int_t icent=0; icent<fCentralityBins+1; icent++) //loop over all centrality bins
  {
  	if(fcentralityArray[icent]<0){return -1;}
  	if(CentralityValue >= fcentralityArray[icent]) { continue; } 
  	else { return icent-1; } 
  }	

  //We went through all centrality edges without returning. This means, that the measured value is bigger than the maximum centrality that we want for our analyis
  return -1;

} //Int_t AliAnalysisSPC::SelectCentrality(AliAODEvent *aAODevent)

//==========================================================================================================================================================================

void AliAnalysisSPC::FisherYatesRandomizing(Int_t Mult, Int_t *RandomIndex)
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

void AliAnalysisSPC::CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights)
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
  if(bUseWeightsNUE || bUseWeightsNUA){wPhi = CalculateQvectors_weights[i];} //Change some point
  for(Int_t h=0;h<113;h++)
  {
   for(Int_t p=0;p<15;p++)
   {
    if(bUseWeightsNUE || bUseWeightsNUA){wPhiToPowerP = pow(wPhi,p);}
    fQvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi2),wPhiToPowerP*TMath::Sin(h*dPhi2));
   } //  for(Int_t p=0;p<kMaxPower;p++)
  } // for(Int_t h=0;h<kMaxHarmonic;h++)
 } //  for(Int_t i=0;i<fParticles;i++) // loop over particles


} // void AliAnalysisSPC::CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights)

//==========================================================================================================================================================================

TComplex AliAnalysisSPC::CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B){

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

 } //TComplex AliAnalysisSPC::CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B)

//==========================================================================================================================================================================

TComplex AliAnalysisSPC::Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][p];} 
 return TComplex::Conjugate(fQvector[-n][p]);
 
} // TComplex AliAnalysisSPC::Q(Int_t n, Int_t p)


//==========================================================================================================================================================================


TComplex AliAnalysisSPC::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 
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

void AliAnalysisSPC::Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Int_t h13, Int_t h14, Double_t *Correlation_Data)
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

void AliAnalysisSPC::MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array)
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

     fResultsAlternativeError[MainTask_CentBin]->Fill(2.*(Float_t)(i)+0.5,CorrelationNum[i],Weight_CorrelationNum[i]); //safe output first set of harmonics
     fResultsAlternativeError[MainTask_CentBin]->Fill(2.*(Float_t)(i)+1.5,CorrelationDenom[i],Weight_CorrelationDenom[i]); //safe output first set of harmonics

     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+0.5,CorrelationNum[i]*CorrelationDenom[i],Weight_CorrelationNum[i]*Weight_CorrelationDenom[i]); //w_D*N*w_D*D
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+1.5,Weight_CorrelationNum[i]*Weight_CorrelationDenom[i],1.); //w_N*w_D
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+2.5,Weight_CorrelationNum[i],1.); //w_N
     fCovResults[MainTask_CentBin]->Fill(4.*(Float_t)(i)+3.5,Weight_CorrelationDenom[i],1.); //w_D 

     fJoinedCovResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+0.5,CorrelationJoinedCov[i],Weight_CorrelationJoinedCov[i]); //Joined Cov Term z 
     fJoinedCovResults[MainTask_CentBin]->Fill(2.*(Float_t)(i)+1.5,Weight_CorrelationJoinedCov[i],1.); //w_z safe output first set of harmonics
	
   } 

	delete [] Data_Correlation; 

  } //if(fParticles>=fMinNumberPart)

} //void AliAnalysisSPC::MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array)

//==========================================================================================================================================================================

  void AliAnalysisSPC::MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)
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
 
} //void AliAnalysisSPC::MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B)

//==========================================================================================================================================================================

void AliAnalysisSPC::ComputeTPCWithEtaGaps(Int_t CentralityBin, Int_t numberOfParticles, Double_t* angles, Double_t* pWeights, Double_t* pseudorapidity)
{
/* Compute the 2-particle correlators using eta gaps for the current event. */
  TComplex  Qminus[8]   = {TComplex(0., 0.)};   // Q-vectors for the negative subset of the eta range, for v_1 to v_8.
  TComplex  Qplus[8]    = {TComplex(0., 0.)};   // Q-vectors for the positive subset of the eta range, for v_1 to v_8.
  Float_t   Mminus[8]   = {0.};                 // Multiplicity in the negative subset of the eta range.
  Float_t   Mplus[8]    = {0.};                 // Multiplicity in the positive subset of the eta range.
  Float_t   iAngle          = 0.;                     // Azimuthal angle of the current particle.
  Float_t   iWeight         = 1.;                     // Particle weight of the current particle (default: unit weight).
  Float_t   iEta            = 0.;                     // Pseudorapidity of the current particle.
  Float_t   iWeightToP      = 1.;                     // Particle weight rised to the power p.
  TComplex  complexCorrel   = TComplex(0., 0.);       // Complex value of the 2-p correlator.
  Double_t  realCorrel      = 0.;                     // Real value of the 2-p correlator.

  fProfileTPCEta[CentralityBin]->Fill(8.5, fEtaGap,1.); //Fill Eta Gap for saving purpose

// Compute the Q-vectors for the negative and positive subsets of the eta range.
  for (Int_t iPart = 0; iPart < numberOfParticles; iPart++)
  {
  // Read the right elements in the provided arrays.
    iAngle  = angles[iPart];
    iWeight = pWeights[iPart];
    iEta    = pseudorapidity[iPart];
    if (bUseWeightsNUE || bUseWeightsNUA) {iWeightToP = iWeight;}   // All weights are multiplied to get the final one.

  // Compute the Q-vectors.
    if (iEta < 0.)    // Negative subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
          if (iEta < ((-0.5)*fEtaGap))    // Compute only if the particle is in the range.
          {
            Qminus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mminus[iHarmo] += iWeightToP;
          }
          else {continue;}
      }   // End of the loop over the harmonics.
    }   // End of the condition "negative subset".
    else if (iEta > 0.)   // Positive subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
          if (iEta > (0.5*fEtaGap))   // Compute only if the particle is in the range.
          {
            Qplus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mplus[iHarmo] += iWeightToP;
          } 
      }   // End of the loop over the harmonics.
    }   // End of the condition "positive subset".
    else {continue;}    // Particle with iEta = 0.
  }   // End of the loop over the particles for the Q-vectors.

// Compute the 2-p correlators using Qminus and Qplus.
  for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
  {

      if (!( (Qminus[iHarmo].TComplex::Rho() > 0.) && (Qplus[iHarmo].TComplex::Rho() > 0.) )) {continue;}
      if (!( (Mminus[iHarmo] > 0.) && (Mplus[iHarmo] > 0.) )) {continue;}

      complexCorrel = Qminus[iHarmo]*TComplex::Conjugate(Qplus[iHarmo]);
      realCorrel    = (complexCorrel.Re())/(Mminus[iHarmo]*Mplus[iHarmo]);
      fProfileTPCEta[CentralityBin]->Fill(iHarmo + 0.5, realCorrel, Mminus[iHarmo]*Mplus[iHarmo]); //GANESHA declare

    // Reset the 2-particle correlator.
      complexCorrel = TComplex(0.,0.);
      realCorrel    = 0.;
    
  }   // End of the loop over the harmonics.

}   // End of "void ComputeTPCWithEtaGaps()".



