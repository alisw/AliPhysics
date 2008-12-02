/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id: $ */

//_________________________________________________________________________
// Steering class for particle (gamma, hadron) identification and correlation analysis
// It is called by the task class AliAnalysisTaskParticleCorrelation and it connects the input 
// (ESD/AOD/MonteCarlo) got with AliCaloTrackReader (produces TClonesArrays of AODs 
// (TParticles in MC case if requested)), with the 
// analysis classes that derive from AliAnaPartCorrBaseClass
//
// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TClonesArray;
class TString ;
//#include "Riostream.h"

//---- AliRoot system ---- 
#include "AliAnaPartCorrBaseClass.h" 
#include "AliAnaPartCorrMaker.h" 
#include "AliCaloTrackReader.h" 
#include "AliLog.h"


ClassImp(AliAnaPartCorrMaker)


//____________________________________________________________________________
AliAnaPartCorrMaker::AliAnaPartCorrMaker() : 
TObject(),
fOutputContainer(new TList ), fAnalysisContainer(new TList ),
fMakeHisto(0), fMakeAOD(0), fAnaDebug(0), 
fReader(0x0), fAODBranchList(new TList )
{
	//Default Ctor
	if(fAnaDebug > 1 ) printf("*** Analysis Maker  Constructor *** \n");
	
	//Initialize parameters, pointers and histograms
	if(!fReader)
		fReader = new AliCaloTrackReader();
	
	InitParameters();
}

//____________________________________________________________________________
AliAnaPartCorrMaker::AliAnaPartCorrMaker(const AliAnaPartCorrMaker & g) :   
TObject(),
fOutputContainer(g. fOutputContainer), fAnalysisContainer(g.fAnalysisContainer), 
fMakeHisto(g.fMakeHisto), fMakeAOD(fMakeAOD), fAnaDebug(g. fAnaDebug),
fReader(g.fReader), fAODBranchList(g.fAODBranchList)
{
	// cpy ctor
	
}

//_________________________________________________________________________
AliAnaPartCorrMaker & AliAnaPartCorrMaker::operator = (const AliAnaPartCorrMaker & source)
{
	// assignment operator
	
	if(this == &source)return *this;
	((TObject *)this)->operator=(source);
	
	fOutputContainer    = source.fOutputContainer ;
	fAnalysisContainer  = source.fAnalysisContainer ;
	fAnaDebug           = source.fAnaDebug;
	fMakeHisto          = source.fMakeHisto;
	fMakeAOD            = source.fMakeAOD;
	
	fReader             = source.fReader ;
	fAODBranchList      = source.fAODBranchList;
	
	return *this;
	
}

//____________________________________________________________________________
AliAnaPartCorrMaker::~AliAnaPartCorrMaker() 
{
	// Remove all pointers.
	
	// Protection added in case of NULL pointers (MG)
	if (fOutputContainer) {
		fOutputContainer->Clear();
		delete fOutputContainer ;
	}   
	
	if (fAnalysisContainer) {
		fAnalysisContainer->Clear();
		delete fAnalysisContainer ;
	}   
	
	if (fReader) delete fReader ;
	
	
	if(fAODBranchList){
//		for(Int_t iaod = 0; iaod < fAODBranchList->GetEntries(); iaod++)
//			fAODBranchList->At(iaod)->Clear();
	
		fAODBranchList->Clear();
		delete fAODBranchList ;
	}
	
}

//________________________________________________________________________
TList * AliAnaPartCorrMaker::GetAODBranchList()
{ 

// Get any new output AOD branches from analysis and put them in a list
// The list is filled in the maker, and new branch passed to the analysis frame
// AliAnalysisTaskPartCorr
 
	for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
			
		AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
		if(ana->NewOutputAOD()) fAODBranchList->Add(ana->GetCreateOutputAODBranch());
	}

	return fAODBranchList ;

}

//________________________________________________________________________
void AliAnaPartCorrMaker::Init()
{  
	//Init container histograms and other common variables
	
	if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0)
			AliFatal("Analysis job list not initialized");
			
	for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
			
		AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
		ana->SetReader(fReader); //SetReader for each analysis
		ana->Init();
		
		if(fMakeHisto){// Analysis with histograms as output on
			//Fill container with appropriate histograms			
			TList * templist =  ana -> GetCreateOutputObjects(); 
				for(Int_t i = 0; i < templist->GetEntries(); i++)
					fOutputContainer->Add(templist->At(i)) ;
		}// Analysis with histograms as output on
	}//Loop on analysis defined
}

//____________________________________________________________________________
void AliAnaPartCorrMaker::InitParameters()
{
	
	//Init data members
	fMakeHisto = kTRUE;
	fMakeAOD = kTRUE; 
	fAnaDebug = 0; // No debugging info displayed by default

}

//__________________________________________________________________
void AliAnaPartCorrMaker::Print(const Option_t * opt) const
{
	
	//Print some relevant parameters set for the analysis
	if(! opt)
		return;
	
	printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
	printf("Debug level                =     %d\n", fAnaDebug) ;
	printf("Produce Histo              =     %d\n", fMakeHisto) ;
	printf("Produce AOD                =     %d\n", fMakeAOD) ;
	
} 


//____________________________________________________________________________
Bool_t AliAnaPartCorrMaker::ProcessEvent(Int_t iEntry){
	//Process analysis for this event
	
	if(fMakeHisto && !fOutputContainer)
		AliFatal("Histograms not initialized");
	
	if(fAnaDebug >= 0 ) printf("***  Event %d   ***  \n",iEntry);
	
	//Each event needs an empty branch	
	for(Int_t iaod = 0; iaod < fAODBranchList->GetEntries(); iaod++)
		fAODBranchList->At(iaod)->Clear();
	
	//Tell the reader to fill the data in the 3 detector lists
	fReader->FillInputEvent();
	
	//Loop on analysis algorithms
	if(fAnaDebug > 0 ) printf("*** Begin analysis *** \n");
	Int_t nana = fAnalysisContainer->GetEntries() ;
	for(Int_t iana = 0; iana <  nana; iana++){
		
		AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ; 
		
		ana->ConnectInputOutputAODBranches(); //Sets branches for each analysis
		
		//Make analysis, create aods in aod branch or AODCaloClusters
		if(fMakeAOD) ana->MakeAnalysisFillAOD()  ;
		//Make further analysis with aod branch and fill histograms
		if(fMakeHisto) ana->MakeAnalysisFillHistograms()  ;
		
	}
		
	fReader->ResetLists();
	
	if(fAnaDebug > 0 ) printf("*** End analysis *** \n");
	
	return kTRUE ;
	
}
