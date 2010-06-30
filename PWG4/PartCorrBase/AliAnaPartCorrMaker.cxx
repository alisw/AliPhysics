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

#include <cstdlib>

// --- ROOT system ---
#include "TClonesArray.h"
#include "TList.h"
#include "TH1I.h"
//#include "Riostream.h"
//#include <TObjectTable.h>

//---- AliRoot system ---- 
#include "AliAnaPartCorrBaseClass.h" 
#include "AliAnaPartCorrMaker.h" 
#include "AliCaloTrackReader.h" 


ClassImp(AliAnaPartCorrMaker)


//____________________________________________________________________________
AliAnaPartCorrMaker::AliAnaPartCorrMaker() : 
TObject(),
fOutputContainer(new TList ), fAnalysisContainer(new TList ),
fMakeHisto(kFALSE), fMakeAOD(kFALSE), fMakeMixing(kFALSE), fAnaDebug(0), 
fReader(new AliCaloTrackReader()), fCaloUtils(new AliCalorimeterUtils()), 
fAODBranchList(new TList ),fCuts(new TList), fhNEvents(0x0)
{
  //Default Ctor
  if(fAnaDebug > 1 ) printf("*** Analysis Maker  Constructor *** \n");
  
  //Initialize parameters, pointers and histograms
  InitParameters();
}

//____________________________________________________________________________
AliAnaPartCorrMaker::AliAnaPartCorrMaker(const AliAnaPartCorrMaker & maker) :   
TObject(),
fOutputContainer(new TList()), fAnalysisContainer(new TList()), 
fMakeHisto(maker.fMakeHisto), fMakeAOD(maker.fMakeAOD), fMakeMixing(maker.fMakeMixing), 
fAnaDebug(maker.fAnaDebug),
fReader(),//new AliCaloTrackReader(*maker.fReader)), 
fCaloUtils(),//(new AliCalorimeterUtils(*maker.fCaloUtils)),
fAODBranchList(new TList()), fCuts(new TList()), 
fhNEvents(maker.fhNEvents)
{
  // cpy ctor
	
}

//_________________________________________________________________________
//AliAnaPartCorrMaker & AliAnaPartCorrMaker::operator = (const AliAnaPartCorrMaker & source)
//{
//  // assignment operator
//  
//  if(this == &source)return *this;
//  ((TObject *)this)->operator=(source);
//  
//  delete fOutputContainer;   fOutputContainer    = source.fOutputContainer ;
//  delete fAnalysisContainer; fAnalysisContainer  = source.fAnalysisContainer ;
//  fAnaDebug           = source.fAnaDebug;
//  fMakeHisto          = source.fMakeHisto;
//  fMakeAOD            = source.fMakeAOD;
//  
//  delete fReader ;        fReader             = new AliCaloTrackReader(*source.fReader) ;
//  delete fAODBranchList; fAODBranchList      = source.fAODBranchList;
//	
//  return *this;
//  
//}

//____________________________________________________________________________
AliAnaPartCorrMaker::~AliAnaPartCorrMaker() 
{
  // Remove all owned pointers.
//printf("======Delete Maker \n");

//  Do not delete it here, already done somewhere else, need to understand where.	
//  if (fOutputContainer) {
//    fOutputContainer->Clear();
//    delete fOutputContainer ;
//  }   
  
  if (fAnalysisContainer) {
    fAnalysisContainer->Delete();
    delete fAnalysisContainer ;
  }   
  
  if (fReader)    delete fReader ;
  if (fCaloUtils) delete fCaloUtils ;

	
  if(fAODBranchList){
//		for(Int_t iaod = 0; iaod < fAODBranchList->GetEntries(); iaod++)
//			fAODBranchList->At(iaod)->Clear();
	
    fAODBranchList->Delete();
    delete fAODBranchList ;
  }
  
  if(fCuts){
	  fCuts->Delete();
	  delete fCuts;
  }
	
//	printf("====== Maker deleted \n");
}

//________________________________________________________________________
TList * AliAnaPartCorrMaker::GetListOfAnalysisCuts()
{ 

	// Get the list of the cuts used for the analysis
	// The list is filled in the maker, called by the task in LocalInit() and posted there
	//printf("GetListOfAnalysis! \n");

	for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
		AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
		TObjString * objstring = ana->GetAnalysisCuts();
		//printf("analysis %d cuts %p\n",iana,objstring);
		if(objstring)fCuts->Add(objstring);
	}
    //printf("Maker::GetEntries() %d\n",fCuts->GetEntries());
	return fCuts ;
  
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
TList *AliAnaPartCorrMaker::GetOutputContainer()
{
// Fill the output list of histograms during the CreateOutputObjects stage.
  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0){
    printf("AliAnaPartCorrMaker::GetOutputContainer() - Analysis job list not initialized!!!\n");
    //abort();
  }

  //Initialize the geometry pointers
  GetCaloUtils()->InitPHOSGeometry();
  GetCaloUtils()->InitEMCALGeometry();
		
  char newname[128];
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
    AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    if(fMakeHisto){// Analysis with histograms as output on
      //Fill container with appropriate histograms			
      TList * templist =  ana ->GetCreateOutputObjects(); 
	  templist->SetOwner(kFALSE); //Owner is fOutputContainer.
      for(Int_t i = 0; i < templist->GetEntries(); i++){

	//Add only  to the histogram name the name of the task
	if(   strcmp((templist->At(i))->ClassName(),"TObjString")   ) {
	  sprintf(newname,"%s%s", (ana->GetAddedHistogramsStringToName()).Data(), (templist->At(i))->GetName());  
	  ((TH1*) templist->At(i))->SetName(newname);
	}
	//Add histogram to general container
	fOutputContainer->Add(templist->At(i)) ;
      }
		delete templist;
    }// Analysis with histograms as output on
  }//Loop on analysis defined
  
	
  fhNEvents        = new TH1I("hNEvents", "Number of analyzed events"   , 1 , 0 , 1  ) ;
  fOutputContainer->Add(fhNEvents);
	
  return fOutputContainer;

}

//________________________________________________________________________
void AliAnaPartCorrMaker::Init()
{  
  //Init container histograms and other common variables
  // Fill the output list of histograms during the CreateOutputObjects stage.
 
  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0){
    printf("AliAnaPartCorrMaker::GetOutputInit() - Analysis job list not initialized!!!\n");
    //abort();
  }

  //Initialize the geometry pointers
  //GetCaloUtils()->InitPHOSGeometry();
  //GetCaloUtils()->InitEMCALGeometry();
	
  //Initialize reader
  fReader->Init();
  fReader->SetCaloUtils(fCaloUtils); // pass the calo utils pointer to the reader
	
  //fCaloUtils->Init();
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
    
    AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    ana->SetReader(fReader); //SetReader for each analysis
	ana->SetCaloUtils(fCaloUtils); //Set CaloUtils for each analysis

    ana->Init();
    
  }//Loop on analysis defined
}

//____________________________________________________________________________
void AliAnaPartCorrMaker::InitParameters()
{	
  //Init data members
  
  fMakeHisto  = kTRUE;
  fMakeAOD    = kTRUE; 
  fMakeMixing = kFALSE;
  fAnaDebug   = 0; // No debugging info displayed by default
	
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
  printf("Mixing Analysis            =     %d\n", fMakeMixing) ;
  printf("Number of analysis tasks   =     %d\n", fAnalysisContainer->GetEntries()) ;
  if(!strcmp("all",opt)){
	  printf("Print analysis Tasks settings :\n") ;
	  for(Int_t iana = 0; iana<fAnalysisContainer->GetEntries(); iana++){
		  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana))->Print("");
	  }
  
	  printf("Print analysis Reader settings :\n") ;
	  fReader->Print("");
	  printf("Print analysis Calorimeter Utils settings :\n") ;
	  fCaloUtils->Print("");
  }
} 


//____________________________________________________________________________
void AliAnaPartCorrMaker::ProcessEvent(const Int_t iEntry, const char * currentFileName){
  //Process analysis for this event
  
  if(fMakeHisto && !fOutputContainer){
    printf("AliAnaPartCorrMaker::ProcessEvent() - Histograms not initialized\n");
    abort();
  }
	
  if(fAnaDebug >= 0 ){ 
		printf("***  Event %d   ***  \n",iEntry);
	  if(fAnaDebug > 1 ) {
		  printf("AliAnaPartCorrMaker::ProcessEvent() - Current File Name : %s\n", currentFileName);
		  //printf("fAODBranchList %p, entries %d\n",fAODBranchList,fAODBranchList->GetEntries());
	  }
  }
  //Each event needs an empty branch
  Int_t nAODBranches = fAODBranchList->GetEntries();
  for(Int_t iaod = 0; iaod < nAODBranches; iaod++){
	  //fAODBranchList->At(iaod)->Clear();
	  TClonesArray *tca = dynamic_cast<TClonesArray*> (fAODBranchList->At(iaod));
	  if(tca) tca->Delete();
  }

  //Tell the reader to fill the data in the 3 detector lists
  Bool_t ok = fReader->FillInputEvent(iEntry, currentFileName);
  if(!ok){
	  if(fAnaDebug >= 1 )printf("*** Skip event *** %d \n",iEntry);
	  return ;
  }
	
  fCaloUtils->SetGeometryTransformationMatrices(fReader->GetInputEvent());	
	
  //printf(">>>>>>>>>> BEFORE >>>>>>>>>>>\n");
  //gObjectTable->Print();
  //Loop on analysis algorithms
  if(fAnaDebug > 0 ) printf("*** Begin analysis *** \n");
  Int_t nana = fAnalysisContainer->GetEntries() ;
  for(Int_t iana = 0; iana <  nana; iana++){
    AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ; 
    
	ana->ConnectInputOutputAODBranches(); //Sets branches for each analysis
    //Make analysis, create aods in aod branch or AODCaloClusters
	if(fMakeAOD   && !fMakeMixing)  ana->MakeAnalysisFillAOD()  ;
    //Make further analysis with aod branch and fill histograms
    if(fMakeHisto && !fMakeMixing)  ana->MakeAnalysisFillHistograms()  ;
    //Make analysis with delta AODs of different events
	if(fMakeMixing)                 ana->MakeMixingAnalysisFillHistograms()  ;

  }
	
  fhNEvents->Fill(0); //Event analyzed
	
  fReader->ResetLists();
  
  //printf(">>>>>>>>>> AFTER >>>>>>>>>>>\n");
  //gObjectTable->Print();
	
  if(fAnaDebug > 0 ) printf("*** End analysis *** \n");
  
}

//________________________________________________________________________
void AliAnaPartCorrMaker::Terminate(TList * outputList)
{  
  //Execute Terminate of analysis
  //Do some final plots.
  
  if (!outputList) {
	  Error("Terminate", "No output list");
	  return;
  }
	
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++){
    
    AliAnaPartCorrBaseClass * ana =  ((AliAnaPartCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    ana->Terminate(outputList);
    
  }//Loop on analysis defined
}
