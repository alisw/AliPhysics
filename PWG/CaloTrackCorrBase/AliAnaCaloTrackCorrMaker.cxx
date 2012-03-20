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

//_____________________________________________________________________________
// Steering class for particle (gamma, hadron) identification and correlation 
// analysis. It is called by the task class AliAnalysisTaskCaloTrackCorrelation 
// and it connects the input (ESD/AOD/MonteCarlo) got with AliCaloTrackReader 
// (produces TClonesArrays of AODs (TParticles in MC case if requested)), with 
// the analysis classes that derive from AliAnaCaloTrackCorrBaseClass
//
// -- Author: Gustavo Conesa (INFN-LNF, LPSC-Grenoble)

#include <cstdlib>

// --- ROOT system ---
#include "TClonesArray.h"
#include "TList.h"
#include "TH1I.h"
//#include <TObjectTable.h>

//---- AliRoot system ---- 
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnaCaloTrackCorrBaseClass.h" 
#include "AliAnaCaloTrackCorrMaker.h" 

ClassImp(AliAnaCaloTrackCorrMaker)


//__________________________________________________
AliAnaCaloTrackCorrMaker::AliAnaCaloTrackCorrMaker() : 
TObject(),
fReader(0),                   fCaloUtils(0),
fOutputContainer(new TList ), fAnalysisContainer(new TList ),
fMakeHisto(kFALSE),           fMakeAOD(kFALSE), 
fAnaDebug(0),                 fCuts(new TList), 
fhNEvents(0x0),               fhTrackMult(0x0),
fhCentrality(0x0)
{
  //Default Ctor
  if(fAnaDebug > 1 ) printf("*** Analysis Maker Constructor *** \n");
  
  //Initialize parameters, pointers and histograms
  InitParameters();
}

//________________________________________________________________________________________
AliAnaCaloTrackCorrMaker::AliAnaCaloTrackCorrMaker(const AliAnaCaloTrackCorrMaker & maker) :   
TObject(),
fReader(),   //(new AliCaloTrackReader(*maker.fReader)),
fCaloUtils(),//(new AliCalorimeterUtils(*maker.fCaloUtils)),
fOutputContainer(new TList()), fAnalysisContainer(new TList()), 
fMakeHisto(maker.fMakeHisto),  fMakeAOD(maker.fMakeAOD),
fAnaDebug(maker.fAnaDebug),    fCuts(new TList()),
fhNEvents(maker.fhNEvents),    fhTrackMult(maker.fhTrackMult),
fhCentrality(maker.fhCentrality)
{
  // cpy ctor
}

//___________________________________________________
AliAnaCaloTrackCorrMaker::~AliAnaCaloTrackCorrMaker() 
{
  // Remove all owned pointers.
  
  //  Do not delete it here, already done somewhere else, need to understand where.	
  //  if (fOutputContainer) {
  //    fOutputContainer->Clear();
  //    delete fOutputContainer ;
  //  }   
  
  if (fAnalysisContainer)
  {
    fAnalysisContainer->Delete();
    delete fAnalysisContainer ;
  }   
  
  if (fReader)    delete fReader ;
  if (fCaloUtils) delete fCaloUtils ;
  
  if(fCuts)
  {
	  fCuts->Delete();
	  delete fCuts;
  }
	
}

//_______________________________________________________
void    AliAnaCaloTrackCorrMaker::AddAnalysis(TObject* ana, Int_t n) 
{
  // Add analysis depending on AliAnaCaloTrackCorrBaseClass to list
  
  if ( fAnalysisContainer)
  { 
    fAnalysisContainer->AddAt(ana,n); 
  }
  else
  { 
    printf("AliAnaCaloTrackCorrMaker::AddAnalysis() - AnalysisContainer not initialized\n");
    abort();
  }
}  

//_________________________________________________________
TList * AliAnaCaloTrackCorrMaker::FillAndGetAODBranchList()
{ 
	
	// Get any new output AOD branches from analysis and put them in a list
	// The list is filled in the maker, and new branch passed to the analysis frame
	// AliAnalysisTaskCaloTrackCorrelation
  
	TList *aodBranchList = fReader->GetAODBranchList() ;
  
	for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
		AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
		if(ana->NewOutputAOD()) aodBranchList->Add(ana->GetCreateOutputAODBranch());
	}
	
	return aodBranchList ;
	
}

//_______________________________________________________
TList * AliAnaCaloTrackCorrMaker::GetListOfAnalysisCuts()
{ 
  
	// Get the list of the cuts used for the analysis
	// The list is filled in the maker, called by the task in LocalInit() and posted there
  
	for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
		AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
		TObjString * objstring = ana->GetAnalysisCuts();
    
		if(objstring)fCuts->Add(objstring);
	}
  
	return fCuts ;
  
}

//___________________________________________________
TList *AliAnaCaloTrackCorrMaker::GetOutputContainer()
{
  // Fill the output list of histograms during the CreateOutputObjects stage.
  
  //Initialize calorimeters  geometry pointers
  //GetCaloUtils()->InitPHOSGeometry();
  //GetCaloUtils()->InitEMCALGeometry();
  
  //General event histograms
  
  fhNEvents        = new TH1I("hNEvents",   "Number of analyzed events"     , 1 , 0 , 1  ) ;
  fhNEvents->SetYTitle("# events");
  fOutputContainer->Add(fhNEvents);
  
  fhTrackMult      = new TH1I("hTrackMult", "Number of tracks per events"   , 2000 , 0 , 2000  ) ;
  fhTrackMult->SetXTitle("# tracks");
  fOutputContainer->Add(fhTrackMult);
  
  fhCentrality     =new TH1F("hCentrality","Number of events in centrality bin",100,0.,100) ;
  fhCentrality->SetXTitle("Centrality bin");
  fOutputContainer->Add(fhCentrality) ;  
  
  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0)
  {
    printf("AliAnaCaloTrackCorrMaker::GetOutputContainer() - Analysis job list not initialized!!!\n");
    return fOutputContainer;
  }
  
  const Int_t buffersize = 255;
  char newname[buffersize];
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    
    if(fMakeHisto) // Analysis with histograms as output on
    {
      
      //Fill container with appropriate histograms			
      TList * templist =  ana ->GetCreateOutputObjects(); 
      templist->SetOwner(kFALSE); //Owner is fOutputContainer.
      
      for(Int_t i = 0; i < templist->GetEntries(); i++)
      {
        
        //Add only  to the histogram name the name of the task
        if(   strcmp((templist->At(i))->ClassName(),"TObjString")   ) 
        {
          snprintf(newname,buffersize, "%s%s", (ana->GetAddedHistogramsStringToName()).Data(), (templist->At(i))->GetName());  
          ((TH1*) templist->At(i))->SetName(newname);
        }
        
        //Add histogram to general container
        fOutputContainer->Add(templist->At(i)) ;
        
      }
      
      delete templist;
      
    }// Analysis with histograms as output on
    
  }//Loop on analysis defined
  
  return fOutputContainer;
  
}

//___________________________________
void AliAnaCaloTrackCorrMaker::Init()
{  
  //Init container histograms and other common variables
  // Fill the output list of histograms during the CreateOutputObjects stage.
  
  //Initialize reader
  GetReader()->Init();
  GetReader()->SetCaloUtils(GetCaloUtils()); // pass the calo utils pointer to the reader
	
  
  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0)
  {
    printf("AliAnaCaloTrackCorrMaker::GetOutputInit() - Analysis job list not initialized!!!\n");
    return;
  }
  
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    
    ana->SetReader(fReader);       // Set Reader for each analysis
    ana->SetCaloUtils(fCaloUtils); // Set CaloUtils for each analysis
    
    ana->Init();
    
  }//Loop on analysis defined
  
}

//_____________________________________________
void AliAnaCaloTrackCorrMaker::InitParameters()
{	
  //Init data members
  
  fMakeHisto  = kTRUE;
  fMakeAOD    = kTRUE; 
  fAnaDebug   = 0; // No debugging info displayed by default
	
}

//______________________________________________________________
void AliAnaCaloTrackCorrMaker::Print(const Option_t * opt) const
{	
  //Print some relevant parameters set for the analysis
	
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  printf("Debug level                =     %d\n", fAnaDebug   ) ;
  printf("Produce Histo              =     %d\n", fMakeHisto  ) ;
  printf("Produce AOD                =     %d\n", fMakeAOD    ) ;
  printf("Number of analysis tasks   =     %d\n", fAnalysisContainer->GetEntries()) ;
  
  if(!strcmp("all",opt))
  {
	  printf("Print analysis Tasks settings :\n") ;
	  for(Int_t iana = 0; iana<fAnalysisContainer->GetEntries(); iana++)
    {
		  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana))->Print("");
	  }
    
	  printf("Print analysis Reader settings :\n") ;
	  fReader->Print("");
	  printf("Print analysis Calorimeter Utils settings :\n") ;
	  fCaloUtils->Print("");
    
  }
  
} 

//_______________________________________________________________________
void AliAnaCaloTrackCorrMaker::ProcessEvent(const Int_t iEntry, 
                                            const char * currentFileName)
{
  //Process analysis for this event
  
  if(fMakeHisto && !fOutputContainer)
  {
    printf("AliAnaCaloTrackCorrMaker::ProcessEvent() - Histograms not initialized\n");
    abort();
  }
	
  if(fAnaDebug >= 0 )
  { 
		printf("***  AliAnaCaloTrackCorrMaker::ProcessEvent() Event %d   ***  \n",iEntry);
	  if(fAnaDebug > 1 ) 
    {
		  printf("AliAnaCaloTrackCorrMaker::ProcessEvent() - Current File Name : %s\n", currentFileName);
		  //printf("fAODBranchList %p, entries %d\n",fAODBranchList,fAODBranchList->GetEntries());
	  }
  }
  
  //Each event needs an empty branch
  TList * aodList = fReader->GetAODBranchList();
  Int_t nAODBranches = aodList->GetEntries();
  for(Int_t iaod = 0; iaod < nAODBranches; iaod++)
  {
	  TClonesArray *tca = dynamic_cast<TClonesArray*> (aodList->At(iaod));
	  if(tca) tca->Clear("C");
  }
  
  //Set geometry matrices before filling arrays, in case recalibration/position calculation etc is needed
  fCaloUtils->AccessGeometry(fReader->GetInputEvent());	
  
  //Set the AODB calibration, bad channels etc. parameters at least once
  fCaloUtils->AccessOADB(fReader->GetInputEvent());	

  
  //Tell the reader to fill the data in the 3 detector lists
  Bool_t ok = fReader->FillInputEvent(iEntry, currentFileName);
  if(!ok)
  {
	  if(fAnaDebug >= 1 )printf("*** Skip event *** %d \n",iEntry);
	  return ;
  }
	
  //Magic line to write events to file
  if(fReader->WriteDeltaAODToFile())AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
  
  //printf(">>>>>>>>>> BEFORE >>>>>>>>>>>\n");
  //gObjectTable->Print();
  
  //Loop on analysis algorithms
  
  if(fAnaDebug > 0 ) printf("*** Begin analysis *** \n");
  
  Int_t nana = fAnalysisContainer->GetEntries() ;
  for(Int_t iana = 0; iana <  nana; iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ; 
    
    ana->ConnectInputOutputAODBranches(); //Sets branches for each analysis
    //Make analysis, create aods in aod branch or AODCaloClusters
    if(fMakeAOD  )  ana->MakeAnalysisFillAOD()  ;
    //Make further analysis with aod branch and fill histograms
    if(fMakeHisto)  ana->MakeAnalysisFillHistograms()  ;
    
  }
	
  // Event control histograms
  fhNEvents   ->Fill(0); //Event analyzed
  fhTrackMult ->Fill(fReader->GetTrackMultiplicity()); 
  fhCentrality->Fill(fReader->GetEventCentrality());
  
  fReader->ResetLists();
  
  //printf(">>>>>>>>>> AFTER >>>>>>>>>>>\n");
  //gObjectTable->Print();
	
  if(fAnaDebug > 0 ) printf("*** End analysis *** \n");
  
}

//__________________________________________________________
void AliAnaCaloTrackCorrMaker::Terminate(TList * outputList)
{  
  //Execute Terminate of analysis
  //Do some final plots.
  
  if (!outputList) 
  {
	  Error("Terminate", "No output list");
	  return;
  }
	
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    if(ana->MakePlotsOn())ana->Terminate(outputList);
    
  }//Loop on analysis defined
  
}

