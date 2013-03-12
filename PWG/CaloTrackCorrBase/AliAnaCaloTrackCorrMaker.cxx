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
#include "TH1F.h"
//#include <TObjectTable.h>

//---- AliRoot system ----
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
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
fScaleFactor(-1),
fhNEvents(0),                 fhNPileUpEvents(0),
fhZVertex(0),                 
fhPileUpClusterMult(0),       fhPileUpClusterMultAndSPDPileUp(0),
fhTrackMult(0),
fhCentrality(0),              fhEventPlaneAngle(0),
fhNMergedFiles(0),            fhScaleFactor(0),
fhEMCalBCEvent(0),            fhEMCalBCEventCut(0),
fhTrackBCEvent(0),            fhTrackBCEventCut(0),
fhPrimaryVertexBC(0),         fhTimeStampFraction(0),
fhNPileUpVertSPD(0),          fhNPileUpVertTracks(0)
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
fScaleFactor(maker.fScaleFactor),
fhNEvents(maker.fhNEvents), 
fhNPileUpEvents(maker.fhNPileUpEvents),
fhZVertex(maker.fhZVertex),    
fhPileUpClusterMult(maker.fhPileUpClusterMult),
fhPileUpClusterMultAndSPDPileUp(maker.fhPileUpClusterMultAndSPDPileUp),
fhTrackMult(maker.fhTrackMult),
fhCentrality(maker.fhCentrality),
fhEventPlaneAngle(maker.fhEventPlaneAngle),
fhNMergedFiles(maker.fhNMergedFiles),          
fhScaleFactor(maker.fhScaleFactor),
fhEMCalBCEvent(maker.fhEMCalBCEvent),
fhEMCalBCEventCut(maker.fhEMCalBCEventCut),
fhTrackBCEvent(maker.fhTrackBCEvent),
fhTrackBCEventCut(maker.fhTrackBCEventCut),
fhPrimaryVertexBC(maker.fhPrimaryVertexBC),
fhTimeStampFraction(maker.fhTimeStampFraction),
fhNPileUpVertSPD(maker.fhNPileUpVertSPD),
fhNPileUpVertTracks(maker.fhNPileUpVertTracks)
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

//__________________________________________________________________
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

//_________________________________________________________
void AliAnaCaloTrackCorrMaker::FillControlHistograms()
{
  // Event control histograms
  
  AliVEvent* event =  fReader->GetInputEvent();
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);
  
  fhNEvents        ->Fill(0); // Number of events analyzed
  
  if( fReader->IsPileUpFromSPD())
    fhNPileUpEvents->Fill(0.5);
  //if( event->IsPileupFromSPDInMultBins())
  //  fhNPileUpEvents->Fill(1.5);
  if( fReader->IsPileUpFromEMCal())
    fhNPileUpEvents->Fill(2.5);
  if( fReader->IsPileUpFromSPDOrEMCal() )
    fhNPileUpEvents->Fill(3.5);
  if( fReader->IsPileUpFromSPDAndEMCal() )
    fhNPileUpEvents->Fill(4.5);
  if( fReader->IsPileUpFromSPDAndNotEMCal() )
    fhNPileUpEvents->Fill(5.5);
  if( fReader->IsPileUpFromEMCalAndNotSPD() )
    fhNPileUpEvents->Fill(6.5);
  if( fReader->IsPileUpFromNotSPDAndNotEMCal() )
    fhNPileUpEvents->Fill(7.5);
  
  if(fReader->IsPileUpFromSPD())
    fhPileUpClusterMultAndSPDPileUp ->Fill(fReader->GetNPileUpClusters());
    
  fhPileUpClusterMult ->Fill(fReader->GetNPileUpClusters  ());
  fhTrackMult         ->Fill(fReader->GetTrackMultiplicity());
  fhCentrality        ->Fill(fReader->GetEventCentrality  ());
  fhEventPlaneAngle   ->Fill(fReader->GetEventPlaneAngle  ());
  
  for(Int_t i = 0; i < 19; i++)
  {
    if(fReader->GetTrackEventBC(i))   fhTrackBCEvent   ->Fill(i);
    if(fReader->GetTrackEventBCcut(i))fhTrackBCEventCut->Fill(i);
    if(fReader->GetEMCalEventBC(i))   fhEMCalBCEvent   ->Fill(i);
    if(fReader->GetEMCalEventBCcut(i))fhEMCalBCEventCut->Fill(i);
  }
  
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v) ;
  fhZVertex->Fill(v[2]);
  
  Int_t bc = fReader->GetVertexBC();
  if(bc!=AliVTrack::kTOFBCNA)fhPrimaryVertexBC->Fill(bc);
  
  
  // N pile up vertices
  Int_t nVerticesSPD    = -1;
  Int_t nVerticesTracks = -1;
  
  if      (esdevent)
  {
    nVerticesSPD    = esdevent->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = esdevent->GetNumberOfPileupVerticesTracks();
    
  }//ESD
  else if (aodevent)
  {
    nVerticesSPD    = aodevent->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = aodevent->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  fhNPileUpVertSPD   ->Fill(nVerticesSPD);
  fhNPileUpVertTracks->Fill(nVerticesTracks);

  // Time stamp
  if(fReader->IsSelectEventTimeStampOn() && esdevent)
  {
    Int_t timeStamp = esdevent->GetTimeStamp();
    Float_t timeStampFrac = 1.*(timeStamp-fReader->GetRunTimeStampMin()) /
                               (fReader->GetRunTimeStampMax()-fReader->GetRunTimeStampMin());
    
    //printf("stamp %d, min %d, max %d, frac %f\n", timeStamp, fReader->GetRunTimeStampMin(), fReader->GetRunTimeStampMax(), timeStampFrac);

    fhTimeStampFraction->Fill(timeStampFrac);
  }
  
  
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
  
  fhNEvents      = new TH1F("hNEvents",   "Number of analyzed events"     , 1 , 0 , 1  ) ;
  fhNEvents->SetYTitle("# events");
  fOutputContainer->Add(fhNEvents);
  
  fhNPileUpEvents      = new TH1F("hNPileUpEvents",   "Number of events considered as pile-up", 8 , 0 , 8 ) ;
  fhNPileUpEvents->SetYTitle("# events");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(1 ,"SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(2 ,"Multi SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(3 ,"EMCal");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(4 ,"EMCal || SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(5 ,"EMCal && SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(6 ,"!EMCal && SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(7 ,"EMCal && !SPD");
  fhNPileUpEvents->GetXaxis()->SetBinLabel(8 ,"!EMCal && !SPD");
  fOutputContainer->Add(fhNPileUpEvents);
  
  fhTrackBCEvent      = new TH1F("hTrackBCEvent",   "Number of events with at least 1 track in a bunch crossing ", 19 , 0 , 19 ) ;
  fhTrackBCEvent->SetYTitle("# events");
  fhTrackBCEvent->SetXTitle("Bunch crossing");
  for(Int_t i = 1; i < 20; i++)
    fhTrackBCEvent->GetXaxis()->SetBinLabel(i ,Form("%d",i-10));
  fOutputContainer->Add(fhTrackBCEvent);
  
  fhTrackBCEventCut      = new TH1F("hTrackBCEventCut",   "Number of events with at least 1 track in a bunch crossing ", 19 , 0 , 19 ) ;
  fhTrackBCEventCut->SetYTitle("# events");
  fhTrackBCEventCut->SetXTitle("Bunch crossing");
  for(Int_t i = 1; i < 20; i++)
    fhTrackBCEventCut->GetXaxis()->SetBinLabel(i ,Form("%d",i-10));
  fOutputContainer->Add(fhTrackBCEventCut);

  fhPrimaryVertexBC      = new TH1F("hPrimaryVertexBC", "Number of primary vertex per bunch crossing ", 41 , -20 , 20  ) ;
  fhPrimaryVertexBC->SetYTitle("# events");
  fhPrimaryVertexBC->SetXTitle("Bunch crossing");
  fOutputContainer->Add(fhPrimaryVertexBC);

  fhEMCalBCEvent      = new TH1F("hEMCalBCEvent",   "Number of events with at least 1 cluster in a bunch crossing ", 19 , 0 , 19 ) ;
  fhEMCalBCEvent->SetYTitle("# events");
  fhEMCalBCEvent->SetXTitle("Bunch crossing");
  for(Int_t i = 1; i < 20; i++)
    fhEMCalBCEvent->GetXaxis()->SetBinLabel(i ,Form("%d",i-10));
  fOutputContainer->Add(fhEMCalBCEvent);
  
  fhEMCalBCEventCut      = new TH1F("hEMCalBCEventCut",   "Number of events with at least 1 cluster in a bunch crossing", 19 , 0 , 19 ) ;
  fhEMCalBCEventCut->SetYTitle("# events");
  fhEMCalBCEventCut->SetXTitle("Bunch crossing");
  for(Int_t i = 1; i < 20; i++)
    fhEMCalBCEventCut->GetXaxis()->SetBinLabel(i ,Form("%d",i-10));
  fOutputContainer->Add(fhEMCalBCEventCut);
  
  fhZVertex      = new TH1F("hZVertex", " Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZVertex);
  
  fhTrackMult    = new TH1F("hTrackMult", "Number of tracks per events"   , 2000 , 0 , 2000  ) ;
  fhTrackMult->SetXTitle("# tracks");
  fOutputContainer->Add(fhTrackMult);
  
  fhPileUpClusterMult    = new TH1F("hPileUpClusterMult", "Number of clusters per event with large time (|t| > 20 ns)" , 100 , 0 , 100  ) ;
  fhPileUpClusterMult->SetXTitle("# clusters");
  fOutputContainer->Add(fhPileUpClusterMult);
  
  fhPileUpClusterMultAndSPDPileUp = new TH1F("hPileUpClusterMultAndSPDPileUp", "Number of clusters per event with large time (|t| > 20 ns, events tagged as pile-up by SPD)" , 100 , 0 , 100 ) ;
  fhPileUpClusterMultAndSPDPileUp->SetXTitle("# clusters");
  fOutputContainer->Add(fhPileUpClusterMultAndSPDPileUp);
  
  fhNPileUpVertSPD  = new TH1F ("hNPileUpVertSPD","N pile-up SPD vertex", 50,0,50);
  fhNPileUpVertSPD->SetYTitle("# vertex ");
  fOutputContainer->Add(fhNPileUpVertSPD);
  
  fhNPileUpVertTracks  = new TH1F ("hNPileUpVertTracks","N pile-up Tracks vertex", 50,0,50);
  fhNPileUpVertTracks->SetYTitle("# vertex ");
  fOutputContainer->Add(fhNPileUpVertTracks);
  
  fhCentrality   = new TH1F("hCentrality","Number of events in centrality bin",100,0.,100) ;
  fhCentrality->SetXTitle("Centrality bin");
  fOutputContainer->Add(fhCentrality) ;  
  
  fhEventPlaneAngle=new TH1F("hEventPlaneAngle","Number of events in event plane",100,0.,TMath::Pi()) ;
  fhEventPlaneAngle->SetXTitle("EP angle (rad)");
  fOutputContainer->Add(fhEventPlaneAngle) ;
  
  if(fReader->IsSelectEventTimeStampOn())
  {
    fhTimeStampFraction = new TH1F("hTimeStampFraction","Fraction of events within a given time stamp range",150, -1, 2) ;
    fhTimeStampFraction->SetXTitle("fraction");
    fOutputContainer->Add(fhTimeStampFraction) ;
  }
  
  if(fScaleFactor > 0)
  {
    fhNMergedFiles = new TH1F("hNMergedFiles",   "Number of merged output files"     , 1 , 0 , 1  ) ;
    fhNMergedFiles->SetYTitle("# files");
    fhNMergedFiles->Fill(1); // Fill here with one entry, while merging it will count the rest
    fOutputContainer->Add(fhNMergedFiles);
    
    fhScaleFactor = new TH1F("hScaleFactor",   "Number of merged output files"     , 1 , 0 , 1  ) ;
    fhScaleFactor->SetYTitle("scale factor");
    fhScaleFactor->SetBinContent(1,fScaleFactor); // Fill here 
    fOutputContainer->Add(fhScaleFactor);    
  }
  
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
          //printf("name %s, new name %s\n",(templist->At(i))->GetName(),newname);
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
    fReader->ResetLists();
	  return ;
  }
	
  //Magic line to write events to file
  if(fReader->WriteDeltaAODToFile())AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
  
  //printf(">>>>>>>>>> BEFORE >>>>>>>>>>>\n");
  //gObjectTable->Print();
  
  //Access pointers, and trigger mask check needed in mixing case
  AliAnalysisManager   *manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  UInt_t isMBTrigger = kFALSE;
  UInt_t isTrigger   = kFALSE;
  if(inputHandler)
  {
    isMBTrigger = inputHandler->IsEventSelected() & fReader->GetMixEventTriggerMask();
    isTrigger   = inputHandler->IsEventSelected() & fReader->GetEventTriggerMask();
  }
  
  //Loop on analysis algorithms
  
  if(fAnaDebug > 0 ) printf("*** Begin analysis *** \n");
  
  Int_t nana = fAnalysisContainer->GetEntries() ;
  for(Int_t iana = 0; iana <  nana; iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ; 
    
    ana->ConnectInputOutputAODBranches(); //Sets branches for each analysis
    
    //Fill pool for mixed event for the analysis that need it
    if(!fReader->IsEventTriggerAtSEOn() && isMBTrigger)
    {
      ana->FillEventMixPool();
      continue; // pool filled do not try to fill AODs or histograms
    }
    
    //Make analysis, create aods in aod branch and in some cases fill histograms
    if(fMakeAOD  )  ana->MakeAnalysisFillAOD()  ;
    
    //Make further analysis with aod branch and fill histograms
    if(fMakeHisto)  ana->MakeAnalysisFillHistograms()  ;
    
  }
	
  fReader->ResetLists();

  // In case of mixing analysis, non triggered events are used,
  // do not fill control histograms for a non requested triggered event
  if(!fReader->IsEventTriggerAtSEOn() && !isTrigger)
  {    
    if(fAnaDebug > 0 ) printf("AliAnaCaloTrackMaker::ProcessEvent() - *** End analysis, MB for mixing *** \n");
    return;
  }
  
  FillControlHistograms();
  
  //printf(">>>>>>>>>> AFTER >>>>>>>>>>>\n");
  //gObjectTable->Print();
	
  if(fAnaDebug > 0 ) printf("AliAnaCaloTrackMaker::ProcessEvent() - *** End analysis *** \n");
  
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

