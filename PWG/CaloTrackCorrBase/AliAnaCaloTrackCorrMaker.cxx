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

#include <cstdlib>

// --- ROOT system ---
#include <TClonesArray.h>
#include <TList.h>
#include <TH1F.h>
//#include <TObjectTable.h>
#include <TGeoGlobalMagField.h>

//---- AliRoot system ----
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAnaCaloTrackCorrBaseClass.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"

/// \cond CLASSIMP
ClassImp(AliAnaCaloTrackCorrMaker) ;
/// \endcond

//__________________________________________________
/// Default constructor.
/// Initialize parameters, pointers and arrays.
//__________________________________________________
AliAnaCaloTrackCorrMaker::AliAnaCaloTrackCorrMaker() :
TObject(),
fReader(0),                   fCaloUtils(0),
fOutputContainer(new TList ), fAnalysisContainer(new TList ),
fMakeHisto(kFALSE),           fMakeAOD(kFALSE),
fAnaDebug(0),                 fCuts(new TList),
fScaleFactor(-1),
fFillDataControlHisto(1),     fSumw2(0),
fCheckPtHard(0),
// Control histograms
fhNEventsIn(0),               fhNEvents(0),
fhNExoticEvents(0),           fhNEventsNoTriggerFound(0),
fhNPileUpEvents(0),           fhNPileUpEventsTriggerBC0(0),
fhXVertex(0),                 fhYVertex(0),                       fhZVertex(0),
fhXVertexExotic(0),           fhYVertexExotic(0),                 fhZVertexExotic(0),
fhPtHard(0),                  fhPtHardWeighted(0),
fhPileUpClusterMult(0),       fhPileUpClusterMultAndSPDPileUp(0),
fhTrackMult(0),
fhCentrality(0),              fhEventPlaneAngle(0),
fhNEventsWeighted(0),         fhTrackMultWeighted(0),
fhCentralityWeighted(0),      fhEventPlaneAngleWeighted(0),
fhNMergedFiles(0),            fhScaleFactor(0),
fhEMCalBCEvent(0),            fhEMCalBCEventCut(0),
fhTrackBCEvent(0),            fhTrackBCEventCut(0),
fhPrimaryVertexBC(0),         fhTimeStampFraction(0),
fhNPileUpVertSPD(0),          fhNPileUpVertTracks(0),
// EMCal trigger control histograms
fhClusterTriggerBC(0),                  fhClusterTriggerBCExotic(0),
fhClusterTriggerBCBadCell(0),           fhClusterTriggerBCBadCellExotic(0),
fhClusterTriggerBCBadCluster(0),        fhClusterTriggerBCBadClusterExotic(0),
fhClusterTriggerBCUnMatch(0),           fhClusterTriggerBCExoticUnMatch(0),
fhClusterTriggerBCBadCellUnMatch(0),    fhClusterTriggerBCBadCellExoticUnMatch(0),
fhClusterTriggerBCBadClusterUnMatch(0), fhClusterTriggerBCBadClusterExoticUnMatch(0),
fhClusterTriggerBCEventBC(0),           fhClusterTriggerBCEventBCUnMatch(0),
fhClusterTriggerBCExoticEventBC(0),     fhClusterTriggerBCExoticEventBCUnMatch(0)
{
  AliDebug(1,"*** Analysis Maker Constructor ***");
  
  for(Int_t i = 0; i < 3; i++)
  {
    fhClusterTriggerBCUnMatchReMatch      [0] = 0;
    fhClusterTriggerBCExoticUnMatchReMatch[0] = 0;
  }
  
  InitParameters();
}

//________________________________________________________________________________________
/// Copy constructor.
//________________________________________________________________________________________
AliAnaCaloTrackCorrMaker::AliAnaCaloTrackCorrMaker(const AliAnaCaloTrackCorrMaker & maker) :
TObject(),
fReader(),   //(new AliCaloTrackReader(*maker.fReader)),
fCaloUtils(),//(new AliCalorimeterUtils(*maker.fCaloUtils)),
fOutputContainer(new TList()), fAnalysisContainer(new TList()),
fMakeHisto(maker.fMakeHisto),  fMakeAOD(maker.fMakeAOD),
fAnaDebug(maker.fAnaDebug),    fCuts(new TList()),
fScaleFactor(maker.fScaleFactor),
fFillDataControlHisto(maker.fFillDataControlHisto),
fSumw2(maker.fSumw2),
fCheckPtHard(maker.fCheckPtHard),
fhNEventsIn(maker.fhNEventsIn),
fhNEvents(maker.fhNEvents),
fhNExoticEvents(maker.fhNExoticEvents),
fhNEventsNoTriggerFound(maker.fhNEventsNoTriggerFound),
fhNPileUpEvents(maker.fhNPileUpEvents),
fhNPileUpEventsTriggerBC0(maker.fhNPileUpEventsTriggerBC0),
fhXVertex(maker.fhXVertex),
fhYVertex(maker.fhYVertex),
fhZVertex(maker.fhZVertex),
fhXVertexExotic(maker.fhXVertexExotic),
fhYVertexExotic(maker.fhYVertexExotic),
fhZVertexExotic(maker.fhZVertexExotic),
fhPtHard(maker.fhPtHard),
fhPtHardWeighted(maker.fhPtHardWeighted),
fhPileUpClusterMult(maker.fhPileUpClusterMult),
fhPileUpClusterMultAndSPDPileUp(maker.fhPileUpClusterMultAndSPDPileUp),
fhTrackMult(maker.fhTrackMult),
fhCentrality(maker.fhCentrality),
fhEventPlaneAngle(maker.fhEventPlaneAngle),
fhNEventsWeighted(maker.fhNEventsWeighted),
fhTrackMultWeighted(maker.fhTrackMultWeighted),
fhCentralityWeighted(maker.fhCentralityWeighted),
fhEventPlaneAngleWeighted(maker.fhEventPlaneAngleWeighted),
fhNMergedFiles(maker.fhNMergedFiles),
fhScaleFactor(maker.fhScaleFactor),
fhEMCalBCEvent(maker.fhEMCalBCEvent),
fhEMCalBCEventCut(maker.fhEMCalBCEventCut),
fhTrackBCEvent(maker.fhTrackBCEvent),
fhTrackBCEventCut(maker.fhTrackBCEventCut),
fhPrimaryVertexBC(maker.fhPrimaryVertexBC),
fhTimeStampFraction(maker.fhTimeStampFraction),
fhNPileUpVertSPD(maker.fhNPileUpVertSPD),
fhNPileUpVertTracks(maker.fhNPileUpVertTracks),
fhClusterTriggerBC(maker.fhClusterTriggerBC),
fhClusterTriggerBCExotic(maker.fhClusterTriggerBCExotic),
fhClusterTriggerBCBadCell(maker.fhClusterTriggerBCBadCell),
fhClusterTriggerBCBadCellExotic(maker.fhClusterTriggerBCBadCellExotic),
fhClusterTriggerBCBadCluster(maker.fhClusterTriggerBCBadCluster),
fhClusterTriggerBCBadClusterExotic(maker.fhClusterTriggerBCBadClusterExotic),
fhClusterTriggerBCUnMatch(maker.fhClusterTriggerBCUnMatch),
fhClusterTriggerBCExoticUnMatch(maker.fhClusterTriggerBCExoticUnMatch),
fhClusterTriggerBCBadCellUnMatch(maker.fhClusterTriggerBCBadCellUnMatch),
fhClusterTriggerBCBadCellExoticUnMatch(maker.fhClusterTriggerBCBadCellExoticUnMatch),
fhClusterTriggerBCBadClusterUnMatch(maker.fhClusterTriggerBCBadClusterUnMatch),
fhClusterTriggerBCBadClusterExoticUnMatch(maker.fhClusterTriggerBCBadClusterExoticUnMatch),
fhClusterTriggerBCEventBC(maker.fhClusterTriggerBCEventBC),
fhClusterTriggerBCEventBCUnMatch(maker.fhClusterTriggerBCEventBCUnMatch),
fhClusterTriggerBCExoticEventBC(maker.fhClusterTriggerBCExoticEventBC),
fhClusterTriggerBCExoticEventBCUnMatch(maker.fhClusterTriggerBCExoticEventBCUnMatch)
{
  for(Int_t i = 0; i < 3; i++)
  {
    fhClusterTriggerBCUnMatchReMatch      [i] = maker.fhClusterTriggerBCUnMatchReMatch      [i];
    fhClusterTriggerBCExoticUnMatchReMatch[i] = maker.fhClusterTriggerBCExoticUnMatchReMatch[i];
  }
}

//___________________________________________________
/// Destructor. Remove only owned pointers.
//___________________________________________________
AliAnaCaloTrackCorrMaker::~AliAnaCaloTrackCorrMaker()
{
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
/// Add analysis depending on AliAnaCaloTrackCorrBaseClass to list.
//__________________________________________________________________
void AliAnaCaloTrackCorrMaker::AddAnalysis(TObject* ana, Int_t n)
{  
  if ( fAnalysisContainer)
  {
    fAnalysisContainer->AddAt(ana,n);
  }
  else
  {
    AliFatal("AnalysisContainer not initialized");
  }
}

//_________________________________________________________
/// \return list with any new output AOD branches from analysis.
/// The list is filled in the maker, and new branch passed to the 
/// main analysis task AliAnalysisTaskCaloTrackCorrelation.
//_________________________________________________________
TList * AliAnaCaloTrackCorrMaker::FillAndGetAODBranchList()
{
  TList *aodBranchList = fReader->GetAODBranchList() ;
  
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    if(ana->NewOutputAOD()) aodBranchList->Add(ana->GetCreateOutputAODBranch());
  }
  
  return aodBranchList ;
}

//____________________________________________________
/// Fill here event control histograms.
//____________________________________________________
void AliAnaCaloTrackCorrMaker::FillControlHistograms()
{  
  AliVEvent* event =  fReader->GetInputEvent();
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);
  
  fhNEvents        ->Fill(0); // Number of events analyzed
  
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v) ;
  fhXVertex->Fill(v[0]);
  fhYVertex->Fill(v[1]);
  fhZVertex->Fill(v[2]);
  
  fhTrackMult      ->Fill(fReader->GetTrackMultiplicity());
  fhCentrality     ->Fill(fReader->GetEventCentrality  ());
  fhEventPlaneAngle->Fill(fReader->GetEventPlaneAngle  ());
      
  if ( fReader->GetWeightUtils()->IsCentralityWeightOn() )
  {      
    Float_t eventWeight = fReader->GetEventWeight();

    fhNEventsWeighted        ->Fill(0.,                             eventWeight);
    fhTrackMultWeighted      ->Fill(fReader->GetTrackMultiplicity(),eventWeight);
    fhCentralityWeighted     ->Fill(fReader->GetEventCentrality  (),eventWeight);
    fhEventPlaneAngleWeighted->Fill(fReader->GetEventPlaneAngle  (),eventWeight);
  }
  
  // Check the pT hard in MC
  if ( fCheckPtHard && fReader->GetGenEventHeader()  )
  {
    if(!strcmp(fReader->GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
    {
      AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) fReader->GetGenEventHeader();
      
      Float_t pTHard = pygeh->GetPtHard();
      
      //printf("pT hard %f, event weight %e\n",pTHard,fReader->GetEventWeight());
      
      fhPtHard->Fill(pTHard);
      
      if ( fReader->GetWeightUtils()->IsMCCrossSectionCalculationOn() && 
          !fReader->GetWeightUtils()->IsMCCrossSectionJustHistoFillOn()  )     
        fhPtHardWeighted->Fill(pTHard,fReader->GetEventWeight());
    }
  }
  
  if(fFillDataControlHisto)
  {
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
    
    Int_t triggerBC = fReader->GetTriggerClusterBC() ;
    if( triggerBC == 0            &&
       !fReader->IsExoticEvent()  &&
       !fReader->IsBadCellTriggerEvent())
    {
      if( fReader->IsPileUpFromSPD())
        fhNPileUpEventsTriggerBC0->Fill(0.5);
      //if( event->IsPileupFromSPDInMultBins())
      //  fhNPileUpEventsTriggerBC0->Fill(1.5);
      if( fReader->IsPileUpFromEMCal())
        fhNPileUpEventsTriggerBC0->Fill(2.5);
      if( fReader->IsPileUpFromSPDOrEMCal() )
        fhNPileUpEventsTriggerBC0->Fill(3.5);
      if( fReader->IsPileUpFromSPDAndEMCal() )
        fhNPileUpEventsTriggerBC0->Fill(4.5);
      if( fReader->IsPileUpFromSPDAndNotEMCal() )
        fhNPileUpEventsTriggerBC0->Fill(5.5);
      if( fReader->IsPileUpFromEMCalAndNotSPD() )
        fhNPileUpEventsTriggerBC0->Fill(6.5);
      if( fReader->IsPileUpFromNotSPDAndNotEMCal() )
        fhNPileUpEventsTriggerBC0->Fill(7.5);
    }
    
    if(fReader->IsPileUpFromSPD())
      fhPileUpClusterMultAndSPDPileUp ->Fill(fReader->GetNPileUpClusters());
    
    fhPileUpClusterMult ->Fill(fReader->GetNPileUpClusters  ());
    
    for(Int_t i = 0; i < 19; i++)
    {
      if(fReader->IsAccessToTrackTimeOn())
      {
        if(fReader->GetTrackEventBC(i))   fhTrackBCEvent   ->Fill(i);
        if(fReader->GetTrackEventBCcut(i))fhTrackBCEventCut->Fill(i);
      }
      if(fReader->GetEMCalEventBC(i))   fhEMCalBCEvent   ->Fill(i);
      if(fReader->GetEMCalEventBCcut(i))fhEMCalBCEventCut->Fill(i);
    }
    
    Int_t bc = fReader->GetVertexBC();
    if(bc!=AliVTrack::kTOFBCNA) fhPrimaryVertexBC->Fill(bc);
    
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
}

//___________________________________________________________
/// Fill here EMCal triggered events control histograms.
//___________________________________________________________
void AliAnaCaloTrackCorrMaker::FillTriggerControlHistograms()
{
  if ( fFillDataControlHisto < 2 ) return;
  
  Int_t  triggerBC   = fReader->GetTriggerClusterBC() ;
  Bool_t exotic      = fReader->IsExoticEvent();
  Bool_t badCluster  = fReader->IsBadCellTriggerEvent();
  Bool_t badCell     = fReader->IsBadMaxCellTriggerEvent();
  Bool_t triggerMatch= fReader->IsTriggerMatched();
  Bool_t triggerBCOK = kTRUE;
  Int_t  triggerId   = fReader->GetTriggerClusterId() ;
  
  Bool_t reMatchOpenTime = fReader->IsTriggerMatchedOpenCuts(0);
  Bool_t reMatchNeigbour = fReader->IsTriggerMatchedOpenCuts(1);
  Bool_t reMatchBoth     = fReader->IsTriggerMatchedOpenCuts(2);
  
  if(triggerId < 0)
  {
    //printf("Trigger id %d\n",triggerId);
    if(triggerId == -2) fhNEventsNoTriggerFound->Fill(0);
    triggerBCOK = kFALSE;
  }
  
  if(exotic)
  {
    fhNExoticEvents->Fill(0) ;
    Double_t v[3];
    fReader->GetInputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
    fhXVertexExotic->Fill(v[0]);
    fhYVertexExotic->Fill(v[1]);
    fhZVertexExotic->Fill(v[2]);
  }
  //if(fReader->IsExoticEvent()) printf("Maker: EXOTIC Cluster trigger\n");
  
  if(!triggerBCOK) return;
  
  Int_t eventBC = fReader->GetInputEvent()->GetBunchCrossNumber();
  if(eventBC%4 < 0 || eventBC%4 > 3 )
    AliWarning(Form("STRANGE: Trigger BC %d - Event BC %d, modulo4 %d",triggerBC,eventBC,eventBC%4));
  
  if(triggerMatch)
  {
    if     (!exotic     && !badCluster) fhClusterTriggerBC->Fill(triggerBC);
    else if( exotic     &&  badCluster)
    {
      fhClusterTriggerBCBadClusterExotic->Fill(triggerBC);
      if(badCell) fhClusterTriggerBCBadCellExotic->Fill(triggerBC);
    }
    else if( exotic     && !badCluster) fhClusterTriggerBCExotic->Fill(triggerBC);
    else if( badCluster && !exotic    )
    {
      fhClusterTriggerBCBadCluster ->Fill(triggerBC);
      if(badCell) fhClusterTriggerBCBadCell->Fill(triggerBC);
    }
    
    if(!exotic) fhClusterTriggerBCEventBC      ->Fill(triggerBC,eventBC%4);
    else        fhClusterTriggerBCExoticEventBC->Fill(triggerBC,eventBC%4);
  }
  else
  {
    if     (!exotic     && !badCluster)
    {
      fhClusterTriggerBCUnMatch->Fill(triggerBC);
      if(reMatchOpenTime) fhClusterTriggerBCUnMatchReMatch[0]->Fill(triggerBC);
      if(reMatchNeigbour) fhClusterTriggerBCUnMatchReMatch[1]->Fill(triggerBC);
      if(reMatchBoth)     fhClusterTriggerBCUnMatchReMatch[2]->Fill(triggerBC);
    }
    else if( exotic     &&  badCluster)
    {
      fhClusterTriggerBCBadClusterExoticUnMatch->Fill(triggerBC);
      if(badCell)   fhClusterTriggerBCBadCellExoticUnMatch   ->Fill(triggerBC);
    }
    else if( exotic     && !badCluster)
    {
      fhClusterTriggerBCExoticUnMatch->Fill(triggerBC);
      if(reMatchOpenTime) fhClusterTriggerBCExoticUnMatchReMatch[0]->Fill(triggerBC);
      if(reMatchNeigbour) fhClusterTriggerBCExoticUnMatchReMatch[1]->Fill(triggerBC);
      if(reMatchBoth)     fhClusterTriggerBCExoticUnMatchReMatch[2]->Fill(triggerBC);
    }
    else if( badCluster && !exotic )
    {
      fhClusterTriggerBCBadClusterUnMatch->Fill(triggerBC);
      if(badCell)fhClusterTriggerBCBadCellUnMatch->Fill(triggerBC);
    }
    
    if(!exotic) fhClusterTriggerBCEventBCUnMatch      ->Fill(triggerBC,eventBC%4);
    else        fhClusterTriggerBCExoticEventBCUnMatch->Fill(triggerBC,eventBC%4);
  }
}

//_______________________________________________________
/// \return the list of the cuts (strings) used for the analysis.
/// The list is filled in the maker, called by the task in LocalInit() and posted there
//_______________________________________________________
TList * AliAnaCaloTrackCorrMaker::GetListOfAnalysisCuts()
{
  // Reader cuts
  TObjString * objstring = fReader->GetListOfParameters();
  fCuts->Add(objstring);
  
  // Analysis wagons cuts
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    objstring = ana->GetAnalysisCuts();
    
    if ( objstring ) fCuts->Add(objstring);
  }
  
  return fCuts ;
}

//___________________________________________________
/// Fill the output list of histograms during the CreateOutputObjects stage.
/// Get/create the histograms of the different analysis tasks and add them to the 
/// common output list. Event control histograms are also created and added here. 
//___________________________________________________
TList *AliAnaCaloTrackCorrMaker::GetOutputContainer()
{
  // General event histograms

  fhNEventsIn      = new TH1F("hNEventsIn",   "Number of input events"     , 1 , 0 , 1  ) ;
  fhNEventsIn->SetYTitle("# events");
  fOutputContainer->Add(fhNEventsIn);
  
  fhNEvents      = new TH1F("hNEvents",   "Number of analyzed events"     , 1 , 0 , 1  ) ;
  fhNEvents->SetYTitle("# events");
  fOutputContainer->Add(fhNEvents);
  
  fhXVertex      = new TH1F("hXVertex", " X vertex distribution"   , 200 , -1 , 1  ) ;
  fhXVertex->SetXTitle("v_{x} (cm)");
  fOutputContainer->Add(fhXVertex);
  
  fhYVertex      = new TH1F("hYVertex", " Y vertex distribution"   , 200 , -1 , 1  ) ;
  fhYVertex->SetXTitle("v_{y} (cm)");
  fOutputContainer->Add(fhYVertex);
  
  fhZVertex      = new TH1F("hZVertex", " Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZVertex);

  fhCentrality   = new TH1F("hCentrality","Number of events in centrality bin", 100, 0., 100) ;
  fhCentrality->SetXTitle("Centrality bin");
  fOutputContainer->Add(fhCentrality) ;
  
  fhEventPlaneAngle = new TH1F("hEventPlaneAngle","Number of events in event plane", 100, 0., TMath::Pi()) ;
  fhEventPlaneAngle->SetXTitle("EP angle (rad)");
  fOutputContainer->Add(fhEventPlaneAngle) ;

  fhTrackMult    = new TH1F
  ("hTrackMult", 
   Form("Number of tracks per event with #it{p}_{T} > %2.2f GeV/#it{c} and |#eta|<%2.2f",
        GetReader()->GetTrackMultiplicityPtCut(),GetReader()->GetTrackMultiplicityEtaCut()), 
   2000 , 0 , 2000) ;
  fhTrackMult->SetXTitle("# tracks");
  fOutputContainer->Add(fhTrackMult);
  
  if(fCheckPtHard)
  {
    fhPtHard  = new TH1F("hPtHard"," #it{p}_{T}-hard for selected triggers",300,0,300); 
    fhPtHard->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
    fOutputContainer->Add(fhPtHard);
    
    if ( fReader->GetWeightUtils()->IsMCCrossSectionCalculationOn() && 
        !fReader->GetWeightUtils()->IsMCCrossSectionJustHistoFillOn()  )     
    {
      fhPtHardWeighted  = new TH1F("hPtHardWeighted"," #it{p}_{T}-hard for selected triggers, weighted by cross section",300,0,300); 
      fhPtHardWeighted->SetXTitle("#it{p}_{T}^{hard} (GeV/#it{c})");
      fOutputContainer->Add(fhPtHardWeighted);
    }
  }
  
  if ( GetReader()->GetWeightUtils()->IsCentralityWeightOn() )
  {
    fhNEventsWeighted         = new TH1F("hNEventsWeighted",   "Number of analyzed events weighted by centrality", 1 , 0 , 1  ) ;
    fhNEventsWeighted->SetYTitle("# events");
    fOutputContainer->Add(fhNEventsWeighted);
      
    fhCentralityWeighted      = new TH1F("hCentralityWeighted","Number of events in centrality bin weighted by centrality", 100, 0.,100) ;
    fhCentralityWeighted->SetXTitle("Centrality bin");
    fOutputContainer->Add(fhCentralityWeighted) ;
      
    fhEventPlaneAngleWeighted = new TH1F("hEventPlaneAngleWeighted","Number of events in event plane weighted by centrality",100, 0., TMath::Pi()) ;
    fhEventPlaneAngleWeighted->SetXTitle("EP angle (rad)");
    fOutputContainer->Add(fhEventPlaneAngleWeighted) ;
      
    fhTrackMultWeighted       = new TH1F
    ("hTrackMultWeighted", 
     Form("Number of tracks per weighted event with #it{p}_{T} > %2.2f GeV/#it{c} and |#eta|<%2.2f",
          GetReader()->GetTrackMultiplicityPtCut(),GetReader()->GetTrackMultiplicityEtaCut()), 
     2000 , 0 , 2000) ;
    fhTrackMultWeighted->SetXTitle("# tracks");
    fOutputContainer->Add(fhTrackMultWeighted);
  }
    
  if(fFillDataControlHisto)
  {
    // Trigger and exoticity related histograms
    if(fFillDataControlHisto > 1)
    {
      fhNExoticEvents      = new TH1F("hNExoticEvents",   "Number of analyzed events triggered by exotic cluster"     , 1 , 0 , 1  ) ;
      fhNExoticEvents->SetYTitle("# exotic events");
      fOutputContainer->Add(fhNExoticEvents);
      
      fhNEventsNoTriggerFound      = new TH1F("hNEventsNoTriggerFound",   "Number of analyzed events triggered but no trigger found"     , 1 , 0 , 1  ) ;
      fhNEventsNoTriggerFound->SetYTitle("# exotic events");
      fOutputContainer->Add(fhNEventsNoTriggerFound);
      
      Int_t   nbin   = 11;
      Float_t minbin =-5.5;
      Float_t maxbin = 5.5;
      Int_t  labelshift = 6;
      
      fhClusterTriggerBCEventBC      = new TH2F("hClusterTriggerBCEventBC", "Found trigger BC and  Event BC",
                                                nbin , minbin ,maxbin,4,0, 4) ;
      fhClusterTriggerBCEventBC->SetXTitle("cluster trigger BC");
      for(Int_t i = 0; i < 4; i++)
        fhClusterTriggerBCEventBC->GetYaxis()->SetBinLabel(i+1 ,Form("BC/4=%d",i));
      fhClusterTriggerBCEventBC->SetXTitle("cluster trigger BC");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCEventBC->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fhClusterTriggerBCEventBC->SetYTitle("Event BC%4");
      fOutputContainer->Add(fhClusterTriggerBCEventBC);
      
      fhClusterTriggerBCExoticEventBC      = new TH2F("hClusterTriggerBCExoticEventBC", "Found exotic trigger BC and  Event BC",
                                                      nbin , minbin ,maxbin,4,1, 4) ;
      for(Int_t i = 0; i < 4; i++)
        fhClusterTriggerBCExoticEventBC->GetYaxis()->SetBinLabel(i+1 ,Form("BC/4=%d",i));
      fhClusterTriggerBCExoticEventBC->SetXTitle("cluster trigger BC");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCExoticEventBC->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fhClusterTriggerBCExoticEventBC->SetYTitle("Event BC%4");
      fOutputContainer->Add(fhClusterTriggerBCExoticEventBC);
      
      fhClusterTriggerBCEventBCUnMatch      = new TH2F("hClusterTriggerBCEventBCUnMatch", "Found unmatched trigger BC and  Event BC",
                                                       nbin , minbin ,maxbin,4,1, 4) ;
      for(Int_t i = 0; i < 4; i++)
        fhClusterTriggerBCEventBCUnMatch->GetYaxis()->SetBinLabel(i+1 ,Form("BC/4=%d",i));
      fhClusterTriggerBCEventBCUnMatch->SetXTitle("cluster trigger BC");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCEventBCUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fhClusterTriggerBCEventBCUnMatch->SetYTitle("Event BC%4");
      fOutputContainer->Add(fhClusterTriggerBCEventBCUnMatch);
      
      fhClusterTriggerBCExoticEventBCUnMatch      = new TH2F("hClusterTriggerExoticBCEventBCUnMatch", "Found unmatched trigger BC and  Event BC",
                                                             nbin , minbin ,maxbin,4,1, 4) ;
      for(Int_t i = 0; i < 4; i++)
        fhClusterTriggerBCExoticEventBCUnMatch->GetYaxis()->SetBinLabel(i+1 ,Form("BC/4=%d",i));
      fhClusterTriggerBCExoticEventBCUnMatch->SetXTitle("cluster trigger BC");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCExoticEventBCUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fhClusterTriggerBCExoticEventBCUnMatch->SetYTitle("Event BC%4");
      fOutputContainer->Add(fhClusterTriggerBCExoticEventBCUnMatch);
      
      fhClusterTriggerBC              = new TH1F("hClusterTriggerBC",
                                                 "Number of analyzed events triggered by a cluster in a given BC",
                                                 nbin , minbin ,maxbin) ;
      fhClusterTriggerBC->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBC->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBC);
      
      fhClusterTriggerBCExotic        = new TH1F("hClusterTriggerBCExotic",
                                                 "Number of analyzed events triggered by a exotic cluster in a given BC",
                                                 nbin , minbin ,maxbin) ;
      fhClusterTriggerBCExotic->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCExotic->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCExotic);
      
      
      fhClusterTriggerBCBadCell         = new TH1F("hClusterTriggerBCBadCell",
                                                   "Number of analyzed events triggered by a bad cell in a given BC",
                                                   nbin , minbin ,maxbin) ;
      
      fhClusterTriggerBCBadCell->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadCell->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadCell);
      
      fhClusterTriggerBCBadCellExotic    = new TH1F("hClusterTriggerBCBadCellExotic",
                                                    "Number of analyzed events triggered by a bad cell & exotic cluster in a given BC",
                                                    nbin , minbin ,maxbin) ;
      fhClusterTriggerBCBadCellExotic->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadCellExotic->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadCellExotic);
      
      fhClusterTriggerBCBadCluster           = new TH1F("hClusterTriggerBCBadCluster",
                                                        "Number of analyzed events triggered by a bad cluster in a given BC",
                                                        nbin , minbin ,maxbin) ;
      
      fhClusterTriggerBCBadCluster->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadCluster->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadCluster);
      
      
      fhClusterTriggerBCBadClusterExotic    = new TH1F("hClusterTriggerBCBadClusterExotic",
                                                       "Number of analyzed events triggered by a bad cluster & exotic cluster in a given BC",
                                                       nbin , minbin ,maxbin) ;
      
      fhClusterTriggerBCBadClusterExotic->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadClusterExotic->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadClusterExotic);
      
      fhClusterTriggerBCUnMatch       = new TH1F("hClusterTriggerBCUnMatch",
                                                 "Number of analyzed events triggered by a cluster (no trigger patch match) in a given BC",
                                                 nbin , minbin ,maxbin) ;
      fhClusterTriggerBCUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCUnMatch);
      
      fhClusterTriggerBCExoticUnMatch = new TH1F("hClusterTriggerBCExoticUnMatch",
                                                 "Number of analyzed events triggered by a exotic cluster (no trigger patch match) in a given BC",
                                                 nbin , minbin ,maxbin) ;
      fhClusterTriggerBCExoticUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCExoticUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCExoticUnMatch);
      
      
      fhClusterTriggerBCBadCellUnMatch    = new TH1F("hClusterTriggerBCBadCellUnMatch",
                                                     "Number of analyzed events triggered by a bad cluster  (no trigger patch match) in a given BC",
                                                     nbin , minbin ,maxbin) ;
      fhClusterTriggerBCBadCellUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadCellUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadCellUnMatch);
      
      
      fhClusterTriggerBCBadCellExoticUnMatch = new TH1F("hClusterTriggerBCBadCellExoticUnMatch",
                                                        "Number of analyzed events triggered by a bad&exotic cluster  (no trigger patch match) in a given BC",
                                                        nbin , minbin ,maxbin) ;
      fhClusterTriggerBCBadCellExoticUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadCellExoticUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadCellExoticUnMatch);
      
      
      fhClusterTriggerBCBadClusterUnMatch    = new TH1F("hClusterTriggerBCBadClusterUnMatch",
                                                        "Number of analyzed events triggered by a bad cluster  (no trigger patch match) in a given BC",
                                                        nbin , minbin ,maxbin) ;
      fhClusterTriggerBCBadClusterUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadClusterUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadClusterUnMatch);
      
      
      fhClusterTriggerBCBadClusterExoticUnMatch = new TH1F("hClusterTriggerBCBadClusterExoticUnMatch",
                                                           "Number of analyzed events triggered by a bad&exotic cluster  (no trigger patch match) in a given BC",
                                                           nbin , minbin ,maxbin) ;
      fhClusterTriggerBCBadClusterExoticUnMatch->SetYTitle("# events");
      for(Int_t i = 1; i < 12; i++)
        fhClusterTriggerBCBadClusterExoticUnMatch->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
      fOutputContainer->Add(fhClusterTriggerBCBadClusterExoticUnMatch);
      
      TString rematch[] = {"OpenTime","CheckNeighbours","Both"};
      for(Int_t j = 0; j < 3; j++)
      {
        fhClusterTriggerBCUnMatchReMatch[j]       = new TH1F(Form("hClusterTriggerBCUnMatch_ReMatch_%s",rematch[j].Data()),
                                                             Form("Number of analyzed events triggered by a cluster (no trigger patch match) in a given BC, re-match %s",rematch[j].Data()),
                                                             nbin , minbin ,maxbin) ;
        fhClusterTriggerBCUnMatchReMatch[j]->SetYTitle("# events");
        for(Int_t i = 1; i < 12; i++)
          fhClusterTriggerBCUnMatchReMatch[j]->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
        fOutputContainer->Add(fhClusterTriggerBCUnMatchReMatch[j]);
        
        fhClusterTriggerBCExoticUnMatchReMatch[j] = new TH1F(Form("hClusterTriggerBCExoticUnMatch_ReMatch_%s",rematch[j].Data()),
                                                             Form("Number of analyzed events triggered by a exotic cluster (no trigger patch match) in a given BC, re-match %s",rematch[j].Data()),
                                                             nbin , minbin ,maxbin) ;
        fhClusterTriggerBCExoticUnMatchReMatch[j]->SetYTitle("# events");
        for(Int_t i = 1; i < 12; i++)
          fhClusterTriggerBCExoticUnMatchReMatch[j]->GetXaxis()->SetBinLabel(i ,Form("BC%d",i-labelshift));
        fOutputContainer->Add(fhClusterTriggerBCExoticUnMatchReMatch[j]);
      }
      
      fhXVertexExotic      = new TH1F("hXVertexExotic", " X vertex distribution in exotic events"   , 200 , -1 , 1  ) ;
      fhXVertexExotic->SetXTitle("v_{x} (cm)");
      fOutputContainer->Add(fhXVertexExotic);
      
      fhYVertexExotic      = new TH1F("hYVertexExotic", " Y vertex distribution in exotic events"   , 200 , -1 , 1  ) ;
      fhYVertexExotic->SetXTitle("v_{y} (cm)");
      fOutputContainer->Add(fhYVertexExotic);
      
      fhZVertexExotic      = new TH1F("hZVertexExotic", " Z vertex distribution in exotic events"   , 200 , -50 , 50  ) ;
      fhZVertexExotic->SetXTitle("v_{z} (cm)");
      fOutputContainer->Add(fhZVertexExotic);
    }
    
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
    
    fhNPileUpEventsTriggerBC0      = new TH1F("hNPileUpEventsTriggerBC0","Number of events considered as pile-up, trigger cluster in BC=0", 8 , 0 , 8 ) ;
    fhNPileUpEventsTriggerBC0->SetYTitle("# events");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(1 ,"SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(2 ,"Multi SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(3 ,"EMCal");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(4 ,"EMCal || SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(5 ,"EMCal && SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(6 ,"!EMCal && SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(7 ,"EMCal && !SPD");
    fhNPileUpEventsTriggerBC0->GetXaxis()->SetBinLabel(8 ,"!EMCal && !SPD");
    fOutputContainer->Add(fhNPileUpEventsTriggerBC0);
    
    if(fReader->IsAccessToTrackTimeOn())
    {
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
    }
    
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
    
    if(fReader->IsSelectEventTimeStampOn())
    {
      fhTimeStampFraction = new TH1F("hTimeStampFraction","Fraction of events within a given time stamp range",150, -1, 2) ;
      fhTimeStampFraction->SetXTitle("fraction");
      fOutputContainer->Add(fhTimeStampFraction) ;
    }
  }
  
  fhNMergedFiles = new TH1F("hNMergedFiles",   "Number of merged output files"     , 1 , 0 , 1  ) ;
  fhNMergedFiles->SetYTitle("# files");
  fhNMergedFiles->Fill(1); // Fill here with one entry, while merging it will count the rest
  fOutputContainer->Add(fhNMergedFiles);
  
  if(fScaleFactor > 0)
  {    
    fhScaleFactor = new TH1F("hScaleFactor",   "Number of merged output files"     , 1 , 0 , 1  ) ;
    fhScaleFactor->SetYTitle("scale factor");
    fhScaleFactor->SetBinContent(1,fScaleFactor); // Fill here
    fOutputContainer->Add(fhScaleFactor);
  }
  
  // Histograms defined and filled in this class, just get the pointers
  // and add them to the list.
  if(GetReader()->GetWeightUtils()->IsMCCrossSectionCalculationOn())
  {
    TList * templist =  GetReader()->GetWeightUtils()->GetCreateOutputHistograms();
      
    if ( templist && templist->GetEntries() == 2 )
    {
      templist->SetOwner(kFALSE); //Owner is fOutputContainer.
    
      fOutputContainer->Add(templist->At(0));
      fOutputContainer->Add(templist->At(1));
    }
  }
  
  // --------------------------------
  // Add control histograms in Reader
  // --------------------------------
  
  TList * templist =  fReader->GetCreateControlHistograms();
  templist->SetOwner(kFALSE); //Owner is fOutputContainer.
    
  for(Int_t ih = 0; ih < templist->GetEntries() ; ih++)
  {        
    //if ( fSumw2 ) ((TH1*) templist->At(ih))->Sumw2();
    
    //printf("histo %d %p %s\n",ih,templist->At(ih), templist->At(ih)->GetName());
    
    //Add histogram to general container
    fOutputContainer->Add(templist->At(ih)) ;
  }
  
  delete templist;
  
  // ------------------------
  // Add analysis histograms
  // ------------------------

  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0)
  {
    AliWarning("Analysis job list not initialized!!!");
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
      templist =  ana->GetCreateOutputObjects();
      templist->SetOwner(kFALSE); //Owner is fOutputContainer.
      
      for(Int_t i = 0; i < templist->GetEntries(); i++)
      {
        //Add only  to the histogram name the name of the task
        if(   strcmp((templist->At(i))->ClassName(),"TObjString")   )
        {
          snprintf(newname,buffersize, "%s%s", (ana->GetAddedHistogramsStringToName()).Data(), (templist->At(i))->GetName());
          //printf("name %s, new name %s\n",(templist->At(i))->GetName(),newname);
          
          ((TH1*) templist->At(i))->SetName(newname);
            
          if ( fSumw2 ) ((TH1*) templist->At(i))->Sumw2();
        }
        
        //Add histogram to general container
        fOutputContainer->Add(templist->At(i)) ;
      }
      
      delete templist;
      
    }// Analysis with histograms as output on
    
  }//Loop on analysis defined
  
  // Initialize calorimeters  geometry pointers
  //GetCaloUtils()->InitPHOSGeometry();
  //GetCaloUtils()->InitEMCALGeometry();
  
  return fOutputContainer;
}

//___________________________________
/// Init container histograms and other common variables
/// Fill the output list of histograms during the CreateOutputObjects stage.
/// Pass to the list of analysis classes the pointers to AliCaloTrackReader and AliCalorimeterUtils, 
/// and call the Init of those analysis.
//___________________________________
void AliAnaCaloTrackCorrMaker::Init()
{  
  // Activate debug level in maker
  if( fAnaDebug >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(this->ClassName(),fAnaDebug);
  
  //Initialize reader
  GetReader()->Init();
  GetReader()->SetCaloUtils(GetCaloUtils()); // pass the calo utils pointer to the reader

  // Activate debug level in reader
  if( fReader->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(fReader->ClassName(), fReader->GetDebug());
  
  // Activate debug level in calo utils
  if( fCaloUtils->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(fCaloUtils->ClassName(), fCaloUtils->GetDebug());
  
  if(!fAnalysisContainer || fAnalysisContainer->GetEntries()==0)
  {
    AliWarning("Analysis job list not initialized");
    return;
  }
  
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    
    ana->SetReader(fReader);       // Set Reader for each analysis
    ana->SetCaloUtils(fCaloUtils); // Set CaloUtils for each analysis
    
    ana->Init();
    ana->InitDebug();
  }//Loop on analysis defined
}

//_____________________________________________
/// Init data members
//_____________________________________________
void AliAnaCaloTrackCorrMaker::InitParameters()
{  
  fMakeHisto  = kTRUE;
  fMakeAOD    = kTRUE;
  fAnaDebug   = 0; // No debugging info displayed by default	
}

//______________________________________________________________
/// Print some relevant parameters set for the analysis
//______________________________________________________________
void AliAnaCaloTrackCorrMaker::Print(const Option_t * opt) const
{	
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

//_____________________________________________________________________________________
/// Main method, analysis are executed here:
/// * 1) Clean-up arrays and stuff, access OADB (once), geometry (once), etc.
/// * 2) Event selection, track filtering and cluster filtering in AliCaloTrackReader. 
/// * 3) Process different defined analysis for this event. 
/// * 4) Event selection for mixing studies.
/// * 5) Control histograms filling
/// * 6) Loop on the defined analysis list and execute the analysis defined steps
/// * 7) Reset lists, fill more control histograms
/// 
//_____________________________________________________________________________________
void AliAnaCaloTrackCorrMaker::ProcessEvent(Int_t iEntry, const char * currentFileName)
{  
  if(fMakeHisto && !fOutputContainer)
    AliFatal("Histograms not initialized");
  
  AliDebug(1,Form("***  AliAnaCaloTrackCorrMaker::ProcessEvent() Event %d   ***",iEntry));
  AliDebug(2,Form("Current File Name : %s", currentFileName));
  //printf("fAODBranchList %p, entries %d\n",fAODBranchList,fAODBranchList->GetEntries());
  
  // Each event needs an empty branch
  TList * aodList = fReader->GetAODBranchList();
  Int_t nAODBranches = aodList->GetEntries();
  for(Int_t iaod = 0; iaod < nAODBranches; iaod++)
  {
	  TClonesArray *tca = dynamic_cast<TClonesArray*> (aodList->At(iaod));
	  if(tca) tca->Clear("C");
  }
  
  // Set geometry matrices before filling arrays, in case recalibration/position calculation etc is needed
  fCaloUtils->AccessGeometry(fReader->GetInputEvent());
  
  // Set the AODB calibration, bad channels etc. parameters at least once
  fCaloUtils->AccessOADB(fReader->GetInputEvent());
  
  // Tell the reader to fill the data in the 3 detector lists
  Bool_t ok = fReader->FillInputEvent(iEntry, currentFileName);
  
  // Access pointers, and trigger mask check needed in mixing case
  AliAnalysisManager   *manager      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
  
  UInt_t isMBTrigger = kFALSE;
  UInt_t isTrigger   = kFALSE;
  if(inputHandler)
  {
    isMBTrigger = inputHandler->IsEventSelected() & fReader->GetMixEventTriggerMask();
    isTrigger   = inputHandler->IsEventSelected() & fReader->GetEventTriggerMask();
  }
  
  //Fill trigger control histograms, make sure it is only for triggered events and
  // not the MB events used for mixing
  if(fReader->IsEventTriggerAtSEOn() || isTrigger)
  {
    fhNEventsIn->Fill(0);
    FillTriggerControlHistograms();
  }
  
  if(!ok)
  {
    AliDebug(1,Form("*** Skip event *** %d",iEntry));
    fReader->ResetLists();
    return ;
  }
  
  // Magic line to write events to file
  if(fReader->WriteDeltaAODToFile())AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
  
  //printf(">>>>>>>>>> BEFORE >>>>>>>>>>>\n");
  //gObjectTable->Print();
  
  // Init mag field for tracks in case of ESDs, not really necessary
  if ( !TGeoGlobalMagField::Instance()->GetField() && fReader->GetInputEvent() )
      (fReader->GetInputEvent())->InitMagneticField();
  
  // Loop on analysis algorithms
  
  AliDebug(1,"*** Begin analysis ***");
  
  Int_t nana = fAnalysisContainer->GetEntries() ;
  for(Int_t iana = 0; iana <  nana; iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana = ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    
    ana->ConnectInputOutputAODBranches(); // Sets branches for each analysis
    
    //Fill pool for mixed event for the analysis that need it
    if(!fReader->IsEventTriggerAtSEOn() && isMBTrigger)
    {
      ana->FillEventMixPool();
      if(!isTrigger) continue; // pool filled do not try to fill AODs or histograms if trigger is not MB
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
    AliDebug(1,"*** End analysis, MB for mixing ***");
    return;
  }
  
  FillControlHistograms();
  
  //printf(">>>>>>>>>> AFTER >>>>>>>>>>>\n");
  //gObjectTable->Print();
	
  AliDebug(1,"*** End analysis ***");
}

//__________________________________________________________
/// Execute Terminate of analysis.
/// Do some final plots.
//__________________________________________________________
void AliAnaCaloTrackCorrMaker::Terminate(TList * outputList)
{
  if (!outputList)
  {
    AliError("No output list");
    return;
  }
  
  for(Int_t iana = 0; iana <  fAnalysisContainer->GetEntries(); iana++)
  {
    AliAnaCaloTrackCorrBaseClass * ana =  ((AliAnaCaloTrackCorrBaseClass *) fAnalysisContainer->At(iana)) ;
    if(ana->MakePlotsOn())ana->Terminate(outputList);
  } // Loop on analysis defined
}



