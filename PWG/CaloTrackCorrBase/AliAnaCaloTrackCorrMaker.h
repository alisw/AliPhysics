#ifndef ALIANACALOTRACKCORRMAKER_H
#define ALIANACALOTRACKCORRMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_____________________________________________________________________________
// Steering class for particle (gamma, hadron) identification and correlation 
// analysis. It is called by the task class AliAnalysisTaskCaloTrackCorrelation 
// and it connects the input (ESD/AOD/MonteCarlo) got with AliCaloTrackReader 
// (produces TClonesArrays of AODs (TParticles in MC case if requested)), with 
// the analysis classes that derive from AliAnaCaloTrackCorrBaseClass
//
// -- Author: Gustavo Conesa (INFN-LNF, LPSC-Grenoble)

// --- ROOT system ---
class TList; 
class TClonesArray;
#include<TObject.h>
class TH1F;

// --- Analysis system ---
#include "AliCaloTrackReader.h" 
#include "AliCalorimeterUtils.h"

class AliAnaCaloTrackCorrMaker : public TObject {

 public: 
  
  AliAnaCaloTrackCorrMaker() ;          // default ctor
  virtual ~AliAnaCaloTrackCorrMaker() ; // virtual dtor
  AliAnaCaloTrackCorrMaker(const AliAnaCaloTrackCorrMaker & maker) ; // cpy ctor
	
  // Setters and Getters
  
  void    AddAnalysis(TObject* ana, Int_t n) ;

  void    FillControlHistograms();
  
  void    FillTriggerControlHistograms();
  
  TList * GetListOfAnalysisContainers() { return fAnalysisContainer ; }
  
  TList * GetListOfAnalysisCuts();
  
  TList * GetOutputContainer() ;
  
  TList * FillAndGetAODBranchList();
  
  Int_t   GetAnaDebug()           const { return fAnaDebug    ; }
  void    SetAnaDebug(Int_t d)          { fAnaDebug = d       ; }
	
  Bool_t  AreHistogramsMade()     const { return fMakeHisto   ; }
  void    SwitchOnHistogramsMaker()     { fMakeHisto = kTRUE  ; }
  void    SwitchOffHistogramsMaker()    { fMakeHisto = kFALSE ; }
 
  Bool_t  AreAODsMade()           const { return fMakeAOD     ; }
  void    SwitchOnAODsMaker()           { fMakeAOD = kTRUE    ; }
  void    SwitchOffAODsMaker()          { fMakeAOD = kFALSE   ; }
  	
  void    SwitchOnDataControlHistograms()  { fFillDataControlHisto = kTRUE  ; }
  void    SwitchOffDataControlHistograms() { fFillDataControlHisto = kFALSE ; }

  
  AliCaloTrackReader  * GetReader()                                   { if(!fReader) fReader = new AliCaloTrackReader ();
                                                                        return fReader    ; }
  void                  SetReader(AliCaloTrackReader * reader)        { fReader = reader  ; }
  	
  AliCalorimeterUtils * GetCaloUtils()                                { if(!fCaloUtils) fCaloUtils = new AliCalorimeterUtils(); 
                                                                        return fCaloUtils      ; }
  void                  SetCaloUtils(AliCalorimeterUtils * caloutils) { fCaloUtils = caloutils ; }
	
  void                  SetScaleFactor(Double_t scale)                { fScaleFactor = scale   ; } 

  
  // Main general methods
  
  void    Init();
  
  void    InitParameters();
  
  void    Print(const Option_t * opt) const;
  
  void    ProcessEvent(Int_t iEntry, const char * currentFileName) ;
  
  void    Terminate(TList * outputList);
  
 private:
  
  //General Data members
  
  AliCaloTrackReader  *  fReader ;    //  Pointer to reader 
  AliCalorimeterUtils *  fCaloUtils ; //  Pointer to CalorimeterUtils
  
  TList *  fOutputContainer ;   //! Output histograms container
  TList *  fAnalysisContainer ; //  List with analysis pointers
  Bool_t   fMakeHisto ;         //  If true makes final analysis with histograms as output
  Bool_t   fMakeAOD ;           //  If true makes analysis generating AODs
  Int_t    fAnaDebug;           //  Debugging info.
  TList *  fCuts ;	            //! List with analysis cuts
  Double_t fScaleFactor ;       //  Scaling factor needed for normalization
  Bool_t   fFillDataControlHisto;//! Fill histograms only interesting with data
  
  // Control histograms
  TH1F *   fhNEventsIn;         //! Number of input events counter histogram
  TH1F *   fhNEvents;           //! Number of acepted events counter histogram
  TH1F *   fhNExoticEvents;     //! Number of events triggered by exotic, counter histogram
  TH1F *   fhNEventsNoTriggerFound; //! number of events where whaetever was done, no trigger is found
  TH1F *   fhNPileUpEvents;     //! N events pasing pile up cut
  TH1F *   fhNPileUpEventsTriggerBC0; //! N events pasing pile up cut
  TH1F *   fhXVertex;           //! X Vertex distribution of accepted event
  TH1F *   fhYVertex;           //! Y Vertex distribution of accepted event
  TH1F *   fhZVertex;           //! Z Vertex distribution of accepted event
  TH1F *   fhXVertexExotic;     //! X Vertex distribution of exotic event
  TH1F *   fhYVertexExotic;     //! Y Vertex distribution of exotic event
  TH1F *   fhZVertexExotic;     //! Z Vertex distribution of exotic event
  
  TH1F *   fhPileUpClusterMult; //! N clusters with high time
  TH1F *   fhPileUpClusterMultAndSPDPileUp; //! N clusters with high time in events tagged as pile-up by SPD
  TH1F *   fhTrackMult;         //! Number of tracks per event histogram
  TH1F *   fhCentrality;        //! Histogram with centrality bins
  TH1F *   fhEventPlaneAngle;   //! Histogram with Event plane angle
  TH1F *   fhNMergedFiles;      //! Number of files merged
  TH1F *   fhScaleFactor;       //! Factor to scale histograms
  TH1F *   fhEMCalBCEvent;      //! N events depending on the existence of a cluster in a given bunch crossing
  TH1F *   fhEMCalBCEventCut;   //! N events depending on the existence of a cluster above acceptance and E cut in a given bunch crossing
  TH1F *   fhTrackBCEvent;      //! N events depending on the existence of a track in a given bunch crossing
  TH1F *   fhTrackBCEventCut;   //! N events depending on the existence of a track above acceptance and pt cut in a given bunch crossing
  TH1F *   fhPrimaryVertexBC;   //! Primary vertex BC
  TH1F *   fhTimeStampFraction; //! event fraction depending on Time Stamp, only if activated on reader
  TH1F *   fhNPileUpVertSPD;    //! number of pile-up vertices from SPD
  TH1F *   fhNPileUpVertTracks; //! number of pile-up vertices from tracks
  
  TH1F *   fhClusterTriggerBC;                       //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCExotic;                 //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadCell;                //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadCellExotic;          //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadCluster;             //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadClusterExotic;       //! number of events triggered, depending on BC of the cluster
  
  TH1F *   fhClusterTriggerBCUnMatch;                //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCExoticUnMatch;          //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadCellUnMatch;         //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadCellExoticUnMatch;   //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadClusterUnMatch;      //! number of events triggered, depending on BC of the cluster
  TH1F *   fhClusterTriggerBCBadClusterExoticUnMatch;//! number of events triggered, depending on BC of the cluster

  TH1F *   fhClusterTriggerBCUnMatchReMatch[3];       //! number of events triggered, depending on BC of the cluster, not matched, open cuts and rematch
  TH1F *   fhClusterTriggerBCExoticUnMatchReMatch[3]; //! number of events triggered by exotic, depending on BC of the clusterm not matched, open cuts and rematch

  
  TH2F *   fhClusterTriggerBCEventBC;                // correlate the found BC in the trigger and the event BC
  TH2F *   fhClusterTriggerBCEventBCUnMatch;         // correlate the found BC in the trigger and the event BC, when there was no match with the trigger BC
  TH2F *   fhClusterTriggerBCExoticEventBC;          // correlate the found BC in the exotic trigger and the event BC
  TH2F *   fhClusterTriggerBCExoticEventBCUnMatch;   // correlate the found BC in the exotic trigger and the event BC, when there was no match with the trigger BC
  
  AliAnaCaloTrackCorrMaker & operator = (const AliAnaCaloTrackCorrMaker & ) ; // cpy assignment
  
  ClassDef(AliAnaCaloTrackCorrMaker,22)
} ;
 

#endif //ALIANACALOTRACKCORRMAKER_H



