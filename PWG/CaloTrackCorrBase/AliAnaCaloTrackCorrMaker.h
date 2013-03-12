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
  
  void    ProcessEvent(const Int_t iEntry, const char * currentFileName) ;
  
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

  // Control histograms
  TH1F *   fhNEvents;           //! Number of events counter histogram
  TH1F *   fhNPileUpEvents;     //! N events pasing pile up cut
  TH1F *   fhZVertex;           //! Vertex of accepted event
  TH1F *   fhPileUpClusterMult; //! N clusters with high time
  //TH1F *   fhPileUpClusterMultAndSPDPileUp; //! N clusters with high time in events tagged as pile-up by SPD
  TH1F *   fhTrackMult;         //! Number of tracks per event histogram
  TH1F *   fhCentrality;        //! Histogram with centrality bins
  TH1F *   fhEventPlaneAngle;   //! Histogram with Event plane angle
  TH1F *   fhNMergedFiles;      //! Number of files merged
  TH1F *   fhScaleFactor;       //! Factor to scale histograms
  TH1F *   fhEMCalBCEvent;      //! N events depending on the existance of a cluster in a given bunch crossing
  TH1F *   fhEMCalBCEventCut;   //! N events depending on the existance of a cluster above acceptance and E cut in a given bunch crossing
  TH1F *   fhTrackBCEvent;      //! N events depending on the existance of a track in a given bunch crossing
  TH1F *   fhTrackBCEventCut;   //! N events depending on the existance of a track above acceptance and pt cut in a given bunch crossing
  TH1F *   fhPrimaryVertexBC;   //! Primary vertex BC
  TH1F *   fhTimeStampFraction; //! event fraction depending on Time Stamp, only if activated on reader
  TH1F *   fhNPileUpVertSPD;    //! number of pile-up vertices from SPD
  TH1F *   fhNPileUpVertTracks; //! number of pile-up vertices from tracks
  
  AliAnaCaloTrackCorrMaker & operator = (const AliAnaCaloTrackCorrMaker & ) ; // cpy assignment
  
  ClassDef(AliAnaCaloTrackCorrMaker,17)
} ;
 

#endif //ALIANACALOTRACKCORRMAKER_H



