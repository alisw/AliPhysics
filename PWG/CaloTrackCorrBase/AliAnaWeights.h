#ifndef ALIANAWEIGHTS_H
#define ALIANAWEIGHTS_H

//_________________________________________________________________________
/// \class AliAnaWeights
/// \brief Calculate the weight to the event to be applied when filling histograms
///
/// This class is intended to get/calculate the global event weight needed to fill histograms.
/// Right now it only gets the weight corresponding to:
///  * the cross section of MC PYTHIA productions produced in pT-hard bins.
///  * in case of centrality dependent weight.
///
/// This task is called by the class AliCaloTrackReader in method
/// AliCaloTrackReader::FillInputEvent().
/// Also, in case of MC pT-hard bins analysis, the cross section and trials used are stored in
/// dedicated histograms via AliAnaCaloTrackCorrMaker, being initialized in
/// AliAnaCaloTrackCorrMaker::GetOutputContainer()
///
/// The different analysis subtasks should recover the weight stored in the reader,
/// and apply it when filling the histogram.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///         
//_________________________________________________________________________

class TH1F;
class TList;
class AliGenPythiaEventHeader;

#include <TObject.h>

class AliAnaWeights : public TObject {
  
 public:
  
  AliAnaWeights();  
    
  /// Destructor
  virtual         ~AliAnaWeights() { if ( fhCentralityWeight ) delete fhCentralityWeight; }
  
  //
  // Centrality weights
  //
    
  void             InitCentralityWeightsHistogram(Int_t nbins = 100, Int_t minCen = 0, Int_t maxCen = 100) ;
    
  TH1F *           GetCentralityWeightsHistogram() ;
    
  void             SetCentralityWeightsHistogram(TH1F* h)  { if ( fhCentralityWeight ) delete fhCentralityWeight ; fhCentralityWeight = h ; }
    
  void             SetCentrality(Float_t cen)              { fCentrality = cen             ; }
   
  Bool_t           IsCentralityWeightOn()            const { return fUseCentralityWeight   ; }
    
  void             SwitchOnCentralityWeight()              { fUseCentralityWeight = kTRUE  ; InitCentralityWeightsHistogram() ;  }
    
  void             SwitchOffCentralityWeight()             { fUseCentralityWeight = kFALSE ; }
    
  //
  // MC weights
  //
    
  TList *          GetCreateOutputHistograms();

  virtual Double_t GetPythiaCrossSection() ;

  static  Bool_t   GetPythiaInfoFromFile(TString currFile, Float_t & xsec, Float_t & trials) ;
    
  virtual Double_t GetWeight() ;
    
  Bool_t           IsMCCrossSectionCalculationOn()   const { return fCheckMCCrossSection   ; }
  Bool_t           IsMCCrossSectionJustHistoFillOn() const { return fJustFillCrossSecHist  ; }
    
  void             SwitchOnMCCrossSectionCalculation()     { fCheckMCCrossSection = kTRUE  ; }
  void             SwitchOffMCCrossSectionCalculation()    { fCheckMCCrossSection = kFALSE ; }

  void             SwitchOnMCCrossSectionHistoFill()       { fCheckMCCrossSection = kTRUE  ;  fJustFillCrossSecHist = kTRUE ; }
  
  void             SetPythiaEventHeader(AliGenPythiaEventHeader* py) { fPyEventHeader = py    ; }

  void             SwitchOnMCCrossSectionFromEventHeader() { fCheckPythiaEventHeader = kTRUE  ; }
  void             SwitchOffMCCrossSectionFromEventHeader(){ fCheckPythiaEventHeader = kFALSE ; }

  
  //
  // General
  //
    
  Bool_t           IsWeightSettingOn()               const { return ( fUseCentralityWeight || fCheckMCCrossSection ) ; }

  Int_t            GetDebug()                        const { return fDebug                 ; }
    
  void             SetDebug(Int_t d)                       { fDebug = d                    ; }

    
 private:
    
  Int_t            fDebug ;               ///< Debug level.
  
  //
  // Centrality weights
  //

  TH1F *           fhCentralityWeight ;   ///<  Container of centrality weights.
  
  Float_t          fCentrality ;          ///<  Container of centrality percentile.
    
  Bool_t           fUseCentralityWeight ; ///<  Return the centratlity weight.
    
  //
  // MC weights
  //
    
  Double_t         fMCWeight ;            ///<  pT-hard bin MC weight. It is used only internally.

  TString          fCurrFileName ;        ///<  Current file path name.

  Bool_t           fCheckMCCrossSection ; ///<  Retrieve from the pyxsec.root file the cross section, only if requested.
  
  Bool_t           fJustFillCrossSecHist; ///< Do not provide a weight, just fill cross section histograms
  
  TH1F *           fhXsec   ;             //!<! Cross section in PYTHIA.
   
  TH1F *           fhTrials ;             //!<! Number of event trials in PYTHIA.
    
  AliGenPythiaEventHeader * fPyEventHeader ; //!<! Pythia event header, needed to retrieve cross section, only in recent MC
  
  Bool_t           fCheckPythiaEventHeader ; ///< Get cross section from pythia event header 
  
  /// Copy constructor not implemented.
  AliAnaWeights(           const AliAnaWeights&); 
  
  /// Assignment operator not implemented.
  AliAnaWeights& operator=(const AliAnaWeights&); 
  
  /// \cond CLASSIMP
  ClassDef(AliAnaWeights, 2) ;
  /// \endcond
  
} ;

#endif //ALIANAWEIGHTS_H
