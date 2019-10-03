#ifndef ALIANALYSISTRACKINGUNCERTAINTIESHMTF_H
#define ALIANALYSISTRACKINGUNCERTAINTIESHMTF_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Analysis task for the systematic study of the uncertainties related to   //
// the tracking and ITS-TPC matching efficiency for different particle      //
// species.                                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TList;
class AliESDEvent;
class AliMCEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpid;


#include "AliAnalysisTaskSE.h"
#include "THn.h"
#include <THnSparse.h>
#include <Rtypes.h>

class AliAnalysisTrackingUncertaintiesHMTF : public AliAnalysisTaskSE {
    
 public:
    
  enum {
    kNumberOfAxes = 8
  };
  enum ESpecies_t {
    kSpecElectron = BIT(0),
    kSpecPion     = BIT(1),
    kSpecKaon     = BIT(2),
    kSpecProton   = BIT(3),
    kAll          = BIT(4)
  };
    
    
  AliAnalysisTrackingUncertaintiesHMTF(const char *name);
  AliAnalysisTrackingUncertaintiesHMTF();
  virtual ~AliAnalysisTrackingUncertaintiesHMTF();
    
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
  void           ProcessTracks(AliStack *stack);
    
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;}
  void           InitializeTrackCutHistograms();
  void           SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void           SetMaxDCAxy(Double_t maxDCA) {fMaxDCAxy = maxDCA;}
  void           SetMaxDCAz(Double_t maxDCA)  {fMaxDCAz  = maxDCA;}
  void           SetEtaRange(Double_t maxEta) {fMaxEta   = maxEta;}
  void           SetCrossRowsOverFndCltTPC(Double_t CrossRowsOverFndClt) {fCrossRowsOverFndCltTPC = CrossRowsOverFndClt;}
  void           SetTriggerClass(TString trigClass) {fTriggerClass = trigClass;}
  void           SetTriggerMask(ULong64_t mask=0)   {fTriggerMask  = mask;}
  void           SetSpecie(ULong64_t specie=0)      {fspecie = specie;}
  void           SetRequireTrackVtx(Bool_t flag)    {fRequireVtxTracks = flag;}
  void           SetUsePtLogScale(Bool_t flag)      {fUsePtLogAxis = flag;}
    
  ULong64_t GetTriggerMask() {return fTriggerMask;}
  ULong64_t GetSpecie() {return fspecie;}
 private:
    
  void   BinLogAxis(const THnSparseF *h, Int_t axisNumber);
  Bool_t IsVertexAccepted(AliESDEvent * esd);
  Bool_t IsElectron(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsPion(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsKaon(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsProton(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsConsistentWithPid(Int_t specie, const AliESDtrack * const tr);//Int_t type, const AliESDtrack * const tr);
    
  Double_t fMaxDCAxy;
  Double_t fMaxDCAz;
  Double_t fMaxEta;
  Double_t fCrossRowsOverFndCltTPC;
    
  TString  fTriggerClass;           /// trigger class
  ULong64_t fTriggerMask;           /// trigger mask
    
  AliESDEvent * fESD;               //! ESD object
  AliESDpid   * fESDpid;            //! basic pid object
  ULong64_t   fspecie;
    
  TH1F *fHistNEvents;               //! histo with number of events
  THnSparse *fHistMC;               //! sparse of the tracks on MC and ITS-TOC matching
  THnSparse *fHistMCTPConly;        //! sparse of the tracks on MC and only TPC request
  THnSparse *fHistData;             //! sparse of the tracks on data and ITS-TPC matching
    
  // histograms for testing
  TH2F *fHistTrackletsTRDvsSPD;     //! histo tracklets: TRD vs SPD (for pile-up checks)

  Bool_t   fMC;                     //flag to switch on the MC analysis for the efficiency estimation
  Bool_t   fRequireVtxTracks;       //flag to require track vertex, if false accepts also SPD
  Bool_t   fUsePtLogAxis;           //flage to use log scale on pt axis in match. eff. sparse
    
  TList           * fListHist;      //! output list for histograms
  AliESDtrackCuts * fESDtrackCuts;  //! cut set which is under study
  AliESDVertex    * fVertex;        //! pointer to ESD vertex
    
  AliAnalysisTrackingUncertaintiesHMTF(const AliAnalysisTrackingUncertaintiesHMTF&);
  AliAnalysisTrackingUncertaintiesHMTF& operator=(const AliAnalysisTrackingUncertaintiesHMTF&);
    
  ClassDef(AliAnalysisTrackingUncertaintiesHMTF, 2);
};

#endif
