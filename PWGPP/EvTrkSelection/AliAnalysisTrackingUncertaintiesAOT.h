#ifndef ALIANALYSISTRACKINGUNCERTAINTIESAOT_H
#define ALIANALYSISTRACKINGUNCERTAINTIESAOT_H

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
//class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "THn.h"
#include <THnSparse.h>
#include <Rtypes.h>

class AliAnalysisTrackingUncertaintiesAOT : public AliAnalysisTaskSE {
    
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
  enum ECentrality {kCentOff,kCentV0M,kCentCL1,kCentZNA,kCentV0A,kCentInvalid};

  // list of possible standard ESD track cuts set
  enum EtrkCuts {
    kDefault=0,
    kStdTPConlyTrkCuts=1,
    kStdITSTPCTrkCuts2009,
    kStdITSTPCTrkCuts2010,
    kStdITSTPCTrkCuts2011,
    kStdITSTPCTrkCuts2015PbPb
    // to be implemented, if needed
    //kStdITSSATrkCuts2009,
    //kStdITSSATrkCuts2010,
    //kStdITSSATrkCutsPbPb2010,
    //kStdITSPureSATrackCuts2009,
    //kStdITSPureSATrackCuts2010
  };

  AliAnalysisTrackingUncertaintiesAOT(const char *name);
  AliAnalysisTrackingUncertaintiesAOT();
  virtual ~AliAnalysisTrackingUncertaintiesAOT();
    
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
  void           ProcessTracks(AliMCEvent *mcEvent);
    
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;}
  void           InitializeTrackCutHistograms();
  void           SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void           SetMaxDCAxy(Double_t maxDCA) {fMaxDCAxy = maxDCA;}
  void           SetMaxDCAz(Double_t maxDCA)  {fMaxDCAz  = maxDCA;}
  void           SetEtaRange(Double_t maxEta) {fMaxEta   = maxEta;}
  void           SetCrossRowsOverFndCltTPC(Double_t CrossRowsOverFndClt) {fCrossRowsOverFndCltTPC = CrossRowsOverFndClt;}
  void           SetTriggerClass(TString trigClass)   {fTriggerClass = trigClass;}
  void           SetTriggerMask(ULong64_t mask=0)     {fTriggerMask  = mask;}
  void           SetSpecie(ULong64_t specie=0)        {fspecie = specie;}
  void           SetRequireTrackVtx(Bool_t flag)      {fRequireVtxTracks = flag;}
  void           SetUseGeneratedPt(Bool_t flag)       {fUseGenPt = flag;}
  void           SetUsePtLogScale(Bool_t flag)        {fUsePtLogAxis = flag;}
  void           SetUseFinePtAxis(Bool_t flag)        {fUseFinePtAxis = flag;}
  void           SetUseCutV0multVsTPCout(Bool_t flag) {fDoCutV0multTPCout=flag;}
  void           SetSPDRequirement(AliESDtrackCuts::ITSClusterRequirement  spdlayreq)   {fSPDlayerReq = spdlayreq;}                   
  void           SetMultSelectionObjectName(TString str){fMultSelectionObjectName=str;}
  void           SetMinCentrality(Double_t val) {fminCent = val;}
  void           SetMaxCentrality(Double_t val) {fmaxCent = val;}
  void           SetUseCentrality(AliAnalysisTrackingUncertaintiesAOT::ECentrality flag);
  void           SetDCAzOn(Bool_t flag = kTRUE) {fDCAz = flag;}
  void           SetTPConly(Bool_t tpconly = kTRUE) {fTPConlyFIT = tpconly;}

  // make the pT binning finer by a factor of 2
  void           SetFinerpTbin(Bool_t flag) {fmakefinerpTbin=flag;}
  
  // geometrical cut for tracks in the TPC (copied from AliRDHFCuts, 2019-May-23rd)
  void SetUseCutGeoNcrNcl(Bool_t opt){fUseCutGeoNcrNcl=opt;}
  void ConfigureCutGeoNcrNcl(Double_t dz, Double_t len, Double_t onept, Double_t fncr, Double_t fncl){
    fDeadZoneWidth=dz;  fCutGeoNcrNclLength=len; fCutGeoNcrNclGeom1Pt=onept;
    fCutGeoNcrNclFractionNcr=fncr; fCutGeoNcrNclFractionNcl=fncl;
  }

  // possibility to modify the ESD track cuts set
  void SetStandardESDtrkCuts(UInt_t whichcuts, UInt_t option_TPCclstcut)  {fWhichCuts=whichcuts; fTPCclstCut=option_TPCclstcut;}

  ULong64_t GetTriggerMask() {return fTriggerMask;}
  ULong64_t GetSpecie() {return fspecie;}

 private:
    
  void   BinLogAxis(const THnSparseF *h, Int_t axisNumber);
  void   BinFinePt(const THnSparseF *h, Int_t axisNumber);
  Bool_t IsVertexAccepted(AliESDEvent * esd);
  Bool_t IsElectron(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsPion(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsKaon(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsProton(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsConsistentWithPid(Int_t specie, const AliESDtrack * const tr);//Int_t type, const AliESDtrack * const tr);
  Bool_t IsEventSelectedInCentrality(AliESDEvent *ESDevent);
    
  AliAnalysisTrackingUncertaintiesAOT::ECentrality   fUseCentrality;
  Double_t fMaxDCAxy;
  Double_t fMaxDCAz;
  Double_t fMaxEta;
  Double_t fCrossRowsOverFndCltTPC;
  Double_t fminCent;
  Double_t fmaxCent;
  AliESDtrackCuts::ITSClusterRequirement fSPDlayerReq; // SPD layers requirement 
    
  TString  fTriggerClass;           /// trigger class
  ULong64_t fTriggerMask;           /// trigger mask
    
  AliESDEvent * fESD;               //! ESD object
  AliESDpid   * fESDpid;            //! basic pid object
  ULong64_t   fspecie;
    
  TH1F *fHistNEvents;               //! histo with number of events
  TH1F *fHistCent;                  //! histo for the centrality 
  THnSparse *fHistMC;               //! sparse of the tracks on MC and ITS-TOC matching
  THnSparse *fHistMCTPConly;        //! sparse of the tracks on MC and only TPC request
  THnSparse *fHistData;             //! sparse of the tracks on data and ITS-TPC matching
  TH2F *fHistAllV0multNTPCout;      //! histo for V0mult vs #tracks TPCout (all)
  TH2F *fHistSelV0multNTPCout;      //! histo for V0mult vs #tracks TPCout (sel)

  Bool_t   fMC;                     //flag to switch on the MC analysis for the efficiency estimation
  Bool_t   fRequireVtxTracks;       //flag to require track vertex, if false accepts also SPD
  Bool_t   fUsePtLogAxis;           //flag to use log scale on pt axis in match. eff. sparse
  Bool_t   fUseFinePtAxis;          //flag to use fine bin width for low pt axis in match. eff. sparse
  Bool_t   fUseGenPt;               //flag to use generated pt in match. eff. sparse
  Bool_t   fDoCutV0multTPCout;      //flag to activate cut on V0mult vs #tracks TPCout
  Bool_t   fDCAz;                   //flag to switch on the DCAz axis
  Bool_t   fTPConlyFIT;             //flag to use only TPC track for DCA fits

  TString fMultSelectionObjectName; /// name of the AliMultSelection object to be considered

  TList           * fListHist;      //! output list for histograms
  AliESDtrackCuts * fESDtrackCuts;  //! cut set which is under study
  AliESDVertex    * fVertex;        //! pointer to ESD vertex
    
  // make the pT binning finer by a factor of 2
  Bool_t fmakefinerpTbin;

  // parameters for geometrical cut for tracks in the TPC (copied from AliRDHFCuts, 2019-May-23rd)
  Bool_t fUseCutGeoNcrNcl; /// flag for enabling/disabling geometrical cut on TPC track
  Double_t fDeadZoneWidth;       /// 1st parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclLength;  /// 2nd parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclGeom1Pt; /// 3rd parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclFractionNcr; /// 4th parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclFractionNcl; /// 5th parameter of GeoNcrNcl cut

  // possibility to modify the ESD track cuts set
  UInt_t fWhichCuts;  ///
  UInt_t fTPCclstCut; /// 0: cut on TPC clusters; 1: cuts on the number of crossed rows and on the ration crossed rows/findable clusters

  AliAnalysisTrackingUncertaintiesAOT(const AliAnalysisTrackingUncertaintiesAOT&);
  AliAnalysisTrackingUncertaintiesAOT& operator=(const AliAnalysisTrackingUncertaintiesAOT&);
    
  ClassDef(AliAnalysisTrackingUncertaintiesAOT, 11);
};

#endif

