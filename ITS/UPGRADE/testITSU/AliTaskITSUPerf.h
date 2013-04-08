#ifndef ALIPERFITSUTASK_H
#define ALIPERFITSUTASK_H

///////////////////////////////////////////////////////////////////////////
// Class AliTaskITSUPerf                                                 //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////

class TH1F; 
class TH2F;
class TH3F;
class AliESDEvent;
class TList;
class TNtuple;

class AliMCParticle;
class AliITSMultRecBg;
class AliESDTrackCuts;
class AliITSUGeomTGeo;
class AliITSURecoDet;

#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include <TMath.h>

class AliTaskITSUPerf : public AliAnalysisTaskSE {
 public:
  enum {kMCPrimBit=31     // flag for primaries
	,kTrCondFail=30   // does not correspond to any tracking condition
	,kITSHitBits=0  // ITS hits pattern stored starting from this bit
  };
  //
  enum {kITS1TPC1,kITS1TPC0,kITS0TPC1,kITS0TPC0,kITSTPCMismatch,kITSTPCNoMatch,kNLabelTypes};
  //
  enum {  // standard histo ID's defined for each centrality and MClabels combination bin
    kHResPTvsPTMC        // pt resolution
    ,kHResDCARvsPTMC     // DCA R resolution
    ,kHResDCAZvsPTMC     // DCA Z resolution
    ,kHNStdHistosCentMC
    //
  };
  //
  enum {  // histos defined for each centrality bin (regardless MClabels combination)
    kHMatchStatus
    ,kHNStdHistosCent
  };
    
  //
  AliTaskITSUPerf(const char *name = "AliTaskITSUPerf");
  virtual ~AliTaskITSUPerf(); 
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  //
  Bool_t     GetUseSpecialOutput()            const {return fUseSpecialOutput;}
  void       SetUseSpecialOutput(Bool_t v=kTRUE)    {fUseSpecialOutput=v;}
  //
  void       CheckTracks();
  void       BuildMCInfo();
  Int_t      GetMCLabType(Int_t labMCTPC,Int_t labMCITS, Int_t nClTPC, Int_t nClITS);
  void       SetTrackingConditions(const TObjArray* arr) {fTrackingCond = arr;}
  //
  //-------------------------------------------------------------
  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}
  void       BookHistos(Int_t bin);
  void       AddHisto(TObjArray* array, TObject* h, Int_t at=-1);
  void       BookStandardHistosCentMCLb(Int_t bin, Int_t mcLb);
  void       BookStandardHistosCent(Int_t bin);
  Int_t      GetHistoID(Int_t htype, Int_t mcStat=-1, Int_t centBin=0) const;
  TH1*       GetHisto(const TObjArray* array, Int_t htype, Int_t mcStat=-1, Int_t centBin=0)  const {return (TH1*)array->At(GetHistoID(htype,mcStat,centBin));}
  //
  Int_t      GetCentralityBin() const;
  //
  void       SetEtaCut(Float_t etaCut)          {fEtaMax = TMath::Abs(etaCut); fEtaMin= -fEtaMax;}
  void       SetEtaMin(Float_t etaMin)          {fEtaMin = etaMin;}
  void       SetEtaMax(Float_t etaMax)          {fEtaMax = etaMax;}
  void       SetPtMin(Float_t  ptMin)           {fPtMin = ptMin;}
  void       SetPtMax(Float_t  ptMax)           {fPtMax = ptMax;}
  void       SetNPtBins(Int_t n=20)             {fNPtBins = n;}
  void       SetZVertexMin(Float_t z)           {fZVertexMin = z;}
  void       SetZVertexMax(Float_t z)           {fZVertexMax = z;}
  void       SetNResBins(Int_t n=100)           {fNResBins = n;}
  //
 protected:
  //
 protected:
  TList*           fOutput;                  // output list send on output slot 1 
  TObjArray        fHistosCentMCLb;          //! local array for histos management, centrality&MClabel selective
  TObjArray        fHistosCent;              //! local array for histos management, centrality selective
  TTree*           fRPTree;                  //! tree of recpoints
  AliStack*        fStack;                   //! MC stack
  AliMCEvent*      fMCEvent;                 //! MC Event
  Float_t          fVtxMC[3];                //! MC gen vertex
  const AliESDVertex* fVtxSPD;               //! SPD vertex
  const AliESDVertex* fVtxTrc;               //! Tracks vertex
  AliESDEvent*     fESDEvent;                //! ESDEvent
  Bool_t           fUseSpecialOutput;        // flag to open special output
  Bool_t           fUseMC;                   // do we use MC info?
  const TObjArray* fTrackingCond;            //! tracking conditions used
  //
  //---------------------------------------------------
  //
  AliITSUGeomTGeo* fGeom;                    //! general interface to ITS geometry
  AliITSURecoDet*  fITS;                     //! interface to ITS reco time (for access to clusters, geometry)
  //
  // MC info
  Int_t            fNSelTracksMC;            //! number of selected MC tracks
  TArrayI          fMCStatus;                //! mc info about every MC particle
  // 
  Int_t            fNPtBins;                 //! N pt bins for histos
  Int_t            fNResBins;                //! N bins for resolution histos
  //
  Double_t         fPtMin;                   //! min pt for histos
  Double_t         fPtMax;                   //! max pt for histos
  Double_t         fEtaMin;                  //! min eta range
  Double_t         fEtaMax;                  //! max eta range
  Double_t         fZVertexMin;              //! min Z vtx to process
  Double_t         fZVertexMax;              //! max Z vtx to process
  //
  Int_t  fCurrCentBin;                     // current centrality bin
  Int_t  fNCentBins;                       // N of mult bins
  Int_t  fUseCentralityVar;                // what is used to determine the centrality
  //
  static const char* fgkLabelTypes[kNLabelTypes]; // label truthness names
 private:    
  AliTaskITSUPerf(const AliTaskITSUPerf&); // not implemented
  AliTaskITSUPerf& operator=(const AliTaskITSUPerf&); // not implemented 
  //  
  ClassDef(AliTaskITSUPerf, 1);  
};


#endif

