/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliAnalysisTaskPhiSA.h*/

#ifndef ALIANALYSISTASKPHISA_H
#define ALIANALYSISTASKPHISA_H


class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TProfile2D;
class THnSparse;
class TList;
class TString;
class TVector2;

class AliQnCorrectionsManager;
class AliQnCorrectionsCutsSet;
class AliQnCorrectionsHistos;
class AliPIDResponse;
class AliTriggerAnalysis;
class AliESDEvent;
class AliEventplane;
class AliFlowTrackCuts;
class AliESDtrackCuts;
class AliESDtrack;
class AliMultiplicity; 
class AliVTrack;
class AliMultSelection;    



//#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsFillEventTask.h"
#include "AliQnCorrectionsQnVector.h"
#include "TMCProcess.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliESDpid.h"

//#endif

class AliAnalysisTaskPhiSA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPhiSA();
  AliAnalysisTaskPhiSA(const char *name,  const Bool_t useshift);
  virtual ~AliAnalysisTaskPhiSA();
  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  
  //Some seters 
  
  void SetAnalysisLevel(Char_t *level){fAnalysisLevel = level;};
  void SetPOIAndRPTrackType(TString tracktype,TString rptracktype);
  void CheckPileUp(Bool_t pile){fCheckPileUp = pile;};
  void AcceptTPCVertex(Bool_t vertpc){fAcceptTPC = vertpc;};    
  void SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  void SetIgnoreTPCzRange( Double_t min, Double_t max ){ fIgnoreTPCzRange=kTRUE; fIgnoreTPCzRangeMin=min; fIgnoreTPCzRangeMax=max; }
  void SetMinimalTPCdedx(Double_t d=10.)        {fMinimalTPCdedx=d; fCutMinimalTPCdedx=kTRUE;}

  void SetPairRapidityCut(Double_t y)           {fPairRapidity = y;}
  void SetDipAngle(Double_t dipangle)           {fDeltaAngle = dipangle; }
  void SetSystemaicCutType(TString systcut)           {fSytCutType = systcut.Data(); fSystematicCuts = kTRUE;}

  //Shifting
  void    SetShiftList(TList* const aShiftList)  {this->fShiftList = (TList*)aShiftList->Clone();}
  TList*  GetShiftList() const                     {return this->fShiftList;}  
  void    SetUseShifting(Bool_t const aShift)       {this->fUseShift = aShift;}
  Bool_t  GetUseShifting() const                   {return this->fUseShift;}


  //getters
  Double_t GetMinimalTPCdedx() const {return fMinimalTPCdedx;}
  const   char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  
 private:
    enum {
      kDim    = 50000,//Maximum number of tracks.
      kDimBuf = 500,//NmixEvent * (NveterxBin + NCentralitiyBin) = 5 * (10+90)
      kCenBin=9,
    };

    //void InitESDcuts() {if (!fAliESDtrackCuts) {fAliESDtrackCuts=new AliESDtrackCuts();}}
    
    TList              *fOutput;                  // Output list
    
    AliESDtrackCuts    *fAliESDtrackCuts;        //alianalysis cuts
    AliESDtrackCuts    *fESDtrackCutsForRP;      //alianalysis cuts for RP
    AliPIDResponse     *fPIDResponse;             // PID response Handler
    AliMultSelection      *fMultSelection;      //! For Centrality 

    AliQnCorrectionsManager *fFlowQnVectorMgr;
    AliTriggerAnalysis *fTriggerAna;              // trigger analysis
    TString  fAnalysisLevel;                      // "ESD"
    TString  fTrackType;                          // "GLOBAL", "TPC"
    TString  fRPTrackType;                          // "GLOBAL", "TPC"
    TString  fSytCutType;
    Bool_t   fSystematicCuts;
    TList    *fShiftList;   // list holding input profile with average sin and cos in shifting
    Bool_t   fUseShift;    // use shift correction
    Bool_t   fAcceptTPC;                          // if kTRUE, the TPC primary vertexes are accepted
    Bool_t   fCheckPileUp;                        // if kTRUE, the SPD pile up is checked
    Bool_t   fIgnoreTPCzRange;                     //ignore tracks going close to central membrane
    Double_t fIgnoreTPCzRangeMax;                  //max z to ignore
    Double_t fIgnoreTPCzRangeMin;                  //min z to ignore
    Bool_t   fCutMinimalTPCdedx;                   //cut on minimal dedx in TPC to reject noise tracks
    Double_t fMinimalTPCdedx;                      //value for minimal TPC dedx
    Bool_t   fUsercuts;	                           // bool, if personal cuts should be used

    Double_t kaonMass;                             //   = 0.49368
    Double_t pionMass;                             //   = 0.49368
    Double_t phiMass;                              //    = 1.01946
    Double_t kstarMass;                            //    = 1.01946

    Int_t    fCurrentEventCentralityBin;//!
    Int_t    fCurrentEventDoughterOne;//!
    Int_t    fCurrentEventDoughterTwo;//!
    Float_t  fCurrentEventVx;//!
    Float_t  fCurrentEventVy;//!
    Float_t  fCurrentEventVz;//!
    Int_t    fBufferPointer;//!
    //Int_t    fRunNumber;//! runnumber

    TVector2* fQVector;	//! Q-Vector of the event
    TVector2* fQVZeroA;//! Q-Vector of VZero A
    TVector2* fQVZeroC;//! Q-Vector of VZero C

    Double_t fPairRapidity; //
    Double_t fDeltaAngle;//

    
    Int_t    fCurrentEventDoughterOneCharge[kDim];//!
    Int_t    fCurrentEventKId[kDim];//!
    Int_t    fCurrentEventDoughterOneIn[kDim];//!
    Int_t    fCurrentEventDoughterOneId[kDim];//!
    Int_t    fCurrentEventDoughterOneTrackNumber[kDim];//!
    Double_t fCurrentEventDoughterOnePx[kDim];//!
    Double_t fCurrentEventDoughterOnePy[kDim];//!
    Double_t fCurrentEventDoughterOnePz[kDim];//!
    Double_t fCurrentEventDoughterOneQx[kDim];//!
    Double_t fCurrentEventDoughterOneQy[kDim];//!

    Int_t    fCurrentEventDoughterTwoCharge[kDim];//!
    Int_t    fCurrentEventDoughterTwoIn[kDim];//!
    Int_t    fCurrentEventDoughterTwoId[kDim];//!
    Int_t    fCurrentEventDoughterTwoTrackNumber[kDim];//!
    Double_t fCurrentEventDoughterTwoPx[kDim];//!
    Double_t fCurrentEventDoughterTwoPy[kDim];//!
    Double_t fCurrentEventDoughterTwoPz[kDim];//!
    Double_t fCurrentEventDoughterTwoQx[kDim];//!
    Double_t fCurrentEventDoughterTwoQy[kDim];//!
  
    Int_t    fBufferEventNEvents[100];//!
    Int_t    fBufferEventFull[100];//!
    Double_t fBufferEventPsi[kDimBuf];//!
    Int_t    fBufferEventDoughterOne[kDimBuf];//!
    Int_t    fBufferEventDoughterOneCharge[kDimBuf][kDim];//!
    Int_t    fBufferEventDoughterOneId[kDimBuf][kDim];//!
    Int_t    fBufferEventDoughterOneIn[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterOnePx[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterOnePy[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterOnePz[kDimBuf][kDim];//!
    Int_t    fBufferEventDoughterTwo[kDimBuf];//!
    Int_t    fBufferEventDoughterTwoCharge[kDimBuf][kDim];//!
    Int_t    fBufferEventDoughterTwoId[kDimBuf][kDim];//!
    Int_t    fBufferEventDoughterTwoIn[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterTwoPx[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterTwoPy[kDimBuf][kDim];//!
    Double_t fBufferEventDoughterTwoPz[kDimBuf][kDim];//!


    // NEW HISTO to be declared here
    TH1F  *fHistZVertex;//!        
    TH1F  *fHistCentralityEvtCount;//!
    TH1F  *fHistEventCount;//!
    TH1F  *fHistCentDist;//!

    
    TH2F  *fHistPionVsKaonMult;//!
    
    /*TH1F  *fHistTPCClusterRP;//!
    TH1F  *fHistITSClusterRP;//!
    TH1F  *fHistChiSqrPerNdfTPCRP;//!
    TH1F  *fHistChiSqrPerNdfITSRP;//!

    TH2F  *fHistdEdxKaonTpc;//!
    TH2F  *fHistdEdxKaonTof;//!
    TH2F  *fHistdEdxPionTpc;//!
    TH2F  *fHistdEdxPionTof;//!
    TH1F  *fHistTPCClusterPOI;//!
    TH1F  *fHistITSClusterPOI;//!
    TH1F  *fHistChiSqrPerNdfTPCPOI;//!
    TH1F  *fHistChiSqrPerNdfITSPOI;//!
    TH1F  *fHistnCrossRowsPOI;//!
    TH1F  *fHistratioCrossedRowsOverFindableClustersTPC;//!
    TH1F  *fHistDCAxyPOI;//!
    TH1F  *fHistDCAzPOI;//!*/
    TProfile *fProCosResThreeSubEventPsiAB;//!
    TProfile *fProCosResThreeSubEventPsiAC;//!
    TProfile *fProCosResThreeSubEventPsiBC;//!
    /*
    TH3F *fhInvMassSAEPvzeroA[9];//!
    TH3F *fhInvMassLikePPSAEPvzeroA[9];//!
    TH3F *fhInvMassLikeMMSAEPvzeroA[9];//!
    TH3F *fhInvMassMixSAEPvzeroA[9];//!
    TH3F *fhInvMassSAEPvzeroC[9];//!
    TH3F *fhInvMassLikePPSAEPvzeroC[9];//!
    TH3F *fhInvMassLikeMMSAEPvzeroC[9];//!
    TH3F *fhInvMassMixSAEPvzeroC[9];//!
    */


    TH1F  *fHistEPTPC[9];//!
    TH1F  *fHistEPV0A[9];//!
    TH1F  *fHistEPV0C[9];//!
    TH1F  *fHistEPFMDA[9];//!
    TH1F  *fHistEPFMDC[9];//!
    TH1F  *fHistEPZDCA[9];//!
    TH1F  *fHistEPZDCC[9];//!
    TH1F  *fHistEPT0A[9];//!
    TH1F  *fHistEPT0C[9];//!

    TH1F  *fHistEPV0A_rw[9];//!
    TH1F  *fHistEPV0C_rw[9];//!
    TH1F  *fHistEPV0A_plain[9];//!
    TH1F  *fHistEPV0C_plain[9];//!
    TH1F  *fHistEPV0A_rec[9];//!
    TH1F  *fHistEPV0C_rec[9];//!
    TH1F  *fHistEPV0A_recShift[9];//!
    TH1F  *fHistEPV0C_recShift[9];//!
    TH1F  *fHistEPV0A_al[9];//!
    TH1F  *fHistEPV0C_al[9];//!

    TH1F  *fHistEPV0A_tw[9];//!
    TH1F  *fHistEPV0C_tw[9];//!
    TH1F  *fHistEPV0A_sc[9];//!
    TH1F  *fHistEPV0C_sc[9];//!



    // Shift correction
    TProfile2D           *fShiftCosTerm_V0A;   // profile holding av. cosine terms for full and sub events
    TProfile2D           *fShiftSinTerm_V0A;   // profile holding av. sine terms for full and sub events
    TProfile2D           *fShiftCosTerm_V0C;   // profile holding av. cosine terms for full and sub events
    TProfile2D           *fShiftSinTerm_V0C;   // profile holding av. sine terms for full and sub events
    
    // profile for shifting
    TProfile2D *fFullCosTerm_V0A;  
    TProfile2D *fFullSinTerm_V0A;
    TProfile2D *fFullCosTerm_V0C;  
    TProfile2D *fFullSinTerm_V0C;

    THnSparse *fhInvMassSAEPvzeroA;//!
    THnSparse *fhInvMassLikePPSAEPvzeroA;//!
    THnSparse *fhInvMassLikeMMSAEPvzeroA;//!
    THnSparse *fhInvMassMixSAEPvzeroA;//!
    THnSparse *fhInvMassSAEPvzeroC;//!
    THnSparse *fhInvMassLikePPSAEPvzeroC;//!
    THnSparse *fhInvMassLikeMMSAEPvzeroC;//!
    THnSparse *fhInvMassMixSAEPvzeroC;//!
    

    Bool_t   PassEvent(AliVEvent *event);
    Int_t    GetCentrality(AliVEvent *evt);
    Float_t  GetCentralityValue(AliVEvent *evt);
    Bool_t   PassTrack(AliVTrack *track);
    Int_t    MakeRealPair(TVector2 *qvA, TVector2 *qvC);
    Int_t    MakeLikePair(TVector2 *qvA, TVector2 *qvC);
    Int_t    MakeLikeSignPair(TVector2 *qv);
    Int_t    MakeRotationalBkgPair(TVector2 *qv);
    Int_t    MakeMixedPair(Int_t bufferPointer,TVector2 *qvA, TVector2 *qvC);
    void     CopyCurrentToBuffer(Int_t bufferPointer,TVector2 *qv);
    Double_t GetDipAngle(Double_t aX,Double_t aY,Double_t aZ,
			 Double_t bX,Double_t bY,Double_t bZ);
    Double_t CosThetaStar(TLorentzVector mother, TLorentzVector daughter0,TVector2& Qvect);
    Double_t CosPhiPsi(TLorentzVector mother, TLorentzVector daughter0,TVector2& Qvect);

    Bool_t   CheckVertex(const AliVVertex *vert);
    Double_t GetNSigmaCut(Double_t pTrack);
    Double_t GetWeight(TObject* track1,TString flag);
    Double_t GetPhiWeight(TObject* track1, TString flag);
    
    Bool_t   IsSelectedPION(AliVTrack *track);//subhash
    Bool_t   IsSelectedKAON(AliVTrack *track);//subhash
    Bool_t   GetStatus(const AliVTrack *track);
    Bool_t   MatchTOF(const AliVTrack *vtrack);
    Double_t GetTOFBeta( AliVTrack *vtrack);
    void     ReSet();   

    AliAnalysisTaskPhiSA(const AliAnalysisTaskPhiSA&); // not implemented
    AliAnalysisTaskPhiSA& operator=(const AliAnalysisTaskPhiSA&); // not implemented
    ClassDef(AliAnalysisTaskPhiSA, 1); // example of analysis
};
//__________________________________________________________________________________________________
inline Bool_t AliAnalysisTaskPhiSA::CheckVertex(const AliVVertex *vertex)
{
  //
  // Checks if a candidate primary vertex is good,
  // which is true if it is not null and has at
  // least one contributor
  //
  if (!vertex) return kFALSE; 
  if (vertex->GetNContributors() < 1) return kFALSE;
  return kTRUE;
}
#endif

