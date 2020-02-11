#ifndef ALIANALYSISTASKTPCTOFCascade_H
#define ALIANALYSISTASKTPCTOFCascade_H

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"

class AliESDEvent;
class AliMCEvent;
//class AliStack;
class AliPhysicsSelection;
class AliESDtrackCuts;
class AliESDpid;
class AliESDtrack;
class AliESDv0;
class AliESDVertex;
class AliAnalysisPIDCascadeV0;
class AliAnalysisPIDCascade;
class AliAODVertex;
class AliTOFcalib;
class AliTOFT0maker;
class TList;
class TH1F;
class TH2F;
class TH1D;
class TObjArray;
class AliAnalysisPIDCascadeEvent;
class AliAnalysisPIDCascadeTrack;
class AliAnalysisPIDCascadeParticle;
class TClonesArray;
class AliCentrality;
class AliPIDResponse;
class AliMultSelection;
class AliAnalysisUtils;
class AliESDtrackCuts;
class TTree;
class AliKFVertex;
class AliKFParticle;
class AliVVZERO;
class AliAnalysisTaskTPCTOFCascade :
public AliAnalysisTaskSE
{

 public:

  AliAnalysisTaskTPCTOFCascade(); // default constructor
  AliAnalysisTaskTPCTOFCascade(Bool_t isMC); // default constructor
  virtual ~AliAnalysisTaskTPCTOFCascade(); // default destructor

  virtual void UserCreateOutputObjects(); // user create output objects
  
  virtual void UserExec(Option_t *option); // user exec
  virtual void Terminate(Option_t *option); // terminate
  virtual Int_t FindCommonMother(Int_t label_1, Int_t label_2);


  /* getters */
  AliESDpid *GetESDpid() const {return fESDpid;}; // get ESD PID
  AliTOFcalib *GetTOFcalib() const {return fTOFcalib;}; // getter
  AliTOFT0maker *GetTOFT0maker() const {return fTOFT0maker;}; // getter
  
  /* setters */
  void SetMCFlag(Bool_t value = kTRUE) {fMCFlag = value;}; // setter
  void SetMCTuneFlag(Bool_t value = kTRUE) {fMCTuneFlag = value;}; // setter
  void SetPbPbFlag(Bool_t value = kTRUE) {fPbPbFlag = value;}; // setter
  void SetVertexSelectionFlag(Bool_t value = kTRUE) {fVertexSelectionFlag = value;}; // setter
  void SetVertexCut(Double_t value) {fVertexCut = value;}; // setter
  void SetRapidityCut(Double_t value) {fRapidityCut = value;}; // setter
  void SetTimeResolution(Double_t value) {fTimeResolution = value;}; // setter
  void ProcessV0s();
  void ProcessCascades();
  void FillHist(Double_t myflag);
  Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex = NULL);
  Bool_t SelectVertex2015pp(AliESDEvent *esd,
			    Bool_t checkSPDres = kTRUE, //enable check on vtx resolution 
			    Bool_t *SPDandTrkExists = NULL, //ask for both trk and SPD vertex
			    Bool_t *checkProximity = NULL); //apply cut on relative position of spd and trk verteces

 protected:

  AliAnalysisTaskTPCTOFCascade(const AliAnalysisTaskTPCTOFCascade &); // copy constructor
  AliAnalysisTaskTPCTOFCascade &operator=(const AliAnalysisTaskTPCTOFCascade &); // operator=
 

  /* methods */
  Bool_t InitRun(); // init run
  Bool_t InitEvent(); // init event
  Bool_t HasPrimaryDCA(AliESDtrack *track); // has primary DCA
  Bool_t MakeTPCPID(AliESDtrack *track, Double_t *nsigma, Double_t *signal); // make TPC PID
  Bool_t MakeTOFPID(AliESDtrack *track, Double_t *nsigma, Double_t *signal); // make TOF PID
  Int_t GetTrackCutsFlag(AliESDtrack *LocalTrack);
  /* flags */
  //AliESDtrackCuts *fESDtrackCuts;
  Bool_t fInitFlag; // init flag
  Bool_t fMCFlag; // MC flag
  Bool_t fMCTuneFlag; // MC tune flag
  Bool_t fPbPbFlag; // PbPb flag
  Bool_t fVertexSelectionFlag; // vertex selection flag
  Bool_t fPrimaryDCASelectionFlag; // primary DCA selection flag
  TTree *fPIDTree;
  TH1D *fEvHist;
  /* ESD analysis */
  AliPIDResponse *fPIDResponse; //! PID object
  AliAnalysisUtils *fAnUtils; //! Analysis Utils
  Int_t fRunNumber; // run number
  UInt_t fStartTime; // start time
  UInt_t fEndTime; // end time
  AliESDEvent *fESDEvent; // ESD event
  AliMCEvent *fMCEvent; // MC event
  //AliStack *fMCStack; // MC stack
  AliESDtrackCuts *fTrackCutsV0;
  AliESDtrackCuts *fTrackCuts2010; //! ITSTPC track cuts 2010
  AliESDtrackCuts *fTrackCuts2011; //! ITSTPC track cuts 2011
  AliESDtrackCuts *fTrackCutsTPCRefit; //! TPC only track cuts + refit
  AliESDtrackCuts *fTrackCuts2011Sys; //! TPC only track cuts + refit 
  AliESDpid *fESDpid; // ESD PID
  Bool_t fIsCollisionCandidate; // is collision candidate
  UInt_t fIsEventSelected; // is event selected
  Bool_t fIsPileupFromSPD; // is pile-up from SPD
  Bool_t fHasVertex; // has vertex
  Float_t fVertexZ; // vertex z
  Float_t fMCTimeZero; // MC time-zero
  AliCentrality *fCentrality; // centrality
  
  AliAnalysisPIDCascadeEvent *fAnalysisEvent; // analysis event
  TClonesArray *fAnalysisTrackArray; // analysis track array
  AliAnalysisPIDCascadeTrack *fAnalysisTrack; // analysis track
  TClonesArray *fAnalysisParticleArray; // analysis particle array
  AliAnalysisPIDCascadeParticle *fAnalysisParticle; // analysis particle
  TClonesArray *fAnalysisV0TrackArray; //V0 track array
  AliAnalysisPIDCascadeV0 *fAnalysisV0Track; //V0 track object
  TClonesArray *fAnalysisCascadeTrackArray; //Cascade track array
  AliAnalysisPIDCascade *fAnalysisCascadeTrack; //Cascade track object

  /* TOF related */
  AliTOFcalib *fTOFcalib; // TOF calib
  AliTOFT0maker *fTOFT0maker; // TOF-T0 maker
  Float_t fTimeResolution; // time resolution

  /*** CUTS ***/

  /* vertex cut */
  Double_t fVertexCut; // vertex cut
  Double_t fRapidityCut; // rapidity cut

  /*** HISTOGRAMS ***/
  
  TH1F *V0MBinCount;

  /* histo lists */
  TList *fHistoList; // histo list
  TList *fMCHistoList; // MC histo list

  
  ClassDef(AliAnalysisTaskTPCTOFCascade, 3);
};

#endif /* ALIANALYSISTASKTPCTOFCascade_H */
