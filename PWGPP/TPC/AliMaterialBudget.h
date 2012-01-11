#ifndef ALIMATERIALBUDGET_H
#define ALIMATERIALBUDGET_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliTPCseed.h>
#include <TString.h>
class AliGenInfoMaker;
class TTreeSRedirector;
class AliMCEventHadnler;
class TParticle;
class AliMCInfo;
class AliMCParticle;
class AliESDRecInfo;
class AliESDEvent;
class AliMCEvent;
class AliComparisonObject;

class AliMaterialBudget : public AliAnalysisTask {
 public:
 AliMaterialBudget();
 AliMaterialBudget(const char *name);
  virtual ~AliMaterialBudget();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void FinishTaskOutput();
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}

  //
  void           FindPairs(AliESDEvent * event);
  Bool_t         IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1);
  //
  void           ProcessMCInfo();
  void           ProcessRefTracker(AliTrackReference* refIn, AliTrackReference* refOut, TParticle*part, Int_t type);
  
  void           FitTrackRefs(TParticle * part, TClonesArray * trefs);

  //
  // debug streamer part
  //
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  //
  static Bool_t PropagateCosmicToDCA(AliExternalTrackParam *param0, AliExternalTrackParam *param1, Double_t mass);
  static AliExternalTrackParam * MakeTrack(const AliTrackReference* ref, TParticle*part);
  static Bool_t  PropagateToPoint(AliExternalTrackParam *param, Double_t *xyz, Double_t mass,  Float_t step);
  //
  AliTrackReference * GetFirstTPCTrackRef(AliMCParticle *mcParticle);
  AliTrackReference * GetAllTOFinfo(AliMCParticle *mcParticle, Int_t & nTrackRef, Int_t &nTrackRefITS, Int_t retValue =0);
 protected:
  void RegisterDebugOutput();
  AliMaterialBudget(const AliMaterialBudget& /*info*/);
  AliMaterialBudget& operator=(const AliMaterialBudget& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;          //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  //
  //
  //
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level 
  Int_t  fDebugLevel;                   //  debug level
  TString      fDebugOutputPath; // debug output path
  //
  // histogran
  //
  TList * fListHist;     // list for histograms
  TH1F  * fHistMult;     // track multiplicity histograms
  //
  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutMaxDz;    // maximal distance in z ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products
  //
  ClassDef(AliMaterialBudget, 1); // Analysis task base class for tracks
};

#endif
