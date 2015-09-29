//d: AliJEfficiencyScanner.h,v 1.5 2012/04/19 15:19:52 jkral Exp $

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIJEFFICIENCYSCANNER_H
#define ALIJEFFICIENCYSCANNER_H

#include <TNamed.h>
#include "AliJRunHeader.h"
#include <iostream>
#include <TClonesArray.h>

#include <AliJConst.h>
#include <TVectorT.h>
#include "AliJMCTrack.h"
#include "AliJTrack.h"
#include <AliJTrackCut.h>
#include <TFile.h>
#include <TF1.h>

//==============================================================

#ifndef AliJMaxDimBuffer
#define AliJMaxDimBuffer
const int kMaxDimBuffer = 300;//max length of a line read to a buffe
#endif

class AliJEventHeader;
//class AliJRunHeader;
class AliJTrack;
class AliAnalysisTaskSE;

class TH1D;
class TH2D;

class AliJEfficiencyScanner : public TNamed  {

 public:
     enum { kJPhysicsPrimary, kJFake, kNPrimaryStatus };
     enum { kJMCTrack, kJGlobal, kJTPCOnly, kJGCG , kNTrackType };
     enum { kNVtxBin=1 };
     enum { kNCentBin=21};
  AliJEfficiencyScanner();
  AliJEfficiencyScanner(const char *name);
  AliJEfficiencyScanner(const AliJEfficiencyScanner& ap);   
  AliJEfficiencyScanner& operator = (const AliJEfficiencyScanner& ap);
  virtual ~AliJEfficiencyScanner();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");

  void SetJTrackList( TClonesArray * l ) { fTrackList  = l ; } 
  void SetJMCTrackList( TClonesArray * l ) { fMCTrackList  = l ; } 
  void SetJEventHeader( AliJEventHeader * h ){ fEventHeader = h ; }
  void SetJRunHeader( AliJRunHeader *h ){ fRunHeader = h; }
  void SetMBTriggMask(int mask ){ fMBTriggMask = mask; }
  
  void SetIsIsolated( Bool_t isIsolatedIn ) { fisIsolated = isIsolatedIn; }
  void SetIsRelative( Bool_t isRelaticeIn ) { fisRelative = isRelaticeIn; }
  void SetIsolParameter( double isolParamIn ) { fisolParam = isolParamIn; }
  void SetIsolCone( double isolConeIn ) { fisolCone = isolConeIn; }
  
  bool IsSelected( AliJTrack * track, int itrigger )const { return AliJTrackCut::GetInstance().IsSelected(track, itrigger); }

  AliJTrack * GetJTrack( int i ){ return (AliJTrack*) fTrackList->At(i); }
  AliJMCTrack * GetJMCTrack( int i ){ return (AliJMCTrack*) fMCTrackList->At(i); }
  TClonesArray * GetJMCTracks() { return fMCTrackList; }
  TClonesArray * GetJTracks() { return fTrackList; }
  AliJEventHeader * GetJEventHeader(){ return fEventHeader; }
  AliJRunHeader * GetJRunHeader(){ return fRunHeader; }

  TH1 * AddTH1(TString name, TH1*h ){
      h->SetName(name);
      h->SetTitle(name);
      h->Sumw2();
      h->SetDirectory(gDirectory);
      return h;
  }
  TH1D * AddTH1D( TString name, TH1D*h){ return (TH1D*)AddTH1(name,h); }
  TH2D * AddTH2D( TString name, TH2D*h){ return (TH2D*)AddTH1(name,h); }

 private:

  Int_t        DebugLevel(){ return 5; }
  inline void   DEBUG(int level, int type, TString msg1, TString msg2=""){
    if(DebugLevel()>level) std::cout<<msg1<<" : "<<msg2<<"\t"<<type<<std::endl;
  }

  void PrintOut() const;

  int fMBTriggMask;
  
  Bool_t fisIsolated; // weather one does isolated or normal efficiency, latter is the default
  Bool_t fisRelative; // true for relative isolation criteria, else absolutive
  double fisolParam;  // when true, isolation threshold = isolParam * pt, else isolParam
  double fisolCone;   // particle is isolated if sum pT inside cone smaller than threshold
  
  TClonesArray *    fTrackList;   //! list of charged track objects
  TClonesArray *    fMCTrackList; //! list of charged track objects
  AliJEventHeader * fEventHeader;
  AliJRunHeader*    fRunHeader; //!  run details (mg field, trigger mask,etc...)

  TH1D *fhChargedPtMC[kNVtxBin][kNCentBin];      //! all MC tracks  filled with MC pt
  TH2D *fh2DChargedPtTrigg[kNVtxBin][kNCentBin];   //! all MC track filled with MC pt in triggered event
  TH2D *fh2DChargedPtTriggVtx[kNVtxBin][kNCentBin];//! all MC track filled with MC pt in trigg event with rec vertex
  TH2D *fh2DChargedPtRec[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts];  //! [centr][cut] well reconstructed MC track filled with MC
  TH2D *fh2DChargedPtAll[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts];  //! [centr][cut] all reconstructed tracks filled with reconstructed pT 

  TH2D* fhVertexZMC; //! 
  TH2D* fhVertexZTrigg; //! 
  TH2D* fhVertexZTriggVtx; //! 

  TH1D* fhVZRawMC; //! 
  TH1D* fhVZRecMC; //! 
  TH1D* fhVZRecAccMC; //! 

  TH1D * fhChargedPtMCTrigg[kNVtxBin][kNCentBin]; //! 
  TH1D * fhChargedPtMCTriggVtx[kNVtxBin][kNCentBin]; //! 
  TH1D * fhChargedPtMCRecoCentVtx[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts][kNPrimaryStatus][kNTrackType]; //! 
  TH2D * fh2VtxCent; //! 

  TH2D * fh2MultGenRawPrimary[AliJTrackCut::kJNTrackCuts]; //! 
  TH2D * fh2MultGenRawAll[AliJTrackCut::kJNTrackCuts]; //! 

  TH1D * fhDCA2VertexXY[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts]; //! 
  TH1D * fhDCA2VertexZ[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts]; //! 

  TH1D * fhL0Input; //! 
  TH1D * fhTriggerAlice; //! 
  TH1D * fhZVtxMCAll; //! 
  TH1D * fhZVtxMCTrigg; //! 
  TH1D * fhZVtxMCTriggVtx; //! 
  TH1D * fhZVtxRecAll; //! 
  TH1D * fhZVtxRecTrigg; //! 
  TH1D * fhZVtxRecTriggVtx; //! 

  TF1 * fVtxReFunc; //! 
  TF1 * fVtxMCFunc; //! 
  TF1 * fVtxRatioFunc; //! 
  double fVtxRatioMax;

  ClassDef(AliJEfficiencyScanner, 1); 
};
#endif // AliJEfficiencyScanner_H

