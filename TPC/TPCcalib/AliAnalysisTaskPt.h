#ifndef AliAnalysisTaskPt_cxx
#define AliAnalysisTaskPt_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
//class AliESDEvent;
class AliVfriendEvent;
class AliVEvent;
class AliESDtrackCuts;
class TList;

#include "AliAnalysisTask.h"

class AliAnalysisTaskPt : public AliAnalysisTask {
 public:
 AliAnalysisTaskPt() : AliAnalysisTask(), fESD(0), fESDfriend(0), fHistPt(0), fCuts(0), fEv(0), fHistQ(0), fListOut(0), fUseFriends(0), fHistNTPCCl(0), fHistNESDtracks(0), fHistNESDfriendtracks(0) 
 {}
  AliAnalysisTaskPt(const char *name);
  virtual ~AliAnalysisTaskPt() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Bool_t GetUseFriends() {return fUseFriends;}
  void   SetUseFriends(Bool_t flag) {fUseFriends = flag;}

 private:
  AliVEvent*       fESD;                      // ESD object
  AliVfriendEvent* fESDfriend;                // ESD friend object
  TH1F*             fHistPt;                   // Pt spectrum
  AliESDtrackCuts*  fCuts;                     // cuts
  Int_t fEv;
  TH1F*             fHistQ;                    // TPC clusters Q spectrum
  TList*            fListOut;                  // output list
  Bool_t            fUseFriends;               // flag to decide whether friends should be use
  TH1F*             fHistNTPCCl;               // histo with the number of TPC clusters
  TH1F*             fHistNESDtracks;           // histo with number of ESD tracks
  TH1F*             fHistNESDfriendtracks;     // histo with number of ESD friend tracks

  AliAnalysisTaskPt(const AliAnalysisTaskPt&); // not implemented
  AliAnalysisTaskPt& operator=(const AliAnalysisTaskPt&); // not implemented
  
  ClassDef(AliAnalysisTaskPt, 1); // example of analysis
};

#endif
