#ifndef MYANTASKSPECTRAINEL_H
#define MYANTASKSPECTRAINEL_H

class TH1I;
class TH1F;
class TH2F;
class TList;
class TString;

class AliESDEvent;
class AliStack;

#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSigLoss : public AliAnalysisTaskSE {

public:
  enum EPileup_Type {kNoPileup, kPileupSPD};
  enum EEvtCut_Type {
    kIsReadable = 1,
    kPassTrig,
    kPassMultSel,
    kIsNotPileup,
    kHasRecVtx,
    kHasGoodVtxZ,
    kNEvtCuts
  };
  enum {
    kNspc =  9,
    kNbin = 36
  };

  AliAnalysisTaskSigLoss();
  AliAnalysisTaskSigLoss(const char *name);
  virtual ~AliAnalysisTaskSigLoss();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  void SetTriggerSel   (UInt_t   tg = AliVEvent::kINT7)  { fTriggerSel   = tg;   }
  void SetMaxVtxZCut            (Float_t  vz = 10)     { fMaxVtxZCut   = vz;   }
  void SetMinRapCut      (Double_t   y = -.5) { fMinRapCut       =   y; }
  void SetMaxRapCut      (Double_t   y =  .5) { fMaxRapCut       =   y; }
  void SetCMSRapFct      (Double_t  dy =  .0) { fCMSRapFct       =  dy; }
  void SetDoMultSel       (Bool_t flag = kTRUE) { fDoMultSel  = flag; }
  void SetMultEst(TString multEst = "V0M") {fMultEstimator = multEst;}
  void SetMultiplicityRange(Float_t low, Float_t high)
  {
    if ((high > low) && (!(low < 0.)) && (!(high > 100.))) {
      fLowMult = low; fHighMult = high;
    }
  }
  
  void SetSPDPileupSel(Int_t cont = 3, Float_t dist = 0.8);
  void SetupStandardEventCutsForRun2();
  
protected:
  Bool_t   IsEventAccepted(EEvtCut_Type& evtSel);
  Bool_t   IsMultSelected();
  Bool_t   IsPileup();
  Bool_t   HasRecVertex();
  Bool_t   IsGoodVtxZ();
  Bool_t   IsRapIn(Double_t y) {return (y > fMinRapCut && y < fMaxRapCut) ? kTRUE : kFALSE;}

  void   DoMCPart(AliStack* lMCstack, EEvtCut_Type step, Bool_t lHasGoodVtxGen);
  Double_t EtaToY(Double_t pt, Double_t m, Double_t eta) const
  {
    Double_t mt = TMath::Sqrt(m * m + pt * pt);
    if(TMath::Abs(eta)<1.) return TMath::ASinH(pt / mt * TMath::SinH(eta));
    else return -1.;
  }

private:
  AliAnalysisTaskSigLoss           (const AliAnalysisTaskSigLoss& source);
  AliAnalysisTaskSigLoss& operator=(const AliAnalysisTaskSigLoss& source);

  AliESDEvent*   fESD;      //! input ESD event
  TList*         fOutput;   //! output list in the root file
  AliEventCuts   fEventCuts;//!<! basic cut variables for events
  AliStack*      fMCstack;  //! MC stack
  AliMCEvent*    fMC;       //! MC Event
  
  TH1I* fHistNEvents;      //! histo w/ number of events
  TH1I* fHistMCEvents;     //! histo w/ number of events ZvtxGen<10 
  TH1F* fHistVtxZ;         //! histo w/ the distribution of the primary vertex Z coordinate
  TH1F* fHistMultBeforeEvtSel;  //! histo w/ event multiplicity before event selection
  TH1F* fHistMultAfterEvtSel;  //! histo w/ event multiplicity after all event selection

  TH2F* fHistPrimMCGenVtxZall[kNspc];   //! histo from events with gen Zvtx cut
  TH2F* fHistPrimMCGenVtxZallCh;        //! histo from events with gen Zvtx cut
  TH2F* fHistPrimMCGenVtxZcut[kNspc];   //! histo from events with rec Zvtx cut
  TH2F* fHistPrimMCGenVtxZcutCh;        //! histo from events with rec Zvtx cut

  Bool_t    fDoMultSel;
  TString   fMultEstimator; // multiplicity framework estimator name  
  Float_t   fLowMult;
  Float_t   fHighMult;
  Float_t   fEvtMult;       // event multiplicity -0.5 by default
  

  UInt_t    fTriggerSel;
  Float_t   fMaxVtxZCut;
  Float_t   fMinRapCut;
  Float_t   fMaxRapCut;
  Float_t   fCMSRapFct;
  
  EPileup_Type  fPlpType;

  Int_t    fMinPlpContribSPD; //minimum contributors to the pilup vertices, SPD
  Float_t  fMinPlpZdistSPD;   //minimum distance for the SPD pileup vertex

  ClassDef(AliAnalysisTaskSigLoss, 1);
};

#endif
