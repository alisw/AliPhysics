#ifndef AliANALYSISTASKHFCORRSIMONFLY_H
#define AliANALYSISTASKHFCORRSIMONFLY_H

/*_____________________________________________________________
 
Class AliAnalysisHFCorrOnFlySim: On the fly Simulation class for
heavy flavour correlations and general event/part properties
 
Authors:
Jitendra Kumar (jitendra.kumar@cern.ch)
Andrea Rossi   (andrea.rossi@cern.ch)
 _____________________________________________________________*/

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

class TH1I;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisHFCorrOnFlySim : public AliAnalysisTaskSE{

 public:
  AliAnalysisHFCorrOnFlySim();
  AliAnalysisHFCorrOnFlySim(const Char_t* name);
  AliAnalysisHFCorrOnFlySim(const AliAnalysisHFCorrOnFlySim& c);
  AliAnalysisHFCorrOnFlySim& operator= (const AliAnalysisHFCorrOnFlySim& c);
  virtual ~AliAnalysisHFCorrOnFlySim();

  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit(){Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
    
  void SetEtaRange(Float_t etamin, Float_t etamax){ fEtaMin=etamin; fEtaMax=etamax; }
  void SetYRange(Float_t ymin, Float_t ymax){ fYMin=ymin; fYMax=ymax; }
  void SetPtRange(Float_t ptmin, Float_t ptmax){ fPtMin=ptmin; fPtMax=ptmax; }
  void SetMultRange(Int_t Mmin, Int_t Mmax){ fMinMultiplicity=Mmin; fMaxMultiplicity=Mmax; }
  void SetQQbarCorrBetween(TString part1, Int_t charge1, TString part2, Int_t charge2)
        {fCorrPart1=part1; fChargeSel1 = charge1; fCorrPart2=part2; fChargeSel2 = charge2;}

    
  void SetPartProperties(Bool_t pYorN0){   fIsPartProp  = pYorN0; }
  void SetEventProperties(Bool_t eYorN0){  fIsEventProp = eYorN0; }
  void SetHFCorrelations(Bool_t YorN2){   fIsCorrOfHeavyFlavor    = YorN2; }
  void SetHHCorrelations(Bool_t YorN3){   fIsCorrOfHadronHadron   = YorN3; }
  void SetQQbarCorrelations(Bool_t YorN1){fIsCorrOfQQbar          = YorN1; }
    
  
 private:
  void CalculateEventProperties(TObject* obj);
  void CalculateParticleProperties(TObject *obj);
  void CalculateHFHadronCorrelations();
  void CalculateHadronHadronCorrelations(TObjArray *ParticleArray);
  void Calculate3PCorrelations();
  void CalculateQQBarCorrelations();
  
    
  TArrayI* CalculateNPartType(TString pname, Int_t &count, Int_t ChargeSel);
  void HeavyFlavourCorrelations(TObject *obj);
  void RemoveNDaughterParticleArray(TObject* obj);
  void DefineHistoNames();

  Double_t AssignCorrectPhiRange(Double_t phi){
    Double_t phiClone = 0.;
    if (phi > +1.5 * TMath::Pi())phi -= TMath::TwoPi();
    if (phi < -0.5 * TMath::Pi())phi += TMath::TwoPi();
    phiClone = phi;
    return phiClone;
  }
  
 protected:
  Bool_t IsMCEventSelected(TObject* obj);
  Bool_t IsMCParticleGenerated(TObject* obj);
  Bool_t IsMCParticleInKineCriteria(TObject* obj);

  AliMCEvent*              fMcEvent;    //! MC event
  AliInputEventHandler*    fMcHandler;  //! MCEventHandler

  TH1I  *fHistEventsProcessed;   //! histo for monitoring the number of events processed slot 1
  TList       *fOutputQA; //! Output list
  TList       *fOutputList; //! Output list
  
  Float_t fEtaMin;   // minimum eta cut
  Float_t fEtaMax;   // maximum eta cut
  Float_t fYMin;     // minimum Y cut
  Float_t fYMax;     // maximum Y cut
  Float_t fPtMin;    // minimum Pt cut
  Float_t fPtMax;    // maximum Pt cut
  Int_t fMinMultiplicity; //Max Mult Limit
  Int_t fMaxMultiplicity; //Min Mult Limit
  TString  fCorrPart1;    //Particle 1 for Corr
  TString  fCorrPart2;    //Particle 2 for Corr
  Int_t fChargeSel1;      //Charge of Part1
  Int_t fChargeSel2;      //Charge of Part2
  AliStack* fStack; 
  TObjArray* fParticleArray;
  Int_t  fcounter;
  Bool_t fIsEventProp;
  Bool_t fIsPartProp;
  Bool_t fIsCorrOfQQbar;
  Bool_t fIsCorrOfHeavyFlavor;
  Bool_t fIsCorrOfHadronHadron;

  TArrayI *fArraySkipDDaugh;//!
  TArrayI *fArrayTrk ;//!
  Int_t flastdaugh;
  ClassDef(AliAnalysisHFCorrOnFlySim,1)
};

#endif
