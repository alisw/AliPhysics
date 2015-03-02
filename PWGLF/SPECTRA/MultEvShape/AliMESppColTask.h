#ifndef ALIMESPPCOLTASK_H
#define ALIMESPPCOLTASK_H

////////////////////////////////////////////////////////////////////////////
//  PP collision task for Multiplicity and Event Shape group                      //
//  Authors:                                                              //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIMESBASETASK_H
#include "AliMESbaseTask.h"
#endif
#define NMAXMULT 100

class AliMEStrackInfo;
class AliEventPool;
class AliEventPoolManager;
class AliMESppColTask : public AliMESbaseTask
{
public:
  class AliMESppColTaskExchange: public TObject
  {
  public:
    AliMESppColTaskExchange();
    void        Add(Float_t deta, Float_t dphi);
  private:
    Int_t   fN;
    Float_t fDEta[NMAXMULT]; //[fN]
    Float_t fDPhi[NMAXMULT]; //[fN]
    ClassDef(AliMESppColTaskExchange, 1)
  };
  class AliMESppColMixEvent: public TNamed
  {
  friend class AliMESppColTask;  
  public:
    AliMESppColMixEvent();
    AliMESppColMixEvent(AliMESppColTask *parent);
    virtual ~AliMESppColMixEvent();
    //setters
    void SetDeltaPhiInterval (Double_t min, Double_t max){
      fPhiMin = min; fPhiMax = max;
      if(TMath::Abs(fPhiMin-fPhiMax) != 2*TMath::Pi()) AliInfo("AliMESppColMixEvent::Warning: the delta phi interval is not set to 2 Pi");
          }
    void SetEventMixing(Bool_t mixON) {fmixing=mixON;}
    void SetTriggerParticleProperties(Double_t ptTrig, Double_t phiTrig, Double_t etaTrig, Int_t idTrig)
        {fPtTrigg = ptTrig; fPhiTrigg = phiTrig; fEtaTrigg = etaTrig; fidTrigg = idTrig;}
    Bool_t DefineEventPool(); // Definition of the Event pool parameters
    Bool_t InitializeEventPool(); // function that initlize everything for the analysis
    Bool_t ProcessEventPool(); // processes the event pool
    Bool_t ProcessAssociatedTracks(Int_t EvLoopIndex/*, const TObjArray* associatedTracks=NULL*/);
    Bool_t Correlate(Int_t loopindex); // function that computes the correlations between the trigger particle and the track n. loopindex
    Bool_t PoolUpdate(/*const TObjArray* associatedTracks=NULL*/);// updates the event pool
    Double_t SetCorrectPhiRange(Double_t phi); // sets all the angles in the correct range
    
//     TH2D* ApplyEventMixingCorrection(TH2D * SEhisto, TH2D * MEhisto);
//     TH2D* NormToPeak(TH2D * inputHisto)
//     TH2D* Divide2DHistos(TH2D * NumHisto, TH2D * DenomHisto)
    
    
    //getters
    AliEventPool* GetPool() {return fPool;}
    TObjArray * GetTrackArray() {return fAssociatedTracks;}
    Int_t GetNofTracks() {return fNTracks;}
    Int_t GetNofEventsInPool() {return fPoolContent;}
    Double_t GetDeltaPhi( ){return fDeltaPhi;} // Delta Phi, needs to be called after the method correlate 
    Double_t GetDeltaEta() {return fDeltaEta;} // Delta Eta
    Double_t GetMultiplicity() {return fMult;} // multiplicity

    // methods to reduce the tracks to correlate with track selection cuts applied here
  TObjArray*  AcceptAndReduceTracks(); // selecting ....

  private:
  AliMESppColMixEvent(const AliMESppColMixEvent&);
  AliMESppColMixEvent& operator=(const AliMESppColMixEvent&);
  
  protected:
  AliMESppColTask     *fParent;           // parent task
  AliEventPoolManager *fPoolManager;      // pool manager running from the nth event onwards
  AliEventPool * fPool;                   // Pool for event mixing
  TObjArray* fAssociatedTracks;          //Array of associated traks
//   AliMEStrackInfo* fMixedTrack;
  
  Bool_t fmixing;                       //switch for event mixing
  Int_t fNTracks;                       //number of tracks in the track array
  Int_t fPoolContent;                   //number of events in pool
  Double_t fMult;                       //multiplicity of event
  Double_t fPhiMin;
  Double_t fPhiMax;
  Double_t fPtTrigg;
  Double_t fEtaTrigg;
  Double_t fPhiTrigg;
  Int_t fidTrigg;
  
  Double_t fDeltaEta;
  Double_t fDeltaPhi;

  ClassDef(AliMESppColMixEvent, 1)
  };

  
  AliMESppColTask();
  AliMESppColTask(const char *name);
  virtual ~AliMESppColTask();
  
  Bool_t         BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess();

private:
  AliMESppColTask(const AliMESppColTask&);
  AliMESppColTask& operator=(const AliMESppColTask&);

  AliMESppColMixEvent *fCorrelator;     // pool manager 
//   Bool_t               fmixing;         // event mixxing switch
  
  ClassDef(AliMESppColTask, 1)            // PP collision task for the Multi Event Shape
};

#endif

