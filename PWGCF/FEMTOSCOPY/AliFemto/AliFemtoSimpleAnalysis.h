/**************************************************************************
 AliFemtoSimpleAnalysis - the most basic analysis there is.
 Most others (e.g. AliFemtoVertexAnalysis) wrap this one.
**************************************************************************/

#ifndef ALIFEMTO_SIMPLE_ANALYSIS_H
#define ALIFEMTO_SIMPLE_ANALYSIS_H
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoAnalysis.h"        // base analysis class
#include "AliFemtoPairCut.h"     
#include "AliFemtoEventCut.h"
#include "AliFemtoParticleCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCorrFctnCollection.h"
#include "AliFemtoPicoEventCollection.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoPicoEvent.h"
#include "AliFemtoV0SharedDaughterCut.h"

class AliFemtoPicoEventCollectionVectorHideAway;

class AliFemtoSimpleAnalysis : public AliFemtoAnalysis {

  // friend class AliFemtoLikeSignAnalysis;

 public:
  AliFemtoSimpleAnalysis();
  AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& OriginalAnalysis);  // copy constructor
  virtual ~AliFemtoSimpleAnalysis();

  AliFemtoSimpleAnalysis& operator=(const AliFemtoSimpleAnalysis& aAna);

  // Gets and Sets

  virtual AliFemtoPairCut*       PairCut();
  virtual AliFemtoEventCut*      EventCut();
  virtual AliFemtoParticleCut*   FirstParticleCut();
  virtual AliFemtoParticleCut*   SecondParticleCut();

  AliFemtoCorrFctnCollection* CorrFctnCollection();
  virtual AliFemtoCorrFctn* CorrFctn(int n);     // Access to CFs within the collection
  void AddCorrFctn(AliFemtoCorrFctn* AnotherCorrFctn);

  void SetPairCut(AliFemtoPairCut* ThePairCut);
  void SetEventCut(AliFemtoEventCut* TheEventCut);
  void SetFirstParticleCut(AliFemtoParticleCut* TheFirstParticleCut);
  void SetSecondParticleCut(AliFemtoParticleCut* TheSecondParticleCut);

  void SetMinSizePartCollection(unsigned int minSize);

  void SetVerboseMode(Bool_t aVerbose);
  
  void SetV0SharedDaughterCut(Bool_t aPerform);
  bool V0SharedDaughterCut();

  void SetEnablePairMonitors(Bool_t aEnable);
  Bool_t EnablePairMonitors();

  unsigned int NumEventsToMix() const;
  void SetNumEventsToMix(const unsigned int& NumberOfEventsToMix);
  AliFemtoPicoEvent* CurrentPicoEvent();
  AliFemtoPicoEventCollection* MixingBuffer();
  bool MixingBufferFull();

  bool AnalyzeIdenticalParticles() const;
  virtual AliFemtoString Report();       //! returns reports of all cuts applied and correlation functions being done
  virtual TList* ListSettings();         // return list of cut settings for the analysis
  virtual TList* GetOutputList();        // Return a TList of objects to be written as output
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin); // startup for EbyE
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);   // cleanup for EbyE
  int GetNeventsProcessed() const;

  virtual void Finish();

#ifdef __ROOT__
  ClassDef(AliFemtoSimpleAnalysis, 0)
#endif 


 protected:

  void AddEventProcessed();
  void MakePairs(const char* type, 
		 AliFemtoParticleCollection* ParticlesPassingCut1,
		 AliFemtoParticleCollection* ParticlesPssingCut2=0,
		 Bool_t enablePairMonitors=kFALSE);

  AliFemtoPicoEventCollectionVectorHideAway* fPicoEventCollectionVectorHideAway; //! Mixing Buffer used for Analyses which wrap this one 

  AliFemtoPairCut*             fPairCut;             //  cut applied to pairs 
  AliFemtoCorrFctnCollection*  fCorrFctnCollection;  //  correlation functions of this analysis 
  AliFemtoEventCut*            fEventCut;            //  cut to select events 
  AliFemtoParticleCut*         fFirstParticleCut;    //  select particles of type #1 
  AliFemtoParticleCut*         fSecondParticleCut;   //  select particles of type #2 
  AliFemtoPicoEventCollection* fMixingBuffer;        //  mixing buffer used in this simplest analysis 
  AliFemtoPicoEvent*           fPicoEvent;           //! The current event, in the small (pico) form 
  unsigned int fNumEventsToMix;                      //  How many "previous" events get mixed with this one, to make background 
  unsigned int fNeventsProcessed;                    //  How many events processed so far 

  unsigned int fMinSizePartCollection;               //  Don't use event if it has fewer than this many particles passing ParticleCuts default 0

  Bool_t fVerbose;
  Bool_t fPerformSharedDaughterCut;
  Bool_t fEnablePairMonitors;

#ifdef __ROOT__
  ClassDef(AliFemtoSimpleAnalysis, 0)
#endif

    };

// Get's
inline AliFemtoPairCut*              AliFemtoSimpleAnalysis::PairCut() {return fPairCut;}
inline AliFemtoEventCut*             AliFemtoSimpleAnalysis::EventCut() {return fEventCut;}
inline AliFemtoParticleCut*          AliFemtoSimpleAnalysis::FirstParticleCut() {return fFirstParticleCut;}
inline AliFemtoParticleCut*          AliFemtoSimpleAnalysis::SecondParticleCut() {return fSecondParticleCut;}
inline AliFemtoCorrFctnCollection*   AliFemtoSimpleAnalysis::CorrFctnCollection() {return fCorrFctnCollection;}
inline unsigned int                  AliFemtoSimpleAnalysis::NumEventsToMix() const {return fNumEventsToMix;}
inline AliFemtoPicoEvent*            AliFemtoSimpleAnalysis::CurrentPicoEvent() {return fPicoEvent;}

inline AliFemtoPicoEventCollection*  AliFemtoSimpleAnalysis::MixingBuffer() {return fMixingBuffer;}

inline bool AliFemtoSimpleAnalysis::AnalyzeIdenticalParticles() const { return (fFirstParticleCut==fSecondParticleCut); }
inline bool AliFemtoSimpleAnalysis::V0SharedDaughterCut() { return fPerformSharedDaughterCut; }
inline bool AliFemtoSimpleAnalysis::EnablePairMonitors() { return fEnablePairMonitors; }

// Set's
inline void AliFemtoSimpleAnalysis::SetPairCut(AliFemtoPairCut* x) { fPairCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
inline void AliFemtoSimpleAnalysis::AddCorrFctn(AliFemtoCorrFctn* cf) {fCorrFctnCollection->push_back(cf); cf->SetAnalysis((AliFemtoAnalysis*)this);}
inline void AliFemtoSimpleAnalysis::SetEventCut(AliFemtoEventCut* x) {fEventCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
inline void AliFemtoSimpleAnalysis::SetFirstParticleCut(AliFemtoParticleCut* x) {fFirstParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}
inline void AliFemtoSimpleAnalysis::SetSecondParticleCut(AliFemtoParticleCut* x) {fSecondParticleCut = x; x->SetAnalysis((AliFemtoAnalysis*)this);}

inline void AliFemtoSimpleAnalysis::SetNumEventsToMix(const unsigned int& nmix){ fNumEventsToMix = nmix;}
inline bool AliFemtoSimpleAnalysis::MixingBufferFull(){return (fMixingBuffer->size() >= fNumEventsToMix);}
inline int AliFemtoSimpleAnalysis::GetNeventsProcessed() const {return fNeventsProcessed;}

inline void AliFemtoSimpleAnalysis::SetMinSizePartCollection(unsigned int minSize){fMinSizePartCollection = minSize;}

inline void AliFemtoSimpleAnalysis::SetVerboseMode(Bool_t aVerbose){fVerbose = aVerbose;}
inline void AliFemtoSimpleAnalysis::SetV0SharedDaughterCut(Bool_t aPerform) { fPerformSharedDaughterCut = aPerform; }
inline void AliFemtoSimpleAnalysis::SetEnablePairMonitors(Bool_t aEnable) { fEnablePairMonitors = aEnable; }
#endif
