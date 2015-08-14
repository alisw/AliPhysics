///
/// \file AliFemtoSimpleAnalysis.h
///

#ifndef ALIFEMTO_SIMPLE_ANALYSIS_H
#define ALIFEMTO_SIMPLE_ANALYSIS_H

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

///
/// \class AliFemtoSimpleAnalysis
/// \brief The most basic (concrete) analysis there is.
///
/// Most other analyses (e.g. AliFemtoVertexAnalysis) inherit from this one. 
/// Provides basic functionality for the analysis. To properly set up the 
/// analysis the following steps should be taken:
///
/// - create particle cuts and add them via SetFirstParticleCut and
///  SetSecondParticleCut. If one analyzes identical particle
///  correlations, the first particle cut must be also the second
///  particle cut.
///
/// - create pair cuts and add them via SetPairCut
///
/// - create one or many correlation functions and add them via
///  AddCorrFctn method.
///
/// - specify how many events are to be strored in the mixing buffer for
///  background construction
///
/// Then, when the analysis is run, for each event, the EventBegin is
/// called before any processing is done, then the ProcessEvent is called
/// which takes care of creating real and mixed pairs and sending them
/// to all the registered correlation functions. At the end of each event,
/// after all pairs are processed, EventEnd is called. After the whole
/// analysis finishes (there is no more events to process) Finish() is
/// called.
///
class AliFemtoSimpleAnalysis : public AliFemtoAnalysis {

// friend class AliFemtoLikeSignAnalysis;

public:
  AliFemtoSimpleAnalysis();
  AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& OriginalAnalysis);  /// Copy all parameters from another analysis. Cuts are cloned, mixing buffer and fNeventsProcessed are NOT copied
  virtual ~AliFemtoSimpleAnalysis();

  AliFemtoSimpleAnalysis& operator=(const AliFemtoSimpleAnalysis& aAna);  /// Copy all parameters from another analysis. Cuts are cloned from the other analysis. Mixing buffer and fNeventsProcessed are NOT copied

  // Gets and Sets

  virtual AliFemtoPairCut*       PairCut();
  virtual AliFemtoEventCut*      EventCut();
  virtual AliFemtoParticleCut*   FirstParticleCut();
  virtual AliFemtoParticleCut*   SecondParticleCut();

  AliFemtoCorrFctnCollection* CorrFctnCollection();     ///< Access to the fCorrFctnCollection
  virtual AliFemtoCorrFctn* CorrFctn(int n);            ///< Access to CFs within the collection
  void AddCorrFctn(AliFemtoCorrFctn* AnotherCorrFctn);  ///< Adds a correlation function to the fCorrFctnCollection member

  void SetPairCut(AliFemtoPairCut* ThePairCut);                         ///< Sets this analysis' pair cut, and the pair cut's analysis.
  void SetEventCut(AliFemtoEventCut* TheEventCut);                      ///< Sets this analysis' event cut, and the event cut's analysis.
  void SetFirstParticleCut(AliFemtoParticleCut* TheFirstParticleCut);   ///< Sets this analysis' particle cut 1, and the particle cut's analysis.
  void SetSecondParticleCut(AliFemtoParticleCut* TheSecondParticleCut); ///< Sets this analysis' particle cut 2, and the particle cut's analysis.

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
  virtual AliFemtoString Report();       //!<! Returns reports of all cuts applied and correlation functions being done
  virtual TList* ListSettings();         ///< return list of cut settings for the analysis
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as output
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin);    ///< Startup code each event - calls EventBegin for each member cut and correlation function
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);   ///< Bulk of analysis code - calls EventBegin, checks cuts, fills correlation functions, then calls EventEnd
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);     ///< Cleanup code each event - calls EventEnd for each member cut and correlation function
  int GetNeventsProcessed() const;      ///< Returns number of events which have been passed to ProcessEvent.

  virtual void Finish();

protected:

  void AddEventProcessed();
  void MakePairs(const char* type, 
        		 AliFemtoParticleCollection* ParticlesPassingCut1,
        		 AliFemtoParticleCollection* ParticlesPssingCut2=NULL,
        		 Bool_t enablePairMonitors=kFALSE);

  AliFemtoPicoEventCollectionVectorHideAway* fPicoEventCollectionVectorHideAway; //!<! Mixing Buffer used for Analyses which wrap this one

  AliFemtoPairCut*             fPairCut;             ///< cut applied to pairs
  AliFemtoCorrFctnCollection*  fCorrFctnCollection;  ///< correlation functions of this analysis
  AliFemtoEventCut*            fEventCut;            ///< cut to select events
  AliFemtoParticleCut*         fFirstParticleCut;    ///< select particles of type #1
  AliFemtoParticleCut*         fSecondParticleCut;   ///< select particles of type #2
  AliFemtoPicoEventCollection* fMixingBuffer;        ///< mixing buffer used in this simplest analysis
  AliFemtoPicoEvent*           fPicoEvent;           //!<! The current event, in the small (pico) form
  unsigned int fNumEventsToMix;                      ///< How many "previous" events get mixed with this one, to make background
  unsigned int fNeventsProcessed;                    ///< How many events processed so far

  unsigned int fMinSizePartCollection;               ///< Don't use event if it has fewer than this many particles passing ParticleCuts default 0

  Bool_t fVerbose;
  Bool_t fPerformSharedDaughterCut;
  Bool_t fEnablePairMonitors;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoSimpleAnalysis, 0)
  /// \endcond
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

