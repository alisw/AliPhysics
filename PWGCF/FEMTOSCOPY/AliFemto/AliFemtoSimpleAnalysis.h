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
#include "AliFemtoV0SharedDaughterCut.h"
#include "AliFemtoXiSharedDaughterCut.h"

class AliFemtoPicoEventCollectionVectorHideAway;
class AliFemtoPicoEvent;

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

  /// Construct with default parameters
  ///
  /// All pointer members are initialized to NULL except for the correlation
  /// function collection (fCorrFctnCollection) and the mixing buffer
  /// (fMixingBuffer) which are created with default parameters.
  AliFemtoSimpleAnalysis();

  /// Copy parameters from another analysis.
  ///
  /// All parameters are copied and cuts & correlation functions are cloned.
  /// A new (empty) mixing buffer is created and the number of events processed
  /// (fNeventsProcessed) is set to 0. The EventCollectionHideAway is NOT
  /// copied, and it's up to the subclass to clone if neccessary.
  AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& OriginalAnalysis);
  AliFemtoSimpleAnalysis& operator=(const AliFemtoSimpleAnalysis& aAna);

  /// Deletes all cuts, correlation functions and events in mixing buffer.
  /// The fPicoEventCollectionVectorHideAway member is left alone.
  virtual ~AliFemtoSimpleAnalysis();

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

  /// Returns whether or not this analysis analyzes identical particles
  ///
  /// This implementation simply returns the equality of the two particle cut
  /// pointers. Even if the cuts are equivalent, this will return false if not
  /// the same memory location.
  bool AnalyzeIdenticalParticles() const;
  virtual AliFemtoString Report();       ///< Returns reports of all cuts applied and correlation functions being done
  virtual TList* ListSettings();         ///< return list of cut settings for the analysis
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as output

  /// Initialization code run at the beginning of processing an event
  ///
  /// This is implemented by calling EventBegin for each member cut
  /// and correlation function
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin);

  /// Bulk of analysis code
  ///
  /// This functions begins by calling EventBegin. If the event passes the
  /// event cut, pairs are made from the particles passing their respective
  /// cuts. The pairs are passed to each correlation function's AddRealPair
  /// method. Pairs made between particles in this event and events in the
  /// mixing buffer, are passed to the correlation functions' AddMixedPair
  /// method. The event is then added to the mixing buffer. The EventEnd() is
  /// called exactly once upon exiting this function.
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);

  /// Cleanup code after processing each event
  ///
  /// Calls EventEnd for each member cut and correlation function.
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);

  /// Returns number of events which have been passed to ProcessEvent.
  int GetNeventsProcessed() const;

  /// Calls Finish method on all correlation functions
  virtual void Finish();

protected:

  /// Increment fNeventsProcessed - is this method neccessary?
  void AddEventProcessed();

  /// Build pairs, check pair cuts, and call CFs' AddRealPair() or
  /// AddMixedPair() methods. If no second particle collection is
  /// specfied, make pairs within first particle collection.
  ///
  /// \param type Either the string "real" or "mixed", specifying which method
  ///             to call (AddRealPair or AddMixedPair)
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
  ClassDef(AliFemtoSimpleAnalysis, 0);
  /// \endcond
#endif

};

// Gets
inline AliFemtoPairCut* AliFemtoSimpleAnalysis::PairCut()
{
  return fPairCut;
}

inline AliFemtoEventCut* AliFemtoSimpleAnalysis::EventCut()
{
  return fEventCut;
}

inline AliFemtoParticleCut* AliFemtoSimpleAnalysis::FirstParticleCut()
{
  return fFirstParticleCut;
}

inline AliFemtoParticleCut* AliFemtoSimpleAnalysis::SecondParticleCut()
{
  return fSecondParticleCut;
}

inline AliFemtoCorrFctnCollection* AliFemtoSimpleAnalysis::CorrFctnCollection()
{
  return fCorrFctnCollection;
}

inline UInt_t AliFemtoSimpleAnalysis::NumEventsToMix() const
{
  return fNumEventsToMix;
}

inline AliFemtoPicoEvent* AliFemtoSimpleAnalysis::CurrentPicoEvent()
{
  return fPicoEvent;
}

inline AliFemtoPicoEventCollection* AliFemtoSimpleAnalysis::MixingBuffer()
{
  return fMixingBuffer;
}

inline bool AliFemtoSimpleAnalysis::AnalyzeIdenticalParticles() const
{
  return (fFirstParticleCut == fSecondParticleCut);
}

inline bool AliFemtoSimpleAnalysis::V0SharedDaughterCut()
{
  return fPerformSharedDaughterCut;
}

inline bool AliFemtoSimpleAnalysis::EnablePairMonitors()
{
  return fEnablePairMonitors;
}

// Sets
inline void AliFemtoSimpleAnalysis::SetPairCut(AliFemtoPairCut* x)
{
  fPairCut = x;
  x->SetAnalysis(this);
}
inline void AliFemtoSimpleAnalysis::AddCorrFctn(AliFemtoCorrFctn* cf)
{
  fCorrFctnCollection->push_back(cf);
  cf->SetAnalysis(this);
}
inline void AliFemtoSimpleAnalysis::SetEventCut(AliFemtoEventCut* x)
{
  fEventCut = x;
  x->SetAnalysis(this);
}
inline void AliFemtoSimpleAnalysis::SetFirstParticleCut(AliFemtoParticleCut* x)
{
  fFirstParticleCut = x;
  x->SetAnalysis(this);
}

inline void AliFemtoSimpleAnalysis::SetSecondParticleCut(AliFemtoParticleCut* x)
{
  fSecondParticleCut = x;
  x->SetAnalysis(this);
}

inline void AliFemtoSimpleAnalysis::SetNumEventsToMix(const unsigned int& nmix)
{
  fNumEventsToMix = nmix;
}

inline bool AliFemtoSimpleAnalysis::MixingBufferFull()
{
  return (fMixingBuffer->size() >= fNumEventsToMix);
}

inline int AliFemtoSimpleAnalysis::GetNeventsProcessed() const
{
  return fNeventsProcessed;
}

inline void AliFemtoSimpleAnalysis::SetMinSizePartCollection(unsigned int minSize)
{
  fMinSizePartCollection = minSize;
}

inline void AliFemtoSimpleAnalysis::SetVerboseMode(Bool_t aVerbose)
{
  fVerbose = aVerbose;
}

inline void AliFemtoSimpleAnalysis::SetV0SharedDaughterCut(Bool_t aPerform)
{
  fPerformSharedDaughterCut = aPerform;
}

inline void AliFemtoSimpleAnalysis::SetEnablePairMonitors(Bool_t aEnable)
{
  fEnablePairMonitors = aEnable;
}

#endif
