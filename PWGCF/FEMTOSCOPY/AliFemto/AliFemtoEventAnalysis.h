///
/// \file AliFemtoEventAnalysis.h
///

#ifndef ALIFEMTO_EVENT_ANALYSIS_H
#define ALIFEMTO_EVENT_ANALYSIS_H

#include "AliFemtoAnalysis.h"        // base analysis class
#include "AliFemtoEventCut.h"
#include "AliFemtoParticleCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCorrFctnCollection.h"
#include "AliFemtoPicoEventCollection.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoV0SharedDaughterCut.h"

class AliFemtoEventAnalysis : public AliFemtoAnalysis
{
public:

  AliFemtoEventAnalysis(double multMin, double multMax);
  AliFemtoEventAnalysis(const AliFemtoEventAnalysis& OriginalAnalysis);
  AliFemtoEventAnalysis& operator=(const AliFemtoEventAnalysis& aAna);

  virtual ~AliFemtoEventAnalysis();

  // Gets and Sets
  virtual AliFemtoEventCut*      EventCut();
  virtual AliFemtoParticleCut*   FirstParticleCut();
  virtual AliFemtoParticleCut*   SecondParticleCut();

  AliFemtoCorrFctnCollection* CorrFctnCollection();     ///< Access to the fCorrFctnCollection
  virtual AliFemtoCorrFctn* CorrFctn(int n);            ///< Access to CFs within the collection

  void AddCorrFctn(AliFemtoCorrFctn* AnotherCorrFctn);  ///< Adds a correlation function to the fCorrFctnCollection member

  void SetEventCut(AliFemtoEventCut* TheEventCut);
  void SetFirstParticleCut(AliFemtoParticleCut* TheFirstParticleCut);
  void SetSecondParticleCut(AliFemtoParticleCut* TheSecondParticleCut);

  void AddParticles(const char* typeIn, AliFemtoParticleCollection *partCollection, bool mixing = false);
  
  void SetNumEventsToMix(const unsigned int& NumberOfEventsToMix);
  void SetV0SharedDaughterCut(Bool_t aPerform);
  bool V0SharedDaughterCut();

  virtual AliFemtoString Report();
  
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin);
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);

  bool MixingBufferFull();
  
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);

  int GetNeventsProcessed() const;

  virtual void Finish();
  virtual TList* ListSettings(){return nullptr;}
  
protected:

  void AddEventProcessed();

  AliFemtoCorrFctnCollection*  fCorrFctnCollection;  ///< correlation functions of this analysis
  AliFemtoEventCut*            fEventCut;            ///< cut to select events
  AliFemtoParticleCut*         fFirstParticleCut;    ///< select particles of type #1
  AliFemtoParticleCut*         fSecondParticleCut;   ///< select particles of type #2

  unsigned int fNeventsProcessed;                    ///< How many events processed so far
  AliFemtoPicoEvent *fPicoEvent;

  AliFemtoPicoEventCollection* fMixingBuffer;        ///< mixing buffer used in this simplest analysis
  unsigned int fNumEventsToMix;                      ///< How many "previous" events get mixed with this one, to make background
  
  double fMultMin;
  double fMultMax;
  
  Bool_t fPerformSharedDaughterCut;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventAnalysis, 0);
  /// \endcond
#endif

};

inline AliFemtoEventCut* AliFemtoEventAnalysis::EventCut()
{
  return fEventCut;
}

inline AliFemtoParticleCut* AliFemtoEventAnalysis::FirstParticleCut()
{
  return fFirstParticleCut;
}

inline AliFemtoParticleCut* AliFemtoEventAnalysis::SecondParticleCut()
{
  return fSecondParticleCut;
}

inline AliFemtoCorrFctnCollection* AliFemtoEventAnalysis::CorrFctnCollection()
{
  return fCorrFctnCollection;
}

inline bool AliFemtoEventAnalysis::V0SharedDaughterCut()
{
  return fPerformSharedDaughterCut;
}

inline void AliFemtoEventAnalysis::AddCorrFctn(AliFemtoCorrFctn* cf)
{
  fCorrFctnCollection->push_back(cf);
  cf->SetAnalysis(this);
}
inline void AliFemtoEventAnalysis::SetEventCut(AliFemtoEventCut* x)
{
  fEventCut = x;
  x->SetAnalysis(this);
}
inline void AliFemtoEventAnalysis::SetFirstParticleCut(AliFemtoParticleCut* x)
{
  fFirstParticleCut = x;
  x->SetAnalysis(this);
}

inline void AliFemtoEventAnalysis::SetSecondParticleCut(AliFemtoParticleCut* x)
{
  fSecondParticleCut = x;
  x->SetAnalysis(this);
}

inline int AliFemtoEventAnalysis::GetNeventsProcessed() const
{
  return fNeventsProcessed;
}

inline void AliFemtoEventAnalysis::SetV0SharedDaughterCut(Bool_t aPerform)
{
  fPerformSharedDaughterCut = aPerform;
}

inline void AliFemtoEventAnalysis::SetNumEventsToMix(const unsigned int& nmix)
{
  fNumEventsToMix = nmix;
}

inline bool AliFemtoEventAnalysis::MixingBufferFull()
{
  return (fMixingBuffer->size() >= fNumEventsToMix);
}

#endif
