/**************************************************************************
 AliFemtoAnalysis - the most basic analysis there is.
 Most others (e.g. AliFemtoVertexAnalysis) wrap this one.
**************************************************************************/

#ifndef AliFemtoAnalysis_hh
#define AliFemtoAnalysis_hh
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoBaseAnalysis.h"        // base analysis class
#include "Base/AliFemtoPairCut.h"     
#include "Base/AliFemtoEventCut.h"
#include "Base/AliFemtoParticleCut.h"
#include "Base/AliFemtoCorrFctn.h"
#include "Infrastructure/AliFemtoCorrFctnCollection.h"
#include "Infrastructure/AliFemtoPicoEventCollection.h"
#include "Infrastructure/AliFemtoParticleCollection.h"
#include "Infrastructure/AliFemtoPicoEvent.h"

class AliFemtoPicoEventCollectionVectorHideAway;


class AliFemtoAnalysis : public AliFemtoBaseAnalysis {

 friend class AliFemtoLikeSignAnalysis;

 public:
  AliFemtoAnalysis();
  AliFemtoAnalysis(const AliFemtoAnalysis& OriginalAnalysis);  // copy constructor
  virtual ~AliFemtoAnalysis();

  AliFemtoAnalysis& operator=(const AliFemtoAnalysis& aAna);

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

  unsigned int NumEventsToMix();
  void SetNumEventsToMix(const unsigned int& NumberOfEventsToMix);
  AliFemtoPicoEvent* CurrentPicoEvent();
  AliFemtoPicoEventCollection* MixingBuffer();
  bool MixingBufferFull();

  bool AnalyzeIdenticalParticles();
  virtual AliFemtoString Report();       //! returns reports of all cuts applied and correlation functions being done

  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin); // startup for EbyE
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);   // cleanup for EbyE
  int GetNeventsProcessed();

  virtual void Finish();

 protected:

  void AddEventProcessed();
  void MakePairs(const char* type, 
		 AliFemtoParticleCollection* ParticlesPassingCut1,
		 AliFemtoParticleCollection* ParticlesPssingCut2=0);

  AliFemtoPicoEventCollectionVectorHideAway* fPicoEventCollectionVectorHideAway;  /* Mixing Buffer used for Analyses which wrap this one */

  AliFemtoPairCut*             fPairCut;             /* cut applied to pairs */
  AliFemtoCorrFctnCollection*  fCorrFctnCollection;  /* correlation functions of this analysis */
  AliFemtoEventCut*            fEventCut;            /* cut to select events */
  AliFemtoParticleCut*         fFirstParticleCut;    /* select particles of type #1 */
  AliFemtoParticleCut*         fSecondParticleCut;   /* select particles of type #2 */
  AliFemtoPicoEventCollection* fMixingBuffer;        /* mixing buffer used in this simplest analysis */
  AliFemtoPicoEvent*           fPicoEvent;           /* The current event, in the small (pico) form */
  unsigned int fNumEventsToMix;                      /* How many "previous" events get mixed with this one, to make background */
  unsigned int fNeventsProcessed;                    /* How many events processed so far */

  unsigned int fMinSizePartCollection;               /* Don't use event if it has fewer than this many particles passing ParticleCuts default 0*/

#ifdef __ROOT__
  ClassDef(AliFemtoAnalysis, 0)
#endif

    };

// Get's
inline AliFemtoPairCut*             AliFemtoAnalysis::PairCut() {return fPairCut;}
inline AliFemtoEventCut*            AliFemtoAnalysis::EventCut() {return fEventCut;}
inline AliFemtoParticleCut*         AliFemtoAnalysis::FirstParticleCut() {return fFirstParticleCut;}
inline AliFemtoParticleCut*         AliFemtoAnalysis::SecondParticleCut() {return fSecondParticleCut;}
inline AliFemtoCorrFctnCollection*  AliFemtoAnalysis::CorrFctnCollection() {return fCorrFctnCollection;}
inline unsigned int              AliFemtoAnalysis::NumEventsToMix(){return fNumEventsToMix;}
inline AliFemtoPicoEvent*           AliFemtoAnalysis::CurrentPicoEvent() {return fPicoEvent;}

inline AliFemtoPicoEventCollection*  AliFemtoAnalysis::MixingBuffer() {return fMixingBuffer;}

// Set's
inline bool AliFemtoAnalysis::AnalyzeIdenticalParticles(){return (fFirstParticleCut==fSecondParticleCut);}
inline void AliFemtoAnalysis::SetPairCut(AliFemtoPairCut* x) { fPairCut = x; x->SetAnalysis((AliFemtoBaseAnalysis*)this);}
inline void AliFemtoAnalysis::AddCorrFctn(AliFemtoCorrFctn* cf) {fCorrFctnCollection->push_back(cf); cf->SetAnalysis((AliFemtoBaseAnalysis*)this);}
inline void AliFemtoAnalysis::SetEventCut(AliFemtoEventCut* x) {fEventCut = x; x->SetAnalysis((AliFemtoBaseAnalysis*)this);}
inline void AliFemtoAnalysis::SetFirstParticleCut(AliFemtoParticleCut* x) {fFirstParticleCut = x; x->SetAnalysis((AliFemtoBaseAnalysis*)this);}
inline void AliFemtoAnalysis::SetSecondParticleCut(AliFemtoParticleCut* x) {fSecondParticleCut = x; x->SetAnalysis((AliFemtoBaseAnalysis*)this);}

inline void AliFemtoAnalysis::SetNumEventsToMix(const unsigned int& nmix){ fNumEventsToMix = nmix;}
inline bool AliFemtoAnalysis::MixingBufferFull(){return (fMixingBuffer->size() >= fNumEventsToMix);}
inline int AliFemtoAnalysis::GetNeventsProcessed() {return fNeventsProcessed;}

inline void AliFemtoAnalysis::SetMinSizePartCollection(unsigned int minSize){fMinSizePartCollection = minSize;}

#endif
