///
/// \file AliFemtoBaryoniaAnalysis.h
/// \author Jeremi Niedziela

#ifndef ALIFEMTO_BARYONIA_ANALYSIS_H
#define ALIFEMTO_BARYONIA_ANALYSIS_H

#include "AliFemtoTrioMinvFctn.h"
#include "AliFemtoTrio.h"

#include "AliFemtoSimpleAnalysis.h"        // base analysis class
#include "AliFemtoEventCut.h"
#include "AliFemtoParticleCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCorrFctnCollection.h"
#include "AliFemtoPicoEventCollection.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoV0SharedDaughterCut.h"

#include <list>

class AliFemtoBaryoniaAnalysis;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef std::list<AliFemtoTrioMinvFctn*, allocator<AliFemtoTrioMinvFctn*> > AliFemtoTrioFctnCollection;
typedef std::list<AliFemtoTrioMinvFctn*, allocator<AliFemtoTrioMinvFctn*> >::iterator  AliFemtoTrioFctnIterator;
#else
typedef std::list<AliFemtoTrioMinvFctn*>            AliFemtoTrioFctnCollection;
typedef std::list<AliFemtoTrioMinvFctn*>::iterator  AliFemtoTrioFctnIterator;
#endif

class AliFemtoBaryoniaAnalysis : public AliFemtoSimpleAnalysis
{
public:
  AliFemtoBaryoniaAnalysis();
  virtual ~AliFemtoBaryoniaAnalysis();
  
  // Gets and Sets
  virtual AliFemtoEventCut*      EventCut(){return fEventCut;}
  virtual AliFemtoParticleCut*   FirstParticleCut(){return fFirstParticleCut;}
  virtual AliFemtoParticleCut*   SecondParticleCut(){return fSecondParticleCut;}
  virtual AliFemtoParticleCut*   ThirdParticleCut(){return fThirdParticleCut;}
  
  // add distribution for the analysis
  void AddDistribution(AliFemtoTrioMinvFctn* function){fTrioFctnCollection->push_back(function);}
  
  void SetEventCut(AliFemtoEventCut* cut){fEventCut = cut;cut->SetAnalysis(this);}
  void SetFirstParticleCut(AliFemtoParticleCut* cut){fFirstParticleCut = cut;cut->SetAnalysis(this);}
  void SetSecondParticleCut(AliFemtoParticleCut* cut){fSecondParticleCut = cut;cut->SetAnalysis(this);}
  void SetThirdParticleCut(AliFemtoParticleCut* cut){fThirdParticleCut = cut;cut->SetAnalysis(this);}
  
  inline void SetCollection1type(AliFemtoTrio::EPart type){fCollection1type=type;}
  inline void SetCollection2type(AliFemtoTrio::EPart type){fCollection2type=type;}
  inline void SetCollection3type(AliFemtoTrio::EPart type){fCollection3type=type;}
  
  inline void SetDoEventMixing(bool mix){fDoEventMixing = mix;}
  
  void AddParticles(AliFemtoParticleCollection *collection1,
                    AliFemtoParticleCollection *collection2,
                    AliFemtoParticleCollection *collection3,
                    bool mixing);
  
  inline void SetV0SharedDaughterCut(bool perform){fPerformSharedDaughterCut = perform;}
  inline bool V0SharedDaughterCut(){return fPerformSharedDaughterCut;}
  
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as
  
  
  virtual AliFemtoString Report(){return "Report";}
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin){};
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp){};
  
  inline int GetNeventsProcessed(){return fNeventsProcessed;}
  
  virtual void Finish(){};
  virtual TList* ListSettings(){return nullptr;}
  
protected:
  AliFemtoTrioFctnCollection*  fTrioFctnCollection;  ///< correlation functions of this analysis
  AliFemtoEventCut*            fEventCut;            ///< cut to select events
  AliFemtoParticleCut*         fFirstParticleCut;    ///< select particles of type #1
  AliFemtoParticleCut*         fSecondParticleCut;   ///< select particles of type #2
  AliFemtoParticleCut*         fThirdParticleCut;    ///< select particles of type #3
  
  unsigned int fNeventsProcessed;                    ///< How many events processed so far
  int fNeventsPassed;
  AliFemtoPicoEvent *fPicoEvent;
  
  AliFemtoPicoEvent *fMixingBuffer[3];        ///< mixing buffer
  Bool_t fPerformSharedDaughterCut;
  
private:
  bool fDoEventMixing;
  
  AliFemtoTrio::EPart fCollection1type;
  AliFemtoTrio::EPart fCollection2type;
  AliFemtoTrio::EPart fCollection3type;
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoBaryoniaAnalysis, 0);
  /// \endcond
#endif
  
};

#endif
