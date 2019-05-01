///
/// \file AliFemtoTrioAnalysis.h
/// \author Jeremi Niedziela

#ifndef ALIFEMTO_TRIO_ANALYSIS_H
#define ALIFEMTO_TRIO_ANALYSIS_H


#include "AliFemtoTrio.h"

#include "AliFemtoSimpleAnalysis.h"        // base analysis class
#include "AliFemtoEventCut.h"
#include "AliFemtoParticleCut.h"
#include "AliFemtoTrioFctn.h"
#include "AliFemtoTrioFctnCollection.h"
#include "AliFemtoPicoEventCollection.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoV0SharedDaughterCut.h"


class AliFemtoTrioAnalysis : public AliFemtoSimpleAnalysis
{
public:
  AliFemtoTrioAnalysis();
  virtual ~AliFemtoTrioAnalysis();
  
  // Gets and Sets
  virtual AliFemtoEventCut*      EventCut();
  virtual AliFemtoParticleCut*   FirstParticleCut();
  virtual AliFemtoParticleCut*   SecondParticleCut();
  virtual AliFemtoParticleCut*   ThirdParticleCut();
  
  // add three particle correlation function for the analysis
  void AddTrioFctn(AliFemtoTrioFctn* function);
  
  void SetEventCut(AliFemtoEventCut* cut);
  void SetFirstParticleCut(AliFemtoParticleCut* cut);
  void SetSecondParticleCut(AliFemtoParticleCut* cut);
  void SetThirdParticleCut(AliFemtoParticleCut* cut);
  
  void SetCollection1type(AliFemtoTrio::EPart type);
  void SetCollection2type(AliFemtoTrio::EPart type);
  void SetCollection3type(AliFemtoTrio::EPart type);
  
  void SetMinSizePart1Collection(unsigned int minSize);
  void SetMinSizePart2Collection(unsigned int minSize);
  void SetMinSizePart3Collection(unsigned int minSize);
  
  void SetDoEventMixing(bool mix);
  
  void AddParticles(AliFemtoParticleCollection *collection1,
                    AliFemtoParticleCollection *collection2=NULL,
                    AliFemtoParticleCollection *collection3=NULL,
                    bool mixing=kFALSE);
  
  void SetV0SharedDaughterCut(bool perform);
  bool V0SharedDaughterCut();
  
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as
  
  
  virtual AliFemtoString Report();
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin);
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp);
  
  int GetNeventsProcessed();
  
  virtual void Finish(){};
  virtual TList* ListSettings();
  
protected:
  AliFemtoTrioFctnCollection*  fTrioFctnCollection;  ///< correlation functions of this analysis
  AliFemtoEventCut*            fEventCut;            ///< cut to select events
  AliFemtoParticleCut*         fFirstParticleCut;    ///< select particles of type #1
  AliFemtoParticleCut*         fSecondParticleCut;   ///< select particles of type #2
  AliFemtoParticleCut*         fThirdParticleCut;    ///< select particles of type #3
  
  int fNeventsPassed;
  AliFemtoPicoEvent *fPicoEvent;
  
  AliFemtoPicoEvent *fMixingBuffer[3];        ///< mixing buffer
  Bool_t fPerformSharedDaughterCut;
  unsigned int fMinSizePart1Collection;
  unsigned int fMinSizePart2Collection;
  unsigned int fMinSizePart3Collection;///< Don't use event if it has fewer than this many particles passing ParticleCuts, for each collection separately, default 1
    
private:
  bool fDoEventMixing;
  
  AliFemtoTrio::EPart fCollection1type;
  AliFemtoTrio::EPart fCollection2type;
  AliFemtoTrio::EPart fCollection3type;
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrioAnalysis, 0);
  /// \endcond
#endif
  
};

#endif
