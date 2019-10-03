///
/// \file AliFemtoMultCorrAnalysis.h
/// \author Jeremi Niedziela

#ifndef ALIFEMTO_MULTCORR_ANALYSIS_H
#define ALIFEMTO_MULTCORR_ANALYSIS_H

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

class AliFemtoMultCorrAnalysis : public AliFemtoSimpleAnalysis
{
public:
  AliFemtoMultCorrAnalysis();
  virtual ~AliFemtoMultCorrAnalysis();
  
  // Gets and Sets
  virtual AliFemtoEventCut*      EventCut(){return fEventCut;}
  
  virtual AliFemtoParticleCut*   MotherCut(int num){return fMotherCut[num];}
  virtual AliFemtoParticleCut*   DaughterCut(int num){return fDaughterCut[num];}
  
  void SetEventCut(AliFemtoEventCut* cut){fEventCut = cut;cut->SetAnalysis(this);}
  
  void SetMotherCut(AliFemtoParticleCut* cut, int num){fMotherCut[num] = cut;cut->SetAnalysis(this);}
  void SetDaughterCut(AliFemtoParticleCut* cut, int num){fDaughterCut[num] = cut;cut->SetAnalysis(this);}
  
  virtual TList* GetOutputList();        ///< Return a TList of objects to be written as
  
  virtual AliFemtoString Report(){return "Report";}
  virtual void ProcessEvent(const AliFemtoEvent* EventToProcess);
  
  virtual void EventBegin(const AliFemtoEvent* TheEventToBegin){};
  virtual void EventEnd(const AliFemtoEvent* TheEventToWrapUp){};
  
  inline int GetNeventsProcessed(){return fNeventsProcessed;}
  
  virtual void Finish(){};
  virtual TList* ListSettings(){return nullptr;}
  
  void FillHbtParticleCollection(AliFemtoParticleCut *cut,const AliFemtoEvent *hbtEvent,
                                 AliFemtoParticleCollection *outputCollection);
  
protected:
  AliFemtoEventCut*            fEventCut;            ///< cut to select events
  
  AliFemtoParticleCut*         fMotherCut[2];    ///< select mother of type #1
  AliFemtoParticleCut*         fDaughterCut[3];    ///< select daughter of type #1
  
  unsigned int fNeventsProcessed;                    ///< How many events processed so far
  int fNeventsPassed;
  
  
private:
  TH2D *fMultCorrFctn;
  TH2D *fNmotherNdaughterCorrFctn;
  
  TH2D *fMultCorrTimesMultFctn;
  TH2D *fNmotherNdaughterRootsCorrFctn;
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoMultCorrAnalysis, 0);
  /// \endcond
#endif
  
};

#endif
