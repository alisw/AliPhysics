////////////////////////////////////////////////////////////////////////////////
/// AliFemtoCorrFctn - the pure virtual base class for correlation function  ///
/// All correlation function classes must inherit from this one              ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCorrFctn_hh
#define AliFemtoCorrFctn_hh

#include "AliFemtoAnalysis.h"
#include "AliFemtoEvent.h"
#include "AliFemtoPair.h"
#include "AliFemtoPairCut.h"

class AliFemtoCorrFctn{

  friend class AliFemtoAnalysis;

public:
  AliFemtoCorrFctn();
  AliFemtoCorrFctn(const AliFemtoCorrFctn& aCorrFctn);
  virtual ~AliFemtoCorrFctn(){/* no-op */};
  AliFemtoCorrFctn& operator=(const AliFemtoCorrFctn& aCorrFctn);

  virtual AliFemtoString Report() = 0;

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void AddFirstParticle(AliFemtoParticle *particle,bool mixing);
  virtual void AddSecondParticle(AliFemtoParticle *particle);
  virtual void CalculateAnglesForEvent();
  
  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual void Finish() = 0;

  virtual TList* GetOutputList() = 0;

  virtual AliFemtoCorrFctn* Clone() { return 0;}

  AliFemtoAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);
  void SetPairSelectionCut(AliFemtoPairCut* aCut);

protected:
  AliFemtoAnalysis* fyAnalysis; //! link to the analysis
  AliFemtoPairCut* fPairCut;    //! this is a PairSelection criteria for this Correlation Function

  private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctn, 1);
  /// \endcond
#endif
};

#endif
