///
/// \file AliFemtoTrioFctn.h
///

#pragma once

#ifndef ALIFEMTOTRIOFCTN_H
#define ALIFEMTOTRIOFCTN_H


#include "AliFemtoEvent.h"
#include "AliFemtoTrio.h"
#include "AliFemtoTrioCut.h"
//#include "AliFemtoTrioAnalysis.h"

/// \class AliFemtoTrioFctn
/// \brief The pure-virtual base class for correlation functions
///
/// All correlation function classes must inherit from this one.
///
/// This class has a optional pointers to the "parent" analysis and
/// a trio cut
///

class AliFemtoTrioAnalysis;

class AliFemtoTrioFctn {

  friend class AliFemtoTrioAnalysis;

public:

  /// Construct with null pointers to "parent" analysis and pair-cut
  AliFemtoTrioFctn();
  /// Copy constructor - copy pointer to parent analysis, clones pair-cut if it exists
  AliFemtoTrioFctn(const AliFemtoTrioFctn&);
  virtual ~AliFemtoTrioFctn(){/* no-op */};
  /// Assignment operator - shallow copy of analysis and cut
  AliFemtoTrioFctn& operator=(const AliFemtoTrioFctn&);

  virtual AliFemtoString Report() = 0;

  /// Add signal trio
  virtual void AddRealTrio(AliFemtoTrio* trio);
  /// Add background trio
  virtual void AddMixedTrio(AliFemtoTrio* trio);

  /// Not Implemented - Add pair with optional
  virtual void AddFirstParticle(AliFemtoParticle *particle, bool mixing);
  virtual void AddSecondParticle(AliFemtoParticle *particle);
  virtual void AddThirdParticle(AliFemtoParticle *particle);
  virtual void CalculateAnglesForEvent();

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  //virtual void Finish() = 0;

  virtual TList* GetOutputList() = 0;

  //virtual AliFemtoTrioFctn* Clone() const = 0;

  AliFemtoTrioAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoTrioAnalysis* aAnalysis);
  void SetTrioCut(AliFemtoTrioCut* aCut);

protected:
  AliFemtoTrioAnalysis* fyAnalysis; //! link to the analysis
  AliFemtoTrioCut* fTrioCut;    //! this is a trio criteria for three particle correlation function

  private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrioFctn, 1);
  /// \endcond
#endif
};






#endif  // ALIFEMTOCORRFCTN_H
