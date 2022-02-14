///
/// \file AliFemtoPairCut.h
///
/// \class AliFemtoPairCut
/// \brief The pure virtual base class for the pair cut
///
/// All pair cuts must inherit from this one
///

#ifndef AliFemtoPairCut_hh
#define AliFemtoPairCut_hh

#include <string>

class AliFemtoAnalysis;

#include "AliFemtoString.h"
#include "AliFemtoEvent.h"
#include "AliFemtoPair.h"
#include "AliFemtoCutMonitorHandler.h"
#include <TList.h>
#include <TObjString.h>

class AliFemtoPairCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoAnalysis;

public:
  enum DataType {kESD=0, kAOD=1, kKine=2};
  typedef enum DataType AliFemtoDataType;

  AliFemtoPairCut();                                 ///< default constructor. - Users should write their own
  AliFemtoPairCut(const AliFemtoPairCut& c);         ///< copy constructor
  virtual ~AliFemtoPairCut();                        ///< destructor
  virtual AliFemtoPairCut* Clone() const { return NULL; }  ///< Clones the object. The default implementation simply returns NULL

  /// Clones the object (non-const)
  /// The default implementation calls the const version
  virtual AliFemtoPairCut* Clone()
    { return const_cast<const AliFemtoPairCut*>(this)->Clone(); }

  AliFemtoPairCut& operator=(const AliFemtoPairCut &aCut);

  virtual bool Pass(const AliFemtoPair* pair) = 0;  ///< true if pair passes, false if not

  virtual AliFemtoString Report() = 0;              ///< user-written method to return string describing cuts
  virtual TList *ListSettings() = 0;                ///< Return a TList of settings

  /// the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  AliFemtoAnalysis* HbtAnalysis() { return fyAnalysis; }

  void SetAnalysis(AliFemtoAnalysis* aAnalysis);    ///< Set back-pointer to Analysis

protected:
  AliFemtoAnalysis* fyAnalysis;                     ///< Link to the base analysis class

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoPairCut, 0);
  /// \endcond
#endif
};


inline AliFemtoPairCut::AliFemtoPairCut(): AliFemtoCutMonitorHandler(), fyAnalysis(NULL) { /* no-op */ }
inline AliFemtoPairCut::AliFemtoPairCut(const AliFemtoPairCut& /* aCut */): AliFemtoCutMonitorHandler(), fyAnalysis(NULL) { /* no-op */ }
inline AliFemtoPairCut::~AliFemtoPairCut(){ /* no-op */ }

inline void AliFemtoPairCut::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoPairCut& AliFemtoPairCut::operator=(const AliFemtoPairCut &aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; return *this; }

#endif
