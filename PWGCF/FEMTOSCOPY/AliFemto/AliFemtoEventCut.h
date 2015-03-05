///
/// \class AliFemtoEventCut
/// \brief The pure virtual base class for the event cut
///
/// All event cuts must inherit from this one and implement the ::Pass and
/// ::Report methods. The ::Clone() function simply returns NULL, so if users
/// want their cuts to behave as expected, they should also write their own.
///

#ifndef AliFemtoEventCut_hh
#define AliFemtoEventCut_hh

class AliFemtoEvent;
class AliFemtoAnalysis;

#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoString.h"

class AliFemtoEventCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoAnalysis;

public:

  AliFemtoEventCut();                          ///< default constructor. - Users should write their own
  AliFemtoEventCut(const AliFemtoEventCut& c); ///< copy constructor
  virtual ~AliFemtoEventCut(){/* no-op */};    ///< destructor
  AliFemtoEventCut& operator=(const AliFemtoEventCut& aCut);  ///< Assignment operator

  virtual bool Pass(const AliFemtoEvent* event) =0;    ///< true if event passes, false if not

  virtual AliFemtoString Report() =0;                  ///< A user-written method to return a string describing cuts
  virtual AliFemtoEventCut* Clone() { return NULL; }   ///< Returns NULL - don't use!

  AliFemtoAnalysis* HbtAnalysis() { return fyAnalysis; };
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);

protected:
  AliFemtoAnalysis* fyAnalysis;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventCut, 0)
  /// \endcond
#endif
};

inline AliFemtoEventCut::AliFemtoEventCut(const AliFemtoEventCut& c) : AliFemtoCutMonitorHandler(), fyAnalysis(c.fyAnalysis) { }
inline void AliFemtoEventCut::SetAnalysis(AliFemtoAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoEventCut::AliFemtoEventCut(): AliFemtoCutMonitorHandler(), fyAnalysis(NULL) { }
inline AliFemtoEventCut& AliFemtoEventCut::operator=(const AliFemtoEventCut& aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; return *this; }

#endif
