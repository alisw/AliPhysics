////////////////////////////////////////////////////////////////////////////////
/// AliFemtoEventCut - the pure virtual base class for the event cut         ///
/// All event cuts must inherit from this one                                ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoEventCut_hh
#define AliFemtoEventCut_hh

class AliFemtoEvent;
class AliFemtoBaseAnalysis;

#include "Infrastructure/AliFemtoCutMonitorHandler.h"
#include "Infrastructure/AliFemtoString.h"

class AliFemtoEventCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoBaseAnalysis;

public:

  AliFemtoEventCut();                // default constructor. - Users should write their own
  AliFemtoEventCut(const AliFemtoEventCut& c); // copy constructor
  virtual ~AliFemtoEventCut(){/* no-op */};       // destructor
  AliFemtoEventCut& operator=(const AliFemtoEventCut& aCut);

  virtual bool Pass(const AliFemtoEvent* event) =0;  // true if passes, false if not

  virtual AliFemtoString Report() =0;    // user-written method to return string describing cuts
  virtual AliFemtoEventCut* Clone() { return 0;}


  AliFemtoBaseAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoBaseAnalysis* aAnalysis);

protected:
  AliFemtoBaseAnalysis* fyAnalysis;

#ifdef __ROOT__
  ClassDef(AliFemtoEventCut, 0)
#endif
};

inline AliFemtoEventCut::AliFemtoEventCut(const AliFemtoEventCut& c) : AliFemtoCutMonitorHandler(), fyAnalysis(0) { }
inline void AliFemtoEventCut::SetAnalysis(AliFemtoBaseAnalysis* analysis) { fyAnalysis = analysis; }
inline AliFemtoEventCut::AliFemtoEventCut(): AliFemtoCutMonitorHandler(), fyAnalysis(0){};                // default constructor. - Users should write their own
inline AliFemtoEventCut& AliFemtoEventCut::operator=(const AliFemtoEventCut& aCut) { if (this == &aCut) return *this; fyAnalysis = aCut.fyAnalysis; return *this; }
#endif
