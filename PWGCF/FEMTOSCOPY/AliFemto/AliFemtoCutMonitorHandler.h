///
/// \file AliFemtoCutMonitorHandler.h
///
/// \class AliFemtoCutMonitorHandler
/// \brief A handler for cut monitors
///
/// You add cut monitors to the collection which are stored in two separate
/// collections - one which stores characteristics of the entities (tracks,
/// particles, pairs, events) that pass the respective cuts and the other for
/// the ones that fail the cut.
///

#ifndef AliFemtoCutMonitorHandler_hh
#define AliFemtoCutMonitorHandler_hh


#include "AliFemtoTypes.h"
#include "AliFemtoEvent.h"
#include "AliFemtoTrack.h"
#include "AliFemtoV0.h"
#include "AliFemtoXi.h"
#include "AliFemtoKink.h"
#include "AliFemtoPair.h" //Gael 12/04/02
#include "AliFemtoParticleCollection.h" // Gael 19/06/02
#include "AliFemtoCutMonitorCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorHandler {
public:

  AliFemtoCutMonitorHandler();
  AliFemtoCutMonitorHandler(const AliFemtoCutMonitorHandler& aHan);
  virtual ~AliFemtoCutMonitorHandler();
  AliFemtoCutMonitorHandler& operator=(const AliFemtoCutMonitorHandler& aHan);

  AliFemtoCutMonitorCollection* PassMonitorColl();
  AliFemtoCutMonitorCollection* FailMonitorColl();
  AliFemtoCutMonitor* PassMonitor(int n);
  AliFemtoCutMonitor* FailMonitor(int n);

  void AddCutMonitor(AliFemtoCutMonitor* cutMoni1, AliFemtoCutMonitor* cutMoni2);
  void AddCutMonitor(AliFemtoCutMonitor* cutMoni);
  void AddCutMonitorPass(AliFemtoCutMonitor* cutMoni);
  void AddCutMonitorFail(AliFemtoCutMonitor* cutMoni);

  void FillCutMonitor(const AliFemtoEvent* event, bool pass);
  void FillCutMonitor(const AliFemtoTrack* track, bool pass);
  void FillCutMonitor(const AliFemtoV0* v0, bool pass);
  void FillCutMonitor(const AliFemtoXi* xi, bool pass);
  void FillCutMonitor(const AliFemtoKink* kink, bool pass);
  void FillCutMonitor(const AliFemtoPair* pair, bool pass);//Gael 11/04/02
  void FillCutMonitor(const AliFemtoParticleCollection* partColl);// Gael 19/06/02
  void FillCutMonitor(const AliFemtoEvent* event, const AliFemtoParticleCollection* partColl);// Gael 19/06/02
  void FillCutMonitor(const AliFemtoParticleCollection* partColl1, const AliFemtoParticleCollection* partColl2);

  void Finish();

  virtual TList *GetOutputList();

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);

private:
  bool fCollectionsEmpty;                  ///< Are the collections empty?
  AliFemtoCutMonitorCollection* fPassColl; ///< Collection of cut monitors for passed entities
  AliFemtoCutMonitorCollection* fFailColl; ///< Collection of cut monitors for failed entities

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCutMonitorHandler, 0);
  /// \endcond
#endif

};

inline AliFemtoCutMonitorCollection* AliFemtoCutMonitorHandler::PassMonitorColl() { return fPassColl;}
inline AliFemtoCutMonitorCollection* AliFemtoCutMonitorHandler::FailMonitorColl() { return fFailColl;}

#endif
