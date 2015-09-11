///
/// \file AliFemtoCutMonitor.h
///
/// \class AliFemtoCutMonitor
/// \brief The base class for cut monitors
///
/// A cut monitor saves attributes of the entities that have passed or failed
/// the given cut. This class is the base class which is present to
/// provide a common interface for storing data.
///
/// Cut monitors are to be used in conjunction with cut monitor handlers
/// (AliFemtoCutMonitorHandler) - of which all standard cuts (e.g.
/// AliFemtoEventCut) inherit from. Your cut monitor objects get added to the
/// monitorhandlers via their AddCutMonitor methods, and the Fill commands get
/// called by the handler upon their FillCutMonitor method, with the particular
/// entity type.
///
/// The default behavior of this base class is to do nothing with the
/// incoming data, as no data members are provided. It is up to the user to use
/// (or write) a subclass with relevant histograms.
///
/// To implement a custom cut monitor, subclass this class and overload the
/// Fill method(s) corresponding to the entity type(s) you wish to monitor.
/// All other 'Fill' methods should be implemented to avoid 'member hiding'
/// compiler warnings.
///
/// All methods of this class are empty except report which returns an empty
/// string and getOutputList which returns a pointer to an empty list.
///

#ifndef ALIFEMTOCUTMONITOR_H
#define ALIFEMTOCUTMONITOR_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoXi;
class AliFemtoKink;
class AliFemtoPair;

#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"

#include <TList.h>

class AliFemtoCutMonitor {
public:
  AliFemtoCutMonitor();
  virtual ~AliFemtoCutMonitor();

  /// Returns an empty string
  virtual AliFemtoString Report();

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);

  /// Returns pointer to empty list
  virtual TList* GetOutputList();

  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack);
  virtual void Fill(const AliFemtoV0* aV0);
  virtual void Fill(const AliFemtoXi* aXi);
  virtual void Fill(const AliFemtoKink* aKink);
  virtual void Fill(const AliFemtoPair* aPair);

  virtual void Fill(const AliFemtoParticleCollection* aCollection);
  virtual void Fill(const AliFemtoEvent* aEvent,
                    const AliFemtoParticleCollection* aCollection);
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,
                    const AliFemtoParticleCollection* aCollection2);

  virtual void Finish();
  virtual void Init();
};

inline
AliFemtoCutMonitor::AliFemtoCutMonitor()
{ /* no-op */
}

inline
AliFemtoCutMonitor::~AliFemtoCutMonitor()
{ /* no-op */
}

inline
AliFemtoString AliFemtoCutMonitor::Report()
{
  AliFemtoString report("");
  return report;
}

inline
TList* AliFemtoCutMonitor::GetOutputList()
{
  TList *tOutputList = new TList();
  return tOutputList;
}

inline
void AliFemtoCutMonitor::Finish()
{ /* no-op */
}

inline
void AliFemtoCutMonitor::Init()
{ /* no-op */
}


#endif
