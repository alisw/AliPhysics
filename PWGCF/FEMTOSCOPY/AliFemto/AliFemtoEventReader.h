///
/// \file AliFemtoEventReader.h
///

#ifndef ALIFEMTOEVENTREADER_H
#define ALIFEMTOEVENTREADER_H

class AliFemtoEvent;
class AliFemtoEventCut;
class AliFemtoTrackCut;
class AliFemtoV0Cut;
class AliFemtoXiCut;
class AliFemtoKinkCut;

#include "AliFemtoString.h"
#include "AliLog.h"
#include <iostream>

using namespace std;


/// \class AliFemtoEventReader
/// \brief The pure virtual base class for femto event readers.
///
/// All event readers must inherit from this one.
///
class AliFemtoEventReader {
public:

  /// Default Constuctor
  ///
  /// Even though it's only a base class and never constructed, if you don't
  /// have an implementation, you get "AliFemtoEventReader type_info node" upon
  /// dynamical loading
  ///
  /// All pointers are set to NULL, the status is set to 0 (good), and debug is
  /// set to 1 (print debug information in methods which run once)
  ///
  AliFemtoEventReader();

  /// Copy Constructor
  ///
  /// This performs a shallow copy, so both the origial and new event readers
  /// point to the same cut objects.
  AliFemtoEventReader(const AliFemtoEventReader& aReader);

  /// Destructor
  ///
  /// No members are deleted - it is up to the entity creating the cuts to
  /// delete them after the event reader has run its course
  virtual ~AliFemtoEventReader() { /* no-op */ };

  /// Assignment Operator
  /// Performs shallow copy of members
  AliFemtoEventReader& operator=(const AliFemtoEventReader& aReader);

  /// Concrete subclasses MUST implement this method, which creates the
  /// AliFemtoEvent.
  virtual AliFemtoEvent* ReturnHbtEvent() = 0;

  /// A user-written method to return a string describing the reader,
  /// including whatever "early" cuts are being done.
  virtual AliFemtoString Report();

  /// This next method does NOT need to be implemented, in which case the
  /// "default" behavior is to print a not-implemented message to stdout.
  virtual int WriteHbtEvent(AliFemtoEvent*);

  /// Init and Finish are optional but would make sense for situations like
  /// opening and closing a file
  virtual int Init(const char* ReadWrite, AliFemtoString& Message);
  virtual void Finish() { /* no-op */ };

  /// AliFemtoManager looks at this for guidance if it gets null pointer from
  /// ReturnHbtEvent
  int Status() const { return fReaderStatus; }

  virtual void SetEventCut(AliFemtoEventCut* ecut);
  virtual void SetTrackCut(AliFemtoTrackCut* pcut);
  virtual void SetV0Cut(AliFemtoV0Cut* pcut);
  virtual void SetXiCut(AliFemtoXiCut* pcut);
  virtual void SetKinkCut(AliFemtoKinkCut* pcut);
  virtual AliFemtoEventCut* EventCut();
  virtual AliFemtoTrackCut* TrackCut();
  virtual AliFemtoV0Cut* V0Cut();
  virtual AliFemtoXiCut* XiCut();
  virtual AliFemtoKinkCut* KinkCut();

  /**
   * Controls the amount of debug information printed.
   * The code indicates which functions should print debug statements:
   *
   * 0: no output at all
   * 1: methods which run once (e.g. constructor, Finsh())
   * 2: once per event
   * 3: once per track
   * 4: once per pair
   */
  int Debug() const { return fDebug; }
  void SetDebug(int d) { fDebug = d; }

protected:
  AliFemtoEventCut* fEventCut;     //!<! link to the front-loaded event cut
  AliFemtoTrackCut* fTrackCut;     //!<! link to the front-loaded track cut
  AliFemtoV0Cut* fV0Cut;           //!<! link to the front-loaded V0 cut
  AliFemtoXiCut* fXiCut;           //!<! link to the front-loaded Xi cut
  AliFemtoKinkCut* fKinkCut;       //!<! link to the front-loaded Kink cut
  int fReaderStatus;               ///< 0="good"
  int fDebug;                      ///< Debug information level

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReader, 0);
  /// \endcond
#endif
};

inline
int AliFemtoEventReader::WriteHbtEvent(AliFemtoEvent*)
{
  AliLog::Message(AliLog::kError,
                  "No WriteHbtEvent implemented",
                  "PWGCF/FEMTOSCOPY",
                  "AliFemtoEventReader",
                  "WriteHbtEvent",
                  __FILE__,
                  __LINE__);
  // Error("AliFemtoEventReader::WriteHbtEvent", "No WriteHbtEvent implemented");
  return 0;
}


#endif
