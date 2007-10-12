////////////////////////////////////////////////////////////////////////////////
/// AliFemtoEventReader - the pure virtual base class for the event reader   ///
/// All event readers must inherit from this one                             ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADER_H
#define ALIFEMTOEVENTREADER_H

class AliFemtoEvent;
class AliFemtoEventCut;
class AliFemtoTrackCut;
class AliFemtoV0Cut;
class AliFemtoXiCut;
class AliFemtoKinkCut;

#include "AliFemtoString.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

class AliFemtoEventReader {
  
 public:
  // even tho it's only a base class and never constructed, if you don't have an implementation,
  // you get "AliFemtoEventReader type_info node" upon dynamical loading
  AliFemtoEventReader() : fEventCut(0), fTrackCut(0), fV0Cut(0), fXiCut(0), fKinkCut(0), fReaderStatus(0), fDebug(1) { /* no-op */ };
  AliFemtoEventReader(const AliFemtoEventReader& aReader);
  virtual ~AliFemtoEventReader(){/* no-op */}
  
  AliFemtoEventReader& operator=(const AliFemtoEventReader& aReader);
  
  virtual AliFemtoEvent* ReturnHbtEvent() =0;

  virtual AliFemtoString Report();    // user-written method to return string describing reader
                                      // Including whatever "early" cuts are being done

  // this next method does NOT need to be implemented, in which case the 
  // "default" method below is executed
  virtual int WriteHbtEvent(AliFemtoEvent*){cout << "No WriteHbtEvent implemented\n"; return (0);}

  // these next two are optional but would make sense for, e.g., opening and closing a file
  virtual int Init(const char* ReadWrite, AliFemtoString& Message){cout << "do-nothing AliFemtoEventReader::Init()\n"; return(0);}
  virtual void Finish(){/*no-op*/};

  int Status() const {return fReaderStatus;} // AliFemtoManager looks at this for guidance if it gets null pointer from ReturnHbtEvent

  virtual void SetEventCut(AliFemtoEventCut* ecut);
  virtual void SetTrackCut(AliFemtoTrackCut* pcut);
  virtual void SetV0Cut(AliFemtoV0Cut* pcut);
  virtual void SetXiCut(AliFemtoXiCut* pcut);
  virtual void SetKinkCut(AliFemtoKinkCut* pcut);
  virtual AliFemtoEventCut* EventCut();
  virtual AliFemtoTrackCut* TrackCut();
  virtual AliFemtoV0Cut*    V0Cut();
  virtual AliFemtoXiCut*    XiCut();
  virtual AliFemtoKinkCut*    KinkCut();

  /* control of debug informations print out, my rule is: */
  /* 0: no output at all                                  */
  /* 1: once (e.g. in constructor, finsh                  */
  /* 2: once per event                                    */
  /* 3: once per track                                    */
  /* 4: once per pair                                     */
  int Debug() const {return fDebug;} 
  void SetDebug(int d){fDebug=d;}

protected:
  AliFemtoEventCut* fEventCut;     //! link to the front-loaded event cut
  AliFemtoTrackCut* fTrackCut;     //! link to the front-loaded track cut
  AliFemtoV0Cut* fV0Cut;           //! link to the front-loaded V0 cut
  AliFemtoXiCut* fXiCut;           //! link to the front-loaded Xi cut
  AliFemtoKinkCut* fKinkCut;       //! link to the front-loaded Kink cut
  int fReaderStatus;               // 0="good"
  int fDebug;                      // Debug information level
#ifdef __ROOT__
  ClassDef(AliFemtoEventReader,0)
#endif
};


#endif

