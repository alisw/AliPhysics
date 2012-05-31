#ifndef ALIESDCOSMIC_H
#define ALIESDCOSMIC_H

#include "TObject.h"

class AliESDEvent;
class TArrayI;
class TTreeSRedirector;

class AliESDcosmic:public TObject {
public:

  AliESDcosmic();
  virtual ~AliESDcosmic();
  void ProcessEvent(AliESDEvent* event);
  void PropagateToAcorde();
  void DumpToTree();
  void SetDebugStreamer(TTreeSRedirector *cstream){fDebugStreamer=cstream;}
public:
  const AliESDEvent  *fESD;            //! associated  ESD event
  TClonesArray * fTracks;        // array  of combined tracks
  TClonesArray * fTracksAcorde;  // array  of combined tracks extrapolated to acorde
  TArrayI      * fPair;           // connected track
  TTreeSRedirector *fDebugStreamer; // debug streamer
private:
  AliESDcosmic&  operator=(const AliESDcosmic&);// not implemented
  AliESDcosmic(const AliESDcosmic&); //not implemented

  ClassDef(AliESDcosmic,1)
};

#endif
