////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDChainKine - the reader class for the Alice ESD and   //
// the model Kinematics information tailored for the Task framework and the   //
// Reads in AliESDfriend to create shared hit/quality information             //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERESDCHAINKINE_H
#define ALIFEMTOEVENTREADERESDCHAINKINE_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include <TTree.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliStack.h>
#include <list>
#include <AliGenEventHeader.h>

#include "AliESDpid.h"

class AliFemtoEvent;

class AliFemtoEventReaderESDChainKine : public AliFemtoEventReader 
{
 public:
  enum TrackType {kGlobal=0, kTPCOnly=1, kITSOnly=2, kSPDTracklet=3};
  typedef enum TrackType ReadTrackType;

  enum EventMult {kTracklet=0, kITSTPC=1, kITSPure=2, kGlobalCount=3, kSPDLayer1=4, kV0Centrality=5 };
  typedef enum EventMult EstEventMult;

  AliFemtoEventReaderESDChainKine();
  AliFemtoEventReaderESDChainKine(const AliFemtoEventReaderESDChainKine& aReader);
  ~AliFemtoEventReaderESDChainKine();

  AliFemtoEventReaderESDChainKine& operator=(const AliFemtoEventReaderESDChainKine& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;
  void SetUseTPCOnly(const bool usetpconly);
  bool GetUseTPCOnly() const;
  void SetReadTPCInner(const bool readinner);
  bool GetReadTPCInner() const;

  void SetUseMultiplicity(EstEventMult aType);

  void SetReadTrackType(ReadTrackType aType);

  void SetESDSource(AliESDEvent *aESD);
  void SetStackSource(AliStack *aStack);
  void SetGenEventHeader(AliGenEventHeader *aGenHeader);
  void SetRotateToEventPlane(short dorotate);

  void SetESDPid(AliESDpid *esdPid) { fESDpid = esdPid; }

 protected:

 private:
  string         fFileName;      // name of current ESD file
  bool           fConstrained;   // flag to set which momentum from ESD file will be use
  bool           fReadInner;     // flag to set if one wants to read TPC-only momentum
  bool           fUseTPCOnly;    // flog to set to read TPC only momentum instead of the full
  int            fNumberofEvent; // number of Events in ESD file
  int            fCurEvent;      // number of current event
  unsigned int   fCurFile;       // number of current file
  AliESDEvent*   fEvent;         // ESD event
  AliStack      *fStack;         // Kinematics stack pointer
  AliGenEventHeader *fGenHeader; // Link to the generator event header
  ReadTrackType  fTrackType;     // Type of track read
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator

  short          fRotateToEventPlane; // Rotate the event so that event plane is at x=0

  Float_t GetSigmaToVertex(double *impact, double *covar);

  AliESDpid *fESDpid;
  Bool_t fIsPidOwner;

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChainKine, 1)
#endif

    };
  
#endif


