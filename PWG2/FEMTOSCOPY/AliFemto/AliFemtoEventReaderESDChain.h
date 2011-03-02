////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDChain - the reader class for the Alice ESD           //
// tailored for the Task framework                                            //
// Reads in AliESDfriend to create shared hit/quality information             //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOEVENTREADERESDCHAIN_H
#define ALIFEMTOEVENTREADERESDCHAIN_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TGraph.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliPhysicsSelection.h"
#include <list>

class AliFemtoEvent;

class AliFemtoEventReaderESDChain : public AliFemtoEventReader 
{
 public:
  enum TrackType {kGlobal=0, kTPCOnly=1, kITSOnly=2, kSPDTracklet=3};
  typedef enum TrackType ReadTrackType;

  enum EventMult {kTracklet=0, kITSTPC=1, kITSPure=2, kGlobalCount=3, kSPDLayer1=4, kV0Centrality=5 };
  typedef enum EventMult EstEventMult;

  AliFemtoEventReaderESDChain();
  AliFemtoEventReaderESDChain(const AliFemtoEventReaderESDChain& aReader);
  ~AliFemtoEventReaderESDChain();

  AliFemtoEventReaderESDChain& operator=(const AliFemtoEventReaderESDChain& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetConstrained(const bool constrained);
  void SetReadTPCInner(const bool readinner);
  void SetUseTPCOnly(const bool usetpconly);

  void SetUsePhysicsSelection(const bool usephysics);
  void SetUseMultiplicity(EstEventMult aType);

  bool GetConstrained() const;
  bool GetReadTPCInner() const;
  bool GetUseTPCOnly() const;
  
  void SetReadTrackType(ReadTrackType aType);

  void SetESDSource(AliESDEvent *aESD);
  //  void SetESDfriendSource(AliESDfriend *aFriend);

 protected:

 private:
  string         fFileName;      //name of current ESD file
  bool           fConstrained;   //flag to set which momentum from ESD file will be use
  bool           fReadInner;     // flag to set if one wants to read TPC-only momentum
                                 // and store it in the hidden info
  bool           fUseTPCOnly;    // flag to set if one wants to replace the global parameters 
                                 // by the TPC only ones
  int            fNumberofEvent; //number of Events in ESD file
  int            fCurEvent;      //number of current event
  unsigned int   fCurFile;       //number of current file
  AliESDEvent*   fEvent;         //ESD event
  //  AliESDfriend*  fEventFriend;
  bool           fUsePhysicsSel; //if true the physics selection class will be used
  AliPhysicsSelection *fSelect;  //Class to select only physics events
  ReadTrackType  fTrackType;     // Type of track read
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator

/*   list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared */
/*   list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow */
		
  Float_t GetSigmaToVertex(double *impact, double *covar);

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChain, 1)
#endif

    };
  
#endif


