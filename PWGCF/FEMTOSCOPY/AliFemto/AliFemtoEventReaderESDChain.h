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
#include "AliFemtoV0.h"
#include "AliESDv0.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TGraph.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include <list>

#include "AliESDpid.h"

class AliFemtoEvent;

class AliFemtoEventReaderESDChain : public AliFemtoEventReader 
{
 public:
  enum TrackType {kGlobal=0, kTPCOnly=1, kITSOnly=2, kSPDTracklet=3};
  typedef enum TrackType ReadTrackType;

  enum EventMult {kCentrality=0, kGlobalCount=1, kReferenceITSTPC=2, kReferenceITSSA=3, kReferenceTracklets=4,  kSPDLayer1=5, kVZERO=6, kCentralityZNA=7, kCentralityCL1=8};
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

  virtual void CopyESDtoFemtoV0(AliESDv0 *tESDv0, AliFemtoV0 *tFemtoV0, AliESDEvent *fESDevent);
  void SetReadV0(bool a);
  void GetGlobalPositionAtGlobalRadiiThroughTPC(AliESDtrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3]);
  void SetMagneticFieldSign(int s);

  void SetUseMultiplicity(EstEventMult aType);
  void SetEventTrigger(UInt_t eventtrig); //trigger

  bool GetConstrained() const;
  bool GetReadTPCInner() const;
  bool GetUseTPCOnly() const;
  
  void SetReadTrackType(ReadTrackType aType);

  void SetESDSource(AliESDEvent *aESD);
  //  void SetESDfriendSource(AliESDfriend *aFriend);
  void SetESDPid(AliESDpid *esdPid) { fESDpid = esdPid; }
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
  ReadTrackType  fTrackType;     // Type of track read
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator
  UInt_t         fEventTrig;     //event trigger


/*   list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared */
/*   list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow */
		
  Float_t GetSigmaToVertex(double *impact, double *covar);


  AliESDpid *fESDpid;
  Bool_t fIsPidOwner;
  bool           fReadV0;           // Read V0 information from the AOD and put it into V0Collection
  int    fMagFieldSign;

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDChain, 1)
#endif

    };
  
#endif


