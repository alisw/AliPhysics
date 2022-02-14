///
/// \file AliFemto/AliFemtoEventReaderESDChain.h
///


#ifndef ALIFEMTOEVENTREADERESDCHAIN_H
#define ALIFEMTOEVENTREADERESDCHAIN_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"
#include "AliFemtoV0.h"
#include "AliESDv0.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"

#include "AliPhysicsSelection.h"

#include <TTree.h>

#include <string>
#include <vector>

class AliFemtoEvent;

/// \class AliFemtoEventReaderESDChain
/// \brief Event reader class for the Alice ESD tailored for the Task framework                                            //
///
/// Reads in AliESDfriend to create shared hit/quality information
///
/// \author Adam Kisiel kisiel@mps.ohio-state.edu
///
class AliFemtoEventReaderESDChain : public AliFemtoEventReader
{
public:
  enum TrackType
  {
    kGlobal = 0,
    kTPCOnly = 1,
    kITSOnly = 2,
    kSPDTracklet = 3
  };
  typedef enum TrackType ReadTrackType;

  enum EventMult
  {
    kCentrality = 0,
    kGlobalCount = 1,
    kReferenceITSTPC = 2,
    kReferenceITSSA = 3,
    kReferenceTracklets = 4,
    kSPDLayer1 = 5,
    kVZERO = 6,
    kCentralityTRK = 7,
    kCentralityZNA = 8,
    kCentralityCL1 = 9,
    kCentralityCND = 10,
    kCentralityV0A = 11,
    kCentralityV0C = 12,
    kCentralityZNC = 13,
    kCentralityCL0 = 14,
    kCentralityFMD = 15,
    kCentralityTKL = 16,
    kCentralityNPA = 17
  };
  typedef enum EventMult EstEventMult;

  AliFemtoEventReaderESDChain();
  AliFemtoEventReaderESDChain(const AliFemtoEventReaderESDChain &aReader);
  ~AliFemtoEventReaderESDChain();

  AliFemtoEventReaderESDChain &operator=(const AliFemtoEventReaderESDChain &aReader);

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

  void CopyESDtoFemtoEvent(AliFemtoEvent *hbtEvent);
  void SetpA2013(Bool_t pa2013);
  void SetUseMVPlpSelection(Bool_t mvplp);
  void SetIsPileUpEvent(Bool_t ispileup);
  void SetMinVtxContr(Int_t contr = 1) { fMinVtxContr = contr; }
  void SetMinPlpContribMV(Int_t minPlpContribMV) { fMinPlpContribMV = minPlpContribMV; }
  void SetMinPlpContribSPD(Int_t minPlpContribSPD) { fMinPlpContribSPD = minPlpContribSPD; }

protected:

private:

  /// name of current ESD file
  std::string fFileName;

  /// flag to set which momentum from ESD file will be use
  bool fConstrained;

  /// flag to set if one wants to read TPC-only momentum
  /// and store it in the hidden info
  bool fReadInner;

  /// flag to set if one wants to replace the global parameters
  /// by the TPC only ones
  bool fUseTPCOnly;

  /// number of Events in ESD file
  int fNumberofEvent;

  /// number of current event
  int fCurEvent;

  /// number of current file
  unsigned int fCurFile;

  /// ESD event
  AliESDEvent *fEvent;

  //  AliESDfriend*  fEventFriend;

  /// Type of track read
  ReadTrackType fTrackType;

  /// Type of the event multiplicity estimator
  EstEventMult fEstEventMult;

  ///event trigger
  UInt_t fEventTrig;

  /*   list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared */
  /*   list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow */

  Float_t GetSigmaToVertex(double *impact, double *covar);

  AliESDpid *fESDpid;
  Bool_t fIsPidOwner;

  /// Read V0 information from the AOD and put it into V0Collection
  bool fReadV0;

  ///
  int fMagFieldSign;

  Bool_t fpA2013;
  Bool_t fisPileUp;
  Bool_t fMVPlp;

  /// no of contributors for pA 2013 data
  Int_t fMinVtxContr;

  /// no of contributors for multivertex pile-up rejection
  Int_t fMinPlpContribMV;

  /// no of contributors for SPD pile-up rejection
  Int_t fMinPlpContribSPD;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderESDChain, 1);
  /// \endcond
#endif
};

#endif
