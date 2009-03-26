#ifndef ALIMUONESDINTERFACE_H
#define ALIMUONESDINTERFACE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONESDInterface
/// \brief Converter between MUON track/cluster/digit and ESDMuon track/cluster/pad
/// 
//  Author Philippe Pillot

#include <AliMpExMap.h>
#include <TObject.h>
#include <TString.h>
#include "AliLog.h"

class AliMUONTrack;
class AliMUONVTrackStore;
class AliMUONVCluster;
class AliMUONVClusterStore;
class AliMUONVDigit;
class AliMUONVDigitStore;
class AliMUONLocalTrigger;
class AliMUONVTriggerStore;
class AliMUONTrackParam;
class AliMUONVTrackReconstructor;
class AliESDEvent;
class AliESDMuonTrack;
class AliESDMuonCluster;
class AliESDMuonPad;
class TIterator;
class AliMUONRecoParam;

class AliMUONESDInterface : public TObject
{
public: // methods to play with internal objects
  
  AliMUONESDInterface();
  virtual ~AliMUONESDInterface();
  
  virtual void Clear(Option_t* = "");
  
  void LoadEvent(AliESDEvent& esdEvent);
  
  /// Return internal track store
  AliMUONVTrackStore* GetTracks() const {return fTracks;}
  /// Return internal track store
  AliMUONVDigitStore* GetDigits() const {return fDigits;}
  /// Return internal trigger store
  AliMUONVTriggerStore* GetTriggers() const {return fTriggers;}
  
  // Return numbers of tracks/clusters/digits
  Int_t GetNTracks() const;
  Int_t GetNClusters() const;
  Int_t GetNClusters(UInt_t trackId) const;
  Int_t GetNDigits() const;
  Int_t GetNDigits(UInt_t trackId) const;
  Int_t GetNDigits(UInt_t trackId, UInt_t clusterId) const;
  Int_t GetNDigitsInCluster(UInt_t clusterId) const;
  Int_t GetNTriggers() const;
  
  // Find internal MUON objects
  AliMUONTrack*        FindTrack(UInt_t trackId) const;
  AliMUONVCluster*     FindCluster(UInt_t clusterId) const;
  AliMUONVCluster*     FindCluster(UInt_t trackId, UInt_t clusterId) const;
  AliMUONVDigit*       FindDigit(UInt_t digitId) const;
  AliMUONLocalTrigger* FindLocalTrigger(Int_t boardNumber) const;
  
  // iterate over internal MUON objects
  TIterator* CreateTrackIterator() const;
  TIterator* CreateClusterIterator() const;
  TIterator* CreateClusterIterator(UInt_t trackId) const;
  TIterator* CreateDigitIterator() const;
  TIterator* CreateDigitIterator(UInt_t trackId) const;
  TIterator* CreateDigitIterator(UInt_t trackId, UInt_t clusterId) const;
  TIterator* CreateDigitIteratorInCluster(UInt_t clusterId) const;
  TIterator* CreateLocalTriggerIterator() const;
  
  
public: // static methods
  
  /// Reset the MUON tracker (using "recoParam" if provided)
  static void ResetTracker(const AliMUONRecoParam* recoParam = 0x0);
  
  /// Set the version of track store
  static void UseTrackStore(TString name) {fgTrackStoreName = name;}
  /// Set the version of cluster store
  static void UseClusterStore(TString name) {fgClusterStoreName = name;}
  /// Set the version of digit store
  static void UseDigitStore(TString name) {fgDigitStoreName = name;}
  /// Set the version of trigger store
  static void UseTriggerStore(TString name) {fgTriggerStoreName = name;}
  
  // Create empty stores (use the version defined in this interface)
  static AliMUONVTrackStore* NewTrackStore();
  static AliMUONVClusterStore* NewClusterStore();
  static AliMUONVDigitStore* NewDigitStore();
  static AliMUONVTriggerStore* NewTriggerStore();
  
  // ESD track parameters --> MUON track parameters
  static void GetParamAtVertex(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam);
  static void GetParamAtDCA(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam);
  static void GetParamAtFirstCluster(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam);
  static void GetParamCov(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam);
  
  // MUON track parameters --> ESD track parameters
  static void SetParamAtVertex(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack);
  static void SetParamAtDCA(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack);
  static void SetParamAtFirstCluster(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack);
  static void SetParamCov(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack);
  
  // ESDMuon objects --> MUON objects conversion
  static void ESDToMUON(const AliESDMuonTrack& esdTrack, AliMUONTrack& track);
  static void ESDToMUON(const AliESDMuonTrack& esdTrack, AliMUONLocalTrigger& locTrg);
  static void ESDToMUON(const AliESDMuonCluster& esdCluster, AliMUONVCluster& cluster);
  static void ESDToMUON(const AliESDMuonPad& esdPad, AliMUONVDigit& digit);
  
  // MUON objects --> ESDMuon objects conversion
  static void MUONToESD(const AliMUONTrack& track, AliESDMuonTrack& esdTrack, const Double_t vertex[3],
			const AliMUONVDigitStore* digits = 0x0, const AliMUONLocalTrigger* locTrg = 0x0);
  static void MUONToESD(const AliMUONLocalTrigger& locTrg, AliESDMuonTrack& esdTrack, UInt_t trackId, UShort_t hitPattern);
  static void MUONToESD(const AliMUONVCluster& cluster, AliESDMuonCluster& esdCluster, const AliMUONVDigitStore* digits = 0x0);
  static void MUONToESD(const AliMUONVDigit& digit, AliESDMuonPad& esdPad);
  
  // Add ESD object into the corresponding MUON store
  // return a pointer to the corresponding MUON object into the store
  static AliMUONTrack*    Add(const AliESDMuonTrack& esdTrack, AliMUONVTrackStore& trackStore);
  static void             Add(const AliESDMuonTrack& esdTrack, AliMUONVTriggerStore& triggerStore);
  static AliMUONVCluster* Add(const AliESDMuonCluster& esdCluster, AliMUONVClusterStore& clusterStore);
  static AliMUONVDigit*   Add(const AliESDMuonPad& esdPad, AliMUONVDigitStore& digitStore);
  
  
protected:
  
  AliMUONESDInterface (const AliMUONESDInterface&); ///< copy constructor
  AliMUONESDInterface& operator=(const AliMUONESDInterface&); ///< assignment operator
  
  
private:
  
  void Reset();
  AliMUONVCluster* FindClusterInTrack(const AliMUONTrack& track, UInt_t clusterId) const;
  
  
private:
  
  static AliMUONRecoParam*           fgRecoParam; ///< reconstruction parameters for refitting
  static AliMUONVTrackReconstructor* fgTracker;   ///< track reconstructor for refitting
    
  static TString fgTrackStoreName;   ///< class name of the track store to use
  static TString fgClusterStoreName; ///< class name of the cluster store to use
  static TString fgDigitStoreName;   ///< class name of the digit store to use
  static TString fgTriggerStoreName; ///< class name of the trigger store to use
  
  // data containers
  AliMUONVTrackStore*   fTracks;   ///< track container
  AliMUONVDigitStore*   fDigits;   ///< digit container
  AliMUONVTriggerStore* fTriggers; ///< trigger container
  
  // maps (to speed up data retrieval)
  AliMpExMap* fClusterMap; ///< map of clusters
  AliMpExMap* fDigitMap;   ///< map of digits
  
  
  ClassDef(AliMUONESDInterface,0)
};

#endif

