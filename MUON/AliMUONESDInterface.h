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
#include <TClonesArray.h>
#include <TString.h>

class AliMUONTrack;
class AliMUONVTrackStore;
class AliMUONVCluster;
class AliMUONVClusterStore;
class AliMUONVDigit;
class AliMUONVDigitStore;
class AliMUONTrackParam;
class AliESDEvent;
class AliESDMuonTrack;
class AliESDMuonCluster;
class AliESDMuonPad;
class TIterator;

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
  
  // Return numbers of tracks/clusters/digits
  Int_t GetNTracks() const;
  Int_t GetNClusters() const;
  Int_t GetNClusters(Int_t iTrack) const;
  Int_t GetNDigits() const;
  Int_t GetNDigits(Int_t iTrack) const;
  Int_t GetNDigits(Int_t iTrack, Int_t iCluster) const;
  Int_t GetNDigitsInCluster(UInt_t clusterId) const;
  
  // Return internal MUON objects (faster than finders)
  // ordering of the MUON objects is the same as in ESD
  AliMUONTrack*    GetTrack(Int_t iTrack) const;
  AliMUONVCluster* GetCluster(Int_t iTrack, Int_t iCluster) const;
  AliMUONVDigit*   GetDigit(Int_t iTrack, Int_t iCluster, Int_t iDigit) const;
  
  // Quickly return internal MUON objects (indices unchecked)
  // ordering of the MUON objects is the same as in ESD
  AliMUONTrack*    GetTrackFast(Int_t iTrack) const;
  AliMUONVCluster* GetClusterFast(Int_t iTrack, Int_t iCluster) const;
  AliMUONVDigit*   GetDigitFast(Int_t iTrack, Int_t iCluster, Int_t iDigit) const;
  
  // Find internal MUON objects (slower than getters)
  AliMUONVCluster* FindCluster(UInt_t clusterId) const;
  AliMUONVDigit*   FindDigit(UInt_t digitId) const;
  
  // iterate over internal MUON objects
  TIterator* CreateTrackIterator() const;
  TIterator* CreateClusterIterator() const;
  TIterator* CreateClusterIterator(Int_t iTrack) const;
  TIterator* CreateDigitIterator() const;
  TIterator* CreateDigitIterator(Int_t iTrack) const;
  TIterator* CreateDigitIterator(Int_t iTrack, Int_t iCluster) const;
  TIterator* CreateDigitIteratorInCluster(UInt_t clusterId) const;
  
  
public: // static methods
  
  /// Set the version of track store
  static void UseTrackStore(TString name) {fgTrackStoreName = name;}
  /// Set the version of cluster store
  static void UseClusterStore(TString name) {fgClusterStoreName = name;}
  /// Set the version of digit store
  static void UseDigitStore(TString name) {fgDigitStoreName = name;}
  
  // Create empty stores (use the version defined in this interface)
  static AliMUONVTrackStore* NewTrackStore();
  static AliMUONVClusterStore* NewClusterStore();
  static AliMUONVDigitStore* NewDigitStore();
  
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
  static void ESDToMUON(const AliESDMuonCluster& esdCluster, AliMUONVCluster& cluster);
  static void ESDToMUON(const AliESDMuonPad& esdPad, AliMUONVDigit& digit);
  
  // MUON objects --> ESDMuon objects conversion
  static void MUONToESD(const AliMUONTrack& track, AliESDMuonTrack& esdTrack, const Double_t vertex[3], const AliMUONVDigitStore* digits = 0x0);
  static void MUONToESD(const AliMUONVCluster& cluster, AliESDMuonCluster& esdCluster, const AliMUONVDigitStore* digits = 0x0);
  static void MUONToESD(const AliMUONVDigit& digit, AliESDMuonPad& esdPad);
  
  // Add ESD object into the corresponding MUON store
  // return a pointer to the corresponding MUON object into the store
  static AliMUONTrack*    Add(const AliESDMuonTrack& esdTrack, AliMUONVTrackStore& trackStore);
  static AliMUONVCluster* Add(const AliESDMuonCluster& esdCluster, AliMUONVClusterStore& clusterStore);
  static AliMUONVDigit*   Add(const AliESDMuonPad& esdPad, AliMUONVDigitStore& digitStore);
  
  
protected:
  
  AliMUONESDInterface (const AliMUONESDInterface&); ///< copy constructor
  AliMUONESDInterface& operator=(const AliMUONESDInterface&); ///< assignment operator
  
  
private:
  
  void Reset();
  AliMUONVCluster* FindClusterInTrack(const AliMUONTrack& track, UInt_t clusterId) const;
  
  
private:
  
  static TString fgTrackStoreName;   ///< class name of the track store to use
  static TString fgClusterStoreName; ///< class name of the cluster store to use
  static TString fgDigitStoreName;   ///< class name of the digit store to use
  
  // data containers
  AliMUONVTrackStore* fTracks; ///< track container
  AliMUONVDigitStore* fDigits; ///< digit container
  
  // maps (to speed up data retrieval)
  AliMpExMap*   fTrackMap;   ///< map of tracks
  AliMpExMap*   fClusterMap; ///< map of clusters
  TClonesArray* fDigitMap;   ///< map of digits
    
    
  ClassDef(AliMUONESDInterface,0)
};


//___________________________________________________________________________
inline AliMUONTrack* AliMUONESDInterface::GetTrackFast(Int_t iTrack) const
{
  /// return MUON track "iTrack" without any check
  return (AliMUONTrack*) fTrackMap->GetObjectFast(iTrack);
}

//___________________________________________________________________________
inline AliMUONVCluster* AliMUONESDInterface::GetClusterFast(Int_t iTrack, Int_t iCluster) const
{
  /// return MUON cluster numbered "iCluster" in track "iTrack" without any check
  return (AliMUONVCluster*) ((AliMpExMap*) fClusterMap->GetObjectFast(iTrack))->GetObjectFast(iCluster);
}

//___________________________________________________________________________
inline AliMUONVDigit* AliMUONESDInterface::GetDigitFast(Int_t iTrack, Int_t iCluster, Int_t iDigit) const
{
  /// return MUON digit numbered "iDigit" in cluster numbered "iCluster" of track "iTrack" without any check
  return (AliMUONVDigit*) ((AliMpExMap*) ((AliMpExMap*) fDigitMap->UncheckedAt(iTrack))->GetObjectFast(iCluster))->GetObjectFast(iDigit);
}


#endif

