#ifndef AliPWG2AODTrack_H
#define AliPWG2AODTrack_H
//-------------------------------------------------------------------------
//     PWG2 specific additional information for the AOD Track
//     Author: Adam Kisiel, OSU, Adam.Kisiel@cern.ch
//-------------------------------------------------------------------------

#include <TBits.h>
#include <TRef.h>

class AliAODTrack;

class AliPWG2AODTrack : public TObject {

 public:
  
  AliPWG2AODTrack();
  AliPWG2AODTrack(Double_t tpcentr[3],
		  Double_t tpcexit[3],
		  TBits tpcshare,
		  TBits tpcclus,
		  AliAODTrack *track);

  virtual ~AliPWG2AODTrack();
  AliPWG2AODTrack(const AliPWG2AODTrack& trk); 
  AliPWG2AODTrack& operator=(const AliPWG2AODTrack& trk);

  void GetTPCNominalEntrancePoint(Double_t *tpce) const;
  void GetTPCNominalExitPoint(Double_t *tpce) const;

  void SetTPCNominalEntrancePoint(Double_t *tpce=0);
  void SetTPCNominalExitPoint(Double_t *tpce=0);

  const TBits &GetTPCSharedMap() const;
  const TBits &GetTPCClusterMap() const;
  
  void SetTPCSharedMap(const TBits &bits);
  void SetTPCClusterMap(const TBits &bits);

  void SetAODTrackRef(AliAODTrack *track);
  AliAODTrack *GetRefAODTrack();

 private :

  // TPC quality and geometrical information
  Double32_t fTPCNominalEntrancePoint[3];  // Nominal entrance point of the track to the TPC
  Double32_t fTPCNominalExitPoint[3];      // Nominal exit point of the track from the TPC
  TBits      fSharedMap;                   // TPC sharing bitmap 
  TBits      fClusterMap;                  // TPC cluster-per-padrow bitmap

  TRef       fAODTrack;                    // pointer to the original AOD track

  ClassDef(AliPWG2AODTrack,1);
};

#endif
