#ifndef ALIESDTOFCLUSTER_H
#define ALIESDTOFCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

//----------------------------------------------------------------------//
//                                                                      //
// AliESDTOFcluster Class                                                //
//                                                                      //
//----------------------------------------------------------------------//

#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "TMath.h"
#include <AliVTOFcluster.h>

class AliESDTOFcluster : public AliVTOFcluster {

 public:
  AliESDTOFcluster();
  AliESDTOFcluster(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency,
		   Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu,
		   Int_t trackIndex,Float_t dX,Float_t dY,Float_t dZ,Float_t length,Double_t expTimes[9]);
  AliESDTOFcluster(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency,
		   Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu);
  AliESDTOFcluster(const AliESDTOFcluster & source);
  AliESDTOFcluster & operator=(const AliESDTOFcluster & source);
  virtual ~AliESDTOFcluster();

  Int_t Update(Int_t trackIndex,Float_t dX,Float_t dY,Float_t dZ,Float_t length,Double_t expTimes[9]);

  Int_t GetClusterIndex() const {return (*fClusterIndex)[0];} // cluster index
  Int_t GetClusterIndex(Int_t ihit) const {return (*fClusterIndex)[ihit];} // cluster index
  Int_t GetTOFchannel() const {return (*fTOFchannel)[0];} // TOF channel
  Int_t GetTOFchannel(Int_t ihit) const {return (*fTOFchannel)[ihit];} // TOF channel
  Float_t GetTime() const {return (*fTime)[0];}; // TOF time
  Float_t GetTime(Int_t ihit) const {return (*fTime)[ihit];}; // TOF time
  Float_t GetTimeRaw() const {return (*fTimeRaw)[0];}; // TOF raw time
  Float_t GetTimeRaw(Int_t ihit) const {return (*fTimeRaw)[ihit];}; // TOF raw time
  Float_t GetTOT() const {return (*fTOT)[0];}; // TOF tot
  Float_t GetTOT(Int_t ihit) const {return (*fTOT)[ihit];}; // TOF tot
  Float_t GetTOFsignalToT() const {return (*fTOT)[0];}; // TOF tot
  Float_t GetTOFsignalToT(Int_t ihit) const {return (*fTOT)[ihit];}; // TOF tot
  Int_t GetLabel(Int_t i=0) const {return i<3 ? (*fTOFlabel)[i] : -999;};
  Int_t GetLabel(Int_t i,Int_t ihit) const {return i<3 ? (*fTOFlabel)[i+3*ihit] : -999;};
  Int_t GetDeltaBC() const { return (*fDeltaBC)[0];};
  Int_t GetDeltaBC(Int_t ihit) const { return (*fDeltaBC)[ihit];};
  Int_t GetL0L1Latency() const { return (*fL0L1Latency)[0];};
  Int_t GetL0L1Latency(Int_t ihit) const { return (*fL0L1Latency)[ihit];};
  Bool_t GetStatus() const {return fStatus;};
  Float_t GetZ() const {return fZ;};
  Float_t GetPhi() const {return fPhi;};
  Float_t GetR() const {return fR;};
  Int_t GetNMatchableTracks() const {return fNmatchableTracks;};
  Int_t GetTrackIndex(Int_t i=0) const {return i<fNmatchableTracks ? fTrackIndex->At(i) : -999;};
  Float_t GetDistanceInStripPlane(Int_t i=0)   const {return i<fNmatchableTracks ? TMath::Sqrt(fDx->At(i)*fDx->At(i)+fDz->At(i)*fDz->At(i)) : -999.;}; // distance
  Float_t GetDx(Int_t i=0)  const {return i<fNmatchableTracks ? fDx->At(i) : -999.;} // distance, X component
  Float_t GetDy(Int_t i=0)  const {return i<fNmatchableTracks ? fDy->At(i) : -999.;} // distance, Y component
  Float_t GetDz(Int_t i=0)  const {return i<fNmatchableTracks ? fDz->At(i) : -999.;} // distance, Z component
  Float_t GetLength(Int_t i=0) const {return i<fNmatchableTracks ? fTrackLength->At(i) : -999.;} // reconstructed track length at TOF
  Double_t GetIntegratedTime(Int_t iPart=0,Int_t i=0) const {return (i<fNmatchableTracks && iPart<9) ? fIntegratedTimes->At(9*i+iPart) : -999.;} // reconstructed track length at TOF
  void SetStatus(Int_t status) {fStatus=status;};

  void AddTOFhit(Int_t clusterIndex,Int_t tofChannel,Float_t tofTime,Float_t timeRaw,Float_t tofTot,Int_t label[3],Int_t deltaBC,Int_t l0l1Latency, Bool_t status,Float_t zClu,Float_t phiClu,Float_t rClu);

  Int_t GetNTOFhits() const {return fNTOFhits;}

 protected:

  Int_t fNTOFhits;        // number of TOF hit in the cluster
  TArrayI *fClusterIndex; // TOF cluster index in the original tree
  TArrayI *fTOFchannel; //  TOF channel
  TArrayF *fTime; // measurement of TOF time
  TArrayF *fTimeRaw; // 
  TArrayF *fTOT; // measurement of time-over-threshould
  TArrayI *fTOFlabel; // labels of tracks (3 at maximum) that contribute to the TOF hit
  TArrayI *fDeltaBC; // 
  TArrayI *fL0L1Latency; // 
  Bool_t fStatus; // !
  Float_t fZ; // !
  Float_t fPhi; // !
  Float_t fR; // !
  Int_t fNmatchableTracks; // number of matchable tracks with the same TOF matchable hit
  TArrayI *fTrackIndex; //  index of the track in the original tree
  TArrayF *fDx; //  X component of track-cluster distance
  TArrayF *fDy; // ! Y component of track-cluster distance
  TArrayF *fDz; //  Z component of track-cluster distance
  TArrayF *fTrackLength; // receonstructed track length
  TArrayD *fIntegratedTimes; // integrated times

  ClassDef(AliESDTOFcluster, 1) // TOF matchable cluster

}; 

#endif
