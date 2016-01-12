#ifndef ALITRACKER_H
#define ALITRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          class AliTracker
//   that is the base for AliTPCtracker, AliITStrackerV2 and AliTRDtracker
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include "AliTrackerBase.h"

#include "AliRecoParam.h"
#include "AliPlaneEff.h"

class TTree;
class AliCluster;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliTrackPoint;
class AliKalmanTrack;
class AliEventInfo;
class TObjArray;

class AliTracker : public AliTrackerBase {
public:
  AliTracker();
  virtual ~AliTracker(){}

  virtual Int_t Clusters2Tracks(AliESDEvent *event)=0;
  virtual Int_t Clusters2TracksHLT(AliESDEvent *event, const AliESDEvent */*hltEvent*/){
    return Clusters2Tracks(event);
  }
  virtual Int_t PropagateBack(AliESDEvent *event)=0;
  virtual Int_t RefitInward(AliESDEvent *event)=0;
  virtual Int_t LoadClusters(TTree *)=0;
  virtual void UnloadClusters()=0;
  virtual AliCluster *GetCluster(Int_t index) const=0;

  virtual Int_t PostProcess(AliESDEvent */*event*/) {return 0;}
  virtual void FillClusterArray(TObjArray* array) const;
  virtual AliPlaneEff *GetPlaneEff() {return NULL;}
  virtual Bool_t GetTrackPoint(Int_t /* index */ , AliTrackPoint& /* p */) const { return kFALSE;}
  virtual Bool_t GetTrackPointTrackingError(Int_t /* index */, 
  	   AliTrackPoint& /* p */, const AliESDtrack* /* t */) { return kFALSE;}
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 

  static void FillResiduals(const AliExternalTrackParam *t,
			   Double_t *p, Double_t *cov, 
                           UShort_t id, Bool_t updated=kTRUE);
  static void FillResiduals(const AliExternalTrackParam *t,
                            const AliCluster *c, Bool_t updated=kTRUE);
  static void SetFillResiduals(AliRecoParam::EventSpecie_t es, Bool_t flag=kTRUE) { fFillResiduals=flag; fEventSpecie = es ;}
  static void SetResidualsArray(TObjArray **arr) { fResiduals=arr; }
  static TObjArray ** GetResidualsArray() { return fResiduals; }

  void                SetEventInfo(AliEventInfo *evInfo) {fEventInfo = evInfo;}
  const AliEventInfo* GetEventInfo() const {return fEventInfo;}
  //
  virtual Bool_t OwnsESDObjects() const {return kFALSE;} //RS query if tracker owns some objects in the ESD/Friends
  virtual void   CleanESDFriendsObjects(AliESDEvent*) {} //RS allow to tracker to clean the objects it ows in the friends
  virtual void   CleanESDObjects(AliESDEvent*) {} //RS allow to tracker to clean the objects it ows in the ESD
  virtual void   CleanESDTracksObjects(TObjArray* trList) {} // RS removes own objects from array of tracks
  //
  Int_t GetNumberOfClusters() const {return fNClusters;}
  //
protected:
  AliTracker(const AliTracker &atr);
private:
  AliTracker & operator=(const AliTracker & atr);
  static Bool_t fFillResiduals;     // Fill residuals flag
  static TObjArray **fResiduals;    //! Array of histograms with residuals

  static AliRecoParam::EventSpecie_t fEventSpecie ; //! event specie, see AliRecoParam
  AliEventInfo*                      fEventInfo;    //! pointer to the event info object

 protected:
  Int_t fNClusters;                                 // number of clusters loaded
  
  ClassDef(AliTracker,6) //abstract tracker
};

#endif
