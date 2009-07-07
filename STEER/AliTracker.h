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

#include <TGeoGlobalMagField.h>
#include <TObject.h>

#include "AliMagF.h"
#include "AliRecoParam.h"
#include "AliPlaneEff.h"

class TTree;
class AliCluster;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliTrackPoint;
class AliKalmanTrack;
class AliTracker : public TObject {
public:
  AliTracker();
  virtual ~AliTracker(){}
  virtual Int_t Clusters2Tracks(AliESDEvent *event)=0;
  virtual Int_t PropagateBack(AliESDEvent *event)=0;
  virtual Int_t RefitInward(AliESDEvent *event)=0;
  virtual Int_t PostProcess(AliESDEvent */*event*/) {return 0;}
  void SetVertex(const Double_t *xyz, const Double_t *ers=0) { 
     fX=xyz[0]; fY=xyz[1]; fZ=xyz[2];
     if (ers) { fSigmaX=ers[0]; fSigmaY=ers[1]; fSigmaZ=ers[2]; } 
  }

//protected:
  virtual Int_t LoadClusters(TTree *)=0;
  virtual void UnloadClusters()=0;
  virtual void FillClusterArray(TObjArray* array) const;
  virtual AliCluster *GetCluster(Int_t index) const=0;
  virtual AliPlaneEff *GetPlaneEff() {return NULL;}
  virtual Bool_t GetTrackPoint(Int_t /* index */ , AliTrackPoint& /* p */) const { return kFALSE;}
  virtual Bool_t GetTrackPointTrackingError(Int_t /* index */, 
  	   AliTrackPoint& /* p */, const AliESDtrack* /* t */) { return kFALSE;}
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetSigmaX() const {return fSigmaX;}
  Double_t GetSigmaY() const {return fSigmaY;}
  Double_t GetSigmaZ() const {return fSigmaZ;}

  static 
  Double_t MeanMaterialBudget(const Double_t *start, const Double_t *end, Double_t *mparam);
  static
  Bool_t PropagateTrackTo(AliExternalTrackParam *track, Double_t x, Double_t m,
	 Double_t maxStep, Bool_t rotateTo=kTRUE, Double_t maxSnp=0.8);  
  Bool_t PropagateTrackToBxByBz(AliExternalTrackParam *track, Double_t x, 
         Double_t m,
	 Double_t maxStep, Bool_t rotateTo=kTRUE, Double_t maxSnp=0.8);  
  //
  static Double_t GetBz(const Double_t *r);
  static void GetBxByBz(const Double_t r[3], Double_t b[3]);
  static Double_t GetBz();
  static Bool_t   UniformField();
  //
  static void FillResiduals(const AliExternalTrackParam *t,
			   Double_t *p, Double_t *cov, 
                           UShort_t id, Bool_t updated=kTRUE);
  static void FillResiduals(const AliExternalTrackParam *t,
                            const AliCluster *c, Bool_t updated=kTRUE);
  static void SetFillResiduals(AliRecoParam::EventSpecie_t es, Bool_t flag=kTRUE) { fFillResiduals=flag; fEventSpecie = es ;}
  static void SetResidualsArray(TObjArray **arr) { fResiduals=arr; }
  static TObjArray ** GetResidualsArray() { return fResiduals; }

protected:
  AliTracker(const AliTracker &atr);
private:
  AliTracker & operator=(const AliTracker & atr);
  static Bool_t fFillResiduals;       // Fill residuals flag
  static TObjArray **fResiduals;    //! Array of histograms with residuals

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex
 
  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z
  
  static AliRecoParam::EventSpecie_t fEventSpecie ; //! event specie, see AliRecoParam
  
  ClassDef(AliTracker,4) //abstract tracker
};

//__________________________________________________________________________
inline Bool_t AliTracker::UniformField()
{
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  return fld ? fld->IsUniform():kTRUE;
}

#endif
