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
#include <TObject.h>
#include <AliPlaneEff.h>

class TTree;
class AliMagF;
class AliCluster;
class AliKalmanTrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliTrackPoint;

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

  static void SetFieldMap(const AliMagF* map, Bool_t uni);
  static const AliMagF *GetFieldMap() {return fgkFieldMap;}
  static Double_t GetBz(const Float_t *r); 
  static Double_t GetBz(const Double_t *r) {
    Float_t rr[]={r[0],r[1],r[2]};
    return GetBz(rr);
  }
  static Double_t GetBz() {return fgBz;}
  static Bool_t UniformField() {return fgUniformField;}

  static void FillResiduals(const AliExternalTrackParam *t,
			   Double_t *p, Double_t *cov, 
                           UShort_t id, Bool_t updated=kTRUE);
  static void SetFillResiduals(Bool_t flag=kTRUE) { fFillResiduals=flag; }
  static void SetResidualsArray(TObjArray *arr) { fResiduals=arr; }

protected:
  AliTracker(const AliTracker &atr);
private:
  AliTracker & operator=(const AliTracker & atr);

  static Bool_t fgUniformField;       // uniform field flag
  static const AliMagF *fgkFieldMap;  //! field map
  static Double_t fgBz;               // Nominal Bz (kG)

  static Bool_t fFillResiduals;       // Fill residuals flag
  static TObjArray *fResiduals;    //! Array of histograms with residuals

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex
 
  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z

  ClassDef(AliTracker,4) //abstract tracker
};

#endif


