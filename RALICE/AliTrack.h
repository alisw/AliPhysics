#ifndef ALITRACK_H
#define ALITRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$


#include "TObject.h"
#include "TObjArray.h"
#include "TArrayD.h"
#include "TArrayI.h"
 
#include "AliSignal.h"
#include "AliBoost.h"
#include "AliPositionObj.h"
#include "AliTimestamp.h"
 
class AliTrack : public TNamed,public Ali4Vector
{
 public:
  AliTrack();                             // Default constructor
  virtual ~AliTrack();                    // Destructor
  AliTrack(const AliTrack& t);            // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  virtual void Reset();                   // Reset all values to 0
  void Set4Momentum(Ali4Vector& p);       // Set track 4-momentum
  void Set3Momentum(Ali3Vector& p);       // Set track 3-momentum
  void SetMass(Double_t m,Double_t dm=0); // Set particle mass and error
  void SetMass();                         // Set mass and error to the values of the hyp. with highest prob.
  void SetCharge(Float_t q);              // Set particle charge
  virtual void Data(TString f="car",TString u="rad"); // Print track information for frame f and ang units u
  virtual void List(TString f="car",TString u="rad"); // Track and decay level 1 info for frame f and ang units u
  virtual void ListAll(TString f="car",TString u="rad");// Track and all decay level info for frame f and ang units u
  Ali3Vector Get3Momentum(Float_t scale=-1) const; // Provide track 3-momentum
  Double_t GetMomentum(Float_t scale=-1); // Provide value of track 3-momentum
  Double_t GetMass(Float_t scale=-1);     // Provide particle mass
  Float_t GetCharge() const;              // Provide particle charge
  Double_t GetEnergy(Float_t scale=-1);   // Provide particle total energy
  void Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms); // Perform 2-body decay
  Int_t GetNdecay() const;                // Provide number of decay products
  AliTrack* GetDecayTrack(Int_t j) const; // Access to decay produced track number j
  void RemoveDecays();                    // Remove all the decay products of this track
  void AddSignal(AliSignal& s,Int_t mode=0);    // Relate an AliSignal to this track
  void RemoveSignal(AliSignal& s,Int_t mode=1); // Remove related AliSignal from this track
  void RemoveSignals(Int_t mode=1);             // Remove all related AliSignals from this track
  Int_t GetNsignals() const;              // Provide number of related AliSignals
  AliSignal* GetSignal(Int_t j) const;    // Access to the related AliSignal number j
  void SetBeginPoint(AliPosition& p);     // Set the track begin-point
  AliPosition* GetBeginPoint();           // Provide the track begin-point
  void SetEndPoint(AliPosition& p);       // Set the track end-point
  AliPosition* GetEndPoint();             // Provide the track end-point
  void SetReferencePoint(AliPosition& p); // Set the track reference-point for the 3-momentum vector
  AliPosition* GetReferencePoint();       // Provide the track reference-point for the 3-momentum vector
  void AddTrackHypothesis(AliTrack& t);   // Add track hypothesis
  void AddTrackHypothesis(Double_t prob,Double_t m,Double_t dm=0); // Add track hypothesis with mass data
  Int_t GetNhypotheses() const;           // Provide number of track hypotheses
  AliTrack* GetTrackHypothesis(Int_t j=0) const; // Provide the j-th track hypothesis 
  void RemoveTrackHypothesis(AliTrack& t);// Remove the specified track hypothesis 
  void RemoveTrackHypotheses();           // Remove all track hypotheses 
  Double_t GetPt(Float_t scale=-1);       // Provide trans. momentum w.r.t. z-axis
  Double_t GetPl(Float_t scale=-1);       // Provide long. momentum w.r.t. z-axis
  Double_t GetEt(Float_t scale=-1);       // Provide trans. energy w.r.t. z-axis
  Double_t GetEl(Float_t scale=-1);       // Provide long. energy w.r.t. z-axis
  Double_t GetMt(Float_t scale=-1);       // Provide trans. mass w.r.t. z-axis
  Double_t GetRapidity();                 // Provide rapidity value w.r.t. z-axis
  void SetImpactPoint(AliPosition& p,TString q); // Set the impact-point in plane "q=0"
  AliPosition* GetImpactPoint(TString q);        // Provide the impact-point in plane "q=0"
  void SetId(Int_t id);                   // Set the user defined unique track identifier
  Int_t GetId() const;                    // Provide the user defined unique track identifier
  void SetClosestPoint(AliPosition& p);   // Set position p as point of closest approach w.r.t. some reference
  AliPosition* GetClosestPoint();         // Provide point of closest approach w.r.t. some reference
  void SetParticleCode(Int_t code);       // Set the user defined particle id code (e.g. the PDF convention)
  Int_t GetParticleCode() const;          // Provide the user defined particle id code
  void SetParentTrack(AliTrack* t);       // Set pointer to the parent track
  AliTrack* GetParentTrack();             // Provide pointer to the parent track
  void SetProb(Double_t prob);            // Set the hypothesis probability for this track
  Float_t GetProb() const;                // Provide the hypothesis probability for this track
  void SetFitDetails(TObject* obj);       // Enter the object containing the fit details
  void SetFitDetails(TObject& obj) { SetFitDetails(&obj); }
  TObject* GetFitDetails();               // Provide pointer to the object containing the fit details
  void SetTimestamp(AliTimestamp& t);     // Set the track timestamp
  AliTimestamp* GetTimestamp();           // Provide the track timestamp
  void RemoveTimestamp();                 // Remove timestamp from this track
  Double_t GetDistance(AliPosition* p,Float_t scale=-1);   // Provide distance to position p
  Double_t GetDistance(AliPosition& p,Float_t scale=-1) { return GetDistance(&p,scale); }
  Double_t GetDistance(AliTrack* t,Float_t scale=-1);      // Provide distance to track t
  Double_t GetDistance(AliTrack& t,Float_t scale=-1) { return GetDistance(&t,scale); }
  void SetEscale(Float_t scale);          // Set the scale of the energy-momentum units of the track
  Float_t GetEscale() const;              // Provide the scale of the energy-momentum units of the track
 
 protected:
  void Init();               // Initialisation of pointers etc...
  Float_t fQ;                // The charge of the particle
  TObjArray* fDecays;        // The array of decay produced tracks
  TObjArray* fSignals;       // The array of related AliSignals
  TObjArray* fHypotheses;    // The array of track hypotheses
  AliPositionObj* fBegin;    // The begin-point of the track 
  AliPositionObj* fEnd;      // The end-point of the track 
  AliPositionObj* fRef;      // The reference-point of the track for the 3-momentum vector
  AliPositionObj* fImpactXY; // The (extrapolated) impact-point in the plane z=0
  AliPositionObj* fImpactXZ; // The (extrapolated) impact-point in the plane y=0
  AliPositionObj* fImpactYZ; // The (extrapolated) impact-point in the plane x=0
  Int_t fUserId;             // The user defined identifier
  AliPositionObj* fClosest;  // The (extrapolated) point of closest approach w.r.t some reference
  Int_t fCode;               // The user defined particle id code
  AliTrack* fParent;         // Pointer to the parent track
  Float_t fProb;             // Probability for this track as a hypothesis
  TObject* fFit;             // Object containing details of the fit
  AliTimestamp* fTstamp;     // The track timestamp
  Float_t fEscale;           // The scale of the energy-momentum units of the track

 private:
  void Dumps(AliTrack* t,Int_t n,TString f,TString u); // Recursively print all decay levels
 
 ClassDef(AliTrack,20) // Handling of the attributes of a reconstructed particle track.
};
#endif
