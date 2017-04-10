//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliFilteredTrack.h
/// @author Matthias Richter
/// @date   2012-04-02
/// @brief  AliVParticle with minimal information
///

#ifndef ALIFILTEREDTRACK_H
#define ALIFILTEREDTRACK_H

#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "TMath.h"

class TMCParticle;

/**
 * @class AliFilteredTrack
 * An AliVTrack/AliVParticle implementation with minimal persistent members.
 * - pt, phi, theta
 * - (x,y,z of initial point) to be implemented
 *
 * Other properties are derived from minimal parameter set.
 * Automatically calculated in the Streamer when loading an object.
 */
class AliFilteredTrack : public AliVTrack {
 public:
  /// default constructor
  AliFilteredTrack();
  /// constructor
  //AliFilteredTrack(float pt, float phi, float theta, int charge);
  /// copy constructor
  AliFilteredTrack(const AliFilteredTrack& other);
  /// copy constructor
  AliFilteredTrack(const AliVParticle& other);
  /// copy constructor
  AliFilteredTrack(const TMCParticle& other);
  /// copy constructor
  AliFilteredTrack(const AliAODTrack& other);
  /// assignment operator
  AliFilteredTrack& operator=(const AliFilteredTrack& other);
  /// assignment operator
  AliFilteredTrack& operator=(const AliVParticle& other);
  /// destructor
  virtual ~AliFilteredTrack();

  enum {
    kX = 0,
    kY,
    kZ,
    kNofDim
  };

  enum {
    kSignBit  		= BIT(14),   // unset -> 'plus', set -> 'minus'
    kGlobalHybrid  	= BIT(15),   // set -> selected by GlobalHybrid
    kBIT4 		= BIT(16),   // set -> selected by FilterBit4
    kBIT5 		= BIT(17),   // set -> selected by FilterBit5
    kBIT6 		= BIT(18),    // set -> selected by FilterBit6
    kMC			= BIT(19),    //set -> is MCparticle
    kLeading		= BIT(20)     //set -> is the Leading pT particle in the event
  };

  // kinematics
  virtual Double_t Px() const {return fP[kX];}
  virtual Double_t Py() const {return fP[kY];}
  virtual Double_t Pz() const {return fP[kZ];}
  Double_t GetPx() const {return fP[kX];}
  Double_t GetPy() const {return fP[kY];}
  Double_t GetPz() const {return fP[kZ];}
  virtual Double_t Pt() const {return fPt;}
  virtual Double_t P()  const {return fPtot;}
  virtual Bool_t   PxPyPz(Double_t p[3]) const {
    p[kX]=fP[kX]; p[kY]=fP[kY]; p[kZ]=fP[kZ]; return kTRUE;
  }

  virtual Double_t Xv() const {return -999.;}
  virtual Double_t Yv() const {return -999.;}
  virtual Double_t Zv() const {return -999.;}
  virtual Bool_t   XvYvZv(Double_t x[3]) const {
    x[kX]=Xv(); x[kY]=Yv(); x[kZ]=Zv(); 
    return kTRUE;
  }

  virtual Double_t OneOverPt()  const {return fOneOverPt;}
  virtual Double_t Phi()        const {return fPhi;}
  virtual Double_t Theta()      const {return fTheta;}


  virtual Double_t E()          const {return -999.;}
  virtual Double_t M()          const {return -999.;}
  
  virtual Double_t Eta()        const {return fEta;}
  virtual Double_t Y()          const {return -999.;}
  
  virtual Short_t Charge()      const {
    return TestBit(kSignBit)?-1:1;
  }
  virtual Double_t GetSign()    const {return Charge();}
  virtual Int_t   GetLabel()    const {return -999;}

  virtual Int_t PdgCode() const {return -999;}
  virtual const Double_t* PID() const {return NULL;}
  virtual Int_t GetID() const {return -1;}
  virtual UChar_t GetITSClusterMap() const {return 0;}
  virtual ULong_t GetStatus() const {return 0;}
  virtual Bool_t GetXYZ(Double_t*) const {return kFALSE;}
  virtual Bool_t GetCovarianceXYZPxPyPz(Double_t*) const {return kFALSE;}
  virtual Bool_t PropagateToDCA(const AliVVertex*, Double_t, Double_t, Double_t*, Double_t*) {return kFALSE;}

  /// overloaded from TObject: cleanup
  virtual void Clear(Option_t * option ="");
  /// overloaded from TObject: print info
  virtual void Print(Option_t *option="") const;

  template<typename T> void Set(const T p[kNofDim]);
  void  Set(const AliVParticle& track);
  void  Set(const AliAODTrack& track);
  void  SetCharge(const AliVParticle& track);
  void  SetEff(float eff){feff = eff;}
  float GetEff(){return feff;}
  bool  IsGlobalHybrid(){return TestBit(kGlobalHybrid);}
  bool  IsBIT4(){return TestBit(kBIT4);}
  bool  IsBIT5(){return TestBit(kBIT5);}
  bool  IsBIT6(){return TestBit(kBIT6);}
  bool  IsLeading(){return TestBit(kLeading);}
  bool 	IsMC(){return TestBit(kMC);}
  void  SetGlobal(){SetBit(kGlobalHybrid);}
  void  SetBIT4(){SetBit(kBIT4);}
  void  SetBIT5(){SetBit(kBIT5);}
  void  SetBIT6(){SetBit(kBIT6);}
  void  SetLeading(){SetBit(kLeading);}
  void 	SetMC(bool isset){
    if(isset){
      SetBit(kMC);
    }
    else{
      ResetBit(kMC);
    }
  }
  void  SetAODFilterBits(const AliAODTrack * t);
  void  Calculate(bool bCalculateMommentumComponents=true);

 protected:
 private:
   
  // persistent members stored for the class
  float fPt;        // transverse momentum
  float fPhi;       // phi
  float fTheta;     // theta
  // all other members are transient and calculated from the
  // momentum vector
  float feff;	    //! efficiency
  float fP[3];      //! momentum vector
  float fOneOverPt; //! one over Pt
  float fPtot;      //! momentum
  float fEta;       //! eta

  ClassDef(AliFilteredTrack, 1)
};

template<typename T>
void AliFilteredTrack::Set(const T p[kNofDim])
{
  // set momentum vector and calculate internal values
  fP[kX]=p[kX];
  fP[kY]=p[kY];
  fP[kZ]=p[kZ];
  fPt = TMath::Sqrt(p[kX]*p[kX] + p[kY]*p[kY]);
  fPhi = TMath::ATan2(p[kY],p[kX]);
  if (fPhi<0) fPhi += 2*TMath::Pi();
  fTheta =  TMath::ATan2(fPt,p[kZ]);

  Calculate(false);
}
#endif
/*
#if defined(__MAKECINT__) 
// define the class to have a custom streamer
#pragma link C++ class AliFilteredTrack-;
#endif*/

