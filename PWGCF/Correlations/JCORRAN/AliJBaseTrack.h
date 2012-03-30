// $Id: AliJBaseTrack.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

///////////////////////////////////////////////////
/*
   \file AliJBaseTrack.h
   \brief
   \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
   \email: djkim@jyu.fi
   \version $Revision: 1.5 $
   \date $Date: 2008/05/08 15:19:52 $
   */
///////////////////////////////////////////////////

#ifndef ALIJBASETRACK_H
#define ALIJBASETRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include  "AliJConst.h"
#include  "TLorentzVector.h"

class AliJBaseTrack : public TLorentzVector {
 public:
  enum { kIsIsolated, kNFlag };
  AliJBaseTrack();
  AliJBaseTrack(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge); // constructor
  AliJBaseTrack(const AliJBaseTrack& a);
  AliJBaseTrack(const TLorentzVector & a);
  virtual ~AliJBaseTrack(){;}		//destructor

  float   GetTwoPiPhi() const {return Phi()>-kJPi/3 ? Phi() : kJTwoPi+Phi();} 
  TLorentzVector GetLorentzVector(){ return TLorentzVector( *this);}

  Int_t         GetID()           const { return fID;}
  Int_t         GetLabel()        const { return fLabel; }
  Short_t       GetParticleType() const { return fParticleType;}
  ULong_t       GetStatus()       const { return fStatus; }
  Short_t       GetCharge()       const { return fCharge; } 
  UInt_t        GetFlags()        const { return fFlags; }
  Bool_t        GetIsIsolated()   const { return IsTrue(kIsIsolated);}


  void SetID      (const int id){fID=id;}
  void SetLabel   (const Int_t label ){ fLabel=label; }
  void SetParticleType(const Short_t ptype){ fParticleType=ptype; }
  void SetStatus  (const ULong_t status){ fStatus=status; }
  void SetCharge  (const Char_t charge){ fCharge=charge; }
  void SetFlags   (const UInt_t bits ){ fFlags=bits; }
  void SetIsIsolated(Bool_t tf){ SetFlag( kIsIsolated, tf); }

  void Print(Option_t* option) const;

  // Handel BitsData
  Bool_t IsTrue(int i ) const { return TESTBIT(fFlags, i); }
  void SetFlag(int i, Bool_t t){ if(t){SETBIT(fFlags,i);}else{CLRBIT(fFlags, i);}}

  // Operators
  AliJBaseTrack& operator=(const AliJBaseTrack& trk);

 protected:
  Int_t         fID;            // Unique track ID
  Int_t         fLabel;         // Unique track label for MC-Data relation
  Short_t       fParticleType;  // ParticleType 
  Char_t        fCharge;        // track charge for real data
  ULong_t       fStatus;        // reconstruction status flags or MC status 
  UInt_t        fFlags;         // store series of any boolen value.

  ClassDef(AliJBaseTrack,1)
};

#endif

