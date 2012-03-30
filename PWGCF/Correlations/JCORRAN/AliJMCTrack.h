// $Id: AliJMCTrack.h,v 1.3 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJMCTrack.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.3 $
  \date $Date: 2008/05/08 15:19:52 $
  */
////////////////////////////////////////////////////

#ifndef ALIJMCTRACK_H
#define ALIJMCTRACK_H

// include ROOT lib
#ifndef ROOT_TObject
#include <TObject.h>
#endif
#include <TDatabasePDG.h>

// include Inharitance
#include "AliJBaseTrack.h"

class AliJMCTrack : public AliJBaseTrack {

 public:
  enum { kPrimary=AliJBaseTrack::kNFlag, kPHOS, kEMCAL, kTPC, kNFlag };//For ALICE
  enum { kFinal=AliJBaseTrack::kNFlag };// for MC
  //usage : this->SetFlag( kPrimary, kTRUE );
  //usage : this->IsTrue( kFinal );

  AliJMCTrack();	    //default constructor
  AliJMCTrack(const AliJMCTrack& a);	//copy constructor

  ~AliJMCTrack(){;}		//destructor

  Int_t  GetPdgCode()    const {return fPdgCode;}
  Int_t  GetMother  (Int_t i) const {return fMother[i];}
  Int_t  GetDaughter(Int_t i) const {return fDaughter[i];}
  Double32_t GetVx()    const{return fVx;}
  Double32_t  GetVy()    const{return fVy;}
  Double32_t  GetVz()    const{return fVz;}

  const TParticlePDG& GetPDGData() const ;
  void SetPdgCode(Int_t icode);// Set PDG and E,charge
  void SetMother  (int m0, int m1){ fMother[0] = m0;fMother[1]=m1; }
  void SetDaughter(int d0, int d1){ fDaughter[0]=d0;fDaughter[1]=d1; }
  void SetProductionVertex(Double_t vx, Double_t vy, Double_t vz)
  {fVx=vx; fVy=vy; fVz=vz;}


  AliJMCTrack& operator=(const AliJMCTrack& trk);


  //Extra
  Bool_t IsFinal() const { return IsTrue( kFinal ); }
  void  SetIsFinal(Bool_t t){ SetFlag( kFinal, t );}
  //TODO
  Bool_t IsHadron() const;
  Bool_t IsCharged() const { return GetCharge(); }
  Bool_t IsParton()  const {return ( fPdgCode < -7 && fPdgCode < 7 && fPdgCode !=0 );}

 private:

  Short_t         fPdgCode;              // PDG code of the particle
  Short_t         fMother[2];            // Index of the mother particles
  Short_t         fDaughter[2];          // Indices of the daughter particles

  Double32_t      fVx;                   //[0.,0.,12] x of production vertex
  Double32_t      fVy;                   //[0.,0.,12] y of production vertex
  Double32_t      fVz;                   //[0.,0.,12] z of production vertex

  ClassDef(AliJMCTrack,2)
};

#endif
