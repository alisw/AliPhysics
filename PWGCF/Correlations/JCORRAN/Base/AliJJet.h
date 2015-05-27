/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

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

#ifndef ALIJJET_H
#define ALIJJET_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TRefArray.h>
#include  "AliJConst.h"
#include "AliJBaseTrack.h"

using namespace std;

class AliJJet : public AliJBaseTrack {
public:
  AliJJet();
  AliJJet(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge); // constructor
  AliJJet(const AliJJet& a);
  AliJJet(const TLorentzVector & a);
  virtual ~AliJJet(){
      ; 
  }    //destructor
  AliJJet& operator=(const AliJJet& trk);
  
  void SetArea(double a){ fArea = a; }
  Double_t GetArea() const{ return fArea; }
  Double_t Area() const{ return fArea; }
  void AddConstituent(TObject* t){ fConstituents.Add(t); }
  void AddConstituentRef(TObject* t){ fConstituentsRef.Add(t); }
  TRefArray* GetConstituentsRef(){ return &fConstituentsRef; }
  TObjArray* GetConstituents(){ GenConstituentsFromRef();return &fConstituents; }
  int GetNConstituentsRef(){ return fConstituentsRef.GetEntriesFast(); }
  int GetNConstituents(){ GenConstituentsFromRef();return fConstituents.GetEntriesFast(); }
  AliJBaseTrack * GetConstituent(int i) { GenConstituentsFromRef();return (AliJBaseTrack*)fConstituents.At(i); }
  AliJBaseTrack * GetConstituentRef(int i) const{ return (AliJBaseTrack*)fConstituentsRef.At(i); }
  void ReSum();
  void ReSum2(){ SetE( E()>fE2?E():fE2 ); }
  int LeadingParticleId(){ return fLeadingTrackId; }
  double LeadingParticlePt(){ return fLeadingTrackPt; }
  double LeadingParticleE(){ return fLeadingTrackE; }
  virtual void    Clear(Option_t* = ""){
      fConstituents.Clear();
      fConstituentsRef.Clear();
  }
  void SetConstituentsOwner(){
      fConstituents.SetOwner(kTRUE);
  }
  double E2nd(){ return fE2; }
  double E2nd2(){ return fE2*fE2; }

  void GenConstituentsFromRef( int force = 0){
      if( force == 1 || fConstituents.GetEntriesFast() < 1 ){
          fConstituents.Clear();
          for( int i=0;i<GetNConstituentsRef();i++ ){
              AddConstituent( GetConstituentRef(i) );
          }
      }
  }


private:
  int      fLeadingTrackId;     //! id of leading track in constituents
  double   fLeadingTrackPt;
  double   fLeadingTrackE;
  double   fNConstituent;
  double   fE2;
  Double_t fArea;              // Area of the jet
  TObjArray fConstituents;     //! Constituent tracks of the jets
  TRefArray fConstituentsRef;     // Constituent tracks of the jets


  ClassDef(AliJJet,3)
};
#endif
