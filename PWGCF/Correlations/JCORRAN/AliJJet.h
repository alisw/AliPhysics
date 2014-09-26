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
#include  "AliJConst.h"
#include "AliJBaseTrack.h"

using namespace std;

class AliJJet : public AliJBaseTrack {
public:
  AliJJet();
  AliJJet(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge); // constructor
  AliJJet(const AliJJet& a);
  AliJJet(const TLorentzVector & a);
  virtual ~AliJJet(){; }    //destructor
  AliJJet& operator=(const AliJJet& trk);
  
  void SetArea(double a){ fArea = a; }
  Double_t GetArea() const{ return fArea; }
  Double_t Area() const{ return fArea; }
  void AddConstituent(TObject* t){ fConstituents.Add(t); }
  TObjArray* GetConstituents(){ return &fConstituents; }
  AliJBaseTrack * GetConstituent(int i) const{ return (AliJBaseTrack*)fConstituents[i]; }
  
private:
  Double_t fArea;              // Area of the jet
  TObjArray fConstituents;     // Constituent tracks of the jets
  
  ClassDef(AliJJet,1)
};
#endif
