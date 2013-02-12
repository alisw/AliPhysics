#ifndef ALIITSUMODULE_H
#define ALIITSUMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSUModule.h 53509 2011-12-10 18:55:52Z masera $ */
///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Class AliITSUModule                                            //
//  is a superclass for AliITSUModuleSSD, SPD and SDD.             //
//  The main function of modules is to simulate DIGITS from          //
//  GEANT HITS and produce POINTS from DIGITS                        //
//  It also make fast simulation without use of DIGITS               //
//                                                                   //
//  created by: A.Boucham, W.Peryt, S.Radomski, P.Skowronski         //
//              R.Barbera, B. Batynia, B. Nilsen                     //
//  ver. 1.0    CERN, 16.09.1999                                     //
//  modified for upgrade: ruben.shahoyan@cern.ch                     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>
#include "AliITSUHit.h"
class AliITSUGeomTGeo;

class AliITSUModule: public TObject {

 public:
  AliITSUModule();             // default constructor
  AliITSUModule(Int_t index, AliITSUGeomTGeo* tg);
  virtual ~AliITSUModule();
  UInt_t     GetIndex()                                 const {return GetUniqueID();}
  void       SetIndex(UInt_t ind)                             {return SetUniqueID(ind);}
  Int_t      GetNHits()                                 const {return fHitsM->GetEntriesFast();}
  TObjArray *GetHits()                                  const {return fHitsM;}
  AliITSUHit *GetHit(Int_t i)                         const {return (AliITSUHit*)fHitsM->UncheckedAt(i);}
  void       AddHit(AliITSUHit *hit)                        {fHitsM->AddLast(hit);}
  void       Clear(Option_t* opt=0);
  //
  Bool_t   MedianHitG(AliITSUHit *h1,AliITSUHit *h2,Float_t &x,Float_t &y,Float_t &z);
  void     MedianHitG(Int_t index, Float_t hitx1,Float_t hity1,Float_t hitz1,Float_t hitx2,Float_t hity2,Float_t hitz2, Float_t &xMg,Float_t &yMg, Float_t &zMg);
  Bool_t   MedianHitL(AliITSUHit *h1,AliITSUHit *h2,Float_t &x,Float_t &y,Float_t &z) const;
  void     MedianHitL(Int_t,AliITSUHit *,AliITSUHit *,Float_t &,Float_t &, Float_t &){};
  Double_t PathLength(const AliITSUHit *itsHit1,const AliITSUHit *itsHit2);
  void     MedianHit(Int_t index, Float_t xg,Float_t yg,Float_t zg,Int_t status,Float_t &xMg, Float_t &yMg, Float_t &zMg,Int_t &flag);
  void     PathLength(Float_t x,Float_t y,Float_t z,Int_t status,Int_t &nseg,Float_t &x1,Float_t &y1,Float_t &z1,Float_t &dx1,Float_t &dy1, Float_t &dz1,Int_t &flag) const;
  Bool_t   LineSegmentL(Int_t hindex,Double_t &a,Double_t &b,Double_t &c,Double_t &d,Double_t &e,Double_t &f,Double_t &de);
  Bool_t   LineSegmentL(Int_t hindex,Double_t &a,Double_t &b,Double_t &c,Double_t &d,Double_t &e,Double_t &f,Double_t &de, Double_t &tof, Int_t &track);
  //
  Bool_t   LineSegmentG(Int_t hindex,Double_t &a,Double_t &b,Double_t &c,Double_t &d,Double_t &e,Double_t &f,Double_t &de);
  Bool_t   LineSegmentG(Int_t hindex,Double_t &a,Double_t &b,Double_t &c,Double_t &d,Double_t &e,Double_t &f,Double_t &de, Double_t &tof, Int_t &track);
  //
 protected:
    AliITSUModule(const AliITSUModule &source); 
    AliITSUModule& operator=(const AliITSUModule &source); 
    TObjArray            *fHitsM;     // Pointer to list of hits on this module
    AliITSUGeomTGeo    *fGeomTG;    // pointed to geometry
    //
    ClassDef(AliITSUModule,1) // Copy the hits into a more useful order
};

inline void AliITSUModule::Clear(Option_t *) {fHitsM->Clear();}

#endif



