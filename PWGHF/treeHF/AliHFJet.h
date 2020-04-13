#ifndef ALIHFJET_H
#define ALIHFJET_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJet
// \helper class to handle jet objects
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include "TObject.h"
class AliHFJet : public TObject
{
  public:
  


  AliHFJet();
  AliHFJet(const AliHFJet &source);
  virtual ~AliHFJet();
  void Reset();

  Float_t GetID() {return fID;}
  Float_t GetHFMeson() {return fHFMeson;}
  Float_t GetPt() {return fPt;}
  Float_t GetEta() {return fEta;}
  Float_t GetPhi() {return fPhi;}
  Float_t GetDeltaEta() {return fDeltaEta;}
  Float_t GetDeltaPhi() {return fDeltaPhi;}
  Float_t GetDeltaR() {return fDeltaR;}
  Float_t GetN() {return fN;}
  Float_t GetZg() {return fZg;}
  Float_t GetRg() {return fRg;}
  Float_t GetNsd() {return fNsd;}
  Float_t GetPt_mother() {return fPt_mother;}
  Float_t Getk0() {return fk0;}
  Float_t Getk1() {return fk1;}
  Float_t Getk2() {return fk2;}
  Float_t GetkT() {return fkT;}

    

  Float_t fID;        //unique (in event) jet ID
  Float_t fHFMeson;   //determines if the jet contains the HF candidtae or particle
  Float_t fPt;        //jet pT
  Float_t fEta;       //jet pseudorapidity
  Float_t fPhi;       //jet phi
  Float_t fDeltaEta;  //pseudorapidity difference of jet axis and HF candidiate or particle
  Float_t fDeltaPhi;  //phi difference of jet axis and HF candidiate or particle
  Float_t fDeltaR;    //pseudorapidity-phi distnace of jet axis and HF candidiate or particle
  Float_t fN;         //number of jet constituents
  Float_t fZg;        //soft dropped splitting momentum fraction
  Float_t fRg;        //soft dropped splitting angle
  Float_t fNsd;       //number of splittings passing soft drop
  Float_t fPt_mother; //total pT going into soft drop splitting
  Float_t fk0;         //dynamical grooming kappa with alpha=0
  Float_t fk1;         //dynamical grooming kappa with alpha=1
  Float_t fk2;         //dynamical grooming kappa with alpha=2
  Float_t fkT;         //Splitting witht the largest kT (following hardest branch)


  /// \cond CLASSIMP
  ClassDef(AliHFJet,1); ///
  /// \endcond
};
#endif
