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
  Float_t GetLeadingPt() {return fLeadingPt;}
  Float_t GetDeltaEta() {return fDeltaEta;}
  Float_t GetDeltaPhi() {return fDeltaPhi;}
  Float_t GetDeltaR() {return fDeltaR;}
  Float_t GetN() {return fN;}
  Float_t GetZ() {return fZ;}
  Float_t GetAngularityk1B1() {return fAngularityk1B1;}
  Float_t GetpTDispersion() {return fpTDispersion;}
  Float_t GetChargek03() {return fChargek03;}
  Float_t GetChargek05() {return fChargek05;}
  Float_t GetChargek07() {return fChargek07;}
  Float_t GetZg() {return fZg;}
  Float_t GetRg() {return fRg;}
  Float_t GetNsd() {return fNsd;}
  Float_t GetPt_splitting() {return fPt_splitting;}
  Float_t Getk0() {return fk0;}
  Float_t GetZk0() {return fZk0;}
  Float_t GetRk0() {return fRk0;}
  Float_t Getk1() {return fk1;}
  Float_t GetZk1() {return fZk1;}
  Float_t GetRk1() {return fRk1;}
  Float_t Getk2() {return fk2;}
  Float_t GetZk2() {return fZk2;}
  Float_t GetRk2() {return fRk2;}
  Float_t GetkT() {return fkT;}
  Float_t GetZkT() {return fZkT;}
  Float_t GetRkT() {return fRkT;}

    

  Float_t fID;        //unique (in event) jet ID
  Float_t fHFMeson;   //determines if the jet contains the HF candidtae or particle
  Float_t fPt;        //jet pT
  Float_t fEta;       //jet pseudorapidity
  Float_t fPhi;       //jet phi
  Float_t fLeadingPt; //leading track pT in jet
  Float_t fDeltaEta;  //pseudorapidity difference of jet axis and HF candidiate or particle
  Float_t fDeltaPhi;  //phi difference of jet axis and HF candidiate or particle
  Float_t fDeltaR;    //pseudorapidity-phi distnace of jet axis and HF candidiate or particle
  Float_t fN;         //number of jet constituents
  Float_t fZ;         //fragmentation function
  Float_t fAngularityk1B1;     //jet angularity with kappa = 1 and beta =1
  Float_t fpTDispersion;       //jet momentum dispersion
  Float_t fChargek03;          //jet charge with kappa = 0.3
  Float_t fChargek05;          //jet charge with kappa = 0.5
  Float_t fChargek07;          //jet charge with kappa = 0.7
  Float_t fZg;        //soft dropped splitting momentum fraction
  Float_t fRg;        //soft dropped splitting angle
  Float_t fNsd;       //number of splittings passing soft drop
  Float_t fPt_splitting; //total pT going into soft drop splitting
  Float_t fk0;         //dynamical grooming kappa with alpha=0
  Float_t fZk0;        //dynamical grooming with alpha=0 splitting momentum fraction
  Float_t fRk0;        //dynamical grooming with alpha=0 splitting angle
  Float_t fk1;         //dynamical grooming kappa with alpha=1
  Float_t fZk1;        //dynamical grooming with alpha=1 splitting momentum fraction
  Float_t fRk1;        //dynamical grooming with alpha=1 splitting angle
  Float_t fk2;         //dynamical grooming kappa with alpha=2
  Float_t fZk2;        //dynamical grooming with alpha=2 splitting momentum fraction
  Float_t fRk2;        //dynamical grooming with alpha=2 splitting angle
  Float_t fkT;         //Splitting witht the largest kT (following hardest branch)
  Float_t fZkT;        //Splitting witht the largest kT (following hardest branch) splitting momentum fraction
  Float_t fRkT;        //Splitting witht the largest kT (following hardest branch) splitting angle


  /// \cond CLASSIMP
  ClassDef(AliHFJet,1); ///
  /// \endcond
};
#endif
