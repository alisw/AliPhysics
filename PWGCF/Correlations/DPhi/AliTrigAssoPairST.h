/*
 *****************************************************************************************/
#ifndef ALITRIGASSOPAIRST_H
#define ALITRIGASSOPAIRST_H

#include "AliVParticle.h"
#include "TString.h"
#include "AliLog.h"


class AliTrigAssoPairST : public AliVParticle {

 public :
   AliTrigAssoPairST(Short_t charge, Float_t eta, Float_t phi, Float_t pt, Float_t pt2,
                       Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate,
                       Double_t multiplicity, Double_t deta_pairs, Double_t dphi_pairs);
   virtual ~AliTrigAssoPairST();

 virtual Double_t Px() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Py() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Pz() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Pt() const { return fpT; }
  virtual Double_t P() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t PxPyPz(Double_t[3]) const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Xv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Yv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Zv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t XvYvZv(Double_t[3]) const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t OneOverPt() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Phi() const { return fPhi; }
  virtual Double_t Theta() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t E() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t M() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Eta() const { return fEta; }
  virtual Double_t Y() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Short_t Charge() const { return fCharge; }
  virtual Int_t GetLabel() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Int_t PdgCode() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual const Double_t *PID() const {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Short_t WhichCandidate() const { return fCandidate; }
  virtual Int_t GetID() const { return fID; }
  virtual Int_t GetIDFirstDaughter() const { return fID1; }
  virtual Int_t GetIDSecondDaughter() const { return fID2; }
  virtual Double_t Multiplicity() const { return fMultiplicity; }
  virtual Double_t Getdeta_pairs() const { return fdeta_pairs; }
  virtual Double_t Getdphi_pairs() const { return fdphi_pairs; }
  virtual Double_t Pt_Asso() const { return fpT_Asso; }

 private:
  // 

  Short_t fCharge;    // Charge
  Float_t fEta;       // Eta
  Float_t fPhi;       // Phi
  Float_t fpT;        // pT
  Float_t fpT_Asso;        // Associate pT
  Float_t fdeta_pairs;        // dEta_Pairs
  Float_t fdphi_pairs;        // dPhi_Pairs
  Int_t fID;          // ID
  Short_t fCandidate; // 1-pi,2-K,3-p
  Double_t fMultiplicity;
  Int_t fID1;
  Int_t fID2;
  ClassDef(AliTrigAssoPairST, 1);

};

#endif
