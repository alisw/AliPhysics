#ifndef ALIANALYSISPIDCASCADEPARTICLE_H
#define ALIANALYSISPIDCASCADEPARTICLE_H

#include "TObject.h"
#include "TLorentzVector.h"

class TParticle;

class AliAnalysisPIDCascadeParticle :
public TObject
{

 public:

  AliAnalysisPIDCascadeParticle(); // default constructor
  AliAnalysisPIDCascadeParticle(const AliAnalysisPIDCascadeParticle &source); // copy constructor
  AliAnalysisPIDCascadeParticle &operator=(const AliAnalysisPIDCascadeParticle &source); // operator=
  virtual ~AliAnalysisPIDCascadeParticle(); // default destructor

  Int_t GetLabel() const {return fLabel;} // get label
  Float_t GetPt() const {return fPt;}; // get pt
  Float_t GetEta() const {return fEta;}; // get eta
  Float_t GetPhi() const {return fPhi;}; // get phi
  Int_t GetPdgCode() const {return fPdgCode;}; // get PDG code
  Int_t GetMotherPdgCode() const {return fMotherPdgCode; }; //get mother PDG code
  Double_t GetY() const; // get Y
  Int_t GetSign() const; // get sign
  Int_t GetPID() const; // get MC PID
  Double_t GetMass() const; // get mass

  void Reset(); // reset
  void Update(TParticle *particle, Int_t label, Int_t MotherPDG); // update

 private:

  Int_t fLabel; // label
  Float_t fPt; // pt
  Float_t fEta; // eta
  Float_t fPhi; // phi
  Int_t fPdgCode; // PDG code
  Int_t fMotherPdgCode; //Mother PDG code

  /*** tools ***/
  static TLorentzVector fgLorentzVector;

  ClassDef(AliAnalysisPIDCascadeParticle, 2);
};

#endif /* ALIANALYSISPIDCASCADEPARTICLE_H */
