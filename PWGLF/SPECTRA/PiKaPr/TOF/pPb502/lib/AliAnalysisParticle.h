#ifndef ALIANALYSISPARTICLE_H
#define ALIANALYSISPARTICLE_H

#include "TObject.h"
#include "TLorentzVector.h"

class TParticle;

class AliAnalysisParticle :
public TObject
{

 public:

  AliAnalysisParticle(); // default constructor
  AliAnalysisParticle(const AliAnalysisParticle &source); // copy constructor
  AliAnalysisParticle &operator=(const AliAnalysisParticle &source); // operator=
  virtual ~AliAnalysisParticle(); // default destructor

  Int_t GetLabel() const {return fLabel;} // get label
  Float_t GetPt() const {return fPt;}; // get pt
  Float_t GetEta() const {return fEta;}; // get eta
  Float_t GetPhi() const {return fPhi;}; // get phi
  Int_t GetPdgCode() const {return fPdgCode;}; // get PDG code

  Double_t GetY() const; // get Y
  Float_t GetSign() const; // get sign
  Int_t GetPID() const; // get MC PID
  Double_t GetMass() const; // get mass

  void Reset(); // reset
  void Update(TParticle *particle, Int_t label); // update
  
 private:

  Int_t fLabel; // label
  Float_t fPt; // pt
  Float_t fEta; // eta
  Float_t fPhi; // phi
  Int_t fPdgCode; // PDG code

  /*** tools ***/
  static TLorentzVector fgLorentzVector;

  ClassDef(AliAnalysisParticle, 1);
};

#endif /* ALIANALYSISPARTICLE_H */
