
#ifndef AliGenParam_H
#define AliGenParam_H
#include "AliGenerator.h"
#include "AliPythia.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TTree.h"
#include "TMCParticle.h"

//-------------------------------------------------------------
// Generators specific to MUON Arm

// Generate heavy mesons - J/Psi, Upsilon, Phi

class AliGenParam : public AliGenerator
{
protected:
    Double_t (*fPtParaFunc)(Double_t*, Double_t*);
    Double_t (*fYParaFunc )(Double_t*, Double_t*);
    Int_t    (*fIpParaFunc )();    
    TF1* fPtPara;
    TF1* fYPara;
    Int_t fIpart;
    Float_t fdNdy0;
    Float_t fYWgt;
    Float_t fPtWgt;
    Float_t fBias;
    Int_t   fTrials;
    Decay_t fForceDecay;
    TArrayI   fChildSelect;
    AliPythia *fPythia;
 private:
    // check if particle is selected as child
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection goes here
    Bool_t KinematicSelection(TMCParticle *particle);
 public:
  AliGenParam();
  AliGenParam(Int_t npart, Int_t ipart);
//		   Double_t (*PtPara)(Double_t*, Double_t*),
//		   Double_t (*YPara )(Double_t*, Double_t*));
  virtual ~AliGenParam();
  virtual void Generate();
  virtual void Init();
  // select particle type
  virtual void SetPart(Int_t part=443) {fIpart=part;}
  // force decay type
  virtual void ForceDecay(Decay_t decay=dimuon) {fForceDecay=decay;}
  ClassDef(AliGenParam,1)
};
#endif










