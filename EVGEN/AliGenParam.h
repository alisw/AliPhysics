
#ifndef AliGenParam_H
#define AliGenParam_H
#include "AliGenerator.h"
#include "AliPythia.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TTree.h"
#include "TParticle.h"

//-------------------------------------------------------------
class AliGenParam : public AliGenerator
{
protected:
  Double_t (*fPtParaFunc)(Double_t*, Double_t*); //! Pointer to Pt parametrisation function
  Double_t (*fYParaFunc )(Double_t*, Double_t*); //! Pointer to Y parametrisation function
  Int_t    (*fIpParaFunc )();    //! Pointer to particle type parametrisation function
    TF1* fPtPara;
    TF1* fYPara;
    Param_t     fParam;
    Float_t     fdNdy0;
    Float_t     fYWgt;
    Float_t     fPtWgt;
    Weighting_t fAnalog;       //Flaf for anolog or pt-weighted generation
    Float_t     fBias;
    Int_t       fTrials;
    Decay_t     fForceDecay;
    Int_t       fCutOnChild;
    TArrayI   fChildSelect;
    AliPythia *fPythia;
 private:
    // check if particle is selected as child
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection goes here
    Bool_t KinematicSelection(TParticle *particle);
 public:
  AliGenParam();
  AliGenParam(Int_t npart, Param_t param);
//		   Double_t (*PtPara)(Double_t*, Double_t*),
//		   Double_t (*YPara )(Double_t*, Double_t*));
  virtual ~AliGenParam();
  virtual void Generate();
  virtual void Init();
  // select particle type
  virtual void SetParam(Param_t param=jpsi_p) {fParam=param;}
  // force decay type
  virtual void ForceDecay(Decay_t decay=dimuon) {fForceDecay=decay;}
  virtual void SetWeighting(Weighting_t flag=analog) {fAnalog=flag;}	
  virtual void SetCutOnChild(Int_t flag=0) {fCutOnChild=flag;}
  ClassDef(AliGenParam,1) // Generator using parameterised pt- and y-distribution
};
#endif










