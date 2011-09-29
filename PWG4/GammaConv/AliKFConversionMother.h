#ifndef ALIKFCONVERSIONMOTHER_H
#define ALIKFCONVERSIONMOTHER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class containing the aod information from conversions
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---

#include "TMath.h"
#include "AliKFParticle.h"
#include "AliKFConversionPhoton.h"

class AliKFConversionMother : public AliKFParticle {

 public: 

  //Constructors
  AliKFConversionMother();    
  //AliKFConversionMother(AliKFParticle &kfparticle);
 // AliKFConversionMother(const AliKFParticle &d1,const AliKFParticle &d2);
  AliKFConversionMother(const AliKFConversionPhoton &d1,const AliKFConversionPhoton &d2);


  //Copy Constructor
  AliKFConversionMother(const AliKFConversionMother & g);           
  //assignment operator
  AliKFConversionMother & operator = (const AliKFConversionMother & g);

  //Destructor
  virtual ~AliKFConversionMother() {;}

  ///Set track or MC labels
  void SetLabel1(Int_t label){fLabel[0] = label;}
  void SetLabel2(Int_t label){fLabel[1] = label;}
  void SetGammaLabels(Int_t label1, Int_t label2){fLabel[0] = label1; fLabel[1] = label2;}

  Int_t GetGammaLabel(Int_t i) const {return fLabel[i];}

  Double_t GetOpeningAngle(){return fOpeningAngle;}
  Double_t GetAlpha(){return fAlpha;}
  Double_t GetRapidity();

  Double_t M(){return GetMass();}

  Double_t Phi() const;

 private:

    Int_t fLabel[2]; // Labels of two decay gammas
    Double_t fOpeningAngle; // of decay gammas
    Double_t fAlpha; // of the meson

  ClassDef(AliKFConversionMother,1)
};


#endif



