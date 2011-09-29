#ifndef ALIAODCONVERSIONMOTHER_H
#define ALIAODCONVERSIONMOTHER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class reconstructing the mother particle of conversion gammas
//---------------------------------------------
////////////////////////////////////////////////

//Author Daniel Lohner (Daniel.Lohner@cern.ch)

#include "TLorentzVector.h"
#include "AliAODConversionParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliKFConversionMother.h"
#include "AliKFParticle.h"

class AliAODConversionMother : public AliAODConversionParticle{

 public: 

  //Default Constructor
     AliAODConversionMother();

  // Constructor for ESD to AOD Conversion
     AliAODConversionMother(AliKFConversionMother *kf);
 
  //Constructor Decay Mother Particle
     AliAODConversionMother(AliAODConversionPhoton *y1,AliAODConversionPhoton *y2);

  //Destructor
     virtual ~AliAODConversionMother();

  ///Set the Chi2 of reconstructed conversion gamma
     void SetChi2(Float_t chi2) {fChi2 = chi2;}

  //Get the Chi2 of particle
     Float_t Chi2() const {return fChi2;}

     ///Set track or MC labels
     void SetLabel1(Int_t label){fLabel[0] = label;}
     void SetLabel2(Int_t label){fLabel[1] = label;}
     void SetLabels(Int_t label1, Int_t label2){fLabel[0] = label1; fLabel[1] = label2;}

     Int_t GetLabel(Int_t i) const {return fLabel[i];}
     Int_t GetLabel1() const {return fLabel[0];}
     Int_t GetLabel2() const {return fLabel[1];}

     Double_t GetOpeningAngle() const { return fOpeningAngle;}

     Double_t GetAlpha() const { return fAlpha;}


private:
    Int_t fLabel[2]; // Labels of the decay photons
    Float_t fChi2; // Chi sq of reconstructed mother
    Double_t fOpeningAngle;
    Double_t fAlpha;

    ClassDef(AliAODConversionMother,1)
};

#endif
