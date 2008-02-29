#ifndef ALIEMCALSURVEY_H
#define ALIEMCALSURVEY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <Rtypes.h>

class TClonesArray;
class TString;

class AliEMCALGeometry;

/*
  Objects of this class read txt file with survey data
  and convert the data into AliAlignObjParams of alignable EMCAL volumes.
  AliEMCALSurvey inherits TObject only to use AliLog "functions".
*/

class AliEMCALSurvey : public TObject {
public:
  AliEMCALSurvey();
  AliEMCALSurvey(const TString &txtFileName);

  virtual ~AliEMCALSurvey();

  //Create AliAlignObjParams for strips.
  void CreateAliAlignObjParams(TClonesArray &array);
  //Create AliAlignObjParams with null shifts and rotations.
  void CreateNullObjects(TClonesArray &alObj, const AliEMCALGeometry *geom)const;

protected:

  struct AliEMCALSuperModuleDelta {
    Float_t fXShift; //x shift
    Float_t fYShift; //y shift
    Float_t fZShift; //z shift
    Float_t fPsi;    //psi
    Float_t fTheta;  //theta
    Float_t fPhi;    //phi
  };

  Int_t 	  fNSuperModule; // Number of supermodules.
  AliEMCALSuperModuleDelta *fSuperModuleData; // Supermodule transformation data

  void InitSuperModuleData(const Double_t *xReal, const Double_t *yReal, const Double_t *zReal);

private:
  //Calculate shifts and rotations for supermodule.
  virtual AliEMCALSuperModuleDelta GetSuperModuleTransformation(Int_t smIndex) const;

  AliEMCALSurvey(const AliEMCALSurvey &);
  AliEMCALSurvey &operator = (const AliEMCALSurvey &);

  ClassDef(AliEMCALSurvey, 1) //Survey data reader
};

#endif
