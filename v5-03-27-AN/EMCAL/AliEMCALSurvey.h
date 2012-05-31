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

  enum SurveyDataType_t { kSurvey = 0, //use real survey parameters
			  kDummy = 1 //use dummy values for testing
  };

  AliEMCALSurvey();
  AliEMCALSurvey(const TString &txtFileName, const SurveyDataType_t dataType=kSurvey);

  virtual ~AliEMCALSurvey();

  //Create AliAlignObjParams for strips.
  void CreateAliAlignObjParams(TClonesArray &array);
  //Create AliAlignObjParams with null shifts and rotations.
  void CreateNullObjects(TClonesArray &alObj, const AliEMCALGeometry *geom)const;

  void  SetDataType(const SurveyDataType_t dataType) { fDataType = dataType; }
  Int_t GetDataType() const { return (Int_t)fDataType; }

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

  void InitSuperModuleData(const Double_t *xReal, const Double_t *yReal, const Double_t *zReal, 
			   const Double_t *psiReal, const Double_t *thetaReal, const Double_t *phiReal);
  void InitSuperModuleData(const TObjArray* surveypoints);

private:
  //Calculate shifts and rotations for supermodule.
  virtual AliEMCALSuperModuleDelta GetSuperModuleTransformation(Int_t smIndex) const;

  AliEMCALSurvey(const AliEMCALSurvey &);
  AliEMCALSurvey &operator = (const AliEMCALSurvey &);

  Int_t  fDataType; //! which date type (survey or dummy) to use

  ClassDef(AliEMCALSurvey, 2) //Survey data reader
};

#endif
