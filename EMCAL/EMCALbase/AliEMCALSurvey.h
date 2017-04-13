#ifndef ALIEMCALSURVEY_H
#define ALIEMCALSURVEY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <Rtypes.h>

class TClonesArray;
class TString;

class AliEMCALGeometry;

////////////////////////////////////////////////////////////////////////////
///
/// \class  AliEMCALSurvey
/// \ingroup EMCALbase
/// \brief Read survey data and create alignement
///
/// Objects of this class read txt file with survey data
/// and convert the data into AliAlignObjParams of alignable EMCAL volumes.
/// AliEMCALSurvey inherits TObject only to use AliLog "functions".
///
/// Dummy functions originally written before EMCAL installation and
/// survey are kept for backward compatibility, but now they are not
/// used.
///
/// Surveyed points on the EMCAL support rails were used with the CATIA
/// 3D graphics program to determine the positions of the bottom
/// corners of the active area for each supermodule.  These numbers are
/// read in from file and converted to position of the center and roll,
/// pitch, yaw angles of each installed SM.
///
/// \author J.L. Klay, Cal Poly
/// \author Adapted for DCAL by M.L. Wang CCNU & Subatech Oct-19-2012
///
////////////////////////////////////////////////////////////////////////////

class AliEMCALSurvey : public TObject 
{

public:
  
  enum SurveyDataType_t 
  { 
    kSurvey = 0, ///< use real survey parameters
    kDummy  = 1  ///< use dummy values for testing
  };

  AliEMCALSurvey();
  AliEMCALSurvey(const TString &txtFileName, const SurveyDataType_t dataType=kSurvey);

  virtual ~AliEMCALSurvey();

  // Create AliAlignObjParams for strips.
  void CreateAliAlignObjParams(TClonesArray &array);
  
  // Create AliAlignObjParams with null shifts and rotations.
  void CreateNullObjects(TClonesArray &alObj, const AliEMCALGeometry *geom) const;

  void  SetDataType(const SurveyDataType_t dataType) { fDataType = dataType    ; }
  Int_t GetDataType()                          const { return (Int_t)fDataType ; }

protected:

  struct AliEMCALSuperModuleDelta {
    Float_t fXShift; //x shift
    Float_t fYShift; //y shift
    Float_t fZShift; //z shift
    Float_t fPsi;    //psi
    Float_t fTheta;  //theta
    Float_t fPhi;    //phi
  };

  Int_t                     fNSuperModule;    ///< Number of supermodules.
  AliEMCALSuperModuleDelta *fSuperModuleData; ///< Supermodule transformation data

  void InitSuperModuleData(const Double_t *xReal  , const Double_t *yReal    , const Double_t *zReal, 
                           const Double_t *psiReal, const Double_t *thetaReal, const Double_t *phiReal);
  void InitSuperModuleData(const TObjArray* surveypoints);

private:
  
  // Calculate shifts and rotations for supermodule.
  virtual AliEMCALSuperModuleDelta GetSuperModuleTransformation(Int_t smIndex) const;

  AliEMCALSurvey             (const AliEMCALSurvey &);
  AliEMCALSurvey &operator = (const AliEMCALSurvey &);

  Int_t  fDataType; //!<! which date type (survey or dummy) to use

  /// \cond CLASSIMP
  ClassDef(AliEMCALSurvey, 2) ;
  /// \endcond

};

#endif // ALIEMCALSURVEY_H
