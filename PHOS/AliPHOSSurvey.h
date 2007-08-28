#ifndef ALIPHOSSURVEY_H
#define ALIPHOSSURVEY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.3  2007/07/10 12:41:38  kharlov
 * Added a new class AliPHOSSurvet1 which read survey data from EDMS files
 *
 * Revision 1.2  2007/05/17 17:13:32  kharlov
 * Coding convensions satisfied (T.Pocheptsov)
 *
 * Revision 1.1  2007/04/19 15:47:20  kharlov
 * Add misalignment of strip units with AliPHOSSurvey class
 *
 */

#include <TObject.h>
#include <Rtypes.h>

class TClonesArray;
class TString;

class AliPHOSGeometry;

/*
  Objects of this class read txt file with survey (photogrammetry) data
  and convert the data into AliAlignObjParams of alignable PHOS volumes.
  It can be used as a base class, you need to override GetStripTransformation.
  AliPHOSSurvey inherits TObject only to use AliLog "functions".
*/

class AliPHOSSurvey : public TObject {
public:
  AliPHOSSurvey();
  AliPHOSSurvey(const TString &txtFileName);

  virtual ~AliPHOSSurvey();

  //Create AliAlignObjParams for strips.
  void CreateAliAlignObjParams(TClonesArray &array);
  //Create AliAlignObjParams with null shifts and rotations.
  void CreateNullObjects(TClonesArray &alObj, const AliPHOSGeometry *geom)const;

protected:

  struct AliPHOSStripDelta {
    Float_t fXShift; //x shift
    Float_t fYShift; //y shift
    Float_t fZShift; //z shift
    Float_t fPsi;    //psi
    Float_t fTheta;  //theta
    Float_t fPhi;    //phi
  };

  Int_t 	           fStrNum; // Number of strips.
  AliPHOSStripDelta *fStripData; // Strip unit transformation data

  void InitStripData(const Double_t *xReal, const Double_t *zReal);

private:
  //Calculate shifts and rotations for strip number stripIndex in a module moduleIndex.
  virtual AliPHOSStripDelta GetStripTransformation(Int_t stripIndex, Int_t moduleIndex)const;

  AliPHOSSurvey(const AliPHOSSurvey &);
  AliPHOSSurvey &operator = (const AliPHOSSurvey &);

  ClassDef(AliPHOSSurvey, 1) //Survey data reader
};

#endif
