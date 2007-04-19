#ifndef ALIPHOSSURVEY_H
#define ALIPHOSSURVEY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

#include <vector>

#include <TObject.h>
#include <Rtypes.h>

class TClonesArray;
class TString;

class AliPHOSGeometry;

/*
  Objects of this class read txt file with survey (photogrammetry) data
  and convert the data into AliAlignObjAngles of alignable PHOS volumes.
  It can be used as a base class, you need to override GetStripTransformation.
  AliPHOSSurvey inherits TObject only to use AliLog "functions".
*/

class AliPHOSSurvey : public TObject {
public:
  AliPHOSSurvey();
  AliPHOSSurvey(const TString &txtFileName);

  //Create AliAlignObjAngles for strips.
  void CreateAliAlignObjAngles(TClonesArray &array);
  //Create AliAlignObjAngles with null shifts and rotations.
  void CreateNullObjects(TClonesArray &, const AliPHOSGeometry *)const;

protected:
  struct Transformation_t {
    Float_t fXShift;
    Float_t fYShift;
    Float_t fZShift;
    Float_t fPsi;
    Float_t fTheta;
    Float_t fPhi;
  };

private:
  //Calculate shifts and rotations for strip number stripIndex in a module moduleIndex.
  virtual Transformation_t GetStripTransformation(Int_t stripIndex, Int_t moduleIndex)const;

  AliPHOSSurvey(const AliPHOSSurvey &);
  AliPHOSSurvey &operator = (const AliPHOSSurvey &);

private:
  std::vector<Transformation_t> fStripData; // Strip unit transformation data

  ClassDef(AliPHOSSurvey, 1) //Survey data reader
};

#endif
