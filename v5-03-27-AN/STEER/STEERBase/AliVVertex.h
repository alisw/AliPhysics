#ifndef AliVVertex_H
#define AliVVertex_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     base class for ESD and AOD vertices
//     Author: A. Dainese
//-------------------------------------------------------------------------

#include <TNamed.h>

class AliVVertex: public TNamed {

public:
  AliVVertex() { }
  virtual ~AliVVertex() { }
  AliVVertex(const AliVVertex& vVert); 
  AliVVertex& operator=(const AliVVertex& vVert);

  // vertex properties
  virtual void     GetXYZ(Double_t position[3]) const = 0;
  virtual Double_t GetX() const = 0;
  virtual Double_t GetY() const = 0;
  virtual Double_t GetZ() const = 0;
  virtual void     GetCovarianceMatrix(Double_t covmatrix[6]) const = 0;
  

  virtual Double_t GetChi2perNDF() const = 0;
  virtual Double_t GetChi2() const = 0;
  virtual Int_t    GetNDF() const = 0;

  virtual Int_t    GetNContributors() const = 0;
  virtual void     PrintIndices() const = 0;
  virtual void     Print(Option_t* option = "") const = 0;

  virtual void Clear(Option_t* option) {TNamed::Clear(option);}


  ClassDef(AliVVertex,1)  // base class for vertices
};

#endif
