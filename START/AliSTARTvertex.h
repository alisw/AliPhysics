#ifndef STARTVERTEX_H
#define STARTVERTEX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */



#include <AliESD.h>

//___________________________________________
class AliSTARTvertex   : public TObject {


////////////////////////////////////////////////////////////////////////
 public:
  AliSTARTvertex():TObject(),fZposition(0) {}
  virtual ~AliSTARTvertex() {}

  void Reconstruct(AliRunLoader* runLoader, AliESD *pESD);

  Float_t GetVertex() const {return fZposition;}
  void SetVertex(Float_t zPosition) {fZposition=zPosition;}

 private:
    Float_t fZposition;        // Z position of vertex (mm)

    ClassDef(AliSTARTvertex,1)  //Reconstructive vertex (Header) object 
};
#endif
