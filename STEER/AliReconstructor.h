#ifndef ALIRECONSTRUCTOR_H
#define ALIRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// base class for reconstruction algorithm
// Derived classes should implement a default constructor and
// the virtual methods
//

#include <TObject.h>
#include <TString.h>

class TTree;
class AliRunLoader;
class AliRawReader;
class AliVertexer;
class AliTracker;
class AliESD;


class AliReconstructor: public TObject {
public:
  AliReconstructor(): TObject(), fOption() {};
  virtual ~AliReconstructor() {};

  virtual void         Init(AliRunLoader* /*runLoader*/) {};

  virtual Bool_t       HasDigitConversion() const {return kFALSE;};
  virtual void         ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;

  virtual Bool_t       HasLocalReconstruction() const {return kFALSE;};
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         Reconstruct(AliRunLoader* runLoader, 
				   AliRawReader* rawReader) const;

  virtual AliVertexer* CreateVertexer(AliRunLoader* /*runLoader*/) const 
    {return NULL;}
  virtual AliTracker*  CreateTracker(AliRunLoader* /*runLoader*/) const 
    {return NULL;}

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESD* esd) const;
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESD* esd) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd) const;

  virtual const char*  GetDetectorName() const;

  void                 SetOption(Option_t* option) {fOption = option;};
  virtual Option_t*    GetOption() const {return fOption.Data();};

private:
  TString              fOption;   //! option for reconstruction

  ClassDef(AliReconstructor, 0)   // base class for reconstruction algorithms
};

#endif
