#ifndef ALITRDQABUILDREFERENCE_H
#define ALITRDQABUILDREFERENCE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaBuildReference.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Build the reference histograms                                        //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TFile;
class TH1D;

class AliTRDqaBuildReference: public TObject {

 public:

  AliTRDqaBuildReference() {}          // ctor
  AliTRDqaBuildReference(const AliTRDqaBuildReference& qadm):TObject(qadm) {}   
  AliTRDqaBuildReference& operator = (const AliTRDqaBuildReference& qadm);
  virtual ~AliTRDqaBuildReference() {;} // dtor
  
  virtual void Copy(TObject &/*qadm*/) const;
  void     BuildRefHistos(TFile *file) const;
  Double_t CalculateQuality(TH1D* /*measured*/, TH1D* /*reference*/) const {return 1.;}

 private:

  ClassDef(AliTRDqaBuildReference,1)   // Creates the TRD QA data

};
#endif
