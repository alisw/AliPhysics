#ifndef ALITRDCLUSTERIZERHLT_H
#define ALITRDCLUSTERIZERHLT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD cluster finder                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDclusterizer.h"

class AliRawReader;
class AliRawReaderMemory;
class AliTRDrawData;
class AliTRDfeeParam;

class AliTRDclusterizerHLT : public AliTRDclusterizer 
{
 public:

  AliTRDclusterizerHLT();
  AliTRDclusterizerHLT(const Text_t* name, const Text_t* title, AliTRDReconstructor* reconstructor);
  AliTRDclusterizerHLT(const Text_t* name, const Text_t* title);
  AliTRDclusterizerHLT(const AliTRDclusterizerHLT &c);
  virtual             ~AliTRDclusterizerHLT();
  AliTRDclusterizerHLT &operator=(const AliTRDclusterizerHLT &c);

  virtual void     Copy(TObject &c) const;
  //virtual Bool_t   MakeClusters();
  virtual Bool_t   ReadDigits(AliRawReaderMemory* rawReader);
  using AliTRDclusterizer::ReadDigits;
  //virtual Bool_t   TreeClusters(Int_t idet);
  virtual Bool_t   InitClusterTree();
  virtual Bool_t   InsertClusters(TClonesArray *tobjarr, Int_t idet);
  virtual Int_t    CountClusters();
  virtual Int_t    GetNclusters();

  virtual Bool_t   ResetTree();
 
  TTree *          GetClusterTree() {return fClusterTree;}
  virtual Bool_t   IsTreeOwner() const {return fTreeCreatedHere;}

  virtual void     SetRawVersion(Int_t iver);

 protected:
  virtual Bool_t   ReadDigits() const {return kFALSE;} //this method not to be used on HLT
  virtual Bool_t   ReadDigits(AliRawReader* /*rawReader*/) const {/*if (!rawReader) */ return kFALSE;} //this method not to be used on HLT

  Bool_t            fTreeCreatedHere; //flag indicating that AliTRDclusterizerHLT has created the cluster tree

 private:

  Int_t fNclusters; //number of clusters found - updated by ::GetNclusters()
  AliTRDrawData  *fRawDataSource; //! pointer to the TRD raw data stream
  AliTRDfeeParam *fFeeParams; //! fee parameters - this is singleton so watch out! with delete

  ClassDef(AliTRDclusterizerHLT,1)           //  TRD-Cluster finder, slow simulator

};

#endif
