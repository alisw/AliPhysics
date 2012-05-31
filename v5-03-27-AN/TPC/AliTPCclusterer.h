#ifndef ALITPCCLUSTERER_H
#define ALITPCCLUSTERER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                  The TPC cluster finder
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include <TObject.h>

class TTree;
class AliTPCParam;
class AliTPCcluster;

class AliTPCclusterer : public TObject {
public:
  AliTPCclusterer():TObject(),fPar(0){};
  AliTPCclusterer(const AliTPCParam *par):TObject(), fPar(par){};
  AliTPCclusterer(const AliTPCclusterer &param); // copy constructor
  AliTPCclusterer &operator = (const AliTPCclusterer & param);
  Int_t Digits2Clusters(TTree *dig, TTree *clu);

private:
  class AliBin {
  public:
    UShort_t GetQ()    const   {return fQ;}
    UInt_t   GetMask() const   {return fMask;}
    void     SetQ(UShort_t q)  {fQ=q;}
    void     SetMask(UInt_t m) {fMask=m;}
  private:
    UShort_t fQ;  //signal
    UInt_t fMask; //peak mask
  };

private:
  static Bool_t IsMaximum(Int_t k, Int_t max, const AliBin *bins); 
  static void FindPeaks(Int_t k,Int_t m,AliBin*b,Int_t*idx,UInt_t*msk,Int_t&n);
  static void MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m);
  static void MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,
   AliTPCcluster &c);

  const AliTPCParam *fPar;      //! pointer to the TPC parameters

  ClassDef(AliTPCclusterer,1)  // the TPC cluster finder
};


inline Bool_t AliTPCclusterer::IsMaximum(Int_t k,Int_t max,const AliBin *bins){
  //is this a local maximum ?
  UShort_t q=bins[k].GetQ();
  if (q==1023) return kFALSE;
  if (bins[k-max].GetQ() > q) return kFALSE;
  if (bins[k-1  ].GetQ() > q) return kFALSE; 
  if (bins[k+max].GetQ() > q) return kFALSE; 
  if (bins[k+1  ].GetQ() > q) return kFALSE; 
  if (bins[k-max-1].GetQ() > q) return kFALSE;
  if (bins[k+max-1].GetQ() > q) return kFALSE; 
  if (bins[k+max+1].GetQ() > q) return kFALSE; 
  if (bins[k-max+1].GetQ() > q) return kFALSE;
  return kTRUE; 
}

//-----------------------------------------------------------------

#endif


