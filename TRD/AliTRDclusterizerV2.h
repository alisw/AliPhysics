#ifndef ALITRDCLUSTERIZERV2_H
#define ALITRDCLUSTERIZERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD cluster finder                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDclusterizer.h"

class AliTRDdataArrayI;
class AliTRDdataArrayF;
class AliTRDdigitsManager;
class AliTRDCalROC;
class AliRawReader;
class AliTRDSignalIndex;
class AliTRDtransform;

class AliTRDclusterizerV2 : public AliTRDclusterizer {

 public:

  AliTRDclusterizerV2();
  AliTRDclusterizerV2(const Text_t* name, const Text_t* title);
  AliTRDclusterizerV2(const AliTRDclusterizerV2 &c);
  virtual             ~AliTRDclusterizerV2();
  AliTRDclusterizerV2 &operator=(const AliTRDclusterizerV2 &c);

  virtual void     Copy(TObject &c) const;
  virtual Bool_t   Raw2Clusters(AliRawReader *rawReader);
  virtual Bool_t   Raw2ClustersChamber(AliRawReader *rawReader);
  virtual Bool_t   MakeClusters();
  virtual Bool_t   MakeClusters(Int_t det);
  virtual Bool_t   ReadDigits();
  virtual Bool_t   ReadDigits(AliRawReader *rawReader);
  virtual Bool_t   ReadDigits(TTree *digitsTree);

  virtual Bool_t   AddLabels(Int_t idet, Int_t firstClusterROC, Int_t nClusterROC);
  virtual Bool_t   SetAddLabels(Bool_t kset) { fAddLabels = kset; 
                                               return fAddLabels;  } // should we assign labels to clusters
  virtual void     SetRawVersion(Int_t iver) { fRawVersion = iver; } // set the expected raw data version

 protected:

          void     DeConvExp(Double_t *source, Double_t *target
                           , Int_t nTimeTotal, Int_t nexp);
	  void     TailCancelation(AliTRDdataArrayI *digitsIn
                                 , AliTRDdataArrayF *digitsOut 
                                 , AliTRDSignalIndex *indexesIn
			         , AliTRDSignalIndex *indexesOut
			         , Int_t nTimeTotal
		                 , Float_t ADCthreshold
		                 , AliTRDCalROC *calGainFactorROC
		                 , Float_t calGainFactorDetValue);
  virtual Double_t Unfold(Double_t eps, Int_t plane, Double_t *padSignal);
          Double_t GetCOG(Double_t signal[5]); 

  virtual void     ResetHelperIndexes(AliTRDSignalIndex *indexesIn);

  AliTRDdigitsManager *fDigitsManager;       //! TRD digits manager

  Bool_t               fAddLabels;           //  Should clusters have MC labels?
  Int_t                fRawVersion;          //  Expected raw version of the data - default is 2

  AliTRDSignalIndex   *fIndexesOut;          //! Helper indexes for clusterization
  AliTRDSignalIndex   *fIndexesMaxima;       //! Helper indexes for clusterization

  AliTRDtransform     *fTransform;           //! Transforms the reconstructed space points

  ClassDef(AliTRDclusterizerV2,3)            //  TRD clusterfinder

};

#endif
