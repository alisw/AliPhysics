#ifndef ALITRDCLUSTERIZER_H
#define ALITRDCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD cluster finder                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TFile;
class TTree;
class TClonesArray;

class AliRunLoader;
class AliRawReader;

class AliTRD;
class AliTRDcluster;

class AliTRDarrayADC;
class AliTRDarraySignal;
class AliTRDdigitsManager;
class AliTRDSignalIndex;
class AliTRDtransform;
class AliTRDCalROC;
class AliTRDReconstructor;

class AliTRDclusterizer : public TNamed 
{

 public:

  // steering flags
  enum{
    kOwner = BIT(14)
  };

  struct MaxStruct
  {
    Int_t       Row;
    Int_t       Col;
    Int_t       Time;
    UChar_t     padStatus;
    Float_t     Signals[3];
    MaxStruct():Row(0),Col(0),Time(0),padStatus(0)
      {}
    MaxStruct &operator=(const MaxStruct &a)
      {Row=a.Row; Col=a.Col; Time=a.Time; padStatus=a.padStatus;
       memcpy(Signals, a.Signals, 3*sizeof(Signals[0])); return *this;}
  };
  
  AliTRDclusterizer(const AliTRDReconstructor *const rec = 0x0);
  AliTRDclusterizer(const Text_t* name, const Text_t* title, const AliTRDReconstructor *const rec = 0x0);
  AliTRDclusterizer(const AliTRDclusterizer &c);
  virtual         ~AliTRDclusterizer();
  AliTRDclusterizer &operator=(const AliTRDclusterizer &c);

  virtual void     Copy(TObject &c) const;

  virtual Bool_t   Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t   OpenInput(Int_t nEvent = 0);
  virtual Bool_t   OpenOutput();
  virtual Bool_t   OpenOutput(TTree *clusterTree);

  virtual Bool_t   ReadDigits();
  virtual Bool_t   ReadDigits(AliRawReader *rawReader);
  virtual Bool_t   ReadDigits(TTree *digitsTree);

  virtual Bool_t   WriteClusters(Int_t det);
          void     ResetRecPoints();
  TClonesArray       *RecPoints();
          Bool_t   WriteTracklets(Int_t det);

  virtual Bool_t   Raw2Clusters(AliRawReader *rawReader);
  virtual Bool_t   Raw2ClustersChamber(AliRawReader *rawReader);

  virtual Bool_t   MakeClusters();
  virtual Bool_t   MakeClusters(Int_t det);

  virtual Bool_t   AddLabels(const Int_t idet, const Int_t firstClusterROC, const Int_t nClusterROC);
  virtual Bool_t   SetAddLabels(const Bool_t kset) { fAddLabels = kset; return fAddLabels;  } // should we assign labels to clusters
  virtual void     SetRawVersion(const Int_t iver) { fRawVersion = iver; } // set the expected raw data version
  void             SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
  static UChar_t   GetStatus(Short_t &signal);

  Bool_t           IsClustersOwner() const {return TestBit(kOwner);}
  void             SetClustersOwner(Bool_t own=kTRUE) {SetBit(kOwner, own); if(!own) fRecPoints = 0x0;}

 protected:

  void             DeConvExp (const Double_t *const source, Double_t *const target
			     ,const Int_t nTimeTotal, const Int_t nexp);
  void             TailCancelation();

  virtual Float_t  Unfold(Double_t eps, Int_t layer, Double_t *padSignal) const;
          Double_t GetCOG(Double_t signal[5]) const; 
  void             FillLUT();
          Double_t LUTposition(Int_t ilayer, Double_t ampL, Double_t ampC, Double_t ampR) const;
  
  void             SetPadStatus(const UChar_t status, UChar_t &encoding);
  UChar_t          GetPadStatus(UChar_t encoding) const;
  Int_t            GetCorruption(UChar_t encoding) const;

  Bool_t           IsMaximum(const MaxStruct &Max, UChar_t &padStatus, Float_t *const Signals);       //for const correctness reasons not const parameters are given separately
  Bool_t           IsFivePadCluster(const MaxStruct &ThisMax, const MaxStruct &NeighbourMax, Float_t &ratio); // ''
  void             CreateCluster(const MaxStruct &Max); 

  const AliTRDReconstructor *fReconstructor;       //! reconstructor
  AliRunLoader        *fRunLoader;           //! Run Loader
  TTree               *fClusterTree;         //! Tree with the cluster
  TClonesArray        *fRecPoints;           //! Array of clusters

  TTree               *fTrackletTree;        //! Tree for tracklets

  AliTRDdigitsManager *fDigitsManager;       //! TRD digits manager

  UInt_t              **fTrackletContainer;  //! tracklet container

  Bool_t               fAddLabels;           //  Should clusters have MC labels?
  Int_t                fRawVersion;          //  Expected raw version of the data - default is 2

  AliTRDtransform     *fTransform;           //! Transforms the reconstructed space points

  Int_t                fLUTbin;              //  Number of bins of the LUT
  Double_t            *fLUT;                 //! The lookup table

  AliTRDarrayADC      *fDigitsIn;
  AliTRDSignalIndex   *fIndexes;
  Float_t              fADCthresh;            // ADC thresholds: There is no ADC threshold anymore, and simParam should not be used in clusterizer. KO
  Float_t              fMaxThresh;            // Threshold value for the maximum
  Float_t              fSigThresh;            // Threshold value for the digit signal
  Float_t              fMinMaxCutSigma;       // Threshold value for the maximum (cut noise)
  Float_t              fMinLeftRightCutSigma; // Threshold value for the sum pad (cut noise)
  Int_t                fLayer;                // Current layer of the detector
  Int_t                fDet;                  // Current detecor
  UShort_t             fVolid;                // Volume ID
  Int_t                fColMax;               // Number of Colums in one detector
  Int_t                fTimeTotal;            // Number of time bins
  AliTRDCalROC        *fCalGainFactorROC;     // Calibration object with pad wise values for the gain factors
  Float_t              fCalGainFactorDetValue;// Calibration value for chamber wise noise
  AliTRDCalROC        *fCalNoiseROC;          // Calibration object with pad wise values for the noise
  Float_t              fCalNoiseDetValue;     // Calibration value for chamber wise noise
  AliTRDarraySignal   *fDigitsOut;
  Int_t                fClusterROC;           // The index to the first cluster of a given ROC
  Int_t                firstClusterROC;       // The number of cluster in a given ROC

  ClassDef(AliTRDclusterizer,6)              //  TRD clusterfinder

};

#endif
