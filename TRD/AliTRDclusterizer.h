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

#include "AliTRDrawStream.h"
#include "AliTRDgeometry.h"

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
class AliTRDCalSingleChamberStatus;
class AliTRDrecoParam;
class AliTRDCalOnlineGainTableROC;

class AliTRDclusterizer : public TNamed 
{

 public:

  // steering flags
  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum{
    kTrOwner = BIT(14)  //  toggle online tracklets ownership
    ,kClOwner= BIT(15)  //  toggle cluster ownership
    ,kLabels = BIT(16)  //  toggle MC labels for clusters
    ,kSkipTrafo = BIT(17)  //  skip the coordinate transformation of clusters?
    ,kLUT    = BIT(18)  //  using look up table for cluster's r-phi position
    ,kGAUS   = BIT(19)  //  using gauss approx. for cluster's r-phi position
    ,knewDM  = BIT(20)  //  was the digitsmanger created by raw2clusters?
    ,kTracksOwner = BIT(21) //  toggle GTU track ownership
  };

  struct MaxStruct
  {
    Int_t       row;           // row of the current cluster candidate
    Int_t       col;           // col of the current cluster candidate
    Int_t       time;          // time -"-
    Short_t     signals[3];    // signals of the maximum pad and it's twon neigbours
    UChar_t     padStatus;     // padStatus of the current cluster candidate
    Bool_t      fivePad;       // is this cluster candidate part of a 5 pad cluster (two overlaping clusters)?
  MaxStruct():row(-1),col(0),time(0),padStatus(0),fivePad(kFALSE)
      {memset(signals, 0, 3*sizeof(Short_t));}
    MaxStruct &operator=(const MaxStruct &a)
    {memcpy(this, &a, sizeof(MaxStruct)); return *this;}
  };
  
  AliTRDclusterizer(const AliTRDReconstructor *const rec = 0x0);
  AliTRDclusterizer(const Text_t* name, const Text_t* title, const AliTRDReconstructor *const rec = 0x0);
  AliTRDclusterizer(const AliTRDclusterizer &c);
  virtual         ~AliTRDclusterizer();
  AliTRDclusterizer &operator=(const AliTRDclusterizer &c);

  void     Copy(TObject &c) const;

  Bool_t   Open(const Char_t *name, Int_t nEvent = 0);
  Bool_t   OpenInput(Int_t nEvent = 0);
  Bool_t   OpenOutput();
  Bool_t   OpenOutput(TTree *const clusterTree);

  Bool_t   ReadDigits();
  Bool_t   ReadDigits(AliRawReader *rawReader);
  Bool_t   ReadDigits(TTree *digitsTree);

  Bool_t   ReadTracklets();
  Bool_t   ReadTracks();

  Bool_t   WriteClusters(Int_t det);
  void     ResetRecPoints();
  virtual TClonesArray    *RecPoints();
  virtual TClonesArray    *TrackletsArray(const TString &trkltype = "");
  virtual TClonesArray    *TracksArray();

  Bool_t   Raw2Clusters(AliRawReader *rawReader);
  Bool_t   Raw2ClustersChamber(AliRawReader *rawReader);

  Bool_t   MakeClusters();
  Bool_t   MakeClusters(Int_t det);

  Bool_t   AddLabels();
  Bool_t   SetUseLabels(const Bool_t kset) { SetBit(kLabels, kset); return TestBit(kLabels);  } // should we assign labels to clusters
  void     SetRawVersion(const Int_t iver) { fRawVersion = iver; } // set the expected raw data version
  void             SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
  static UChar_t   GetStatus(Short_t &signal);
  Int_t            GetAddedClusters() const {return fNoOfClusters;}
  Int_t            GetNTimeBins() const {return fTimeTotal;}

  Bool_t           IsClustersOwner() const {return TestBit(kClOwner);}
  virtual void     SetClustersOwner(Bool_t own=kTRUE) {SetBit(kClOwner, own); if(!own) {fRecPoints = 0x0; fNoOfClusters=0;} }
  void             SetTrackletsOwner(Bool_t own=kTRUE) {SetBit(kTrOwner, own); if(!own) {fTracklets = 0x0; } }
  void             SetTracksOwner(Bool_t own=kTRUE) {SetBit(kTracksOwner, own); if(!own) {fTracks = 0x0; } }
  void             SetSkipTransform(Bool_t b=kTRUE) {SetBit(kSkipTrafo, b); }

  UInt_t   GetTriggerFlags(const Int_t sector) const { return fTrgFlags[sector]; }

protected:

  void             ApplyTCTM(Short_t *const arr, const Int_t nTime, const Int_t nexp);
  void             DeConvExp (Short_t *const arr, const Int_t nTime, const Int_t nexp);
  void             ConvExp(Short_t *const arr, const Int_t nTime);
  void             TailCancelation(const AliTRDrecoParam* const recoParam);

  Float_t          Unfold(Double_t eps, Int_t layer, const Double_t *const padSignal) const;
  
  void             SetPadStatus(const UChar_t status, UChar_t &encoding) const;
  UChar_t          GetPadStatus(UChar_t encoding) const;
  Int_t            GetCorruption(UChar_t encoding) const;

  Bool_t           IsMaximum(const MaxStruct &Max, UChar_t &padStatus, Short_t *const Signals);       //for const correctness reasons not const parameters are given separately
  Bool_t           FivePadCluster(MaxStruct &ThisMax, MaxStruct &NeighbourMax);
  void             CreateCluster(const MaxStruct &Max); 

  virtual void     AddClusterToArray(AliTRDcluster* cluster);

private:
  inline void      CalcAdditionalInfo(const MaxStruct &Max, Short_t *const signals, Int_t &nPadCount);

protected:
  const AliTRDReconstructor *fReconstructor; //! reconstructor
  AliRunLoader        *fRunLoader;           //! Run Loader
  TTree               *fClusterTree;         //! Tree with the cluster
  TClonesArray        *fRecPoints;           //! Array of clusters
  TClonesArray        *fTracklets;           //! Array of online tracklets
  TClonesArray        *fTracks;              //! Array of GTU tracks

  TTree               *fTrackletTree;        //! Tree for tracklets

  AliTRDdigitsManager *fDigitsManager;       //! TRD digits manager

  UInt_t              **fTrackletContainer;  //! tracklet container
					     // legacy code to avoid breakint AliHLTTRDClusterizer
					     // but it's useless

  Int_t                fRawVersion;          //  Expected raw version of the data - default is 2

  AliTRDtransform     *fTransform;           //! Transforms the reconstructed space points

  AliTRDarrayADC      *fDigits;               // Array holding the digits
  AliTRDSignalIndex   *fIndexes;              // Array holding the indexes to the digits
  Short_t              fMaxThresh;            // Threshold value for the maximum
  Short_t              fMaxThreshTest;        // Threshold value for the maximum to test agains
  Short_t              fSigThresh;            // Threshold value for the digit signal
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
  AliTRDCalSingleChamberStatus *fCalPadStatusROC; // Calibration object with the pad status
  AliTRDCalOnlineGainTableROC *fCalOnGainROC; // Calibration table of online gain factors
  Int_t                fClusterROC;           // The index to the first cluster of a given ROC
  Int_t                firstClusterROC;       // The number of cluster in a given ROC
  Int_t                fNoOfClusters;         // Number of Clusters already processed and still owned by the clusterizer
  Int_t                fBaseline;             // Baseline of the ADC values
  AliTRDrawStream     *fRawStream;            // Raw data streamer
  UInt_t               fTrgFlags[AliTRDgeometry::kNsector]; // trigger flags

  ClassDef(AliTRDclusterizer,6)               //  TRD clusterfinder

};

#endif
