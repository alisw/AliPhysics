#ifndef ALITOFCLUSTERFINDERV1_H
#define ALITOFCLUSTERFINDERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// AliTOFClusterFinderV1 Class
// Task: Transform digits/raw data to TOF Clusters, to fill TOF RecPoints
// and feed TOF tracking 

#include "TObject.h"
#include "TNamed.h"

#include "AliTOFGeometry.h"
#include "AliTOFRawStream.h"

class TClonesArray;
class TTree;

class AliRunLoader;
class AliRawReader;

class AliTOFcluster;
class AliTOFcalib;
class AliTOFDigitMap;
class AliTOFRecoParam;

class  AliTOFselectedDigit : public TObject { 
 public:
  AliTOFselectedDigit() :
    fTDC(0.),fADC(0.),fTOT(0.),fWeight(0.),fIndex(-1) {
    for (Int_t ii=0; ii<5; ii++) fDetectorIndex[ii]=-1;
    for (Int_t ii=0; ii<3; ii++) fTrackLabel[ii]=-1;
  };
  AliTOFselectedDigit(Int_t * const ind, Double_t h1, Double_t h2, Double_t h3, Double_t h4, Int_t idx, Int_t * const l):
    TObject(),
    fTDC(h1),fADC(h2),fTOT(h3),fWeight(h4),fIndex(idx) {
    for (Int_t ii=0; ii<5; ii++) fDetectorIndex[ii]=ind[ii];
    for (Int_t ii=0; ii<3; ii++) fTrackLabel[ii]=l[ii];
  };
  AliTOFselectedDigit(const AliTOFselectedDigit & source) :
    TObject(source),
    fTDC(source.fTDC),fADC(source.fADC),fTOT(source.fTOT),fWeight(source.fWeight),fIndex(source.fIndex)
    {
      for (Int_t ii=0; ii<5; ii++) fDetectorIndex[ii]=source.fDetectorIndex[ii];
      for (Int_t ii=0; ii<3; ii++) fTrackLabel[ii]=source.fTrackLabel[ii];
    };
  AliTOFselectedDigit & operator=(const AliTOFselectedDigit & source)
    { if (this == &source) return *this;
      TObject::operator=(source);
      for (Int_t ii=0; ii<5; ii++) fDetectorIndex[ii]=source.fDetectorIndex[ii];
      fTDC=source.fTDC;fADC=source.fADC;fTOT=source.fTOT;fWeight=source.fWeight;fIndex=source.fIndex;
      for (Int_t ii=0; ii<3; ii++) fTrackLabel[ii]=source.fTrackLabel[ii];
      return *this; };

  Double_t GetTDC()    const {return fTDC;} // digit TOF
  Double_t GetADC()    const {return fADC;} // digit ADC
  Double_t GetTOT()    const {return fTOT;} // digit TOT
  Double_t GetWeight() const {return fWeight;} // digit weight
  Int_t GetTrackLabel(Int_t n)    const {return fTrackLabel[n];} // Labels of tracks in digit
  Int_t GetDetectorIndex(Int_t n) const {return fDetectorIndex[n];} // Digit Detector Index n
  Int_t GetIndex()     const {return fIndex;} // Digit Index

 private: 

  Int_t fDetectorIndex[5]; //digit detector indices (sector, plate, strip, padX, padZ) 
  Double_t fTDC; //TDC count
  Double_t fADC; //ADC count
  Double_t fTOT; //TOT count
  Double_t fWeight; //weight
  Int_t fIndex; //index of the digit in the TOF digit tree
  Int_t fTrackLabel[3]; //track labels

}; 


class AliTOFClusterFinderV1 : public TNamed
{

  enum {kTofMaxCluster=77777}; //maximal number of the TOF clusters

 public:

  AliTOFClusterFinderV1(AliTOFcalib *calib);
  AliTOFClusterFinderV1(const AliTOFClusterFinderV1 &source); // copy constructor
  AliTOFClusterFinderV1& operator=(const AliTOFClusterFinderV1 &source); // ass. op.
  virtual ~AliTOFClusterFinderV1();

  void Digits2RecPoints(TTree* digitsTree, TTree* clusterTree);
  void Digits2RecPoints(AliRawReader *rawReader, TTree *clustersTree);
  void Raw2Digits(AliRawReader *rawReader, TTree* digitsTree); 

  AliTOFClusterFinderV1(AliRunLoader* runLoader, AliTOFcalib *calib);
  void Digits2RecPoints(Int_t iEvent);
  void Digits2RecPoints(Int_t ievt, AliRawReader *rawReader);
  void Raw2Digits(Int_t ievt, AliRawReader *rawReader);

  void FillRecPoint();
  void ResetRecpoint();
  void ResetDigits();
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;} // To set the verbose level
  void SetDecoderVersion(Int_t version){fDecoderVersion=version;} // To set the decoder version
  Bool_t GetDecoderVersion() const {return fDecoderVersion;} // To get the decoder version
  //UShort_t GetClusterVolIndex(Int_t *ind) const; //Volume Id getter
  void GetClusterPars(Int_t *ind, Double_t *pos, Double_t *cov) const; //cluster par getter
  void GetClusterPars(/*Bool_t check,*/ Int_t counter, Int_t **ind, Double_t *weight,
		      Double_t *pos, Double_t *cov) const; //cluster par getter

  void FindOnePadClusterPerStrip(Int_t nSector, Int_t nPlate, Int_t nStrip);
  void FindClustersWithoutTOT(Int_t nSector, Int_t nPlate, Int_t nStrip);
  void FindClustersPerStrip(Int_t nSector, Int_t nPlate, Int_t nStrip, Int_t group);

  void FindClusters34(Int_t nSector, Int_t nPlate, Int_t nStrip);
  void FindClusters23(Int_t nSector, Int_t nPlate, Int_t nStrip);
  void FindClusters24(Int_t nSector, Int_t nPlate, Int_t nStrip);

  void  SetMaxDeltaTime(Int_t a)   {fMaxDeltaTime = a;}; // to set deltaTime [bin number]
  void  SetMaxDeltaTime(Float_t a) {fMaxDeltaTime = (Int_t)(a/AliTOFGeometry::TdcBinWidth());}; // to set deltaTime [ps]
  Int_t GetMaxDeltaTime()     const {return fMaxDeltaTime;};


  void SetCalibrateFlag(Bool_t dummy) {fCalibrateTOFtimes = dummy;};
  Bool_t GetCalibrateFlag() const {return fCalibrateTOFtimes;};

 protected:

  AliRunLoader *fRunLoader;   // Pointer to Run Loader
  AliTOFcluster *fTofClusters[kTofMaxCluster];  // pointers to the TOF clusters
  TClonesArray *fDigits;      // List of digits
  TClonesArray *fRecPoints;   // List of reconstructed points
  Int_t fNumberOfTofClusters; // Number of TOF Clusters
  Int_t fNumberOfTofDigits;   // Number of TOF Digits

 private:

  const AliTOFRecoParam* fkRecoParam; // pointer to TOF reconstruction parameters

  Int_t fMaxDeltaTime; // max time difference in between two tof
                       // measurements for two neighbouring pads

  Int_t InsertCluster(AliTOFcluster *tofCluster);    // Fills TofClusters Array
  Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
  Bool_t MakeSlewingCorrection(Int_t *detectorIndex, Int_t tofDigitToT, Int_t tofDigitTdc,
			       Int_t &tdcCorr);
  void TOFclusterError(/*Bool_t check,*/ Int_t counter, Int_t **ind, Double_t *weight,
		       Double_t ppos[], Double_t cov[]) const;

  void AverageCalculations(Int_t number, Float_t *interestingX,
			   Float_t *interestingY, Float_t *interestingZ,
			   Double_t *interestingTOF, Double_t *interestingTOT,
			   Double_t *interestingADC, Double_t *interestingWeight,
			   Int_t *parTOF, Double_t *posClus, Bool_t &check);

  Int_t fVerbose;           // Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  Bool_t fDecoderVersion;   // setting whether to use the new decoder version 
                            //  -true -> new version
                            //  -false ->old version  (default value!!)
  AliTOFcalib *fTOFcalib;        // pointer to the TOF calibration info
  AliTOFDigitMap* fTOFdigitMap;  // TOF digit map pointer
  AliTOFGeometry *fTOFGeometry;  // pointer to the TOF geometry
  TTree *fTOFdigits;             // pointer to the TOF digit tree

  AliTOFRawStream fTOFRawStream; // AliTOFRawStream variable

  Bool_t fCalibrateTOFtimes;     // used for check

  ClassDef(AliTOFClusterFinderV1,5) // To run TOF clustering
};
#endif

