#ifndef ALITPCCLUSTERERMI_H
#define ALITPCCLUSTERERMI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */



//-------------------------------------------------------
//                       TPC clusterer
//
//   Origin: Marian Ivanov  
//-------------------------------------------------------
#include <Rtypes.h>
#include <TObject.h>
#include <AliTPCRecoParam.h>
#define kMAXCLUSTER 2500

class TFile;
class TClonesArray;
class AliTPCParam;
class AliTPCRecoParam;
class AliTPCclusterMI;
class AliTPCClustersRow;
class AliRawReader;
class AliSimDigits;
class TTree;
class TTreeSRedirector;
class  AliRawEventHeaderBase;
class AliTPCCalROC;

class AliTPCclustererMI : public TObject{
public:
  AliTPCclustererMI(const AliTPCParam* par, const AliTPCRecoParam * recoParam = 0);
  AliTPCclustererMI(const AliTPCclustererMI &param); // copy constructor
  AliTPCclustererMI &operator = (const AliTPCclustererMI & param); //assignment
  virtual ~AliTPCclustererMI();
  virtual void Digits2Clusters();
  virtual void Digits2ClustersOld(AliRawReader* rawReader);
  virtual void Digits2Clusters(AliRawReader* rawReader);
  virtual void SetInput(TTree * tree);  // set input tree with digits
  virtual void SetOutput(TTree * tree); // set output tree with 
  virtual void FillRow();               // fill the output container - Tree or TObjArray
  TObjArray * GetOutputArray(){return fOutputArray;}
  TClonesArray * GetOutputClonesArray(){return fOutputClonesArray;}

  void StoreInClonesArray(Bool_t bOutput = kTRUE) {fBClonesArray = bOutput;} // store output clusters in TClonesArray

private:
  Bool_t IsMaximum(Float_t k, Int_t max, const Float_t *bins) const; 
  void MakeCluster2(Int_t k,Int_t max,Float_t *bins,UInt_t m,
   AliTPCclusterMI &c);  
  void MakeCluster(Int_t k,Int_t max,Float_t *bins,UInt_t m,
   AliTPCclusterMI &c); 
  Float_t  GetSigmaY2(Int_t iz);
  Float_t  GetSigmaZ2(Int_t iz);
  Float_t  FitMax(Float_t vmatrix[5][5], Float_t y, Float_t z, Float_t sigmay, Float_t sigmaz);
  void AddCluster(AliTPCclusterMI &c, Float_t *matrix, Int_t pos);  // add the cluster to the array
  void UnfoldCluster(Float_t * matrix[7], Float_t recmatrix[5][5], 
		     Float_t & meani, Float_t & meanj, Float_t & sum, Float_t &overlap );
  void FindClusters(AliTPCCalROC * noiseROC);
  Bool_t AcceptCluster(AliTPCclusterMI*c);
  Double_t  ProcesSignal(Float_t * signal, Int_t nchannels, Int_t id[3], Double_t &rms, Double_t &pedestalCalib);
  void ProcessSectorData();
  Int_t ReadHLTClusters();
  
  Float_t * fBins;       //!digits array
  Int_t   * fSigBins; //!digits array containg only timebins above threshold
  Int_t     fNSigBins;//!size of fSigBins
  Int_t fLoop;         //loop - cf in 2 loops
  Int_t fMaxBin;       //current ( for current sector)  maximal bin
  Int_t fMaxTime;      //current ( for current sector)  maximal time
  Int_t fMaxPad;       //current ( for current sector)  maximal pad
  Int_t fSector;      //!current sector
  Int_t fRow;         //!current row
  Float_t fSign;      //!current sign 
  Float_t fRx;        // current radius
  Float_t fPadWidth;  // the width of the pad
  Float_t fPadLength;  // the width of the pad
  Float_t fZWidth;     //the z bin width
  Bool_t  fPedSubtraction; // perform pedestal subtraction or not
  AliRawEventHeaderBase *fEventHeader; //! event header information
  UInt_t  fTimeStamp;   // Time Stamp
  UInt_t  fEventType;   // Event Type
  TTree * fInput;   //!input  tree with digits - object not owner
  TTree * fOutput;   //!output tree with digits - object not owner
  TObjArray *fOutputArray;     //! output TObjArray with pointers arrays of cluster
  TClonesArray *fOutputClonesArray; //! output TClonesArray with clusters
  AliTPCClustersRow * fRowCl;  //! current cluster row
  AliSimDigits * fRowDig;      //! current digits row
  const AliTPCParam * fParam;        //! tpc parameters
  Int_t fNcluster;             // number of clusters - for given row
  Int_t fNclusters;            // tot number of clusters
  TTreeSRedirector *fDebugStreamer;     //!debug streamer
  const AliTPCRecoParam  * fRecoParam;        //! reconstruction parameters
  Bool_t  fBDumpSignal; // dump signal flag
  Bool_t  fBClonesArray; // output clusters stored in TClonesArray 
  Bool_t  fBUseHLTClusters; // use HLT clusters instead of offline clusters

  // Non-persistent arrays

  Float_t** fAllBins; //! All sector bins
  Int_t** fAllSigBins;//! All signal bins in a sector
  Int_t*  fAllNSigBins;//! Number of signal bins in a sector


  ClassDef(AliTPCclustererMI,2)  // Time Projection Chamber digits
};

inline Bool_t AliTPCclustererMI::IsMaximum(Float_t q,Int_t max,const Float_t *bins) const {
  //is this a local maximum ?
  if (bins[-max] >= q) return kFALSE;
  if (bins[-1  ] >= q) return kFALSE; 
  if (bins[+max] > q) return kFALSE; 
  if (bins[+1  ] > q) return kFALSE; 
  if (bins[-max-1] >= q) return kFALSE;
  if (bins[+max-1] >= q) return kFALSE; 
  if (bins[+max+1] > q) return kFALSE; 
  if (bins[-max+1] >= q) return kFALSE;
  //
  //
  if (fRecoParam->GetClusterMaxRange(1)>0){  //local maxima in z on more than 1 time bin 
    if (bins[-2] > q) return kFALSE;
    if (bins[ 2] > q) return kFALSE;
  }
  if (fRecoParam->GetClusterMaxRange(0)>0){  //local maxima in y on more than pad 
    if (bins[-2*max] > q) return kFALSE;
    if (bins[ 2*max] > q) return kFALSE;
  }

  return kTRUE; 
}



//-----------------------------------------------------------------

#endif


