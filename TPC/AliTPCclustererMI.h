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
#define kMAXCLUSTER 2500

class TFile;
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
  virtual void Digits2Clusters(AliRawReader* rawReader);
  virtual void SetOldRCUFormat(Bool_t rcuFormat = kFALSE)
    { fIsOldRCUFormat = rcuFormat; };
  virtual void SetInput(TTree * tree);  // set input tree with digits    
  virtual void SetOutput(TTree * tree); //set output tree with 
private:
  Bool_t IsMaximum(Float_t k, Int_t max, const Float_t *bins) const; 
  void MakeCluster2(Int_t k,Int_t max,Float_t *bins,UInt_t m,
   AliTPCclusterMI &c);  
  void MakeCluster(Int_t k,Int_t max,Float_t *bins,UInt_t m,
   AliTPCclusterMI &c); 
  Float_t  GetSigmaY2(Int_t iz);
  Float_t  GetSigmaZ2(Int_t iz);
  Float_t  FitMax(Float_t vmatrix[5][5], Float_t y, Float_t z, Float_t sigmay, Float_t sigmaz);
  void AddCluster(AliTPCclusterMI &c);  // add the cluster to the array
  void UnfoldCluster(Float_t * matrix[7], Float_t recmatrix[5][5], 
		     Float_t & meani, Float_t & meanj, Float_t & sum, Float_t &overlap );
  void FindClusters(AliTPCCalROC * noiseROC);
  Double_t  ProcesSignal(Float_t * signal, Int_t nchannels, Int_t id[3], Double_t &rms, Double_t &pedestalCalib);
  void DumpHistos();


  Float_t * fBins;       //!digits array
  Float_t * fResBins;    //!digits array with res. after 1 finder
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
  Bool_t  fIsOldRCUFormat; // assume old RCU raw data format
  AliRawEventHeaderBase *fEventHeader; //! event header information
  UInt_t  fTimeStamp;   // Time Stamp
  UInt_t  fEventType;   // Event Type
  TTree * fInput;   //!input  tree with digits - object not owner
  TTree * fOutput;   //!output tree with digits - object not owner
  AliTPCClustersRow * fRowCl;  //! current cluster row
  AliSimDigits * fRowDig;      //! current digits row
  const AliTPCParam * fParam;        //! tpc parameters
  Int_t fNcluster;             // number of clusters - for given row
  TObjArray * fAmplitudeHisto;          //! array of histograms of amplitudes
  TTreeSRedirector *fDebugStreamer;     //!debug streamer
  const AliTPCRecoParam  * fRecoParam;        //! reconstruction parameters
  ClassDef(AliTPCclustererMI,1)  // Time Projection Chamber digits
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
  return kTRUE; 
}



//-----------------------------------------------------------------

#endif


