#ifndef ALIITSSIMULATIONSDD_H
#define ALIITSSIMULATIONSDD_H



#include <TNtuple.h>
#include <TArrayF.h>


#include "AliITSsimulation.h"

//___________________________________________________


class TH1F;
class TFile;
class TVector;
class TArrayI;
class TArrayF;
class AliITS;
class AliITSMap;
class AliITSMapA1;
class AliITSMapA2;
class AliITSetfSDD;
class AliITSInStream;
class AliITSresponse;

//___________________________________________________

class AliITSsimulationSDD : public AliITSsimulation {

public:

  AliITSsimulationSDD();
  AliITSsimulationSDD(AliITSsegmentation *seg, AliITSresponse *res);
  AliITSsimulationSDD(AliITSsimulationSDD &source);
  virtual ~AliITSsimulationSDD();
  AliITSsimulationSDD& operator=(AliITSsimulationSDD &source);

  // get the address of the array mapping the signal or pointers to arrays
  virtual AliITSMap*  HitMap(Int_t i);

  // set the scale size factor for the smples in FFT
  virtual void SetScaleFourier(Int_t scale=4) {fScaleSize=scale;}
  Int_t ScaleFourier() {return fScaleSize;}
  // set perpendicular tracks flag
  virtual void SetPerpendTracksFlag(Bool_t flag=1) {fFlag=1;}
  Bool_t PerpendTracksFlag() {return fFlag;} 
  // set compression parameters for 2D or 1D via response functions
  void SetCompressParam();
  // retrieve compression parameters for 2D or 1D
  void CompressionParam(Int_t i, Int_t &db, Int_t &tl, Int_t &th);
  void CompressionParam(Int_t i, Int_t &db, Int_t &tl);

  virtual Int_t Convert10to8(Int_t signal);
  virtual Int_t Convert8to10(Int_t signal);
  virtual void ZeroSuppression(const char *opt);
  virtual void Init2D();
  virtual void Compress2D();
  virtual void Init1D();
  virtual void Compress1D();
  virtual void StoreAllDigits();
  virtual void ReadBaseline();
  virtual void GetAnodeBaseline(Int_t i, Float_t &baseline, Float_t &noise);
  virtual void AddDigit(Int_t i, Int_t j, Int_t signal);
  virtual void  FindCluster
       (Int_t i, Int_t j,Int_t signal,Int_t minval,Bool_t &cond);


  // get parameters for 1D - this could be changed when we get more
  // input from Torino after they have a look at the code 
  virtual Int_t Tolerance(Int_t i) {return fTol[i];}
  virtual Int_t Disable(Int_t i) {return fT2[i];}
  // Set the output file name - for 1D encoding 
  virtual void SetFileName(const char *filnam) {fFileName=filnam;}

  void ChargeToSignal();
  void DigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev);
  void SortTracks(Int_t *tracks,Float_t *charges,Int_t *hits,Int_t ntracks);
  void ListOfFiredCells(Int_t *arg,Double_t timeAmplitude,TObjArray *list,
                        TClonesArray *padr);

  void CreateHistograms(Int_t scale);
  void FillHistograms();
  void ResetHistograms();
  // Get the pointer to the array of histograms
  TObjArray*  GetHistArray() {return fHis;}

  // create a separate tree for background monitoring (2D) 
  virtual  void  MakeTreeB(Option_t *option="B") 
    { fTreeB = new TNtuple("ntuple","2D backgr","nz:nl:nh:low:anode");}
  void           GetTreeB(Int_t) { }

  // Return pointer to TreeB
  TNtuple      *TreeB() {return fTreeB;} 

  void WriteToFile(TFile *fp);
  TH1F *GetAnode(Int_t wing, Int_t anode); 
  void SetCheckNoise(Bool_t check=kFALSE) {fCheckNoise=check;}
  Float_t GetNoise();
  void SetDoFFT(Int_t doFFT=1) {fDoFFT=doFFT;}


private:
  AliITS              *fITS;  //! local pointer to ITS

  AliITSMapA1         *fHitMap1; //! local pointer to map of digits
  AliITSMapA2         *fHitMap2; //! local pointer to map of signals
  AliITSInStream      *fStream;  //! input file stream
  AliITSetfSDD        *fElectronics; //! local pointer to electronics simulation
  
  TArrayI             fD;            // decrease values for baseline eq.
  TArrayI             fT1;           // low thresholds
  TArrayI             fT2;           // high thresholds(2D) or disable (1D) 
  TArrayI             fTol;          // tolerance
  TArrayF             fBaseline;     // Baseline
  TArrayF             fNoise;        // Noise value
  TNtuple            *fTreeB;        // Background info tree for 2D
  TString             fParam;        // Compresion algorithm options
  TString             fFileName;     // File name for possible options above

  Bool_t fFlag;         // Flag used to simulate perpendicular tracks
  Bool_t fCheckNoise;   // Flag used to check the simulated noise
  Int_t  fDoFFT;        // Flag used to switch off electronics when 0
  Int_t fNofMaps;       // Number of anodes used ( 1 - 2*nanodes per wing )
  Int_t fMaxNofSamples; // Number of time samples
  Int_t fScaleSize;     // scale size factor for the samples in FFT
  Int_t fModule;  // in case bgr, noise, param change module-by-module
  Int_t fEvent;   // solely for output from bgr monitoring of 2D
  
  TObjArray *fHis;             // just in case for histogramming

  Double_t            *fInZR;  // ! [fScaleSize*fMaxNofSamples]  
                               // input of the real part of FFT
  Double_t            *fInZI;  // ! [fScaleSize*fMaxNofSamples] 
                               // input of the imaginary part of FFT
  Double_t            *fOutZR; // ! [fScaleSize*fMaxNofSamples] 
                               // output of the real part of FFT
  Double_t            *fOutZI; // ! [fScaleSize*fMaxNofSamples] 
                               // output of the imaginary part of FFT

  ClassDef(AliITSsimulationSDD,1)  // Simulation of SDD clusters
    
};
#endif
