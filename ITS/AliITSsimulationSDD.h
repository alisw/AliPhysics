#ifndef ALIITSSIMULATIONSDD_H
#define ALIITSSIMULATIONSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

#include <TNtuple.h>
#include <TArrayF.h>

#include "AliITSsimulation.h"

class TH1F;
class TFile;
class TVector;
class TArrayI;
class TArrayF;
class AliITS;
class AliITSpList;
class AliITSMap;
class AliITSMapA1;
class AliITSMapA2;
class AliITSetfSDD;
class AliITSsegmentationSDD;
class AliITSInStream;
class AliITSresponse;
class AliITSresponseSDD;

class AliITSsimulationSDD : public AliITSsimulation {
  public:
    AliITSsimulationSDD(); // default constructor
    //Standard Constructor
    AliITSsimulationSDD(AliITSsegmentation *seg, AliITSresponse *res);
    // Copy opporator
    AliITSsimulationSDD(AliITSsimulationSDD &source);
    virtual ~AliITSsimulationSDD(); // Destructor
    // = opporator
    AliITSsimulationSDD& operator=(const AliITSsimulationSDD &source);
    virtual AliITSsimulationSDD& operator=(const AliITSsimulation &source);
    // Initilize variables for this simulation
    void Init();

    // get the address of the array mapping the signal or pointers to arrays
    //virtual AliITSMap*  HitMap(Int_t i);

    // Return The Responce class
    AliITSresponseSDD* GetResp(){return (AliITSresponseSDD*)fResponse;}
    // Return The Segmentation class
    AliITSsegmentationSDD* GetSeg(){
        return (AliITSsegmentationSDD*)fSegmentation;}
    // set the scale size factor for the smples in FFT
    virtual void SetScaleFourier(Int_t scale=4) {fScaleSize=scale;}
    Int_t ScaleFourier() const {return fScaleSize;} // returns the scale factor
    // set perpendicular tracks flag
    virtual void SetPerpendTracksFlag(Bool_t flag=kFALSE) {fFlag=flag;}
    // returns perpendicular track flag.
    Bool_t PerpendTracksFlag() const {return fFlag;} 
    // set crosstalk flag
    virtual void SetCrosstalkFlag(Bool_t flag=kFALSE) {fCrosstalkFlag=flag;}
    // return crosstalk flag
    Bool_t CrosstalkFlag() const {return fCrosstalkFlag;}
    // set compression parameters for 2D or 1D via response functions
    void SetCompressParam();
    // retrieve compression parameters for 2D or 1D
    void CompressionParam(Int_t i, Int_t &db, Int_t &tl, Int_t &th);
    // retrieve compression parameters for 2D or 1D
    void CompressionParam(Int_t i, Int_t &db, Int_t &tl);

    virtual Int_t Convert10to8(Int_t signal) const;//10 to 8 bit SDD compresion
    virtual void ZeroSuppression(const char *opt); // Apply zero suppresion
    virtual void Init2D();   // initiilzes 2D compresion algorithm
    virtual void Compress2D(); // Applies 2D compresion algorithm
    virtual void Init1D();   // initilizes 1D compresion algorithm
    virtual void Compress1D(); // Applies 1D compresion algorithm
    virtual void StoreAllDigits(); // if No compresion run this.
    virtual void ReadBaseline();  // read baseline values from a file
    // returns baseline and noise for a given anode i.
    virtual void GetAnodeBaseline(Int_t i,Double_t &baseline,Double_t &noise);
    // local implementation of ITS->AddDigit. Specific for SDD
    virtual void AddDigit(Int_t i, Int_t j, Int_t signal);
    // Finds clulsters of signals. Use with regards to Compresion algorithms
    virtual void  FindCluster(Int_t i, Int_t j,Int_t signal,
                              Int_t minval,Bool_t &cond);

    // get parameters for 1D - this could be changed when we get more
    // input from Torino after they have a look at the code 
    virtual Int_t Tolerance(Int_t i) {return fTol[i];}//returns tolerance
    virtual Int_t Disable(Int_t i)  {return fT2[i];}//high threshold  2D
    // Set the output file name - for 1D encoding 
    virtual void SetFileName(const char *filnam) {fFileName=filnam;}

    // add baseline, noise, electronics and ADC saturation effects
    void ChargeToSignal(Bool_t bAddNoise=kFALSE);
    // add dead channels
    void ApplyDeadChannels();
    // add crosstalk effect
    void ApplyCrosstalk();
    
    // create maps to build the lists of tracks for each summable digit
    void InitSimulationModule( Int_t module, Int_t event );
    // clear maps
    void ClearMaps();
    // Summable Digitses a SDD module
    void SDigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev);
    // Add Summable digits to module maps.
    Bool_t AddSDigitsToModule( TClonesArray *pItemArray, Int_t mask );
    // digitize module from the sum of summable digits.
    void FinishSDigitiseModule();
    // Writes summable digits
    void WriteSDigits();
    // Introduces electronics effects and does zero-suppresion if required
    void FinishDigits();
    // Digitses a SDD module
    void DigitiseModule(AliITSmodule *mod,Int_t md,Int_t ev);
    // Spread charge in a SDD module
    void HitsToAnalogDigits(AliITSmodule *mod);
    // Sorts tracks for the 3 most highly contributed one to be added to digit.
    //void SortTracks(Int_t *tracks,Float_t *charges,Int_t *hits
    //                Int_t ntracks);
    // collects and returns the fired SDD cells (uses AliITSMapA2...).
    //void ListOfFiredCells(Int_t *arg,Double_t timeAmplitude,TObjArray *list,
    //		  TClonesArray *padr);

    // Creates histograms of maps for debugging
    void CreateHistograms(Int_t scale);
    // Fills histograms of maps for debugging
    void FillHistograms();
    // Resets histograms of maps for debugging
    void ResetHistograms();
    // Get the pointer to the array of histograms
    TObjArray*  GetHistArray() {return fHis;}
    // create a separate tree for background monitoring (2D) 
    virtual  void  MakeTreeB(Option_t *option="B") 
        { if(strstr(option,"B"))
            fTreeB = new TNtuple("ntuple", "2D backgr","nz:nl:nh:low:anode");}
    // presently a dummy routine use TreeB() instead
    void GetTreeB(Int_t) { }
    // Return pointer to TreeB
    TNtuple *TreeB() {return fTreeB;}
    void WriteToFile(TFile *fp);// Writes the histograms to a file
    // Get's histogram of a particular anode.
    TH1F *GetAnode(Int_t wing, Int_t anode);

    // Sets the check noise flag.
    void SetCheckNoise(Bool_t check=kFALSE) {fCheckNoise=check;}
    // Returns the noise value
    Float_t GetNoise();
    // sets DoFFT value.
    void SetDoFFT(Int_t doFFT=1) {fDoFFT=doFFT;}

    // Print SSD simulation Parameters
    virtual void PrintStatus() const;

  private:
    // Variables and pointers for local use only. Not Streamed out.
    AliITS         *fITS;          //! local pointer to ITS
    AliITSMapA2    *fHitMap2;      //! local pointer to map of signals
    AliITSMapA2    *fHitSigMap2;   //! local pointer to map of signals
    AliITSMapA2    *fHitNoiMap2;   //! local pointer to map of signals
    AliITSInStream *fStream;       //! input file stream
    AliITSetfSDD   *fElectronics;  //! local pointer to electronics simulation
    Double_t       *fInZR;         //! [fScaleSize*fMaxNofSamples] input of the
                                   // real part of FFT
    Double_t       *fInZI;         //! [fScaleSize*fMaxNofSamples] 
                                   // input of the imaginary part of FFT
    Double_t       *fOutZR;        //! [fScaleSize*fMaxNofSamples] 
                                   // output of the real part of FFT
    Double_t       *fOutZI;        //! [fScaleSize*fMaxNofSamples] 
                                   // output of the imaginary part of FFT
    Bool_t         *fAnodeFire;     //! [#of anodes] Flag if there is a signal

    TObjArray *fHis;          // just in case for histogramming
    TArrayI    fD;            // decrease values for baseline eq.
    TArrayI    fT1;           // low thresholds
    TArrayI    fT2;           // high thresholds(2D) or disable (1D) 
    TArrayI    fTol;          // tolerance
    TArrayF    fBaseline;     // Baseline
    TArrayF    fNoise;        // Noise value
    TNtuple   *fTreeB;        // Background info tree for 2D
    TString    fParam;        // Compresion algorithm options
    TString    fFileName;     // File name for possible options above
    Bool_t     fFlag;         // Flag used to simulate perpendicular tracks
    Bool_t     fCheckNoise;   // Flag used to check the simulated noise
    Bool_t     fCrosstalkFlag; // Flag used to apply the crosstalk effect
    Int_t      fDoFFT;        // Flag used to switch off electronics when 0
    Int_t      fNofMaps;      // Number of anodes used ( 1-2*nanodes per wing )
    Int_t      fMaxNofSamples;// Number of time samples
    Int_t      fScaleSize;    // scale size factor for the samples in FFT

    ClassDef(AliITSsimulationSDD,1)  // Simulation of SDD clusters

};
#endif
