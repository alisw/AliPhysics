#ifndef ALIITSSIMULATIONSDD_H
#define ALIITSSIMULATIONSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */


/* $Id$ */

////////////////////////////////////////////////////////////
// Simulation class for SDD                               //
////////////////////////////////////////////////////////////

#include <TNtuple.h>
#include "AliITSsimulation.h"
#include "AliITSsegmentationSDD.h"
class TH1F;
class TFile;
class TArrayI;
class TArrayF;
class AliITS;
class AliITSpList;
class AliITSMap;
class AliITSMapA1;
class AliITSMapA2;
class AliITSetfSDD;
class AliITSCalibration;
class AliITSCalibrationSDD;

class AliITSsimulationSDD : public AliITSsimulation {
  public:
    AliITSsimulationSDD(); // default constructor
    //Standard Constructor
    AliITSsimulationSDD(AliITSDetTypeSim* dettyp);

    virtual ~AliITSsimulationSDD(); // Destructor

    //    virtual AliITSsimulation& operator=(const AliITSsimulation &source);
    // Initilize variables for this simulation
    void Init();

    // Get a pointer to the segmentation object
    virtual AliITSsegmentation* GetSegmentationModel(Int_t /*dt*/){return fDetType->GetSegmentationModel(1);}
    // set pointer to segmentation object
    virtual void SetSegmentationModel(Int_t /*dt*/, AliITSsegmentation *seg){fDetType->SetSegmentationModel(1,seg);}

    static Int_t ScaleFourier(const AliITSsegmentationSDD* seg) 
    {if(seg->Npx()==128) {return 8;} else {return 4;}} // returns the scale factor
    // set perpendicular tracks flag
    virtual void SetPerpendTracksFlag(Bool_t flag=kFALSE) {fFlag=flag;}
    // returns perpendicular track flag.
    Bool_t PerpendTracksFlag() const {return fFlag;} 
    // set crosstalk flag
    virtual void SetCrosstalkFlag(Bool_t flag=kFALSE) {fCrosstalkFlag=flag;}
    // return crosstalk flag
    Bool_t CrosstalkFlag() const {return fCrosstalkFlag;}
    void FastFourierTransform(Double_t *real, Double_t *imag, Int_t direction);
    virtual Int_t Convert10to8(Int_t signal) const;//10 to 8 bit SDD compresion
    virtual Int_t Convert8to10(Int_t signal) const;//8 to 10 bit decompresion
    virtual void Compress2D(); // Applies 2D compresion algorithm
    virtual void StoreAllDigits(); // if No compresion run this.
    // returns baseline and noise for a given anode i.
    //virtual void GetAnodeBaseline(Int_t i,Double_t &baseline,Double_t &noise) const;
    // local implementation of ITS->AddDigit. Specific for SDD
    virtual void AddDigit(Int_t i, Int_t j, Int_t signalc, Int_t signale);

    // add baseline, noise, gain, electronics and ADC saturation effects
    void ChargeToSignal(Int_t mod,Bool_t bAddNoise=kFALSE, Bool_t bAddGain=kTRUE);
    // add crosstalk effect
    void ApplyCrosstalk(Int_t mod);
    
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
    void WriteToFile(TFile *fp);// Writes the histograms to a file
    // Get's histogram of a particular anode.
    TH1F *GetAnode(Int_t wing, Int_t anode);

    // sets DoFFT value.
    void SetDoFFT(Int_t doFFT=1) {fDoFFT=doFFT;}

    // Print SSD simulation Parameters
    virtual void PrintStatus() const;

  private:
    AliITSsimulationSDD(const AliITSsimulationSDD &source);
    AliITSsimulationSDD& operator=(const AliITSsimulationSDD &source);

    // virtual void GetBaseline(Int_t mod);  // read baseline values from a file
    // Variables and pointers for local use only. Not Streamed out.
    AliITS         *fITS;          //! local pointer to ITS
    AliITSMapA2    *fHitMap2;      //! local pointer to map of signals
    AliITSMapA2    *fHitSigMap2;   //! local pointer to map of signals
    AliITSMapA2    *fHitNoiMap2;   //! local pointer to map of signals
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
    Bool_t     fFlag;         // Flag used to simulate perpendicular tracks
    Bool_t     fCrosstalkFlag; // Flag used to apply the crosstalk effect
    Int_t      fDoFFT;        // Flag used to switch off electronics when 0
    Int_t      fNofMaps;      // Number of anodes used ( 1-2*nanodes per wing )
    Int_t      fMaxNofSamples;// Number of time samples
    Int_t      fScaleSize;    // scale size factor for the samples in FFT

    ClassDef(AliITSsimulationSDD,3)  // Simulation of SDD clusters

};
#endif
