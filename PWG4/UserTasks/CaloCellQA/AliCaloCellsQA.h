/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Class for bad channels & bad runs identification
// Author: Olga Driga (SUBATECH)

#ifndef ALICALOCELLSQA_H
#define ALICALOCELLSQA_H

// --- ROOT system ---
#include <TObjArray.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>

class AliCaloCellsQA : public TObject {

public:

  // detectors
  enum {
    kEMCAL = 0,
    kPHOS  = 1
// ,kDCAL  = 2      // not yet implemented
  };

  AliCaloCellsQA();
  AliCaloCellsQA(Int_t nmods, Int_t det = kEMCAL, Int_t startRunNumber = 100000, Int_t endRunNumber = 200000);

  virtual ~AliCaloCellsQA();

  virtual void ActivateFullAnalysis();
  virtual void InitSummaryHistograms(Int_t nbins = 400, Double_t emax = 4.,
                                     Int_t nbinsh = 100, Double_t emaxh = 300.,
                                     Int_t nbinst = 250, Double_t tmin = 0.4e-6, Double_t tmax = 0.85e-6);
  virtual void InitTransientFindCurrentRun(Int_t runNumber);
  virtual void Fill(Int_t runNumber, TObjArray *clusArray, AliVCaloCells *cells, Double_t vertexXYZ[3]);  // main method

  // getters
  virtual Int_t      GetDetector()     { return fDetector; }
  virtual Int_t      GetNMods()        { return fNMods; }
  virtual Double_t   GetClusElowMin()  { return fClusElowMin; }
  virtual Double_t   GetClusEhighMin() { return fClusEhighMin; }
  virtual Double_t   GetPi0EClusMin()  { return fPi0EClusMin; }
  virtual Bool_t     GetFullAnalysis() { return fkFullAnalysis; }

  virtual TObjArray* GetListOfHistos() { return fListOfHistos; }

  // setters
  virtual void SetClusterEnergyCuts(Double_t pi0EClusMin = 0.5, Double_t ElowMin = 0.3, Double_t EhighMin = 1.0);
  virtual void SetBinningParameters(Int_t nbins1 = 100, Int_t nbins2 = 250, Int_t nbins3x = 200, Int_t nbins3y = 30,
                                    Double_t xmax1 = 5., Double_t xmax2 = 0.5, Double_t xmax3 = 20.);

protected:

  virtual void    Init(Int_t nmods, Int_t det, Int_t startRunNumber, Int_t endRunNumber);
  virtual void    InitTransientMembers(Int_t run);
  virtual void    InitHistosForRun(Int_t run);

  virtual void    FillCellsInCluster(TObjArray *clusArray, AliVCaloCells *cells);
  virtual void    FillJustCells(AliVCaloCells *cells);
  virtual void    FillPi0Mass(TObjArray *clusArray, Double_t vertexXYZ[3]);

  virtual Int_t   CheckClusterGetSM(AliVCluster* clus);
  virtual Int_t   GetSM(Int_t absId);
  virtual Bool_t  IsCellLocalMaximum(Int_t c, AliVCluster* clus, AliVCaloCells* cells);
  virtual Bool_t  IsCellLocalMaximum(Int_t absId, AliVCaloCells* cells);
  virtual void    AbsIdToSMEtaPhi(Int_t absId, Int_t &sm, Int_t &eta, Int_t &phi);

private:

  AliCaloCellsQA(const AliCaloCellsQA &);
  AliCaloCellsQA & operator = (const AliCaloCellsQA &);

private:

  // changeable parameters
  Int_t     fDetector;                  // kEMCAL or kPHOS
  Int_t     fNMods;                     // maximum supermodule number + 1 (4 or 10)
  Double_t  fClusElowMin;               // minimum cluster energy cut for low energy per run histograms
  Double_t  fClusEhighMin;              // minimum cluster energy cut for high energy per run histograms
  Double_t  fPi0EClusMin;               // minimum cluster energy cut for pi0 mass histograms
  Bool_t    fkFullAnalysis;             // flag to activate all available histograms

  // changeable binning parameters
  Int_t     fNBinsECells;               // number of bins in hECells
  Int_t     fNBinsPi0Mass;              // number of bins in hPi0Mass
  Int_t     fNBinsXNCellsInCluster;     // number of bins in hNCellsInCluster, X axis
  Int_t     fNBinsYNCellsInCluster;     // number of bins in hNCellsInCluster, Y axis
  Double_t  fXMaxECells;                // X axis maximum in hECells
  Double_t  fXMaxPi0Mass;               // X axis maximum in hPi0Mass
  Double_t  fXMaxNCellsInCluster;       // X axis maximum in hNCellsInCluster

  // internal parameters
  Int_t     fRunNumbers[1000];          // already encountered runs
  Int_t     fNRuns;                     // number of encountered runs
  Int_t     fRI;                        //! current (cached) run index

  // internal parameters, used for coding convenience
  Int_t     fAbsIdMin;                  // minimum absId number (0/EMCAL, 1/PHOS)
  Int_t     fAbsIdMax;                  // maximum absId number + 1

  TObjArray  *fListOfHistos;            // array with all the histograms


  /* All the histograms below are present in fListOfHistos.
   *
   * NOTE: strictly speaking, the transient members below are not necessary, but the
   * caching mechanism used in this class boosts the analysis speed considerably
   * (no searching for a histogram by its name is performed in each event).
   */

  TH1D *fhNEventsProcessedPerRun;             //! number of processed events per run

  // current run histograms; X axis -- cell absId number;
  TH1F *fhCellLocMaxNTimesInClusterElow;      //! number of times cell was local maximum in a low energy cluster
  TH1F *fhCellLocMaxNTimesInClusterEhigh;     //! number of times cell was local maximum in a high energy cluster
  TH1F *fhCellLocMaxETotalClusterElow;        //! total cluster energy for local maximum cell, low energy
  TH1F *fhCellLocMaxETotalClusterEhigh;       //! total cluster energy for local maximum cell, high energy
  TH1F *fhCellNonLocMaxNTimesInClusterElow;   //! number of times cell wasn't local maximum in a low energy cluster
  TH1F *fhCellNonLocMaxNTimesInClusterEhigh;  //! number of times cell wasn't local maximum in a high energy cluster
  TH1F *fhCellNonLocMaxETotalClusterElow;     //! total cluster energy for not local maximum cell, low energy
  TH1F *fhCellNonLocMaxETotalClusterEhigh;    //! total cluster energy for not local maximum cell, high energy

  // current run, per supermodule histograms; the maximum number of supermodules is 10
  TH1F *fhECells[10];                         //! cell amplitude distribution
  TH1F *fhPi0Mass[10][10];                    //! pi0 mass spectrum
  TH2F *fhNCellsInCluster[10];                //! distribution of number of cells in cluster vs cluster energy

  // summary histograms: cells spectra at different conditions; X axis -- cell absId number
  TH2F *fhCellAmplitude;                      //! amplitude distribution per cell
  TH2F *fhCellAmplitudeEhigh;                 //! amplitude distribution per cell, high energies
  TH2F *fhCellAmplitudeNonLocMax;             //! amplitude distribution per not a local maximum cell
  TH2F *fhCellAmplitudeEhighNonLocMax;        //! amplitude distribution per not a local maximum cell, high energies
  TH2F *fhCellTime;                           //! time distribution per cell

  ClassDef(AliCaloCellsQA,2)
};

#endif
