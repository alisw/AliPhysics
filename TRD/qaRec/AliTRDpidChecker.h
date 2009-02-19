#ifndef ALITRDPIDCHECKER_H
#define ALITRDPIDCHECKER_H

//////////////////////////////////////////////////////
//
// Task to check PID performance of the TRD
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//          Alex Bercuci <A.Bercuci@gsi.de>
//          Markus Fasel <M.Fasel@gsi.de>
//
///////////////////////////////////////////////////////

#ifndef ROOT_TAxis
#include "TAxis.h"
#endif

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDReconstructor;
class AliTRDpidUtil;
class AliTRDpidChecker : public AliTRDrecoTask 
{
  // Plots registered for this task
  enum{
     kEfficiency     =  0     // pi Efficiency plot
    ,kdEdx           =  1     // dE/dx spectra
    ,kdEdxSlice      =  2     // dE/dx spectra
    ,kPH             =  3     // pulse height spectra
    ,kNClus          =  4     //  number of clusters per track
    ,kMomentum       =  5     // momentum distribution
    ,kMomentumBin    =  6     // momentum distribution
    ,kThresh         =  7     // threshold in efficiency
    ,kNTracklets     =  8     // Number of tracklets per track
  };
public:
  AliTRDpidChecker();
  virtual ~AliTRDpidChecker();
  
  virtual void    CreateOutputObjects();
  virtual Bool_t  GetRefFigure(Int_t ifig);
  virtual Bool_t  PostProcess();
  virtual void    Terminate(Option_t *);

  TH1 *PlotLQ(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNN(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotESD(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdx(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdxSlice(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotPH(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNClus(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNTracklets(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMom(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMomBin(const AliTRDtrackV1 *track = 0x0);

  void SetRequireMinNTracklets(Int_t mintracklets) { fMinNTracklets = mintracklets; }
  void SetRequireMaxNTracklets(Int_t maxtracklets) { fMaxNTracklets = maxtracklets; }

  TObjArray *GetGraphs() { return fGraph; };
  //TObjArray *GetHistos() { return fContainer; };
  virtual TObjArray *Histos();
  void EvaluatePionEfficiency(TObjArray *histoContainer, TObjArray *results, Float_t electron_efficiency);
  inline void SetMomentumBinning(Int_t nBins, Double_t *bins);
  inline Int_t FindBin(Int_t species, Double_t momentum);
  inline Bool_t IsInRange(Double_t momentum);

private:
  AliTRDpidChecker(const AliTRDpidChecker&);               // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&);    // not implemented

  Int_t  CalcPDG(AliTRDtrackV1* track = 0x0);
  Bool_t CheckTrackQuality(const AliTRDtrackV1* track = 0x0);
  
  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID
  AliTRDpidUtil       *fUtil;              //! utility class for PID calculations
  TObjArray           *fGraph;             //! array of graphs filled in PostProcess
  TObjArray           *fEfficiency;        //! array of histograms with efficiency
  TAxis               *fMomentumAxis;      //! helper mementum binning
  Int_t                fMinNTracklets;     // minimum number of required Tracklets (for systematic studies)
  Int_t                fMaxNTracklets;     // maximum number of required Tracklets (for systematic studies) 
  ClassDef(AliTRDpidChecker, 1); // TRD PID checker
};

//________________________________________________________________________
inline void AliTRDpidChecker::SetMomentumBinning(Int_t nBins, Double_t *bins){
  //
  // Set the Momentum Bins
  //
  if(fMomentumAxis) delete fMomentumAxis;
  fMomentumAxis = new TAxis(nBins, bins);
}

//________________________________________________________________________
inline Int_t AliTRDpidChecker::FindBin(Int_t species, Double_t momentum){
  //
  // Find the Bin in the 2D Histogram
  //
  return species * fMomentumAxis->GetNbins() + (fMomentumAxis->FindBin(momentum) -1);
}

//________________________________________________________________________
inline Bool_t AliTRDpidChecker::IsInRange(Double_t momentum){
  //
  // Check Whether momentum is in the defined Range
  //
  return (momentum >= fMomentumAxis->GetXmin() && momentum <= fMomentumAxis->GetXmax());
}

#endif
