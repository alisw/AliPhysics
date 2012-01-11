#ifndef ALITRDCHECKPID_H
#define ALITRDCHECKPID_H

//////////////////////////////////////////////////////
//
// PID performance checker of the TRD
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

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class AliTRDpidUtil;
class AliTRDcheckPID : public AliTRDrecoTask 
{
public:
  // Plots registered for this task
  enum{
    kEfficiency        =  0     // e Efficiency plot
    ,kdEdx             =  1     // dE/dx spectra
    ,kdEdxSlice        =  2     // dE/dx spectra
    ,kPH               =  3     // pulse height spectra
    ,kNClus            =  4     //  number of clusters per track
    ,kMomentum         =  5     // momentum distribution
    ,kMomentumBin      =  6     // momentum distribution
    ,kNTracklets       =  7     // Number of tracklets per track
    ,kEfficiencyMu     =  8     // mu Efficiency plot
    ,kEfficiencyPi     =  9     // pi Efficiency plot
    ,kEfficiencyKa     =  10    // K Efficiency plot
    ,kEfficiencyPr     =  11    // pr Efficiency plot
    ,kV0               =  12    // V0 performance
    ,kdQdl             =  13  // Plot total charge used for the 1D-Likelihood calculation
    ,kNPlots           =  14    // Number of plots for this tasks
 };
  AliTRDcheckPID();
  AliTRDcheckPID(char* name);
  virtual ~AliTRDcheckPID();
  
  virtual void    UserCreateOutputObjects();
  virtual Bool_t  GetRefFigure(Int_t ifig);
  virtual void    UserExec(Option_t *opt);
  virtual Bool_t  PostProcess();

  TH1 *PlotLQ(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNN(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotESD(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdQdl(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdx(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdxSlice(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotPH(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNClus(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNTracklets(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMom(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMomBin(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotV0(const AliTRDtrackV1 *track = 0x0);

  void SetRequireMinNTracklets(Int_t mintracklets) { fMinNTracklets = mintracklets; }
  void SetRequireMaxNTracklets(Int_t maxtracklets) { fMaxNTracklets = maxtracklets; }

  TObjArray *GetGraphs() const { return fGraph; };
  static Char_t const* MethodName(Int_t id) { return fgMethod[id]; };
  //TObjArray *GetHistos() { return fContainer; };
  virtual TObjArray *Histos();
  void EvaluateEfficiency(const TObjArray* const histoContainer, TObjArray *results, Int_t species, Float_t electronEfficiency);
  inline void SetMomentumBinning(Int_t nBins, Double_t *bins);
  inline Int_t FindBin(Int_t species, Double_t momentum);
  inline Bool_t IsInRange(Double_t momentum);

private:
  AliTRDcheckPID(const AliTRDcheckPID&);               // not implemented
  AliTRDcheckPID& operator=(const AliTRDcheckPID&);    // not implemented

  Int_t  CalcPDG(AliTRDtrackV1* track = 0x0);
  Bool_t CheckTrackQuality(const AliTRDtrackV1* track = 0x0) const;
  void   LocalInit();

  static Char_t const *fgMethod[3];        // PID method name
  AliTRDpidUtil       *fUtil;              // utility class for PID calculations
  TObjArray           *fGraph;             //! array of graphs filled in PostProcess
  TObjArray           *fPID;               //! array of PID info/track for calibration
  TObjArray           *fV0s;               //! array of V0 info
  TObjArray           *fEfficiency[AliPID::kSPECIES];      //! array of histograms with efficiency
  TAxis               *fMomentumAxis;      // helper mementum binning
  Int_t                fMinNTracklets;     // minimum number of required Tracklets (for systematic studies)
  Int_t                fMaxNTracklets;     // maximum number of required Tracklets (for systematic studies) 
  ClassDef(AliTRDcheckPID, 3); // TRD PID checker
};

//________________________________________________________________________
inline void AliTRDcheckPID::SetMomentumBinning(Int_t nBins, Double_t *bins){
  //
  // Set the Momentum Bins
  //
  if(fMomentumAxis) delete fMomentumAxis;
  fMomentumAxis = new TAxis(nBins, bins);
}

//________________________________________________________________________
inline Int_t AliTRDcheckPID::FindBin(Int_t species, Double_t momentum){
  //
  // Find the Bin in the 2D Histogram
  //
  return species * fMomentumAxis->GetNbins() + (fMomentumAxis->FindBin(momentum) -1);
}

//________________________________________________________________________
inline Bool_t AliTRDcheckPID::IsInRange(Double_t momentum){
  //
  // Check Whether momentum is in the defined Range
  //
  return (momentum >= fMomentumAxis->GetXmin() && momentum <= fMomentumAxis->GetXmax());
}

#endif
