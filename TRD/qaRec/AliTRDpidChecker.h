#ifndef ALITRDPIDCHECKER_H
#define ALITRDPIDCHECKER_H

//////////////////////////////////////////////////////
//
// Task to check PID performance of the TRD
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//
///////////////////////////////////////////////////////

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
  };
  // PID methods
  enum {
     kLQ   = 0 // 2D likelihood method
    ,kNN   = 1 // Neural network method
    ,kESD  = 2 // ESD results - check offline
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
  TH1 *PlotMom(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMomBin(const AliTRDtrackV1 *track = 0x0);

  TObjArray *GetGraphs() { return fGraph; };
  virtual TObjArray *Histos();

private:
  AliTRDpidChecker(const AliTRDpidChecker&);               // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&);    // not implemented

  Int_t  CalcPDG(AliTRDtrackV1* track = 0x0);
  Bool_t CheckTrackQuality(const AliTRDtrackV1* track = 0x0);
  
  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID
  AliTRDpidUtil       *fUtil;              //! utility class for PID calculations
  TObjArray           *fGraph;             //! array of graphs filled in PostProcess
  TObjArray           *fEfficiency;        //! array of histograms with efficiency

  ClassDef(AliTRDpidChecker, 1); // TRD PID checker
};

#endif
