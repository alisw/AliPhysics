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

  enum{
    kLQlikelihood    = 0     // place for 2-dim LQ electron likelihood distributions
    ,kNNlikelihood   = 1     // place for NN electron likelihood distributions
    ,kdEdx           = 2     // place for the dE/dx spectra
    ,kdEdxSlice      = 3     // place for the dE/dx spectra
    ,kPH             = 4     // place for pulse height spectra
    ,kNClus          = 5     // place for the number of clusters per track
    ,kMomentum       = 6     // place for the momentum distribution
    ,kMomentumBin    = 7     // place for the momentum distribution
    ,kGraphLQ        = 8     // place for the 2-dim LQ pion efficiencies
    ,kGraphNN        = 9     // place for the NN pion efficiencies
  };

  enum{
    kGraphStart = kGraphLQ
  };

public:
  AliTRDpidChecker();
  virtual ~AliTRDpidChecker();
  
  virtual void    CreateOutputObjects();
  virtual void    GetRefFigure(Int_t ifig);
  virtual Bool_t  PostProcess();
  virtual void    Terminate(Option_t *);

  TH1 *PlotLQ(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNN(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdx(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotdEdxSlice(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotPH(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNClus(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMom(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotMomBin(const AliTRDtrackV1 *track = 0x0);

  virtual TObjArray *Histos();

private:
  AliTRDpidChecker(const AliTRDpidChecker&);               // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&);    // not implemented

  Int_t  CalcPDG(AliTRDtrackV1* track = 0x0);
  Bool_t CheckTrackQuality(const AliTRDtrackV1* track = 0x0);
  
  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID
  AliTRDpidUtil       *fUtil;              //! utility class for PID calculations

  ClassDef(AliTRDpidChecker, 1); // TRD PID checker
};

#endif
