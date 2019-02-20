#ifndef ALITRDV0MONITOR_H
#define ALITRDV0MONITOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
//
//  Monitor V0 for TRD
//
//  Authors:                                          
//  Markus Heide <mheide@uni-muenster.de> 
//////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif
#ifndef ALIPID_H
#include "AliPID.h"
#endif
#ifndef ALITRDV0INFO_H
#include "AliTRDv0Info.h"
#endif

class TTree;
class TH2F;
class TH1I;
class TObjArray;
class AliTRDv0Monitor : public AliTRDrecoTask
{
public:
  enum ETRDv0Monitor {
    kNSamples = 3
    ,kNPlots           =  13    // Number of plots for this task 
    ,kNCutSteps = 3
    ,kNDets = 3
  };

  AliTRDv0Monitor();
  AliTRDv0Monitor(const char *name);
  virtual ~AliTRDv0Monitor(){};

  Bool_t      GetRefFigure(Int_t ifig); 
  TObjArray*  Histos();
  void        UserExec(Option_t *option);
  //void        MakeSummary();

private:
  AliTRDv0Monitor(const AliTRDv0Monitor&);              // not implemented
  AliTRDv0Monitor& operator=(const AliTRDv0Monitor&);   // not implemented

  TH1I *fhCutReductions[AliPID::kSPECIES];//!histo for sample reductions by each ID cut
  TH1I *fhQualityReductions;//!histo for sample reductions by each quality cut
  TH2F *fhV0Chi2ndf[AliTRDv0Info::kNDecays][kNCutSteps];//!Chi2/ndf distributions before cuts, after inv. mass cut, after all cuts (same for all arrays below!!!)
  TH2F *fhInvMass[AliTRDv0Info::kNDecays];//!invariant mass distributions
  TH2F *fhPsiPair[AliTRDv0Info::kNDecays][kNCutSteps];//!Psi_pair angle distributions
  TH2F *fhPointAngle[AliTRDv0Info::kNDecays][kNCutSteps];//!pointing angle
  TH2F *fhDCA[AliTRDv0Info::kNDecays][kNCutSteps];//!Distance of closest approach between daughters
  TH2F *fhOpenAngle[AliTRDv0Info::kNDecays][kNCutSteps];//!opening angle between daughters
  TH2F *fhDetPID[kNDets][AliPID::kSPECIES];//!likelihood outputs from different detectors
  TH2F *fhComPID[AliPID::kSPECIES];//!combined PID from TPC and TOF
  TH2F *fhRadius[AliTRDv0Info::kNDecays][kNCutSteps];//!radial distance of secondary vertex to primary vertex
  TH2F *fhTPCdEdx[AliPID::kSPECIES][kNCutSteps];//!energy deposition in TPC
 

 
  TObjArray     *fV0s;                  //! v0 array
  TTree         *fData;                 //! dEdx-P data
  TObjArray     *fInfo;                 //! list of PID info
  Float_t       fP;                     //! momentum
  Float_t       fPID[AliPID::kSPECIES]; //! pid from v0s

  ClassDef(AliTRDv0Monitor, 3); // V0 Monitor
};

#endif
