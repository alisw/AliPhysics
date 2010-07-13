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
#include "info/AliTRDv0Info.h"
#endif

class TTree;
class TList;
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
  AliTRDv0Monitor(const char *name, const char *title);

  virtual ~AliTRDv0Monitor();
  
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *option);
 
  TList *fOutput;         //! Container for output histos
  
  TH1I *hCutReductions[AliPID::kSPECIES];
  TH1I *hQualityReductions;
  TH2F *hV0Chi2ndf[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hInvMass[AliTRDv0Info::kNDecays];
  TH2F *hPsiPair[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hPointAngle[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hDCA[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hOpenAngle[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hDetPID[kNDets][AliPID::kSPECIES];
  TH2F *hComPID[AliPID::kSPECIES];
  TH2F *hRadius[AliTRDv0Info::kNDecays][kNCutSteps];
  TH2F *hTPCdEdx[AliPID::kSPECIES][kNCutSteps];
 


protected:
 
  TObjArray     *fV0s;                  //! v0 array
  TTree         *fData;                 //! dEdx-P data
  TObjArray     *fInfo;                 //! list of PID info
  
  Float_t       fP;                     // momentum
 
  Float_t       fPID[AliPID::kSPECIES]; // pid from v0s
  //Int_t GetPDG(Int_t index);            // PDG information from MC

private:
  AliTRDv0Monitor(const AliTRDv0Monitor&);              // not implemented
  AliTRDv0Monitor& operator=(const AliTRDv0Monitor&);   // not implemented

  ClassDef(AliTRDv0Monitor, 3); // V0 Monitor
};

#endif
