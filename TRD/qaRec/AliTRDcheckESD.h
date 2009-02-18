#ifndef AliTRDcheckESD_H
#define AliTRDcheckESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDcheckESD.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

class AliESDEvent;
class AliMCEvent;
class TObjArray;
class AliTRDcheckESD : public AliAnalysisTask {
public:
  enum ETRDcheckESDstatus {
    kMC = BIT(0)
  };
  enum ETRDcheckESDhistos {
    kNCl     = 0   // number of clusters per track
   ,kTRDstat = 2   // TRD tracks status
  };
  AliTRDcheckESD();
  virtual ~AliTRDcheckESD();
  
  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();

  Bool_t  HasMC() const { return TESTBIT(fStatus, kMC);}

  void    Exec(Option_t *);
  void    SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  void    Terminate(Option_t *);

  static const Float_t xTPC;
  static const Float_t xTOF;

private:
  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  Int_t            fStatus;            // bit mask for controlling the task
  AliESDEvent      *fESD;              // ESD event
  AliMCEvent       *fMC;               // MC event
  TObjArray        *fHistos;           // QA histos
  ClassDef(AliTRDcheckESD, 1)          // user oriented TRD analysis based on ESD-MC data
};
#endif
