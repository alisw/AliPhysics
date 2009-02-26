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
class TGraphErrors;
class AliTRDcheckESD : public AliAnalysisTask {
public:
  enum ETRDcheckESDstatus {
    kMC = BIT(0)
  };
  enum ETRDcheckESDhistos {
    kNCl  = 0 // number of clusters per track
   ,kTRDstat  // TRD tracks status
   ,kResults  // graphs as results
   ,kNhistos = 3 // number of histograms
  };
  enum ETRDcheckESDbits {
    kTPCout = 1 // track left TPC
   ,kTRDin      // track reach TRD fiducial volume
   ,kTRDout     // track reconstructed in TRD
   ,kTRDpid     // PID calculated in TRD
   ,kTRDref     // track refitted in TRD
   ,kNbits  = 5 // number of check bits
  };
  AliTRDcheckESD();
  virtual ~AliTRDcheckESD();
  
  void          ConnectInputData(Option_t *);
  void          CreateOutputObjects();
  TGraphErrors* GetGraph(Int_t id, Option_t *opt=0x0);
  void          Exec(Option_t *);

  Bool_t        HasMC() const { return TESTBIT(fStatus, kMC);}
  TObjArray*    Histos();
  Bool_t        Load(const Char_t *fn, const Char_t *name=0x0);
  void          SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  void          Terminate(Option_t *);

  static const Int_t   fgkNgraphs;
  static const Float_t fgkxTPC;
  static const Float_t fgkxTOF;

private:
  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  void          Process(TH1 **h, TGraphErrors *g);

  Int_t            fStatus;            // bit mask for controlling the task
  AliESDEvent      *fESD;              // ESD event
  AliMCEvent       *fMC;               // MC event
  TObjArray        *fHistos;           // QA histos
  ClassDef(AliTRDcheckESD, 1)          // user oriented TRD analysis based on ESD-MC data
};
#endif
