#ifndef ALITRDCHECKESD_H
#define ALITRDCHECKESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDcheckESD.h 27496 2008-07-22 08:35:45Z cblume $ */

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

class AliESDEvent;
class AliMCEvent;
class TH1;
class TObjArray;
class TGraph;
class TGraphErrors;
class AliTRDcheckESD : public AliAnalysisTask {
public:
  enum ETRDcheckESDstatus {
    kMC = BIT(0)
  };
  enum ETRDcheckESDhistos {
    kNCl  = 0    // number of clusters per track
   ,kTRDstat     // TRD tracks status
   ,kTRDmom      // TRD track momentum
   ,kNhistos = 3 // number of histograms
   ,kNgraphs = 6 // number of graphs
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
  TGraph*       GetGraph(Int_t id, Option_t *opt="bc");
  void          Exec(Option_t *);

  Bool_t        HasMC() const { return TESTBIT(fStatus, kMC);}
  TObjArray*    Histos();
  Bool_t        Load(const Char_t *fn, const Char_t *name=0x0);
  void          SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  Bool_t        PutTrendValue(const Char_t *name, Double_t val);
  void          Terminate(Option_t *);

private:
  static const Float_t fgkxTPC; // end radial position of TPC
  static const Float_t fgkxTOF; // start radial position of TOF

  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  void          Process(TH1 **h, TGraphErrors *g);
  void          PrintStatus(ULong_t s);

  Int_t            fStatus;            // bit mask for controlling the task
  AliESDEvent      *fESD;              // ESD event
  AliMCEvent       *fMC;               // MC event
  TObjArray        *fHistos;           // QA histos
  TObjArray        *fResults;          // QA graphs
  static FILE *fgFile;                 //! trend file streamer
  ClassDef(AliTRDcheckESD, 3)          // user oriented TRD analysis based on ESD-MC data
};
#endif
