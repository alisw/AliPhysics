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
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliMCEvent;
class TH1;
class TH2;
class TObjArray;
class TGraph;
class TGraphErrors;
class AliTRDcheckESD : public AliAnalysisTaskSE {
public:
  enum ETRDcheckESDstatus {
     kMC        = BIT(0)  // use MC info
    ,kCollision = BIT(1)  // 
  };
  enum ETRDcheckESDhistos {
    kNCl  = 0    // number of clusters per track
   ,kTRDstat     // TRD tracks status
   ,kTRDmom      // TRD track momentum
   ,kPtRes       // Pt resolution @ vertex for TRD
   ,kNhistos = 4 // number of histograms
   ,kNrefs   = 4 // number of reference plots
  };
  enum ETRDcheckESDbits {
    kTPCout = 1 // track left TPC
   ,kTRDin      // track reach TRD fiducial volume
   ,kTRDout     // track reconstructed in TRD
   ,kTRDpid     // PID calculated in TRD
   ,kTRDref     // track refitted in TRD
  };
  AliTRDcheckESD();
  AliTRDcheckESD(char* name);
  virtual ~AliTRDcheckESD();
  
  void          UserCreateOutputObjects();
  Bool_t        GetRefFigure(Int_t ifig);
  Int_t         GetNRefFigures() const  { return fNRefFigures; } 
  void          UserExec(Option_t *);

  Bool_t        HasMC() const { return TESTBIT(fStatus, kMC);}
  Bool_t        IsCollision() const {return TESTBIT(fStatus, kCollision);}
  void          SetCollision(Bool_t set=kTRUE) {set ? SETBIT(fStatus, kCollision) : CLRBIT(fStatus, kCollision);}
  TObjArray*    Histos();
  Bool_t        Load(const Char_t *fn="AnalysisResults.root", const Char_t *dir="TRD_Performance", const Char_t *name=NULL);
  void          SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  Bool_t        PutTrendValue(const Char_t *name, Double_t val);
  void          Terminate(Option_t *);

private:
  static const Float_t fgkxTPC; // end radial position of TPC
  static const Float_t fgkxTOF; // start radial position of TOF
  static const UChar_t fgkNgraph[kNrefs]; // number of graphs/ref plot

  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  Int_t         Pdg2Idx(Int_t pdg);
  void          Process(TH1 **h, TGraphErrors *g);
  void          Process2D(TH2 * const h, TGraphErrors **g);
  void          PrintStatus(ULong_t s);

  Int_t            fStatus;            // bit mask for controlling the task
  Int_t            fNRefFigures;       // number of current ref plots
  AliESDEvent      *fESD;              //! ESD event
  AliMCEvent       *fMC;               //! MC event
  TObjArray        *fHistos;           //! QA histos
  TObjArray        *fResults;          // QA graphs
  static FILE      *fgFile;            //! trend file streamer
  // Vertex selection
  static const Float_t fgkEvVertexZ;// cm
  static const Int_t   fgkEvVertexN;// cm
  // Track selection
  static const Float_t fgkTrkDCAxy; // cm
  static const Float_t fgkTrkDCAz;  // cm
  static const Int_t   fgkNclTPC;   // N clusters TPC
  static const Float_t fgkPt;       // min. pt
  static const Float_t fgkEta;      // eta range

  ClassDef(AliTRDcheckESD, 5)          // user oriented TRD analysis based on ESD-MC data
};
#endif
