/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Container class for bad channels & bad runs identification
// Author: Olga Driga (SUBATECH)

#ifndef ALIANALYSISTASKCALOCELLSQA_H
#define ALIANALYSISTASKCALOCELLSQA_H

// --- ROOT system ---
#include <TString.h>

// --- AliRoot header files ---
#include <AliAnalysisTaskSE.h>
#include <AliCaloCellsQA.h>

class AliAnalysisTaskCaloCellsQA : public AliAnalysisTaskSE {

public:

  // detectors
  enum {
    kEMCAL = 0,
    kPHOS  = 1
// ,kDCAL  = 2      // not implemented
  };

  AliAnalysisTaskCaloCellsQA();
  AliAnalysisTaskCaloCellsQA(const char *name, Int_t nmods = 10, Int_t det = kEMCAL, const char *outfile = NULL);
  virtual ~AliAnalysisTaskCaloCellsQA();

  void   UserCreateOutputObjects();
  void   UserExec(Option_t *);
  void   Terminate(Option_t *);

  void   SetBadCells(Int_t badcells[], Int_t nbad);

  // getters and setters
  AliCaloCellsQA*  GetCaloCellsQA()                { return fCellsQA; }
  Bool_t           GetAvoidPileup()                { return fkAvoidPileup; }
  const char*      GetOutputFileName()             { return fOutfile.Data(); }
  void             SetAvoidPileup(Bool_t flag)     { fkAvoidPileup = flag; }
  void             SetOutputFileName(char* fname)  { fOutfile = fname; }

protected:
  Bool_t IsClusterBad(AliVCluster *clus);

private:
  AliAnalysisTaskCaloCellsQA(const AliAnalysisTaskCaloCellsQA &);
  AliAnalysisTaskCaloCellsQA & operator = (const AliAnalysisTaskCaloCellsQA &);

protected:
  Bool_t           fkAvoidPileup;     // flag not to process pileup events
  AliCaloCellsQA*  fCellsQA;          // analysis instance
  TString          fOutfile;          // output file name
  Int_t            fNBad;             // number of entries in fBadCells
  Int_t*           fBadCells;         //[fNBad] bad cells array

  ClassDef(AliAnalysisTaskCaloCellsQA, 2);
};

#endif
