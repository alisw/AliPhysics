#ifndef ALICOMPARISONRESTASK_H
#define ALICOMPARISONRESTASK_H

//------------------------------------------------------------------------------
// Class to compare properties of reconstructed and MC particle tracks. 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class AliComparisonRes;
class AliComparisonEff;
class AliComparisonDEdx;
class AliComparisonDCA;
class AliMagFMaps;
class TList;

#include "AliAnalysisTask.h"

class AliComparisonTask : public AliAnalysisTask {
 public:
  AliComparisonTask(const char *name = "AliComparisonTask");
  virtual ~AliComparisonTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // Read TTree entry (event by event)
  Bool_t  ReadEntry(Int_t evt);

  // Set comparison objects
  void SetAliComparisonRes(AliComparisonRes* comp) {fCompRes = comp;}
  void SetAliComparisonEff(AliComparisonEff* comp) {fCompEff = comp;}
  void SetAliComparisonDEdx(AliComparisonDEdx* comp) {fCompDEdx = comp;}
  void SetAliComparisonDCA(AliComparisonDCA* comp) {fCompDCA = comp;}

  void SetMagField(Int_t mag = 2) {fMagField = mag;}
  void SetGeometry(char* geom = "/d/alice12/jacek/sim/v4-10-Release/pp/0/geometry.root")  {fGeom = geom;}

 private:
  TTree* fTree;                   //! input tree
  AliMCInfo *fInfoMC;             //! AliMCInfo object
  AliESDRecInfo *fInfoRC;         //! AliESDRecInfo object
  AliComparisonRes* fCompRes;     // TPC resolution comparison object
  AliComparisonEff* fCompEff;     // TPC efficiency comparison object
  AliComparisonDEdx* fCompDEdx;   // TPC DEdx comparison object
  AliComparisonDCA* fCompDCA;     // TPC DCA comparison object

  TList* fOutput;                 //! list send on output slot 0
  static Int_t evtNumber;         //! event number
  Int_t  fMagField;               //! mag. field (0 - 0.2 T, 1 - 0.4 T, 2 - 0.5 T) 
  AliMagFMaps *fMagFMap;          //! mag. field map 
  const char *fGeom;              //! ROOT file with detector geometry

  AliComparisonTask(const AliComparisonTask&); // not implemented
  AliComparisonTask& operator=(const AliComparisonTask&); // not implemented
  
  ClassDef(AliComparisonTask, 1); // example of analysis
};

#endif
