#ifndef ALIANALYSISGOODIES_H
#define ALIANALYSISGOODIES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Various utilities usefull for analysis
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h> 
#include <TStopwatch.h> 
#include <TString.h> 

#include "AliLog.h"

class AliAnalysisTask ; 
class TChain ; 
class TClass ; 
class AliEventTagCuts ;  
class AliRunTagCuts ;  

class AliAnalysisGoodies : public TObject {

public:
  AliAnalysisGoodies() ; 
  ~AliAnalysisGoodies() {;} 

  virtual void Help() const; 
  const Bool_t Alien2Local(const TString collectionNameIn, const TString localDir) ; 
  const Bool_t Make(AliRunTagCuts *runCuts, AliEventTagCuts *evtCuts, const char * in, const char * out) const  ; 
  const Bool_t Merge(const char * collection, const char * subFile = 0, const char * outFile = 0) ; 
  const Bool_t Register( const char * lfndir, const char * pfndir, const char * file)  ;   
  const Bool_t Process(TChain * chain) ;  
  const Bool_t Process(const char * esdFile)  ;  
  const Bool_t Process(const char * esdFile, AliRunTagCuts * runCuts, AliEventTagCuts * evtCuts)  ;  
  const Bool_t Process(const char * esdFile, const char * runCuts, const char * evtCuts)  ;  
  void         SetESDTreeName(const char * name) { fESDTreeName = name ; }
  void         SetTasks(Int_t nb, AliAnalysisTask ** taskList, TClass ** inputType, TClass ** outputType) ;
  const Bool_t MakeEsdCollectionFromTagFile(AliRunTagCuts * runCuts, AliEventTagCuts * evCuts, const char * in, const char * out) const  ; 

private:
  const Bool_t MakeEsdCollectionFromTagFile(const char * , const char * , const char * , const char *) const 
  { AliError("Not implemented") ; return 0 ;}
  const Bool_t MakeEsdCollectionFromTagCollection(AliRunTagCuts * runCuts, AliEventTagCuts * evtCuts, const char * in, const char * out) const ; 
  const Bool_t MakeEsdCollectionFromTagCollection(const char * , const char * , const char * , const char * ) const 
  { AliError("Not implemented") ; return 0 ;}
  const Bool_t ProcessChain(TChain * chain) const ; 
  const Bool_t ProcessEsdFile(const char * esdFile) const ;
  const Bool_t ProcessTagFile(const char * tagFile, AliRunTagCuts *runCuts, AliEventTagCuts *evtCuts) const ;
  const Bool_t ProcessTagFile(const char * tagFile, const char * runCuts, const char * evtCuts) const ;
  const Bool_t ProcessEsdXmlCollection(const char * esdFile) const ;
  const Bool_t ProcessTagXmlCollection(const char * esdFile, AliRunTagCuts * runCuts, AliEventTagCuts * evtCuts) const ;
  const Bool_t ProcessTagXmlCollection(const char * esdFile, const char * runCuts, const char * evtCuts) const ;

  TStopwatch        fTimer         ;   //! stopwatch
  TString           fESDTreeName   ;   //! name of the ESD TTree
  UShort_t          fnumberOfTasks ;   //! number of tasks
  AliAnalysisTask ** fTaskList      ;  //! list of tasks
  TClass          ** fTaskInType    ;  //! list of tasks input
  TClass          ** fTaskOuType    ;  //! list of tasks output

  ClassDef(AliAnalysisGoodies, 0); // an analysis utilities class
};
#endif // ALIANALYSISGOODIES_H
