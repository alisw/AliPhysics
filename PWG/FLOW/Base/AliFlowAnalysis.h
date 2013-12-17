/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWANALYSISBASE_H
#define ALIFLOWANALYSISBASE_H

class AliFlowEventSimple;
class AliFlowEventSimpleCuts;
class TList;
class TDirectoryFile;
#include "TNamed.h"

/////////////////////////////////////////////////////////////////////////////
// Description: base class for flow analysis classes
// origin:      Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
/////////////////////////////////////////////////////////////////////////////
 
class AliFlowAnalysis : public TNamed {
 public:
 
   AliFlowAnalysis(const char* name="AliFlowAnalysis");            //default constructor
   virtual  ~AliFlowAnalysis();  //destructor
 
   virtual void Init() {}                                       //Define output objects
   virtual void ProcessEvent(AliFlowEventSimple* /*anEvent*/);  //Main routine executed by the framework
   virtual void Make(AliFlowEventSimple* /*anEvent*/) {}        //Main routine to be implemened by user
   virtual void Finish() {}                                           //Fill results
   virtual void GetOutputHistograms(TList* /* outputListHistos */) {} //Copy output objects from TList
   virtual void WriteHistograms(TDirectoryFile* /* outputFileName */) const {} //writes histograms locally (for OnTheFly)
   
   virtual TList*    GetHistList()      const {return NULL;}
   virtual void      SetHistList(TList*) {}
 
   void SetPOItype(Int_t i) {fPOItype=i;}
   Int_t GetPOItype() const {return fPOItype;}

   void SetEventCuts(AliFlowEventSimpleCuts* cuts) {fEventCuts=cuts;}
   AliFlowEventSimpleCuts* GetEventCuts() const {return fEventCuts;}

 protected:
   AliFlowEventSimpleCuts* fEventCuts;  //some analysis level event cuts
   Int_t fPOItype;                       //which POI type are we processing in this analysis?
 
 private:
   AliFlowAnalysis(const AliFlowAnalysis& anAnalysis);            //copy constructor
   AliFlowAnalysis& operator=(const AliFlowAnalysis& anAnalysis); //assignment operator

   ClassDef(AliFlowAnalysis,1)  // class version
};
 
#endif
