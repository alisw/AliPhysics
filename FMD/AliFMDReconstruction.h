//   Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
//  See cxx source for full Copyright notice                               
//  AliFMDReconstruction.h 
//  Task Class for making TreeR in FMD                        
//-- Authors: Evgeny Karpechev (INR) and Alla Maevskaia (INR)


#ifndef AliFMDReconstruction_h
#define AliFMDReconstruction_h

#include "TTask.h"
#include "TString.h"
#include "AliFMD.h"

class AliFMDReconstruction: public TTask 
{
 public:
  AliFMDReconstruction() ; 
  AliFMDReconstruction(char* HeaderFile,char *SdigitsFile = 0) ; 
  virtual ~AliFMDReconstruction();
  char *GetReconstParticlesFile()const{return (char*) fReconstParticlesFile.Data();}  
  virtual void  Exec(Option_t *option); 
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  Stat_t GetNEvents(){return fNevents;}
  void SetReconstParticlesFile(char * file ) ;
  virtual void Print(Option_t* option) const ;   

 private:
  Int_t   fNevents ;                         // Number of events
  TString fReconstParticlesFile;             //output file 
  TString fHeadersFile ;                     //input file
  ClassDef(AliFMDReconstruction,2) 
};
#endif









