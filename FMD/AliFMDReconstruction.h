#ifndef AliFMDReconstruction_h
//#define AliFMDReconstruction_h 1
#define AliFMDReconstruction_h

// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
#include "AliFMD.h"
#include "AliDetector.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliFMDReconstruction: public TTask 
{
 public:
  AliFMDReconstruction() ; 
  AliFMDReconstruction(char* HeaderFile,char *SdigitsFile = 0) ; 
  virtual ~AliFMDReconstruction();
  char *GetReconstParticlesFile()const{return (char*) fReconstParticlesFile.Data();}  
  virtual void  Exec(TClonesArray *fReconParticles,Option_t *option); 
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  Stat_t GetNEvents(){return fNevents;}
  void SetReconstParticlesFile(char * file ) ;
  virtual void Print(Option_t* option) const ;
  //TClonesArray *ReconParticles() const {return fReconParticles;}   

 private:
  Int_t   fNevents ;                         // Number of events to digitize
  TString fReconstParticlesFile;    //output file 
  TString fHeadersFile ;                     //input file
  //  TClonesArray *fReconParticles;
  ClassDef(AliFMDReconstruction,1) 
};
#endif









