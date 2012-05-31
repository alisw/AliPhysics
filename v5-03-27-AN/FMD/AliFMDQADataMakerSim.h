#ifndef ALIFMDQADATAMAKERSIM_H
#define ALIFMDQADATAMAKERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */

#include "AliQADataMakerSim.h"
// #include "TClonesArray.h"
class TH1F ; 
class TH1I ; 
class TList ; 
class TClonesArray;
//_____________________________________________________________________
// This class implements the AliQADataMakerSim for the FMD. Some
// functions are not implemented yet.
// Author : Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//_____________________________________________________________________



class AliFMDQADataMakerSim: public AliQADataMakerSim {

 public:
  AliFMDQADataMakerSim() ;          // ctor
  AliFMDQADataMakerSim(const AliFMDQADataMakerSim& qadm) ;   
  AliFMDQADataMakerSim& operator = (const AliFMDQADataMakerSim& qadm) ;
  virtual ~AliFMDQADataMakerSim();  // dtor
  
 private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list);
  virtual void   InitHits(); 
  virtual void   InitDigits(); 
  // virtual void   InitRaws() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeHits() ;
  virtual void   MakeHits(TTree* hitTree) ;
  virtual void   MakeDigits() ; 
  virtual void   MakeDigits(TTree* digitTree) ; 
  // virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   MakeSDigits() ; 
  virtual void   MakeSDigits(TTree * sigitTree) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliFMDQADataMakerSim,0)  // description 
    };

#endif // AliFMDQADataMakerSim_H
//____________________________________________________________________
//
// Local Variables: 
//  mode: c++
// End:
//
