#ifndef ALIANALYSISTASKPTMOTHFROMPTDAUGH_H
#define ALIANALYSISTASKPTMOTHFROMPTDAUGH_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//             Class AnalysisTaskAliPtMothFromPtDaugh                    //
//   AnalysisTaskSE used for the reconstruction of mothers particles     //
//   spectra (pT and pTMin) starting from the pT-spectra of              //
//   daughters particles.                                                //
//                                                                       //
//   Contact: Giuseppe.Bruno@ba.infn.it & Fiorella.Fionda@ba.infn.it     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

class AliPtMothFromPtDaugh;
class TNtuple;
class TH1F;
class TList;

class AliAnalysisTaskPtMothFromPtDaugh : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskPtMothFromPtDaugh();
  AliAnalysisTaskPtMothFromPtDaugh(Bool_t IsNtuplaCreated);
  virtual ~AliAnalysisTaskPtMothFromPtDaugh();

  virtual void  UserExec(Option_t *option);
  virtual void  UserCreateOutputObjects();
  virtual void  Terminate(Option_t *option); 
 
  void SetPtMothFromPtDaugh(AliPtMothFromPtDaugh * const ptExtr)
              { fPtMothDaugh = ptExtr; }             // set AliPtMothFromPtDaugh object
  void SetReadKineFromNtupla(Bool_t ReadKinematic)
              {fReadKineFromNtupla = ReadKinematic;} // set flag to read kinematics from Ntupla  
  void SetNtuplaFileName(char *fileNtuplaName)
              {fFileNtuplaName=fileNtuplaName;}      // set file name from which Ntupla is read
  TNtuple *ReadNtuplaFromFile(char * inFileName);    // get the Ntupla from the file 
  AliPtMothFromPtDaugh *GetPtMothFromPtDaugh(){return fPtMothDaugh;}   

private:
  
  AliPtMothFromPtDaugh    *fPtMothDaugh;        //Pointer to AliPtMothFromPtDaugh object   
  TNtuple                 *fDecayKine;          //Ntupla to store kinematic information of Decay (optional output)
  Bool_t                   fReadKineFromNtupla; //kTRUE: read kinematics from Ntupla
                                                //kFALSE: loops on events to evaluate Ntupla
  char                    *fFileNtuplaName;     //file name from which Ntupla is read 
  TList                   *fList;               //List of mothers histograms (standard output)

  AliAnalysisTaskPtMothFromPtDaugh(const AliAnalysisTaskPtMothFromPtDaugh &c);
  AliAnalysisTaskPtMothFromPtDaugh& operator= (const AliAnalysisTaskPtMothFromPtDaugh &c);
  
  ClassDef(AliAnalysisTaskPtMothFromPtDaugh,1);  // task for analysis of mother pt spectrum from daughter pt spectrum
};
#endif
