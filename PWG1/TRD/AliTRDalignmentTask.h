#ifndef ALITRDALIGNMENTTASK_H
#define ALITRDALIGNMENTTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD Alignment QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TH1;
class TTree;
class AliTrackPoint;
class AliTrackPointArray;
class AliTRDtrackV1;
class AliTRDalignmentTask : public AliTRDrecoTask
{
public:

  AliTRDalignmentTask();
  AliTRDalignmentTask(char* name);
  virtual ~AliTRDalignmentTask();
  
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *opt);
  TH1*    PlotTrackPoints(const AliTRDtrackV1 *track=NULL);
  //Bool_t  PostProcess() { return kTRUE;}
  
private:
  Bool_t IsIdenticalWithOneOf(AliTrackPoint * const p, AliTrackPointArray *parray, int nmax);
  AliTRDalignmentTask(const AliTRDalignmentTask&);
  AliTRDalignmentTask& operator=(const AliTRDalignmentTask&);

private:
  TTree          *fTree;    //! pointer to the output TTree 
  AliTrackPointArray *fArray; // pointer to the track points
 
  ClassDef(AliTRDalignmentTask, 1) // tracking resolution task
};
#endif
