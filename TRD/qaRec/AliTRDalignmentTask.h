#ifndef AliTRDALIGNMENTTASK_H
#define AliTRDALIGNMENTTASK_H
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
  virtual ~AliTRDalignmentTask();
  
  void    CreateOutputObjects();
  void    Exec(Option_t *opt);
  void    GetRefFigure(Int_t ifig);
  TObjArray*  Histos(); 
  TH1*    PlotTrackPoints(const AliTRDtrackV1 *track=0x0);
  Bool_t  PostProcess(){return kTRUE;}
  void    Terminate(Option_t *);
  
private:
  Bool_t IsIdenticalWithOneOf(AliTrackPoint *p, AliTrackPointArray *parray, int nmax);
  AliTRDalignmentTask(const AliTRDalignmentTask&);
  AliTRDalignmentTask& operator=(const AliTRDalignmentTask&);

private:
  TTree          *fTree;    //! pointer to the output TTree 
  AliTrackPointArray *fArray; // pointer to the track points
 
  ClassDef(AliTRDalignmentTask, 1) // tracking resolution task
};
#endif
