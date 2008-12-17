#ifndef ALITRDTRACKINGEFFICIENCYCOMBINED_H
#define ALITRDTRACKINGEFFICIENCYCOMBINED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingEfficiencyCombined.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDtrackingEfficiencyCombined : public AliTRDrecoTask{
public:
  AliTRDtrackingEfficiencyCombined();
  virtual ~AliTRDtrackingEfficiencyCombined(){;}
  
  void        CreateOutputObjects();
  void        Exec(Option_t *);
  void        Terminate(Option_t *);
  
  Bool_t      PostProcess();
  TObjArray*  Histos();
  Bool_t      GetRefFigure(Int_t ifig);
    
private:
  enum{
    kEfficiencyHistogram = 0,
    kContaminationHistogram = 1,
    kEfficiencySpeciesHistogram = 2,
    kContaminationSpeciesHistogram = 7,
    kEfficiencyNoPID = 12,
    kContaminationNoPID = 13
  };
  typedef enum{
    kAccepted = 0,
    kRejected = 1,
    kContamination = 2
  } FillingMode_t;
  AliTRDtrackingEfficiencyCombined(const AliTRDtrackingEfficiencyCombined &);
  AliTRDtrackingEfficiencyCombined& operator=(const AliTRDtrackingEfficiencyCombined &);
  
  void    FillHistograms(Int_t ntracks, Int_t *indices, FillingMode_t mode);
  void    FillStreamTrackWOMC(AliTRDtrackInfo *trkInf);

  Bool_t  IsFindable(AliTRDtrackInfo *trkInf);
  Bool_t  IsRegistered(AliTRDtrackInfo *trkInf, Int_t *indices, Int_t nTracks);
    
  ClassDef(AliTRDtrackingEfficiencyCombined, 1); // Combined tracking efficiency
};
    
#endif
