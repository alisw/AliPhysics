#ifndef ALITRDEFFICIENCYMC_H
#define ALITRDEFFICIENCYMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDefficiencyMC.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDefficiencyMC : public AliTRDrecoTask{
public:
  AliTRDefficiencyMC();
  virtual ~AliTRDefficiencyMC(){;}
  
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
  AliTRDefficiencyMC(const AliTRDefficiencyMC &);
  AliTRDefficiencyMC& operator=(const AliTRDefficiencyMC &);
  
  void    FillHistograms(Int_t ntracks, Int_t *indices, FillingMode_t mode);
  void    FillStreamTrackWOMC(AliTRDtrackInfo *trkInf);

  Bool_t  IsFindable(AliTRDtrackInfo *trkInf);
  Bool_t  IsRegistered(AliTRDtrackInfo *trkInf, Int_t *indices, Int_t nTracks);
    
  ClassDef(AliTRDefficiencyMC, 1); // Combined tracking efficiency
};
    
#endif
