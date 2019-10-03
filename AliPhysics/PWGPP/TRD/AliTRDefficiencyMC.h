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
  enum ETRDefficiencyMC{
    kFindable= 0
   ,kNoChmb
   ,kCurved
   ,kPCut
   ,kPhiCut
   ,kThtCut
   ,kLayer
   ,kPrimary
  };
  enum ETRDefficiencyMCstatus{
    kAccept = 0
   ,kMiss
   ,kFake
  };
  AliTRDefficiencyMC();
  AliTRDefficiencyMC(char* name);
  virtual ~AliTRDefficiencyMC(){;}
  
  void        UserCreateOutputObjects();
  void        UserExec(Option_t *);
  
  Bool_t      PostProcess();
  TObjArray*  Histos();
  Bool_t      GetRefFigure(Int_t ifig);
    
private:
  enum{
    kEfficiencyHistogram = 0
   ,kContaminationHistogram = 1
   ,kEfficiencySpeciesHistogram = 2
  };
  AliTRDefficiencyMC(const AliTRDefficiencyMC &);
  AliTRDefficiencyMC& operator=(const AliTRDefficiencyMC &);
  
  void    FillHistograms(Int_t ntracks, Int_t *indices, ETRDefficiencyMCstatus mode);
  void    FillStreamTrackWOMC(AliTRDtrackInfo * const trkInf);

  Int_t   IsFindableNot(AliTRDtrackInfo * const trkInf);
  Int_t   IsRegistered(const AliTRDtrackInfo * const trkInf, Int_t *indices, Int_t nTracks);
    
  static Float_t      fgPCut;   // lower momentum cut
  static Float_t      fgPhiCut; // higher phi cut
  static Float_t      fgThtCut; // higher theta cut


  ClassDef(AliTRDefficiencyMC, 2); // Combined tracking efficiency
};
    
#endif
