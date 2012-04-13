#ifndef ALITRDEFFICIENCY_H
#define ALITRDEFFICIENCY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDefficiency.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDtrackV1;
class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDefficiency : public AliTRDrecoTask
{
public:
  AliTRDefficiency();
  AliTRDefficiency(char* name);
  virtual       ~AliTRDefficiency();
//  void        UserCreateOutputObjects();
  void          LocalUserExec(Option_t *);
  Bool_t        GetRefFigure(Int_t ifig);
  static Int_t  GetPtBin(Float_t pt);
  static Int_t  GetPtBinSignificant(Float_t pt);
  TObjArray*    Histos();
  TH1*          PlotBasicEff(const AliTRDtrackV1 *t=NULL);
//  TH1*          PlotMC(const AliTRDtrackV1 *t=NULL);
  void          MakeSummary();
  Bool_t        PostProcess();
  TObjArray*    Results() const {return fProj;}
protected:
  Bool_t        MakeProjectionBasicEff();

private:
  AliTRDefficiency(const AliTRDefficiency&);
  AliTRDefficiency& operator=(const AliTRDefficiency&);

  TClonesArray     *fMissed;          // Missed ?
  TObjArray        *fProj;            //! result holder - sigma values

  ClassDef(AliTRDefficiency, 2) // TRD tracking efficiency
};

#endif

