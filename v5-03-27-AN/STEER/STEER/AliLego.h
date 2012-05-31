#ifndef ALILEGO_H
#define ALILEGO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
class TH2F;
class AliLegoGenerator;
class TClonesArray;

class AliLego : public TNamed  {

public:
  AliLego();
  AliLego(const char *title, Int_t ntheta,Float_t themin,
	  Float_t themax, Int_t nphi, Float_t phimin,
	  Float_t phimax,Float_t rmin,Float_t rmax,Float_t zmax);
  AliLego(const char *title, AliLegoGenerator* generator);
  AliLego(const AliLego &lego);
  virtual ~AliLego();
  virtual void  StepManager();
  virtual void  BeginEvent();
  virtual void  FinishEvent();
  virtual void  FinishRun();
  virtual AliLego &operator=(const AliLego &lego) 
  {lego.Copy(*this);return(*this);}

private:
  void Copy(TObject &lego) const;
  void DumpVolumes();
  
private:
   AliLegoGenerator *fGener;     //Lego generator
   Float_t    fTotRadl;          //!Total Radiation length
   Float_t    fTotAbso;          //!Total absorption length
   Float_t    fTotGcm2;          //!Total G/CM2 traversed
   TH2F      *fHistRadl;         //Radiation length map 
   TH2F      *fHistAbso;         //Interaction length map
   TH2F      *fHistGcm2;         //g/cm2 length map
   TH2F      *fHistReta;         //Radiation length map as a function of eta
   TClonesArray *fVolumesFwd;    //!Volume sequence forward
   TClonesArray *fVolumesBwd;    //!Volume sequence backward   
   Int_t      fStepBack;         //!Flag for backstepping
   Int_t      fStepsBackward;    //!Counts steps forward
   Int_t      fStepsForward;     //!Counts steps backward
   Int_t      fErrorCondition;   //!Error condition flag
   Int_t      fDebug;            // Debug Flag
   Bool_t     fStopped;          //!Scoring has been stopped 
   
  ClassDef(AliLego,1) //Utility class to compute and draw Radiation Length Map

};


#endif











