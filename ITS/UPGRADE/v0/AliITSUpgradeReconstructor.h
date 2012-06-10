#ifndef ALIITSUPGRADERECONSTRUCTOR_H
#define ALIITSUPGRADERECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//.
// ITS upgrade  base class to reconstruct an event
//.
#include "AliITSReconstructor.h"        //base class
#include "AliITSDigitUpgrade.h"           
#include "AliITSsegmentationUpgrade.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliITSUpgradeClusterFinder.h"
#include "AliITSRecoParam.h"
#include <TMatrixF.h>                //UseDig()
#include <TFile.h>
#include <TNtupleD.h>
#include <TClonesArray.h>            //UseDig()
#include <TObjArray.h>               //SigConv()
class AliRawReader;                  //Reconstruct() with raw data   
class AliITSRecPoint;
class AliITSUpgradeReconstructor: public AliReconstructor
{ 
 public:
  AliITSUpgradeReconstructor();               
  virtual ~AliITSUpgradeReconstructor();                    //dtor  
  virtual void Init();

  virtual void ResetDigits(); 
  virtual void ResetDigits(Int_t branch);
  virtual AliTracker*  CreateTracker() const;

  virtual void  Reconstruct(TTree* digitsTree, TTree* clustersTree) const; 
  virtual void  Reconstruct(AliRawReader * /*rawdata*/, TTree* /*clustersTree*/) const {AliInfo("Not implemented");} 

  static const AliITSRecoParam* GetRecoParam() { return dynamic_cast<const AliITSRecoParam*>(AliReconstructor::GetRecoParam(0)); }

 private:
  AliITSUpgradeReconstructor(const AliITSUpgradeReconstructor&);              //Not implemented
  AliITSUpgradeReconstructor &operator=(const AliITSUpgradeReconstructor&);   //Not implemented
  TObjArray  *fDigits;     
  Int_t      fNlayers;

 
  ClassDef(AliITSUpgradeReconstructor, 1)        // class for the ITS reconstruction
    };


#endif


