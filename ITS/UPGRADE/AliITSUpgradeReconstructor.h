#ifndef ALIITSUPGRADERECONSTRUCTOR_H
#define ALIITSUPGRADERECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//.
// ITS upgrade  base class to reconstruct an event
//.
#include "AliITSReconstructor.h"        //base class
#include "AliITSDigitUpgrade.h"           //Dig2Clu(), UseDig()
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
class AliITSUpgradeReconstructor: public AliITSReconstructor
{ 
 public:
  AliITSUpgradeReconstructor();               
  virtual ~AliITSUpgradeReconstructor();                    //dtor  
  virtual void Init();

  virtual void SetTreeAddressD(TTree* const treeD); 
  void SetTreeAddressR(TTree* const treeR);
  void AddRecPoint(const AliITSRecPoint &p);
  TBranch*  MakeBranchInTree(TTree* const tree,const char* name, const char *classname,void* address,Int_t size,Int_t splitlevel);
		
  virtual void ResetDigits(); 
  virtual void ResetDigits(Int_t branch);
  void MakeBranchRF(TTree *treeR){MakeBranchR(treeR,"Fast");} 
  void DigitsToRecPoints(TTree *treeD,TTree *treeR);
  void MakeBranchR(TTree *treeR,Option_t *opt=" ");
  void ResetRecPoints(){if(fRecPoints) fRecPoints->Clear();fNRecPoints = 0;};
  virtual void MakeBranch(TTree *tree,Option_t *opt);
  virtual AliTracker*  CreateTracker() const;
 private:
  AliITSUpgradeReconstructor(const AliITSUpgradeReconstructor&);              //Not implemented
  AliITSUpgradeReconstructor &operator=(const AliITSUpgradeReconstructor&);   //Not implemented
  TClonesArray *fRecPoints;  //! List of reconstructed points
  Int_t         fNRecPoints; // Number of rec points
  TObjArray    *fDigits;     


  //TString fmemberToWriteLater;
  //  
  ClassDef(AliITSUpgradeReconstructor, 1)        // class for the ITS reconstruction
    };


#endif


