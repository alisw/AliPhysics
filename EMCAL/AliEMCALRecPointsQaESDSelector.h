#ifndef ALIEMCALRECPOINTSQAESDSELECTOR_H
#define ALIEMCALRECPOINTSQAESDSELECTOR_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
 
//*--  Authors: Aleksei Pavlinov (WSU)

#include "AliSelector.h"

class AliEMCALFolder;

class TList;
class TBrowser;
class TChain;
class TObjectSet;

class AliEMCALRecPointsQaESDSelector :  public AliSelector {
  public:
    AliEMCALRecPointsQaESDSelector();
    virtual ~AliEMCALRecPointsQaESDSelector();

    virtual void    Begin(TTree* tree);
    virtual void    Init(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();
  //
  void InitStructure(Int_t it);
  static TList *DefineHistsOfRP(const char *name="RP");
  // static TList *DefineHistsOfTowers(const char *name="towers");
  // 
  void FitEffMassHist(); //*MENU*  
  void PrintInfo();      // *MENU*
  //
  void    SetChain(TChain *chain)  {fChain = chain;}
  TChain* GetChain()               {return fChain;}

  AliEMCALFolder* CreateEmcalFolder(const Int_t it);
  void SetEmcalFolder(AliEMCALFolder* folder); 
  void SetEmcalOldFolder(AliEMCALFolder* folder); 
  AliEMCALFolder* GetEmcalFolder() {return fEMCAL;}
  AliEMCALFolder* GetEmcalOldFolder() {return fEMCALOld;}
  AliEMCALFolder* GetEmcalOldFolder(const Int_t nsm);
  // 
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const;

  //
  //// Pictures staf - Jun 26, 2007
  //
  void ReadAllEmcalFolders();
  void PictVsIterNumber(const Int_t ind=0, const Int_t nsm=0);

 protected:
  //
  TChain* fChain; //! chain if ESD files
  TList* fLofHistsPC; //! list of histograms of pseudo clusters 
  TList* fLofHistsRP; //! list of histograms of rec.points 
  //
  TObjectSet*     fEmcalPool;
  AliEMCALFolder* fEMCAL;    //! current  EMCAL object
  AliEMCALFolder* fEMCALOld; //! previous EMCAL object

 private:
  AliEMCALRecPointsQaESDSelector(const AliEMCALRecPointsQaESDSelector&);
  AliEMCALRecPointsQaESDSelector& operator=(const AliEMCALRecPointsQaESDSelector&);

  ClassDef(AliEMCALRecPointsQaESDSelector, 1);
};

#endif
