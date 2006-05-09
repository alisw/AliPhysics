// $Header$

#ifndef REVE_VSD_H
#define REVE_VSD_H

#include "Reve.h"
#include "PODs.h"
#include <TTree.h>

namespace Reve {

class VSDTree : public TTree
{
public:
  ClassDef(VSDTree, 1);
};

class VSD : public TObject
{
protected:
  Int_t        fBuffSize;

  TFile*       mFile;        //!
  TDirectory*  mDirectory;   //!

public:
  TTree*       mTreeK;       //! X{g}
  //TTree*       mTreeTR;      //! X{g}
  TTree*       mTreeH;       //! X{g}
  TTree*       mTreeC;       //! X{g}
  TTree*       mTreeR;       //! X{g}
  TTree*       mTreeKK;      //! X{g}
  TTree*       mTreeV0;      //! X{g}
  TTree*       mTreeGI;      //! X{g}

  MCTrack      mK,  *mpK;    //!
  //MCTrackRef   mTR, *mpTR;   //!
  Hit          mH,  *mpH;    //!
  Cluster      mC,  *mpC;    //!
  RecTrack     mR,  *mpR;    //!
  RecKink      mKK, *mpKK;   //!
  RecV0        mV0, *mpV0;   //!
  GenInfo      mGI, *mpGI;   //!

public:
  VSD();
  VSD(const Text_t* name, const Text_t* title="");
  virtual void InitTreeVars();

  virtual void SetDirectory(TDirectory* dir);

  virtual void CreateTrees();
  virtual void DeleteTrees();

  virtual void CreateBranches();
  virtual void SetBranchAddresses();

  virtual void WriteTrees();
  virtual void LoadTrees();

  virtual void LoadVSD(const Text_t* vsd_file_name,
		       const Text_t* dir_name="Event0");

  ClassDef(VSD, 1);
}; // endclass VSD

}

#endif
