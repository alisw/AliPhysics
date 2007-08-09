/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD cluster finder                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h"
#include "AliLog.h"

#include "AliTRDclusterizerV1HLT.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDpadPlane.h"
#include "AliTRDrawData.h"
#include "AliTRDcalibDB.h"
#include "AliTRDSimParam.h"
#include "AliTRDRecParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcluster.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"

ClassImp(AliTRDclusterizerV1HLT)

//_____________________________________________________________________________
AliTRDclusterizerV1HLT::AliTRDclusterizerV1HLT() 
  : AliTRDclusterizerV1(),
    fTreeCreatedHere(kFALSE),
    fNclusters(-1),
    fRawDataVersion(2),
    fRawDataSource(0)
{
  //
  // AliTRDclusterizerV1HLT default constructor
  //
}

//_____________________________________________________________________________
AliTRDclusterizerV1HLT::AliTRDclusterizerV1HLT(const Text_t *name, const Text_t *title) 
  : AliTRDclusterizerV1(name,title),
    fTreeCreatedHere(kFALSE),
    fNclusters(-1),
    fRawDataVersion(2),
    fRawDataSource(0)
{
  //
  // AliTRDclusterizerV1HLT constructor
  //
}

//_____________________________________________________________________________
AliTRDclusterizerV1HLT::AliTRDclusterizerV1HLT(const AliTRDclusterizerV1HLT &c)
  : AliTRDclusterizerV1(c),
    fTreeCreatedHere(kFALSE),
    fNclusters(-1),
    fRawDataVersion(2),
    fRawDataSource(0)
{
  //
  // AliTRDclusterizerV1HLT copy constructor
  //
}

//_____________________________________________________________________________
AliTRDclusterizerV1HLT::~AliTRDclusterizerV1HLT()
{
  //
  // AliTRDclusterizerV1HLT destructor
  //
  if (fTreeCreatedHere == kTRUE)
    delete fClusterTree;
  
  delete fRawDataSource;
  
}

//_____________________________________________________________________________
AliTRDclusterizerV1HLT &AliTRDclusterizerV1HLT::operator=(const AliTRDclusterizerV1HLT &c)
{
  //
  // Assignment operator
  //

  this->fRawDataSource = 0;
  if (this != &c) ((AliTRDclusterizerV1HLT &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizerV1HLT::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliFatal("Not implemented");

//   ((AliTRDclusterizerV1HLT &) c).fDigitsManager = 0;
//   ((AliTRDclusterizerV1HLT &) c).fTreeCreatedHere = kFALSE;
//   AliTRDclusterizerV1::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1HLT::ReadDigits(AliRawReaderMemory *rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  if (fRawDataSource == 0)
    fRawDataSource = new AliTRDrawData;

  //PH  fRawDataSource->SetRawVersion(fRawDataVersion);
  fDigitsManager = fRawDataSource->Raw2Digits((AliRawReader*)rawReader);
  //AliInfo(Form("Digits manager at 0x%x", fDigitsManager));
  AliDebug(1, Form("Digits manager at 0x%x", fDigitsManager));
  if (fDigitsManager)
    return kTRUE;
  else
    return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1HLT::InitClusterTree()
{
  //
  // This has to be called on HLT - creation of the cluster Tree used by the offline clusterizer (base class)
  //
  Bool_t kReturn = kFALSE;
  if (fClusterTree == 0)
    {
      fClusterTree = new TTree("TRDclusterTree", "TRDclusterTree");
      fTreeCreatedHere = kTRUE;
    }

  if (fClusterTree != 0)
    kReturn = kTRUE;    

  return kReturn;
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1HLT::InsertClusters(TObjArray *tobjarr, Int_t idet)
{
  //
  // Fill the tree with clusters - from a different component for example
  //

  //clear the current
  ResetRecPoints();
  delete fRecPoints;
  fRecPoints = 0;

  //set the pointer used in WriteClusters
  fRecPoints = tobjarr;  
  Bool_t kRet = kFALSE;
  if (InitClusterTree())
    kRet = WriteClusters(idet);  

  fRecPoints = 0;

  return kRet;
}
//_____________________________________________________________________________
Int_t AliTRDclusterizerV1HLT::GetNclusters()
{
  //
  // Returns the number of clusters AliTRDclusterizerV1HLT::fNclusters
  // Count them first if fNclusters < 0
  //

  if (fNclusters < 0)
    {
      CountClusters();
    }
  return fNclusters;
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1HLT::ResetTree()
{
  //
  // Recreate the cluster tree
  // 

  //   if (fClusterTree != 0)
  //     fClusterTree->Reset();
  // well we'd better delete the whole tree and branches
  delete fClusterTree;
  fClusterTree = NULL;
  fClusterTree = new TTree("TRDclusterTree", "TRDclusterTree");
  if (fClusterTree)
    {
      fTreeCreatedHere = kTRUE;
      fNclusters = -1;
      //AliInfo("Tree Reset Successful");
      AliDebug(1,Form("Tree Reset Successful"));
    }
  else
    {
      fTreeCreatedHere = kFALSE;
      AliError("Reset Tree Error!\n");
    }
  
  return fTreeCreatedHere;
}

//_____________________________________________________________________________
Int_t AliTRDclusterizerV1HLT::CountClusters()
{
  //
  // Count the clusters - runs over the cluster tree
  //

  fNclusters = -1;
  if (fClusterTree == 0)
    {
      AliError("No tree to count clusters!\n");
      return -1;
    }
  TList *lt = (TList*)fClusterTree->GetListOfBranches();
  TIter it(lt);
  it.Reset();
  TBranch *tb = 0;
  Int_t icount = 0;
  while ((tb = (TBranch*)it.Next()) != 0)
    {
      TObjArray *clusters = 0;
      tb->SetAddress(&clusters);
      for (Int_t icb = 0; icb < tb->GetEntries(); icb++)
	{
	  tb->GetEntry(icb);
	  icount += clusters->GetEntries();
	}
    }
  fNclusters = icount;
  //AliInfo(Form("Recounted clusters %d", fNclusters));
  AliDebug(2, Form("Recounted clusters %d", fNclusters));

  return fNclusters;
}
