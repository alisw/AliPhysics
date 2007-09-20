
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD cluster finder                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliAlignObj.h"

#include "AliTRDclusterizer.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDrawData.h"
#include "AliTRDcalibDB.h"
#include "AliTRDSimParam.h"
#include "AliTRDRecParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDtransform.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDRawStream.h"
#include "AliTRDRawStreamV2.h"
#include "AliTRDfeeParam.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer()
  :TNamed()
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fDigitsManager(NULL)
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(NULL)
{
  //
  // AliTRDclusterizer default constructor
  //

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fDigitsManager(new AliTRDdigitsManager())
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(new AliTRDtransform(0))
{
  //
  // AliTRDclusterizer constructor
  //

  fDigitsManager->CreateArrays();

  fRawVersion = AliTRDfeeParam::Instance()->GetRAWversion();

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const AliTRDclusterizer &c)
  :TNamed(c)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
  ,fDigitsManager(NULL)
  ,fAddLabels(kTRUE)
  ,fRawVersion(2)
  ,fIndexesOut(NULL)
  ,fIndexesMaxima(NULL)
  ,fTransform(NULL)
{
  //
  // AliTRDclusterizer copy constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizer::~AliTRDclusterizer()
{
  //
  // AliTRDclusterizer destructor
  //

  if (fRecPoints) 
    {
      fRecPoints->Delete();
      delete fRecPoints;
    }

  if (fDigitsManager) 
    {
      delete fDigitsManager;
      fDigitsManager = NULL;
    }

  if (fIndexesOut)
    {
      delete fIndexesOut;
      fIndexesOut    = NULL;
    }

  if (fIndexesMaxima)
    {
      delete fIndexesMaxima;
      fIndexesMaxima = NULL;
    }

  if (fTransform)
    {
      delete fTransform;
      fTransform     = NULL;
    }

}

//_____________________________________________________________________________
AliTRDclusterizer &AliTRDclusterizer::operator=(const AliTRDclusterizer &c)
{
  //
  // Assignment operator
  //

  if (this != &c) 
    {
      ((AliTRDclusterizer &) c).Copy(*this);
    }
  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizer::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDclusterizer &) c).fClusterTree   = NULL;
  ((AliTRDclusterizer &) c).fRecPoints     = NULL;  
  ((AliTRDclusterizer &) c).fDigitsManager = NULL;
  ((AliTRDclusterizer &) c).fAddLabels     = fAddLabels;
  ((AliTRDclusterizer &) c).fRawVersion    = fRawVersion;
  ((AliTRDclusterizer &) c).fIndexesOut    = NULL;
  ((AliTRDclusterizer &) c).fIndexesMaxima = NULL;
  ((AliTRDclusterizer &) c).fTransform     = NULL;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file. Output and input are in the same file
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader         = AliRunLoader::GetRunLoader(evfoldname);

  if (!fRunLoader) {
    fRunLoader = AliRunLoader::Open(name);
  }

  if (!fRunLoader) {
    AliError(Form("Can not open session for file %s.",name));
    return kFALSE;
  }

  OpenInput(nEvent);
  OpenOutput();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput()
{
  //
  // Open the output file
  //

  TObjArray *ioArray = 0;

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  loader->MakeTree("R");

  fClusterTree = loader->TreeR();
  fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput(TTree *clusterTree)
{
  //
  // Connect the output tree
  //

  TObjArray *ioArray = 0;

  fClusterTree = clusterTree;
  fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenInput(Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the digits-tree
  //

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(nEvent);
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteClusters(Int_t det)
{
  //
  // Fills TRDcluster branch in the tree with the clusters 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det <                      -1) || 
      (det >= AliTRDgeometry::Ndet())) {
    AliError(Form("Unexpected detector index %d.\n",det));
    return kFALSE;
  }
 
  TBranch *branch = fClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    TObjArray *ioArray = 0;
    branch = fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);
  }

  if ((det >=                      0) && 
      (det <  AliTRDgeometry::Ndet())) {

    Int_t nRecPoints = RecPoints()->GetEntriesFast();
    TObjArray *detRecPoints = new TObjArray(400);

    for (Int_t i = 0; i < nRecPoints; i++) {
      AliTRDcluster *c = (AliTRDcluster *) RecPoints()->UncheckedAt(i);
      if (det == c->GetDetector()) {
        detRecPoints->AddLast(c);
      }
      else {
        AliError(Form("Attempt to write a cluster with unexpected detector index: got=%d expected=%d\n"
                     ,c->GetDetector()
                     ,det));
      }
    }

    branch->SetAddress(&detRecPoints);
    fClusterTree->Fill();

    delete detRecPoints;
    
    return kTRUE;

  }

  if (det == -1) {

    AliInfo(Form("Writing the cluster tree %s for event %d."
	        ,fClusterTree->GetName(),fRunLoader->GetEventNumber()));

    if (fRecPoints) {

      branch->SetAddress(&fRecPoints);

      AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
      loader->WriteRecPoints("OVERWRITE");
  
    }
    else {

      AliError("Cluster tree does not exist. Cannot write clusters.\n");
      return kFALSE;

    }

    return kTRUE;  

  }

  AliError(Form("Unexpected detector index %d.\n",det));
 
  return kFALSE;  
  
}

//_____________________________________________________________________________
void AliTRDclusterizer::ResetHelperIndexes(AliTRDSignalIndex *indexesIn)
{
  // 
  // Reset the helper indexes
  //

  if (fIndexesOut)
    {
      // carefull here - we assume that only row number may change - most probable
      if (indexesIn->GetNrow() <= fIndexesOut->GetNrow())
 	fIndexesOut->ResetContent();
      else
	fIndexesOut->ResetContentConditional(indexesIn->GetNrow()
                                           , indexesIn->GetNcol()
                                           , indexesIn->GetNtime());
    }
  else
    {
      fIndexesOut = new AliTRDSignalIndex(indexesIn->GetNrow()
                                        , indexesIn->GetNcol() 
                                        , indexesIn->GetNtime());
    }
  
  if (fIndexesMaxima)
    {
      // carefull here - we assume that only row number may change - most probable
      if (indexesIn->GetNrow() <= fIndexesMaxima->GetNrow()) 
        {
 	  fIndexesMaxima->ResetContent();
	}
      else
	{
	  fIndexesMaxima->ResetContentConditional(indexesIn->GetNrow()
                                                , indexesIn->GetNcol()
                                                , indexesIn->GetNtime());
	}
    }
  else
    {
      fIndexesMaxima = new AliTRDSignalIndex(indexesIn->GetNrow()
                                           , indexesIn->GetNcol()
                                           , indexesIn->GetNtime());
    }

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits()
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fRunLoader) {
    AliError("No run loader available");
    return kFALSE;
  }

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader->TreeD()) {
    loader->LoadDigits();
  }

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits(loader->TreeD()));

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits(TTree *digitsTree)
{
  //
  // Reads the digits arrays from the input tree
  //

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits(digitsTree));

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::ReadDigits(AliRawReader *rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData raw;
  fDigitsManager = raw.Raw2Digits(rawReader);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::MakeClusters()
{
  //
  // Creates clusters from digits
  //

  // Propagate info from the digits manager
  if (fAddLabels == kTRUE)
    {
      fAddLabels = fDigitsManager->UsesDictionaries();
    }

  Bool_t fReturn = kTRUE;
  for (Int_t i = 0; i < AliTRDgeometry::kNdet; i++)
    {

      AliTRDdataArrayI *digitsIn = fDigitsManager->GetDigits(i);      
      // This is to take care of switched off super modules
      if (digitsIn->GetNtime() == 0) 
        {
	  continue;
        }
      digitsIn->Expand();
      AliTRDSignalIndex* indexes = fDigitsManager->GetIndexes(i);
      if (indexes->IsAllocated() == kFALSE)
	{
	  fDigitsManager->BuildIndexes(i);
	}

      Bool_t fR = kFALSE;
      if (indexes->HasEntry())
	{
	  if (fAddLabels)
	    {
	      for (Int_t iDict = 0; iDict < AliTRDdigitsManager::kNDict; iDict++) 
		{
		  AliTRDdataArrayI *tracksIn = 0;
		  tracksIn = fDigitsManager->GetDictionary(i,iDict);
		  tracksIn->Expand();
		}
	    }
	  fR = MakeClusters(i);
	  fReturn = fR && fReturn;
	}

      if (fR == kFALSE)
	{
	  WriteClusters(i);
	  ResetRecPoints();
	}

      // No compress just remove
      fDigitsManager->RemoveDigits(i);
      fDigitsManager->RemoveDictionaries(i);      
      fDigitsManager->ClearIndexes(i);

    }

  return fReturn;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Raw2Clusters(AliRawReader *rawReader)
{
  //
  // Creates clusters from raw data
  //

  AliTRDdataArrayI *digits = 0;
  AliTRDdataArrayI *track0 = 0;
  AliTRDdataArrayI *track1 = 0;
  AliTRDdataArrayI *track2 = 0; 

  AliTRDSignalIndex *indexes = 0;

  // Create the digits manager
  if (!fDigitsManager)
    {
      fDigitsManager = new AliTRDdigitsManager();
      fDigitsManager->CreateArrays();
    }

  AliTRDRawStreamV2 input(rawReader);
  input.SetRawVersion( fRawVersion );
  input.Init();

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));

  // Loop through the digits
  Int_t lastdet = -1;
  Int_t det     =  0;
  Int_t it      =  0;
  while (input.Next()) 
    {

      det = input.GetDet();

      if (det != lastdet) 
	{
	
	  if (lastdet != -1)
	    {
	      digits = fDigitsManager->GetDigits(lastdet);
	      Bool_t iclusterBranch = kFALSE;
	      if (indexes->HasEntry())
		iclusterBranch = MakeClusters(lastdet);
	      if (iclusterBranch == kFALSE)
		{
		  WriteClusters(lastdet);
		  ResetRecPoints();
		}
	    }

	  if (digits)
	    {
	      fDigitsManager->RemoveDigits(lastdet);
	      fDigitsManager->RemoveDictionaries(lastdet);
	      fDigitsManager->ClearIndexes(lastdet);
	    }

	  lastdet = det;

	  // Add a container for the digits of this detector
	  digits = fDigitsManager->GetDigits(det);
	  track0 = fDigitsManager->GetDictionary(det,0);
	  track1 = fDigitsManager->GetDictionary(det,1);
	  track2 = fDigitsManager->GetDictionary(det,2);

	  // Allocate memory space for the digits buffer
	  if (digits->GetNtime() == 0) 
	    {
	      //AliDebug(5, Form("Alloc digits for det %d", det));
	      digits->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track0->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track1->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track2->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	    }
	  
	  indexes = fDigitsManager->GetIndexes(det);
	  indexes->SetSM(input.GetSM());
	  indexes->SetStack(input.GetStack());
	  indexes->SetLayer(input.GetLayer());
	  indexes->SetDetNumber(det);
	  if (indexes->IsAllocated() == kFALSE)
	    {
	      indexes->Allocate(input.GetMaxRow(), input.GetMaxCol(), input.GetNumberOfTimeBins());
	    }

	}
      
      for (it = 0; it < 3; it++)
	{
	  if ( input.GetTimeBin() + it < input.GetNumberOfTimeBins() )
	    {
	      if (input.GetSignals()[it] > 0)
		{
		  digits->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, input.GetSignals()[it]);

		  indexes->AddIndexTBin(input.GetRow(), input.GetCol(),
					input.GetTimeBin() + it);
		  track0->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		  track1->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		  track2->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		}
	    }
	}

    }

  if (lastdet != -1)
    {
      Bool_t iclusterBranch = kFALSE;
      if (indexes->HasEntry()) 
        {
	  iclusterBranch = MakeClusters(lastdet);
        }
      if (iclusterBranch == kFALSE)
	{
	  WriteClusters(lastdet);
	  ResetRecPoints();
	}
      //MakeClusters(lastdet);
      if (digits)
	{
	  fDigitsManager->RemoveDigits(lastdet);
	  fDigitsManager->RemoveDictionaries(lastdet);
	  fDigitsManager->ClearIndexes(lastdet);
	}
    }

  delete fDigitsManager;
  fDigitsManager = NULL;
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Raw2ClustersChamber(AliRawReader *rawReader)
{
  //
  // Creates clusters from raw data
  //

  // Create the digits manager
  if (!fDigitsManager)
    {
      fDigitsManager = new AliTRDdigitsManager();
      fDigitsManager->CreateArrays();
    }

  fDigitsManager->SetUseDictionaries(fAddLabels);

  AliTRDRawStreamV2 input(rawReader);
  input.SetRawVersion( fRawVersion );
  input.Init();

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));
  
  Int_t det    = 0;
  while ((det = input.NextChamber(fDigitsManager)) >= 0)
    {
      Bool_t iclusterBranch = kFALSE;
      if (fDigitsManager->GetIndexes(det)->HasEntry())
	{
	  iclusterBranch = MakeClusters(det);
	}
      if (iclusterBranch == kFALSE)
	{
	  WriteClusters(det);
	  ResetRecPoints();
	}
      fDigitsManager->RemoveDigits(det);
      fDigitsManager->RemoveDictionaries(det);      
      fDigitsManager->ClearIndexes(det);
    }

  delete fDigitsManager;
  fDigitsManager = NULL;
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::MakeClusters(Int_t det)
{
  //
  // Generates the cluster.
  //

  // Get the digits
  //   digits should be expanded beforehand!
  //   digitsIn->Expand();
  AliTRDdataArrayI *digitsIn = fDigitsManager->GetDigits(det);      

  // This is to take care of switched off super modules
  if (digitsIn->GetNtime() == 0) 
    {
      return kFALSE;
    }

  AliTRDSignalIndex *indexesIn = fDigitsManager->GetIndexes(det);
  if (indexesIn->IsAllocated() == kFALSE)
    {
      AliError("Indexes do not exist!");
      return kFALSE;      
    }
    
  AliTRDcalibDB  *calibration = AliTRDcalibDB::Instance();
  if (!calibration) 
    {
      AliFatal("No AliTRDcalibDB instance available\n");
      return kFALSE;  
    }

  AliTRDRecParam *recParam    = AliTRDRecParam::Instance();
  if (!recParam) 
    {
      AliError("No AliTRDRecParam instance available\n");
      return kFALSE;  
    }

  // ADC thresholds
  // There is no ADC threshold anymore, and simParam should not be used in clusterizer. KO
  Float_t ADCthreshold   = 0; 

  // Threshold value for the maximum
  Float_t maxThresh      = recParam->GetClusMaxThresh();
  // Threshold value for the digit signal
  Float_t sigThresh      = recParam->GetClusSigThresh();

  // Iteration limit for unfolding procedure
  const Float_t kEpsilon = 0.01;             
  const Int_t   kNclus   = 3;  
  const Int_t   kNsig    = 5;

  Int_t    iUnfold       = 0;  
  Double_t ratioLeft     = 1.0;
  Double_t ratioRight    = 1.0;

  Double_t padSignal[kNsig];   
  Double_t clusterSignal[kNclus];

  Int_t icham = indexesIn->GetChamber();
  Int_t iplan = indexesIn->GetPlane();
  Int_t isect = indexesIn->GetSM();

  // Start clustering in the chamber

  Int_t idet  = AliTRDgeometry::GetDetector(iplan,icham,isect);
  if (idet != det)
    {
      AliError("Strange Detector number Missmatch!");
      return kFALSE;
    }

  // TRD space point transformation
  fTransform->SetDetector(det);

  Int_t    ilayer  = AliGeomManager::kTRD1 + iplan;
  Int_t    imodule = icham + AliTRDgeometry::Ncham() * isect;
  UShort_t volid   = AliGeomManager::LayerToVolUID(ilayer,imodule); 

  Int_t nColMax    = digitsIn->GetNcol();
  Int_t nTimeTotal = digitsIn->GetNtime();

  // Detector wise calibration object for the gain factors
  const AliTRDCalDet *calGainFactorDet      = calibration->GetGainFactorDet();
  // Calibration object with pad wise values for the gain factors
  AliTRDCalROC       *calGainFactorROC      = calibration->GetGainFactorROC(idet);
  // Calibration value for chamber wise gain factor
  Float_t             calGainFactorDetValue = calGainFactorDet->GetValue(idet);

  Int_t nClusters = 0;

  AliTRDdataArrayF *digitsOut = new AliTRDdataArrayF(digitsIn->GetNrow()
                                                    ,digitsIn->GetNcol()
			                            ,digitsIn->GetNtime()); 

  ResetHelperIndexes(indexesIn);

  // Apply the gain and the tail cancelation via digital filter
  TailCancelation(digitsIn
	         ,digitsOut  
	         ,indexesIn
	         ,fIndexesOut
	         ,nTimeTotal
                 ,ADCthreshold
                 ,calGainFactorROC
                 ,calGainFactorDetValue);	
	
  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;
  Int_t iPad = 0;
    
  fIndexesOut->ResetCounters();
  while (fIndexesOut->NextRCTbinIndex(row, col, time))
    {

      Float_t signalM = TMath::Abs(digitsOut->GetDataUnchecked(row,col,time));
	    
      // Look for the maximum
      if (signalM >= maxThresh) 
	{
		
	  if (col + 1 >= nColMax || col-1 < 0)
	    continue;

	  Float_t signalL = TMath::Abs(digitsOut->GetDataUnchecked(row,col+1,time));
	  Float_t signalR = TMath::Abs(digitsOut->GetDataUnchecked(row,col-1,time));

	  if ((TMath::Abs(signalL) <= signalM) && 
	      (TMath::Abs(signalR) <  signalM)) 
	    {
	      if ((TMath::Abs(signalL) >= sigThresh) ||
		  (TMath::Abs(signalR) >= sigThresh)) 
		{
		  // Maximum found, mark the position by a negative signal
		  digitsOut->SetDataUnchecked(row,col,time,-signalM);
		  fIndexesMaxima->AddIndexTBin(row,col,time);
		}
	    }

	}

    }
	       
  // The index to the first cluster of a given ROC
  Int_t firstClusterROC = -1;
  // The number of cluster in a given ROC
  Int_t nClusterROC     =  0;

  // Now check the maxima and calculate the cluster position
  fIndexesMaxima->ResetCounters();
  while (fIndexesMaxima->NextRCTbinIndex(row, col, time)) 
    {

      // Maximum found ?             
      if (digitsOut->GetDataUnchecked(row,col,time) < 0.0) 
        {

	  for (iPad = 0; iPad < kNclus; iPad++) 
            {
	      Int_t iPadCol = col - 1 + iPad;
	      clusterSignal[iPad] = TMath::Abs(digitsOut->GetDataUnchecked(row,iPadCol,time));
	    }

	  // Count the number of pads in the cluster
	  Int_t nPadCount = 0;
	  Int_t ii;
	  // Look to the left
	  ii = 0;
	  while (TMath::Abs(digitsOut->GetDataUnchecked(row,col-ii  ,time)) >= sigThresh) 
            {
	      nPadCount++;
	      ii++;
	      if (col-ii   <        0) break;
	    }
	  // Look to the right
	  ii = 0;
	  while (TMath::Abs(digitsOut->GetDataUnchecked(row,col+ii+1,time)) >= sigThresh) 
            {
	      nPadCount++;
	      ii++;
	      if (col+ii+1 >= nColMax) break;
	    }
	  nClusters++;

	  // Look for 5 pad cluster with minimum in the middle
	  Bool_t fivePadCluster = kFALSE;
	  if (col < (nColMax - 3)) 
            {
	      if (digitsOut->GetDataUnchecked(row,col+2,time) < 0) 
                {
	          fivePadCluster = kTRUE;
	        }
	      if ((fivePadCluster) && (col < (nColMax - 5))) 
                {
	          if (digitsOut->GetDataUnchecked(row,col+4,time) >= sigThresh) 
                    {
	              fivePadCluster = kFALSE;
	            }
	        }
	      if ((fivePadCluster) && (col >             1)) 
                {
	          if (digitsOut->GetDataUnchecked(row,col-2,time) >= sigThresh) 
                    {
	              fivePadCluster = kFALSE;
	            }
	        }
	    }

	  // 5 pad cluster
	  // Modify the signal of the overlapping pad for the left part 
	  // of the cluster which remains from a previous unfolding
	  if (iUnfold) 
            {
	      clusterSignal[0] *= ratioLeft;
	      iUnfold = 0;
	    }

	  // Unfold the 5 pad cluster
	  if (fivePadCluster) 
            {
	      for (iPad = 0; iPad < kNsig; iPad++) 
                {
	          padSignal[iPad] = TMath::Abs(digitsOut->GetDataUnchecked(row
			   					          ,col-1+iPad
								          ,time));
	        }
	      // Unfold the two maxima and set the signal on 
	      // the overlapping pad to the ratio
	      ratioRight        = Unfold(kEpsilon,iplan,padSignal);
	      ratioLeft         = 1.0 - ratioRight; 
	      clusterSignal[2] *= ratioRight;
	      iUnfold = 1;
	    }

	  // The position of the cluster in COL direction relative to the center pad (pad units)
          Double_t clusterPosCol = 0.0;
	  if (recParam->LUTOn()) 
            {
	      // Calculate the position of the cluster by using the
	      // lookup table method
	      clusterPosCol = recParam->LUTposition(iplan
                                                   ,clusterSignal[0]
                                                   ,clusterSignal[1]
 	                                           ,clusterSignal[2]);
	    }
	  else 
            {
	      // Calculate the position of the cluster by using the
	      // center of gravity method
	      for (Int_t i = 0; i < kNsig; i++) 
                {
	          padSignal[i] = 0.0;
	        }
	      padSignal[2] = TMath::Abs(digitsOut->GetDataUnchecked(row,col  ,time)); // Central pad
	      padSignal[1] = TMath::Abs(digitsOut->GetDataUnchecked(row,col-1,time)); // Left    pad
	      padSignal[3] = TMath::Abs(digitsOut->GetDataUnchecked(row,col+1,time)); // Right   pad
	      if ((col >           2) && 
	          (TMath::Abs(digitsOut->GetDataUnchecked(row,col-2,time)) < padSignal[1])) 
                {
	          padSignal[0] = TMath::Abs(digitsOut->GetDataUnchecked(row,col-2,time));
	        }
	      if ((col < nColMax - 3) &&
	          (TMath::Abs(digitsOut->GetDataUnchecked(row,col+2,time)) < padSignal[3])) 
                {
	          padSignal[4] = TMath::Abs(digitsOut->GetDataUnchecked(row,col+2,time));
	        }  
	      clusterPosCol = GetCOG(padSignal);
	    }

	  // Store the amplitudes of the pads in the cluster for later analysis
	  Short_t signals[7] = { 0, 0, 0, 0, 0, 0, 0 };
	  for (Int_t jPad = col-3; jPad <= col+3; jPad++) 
            {
	      if ((jPad <          0) || 
	          (jPad >= nColMax-1)) 
                {
	          continue;
	        }
	      signals[jPad-col+3] = TMath::Nint(TMath::Abs(digitsOut->GetDataUnchecked(row,jPad,time)));
	    }

          // Transform the local cluster coordinates into calibrated 
          // space point positions defined in the local tracking system.
          // Here the calibration for T0, Vdrift and ExB is applied as well.
	  Double_t clusterXYZ[6];
	  clusterXYZ[0] = clusterPosCol;
	  clusterXYZ[1] = clusterSignal[0];
	  clusterXYZ[2] = clusterSignal[1];
	  clusterXYZ[3] = clusterSignal[2];
	  clusterXYZ[4] = 0.0;
	  clusterXYZ[5] = 0.0;
          Int_t    clusterRCT[3];
          clusterRCT[0] = row;
          clusterRCT[1] = col;
          clusterRCT[2] = 0;
	  fTransform->Transform(clusterXYZ,clusterRCT,((UInt_t) time),0);

	  // Add the cluster to the output array
	  // The track indices will be stored later 
          Float_t clusterPos[3];
          clusterPos[0] = clusterXYZ[0];
          clusterPos[1] = clusterXYZ[1];
          clusterPos[2] = clusterXYZ[2];
	  Float_t clusterSig[2];
	  clusterSig[0] = clusterXYZ[4];
	  clusterSig[1] = clusterXYZ[5];
          Double_t clusterCharge  = clusterXYZ[3];
	  Char_t   clusterTimeBin = ((Char_t) clusterRCT[2]);
	  AliTRDcluster *cluster = new AliTRDcluster(idet
						    ,clusterCharge
						    ,clusterPos
						    ,clusterSig
						    ,0x0
						    ,((Char_t) nPadCount)
						    ,signals
						    ,((UChar_t) col)
						    ,((UChar_t) row)
						    ,((UChar_t) time)
						    ,clusterTimeBin
						    ,clusterPosCol
						    ,volid);

	  // Temporarily store the row, column and time bin of the center pad
	  // Used to later on assign the track indices
	  cluster->SetLabel( row,0);
	  cluster->SetLabel( col,1);
	  cluster->SetLabel(time,2);

	  RecPoints()->Add(cluster);

	  // Store the index of the first cluster in the current ROC
	  if (firstClusterROC < 0) 
            {
	      firstClusterROC = RecPoints()->GetEntriesFast() - 1;
	    }

	  // Count the number of cluster in the current ROC
	  nClusterROC++;

        } // if: Maximum found ?

    }

  delete digitsOut;

  if (fAddLabels) 
    {
      AddLabels(idet, firstClusterROC, nClusterROC);
    }

  // Write the cluster and reset the array
  WriteClusters(idet);
  ResetRecPoints();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::AddLabels(Int_t idet, Int_t firstClusterROC, Int_t nClusterROC)
{
  //
  // Add the track indices to the found clusters
  //
  
  const Int_t   kNclus  = 3;  
  const Int_t   kNdict  = AliTRDdigitsManager::kNDict;
  const Int_t   kNtrack = kNdict * kNclus;

  Int_t iClusterROC = 0;

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;
  Int_t iPad = 0;

  // Temporary array to collect the track indices
  Int_t *idxTracks = new Int_t[kNtrack*nClusterROC];

  // Loop through the dictionary arrays one-by-one
  // to keep memory consumption low
  AliTRDdataArrayI *tracksIn = 0;
  for (Int_t iDict = 0; iDict < kNdict; iDict++) {

    // tracksIn should be expanded beforehand!
    tracksIn = fDigitsManager->GetDictionary(idet,iDict);

    // Loop though the clusters found in this ROC
    for (iClusterROC = 0; iClusterROC < nClusterROC; iClusterROC++) {
 
      AliTRDcluster *cluster = (AliTRDcluster *)
	RecPoints()->UncheckedAt(firstClusterROC+iClusterROC);
      row  = cluster->GetLabel(0);
      col  = cluster->GetLabel(1);
      time = cluster->GetLabel(2);

      for (iPad = 0; iPad < kNclus; iPad++) {
	Int_t iPadCol = col - 1 + iPad;
	Int_t index   = tracksIn->GetDataUnchecked(row,iPadCol,time) - 1;
	idxTracks[3*iPad+iDict + iClusterROC*kNtrack] = index;     
      }

    }

  }

  // Copy the track indices into the cluster
  // Loop though the clusters found in this ROC
  for (iClusterROC = 0; iClusterROC < nClusterROC; iClusterROC++) {
 
    AliTRDcluster *cluster = (AliTRDcluster *)
      RecPoints()->UncheckedAt(firstClusterROC+iClusterROC);
    cluster->SetLabel(-9999,0);
    cluster->SetLabel(-9999,1);
    cluster->SetLabel(-9999,2);
  
    cluster->AddTrackIndex(&idxTracks[iClusterROC*kNtrack]);

  }

  delete [] idxTracks;

  return kTRUE;

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::GetCOG(Double_t signal[5])
{
  //
  // Get COG position
  // Used for clusters with more than 3 pads - where LUT not applicable
  //

  Double_t sum = signal[0]
               + signal[1]
               + signal[2] 
               + signal[3]
               + signal[4];

  Double_t res = (0.0 * (-signal[0] + signal[4])
                      + (-signal[1] + signal[3])) / sum;

  return res;		  

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::Unfold(Double_t eps, Int_t plane, Double_t *padSignal)
{
  //
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.
  //

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No AliTRDcalibDB instance available\n");
    return kFALSE;  
  }
  
  Int_t   irc                = 0;
  Int_t   itStep             = 0;                 // Count iteration steps

  Double_t ratio             = 0.5;               // Start value for ratio
  Double_t prevRatio         = 0.0;               // Store previous ratio

  Double_t newLeftSignal[3]  = { 0.0, 0.0, 0.0 }; // Array to store left cluster signal
  Double_t newRightSignal[3] = { 0.0, 0.0, 0.0 }; // Array to store right cluster signal
  Double_t newSignal[3]      = { 0.0, 0.0, 0.0 };

  // Start the iteration
  while ((TMath::Abs(prevRatio - ratio) > eps) && (itStep < 10)) {

    itStep++;
    prevRatio = ratio;

    // Cluster position according to charge ratio
    Double_t maxLeft  = (ratio*padSignal[2] - padSignal[0]) 
                      / (padSignal[0] + padSignal[1] + ratio*padSignal[2]);
    Double_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) 
                      / ((1.0 - ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // Set cluster charge ratio
    irc = calibration->PadResponse(1.0,maxLeft ,plane,newSignal);
    Double_t ampLeft  = padSignal[1] / newSignal[1];
    irc = calibration->PadResponse(1.0,maxRight,plane,newSignal);
    Double_t ampRight = padSignal[3] / newSignal[1];

    // Apply pad response to parameters
    irc = calibration->PadResponse(ampLeft ,maxLeft ,plane,newLeftSignal );
    irc = calibration->PadResponse(ampRight,maxRight,plane,newRightSignal);

    // Calculate new overlapping ratio
    ratio = TMath::Min((Double_t) 1.0
                      ,newLeftSignal[2] / (newLeftSignal[2] + newRightSignal[0]));

  }

  return ratio;

}

//_____________________________________________________________________________
void AliTRDclusterizer::TailCancelation(AliTRDdataArrayI *digitsIn
	       		              , AliTRDdataArrayF *digitsOut
				      , AliTRDSignalIndex *indexesIn
				      , AliTRDSignalIndex *indexesOut
                                      , Int_t nTimeTotal
			              , Float_t ADCthreshold
			              , AliTRDCalROC *calGainFactorROC
			              , Float_t calGainFactorDetValue)
{
  //
  // Applies the tail cancelation and gain factors: 
  // Transform digitsIn to digitsOut
  //

  Int_t iRow  = 0;
  Int_t iCol  = 0;
  Int_t iTime = 0;

  AliTRDRecParam *recParam = AliTRDRecParam::Instance();
  if (!recParam) 
    {
      AliError("No AliTRDRecParam instance available\n");
      return;
    }

  Double_t *inADC  = new Double_t[nTimeTotal];  // ADC data before tail cancellation
  Double_t *outADC = new Double_t[nTimeTotal];  // ADC data after tail cancellation
  indexesIn->ResetCounters();
  while (indexesIn->NextRCIndex(iRow, iCol))
    {
      Float_t  calGainFactorROCValue = calGainFactorROC->GetValue(iCol,iRow);
      Double_t gain                  = calGainFactorDetValue 
                                     * calGainFactorROCValue;

      for (iTime = 0; iTime < nTimeTotal; iTime++) 
	{	  
	  // Apply gain gain factor
	  inADC[iTime]   = digitsIn->GetDataUnchecked(iRow,iCol,iTime);
	  inADC[iTime]  /= gain;
	  outADC[iTime]  = inADC[iTime];
	}

      // Apply the tail cancelation via the digital filter
      if (recParam->TCOn()) 
        {
	  DeConvExp(inADC,outADC,nTimeTotal,recParam->GetTCnexp());
        }

      indexesIn->ResetTbinCounter();
      while (indexesIn->NextTbinIndex(iTime))
	{
	  // Store the amplitude of the digit if above threshold
	  if (outADC[iTime] > ADCthreshold) 
	    {
	      digitsOut->SetDataUnchecked(iRow,iCol,iTime,outADC[iTime]);
	      indexesOut->AddIndexTBin(iRow,iCol,iTime);
	    }	  
	} // while itime

    } // while irow icol
  
  delete [] inADC;
  delete [] outADC;

  return;

}

//_____________________________________________________________________________
void AliTRDclusterizer::DeConvExp(Double_t *source, Double_t *target
				 , Int_t n, Int_t nexp) 
{
  //
  // Tail cancellation by deconvolution for PASA v4 TRF
  //

  Double_t rates[2];
  Double_t coefficients[2];

  // Initialization (coefficient = alpha, rates = lambda)
  Double_t r1 = 1.0;
  Double_t r2 = 1.0;
  Double_t c1 = 0.5;
  Double_t c2 = 0.5;

  if (nexp == 1) {   // 1 Exponentials
    r1 = 1.156;
    r2 = 0.130;
    c1 = 0.066;
    c2 = 0.000;
  }
  if (nexp == 2) {   // 2 Exponentials
    r1 = 1.156;
    r2 = 0.130;
    c1 = 0.114;
    c2 = 0.624;
  }

  coefficients[0] = c1;
  coefficients[1] = c2;

  Double_t dt = 0.1;

  rates[0] = TMath::Exp(-dt/(r1));
  rates[1] = TMath::Exp(-dt/(r2));
  
  Int_t i = 0;
  Int_t k = 0;

  Double_t reminder[2];
  Double_t correction = 0.0;
  Double_t result     = 0.0;

  // Attention: computation order is important
  for (k = 0; k < nexp; k++) {
    reminder[k] = 0.0;
  }

  for (i = 0; i < n; i++) {

    result    = (source[i] - correction);    // No rescaling
    target[i] = result;

    for (k = 0; k < nexp; k++) {
      reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
    }

    correction = 0.0;
    for (k = 0; k < nexp; k++) {
      correction += reminder[k];
    }

  }

}

//_____________________________________________________________________________
void AliTRDclusterizer::ResetRecPoints() 
{
  //
  // Resets the list of rec points
  //

  if (fRecPoints) {
    fRecPoints->Delete();
  }

}

//_____________________________________________________________________________
TObjArray *AliTRDclusterizer::RecPoints() 
{
  //
  // Returns the list of rec points
  //

  if (!fRecPoints) {
    fRecPoints = new TObjArray(400);
  }
 
  return fRecPoints;

}
