/***************************************************************************
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

//-----------------------------------------------------//
//                                                     //
//           Date   : August 05 2003                   //
//  This reads the file PMD.digits.root(TreeD),        //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h"

#include "AliPMDdigit.h"
#include "AliPMDClusterFinder.h"
#include "AliPMDClustering.h"
#include "AliPMDClusteringV1.h"
#include "AliPMDcluster.h"
#include "AliPMDrecpoint1.h"
#include "AliPMDrechit.h"
#include "AliPMDRawStream.h"
#include "AliPMDCalibData.h"
#include "AliPMDddldata.h"

#include "AliDAQ.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"



ClassImp(AliPMDClusterFinder)

AliPMDClusterFinder::AliPMDClusterFinder():
  fRunLoader(0),
  fPMDLoader(0),
  fCalibData(GetCalibData()),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fRecpoints(new TClonesArray("AliPMDrecpoint1", 1000)),
  fRechits(new TClonesArray("AliPMDrechit", 1000)),
  fNpoint(0),
  fNhit(0),
  fDetNo(0),
  fEcut(0.)
{
//
// Constructor
//
}
// ------------------------------------------------------------------------- //
AliPMDClusterFinder::AliPMDClusterFinder(AliRunLoader* runLoader):
  fRunLoader(runLoader),
  fPMDLoader(runLoader->GetLoader("PMDLoader")),
  fCalibData(GetCalibData()),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fRecpoints(new TClonesArray("AliPMDrecpoint1", 1000)),
  fRechits(new TClonesArray("AliPMDrechit", 1000)),
  fNpoint(0),
  fNhit(0),
  fDetNo(0),
  fEcut(0.)
{
//
// Constructor
//
}
// ------------------------------------------------------------------------- //
AliPMDClusterFinder::AliPMDClusterFinder(const AliPMDClusterFinder & finder):
  TObject(finder),
  fRunLoader(0),
  fPMDLoader(0),
  fCalibData(GetCalibData()),
  fTreeD(0),
  fTreeR(0),
  fDigits(NULL),
  fRecpoints(NULL),
  fRechits(NULL),
  fNpoint(0),
  fNhit(0),
  fDetNo(0),
  fEcut(0.)
{
  // copy constructor
  AliError("Copy constructor not allowed");
}
// ------------------------------------------------------------------------- //
AliPMDClusterFinder &AliPMDClusterFinder::operator=(const AliPMDClusterFinder & /*finder*/)
{
 // assignment op
  AliError("Assignment Operator not allowed");
  return *this;
}
// ------------------------------------------------------------------------- //
AliPMDClusterFinder::~AliPMDClusterFinder()
{
  // Destructor
  if (fDigits)
    {
      fDigits->Delete();
      delete fDigits;
      fDigits=0;
    }
  if (fRecpoints)
    {
      fRecpoints->Delete();
      delete fRecpoints;
      fRecpoints=0;
    }
  if (fRechits)
    {
      fRechits->Delete();
      delete fRechits;
      fRechits=0;
    }
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(Int_t ievt)
{
  // Converts digits to recpoints after running clustering
  // algorithm on CPV plane and PREshower plane
  //

  Int_t    det = 0,smn = 0;
  Int_t    xpos,ypos;
  Float_t  adc;
  Int_t    ismn;
  Int_t    idet;
  Float_t  clusdata[6];

  TObjArray *pmdcont = new TObjArray();
  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  pmdclust->SetEdepCut(fEcut);

  fRunLoader->GetEvent(ievt);


  fTreeD = fPMDLoader->TreeD();
  if (fTreeD == 0x0)
    {
      AliFatal("AliPMDClusterFinder: Can not get TreeD");

    }
  AliPMDdigit  *pmddigit;
  TBranch *branch = fTreeD->GetBranch("PMDDigit");
  branch->SetAddress(&fDigits);

  ResetRecpoint();

  fTreeR = fPMDLoader->TreeR();
  if (fTreeR == 0x0)
    {
      fPMDLoader->MakeTree("R");
      fTreeR = fPMDLoader->TreeR();
    }

  Int_t bufsize = 16000;
  TBranch * branch1 = fTreeR->Branch("PMDRecpoint", &fRecpoints, bufsize); 
  TBranch * branch2 = fTreeR->Branch("PMDRechit", &fRechits, bufsize); 

  Int_t nmodules = (Int_t) fTreeD->GetEntries();

  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      ResetCellADC();
      fTreeD->GetEntry(imodule); 
      Int_t nentries = fDigits->GetLast();
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
	  
	  det    = pmddigit->GetDetector();
	  smn    = pmddigit->GetSMNumber();
	  xpos   = pmddigit->GetRow();
	  ypos   = pmddigit->GetColumn();
	  adc    = pmddigit->GetADC();
	  
	  // CALIBRATION
	  Float_t gain = fCalibData->GetGainFact(det,smn,xpos,ypos);
	  // printf("adc = %d gain = %f\n",adc,gain);
	  
	  adc = adc*gain;

	  //Int_t trno   = pmddigit->GetTrackNumber();
	  fCellADC[xpos][ypos] = (Double_t) adc;
	}

      idet = det;
      ismn = smn;
      pmdclust->DoClust(idet,ismn,fCellADC,pmdcont);
      
      Int_t nentries1 = pmdcont->GetEntries();

      AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

      for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	{
	  AliPMDcluster *pmdcl = (AliPMDcluster*)pmdcont->UncheckedAt(ient1);
	  idet        = pmdcl->GetDetector();
	  ismn        = pmdcl->GetSMN();
	  clusdata[0] = pmdcl->GetClusX();
	  clusdata[1] = pmdcl->GetClusY();
	  clusdata[2] = pmdcl->GetClusADC();
	  clusdata[3] = pmdcl->GetClusCells();
	  clusdata[4] = pmdcl->GetClusSigmaX();
	  clusdata[5] = pmdcl->GetClusSigmaY();

	  AddRecPoint(idet,ismn,clusdata);

	  Int_t ncell = (Int_t) clusdata[3];
	  for(Int_t ihit = 0; ihit < ncell; ihit++)
	    {
	      Int_t celldataX = pmdcl->GetClusCellX(ihit);
	      Int_t celldataY = pmdcl->GetClusCellY(ihit);
	      AddRecHit(celldataX, celldataY);
	    }
	  branch2->Fill();
	  ResetRechit();
	}
      pmdcont->Clear();
      
      branch1->Fill();
      ResetRecpoint();

    } // modules

  ResetCellADC();
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");  
  fPMDLoader->WriteRecPoints("OVERWRITE");

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
    
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(AliRawReader *rawReader,
					   TTree *clustersTree)
{
  // Converts RAW data to recpoints after running clustering
  // algorithm on CPV and PREshower plane
  //
  // This method is called at the time of reconstruction


  Float_t  clusdata[6];
  TObjArray pmdddlcont;

  TObjArray *pmdcont = new TObjArray();
  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  pmdclust->SetEdepCut(fEcut);

  ResetRecpoint();

  Int_t bufsize = 16000;
  TBranch *branch1 = clustersTree->Branch("PMDRecpoint", &fRecpoints, bufsize); 

  TBranch * branch2 = clustersTree->Branch("PMDRechit", &fRechits, bufsize); 

  const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");
  const Int_t kRow = 48;
  const Int_t kCol = 96;

  Int_t idet = 0;
  Int_t iSMN = 0;

  
  for (Int_t indexDDL = 0; indexDDL < kDDL; indexDDL++)
    {
      if (indexDDL < 4)
	{
	  iSMN = 6;
	}
      else if (indexDDL >= 4)
	{
	  iSMN = 12;
	}
      Int_t ***precpvADC;
      precpvADC = new int **[iSMN];
      for (Int_t i=0; i<iSMN; i++) precpvADC[i] = new int *[kRow];
      for (Int_t i=0; i<iSMN;i++)
	{
	  for (Int_t j=0; j<kRow; j++) precpvADC[i][j] = new int [kCol];
	}
      for (Int_t i = 0; i < iSMN; i++)
	{
	  for (Int_t j = 0; j < kRow; j++)
	    {
	      for (Int_t k = 0; k < kCol; k++)
		{
		  precpvADC[i][j][k] = 0;
		}
	    }
	}
      ResetCellADC();
      rawReader->Reset();
      AliPMDRawStream pmdinput(rawReader);

      rawReader->Select("PMD", indexDDL, indexDDL);

      pmdinput.DdlData(indexDDL,&pmdddlcont);

      Int_t indexsmn = 0;
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();

	  Float_t sig1 = (Float_t) sig;
	  // CALIBRATION
	  Float_t gain = fCalibData->GetGainFact(det,smn,row,col);
	  //printf("sig = %d gain = %f\n",sig,gain);
	  sig = (Int_t) (sig1*gain);

	  if (indexDDL < 4)
	    {
	      if (det != 0)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn - indexDDL * 6;
	    }
	  else if (indexDDL == 4)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      if (smn < 6)
		{
		  indexsmn = smn;
		}
	      else if (smn >= 18 && smn < 24)
		{
		  indexsmn = smn - 12;
		}
	    }
	  else if (indexDDL == 5)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      if (smn >= 6 && smn < 12)
		{
		  indexsmn = smn - 6;
		}
	      else if (smn >= 12 && smn < 18)
		{
		  indexsmn = smn - 6;
		}
	    }	      
	  precpvADC[indexsmn][row][col] = sig;
	}
      
      pmdddlcont.Clear();

      Int_t ismn = 0;
      for (Int_t indexsmn = 0; indexsmn < iSMN; indexsmn++)
	{
	  ResetCellADC();
	  for (Int_t irow = 0; irow < kRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kCol; icol++)
		{
		  fCellADC[irow][icol] = 
		    (Double_t) precpvADC[indexsmn][irow][icol];
		} // row
	    }     // col
	  
	  if (indexDDL < 4)
	    {
	      ismn = indexsmn + indexDDL * 6;
	      idet = 0;
	    }
	  else if (indexDDL == 4)
	    {
	      if (indexsmn < 6)
		{
		  ismn = indexsmn;
		}
	      else if (indexsmn >= 6 && indexsmn < 12)
		{
		  ismn = indexsmn + 12;
		}
	      idet = 1;
	    }
	  else if (indexDDL == 5)
	    {
	      if (indexsmn < 6)
		{
		  ismn = indexsmn + 6;
		}
	      else if (indexsmn >= 6 && indexsmn < 12)
		{
		  ismn = indexsmn + 6;
		}
	      idet = 1;
	    }

	  pmdclust->DoClust(idet,ismn,fCellADC,pmdcont);
	  Int_t nentries1 = pmdcont->GetEntries();

	  AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      AliPMDcluster *pmdcl = 
		(AliPMDcluster*)pmdcont->UncheckedAt(ient1);
	      idet        = pmdcl->GetDetector();
	      ismn        = pmdcl->GetSMN();
	      clusdata[0] = pmdcl->GetClusX();
	      clusdata[1] = pmdcl->GetClusY();
	      clusdata[2] = pmdcl->GetClusADC();
	      clusdata[3] = pmdcl->GetClusCells();
	      clusdata[4] = pmdcl->GetClusSigmaX();
	      clusdata[5] = pmdcl->GetClusSigmaY();

	      AddRecPoint(idet,ismn,clusdata);

	      Int_t ncell = (Int_t) clusdata[3];
	      for(Int_t ihit = 0; ihit < ncell; ihit++)
		{
		  Int_t celldataX = pmdcl->GetClusCellX(ihit);
		  Int_t celldataY = pmdcl->GetClusCellY(ihit);
		  AddRecHit(celldataX, celldataY);
		}
	      branch2->Fill();
	      ResetRechit();

	    }
	  pmdcont->Clear();
	  
	  branch1->Fill();
	  ResetRecpoint();


	} // smn

      for (Int_t i=0; i<iSMN; i++)
	{
	  for (Int_t j=0; j<kRow; j++) delete [] precpvADC[i][j];
	}
      for (Int_t i=0; i<iSMN; i++) delete [] precpvADC[i];
      delete precpvADC;
    } // DDL Loop
  
  ResetCellADC();
  
  //   delete the pointers
  delete pmdclust;
  delete pmdcont;

}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(Int_t ievt, AliRawReader *rawReader)
{
  // Converts RAW data to recpoints after running clustering
  // algorithm on CPV and PREshower plane
  //

  Float_t  clusdata[6];
  TObjArray pmdddlcont;
  TObjArray *pmdcont = new TObjArray();

  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  pmdclust->SetEdepCut(fEcut);

  fRunLoader->GetEvent(ievt);

  ResetRecpoint();

  fTreeR = fPMDLoader->TreeR();
  if (fTreeR == 0x0)
    {
      fPMDLoader->MakeTree("R");
      fTreeR = fPMDLoader->TreeR();
    }
  Int_t bufsize = 16000;
  TBranch *branch1 = fTreeR->Branch("PMDRecpoint", &fRecpoints, bufsize); 
  TBranch *branch2 = fTreeR->Branch("PMDRechit", &fRechits, bufsize); 

  const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");
  const Int_t kRow = 48;
  const Int_t kCol = 96;

  Int_t idet = 0;
  Int_t iSMN = 0;
  
  for (Int_t indexDDL = 0; indexDDL < kDDL; indexDDL++)
    {
      if (indexDDL < 4)
	{
	  iSMN = 6;
	}
      else if (indexDDL >= 4)
	{
	  iSMN = 12;
	}
      Int_t ***precpvADC;
      precpvADC = new int **[iSMN];
      for (Int_t i=0; i<iSMN; i++) precpvADC[i] = new int *[kRow];
      for (Int_t i=0; i<iSMN;i++)
	{
	  for (Int_t j=0; j<kRow; j++) precpvADC[i][j] = new int [kCol];
	}
      for (Int_t i = 0; i < iSMN; i++)
	{
	  for (Int_t j = 0; j < kRow; j++)
	    {
	      for (Int_t k = 0; k < kCol; k++)
		{
		  precpvADC[i][j][k] = 0;
		}
	    }
	}
      ResetCellADC();
      rawReader->Reset();
      rawReader->Select("PMD", indexDDL, indexDDL);

      AliPMDRawStream pmdinput(rawReader);
      pmdinput.DdlData(indexDDL,&pmdddlcont);
    
      Int_t indexsmn = 0;
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();

	  Float_t sig1 = (Float_t) sig;
	  // CALIBRATION
	  Float_t gain = fCalibData->GetGainFact(det,smn,row,col);

	  //printf("sig = %d gain = %f\n",sig,gain);
	  sig = (Int_t) (sig1*gain);


	  if (indexDDL < 4)
	    {
	      if (det != 0)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn - indexDDL * 6;
	    }
	  else if (indexDDL == 4)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      if (smn < 6)
		{
		  indexsmn = smn;
		}
	      else if (smn >= 18 && smn < 24)
		{
		  indexsmn = smn - 12;
		}
	    }
	  else if (indexDDL == 5)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      if (smn >= 6 && smn < 12)
		{
		  indexsmn = smn - 6;
		}
	      else if (smn >= 12 && smn < 18)
		{
		  indexsmn = smn - 6;
		}
	    }	      
	  precpvADC[indexsmn][row][col] = sig;

	}
      
      pmdddlcont.Clear();

      Int_t ismn = 0;
      for (Int_t indexsmn = 0; indexsmn < iSMN; indexsmn++)
	{
	  ResetCellADC();
	  for (Int_t irow = 0; irow < kRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kCol; icol++)
		{
		  fCellADC[irow][icol] = 
		    (Double_t) precpvADC[indexsmn][irow][icol];
		} // row
	    }     // col

	  
	  if (indexDDL < 4)
	    {
	      ismn = indexsmn + indexDDL * 6;
	      idet = 0;
	    }
	  else if (indexDDL == 4)
	    {
	      if (indexsmn < 6)
		{
		  ismn = indexsmn;
		}
	      else if (indexsmn >= 6 && indexsmn < 12)
		{
		  ismn = indexsmn + 12;
		}
	      idet = 1;
	    }
	  else if (indexDDL == 5)
	    {
	      if (indexsmn < 6)
		{
		  ismn = indexsmn + 6;
		}
	      else if (indexsmn >= 6 && indexsmn < 12)
		{
		  ismn = indexsmn + 6;
		}
	      idet = 1;
	    }

	  pmdclust->DoClust(idet,ismn,fCellADC,pmdcont);
	  Int_t nentries1 = pmdcont->GetEntries();

	  AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      AliPMDcluster *pmdcl = 
		(AliPMDcluster*)pmdcont->UncheckedAt(ient1);
	      idet        = pmdcl->GetDetector();
	      ismn        = pmdcl->GetSMN();
	      clusdata[0] = pmdcl->GetClusX();
	      clusdata[1] = pmdcl->GetClusY();
	      clusdata[2] = pmdcl->GetClusADC();
	      clusdata[3] = pmdcl->GetClusCells();
	      clusdata[4] = pmdcl->GetClusSigmaX();
	      clusdata[5] = pmdcl->GetClusSigmaY();

	      AddRecPoint(idet,ismn,clusdata);

	      Int_t ncell = (Int_t) clusdata[3];
	      for(Int_t ihit = 0; ihit < ncell; ihit++)
		{
		  Int_t celldataX = pmdcl->GetClusCellX(ihit);
		  Int_t celldataY = pmdcl->GetClusCellY(ihit);
		  AddRecHit(celldataX, celldataY);
		}
	      branch2->Fill();
	      ResetRechit();

	    }
	  pmdcont->Clear();
	  
	  branch1->Fill();
	  ResetRecpoint();


	} // smn

      for (Int_t i=0; i<iSMN; i++)
	{
	  for (Int_t j=0; j<kRow; j++) delete [] precpvADC[i][j];
	}
      for (Int_t i=0; i<iSMN; i++) delete [] precpvADC[i];
      delete precpvADC;
    } // DDL Loop
  
  ResetCellADC();
  
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");  
  fPMDLoader->WriteRecPoints("OVERWRITE");

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;

}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::SetCellEdepCut(Float_t ecut)
{
  fEcut = ecut;
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::AddRecPoint(Int_t idet,Int_t ismn,Float_t *clusdata)
{
  // Add Reconstructed points
  //
  TClonesArray &lrecpoints = *fRecpoints;
  AliPMDrecpoint1 *newrecpoint;
  newrecpoint = new AliPMDrecpoint1(idet, ismn, clusdata);
  new(lrecpoints[fNpoint++]) AliPMDrecpoint1(newrecpoint);
  delete newrecpoint;
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::AddRecHit(Int_t celldataX,Int_t celldataY)
{
  // Add associated cell hits to the Reconstructed points
  //
  TClonesArray &lrechits = *fRechits;
  AliPMDrechit *newrechit;
  newrechit = new AliPMDrechit(celldataX, celldataY);
  new(lrechits[fNhit++]) AliPMDrechit(newrechit);
  delete newrechit;
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::ResetCellADC()
{
  // Reset the individual cell ADC value to zero
  //
  for(Int_t irow = 0; irow < fgkRow; irow++)
    {
      for(Int_t icol = 0; icol < fgkCol; icol++)
	{
	  fCellADC[irow][icol] = 0.;
	}
    }
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::ResetRecpoint()
{
  // Clear the list of reconstructed points
  fNpoint = 0;
  if (fRecpoints) fRecpoints->Clear();
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::ResetRechit()
{
  // Clear the list of reconstructed points
  fNhit = 0;
  if (fRechits) fRechits->Clear();
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::Load()
{
  // Load all the *.root files
  //
  fPMDLoader->LoadDigits("READ");
  fPMDLoader->LoadRecPoints("recreate");
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::LoadClusters()
{
  // Load all the *.root files
  //
  fPMDLoader->LoadRecPoints("recreate");
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::UnLoad()
{
  // Unload all the *.root files
  //
  fPMDLoader->UnloadDigits();
  fPMDLoader->UnloadRecPoints();
}
// ------------------------------------------------------------------------- //
void AliPMDClusterFinder::UnLoadClusters()
{
  // Unload all the *.root files
  //
  fPMDLoader->UnloadRecPoints();
}
// ------------------------------------------------------------------------- //

AliPMDCalibData* AliPMDClusterFinder::GetCalibData() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  // Added by ZA
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Data");
  
  if(!entry){
    AliWarning("Calibration object retrieval failed! Dummy calibration will be used.");
    
    // this just remembers the actual default storage. No problem if it is null.
    AliCDBStorage *origStorage = AliCDBManager::Instance()->GetDefaultStorage();
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    
    entry = AliCDBManager::Instance()->Get("PMD/Calib/Data");
    
    // now reset the original default storage to AliCDBManager...
    AliCDBManager::Instance()->SetDefaultStorage(origStorage);  
  }
  
  AliPMDCalibData *calibdata=0;
  if (entry) calibdata = (AliPMDCalibData*) entry->GetObject();
  
  if (!calibdata)  AliError("No calibration data from calibration database !");
  
  return calibdata;
}
