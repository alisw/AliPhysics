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
#include "AliPMDisocell.h"
#include "AliPMDRawStream.h"
#include "AliPMDCalibData.h"
#include "AliPMDPedestal.h"
#include "AliPMDddldata.h"

#include "AliDAQ.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"



ClassImp(AliPMDClusterFinder)

AliPMDClusterFinder::AliPMDClusterFinder():
  fRunLoader(0),
  fPMDLoader(0),
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
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
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
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
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
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
	fDigits->Clear();
    }
  if (fRecpoints)
    {
      fRecpoints->Clear();
    }
  if (fRechits)
    {
      fRechits->Clear();
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

  AliPMDisocell *pmdiso = 0x0;

  TObjArray *pmdcont = new TObjArray();
  TObjArray *pmdisocell = new TObjArray();

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
	  if(xpos < 0 || xpos > 48 || ypos < 0 || ypos > 96)
	    {
	      AliError(Form("*Row %d and Column NUMBER %d NOT Valid *",
			      xpos, ypos));
	      continue; 
	    }
	  // CALIBRATION
	  Float_t gain = fCalibGain->GetGainFact(det,smn,xpos,ypos);
	  // printf("adc = %d gain = %f\n",adc,gain);
	  
	  adc = adc*gain;

	  //Int_t trno   = pmddigit->GetTrackNumber();
	  fCellTrack[xpos][ypos] = pmddigit->GetTrackNumber();
	  fCellPid[xpos][ypos]   = pmddigit->GetTrackPid();
	  fCellADC[xpos][ypos]   = (Double_t) adc;
	}

      idet = det;
      ismn = smn;
      pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,
			pmdisocell,pmdcont);

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
	      Int_t celldataTr = pmdcl->GetClusCellTrack(ihit);
	      Int_t celldataPid   = pmdcl->GetClusCellPid(ihit);
	      Float_t celldataAdc = pmdcl->GetClusCellAdc(ihit);
	      AddRecHit(celldataX, celldataY, celldataTr, celldataPid, celldataAdc);
	    }
	  branch2->Fill();
	  ResetRechit();
	}
      pmdcont->Delete();

	  // Added single isolated cell for offline gain calibration
	  nentries1 = pmdisocell->GetEntries();
	  AliDebug(1,Form("Total number of isolated single cell clusters = %d",nentries1));
	  
	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      pmdiso = (AliPMDisocell*)pmdisocell->UncheckedAt(ient1);
	      idet = pmdiso->GetDetector();
	      ismn = pmdiso->GetSmn();
	      clusdata[0] = (Float_t) pmdiso->GetRow();
	      clusdata[1] = (Float_t) pmdiso->GetCol();
	      clusdata[2] = pmdiso->GetADC();
	      clusdata[3] = 1.;
	      clusdata[4] = -99.;
	      clusdata[5] = -99.;
      
	      AddRecPoint(idet,ismn,clusdata);
	    }
	  pmdisocell->Delete();
      
      branch1->Fill();
      ResetRecpoint();

    } // modules

  ResetCellADC();
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");  
  fPMDLoader->WriteRecPoints("OVERWRITE");

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
  delete pmdisocell;
    
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(TTree *digitsTree,
					   TTree *clustersTree)
{
  // Converts digits to recpoints after running clustering
  // algorithm on CPV plane and PREshower plane
  //
  // This algorithm is called during the reconstruction from digits

  Int_t    det = 0,smn = 0;
  Int_t    xpos,ypos;
  Float_t  adc;
  Int_t    ismn;
  Int_t    idet;
  Float_t  clusdata[6];

  AliPMDcluster *pmdcl = 0x0;
  AliPMDisocell *pmdiso = 0x0;

  TObjArray *pmdcont = new TObjArray();
  TObjArray *pmdisocell = new TObjArray();
  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  pmdclust->SetEdepCut(fEcut);

  AliPMDdigit  *pmddigit;
  TBranch *branch = digitsTree->GetBranch("PMDDigit");
  branch->SetAddress(&fDigits);

  ResetRecpoint();

  Int_t bufsize = 16000;
  TBranch * branch1 = clustersTree->Branch("PMDRecpoint", &fRecpoints, bufsize); 
  TBranch * branch2 = clustersTree->Branch("PMDRechit", &fRechits, bufsize); 

  Int_t nmodules = (Int_t) digitsTree->GetEntries();

  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {

      Int_t totADCMod = 0;
      ResetCellADC();
      digitsTree->GetEntry(imodule); 
      Int_t nentries = fDigits->GetLast();
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
	  
	  det    = pmddigit->GetDetector();
	  smn    = pmddigit->GetSMNumber();
	  xpos   = pmddigit->GetRow();
	  ypos   = pmddigit->GetColumn();
	  adc    = pmddigit->GetADC();
	  if(xpos < 0 || xpos > 48 || ypos < 0 || ypos > 96)
	    {
	      AliError(Form("*Row %d and Column NUMBER %d NOT Valid *",
			    xpos, ypos));
	      continue; 
	    }
	  
	  // Pedestal Subtraction
	  Int_t   pedmeanrms = fCalibPed->GetPedMeanRms(det,smn,xpos,ypos);
	  Int_t   pedrms1    = (Int_t) pedmeanrms%100;
	  Float_t pedrms     = (Float_t)pedrms1/10.;
	  Float_t pedmean    = (Float_t) (pedmeanrms - pedrms1)/1000.0;
	  //printf("%f %f\n",pedmean, pedrms);

	  Float_t adc1 = adc - (pedmean + 3.0*pedrms);

	  // CALIBRATION
	  Float_t gain = fCalibGain->GetGainFact(det,smn,xpos,ypos);
	  // printf("adc = %d gain = %f\n",adc,gain);

	  adc = adc1*gain;

	  fCellTrack[xpos][ypos] = pmddigit->GetTrackNumber();
	  fCellPid[xpos][ypos] = pmddigit->GetTrackPid();
	  fCellADC[xpos][ypos] = (Double_t) adc;

	  totADCMod += (Int_t) adc;

	}

      idet = det;
      ismn = smn;

      if (totADCMod <= 0) continue;

      pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,
			pmdisocell,pmdcont);
      
      Int_t nentries1 = pmdcont->GetEntries();

      AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

      for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	{
	  pmdcl = (AliPMDcluster*)pmdcont->UncheckedAt(ient1);
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
	      Int_t celldataTr = pmdcl->GetClusCellTrack(ihit);
	      Int_t celldataPid   = pmdcl->GetClusCellPid(ihit);
	      Float_t celldataAdc = pmdcl->GetClusCellAdc(ihit);
	      AddRecHit(celldataX, celldataY, celldataTr, celldataPid, celldataAdc);
	    }
	  branch2->Fill();
	  ResetRechit();
	}
      pmdcont->Delete();

      // Added single isolated cell for offline gain calibration
      nentries1 = pmdisocell->GetEntries();
      AliDebug(1,Form("Total number of isolated single cell clusters = %d",nentries1));
      
      for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	{
	  pmdiso = (AliPMDisocell*)pmdisocell->UncheckedAt(ient1);
	  idet = pmdiso->GetDetector();
	  ismn = pmdiso->GetSmn();
	  clusdata[0] = (Float_t) pmdiso->GetRow();
	  clusdata[1] = (Float_t) pmdiso->GetCol();
	  clusdata[2] = pmdiso->GetADC();
	  clusdata[3] = 1.;
	  clusdata[4] = -99.;
	  clusdata[5] = -99.;
	  
	  AddRecPoint(idet,ismn,clusdata);
	}
      pmdisocell->Delete();
      
      branch1->Fill();
      ResetRecpoint();

    } // modules


  ResetCellADC();

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
  delete pmdisocell;
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(AliRawReader *rawReader,
					   TTree *clustersTree)
{
  // Converts RAW data to recpoints after running clustering
  // algorithm on CPV and PREshower plane
  //
  // This method is called at the time of reconstruction from RAW data


  AliPMDddldata *pmdddl = 0x0;
  AliPMDcluster *pmdcl  = 0x0;
  AliPMDisocell *pmdiso = 0x0;


  Float_t  clusdata[6];
  TObjArray pmdddlcont;

  TObjArray *pmdcont = new TObjArray();
  TObjArray *pmdisocell = new TObjArray();

  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  pmdclust->SetEdepCut(fEcut);

  ResetRecpoint();

  Int_t bufsize = 16000;
  TBranch *branch1 = clustersTree->Branch("PMDRecpoint", &fRecpoints, bufsize); 

  TBranch * branch2 = clustersTree->Branch("PMDRechit", &fRechits, bufsize); 

  const Int_t kRow = 48;
  const Int_t kCol = 96;

  Int_t idet = 0;
  Int_t iSMN = 0;

  Int_t indexDDL = -1;
  AliPMDRawStream pmdinput(rawReader);

  while ((indexDDL = pmdinput.DdlData(&pmdddlcont)) >=0)
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

      Int_t indexsmn = 0;
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();

	  if(smn == -1)
	    {
	      AliError(Form("*MODULE NUMBER WRONG %d *",smn));
	      continue; 
	    }
	  if(row < 0 || row > 48 || col < 0 || col > 96)
	    {
	      AliError(Form("*Row %d and Column NUMBER %d NOT Valid *",
			    row, col));
	      continue; 
	    }

	  // Pedestal Subtraction
	  Int_t   pedmeanrms = fCalibPed->GetPedMeanRms(det,smn,row,col);
	  Int_t   pedrms1    = (Int_t) pedmeanrms%100;
	  Float_t pedrms     = (Float_t)pedrms1/10.;
	  Float_t pedmean    = (Float_t) (pedmeanrms - pedrms1)/1000.0;

	  //printf("%f %f\n",pedmean, pedrms);

	  // Float_t sig1 = (Float_t) sig;
	  Float_t sig1 = (Float_t) sig - (pedmean + 3.0*pedrms);

	  // CALIBRATION
	  Float_t gain = fCalibGain->GetGainFact(det,smn,row,col);
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
      
      pmdddlcont.Delete();

      Int_t totAdcMod = 0;

      Int_t ismn = 0;
      for (indexsmn = 0; indexsmn < iSMN; indexsmn++)
	{
	  ResetCellADC();
	  totAdcMod = 0;
	  for (Int_t irow = 0; irow < kRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kCol; icol++)
		{
		  fCellTrack[irow][icol] = -1;
		  fCellPid[irow][icol]   = -1;

		  fCellADC[irow][icol] = 
		    (Double_t) precpvADC[indexsmn][irow][icol];
		  totAdcMod += precpvADC[indexsmn][irow][icol];
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

	  if (totAdcMod <= 0) continue;


	  pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,
			    pmdisocell,pmdcont);

	  Int_t nentries1 = pmdcont->GetEntries();

	  AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      pmdcl = (AliPMDcluster*)pmdcont->UncheckedAt(ient1);
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
		  Int_t celldataTr = pmdcl->GetClusCellTrack(ihit);
		  Int_t celldataPid   = pmdcl->GetClusCellPid(ihit);
		  Float_t celldataAdc = pmdcl->GetClusCellAdc(ihit);
		  AddRecHit(celldataX, celldataY, celldataTr, celldataPid, celldataAdc);
		}
	      branch2->Fill();
	      ResetRechit();

	    }
	  pmdcont->Delete();


	  // Added single isolated cell for offline gain calibration
	  nentries1 = pmdisocell->GetEntries();
	  AliDebug(1,Form("Total number of isolated single cell clusters = %d",nentries1));
	  
	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      pmdiso = (AliPMDisocell*)pmdisocell->UncheckedAt(ient1);
	      idet = pmdiso->GetDetector();
	      ismn = pmdiso->GetSmn();
	      clusdata[0] = (Float_t) pmdiso->GetRow();
	      clusdata[1] = (Float_t) pmdiso->GetCol();
	      clusdata[2] = pmdiso->GetADC();
	      clusdata[3] = 1.;
	      clusdata[4] = -99.;
	      clusdata[5] = -99.;
      
	      AddRecPoint(idet,ismn,clusdata);
	    }
	  pmdisocell->Delete();
	  
	  branch1->Fill();
	  ResetRecpoint();


	} // smn

      for (Int_t i=0; i<iSMN; i++)
	{
	  for (Int_t j=0; j<kRow; j++) delete [] precpvADC[i][j];
	}
      for (Int_t i=0; i<iSMN; i++) delete [] precpvADC[i];
      delete [] precpvADC;

    } // DDL Loop

  
  ResetCellADC();
  
  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
  delete pmdisocell;

}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(Int_t ievt, AliRawReader *rawReader)
{
  // Converts RAW data to recpoints after running clustering
  // algorithm on CPV and PREshower plane
  //

  Float_t  clusdata[6];

  TObjArray pmdddlcont;

  AliPMDcluster *pmdcl  = 0x0;
  AliPMDisocell *pmdiso  = 0x0;


  TObjArray *pmdcont = new TObjArray();
  TObjArray *pmdisocell = new TObjArray();

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

  const Int_t kRow = 48;
  const Int_t kCol = 96;

  Int_t idet = 0;
  Int_t iSMN = 0;

  AliPMDRawStream pmdinput(rawReader);
  Int_t indexDDL = -1;

  while ((indexDDL = pmdinput.DdlData(&pmdddlcont)) >=0) {

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
	  if(row < 0 || row > 48 || col < 0 || col > 96)
	    {
	      AliError(Form("*Row %d and Column NUMBER %d NOT Valid *",
			      row, col));
	      continue; 
	    }
	  // Pedestal Subtraction
	  Int_t   pedmeanrms = fCalibPed->GetPedMeanRms(det,smn,row,col);
	  Int_t   pedrms1    = (Int_t) pedmeanrms%100;
	  Float_t pedrms     = (Float_t)pedrms1/10.;
	  Float_t pedmean    = (Float_t) (pedmeanrms - pedrms1)/1000.0;

	  //printf("%f %f\n",pedmean, pedrms);

	  //Float_t sig1 = (Float_t) sig;
	  Float_t sig1 = (Float_t) sig - (pedmean + 3.0*pedrms);
	  // CALIBRATION
	  Float_t gain = fCalibGain->GetGainFact(det,smn,row,col);

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
      
      pmdddlcont.Delete();

      Int_t ismn = 0;
      for (indexsmn = 0; indexsmn < iSMN; indexsmn++)
	{
	  ResetCellADC();
	  for (Int_t irow = 0; irow < kRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kCol; icol++)
		{
		  fCellTrack[irow][icol] = -1;
		  fCellPid[irow][icol]   = -1;
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

	  pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,
			    pmdisocell,pmdcont);

	  Int_t nentries1 = pmdcont->GetEntries();

	  AliDebug(1,Form("Total number of clusters/module = %d",nentries1));

	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      pmdcl       = (AliPMDcluster*)pmdcont->UncheckedAt(ient1);
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
		  Int_t celldataTr = pmdcl->GetClusCellTrack(ihit);
		  Int_t celldataPid = pmdcl->GetClusCellPid(ihit);
		  Float_t celldataAdc = pmdcl->GetClusCellAdc(ihit);
		  AddRecHit(celldataX, celldataY, celldataTr, celldataPid, celldataAdc);
		}
	      branch2->Fill();
	      ResetRechit();

	    }
	  pmdcont->Delete();

	  // Added single isolated cell for offline gain calibration
	  nentries1 = pmdisocell->GetEntries();
	  AliDebug(1,Form("Total number of isolated single cell clusters = %d",nentries1));
	  
	  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	    {
	      pmdiso = (AliPMDisocell*)pmdisocell->UncheckedAt(ient1);
	      idet = pmdiso->GetDetector();
	      ismn = pmdiso->GetSmn();
	      clusdata[0] = (Float_t) pmdiso->GetRow();
	      clusdata[1] = (Float_t) pmdiso->GetCol();
	      clusdata[2] = pmdiso->GetADC();
	      clusdata[3] = 1.;
	      clusdata[4] = -99.;
	      clusdata[5] = -99.;
      
	      AddRecPoint(idet,ismn,clusdata);
	    }
	  pmdisocell->Delete();
	  
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
  delete pmdisocell;
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
void AliPMDClusterFinder::AddRecHit(Int_t celldataX,Int_t celldataY,
				    Int_t celldataTr, Int_t celldataPid,
				    Float_t celldataAdc)
{
  // Add associated cell hits to the Reconstructed points
  //
  TClonesArray &lrechits = *fRechits;
  AliPMDrechit *newrechit;
  newrechit = new AliPMDrechit(celldataX, celldataY, celldataTr, celldataPid, celldataAdc);
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
	  fCellTrack[irow][icol] = -1;
	  fCellPid[irow][icol]   = -1;
	  fCellADC[irow][icol]   = 0.;
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
AliPMDCalibData* AliPMDClusterFinder::GetCalibGain() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  // Added by ZA
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Gain");
  
  if(!entry)  AliFatal("Calibration object retrieval failed! ");
  
  AliPMDCalibData *calibdata=0;
  if (entry) calibdata = (AliPMDCalibData*) entry->GetObject();
  
  if (!calibdata)  AliFatal("No calibration data from calibration database !");
  
  return calibdata;
}
// ------------------------------------------------------------------------- //
AliPMDPedestal* AliPMDClusterFinder::GetCalibPed() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Ped");
  
  if(!entry) AliFatal("Pedestal object retrieval failed!");
    
  AliPMDPedestal *pedestal = 0;
  if (entry) pedestal = (AliPMDPedestal*) entry->GetObject();
  
  if (!pedestal)  AliFatal("No pedestal data from pedestal database !");
  
  return pedestal;
}
