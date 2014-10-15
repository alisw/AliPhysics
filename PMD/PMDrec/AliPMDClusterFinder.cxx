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
#include <TSystem.h>

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
#include "AliPMDPedestal.h"
#include "AliPMDddldata.h"
#include "AliPMDHotData.h"
#include "AliPMDNoiseCut.h"
#include "AliPMDddlinfoData.h"
#include "AliPMDRecoParam.h"
#include "AliRecoParam.h"
#include "AliPMDReconstructor.h"

#include "AliDAQ.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"



ClassImp(AliPMDClusterFinder)

AliPMDClusterFinder::AliPMDClusterFinder():
  fRunLoader(0),
  fPMDLoader(0),
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
  fCalibHot(GetCalibHot()),
  fNoiseCut(GetNoiseCut()),
  fDdlinfo(GetDdlinfoData()),
  fRecoParam(0x0),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fRecpoints(new TClonesArray("AliPMDrecpoint1", 1000)),
  fRechits(new TClonesArray("AliPMDrechit", 1000)),
  fNpoint(0),
  fNhit(0),
  fDetNo(0)
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
  fCalibHot(GetCalibHot()),
  fNoiseCut(GetNoiseCut()),
  fDdlinfo(GetDdlinfoData()),
  fRecoParam(0x0),
  fTreeD(0),
  fTreeR(0),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fRecpoints(new TClonesArray("AliPMDrecpoint1", 1000)),
  fRechits(new TClonesArray("AliPMDrechit", 1000)),
  fNpoint(0),
  fNhit(0),
  fDetNo(0)
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
  fCalibHot(GetCalibHot()),
  fNoiseCut(GetNoiseCut()),
  fDdlinfo(GetDdlinfoData()),
  fRecoParam(0x0),
  fTreeD(0),
  fTreeR(0),
  fDigits(NULL),
  fRecpoints(NULL),
  fRechits(NULL),
  fNpoint(0),
  fNhit(0),
  fDetNo(0)
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

void AliPMDClusterFinder::Digits2RecPoints(TTree *digitsTree,
					   TTree *clustersTree, Int_t gRecoMode)
{
  // Converts digits to recpoints after running clustering
  // algorithm on CPV plane and PREshower plane
  //
  // This algorithm is called during the reconstruction from digits

  Int_t    det  = 0, smn = 0;
  Int_t    xpos = 0, ypos = 0;
  Int_t    ismn = 0;
  Int_t    idet = 0;
  Float_t  adc  = 0.;
  Float_t  clusdata[6] = {0.,0.,0.,0.,0.,0.};

  AliPMDcluster *pmdcl = 0x0;

  TObjArray *pmdcont = new TObjArray();

  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  // Fetch the reco param object

  fRecoParam = AliPMDReconstructor::GetRecoParam();
  if(fRecoParam == 0x0)
    {
       AliFatal("No Reco Param found for PMD!!!");
    }


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

	  if(det < 0 || det > 1)
	    {
	      AliError(Form("*CPV/PRE NUMBER WRONG %d *",det));
	      continue; 
	    }
	  if(smn == -1 || smn > 23)
	    {
	      AliError(Form("*MODULE NUMBER WRONG %d *",smn));
	      continue; 
	    }

	  if(xpos < 0 || xpos > 47 || ypos < 0 || ypos > 95)
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

	  // Hot cell - set the cell adc = 0
	  Float_t hotflag = fCalibHot->GetHotChannel(det,smn,xpos,ypos);
	  if (hotflag == 1.) adc1 = 0;

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

      // Set the minimum noise cut per module before clustering

      // Int_t imod = idet*24 + ismn;


      // Int_t cluspar = fRecoParam->GetPbPbParam()->GetClusteringParam();
      AliPMDRecoParam * par = fRecoParam->GetPPParam();
      Int_t cluspar = par->GetClusteringParam();
      delete par;

      // Int_t cluspar = fRecoParam->GetCosmicParam()->GetClusteringParam();
      
      //_______________________________________________________// 
      //Added to switch Refine and crude Clustering - satya//
      // temporary solution - will be sorted out later
      /*cluspar = 1;
      static AliPMDRecoParam *reconp = NULL;
      reconp = (AliPMDRecoParam*)AliPMDReconstructor::GetRecoParam();
      if(!reconp) {
	cluspar = 1;
      } 
      else { 
	
      if( reconp->GetClusteringParam() == 1) 
	cluspar = 1;
      if( reconp->GetClusteringParam() == 2) 
	cluspar = 2;
	}

      */
      cluspar = gRecoMode;
      //_______________________________________________________// 
      
      pmdclust->SetClusteringParam(cluspar);

      Float_t encut = 4.;
      pmdclust->SetEdepCut(encut);
      pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,pmdcont);
      
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
	  if (ncell > 19) ncell = 19;
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

      branch1->Fill();
      ResetRecpoint();

    } // modules


  ResetCellADC();

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
}
// ------------------------------------------------------------------------- //

void AliPMDClusterFinder::Digits2RecPoints(AliRawReader *rawReader,
					   TTree *clustersTree, Int_t gRecoMode)
{
  // Converts RAW data to recpoints after running clustering
  // algorithm on CPV and PREshower plane
  //
  // This method is called at the time of reconstruction from RAW data


  AliPMDddldata *pmdddl = 0x0;
  AliPMDcluster *pmdcl  = 0x0;

  Float_t  clusdata[6];
  TObjArray pmdddlcont;

  TObjArray *pmdcont = new TObjArray();

  AliPMDClustering *pmdclust = new AliPMDClusteringV1();

  // access the ddlinfo database to fetch  the no of modules per DDL

  Int_t moduleddl[6] = {0,0,0,0,0,0};

  for(Int_t jddl = 0; jddl < 6; jddl++)
    {
      moduleddl[jddl] = fDdlinfo->GetNoOfModulePerDdl(jddl);
    }

  // Set the minimum noise cut per module before clustering

  fRecoParam = AliPMDReconstructor::GetRecoParam();

  if(fRecoParam == 0x0)
    {
       AliFatal("No Reco Param found for PMD!!!");
    }

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
      iSMN = moduleddl[indexDDL];

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


	  if(det < 0 || det > 1)
	    {
	      AliError(Form("*CPV/PRE NUMBER WRONG %d *",det));
	      continue; 
	    }
	  if(smn < 0 || smn > 23)
	    {
	      AliError(Form("*MODULE NUMBER WRONG %d *",smn));
	      continue; 
	    }
	  if(row < 0 || row > 47 || col < 0 || col > 95)
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

	  // Hot cell - set the cell adc = 0
	  Float_t hotflag = fCalibHot->GetHotChannel(det,smn,row,col);
	  if (hotflag == 1.) sig1 = 0;

	  // CALIBRATION
	  Float_t gain = fCalibGain->GetGainFact(det,smn,row,col);
	  //printf("sig = %d gain = %f\n",sig,gain);
	  sig = (Int_t) (sig1*gain);

	  if (indexDDL == 0)
	    {
	      if (det != 0)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      if (iSMN == 6)
		{
		  indexsmn = smn;
		}
	      else if (iSMN == 12)
		{
		  if (smn < 6)
		    indexsmn = smn;
		  else if (smn >= 18 && smn < 24)
		    indexsmn = smn-12;
		}
	    }
	  else if (indexDDL >= 1 && indexDDL < 4)
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
	      if (smn >= 6 && smn < 18)
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
	  
	  if (indexDDL == 0)
	    {
	      if (iSMN == 6)
		{
		  ismn = indexsmn;
		}
	      else if (iSMN == 12)
		{
		  
		  if (indexsmn < 6)
		    ismn = indexsmn;
		  else if (indexsmn >= 6 && indexsmn < 12)
		    ismn = indexsmn + 12;
		}
	      idet = 0;
	    }
	  else if (indexDDL >= 1 && indexDDL < 4)
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
	      ismn = indexsmn + 6;
	      idet = 1;
	    }

	  if (totAdcMod <= 0) continue;

	  Int_t imod = idet*24 + ismn;

	  // Int_t cluspar = fRecoParam->GetPbPbParam()->GetClusteringParam();
	  AliPMDRecoParam * par = fRecoParam->GetPPParam();
	   Int_t cluspar = par->GetClusteringParam();
	   delete par;
	  // Int_t cluspar = fRecoParam->GetCosmicParam()->GetClusteringParam();

	  //_______________________________________________________// 
	  //Added to switch Refine and crude Clustering - satya//
	  // temporary solution - will be sorted out later
	  /* cluspar = 1;
	  static AliPMDRecoParam *reconp = NULL;
	  reconp = (AliPMDRecoParam*)AliPMDReconstructor::GetRecoParam();
	  if(!reconp) {
	    cluspar = 1;
	  } 
	  else { 
	    if( reconp->GetClusteringParam() == 1) 
	      cluspar = 1;
	    if( reconp->GetClusteringParam() == 2) 
	      cluspar = 2;
	  }*/


	  cluspar = gRecoMode; // permanent solution

	  //_______________________________________________________// 

	  pmdclust->SetClusteringParam(cluspar);
	  Float_t encut = fNoiseCut->GetNoiseCut(imod);

	  pmdclust->SetEdepCut(encut);
	  pmdclust->DoClust(idet,ismn,fCellTrack,fCellPid,fCellADC,pmdcont);

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
	      if (ncell > 19) ncell = 19;
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
//--------------------------------------------------------------------//
AliPMDHotData* AliPMDClusterFinder::GetCalibHot() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Hot");
  
  if(!entry) AliFatal("HotData object retrieval failed!");
  
  AliPMDHotData *hot = 0;
  if (entry) hot = (AliPMDHotData*) entry->GetObject();
  
  if (!hot)  AliFatal("No hot data from  database !");
  
  return hot;
}
//--------------------------------------------------------------------//
AliPMDNoiseCut* AliPMDClusterFinder::GetNoiseCut() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/NoiseCut");
  
  if(!entry) AliFatal("Noisecut object retrieval failed!");
  
  AliPMDNoiseCut *ncut = 0;
  if (entry) ncut = (AliPMDNoiseCut*) entry->GetObject();
  
  if (!ncut)  AliFatal("No noise cut data from  database !");
  
  return ncut;
}
//--------------------------------------------------------------------//
AliPMDddlinfoData* AliPMDClusterFinder::GetDdlinfoData() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Ddlinfo");
  
  if(!entry) AliFatal("ddlinfo object retrieval failed!");
  
  AliPMDddlinfoData *ddlinfo = 0;
  if (entry) ddlinfo = (AliPMDddlinfoData*) entry->GetObject();
  
  if (!ddlinfo)  AliFatal("No ddl info data from  database !");
  
  return ddlinfo;
}
  
