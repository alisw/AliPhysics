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


//Root includes
#include <TObjArray.h>
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
#include "TTreeStream.h"

//AliRoot includes
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"
#include "AliBitPacking.h"
//header file
#include "AliPMDCalibPedestal.h"


ClassImp(AliPMDCalibPedestal)


AliPMDCalibPedestal::AliPMDCalibPedestal() :
  TObject()
{
    //
    // default constructor
    //

    for (int i = 0; i < kDet; i++)
    {
	for (int j = 0; j < kMaxSMN; j++)
	{
	    for (int k = 0; k < kMaxRow; k++)
	    {
		for (int l = 0; l < kMaxCol; l++)
		{
		    fPedVal[i][j][k][l]   = 0.;
		    fPedValSq[i][j][k][l] = 0.;
		    fPedCount[i][j][k][l] = 0.;
		    fPedChain[i][j][k][l] = 0;
		}
	    }
	}
    }


}
//_____________________________________________________________________
AliPMDCalibPedestal::AliPMDCalibPedestal(const AliPMDCalibPedestal &ped) :
  TObject(ped)
{
    //
    // copy constructor
    //
    for (int i = 0; i < kDet; i++)
    {
	for (int j = 0; j < kMaxSMN; j++)
	{
	    for (int k = 0; k < kMaxRow; k++)
	    {
		for (int l = 0; l < kMaxCol; l++)
		{
		    fPedVal[i][j][k][l]   = ped.fPedVal[i][j][k][l];
		    fPedValSq[i][j][k][l] = ped.fPedValSq[i][j][k][l];
		    fPedCount[i][j][k][l] = ped.fPedCount[i][j][k][l];
		    fPedChain[i][j][k][l] = ped.fPedChain[i][j][k][l];
		}
	    }
	}
    }
    
}
//_____________________________________________________________________
AliPMDCalibPedestal& AliPMDCalibPedestal::operator = (const  AliPMDCalibPedestal &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliPMDCalibPedestal(source);

  return *this;
}
//_____________________________________________________________________
AliPMDCalibPedestal::~AliPMDCalibPedestal()
{
    //
    // destructor
    //
}
//_____________________________________________________________________
Bool_t AliPMDCalibPedestal::ProcessEvent(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //

    const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");

    UInt_t busmcmch;
    UInt_t pbus, mcm, chno;

    AliPMDRawStream rawStream(rawReader);

    TObjArray pmdddlcont;
    Int_t iddl = -1;
    Int_t numberofDDLs = 0;

    while ((iddl = rawStream.DdlData(&pmdddlcont)) >=0) {
      numberofDDLs++;
      Int_t ientries = pmdddlcont.GetEntries();
      //printf("iddl = %d ientries = %d\n",iddl, ientries);
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	    AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	    
	    Int_t det = pmdddl->GetDetector();
	    Int_t smn = pmdddl->GetSMN();
	    Int_t row = pmdddl->GetRow();
	    Int_t col = pmdddl->GetColumn();
	    Float_t sig = (Float_t) pmdddl->GetSignal();

	    pbus = (UInt_t) pmdddl->GetPatchBusId();
	    mcm  = (UInt_t) pmdddl->GetMCM();
	    chno = (UInt_t) pmdddl->GetChannel();

	    busmcmch = 0;
	    AliBitPacking::PackWord(chno,busmcmch,0,7);
	    AliBitPacking::PackWord(mcm,busmcmch,8,15);
	    AliBitPacking::PackWord(pbus,busmcmch,16,23);


	    if (fPedChain[det][smn][row][col] == 0)
		fPedChain[det][smn][row][col]   = busmcmch;

	    fPedVal[det][smn][row][col]   += sig;
	    fPedValSq[det][smn][row][col] += sig*sig;
	    fPedCount[det][smn][row][col]++;
	}
      pmdddlcont.Delete();
    }
    if (numberofDDLs < kDDL)
      return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________

void AliPMDCalibPedestal::Analyse(TTree *pedtree)
{
    //
    //  Calculate pedestal Mean and RMS
    //

    FILE *fpw0 = fopen("pedestal2304.ped","w");
    FILE *fpw1 = fopen("pedestal2305.ped","w");
    FILE *fpw2 = fopen("pedestal2306.ped","w");
    FILE *fpw3 = fopen("pedestal2307.ped","w");
    FILE *fpw4 = fopen("pedestal2308.ped","w");
    FILE *fpw5 = fopen("pedestal2309.ped","w");

    UInt_t  busmcmch;
    UInt_t  pbus, mcm, chno;
    Int_t   ddlno;
    Int_t   det, sm, row, col;
    Float_t mean, rms;
    Float_t meansq, diff;



    pedtree->Branch("det",&det,"det/I");
    pedtree->Branch("sm",&sm,"sm/I");
    pedtree->Branch("row",&row,"row/I");
    pedtree->Branch("col",&col,"col/I");
    pedtree->Branch("mean",&mean,"mean/F");
    pedtree->Branch("rms",&rms,"rms/F");

    for (int idet = 0; idet < kDet; idet++)
    {
	for (int ism = 0; ism < kMaxSMN; ism++)
	{
	    ConvertDDL(idet,ism,ddlno);
	    for (int irow = 0; irow < kMaxRow; irow++)
	    {
		for (int icol = 0; icol < kMaxCol; icol++)
		{
		    det  = idet;
		    sm   = ism;
		    row  = irow;
		    col  = icol;
		    mean = 0.;
		    rms  = 0.;

		    if (fPedCount[idet][ism][irow][icol] > 0)
		    {
		      mean = fPedVal[idet][ism][irow][icol]/fPedCount[idet][ism][irow][icol];

		      meansq = fPedValSq[idet][ism][irow][icol]/fPedCount[idet][ism][irow][icol];

			diff = meansq - mean*mean;
			if (diff > 0.)
			{
			    rms  = sqrt(diff);
			}
			else
			{
			    rms = 0.;
			}
			pedtree->Fill();

			busmcmch = fPedChain[idet][ism][irow][icol];

			chno = busmcmch & 0x00FF;
			mcm  = (busmcmch >> 8) & 0x00FF;
			pbus = (busmcmch >> 16) & 0x00FF;
			
			if (ddlno == 0)
			{
			    fprintf(fpw0,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			else if (ddlno == 1)
			{
			    fprintf(fpw1,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			else if (ddlno == 2)
			{
			    fprintf(fpw2,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			else if (ddlno == 3)
			{
			    fprintf(fpw3,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			else if (ddlno == 4)
			{
			    fprintf(fpw4,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			else if (ddlno == 5)
			{
			    fprintf(fpw5,"%d %d %d %f %f\n",
				    pbus, mcm, chno, mean, rms);
			}
			
		    }
		    
		}
	    }
	}
    }
    
    fclose(fpw0);
    fclose(fpw1);
    fclose(fpw2);
    fclose(fpw3);
    fclose(fpw4);
    fclose(fpw5);
}


// -------------------------------------------------------------------

void AliPMDCalibPedestal::ConvertDDL(Int_t det, Int_t smn, Int_t &ddlno)
{
// Given the plane number and serial module number, ddlno is calculated

    if (det == 0)
    {
	if (smn <= 5)
	{
	    ddlno = 0;
	}
	else if (smn > 5 && smn <= 11)
	{
	    ddlno = 1;
	}
	else if (smn > 11 && smn <= 17)
	{
	    ddlno = 2;
	}
	else if (smn > 17 && smn <= 23)
	{
	    ddlno = 3;
	}
    }
    else if (det == 1)
    {
	if (smn <= 5 || (smn >= 18 && smn <=23))
	{
	    ddlno = 4;
	}
	else if (smn >= 6 && smn <= 17)
	{
	    ddlno = 5;
	}
    }

}

