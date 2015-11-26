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
  TObject(),
  fRunNumber(-1),
  fEventNumber(0)
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
		}
	    }
	}
    }

    for (int i = 0; i < 6; i++)
    {
	for (int j = 0; j < 51; j++)
	{
	    for (int k = 0; k < 25; k++)
	    {
		for (int l = 0; l < 64; l++)
		{
		    fPedChain[i][j][k][l] = -1;
		}
	    }
	}
    }


}
//_____________________________________________________________________
AliPMDCalibPedestal::AliPMDCalibPedestal(const AliPMDCalibPedestal &ped) :
  TObject(ped),
  fRunNumber(ped.fRunNumber),
  fEventNumber(ped.fEventNumber)
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
		}
	    }
	}
    }

    for (int i = 0; i < 6; i++)
    {
	for (int j = 0; j < 51; j++)
	{
	    for (int k = 0; k < 25; k++)
	    {
		for (int l = 0; l < 64; l++)
		{
		    fPedChain[i][j][k][l] = -1;
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
Bool_t AliPMDCalibPedestal::ProcessEvent(AliRawReader *rawReader, TObjArray *pmdddlcont)
{
  //
  //  Event processing loop - AliRawReader
  //

    const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");

    UInt_t detsmnrowcol = 0;
    UInt_t pbus = 0, mcm = 0, chno = 0;

    fRunNumber = rawReader->GetRunNumber();

    AliPMDRawStream rawStream(rawReader);

    fEventNumber++;

    Int_t iddl = -1;
    Int_t numberofDDLs = 0;

    while ((iddl = rawStream.DdlData(pmdddlcont)) >=0) {
      numberofDDLs++;
      Int_t ientries = pmdddlcont->GetEntries();
      //printf("iddl = %d ientries = %d\n",iddl, ientries);
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	    AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont->UncheckedAt(ient);
	    
	    Int_t det = pmdddl->GetDetector();
	    Int_t smn = pmdddl->GetSMN();
	    Int_t row = pmdddl->GetRow();
	    Int_t col = pmdddl->GetColumn();
	    Float_t sig = (Float_t) pmdddl->GetSignal();

	    pbus = (UInt_t) pmdddl->GetPatchBusId();
	    mcm  = (UInt_t) pmdddl->GetMCM();
	    chno = (UInt_t) pmdddl->GetChannel();

	    detsmnrowcol = 0;
	    AliBitPacking::PackWord(det,detsmnrowcol,0,7);
	    AliBitPacking::PackWord(smn,detsmnrowcol,8,15);
	    AliBitPacking::PackWord(row,detsmnrowcol,16,23);
	    AliBitPacking::PackWord(col,detsmnrowcol,24,31);

	    if (fPedChain[iddl][pbus][mcm][chno] == -1)
		fPedChain[iddl][pbus][mcm][chno] = (Int_t)detsmnrowcol;


	    fPedVal[det][smn][row][col]   += sig;
	    fPedValSq[det][smn][row][col] += sig*sig;
	    fPedCount[det][smn][row][col]++;
	}
      pmdddlcont->Delete();
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

    UInt_t  detsmnrowcol = 0;
    Int_t   det = 0, sm = 0, row = 0, col = 0;
    Int_t   idet = 0, ism = 0, irow = 0, icol = 0;
    Float_t mean = 0., rms = 0., meanoff = 0.;
    Double_t meansq = 0., diff = 0.;

    FILE *fpw0 = fopen("pedestal2304.ped","w");
    FILE *fpw1 = fopen("pedestal2305.ped","w");
    FILE *fpw2 = fopen("pedestal2306.ped","w");
    FILE *fpw3 = fopen("pedestal2307.ped","w");
    FILE *fpw4 = fopen("pedestal2308.ped","w");
    FILE *fpw5 = fopen("pedestal2309.ped","w");

    fprintf(fpw0,"//=============================================\n");
    fprintf(fpw0,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw0,"//=============================================\n");
    fprintf(fpw0,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw0,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw0,"//---------------------------------------------\n");
    fprintf(fpw0,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");
    fprintf(fpw0,"//---------------------------------------------\n");

    fprintf(fpw1,"//=============================================\n");
    fprintf(fpw1,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw1,"//=============================================\n");
    fprintf(fpw1,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw1,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw1,"//---------------------------------------------\n");
    fprintf(fpw1,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");

    fprintf(fpw2,"//=============================================\n");
    fprintf(fpw2,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw2,"//=============================================\n");
    fprintf(fpw2,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw2,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw2,"//---------------------------------------------\n");
    fprintf(fpw2,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");
    fprintf(fpw2,"//---------------------------------------------\n");

    fprintf(fpw3,"//=============================================\n");
    fprintf(fpw3,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw3,"//=============================================\n");
    fprintf(fpw3,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw3,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw3,"//---------------------------------------------\n");
    fprintf(fpw3,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");
    fprintf(fpw3,"//---------------------------------------------\n");

    fprintf(fpw4,"//=============================================\n");
    fprintf(fpw4,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw4,"//=============================================\n");
    fprintf(fpw4,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw4,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw4,"//---------------------------------------------\n");
    fprintf(fpw4,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");
    fprintf(fpw4,"//---------------------------------------------\n");

    fprintf(fpw5,"//=============================================\n");
    fprintf(fpw5,"//     Pedestal file Calculated by Online DA\n");
    fprintf(fpw5,"//=============================================\n");
    fprintf(fpw5,"//     RUN             :%d\n",fRunNumber);
    fprintf(fpw5,"//     Statistics      :%d\n",fEventNumber);
    fprintf(fpw5,"//---------------------------------------------\n");
    fprintf(fpw5,"//format:CHAIN_NO  FEE_ID  CHANNEL  MEAN  SIGMA\n");
    fprintf(fpw5,"//---------------------------------------------\n");


    for(Int_t iddl = 0; iddl < 6; iddl++)
      {
	for(Int_t ibus = 1; ibus < 51; ibus++)
	  {
	    for(Int_t imcm = 1; imcm < 25; imcm++)
	      {
		for(Int_t ich = 0; ich < 64; ich++)
		  {

		    if (fPedChain[iddl][ibus][imcm][ich] != -1)
		      {
			detsmnrowcol = (UInt_t)fPedChain[iddl][ibus][imcm][ich];

			idet =  detsmnrowcol & 0x00FF;
			ism  = (detsmnrowcol >> 8) & 0x00FF;
			irow = (detsmnrowcol >> 16) & 0x00FF;
			icol = (detsmnrowcol >> 24) & 0x00FF;
			
			mean = 0.;
			rms  = 0.;
			if (fPedCount[idet][ism][irow][icol] > 0)
			  {
			    mean = fPedVal[idet][ism][irow][icol]/fPedCount[idet][ism][irow][icol];
			    // Fix for 2015 data sjena
                            // should be used for pedestal
			    meanoff = mean; // Variable kept to add an offset - it is off now   
			    meansq  = fPedValSq[idet][ism][irow][icol]/fPedCount[idet][ism][irow][icol];
			    
			    diff = meansq - mean*mean;
			    if (diff > 0.)
			      {
				rms  = sqrt(diff);
			      }
			    else
			      {
				rms = 0.;
			      }

			    if (iddl == 0)
			      {
				fprintf(fpw0,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
			    else if (iddl == 1)
			      {
				fprintf(fpw1,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
			    else if (iddl == 2)
			      {
				fprintf(fpw2,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
			    else if (iddl == 3)
			      {
				fprintf(fpw3,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
			    else if (iddl == 4)
			      {
				fprintf(fpw4,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
			    else if (iddl == 5)
			      {
				fprintf(fpw5,"%d %d %d %f %f\n",
					ibus, imcm, ich, meanoff, rms);
			      }
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

    pedtree->Branch("det",&det,"det/I");
    pedtree->Branch("sm",&sm,"sm/I");
    pedtree->Branch("row",&row,"row/I");
    pedtree->Branch("col",&col,"col/I");
    pedtree->Branch("mean",&mean,"mean/F");
    pedtree->Branch("rms",&rms,"rms/F");

    for (idet = 0; idet < kDet; idet++)
      {
	for (ism = 0; ism < kMaxSMN; ism++)
	  {
	    for (irow = 0; irow < kMaxRow; irow++)
	      {
		for (icol = 0; icol < kMaxCol; icol++)
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
		      }
		    
		  }
	      }
	  }
      }
}
// -------------------------------------------------------------------
