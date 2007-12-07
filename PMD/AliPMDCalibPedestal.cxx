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

    AliPMDRawStream rawStream(rawReader);

    TObjArray pmdddlcont;
    Bool_t streamout = kTRUE;

    for (Int_t iddl = 0; iddl < kDDL; iddl++)
    {
	rawReader->Select("PMD", iddl, iddl);
	//cout << reader.GetDataSize() << endl;
	streamout = rawStream.DdlData(iddl, &pmdddlcont);
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
	    Float_t sig = (Float_t) pmdddl->GetSignal();

	    fPedVal[det][smn][row][col]   += sig;
	    fPedValSq[det][smn][row][col] += sig*sig;
	    fPedCount[det][smn][row][col]++;
	}
	pmdddlcont.Clear();
    }
    return streamout;
}
//_____________________________________________________________________

void AliPMDCalibPedestal::Analyse(TTree *pedtree)
{
    //
    //  Calculate pedestal Mean and RMS
    //
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
	    for (int irow = 0; irow < kMaxRow; irow++)
	    {
		for (int icol = 0; icol < kMaxCol; icol++)
		{
		    det  = idet;
		    sm   = ism;
		    row  = irow;
		    col  = icol;
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
		    }

		    pedtree->Fill();
		}
	    }
	}
    }
}
//_____________________________________________________________________
