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

    for (int i = 0; i < 2; i++)
    {
	for (int j = 0; j < 24; j++)
	{
	    for (int k = 0; k < 96; k++)
	    {
		for (int l = 0; l < 96; l++)
		{

		    fPedHisto[i][j][k][l] = new TH1F("","",300,0.,300.);
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
    for (int i = 0; i < 2; i++)
    {
	for (int j = 0; j < 24; j++)
	{
	    for (int k = 0; k < 96; k++)
	    {
		for (int l = 0; l < 96; l++)
		{
		    
		    fPedHisto[i][j][k][l] = ped.fPedHisto[i][j][k][l];
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
    delete fPedHisto;
}
//_____________________________________________________________________
Bool_t AliPMDCalibPedestal::ProcessEvent(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
    AliPMDRawStream rawStream(rawReader);

    TObjArray pmdddlcont;
    Bool_t streamout = kTRUE;

    for (Int_t iddl = 0; iddl < 6; iddl++)
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
	    Int_t sig = pmdddl->GetSignal();

	    fPedHisto[det][smn][row][col]->Fill((Float_t) sig);
	    
	}
	pmdddlcont.Clear();
    }
    return streamout;
}
//_____________________________________________________________________

void AliPMDCalibPedestal::Analyse()
{
    //
    //  Calculate pedestal Mean and RMS
    //
    for (int i = 0; i < 2; i++)
    {
	for (int j = 0; j < 24; j++)
	{
	    for (int k = 0; k < 96; k++)
	    {
		for (int l = 0; l < 96; l++)
		{

		    Float_t mean = fPedHisto[i][j][k][l]->GetMean();
		    Float_t rms  = fPedHisto[i][j][k][l]->GetRMS();
		}
	    }
	}
    }
}
//_____________________________________________________________________
