//////////////////////////////////////////////////////////////////////////////
//
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.  	    *
// *									    *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is	    *
// * provided "as is" without express or implied warranty.		    *
// **************************************************************************/
//
//////////////////////////////////////////////////////////////////////////////

#include "TF1.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"
#include "AliPMDCalibGain.h"

ClassImp(AliPMDCalibGain)

AliPMDCalibGain::AliPMDCalibGain(): TObject()
{
  // Standard Constructor
    for(Int_t idet = 0; idet < kDet; idet++)
    {
	for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	    fHsmIso[idet][ismn] = new TH1F(Form("HmsIso_%d_%d",idet,ismn),"",100,0.,1000.);
	    for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
		for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		    fHadcIso[idet][ismn][jrow][kcol]  = new TH1F(Form("HadcIso_%d_%d_%d_%d",idet,ismn,jrow,kcol),"",100,0.,1000.);
		}
	    }
	}
    }
}
// ------------------------------------------------------------------------ //
AliPMDCalibGain::AliPMDCalibGain(const AliPMDCalibGain &pmdcalibgain):
    TObject(pmdcalibgain)
{
    for(Int_t idet = 0; idet < kDet; idet++)
    {
	for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	    fHsmIso[idet][ismn] = pmdcalibgain.fHsmIso[idet][ismn] ;
	    for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
		for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		    fHadcIso[idet][ismn][jrow][kcol]  = pmdcalibgain.fHadcIso[idet][ismn][jrow][kcol];
		}
	    }
	}
    }
    
}
// ------------------------------------------------------------------------ //
AliPMDCalibGain &AliPMDCalibGain::operator=(const AliPMDCalibGain &pmdcalibgain)
{
    if(this != &pmdcalibgain)
    {
	for(Int_t idet = 0; idet < kDet; idet++)
	{
	    for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	    {
		fHsmIso[idet][ismn] = pmdcalibgain.fHsmIso[idet][ismn] ;
		for(Int_t jrow = 0; jrow < kMaxRow;jrow++)
		{
		    for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		    {
			fHadcIso[idet][ismn][jrow][kcol]  =
			    pmdcalibgain.fHadcIso[idet][ismn][jrow][kcol];
		    }
		}
	    }
	}
    }
    return *this;
}
// ------------------------------------------------------------------------ //
AliPMDCalibGain::~AliPMDCalibGain()
{
    // dtor
    for(Int_t idet = 0; idet < kDet; idet++)
    {
	for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  delete fHsmIso[idet][ismn];
	    for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
		for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		  delete fHadcIso[idet][ismn][jrow][kcol];
		}
	    }
	}
    }
}
// ------------------------------------------------------------------------ //
Bool_t AliPMDCalibGain::ProcessEvent(AliRawReader *rawReader)
{
  // Calculates the ADC of isolated cell

  TObjArray pmdddlcont;

  const Int_t kDDL           = AliDAQ::NumberOfDdls("PMD");
  const Int_t kCellNeighbour = 6;

  Int_t neibx[6] = {1,0,-1,-1,0,1};
  Int_t neiby[6] = {0,1,1,0,-1,-1};
  
  Int_t id1,jd1;  //neighbour row/col
  Int_t isocount; //number of neighbours with 0 signal

  Int_t d1[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Bool_t streamout = kFALSE;

  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
        {
          for(Int_t irow = 0; irow < kMaxRow; irow++)
            {
              for(Int_t icol = 0; icol < kMaxCol; icol++)
                {
                  d1[idet][ismn][irow][icol] = 0;
                }
	    }
	}
    }

  AliPMDRawStream rawStream(rawReader);

  for (Int_t iddl = 0; iddl < kDDL; iddl++)
  {
      rawReader->Select("PMD", iddl, iddl);
      streamout = rawStream.DdlData(iddl,&pmdddlcont);
      
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
      {
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t idet = pmdddl->GetDetector();
	  Int_t ismn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t ichno = pmdddl->GetChannel();
	  Int_t irow = pmdddl->GetRow();
	  Int_t icol = pmdddl->GetColumn();
	  Int_t isig = pmdddl->GetSignal();
	  
	  if (isig>0)
	  {
	      d1[idet][ismn][irow][icol] = isig;
	  }
      }
      pmdddlcont.Clear();
  }
  
  for(Int_t idet=0; idet < kDet; idet++)
  {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
      {
	  for(Int_t irow = 0; irow < kMaxRow; irow++)
	  {
	      for(Int_t icol = 0; icol < kMaxCol; icol++)
	      {
		  if(d1[idet][ismn][irow][icol] > 0)
		  {
		      isocount = 0;
		      for(Int_t ii = 0; ii < kCellNeighbour; ii++)
		      {
			  id1 = irow + neibx[ii];
			  jd1 = icol + neiby[ii];
			  if(d1[idet][ismn][id1][jd1] == 0)
			  {
			      isocount++;
			      if(isocount == kCellNeighbour)
			      {
				  fHsmIso[idet][ismn]->Fill(d1[idet][ismn][irow][icol]);
				  fHadcIso[idet][ismn][irow][icol]->Fill(d1[idet][ismn][irow][icol]);
			      }
			  }
		      }  // neigh cell cond.
		  }     // d>0 cond.
	      }
	  }
      }
      
  }
  return streamout;
}
// ------------------------------------------------------------------------ //
void AliPMDCalibGain::Analyse(TTree *gaintree)
{
    // Calculates the mean
    Int_t   DET, SM, ROW, COL;
    Float_t GAIN;
    Float_t modmean, cellmean;

    gaintree->Branch("DET",&DET,"DET/I");
    gaintree->Branch("SM",&SM,"SM/I");
    gaintree->Branch("ROW",&ROW,"ROW/I");
    gaintree->Branch("COL",&COL,"COL/I");
    gaintree->Branch("GAIN",&GAIN,"GAIN/F");

    for(Int_t idet = 0; idet < kDet; idet++)
    {
	for(Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	    modmean = fHsmIso[idet][ism]->GetMean();
	    for(Int_t irow = 0; irow < kMaxRow; irow++)
	    {
		for(Int_t icol = 0; icol < kMaxCol; icol++)
		{
		    cellmean = fHadcIso[idet][ism][irow][icol]->GetMean();
		    DET      = idet;
		    SM       = ism;
		    ROW      = irow;
		    COL      = icol;
		    GAIN     = 0.;

		    if(modmean > 0.0)
		    {
			GAIN = cellmean/modmean;
		    }
		    gaintree->Fill();
		}
	    }
	}
    }
    
}
// ------------------------------------------------------------------------ //
