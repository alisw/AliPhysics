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
// Author : Z. Ahamed

#include "TF1.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TObjArray.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRawReaderFile.h"
#include "AliPMDCalibrator.h"
#include "AliRawReaderDate.h"
#include "AliPMDRawStream.h"
#include "AliPMDCalibData.h"
#include "AliPMDddldata.h"
#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliDAQ.h"

ClassImp(AliPMDCalibrator)

const Int_t kDet =2;
const Int_t kMaxSMN = 24;
const Int_t kMaxRow =96;
const Int_t kMaxCol =96;

AliPMDCalibrator::AliPMDCalibrator():
  fCalibData(new AliPMDCalibData())
{
  // Standard Constructor
  for(Int_t d=0;d<2;d++)
    {
      for(Int_t i=0;i<24;i++)
	{
	  fHsmIso[d][i] = NULL ;
	  for(Int_t j=0;j<96;j++)
	    {
	      for(Int_t k=0;k<96;k++)
		{
		  fGainFact[d][i][j][k] = 0.0;
		  fHadcIso[d][i][j][k]  = NULL;
		}
	    }
	}
    }
}
// ------------------------------------------------------------------------ //
AliPMDCalibrator::AliPMDCalibrator(const AliPMDCalibrator &pmdcalibrator):
  fCalibData(new AliPMDCalibData())
{
  for(Int_t d=0;d<2;d++)
    {
      for(Int_t i=0;i<24;i++)
	{
	  fHsmIso[d][i] = pmdcalibrator.fHsmIso[d][i] ;
	  for(Int_t j=0;j<96;j++)
	    {
	      for(Int_t k=0;k<96;k++)
		{
		  fGainFact[d][i][j][k] = pmdcalibrator.fGainFact[d][i][j][k];
		  fHadcIso[d][i][j][k]  = pmdcalibrator.fHadcIso[d][i][j][k];
		}
	    }
	}
    }

}
// ------------------------------------------------------------------------ //
AliPMDCalibrator &AliPMDCalibrator::operator=(const AliPMDCalibrator &pmdcalibrator)
{
  if(this != &pmdcalibrator)
    {
      for(Int_t d=0;d<2;d++)
	{
	  for(Int_t i=0;i<24;i++)
	    {
	      fHsmIso[d][i] = pmdcalibrator.fHsmIso[d][i] ;
	      for(Int_t j=0;j<96;j++)
		{
		  for(Int_t k=0;k<96;k++)
		    {
		      fGainFact[d][i][j][k] =
			pmdcalibrator.fGainFact[d][i][j][k];
		      fHadcIso[d][i][j][k]  =
			pmdcalibrator.fHadcIso[d][i][j][k];
		    }
		}
	    }
	}
    }
  return *this;
}
// ------------------------------------------------------------------------ //
AliPMDCalibrator::~AliPMDCalibrator()
{
  // dtor
  if(fHsmIso)  delete fHsmIso ;
  if(fHadcIso) delete fHadcIso ;
  delete fCalibData;
}
// ------------------------------------------------------------------------ //

void AliPMDCalibrator::Exec()
{
  // reads parameters and does the calibration
  CalculateIsoCell() ;

}
// ------------------------------------------------------------------------ //

void AliPMDCalibrator::Init()
{
  // intializes everything
  char hname[kMaxSMN];
  char hname24[kMaxSMN];
  char hnameiso[120];
  char htitle1[120];

  for(Int_t d=0;d<2;d++) {
    for(Int_t i1=0; i1<kMaxSMN;i1++) {
      sprintf(hname,"det_%d_iso_sm_%2d",d,i1);
      sprintf(hname24,"det_%d_iso_sm_%2d",d,i1);
      fHsmIso[d][i1]= new TH1F(hname,hname24,100,0,1000);
      for(Int_t j1 = 0; j1 < kMaxRow; j1++) {
	for(Int_t k1 = 0; k1 < kMaxCol; k1++) {
	  sprintf(hnameiso,"Isolated Cell ADC for det_%d_cell_sm%d_row%d_col%d"
		  ,d,i1,j1,k1);
	  sprintf(htitle1,"Isolated Cell ADC for det_%d_cell_sm%d_row%d_col%d"
		  ,d,i1,j1,k1);
	  
	  TObject *old=gDirectory->GetList()->FindObject(hnameiso);
	  if (old) gDirectory->GetList()->Remove(old);
	  fHadcIso[d][i1][j1][k1] = new TH1F(hnameiso,htitle1,100,0.,4000.);
	}
      }
    }
  }
  
}

// ------------------------------------------------------------------------ //

void AliPMDCalibrator::CalculateIsoCell()
{
  // Calculates the ADC of isolated cell

  TObjArray pmdddlcont;
  const Int_t kDDL           = AliDAQ::NumberOfDdls("PMD");
  const Int_t kMaxHit        = 60000;
  const Int_t kCellNeighbour = 6;

  Int_t neibx[6] = {1,0,-1,-1,0,1};
  Int_t neiby[6] = {0,1,1,0,-1,-1};
  
  Int_t id1,jd1; //neighbour row/col
  Int_t countisocell = 0 ;//number of isilated cell
  Int_t isocount; //number of neighbours with 0 signal
  Int_t d1[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Int_t ch[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Int_t maxhit;
  
  
  Int_t det[kMaxHit],smn[kMaxHit];
  Int_t row[kMaxHit],col[kMaxHit],sig[kMaxHit],chno[kMaxHit];

  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
        {
          for(Int_t irow = 0; irow < kMaxRow; irow++)
            {
              for(Int_t icol = 0; icol < kMaxCol; icol++)
                {
                  d1[idet][ismn][irow][icol] = 0;
                  ch[idet][ismn][irow][icol] = 0;
                }
	    }
	}
    }
  //accessing raw data
  AliRawReaderFile reader(".");
  AliPMDRawStream stream(&reader);
  while(reader.NextEvent())
    {
      // printf("In CalculateIsoCell before while(stream.Next()), ...\n");
      
      /*
      while(stream.Next())
	{
	  Int_t idet = stream.GetDetector();
	  Int_t ismn = stream.GetSMN();
	  Int_t ichno = stream.GetChannel();
	  Int_t irow = stream.GetRow();
	  Int_t icol = stream.GetColumn();
	  Int_t isig = stream.GetSignal();
	  
	  if (isig>0)
	    {
	      d1[idet][ismn][irow][icol] = isig;
	      ch[idet][ismn][irow][icol] = ichno;
	    }
	}
      */
      // New PMD Reader is plugged in
      
      for (Int_t iddl = 0; iddl < kDDL; iddl++)
	{
	  reader.Select("PMD", iddl, iddl);
	  stream.DdlData(iddl,&pmdddlcont);
	  
	  Int_t ientries = pmdddlcont.GetEntries();
	  for (Int_t ient = 0; ient < ientries; ient++)
	    {
	      AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	      
	      Int_t idet = pmdddl->GetDetector();
	      Int_t ismn = pmdddl->GetSMN();
	      //Int_t mcm = pmdddl->GetMCM();
	      Int_t ichno = pmdddl->GetChannel();
	      Int_t irow = pmdddl->GetRow();
	      Int_t icol = pmdddl->GetColumn();
	      Int_t isig = pmdddl->GetSignal();
	      
	      if (isig>0)
		{
		  d1[idet][ismn][irow][icol] = isig;
		  ch[idet][ismn][irow][icol] = ichno;
		}
	    }
	  pmdddlcont.Clear();
	}


      maxhit = 0;
      
      for(Int_t idet=0; idet < kDet; idet++)
	{
	  for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	    {
	      for(Int_t irow = 0; irow < kMaxRow; irow++)
		{
		  for(Int_t icol = 0; icol < kMaxCol; icol++)
		    {
		      
		      //printf("d1[%d][%d][%d][%d]=%d\n",
		      //det,ksmn,irow,jcol, d1[det][ksmn][irow][jcol] ); 
		      //printf("ch[%d][%d][%d][%d]=%d\n",
		      //det,ksmn,irow,jcol, ch[det][ksmn][irow][jcol] ); 
		      
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
				      countisocell++;
				      det[maxhit]  = idet;
				      smn[maxhit]  = ismn;
				      row[maxhit]  = irow;
				      col[maxhit]  = icol;
				      sig[maxhit]  = d1[idet][ismn][irow][icol];
				      chno[maxhit] = ch[idet][ismn][irow][icol];
				      maxhit++;
				      fHsmIso[idet][ismn]->Fill(d1[idet][ismn][irow][icol]);
				      fHadcIso[idet][ismn][irow][icol]->Fill(d1[idet][ismn][irow][icol]);
				    }
				}
			    }  // neigh cell cond.
			}     // d>0 cond.
		    }
		}
	    }
	} //event loop
    }
  Double_t histMean[2][24];
  Double_t isoMean[2][24][96][96];
  for(Int_t d1=0;d1<2;d1++)
    {
      for(Int_t i1=0;i1<24;i1++)
	{
	  histMean[d1][i1]= fHsmIso[d1][i1]->GetMean();
	  for(Int_t j1=0;j1<96;j1++)
	    {
	      for(Int_t k1=0;k1<96;k1++)
		{
		  isoMean[d1][i1][j1][k1]=fHadcIso[d1][i1][j1][k1]->GetMean();
		  if(isoMean[d1][i1][j1][k1]>0.0 && histMean[d1][i1]>0.0)
		    {
		      fGainFact[d1][i1][j1][k1]=isoMean[d1][i1][k1][j1]/histMean[d1][i1];
		      Float_t gain=fGainFact[d1][i1][j1][k1];
		      fCalibData->SetGainFact(d1,i1,j1,k1,gain);
		    }                              
		}
	    }
	}
    }
  
}
// ------------------------------------------------------------------------ //
Bool_t AliPMDCalibrator::Store()
{
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBId id("PMD/Calib/Data",0,0);
  AliCDBMetaData md;
  md.SetResponsible("Zubayer");
  md.SetBeamPeriod(0);
  md.SetAliRootVersion("28.02.2006");
  md.SetComment("Test");
  
  printf("\n\n\n fCalibData\n");
  //fCalibData->Print(0);
  //printf("\n\n\n fCalibData\n");
  
  Bool_t result = man->Put(fCalibData,id,&md);

  return result;
}

