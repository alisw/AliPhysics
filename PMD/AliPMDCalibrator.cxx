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
#include "TMath.h"

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

//#include "AliPMDCleanNoise.h"
//#include "AliPMDCleaner.h"

#include "AliPMDPedestal.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

ClassImp(AliPMDCalibrator)


AliPMDCalibrator::AliPMDCalibrator():
  fCalibGain(new AliPMDCalibData()), fCalibPed(new AliPMDPedestal())
{
  // Standard Constructor
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      fHdetIso[idet] = NULL ;
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  fHsmIso[idet][ismn] = NULL ;
	  for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
	      for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		  fGainFact[idet][ismn][jrow][kcol] = 0.0;
		  fHadcIso[idet][ismn][jrow][kcol]  = NULL;
		}
	    }
	}
    }
}
// ------------------------------------------------------------------------ //
AliPMDCalibrator::AliPMDCalibrator(const AliPMDCalibrator &pmdcalibrator):
  fCalibGain(new AliPMDCalibData()), fCalibPed(new AliPMDPedestal())
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      fHdetIso[idet] = pmdcalibrator.fHdetIso[idet] ;
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  fHsmIso[idet][ismn] = pmdcalibrator.fHsmIso[idet][ismn] ;
	  for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
	      for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		  fGainFact[idet][ismn][jrow][kcol] = pmdcalibrator.fGainFact[idet][ismn][jrow][kcol];
		  fHadcIso[idet][ismn][jrow][kcol]  = pmdcalibrator.fHadcIso[idet][ismn][jrow][kcol];
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
      for(Int_t idet = 0; idet < kDet; idet++)
	{
	  fHdetIso[idet] = pmdcalibrator.fHdetIso[idet] ;
	  for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	    {
	      fHsmIso[idet][ismn] = pmdcalibrator.fHsmIso[idet][ismn] ;
	      for(Int_t jrow = 0; jrow < kMaxRow;jrow++)
		{
		  for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		    {
		      fGainFact[idet][ismn][jrow][kcol] =
			pmdcalibrator.fGainFact[idet][ismn][jrow][kcol];
		      fHadcIso[idet][ismn][jrow][kcol]  =
			pmdcalibrator.fHadcIso[idet][ismn][jrow][kcol];
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
  // destructor
  if(fHdetIso) delete fHdetIso ;
  if(fHsmIso)  delete fHsmIso ;
  if(fHadcIso) delete fHadcIso ;
  delete fCalibGain;
  delete fCalibPed;
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
  char hname2[kDet];
  char htitle2[2];
  char hname[kMaxSMN];
  char hname24[kMaxSMN];
  char hnameiso[120];
  char htitle1[120];
  
  for(Int_t d = 0; d < kDet; d++) {
    sprintf(hname2,"Isolated cell adc for Det Plane %d",d);
    fHdetIso[d]= new TH1F(hname2,htitle2,100,0,1000);
    for(Int_t i1 = 0; i1 < kMaxSMN; i1++) {
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
	  fHadcIso[d][i1][j1][k1] = new TH1F(hnameiso,htitle1,100,0.,1000.);
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

  const Int_t kCellNeighbour = 6;

  Int_t neibx[6] = {1,0,-1,-1,0,1};
  Int_t neiby[6] = {0,1,1,0,-1,-1};

  Int_t id1,jd1;            //neighbour row/col
  Int_t countisocell = 0 ;  //number of isilated cell
  Int_t isocount;           //number of neighbours with 0 signal
  Int_t d1[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Int_t maxhit;
  Int_t nhit[kDet][kMaxSMN];
  Int_t nhitcell[kDet][kMaxSMN][kMaxRow][kMaxCol];
  
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
        {
	  nhit[idet][ismn] = 0;
	  for(Int_t irow = 0; irow < kMaxRow; irow++)
	    {
	      for(Int_t icol = 0; icol < kMaxCol; icol++)
		{
		  d1[idet][ismn][irow][icol] = 0;
		  nhitcell[idet][ismn][irow][icol] = 0;
		}
	    }
	}
    }

  Float_t tempnhit1[kDet][kMaxSMN];
  Float_t tempnhit2[kDet][kMaxSMN];
  Float_t meannhit[kDet][kMaxSMN];
  Float_t meanSqnhit[kDet][kMaxSMN];
  Float_t sigmanhit[kDet][kMaxSMN];
  Float_t count[kDet][kMaxSMN];
  Float_t nhitcut[kDet][kMaxSMN];

  for (Int_t idet = 0; idet < kDet; idet++)
  {
      for (Int_t ismn = 0; ismn < kMaxSMN; ismn++)
      {
	  tempnhit1[idet][ismn]  = 0.;
	  tempnhit2[idet][ismn]  = 0.;
	  meannhit[idet][ismn]   = 0.;
	  meanSqnhit[idet][ismn] = 0.;
	  sigmanhit[idet][ismn]  = 0.;
	  count[idet][ismn]      = 0.;
	  nhitcut[idet][ismn]    = 0.;
      }
  }


  //accessing raw data
  AliRawReaderFile reader(".");
  AliPMDRawStream stream(&reader);
  while(reader.NextEvent())
    { 
      // New PMD Reader is plugged in
	Int_t iddl = -1;
	while ((iddl = stream.DdlData(&pmdddlcont)) >=0) {

	    //reader.Select("PMD", iddl, iddl);
	    //stream.DdlData(iddl,&pmdddlcont);
	  Int_t ientries = pmdddlcont.GetEntries();
	  for (Int_t ient = 0; ient < ientries; ient++)
	    {
	      AliPMDddldata *pmdddl = 
		(AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	      Int_t idet = pmdddl->GetDetector();
	      Int_t ismn = pmdddl->GetSMN();
	      Int_t irow = pmdddl->GetRow();
	      Int_t icol = pmdddl->GetColumn();
	      Float_t isig1 = pmdddl->GetSignal();
	      // Pedestal Subtraction
	      Int_t   pedmeanrms = 
		fCalibPed->GetPedMeanRms(idet,ismn,irow,icol);
	      Int_t   pedrms1    = (Int_t) pedmeanrms%100;
	      Float_t pedrms     = (Float_t)pedrms1/10.;
	      Float_t pedmean    = (Float_t) (pedmeanrms - pedrms1)/1000.0;
	      Float_t isig = isig1 - (pedmean + 3.0*pedrms);
	      if (isig>0)
		{
		  d1[idet][ismn][irow][icol] = (Int_t)isig;
		  nhitcell[idet][ismn][irow][icol] += 1; 
		}
	    }//ient loop
	  pmdddlcont.Clear();
	}//iddl loop
      maxhit = 0;
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
				      countisocell++;
				      maxhit++;
                                      fHdetIso[idet]->
					Fill(d1[idet][ismn][irow][icol]);
				      fHsmIso[idet][ismn]->
					Fill(d1[idet][ismn][irow][icol]);
				      fHadcIso[idet][ismn][irow][icol]->
					Fill(d1[idet][ismn][irow][icol]);
				    }
				}
			    }  // neigh cell cond.
			}     // d>0 cond.
		    }
		}
	    }
	} //det
    }//event loop

  // Mean and Sigma Calculations
  for(Int_t idet=0; idet < kDet; idet++){
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++){
	  for(Int_t irow = 0; irow < kMaxRow; irow++){
	      for(Int_t icol = 0; icol < kMaxCol; icol++){
		  if(nhitcell[idet][ismn][irow][icol]>0){
		      count[idet][ismn] += 1;
		      tempnhit1[idet][ismn] += nhitcell[idet][ismn][irow][icol];
		      tempnhit2[idet][ismn] += nhitcell[idet][ismn][irow][icol]
			  *nhitcell[idet][ismn][irow][icol];
		  }
	      }
	  }
      }
  }//det loop

  //cout<<"nhit cell = "<<idet<<"  "<<ismn<<"  "
  //  <<irow<<"  "<<icol<<"  "<<nhitcell[idet][ismn][irow][icol]<<endl;
  //count[idet][ismn] += 1;

  for(Int_t i=0; i < kDet; i++)
    {
      for(Int_t j=0; j < kMaxSMN; j++)
	{
	  if(count[i][j] > 0.)
	    {
		meannhit[i][j]   = tempnhit1[i][j]/count[i][j];
		meanSqnhit[i][j] = tempnhit2[i][j]/count[i][j]; 
		sigmanhit[i][j]  = sqrt(meanSqnhit[i][j]-
				      (meannhit[i][j]*meannhit[i][j]));
		nhitcut[i][j]    = 3*sigmanhit[i][j] + meannhit[i][j];
	    }
	}
    }


  Double_t histdetMean[kDet];
  Double_t histMean[kDet][kMaxSMN];
  Double_t isoMean[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Double_t smNormFactor[kDet][kMaxSMN];
  
  for(Int_t det1 = 0; det1 < kDet; det1++)
    {
      histdetMean[det1]= fHdetIso[det1]->GetMean();
      for(Int_t i1 = 0; i1 < kMaxSMN; i1++)
	{
	  histMean[det1][i1]= fHsmIso[det1][i1]->GetMean();
	  if(histMean[det1][i1]>0.0 && histdetMean[det1]>0.0)
	    {
	      smNormFactor[det1][i1]= histdetMean[det1]/histMean[det1][i1];
	    }
	  for(Int_t j1 = 0; j1 < kMaxRow; j1++)
	    {
	      for(Int_t k1 = 0; k1 < kMaxCol; k1++)
		{
		 if(nhitcell[det1][i1][j1][k1]> nhitcut[det1][i1]) fGainFact[det1][i1][j1][k1]=-99.0;
		   if(nhitcell[det1][i1][j1][k1]< nhitcut[det1][i1])
		     {
		      isoMean[det1][i1][j1][k1]=fHadcIso[det1][i1][j1][k1]->
			GetMean();
		      if(isoMean[det1][i1][j1][k1]>0.0 && histMean[det1][i1]>0.0)
			{
			  fGainFact[det1][i1][j1][k1]=
			    isoMean[det1][i1][j1][k1]/(histMean[det1][i1]*
						     smNormFactor[det1][i1]);

			}
		    }   
 		   Float_t gain = fGainFact[det1][i1][j1][k1];
 		   fCalibGain->SetGainFact(det1,i1,j1,k1,gain);
		}
                          
 			
	    }
	}
    }
}//CalculateIsoCell()

// ------------------------------------------------------------------------ //

Bool_t AliPMDCalibrator::Store()
{
  AliCDBManager *man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if(!man->IsDefaultStorageSet()) return kFALSE;
  AliCDBId id("PMD/Calib/Gain",0,999999999);
  AliCDBMetaData md;
  md.SetBeamPeriod(0);
  md.SetComment("Test");
  
  printf("\n\n\n fCalibData\n");
  //fCalibData->Print(0);
  //printf("\n\n\n fCalibData\n");
  
  Bool_t result = man->Put(fCalibGain,id,&md);

  return result;
}

