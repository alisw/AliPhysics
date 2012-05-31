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

#include "Riostream.h"
#include "TF1.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TMath.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"
#include "AliPMDCalibGain.h"

ClassImp(AliPMDCalibGain)

AliPMDCalibGain::AliPMDCalibGain():
  TObject(),
  fpw(NULL)
{
  // Standard Constructor
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  fSMIso[idet][ismn]   = 0.;
	  fSMCount[idet][ismn] = 0.;
	  fCountSm[idet][ismn]=0.;
	  fTempnhit[idet][ismn]=0.;
	  fTempnhitSq[idet][ismn]=0.;
	  for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
	      for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		  fCellIso[idet][ismn][jrow][kcol]    = 0.;
		  fCellCount[idet][ismn][jrow][kcol]  = 0.;
		  fNhitCell[idet][ismn][jrow][kcol]   = 0.;
		  fPedMeanRMS[idet][ismn][jrow][kcol] = 0.;
		  fHotFlag[idet][ismn][jrow][kcol]    = 0.;
		  
		}
	    }
	}
    }
  
}
// ------------------------------------------------------------------------ //
AliPMDCalibGain::AliPMDCalibGain(const AliPMDCalibGain &pmdcalibgain):
  TObject(pmdcalibgain),
  fpw(NULL)
{
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  fSMIso[idet][ismn] = pmdcalibgain.fSMIso[idet][ismn] ;
	  fSMCount[idet][ismn] = pmdcalibgain.fSMCount[idet][ismn] ;
	  fCountSm[idet][ismn] = pmdcalibgain.fCountSm[idet][ismn];
	  fTempnhit[idet][ismn] = pmdcalibgain.fTempnhit[idet][ismn];
	  fTempnhitSq[idet][ismn] = pmdcalibgain.fTempnhitSq[idet][ismn];
	  for(Int_t jrow = 0; jrow < kMaxRow; jrow++)
	    {
	      for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		{
		  fCellIso[idet][ismn][jrow][kcol] = pmdcalibgain.fCellIso[idet][ismn][jrow][kcol];
		  fCellCount[idet][ismn][jrow][kcol] = pmdcalibgain.fCellCount[idet][ismn][jrow][kcol];
		  fNhitCell[idet][ismn][jrow][kcol] = pmdcalibgain.fNhitCell[idet][ismn][jrow][kcol];
		  fPedMeanRMS[idet][ismn][jrow][kcol] = pmdcalibgain.fPedMeanRMS[idet][ismn][jrow][kcol];
		  fHotFlag[idet][ismn][jrow][kcol] = pmdcalibgain.fHotFlag[idet][ismn][jrow][kcol];
		  
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
      this->fpw = pmdcalibgain.fpw;
      for(Int_t idet = 0; idet < kDet; idet++)
	{
	  for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	    {
	      fSMIso[idet][ismn]      = pmdcalibgain.fSMIso[idet][ismn];
	      fSMCount[idet][ismn]    = pmdcalibgain.fSMCount[idet][ismn];
	      fCountSm[idet][ismn]    = pmdcalibgain.fCountSm[idet][ismn];
	      fTempnhit[idet][ismn]   = pmdcalibgain.fTempnhit[idet][ismn];
	      fTempnhitSq[idet][ismn] = pmdcalibgain.fTempnhitSq[idet][ismn];
	      for(Int_t jrow = 0; jrow < kMaxRow;jrow++)
		{
		  for(Int_t kcol = 0; kcol < kMaxCol; kcol++)
		    {
		      fCellIso[idet][ismn][jrow][kcol] = pmdcalibgain.fCellIso[idet][ismn][jrow][kcol];
		      fCellCount[idet][ismn][jrow][kcol] = pmdcalibgain.fCellCount[idet][ismn][jrow][kcol];
		      fNhitCell[idet][ismn][jrow][kcol] = pmdcalibgain.fNhitCell[idet][ismn][jrow][kcol]; //za 
		      fPedMeanRMS[idet][ismn][jrow][kcol] = pmdcalibgain.fPedMeanRMS[idet][ismn][jrow][kcol];
		      fHotFlag[idet][ismn][jrow][kcol] = pmdcalibgain.fHotFlag[idet][ismn][jrow][kcol];
		      
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
}

// ------------------------------------------------------------------------ //

Int_t AliPMDCalibGain::ExtractPedestal(const Char_t *rootFile)
{
  // Pedestal extraction from the PMD_PED.root file
  // To be called once at the beginning
  
  Int_t   det=0, sm=0, row=0, col=0;
  Float_t mean=0., rms=0.;
  
  TFile *pedfile = new TFile(rootFile);
  
  if(!pedfile)
    {
      printf("ERROR --- NO PEDESTAL (PMD_PED1.root) FILE IS FOUND --- STOP GAIN DA\n");
      return -3;
    }


  TTree *ped =(TTree*)pedfile->Get("ped");

  ped->SetBranchAddress("det",&det);
  ped->SetBranchAddress("sm",&sm);
  ped->SetBranchAddress("row",&row);
  ped->SetBranchAddress("col",&col);
  ped->SetBranchAddress("mean",&mean);
  ped->SetBranchAddress("rms",&rms);

  Int_t nentries = (Int_t)ped->GetEntries();

  for (Int_t ient = 0; ient < nentries; ient++)
    {
      ped->GetEntry(ient);
      fPedMeanRMS[det][sm][row][col] = mean + 3.*rms;
      //printf("Mean= %f, RMS= %f, PedMeanRMS=%f\n",mean,rms,fPedMeanRMS[det][sm][row][col]);

    }

  pedfile->Close();
  delete pedfile;
  pedfile = 0x0;

  return 1;
}
//------------------------------------------------------------------------------------------------

Int_t AliPMDCalibGain::ExtractHotChannel(const Char_t *rootFile)
{
  // HotChannel extraction from the PMD_HOT.root file
  // To be called once at the beginning

  Int_t   det=0, sm=0, row=0, col=0;
  Float_t flag=0.;

  TFile *hotmapfile = new TFile(rootFile);

  if(!hotmapfile)
    {
      printf(" NO HOTCHANNEL MAP (PMD_HOT.root) FILE IS FOUND \n");

      for (Int_t idet = 0; idet < kDet; idet++)
	{
	  for (Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	    {
	      for (Int_t irow = 0; irow < kMaxRow; irow++)
		{
		  for (Int_t icol = 0; icol < kMaxCol; icol++)
		    {
		      fHotFlag[idet][ismn][irow][icol] = 0.;
		    }
		}
	    }
	}
    }


  TTree *hot =(TTree*)hotmapfile->Get("hot");
  
  hot->SetBranchAddress("det",&det);
  hot->SetBranchAddress("sm",&sm);
  hot->SetBranchAddress("row",&row);
  hot->SetBranchAddress("col",&col);
  hot->SetBranchAddress("flag",&flag);
  
  Int_t nentries = (Int_t)hot->GetEntries();
  
  for (Int_t ient = 0; ient < nentries; ient++)
    {
      hot->GetEntry(ient);
      fHotFlag[det][sm][row][col] = flag;

      //printf(" HotFlag=%f\n",fHotFlag[det][sm][row][col]);
    }

  hotmapfile->Close();
  delete hotmapfile;
  hotmapfile = 0x0;
  
  return 1;
}


// ------------------------------------------------------------------------ //

void AliPMDCalibGain::ReadTempFile(const Char_t *tempFile)
{
  // Read the variables from the file
  
  ifstream intmpfile;
  intmpfile.open(tempFile);

  Int_t iddet = 0, issm = 0, irrow = 0, iccol = 0;
  Float_t smcount = 0., smiso = 0.;
  Float_t cellcount = 0., celliso = 0.;


  for (Int_t idet = 0; idet < kDet; idet++)
    {
      for (Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  intmpfile >> iddet >> issm >> smcount >> smiso;
	  fSMCount[idet][ism] = smcount;
	  fSMIso[idet][ism]   = smiso;
	}
    }

  for (Int_t idet = 0; idet < kDet; idet++)
    {
      for (Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  for (Int_t irow = 0; irow < kMaxRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kMaxCol; icol++)
		{
		  intmpfile >> iddet >> issm >> irrow >> iccol
			    >> cellcount >> celliso;
		  fCellCount[idet][ism][irow][icol] = cellcount;
		  fCellIso[idet][ism][irow][icol]   = celliso;
		}
	    }
	}
    }
  

  intmpfile.close();

}
// ------------------------------------------------------------------------ //
void AliPMDCalibGain::WriteTempFile(const Char_t *tempFile)
{
  // Write the Temporary file if the required statics is not achieved


  /*
    Following variables to be written to a file
    fDetIso[idet] ;
    fSMIso[idet][ismn]; 
    fCellIso[idet][ismn][irow][icol]; 
    
    fDetCount[idet];
    fSMCount[idet][ismn];
    fCellCount[idet][ismn][irow][icol];
  */				  


  fpw = fopen(tempFile,"w+");
  
  for (Int_t idet = 0; idet < kDet; idet++)
    {
     // fprintf(fpw,"%d %f %f\n",idet,fDetCount[idet],fDetIso[idet]);
    }
  
  for (Int_t idet = 0; idet < kDet; idet++)
    {
      for (Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  fprintf(fpw,"%d %d %f %f\n",idet,ism, fSMCount[idet][ism],fSMIso[idet][ism]);
	}
    }
  
  for (Int_t idet = 0; idet < kDet; idet++)
    {
      for (Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  for (Int_t irow = 0; irow < kMaxRow; irow++)
	    {
	      for (Int_t icol = 0; icol < kMaxCol; icol++)
		{
		  fprintf(fpw,"%d %d %d %d %f %f\n",idet,ism,irow,icol,
			  fCellCount[idet][ism][irow][icol],
			  fCellIso[idet][ism][irow][icol]);
		}
	    }
	}
    }
  
  fclose(fpw);
  
}

// ------------------------------------------------------------------------ //

Bool_t AliPMDCalibGain::ProcessEvent(AliRawReader *rawReader, TObjArray *pmdddlcont)
{
  // Calculates the ADC of isolated cell
  
  const Int_t kDDL           = AliDAQ::NumberOfDdls("PMD");
  const Int_t kCellNeighbour = 6;
  Int_t neibx[6] = {1,0,-1,-1,0,1};
  Int_t neiby[6] = {0,1,1,0,-1,-1};
  
  Int_t id1 = 0,jd1 = 0;  //neighbour row/col
  Int_t isocount = 0;     //number of neighbours with 0 signal

  Float_t d1[kDet][kMaxSMN][kMaxRow][kMaxCol];
  
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
        {
	  for(Int_t irow = 0; irow < kMaxRow; irow++)
            {
              for(Int_t icol = 0; icol < kMaxCol; icol++)
                {
		  d1[idet][ismn][irow][icol] = 0.;
                }
	    }
	}
    }
  
  AliPMDRawStream rawStream(rawReader);

  Int_t iddl = -1;
  
  Int_t numberofDDLs = 0;

  while ((iddl = rawStream.DdlData(pmdddlcont)) >=0) {
    numberofDDLs++;
    
    Int_t ientries = pmdddlcont->GetEntries();
    
    for (Int_t ient = 0; ient < ientries; ient++)
      {
	AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont->UncheckedAt(ient);
	
	Int_t idet = pmdddl->GetDetector();
	Int_t ismn = pmdddl->GetSMN();
	Int_t mcm = pmdddl->GetMCM();
	//Int_t ichno = pmdddl->GetChannel();
	Int_t irow = pmdddl->GetRow();
	Int_t icol = pmdddl->GetColumn();
	Int_t isig = pmdddl->GetSignal();
	
	// This is the protection not to crash the code 
	
	if(mcm == 0) continue;
	if (irow < 0 || icol < 0 || irow > 47 || icol > 95) continue;
	
	// Pedestal subtraction
	if(fHotFlag[idet][ismn][irow][icol] == 1.0) isig = 0;
	
	if (isig>0)
	  {
	    d1[idet][ismn][irow][icol] =
	      (Float_t) isig - fPedMeanRMS[idet][ismn][irow][icol];
	    
	    //printf("Signal_ped_subtracted=%f, pedestal=%f\n",d1[idet][ismn][irow][icol]),fPedMeanRMS[idet][ismn][irow][icol];
	    
	    fNhitCell[idet][ismn][irow][icol]++;     //cell hit frequency
	  }
	

      }
    pmdddlcont->Delete();
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
			  if (id1 < 0) id1 = 0;
			  if (id1 > kMaxRow-1) id1 = kMaxRow - 1;
			  if (jd1 < 0) jd1 = 0;
			  if (jd1 > kMaxCol-1) jd1 = kMaxCol - 1;
			  if(d1[idet][ismn][id1][jd1] == 0)
			    {
			      isocount++;
			      if(isocount == kCellNeighbour)
				{
				  //fDetIso[idet] += d1[idet][ismn][irow][icol];
				  fSMIso[idet][ismn] += d1[idet][ismn][irow][icol];
				  fCellIso[idet][ismn][irow][icol] += d1[idet][ismn][irow][icol];
				  //fDetCount[idet]++;
				  fSMCount[idet][ismn]++;
				  fCellCount[idet][ismn][irow][icol]++;
				  
				}
			    }
			}  // neigh cell cond.
		    }     // d>0 cond.
		}
	    }
	}
    }
  
  for(Int_t idet=0; idet < kDet; idet++)
    {
      for(Int_t ismn = 0; ismn < kMaxSMN; ismn++)
	{
	  for(Int_t irow = 0; irow < kMaxRow; irow++)
	    {
	      for(Int_t icol = 0; icol < kMaxCol; icol++)
		{
		  if(fNhitCell[idet][ismn][irow][icol]>0)
		    {
		      fCountSm[idet][ismn]  += 1;
		      fTempnhit[idet][ismn] += fNhitCell[idet][ismn][irow][icol];
		      fTempnhitSq[idet][ismn] += fNhitCell[idet][ismn][irow][icol]
			*fNhitCell[idet][ismn][irow][icol];
		    }
		}
	    }
	}
    }
  
  if (numberofDDLs < kDDL) return kFALSE;
  return kTRUE;
  
}
// ------------------------------------------------------------------------ //
void AliPMDCalibGain::Analyse(TTree *gaintree, TTree *meantree)
{
  // Calculates the mean
  Int_t   det = 0, sm = 0, row = 0, col = 0;
  Float_t gain = 0.;
  Float_t cellmean = 0.;

  Float_t modmean[2][24];

  for (Int_t idet=0; idet < 2; idet++)
    {
      for (Int_t ism = 0; ism < 24; ism++)
	{
	  modmean[idet][ism] = 0.;
	}
    }

  gaintree->Branch("det",&det,"det/I");
  gaintree->Branch("sm",&sm,"sm/I");
  gaintree->Branch("row",&row,"row/I");
  gaintree->Branch("col",&col,"col/I");
  gaintree->Branch("gain",&gain,"gain/F");
  
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  if (fSMCount[idet][ism] > 0)
	    modmean[idet][ism] = fSMIso[idet][ism]/fSMCount[idet][ism];
	  for(Int_t irow = 0; irow < kMaxRow; irow++)
	    {
	      for(Int_t icol = 0; icol < kMaxCol; icol++)
		{
		  if (fCellCount[idet][ism][irow][icol] > 0.)
		    {
		      cellmean = fCellIso[idet][ism][irow][icol]/fCellCount[idet][ism][irow][icol];
		    }
		  det      = idet;
		  sm       = ism;
		  row      = irow;
		  col      = icol;
		  if (cellmean > 0.0 && fCellCount[idet][ism][irow][icol]>0.)
		    {
		      gain = cellmean/modmean[idet][ism];
		    }
		  else
		    {
		      gain = 0.;
		    }
		  //if(fCellCount[idet][ism][irow][icol]>0.) printf("CellCount =%f, gain= %f\n",fCellCount[idet][ism][irow][icol],gain);
		  gaintree->Fill();
		}
	    }
	}
    }

  Float_t smmean;

  // Writing each module mean value
  meantree->Branch("det",&det,"det/I");
  meantree->Branch("sm",&sm,"sm/I");
  meantree->Branch("smmean",&smmean,"row/F");
  
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for (Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  det    = idet;
	  sm     = ism;
	  smmean = modmean[idet][ism];
	  meantree->Fill();
	}
    }

}
// ------------------------------------------------------------------------ //
void AliPMDCalibGain::FindHotCell(TTree *hottree, Float_t xvar)
{
  // Calculates the mean
  Int_t   det = 0, sm = 0, row = 0, col = 0;
  Float_t flag = 0.;
  Float_t meannhit = 0.;
  Float_t meanSqnhit = 0.;
  Float_t sigmanhit = 0.,nhitcut = 0.;

  //Float_t xvar = 5.;

  hottree->Branch("det",&det,"det/I");
  hottree->Branch("sm",&sm,"sm/I");
  hottree->Branch("row",&row,"row/I");
  hottree->Branch("col",&col,"col/I");
  hottree->Branch("flag",&flag,"flag/F");
  
  for(Int_t idet = 0; idet < kDet; idet++)
    {
      for(Int_t ism = 0; ism < kMaxSMN; ism++)
	{
	  if (fCountSm[idet][ism]> 0)
	    {
	      meannhit   = fTempnhit[idet][ism]/fCountSm[idet][ism];
	      meanSqnhit = fTempnhitSq[idet][ism]/fCountSm[idet][ism];
	      sigmanhit  = sqrt(meanSqnhit-(meannhit*meannhit));
	      nhitcut    = meannhit + xvar*sigmanhit;

	      for(Int_t irow = 0; irow < kMaxRow; irow++)
		{
		  for(Int_t icol = 0; icol < kMaxCol; icol++)
		    {
		      det      = idet;
		      sm       = ism;
		      row      = irow;
		      col      = icol;
		      
		      if(fNhitCell[idet][ism][irow][icol] > nhitcut)
			{
			  flag  = 1.0;
			}
		      else
			{
			  flag = 0.;
			}
		      hottree->Fill();
		    }
		  
		}
	    }
	}
    }
}

