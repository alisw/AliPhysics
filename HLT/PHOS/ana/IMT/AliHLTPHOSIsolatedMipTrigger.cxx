/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSIsolatedMipTrigger.h"
#include "stdio.h"
#include "AliHLTPHOSDigit.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include <iostream>
#include "AliHLTPHOSMipCluster.h"
#include "AliHLTPHOSMipClusterManager.h"


#define Q_MIN 27
#define Q_MAX 33

#define XBIN_LOW  0
#define XBIN_UP   1023
#define N_BINS    1023

AliHLTPHOSIsolatedMipTrigger::AliHLTPHOSIsolatedMipTrigger() :  fMipLowCut(3),  
								fMipHighCut(150),
								fModuleId(2)
{
  fClusterManager = new AliHLTPHOSMipClusterManager();
  Init();
}


AliHLTPHOSIsolatedMipTrigger::~AliHLTPHOSIsolatedMipTrigger()
{

}

/*
void
AliHLTPHOSIsolatedMipTrigger::SetFilePath(char *path)
{
  sprintf(fCurrentFilePath, "%s", path);
  printf("\n AliHLTPHOSIsolatedMipTrigger::SetFilePath, filepath was set to %s \n", fCurrentFilePath);
}
*/

void
AliHLTPHOSIsolatedMipTrigger::SetMipLowCut(Float_t lCut)
{
  fMipLowCut = lCut;
}

void
AliHLTPHOSIsolatedMipTrigger::SetMipHighCut(Float_t hCut)
{
  fMipHighCut = hCut;
}


int
AliHLTPHOSIsolatedMipTrigger::ScanFileList()
{
  int iRet = 0;
  FILE *fp = fopen("filelist.txt", "r");
  
  if(fp != 0)
    {
      fscanf(fp, "%d\n", &fNFiles);
      cout << "there are " << fNFiles << "  files " <<endl;
 
      for(int i=0; i< fNFiles; i++)
	{
	  fscanf(fp, "%s\n"  ,fFileList[i]);
	}
      
      for(int i=0; i< fNFiles; i++)
	{
	  cout << fFileList[i] << endl;
	}
    }
  else
    {
      cout << "ERROR opening file, filelist.txt" << endl;
      iRet = -1;
    } 
  
  return iRet;
}



int
AliHLTPHOSIsolatedMipTrigger::Analyze()
{
  //  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDigit" , 500); 

  ScanFileList();
  cout << "AliHLTPHOSIsolatedMipTrigger::Analyze fNFiles = " << fNFiles  <<endl;

  for(int file=0; file < fNFiles; file ++)
    {
 
     printf("\n AliHLTPHOSIsolatedMipTrigger::Analyze, analyzing file %s\n", fFileList[file]);   
      
      TFile *infile = infile->Open(fFileList[file],"read");
      //      TFile *infile = infile->Open("/tmp/data2/digits151007/run8340_digitTree_5.root","read");

      AliHLTPHOSDigit *tmpDigit = 0;  

 
      if(infile == 0)
	{
	  cout << "ERROR could not open file" << endl;
	}
      else
	{
	  
	  infile->cd();
	  cout << "File successfully opened, creating" << endl;
	  TChain *cain= new TChain("digitTree");
	  TClonesArray *digArray = new TClonesArray("AliHLTPHOSDigit" , 100);

	  cain->Add(fFileList[file]);

	  int nThrees =  cain->GetEntries() ;
	  cout << "there are "  <<  nThrees  <<  " entries in this file" << endl;
	  cain->SetBranchAddress("Digit", &digArray);

	  for(int i=0; i< nThrees; i++)
	    {
	      fClusterManager->ResetMipManager(); 
	      cain->GetEntry(i);
	      int nArrays = digArray->GetEntriesFast();

	      if(i%100 == 0)
		{
		  cout << "processing event " << i << "  of file  " << fFileList[file] <<endl;
		}
	      
	 
	      for(int j=0; j < nArrays; j++)
		{

		  tmpDigit = (AliHLTPHOSDigit*)digArray->At(j);
		      
		  if( IsMipCandidate(tmpDigit) == true)
		    {
		      fClusterManager->AddDigit(tmpDigit);
		    }
		}
	      
	      FillClusterHistograms();
	    }
 
	  delete digArray;
	  delete cain;
	  infile->Close(); 
	  delete infile;
	}
      
      WriteHistograms();

    }
 
  return 0;
}


bool
AliHLTPHOSIsolatedMipTrigger::IsMipCandidate(AliHLTPHOSDigit  *digit)
{
  bool ret = false;
  Float_t amplitude = digit->GetAmplitude();
  Int_t *data =  digit->GetRawData();
  Float_t sum = 0;
  Float_t q;

  if( amplitude > fMipLowCut && amplitude < fMipHighCut)
    {
      int start = digit->GetPreSamples();
      int end = digit->GetTotalSamples() - start;
 
      for(int  i = start; i < end; i++)
	{
	  sum+=data[i];
	}
      
      q = ((Float_t)sum)/amplitude;

      if( (q > Q_MIN) && (q < Q_MAX))
	{
	  ret = true;
	}
    }
  return ret;

}



void
AliHLTPHOSIsolatedMipTrigger::Init()
{
  
  char tmpHistoName[256];

  for(int z = 0; z < N_ZROWS_MOD; z ++)
    {
      for(int x = 0; x < N_XCOLUMNS_MOD; x ++)
	{
	  sprintf(tmpHistoName, "ClusterEnergies3x3_%d_%d_%d_%d",(int)fModuleId,  z, x, 1);
	  fClusterHistograms[z][x] = new TH1F( tmpHistoName, tmpHistoName, N_BINS, XBIN_LOW, XBIN_UP);
	}
    } 
}



void 
AliHLTPHOSIsolatedMipTrigger::FillClusterHistograms()
{
  int cnt = 0;
  AliHLTPHOSMipCluster* tmp;

  for(int z = 0; z < N_ZROWS_MOD; z ++)
    {
      for(int x = 0; x < N_XCOLUMNS_MOD; x ++)
	{
	  tmp = fClusterManager->GetCluster(z, x);
	  if(tmp != 0)
	    {
	      fClusterHistograms[z][x]->Fill(tmp->Get3x3Sum());
	      cnt ++;
	    }
	}
    } 
  
  //  cout << "filled  " << cnt <<"  clusters"  <<endl;
}

void
AliHLTPHOSIsolatedMipTrigger::WriteHistograms()
{
  char tmpFileName[256];
  sprintf(tmpFileName, "ClusterEnergies.root");

  TFile *clusterFile =  new TFile(tmpFileName, "recreate");

  for(int z = 0; z < N_ZROWS_MOD; z ++)
    {
      for(int x = 0; x < N_XCOLUMNS_MOD; x ++)
	{
	  if(fClusterHistograms[z][x]->GetEntries() > 10 )
	    {
	      fClusterHistograms[z][x]->Write();
	    }
	}
    }  

  clusterFile->Close();

  delete clusterFile;
}
