// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Boris Polichtchouk & Per Thomas Hille for the ALICE           *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRcuHistogramProducer.h"
#include <iostream>
//#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "TFile.h"
#include "unistd.h"
#include <time.h>
#include "AliHLTCaloUtilities.h" 

#define THRESHOLD 30

using  namespace std;


/*************************************************************************
* Class AliHLTPHOSRcuHistogramProducer accumulating histograms           *
* with amplitudes per PHOS channel                                       *
* It is intended to run at the HLT farm                                  *
* and it fills the histograms with amplitudes per channel.               * 
* Usage example see in PHOS/macros/Shuttle/AliPHOSCalibHistoProducer.C   *
**************************************************************************/


AliHLTPHOSRcuHistogramProducer:: AliHLTPHOSRcuHistogramProducer(): // AliHLTPHOSBase(), 
  //		    AliHLTPHOSRcuProcessor(),
								    fCellAccEnergy(), 
								    fModuleID(0), 
								    fRcuX(0), 
								    fRcuZ(0),
								    fUtilitiesPtr(0)
{
  fUtilitiesPtr = new  AliHLTPHOSUtilities;

  //Default constructor
  //  cout << "WARNING: You cannot invoke the AliHLTPHOSRcuHistogramProducer without arguments" << endl;
  //  cout << "Usage AliHLTPHOSRcuHistogramProducer(ModuleID, X. Z)" << endl;
} 


AliHLTPHOSRcuHistogramProducer::AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ): //AliHLTPHOSBase(), 
  //	AliHLTPHOSRcuProcessor(),
																fCellAccEnergy(),
																fModuleID(moduleID), 
																fRcuX(rcuX), 
																fRcuZ(rcuZ),
																fUtilitiesPtr(0)
{

  //Se header file for documentation
  fUtilitiesPtr = new AliHLTPHOSUtilities();

  char *tmp = getenv("HOME");
  if(tmp == 0)
    {
      //      cout << "ERROR, environment vriable HOME is not set" << endl;
    }
  else
    {
      //  sprintf(fHistoOutDir,"%s/rundir/output/histograms",tmp);
      sprintf(fHistoOutDir,"%s/rundir/output/histograms/",tmp);
    }
  Init();
}


AliHLTPHOSRcuHistogramProducer::~AliHLTPHOSRcuHistogramProducer()
{
  //Destructor
}


void
AliHLTPHOSRcuHistogramProducer::SetHistoOutDir(char *outDir)
{
  //comment
  sprintf(fHistoOutDir,"%s", outDir);
}

void 
AliHLTPHOSRcuHistogramProducer::SetDefaultHistoOutDir()
{
  //comment
  char *tmp = getenv("HOME/rundir"); 
  sprintf(fHistoOutDir,"%s/rundir/output/histograms", tmp);
  //testing wether or not directry exist
  FILE *fp = fopen(fHistoOutDir, "w");
  
  if(fp == 0)
    {
//       cout << "ERROR, directory =" << fHistoOutDir << "  Doesnt exist, or you don have write permissions to the directory" << endl;
//       cout << "WARNING: Histograms will not bew written to files at end of run unless a valid directory is set" << endl;
//       cout << "INFO, You must either" << endl;
//       cout << "1) Create the directory " << fHistoOutDir <<  "Manually" <<endl;
//       cout << "OR "<< endl;
//       cout << "2)  Se a valid output directory with the function AliHLTPHOSRcuHistogramProducer::SetHistoOutDir(*char outdir) "<< endl;
   }
  else
    {
      //      cout << "INFO: Output ddirectory for Histograms was set tot" << fHistoOutDir << endl;
      //      cout << "INFO: if you want another output directory use the AliHLTPHOSRcuHistogramProducer::SetHistoOutDir(*char outdir)" << endl;
    }

} 

void
AliHLTPHOSRcuHistogramProducer::Init()
{
  //See header file for documentation
  char tmpHistoName[256];
  int geomx;
  int geomz;

  for(int gain=0; gain< NGAINS; gain++)
    {
      sprintf(tmpHistoName, "DeadChanneMap_Module%d_rcuz%d_rcux%d_gain%d",(int)fModuleID,  fRcuZ, fRcuX, gain);
      //    fDeadChannelMapHistogramPtrs[gain] = new TH2D( tmpHistoName, tmpHistoName, NBINS, XBINLOW, XBINUP);
      fDeadChannelMapHistogramPtrs[gain] = new TH2D(tmpHistoName, tmpHistoName,  
				      NXCOLUMNSRCU , 0, NXCOLUMNSRCU , 
				      NZROWSRCU,         0, NZROWSRCU);
      fDeadChannelMapHistogramPtrs[gain]->Reset(); 
      //     fgCalibHistPtr[gain]->GetXaxis()->SetRange(128, 128 + 64);
    }

  for(int x = 0; x < NXCOLUMNSRCU; x ++)
    {
      for(int z = 0; z < NZROWSRCU; z ++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      geomx = x + NXCOLUMNSRCU*fRcuX;
	      geomz = z + NZROWSRCU*fRcuZ;
	      fEnergyAverageValues[x][z][gain] = 0; 
	      fAccumulatedValues[x][z][gain]   = 0;
	      fTimingAverageValues[x][z][gain] = 0; 
	      fHits[x][z][gain]                = 0;
	      fDeadChannelMap[x][z][gain]      = 0;
	      sprintf(tmpHistoName, "Edistribution_%d_%d_%d_%d",(int)fModuleID,  geomx, geomz, gain);
	      //	      fEnergyHistogramPtrs[x][z][gain] = 0;
	      fEnergyHistogramPtrs[x][z][gain] = new TH1F( tmpHistoName, tmpHistoName, NBINS, XBINLOW, XBINUP);
	      sprintf(tmpHistoName, "TOFdistribution_module%d_x%d_z%d_gain%d", (int)fModuleID,  geomx, geomz, gain);
	      
	      //	      fTimingHistogramPtrs[x][z][gain] = 0;
	      
	      //	      sprintf(tmpHistoName, "DeadChanneMap_%d_%d_%d_%d",(int)fModuleID,  geomx, geomz, gain);
	      //	      fDeadChannelMapHistogramPtrs[x][z][gain] = new TH1D( tmpHistoName, tmpHistoName, N_BINS, XBIN_LOW, XBIN_UP);


	      fTimingHistogramPtrs[x][z][gain] = new TH1F(tmpHistoName , tmpHistoName, NBINS, XBINLOW, XBINUP);
	      fCellAccEnergy.fAccumulatedEnergies[x][z][gain] = 0;
	      fCellAccEnergy.fHits[x][z][gain] = 0;
	      fCellAccEnergy.fDeadChannelMap[x][z][gain] = 0;
	      fCellAccEnergy.fModuleID = 0;
	      fCellAccEnergy.fRcuX = 0;
	      fCellAccEnergy.fRcuZ = 0; 
	    }
	} 
    } 
}



void 
AliHLTPHOSRcuHistogramProducer::FillEnergy(AliHLTUInt8_t x, AliHLTUInt8_t z,  AliHLTUInt8_t gain, float energy)
{
  //comment
  fCellAccEnergy.fAccumulatedEnergies[x][z][gain] += energy;
  fCellAccEnergy.fHits[x][z][gain] ++;
  fEnergyHistogramPtrs[x][z][gain]->Fill(energy); 
}


void 
AliHLTPHOSRcuHistogramProducer::FillTime(AliHLTUInt8_t x, AliHLTUInt8_t z, AliHLTUInt8_t gain, float time)
{
  //See header file for documentation
  fTimingHistogramPtrs[x][z][gain]->Fill(time);
}

const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct& 
AliHLTPHOSRcuHistogramProducer::GetCellAccumulatedEnergies()
{
  //  return &fCellAccEnergy ;
  return fCellAccEnergy ;
}

void 
AliHLTPHOSRcuHistogramProducer::FillLiveChannels(Int_t data[], int size,  Int_t x, Int_t z, Int_t gain)
{
  //comment
  for(Int_t i = 0; i < size; i++)
    {
      if(data[i] > THRESHOLD)
	{
	  if(data[i+1] > THRESHOLD) 
	    {
	      if(data[i+2] > THRESHOLD)
		{
		  if(data[i+3] > THRESHOLD)
		    {
		      
		      fCellAccEnergy.fDeadChannelMap[x][z][gain] = 10;
		      
		      return;
		    }
		}
	    }
	}
    }
}

void 
AliHLTPHOSRcuHistogramProducer::FillLiveChannelHistograms()
{
  //comment
  for(int x = 0; x <  NXCOLUMNSRCU; x ++)
    {
      for(int z = 0; z < NZROWSRCU; z ++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      fDeadChannelMapHistogramPtrs[gain]->SetBinContent(x ,z , fCellAccEnergy.fDeadChannelMap[x][z][gain]);
	    }
  	} 
    }

}


void
AliHLTPHOSRcuHistogramProducer::Reset()
{
  //See header file for documentation
  for(int x = 0; x < NXCOLUMNSRCU; x ++)
    {
      for(int z = 0; z < NZROWSRCU; z ++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      fEnergyAverageValues[x][z][gain] = 0; 
	      fAccumulatedValues[x][z][gain]   = 0;
	      fTimingAverageValues[x][z][gain] = 0; 
	      fHits[x][z][gain]                = 0;
	      fDeadChannelMap[x][z][gain]      = 0;
	    }
	} 
    }
  
  for(int i = 0; i <ALTROMAXSAMPLES;  i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSRcuHistogramProducer::WriteAllHistograms(char *opt)
{
  //comment
  printf("\nAliHLTPHOSRcuHistogramProducer::WriteAllHistogram, opt = %s\n", opt);

  FillLiveChannelHistograms();
  //  cout <<<< endl;

  int runNumber = 0;
  char tmpEFileName[256];
  char tmpDeadFileName_gain0[256];
  char tmpDeadFileName_gain1[256]; 
 // char *tmpRundir = getenv("HOME");
  char runNumberFile[256]; 
  char timeString[256];

  fUtilitiesPtr->ResetArray(runNumberFile, 256);
  fUtilitiesPtr->ResetArray(tmpEFileName, 256);
  fUtilitiesPtr->ResetArray(timeString, 256);
  

  sprintf(runNumberFile, "%s/rundir/runNumber.txt", getenv("HOME"));

  FILE *fp = fopen(runNumberFile, "r");

  Int_t res = 0; //OD to get rid of warnings
  if(fp == 0)
    {
      ScanTimeString(timeString);  
      //      cout << "WARNING, could not find file "<< runNumberFile  <<endl;
      //      cout <<"Filename will be stamped with data and time instead " << endl;
      sprintf(tmpEFileName, "%s/Energy/EnergyHistograms_%s_mod%d_rcuZ%d_rcuX%d.root", fHistoOutDir, timeString, (int)fModuleID, (int)fRcuZ, (int)fRcuZ);
      sprintf(tmpDeadFileName_gain0, "%s/DeadMap/DeadChannelHistograms_%s_mod%d_rcuZ%d_rcuX%d_LG.root", fHistoOutDir, timeString, (int)fModuleID, (int)fRcuZ, (int)fRcuZ); 
      sprintf(tmpDeadFileName_gain1, "%s/DeadMap/DeadChannelHistograms_%s_mod%d_rcuZ%d_rcuX%d_HG.root", fHistoOutDir, timeString, (int)fModuleID, (int)fRcuZ, (int)fRcuZ); 
    } 
  else
    {
      res = fscanf(fp, "%d", &runNumber);
      sprintf(tmpEFileName, "%s/Energy/EnergyHisttograms_run%d_mod%d_rcuZ%d_rcuX%d.root", fHistoOutDir, runNumber, (int)fModuleID, (int)fRcuZ, (int)fRcuX);
      sprintf(tmpDeadFileName_gain0, "%s/DeadMap/DeadChannleHistograms_run%d_mod%d_rcuZ%d_rcuX%d_LG.root", fHistoOutDir, runNumber, (int)fModuleID, (int)fRcuZ, (int)fRcuX);
      sprintf(tmpDeadFileName_gain1, "%s/DeadMap/DeadChannleHistograms_run%d_mod%d_rcuZ%d_rcuX%d_HG.root", fHistoOutDir, runNumber, (int)fModuleID, (int)fRcuZ, (int)fRcuX);
      fclose(fp);
    }

  //  cout << "tmpEFileName = "<< tmpEFileName  << endl;

  TFile *energyHistoFile =  new TFile(tmpEFileName, opt);
  if(!energyHistoFile) return;
  if(!energyHistoFile->IsOpen()) return;

  //  cout <<"printing histograms"<< endl;

  for(int x = 0; x <  NXCOLUMNSRCU; x ++)
    {
      for(int z = 0; z < NZROWSRCU; z ++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      //     cout << "the number of entries is " <<fEnergyHistogramPtrs[x][z][gain]->GetEntries()<< endl;
	      fEnergyHistogramPtrs[x][z][gain]->Write();
	    }
	} 
    }
  energyHistoFile->Close();


  //LOW GAIN
  //TFile *deadHistoFile_gain0 =  new TFile(tmpDeadFileName_gain0,"update");
  TFile *deadHistoFile_gain0 =  new TFile(tmpDeadFileName_gain0, opt);

  if(!deadHistoFile_gain0) return;
  if(!deadHistoFile_gain0->IsOpen()) return;
  fDeadChannelMapHistogramPtrs[0]->Write();


  //HIGH GAIN
  //  TFile *deadHistoFile_gain1 =  new TFile(tmpDeadFileName_gain1,"update");
  TFile *deadHistoFile_gain1 =  new TFile(tmpDeadFileName_gain1, opt);
  if(!deadHistoFile_gain1) return;
  if(!deadHistoFile_gain1->IsOpen()) return;
  fDeadChannelMapHistogramPtrs[1]->Write();


  deadHistoFile_gain0->Close();
  deadHistoFile_gain1->Close(); 
  
  //  cout << "printing histograms, finished"<< endl;
}

void
AliHLTPHOSRcuHistogramProducer::ScanTimeString(char *timeString)
{
  //comment
  time_t timePtr;
  tm *tmPtr;
  time(&timePtr); 
  tmPtr=localtime(&timePtr);
  timeString=asctime(tmPtr);
  char day[10];
  char month[10];
  int date;
  int hour;
  int min;
  int sec;
  int year;
  sscanf(timeString, "%s %s %d %d:%d:%d %d\n", day, month, &date, &hour, &min, &sec, &year);
  
}
