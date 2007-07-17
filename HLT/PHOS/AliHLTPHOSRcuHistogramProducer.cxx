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
//#include "stdio.h"
#//include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "TFile.h"

using  namespace std;


/*************************************************************************
* Class AliHLTPHOSRcuHistogramProducer accumulating histograms           *
* with amplitudes per PHOS channel                                       *
* It is intended to run at the HLT farm                                  *
* and it fills the histograms with amplitudes per channel.               * 
* Usage example see in PHOS/macros/Shuttle/AliPHOSCalibHistoProducer.C   *
**************************************************************************/


AliHLTPHOSRcuHistogramProducer:: AliHLTPHOSRcuHistogramProducer(): fModuleID(0), fRcuX(0), fRcuZ(0)

{
  //Default constructor
  cout << "WARNING: You cannot invoke the AliHLTPHOSRcuHistogramProducer without arguments" << endl;
  cout << "Usage AliHLTPHOSRcuHistogramProducer(ModuleID, X. Z)" << endl;
} 

//AliHLTPHOSRcuHistogramProducer::AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ)
AliHLTPHOSRcuHistogramProducer::AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ)
{
  //Se header file for documentation
  SetModuleID(moduleID);
  SetRcuX(rcuX); 
  SetRcuZ(rcuZ); 
  Init();
}

AliHLTPHOSRcuHistogramProducer::~ AliHLTPHOSRcuHistogramProducer()
{
  //Destructor
}


void
AliHLTPHOSRcuHistogramProducer::Init()
{
  //See header file for documentation
  char tmpHistoName[256];
  int geomx;
  int geomz;
    

  for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
    {
      for(int z = 0; z < N_ZROWS_RCU; z ++)
	{
	  for(int gain = 0; gain < N_GAINS; gain ++)
	    {
	      geomx = x + N_XCOLUMNS_RCU*fRcuX;
	      geomz = z + N_ZROWS_RCU*fRcuZ;

	      fEnergyAverageValues[x][z][gain] = 0; 
	      fAccumulatedValues[x][z][gain]   = 0;
	      fTimingAverageValues[x][z][gain] = 0; 
	      fHits[x][z][gain]                = 0;
	      sprintf(tmpHistoName, "Edistribution_%d_%d_%d_%d",(int)fModuleID,  geomx, geomz, gain);
	      fEnergyHistogramPtrs[x][z][gain] = 0;
	      fEnergyHistogramPtrs[x][z][gain] = new TH1F( tmpHistoName, tmpHistoName, N_BINS, XBIN_LOW, XBIN_UP);
	      sprintf(tmpHistoName, "TOFdistribution_%d_%d_%d_%d",(int)fModuleID,  geomx, geomz, gain);
	      fTimingHistogramPtrs[x][z][gain] = 0;
	      fTimingHistogramPtrs[x][z][gain] = new TH1F(tmpHistoName , tmpHistoName, N_BINS, XBIN_LOW, XBIN_UP);
	      fCellAccEnergy.fAccumulatedEnergies[x][z][gain] = 0;
	      fCellAccEnergy.fHits[x][z][gain] = 0;
	      fCellAccEnergy.fModuleID = 0;
	      fCellAccEnergy.fRcuX = 0;
	      fCellAccEnergy.fRcuZ = 0; 
	    }
	} 
    } 
}


void 
AliHLTPHOSRcuHistogramProducer::SetRcuX(AliHLTUInt8_t X)
{
  //See header file for documentation
  fRcuX = X; 
  fCellAccEnergy.fRcuX = X;
}



void 
AliHLTPHOSRcuHistogramProducer::SetRcuZ(AliHLTUInt8_t Z)
{
  //See header file for documentation
  fRcuZ = Z; 
  fCellAccEnergy.fRcuZ = Z;
}




void 
AliHLTPHOSRcuHistogramProducer::SetModuleID(AliHLTUInt8_t moduleID)
{
  //See header file for documentation
 fModuleID = moduleID;
}


void 
AliHLTPHOSRcuHistogramProducer::FillEnergy(AliHLTUInt8_t x, AliHLTUInt8_t z,  AliHLTUInt8_t gain, float energy)
{
  //  cout << "AliHLTPHOSRcuHistogramProducer::FillEnergy x = " << (int)x<< "  z=  " << (int)z;
  //  cout << " Energy =" <<  energy;
  //  cout << " Accumulated energy" <<fCellAccEnergy.fAccumulatedEnergies[x][z][gain] << endl; 
  //See header file for documentation
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
AliHLTPHOSRcuHistogramProducer::Reset()
{
  //See header file for documentation
  for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
    {
      for(int z = 0; z < N_ZROWS_RCU; z ++)
	{
	  for(int gain = 0; gain < N_GAINS; gain ++)
	    {
	      fEnergyAverageValues[x][z][gain] = 0; 
	      fAccumulatedValues[x][z][gain]   = 0;
	      fTimingAverageValues[x][z][gain] = 0; 
	      fHits[x][z][gain]                = 0;
	    }
	} 
    }
  
  for(int i = 0; i <ALTRO_MAX_SAMPLES;  i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSRcuHistogramProducer::WriteEnergyHistograms()
{
  //See header file for documentation
  char tmpFileName[256];
  sprintf(tmpFileName,"/home/aliphoshlt/rundir/outdata/calibHisto_%d_%d_%d.root", (int)fModuleID, (int)fRcuX, (int)fRcuZ);
  TFile *histoFile =  new TFile(tmpFileName,"update");
  char hname[128];
  if(!histoFile) return;
  if(!histoFile->IsOpen()) return;

  cout <<"printing histograms"<< endl;
  cout <<"histofile-Getname() =" << histoFile->GetName() << endl;

    for(int x = 0; x <  N_XCOLUMNS_RCU; x ++)
    {
      for(int z = 0; z < N_ZROWS_RCU; z ++)
	{
	  for(int gain = 0; gain < N_GAINS; gain ++)
	    {
	      fEnergyHistogramPtrs[x][z][gain]->Write();
	    }
	} 
    }

    cout << "printing histograms, finished"<< endl;
    histoFile->Close();

}
