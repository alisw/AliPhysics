// $Id: AliHLTEMCALOnlineDisplay.cxx 35108 2009-09-30 01:58:37Z phille $

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


#ifndef __CINT__
# include <stdexcept>
# include <TSystem.h>
# include <TApplication.h>
# include "TStyle.h" 
#endif



#include  "AliHLTEMCALOnlineDisplay.h"
#include  "AliHLTDataTypes.h"
//#include  "AliHLTEMCALRcuCellEnergyDataStruct.h"
#include  <vector>
#include  "stdio.h"
#include <string>
#include <sys/ipc.h>
#include <errno.h>
#include "TH2.h"
#include "TCanvas.h"

//#include "AliHLTCaloRcuCellAccumulatedEnergyDataStruct.h"

//#include "AliHLTEMCALCommonDefs.h"
#include <iostream>
#include "AliHLTEMCALOnlineDisplayEventTab.h"
#include "AliHLTEMCALOnlineDisplayCalibTab.h"
#include "AliHLTEMCALOnlineDisplayFourierTab.h"

//#include "AliHLTEMCALFourier.h"


AliHLTEMCALOnlineDisplayEventTab*     AliHLTEMCALOnlineDisplay::fgEventTabPtr       = 0;
AliHLTEMCALOnlineDisplayFourierTab*   AliHLTEMCALOnlineDisplay::fgFourierTabPtr     = 0;
AliHLTEMCALOnlineDisplayCalibTab*     AliHLTEMCALOnlineDisplay::fgCalibTabPtr       = 0;
AliHLTEMCALOnlineDisplay*             AliHLTEMCALOnlineDisplay::fgInstancePtr       = 0;          /**<The one an only instance of PhosOnlineDisplay*/
HOMERReader*                         AliHLTEMCALOnlineDisplay::fgHomerReaderPtr    = 0;          /**<Homer reader that fetches events from the HLT online stream*/
HOMERReader*                         AliHLTEMCALOnlineDisplay::fgHomerReadersPtr[MAXHOSTS];     /**<Homer reader that fetches events from the HLT online stream*/
Bool_t                               AliHLTEMCALOnlineDisplay::fgAccumulate        = kFALSE ;    /**<If set to kFALSE reset fgLegoplot between event, kTRUE adds current energies to previous plot*/
Bool_t                               AliHLTEMCALOnlineDisplay::fgSyncronize        = kFALSE ;
unsigned int                         AliHLTEMCALOnlineDisplay::fgNHosts            = 0;
unsigned int                         AliHLTEMCALOnlineDisplay::fgNPorts            = 0;
char*                                AliHLTEMCALOnlineDisplay::fgHosts[MAXHOSTS];
short unsigned int*                  AliHLTEMCALOnlineDisplay::fgPorts             = 0; 
TGTab*                               AliHLTEMCALOnlineDisplay::fgTab                = 0;

using namespace std;

//gStyle->SetOptStat(false);


AliHLTEMCALOnlineDisplay*
AliHLTEMCALOnlineDisplay::Instance(int argc, char** argv) 
{
  // See header file for documentation
  if (!fgInstancePtr) fgInstancePtr = new AliHLTEMCALOnlineDisplay(argc, argv);
  return fgInstancePtr;
}


//AliHLTEMCALOnlineDisplay::AliHLTEMCALOnlineDisplay() : AliHLTEMCALBase(), fRunNumber(-1)
AliHLTEMCALOnlineDisplay::AliHLTEMCALOnlineDisplay() : fRunNumber(-1)

{
  // See header file for documentation
  cout << "ERROR ! level: FATAL, you cannot invoke the onlinedisplay without arguments" << endl;
}


//AliHLTEMCALOnlineDisplay::AliHLTEMCALOnlineDisplay(int argc, char** argv) : AliHLTEMCALBase()
AliHLTEMCALOnlineDisplay::AliHLTEMCALOnlineDisplay(int argc, char** argv) 
{
  // See header file for documentation
  gStyle->SetOptStat(false);
  
  ScanArguments(argc, argv);
  char **tmp;
  cout << "creating new EMCAL Onlinedisplay" << endl;
  fgHomerReaderPtr = new  HOMERReader(fgNHosts, (const char**)fgHosts,  fgPorts);
  cout << "AliHLTEMCALOnlineDisplay::AliHLTEMCALOnlineDisplay: fgHomerReaderPtr =  " <<  fgHomerReaderPtr << endl;

  for(int i = 0; i <fgNHosts; i++)
    {
      fgHomerReadersPtr[i] =      new  HOMERReader(fgHosts[i], fgPorts[i]); 
    }
  InitDisplay();
}


AliHLTEMCALOnlineDisplay::~AliHLTEMCALOnlineDisplay()
{
  // See header file for documentation
}



void
AliHLTEMCALOnlineDisplay::InitDisplay()
{
  // See header file for documentation
  char tmpHistoName[256];
  char tmpChDtaName[256];

  gStyle->SetPalette(1);
  fgTab = new TGTab(this, 100, 100);
  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);
  fgEventTabPtr =    new  AliHLTEMCALOnlineDisplayEventTab(this, fgTab, fgHomerReaderPtr, fgHomerReadersPtr, fgNHosts, fRunNumber);
  fgFourierTabPtr =  new  AliHLTEMCALOnlineDisplayFourierTab(this, fgTab, fgHomerReaderPtr, fgHomerReadersPtr, fgNHosts); 
  fgCalibTabPtr =    new  AliHLTEMCALOnlineDisplayCalibTab(fgTab, fgHomerReaderPtr, fgHomerReadersPtr, fgNHosts);
  
  //  fgEventTabPtr->SetRunNumber(fRunNumber);
  // fgFourierTabPtr->SetRunNumber(fRunNumber);
  // fgCalibTabPtr->SetRunNumber(fRunNumber);

  AddFrame(fgTab, fL1);
  MapSubwindows();
  Resize();
  SetWindowName("EMCAL HLT OnlineDisplay");
  MapWindow();
  
  // MoveResize(100,100, 800,1000);
  //  MoveResize(100,100, 500,700);
  
  MoveResize(100,100, 700, 800 );

}



int
AliHLTEMCALOnlineDisplay::GetNextEvent()
{
  // See header file for documentation
  cout << "AliHLTEMCALOnlineDisplay::GetNextEvent()" << endl;
  fgEventTabPtr->GetNextEvent();
  fgFourierTabPtr->GetNextEvent();
}


int
AliHLTEMCALOnlineDisplay::GetHistogram()
{
  // See header file for documentation
  fgCalibTabPtr->GetNextEvent(); 
}



int
AliHLTEMCALOnlineDisplay::ScanArguments(int argc, char** argv)
{
  // See header file for documentation
  for(int i=0; i< MAXHOSTS; i++)
    {
      fgHosts[i] = new char[256];
    }

  fgPorts = new short unsigned[100];
  Bool_t hostIsSet = kFALSE;
  Bool_t portIsSet = kFALSE;
  int iResult=0;
  TString argument="";

  for (int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      
      if (argument.IsNull()) 
	{
	  continue;
	}
      
      if (argument.CompareTo("-sync")==0)
	{
	  cout << "setting Synchronize to true" << endl;
	  fgSyncronize = kTRUE;
	}  
      
      if (argument.CompareTo("-acc")==0)
	{
	  cout << "setting Accumulate to true" << endl;
	  fgAccumulate = kTRUE;
	}  
      
      if (argument.CompareTo("-run")==0) 
	{
	  if(i+1 <= argc)
	    {
	      i++;
	      fRunNumber =  atoi(argv[i]);
	      cout << __FILE__ <<":" <<__LINE__ << ", !!!!!!!!!!!!! setting runnumber too   " <<  fRunNumber   <<endl;;
	      //	      fIsSetRunNumber = true;
	    }

	}


      if (argument.CompareTo("-host")==0) 
	{
	  if(i+1 <= argc)
	    {
	      i++;
	      sprintf(fgHosts[fgNHosts],"%s", argv[i]);
	      fgNHosts ++; 
	      cout <<"fgNHosts set to "<< fgNHosts <<endl;
	      hostIsSet = kTRUE; 
	      if(i+1 <= argc)
		{
		  argument=argv[i+1];
		  if(argument.CompareTo("-port")==0)
		    {
		      i++;
		      if(i+1 <= argc)
			{
			  i++;
			  fgPorts[fgNPorts] = atoi(argv[i]);	
			  cout << "A setting port to   " << fgPorts[fgNPorts]  <<endl; 
			  fgNPorts ++;
			  portIsSet = kTRUE;
			}
		    }
		  else
		    {
		      fgPorts[fgNPorts] =  DEFAULTEVENTPORT;	
		      cout << "B setting port to   " << fgPorts[fgNPorts]  <<endl; 
		      fgNPorts ++;
		      portIsSet = kTRUE;
		    }
		}
	    }
	}
    }

  if(hostIsSet != kTRUE ||  portIsSet != kTRUE)
    {
      if(hostIsSet == kFALSE)
	{
	  printf("\nERROR: no hostname is specified\n");
	}
      
      if( portIsSet == kFALSE)
	{
	  printf("ERROR: no port spcified\n");
	}
      printf("\nYou must specify at least one host \n\n");
      printf("*****************************************************************\n");
      printf("\nUsage: ./onlinedisplay  -host  <hostname>   -port  <port>");
      printf("\n-port is optional, if not set  port 42001 will be used\n");
      printf("*****************************************************************\n\n\n");
      iResult = -1;
    }
 
 
  else
    {
      iResult = 0;
    }

  return iResult;
}//end ScanArguments



void 
AliHLTEMCALOnlineDisplay::Gain2Text(const int gain,  char *txt) const
{
  // See header file for documentation
  if(gain == LOWGAIN)
    {
      sprintf(txt,"High Gain");
      //   sprintf(txt,"Low Gain");
    }

  else if(gain == HIGHGAIN)
    {
      sprintf(txt,"Low Gain");
      //  sprintf(txt,"High Gain");
    }

  else
    {
      sprintf(txt,"Error!! invalid gain %d", gain);
    }
  
}
