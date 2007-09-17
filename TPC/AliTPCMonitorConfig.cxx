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

/*
$Log$
*/ 

#include "AliTPCMonitorConfig.h"
ClassImp(AliTPCMonitorConfig)

//_______________________________________________________________________________________________________________
AliTPCMonitorConfig::AliTPCMonitorConfig(Char_t* name, Char_t* title) : TNamed(name,title)
{

  // Constructor 
  // Set default values for Size of Window 
  
  fRangeMaxAdcMin       = 50  ;
  fRangeMaxAdcMax       = 100 ;
  
  fRangeBaseMin         = 300 ;
  fRangeBaseMax         = 600 ;
  
  fRangeSumMin          = 50  ; 
  fRangeSumMax          = 100 ;
  
  
  fFormat               = -1;  
  fSector               = 0;
  fEventNext            = 1;
  fEventNextID          = 1;

  SetMainSize(130,720,10,26 );

  fFileLast             = new Char_t[256]       ;
  fFileLastSet          = 0;
  fFileCurrent          = new Char_t[256];
  
  fSectorLast           =-1;
  fSectorLastDisplayed  =-1;
  
  fWrite10Bit           = 0;
  
  fSectorArr            = new Int_t[36];    for(Int_t i =0; i<36; i++) { fSectorArr[i]  =  0;}
  fComponents           = new Float_t[10];  for(Int_t i =0; i<10;i++)  { fComponents[i] =0.0;}
  
  fSamplingFreq         = 1000000;
  fTimeBins             = 1024;
  fPedestals            = 1;
  fNumOfChannels        = 16000;
  fMaxHwAddr            = 24000;

  fFitPulse             = 1;
  fEventProcessed       = 0;
  
} 
 

//_______________________________________________________________________________________________________________
AliTPCMonitorConfig::~AliTPCMonitorConfig() 
{
  // Destructor
} 


//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::SetMainSize(Int_t mainx, Int_t mainy, Int_t borderx=10, Int_t bordery=26) 
{
  // Set Size of Main window, buttons and canvases

  
  // size of main frame ///////////////////////
  fMainXSize            = mainx;
  fMainYSize            = mainy;
  
  // border size depending on window manager///
  fBorderXSize          = borderx;
  fBorderYSize          = bordery;
  
  // y-pos of first sector button /////////////
  fButtonFirstY         = (int) (0.50*fMainYSize);
 
  // canvas size for Minitor Canvases ///////////////////////
  fCanvasMainSize       = (Int_t)(fMainYSize/3);
  
  fCanvasXSize          = (int)fCanvasMainSize;
  fCanvasYSize          = (int)fCanvasMainSize;
  fCanvasXSpace         = (int)(-1*(fCanvasMainSize+ fBorderXSize ));
  fCanvasYSpace         = (int)(    fCanvasMainSize+ fBorderYSize );
  fCanvasXOffset        = (int)(-1*(fMainXSize     + fBorderXSize ));
  
  // adjust main size in y ////////////////////
  fMainYSize            = fMainYSize+fBorderYSize;
  
  // y-size of the main window also defines ///
  fButtonXSize          = (Int_t)(fMainXSize/2-10);
  fButtonYSize          = (Int_t)fMainYSize/47;
  fButtonFirstX1        =         5;
  fButtonFirstX2        = fMainXSize/2+5;

}
 
//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::SetBaseConfig(Float_t* conf_arr)
{
  // Set base configuration stored in array conf_arr 
  
  fRangeBaseMin   = (int)conf_arr[0];
  fRangeBaseMax   = (int)conf_arr[1];
  
  fRangeMaxAdcMin = (int)conf_arr[2];
  fRangeMaxAdcMax = (int)conf_arr[3];
  
  fRangeSumMin    = (int)conf_arr[4];
  fRangeSumMax    = (int)conf_arr[5];
  
  cout << " Set Ranges to : " << endl;
  cout << " range  base       :: " <<  fRangeBaseMin   << "\t :  " << fRangeBaseMax   << endl;
  cout << " range  adc max    :: " <<  fRangeMaxAdcMin << "\t :  " << fRangeMaxAdcMax << endl;
  cout << " range  sum        :: " <<  fRangeSumMin    << "\t :  " << fRangeSumMax    << endl;
  
}

//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::ReadConfig(Char_t* nameconf) 
{
  // Read base configuration from file
  // Update main window size

  string line;
  ifstream datin;
  datin.open(nameconf);
  
  if(!datin) {  AliWarning("Could not read configfile");}
  
  cout << "////  Read Configuration ///////// " << endl;
  while(!datin.eof())
    {
      string line;
      getline(datin,line);
      if(line.find("max adc")!=string::npos)
	{ 
	  datin >> fRangeMaxAdcMin     ; 
	  datin >> fRangeMaxAdcMax     ;
	  cout << " range max          :: " << fRangeMaxAdcMin << " : " << fRangeMaxAdcMax << endl;
	}
      
      if(line.find("baseline")!=string::npos)
	{ 
	  datin >> fRangeBaseMin  ;
	  datin >> fRangeBaseMax  ;
	  cout << " range  base        :: " <<  fRangeBaseMin << " : " <<  fRangeBaseMax  << endl;
	}
      if(line.find("adc sum")!=string::npos)
	{ 
	  datin >> fRangeSumMin  ;
	  datin >> fRangeSumMax  ;
	  cout << " range sum          :: " <<  fRangeSumMin << " : " << fRangeSumMax  << endl;
	}
      if(line.find("frequency")!=string::npos)
	{ 
	  datin >> fSamplingFreq  ;
	  cout << " sampling frequency :: " <<  fSamplingFreq << endl;
	}
      if(line.find("timebins")!=string::npos)
	{ 
	  datin >> fTimeBins  ;
	  cout << " timebins           :: " <<  fTimeBins << endl;
	}
      if(line.find("pedestal")!=string::npos)
	{ 
	  datin >> fPedestals  ;
	  cout << " pedestal scheme    :: " <<  fPedestals << endl;
	}
      if(line.find("main window size")!=string::npos)
	{ 
	  datin >> fMainXSize  ;
	  datin >> fMainYSize  ;
	  cout << " main window size   :: " <<  fMainXSize  << "  , " << fMainYSize <<  endl;
	}
      if(line.find("border size")!=string::npos)
	{ 
	  datin >> fBorderXSize  ; 
	  datin >> fBorderYSize  ;
	  cout << " border size        :: " <<  fBorderXSize  << "  , " << fBorderYSize <<  endl;
	}
    }
  cout << "////  Read Configuration done //// " << endl;
  SetMainSize(fMainXSize,fMainYSize,fBorderXSize,fBorderYSize) ;
}

//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::PrintConfig()
{
  // Print base configuration 
  
  cout << " /////// Configuration /////////////////////////////// " << endl;
  cout << " Timebins                :: " <<  fTimeBins                                   << endl;
  cout << " Range to det ADC max    :: " <<  fRangeMaxAdcMin << " : " << fRangeMaxAdcMax << endl;
  cout << " Range to det Baseline   :: " <<  fRangeBaseMin   << " : " << fRangeBaseMax   << endl;
  cout << " Range to det ADC sum    :: " <<  fRangeSumMin    << " : " << fRangeSumMax    << endl;
  cout << " Pedestal setting        :: " <<  fPedestals                                  << endl;
  cout << " Sampling Frequency      :: " <<  fSamplingFreq                               << endl;
  cout << " Sector (Last)           :: " <<  fSector <<"  ("<< fSectorLast   <<")"       << endl; 
  cout << " Data Format             :: " <<  fFormat                                     << endl;
  cout << " File current            :: " <<  fFileCurrent                                << endl;
  
  if(fFileLastSet)  {      cout << " File last           :: " <<  fFileLast << endl;    }
  return;										  
  
}

//______________________________________________________________________________________________________________
Char_t* AliTPCMonitorConfig::GetLastProcFile()
{
  // Return name of last processed file
  // Get from file if written  
  if(!fFileLastSet)
    {
      Char_t fnlast[256];
      sprintf(fnlast,"AliTPCMonitorLastFile.txt");
      ifstream datin(fnlast);
      if(!datin) {  AliWarning("Could not read file containing  name of last processed file: Check permissions and path");}
      datin >> fFileLast;
      datin.close(); 
      fFileLastSet=1;
    }
  return fFileLast;
}

//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::SetLastProcFile(Char_t* file)
{
  // Set name of last processed file 
  // Write name to a file 
  
  sprintf(fFileLast,file);
  Char_t fnlast[256];
  sprintf(fnlast,"AliTPCMonitorLastFile.txt");
  ofstream datout(fnlast);
  if(!datout) {  AliWarning("Could not write file containing name of last processed file: Check permissions and path");}
  datout << fFileLast << endl;
  datout.close();
  fFileLastSet =1;
}
