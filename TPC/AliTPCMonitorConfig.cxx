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
Revision 1.2  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

//////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorConfig class
////
//// Configuration handler class for AliTPCMonitor
////
//// The basic configuration will be read from the file AliTPCMonitorConfig.txt
//// and can be changed online via the Button "Conf. Ranges"
//// Basic configuration settings are e.g. the range for the determination 
//// of the baseline, maximum adc value and the settings for the pedestal calculation.
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////


#include "AliTPCMonitorConfig.h"
#include "AliLog.h" 
#include <Riostream.h>

ClassImp(AliTPCMonitorConfig)

// _______________________________________________________________________________________________________________
AliTPCMonitorConfig::AliTPCMonitorConfig(const Char_t* name, const Char_t* title) : 
  TNamed(name,title),
  fFormat(-1),
  fSector(0),
  fSectorLast(-1),
  fSectorLastDisplayed(-1),
  fSectorArr(new Int_t[36]),
  fFileLast(),
  fFileLastSet(0),
  fFileCurrent(),
  fEventNext(1),
  fEventNextID(1),
  fEventProcessed(0),
  fRangeMaxAdcMin(50),
  fRangeMaxAdcMax(100),
  fRangeBaseMin(300),
  fRangeBaseMax(600),
  fRangeSumMin(50), 
  fRangeSumMax(100),
  fCanvasXSize(100),
  fCanvasYSize(100), 
  fCanvasXSpace(10), 
  fCanvasYSpace(10), 
  fCanvasXOffset(130),
  fCanvasMainSize(200),
  fMainXSize(100),
  fMainYSize(600),
  fBorderXSize(10),
  fBorderYSize(10),
  fButtonXSize(100), 
  fButtonYSize(20), 
  fButtonFirstX1(10),
  fButtonFirstX2(50),
  fButtonFirstY(300),  
  fWrite10Bit(0),
  fComponents(new Float_t[10]), 
  fSamplingFreq(1000000),
  fPedestals(1),
  fNumOfChannels(16000),
  fTimeBins(1024),
  fMaxHwAddr(24000),
  fFitPulse(1),
  fProcOneSector(0)
{
  // Constructor 
  for(Int_t i =0; i<36; i++) { fSectorArr[i]  =  0;}
  for(Int_t i =0; i<10;i++)  { fComponents[i] =0.0;}
  SetMainSize(130,720,10,26 );
  
} 
 

//_______________________________________________________________________________________________________________
AliTPCMonitorConfig::AliTPCMonitorConfig(const AliTPCMonitorConfig &config) :
  TNamed(config.GetName(),config.GetTitle()),
  fFormat(config.fFormat),
  fSector(config.fSector),
  fSectorLast(config.fSectorLast),
  fSectorLastDisplayed(config.fSectorLastDisplayed),
  fSectorArr(new Int_t[36]),
  fFileLast(config.fFileLast),
  fFileLastSet(config.fFileLastSet),
  fFileCurrent(config.fFileCurrent),
  fEventNext(config.fEventNext),
  fEventNextID(config.fEventNextID),
  fEventProcessed(config.fEventProcessed),
  fRangeMaxAdcMin(config.fRangeMaxAdcMin),
  fRangeMaxAdcMax(config.fRangeMaxAdcMax),
  fRangeBaseMin(config.fRangeBaseMin),
  fRangeBaseMax(config.fRangeBaseMax),
  fRangeSumMin(config.fRangeSumMin),
  fRangeSumMax(config.fRangeSumMax),
  fCanvasXSize(config.fCanvasXSize),
  fCanvasYSize(config.fCanvasYSize),
  fCanvasXSpace(config.fCanvasXSpace),
  fCanvasYSpace(config.fCanvasYSpace),
  fCanvasXOffset(config.fCanvasXOffset),
  fCanvasMainSize(config.fCanvasMainSize),
  fMainXSize(config.fMainXSize),
  fMainYSize(config.fMainYSize),
  fBorderXSize(config.fBorderXSize),
  fBorderYSize(config.fBorderYSize),
  fButtonXSize(config.fBorderXSize),
  fButtonYSize(config.fButtonYSize),
  fButtonFirstX1(config.fButtonFirstX1),
  fButtonFirstX2(config.fButtonFirstX2),
  fButtonFirstY(config.fButtonFirstY),
  fWrite10Bit(config.fWrite10Bit),
  fComponents(new Float_t[10]),
  fSamplingFreq(config.fSamplingFreq),
  fPedestals(config.fPedestals),
  fNumOfChannels(config.fNumOfChannels),
  fTimeBins(config.fTimeBins),
  fMaxHwAddr(config.fMaxHwAddr),
  fFitPulse(config.fFitPulse),
  fProcOneSector(config.fProcOneSector)
{
  // copy constructor
  
  for(Int_t i =0; i<36; i++) { fSectorArr[i]  =  0;}
  for(Int_t i =0; i<10;i++)  { fComponents[i] =0.0;}


}
//_______________________________________________________________________________________________________________
AliTPCMonitorConfig &AliTPCMonitorConfig::operator =(const AliTPCMonitorConfig& config)
{
  // assignement operator
  if(this!=&config){ 
    ((TNamed *)this)->operator=(config);
    fFormat=config.fFormat;
    fSector=config.fSector;
    fSectorLast=config.fSectorLast;
    fSectorLastDisplayed=config.fSectorLastDisplayed;
    fFileLastSet=config.fFileLastSet;
    fEventNext=config.fEventNext;
    fEventNextID=config.fEventNextID;
    fEventProcessed=config.fEventProcessed;
    fRangeMaxAdcMin=config.fRangeMaxAdcMin;
    fRangeMaxAdcMax=config.fRangeMaxAdcMax;
    fRangeBaseMin=config.fRangeBaseMin;
    fRangeBaseMax=config.fRangeBaseMax;
    fRangeSumMin=config.fRangeSumMin;
    fRangeSumMax=config.fRangeSumMax;
    fCanvasXSize=config.fCanvasXSize;
    fCanvasYSize=config.fCanvasYSize;
    fCanvasXSpace=config.fCanvasXSpace;
    fCanvasYSpace=config.fCanvasYSpace;
    fCanvasXOffset=config.fCanvasXOffset;
    fCanvasMainSize=config.fCanvasMainSize;
    fMainXSize=config.fMainXSize;
    fMainYSize=config.fMainYSize;
    fBorderXSize=config.fBorderXSize;
    fButtonYSize=config.fButtonYSize;
    fButtonFirstX1=config.fButtonFirstX1;
    fButtonFirstX2=config.fButtonFirstX2;
    fButtonFirstY=config.fButtonFirstY;
    fWrite10Bit=config.fWrite10Bit;
    fPedestals=config.fPedestals;
    fNumOfChannels=config.fNumOfChannels;
    fTimeBins=config.fTimeBins;
    fMaxHwAddr=config.fMaxHwAddr;
    fFitPulse=config.fFitPulse; 
    fProcOneSector=config.fProcOneSector;

    fFileLast             = config.fFileLast;
    fFileCurrent          = config.fFileCurrent;
    
    fSectorArr            = new Int_t[36];    for(Int_t i =0; i<36; i++) { fSectorArr[i]  =  0;}
    fComponents           = new Float_t[10];  for(Int_t i =0; i<10;i++)  { fComponents[i] =0.0;}
    
      
  }
  return *this;
}


//_______________________________________________________________________________________________________________
AliTPCMonitorConfig::~AliTPCMonitorConfig() 
{
  // Destructor
  delete[] fSectorArr;
  delete[] fComponents;
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
void AliTPCMonitorConfig::SetBaseConfig(Float_t* confarr)
{
  // Set base configuration stored in array confarr 
  
  fRangeBaseMin   = (int)confarr[0];
  fRangeBaseMax   = (int)confarr[1];
  
  fRangeMaxAdcMin = (int)confarr[2];
  fRangeMaxAdcMax = (int)confarr[3];
  
  fRangeSumMin    = (int)confarr[4];
  fRangeSumMax    = (int)confarr[5];
  
  cout << " Set Ranges to : " << endl;
  cout << " range  base       :: " <<  fRangeBaseMin   << "\t :  " << fRangeBaseMax   << endl;
  cout << " range  adc max    :: " <<  fRangeMaxAdcMin << "\t :  " << fRangeMaxAdcMax << endl;
  cout << " range  sum        :: " <<  fRangeSumMin    << "\t :  " << fRangeSumMax    << endl;
  
}

//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::ReadConfig(const Char_t* nameconf) 
{
  // Read base configuration from file
  // Update main window size

  //  string line;
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
  
  if(fFileLastSet)  {      cout << " File last           :: " <<  fFileLast.Data() << endl;    }
  return;										  
  
}

//______________________________________________________________________________________________________________
const Char_t* AliTPCMonitorConfig::GetLastProcFile()
{
  // Return name of last processed file
  // Get from file if written  
  if(!fFileLastSet)
    {
      ifstream datin("AliTPCMonitorLastFile.txt");
      if(!datin.is_open()) {
        AliWarning("Could not read file containing  name of last processed file: Check permissions and path");
        fFileLast="";
        fFileLastSet=0;
        return fFileLast.Data();
      }
      datin >> fFileLast;
      datin.close(); 
      fFileLastSet=1;
    }
  return fFileLast;
}

//_______________________________________________________________________________________________________________
void AliTPCMonitorConfig::SetLastProcFile(const Char_t* file)
{
  // Set name of last processed file 
  // Write name to a file 
  
  fFileLast=file;
  ofstream datout("AliTPCMonitorLastFile.txt");
  if(!datout) {  AliWarning("Could not write file containing name of last processed file: Check permissions and path");}
  datout << fFileLast << endl;
  datout.close();
  fFileLastSet =1;
}
