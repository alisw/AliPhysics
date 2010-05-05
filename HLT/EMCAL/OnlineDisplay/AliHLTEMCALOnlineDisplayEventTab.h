
//-*- Mode: C++ -*-
// $Id: AliHLTEMCALOnlineDisplayEventTab.h 35071 2009-09-29 05:26:09Z phille $

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
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

#ifndef ALIHLTEMCALONLINEDISPLAYEVENTTAB_H
#define ALIHLTEMCALONLINEDISPLAYEVENTTAB_H


#include "AliHLTEMCALOnlineDisplayTab.h"
#include "AliHLTEMCALConstants.h"
#include "AliHLTEMCALMapper.h"


#define NZRCUCOORD 2
#define NXRCUCOORD 2


using namespace EmcalHLTConst;


class TGTab;
class TRootEmbeddedCanvas;
class TCanvas;
class TH2D;
class TH1D;
class AliHLTEMCALOnlineDisplayTH2D;
class AliHLTEMCALGetEventButton;
class AliHLTHOMERReader;

//class AliHLTEMCALRcuCellEnergyDataStruct;

class AliHLTCaloRcuCellEnergyDataStruct;


class AliHLTEMCALOnlineDisplay;
//class AliHLTEMCALSharedMemoryInterface;

//class AliHLTEMCALSharedMemoryInterfacev2;
 
//class  AliHLTCaloSharedMemoryInterfacev2;

class  AliHLTEMCALSharedMemoryInterface;


//class AliHLTEMCALChannelRawDataStruct;


class AliHLTCaloChannelRawDataStruct;


class AliHLTEMCALOnlineDisplayEventTab : public AliHLTEMCALOnlineDisplayTab
{
 public:
 
  virtual ~AliHLTEMCALOnlineDisplayEventTab();


  AliHLTEMCALOnlineDisplayEventTab(AliHLTEMCALOnlineDisplay * onlineDisplayPtr, TGTab  *tabPtr, 
				  AliHLTHOMERReader * homerSyncPtr, 
				  AliHLTHOMERReader * homerPtrs[MAXHOSTS], 
				  int nHosts,  int runnumber = -1);
//   AliHLTEMCALOnlineDisplayEventTab(AliHLTEMCALOnlineDisplay *onlineDisplayPtr, TGTab *tabPtr, 
// 				  AliHLTHOMERReader *fgHomerReaderPtr, 
// 				  AliHLTHOMERReader *fgHomerReadersPtr[MAXHOSTS], 
// 				  int nHosts, const int runnumber = -1);
    //    {

 
  
 

/* 
  void SetRunNumber(const int runnumber) 
  {
    
    fRunNumber = runnumber ;
    cout << __FILE__ <<":"<< __LINE__ << "RunNumber was set to "<< fRunNumber  <<endl; ;
  };
  */

  Int_t GetRawData(TH1D *histPtr, int x, int z, int gain);

  void UpdateDisplay();
  int GetNextEvent();
  virtual void ReadBlockData(AliHLTHOMERReader *homeReaderPtr);
  void FindFourierBlocks(AliHLTHOMERReader *homeReaderPtr) const;

  void ResetDisplay();
  TGTab               *fTab;
  TGTab               *fSubTab1;
  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6;
  TGCompositeFrame    *fSubF1, *fSubF2, *fSubF3;
  TCanvas *fgCanvasPtr[NGAINS];
  AliHLTEMCALOnlineDisplayTH2D *fgLegoPlotPtr[NGAINS];

  /*
  int *fChannelData[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  Int_t fNChannelSamples[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  Int_t fChannelEnergy[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  */  

  int *fChannelData[NMODULES][NZROWSMOD][NXCOLUMNSMOD][NGAINS];
  Int_t fNChannelSamples[NMODULES][NZROWSMOD][NXCOLUMNSMOD][NGAINS];
  Int_t fChannelEnergy[NMODULES][NZROWSMOD][NXCOLUMNSMOD][NGAINS];
  
 private:
  AliHLTEMCALOnlineDisplayEventTab();
  
  
  // void FillRawData(const AliHLTEMCALChannelRawDataStruct &rawStr);
  void FillRawData(const AliHLTCaloChannelRawDataStruct &rawStr);
  AliHLTEMCALGetEventButton* fgEventButtPtr; 
  void InitDisplay(TGTab *tabPtr){};
  void InitDisplay(TGTab * tabPtr, int runnumber);
  AliHLTEMCALOnlineDisplay *fOnlineDisplayPtr;
  // AliHLTCaloSharedMemoryInterfacev2  *fShmPtr;   
  AliHLTEMCALSharedMemoryInterface  *fShmPtr;   

};


#endif
