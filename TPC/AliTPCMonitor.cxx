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
Revision 1.3  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.2  2007/09/18 09:44:45  cvetan
Sorting out some issues concerning the compilation with and without DATE support

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/   

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitor class
//// 
//// Main class for TPC Monitor
//// Monitor can handle rootified data, files and online streams in DATE format.
//// The monitor GUI is started by the macro TPCMonitor.C
//// 
//// In the The data are read in in two read cycles. 
//// If no sector is specified (sectorid==-1) all sectors are read and only the max value 
//// for each channel is stored in a global histogram for Side A and C.
//// In this way the whole TPC can be read in at once.
////
//// If the sector is specified  only one sector is read in and additional quantities
//// e.g baseline and baseline rms are calculated and stored. 
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////         Jens Wiechula, Uni Tuebingen (Jens.Wiechula@cern.ch)
////
/////////////////////////////////////////////////////////////////////////



#include "AliLog.h"
#include "AliTPCMonitor.h"   
#include "AliTPCMonitorMappingHandler.h"
#include "AliTPCMonitorFFT.h"
#include "AliRawReader.h"
#include "AliRawReaderDateOnline.h"
#include "AliRawReaderRoot.h"
#include "AliRawEventHeaderBase.h"
#include "AliAltroRawStreamV3.h"
#include "TH2F.h" 
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3S.h"
#include "TLegend.h"
#include "TROOT.h" 
#include "TDirectory.h"
#include "TSystem.h"
#include "TString.h"
#include "TPaveText.h"   
#include "TFile.h"
#include <Riostream.h>

ClassImp(AliTPCMonitor)

const Int_t AliTPCMonitor::fgkHwMaskFEC              = 0x0780;                          
const Int_t AliTPCMonitor::fgkHwMaskBranch           = 0x0800;                          
const Int_t AliTPCMonitor::fgkHwMaskFECChannel       = 0x007f;                          
const Int_t AliTPCMonitor::fgkHwMaskAltroChannel     = 0x000f;                          
const Int_t AliTPCMonitor::fgkHwMaskAltroChip        = 0x0070;                          
const Int_t AliTPCMonitor::fgkHwMaskRCU              = 0x7000;                          

//____________________________________________________________________________
AliTPCMonitor::AliTPCMonitor(const char* name, const char* title) : 
AliTPCMonitorConfig(name,title),
fPad(new Int_t*[GetMaxHwAddr()]),
fPadMapHw(new Float_t[GetMaxHwAddr()]),
fPadMapRCU(new Int_t*[GetMaxHwAddr()]),
fHistIROC(0),
fHistOROC(0),
fHistIROCIndex(0),
fHistOROCIndex(0),
fHistIROCTime(0),
fHistOROCTime(0),
fHistIROCClone(0),
fHistOROCClone(0),
fHistIROCRMS(0),
fHistOROCRMS(0),
fHistIROCBASE(0),
fHistOROCBASE(0),
fHistIROCSUM(0),
fHistOROCSUM(0),
fHistChannelTime(0),
fHistAddrMapIndex(0),
fHistAddrMaxAdc(0),
fHistAddrBaseMean(0),
fHistAddrMaxAdcX(0),
fHistAddrAdcSum(0),
fHistAddrBaseRms(0),
fHistDistrSumIROC(0),
fHistDistrMaxIROC(0),
fHistDistrSumOROC(0),
fHistDistrMaxOROC(0),
fHistDistrBase2dIROC(0),
fHistDistrBase2dOROC(0),
fHistDistrBaseRmsIROC(0),
fHistDistrBaseMeanIROC(0),
fHistDistrBaseRmsOROC(0),
fHistDistrBaseMeanOROC(0),
fHistGlobalMaxA(0),
fHistGlobalMaxC(0),
fHistList(new TObjArray()),

fkNRowsIroc(63),
fkNRowsOroc(96),

fkNPadsIroc(110),
fkNPadsOroc(140),
fkNPadMinIroc(-55),
fkNPadMinOroc(-70),
fkNPadMaxIroc(55),
fkNPadMaxOroc(70),
fVerb(0),
fLastEv(0),
fEventNumber(0),
fEventNumberOld(0),
fDisableFit(0),
fExecGlob(0),
fExecPlaneMax(0),
fExecPadIrocRms(0),
fExecPadOrocRms(0),
fRunId(0),
fEqId(0),
fPadUsedRoc(-1),
fPadUsedHwAddr(-1),
fGdcId(0),
fLdcId(0),
fLdcIdOld(1),
fMapEqidsSec(new Int_t*[36]),
fMapEqidsRcu(new Int_t[1000]),
fMirror(1),
fChannelIter(0),
fMapHand(0),
fRawReader(0),
fkMonTable(0x0),
fMonTableString(""),
fMonTableArray(0x0),
fMonTableChanged(kFALSE)

{
  // Constructor
  
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPad[i]       = new Int_t[GetTimeBins()];}
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPadMapRCU[i] = new Int_t[6];}
  for(Int_t i = 0; i<36; i++) { fMapEqidsSec[i]       = new Int_t[6];}
  
  if (gDirectory)
  {
    if (!gDirectory->GetList())
    {
      Warning("Build","Current directory is not a valid directory");
      return;
    }
    AliTPCMonitor *hold = (AliTPCMonitor*)gDirectory->GetList()->FindObject(GetName());
    if(hold)
    {
      Warning("Build","Replacing existing histogram: %s (Potential memory leak).",GetName());
      gDirectory->GetList()->Remove(hold);
    }
    gDirectory->Append(this);
  }
  CreateHistos();
  SetEqIds();
  
}

//____________________________________________________________________________
AliTPCMonitor::AliTPCMonitor(const AliTPCMonitor &monitor):
AliTPCMonitorConfig(monitor.GetName(),monitor.GetTitle()),
fPad(new Int_t*[GetMaxHwAddr()]),
fPadMapHw(new Float_t[GetMaxHwAddr()]),
fPadMapRCU(new Int_t*[GetMaxHwAddr()]),
fHistIROC(0),
fHistOROC(0),
fHistIROCIndex(0),
fHistOROCIndex(0),
fHistIROCTime(0),
fHistOROCTime(0),
fHistIROCClone(0),
fHistOROCClone(0),
fHistIROCRMS(0),
fHistOROCRMS(0),
fHistIROCBASE(0),
fHistOROCBASE(0),
fHistIROCSUM(0),
fHistOROCSUM(0),
fHistChannelTime(0),
fHistAddrMapIndex(0),
fHistAddrMaxAdc(0),
fHistAddrBaseMean(0),
fHistAddrMaxAdcX(0),
fHistAddrAdcSum(0),
fHistAddrBaseRms(0),
fHistDistrSumIROC(0),
fHistDistrMaxIROC(0),
fHistDistrSumOROC(0),
fHistDistrMaxOROC(0),
fHistDistrBase2dIROC(0),
fHistDistrBase2dOROC(0),
fHistDistrBaseRmsIROC(0),
fHistDistrBaseMeanIROC(0),
fHistDistrBaseRmsOROC(0),
fHistDistrBaseMeanOROC(0),
fHistGlobalMaxA(0),
fHistGlobalMaxC(0),
fHistList(new TObjArray()),
fkNRowsIroc(monitor.fkNRowsIroc),
fkNRowsOroc(monitor.fkNRowsOroc),
fkNPadsIroc(monitor.fkNPadsIroc),
fkNPadsOroc(monitor.fkNPadsOroc),
fkNPadMinIroc(monitor.fkNPadMinIroc),
fkNPadMinOroc(monitor.fkNPadMinOroc),
fkNPadMaxIroc(monitor.fkNPadMaxIroc),
fkNPadMaxOroc(monitor.fkNPadMaxOroc),
fVerb(monitor.fVerb),
fLastEv(monitor.fLastEv),
fEventNumber(monitor.fEventNumber),
fEventNumberOld(monitor.fEventNumberOld),
fDisableFit(monitor.fDisableFit),
fExecGlob(monitor.fExecGlob),
fExecPlaneMax(monitor.fExecPlaneMax),
fExecPadIrocRms(monitor.fExecPadIrocRms),
fExecPadOrocRms(monitor.fExecPadOrocRms),
fRunId(monitor.fRunId),
fEqId(monitor.fEqId),
fPadUsedRoc(monitor.fPadUsedRoc),
fPadUsedHwAddr(monitor.fPadUsedHwAddr),
fGdcId(monitor.fGdcId),
fLdcId(monitor.fLdcId),
fLdcIdOld(monitor.fLdcIdOld),
fMapEqidsSec(new Int_t*[36]),
fMapEqidsRcu(new Int_t[1000]),
fMirror(monitor.fMirror),
fChannelIter(monitor.fChannelIter),
fMapHand(monitor.fMapHand),
fRawReader(monitor.fRawReader),
fkMonTable(0x0),
fMonTableString(""),
fMonTableArray(0x0),
fMonTableChanged(kFALSE)
{
  // copy constructor
  
  fHistIROC=(TH2F*)monitor.fHistIROC->Clone(); fHistList->Add(fHistIROC);
  fHistOROC=(TH2F*)monitor.fHistOROC->Clone(); fHistList->Add(fHistOROC);
  fHistIROCIndex=(TH2S*)monitor.fHistIROCIndex->Clone(); fHistList->Add(fHistIROCIndex);
  fHistOROCIndex=(TH2S*)monitor.fHistOROCIndex->Clone(); fHistList->Add(fHistOROCIndex);
  fHistIROCTime=(TH2F*)monitor.fHistIROCTime->Clone(); fHistList->Add(fHistIROCTime);
  fHistOROCTime=(TH2F*)monitor.fHistOROCTime->Clone(); fHistList->Add(fHistOROCTime);
  fHistIROCClone=(TH2F*)monitor.fHistIROCClone->Clone(); fHistList->Add(fHistIROCClone);
  fHistOROCClone=(TH2F*)monitor.fHistOROCClone->Clone(); fHistList->Add(fHistOROCClone);
  fHistIROCRMS=(TH2F*)monitor.fHistIROCRMS->Clone(); fHistList->Add(fHistIROCRMS);
  fHistOROCRMS=(TH2F*)monitor.fHistOROCRMS->Clone(); fHistList->Add(fHistOROCRMS);
  fHistIROCBASE=(TH2F*)monitor.fHistIROCBASE->Clone(); fHistList->Add(fHistIROCBASE);
  fHistOROCBASE=(TH2F*)monitor.fHistOROCBASE->Clone(); fHistList->Add(fHistOROCBASE);
  fHistIROCSUM=(TH2F*)monitor.fHistIROCSUM->Clone(); fHistList->Add(fHistIROCSUM);
  fHistOROCSUM=(TH2F*)monitor.fHistOROCSUM->Clone(); fHistList->Add(fHistOROCSUM);
  
  fHistChannelTime=(TH2F*)monitor.fHistChannelTime->Clone(); fHistList->Add(fHistChannelTime);
  
  fHistAddrMapIndex=(TH1F*)monitor.fHistAddrMapIndex->Clone(); fHistList->Add(fHistAddrMapIndex);
  fHistAddrMaxAdc=(TH1F*)monitor.fHistAddrMaxAdc->Clone(); fHistList->Add(fHistAddrMaxAdc);
  fHistAddrBaseMean=(TH1F*)monitor.fHistAddrBaseMean->Clone(); fHistList->Add(fHistAddrBaseMean);
  fHistAddrMaxAdcX=(TH1F*)monitor.fHistAddrMaxAdcX->Clone(); fHistList->Add(fHistAddrMaxAdcX);
  fHistAddrAdcSum=(TH1F*)monitor.fHistAddrAdcSum->Clone(); fHistList->Add(fHistAddrAdcSum);
  fHistAddrBaseRms=(TH1F*)monitor.fHistAddrBaseRms->Clone(); fHistList->Add(fHistAddrBaseRms);
  
  
  fHistDistrSumIROC=(TH1F*)monitor.fHistDistrSumIROC->Clone(); fHistList->Add(fHistDistrSumIROC);
  fHistDistrMaxIROC=(TH1F*)monitor.fHistDistrMaxIROC->Clone(); fHistList->Add(fHistDistrMaxIROC);
  fHistDistrSumOROC=(TH1F*)monitor.fHistDistrSumOROC->Clone(); fHistList->Add(fHistDistrSumOROC);
  fHistDistrMaxOROC=(TH1F*)monitor.fHistDistrMaxOROC->Clone(); fHistList->Add(fHistDistrMaxOROC);
  
  fHistDistrBase2dIROC=(TH2F*)monitor.fHistDistrBase2dIROC->Clone(); fHistList->Add(fHistDistrBase2dIROC);
  fHistDistrBase2dOROC=(TH2F*)monitor.fHistDistrBase2dOROC->Clone(); fHistList->Add(fHistDistrBase2dOROC);
  fHistDistrBaseRmsIROC=(TH1D*)monitor.fHistDistrBaseRmsIROC->Clone(); fHistList->Add(fHistDistrBaseRmsIROC);
  fHistDistrBaseMeanIROC=(TH1D*)monitor.fHistDistrBaseMeanIROC->Clone(); fHistList->Add(fHistDistrBaseMeanIROC);
  fHistDistrBaseRmsOROC=(TH1D*)monitor.fHistDistrBaseRmsOROC->Clone(); fHistList->Add(fHistDistrBaseRmsOROC);
  fHistDistrBaseMeanOROC=(TH1D*)monitor.fHistDistrBaseMeanOROC->Clone(); fHistList->Add(fHistDistrBaseMeanOROC);
  
  fHistGlobalMaxA=(TH2S*)monitor.fHistGlobalMaxA->Clone(); fHistList->Add(fHistGlobalMaxA);
  fHistGlobalMaxC=(TH2S*)monitor.fHistGlobalMaxC->Clone(); fHistList->Add(fHistGlobalMaxC);
  
  // fPad       = new Int_t*[monitor.GetMaxHwAddr()];
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPad[i]       = new Int_t[monitor.GetTimeBins()];}
  
  //fPadMapRCU = new Int_t*[monitor.GetMaxHwAddr()];
  for(Int_t i = 0; i<monitor.GetMaxHwAddr(); i++) { fPadMapRCU[i] = new Int_t[6];}
  
  //fPadMapHw  = new Float_t[monitor.GetMaxHwAddr()];
  
  //fMapEqidsRcu   = new Int_t[1000];
  //fMapEqidsSec   = new Int_t*[36];
  for(Int_t i = 0; i<36; i++) { fMapEqidsSec[i]       = new Int_t[6];}
  SetEqIds();
  
  if (gDirectory)
  {
    if (!gDirectory->GetList())
    {
      Warning("Build","Current directory is not a valid directory");
      
    }
    else
    {
      AliTPCMonitor *hold = (AliTPCMonitor*)gDirectory->GetList()->FindObject(monitor.GetName());
      if(hold)
      {
        Warning("Build","Replacing existing histogram: %s (Potential memory leak).",monitor.GetName());
        gDirectory->GetList()->Remove(hold);
      }
      gDirectory->Append(this);
    }
  }
}


//____________________________________________________________________________

AliTPCMonitor &AliTPCMonitor:: operator= (const AliTPCMonitor& monitor)
{
  // assigment operator
  if(this!=&monitor)
  {
    ((AliTPCMonitorConfig *)this)->operator=(monitor);
    fkNRowsIroc=monitor.fkNRowsIroc;
    fkNRowsOroc=monitor.fkNRowsOroc;
    fkNPadsIroc=monitor.fkNPadsIroc;
    fkNPadsOroc=monitor.fkNPadsOroc;
    fkNPadMinIroc=monitor.fkNPadMinIroc;
    fkNPadMinOroc=monitor.fkNPadMinOroc;
    fkNPadMaxIroc=monitor.fkNPadMaxIroc;
    fkNPadMaxOroc=monitor.fkNPadMaxOroc;
    fVerb=monitor.fVerb;
    fLastEv=monitor.fLastEv;
    fEventNumber=monitor.fEventNumber;
    fEventNumberOld=monitor.fEventNumberOld;
    fDisableFit=monitor.fDisableFit;
    fExecGlob=monitor.fExecGlob;
    fExecPlaneMax=monitor.fExecPlaneMax;
    fExecPadIrocRms=monitor.fExecPadIrocRms;
    fExecPadOrocRms=monitor.fExecPadOrocRms;
    fRunId=monitor.fRunId;
    fEqId=monitor.fEqId;
    fPadUsedRoc=monitor.fPadUsedRoc;
    fPadUsedHwAddr=monitor.fPadUsedHwAddr;
    fGdcId=monitor.fGdcId;
    fLdcId=monitor.fLdcId;
    fLdcIdOld=monitor.fLdcIdOld;
    fMapHand=monitor.fMapHand;
    fRawReader=monitor.fRawReader;
    fkMonTable=0x0;
    fMonTableString="";
    fMonTableArray=0x0;
    fMonTableChanged=kFALSE;
    
    fHistList = new TObjArray();
    fHistIROC=(TH2F*)monitor.fHistIROC->Clone(); fHistList->Add(fHistIROC);
    fHistOROC=(TH2F*)monitor.fHistOROC->Clone(); fHistList->Add(fHistOROC);
    fHistIROCIndex=(TH2S*)monitor.fHistIROCIndex->Clone(); fHistList->Add(fHistIROCIndex);
    fHistOROCIndex=(TH2S*)monitor.fHistOROCIndex->Clone(); fHistList->Add(fHistOROCIndex);
    fHistIROCTime=(TH2F*)monitor.fHistIROCTime->Clone(); fHistList->Add(fHistIROCTime);
    fHistOROCTime=(TH2F*)monitor.fHistOROCTime->Clone(); fHistList->Add(fHistOROCTime);
    fHistIROCClone=(TH2F*)monitor.fHistIROCClone->Clone(); fHistList->Add(fHistIROCClone);
    fHistOROCClone=(TH2F*)monitor.fHistOROCClone->Clone(); fHistList->Add(fHistOROCClone);
    fHistIROCRMS=(TH2F*)monitor.fHistIROCRMS->Clone(); fHistList->Add(fHistIROCRMS);
    fHistOROCRMS=(TH2F*)monitor.fHistOROCRMS->Clone(); fHistList->Add(fHistOROCRMS);
    fHistIROCBASE=(TH2F*)monitor.fHistIROCBASE->Clone(); fHistList->Add(fHistIROCBASE);
    fHistOROCBASE=(TH2F*)monitor.fHistOROCBASE->Clone(); fHistList->Add(fHistOROCBASE);
    fHistIROCSUM=(TH2F*)monitor.fHistIROCSUM->Clone(); fHistList->Add(fHistIROCSUM);
    fHistOROCSUM=(TH2F*)monitor.fHistOROCSUM->Clone(); fHistList->Add(fHistOROCSUM);
    
    fHistChannelTime=(TH2F*)monitor.fHistChannelTime->Clone(); fHistList->Add(fHistChannelTime);
    
    fHistAddrMapIndex=(TH1F*)monitor.fHistAddrMapIndex->Clone(); fHistList->Add(fHistAddrMapIndex);
    fHistAddrMaxAdc=(TH1F*)monitor.fHistAddrMaxAdc->Clone(); fHistList->Add(fHistAddrMaxAdc);
    fHistAddrBaseMean=(TH1F*)monitor.fHistAddrBaseMean->Clone(); fHistList->Add(fHistAddrBaseMean);
    fHistAddrMaxAdcX=(TH1F*)monitor.fHistAddrMaxAdcX->Clone(); fHistList->Add(fHistAddrMaxAdcX);
    fHistAddrAdcSum=(TH1F*)monitor.fHistAddrAdcSum->Clone(); fHistList->Add(fHistAddrAdcSum);
    fHistAddrBaseRms=(TH1F*)monitor.fHistAddrBaseRms->Clone(); fHistList->Add(fHistAddrBaseRms);
    
    
    fHistDistrSumIROC=(TH1F*)monitor.fHistDistrSumIROC->Clone(); fHistList->Add(fHistDistrSumIROC);
    fHistDistrMaxIROC=(TH1F*)monitor.fHistDistrMaxIROC->Clone(); fHistList->Add(fHistDistrMaxIROC);
    fHistDistrSumOROC=(TH1F*)monitor.fHistDistrSumOROC->Clone(); fHistList->Add(fHistDistrSumOROC);
    fHistDistrMaxOROC=(TH1F*)monitor.fHistDistrMaxOROC->Clone(); fHistList->Add(fHistDistrMaxOROC);
    
    fHistDistrBase2dIROC=(TH2F*)monitor.fHistDistrBase2dIROC->Clone(); fHistList->Add(fHistDistrBase2dIROC);
    fHistDistrBase2dOROC=(TH2F*)monitor.fHistDistrBase2dOROC->Clone(); fHistList->Add(fHistDistrBase2dOROC);
    fHistDistrBaseRmsIROC=(TH1D*)monitor.fHistDistrBaseRmsIROC->Clone(); fHistList->Add(fHistDistrBaseRmsIROC);
    fHistDistrBaseMeanIROC=(TH1D*)monitor.fHistDistrBaseMeanIROC->Clone(); fHistList->Add(fHistDistrBaseMeanIROC);
    fHistDistrBaseRmsOROC=(TH1D*)monitor.fHistDistrBaseRmsOROC->Clone(); fHistList->Add(fHistDistrBaseRmsOROC);
    fHistDistrBaseMeanOROC=(TH1D*)monitor.fHistDistrBaseMeanOROC->Clone(); fHistList->Add(fHistDistrBaseMeanOROC);
    
    fHistGlobalMaxA=(TH2S*)monitor.fHistGlobalMaxA->Clone(); fHistList->Add(fHistGlobalMaxA);
    fHistGlobalMaxC=(TH2S*)monitor.fHistGlobalMaxC->Clone(); fHistList->Add(fHistGlobalMaxC);
    
    fPad       = new Int_t*[monitor.GetMaxHwAddr()];
    for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPad[i]       = new Int_t[monitor.GetTimeBins()];}
    
    fPadMapRCU = new Int_t*[monitor.GetMaxHwAddr()];
    for(Int_t i = 0; i<monitor.GetMaxHwAddr(); i++) { fPadMapRCU[i] = new Int_t[6];}
    
    fPadMapHw  = new Float_t[monitor.GetMaxHwAddr()];
    
    fMapEqidsRcu   = new Int_t[1000];
    fMapEqidsSec   = new Int_t*[36];
    for(Int_t i = 0; i<36; i++) { fMapEqidsSec[i]       = new Int_t[6];}
    SetEqIds();
    
    if (gDirectory)
    {
      if (!gDirectory->GetList())
      {
        Warning("Build","Current directory is not a valid directory");
        return  *this;
      }
      AliTPCMonitor *hold = (AliTPCMonitor*)gDirectory->GetList()->FindObject(monitor.GetName());
      if(hold)
      {
        Warning("Build","Replacing existing histogram: %s (Potential memory leak).",monitor.GetName());
        gDirectory->GetList()->Remove(hold);
      }
      gDirectory->Append(this);
    }
  }
  return *this;
}


//____________________________________________________________________________
AliTPCMonitor::~AliTPCMonitor() 
{
  // Destructor
  
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { delete[] fPad[i] ;}
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { delete[] fPadMapRCU[i];}
  delete[] fPadMapHw ;
  DeleteHistos();
}

//____________________________________________________________________________
void AliTPCMonitor::CreateHistos() 
{
  // Create histograms to be displayed
  
  if(fVerb) cout << " create new ones " << endl;
  fHistIROC            = new TH2F("fHistIROC"            ,"fHistIROC"           ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROC);
  fHistOROC            = new TH2F("fHistOROC"            ,"fHistOROC"           ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROC);
  
  fHistIROCIndex       = new TH2S("fHistIROCIndex"       ,"fHistIROCIndex"      ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROCIndex);
  fHistOROCIndex       = new TH2S("fHistOROCIndex"       ,"fHistOROCIndex"      ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROCIndex);
  
  fHistIROCTime        = new TH2F("fHistIROCTime"        ,"fHistIROCTime"       ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROCTime);
  fHistOROCTime        = new TH2F("fHistOROCTime"        ,"fHistOROCTime"       ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROCTime);
  
  fHistIROCRMS         = new TH2F("fHistIROCRMS"         ,"fHistIROCRMS"        ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROCRMS);
  fHistOROCRMS         = new TH2F("fHistOROCRMS"         ,"fHistOROCRMS"        ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROCRMS);
  
  fHistIROCSUM         = new TH2F("fHistIROCSUM"         ,"fHistIROCSUM"        ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROCSUM);
  fHistOROCSUM         = new TH2F("fHistOROCSUM"         ,"fHistOROCSUM"        ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROCSUM);
  
  fHistIROCBASE        = new TH2F("fHistIROCBASE"        ,"fHistIROCBASE"       ,fkNRowsIroc,0,fkNRowsIroc,fkNPadsIroc, fkNPadMinIroc, fkNPadMaxIroc);       fHistList->Add(fHistIROCBASE);
  fHistOROCBASE        = new TH2F("fHistOROCBASE"        ,"fHistOROCBASE"       ,fkNRowsOroc,0,fkNRowsOroc,fkNPadsOroc, fkNPadMinOroc, fkNPadMaxOroc);       fHistList->Add(fHistOROCBASE);
  
  
  fHistChannelTime     = new TH2F("fHistChannelTime"     ,"fHistChannelTime"    ,GetNumOfChannels(),0,GetNumOfChannels(),GetTimeBins(),0,GetTimeBins());fHistList->Add(fHistChannelTime);
  fHistAddrMapIndex    = new TH1F("fHistAddrMapIndex"    ,"fHistAddrMapIndex"   ,GetMaxHwAddr()    ,0,GetMaxHwAddr());                                  fHistList->Add(fHistAddrMapIndex);
  fHistAddrMaxAdc      = new TH1F("fHistAddrMaxAdc"      ,"fHistAddrMaxAdc"     ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrMaxAdc);
  fHistAddrMaxAdcX     = new TH1F("fHistAddrMaxAdcX"     ,"fHistAddrMaxAdcX"    ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrMaxAdcX);
  fHistAddrBaseMean    = new TH1F("fHistAddrBaseMean"    ,"fHistAddrBaseMean"   ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrBaseMean);
  fHistAddrAdcSum      = new TH1F("fHistAddrAdcSum"      ,"fHistAddrAdcSum"     ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrAdcSum);
  fHistAddrBaseRms     = new TH1F("fHistAddrBaseRms"     ,"fHistAddrBaseRms"    ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrBaseRms);
  fHistDistrSumIROC    = new TH1F("fHistDistrSumIROC"    ,"fHistDistrSumIROC"   ,400,0.0,4000.0);                                                       fHistList->Add(fHistDistrSumIROC);
  fHistDistrMaxIROC    = new TH1F("fHistDistrMaxIROC"    ,"fHistDistrMaxIROC"   ,500,0.0,1000.0);                                                       fHistList->Add(fHistDistrMaxIROC);
  fHistDistrSumOROC    = new TH1F("fHistDistrSumOROC"    ,"fHistDistrSumOROC"   ,400,0.0,4000.0);                                                       fHistList->Add(fHistDistrSumOROC);
  fHistDistrMaxOROC    = new TH1F("fHistDistrMaxOROC"    ,"fHistDistrMaxOROC"   ,500,0.0,1000.0);                                                       fHistList->Add(fHistDistrMaxOROC);
  
  fHistDistrBase2dIROC = new TH2F("fHistDistrBase2dIROC" ,"fHistDistrBase2dIROC",100,0.0,100.0,100,0.0,10.0);                                           fHistList->Add(fHistDistrBase2dIROC);
  fHistDistrBase2dOROC = new TH2F("fHistDistrBase2dOROC" ,"fHistDistrBase2dOROC",100,0.0,100.0,100,0.0,10.0);                                           fHistList->Add(fHistDistrBase2dOROC);
  
  fHistGlobalMaxA      = new TH2S("SIDE A"               ,"SIDE A"              ,500,-3000,3000,500,-3000,3000);                                        fHistList->Add(fHistGlobalMaxA);
  fHistGlobalMaxC      = new TH2S("SIDE C"               ,"SIDE C"              ,500,-3000,3000,500,-3000,3000);                                        fHistList->Add(fHistGlobalMaxC);
  
  ResetArrays();
}
//____________________________________________________________________________
Int_t AliTPCMonitor::ProcessEvent()
{
  // Process Event
  // Depending on the value of the sector id all sectors (sectorid == -1) are processed.
  //
  // In this case only the maximum values are calculated per pad and filled to the global histograms
  // In a second loop the last processed(displayed) sector will be processed (sectorid!=-1)
  // again and the baseline rms and further quantities are calculated
  //
  // If only one sector should be processed SetProcOneSector(1) should be set.
  // In this case only the specified (last/last displayed) sector will be processed.
  //
  // If GetProcNextEvent()==0 the same event will be processed again
  
  
  Int_t sectorid = 0;
  Int_t retflag  = 0; // id of last sector + 1000,  or error flag
  if(GetProcNextEvent()==1 && fLastEv) { AliInfo("Last event already processed"); }
  if(GetProcNextEvent()==1) ResetSectorArray();
  
  
  if(GetProcNextEvent()==0 || GetProcOneSector()==1 ) sectorid = GetLastSector();
  else                                                sectorid = -1;
  
  // first iteration
  retflag = ReadDataNew(sectorid);
  
  SetLastProcFile(GetFile());
  
  if(retflag>=10 && retflag<1000){ AliError("Could not read event properly: Check file name and format or try next event");  return 0  ;}
  
  DrawHists(3);
  
  // second iteration
  if(sectorid==-1 && retflag >1000)
  {
    AliInfo("Second read cycle");
    SetProcNextEvent(0);
    if(GetLastSectorDisplayed()==-1) {sectorid =  GetLastSector()         ;    }
    else                             {sectorid =  GetLastSectorDisplayed();  SetLastSector(sectorid)           ;  }
    retflag = ReadDataNew(sectorid);
  }
  
  SetLastSectorDisplayed(sectorid)  ;
  fMapHand->ReadfecHwMap(GetLastSector());
  FillHistsPadPlane();
  DrawHists(1);
  SetEventProcessed(1);
  return retflag;
}
//__________________________________________________________________
Int_t AliTPCMonitor::ReadDataNew(Int_t secid)
{
  // Read Data File/Stream  for specified Format.
  // Payload will be extracted from either ROOT or DATE format
  // and passed to FillHistsDecode for decoding of the adc information

  if (TString(GetLastProcFile())!=GetFile() || fMonTableChanged) {
    delete fRawReader;
    fRawReader=0x0;
  }
  
  if (!fRawReader){
    TString file(GetFile());
    if (file.BeginsWith("mem://")){
      file.ReplaceAll("mem://","");
      fRawReader = new AliRawReaderDateOnline(file.Data(),fkMonTable);
    }else{
      fRawReader = AliRawReader::Create(GetFile());
    }
    SetLastProcFile(GetFile());
  }
  
  if (!fRawReader){
    AliWarning("Coult not initialize raw reader");
    return 11;
  }
  
  Bool_t skip=kTRUE;
  while(skip && GetProcNextEvent())
  {
    if(fVerb) cout << "AliTPCMonitor::ReadDataNew get event " << endl;
    if(fRawReader->IsA()==AliRawReaderRoot::Class()){
      AliInfo(Form("Root, NextEvent: %d\n",GetEventID()));
      if (!fRawReader->GotoEvent(GetEventID())){AliError("Could not get next Event"); return 11 ;}
    } else if(!fRawReader->NextEvent()) { AliError("Could not get next Event"); return 11 ;}
    // skip all events but physics, calibration and software trigger events!
    UInt_t eventType=fRawReader->GetType();
    if ( !(eventType==AliRawEventHeaderBase::kPhysicsEvent ||
           eventType==AliRawEventHeaderBase::kCalibrationEvent ||
           eventType==AliRawEventHeaderBase::kSystemSoftwareTriggerEvent ||
           eventType==AliRawEventHeaderBase::kDetectorSoftwareTriggerEvent) ) {
             if (fVerb) cout<< "Skipping event! Its neither of 'physics, calibration and software trigger event'" << endl;
             if(fRawReader->IsA()==AliRawReaderRoot::Class()){
               if (fEventNumber<fRawReader->GetNumberOfEvents()-1) ++fEventNumber;
               else {AliError("No more events");return 11;}
             }
             continue;
           }
    skip=kFALSE;
    //test if the TPC has data
    UChar_t *data=0;
    fRawReader->Select("TPC");
    Int_t eventNr=fRawReader->GetEventIndex();
    if (!fRawReader->ReadNextData(data)) {
      skip=kTRUE;
      AliInfo(Form("%d / %d",fEventNumber,fRawReader->GetNumberOfEvents()));
      if(fRawReader->IsA()==AliRawReaderRoot::Class()){
        if (fEventNumber<fRawReader->GetNumberOfEvents()-1){
          ++eventNr;
//           printf("inc conter\n");
        }
        else {
          AliError("No more events");
          return 11;
        }
      }
    }
    fEventNumber = eventNr;
    fEventNumberOld = eventNr;
  }
  
  AliInfo(Form("secid: %d",secid));
  
  //========================== Histogram filling ======================
  ResetHistos() ;
  AliAltroRawStreamV3 *altro=new AliAltroRawStreamV3(fRawReader);
  altro->SelectRawData("TPC");
  altro->Reset();
  fChannelIter=0;
  
  Int_t     hw                = 0;
  Int_t     nextHwAddress     = 0;
  Int_t     rcupatch          = 0;
  Int_t     maxADC            = 0;
  Int_t     maxx              = 0;
  Int_t     sum               = 0;
  Int_t     sumn              = 0;
  Int_t           lastrcuid   = 0;
  
  while ( altro->NextDDL() ){
    fGdcId        = fRawReader->GetGDCId() ;
    fLdcId        = fRawReader->GetLDCId() ;
    fRunId        = fRawReader->GetRunNumber() ;
    fEqId         = fRawReader->GetEquipmentId();
    rcupatch      = GetRCUPatch(fRunId, fEqId);
    Int_t rcupatchSector=rcupatch%6;
    lastrcuid     = (rcupatch+1000);

    Int_t currentSector=rcupatch/6;
    
    if(fLdcIdOld!=fLdcId &&  fChannelIter!=0) {
      if(secid==-1)
      {
        FillGlobal(GetLastSector());
        ResetArrays();
        if (fVerb) printf("filled sector: %d\n",GetLastSector());
        fChannelIter =0;
      }
//       else
//       {
//         printf("RET: filled sector: %d\n",GetLastSector());
//         return lastrcuid;
//         if (rcupatch/6!=secid) continue;
//       }
//       fChannelIter=0;
    }
    
    if (!CheckEqId(secid,fEqId)) continue;
    if (fVerb) printf("Sector: %d, RCU patch: %d, LDC: %d, EqId: %d\n",currentSector,rcupatch, fLdcId, fEqId);
    
    while ( altro->NextChannel() ){
      hw=altro->GetHWAddress();
      nextHwAddress         = ( hw + (rcupatchSector<<12) );
      fPad[fChannelIter][0] = nextHwAddress ;
      fPadMapHw[nextHwAddress]    =  fChannelIter  ;
      maxADC=0;
      while ( altro->NextBunch() ){
        Int_t  startTbin    = (Int_t)altro->GetStartTimeBin();
        Int_t  bunchlength  = (Int_t)altro->GetBunchLength();
        const UShort_t *sig = altro->GetSignals();
        for (Int_t iTimeBin = 0; iTimeBin<bunchlength; iTimeBin++){
          Int_t adc=(Int_t)sig[iTimeBin];
          Int_t ntime=startTbin-iTimeBin;
          //fill channel information
          fPad[fChannelIter][ntime]  = adc;
          if( (adc>maxADC) && (ntime>=GetRangeMaxAdcMin()) && (ntime<GetRangeMaxAdcMax()  ))  {maxADC = adc;maxx = ntime ;}
          if(              (ntime>=GetRangeSumMin())    && (ntime<GetRangeSumMax()     ))  {sum+=adc; sumn++;}
          
        }
      }
      //get pedestal, noise
      Float_t pedestal=TMath::Mean(GetRangeBaseMax()-GetRangeBaseMin(),fPad[fChannelIter]+GetRangeBaseMin());
      Float_t noise   =TMath::RMS(GetRangeBaseMax()-GetRangeBaseMin(),fPad[fChannelIter]+GetRangeBaseMin());
      fHistAddrMaxAdc->SetBinContent(nextHwAddress,maxADC-GetPedestals()*pedestal);
      
      if(secid!=-1)
      {
        if(rcupatchSector<2)
        {
          fHistDistrBase2dIROC->Fill(pedestal,noise);
          fHistDistrSumIROC->Fill(sum);
          fHistDistrMaxIROC->Fill(maxADC-pedestal*GetPedestals());
          fHistDistrSumIROC->Fill(sum -sumn*pedestal*GetPedestals());
        }
        else
        {
          fHistDistrBase2dOROC->Fill(pedestal,noise);
          fHistDistrSumOROC->Fill(sum);
          fHistDistrMaxOROC->Fill(maxADC-pedestal*GetPedestals());
          fHistDistrSumOROC->Fill(sum -sumn*pedestal*GetPedestals());
        }
        
        fHistAddrAdcSum->SetBinContent(  nextHwAddress,sum);
        fHistAddrMapIndex->SetBinContent(nextHwAddress,fChannelIter);
        fHistAddrBaseMean->SetBinContent(nextHwAddress,pedestal);
        fHistAddrMaxAdcX->SetBinContent( nextHwAddress,maxx);
        fHistAddrBaseRms->SetBinContent( nextHwAddress,noise);
      }
      
      ++fChannelIter;
    }
    SetLastSector(rcupatch/6);
    if(fChannelIter!=0 && secid==-1 ) SetSectorFilled(rcupatch/6);
    fLdcIdOld = fLdcId ;
  }
  delete altro;
  if(fChannelIter!=0 && secid==-1) { FillGlobal(GetLastSector());}
  return lastrcuid;
}











//____________________________________________________________________________
void AliTPCMonitor::FillHistsPadPlane() 
{
  // Fill 2Dim histograms for IROC and OROC (max , rms and sum)
  
  if(fVerb)cout << "AliTPCMonitor::FillHistsPadPlane() Start  " << endl;
  if(fVerb)PrintConfig();
  
  Int_t pad    = 0;
  Int_t row    = 0;
  Int_t padmax = 0;
  Int_t hwadd  = 0;
  
  for(Int_t ch = 0; ch<fChannelIter; ch++)
  {
    hwadd=  fPad[ch][0];
    fHistChannelTime->SetCellContent(ch,0,hwadd);
    
    for(Int_t bin = 1; bin <GetTimeBins(); bin++)
    {
      if( fHistChannelTime->GetCellContent(ch,bin)!=0)  cout << " cellcontent already set " << endl;
      if(     GetPedestals()==1 ) fHistChannelTime->SetCellContent(ch,bin,(fPad[ch][bin]- fHistAddrBaseMean->GetBinContent(hwadd)));
      else                        fHistChannelTime->SetCellContent(ch,bin,(fPad[ch][bin]));
    }
    
    pad    = fMapHand->GetPad(   hwadd);
    row    = fMapHand->GetPadRow(hwadd);
    padmax = fMapHand->GetNumofPads(row);
    
    if(row<63)
    {
      fHistIROC->SetCellContent(     row    +1 ,pad +55 -padmax/2 +1,fHistAddrMaxAdc->GetBinContent(  hwadd));
      fHistIROCIndex->SetCellContent(row    +1 ,pad +55 -padmax/2 +1,ch);
      fHistIROCRMS->SetCellContent(  row    +1 ,pad +55 -padmax/2 +1,fHistAddrBaseRms->GetBinContent( hwadd));
      fHistIROCBASE->SetCellContent( row    +1 ,pad +55 -padmax/2 +1,fHistAddrBaseMean->GetBinContent(hwadd));
      fHistIROCSUM->SetCellContent(  row    +1 ,pad +55 -padmax/2 +1,fHistAddrAdcSum->GetBinContent(  hwadd));
    }
    else
    {
      fHistOROC->SetCellContent(     row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrMaxAdc->GetBinContent(  hwadd));
      fHistOROCIndex->SetCellContent(row-63 +1 ,pad +70 -padmax/2 +1,ch);
      fHistOROCRMS->SetCellContent(  row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrBaseRms->GetBinContent( hwadd));
      fHistOROCBASE->SetCellContent( row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrBaseMean->GetBinContent(hwadd));
      fHistOROCSUM->SetCellContent(  row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrAdcSum->GetBinContent(  hwadd));
    }
  }
  
  fHistChannelTime->GetXaxis()->SetRange(0,fChannelIter);
  fHistChannelTime->GetYaxis()->SetRange(0,GetTimeBins());
}



//____________________________________________________________________________
void AliTPCMonitor::ResetArrays() 
{
  // Reset data arrays
  for(Int_t row = 0 ; row < fkNRowsIroc; row++)
  {
    for(Int_t pad = 0 ; pad <  fkNPadsIroc ; pad++)
    {
      fHistIROCIndex->SetCellContent(row+1,pad+1,-1);
    }
  }
  for(Int_t row = 0 ; row < fkNRowsOroc; row++)
  {
    for(Int_t pad = 0 ; pad <  fkNPadsOroc ; pad++)
    {
      fHistOROCIndex->SetCellContent(row+1,pad+1,-1);
    }
  }
  
  for(Int_t ch= 0; ch<GetMaxHwAddr(); ch++)
  {
    fHistAddrMaxAdcX->SetBinContent(ch,-1);
    fHistAddrMapIndex->SetBinContent(ch,-1);
    fHistAddrMaxAdc->SetBinContent(  ch, 0);
    fHistAddrBaseMean->SetBinContent(  ch, 0);
    fHistAddrAdcSum->SetBinContent(  ch, 0);
    fHistAddrBaseRms->SetBinContent(ch, 0);
  }
  
  for(Int_t ch = 0; ch< GetNumOfChannels(); ch++)
  {
    for(Int_t bin = 0; bin< GetTimeBins(); bin++)
    {
      fPad[ch][bin] = 0;
    }
  }
  for(Int_t ch = 0; ch< GetMaxHwAddr(); ch++)
  {
    fPadMapHw[ch]=-1;
    fPadMapRCU[ch][0]=-1;
    fPadMapRCU[ch][1]=-1;
    fPadMapRCU[ch][2]=-1;
    fPadMapRCU[ch][3]=-1;
    fPadMapRCU[ch][4]=-1;
    fPadMapRCU[ch][5]=-1;
  }
  fPadUsedHwAddr=-1;
}


//____________________________________________________________________________
void AliTPCMonitor::ResetHistos() 
{
  // Reset all but
  for(Int_t i =0; i<fHistList->GetEntries(); i++)
  {
    TString name=fHistList->At(i)->GetName();
    if(GetProcNextEvent()==0 && (name=="SIDE A" || name=="SIDE C")) continue;
    ((TH1*)fHistList->At(i))->Reset();
  }
  ResetArrays();
}

//____________________________________________________________________________
void AliTPCMonitor::DeleteHistos() 
{
  // Delete histograms
  for(Int_t i =0; i<fHistList->GetEntries(); i++)
  {
    delete (TH1*)fHistList->At(i);
  }
}


//__________________________________________________________________
Int_t AliTPCMonitor::CheckEqId(Int_t secid,Int_t eqid)
{
  // Check if equipment id corresponds to any  rcu patch in sector
  // Equipment ids changed during commisioning in 2006 (starting from run 704)
  // However Runids started from 0 again in 2007
  // Now only runids from commissioning in 2006 after runid 704 and all new ones are supported.
  // Comment in equipment check for runids < 704 if old runs should be processed
  
  if(fVerb) cout << "AliTPCMonitor::CheckEqId  : SectorId  " << secid << " EquipmentId " << eqid  << " runid " << fRunId <<  endl;
  Int_t passed =1;
  //skip all eqids which do not belong to the TPC
  if ( eqid<768||eqid>983 ) return 0;
  //
//   if(fRunId<704 && 0) // commented out --> runs with runid < 704 in 2006 are not recognized anymore
//   {
//     if( (secid>-1) && (secid<36) )   // if ( secid is in range) { take only specific eqids}  else { take all }
//     {
//       if(      (secid==13) && ( eqid!=408 && eqid!=409  &&  eqid!=509 && eqid!=512  && eqid!=513 && eqid!=517 )) {passed=0;}
//       else if( (secid==4)  && ( eqid!=404 && eqid!=504  &&  eqid!=407 && eqid!=503  && eqid!=508 && eqid!=506 )) {passed=0;}
//     }
//     else                                                                   {if(fVerb) cout << "passed check "<< endl; }
//   }
//   else
//   {
    if( (secid>-1) && (secid<36) )   // if ( secid is in range) { take only specific eqids}  else { take all }
    {
      if(eqid!=fMapEqidsSec[secid][0] && eqid!= fMapEqidsSec[secid][1] && eqid!=fMapEqidsSec[secid][2] &&
         eqid!=fMapEqidsSec[secid][3] && eqid!= fMapEqidsSec[secid][4] && eqid!=fMapEqidsSec[secid][5] )  {passed=0;}
    }
    else                                                                   {if(fVerb) cout << "passed check "<< endl;}
//   }
  
  return passed;
}

//__________________________________________________________________
void AliTPCMonitor::SetEqIds()
{
  // Set mapping for equipment ids
  for(Int_t i = 0; i<36 ; i++)
  {
    for(Int_t j = 0; j<6; j++)
    {
      if(j<2) fMapEqidsSec[i][j]= 768+i*2+j;
      else    fMapEqidsSec[i][j]= 840+i*4+j-2;
    }
  }
  
  for(Int_t i = 0; i<36 ; i++)
  {
    for(Int_t j = 0; j<6; j++)
    {
      if(j<2) fMapEqidsRcu[768+i*2+j]   = i*6 +j;
      else    fMapEqidsRcu[840+i*4+j-2] = i*6 +j;
    }
  }
}

//__________________________________________________________________
void AliTPCMonitor::FillGlobal(Int_t sector) 
{
  
  // Fill global histograms with max adc for each channel
  
  TH2S* hglob =0;
  if((sector/18) ==0) hglob =  fHistGlobalMaxA;
  else                hglob =  fHistGlobalMaxC;
  
  Float_t  rotsec  = (2*TMath::Pi()/18.0);
  Float_t  rot     = (-rotsec*(sector%18)  +4*rotsec);
  
  Float_t  m11     =    TMath::Cos(rot);
  Float_t  m12     =    TMath::Sin(rot);
  Float_t  m21     = -1*TMath::Sin(rot);
  Float_t  m22     =    TMath::Cos(rot);
  
  Int_t    max     = 0; // use integer for global view
  
  Double_t xval    = 0.0;
  Double_t yval    = 0.0;
  
  Int_t    pad     = 0;
  Int_t    row     = 0;
  Int_t    padmax  = 0;
  
  Float_t  xdr     = 0;
  Float_t  ydr     = 0;
  
  for(Int_t hw = 0; hw<fHistAddrMaxAdc->GetNbinsX(); hw++)
  {
    max = (Int_t)fHistAddrMaxAdc->GetBinContent(hw);
    if(max!=-1)
    {
      pad    = fMapHand->GetPad(   hw);
      row    = fMapHand->GetPadRow(hw);
      if (row==-1) continue;
      padmax = fMapHand->GetNumofPads(row);
      if (sector%36>17) fMirror=-1;
      else fMirror=1;
      GetXY(xval ,yval , padmax,row ,pad);
      xdr    =  xval*m11 +yval*m12;
      ydr    =  xval*m21 +yval*m22;
      if(hglob->GetBinContent(hglob->GetXaxis()->FindBin(xdr),hglob->GetYaxis()->FindBin(ydr))==0)   hglob->Fill(xdr,ydr,(Int_t)max);
    }
  }
}


//__________________________________________________________________
void AliTPCMonitor::GetXY( Double_t& xval , Double_t& yval , Int_t padmax, Int_t row , Int_t pad) const 
{
  // Get x and y position of pad
  
  if(row<63)
  {
    xval = fMirror*( 2*padmax -4*pad -2);
    yval = 852.25 +7.5*row;
  }
  else
  {
    xval = fMirror*( 3*padmax -6*pad -3);
    if((row-63)<63) {	  yval = 10*(row-63) +1351;	}
    else 	      {	  yval = 15*(row-63-64)+1993.5;	}
  }
  
} 

//__________________________________________________________________
Int_t AliTPCMonitor::GetPadAtX(Float_t xval, Int_t row, Int_t padmax) const
{
  // Get pad number at given position in x
  Int_t pad    = 0;
  
  if(row<63) {pad = (Int_t)( ( (xval/fMirror) +2 -2*padmax)/-4);}
  else       {pad = (Int_t)( ( (xval/fMirror) +3 -3*padmax)/-6);}
  
  if(pad>=padmax) return -1;
  else            return pad ;
  
}

//__________________________________________________________________
Int_t AliTPCMonitor::GetPadAtX(Float_t xval, Int_t row) const
{
  
 // Get pad number at given position in x
  
  Int_t padmax = fMapHand->GetNumofPads(row);
  Int_t pad    = 0;
  
  if(row<63) {pad = (Int_t)( ( (xval/fMirror) +2 -2*padmax)/-4);}
  else       {pad = (Int_t)( ( (xval/fMirror) +3 -3*padmax)/-6);}
  
  if(pad>=padmax) return -1;
  else            return pad ;
  
}

//__________________________________________________________________
void AliTPCMonitor::DrawHists(Int_t histos) 
{
  
  // Draw sets of histograms
  // histos==1 : 2Dim histos for MAX adc  and add executables
  // histos==2 : distributions max/rms/sum
  // histos==3 : global max adc for specified SideA/C
  
  
  if(fVerb)    cout << " Draw histos " << endl;
  TString cside;
  if(GetLastSector()/18==0 ) cside="A";
  else                       cside="C";
  
  TString titleSEC  = Form("Sector %i Side %s Run : %05i EventID %i "       ,GetLastSector()%18,cside.Data(),fRunId, fEventNumber);
  TString titleEvent= Form("Time <-> Channles  %s"                          ,titleSEC.Data());
  TString titleIROC = Form("IROC %s"                                        ,titleSEC.Data());
  TString titleOROC = Form("OROC %s"                                        ,titleSEC.Data());
  
  TString titleMAX  = Form("Max (timebin: %i,%i) %s"                        ,GetRangeMaxAdcMin(),GetRangeMaxAdcMax(),titleSEC.Data());
  TString titleSUM  = Form("Sum (timebin: %i,%i) %s"                        ,GetRangeSumMin()   ,GetRangeSumMax()   ,titleSEC.Data());
  TString titleBASE = Form("Baseline RMS<->Mean  (timebin: %i-%i) %s"       ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC.Data());
  TString titleMEAN = Form("Baseline Mean (timebin: %i-%i) %s"              ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC.Data());
  TString titleRMS  = Form("Baseline RMS (timebin: %i-%i) %s"               ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC.Data());
  
  if(histos==1)
  {
      // IROC _______________________________________________________________
    TCanvas* ciroc = 0;
    ciroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("ciroc");
    if(!ciroc)
    {
      ciroc = CreateCanvas("ciroc");
      fExecPlaneMax=0;
    }
    ciroc->cd();
    
    fHistIROC->SetXTitle("row");
    fHistIROC->SetYTitle("pad");
    if(GetPedestals()) fHistIROC->SetZTitle("max ADC (baseline sub)");
    else               fHistIROC->SetZTitle("max ADC ");
    fHistIROC->SetTitle(titleIROC.Data());
    fHistIROC->SetMinimum(0.01);
    fHistIROC->Draw("COLZ");
    ciroc->UseCurrentStyle();
    
    
    fHistIROCTime->SetXTitle("row"); fHistIROCTime->SetZTitle("peak time (fit)");       fHistIROCTime->SetYTitle("pad");      fHistIROCTime->SetTitle(titleIROC.Data());
    fHistIROCRMS->SetXTitle("row");  fHistIROCRMS->SetZTitle( "baseline rms (ADC)");    fHistIROCRMS->SetYTitle("pad");       fHistIROCRMS->SetTitle(titleIROC.Data());
    
      // OROC
    TCanvas* coroc = 0;
    coroc =(TCanvas*)gROOT->GetListOfCanvases()->FindObject("coroc");
    if(!coroc) {
      coroc = CreateCanvas("coroc");
      fExecPlaneMax=0;
    }
    coroc->cd();
    
    fHistOROC->SetXTitle("row");
    fHistOROC->SetYTitle("pad");
    if(GetPedestals()) fHistOROC->SetZTitle("max ADC (baseline sub)");
    else               fHistOROC->SetZTitle("max ADC ");
    fHistOROC->SetTitle(titleOROC.Data());
    fHistOROC->SetMinimum(0.01);
    fHistOROC->Draw("COLZ");
    coroc->UseCurrentStyle();
    
    
    fHistOROCTime->SetXTitle("row"); fHistOROCTime->SetZTitle("peak time (fit) (timebins)"); fHistOROCTime->SetYTitle("pad"); fHistOROCTime->SetTitle(titleOROC.Data());
    fHistOROCRMS->SetXTitle("row");  fHistOROCRMS->SetZTitle("baseline rms (ADC)");          fHistOROCRMS->SetYTitle("pad");  fHistOROCRMS->SetTitle(titleOROC.Data());
    
    // SUM
    TString namesum=Form("ADC sum (bins: %i, %i)",GetRangeSumMin() ,GetRangeSumMax() );
    fHistIROCSUM->SetXTitle("row");
    fHistIROCSUM->SetZTitle(namesum.Data());
    fHistIROCSUM->SetYTitle("pad");
    fHistIROCSUM->SetTitle(titleIROC.Data());
    
    fHistOROCSUM->SetXTitle("row");
    fHistOROCSUM->SetZTitle(namesum);
    fHistOROCSUM->SetYTitle("pad");
    fHistOROCSUM->SetTitle(titleOROC.Data());
    
    // BASE
    TString namebase=Form("base mean (timbebin: %i, %i )",GetRangeBaseMin(),GetRangeBaseMax());
    fHistIROCBASE->SetXTitle("row");
    fHistIROCBASE->SetZTitle(namebase.Data());
    fHistIROCBASE->SetYTitle("pad");
    fHistIROCBASE->SetTitle(titleIROC.Data());
    
    fHistOROCBASE->SetXTitle("row");
    fHistOROCBASE->SetZTitle(namebase);
    fHistOROCBASE->SetYTitle("pad");
    fHistOROCBASE->SetTitle(titleOROC.Data());
    
    if(fHistIROCClone) fHistIROCClone->Delete();
    if(fHistOROCClone) fHistOROCClone->Delete();
    fHistIROCClone = (TH2F*)fHistIROC->Clone("fHistIROCClone");
    fHistOROCClone = (TH2F*)fHistOROC->Clone("fHistOROCClone");
    
      // Executables
    if(fExecPlaneMax==0)
    {
      TString carry1=Form(".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
      ciroc->AddExec("pad",carry1.Data());
      coroc->AddExec("pad",carry1.Data());
      
      TString carry2=Form(".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
      ciroc->AddExec("row",carry2.Data());
      coroc->AddExec("row",carry2.Data());
      fExecPlaneMax=1;
    }
    coroc->Update();
    ciroc->Update();
  }
  else if(histos==2)
  {
      // MAX ADC distribution  ____________________________________________
    TCanvas* cmax = 0;
    cmax = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cmax");
    if(!cmax)  cmax = CreateCanvas("cmax");
    
    cmax->cd();
    fHistDistrMaxIROC->GetXaxis()->SetRangeUser(0.0,1000.0);
    fHistDistrMaxIROC->SetXTitle("max ADC (ADC)");
    fHistDistrMaxIROC->SetYTitle("counts");
    fHistDistrMaxIROC->SetTitle(titleMAX.Data());
    fHistDistrMaxIROC->Draw("");
    fHistDistrMaxOROC->SetLineColor(2);
    fHistDistrMaxOROC->Draw("same");
    
    if(fHistDistrMaxOROC->GetMaximum()> fHistDistrMaxIROC->GetMaximum())  fHistDistrMaxIROC->SetMaximum(fHistDistrMaxOROC->GetMaximum()*1.1);
    
    TLegend* legio = new TLegend(0.6,0.6,0.8,0.8);
    legio->SetFillColor(0);
    legio->AddEntry(fHistDistrMaxIROC,"IROC","l");
    legio->AddEntry(fHistDistrMaxOROC,"OROC","l");
    legio->Draw("same");
    
      // ADC sum distribution
    TCanvas* csum = 0;
    csum = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("csum");
    if(!csum)  csum = CreateCanvas("csum") ;
    csum->cd();
    
    fHistDistrSumIROC->SetXTitle("sum ADC (ADC)");
    fHistDistrSumIROC->SetYTitle("counts");
    fHistDistrSumIROC->SetTitle(titleSUM.Data());
    fHistDistrSumIROC->Draw("");
    fHistDistrSumOROC->SetLineColor(2);
    fHistDistrSumOROC->Draw("same");
    if(fHistDistrSumOROC->GetMaximum()> fHistDistrSumIROC->GetMaximum())  fHistDistrSumIROC->SetMaximum(fHistDistrSumOROC->GetMaximum()*1.1);
    legio->Draw("same");
    
      // BASELINE MEAN distribution
    TCanvas* cbasemean = 0;
    cbasemean = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cbasemean");
    if(!cbasemean)  cbasemean = CreateCanvas("cbasemean");
    cbasemean->cd();
    
    fHistDistrBaseMeanIROC = fHistDistrBase2dIROC->ProjectionX("fHistDistrBaseMeanIROC");
    fHistDistrBaseMeanIROC->SetXTitle("base mean (ADC)");
    fHistDistrBaseMeanIROC->SetYTitle("counts");
    fHistDistrBaseMeanIROC->SetTitle(titleMEAN.Data());
    fHistDistrBaseMeanIROC->Draw("");
    
    fHistDistrBaseMeanOROC = fHistDistrBase2dOROC->ProjectionX("fHistDistrBaseMeanOROC");
    fHistDistrBaseMeanOROC->SetLineColor(2);
    fHistDistrBaseMeanOROC->Draw("same");
    if(fHistDistrBaseMeanOROC->GetMaximum()>fHistDistrBaseMeanIROC->GetMaximum())   fHistDistrBaseMeanIROC->SetMaximum(fHistDistrBaseMeanOROC->GetMaximum()*1.1);
    legio->Draw("same");
    
    TCanvas* cbaserms = 0;
    cbaserms = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cbaserms");
    if(!cbaserms) cbaserms = CreateCanvas("cbaserms") ;
    cbaserms->cd();
    
      // BASELINE RMS distribution
    fHistDistrBaseRmsIROC = fHistDistrBase2dIROC->ProjectionY("fHistDistrBaseRmsIROC");
    fHistDistrBaseRmsIROC->SetXTitle("base rms (ADC)");
    fHistDistrBaseRmsIROC->SetYTitle("counts");
    fHistDistrBaseRmsIROC->SetTitle(titleRMS.Data());
    fHistDistrBaseRmsIROC->Draw("");
    
    fHistDistrBaseRmsOROC = fHistDistrBase2dOROC->ProjectionY("fHistDistrBaseRmsOROC");
    fHistDistrBaseRmsOROC->SetLineColor(2);
    fHistDistrBaseRmsOROC->Draw("same");
    if(fHistDistrBaseRmsOROC->GetMaximum()>fHistDistrBaseRmsIROC->GetMaximum())  fHistDistrBaseRmsIROC->SetMaximum(fHistDistrBaseRmsOROC->GetMaximum()*1.1);
    legio->Draw("same");
    
    cmax->Update();
    csum->Update();
    cbasemean->Update();
    cbaserms->Update();
  }
  else  if(histos==3)
  {
      // GLOBAL MAX ADC _________________________________
    if(GetProcNextEvent()==1)
    {
      TCanvas* cglobA =0;
      TCanvas* cglobC =0;
      
      if(!(cglobC=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("SIDE C all"))) cglobC = CreateCanvas("SIDE C all");
      if(!(cglobA=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("SIDE A all"))) cglobA = CreateCanvas("SIDE A all");
      
      TString globtitle1=Form("SIDE A Run %05i (EventID %i)",fRunId,fEventNumber);
      TString globtitle2=Form("SIDE C Run %05i (EventID %i)",fRunId,fEventNumber);
      
      fHistGlobalMaxA->SetTitle(globtitle1.Data());
      fHistGlobalMaxC->SetTitle(globtitle2.Data());
      fHistGlobalMaxA->SetXTitle("x/mm");
      fHistGlobalMaxA->SetYTitle("y/mm");
      fHistGlobalMaxC->SetXTitle("x/mm");
      fHistGlobalMaxC->SetYTitle("y/mm");
      
      if(GetPedestals()==0)     {	      fHistGlobalMaxA->SetZTitle("max adc (not baseline sub)");   fHistGlobalMaxC->SetZTitle("max adc (not baseline sub)");  }
      else 	                    {	      fHistGlobalMaxA->SetZTitle("max adc ");	                  fHistGlobalMaxC->SetZTitle("max adc ");	                   }
      
      fHistGlobalMaxA->SetMinimum(0.01);
      fHistGlobalMaxC->SetMinimum(0.01);
      
      cglobC->cd() ; fHistGlobalMaxC->Draw("COLZ");
      cglobA->cd() ; fHistGlobalMaxA->Draw("COLZ");
      
      TString nameom=Form(".x  %s/TPC/AliTPCMonitorExec.C(3)",gSystem->Getenv("ALICE_ROOT"));
      
      if(fExecGlob==0)
      {
        if(fVerb)cout << " set exec " << nameom.Data() << endl;
        cglobC->AddExec("glob",nameom.Data());
        cglobA->AddExec("glob",nameom.Data());
        fExecGlob = 1;
      }
      else
      {
        cglobC->DeleteExec("glob");
        cglobA->DeleteExec("glob");
        
        if(fVerb)  cout << " set exec " << nameom.Data() << endl;
        cglobC->AddExec("glob",nameom.Data());
        cglobA->AddExec("glob",nameom.Data());
        
      }
      cglobC->Update();
      cglobA->Update();
    }
    
  }
}



//__________________________________________________________________
void AliTPCMonitor::DrawRMSMap() 
{
  // Draw 2Dim rms histos for IROC and OROC
  // and set executables for canvases
  
  TCanvas* crmsoroc =0;
  TCanvas* crmsiroc =0;
  if(!(crmsoroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crmsoroc")))    crmsoroc    = CreateCanvas("crmsoroc");
  if(!(crmsiroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crmsiroc")))    crmsiroc    = CreateCanvas("crmsiroc");
  
  crmsiroc->cd();  fHistIROCRMS->Draw("COLZ");
  crmsoroc->cd();  fHistOROCRMS->Draw("COLZ");
  
  TString carry1=Form(".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
  TString carry2=Form(".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
  
  if(fExecPadIrocRms==0)
  {
    crmsiroc->AddExec("pad",carry1.Data());
    crmsiroc->AddExec("row",carry2.Data());
    fExecPadIrocRms=1;
  }
  
  if(fExecPadOrocRms==0)
  {
    crmsoroc->AddExec("pad",carry1.Data());
    crmsoroc->AddExec("row",carry2.Data());
    fExecPadOrocRms=1;
  }
  
  crmsiroc->Update();
  crmsoroc->Update();
  
  DrawHists(2);
  
}

//__________________________________________________________________
void AliTPCMonitor::ExecPad() 
{
  
  // Executable for Pad
  // Show time profile for channel the mouse is pointing at
  
  Int_t event = gPad->GetEvent();
  if (event != 51)   return;
  
  TObject *select = gPad->GetSelected();
  if(!select)    return;
  if(!select->InheritsFrom("TH2")) { return;  }
  gPad->GetCanvas()->FeedbackMode(kTRUE);
  
  // get position
  Int_t    px        = gPad->GetEventX();
  Int_t    py        = gPad->GetEventY();
  Float_t  upy       = gPad->AbsPixeltoY(py);
  Float_t  upx       = gPad->AbsPixeltoX(px);
  Float_t  y         = gPad->PadtoY(upy);
  Float_t  x         = gPad->PadtoX(upx);
  
  Int_t    setrange  = 0;
  
  TCanvas* cpad      = 0;
  //  Char_t   namehist[50];
  TString   projhist;
  TString   namesel;
  TString   namecanv;
  
  Int_t    xbinmin  = 0;
  Int_t    xbinmax  = 0;
  Float_t  ybinmin  = 0;
  Float_t  ybinmax  = 0;
  Int_t    rocid     = 0;
  // Check wich Canvas executed the event
  TH2S* fHistIndex=0;
  namesel=select->GetName();
  if(namesel=="fHistOROC" || namesel=="fHistOROCRMS" || namesel=="fHistOROCTime" )
  {
    rocid = 1;
    fPadUsedRoc =1;
    projhist=Form("ProjectionOROC");
    namecanv=Form("coroc_ch");
    fHistIndex = fHistOROCIndex;
  }
  if(namesel=="fHistIROC" || namesel=="fHistIROCRMS" || namesel=="fHistIROCTime" )
  {
    rocid = 0;
    fPadUsedRoc=0;
    projhist=Form("ProjectionIROC");
    namecanv=Form("ciroc_ch");
    fHistIndex = fHistIROCIndex;
  }
  
  // Check if Canvas already existed and get Ranges from former Prjection histogram
  if((cpad=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(namecanv)))
  {
    cpad->cd();
    if(gROOT->Get(projhist))
    {
      setrange = 1;
      xbinmin = ((TH1D*)gROOT->Get(projhist))->GetXaxis()->GetFirst();
      xbinmax = ((TH1D*)gROOT->Get(projhist))->GetXaxis()->GetLast();
      ybinmin = ((TH1D*)gROOT->Get(projhist))->GetMinimum();
      ybinmax = ((TH1D*)gROOT->Get(projhist))->GetMaximum();
      delete gROOT->Get("legfit");
      delete gROOT->Get("fg");
    }
  }
  else
  {
    cpad = CreateCanvas(namecanv); cpad->cd();
  }
  
  // Get Bin
  if (!fHistIndex) return;
  Int_t testy = fHistIndex->GetYaxis()->FindBin(y);
  Int_t testx = fHistIndex->GetXaxis()->FindBin(x);
  Int_t binchannel = (Int_t)fHistIndex->GetCellContent(testx,testy);
  if(binchannel>30000 || binchannel<0) return;
  
  TH1D *hp=(TH1D*)gROOT->Get(projhist);
  if(!hp)  hp=new TH1D(projhist,projhist,GetTimeBins(),0,GetTimeBins());//delete gROOT->Get(projhist);
  // Get Projection
//   TH1D *hp = (TH1D*)(((TH1D*)fHistChannelTime->ProjectionY("hp",binchannel,binchannel))->Clone(projhist));
//   TH1D *hp = fHistChannelTime->ProjectionY(projhist,binchannel,binchannel);
//   hp->Reset();
  Int_t nbinsx=fHistChannelTime->GetNbinsX();
  for (Int_t itb=0;itb<GetTimeBins();++itb){
    hp->GetArray()[itb]=fHistChannelTime->GetArray()[binchannel+itb*(nbinsx+2)];
  }
//   printf("hi\n");
//   return;
  
  // Make title and Pave for channel Info
  Int_t npadRow , npad  , nhw , nmax , hwadd;
  
  hwadd   = (Int_t)fHistChannelTime->GetCellContent(binchannel,0);
  fPadUsedHwAddr = hwadd;
  
  if(rocid==0)npadRow = fMapHand->GetPadRow(hwadd);
  else        npadRow = fMapHand->GetPadRow(hwadd)-63;
  npad                = fMapHand->GetPad(hwadd);
  nhw                 = hwadd;
  nmax                = (Int_t)hp->GetMaximum();
  
  
  TPaveText*  legstat = new TPaveText(0.18,0.65,0.3,0.8,"NDC");
  
  Int_t   connector   =  fMapHand->GetFECconnector(hwadd);
  Int_t   fecnr       =  fMapHand->GetFECfromHw(hwadd);
  Int_t   fecch       =  fMapHand->GetFECchannel(hwadd);
  Int_t   altrochip   =  fMapHand->GetAltro(hwadd);
  Int_t   altrochannel=  (fMapHand->GetAltroChannel(hwadd))%16;
  Int_t   fecloc      =  fMapHand->U2fGetFECinRCU(fecnr) ;
  Int_t   feclocbran  =  fMapHand->U2fGetFECinBranch(fecnr);
  Int_t   branch      =  fMapHand->U2fGetBranch(fecnr);
  
  
  Short_t fecget      = (hwadd & fgkHwMaskFEC)   >> 7;
  Short_t branchget   = (hwadd & fgkHwMaskBranch)>> 11;
  
  legstat->AddText(Form("Branch (map) \t %i (%i) \n",branchget,branch));
  legstat->AddText(Form("Fec in patch \t %i \n",fecloc));
  legstat->AddText(Form("Fec in branch (map)\t %i (%i)\n",fecget,feclocbran));
  legstat->AddText(Form("Connector  \t %i \n",connector));
  legstat->AddText(Form("Fec No.   \t %i \n",fecnr));
  legstat->AddText(Form("Fec chan  \t %i \n",fecch));
  legstat->AddText(Form("Altro chip\t %i \n",altrochip));
  legstat->AddText(Form("Altro chan\t %i \n",altrochannel));
  
  TString title=Form("Row=%d Pad=%d Hw =%d maxADC =%d count =%d",npadRow,npad,nhw,nmax,binchannel);
  
//   hp->SetName(projhist);
  hp->SetTitleSize(0.04);
  hp->SetTitle(title.Data());
  hp->SetYTitle("ADC");
  hp->SetXTitle("Timebin");
  hp->GetXaxis()->SetTitleColor(1);
  
  if(setrange)
  {
    hp->GetXaxis()->SetRange(xbinmin,xbinmax);
    hp->SetMinimum(ybinmin);
    hp->SetMaximum(ybinmax);
  }
  else
  {
    hp->SetMinimum(0.0);
    hp->SetMaximum(1000.0);
  }
  
  cpad->cd();
  hp->Draw();
  
  // Make Fit to peak
  if(GetPedestals() && fDisableFit==0)
  {
    Int_t maxx  =    (Int_t)fHistAddrMaxAdcX->GetBinContent(hwadd);
    Float_t max  =  (Float_t)fHistAddrMaxAdc->GetBinContent(hwadd);
    Float_t base =  (Float_t)fHistAddrBaseMean->GetBinContent(hwadd);
    if(base!=0)
    {
      if( ((max+base)/base)>1.2)
      {
        TF1* fg = new TF1("fg",AliTPCMonitor::Gamma4,maxx-5,maxx+5,4);
        fg->SetParName(0,"Normalisation");
        fg->SetParName(1,"Minimum");
        fg->SetParName(2,"Width");
        fg->SetParName(3,"Base");
        fg->SetParameter(0,max);
        fg->SetParameter(1,maxx-2);
        fg->SetParameter(2,1.5);
        fg->FixParameter(3,0);
        fg->SetLineColor(4);
        fg->SetLineWidth(1);
        hp->Fit("fg","RQ");
        
        TLegend* legfit = new TLegend(0.6,0.7,0.7,0.8);
        legfit->AddEntry("fg","#Gamma 4 fit","l");
        legfit->SetFillColor(0);
        legfit->SetName("legfit");
        legfit->Draw("same");
      }
    }
  }
  legstat->SetFillColor(0);
  legstat->Draw("same");
  cpad->Update();
  return;
}

//__________________________________________________________________
void AliTPCMonitor::ExecRow() 
{
  
  // Executable for Pad
  // Show profile of max adc over given pad row
  // and 2dim histo  adc(pad-in-row,time bin)
  
  Int_t event = gPad->GetEvent();
  if (event != 61)    return;
  gPad->cd();
  TObject *select = gPad->GetSelected();
  if(!select)    return;
  if(!select->InheritsFrom("TH2")) {   return;  }
  
  Int_t rocid = 0;
  
  const Char_t *rowhist="hrowtime";
  const Char_t *rowhistxmax="hxmax";
  const Char_t *rowhistmax="hrowmax";
  
  // get position
  Int_t   px  = gPad->GetEventX();
  Int_t   py  = gPad->GetEventY();
  Float_t upy = gPad->AbsPixeltoY(py);
  Float_t upx = gPad->AbsPixeltoX(px);
  Float_t y   = gPad->PadtoY(upy);
  Float_t x   = gPad->PadtoX(upx);
  
  TCanvas*crowtime   = 0;
  TCanvas*crowmax    = 0;
  TCanvas*cxmax      = 0;
  
  TH2S*  fHistIndex  = 0;
  
  // ranges from already existing histos
  Int_t    rowtimexmin = 0;
  Int_t    rowtimexmax = 0;
  Int_t    rowtimeymin = 0;
  Int_t    rowtimeymax = 0;
  Float_t  rowtimezmin = 0;
  Float_t  rowtimezmax = 0;
  
  Int_t    profrowxmin = 0;
  Int_t    profrowxmax = 0;
  Double_t profrowymin = 0;
  Double_t profrowymax = 0;
  
  Int_t    profxxmin   = 0;
  Int_t    profxxmax   = 0;
  Double_t profxymin   = 0;
  Double_t profxymax   = 0;
  
  
  Int_t    setrange      = 0;
  
  TString name=select->GetName();
  if(     name=="fHistIROC" || name=="fHistIROCRMS" ) { fHistIndex = fHistIROCIndex;     rocid =1;   }
  else if(name=="fHistOROC" || name=="fHistOROCRMS" ) { fHistIndex = fHistOROCIndex;     rocid =2;   }
  else  { cout << " not implemented for this histo " << endl; return; }
  
  gPad->GetCanvas()->FeedbackMode(kTRUE);
  
  
  
  // check if canvases exist //
  crowtime  = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crowtime");
  crowmax   = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crowmax");
  cxmax     = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cxmax");
  
  if(!crowtime)   crowtime   = CreateCanvas("crowtime") ;
  if(!crowmax)    crowmax    = CreateCanvas("crowmax")  ;
  if(!cxmax  )    cxmax      = CreateCanvas("cxmax")    ;
  
  // check ranges of already existing histos
  if(gROOT->Get(rowhist))
  {
    rowtimexmin  = ((TH2F*)gROOT->Get(rowhist))->GetXaxis()->GetFirst();
    rowtimexmax  = ((TH2F*)gROOT->Get(rowhist))->GetXaxis()->GetLast();
    rowtimeymin  = ((TH2F*)gROOT->Get(rowhist))->GetYaxis()->GetFirst();
    rowtimeymax  = ((TH2F*)gROOT->Get(rowhist))->GetYaxis()->GetLast();
    rowtimezmin  = ((TH2F*)gROOT->Get(rowhist))->GetMinimum();
    rowtimezmax  = ((TH2F*)gROOT->Get(rowhist))->GetMaximum();
    
    profrowxmin  = ((TH1F*)gROOT->Get(rowhistmax))->GetXaxis()->GetFirst();
    profrowxmax  = ((TH1F*)gROOT->Get(rowhistmax))->GetXaxis()->GetLast();
    profrowymin  = ((TH1F*)gROOT->Get(rowhistmax))->GetMinimum();
    profrowymax  = ((TH1F*)gROOT->Get(rowhistmax))->GetMaximum();
    
    profxxmin    = ((TH1F*)gROOT->Get(rowhistxmax))->GetXaxis()->GetFirst();
    profxxmax    = ((TH1F*)gROOT->Get(rowhistxmax))->GetXaxis()->GetLast();
    profxymin    = ((TH1F*)gROOT->Get(rowhistxmax))->GetMinimum();
    profxymax    = ((TH1F*)gROOT->Get(rowhistxmax))->GetMaximum();
    
    setrange =1;
    
    delete gROOT->Get(rowhist);
    delete gROOT->Get(rowhistmax);
    delete gROOT->Get("hxmax");
    delete gROOT->Get("legrow");
  }
  
  // get channel for xy bin -> getRow -> getNrows -> getHw for each Pad in Row -> get channel for each hw -> make proj
  Int_t testy      = fHistIndex->GetYaxis()->FindBin(y);
  Int_t testx      = fHistIndex->GetXaxis()->FindBin(x);
  Int_t binchannel = (Int_t)fHistIndex->GetCellContent(testx,testy);
  
  if(binchannel>30000)    return;
  if(binchannel<=0 ) { crowtime->Update() ;    crowmax->Update() ;    return ;  }
  
  // get hwaddress
  Int_t hwadd     = (Int_t)fHistChannelTime->GetCellContent(binchannel,0);
  Int_t row       = fMapHand->GetPadRow(hwadd);
  Int_t pad       = fMapHand->GetPad(hwadd)   ;
  Int_t numofpads =  fMapHand->GetNumofPads(row);
  
  // create histos
  TH2F *hrowtime = new TH2F(rowhist     , ""  ,numofpads,0,numofpads,GetTimeBins(),0.0,GetTimeBins());
  TH1F *hrowmax  = new TH1F(rowhistmax , ""  ,numofpads,0,numofpads);
  TH1F *hxmax    = new TH1F(rowhistxmax, ""  ,159,0,159      );
  
  // Row profile ///////////
  if(fVerb) cout << " Number of pads " << numofpads << endl;
  for(Int_t padnr = 0; padnr<numofpads;padnr++)
  {
    Int_t addrinrow = fMapHand->GetPadAddInRow(row,padnr );
    Int_t channel   = (Int_t)fHistAddrMapIndex->GetBinContent(addrinrow);
    if(channel==-1) continue;
    
    hrowmax->SetBinContent(padnr+1,fHistAddrMaxAdc->GetBinContent(addrinrow));
    TH1D *hp = fHistChannelTime->ProjectionY("hp",channel,channel);
    for(Int_t time = 0;time<GetTimeBins();time++) {
      
      Float_t val = hp->GetBinContent(time);
      hrowtime->SetCellContent(padnr+1,time+1,val);
    }
  }
  
  // X profile  /////////////
  Double_t xval  = 0.0;
  Double_t yval  = 0.0;
  GetXY(xval,yval,numofpads,row,pad);
  
  Int_t padnr = 0;
  Int_t hw    = 0;
  for(Int_t nrow = 0; nrow<159; nrow++)
  {
    padnr = GetPadAtX(xval,nrow);
    if(padnr>=0)
    {
      hw = fMapHand->GetPadAddInRow(nrow,padnr);
      if(fPadMapHw[hw]==-1){ continue                                                      ; }
      else                { hxmax->SetBinContent(nrow+1,fHistAddrMaxAdc->GetBinContent(hw))   ;	}
    }
  }
  
  cxmax->cd();
  hxmax->SetTitle(Form("max adc in pads at x=%5.1f mm",xval));
  hxmax->SetXTitle("row");
  if(!GetPedestals()) hxmax->SetYTitle("max adc (baseline sub.)");
  else                hxmax->SetYTitle("max adc ");
  hxmax->SetMinimum(0.01);
  hxmax->Draw("l");
  
  if(setrange)
  {
    hxmax->GetXaxis()->SetRange(profxxmin,profxxmax);
    hxmax->SetMinimum(profxymin);
    hxmax->SetMaximum(profxymax);
  }
  
  cxmax->Update();
  
  crowtime->cd();
  TString title;
  TString titlemax;
  if(rocid==1) {
    title.Form("%s Row=%d",((TH2*)select)->GetTitle(),row)   ;
    titlemax.Form("IROC  max/sum Row=%d",row   );
  }  else  {
    title.Form("%s Row=%d",((TH2*)select)->GetTitle(),row-63);
    titlemax.Form("OROC  max/sum Row=%d",row-63);
  }
  
  if(fVerb) cout << " set name " << endl;
  
  
  // row vs time
  crowtime->cd();
  hrowtime->SetTitleSize(0.04);
  hrowtime->SetTitle(title.Data());
  hrowtime->SetYTitle("timbin");
  hrowtime->SetXTitle("pad in row");
  hrowtime->SetZTitle("signal (ADC)");
  
  hrowtime->GetXaxis()->SetTitleColor(1);
  hrowtime->SetMaximum(1000.0);
  hrowtime->SetMinimum(0.0);
  
  if(setrange)
  {
    hrowtime->GetXaxis()->SetRange(rowtimexmin,rowtimexmax);
    hrowtime->GetYaxis()->SetRange(rowtimeymin,rowtimeymax);
    hrowtime->SetMinimum(rowtimezmin);
    hrowtime->SetMaximum(rowtimezmax);
  }
  
  hrowtime->Draw("COLZ");
  crowtime->UseCurrentStyle();
  crowtime->Update();
  
  // max and sum /////////////////////////
  crowmax->cd();
  if(setrange) {
    hrowmax->GetXaxis()->SetRange(profrowxmin,profrowxmax);
    hrowmax->SetMinimum(profrowymin);
    hrowmax->SetMaximum(profrowymax);
  }
  hrowmax->SetTitleSize(0.04);
  hrowmax->SetTitle(title.Data());
  hrowmax->SetYTitle("max adc");
  hrowmax->SetXTitle("pad in row");
  hrowmax->GetXaxis()->SetTitleColor(1);
  
  hrowmax->SetLineColor(2);
  hrowmax->Draw("l");
  crowmax->Update();
  
  return;
}

//__________________________________________________________________
void AliTPCMonitor::Write10bitChannel()
{
  
  // Write 10 bit words form histogram for active(last pointed)  channel
  
  if(fPadUsedHwAddr==-1){ AliWarning(" No active pad "); return ;}
  
  Int_t  pad     = (Int_t)fMapHand->GetPad(   fPadUsedHwAddr);
  Int_t  row     = (Int_t)fMapHand->GetPadRow(fPadUsedHwAddr);
  Int_t  channel = (Int_t)fPadMapHw[fPadUsedHwAddr];
  
  TString filenameroot;
  TString filenamedat;
  TString projhist;
  
  if(fPadUsedRoc==1) { projhist.Form("ProjectionOROC"); }
  if(fPadUsedRoc==0) { projhist.Form("ProjectionIROC"); }
  
  filenamedat.Form("Channel_Run%05i_EventID_%i_Pad_%i_Row_%i.dat"      ,fRunId,fEventNumber,pad,row);
  filenameroot.Form("Channel_Run%05i_EventID_%i_Pad_%i_Row_%i.root"     ,fRunId,fEventNumber,pad,row);
  
  TH1D* hpr = 0;
  if((hpr=(TH1D*)gROOT->Get(projhist)))
  {
      // root histo
    TFile f(filenameroot.Data(),"recreate");
    hpr->Write();
    f.Close();
    
      // raw singal
    ofstream datout(filenamedat.Data(),ios::out);
    datout <<"Timebin \t ADC value " << endl;
    for(Int_t i = 1; i <GetTimeBins(); i++)
    {
      datout << i << " \t \t " << fPad[channel][i] << endl;
    }
    datout.close();
  }
  else
  {
    AliWarning("No projection histo found ");
  }
}

//__________________________________________________________________
void AliTPCMonitor::ExecTransform() 
{
  
  // Make Fast Fourier Transformation for active pad
  // fft is only performed for a data sample of size 2^n
  // reduce window according to largest  power of 2 which is smaller than the viewing  range
  
  TString namecanv;
  TString namecanv2;
  TString projhist;
  TString namehtrimag;
  TString namehtrreal;
  TString namehtrmag;
  
  if(fPadUsedRoc==1) {    namecanv="coroc_ch_trans"; namecanv2="coroc_ch_trans2";  projhist="ProjectionOROC";  }
  if(fPadUsedRoc==0) {    namecanv="ciroc_ch_trans"; namecanv2="ciroc_ch_trans2";  projhist="ProjectionIROC";  }
  
  TH1D*  hproj = 0;
  
  if((TH1D*)gROOT->Get(projhist)==0){AliWarning("Proj histo does not exist \n Move mouse over 2d histo choose channel \n and drag mouse form histo again!");  return ;}
  else      hproj = (TH1D*)gROOT->Get(projhist.Data()) ;
  
  
  if(fPadUsedRoc==1) {  namehtrimag="htransimagfreq_oroc";    namehtrreal="htransrealfreq_oroc"; namehtrmag="htransmagfreq_oroc";  }
  else               {  namehtrimag="htransimagfreq_iroc";    namehtrreal="htransrealfreq_iroc"; namehtrmag="htransmagfreq_iroc";  }
  
  if( gROOT->Get(namehtrimag.Data()))  delete  gROOT->Get(namehtrimag.Data());
  if( gROOT->Get(namehtrreal.Data()))  delete  gROOT->Get(namehtrreal.Data());
  if( gROOT->Get(namehtrmag.Data()))  delete  gROOT->Get(namehtrmag.Data());
  
  TCanvas *ctrans = 0;
  if(!(ctrans = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(namecanv.Data())))
  {
    ctrans = CreateCanvas(namecanv.Data());
    ctrans->Divide(1,2);
  }
  TCanvas *ctrans2 = 0;
  if(!(ctrans2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(namecanv2.Data())))
  {
    ctrans2 = CreateCanvas(namecanv2.Data());
//      ctrans2->Divide(1,2);
  }
  
  Int_t binfirst  =  hproj->GetXaxis()->GetFirst();
  Int_t binlast   =  hproj->GetXaxis()->GetLast();
  Int_t bins       =  binlast -binfirst +1;
  
  Int_t power = 0;
  for(Int_t pot = 0; pot<=10 ; pot++)
  {
    Int_t comp =  (Int_t)TMath::Power(2,pot);
    if(bins>=comp)power = pot;
  }
  
  bins = (Int_t)TMath::Power(2,power);
  
  // sampling frequency ;
  Double_t  deltat = 1.0/(Float_t)GetSamplingFrequency();
  
  // output histo
  TH1D* htransrealfreq = new TH1D(namehtrreal.Data(),namehtrreal.Data(),10000,-1/(2*deltat),1/(2*deltat));
  TH1D* htransimagfreq = new TH1D(namehtrimag.Data(),namehtrimag.Data(),10000,-1/(2*deltat),1/(2*deltat));
  TH1D* htransmag      = new TH1D(namehtrmag.Data(),namehtrmag.Data(),10000,-1/(2*deltat),1/(2*deltat));
  
  TString titlereal;
  TString titleimag;
  TString titlemag;
  if(fPadUsedRoc==1) {    titlereal="OROC DFT real part";  titleimag="OROC DFT imag part";  titlemag="OROC DFT magnitude";  }
  else {                  titlereal="IROC DFT real part";  titleimag="IROC DFT imag part";  titlemag="IROC DFT magnitude";  }
  
  htransrealfreq->SetTitle(titlereal.Data());  htransrealfreq->SetXTitle("f/hz");  htransrealfreq->SetYTitle("z_{real}(f)");
  htransimagfreq->SetTitle(titleimag.Data());  htransimagfreq->SetXTitle("f/hz");  htransimagfreq->SetYTitle("z_{imag}(f)");
  htransmag->SetTitle(titlemag.Data());  htransmag->SetXTitle("f/hz");  htransmag->SetYTitle("mag(f)");
  
  // create complex packed data array
  const Int_t kdatasiz = 2*bins;
  Double_t*  data = new Double_t[kdatasiz];
  for(Int_t i=0;i<2*bins;i++)  { data[i]   =  0.0;}
  for(Int_t i=0;i<bins;i++)    { data[2*i] = (Double_t)hproj->GetBinContent(binfirst+i); }
  
  // make fourier transformation
  AliTPCMonitorFFT* four = new AliTPCMonitorFFT();
  four->ComplexRadix2ForwardWrap(data,1,bins);
  
  // write output  and fill histos forward
  Double_t freq =  0.0;
  for(Int_t i=0;i<2*bins;i++)
  {
    if(i<bins)
    {
      if(i<(bins/2))  { freq = i/(bins*deltat)            ; }
      else            { freq = -1*((bins-i)/(bins*deltat)); }
      htransrealfreq->Fill( freq,data[2*i]  );
      htransimagfreq->Fill( freq,data[2*i+1]);
      htransmag->Fill( freq, TMath::Sqrt(data[2*i]*data[2*i]+data[2*i+1]*data[2*i+1]) );
    }
  }
  
  ctrans->cd(1);
  htransrealfreq->Draw();
  ctrans->cd(2);
  htransimagfreq->Draw();
  ctrans->Update();
  ctrans2->cd();
  htransmag->Draw();
  ctrans2->Update();
  delete four;
  delete [] data;
}

//__________________________________________________________________
void AliTPCMonitor::ShowSel(const Int_t* compval)
{
  
  // Show only selected components
  // First restore original histogram from clone
  // Than remove all not matching pads form histos
  
  Int_t   connector   =  0;
  Int_t   fecnr       =  0;
  Int_t   altrochip   =  0;
  Int_t   feclocbran  =  0;
  Int_t   branch      =  0;
  Short_t rcuget      =  0;
  Int_t   emptyI      =  1;
  Int_t   index       = -1;
  Int_t   hwadd       =  0;
  
  Float_t maxiroc     =  fHistIROCClone->GetMaximum();
  Float_t maxoroc     =  fHistOROCClone->GetMaximum();
  
  
  //  restore original histos
  for(Int_t row = 0; row<fkNRowsIroc;  row++)
  {
    for(Int_t pad = 0; pad<fkNPadsIroc;  pad++)
    {
      index = (Int_t)fHistIROCIndex->GetCellContent(row+1,pad+1);
      if(index==-1)continue;
      else  fHistIROC->SetCellContent(row+1,pad+1,fHistIROCClone->GetCellContent(row+1,pad+1));
    }
  }
  for(Int_t row = 0; row<fkNRowsOroc;  row++)
  {
    for(Int_t pad = 0; pad<fkNPadsOroc;  pad++)
    {
      index = (Int_t)fHistOROCIndex->GetCellContent(row+1,pad+1);
      if(index==-1)continue;
      else    fHistOROC->SetCellContent(row+1,pad+1,fHistOROCClone->GetCellContent(row+1,pad+1));
    }
  }
  
  
  // remove not matching entries from fHistIROC/fHistOROC
  
  TH2F* fHist       =0;
  TH2S* fHistIndex  =0;
  Int_t npads       =0;
  Int_t subrows     =0;
  
  for(Int_t row = 0; row< (fkNRowsIroc + fkNRowsOroc);  row++)
  {
    if(row<fkNRowsIroc) {  fHist=fHistIROC ; fHistIndex = fHistIROCIndex; npads = fkNPadsIroc; subrows =0         ;}
    else                {  fHist=fHistOROC ; fHistIndex = fHistOROCIndex; npads = fkNPadsOroc; subrows =fkNRowsIroc;}
    
    for(Int_t pad = 0; pad<npads;  pad++)
    {
      index    = (Int_t)fHistIndex->GetCellContent(row -subrows +1,pad+1);
      if(index==-1)  continue ;
      hwadd    = (Int_t)fHistChannelTime->GetCellContent(index,0);
      
    // global fecnr
      fecnr     =  fMapHand->GetFECfromHw(hwadd);
      if(compval[0]!=-1 && fecnr!=compval[0])      {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
      
    // rcu
      rcuget      = (hwadd & fgkHwMaskRCU)>> 12;
      if(compval[1]!=-1 && rcuget!=compval[1])     {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
      
    // branch
      branch    =  fMapHand->U2fGetBranch(fecnr) ;
      if(compval[2]!=-1 && branch!=compval[2])     {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
      
    // local fec
      feclocbran=   fMapHand->U2fGetFECinBranch(fecnr) ;
      if(compval[3]!=-1 && feclocbran!=compval[3]) { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
      
    // connector
      connector =  fMapHand->GetFECconnector(hwadd);
      if(compval[4]!=-1 && connector!=compval[4])  { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
      
    // Altro chip
      altrochip     =  fMapHand->GetAltro(hwadd);
      if(compval[5]!=-1 && altrochip!=compval[5])      { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
      emptyI =0;
    }
  }
  
  TCanvas* c1 =0;
  TCanvas* c2 =0;
  if(gROOT->GetListOfCanvases()->FindObject("ciroc"))
  {
    c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("ciroc");
    c1->cd() ;
    fHistIROC->Draw("COLZ");
    fHistIROC->SetMaximum(maxiroc);
    fHistIROC->SetMinimum(0.0);
    c1->Update();
  }
  if(gROOT->GetListOfCanvases()->FindObject("coroc"))
  {
    c2 =  (TCanvas*)gROOT->GetListOfCanvases()->FindObject("coroc");
    c2->cd() ;
    fHistOROC->Draw("COLZ");
    fHistOROC->SetMaximum(maxoroc);
    fHistOROC->SetMinimum(0.0);
    c2->Update();
  }
  return ;
}

//__________________________________________________________________
void AliTPCMonitor::ResizeCanv() 
{
  // Resize canvases and delete some of them
  
  TString carry1=Form(".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
  TString carry3=Form(".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
  if(fVerb) cout <<  " canv 1 " << endl;
  
  if(gROOT->GetListOfCanvases()->FindObject(        "coroc_ch")) {  delete gROOT->GetListOfCanvases()->FindObject("coroc_ch") ; }
  if(gROOT->GetListOfCanvases()->FindObject(        "ciroc_ch")) {  delete gROOT->GetListOfCanvases()->FindObject("ciroc_ch") ; }
  
  // for 2dim plots delete create and draw again
  if(gROOT->GetListOfCanvases()->FindObject("ciroc"))
  {
    delete  gROOT->GetListOfCanvases()->FindObject("ciroc");
    TCanvas* ciroc = CreateCanvas("ciroc");
    ciroc->cd();
    fHistIROC->Draw("COLZ");
    ciroc->AddExec("pad",carry1.Data());
    ciroc->AddExec("row",carry3.Data());
    fExecPlaneMax=1;
    ciroc->Update();
  }
  // for 2dim plots delete create and draw again
  if(gROOT->GetListOfCanvases()->FindObject("coroc"))
  {
    delete gROOT->GetListOfCanvases()->FindObject("coroc");
    TCanvas* coroc = CreateCanvas("coroc");
    coroc->cd();
    fHistOROC->Draw("COLZ");
    
    coroc->AddExec("pad",carry1.Data());
    coroc->AddExec("row",carry3.Data());
    coroc->Update();
    fExecPlaneMax=1;
  }
  
  if(gROOT->GetListOfCanvases()->FindObject(       "cbasemean")) {    delete gROOT->GetListOfCanvases()->FindObject("cbasemean"); }
  if(gROOT->GetListOfCanvases()->FindObject(           "cbase")) {    delete gROOT->GetListOfCanvases()->FindObject("cbase");}
  if(gROOT->GetListOfCanvases()->FindObject(        "cbaserms")) {    delete gROOT->GetListOfCanvases()->FindObject("cbaserms");  }
  if(gROOT->GetListOfCanvases()->FindObject(            "cmax")) {    delete gROOT->GetListOfCanvases()->FindObject("cmax");      }
  if(gROOT->GetListOfCanvases()->FindObject(            "csum")) {    delete gROOT->GetListOfCanvases()->FindObject("csum");      }
  if(gROOT->GetListOfCanvases()->FindObject(  "ciroc_ch_trans")) {    delete gROOT->GetListOfCanvases()->FindObject("ciroc_ch_trans");}
  if(gROOT->GetListOfCanvases()->FindObject(  "coroc_ch_trans")) {    delete gROOT->GetListOfCanvases()->FindObject("coroc_ch_trans");}
  if(gROOT->GetListOfCanvases()->FindObject(       "crowtime")) {    delete gROOT->GetListOfCanvases()->FindObject("crowtime"); }
  if(gROOT->GetListOfCanvases()->FindObject(        "crowmax")) {    delete gROOT->GetListOfCanvases()->FindObject("crowmax");  }
  if(gROOT->GetListOfCanvases()->FindObject(          "cxmax")) {    delete gROOT->GetListOfCanvases()->FindObject("cxmax");    }
  if(gROOT->GetListOfCanvases()->FindObject(        "crmsoroc")) {    delete gROOT->GetListOfCanvases()->FindObject("crmsoroc");         fExecPadOrocRms = 0;   }
  if(gROOT->GetListOfCanvases()->FindObject(        "crmsiroc")) {    delete gROOT->GetListOfCanvases()->FindObject("crmsiroc");         fExecPadIrocRms = 0;   }
  if(gROOT->GetListOfCanvases()->FindObject(        "crowmax")) {    delete gROOT->GetListOfCanvases()->FindObject("crowmax");  }
  if(gROOT->GetListOfCanvases()->FindObject(        "crowmax")) {    delete gROOT->GetListOfCanvases()->FindObject("crowmax");  }
  
}




//__________________________________________________________________
Int_t AliTPCMonitor::ExecProcess() 
{
  // Executable for global Histogram
  // Will be called from /TPC/AliTPCMonitorExec.C(3)
  // Call ProcessEvent for same event and sector pointed at
  
  Int_t side   = 0;
  Int_t sector = 0;
  
  Int_t event = gPad->GetEvent();
  if(event != 61)  return -1;
  
  TObject *select = gPad->GetSelected();
  if(!select)  return -1;
  TString name=select->GetName();
  if(!select->InheritsFrom("TH2")) {gPad->SetUniqueID(0);    return -1;  }
  if(       name=="hglobal"  || name=="SIDE A"  ) side = 0;
  else  if( name=="hglobal2" || name=="SIDE C"  ) side = 1;
  
  // get position
  Int_t   px    = gPad->GetEventX();
  Int_t   py    = gPad->GetEventY();
  Float_t upy   = gPad->AbsPixeltoY(py);
  Float_t upx   = gPad->AbsPixeltoX(px);
  Float_t y     = gPad->PadtoY(upy);
  Float_t x     = gPad->PadtoX(upx);
  
  Int_t testy = ((TH2*)select)->GetYaxis()->FindBin(y);
  Int_t testx = ((TH2*)select)->GetXaxis()->FindBin(x);
  if(((TH2*)select)->GetCellContent(testx,testy)==0) return -1 ;
  
  Float_t alpha = 0.0;
  if(x>0.0 && y > 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(x/y)));
  if(x>0.0 && y < 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(y/x)));
  if(x<0.0 && y < 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(x/y)));
  if(x<0.0 && y > 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(y/x)));
  
  if(x>0.0 && y < 0.0)    alpha += (    TMath::Pi()/2);
  if(x<0.0 && y < 0.0)    alpha += (    TMath::Pi());
  if(x<0.0 && y > 0.0)    alpha += (1.5*TMath::Pi());
  
  sector =   (Int_t)(alpha/(2*TMath::Pi()/18.0));
  if(alpha> (sector+0.5)*(2*TMath::Pi()/18.0))  sector+=1;
  
  if(sector==18 && side ==0 ) {
    AliWarning("There was a wromg assignment of sector 0 with sector 18. Check sectors");
    sector =0;
  }
  
  sector = (18-sector +4)%18;
  SetLastSector(sector+ side*18);
  SetProcNextEvent(0);
  
  if(fVerb) cout << "AliTPCMonitor::ExecProcess()  next side          " <<   side    << " next sector        " <<    sector  << endl;
  
  return (Int_t)ProcessEvent();
  
}

//__________________________________________________________________
Int_t AliTPCMonitor::GetRCUPatch(Int_t runid, Int_t eqid) const
{
  
  // Return RCU patch index for given equipment id eqid
  Int_t patch = 0;
  //if(runid>=704)
  if(runid>=0)
  {
    if ( eqid<768 || eqid>983 ) return 0; //no TPC eqid
//     if(eqid>=1000) return 0;
    patch = fMapEqidsRcu[eqid] ;
  }
  else
  {
    if(eqid==408) {patch =  13*6+4 +0;  }
    if(eqid==409) {patch =  13*6+5 +0;  }
    if(eqid==509) {patch =  13*6+0 +0;  }
    if(eqid==512) {patch =  13*6+3 +0;  }
    if(eqid==513) {patch =  13*6+1 +0;  }
    if(eqid==517) {patch =  13*6+2 +0;  }
    
    if(eqid==404) {patch =  4*6+5 +0;   }
    if(eqid==504) {patch =  4*6+4 +0;   }
    if(eqid==407) {patch =  4*6+3 +0;   }
    if(eqid==503) {patch =  4*6+2 +0;   }
    if(eqid==508) {patch =  4*6+1 +0;   }
    if(eqid==506) {patch =  4*6+0 +0;   }
  }
  return patch;
}

//__________________________________________________________________
void AliTPCMonitor::DumpHeader(AliRawReader * reader) const
{
  // Dump Event header for  format
  
  cout << "EventHeader     : fReader->GetEquipmentSize()            :" << reader->GetEquipmentSize()        << endl;
  cout << "EventHeader     : fReader->GetType()                     :" << reader->GetType()                 << endl;
  cout << "EventHeader     : fReader->GetRunNumber()                :" << reader->GetRunNumber()            << endl;
  cout << "EventHeader     : fReader->GetEventId()                  :" << *(reader->GetEventId())           << endl;
  cout << "EventHeader     : fReader->GetLDCId()                    :" << reader->GetLDCId()                << endl;
  cout << "EventHeader     : fReader->GetGDCId()                    :" << reader->GetGDCId()                << endl;
}


//__________________________________________________________________
Double_t AliTPCMonitor::Gamma4(const Double_t* x, const Double_t* par) {
  
  // Gamma4 function used to fit signals
  // Defined in sections: diverging branch set to 0
  
  Double_t val  = 0.0;
  if(x[0] > par[1])
    val = par[0]*exp(4.0)* pow((x[0]-par[1])/par[2],4)*exp(-4.0*(x[0]-par[1])/par[2])+ par[3];
  else
    val = 0;
  return val;
}

//__________________________________________________________________
TCanvas* AliTPCMonitor::CreateCanvas(const Char_t* name)
{
  // Create Canvases
  
  TCanvas* canv =0;
  
  Int_t xoffset  = GetCanvasXOffset();
  Int_t xsize    = GetCanvasXSize();
  Int_t ysize    = GetCanvasYSize();
  Int_t xspace   = GetCanvasXSpace();
  Int_t yspace   = GetCanvasYSpace();

  TString sname=name;
  // ROC 2dim max distribution
  if(     sname=="coroc") {    canv   = new TCanvas("coroc"         ,"coroc"      ,                   -1+xoffset,(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="ciroc") {    canv   = new TCanvas("ciroc"         ,"ciroc"      ,                   -1+xoffset,                     0    ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // ROC  2dim rms distribution
  else if(sname=="crmsoroc") {    canv   = new TCanvas("crmsoroc"      ,"crmsoroc"   ,                   -1+xoffset,(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="crmsiroc") {    canv   = new TCanvas("crmsiroc"      ,"crmsiroc"   ,                   -1+xoffset,                        0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // Global ADC max Histos
  else if(sname=="SIDE C all") {    canv   = new TCanvas("SIDE C all"    ,"SIDE C all" ,   (Int_t)(3*xspace+ xoffset),(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="SIDE A all") {    canv   = new TCanvas("SIDE A all"    ,"SIDE A all" ,   (Int_t)(3*xspace+ xoffset),                        0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // 1 dim max sum basekine distribution
  else if(sname=="cmax") {    canv   = new TCanvas("cmax"          ,"cmax"       ,                   -1+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(sname=="csum") {    canv   = new TCanvas("csum"          ,"csum"       ,               xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(sname=="cbasemean") {    canv   = new TCanvas("cbasemean"     ,"cbasemean"  ,             2*xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(sname=="cbaserms") {    canv   = new TCanvas("cbaserms"      ,"cbaserms"   ,             3*xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  // Projections of single channel
  else if(sname=="coroc_ch") {    canv   = new TCanvas("coroc_ch"      ,"coroc_ch"   ,   (Int_t)(1.5*xspace+xoffset),(Int_t)(yspace+0.5*ysize),(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="ciroc_ch") {    canv   = new TCanvas("ciroc_ch"      ,"ciroc_ch"   ,   (Int_t)(1.5*xspace+xoffset),                       0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // FFT for single channel
  else if(sname=="coroc_ch_trans") {    canv   = new TCanvas("coroc_ch_trans","coroc_ch_trans",(Int_t)(3.0*xspace+xoffset),(Int_t)(yspace+0.5*ysize),(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="ciroc_ch_trans") {    canv   = new TCanvas("ciroc_ch_trans","ciroc_ch_trans",(Int_t)(3.0*xspace+xoffset),                       0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="coroc_ch_trans2") {    canv   = new TCanvas("coroc_ch_trans2","coroc_ch_trans2",(Int_t)(3.0*xspace+xoffset),(Int_t)(yspace+0.5*ysize),(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(sname=="ciroc_ch_trans2") {    canv   = new TCanvas("ciroc_ch_trans2","ciroc_ch_trans2",(Int_t)(3.0*xspace+xoffset),                       0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // row profile histograms
  else if(sname=="crowtime") {    canv   = new TCanvas("crowtime"     ,"crowtime"  ,              1*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(sname=="crowmax") {    canv   = new TCanvas("crowmax"      ,"crowmax"   ,              2*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(sname=="cxmax") {    canv   = new TCanvas("cxmax"        ,"cxmax"     ,              3*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else                                      {    cout   << " Warning Canvas name unknown "  << endl;                                                                                                  return 0   ;  }
}

//__________________________________________________________________
void AliTPCMonitor::WriteHistos()  
{
  // Writes all available histograms to a file in current directory
  // File name will be specified by : sector, side, runid and eventnumber
  
  if(GetEventProcessed())
  {
    AliInfo("Write histos to file");
    TString name=Form("SIDE_%i_SECTOR_%02i_RUN_%05i_EventID_%06i.root",(GetLastSector()/18),(GetLastSector()%18),fRunId,fEventNumber);
    TFile* f = new TFile(name.Data(),"recreate");
    for(Int_t i =0; i<fHistList->GetEntries(); i++)
    {
      if(((TH1*)fHistList->At(i))!=0)
      {
        ((TH1*)fHistList->At(i))->Write();
      }
    }
    f->ls();
    f->Close();
  }
  else
  {
    AliError("No Event Processed : Chose Format , File and push 'Next Event' ");
  }
}


//__________________________________________________________________
TH1* AliTPCMonitor::GetHisto(char* histname) 
{
  
  // Returns histogram specified by histname
  // check available names for histos in CreateHistos()
  
  TH1* hist = 0;
  if((TH1*)fHistList->FindObject(histname))
  {
    hist = (TH1*)fHistList->FindObject(histname);
  }
  else
  {
    cout << " AliTPCMonitor::GetHisto :: Can not find histo with name " << histname << endl;
  }
  return hist ;
}

//_________________________________
void AliTPCMonitor::SetupMonitoringTable(const char* table)
{
  //
  // Setup the monitoring table
  //
  fMonTableChanged=kFALSE;
  TString newTable(table);

  //check if we have a new table
  if (newTable==fMonTableString) return;
  fMonTableString=newTable;
  
  //delete old data
  delete fkMonTable;
  delete fMonTableArray;
  fkMonTable=0x0;

  //parse table request
  fMonTableArray=fMonTableString.Tokenize(",;");
  fMonTableArray->SetOwner();

  //consistency check
  Int_t entries=fMonTableArray->GetEntries();
  if (entries%4){
    AliError(Form("Monitoring table has to be given in muliples of 4\nEntries need to be sparated by ',' or ';'\nCannot parse the current table request: %s",fMonTableString.Data()));
    delete fMonTableArray;
    fMonTableArray=0x0;
    return;
  }

  fkMonTable=new const Char_t*[entries+1];
  for (Int_t i=0;i<entries;++i) fkMonTable[i]=fMonTableArray->At(i)->GetName();
  fkMonTable[entries]=0x0;

  fMonTableChanged=kTRUE;
}
