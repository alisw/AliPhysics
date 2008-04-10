/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

/** @file    AliFMDGainDA.cxx
    @author  Hans Hjersing Dalsgaard <canute@nbi.dk>
    @date    Mon Mar 13 13:46:05 2008
    @brief   Derived class for the pulse gain detector algorithm.
*/
// This class contains the implementation of the gain detector algorithms (DA) for the FMD.
// The data is collected in histograms that are reset for each pulse length after the 
// mean and standard deviation are put into a TGraphErrors object. After a certain number of pulses
// (usually 8) the graph is fitted to a straight line. The gain is then slope of this line as
// it combines the known pulse and the response of the detector.


#include "AliFMDGainDA.h"
#include "iostream"
#include "fstream"
#include "AliLog.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "AliFMDParameters.h"

//_____________________________________________________________________
ClassImp(AliFMDGainDA)

//_____________________________________________________________________
AliFMDGainDA::AliFMDGainDA() : AliFMDBaseDA(),
  fGainArray(),
  fPulseSize(32),
  fHighPulse(256), 
  fPulseLength(100),
  fEventsPerChannel(0),
  fCurrentPulse(0),
  fCurrentChannel(0),
  fNumberOfStripsPerChip(128)
{
  fOutputFile.open("gains.csv");
  fGainArray.SetOwner();
  
}

//_____________________________________________________________________
AliFMDGainDA::AliFMDGainDA(const AliFMDGainDA & gainDA) :  
  AliFMDBaseDA(gainDA),
  fGainArray(gainDA.fGainArray),
  fPulseSize(gainDA.fPulseSize),
  fHighPulse(gainDA.fHighPulse),
  fPulseLength(gainDA.fPulseLength),
  fEventsPerChannel(gainDA.fEventsPerChannel),
  fCurrentPulse(gainDA.fCurrentPulse),
  fCurrentChannel(gainDA.fCurrentChannel),
  fNumberOfStripsPerChip(gainDA.fNumberOfStripsPerChip)
{
  
}

//_____________________________________________________________________
AliFMDGainDA::~AliFMDGainDA() {
  
 

}

//_____________________________________________________________________
void AliFMDGainDA::Init() {
    
  fEventsPerChannel = (fPulseLength*fHighPulse) / fPulseSize ;
  
  SetRequiredEvents(fEventsPerChannel*fNumberOfStripsPerChip); //8 pulser values * 128 strips * 100 samples
  TObjArray* detArray;
  TObjArray* ringArray;
  TObjArray* sectorArray;
  
  for(UShort_t det=1;det<=3;det++) {
    detArray = new TObjArray();
    detArray->SetOwner();
    fGainArray.AddAtAndExpand(detArray,det);
    for (UShort_t ir = 0; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 2 : 4);
      ringArray = new TObjArray();
      ringArray->SetOwner();
      detArray->AddAtAndExpand(ringArray,ir);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	sectorArray = new TObjArray();
	sectorArray->SetOwner();
	ringArray->AddAtAndExpand(sectorArray,sec);
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  TH1S* hChannel = new TH1S(Form("hFMD%d%c_%d_%d",det,ring,sec,strip),Form("hFMD%d%c_%d_%d",det,ring,sec,strip),1024,0,1023);
	  hChannel->SetDirectory(0);
	  sectorArray->AddAtAndExpand(hChannel,strip);
	}
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDGainDA::AddChannelContainer(TObjArray* sectorArray, 
				       UShort_t det, 
				       Char_t ring, 
				       UShort_t sec, 
				       UShort_t strip) {
  
  TGraphErrors* hChannel  = new TGraphErrors();
  sectorArray->AddAtAndExpand(hChannel,strip);
}

//_____________________________________________________________________
void AliFMDGainDA::FillChannels(AliFMDDigit* digit) {

  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();
  
  
  if(strip%fNumberOfStripsPerChip) return;
  Int_t VAchip = strip / fNumberOfStripsPerChip; 
  TH1S* hChannel = GetChannelHistogram(det,ring,sec,VAchip);
  hChannel->Fill(digit->Counts());
  UpdatePulseAndADC(det,ring,sec,strip);
}

//_____________________________________________________________________
void AliFMDGainDA::Analyse(UShort_t det, 
			   Char_t ring, 
			   UShort_t sec, 
			   UShort_t strip) {
  TGraphErrors* grChannel = GetChannel(det,ring,sec,strip);
  if(!grChannel->GetN()) {
    // AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",det, ring , sec, strip));
    return;
  }
  TF1 fitFunc("fitFunc","pol1",-10,280); 
  fitFunc.SetParameters(100,3);
  grChannel->Fit("fitFunc","Q","Q",0,fHighPulse);
  AliFMDParameters* pars = AliFMDParameters::Instance();
  UInt_t ddl, board,chip,channel;
  pars->Detector2Hardware(det,ring,sec,strip,ddl,board,chip,channel);
    
  Float_t chi2ndf = 0;
  if(fitFunc.GetNDF())
    chi2ndf = fitFunc.GetChisquare() / fitFunc.GetNDF();
  ddl = ddl + kBaseDDL;
  
  Int_t relStrip = strip%fNumberOfStripsPerChip;
  
  fOutputFile << ddl                         << ','
	      << board                       << ','
	      << chip                        << ','
	      << channel                     << ','
	      << relStrip                    << ','
	      << fitFunc.GetParameter(1)     << ','
	      << fitFunc.GetParError(1)      << ','
	      << chi2ndf                     <<"\n";
  
  
  if(fSaveHistograms) {
    
    gDirectory->cd(Form("%s:FMD%d%c/sector_%d/strip_%d",fDiagnosticsFilename,det,ring,sec,strip));
    grChannel->Write(Form("grFMD%d%c_%d_%d",det,ring,sec,strip));
  }
  
 
  
}

//_____________________________________________________________________
void AliFMDGainDA::WriteHeaderToFile() {
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fOutputFile.write(Form("# %s \n",pars->GetGainShuttleID()),9);
  fOutputFile.write("# Rcu, Board, Chip, Channel, Strip, Gain, Error, Chi2/NDF \n",59);
  
}

//_____________________________________________________________________
TH1S* AliFMDGainDA::GetChannelHistogram(UShort_t det, 
					Char_t ring, 
					UShort_t sec, 
					UShort_t strip) {
  
  UShort_t  Ring = 1;
  if(ring == 'O')
    Ring = 0;
  
  
  TObjArray* detArray  = static_cast<TObjArray*>(fGainArray.At(det));
  TObjArray* ringArray = static_cast<TObjArray*>(detArray->At(Ring));
  TObjArray* secArray  = static_cast<TObjArray*>(ringArray->At(sec));
  TH1S* hChannel       = static_cast<TH1S*>(secArray->At(strip));
  
  return hChannel;
}

//_____________________________________________________________________
TGraphErrors* AliFMDGainDA::GetChannel(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip) {
  
  UShort_t  Ring = 1;
  if(ring == 'O')
    Ring = 0;
  
  
  TObjArray* detArray          = static_cast<TObjArray*>(fDetectorArray.At(det));
  TObjArray* ringArray         = static_cast<TObjArray*>(detArray->At(Ring));
  TObjArray* secArray          = static_cast<TObjArray*>(ringArray->At(sec));
  TGraphErrors* hChannel       = static_cast<TGraphErrors*>(secArray->At(strip));
  
  return hChannel;
}

//_____________________________________________________________________
void AliFMDGainDA::UpdatePulseAndADC(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip) {
  
  if(strip%fNumberOfStripsPerChip) return;
  if(((GetCurrentEvent())%fPulseLength) && GetCurrentEvent()>0) return;
    
  Int_t VAchip = strip/fNumberOfStripsPerChip; 
  TH1S* hChannel = GetChannelHistogram(det,ring,sec,VAchip);
  
  if(!hChannel->GetEntries()) {
    AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",det, ring , sec, strip));
    return;
  }
  Double_t mean      = hChannel->GetMean();
  Double_t rms       = hChannel->GetRMS();
  Double_t pulse     = (Double_t)fCurrentPulse*fPulseSize;
  
  hChannel->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
  
  mean      = hChannel->GetMean();
  rms       = hChannel->GetRMS();
    
  Int_t channelNumber = strip + (GetCurrentEvent()-1)/800; 
  
  TGraphErrors* channel = GetChannel(det,ring,sec,channelNumber);
  
  channel->SetPoint(fCurrentPulse,pulse,mean);
  channel->SetPointError(fCurrentPulse,0,rms);
  
  if(fSaveHistograms) {
    gDirectory->cd(Form("%s:FMD%d%c/sector_%d/strip_%d",fDiagnosticsFilename,det,ring,sec,channelNumber));
    hChannel->GetXaxis()->SetRange(0,1023);
    hChannel->Write(Form("hFMD%d%c_%d_%d_pulse_%d",det,ring,sec,channelNumber,fCurrentPulse));
  }
  
  hChannel->Reset();
  
}

//_____________________________________________________________________
void AliFMDGainDA::ResetPulseAndUpdateChannel() {
  
  fCurrentPulse   = 0;
  
}

//_____________________________________________________________________
void AliFMDGainDA::FinishEvent() {
  
  if(GetCurrentEvent()>0 && (GetCurrentEvent()%fPulseLength==0))
    fCurrentPulse++;
  
  if(GetCurrentEvent()>0 && (GetCurrentEvent())%fEventsPerChannel==0)
    fCurrentPulse = 0;


}
//_____________________________________________________________________
//
//EOF
//
