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
//____________________________________________________________________
//
// This class contains the implementation of the gain detector
// algorithms (DA) for the FMD.  The data is collected in histograms
// that are reset for each pulse length after the mean and standard
// deviation are put into a TGraphErrors object. After a certain
// number of pulses (usually 8) the graph is fitted to a straight
// line. The gain is then slope of this line as it combines the known
// pulse and the response of the detector.
//
#include "AliFMDGainDA.h"
#include "AliFMDAltroMapping.h"
#include "AliFMDParameters.h"
#include "AliFMDCalibGain.h"
#include "AliFMDDigit.h"
#include "AliLog.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1S.h>
#include <TGraphErrors.h>
#include <TDatime.h>
#include <TH2.h>
#include <TROOT.h>

//_____________________________________________________________________
ClassImp(AliFMDGainDA)
#if 0 // Do not delete - here to let Emacs indent properly
;
#endif

//_____________________________________________________________________
AliFMDGainDA::AliFMDGainDA() 
  : AliFMDBaseDA(),
    fHighPulse(256), 
    fEventsPerChannel(10),
    fCurrentPulse(10),
    fCurrentChannel(10),
    fNumberOfStripsPerChip(128),
    fSummaryGains("GainsSummary","Summary of gains",51200,0,51200),
    fCurrentSummaryStrip(1),
    fGainFMD1i(0),
    fGainFMD2i(0),
    fGainFMD2o(0),
    fGainFMD3i(0),
    fGainFMD3o(0),
    fChi2FMD1i(0),
    fChi2FMD2i(0),
    fChi2FMD2o(0),
    fChi2FMD3i(0),
    fChi2FMD3o(0)
{
  // Constructor 
  // 
  // Parameters: 
  //   None
  fCurrentPulse.Reset(0);
  fCurrentChannel.Reset(0);
  fDiagnosticsFilename = "diagnosticsGain.root";
}

//_____________________________________________________________________
AliFMDGainDA::AliFMDGainDA(const AliFMDGainDA & gainDA) 
  : AliFMDBaseDA(gainDA),
    fHighPulse(gainDA.fHighPulse),
    fEventsPerChannel(gainDA.fEventsPerChannel),
    fCurrentPulse(gainDA.fCurrentPulse),
    fCurrentChannel(gainDA.fCurrentChannel),
    fNumberOfStripsPerChip(gainDA.fNumberOfStripsPerChip),
    fSummaryGains(gainDA.fSummaryGains),
    fCurrentSummaryStrip(gainDA.fCurrentSummaryStrip),
    fGainFMD1i(gainDA.fGainFMD1i),
    fGainFMD2i(gainDA.fGainFMD2i),
    fGainFMD2o(gainDA.fGainFMD2o),
    fGainFMD3i(gainDA.fGainFMD3i),
    fGainFMD3o(gainDA.fGainFMD3o),
    fChi2FMD1i(gainDA.fChi2FMD1i),
    fChi2FMD2i(gainDA.fChi2FMD2i),
    fChi2FMD2o(gainDA.fChi2FMD2o),
    fChi2FMD3i(gainDA.fChi2FMD3i),
    fChi2FMD3o(gainDA.fChi2FMD3o)
{  
  // Copy Constructor 
  // 
  // Parameters: 
  //   gainDA   Object to copy from 
  fCurrentPulse.Reset(0);
  fCurrentChannel.Reset(0);
}

//_____________________________________________________________________
AliFMDGainDA::~AliFMDGainDA() 
{
  // Destructor 
  // 
  // Parameters: 
  //   None
}


//_____________________________________________________________________
Bool_t
AliFMDGainDA::OpenFiles(Bool_t appendRun) 
{
  if (!AliFMDBaseDA::OpenFiles(appendRun)) return false;
  if (!appendRun || fRunno == 0) { 
    Rotate("gains.csv", 3);
    fOutputFile.open("gains.csv");
  }
  else 
    fOutputFile.open(Form("gains_%09d.csv", fRunno));
  if (!fOutputFile) { 
    Error("OpenFiles", "Failed to open gains file");
    return false;
  }
  return true;
}

//_____________________________________________________________________
void AliFMDGainDA::Init() 
{
  // Initialize 
  // 
  // Parameters: 
  //   None
  Int_t nEventsRequired = 0;
  
  for(Int_t idx = 0;idx<fEventsPerChannel.GetSize();idx++) {
    Int_t nEvents = 0;
    if(fPulseSize.At(idx))
      nEvents = (fPulseLength.At(idx)*fHighPulse) / fPulseSize.At(idx);
    fEventsPerChannel.AddAt(nEvents,idx);
    if(nEvents>nEventsRequired) 
      nEventsRequired = nEvents * fNumberOfStripsPerChip;    
  }
  SetRequiredEvents(nEventsRequired); 
}

//_____________________________________________________________________
void AliFMDGainDA::InitContainer(TDirectory* dir) 
{
  AliFMDBaseDA::InitContainer(dir);
  
  for(UShort_t det=1;det<=3;det++) {
    UShort_t sr = (det == 1 ? 1 : 0);
    for (UShort_t ir = sr; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nva  = (ir == 0 ? 2 : 4);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	Array* sectorArray = GetSectorArray(det, ring, sec);
	Array* cache = new Array(nva);
	cache->SetName("Cache");
	cache->SetOwner();
	Int_t n = sectorArray->GetEntriesFast();
	sectorArray->AddAtAndExpand(cache, n);
	for(UShort_t va = 0; va < nva; va++) {
	  TH1S* hChannel = new TH1S(Form("FMD%d%c[%02d]_va%d",
					 det,ring,sec,va),
				    Form("FMD%d%c[%02d] VA%d Cache",
					 det,ring,sec,va),
				    1024,-.5,1023.5);
	  hChannel->SetDirectory(0);
	  cache->AddAtAndExpand(hChannel,va);
	}
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDGainDA::AddChannelContainer(Array* stripArray, 
				       UShort_t det  , 
				       Char_t   ring,  
				       UShort_t sec, 
				       UShort_t strip) 
{  
  // Make a channel container 
  // 
  // Parameters: 
  //  sectorArray    Sectors 
  //  det            Detector number
  //  ring           Ring identifier 
  //  sec            Sector number
  //  strip          Strip number

  AliFMDParameters* pars     = AliFMDParameters::Instance();
  UShort_t          board    = pars->GetAltroMap()->Sector2Board(ring, sec);
  Int_t             halfring = GetHalfringIndex(det,ring,board/16);
  Int_t             dPulse   = fPulseSize.At(halfring);
  Int_t             nPulses  = dPulse > 0 ? fHighPulse / dPulse : 0;

  TGraphErrors* gChannel  = new TGraphErrors(nPulses);
  gChannel->SetName(Form("FMD%d%c[%02d,%03d]", det, ring, sec, strip));
  gChannel->SetTitle(Form("FMD%d%c[%02d,%03d] ADC vs DAC", 
			  det, ring, sec, strip));
  stripArray->AddAtAndExpand(gChannel,0);
}

//_____________________________________________________________________
void AliFMDGainDA::AddSectorSummary(Array* sectorArray, 
				    UShort_t det, 
				    Char_t   ring, 
				    UShort_t sec, 
				    UShort_t nStr) 
{
  TH1F* summary = new TH1F("Summary", Form("Summary of gains in FMD%d%c[%02d]", 
					   det, ring, sec), 
			   nStr, -.5, nStr-.5);
  summary->SetXTitle("Strip");
  summary->SetYTitle("Gain [ADC/DAC]");
  summary->SetDirectory(0);

  Int_t n = sectorArray->GetEntriesFast();
  sectorArray->AddAtAndExpand(summary, n);
}

//_____________________________________________________________________
void AliFMDGainDA::FillChannels(AliFMDDigit* digit) 
{
  // Fill data into histogram
  // 
  // Parameters: 
  //   digit  Digit to get the data from

  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();
  
  //Strip is always seen as the first in a VA chip. All other strips are junk.
  //Strips are counted from zero on even sectors and from 511 on odd sectors...
   
  if((sec%2)     && ((strip+1) % fNumberOfStripsPerChip)) return;
  if(((sec+1)%2) && (strip % fNumberOfStripsPerChip)) return;
  
  UShort_t vaChip   = strip / fNumberOfStripsPerChip; 
  TH1S*    hChannel = GetChannelHistogram(det, ring, sec, vaChip);
  hChannel->Fill(digit->Counts());
  UpdatePulseAndADC(det,ring,sec,strip);
}

//_____________________________________________________________________
void AliFMDGainDA::MakeSummary(UShort_t det, Char_t ring)
{
  // 
  // Create summary hists for FMD gains and chi2 of the fits
  //
  switch (det) { 
  case 1: 
    fGainFMD1i = MakeSummaryHistogram("gain", "Gains", det, ring);
    fChi2FMD1i = MakeSummaryHistogram("chi2", "#Chi^{2}/NDF", det, ring);
    break;
  case 2:
    switch (ring) { 
    case 'I': case 'i':
      fGainFMD2i = MakeSummaryHistogram("gain", "Gains", det, ring);
      fChi2FMD2i = MakeSummaryHistogram("chi2", "#Chi^{2}/NDF", det, ring);
      break;
    case 'O': case 'o':
      fGainFMD2o = MakeSummaryHistogram("gain", "Gains", det, ring);
      fChi2FMD2o = MakeSummaryHistogram("chi2", "#Chi^{2}/NDF", det, ring);
      break;
    }
    break;
  case 3:
    switch (ring) { 
    case 'I': case 'i':
      fGainFMD3i = MakeSummaryHistogram("gain", "Gains", det, ring);
      fChi2FMD3i = MakeSummaryHistogram("chi2", "#Chi^{2}/NDF", det, ring);
      break;
    case 'O': case 'o':
      fGainFMD3o = MakeSummaryHistogram("gain", "Gains", det, ring);
      fChi2FMD3o = MakeSummaryHistogram("chi2", "#Chi^{2}/NDF", det, ring);
      break;
    }
    break;
  }
}

//_____________________________________________________________________
void AliFMDGainDA::Analyse(UShort_t det, 
			   Char_t   ring, 
			   UShort_t sec, 
			   UShort_t strip) 
{
  // Analyse result of a single strip
  // 
  // Parameters: 
  //  det            Detector number
  //  ring           Ring identifier 
  //  sec            Sector number
  //  strip          Strip number
  TGraphErrors* grChannel = GetChannel(det,ring,sec,strip);
  if(!grChannel->GetN()) {
    AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",
                     det, ring , sec, strip));
    return;
  }
  TF1* fitFunc = new TF1("fitFunc","pol1",-10,280); 
  fitFunc->SetParameters(100,3);
  grChannel->Fit(fitFunc,"Q0","",0,fHighPulse);
  if (grChannel->GetListOfFunctions())
    grChannel->GetListOfFunctions()->Remove(fitFunc);
  gROOT->GetListOfFunctions()->Remove(fitFunc);

     
  Float_t gain    = -1;
  Float_t error   = -1; 
  Float_t chi2ndf = -1;
  if((fitFunc->GetParameter(1)) == (fitFunc->GetParameter(1))) {
    gain    = fitFunc->GetParameter(1);
    error   = fitFunc->GetParError(1);
    if(fitFunc->GetNDF())
      chi2ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();
  }

  fOutputFile << det                         << ','
	      << ring                        << ','
	      << sec                         << ','
	      << strip                       << ','
	      << gain                        << ','
	      << error                       << ','
	      << chi2ndf                     <<"\n";
  
  //due to RCU trouble, first strips on VAs are excluded
  // if(strip%128 != 0) {
    
  fSummaryGains.SetBinContent(fCurrentSummaryStrip,fitFunc->GetParameter(1));
  fSummaryGains.SetBinError(fCurrentSummaryStrip,fitFunc->GetParError(1));
  
  fCurrentSummaryStrip++;

  TH2* hGain = 0;
  TH2* hChi2 = 0;
  switch (det) { 
  case 1: hGain = fGainFMD1i; hChi2 = fChi2FMD1i; break;
  case 2: 
    switch (ring) { 
    case 'I':  hGain = fGainFMD2i; hChi2 = fChi2FMD2i; break;
    case 'O':  hGain = fGainFMD2o; hChi2 = fChi2FMD2o; break;
    }
    break;
  case 3:
    switch (ring) { 
    case 'I':  hGain = fGainFMD3i; hChi2 = fChi2FMD3i; break;
    case 'O':  hGain = fGainFMD3o; hChi2 = fChi2FMD3o; break;
    }
    break;
  }
  if (hGain && hChi2) {
    Int_t bin = hGain->FindBin(sec, strip);
    hGain->SetBinContent(bin, gain);
    hGain->SetBinError(bin, error);
    hChi2->SetBinContent(bin, chi2ndf);
  }

    // }
  if(fSaveHistograms) {
    TH1F* summary = GetSectorSummary(det, ring, sec);
    summary->SetBinContent(strip+1, fitFunc->GetParameter(1));
    summary->SetBinError(strip+1, fitFunc->GetParError(1));
  }  
  delete fitFunc;
}

//_____________________________________________________________________
void AliFMDGainDA::Terminate(TFile* diagFile)
{
  // End of file 
  // 
  // Parameters: 
  //   None
  if(diagFile) {
    diagFile->cd();
    fSummaryGains.Write();
  }
}

//_____________________________________________________________________
void AliFMDGainDA::WriteHeaderToFile() 
{
  // Write header to the output file
  // 
  // Parameters: 
  //   None
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fOutputFile.write(Form("# %s \n",pars->GetGainShuttleID()),9);
  TDatime now;
  fOutputFile << "# This file created from run # " << fRunno 
	      << " @ " << now.AsString() << std::endl;
  fOutputFile.write("# Detector, "
		    "Ring, "
		    "Sector, "
		    "Strip, "
		    "Gain, "
		    "Error, "
		    "Chi2/NDF \n",56);  
}

//_____________________________________________________________________
TH1S* AliFMDGainDA::GetChannelHistogram(UShort_t det, 
					Char_t   ring, 
					UShort_t sec, 
					UShort_t va) 
{
  // Get the current histogram of a single strip
  // 
  // Parameters: 
  //  det            Detector number
  //  ring           Ring identifier 
  //  sec            Sector number
  //  va             VA number
  Array* secArray  = GetSectorArray(det, ring, sec);
  Int_t      n         = secArray->GetEntriesFast();
  Array* cache     = static_cast<Array*>(secArray->At(n-1));
  TH1S* hChannel       = static_cast<TH1S*>(cache->At(va));
  
  return hChannel;
}

//_____________________________________________________________________
TGraphErrors* AliFMDGainDA::GetChannel(UShort_t det, 
				       Char_t   ring, 
				       UShort_t sec, 
				       UShort_t strip) 
{  
  // Get the graph of a single strip
  // 
  // Parameters: 
  //  det            Detector number
  //  ring           Ring identifier 
  //  sec            Sector number
  //  strip          Strip number
  Array*    stripArray = GetStripArray(det, ring, sec, strip);
  TGraphErrors* hChannel   = static_cast<TGraphErrors*>(stripArray->At(0));
  
  return hChannel;
}

//_____________________________________________________________________
TH1F* AliFMDGainDA::GetSectorSummary(UShort_t det, 
				     Char_t   ring, 
				     UShort_t sec) 
{
  Array* secArray    = GetSectorArray(det, ring, sec);
  Int_t      n           = secArray->GetEntriesFast();
  return static_cast<TH1F*>(secArray->At(n-2)); // Cache added later
}

//_____________________________________________________________________
void AliFMDGainDA::UpdatePulseAndADC(UShort_t det, 
				     Char_t ring, 
				     UShort_t sec, 
				     UShort_t strip) 
{
  // Go to next pulse size
  // 
  // Parameters: 
  //  det            Detector number
  //  ring           Ring identifier 
  //  sec            Sector number
  //  strip          Strip number
  
  AliFMDParameters* pars     = AliFMDParameters::Instance();
  UShort_t          board    = pars->GetAltroMap()->Sector2Board(ring, sec);
  Int_t             halfring = GetHalfringIndex(det,ring,board/16);
  
  if(GetCurrentEvent()> (fNumberOfStripsPerChip*fEventsPerChannel.At(halfring)))
    return;
  
  if ((sec       % 2) && ((strip+1) % fNumberOfStripsPerChip)) return;
  if (((sec + 1) % 2) && (strip     % fNumberOfStripsPerChip)) return;
  
  if(((GetCurrentEvent()) % fPulseLength.At(halfring)) 
     && GetCurrentEvent() > 0) return;
     
  Int_t vaChip   = strip/fNumberOfStripsPerChip; 
  TH1S* hChannel = GetChannelHistogram(det,ring,sec,vaChip);
  
  if(!hChannel->GetEntries()) {
    AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",
		    det, ring , sec, strip));
    return;
  }
  Double_t mean      = hChannel->GetMean();
  Double_t rms       = hChannel->GetRMS();
  Double_t pulse     = (Double_t(fCurrentPulse.At(halfring)) 
			* fPulseSize.At(halfring));
  Int_t    firstBin  = hChannel->GetXaxis()->GetFirst();
  Int_t    lastBin   = hChannel->GetXaxis()->GetLast();
  hChannel->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
  
  mean               = hChannel->GetMean();
  rms                = hChannel->GetRMS();
  
  hChannel->GetXaxis()->SetRange(firstBin,lastBin);
  
  Int_t         channelOff = ((GetCurrentEvent()-1)
			      / ((fPulseLength.At(halfring)*fHighPulse)
				 / fPulseSize.At(halfring)));
  Int_t         channelNo  = (strip + ((sec % 2 == 1) ? -1 : 1) * channelOff);
  TGraphErrors* channel    = GetChannel(det,ring,sec,channelNo);

  channel->SetPoint(fCurrentPulse.At(halfring),pulse,mean);
  channel->SetPointError(fCurrentPulse.At(halfring),0,rms);
  
  if(fSaveHistograms) {
    TH1S* out = 
      static_cast<TH1S*>(hChannel->Clone(Form("FMD%d%c[%02d,%03d]_0x%02x",
					      det, ring, sec, channelNo, 
					      int(pulse))));
    out->SetTitle(Form("FMD%d%c[%02d,%03d] DAC=0x%02x (%3d)",
		       det, ring, sec, channelNo, int(pulse), int(pulse)));
    out->SetDirectory(0);
    Array* arr = GetStripArray(det, ring, sec, channelNo);
    arr->AddAtAndExpand(out, fCurrentPulse.At(halfring)+1);
  }
    
  hChannel->Reset();
  
}

//_____________________________________________________________________
void AliFMDGainDA::ResetPulseAndUpdateChannel() 
{  
  // Reset all
  // 
  // Parameters: 
  //   None
  fCurrentPulse.Reset(0); 
}

//_____________________________________________________________________
void AliFMDGainDA::FinishEvent() 
{
  // End of event
  // 
  // Parameters: 
  //   None
  for(UShort_t det=1; det<=3;det++) {
    UShort_t firstring = (det == 1 ? 1 : 0);
    for(UShort_t iring = firstring; iring <=1;iring++) {
      Char_t ring = (iring == 1 ? 'I' : 'O');
      for(UShort_t board =0 ; board <=1; board++) {
	Int_t idx = GetHalfringIndex(det,ring,board);
	
	if(!fPulseLength.At(idx) || !fEventsPerChannel.At(idx))
	  continue;
	if(GetCurrentEvent()>0 && 
	   ((GetCurrentEvent() % fPulseLength.At(idx)) == 0))
	  fCurrentPulse.AddAt(fCurrentPulse.At(idx)+1,idx);
	
	if(GetCurrentEvent()>0 && 
	   ((GetCurrentEvent()) % fEventsPerChannel.At(idx)) == 0)
	  fCurrentPulse.AddAt(0,idx);
      }
    }
  }
}
//_____________________________________________________________________
//
//EOF
//
