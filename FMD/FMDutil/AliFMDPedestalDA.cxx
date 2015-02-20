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

/** @file    AliMDPedestalDA.cxx
    @author  Hans Hjersing Dalsgaard <canute@nbi.dk>
    @date    Mon Mar 10 09:46:05 2008
    @brief   Derived class for the pedestal detector algorithm.
*/
//
// This class implements the virtual functions of the AliFMDBaseDA
// class.  The most important of these functions, FillChannels(..) and
// Analyse(...) collect and analyse the data of each channel. The
// resulting pedestal and noise values are written to a comma
// separated values (csv) file on the go. The csv files produced in
// this way are the basic input to the AliFMDPreprocessor.
//

#include "AliFMDPedestalDA.h"
#include "AliFMDAltroMapping.h"
#include "AliFMDParameters.h"
#include "AliFMDCalibPedestal.h"
#include "AliFMDDigit.h"
#include "AliLog.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TFile.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TDatime.h>
#include <TH2.h>
#include <TROOT.h>

//_____________________________________________________________________
ClassImp(AliFMDPedestalDA)

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA() 
  : AliFMDBaseDA(),
    fCurrentChannel(1),
    fPedSummary("PedestalSummary","pedestals",51200,0,51200),
    fNoiseSummary("NoiseSummary","noise",51200,0,51200),
    fZSfileFMD1(),
    fZSfileFMD2(),
    fZSfileFMD3(), 
    fMinTimebin(3 * 4 * 3 * 16), // 3 ddls, 4 FECs, 3 Altros, 16 channels
    fMaxTimebin(3 * 4 * 3 * 16), // 3 ddls, 4 FECs, 3 Altros, 16 channels
    fSummaryFMD1i(0),
    fSummaryFMD2i(0),
    fSummaryFMD2o(0),
    fSummaryFMD3i(0),
    fSummaryFMD3o(0)
{
  // Default constructor 
  fDiagnosticsFilename = "diagnosticsPedestal.root";
}

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA(const AliFMDPedestalDA & pedDA) : 
  AliFMDBaseDA(pedDA),
  fCurrentChannel(1),
  fPedSummary("PedestalSummary","pedestals",51200,0,51200),
  fNoiseSummary("NoiseSummary","noise",51200,0,51200),
  fZSfileFMD1(),
  fZSfileFMD2(),
  fZSfileFMD3(),
  fMinTimebin(pedDA.fMinTimebin),
  fMaxTimebin(pedDA.fMaxTimebin),
  fSummaryFMD1i(pedDA.fSummaryFMD1i),
  fSummaryFMD2i(pedDA.fSummaryFMD2i),
  fSummaryFMD2o(pedDA.fSummaryFMD2o),
  fSummaryFMD3i(pedDA.fSummaryFMD3i),
  fSummaryFMD3o(pedDA.fSummaryFMD3o)
{
  // Copy constructor 
}

//_____________________________________________________________________
AliFMDPedestalDA::~AliFMDPedestalDA() 
{ 
  // Destructor.
}

//_____________________________________________________________________
Bool_t 
AliFMDPedestalDA::OpenFiles(Bool_t appendRun)
{
  if (!AliFMDBaseDA::OpenFiles(appendRun)) return false;
  if (!appendRun || fRunno == 0) {
    Rotate("peds.csv", 3);
    fOutputFile.open("peds.csv");
  }
  else 
    fOutputFile.open(Form("peds_%09d.csv",fRunno));
  
  if (!fOutputFile) { 
    Error("OpenFiles", "Failed to open pedestal file");
    return false;
  }
  Rotate("ddl3072.csv", 10);
  fZSfileFMD1.open("ddl3072.csv");
  Rotate("ddl3073.csv", 10);
  fZSfileFMD2.open("ddl3073.csv");
  Rotate("ddl3074.csv", 10);
  fZSfileFMD3.open("ddl3074.csv");  
  return true;
}

//_____________________________________________________________________
void AliFMDPedestalDA::Init() 
{ 
  // Initialise 
  SetRequiredEvents(1000);
  fMinTimebin.Reset(1024);
  fMaxTimebin.Reset(-1);


}

//_____________________________________________________________________
void AliFMDPedestalDA::AddChannelContainer(Array* sampleArray, 
					   UShort_t det, 
					   Char_t   ring, 
					   UShort_t sec, 
					   UShort_t strip) 
{
  // Add a channel to the containers.
  //
  // Parameters:
  //     sectorArray  Array of sectors
  //     det          Detector 
  //     ring         Ring
  //     sec          Sector 
  //     strip        Strip
  AliFMDParameters* pars        = AliFMDParameters::Instance();
  UInt_t            samples     = pars->GetSampleRate(det, ring, sec, strip);
  for (UInt_t sample = 0; sample < samples; sample++) {
    TString name(Form("FMD%d%c[%02d,%03d]_%d", det,ring,sec,strip,sample));
    TH1S* hSample = new TH1S(name.Data(),name.Data(), 1024,-.5,1023.5);
    hSample->SetXTitle("ADC");
    hSample->SetYTitle("Events");
    hSample->SetDirectory(0);
    hSample->ResetBit(TObject::kMustCleanup);
    sampleArray->AddAtAndExpand(hSample, sample);
  }
}

//_____________________________________________________________________
void AliFMDPedestalDA::AddSectorSummary(Array* sectorArray, 
					UShort_t det, 
					Char_t   ring, 
					UShort_t sec, 
					UShort_t nStr) 
{
  TH1F* sumPed = new TH1F("Pedestals", 
			  Form("Summary of pedestals in FMD%d%c[%02d]", 
			       det, ring, sec), 
			  nStr, -.5, nStr-.5);
  sumPed->SetXTitle("Strip");
  sumPed->SetYTitle("Pedestal [ADC]");
  sumPed->SetDirectory(0);
  sumPed->ResetBit(TObject::kMustCleanup);
  
  TH1F* sumNoise = static_cast<TH1F*>(sumPed->Clone("Noise"));
  sumNoise->SetYTitle("Noise [ADC]");
  sumNoise->SetDirectory(0);
  sumNoise->ResetBit(TObject::kMustCleanup);
  
  Int_t n = sectorArray->GetEntriesFast();
  sectorArray->AddAtAndExpand(sumPed,   n + kPedestalOffset - 1);
  sectorArray->AddAtAndExpand(sumNoise, n + kNoiseOffset - 1);
}

//_____________________________________________________________________
void AliFMDPedestalDA::FillChannels(AliFMDDigit* digit) 
{
  // Fill ADC values from a digit into the corresponding histogram.
  //
  // Parameters:
  //    digit Digit to fill ADC values for.
  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();
  
  AliFMDParameters* pars     = AliFMDParameters::Instance();
  UInt_t            samples  = pars->GetSampleRate(det, ring, sec, strip);
  for (UInt_t sample = 0; sample < samples; sample++) {
    TH1S* hSample = GetChannel(det, ring, sec, strip, sample);
    if (!hSample) continue;
    
    hSample->Fill(digit->Count(sample));
  }
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::MakeSummary(UShort_t det, Char_t ring)
{
  //Create summary hists for FMD pedestals
  // std::cout << "Making summary for FMD" << det << ring << " ..." 
  //           << std::endl;
  switch (det) { 
  case 1: 
    fSummaryFMD1i = MakeSummaryHistogram("ped", "Pedestals", det, ring);
    break;
  case 2:
    switch (ring) { 
    case 'I': case 'i':
      fSummaryFMD2i = MakeSummaryHistogram("ped", "Pedestals", det, ring);
      break;
    case 'O': case 'o':
      fSummaryFMD2o = MakeSummaryHistogram("ped", "Pedestals", det, ring);
      break;
    }
    break;
  case 3:
    switch (ring) { 
    case 'I': case 'i':
      fSummaryFMD3i = MakeSummaryHistogram("ped", "Pedestals", det, ring);
      break;
    case 'O': case 'o':
      fSummaryFMD3o = MakeSummaryHistogram("ped", "Pedestals", det, ring);
      break;
    }
    break;
  }
}

//_____________________________________________________________________
void AliFMDPedestalDA::Analyse(UShort_t det, 
			       Char_t   ring, 
			       UShort_t sec, 
			       UShort_t strip) 
{
  // Analyse a strip.  That is, compute the mean and spread of the ADC
  // spectra for all strips.  Also output on files the values. 
  //
  // Parameters:
  //     det   Detector
  //     ring  Ring
  //     sec   Sector 
  //     strip Strip.
  AliFMDParameters* pars     = AliFMDParameters::Instance();
  TH2* summary = 0;
  switch (det) { 
  case 1: summary = fSummaryFMD1i; break;
  case 2: 
    switch (ring) { 
    case 'I':  summary = fSummaryFMD2i; break;
    case 'O':  summary = fSummaryFMD2o; break;
    }
    break;
  case 3:
    switch (ring) { 
    case 'I':  summary = fSummaryFMD3i; break;
    case 'O':  summary = fSummaryFMD3o; break;
    }
    break;
  }

  UInt_t samples  = pars->GetSampleRate(det, ring, sec, strip);
  for (UShort_t sample = 0; sample < samples; sample++) {
    TH1S* hChannel = GetChannel(det, ring, sec, strip,sample);
    if(!hChannel || hChannel->GetEntries() == 0) {
      //AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",
      //	    det,ring,sec,strip));
      return;
    }
    
    AliDebugF(50, "Fitting FMD%d%c[%02d,%03d] with %d entries",
	      det,ring,sec,strip, int(hChannel->GetEntries()));
    TF1* fitFunc = new TF1("fitFunc","gausn",0,300);
    fitFunc->ResetBit(TObject::kMustCleanup);
    fitFunc->SetParameters(100,100,1);
    hChannel->Fit(fitFunc,"Q0","",10,200);
    hChannel->GetListOfFunctions()->Remove(fitFunc);
    gROOT->GetListOfFunctions()->Remove(fitFunc);

    Float_t mean = hChannel->GetMean();
    Float_t rms  = hChannel->GetRMS();
    
    hChannel->GetXaxis()->SetRangeUser(mean-5*rms,mean+5*rms);
    mean = hChannel->GetMean();
    rms  = hChannel->GetRMS();
  
    
    UShort_t ddl, board, altro, channel;
    UShort_t timebin;
    
    pars->Detector2Hardware(det,ring,sec,strip,sample,
			    ddl,board,altro,channel,timebin);
    Int_t idx = HWIndex(ddl, board, altro, channel);
    if (idx >= 0) { 
      fMinTimebin[idx] = TMath::Min(Short_t(timebin),   fMinTimebin[idx]);
      fMaxTimebin[idx] = TMath::Max(Short_t(timebin+1), fMaxTimebin[idx]);
    }
    
    std::ostream* zsFile = 0;
    switch(det) {
    case 1:  zsFile = &fZSfileFMD1; break;
    case 2:  zsFile = &fZSfileFMD2; break;
    case 3:  zsFile = &fZSfileFMD3; break;
    default: AliWarning("Unknown sample!"); break;
      
    }
    *zsFile << board   << ',' 
	    << altro   << ',' 
	    << channel << ',' 
	    << timebin << ','
	    << mean    << ',' 
	    << rms     << "\n";
    
    Float_t chi2ndf = 0;
    
    
    if(fitFunc->GetNDF())
      chi2ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();
    
    
    Int_t sampleToWrite = 2;
    if      (samples == 2) sampleToWrite = 1;
    else if (samples <  2) sampleToWrite = 0;
    
    hChannel->GetXaxis()->SetRange(1,1024);

    if(sample != sampleToWrite) continue;
    
    
    fOutputFile << det                         << ','
		<< ring                        << ','
		<< sec                         << ','
		<< strip                       << ','
		<< mean                        << ','
		<< rms                         << ','
		<< fitFunc->GetParameter(1)    << ','
		<< fitFunc->GetParameter(2)    << ','
		<< chi2ndf                     <<"\n";

    delete fitFunc;

    if (summary) { 
      Int_t bin = summary->FindBin(sec, strip);
      summary->SetBinContent(bin, mean);
      summary->SetBinError(bin, rms);
    }

    if(fSaveHistograms  ) {
      TH1F* sumPed   = GetSectorSummary(det, ring, sec, true);
      TH1F* sumNoise = GetSectorSummary(det, ring, sec, false);
      sumPed->SetBinContent(strip+1, mean);
      sumPed->SetBinError(strip+1, rms);
      sumNoise->SetBinContent(strip+1, rms);

      fPedSummary.SetBinContent(fCurrentChannel,mean);
      fNoiseSummary.SetBinContent(fCurrentChannel,rms);
      fCurrentChannel++;
    }
  }
}

//_____________________________________________________________________
void AliFMDPedestalDA::Terminate(TFile* diagFile) 
{
  // Called at the end of a job.  Fills in missing time-bins and
  // closes output files
  if(fSaveHistograms && diagFile) {
    diagFile->cd();
    
    fPedSummary.Write();
    fNoiseSummary.Write();
  }
  AliFMDAltroMapping* map = AliFMDParameters::Instance()->GetAltroMap();
  for (Int_t i = 0; i < 3; i++) { 
    std::ofstream& out = (i == 0 ? fZSfileFMD1 : 
			  i == 1 ? fZSfileFMD2 : 
			  fZSfileFMD3);
    if (out.is_open() && fSeenDetectors[i]) { 
      FillinTimebins(out, map->Detector2DDL(i+1));
    }
    if (!fSeenDetectors[i]) {
      TString n(Form("ddl%d.csv",3072+map->Detector2DDL(i+1)));
      gSystem->Unlink(n.Data());
    }
  }
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::FillinTimebins(std::ofstream& out, UShort_t /*ddl*/)
{
  // 
  // Fill missing timebins
  //
#if 0
  unsigned short  boards[] = { 0x0, 0x1, 0x10, 0x11, 0xFFFF };
  unsigned short* board    = boards;
  while ((*boards) != 0xFFFF) { 
    for (UShort_t altro = 0; altro < 3; altro++) { 
      for (UShort_t channel = 0; channel < 16; channel++) { 
	Int_t idx = HWIndex(ddl, *board, altro, channel);
	if (idx < 0) { 
	  AliWarningF("Invalid index for %4d/0x%02x/0x%x/0x%x: %d", 
		      ddl, *board, altro, channel, idx);
	  continue;
	}
	Short_t min = fMinTimebin[idx];
	Short_t max = fMaxTimebin[idx];
	
	// Channel not seen at all. 
	if (min > 1023 || max < 0) continue;

	out << "# Extra timebins for 0x" << std::hex 
	    << board << ',' << altro << ',' << channel 
	    << " got time-bins " << min << " to " << max-1 
	    << std::dec << std::endl;
	
	for (UShort_t t = 15; t < min; t++) 
	  // Write a phony line 
	  out << board << "," << altro << "," << channel << "," 
	      << t << "," << 1023 << "," << 0 << std::endl;

	for (UShort_t t = max; t < 1024; t++) 
	  // Write a phony line 
	  out << board << "," << altro << "," << channel << "," 
	      << t << "," << 1023 << "," << 0 << std::endl;
      } // channel loop
    } // altro loop
  } // board loop
  // Write trailer, and close 
#endif
  out.write("# EOF\n", 6);
  out.close();
}

//_____________________________________________________________________
void AliFMDPedestalDA::WriteHeaderToFile() 
{
  //
  // Write headers to output files 
  //
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fOutputFile.write(Form("# %s \n",pars->GetPedestalShuttleID()),13);
  TDatime now;
  fOutputFile << "# This file created from run # " << fRunno 
	      << " @ " << now.AsString() << std::endl;
  fOutputFile.write("# Detector, "
		    "Ring, "
		    "Sector, "
		    "Strip, "
		    "Pedestal, "
		    "Noise, "
		    "Mu, "
		    "Sigma, "
		    "Chi2/NDF \n", 71);

  std::ostream* zss[] = { &fZSfileFMD1, &fZSfileFMD2, &fZSfileFMD3, 0 };
  for (size_t i = 0; i < 3; i++)  {
    *(zss[i]) << "# FMD " << (i+1) << " pedestals \n"
	      << "# board, "
	      << "altro, "
	      << "channel, "
	      << "timebin, "
	      << "pedestal, "
	      << "noise\n";
    *(zss[i]) << "# This file created from run # " << fRunno 
	      << " @ " << now.AsString() << std::endl;
  }
}

//_____________________________________________________________________
TH1S* AliFMDPedestalDA::GetChannel(UShort_t det, 
				   Char_t   ring, 
				   UShort_t sec, 
				   UShort_t strip,
				   UInt_t   sample) 
{
  // Get the histogram corresponding to a strip sample.
  // 
  // Parameters:
  //     det    Detector
  //     ring   Ring
  //     sec    Sector
  //     strip  Strip
  //     sample Sample
  //
  // Return:
  //     ADC spectra of a strip.
  Array* sampleArray = GetStripArray(det, ring, sec, strip);
  if (!sampleArray) return 0;
  TH1S*      hSample     = static_cast<TH1S*>(sampleArray->At(sample));
  if (!hSample) {
    AliErrorF("No channel histogram for FMD%d%c[%02d,%03d]_%d", 
	      det, ring, sec, strip, sample);
    sampleArray->ls();
    AliErrorF("Path is %s <- %s <- %s <- %s", 
	      sampleArray->GetName(),
	      GetSectorArray(det, ring, sec)->GetName(), 
	      GetRingArray(det, ring)->GetName(),
	      GetDetectorArray(det)->GetName());
    
  }
  return hSample;
  
}
//_____________________________________________________________________
TH1F* AliFMDPedestalDA::GetSectorSummary(UShort_t det, 
					 Char_t   ring, 
					 UShort_t sec, 
					 Bool_t   pedNotNoise) 
{
  Array*     secArray    = GetSectorArray(det, ring, sec);
  Int_t      n           = secArray->GetEntriesFast();
  Int_t      i           = n - (pedNotNoise ? kNoiseOffset : kPedestalOffset);
  return static_cast<TH1F*>(secArray->At(i));
}

//_____________________________________________________________________
//
//EOF
//
