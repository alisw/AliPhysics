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

/** @file    AliFMDPedestalDA.cxx
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
#include "iostream"
#include "fstream"
#include "AliLog.h"
#include "TF1.h"
#include "TObject.h"
#include "TMath.h"

//_____________________________________________________________________
ClassImp(AliFMDPedestalDA)

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA() : AliFMDBaseDA(),
  fCurrentChannel(1),
  fPedSummary("PedestalSummary","pedestals",51200,0,51200),
  fNoiseSummary("NoiseSummary","noise",51200,0,51200)
{
  fOutputFile.open("peds.csv");
  
}

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA(const AliFMDPedestalDA & pedDA) : 
  AliFMDBaseDA(pedDA),
  fCurrentChannel(1),
  fPedSummary("PedestalSummary","pedestals",51200,0,51200),
  fNoiseSummary("NoiseSummary","noise",51200,0,51200)
{
  
}

//_____________________________________________________________________
AliFMDPedestalDA::~AliFMDPedestalDA() 
{ 
}

//_____________________________________________________________________
void AliFMDPedestalDA::Init() 
{ 
  SetRequiredEvents(1000);
}

//_____________________________________________________________________
void AliFMDPedestalDA::AddChannelContainer(TObjArray* sectorArray, 
					   UShort_t det, 
					   Char_t   ring, 
					   UShort_t sec, 
					   UShort_t strip) 
{
#ifdef USE_SAMPLES
  AliFMDParameters* pars        = AliFMDParameters::Instance();
  Int_t             samples     = pars->GetSampleRate(det, ring, sec, strip);
  TObjArray*        sampleArray = new TObjArray(samples);
  for (size_t sample = 0; sample < samples; sample++) {
    TH1S* hSample = new TH1S(Form("FMD%d%c[%02d,03%d]_%d",
				  det,ring,sec,strip,sample),
			     Form("FMD%d%c[%02d,%03%d]_%d",
				  det,ring,sec,strip),
			     1024,-.5,1023.5);
    hSample->SetXTitle("ADC");
    hSample->SetYTitle("Events");
    sampleArray->AddAt(hSample, sample);
  }
  sectorArray->AddAtAndExpand(sampleArray, strip);
#else
  TH1S* hChannel = new TH1S(Form("hFMD%d%c_%d_%d",det,ring,sec,strip),
			    Form("hFMD%d%c_%d_%d",det,ring,sec,strip),
			    1024,-.5,1023.5);
  
  hChannel->SetDirectory(0);
  sectorArray->AddAtAndExpand(hChannel,strip);
#endif

}

//_____________________________________________________________________
void AliFMDPedestalDA::FillChannels(AliFMDDigit* digit) 
{
  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();

#ifdef USE_SAMPLES    
  AliFMDParameters* pars     = AliFMDParameters::Instance();
  Int_t             samples  = pars->GetSampleRate(det, ring, sec, strip);
  for (size_t sample = 0; sample < samples; sample++) {
    TH1S* hSample = GetChannel(det, ring, sec, strip, sample);
    hSample->Fill(digit->Count(sample));
  }
#else
  TH1S* hChannel = GetChannel(det, ring, sec, strip);
  hChannel->Fill(digit->Counts());
#endif
}

//_____________________________________________________________________
void AliFMDPedestalDA::Analyse(UShort_t det, 
			       Char_t   ring, 
			       UShort_t sec, 
			       UShort_t strip) 
{

  TH1S* hChannel       = GetChannel(det, ring, sec, strip);
  if(hChannel->GetEntries() == 0) {
    //  AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",
    //                  det,ring,sec,strip));
    return;
  }
 
  AliDebug(50, Form("Fitting FMD%d%c_%d_%d with %d entries",det,ring,sec,strip,
		   hChannel->GetEntries()));
  TF1 fitFunc("fitFunc","gausn",0,300);
  fitFunc.SetParameters(100,100,1);
  hChannel->Fit("fitFunc","Q0+","",10,200);
  
  Float_t mean = hChannel->GetMean();
  Float_t rms  = hChannel->GetRMS();
  
  
  
  hChannel->GetXaxis()->SetRangeUser(mean-5*rms,mean+5*rms);
  
  mean = hChannel->GetMean();
  rms  = hChannel->GetRMS();
  


  Float_t chi2ndf = 0;
  if(fitFunc.GetNDF())
    chi2ndf = fitFunc.GetChisquare() / fitFunc.GetNDF();
 
  fOutputFile << det                         << ','
	      << ring                        << ','
	      << sec                         << ','
              << strip                       << ','
	      << mean                        << ','
	      << rms                         << ','
	      << fitFunc.GetParameter(1)     << ','
	      << fitFunc.GetParameter(2)     << ','
	      << chi2ndf                     <<"\n";
  
  if(fSaveHistograms) {
    gDirectory->cd(GetSectorPath(det, ring, sec, kTRUE));
    TH1F* sumPed   = dynamic_cast<TH1F*>(gDirectory->Get("Pedestals"));
    TH1F* sumNoise = dynamic_cast<TH1F*>(gDirectory->Get("Noise"));
    Int_t nStr = (ring == 'I' ? 512 : 256);
    if (!sumPed) {
      sumPed = new TH1F("Pedestals", 
			Form("Summary of pedestals in FMD%d%c[%02d]", 
			     det, ring, sec), 
			nStr, -.5, nStr-.5);
      sumPed->SetXTitle("Strip");
      sumPed->SetYTitle("Pedestal [ADC]");
      sumPed->SetDirectory(gDirectory);
    }
    if (!sumNoise) { 
      sumNoise = new TH1F("Noise", 
			  Form("Summary of noise in FMD%d%c[%02d]", 
			       det, ring, sec), 
			  nStr, -.5, nStr-.5);
      sumNoise->SetXTitle("Strip");
      sumNoise->SetYTitle("Noise [ADC]");
      
      sumNoise->SetDirectory(gDirectory);
    }
    sumPed->SetBinContent(strip+1, mean);
    sumPed->SetBinError(strip+1, rms);
    sumNoise->SetBinContent(strip+1, rms);
    
    if(sumNoise->GetEntries() == nStr)
      sumNoise->Write(sumNoise->GetName(),TObject::kOverwrite);
    if(sumPed->GetEntries() == nStr)
      sumPed->Write(sumPed->GetName(),TObject::kOverwrite);
    
    fPedSummary.SetBinContent(fCurrentChannel,mean);
    fNoiseSummary.SetBinContent(fCurrentChannel,rms);
    fCurrentChannel++;
    
    gDirectory->cd(GetStripPath(det, ring, sec, strip, kTRUE));
    hChannel->GetXaxis()->SetRange(1,1024);
    
    hChannel->Write();
  }  
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::Terminate(TFile* diagFile) 
{
  diagFile->cd();
  
  fPedSummary.Write();
  fNoiseSummary.Write();
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::WriteHeaderToFile() 
{
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fOutputFile.write(Form("# %s \n",pars->GetPedestalShuttleID()),13);
  fOutputFile.write("# Detector, "
		    "Ring, "
		    "Sector, "
		    "Strip, "
		    "Pedestal, "
		    "Noise, "
		    "Mu, "
		    "Sigma, "
		    "Chi2/NDF \n", 71);
}

//_____________________________________________________________________
TH1S* AliFMDPedestalDA::GetChannel(UShort_t det, 
				   Char_t   ring, 
				   UShort_t sec, 
				   UShort_t strip) 
{
  UShort_t   iring     = (ring == 'O' ? 0 : 1);
  TObjArray* detArray  = static_cast<TObjArray*>(fDetectorArray.At(det));
  TObjArray* ringArray = static_cast<TObjArray*>(detArray->At(iring));
  TObjArray* secArray  = static_cast<TObjArray*>(ringArray->At(sec));
#ifdef USE_SAMPLES
  AliFMDParameters* pars        = AliFMDParameters::Instance();
  Int_t             samples     = pars->GetSampleRate(det, ring, sec, strip);
  TObjArray*        sampleArray = static_cast<TObjArray*>(secArray->At(strip));
  TH1S*      hSample = static_cast<TH1S*>(sampleArray->At(sample));
  return hSample;
#else
  TH1S*      hChannel  = static_cast<TH1S*>(secArray->At(strip));
  return hChannel;
#endif  
}

//_____________________________________________________________________
//
//EOF
//
