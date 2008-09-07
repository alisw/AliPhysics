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
  fNoiseSummary("NoiseSummary","noise",51200,0,51200),
  fZSfileFMD1(),
  fZSfileFMD2(),
  fZSfileFMD3()
{
  fOutputFile.open("peds.csv");
  fZSfileFMD1.open("ddl3072.csv");
  fZSfileFMD2.open("ddl3073.csv");
  fZSfileFMD3.open("ddl3074.csv");  
}

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA(const AliFMDPedestalDA & pedDA) : 
  AliFMDBaseDA(pedDA),
  fCurrentChannel(1),
  fPedSummary("PedestalSummary","pedestals",51200,0,51200),
  fNoiseSummary("NoiseSummary","noise",51200,0,51200),
  fZSfileFMD1(),
  fZSfileFMD2(),
  fZSfileFMD3()
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
  AliFMDParameters* pars        = AliFMDParameters::Instance();
  UInt_t             samples     = pars->GetSampleRate(det, ring, sec, strip);
  TObjArray*        sampleArray = new TObjArray(samples);
  sampleArray->SetOwner();
  for (UInt_t sample = 0; sample < samples; sample++) {
    TH1S* hSample = new TH1S(Form("FMD%d%c[%02d,03%d]_%d",
				  det,ring,sec,strip,sample),
			     Form("FMD%d%c[%02d,%03%d]_%d",
				  det,ring,sec,strip),
			     1024,-.5,1023.5);
    hSample->SetXTitle("ADC");
    hSample->SetYTitle("Events");
    hSample->SetDirectory(0);
    sampleArray->AddAt(hSample, sample);
  }
  sectorArray->AddAtAndExpand(sampleArray, strip);
  

}

//_____________________________________________________________________
void AliFMDPedestalDA::FillChannels(AliFMDDigit* digit) 
{
  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();

  AliFMDParameters* pars     = AliFMDParameters::Instance();
  UInt_t             samples  = pars->GetSampleRate(det, ring, sec, strip);
  for (UInt_t sample = 0; sample < samples; sample++) {
    TH1S* hSample = GetChannel(det, ring, sec, strip, sample);
    hSample->Fill(digit->Count(sample));
  }
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::Analyse(UShort_t det, 
			       Char_t   ring, 
			       UShort_t sec, 
			       UShort_t strip) {
  
  AliFMDParameters* pars     = AliFMDParameters::Instance();
  UInt_t             samples  = pars->GetSampleRate(det, ring, sec, strip);
  for (UShort_t sample = 0; sample < samples; sample++) {
  
    TH1S* hChannel       = GetChannel(det, ring, sec, strip,sample);
    if(hChannel->GetEntries() == 0) {
      //AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",
      //	    det,ring,sec,strip));
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
  
    
    UShort_t ddl, board, altro, channel;
    UShort_t timebin;
    
    pars->Detector2Hardware(det,ring,sec,strip,sample,ddl,board,altro,channel,timebin);
    
    switch(det) {
    case 1:
      fZSfileFMD1 << board << ',' << altro << ',' << channel << ',' << timebin << ','
		  << mean  << ',' << rms << "\n"; break;
    case 2:
      fZSfileFMD2 << board << ',' << altro << ',' << channel << ',' << timebin << ','
		  << mean  << ',' << rms << "\n"; break;
    case 3:
      fZSfileFMD3 << board << ',' << altro << ',' << channel << ',' << timebin << ','
		  << mean  << ',' << rms << "\n"; break;
    default:
      AliWarning("Unknown sample!"); break;
      
    }
    
    Float_t chi2ndf = 0;
    if(fitFunc.GetNDF())
      chi2ndf = fitFunc.GetChisquare() / fitFunc.GetNDF();
    
    if(sample==samples-1) {
    
      fOutputFile << det                         << ','
		  << ring                        << ','
		  << sec                         << ','
		  << strip                       << ','
		  << mean                        << ','
		  << rms                         << ','
		  << fitFunc.GetParameter(1)     << ','
		  << fitFunc.GetParameter(2)     << ','
		  << chi2ndf                     <<"\n";
      
      if(fSaveHistograms  ) {
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
  }
}

//_____________________________________________________________________
void AliFMDPedestalDA::Terminate(TFile* diagFile) 
{
  if(fSaveHistograms) {
    diagFile->cd();
    
    fPedSummary.Write();
    fNoiseSummary.Write();
  }
  
  if(fZSfileFMD1.is_open()) { 
    fZSfileFMD1.write("# EOF\n",6);
    fZSfileFMD1.close();  }
  if(fZSfileFMD2.is_open()) {
    fZSfileFMD2.write("# EOF\n",6);
    fZSfileFMD2.close(); }
  if(fZSfileFMD3.is_open()) {
    fZSfileFMD3.write("# EOF\n",6);
    fZSfileFMD3.close(); }
  
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
  fZSfileFMD1.write("# FMD 1 pedestals \n",19);
  fZSfileFMD1.write("# board, "
		    "altro, "
		    "channel, "
		    "timebin, "
		    "pedestal, "
		    "noise \n", 51);
  fZSfileFMD2.write("# FMD 2 pedestals \n",19);
  fZSfileFMD2.write("# board, "
		    "altro, "
		    "channel, "
		    "timebin, "
		    "pedestal, "
		    "noise \n", 51);
  fZSfileFMD3.write("# FMD 3 pedestals \n",19);
  fZSfileFMD3.write("# board, "
		    "altro, "
		    "channel, "
		    "timebin, "
		    "pedestal, "
		    "noise \n", 51);
  
}

//_____________________________________________________________________
TH1S* AliFMDPedestalDA::GetChannel(UShort_t det, 
				   Char_t   ring, 
				   UShort_t sec, 
				   UShort_t strip,
				   UInt_t    sample) 
{
  UShort_t   iring     = (ring == 'O' ? 0 : 1);
  TObjArray* detArray  = static_cast<TObjArray*>(fDetectorArray.At(det));
  TObjArray* ringArray = static_cast<TObjArray*>(detArray->At(iring));
  TObjArray* secArray  = static_cast<TObjArray*>(ringArray->At(sec));
  TObjArray*        sampleArray = static_cast<TObjArray*>(secArray->At(strip));
  TH1S*      hSample = static_cast<TH1S*>(sampleArray->At(sample));
  return hSample;
  
}

//_____________________________________________________________________
//
//EOF
//
