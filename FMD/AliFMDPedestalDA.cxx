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
// This class implements the virtual functions of the AliFMDBaseDA class.
// The most important of these functions, FillChannels(..) and Analyse(...) collect 
// and analyse the data of each channel. The resulting pedestal and noise values are
// written to a comma separated values (csv) file on the go. The csv files produced 
// in this way are the basic input to the AliFMDPreprocessor.
//

#include "AliFMDPedestalDA.h"
#include "iostream"
#include "fstream"
#include "AliLog.h"
#include "TF1.h"

//_____________________________________________________________________
ClassImp(AliFMDPedestalDA)

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA() : AliFMDBaseDA()
{
  fOutputFile.open("peds.csv");
  
}

//_____________________________________________________________________
AliFMDPedestalDA::AliFMDPedestalDA(const AliFMDPedestalDA & pedDA) : 
  AliFMDBaseDA(pedDA)
{
  
}

//_____________________________________________________________________
AliFMDPedestalDA::~AliFMDPedestalDA() {
 
}

//_____________________________________________________________________
void AliFMDPedestalDA::Init() {
 
  SetRequiredEvents(1000);
}

//_____________________________________________________________________
void AliFMDPedestalDA::AddChannelContainer(TObjArray* sectorArray, 
					   UShort_t det, 
					   Char_t ring, 
					   UShort_t sec, 
					   UShort_t strip) {

  TH1S* hChannel = new TH1S(Form("hFMD%d%c_%d_%d",det,ring,sec,strip),
			    Form("hFMD%d%c_%d_%d",det,ring,sec,strip),
			    1024,0,1023);
  
  hChannel->SetDirectory(0);
  sectorArray->AddAtAndExpand(hChannel,strip);
}

//_____________________________________________________________________
void AliFMDPedestalDA::FillChannels(AliFMDDigit* digit) {
  UShort_t det   = digit->Detector();
  Char_t   ring  = digit->Ring();
  UShort_t sec   = digit->Sector();
  UShort_t strip = digit->Strip();
    
  TH1S* hChannel       = GetChannel(det, ring, sec, strip);
  hChannel->Fill(digit->Counts());
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::Analyse(UShort_t det, 
			       Char_t ring, 
			       UShort_t sec, 
			       UShort_t strip) {

  TH1S* hChannel       = GetChannel(det, ring, sec, strip);
  if(hChannel->GetEntries() == 0) {
    //  AliWarning(Form("No entries for FMD%d%c, sector %d, strip %d",det,ring,sec,strip));
    return;
  }
  
  //std::cout<<"Fitting Channel #"<<Form("FMD%d%c_%d_%d",det,ring,sec,strip)<<" with "<<hChannel->GetEntries()<<" entries     \r"<<std::flush;
  
  Float_t mean = hChannel->GetMean();
  Float_t rms  = hChannel->GetRMS();
  
  hChannel->GetXaxis()->SetRangeUser(mean-4*rms,mean+4*rms);
  
  mean = hChannel->GetMean();
  rms  = hChannel->GetRMS();
  
  hChannel->Fit("gaus","QO","QO",mean-5*rms,mean+5*rms);
  TF1* fitFunc = hChannel->GetFunction("gaus");
    
  UInt_t ddl, board, chip, channel;
  
  UShort_t relStrip = strip%128;
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Detector2Hardware(det,ring,sec,strip,ddl,board,chip,channel);
  Float_t chi2ndf = 0;
  if(fitFunc->GetNDF())
    chi2ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();
  ddl = ddl + kBaseDDL;
  
  fOutputFile << ddl                         << ','
	      << board                       << ','
	      << chip                        << ','
	      << channel                     << ','
	      << relStrip                    << ','
	      << 0                           << ','
	      << 0                           << ','
	      << mean                        << ','
	      << rms                         << ','
	      << fitFunc->GetParameter(1)    << ','
	      << fitFunc->GetParameter(2)    << ','
	      << chi2ndf                     <<"\n";
  
  if(fSaveHistograms) {
    gDirectory->cd(Form("%s:FMD%d%c/sector_%d/strip_%d",fDiagnosticsFilename.Data(),det,ring,sec,strip));
    hChannel->GetXaxis()->SetRange(0,1023);
    hChannel->Write();
    
  }
  
}

//_____________________________________________________________________
void AliFMDPedestalDA::WriteHeaderToFile() {
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fOutputFile.write(Form("# %s \n",pars->GetPedestalShuttleID()),13);
  fOutputFile.write("# Rcu, Board, Chip, Channel, Strip, Sample, TimeBin, Pedestal, Noise, Mu, Sigma, Chi2/NDF \n",91);
  
}

//_____________________________________________________________________
TH1S* AliFMDPedestalDA::GetChannel(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip) {
  
  UShort_t  Ring = 1;
  if(ring == 'O')
    Ring = 0;
  
  
  TObjArray* detArray  = static_cast<TObjArray*>(fDetectorArray.At(det));
  TObjArray* ringArray = static_cast<TObjArray*>(detArray->At(Ring));
  TObjArray* secArray  = static_cast<TObjArray*>(ringArray->At(sec));
  TH1S* hChannel       = static_cast<TH1S*>(secArray->At(strip));
  
  return hChannel;
}
//_____________________________________________________________________
//
//EOF
//
