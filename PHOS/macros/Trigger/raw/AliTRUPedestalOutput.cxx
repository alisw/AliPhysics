/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Henrik Qvigstad <henrik.qvigstad@cern.ch>
/* $Id$ */


#include "AliTRUPedestalOutput.h"

ClassImp(AliTRUPedestalOutput)

#include "TMath.h"
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include "TCanvas.h"

#include <iostream>
using namespace std;

AliTRUPedestalOutput::AliTRUPedestalOutput()
: fRun(-1),
  fNEvents(0),
  fPedestals(0),
  fPedestalRMS(0),
  fPedestalSamples(0),
  fPedestals2d({0}),
  fPedestalRMS2d({0}),
  fPedestalsId(0),
  fPedestals_branch({{{0}}}),
  fTRUSignals({{{{{0}}}}})
{

}

AliTRUPedestalOutput::~AliTRUPedestalOutput()
{
  // All the class contains is Histograms, they are owned by current ROOT folder.
  // new TCanvas;
  // fTRUSignals[2][0][1][0][1]->DrawCopy();
}

void AliTRUPedestalOutput::SetRun ( Int_t run )
{
  if( fRun >= 0 ) // fRun is initilized to -1
    Error("SetRun", "Only one run allowed, i.e. run should not be set, then set to something different.");
  fRun = run;
}

void AliTRUPedestalOutput::EventAdded()
{
  fNEvents++;
}

TH1F* AliTRUPedestalOutput::GetPedestals()
{
  // Delete Histogram if exists
  if( fPedestals )
    fPedestals->Delete();
  fPedestals = 0;

  // Create Histogram
  fPedestals = new TH1F("fPedestals", "TRU Pedestals", 10*1024, 0, 1024);
  fPedestals->GetXaxis()->SetTitle("Pedestals, <signal>_c");
  fPedestals->GetYaxis()->SetTitle("Count");

  // Fill Histogram
  for (unsigned int mod=0; mod<kNMods; ++mod) {
    for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    double pedestal = GetPedestal(mod, row, branch, xrow, zcol);
	    fPedestals->Fill(pedestal);
	  }
	}
      }
    }
  }

  return fPedestals;
}

TH1F* AliTRUPedestalOutput::GetPedestalRMS()
{
  // Delete Histogram if exists
  if( fPedestalRMS )
    fPedestalRMS->Delete();
  fPedestalRMS = 0;

  fPedestalRMS = new TH1F("fPedestalRMS", "TRU Pedestal RMS", 1020, 0, 102);
  fPedestalRMS->GetXaxis()->SetTitle("RMS");
  fPedestalRMS->GetYaxis()->SetTitle("Count");

  for (unsigned int mod=0; mod<kNMods; ++mod) {
    for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    double rms = GetRMS(mod, row, branch, xrow, zcol);
	    fPedestalRMS->Fill(rms);
	    if ( 2 < rms )
	      printf("AliTRUPedestalOutput::GetPedestalRMS: ped. RMS:%f, mod:%d row:%d branch:%d x:%02d z:%02d\n", rms, mod, row, branch, xrow, zcol);

	  }
	}
      }
    }
  }

  return fPedestalRMS;
}

TH1I* AliTRUPedestalOutput::GetPedestalSamples()
{
  if( fPedestalSamples )
    fPedestalSamples->Delete();
  fPedestalSamples = 0;

  fPedestalSamples = new TH1I("fPedestalSamples", "TRU Pedestal Samples", 1000, 0, kNTRUTimeBins*fNEvents);
  fPedestalSamples->GetXaxis()->SetTitle("Samples");
  fPedestalSamples->GetYaxis()->SetTitle("Count");

  int smallThreshold = 2600;
  for (unsigned int mod=0; mod<kNMods; ++mod) {
    for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	int smallCount = 0;
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    double entries = GetSamples(mod, row, branch, xrow, zcol);
	    fPedestalSamples->Fill(entries);
	    if ( entries < smallThreshold )
	      ++smallCount;
	  } // z
	} // x 
	if(smallCount)
	  printf("AliTRUPedestalOutput::GetPedestalSamples: Small Sample! %d cells in TRU: mod:%d row:%d branch:%d\n", smallCount, mod, row, branch);
      } // branch 
    } // row
  } // mod

  return fPedestalSamples;
}

TH2F* AliTRUPedestalOutput::GetPedestals2d(UInt_t mod)
{
  if( fPedestals2d[mod] )
    fPedestals2d[mod]->Delete();
  fPedestals2d[mod] = 0;

  const int nhr = kN2x2X +kNTRURows -1; // number of histogram rows, x
  const int nhc = kN2x2Z +kNBranches -1; // number of histogram columns, z
  TString name("fPedestals2d_mod:"); name += mod;
  TString title("fPedestals, Mod: "); title += mod;
  fPedestals2d[mod] = new TH2F(name, title, nhr, 0, nhr,   nhc, 0, nhc);
  fPedestals2d[mod]->GetXaxis()->SetTitle("x = x_tru + 9*tru");
  fPedestals2d[mod]->GetYaxis()->SetTitle("z = z_tru + 15*branch");
  fPedestals2d[mod]->SetStats(0);

  for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    int xindex = xrow + kN2x2XPrTRURow*row +row;
	    int zindex = zcol + kN2x2ZPrBranch*branch +branch;
	    double pedestal = GetPedestal(mod, row, branch, xrow, zcol);

    	    fPedestals2d[mod]->SetBinContent(xindex+1, zindex+1, pedestal);
	  } // z col
	} // x row
      } // branch, tru
    } // row, tru
  fPedestals2d[mod]->SetDrawOption("colz");
  return fPedestals2d[mod];
}

TH2F* AliTRUPedestalOutput::GetPedestalRMS2d(UInt_t mod)
{
  if( fPedestalRMS2d[mod] )
    fPedestalRMS2d[mod]->Delete();
  fPedestalRMS2d[mod] = 0;

  const int nhr = kN2x2X +kNTRURows -1; // number of histogram rows, x
  const int nhc = kN2x2Z +kNBranches -1; // number of histogram columns, z
  TString name("fPedestalRMS2d_mod:"); name += mod;
  TString title("RMS, Mod: "); title += mod;
  fPedestalRMS2d[mod] = new TH2F(name, title, nhr, 0, nhr,   nhc, 0, nhc);
  fPedestalRMS2d[mod]->GetXaxis()->SetTitle("x = x_tru + 9*tru");
  fPedestalRMS2d[mod]->GetYaxis()->SetTitle("z = z_tru + 15*branch");
  fPedestalRMS2d[mod]->SetStats(0);

  for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    int xindex = xrow + kN2x2XPrTRURow*row +row;
	    int zindex = zcol + kN2x2ZPrBranch*branch +branch;
	    double rms = GetRMS(mod, row, branch, xrow, zcol);

    	    fPedestalRMS2d[mod]->SetBinContent(xindex+1, zindex+1, rms);
	  } // z col
	} // x row
      } // branch, tru
    } // row, tru
  fPedestalRMS2d[mod]->SetDrawOption("colz");
  return fPedestalRMS2d[mod];
}

TH1F* AliTRUPedestalOutput::GetPedestalsId()
{
  if( fPedestalsId )
    fPedestalsId->Delete();
  fPedestalsId = 0;

  // else Initilize
  int nhc = kNMods*kNTRURows*kNBranches*kN2x2XPrTRURow*kN2x2ZPrBranch;
  TString title("Id v Pedestal");
  fPedestalsId = new TH1F("fPedestalsId", title.Data(), nhc, 0, nhc);
  fPedestalsId->GetXaxis()->SetTitle("Id ([mod][row][branch][x][z], Row-major order)");
  fPedestalsId->GetYaxis()->SetTitle("Pedestal");
  fPedestalsId->SetStats(0);

  for (unsigned int mod=0; mod<kNMods; ++mod) {
    for (unsigned int row=0; row<kNTRURows; ++row) {
      for (unsigned int branch=0; branch<kNBranches; ++branch) {
	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	      continue;

	    int index =
	      mod*kNTRURows*kNBranches*kN2x2XPrTRURow*kN2x2ZPrBranch
	    + row*kNBranches*kN2x2XPrTRURow*kN2x2ZPrBranch
	    + branch*kN2x2XPrTRURow*kN2x2ZPrBranch
	    + xrow*kN2x2ZPrBranch
	    + zcol;
	    double pedestal = GetPedestal(mod, row, branch, xrow, zcol);
	    double pederror = GetPedestalError(mod, row, branch, xrow, zcol);
	    fPedestalsId->SetBinContent(index+1, pedestal);
	    fPedestalsId->SetBinError(index+1, pederror);
	  }
	}
      }
    }
  }
  fPedestalsId->SetDrawOption("E");
  return fPedestalsId;
}

TH1F* AliTRUPedestalOutput::GetPedestals_branch ( UInt_t mod, UInt_t row, UInt_t branch )
{
  if( fPedestals_branch[mod][row][branch] )
    fPedestals_branch[mod][row][branch]->Delete();
  fPedestals_branch[mod][row][branch] = 0;
  
  TH1F* hist = new TH1F("fPedestals", "Pedestals", 10*1024, 0, 1024);
  hist->GetXaxis()->SetTitle("Pedestals, <signal>");
  hist->GetYaxis()->SetTitle("Count");

  for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
    for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
      if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
	continue;
      double pedestal = GetPedestal(mod, row, branch, xrow, zcol);
      hist->Fill(pedestal);
    }
  }

  return fPedestals_branch[mod][row][branch] = hist;
}

TH1I* AliTRUPedestalOutput::GetTRUSignals ( UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z )
{
  if( fTRUSignals[mod][row][branch][x][z] )
    return fTRUSignals[mod][row][branch][x][z];
  else
  {
    // else Initilize
    char name[256];
    sprintf(name, "fTRUSignals_m:%d_r:%d_b:%d_x:%02d_z:%02d", mod, row, branch, x, z);
    char title[256];
    sprintf(title, "TRU Signal mod:%d row:%d branch:%d x:%02d z:%02d", mod, row, branch, x, z);
  
    TH1I* hist = new TH1I(name, title, 1024, 0, 1024);
    hist->GetXaxis()->SetTitle("Bin");
    hist->GetYaxis()->SetTitle("Count");
    return fTRUSignals[mod][row][branch][x][z] = hist;
  }
}

Double_t AliTRUPedestalOutput::GetPedestal ( UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z )
{
  
  if( fTRUSignals[mod][row][branch][x][z] )
    return fTRUSignals[mod][row][branch][x][z]->GetMean();
  else
    return -1;
}

Double_t AliTRUPedestalOutput::GetPedestalError ( UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z )
{
  if( fTRUSignals[mod][row][branch][x][z] )
    return fTRUSignals[mod][row][branch][x][z]->GetMeanError();
  else
    return -1;

}

Double_t AliTRUPedestalOutput::GetRMS ( UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z )
{
  if( fTRUSignals[mod][row][branch][x][z] )
    return fTRUSignals[mod][row][branch][x][z]->GetRMS();
  else
    return -1;
}

Double_t AliTRUPedestalOutput::GetSamples(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z)
{
  if( fTRUSignals[mod][row][branch][x][z] )
    return fTRUSignals[mod][row][branch][x][z]->GetEntries();
  else
    return -1;
}

// TH2F* AliTRUPedestalOutput::GetPedestalRMS2d_old(UInt_t mod)
// {
//   if( fPedestalRMS2d[mod] )
//     fPedestalRMS2d[mod]->Delete();
//   fPedestalRMS2d[mod] = 0;
//
//   const int nhr = 2*kN2x2X +2*kNTRURows; // number of histogram rows, x
//   const int nhc = 2*kN2x2Z +2*kNBranches; // number of histogram columns, z
//   TString name("fPedestalRMS2d_mod:"); name += mod;
//   TString title("RMS, Mod: "); title += mod;
//   fPedestalRMS2d[mod] = new TH2F(name, title, nhr, 0, nhr,   nhc, 0, nhc);
//   fPedestalRMS2d[mod]->GetXaxis()->SetTitle("x");
//   fPedestalRMS2d[mod]->GetYaxis()->SetTitle("z");
//
//   for (unsigned int row=0; row<kNTRURows; ++row) {
//       for (unsigned int branch=0; branch<kNBranches; ++branch) {
// 	for (unsigned int xrow=0; xrow<kN2x2XPrTRURow; ++xrow) {
// 	  for (unsigned int zcol=0; zcol<kN2x2ZPrBranch; ++zcol) {
// 	    if( ! fTRUSignals[mod][row][branch][xrow][zcol] )
// 	      continue;
//
// 	    int xindex = 2*xrow + 2*kN2x2XPrTRURow*row +2*row;
// 	    int zindex = 2*zcol + 2*kN2x2ZPrBranch*branch +2*branch;
// 	    double rms = GetRMS(mod, row, branch, xrow, zcol);
//
//     	    fPedestalRMS2d[mod]->SetBinContent(xindex+1, zindex+1, rms);
// 	  } // z col
// 	} // x row
//       } // branch, tru
//     } // row, tru
//   fPedestalRMS2d[mod]->SetDrawOption("colz");
//   return fPedestalRMS2d[mod];
// }
