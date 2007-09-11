/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
*/ 

//_________________________________________________________________________
// Cell folder which will keep all information about cell(tower) itself
//  Initial version was created with TDataSet staf
//  TObjectSet -> TFolder; Sep 6, 2007
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALCell.h"
#include "AliEMCALHistoUtilities.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALFolder.h"
#include "AliEMCALSuperModule.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALRecPointsQaESDSelector.h"

#include "AliEMCALCalibCoefs.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include <TH1.h>
#include <TF1.h>
#include <TNtuple.h>

typedef  AliEMCALHistoUtilities u;

ClassImp(AliEMCALCell)

Double_t ADCCHANNELEC = 0.0153;  // Update 24 Apr 2007: 250./16/1024 - width of one ADC channel in GeV
Double_t MPI0         = 0.13498; // mass of pi0
Double_t MPI02        = MPI0*MPI0; // mass**2

AliEMCALCell::AliEMCALCell() : 
TFolder(), 
 fParent(0),fLh(0),
fAbsId(0),fSupMod(0),fModule(0),fPhi(0),fEta(0),fPhiCell(0),fEtaCell(0),fCcIn(0),fCcOut(0),
fFun(0)
{
}

AliEMCALCell::AliEMCALCell(const Int_t absId, const char* title) : 
  TFolder(Form("Cell%4.4i",absId),title), 
 fParent(0),fLh(0),
fAbsId(absId),fSupMod(0),fModule(0),fPhi(0),fEta(0),fPhiCell(0),fEtaCell(0),fCcIn(0),fCcOut(0),
fFun(0)
{
  
  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();
  g->GetCellIndex(fAbsId, fSupMod, fModule, fPhi, fEta);
  g->GetCellPhiEtaIndexInSModule(fSupMod, fModule, fPhi, fEta, fPhiCell, fEtaCell);

} 

AliEMCALCell::~AliEMCALCell()
{
  // dtor
}
//-------------------------------------------------------------------------------------

void AliEMCALCell::SetCCfromDB(AliEMCALCalibData *ccDb)
{
  if(ccDb == 0) return;
  // fADCchannelEC = fCalibData->GetADCchannel(iSupMod,ieta,iphi);
  // fetaCel-column; fPhiCell- row
  fCcIn = ccDb->GetADCchannel(fSupMod, fEtaCell, fPhiCell); // in GeV

  TH1* h = (TH1*)GetHists()->At(0);
  u::AddToNameAndTitle(h, 0, Form(", cc %5.2f MeV", fCcIn*1.e+3));
}

void AliEMCALCell::SetCCfromCCTable(AliEMCALCalibCoefs *t)
{
  if(t == 0) {
    //    Dump();
    return;
  }
  if(fLh == 0) {
    fLh = BookHists();
    Add(fLh);
  }

  calibCoef *r = t->GetTable(fAbsId);
  if(r && r->absId == fAbsId) {
    fCcIn = r->cc;
  } else { // something wrong
    if(r) printf(" fAbsId %i : r->absId %i \n", fAbsId, r->absId);
    assert(0);
  }

  TH1* h = (TH1*)GetHists()->At(0);
  u::AddToNameAndTitle(h, 0, Form(", cc %5.2f MeV", fCcIn*1.e+3));
}

void AliEMCALCell::FillEffMass(const Double_t mgg)
{
  u::FillH1(GetHists(), 0, mgg);
}

void AliEMCALCell::FillCellNtuple(TNtuple *nt)
{
  if(nt==0) return;
  nt->Fill(fAbsId,fSupMod,fModule,fPhi,fEta,fPhiCell,fEtaCell,fCcIn,fCcOut);
}

void AliEMCALCell::FitHist(TH1* h, const char* name, const char* opt)
{
  TString optFit(""), OPT(opt);
  OPT.ToUpper();
  if(h==0) return; 
  printf("<I> AliEMCALCell::FitHist : |%s| is started : opt %s\n", h->GetName(), opt);
  TString tit(h->GetTitle());

  TF1 *GausPol2 = 0, *g=0, *bg=0;
  if(h->GetListOfFunctions()->GetSize() == 0 || 1) {
    g = u::Gausi(name, 0.0, 0.4, h); // gaus estimation

    g->SetParLimits(0, h->Integral()/20., h->Integral());
    g->SetParLimits(1, 0.08, 0.20);
    g->SetParLimits(2, 0.001, 0.02);

    g->SetParameter(0, 1200.);
    g->SetParameter(1, h->GetBinLowEdge(h->GetMaximumBin()));
    g->SetParameter(2, 0.010);

    g->SetLineColor(kRed);
    g->SetLineWidth(2);

    optFit = "0NQ";
    h->Fit(g, optFit.Data(),"", 0.001, 0.3);

    optFit = "0NQIME";
    h->Fit(g, optFit.Data(),"", 0.001, 0.3);

    bg = new TF1(Form("bg%s",name), "pol2", 0.0, 0.3);  
    optFit = "0NQ";
    h->Fit(bg, optFit.Data(),"", 0.0, 0.3);

    GausPol2 = u::GausiPol2(name, 0.00, 0.3, g, bg);
    optFit = "0Q";
    h->Fit(GausPol2, optFit.Data(),"", 0.03, 0.28);
  // Clean up
    delete g;
    delete bg;
    optFit = "0IME+"; // no drwaing at all
    if(tit.Contains("SM") || OPT.Contains("DRAW")) optFit = "IME+";
  } else {
    GausPol2 = (TF1*)h->GetListOfFunctions()->At(0);
    optFit = "IME+";
    printf("<I> Function is defined alredy : %s optFit %s \n", GausPol2->GetTitle(), optFit.Data());
  }
  //  optFit = "IME+";
  h->Fit(GausPol2, optFit.Data(),"", 0.01, 0.28);

  if(optFit.Contains("0") == 0) {
    gStyle->SetOptFit(111);
    u::DrawHist(h,2);
  }
  printf("<I> AliEMCALCell::FitHist : |%s| is ended \n\n", h->GetName());
}

void AliEMCALCell::FitEffMassHist(const char* opt)
{
  AliEMCALFolder* EMCAL = AliEMCALRecPointsQaESDSelector::GetEmcalFolder();
  Int_t it = EMCAL->GetIterationNumber();

  TH1* h = (TH1*)GetHists()->At(0);

  FitHist(h, GetName(), opt);

  fFun = (TF1*)h->GetListOfFunctions()->At(0);
  if(fFun) {
    Double_t mpi = fFun->GetParameter(1), mpi2 = mpi*mpi;
    Double_t ccTmp = fCcIn * MPI02 / mpi2;
    if(it<=1) { // Jul 16, 2007
      fCcOut = ccTmp;
    } else {
      fCcOut = (ccTmp + fCcIn)/2.;
    }
    printf(" fFun %s | %s : iet %i\n", fFun->GetName(), fFun->GetTitle(), it);
  }
  printf(" %s | fCcIn %6.5f -> % 6.5f <- fCcOut \n", GetTitle(), fCcIn , fCcOut);
}

void AliEMCALCell::PrintInfo()
{
  printf(" %s %s \n", GetName(), GetTitle());
  if(fLh == 0 ) return;
  TH1* h = (TH1*)GetHists()->At(0);
  TF1 *f = (TF1*)h->GetListOfFunctions()->At(0);
  if(fFun) printf(" fFun : %s | %s \n", fFun->GetName(),  fFun->GetTitle());
  else fFun = f;
  // if(f) f->Dump();
}

TList* AliEMCALCell::BookHists()
{
  gROOT->cd();
  TH1::AddDirectory(1);

  AliEMCALFolder* EMCAL = AliEMCALRecPointsQaESDSelector::GetEmcalFolder();
  Int_t it = EMCAL->GetIterationNumber();

  new TH1F("01_EffMass", "effective mass of #gamma,#gamma(m_{#pi^{0}}=134.98 MeV) ", 60,0.0,0.3);

  TList *l = AliEMCALHistoUtilities::MoveHistsToList(Form("HistsOfCell%4.4i",fAbsId), kFALSE);
  AliEMCALHistoUtilities::AddToNameAndTitleToList(l, Form("4.4i_It%i",fAbsId,it),
						  Form(" Cell %4.4i, iter. %i",fAbsId, it));

  TH1::AddDirectory(0);
  return l;
}
