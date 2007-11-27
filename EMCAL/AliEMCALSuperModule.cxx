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
Revision 1.3  2007/11/14 15:34:05  gustavo
Coding violations corrected

Revision 1.2  2007/09/11 19:38:15  pavlinov
added pi0 calibration, linearity, shower profile
 
*/

//_________________________________________________________________________
// Top EMCAL folder which will keep all information about EMCAL itself,
// super Modules (SM), modules, towers, set of hists and so on.
//  TObjectSet -> TFolder; Sep 6, 2007
//
//
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALSuperModule.h"
#include "AliEMCALFolder.h"
#include "AliEMCALCell.h"
#include "AliEMCALHistoUtilities.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

typedef  AliEMCALHistoUtilities u;

ClassImp(AliEMCALSuperModule)

AliEMCALSuperModule::AliEMCALSuperModule() : TFolder()
, fParent(0),fLh(0),fSMNumber(0)
{
  //default ctor
}

AliEMCALSuperModule::AliEMCALSuperModule(const Int_t m, const char* title) : 
TFolder(Form("SM%2.2i",m), title)
, fParent(0),fLh(0),fSMNumber(m)
{ 
  //ctor
} 

AliEMCALSuperModule::AliEMCALSuperModule(const AliEMCALSuperModule &sm) : 
TFolder(), fParent(sm.fParent),fLh(sm.fLh),fSMNumber(sm.fSMNumber)
{ 
  //copy ctor
} 

AliEMCALSuperModule & AliEMCALSuperModule::operator =(const AliEMCALSuperModule &sm)
{ 
    // assignment operator
//   if(this == &sm)return *this;
//   ((TObject *)this)->operator=(sm);

  fParent = sm.fParent; 
  fLh = sm.fLh;
  fSMNumber = sm.fSMNumber;
  return *this;
} 

AliEMCALSuperModule::~AliEMCALSuperModule()
{
  // dtor
}

void AliEMCALSuperModule::Init()
{
  //Initialization method

  if(GetHists()==0) {
    fLh = BookHists();
    Add(fLh);
  }
}

void  AliEMCALSuperModule::AddCellToEtaRow(AliEMCALCell *cell, const Int_t etaRow)
{
  //Adds cells to corresponding Super module Eta Row
  static TFolder *set;
  set = dynamic_cast<TFolder*>(FindObject(Form("ETA%2.2i",etaRow))); 
  if(set==0) {
    set = new  TFolder(Form("ETA%2.2i",etaRow),"eta row");
    Add(set);
  }
  set->Add(cell);
}

void AliEMCALSuperModule::FitForAllCells()
{
  //Fit histograms of each cell

  Int_t ncells=0;
  for(int eta=0; eta<48; eta++) { // eta row
    TFolder *setEta = dynamic_cast<TFolder*>(FindObject(Form("ETA%2.2i",eta)));
    if(setEta) {
      TList* l = (TList*)setEta->GetListOfFolders();
      printf(" eta %i : %s : cells %i \n", eta, setEta->GetName(), l->GetSize());
      for(int phi=0; phi<l->GetSize(); phi++) { // cycle on cells (phi directions)
        AliEMCALCell* cell = dynamic_cast<AliEMCALCell*>(l->At(phi));
        if(cell == 0) continue; 

        cell->FitEffMassHist();

        u::FillH1(GetHists(), 1, cell->GetCcIn()*1.e+3);
        u::FillH1(GetHists(), 2, cell->GetCcOut()*1.e+3);

        TF1 *f = cell->GetFunction();
        u::FillH1(GetHists(), 3, f->GetParameter(1));
        u::FillH1(GetHists(), 4, f->GetParameter(2));
        u::FillH1(GetHists(), 5, f->GetChisquare()/f->GetNDF());
        u::FillH1(GetHists(), 6, f->GetParameter(0));

        ncells++;
      }
    }
  }
  printf(" <I> AliEMCALSuperModule::FitForAllCells() : ncells %i with fit \n", ncells);
}

void AliEMCALSuperModule::FitEffMassHist()
{
  //Fit effective mass histogram
  TH1* h = (TH1*)GetHists()->At(0);
  AliEMCALCell::FitHist(h, GetName());
}


void AliEMCALSuperModule::PrintInfo()
{
  //Print
  printf(" Super Module :   %s    :   %i \n", GetName(), fSMNumber);
  printf(" # of active cells                %i \n", GetNumberOfCells());
  TH1* h = (TH1*)GetHists()->At(0);
  printf(" # h->Integral() of eff.mass hist %i \n", int(h->Integral()));
}

void AliEMCALSuperModule::DrawCC(int iopt)
{
  //Draw different cell histograms
  TCanvas *c=0; 
  if(iopt==1) c = new TCanvas("ccInOut","ccInOut");
 
  gStyle->SetOptStat(0);

  TH1 *hCCIn  = (TH1*)GetHists()->At(1);
  TH1 *hCCOut = (TH1*)GetHists()->At(2);

  if(hCCIn == 0)              return;
  if(hCCIn->GetEntries()<10.) return;

  hCCIn->SetStats(kFALSE);
  hCCOut->SetStats(kFALSE);
  hCCOut->SetTitle("CC in and out; Iter 1; Jul 26; All Statistics");
  hCCOut->SetXTitle("cc in MeV");
  hCCOut->SetYTitle("  N  ");

  u::DrawHist(hCCOut,2);
  hCCOut->SetAxisRange(10., 24.);
  u::DrawHist(hCCIn,2, kRed, "same");

  TLegend *leg = new TLegend(0.5,0.36, 0.97,0.80);
  TLegendEntry *leIn = leg->AddEntry(hCCIn, Form("input cc : %6.2f #pm %6.2f", hCCIn->GetMean(),hCCIn->GetRMS()), "L");
  leIn->SetTextColor(hCCIn->GetLineColor());

  if(hCCOut->GetEntries()>10.) 
  leg->AddEntry(hCCOut, Form("output cc : %6.2f #pm %6.2f", hCCOut->GetMean(),hCCOut->GetRMS()), "L");

  leg->Draw();

  if(c) c->Update();
}

Int_t AliEMCALSuperModule::GetNumberOfCells()
{
  //Returns number of cells in SM
  Int_t ncells=0;
  TList* l = (TList*)GetListOfFolders();
  for(int eta=0; eta<l->GetSize(); eta++) { // cycle on eta row
    TFolder *setEta = dynamic_cast<TFolder*>(l->At(eta));
    if(setEta==0) continue;

    TList* le = (TList*)setEta->GetListOfFolders();
    ncells += le->GetSize();
  }
  return ncells;
}

TList* AliEMCALSuperModule::BookHists()
{
  //Initializes histograms

  gROOT->cd();
  TH1::AddDirectory(1);

  AliEMCALFolder* emcal = (AliEMCALFolder*)GetParent(); 
  Int_t it = emcal->GetIterationNumber();

  new TH1F("00_EffMass",  "effective mass of #gamma,#gamma(m_{#pi^{0}}=134.98 MeV) ", 250,0.0,0.5);
  new TH1F("01_CCInput",  "input CC dist.(MEV) ", 200, 5., 25.);
  new TH1F("02_CCOutput", "output CC dist.(MEV) ", 200, 5., 25.);
  new TH1F("03_MPI0", "mass of #pi_{0} dist. ", 170, 0.05, 0.22);
  new TH1F("04_RESPI0", "resolution at #pi_{0} dist. ", 50, 0.0, 0.05);
  new TH1F("05_XI2/NDF", "#chi^{2} / ndf", 50, 0.0, 5.0);
  new TH1F("06_NPI0", "number of #pi_{0}", 150, 0.0, 1500.);

  TList *l = AliEMCALHistoUtilities::MoveHistsToList(Form("HistsOfSM_%2.2i",fSMNumber), kFALSE);
  AliEMCALHistoUtilities::AddToNameAndTitleToList(l, Form("_%2.2i_It%i",fSMNumber, it), 
						  Form(" SM %2.2i, Iter %i",fSMNumber, it));

  TH1::AddDirectory(0);
  return l;
}
