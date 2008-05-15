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
Revision 1.3  2007/09/12 17:44:22  pavlinov
fixed compilation problem under SuSe Linux

Revision 1.2  2007/09/11 19:38:15  pavlinov
added pi0 calibration, linearity, shower profile

*/ 

//_________________________________________________________________________
// Cell folder which will keep all information 
// about cell(tower) itself
// Initial version was created with TDataSet staf
// TObjectSet -> TFolder; Sep 6, 2007
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALCell.h"
#include "AliEMCALHistoUtilities.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALFolder.h"
#include "AliEMCALSuperModule.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibCoefs.h"
#include "AliEMCALPi0Calibration.h"

#include <cassert>

#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include <TH1.h>
#include <TF1.h>
#include <TNtuple.h>

typedef  AliEMCALHistoUtilities u;

ClassImp(AliEMCALCell)

//______________________________________________________________
AliEMCALCell::AliEMCALCell() : 
TFolder(), 
 fParent(0),fLh(0),
fAbsId(0),fSupMod(0),fModule(0),fPhi(0),fEta(0),fPhiCell(0),fEtaCell(0),fCcIn(0),fCcOut(0),
fFun(0)
{
  //default ctor
}

//______________________________________________________________
AliEMCALCell::AliEMCALCell(const AliEMCALCell& cell) : 
  TFolder(cell.GetName(),cell.GetTitle()), 
  fParent(cell.fParent),fLh(cell.fLh),
  fAbsId(cell.fAbsId),fSupMod(cell.fSupMod),
  fModule(cell.fModule),fPhi(cell.fPhi),
  fEta(cell.fEta),fPhiCell(cell.fPhiCell),
  fEtaCell(cell.fEtaCell),fCcIn(cell.fCcIn),
  fCcOut(cell.fCcOut),fFun(cell.fFun)
{
  //copy ctor
}

//______________________________________________________________
AliEMCALCell::AliEMCALCell(const Int_t absId, const char* title) : 
  TFolder(Form("Cell%4.4i",absId),title), 
 fParent(0),fLh(0),
fAbsId(absId),fSupMod(0),fModule(0),fPhi(0),fEta(0),fPhiCell(0),fEtaCell(0),fCcIn(0),fCcOut(0),
fFun(0)
{
  // Oct 15, 2007  
  AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();
  g->GetCellIndex(fAbsId, fSupMod, fModule, fPhi, fEta);
  g->GetCellPhiEtaIndexInSModule(fSupMod, fModule, fPhi, fEta, fPhiCell, fEtaCell);

} 

//______________________________________________________________
AliEMCALCell::~AliEMCALCell()
{
  // dtor
}

//-------------------------------------------------------------------------------------
void AliEMCALCell::SetCCfromDB(AliEMCALCalibData *ccDb)
{
  // Oct 15, 2007
  if(ccDb == 0) return;
  // fADCchannelEC = fCalibData->GetADCchannel(iSupMod,ieta,iphi);
  // fetaCel-column; fPhiCell- row
  fCcIn = ccDb->GetADCchannel(fSupMod, fEtaCell, fPhiCell); // in GeV

  TH1* h = (TH1*)GetHists()->At(0);
  u::AddToNameAndTitle(h, 0, Form(", cc %5.2f MeV", fCcIn*1.e+3));
}

//______________________________________________________________
void AliEMCALCell::SetCCfromCCTable(AliEMCALCalibCoefs *t)
{
  // Oct 15, 2007
  if(t == 0) return;

  if(fLh == 0) {
    fLh = BookHists();
    Add(fLh);
  }

  AliEMCALCalibCoef *r = t->GetTable(fAbsId);
  if(r && r->fAbsId == fAbsId) {
    fCcIn = r->fCc;
  } else { // something wrong
    if(r) printf(" fAbsId %i : r->absId %i \n", fAbsId, r->fAbsId);
    assert(0);
  }

  TH1* h = (TH1*)GetHists()->At(0);
  u::AddToNameAndTitle(h, 0, Form(", cc %5.2f MeV", fCcIn*1.e+3));
}

//______________________________________________________________
void AliEMCALCell::FillEffMass(const Double_t mgg)
{
  u::FillH1(GetHists(), 0, mgg);
}

//______________________________________________________________
void AliEMCALCell::FillCellNtuple(TNtuple *nt)
{
  if(nt==0) return;
  nt->Fill(fAbsId,fSupMod,fModule,fPhi,fEta,fPhiCell,fEtaCell,fCcIn,fCcOut);
}

//______________________________________________________________
void AliEMCALCell::FitHist(TH1* h, const char* name, const char* opt)
{
  // Oct 15, 2007
  TString optFit(""), sopt(opt);
  sopt.ToUpper();
  if(h==0) return; 
  printf("<I> AliEMCALCell::FitHist : |%s| is started : opt %s\n", h->GetName(), opt);
  TString tit(h->GetTitle());

  TF1 *gausPol2 = 0, *g=0, *bg=0;
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

    gausPol2 = u::GausiPol2(name, 0.00, 0.3, g, bg);
    optFit = "0Q";
    h->Fit(gausPol2, optFit.Data(),"", 0.03, 0.28);
  // Clean up
    delete g;
    delete bg;
    optFit = "0IME+"; // no drwaing at all
    if(tit.Contains("SM") || sopt.Contains("DRAW")) optFit = "IME+";
  } else {
    gausPol2 = (TF1*)h->GetListOfFunctions()->At(0);
    optFit = "IME+";
    printf("<I> Function is defined alredy : %s optFit %s \n", gausPol2->GetTitle(), optFit.Data());
  }
  //  optFit = "IME+";
  h->Fit(gausPol2, optFit.Data(),"", 0.01, 0.28);

  if(optFit.Contains("0") == 0) {
    gStyle->SetOptFit(111);
    u::DrawHist(h,2);
  }
  printf("<I> AliEMCALCell::FitHist : |%s| is ended \n\n", h->GetName());
}

//______________________________________________________________
void AliEMCALCell::FitEffMassHist(const char* opt)
{
  // Oct 15, 2007
  static Double_t mPI0  = 0.13498; // mass of pi0
  static Double_t mPI02 = mPI0*mPI0; // mass**2

  AliEMCALFolder* emcal = AliEMCALPi0Calibration::GetEmcalFolder();
  Int_t it = emcal->GetIterationNumber();

  TH1* h = (TH1*)GetHists()->At(0);

  FitHist(h, GetName(), opt);

  fFun = (TF1*)h->GetListOfFunctions()->At(0);
  if(fFun) {
    Double_t mpi = fFun->GetParameter(1), mpi2 = mpi*mpi;
    Double_t ccTmp = fCcIn * mPI02 / mpi2;
    if(it<=1) { // Jul 16, 2007
      fCcOut = ccTmp;
    } else {
      fCcOut = (ccTmp + fCcIn)/2.;
    }
    printf(" fFun %s | %s : iet %i\n", fFun->GetName(), fFun->GetTitle(), it);
  }
  printf(" %s | fCcIn %6.5f -> % 6.5f <- fCcOut \n", GetTitle(), fCcIn , fCcOut);
}

//______________________________________________________________
void AliEMCALCell::PrintInfo()
{
  // Oct 15, 2007
  printf(" %s %s \n", GetName(), GetTitle());
  if(fLh == 0 ) return;
  TH1* h = (TH1*)GetHists()->At(0);
  TF1 *f = (TF1*)h->GetListOfFunctions()->At(0);
  if(fFun) printf(" fFun : %s | %s \n", fFun->GetName(),  fFun->GetTitle());
  else fFun = f;
  // if(f) f->Dump();
}

//______________________________________________________________
TList* AliEMCALCell::BookHists()
{
  // Oct 15, 2007
  gROOT->cd();
  TH1::AddDirectory(1);

  AliEMCALFolder* emcal = AliEMCALPi0Calibration::GetEmcalFolder();
  Int_t it = emcal->GetIterationNumber();

  new TH1F("01_EffMass", "effective mass of #gamma,#gamma(m_{#pi^{0}}=134.98 MeV) ", 60,0.0,0.3);

  TList *l = AliEMCALHistoUtilities::MoveHistsToList(Form("HistsOfCell%4.4i",fAbsId), kFALSE);
  AliEMCALHistoUtilities::AddToNameAndTitleToList(l, Form("4.4i_It%i",fAbsId,it),
						  Form(" Cell %4.4i, iter. %i",fAbsId, it));

  TH1::AddDirectory(0);
  return l;
}
