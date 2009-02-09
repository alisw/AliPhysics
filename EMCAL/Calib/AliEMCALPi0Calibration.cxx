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
/* History of svn commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
//    This is analysis task for doing EMCAL pi0 calibration
//--  Authors: Aleksei Pavlinov (WSU) 
//--  Feb 17, 2007 - Sep 11, 2007
//--  Pi0 calibration
//--  Recalibration, linearity, shower profile
//--  May 2008 - move to AliAnalysisTaskSE
//--  You should load next library before running:
//    root [0] gSystem->Load("libANALYSIS");
//    root [1] gSystem->Load("libANALYSISalice");
//    root [2] gSystem->Load("libEMCALcalib");
//_________________________________________________________________________

#include "AliEMCALPi0Calibration.h"
#include "AliEMCALFolder.h"
#include "AliEMCALSuperModule.h"
#include "AliEMCALCell.h"
#include "AliEMCALCellInfo.h"
#include "AliEMCALPi0SelectionParam.h"
#include "AliEMCALHistoUtilities.h"
#include "AliEMCALCalibCoefs.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"
#include "AliESD.h"
#include "AliEMCALRecPoint.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliEMCALDigit.h"
#include "AliESDEvent.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <map>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <TBrowser.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TObjString.h>
#include <TArrayI.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TNtuple.h>

using namespace std;

typedef  AliEMCALHistoUtilities u;

AliEMCALGeometry* AliEMCALPi0Calibration::fgEmcalGeo=0;
Int_t             AliEMCALPi0Calibration::fgNmaxCell = 4*12*24*11;

Double_t AliEMCALPi0Calibration::fgDistEff = 3.;
Double_t AliEMCALPi0Calibration::fgW0      = 4.5;
Double_t AliEMCALPi0Calibration::fgSlopePhiShift = 0.01; // ?? just guess

AliEMCALFolder* AliEMCALPi0Calibration::fgEMCAL = 0;
AliEMCALFolder* AliEMCALPi0Calibration::fgEMCALOld = 0;

const Char_t **AliEMCALPi0Calibration::fgAnaOpt=0;
Int_t AliEMCALPi0Calibration::fgNanaOpt = 0; 
enum keyOpt{
  kCORR1,
  kRECALIB,
  kIDEAL,
  kPI0,
  kGAMMA,
  kKINE,
  kPROF,
  kFIT,
 //-- 
  kEND
};

// Enumeration variables
// enum {kUndefined=-1, kLed=0, kBeam=1};
// enum {kPatchResponse=0, kClusterResponse=1};
// enum {kSkip=0,  kSaveToMemory=1};

ClassImp(AliEMCALPi0Calibration)

AliEMCALPi0Calibration::AliEMCALPi0Calibration() :
  AliAnalysisTaskSE(),
  fPmom(0),
  fChain(0),
  fLofHistsPC(0),
  fLofHistsRP(0),
  fLKineVsRP(0),
  fLShowerProfile(0),
  fCellsInfo(0),
  fEmcalPool(0),
  fRunOpts(0),
  fArrOpts(0),
  fKeyOpts(0)
{ 
  // Default constructor - for reading 
}

AliEMCALPi0Calibration::AliEMCALPi0Calibration(const char* name) :
  AliAnalysisTaskSE(name),
  fPmom(0),
  fChain(0),
  fLofHistsPC(0),
  fLofHistsRP(0),
  fLKineVsRP(0),
  fLShowerProfile(0),
  fCellsInfo(0),
  fEmcalPool(0),
  fRunOpts(0),
  fArrOpts(0),
  fKeyOpts(0)
{
  //
  // Constructor. Initialization of pointers
  //
  const Char_t *anaOpt[]={
  "CORR1",   // GetCorrectedEnergyForGamma1(Double_t eRec);
  "RECALIB",
  "IDEAL",
  "PI0",
  "GAMMA",
  "KINE",   // reading kine file
  "PROF",   // Shower profile: phi direction now
  "FIT"     // define parameters : deff, w0 and phislope
  };
  
  fgNanaOpt = sizeof(anaOpt) / sizeof(Char_t*); 
  fgAnaOpt = new const Char_t*[fgNanaOpt];
  for(int i=0; i<fgNanaOpt; i++) fgAnaOpt[i] = anaOpt[i];
}

//AliEMCALPi0Calibration::AliEMCALPi0Calibration(const AliAnalysisTaskSE& obj) :
//  AliAnalysisTaskSE(),
//  fPmom(0),
//  fChain(0),
//  fLofHistsPC(0),
//  fLofHistsRP(0),
//  fLKineVsRP(0),
//  fLShowerProfile(0),
//  fCellsInfo(0),
//  fEmcalPool(0),
//  fRunOpts(0),
//  fArrOpts(0),
//  fKeyOpts(0)
//{ 
//  // Copy constructor - unused
//}
//
//AliEMCALPi0Calibration& AliEMCALPi0Calibration::operator=(const AliAnalysisTaskSE& other)
//{ 
//  // Assignment
// 
//  return *this;
//}

AliEMCALPi0Calibration::~AliEMCALPi0Calibration()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliEMCALPi0Calibration::InitStructure(Int_t it)
{
  //
  // Initialize the common structure of selector
  //
  fgEmcalGeo = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE"); // initialize geometry just once
  fCellsInfo = AliEMCALCellInfo::GetTableForGeometry(fgEmcalGeo);

  if(fRunOpts.Length()>0) CheckRunOpts();

  Int_t key = 0;
  if(GetKeyOptsValue(kGAMMA)) key = 1;
  fLofHistsRP = DefineHistsOfRP("RP", fPmom, key);
  fLofHistsPC = DefineHistsOfRP("PseudoCl", fPmom);

  fEmcalPool = new TFolder("PoolOfEMCAL","");
  if(it <= 1) {
    fgEMCAL     = new AliEMCALFolder(it); // folder for first itteration   
    fEmcalPool->Add(fgEMCAL);
  }
  //if(it<=0) SetName("GammaSel"); // For convinience
  //else      SetName("Pi0Sel");
  if(GetKeyOptsValue(kKINE)) fLKineVsRP = DefineHistsOfKineVsRP("KineVsRP",fPmom, key);
  if(GetKeyOptsValue(kPROF)) fLShowerProfile = DefineHistsForShowerProfile("ProfY", fPmom);
}

void AliEMCALPi0Calibration::CheckRunOpts()
{
  // Check run options
  fRunOpts.ToUpper();
  int nopt = u::ParseString(fRunOpts, fArrOpts);
  printf("<I> AliEMCALPi0Calibration::CheckRunOpts() analyze %i(%i) options : fgNanaOpt %i\n", 
         nopt, fArrOpts.GetEntries(), fgNanaOpt);  
  if(nopt <=0) return;

  fKeyOpts = new TArrayI(fgNanaOpt);
  for(Int_t i=0; i<fgNanaOpt; i++) (*fKeyOpts)[i] = 0;

  for(Int_t i=0; i<fArrOpts.GetEntries(); i++ ) {
    TObjString *o = (TObjString*)fArrOpts.At(i); 

    TString runOpt = o->String();
    Int_t indj=-1;

    for(Int_t j=0; j<fgNanaOpt; j++) {
      TString opt = fgAnaOpt[j];
      if(runOpt.Contains(opt,TString::kIgnoreCase)) {
	indj = j;
        break;
      }
    }
    if(indj<0) {
      printf("<E> option |%s| unavailable **\n", 
      runOpt.Data());
      assert(0);
    } else {
      (*fKeyOpts)[indj] = 1;
      printf("<I> option |%s| is valid : number %i : |%s| \n", 
	     runOpt.Data(), indj, fgAnaOpt[indj]);
    }
  }
}

Int_t AliEMCALPi0Calibration::GetKeyOptsValue(Int_t key)
{
  // Oct 14, 2007
  static Int_t val=0;
  val = 0;
  if(fKeyOpts && key>=0 && key<fKeyOpts->GetSize()) {
    val = fKeyOpts->At(key);
  }
  // printf(" key %i : val %i : opt %s \n", key, val, fgAnaOpt[key]);
  return val;
}

void AliEMCALPi0Calibration::UserExec(const Option_t* /* */)
{
  AliESDEvent * esd = dynamic_cast<AliESDEvent *>(InputEvent());

  if (!esd)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return;
  }
  
  static AliEMCALPi0SelectionParRec* rPar = GetEmcalFolder()->GetPi0SelectionParRow(0);

  static Int_t nEmcalClusters, indOfFirstEmcalRP, nEmcalRP,nEmcalPseudoClusters;
  nEmcalClusters    = esd->GetNumberOfEMCALClusters();
  indOfFirstEmcalRP = esd->GetFirstEMCALCluster();
  u::FillH1(fLofHistsRP, 1, double(indOfFirstEmcalRP));

  static AliRunLoader* rl = 0;
  static Int_t nev = 0; // Temporary - 0nly for reading one file now !!
  
  static AliESDCaloCluster *cl = 0; 
  nEmcalRP = nEmcalPseudoClusters = 0;
  TList *l=0;
  double eDigi=0;   

  TClonesArray lvM1("TLorentzVector", 100);  // for convenience; 100 is enough now 
  TArrayI      indLv(100);                   // index of RP

  static TLorentzVector v, vcl;
  int nrp = 0; // # of RP for gg analysis
  static Double_t erec=0., ecorr=0.0;
  for(int i=indOfFirstEmcalRP; i<indOfFirstEmcalRP+nEmcalClusters; i++) {
    cl = esd->GetCaloCluster(i);
    if(cl->GetClusterType() == AliESDCaloCluster::kEMCALPseudoCluster) {
      nEmcalPseudoClusters++;
      l = fLofHistsPC;
    } else if(cl->GetClusterType() == AliESDCaloCluster::kEMCALClusterv1){
      nEmcalRP++;
      if(fgEMCAL->GetIterationNumber()>1||GetKeyOptsValue(kIDEAL)||GetKeyOptsValue(kRECALIB)||GetKeyOptsValue(kFIT)) {
        AliEMCALRecPoint *rp=0;
        if(GetKeyOptsValue(kFIT) == kFALSE) fgDistEff = -1.; // No fitting ; Sep 4, 2007
        if(GetKeyOptsValue(kIDEAL)) {
          rp = AliEMCALFolder::GetRecPoint(cl, fgEMCAL->GetCCFirst(), 0, fLofHistsRP, fgDistEff, fgW0, fgSlopePhiShift);
	} else {
          rp = AliEMCALFolder::GetRecPoint(cl, fgEMCAL->GetCCFirst(), fgEMCAL->GetCCIn(), fLofHistsRP, fgDistEff, fgW0, fgSlopePhiShift);
	}
        if(GetKeyOptsValue(kPROF)) {
          FillHistsForShowerProfile( GetListShowerProfile(), rp, GetCellsInfo());
	}
	//if(rp->GetPointEnergy()>=rPar->fEOfRpMin && u::GetLorentzVectorFromRecPoint(v, rp)) {
        if(u::GetLorentzVectorFromRecPoint(v, rp)) { // comparing with RP
          if(GetKeyOptsValue(kCORR1)) {
	    erec  = v.Rho();
            ecorr = u::GetCorrectedEnergyForGamma1(erec);
            v.SetRho(ecorr);
            v.SetE(ecorr); // This is gamma
	    //printf("<1> erec %f | ecorr %f \n", erec, ecorr);
	  }
          new(lvM1[nrp]) TLorentzVector(v);
          indLv[nrp] = i;
          nrp++;
	  // Conroling of recalibration
          u::GetLorentzVectorFromESDCluster(vcl, cl);
	  u::FillH1(fLofHistsRP, 11, vcl.P()-v.P());
	  u::FillH1(fLofHistsRP, 12, TMath::RadToDeg()*(vcl.Theta()-v.Theta()));
	  u::FillH1(fLofHistsRP, 13, TMath::RadToDeg()*(vcl.Phi()-v.Phi()));

	  u::FillH1(fLofHistsRP, 16, v.P());
	  u::FillH2(fLofHistsRP, 17, vcl.P(), vcl.P()-v.P());
          l = 0; // no filling
          if(GetKeyOptsValue(kIDEAL) || GetKeyOptsValue(kRECALIB)) l = fLofHistsRP;
        }
        if(rp) delete rp;
      } else { // first iteration
	//        if(cl->E()>=rPar->fEOfRpMin && u::GetLorentzVectorFromESDCluster(v, cl)) {
        if(u::GetLorentzVectorFromESDCluster(v, cl)) { // comparing with RP
	// cut 0.4 GeV may be high ! 
          if(GetKeyOptsValue(kCORR1)) {
	    erec  = v.Rho();
            ecorr = u::GetCorrectedEnergyForGamma1(erec);
            v.SetRho(ecorr);
            v.SetE(ecorr); // This is gamma now
            // printf("<2> erec %f | ecorr %f \n", erec, ecorr);
            // printf("    v.Rho()  %f \n", v.Rho());
	  }
          new(lvM1[nrp]) TLorentzVector(v);
          indLv[nrp] = i;
          nrp++;
          l = fLofHistsRP;
        }
      }
      
    } else {
      printf(" wrong cluster type : %i\n", cl->GetClusterType());
      assert(0);
    }
    u::FillH1(l, 2, double(cl->GetClusterType()));

    u::FillH1(l, 3, double(cl->GetNumberOfDigits()));  
    u::FillH1(l, 4, double(cl->E()));  
    // Cycle on digits (towers)
    Short_t *digiAmpl  = cl->GetDigitAmplitude()->GetArray();
    Short_t *digiTime  = cl->GetDigitTime()->GetArray();
    Short_t *digiAbsId = cl->GetDigitIndex()->GetArray();
    for(int id=0; id<cl->GetNumberOfDigits(); id++) {
      eDigi = double(digiAmpl[id]) / 500.; // See AliEMCALClusterizerv1
      //      if(eDigi <= 0.0) { // sometimes it is happen
      //if(eDigi > 10.0 && cl->GetClusterType() == AliESDCaloCluster::kEMCALClusterv1) {
      // printf(" %i digiAmpl %i : %f \n", id, int(digiAmpl[id]), eDigi);
      //}
      u::FillH1(l, 5, eDigi);
      u::FillH1(l, 6, double(digiTime[id]));
      u::FillH1(l, 7, double(digiAbsId[id]));
      if(int(digiAbsId[id]) >= fgNmaxCell) {
        printf(" id %i :  digiAbsId[id] %i (%i) : %s \n", 
	       id, int(digiAbsId[id]), fgNmaxCell, l->GetName());
      }
    }
  }
  u::FillH1(fLofHistsRP, 0, double(nEmcalRP));
  u::FillH1(fLofHistsPC, 0, double(nEmcalPseudoClusters));

  static TLorentzVector *lv1=0, *lv2=0, lgg;
  for(int i1=0; i1<nrp; i1++){
    lv1 = (TLorentzVector*)lvM1.At(i1);
    u::FillH1(fLofHistsRP, 18, lv1->P());
  }
  static Double_t mgg, pgg;
  mgg = pgg = 0.;
  nrp = lvM1.GetEntriesFast(); 
  if(nrp >= 2) {
    // eff.mass analysis
    for(int i1=0; i1<nrp-1; i1++){
      lv1 = (TLorentzVector*)lvM1.At(i1);
      for(int i2=i1+1; i2<nrp; i2++){
        lv2 = (TLorentzVector*)lvM1.At(i2);
        lgg = (*lv1) + (*lv2);
        mgg = lgg.M();  // eff.mass
        pgg = lgg.P();  // momentum
        u::FillH1(fLofHistsRP, 8, mgg);
 
	if((mgg>=rPar->fMassGGMin && mgg<=rPar->fMassGGMax)) {// pi0 candidates
	  if((pgg>=rPar->fMomPi0Min && pgg>=rPar->fMomPi0Min)) {
            if(fgEMCAL && fgEMCAL->GetIterationNumber()>=1) {
              fgEMCAL->FillPi0Candidate(mgg,esd->GetCaloCluster(indLv[i1]),esd->GetCaloCluster(indLv[i2]));
              u::FillH1(fLofHistsRP, 9, pgg); 
              u::FillH1(fLofHistsRP,10, lv1->P());
              u::FillH1(fLofHistsRP,10, lv2->P());
	    }
	  }
	}
      }
    }
  }

  //  static Int_t fileNumber = 0;
  static TString curFileName;
  //  if(GetKeyOptsValue(kKINE) && nev<fChain->GetEntries()) {
  if(GetKeyOptsValue(kKINE)) {
  // Get galice.root file in current directory
    if(nev%1000==0) {
      printf(" current file |%s|\n", fChain->GetCurrentFile()->GetName());
      curFileName = fChain->GetCurrentFile()->GetName();
      curFileName.ReplaceAll("AliESDs.","galice."); 
    }
    rl = u::InitKinematics(nev, curFileName.Data());
  // Compare kineamtics vs EMCal clusters
    FillHistsOfKineVsRP(fLKineVsRP, rl, lvM1);
  }

  lvM1.Delete();

  if(nEmcalClusters != (nEmcalRP+nEmcalPseudoClusters))
    Info("Process","nEmcalClusters %i : RP %i + PC %i ",nEmcalClusters, nEmcalRP, nEmcalPseudoClusters); 

  nev++;

  return;
}

//
TList *AliEMCALPi0Calibration::DefineHistsOfRP(const char *name,Double_t p,Int_t keyOpt)
{
  //
  // Define histogramms of rec.points
  //
  printf("<I> DefineHistsOfRP :%s : p %f : keyOpt %i \n", name, p, keyOpt);
  Double_t adcChannelEC = 0.0153; // ~15mev per adc count
  Double_t xma = p*1.4, xmi=0.0, step=0.0, xmic=xmi, xmac = xma;
  if(xma<0) xma = 20.;
  Int_t nmax=1000, scale=4, nmaxc = nmax;

  gROOT->cd();
  TH1::AddDirectory(1);
  new TH1F("00_EmcalMultiplicity", "multiplicity of EMCAL  ", 201, -0.5, 200.5); // real and pseudo
  new TH1F("01_IndexOfFirstEmcal", "index of first emcal rec.points ", 201, -0.5, 200.5);

  new TH1F("02_NumberOf", "number of  ", 6, -0.5, 5.5);
  new TH1F("03_NumberOfDigitsIn", "number of digits(towers) in rec.points ", 101, -0.5, 100.5);

  if(keyOpt==1 && p>0.1) {
    if   (p>=100.)  scale=12;
    else if(p>=50.) scale=10;
    else if(p>=30.) scale=8;  
    xma = p + scale*(0.15*TMath::Sqrt(p));
    xmi = p - scale*(0.15*TMath::Sqrt(p));
    step = (xma-xmi)/nmax;
  }

  if(step < 0.0153) {
    nmax = int((xma-xmi) / adcChannelEC)+1;
    xma =  xmi + adcChannelEC*nmax;
  } 
  new TH1F("04_EnergyOf", "energy of ", nmax, xmi, xma);
  nmaxc = nmax; xmic=xmi; xmac = xma;

  nmax = 10000;
  xmi  = adcChannelEC/2.; xma = xmi + adcChannelEC*nmax;
  // All energy(momentum) unit is GeV if don't notice
  new TH1F("05_DigitEnergyIn", "digit energy in ", nmaxc, xmic, xmac);
  new TH1F("06_DigitTimeIn", "digit time in 10ps(0.01ns) ", 1000, 0.0, 3.e+3); // ns/100 = 10 ps
  new TH1F("07_DigitAbsIdIn", "digit abs id in ", fgNmaxCell, -0.5, double(fgNmaxCell)-0.5);
  new TH1F("08_EffMass", "effective mass of #gamma,#gamma(m_{#pi^{0}}=134.9766 MeV)", 100, 0.0, 0.5);
  new TH1F("09_MomOfPi0Candidate", "momentum of #pi^{0} candidates (0.085 <mgg<0.185)", 600, 0.0, 30.0);
  new TH1F("10_MomOfRpPi0Candidate", "momentum of RP for #pi^{0} candidates (0.085 <mgg<0.185)", 600, 0.0, 30.0);
  // Recalibration staf
  Double_t pmax=11., dpmax=0.20;
  if(p>0.1) {
    pmax=p*1.2;
    dpmax *= p;
  }
  new TH1F("11_MomClESD-RpRecalib", "difference of momentum cl(ESD) - rp(Recalib)", 100, -.01, dpmax);
  new TH1F("12_ThetaClESD-RpRecalib", "difference of #theta cl(ESD) - rp(Recalib) in degree", 100, -0.05, +0.05);
  new TH1F("13_PhiClESD-RpRecalib", "difference of #phi cl(ESD) - rp(Recalib) in degree ", 100, -0.05, +0.05);
  // Digi
  new TH1F("14_EDigiRecalib", "energy of digits after recalibration", 2000, 0.0, 20.);
  //  AliEMCALGeometry* g = AliEMCALGeometry::GetInstance();
  new TH1F("15_AbsIdRecalib", "abs Id of digits after recalibration", fgEmcalGeo->GetNCells(),-0.5,Double_t(fgEmcalGeo->GetNCells())-0.5);
  new TH1F("16_EnergyOfRecalibRp_", "energy of recalibrated rec.points", nmaxc, xmic, xmac); // Jul 12, 2007
  new TH2F("17_ShiftRecalib_", "E(clESD) - E(recalib)", 110,0.0, pmax, 50,0.0,dpmax); // Jul 13, 2007
  
  // Corrected staff
  new TH1F("18_EnergyOfCorrectedLV", "energy of corrected LV", nmaxc, xmic, xmac);

  TString st = Form("ListOfHists%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  TList *l = u::MoveHistsToList(st.Data(), kFALSE);
  st = Form("%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  u::AddToNameAndTitleToList(l, st.Data(), st.Data());

  return l;
}

TList* AliEMCALPi0Calibration::DefineHistsOfKineVsRP(const char *name,  Double_t p, Int_t keyOpt)
{
  //
  // Define histogramms for comparing a initial kinematics with rec.points
  //
 printf("<I>  DefineHistsOfKineVsRP :%s : p %f : keyOpt %i \n", name, p, keyOpt);

  gROOT->cd();
  TH1::AddDirectory(1);
  new TH1F("00_hVx",Form("Vx of primary vertex"), 100, -5., +5.);   // 00
  new TH1F("01_hVy",Form("Vy of primary vertex"), 100, -5., +5.);   // 01
  new TH1F("02_hVz",Form("Vz of primary vertex"), 100, -50., +50.); // 02

  //  Double_t adcChannelEC = 0.0153; // ~15mev per adc count
  Double_t xma = p*1.4, xmi=0.0, sig=0.15*TMath::Sqrt(p);
  //  Double_t step=0.0, xmic=xmi, xmac = xma;
  if(xma<0) xma = 20.;
  //Int_t nmax=1000;
  // scale=4, nmaxc = nmax;

  new TH1F("03_hGidPrimar", "Geant Id of primary particle ", 101, 0.5, 100.5);   // 03
  new TH1F("04_hPmomPrim","momentum of primary particle", 100, xmi, xma);        // 04
  new TH1F("05_hEtaPrim","#eta of primary particle ", 200, -1., 1.);             // 05
  new TH1F("06_hPhiPrim","#phi of primary particle ", 63*2, 0.0,TMath::Pi()*2.); // 06
  new TH1F("07_hThetaPrim","#theta of primary particle", 63, 0.0,TMath::Pi());   // 07 
  new TH1F("08_hPhiPrimInDegree","#phi of primary particle in degree", 120, 75.0, 195.0);  // 08
  new TH1F("09_hThetaPrimInDegree","#theta of primary particle in degree", 90, 45., 135.); // 09

  // Particle vs cluster
  new TH1F("10_hE(P)-E(RP)"," E(p) - E(RP) ", 100, -10.*sig, 10.*sig);   // 10

  new TH1F("11_hPvsRpAngleInDegree","angle between P and RP (in degree)", 110, 0.0, 1.0); // 11
  double dphi=0.5; // for fitting
  new TH1F("12_hPvsRpDPhiInDegree","dif(#phi) between P and RP (in degree)", 100, -dphi, dphi); // 12
  new TH1F("13_hPvsRpDThetaInDegree","dif(#theta) between P and RP (in degree)", 100, -0.5, 0.5); // 13

  new TH1F("14_hPvsRpDPhiInDegree","dif(#phi) between P and RP (in degree) right part of SM", 100, -dphi, +dphi); // 14
  new TH1F("15_hPvsRpDPhiInDegree","dif(#phi) between P and RP (in degree) left part of SM", 100, -dphi, dphi); // 15

  new TH1F("16_hPvsRpDThetaInDegree","dif(#theta) between P and RP (even index)", 100, -0.5, 0.5); // 16
  new TH1F("17_hPvsRpDThetaInDegree","dif(#theta) between P and RP (odd  index)", 100, -0.5, 0.5); // 17

  TString st = Form("ListOfHists%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  TList *l = u::MoveHistsToList(st.Data(), kFALSE);
  st = Form("%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  u::AddToNameAndTitleToList(l, st.Data(), st.Data());
  return l;
}

TList *AliEMCALPi0Calibration::DefineHistsForShowerProfile(const char *name, Double_t p)
{
  // Aug 1, 2007 - shower profile business as in testbeam analysis
  //               Phi direction here is Y direction in testbeam 
  printf("<I> DefineHistsForShowerProfile: %s : p %f \n", 
  name, p);

  gROOT->cd();
  TH1::AddDirectory(1);

  // Cell unit
  new TH1F("00_hPHiMean"," mean in #phi direction ", 24*10, 0.0, 24.);
  new TH1F("01_hPHiLog"," mean with log in #phi direction ", 24*10, 0.0, 24.);
  new TH2F("02_XmeanVsXlog", "xmean vs xlog", 240,0.0, 24.,240,0.0, 24.);
  new TH2F("03_XmeanVsXlog", "xmean vs xlog (system of cell with max energy, first half on #eta)", 
  50,-0.5,+0.5, 100,-1.,+1.);
  new TH2F("04_XmeanVsXlog", "xmean vs xlog (system of cell with max energy, second half on #eta)", 
  50,-0.5,+0.5, 100,-1.,+1.);

  new TH2F("05_PhiShowerProfile", " #phi shower profile - cumulative distribution",
  6*25,-3.,3., 201, -0.0025, 1.0025);   // 00; x direction in cell unit

  TString st = Form("L%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  TList *l = u::MoveHistsToList(st.Data(), kFALSE);
  st = Form("%s_P=%5.1f",name,p);
  st.ReplaceAll(" ","");
  u::AddToNameAndTitleToList(l, st.Data(), st.Data());

  return l;
}

void AliEMCALPi0Calibration::FillHistsOfKineVsRP(TList *l, AliRunLoader* rl,  TClonesArray &lvM)
{
  //
  // lvM - array of TLorentzVector's which was cretaef from AliESDCaloCluster's
  //

  if(l==0 || rl==0) return;

  // TNtuple for qucik analysis
  static TNtuple *nt=0;
  if(nt==0) {
    gROOT->cd();
    Int_t bsize = int(1.e+7);
    nt = new TNtuple("angle","angle ntuple for quick analysis",
    "dPhi:dTheta:phiP:thetaP", bsize);
    gROOT->GetListOfBrowsables()->Add(nt);
  }
  static AliStack* st=0;
  static TParticle *p=0;
  static Int_t gid=0, ic=0, pdg=0, i=0;
  gid = ic = pdg = 0;

  st = rl->Stack();
  if(st == 0) return;
  // first primary particle
  p = st->Particle(0);
  if(p == 0) return;

  u::FillH1(l, ic++, p->Vx());
  u::FillH1(l, ic++, p->Vy());
  u::FillH1(l, ic++, p->Vz());

  pdg = p->GetPdgCode();
  //gid = gMC->IdFromPDG(pdg); // gMc should be defined

  u::FillH1(l, ic++, Double_t(gid));
  u::FillH1(l, ic++, p->P());
  u::FillH1(l, ic++, p->Eta());
  u::FillH1(l, ic++, TVector2::Phi_0_2pi(p->Phi()) );
  u::FillH1(l, ic++, p->Theta());
  // In degree
  u::FillH1(l, ic++, TVector2::Phi_0_2pi(p->Phi())*TMath::RadToDeg());
  u::FillH1(l, ic++, p->Theta()*TMath::RadToDeg()); // 09

  if(lvM.GetEntriesFast() == 0) return;
  //
  // Initial kinematics vs Calo clusters
  //
  static TLorentzVector *lv = 0, lvp;
  lvp.SetPxPyPzE(p->Px(),p->Py(),p->Pz(),p->Energy()); // Particle 

  Double_t angle = 0.0, angleMin = 180.0, eDiff=0.0, dPhi=0., dTheta=0., phiP=0.0, thetaP=0.0;
  Int_t indMin = 0;
  phiP   = lvp.Phi() * TMath::RadToDeg();
  if (phiP<81. || phiP>99.)  return;     // cut phi boundaries
  thetaP = lvp.Theta() * TMath::RadToDeg();
  if (thetaP<56. || thetaP>89.)  return; // cut theta boundaries

  for(i=0; i<lvM.GetEntriesFast(); i++) {
    lv = (TLorentzVector*)lvM.At(i);
    angle = lv->Angle(lvp.Vect())*TMath::RadToDeg();
    if(angleMin > angle) {
      angleMin = angle;
      indMin   = i;
    }
  }
  lv = (TLorentzVector*)lvM.At(indMin);
  eDiff   = lvp.E() - lv->E();
  u::FillH1(l, ic++, eDiff); 

  dPhi    = TVector2::Phi_mpi_pi(lvp.Phi()-lv->Phi())*TMath::RadToDeg();
  dTheta  = TVector2::Phi_mpi_pi(lvp.Theta()-lv->Theta())*TMath::RadToDeg();
  u::FillH1(l, ic++, angleMin);
  u::FillH1(l, ic++, dPhi); 
  u::FillH1(l, ic++, dTheta); 

  if       (phiP>=81. && phiP<=90.) {
    u::FillH1(l, 14, dPhi); 
  } else if(phiP> 90. && phiP<=99.) {
    u::FillH1(l, 15, dPhi); 
  }
  // Unfinished - have to get absid of digit with max energy 
  u::FillH1(l, 16, dTheta); 
  u::FillH1(l, 17, dTheta); 
  if(nt) {
    nt->Fill(dPhi,dTheta,phiP,thetaP);
    /*
     tv__tree = (TTree *) gROOT->FindObject("angle");
     tv__tree->Draw("dPhi:phiP","abs(dPhi)<0.5","", 9182, 0);
     tv__tree->Draw("dTheta:thetaP","abs(dTheta)<0.5","", 9182, 0);
     */
  }
}

void AliEMCALPi0Calibration::FillHistsForShowerProfile
(TList *l, AliEMCALRecPoint *rp, AliEMCALCellInfo * t)
{
  // Aug 1, 2007
  if(l==0 || rp==0 || t==0) return;
  // --
  static Double_t xmean=0., xlog=0.;
  static AliEMCALCellIndexes rMax;
  static Int_t phiSize;

  if(rp->GetPointEnergy() < 1.0) return;
  //if(rp->GetPointEnergy() > 8.0) return; // discard merged clusters for pi0 case

  EvalLocalPhiPosition(1.,  rp, t, xmean, phiSize, rMax);
  if(phiSize == 0) return; // just one row in cell directions

  EvalLocalPhiPosition(5.5, rp, t, xlog, phiSize, rMax);
  if(rMax.fIPhi>1.5&&rMax.fIPhi<21.5 && rMax.fIEta>1.5&&rMax.fIEta<46.5) { 
    u::FillH1(l, 0, xmean); 
    u::FillH1(l, 1, xlog); 
    u::FillH2(l, 2, xlog, xmean);
    // Select two central modules
    if((rMax.fIPhim==5 || rMax.fIPhim==6)){
    // Transition to system of cell with max energy
      xmean -= (double(rMax.fIPhi)+0.5);  
      xlog  -= (double(rMax.fIPhi)+0.5);  
      if(rMax.fIEtam>=2 && rMax.fIEtam<=12){ // approximatively first half on eta 
        u::FillH2(l, 3, xlog, xmean);
      } else {// approximatively second half on eta 
        u::FillH2(l, 4, xlog, xmean);
      }
    }
  }
}

void   AliEMCALPi0Calibration::EvalLocalPhiPosition(const Double_t wlog, const AliEMCALRecPoint *rp, const AliEMCALCellInfo* t, Double_t &xcog, Int_t &phiSize, AliEMCALCellIndexes &rMax)
{
  // wlog = 1 - usual center of gravity; >1 - with logarithmic weight.
  // digits - array of digits
  // t      - 
  //==
  // xcog - position of center of gravity in cell unit 
  if(wlog<1.0 || rp==0 || t==0) return;
  xcog = 0;

  static Double_t wtot=0., w=0., edigi=0., e=0., edigiMax=0.;
  static Int_t absid = 0, phiMin=0, phiMax=0;
  static AliEMCALCellIndexes* r=0;

  e = rp->GetPointEnergy();
  wtot = 0.0;
  edigiMax=0.;

  phiMin = 23; phiMax = 0;
  for(Int_t iDigit=0; iDigit<rp->GetMultiplicity(); iDigit++) {
    absid = rp->GetAbsId()[iDigit];
    edigi = rp->GetEnergiesList()[iDigit];
    if(wlog > 1.0)  w = TMath::Max( 0., wlog + TMath::Log(edigi/e));
    else            w = edigi; // just energy

    r     = t->GetTable(absid);
    xcog += w*(Double_t(r->fIPhi) + 0.5);
    wtot += w;
    if(edigi > edigiMax) {
      edigiMax = edigi;
      rMax = (*r);
    }
    if(phiMin > r->fIPhi) phiMin = r->fIPhi; 
    if(phiMax < r->fIPhi) phiMax = r->fIPhi; 
  }
  xcog /= wtot;
  phiSize = phiMax - phiMin;
  //  printf("phiSize %i \n", phiSize); 
}
/* unused now
TList *AliEMCALPi0Calibration::DefineHistsOfTowers(const char *name)
{
  //
  // ESD: towers information was saved to pseudo clusters
  // 
  gROOT->cd();
  TH1::AddDirectory(1);
  new TH1F("00_EmcalMultiplicity", "multiplicity of EMCAL  ", 201, -0.5, 200.5); // number of pseudo RP

  new TH1F("01_EnergyOf", "energy of ", 1000, 0.0, 100.);

  TList *l = u::MoveHistsToList(Form("ListOfHists%s",name, " - ESD"), kFALSE);
  u::AddToNameAndTitleToList(l, name, name);

  return l;
}
*/

void AliEMCALPi0Calibration::FitEffMassHist()
{
  TH1* h = (TH1*)fLofHistsRP->At(8);
  AliEMCALCell::FitHist(h, GetName(), "draw");
}

void AliEMCALPi0Calibration::PrintInfo()
{
  // Service routine
  printf("\n %i Entrie(s) | Option(s) |%s| \n", GetOptsArray().GetEntries(), fRunOpts.Data());  
  for(int i=0; i<fgNanaOpt; i++) {
    if(GetKeyOptsValue(i)) printf(" %i |%s| \n", i, fgAnaOpt[i]);
  }

  TList *l[2] = {fLofHistsPC, fLofHistsRP};
  printf("\n");
  for(int i=0; i<2; i++){
    TH1F *h = (TH1F*)l[i]->At(2);
    printf(" %s \t: %i \n", h->GetTitle(), int(h->GetEntries()));
  }
  printf(" fgDistEff %f fgW0 %f fgSlopePhiShift %f \n", fgDistEff, fgW0, fgSlopePhiShift);
}

void  AliEMCALPi0Calibration::SetMomentum(Double_t p) 
{ // Jul 9, 2007
  fPmom = p;
  //  SetName(Form("%s_p_%f5.1", GetName(), p));
}


AliEMCALFolder*  AliEMCALPi0Calibration::CreateEmcalFolder(const Int_t it)
{
  //
  // Create emcal folder for iteration number it
  //
  AliEMCALFolder* newFolder = new AliEMCALFolder(it); // folder for iteration #it   
  if(it>1) {
    fgEMCALOld = fgEMCAL; 
    AliEMCALCalibCoefs* tabOldOut = fgEMCALOld->GetCCOut();
    AliEMCALCalibCoefs* tabNewIn = new AliEMCALCalibCoefs(*tabOldOut);
    tabNewIn->SetName(AliEMCALFolder::GetCCinName().Data());
    newFolder->Add(tabNewIn);
  } 
  fEmcalPool->Add(newFolder);
  fgEMCAL = newFolder;

  return fgEMCAL;
}

AliEMCALFolder* AliEMCALPi0Calibration::GetEmcalOldFolder(const Int_t nsm)
{
  // Return emcal folder with number nsm
  AliEMCALFolder* folder=0;
  if(fEmcalPool) folder =  (AliEMCALFolder*)fEmcalPool->FindObject(Form("EMCAL_%2.2i",nsm));
  return folder;
}


void AliEMCALPi0Calibration::SetEmcalFolder(AliEMCALFolder* folder)
{
  fgEMCAL = folder;
  fEmcalPool->Add(fgEMCAL);
}

void AliEMCALPi0Calibration::SetEmcalOldFolder(AliEMCALFolder* folder)
{
  fgEMCALOld = folder;
  fEmcalPool->Add(fgEMCALOld);
}

void AliEMCALPi0Calibration::Browse(TBrowser* b)
{
  // What we see at browser
  // if(esd)        b->Add(esd);
  if(fChain)      b->Add(fChain);
  if(fEmcalPool)  b->Add(fEmcalPool);
  if(fgEmcalGeo) b->Add(fgEmcalGeo);
  if(fCellsInfo)  b->Add(fCellsInfo);
  //
  if(fLofHistsPC) b->Add(fLofHistsPC);
  if(fLofHistsRP) b->Add(fLofHistsRP);
  if(fLKineVsRP)  b->Add(fLKineVsRP);
  if(fLShowerProfile) b->Add(fLShowerProfile);
  //  if(u) b->Add(u);
}

Bool_t AliEMCALPi0Calibration::IsFolder() const
{
  if(fLofHistsRP || fEmcalPool) return kTRUE;
  return kFALSE;
}

void AliEMCALPi0Calibration::Save(Int_t ver, const char *optIO)
{ 
  // Aug 3, 2007
  // Save selector to file
  TString dir("/home/pavlinov/ALICE/SHISHKEBAB/RF/CALIB/"); // Root directory for saving
  TString nf=dir;
  if(GetKeyOptsValue(kPROF)) {
    nf += "PROF/PROFILE_";
    nf += ver;
    nf += ".root";
    TFile f(nf.Data(), optIO);
    if(f.IsOpen()) {
      this->Write();
      f.Close();
      printf("<I> Save selectort to file |%s| : optIO %s \n",nf.Data(), optIO);
    } else {
      printf("<W> File %s exits ! Increase version : ver %i !\n", nf.Data(), ver);
    }
  } else {
    printf("No PROF option -> no saving now !\n");
  }
}

AliEMCALPi0Calibration* AliEMCALPi0Calibration::ReadSelector(const char* nf)
{
  // Read selector to file
  AliEMCALPi0Calibration* selector=0;

  TH1::AddDirectory(0);
  TFile f(nf,"READ");
  if(f.IsOpen()) {
    TObject *o = f.Get("AliEMCALPi0Calibration");
    if(o) selector = dynamic_cast<AliEMCALPi0Calibration *>(o);
  }
  // printf("<I> read selector %p : file |%s| \n", selector, nf);
  return selector;
}

void AliEMCALPi0Calibration::ReadAllEmcalFolders()
{
  // Oct 14, 2007
  if(fEmcalPool==0) {
    fEmcalPool = new TFolder("PoolOfEMCAL","");
    for(Int_t it=1; it<=10; it++){
      AliEMCALFolder* fold = AliEMCALFolder::ReadFolder(Form("EMCALFOLDER_It%i_fit.root",it), "READ");
      //      AliEMCALFolder* fold = AliEMCALFolder::Read(Form("EMCALFOLDER_It%i_fit.root",it), "READ");
      if(fold) fEmcalPool->Add(fold);
    }
  }
}

void AliEMCALPi0Calibration::PictVsIterNumber(const Int_t ind, const Int_t nsm)
{
  // Jun 27, 2007 - unfinished; which picture is the best
  if(ind<0 || ind>5) return;
  gROOT->cd();
  TH1::AddDirectory(1);

  Int_t itMax = 10, it=0;
  map <int, const char*> indName;
  indName[0] = "eff.mass";
  indName[3] = "mass of #pi_{0}";
  indName[4] = "resolution of #pi_{0}";
  indName[5] = "chi^{2}/ndf";

  TH1F *hout = new TH1F(indName[ind], indName[ind], itMax, 0.5, double(itMax)+0.5); 
   
  TH1::AddDirectory(0);
  Double_t content, error;
  TList* col = (TList*)fEmcalPool->GetListOfFolders();
  for(Int_t i=0; i<col->GetSize(); i++) { // cycle on EMCAL folders
    AliEMCALFolder* folder = static_cast<AliEMCALFolder*>(col->At(i));
    if(folder==0) continue;
    it = folder->GetIterationNumber();

    AliEMCALSuperModule* sm = folder->GetSuperModule(nsm);
    if(sm==0) continue;

    TList* l = sm->GetHists();
    if(l==0) continue;

    TH1F *hin = (TH1F*)l->At(ind);
    if(hin==0) continue;

    if(ind !=0 ) {
      content = hin->GetMean();
      error   = hin->GetRMS();
    } else {
      sm->FitEffMassHist();
      TF1 *f = (TF1*)hin->GetListOfFunctions()->At(0);
      content = error = -1.;
      if(f) {
	//        content = f->GetParameter(1);
        //error   = f->GetParameter(2);
        content = f->GetParameter(2);
        error   = f->GetParError(2);
      }
    }

    if(content > 0.0) {
      hout->SetBinContent(it, content);
      hout->SetBinError(it, error);
      printf(" it %i content %f +/- %f \n", it, content, error);
    }
  }

  u::DrawHist(hout,2);
  hout->SetMinimum(0.0);

}

TH1F* AliEMCALPi0Calibration::FitHistOfRecPointEnergy(const char *opt)
{
  // Fit hist of rec.point energy
  TH1::AddDirectory(0);

  TString sopt(opt);
  sopt.ToUpper();

  Int_t ind = 4, ind2 = 16;
  if(GetKeyOptsValue(kIDEAL)) {
    ind  = 16; // Jul 12, 2007
    ind2 =  4;
  } else if(GetKeyOptsValue(kCORR1)) {
    ind  = 18;
    ind2 =  4;
  }
  TH1F *hold = (TH1F*)fLofHistsRP->At(ind), *h=0;
  if(hold == 0) return 0;
  if(hold->GetEntries() <10.) return hold;

  if(sopt.Contains("CLONE")) {
    TString newName(Form("C_%s",hold->GetName()));
    h = (TH1F*)hold->Clone(newName.Data());
    printf(" Clone hist %s -> |%s|%s| \n",hold->GetName(),h->GetName(),h->GetTitle()); 
  } else {
    h = hold;
  }
  TH1F* h2 = (TH1F*)fLofHistsRP->At(ind2);

  Double_t xmax = h->GetXaxis()->GetXmax(), xmin = 0.4*xmax;;
  TF1 *g = u::Gausi("fRecPointE", xmin, xmax, h);
  g->SetLineColor(kRed);
  gStyle->SetOptFit(111);

  h->Fit(g,"N","", xmin, xmax);
  printf(" (1) xmin %f : xmax %f \n", xmin, xmax);

  xmin = g->GetParameter(1) - 4.*g->GetParameter(2);
  xmax = g->GetParameter(1) + 4.*g->GetParameter(2);
  h->Fit(g,"Q+","", xmin, xmax);
  u::DrawHist(h2,1, 1, "same",2);
  printf(" (2) xmin %f : xmax %f \n", xmin, xmax);

  return h;
}

TCanvas *AliEMCALPi0Calibration::Linearity(TList *l, int ifun)
{ 
  // Jul 10, 2007
  // Draw picture of EMCal linearity 
  if(l==0) {
    printf("<E> AliEMCALPi0Calibration::Linearity :TList is zero ! Bye ! \n");
    return 0;
  }
  Double_t p[9]={0.5, 1., 3., 5., 10., 20., 30., 50., 100.}, ep[9];
  // see macro rHits.C ->defineSampleFraction
  TH1F* hErecOverEin = new TH1F("hErecOverEin","Ratio E_{rec}/E_{#gamma} vs E_{#gamma}", 101,0.5,101.5);
  // Fill hist
  TArrayD invRat(9), eInvRat(9), erec(9),derec(9), residual(9);
  TArrayD res(9), eres(9); // resolution
  Int_t np=9;
  for(Int_t var=0; var<9; var++){
    TH1F *h = (TH1F*)l->At(var);
    TF1  *f = (TF1*)h->GetListOfFunctions()->At(0);
    Double_t  mean = f->GetParameter(1);
    Double_t emean = f->GetParError(1);

    Int_t bin = hErecOverEin->FindBin(p[var]);
    hErecOverEin->SetBinContent(bin, mean/p[var]);
    hErecOverEin->SetBinError(bin,  emean/p[var]); 
    //
    invRat[var]  = p[var]/mean;
    eInvRat[var] = emean/p[var]*invRat[var];
    erec[var]    = mean;
    derec[var]   = emean;
    // Resolution in %
    res[var]  = 100.*f->GetParameter(2)/p[var];
    eres[var] = 100.*f->GetParError(2)/p[var];
    ep[var] = 0.0;
  }

  TCanvas *c = new TCanvas("Linearity","Linearity", 20,20, 700, 500);
  gStyle->SetOptStat(0);
  if(0) { // E(rec) / E(gamma)
    u::DrawHist(hErecOverEin,2);
    hErecOverEin->SetMinimum(0.984);
    hErecOverEin->SetMaximum(1.001);
  }
  Int_t markerColor=1;
  const char *fun="", *optFit="";
  TF1 *f = 0;
  if(0) {
    if(ifun==-5) {
      fun = "fun5"; optFit="R+";
      f = new TF1(fun,"([0]*(x-7.5)*(x-7.5))*exp([1]+[2]*x+[3]*x*x)+1.", 0.0, 100.);
      f->FixParameter(0, 1.95380e-05);
      //      f->SetParameter(0, 1.95380e-05);
      //      f->SetParameter(1, 1.0);
    } else if(ifun==-4) {
      fun = "fun4"; optFit="R+";
      f = new TF1(fun,"([0]*(x-7.5)*(x-7.5))*([1]+[2]*abs(x)+[3]*(x-[4])*(x-[4]))+1.", 0.0, 100.);
      f->FixParameter(0, 1.95380e-05);
      f->SetParameter(4, 50.);
      //      f = new TF1(fun,"exp(([0]+[1]*x)*(x-[2])*(x-[3]))", 0.0, 100.);
      //f = new TF1(fun,"([0]+[1]*x)*(x-[2])*(x-[3])", 0.0, 100.);
      //f->SetParameter(0, 1.);
      //f->SetParameter(1, 1.);

      //      f->SetParameter(2, 5.);
      //f->SetParLimits(2, 3.,8.);
      //f->SetParameter(3, 10.);
      //f->SetParLimits(3, 9.,15.);
      //f->SetParameter(2, 2.e-4);
    } else if(ifun==-3) {
      fun = "fun3"; optFit="R+";
      f = new TF1(fun,"[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0.0, 11.);
    } else if(ifun==-2) {
      fun = "fun2"; optFit="R+";

      /*
      f = new TF1(fun,"[0]*(x-7.5)+[1]*(x-7.5)*(x-7.5)+[2]", 0.0, 10.1); 
      f->SetParameter(0, 5.98727e-04);
      f->SetParameter(1, 2.12045e-04);
      f->SetParameter(2, 1.);
      */
      f = new TF1(fun,"[0]+[1]*x+[2]*x*x", 9.0, 100.1); 
    } else if(ifun==-1) {
      fun = "fun0"; optFit="R+";
      f = new TF1(fun,"[0]", 0.0, 100.1); 
    } else if(ifun>=2) {
      fun = Form("pol%i",ifun); 
      optFit="+";
    }
    TGraphErrors *gr = u::DrawGraphErrors(np, erec.GetArray(),invRat.GetArray(),derec.GetArray(),eInvRat.GetArray(),
   //    markerColor,21+markerColor,"AP", " Ratio E_{#gamma}/E_{Rec} ", "E_{Rec}  ","   E_{#gamma}/E_{Rec}",
 markerColor,21+markerColor,"AP", " Ratio E_{#gamma}/E_{Rec}^{corr} ", "E_{Rec}^{corr}  ","   E_{#gamma}/E_{Rec}^{corr}",
					  ifun, optFit, fun);
    gr->GetHistogram()->SetAxisRange(0.0,100.);
    //    double xmi=0.999, xma=1.017;
    double xmi=0.995, xma=1.005;
    gr->GetHistogram()->SetMaximum(xma);
    gr->GetHistogram()->SetMinimum(xmi);
    gr->GetHistogram()->SetTitleOffset(1.4,"y");
    if(ifun==0) {
      f = new TF1("fres", "AliEMCALHistoUtilities::EnergyCorrectionForGamma1(x)", 0., 101.); 
      f->Draw("same");
    }
  }
  TLine *line = new TLine(0.0,1.0, 100.,1.0);
  line->Draw();

  if(0) {
    c->Clear();
    for(int i=0; i<9; i++) {
      residual[i] =  100.*(invRat[i] - u::GetCorrectionCoefficientForGamma1(erec[i])); // in percent
      printf(" erec %f : residual %5.3f \n", erec[i], residual[i]);
    }
    markerColor = 2;
    TGraphErrors *gr2=u::DrawGraphErrors(np, erec.GetArray(),residual.GetArray(),derec.GetArray(),eInvRat.GetArray(),
		  markerColor,21+markerColor,"AP"," residual in %, rec.point level", "E_{Rec}  ","  residual in %",
    -1, "", 0);
    gr2->GetHistogram()->SetAxisRange(0.0,100.);
    gr2->GetHistogram()->SetMaximum(0.2);
    gr2->GetHistogram()->SetMinimum(-0.1);
    line = new TLine(0.0,0.0, 101.,0.0);
    line->Draw();
    //TLatex *lat = 
	u::Lat("linearity better 0.2% after correction",20., 0.15, 12, 0.06, 1); 
    //if(lat); //For what is this?, commented due to compilation warnings
  }
  if(1) {
    TString latexName;
    c->Clear();
    gStyle->SetOptFit(111);
    markerColor = 1;
    ifun        = -11;
    f = u::GetResolutionFunction("FRES2", latexName);  
    f->SetNpx(10000);
    f->SetLineColor(kBlack);

    TGraphErrors *gr3=u::DrawGraphErrors(np, p,res.GetArray(), ep, eres.GetArray(),
    markerColor,21+markerColor,"AP"," resolution in %, rec.point level","E_{#gamma} "," resolution in %",
					 ifun, "+", f->GetName());
    gr3->GetHistogram()->SetAxisRange(0.0,101.);
    gr3->GetHistogram()->SetMaximum(14.);
    gr3->GetHistogram()->SetMinimum(0.0);
    gr3->SetMarkerSize(1.5);

    //TLatex *lat = 
	u::Lat(latexName.Data(),82., 11., 12, 0.06, 1);
    //if(lat); //For what is this?, commented due to compilation warnings
    // Exp. data
    TF1 *fexp = new TF1(*f);
    fexp->SetName("fexp");
    fexp->SetParameter(0, 2.168);
    fexp->SetParameter(1, 8.818);
    fexp->SetLineWidth(1);
    fexp->SetLineColor(kRed);
    fexp->SetLineStyle(1);
    fexp->Draw("same");

    TLegend *leg = new TLegend(0.21,0.36, 0.68,0.82);
    leg->AddEntry(gr3,  "MC data", "P");
    leg->AddEntry(f,    "fit of MC data", "L");
    TLegendEntry *ent3 = leg->AddEntry(fexp, "#splitline{fit of exp.data}{FNAL, Nov 2005}", "L");
    ent3->SetTextColor(fexp->GetLineColor());
    leg->Draw();
  }
  c->Update();

  return c;
}

TCanvas *AliEMCALPi0Calibration::DrawKineVsRP(TList *l)
{ 
  //Jul 25, 2007
  if(l==0) {
    printf("<W> AliEMCALPi0Calibration::DrawKineVsRP : TList is zero ! \n");
    return 0;
  }
  TCanvas *c = new TCanvas("KineVsRP","KineVsRP", 20,20, 700, 500);
  gStyle->SetOptStat(1110);
  c->Divide(2,2);

  c->cd(1);
  TH1F* h1 = (TH1F*)l->At(10);
  u::DrawHist(h1,2);

  c->cd(2);
  TH1F* h2 = (TH1F*)l->At(11);
  u::DrawHist(h2,2);

  c->cd(3);
  TH1F* h3 = (TH1F*)l->At(12);
  u::DrawHist(h3,2);
  u::DrawHist((TH1F*)l->At(14), 1,kRed,"same");
  u::DrawHist((TH1F*)l->At(15), 1,kGreen,"same");

  c->cd(4);
  TH1F* h4 = (TH1F*)l->At(13);
  u::DrawHist(h4, 2);
  u::DrawHist((TH1F*)l->At(16), 1,kRed,"same");
  u::DrawHist((TH1F*)l->At(17), 1,kGreen,"same");

  /*
  TH1F* h2 = (TH1F*)l->At(11);
  u::DrawHist(h2,2);
  */

  c->Update();
  return c;
}

TCanvas* AliEMCALPi0Calibration::DrawMeanVsLog(TH2F *h2)
{
  // h - input histogramm : mean vds log coordinates 
  if(h2==0) return 0;

  TCanvas *c = new TCanvas("KineVsRP","KineVsRP", 20,20, 700, 500);

  TH1F *hid1 = new TH1F("h1","h1 title ", h2->GetNbinsX(), 
  h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());  

  gROOT->cd();
  TH1::AddDirectory(1);
  TString newName;
  for(int ix=1; ix<=h2->GetNbinsX();ix++) {
    newName = "hd_"; newName += ix;
    TH1D *hd = h2->ProjectionY(newName.Data(),ix,ix);
    if(hd->Integral()>=4.) {
      //      lhd->Add(hd);
      hid1->SetBinContent(ix, hd->GetMean());
      hid1->SetBinError(ix, hd->GetRMS());
    }
  }
  TF1 *f1 = new TF1("fcg", "0.5*TMath::SinH(x/[0])/TMath::SinH(0.5/[0])", -0.5, 0.5);
  f1->SetParameter(0,0.13);
  f1->SetLineColor(kRed);
  f1->SetLineWidth(1);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  hid1->Fit(f1,"R+");

  u::DrawHist(hid1,2);

//u::DrawHist(h2,2);

  c->Update();
  return c;
}

TCanvas* AliEMCALPi0Calibration::DrawPhiEtaAnglesDistribution(const char* gn)
{ 
  // Aug 6, 2007
  // Proper geometry should be defined already 
  // dTheta distibution has two peaks which coresponding different cells inside module ! 

  TCanvas *c = new TCanvas("Geometry","Geometry", 20,20, 700, 500);
  c->Divide(2,2);

  if(fgEmcalGeo==0) fgEmcalGeo = AliEMCALGeometry::GetInstance(gn);

  gROOT->cd();
  TH1::AddDirectory(1);
  TH1F *hDtheta  = new TH1F("hDtheta","#Delta#theta in one SM", 60, -2.0, +1.0); // in degree
  TH2F *hDtheta2 = new TH2F("hDtheta2","#Delta#theta vs fIEta of cell", 48, -0.5, 47.5, 60,-2.0,+1.0);

  TH1F *hDphi  = new TH1F("hDphi","#Delta#ph in one SM", 2000, -10.0, +10.0); // in degree

  TVector3 vg3;
  AliEMCALCellInfo* t = GetCellsInfo();
  Double_t thetaCell=0., thetaModule=0.;
  Double_t phiCell=0., dphi=0.;

  for(int absid=0; absid<12*24*4; absid++){
    AliEMCALCellIndexes *r = t->GetTable(absid);
    fgEmcalGeo->GetGlobal(absid, vg3);

    thetaCell   = vg3.Theta()*TMath::RadToDeg();
    thetaModule = 90. - 1.5*r->fIEtam;
    hDtheta->Fill(thetaCell - thetaModule);
    hDtheta2->Fill(double(r->fIEta), thetaCell - thetaModule);

    phiCell   = vg3.Phi()*TMath::RadToDeg();
    dphi      = phiCell - 90.;
    hDphi->Fill(dphi);
  }

  c->cd(1);
  u::DrawHist(hDtheta,2);
  c->cd(2);
  u::DrawHist(hDtheta2,2);
  c->cd(3);
  u::DrawHist(hDphi,2);

  c->Update();
  return c;
}

TCanvas* AliEMCALPi0Calibration::DrawDeffVsEnergy()
{
  // Aug 31, 2007 - obsolete
  Double_t p[]={0.5, 1.0, 3., 5., 10., 20.,30.,50.,100.};
  Double_t ep[]={0., 0., 0., 0., 0., 0., 0., 0., 0.};
  // 2 pars
  Double_t  deff[]={9.07515, 10.0835, 11.2032, 11.4229, 12.3578, 13.0332, 13.3281, 13.7910, 14.3319};
  Double_t edeff[]={2.69645e-03, 3.52180e-03, 4.19236e-02, 1.21201e-01, 1.58886e-01, 3.96680e-01, 3.29985e-02, 1.17113e-01, 5.22763e-01};
  //
  Double_t w0[]={3.44679e+00, 3.82714e+00, 4.15035e+0, 4.36650e+00, 4.51511e+00, 4.65590, 4.63289e+00, 4.66568, 4.68125};
  Double_t ew0[]={7.58982e-01, 1.26420e-02, 2.36129e-10, 1.21201e-01, 2.12999e-01, 7.95650e-02, 6.15307e-03, 1.88803e-01, 5.18022e-05};
  int np=sizeof(p)/sizeof(Double_t);
  printf("<I> AliEMCALPi0Calibration::DrawDeffVsEnergy() | np %i \n", np);

  TCanvas *c = new TCanvas("Deff","Deff", 20,20, 700, 500);
  c->Divide(2,1);

  Int_t markerColor=1;
  c->cd(1);
  //TF1 *fdeff= new TF1("fdeff","[0]+[1]*log(x)",0.4, 100.4);
  //if(fdeff); //For what is this?, commented due to compilation warnings
  TGraphErrors *gr = u::DrawGraphErrors(np, p,deff, ep, edeff,
    markerColor,21+markerColor,"AP"," D_{eff} vs E_{#gamma} ","E_{#gamma}         "," D_{eff} in cm ",
					-1, "", 0);
  //					-1, "+", fdeff->GetName());
  gr->GetHistogram()->SetMaximum(15.);
  gPad->SetLogx(1);

  c->cd(2);
  TGraphErrors *gr2 = u::DrawGraphErrors(np, p,w0, ep, ew0,
    markerColor,21+markerColor,"AP"," w_{0} vs E_{#gamma} ","E_{#gamma}         "," w_{0}    ",
    -1, "", 0);

  gr2->GetHistogram()->SetMaximum(5.);
  gPad->SetLogx(1); 

  c->Update();
  return c;
}

TCanvas* AliEMCALPi0Calibration::DrawDeffVsEnergy2(const char *opt)
{
  // Aug 28, 2008 - read pars and pars errors from files
  Double_t p[]={0.5, 1.0, 3., 5., 10., 20.,30.,50.,100.};
  Double_t ep[]={0., 0., 0., 0., 0., 0., 0., 0., 0.};
  // 2 pars
  Double_t deff[9], edeff[9], w0[9], ew0[9]; // max size now 
  TString sopt(opt);

  int np = sizeof(p)/sizeof(Double_t);
  printf("<I> AliEMCALPi0Calibration::DrawDeffVsEnergy2() | np %i \n", np);
  ReadParsDeffAndW0("/data/r22b/ALICE/CALIB/FIT/", deff, edeff, w0, ew0, 1);

  TCanvas *c = new TCanvas("Deff","Deff", 20,20, 700, 500);
  c->Divide(2,1);

  TF1 *fdeff = 0, *fw0 = 0;
  TString optFit(""), funName("");
  if(sopt.Contains("fit1")) {
    fdeff= new TF1("fdeff","[0]+[1]*log(x)",0.1, 101.);
    fdeff->SetLineColor(kRed);
    fdeff->SetLineWidth(1);

    /* good description - "[0]/(1.+exp([1]*x))"
   1  p0           4.82208e+00   9.93617e-04  -0.00000e+00   7.16648e-06
   2  p1          -7.79655e-01   4.07500e-03  -0.00000e+00   2.01009e-03
      better description - "[0]/(1.+exp([1]*(x-[2])))" (like the Woods-Saxon potential)
   1  p0           4.83713e+00   1.00437e-03  -6.98984e-10   4.90371e-07
   2  p1          -2.77970e-01   3.35587e-03  -1.79239e-09   2.67312e-07
   3  p2           4.41116e+00   8.72191e-02   1.82791e-04   1.55643e-05
     */
    fw0= new TF1("fw0","[0]/(1.+exp([1]*(x+[2])))",0.1, 101.);
    fw0->SetLineColor(kRed);
    fw0->SetLineWidth(1);
    fw0->SetParameter(0, 4.8);
    fw0->SetParameter(1, -2.77970e-01);
    fw0->SetParameter(2,  4.41116);
    optFit = "+";
  }
  if(fdeff) funName =  fdeff->GetName();

  Int_t markerColor=1;
  c->cd(1);
  gStyle->SetOptFit(111);
  TGraphErrors *gr = u::DrawGraphErrors(np, p,deff, ep, edeff,
    markerColor,21+markerColor,"AP"," D_{eff} vs E_{#gamma} ","E_{#gamma}         "," D_{eff} in cm ",
					-1, optFit.Data(), funName.Data());
  gr->GetHistogram()->SetMaximum(15.);
  gPad->SetLogx(1);
  TLegend *leg1 = new TLegend(0.12,0.76, 0.70,0.90);
  TLegendEntry *le1 = leg1->AddEntry(fdeff, Form("%s",fdeff->GetTitle()), "lp");
  //TLegendEntry *le1 = leg1->AddEntry(fdeff, Form("%4.2f+%4.2f*log(E_{#gamma})",
  //fdeff->GetParameter(0),fdeff->GetParameter(1)), "lp");
  le1->SetTextColor(fdeff->GetLineColor());
  leg1->Draw();

  c->cd(2);
  gStyle->SetOptFit(111);
  if(fw0) funName =  fw0->GetName();
  TGraphErrors *gr2 = u::DrawGraphErrors(np, p,w0, ep, ew0,
    markerColor,21+markerColor,"AP"," w_{0} vs E_{#gamma} ","E_{#gamma}         "," w_{0}    ",
					-1, optFit.Data(), funName.Data());

  gr2->GetHistogram()->SetMaximum(5.);

  TLegend *leg2 = new TLegend(0.17,0.6, 0.99,0.72);
  TLegendEntry *le2 = leg2->AddEntry(fw0, Form("%s",fw0->GetTitle()), "lp");
  //TLegendEntry *le2 = leg2->AddEntry(fw0, Form("#frac{%4.2f}{1.+exp(%4.2f*(x+%4.2f)}",
  //fw0->GetParameter(0),fw0->GetParameter(1),fw0->GetParameter(2)), "lp");
  le2->SetTextColor(fw0->GetLineColor());
  leg2->Draw();
  //gPad->SetLogx(1); 

  c->Update();
  return c;
}

void AliEMCALPi0Calibration::ReadParsDeffAndW0
(const char *dirName, double *deff, double *edeff, double *w0, double *ew0, const Int_t pri)
{
  // read pars and W0
  int strategy = 0, itmp=0;
  char line[100];
  for(int var=11; var<=19; var++){
    int ind = var -11;
    ifstream fin;
    TString fname = Form("%s/fitVar%iStrategy%i.txt", dirName, var, strategy);
    //printf(" open file %s \n", fname.Data());
    fin.open(fname.Data());
   
    for(int i=1; i<=2; i++) { // skip to lines
      fin.getline(line,100).eof();
      if(pri>=2) printf("%s \n", line);
    }

    fin >> itmp >> deff[ind] >> edeff[ind];
    fin >> itmp >> w0[ind]   >> ew0[ind];
    if(pri>=1) 
    printf(" %i def %f+/-%f : def %f+/-%f\n", ind, deff[ind],edeff[ind],w0[ind],ew0[ind]); 
    fin.getline(line,100).eof(); // skip last line
    fin.close();
  }
}

TCanvas* AliEMCALPi0Calibration::DrawSpaceResolution()
{
  // Sep 4, 2007;
  // Space resolution vs optimised space resolution
  Double_t p[]={0.5, 1.0, 3., 5., 10., 20.,30.,50.,100.};
  Double_t ep[]={0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t  spAng[]={0.1877, 0.1351, 0.08144, 0.06331, 0.05384, 0.04876, 0.04706, 0.04656, 0.04726};
  Double_t rmsSpAng[]={0.1033, 0.07685, 0.05235, 0.03992, 0.04012, 0.04257, 0.03472, 0.02814, 0.02784};
  Double_t  spAngOpt[]={0.1733, 0.1311, 0.07961, 0.06401, 0.05347, 0.04618, 0.04288, 0.04, 0.03802};
  Double_t rmsSpAngOpt[]={0.09476, 0.07472, 0.05011, 0.04242, 0.04075, 0.04304, 0.03545, 0.02744, 0.02593};

  int np=sizeof(p)/sizeof(Double_t);
  printf("<I> AliEMCALPi0Calibration::DrawSpaceResolution() | np %i \n", np);

  Double_t* eSpAng    = new Double_t[np];
  Double_t* eSpAngOpt = new Double_t[np];
  Double_t cc=TMath::Sqrt(8000.), cc2 = 1000.*TMath::DegToRad();
  for(int i=0; i<np; i++){
    spAng[i]       *= cc2;
    rmsSpAng[i]    *= cc2;
    spAngOpt[i]    *= cc2;
    rmsSpAngOpt[i] *= cc2;

    eSpAng[i]    = spAng[i]/cc;
    eSpAngOpt[i] = spAngOpt[i]/cc;
  }

  TCanvas *c = new TCanvas("Deff","Deff", 20,20, 700, 500);
  //c->Divide(2,1);

  c->cd(1);
  gStyle->SetOptFit(111);

  TGraphErrors *gr = u::DrawGraphErrors(np, p,spAng, ep, eSpAng,
  kBlack, 21,"AP","Angle resolution (mrad) vs E_{#gamma} ","E_{#gamma}       "," anlgle resolution (mrad) ",
					-1, "", 0);
  gr->GetHistogram()->SetMaximum(4.);
  gr->GetHistogram()->SetMinimum(0.6);
  gPad->SetLogx(1);
  gPad->SetLogy(1);

  TF1 *fang = new TF1("fang","[0]+[1]/sqrt(x)",0.1, 101.);
  fang->SetLineColor(kRed);
  fang->SetLineWidth(1);

  TGraphErrors *gr2 = u::DrawGraphErrors(np, p,spAngOpt, ep, eSpAngOpt,
  kRed,22,"P"," space vs E_{#gamma} ","E_{#gamma}         "," anlgle resolution (mrad) ",
					 -1, "+", fang->GetName());
  TLegend *leg1 = new TLegend(0.40,0.57, 0.92,0.87);
  TLegendEntry *le1 = leg1->AddEntry(gr, "Initial angle resolution", "p");
  le1->SetTextColor(gr->GetLineColor());

  TLegendEntry *le2 = leg1->AddEntry(gr2, "Optimized angle resolution", "p");
  le2->SetTextColor(gr2->GetLineColor());

  TLegendEntry *le3 = leg1->AddEntry(fang, 
  Form("%3.2f+#frac{%4.2f}{#sqrt{E}}",fang->GetParameter(0),fang->GetParameter(1)), "lp");
  le3->SetTextColor(fang->GetLineColor());

  leg1->Draw();

  c->Update();

  return c;
}

void AliEMCALPi0Calibration::ResetAllListOfHists()
{
  // Reset all list of hits
  u::ResetListOfHists(fLofHistsPC);
  u::ResetListOfHists(fLofHistsRP);
  u::ResetListOfHists(fLKineVsRP);
  u::ResetListOfHists(fLShowerProfile);
}

void AliEMCALPi0Calibration::ReloadChain(Long64_t entry)
{ 
  // Oct 14, 2007 - unused now
  if(fChain) {
    fChain->LoadTree(entry);
  }
}

void AliEMCALPi0Calibration::GetInitialParsForFit(const Int_t var, Double_t &deff, Double_t &w0, Double_t &phislope, const int phiCase)
{
  //  int phiCase=0; // 0 or 1
  phislope = 0.001;
  if(var==11)         { // stay here
    deff = 9.07515;
    w0   = 3.44679;
  }  else if(var==12) {
    deff = 1.00835e+01;
    w0   = 3.82714e+00;
  }  else if(var==13) {
    deff = 1.12032e+01;
    w0   = 4.15035e+00;
  }  else if(var==14) {
    deff = 1.14229e+01;
    w0   = 4.36650;
  }  else if(var==15) {
    deff = 1.21361e+01; 
    w0   = 4.72875;
  }  else if(var==16) {
    deff = 1.28096e+01;
    w0   = 4.81753e+00;
  }  else if(var==17) {
    deff = 1.33281e+01; 
    w0   = 4.63289e+00;
  }  else if(var==18) {
    deff = 1.37910e+01; 
    w0   = 4.66568;
  }  else if(var==19) {// Aug 15, 2007
    switch (phiCase) {
    case 0: // phislope was defined than deff and w0 were fixed
      deff = 1.43319e+01;
      w0   = 4.69279;
      phislope = 2.41419e-04;
      break;
    case 1: // all parameters free
      deff = 1.46327e+01;
      w0   = 4.68125;
      phislope = 8.95559e-04;
      break;
    }
  } else { 
    printf("<E> var % i -> no definition ! \n", var);
    assert(0);
  }
  printf("<I> AliEMCALPi0Calibration::GetInitialParForFit()\n deff %f w0 %f phislope %f \n", 
	 deff,w0, phislope);
}
