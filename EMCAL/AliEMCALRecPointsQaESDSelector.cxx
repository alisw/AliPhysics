/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

/* $Log$ */

//_________________________________________________________________________
// This is selector for making a set of histogramms from ESD for rec.points
//*--  Authors: Aleksei Pavlinov (WSU) 
//*--  Feb 17, 2007 - May 23, 2007
//*--  Pi0 calibration
//_________________________________________________________________________

#include "AliEMCALRecPointsQaESDSelector.h"
#include "AliEMCALFolder.h"
#include "AliEMCALSuperModule.h"
#include "AliEMCALCell.h"
#include "AliEMCALPi0SelectionParam.h"
#include "AliEMCALHistoUtilities.h"
#include "AliEMCALCalibCoefs.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"
#include "AliESD.h"
#include "AliEMCALRecPoint.h"

#include <assert.h>
#include <map>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TArrayS.h>
#include <TBrowser.h>
#include <TDataSet.h>

using namespace std;

typedef  AliEMCALHistoUtilities u;

ClassImp(AliEMCALRecPointsQaESDSelector)

Int_t  nmaxCell = 4*12*24*11;

AliEMCALRecPointsQaESDSelector::AliEMCALRecPointsQaESDSelector() :
  AliSelector(),
  fChain(0),
  fLofHistsPC(0),
  fLofHistsRP(0),
  fEmcalPool(0),
  fEMCAL(0),
  fEMCALOld(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliEMCALRecPointsQaESDSelector::~AliEMCALRecPointsQaESDSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliEMCALRecPointsQaESDSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelector::Begin(tree);
}

void AliEMCALRecPointsQaESDSelector::Init(TTree* tree)
{
  // Init function

  AliSelector::Init(tree);
}

void AliEMCALRecPointsQaESDSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  printf("**** <I> Slave Begin \n **** ");

  AliSelector::SlaveBegin(tree);
}

void AliEMCALRecPointsQaESDSelector::InitStructure(Int_t it)
{
  AliEMCALGeometry::GetInstance("SHISH_TRD1_CURRENT_2X2"); // initialize geometry just once
  //AliEMCALGeometry::GetInstance(""); // initialize geometry just once

  fLofHistsPC = DefineHistsOfRP("PseudoCl");
  fLofHistsRP = DefineHistsOfRP("RP");

  fEmcalPool = new TObjectSet("PoolOfEMCAL");
  if(it == 1) {
    fEMCAL     = new AliEMCALFolder(it); // folder for first itteration   
    fEmcalPool->Add(fEMCAL);
  }
}

Bool_t AliEMCALRecPointsQaESDSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }
  static pi0SelectionParam* rPar = GetEmcalFolder()->GetPi0SelectionParRow(0);

  static Int_t nEmcalClusters, indOfFirstEmcalRP, nEmcalRP,nEmcalPseudoClusters;
  nEmcalClusters    = fESD->GetNumberOfEMCALClusters();
  indOfFirstEmcalRP = fESD->GetFirstEMCALCluster();
  u::FillH1(fLofHistsRP, 1, double(indOfFirstEmcalRP));

  static AliESDCaloCluster *cl = 0; 
  nEmcalRP = nEmcalPseudoClusters = 0;
  TList *l=0;
  double eDigi=0;   

  TClonesArray lvM1("TLorentzVector", 100);  // for convenience; 100 is enough now 
  TArrayI      indLv(100);                   // index of RP

  static TLorentzVector v, vcl;
  int nrp = 0; // # of RP for gg analysis
  for(int i=indOfFirstEmcalRP; i<indOfFirstEmcalRP+nEmcalClusters; i++) {
    cl = fESD->GetCaloCluster(i);
    if(cl->GetClusterType() == AliESDCaloCluster::kPseudoCluster) {
      nEmcalPseudoClusters++;
      l = fLofHistsPC;
    } else if(cl->GetClusterType() == AliESDCaloCluster::kClusterv1){
      nEmcalRP++;
      if(fEMCAL->GetIterationNumber()>1) {
        AliEMCALRecPoint *rp = AliEMCALFolder::GetRecPoint(cl, fEMCAL->GetCCFirst(), fEMCAL->GetCCIn(), fLofHistsRP);
	//        if(rp->GetPointEnergy()>=rPar->eOfRpMin && u::GetLorentzVectorFromRecPoint(v, rp)) {
        if(u::GetLorentzVectorFromRecPoint(v, rp)) { // comparing with RP
          new(lvM1[nrp]) TLorentzVector(v);
          indLv[nrp] = i;
          nrp++;
	  // Conroling of recalibration
          u::GetLorentzVectorFromESDCluster(vcl, cl);
	  u::FillH1(fLofHistsRP, 11, vcl.P()-v.P());
	  u::FillH1(fLofHistsRP, 12, TMath::RadToDeg()*(vcl.Theta()-v.Theta()));
	  u::FillH1(fLofHistsRP, 13, TMath::RadToDeg()*(vcl.Phi()-v.Phi()));
          l = 0; // no filling
        }
        delete rp;
      } else { // first iteration
	//        if(cl->E()>=rPar->eOfRpMin && u::GetLorentzVectorFromESDCluster(v, cl)) {
        if(u::GetLorentzVectorFromESDCluster(v, cl)) { // comparing with RP
	// cut 0.4 GeV may be high ! 
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
      //if(eDigi > 10.0 && cl->GetClusterType() == AliESDCaloCluster::kClusterv1) {
      // printf(" %i digiAmpl %i : %f \n", id, int(digiAmpl[id]), eDigi);
      //}
      u::FillH1(l, 5, eDigi);
      u::FillH1(l, 6, double(digiTime[id]));
      u::FillH1(l, 7, double(digiAbsId[id]));
      if(int(digiAbsId[id]) >= nmaxCell) {
        printf(" id %i :  digiAbsId[id] %i (%i) : %s \n", 
	       id, int(digiAbsId[id]), nmaxCell, l->GetName());
      }
    }
  }
  u::FillH1(fLofHistsRP, 0, double(nEmcalRP));
  u::FillH1(fLofHistsPC, 0, double(nEmcalPseudoClusters));

  static TLorentzVector *lv1=0, *lv2=0, lgg;
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
 
	if((mgg>=rPar->massGGMin && mgg<=rPar->massGGMax)) {// pi0 candidates
	  if((pgg>=rPar->momPi0Min && pgg>=rPar->momPi0Min)) {
            fEMCAL->FillPi0Candidate(mgg,fESD->GetCaloCluster(indLv[i1]),fESD->GetCaloCluster(indLv[i2]));
            u::FillH1(fLofHistsRP, 9, pgg); 
            u::FillH1(fLofHistsRP,10, lv1->P());
            u::FillH1(fLofHistsRP,10, lv2->P());
	  }
	}
      }
    }
  }

  lvM1.Delete();

  if(nEmcalClusters != (nEmcalRP+nEmcalPseudoClusters))
    Info("Process","nEmcalClusters %i : RP %i + PC %i ",nEmcalClusters, nEmcalRP, nEmcalPseudoClusters); 

  return kTRUE;
}

void AliEMCALRecPointsQaESDSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  // fOutput->Add(fMultiplicity);
}

void AliEMCALRecPointsQaESDSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();
  /*
  fMultiplicity = dynamic_cast<TH1F*> (fOutput->FindObject("fMultiplicity"));

  if (!fMultiplicity)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histogram not available %p", (void*) fMultiplicity));
    return;
  }

  new TCanvas;
  fMultiplicity->DrawCopy();
  */
  PrintInfo();
}
//
TList *AliEMCALRecPointsQaESDSelector::DefineHistsOfRP(const char *name)
{
  double ADCchannelEC = 0.00305; // ~3mev per adc count

  gROOT->cd();
  TH1::AddDirectory(1);
  new TH1F("00_EmcalMultiplicity", "multiplicity of EMCAL  ", 201, -0.5, 200.5); // real and pseudo
  new TH1F("01_IndexOfFirstEmcal", "index of first emcal rec.points ", 201, -0.5, 200.5);

  new TH1F("02_NumberOf", "number of  ", 6, -0.5, 5.5);
  new TH1F("03_NumberOfDigitsIn", "number of digits(towers) in rec.points ", 101, -0.5, 100.5);
  new TH1F("04_EnergyOf", "energy of ", 1000, 0.0, 100.);

  int nmax=10000;
  double xmi = ADCchannelEC/2., xma = xmi + ADCchannelEC*nmax;
  // All energy(momentum) unit is GeV if don't notice
  new TH1F("05_DigitEnergyIn", "digit energy in ", nmax, xmi, xma);
  new TH1F("06_DigitTimeIn", "digit time in 10ps(0.01ns) ", 1000, 0.0, 3.e+3); // ns/100 = 10 ps
  new TH1F("07_DigitAbsIdIn", "digit abs id in ", nmaxCell, -0.5, double(nmaxCell)-0.5);
  new TH1F("08_EffMass", "effective mass of #gamma,#gamma(m_{#pi^{0}}=134.9766 MeV)", 100, 0.0, 0.5);
  new TH1F("09_MomOfPi0Candidate", "momentum of #pi^{0} candidates (0.085 <mgg<0.185)", 600, 0.0, 30.0);
  new TH1F("10_MomOfRpPi0Candidate", "momentum of RP for #pi^{0} candidates (0.085 <mgg<0.185)", 600, 0.0, 30.0);
  // Recalibration staf
  new TH1F("11_MomClESD-RpRecalib", "difference of momentum cl(ESD) - rp(Recalib)", 200, -1., +1.);
  new TH1F("12_ThetaClESD-RpRecalib", "difference of #theta cl(ESD) - rp(Recalib) in degree", 100, -0.05, +0.05);
  new TH1F("13_PhiClESD-RpRecalib", "difference of #phi cl(ESD) - rp(Recalib) in degree ", 100, -0.05, +0.05);
  // Digi
  new TH1F("14_EDigiRecalib", "energy of digits after recalibration", 2000, 0.0, 20.);
  AliEMCALGeometry* g = AliEMCALGeometry::GetInstance();
  new TH1F("15_AbsIdRecalib", "abs Id of digits after recalibration", g->GetNCells(),-0.5,Double_t(g->GetNCells())-0.5);
  
  TList *l = u::MoveHistsToList(Form("ListOfHists%s",name), kFALSE);
  u::AddToNameAndTitleToList(l, name, name);

  return l;
}

/* unused now
TList *AliEMCALRecPointsQaESDSelector::DefineHistsOfTowers(const char *name)
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

void AliEMCALRecPointsQaESDSelector::FitEffMassHist()
{
  TH1* h = (TH1*)fLofHistsRP->At(8);
  AliEMCALCell::FitHist(h, GetName(), "draw");
}

void AliEMCALRecPointsQaESDSelector::PrintInfo()
{
  // Service routine
  TList *l[2] = {fLofHistsPC, fLofHistsRP};
  printf("\n");
  for(int i=0; i<2; i++){
    TH1F *h = (TH1F*)l[i]->At(2);
    printf(" %s \t: %i \n", h->GetTitle(), int(h->GetEntries()));
  }
}

AliEMCALFolder*  AliEMCALRecPointsQaESDSelector::CreateEmcalFolder(const Int_t it)
{
  AliEMCALFolder* newFolder = new AliEMCALFolder(it); // folder for iteration #it   
  if(it>1) {
    fEMCALOld = fEMCAL; 
    AliEMCALCalibCoefs* tabOldOut = fEMCALOld->GetCCOut();
    AliEMCALCalibCoefs* tabNewIn = new AliEMCALCalibCoefs(*tabOldOut);
    tabNewIn->SetName(AliEMCALFolder::fgkCCinName.Data());
    newFolder->Add(tabNewIn);
  } 
  fEmcalPool->Add(newFolder);
  fEMCAL = newFolder;

  return fEMCAL;
}

AliEMCALFolder* AliEMCALRecPointsQaESDSelector::GetEmcalOldFolder(const Int_t nsm)
{
  AliEMCALFolder* folder=0;
  if(fEmcalPool) folder =  (AliEMCALFolder*)fEmcalPool->FindByName(Form("EMCAL_%2.2i",nsm));
  return folder;
}


void AliEMCALRecPointsQaESDSelector::SetEmcalFolder(AliEMCALFolder* folder)
{
  fEMCAL = folder;
  fEmcalPool->Add(fEMCAL);
}

void AliEMCALRecPointsQaESDSelector::SetEmcalOldFolder(AliEMCALFolder* folder)
{
  fEMCALOld = folder;
  fEmcalPool->Add(fEMCALOld);
}

void AliEMCALRecPointsQaESDSelector::Browse(TBrowser* b)
{
  if(fESD)        b->Add(fESD);
  if(fLofHistsPC) b->Add(fLofHistsPC);
  if(fLofHistsRP) b->Add(fLofHistsRP);
  if(fChain)      b->Add(fChain);
  if(fEmcalPool)  b->Add(fEmcalPool);
  //  if(u) b->Add(u);
}

Bool_t AliEMCALRecPointsQaESDSelector::IsFolder() const
{
  if(fESD || fLofHistsRP || fEmcalPool) return kTRUE;
  return kFALSE;
}

void AliEMCALRecPointsQaESDSelector::ReadAllEmcalFolders()
{
  if(fEmcalPool==0) {
    fEmcalPool = new TObjectSet("PoolOfEMCAL");
    for(Int_t it=1; it<=10; it++){
      AliEMCALFolder* fold = AliEMCALFolder::Read(Form("EMCALFOLDER_It%i_fit.root",it), "READ");
      if(fold) fEmcalPool->Add(fold);
    }
  }
}

void AliEMCALRecPointsQaESDSelector::PictVsIterNumber(const Int_t ind, const Int_t nsm)
{
  // Jun 27, 2007 - unfinished; which picture is the best
  if(ind<0 || ind>5) return;
  gROOT->cd();
  TH1::AddDirectory(1);

  Int_t itMax = 10, it=0;
  map <int, char*> indName;
  indName[0] = "eff.mass";
  indName[3] = "mass of #pi_{0}";
  indName[4] = "resolution of #pi_{0}";
  indName[5] = "chi^{2}/ndf";

  TH1F *hout = new TH1F(indName[ind], indName[ind], itMax, 0.5, double(itMax)+0.5); 
   
  TH1::AddDirectory(0);
  Double_t content, error;
  for(Int_t i=0; i<fEmcalPool->GetListSize(); i++) { // cycle on EMCAL folders
    AliEMCALFolder* folder = static_cast<AliEMCALFolder*>(fEmcalPool->At(i));
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
