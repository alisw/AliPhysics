//____________________________________________________________________
//
// $Id$
//
// Script to get the various cross sections, energy loss, and ranges
// of a particular particle type in a particular medium. 
//
// This script should be compiled to speed it up.  
// 
// It creates a tree on the current output file, with the relevant
// information. 
//
// Note, that VMC _must_ be the TGeant3TGeo VMC. 
//
#include <TArrayF.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoMedium.h>
#include <TGeant3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH2.h>
#include <TLegend.h>
#include <AliRun.h>
#include <TCanvas.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
/** @defgroup FMD_xsec_script X-section script 
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** @ingroup FMD_xsec_script
 */
struct Mech 
{
  char* fName;
  char* fTitle;
  char* fUnit;
  int   fColor;
  int   fStyle;
};

Mech fgMechs[] = 
  {{ "HADF","FLUKA total hadronic x-section",        "cm^{-1}",1,1},
   { "INEF","FLUKA hadronic inelastic x-section",    "cm^{-1}",2,1},
   { "ELAF","FLUKA hadronic elastic x-section",      "cm^{-1}",3,1},
   { "HADG","GHEISHA total hadronic x-section",      "cm^{-1}",4,1},
   { "INEG","GHEISHA hadronic inelastic x-section",  "cm^{-1}",6,1},
   { "ELAG","GHEISHA hadronic elastic x-section",    "cm^{-1}",7,1},
   { "FISG","GHEISHA nuclear fission x-section",     "cm^{-1}",1,2},
   { "CAPG","GHEISHA neutron capture x-section",     "cm^{-1}",2,2},
   { "LOSS","stopping power",                        "MeV/cm", 3,2},
   { "PHOT","photoelectric x-section",               "cm^{-1}",4,2},
   { "ANNI","positron annihilation x-section",       "cm^{-1}",6,2},
   { "COMP","Compton effect x-section",              "cm^{-1}",7,2},
   { "BREM","bremsstrahlung x-section",              "cm^{-1}",1,3},
   { "PAIR","photon and muon direct- pair x-section","cm^{-1}",2,3},
   { "DRAY","delta-rays x-section",                  "cm^{-1}",3,3},
   { "PFIS","photo-fission x-section",               "cm^{-1}",4,3},
   { "RAYL","Rayleigh scattering x-section",         "cm^{-1}",6,3},
   { "MUNU","muon-nuclear interaction x-section",    "cm^{-1}",7,3},
   { "RANG","range",                                 "cm",     1,4},
   { "STEP","maximum step",                          "cm",     2,4},
   { 0, 0, 0, 0, 0}};

//____________________________________________________________________
/** @ingroup FMD_xsec_script
 */
struct MechValue 
{
  MechValue(const Mech& m, size_t n) 
    : fMech(m), fValues(n), fStatus(0), fELoss(0)
  {
    fGraph.SetName(fMech.fName);
    fGraph.SetTitle(fMech.fTitle);
    fGraph.GetXaxis()->SetTitle("#beta#gamma");
    fGraph.GetYaxis()->SetTitle(fMech.fUnit);
    fGraph.SetLineColor(fMech.fColor);
    fGraph.SetLineStyle(fMech.fStyle);
    fGraph.SetLineWidth(2);
    std::string name(fMech.fName);
    if (name != "LOSS") return;
    
    fELoss = new TGraphErrors(n);
    fELoss->SetName("ELOSS");
    fELoss->SetTitle(fMech.fTitle);
    fELoss->GetXaxis()->SetTitle("#beta#gamma");
    fELoss->GetYaxis()->SetTitle(fMech.fUnit);
    fELoss->SetLineColor(fMech.fColor);
    fELoss->SetLineStyle(fMech.fStyle);
    fELoss->SetLineWidth(2);
    fELoss->SetFillColor(fMech.fColor+1);
    fELoss->SetFillStyle(3001);
  }
  bool Get(TGeant3* mc, std::vector<float>& t, 
	   std::vector<float>& c, int medNo, int pidNo, 
	   float mass) 
  {
    int ixst;
    mc->Gftmat(medNo, pidNo, fMech.fName, t.size(), 
	       &(t[0]), &(fValues[0]), &(c[0]), ixst);
    fStatus = ixst;
    if (!fStatus) return false;
    fGraph.Set(t.size());
    for (size_t i = 0; i < t.size(); i++) {
      fGraph.SetPoint(i, t[i]/mass, fValues[i]);
      if (!fELoss) continue;
      fELoss->SetPoint(i, t[i]/mass, fValues[i]);
      // ~ 5 sigma
      fELoss->SetPointError(i, 0, .5 * fValues[i]);
    }
    fGraph.Write();
    if (fELoss) fELoss->Write();
    return true;
  }
  TGraph& Draw()
  {
    if (fELoss) fELoss->DrawClone("4 same");
    fGraph.DrawClone("l same");
    return fGraph;
  }
  const Mech&        fMech;
  std::vector<float> fValues;
  bool               fStatus;
  TGraph             fGraph;
  TGraphErrors*      fELoss;
};

//____________________________________________________________________
/** @ingroup FMD_xsec_script
 */
struct XSections 
{
  XSections(size_t n=91, float emin=1e-5, float emax=1e4) 
    : fTKine(n), 
      fCuts(5)
  {
    float dp   = 1. / log10(emax/emin);
    float pmin = log10(emin);
    fTKine[0]  = emin;
    for (size_t i = 1; i < fTKine.size(); i++) {
      float el  = pmin + i * dp;
      fTKine[i] = pow(10, el);
    }
    for (float_array::iterator i = fCuts.begin(); i != fCuts.end(); ++i) 
      *i = 1e-4;
    Mech* mech = &(fgMechs[0]);
    while (mech->fName) {
      fMechs.push_back(new MechValue(*mech, n));
      mech++;
    }
  }
  void Run(const std::string& medName, const std::string& pdgName) 
  {
    TGeant3* mc = static_cast<TGeant3*>(gMC);
    if (!mc) {
      std::cerr << "Couldn't get VMC" << std::endl;
      return;
    }
    TGeoMedium* medium = gGeoManager->GetMedium(medName.c_str());
    if (!medium) {
      std::cerr << "Couldn't find medium " << medName << std::endl;
      return;
    }
    int medNo = medium->GetMaterial()->GetUniqueID();
    TDatabasePDG* pdgDb = TDatabasePDG::Instance();
    fPDG                = pdgDb->GetParticle(pdgName.c_str());
    if (!fPDG) {
      std::cerr << "Couldn't find particle " << pdgName << std::endl;
      return;
    }
    int pdgNo = fPDG->PdgCode();
    int pidNo = mc->IdFromPDG(pdgNo);

    std::stringstream vars;
    vars << "betagamma/F";

    size_t nOk   = 0;
    // Loop over defined mechanisms 
    for (mech_array::iterator i = fMechs.begin(); i != fMechs.end(); ++i) {
      if (!(*i)->Get(mc, fTKine, fCuts, medNo, pidNo, fPDG->Mass()))continue;
      vars << ":" << (*i)->fMech.fName;
      nOk ++;
    }
    
    std::stringstream tName;
    tName << medName << "_" << pdgName;
    TTree* tree = new TTree(tName.str().c_str(), tName.str().c_str());

    float_array cache(nOk+1);
    tree->Branch("xsec", &(cache[0]), vars.str().c_str());
    for (size_t i = 0; i < fTKine.size(); i++) {
      cache[0] = fTKine[i] / fPDG->Mass();
      int k = 0;
      for (mech_array::iterator j = fMechs.begin(); j != fMechs.end(); ++j) {
	if (!(*j)->fStatus) continue;
	cache[++k] = (*j)->fValues[i];
      }
      tree->Fill();
    }
    tree->Write();
  }
  void Draw() 
  {
    float min = 100000;
    float max = 0;
    std::vector<TGraph*> gs;
    float_array bg(fTKine.size());
    for (size_t i = 0; i < fTKine.size(); i++) bg[i] = fTKine[i]/fPDG->Mass();
    for (mech_array::iterator j = fMechs.begin(); j != fMechs.end(); ++j) {
      if (!(*j)->fStatus) continue;
      for (size_t i = 0; i < fTKine.size(); i++) {
	if ((*j)->fValues[i] == 0) continue;
	min = std::min(min, (*j)->fValues[i]);
	max = std::max(max, (*j)->fValues[i]);
      }
    }

    TCanvas* c = new TCanvas("c", "C", 700, 700);
    c->SetFillColor(0);
    c->SetLogy();
    c->SetLogx();
    c->SetGridy();
    c->SetGridx();
    float_array y(101);
    float ymin = log10(min);
    float dy   = (log10(max)+.5 - log10(min)) / y.size();
    for (size_t i = 1; i < y.size(); i++)  y[i] = pow(10, ymin + i * dy);
    TH2* f = new TH2F("x", "X-sec",bg.size()-1,&(bg[0]),y.size()-1,&(y[0]));
    f->SetXTitle("#beta#gamma");
    f->SetDirectory(0);
    f->SetStats(kFALSE);
    f->Draw();
    TLegend* l = new TLegend(0.45, 0.125, 0.90, 0.45);
    l->SetFillColor(0);
    // l->SetFillStyle(0);
    for (mech_array::iterator j = fMechs.begin(); j != fMechs.end(); ++j) {
      if (!(*j)->fStatus) continue;
      TGraph& g = (*j)->Draw();
      l->AddEntry(&g, g.GetTitle(), "l");
    }
    l->Draw("same");
  }
protected: 
  typedef std::vector<float> float_array;
  typedef std::vector<MechValue*> mech_array;
  float_array   fTKine;
  float_array   fCuts;
  mech_array    fMechs;
  TParticlePDG* fPDG;
};

bool init = false;

void XSection(const char* pdgName="pi-")
{
  if (!init) {
    gAlice->InitMC("$(ALICE_ROOT)/FMD/Config.C");
    init = true;
  }
  TFile* file = TFile::Open("xsec.root", "RECREATE");
  XSections xs;
  xs.Run("FMD_Si$", pdgName);
  xs.Draw();
  file->Close();
}
//____________________________________________________________________
//
// EOF
//
