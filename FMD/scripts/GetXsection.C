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
#ifdef __CINT__
XSection()
{
  gROOT->ProcessLine(".x Compile.C(\"$(ALICE_ROOT)/FMD/scripts/XSection.C\"");
  gAlice->InitMC("$(ALICE_ROOT)/FMD/Config.C");
  TFile* file = TFile::Open("xsec.root", "RECREATE");
  GetXsection("FMD_Si$", "pi+");
  file->Close();
}
#else  
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
#include <TAxis.h>
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
};

Mech fgMechs[] = 
  {{ "HADF","total hadronic x-section according to FLUKA",      "cm^{1}"},
   { "INEF","hadronic inelastic x-section according to FLUKA",  "cm^{1}"},
   { "ELAF","hadronic elastic x-section according to FLUKA",    "cm^{1}"},
   { "HADG","total hadronic x-section according to GHEISHA",    "cm^{1}"},
   { "INEG","hadronic inelastic x-section according to GHEISHA","cm^{1}"},
   { "ELAG","hadronic elastic x-section according to GHEISHA",  "cm^{1}"},
   { "FISG","nuclear fission x-section according to GHEISHA",   "cm^{1}"},
   { "CAPG","neutron capture x-section according to GHEISHA",   "cm^{1}"},
   { "LOSS","stopping power",                                   "cm^{1}"},
   { "PHOT","photoelectric x-section",                          "MeV/cm"},
   { "ANNI","positron annihilation x-section",                  "cm^{1}"},
   { "COMP","Compton effect x-section",                         "cm^{1}"},
   { "BREM","bremsstrahlung x-section",                         "cm^{1}"},
   { "PAIR","photon and muon direct- pair x-section",           "cm^{1}"},
   { "DRAY","delta-rays x-section",                             "cm^{1}"},
   { "PFIS","photo-fission x-section",                          "cm^{1}"},
   { "RAYL","Rayleigh scattering x-section",                    "cm^{1}"},
   { "MUNU","muon-nuclear interaction x-section",               "cm^{1}"},
   { "RANG","range",                                            "cm"},
   { "STEP","maximum step",                                     "cm"},
   { 0, 0, 0,}};

//____________________________________________________________________
/** @ingroup FMD_xsec_script
 */
struct MechValue 
{
  MechValue(const Mech& m, size_t n) 
    : fMech(m), fValues(n), fStatus(0)
  {}
  bool Get(TGeant3* mc, std::vector<float>& t, 
	   std::vector<float>& c, int medNo, int pidNo) 
  {
    TGraph g(t.size());
    g.SetName(fMech.fName);
    g.SetTitle(fMech.fTitle);
    g.GetXaxis()->SetTitle("p [GeV]");
    g.GetYaxis()->SetTitle(fMech.fUnit);
    int ixst;
    mc->Gftmat(medNo, pidNo, fMech.fName, t.size(), 
	       &(t[0]), &(fValues[0]), &(c[0]), ixst);
    fStatus = ixst;
    if (!fStatus) return false;
    for (size_t i = 0; i < t.size(); i++) g.SetPoint(i, t[i], fValues[i]);
    g.Write();
    return true;
  }
  const Mech&        fMech;
  std::vector<float> fValues;
  bool               fStatus;
};

//____________________________________________________________________
/** @ingroup FMD_xsec_script
 */
struct XSection 
{
  XSection(size_t n=91, float emin=1e-5, float emax=1e4) 
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
    size_t i = 0;
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
    TParticlePDG* pdgP  = pdgDb->GetParticle(pdgName.c_str());
    if (!pdgP) {
      std::cerr << "Couldn't find particle " << pdgName << std::endl;
      return;
    }
    int pdgNo = pdgP->PdgCode();
    int pidNo = mc->IdFromPDG(pdgNo);

    std::stringstream vars;
    vars << "T/F";

    size_t nOk   = 0;
    // Loop over defined mechanisms 
    for (mech_array::iterator i = fMechs.begin(); i != fMechs.end(); ++i) {
      if (!(*i)->Get(mc, fTKine, fCuts, medNo, pidNo))continue;
      vars << ":" << (*i)->fMech.fName;
      nOk ++;
    }
    
    std::stringstream tName;
    tName << medName << "_" << pdgName;
    TTree* tree = new TTree(tName.str().c_str(), tName.str().c_str());

    float_array cache(nOk+1);
    tree->Branch("xsec", &(cache[0]), vars.str().c_str());
    for (size_t i = 0; i < fTKine.size(); i++) {
      cache[0] = fTKine[i];
      int k = 0;
      for (mech_array::iterator j = fMechs.begin(); j != fMechs.end(); ++j) {
	if (!(*j)->fStatus) continue;
	cache[++k] = (*j)->fValues[i];
      }
      tree->Fill();
    }
    tree->Write();
  }
protected: 
  typedef std::vector<float> float_array;
  typedef std::vector<MechValue*> mech_array;
  float_array  fTKine;
  float_array  fCuts;
  mech_array   fMechs;
};

#if 0
//____________________________________________________________________
/** @ingroup FMD_xsec_script
    @param medName 
    @param pdgName 
    @param n 
    @param emin 
    @param emax 
*/
void
GetXsection(const char* medName, const char* pdgName,
	    Int_t n=91, Float_t emin=1e-5, Float_t emax=1e4)
{
  TArrayF tkine(n);
  Float_t dp   = 1/TMath::Log10(emax/emin);
  Float_t pmin = TMath::Log10(emin);
  tkine[0]     = emin;
  for (Int_t i=1; i < tkine.fN; i++) {
    Float_t el = pmin + i * dp;
    tkine[i]   = TMath::Power(10, el);
  }
  TArrayF cuts(5);
  cuts.Reset(1e-4);

  Mech mechs[] = 
    {{ "HADF","total hadronic x-section according to FLUKA","cm^{1}",n,0},
     { "INEF","hadronic inelastic x-section according to FLUKA","cm^{1}",n,0},
     { "ELAF","hadronic elastic x-section according to FLUKA","cm^{1}",n,0},
     { "HADG","total hadronic x-section according to GHEISHA","cm^{1}",n,0},
     { "INEG","hadronic inelastic x-section according to GHEISHA","cm^{1}",n,0},
     { "ELAG","hadronic elastic x-section according to GHEISHA","cm^{1}",n,0},
     { "FISG","nuclear fission x-section according to GHEISHA","cm^{1}",n,0},
     { "CAPG","neutron capture x-section according to GHEISHA","cm^{1}",n,0},
     { "LOSS","stopping power","cm^{1}",n,0},
     { "PHOT","photoelectric x-section","MeV/cm",n,0},
     { "ANNI","positron annihilation x-section","cm^{1}",n,0},
     { "COMP","Compton effect x-section","cm^{1}",n,0},
     { "BREM","bremsstrahlung x-section","cm^{1}",n,0},
     { "PAIR","photon and muon direct- pair x-section","cm^{1}",n,0},
     { "DRAY","delta-rays x-section","cm^{1}",n,0},
     { "PFIS","photo-fission x-section","cm^{1}",n,0},
     { "RAYL","Rayleigh scattering x-section","cm^{1}",n,0},
     { "MUNU","muon-nuclear interaction x-section","cm^{1}",n,0},
     { "RANG","range","cm",n,0},
     { "STEP","maximum step","cm",n,0},
     { 0, 0, 0, 0, 0}};
  TGeant3* mc = (TGeant3*)gMC;
  if (!mc) {
    std::cerr << "Couldn't get VMC" << std::endl;
    return;
  }
  TGeoMedium* medium = gGeoManager->GetMedium(medName);
  if (!medium) {
    std::cerr << "Couldn't find medium " << medName << std::endl;
    return;
  }
  Int_t medNo = medium->GetMaterial()->GetUniqueID();
  TDatabasePDG* pdgDb = TDatabasePDG::Instance();
  TParticlePDG* pdgP  = pdgDb->GetParticle(pdgName);
  if (!pdgP) {
    std::cerr << "Couldn't find particle " << pdgName << std::endl;
    return;
  }
  Int_t pdgNo = pdgP->PdgCode();
  Int_t pidNo = mc->IdFromPDG(pdgNo);
    
  Mech* mech = &(mechs[0]);
  Int_t nMech = 0;
  Int_t nOk = 0;
  TString vars("T/F");
  while (mech->name) {
    cout << mech->name << ": " << mech->title << " ... " << std::flush;
    nMech++;
    Int_t ixst;
    mc->Gftmat(medNo, pidNo, mech->name, n, 
	       tkine.fArray, mech->values.fArray, cuts.fArray, ixst);
    mech->status = ixst;
    if (ixst) {
      nOk++;
      vars.Append(Form(":%s", mech->name));
      if (!strcmp("LOSS", mech->name)) {
	for (Int_t i = 0; i < n; i++) 
	  std::cout << i << "\t" << tkine[i] << "\t" 
		    << mech->values[i] << std::endl;
      }
    }
    std::cout << (ixst ? "ok" : "failed") << std::endl;
    mech++;
  }
  // TFile* file = TFile::Open(Form("xsec-%d.root", pdgNo),
  // "RECREATE");
  TArrayF cache(nOk+1);
  TTree* tree = new TTree(Form("%s_%s", medName, pdgName), 
			  Form("%s_%s", medName, pdgName));
  tree->Branch("xsec", cache.fArray, vars.Data());
  for (Int_t i = 0; i < n; i++) {
    cache[0] = tkine[i];
    Int_t k = 0;
    for (Int_t j = 0; j < nMech; j++) {
      if (mechs[j].status) {
	if (!strcmp(mechs[j].name, "LOSS")) 
	  std::cout << tkine[i] << "\t" << mechs[j].values[i] << std::endl;
	cache[k+1] = mechs[j].values[i];
	k++;
      }
    }
    std::cout << k << "\t" << (k == nOk) << std::endl;
    tree->Fill();
  }
  tree->Write();
}
#endif
#endif
//____________________________________________________________________
//
// EOF
//
