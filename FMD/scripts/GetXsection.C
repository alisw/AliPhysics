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
/** @defgroup xsec_script X-section script 
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** @ingroup xsec_script
 */
struct Mech 
{
  char*    name;
  char*    title;
  char*    unit;
  TArrayF  values;
  int      status;
};


//____________________________________________________________________
/** @ingroup xsec_script
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
//____________________________________________________________________
//
// EOF
//
