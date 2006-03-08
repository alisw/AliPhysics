#ifndef __CINT__
#include <AliESD.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TChain.h>
#include <iostream>
#endif

void
CompareESD()
{
  TChain*           chain = new TChain("esdTree");
  TSystemDirectory* dir   = new TSystemDirectory(".", ".");
  TList*            files = dir->GetListOfFiles();
  files->Sort();
  TIter next(files);
  TSystemFile*      file  = 0;
  TObjArray*        other = new TObjArray;
  while ((file = static_cast<TSystemFile*>(next()))) {
    TString fname(file->GetName());
    if (fname.Contains("AliESDs")) 
      chain->AddFile(fname.Data());
    if (fname.Contains("FMD.ESD.")) {
      TFile* o = TFile::Open(fname.Data());
      if (!o) {
	Warning("CompareESD", "Failed to open file %s", fname.Data());
	other->Add(0);
	continue;
      }
      AliESDFMD* fmdEsd = static_cast<AliESDFMD*>(o->Get("AliESDFMD"));
      if (!fmdEsd) {
	Warning("CompareESD", "Failed to get FMD object from reference %s", 
		fname.Data());
	other->Add(0);
	continue;
      }
      other->Add(fmdEsd);
    }
  }
  delete dir;
  
  AliESD* esd = 0;
  chain->SetBranchAddress("ESD", &esd);
  Int_t   n   = chain->GetEntries();
  for (Int_t i = 0; i < n; i++) {
    Int_t read = chain->GetEntry(i);
    if (read <= 0) break;
    std::cout << "Entry # " << i << std::endl;
    if (!esd) break;
    esd->Print();
    AliESDFMD* esdObj = esd->GetFMDData();
    AliESDFMD* refObj = static_cast<AliESDFMD*>(other->At(i));
    if (!esdObj) {
      Warning("CompareESD", "no FMD object in ESD");
      continue;
    }
    if (!refObj) {
      Warning("CompareESD", "no reference FMD object");
      continue;
    }
    std::cout << " Comparing ... " << std::endl;
    for (UShort_t det = 1; det <= esdObj->MaxDetectors(); det++) {
      std::cout << "FMD" << det << std::endl;
      for (UShort_t ir = 0; ir < esdObj->MaxRings(); ir++) {
	Char_t ring = (ir == 0 ? 'I' : 'O');
	std::cout << " Ring " << ring << std::endl;
	for (UShort_t sec = 0; sec < esdObj->MaxSectors(); sec++) {
	  std::cout << "  Sector # " << sec << std::endl;
	  for (UShort_t str = 0; str < esdObj->MaxStrips(); str++) {
	    if (esdObj->Multiplicity(det, ring, sec, str) != 
		refObj->Multiplicity(det, ring, sec, str))
	      Warning("CompareESD", 
		      "Mult for FMD%d%c[%2d,%3d]",det,ring,sec,str);
	    if (esdObj->Eta(det, ring, sec, str) != 
		refObj->Eta(det, ring, sec, str))
	      Warning("CompareESD", 
		      "Eta for FMD%d%c[%2d,%3d]", det,ring,sec,str);
	    if (esdObj->Multiplicity(det, ring, sec, str) > 0) 
	      Info("CompareESD", "Mult in FMD%c%d[%2d,%3d] is %f",
		   det,ring,sec,str,
		   esdObj->Multiplicity(det, ring, sec, str));
	  }
	}
      }
    }
  }
}

      
  
