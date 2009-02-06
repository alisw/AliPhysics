void
FillEtas(AliFMDFloatMap& etas)
{
  AliCDBManager*  cdb  = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  AliGeomManager::LoadGeometry("geometry.root");
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  
  etas.Clear();
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRing = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRing; q++) { 
      Char_t   r    = (q == 0 ? 'I' : 'O');
      UShort_t nStr = (q == 0 ? 512 : 256);
      for (UShort_t t = 0; t < nStr; t++) { 
	Double_t x, y, z;
	geom->Detector2XYZ(d, r, 0, t, x, y, z);
	Double_t phi   =  TMath::ATan2(y, x);
	Double_t rr    =  TMath::Sqrt(y * y + x * x);
	Double_t theta =  TMath::ATan2(rr, z);
	Double_t eta   = -TMath::Log(TMath::Tan(theta / 2));
	etas(d, r, 0, t) =  eta;
      }
    }
  }
}

Double_t
MakeMult()
{
  return TMath::Min(20., gRandom->Landau(1, .1));
}

void
FillESDFMD(AliESDFMD& fmd, AliFMDFloatMap& etas, Double_t ratio)
{
  fmd.Clear();
  
  Int_t nTotal  = 51200;
  Int_t nFilled = 0;

  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRing = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRing; q++) { 
      Char_t   r    = (q == 0 ? 'I' : 'O');
      UShort_t nSec = (q == 0 ?  20 :  40);
      UShort_t nStr = (q == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nSec; s++) {
	for (UShort_t t = 0; t < nStr; t++) { 
	  if (s == 0) { 
	    Double_t eta = etas(d, r, s, t);
	    fmd.SetEta(d, r, s, t, eta);
	  }
	  Double_t now  = Double_t(nFilled) / nTotal;
	  Bool_t   fill = now < ratio;
	  fmd.SetMultiplicity(d, r, s, t, (!fill ? 0 : MakeMult()));
	  if (fill) nFilled++;
	}
      }
    }
  }
}

Float_t
MakeTree(Int_t ratio, Int_t split, AliFMDFloatMap& etas)
{
  // Split is 
  // 
  //  0 - No split of members
  //  1 - Split on data and sub-objects.
  //  2 - Same as 1
  // 99 - Default - same as 2
  
  AliESDFMD* fmd = new AliESDFMD;
  TString    fname(Form("test_%03d_%02d.root", ratio, split));
  TString    tname(Form("T_%03d_%02d", ratio, split));
  TFile*     file = TFile::Open(fname.Data(),"RECREATE");
  TTree*     tree = new TTree(tname.Data(), 
			      Form("Ratio=%3d%%, Split=%2d", ratio, split));
  TBranch* branch = tree->Branch(split == 0 ? "FMDESD" : "FMDESD.", 
				 &fmd, 32000, split);
  // branch->SetCompressionLevel(9);
  
  std::cout << "  Doing " << tree->GetTitle() << " " << std::flush;
  Int_t nEvents = 10;
  for (Int_t i = 0; i < nEvents; i++) { 
    std::cout << "." << std::flush;
    FillESDFMD(*fmd, etas, Double_t(ratio) / 100);
    tree->Fill();
  }
  
  
  // branch->Print();
  tree->GetCurrentFile()->cd();
  tree->Write();
  tree->GetCurrentFile()->Close();
  
  file = TFile::Open(fname.Data(), "UPDATE");
  tree = static_cast<TTree*>(file->Get(tname.Data()));
  if (split == 0) 
    branch = tree->FindBranch("FMDESD");
  else if (split == 1) { 
    TLeaf* leaf = tree->FindLeaf("FMDESD.fMultiplicity");
    if (leaf) branch = leaf->GetBranch();
    else Warning("FillTree", "Failed to find FMDESD.fMultiplicity");
  }
  else if (split > 1) { 
    TLeaf* leaf = tree->FindLeaf("FMDESD.fMultiplicity.fData");
    if (leaf) branch = leaf->GetBranch();
    else Warning("FillTree", "Failed to find FMDESD.fMultiplicity.fData");
  }
  
  TH1* h = new TH1F("h", tree->GetTitle(), 100, 0, 20);
  tree->Draw("FMDESD.fMultiplicity.fData>>h", 
	     "FMDESD.fMultiplicity.fData<20&&FMDESD.fMultiplicity.fData>0");
  Int_t   nEntries = h->GetEntries();
  Float_t fratio   = Float_t(nEntries) / 51200 / nEvents;
  
  Int_t    comp    = branch->GetCompressionLevel();
  Long64_t totSize = branch->GetTotBytes();
  Long64_t zipSize = branch->GetZipBytes();
  Float_t  compz   = (zipSize ? Float_t(totSize+0.00001)/zipSize : 1);
  std::cout << " done - " 
	    << "Compression: " << compz << " (" << comp << ", " 
	    << Int_t(fratio*1000+.5)/10. << "%)" << std::endl;
  // tree->GetCurrentFile()->Write();
  // tree->GetCurrentFile()->Close();

  // delete tree;
  // delete fmd;
  
  return compz;
}

void
ESDDataSize()
{
  gStyle->SetPalette(1);
  
  AliFMDFloatMap etas(3, 2, 1, 512);
  FillEtas(etas);
  
  Int_t  ratios[] = { 1, 10, 25, 50, 100, -1 };
  Int_t  splits[] = { 0, 1, 99, -1 };
  Int_t* pratio   = ratios;

  Double_t arat[] = { 0.5, 1.5, 9.5, 10.5, 24.5, 25.5, 49.5, 50.5, 99.5, 100.5};
  Double_t aspl[] = { -0.5, 0.5, 1.5, 98.5, 99.5 }; 
  TH2* h = new TH2F("size", "ESD size", 9, arat, 4, aspl);
  h->GetXaxis()->SetTitle("Fill ratio");
  h->GetYaxis()->SetTitle("Split level");
  h->SetDirectory(0);
  h->Draw("COLZ");
  std::cout << "Loop over ratios" << std::endl;
  while (*pratio > 0) { 
    Int_t  ratio  = *pratio++;
    Int_t* psplit = splits;

    std::cout << " Loop over splits" << std::endl;
    while (*psplit >= 0) { 
      Int_t split = *psplit++;
      
      Double_t comp = MakeTree(ratio, split, etas);
      h->Fill(ratio, split, comp);
      // h->Draw("COLZ");
    }
  }
  h->Draw("LEGO2");
}

  
