void MakeVZEROEqualizationFactorsFor2010pp(const char *cdburi = "local://$ALICE_ROOT/OCDB")
{
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdburi);

  TFile *finput = TFile::Open("../PWGCF/EBYE/BalanceFunctions/VZEROEqualization/pp/CalibrationVZERO.root");
  const char *periods[4] = {"b","c","d","e"};

  for(Int_t i = 0; i < 4; ++i) {
    TList *list = (TList*)finput->Get(Form("LHC10%s",periods[i]));
    TH2F *histo = (TH2F*)list->FindObject("gHistVZEROChannelGainEqualizationMap");
    for(Int_t j = 1; j <= histo->GetNbinsX(); ++j) {
      TString str = histo->GetXaxis()->GetBinLabel(j);
      Int_t runN = str.Atoi();
      Double_t sum = 0.;
      for(Int_t k = 0; k < 64; ++k) {
	sum += histo->GetBinContent(j,k+1);
      }
      if (sum < 1e-6) continue;
      Double_t factors[66];
      factors[0] = factors[65] = 0.;
      for(Int_t k = 0; k < 64; ++k) {
	factors[k+1] = histo->GetBinContent(j,k+1)/sum*64.;
      }

      TH1F *eqFactors = new TH1F("VZEROEqualizationFactors","VZERO Equalization Factors for 2010 p-p data",64,-0.5,63.5);
      eqFactors->SetContent(factors);

      AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
      md->SetResponsible("Brigitte Cheynis");
      md->SetBeamPeriod(0);
      md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
      md->SetComment("Entry for VZERO Equalization Factors object (2010 p-p)");
      md->PrintMetaData();

      AliCDBStorage *storLoc = man->GetDefaultStorage();
      AliCDBId id("VZERO/Calib/EqualizationFactors",runN,runN);

      storLoc->Put(eqFactors, id, md);

      delete eqFactors;
      //      storLoc->Delete();
      //     delete md;
    }
  }
}
