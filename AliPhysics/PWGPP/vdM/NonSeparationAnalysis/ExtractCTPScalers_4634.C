// -*- C++ -*-

TString classNames[] = {
  "CINT7-B-NOPF-FASTNOTRD",
  "CINT7-I1-NOPF-FASTNOTRD",
  "CINT7-I2-NOPF-FASTNOTRD",
  "CINT7-I3-NOPF-FASTNOTRD",
  "CINT7-I4-NOPF-FASTNOTRD",
  "CINT7-I5-NOPF-FASTNOTRD",
  "CINT7-I6-NOPF-FASTNOTRD",
  "CINT7-I7-NOPF-FASTNOTRD",
  "CINT7-I8-NOPF-FASTNOTRD",
  "C0TVX-B-NOPF-FASTNOTRD",
  "C0TVX-I1-NOPF-FASTNOTRD",
  "C0TVX-I2-NOPF-FASTNOTRD",
  "C0TVX-I3-NOPF-FASTNOTRD",
  "C0TVX-I4-NOPF-FASTNOTRD",
  "C0TVX-I5-NOPF-FASTNOTRD",
  "C0TVX-I6-NOPF-FASTNOTRD",
  "C0TVX-I7-NOPF-FASTNOTRD",
  "C0TVX-I8-NOPF-FASTNOTRD",
  "CADAND-B-NOPF-FASTNOTRD",
  "CADAND-I1-NOPF-FASTNOTRD",
  "CADAND-I2-NOPF-FASTNOTRD",
  "CADAND-I3-NOPF-FASTNOTRD",
  "CADAND-I4-NOPF-FASTNOTRD",
  "CADAND-I5-NOPF-FASTNOTRD",
  "CADAND-I6-NOPF-FASTNOTRD",
  "CADAND-I7-NOPF-FASTNOTRD",
  "CADAND-I8-NOPF-FASTNOTRD",
  "CVBANOTC-I1-NOPF-FASTNOTRD",
  "CVBANOTC-I2-NOPF-FASTNOTRD",
  "CVBANOTC-I3-NOPF-FASTNOTRD",
  "CVBANOTC-I4-NOPF-FASTNOTRD",
  "CVBANOTC-I5-NOPF-FASTNOTRD",
  "CVBANOTC-I6-NOPF-FASTNOTRD",
  "CVBANOTC-I7-NOPF-FASTNOTRD",
  "CVBANOTC-I8-NOPF-FASTNOTRD",
  "CVBCNOTA-I1-NOPF-FASTNOTRD",
  "CVBCNOTA-I2-NOPF-FASTNOTRD",
  "CVBCNOTA-I3-NOPF-FASTNOTRD",
  "CVBCNOTA-I4-NOPF-FASTNOTRD",
  "CVBCNOTA-I5-NOPF-FASTNOTRD",
  "CVBCNOTA-I6-NOPF-FASTNOTRD",
  "CVBCNOTA-I7-NOPF-FASTNOTRD",
  "CVBCNOTA-I8-NOPF-FASTNOTRD",
  "CADANOTC-I1-NOPF-FASTNOTRD",
  "CADANOTC-I2-NOPF-FASTNOTRD",
  "CADANOTC-I3-NOPF-FASTNOTRD",
  "CADANOTC-I4-NOPF-FASTNOTRD",
  "CADANOTC-I5-NOPF-FASTNOTRD",
  "CADANOTC-I6-NOPF-FASTNOTRD",
  "CADANOTC-I7-NOPF-FASTNOTRD",
  "CADANOTC-I8-NOPF-FASTNOTRD",
  "CADCNOTA-I1-NOPF-FASTNOTRD",
  "CADCNOTA-I2-NOPF-FASTNOTRD",
  "CADCNOTA-I3-NOPF-FASTNOTRD",
  "CADCNOTA-I4-NOPF-FASTNOTRD",
  "CADCNOTA-I5-NOPF-FASTNOTRD",
  "CADCNOTA-I6-NOPF-FASTNOTRD",
  "CADCNOTA-I7-NOPF-FASTNOTRD",
  "CADCNOTA-I8-NOPF-FASTNOTRD",
  "CT0ANOTC-I1-NOPF-FASTNOTRD",
  "CT0ANOTC-I2-NOPF-FASTNOTRD",
  "CT0ANOTC-I3-NOPF-FASTNOTRD",
  "CT0ANOTC-I4-NOPF-FASTNOTRD",
  "CT0ANOTC-I5-NOPF-FASTNOTRD",
  "CT0ANOTC-I6-NOPF-FASTNOTRD",
  "CT0ANOTC-I7-NOPF-FASTNOTRD",
  "CT0ANOTC-I8-NOPF-FASTNOTRD",
  "CT0CNOTA-I1-NOPF-FASTNOTRD",
  "CT0CNOTA-I2-NOPF-FASTNOTRD",
  "CT0CNOTA-I3-NOPF-FASTNOTRD",
  "CT0CNOTA-I4-NOPF-FASTNOTRD",
  "CT0CNOTA-I5-NOPF-FASTNOTRD",
  "CT0CNOTA-I6-NOPF-FASTNOTRD",
  "CT0CNOTA-I7-NOPF-FASTNOTRD",
  "CT0CNOTA-I8-NOPF-FASTNOTRD"
};

const Int_t nClasses = sizeof(classNames)/sizeof(TString);

TTree *TT = new TTree;
void ExtractCTPScalers_4634()
{
  gROOT->LoadMacro("AnalyzeCTPRecords.C+");
  gROOT->LoadMacro("ExtractRateFromCTPScalers.C+");

  // AnalyzeCTPRecords(244369, 2015, classNames, nClasses);
  // AnalyzeCTPRecords(244375, 2015, classNames, nClasses);

  const char *suffix[5] = { "", "_PU1", "_PU2", "_BB1", "_BB2" };
  const Int_t     bb[5] = { 0,    0,      0,      -1,    1     };
  const Double_t  pu[5] = { 0.0,  0.9,    1.0,   0.0,    0.0   };

  const TString classIDs[] = { "CADAND" };
  const TString bunchIDs[] = { "B", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8" };
  const Int_t   nBunches[] = {  22,    1,    1,    1,    1,    1,    1,    1,    1 };

  const Double_t ratioA = 0.129; // AnotC/AandC (AD)
  const Double_t ratioC = 0.287; // CnotA/AandC (AD)

  for (Int_t l=0; l<5; ++l) {
    for (Int_t i=0; i<sizeof(classIDs)/sizeof(TString); ++i) {
      for (Int_t j=0; j<sizeof(bunchIDs)/sizeof(TString); ++j) {
	ExtractRateFromCTPScalers(classIDs[i], bunchIDs[j], 244369, "txt/4634/ScanXY.txt",     nBunches[j], 4634, pu[l]*ratioA, pu[l]*ratioC,
				  "root/4634/IntensityDecay.root", "root/4634/SatelliteContaminationAD.root", bb[l], TT);
	ExtractRateFromCTPScalers(classIDs[i], bunchIDs[j], 244375, "txt/4634/ScanOffset.txt", nBunches[j], 4634, pu[l]*ratioA, pu[l]*ratioC,
				  "", "root/4634/SatelliteContaminationAD.root", bb[l], TT);
      }
    }
    TFile *f = TFile::Open(Form("root/4634/RatesAndSep_4634%s.root", suffix[l]), "RECREATE");
    TT->Write("TT");
    f->Write();
    f->Close();
  }
}
