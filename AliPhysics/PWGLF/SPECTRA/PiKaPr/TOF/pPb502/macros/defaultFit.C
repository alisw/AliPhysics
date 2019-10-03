defaultFit(const Char_t *filename, Int_t ipart, Int_t icharge, Int_t iipart, Int_t icent)
{
  gROOT->LoadMacro("TOFpid.C");
  TOFspectrum(filename, ipart, icharge, iipart, icent);
  TOFpid_rawSpectra();
  TOFpid_rawSpectra_mismatchCorrected();
}
