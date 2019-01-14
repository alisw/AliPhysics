
int fitPolynomialField()
{
  gSystem->Load("libAliHLTTPC.so");
  AliGPUTPCGMPolynomialField polyField;
  AliMagF* fld = new AliMagF("Fit", "Fit", 1., 1., AliMagF::k2kG);
  AliGPUTPCGMPolynomialFieldCreator::FitField( fld,  polyField );
  return 0;
}
