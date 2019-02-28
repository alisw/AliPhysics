
int fitPolynomialField()
{
  gSystem->Load("libAliHLTTPC.so");
  GPUTPCGMPolynomialField polyField;
  AliMagF* fld = new AliMagF("Fit", "Fit", 1., 1., AliMagF::k2kG);
  GPUTPCGMPolynomialFieldCreator::FitField( fld,  polyField );
  return 0;
}
