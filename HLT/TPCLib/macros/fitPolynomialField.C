
int fitPolynomialField()
{
  gSystem->Load("libAliHLTTPC.so");
  AliHLTTPCGMPolynomialField polyField;
  AliMagF* fld = new AliMagF("Fit", "Fit", 1., 1., AliMagF::k2kG);
  AliHLTTPCGMPolynomialFieldCreator::FitField( fld,  polyField );
  return 0;
}
