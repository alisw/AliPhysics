void AliTRDdigits2raw()
{

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(2);
  raw->OpenInput("galice.root");
  raw->Digit2Raw();

}
