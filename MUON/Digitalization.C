void Digitalization()
{
Int_t nInputStreams = 1;
Int_t sperb=1;
gAlice->Delete();
AliRunDigitizer   * manager = new AliRunDigitizer(nInputStreams,sperb);
manager->SetInputStream(0,"galice.root");
manager->SetOutputFile("galice_digits.root");
AliMUONDigitizerv1* dMUON   = new AliMUONDigitizerv1(manager);
//dMUON->SetDebug(3);
manager->AddDigitizer(dMUON);
manager->Exec("deb");
}
