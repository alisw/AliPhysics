void Config() 
{
//
// A simple macro to test the AliFluka functionality
//

//
//  Fluka
//  
    new AliFluka("C++ Interface to FLUKA");
    AliFluka* fluka = (AliFluka*) gMC;
    printf("Pointer to Fluka Interface %p\n", gMC);
//
//  Output file
//
    TFile  *rootfile = new TFile("galice.root", "recreate");
    rootfile->SetCompressionLevel(2);
//
//  Generator
//
    AliGenHIJINGpara *gener = new AliGenHIJINGpara(50);
    gener->SetMomentumRange(0, 999);
    gener->SetPhiRange(0, 360);
    gener->SetThetaRange(45.,135.);
    gener->Init();

}
