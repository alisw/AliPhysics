AliGenerator *AddMCGenAmpt()
{
// User defined generator

  gSystem->Load("libampt.so");       
  gSystem->Load("libTAmpt.so");

  AliGenAmpt *genAMPT = new AliGenAmpt(-1);

  // will be made optional later
  genAMPT->SetEnergyCMS(2760);
  genAMPT->SetReferenceFrame("CMS");
  genAMPT->SetProjectile("A", 208, 82);
  genAMPT->SetTarget    ("A", 208, 82);
  genAMPT->SetPtHardMin (2);
  genAMPT->SetImpactParameterRange(0.00,20.00);
  genAMPT->SetJetQuenching(0); // enable jet quenching
  genAMPT->SetShadowing(1);    // enable shadowing
  genAMPT->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
  genAMPT->SetSpectators(0);   // track spectators 
  genAMPT->SetIsoft(4);        // 4=string melting, 1=standard AMPT
  genAMPT->SetXmu(3.2264);     // parton xsection
  genAMPT->SetNtMax(150);      // time bins
  
  genAMPT->SetAlpha(1./3.);    //alpha =0.333
  genAMPT->SetStringFrag(0.5,0.9); //string fragmentation parameters
  genAMPT->SetIpop(1); //enable popcorn mechanism (net-baryon stopping)
  // This particular choice of gives scattering cross section to be 1.5 mb

  genAMPT->SetRandomReactionPlane(kTRUE);

 return genAMPT;




}
