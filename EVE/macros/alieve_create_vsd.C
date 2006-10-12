void alieve_create_vsd()
{
  // Invoke as: aliroot alieve_create_vsd.C

  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libTreePlayer");
  gSystem->Load("libGed");
  gSystem->Load("libRGL");

  gSystem->Load("libReve");
  gSystem->Load("libAlieve");

  Reve::DisablePODTObjectStreamers();

  Alieve::VSDCreator vc;
  vc.SetDebugLevel(2);
  vc.CreateVSD(".", 0, "AliVSD.root");
}
