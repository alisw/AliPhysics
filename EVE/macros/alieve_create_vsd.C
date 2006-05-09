void alieve_create_vsd()
{
  // Invoke as: aliroot alieve_create_vsd.C

  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libTreePlayer.so");
  gSystem->Load("libGed.so");

  gSystem->Load("libReve.so");
  gSystem->Load("libAlieve.so");

  Reve::DisablePODTObjectStreamers();

  Alieve::VSDCreator vc;
  vc.SetDebugLevel(2);
  vc.CreateVSD(".", 0, "AliVSD.root");
}
