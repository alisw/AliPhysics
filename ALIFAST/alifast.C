void alifast()
{

  AliFast * fast =new AliFast("blabla","blabla");
  gAliFast->Init("Config.C");

  TFile file("alifast.root","recreate","ALIFast root file",2);
  gAliFast.Init();
  gAliFast.MakeTree();
  gAliFast->Generator()->Generate();

  gAliFast.Make(0);       // Generate and Reconstruct event
  gAliFast.FillTree();    
  gAliFast.Clear();       // Clear reconstructed event lists

  gAliFast.Finish();
  gAliFast.Write(); 
}
