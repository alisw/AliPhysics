void alifasttest()
{

  AliFast * fast =new AliFast("blabla","blabla");
  gAliFast->Init("Config.C");

  TFile file("alifasttest.root","recreate","ALIFast root file",2);
  gAliFast.Init();
  gAliFast.SetTestTrack(1);
  gAliFast.MakeTree();
  gAliFast->Generator()->Generate();

  gAliFast.Make();       // Generate and Reconstruct event
  gAliFast.FillTree();    
  gAliFast.Clear();       // Clear reconstructed event lists

  gAliFast.Finish();
  gAliFast.Write(); 
}
