Int_t eval(){

  AliL3Logger l;
  l.UseStdout();
  l.UseStream();

  Char_t *mcfile = "/heim/franken/data/V3.04/fast/hg_8k_v0_s1-3_e0.root";
//  Char_t *mcfile = "/heim/franken/data/V3.04/slow/hg_8k_v0_s1-3_e0.root";
  Char_t *mcClusterfile= "hg_8k_v0_s1-3_e0_cl.root";
  Char_t *trackfile ="tracks.raw";



  Int_t slice[2]={1,3};
  AliL3Evaluate *eval = new AliL3Evaluate(mcfile,slice);
  eval->SetupFast(trackfile,mcClusterfile);
//  eval->SetupSlow(trackfile);
//  eval->EvaluateSlice(1,10,10);
  eval->EvaluateGlobal(1,1);
  eval->Write2File("eval.root");
  
  return 0;
}
