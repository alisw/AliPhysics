void 
TestTaskIO(bool read=false) 
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGLFforward.so");

  TFile* file = TFile::Open("task.root", (read ? "READ" : "RECREATE"));

  if (!read) TestTaskIOWrite(file);
  else       TestTaskIORead(file);

  file->Close();
}

void
TestTaskIOWrite(TFile* f)
{
  AliFMDAnaParameters* p = AliFMDAnaParameters::Instance();
  p->SetEnergy(AliFMDAnaParameters::k900);
  p->Init();
  p->Dump();

  f->cd();

  AliFMDAnalysisTaskSE* t = new AliFMDAnalysisTaskSE("FMD");
  t->Write();
  t->Print("p");
}


void
TestTaskIORead(TFile* f)
{
  AliFMDAnalysisTaskSE* t = static_cast<AliFMDAnalysisTaskSE*>(f->Get("FMD"));
  t->Print("p");

  AliFMDAnaParameters* p = AliFMDAnaParameters::Instance();
  p->Dump();
}
