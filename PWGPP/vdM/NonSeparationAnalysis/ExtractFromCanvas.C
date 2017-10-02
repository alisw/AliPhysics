// -*- C++ -*-


TString ExtractFromCanvas(TString fn) {
  TString line;
  TFile::Open(fn);
  TCanvas *c = (TCanvas*)canvas_par->GetListOfPrimitives()->At(0);
  TList *l = (TList*)c->GetListOfPrimitives();
  for (Int_t i=0; i<l->GetEntries()-1; ++i) {
    TString s = l->At(i)->GetTitle();
    //    Printf("%s", s.Data());
    s=s(s.Index("=")+1, 123456789);
    s.ReplaceAll("_", "");
    s.ReplaceAll("{", " ");
    s.ReplaceAll("^", " ");
    s.ReplaceAll("}", " ");
    s.ReplaceAll("#pm", " ");
    s.ReplaceAll("#mum", " ");
    s.ReplaceAll("#murad", " ");
    s.ReplaceAll("cm", " ");
    s.ReplaceAll("/", " ");
    s.ReplaceAll("=", " ");
    line += s.Data();
    line += " ";
  }
  Printf("%s", line.Data());
  return line;
}
