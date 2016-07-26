#include <cstdarg>
#include <cstring>
#include <cstdio>

void MyPrint(const char* msg)
{
  gROOT->IndentLevel();
  Printf(msg);
}
void MyIncr()
{
  gROOT->IncreaseDirLevel();
}
void MyDecr()
{
  gROOT->DecreaseDirLevel();
}
struct MyGuard
{
  TString fMsg;
  MyGuard(const char* msg)
    : fMsg(msg)
  {
    MyPrint(msg);
    MyIncr();
  }
  ~MyGuard() { MyDecr(); MyPrint(Form("End %s", fMsg.Data())); }
};
  


void AddPath(const TString& dir, Bool_t pre=true)
{
  TString mPath(gROOT->GetMacroPath());
  if (pre) {
    mPath.Prepend(":");
    mPath.Prepend(dir);
  }
  else {
    mPath.Append(":");
    mPath.Append(dir);
  }
  gROOT->SetMacroPath(mPath);  
  gSystem->AddIncludePath(Form("-I%s", dir.Data()));
}


void ProcessCandidates(const TSeqCollection* inp, TList& toDraw)
{
  TIter next(inp);
  TObject* o = 0;
  while ((o = next())) {
    GraphSysErr* g = static_cast<GraphSysErr*>(o);
    TString obsKey = g->GetKey("obskey");
    TString det    = g->GetKey("detector");
    MyPrint(Form("ObsKey=%s Detector=%s", obsKey.Data(), det.Data()));
    // g->Print("key");
    TString opt("");
    if (det.Contains("TRACKLET", TString::kIgnoreCase)) {
      opt="XBASE=-4";
      Color_t col = g->GetMarkerColor();
      g->SetDataOption(GraphSysErr::kNoTick);
      g->SetSumOption(GraphSysErr::kBox);
      g->SetSumFillStyle(0);
      g->SetSumFillColor(col);
      g->SetSumLineColor(col);
      g->SetCommonSumOption(GraphSysErr::kRect);
      g->SetCommonSumFillStyle(1001);
      g->SetCommonSumFillColor(col);
      g->SetCommonSumLineColor(col);
    }
    if (obsKey.EqualTo("DN/DETARAP", TString::kIgnoreCase))
      toDraw.Add(g, opt);
    
  }
}

void ProcessHEPData(const TString& filename, TList& toDraw)
{
  MyGuard g(Form("HEPData %s", filename.Data()));
  TSeqCollection* list = GraphSysErr::Import(filename);
  ProcessCandidates(list,toDraw);
}

void ProcessFile(const TString& filename, TList& toDraw)
{
  MyGuard g(Form("File %s", filename.Data()));
  TFile* file = TFile::Open(filename, "READ");
  if (!file) return;

  TList candidates;
  ProcessFileDirectory(file, candidates);
  ProcessCandidates(&candidates, toDraw);
  file->Close();
}

void ProcessCollection(const TCollection* d, TList& toInspect)
{
  MyGuard g(Form("Collection: %s", d->GetName()));
  TIter    next(d);
  TObject* obj = 0;
  while ((obj = next())) {
    TClass* cls = obj->IsA();
    if (!cls) {
      MyPrint(Form("W: Object %s in %s has unknown class: %s",
		   obj->GetName(), d->GetName(), obj->ClassName()));
      continue;
    }
    if (cls->InheritsFrom(TDirectory::Class())) {
      ProcessFileDirectory(static_cast<TDirectory*>(obj),toInspect);
      continue;
    }
    if (cls->InheritsFrom(TCollection::Class())) {
      ProcessCollection(static_cast<TCollection*>(obj), toInspect);
      continue;
    }       
    if (cls->InheritsFrom(GraphSysErr::Class())) {
      MyPrint(Form("Adding %s/%s to candidate list",
		   d->GetName(), obj->GetName()));
      toInspect.Add(obj);
    }
  }
}

void ProcessFileDirectory(const TDirectory* d, TList& toInspect)
{
  MyGuard g(Form("File directory: %s", d->GetPath()));
  TIter next(d->GetListOfKeys());
  TKey* key = 0;
  while ((key = static_cast<TKey*>(next()))) {
    TClass* cls = gROOT->GetClass(key->GetClassName());
    if (!cls) {
      MyPrint(Form("W: Object %s in %s has unknown class: %s",
		   key->GetName(), d->GetPath(), key->GetClassName()));
      continue;
    }
    if (cls->InheritsFrom(TDirectory::Class())) {
      ProcessFileDirectory(static_cast<TDirectory*>(key->ReadObj()),toInspect);
      continue;
    }
    if (cls->InheritsFrom(TCollection::Class())) {
      ProcessCollection(static_cast<TCollection*>(key->ReadObj()), toInspect);
      continue;
    }       
    if (cls->InheritsFrom(GraphSysErr::Class())) {
      MyPrint(Form("Adding %s/%s to candidate list",
		   d->GetPath(), key->GetName()));
      toInspect.Add(key->ReadObj());
    }
  }
}

void ProcessOne(const TString& filename, TList& toDraw)
{
  if (filename.EndsWith(".input"))
    ProcessHEPData(filename, toDraw);
  else
    ProcessFile(filename, toDraw);
}

void QuickDraw()
{
  QuickDraw(&(gApplication->Argv()[1]));
}
void QuickDraw(const char** args)
{
  if (!gROOT->GetClass("GraphSysErr")) {
    AddPath("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/gse");
    AddPath("$HOME/GraphSysErr");
    gROOT->LoadMacro("GraphSysErr.C+g");
  }    
  TList list;
  {
    // Int_t  argc = gApplication->Argc();
    // gApplication->ClearInputFiles();
    MyGuard* g1 = new MyGuard("Processing arguments");
    const char** ptr = args;
    while (*ptr) {
      TString argi = *ptr;
      ptr++;
      if (argi.BeginsWith("-")) continue;
      if (argi.EndsWith("root.exe")) continue;
      if (argi.EndsWith("QuickDraw.C")) continue;
      if (!argi.EndsWith(".input") && !argi.EndsWith(".root")) continue;
      gSystem->ExpandPathName(argi);
      MyPrint(Form("Adding file %s", argi.Data()));
      list.Add(new TObjString(argi));
      // gApplication->Argv(i) = 0;
    }
    delete g1;
  }
  

  TList toDraw;
  {
    MyGuard* g2 = new MyGuard("Processing found files");
    TIter next(&list);
    TObject* obj = 0;
    while ((obj = next())) {
      ProcessOne(obj->GetName(), toDraw);
    }
    delete g2;
  }

  // toDraw.ls();
  TString opt("COMBINED QUAD");
  Double_t     ymin = 1e9, ymax = -1e9;
  Double_t     xmin = 1e9, xmax = -1e9;
  {
    MyGuard* g3 = new MyGuard("Figuring out min/max");
    TIter        next(&toDraw);
    GraphSysErr* g = 0;
    while ((g = static_cast<GraphSysErr*>(next()))) {
      Double_t mn, mx;
      g->GetMinMax(opt, mn, mx);
      ymax        = TMath::Max(mx,ymax);
      ymin        = TMath::Min(mn,ymin);
      Int_t    nx = g->GetN()-1;
      Double_t x1 = g->GetX(0) -2*g->GetErrorXLeft(0);
      Double_t x2 = g->GetX(nx)+2*g->GetErrorXRight(nx);
      xmin        = TMath::Min(x1,xmin);
      xmax        = TMath::Max(x2,xmax);

    }
    MyPrint(Form("x-range: [%f,%f]   y-range: [%f,%f]",xmin,xmax,ymin,ymax));
    delete g3;
  }

  
  TCanvas* c = new TCanvas("c","c", 1000, 1000);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->SetLeftMargin(ymax > 100 ? 0.14 : 0.1);
  TH1*     f = new TH2D("frame","",300,
			(xmin < 0 ? 1.3 : 0.7)*xmin,
			(xmax < 0 ? 0.7 : 1.3)*xmax,
			300, 0.9*ymin, 1.1*ymax);
  f->SetStats(0);
  f->SetXTitle("\\eta");
  f->SetYTitle("\\hbox{d}N_{\\hbox{ch}}/\\hbox{d}\\eta");
  f->GetYaxis()->SetTitleOffset(ymax > 1000 ? 1.7 : 1.1);
  f->GetYaxis()->SetNdivisions(210);
  f->GetXaxis()->SetNdivisions(210);
  f->Draw();
  {
    MyGuard*     g4 = new MyGuard("Drawing stuff");
    TIter        next(&toDraw);
    GraphSysErr* g = 0;
    TH1*         f = 0;
    while ((g = static_cast<GraphSysErr*>(next()))) {
      TString lOpt(opt); lOpt.Append(" "); lOpt.Append(next.GetOption());
      MyPrint(Form("%s w/options \"%s\" (%d)",
		   g->GetName(),lOpt.Data(),
		   g->GetSumOption()));
      g->Draw(lOpt);
    }
    delete g4;
  }

  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs("QuickDraw.png");
  gPad = c;
}
