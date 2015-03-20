#ifndef __CINT__
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TString.h>
#include <TError.h>
#include <TClass.h>
#include <fstream>
#include <iomanip>
#else 
class TGraphAsymmErrors;
#endif

TGraphAsymmErrors* GetGraph(UShort_t sNN, const TString& trig) {
  TString fname(Form("%s_%04d.root", trig.Data(), sNN));
  TFile* file = TFile::Open(fname.Data(), "READ");
  if (!file) { 
    Error("", "Failed to open %s", fname.Data());
    return 0;
  }
  
  TString objName(trig);
  if (objName == "INELGt0") objName = "INELg0";
  
  TObject* o = file->Get(objName);
  if (!o) { 
    Error("", "Failed to find %s in %s", objName.Data(), fname.Data());
    file->Close();
    return 0;
  }

  if (!o->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) { 
    Error("", "%s is not a TGraphAsymmErrors, but a %s", 
	  objName.Data(), o->ClassName());
    file->Close();
    return 0;
  }
  TGraphAsymmErrors* ret = static_cast<TGraphAsymmErrors*>(o);
  file->Close();

  return ret;
}
void WriteArray(const char* name, Double_t* a, Int_t first, Int_t last, std::ostream& o) 
{

  o << "    Double_t " << name << "[] = {";
  Int_t n = last-first+1;
  for (Int_t i = 0; i < n; i++) {
    if (i != 0) {
      o << ",";
      if (i % 5 == 0) o << "\n    ";
    }
    o << a[first+i];
  }
  o << "};" << std::endl;
} 
      
void WriteGraph(UShort_t sNN, const TString& trig, const TString& date, std::ostream& o) 
{
  TGraphAsymmErrors* src = GetGraph(sNN, trig);
  if (!src) return;
  
  o << "// sNN=" << sNN << " trig=" << trig << std::endl;
  
  Int_t     maxN  = src->GetN();
  Double_t* px    = src->GetX();
  Double_t* pxel  = src->GetEXlow();
  Double_t* pxeh  = src->GetEXhigh();
  Double_t* py    = src->GetY();
  Double_t* pyel  = src->GetEYlow();
  Double_t* pyeh  = src->GetEYhigh();
  Int_t     first = -1;
  Int_t     last  = -1;
  for (Int_t i = 0; i < maxN; i++) { 
    if (py[i] > 1e-9) { 
      last = i;
      if (first < 0) first = i;
    }
  }
  if (first < 0 || last < 0) { 
    Warning("", "No points found for %d/%s", sNN, trig.Data());
    return;
  }
  Int_t  n = last-first+1;
  Printf("Got %d-%d+1=%d non-zero points for %s",
	 last,first,last-first+1, src->GetName());

  o << "  if (sNN==" << sNN << " && trg.EqualTo(\"" << trig << "\",TString::kIgnoreCase))  {\n"
    << std::scientific 
    << "    const Int_t n = " << n << ";" << std::endl;
  WriteArray("x",   px,   first, last, o);
  WriteArray("xel", pxel, first, last, o);
  WriteArray("xeh", pxeh, first, last, o);
  WriteArray("y",   py,   first, last, o);
  WriteArray("yel", pyel, first, last, o);
  WriteArray("yeh", pyeh, first, last, o);

  TString otrg(trig);
  if      (otrg.EqualTo("INEL"))   otrg = "Inel";
  else if (otrg.EqualTo("NSD"))    otrg = "Nsd";
  else if (otrg.EqualTo("INELGt0")) otrg = "InelGt0";

  o << "    TGraphAsymmErrors* g = new TGraphAsymmErrors(" << n << ",x,y,xel,xeh,yel,yeh);\n"
    << "    SetGraphAttributes(g," << trig << ",WIP,false,\n"
    << "                       \"alice_pp" << sNN << otrg << "Work\",\n"
    << "                       \"PWG-UD/MULT - " << date << "\");\n"
    << "    return g;\n"
    << "  }"
    << std::endl;
}

void Extract()
{
  TString date("2014/05/01");
  std::ofstream o("PWGUD.C");
  
  o << "Int_t INEL = kRed+2;\n"
    << "Int_t NSD = kGreen+2;\n"
    << "Int_t INELGt0 = kBlue+2;\n"
    << "Int_t WIP = 20;\n"
    << "void SetGraphAttributes(TGraph* g, Int_t t, Int_t exp, bool b, const char* name, const char* title) {\n"
    << "  g->SetMarkerColor(t);\n"
    << "  g->SetMarkerStyle(exp);\n"
    << "  g->SetLineColor(t);\n"
    << "  g->SetFillColor(t);\n"
    << "  g->SetName(name);\n"
    << "  g->SetTitle(title);\n"
    << "  g->SetFillStyle(0);\n"
    << "  g->GetHistogram()->SetStats(kFALSE);\n"
    << "  g->GetHistogram()->SetXTitle(\"#eta\");\n"
    << "  g->GetHistogram()->SetYTitle(\"#frac{1}{N} #frac{dN_{ch}}{#eta}\");\n"
    << "}\n\n"
    << "TGraphAsymmErrors* GetPWGUD(UShort_t sNN, const TString& trg)\n"
    << "{\n";
  UShort_t aSNN[] = { 900, 2760, 7000, 8000, 0 };
  const Char_t* aTrig[] = { "INEL", "NSD", "INELGt0", 0 };
  UShort_t* pSNN = aSNN;

  while ((*pSNN > 0)) { 
    const char** pTrig = aTrig;
    while ((*pTrig)) {
      WriteGraph(*pSNN, *pTrig, date, o);
      pTrig++;
    }
    pSNN++;
  }
  
  o << "  Warning(\"\",\"No data for sNN=%d and trig=%s\",sNN,trg.Data());\n" 
    << "  return 0;\n"
    << "}\n"
    << "\n"
    << "TMultiGraph* PWGUD() {\n"
    << "  TMultiGraph* mg = new TMultiGraph();\n";
  pSNN = aSNN;
  while ((*pSNN > 0)) { 
    const char** pTrig = aTrig;
    while ((*pTrig)) {
      o << "  mg->Add(GetPWGUD(" << *pSNN << ",\"" << *pTrig << "\"));\n";
      pTrig++;
    }
    pSNN++;
  }
  o << "\n"
    << "  mg->Draw(\"APL\");\n"
    << "  return mg;\n"
    << "}\n"
    << std::endl;
}
