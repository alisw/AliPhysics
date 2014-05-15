TMap* qaConfig()
{
  TMap* c = new TMap();
  
  ConfigEntryNew(c, "meanTPCncl:run", 130., 140., "description");
  ConfigEntryNew(c, "meanTPCnclF:run", 0.8, 1.0, "description");
  ConfigEntryNew(c, "meanMIP:run", 47., 53., "description");
  ConfigEntryNew(c, "resolutionMIP:run", 0.002, 0.1, "description");
  
  ConfigEntryNew(c, "meanVertX:run", 0., 0.2, "description");
  ConfigEntryNew(c, "meanVertY:run", 0.2, 0.4, "description");
  ConfigEntryNew(c, "meanVertZ:run", -2., 2., "description");
  
  ConfigEntryNew(c, "offsetdRA:run", -0.2, 0.6, "description");
  ConfigEntryNew(c, "meanMultPos:run", 50., 100., "description");
  ConfigEntryNew(c, "tpcItsMatchA:run", 0.7, 1., "description");// TPC-ITS matching eff.
  ConfigEntryNew(c, "lambdaPull:run", -1., 1., "description"); // ITS-TPC matching eff.
  ConfigEntryNew(c, "tpcConstrainPhiA:run", -0.4, 0.4, "description");
  ConfigEntryNew(c, "deltaPt:run", -0.006, 0.006, "description");
  
  ConfigEntryNew(c, "dcarAP0:run", -0.6, 0.6, "description");
  ConfigEntryNew(c, "dcar_0:run", -0.02, 0.125, "description");
  ConfigEntryNew(c, "dcar_1:run", -0.08, 0.04, "description");
  ConfigEntryNew(c, "dcar_2:run", -0.04, 0.06, "description");
  ConfigEntryNew(c, "dcaz_0:run", -0.04, 0.04, "description");
  ConfigEntryNew(c, "dcaz_1:run", -0.01, 0.15, "description");
  ConfigEntryNew(c, "dcaz_2:run", -0.04, 0.04, "description");
   
  return c;
}

ConfigEntryNew(TMap* m, TString keyName, Float_t min, Float_t max, TString desc)
{
  TObjString* key = new TObjString(keyName);
  TObjString* description = new TObjString(desc);
  TVectorF* values = new TVectorF(2);
  values(0) = min; values(1) = max;
  TList* list = new TList();
  list->Add(values); //0
  list->Add(description); //1
  m->Add(key,list);
}

Float_t ConfigEntryMin(TMap* m, TNamed* h)
{
  TPair* p = dynamic_cast<TPair*>(m->FindObject(h->GetName()));
  TList* l = dynamic_cast<TList*>(p->Value());
  l->Print();
  if (!l) return -9999999;
  TVectorF* v = dynamic_cast<TVectorF*>(l->At(0));
  if (!v) return -1111111;  
  return v(0);
}

Float_t ConfigEntryMax(TMap* m, TNamed* h)
{
  TPair* p = dynamic_cast<TPair*>(m->FindObject(h->GetName()));
  TList* l = dynamic_cast<TList*>(p->Value());
  if (!l) return -9999999;
  TVectorF* v = dynamic_cast<TVectorF*>(l->At(0));
  if (!v) return 11111111;  
  return v(1);
}

TString ConfigEntryDescription(TMap* m, TNamed* h)
{
  TString dummy;
  TPair* p = dynamic_cast<TPair*>(m->FindObject(h->GetName()));
  TList* l = dynamic_cast<TList*>(p->Value());
  if (!l) return dummy;
  TObjString* os = dynamic_cast<TObjString*>(l->At(1));
  return os->GetString();
}

/*
example()
{
  gROOT->LoadMacro("qaConfig.C");
  TMap* configMap = qaConfig();

  TH1F* histi = new TH1F("test","test",1,0,1);
  Float_t min = ConfigEntryMin(configMap,hist);
  Float_t max = ConfigEntryMax(configMap,hist);
  TString desc = ConfigEntryDescription(configMap,hist);
  hist->Draw();
  TAxis* xaxis=hist->GetXaxis();
  Float_t x1 = xaxis->GetBinLowEdge(1);
  Float_t x2 = xaxis->GetBinUpEdge(xaxis->GetLast());
  TLine* lineMin = TLine(x1,min,x2,min); lineMin->SetLineColor(kRed);
  lineMin->Draw();
  TLine* lineMax = TLaxe(x1,max,x2,max); lineMax->SetLineColor(kRed);
  lineMax->Draw();
}
*/
