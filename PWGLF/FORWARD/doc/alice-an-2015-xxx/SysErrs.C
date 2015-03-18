void
SysErrs(const TString& sys="pp", UShort_t sNN=900, const TString& trig="INEL")
{
  gROOT->LoadMacro("~/GSE-0.1/GraphSysErr.C+g");
  
  GraphSysErr* g = new GraphSysErr("sysErr", "Systematic errors");
  
  Double_t eventSelectionLow  = 0;
  Double_t eventSelectionHigh = 0;
  
  if (sys.EqualTo("pp")) { 
    if (trig.EqualTo("INEL")) { 
      switch (sNN) { 
      case 900:  eventSelectionLow = 0.1;  eventSelectionHigh = 0.3; break;
      case 2760: eventSelectionLow = 0.35; eventSelectionHigh = 0.6; break;
      case 7000: eventSelectionLow = 0.3;  eventSelectionHigh = 0.6; break;
      case 8000: eventSelectionLow = 0.3;  eventSelectionHigh = 0.6; break;
      }
    }
    else if (trig.EqualTo("NSD")) { 
      eventSelectionLow = eventSelectionHigh = (sNN == 2760 ? 3 : 2);
    }
  }
  else if (sys.EqualTo("pPb")) { 
    if (trig.EqualTo("MB")) 
      eventSelectionLow = eventSelectionHigh = 2;
    else 
      eventSelectionLow = eventSelectionLow = 4;
  }
  else if (sys.EqualTo("Pbp")) { 
    if (trig.EqualTo("MB")) 
      eventSelectionLow = eventSelectionHigh = 2;
    else 
      eventSelectionLow = eventSelectionLow = 6;
  }
  else if (sys.EqualTo("PbPb"))
    eventSelectionLow = eventSelectionHigh = 2;

  Int_t es = g->DefineCommon("Event selection", false, eventSelectionLow, 
			     eventSelectionHigh);
  Int_t sf = g->DeclarePoint2Point("Merging", false);
  Int_t dc = g->DeclarePoint2Point("Density", false);
  Int_t sa = g->DeclarePoint2Point("Satellite", false);
  Int_t re = g->DeclarePoint2Point("Reference", false);

  g->SetSysFillColor(es, kRed+2);
  g->SetSysFillStyle(es, 3001);
  g->SetSysLineColor(es, kRed+2);
  g->SetSysOption(es, GraphSysErr::kFill);

  g->SetSysFillColor(sf, kCyan+2);
  g->SetSysFillStyle(sf, 3001);
  g->SetSysOption(sf, GraphSysErr::kFill);

  g->SetSysLineColor(dc, kMagenta+2);
  g->SetSysLineWidth(dc, 2);
  g->SetSysFillColor(dc, kMagenta+2);
  g->SetSysFillStyle(dc, 3001);
  g->SetSysOption(dc, GraphSysErr::kFill);


  g->SetSysLineColor(sa, kBlue+2);
  g->SetSysFillColor(sa, kBlue+2);
  g->SetSysLineWidth(sa, 2);
  g->SetSysFillStyle(sa, 3001);
  g->SetSysOption(sa, GraphSysErr::kFill);

  g->SetSysLineColor(re, kGreen+2);
  g->SetSysFillColor(re, kGreen+2);
  g->SetSysLineWidth(re, 3);
  g->SetSysFillStyle(re, 3001);
  g->SetSysOption(re, GraphSysErr::kFill);


  TAxis a(40, -4, 6);
  for (Int_t i = 0; i < 40; i++) { 
    Double_t eta  = a.GetBinCenter(i+1);
    Double_t dEta = a.GetBinWidth(i+1);
    
    g->SetPoint(i, eta, 0);
    g->SetPointError(i, dEta);

    g->SetSysError(sf, i, 0, 1);
    g->SetSysError(dc, i, 0, 1);
    g->SetSysError(sa, i, 0, 5.3);
    g->SetSysError(re, i, 0, 5);
  }

  g->Draw("STACK COMMON STAT quad axis");
  
  g->Export();
}

    
    
  
