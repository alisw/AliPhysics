#include "AliAnalysisMuMuGraphUtil.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <vector>
#include <map>
#include "TMath.h"
#include "TObjArray.h"
#include "AliLog.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TStyle.h"
#include "AliAnalysisTriggerScalers.h"
#include <set>
#include "AliAnalysisMuMuResult.h"

ClassImp(AliAnalysisMuMuGraphUtil)

//____________________________________________________________________________
AliAnalysisMuMuGraphUtil::AliAnalysisMuMuGraphUtil(const char* ocdbpath) : TObject(),
fOCDBPath(ocdbpath),
fAttLine(),
fAttMarker(),
fAttFill(),
fAttXaxis(),
fAttYaxis(),
fDrawOptions(),
fShouldDrawPeriods(kFALSE)
{
  // default ctor
  DefaultStyle();
}

//____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuGraphUtil::Combine(TObjArray& graphs, Bool_t compact)
{
  // make one graph out of several
  // x axis is supposed to be run numbers and will end up ordered in the
  // returned graph
  
  std::map<int, std::vector<double> > values;
  std::map<int, std::vector<double> >::const_iterator it;

  TIter next(&graphs);
  TGraph* g;
  
  while ( ( g = static_cast<TGraph*>(next())) )
  {
    TGraphErrors* ge = dynamic_cast<TGraphErrors*>(g);
    
    for ( Int_t i = 0; i < g->GetN(); ++i )
    {
      Int_t runNumber = GetRunNumber(*g,i); // by doing this we "de-compact" the graph

      it = values.find(runNumber);
      if ( it != values.end() )
      {
        AliErrorClass(Form("Already got values for run %d !",runNumber));
        StdoutToAliErrorClass(graphs.Print(););
        return 0x0;
      }

      std::vector<double> quartet;
    
      quartet.push_back(runNumber);
      quartet.push_back(g->GetY()[i]);
      
      if ( ge )
      {
        quartet.push_back(ge->GetEX()[i]);
        quartet.push_back(ge->GetEY()[i]);
      }
      else
      {
        quartet.push_back(0.0);
        quartet.push_back(0.0);
      }
      
      values.insert( std::make_pair(runNumber,quartet));      
    }
  }
  
  TGraphErrors* rv(0x0);
  
  if ( values.size() )
  {
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vxerr;
    std::vector<double> vyerr;
    
    for ( it = values.begin(); it != values.end(); ++it )
    {
      const std::vector<double>& q = it->second;
      
      vx.push_back(q[0]);
      vy.push_back(q[1]);
      vxerr.push_back(q[2]);
      vyerr.push_back(q[3]);
    }
    
    rv = new TGraphErrors(values.size(),&vx[0],&vy[0],&vxerr[0],&vyerr[0]);
    rv->GetXaxis()->SetNoExponent();
    
    g = static_cast<TGraph*>(graphs.At(0));
    
    rv->SetName(g->GetName());
    rv->SetTitle(g->GetTitle());
    
    if ( compact || IsCompact(*g) )
    {
      Compact(*rv);
    }
  }
  
  return rv;
}

//____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::DefaultStyle()
{
  // Define default color/styles to be used, for at least 2 graphs
  // (here 5)
  
  Int_t colors[] = { 1, kGray+1, 4, 2, 6 };
  
  for ( Int_t i = 0; i < 5; ++i )
  {
    Int_t color = colors[i];
    
    fAttLine.push_back(TAttLine(color,1,1));
    fAttFill.push_back(TAttFill(color,1001));
    fAttMarker.push_back(TAttMarker(color,20+i,1));
    fAttXaxis.push_back(TAttAxis());
    
    TAttAxis a;
    
    a.ResetAttAxis();
    
    a.SetLabelColor(color);
    a.SetTitleColor(color);
    
    fAttYaxis.push_back(a);
    
    fDrawOptions.push_back("LP");
  }
}

//____________________________________________________________________________
TCanvas* AliAnalysisMuMuGraphUtil::DrawWith2Scales(TGraph& g1, TGraph& g2, const char* canvasName)
{
  TCanvas* c1 = new TCanvas(canvasName,canvasName);
  c1->Draw();
  
  TPad* pad1 = new TPad("pad1","",0,0,1,1);
  TPad* pad2 = new TPad("pad2","",0,0,1,1);
  
  g2.GetYaxis()->SetTitle(g2.GetTitle());
  g1.GetYaxis()->SetTitle(g1.GetTitle());
  
  pad1->SetFillStyle(4000);
  pad1->SetFrameFillStyle(0); //  transparent pad
  
  pad2->Draw();
  pad2->cd();
  
  StyleGraph(g2,1);
  
  g2.Draw("abxy+");
  
  pad1->Draw();
  pad1->cd();
  
  StyleGraph(g1,0);
  
  g1.Draw("alp");
  
  return c1;
}

//____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::Compact(TGraph& g)
{
  /// Compact (i.e. get the equivalent of 1 bin = 1 run number for an histogram)
  /// the graph.
  /// Only works if the x content of this graph represents run numbers. Otherwise
  /// result is unpredictable ;-)
  
  if ( !g.GetN() ) return;
  
  Double_t x,y;
  
  std::vector<int> runs;
  std::vector<double> bins;
  
  Int_t i(0);
  
  for ( i = 0; i < g.GetN(); ++i )
  {
    g.GetPoint(i,x,y);
    runs.push_back(TMath::Nint(x));
    bins.push_back(i);
    g.SetPoint(i,i+0.5,y);
  }
  
  bins.push_back(i);
  
  TAxis* axis = g.GetXaxis();
  
  axis->Set(g.GetN(),&bins[0]);
  
  for ( std::vector<int>::size_type j = 0; j < runs.size(); ++j )
  {
    axis->SetBinLabel(j+1,TString::Format("%d",runs[j]).Data());
  }
  
  axis->LabelsOption("v");
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuGraphUtil::GetRunNumber(const TGraph& g, Int_t i)
{
  // get the run number associated with bin i+1
  // if graph is not compacted then run number = x-value
  // otherwise we get it from the axis label
 
  Int_t runNumber = TMath::Nint(g.GetX()[i]);

  TString runLabel = g.GetXaxis()->GetBinLabel(i+1);
  
  if ( runLabel.Length() )
  {
    runNumber = runLabel.Atoi();
  }

  return runNumber;
}

//_____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::GetRuns(std::set<int>& runs, TGraph& graph) const
{
  // extract the list of runs in graph's x-axis
  
  for ( Int_t i = 0; i < graph.GetN(); ++i )
  {
    runs.insert(GetRunNumber(graph,i));
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::GetYMinAndMax(TGraph& graph, Double_t& ymin, Double_t& ymax)
{
  // find graph y-range
  // note that ymin and ymax *must* be initialized correctly outside
  // (this is done this way so that this method can be used easily
  // to get the range of a set of graphs)
  
  Double_t x,y;
  
  for ( Int_t i = 0; i < graph.GetN(); ++i )
  {
    graph.GetPoint(i,x,y);
    ymin = TMath::Min(ymin,y);
    ymax = TMath::Max(ymax,y);
  }
  
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuGraphUtil::IsCompact(TGraph& g)
{
  // whether the graph is compact or not
  Double_t delta(0.0);
  
  for ( Int_t i = 1; i < g.GetN(); ++i )
  {
    delta = TMath::Max(delta,g.GetX()[i] - g.GetX()[i-1]);
  }
  
  Bool_t hasLabels(kFALSE);

  for ( Int_t i = 1; ( i <= g.GetN() ) && ( !hasLabels ); ++i )
  {
    TString label(g.GetXaxis()->GetBinLabel(i));
    if ( label.Length() ) hasLabels = kTRUE;
  }
  
  return hasLabels && delta == 1.0;
}

//_____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::PlotSameWithLegend(TObjArray& a,
                                                  Double_t ymin, Double_t ymax) const
{
  // plot on same canvas
  if (!gPad) new TCanvas;
  
  Double_t xmin = TMath::Limits<Double_t>::Max();
  Double_t xmax = TMath::Limits<Double_t>::Min();
  
  TIter next(&a);
  TGraph* g;
  
  while ( ( g = static_cast<TGraph*>(next())))
  {
    xmin = TMath::Min(xmin,g->GetX()[0]);

    xmax = TMath::Max(xmax,g->GetX()[g->GetN()-1]);

  }

  TH2* hframe = new TH2F("hframe","hframe",100,xmin,xmax,100,ymin,ymax);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  hframe->GetXaxis()->SetNoExponent();

//  if ( IsCompact(g1) )
//  {
//    (*(hframe->GetXaxis()))=(*(g1.GetXaxis()));
//  }

  hframe->Draw();

  if ( fShouldDrawPeriods )
  {
    std::set<int> runs;
  
    next.Reset();
    
    while ( ( g = static_cast<TGraph*>(next())))
    {
      GetRuns(runs,*g);
    }
    
    AliAnalysisTriggerScalers ts(runs,fOCDBPath);
  
    ts.DrawPeriods(ymin,ymax,kGray);
    
    hframe->Draw("axissame");
  }

  next.Reset();
  
  Int_t i(0);
  TLegend* l = new TLegend(0.5,0.7,0.9,0.9); // fixme: how to get the legend position/size ?
  l->SetFillColor(0);
  
  while ( ( g = static_cast<TGraph*>(next())))
  {
    StyleGraph(*g,i);
    g->Draw(fDrawOptions[i].c_str());
    ++i;
    l->AddEntry(g,g->GetName(),fDrawOptions[0].c_str());
  }
  
  l->Draw();  

}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMuGraphUtil::RelDif(TGraph& ga, TGraph& gb)
{
  // compute the relative difference between two graphs
  
  std::vector<double> vx;
  std::vector<double> vxerr;
  std::vector<double> vy;
  std::vector<double> vyerr;

  for ( Int_t i = 0; i < ga.GetN(); ++i )
  {
    Double_t xa,xb,ya,yb;
  
    ga.GetPoint(i,xa,ya);
    gb.GetPoint(i,xb,yb);
  
    if ( xa != xb )
    {
      AliErrorClass(Form("Incompatible graphs : got xa=%e and xb=%e",xa,xb));
      return 0x0;
    }
  
    Double_t newvalue = 0.0;
  
    if ( TMath::Abs(xa) > 1E-12 )
    {
      newvalue = 100.0*( yb - ya ) / ya;
    }
  
    Double_t yerr = 0.0;
    
    if ( dynamic_cast<TGraphErrors*>(&ga) && dynamic_cast<TGraphErrors*>(&gb) )
    {
        yerr = newvalue*AliAnalysisMuMuResult::ErrorAB(ya,ga.GetEY()[i],
                                                       yb,gb.GetEY()[i]);
    }
  
    vx.push_back(xa);
    vxerr.push_back(0.5);
    vy.push_back(newvalue);
    vyerr.push_back(yerr);
  }
  
  
  return new TGraphErrors(vx.size(),&vx[0],&vy[0],&vxerr[0],&vyerr[0]);
}

//_____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::StyleGraph(TGraph& g, UInt_t index) const
{
  if ( index >= fAttFill.size() ) index = 0;
  
  static_cast<TAttFill&>(g) = fAttFill[index];
  static_cast<TAttLine&>(g) = fAttLine[index];
  static_cast<TAttMarker&>(g) = fAttMarker[index];
  
  g.GetYaxis()->SetLabelColor(fAttYaxis[index].GetLabelColor());
  g.GetYaxis()->SetTitleColor(fAttYaxis[index].GetTitleColor());
  
  //static_cast<TAttAxis&>((*g.GetYaxis())) = fAttYaxis[index];
}

//_____________________________________________________________________________
void AliAnalysisMuMuGraphUtil::UnCompact(TGraph& g)
{
  /// Reverse operation of the Compact method
  /// Only works if the labels of this graph represents run numbers. Otherwise
  /// result is unpredictable ;-)
  
  if ( !g.GetN() ) return;
  
  //  Int_t run1 = TString(g.GetXaxis()->GetBinLabel(1)).Atoi();
  //  Int_t run2 = TString(g.GetXaxis()->GetBinLabel(g.GetN())).Atoi();
  
  std::vector<double> runs;
  Int_t runNumber(-1);
  
  for ( Int_t i = 0; i < g.GetN(); ++i )
  {
    runNumber = TString(g.GetXaxis()->GetBinLabel(i+1)).Atoi();
    runs.push_back(runNumber*1.0);
  }
  
  runs.push_back(runNumber+1);
  
  g.GetXaxis()->Set(g.GetN(),&runs[0]);
  
  for ( Int_t i = 0; i < g.GetN(); ++i )
  {
    g.SetPoint(i,runs[i],g.GetY()[i]);
  }
  
  g.GetXaxis()->SetNoExponent();
  
}


