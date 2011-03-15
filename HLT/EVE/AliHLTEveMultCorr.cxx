#include "AliHLTEveMultCorr.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <TCanvas.h>
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTEveHistoMerger.h"
#include "TList.h"
#include "AliLog.h"

AliHLTEveMultCorr::AliHLTEveMultCorr(const char* name): AliHLTEveBase(name)
        ,fVzeroMultCanvas(0)
        ,fZdcMultCanvas(0)
        ,fTrackMultCanvas(0)
        ,fCorrCanvas(0)
        ,fEtCorrCanvas(0)
        ,fZdcVzeroSpdCorrCanvas(0)
        ,fMerger(0)
        ,fMyList(0)
{
  fMerger = new AliHLTEveHistoMerger();
}

AliHLTEveMultCorr::~AliHLTEveMultCorr()
{
  if(fMerger)
  {
    delete fMerger;
  }
  if(fMyList)
  {
    delete fMyList;
  }
}

void AliHLTEveMultCorr::ResetElements()
{

}

void AliHLTEveMultCorr::UpdateElements()
{
    fVzeroMultCanvas->Update();
    fZdcMultCanvas->Update();
    fTrackMultCanvas->Update();
    fCorrCanvas->Update();
    fEtCorrCanvas->Update();
    fZdcVzeroSpdCorrCanvas->Update();

}

void AliHLTEveMultCorr::ProcessBlock(AliHLTHOMERBlockDesc* block)
{
    TString msg;
    
    TList *hlist = dynamic_cast<TList*>(block->GetTObject());
    if(!hlist)
    {
      hlist = dynamic_cast<TList*>(fMerger->Process(hlist, block->GetSpecification()));
    }
    else
    {
      if(fMyList)
      {
	delete fMyList;
	fMyList = 0;
      }
      fMyList = dynamic_cast<TList*>(hlist->Clone());
      hlist = fMyList;
    }
    if (hlist)
    {

        if (hlist->Contains("fVzeroMult")) // These are the correlation histograms
        {
            TH1F *oneDf = 0;
            TH2F *twoDf = 0;
	    int cd = 1;

            //VZERO multiplicity hists
            if (!fVzeroMultCanvas) 
	    {
	      fVzeroMultCanvas = CreateCanvas("V0 M", "V0 Multiplicities");
	      fVzeroMultCanvas->Divide(2, 2);
	    }
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fVzeroMult"));
            AddHistogramToCanvas(oneDf, fVzeroMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fVzeroMultAC"));
            AddHistogramToCanvas(twoDf, fVzeroMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fVzeroFlaggedMult"));
            AddHistogramToCanvas(oneDf, fVzeroMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fVzeroFlaggedMultAC"));
            AddHistogramToCanvas(twoDf, fVzeroMultCanvas, cd);

	    cd = 1;
            //ZDC multiplicity hists
            if (!fZdcMultCanvas) 
	    {
	      fZdcMultCanvas = CreateCanvas("ZDC M", "ZDC Multiplicities");
	      fZdcMultCanvas->Divide(3, 3);
	    }
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcEzdc"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcEzn"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcEzp"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcEzem"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcNpart"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fZdcB"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fZdcEzemEzdc"));
            AddHistogramToCanvas(twoDf, fZdcMultCanvas, cd);

	    cd = 1;
            // TPC multiplicity
            if (!fTrackMultCanvas) 
	    {
	      fTrackMultCanvas = CreateCanvas("Track M", "Track Multiplicities");
	      fTrackMultCanvas->Divide(2, 2);
	    }
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fTpcNch2"));
            AddHistogramToCanvas(oneDf, fTrackMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fTpcNch3"));
            AddHistogramToCanvas(oneDf, fTrackMultCanvas, cd);
	    cd++;
	    oneDf = dynamic_cast<TH1F*>(FindHistogram(hlist, "fSpdNClusters"));
	    
            AddHistogramToCanvas(oneDf, fTrackMultCanvas, cd);

	    cd = 1;
            // Correlations
            if (!fCorrCanvas) 
	    {
	      fCorrCanvas = CreateCanvas("Corr", "Multiplicity Correlations");
	      fCorrCanvas->Divide(3, 2);
	    }
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzdcNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzeroNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrSpdTpcNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzemNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzdcVzero"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzemVzero"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    
	    cd = 1;
            // ET Correlations
            if (!fEtCorrCanvas) 
	    {
	      fEtCorrCanvas = CreateCanvas("ET", "E_{T} Correlations");
	      fEtCorrCanvas->Divide(3, 2);
	    }
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrZdcTotEvsPhosTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrZdcTotEvsEmcalTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrZdcTotEvsTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzerovsPhosTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzerovsEmcalTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzerovsTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);

	    cd = 1;
            // ZDC,V0 vs SPD
            if (!fZdcVzeroSpdCorrCanvas) 
	    {
	      fZdcVzeroSpdCorrCanvas = CreateCanvas("ZDC/V0 vs SPD", "ZDC/V0 vs. SPD");
	      fZdcVzeroSpdCorrCanvas->Divide(3, 2);
	      
	    }
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzeroSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzeroASpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrVzeroCSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzdcSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzdcASpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(FindHistogram(hlist, "fCorrEzdcCSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
        }
        else
        {
            AliWarning("This block does not contain what you think it contains!");
        }
    }
    else
    {
      AliWarning("There is no TList object in the block");
    }


}

void AliHLTEveMultCorr::AddHistogramsToCanvas(AliHLTHOMERBlockDesc* /*block*/, TCanvas* /*canvas*/, Int_t& /*cdCount*/)
{
}

void AliHLTEveMultCorr::AddHistogramToCanvas(TH1* hist, TCanvas* canvas, Int_t& cdCount, Bool_t zoom)
{
  
  TString msg;
    if (hist)
    {
        if (!strcmp(hist->ClassName(), "TH1F"))
        {
	    canvas->cd(cdCount);
	    TPad * pad = dynamic_cast<TPad*>(canvas->cd(cdCount));
	    if(pad) pad->SetLogy();
	    if(zoom) 
	    {
	      TH1F *h = dynamic_cast<TH1F*>(hist);
	      if(h) 
		{
		  h->GetXaxis()->SetRange(0, (Int_t) (h->GetMaximumBin() + h->GetMaximumBin()*0.2));
		  h->Draw();
		}
	    }
        }
        else if (!strcmp(hist->ClassName(), "TH2F"))
        {
	    canvas->cd(cdCount);
	    if(zoom)
	    {
	      TH2F *h = dynamic_cast<TH2F*>(hist);
	      if(h)
		{
		  h->GetXaxis()->SetRange(0, (Int_t) (h->GetMaximumBin() + h->GetMaximumBin()*0.2));
		  h->GetYaxis()->SetRange(0, (Int_t) (h->GetMaximumBin() + h->GetMaximumBin()*0.2));
		}
	    }
            dynamic_cast<TH2F*>(hist)->Draw("COLZ");
        }
        else
        {	
	  msg.Form("I don't want histograms of type %s", hist->ClassName());
	  AliWarning(msg.Data());
        }
    }
    else
    {
      AliWarning("You gave me a null pointer!");
    }
}

TH1* AliHLTEveMultCorr::FindHistogram(TCollection* coll, const char* name)
{
  TString msg;
  TH1 *hist = dynamic_cast<TH1*>(coll->FindObject(name));
  if(!hist)
  {
    msg.Form("Could not find object %s", name);
    AliWarning(msg);
  }
  return hist;
}
