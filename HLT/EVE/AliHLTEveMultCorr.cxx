#include "AliHLTEveMultCorr.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <TCanvas.h>
#include "AliHLTHOMERBlockDesc.h"
#include "TList.h"

AliHLTEveMultCorr::AliHLTEveMultCorr(const char* name): AliHLTEveBase(name)
        ,fVzeroMultCanvas(0)
        ,fZdcMultCanvas(0)
        ,fTpcMultCanvas(0)
        ,fCorrCanvas(0)
        ,fEtCorrCanvas(0)
        ,fZdcVzeroSpdCorrCanvas(0)
{

}

AliHLTEveMultCorr::~AliHLTEveMultCorr()
{
}

void AliHLTEveMultCorr::ResetElements()
{

}

void AliHLTEveMultCorr::UpdateElements()
{
    fVzeroMultCanvas->Update();
    fZdcMultCanvas->Update();
    fTpcMultCanvas->Update();
    fCorrCanvas->Update();
    fEtCorrCanvas->Update();
    fZdcVzeroSpdCorrCanvas->Update();

}

void AliHLTEveMultCorr::ProcessBlock(AliHLTHOMERBlockDesc* block)
{
    TList *hlist = dynamic_cast<TList*>(block->GetTObject());
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
	      fVzeroMultCanvas = CreateCanvas("V0 Mult", "V0 Multiplicities");
	      fVzeroMultCanvas->Divide(2, 2);
	    }
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fVzeroMult"));
            AddHistogramToCanvas(oneDf, fVzeroMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fVzeroMultAC"));
            AddHistogramToCanvas(twoDf, fVzeroMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fVzeroFlaggedMult"));
            AddHistogramToCanvas(oneDf, fVzeroMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fVzeroFlaggedMultAC"));
            AddHistogramToCanvas(twoDf, fVzeroMultCanvas, cd);

	    cd = 1;
            //ZDC multiplicity hists
            if (!fZdcMultCanvas) 
	    {
	      fZdcMultCanvas = CreateCanvas("ZDC Mult", "ZDC Multiplicities");
	      fZdcMultCanvas->Divide(2, 4);
	    }
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcEzdc"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcEzn"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcEzp"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcEzem"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcNpart"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fZdcB"));
            AddHistogramToCanvas(oneDf, fZdcMultCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fZdcEzemEzdc"));
            AddHistogramToCanvas(twoDf, fZdcMultCanvas, cd);

	    cd = 1;
            // TPC multiplicity
            if (!fTpcMultCanvas) 
	    {
	      fTpcMultCanvas = CreateCanvas("TPC Mult", "TPC Multiplicities");
	      fTpcMultCanvas->Divide(2);
	    }
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fTpcNch2"));
            AddHistogramToCanvas(oneDf, fTpcMultCanvas, cd);
	    cd++;
            oneDf = dynamic_cast<TH1F*>(hlist->FindObject("fTpcNch3"));
            AddHistogramToCanvas(oneDf, fTpcMultCanvas, cd);

	    cd = 1;
            // Correlations
            if (!fCorrCanvas) 
	    {
	      fCorrCanvas = CreateCanvas("Correlations", "Multiplicity Correlations");
	      fCorrCanvas->Divide(2, 3);
	    }
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzdcNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzeroNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrSpdTpcNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzemNch"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzdcVzero"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzemVzero"));
            AddHistogramToCanvas(twoDf, fCorrCanvas, cd);
	    
	    cd = 1;
            // ET Correlations
            if (!fEtCorrCanvas) 
	    {
	      fEtCorrCanvas = CreateCanvas("E_{T} corr", "E_{T} Correlations");
	      fEtCorrCanvas->Divide(2, 3);
	    }
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrZdcTotEvsPhosTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrZdcTotEvsEmcalTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrZdcTotEvsTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzerovsPhosTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzerovsEmcalTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzerovsTotEt"));
            AddHistogramToCanvas(twoDf, fEtCorrCanvas, cd);

	    cd = 1;
            // ZDC,V0 vs SPD
            if (!fZdcVzeroSpdCorrCanvas) 
	    {
	      fZdcVzeroSpdCorrCanvas = CreateCanvas("ZDC/V0 vs SPD", "ZDC/V0 vs. SPD");
	      fZdcVzeroSpdCorrCanvas->Divide(2, 3);
	      
	    }
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzeroSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzeroASpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrVzeroCSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzdcSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzdcASpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
	    cd++;
            twoDf = dynamic_cast<TH2F*>(hlist->FindObject("fCorrEzdcCSpd"));
            AddHistogramToCanvas(twoDf, fZdcVzeroSpdCorrCanvas, cd);
        }
        else
        {
            std::cout << "ERROR: This block does not contain what you think it contains!" << std::endl;
        }
    }
    else
    {
      std::cout << "ERROR: There is no TList object in the block" << std::endl;
    }


}

void AliHLTEveMultCorr::AddHistogramsToCanvas(AliHLTHOMERBlockDesc* /*block*/, TCanvas* /*canvas*/, Int_t& /*cdCount*/)
{
}

void AliHLTEveMultCorr::AddHistogramToCanvas(TH1* hist, TCanvas* canvas, Int_t& cdCount)
{
    if (hist)
    {
        if (!strcmp(hist->ClassName(), "TH1F"))
        {
	    canvas->cd(cdCount);
            dynamic_cast<TH1F*>(hist)->Draw();
        }
        else if (!strcmp(hist->ClassName(), "TH2F"))
        {
	    canvas->cd(cdCount);
            dynamic_cast<TH2F*>(hist)->Draw("COLZ");
        }
        else
        {
            std::cout << "I don't want histograms of type " << hist->ClassName() << std::endl;
        }
    }
    else
    {
      std::cout << "ERROR (AddHistogramsToCanvas): You gave me a null pointer" << std::endl;
    }
}

