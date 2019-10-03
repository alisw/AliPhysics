#include "AliMuonCompactQuickAccEffChecker.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TParameter.h"
#include "TString.h"
#include <cassert>
#include <iostream>
#include <vector>

/// \ingroup compact
void AliMuonCompactQuickAccEffChecker::CompareWithFullAccEff(const char* fullacceff, const char* quickacceff, const char* quickAccEffWithoutMonoCathodeClusters)
{
    // mind your steps here.
    // this method is far for being bullet-proof.
    // some checks are missing, in particular in the case quickAccEffWithoutMonoCathodeCluster is provided,
    // one should check the reference number of jpsi, etc... are the same as the ones in
    // quickacceff file...
    
    TFile* f = TFile::Open(fullacceff);
    TH1* hfullacceff = static_cast<TH1*>(f->Get("Eff"));
    hfullacceff->SetDirectory(0);
    delete f;

    f = TFile::Open(quickacceff);

    std::vector<TGraphErrors*> vg;
    
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_occ_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_ped_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_hv_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_hv_occ_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_ped_hv_occ_config")));
    vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_ped_hv_occ_config_reject")));

    TParameter<Double_t>* RefAccEff = static_cast<TParameter<Double_t>*>(f->Get("RefAccEff"));
    TParameter<Double_t>* RefAccEffError = static_cast<TParameter<Double_t>*>(f->Get("RefAccEffError"));

    delete f;

    if  (strlen(quickAccEffWithoutMonoCathodeClusters))
    {
        f = TFile::Open(quickAccEffWithoutMonoCathodeClusters);
        vg.push_back(static_cast<TGraphErrors*>(f->Get("acceffdrop_ped_hv_occ_config_reject")));
        vg[vg.size()-1]->SetName("acceffdrop_ped_hv_occ_config_reject-no-mono");
        delete f;
    }
    std::vector<TH1*> vquick;
    std::vector<TH1*> vratio;

    TAxis* x = hfullacceff->GetXaxis();

    std::vector<int> colors = { 2,4,6,8, kRed-7, kRed-10, kBlue-3, kBlue-8 };

    TCanvas* ccompare = new TCanvas("ccompare","ccompare");

    int icolor = 0;

    hfullacceff->Draw("histe");

    TLegend* legend = new TLegend(0.6,0.1,0.9,0.3);

    std::vector<TH1*> vratiomean;

    double refAccEff = RefAccEff->GetVal(); // 0.187829; 
    double refAccEffRelativeError = RefAccEffError->GetVal()/refAccEff;

    for ( auto g : vg )
    {
        if (!g) continue;

        TH1* hQuickAccEff = static_cast<TH1*>(hfullacceff->Clone(Form("h_%s",g->GetName())));
        hQuickAccEff->Reset();
        hQuickAccEff->SetMarkerStyle(0);

        TH1* hratiomean = new TH1F(Form("hratiomean_%s",g->GetName()),"",200,0.5,1.5);

        hratiomean->SetLineColor(colors[icolor]);

        vratiomean.push_back(hratiomean);

        for ( Int_t i = 1; i <= x->GetNbins(); ++i )
        {
            TString srn = x->GetBinLabel(i);
            Int_t rn = static_cast<Int_t>(g->GetX()[i-1]);
            assert(rn=srn.Atoi());
            double drop = g->GetY()[i-1];
            double dropRelativeError = g->GetEY()[i-1]/drop;
            double quick = refAccEff*(1 - drop/100.0 );
            double quickRelativeError = TMath::Sqrt(dropRelativeError*dropRelativeError + refAccEffRelativeError*refAccEffRelativeError);
            hQuickAccEff->SetBinContent(i,quick);
            hQuickAccEff->SetBinError(i,quick*quickRelativeError);
            hratiomean->Fill(quick/hfullacceff->GetBinContent(i));
            std::cout << Form("RUN %6d %30s DROP %7.2f %% +- %5.2f %% QuickAccEff %7.2f %% +- %5.2f %% FullAccEff %7.2f %% +- %5.2f %%",
                    rn,
                    g->GetName(),
                    drop,dropRelativeError*drop,
                    100*quick,100*quickRelativeError*quick,
                    100*hfullacceff->GetBinContent(i),100*hfullacceff->GetBinError(i)) << std::endl;
        }

        legend->AddEntry(hQuickAccEff,g->GetName(),"l");
        hQuickAccEff->SetLineColor(colors[icolor]);
        hQuickAccEff->Draw("histesame");
        TH1* h = static_cast<TH1*>(hfullacceff->Clone(Form("hratio_%s",g->GetName())));
        h->Divide(hQuickAccEff);
        vratio.push_back(h);
        h->SetLineColor(colors[icolor]);
        ++icolor;
    }

    legend->Draw();

    TCanvas* canvas = new TCanvas("cratio","cratio");

    bool first = true;

    for ( auto h : vratio)
    {
        if ( !TString(h->GetName()).Contains("acceffdrop_ped_hv_occ_config_reject")) continue;

        h->SetMinimum(0.50);
        h->SetMaximum(1.15);
        if (first)
        {
            h->Draw("hist");
            first=false;
        }
        else
        {
            h->Draw("histsame");
        }
    }

    TCanvas* cratiomean = new TCanvas("cratiomean","cratiomean");

    first = false;

    for ( auto h : vratiomean)
    {
        if ( !TString(h->GetName()).Contains("acceffdrop_ped_hv_occ_config_reject")) continue;

        std::cout << h->GetName() << " Mean=" << h->GetMean() << std::endl;

        if (first)
        {
            h->Draw("hist");
            first=false;
        }
        else
        {
            h->Draw("histsame");
        }
    }
}
