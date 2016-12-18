#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <Riostream.h>
#endif



void PlotMatchEffSystUncGraph(){
    
    //good to plot vs pt
    //to be updated to plot vs phi, eta
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111111);
    
    Double_t Eff_mc_prim[14];
    Double_t Eff_mc_sec[14];
    Double_t Eff_mc_tot[14];
    Double_t Eff_dati[14];
    
    Double_t frazPrim[14];
    Double_t TPCOnly[14];
    
    Double_t erD[14];
    Double_t erP[14];
    Double_t erS[14];
    Double_t erfP[14];
    Double_t erTPC[14];

    
    //pidfix, golden
    ifstream primaries("MatchEff_MCprim_15o_pidfix_Golden.txt",ios::in);
    ifstream primErr("MatchEff_MCprim_err_15o_pidfix_Golden.txt",ios::in);
    ifstream secondaries("MatchEff_MCsec_15o_pidfix_Golden.txt",ios::in);
    ifstream secErr("MatchEff_MCsec_err_15o_pidfix_Golden.txt",ios::in);
    ifstream inclusive("MatchEff_MCall_15o_pidfix_Golden.txt",ios::in);
    ifstream data("MatchEff_data_15o_pidfix_Golden.txt",ios::in);
    ifstream dataErr("MatchEff_data_err_15o_pidfix_Golden.txt",ios::in);
    
    ifstream frazP("fractions__PiK__15o_pidfix_Golden_VsPt.txt",ios::in);
    ifstream TPCOnlyCorr("corrTPC_VsPt_15o_pidfix_Golden.txt",ios::in);
    ifstream frazPerr("fractionErrs__PiK__15o_pidfix_Golden_VsPt.txt",ios::in);
    ifstream TPCOnlyCorrErr("corrTPCErr_VsPt_15o_pidfix_Golden.txt",ios::in);

    if (primaries.is_open())
        for(Int_t i=0;i<14;i++){
            
            primaries >> Eff_mc_prim[i];
            primErr >> erP[i];
            secondaries >> Eff_mc_sec[i];
            secErr >> erS[i];
            inclusive >> Eff_mc_tot[i];
            data >> Eff_dati[i];
            dataErr >> erD[i];
            
            frazP >> frazPrim[i];
            TPCOnlyCorr >> TPCOnly[i];
            frazPerr >> erfP[i];
            TPCOnlyCorrErr >> erTPC[i];
            
            //    primariesErr >> erP[i];
            //    secondariesErr >> erS[i];
        }
    primaries.close();
    primErr.close();
    secondaries.close();
    secErr.close();
    inclusive.close();
    data.close();
    dataErr.close();
    frazP.close();
    TPCOnlyCorr.close();
    frazPerr.close();
    TPCOnlyCorrErr.close();
    //primariesErr.close();
    //secondariesErr.close();
    
    
    Double_t erTmc[14];
    Double_t erSys[14];
    Double_t x[14]={0.75,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5};
    Double_t erx[14]={0.1,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
    
    TGraphErrors *gr_Eff_dati=new TGraphErrors(14,x,Eff_dati,erx,erD);
    
    TGraphErrors *gr_Eff_mc_prim=new TGraphErrors(14,x,Eff_mc_prim,erx,erP);
    
    TGraphErrors *gr_Eff_mc_sec=new TGraphErrors(14,x,Eff_mc_sec,erx,erS);
    
    TGraph *gr_Eff_mc_tot=new TGraph(14,x,Eff_mc_tot);
    
    gr_Eff_dati->SetTitle("effData");
    gr_Eff_mc_prim->SetTitle("eff_mc_prim");
    gr_Eff_mc_sec->SetTitle("eff_mc_sec");
    gr_Eff_mc_tot->SetTitle("eff_mc_tot");
    
    gr_Eff_dati->SetMarkerColor(kBlue);
    gr_Eff_mc_prim->SetMarkerColor(kGreen);
    gr_Eff_mc_sec->SetMarkerColor(kRed);
    gr_Eff_mc_tot->SetMarkerColor(kGreen+2);
    gr_Eff_dati->SetLineColor(kBlue);
    gr_Eff_mc_prim->SetMarkerColor(kGreen);
    gr_Eff_mc_prim->SetLineColor(kGreen);
    gr_Eff_mc_sec->SetMarkerColor(kRed);
    gr_Eff_mc_sec->SetLineColor(kRed);
    gr_Eff_mc_tot->SetMarkerColor(kGreen+2);
    gr_Eff_mc_tot->SetLineColor(kGreen+2);
    
    
    gr_Eff_dati->SetFillColor(0);
    gr_Eff_mc_prim->SetFillColor(0);
    gr_Eff_mc_sec->SetFillColor(0);
    gr_Eff_mc_tot->SetFillColor(0);
    
    Double_t effMCcorr[14];
    Double_t systMCcorr[14];
    Double_t systMC[14];
    for(Int_t i=0;i<14;i++){
        effMCcorr[i]=TPCOnly[i]*frazPrim[i]*Eff_mc_prim[i]+(1-TPCOnly[i]*frazPrim[i])*Eff_mc_sec[i];
        erTmc[i] = TMath::Sqrt((TPCOnly[i]*Eff_mc_prim[i]-TPCOnly[i]*Eff_mc_sec[i])*(TPCOnly[i]*Eff_mc_prim[i]-TPCOnly[i]*Eff_mc_sec[i])*erfP[i]*erfP[i] +
                               ((frazPrim[i]*Eff_mc_prim[i]-frazPrim[i]*Eff_mc_sec[i])*(frazPrim[i]*Eff_mc_prim[i]-frazPrim[i]*Eff_mc_sec[i]))*erTPC[i]*erTPC[i] +
                               (frazPrim[i]*TPCOnly[i]*frazPrim[i]*TPCOnly[i])*erP[i]*erP[i] +
                               (1 - TPCOnly[i]*frazPrim[i])*(1 - TPCOnly[i]*frazPrim[i])*erS[i]*erS[i]);
        
        systMCcorr[i]=effMCcorr[i]/Eff_dati[i];
        
        
        systMC[i] = Eff_mc_tot[i]/Eff_dati[i];
        erSys[i]=TMath::Sqrt(1/(Eff_dati[i]*Eff_dati[i])*erTmc[i]*erTmc[i]+(effMCcorr[i]/(Eff_dati[i]*Eff_dati[i]))*(effMCcorr[i]/(Eff_dati[i]*Eff_dati[i]))*erD[i]*erD[i]);
    }
    
    TLine *l = new TLine(0.5,1.02,14,1.02);
    TLine *l2= new TLine(0.5,0.98,14,0.98);
    TCanvas *c=new TCanvas("syst","syst");
    TGraphErrors *gr_Eff_systCorr=new TGraphErrors(14,x,systMCcorr,erx,erSys);
    TGraph *gr_Eff_syst=new TGraph(14,x,systMC);
    TF1 *fitfun = new TF1("fitfun","pol0",0,13.5);
    gr_Eff_systCorr->SetTitle("systMCcorr");
    gr_Eff_syst->SetTitle("systMC");
    gr_Eff_systCorr->SetMarkerColor(kBlue);
    gr_Eff_syst->SetMarkerColor(kRed);
    gr_Eff_systCorr->Fit(fitfun,"R");
    gr_Eff_systCorr->SetFillColor(0);
    gr_Eff_syst->SetFillColor(0);
    gr_Eff_systCorr->Draw("A* ");
    gr_Eff_syst->Draw("* same");
    l->Draw("same");
    l2->Draw("same");
    
    TCanvas *ci3=new TCanvas("eff","eff");
    TGraphErrors *gr_Eff_Corr=new TGraphErrors(14,x,effMCcorr,erx,erTmc);
    TGraph *gr_Eff_NOTCORR=new TGraph(14,x,Eff_mc_tot);
    gr_Eff_Corr->SetTitle("effMCcorr");
    gr_Eff_NOTCORR->SetTitle("effMC");
    gr_Eff_Corr->SetMarkerColor(kBlack);
    gr_Eff_Corr->SetLineColor(kBlack);
    gr_Eff_NOTCORR->SetLineColor(kMagenta);
    gr_Eff_NOTCORR->SetMarkerColor(kMagenta);
    gr_Eff_Corr->SetFillColor(0);
    gr_Eff_NOTCORR->SetFillColor(0);
    
    gr_Eff_Corr->GetYaxis()->SetRangeUser(0.,1.);
    gr_Eff_Corr->Draw("A* ");
    gr_Eff_NOTCORR->Draw("* same");
    gr_Eff_dati->Draw("* same");
    gr_Eff_mc_prim->Draw("* same");
    gr_Eff_mc_sec->Draw("* same");
}
