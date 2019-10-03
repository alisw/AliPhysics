#include "extractV2.C"
#include "TMath.h"
Bool_t CalculateSystError=1;
Bool_t CalculateSystErrorOnV2=1;
Bool_t CalculateSystErrorOnRatio=1;
Bool_t CalculateSystErrorOnDiff=1;
const Int_t nc = 3;

void syst_mu_tkl(){
  gStyle->SetOptTitle(0);
  Int_t     nbinst;
  Double_t* binst;
  Double_t binst_centers[100];

  for (Int_t im=0;im<nm;im++){ // cent method
//                if(im <= kNtkl)continue; //only AMPT
//                if(im<=kDPM)continue; //single estimator
    if (im!=kDPM) continue;
//                if (im!=kAMM) continue;
    //    if(im != kLEE && im != kV0S)continue; //single estimator
    TString centMethod = centMethods[im];
    Bool_t isTkl=im>kNtrk && im<=kNtkl? 1 : 0;
    Bool_t isAMPTgen=im>kNtkl? 1 : 0;
    Int_t assBin = isTkl ? 4 : 0;
    Int_t     nbinsa = isTkl ? nbins_tkl      : nbins_trk;
    Double_t* binsa  = isTkl ? bins_tkl       : bins_trk;

    for (Int_t ip=0;ip<2;ip++){ // period

      if(CalculateSystError){ //calculate syst check
        for (Int_t is=0;is<nsyst;is++){
          printf("ip=%i im=%i is=%i\n",ip,im,is);
          //          if(is!=k_LowMultScale)continue;
          gStudySystematic=isystCheck[is];
          TString inTkFileName =Form("dphi/dphi_corr_%s_%s_%s.root","tkl_tkl",centMethod.Data(),period[ip].Data());
          TString inMuFileName =Form("dphi/dphi_corr_%s_%s_%s.root","mu_tkl",centMethod.Data(),period[ip].Data());
          TString outTkFileName=Form("v2/v2_%s_%s_%s%s.root","tkl_tkl",centMethod.Data(),period[ip].Data(),systCheck[is].Data());//for the syst on the DetaRegion we want to use the default reference
          TString outMuFileName=Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[ip].Data(),systCheck[is].Data());
          // tracklet v2
          if(gStudySystematic==k_vertex05 || gStudySystematic==k_vertex01 || gStudySystematic==k_mixed_norm || gStudySystematic==k_mixed_integrated) {
            inTkFileName.ReplaceAll(".root",Form("%s.root",systCheck[is].Data()));
            inMuFileName.ReplaceAll(".root",Form("%s.root",systCheck[is].Data()));
          }
          nbinst = isTkl ? nbins_tkl : nbins_trk;
          binst  = isTkl ? bins_tkl  : bins_trk;
          for (Int_t ibin=0;ibin<nbinst;ibin++) binst_centers[ibin]=(binst[ibin]+binst[ibin+1])/2.;
          extractV2(inTkFileName,outTkFileName,nc,nbinst,nbinsa,binst_centers,0);

          nbinst = isTkl ? nbins_muon_tkl : nbins_muon_trk;
          binst  = isTkl ? bins_muon_tkl  : bins_muon_trk;
          if(isAMPTgen){
            nbinst = nbins_muon_ampt_gen;
            binst  = bins_muon_ampt_gen;
          }
          for (Int_t ibin=0;ibin<nbinst;ibin++) binst_centers[ibin]=(binst[ibin]+binst[ibin+1])/2.;
          extractV2(inMuFileName,outMuFileName,nc,nbinst,nbinsa,binst_centers,outTkFileName.Data());
        }//syst check
      }

      if(CalculateSystErrorOnV2){//draw syst for mu_tkl

        TCanvas *c=new TCanvas(Form("csyst_%s_%s",centMethod.Data(),period[ip].Data()),Form("csyst_%s_%s",centMethod.Data(),period[ip].Data()),600,1200);
        c->Divide(1,nc);
        TLegend *leg=new TLegend(0.01,0.01,0.99,0.99);
        //output file name
        TString outFileNameSyst = Form("syst/syst_v2_%s_%s_%s.root","mu_tkl",centMethod.Data(),period[ip].Data());

        //create graphs
        TGraphErrors*** g = new TGraphErrors**[nsyst];
        for (Int_t is=0; is<nsyst; is++) g[is] = new TGraphErrors*[nc];

        //get graphs
        for(Int_t is=0;is<nsyst;is++){
          TFile * f=new TFile(Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[ip].Data(),systCheck[is].Data()));
          for(Int_t ic=0;ic<nc;ic++)g[is][ic]=(TGraphErrors*)f->Get(Form("gv2_cent%d_ass%d",ic,assBin))->Clone();
          f->Close();
          delete f;
        }

        DrawSystError(g,c,leg,centMethod,outFileNameSyst,period[ip].Data(),"pol1");
        //        DrawSystError(g,c,leg,centMethod,outFileNameSyst,period[ip].Data(),"");

        TCanvas* cleg = new TCanvas("cleg","cleg");
        leg->Draw();
      }
      //draw syst on v2
    }//period


    if(CalculateSystErrorOnRatio){//draw syst for mu_tkl
      TCanvas *c=new TCanvas(Form("csystRatio_%s",centMethod.Data()),Form("csystRatio_%s",centMethod.Data()),600,1200);
      c->Divide(1,nc);
      TLegend *leg=new TLegend(0.2,0.2,0.8,0.8);
      //output file name
      TString outFileNameSyst=Form("syst/syst_v2_%s_%s_%s.root","mu_tkl",centMethod.Data(),"Ratio");

      //p going direction
      //create graphs
      TGraphErrors*** g1 = new TGraphErrors**[nsyst];
      for (Int_t is=0; is<nsyst; is++) g1[is] = new TGraphErrors*[nc];
      //get graphs
      for(Int_t is=0;is<nsyst;is++){
        TFile* f = new TFile(Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[0].Data(),systCheck[is].Data()));
        for(Int_t ic=0;ic<nc;ic++) g1[is][ic]=(TGraphErrors*)f->Get(Form("gv2_cent%d_ass%d",ic,assBin))->Clone();
        f->Close();
        delete f;
      }

      //Pb going direction
      //create graphs
      TGraphErrors*** g2 = new TGraphErrors**[nsyst];
      for (Int_t is=0; is<nsyst; is++) g2[is] = new TGraphErrors*[nc];
      //get graphs
      for(Int_t is=0;is<nsyst;is++){
        TFile* f = new TFile(Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[1].Data(),systCheck[is].Data()));
        for(Int_t ic=0;ic<nc;ic++)g2[is][ic]=(TGraphErrors*)f->Get(Form("gv2_cent%d_ass%d",ic,assBin))->Clone();
        f->Close();
        delete f;
      }

      //calculate ratio between periods
      for(Int_t is=0;is<nsyst;is++)
        for(Int_t ic=0;ic<nc;ic++){
          Double_t x1,y1,ey1;
          Double_t x2,y2,ey2;
          for(Int_t it=0;it<g1[is][ic]->GetN();it++){
            g1[is][ic]->GetPoint(it,x1,y1);
            ey1=g1[is][ic]->GetErrorY(it);
            g2[is][ic]->Print();
            g2[is][ic]->GetPoint(it,x2,y2);
            ey2=g2[is][ic]->GetErrorY(it);
            //set the ratio
            g1[is][ic]->SetPoint(it,x2,y2/y1);
            g1[is][ic]->SetPointError(it,0,y2/y1 *(TMath::Sqrt(TMath::Abs(ey2*ey2/y2/y2-ey1*ey1/y1/y1))));//difference in quadrature of stat errors
          }
        }
      DrawSystError(g1,c,leg,centMethod,outFileNameSyst,"ratio","pol1");
      //      DrawSystError(g1,c,leg,centMethod,outFileNameSyst,"ratio","");

      TCanvas *c=new TCanvas("cleg","cleg");
      leg->Draw();
    }

    if(CalculateSystErrorOnDiff){//draw syst for mu_tkl
      TCanvas *c=new TCanvas(Form("csystDiff_%s",centMethod.Data()),Form("csystDiff_%s",centMethod.Data()),600,1200);
      c->Divide(1,nc);
      TLegend *leg=new TLegend(0.2,0.2,0.8,0.8);
      //output file name
      TString outFileNameSyst=Form("syst/syst_v2_%s_%s_%s.root","mu_tkl",centMethod.Data(),"Diff");

      //p going direction
      //create graphs
      TGraphErrors*** g1 = new TGraphErrors**[nsyst];
      for (Int_t is=0; is<nsyst; is++) g1[is] = new TGraphErrors*[nc];
      //get graphs
      for(Int_t is=0;is<nsyst;is++){
        TFile* f = new TFile(Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[0].Data(),systCheck[is].Data()));
        for(Int_t ic=0;ic<nc;ic++) g1[is][ic]=(TGraphErrors*)f->Get(Form("gv2_cent%d_ass%d",ic,assBin))->Clone();
        f->Close();
        delete f;
      }

      //Pb going direction
      //create graphs
      TGraphErrors*** g2 = new TGraphErrors**[nsyst];
      for (Int_t is=0; is<nsyst; is++) g2[is] = new TGraphErrors*[nc];
      //get graphs
      for(Int_t is=0;is<nsyst;is++){
        TFile* f = new TFile(Form("v2/v2_%s_%s_%s%s.root","mu_tkl",centMethod.Data(),period[1].Data(),systCheck[is].Data()));
        for(Int_t ic=0;ic<nc;ic++)g2[is][ic]=(TGraphErrors*)f->Get(Form("gv2_cent%d_ass%d",ic,assBin))->Clone();
        f->Close();
        delete f;
      }

      //calculate ratio between periods
      for(Int_t is=0;is<nsyst;is++)
        for(Int_t ic=0;ic<nc;ic++){
          Double_t x1,y1,ey1;
          Double_t x2,y2,ey2;
          for(Int_t it=0;it<g1[is][ic]->GetN();it++){
            g1[is][ic]->GetPoint(it,x1,y1);
            ey1=g1[is][ic]->GetErrorY(it);
            g2[is][ic]->Print();
            g2[is][ic]->GetPoint(it,x2,y2);
            ey2=g2[is][ic]->GetErrorY(it);
            //set the ratio
            g1[is][ic]->SetPoint(it,x2,y2-y1);
            g1[is][ic]->SetPointError(it,0,sqrt(ey1*ey1+ey2*ey2));
          }
        }
      DrawSystError(g1,c,leg,centMethod,outFileNameSyst,"diff","pol1");

      TCanvas *c=new TCanvas("cleg","cleg");
      leg->Draw();
    }

  } // multiplicity estimator
}

void DrawSystError(TGraphErrors*** g,TCanvas *c,TLegend *leg,TString centMethod,TString outFileNameSyst,TString period, TString fitFunction=""){
  Int_t DrawFirst=-1;
  for(Int_t is=0;is<nsyst;is++){ // syst check
    Printf("\033[1;31m %s  period %s \033[m",systCheck[is].Data(),period.Data());
    //    remove syst not considered
    //    if(is==k_mixed_norm || is==k_parabolic_fit  || is==k_fit_opt  || is==k_mixed_integrated)continue;
    //    if(is!=kDefault && is!=k_vertex05 && is!=k_vertex01 && is!=k_60_to_70)continue;
    //    if(is!=kDefault && is!=k_hist_integral && is!=k_constant_fit)continue;
    //        if(is==k_mixed_norm || is==k_hist_integral || is==k_parabolic_fit || is==k_fit_opt || is==k_v2_analytical || is==k_no_subtraction || is==k_vertex01)continue;
    DrawFirst++;
    for(Int_t ic=0;ic<nc;ic++){
      Printf("\033[1;31m cent: %d \033[m",ic);
      c->cd(1);
      gPad->SetGridy();
      g[is][ic]->SetMarkerSize(1.2);
      g[is][ic]->SetMarkerStyle(marker[is]);
      g[is][ic]->SetMarkerColor(color[is]);
      g[is][ic]->SetLineColor(color[is]);
      g[is][ic]->SetLineWidth(1.2);
      if(ic==0)leg->AddEntry(g[is][ic],systCheck[is].Data(),"lp");
      g[is][ic]->GetYaxis()->SetTitle("v_{2,#mul}");
      g[is][ic]->GetXaxis()->SetTitle("#it{p}_{T}");
      g[is][ic]->GetYaxis()->SetTitleSize(0.07);
      g[is][ic]->GetYaxis()->SetTitleOffset(0.6);
      g[is][ic]->GetYaxis()->SetLabelSize(0.05);
      g[is][ic]->GetXaxis()->SetTitleSize(0.07);
      g[is][ic]->GetXaxis()->SetTitleOffset(0.6);
      g[is][ic]->GetXaxis()->SetLabelSize(0.05);

      if(ic==0)g[is][ic]->Clone()->Draw((is==0)?"ap":"psame");
      if(is==0)continue;//no ratio for the default case
      TGraphErrors* gratio=(TGraphErrors*)g[is][ic]->Clone("gratio");
      gratio->SetName(Form("gUnc_%s%s_Cent%d",centMethod.Data(),systCheck[is].Data(),ic));
      TGraphErrors* gdelta=(TGraphErrors*)g[is][ic]->Clone("gdelta");
      for(Int_t ipoint=0;ipoint<g[is][ic]->GetN();ipoint++){
        Double_t x,y,ey;
        Double_t xref,yref,eyref;
        gratio->GetPoint(ipoint,x,y);
        ey=gratio->GetErrorY(ipoint);
        g[0][ic]->GetPoint(ipoint,xref,yref);
        eyref=g[0][ic]->GetErrorY(ipoint);
        //set the ratio
        //        gratio->SetPoint(ipoint,x+0.05*is,y/yref);
        gratio->SetPoint(ipoint,x,y/yref);
        gratio->SetPointError(ipoint,0,y/yref *(TMath::Sqrt(TMath::Abs(ey*ey/y/y-eyref*eyref/yref/yref))));//difference in quadrature of stat errors
        //set the Delta/sigma
        Double_t err=TMath::Sqrt(TMath::Abs(ey*ey-eyref*eyref));//difference in quadrature as worst case, as in arXiv:hep-ex/0207026v1
        //        gdelta->SetPoint(ipoint,x+0.05*is,TMath::Abs(y-yref)/err);
        gdelta->SetPoint(ipoint,x,TMath::Abs(y-yref)/err);
        gdelta->SetPointError(ipoint,0,0);
      }
      c->cd(2);
      gPad->SetGridy();
      gratio->GetYaxis()->SetTitle("v_{2,#mu}/v_{2,#mu,def}");
      gratio->GetXaxis()->SetTitle("#it{p}_{T}");
      gratio->GetYaxis()->SetRangeUser(0.5,1.5);
      gratio->GetYaxis()->SetTitleSize(0.07);
      gratio->GetYaxis()->SetTitleOffset(0.6);
      gratio->GetYaxis()->SetLabelSize(0.05);
      gratio->GetXaxis()->SetTitleSize(0.07);
      gratio->GetXaxis()->SetTitleOffset(0.6);
      gratio->GetXaxis()->SetLabelSize(0.05);
      if(ic==0)gratio->Clone()->Draw((DrawFirst==1)?"ap":"psame");

      TF1* fitfunc=0x0;
      if(fitFunction==""){//fit function not defined, selected according to the smaller chi2/ndf
        TF1* pol0=new TF1(Form("Unc_%s%s_Cent%d",centMethod.Data(),systCheck[is].Data(),ic),"pol0",0.5,4.);
        TF1* pol1=new TF1(Form("Unc_%s%s_Cent%d",centMethod.Data(),systCheck[is].Data(),ic),"pol1",0.5,4.);
        TF1* pol2=new TF1(Form("Unc_%s%s_Cent%d",centMethod.Data(),systCheck[is].Data(),ic),"pol2",0.5,4.);
        gratio->Fit(pol0,"","MN",0.50001,4);
        gratio->Fit(pol1,"","MN",0.50001,4);
        gratio->Fit(pol2,"","MN",0.50001,4);
        Double_t chi2ndf0=pol0->GetChisquare()/pol0->GetNDF();
        Double_t chi2ndf1=pol1->GetChisquare()/pol1->GetNDF();
        Double_t chi2ndf2=pol2->GetChisquare()/pol2->GetNDF();
        Printf("pol0: chi2: %f  ndf: %d  chi2/ndf: %f",pol0->GetChisquare(),pol0->GetNDF(),chi2ndf0);
        Printf("pol1: chi2: %f  ndf: %d  chi2/ndf: %f",pol1->GetChisquare(),pol1->GetNDF(),chi2ndf1);
        Printf("pol2: chi2: %f  ndf: %d  chi2/ndf: %f",pol2->GetChisquare(),pol2->GetNDF(),chi2ndf2);
        if(chi2ndf0<=chi2ndf1 && chi2ndf0<=chi2ndf2){
          fitfunc=pol0->Clone();
        }else if(chi2ndf1<=chi2ndf2){
          fitfunc=pol1->Clone();
        }else{
          fitfunc=pol2->Clone();
        }
        delete pol0;
        delete pol1;
        delete pol2;
      } else fitfunc=new TF1(Form("Unc_%s%s_Cent%d",centMethod.Data(),systCheck[is].Data(),ic),fitFunction.Data(),0.5,4.);
      gratio->Fit(fitfunc,"","MN",0.50001,4);
      fitfunc->SetLineColor(gratio->GetLineColor());
      if(ic==0)fitfunc->DrawCopy("same");
      fitfunc->SetParameter(0,fitfunc->GetParameter(0)-1.);
      //save TF1 with errors
      TFile * fsys=new TFile(outFileNameSyst,(is==1 && ic==0)?"recreate":"update");
      fitfunc->Write();
      //transform the ratio in order to use directly in the syst error combination
      for(Int_t ipoint=0;ipoint<gratio->GetN();ipoint++){
        // remove non-significant points
        // if(gdelta->GetY()[ipoint]<1.)gratio->SetPoint(ipoint,gratio->GetX()[ipoint],0.);
        // else gratio->SetPoint(ipoint,gratio->GetX()[ipoint],gratio->GetY()[ipoint]-1);
        gratio->SetPoint(ipoint,gratio->GetX()[ipoint],gratio->GetY()[ipoint]-1);
      }
      gratio->Write();
      fsys->Close();
      //draw significance of the check
      c->cd(3);
      gPad->SetGridy();
      gdelta->GetYaxis()->SetRangeUser(0,5);
      gdelta->GetYaxis()->SetTitle("#Deltav_{2,#mu}/#sigma_{#Deltav_{2,#mu}}");
      gdelta->GetXaxis()->SetTitle("#it{p}_{T}");
      gdelta->GetYaxis()->SetTitleSize(0.07);
      gdelta->GetYaxis()->SetTitleOffset(0.6);
      gdelta->GetYaxis()->SetLabelSize(0.05);
      gdelta->GetXaxis()->SetTitleSize(0.07);
      gdelta->GetXaxis()->SetTitleOffset(0.6);
      gdelta->GetXaxis()->SetLabelSize(0.05);
      if(ic==0)gdelta->Clone()->Draw((DrawFirst==1)?"alp":"plsame");
      if(is==1){
        TLine *l =new TLine(0,1,4.,1);
        if(ic==0)l->Draw();
      }
      delete fsys;
      delete fitfunc;
      delete gratio;
      delete gdelta;
    }//mult class
  }//nsyst
}
