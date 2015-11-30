#include "subtraction.C"

void extractV2(TString inFileName, TString outFileName, const Int_t nc, const Int_t nt, const Int_t na, Double_t* pt, Char_t* trackletFileName=0){
  Double_t res[14];

  Double_t V2val[nc][na][nt];
  Double_t V2err[nc][na][nt];
  Double_t fval[nc][na][nt];
  Double_t ferr[nc][na][nt];
  Double_t sval1[nc][na][nt];
  Double_t sval2[nc][na][nt];

  TString outFileNameFactors=outFileName.Data();
  outFileNameFactors.ReplaceAll("/v2","/Factors");
  TFile *fFactors=0x0;

  for (Int_t ic=0;ic<nc;ic++){
    for (Int_t it=0;it<nt;it++) {
      for (Int_t ia=0;ia<na;ia++) {
        printf("ic=%i it=%i ia=%i n=%p\n",ic,it,ia,trackletFileName);
        if (!trackletFileName){
          Double_t ExclMax=1.2;
          if(gStudySystematic==k_exclusion10) ExclMax=1.0;
          if(gStudySystematic==k_exclusion08) ExclMax=0.8;
          subtraction(inFileName,it,ia,ic,(gStudySystematic==k_60_to_70)?4:3,-2.,2.,-ExclMax,ExclMax,0,res); //tkl-tkl
        }else {
          subtraction(inFileName,it,ia,ic,(gStudySystematic==k_60_to_70)?4:3,(gStudySystematic==k_hist_integral)?-4.5:-5.,(gStudySystematic==k_hist_integral)?-2:-1.5,-10,-10,0,res); //mu-tkl
        }
        V2val[ic][ia][it]=res[2];
        V2err[ic][ia][it]=res[6];
        fval[ic][ia][it]=res[10];
        ferr[ic][ia][it]=res[11];
        sval1[ic][ia][it]=res[12];
        sval2[ic][ia][it]=res[13];
      }
    }
  }

  // diagonal (calculated before non-diagonal)
  Double_t v2v[nc][na];
  Double_t v2e[nc][na];
  TGraphErrors* gv2[nc];
  // full matrix ia vs it (calculated with diagonal values)
  Double_t v2val[nc][na][nt];
  Double_t v2err[nc][na][nt];
  TGraphErrors* gv2ass[nc][na];

  if (trackletFileName) {
    // muon case: we take v2v[ic][ia] - reference v2 from track(lets) we want to devide out
    TFile* fv2 = new TFile(trackletFileName);
    for (Int_t ic=0;ic<nc;ic++){
      gv2[ic] = (TGraphErrors*) fv2->Get(Form("gv2_cent%i",ic));//in the mu-tkl file we have the gv2[] of the tkl-tkl
      for (Int_t ia=0;ia<na;ia++){
        v2v[ic][ia] = gv2[ic]->GetY()[ia];
        v2e[ic][ia] = gv2[ic]->GetEY()[ia];
      }
    }
    fv2->Close();
    delete fv2;
  }
  else {
    // track(let) case: we calculate diagonal v2v[ic][ia] with symmetric bins
    for (Int_t ic=0;ic<nc;ic++){
      for (Int_t ia=0;ia<na;ia++){
        v2v[ic][ia]=sqrt(V2val[ic][ia][ia]);
        v2e[ic][ia]=0.5*v2v[ic][ia]*(V2err[ic][ia][ia]/V2val[ic][ia][ia]);
      }
      gv2[ic] = new TGraphErrors(na,pt,v2v[ic],0,v2e[ic]);//gv2[] contains symmetric bins in dphi
      gv2[ic]->SetName(Form("gv2_cent%i",ic));
    }
  }


  for (Int_t ic=0;ic<nc;ic++){
    for (Int_t ia=0;ia<na;ia++){
      for (Int_t it=0;it<nt;it++){
        v2val[ic][ia][it]=V2val[ic][ia][it]/v2v[ic][ia];//v2v[ic][ia] calculated with symmetric bins
        v2err[ic][ia][it]=v2val[ic][ia][it]*sqrt(pow(V2err[ic][ia][it]/V2val[ic][ia][it],2)+pow(v2e[ic][ia]/v2v[ic][ia],2));
      }
      gv2ass[ic][ia] = new TGraphErrors(nt,pt,v2val[ic][ia],0,v2err[ic][ia]);
      gv2ass[ic][ia]->SetName(Form("gv2_cent%i_ass%i",ic,ia));
      gv2ass[ic][ia]->SetLineColor(color[ia]);
      gv2ass[ic][ia]->SetMarkerStyle(marker[ia]);
      gv2ass[ic][ia]->SetMarkerColor(color[ia]);
    }
  }

  if(gStudySystematic==k_LowMultScale){//saving scaling factors in file
    TGraphErrors* gfass[nc][na];
    TGraphErrors* gfassS[nc][na];
    TGraphErrors* gs1[nc][na];
    TGraphErrors* gs2[nc][na];
    for (Int_t ic=0;ic<nc;ic++){
      for (Int_t ia=0;ia<na;ia++){
        gfass[ic][ia] = new TGraphErrors(nt,pt,fval[ic][ia],0,ferr[ic][ia]);
        gfass[ic][ia]->SetName(Form("gf_cent%i_ass%i",ic,ia));
        gfass[ic][ia]->SetLineColor(color[ia]);
        gfass[ic][ia]->SetMarkerStyle(marker[ia]);
        gfass[ic][ia]->SetMarkerColor(color[ia]);
        //fit and save the scaling factors smoothed
        gfassS[ic][ia] = (TGraphErrors*)gfass[ic][ia]->Clone(Form("%s_smoothed",gfass[ic][ia]->GetName()));
        //gfassS[ic][ia]->SetName(Form("gf_cent%i_ass%i_smoothed",ic,ia));
        TF1 *pol1=new TF1("pol1","pol1",0,5);
        gfassS[ic][ia]->Fit(pol1,"0","N");
        for(Int_t it=0;it<gfassS[ic][ia]->GetN();it++){
          gfassS[ic][ia]->SetPoint(it,gfassS[ic][ia]->GetX()[it],pol1->Eval(gfassS[ic][ia]->GetX()[it]));
        }
        gfassS[ic][ia]->SetLineColor(color[ia]);
        gfassS[ic][ia]->SetMarkerStyle(marker[ia]);
        gfassS[ic][ia]->SetMarkerColor(color[ia]);
        //save sigma
        gs1[ic][ia] = new TGraphErrors(nt,pt,sval1[ic][ia],0,0);
        gs1[ic][ia]->SetName(Form("gsHM_cent%i_ass%i",ic,ia));
        gs1[ic][ia]->SetLineColor(color[ia]);
        gs1[ic][ia]->SetMarkerStyle(marker[ia]);
        gs1[ic][ia]->SetMarkerColor(color[ia]);
        gs2[ic][ia] = new TGraphErrors(nt,pt,sval2[ic][ia],0,0);
        gs2[ic][ia]->SetName(Form("gsLM_cent%i_ass%i",ic,ia));
        gs2[ic][ia]->SetLineColor(color[ia]);
        gs2[ic][ia]->SetMarkerStyle(marker[ia]);
        gs2[ic][ia]->SetMarkerColor(color[ia]);

      }
    }
    TFile* fout = new TFile(outFileNameFactors,"recreate");
    for (Int_t ic=0;ic<1;ic++) {//save only 0-20%
//    for (Int_t ic=0;ic<nc;ic++) {
      for (Int_t ia=0;ia<na;ia++) {
        gfass[ic][ia]->Write();
        gfassS[ic][ia]->Write();
        gs1[ic][ia]->Write();
        gs2[ic][ia]->Write();

      }
    }
    fout->Close();
    delete fout;
  }

  //add trigger stat file to the output
  TFile *infile=new TFile(inFileName);
  TFile* fout = new TFile(outFileName,"recreate");
  for (Int_t ic=0;ic<nc;ic++) {
    gv2[ic]->Write();
    TH1D* hTriggerStat = (TH1D*) infile->Get(Form("TriggerStat_%d",ic));
    if(hTriggerStat)hTriggerStat->Write();
    for (Int_t ia=0;ia<na;ia++) {
      gv2ass[ic][ia]->Write();
    }
  }
  infile->Clone();
  delete infile;
  fout->Close();
  delete fout;
}
