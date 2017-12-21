#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#endif

// Macro for automatized checks on the output of AliMultSelectionTask
//    -> aimed at fast checks on the output of the AnalysisQA trains
//
// 1) checks that the QA histograms are filled for the proper estimators
//    NOTE: the enabled estimators are different in pp, p-Pb and Pb-Pb
// 2) check that the percentile distributions are flat
//    NOTE: the range in which the percentile distribution has to flat
//          is different for  pp, p-Pb and Pb-Pb
//



void ProcessOutputMultSelTask(TString filename="AnalysisResults.root", TString collsyst="pp"){

  TFile* fil=new TFile(filename.Data());
  if(!fil || !fil->IsOpen()){
    printf("ERROR: file %s not open\n",filename.Data());
    return;
  }
  TDirectoryFile* mso=(TDirectoryFile*)fil->Get("MultSelection");
  Bool_t isMC=kFALSE;
  if(!mso){
    mso=(TDirectoryFile*)fil->Get("MultSelection_MC");
    isMC=kTRUE;
  }
  if(!mso){
    printf("ERROR: TDirectoryFile MultSelection not present in %s\n",filename.Data());
    return;
  }
  TList* msl=(TList*)mso->Get("cListMultSelection");
  if(!msl){
    printf("TList cListMultSelection not present\n");
    return;
  }

  // Extract collsyst from histogram title
  TH1D* hev=(TH1D*)msl->FindObject("fHistEventCounter");
  Bool_t readFromHisto=kFALSE;
  if(hev){
    TString htit=hev->GetTitle();
    cout<<htit.Data()<<endl;
    if(htit.Contains("Event type: PbPb")){
      collsyst="Pb-Pb";
      readFromHisto=kTRUE;
    }else if(htit.Contains("Event type: XeXe")){
      collsyst="Xe-Xe";
      readFromHisto=kTRUE;
    }else if(htit.Contains("Event type: pA")){
      collsyst="p-Pb";
      readFromHisto=kTRUE;
    }else if(htit.Contains("Event type: pp")){
      collsyst="pp";
      readFromHisto=kTRUE;
    }
  }


  TString estimNames[11]={"V0M","V0A","V0C","CL0","CL1","SPDClusters","SPDTracklets","ZNA","ZNC","ZNApp","ZNCpp"};
  TCanvas* ce=new TCanvas("ce","Estimators",1200,800);
  ce->Divide(4,3);
  ce->cd(1);
  TString cmc="Monte Carlo";
  if(!isMC) cmc="Data";
  TLatex* tt1=new TLatex(0.1,0.76,cmc.Data());
  tt1->Draw();
  tt1->SetTextFont(43);
  tt1->SetTextSize(22);
  TString ccs=Form("Collision system: %s",collsyst.Data());
  TLatex* tt2=new TLatex(0.1,0.6,ccs.Data());
  tt2->SetTextFont(43);
  tt2->SetTextSize(22);
  tt2->Draw();
  TString chow="  (from macro argument)";
  if(readFromHisto) chow="  (from fHistEventCounter title)";
  TLatex* tt3=new TLatex(0.1,0.5,chow.Data());
  tt3->SetTextFont(43);
  tt3->SetTextSize(18);
  tt3->Draw();

  for(Int_t je=0; je<11; je++){
    Bool_t shouldBeFilled=kFALSE;
    // maxCent is the upper limit of the percentile interval in which the distribution should be flat
    // depends on collision system and estimator
    Double_t maxCent=100.;
    if(collsyst=="pp"){
      // enables the estimators which are calibrated for pp
      if(estimNames[je]=="V0M" || estimNames[je]=="V0A" || estimNames[je]=="V0C" || 
	 estimNames[je]=="ZNApp" ||  estimNames[je]=="ZNCpp" || 
	 estimNames[je]=="SPDClusters" || estimNames[je]=="SPDTracklets") shouldBeFilled=kTRUE;
    }else if(collsyst=="Pb-Pb" || collsyst=="Xe-Xe"){
      if(estimNames[je]=="V0M" || estimNames[je]=="SPDTracklets" ||
	 estimNames[je]=="CL0" || estimNames[je]=="CL1") shouldBeFilled=kTRUE;
      maxCent=90.;
    }else{
      // enables the estimators which are calibrated for p-Pb 
      if(estimNames[je]=="V0M" || estimNames[je]=="V0A" || estimNames[je]=="V0C" || 
	 estimNames[je]=="ZNA" ||  estimNames[je]=="ZNC" || 
	 estimNames[je]=="CL0" || estimNames[je]=="CL1") shouldBeFilled=kTRUE;
      if(collsyst=="p-Pb" && (estimNames[je]=="ZNA" ||  estimNames[je]=="ZNC")) maxCent=95.;
    }
    TString histoname=Form("fHistQASelected_%s",estimNames[je].Data());
    printf("--- Retrieving histogram %s ---\n",histoname.Data());
    TH1D* h=(TH1D*)msl->FindObject(histoname.Data());
    if(h){
      ce->cd(je+2);
      h->SetLineWidth(2);
      h->Draw("e");
      Int_t ib0=h->FindBin(0.01);
      Int_t ibMax=h->FindBin(maxCent-0.1);
      Double_t entries0Max=h->Integral(ib0,ibMax);
      printf("Estimator %s  events in 0-%d = %.0f\n",estimNames[je].Data(),TMath::Nint(maxCent),entries0Max);
      TString check="OK";
      if(shouldBeFilled && entries0Max<0.5) check="BAD: should fe filled";
      if(!shouldBeFilled && entries0Max>0.5) check="BAD: should be empty";
      if(entries0Max>300){
	h->Fit("pol0","Q","",0.,maxCent);
	TF1* ffit=(TF1*)h->GetListOfFunctions()->FindObject("pol0");
	if(ffit){
	  ffit->SetLineWidth(2);
	  ffit->SetLineColor(kGray+1);
	  Double_t mfit=ffit->GetParameter(0);
	  Double_t ifflat=entries0Max*h->GetBinWidth(1)/maxCent;
	  Double_t emfit=ffit->GetParError(0);
	  Double_t delta=TMath::Abs(mfit-entries0Max*h->GetBinWidth(1)/maxCent);
	  if(TMath::Abs(mfit-ifflat)>3.*emfit) check="BAD: not flat";
	}
	if(check=="OK"){
	  Double_t devsig=0.;
	  Double_t chi2=0.;
	  Int_t npt=0;
	  for(Int_t ib=ib0; ib<=ibMax; ib++){
	    Double_t bw=h->GetBinWidth(ib);
	    Double_t bc=h->GetBinContent(ib);
	    Double_t ebc=TMath::Sqrt(bc);
	    Double_t exp=entries0Max*bw/maxCent;
	    if(bc==0) ebc=1;
	    devsig+=(bc-exp)/ebc;
	    chi2+=devsig*devsig;
	    npt++;
	  }
	  if(TMath::Abs(devsig/npt)>3. || TMath::Sqrt(chi2)/npt>3.) check="BAD: not flat";
	}
      }
      TLatex* tes=new TLatex(0.17,0.78,estimNames[je].Data());
      tes->SetNDC();
      tes->SetTextFont(43);
      tes->SetTextSize(22);
      tes->Draw();
      TLatex* tqa=new TLatex(0.17,0.70,check.Data());
      tqa->SetNDC();
      if(check=="OK") tqa->SetTextColor(kGreen+1);
      else tqa->SetTextColor(2);
      tqa->SetTextFont(63);
      tqa->SetTextSize(22);
      tqa->Draw();
    }
  }



}
