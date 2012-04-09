TF1 *CreateFuncGausTail(TH1F *h,TString funcname);
TF1 *CreateFuncGaussTail(TH1F *h,TString funcname,Double_t parone,Double_t partwo);

void PlotImpParResVsPt()
{
  gStyle->SetOptStat(0);


  TFile *file1=TFile::Open("QAresults_LHC11h_p2_170593.root","READ");
  TDirectoryFile *dir1=(TDirectoryFile*)file1->GetDirectory("ImpParRes_Performance");
  
  //TFile *d0dist1=TFile::Open("/home/luojb/2011/ImpactParameter/OutputFile/PP/LHC10b/ImpParRes.Performance_LHC10bcde_pass2_ESDdiamond_2011_01_26.root","READ");
  
  //if (!d0dist1->IsOpen()) return;
  
  

  TList *d0allpointRec1 = (TList*)dir1->Get("coutputd0allPointRec_0_1000000");
   
  TList *d0allpointSkip1 = (TList*)dir1->Get("coutputd0allPointSkip_0_1000000");
  
  TList *d0particlePID1 = (TList*)dir1->Get("coutputd0particlePID_0_1000000");
   
  TList *d0pt = (TList*)dir1->Get("coutputd0Pt_0_1000000");
  //TFile *d0ppr=TFile::Open("../ResolutionsAnalysis_def.root","READ"); 
  //if(!d0ppr->IsOpen()) return;
 
  //define the histgram
  
  TH1F **d0AllpointrphiSkip1_=new TH1F*[26]; 
  
  TH1F **d0AllpointrphiSkipTail1_=new TH1F*[26]; 
  TH1F **d0AllpointrphiSkipGaus1_=new TH1F*[26];  
 
  TH1F **d0Pt_=new TH1F*[26];
  Int_t nbins = 51;  
  Double_t xbins[52]={0.22,0.23,0.26,0.27,0.35,0.36,0.45,0.46,0.55,0.56,0.65,0.66,0.75,0.76,0.85,0.86,1.05,1.06,1.25,1.27,1.45,1.47,1.65,1.67,1.85,1.87,2.15,2.17,2.45,2.47,2.65,2.67,2.85,2.87,3.25,3.27,3.75,3.77,4.15,4.20,4.95,5.15,5.35,5.55,6.0,6.8,8.5,10.5,12,19.,21.,32};
 
  //TGraph *Resrphi_allpoint =new TGraph("Resrphi_allpot","d_0 impact parameter resolution_rphi [#mu m];    p_{t} [GeV/c]",nbins,xbins); 
  TGraphErrors *Resrphi_allpointRec =new TGraphErrors(26);
  //Resrphi_allpointRec->SetMaximum(400); 
  //Resrphi_allpointRec->SetMinimum(0); 
  //Resrphi_allpoint->SetXTitle("p_{t} [GeV/c]");
  //Resrphi_Rec->SetYTitle("d^{0} impact parameter resolution_rphi [#mu m]");
  Resrphi_allpointRec->SetLineWidth(1);
  Resrphi_allpointRec->SetMarkerStyle(20);  
  
  TGraphErrors *Resrphi_allpointSkip1 =(TGraphErrors*)Resrphi_allpointRec->Clone("Resrphi_allpointSkip1");
  
  TGraphErrors *Resrphi_allpointSkipmean1 =(TGraphErrors*)Resrphi_allpointRec->Clone("Resrphi_allpointSkipmean1");
   
 
  // fit gaussian
  /*
  TF1 *allpointRec1 = new TF1("allpointRec1","gaus",-10000.,10000.);
  TF1 *allpointRec2 = new TF1("allpointRec2","gaus",-10000.,10000.);
  TF1 *allpointRec3 = new TF1("allpointRec3","gaus",-10000.,10000.);
  TF1 *allpointRec4 = new TF1("allpointRec4","gaus",-10000.,10000.);
  TF1 *allpointRec5 = new TF1("allpointRec5","gaus",-10000.,10000.);
  TF1 *allpointRec6 = new TF1("allpointRec6","gaus",-10000.,10000.);
 */

  
  TF1 *f1 = new TF1("f1","gaus",-10000.,10000.);
  TF1 *wn1 = new TF1("wn1","gaus",-10000.,10000.);
  //fitting  histogram
  for (Int_t i=0;i<26;i++){
    TString getstr;
    Int_t jj=i+1;
    Double_t j =3.;  
    getstr="d0allpointrphiSkip_"; 
    getstr += jj;
    d0AllpointrphiSkip1_[i] = (TH1F*)d0allpointSkip1->FindObject(getstr);//->Clone(getstr+"d0allrphiRec");

    getstr="d0pt_";
    getstr += jj;
    d0Pt_[i] = (TH1F*)d0pt->FindObject(getstr);

   
    //Set fitting function
    //normalization
    TF1 *h = new TF1("h","gaus",-10000.,10000.);
    d0AllpointrphiSkip1_[i]->Fit("h","NR,Q");
    Double_t d0rphirange_allpointSkip1 = h->GetParameter(2);   
    Double_t Entries = d0AllpointrphiSkip1_[i]->GetEntries();
    Double_t Integral = d0AllpointrphiSkip1_[i]->Integral();
    d0AllpointrphiSkip1_[i]->Scale(1./Integral);

    //cp histogram
    d0AllpointrphiSkipTail1_[i] = new TH1F(*d0AllpointrphiSkip1_[i]);   // = new TH1F(*(TH1F*)d0AllpointrphiSkip1_[i]);               
    d0AllpointrphiSkipGaus1_[i] = new TH1F(*d0AllpointrphiSkip1_[i]);
    Double_t wholerangeint=d0AllpointrphiSkip1_[i]->Integral();
    Double_t cutleft1= -j*(d0rphirange_allpointSkip1);
    Double_t cutright1 =j*d0rphirange_allpointSkip1;
    d0AllpointrphiSkipTail1_[i]->Reset(0);
    d0AllpointrphiSkipGaus1_[i]->Reset(0);
    for (Int_t bin=1;bin<d0AllpointrphiSkip1_[i]->GetNbinsX();bin++){
      Int_t bincontent = d0AllpointrphiSkip1_[i]->GetBinCenter(bin);
      if(bincontent<cutleft1 || bincontent>cutright1) {
	d0AllpointrphiSkipTail1_[i]->SetBinContent(bin,d0AllpointrphiSkip1_[i]->GetBinContent(bin));
	d0AllpointrphiSkipTail1_[i]->SetBinError(bin,d0AllpointrphiSkip1_[i]->GetBinError(bin));
	d0AllpointrphiSkipGaus1_[i]->SetBinContent(bin,0.);  
	//This sentence is very important,otherwise we will get the information the data is empty when we fit it .
      }
      else if(bincontent>=cutleft1 && bincontent<=cutright1){
	d0AllpointrphiSkipTail1_[i]->SetBinContent(bin,0);
	d0AllpointrphiSkipGaus1_[i]->SetBinContent(bin,d0AllpointrphiSkip1_[i]->GetBinContent(bin));
      	d0AllpointrphiSkipGaus1_[i]->SetBinError(bin,d0AllpointrphiSkip1_[i]->GetBinError(bin));  
      }
      //d0AllpointrphiSkipGaus_[i]->SetBinContent(bin,d0AllpointrphiSkip1_[i]->GetBinContent(bin));
      //d0AllpointrphiSkipGaus_[i]->SetBinError(bin,d0AllpointrphiSkip1_[i]->GetBinError(bin));  
    }
    //return;
    TF1 *hh;
    hh =CreateFuncTail(d0AllpointrphiSkipTail1_[i],"hh");
    hh->SetLineColor(4);
    d0AllpointrphiSkipTail1_[i]->SetLineColor(5);
    d0AllpointrphiSkipTail1_[i]->Fit("hh","NR,Q","",-2500,2500);  
    Double_t Sigmatail_allpointSkip1 = hh->GetParameter(2);
    d0AllpointrphiSkipGaus1_[i]->Fit("h","NR,Q","",-2500,2500);
    Double_t Sigmagaus_allpointSkip1 = h->GetParameter(2);
    TF1 *allpointSkip1;
    allpointSkip1=CreateFuncGaussTail(d0AllpointrphiSkip1_[i],"allpointSkip1",Sigmagaus_allpointSkip1,Sigmatail_allpointSkip1);
   
   

    Double_t pt =d0Pt_[i]->GetMean();
    Double_t RMS=d0Pt_[i]->GetRMS();
    cout<<"i:"<<i<<endl;  
    cout<<"pt:"<<pt<<endl;

    //fill the histogram 
    //Double_t ptpion=ptt*ptt/TMath::Sqrt(ptt*ptt+mass*mass);   
    //allpointSkip1->SetRange(-j*d0rphirange_allpointSkip1,j*d0rphirange_allpointSkip1);
    d0AllpointrphiSkip1_[i]->Fit("allpointSkip1","NR,Q");
    Resrphi_allpointSkip1->SetPoint(jj,pt,allpointSkip1->GetParameter(2));
    Resrphi_allpointSkip1->SetPointError(jj,RMS,allpointSkip1->GetParError(2));
    Resrphi_allpointSkipmean1->SetPoint(jj,pt,allpointSkip1->GetParameter(1));
    Resrphi_allpointSkipmean1->SetPointError(jj,RMS,allpointSkip1->GetParError(1));
    
  }
 
  //********************************************************************************************** 
  //************************************  Plot figures  ******************************************
  //**********************************************************************************************

  TCanvas *d0ResrphiSkip1 = new TCanvas("d0ResrphiSkip1","d0ResrphiSkip1",0,0,540,420);       
  d0ResrphiSkip1->cd();
  gPad->SetLogx(); 
  TH2F *hFrame111 = new TH2F("hFrame","; p_{t} [GeV/c]; d0_resolution [#mum]",2,0.15,38,2,0,350);   
  hFrame111->Draw();
   
  Resrphi_allpointSkip1->SetMarkerStyle(20);
  Resrphi_allpointSkip1->SetMarkerColor(2);
  Resrphi_allpointSkip1->Draw("p"); 
  /*
  TLegend *legend = new TLegend(0.28,0.75,0.65,0.88);
  legend-> SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(Resrphi_allpointSkip1,"d0ResrphiSkip_LHCbcd","pl");
  //legend->AddEntry(Resrphi_allpointRec4,"d0ResrphiSkip_LHC11a10a","pl");
  //legend->AddEntry(Resrphi_allpointRec2,"d0ResrphiSkip_LHC11a10a_m","pl");
  legend->Draw();
  */
  d0ResrphiSkip1->Update();

  TCanvas *d0MeanrphiSkip1 = new TCanvas("d0MeanrphiSkip1","d0MeanrphiSkip1",0,0,540,420); 
  d0MeanrphiSkip1->cd();
  gPad->SetLogx(); 
  TH2F *hFrame1 = new TH2F("hFrame1","; p_{t} [GeV/c]; d0_mean [#mum]",2,0.15,38,2,-35,30);
  hFrame1->Draw();  
  Resrphi_allpointSkipmean1->SetMarkerStyle(20);
  Resrphi_allpointSkipmean1->SetMarkerColor(2);
  Resrphi_allpointSkipmean1->Draw("p");  
  //Resrphi_allpointSkipmean2->SetMarkerStyle(20);
  //Resrphi_allpointSkipmean2->SetMarkerColor(3);
  //Resrphi_allpointSkipmean2->Draw("p,sames"); 
 
  /*
  TLegend *legend = new TLegend(0.2,0.7,0.7,0.88);
  legend-> SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(Resrphi_allpointSkipmean1,"d0MeanrphiSkip_LHC11a_partpointSkip_pass2_with_SDD","pl");
  //legend->AddEntry(Resrphi_allpointSkipmean2,"d0MeanrphiSkip_LHC11a_partpointSkip_pass3_with_SDD","pl");
  //legend->AddEntry(Resrphi_allpointSkipmean3,"d0MeanrphiSkip_LHC10d_partpointSkip_diamond","pl");
  legend->Draw();
  */
  d0MeanrphiSkip1->Update();

  return;

}
//---------------------------------------------------------------------------
TF1 *CreateFuncTail(TH1F *hh,TString funcname,Double_t wholeRangeInt=-1.)
{
  /*  TF1 *f3 = new TF1("f3","gaus",300,10000.);
  cout<<"here:"<<2222222<<endl;
  hh->Fit("f3","NR,Q");
  Double_t value =f3->GetParameter(2);
  cout<<"value:"<<value<<endl;*/
  TF1 *tail=new TF1(funcname.Data(),"[0]*(1./(2.*[2])*TMath::Exp(-TMath::Abs(x-[1])/[2]))",-2500.,2500.);
  //tail->SetParLimits(0,0.,250.);//Set the first parameter [0] range
  tail->SetParLimits(1,-50,50.);//Set the second parameter [1] range
  tail->SetParLimits(2,0,15000.);
  Double_t binwidth=hh->GetBinWidth(10);
  Double_t integral=hh->Integral();
  if(wholeRangeInt!=-1.)tail->SetParLimits(0,(0.2)*wholeRangeInt*binwidth,(0.5)*wholeRangeInt*binwidth);
  Double_t RMS1=TMath::Abs(hh->GetRMS());
  Double_t firstvalue1=binwidth*integral;
  //cout<<"binwidth:"<<binwidth<<endl;
  //cout<<"integral:"<<integral<<endl;
  //cout<<"RMS:"<<RMS<<endl;
  //cout<<"firstvalue1:"<<firstvalue1<<endl;
  tail->SetParameters(1.,0,100.);//Set the initial value of parameter
  //  tail->SetParameters(1,0,RMS1);
  
  return tail;

} 
//----------------------------------------------------------------------------
TF1 *CreateFuncGausTail(TH1F *h,TString funcname)
{
  TF1 *gaustail=new TF1(funcname.Data(),"[0]*([4]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))+(1.-[4])/(2.*[3])*TMath::Exp(-TMath::Abs(x-[1])/[3]))",-2500.,2500.);
  gaustail->SetParLimits(0,0.,250000000.);//Set the first parameter [0] range
  gaustail->SetParLimits(1,-20,20.);//Set the second parameter [1] range
  gaustail->SetParLimits(2,0.,500.);
  gaustail->SetParLimits(3,0.,500.);
  gaustail->SetParLimits(4,0.,1.);
  Double_t binwidth=h->GetBinWidth(10);
  Double_t integral=h->Integral();
  Double_t RMS1=TMath::Abs(h->GetRMS());
  Double_t firstvalue=binwidth*integral;
  //cout<<"binwidth:"<<binwidth<<endl;
  //cout<<"integral:"<<integral<<endl;
  //cout<<"RMS:"<<RMS<<endl;
  cout<<"firstvalue:"<<firstvalue<<endl;
  gaustail->SetParameters(firstvalue,0,RMS1,RMS1/20,0.82);//Set the initial value of parameter
  //gaustail->FixParameter(0,firstvalue);
  return gaustail;

} 
//----------------------------------------------------------------------------
TF1 *CreateFuncGaussTail(TH1F *h,TString funcname,Double_t parone,Double_t partwo)
{
  TF1 *f = new TF1("f","gaus",-10000.,10000.);
  h->Fit("f","NR,Q");
  Double_t value =f->GetParameter(2);
  Double_t par1=1.1*parone;
  Double_t par2=partwo*1.2;//par2*(1.-0.2),par2*(1.+0.2));
  TF1 *gaustail=new TF1(funcname.Data(),"[0]*([4]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))+(1.-[4])/(2.*[3])*TMath::Exp(-TMath::Abs(x-[1])/[3]))",-2500.,2500.);
  //gaustail->SetParLimits(0,0.,10.);//Set the first parameter [0] range
  gaustail->SetParLimits(1,-50.,50.);//Set the second parameter [1] range
  //gaustail->SetParLimits(2,0.,par1);
  gaustail->SetParLimits(3,0,par2);
  gaustail->SetParLimits(4,0.,1.);
  Double_t binwidth=h->GetBinWidth(10);
  Double_t integral=h->Integral();
  Double_t RMS1=TMath::Abs(h->GetRMS());
  Double_t firstvalue=binwidth*integral;
  //cout<<"binwidth:"<<binwidth<<endl;
  //cout<<"integral:"<<integral<<endl;
  //cout<<"RMS:"<<RMS<<endl;
  //cout<<"firstvalue:"<<firstvalue<<endl;
  gaustail->SetParameters(1,0,value,partwo,0.82);//Set the initial value of parameter
 
  delete f;
  return gaustail;

} 
