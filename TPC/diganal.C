TH1F * h2max = 0;
void tpcstart(const  char * param= "Param2",const char * out="out.root",
	      Int_t nevent=0, Float_t th=3)
{


  gtpc.SetIO("galice.root", "out.root");

  AliTPCD * dig=new AliTPCD;
  dig.SetName(param);
  gtpc.GetIn()->cd();
  dig.Read(param);
  gtpc.SetDParam(dig);
  gtpc.SetTree(nevent);

  gtpc.SetbDelHisto(kTRUE);
  gtpc.SetThreshold(3);
}
  


void tpccfind(Int_t threshold, Float_t thr = 2,
	      Int_t i1 = 2, Int_t i2 = 2, 
	      Int_t tharea =20, Int_t thmax=100)
{  
  ///////////////////////GRAPHICS//DECLARATION//////////////////////
  TCanvas  * c1 = new TCanvas("padcluster","Cluster finder 1",700,900);
  c1->cd();
  TCanvas  * c2 = new TCanvas("padcluster2","Cluster finder 2",700,900);
  c2->cd();
  c1->cd();
  TPad * pad11 = new TPad("pad11","",0.01,0.76,0.48,0.95,21);
  pad11->Draw();
  TPad * pad12 = new TPad("pad12","",0.51,0.76,0.95,0.95,21);
  pad12->Draw();
  TPad * pad21 = new TPad("pad21","",0.01,0.56,0.49,0.74,21);
  pad21->Draw();
  TPad * pad22 = new TPad("pad22","",0.51,0.56,0.95,0.74,21);
  pad22->Draw();
  TPad * pad31 = new TPad("pad31","",0.01,0.36,0.49,0.54,21);
  pad31->Draw();
  TPad * pad32 = new TPad("pad32","",0.51,0.36,0.95,0.54,21);
  pad32->Draw();
  TPad * pad41 = new TPad("pad41","",0.01,0.16,0.49,0.34,21);
  pad41->Draw();
  TPad * pad42 = new TPad("pad42","",0.51,0.16,0.95,0.34,21);
  pad42->Draw();

  c2->cd();
  TPad * pad11_2 = new TPad("pad11_2","",0.01,0.76,0.48,0.95,21);
  pad11_2->Draw();
  TPad * pad12_2 = new TPad("pad12_2","",0.51,0.76,0.95,0.95,21);
  pad12_2->Draw();
  /////////////////////HISTOGRAMS///DECLARATION///////////////////////
  pad11->cd();
  TH1F * hsx = new TH1F("hsx","Sigma of distribution in time",40,0,2);
  pad12->cd();
  TH1F * hsy = new TH1F("hsy","Sigma of distribution in pads",40,0,2);

  pad21->cd();
  TProfile * hsx2 = new TProfile("hsx2","Sigma of distribution in time",
			 20,100,500);
  pad22->cd();
  TProfile * hsy2 = new TProfile("hsy2","Sigma of distribution in pads",
			 20,100,500);

  pad31->cd();
  TH1F * harea = new TH1F("harea","Area of the peak",26,0,52);
  pad32->cd();
  TH1F * hmax = new TH1F("hmax","Maximal amplitude in peak",30,0,150);
 

  pad41->cd();  
  TProfile * harea2= new TProfile("harea2","Area dependence z coordinata",
			 20,100,500);
  pad42->cd();
  TProfile * hmax2 = new TProfile("hmax2","Maximal amplitude dependence",
			 20,100,500);
  pad41->cd();  

  pad11_2->cd();
  TProfile * harea2p= new TProfile("harea2p","Area dependence on pad coordinata",
			 20,0,100);
  pad12_2->cd();
  TProfile * hmax2p = new TProfile("hmax2p","Maximal amplitude dependence on pad",
			 20,0,50);
  //////////////////CALCULATION//////////////////////////////////////////

  for (Int_t k = i1;k <=i2; k++)
    {
      tpcanal(1,k,10,0,kFALSE);
      TClusterFinder * cf=new TClusterFinder(0,0,threshold,1);
      cf->GetHisto(&gtpc.GetHis1());
      TClonesArray * arr = cf->FindClusters();
      cf->Delete();
      Int_t size = arr->GetEntries();
      
      if ( size>0 )   
	for (Int_t i=0 ; i<size;i++)
	  {
	    Int_t index;
	    TCluster *c=(TCluster*)arr->UncheckedAt(i);
	    hsx->Fill(TMath::Sqrt(c.fSigmaX2));
	    hsy->Fill(TMath::Sqrt(c.fSigmaY2));    
	    if  (TMath::Sqrt(c.fSigmaX2)<thr)
	      hsx2->Fill(c.fX,TMath::Sqrt(c.fSigmaX2),1);
	    if  (TMath::Sqrt(c.fSigmaY2)<thr)
	      hsy2->Fill(c.fX,TMath::Sqrt(c.fSigmaY2),1);       
	    hmax->Fill(c.fMax);
	    harea->Fill(c.fArea);
	    if (c.fArea<tharea)  harea2->Fill(c.fX,c.fArea,1);
            if (c.fMax<thmax) hmax2->Fill(c.fX,c.fMax,1);
	    if (c.fArea<tharea)  harea2p->Fill(c.fY,c.fArea,1);
            if (c.fMax<thmax) hmax2p->Fill(c.fY,c.fMax,1);
	    
	  }
    } 
  gStyle->SetOptStat(1);
  pad11->cd();
  hsx->Draw();
  pad12->cd();
  hsy->Draw();

  pad21->cd();
  hsx2->Draw();
  pad22->cd();
  hsy2->Draw();

  pad31->cd();
  harea->Draw();
  pad32->cd();
  hmax->Draw();

  pad41->cd();
  harea2->Draw();
  pad42->cd();
  hmax2->Draw();


  c1->cd();
  TPaveText * comment = new TPaveText(0.05,0.02,0.95,0.14,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);   
  comment->ReadFile("clusters.txt");
  comment->Draw();
  c2->cd();
  pad11_2->cd();
  harea2p->Draw();
  pad12_2->cd();
  hmax2p->Draw();


}


void tpcanal(Int_t sec, Int_t row, Int_t pad, Float_t * res=0, 
	     Bool_t bdraw = kTRUE,Float_t xmin=400,Float_t xmax=500)
{
  //calculate occupancy
  Double_t par[3];
  gtpc.SetSecRowTime(sec,row);
  gtpc.SetHisto(pad);
  //fit occupancy dependence

  gStyle->SetOptStat(0); 
  g1 = new TF1("pol0_r","pol0",xmin,xmax);  
  gtpc.GetHis3()->Fit("pol0_r","R0Q");  
  char s[100];
  
  if (bdraw == kTRUE) 
    {       
      gtpc.Draw("box");
      gtpc.GetPad3().cd();
      gtpc.GetPad3().SetGridx();
      gtpc.GetPad3().SetGridy();       
      gtpc.GetPad3().SetLogy();          
      gtpc.GetPad3().Draw();
      gtpc.GetPad2().cd();
      gtpc.GetPad2().SetGridx();
      gtpc.GetPad2().SetGridy();   
      fitText = new TPaveText(0.1,0.7,0.4,0.9,"NDC");
      gtpc.GetHis3()->Draw();
      sprintf(s,"p0 fit on interval %3.0f+- %3.0f",xmin,xmax);
      fitText->AddText(s);
    }
  g1->GetParameters(&par[0]);
  Float_t error = g1->GetParError(0);

  sprintf(s,"%0.3f+- %0.3f",par[0],error);
  if (bdraw == kTRUE)
    {
      gtpc.GetHis3()->Fit("pol0_r","R0Q");  
      fitText->AddText(s);
      fitText->Draw();
      gtpc.GetPad2().Update();  

      //plot histograms with specified options
      //move pads to another position be possible add text
      gtpc.GetPad1().SetPad(0.05,0.72,0.95,0.95);
      gtpc.GetPad2().SetPad(0.05,0.47,0.95,0.70); 
      gtpc.GetPad3().SetPad(0.05,0.22,0.95,0.45);  
      //add comments to the histograms 
      gtpc.GetCanvas().cd();
      TPaveText * comment = new TPaveText(0.05,0.03,0.95,0.2,"NDC");
      comment->SetTextAlign(12);
      comment->SetFillColor(42);   
      comment->ReadFile("comment.txt");
      comment->Draw();
    }  
  if (res != 0)
    {
      res[0] = par[0];
      res[1] = error;
      cout<<s<<"   "<<res[0]<<"   "<<res[1]<<"\n";
    }
}


void tpcanalall(Int_t isec =1, Int_t lstep = 1,Float_t tmin=400)
{
   AliTPCParam * tpcparam = gtpc->GetParam();
 //make window for displaying results
  TCanvas  * c_occu = new TCanvas("coccu","Occupancy dependence",700,900);
  c_occu->Update();
  TPad * pad1 = new TPad("occupancy","occupancy",0.05,0.25,0.95,0.95,21);
  pad1->Draw();
  //add comments to the histograms 
  TPaveText * comment = new TPaveText(0.05,0.03,0.95,0.2,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  comment->ReadFile("comment.txt");
  comment->Draw();
  //prepare histogram 
  Int_t irow = tpcparam->GetNRow(isec);
  Float_t xmin = tpcparam->GetPadRowRadii(isec,1);
  Float_t xmax = tpcparam->GetPadRowRadii(isec,irow);  
  pad1->cd();
  char s[220];
  char sh[220];
  sprintf(s,"occu_sector%d",isec);
  sprintf(sh,"Occupancy in sector %d as function of pad raw",isec);  
  TH1F * occudep = new TH1F(s,sh,300,xmin,xmax);
  Float_t   res[20];
  Float_t   x;
  for (Int_t i=2;i<irow;i+=lstep)
    { 
      tpcanal(isec,i,10,&res[0],kFALSE,tmin);
      x = tpcparam->GetPadRowRadii(isec,i) ;
      Int_t index = (300*(x-xmin))/(xmax-xmin);
      cout<<i<<"  "<<index<<"   "<<x<<"   "<<res[0]<<"   "<<res[1]<<"\n";  
      occudep->SetBinContent(index,res[0]);
      occudep->SetBinError(index,res[1]);      
    }
  //plot occupancy histogram
  
  pad1->SetGridx();
  pad1->SetGridy();
  gStyle->SetOptFit(0);       
  occudep->Draw("error");
  occudep->SetXTitle("pad row center position [cm]");
  occudep->SetYTitle("occupancy"); 
 
  //fit occupancy dependence
  //linear fit
  TF1 * g1 = new TF1("pol1_r","pol1");  
  occudep->Fit("pol1_r","+Q"); 
  Double_t par[3];
  Float_t error[3]; 
  Float_t chi;
  g1->GetParameters(&par[0]);
  error[0]=g1->GetParError(0);
  error[1]=g1->GetParError(1);   
  Float_t  chi = g1->GetChisquare();
  sprintf(s,"Linear fit     ocupancy = (%2.3f - %2.3f) +(%2.1f+- %2.1f).r   chi2 = %2.2f",
	  par[0],error[0],1000*par[1],1000*error[1],chi);
  comment->AddText(s);
 //(1-exp([0]1/(r*2+[1]**2)  fit
   TF1 * g1 = new TF1("polm1",occur,1,00,1);  
    occudep->Fit("polm1","+Q"); 
    Double_t par[3];
    Float_t error[3]; 
    g1->GetParameters(&par[0]);
    error[0]=g1->GetParError(0);
    //    error[1]=g1->GetParError(1);
    chi = g1->GetChisquare();
    sprintf(s,"(1-exp(P1/(x^2) fit   P1=(%2.3f+- %2.3f)    chi2=%2.2f ",
  	  par[0],error[0],chi);
    comment->AddText(s);
  c_occu->Update();
    
}

void tpcanalspectra(Int_t isec =1, Int_t r1= 2)
{
   AliTPCParam * tpcparam = gtpc->GetParam();
 //make window for displaying results
  TCanvas  * c_occu = new TCanvas("occuhis","Occupancy dependence",700,900);
  c_occu->Update();
  TPad * pad1 = new TPad("ocpad1","occupancy1",0.05,0.61,0.95,0.95,21);
  pad1->Draw();
  TPad * pad2 = new TPad("ocpad2","occupancy",0.05,0.61,0.95,0.95,21);
  pad2->Draw();

  TPad * pad3 = new TPad("ocpad3","occupancy1",0.05,0.25,0.95,0.60,21);
  pad3->Draw();
  TPad * pad4 = new TPad("ocpad4","occupancy",0.05,0.25,0.95,0.60,21);
  pad4->Draw();

  //add comments to the histograms 
  TPaveText * comment = new TPaveText(0.05,0.03,0.95,0.2,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  comment->ReadFile("comment.txt");
  comment->Draw();
  TH2F  his
  //prepare histogram 
  for (Int_t i=1;i<irow;i+=lstep)
    { 
      tpcanal(isec,i,10,&res[0],kFALSE,tmin);
    }
  //plot occupancy histogram
  
  pad1->SetGridx();
  pad1->SetGridy();
}




void tpcdraw(Int_t sec, Int_t row, Int_t pad)
{
   gStyle->SetOptStat(0); 
  //calculate occupancy for selected sector and pad row 
  //for selected pad is obtained signal shape 
  Double_t par[3];
  gtpc.SetSecRowTime(sec,row);
  gtpc.SetHisto(pad);
  gtpc.Draw("box");  
  //plot histograms with specified options
  //move pads to another position be possible add text
  gtpc.GetPad1().SetPad(0.05,0.72,0.95,0.95);
  gtpc.GetPad2().SetPad(0.05,0.47,0.95,0.70); 
  gtpc.GetPad3().SetPad(0.05,0.22,0.95,0.45);  
  //fit histogram of occupancy on specified range <150,500> 
  gtpc.GetPad2().cd();
  g1 = new TF1("pol0_r","pol0",150,500); 
  gtpc.GetHis3()->Fit("pol0_r","R0Q");    
  g1->GetParameters(&par[0]);
  Float_t error = g1->GetParError(0);
  fitText = new TPaveText(0.15,0.7,0.3,0.9,"NDC");
  fitText->AddText("p0 fit on interval <150-500>");
  char s[100];
  sprintf(s,"%0.3f+- %0.3f",par[0],error);
  fitText->AddText(s);
  fitText->Draw();
  gtpc.GetPad2().Update();     
  //set logarithmic 
  gtpc.GetPad3().cd();
  gtpc.GetPad3().SetLogy();  
   gtpc.GetPad3().Draw();    
  //add comments to the histograms 
  gtpc.GetCanvas().cd();
  TPaveText * comment = new TPaveText(0.05,0.03,0.95,0.2,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  comment->ReadFile("comment.txt");
  comment->Draw();
  gtpc.GetCanvas().Update();
  

}


void oDependence()
{
  //set plot options
  gStyle->SetOptFit(1); 
  gStyle->SetOptStat(1);  
  TCanvas  * c1 = new TCanvas("canPRF","Pad response function",700,900);
  TPad * pad1 = new TPad("pad1THR","",0.05,0.55,0.45,0.95,21);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2PRF","",0.55,0.55,0.95,0.95,21);
  pad2->Draw(); 
  TPad * pad3 = new TPad("pad3PRF","",0.55,0.05,0.95,0.45,21);
  pad3->Draw(); 
 
  pad1->cd();
  pad1->SetGridx();
  pad1->SetGridy();
  pad2->SetGridx();
  pad2->SetGridy();
  pad3->SetGridx();
  pad3->SetGridy();

  //make histogram of threshold dependence
  TH1F * hotd =new TH1F("Occupancy dependence on threshold",
			"Ocupancy at first pad row as function of threshold",
                        25,0.,25.);

  hotd->SetBinContent(5,0.625);
  hotd->SetBinError(5,0.02);
  hotd->SetBinContent(10,0.559);
  hotd->SetBinError(10,0.02); 
  hotd->SetBinContent(20,0.478);
  hotd->SetBinError(20,0.02);
  hotd->SetXTitle("Threshold   [channels]");
  hotd->SetYTitle("occupancy");
  hotd->Fit("pol1","+");  
  hotd->Draw("error");
  //make histogram of PRF  dependence
  TH1F * hoprfd =new TH1F("Occupancy dependence on PRF width",
			"Occupancy at first pad row as function of generic PRF sigma for  2.05x0.35 cm pad size ",
                        65, 0.,6.5);
  hoprfd->SetBinContent(10,0.492);
  hoprfd->SetBinError(10,0.02);

  hoprfd->SetBinContent(20,0.524);
  hoprfd->SetBinError(20,0.02); 

  hoprfd->SetBinContent(30,0.559);
  hoprfd->SetBinError(30,0.02);
  hoprfd->SetXTitle("Sigma of PRF   [mm]");
  hoprfd->SetYTitle("occupancy");
  pad2->cd();
  hoprfd->Fit("pol1","+");  
  hoprfd->Draw("error");
  pad2->Draw();
  //pad 3 histogram  
  pad3->cd();
   TH1F * hoprfd88 =new TH1F("Occupancy dependence on PRF width 08x08",
			"Occupancy at first pad row as function of generic PRF sigma for  0.8x0.8 cm pad size ",
                        65, 0.,6.5);

  hoprfd88->SetBinContent(20,0.322);
  hoprfd88->SetBinError(20,0.02);

  hoprfd88->SetBinContent(30,0.344);
  hoprfd88->SetBinError(30,0.02); 

  hoprfd88->SetBinContent(40,0.369);
  hoprfd88->SetBinError(40,0.02);
 
  hoprfd88->SetBinContent(60,0.416);
  hoprfd88->SetBinError(60,0.02);
  hoprfd88->SetXTitle("Sigma of PRF   [mm]");
  hoprfd88->SetYTitle("occupancy");
  hoprfd88->Fit("pol1","+");  
  hoprfd88->Draw("error");
  c1->cd();
  TPaveText * comment = new TPaveText(0.05,0.15,0.45,0.35,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);  
  comment->ReadFile("commentdep.txt");  
  comment->Draw();

}

void produceHisto()
{ 
  TH1F * hmostimp =new TH1F("number most important",
			"Mean value number of over threshold digits produced by \
most important  particle",
                        20,0.,23.);
  gStyle->SetOptStat(0); 
  hmostimp->Fill(1.,33.51);
  hmostimp->Fill(5.,34.32);
  hmostimp->Fill(10.,35.78);
  hmostimp->Fill(15.,36.85);
  hmostimp->Fill(20.,37.61);
  hmostimp->Write();
  hmostimp->SetMinimum(32.);
  delete hmostimp;
  
  TH1F * hall =new TH1F("number all3",
			"Mean value Number of over threshold digits produced by \
alll three stored particles",
                        20,0.,23.);
  gStyle->SetOptStat(0); 
  hall->Fill(1.,77.05);
  hall->Fill(5.,76.905);
  hall->Fill(10.,71.9);
  hall->Fill(15.,69.7);
  hall->Fill(20.,65.15);
  hall->Write();
  hall->SetMinimum(32.);    
}

void create6()
{
  TCanvas *cnumber =new TCanvas("Number","Number",600,900);    
  TPad *pad1 =new TPad("pad1","pad1",0.05,0.7,0.45,0.95);  
  pad1->cd();    
  pad1->Draw();
  gtpc.fout->Get("His_1_1").Draw(); 
  pad1->Update();
  TPad *pad2 =new TPad("pad2","pad2",0.55,0.7,0.95,0.95); 
  TPad *pad3 =new TPad("pad3","pad3",0.05,0.35,0.45,0.66);      
  TPad *pad4 =new TPad("pad4","pad4",0.55,0.35,0.95,0.66);
  TPad *pad5 =new TPad("pad5","pad5",0.05,0.05,0.45,0.32);      
  TPad *pad6 =new TPad("pad5","pad6",0.55,0.05,0.95,0.32);   
}

Double_t xtom1(Double_t *x, Double_t *par)
{ 
  Double_t xm=x[0]/100.;
  return (par[0]/(xm*xm));
}

Double_t xtomo(Double_t *x, Double_t *par)
{ 
  Double_t xm=x[0]/100.;
  return (par[0]/(xm**par[1]));
}

Double_t occur(Double_t *x, Double_t *par)
{
  Double_t xm=x[0]/100.;
  return  (1-exp(-par[0]/(xm**2)));
} 

void probability(Float_t param = 1, Float_t x1 = 0.9, Float_t x2 = 1.3, 
		 Float_t over = 2,const char * com ="", Int_t N=20)
{  

  //create canvas for drawing
  TCanvas  * c1 = new TCanvas("canprob","Pad response function",700,900);
  TPad * pad1 = new TPad("Theoretical probability","",0.05,0.22,0.95,0.95,21);
  pad1->Draw();
   
  //create histogram with estimated occupancy 
  //normalised to obtained occupancy at first pad
  Float_t y1=0;
  Float_t y2;
  char s[120];
  sprintf(s,"1-exp(-[1]/x**%f)",over);
  cout<<s<<"\n";
  TF1 *funr1 = new TF1("funr1",s,x1,x2);
  funr1->SetParameters(1,param);
  sprintf(s,"Probability  according 1-exp(-%1.3f/x**%1.1f distribution)",
	  param,over);
  pad1->cd();
  TH1F * hOccuT = new TH1F("hOccuT",s,5*N,x1,x2);
  Float_t x=x1;
  Float_t y;

  for (Float_t i = 0;i<N+1;i++)
    {     
      y = funr1->Eval(x);
      hOccuT->Fill(x,y);
      x+=(x2-x1)/Float_t(N);
    };
  //fitting calculated dependence with linear fit and with
  //generic function
  pad1->cd();
  sprintf(s,"[1]/(x**%1.1f)",over);
  TF1 *lin1 = new TF1("lin1","pol1",x1,x2);
  lin1->SetLineColor(2);
  lin1->SetLineWidth(5);
  lin1->SetLineStyle(1);
  hOccuT->Draw();
  hOccuT->Fit("lin1","S+");

  sprintf(s,"[1]/(x**%1.1f)",over);
  TF1 *funorig = new TF1("funorig",s,x1,x2);
  funorig->SetLineColor(3);
  funorig->SetLineWidth(5);
  funorig->SetLineStyle(2);
  hOccuT->Fit("funorig","S+");
  
  //find minimum and maximum and scale histo  
  if (y1 == 0)  
    {
      Float_t ymin,ymax;
      y1=lin1->Eval(x2);
      y2=lin1->Eval(x1);
      ymin= funorig->Eval(x2);
      ymax= funorig->Eval(x1);
      if (ymin<y1) y1=ymin;
      if (ymax>y2) y2=ymax;
    }
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0); 
  hOccuT->SetMaximum(y2);
  hOccuT->SetMinimum(y1); 
  hOccuT->SetXTitle("r position [m]");
  hOccuT->SetYTitle("probability");  
  //add comments to the histograms 
  c1->cd();
  TPaveText * comment = new TPaveText(0.05,0.03,0.95,0.20,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  TText *title = comment->AddText("Estimation of occupancy dependence on R position");
  title->SetTextSize(0.04);
  comment->AddText("Observed efect of probability saturation");
  sprintf(s,"Supposed generic flux dependence : %1.3f/(x**%1.1f)",param,over);
  comment->AddText(s);
  comment->AddText("Probility : 1-exp(-flux*mean particle \"pad x time\" area)");
  comment->AddText("Full line  linear fit ");
  comment->AddText("Dashed line : fit by generic flux function ");
  
  comment->AddText(com);
  comment->Draw();
}



void digamp(Float_t t1, Float_t t2, Float_t p1, Float_t p2, Float_t th)
{
  gStyle->SetOptStat(1); 
  TH1F * h2max = new TH1F("all amplitudes","all amplitude", 20, 1,600);
  for (Int_t itime = t1;itime<t2;itime++)
  for (Int_t ipad = p1;ipad<p2;ipad++)  
    {
      Int_t index = gtpc.GetHis1().GetBin(itime,ipad);
      Float_t weight = gtpc.GetHis1().GetBinContent(index);
      //      cout<<itime<<"\t"<<ipad<<"\t"<<weight<<"\n";
      if (weight > th) h2max->Fill(weight);
    };  
   h2max->Draw();
}            

void digampM(Float_t t1, Float_t t2, Float_t p1, Float_t p2, Float_t th)
{
  gStyle->SetOptStat(1); 

  //create canvas for drawing
  //  TCanvas  * c1 = new TCanvas("dh","Amplitude",700,900);
  //TPad * pad1 = new TPad("amplitude","",0.05,0.22,0.95,0.95,21);
  //pad1->Draw();

  if (h2max == 0) h2max = new TH1F("max amplitudes","max amplitude", 25, 1,250);
  for (Int_t itime = t1;itime<t2;itime++)
  for (Int_t ipad = p1;ipad<p2;ipad++)  
    {
      Bool_t bmax = kTRUE;
      Int_t index = gtpc.GetHis1().GetBin(itime,ipad);
      Float_t weight = gtpc.GetHis1().GetBinContent(index);
      if (weight>th)
	{
	  for (Int_t i = -1;i<2;i++) 
	    for (Int_t j = -1;j<2;j++) 	  
	      if (!((0==i)&&(0==j)))
		{
		  index = gtpc.GetHis1().GetBin(itime+i,ipad+j);
		  Float_t weightl = gtpc.GetHis1().GetBinContent(index);
		  if (!(weightl<weight)) bmax = kFALSE;
		}     
	  if (kTRUE==bmax) h2max->Fill(weight);
	}
    };  
   h2max->Draw("error");

}            



void digampMALL(Float_t t1, Float_t t2, Float_t p1, Float_t p2, Float_t th )
{
  for (Int_t i=1;i<20;i++)
    {
      tpcanal(1,i,10);
      digampM(t1,t2,p1,p2,th);
    }
}
