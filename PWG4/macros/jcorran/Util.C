/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//Functions used by the other macros
//Author: Jason Glyndwr Ulery, ulery@uni-frankfurt.de

void MakeProjections(Float_t &xTPt1, Float_t &xTPt2, Float_t &xAPt1, Float_t &xAPt2, Int_t xCent, TList *xList, TH1F *xOutHistPhi, TH1F *xOutHistPhiMix, TH1F *xOutHistEtaN, TH1F *xOutHistEtaNMix, TH1F *xOutHistEtaA, TH1F *xOutHistEtaAMix, TH2F *xOutHistPhiEta, TH2F *xOutHistPhiEtaMix, Float_t xMPt[4], Int_t xptweighted=0, Int_t xMC=0, Int_t xNotDelta=0, Int_t xSign=0){//(Trigger Pt Low, Trigger Pt High, Assoc Pt Low, Assoc Pt High, Centrality Bin, List, Hists, mean pt array, Pt weighted?, Monty Carlo?, gives phi from triggered and from all for delta phi and delta phi mixed respectively) 
  // TFile *xFile=new TFile(xFileName);
  //TList *xList=xFile->Get("julery_DiHadron");
  //This needs changed if Input binning is changed
  //cout << xptweighted << " " << xMC << " " << xNotDelta << " " << xSign << endl;
  const int xnpt=12;
  //float xptbins[(xnpt+1)]={2.5,3,4,6,8,10,15,20,30,40,50,75,100};
  //IF you change also change in EffFit
  float xptbins[(xnpt+1)]={2,2.5,3,4,5,6,8,10,15,20,30,40,50};
  const int xapt=30;
  float xaptbins[(xapt+1)]={0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100};
  float xPi=3.1415926535898;
  
  //After this no changes should be needed
  float xsum1,xsum2,xerr1,xerr2;
  int xtpt1=-1,xtpt2=-1;
  TH2F *xtemphist1;
  TH2F *xtemphist2;
  TH3F *xtemphist3;
  TH3F *xtemphist4;
  TH2F *xtemphistEN;
  TH2F *xtemphistENM;
  TH2F *xtemphistEA;
  TH2F *xtemphistEAM;
  int xapt1,xapt2,xapt3,xapt4;
  
  for(int xi=0;xi<=xnpt;xi++){
    if(fabs(xTPt1-xptbins[xi])<0.2)xtpt1=xi;
    if(fabs(xTPt2-xptbins[xi])<0.2)xtpt2=xi;
  }
  //cout << xtpt1 << " " << xtpt2 << endl;
  if(xtpt1<0)cout << "Invalid Trigger Pt1" << endl;
  if(xtpt2<0)cout << "Invalid Trigger Pt2" << endl;
  if(xtpt1<0||xtpt2<0)break;
  
  
  if(xAPt2>xTPt1){
    cout << "Associateds must be less then trigger (Automatically forced)" << endl;
    xAPt2=xTPt1;
  }
  for(int xi=0;xi<=xapt;xi++){
    if(fabs(xAPt1-xaptbins[xi])<0.1)xapt1=xi+1;
    if(fabs(xAPt2-xaptbins[xi])<0.1)xapt2=xi;
  }
  // cout << xapt1 << " " << xapt2 << endl;
  //since there are some tolerences reset the values to what are actually used
  // xAPt1=(xapt1-1)*xasswidth;
  // xAPt2=(xapt2-1)*xasswidth;
  xAPt1=xaptbins[xapt1-1];
  xAPt2=xaptbins[xapt2];
  xTPt1=xptbins[xtpt1];
  xTPt2=xptbins[xtpt2];
  cout << "xTPt: "<< xTPt1 <<  " - " << xTPt2 << endl; 
  cout << "xAPt: " << xAPt1 << " - " << xAPt2 << endl;
  float xntrig=0,xnmix=0,xntrigpt=0;
  float exntrig=0,exntrigpt=0;
  char xname[150];
  char *cpt1[2]={"","Pt"};
  char *cdelta2[2]={"DeltaPhiEta","PhiEtaTrig"};
  char *cdelta1[2]={"DeltaPhi","PhiTrig"};
  char *cmc1[2]={"","_MC"};
  char *sign1[3]={"","_LS","_ULS"};
  char *sign2[3]={""," LikeSign"," UnlikeSign"};
  char *xEtaN[2]={"DeltaEtaN","EtaTrig"};
  char *xEtaA[2]={"DeltaEtaA","EtaTrig"};
  
  
  sprintf(xname,"fHistNTrigger_C%d%s",xCent,cmc1[xMC]);
  TH1F *xNTrig=(TH1F*)xList->FindObject(xname);
  sprintf(xname,"fHistNMix_C%d%s",xCent,cmc1[xMC]);
  TH1F *xNMix=(TH1F*)xList->FindObject(xname);
  sprintf(xname,"fHistNTriggerPt_C%d%s",xCent,cmc1[xMC]);
  TH1F *xNTrigPt=(TH1F*)xList->FindObject(xname);
  sprintf(xname,"fHistMult%s",cmc1[xMC]);
  TH1F *xMult=(TH1F*)xList->FindObject(xname);
  int xnevent=xMult->GetEntries();
  
  for(int xi=xtpt1;xi<xtpt2;xi++){//Do for all trigger pt ranges we are adding
    sprintf(xname,"fHist%s%s_P%d_C%d%s%s",cdelta1[xNotDelta],cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    xtemphist1=(TH2F*)xList->FindObject(xname);
    sprintf(xname,"fHistDeltaPhiMix%s_P%d_C%d%s%s",cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    if(xNotDelta) sprintf(xname,"fHistPhi%s_C%d%s",cpt1[xptweighted],xCent,cmc1[xMC]);
    xtemphist2=(TH2F*)xList->FindObject(xname);
    
    sprintf(xname,"fHist%s%s_P%d_C%d%s%s",xEtaN[xNotDelta],cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    xtemphistEN=(TH2F*)xList->FindObject(xname);
    
    sprintf(xname,"fHistDeltaEtaNMix%s_P%d_C%d%s%s",cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    if(xNotDelta) sprintf(xname,"fHistEta%s_C%d%s",cpt1[xptweighted],xCent,cmc1[xMC]);
    xtemphistENM=(TH2F*)xList->FindObject(xname);
    
    
    sprintf(xname,"fHist%s%s_P%d_C%d%s%s",xEtaA[xNotDelta],cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    xtemphistEA=(TH2F*)xList->FindObject(xname);
    sprintf(xname,"fHistDeltaEtaAMix%s_P%d_C%d%s%s",cpt1[xptweighted],xi,xCent,sign1[xSign],cmc1[xMC]);
    if(xNotDelta) sprintf(xname,"fHistEta%s_C%d%s",cpt1[xptweighted],xCent,cmc1[xMC]);
    xtemphistEAM=(TH2F*)xList->FindObject(xname);
    
    sprintf(xname,"fHist%s%s_P%d_C%d%s",cdelta2[xNotDelta],cpt1[xptweighted],xi,xCent,cmc1[xMC]);
    xtemphist3=(TH3F*)xList->FindObject(xname);
    sprintf(xname,"fHistDeltaPhiEtaMix%s_P%d_C%d%s",cpt1[xptweighted],xi,xCent,cmc1[xMC]);
    if(xNotDelta) sprintf(xname,"fHistPhiEta%s_C%d%s",cpt1[xptweighted],xCent,cmc1[xMC]);
    xtemphist4=(TH3F*)xList->FindObject(xname);
    //c12355=new TCanvas("c12355","c12355");
    // xtemphist4->Draw();
    
    xntrig+=xNTrig->GetBinContent(xi+1);
    xnmix+=xNMix->GetBinContent(xi+1);
    exntrig+=pow(xNTrig->GetBinError(xi+1),2);
    xntrigpt+=xNTrigPt->GetBinContent(xi+1);
    exntrigpt+=pow(xNTrigPt->GetBinError(xi+1),2);
    if(xi==xtpt1){
      TH2F *xDeltaPhi=(TH2F*)xtemphist1->Clone();
      xDeltaPhi->SetName("xDeltaPhi");
      TH2F *xDeltaPhiMix=(TH2F*)xtemphist2->Clone();
      xDeltaPhiMix->SetName("xDeltaPhiMix");
      int xnbins=xDeltaPhi->GetNbinsX();
      float xmin=xDeltaPhi->GetBinCenter(1)-(xDeltaPhi->GetBinWidth(1))/2.;
      float xmax=xDeltaPhi->GetBinCenter(xnbins)+(xDeltaPhi->GetBinWidth(xnbins))/2.;
      
      TH2F *xDeltaEtaN=(TH2F*)xtemphistEN->Clone();
      xDeltaEtaN->SetName("xDeltaEtaN");
      
      TH2F *xDeltaEtaNMix=(TH2F*)xtemphistENM->Clone();
      xDeltaEtaNMix->SetName("xDeltaEtaNMix");
      int xnEbins=xDeltaEtaN->GetNbinsX();
      //cout << xnEbins << "EtaBins" << endl;
      float xEmin=xDeltaEtaN->GetBinCenter(1)-(xDeltaEtaN->GetBinWidth(1))/2.;
      float xEmax=xDeltaEtaN->GetBinCenter(xnEbins)+(xDeltaEtaN->GetBinWidth(xnEbins))/2.; 
      
      TH2F *xDeltaEtaA=(TH2F*)xtemphistEA->Clone();
      xDeltaEtaA->SetName("xDeltaEtaA");
      TH2F *xDeltaEtaAMix=(TH2F*)xtemphistEAM->Clone();
      xDeltaEtaAMix->SetName("xDeltaEtaAMix");
      
      TH3F *xDeltaPhiEta=(TH3F*)xtemphist3->Clone();
      xDeltaPhiEta->SetName("xDeltaPhiEta");
      TH3F *xDeltaPhiEtaMix=(TH3F*)xtemphist4->Clone();
      xDeltaPhiEtaMix->SetName("xDeltaPhiEtaMix");
      TAxis *xAxis=xDeltaPhiEta->GetXaxis();
      TAxis *xAxis2=xDeltaPhiEta->GetXaxis();
      TAxis *yAxis2=xDeltaPhiEta->GetYaxis();
      int xnbins2=xAxis2->GetNbins();
      int ynbins2=yAxis2->GetNbins();
      float xmin2=xAxis2->GetBinCenter(1)-(xAxis2->GetBinWidth(1))/2.;
      float xmax2=xAxis2->GetBinCenter(xnbins2)+(xAxis2->GetBinWidth(xnbins2))/2.;
      float ymin2=yAxis2->GetBinCenter(1)-(yAxis2->GetBinWidth(1))/2.;
      float ymax2=yAxis2->GetBinCenter(ynbins2)+(yAxis2->GetBinWidth(ynbins2))/2.;
      //cout << "Sum Phi " << xDeltaPhi->GetSum() << endl;
    }
    else{
      for(int xj=1;xj<=xnbins;xj++){
	for(int xk=xapt1;xk<=xapt2;xk++){
	  xDeltaPhi->SetBinContent(xj,xk,(xDeltaPhi->GetBinContent(xj,xk)+xtemphist1->GetBinContent(xj,xk)));
	  xDeltaPhi->SetBinError(xj,xk,sqrt(pow(xDeltaPhi->GetBinError(xj,xk),2)+pow(xtemphist1->GetBinError(xj,xk),2)));
	  if(!xNotDelta)xDeltaPhiMix->SetBinContent(xj,xk,(xDeltaPhiMix->GetBinContent(xj,xk)+xtemphist2->GetBinContent(xj,xk)));
	  if(!xNotDelta)xDeltaPhiMix->SetBinError(xj,xk,sqrt(pow(xDeltaPhiMix->GetBinError(xj,xk),2)+pow(xtemphist2->GetBinError(xj,xk),2)));			
	}
      }//end phi loop
      
      for(int xj=1;xj<=xnEbins;xj++){
	for(int xk=xapt1;xk<=xapt2;xk++){
	  xDeltaEtaN->SetBinContent(xj,xk,(xDeltaEtaN->GetBinContent(xj,xk)+xtemphistEN->GetBinContent(xj,xk)));
	  xDeltaEtaN->SetBinError(xj,xk,sqrt(pow(xDeltaEtaN->GetBinError(xj,xk),2)+pow(xtemphistEN->GetBinError(xj,xk),2)));
	  if(!xNotDelta)xDeltaEtaNMix->SetBinContent(xj,xk,(xDeltaEtaNMix->GetBinContent(xj,xk)+xtemphistENM->GetBinContent(xj,xk)));
	  if(!xNotDelta)xDeltaEtaNMix->SetBinError(xj,xk,sqrt(pow(xDeltaEtaNMix->GetBinError(xj,xk),2)+pow(xtemphistENM->GetBinError(xj,xk),2)));
	  
	  xDeltaEtaA->SetBinContent(xj,xk,(xDeltaEtaA->GetBinContent(xj,xk)+xtemphistEA->GetBinContent(xj,xk)));
	  xDeltaEtaA->SetBinError(xj,xk,sqrt(pow(xDeltaEtaA->GetBinError(xj,xk),2)+pow(xtemphistEA->GetBinError(xj,xk),2)));
	  if(!xNotDelta)xDeltaEtaAMix->SetBinContent(xj,xk,(xDeltaEtaAMix->GetBinContent(xj,xk)+xtemphistEAM->GetBinContent(xj,xk)));
	  if(!xNotDelta)xDeltaEtaAMix->SetBinError(xj,xk,sqrt(pow(xDeltaEtaAMix->GetBinError(xj,xk),2)+pow(xtemphistEAM->GetBinError(xj,xk),2)));
	}
      }//end eta loop
      
      for(int xj=1;xj<=xnbins2;xj++){
	for(int xk=1;xk<=ynbins2;xk++){
	  for(int xl=xapt1;xl<=xapt2;xl++){
	    xDeltaPhiEta->SetBinContent(xj,xk,xl,(xDeltaPhiEta->GetBinContent(xj,xk,xl)+xtemphist3->GetBinContent(xj,xk,xl)));
	    xDeltaPhiEta->SetBinError(xj,xk,xl,sqrt(pow(xDeltaPhiEta->GetBinError(xj,xk,xl),2)+pow(xtemphist3->GetBinError(xj,xk,xl),2)));
	    if(!xNotDelta)xDeltaPhiEtaMix->SetBinContent(xj,xk,xl,(xDeltaPhiEtaMix->GetBinContent(xj,xk,xl)+xtemphist4->GetBinContent(xj,xk,xl)));
	    if(!xNotDelta)xDeltaPhiEtaMix->SetBinError(xj,xk,xl,sqrt(pow(xDeltaPhiEtaMix->GetBinError(xj,xk,xl),2)+pow(xtemphist4->GetBinError(xj,xk,xl),2)));
	  }
	}
      }//end eta-phi loop
    }//else
  }//xi 
  
  /*
    c124=new TCanvas("c124");
    xtemphist4->Draw();
    cout << xtemphist4->GetBinContent(10,10,2) << " " << xDeltaPhiEtaMix->GetBinContent(10,10,2) << endl;
    cout << xtemphist4->GetBinError(10,10,2) << " " << xDeltaPhiEtaMix->GetBinError(10,10,2) << endl;
    cout << xtemphist4->GetBinContent(10,1,2) << " " << xDeltaPhiEtaMix->GetBinContent(10,1,2) << endl;
  */
  
  
  cout << "Number of Triggers: " << xntrig << endl;
  cout << "Number of Mixed Events: " << xnmix << endl;
  sprintf(xname,"fHistNEvents_C%d",xCent);
  TH1F *xNEvents=(TH1F*)xList->FindObject(xname);
  cout << "Number of Events: " << xNEvents->GetBinContent(1) << endl;
  if(xptweighted)cout << "<P_T Trigger>: " << xntrigpt/xntrig << endl;
  //cout << xntrigpt << " " << exntrigpt << " " << xntrig << " " << exntrig << endl;
  if(xptweighted)xMPt[0]=xntrigpt/xntrig;
  exntrigpt=sqrt(exntrigpt);
  exntrig=sqrt(exntrig);
  if(xptweighted)xMPt[1]=xMPt[0]*sqrt(pow(exntrigpt/xntrigpt,2)+pow(exntrig/xntrig,2));
  xMPt[2]=xntrig;
  xMPt[3]=xnmix;
  xDeltaPhi->Scale(xnbins/(2.*xPi)/xntrig);
  xDeltaEtaN->Scale(xnEbins/(xEmax-xEmin)/xntrig);
  xDeltaEtaA->Scale(xnEbins/(xEmax-xEmin)/xntrig);
  xDeltaPhiEta->Scale(xnbins2/(2.*xPi)*ynbins2/(ymax2-ymin2)/xntrig);
  if(!xNotDelta){
    xDeltaPhiMix->Scale(xnbins/(2.*xPi)/xnmix);
    xDeltaEtaNMix->Scale(xnEbins/(xEmax-xEmin)/xnmix);
    xDeltaEtaAMix->Scale(xnEbins/(xEmax-xEmin)/xnmix);
    xDeltaPhiEtaMix->Scale(xnbins2/(2.*xPi)*ynbins2/(ymax2-ymin2)/xnmix);
  }
  else{
    xDeltaPhiMix->Scale(xnbins/(2.*xPi)/xnevent);
    xDeltaEtaNMix->Scale(xnEbins/(xEmax-xEmin)/xnevent);
    xDeltaEtaAMix->Scale(xnEbins/(xEmax-xEmin)/xnevent);
    xDeltaPhiEtaMix->Scale(xnbins2/(2.*xPi)*ynbins2/(ymax2-ymin2)/xnevent);
  }
  
  char *tit1="#Delta#phi";
  if(xNotDelta)tit1="#phi";
  char *tit11="";
  if(xptweighted)tit11="p_{T} Weighted";
  char *tit12="";
  if(xNotDelta)tit12="Triggered";
  char *tit2="1";
  if(xptweighted)tit2="p_{T}";
  char *tit3="Mixed #Delta#phi";
  if(xNotDelta)tit3="#phi";
  
  sprintf(xname,"%s %s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit12,tit11,tit1,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistPhi=new TH1F("dPhi",xname,xnbins,xmin,xmax);
  xOutHistPhi->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistPhi->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit1);
  xOutHistPhi->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistPhi);
  xOutHistPhi->SetMarkerStyle(20);
  xOutHistPhi->SetMarkerColor(2);
  xOutHistPhi->SetLineColor(2);
  
  sprintf(xname,"%s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit11,tit3, xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistPhiMix=new TH1F("dMix",xname,xnbins,xmin,xmax);
  xOutHistPhiMix->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistPhiMix->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit1);
  xOutHistPhiMix->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistPhiMix);
  xOutHistPhiMix->SetMarkerStyle(21);
  xOutHistPhiMix->SetMarkerColor(1);
  xOutHistPhiMix->SetLineColor(1);
  
  char *tit14="#Delta#eta";
  if(xNotDelta) tit14="#eta";
  char *tit15="Mixed #Delta#eta";
  if(xNotDelta) tit15="#eta";
  char *titNear="Near-Side";
  if(xNotDelta)titNear="";
  
  sprintf(xname,"%s %s %s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit12,titNear,tit11,tit14,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistEtaN=new TH1F("dEtaN",xname,xnEbins,xEmin,xEmax);
  xOutHistEtaN->Sumw2();
  sprintf(xname,"%s     ",tit14);
  xOutHistEtaN->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit14);
  xOutHistEtaN->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistEtaN);
  xOutHistEtaN->SetMarkerStyle(20);
  xOutHistEtaN->SetMarkerColor(2);
  xOutHistEtaN->SetLineColor(2);
  
  sprintf(xname,"%s %s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit11,titNear,tit15,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistEtaNMix=new TH1F("dEtaMixN",xname,xnEbins,xEmin,xEmax);
  xOutHistEtaNMix->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistEtaNMix->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit14);
  xOutHistEtaNMix->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistEtaNMix);
  xOutHistEtaNMix->SetMarkerStyle(21);
  xOutHistEtaNMix->SetMarkerColor(1);
  xOutHistEtaNMix->SetLineColor(1);
  
  char *titAway="Away-Side";
  if(xNotDelta)titAway="";
  
  sprintf(xname,"%s %s %s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit12,titAway,tit11,tit14,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistEtaA=new TH1F("dEtaA",xname,xnEbins,xEmin,xEmax);
  xOutHistEtaA->Sumw2();
  sprintf(xname,"%s     ",tit14);
  xOutHistEtaA->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit14);
  xOutHistEtaA->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistEtaA);
  xOutHistEtaA->SetMarkerStyle(20);
  xOutHistEtaA->SetMarkerColor(2);
  xOutHistEtaA->SetLineColor(2);
  
  sprintf(xname,"%s %s %s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit11,titAway,tit15,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistEtaAMix=new TH1F("dEtaMixA",xname,xnEbins,xEmin,xEmax);
  xOutHistEtaAMix->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistEtaAMix->GetXaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%s}   ",tit2,tit14);
  xOutHistEtaAMix->GetYaxis()->SetTitle(xname);
  SetTitles1D(xOutHistEtaAMix);
  xOutHistEtaAMix->SetMarkerStyle(21);
  xOutHistEtaAMix->SetMarkerColor(1);
  xOutHistEtaAMix->SetLineColor(1);
  
  sprintf(xname,"%s %s %s-%s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit12,tit11,tit1,tit14,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistPhiEta=new TH2F("dPhiEta",xname,xnbins2,xmin2,xmax2,ynbins2,ymin2,ymax2);
  xOutHistPhiEta->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistPhiEta->GetXaxis()->SetTitle(xname);
  sprintf(xname,"%s     ",tit14);
  xOutHistPhiEta->GetYaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%sd%s}   ",tit2,tit1,tit14);
  xOutHistPhiEta->GetZaxis()->SetTitle(xname);
  SetTitles2D(xOutHistPhiEta);
  
  sprintf(xname,"%s %s-%s Distribution for %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",tit11,tit3,tit14,xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
  *xOutHistPhiEtaMix=new TH2F("dPhiEtaMix",xname,xnbins2,xmin2,xmax2,ynbins2,ymin2,ymax2);
  xOutHistPhiEtaMix->Sumw2();
  sprintf(xname,"%s (radians)     ",tit1);
  xOutHistPhiEtaMix->GetXaxis()->SetTitle(xname);
  sprintf(xname,"%s     ",tit14);
  xOutHistPhiEtaMix->GetYaxis()->SetTitle(xname);
  sprintf(xname,"#frac{%s}{N_{Trig}}#frac{dN}{d%sd%s}   ",tit2,tit1,tit14);
  xOutHistPhiEtaMix->GetZaxis()->SetTitle(xname);
  SetTitles2D(xOutHistPhiEtaMix);
  
  float xsum1,xsum2,xerr1,xerr2,xsum3,xsum4,xerr3,xerr4;
  for(int xi=1;xi<=xnbins;xi++){
    xsum1=0; xerr1=0; xsum2=0; xerr2=0;
    for(int xj=xapt1; xj<=xapt2;xj++){
      xsum1+=xDeltaPhi->GetBinContent(xi,xj);
      xerr1+=pow(xDeltaPhi->GetBinError(xi,xj),2);
      xsum2+=xDeltaPhiMix->GetBinContent(xi,xj);
      xerr2+=pow(xDeltaPhiMix->GetBinError(xi,xj),2);
    }
    
    xOutHistPhi->SetBinContent(xi,xsum1);
    xOutHistPhi->SetBinError(xi,sqrt(xerr1));
    xOutHistPhiMix->SetBinContent(xi,xsum2);
    xOutHistPhiMix->SetBinError(xi,sqrt(xerr2));
  }//philoop
  for(int xi=1;xi<=xnEbins;xi++){
    xsum1=0; xerr1=0; xsum2=0; xerr2=0; 
    xsum3=0; xerr3=0; xsum4=0; xerr4=0;
    //cout << "apt" << xapt1 << " " << xapt2 << endl;
    for(int xj=xapt1; xj<=xapt2;xj++){
      xsum1+=xDeltaEtaN->GetBinContent(xi,xj);
      xerr1+=pow(xDeltaEtaN->GetBinError(xi,xj),2);
      xsum2+=xDeltaEtaNMix->GetBinContent(xi,xj);
      xerr2+=pow(xDeltaEtaNMix->GetBinError(xi,xj),2);
      xsum3+=xDeltaEtaA->GetBinContent(xi,xj);
      xerr3+=pow(xDeltaEtaA->GetBinError(xi,xj),2);
      xsum4+=xDeltaEtaAMix->GetBinContent(xi,xj);
      xerr4+=pow(xDeltaEtaAMix->GetBinError(xi,xj),2);
    }
    xOutHistEtaN->SetBinContent(xi,xsum1);
    xOutHistEtaN->SetBinError(xi,sqrt(xerr1));
    xOutHistEtaNMix->SetBinContent(xi,xsum2);
    xOutHistEtaNMix->SetBinError(xi,sqrt(xerr2));
    xOutHistEtaA->SetBinContent(xi,xsum3);
    xOutHistEtaA->SetBinError(xi,sqrt(xerr3));
    xOutHistEtaAMix->SetBinContent(xi,xsum4);
    xOutHistEtaAMix->SetBinError(xi,sqrt(xerr4));
  }//etaloop
  
  for(int xi=1;xi<=xnbins2;xi++){
    for(int xj=1;xj<=ynbins2;xj++){
      xsum1=0; xerr1=0; xsum2=0; xerr2=0;
      for(int xk=xapt1; xk<=xapt2;xk++){
        xsum1+=xDeltaPhiEta->GetBinContent(xi,xj,xk);
	xerr1+=pow(xDeltaPhiEta->GetBinError(xi,xj,xk),2);
	xsum2+=xDeltaPhiEtaMix->GetBinContent(xi,xj,xk);
	xerr2+=pow(xDeltaPhiEtaMix->GetBinError(xi,xj,xk),2);
      }
      //cout << xDeltaPhiEtaMix->GetYaxis()->GetBinCenter(xj) << " " << xOutHistPhiEtaMix->GetYaxis()->GetBinCenter(xj) << " " << xsum2 << endl;
      xOutHistPhiEta->SetBinContent(xi,xj,xsum1);
      xOutHistPhiEta->SetBinError(xi,xj,sqrt(xerr1));
      xOutHistPhiEtaMix->SetBinContent(xi,xj,xsum2);
      xOutHistPhiEtaMix->SetBinError(xi,xj,sqrt(xsum2));
    }
  }//phi-eta loop
  
  //c10=new TCanvas("c10","",800,600);
  //xOutHistPhi->Draw();
  //delete xDeltaPhi;
  // delete xDeltaMix;
  //delete xtemphist1;
  //delete xtemphist2;
  //xFile->Close();
  
}

//--------------------------------------------------------------------------------------------------------		 

void ZYA1(TH1F *yHistSig, TH1F *yHistBg, float scalea[3], float ZYAMCent=1.,float ZYAMWidth=0.2){
  //yerr=0 no error
  //yerr=1 no mixed event
  //yerr=2 no signal
  //yerr=3 no signal and no mixed
  
  //float ZYAMCent=1.3;
  //float ZYAMWidth=0.2;//window is cent +/- width so its really a 1/2 width of the window
  float a,yea,yerr;
  //static float scalea[3];
  float ysumsig=0,ysumbg=0,yerrsig=0,yerrbg=0;
  float ycent;
  int ybins=yHistSig->GetNbinsX();
  float yPi=3.1415926535898;
  //cout << ybins << endl;
  //cout << yHistSig->GetSum() << endl;
  //c20=new TCanvas("c20","",800,600);
  //yHistSig->Draw();
  cout << "ZYAMRegion: " << (ZYAMCent-ZYAMWidth) << "-" << (ZYAMCent+ZYAMWidth) << endl;
  
  for(int yi=1;yi<=ybins;yi++){
    ycent=yHistSig->GetBinCenter(yi);
    if((fabs(fabs(ycent)-ZYAMCent)<ZYAMWidth)||(fabs(fabs(ycent)-(2*yPi-ZYAMCent))<ZYAMWidth)){
      ysumsig+=yHistSig->GetBinContent(yi);
      yerrsig+=pow(yHistSig->GetBinError(yi),2);
      ysumbg+=yHistBg->GetBinContent(yi);
      yerrbg+=pow(yHistBg->GetBinError(yi),2);
    }
  }
  cout << ysumsig << " " << yerrsig << " " << ysumbg << " " << yerrbg << endl;
  if(ysumbg==0&&ysumsig==0){
    a=0;
    yea=0;
    yerr=3;
  }
  else if(ysumbg==0){
    a=0;
    yea=0;
    yerr=1;
    
  }
  else if(ysumsig==0){
    a=0;
    yea=0;
    yerr=2;
  }
  else{
    a=ysumsig/ysumbg;
    yea=a*sqrt(yerrsig/ysumsig/ysumsig+yerrbg/ysumbg/ysumbg);
    yerr=0;
  }
  cout << a << endl;
  scalea[0]=a;
  scalea[1]=yea;
  scalea[2]=yerr;
  if(yerr) cout << "ZYA1 Error:" << yerr << endl;
  //return scalea;
}
//----------------------------------------------------------------------------------------------------
void ZYAM2D(TH2F *yHistSig, TH2F *yHistBg, Float_t scalea[3], Int_t nMin=10, Int_t nIter=10){
  //for pp nIter=1 might be fine, for heavy-ions where the background has shape it shouldn't

  const Int_t NMin=100;
 
  float sumSig, sumBg;
  float esumSig, esumBg;
  float sumMin, esumMin;
  scalea[0]=1;
  scalea[1]=0;
  float MinArray[NMin][2];
  float MinArraySig[NMin][2];
  float MinArrayBg[NMin][2];
 
  float maxBinVal;
  int maxBin;
  float bin, err;
  float s,b,es,eb;
  for(int i=0;i<nIter;i++){
    for(int q=0;q<nMin;q++){
      MinArray[q][0]=10000;
    }
    for(int x=1;x<=yHistSig->GetNbinsX();x++){
      for(int y=2;y<=(yHistSig->GetNbinsY()-1);y++){
	maxBinVal=-10000;
	for(int j=0;j<nMin;j++){//Find highest bin
	  if(MinArray[j][0]>maxBinVal){
	    maxBinVal=MinArray[j][0];
	    maxBin=j;
	  }
	}
	s=yHistSig->GetBinContent(x,y);
	b=yHistBg->GetBinContent(x,y);
	es=yHistSig->GetBinError(x,y);
	eb=yHistBg->GetBinError(x,y);
	bin=s-scalea[0]*b;
	//float a11=scalea[0];
	//bin2=s-a11*b;
	//cout << "s" << s << " " << "b"  << b << " " << bin << endl;
	//cout << bin << " " << maxBinVal << " " << maxBin << endl;
	if(bin<maxBinVal){//replce highest bin with current if lower
	  MinArray[maxBin][0]=bin;
	  // cout << maxBin << " " << maxBinVal << " " << bin << " " << MinArray[maxBin][0] << endl;
	  //divide by 0 protections
	  if(scalea[0]==0)scalea[0]=1E-10;
	  if(scalea[1]==0)scalea[1]=1E-10;
	  if(es==0)es=1E-10;
	  if(eb==0)eb=1E-10;
	  if(s==0)s=1E-10;
	  if(b==0)b=1E-10;
	  MinArray[maxBin][1]=pow(es,2)+pow(scalea[0]*b,2)*(pow(scalea[1]/scalea[0],2)+pow(eb/b,2));
	  MinArraySig[maxBin][0]=s;
	  MinArraySig[maxBin][1]=es*es;
	  MinArrayBg[maxBin][0]=b;
	  MinArrayBg[maxBin][1]=eb*eb;
	  //cout << s << " " << b << " " << bin  << " " << maxBin << endl;
	}
      }
    }
    sumSig=0, sumBg=0, esumSig=0, esumBg=0, sumMin=0, esumMin=0;
    for(int j=0;j<nMin;j++){
      //cout << MinArray[j][0] << endl;
      if(MinArray[j][0]!=10000&&MinArraySig[j][0]>1E-9&&MinArrayBg[j][0]>1E-9){
	//	cout << MinArraySig[j][0] << endl;
      sumMin+=MinArray[j][0];
      esumMin+=MinArray[j][1];
      sumSig+=MinArraySig[j][0];
      esumSig+=MinArraySig[j][1];
      sumBg+=MinArrayBg[j][0];
      esumBg+=MinArrayBg[j][1];
      }
    }
    if(sumSig==0){sumSig=1E-10; cout << "Warning: Zero sumSig" << endl;}
    if(sumBg==0){sumBg=1E-10; cout << "Warning: Zero sumBg" << endl;}
    scalea[0]=sumSig/sumBg;
    //cout << sumSig << endl;
    if(sumSig==1E-10)scalea[0]=0;
    scalea[1]=scalea[0]*pow(esumSig/sumSig/sumSig+esumBg/sumBg/sumBg,0.5);
    esumMin=pow(esumMin,0.5);
    esumSig=pow(esumSig,0.5);
    esumBg=pow(esumBg,0.5);
    cout << "Iter:  " << i << "   Min=" << sumMin << " +/- " << esumMin << "  a=" << scalea[0] << endl;  
  }
}
//--------------------------------------------------------------------------------------------------
void ZYAM2D2(TH2F *yHistSig, TH2F *yHistBg, float scalea[3], float ZYAMCent=1.,float ZYAMWidth=0.2){
  //yerr=0 no error
  //yerr=1 no mixed event
  //yerr=2 no signal
  //yerr=3 no signal and no mixed
  
  //float ZYAMCent=1.3;
  //float ZYAMWidth=0.2;//window is cent +/- width so its really a 1/2 width of the window
  float a,yea,yerr;
  //static float scalea[3];
  float ysumsig=0,ysumbg=0,yerrsig=0,yerrbg=0;
  float xcent,ycent;
  int xbins=yHistSig->GetNbinsX();
  int ybins=yHistSig->GetNbinsY();
  float yPi=3.1415926535898;
  //cout << ybins << endl;
  //cout << yHistSig->GetSum() << endl;
  //c20=new TCanvas("c20","",800,600);
  //yHistSig->Draw();
  cout << "ZYAMRegion: " << (ZYAMCent-ZYAMWidth) << "-" << (ZYAMCent+ZYAMWidth) << endl;
  for(int xi=1;xi<=xbins;xi++){
   xcent=yHistSig->GetBinCenter(xi);
    if((fabs(fabs(xcent)-ZYAMCent)<ZYAMWidth)||(fabs(fabs(xcent)-(2*yPi-ZYAMCent))<ZYAMWidth)){
    for(int yi=1;yi<=ybins;yi++){
      ysumsig+=yHistSig->GetBinContent(xi,yi);
      yerrsig+=pow(yHistSig->GetBinError(xi,yi),2);
      ysumbg+=yHistBg->GetBinContent(xi,yi);
      yerrbg+=pow(yHistBg->GetBinError(xi,yi),2);
    }
    }
  }
  cout << ysumsig << " " << yerrsig << " " << ysumbg << " " << yerrbg << endl;
  if(ysumbg==0&&ysumsig==0){
    a=0;
    yea=0;
    yerr=3;
  }
  else if(ysumbg==0){
    a=0;
    yea=0;
    yerr=1;
    
  }
  else if(ysumsig==0){
    a=0;
    yea=0;
    yerr=2;
  }
  else{
    a=ysumsig/ysumbg;
    yea=a*sqrt(yerrsig/ysumsig/ysumsig+yerrbg/ysumbg/ysumbg);
    yerr=0;
  }
  cout << a << endl;
  scalea[0]=a;
  scalea[1]=yea;
  scalea[2]=yerr;
  if(yerr) cout << "ZYA1 Error:" << yerr << endl;
  //return scalea;
}
//-------------------------------------------------------------------------------------------------------
void EffCorr(Float_t &yTPt1, Float_t &yTPt2, Float_t &yAPt1, Float_t &yAPt2, Int_t  yCent, char *yFileName, char *yEffFileName, TH1F *yPhiRaw, TH1F *yPhiCorr, TH1F *yPhiEff,  TH1F *yPhiMixRaw, TH1F *yPhiMixCorr, TH1F *yPhiMixEff, TH2F *yPhiEtaRaw, TH2F *yPhiEtaCorr, TH2F *yPhiEtaEff,  TH2F *yPhiEtaMixRaw, TH2F *yPhiEtaMixCorr, TH2F *yPhiEtaMixEff, Float_t yMPt[3], Int_t yptweighted=0, Int_t yNotDelta=0, Int_t yMethod=1){
  char name[120], name2[105];
  
  float yMPt2[2];
  
  MakeProjections(yTPt1,yTPt2,yAPt1,yAPt2,yCent,yFileName,yPhiRaw,yPhiMixRaw,yPhiEtaRaw,yPhiEtaMixRaw,yMPt,yptweighted,0,yNotDelta);
  
  sprintf(name,"yPhiRaw_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiRaw->SetName(name);
  
  sprintf(name,"yPhiMixRaw_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiMixRaw->SetName(name);
  
  sprintf(name,"yPhiEtaRaw_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaRaw->SetName(name);
  
  sprintf(name,"yPhiEtaMixRaw_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaMixRaw->SetName(name);
  
  
  //Projections from MC (from eff file)
  TH1F *yPhiMC=new TH1F("yPhiMC","",1,0,1);
  TH1F *yPhiMixMC=new TH1F("yPhiMixMC","",1,0,1);
  TH2F *yPhiEtaMC=new TH2F("yPhiEtaMC","",1,0,1,1,0,1);
  TH2F *yPhiEtaMixMC=new TH2F("yPhiEtaMixMC","",1,0,1,1,0,1);
  MakeProjections(yTPt1,yTPt2,yAPt1,yAPt2,yCent,yEffFileName,yPhiMC,yPhiMixMC,yPhiEtaMC,yPhiEtaMixMC,yMPt2,yptweighted,1,yNotDelta);
  sprintf(name,"yPhiMC_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiMC->SetName(name);
  sprintf(name,"yPhiMixMC_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiMixMC->SetName(name);
  sprintf(name,"yPhiEtaMC_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaMC->SetName(name);
  sprintf(name,"yPhiEtaMixMC_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaMixMC->SetName(name);
  
  //Efficiency Histograms
  MakeProjections(yTPt1,yTPt2,yAPt1,yAPt2,yCent,yEffFileName,yPhiEff,yPhiMixEff,yPhiEtaEff,yPhiEtaMixEff,yMPt2,yptweighted,0,yNotDelta);
  sprintf(name,"yPhiEff_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEff->SetName(name);
  sprintf(name2,"Efficiency %s",yPhiEff->GetTitle());
  yPhiEff->SetTitle(name2);
  
  sprintf(name,"yPhiMixEff_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiMixEff->SetName(name);
  sprintf(name2,"Efficiency %s",yPhiMixEff->GetTitle());
  yPhiMixEff->SetTitle(name2);
  
  sprintf(name,"yPhiEtaEff_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaEff->SetName(name);
  sprintf(name2,"Efficiency %s",yPhiEtaEff->GetTitle());
  yPhiEtaEff->SetTitle(name2);
  
  sprintf(name,"yPhiEtaMixEff_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaMixEff->SetName(name);
  sprintf(name2,"Efficiency %s",yPhiEtaMixEff->GetTitle());
  yPhiEtaMixEff->SetTitle(name2);
  
  
  //Calculate Eff
  yPhiEff->Divide(yPhiMC);
  yPhiMixEff->Divide(yPhiMixMC);
  yPhiEtaEff->Divide(yPhiEtaMC);
  yPhiEtaMixEff->Divide(yPhiEtaMixMC);
  
  // if(yMethod==0){
  // yPhiEff
  
  //Corrected Plots
  *yPhiCorr=(TH1F*)yPhiRaw->Clone();
  sprintf(name,"yPhiCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiCorr->SetName(name);
  yPhiCorr->Divide(yPhiEff);
  
  *yPhiMixCorr=(TH1F*)yPhiMixRaw->Clone();
  sprintf(name,"yPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiMixCorr->SetName(name);
  yPhiMixCorr->Divide(yPhiMixEff);
  
  *yPhiEtaCorr=(TH2F*)yPhiEtaRaw->Clone();
  sprintf(name,"yPhiEtaCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaCorr->SetName(name);
  yPhiEtaCorr->Divide(yPhiEtaEff);
  
  *yPhiEtaMixCorr=(TH2F*)yPhiEtaMixRaw->Clone();
  sprintf(name,"yPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",yTPt1,yTPt2,yAPt1,yAPt2,yCent);
  yPhiEtaMixCorr->SetName(name);
  yPhiEtaMixCorr->Divide(yPhiEtaMixEff);
  
  //Change some titles
  sprintf(name2,"Uncorrected %s",yPhiRaw->GetTitle());
  yPhiRaw->SetTitle(name2);
  sprintf(name2,"Uncorrected %s",yPhiMixRaw->GetTitle());
  yPhiMixRaw->SetTitle(name2);
  sprintf(name2,"Uncorrected %s",yPhiEtaRaw->GetTitle());
  yPhiEtaRaw->SetTitle(name2);
  sprintf(name2,"Uncorrected %s",yPhiEtaMixRaw->GetTitle());
  yPhiEtaMixRaw->SetTitle(name2);
  
  //Geometry Correction for the eta-phi plots
  if(!yNotDelta){
    int xbins=yPhiEtaCorr->GetNbinsX();
    int ybins=yPhiEtaCorr->GetNbinsY();
    float corr;
    float center;
    float scale=0;
    float sumy;
    //uses the mc for the geomery correction
    for(int yx=1;yx<=xbins;yx++){
      scale+=yPhiEtaMixMC->GetBinContent(yx,ybins/2);
      scale+=yPhiEtaMixMC->GetBinContent(yx,ybins/2+1);
    }
    scale/=2;
    for(int yy=1;yy<=ybins;yy++){
      sumy=0;
      for(int yx=1;yx<=xbins;yx++){
	sumy+=yPhiEtaMixMC->GetBinContent(yx,yy);
      }
      if(sumy==0)sumy=0.00001;
      corr=scale/sumy;
      for(int yx=1;yx<=xbins;yx++){
	yPhiEtaCorr->SetBinContent(yx,yy,corr*yPhiEtaCorr->GetBinContent(yx,yy));
	yPhiEtaCorr->SetBinError(yx,yy,corr*yPhiEtaCorr->GetBinError(yx,yy));
	yPhiEtaMixCorr->SetBinContent(yx,yy,corr*yPhiEtaMixCorr->GetBinContent(yx,yy));
	yPhiEtaMixCorr->SetBinError(yx,yy,corr*yPhiEtaMixCorr->GetBinError(yx,yy));
      }
    }
  }
  /*
  //uses equation for the geometry correct, must be changed if trigger or associtated eta cuts change
  TAxis *yAxis=yPhiEtaCorr->GetYaxis();
  for(int yx=1;yx<=xbins;yx++){
  for(int yy=1;yy<=xbins;yy++){
  center=yAxis->GetBinCenter(yy);
  if(center<0)corr=1/(1+center/1.8);
  else corr=1/(1-center/1.8);
  yPhiEtaCorr->SetBinContent(yx,yy,corr*yPhiEtaCorr->GetBinContent(yx,yy));
  yPhiEtaCorr->SetBinError(yx,yy,corr*yPhiEtaCorr->GetBinError(yx,yy));
  yPhiEtaMixCorr->SetBinContent(yx,yy,corr*yPhiEtaMixCorr->GetBinContent(yx,yy));
  yPhiEtaMixCorr->SetBinError(yx,yy,corr*yPhiEtaMixCorr->GetBinError(yx,yy));
  }
  }
  */
}
//----------------------------------------------------------------------------------------------------------
void EffCorr2(Float_t yTPt1, Float_t yTPt2, Float_t yAPt1, Float_t yAPt2, Int_t yCent, TH1F *yPhiEff, TH1F *yPhiMC, TH1F *yPhiMixEff, TH1F *yPhiMixMC, TH1F *yEtaNEff, TH1F *yEtaNMC, TH1F *yEtaNMixEff, TH1F *yEtaNMixMC, TH1F *yEtaAEff, TH1F *yEtaAMC, TH1F *yEtaAMixEff, TH1F *yEtaAMixMC, TH2F *yPhiEtaEff,  TH2F *yPhiEtaMC, TH2F *yPhiEtaMixEff, TH2F *yPhiEtaMixMC, Int_t yMethod=1){//MovingMax is special for yMethod==4
  char effname[120], effname2[105];
  if(yMethod==0){
    for(int i=1;i<=yPhiEff->GetNbinsX();i++){
      yPhiEff->SetBinContent(i,1);
      yPhiEff->SetBinError(i,0);
      yPhiMixEff->SetBinContent(i,1);
      yPhiMixEff->SetBinError(i,0);
    }
    for(int i=1;i<=yEtaNEff->GetNbinsX();i++){
      yEtaNEff->SetBinContent(i,1);
      yEtaNEff->SetBinError(i,0);
      yEtaAEff->SetBinContent(i,1);
      yEtaAEff->SetBinError(i,0);
      yEtaNMixEff->SetBinContent(i,1);
      yEtaNMixEff->SetBinError(i,0);
      yEtaAMixEff->SetBinContent(i,1);
      yEtaAMixEff->SetBinError(i,0);
    }
    for(int i=1;i<=yPhiEtaEff->GetXaxis()->GetNbins();i++){
      for(int j=1;j<=yPhiEtaEff->GetYaxis()->GetNbins();j++){
	yPhiEtaEff->SetBinContent(i,j,1);
	yPhiEtaEff->SetBinError(i,j,0);
	yPhiEtaMixEff->SetBinContent(i,j,1);
	yPhiEtaMixEff->SetBinError(i,j,0);
      }
    }
  }//yMethod==0
  else{
    //Calculate Eff (Method 1)
    yPhiEff->Divide(yPhiMC);
    yPhiMixEff->Divide(yPhiMixMC);
    yEtaNEff->Divide(yEtaNMC);
    yEtaNMixEff->Divide(yEtaNMixMC);
    yEtaAEff->Divide(yEtaAMC);
    yEtaAMixEff->Divide(yEtaAMixMC);
    yPhiEtaEff->Divide(yPhiEtaMC);
    yPhiEtaMixEff->Divide(yPhiEtaMixMC);
    
    if(yMethod==2){//Use mixed event efficiencys for triggered events
      /*
	yPhiEff=(TH1F*)yPhiMixEff->Clone();
	yPhiEff->SetName("yPhiEff");
      yEtaNEff=(TH1F*)yEtaNMixEff->Clone();
      yEtaNEff->SetName("yEtaNEff");
      yEtaAEff=(TH1F*)yEtaAMixEff->Clone();
      yEtaAEff->SetName("yEtaAEff");
      yPhiEtaEff=(TH1F*)yPhiEtaMixEff->Clone();
      yPhiEtaEff->SetName("yPhiEtaEff");
      */
      for(int i=1;i<=yPhiEff->GetNbinsX();i++){
	yPhiEff->SetBinContent(i,yPhiMixEff->GetBinContent(i));
	yPhiEff->SetBinError(i,yPhiMixEff->GetBinError(i));
      }
      for(int i=1;i<=yEtaNEff->GetNbinsX();i++){
	yEtaNEff->SetBinContent(i,yEtaNMixEff->GetBinContent(i));
	yEtaNEff->SetBinError(i,yEtaNMixEff->GetBinError(i));
	yEtaAEff->SetBinContent(i,yEtaAMixEff->GetBinContent(i));
	yEtaAEff->SetBinError(i,yEtaAMixEff->GetBinError(i));
      }
      for(int i=1;i<=yPhiEtaEff->GetXaxis()->GetNbins();i++){
	for(int j=1;j<=yPhiEtaEff->GetYaxis()->GetNbins();j++){
	  yPhiEtaEff->SetBinContent(i,j,yPhiEtaMixEff->GetBinContent(i,j));
	  yPhiEtaEff->SetBinError(i,j,yPhiEtaMixEff->GetBinError(i,j));
	}
      }
    }
    if(yMethod==3){//Sum one distribution and use it for all
      float eff=0, effE=0, effM=0, effME=0;
      for(int i=1;i<=yPhiEff->GetNbinsX();i++){
	eff+=yPhiEff->GetBinContent(i);
	effE+=pow(yPhiEff->GetBinError(i),2);
	effM+=yPhiMixEff->GetBinContent(i);
	effME+=pow(yPhiMixEff->GetBinError(i),2);
      }
      eff/=yPhiEff->GetNbinsX();
      effM/=yPhiEff->GetNbinsX();
      effE=sqrt(effE)/yPhiEff->GetNbinsX();
      effME=sqrt(effME)/yPhiEff->GetNbinsX();
      cout << "Triggered Event Eff: " << eff << " +/- " << effE << endl;
      cout << "Mixed Event Efficie: " << effM << " +/- " << effME << endl;
      
      for(int i=1;i<=yPhiEff->GetNbinsX();i++){
	yPhiEff->SetBinContent(i,eff);
	yPhiEff->SetBinError(i,effE);
	yPhiMixEff->SetBinContent(i,effM);
	yPhiMixEff->SetBinError(i,effME);
      }
      for(int i=1;i<=yEtaNEff->GetNbinsX();i++){
	yEtaNEff->SetBinContent(i,eff);
	yEtaNEff->SetBinError(i,effE);
	yEtaNMixEff->SetBinContent(i,effM);
	yEtaNMixEff->SetBinError(i,effME);
	yEtaAEff->SetBinContent(i,eff);
	yEtaAEff->SetBinError(i,effE);
	yEtaAMixEff->SetBinContent(i,effM);
	yEtaAMixEff->SetBinError(i,effME);
      }
      for(int i=1;i<=yPhiEtaEff->GetXaxis()->GetNbins();i++){
	for(int j=1;j<=yPhiEtaEff->GetYaxis()->GetNbins();j++){
	  yPhiEtaEff->SetBinContent(i,j,eff);
	  yPhiEtaEff->SetBinError(i,j,effE);
	  yPhiEtaMixEff->SetBinContent(i,j,effM);
	  yPhiEtaMixEff->SetBinError(i,j,effME);
	}
      }
    }
  }
  
  
  yPhiEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yPhiMixEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yEtaNEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yEtaNMixEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yEtaAEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yEtaAMixEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  yPhiEtaEff->GetZaxis()->SetTitle("Efficiency+Contamination");
  yPhiEtaMixEff->GetZaxis()->SetTitle("Efficiency+Contamination");
  sprintf(effname,"Efficiency Triggered %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f",yTPt1,yTPt2,yAPt1,yAPt2);
  yPhiEff->SetTitle(effname);
  sprintf(effname,"Efficiency Mixed %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f", yTPt1,yTPt2,yAPt1,yAPt2);
  yPhiMixEff->SetTitle(effname);
  
  sprintf(effname,"Near-Side Efficiency Triggered %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f",yTPt1,yTPt2,yAPt1,yAPt2);
  yEtaNEff->SetTitle(effname);
  sprintf(effname,"Near-Side Efficiency Mixed %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f", yTPt1,yTPt2,yAPt1,yAPt2);
  yEtaNMixEff->SetTitle(effname);
  
  sprintf(effname,"Away-Side Efficiency Triggered %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f",yTPt1,yTPt2,yAPt1,yAPt2);
  yEtaAEff->SetTitle(effname);
  sprintf(effname,"Away-Side Efficiency Mixed %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f", yTPt1,yTPt2,yAPt1,yAPt2);
  yEtaAMixEff->SetTitle(effname);
  
  sprintf(effname,"Efficiency Triggered %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f",yTPt1,yTPt2,yAPt1,yAPt2);
  yPhiEtaEff->SetTitle(effname);
  sprintf(effname,"Efficiency Mixed %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f",yTPt1,yTPt2,yAPt1,yAPt2);
  yPhiEtaMixEff->SetTitle(effname);
  
  
  //Corrected Plots
  //if still problems try shifting this to the other code
  // yPhiCorr->Divide(yPhiEff);
  /*
  //Change some titles
  sprintf(effname2,"Uncorrected %s",yPhiRaw->GetTitle());
  yPhiRaw->SetTitle(effname2);
  sprintf(effname2,"Uncorrected %s",yPhiMixRaw->GetTitle());
  yPhiMixRaw->SetTitle(effname2);
  sprintf(effname2,"Uncorrected %s",yPhiEtaRaw->GetTitle());
  yPhiEtaRaw->SetTitle(effname2);
  sprintf(effname2,"Uncorrected %s",yPhiEtaMixRaw->GetTitle());
  yPhiEtaMixRaw->SetTitle(effname2);
  */
}

//Efficiency Fits-----------------------------------------------------------------------
void EffFit(Float_t yAPt1,Float_t yAPt2,Int_t Cent,TList *effList, TH1F *yPhiEff, TH1F *yPhiMixEff, TH1F *yEtaNEff, TH1F *yEtaNMixEff, TH1F* yEtaAEff, TH1F *yEtaAMixEff, TH2F *yPhiEtaEff, TH2F *yPhiEtaMixEff, Int_t LSign=0, Int_t VariablePt=0){

  const int xnpt=12;
  //float xptbins[(xnpt+1)]={2.5,3,4,6,8,10,15,20,30,40,50,75,100};
  float xptbins[(xnpt+1)]={2,2.5,3,4,5,6,8,10,15,20,30,40,50};
  
  int maxPtTrig=10;//we don't want to go beyond our resolution/statistics
  //if(yAPt1>0.5)maxPtTrig=15;
  if(yAPt1>1)maxPtTrig=8;
  if(yAPt1>1.5)maxPtTrig=6;

  int minPtTrig=yAPt1;
  //set up trigger array
  int size1=0;
  for(int i=0;i<xnpt;i++)if(xptbins[i]>(yAPt1+0.0001)&&xptbins[i]<(maxPtTrig-0.0001))size1++;
  const int NTPt=size1;
  float yPi=3.1415962;
  Float_t X[NTPt];
  Float_t T0[2][NTPt];
  Float_t T1[2][NTPt];
  Float_t T2[2][NTPt];
  Float_t T3[2][NTPt];
  Float_t T4[2][NTPt];
  Float_t M0[2][NTPt];
  Float_t NT0[2][NTPt];
  Float_t NT1[2][NTPt];
  Float_t NT2[2][NTPt];
  Float_t NT3[2][NTPt];
  Float_t NT4[2][NTPt];
  Float_t NM0[2][NTPt];
  Float_t AT0[2][NTPt];
  Float_t AM0[2][NTPt];
  
  TH1F *hPhiEff=new TH1F("hPhiEff","",1,0,1);
  TH1F *hPhiMC=new TH1F("hPhiMC","",1,0,1);
  TH1F *hPhiMixEff=new TH1F("hPhiMixEff","",1,0,1);
  TH1F *hPhiMixMC=new TH1F("hPhiMixMC","",1,0,1);
  
  TH1F *hEtaNEff=new TH1F("hEtaNEff","",1,0,1);
  TH1F *hEtaNMC=new TH1F("hEtaNMC","",1,0,1);
  TH1F *hEtaNMixEff=new TH1F("hEtaNMixEff","",1,0,1);
  TH1F *hEtaNMixMC=new TH1F("hEtaNMixMC","",1,0,1);
  
  TH1F *hEtaAEff=new TH1F("hEtaAEff","",1,0,1);
  TH1F *hEtaAMC=new TH1F("hEtaAMC","",1,0,1);
  TH1F *hEtaAMixEff=new TH1F("hEtaAMixEff","",1,0,1);
  TH1F *hEtaAMixMC=new TH1F("hEtaAMixMC","",1,0,1);
  
  TH2F *hPhiEtaEff=new TH2F("hPhiEtaEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMC=new TH2F("hPhiEtaMC","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixEff=new TH2F("hPhiEtaMixEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixMC=new TH2F("hPhiEtaMixMC","",1,0,1,1,0,1);
  
  TF1 *fit1=new TF1("fit1","[0]+[1]/[2]*exp(-0.5*pow(x/[2],2))+[3]/[4]*exp(-0.5*pow((x-3.1415962)/[4],2))",-1.5,6);
  TF1 *fit2=new TF1("fit2","[0]",-2,6);
 TF1 *fit3=new TF1("fit3","[0]+[1]/[2]*exp(-0.5*pow(x/[2],2))",-2,6);
  fit1->SetParameters(1,0.1,0.1,0.1,0.1);
  fit2->SetParameter(0,1);
  fit3->SetParameters(1,0.1,0.1);
  fit1->SetParLimits(0,0,1000);
  fit1->SetParLimits(1,0,10);
  fit1->SetParLimits(2,0,10);
  fit1->SetParLimits(3,0,10);
  fit1->SetParLimits(4,0,10);
  fit3->SetParLimits(0,0,1000);
  fit3->SetParLimits(1,0,10);
  fit3->SetParLimits(2,0,10);

  float MPt[4];
  float TPt1,TPt2,APt1,APt2;
  for(int i=0;i<NTPt;i++){
    TPt1=xptbins[i];
    TPt2=xptbins[i+1];
    APt1=yAPt1;
    APt2=yAPt2;
    if(VariablePt)APt2=maxPtTrig;
    
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,MPt,0,0,0,LSign);
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiMC,hPhiMixMC,hEtaNMC,hEtaNMixMC,hEtaAMC,hEtaAMixMC,hPhiEtaMC,hPhiEtaMixMC,MPt,0,1,0,LSign);
    // EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEff,hPhiMC,hPhiMixEff,hPhiMixMC,hEtaNEff,hEtaNMC,hEtaNMixEff,hEtaNMixMC,hEtaAEff,hEtaAMC,hEtaAMixEff,hEtaAMixMC,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixEff,hPhiEtaMixMC,1);
    EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEff,hPhiMC,hPhiMixEff,hPhiMixMC,hEtaNEff,hEtaNMC,hEtaNMixEff,hEtaNMixMC,hEtaAEff,hEtaAMC,hEtaAMixEff,hEtaAMixMC,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixEff,hPhiEtaMixMC,1);
    X[i]=(TPt1+TPt2)*0.5;
    
    hPhiEff->Fit("fit1");
    T0[0][i]=fit1->GetParameter(0);
    T0[1][i]=fit1->GetParError(0);
    T1[0][i]=fit1->GetParameter(1);
    T1[1][i]=fit1->GetParError(1);
    T2[0][i]=fit1->GetParameter(2);
    T2[1][i]=fit1->GetParError(2);
    T3[0][i]=fit1->GetParameter(3);
    T3[1][i]=fit1->GetParError(3);
    T4[0][i]=fit1->GetParameter(4);
    T4[1][i]=fit1->GetParError(4);

    hPhiMixEff->Fit("fit2");
    M0[0][i]=fit2->GetParameter(0);
    M0[1][i]=fit2->GetParError(0);
    
    hEtaNEff->Fit("fit3");
    NT0[0][i]=fit3->GetParameter(0);
    NT0[1][i]=fit3->GetParError(0);
    NT1[0][i]=fit3->GetParameter(1);
    NT1[1][i]=fit3->GetParError(1);
    NT2[0][i]=fit3->GetParameter(2);
    NT2[1][i]=fit3->GetParError(2);
 
    
    hEtaNMixEff->Fit("fit2");
    NM0[0][i]=fit2->GetParameter(0);
    NM0[1][i]=fit2->GetParError(0);
    
    hEtaAEff->Fit("fit2");
    AT0[0][i]=fit2->GetParameter(0);
    AT0[1][i]=fit2->GetParError(0);
    
    hEtaNMixEff->Fit("fit2");
    AM0[0][i]=fit2->GetParameter(0);
    AM0[1][i]=fit2->GetParError(0);
  }
  
  gT0=new TGraphErrors(NTPt,X,T0[0],0,T0[1]);
  gT1=new TGraphErrors(NTPt,X,T1[0],0,T1[1]);
  gT2=new TGraphErrors(NTPt,X,T2[0],0,T2[1]);
  gT3=new TGraphErrors(NTPt,X,T3[0],0,T3[1]);
  gT4=new TGraphErrors(NTPt,X,T4[0],0,T4[1]);

  gM0=new TGraphErrors(NTPt,X,M0[0],0,M0[1]);
  
  gNT0=new TGraphErrors(NTPt,X,NT0[0],0,NT0[1]);
  gNT1=new TGraphErrors(NTPt,X,NT1[0],0,NT1[1]);
  gNT2=new TGraphErrors(NTPt,X,NT2[0],0,NT2[1]);

  gNM0=new TGraphErrors(NTPt,X,NM0[0],0,NM0[1]);

  gAT0=new TGraphErrors(NTPt,X,AT0[0],0,AT0[1]);
  gAM0=new TGraphErrors(NTPt,X,AM0[0],0,AM0[1]);
  
  float pT0[2],pT1[2],pT2[2],pT3[2],pT4[2];
  float pM0[2];
  float pNT0[2],pNT1[2],pNT2[2],pNT3[2],pNT4[2];
  float pNM0[2];
  float pAT0[2];
  float pAM0[2];
  
  gT0->Fit("fit2");
  pT0[0]=fit2->GetParameter(0);
  pT0[1]=fit2->GetParError(0);
  gT1->Fit("fit2");
  pT1[0]=fit2->GetParameter(0);
  pT1[1]=fit2->GetParError(0);
  gT2->Fit("fit2");
  pT2[0]=fit2->GetParameter(0);
  pT2[1]=fit2->GetParError(0);
  gT3->Fit("fit2");
  pT3[0]=fit2->GetParameter(0);
  pT3[1]=fit2->GetParError(0);
  gT4->Fit("fit2");
  pT4[0]=fit2->GetParameter(0);
  pT4[1]=fit2->GetParError(0);
  
  gNT0->Fit("fit2");
  pNT0[0]=fit2->GetParameter(0);
  pNT0[1]=fit2->GetParError(0);
  gNT1->Fit("fit2");
  pNT1[0]=fit2->GetParameter(0);
  pNT1[1]=fit2->GetParError(0);
  gNT2->Fit("fit2");
  pNT2[0]=fit2->GetParameter(0);
  pNT2[1]=fit2->GetParError(0);
  
  gM0->Fit("fit2");
  pM0[0]=fit2->GetParameter(0);
  pM0[1]=fit2->GetParError(0);
  
  gNM0->Fit("fit2");
  pNM0[0]=fit2->GetParameter(0);
  pNM0[1]=fit2->GetParError(0);
  
  gAM0->Fit("fit2");
  pAM0[0]=fit2->GetParameter(0);
  pAM0[1]=fit2->GetParError(0);
  
  gAT0->Fit("fit2");
  pAT0[0]=fit2->GetParameter(0);
  pAT0[1]=fit2->GetParError(0);
  
  float eff,gaus,err,cent;
  float gausA,errA;
  for(int x=1;x<=yPhiEff->GetNbinsX();x++){
    cent=yPhiEff->GetBinCenter(x);
    if(pT2[0]==0)cout << "Error pT2[0]" << endl;
    gaus=pT1[0]/pT2[0]*exp(-0.5*pow(cent/pT2[0],2));
    gausA=pT3[0]/pT4[0]*exp(-0.5*pow((cent-yPi)/pT4[0],2));
    eff=gaus+gausA+pT0[0];
    if(pT1[0]==0) cout << "Error pT1[0]" << endl;
    if(gaus==0)gaus+=1E-10;//divide by 0 protection
    if(gausA==0)gausA+=1E-10;
    err=pow(pow(pT0[1],2)+pow(pT1[1]*gaus/pT1[0],2)+pow(pT2[1]*gaus/pT2[0]*(pow(cent/pT2[0],2)-1),2)+pow(pT3[1]*gausA/pT3[0],2)+pow(pT4[1]*gausA/pT4[0]*(pow((cent-yPi)/pT4[0],2)-1),2),0.5);
    yPhiEff->SetBinContent(x,eff);
    yPhiEff->SetBinError(x,err);
    
    yPhiMixEff->SetBinContent(x,pM0[0]);
    yPhiMixEff->SetBinError(x,pM0[1]);
  }
  for(int x=1;x<=yEtaNEff->GetNbinsX();x++){
    cent=yEtaNEff->GetBinCenter(x);
    if(pNT2[0]==0)cout << "Error pNT2[0]" << endl;
    gaus=pNT1[0]/pNT2[0]*exp(-0.5*pow(cent/pNT2[0],2));
    eff=gaus+pNT0[0];
    if(pNT1[0]==0) cout << "Error pNT1[0]" << endl;
    err=pow(pow(pNT0[1],2)+pow(pNT1[1]*gaus/pNT1[0],2)+pow(pNT2[1]*gaus/pNT2[0]*(pow(cent/pNT2[0],2)-1),2),0.5);
    yEtaNEff->SetBinContent(x,eff);
    yEtaNEff->SetBinError(x,err);

    yEtaNMixEff->SetBinContent(x,pNM0[0]);
    yEtaNMixEff->SetBinError(x,pNM0[1]);

    yEtaAEff->SetBinContent(x,pAT0[0]);
    yEtaAEff->SetBinError(x,pAT0[1]);
    
    yEtaAMixEff->SetBinContent(x,pAM0[0]);
    yEtaAMixEff->SetBinError(x,pAM0[1]);
  }

  float gaus2,eff2,err2,cent2;
  float eff3,err3;
  for(int x=1;x<=yPhiEtaEff->GetXaxis()->GetNbins();x++){
    cent=yPhiEtaEff->GetXaxis()->GetBinCenter(x);
    for(int y=1;y<=yPhiEtaEff->GetYaxis()->GetNbins();y++){
      if(fabs(cent)<1.5||fabs(cent+2*yPi)<1.5||fabs(cent-2*yPi)<1.5){
      gaus=pT1[0]/pT2[0]*exp(-0.5*pow(cent/pT2[0],2));
      eff=gaus;
      if(eff==0)eff+=1E-10;//divide by 0 protection
      err=pow(pow(pT1[1]*gaus/pT1[0],2)+pow(pT2[1]*gaus/pT2[0]*(pow(cent/pT2[0],2)-1),2),0.5); 
      cent2=yPhiEtaEff->GetYaxis()->GetBinCenter(y);
      gaus2=pNT1[0]/pNT2[0]*exp(-0.5*pow(cent2/pNT2[0],2));
      eff2=gaus2;
      if(eff2==0)eff2+=1E-10;//divide by 0 protection
      err2=pow(pow(pNT1[1]*gaus2/pNT1[0],2)+pow(pNT2[1]*gaus2/pNT2[0]*(pow(cent2/pNT2[0],2)-1),2),0.5);
      // eff3=pow(eff*eff2,0.5);
      eff3=pow(gaus*gaus2,0.5)+pow((pT0[0])*pNT0[0],0.5);
      // eff=pow(pT0[0]*pNT0[0],0.5);
      err3=pow(gaus*gaus2*(pow(err/eff,2)+pow(err2/eff2,2))+pT0[0]*pNT0[0]*(pow(pT0[1]/pT0[0],2)+pow(pNT0[1]/pNT0[0],2)),0.5);
    
      }
      else{
	gausA=pT3[0]/pT4[0]*exp(-0.5*pow((cent-yPi)/pT4[0],2));
	if(gausA==0)gausA+=1E-10;
	eff3=pow((pT0[0]+gausA)*pAT0[0],0.5);
	//err3=pow(pow(pT0[1]/pT0[0],2)+pow(pAT0[1]/pAT0[0],2),0.5);
	err3=eff3*pow((pow(pT0[1],2)+pow(pT3[1]*gausA/pT3[0],2)+pow(pT4[1]*gausA/pT4[0]*(pow((cent-yPi)/pT4[0],2)-1),2))/pow(pT0[0]+gausA,2)+pow(pAT0[1]/pAT0[1],2),0.5);
      }
       yPhiEtaEff->SetBinContent(x,y,eff3);
      yPhiEtaEff->SetBinError(x,y,err3);
      
      if(fabs(cent)<1.5){
	eff3=pow(pM0[0]*pNM0[0],0.5);
	err3=pow(pow(pM0[1]/pM0[0],2)+pow(pNM0[1]/pNM0[0],2),0.5);
      }
      else{
	eff3=pow(pM0[0]*pAM0[0],0.5);
	err3=pow(pow(pM0[1]/pM0[0],2)+pow(pAM0[1]/pAM0[0],2),0.5);
      }
      yPhiEtaMixEff->SetBinContent(x,y,eff3);
      yPhiEtaMixEff->SetBinError(x,y,err3);
    }
  }
  }

 //Geometry Correction for the eta-phi plots
 //\-------------------------------------------------------------------
void GeoCorr(TH1F *yEtaNCorr, TH1F *yEtaNMixCorr, TH1F *yEtaNMixMC, TH1F *yEtaACorr, TH1F *yEtaAMixCorr, TH1F *yEtaAMixMC, TH2F *yPhiEtaCorr, TH2F *yPhiEtaMixCorr, TH2F *yPhiEtaMixMC){
  
  int xbins=yEtaNCorr->GetNbinsX();
   
  float corr;
  float center;
  float scale=0,scaleN,scaleA;
  float sumy;
  float escaleN,escaleA, ecorr;
  float con,err;
  scaleN=(yEtaNMix->GetBinContent(xbins/2)+yEtaNMix->GetBinContent(xbins/2+1))/2.;
  escaleN=pow(yEtaNMix->GetBinError(xbins/2),2)+pow(yEtaNMix->GetBinError(xbins/2+1),2)/4.;
  scaleA=(yEtaAMix->GetBinContent(xbins/2)+yEtaAMix->GetBinContent(xbins/2+1))/2.;
  escaleA=pow(yEtaAMix->GetBinError(xbins/2),2)+pow(yEtaAMix->GetBinError(xbins/2+1),2)/4.;
  if(scaleN==0)scaleN=1E-5;
  if(scaleA==0)scaleA=1E-5;
  for(int yx=1;yx<=xbins;yx++){
    con=yEtaNMixMC->GetBinContent(yx);
    if(con==0){con=1E-5;corr=1E-5;}//divide by 0 protection
    else corr=scaleN/con;
    ecorr=corr*corr*(escaleN/scaleN/scaleN+pow(yEtaNMixMC->GetBinError(yx),2)/con/con);
    con=yEtaNCorr->GetBinContent(yx);
    err=pow(yEtaNCorr->GetBinError(yx),2);
    yEtaNCorr->SetBinContent(yx,corr*con);
    if(con==0)con=1E-5;//Divide by 0 protection
    if(corr==0)corr=1E-5;
    yEtaNCorr->SetBinError(yx,corr*con*pow(ecorr/corr/corr+err/con/con,0.5));
    con=yEtaNMixCorr->GetBinContent(yx);
    err=pow(yEtaNMixCorr->GetBinError(yx),2);
    yEtaNMixCorr->SetBinContent(yx,corr*con);
    if(con==0)con=1E-5;
    if(corr==0)corr=1E-5;
    yEtaNMixCorr->SetBinError(yx,corr*con*pow(ecorr/corr/corr+err/con/con,0.5));

    con=yEtaAMixMC->GetBinContent(yx);
    if(con==0){con=1E-5;corr=1E-5;}
    else corr=scaleA/con;
    ecorr=corr*corr*(escaleA/scaleA/scaleA+pow(yEtaAMixMC->GetBinError(yx),2)/con/con);
    con=yEtaACorr->GetBinContent(yx);
    err=pow(yEtaACorr->GetBinError(yx),2);
    yEtaACorr->SetBinContent(yx,corr*con);
    if(con==0)con=1E-5;
    if(corr==0)corr=1E-5;
    yEtaACorr->SetBinError(yx,corr*con*pow(ecorr/corr/corr+err/con/con,0.5));
    con=yEtaAMixCorr->GetBinContent(yx);
    err=pow(yEtaAMixCorr->GetBinError(yx),2);
    yEtaAMixCorr->SetBinContent(yx,corr*con);
    if(con==0)con=1E-5;
    if(corr==0)corr=1E-5;
    yEtaAMixCorr->SetBinError(yx,corr*con*pow(ecorr/corr/corr+err/con/con,0.5));
  }//end eta loop
  
  xbins=yPhiEtaCorr->GetNbinsX();
  int ybins=yPhiEtaCorr->GetNbinsY();
  scale=0;
  //uses the mc for the geomery correction
  for(int yx=1;yx<=xbins;yx++){
    scale+=yPhiEtaMixMC->GetBinContent(yx,ybins/2);
    scale+=yPhiEtaMixMC->GetBinContent(yx,ybins/2+1);
  }
  scale/=2.;
  for(int yy=1;yy<=ybins;yy++){
    sumy=0;
    for(int yx=1;yx<=xbins;yx++){
      sumy+=yPhiEtaMixMC->GetBinContent(yx,yy);
    }
    if(sumy==0)sumy=0.00001;
    corr=scale/sumy;
    for(int yx=1;yx<=xbins;yx++){
      yPhiEtaCorr->SetBinContent(yx,yy,corr*yPhiEtaCorr->GetBinContent(yx,yy));
      yPhiEtaCorr->SetBinError(yx,yy,corr*yPhiEtaCorr->GetBinError(yx,yy));
      yPhiEtaMixCorr->SetBinContent(yx,yy,corr*yPhiEtaMixCorr->GetBinContent(yx,yy));
      yPhiEtaMixCorr->SetBinError(yx,yy,corr*yPhiEtaMixCorr->GetBinError(yx,yy));
    }
  }
}
//-----------------------------------------------------------------------------------------------
void GeoCorr2(TH1F *yEtaNCorr, TH1F *yEtaNMixCorr, TH1F *yEtaACorr, TH1F *yEtaAMixCorr, TH2F *yPhiEtaCorr, TH2F *yPhiEtaMixCorr){
  int nbins=yEtaNCorr->GetNbinsX();
  int nbinsPhi=yPhiEtaCorr->GetNbinsX();
  int nbinsEta=yPhiEtaCorr->GetNbinsY();
  float EtaCut=yEtaNCorr->GetBinCenter(nbins)+yEtaNCorr->GetBinWidth(nbins)/2.;
  float corr;
  float cent;
  for(int i=1;i<=nbins;i++){
    cent=yEtaNCorr->GetBinCenter(i);
    corr=1-fabs(cent)/EtaCut;
    yEtaNCorr->SetBinContent(i,yEtaNCorr->GetBinContent(i)/corr);
    yEtaNCorr->SetBinError(i,yEtaNCorr->GetBinError(i)/corr);
    yEtaNMixCorr->SetBinContent(i,yEtaNMixCorr->GetBinContent(i)/corr);
    yEtaNMixCorr->SetBinError(i,yEtaNMixCorr->GetBinError(i)/corr);
    yEtaACorr->SetBinContent(i,yEtaACorr->GetBinContent(i)/corr);
    yEtaACorr->SetBinError(i,yEtaACorr->GetBinError(i)/corr);
    yEtaAMixCorr->SetBinContent(i,yEtaAMixCorr->GetBinContent(i)/corr);
    yEtaAMixCorr->SetBinError(i,yEtaAMixCorr->GetBinError(i)/corr);
  }
  for(int i=1;i<=nbinsEta;i++){
    cent=yPhiEtaCorr->GetYaxis()->GetBinCenter(i);
    corr=1-fabs(cent)/EtaCut;
    for(int j=1;j<=nbinsPhi;j++){
      yPhiEtaCorr->SetBinContent(j,i,yPhiEtaCorr->GetBinContent(j,i)/corr);
      yPhiEtaCorr->SetBinError(j,i,yPhiEtaCorr->GetBinError(j,i)/corr);
      yPhiEtaMixCorr->SetBinContent(j,i,yPhiEtaMixCorr->GetBinContent(j,i)/corr);
      yPhiEtaMixCorr->SetBinError(j,i,yPhiEtaMixCorr->GetBinError(j,i)/corr);
    }
  }
}
//---------------------------------------------------------------------------------
void Load3Particle(Float_t &xTPt1, Float_t &xTPt2, Float_t &xAPt1, Float_t &xAPt2,Int_t xCent, TList *xList, TH2F *xPhiPhiRaw,TH2F *xPhiPhiSS,TH2F *xPhiPhiMix, TH2F *xEtaEtaRaw, TH2F *xEtaEtaSS, TH2F *xEtaEtaMix,Int_t xMC=0, Int_t xSign=0){
  const int xnAPt=7;
  float xAPt31[xnAPt]={0.5,1.0,1.5,2.0,3,4,1};
  float xAPt32[xnAPt]={1.0,1.5,2.0,3.0,4,5,2};
  const int xnTPt=12;
  float xTPt[(xnTPt+1)]={2,2.5,3,4,5,6,8,10,15,20,30,40,50};//any element in this array should also been in the one 2 lines down
 float xPi=3.1415962;
  int APtBin=-1;
  int TPtBin1=-1;
  int TPtBin2=-1;
  for(int i=0;i<xnAPt;i++){
    if(fabs(xAPt1-xAPt31[i])<0.1&&fabs(xAPt2-xAPt32[i])<0.1){
      xAPt1=xAPt31[i];
      xAPt2=xAPt32[i];
      APtBin=i;
    }
  }
  if(APtBin==-1){
    cout << "Invalid Associated Pt: " << xAPt1 << " " << xAPt2 << endl;
    cout << "Valid values are:" << endl;
    for(int i=0;i<xnAPt;i++){
      cout << "APt1: " << xAPt31[i] << "  APt2: " << xAPt32[i] << endl;
    }
    break;
  }
  for(int i=0;i<xnTPt;i++){
    if(fabs(xTPt1-xTPt[i])<0.1){
     xTPt1=xTPt[i];
     TPtBin1=i;
    }
    if(fabs(xTPt2-xTPt[i+1])<0.1){
      xTPt2=xTPt[i+1];
      TPtBin2=i;
    }
  }
  if(TPtBin1==-1||TPtBin2==-1){
    cout << "Invalid Trigger Pt: " << xTPt1 << " " << xTPt2 << end;
     cout << "Valid values are:" << endl;
    for(int i=0;i<=xnTPt;i++){
      cout <<  xTPt[i] << "  ";
    }
    cout << endl;
    break;
  }
  char xname[100];
  char xtit[200];
 char *cmc1[2]={"","_MC"};
  sprintf(xname,"fHistNTrigger_C%d%s",xCent,cmc1[xMC]);
  TH1F *xNTrig=(TH1F*)xList->FindObject(xname);
  sprintf(xname,"fHistNMix_C%d%s",xCent,cmc1[xMC]);
  TH1F *xNMix=(TH1F*)xList->FindObject(xname);
  int ntriggers=xNTrig->GetBinContent(1);
  int nmix=xNMix->GetBinContent(1);
  // c100=new TCanvas("c100");
  // xNTrig->Draw();
  TH2F *xtemphistP;
  TH2F *xtemphistPSS;
  TH2F *xtemphistPMix;
  TH2F *xtemphistE;
  TH2F *xtemphistESS;
  TH2F *xtemphistEMix;
  //char *cmc1[2]={"","_MC"};
  char *sign1[3]={"","_LS","_ULS"};
  char *sign2[3]={""," LikeSign"," UnlikeSign"};
  int nbins;
  float min,max;
  float conP, conPSS, conPMix, conE, conESS, conEMix;
  float errP, errPSS, errPMix, errE, errESS, errEMix;
  float xntrig=0;
  float xnmix=0;
  

  for(int xi=TPtBin1;xi<=TPtBin2;xi++){//Do for all trigger pt ranges we are adding
    xntrig+=xNTrig->GetBinContent(xi+1);
    //cout << xi << " " << xNTrig->GetBinContent(xi+1) << endl;
    xnmix+=xNMix->GetBinContent(xi+1);
    sprintf(xname,"fHistDeltaPhiPhi_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistP=(TH2F*)xList->FindObject(xname);

    sprintf(xname,"fHistDeltaPhiPhiSS_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistPSS=(TH2F*)xList->FindObject(xname);

    sprintf(xname,"fHistDeltaPhiPhiMix_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistPMix=(TH2F*)xList->FindObject(xname);
   
    sprintf(xname,"fHistDeltaEtaEta_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistE=(TH2F*)xList->FindObject(xname);

    sprintf(xname,"fHistDeltaEtaEtaSS_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistESS=(TH2F*)xList->FindObject(xname);

    sprintf(xname,"fHistDeltaEtaEtaMix_P%dp%d_C%d%s%s",xi,APtBin,xCent,cmc1[xMC],sign1[xSign]);
    xtemphistEMix=(TH2F*)xList->FindObject(xname);
  
    if(xi==TPtBin1){
      nbins=xtemphistP->GetNbinsX();
      min=xtemphistP->GetBinCenter(1)-xtemphistP->GetBinWidth(1)/2.;
      max=xtemphistP->GetBinCenter(nbins)+xtemphistP->GetBinWidth(nbins)/2.;
      
      sprintf(xtit,"#Delta#phi-#Delta#phi Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
      *xPhiPhiRaw=new TH2F("dPhiPhiRaw",xtit,nbins,min,max,nbins,min,max);
      xPhiPhiRaw->Sumw2();
      xPhiPhiRaw->GetXaxis()->SetTitle("#Delta#phi");
      xPhiPhiRaw->GetYaxis()->SetTitle("#Delta#phi");
      xPhiPhiRaw->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#phid#Delta#phi}");
      SetTitles2D(xPhiPhiRaw);
      
      sprintf(xtit,"Soft-Soft #Delta#phi-#Delta#phi Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
      *xPhiPhiSS=new TH2F("dPhiPhiSS",xtit,nbins,min,max,nbins,min,max);
      xPhiPhiSS->Sumw2();
      xPhiPhiSS->GetXaxis()->SetTitle("#Delta#phi");
      xPhiPhiSS->GetYaxis()->SetTitle("#Delta#phi");
      xPhiPhiSS->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#phid#Delta#phi}");
      SetTitles2D(xPhiPhiSS);

       sprintf(xtit,"Mixed #Delta#phi-#Delta#phi Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
      *xPhiPhiMix=new TH2F("dPhiPhiMix",xtit,nbins,min,max,nbins,min,max);
      xPhiPhiMix->Sumw2();
      xPhiPhiMix->GetXaxis()->SetTitle("#Delta#phi");
      xPhiPhiMix->GetYaxis()->SetTitle("#Delta#phi");
      xPhiPhiMix->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#phid#Delta#phi}");
      SetTitles2D(xPhiPhiMix);
      nbins=xtemphistE->GetNbinsX();
      min=xtemphistE->GetBinCenter(1)-xtemphistE->GetBinWidth(1)/2.;
      max=xtemphistE->GetBinCenter(nbins)+xtemphistE->GetBinWidth(nbins)/2.;
      sprintf(xtit,"#Delta#eta-#Delta#eta Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
      *xEtaEtaRaw=new TH2F("dEtaEtaRaw",xtit,nbins,min,max,nbins,min,max);
      xEtaEtaRaw->Sumw2();
      xEtaEtaRaw->GetXaxis()->SetTitle("#Delta#eta");
      xEtaEtaRaw->GetYaxis()->SetTitle("#Delta#eta");
      xEtaEtaRaw->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#etad#Delta#eta}");
      SetTitles2D(xEtaEtaRaw);
      
       sprintf(xtit,"Soft-Soft #Delta#eta-#Delta#eta Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
       *xEtaEtaSS=new TH2F("dEtaEtaSS",xtit,nbins,min,max,nbins,min,max);
      xEtaEtaSS->Sumw2();
      xEtaEtaSS->GetXaxis()->SetTitle("#Delta#eta");
      xEtaEtaSS->GetYaxis()->SetTitle("#Delta#eta");
      xEtaEtaSS->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#etad#Delta#eta}");
      SetTitles2D(xEtaEtaSS);

       sprintf(xtit,"Mixed #Delta#eta-#Delta#eta Distribution for  %3.1f<p_{T}^{Trig}<%3.1f %2.2f<p_{T}^{Assoc}<%2.2f%s",xTPt1,xTPt2,xAPt1,xAPt2,sign2[xSign]);
      *xEtaEtaMix=new TH2F("dEtaEtaMix",xtit,nbins,min,max,nbins,min,max);
      xEtaEtaMix->Sumw2();
      xEtaEtaMix->GetXaxis()->SetTitle("#Delta#eta");
      xEtaEtaMix->GetYaxis()->SetTitle("#Delta#eta");
      xEtaEtaMix->GetZaxis()->SetTitle("#frac{1}{N_{Trig}}#frac{d^{2}N_{pairs}}{d#Delta#etad#Delta#eta}");
      SetTitles2D(xEtaEtaMix);

    }
    for(int x=1;x<=xtemphistP->GetNbinsX();x++){
      for(int y=1;y<=xtemphistP->GetNbinsY();y++){
	if(xi==TPtBin1){
	  conP=0;
	  errP=0;
	  conPSS=0;
	  errPSS=0;
	  conPMix=0;
	  errPMix=0;
	}
	else{
	  conP=xPhiPhiRaw->GetBinContent(x,y);
	  errP=xPhiPhiRaw->GetBinError(x,y);
	  conPSS=xPhiPhiSS->GetBinContent(x,y);
	  errPSS=xPhiPhiSS->GetBinError(x,y);
	  conPMix=xPhiPhiMix->GetBinContent(x,y);
	  errPMix=xPhiPhiMix->GetBinError(x,y);
	}
	xPhiPhiRaw->SetBinContent(x,y,conP+xtemphistP->GetBinContent(x,y));
	xPhiPhiRaw->SetBinError(x,y,pow(errP*errP+pow(xtemphistP->GetBinError(x,y),2),0.5));
	xPhiPhiSS->SetBinContent(x,y,conPSS+xtemphistPSS->GetBinContent(x,y));
	xPhiPhiSS->SetBinError(x,y,pow(errPSS*errPSS+pow(xtemphistPSS->GetBinError(x,y),2),0.5));
	xPhiPhiMix->SetBinContent(x,y,conPMix+xtemphistPMix->GetBinContent(x,y));
	xPhiPhiMix->SetBinError(x,y,pow(errPMix*errPMix+pow(xtemphistPMix->GetBinError(x,y),2),0.5));
      }
    }
    
    for(int x=1;x<=xtemphistE->GetNbinsX();x++){
      for(int y=1;y<=xtemphistE->GetNbinsY();y++){
	if(xi==TPtBin1){
	  conE=0;
	  errE=0;
	  conESS=0;
	  errESS=0;
	  conEMix=0;
	  errEMix=0;
	}
	else{
	  conE=xEtaEtaRaw->GetBinContent(x,y);
	  errE=xEtaEtaRaw->GetBinError(x,y);
	  conESS=xEtaEtaSS->GetBinContent(x,y);
	  errESS=xEtaEtaSS->GetBinError(x,y);
	  conEMix=xEtaEtaMix->GetBinContent(x,y);
	  errEMix=xEtaEtaMix->GetBinError(x,y);
	}
	xEtaEtaRaw->SetBinContent(x,y,conE+xtemphistE->GetBinContent(x,y));
	xEtaEtaRaw->SetBinError(x,y,pow(errE*errE+pow(xtemphistE->GetBinError(x,y),2),0.5));
	xEtaEtaSS->SetBinContent(x,y,conESS+xtemphistESS->GetBinContent(x,y));
	xEtaEtaSS->SetBinError(x,y,pow(errESS*errESS+pow(xtemphistESS->GetBinError(x,y),2),0.5));
	xEtaEtaMix->SetBinContent(x,y,conEMix+xtemphistEMix->GetBinContent(x,y));
	xEtaEtaMix->SetBinError(x,y,pow(errEMix*errEMix+pow(xtemphistEMix->GetBinError(x,y),2),0.5));
      
      }
    }
  }

//Normalize the histograms
  //cout << xPi << " " << max << " " << xntrig << " " << xnmix << endl;
  float phibins=xPhiPhiRaw->GetNbinsX();
  float etabins=xEtaEtaRaw->GetNbinsX();
  float scalephi=pow(phibins/(2*xPi),2);
  float scaleeta=pow(etabins/(max-min),2);//max still set for eta from making histograms
  xPhiPhiRaw->Scale(scalephi/xntrig);
  xPhiPhiSS->Scale(scalephi/xnmix);
  xPhiPhiMix->Scale(scalephi/(2.*xnmix));
  xEtaEtaRaw->Scale(scaleeta/xntrig);
  xEtaEtaSS->Scale(scaleeta/xnmix);
  xEtaEtaMix->Scale(scaleeta/(2.*xnmix));
}

//------------------------------------------------------------------------------------
void EffCorr3Part(Float_t TPt1, Float_t TPt2, Float_t APt1, Float_t APt2, Int_t Cent, TH2F *yPhiPhiEff, TH2F *yPhiPhiMC, TH2F *yPhiPhiSSEff, TH2F *yPhiPhiSSMC, TH2F *yPhiPhiMixEff, TH2F *yPhiPhiMixMC, TH2F *yEtaEtaEff, TH2F *yEtaEtaMC, TH2F *yEtaEtaSSEff, TH2F *yEtaEtaSSMC, TH2F *yEtaEtaMixEff, TH2F *yEtaEtaMixMC, Int_t yMethod){
  //Right now only putting in for case where efficiency is taken care of before hand.  I will add others if used
  if(yMethod==0){
    int nBinPhi=yPhiPhiEff->GetNbinsX();
    int nBinEta=yEtaEtaEff->GetNbinsX();
    for(int x=1;x<=nBinPhi;x++){
      for(int y=1;y<=nBinPhi;y++){
	yPhiPhiEff->SetBinContent(x,y,1);
	yPhiPhiEff->SetBinError(x,y,0);
	yPhiPhiSSEff->SetBinContent(x,y,1);
	yPhiPhiSSEff->SetBinError(x,y,0);
	yPhiPhiMixEff->SetBinContent(x,y,1);
	yPhiPhiMixEff->SetBinError(x,y,0);
      }
    }
  for(int x=1;x<=nBinEta;x++){
      for(int y=1;y<=nBinEta;y++){
	yEtaEtaEff->SetBinContent(x,y,1);
	yEtaEtaEff->SetBinError(x,y,0);
	yEtaEtaSSEff->SetBinContent(x,y,1);
	yEtaEtaSSEff->SetBinError(x,y,0);
	yEtaEtaMixEff->SetBinContent(x,y,1);
	yEtaEtaMixEff->SetBinError(x,y,0);
      }
    }
  }
  else{
    Printf("Method needs to be inserted into EffCorr3Part");
    break;
  }
  

}

//-------------------------------------------------------------------------------------------------------
void GeoCorr3Part2(TH2F *yEtaEtaCorr, TH2F *yEtaEtaSSCorr, TH2F *yEtaEtaMixCorr, TH2F *yEtaEtaHSCorr){//using equation to correct
  int bin=yEtaEtaCorr->GetNbinsX();
  float centx, centy;
  float corr,corrx,corry,corrz;
  float EtaCut=yEtaEtaCorr->GetBinCenter(bin)+yEtaEtaCorr->GetBinWidth(bin)/2.;
  cout << "EtaCut " << EtaCut << endl;
  for(int x=1;x<=bin;x++){
    centx=yEtaEtaCorr->GetXaxis()->GetBinCenter(x);
    corrx=1-fabs(centx)/EtaCut;
    for(int y=1;y<=bin;y++){
      centy=yEtaEtaCorr->GetYaxis()->GetBinCenter(y);
      corry=1-fabs(centy)/EtaCut;
      corrz=1-fabs(centx-centy)/EtaCut;
      if(fabs(centx-centy)>EtaCut)corrz=0;
      corr=corrz*pow(corrx*corry,1/2.);
      if(corr==0)corr=1;
      corr=1;
      yEtaEtaCorr->SetBinContent(x,y,yEtaEtaCorr->GetBinContent(x,y)/corr);
      yEtaEtaCorr->SetBinError(x,y,yEtaEtaCorr->GetBinError(x,y)/corr);
      yEtaEtaMixCorr->SetBinContent(x,y,yEtaEtaMixCorr->GetBinContent(x,y)/corr);
      yEtaEtaMixCorr->SetBinError(x,y,yEtaEtaMixCorr->GetBinError(x,y)/corr);
      yEtaEtaSSCorr->SetBinContent(x,y,yEtaEtaSSCorr->GetBinContent(x,y)/corr);
      yEtaEtaSSCorr->SetBinError(x,y,yEtaEtaSSCorr->GetBinError(x,y)/corr);
      yEtaEtaHSCorr->SetBinContent(x,y,yEtaEtaHSCorr->GetBinContent(x,y)/corr);
      yEtaEtaHSCorr->SetBinError(x,y,yEtaEtaHSCorr->GetBinError(x,y)/corr);
    }
  }
}
//--------------------------------------------------------------------------------------------------------------------
void ZYAM3(TH1F *yPhi, TH1F *yPhiMix, TH2F *yPhiPhi, TH2F *yPhiPhiSS, TH2F *yPhiPhiMix, Float_t a[3], Float_t b, Int_t nMin=10, Int_t nIter=10){//xxx scale the mixed event to the SS before it comes in
  const Int_t NMin=100;
  
  // float sumSig2, sumMix2, sumSig3, sumSS3, sumMix3;
  //float esumSig2, esumMix2, esumSig3, esumSS3, esumMix3;
  //float sumMin, esumMin;
  float itersize=0.1;
  float oldmin=1000;
  a[0]=1;
  a[1]=0;
  float MinArray[NMin][2];
  float Phi[NMin][2];
  float SS[NMin][2];
  float Mix[NMin][2];
  float PhiX[NMin][2];
  float PhiY[NMin][2];
   float MixX[NMin][2];
  float MixY[NMin][2];
  int nBins2=yPhi->GetNbinsX();
  int nBins3=yPhiPhi->GetNbinsX();
  int ReBin=nBins2/nBins3;
    if(fabs(1.0*nBins2/nBins3-ReBin)>0.01)Printf("Please use a divisor of the 2-particle binning for the 3-particle");

  float maxBinVal;
  int maxBin;
  float bin, err;
  float phi, ss, mix, phiX, mixX, phiY, mixY;
  float ephi, ess, emix, ephiX, emixX, ephiY, emixY;
  //float b1=yPhiPhi->GetSum()/pow(yPhi->GetSum(),2);
  //float b2=yPhiPhiSS->GetSum()/pow(yPhiMix->GetSum(),2);
  // b=b1/b2;
  for(int i=0;i<nIter;i++){
    for(int q=0;q<nMin;q++){
      MinArray[q][0]=10000;
    }
    for(int x=1;x<=yPhiPhi->GetNbinsX();x++){
      phiX=0, ephiX=0, mixX=0, emixX=0;
      for(int j=(ReBin*(x-1)+1);j<ReBin*x;j++){
	phiX+=yPhi->GetBinContent(j);
	ephiX+=pow(yPhi->GetBinError(j),2);
	mixX+=yPhiMix->GetBinContent(j);
	emixX+=pow(yPhiMix->GetBinError(j),2);
      }
      for(int y=1;y<=yPhiPhi->GetNbinsY();y++){
	phiY=0, ephiY=0, mixY=0, emixY=0;
	for(int j=(ReBin*(y-1)+1);j<ReBin*y;j++){
	phiY+=yPhi->GetBinContent(j);
	ephiY+=pow(yPhi->GetBinError(j),2);
	mixY+=yPhiMix->GetBinContent(j);
	emixY+=pow(yPhiMix->GetBinError(j),2);
      }
	maxBinVal=-10000;
	for(int j=0;j<nMin;j++){//Find highest bin
	  if(MinArray[j][0]>maxBinVal){
	    maxBinVal=MinArray[j][0];
	    maxBin=j;
	  }
	}
	phi=yPhiPhi->GetBinContent(x,y);
	ephi=yPhiPhi->GetBinError(x,y);
	ss=yPhiPhiSS->GetBinContent(x,y);
	ess=yPhiPhiSS->GetBinError(x,y);
	mix=yPhiPhiMix->GetBinContent(x,y);
	emix=yPhiPhiMix->GetBinError(x,y);
	//xxx
	bin=phi-a[0]*mixY*(phiX-a[0]*mixX)-a[0]*mixX*(phiY-a[0]*mixY)-a[0]*a[0]*b*ss;
	//cout << bin << " " << maxBinVal << " " << maxBin << endl;
	if(bin<maxBinVal){//replce highest bin with current if lower
	  MinArray[maxBin][0]=bin;
	  //divide by 0 protections
	  if(a[0]==0)a[0]=1E-10;
	  if(a[1]==0)a[1]=1E-10;
	  if(ephi==0)ephi=1E-10;
	  if(ess==0)ess=1E-10;
	  if(emix==0)emix=1E-10;
	  if(phi==0)phi=1E-10;
	  if(ss==0)ss=1E-10;
	  if(mix==0)mix=1E-10;
	  MinArray[maxBin][1]=ephi*ephi+a[0]*mixY*phiX*(pow(emixY/mixY,2)+pow(ephiX/phiX,2))+2*a[0]*a[0]*mixY*mixX*(pow(emixY/mixY,2)+pow(emixX/mixX,2))+a[0]*mixX*phiY*(pow(emixX/mixX,2)+pow(ephiY/phiY,2))+pow(a[0]*a[0]*b*ess,2);
	  Phi[maxBin][0]=phi;
	  SS[maxBin][0]=ss;
	  Mix[maxBin][0]=mix;
	  PhiX[maxBin][0]=phiX;
	  PhiY[maxBin][0]=phiY;
	  MixX[maxBin][0]=MixX;
	  MixY[maxBin][0]=MixY;
	 
	  //cout << s << " " << b << " " << bin  << " " << maxBin << endl;
	}
      }
    }
    //sumSig=0, sumBg=0, esumSig=0, esumBg=0, sumMin=0, esumMin=0;
    float sumPhi=0,sumHS=0,sumSS=0,sumMix=0,sumMixXY=0, sumMin=0, esumMin=0;
    for(int j=0;j<nMin;j++){
      sumMin+=MinArray[j][0];
      esumMin+=MinArray[j][1];
      sumPhi+=Phi[j][0];
      sumHS+=(PhiX[j][0]*MixY[j][0]+PhiY[j][0]*MixX[j][0]);
      sumMix+=mix;
      sumSS+=b*ss;
      sumMixXY+=(2*MixY[j][0]*MixX[j][0]);
    }
    // if(sumSig==0){sumSig=1E-10; cout << "Warning: Zero sumSig" << endl;}
    // if(sumBg==0){sumBg=1E-10; cout << "Warning: Zero sumBg" << endl;}
    float a1,a2;
    // cout << pow(sumHS*sumHS-4*sumPhi*(sumMixXY-sumSS),0.5) << endl;
    a1=(sumHS+pow(sumHS*sumHS-4*sumPhi*(sumMixXY-sumSS),0.5))/(2*sumPhi*(sumMixXY-sumSS));
    a2=(sumHS-pow(sumHS*sumHS-4*sumPhi*(sumMixXY-sumSS),0.5))/(2*sumPhi*(sumMixXY-sumSS));
    if(fabs(oldmin)<fabs(sumMin))itersize/=2;
    oldmin=sumMin;
    if(a1>0&&a2>0){
      if(fabs(a1-1)<fabs(a2-1)){
	a[0]=a1;
      }
      else a[0]=a2;
    }
    else if(a1<0)a[0]=a2;
    else if(a2<0)a[0]=a1;
    else if(sumMin<0) a[0]-=itersize;
    else a[0]+=itersize;

   
  
    cout << "a1: " << a1 << "    a2: " << a2 << endl;
    //scalea[1]=scalea[0]*pow(esumSig/sumSig/sumSig+esumBg/sumBg/sumBg,0.5);
    esumMin=pow(esumMin,0.5);
    //esumSig=pow(esumSig,0.5);
    //esumBg=pow(esumBg,0.5);
    cout << "Iter:  " << i << "   Min=" << sumMin/nMin << " +/- " << esumMin/nMin << "  a=" << a[0] << endl;  
  }
}
  

//----------------------------------------------------------------------------------------------
void Add1D(TH1F *Histo, TH1F *HistAdd, float xscale){
  int nbinsx=Histo->GetNbinsX();
  for(int xi=1;xi<=nbinsx;xi++){
    if(fabs(Histo->GetBinContent(xi)>1E-10)&&fabs(HistAdd->GetBinContent(xi)>1E-10)){//only add if contents in both
      Histo->SetBinContent(xi,Histo->GetBinContent(xi)+xscale*HistAdd->GetBinContent(xi));
      Histo->SetBinError(xi,sqrt(pow(Histo->GetBinError(xi),2)+pow(xscale*HistAdd->GetBinError(xi),2)));
    }
  }
}
      
//----------------------------------------------------------------------------------
void Add2D(TH2F *Histo, TH2F *HistAdd,float xscale){
  int nbinsx=Histo->GetNbinsX();
  int nbinsy=Histo->GetNbinsY();
  //histOut=(TH2F*)Histo->Clone();
  // histOut->SetName("Added Histo");	
  
  for(int xi=1;xi<=nbinsx;xi++){
    for(int yi=1;yi<=nbinsy;yi++){
      //cout << Histo->GetBinContent(xi,yi) << endl;
      if(fabs(Histo->GetBinContent(xi,yi)>1E-10)&&fabs(HistAdd->GetBinContent(xi,yi))>1E-10){
      Histo->SetBinContent(xi,yi,Histo->GetBinContent(xi,yi)+xscale*HistAdd->GetBinContent(xi,yi));
      Histo->SetBinError(xi,yi,sqrt(pow(Histo->GetBinError(xi,yi),2)+pow(xscale*HistAdd->GetBinError(xi,yi),2)));
      }
    }
  }
}
//----------------------------------------------------------------------------------
void DivideMean1D(TH1F *Hist, TH1F *Hist2,float a1, float a2,float ntrig){
    int nbinsx=Hist->GetNbinsX();
    float bin,err,bin2,err2;
    float sigma=(a2-a1)/2;
    for(int xi=1;xi<=nbinsx;xi++){
      bin=Hist->GetBinContent(xi);
      err=Hist->GetBinError(xi);
      bin2=Hist2->GetBinContent(xi);
      err2=Hist2->GetBinError(xi);
      if(bin2==0){
	bin2=1;
	bin=0;
	err2=1;
      }
      Hist->SetBinContent(xi,bin/bin2);
      Hist->SetBinError(xi,sigma/(err2*ntrig*Hist->GetBinWidth(1)));
    }

}
//---------------------------------------------------------------------------------
//Correction factor screws with the statistics
void DivideMean1DCorr(TH1F *Hist, TH1F *Hist2, TH1F *Hist3, float a1, float a2,float ntrig){
    int nbinsx=Hist->GetNbinsX();
    float bin,err,bin2,err2, bin3, err3;
    float sigma=(a2-a1)/2;
    for(int xi=1;xi<=nbinsx;xi++){
      bin=Hist->GetBinContent(xi);
      err=Hist->GetBinError(xi);
      bin2=Hist2->GetBinContent(xi);
      err2=Hist2->GetBinError(xi);
      bin3=Hist3->GetBinContent(xi);
      err3=Hist3->GetBinError(xi);
      if(bin2==0){
	bin2=1;
	bin=0;
	err2=1;
	err3=1;
      }
      Hist->SetBinContent(xi,bin/bin2);
      Hist->SetBinError(xi,sigma/(err3*ntrig*Hist->GetBinWidth(1)));//error should go with uncorrected counts
    }
}

//-----------------------------
//Background subtraction screws with them more
void DivideMean1DBgSub(TH1F *Hist, TH1F *Hist2,TH1F *Hist3,TH1F *Hist4,float a1, float a2,float ntrig, float nmix){
    int nbinsx=Hist->GetNbinsX();
    float bin,err,bin2,err2,bin3,err3,bin4,err4;
    float sigma=(a2-a1);
    float e1,e2;
    for(int xi=1;xi<=nbinsx;xi++){
      bin=Hist->GetBinContent(xi);
      err=Hist->GetBinError(xi);
      bin2=Hist2->GetBinContent(xi);
      err2=Hist2->GetBinError(xi);
      bin3=Hist3->GetBinContent(xi);
      err3=Hist3->GetBinError(xi);
      bin4=Hist4->GetBinContent(xi);
      err4=Hist4->GetBinError(xi);
      if(bin2==0){
	bin2=1;
	bin=0;
	err2=1;
	err3=1;
	err4=1;
      }
      Hist->SetBinContent(xi,bin/bin2);
      cout << sigma << endl;
      e1=pow(sigma/(err3*ntrig*Hist->GetBinWidth(xi)),2);
      e2=pow(sigma/(err4*nmix*Hist->GetBinWidth(xi)),2);
      cout << e1 << endl;
      
      Hist->SetBinError(xi,sqrt(e1));
    }
}
//----------------------------------------------------------------------------------
void DivideMean2D(TH2F *Hista, TH2F *Hista2, float a1, float a2,float ntrig){
  int nbinsx=Hista->GetXaxis()->GetNbins();
  int nbinsy=Hista->GetXaxis()->GetNbins();
    float bin,err,bin2,err2;
    float sigma=(a2-a1)/2;
    for(int xi=1;xi<=nbinsx;xi++){
      for(int yi=1;yi<=nbinsy;yi++){
	bin=Hista->GetBinContent(xi,yi);
	err=Hista->GetBinError(xi,yi);
	bin2=Hista2->GetBinContent(xi,yi);
	err2=Hista2->GetBinError(xi,yi);
	if(bin2==0){
	  bin2=1;
	  bin=0;
	  err2=1;
	}

	Hista->SetBinContent(xi,yi,bin/bin2);
	Hista->SetBinError(xi,yi,sigma/(err2*ntrig*HistGetXaxis()->GetBinWidth(1)*HistGetYaxis()->GetBinWidth(1)));
      }
    }
}
//----------------------------------------------------------------------------------------------------------------
void SetMargins1D(TCanvas *xc){
  xc->SetRightMargin(0.02);
  xc->SetLeftMargin(0.16);
  xc->SetBottomMargin(0.15);
}

void SetMargins1DCyl(TCanvas *xc){
  xc->SetRightMargin(0.1);
  xc->SetLeftMargin(0.1);
  xc->SetTopMargin(0.05);
  xc->SetBottomMargin(0.05);
  xc->SetTheta(90);
  xc->SetPhi(209);//151,209
}

void SetMargins2DSurf(TCanvas *xc){
  xc->SetBottomMargin(0.15);
  xc->SetLeftMargin(0.21);
  xc->SetRightMargin(0.05);
  xc->SetTopMargin(0.05);
}
void SetMargins2D(TCanvas *xc){
  xc->SetBottomMargin(0.15);
  xc->SetLeftMargin(0.12);
  xc->SetRightMargin(0.25);
}
void SetTitles1D(TH1 *xhisto){
  xhisto->GetXaxis()->SetTitleColor(1);
  xhisto->GetXaxis()->SetTitleSize(0.06);
  xhisto->GetYaxis()->SetTitleSize(0.06);
  xhisto->GetXaxis()->SetLabelSize(0.05);
  xhisto->GetYaxis()->SetLabelSize(0.05);
  xhisto->GetYaxis()->SetTitleOffset(1.2);
}
void SetTitles2D(TH1 *xhisto){
  xhisto->GetXaxis()->SetTitleColor(1);
  xhisto->GetXaxis()->SetTitleSize(0.06);
  xhisto->GetYaxis()->SetTitleSize(0.06);
  xhisto->GetZaxis()->SetTitleSize(0.06);
  xhisto->GetXaxis()->SetLabelSize(0.05);
  xhisto->GetYaxis()->SetLabelSize(0.05);
  xhisto->GetZaxis()->SetLabelSize(0.05);
  xhisto->GetZaxis()->SetTitleOffset(1.43);
}

void Style(int i) 
{
  switch (i) {
  case 1:
    gStyle->SetStatColor(kWhite);
    gStyle->SetTitleColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(1);
    break;
  case 2:  //paper figure style
    gStyle->SetStatColor(kWhite);
    gStyle->SetTitleColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetPalette(1);
    break;
  default:
    break;
  }
}
//------------------------------------------------------------------
void keySymbol(Float_t x, Float_t y, const Char_t* tt, 
               Int_t color=1, Int_t marker=1, Float_t tSize=0.04, Float_t Msize=1.0){

  gPad->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(12);
  t->SetTextSize(tSize);
  t->SetTextColor(color);
  t->DrawLatex(x+0.025,y,tt);
  
      
  TMarker *l = new TMarker();
  l->SetNDC(kTRUE);
  l->SetMarkerStyle(marker);
  l->SetMarkerColor(color);
  l->SetMarkerSize(Msize);
  l->PaintMarkerNDC(x-0.0025,y);
  l->SetX(x-0.0025);
  l->SetY(y);
  l->Draw();
}
//------------------------------------------------
void keyText(Float_t x, Float_t y, const Char_t* tt, 
	     Int_t color=1,  Float_t tSize=0.04, Float_t tangle=0){
 
  gPad->Update();
  Float_t x1 = gPad->GetFrame()->GetX1();
  Float_t x2 = gPad->GetFrame()->GetX2();
  Float_t dx = x2-x1;
  x = x1 + dx*x;
  Float_t y1 = gPad->GetFrame()->GetY1();
  Float_t y2 = gPad->GetFrame()->GetY2();
  Float_t dy = y2-y1;
  y = y1 + dy*y; 
TLatex *t = new TLatex(x,y,tt);
  t->SetTextAlign(12);
  t->SetTextSize(tSize);
  t->SetTextColor(color);
  t->SetTextAngle(tangle);
  t->Draw();
}
////////////////////////////////////////////////
void rebin1D(TH1 *&temp1D, int rbins){
  int nbinxx1=temp1D->GetNbinsX();
  int nbinxx2=nbinxx1/rbins;
  float widxx=temp1D->GetBinWidth(1);
  float minxx=temp1D->GetBinCenter(1)-widxx/2;
  float maxxx=temp1D->GetBinCenter(nbinxx1)+widxx/2;
  char *titxx=temp1D->GetTitle();
 
  TH1D *temp1D2=new TH1D("temp1D2",titxx,nbinxx2,minxx,maxxx);
  for(int xx1=1;xx1<=nbinxx2;xx1++){
    float conx1=0,errx1=0;
    for(int xx2=1;xx2<=rbins;xx2++){
      conx1+=temp1D->GetBinContent(rbins*(xx1-1)+xx2);
   
      errx1+=pow(temp1D->GetBinError(rbins*(xx1-1)+xx2),2);

    }

    conx1/=rbins;
    errx1=sqrt(errx1)/rbins;
    temp1D2->SetBinContent(xx1,conx1);
    temp1D2->SetBinError(xx1,errx1);
  }
  char *namex=temp1D->GetName();
  char *xtit=temp1D->GetXaxis()->GetTitle();
  char *ytit=temp1D->GetYaxis()->GetTitle();
  int markx=temp1D->GetMarkerStyle();
  int colx=temp1D->GetMarkerColor();
  temp1D=(TH1D*)temp1D2->Clone();
  temp1D->SetName(namex);
  temp1D->GetXaxis()->SetTitle(xtit);
  temp1D->GetYaxis()->SetTitle(ytit);
  temp1D->SetMarkerStyle(markx);
  temp1D->SetLineColor(colx);
  temp1D->SetMarkerColor(colx);
  SetTitles1D(temp1D);
 
}
////////////////////////////////
//Broken need to fix before I use
void rebin2D(TH2D *temp2D,int rbins){
  int nbinxx1=temp2D->GetNbinsX();
  int nbinxx2=nbinxx2/rbins;
  float widxx=temp2D->GetBinWidth(1);
  float minxx=temp2D->GetBinCenter(1)-widxx/2;
  float maxxx=temp2D->GetBinCenter(nbinxx1)+widxx/2;
  char titxx=temp2D->GetTitle();
  TH2D *temp2D2=new TH2D("temp2D2",titxx,nbinxx2,minxx,maxxx,nbinxx2,minxx,maxxx);
  for(int xx1=1;xx1<=nbinxx2;xx1++){
    for(int yy1=1;yy1<=nbinxx2;yy1++){
      float conx1=0,errx1=0;
      for(int xx2=1;xx2<=rbins;xx2++){
	for(int yy2=1;yy2<=rbins;yy2++){
	  conx1+=temp2D->GetBinContent(rbins*(xx1-1)+xx2,rebins*(yy1-1)+yy2);
	  errx1+=pow(temp2D->GetBinError(rbins(xx1-1)+xx2,rebins*(yy1-1)+yy2),2);
	}
      }
      conx1/=pow(rbins,2);
      errx1=sqrt(errx1)/pow(rbins,2);
      temp2D2->SetBinContent(xx1,conx1);
      temp2D2->SetBinError(xx1,errx1);
    }
  }
  return temp2D2;
}
//----------------------------------
void DrawLego1D(TCanvas *&tc, TH1 *&thist){

  //cout << "off" << thist->GetXaxis()->GetLabelOffset() << endl;
  if(thist->GetMinimum()>0)thist->SetMinimum(0);
  thist->GetXaxis()->SetLabelOffset(-0.1);
  tc1=new TPad("tc1","ytit",0,0.0,0.1,1);
  tc2=new TPad("tc2","main",0.1,0,1,1);
  tc1->Draw();
  tc2->Draw();
  tc2->cd();
  tc2->SetTheta(0);
  tc2->SetPhi(0);
  tc2->SetRightMargin(0.02);
  //xc->SetLeftMargin(0.16);
  tc2->SetBottomMargin(0.15);
  thist->Draw("lego2");
  tc1->cd();
  keyText(.45,.65,"#frac{1}{N_{Trig}}#frac{dN}{d#Delta#phi}",1,0.4,90);
  // thist->GetXaxis()->SetLabelOffset(0);
}
//---------------------------------------
void SaveAsText(TH1 *histo, char *filename){
  
  //I should also make it works with TH2
  FILE *fp;
  fp=fopen(filename,"w");
  fprintf(fp,"Bin\t%s\tWidth\t%s\tError\n",histo->GetXaxis()->GetTitle(),histo->GetYaxis()->GetTitle());
  for(int x=1;x<=histo->GetNbinsX();x++){
    fprintf(fp,"%d\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n",x,histo->GetBinCenter(x),histo->GetBinWidth(x),histo->GetBinContent(x),histo->GetBinError(x));
  }
  fclose(fp);
}
  
//==========================================
void MixedCorrect(TH1F *hPhi, TH1F *hPhiMix, TH1F *hEtaN, TH1F *hEtaNMix, TH1F *hEtaA, TH1F *hEtaAMix, TH2F *hPhiEta, TH2F *hPhiEtaMix){

//Do not use GeoCorr if using this
  // float avePhiMix=hPhiMix->GetSum()/hPhiMix->GetNbinsX();
  float avePhiMix=0, eavePhiMix=0;
  for(int i=0;i<=hPhiMix->GetNbinsX();i++){
    avePhiMix+=hPhiMix->GetBinContent(i);
    eavePhiMix+=pow(hPhiMix->GetBinError(i),2);
  }
  avePhiMix/=hPhiMix->GetNbinsX();
  eavePhiMix=pow(eavePhiMix,0.5)/hPhiMix->GetNbinsX();
  hPhiMix->Scale(1/avePhiMix);
  hPhi->Divide(hPhiMix);

  
  float maxEtaNMix=hEtaNMix->GetBinContent(hEtaNMix->GetNbinsX()/2);
  // float maxEtaNMix=0, emaxEtaNMix=0;
  // for(int i=0;i<=hEtaNMix->GetNbinsX();i++){
  //   aveEtaMix+=hEtaMix->GetBinContent(i);
  //   eaveEtaMix+=pow(hEtaMix->GetBinError(i),2);
  //  }
  //aveEtaMix/=hEtaMix->GetNbinsX();
  // eaveEtaMix=pow(eaveEtaMix,0.5)/hEtaMix->GetNbinsX();
  hEtaNMix->Scale(1/maxEtaNMix);
  hEtaN->Divide(hEtaNMix);
  
  float maxEtaAMix=hEtaAMix->GetBinContent(hEtaAMix->GetNbinsX()/2);
  hEtaAMix->Scale(1/maxEtaAMix);
  hEtaA->Divide(hEtaAMix);

  float maxavePhiEta=0,emaxavePhiEta=0;
  float nbin=0;
  for(int i=1;i<=hPhiEtaMix->GetNbinsX();i++){
    maxavePhiEta+=hPhiEtaMix->GetBinContent(i,hPhiEtaMix->GetNbinsY()/2);
    emaxavePhiEta+=pow(hPhiEtaMix->GetBinError(i,hPhiEtaMix->GetNbinsY()/2),2);
  }
  maxavePhiEta/=hPhiEtaMix->GetNbinsX();
  emaxavePhiEta=pow(emaxavePhiEta,0.5)/hPhiEtaMix->GetNbinsX();
  hPhiEtaMix->Scale(1/maxavePhiEta);
  hPhiEta->Divide(hPhiEtaMix);

  for(int i=1;i<=hPhiMix->GetNbinsX();i++){
    hPhiMix->SetBinContent(i,avePhiMix);
    hPhiMix->SetBinError(i,0);
    }
 
  for(int i=1;i<=hEtaNMix->GetNbinsX();i++){
    hEtaNMix->SetBinContent(i,maxEtaNMix);
    hEtaAMix->SetBinContent(i,maxEtaAMix);
    hEtaNMix->SetBinError(i,0);
    hEtaAMix->SetBinError(i,0);
  }

  for(int i=1;i<=hPhiEta->GetNbinsX();i++){
    for(int j=1;j<=hPhiEta->GetNbinsY();j++){
      hPhiEtaMix->SetBinContent(i,j,maxavePhiEta);
      hPhiEtaMix->SetBinError(i,j,0);
    }
  }
}

void MixedCorrect3(TH2F *hPhiPhi, TH2F *hPhiPhiSS, TH2F *hPhiPhiMix,TH2F *hEtaEta, TH2F *hEtaEtaSS, TH2F *hEtaEtaMix){
  float avePhiPhiMix=hPhiPhiMix->GetSum()/pow(hPhiPhiMix->GetNbinsX(),2);
  cout << "hPhiPhiMix: " << hPhiPhiMix->GetBinContent(1,1) << " " << hPhiPhiMix->GetBinError(1,1) << endl;
  hPhiPhiMix->Scale(1/avePhiPhiMix);
  cout << "hPhiPhiMix: " << hPhiPhiMix->GetBinContent(1,1) << " " << hPhiPhiMix->GetBinError(1,1) << endl;
  hPhiPhi->Divide(hPhiPhiMix);
  hPhiPhiSS->Divide(hPhiPhiMix);

  float maxEtaEtaMix=hEtaEtaMix->GetBinContent(hEtaEtaMix->GetNbinsX()/2,hEtaEtaMix->GetNbinsX()/2);
  hEtaEtaMix->Scale(1/maxEtaEtaMix);
  hEtaEta->Divide(hEtaEtaMix);
  hEtaEtaSS->Divide(hEtaEtaMix);
  
  float conX, conY;
  int nbin=hEtaEta->GetNbinsX();
  float EtaCut=hEtaEta->GetBinCenter(nbin)+hEtaEta->GetBinWidth(1)/2.;

  for(int i=1;i<=hPhiPhi->GetNbinsX();i++){
 
    for(int j=0;j<=hPhiPhi->GetNbinsX();j++){  
      hPhiPhiMix->SetBinContent(i,j,avePhiPhiMix);
      hPhiPhiMix->SetBinError(i,j,0);
    }
  }
  for(int i=1;i<=hEtaEta->GetNbinsX();i++){
    conX=hEtaEta->GetBinCenter(i);
    for(int j=0;j<=hEtaEta->GetNbinsX();j++){
      conY=hEtaEta->GetBinCenter(j);
      if(hEtaEtaMix->GetBinContent(i,j)){
	hEtaEtaMix->SetBinContent(i,j,maxEtaEtaMix);
	hEtaEtaMix->SetBinError(i,j,0);
      }
    }
  }
}
//=================================================================

void ProjPhiPhi(TH2F *shist, TH1F *Non1D, TH1F *Noff1D, TH1F *Aon1D, TH1F *Aoff1D , Float_t NearCut=1.5, float DiagWidth=0.35,int nrand=1000){
  delete gRandom;
  float pi=3.14159265;
  int binxy=shist->GetNbinsX();
  float widthxy=shist->GetXaxis()->GetBinWidth(1);
  gRandom=new TRandom3();
  gRandom->SetSeed(1);
  float s,d;
  int Nbinxy=2.8*NearCut/widthxy;
  float AwayCut=pi-NearCut;
  int Abinxy=2.8*AwayCut/widthxy;
  *Non1D=new TH1F("Non1D","Near On-Diag",Nbinxy,-NearCut,NearCut);
  *Noff1D=new TH1F("Noff1D","Near Off-Diag",Nbinxy,-NearCut,NearCut);
  *Aon1D=new TH1F("Aon1D","Away On-Diag",Abinxy,-AwayCut,AwayCut);
  *Aoff1D=new TH1F("Aoff1D","Away Off-Diag",Abinxy,-AwayCut,AwayCut);

  Non1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2  or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Noff1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2 or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Aon1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2-#pi or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Aoff1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2-#pi or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Non1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  Noff1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  Aon1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  Aoff1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  SetTitles1D(Non1D);
  SetTitles1D(Noff1D);
  SetTitles1D(Aon1D);
  SetTitles1D(Aoff1D);

  TH1F *eNon1D=(TH1F*)Non1D->Clone();
  eNon1D->SetName("eNon1D");
  TH1F *eNoff1D=(TH1F*)Noff1D->Clone();
  eNoff1D->SetName("eNoff1D");
  TH1F *eAon1D=(TH1F*)Aon1D->Clone();
  eAon1D->SetName("eAon1D");
  TH1F *eAoff1D=(TH1F*)Aoff1D->Clone();
  eAoff1D->SetName("eAoff1D");
  
  //float binsw=2*pi/binxy*2*pi/(2*pi-2);
  float binsw;
  float NbinsW=Non1D->GetBinWidth(1);
  float AbinsW=Aon1D->GetBinWidth(1);
  float x10,y10,x20,y20,conp1,errp1;
  int nears=10;
  for(int xx1=1;xx1<=binxy;xx1++){
    x10=shist->GetBinCenter(xx1);
    for(int yy1=1;yy1<=binxy;yy1++){
      y10=shist->GetBinCenter(yy1);
      if((x10<NearCut&&y10<NearCut)){binsw=NbinsW/widthxy; nears=1;}
      else if((x10>NearCut&&y10>NearCut)){binsw=AbinsW/widthxy;nears=0;}
      else continue;
      conp1=shist->GetBinContent(xx1,yy1)*binsw/nrand;
      errp1=pow(shist->GetBinError(xx1,yy1)*binsw,2)/nrand;
      //     cout << conp1 << " " << errp1 << " " << endl; 
      // cout << shist->GetBinContent(xx1,yy1) << " " << shist->GetBinError(xx1,yy1) << endl; 
      for(int r1=0;r1<nrand;r1++){
	x20=x10+binsw*(gRandom->Rndm()-0.5);
	y20=y10+binsw*(gRandom->Rndm()-0.5);
	s=(x20+y20)/2;
	d=(x20-y20)/2;
	if(fabs(s-pi)<DiagWidth&&nears==0){
	  Aoff1D->Fill(d,conp1);
	  eAoff1D->Fill(d,errp1);
	  
	}
	
	if(fabs(s)<DiagWidth&&nears==1){
	  Noff1D->Fill(d,conp1);
	  eNoff1D->Fill(d,errp1);
	}
	if(d<DiagWidth&&d>0){
	  if(nears==0){
	    Aon1D->Fill((s-pi),2*conp1);//conp1*2
	    eAon1D->Fill((s-pi),4*errp1);//errp1*4
	  }
	  else if(nears==1){
	    Non1D->Fill(s,2*conp1);//conp1*2
	    eNon1D->Fill(s,4*errp1);//errp1*4
	  }
	}
      }
    }
  }
  for(int xx1=1;xx1<=Nbinxy;xx1++){
    Noff1D->SetBinError(xx1,sqrt(eNoff1D->GetBinContent(xx1)));
    Non1D->SetBinError(xx1,sqrt(eNon1D->GetBinContent(xx1)));
  }
  for(int xx1=1;xx1<=Abinxy;xx1++){
    Aoff1D->SetBinError(xx1,sqrt(eAoff1D->GetBinContent(xx1)));
    Aon1D->SetBinError(xx1,sqrt(eAon1D->GetBinContent(xx1)));
  }
  Non1D->SetMarkerStyle(25);
  Non1D->SetMarkerColor(4);
  Non1D->SetLineColor(4);
  Noff1D->SetMarkerStyle(20);
  Noff1D->SetMarkerColor(2);
  Noff1D->SetLineColor(2);
 Aon1D->SetMarkerStyle(25);
 Aon1D->SetMarkerColor(4);
 Aon1D->SetLineColor(4);
 Aoff1D->SetMarkerStyle(20);
  Aoff1D->SetMarkerColor(2);
  Aoff1D->SetLineColor(2);
}
//========================================
ProjEtaEta(TH2F *shist, TH1F *Non1D, TH1F *Noff1D, float DiagWidth=0.35,int nrand=1000){
  float pi=3.14159265;
  int binxy=shist->GetNbinsX();
  float widthxy=shist->GetXaxis()->GetBinWidth(1);
  gRandom=new TRandom3();
  gRandom->SetSeed(1);
  float NearCut=shist->GetBinCenter(binxy)+shist->GetBinWidth(binxy)/2.;
  float s,d;
  int Nbinxy=2.8*NearCut/widthxy;

  *Non1D=new TH1F("NEon1D","Near On-Diag",Nbinxy,-NearCut,NearCut);
  *Noff1D=new TH1F("NEoff1D","Near Off-Diag",Nbinxy,-NearCut,NearCut);

  Non1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2  or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Noff1D->GetXaxis()->SetTitle("#Sigma=(#Delta#phi_{1}+#Delta#phi_{2})/2 or #Delta=(#Delta#phi_{1}-#Delta#phi_{2})/2     ");
  Non1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  Noff1D->GetYaxis()->SetTitle("1/N_{Trig} dN/d(#Sigma or #Delta)   ");
  SetTitles1D(Non1D);
  SetTitles1D(Noff1D);

 TH1F *eNon1D=(TH1F*)Non1D->Clone();
  eNon1D->SetName("eNEon1D");
  TH1F *eNoff1D=(TH1F*)Noff1D->Clone();
  eNoff1D->SetName("eNEoff1D");

  float binsw=Non1D->GetBinWidth(1);
  float x10,y10,x20,y20,conp1,errp1;

  for(int xx1=1;xx1<=binxy;xx1++){
    x10=shist->GetBinCenter(xx1);
    for(int yy1=1;yy1<=binxy;yy1++){
      y10=shist->GetBinCenter(yy1);
      conp1=shist->GetBinContent(xx1,yy1)*binsw/nrand;
      errp1=pow(shist->GetBinError(xx1,yy1)*binsw,2)/nrand;
      for(int r1=0;r1<nrand;r1++){
	x20=x10+binsw*(gRandom->Rndm()-0.5);
	y20=y10+binsw*(gRandom->Rndm()-0.5);
	s=(x20+y20)/2;
	d=(x20-y20)/2;

	if(fabs(s)<DiagWidth){
	  Noff1D->Fill(d,conp1);
	  eNoff1D->Fill(d,errp1);
	}
	if(d<DiagWidth&&d>0){
	  Non1D->Fill(s,2*conp1);//conp1*2
	  eNon1D->Fill(s,4*errp1);//errp1*4
	}
      }
    }
  }
 for(int xx1=1;xx1<=Nbinxy;xx1++){
    Noff1D->SetBinError(xx1,sqrt(eNoff1D->GetBinContent(xx1)));
    Non1D->SetBinError(xx1,sqrt(eNon1D->GetBinContent(xx1)));
  }
Non1D->SetMarkerStyle(25);
  Non1D->SetMarkerColor(4);
  Non1D->SetLineColor(4);
  Noff1D->SetMarkerStyle(20);
  Noff1D->SetMarkerColor(2);
  Noff1D->SetLineColor(2);
}
