#include <iostream>
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TNamed.h"
#include "TList.h"
#include "TDirectory.h"
#include <X3DDefs.h>
#include "TPaveText.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TObjString.h"
#include "THn.h"
#include "TArrayD.h"
#include <TMatrixTBase.h>
#include "TParameter.h"


using namespace std;

THnF * SumCent(Double_t Mincent , Double_t Maxcent, THnF* source, const char * name){
  THnF * output = (THnF*)source->Clone(name);
  output->Reset();
  for(int mult = 0; mult<=source->GetAxis(0)->GetNbins();mult++){
    double lowedge = source->GetAxis(0)->GetBinLowEdge(mult);
    double upedge = source->GetAxis(0)->GetBinUpEdge(mult);
    double multval = source->GetAxis(0)->GetBinCenter(mult);
    if(lowedge>=Mincent){
      if(upedge<=Maxcent){
      for(int vz = 0; vz<=source->GetAxis(1)->GetNbins();vz++){
	double vzval = source->GetAxis(1)->GetBinCenter(vz);
	for(int eta = 0;eta<=source->GetAxis(3)->GetNbins();eta++){
	  //find the corresponding bin:
	  double etaval = source->GetAxis(3)->GetBinCenter(eta);
	  for(int pT = 0;pT<=source->GetAxis(4)->GetNbins();pT++){
	    //find the corresponding bin:
	    double pTval = source->GetAxis(4)->GetBinCenter(pT);
	    for(int phi = 0;phi<=source->GetAxis(2)->GetNbins();phi++){
	    //find the corresponding bin:
	    double phival = source->GetAxis(2)->GetBinCenter(phi);
	    int bin[5] = {mult,vz,phi,eta,pT};
	      Double_t value[5] = {multval,vzval,phival,etaval,pTval};
	      for(int i=1;i<=source->GetBinContent(bin);i++){
		output->Fill(value);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return output;
}

TH2D * SumCent(Double_t Mincent , Double_t Maxcent, TH2D* source, const char * name){
  TH2D * output = (TH2D*)source->Clone(name);
  output->Reset();
  for(int mult = 0; mult<=source->GetXaxis()->GetNbins();mult++){
    double lowedge = source->GetXaxis()->GetBinLowEdge(mult);
    double upedge = source->GetXaxis()->GetBinUpEdge(mult);
    double multval = source->GetXaxis()->GetBinCenter(mult);
    if(lowedge>=Mincent){
      if(upedge<=Maxcent){
      for(int vz = 0; vz<=source->GetYaxis()->GetNbins();vz++){
	double vzval = source->GetYaxis()->GetBinCenter(vz);
	for(int i=1;i<=source->GetBinContent(mult,vz);i++){
	  output->Fill(multval,vzval);
	  }      
	}	
      }
    }
  }
  return output;
}

THnF * Rebin(Int_t Nmult,const Double_t * MultAxis,Int_t nvz, const Double_t * VZAxis,Int_t nphi, Double_t phimin, Double_t phimax,Int_t neta, Double_t etamin,Double_t etamax,Int_t npt, const Double_t * pTAxis,THnF * input, const char* name){
  //Function that takes a THnF, and rebins it to the given axes.
  THnF * Resulthist;
  if(nphi>1){
    Int_t    nbins[5] 	= {Nmult,nvz,nphi,neta,npt};
    Double_t xmin[5] 	= {MultAxis[0],VZAxis[0],phimin,etamin,pTAxis[0]};
    Double_t xmax[5] 	= {MultAxis[Nmult],VZAxis[nvz],phimax,etamax,pTAxis[npt]};
    Resulthist = new THnF(name,input->GetTitle(), 5,nbins,xmin,xmax);
    Resulthist->GetAxis(2)->SetTitle(input->GetAxis(2)->GetTitle());
    Resulthist->GetAxis(3)->SetTitle(input->GetAxis(3)->GetTitle());
    Resulthist->GetAxis(4)->SetTitle(input->GetAxis(4)->GetTitle());
    Resulthist->GetAxis(4)->Set(npt,pTAxis);
  }
  else{
    Int_t    nbins[4] 	= {Nmult,nvz,neta,npt};
    Double_t xmin[4] 	= {MultAxis[0],VZAxis[0],etamin,pTAxis[0]};
    Double_t xmax[4] 	= {MultAxis[Nmult],VZAxis[nvz],etamax,pTAxis[npt]};
    Resulthist = new THnF(name,input->GetTitle(), 4,nbins,xmin,xmax);
    Resulthist->GetAxis(2)->SetTitle(input->GetAxis(3)->GetTitle());
    Resulthist->GetAxis(3)->SetTitle(input->GetAxis(4)->GetTitle());
    Resulthist->GetAxis(3)->Set(npt,pTAxis);    
  }
  Resulthist->GetAxis(0)->SetTitle(input->GetAxis(0)->GetTitle());
  Resulthist->GetAxis(0)->Set(Nmult,MultAxis);
  Resulthist->GetAxis(1)->SetTitle(input->GetAxis(1)->GetTitle());
  Resulthist->GetAxis(1)->Set(nvz,VZAxis);

  Resulthist->Sumw2();

  for(int mult = 0;mult<=input->GetAxis(0)->GetNbins();mult++){
    //find the corresponding bin:
    double multval = input->GetAxis(0)->GetBinCenter(mult);
    for(int vz = 0;vz<=input->GetAxis(1)->GetNbins();vz++){
      //find the corresponding bin:
      double vzval = input->GetAxis(1)->GetBinCenter(vz);
	for(int eta = 0;eta<=input->GetAxis(3)->GetNbins();eta++){
	  //find the corresponding bin:
	  double etaval = input->GetAxis(3)->GetBinCenter(eta);
	  for(int pT = 0;pT<=input->GetAxis(4)->GetNbins();pT++){
	    //find the corresponding bin:
	    double pTval = input->GetAxis(4)->GetBinCenter(pT);
	    for(int phi = 0;phi<=input->GetAxis(2)->GetNbins();phi++){
	    //find the corresponding bin:
	    double phival = input->GetAxis(2)->GetBinCenter(phi);
	    int bin[5] = {mult,vz,phi,eta,pT};
	    if(nphi>1){
	      Double_t value[5] = {multval,vzval,phival,etaval,pTval};
	      for(int i=1;i<=input->GetBinContent(bin);i++){
		Resulthist->Fill(value);
	      }
	    }
	    else{
	      Double_t value[4] = {multval,vzval,etaval,pTval};
	      for(int i=1;i<=input->GetBinContent(bin);i++){
		Resulthist->Fill(value);
	      }
	    }
	  }	      
	}
      }
    }
  }

  return Resulthist;
}

THnF * Rebin(Int_t Nmult,const Double_t * MultAxis,Int_t nvz, const Double_t * VZAxis,Int_t nphi, Double_t * phiaxis,Int_t neta, Double_t * etaaxis,Int_t npt, const Double_t * pTAxis,THnF * input, const char* name){
  //Function that takes a THnF, and rebins it to the given axes.
  Int_t    nbins[5] 	= {Nmult,nvz,nphi,neta,npt};
  Double_t xmin[5] 	= {MultAxis[0],VZAxis[0],phiaxis[0],etaaxis[0],pTAxis[0]};
  Double_t xmax[5] 	= {MultAxis[Nmult],VZAxis[nvz],phiaxis[nphi],etaaxis[neta],pTAxis[npt]};
  
  THnF * Resulthist = new THnF(name,input->GetTitle(), 5,nbins,xmin,xmax);
  Resulthist->GetAxis(0)->SetTitle(input->GetAxis(0)->GetTitle());
  Resulthist->GetAxis(0)->Set(Nmult,MultAxis);
  Resulthist->GetAxis(1)->SetTitle(input->GetAxis(1)->GetTitle());
  Resulthist->GetAxis(1)->Set(nvz,VZAxis);
  Resulthist->GetAxis(2)->SetTitle(input->GetAxis(2)->GetTitle());
  Resulthist->GetAxis(2)->Set(nphi,phiaxis);
  Resulthist->GetAxis(3)->SetTitle(input->GetAxis(3)->GetTitle());
  Resulthist->GetAxis(3)->Set(neta,etaaxis);
  Resulthist->GetAxis(4)->SetTitle(input->GetAxis(4)->GetTitle());
  Resulthist->GetAxis(4)->Set(npt,pTAxis);
  Resulthist->Sumw2();

  for(int mult = 0;mult<=input->GetAxis(0)->GetNbins();mult++){
    //find the corresponding bin:
    double multval = input->GetAxis(0)->GetBinCenter(mult);
    for(int vz = 0;vz<=input->GetAxis(1)->GetNbins();vz++){
      //find the corresponding bin:
      double vzval = input->GetAxis(1)->GetBinCenter(vz);
      for(int phi = 0;phi<=input->GetAxis(2)->GetNbins();phi++){
	//find the corresponding bin:
	double phival = input->GetAxis(2)->GetBinCenter(phi);
	for(int eta = 0;eta<=input->GetAxis(3)->GetNbins();eta++){
	  //find the corresponding bin:
	  double etaval = input->GetAxis(3)->GetBinCenter(eta);
	  for(int pT = 0;pT<=input->GetAxis(4)->GetNbins();pT++){
	    //find the corresponding bin:
	    double pTval = input->GetAxis(4)->GetBinCenter(pT);
	    int bin[5] = {mult,vz,phi,eta,pT};
	    Double_t value[5] = {multval,vzval,phival,etaval,pTval};
	    for(int i=1;i<=input->GetBinContent(bin);i++){
	      Resulthist->Fill(value);
	    }
	  }	
	}
      }
    }
  }

  return Resulthist;
}

THnF* Rebin(THnF * input, const char* name,TArrayD * multax,TArrayD* vzaxis, TArrayD* phiaxis, TArrayD* etaaxis, TArrayD* ptaxis){
  return Rebin(multax->GetSize()-1,multax->GetArray(),vzaxis->GetSize()-1,vzaxis->GetArray(),phiaxis->GetSize()-1,phiaxis->GetArray(),etaaxis->GetSize()-1,etaaxis->GetArray(),ptaxis->GetSize()-1,ptaxis->GetArray(),input,name);
}

Float_t GetBinAx(THnF * hist,Int_t Ax1,Int_t Bin1,Int_t Ax2,Int_t Bin2,Int_t Ax3,Int_t Bin3,Int_t Ax4,Int_t Bin4,Int_t Ax5,Int_t Bin5){
  int Bin[5] = {0,0,0,0,0};
  Bin[Ax1] = Bin1;
  Bin[Ax2] = Bin2;
  Bin[Ax3] = Bin3;
  Bin[Ax4] = Bin4;
  Bin[Ax5] = Bin5;
  return hist->GetBinContent(Bin);
}

TH3D * Rebin(THnF * input, const char* name, Int_t Ax1, TArrayD * axis1,Int_t Ax2, TArrayD * axis2, Int_t Ax3, TArrayD* axis3, Int_t Ax4, Double_t Ax4min,Int_t Ax5){
  TH3D * hist3d = new TH3D(name,"",axis1->GetSize()-1,axis1->GetArray(),axis2->GetSize()-1,axis2->GetArray(),axis3->GetSize()-1,axis3->GetArray());
  TAxis* Axis1 = input->GetAxis(Ax1);
  hist3d->GetXaxis()->SetTitle(Axis1->GetTitle());
  TAxis* Axis2 = input->GetAxis(Ax2);
  hist3d->GetYaxis()->SetTitle(Axis2->GetTitle());
  TAxis* Axis3 = input->GetAxis(Ax3);
  hist3d->GetZaxis()->SetTitle(Axis3->GetTitle());
  TAxis* Axis4 = input->GetAxis(Ax4);
//   cout << Axis4->GetBinCenter(Axis4->GetNbins())<<endl;
  TAxis* Axis5 = input->GetAxis(Ax5);
  hist3d->Sumw2();
  for(int x = 0;x<=Axis1->GetNbins()+1;x++){
    Double_t Ax1val = Axis1->GetBinCenter(x);
    for(int y = 0;y<=Axis2->GetNbins()+1;y++){
      Double_t Ax2val = Axis2->GetBinCenter(y);
      for(int z = 0;z<=Axis3->GetNbins()+1;z++){
	Double_t Ax3val = Axis3->GetBinCenter(z);
	for(int x1=1;x1<=Axis4->GetNbins();x1++){
	  Double_t Ax4val = Axis4->GetBinLowEdge(x1);
	  if(Ax4val<Ax4min){ continue;}
	  for(int x2=1;x2<=Axis5->GetNbins();x2++){
	    Double_t content = GetBinAx(input,Ax1,x,Ax2,y,Ax3,z,Ax4,x1,Ax5,x2);
	    for(int i=1;i<=content;i++){
	      hist3d->Fill(Ax1val,Ax2val,Ax3val);
	    }
	  }
	}
      }
    }
  }
  return hist3d;
}

TH3D * Rebin(THnF * input, const char* name, Int_t Ax1, TArrayD * axis1,Int_t Ax2, TArrayD * axis2, Int_t Ax3, TArrayD* axis3, Int_t Ax4,Int_t Ax5){
  TH3D * hist3d = new TH3D(name,"",axis1->GetSize()-1,axis1->GetArray(),axis2->GetSize()-1,axis2->GetArray(),axis3->GetSize()-1,axis3->GetArray());
  TAxis* Axis1 = input->GetAxis(Ax1);
  hist3d->GetXaxis()->SetTitle(Axis1->GetTitle());
  TAxis* Axis2 = input->GetAxis(Ax2);
  hist3d->GetYaxis()->SetTitle(Axis2->GetTitle());
  TAxis* Axis3 = input->GetAxis(Ax3);
  hist3d->GetZaxis()->SetTitle(Axis3->GetTitle());
  TAxis* Axis4 = input->GetAxis(Ax4);
  TAxis* Axis5 = input->GetAxis(Ax5);
  hist3d->Sumw2();
  for(int x = 0;x<=Axis1->GetNbins()+1;x++){
    Double_t Ax1val = Axis1->GetBinCenter(x);
    for(int y = 0;y<=Axis2->GetNbins()+1;y++){
      Double_t Ax2val = Axis2->GetBinCenter(y);
      for(int z = 0;z<=Axis3->GetNbins()+1;z++){
	Double_t Ax3val = Axis3->GetBinCenter(z);
	for(int x1=1;x1<=Axis4->GetNbins();x1++){
	  for(int x2=1;x2<=Axis5->GetNbins();x2++){
	    Double_t content = GetBinAx(input,Ax1,x,Ax2,y,Ax3,z,Ax4,x1,Ax5,x2);
	    for(int i=1;i<=content;i++){
	      hist3d->Fill(Ax1val,Ax2val,Ax3val);
	    }
	  }
	}
      }
    }
  }
  return hist3d;
}

TH2D * Rebin(THnF * input, const char* name, Int_t Ax1, TArrayD * axis1,Int_t Ax2, TArrayD * axis2, Int_t Ax3, Int_t Ax5,Int_t Ax4, Double_t Ax4min){
  TH2D * hist2d = new TH2D(name,"",axis1->GetSize()-1,axis1->GetArray(),axis2->GetSize()-1,axis2->GetArray());
  TAxis* Axis1 = input->GetAxis(Ax1);
  hist2d->GetXaxis()->SetTitle(Axis1->GetTitle());
  TAxis* Axis2 = input->GetAxis(Ax2);
  hist2d->GetYaxis()->SetTitle(Axis2->GetTitle());
  TAxis* Axis3 = input->GetAxis(Ax3);
  TAxis* Axis4 = input->GetAxis(Ax4);
//   cout << Axis4->GetBinCenter(Axis4->GetNbins())<<endl;
  TAxis* Axis5 = input->GetAxis(Ax5);
  hist2d->Sumw2();
  for(int x = 0;x<=Axis1->GetNbins()+1;x++){
    Double_t Ax1val = Axis1->GetBinCenter(x);
    for(int y = 0;y<=Axis2->GetNbins()+1;y++){
      Double_t Ax2val = Axis2->GetBinCenter(y);
      for(int z = 0;z<=Axis3->GetNbins()+1;z++){
	for(int x1=1;x1<=Axis4->GetNbins();x1++){
	  Double_t Ax4val = Axis4->GetBinLowEdge(x1);
	  if(Ax4val<Ax4min){ continue;}
	  for(int x2=1;x2<=Axis5->GetNbins();x2++){
	    Double_t content = GetBinAx(input,Ax1,x,Ax2,y,Ax3,z,Ax4,x1,Ax5,x2);
	    for(int i=1;i<=content;i++){
	      hist2d->Fill(Ax1val,Ax2val);
	    }
	  }
	}
      }
    }
  }
  return hist2d;
}
  


TH1D * Errors(THnF * hist,const char* name,bool err = false){
  TH1D * errors = new TH1D(name,"distribution of errors in histogram in percent.",1000,0,100);
  errors->GetYaxis()->SetTitle("Relative error (%%)");
  for(int i = 1; i<=hist->GetAxis(0)->GetNbins();i++){
    for(int j = 1; j<=hist->GetAxis(1)->GetNbins();j++){
      for(int k = 1; k<=hist->GetAxis(2)->GetNbins();k++){
	for(int l = 1; l<=hist->GetAxis(3)->GetNbins();l++){
	  for(int m = 1; m<=hist->GetAxis(4)->GetNbins();m++){
	    Int_t index[5] = {i,j,k,l,m};
	    Double_t BinContent = hist->GetBinContent(index) ;
	    Double_t BinError =0.0;
	    if(err) BinError = hist->GetBinError2(hist->GetBin(index));
	    else BinError = TMath::Sqrt(BinContent);
	    if(BinError>10e-10){
	      Double_t RelErr = BinError/BinContent;
	      errors->Fill(RelErr*100);
	    }
	  }
	}
      }
    }
  }
  return errors;
}


TH2D * Slice(Int_t Ax1, Int_t Ax2,Int_t Ax3,Int_t Bin3,Int_t Ax4, Int_t Bin4, Int_t Ax5, Int_t Bin5,THnF * hist){
  //creates a 3d slice of a THnF:
  TH2D * slice = hist->Projection(Ax2,Ax1);
  slice->Sumw2(kFALSE);
  slice->Clear();
  for(int i = 0;i<=slice->GetNbinsX();i++){
    for(int j = 0;j<=slice->GetNbinsY();j++){
      slice->SetBinContent(i,j,GetBinAx(hist,Ax1,i,Ax2,j,Ax3,Bin3,Ax4,Bin4,Ax5,Bin5));
    }
  }
  slice->Sumw2();
  return slice;
}

TH2D * Slice(Int_t Ax1, Int_t Ax2,Int_t Ax3,Int_t Bin3,Int_t Ax4, Int_t Ax5,THnF * hist){
  //creates a 3d slice of a THnF:
  TH2D * slice = hist->Projection(Ax2,Ax1);
  slice->Sumw2(kFALSE);
  slice->Clear();
  for(int i = 0;i<=slice->GetNbinsX();i++){
    for(int j = 0;j<=slice->GetNbinsY();j++){
      Double_t content = 0.0;
      for(int k =1;k<=hist->GetAxis(Ax4)->GetNbins();k++){
	for(int l = 1;l<=hist->GetAxis(Ax5)->GetNbins();l++){
	  content += GetBinAx(hist,Ax1,i,Ax2,j,Ax3,Bin3,Ax4,k,Ax5,l);
	}
      }
    slice->SetBinContent(i,j,content);
    }
  }
  slice->Sumw2();
  return slice;
}


void MakeFoldersFill2dSlices(TFile* outfile, THnF * histpp, THnF * histmc){
  if(!outfile||!histpp||!histmc) return;
  TAxis * axis1 = histpp->GetAxis(0);
//   TString * axis1name = new TString("Cent");
//   TAxis * axis2 = histpp->GetAxis(1);
//   TString * axis1name = new TString("ZVert");
  TAxis * axis3 = histpp->GetAxis(2);
//   TString * axis1name = new TString("Phi");
//   TAxis * axis4 = histpp->GetAxis(3);
//   TString * axis1name = new TString("Eta");  
  TAxis * axis5 = histpp->GetAxis(4);
//   TString * axis1name = new TString("Pt");
  TDirectory * slicedir = outfile->mkdir("Slices");
  //z vertex vs eta vs pT:
  TDirectory * vzetadir =   slicedir->mkdir("ZVert_vs_Eta");
  TCanvas * canvas = new TCanvas("eff");
  for(int Cent = 1;Cent<= axis1->GetNbins();Cent++){
    TDirectory * CentDir = vzetadir->mkdir(Form("%4.2f_Cent_%4.2f",axis1->GetBinLowEdge(Cent),axis1->GetBinUpEdge(Cent)));
    for(int Phi = 1;Phi<=axis3->GetNbins();Phi++){
      TDirectory * PhiDir = CentDir->mkdir(Form("%4.2f_Phi_%4.2f",axis3->GetBinLowEdge(Phi),axis3->GetBinUpEdge(Phi)));
      for(int pT = 1;pT<=axis5->GetNbins();pT++){
	TDirectory * pTDir = PhiDir->mkdir(Form("%4.2f_pT_%4.2f",axis5->GetBinLowEdge(pT),axis5->GetBinUpEdge(pT)));
	pTDir->cd();
	TH2D * histpp2d = Slice(1,3,4,pT,0,Cent,2,Phi,histpp);
	histpp2d->Write("z_vert_vs_eta_dat");
	TH2D * histmc2d = Slice(1,3,4,pT,0,Cent,2,Phi,histmc);    
	histmc2d->Write("z_vert_vs_eta_vs_pT_mc");
	histmc2d->Divide(histpp2d);
	histmc2d->Write("z_vert_vs_eta_vs_pT_mc_div_dat");
	canvas->cd();
	histmc2d->SetStats(false);
	histmc2d->Draw("colz");
	pTDir->cd();
	canvas->Write();
	canvas->Clear();
	delete histpp2d;delete histmc2d;
      }
    }
  }
  
}

void CompareSliceInt(TDirectory * savedir,THnF * histrec,THnF * histmc){
  savedir->cd();
  if(TString(savedir->GetName()).CompareTo("phiptvscent")==0){
    TH2D * integrated = histrec->Projection(4,2);
    integrated->Divide(histmc->Projection(4,2));
    TCanvas * canvas = new TCanvas("phiptvscent");
    canvas->Divide(5,2);
    for(int mult = 1;mult<=histrec->GetAxis(0)->GetNbins();mult++){
      TH2D* slice = Slice(2,4,0,mult,1,3,histrec);
      slice->Divide(Slice(2,4,0,mult,1,3,histmc));
      slice->Divide(integrated);
      slice->SetTitle(Form("%4.2f_Cent_%4.2f",histrec->GetAxis(0)->GetBinLowEdge(mult),histrec->GetAxis(0)->GetBinUpEdge(mult)));
      slice->Write(Form("%4.2f_Cent_%4.2f",histrec->GetAxis(0)->GetBinLowEdge(mult),histrec->GetAxis(0)->GetBinUpEdge(mult)));
      canvas->cd(mult);
      slice->Draw("colz");
    }
    canvas->Update();
    canvas->SetTitle("Phi vs pT in different slices in Centrality divided by the integral over all centralities.");
    canvas->Write();
  }
  if(TString(savedir->GetName()).CompareTo("etavzvscent")==0){
    TH2D * integrated = histrec->Projection(3,1);
    integrated->Divide(histmc->Projection(3,1));
    TCanvas * canvas = new TCanvas("etavzvscent");
    canvas->Divide(5,2);
    for(int mult = 1;mult<=histrec->GetAxis(0)->GetNbins();mult++){
      TH2D* slice = Slice(1,3,0,mult,4,2,histrec);
      slice->Divide(Slice(1,3,0,mult,4,2,histmc));
      slice->Divide(integrated);
      slice->SetTitle(Form("%4.2f_Cent_%4.2f",histrec->GetAxis(0)->GetBinLowEdge(mult),histrec->GetAxis(0)->GetBinUpEdge(mult)));
      slice->Write(Form("%4.2f_Cent_%4.2f",histrec->GetAxis(0)->GetBinLowEdge(mult),histrec->GetAxis(0)->GetBinUpEdge(mult)));
      canvas->cd(mult);
      slice->Draw("colz");
    }
    canvas->Update();
    canvas->SetTitle("Eta vs Vz in different slices in Centrality divided by the integral over all centralities.");
    canvas->Write();
    
  }
}

void TestEffDependence(TDirectory * savedir,THnF * histrec, THnF * histmc,TArrayD * multax,TArrayD* vzaxis, TArrayD* phiaxis, TArrayD* etaaxis, TArrayD* ptaxis){
  //Test the dependence of variables on each other for different binnings.
  savedir->cd();
  THnF* histrebinrec = Rebin(histrec, Form("%s_recthn",savedir->GetName()),multax,vzaxis,phiaxis,etaaxis,ptaxis);
  histrebinrec->Write();
  THnF* histrebinmc  = Rebin(histmc, Form("%s_mcthn",savedir->GetName()),multax,vzaxis,phiaxis,etaaxis,ptaxis);
  histrebinmc->Write();
  TDirectory * phiptvzcent = savedir->mkdir("phiptvscent");
  CompareSliceInt(phiptvzcent,histrebinrec,histrebinmc);

  TDirectory * etavzvscent = savedir->mkdir("etavzvscent");
  CompareSliceInt(etavzvscent,histrebinrec,histrebinmc);
  
}

void MakeTestHists(){
  //open the input file:
  TFile * infile = TFile::Open("AnalysisResults.root","READ");
  THnF * hnTracksInBins;
  THnF * hnTracksInBinsMC;
  TList * dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")->Get("ThreePartTrackEfficienciesPbPb_0_16_0_0Coutput1"));
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBins = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMC = dynamic_cast<THnF*>(dir->At(i));
  }
  TFile* outfile = TFile::Open("eff.root","RECREATE");
  TAxis * multaxis = hnTracksInBins->GetAxis(0);
  multaxis->SetTitle("Centrality [%]");
  hnTracksInBinsMC->GetAxis(0)->SetTitle("Centrality [%]");
  TAxis * vzAxis = hnTracksInBins->GetAxis(1);
  vzAxis->SetTitle("Z-Vertex [cm]");  
  hnTracksInBinsMC->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  TAxis * phiaxis = hnTracksInBins->GetAxis(2);
  phiaxis->SetTitle("#Phi [radians]");    
  hnTracksInBinsMC->GetAxis(2)->SetTitle("#Phi [radians]");
  TAxis * Etaaxis = hnTracksInBins->GetAxis(3);
  Etaaxis->SetTitle("#eta []");  
  hnTracksInBinsMC->GetAxis(3)->SetTitle("#eta []");
  TAxis * pTaxis = hnTracksInBins->GetAxis(4);
  pTaxis->SetTitle("pT [GeV/c]");    
  hnTracksInBinsMC->GetAxis(4)->SetTitle("pT [GeV/c]");
  outfile->cd();  
  //Binning 1:
  Double_t multaxisd[11] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,80.0,95.0} ;
  TArrayD * multaxisar = new TArrayD(11,multaxisd);
  Double_t  vzaxisd[8] = {-10.0,-8.0,-6.0,-2.0,2.0,6.0,8.0,10.0};
  TArrayD * vzaxisar = new TArrayD(8,vzaxisd);
  TArrayD * phiar = new TArrayD(phiaxis->GetNbins()+1);
  for(int i = 1;i<=phiaxis->GetNbins();i++){phiar->AddAt(phiaxis->GetBinLowEdge(i),i-1);}
  phiar->AddAt(phiaxis->GetBinUpEdge(phiaxis->GetNbins()),phiaxis->GetNbins());
  TArrayD * etaar = new TArrayD(Etaaxis->GetNbins()+1);
  for(int i = 1;i<=Etaaxis->GetNbins();i++){etaar->AddAt(Etaaxis->GetBinLowEdge(i),i-1);}
  etaar->AddAt(Etaaxis->GetBinUpEdge(Etaaxis->GetNbins()),Etaaxis->GetNbins());
  Double_t  pTaxisd[31] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,5.0,8.0,16.0};
  TArrayD * pTaxisar = new TArrayD(31,pTaxisd);
 
  TestEffDependence(outfile->mkdir("binning1"),hnTracksInBins,hnTracksInBinsMC,multaxisar,vzaxisar,phiar,etaar,pTaxisar);

}

void MakeEffHistspp(const char* file = "AnalysisResults.root"){
   //open the input file:
  TString filename = TString(file);
  //test if local or global
  bool global = false;
  if(filename.CompareTo("AnalysisResults.root")==0) global = true;
  TFile * infile = TFile::Open(filename.Data(),"READ");
  THnF * hnTracksInBins;
  THnF * hnTracksInBinsMC;
  TH2D * hnCentVsVertex;
  TList * dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficiencies_0_16")->Get("ThreePartTrackEfficiencies_0_16_0_0Coutput1"));
  TFile* outfile = TFile::Open("eff.root","RECREATE");
  outfile->cd();
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBins = dynamic_cast<THnF*>(dir->At(i)->Clone("hnTracksinBins_eff"));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMC = dynamic_cast<THnF*>(dir->At(i)->Clone("hnTracksinBinsMC_eff"));
    if(TString(dir->At(i)->GetName()).CompareTo("centVsZVertex")==0) hnCentVsVertex = dynamic_cast<TH2D*>(dir->At(i)->Clone("centVsZVertex_eff"));
  }
  infile->Close();
  //Axis: 0 - Centrality or multiplicity, 1 - Vertex , 2 - phi , 3 - eta , 4 - pT
  TAxis * multaxis = hnTracksInBins->GetAxis(0);
  multaxis->SetTitle("Centrality [%]");
  hnTracksInBinsMC->GetAxis(0)->SetTitle("Centrality [%]");
  TAxis * vzAxis = hnTracksInBins->GetAxis(1);
  vzAxis->SetTitle("Z-Vertex [cm]");  
  hnTracksInBinsMC->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  TAxis * phiaxis = hnTracksInBins->GetAxis(2);
  phiaxis->SetTitle("#Phi [radians]");    
  hnTracksInBinsMC->GetAxis(2)->SetTitle("#Phi [radians]");
  TAxis * Etaaxis = hnTracksInBins->GetAxis(3);
  Etaaxis->SetTitle("#eta []");  
  hnTracksInBinsMC->GetAxis(3)->SetTitle("#eta []");
  TAxis * pTaxis = hnTracksInBins->GetAxis(4);
  pTaxis->SetTitle("pT [GeV/c]");    
  hnTracksInBinsMC->GetAxis(4)->SetTitle("pT [GeV/c]");
  outfile->cd();
  hnTracksInBins->Write();
  hnTracksInBinsMC->Write();
  Errors(hnTracksInBins,"Errorsraw")->Write();
  Errors(hnTracksInBinsMC,"ErrorsMC")->Write();
  hnCentVsVertex->Write();
  
  //1D projections
  TDirectory * dim1 = outfile->mkdir("Projections_1d");
  dim1->cd();
  
  TH1D * multRec = hnTracksInBins->Projection(0);
  TH1D * multMC = hnTracksInBinsMC->Projection(0);
  multRec->Sumw2();
  multRec->Write("Multiplicity_All_Reconstructed");
  multMC->Sumw2();
  multMC->Write("Multiplicity_All_Produced");
  TH1D * multEffRec = dynamic_cast<TH1D*>(multRec->Clone("Efficiency_Mult_Rec"));
  multEffRec->Divide(multMC);
  multEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRec->Write();
  delete multRec; delete multMC; delete multEffRec; 

  TH1D * vzRec = hnTracksInBins->Projection(1);
  TH1D * vzMC = hnTracksInBinsMC->Projection(1);
  vzRec->Sumw2();
  vzRec->Write("VZ_All_Reconstructed");
  vzMC->Sumw2();
  vzMC->Write("VZ_All_Produced");
  TH1D * vzEffRec = dynamic_cast<TH1D*>(vzRec->Clone("Efficiency_VZ_Rec"));
  vzEffRec->Divide(vzMC);
  vzEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRec->Write();
  delete vzRec; delete vzMC;delete vzEffRec;

  TH1D * phiRec = hnTracksInBins->Projection(2);
  TH1D * phiMC = hnTracksInBinsMC->Projection(2);
  phiRec->Sumw2();
  phiRec->Write("phi_All_Reconstructed");
  phiMC->Sumw2();
  phiMC->Write("phi_All_Produced");
  TH1D * phiEffRec = dynamic_cast<TH1D*>(phiRec->Clone("Efficiency_phi_Rec"));
  phiEffRec->Divide(phiMC);
  phiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiEffRec->Write();
  delete phiRec; delete phiMC;delete phiEffRec;

  TH1D * etaRec = hnTracksInBins->Projection(3);
  TH1D * etaMC = hnTracksInBinsMC->Projection(3);
  etaRec->Sumw2();
  etaRec->Write("eta_All_Reconstructed");
  etaMC->Sumw2();
  etaMC->Write("eta_All_Produced");
  TH1D * etaEffRec = dynamic_cast<TH1D*>(etaRec->Clone("Efficiency_eta_Rec"));
  etaEffRec->Divide(etaMC);
  etaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaEffRec->Write();
  delete etaRec; delete etaMC;delete etaEffRec;

  TH1D * pTRec = hnTracksInBins->Projection(4);
  TH1D * pTMC = hnTracksInBinsMC->Projection(4);
  pTRec->Sumw2();
  pTRec->Write("pT_All_Reconstructed");
  pTMC->Sumw2();
  pTMC->Write("pT_All_Produced");
  TH1D * pTEffRec = dynamic_cast<TH1D*>(pTRec->Clone("Efficiency_pT_Rec"));
  pTEffRec->Divide(pTMC);
  pTEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  pTEffRec->Write();
  TH1D * pTWeightRec = dynamic_cast<TH1D*>(pTMC->Clone("Weight_pT_Rec"));
  pTWeightRec->Divide(pTRec);
  pTWeightRec->SetTitle(" produced MC particles divided by Reconstructed tracks");
  delete pTRec; 
  delete pTMC;delete pTEffRec;

  
  //2D projections
  TDirectory * dim2 = outfile->mkdir("Projections_2d");
  dim2->cd();
  TCanvas * can2d = new TCanvas("canvas");
  
  TH2D * multvzrec = hnTracksInBins->Projection(1,0);
  multvzrec->Sumw2();
  multvzrec->Write("mult_VZ_All_Reconstructed");
  can2d->cd();
  multvzrec->Draw("colz");
  can2d->Write("mult_VZ_All_Reconstructed_can");
  can2d->Clear();
  TH2D * multvzMC = hnTracksInBinsMC->Projection(1,0);
  multvzMC->Sumw2();
  multvzMC->Write("mult_VZ_All_Produced");    
  can2d->cd();
  multvzMC->Draw("colz");
  can2d->Write("mult_VZ_All_Produced_can");
  can2d->Clear();
  TH2D * multvzeffrec = dynamic_cast<TH2D*>(multvzrec->Clone("Efficiency_mult_vz_Rec"));
  multvzeffrec->Divide(multvzMC);
  multvzeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multvzeffrec->Write();
  can2d->cd();
  multvzeffrec->Draw("colz");
  can2d->Write("Efficiency_mult_vz_Rec_can");
  can2d->Clear();
  delete multvzrec; delete multvzMC;delete multvzeffrec;

  TH2D * multphirec = hnTracksInBins->Projection(2,0);
  multphirec->Sumw2();
  multphirec->Write("mult_phi_All_Reconstructed");
  can2d->cd();
  multphirec->Draw("colz");
  can2d->Write("mult_phi_All_Reconstructed_can");
  can2d->Clear();
  TH2D * multphiMC = hnTracksInBinsMC->Projection(2,0);
  multphiMC->Sumw2();
  multphiMC->Write("mult_phi_All_Produced");    
  can2d->cd();
  multphiMC->Draw("colz");
  can2d->Write("mult_phi_All_Produced_can");
  can2d->Clear();
  TH2D * multphieffrec = dynamic_cast<TH2D*>(multphirec->Clone("Efficiency_mult_phi_Rec"));
  multphieffrec->Divide(multphiMC);
  multphieffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multphieffrec->Write();
  can2d->cd();
  multphieffrec->Draw("colz");
  can2d->Write("Efficiency_mult_phi_Rec_can");
  can2d->Clear();
  delete multphirec; delete multphiMC;delete multphieffrec;

  TH2D * multetarec = hnTracksInBins->Projection(3,0);
  multetarec->Sumw2();
  multetarec->Write("mult_eta_All_Reconstructed");
  can2d->cd();
  multetarec->Draw("colz");
  can2d->Write("mult_eta_All_Reconstructed_can");
  can2d->Clear();
  TH2D * multetaMC = hnTracksInBinsMC->Projection(3,0);
  multetaMC->Sumw2();
  multetaMC->Write("mult_eta_All_Produced");    
  can2d->cd();
  multetaMC->Draw("colz");
  can2d->Write("mult_eta_All_Produced_can");
  can2d->Clear();
  TH2D * multetaeffrec = dynamic_cast<TH2D*>(multetarec->Clone("Efficiency_mult_eta_Rec"));
  multetaeffrec->Divide(multetaMC);
  multetaeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multetaeffrec->Write();
  can2d->cd();
  multetaeffrec->Draw("colz");
  can2d->Write("Efficiency_mult_eta_Rec_can");
  can2d->Clear();
  delete multetarec; delete multetaMC;delete multetaeffrec;

  TH2D * multptRec = hnTracksInBins->Projection(4,0);
  multptRec->Sumw2();
  multptRec->Write("mult_pT_All_Reconstructed");
  can2d->cd();
  multptRec->Draw("colz");
  can2d->Write("mult_pT_All_Reconstructed_can");
  can2d->Clear();
  TH2D * multptMC = hnTracksInBinsMC->Projection(4,0);
  multptMC->Sumw2();
  multptMC->Write("mult_pT_All_Produced");
  can2d->cd();
  multptMC->Draw("colz");
  can2d->Write("mult_pT_All_Produced_can");
  can2d->Clear();
  TH2D * multptEffRec = dynamic_cast<TH2D*>(multptRec->Clone("Efficiency_mult_pT_Rec"));
  multptEffRec->Divide(multptMC);
  multptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multptEffRec->Write();
  can2d->cd();
  multptEffRec->Draw("colz");
  can2d->Write("Efficiency_mult_pT_Rec_can");
  can2d->Clear();
  delete multptRec; delete multptMC;delete multptEffRec;

  TH2D * vzphiRec = hnTracksInBins->Projection(2,1);
  vzphiRec->Sumw2();
  vzphiRec->Write("VZ_phi_All_Reconstructed");
  can2d->cd();
  vzphiRec->Draw("colz");
  can2d->Write("VZ_phi_All_Reconstructed_can");
  can2d->Clear();
  TH2D * vzphiMC = hnTracksInBinsMC->Projection(2,1);
  vzphiMC->Sumw2();
  vzphiMC->Write("VZ_phi_All_Produced");
  can2d->cd();
  vzphiMC->Draw("colz");
  can2d->Write("VZ_phi_All_Produced_can");
  can2d->Clear();
  TH2D * vzphiEffRec = dynamic_cast<TH2D*>(vzphiRec->Clone("Efficiency_vz_phi_Rec"));
  vzphiEffRec->Divide(vzphiMC);
  vzphiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzphiEffRec->Write();
  can2d->cd();
  vzphiEffRec->Draw("colz");
  can2d->Write("Efficiency_vz_phi_Rec_can");
  can2d->Clear();
  delete vzphiRec;delete vzphiMC;delete vzphiEffRec;
  
  TH2D * vzetaRec = hnTracksInBins->Projection(3,1);
  vzetaRec->Sumw2();
  vzetaRec->Write("VZ_Eta_All_Reconstructed");
  TH2D * vzetaMC = hnTracksInBinsMC->Projection(3,1);
  vzetaMC->Sumw2();
  vzetaMC->Write("VZ_Eta_All_Produced");    
  TH2D * vzetaEffRec = dynamic_cast<TH2D*>(vzetaRec->Clone("Efficiency_vz_eta_Rec"));
  vzetaEffRec->Divide(vzetaMC);
  vzetaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzetaEffRec->Write();
  delete vzetaRec;delete vzetaMC;delete vzetaEffRec;
  
  TH2D * vzptRec = hnTracksInBins->Projection(4,1);
  vzptRec->Sumw2();
  vzptRec->Write("VZ_pT_All_Reconstructed");
  TH2D * vzptMC = hnTracksInBinsMC->Projection(4,1);
  vzptMC->Sumw2();
  vzptMC->Write("VZ_pT_All_Produced");    
  TH2D * vzptEffRec = dynamic_cast<TH2D*>(vzptRec->Clone("Efficiency_vz_pT_Rec"));
  vzptEffRec->Divide(vzptMC);
  vzptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzptEffRec->Write();
  delete vzptRec;  delete vzptMC;delete vzptEffRec;

  TH2D * phietaRec = hnTracksInBins->Projection(3,2);
  phietaRec->Sumw2();
  phietaRec->Write("phi_Eta_All_Reconstructed");
  TH2D * phietaMC = hnTracksInBinsMC->Projection(3,2);
  phietaMC->Sumw2();
  phietaMC->Write("phi_Eta_All_Produced");    
  TH2D * phietaEffRec = dynamic_cast<TH2D*>(phietaRec->Clone("Efficiency_phi_eta_Rec"));
  phietaEffRec->Divide(phietaMC);
  phietaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phietaEffRec->Write();
  delete phietaRec;  delete phietaMC;delete phietaEffRec;

  TH2D * phiptRec = hnTracksInBins->Projection(4,2);
  phiptRec->Sumw2();
  phiptRec->Write("phi_pT_All_Reconstructed");
  TH2D * phiptMC = hnTracksInBinsMC->Projection(4,2);
  phiptMC->Sumw2();
  phiptMC->Write("phi_pT_All_Produced");    
  TH2D * phiptEffRec = dynamic_cast<TH2D*>(phiptRec->Clone("Efficiency_phi_pT_Rec"));
  phiptEffRec->Divide(phiptMC);
  phiptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiptEffRec->Write();
  delete phiptRec; delete phiptMC;delete phiptEffRec;
  
  TH2D * etaptRec = hnTracksInBins->Projection(4,3);
  etaptRec->Sumw2();
  etaptRec->Write("eta_pT_All_Reconstructed");
  TH2D * etaptMC = hnTracksInBinsMC->Projection(4,3);
  etaptMC->Sumw2();
  etaptMC->Write("eta_pT_All_Produced");    
  TH2D * etaptEffRec = dynamic_cast<TH2D*>(etaptRec->Clone("Efficiency_eta_pT_Rec"));
  etaptEffRec->Divide(etaptMC);
  etaptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaptEffRec->Write();  
  delete etaptRec; delete etaptMC;delete etaptEffRec;

  //make more meaningful axes:
  outfile->cd();
  //mult: keep for now, unsure how it actually should look.
  Double_t multaxisArray[2] = {0.0,300};
  TArrayD * multaxisA = new TArrayD(2,multaxisArray);
  //vz: 4 cm bins from -6 to 6,6-8,8-10 and the same on the other side:
  Double_t  vzaxisArray[8] = {-10.0,-7.5,-5.0,-2.0,2.0,5.0,7.5,10.0};
  TArrayD * vzaxisA = new TArrayD(8,vzaxisArray);//   
  
  Double_t etamin = -0.9;
  Double_t etamax = 0.9;
  Double_t deta = (etamax - etamin)/63;
  TArrayD * etaaxisAD = new TArrayD(63+1);
  for(int i = 0; i<=63;i++){
    etaaxisAD->AddAt(etamin+i*deta,i);
  }
  
  //make more meaningful axes:
  Double_t  pTaxisArray[16] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  TArrayD * pTaxisA = new TArrayD(16,pTaxisArray);

  TH3D * histrec = Rebin(hnTracksInBins,"hist3drec",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histrec->Write();
  TH3D * histMC = Rebin(hnTracksInBinsMC,"hist3dMC",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histMC->Write();
  histMC->Divide(histrec);
  histMC->Write("Eff3d");
  TFile* outfile2 = TFile::Open("LHC11aWeight.root","RECREATE");
  outfile2->cd();
  histMC->Write("hnWeight");
//   outfile2->Close();
  outfile->cd();
  TH1D * err = new TH1D("err","Errors in the 3d histogram",100,0.0,20.0);
  for(int x=1;x<=histrec->GetNbinsX();x++){
    for(int y=1;y<=histrec->GetNbinsY();y++){
      for(int z=1;z<=histrec->GetNbinsZ();z++){
	if(histrec->GetBinContent(x,y,z)>1.0e-10)err->Fill(100*histrec->GetBinError(x,y,z)/histrec->GetBinContent(x,y,z));
      }
    }
  }
  err->Write();
  TH2D * histrechpT = Rebin(hnTracksInBins,"hist3drechpT",3,etaaxisAD,1,vzaxisA,0,2,4,2.0);
  histrechpT->Write();
  TH2D * histhpTMC = Rebin(hnTracksInBinsMC,"hist3dMChpT",3,etaaxisAD,1,vzaxisA,0,2,4,2.0);
  histhpTMC->Write();  
  histhpTMC->Divide(histrechpT);
  Double_t nbins = 0.0;
  Double_t content = 0.0;
  for(int x = 1;x<=histhpTMC->GetXaxis()->GetNbins();x++){
     for(int y = 1;y<=histhpTMC->GetYaxis()->GetNbins();y++){
	if(histhpTMC->GetBinContent(x,y)>1.0E-10){
	  nbins +=1.0;
	  content+=histhpTMC->GetBinContent(x,y);       
	}
      }
    }
  histhpTMC->Write("Eff2dhpT");
  TF1* ptfunc;
  outfile2->cd();
  histhpTMC->Write("hnWeight_highpt");
  ptfunc = new TF1("pT_function","[0]+[1]*x+[2]*x*x+[3]*exp(-x)",2,100);
  pTWeightRec->Scale(nbins/content);
  pTWeightRec->Fit(ptfunc,"","",2,100);
  outfile->cd();
  pTWeightRec->Write();
  outfile2->cd();
  ptfunc->Write();
  outfile2->Close();  
  delete hnTracksInBins; delete hnTracksInBinsMC; delete hnCentVsVertex;    
  delete pTWeightRec;
  delete histrec; delete histMC;delete err;delete histrechpT;delete histhpTMC;
  delete ptfunc;
  delete multaxisA; delete vzaxisA; delete pTaxisA;  
  

  outfile->Close();

  
//   infile->Close();
// //   delete multEffRec;delete multEffRecPP;
}

void MakeEffHistsPbPb(const char* file = "AnalysisResults.root"){
  //open the input file:
  TString filename = TString(file);
  //test if local or global
  bool global = false;
  if(filename.CompareTo("AnalysisResults.root")==0) global = true;
  TFile * infile = TFile::Open(filename.Data(),"READ");
  if(!infile)return;
  
  THnF * hnTracksInBinstmp;  THnF * hnTracksInBins;
//   THnF * hnTracksInBinsRecPPtmp;THnF * hnTracksInBinsRecPP;
  THnF * hnTracksInBinsMCtmp;  THnF * hnTracksInBinsMC;
  TH2D * hnCentVsVertextmp;  TH2D * hnCentVsVertex;
  TList * dir;
  if(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")){
   dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")->Get("ThreePartTrackEfficienciesPbPb_0_16_0_0Coutput1"));
  }
  else if(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_100")){
   dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_100")->Get("ThreePartTrackEfficienciesPbPb_0_100_0_0Coutput1"));    
  }
  else return;
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBinstmp = dynamic_cast<THnF*>(dir->At(i));
//     if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsRecPP")==0)hnTracksInBinsRecPPtmp = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMCtmp = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("centVsZVertex")==0) hnCentVsVertextmp = (dynamic_cast<TH2D*>(dir->At(i)));
  }

  TFile* outfile = TFile::Open(filename.ReplaceAll("AnalysisResults.root","eff.root"),"RECREATE");
  outfile->cd();
  if(filename.Contains("LHC12a17g")){
    hnTracksInBins = SumCent(0.0,10.0,hnTracksInBinstmp,"hnTracksReconstruced");
//     hnTracksInBinsRecPP = SumCent(0.0,10.0,hnTracksInBinsRecPPtmp,"hnTracksReconstruced_PP");
    hnTracksInBinsMC = SumCent(0.0,10.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(0.0,10.0,hnCentVsVertextmp,"hnCentralityvsVertex");
    cout << "hijhoh5"<<endl;
  }
  if(filename.Contains("LHC12a17h")){
    hnTracksInBins = SumCent(10.0,50.0,hnTracksInBinstmp,"hnTracksReconstruced");
//     hnTracksInBinsRecPP = SumCent(10.0,50.0,hnTracksInBinsRecPPtmp,"hnTracksReconstruced_PP");
    hnTracksInBinsMC = SumCent(10.0,50.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(10.0,50.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  if(filename.Contains("LHC12a17i")){
    hnTracksInBins = SumCent(50.0,100.0,hnTracksInBinstmp,"hnTracksReconstruced");
//     hnTracksInBinsRecPP = SumCent(50.0,100.0,hnTracksInBinsRecPPtmp,"hnTracksReconstruced_PP");
    hnTracksInBinsMC = SumCent(50.0,100.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(50.0,100.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  if(filename.Contains("LHC11a10a")){
    hnTracksInBins = SumCent(0.0,100.0,hnTracksInBinstmp,"hnTracksReconstruced");
//     hnTracksInBinsRecPP = SumCent(0.0,100.0,hnTracksInBinsRecPPtmp,"hnTracksReconstruced_PP");
    hnTracksInBinsMC = SumCent(0.0,100.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(0.0,100.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  infile->Close();
  cout << "hoiho2"<<endl;
  delete hnTracksInBinstmp; delete hnTracksInBinsMCtmp;
//   delete hnTracksInBinsRecPPtmp;
  delete hnCentVsVertextmp;
  cout << "tres"<<endl;

  //Axis: 0 - Centrality or multiplicity, 1 - Vertex , 2 - phi , 3 - eta , 4 - pT

  TAxis * multaxis = hnTracksInBins->GetAxis(0);
  multaxis->SetTitle("Centrality [%]");
//   hnTracksInBinsRecPP->GetAxis(0)->SetTitle("Centrality [%]");
  hnTracksInBinsMC->GetAxis(0)->SetTitle("Centrality [%]");
  TAxis * vzAxis = hnTracksInBins->GetAxis(1);
  vzAxis->SetTitle("Z-Vertex [cm]");  
//   hnTracksInBinsRecPP->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  hnTracksInBinsMC->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  TAxis * phiaxis = hnTracksInBins->GetAxis(2);
  phiaxis->SetTitle("#Phi [radians]");    
//   hnTracksInBinsRecPP->GetAxis(2)->SetTitle("#Phi [radians]");
  hnTracksInBinsMC->GetAxis(2)->SetTitle("#Phi [radians]");
  TAxis * Etaaxis = hnTracksInBins->GetAxis(3);
  Etaaxis->SetTitle("#eta []");  
//   hnTracksInBinsRecPP->GetAxis(3)->SetTitle("#eta []");
  hnTracksInBinsMC->GetAxis(3)->SetTitle("#eta []");
  TAxis * pTaxis = hnTracksInBins->GetAxis(4);
  pTaxis->SetTitle("pT [GeV/c]");    
//   hnTracksInBinsRecPP->GetAxis(4)->SetTitle("pT [GeV/c]");
  hnTracksInBinsMC->GetAxis(4)->SetTitle("pT [GeV/c]");
  outfile->cd();
  hnTracksInBins->Write();
  hnTracksInBinsMC->Write();
  Errors(hnTracksInBins,"Errorsraw")->Write();
  Errors(hnTracksInBinsMC,"ErrorsMC")->Write();
  hnCentVsVertex->Write();
  //1D projections
  TDirectory * dim1 = outfile->mkdir("Projections_1d");
  dim1->cd();
  
  TH1D * multRec = hnTracksInBins->Projection(0);
//   TH1D * multRecPP = hnTracksInBinsRecPP->Projection(0);
  TH1D * multMC = hnTracksInBinsMC->Projection(0);
  multRec->Sumw2();
  multRec->Write("Multiplicity_All_Reconstructed");
//   multRecPP->Sumw2();
//   multRecPP->Write("Multiplicity_Reconstructed_from_Physical_Primary");  
  multMC->Sumw2();
  multMC->Write("Multiplicity_All_Produced");  
  TH1D * multEffRec = dynamic_cast<TH1D*>(multRec->Clone("Efficiency_Mult_Rec"));
  multEffRec->Divide(multMC);
  multEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRec->Write();
//   TH1D * multEffRecPP = dynamic_cast<TH1D*>(multRecPP->Clone("Efficiency_Mult_Rec_pp"));
//   multEffRecPP->Divide(multMC);
//   multEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   multEffRecPP->Write();
  delete multRec; 
//   delete multRecPP; 
  delete multMC;delete multEffRec;
//   delete multEffRecPP;


  TH1D * vzRec = hnTracksInBins->Projection(1);
//   TH1D * vzRecPP = hnTracksInBinsRecPP->Projection(1);
  TH1D * vzMC = hnTracksInBinsMC->Projection(1);
  vzRec->Sumw2();
  vzRec->Write("VZ_All_Reconstructed");
//   vzRecPP->Sumw2();
//   vzRecPP->Write("VZ_Reconstructed_from_Physical_Primary");  
  vzMC->Sumw2();
  vzMC->Write("VZ_All_Produced");
  TH1D * vzEffRec = dynamic_cast<TH1D*>(vzRec->Clone("Efficiency_VZ_Rec"));
  vzEffRec->Divide(vzMC);
  vzEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRec->Write();
//   TH1D * vzEffRecPP = dynamic_cast<TH1D*>(vzRecPP->Clone("Efficiency_VZ_Rec_pp"));
//   vzEffRecPP->Divide(vzMC);
//   vzEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   vzEffRecPP->Write();
  delete vzRec; 
//   delete vzRecPP; 
  delete vzMC;delete vzEffRec;
//   delete   vzEffRecPP;
  
  TH1D * phiRec = hnTracksInBins->Projection(2);
//   TH1D * phiRecPP = hnTracksInBinsRecPP->Projection(2);
  TH1D * phiMC = hnTracksInBinsMC->Projection(2);
  phiRec->Sumw2();
  phiRec->Write("phi_All_Reconstructed");
//   phiRecPP->Sumw2();
//   phiRecPP->Write("phi_Reconstructed_from_Physical_Primary");  
  phiMC->Sumw2();
  phiMC->Write("phi_All_Produced");
  TH1D * phiEffRec = dynamic_cast<TH1D*>(phiRec->Clone("Efficiency_phi_Rec"));
  phiEffRec->Divide(phiMC);
  phiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiEffRec->Write();
//   TH1D * phiEffRecPP = dynamic_cast<TH1D*>(phiRecPP->Clone("Efficiency_phi_Rec_pp"));
//   phiEffRecPP->Divide(phiMC);
//   phiEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   phiEffRecPP->Write();
  delete phiRec; 
//   delete phiRecPP; 
  delete phiMC;delete phiEffRec;
//   delete phiEffRecPP;
  
  TH1D * etaRec = hnTracksInBins->Projection(3);
//   TH1D * etaRecPP = hnTracksInBinsRecPP->Projection(3);
  TH1D * etaMC = hnTracksInBinsMC->Projection(3);
  etaRec->Sumw2();
  etaRec->Write("eta_All_Reconstructed");
//   etaRecPP->Sumw2();
//   etaRecPP->Write("eta_Reconstructed_from_Physical_Primary");  
  etaMC->Sumw2();
  etaMC->Write("eta_All_Produced");  
  TH1D * etaEffRec = dynamic_cast<TH1D*>(etaRec->Clone("Efficiency_eta_Rec"));
  etaEffRec->Divide(etaMC);
  etaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaEffRec->Write();
//   TH1D * etaEffRecPP = dynamic_cast<TH1D*>(etaRecPP->Clone("Efficiency_eta_Rec_pp"));
//   etaEffRecPP->Divide(etaMC);
//   etaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   etaEffRecPP->Write();
  delete etaRec; 
//   delete etaRecPP; 
  delete etaMC;delete etaEffRec;
//   delete etaEffRecPP;

  TH1D * pTRec = hnTracksInBins->Projection(4);
//   TH1D * pTRecPP = hnTracksInBinsRecPP->Projection(4);
  TH1D * pTMC = hnTracksInBinsMC->Projection(4);
  pTRec->Sumw2();
  pTRec->Write("pT_All_Reconstructed");
//   pTRecPP->Sumw2();
//   pTRecPP->Write("pT_Reconstructed_from_Physical_Primary");  
  pTMC->Sumw2();
  pTMC->Write("pT_All_Produced");
  TH1D * pTEffRec = dynamic_cast<TH1D*>(pTRec->Clone("Efficiency_pT_Rec"));
  pTEffRec->Divide(pTMC);
  pTEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  pTEffRec->Write();
//   TH1D * pTEffRecPP = dynamic_cast<TH1D*>(pTRecPP->Clone("Efficiency_pT_Rec_pp"));
//   pTEffRecPP->Divide(pTMC);
//   pTEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   pTEffRecPP->Write();
  TH1D * pTWeightRec = dynamic_cast<TH1D*>(pTMC->Clone("Weight_pT_Rec"));
  pTWeightRec->Divide(pTRec);
  pTWeightRec->SetTitle(" produced MC particles divided by Reconstructed tracks");
  delete pTRec; 
//   delete pTRecPP; 
  delete pTMC;delete pTEffRec;
//   delete pTEffRecPP;
  

  //2D projections
  TDirectory * dim2 = outfile->mkdir("Projections_2d");
  dim2->cd();
  
  TH2D * multvzrec = hnTracksInBins->Projection(1,0);
  multvzrec->Sumw2();
  multvzrec->Write("mult_VZ_All_Reconstructed");
/*  TH2D * multvzrecPP = hnTracksInBinsRecPP->Projection(1,0);
  multvzrecPP->Sumw2();
  multvzrecPP->Write("mult_VZ_All_Reconstructed_from_Physical_Primary"); */ 
  TH2D * multvzMC = hnTracksInBinsMC->Projection(1,0);
  multvzMC->Sumw2();
  multvzMC->Write("mult_VZ_All_Produced");    
  TH2D * multvzeffrec = dynamic_cast<TH2D*>(multvzrec->Clone("Efficiency_mult_vz_Rec"));
  multvzeffrec->Divide(multvzMC);
  multvzeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multvzeffrec->Write();
//   TH2D * multvzeffrecPP = dynamic_cast<TH2D*>(multvzrecPP->Clone("Efficiency_mult_vz_Rec_pp"));
//   multvzeffrecPP->Divide(multvzMC);
//   multvzeffrecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   multvzeffrecPP->Write();
  delete multvzrec; /*delete multvzrecPP;*/ delete multvzMC;delete multvzeffrec;/*delete multvzeffrecPP;*/

  TH2D * multphirec = hnTracksInBins->Projection(2,0);
  multphirec->Sumw2();
  multphirec->Write("mult_phi_All_Reconstructed");
//   TH2D * multphirecPP = hnTracksInBinsRecPP->Projection(2,0);
//   multphirecPP->Sumw2();
//   multphirecPP->Write("mult_phi_All_Reconstructed_from_Physical_Primary");  
  TH2D * multphiMC = hnTracksInBinsMC->Projection(2,0);
  multphiMC->Sumw2();
  multphiMC->Write("mult_phi_All_Produced");    
  TH2D * multphieffrec = dynamic_cast<TH2D*>(multphirec->Clone("Efficiency_mult_phi_Rec"));
  multphieffrec->Divide(multphiMC);
  multphieffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multphieffrec->Write();
//   TH2D * multphieffrecPP = dynamic_cast<TH2D*>(multphirecPP->Clone("Efficiency_mult_phi_Rec_pp"));
//   multphieffrecPP->Divide(multphiMC);
//   multphieffrecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   multphieffrecPP->Write();
  delete multphirec; /*delete multphirecPP;*/ delete multphiMC;delete multphieffrec;/*delete multphieffrecPP;*/

  TH2D * multetarec = hnTracksInBins->Projection(3,0);
  multetarec->Sumw2();
  multetarec->Write("mult_eta_All_Reconstructed");
//   TH2D * multetarecPP = hnTracksInBinsRecPP->Projection(3,0);
//   multetarecPP->Sumw2();
//   multetarecPP->Write("mult_eta_All_Reconstructed_from_Physical_Primary");  
  TH2D * multetaMC = hnTracksInBinsMC->Projection(3,0);
  multetaMC->Sumw2();
  multetaMC->Write("mult_eta_All_Produced");    
  TH2D * multetaeffrec = dynamic_cast<TH2D*>(multetarec->Clone("Efficiency_mult_eta_Rec"));
  multetaeffrec->Divide(multetaMC);
  multetaeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multetaeffrec->Write();
//   TH2D * multetaeffrecPP = dynamic_cast<TH2D*>(multetarecPP->Clone("Efficiency_mult_eta_Rec_pp"));
//   multetaeffrecPP->Divide(multetaMC);
//   multetaeffrecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   multetaeffrecPP->Write();
  delete multetarec; /*delete multetarecPP;*/ delete multetaMC;delete multetaeffrec;/*delete multetaeffrecPP;*/

  TH2D * multptRec = hnTracksInBins->Projection(4,0);
  multptRec->Sumw2();
  multptRec->Write("mult_pT_All_Reconstructed");
//   TH2D * multptRecPP = hnTracksInBinsRecPP->Projection(4,0);
//   multptRecPP->Sumw2();
//   multptRecPP->Write("mult_pT_All_Reconstructed_from_Physical_Primary");  
  TH2D * multptMC = hnTracksInBinsMC->Projection(4,0);
  multptMC->Sumw2();
  multptMC->Write("mult_pT_All_Produced");    
  TH2D * multptEffRec = dynamic_cast<TH2D*>(multptRec->Clone("Efficiency_mult_pT_Rec"));
  multptEffRec->Divide(multptMC);
  multptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multptEffRec->Write();
//   TH2D * multptEffRecPP = dynamic_cast<TH2D*>(multptRecPP->Clone("Efficiency_mult_pT_Rec_pp"));
//   multptEffRecPP->Divide(multptMC);
//   multptEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   multptEffRecPP->Write();
  delete multptRec; /*delete multptRecPP;*/ delete multptMC;delete multptEffRec;/*delete multptEffRecPP;*/

  TH2D * vzphiRec = hnTracksInBins->Projection(2,1);
  vzphiRec->Sumw2();
  vzphiRec->Write("VZ_phi_All_Reconstructed");
//   TH2D * vzphiRecPP = hnTracksInBinsRecPP->Projection(2,1);
//   vzphiRecPP->Sumw2();
//   vzphiRecPP->Write("VZ_phi_All_Reconstructed_from_Physical_Primary");  
  TH2D * vzphiMC = hnTracksInBinsMC->Projection(2,1);
  vzphiMC->Sumw2();
  vzphiMC->Write("VZ_phi_All_Produced");    
  TH2D * vzphiEffRec = dynamic_cast<TH2D*>(vzphiRec->Clone("Efficiency_vz_phi_Rec"));
  vzphiEffRec->Divide(vzphiMC);
  vzphiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzphiEffRec->Write();
//   TH2D * vzphiEffRecPP = dynamic_cast<TH2D*>(vzphiRecPP->Clone("Efficiency_vz_phi_Rec_pp"));
//   vzphiEffRecPP->Divide(vzphiMC);
//   vzphiEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   vzphiEffRecPP->Write();  
  delete vzphiRec; /*delete vzphiRecPP;*/ delete vzphiMC;delete vzphiEffRec;/*delete vzphiEffRecPP;*/
  
  TH2D * vzetaRec = hnTracksInBins->Projection(3,1);
  vzetaRec->Sumw2();
  vzetaRec->Write("VZ_Eta_All_Reconstructed");
//   TH2D * vzetaRecPP = hnTracksInBinsRecPP->Projection(3,1);
//   vzetaRecPP->Sumw2();
//   vzetaRecPP->Write("VZ_Eta_All_Reconstructed_from_Physical_Primary");  
  TH2D * vzetaMC = hnTracksInBinsMC->Projection(3,1);
  vzetaMC->Sumw2();
  vzetaMC->Write("VZ_Eta_All_Produced");    
  TH2D * vzetaEffRec = dynamic_cast<TH2D*>(vzetaRec->Clone("Efficiency_vz_eta_Rec"));
  vzetaEffRec->Divide(vzetaMC);
  vzetaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzetaEffRec->Write();
//   TH2D * vzetaEffRecPP = dynamic_cast<TH2D*>(vzetaRecPP->Clone("Efficiency_vz_eta_Rec_pp"));
//   vzetaEffRecPP->Divide(vzetaMC);
//   vzetaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   vzetaEffRecPP->Write();
  delete vzetaRec; /*delete vzetaRecPP;*/ delete vzetaMC;delete vzetaEffRec;/*delete vzetaEffRecPP;*/
  
  TH2D * vzptRec = hnTracksInBins->Projection(4,1);
  vzptRec->Sumw2();
  vzptRec->Write("VZ_pT_All_Reconstructed");
/*  TH2D * vzptRecPP = hnTracksInBinsRecPP->Projection(4,1);
  vzptRecPP->Sumw2();
  vzptRecPP->Write("VZ_pT_All_Reconstructed_from_Physical_Primary"); */ 
  TH2D * vzptMC = hnTracksInBinsMC->Projection(4,1);
  vzptMC->Sumw2();
  vzptMC->Write("VZ_pT_All_Produced");    
  TH2D * vzptEffRec = dynamic_cast<TH2D*>(vzptRec->Clone("Efficiency_vz_pT_Rec"));
  vzptEffRec->Divide(vzptMC);
  vzptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzptEffRec->Write();
//   TH2D * vzptEffRecPP = dynamic_cast<TH2D*>(vzptRecPP->Clone("Efficiency_vz_pT_Rec_pp"));
//   vzptEffRecPP->Divide(vzptMC);
//   vzptEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   vzptEffRecPP->Write();
  delete vzptRec;/* delete vzptRecPP;*/ delete vzptMC;delete vzptEffRec;/*delete vzptEffRecPP;*/

  TH2D * phietaRec = hnTracksInBins->Projection(3,2);
  phietaRec->Sumw2();
  phietaRec->Write("phi_Eta_All_Reconstructed");
//   TH2D * phietaRecPP = hnTracksInBinsRecPP->Projection(3,2);
//   phietaRecPP->Sumw2();
//   phietaRecPP->Write("phi_Eta_All_Reconstructed_from_Physical_Primary");  
  TH2D * phietaMC = hnTracksInBinsMC->Projection(3,2);
  phietaMC->Sumw2();
  phietaMC->Write("phi_Eta_All_Produced");    
  TH2D * phietaEffRec = dynamic_cast<TH2D*>(phietaRec->Clone("Efficiency_phi_eta_Rec"));
  phietaEffRec->Divide(phietaMC);
  phietaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phietaEffRec->Write();
//   TH2D * phietaEffRecPP = dynamic_cast<TH2D*>(phietaRecPP->Clone("Efficiency_phi_eta_Rec_pp"));
//   phietaEffRecPP->Divide(phietaMC);
//   phietaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   phietaEffRecPP->Write();
  delete phietaRec; /*delete phietaRecPP; */delete phietaMC;delete phietaEffRec;/*delete phietaEffRecPP;*/

  TH2D * phiptRec = hnTracksInBins->Projection(4,2);
  phiptRec->Sumw2();
  phiptRec->Write("phi_pT_All_Reconstructed");
//   TH2D * phiptRecPP = hnTracksInBinsRecPP->Projection(4,2);
//   phiptRecPP->Sumw2();
//   phiptRecPP->Write("phi_pT_All_Reconstructed_from_Physical_Primary");  
  TH2D * phiptMC = hnTracksInBinsMC->Projection(4,2);
  phiptMC->Sumw2();
  phiptMC->Write("phi_pT_All_Produced");    
  TH2D * phiptEffRec = dynamic_cast<TH2D*>(phiptRec->Clone("Efficiency_phi_pT_Rec"));
  phiptEffRec->Divide(phiptMC);
  phiptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiptEffRec->Write();
//   TH2D * phiptEffRecPP = dynamic_cast<TH2D*>(phiptRecPP->Clone("Efficiency_phi_pT_Rec_pp"));
//   phiptEffRecPP->Divide(phiptMC);
//   phiptEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
//   phiptEffRecPP->Write();
  delete phiptRec; /*delete phiptRecPP; */delete phiptMC;delete phiptEffRec;/*delete phiptEffRecPP;*/
  
  TH2D * etaptRec = hnTracksInBins->Projection(4,3);
  etaptRec->Sumw2();
  etaptRec->Write("eta_pT_All_Reconstructed");
  TH2D * etaptMC = hnTracksInBinsMC->Projection(4,3);
  etaptMC->Sumw2();
  etaptMC->Write("eta_pT_All_Produced");    
  TH2D * etaptEffRec = dynamic_cast<TH2D*>(etaptRec->Clone("Efficiency_eta_pT_Rec"));
  etaptEffRec->Divide(etaptMC);
  etaptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaptEffRec->Write();  
  delete etaptRec; delete etaptMC;delete etaptEffRec;

  //make more meaningful axes:
  outfile->cd();
  //mult: keep for now, unsure how it actually should look.
  Double_t multaxisArray[9] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,95.0};
  TArrayD * multaxisA = new TArrayD(9,multaxisArray);
  //vz: 4 cm bins from -6 to 6,6-8,8-10 and the same on the other side:
  Double_t  vzaxisArray[8] = {-10.0,-7.5,-5.0,-2.0,2.0,5.0,7.5,10.0};
  TArrayD * vzaxisA = new TArrayD(8,vzaxisArray);
  Double_t etamin = -0.9;
  Double_t etamax = 0.9;
  Double_t deta = (etamax - etamin)/63;
  TArrayD * etaaxisAD = new TArrayD(63+1);
  for(int i = 0; i<=63;i++){
    etaaxisAD->AddAt(etamin+i*deta,i);
  }
  
  //phi need all the bins I can get... 
//   const Double_t * phibins = hnTracksInBins->GetAxis(2)->GetXbins()->GetArray();
  //eta: same
//   const Double_t * etabins = hnTracksInBins->GetAxis(3)->GetXbins()->GetArray();
  //pt: 0.1GeV/c up to 3 GeV/c = 25 bins, 0.5GeV/c up to 5 GeV/c = 4 bins, then 1 5GeV/c bin and one 6GeV/c bin: total 30 bins
  Double_t  pTaxisArray[16] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  TArrayD * pTaxisA = new TArrayD(16,pTaxisArray);
  
  TH3D * histrec = Rebin(hnTracksInBins,"hist3drec",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histrec->Write();
  TH3D * histMC = Rebin(hnTracksInBinsMC,"hist3dMC",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histMC->Write();
  histMC->Divide(histrec);
  histMC->Write("Eff3d");
  TFile* outfile2;
  if(global){
    outfile2 = TFile::Open(filename.ReplaceAll("eff.root","LHC10hWeight.root"),"RECREATE");
    outfile2->cd();
    histMC->Write("hnWeight");
  }
  outfile->cd();
  TH1D * err = new TH1D("err","Errors in the 3d histogram",100,0.0,20.0);
  for(int x=1;x<=histrec->GetNbinsX();x++){
    for(int y=1;y<=histrec->GetNbinsY();y++){
      for(int z=1;z<=histrec->GetNbinsZ();z++){
	if(histrec->GetBinContent(x,y,z)>1.0e-10)err->Fill(100*histrec->GetBinError(x,y,z)/histrec->GetBinContent(x,y,z));
      }
    }
  }
  err->Write();
  
//   Double_t etamin = Etaaxis->GetBinLowEdge(1);
//   Double_t etamax = Etaaxis->GetBinUpEdge(Etaaxis->GetNbins());
//   Double_t deta = (etamax - etamin)/Etaaxis->GetNbins();
//   TArrayD * etaaxisAD = new TArrayD(Etaaxis->GetNbins()/3.0+1);
//   for(int i = 0; i<=Etaaxis->GetNbins()/3.0;i++){
//     etaaxisAD->AddAt(etamin+i*3.0*deta,i);
//   }
//   const TArrayD * etabins = hnTracksInBins->GetAxis(3)->GetXbins();
//   cout << etabins->GetSize()<<endl;
  TH2D * histrechpT = Rebin(hnTracksInBins,"hist3drechpT",3,etaaxisAD,1,vzaxisA,0,2,4,2.0);
  histrechpT->Write();
  TH2D * histhpTMC = Rebin(hnTracksInBinsMC,"hist3dMChpT",3,etaaxisAD,1,vzaxisA,0,2,4,2.0);
  histhpTMC->Write();  
  histhpTMC->Divide(histrechpT);
  Double_t nbins = 0.0;
  Double_t content = 0.0;
  for(int x = 1;x<=histhpTMC->GetXaxis()->GetNbins();x++){
     for(int y = 1;y<=histhpTMC->GetYaxis()->GetNbins();y++){
	if(histhpTMC->GetBinContent(x,y)>1.0E-10){
	  nbins +=1.0;
	  content+=histhpTMC->GetBinContent(x,y);       
	}
      }
    }
//   cout << content/nbins<<endl;  
  outfile->cd();

  histhpTMC->Write("Eff2dhpT");
  TF1* ptfunc;
  if(global){
    outfile2->cd();
    histhpTMC->Write("hnWeight_highpt");
    ptfunc = new TF1("pT_function","[0]+[1]*x+[2]*x*x+[3]*exp(-x)",2,100);
    pTWeightRec->Scale(nbins/content);
    pTWeightRec->Fit(ptfunc,"","",2,100);
  }
  outfile->cd();
  pTWeightRec->Write();
  if(global){
    outfile2->cd();
    ptfunc->Write();
    outfile2->Close();
  }
  delete hnTracksInBins; delete hnTracksInBinsMC; /*delete hnTracksInBinsRecPP;*/delete hnCentVsVertex;    
  delete pTWeightRec;
  delete histrec; delete histMC;delete err;delete histrechpT;delete histhpTMC;
  if(global)delete ptfunc;
//   delete etaaxisAD; 
  delete multaxisA; delete vzaxisA; delete pTaxisA;
  outfile->Close();
  
  /*
  
  
  
  

  




//   cout << "beefopre"<<endl;
//   THnF * RebinRec = Rebin(7,multaxisArray,7,vzaxisArray,1,hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins()/3.0,hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBins,"hnRebinRec");
//   RebinRec->Write();
//   THnF * RebinMC = Rebin(7,multaxisArray,7,vzaxisArray,1,hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins()/3.0,hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBinsMC,"hnRebinMC");
//   RebinMC->Write();  
//   
//   THnF * weights = dynamic_cast<THnF*>(RebinRec->Clone("hnWeights"));
//   THnF * weightsnoerr = dynamic_cast<THnF*>(RebinRec->Clone("hnWeight"));
//   weights->Clear();
//   weightsnoerr->Clear();
//   cout << "beefopreas"<<endl;
// 
//   int nerr = 0;
//   int nnoerr = 0;
//   for(int mult = 1;   mult<=weights->GetAxis(0)->GetNbins();mult++){
//     for(int vz = 1;   vz<=weights->GetAxis(1)->GetNbins();vz++){
// //       for(int phi = 1;phi<=weights->GetAxis(2)->GetNbins();phi++){
// 	for(int eta=1;eta<=weights->GetAxis(2)->GetNbins();eta++){
// 	  for(int pT=1;pT<=weights->GetAxis(3)->GetNbins();pT++){
// 	    Int_t index[5] = {mult,vz,eta,pT};
// 	    Double_t BinContentRec = RebinRec->GetBinContent(index) ;
// 	    Double_t BinContentMC = RebinMC->GetBinContent(index) ;
// 	    Double_t BinContent = 0.0;
// 	    Double_t BinError = 0.0;
// 	    if(!BinContentRec<0.5){
// 	      BinContent = BinContentMC/BinContentRec;
// 	      BinError = BinContent*BinContent*(RebinRec->GetBinError2(RebinRec->GetBin(index))/(BinContentRec*BinContentRec) + RebinMC->GetBinError2(RebinMC->GetBin(index))/(BinContentMC*BinContentMC));
// // 	      cout << 1 <<" "<< RebinRec->GetBinError2(RebinRec->GetBin(index)) <<" "<<BinContentRec<<endl;
// // 	      cout << 2 <<" "<< RebinMC->GetBinError2(RebinMC->GetBin(index)) <<" "<<BinContentMC<<endl;
// // 	      cout << 3 <<" "<< BinError <<" "<<BinContent<<endl;
// 	    }
// 	    weights->SetBinContent(index,BinContent);
// 	    weightsnoerr->SetBinContent(index,BinContent);
// 	    if(BinError >1.0E-10)weights->SetBinError2(weights->GetBin(index),BinError);
// 	    
// 	    if(TMath::Sqrt(BinError)*100>BinContent) {nerr+=1;}
// 	    else nnoerr+=1;
// 	    
// 	  }
// 	}
// //       }
//     }
//   }
//   cout << nerr<<" "<<nnoerr<<endl;
//   weights->Write("hnWeights");
//   Errors(weights, "errorsweights",true)->Write();
  /*
  //Canvases to compare efficiencies in different binnings:
  TCanvas * plotcanvas = new TCanvas("plotcanvas");
  multEffRec->Draw("E");
  TH1D * multRecrb = RebinRec->Projection(0);
  TH1D * multMCrb = RebinMC->Projection(0);
  multRecrb->Sumw2();
  multMCrb->Sumw2();
  TH1D * multEffRecrb = dynamic_cast<TH1D*>(multRecrb->Clone("Efficiency_Mult_Rec_reb"));
  multEffRecrb->Divide(multMCrb);
  multEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRecrb->Write();
  multEffRecrb->Draw("sameE");
  plotcanvas->Write("MultEffComparison");
  plotcanvas->Clear();
  
  vzEffRec->Draw("E");
  TH1D * vzRecrb = RebinRec->Projection(1);
  TH1D * vzMCrb = RebinMC->Projection(1);
  vzRecrb->Sumw2();
  vzMCrb->Sumw2();
  TH1D * vzEffRecrb = dynamic_cast<TH1D*>(vzRecrb->Clone("Efficiency_VZ_Rec_rb"));
  vzEffRecrb->Divide(vzMCrb);
  vzEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRecrb->Write();
  vzEffRecrb->Draw("sameE");
  plotcanvas->Write("VZEffComparison");
  plotcanvas->Clear();

  phiEffRec->Draw("E");
  TH1D * phiRecrb = RebinRec->Projection(2);
  TH1D * phiMCrb = RebinMC->Projection(2);
  phiRecrb->Sumw2();
  phiMCrb->Sumw2();
  TH1D * phiEffRecrb = dynamic_cast<TH1D*>(phiRecrb->Clone("Efficiency_phi_Rec_rb"));
  phiEffRecrb->Divide(phiMCrb);
  phiEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiEffRecrb->Write();
  phiEffRecrb->Draw("sameE");
  plotcanvas->Write("phiEffComparison");
  plotcanvas->Clear();
  
  etaEffRec->Draw("E");
  TH1D * etaRecrb = RebinRec->Projection(3);
  TH1D * etaMCrb = RebinMC->Projection(3);
  etaRecrb->Sumw2();
  etaMCrb->Sumw2();
  TH1D * etaEffRecrb = dynamic_cast<TH1D*>(etaRecrb->Clone("Efficiency_eta_Rec_rb"));
  etaEffRecrb->Divide(etaMCrb);
  etaEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaEffRecrb->Write();
  etaEffRecrb->Draw("sameE");
  plotcanvas->Write("etaEffComparison");
  plotcanvas->Clear();

  pTEffRec->Draw("E");
  TH1D * pTRecrb = RebinRec->Projection(4);
  TH1D * pTMCrb = RebinMC->Projection(4);
  pTRecrb->Sumw2();
  pTMCrb->Sumw2();
  TH1D * pTEffRecrb = dynamic_cast<TH1D*>(pTRecrb->Clone("Efficiency_pT_Rec_rb"));
  pTEffRecrb->Divide(pTMCrb);
  pTEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  pTEffRecrb->Write();
  pTEffRecrb->Draw("sameE");
  plotcanvas->Write("pTEffComparison");
  plotcanvas->Clear();
  
  //2d projections in correlated variables: 
  TH2D * phiptweight = RebinRec->Projection(4,2);
  phiptweight->Divide(RebinMC->Projection(4,2));
  phiptweight->Write("Weights_Phi_pT");
//   phiptweight->Scale(1.0/phiptweight->Integral());
//   phiptweight->Write("Weights_Phi_pT_norm");
  TH2D * vzetaweight = RebinRec->Projection(3,1);
  vzetaweight->Divide(RebinMC->Projection(3,1));
  vzetaweight->Write("Weights_Eta_VZ");
  vzetaweight->Scale(1.0/vzetaweight->GetMaximum());
  vzetaweight->Write("Weights_Eta_VZ_norm");
  //Centnorm:
  multEffRecrb->Scale(1.0/multEffRecrb->GetMaximum());
  multEffRecrb->Write("Weights_Cent_norm");

  
 // TH2D * multptweight = weights->Projection(4,0);
 // multptweight->Write("weights_pT_All_Reconstructed");
  
  TH1D* differences = new TH1D("differencesinbins","Differences between N-dim and factorized",100,0.0,1.0);
  for(int mult = 1;mult<=weightsnoerr->GetAxis(0)->GetNbins();mult++){
    for(int vz = 1;vz<=weightsnoerr->GetAxis(1)->GetNbins();vz++){
      for(int phi =1;phi<=weightsnoerr->GetAxis(2)->GetNbins();phi++){
	for(int eta = 1;eta<=weightsnoerr->GetAxis(3)->GetNbins();eta++){
	  for(int pt = 1;pt<=weightsnoerr->GetAxis(4)->GetNbins();pt++){
	    Int_t index[5] = {mult,vz,phi,eta,pt};
	    differences->Fill((TMath::Abs(weightsnoerr->GetBinContent(index)- phiptweight->GetBinContent(phi,pt)*vzetaweight->GetBinContent(vz,eta)*multEffRecrb->GetBinContent(mult))/TMath::Abs(weightsnoerr->GetBinContent(index))));
	  }
	}
      }
    }	
  }
  differences->Write();
  
  outfile->Close();
//   TFile* outfile2 = TFile::Open("LHC10hWeight.root","RECREATE");
//   weightsnoerr->Write("hnWeight");
//   
//   
//   
//   outfile2->Close();
  
//   infile->Close();
  
//   delete hnTracksInBins; delete hnTracksInBinsRecPP; delete hnTracksInBinsMC;
// //   delete multEffRec;delete multEffRecPP;
//   
*/
}

