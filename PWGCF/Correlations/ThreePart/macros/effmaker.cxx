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
double glimit = 4.5;
double gMaxWeight = 1000.0;
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

Double_t myfunc(Double_t* x,Double_t* par){
  double limit = glimit;
  if(x[0]<limit){
    return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
  }
  else return par[0]+limit*par[1]+limit*limit*par[2]+limit*limit*limit*par[3];
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
  TList * dir = NULL;
  TDirectory * outdir = infile->GetDirectory("ThreePartTrackEfficiencies_0_100");
  if(outdir){
    dir = dynamic_cast<TList*>(outdir->Get("ThreePartTrackEfficiencies_0_100_0_0Coutput1"));
  }
  else{
    outdir = infile->GetDirectory("AddTaskThreePartTrackEfficiencies_0_100");
    if(outdir){
      dir = dynamic_cast<TList*>(outdir->Get("AddTaskThreePartTrackEfficiencies_0_100Coutput1"));
    }
    else return;
  }
  TFile* outfile = TFile::Open("eff.root","RECREATE");
  outfile->cd();
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBins = dynamic_cast<THnF*>(dir->At(i)->Clone("hnTracksinBins_eff"));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMC = dynamic_cast<THnF*>(dir->At(i)->Clone("hnTracksinBinsMC_eff"));
    if(TString(dir->At(i)->GetName()).CompareTo("centVsZVertex")==0) hnCentVsVertex = dynamic_cast<TH2D*>(dir->At(i)->Clone("centVsZVertex_eff"));
  }
  infile->Close();
  
  TCanvas * canvas = new TCanvas("canvas"); 
  canvas->cd();
  
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
  multEffRec->SetMinimum(0.00);
  multEffRec->GetYaxis()->SetTitle("#varepsilon []");
  multEffRec->SetStats(false);
  multEffRec->Draw();
  canvas->SaveAs("Efficiency_Mult.eps");
  canvas->Clear();
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
  vzEffRec->SetMinimum(0.00);
  vzEffRec->GetYaxis()->SetTitle("#varepsilon []");
  vzEffRec->SetStats(false);
  vzEffRec->Draw();
  canvas->SaveAs("Efficiency_VZ.eps");
  canvas->Clear();
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
  etaEffRec->SetMinimum(0.00);  
  etaEffRec->GetYaxis()->SetTitle("#varepsilon []");
  etaEffRec->SetStats(false);
  etaEffRec->Draw();
  canvas->SaveAs("Efficiency_eta.eps");
  canvas->Clear();
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
  pTEffRec->GetXaxis()->SetRangeUser(0.0,50.0);
  pTEffRec->SetMinimum(0.00);
  pTEffRec->GetYaxis()->SetTitle("#varepsilon []");
  pTEffRec->SetStats(false);
  pTEffRec->Draw();
  canvas->SaveAs("Efficiency_pT.eps");
  canvas->Clear();
  TH1D * pTWeightRec = dynamic_cast<TH1D*>(pTMC->Clone("Weight_pT_Rec"));
  
  double averagept = 0;
  double ntracks = 0;
  for(int i = pTRec->FindBin(2.0);i<pTRec->GetNbinsX();i++){
    averagept += pTRec->GetBinCenter(i)*pTRec->GetBinContent(i);
    ntracks += pTRec->GetBinContent(i);
  }
  averagept = averagept/ntracks;
  pTWeightRec->Divide(pTRec);
  double weightataverage = pTWeightRec->GetBinContent(pTWeightRec->FindBin(averagept));
  
  
  pTWeightRec->SetTitle(" produced MC particles divided by Reconstructed tracks");
  delete pTRec; 
  delete pTMC;delete pTEffRec;

  
  //2D projections
  TDirectory * dim2 = outfile->mkdir("Projections_2d");
  dim2->cd();
  TCanvas * can2d = new TCanvas("canvas2");
  
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
//   can2d->SaveAs("Efficiency_Mult_pT.eps");
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
  
//   can2d->SaveAs("Efficiency_Mult_pT.eps");

  
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

//   outfile2->Close();
  outfile->cd();
  TH1D * err = new TH1D("err","Errors in the 3d histogram",100,0.0,20.0);
  for(int x=1;x<=histMC->GetNbinsX();x++){
    for(int y=1;y<=histMC->GetNbinsY();y++){
      for(int z=1;z<=histMC->GetNbinsZ();z++){
	if(histMC->GetBinContent(x,y,z)>1.0e-10){
	  if(histMC->GetBinContent(x,y,z)>gMaxWeight){
	    histMC->SetBinContent(x,y,z,0.0);
	    histMC->SetBinError(x,y,z,0.0);
	  }
	  else err->Fill(100*histMC->GetBinError(x,y,z)/histMC->GetBinContent(x,y,z));
	}
      }
    }
  }
  histMC->Write("Eff3d");
  TFile* outfile2 = TFile::Open("LHC11aWeight.root","RECREATE");
  outfile2->cd();
  histMC->Write("hnWeight");  
  outfile->cd();
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
	  if(histhpTMC->GetBinContent(x,y)>gMaxWeight){
	    histhpTMC->SetBinContent(x,y,0.0);
	    histhpTMC->SetBinError(x,y,0.0);
	  }
	  else{
	    nbins +=1.0;
	    content+=histhpTMC->GetBinContent(x,y);
	  }
	}
      }
    }
  histhpTMC->Write("Eff2dhpT");
  TF1* ptfunc;
  outfile2->cd();
  histhpTMC->Write("hnWeight_highpt");
  can2d->cd();
  histhpTMC->SetStats(kFALSE);
  histhpTMC->SetTitle("Weight (N_{MC}/N_{reconstructed}).");
  histhpTMC->GetZaxis()->SetLabelSize(0.02);
  histhpTMC->Draw("colz");
  can2d->SaveAs("Weight_eta_vz.eps");
  ptfunc = new TF1("pT_function",myfunc,2,50,4);
  ptfunc->SetNpx(10000);
  glimit = 6.0;
  pTWeightRec->Scale(1.0/weightataverage);
  pTWeightRec->Fit(ptfunc,"","",2,50);
  outfile->cd();
  pTWeightRec->Write();
  canvas->cd();
  pTWeightRec->Draw();
  pTWeightRec->GetXaxis()->SetRangeUser(0.0,50.0);
  pTWeightRec->SetMaximum(1.2);
  pTWeightRec->SetMinimum(0.6);
  pTWeightRec->SetTitle("Weight in pT scaled to 1 at average pT.");
  canvas->Update();
  canvas->SaveAs("Weight_pT.eps");
  canvas->Clear();
  
  outfile2->cd();
  ptfunc->Write();
  outfile2->Close();  
  delete hnTracksInBins; delete hnTracksInBinsMC; delete hnCentVsVertex;    
  delete pTWeightRec;
  delete histrec; delete histMC;delete err;delete histrechpT;delete histhpTMC;
  delete ptfunc;
  delete multaxisA; delete vzaxisA; delete pTaxisA;  
  delete canvas;delete can2d;

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
  THnF * hnTracksInBinsMCtmp;  THnF * hnTracksInBinsMC;
  TH2D * hnCentVsVertextmp;  TH2D * hnCentVsVertex;
  TList * dir;
  if(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")){
   dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")->Get("ThreePartTrackEfficienciesPbPb_0_16_0_0Coutput1"));
  }
  else if(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_100")){
   dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_100")->Get("ThreePartTrackEfficienciesPbPb_0_100_0_0Coutput1"));    
  }
  else{ 
    infile->Close();
    return;
  }
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBinstmp = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMCtmp = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("centVsZVertex")==0) hnCentVsVertextmp = (dynamic_cast<TH2D*>(dir->At(i)));
  }

  TFile* outfile = TFile::Open(filename.ReplaceAll("AnalysisResults.root","eff.root"),"RECREATE");
  cout << filename.Data()<<endl;
  outfile->cd();
  if(filename.Contains("LHC12a17g")){
    hnTracksInBins = SumCent(0.0,10.0,hnTracksInBinstmp,"hnTracksReconstruced");
    hnTracksInBinsMC = SumCent(0.0,10.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(0.0,10.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  if(filename.Contains("LHC12a17h")){
    hnTracksInBins = SumCent(10.0,50.0,hnTracksInBinstmp,"hnTracksReconstruced");
    hnTracksInBinsMC = SumCent(10.0,50.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(10.0,50.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  if(filename.Contains("LHC12a17i")){
    hnTracksInBins = SumCent(50.0,100.0,hnTracksInBinstmp,"hnTracksReconstruced");
    hnTracksInBinsMC = SumCent(50.0,100.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(50.0,100.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  if(filename.Contains("LHC11a10a")||filename.CompareTo("eff.root")==0){
    hnTracksInBins = SumCent(0.0,100.0,hnTracksInBinstmp,"hnTracksReconstruced");
    hnTracksInBinsMC = SumCent(0.0,100.0,hnTracksInBinsMCtmp,"hnTracksProduced");
    hnCentVsVertex = SumCent(0.0,100.0,hnCentVsVertextmp,"hnCentralityvsVertex");
  }
  infile->Close();
  dir->Delete();
  delete dir;
  TCanvas * canvas = new TCanvas("canvas"); 
  canvas->cd();
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
  multEffRec->SetMinimum(0.0);
  multEffRec->GetYaxis()->SetTitle("#varepsilon []");
  multEffRec->SetStats(false);
  multEffRec->Draw();
  canvas->SaveAs("Efficiency_Mult.eps");
  canvas->Clear();
  delete multRec; 
  delete multMC;delete multEffRec;


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
  vzEffRec->SetMinimum(0.00);
  vzEffRec->GetYaxis()->SetTitle("#varepsilon []");
  vzEffRec->SetStats(false);
  vzEffRec->Draw();
  canvas->SaveAs("Efficiency_VZ.eps");
  canvas->Clear();
  delete vzRec; 
  delete vzMC;delete vzEffRec;
  
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
  delete phiRec; 
  delete phiMC;delete phiEffRec;
  
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
  etaEffRec->SetMinimum(0.00);
  etaEffRec->GetYaxis()->SetTitle("#varepsilon []");
  etaEffRec->SetStats(false);
  etaEffRec->Draw();
  canvas->SaveAs("Efficiency_eta.eps");
  canvas->Clear();
  delete etaRec; 
  delete etaMC;delete etaEffRec;

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
  pTEffRec->GetXaxis()->SetRangeUser(0.0,50.0);
  pTEffRec->SetMinimum(0.00);
  pTEffRec->GetYaxis()->SetTitle("#varepsilon []");
  pTEffRec->SetStats(false);
  pTEffRec->Draw();
  canvas->SaveAs("Efficiency_pT.eps");
  canvas->Clear();
  TH1D * pTWeightRec = dynamic_cast<TH1D*>(pTMC->Clone("Weight_pT_Rec"));
  double averagept = 0;
  double ntracks = 0;
  for(int i = pTRec->FindBin(2.0);i<pTRec->GetNbinsX();i++){
    averagept += pTRec->GetBinCenter(i)*pTRec->GetBinContent(i);
    ntracks += pTRec->GetBinContent(i);
  }
  averagept = averagept/ntracks;
  pTWeightRec->Divide(pTRec);
  double weightataverage = pTWeightRec->GetBinContent(pTWeightRec->FindBin(averagept));
  pTWeightRec->SetTitle(" produced MC particles divided by Reconstructed tracks");
  delete pTRec; 
  delete pTMC;delete pTEffRec;
  

  //2D projections
  TDirectory * dim2 = outfile->mkdir("Projections_2d");
  dim2->cd();
  
  TH2D * multvzrec = hnTracksInBins->Projection(1,0);
  multvzrec->Sumw2();
  multvzrec->Write("mult_VZ_All_Reconstructed");
  TH2D * multvzMC = hnTracksInBinsMC->Projection(1,0);
  multvzMC->Sumw2();
  multvzMC->Write("mult_VZ_All_Produced");    
  TH2D * multvzeffrec = dynamic_cast<TH2D*>(multvzrec->Clone("Efficiency_mult_vz_Rec"));
  multvzeffrec->Divide(multvzMC);
  multvzeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multvzeffrec->Write();
  delete multvzrec;  delete multvzMC;delete multvzeffrec;

  TH2D * multphirec = hnTracksInBins->Projection(2,0);
  multphirec->Sumw2();
  multphirec->Write("mult_phi_All_Reconstructed");
  TH2D * multphiMC = hnTracksInBinsMC->Projection(2,0);
  multphiMC->Sumw2();
  multphiMC->Write("mult_phi_All_Produced");    
  TH2D * multphieffrec = dynamic_cast<TH2D*>(multphirec->Clone("Efficiency_mult_phi_Rec"));
  multphieffrec->Divide(multphiMC);
  multphieffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multphieffrec->Write();
  delete multphirec;  delete multphiMC;delete multphieffrec;

  TH2D * multetarec = hnTracksInBins->Projection(3,0);
  multetarec->Sumw2();
  multetarec->Write("mult_eta_All_Reconstructed");
  TH2D * multetaMC = hnTracksInBinsMC->Projection(3,0);
  multetaMC->Sumw2();
  multetaMC->Write("mult_eta_All_Produced");    
  TH2D * multetaeffrec = dynamic_cast<TH2D*>(multetarec->Clone("Efficiency_mult_eta_Rec"));
  multetaeffrec->Divide(multetaMC);
  multetaeffrec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multetaeffrec->Write();
  delete multetarec;  delete multetaMC;delete multetaeffrec;
  
  TH2D * multptRec = hnTracksInBins->Projection(4,0);
  multptRec->Sumw2();
  multptRec->Write("mult_pT_All_Reconstructed");
  TH2D * multptMC = hnTracksInBinsMC->Projection(4,0);
  multptMC->Sumw2();
  multptMC->Write("mult_pT_All_Produced");    
  TH2D * multptEffRec = dynamic_cast<TH2D*>(multptRec->Clone("Efficiency_mult_pT_Rec"));
  multptEffRec->Divide(multptMC);
  multptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multptEffRec->Write();
  delete multptRec;  delete multptMC;delete multptEffRec;

  TH2D * vzphiRec = hnTracksInBins->Projection(2,1);
  vzphiRec->Sumw2();
  vzphiRec->Write("VZ_phi_All_Reconstructed");
  TH2D * vzphiMC = hnTracksInBinsMC->Projection(2,1);
  vzphiMC->Sumw2();
  vzphiMC->Write("VZ_phi_All_Produced");    
  TH2D * vzphiEffRec = dynamic_cast<TH2D*>(vzphiRec->Clone("Efficiency_vz_phi_Rec"));
  vzphiEffRec->Divide(vzphiMC);
  vzphiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzphiEffRec->Write();
  delete vzphiRec; delete vzphiMC;delete vzphiEffRec;
  
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
  delete vzetaRec;  delete vzetaMC;delete vzetaEffRec;
  
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
  delete vzptRec; delete vzptMC;delete vzptEffRec;

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
  delete phietaRec; delete phietaMC;delete phietaEffRec;

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

  //pt: 0.1GeV/c up to 3 GeV/c = 25 bins, 0.5GeV/c up to 5 GeV/c = 4 bins, then 1 5GeV/c bin and one 6GeV/c bin: total 30 bins
  Double_t  pTaxisArray[16] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  TArrayD * pTaxisA = new TArrayD(16,pTaxisArray);
  
  TH3D * histrec = Rebin(hnTracksInBins,"hist3drec",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histrec->Write();
  TH3D * histMC = Rebin(hnTracksInBinsMC,"hist3dMC",3,etaaxisAD,1,vzaxisA,4,pTaxisA,2,0);
  histMC->Write();
  histMC->Divide(histrec);
  outfile->cd();
  TH1D * err = new TH1D("err","Errors in the 3d histogram",100,0.0,20.0);
  for(int x=1;x<=histMC->GetNbinsX();x++){
    for(int y=1;y<=histMC->GetNbinsY();y++){
      for(int z=1;z<=histMC->GetNbinsZ();z++){
	if(histMC->GetBinContent(x,y,z)>1.0e-10){
	  if(histMC->GetBinContent(x,y,z)>gMaxWeight){
	    histMC->SetBinContent(x,y,z,0.0);
	    histMC->SetBinError(x,y,z,0.0);
	  }
	  else err->Fill(100*histMC->GetBinError(x,y,z)/histMC->GetBinContent(x,y,z));
	}
      }
    }
  }
  err->Write();
  histMC->Write("Eff3d");

  TFile* outfile2;
  if(global){
    outfile2 = TFile::Open(filename.ReplaceAll("eff.root","LHC10hWeight.root"),"RECREATE");
    outfile2->cd();
    histMC->Write("hnWeight");
  }
  outfile->cd();
  
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
	  if(histhpTMC->GetBinContent(x,y)>gMaxWeight){
	    histhpTMC->SetBinContent(x,y,0.0);
	    histhpTMC->SetBinError(x,y,0.0);
	  }
	  else{
	    nbins +=1.0;
	    content+=histhpTMC->GetBinContent(x,y);
	  }
	}
      }
    }
  outfile->cd();

  histhpTMC->Write("Eff2dhpT");
  TF1* ptfunc;
  if(global){
    outfile2->cd();
    canvas->cd();
    histhpTMC->Write("hnWeight_highpt");
    histhpTMC->SetStats(kFALSE);
    histhpTMC->SetTitle("Weight (N_{MC}/N_{reconstructed}).");
    histhpTMC->GetZaxis()->SetLabelSize(0.02);
    histhpTMC->Draw("colz");
    canvas->SaveAs("Weight_eta_vz.eps");
    canvas->Clear();
    ptfunc = new TF1("pT_function",myfunc,2,50,4);
    ptfunc->SetNpx(10000);
    pTWeightRec->Scale(1.0/weightataverage);
    glimit = 4.5;
    pTWeightRec->Fit(ptfunc,"","",2,50);
    pTWeightRec->Draw();
    pTWeightRec->GetXaxis()->SetRangeUser(0.0,50.0);
    pTWeightRec->SetTitle("Weight in pT scaled to 1 at average pT.");
//     pTWeightRec->SetMaximum(1.2);
    pTWeightRec->SetMinimum(0.6);
    canvas->Update();
    canvas->SaveAs("Weight_pT.eps");
    canvas->Clear();
    
  }
  outfile->cd();
  pTWeightRec->Write();
  if(global){
    outfile2->cd();
    ptfunc->Write();
    outfile2->Close();
  }
  delete hnTracksInBins; delete hnTracksInBinsMC; delete hnCentVsVertex;    
  delete pTWeightRec;
  delete histrec; delete histMC;delete err;delete histrechpT;delete histhpTMC;
  if(global)delete ptfunc;
  delete etaaxisAD; 
  delete multaxisA; delete vzaxisA; delete pTaxisA;
  outfile->Close();
  

}

