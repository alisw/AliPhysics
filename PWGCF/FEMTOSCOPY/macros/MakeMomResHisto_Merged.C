#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"

using namespace std;

void MakeMomResHisto_Merged(){

  int MBmerged = 0;
  //

  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R11.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R10.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_lego_L0p68R9andR8.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R7.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R6.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_lego_L0p68R11_PID1p5_and_FB5.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R11_chi3p1.root","READ");
  TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R11_highKt.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17d_fix_genSignal_Rinv11.root","READ");
  //TFile *infile1 = new TFile("Results/PDC_HIJING_12a17a_myRun_L0p68R11_ncls.root","READ");
  //TFile *infile1 = new TFile("MyOutput.root","READ");
  //
  //TFile *infile2 = new TFile("Results/PDC_HIJING_12a17b_12a17e_noPadCut_lego.root","READ");
  //TFile *infile2 = new TFile("Results/PDC_HIJING_12a17b_newr3_lego.root","READ");
  //TFile *infile2 = new TFile("Results/PDC_HIJING_12a17b_GRS_lego.root","READ");
  //TFile *infile2 = new TFile("Results/PDC_HIJING_12a17b_myRun_L0p4R8.root","READ");

  //TDirectoryFile *tdir1 = (TDirectoryFile*)infile1->Get("PWGCF.outputChaoticityAnalysis.root");
  //TList *list1=(TList*)tdir1->Get("ChaoticityOutput_2");
  TList *list1=(TList*)infile1->Get("MyList");
  infile1->Close();
  //TDirectoryFile *tdir2 = (TDirectoryFile*)infile2->Get("PWGCF.outputChaoticityAnalysis.root");
  //TList *list2=(TList*)tdir2->Get("ChaoticityOutput");
  //TList *list2=(TList*)infile2->Get("MyList");
  //infile2->Close();

  TList *list;

  const int Mbins=2;

  TH2D *num_i_2[2];// 2-particle
  TH2D *den_i_2[2];// 2-particle
  TH2D *num_s_2[2];// 2-particle
  TH2D *den_s_2[2];// 2-particle
  //
  TH2D *merged_num_i_2[2];
  TH2D *merged_den_i_2[2];
  TH2D *merged_num_s_2[2];
  TH2D *merged_den_s_2[2];
  //

 
  TString *names[2][4];
  // 3d
  TString *names_terms[2][5][2];// charge-combo (ss or os), terms, ideal/smeared, Mbins
  TH3D *merged_terms_i[2][5];// charge-combo (ss or os), terms, Mbins
  TH3D *merged_terms_s[2][5];// charge-combo (ss or os), terms, Mbins
  TH1D *merged_terms1D_i[2][5];
  TH1D *merged_terms1D_s[2][5];
  //
  TString *names_termsK3[2][2];// charge-combo (ss or os), Sum/En, Mbins
  TH3D *merged_terms_SumK3[2];// charge-combo (ss or os), Mbins
  TH3D *merged_terms_EnK3[2];// charge-combo (ss or os), Mbins
  //
  TString *names_terms_FVP[5][2][2];// terms, ideal/smeared, Projection 1 or 2, Mbins
  TH3D *merged_terms_FVP_i[5][2];// terms, Projection 1 or 2, Mbins
  TH3D *merged_terms_FVP_s[5][2];// terms, Projection 1 or 2, Mbins
  //
  TString *names_terms_FVP_K[5][2][2];// terms, Sum/En, Projection 1 or 2, Mbins
  TH3D *merged_terms_FVP_SumK[5][2];// terms, Projection 1 or 2, Mbins
  TH3D *merged_terms_FVP_EnK[5][2];// terms, Projection 1 or 2, Mbins
  //
  TString *names_terms_FVPTPN[2][2];// ideal/smeared, Projection 1 or 2, Mbins
  TH3D *merged_terms_FVPTPN_i[2];// Projection 1 or 2, Mbins
  TH3D *merged_terms_FVPTPN_s[2];// Projection 1 or 2, Mbins

  
  for(int mb=0; mb<Mbins; mb++){// Mbin loop
    cout<<"Cent bin "<<mb<<endl;
    for(int ch=0; ch<2; ch++){// + or - loop
      
      int MixedCbin1=0, MixedCbin2=0, MixedCbin3=0;
      // 3d histos
      for(int term=0; term<5; term++){// term loop
	names_terms[0][term][0]=new TString("PairCut3_Charge1_"); *names_terms[0][term][0] += ch;
	names_terms[0][term][0]->Append("_Charge2_"); *names_terms[0][term][0] += ch;
	names_terms[0][term][0]->Append("_Charge3_"); *names_terms[0][term][0] += ch;
	names_terms[0][term][0]->Append("_SC_0_M_");
	*names_terms[0][term][0] += mb;
	names_terms[0][term][0]->Append("_ED_0_");
	TString *TPNnameBase = new TString(names_terms[0][term][0]->Data());
	names_terms[0][term][0]->Append("Term_"); *names_terms[0][term][0] += term+1; 
	names_terms[0][term][1]=new TString(names_terms[0][term][0]->Data());
	if(term==0) {
	  names_termsK3[0][0]=new TString(names_terms[0][term][0]->Data());
	  names_termsK3[0][1]=new TString(names_terms[0][term][0]->Data());
	  names_termsK3[0][0]->Append("_SumK3");
	  names_termsK3[0][1]->Append("_EnK3");
	}	
	//
	
	names_terms_FVP[term][0][0] = new TString(names_terms[0][term][0]->Data());
	names_terms_FVP[term][1][0] = new TString(names_terms[0][term][0]->Data());
	names_terms_FVP[term][0][0]->Append("_4VectProd1Ideal");
	names_terms_FVP[term][1][0]->Append("_4VectProd1Smeared");
	//
	names_terms_FVP[term][0][1] = new TString(names_terms[0][term][1]->Data());
	names_terms_FVP[term][1][1] = new TString(names_terms[0][term][1]->Data());
	names_terms_FVP[term][0][1]->Append("_4VectProd2Ideal");
	names_terms_FVP[term][1][1]->Append("_4VectProd2Smeared");
	//
	names_terms_FVP_K[term][0][0] = new TString(names_terms[0][term][1]->Data());// Prod1 SumK
	names_terms_FVP_K[term][1][0] = new TString(names_terms[0][term][1]->Data());// Prod1 EnK
	names_terms_FVP_K[term][0][1] = new TString(names_terms[0][term][1]->Data());// Prod2 SumK
	names_terms_FVP_K[term][1][1] = new TString(names_terms[0][term][1]->Data());// Prod2 EnK
	if(term==0) {
	  names_terms_FVP_K[term][0][0]->Append("_4VectProd1SumK3");
	  names_terms_FVP_K[term][1][0]->Append("_4VectProd1EnK3");
	  names_terms_FVP_K[term][0][1]->Append("_4VectProd2SumK3");
	  names_terms_FVP_K[term][1][1]->Append("_4VectProd2EnK3");
	  names_terms_FVPTPN[0][0] = new TString(TPNnameBase->Data());
	  names_terms_FVPTPN[0][1] = new TString(TPNnameBase->Data());
	  names_terms_FVPTPN[1][0] = new TString(TPNnameBase->Data());
	  names_terms_FVPTPN[1][1] = new TString(TPNnameBase->Data());
	  names_terms_FVPTPN[0][0]->Append("TPN_0_4VectProd1Ideal");
	  names_terms_FVPTPN[0][1]->Append("TPN_0_4VectProd2Ideal");
	  names_terms_FVPTPN[1][0]->Append("TPN_0_4VectProd1Smeared");
	  names_terms_FVPTPN[1][1]->Append("TPN_0_4VectProd2Smeared");
	}
	if(term>0 && term<4){
	  names_terms_FVP_K[term][0][0]->Append("_4VectProd1SumK2");
	  names_terms_FVP_K[term][1][0]->Append("_4VectProd1EnK2");
	  names_terms_FVP_K[term][0][1]->Append("_4VectProd2SumK2");
	  names_terms_FVP_K[term][1][1]->Append("_4VectProd2EnK2");
	}
	
	names_terms[0][term][0]->Append("_Ideal");
	names_terms[0][term][1]->Append("_Smeared");
	
	//
	if(ch==0) {MixedCbin1=0; MixedCbin2=0; MixedCbin3=1;}
	else {MixedCbin1=0; MixedCbin2=1; MixedCbin3=1;}
	names_terms[1][term][0]=new TString("PairCut3_Charge1_"); *names_terms[1][term][0] += MixedCbin1;
	names_terms[1][term][0]->Append("_Charge2_"); *names_terms[1][term][0] += MixedCbin2;
	names_terms[1][term][0]->Append("_Charge3_"); *names_terms[1][term][0] += MixedCbin3;
	names_terms[1][term][0]->Append("_SC_0_M_");
	*names_terms[1][term][0] += mb;
	names_terms[1][term][0]->Append("_ED_0_Term_"); *names_terms[1][term][0] += term+1;
	names_terms[1][term][1]=new TString(names_terms[1][term][0]->Data());
	if(term==0) {
	  names_termsK3[1][0]=new TString(names_terms[1][term][0]->Data());
	  names_termsK3[1][1]=new TString(names_terms[1][term][0]->Data());
	  names_termsK3[1][0]->Append("_SumK3");
	  names_termsK3[1][1]->Append("_EnK3");
	}	
	names_terms[1][term][0]->Append("_Ideal");
	names_terms[1][term][1]->Append("_Smeared");
      }// term loop
    
      
      TString *base1 = new TString("Explicit2_Charge1_"); *base1 += ch; 
      base1->Append("_Charge2_"); *base1 += ch; base1->Append("_SC_0_M_");
      names[0][0]=new TString(base1->Data());
      *names[0][0] += mb;
      names[0][0]->Append("_ED_0_Term_1_Ideal");
      names[0][1]=new TString(base1->Data());
      *names[0][1] += mb;
      names[0][1]->Append("_ED_0_Term_2_Ideal");
      names[0][2]=new TString(base1->Data());
      *names[0][2] += mb;
      names[0][2]->Append("_ED_0_Term_1_Smeared");
      names[0][3]=new TString(base1->Data());
      *names[0][3] += mb;
      names[0][3]->Append("_ED_0_Term_2_Smeared");
      // mixed charge 2-particle case is double counted here.  It's ok (only the stat errors are wrong...not used anyway).
      TString *base2 = new TString("Explicit2_Charge1_0_Charge2_1"); base2->Append("_SC_0_M_");
      names[1][0]=new TString(base2->Data());
      *names[1][0] += mb;
      names[1][0]->Append("_ED_0_Term_1_Ideal");
      names[1][1]=new TString(base2->Data());
      *names[1][1] += mb;
      names[1][1]->Append("_ED_0_Term_2_Ideal");
      names[1][2]=new TString(base2->Data());
      *names[1][2] += mb;
      names[1][2]->Append("_ED_0_Term_1_Smeared");
      names[1][3]=new TString(base2->Data());
      *names[1][3] += mb;
      names[1][3]->Append("_ED_0_Term_2_Smeared");
      //
      if(mb<2) list=(TList*)list1->Clone();
      //else list=(TList*)list2->Clone();
      
      
      for(int ChProd=0; ChProd<2; ChProd++){// charge product loop
	//
	num_i_2[ChProd]=(TH2D*)list->FindObject(names[ChProd][0]->Data());
	den_i_2[ChProd]=(TH2D*)list->FindObject(names[ChProd][1]->Data());
	num_s_2[ChProd]=(TH2D*)list->FindObject(names[ChProd][2]->Data());
	den_s_2[ChProd]=(TH2D*)list->FindObject(names[ChProd][3]->Data());
	//
	
	if(ch==0 && mb==0){
	  for(int term=0; term<5; term++){
	    merged_terms_i[ChProd][term] = (TH3D*)(list->FindObject(names_terms[ChProd][term][0]->Data()))->Clone();
	    merged_terms_s[ChProd][term] = (TH3D*)(list->FindObject(names_terms[ChProd][term][1]->Data()))->Clone();

	    if(term==0){
	      merged_terms_SumK3[ChProd] = (TH3D*)(list->FindObject(names_termsK3[ChProd][0]->Data()))->Clone();
	      merged_terms_EnK3[ChProd] = (TH3D*)(list->FindObject(names_termsK3[ChProd][1]->Data()))->Clone();
	    }
	    if(ChProd==0) {
	      //merged_terms_FVP_i[term][0] = (TH3D*)(list->FindObject(names_terms_FVP[term][0][0]->Data()))->Clone();// ideal, Prod 1
	      //merged_terms_FVP_s[term][0] = (TH3D*)(list->FindObject(names_terms_FVP[term][1][0]->Data()))->Clone();// smeared, Prod 1
	      //merged_terms_FVP_i[term][1] = (TH3D*)(list->FindObject(names_terms_FVP[term][0][1]->Data()))->Clone();// ideal, Prod 2
	      //merged_terms_FVP_s[term][1] = (TH3D*)(list->FindObject(names_terms_FVP[term][1][1]->Data()))->Clone();// smeared, Prod 2
	      if(term !=4){
		//merged_terms_FVP_SumK[term][0] = (TH3D*)(list->FindObject(names_terms_FVP_K[term][0][0]->Data()))->Clone();// Prod 1 Sum K
		//merged_terms_FVP_EnK[term][0] = (TH3D*)(list->FindObject(names_terms_FVP_K[term][1][0]->Data()))->Clone();// Prod 1 En K
		//merged_terms_FVP_SumK[term][1] = (TH3D*)(list->FindObject(names_terms_FVP_K[term][0][1]->Data()))->Clone();// Prod 2 Sum K
		//merged_terms_FVP_EnK[term][1] = (TH3D*)(list->FindObject(names_terms_FVP_K[term][1][1]->Data()))->Clone();// Prod 2 En K
	      }
	      if(term==0){
		//merged_terms_FVPTPN_i[0] = (TH3D*)(list->FindObject(names_terms_FVPTPN[0][0]->Data()))->Clone();// ideal, Prod 1 
		//merged_terms_FVPTPN_s[0] = (TH3D*)(list->FindObject(names_terms_FVPTPN[1][0]->Data()))->Clone();// smeared, Prod 1 
		//merged_terms_FVPTPN_i[1] = (TH3D*)(list->FindObject(names_terms_FVPTPN[0][1]->Data()))->Clone();// ideal, Prod 2
		//merged_terms_FVPTPN_s[1] = (TH3D*)(list->FindObject(names_terms_FVPTPN[1][1]->Data()))->Clone();// smeared, Prod 2
	      }
	    }
	  }
	}else{// ch==1 || mb!=0
	  for(int term=0; term<5; term++){
	    if(ChProd==0) {
	      merged_terms_i[0][term]->Add((TH3D*)(list->FindObject(names_terms[0][term][0]->Data())));
	      merged_terms_s[0][term]->Add((TH3D*)(list->FindObject(names_terms[0][term][1]->Data())));
	      if(term==0) {
		merged_terms_SumK3[0]->Add((TH3D*)(list->FindObject(names_termsK3[0][0]->Data())));
		merged_terms_EnK3[0]->Add((TH3D*)(list->FindObject(names_termsK3[0][1]->Data())));
		/*merged_terms_FVPTPN_i[0]->Add((TH3D*)(list->FindObject(names_terms_FVPTPN[0][0]->Data())));
		merged_terms_FVPTPN_i[1]->Add((TH3D*)(list->FindObject(names_terms_FVPTPN[0][1]->Data())));
		merged_terms_FVPTPN_s[0]->Add((TH3D*)(list->FindObject(names_terms_FVPTPN[1][0]->Data())));
		merged_terms_FVPTPN_s[1]->Add((TH3D*)(list->FindObject(names_terms_FVPTPN[1][1]->Data())));*/
	      }
	      /* merged_terms_FVP_i[term][0]->Add((TH3D*)(list->FindObject(names_terms_FVP[term][0][0]->Data())));
	      merged_terms_FVP_s[term][0]->Add((TH3D*)(list->FindObject(names_terms_FVP[term][1][0]->Data())));
	      merged_terms_FVP_i[term][1]->Add((TH3D*)(list->FindObject(names_terms_FVP[term][0][1]->Data())));
	      merged_terms_FVP_s[term][1]->Add((TH3D*)(list->FindObject(names_terms_FVP[term][1][1]->Data())));*/
	      
	      if(term<4){
		/*merged_terms_FVP_SumK[term][0]->Add((TH3D*)(list->FindObject(names_terms_FVP_K[term][0][0]->Data())));
		merged_terms_FVP_EnK[term][0]->Add((TH3D*)(list->FindObject(names_terms_FVP_K[term][1][0]->Data())));
		merged_terms_FVP_SumK[term][1]->Add((TH3D*)(list->FindObject(names_terms_FVP_K[term][0][1]->Data())));
		merged_terms_FVP_EnK[term][1]->Add((TH3D*)(list->FindObject(names_terms_FVP_K[term][1][1]->Data())));*/
	      }
	    }else {// different placement of same-charge pair
	      int termOther=term;
	      if(ch==1){
		if(term==0) termOther=0;
		else if(term==1) termOther=3;
		else if(term==2) termOther=2;
		else if(term==3) termOther=1;
		else termOther=4;
	      }
	      TH3D *temphisto_i=((TH3D*)(list->FindObject(names_terms[1][termOther][0]->Data())));
	      TH3D *temphisto_s=((TH3D*)(list->FindObject(names_terms[1][termOther][1]->Data())));
	      TH3D *temphisto_SumK3=((TH3D*)(list->FindObject(names_termsK3[1][0]->Data())));
	      TH3D *temphisto_EnK3=((TH3D*)(list->FindObject(names_termsK3[1][1]->Data())));
	      for(int bin1=1; bin1<merged_terms_i[1][term]->GetNbinsX(); bin1++){
		for(int bin2=1; bin2<merged_terms_i[1][term]->GetNbinsY(); bin2++){
		  for(int bin3=1; bin3<merged_terms_i[1][term]->GetNbinsZ(); bin3++){
		    
		    if(ch==0){
		      merged_terms_i[1][term]->SetBinContent(bin1,bin2,bin3, merged_terms_i[1][term]->GetBinContent(bin1,bin2,bin3)+temphisto_i->GetBinContent(bin1,bin2,bin3));
		      merged_terms_s[1][term]->SetBinContent(bin1,bin2,bin3, merged_terms_s[1][term]->GetBinContent(bin1,bin2,bin3)+temphisto_s->GetBinContent(bin1,bin2,bin3));
		    }else {
		      merged_terms_i[1][term]->SetBinContent(bin1,bin2,bin3, merged_terms_i[1][term]->GetBinContent(bin1,bin2,bin3)+temphisto_i->GetBinContent(bin3,bin2,bin1));
		      merged_terms_s[1][term]->SetBinContent(bin1,bin2,bin3, merged_terms_s[1][term]->GetBinContent(bin1,bin2,bin3)+temphisto_s->GetBinContent(bin3,bin2,bin1));
		    }
		    
		    //
		    if(term==0){
		      //merged_terms_SumK3[1]->SetBinContent(bin1,bin2,bin3, merged_terms_SumK3[1]->GetBinContent(bin1,bin2,bin3)+temphisto_SumK3->GetBinContent(bin3,bin2,bin1));
		      //merged_terms_EnK3[1]->SetBinContent(bin1,bin2,bin3, merged_terms_EnK3[1]->GetBinContent(bin1,bin2,bin3)+temphisto_EnK3->GetBinContent(bin3,bin2,bin1));
		    }
		  }
		}
	      }
	    }
	  }
	}// else for 3-particle part
	
	// 2-particles
	if(ch==0 && mb==0) {
	  merged_num_i_2[ChProd] = (TH2D*)num_i_2[ChProd]->Clone();
	  merged_den_i_2[ChProd] = (TH2D*)den_i_2[ChProd]->Clone();
	  merged_num_s_2[ChProd] = (TH2D*)num_s_2[ChProd]->Clone();
	  merged_den_s_2[ChProd] = (TH2D*)den_s_2[ChProd]->Clone();
	}else {
	  merged_num_i_2[ChProd]->Add(num_i_2[ChProd]);
	  merged_den_i_2[ChProd]->Add(den_i_2[ChProd]);
	  merged_num_s_2[ChProd]->Add(num_s_2[ChProd]);
	  merged_den_s_2[ChProd]->Add(den_s_2[ChProd]);
	  
	}
      }// ch loop
      
      
    }// mb loop
    
  }// ChProd loop
  
  
  cout<<"Start Writing to Output File"<<endl;
  TFile *outfile = new TFile("MomResFile_temp.root","RECREATE");
  TFile *outfileOffline = new TFile("MomResFile_Offline_temp.root","RECREATE");
  int NbinsX = merged_num_i_2[0]->GetNbinsX();
  int NbinsY = merged_num_i_2[0]->GetNbinsY();
  float Xlow = merged_num_i_2[0]->GetXaxis()->GetBinLowEdge(1);
  float Xhigh = merged_num_i_2[0]->GetXaxis()->GetBinUpEdge(NbinsX);
  float Ylow = merged_num_i_2[0]->GetYaxis()->GetBinLowEdge(1);
  float Yhigh = merged_num_i_2[0]->GetYaxis()->GetBinUpEdge(NbinsY);
  TH2D *MomResHisto_pp = new TH2D("MomResHisto_pp","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  TH2D *MomResHisto_mp = new TH2D("MomResHisto_mp","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  //
  TH2D *MomResHisto_pp_term1 = new TH2D("MomResHisto_pp_term1","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  TH2D *MomResHisto_pp_term2 = new TH2D("MomResHisto_pp_term2","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  TH2D *MomResHisto_mp_term1 = new TH2D("MomResHisto_mp_term1","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  TH2D *MomResHisto_mp_term2 = new TH2D("MomResHisto_mp_term2","",NbinsX,Xlow,Xhigh, NbinsY,Ylow,Yhigh);
  //
  // 3d part
  int NbinsX_3d = merged_terms_i[0][0]->GetNbinsX();
  float Xlow_3d = merged_terms_i[0][0]->GetXaxis()->GetBinLowEdge(1);
  float Xhigh_3d = merged_terms_i[0][0]->GetXaxis()->GetBinUpEdge(NbinsX_3d);
  TH3D *MomResHisto_3d_SC[5];
  TH3D *MomResHisto_3d_MC[5];
  TH1D *MomResHisto_1d_SC[5];
  TH1D *MomResHisto_1d_MC[5];
  TH3D *MomResHisto_3d_SC_K3;
  TH3D *MomResHisto_3d_MC_K3;
  // 4vectProd
  /*const int NEdges_FVP_temp = merged_terms_FVP_i[0][0]->GetNbinsX() + 1;
  const int NEdges_FVP = NEdges_FVP_temp;
  double Edges[NEdges_FVP]={0};
  for(int edg=0; edg<NEdges_FVP; edg++){
    Edges[edg] = merged_terms_FVP_i[0][0]->GetXaxis()->GetBinLowEdge(edg+1);
    //cout<<Edges[edg]<<endl;
  }
  TH3D *MomResHisto_3d_FVP1[5];
  TH3D *MomResHisto_3d_FVP2[5];
  TH3D *MomResHisto_3d_FVP1K[5];
  TH3D *MomResHisto_3d_FVP2K[5];
  TH3D *MomResHisto_TPN_FVP1;
  TH3D *MomResHisto_TPN_FVP2;*/


  for(int term=0; term<5; term++){
    TString *nameSC = new TString("MomResHisto3_SC_term");
    TString *nameMC = new TString("MomResHisto3_MC_term");
    TString *nameSC1D = new TString("MomResHisto1D_SC_term");
    TString *nameMC1D = new TString("MomResHisto1D_MC_term");
    TString *nameFVP1 = new TString("MomResHisto3_FVP1_term");
    TString *nameFVP2 = new TString("MomResHisto3_FVP2_term");
    TString *nameFVP1K = new TString("MomResHisto3_FVP1K_term");
    TString *nameFVP2K = new TString("MomResHisto3_FVP2K_term");
    *nameSC += term+1;
    *nameMC += term+1;
    *nameSC1D += term+1;
    *nameMC1D += term+1;
    *nameFVP1 += term+1;
    *nameFVP2 += term+1;
    *nameFVP1K += term+1;
    *nameFVP2K += term+1;
    nameSC->Append("_M");
    nameMC->Append("_M");
    nameSC1D->Append("_M");
    nameMC1D->Append("_M");
    nameFVP1->Append("_M");
    nameFVP2->Append("_M");
    nameFVP1K->Append("_M");
    nameFVP2K->Append("_M");
    *nameSC += MBmerged;
    *nameMC += MBmerged;
    *nameSC1D += MBmerged;
    *nameMC1D += MBmerged;
    *nameFVP1 += MBmerged;
    *nameFVP2 += MBmerged;
    *nameFVP1K += MBmerged;
    *nameFVP2K += MBmerged;
    MomResHisto_3d_SC[term] = new TH3D(nameSC->Data(),"",NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d);
    MomResHisto_3d_MC[term] = new TH3D(nameMC->Data(),"",NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d);
    MomResHisto_1d_SC[term] = new TH1D(nameSC1D->Data(),"",20,0,0.2);
    MomResHisto_1d_MC[term] = new TH1D(nameMC1D->Data(),"",20,0,0.2);
    /*MomResHisto_3d_FVP1[term] = new TH3D(nameFVP1->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);
      MomResHisto_3d_FVP2[term] = new TH3D(nameFVP2->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);
      MomResHisto_3d_FVP1K[term] = new TH3D(nameFVP1K->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);
      MomResHisto_3d_FVP2K[term] = new TH3D(nameFVP2K->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);*/
  }
  //TString *nameTPN_FVP1 = new TString("MomResHisto3_TPN_FVP1_M");
  //TString *nameTPN_FVP2 = new TString("MomResHisto3_TPN_FVP2_M");
  TString *nameSC_K3 = new TString("AvgK3_SC_M");
  TString *nameMC_K3 = new TString("AvgK3_MC_M");
  *nameSC_K3 += MBmerged;
  *nameMC_K3 += MBmerged;
  //*nameTPN_FVP1 += MBmerged;
  //*nameTPN_FVP2 += MBmerged;
  MomResHisto_3d_SC_K3 = new TH3D(nameSC_K3->Data(),"",NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d);
  MomResHisto_3d_MC_K3 = new TH3D(nameMC_K3->Data(),"",NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d, NbinsX_3d,Xlow_3d,Xhigh_3d);
  //MomResHisto_TPN_FVP1 = new TH3D(nameTPN_FVP1->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);
  //MomResHisto_TPN_FVP2 = new TH3D(nameTPN_FVP2->Data(),"",NEdges_FVP-1,Edges, NEdges_FVP-1,Edges, NEdges_FVP-1,Edges);
  
  
  
  for(int term=0; term<5; term++){
    
    // 1D Q3 projection
    double Sum_i_SC[20]={0};
    double Sum_i_MC[20]={0};
    double Sum_s_SC[20]={0};
    double Sum_s_MC[20]={0};
    
    for(int i=1; i<NbinsX_3d; i++){
      for(int j=1; j<NbinsX_3d; j++){
	for(int l=1; l<NbinsX_3d; l++){
	  //if((i+j)==2 || (i+l)==2 || (j+l)==2) continue;
	  double q3 = sqrt(pow(merged_terms_i[0][term]->GetXaxis()->GetBinCenter(i+1),2)+pow(merged_terms_i[0][term]->GetYaxis()->GetBinCenter(j+1),2)+pow(merged_terms_i[0][term]->GetZaxis()->GetBinCenter(l+1),2));
	  int q3bin = MomResHisto_1d_SC[term]->GetXaxis()->FindBin(q3);
	  if(q3bin>20) continue;
	  Sum_i_SC[q3bin-1] += merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1);
	  Sum_i_MC[q3bin-1] += merged_terms_i[1][term]->GetBinContent(i+1,j+1,l+1);
	  Sum_s_SC[q3bin-1] += merged_terms_s[0][term]->GetBinContent(i+1,j+1,l+1);
	  Sum_s_MC[q3bin-1] += merged_terms_s[1][term]->GetBinContent(i+1,j+1,l+1);
	}
      }
    }
    for(int i=0; i<20; i++){
      
      if(Sum_s_SC[i]>0) MomResHisto_1d_SC[term]->SetBinContent(i+1, Sum_i_SC[i]/Sum_s_SC[i]);
      if(Sum_s_MC[i]>0) MomResHisto_1d_MC[term]->SetBinContent(i+1, Sum_i_MC[i]/Sum_s_MC[i]);
    }

    // full 3d method
    merged_terms_i[0][term]->Divide(merged_terms_s[0][term]);
    merged_terms_i[1][term]->Divide(merged_terms_s[1][term]);
    //merged_terms_FVP_i[term][0]->Divide(merged_terms_FVP_s[term][0]);// Prod 1
    //merged_terms_FVP_i[term][1]->Divide(merged_terms_FVP_s[term][1]);// Prod 2
    if(term < 4){
      //merged_terms_FVP_SumK[term][0]->Divide(merged_terms_FVP_EnK[term][0]);// Prod 1
      //merged_terms_FVP_SumK[term][1]->Divide(merged_terms_FVP_EnK[term][1]);// Prod 2
    }
  }
  //merged_terms_FVPTPN_i[0]->Divide(merged_terms_FVPTPN_s[0]);// Prod 1
  //merged_terms_FVPTPN_i[1]->Divide(merged_terms_FVPTPN_s[1]);// Prod 2
  merged_terms_SumK3[0]->Divide(merged_terms_EnK3[0]);// Avg SC K3 for Q12,Q13,Q23
  merged_terms_SumK3[1]->Divide(merged_terms_EnK3[1]);// Avg MC K3 for Q12,Q13,Q23
  
  
  
  for(int term=0; term<5; term++){
    for(int i=0; i<NbinsX_3d; i++){
      for(int j=0; j<NbinsX_3d; j++){
	for(int l=0; l<NbinsX_3d; l++){
	  
	  double MomResCorr=0; double terms=0;
	  if(term==0 || term==4){
	    if(merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(i+1,l+1,j+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,l+1,j+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(j+1,i+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(j+1,i+1,l+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(j+1,l+1,i+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(j+1,l+1,i+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(l+1,i+1,j+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(l+1,i+1,j+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(l+1,j+1,i+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(l+1,j+1,i+1); terms++;}
	  }else if(term==1){
	    if(merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(i+1,l+1,j+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,l+1,j+1); terms++;}
	  }else if(term==2){
	    if(merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(l+1,j+1,i+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(l+1,j+1,i+1); terms++;}
	  }else {
	    if(merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	    if(merged_terms_i[0][term]->GetBinContent(j+1,i+1,l+1)>0) {MomResCorr += merged_terms_i[0][term]->GetBinContent(j+1,i+1,l+1); terms++;}
	  }
	  
	  if(terms > 0) MomResCorr /= terms;
	  MomResHisto_3d_SC[term]->SetBinContent(i+1,j+1,l+1, MomResCorr);

	  // Mixed-Charge
	  MomResCorr=0; terms=0;
	  if(term==0 || term==1 || term==4){
	    if(merged_terms_i[1][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[1][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	    if(merged_terms_i[1][term]->GetBinContent(i+1,l+1,j+1)>0) {MomResCorr += merged_terms_i[1][term]->GetBinContent(i+1,l+1,j+1); terms++;}
	  }else{
	    if(merged_terms_i[1][term]->GetBinContent(i+1,j+1,l+1)>0) {MomResCorr += merged_terms_i[1][term]->GetBinContent(i+1,j+1,l+1); terms++;}
	  }
	  
	  if(terms > 0) MomResCorr /= terms;
	  MomResHisto_3d_MC[term]->SetBinContent(i+1,j+1,l+1, MomResCorr);
	  
	  if(term==0){
	    MomResHisto_3d_SC_K3->SetBinContent(i+1,j+1,l+1, merged_terms_SumK3[0]->GetBinContent(i+1,j+1,l+1));
	    MomResHisto_3d_MC_K3->SetBinContent(i+1,j+1,l+1, merged_terms_SumK3[1]->GetBinContent(i+1,j+1,l+1));
	  }
	  // condition edge effects
	  if(MomResHisto_3d_SC[term]->GetBinContent(i+1,j+1,l+1) > 10) MomResHisto_3d_SC[term]->SetBinContent(i+1,j+1,l+1, 1.0);
	  if(MomResHisto_3d_MC[term]->GetBinContent(i+1,j+1,l+1) > 10) MomResHisto_3d_MC[term]->SetBinContent(i+1,j+1,l+1, 1.0);
	}
      }
    }
  }// term
  
  // FVP
  /*for(int mb=0; mb<1; mb++){
    for(int term=0; term<5; term++){
      for(int i=0; i<NEdges_FVP-1; i++){
	for(int j=0; j<NEdges_FVP-1; j++){
	  for(int l=0; l<NEdges_FVP-1; l++){
	    MomResHisto_3d_FVP1[term]->SetBinContent(i+1,j+1,l+1, merged_terms_FVP_i[term][0]->GetBinContent(i+1,j+1,l+1));
	    MomResHisto_3d_FVP2[term]->SetBinContent(i+1,j+1,l+1, merged_terms_FVP_i[term][1]->GetBinContent(i+1,j+1,l+1));
	    if(term < 4){
	      MomResHisto_3d_FVP1K[term]->SetBinContent(i+1,j+1,l+1, merged_terms_FVP_SumK[term][0]->GetBinContent(i+1,j+1,l+1));
	      MomResHisto_3d_FVP2K[term]->SetBinContent(i+1,j+1,l+1, merged_terms_FVP_SumK[term][1]->GetBinContent(i+1,j+1,l+1));
	    }
	    if(term==0){
	      MomResHisto_TPN_FVP1->SetBinContent(i+1,j+1,l+1, merged_terms_FVPTPN_i[0]->GetBinContent(i+1,j+1,l+1));
	      MomResHisto_TPN_FVP2->SetBinContent(i+1,j+1,l+1, merged_terms_FVPTPN_i[1]->GetBinContent(i+1,j+1,l+1));
	    }
	    // condition edge effects
	    if(MomResHisto_3d_FVP1[term]->GetBinContent(i+1,j+1,l+1) > 10) MomResHisto_3d_FVP1[term]->SetBinContent(i+1,j+1,l+1, 1.0);
	    if(MomResHisto_3d_FVP2[term]->GetBinContent(i+1,j+1,l+1) > 10) MomResHisto_3d_FVP2[term]->SetBinContent(i+1,j+1,l+1, 1.0);
	  }
	}
      }
    }
    }*/
  
  // 2-particles
  for(int ChProd=0; ChProd<2; ChProd++){
    for(int i=0; i<NbinsX; i++){
      for(int j=0; j<NbinsY; j++){
	
	double weight = 1.0, weight_term1=1.0, weight_term2=1.0;
	double weight_e = 0, weight_term1_e=0, weight_term2_e=0;
	bool goodbin=kTRUE;
	if(merged_num_i_2[ChProd]->GetBinContent(i+1,j+1) == 0) goodbin=kFALSE;
	if(merged_den_i_2[ChProd]->GetBinContent(i+1,j+1) == 0) goodbin=kFALSE;
	if(merged_num_s_2[ChProd]->GetBinContent(i+1,j+1) == 0) goodbin=kFALSE;
	if(merged_den_s_2[ChProd]->GetBinContent(i+1,j+1) == 0) goodbin=kFALSE;
	
	if(goodbin==kTRUE){
	  double NUM_i=merged_num_i_2[ChProd]->GetBinContent(i+1,j+1);
	  double NUM_s=merged_num_s_2[ChProd]->GetBinContent(i+1,j+1);
	  double DEN_i=merged_den_i_2[ChProd]->GetBinContent(i+1,j+1);
	  double DEN_s=merged_den_s_2[ChProd]->GetBinContent(i+1,j+1);
	  weight = (NUM_i/DEN_i)/(NUM_s/DEN_s);
	  weight_e = pow((sqrt(NUM_s)*NUM_i/DEN_i)/(pow(NUM_s,2)/DEN_s),2);// approximate. does not take lambda into account!!
	  weight_e += pow((sqrt(DEN_s)*NUM_i/DEN_i)/(NUM_s),2);
	  weight_e = sqrt(weight_e); 
	  //
	  weight_term1 = NUM_i/NUM_s;
	  weight_term1_e = pow(sqrt(NUM_i)/NUM_s,2);
	  weight_term1_e += pow(sqrt(NUM_s)*NUM_i/pow(NUM_s,2),2);
	  weight_term1_e = sqrt(weight_term1_e);
	  weight_term2 = DEN_i/DEN_s;
	  weight_term2_e = pow(sqrt(DEN_i)/DEN_s,2);
	  weight_term2_e += pow(sqrt(DEN_s)*DEN_i/pow(DEN_s,2),2);
	  weight_term2_e = sqrt(weight_term2_e);
	}
	
	if(ChProd==0) {
	  MomResHisto_pp->SetBinContent(i+1,j+1, weight);
	  MomResHisto_pp->SetBinError(i+1,j+1, weight_e);
	  MomResHisto_pp_term1->SetBinContent(i+1,j+1, weight_term1);
	  MomResHisto_pp_term2->SetBinContent(i+1,j+1, weight_term2);
	  MomResHisto_pp_term1->SetBinError(i+1,j+1, weight_term1_e);
	  MomResHisto_pp_term2->SetBinError(i+1,j+1, weight_term2_e);
	}else {
	  MomResHisto_mp->SetBinContent(i+1,j+1, weight); 
	  MomResHisto_mp->SetBinError(i+1,j+1, weight_e);
	  MomResHisto_mp_term1->SetBinContent(i+1,j+1, weight_term1);
	  MomResHisto_mp_term2->SetBinContent(i+1,j+1, weight_term2);
	  MomResHisto_mp_term1->SetBinError(i+1,j+1, weight_term1_e);
	  MomResHisto_mp_term2->SetBinError(i+1,j+1, weight_term2_e);
	}
      }
    }
  }
  

 
  outfile->cd();
  MomResHisto_pp->Write("MomResHisto_pp");
  MomResHisto_mp->Write("MomResHisto_mp");
  MomResHisto_pp_term1->Write("MomResHisto_pp_term1");
  MomResHisto_pp_term2->Write("MomResHisto_pp_term2");
  MomResHisto_mp_term1->Write("MomResHisto_mp_term1");
  MomResHisto_mp_term2->Write("MomResHisto_mp_term2");
  outfile->Close();

  outfileOffline->cd();
  if(MBmerged==0){
    MomResHisto_pp->Write("MomResHisto_pp");
    MomResHisto_mp->Write("MomResHisto_mp");
    MomResHisto_pp_term1->Write("MomResHisto_pp_term1");
    MomResHisto_pp_term2->Write("MomResHisto_pp_term2");
    MomResHisto_mp_term1->Write("MomResHisto_mp_term1");
    MomResHisto_mp_term2->Write("MomResHisto_mp_term2");
  }    
  //
  
  for(int term=0; term<5; term++){
    MomResHisto_3d_SC[term]->Write();
    MomResHisto_3d_MC[term]->Write();
    MomResHisto_1d_SC[term]->Write();
    MomResHisto_1d_MC[term]->Write();
    //MomResHisto_3d_FVP1[term]->Write();
    //MomResHisto_3d_FVP2[term]->Write();
    if(term < 4){
      //MomResHisto_3d_FVP1K[term]->Write();
      //MomResHisto_3d_FVP2K[term]->Write();
    }
  }
  MomResHisto_3d_SC_K3->Write();
  MomResHisto_3d_MC_K3->Write();
  //MomResHisto_TPN_FVP1->Write();
  //MomResHisto_TPN_FVP2->Write();
  
  
  outfileOffline->Close();
   
}
