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


void MakeWeightFile()
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);

  const int NumEDbins=2;
  bool EDBinning=kTRUE;
  //
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h.root","READ");
  //TFile *InputFile = new TFile("Results/PDC_11h_pT_0p2to1p0_FullRunWrongWeightsNoPadRowTTCandMCTTC_RawWeightFile.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h_0p02eta0p045phi_0p03eta0p067phi.root","READ");
  TFile *InputFile = new TFile("Results/RawWeightFile_11h_q2Binning_LowPtMultBinning.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h_81EMbins.root","READ");
  

  TDirectoryFile *tdir = (TDirectoryFile*)InputFile->Get("PWGCF.outputFourPionAnalysis.root");
  TList *MyList=(TList*)tdir->Get("FourPionOutput_2");
  //TList *MyList=(TList*)InputFile->Get("MyList");
  InputFile->Close();
  
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  cout<<Events->GetBinContent(1)<<endl;
  
 

  const int KtBins=4;//4
  const int KyBins=1;
  const int MBins=10;
  const int MBLimit=10;// 0-10. 1 for GenSignal, 10 for the rest!!!!!!!!!!!!!!!!!!!!!!!!!!

  TFile *OutFile = new TFile("WeightFile_temp.root","RECREATE");
  TH3F *WeightHistos[KtBins][MBins][NumEDbins];

 
  for(int ktB=1; ktB<=KtBins; ktB++){
    cout<<"kT Bin "<<ktB<<endl;
    for(int MB=1; MB<=MBins; MB++){
      for(int q2B=1; q2B<=NumEDbins; q2B++){
	TString *InNameNum = new TString("TPN_num_Kt_");
	*InNameNum += ktB-1;
	InNameNum->Append("_Ky_0_M_");
	if(MB<=MBLimit) *InNameNum += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *InNameNum += 1;
	InNameNum->Append("_ED_");
	if(EDBinning) *InNameNum += q2B-1;
	else *InNameNum += 0;// for non q2 binning
	//
	TString *InNameDen = new TString("TPN_den_Kt_");
	*InNameDen += ktB-1;
	InNameDen->Append("_Ky_0_M_");
	if(MB<=MBLimit) *InNameDen += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *InNameDen += 1;
	InNameDen->Append("_ED_");
	if(EDBinning) *InNameDen += q2B-1;
	else *InNameDen += 0;// for non q2 binning
	//
	TH3D *tempNum = (TH3D*)MyList->FindObject(InNameNum->Data());
	TH3D *tempDen = (TH3D*)MyList->FindObject(InNameDen->Data());
	//
	int NormBinStartOut = tempNum->GetXaxis()->FindBin(0.06);// 0.135.  0.06 as a systematic
	int NormBinStartSideLong = tempNum->GetXaxis()->FindBin(0.06);// 0.135.  0.06 as a systematic
	int NormBinEnd = tempNum->GetXaxis()->FindBin(0.08);// 0.2.  0.08 as a systematic
	double Norm = tempNum->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd);
	if(tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd) > 0){
	  Norm /= tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd);
	}else Norm=0;
	cout<<"Normalization = "<<Norm<<endl;
	cout<<tempNum->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd)<<"  "<<tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd)<<endl;
	//
	TString *OutNameWeight = new TString("Weight_Kt_");
	*OutNameWeight += ktB-1;
	OutNameWeight->Append("_Ky_0_M_");
	*OutNameWeight += MB-1;
	OutNameWeight->Append("_ED_");
	*OutNameWeight += q2B-1;
	
	int Nbins=tempNum->GetNbinsX();
	double QLimit = tempNum->GetXaxis()->GetBinUpEdge(Nbins);
	WeightHistos[ktB-1][MB-1][q2B-1] = new TH3F(OutNameWeight->Data(),"r3 Weights", Nbins,0,QLimit, Nbins,0,QLimit, Nbins,0,QLimit);
	WeightHistos[ktB-1][MB-1][q2B-1]->GetXaxis()->SetTitle("out");
	WeightHistos[ktB-1][MB-1][q2B-1]->GetYaxis()->SetTitle("side");
	WeightHistos[ktB-1][MB-1][q2B-1]->GetZaxis()->SetTitle("long");
	WeightHistos[ktB-1][MB-1][q2B-1]->GetXaxis()->SetTitleOffset(1.8);
	WeightHistos[ktB-1][MB-1][q2B-1]->GetYaxis()->SetTitleOffset(1.8);
	WeightHistos[ktB-1][MB-1][q2B-1]->GetZaxis()->SetTitleOffset(1.8);
	double LowQcount=0, HighQcount=0;
	for(int outB=1; outB<=Nbins; outB++){
	  for(int sideB=1; sideB<=Nbins; sideB++){
	    for(int longB=1; longB<=Nbins; longB++){
	      if(Norm==0) continue;
	      double weight=1, weight_e=0;
	      if(tempDen->GetBinContent(outB,sideB,longB) > 0 && tempNum->GetBinContent(outB,sideB,longB) > 0) {
		//if(outB==sideB && outB==longB) cout<<outB<<"  "<<tempNum->GetBinContent(outB,sideB,longB)<<"  "<<tempDen->GetBinContent(outB,sideB,longB)<<endl;
		weight = double(tempNum->GetBinContent(outB,sideB,longB))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm;
		weight_e = pow(sqrt(double(tempNum->GetBinContent(outB,sideB,longB)))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm,2);
		weight_e += pow(sqrt(double(tempDen->GetBinContent(outB,sideB,longB)))*double(tempNum->GetBinContent(outB,sideB,longB))/pow(double(tempDen->GetBinContent(outB,sideB,longB)),2) / Norm,2);
		weight_e = sqrt(weight_e);
		//if(weight < 0.8 && weight > 0.1) cout<<outB<<"  "<<sideB<<"  "<<longB<<"  "<<tempNum->GetBinContent(outB,sideB,longB)<<"  "<<tempDen->GetBinContent(outB,sideB,longB)<<endl;
		double qo = tempNum->GetXaxis()->GetBinCenter(outB);
		double qs = tempNum->GetYaxis()->GetBinCenter(sideB);
		double ql = tempNum->GetZaxis()->GetBinCenter(longB);
		double qmag = sqrt(pow(qo,2) + pow(qs,2) + pow(ql,2));
		if(qmag > 0.01 && qmag < 0.04) LowQcount++;
		if(qmag > 0.04 && qmag < 0.07) HighQcount++;
	      }
	      if(weight==0){
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, 0);
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinError(outB,sideB,longB, 0);
	      }else {
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, weight-1.0);// difference from unity
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinError(outB,sideB,longB, weight_e);
	      }
	      
	      
	    }
	  }
	}
	cout<<"PairCount   "<<LowQcount<<"  "<<HighQcount<<endl;
	WeightHistos[ktB-1][MB-1][q2B-1]->Write();
	
      }// q2B
    }// MB
  }// ktB

 

  int BOI1=2;
  int BOI2=2;
  TH3D *histoOI = WeightHistos[0][0][0]->Clone();
  histoOI->SetDirectory(0);
  //OutFile->Close();
  
  TH1D *pro = WeightHistos[0][0][0]->ProjectionX("pro",BOI1,BOI1,BOI2,BOI2);
  pro->Draw();

  //for(int bin=1; bin<=40; bin++) cout<<pro->GetBinContent(bin)<<", ";
  //cout<<endl;
  //cout<<endl;
  //for(int bin=1; bin<=40; bin++) cout<<pro->GetBinError(bin)<<", ";
  //cout<<endl;

  // qout ,2--2 other bin integration (lowest kT, M0)
  //double Ref[40]={0.186581, 0.180099, 0.16335, 0.137551, 0.103902, 0.0754772, 0.0381106, 0.0250045, 0.0194101, 0.00606144, 0.00175905, -0.00335014, 0.00347716, -0.0020038, 0.00276897, -0.0114744, -0.00412381, 0.00880619, -0.00533615, 0.000337426, -0.00603902, 0.00860179, -0.00154849, 0.00605317, -0.00328388, -0.0257134, 0.00816, 0.00588766, 0.00238239, -0.0208381, 0.00262856, 0.000531408, -9.19566e-05, -0.011538, -0.0216935, -0.0439466, 0, 0, 0, 0};
  //double Ref_e[40]={0.00258484, 0.00262131, 0.00265747, 0.00268686, 0.00271356, 0.00276084, 0.00277321, 0.00281418, 0.00285723, 0.00284373, 0.00286302, 0.00289162, 0.00297097, 0.00304969, 0.00316452, 0.00323752, 0.00339239, 0.0035809, 0.00369192, 0.00389189, 0.00406313, 0.00436902, 0.00457948, 0.00493212, 0.00524766, 0.00552003, 0.00621051, 0.0068126, 0.00754255, 0.00825506, 0.00960699, 0.0111344, 0.0131076, 0.0162349, 0.021695, 0.0372195, 0, 0, 0, 0};
  // qside, 2--2 other bin integration
  //double Ref[40]={0.180498, 0.180099, 0.185806, 0.163894, 0.132356, 0.098442, 0.0732086, 0.0419, 0.0273973, 0.0140072, 0.00491044, 0.00828187, -0.00346015, 0.00556346, 0.00280363, 0.000493745, 0.00483145, -0.000884939, 0.00825155, -0.00551758, -0.00208245, -0.00369946, -0.000581864, -0.000484325, -0.000560278, 0.00125169, 0.00402387, 0.00369215, 0.000496307, 0.00343727, 0.00382731, -0.00107187, -0.00187021, -0.00498415, 0.00312242, -0.000784691, 0.00153594, 0.00514852, -0.00480193, 0.00576707};
  //double Ref_e[40]={0.00263533, 0.00262131, 0.00262328, 0.00255935, 0.00250353, 0.00244688, 0.00240138, 0.00234001, 0.00231244, 0.00228984, 0.00226491, 0.00227351, 0.00224481, 0.00226501, 0.00224979, 0.00223847, 0.00224863, 0.0022289, 0.00224206, 0.00220784, 0.0022124, 0.00219928, 0.00219753, 0.0021926, 0.00218304, 0.00218404, 0.00217921, 0.00217676, 0.00216391, 0.00216305, 0.00215278, 0.00213398, 0.00212274, 0.00211043, 0.0021198, 0.00210246, 0.00210134, 0.00209683, 0.00207121, 0.00208596};
  // qlong, 2--2 other bin integration
  //double Ref[40]={0.19068, 0.180099, 0.16035, 0.118137, 0.0800825, 0.0484194, 0.0314626, 0.0140836, 0.015328, 0.00456526, 0.00461427, 0.00444591, -0.00385123, 0.0076453, 0.00418319, 0.00385833, 0.00192732, 0.0015057, 0.00151401, 0.00441839, 0.00255238, -0.00562173, -0.00434487, 0.00327595, -0.000923301, 0.00340199, -0.00116685, -0.00639669, -0.000535647, 0.000632789, -0.00431929, -0.00590313, -0.00334011, -0.00415945, 0.0113907, -0.00252384, 0.00152043, 0.000260006, -0.00476303, -0.00232199};
  //double Ref_e[40]={0.00394138, 0.00262131, 0.00257243, 0.00251757, 0.00247518, 0.00244153, 0.00243551, 0.00242434, 0.00245618, 0.00245905, 0.00248088, 0.00250852, 0.00251633, 0.00257305, 0.0025939, 0.00261776, 0.0026447, 0.00267184, 0.00269996, 0.00273682, 0.00275901, 0.00277727, 0.00280127, 0.00286128, 0.00288497, 0.00292533, 0.00294705, 0.00296652, 0.00302051, 0.00306191, 0.00308495, 0.00311484, 0.00315736, 0.00320054, 0.00328945, 0.00328307, 0.00333872, 0.00337698, 0.0034125, 0.00346363};
  //
  // qout ,2--2 other bin integration (highest kT, M0)
  //double Ref[40]={0.347024, 0.340981, 0.207019, 0.204354, 0.235698, 0.113875, 0.103287, 0.0390906, 0.049902, 0.0632777, 0.0196594, 0.0168798, -0.0258711, -0.0173084, -0.0255809, -0.03879, -0.0467173, -0.0285905, -0.0280937, -0.0144413, -0.0242922, -0.0183046, -0.0232596, -0.0242507, -0.0322753, -0.0337653, -0.0270012, -0.00878507, -0.00345743, -0.0201966, -0.023139, -0.0138009, -0.0236573, -0.0212273, -0.0261684, -0.014746, 0.00113992, -0.00296493, -0.00872464, -0.0149615};
  //double Ref_e[40]={0.0291587, 0.028769, 0.0247413, 0.0233308, 0.0228364, 0.0193632, 0.0182532, 0.0163304, 0.0156345, 0.0151442, 0.01384, 0.0133412, 0.0123078, 0.0116285, 0.0108021, 0.00981288, 0.00906128, 0.00856431, 0.00804651, 0.00772649, 0.00733884, 0.00705491, 0.0068498, 0.00661994, 0.00638918, 0.00626731, 0.00617582, 0.00618696, 0.00611707, 0.00593496, 0.00581311, 0.00576991, 0.00562799, 0.00553088, 0.00539531, 0.00533839, 0.00533781, 0.00523377, 0.00512349, 0.00503239};
  // qout ,2--2 other bin integration (kT1, M9)
  //double Ref[40]={0.274214, 0.266342, 0.262087, 0.258269, 0.342804, 0.271642, 0.262439, 0.191731, 0.161266, 0.0480234, 0.0678067, 0.0298723, 0.140997, 0.0225363, 0.0369282, 0.0474108, 0.00666558, -0.00719739, 0.0397449, -0.0154514, 0.0563656, -0.0283267, -0.0233019, 0.00502379, -0.0142704, 0.00171802, 0.00224826, -0.0103424, -0.00215456, -0.00290923, 0.0294002, 0.00502061, -0.0311593, -0.00107628, 0.0157102, 0.0242639, -0.00358242, 0.0273829, -0.00498812, -0.0191332};
  //double Ref_e[40]={0.0219786, 0.0219724, 0.0219296, 0.021582, 0.023235, 0.0219217, 0.0219606, 0.0209754, 0.0203662, 0.0185691, 0.0189119, 0.0183129, 0.0204752, 0.017965, 0.0179198, 0.0177246, 0.0170545, 0.0166412, 0.0172614, 0.0163279, 0.0174978, 0.016035, 0.0160211, 0.0165804, 0.0164245, 0.016433, 0.0167745, 0.0163847, 0.0166272, 0.0166435, 0.0172967, 0.0171336, 0.0165754, 0.0172471, 0.0177293, 0.0181375, 0.018193, 0.0189333, 0.0187681, 0.0189596};
  
  /*TH1D *Ratio = (TH1D*)pro->Clone();
  Ratio->SetTitle("C_{2}^{#pm#pm} ratio of nanoAOD160 to AOD160");
  Ratio->GetXaxis()->SetTitleOffset(0.9);
  Ratio->GetXaxis()->SetTitle("q_{out} GeV/c");
  for(int bin=1; bin<=40; bin++) {
    double value = (pro->GetBinContent(bin)+1) / (Ref[bin-1]+1);
    double value_e = sqrt( fabs(pow(pro->GetBinError(bin),2) - pow(Ref_e[bin-1],2)));
    Ratio->SetBinContent(bin, value);
    Ratio->SetBinError(bin, value_e);
  }
  Ratio->Draw();
  Ratio->Fit("pol0","","",0,0.2);
  */

  /*TFile *projections=new TFile("projections.root","RECREATE");
  histoOI->GetXaxis()->SetRange(BOI1,BOI1);
  TH2D *proYZ = (TH2D*)histoOI->Project3D("yz");
  histoOI->GetXaxis()->SetRange(1,40);
  histoOI->GetYaxis()->SetRange(BOI1,BOI1);
  TH2D *proXZ = (TH2D*)histoOI->Project3D("xz");
  histoOI->GetYaxis()->SetRange(1,40);
  histoOI->GetZaxis()->SetRange(BOI1,BOI1);
  TH2D *proXY = (TH2D*)histoOI->Project3D("xy");
  proYZ->GetXaxis()->SetRangeUser(0,0.1);
  proYZ->GetYaxis()->SetRangeUser(0,0.1);
  proXZ->GetXaxis()->SetRangeUser(0,0.1);
  proXZ->GetYaxis()->SetRangeUser(0,0.1);
  proXY->GetXaxis()->SetRangeUser(0,0.1);
  proXY->GetYaxis()->SetRangeUser(0,0.1);
  proXY->Draw("lego");
  
  proYZ->Write();
  proXZ->Write();
  proXY->Write();
  */

}
