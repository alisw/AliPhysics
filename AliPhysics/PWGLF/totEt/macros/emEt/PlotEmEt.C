
#ifndef __CINT__
#include "TGraphErrors.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <iostream>
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#endif
void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
}

void PlotEmEt(TString filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root", TString det = "Emcal", bool sim = false, double scaleFactor = 1, double pionScale = 1.0)
{
  Int_t colors[] = {TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
  Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};
  Float_t pionPlusEt[10] = {360.7,298.3,223.8,149.9,96.1, 58.1,32.4,16.4,7.3,2.7};
  Float_t pionMinusEt[10] ={363.7,300.4,225.4,150.5,96.6, 58.4,32.5,16.5,7.4,2.8};
  Float_t pionEtError[10] = {19.3, 15.3,11.3 ,7.5  , 4.8,  2.9, 1.6, 0.8,0.4,0.1};
  Float_t pionEt[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t emEtPerNpart[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t emEtPerNpartErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t emEt[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t emEtErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  for(int i=0;i<10;i++){
    pionEt[i] = (pionPlusEt[i]+pionMinusEt[i])/2.0;
    emEt[i] = pionEt[i]*1.085;
    emEtErr[i] = emEt[i]*TMath::Sqrt(TMath::Power(0.030/1.085,2)+TMath::Power(pionEtError[i]/pionEt[i],2));
  }


    Float_t x[20], xerr[20],  y_rel_error[20] ,yerr[20], xpion[20], xpionerr[20], ypion[20], ypionerr[20], ysim[20], yrec[20];
    Float_t nonLinError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t efficiencyError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t minEtErr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t minEtCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEtErr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEtValue[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEt500MeVErr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEt500MeVValue[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtErr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtValue[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtErr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtValue[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    TString centLabels[] = {"0-5%","5-10%","10-15%","15-20%","20-25%",  "25-30%","30-35%","35-40%","40-45%","45-50%",  "50-55%","55-60%","60-65%","65-70%","70-75%",  "75-80%","80-85%","85-90%","90-95%","95-100%"};

    Float_t nonLinErrorShort[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t efficiencyErrorShort[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t nonLinErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t efficiencyErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadErrorShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t hadCorrShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonErrorShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t kaonCorrShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t minEtErrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t minEtCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEtErrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEtValueShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEt500MeVErrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t rawEt500MeVValueShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtErrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtValueShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtErrShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t partialCorrEtValueShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtErrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtValueShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtErrShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtValueShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtErrPerNPartShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtValuePerNPartShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtErrPerNPartShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    Float_t corrEtValuePerNPartShort500MeV[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
    TString centLabelsShort[] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%"};

    Float_t hadErrorEmcal[20] = {9.96941,4.34033,3.2544,2.46096,1.85974,1.39188,1.04314,0.768989,0.566392,0.416278,0.303809,0.219583,0.169077,0.137327,0.115232,0.099452,0.095027,0.0887826,0.0645363,};
    Float_t hadCorrEmcal[20] = {58.0371,29.4045,23.0345,18.0506,14.0786,10.8591,8.31337,6.26566,4.66497,3.46897,2.5448,1.90246,1.5164,1.2924,1.15985,1.07691,1.03915,1.00897,0.944257,};
    Float_t hadError500MeVEmcal[20] = {5.19955,2.91943,2.30763,1.82901,1.44725,1.13405,0.880527,0.673811,0.506236,0.380598,0.278777,0.207644,0.158179,0.124682,0.101618,0.0854513,0.074277,0.0710671,0.0628612,};
    Float_t hadCorr500MeVEmcal[20] = {41.4917,18.5258,14.5182,11.3836,8.87434,6.84412,5.23258,3.94014,2.9284,2.17168,1.5966,1.2164,1.01307,0.92015,0.886337,0.876897,0.896635,0.912078,0.890611,};


    Float_t hadErrorPhos[20] = {2.17752,1.02329,0.808674,0.639466,0.507094,0.413639,0.333212,0.272335,0.229312,0.192292,0.16624,0.147969,0.132668,0.118544,0.117667,0.111382,0.109401,0.115973,0.111227,};
    Float_t hadCorrPhos[20] = {16.8125,8.3289,6.67917,5.38945,4.36484,3.61504,3.02595,2.6012,2.28153,2.03373,1.85627,1.72462,1.62872,1.54236,1.52404,1.48698,1.47255,1.48271,1.46966,};
    //Float_t hadError500MeVPhos[20] = {1.78823,0.965582,0.768239,0.61382,0.491313,0.398799,0.323919,0.267179,0.223729,0.188801,0.163673,0.143605,0.129474,0.117266,0.120714,0.124057,0.119597,0.162631,0.121791,};
    //Float_t hadCorr500MeVPhos[20] = {9.91472,4.0766,3.29702,2.68898,2.20738,1.86918,1.6142,1.44665,1.3227,1.23356,1.1716,1.13759,1.11538,1.08378,1.1192,1.12753,1.12033,1.18281,1.12409,};
    Float_t hadError500MeVPhos[20] = {1.63912,0.986321,0.777252,0.627139,0.507912,0.409661,0.337173,0.281528,0.237638,0.200728,0.176946,0.153582,0.135979,0.128235,0.114839,0.105581,0.103842,0.110116,0.127558,};
    Float_t hadCorr500MeVPhos[20] = {10.3323,4.33103,3.44988,2.82032,2.31597,1.95055,1.70481,1.51718,1.38961,1.29369,1.22474,1.18656,1.16962,1.13756,1.19274,1.16588,1.17661,1.19598,1.23654,};

    //PHOS kaon deposits corresponding to minimum energy cut of 250 MeV
    Float_t kaonCorrPhos[10] = {0.522969,1.0701,0.33745,0.386312,0.091973,0.0874692,0.0599575,0.0425214,0.0178472,0.0543614};
    Float_t kaonErrorPhos[10] = {0.0443465,0.112813,0.0256754,0.0344671,0.00641868,0.00744562,0.00488688,0.00348225,0.00161517,0.00590774};
    Float_t kaonCorrPhos500MeV[10] = {0.33013,0.626137,0.20589,0.28063,0.0917355,0.0555848,0.0490011,0.0282463,0.0122385,0.0325976};
    Float_t kaonErrorPhos500MeV[10] = {0.0279942,0.0499577,0.0156654,0.0250381,0.00640211,0.00473154,0.00399387,0.0023132,0.00110758,0.00354255};

    //EMCal kaon deposits corresponding to minimum energy cut of 300 MeV
    Float_t kaonCorrEmcal[10] = {14.8391,14.6397,4.80481,3.72969,1.76575,1.49971,0.929531,0.615978,0.223484,0.668516};
    Float_t kaonErrorEmcal[10] = {1.76834,1.89738,0.460842,0.398095,0.161214,0.155527,0.0938913,0.0626049,0.0243395,0.083329};
    Float_t kaonCorrEmcal500MeV[10] = {12.586,13.8596,4.30586,3.17298,1.45638,1.21947,0.775579,0.516616,0.208884,0.626861};
    Float_t kaonErrorEmcal500MeV[10] = {1.56637,1.82579,0.412986,0.338674,0.132968,0.126465,0.0783407,0.0525063,0.0227494,0.0781368};

    //Cuts for 250 MeV min ET cut
    Float_t minEtCorrPhos[10] = {0.757,0.756,0.756,0.753,0.747,0.739,0.729,0.718,0.706,0.690};
    Float_t minEtErrPhos[10] = {0.027,0.027,0.028,0.028,0.029,0.030,0.031,0.034,0.034,0.034};
    //Cuts for 300 MeV min ET cut
    Float_t minEtCorrEmcal[10] = {0.691,0.691,0.691,0.688,0.681,0.673,0.662,0.650,0.637,0.619};
    Float_t minEtErrEmcal[10] = {0.031,0.031,0.032,0.032,0.032,0.034,0.035,0.037,0.037,0.038};
    //Cuts for 500 MeV min ET cut
    Float_t minEtCorr500MeV[10] = {0.465,0.466,0.467,0.465,0.460,0.452,0.441,0.430,0.418,0.398};
    Float_t minEtErr500MeV[10] = {0.036,0.036,0.036,0.036,0.036,0.037,0.038,0.038,0.038,0.039};
    //10 centrality bins are 0-5,5-10,10-20,20-30,30-40,  40-50,50-60,60-70,70-80,80-90
    int j=0;
    float scale = 1.0;//Scale for eta acceptance and geometric acceptance
    if(det.Contains("Emc")){//Output from the PlotMatchedTracksDeposts.C
      scale = 360.0/40.0/1.2/2;//Azimuthal acceptance over eta range
      for(int i=0;i<10;i++){
	if(i<2){//These bins are exactly what they should bin in the 20 bin binning
	  hadCorrShort[i] = hadCorrEmcal[i];
	  hadErrorShort[i] = hadErrorEmcal[i];
	  hadCorrShort500MeV[i] = hadCorr500MeVEmcal[i];
	  hadErrorShort500MeV[i] = hadError500MeVEmcal[i];
	  j++;//i=0 j=0; i=1 j=1
	}
	else{
	  //i=2 j=2
	  //cout<<"adding bins "<<j<<" and "<<j+1<<endl;
	  hadCorrShort[i] = (hadCorrEmcal[j]+hadCorrEmcal[j+1])/2.0;
	  hadErrorShort[i] = (hadErrorEmcal[j]+hadErrorEmcal[j+1])/2.0;
	  hadCorrShort500MeV[i] = (hadCorr500MeVEmcal[j]+hadCorr500MeVEmcal[j+1])/2.0;
	  hadErrorShort500MeV[i] = (hadError500MeVEmcal[j]+hadError500MeVEmcal[j+1])/2.0;
	  j+=2;
	}
	//The corrections below are already in 10 bin binning
	kaonCorrShort500MeV[i] = kaonCorrEmcal500MeV[i];
	kaonErrorShort500MeV[i] = kaonErrorEmcal500MeV[i];
	kaonCorrShort[i] = kaonCorrEmcal[i];
	kaonErrorShort[i] = kaonErrorEmcal[i];
	minEtCorrShort[i] = minEtCorrEmcal[i];
	minEtErrShort[i] = minEtErrEmcal[i];
	//cout<<"cent bin "<<centLabelsShort[i].Data()<<" hadronic corr "<<hadCorrShort[i]<<endl;
      }
    }
    else{
      scale = 360.0/60/0.24/2;//Azimuthal acceptance over eta range
      for(int i=0;i<10;i++){
	if(i<2){//These bins are exactly what they should bin in the 20 bin binning
	  hadCorrShort[i] = hadCorrPhos[i];
	  hadErrorShort[i] = hadErrorPhos[i];
	  hadCorrShort500MeV[i] = hadCorr500MeVPhos[i];
	  hadErrorShort500MeV[i] = hadError500MeVPhos[i];
	  j++;//i=0 j=0; i=1 j=1
	}
	else{
	  //i=2 j=2
	  hadCorrShort[i] = (hadCorrPhos[j]+hadCorrPhos[j+1])/2.0;
	  hadErrorShort[i] = (hadErrorPhos[j]+hadErrorPhos[j+1])/2.0;
	  hadCorrShort500MeV[i] = (hadCorr500MeVPhos[j]+hadCorr500MeVPhos[j+1])/2.0;
	  hadErrorShort500MeV[i] = (hadError500MeVPhos[j]+hadError500MeVPhos[j+1])/2.0;
	  j+=2;
	}
	//The corrections below are already in 10 bin binning
	kaonCorrShort[i] = kaonCorrPhos[i];
	kaonErrorShort[i] = kaonErrorPhos[i];
	kaonCorrShort500MeV[i] = kaonCorrPhos500MeV[i];
	kaonErrorShort500MeV[i] = kaonErrorPhos500MeV[i];
	minEtCorrShort[i] = minEtCorrPhos[i];
	minEtErrShort[i] = minEtErrPhos[i];
      }
    }
    cout<<"scale:  "<<scale<<endl;

    int j = -1;
    if(det.Contains("Emc")){//Output from the PlotMatchedTracksDeposts.C
      for(int i=0;i<20;i++){
	if(i%2==0){//even numbers
	  j++;//i=0,j=0; i=1,j=0; i=2,j=1
	}
	//cout<<"i "<<i<<" j "<<j<<endl;
	hadCorr[i] = hadCorrEmcal[i];
	hadError[i] = hadErrorEmcal[i];
	kaonCorr[i] = kaonCorrEmcal[j];
	kaonError[i] = kaonErrorEmcal[j];
	minEtCorr[i] = minEtCorrEmcal[j];
	minEtErr[i] = minEtErrEmcal[j];
      }
    }
    else{
      for(int i=0;i<20;i++){
	if(i%2==0){//even numbers
	  j++;//i=0,j=0; i=1,j=0; i=2,j=1
	}
	hadCorr[i] = hadCorrPhos[i];
	hadError[i] = hadErrorPhos[i];
	kaonCorr[i] = kaonCorrPhos[j];
	kaonError[i] = kaonErrorPhos[j];
	minEtCorr[i] = minEtCorrPhos[j];
	minEtErr[i] = minEtErrPhos[j];
      }
    }

    TCanvas *cents = new TCanvas("cents", "cents",700, 600);
    TCanvas *centsAlt = new TCanvas("centsAlt", "centsAlt",700, 600);
    cents->Divide(6,3);
    cents->cd(1);

    int n = 0;
    xpion[n] = 382.8;
    xpionerr[n] = 6;
    n++;
    xpion[n] = 329.7;//n=1
    xpionerr[n] = 6;
    x[n] = 356;
    xerr[n] = 6;
    n++;
    x[n] = xpion[n] = 260.5;//n=2
    xerr[n] = xpionerr[n] = 4.4;
    n++;
    x[n] = xpion[n] = 186.4;//n=3
    xerr[n] = xpionerr[n] = 3.9;
    n++;
    x[n] = xpion[n] = 128.9;//n=4
    xerr[n] = xpionerr[n] = 3.3;
    n++;
    x[n] = xpion[n] = 85;//n=5
    xerr[n] = xpionerr[n] = 2.6;
    n++;
    x[n] = xpion[n] = 52.8;//n=6
    xerr[n] = xpionerr[n] = 2;
    n++;
    x[n] = xpion[n] = 30.0;//n=7
    xerr[n] = xpionerr[n] = 1.3;
    n++;
    x[n] = xpion[n] = 15.8;//n=8
    xerr[n] = xpionerr[n] = 0.6;
    n++;
    x[n] = xpion[n] = 7.48;//n=9
    xerr[n] = xpionerr[n] = 0.29;
    n++;

    for(int i=0;i<10;i++){
      ypion[i] = pionEt[i]/(xpion[i]/2);
      ypionerr[i] = pionEtError[i]/(xpion[i]/2);
      emEtPerNpart[i] = emEt[i]/(xpion[i]/2);
      emEtPerNpartErr[i] = emEtErr[i]/(xpion[i]/2);
    }

    int nmax = 18;

    TFile *f = TFile::Open(filename, "READ");
    if (!f)
    {
        std::cerr << "Could not open file: " << filename << std::endl;
    }

    TList *l = dynamic_cast<TList*>(f->Get("out1"));

    if (!l)
    {
        std::cerr << "Could not get object list from: " << filename << std::endl;
    }
    TString treename = "fEventSummaryTree"+det+"Rec";
    TTree *t = dynamic_cast<TTree*>(l->FindObject(treename.Data()));

    if (!t)
    {
        std::cerr << "Could not get tree from: " << filename << std::endl;
	return;
    }
    TString prefix = "fHistNominal";
    TH2F *fHistNominalRawEt = l->FindObject((prefix+"RawEt"+det+"Rec").Data());
    TH2F *fHistNominalNonLinLowEt = l->FindObject((prefix+"NonLinLowEt"+det+"Rec").Data());
    TH2F *fHistNominalNonLinHighEt = l->FindObject((prefix+"NonLinHighEt"+det+"Rec").Data());
    TH2F *fHistNominalEffLowEt = l->FindObject((prefix+"EffLowEt"+det+"Rec").Data());
    TH2F *fHistNominalEffHighEt = l->FindObject((prefix+"EffHighEt"+det+"Rec").Data());
    TH2F *fHistTotRawEt = l->FindObject("fHistTotRawEt");
    TH2F *fHistTotRawEt500MeV = l->FindObject("fHistTotRawEt500MeV");
    TObjArray rawEt(20);
    TObjArray rawEt500MeV(20);
    TLegend *centLegends[20];
    for(int i = 0; i < nmax; i++)
    {
      //cent bin is nmax-i-1
      int centbin = i;
      //nmax = 18 so this is 17, 16, 15... 0
      rawEt[i]= fHistTotRawEt->ProjectionX(Form("RawEt%i",i),centbin+1,centbin+1);
      rawEt500MeV[i]= fHistTotRawEt500MeV->ProjectionX(Form("RawEt%i",i),centbin+1,centbin+1);
      SetStyles((TH1*)rawEt[i],markers[i],colors[i],"E_{T}","N_{eve}");
      SetStyles((TH1*)rawEt500MeV[i],markers[i],colors[i],"E_{T}","N_{eve}");
      if(((TH1*)rawEt[0])->GetMaximum() < ((TH1*)rawEt[i])->GetMaximum()){
	((TH1*)rawEt[0])->SetMaximum(((TH1*)rawEt[i])->GetMaximum() );
      }
      rawEtValue[centbin] = ((TH1D*)rawEt[i])->GetMean();
      rawEtErr[centbin] = ((TH1D*)rawEt[i])->GetMeanError();
      rawEt500MeVValue[centbin] = ((TH1D*)rawEt[i])->GetMean();
      rawEt500MeVErr[centbin] = ((TH1D*)rawEt[i])->GetMeanError();
      cout<<"Cent bin "<<centbin<<" "<<centLabels[centbin].Data()<<" Raw ET "<<rawEtValue[i]<<" +/- "<<rawEtErr[i]<<" Raw ET> 500 MeV "<<rawEt500MeVValue[i]<<" +/- "<<rawEt500MeVErr[i]<<endl;
      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",centbin+1,centbin+1);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",centbin+1,centbin+1);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",centbin+1,centbin+1);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",centbin+1,centbin+1);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",centbin+1,centbin+1);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
      float efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      nonLinError[centbin] = nonlinfracerr;
      efficiencyError[centbin] = efffracerr;
      cents->cd(i+1);
      cents->cd(i+1)->SetLogy(1);
      ((TH1D*)rawEt[i])->Draw("same");
      cents->cd(i+1)->Update();
      centLegends[i] = new TLegend(0.55,0.7,0.9,0.8);
      
      centLegends[i]->SetFillColor(kWhite);
      centLegends[i]->AddEntry((TH1D*)rawEt[i],"Reconstructed E_{T} distribution", "l");
      centsAlt->cd();
      ((TH1D*)rawEt[i])->Draw("same");
    }
    nmax = 10;
    TObjArray rawEtShort(20);
    TObjArray rawEt500MeVShort(20);
    TLegend *centLegendsShort[20];
    int cbMin = 0;
    int cbMax = 0;
    TCanvas *centsShort = new TCanvas("centsShort", "centsShort",700, 600);
    centsShort->Divide(3,3);
    centsShort->cd(1);
    TCanvas *centsAltShort = new TCanvas("centsAltShort", "centsAltShort",700, 600);
    for(int i = 0; i < nmax; i++)
    {
      if(i>1) cbMax = cbMin+1;//For everything but the first two bins, rebin
      else{cbMax=cbMin;}
      rawEtShort[i]= fHistTotRawEt->ProjectionX(Form("RawEtShort%i",i),cbMin+1,cbMax+1);
      rawEtValueShort[i] = ((TH1D*)rawEtShort[i])->GetMean();
      rawEtErrShort[i] = ((TH1D*)rawEtShort[i])->GetMeanError();
      rawEt500MeVShort[i]= fHistTotRawEt500MeV->ProjectionX(Form("RawEtShort500MeV%i",i),cbMin+1,cbMax+1);
      rawEt500MeVValueShort[i] = ((TH1D*)rawEt500MeVShort[i])->GetMean();
      rawEt500MeVErrShort[i] = ((TH1D*)rawEt500MeVShort[i])->GetMeanError();

      SetStyles((TH1*)rawEtShort[i],markers[i],colors[i],"E_{T}","N_{eve}");
      SetStyles((TH1*)rawEt500MeVShort[i],markers[i],colors[i],"E_{T}","N_{eve}");
      TLine *line = new TLine(rawEtValueShort[i],0,rawEtValueShort[i],((TH1*)rawEtShort[i])->GetMaximum());
      line->SetLineColor(colors[i]);
      if(((TH1*)rawEtShort[0])->GetMaximum() < ((TH1*)rawEtShort[i])->GetMaximum()){
	((TH1*)rawEtShort[0])->SetMaximum(((TH1*)rawEtShort[i])->GetMaximum() );
      }

      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",cbMin+1,cbMax+1);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",cbMin+1,cbMax+1);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",cbMin+1,cbMax+1);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",cbMin+1,cbMax+1);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",cbMin+1,cbMax+1);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
      float efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      nonLinError[nmax-i-1] = nonlinfracerr;
      efficiencyError[nmax-i-1] = efffracerr;
      centsShort->cd(i+1);
      centsShort->cd(i+1)->SetLogy(1);
      ((TH1D*)rawEtShort[i])->Draw("same");
      centsShort->cd(i+1)->Update();
      centsAltShort->cd();
      ((TH1D*)rawEtShort[i])->Draw("same");
      line->Draw();
      centLegendsShort[i] = new TLegend(0.55,0.7,0.9,0.8);
      
      centLegendsShort[i]->SetFillColor(kWhite);
      centLegendsShort[i]->AddEntry((TH1D*)rawEt[i],"Reconstructed E_{T} distribution", "l");
      cbMin = cbMax+1;
      partialCorrEtValueShort[i] = rawEtValueShort[i] - kaonCorrShort[i] - hadCorrShort[i];
      partialCorrEtErrShort[i] = TMath::Sqrt(TMath::Power(kaonErrorShort[i],2)+TMath::Power(hadErrorShort[i],2)+TMath::Power(partialCorrEtValueShort[i]*efffracerr,2)+TMath::Power(partialCorrEtValueShort[i]*nonlinfracerr,2));
      corrEtValueShort[i] = scale*partialCorrEtValueShort[i]/minEtCorr[i];
      corrEtErrShort[i] = corrEtValueShort[i]*TMath::Sqrt(TMath::Power(minEtErr[i]/minEtCorr[i],2)+TMath::Power(partialCorrEtErrShort[i]/partialCorrEtValueShort[i],2));
      corrEtValuePerNPartShort[i] = corrEtValueShort[i]/(xpion[i]/2);
      corrEtErrPerNPartShort[i] = corrEtErrShort[i]/(xpion[i]/2);
      y_rel_error[i] = (corrEtValueShort[i]-emEt[i])/emEt[i];
      cout<<y_rel_error[i]<<" = ("<<corrEtValueShort[i]<<"-"<<emEt[i]<<")/"<<emEt[i]<<endl;
      yerr[i] =  corrEtErrShort[i]/emEt[i];


      partialCorrEtValueShort500MeV[i] = rawEt500MeVValueShort[i] - kaonCorrShort500MeV[i] - hadCorrShort500MeV[i];
      partialCorrEtErrShort500MeV[i] = TMath::Sqrt(TMath::Power(kaonErrorShort500MeV[i],2)+TMath::Power(hadErrorShort500MeV[i],2)+TMath::Power(partialCorrEtValueShort500MeV[i]*efffracerr,2)+TMath::Power(partialCorrEtValueShort500MeV[i]*nonlinfracerr,2));
      corrEtValueShort500MeV[i] = scale*partialCorrEtValueShort500MeV[i]/minEtCorr500MeV[i];
      corrEtErrShort500MeV[i] = corrEtValueShort500MeV[i]*TMath::Sqrt(TMath::Power(minEtErr500MeV[i]/minEtCorr500MeV[i],2)+TMath::Power(partialCorrEtErrShort500MeV[i]/partialCorrEtValueShort500MeV[i],2));
      corrEtValuePerNPartShort500MeV[i] = corrEtValueShort500MeV[i]/(xpion[i]/2);
      corrEtErrPerNPartShort500MeV[i] = corrEtErrShort500MeV[i]/(xpion[i]/2);

      cout<<" partial corr ET "<<partialCorrEtValueShort[i]<<" +/- "<<partialCorrEtErrShort[i]<<" = "<<rawEtValueShort[i]<<"+/-"<<TMath::Power(partialCorrEtValueShort[i]*nonlinfracerr,2)<<"(NL)+/-"<<TMath::Power(partialCorrEtValueShort[i]*efffracerr,2)<<"(eff) - "<<kaonCorrShort[i]<<"+/-"<<kaonErrorShort[i]<<"(kaon) - "<<hadCorrShort[i]<<"+/-"<<hadErrorShort[i]<<"(Had)"<<" frac had ET "<<hadCorrShort[i]/rawEtValueShort[i]<<endl;
      cout<<" corr et "<< corrEtValueShort[i]<<" +/- "<<corrEtErrShort[i]<<" = "<< scale<<"/"<<partialCorrEtValueShort[i]<<"*"<<minEtCorr[i]<<" +/- "<<corrEtValueShort[i]<<"*("<<minEtErr[i]/minEtCorr[i]<<"(Min ET) * "<<partialCorrEtErrShort[i]/partialCorrEtValueShort[i]<<"(others) )"<<endl;
      cout<<"Cent bin "<<centLabelsShort[i].Data()<<" cb min "<<cbMin<<" cb max "<<cbMax<<" Raw ET "<<rawEtValueShort[i]<<" +/- "<<rawEtErrShort[i];//<<endl;
      cout<<" partial corr ET "<<partialCorrEtValueShort[i]<<" +/- "<<partialCorrEtErrShort[i];
      cout<<endl;
      cout<<" corr ET "<<corrEtValueShort[i]<<" +/- "<<corrEtErrShort[i];
      cout<<" pion ET "<<pionEt[i]<<" +/- "<<pionEtError[i]<<" Em Et "<<emEt[i]<<" +/- "<<emEtErr[i]<<endl;
      float correcthadet = rawEtValueShort[i] - pionEt[i]/scale*minEtCorr[i] - kaonCorrShort[i];
      cout<<"Partial corr ET I should measure "<<pionEt[i]/scale*minEtCorr[i]<<" \"correct\" had ET "<<correcthadet<<" correct had et frac "<<correcthadet/rawEtValueShort[i]<<endl;

      cout<<" partial corr ET 500 MeV "<<partialCorrEtValueShort500MeV[i]<<" +/- "<<partialCorrEtErrShort500MeV[i]<<" = "<<rawEt500MeVValueShort[i]<<"+/-"<<TMath::Power(partialCorrEtValueShort500MeV[i]*nonlinfracerr,2)<<"(NL)+/-"<<TMath::Power(partialCorrEtValueShort500MeV[i]*efffracerr,2)<<"(eff) - "<<kaonCorrShort500MeV[i]<<"+/-"<<kaonErrorShort500MeV[i]<<"(kaon) - "<<hadCorrShort500MeV[i]<<"+/-"<<hadErrorShort500MeV[i]<<"(Had)"<<" frac had ET "<<hadCorrShort500MeV[i]/rawEt500MeVValueShort[i]<<endl;
      cout<<" corr et "<< corrEtValueShort500MeV[i]<<" +/- "<<corrEtErrShort500MeV[i]<<" = "<< scale<<"/"<<partialCorrEtValueShort500MeV[i]<<"*"<<minEtCorr500MeV[i]<<" +/- "<<corrEtValueShort500MeV[i]<<"*("<<minEtErr500MeV[i]/minEtCorr500MeV[i]<<"(Min ET) * "<<partialCorrEtErrShort500MeV[i]/partialCorrEtValueShort500MeV[i]<<"(others) )"<<endl;
      cout<<"Cent bin "<<centLabelsShort[i].Data()<<" cb min "<<cbMin<<" cb max "<<cbMax<<" Raw ET "<<rawEt500MeVValueShort[i]<<" +/- "<<rawEt500MeVErrShort[i];//<<endl;
      cout<<" partial corr ET "<<partialCorrEtValueShort500MeV[i]<<" +/- "<<partialCorrEtErrShort500MeV[i];
      cout<<endl;
      cout<<" corr ET "<<corrEtValueShort500MeV[i]<<" +/- "<<corrEtErrShort500MeV[i];
      cout<<endl<<endl;
    }

    TCanvas *c2 = new TCanvas("c2", "dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]",700, 600);
   
    for(int i=0;i<10;i++){
      cout<<"x "<<xpion[i]<<" y "<<corrEtValuePerNPartShort[i]<<" +/- "<<corrEtErrPerNPartShort[i]<<endl;
    }
    TGraphErrors *gr = new TGraphErrors(10,xpion,corrEtValuePerNPartShort,xerr,corrEtErrPerNPartShort);
    gr->SetMarkerStyle(20);

    TLegend *leg = new TLegend(0.363422,0.157061,0.813187,0.357349);
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kWhite);
    leg->AddEntry(gr,"Reconstructed EM E_{T}", "lp");

    TH1F *frame = new TH1F("frame","frame",1,0,2);
    frame->GetYaxis()->SetTitle("dE_{T}/d#eta");
    frame->GetXaxis()->SetTitle("N_{part}");
    //fPion->SetRange(0,2);
    frame->SetMinimum(0);
    frame->SetMaximum(5);
    //frame->Draw();

    TGraphErrors *gr2 = new TGraphErrors(10,xpion,ypion,xpionerr,ypionerr);
    gr2->GetYaxis()->SetTitle("dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]");
    gr2->GetXaxis()->SetTitle("N_{part}");
    gr2->SetTitle("");
    gr2->GetXaxis()->SetRangeUser(0, 400);
    gr2->GetYaxis()->SetRangeUser(0, 2.5);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("AP same");
    leg->AddEntry(gr2,"Reconstructed charged pion E_{T} (scaled by 0.5)", "lp");

    TGraphErrors *gr3 = new TGraphErrors(10,xpion,emEtPerNpart,xpionerr,emEtPerNpartErr);
    gr3->SetMarkerStyle(29);
    gr3->SetMarkerColor(kBlue);
    gr3->Draw("P same");
    leg->AddEntry(gr3,"E_{T}^{em} from charged pions", "lp");
    
    gr->Draw("P same");
    
    leg->Draw();
    
    
    TString title;
    TString ytitle;
    TString xtitle = "N_{part}";
    if(sim)
    {
	title = "Relative Discrepancy for reconstructed vs true";
	ytitle = "#frac{E_{T,rec}-E_{T,true}}{E_{T,true}}";
    }
    else
    {
      title = "Relative Discrepancy for EM E_{T} vs EM #pi^{+/-}";
      ytitle = "#frac{E_{T,rec}-E_{T,#pi^{+/-}}}{E_{T,#pi^{+/-}}}";
    }
     
    TCanvas *c3 = new TCanvas("c3", title,700, 600);
    TGraphErrors *gr4 = 0;
    if(sim)
    {
      gr4 = new TGraphErrors(10,xpion,y_rel_error,xpionerr,yerr);
    }
    else
    {
      gr4 = new TGraphErrors(10,xpion,y_rel_error,xpionerr,yerr);
    }
    gr4->SetTitle(title);
    gr4->GetYaxis()->SetTitle(ytitle);
    gr4->GetXaxis()->SetTitle(xtitle);
    gr4->GetXaxis()->SetRangeUser(0, 400);
    if(sim) gr4->GetYaxis()->SetRangeUser(-0.2,0.2);
    else gr4->GetYaxis()->SetRangeUser(-1,0);
    gr4->SetMarkerStyle(20);
    gr4->Draw("AP");
    
    TString name = "dn_deta_npart_scale_factor_" + TString::Format("%.2f", scaleFactor);
    if(sim) name += "_sim.png";
    else name += "_real.png";
    c2->Print(name, "png");

    name = "rel_error_";
    if(sim) name += "_sim.png";
    else name += "_real.png";
    c3->Print(name, "png");

    if(sim) name = "tot_neutral_et_sim.png";
    else name = "tot_neutral_et_real.png";
    cents->Print(name, "png");



}
