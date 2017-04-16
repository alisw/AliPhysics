/*
 Small macro for testing the raw response function used in
 AliEMCALRawUtils

 J.L. Klay (Cal Poly)

*/
const Double_t fgTimeTrigger = 1.5E-06;
const Double_t fTau = 2.35;
const Double_t fOrder = 2.;
const Double_t fgPedestal = 32;
const Double_t fTimeBinWidth = 100E-9;
const Int_t fgOverflow = 0x3FF;
const Double_t fHighLowGain = 16.;
const Double_t fRawFormatTimeMax = 256*fTimeBinWidth;
const Double_t fgFEENoise = 3.;


void TestRawResponse(const Double_t damp = 20, const Double_t dtime = 1e-09) {

  TH1I* adcHigh = new TH1I("adcHigh","adcHigh",256,0,255);
  TH1I* adcLow = new TH1I("adcLow","adcLow",256,0,255);

  TF1* signalF = new TF1("signal",RawResponseFunction, 0, 256, 5);
  signalF->SetParameter(0,damp);
  signalF->SetParameter(1,(dtime+fgTimeTrigger)/fTimeBinWidth);
  signalF->SetParameter(2,fTau);
  signalF->SetParameter(3,fOrder);
  signalF->SetParameter(4,fgPedestal); 

  for(Int_t itime = 0; itime < 256; itime++) {
    Double_t signal = signalF->Eval(itime);

    Double_t noise = gRandom->Gaus(0.,fgFEENoise);
    signal = sqrt(signal*signal + noise*noise);

    cout << "itime = " << itime << " highgain = " << signal << " signalI = " << static_cast<Int_t>(signal +0.5);

    if(static_cast<Int_t>(signal +0.5) > fgOverflow) {
      adcHigh->SetBinContent(itime+1,fgOverflow);
    } else {
      adcHigh->SetBinContent(itime+1,static_cast<Int_t>(signal +0.5));
    }

    signal /= fHighLowGain;

    cout << " lowgain = " << signal << " signalI = " << static_cast<Int_t>(signal +0.5) << endl;

    if(static_cast<Int_t>(signal +0.5) > fgOverflow) {
      adcLow->SetBinContent(itime+1,fgOverflow);
    } else {
      adcLow->SetBinContent(itime+1,static_cast<Int_t>(signal +0.5));      
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",20,20,600,1000);
  c1->Divide(1,3);
  c1->cd(1);
  signalF->Draw();
  c1->cd(2);
  adcHigh->Draw();
  c1->cd(3);
  adcLow->Draw();
  
}

Double_t RawResponseFunction(Double_t *x, Double_t *par)
{
  Double_t signal ;
  Double_t tau =par[2];
  Double_t N =par[3];
  Double_t ped = par[4];
  Double_t xx = ( x[0] - par[1] + tau ) / tau ;

  if (xx <= 0)
    signal = ped ;
  else {
    signal = ped + par[0] * TMath::Power(xx , N) * TMath::Exp(N * (1 - xx )) ;
  }
  return signal ;
}
