#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"

#include "TRandom.h"
#include "TPad.h"
#include "TCanvas.h"

/// \file LandauTest.C

class TLandauMean: public TObject {
public:
  void Init(Int_t n, Float_t mean, Float_t sigma);  // initial parameters
  void Gener();          // gener sample 

  //void Anal();

  Int_t fNSample;      // number of samples
  Float_t fLMean;        // landau mean
  Float_t fLSigma;       // landau sigma
  //
  Float_t fTM_0_6[3];    // truncated method  - first 3 momenta
  Float_t fTM_0_7[3];    // truncated method  - first 3 momenta
  Float_t fTM_0_8[3];    // truncated method  - first 3 momenta
  Float_t fTM_0_10[3];   // truncated method  - first 3 momenta
  //
  Float_t fLM_0_6[3];    // truncated log.  method  - first 3 momenta
  Float_t fLM_0_7[3];    // truncated log.  method  - first 3 momenta
  Float_t fLM_0_8[3];    // truncated log.  method  - first 3 momenta
  Float_t fLM_0_10[3];   // truncated log.  method  - first 3 momenta

  Float_t fMedian3;      // median 3 value
private:
  Float_t Moment3(Float_t sum1, Float_t sum2, Float_t sum3, Int_t n, Float_t m[3]);
  ClassDef(TLandauMean,1)
};

/// \cond CLASSIMP
ClassImp(TLandauMean)
/// \endcond

void TLandauMean::Init(Int_t n, Float_t mean, Float_t sigma)
{
  /// init parameters

  fNSample = n;
  fLMean   = mean;
  fLSigma  = sigma;
}

Float_t TLandauMean::Moment3(Float_t sumi1, Float_t sumi2, Float_t sumi3, Int_t sum, Float_t m[3])
{
  Float_t m3=0;

  //  m3 = (sumi3-3*pos*sumi2+3*pos*pos*sumi-pos*pos*pos*sum)/sum; 
  Float_t pos = sumi1/sum;
  m[0] = pos;
  m[1] = sumi2/sum-pos*pos;
  if (m[1]<=0){
    printf("pici pici\n");
  }
  else
    m[1] = TMath::Sqrt(m[1]);
  m3 = (sumi3-3*pos*sumi2+3*pos*pos*sumi1-pos*pos*pos*sum)/sum; 
  Float_t sign = m3/TMath::Abs(m3);
  m3 = TMath::Power(sign*m3,1/3.);
  m3*=sign;

  m[2] = m3;  
  return m3;
}

void TLandauMean::Gener()
{
  /// generate sample

  Float_t * buffer = new Float_t[fNSample];
  
  for (Int_t i=0;i<fNSample;i++) {
    buffer[i] = gRandom->Landau(fLMean,fLSigma); 
    if (buffer[i]>1000) buffer[i]=1000;
  }

  Int_t *index = new Int_t[fNSample];
  TMath::Sort(fNSample,buffer,index,kFALSE);

  //
  Float_t median = buffer[index[fNSample/3]];
  //
  Float_t sum06[4]  = {0.,0.,0.,0.};
  Float_t sum07[4]  = {0.,0.,0.,0.};
  Float_t sum08[4]  = {0.,0.,0.,0.};
  Float_t sum010[4]  = {0.,0.,0.,0.};
  //
  Float_t suml06[4] = {0.,0.,0.,0.};
  Float_t suml07[4] = {0.,0.,0.,0.};
  Float_t suml08[4] = {0.,0.,0.,0.};
  Float_t suml010[4] = {0.,0.,0.,0.};
  //

  for (Int_t i =0; i<fNSample; i++){
    Float_t amp  = buffer[index[i]];
    Float_t lamp = median*TMath::Log(1.+amp/median);
    if (i<0.6*fNSample){
      sum06[0]+= amp;
      sum06[1]+= amp*amp;
      sum06[2]+= amp*amp*amp;
      sum06[3]++;
      suml06[0]+= lamp;
      suml06[1]+= lamp*lamp;
      suml06[2]+= lamp*lamp*lamp;
      suml06[3]++;
    }

    if (i<0.7*fNSample){
      sum07[0]+= amp;
      sum07[1]+= amp*amp;
      sum07[2]+= amp*amp*amp;
      sum07[3]++;
      suml07[0]+= lamp;
      suml07[1]+= lamp*lamp;
      suml07[2]+= lamp*lamp*lamp;
      suml07[3]++;
    }
    if (i<0.8*fNSample){
      sum08[0]+= amp;
      sum08[1]+= amp*amp;
      sum08[2]+= amp*amp*amp;
      sum08[3]++;
      suml08[0]+= lamp;
      suml08[1]+= lamp*lamp;
      suml08[2]+= lamp*lamp*lamp;
      suml08[3]++;
    }
    if (i<1*fNSample){
      sum010[0]+= amp;
      sum010[1]+= amp*amp;
      sum010[2]+= amp*amp*amp;
      sum010[3]++;
      suml010[0]+= lamp;
      suml010[1]+= lamp*lamp;
      suml010[2]+= lamp*lamp*lamp;
      suml010[3]++;

    }
  }
  //  
  fMedian3 = median;
  //
  Moment3(sum06[0],sum06[1],sum06[2],sum06[3],fTM_0_6);  
  Moment3(sum07[0],sum07[1],sum07[2],sum07[3],fTM_0_7);  
  Moment3(sum08[0],sum08[1],sum08[2],sum08[3],fTM_0_8);  
  Moment3(sum010[0],sum010[1],sum010[2],sum010[3],fTM_0_10);  
  //

  Moment3(suml06[0],suml06[1],suml06[2],suml06[3],fLM_0_6);  
  Moment3(suml07[0],suml07[1],suml07[2],suml07[3],fLM_0_7);  
  Moment3(suml08[0],suml08[1],suml08[2],suml08[3],fLM_0_8);  
  Moment3(suml010[0],suml010[1],suml010[2],suml010[3],fLM_0_10);  
  //
  fLM_0_6[0] = (TMath::Exp(fLM_0_6[0]/median)-1.)*median;
  fLM_0_7[0] = (TMath::Exp(fLM_0_7[0]/median)-1.)*median;
  fLM_0_8[0] = (TMath::Exp(fLM_0_8[0]/median)-1.)*median;
  fLM_0_10[0] = (TMath::Exp(fLM_0_10[0]/median)-1.)*median;
  //
  delete [] buffer;
}   


void GenerLandau(Int_t nsamples)
{
  TLandauMean * landau = new TLandauMean;
  TFile f("Landau.root","recreate");
  TTree * tree = new TTree("Landau","Landau");
  tree->Branch("Landau","TLandauMean",&landau); 
  
  for (Int_t i=0;i<nsamples;i++){
    Int_t   n     = 20 + Int_t(gRandom->Rndm()*150);
    Float_t mean  = 40. +gRandom->Rndm()*50.;
    Float_t sigma = 5.  +gRandom->Rndm()*15.;
    landau->Init(n, mean, sigma);
    landau->Gener();
    tree->Fill();
  }
  tree->Write();
  f.Close();

}





TH1F *  LandauTest(Float_t meano,  Float_t sigma, Float_t meanlog0, Int_t n,Float_t ratio)
{ 
  /// test for different approach of de dx resolution
  /// meano, sigma  - mean value of Landau distribution and sigma
  /// meanlog0      - scaling factor for logarithmic mean value
  /// n             - number of used layers
  /// ratio         - ratio of used amplitudes for truncated mean


  TCanvas * pad = new TCanvas("Landau test");
  pad->Divide(2,2);
  TH1F * h1 = new TH1F("h1","Logarithmic mean",300,0,4*meano);
  TH1F * h2 = new TH1F("h2","Logarithmic amplitudes",300,0,8*meano);
  TH1F * h3 = new TH1F("h3","Mean",300,0,4*meano);
  TH1F * h4 = new TH1F("h4","Amplitudes",300,0,8*meano);

  for(Int_t j=0;j<10000;j++){
    //generate sample and sort it
    Float_t * buffer = new Float_t[n];
    Float_t * buffer2= new Float_t[n];
    
    for (Int_t i=0;i<n;i++) {
      buffer[i] = gRandom->Landau(meano,sigma); 
      buffer2[i] = buffer[i];    
    }
    //add crosstalk
    for (Int_t i=1;i<n-1;i++) {
      buffer[i] =    buffer2[i]*1.0+ buffer2[i-1]*0.0+ buffer2[i+1]*0.0;
      buffer[i] = TMath::Min(buffer[i],1000.); 
    }
    Int_t *index = new Int_t[n];
    TMath::Sort(n,buffer,index,kFALSE);

    //calculate mean
    Float_t sum;
    sum=0;
    Float_t mean;
    Float_t used = 0;
    for (Int_t i=0;i<n*ratio;i++) {
      if (buffer[index[i]]<1000.){
	Float_t amp = meanlog0*TMath::Log(1+buffer[index[i]]/meanlog0);
	sum += amp;            
	used++;
      }
    }
    mean = sum/used;
    //
    sum=0;
    used=0;
    Float_t sum2=0;
    Float_t meanlog =meanlog0;
    for (Int_t i=0;i<n*ratio;i++) {
      if (buffer[index[i]]<1000.){
	Float_t amp = meanlog*TMath::Log(1.+buffer[index[i]]/(meanlog));
	sum +=amp;
	sum2+=buffer[index[i]];
	used++;
	h2->Fill(amp);
	h4->Fill(buffer[index[i]]);
      }
    }
    mean = sum/used;
    mean = (TMath::Exp(mean/meanlog)-1)*meanlog;
    Float_t mean2 = sum2/used;
    //mean2 = (mean+mean2)/2.;
    h1->Fill(mean);    
    h3->Fill(mean2);
  }

  pad->cd(1);
  h1->Draw();
  pad->cd(2);
  h2->Draw();
  pad->cd(3);
  h3->Draw();
  pad->cd(4);
  h4->Draw();


  return h1;

}

