#include "AliForwardQCumulantRun2.h"
#include "TMath.h"
#include <iostream>
#include "TRandom.h"
#include "AliForwardFlowRun2Settings.h"
#include "TH1D.h"
#include <complex>
#include <cmath>
using namespace std;


//_____________________________________________________________________
AliForwardQCumulantRun2::AliForwardQCumulantRun2():
fMaxMoment(5),
nHBins(fMaxMoment*4+2),
nRefBins(1),
nEtaBins(24),
fCumuRef("fCumuRef","fCumuRef", nRefBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge, nHBins, 0.5, nHBins+0.5),
fCumuDiff("fCumuDiff","fCumuDiff", nEtaBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge, nHBins, 0.5, nHBins+0.5),
fSumOfWeights(0),
fSumOfWeightsSquared(0),
useEvent(kTRUE)
{
}


//_____________________________________________________________________
void AliForwardQCumulantRun2::CumulantsAccumulate(TH2D& dNdetadphi, TList* outputList, double cent, double vertexpos,UInt_t r, TString detType)
{
  TList* eventList = static_cast<TList*>(outputList->FindObject("EventInfo"));
  TH2F* fOutliers = static_cast<TH2F*>(eventList->FindObject("hOutliers"));
  TH1D* fFMDHits = static_cast<TH1D*>(eventList->FindObject("FMDHits"));

  Int_t nBadBins = 0;
  //Double_t limit = 9999.;
  Int_t phibins = dNdetadphi.GetNbinsY();
    Double_t sumOfWeights = 0;

  for (Int_t etaBin = 1; etaBin <= dNdetadphi.GetNbinsX(); etaBin++) {
  
    Double_t acceptance = 1.;

    Double_t eta = dNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    Double_t runAvg = 0;
    Double_t avgSqr = 0;
    Double_t max = 0;
    Int_t nInAvg = 0;

    // Check for acceptance
    //if ( fabs(eta) > 1.7) {
    //  if (dNdetadphi.GetBinContent(etaBin, 0) == 0) break;
    //}
    /*
    if (fabs(eta) > 1.7 && detType == "forward"){
      acceptance = dNdetadphi.GetBinContent(etaBin, kphiAcceptanceBin);
      
      if (acceptance == 0 || acceptance > 2.0) continue;
    }*/

    for (Int_t phiBin = 1; phiBin <= phibins; phiBin++) {
  

      Double_t phi = dNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      Double_t refEtaBin = fCumuRef.GetXaxis()->FindBin(eta);
      Double_t weight = dNdetadphi.GetBinContent(etaBin, phiBin);
      Double_t refEta = fCumuRef.GetXaxis()->GetBinCenter(refEtaBin);

      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) max = weight;
      
      // Fill into Cos() and Sin() hists
      fFMDHits->Fill(weight);
      if ((fSettings.fFlowFlags & fSettings.kEtaGap) && (fabs(eta) > 0.4 && fabs(eta) < 0.8)) { 
        fCumuRef.Fill(refEta, 0., weight);// mult goes in underflowbin kWA2
        //fCumuRef.Fill(refEta, 1., weight*weight);// mult goes in underflowbin kWA2
      }
      //else if (!(fSettings.fFlowFlags & fSettings.kEtaGap)) fCumuRef.Fill(refEta, 0., weight);// mult goes in underflowbin
      else if (detType == "central") fCumuRef.Fill(refEta, 0., weight);// mult goes in underflowbin

      //if (fabs(eta) >= 1.7 && detType == "forward"){
      //  fCumuDiff.Fill(eta, 0., weight);
        //fCumuDiff.Fill(eta, 1., weight*weight);
      //}
      //else if (fabs(eta) < 1.7 && detType == "central"){
        fCumuDiff.Fill(eta, 0., weight);
        //fCumuDiff.Fill(eta, 1., weight*weight);
      //}
      
      for (Int_t n = 1; n <= 2*fMaxMoment; n++) {         
        Double_t cosBin = fCumuDiff.GetYaxis()->GetBinCenter(GetBinNumberCos(n));         
        Double_t sinBin = fCumuDiff.GetYaxis()->GetBinCenter(GetBinNumberSin(n));   
        Double_t cosnPhi = weight*TMath::Cos(n*phi);         
        Double_t sinnPhi = weight*TMath::Sin(n*phi);        
        sumOfWeights += weight;

        fSumOfWeights += weight;
        fSumOfWeightsSquared += weight*weight;

        // fill diff         
        //if (fabs(eta) >= 1.7 && detType == "forward"){
        //  fCumuDiff.Fill(eta, cosBin, cosnPhi); 
        //  fCumuDiff.Fill(eta, sinBin, sinnPhi);
        //}
        //else if (fabs(eta) < 1.7 && detType == "central"){
        //  fCumuDiff.Fill(eta, cosBin, cosnPhi);         
          fCumuDiff.Fill(eta, sinBin, sinnPhi); 
        //}

        // fill ref         
        if ((fSettings.fFlowFlags & fSettings.kEtaGap) &&  ((fabs(eta) > 0.4) && (fabs(eta) < 0.8) ) ){
            fCumuRef.Fill(refEta, cosBin, cosnPhi);         
            fCumuRef.Fill(refEta, sinBin, sinnPhi);
        }
        else if (detType == "central") {
        //else if (!(fSettings.fFlowFlags & fSettings.kEtaGap) && detType == "central") {
          fCumuRef.Fill(refEta, cosBin, cosnPhi);         
          fCumuRef.Fill(refEta, sinBin, sinnPhi);
        }
      } // End of n loop
    } // End of phi loop
    // Outlier cut calculations
    double fSigmaCut = 4.0;
    if (nInAvg > 0) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = (nInAvg > 1 ? TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg) : 0);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent < 60) nBadBins++;
      else nBadBins = 0;
      fOutliers->Fill(cent, nSigma);
      // We still finish the loop, for fOutliers to make sense, 
      // but we do no keep the event for analysis 
      if (nBadBins > 3) useEvent = kFALSE;
    }
  } // End of eta bin


  if (sumOfWeights < 10) useEvent = kFALSE;
  if (!useEvent) std::cout << "BAD EVENT" << std::endl;
  return;
}

void AliForwardQCumulantRun2::saveEvent(TH2D& dNdetadphi, TList* outputList, double cent, double vertexpos,UInt_t r, TString detType){
  TList* stdQCList = static_cast<TList*>(outputList->FindObject("StdQC"));
  TList* refList = static_cast<TList*>(stdQCList->FindObject("Reference"));
  TList* difList = static_cast<TList*>(stdQCList->FindObject("Differential"));

  THnD* cumuRef = 0;
  THnD* cumuDiff = 0;

  // REFERENCE FLOW
  // For each n we loop over the hists
    Double_t noSamples = static_cast<Double_t>(r);
    for (Int_t n = 2; n <= fMaxMoment; n++) {
      Int_t prevRefEtaBin = 0;

      cumuRef = static_cast<THnD*>(refList->FindObject(Form("cumuRef_v%d", n)));
      cumuDiff = static_cast<THnD*>(difList->FindObject(Form("cumuDiff_v%d", n)));

      // Per mom. quantities
      Double_t dQnReA = 0, dQnImA = 0, multA = 0; 
      Double_t dQ2nReA = 0, dQ2nImA = 0;
      Double_t two = 0, w2 = 0, four = 0, w4 = 0;
      Double_t dQnReB = 0, dQnImB = 0, multB = 0; 
      //Double_t dQ2nReB = 0, dQ2nImB = 0;

      for (Int_t etaBin = 1; etaBin <= fCumuDiff.GetNbinsX(); etaBin++) {
        Double_t eta = fCumuDiff.GetXaxis()->GetBinCenter(etaBin);
        std::cout << "eta = " << eta << std::endl;

        Double_t refEtaBinA = fCumuRef.GetXaxis()->FindBin(eta);
        Double_t refEtaA = fCumuRef.GetXaxis()->GetBinCenter(refEtaBinA);
        Double_t refEtaBinB = fCumuRef.GetXaxis()->FindBin(eta);

        if ((fSettings.fFlowFlags & fSettings.kEtaGap)) refEtaBinB = fCumuRef.GetXaxis()->FindBin(-eta);
        //if (kEtaGap && prevRefEtaBin > 0)prevRefEtaBin = refEtaBinA;

        //Double_t refEtaB = fCumuRef.GetXaxis()->GetBinCenter(refEtaBinB);
        //Double_t sumOfWeights  = fCumuRef.GetBinContent(refEtaBinA, 0);

        //complex<double> Q (dQnReA, dQnImA);

        multA   = fCumuRef.GetBinContent(refEtaBinA, 0);
        dQnReA  = fCumuRef.GetBinContent(refEtaBinA, GetBinNumberCos(n));
        dQnImA  = fCumuRef.GetBinContent(refEtaBinA, GetBinNumberSin(n));
        dQ2nReA = fCumuRef.GetBinContent(refEtaBinA, GetBinNumberCos(2*n));
        dQ2nImA = fCumuRef.GetBinContent(refEtaBinA, GetBinNumberSin(2*n));

        multB   = fCumuRef.GetBinContent(refEtaBinB, 0);
        dQnReB  = fCumuRef.GetBinContent(refEtaBinB, GetBinNumberCos(n));
        dQnImB  = fCumuRef.GetBinContent(refEtaBinB, GetBinNumberSin(n));
        //dQ2nReB = fCumuRef.GetBinContent(refEtaBinB, GetBinNumberCos(2*n));
        //dQ2nImB = fCumuRef.GetBinContent(refEtaBinB, GetBinNumberSin(2*n));

        if (refEtaBinA != prevRefEtaBin) {
            if (multA <= 3) continue; 

        // The reference flow is calculated 
        // 2-particle
            if ((fSettings.fFlowFlags & fSettings.kEtaGap)){
              w2 = multA * multB;
              two = (dQnReA*dQnReB + dQnImA*dQnImB);   
            }
            else{
              w2 = multA * (multA - 1.);
              two = (dQnReA*dQnReA + dQnImA*dQnImA - multA);
            }

            //if (w2 <= 0) continue;
            //if (two/w2 < 0) continue;
            
            Double_t x[5] = {noSamples, vertexpos, refEtaA, cent, fSettings.kW4Four};

            x[4] = fSettings.kW2Two;
            cumuRef->Fill(x, two);
            x[4] = fSettings.kW2;
            cumuRef->Fill(x, w2);
            x[4] = fSettings.kWA;
            cumuRef->Fill(x, multA);
            x[4] = fSettings.kWB;
            cumuRef->Fill(x, multB);


            // NUA
            x[4] = fSettings.kCosphi1A;
            cumuRef->Fill(x, dQnReA);
            x[4] = fSettings.kSinphi1A;
            cumuRef->Fill(x, dQnImA);
            x[4] = fSettings.kCosphi1B;
            cumuRef->Fill(x, dQnReB);
            x[4] = fSettings.kSinphi1B;
            cumuRef->Fill(x, dQnImB);

            if (!(fSettings.fFlowFlags & fSettings.kEtaGap)) {
              //Double_t w3 = (mp*multA-2.*mq)*(multA-1.);

              w4 = multA * (multA - 1.) * (multA - 2.) * (multA - 3.);
              four = 
              2.*multA*(multA-3.) + TMath::Power((TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.)),2.)
              -4.*(multA-2.)*(TMath::Power(dQnReA,2.) + TMath::Power(dQnImA,2.))
              -2.*(TMath::Power(dQnReA,2.)*dQ2nReA+2.*dQnReA*dQnImA*dQ2nImA-TMath::Power(dQnImA,2.)*dQ2nReA)
              +(TMath::Power(dQ2nReA,2.)+TMath::Power(dQ2nImA,2.));


              x[4] = fSettings.kW4Four;
              cumuRef->Fill(x, four);
              x[4] = fSettings.kW4;
              cumuRef->Fill(x, w4);

              x[4] = fSettings.k3pWeight;
              cumuRef->Fill(x, multA*(multA-1.)*(multA-2.));

            // NUA
            x[4] = fSettings.kCosphi1phi2p;
            cumuRef->Fill(x, dQnReA*dQnReA - dQ2nReA);
            x[4] = fSettings.kCosphi1phi2m;
            cumuRef->Fill(x, two);
            x[4] = fSettings.kSinphi1phi2p;
            cumuRef->Fill(x, dQnImA*dQnImA - dQ2nImA);
            x[4] = fSettings.kCosphi1phi2phi3m;
            cumuRef->Fill(x, calcCosPhi1Phi2Phi3m(dQnReA, dQnImA, dQ2nReA, dQ2nImA, multA) );
            x[4] = fSettings.kSinphi1phi2phi3m;
            cumuRef->Fill(x, calcSinPhi1Phi2Phi3m(dQnReA, dQnImA, dQ2nReA, dQ2nImA, multA) );
            }
            prevRefEtaBin = refEtaBinA;
          
        }

    // DIFFERENTIAL FLOW -----------------------------------------------------------------------------

    // For each etaBin bin the necessary values for differential flow is calculated

          Double_t mp = fCumuDiff.GetBinContent(etaBin, 0);
        Double_t mq = 0;
        Double_t pnRe  = fCumuDiff.GetBinContent(etaBin, GetBinNumberCos(n));
        Double_t pnIm  = fCumuDiff.GetBinContent(etaBin, GetBinNumberSin(n));
        //Double_t p2nRe = fCumuDiff.GetBinContent(etaBin, GetBinNumberCos(2*n));
        //Double_t p2nIm = fCumuDiff.GetBinContent(etaBin, GetBinNumberSin(2*n));
        Double_t qnRe  = 0;
        Double_t qnIm  = 0;
        Double_t q2nRe = 0;
        Double_t q2nIm = 0;
        if (!(fSettings.fFlowFlags & fSettings.kEtaGap) && detType == "central"){
          //Double_t qnRe  = pnRe;
          //Double_t qnIm  = pnIm;
          //Double_t q2nRe = p2nRe;
          //Double_t q2nIm = p2nIm;
          mq = mp;
        }


        if (mp == 0) continue;
        if (multB <= 3) continue; 

        Double_t twoPrime = 0.;

        Double_t w2p = mp * multB - mq;
        if (w2p == 0) continue;
          twoPrime = (pnRe*dQnReB + pnIm*dQnImB - mq);

        //if (pnRe == 0) continue;
        if (twoPrime < 0){
           std::cout << "two = " << two << std::endl;
           std::cout << "multB = " << multB << std::endl;
           std::cout << "mp = " << mp << std::endl;
           std::cout << "pnRe = " << pnRe << std::endl;
           std::cout << "dQnReB = " << dQnReB << std::endl;
                      std::cout << "pnIm = " << pnIm << std::endl;
           std::cout << "dQnImB = " << dQnImB << std::endl;
}
        Double_t x[5] = {noSamples,vertexpos, eta, cent, fSettings.kW2Two};
        cumuDiff->Fill(x, twoPrime);

        x[4] = fSettings.kW2;
        cumuDiff->Fill(x, w2p);

        x[4] = fSettings.kWA;
        cumuDiff->Fill(x, mp);

        if ((fSettings.fFlowFlags & fSettings.kEtaGap)) continue;

        Double_t w4p = (mp * multA - 3.*mq)*(multA - 1.)*(multA - 2.);
        Double_t fourPrime = calcFourPrime(dQnReA, dQnImA, qnRe, qnIm, pnRe, pnIm, q2nRe, q2nIm, dQ2nReA, dQ2nImA, mq, multA); 

        x[4] = fSettings.k3pWeight;
        cumuDiff->Fill(x, (mp*multA-2.*mq)*(multA-1.));       
        x[4] = fSettings.kW4Four;
        cumuDiff->Fill(x, fourPrime);
        x[4] = fSettings.kW4;
        cumuDiff->Fill(x, w4p);

        // NUA
        x[4] = fSettings.kSinphi1A;
        cumuDiff->Fill(x,pnIm);
        x[4] = fSettings.kCosphi1A;
        cumuDiff->Fill(x,pnRe);     
        x[4] = fSettings.kCosphi1phi2m;
        cumuDiff->Fill(x, twoPrime);
        x[4] = fSettings.kCosphi1phi2p;
        cumuDiff->Fill(x, pnRe*dQnReB - q2nRe);
        x[4] = fSettings.kSinphi1phi2p;
        cumuDiff->Fill(x, pnIm*dQnImB - q2nIm);
        x[4] = fSettings.kCosphi1phi2phi3m;
        cumuDiff->Fill(x, calcCosPsi1Phi2Phi3m(pnRe, dQnReA, dQnImA, pnIm, dQ2nReA, dQ2nImA, mq, qnRe));
        x[4] = fSettings.kSinphi1phi2phi3m;
        cumuDiff->Fill(x, calcSinPsi1Phi2Phi3m(pnIm, dQnReA, dQnImA, pnRe, dQ2nReA, dQ2nImA, mq, qnIm));
        x[4] = fSettings.kCosphi1phi2phi3p;
        cumuDiff->Fill(x, calcCosPsi1Phi2Phi3p(pnRe, dQnImA, dQnReA, multA, q2nRe, q2nIm, mq, qnRe));
        x[4] = fSettings.kSinphi1phi2phi3p;
        cumuDiff->Fill(x, calcSinPsi1Phi2Phi3p(pnIm, dQnImA, dQnReA, multA, mq, qnIm, q2nIm, q2nRe));


      }
    } // End of moment loop
  


  return;
}


//_____________________________________________________________________
Int_t AliForwardQCumulantRun2::GetBinNumberCos(Int_t n) const  
{
  //
  //  Get the bin number of <<cos(nphi)>>
  //
  //  Parameters:
  //   n: moment
  //
  //  Return: bin number
  //

  return TMath::Abs(n*2+1);
}

//_____________________________________________________________________
Int_t AliForwardQCumulantRun2::GetBinNumberSin(Int_t n) const  
{
  //
  //  Get the bin number of <<sin(nphi)>>
  //
  //  Parameters:
  //   n: moment
  //
  //  Return: bin number
  //

  return TMath::Abs(n*2);
}


void AliForwardQCumulantRun2::reset() {
  fCumuRef.Reset();
  fCumuDiff.Reset();
  fCumuNUARef.Reset();
  fCumuNUADiff.Reset();
}
