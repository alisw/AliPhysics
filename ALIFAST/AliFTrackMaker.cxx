//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast TrackMaker class.                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------//
//                                                                      //
// origin: "res.f" fortran by Karel Safarik which was used to           //
//         calculate the track resolution for TP.                       //
//         Different detectors and material can be selected.            //
//         The basic routines compute information and error matrices    //
//         used for the calculation of momentum resolution.             //
//         see references: ASK KAREL??                                  //
//                                                                      //
// C++ in AliFast framework: Elzbieta Richter-Was and Yiota Foka        //
//                           following general structure od Makers in   //
//                           ATLFast by R. Brun.                        //
//                                                                      //
// purpose: provide a Maker which by using general basic routines of    //
//          "res.f" computes the necessary elements of covariance matrix// 
//          for the calculation of Track Resolution.                    //
//          Those elements are the product of the TrackResolMaker and   //
//          are hold in TrackResol class. They are expected to be used  //
//          together with additional information for the calculation of //
//          the smeared momenta.                                        //
//          Additional information necessary for this calculation       //
//          will be provided via classes or functions specific to the   //
//          specific study and/or detectors.                            //
//          One can select the detector and/or material for a specific  //
//          study.                                                      //
//                                                                      //
// starting point: res.f will be initialy partially contained in        // 
//                 AliFTrackResolMaker and in AliFDet                   // 
//                 It will be reorganised further in the framework of   //
//                 AliFast according to the needs.                      //
//                 Names of variables are kept as in fortran code.      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifdef WIN32
// there is a bug in the Microsoft VisualC++ compiler
// this class must be compiled with optimization off on Windows
# pragma optimize( "", off )
#endif

#include <TParticle.h>
#include <TFile.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "AliFast.h"
//#include "AliFMCMaker.h"
#include "AliFTrackMaker.h"
#include "AliFTrack.h"
#include "AliFDet.h"

const Double_t kPi       = TMath::Pi();
const Double_t k2Pi      = 2*kPi;
const Double_t kPiHalf   = kPi/2.;
extern  AliFast * gAliFast;
ClassImp(AliFTrackMaker)

//_____________________________________________________________________________
AliFTrackMaker::AliFTrackMaker()
{
   fNTracks = 0;
}

//_____________________________________________________________________________
AliFTrackMaker::AliFTrackMaker(const char *name, const char *title)
                 :AliFMaker(name,title)
{
//    Default Setters for tracks

   fFruits     = new TClonesArray("AliFTrack",100, kFALSE);
   fBranchName = "Tracks";
   fNTracks    = 0;
// Please, how to do this optionally ??!!!
   Save();
}

//_____________________________________________________________________________
AliFTrackMaker::~AliFTrackMaker()
{
   //dummy
}

//_____________________________________________________________________________
AliFTrack *AliFTrackMaker::AddTrack(Int_t code, Double_t charge, 
                                    Double_t pT, Double_t eta,Double_t phi,
                                    Double_t v11, Double_t v22, Double_t v33,
                                    Double_t v12, Double_t v13, Double_t v23, Int_t iFlag)
{
//            Add a new track to the list of tracks

 //Note the use of the "new with placement" to create a new track object.
 //This complex "new" works in the following way:
 //   tracks[i] is the value of the pointer for track number i in the TClonesArray
 //   if it is zero, then a new track must be generated. This typically
 //   will happen only at the first events
 //   If it is not zero, then the already existing object is overwritten
 //   by the new track parameters.
 // This technique should save a huge amount of time otherwise spent
 // in the operators new and delete.

   TClonesArray &tracks = *(TClonesArray*)fFruits;
   return new(tracks[fNTracks++]) AliFTrack(code,charge,pT,eta,phi,
                                  v11, v22, v33, v12, v13, v23, iFlag);
}

//_____________________________________________________________________________
void AliFTrackMaker::Clear(Option_t *option)
{
   //Reset Track Maker

   fNTracks = 0;
   AliFMaker::Clear(option);
}

//_____________________________________________________________________________
void AliFTrackMaker::Draw(Option_t *)
{
//    Dummy Draw

}

//_____________________________________________________________________________
void AliFTrackMaker::Init()
{
  //Create control histograms 
  if(gAliFast->TestTrack() == 0){

     fResID11 = new TH1D("ResID11","Elec: delta(1/pTotal)*pTotal",1000,-0.5,0.5); 
     fResID12 = new TH1D("ResID12","Elec: delta(lambda)/lambda",1000,-0.01,0.01); 
     fResID13 = new TH1D("ResID13","Elec: delta(phi)/phi",1000,-0.01,0.01); 

     fResID21 = new TH1D("ResID21","Pion: delta(1/pTotal)*pTotal",1000,-1.0,1.0); 
     fResID22 = new TH1D("ResID22","Pion: delta(lambda)/lambda",1000,-1.0,1.0); 
     fResID23 = new TH1D("ResID23","Pion: delta(phi)/phi",1000,-1.0,1.0); 

     fResID31 = new TH1D("ResID31","Kaon: delta(1/pTotal)*pTotal",1000,-1.0,1.0); 
     fResID32 = new TH1D("ResID32","Kaon: delta(lambda)/lambda",1000,-1.0,1.0); 
     fResID33 = new TH1D("ResID33","Kaon: delta(phi)/phi",1000,-1.0,1.0); 

     fResID41 = new TH1D("ResID41","Proton: delta(1/pTotal)*pTotal",1000,-1.0,1.0); 
     fResID42 = new TH1D("ResID42","Proton: delta(lambda)/lambda",1000,-1.0,1.0); 
     fResID43 = new TH1D("ResID43","Proton: delta(phi)/phi",1000,-1.0,1.0); 

  }   
  //Create test histograms for TestJob only  
  if(gAliFast->TestTrack() == 1){
     fResID1Test  = new TH1D("ResID1Test","histogram21 from res.f",1000,0.075,10.075); 
     fResID2Test  = new TH1D("ResID2Test","histogram21 from res.f",1000,0.075,10.075); 
     fResID3Test  = new TH1D("ResID3Test","histogram21 from res.f",1000,0.075,10.075); 
     fResID4Test  = new TH1D("ResID4Test","histogram21 from res.f",1000,0.075,10.075); 
     fResID5Test  = new TH1D("ResID5Test","histogram21 from res.f",1000,0.075,10.075); 
  }

  //Set particle masses
   SetPionMass();
   SetKaonMass();
   SetElectronMass();
   SetProtonMass();

  //Switch on/off tracks reconstruction
   SetRecTrack();

}

//_____________________________________________________________________________
// Calculate track and its resolution
//_____________________________________________________________________________
void AliFTrackMaker::Make()
{
  Double_t v11, v22, v33, v12, v13, v23;
  Int_t iFlag;

  fNTracks = 0; 

  // Check if it is a TestJob
  if(gAliFast->TestTrack() == 1){
     // Run test job
     MakeTest(10);
  }else{
     // Run production job  
     // Get pointers to Particles arrays and TClonesArray

     TClonesArray *particles = gAliFast->Particles();	
     Int_t idPart, idTrack;
     Double_t  charge, pT, eta, phi;
     TParticle *part;
     Int_t  nparticles = particles->GetEntriesFast();
     printf("%10s%10d\n","nparticles",nparticles);
     for(Int_t ind=0;ind<nparticles;ind++) {       
       part = (TParticle*)particles->UncheckedAt(ind);
       idPart  = part->GetPdgCode();
       charge  = part->GetPDG()->Charge();
       pT      = part->Pt();  
       eta     = part->Eta();
       phi     = part->Phi();
       printf("%10s%10d%20.5e%20.5e%20.5e%20.5e\n","Particle",idPart,charge,pT,eta,phi);
       // Check convention for tracks reconstruction
       idTrack = 0;
       if(TMath::Abs(idPart) ==   11)                 idTrack = 1;
       if(TMath::Abs(idPart) == 111 || TMath::Abs(idPart) == 211) idTrack = 2;
       if(TMath::Abs(idPart) == 311 || TMath::Abs(idPart) == 321) idTrack = 3;
       if(TMath::Abs(idPart) == 2212)                 idTrack = 4;
       
       if(idTrack > 0 && fRecTrack > 0){
	 // Check if track should be reconstructed
	 if((fRecTrack == 1 && idTrack == 1) ||
	    (fRecTrack == 2 && idTrack == 2) ||
	    (fRecTrack == 3 && idTrack == 3) ||
	    (fRecTrack == 4 && idTrack == 4) ||
	    fRecTrack == 100 ) {
	   // Tracks are  reconstructed
	   ErrorMatrix(idTrack,pT,eta, v11, v22, v33, v12, v13, v23, iFlag);
	   
	   // Calculate and smear track parameters
	   Double_t lambda, cosLambda, pTotal,pInverse;
	   Double_t pInverseSmea, lambdaSmea, phiSmea;
	   Double_t a1, a2, a3, b2, b3, c3;
	   Double_t rn1, rn2, rn3;
	   
	   lambda    = kPiHalf -2.0*TMath::ATan(TMath::Exp(-eta));
	   cosLambda = TMath::Cos(lambda);
	   pTotal    = pT/cosLambda;
	   pInverse  = 1.0/pTotal;
	   
	   a1  = TMath::Sqrt(v11);
	   if(a1 == 0.){
	     a2 = 0;
	     a3 = 0;
	   }else{
	     a2  = v12/a1;
	     a3  = v13/a1;
	   }
	   b2  = TMath::Sqrt(v22-a2*a2);
	   if(b2 == 0.){
	     b3 = 0;
	   }else{
	     b3  = (v23 - a2*a3)/b2;
	   }
	   c3  = TMath::Sqrt(v33 - a3*a3 -b3*b3);
	   rn1 = gRandom->Gaus(0,1);
	   rn2 = gRandom->Gaus(0,1);
	   rn3 = gRandom->Gaus(0,1);
	   
	   pInverseSmea  = pInverse + a1*rn1;
	   lambdaSmea    = lambda + a2*rn1 + b2*rn2;
	   phiSmea       = phi + a3*rn1 + b3*rn2 + c3*rn3; 
	   
	   // Fill control histograms
	   if(idTrack == 1){
	     fResID11->Fill((pInverseSmea-pInverse)/pInverse);
	     fResID12->Fill((lambdaSmea-lambda)/lambda);
	     fResID13->Fill((phiSmea-phi)/phi);
	   }
	   else if(idTrack == 2){
	     fResID21->Fill((pInverseSmea-pInverse)/pInverse);
	     fResID22->Fill((lambdaSmea-lambda)/lambda);
	     fResID23->Fill((phiSmea-phi)/phi);
	   }
	   else if(idTrack == 3){
	     fResID31->Fill((pInverseSmea-pInverse)/pInverse);
	     fResID32->Fill((lambdaSmea-lambda)/lambda);
	     fResID33->Fill((phiSmea-phi)/phi);
	   }
	   else if(idTrack == 4){
	     fResID41->Fill((pInverseSmea-pInverse)/pInverse);
	     fResID42->Fill((lambdaSmea-lambda)/lambda);
	     fResID43->Fill((phiSmea-phi)/phi);
	   }
	 }else{
	   // Tracks are not reconstructed
	   v11=0.;
	   v12=0.;
	   v13=0.;
	   v22=0.;
	   v23=0.;
	   v33=0.;
	   iFlag=0;
	 }
	 // Store resolution variables  to AliFTrack  ClonesArray
	 AddTrack(idTrack, charge, pT, eta, phi, v11, v22, v33, v12, v13, v23, iFlag);
         printf("%10s%10d%20.5e%20.5e%20.5e%20.5e%20.5e%20.5e%20.5e%20.5e%20.5e%20.5e%10d\n",
              "Track",idTrack,charge,pT,eta,phi,v11, v22, v33, v12, v13, v23, iFlag);
       }
       
     }
  }
}

//_____________________________________________________________________________
void AliFTrackMaker::Finish()
{
  // For TestJob only  
  if(gAliFast->TestTrack() == 1){
    /*
    // Draw test histograms
    TCanvas *c1 = new TCanvas("c1"," ",200,10,600,480);
    c1->Divide(2,3);
    c1->cd(1);   fResID1Test->Draw();
    c1->cd(2);   fResID2Test->Draw();
    c1->cd(3);   fResID3Test->Draw();
    c1->cd(4);   fResID4Test->Draw();
    c1->cd(5);   fResID5Test->Draw();
    c1->Update();
    // Store TestRes.eps file
    c1->Print("TestRes.eps");
    */
    // Store histograms on file
    TFile f2("TestRes.root","RECREATE","Test Res.f");
    fResID1Test->Write();
    fResID2Test->Write();
    fResID3Test->Write();
    fResID4Test->Write();
    fResID5Test->Write();
    f2.Close();
  } 
}
//_____________________________________________________________________________
void AliFTrackMaker::ErrorMatrix(Int_t idTrack, Double_t pT,  Double_t eta,
  Double_t &v11, Double_t &v22, Double_t &v33, Double_t &v12, Double_t &v13, Double_t &v23,
  Int_t &iFlag)
{
  ///////////////////////////////////////////////
  //idTrack      track type            input   //
  //pT           transverse mom        input   //
  //lambda       deep angle            input   //
  //v11,v22,v23  error matrix          output  //
  //v12,v13,v23                        output  //
  //iFlag                              output  //
  ///////////////////////////////////////////////
 
  AliFDet *detector = gAliFast->Detector();
  Int_t nDet = detector->NDet();
  Int_t nDetActive = detector->NDetActive();
  Int_t nTwice = nDetActive + nDetActive;

  Double_t rTrack, rTrackInverse, pTotal, pInverse, diffPInverse;
  Double_t safety;
  Double_t cosLambda, tanLambda, diffLambda;
  Double_t rDet;
 
  Double_t hh0[kNMaxDet2][kNMaxDet2], hhi0[kNMaxDet2][kNMaxDet2];  
  Double_t hh1[kNMaxDet2][kNMaxDet2], hhi1[kNMaxDet2][kNMaxDet2];  
  Double_t dhhiOverPInverse[kNMaxDet2][kNMaxDet2];  
  Double_t dhhiOverLambda[kNMaxDet2][kNMaxDet2];  
  Double_t a1[kNMaxDet2][kNMaxDet2], a2[kNMaxDet2][kNMaxDet2];  
  Double_t a0PInverse[kNMaxDet2];  
  Double_t a0Lambda[kNMaxDet2];  
  Double_t a0Phi[kNMaxDet2]; 

  Double_t vF11, vF12, vF13, vF22, vF23, vF33, d1, d2, d3, det; 
  Int_t   idet, icyl, im, in;
  Double_t phiHalf;
  Double_t lambda;

  lambda = kPiHalf -2.0*TMath::ATan(TMath::Exp(-eta));
  rTrack    = detector->ConstMag()*pT;
  safety    = 10.0;
  if(2.0*rTrack < (detector->RDet(nDet) + safety)){
      iFlag     = 0;
      v11 = 0;
      v22 = 0;
      v33 = 0;
      v12 = 0;
      v13 = 0;
      v23 = 0;
      return;
  }
  iFlag        = 1;
  cosLambda    = TMath::Cos(lambda);
  pTotal       = pT/cosLambda;
  pInverse     = 1.0/pTotal;
  diffPInverse = pInverse*1.0e-5;
  diffLambda   = 1.0e-4; 

  // Compute likelihood and derivatives

  LogLikelyhood(idTrack, pInverse, lambda);
  for(icyl=1; icyl<nTwice+1; icyl++){
     for(im=1; im<nTwice+1; im++){
      hh0[icyl][im]  = HH(icyl,im); 
      hhi0[icyl][im] = HHI(icyl,im);
     }
  }
  LogLikelyhood(idTrack, pInverse+diffPInverse,lambda);
  for(icyl=1; icyl<nTwice+1; icyl++){   
     for(im=1; im<nTwice+1; im++){
      hh1[icyl][im]  = HH(icyl,im);
      hhi1[icyl][im] = HHI(icyl,im);
     }
  }  
  for(icyl=1; icyl<nTwice+1; icyl++){
     for(im=1; im<icyl+1; im++){
        dhhiOverPInverse[icyl][im] = (hhi1[icyl][im]-hhi0[icyl][im])/diffPInverse; 
     }
  }
  LogLikelyhood(idTrack, pInverse, lambda+diffLambda);
  for(icyl=1; icyl<nTwice+1; icyl++){
     for(im=1; im<nTwice+1; im++){
      hh1[icyl][im]  = HH(icyl,im);
      hhi1[icyl][im] = HHI(icyl,im);
     }
  }
  for(icyl=1; icyl<nTwice+1; icyl++){
     for(im=1; im<icyl+1; im++){
        dhhiOverLambda[icyl][im] = (hhi1[icyl][im]-hhi0[icyl][im])/diffLambda; 
     }
  }

  // Compute additional derivatives
  rTrackInverse = 1.0/rTrack;
  tanLambda    = TMath::Tan(lambda);
  icyl = 0;
  for(idet=1; idet<nDet+1;idet++){
     if(detector->IFlagDet(idet) > 0){
        icyl = icyl + 1;
        rDet = detector->RDet(idet);
        phiHalf = TMath::ASin(0.5*rDet*rTrackInverse);
        Double_t rHelp   = rDet /
                          (2.0 * TMath::Sqrt(1.0-(0.5 *rDet*rTrackInverse)*
                                                 (0.5 *rDet*rTrackInverse)));
        a0PInverse[icyl] = - rDet* rHelp
                           /(detector->ConstMag()*cosLambda);
        a0Lambda[icyl]   = - rDet* rHelp
                           * tanLambda * rTrackInverse;
        a0Phi[icyl]      =   rDet;
        a0PInverse[nDetActive+icyl] = 2.0 * tanLambda
                           *rTrack*(rHelp-rTrack*phiHalf)
                           /(detector->ConstMag()*cosLambda);
        a0Lambda[nDetActive+icyl]   = 2.0 * (  rHelp*tanLambda*tanLambda
                                             + rTrack*phiHalf);
        a0Phi[nDetActive+icyl] = 0.0 ;
    }
  }
 
  // Compute information matrix

    vF11=0.0;
    vF12=0.0;
    vF13=0.0;
    vF22=0.0;
    vF23=0.0;
    vF33=0.0;
    for(icyl=1; icyl<nTwice+1; icyl++){
       d1=0.0;     
       d2=0.0;     
       d3=0.0; 
       for(im=1; im < icyl+1; im++){
          d1 = d1 + hhi0[icyl][im]*a0PInverse[im];
          d2 = d2 + hhi0[icyl][im]*a0Lambda[im];
          d3 = d3 + hhi0[icyl][im]*a0Phi[im];
       }
       vF11 =vF11 + d1*d1;
       vF12 =vF12 + d1*d2;
       vF13 =vF13 + d1*d3;
       vF22 =vF22 + d2*d2;
       vF23 =vF23 + d2*d3;
       vF33 =vF33 + d3*d3;
    }
    for(icyl=1; icyl<nTwice+1; icyl++){
       for(im=1; im<icyl+1; im++){
          a1[icyl][im] = 0;
          a2[icyl][im] = 0;
          for(in=im; in<icyl+1;in++){
             a1[icyl][im]=a1[icyl][im]+dhhiOverPInverse[icyl][in]*hh0[im][in];
             a2[icyl][im]=a2[icyl][im]+dhhiOverLambda[icyl][in]*hh0[im][in];
	  }
          vF11=vF11+a1[icyl][im]*a1[icyl][im];
          vF12=vF12+a1[icyl][im]*a2[icyl][im];
          vF22=vF22+a2[icyl][im]*a2[icyl][im];
       }
       vF11=vF11+a1[icyl][icyl]*a1[icyl][icyl];
       vF12=vF12+a1[icyl][icyl]*a2[icyl][icyl];
       vF22=vF22+a2[icyl][icyl]*a2[icyl][icyl];
       }
 
  // Invert information matrix

    det=( vF11*vF22 - vF12*vF12 ) *vF33 + (vF12*vF23 - vF13*vF22)*vF13
                                        + (vF12*vF13 - vF11*vF23)*vF23;

    v11 = (vF22*vF33 - vF23*vF23)/det;
    v22 = (vF11*vF33 - vF13*vF13)/det;
    v33 = (vF11*vF22 - vF12*vF12)/det;
    v12 = (vF13*vF23 - vF12*vF33)/det;
    v13 = (vF12*vF23 - vF13*vF22)/det;
    v23 = (vF12*vF13 - vF11*vF23)/det;
  
    }
//_____________________________________________________________________________//
void AliFTrackMaker::LogLikelyhood(Int_t idTrack, Double_t pInverse,Double_t lambda)
{
  ///////////////////////////////////////////////
  //hh           ??                    output  //
  //hhi          ??                    output  //
  //idTrack       track type           input   //
  //pInverse      inverse  momentum    input   //
  //lambda       polar angle of track  input   //
  ///////////////////////////////////////////////

 
  AliFDet *detector = gAliFast->Detector();
  Int_t nDet = detector->NDet();
  Int_t nDetActive = detector->NDetActive();
  Int_t nTwice = nDetActive + nDetActive;

  Double_t    rDet, rDetSQ;
  Int_t      idet, icyl, im, imc;
  Double_t    cosLambda, tanLambda, pTotal, pT, rTrack, rTrackSQ;
  Double_t    beta, overPBeta, rTrackInv, thickCorr, temp1, temp2;
  Double_t    partMassSQ;
  Double_t    aShelp[kNMaxDet2], dShelp[kNMaxDet2];
  Double_t    projXVXT[kNMaxDet2],projYVXT[kNMaxDet2], projZVXT[kNMaxDet2]; 
  Double_t    proj[kNMaxDet2][kNMaxDet2]; 
  Double_t    erroScatt[kNMaxDet2], variance[kNMaxDet2][kNMaxDet2];
  Double_t    erroSQ[kNMaxDet2];
  Double_t    hh[kNMaxDet2][kNMaxDet2];
  Double_t    hhi[kNMaxDet2][kNMaxDet2];
  Double_t    errorVX, errorVY, errorVZ;

  cosLambda = TMath::Cos(lambda);
  tanLambda = TMath::Tan(lambda);
  pTotal    = 1.0/pInverse;
  pT        = pTotal * cosLambda;
  rTrack    = detector->ConstMag() * pTotal * cosLambda;
  rTrackSQ  = rTrack * rTrack;
  partMassSQ= ParticleMass(idTrack)*ParticleMass(idTrack);
  beta      = pTotal / TMath::Sqrt(partMassSQ+pTotal*pTotal);
  overPBeta = 1./(pTotal*beta);
  rTrackInv = 1./rTrack;
  errorVX   = detector->ErrorVertexX();
  errorVY   = detector->ErrorVertexY();
  errorVZ   = detector->ErrorVertexZ();

    
  erroScatt[0]=0.0;
  erroScatt[1]=0.0;
  for(idet=1; idet < nDet; idet++){
     thickCorr = detector->ThickDet(idet)/TMath::Sqrt(cosLambda*
                 TMath::Sqrt(1.0-0.25*(detector->RDetSQ(idet)/rTrackSQ)));
     if(detector->IFlagGas(idet) == 0){
         thickCorr = thickCorr * (1.3266 + 0.076 * TMath::Log(thickCorr));}
     thickCorr = overPBeta * thickCorr;
     erroScatt[idet+1]=thickCorr*thickCorr;
  }


  icyl = 0;
  for(idet=1; idet<nDet+1; idet++){
    rDet   = detector->RDet(idet);
    rDetSQ = rDet*rDet;
    dShelp[idet] = TMath::Sqrt(4.0*rTrackSQ-rDetSQ);
    aShelp[idet] = TMath::ASin(rDet/(rTrack+rTrack));
    if(detector->IFlagDet(idet) > 0) {
       icyl = icyl + 1;
       projXVXT[icyl] = rDet * rTrackInv;
       projXVXT[nDetActive+icyl] = -tanLambda;
       temp1 = (rTrackSQ + rTrackSQ - rDetSQ)/dShelp[idet];
       temp2 = rDet/dShelp[idet];
       projYVXT[icyl] = temp1*rTrackInv;
       projYVXT[nDetActive+icyl] = tanLambda * temp2;
       projZVXT[icyl] = 0.0;
       projZVXT[nDetActive+icyl] = 1.0;
       proj[icyl][1] = 0.0;
       proj[nDetActive+icyl][0] = 0.0;
       proj[nDetActive+icyl][nDet] = 0.0;
       proj[icyl][nDet] = 0.0;
       for(im=2; im<idet+1; im++){
          proj[icyl][im]= (( rDet
                           *(rTrackSQ+rTrackSQ-detector->RDetSQ(im-1))
                           - detector->RDet(im-1)*temp1*dShelp[im-1])
                         /((rTrackSQ + rTrackSQ)*cosLambda));
          proj[nDetActive+icyl][im]= 0.5 * detector->RDet(im-1)
                                   * rTrackInv
                                   * tanLambda * (detector->RDet(im-1)
                                   - dShelp[im-1]*temp2)/cosLambda;
          proj[nDetActive+icyl][nDet+im]= (rTrack+rTrack) * (aShelp[idet] - aShelp[im-1])
	     + ( rDet*dShelp[im-1]-detector->RDet(im-1)*dShelp[idet])
	     *  dShelp[im-1] * tanLambda * tanLambda
	       / (dShelp[idet] * (rTrack+rTrack));
          proj[icyl][nDet+im]= tanLambda 
                        *  (rDet*detector->RDet(im-1)*dShelp[im-1]
                        / (rTrackSQ+rTrackSQ)
                        - (rDetSQ + detector->RDetSQ(im-1)
                        -  rDetSQ * detector->RDetSQ(im-1)
                        / (rTrackSQ+rTrackSQ))
                        / dShelp[idet]);
       }
       for(im=idet+1; im < nDet+1; im++){
          proj[icyl][im] = 0.0;
          proj[nDetActive+icyl][im] = 0.0;
          proj[nDetActive+icyl][nDet+im] = 0.0;
          proj[icyl][nDet+im] = 0.0;
       }
       if(detector->IFlagDet(idet) == 1){
          erroSQ[icyl] = detector->ErrorRPhi(idet);
          erroSQ[nDetActive+icyl] = detector->ErrorZ(idet);
       }else{
       	  TPCResolution(pT, rDet, lambda);
          erroSQ[icyl]            = SigmaRPhiSQ();
          erroSQ[nDetActive+icyl] = SigmaZSQ();
       }
       erroSQ[icyl] = erroSQ[icyl] + detector->ErrorR(idet)*temp2*temp2;
       erroSQ[nDetActive+icyl] =  erroSQ[nDetActive+icyl] 
                                + detector->ErrorR(idet)*tanLambda*tanLambda;
    }
  }
    for(icyl=1; icyl<nTwice+1; icyl++){
      for(im=1; im<icyl+1; im++){
          variance[icyl][im]=
               projXVXT[icyl]*projXVXT[im]*errorVX
              +projYVXT[icyl]*projYVXT[im]*errorVY
              +projZVXT[icyl]*projZVXT[im]*errorVZ;
	  for(imc=1; imc<nDet+1; imc++){
             variance[icyl][im]=variance[icyl][im]
                    +(proj[icyl][imc]*proj[im][imc]
                    + proj[icyl][nDet+imc]*proj[im][nDet+imc])
                    * erroScatt[imc];
	  }
      }
      variance[icyl][icyl] = variance[icyl][icyl]+erroSQ[icyl];
    }
   
    for(icyl=1; icyl<nTwice+1; icyl++){
      for(im=icyl; im<nTwice+1; im++){
          hh[icyl][im]=variance[im][icyl];
          for(imc=1; imc<icyl;imc++){
               hh[icyl][im]=hh[icyl][im]-hh[imc][icyl]*hh[imc][im];
	    }
            if(im == icyl){
               hh[icyl][im] = TMath::Sqrt(hh[icyl][im]);
            }  else {
	       hh[icyl][im] = hh[icyl][im]/hh[icyl][icyl];
            }
          }
    }
       
    for(icyl=1; icyl<nTwice+1; icyl++){
        hhi[icyl][icyl] = 1.0 / hh[icyl][icyl];
        for(im=1; im<icyl; im++){
           hhi[icyl][im] = 0.0;
           for(imc=im; imc<icyl; imc++){
	      hhi[icyl][im] = hhi[icyl][im]-hh[imc][icyl]*hhi[imc][im];
           }
	   hhi[icyl][im] = hhi[icyl][im]*hhi[icyl][icyl];
       }
    }
    
    for(icyl=1; icyl<nTwice+1; icyl++){
      for(im=1; im<nTwice+1; im++){
	  SetHH(icyl,im,hh[icyl][im]);
          SetHHI(icyl,im,hhi[icyl][im]);
      }
    }

  
}
//_____________________________________________________________________________
// translation of routine tpc_resolution of res.f
//_____________________________________________________________________________
void AliFTrackMaker::TPCResolution(Double_t pTransv, Double_t radiPad, Double_t lambda)
{
  ///////////////////////////////////////////////
  //sigmaRPhiSQ  resolution in r-phi   output  //
  //sigmaZSQ     resolution in z       output  //
  //pTransv      transverse momentum   input   //
  //radiPad      radius of pad row     input   //
  //lambda       polar angle of track  input   //
  //
  //units: cm, GeV/c, radian                   //
  //parametrisation of TPC resolution          //
  //version of 03.07.1995                      //
  //source: Marek Kowalski, Karel Safarik      //
  ///////////////////////////////////////////////

  Double_t aRCoeff=0.41818e-2;
  Double_t bRCoeff=0.17460e-4;
  Double_t cRCoeff=0.30993e-8;
  Double_t dRCoeff=0.41061e-6;
  Double_t aZCoeff=0.39610e-2;
  Double_t bZCoeff=0.22443e-4;
  Double_t cZCoeff=0.51504e-1;

  Double_t sigmaRPhiSQ;
  Double_t sigmaZSQ;

  sigmaRPhiSQ = aRCoeff - bRCoeff * radiPad * TMath::Tan(lambda)+
               (cRCoeff * (radiPad/pTransv) + dRCoeff) * radiPad/pTransv;

  sigmaZSQ    = aZCoeff - bZCoeff * radiPad * TMath::Tan(lambda)+
                cZCoeff * TMath::Tan(lambda)*TMath::Tan(lambda);

  if(sigmaRPhiSQ < 1.0e-6 ) sigmaRPhiSQ = 1.0e-6;
  if(sigmaZSQ    < 1.0e-6 ) sigmaZSQ    = 1.0e-6;

  sigmaRPhiSQ =   (TMath::Sqrt(sigmaRPhiSQ) + 0.005)
                * (TMath::Sqrt(sigmaRPhiSQ) + 0.005);
  sigmaZSQ    =   (TMath::Sqrt(sigmaZSQ)    + 0.005)
                * (TMath::Sqrt(sigmaZSQ)    + 0.005);

  SetSigmaRPhiSQ(sigmaRPhiSQ);
  SetSigmaZSQ(sigmaZSQ); 

  
}

//_____________________________________________________________________________
// returns the mass given particle ID 
//-----------------------------------------------------------------------------
Double_t AliFTrackMaker::ParticleMass(Int_t idTrack)
{
  Double_t mass = 0.0;

       if(idTrack == 2){ mass = fPionMass;}
  else if(idTrack == 3){ mass = fKaonMass;}
  else if(idTrack == 1) {mass = fElectronMass;}
  else if(idTrack == 4) {mass = fProtonMass;}

  return mass;

}

//_____________________________________________________________________________
// returns the rapidity given particle pT, pz 
//-----------------------------------------------------------------------------
Double_t AliFTrackMaker::Rapidity(Double_t pt, Double_t pz)
{
//   Compute rapidity

   Double_t etalog = TMath::Log((TMath::Sqrt(pt*pt + pz*pz) + TMath::Abs(pz))/pt);
   if (pz < 0 ) return -TMath::Abs(etalog);
   else         return  TMath::Abs(etalog);
}

//_____________________________________________________________________________
// returns the phi angle given particle px, py 
//-----------------------------------------------------------------------------
Double_t AliFTrackMaker::Angle(Double_t x, Double_t y)
{
//   Compute phi angle of particle
// ... this is a copy of function ULANGL
//  .. sign(a,b) = -abs(a) if b <  0
//               =  abs(a) if b >= 0

   Double_t angle = 0;
   Double_t r = TMath::Sqrt(x*x + y*y);
   if (r < 1e-20) return angle;
   if (TMath::Abs(x)/r < 0.8) {
      angle = TMath::Sign((Double_t)TMath::Abs(TMath::ACos(x/r)), y);
   } else {
      angle = TMath::ASin(y/r);
      if (x < 0 ) {
         if(angle >= 0) angle =  kPi - angle;
         else           angle = -kPi - angle;
      }
   }
   return angle;
}
//_____________________________________________________________________________
Int_t AliFTrackMaker::Charge(Int_t kf)
{
//...this is copy of function LUCHGE 
//...Purpose: to give three times the charge for a particle/parton. 

  static Int_t kchg[500] = { -1,2,-1,2,-1,2,-1,2,0,0,-3,0,-3,0,-3,0,-3,0,
        0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,3,0,0,3,0,-1,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       2,-1,2,-1,2,3,0,0,0,0,0,0,0,0,0,0,0,3,0,3,3,0,3,0,3,0,3,0,0,0,0,0,
       0,0,0,0,0,3,0,3,3,0,3,0,3,0,3,0,0,0,0,0,0,0,0,0,0,3,0,3,3,0,3,0,3,
       0,3,0,0,0,0,0,0,0,0,0,0,3,0,3,3,0,3,0,3,0,3,0,0,0,0,0,0,0,0,0,0,3,
       0,3,3,0,3,0,3,0,3,0,0,0,0,0,0,0,0,0,0,3,0,3,3,0,3,0,3,0,3,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,
       0,3,0,0,0,0,0,0,0,0,-3,0,0,0,0,0,0,0,0,3,0,-3,0,3,-3,0,0,0,3,6,0,
       3,0,0,0,0,0,-3,0,3,-3,0,-3,0,0,0,0,-3,0,3,6,-3,0,3,-3,0,-3,0,3,6,
       0,3,0,0,0,0,0,-3,0,3,-3,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

//    extern integer kfcomp_(integer *);
  Int_t ipower;
  Int_t ret = 0;
  Int_t kfa = TMath::Abs(kf);
  Int_t kc  = Compress(kfa);

//...Initial values. Simple case of direct readout. 
  if (kc == 0) {
  } else if (kfa <= 100 || kc <= 80 || kc > 100) {
     ret = kchg[kc-1];

// ...Construction from quark content for heavy meson, diquark, baryon.
  } else if (kfa/1000 % 10 == 0) {
	ipower = kfa/100 % 10;
	ret = (kchg[kfa / 100 % 10 - 1] - kchg[kfa/10 % 10 - 1])*Int_t(TMath::Power(-1, ipower));
  } else if (kfa / 10 % 10 == 0) {
	ret = kchg[kfa/1000 % 10 - 1] + kchg[kfa/100 % 10 - 1];
  } else {
	ret = kchg[kfa/1000 % 10 - 1] + kchg[kfa/100 % 10 - 1] + kchg[kfa/10 % 10 - 1];
  }

// ...Add on correct sign.
  if (kf > 0) return ret;
  else        return -ret;
}
//_____________________________________________________________________________
Int_t AliFTrackMaker::Compress(Int_t kf)
{
//...this is copy of function LUCOMP 
//...Purpose: to compress the standard KF codes for use in mass and decay 
//...arrays; also to check whether a given code actually is defined.
//     from BLOCK LUDATA
  static Int_t  kchg[500] = { 1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,
        0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,
        1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,
        1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,
        0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,
        1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,
        0,0,1,0,1,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,
        1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,
        1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static Int_t kftab[25] = { 211,111,221,311,321,130,310,213,113,223,313,
	    323,2112,2212,210,2110,2210,110,220,330,440,30443,30553,0,0 };
  static Int_t kctab[25] = { 101,111,112,102,103,221,222,121,131,132,122,
	    123,332,333,281,282,283,284,285,286,287,231,235,0,0 };

  Int_t ret = 0;
  Int_t kfla, kflb, kflc, kflr, kfls, kfa, ikf;

  kfa = TMath::Abs(kf);
//...Simple cases: direct translation or table.
  if (kfa == 0 || kfa >= 100000) {
     return ret;
  } else if (kfa <= 100) {
     ret = kfa;
     if (kf < 0 && kchg[kfa - 1] == 0)  ret = 0;
     return ret;
  } else {
     for (ikf = 1; ikf <= 23; ++ikf) {
        if (kfa == kftab[ikf-1]) {
           ret = kctab[ikf-1];
           if (kf < 0 && kchg[ret-1] == 0)  ret = 0;
           return ret;
        }
     }
  }
// ...Subdivide KF code into constituent pieces.
  kfla = kfa / 1000%10;
  kflb = kfa / 100%10;
  kflc = kfa / 10%10;
  kfls = kfa%10;
  kflr = kfa / 10000%10;
// ...Mesons.
  if (kfa - kflr*10000 < 1000) {
     if (kflb == 0 || kflb == 9 || kflc == 0 || kflc == 9) {
     } else if (kflb < kflc) {
     } else if (kf < 0 && kflb == kflc) {
     } else if (kflb == kflc) {
        if (kflr == 0 && kfls == 1) {        ret = kflb + 110;
        } else if (kflr == 0 && kfls == 3) { ret = kflb + 130;
        } else if (kflr == 1 && kfls == 3) { ret = kflb + 150;
        } else if (kflr == 1 && kfls == 1) { ret = kflb + 170;
        } else if (kflr == 2 && kfls == 3) { ret = kflb + 190;
        } else if (kflr == 0 && kfls == 5) { ret = kflb + 210;
        }
     } else if (kflb <= 5) {
        if (kflr == 0 && kfls == 1) {        ret = (kflb-1)*(kflb-2)/2 + 100 + kflc;
        } else if (kflr == 0 && kfls == 3) { ret = (kflb-1)*(kflb-2)/2 + 120 + kflc;
        } else if (kflr == 1 && kfls == 3) { ret = (kflb-1)*(kflb-2)/2 + 140 + kflc;
        } else if (kflr == 1 && kfls == 1) { ret = (kflb-1)*(kflb-2)/2 + 160 + kflc;
        } else if (kflr == 2 && kfls == 3) { ret = (kflb-1)*(kflb-2)/2 + 180 + kflc;
        } else if (kflr == 0 && kfls == 5) { ret = (kflb-1)*(kflb-2)/2 + 200 + kflc;
        }
     } else if (kfls == 1 && kflr <= 1 || kfls == 3 && kflr <= 2 || kfls == 5 && kflr == 0) {
        ret = kflb + 80;
     }
// ...Diquarks.
  } else if ((kflr == 0 || kflr == 1) && kflc == 0) {
     if (kfls != 1 && kfls != 3) {
     } else if (kfla == 9 || kflb == 0 || kflb == 9) {
     } else if (kfla < kflb) {
     } else if (kfls == 1 && kfla == kflb) {
     } else { ret = 90;
     }
// ...Spin 1/2 baryons.
  } else if (kflr == 0 && kfls == 2) {
     if (kfla == 9 || kflb == 0 || kflb == 9 || kflc == 9) {
     } else if (kfla <= kflc || kfla < kflb) {
     } else if (kfla >= 6 || kflb >= 4 || kflc >= 4) {
         ret = kfla + 80;
     } else if (kflb < kflc) {
         ret = (kfla+1)*kfla*(kfla-1)/6 + 300 + kflc*(kflc-1)/2 + kflb;
     } else {
         ret = (kfla+1)*kfla*(kfla-1)/6 + 330 + kflb*(kflb-1)/2 + kflc;
     }
// ...Spin 3/2 baryons.
  } else if (kflr == 0 && kfls == 4) {
     if (kfla == 9 || kflb == 0 || kflb == 9 || kflc == 9) {
     } else if (kfla < kflb || kflb < kflc) {
     } else if (kfla >= 6 || kflb >= 4) {
         ret = kfla + 80;
     } else {
         ret = (kfla+1)*kfla*(kfla-1) / 6 + 360 + kflb*(kflb -1) / 2 + kflc;
     }
  }
    return ret;
}

//_____________________________________________________________________________
// TEST JOB: Calculate tracks resolution
//_____________________________________________________________________________
void AliFTrackMaker::MakeTest(Int_t n)
{
  Double_t v11, v22, v33, v12, v13, v23;
  Int_t iFlag;
  Int_t idTrack;
  Double_t  pTStart, pT, eta;

  Double_t sumDPop,sumDDip,sumDPhi;
  Double_t isum,fm;
  Double_t pTotal,partMassSQ,beta,lambda;
  Double_t dPop,dLop,dDip,dPhi,rho12,rho13,rho23;
  Double_t dPPStrag,dPPTot;
  //  Double_t resol1[1001][11],resol2[10001][11],resol3[1001][11],
  //           resol4[1001][11],resol5[10001][11]
  Double_t store1[1001],store2[10001],store3[1001],
           store4[1001],store5[10001];


  idTrack  = 2;
  pTStart = 0.07;
  for(Int_t istep=1; istep<n; istep++){
      if(istep < 100 && istep >  20) istep = istep -1 +  5;
      if(istep < 500 && istep > 100) istep = istep -1 + 25;
      if(istep <1000 && istep > 500) istep = istep -1 + 100;
      pT = pTStart + 0.01*istep;
      eta = - 0.044;
      sumDPop = 0;
      sumDDip = 0;
      sumDPhi = 0;
      isum    = 0;
      for(Int_t in=1; in<11; in++){
         eta    = eta + 0.088;
         lambda = kPiHalf -2.0*TMath::ATan(TMath::Exp(-eta));
         pTotal = pT / TMath::Cos(lambda);
         if(idTrack == 1){
           dPPStrag = 0.055 /pT;}
         else{
           partMassSQ = ParticleMass(idTrack)*ParticleMass(idTrack);
           beta       = pTotal/ TMath::Sqrt(pTotal*pTotal + partMassSQ); 
           dPPStrag   = 0.04/(pT*TMath::Power(beta,2.6666666667));
         }
 	 ErrorMatrix(idTrack,pT,eta, v11, v22, v33, v12, v13, v23, iFlag);
         if(iFlag == 1){
	    dLop   = TMath::Sqrt(v11);
            dDip   = TMath::Sqrt(v22);
            dPhi   = TMath::Sqrt(v33);
            rho12  = v12/(dLop*dDip);
            rho13  = v13/(dLop*dPhi);
            rho23  = v23/(dDip*dPhi);
            dPop   = 100. *dLop * pTotal;
            dDip   = 1000. * dDip;
            dPhi   = 1000. * dPhi;
            dPPTot = TMath::Sqrt(dPop*dPop + dPPStrag*dPPStrag);
	    //            resol1[istep][in] = dPop;
	    //            resol2[istep][in] = dDip;
	    //            resol3[istep][in] = dPhi;
	    //            resol4[istep][in] = dPPTot;
	    //            resol5[istep][in] = dPPStrag;
            sumDPop = sumDPop + dPop;
            sumDDip = sumDDip + dDip;
            sumDPhi = sumDPhi + dPhi;
            isum    = isum + 1;}
         else{
            printf("%20s %10.5f %10.5f %20s\n","pT,eta",pT,eta,"cannot smear");
         }
       }
       if(isum > 0){
         dPop   = sumDPop/isum;
         dDip   = sumDDip/isum;
         dPhi   = sumDPhi/isum;
         dPPTot = TMath::Sqrt(dPop*dPop + dPPStrag*dPPStrag);}
       else{
         dPop   = 0;
         dDip   = 0;
         dPhi   = 0;
       }
       store1[istep] = dPop;
       store2[istep] = dDip;
       store3[istep] = dPhi;
       store4[istep] = dPPTot;
       store5[istep] = dPPStrag;
       if(istep > 20 ){
          Int_t im = 5;
          if(istep > 100) {im =  25;}
          if(istep > 500) {im = 100;}
          fm = 1./(1.*im);
          for(Int_t ist=1; ist<im; ist++){
	      //	       for(Int_t in=1; in < 11; in++){
	      //                   resol1[istep-im+ist][in] = resol1[istep-im][in]+
	      //                         ist*fm*(resol1[istep][in]-resol1[istep-im][in]);
	      //                   resol2[istep-im+ist][in] = resol2[istep-im][in]+
	      //                         ist*fm*(resol2[istep][in]-resol2[istep-im][in]);
	      //                   resol3[istep-im+ist][in] = resol3[istep-im][in]+
	      //                         ist*fm*(resol3[istep][in]-resol3[istep-im][in]);
	      //                   resol4[istep-im+ist][in] = resol4[istep-im][in]+
	      //                         ist*fm*(resol4[istep][in]-resol4[istep-im][in]);
	      //                   resol5[istep-im+ist][in] = resol5[istep-im][in]+
	      //                         ist*fm*(resol5[istep][in]-resol5[istep-im][in]);
	      //	     }
	    store1[istep-im+ist]=store1[istep-im]+
	                         ist*fm*(store1[istep]-store1[istep-im]);
	    store2[istep-im+ist]=store2[istep-im]+
	                         ist*fm*(store2[istep]-store2[istep-im]);
	    store3[istep-im+ist]=store3[istep-im]+
	                         ist*fm*(store3[istep]-store3[istep-im]);
            store4[istep-im+ist]=store4[istep-im]+
	                         ist*fm*(store4[istep]-store4[istep-im]);
            store5[istep-im+ist]=store5[istep-im]+
	                         ist*fm*(store5[istep]-store5[istep-im]);
            // Fill control histograms
	    fResID1Test->Fill(pTStart + 0.01*(istep-im+ist),store1[istep-im+ist]);
	    fResID2Test->Fill(pTStart + 0.01*(istep-im+ist),store2[istep-im+ist]);
	    fResID3Test->Fill(pTStart + 0.01*(istep-im+ist),store3[istep-im+ist]);
	    fResID4Test->Fill(pTStart + 0.01*(istep-im+ist),store4[istep-im+ist]);
	    fResID5Test->Fill(pTStart + 0.01*(istep-im+ist),store5[istep-im+ist]);
          }
	  printf("%10s %10d %20.15f %20.15f %20.15f %20.15f %20.15f \n", 
	               "TestTrack:",istep,store1[istep],store2[istep],store3[istep],
                                    store4[istep],store5[istep]);
       } else {     
	  printf("%10s %10d %20.15f %20.15f %20.15f %20.15f %20.15f \n", 
	               "TestTrack:",istep,store1[istep],store2[istep],store3[istep],
                                    store4[istep],store5[istep]);
	  fResID1Test->Fill(pT,store1[istep]);
	  fResID2Test->Fill(pT,store2[istep]);
	  fResID3Test->Fill(pT,store3[istep]);
	  fResID4Test->Fill(pT,store4[istep]);
	  fResID5Test->Fill(pT,store5[istep]);
       }
  }
}
//_____________________________________________________________________________
