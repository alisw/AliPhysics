/* $Id $ */

//-----------------------------------------------------------
//This class introduces the weights calculated according 
//with functions of efficiency of identification (TPC+TOF) 
//(calculated by B.V. Batyunia).
//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//-----------------------------------------------------------

#include "AliHBTWeightsPID.h"
#include "AliHBTPair.h"
#include "AliHBTParticle.h"
#include <TRandom.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>

ClassImp(AliHBTWeightsPID)  
 
AliHBTWeightsPID* AliHBTWeightsPID::fgWeightsPID=NULL; 

AliHBTWeightsPID::AliHBTWeightsPID()
{
  //ctor
  //initial parameters of model

  fPtK   = new TH1F("fPtK"," pt of K+- ",40,0,4);
  fPtKefftpc   = new TH1F("fPtEfficKpos"," pt of K+ after efficiency cut ",40,0,4);
  fPtKefftpcboth   = new TH1F("fPtEfficKneg"," pt of K- after efficiency cut ",40,0,4);

  // efficiency functions (from pt) for K+- in the TPC
  fEffic1pol = new TF1("fEffic1pol4","pol4",100.,480.); 
  // for -0.9 < eta < -0.7
  // for  0.7 < eta <  0.9
  fEffic1pol->SetParameter(0,-0.266362);
  fEffic1pol->SetParameter(1,0.00565461);
  fEffic1pol->SetParameter(2,-4.06686e-05);
  fEffic1pol->SetParameter(3,1.39387e-07);
  fEffic1pol->SetParameter(4,-1.59674e-10);

  fEffic2pol = new TF1("fEffic2pol4","pol4",100.,540.); 

  fEffic2pol->SetParameter(0,-0.324881);
  fEffic2pol->SetParameter(1,0.00565381);
  fEffic2pol->SetParameter(2,-3.23633e-05);
  fEffic2pol->SetParameter(3,9.72523e-08);
  fEffic2pol->SetParameter(4,-1.01013e-10);
                                                                               
  fEffic3pol = new TF1("fEffic3pol4","pol4",100.,585.); 
  // for -0.5 < eta < -0.3
  // for  0.3 < eta <  0.5
  fEffic3pol->SetParameter(0,-0.306572);
  fEffic3pol->SetParameter(1,0.00557472);
  fEffic3pol->SetParameter(2,-3.33752e-05);
  fEffic3pol->SetParameter(3,9.83241e-08);
  fEffic3pol->SetParameter(4,-9.5827e-11);
  fEffic4pol = new TF1("fEffic4pol4","pol4",100.,600.);
  // for -0.3 < eta < 0.3
  fEffic4pol->SetParameter(0,-0.168648);
  fEffic4pol->SetParameter(1,0.00252021);
  fEffic4pol->SetParameter(2,-1.09113e-05);
  fEffic4pol->SetParameter(3,3.34871e-08);
  fEffic4pol->SetParameter(4,-3.31691e-11);

  // efficiency functions (from pt) for K+- in the TOF
  fEffic1polTOF = new TF1("fEffic1pol4TOF","pol4",0.2,2.0); 
  // for -0.9 < eta < -0.7
  // for  0.7 < eta <  0.9
  fEffic1polTOF->SetParameter(0,-0.165885);
  fEffic1polTOF->SetParameter(1,0.717459);
  fEffic1polTOF->SetParameter(2,-0.457131);
  fEffic1polTOF->SetParameter(3,0.284753);
  fEffic1polTOF->SetParameter(4,-0.105215);
                                                                               
  fEffic2polTOF = new TF1("fEffic2pol4TOF","pol4",0.2,2.4); 
  // for -0.7 < eta < -0.5
  // for  0.5 < eta <  0.7

  fEffic2polTOF->SetParameter(0,-0.165947);
  fEffic2polTOF->SetParameter(1,0.702475);
  fEffic2polTOF->SetParameter(2,-0.300313);
  fEffic2polTOF->SetParameter(3,0.127047);
  fEffic2polTOF->SetParameter(4,-0.0489395);
                                                                               

  fEffic3polTOF = new TF1("fEffic3pol4TOF","pol4",0.2,2.4); 
  // for -0.5 < eta < -0.3
  // for  0.3 < eta <  0.5
  fEffic3polTOF->SetParameter(0,-0.339516);
  fEffic3polTOF->SetParameter(1,1.56942);
  fEffic3polTOF->SetParameter(2,-1.43132);
  fEffic3polTOF->SetParameter(3,0.727148);
  fEffic3polTOF->SetParameter(4,-0.158444);
  fEffic4polTOF = new TF1("fEffic4pol4TOF","pol4",0.2,2.6);
  // for -0.3 < eta < 0.3
  fEffic4polTOF->SetParameter(0,-0.243435);
  fEffic4polTOF->SetParameter(1,1.00928);
  fEffic4polTOF->SetParameter(2,-0.594597);
  fEffic4polTOF->SetParameter(3,0.212601);
  fEffic4polTOF->SetParameter(4,-0.0453419);

}                                            

AliHBTWeightsPID* AliHBTWeightsPID::Instance()
{
  //Creates an instance of the class 
  //or returns pointer to already existing one
  if (fgWeightsPID) {
    return fgWeightsPID;
  } else {
    fgWeightsPID = new AliHBTWeightsPID();
    return fgWeightsPID;
  }
}


Double_t AliHBTWeightsPID::GetWeightPID(const AliHBTPair* trackpair)
{
  //Calculates the weight of "trackpair"
  AliVAODParticle *track1 = trackpair->Particle1();
  AliVAODParticle *track2 = trackpair->Particle2();
    
  Double_t pt1=track1->Pt();
  Double_t eta1=track1->Eta();

  fPtK->Fill(pt1); 

  if(TMath::Abs(eta1) > 0.7 && TMath::Abs(eta1) < 0.9 && pt1 < 0.48) {
    fEfficTPC1 = fEffic1pol->Eval(pt1*1000.);
  }else if(TMath::Abs(eta1) > 0.5 && TMath::Abs(eta1) < 0.7 && pt1 < 0.54) {
    fEfficTPC1 = fEffic2pol->Eval(pt1*1000.);
  }else if(TMath::Abs(eta1) > 0.3 && TMath::Abs(eta1) < 0.5 && pt1 < 0.585) {
    fEfficTPC1 = fEffic3pol->Eval(pt1*1000.);
  }else if(eta1 > -0.3 && eta1 < 0.3 && pt1 < 0.6) {
    fEfficTPC1 = fEffic4pol->Eval(pt1*1000.);
  }
  
  // TOF efficiency
  
  if(pt1 > 0.2) {
    if(TMath::Abs(eta1) > 0.7 && TMath::Abs(eta1) < 0.9 && pt1 < 2.0) {
      fEfficTOF1 = fEffic1polTOF->Eval(pt1);
    }else if(TMath::Abs(eta1) > 0.5 && TMath::Abs(eta1) < 0.7 && pt1 < 2.4) {
      fEfficTOF1 = fEffic2polTOF->Eval(pt1);
    }else if(TMath::Abs(eta1) > 0.3 && TMath::Abs(eta1) < 0.5 && pt1 < 2.4) {
      fEfficTOF1 = fEffic3polTOF->Eval(pt1);
    }else if(eta1 > -0.3 && eta1 < 0.3 && pt1 < 2.6) {
      fEfficTOF1 = fEffic4polTOF->Eval(pt1);
    }
  }
  
  Double_t rndmtpc=gRandom->Rndm();
  Double_t rndmtof=gRandom->Rndm();
  Double_t weightPID1=1.;
  if(fEfficTPC1 < rndmtpc && fEfficTOF1 < rndmtof) { weightPID1=0.;}
  fPtKefftpc->Fill(pt1,weightPID1); 
  
  
  Double_t pt2=track2->Pt();
  Double_t eta2=track2->Eta();
  
  // TPC efficiency

  if(TMath::Abs(eta2) > 0.7 && TMath::Abs(eta2) < 0.9 && pt2 < 0.48) {
    fEfficTPC1 = fEffic1pol->Eval(pt2*1000.);
  }else if(TMath::Abs(eta2) > 0.5 && TMath::Abs(eta2) < 0.7 && pt2 < 0.54) {
    fEfficTPC1 = fEffic2pol->Eval(pt2*1000.);
  }else if(TMath::Abs(eta2) > 0.3 && TMath::Abs(eta2) < 0.5 && pt2 < 0.585) {
    fEfficTPC1 = fEffic3pol->Eval(pt2*1000.);
  }else if(eta2 > -0.3 && eta2 < 0.3 && pt2 < 0.6) {
    fEfficTPC1 = fEffic4pol->Eval(pt2*1000.);
  }
  
  // TOF efficiency
  
  if(pt2 > 0.2) {
    if(TMath::Abs(eta2) > 0.7 && TMath::Abs(eta2) < 0.9 && pt2 < 2.0) {
      fEfficTOF1 = fEffic1polTOF->Eval(pt2);
    }else if(TMath::Abs(eta2) > 0.5 && TMath::Abs(eta2) < 0.7 && pt2 < 2.4) {
      fEfficTOF1 = fEffic2polTOF->Eval(pt2);
    }else if(TMath::Abs(eta2) > 0.3 && TMath::Abs(eta2) < 0.5 && pt2 < 2.4) {
      fEfficTOF1 = fEffic3polTOF->Eval(pt2);
    }else if(eta2 > -0.3 && eta2 < 0.3 && pt2 < 2.6) {
      fEfficTOF1 = fEffic4polTOF->Eval(pt2);
    }
  }
  
  rndmtpc=gRandom->Rndm();
  rndmtof=gRandom->Rndm();
  
  Double_t weightPID2=1.;
  if(fEfficTPC1 < rndmtpc && fEfficTOF1 < rndmtof) { weightPID2=0.;}
  
  Double_t weightPID=weightPID1*weightPID2;
  fPtKefftpcboth->Fill(pt1,weightPID); 
  
  return weightPID;
}
