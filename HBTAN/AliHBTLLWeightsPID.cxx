//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
#include "AliHBTLLWeightsPID.h"
#include "AliPDG.h"
#include "AliHBTPair.h"
#include "AliHBTParticle.h"
#include <TList.h>
#include <TRandom.h>                                                                     
#include <TMath.h>                                                                       


ClassImp(AliHBTLLWeightsPID)  
 
AliHBTLLWeightsPID* AliHBTLLWeightsPID::fWeightsPID=NULL; 

AliHBTLLWeightsPID::AliHBTLLWeightsPID()
{                                                                                           
//ctor
//initial parameters of model

  ptK   = new TH1F("ptK"," pt of K+- ",40,0,4);
  ptKefftpc   = new TH1F("ptEfficKpos"," pt of K+ after efficiency cut ",40,0,4);
  ptKefftpcboth   = new TH1F("ptEfficKneg"," pt of K- after efficiency cut ",40,0,4);

// efficiency functions (from pt) for K+- in the TPC
  effic1pol = new TF1("effic1pol4","pol4",100.,480.); 
                                                      // for -0.9 < eta < -0.7
                                                      // for  0.7 < eta <  0.9
  effic1pol->SetParameter(0,-0.266362);
  effic1pol->SetParameter(1,0.00565461);
  effic1pol->SetParameter(2,-4.06686e-05);
  effic1pol->SetParameter(3,1.39387e-07);
  effic1pol->SetParameter(4,-1.59674e-10);

  effic2pol = new TF1("effic2pol4","pol4",100.,540.); 

  effic2pol->SetParameter(0,-0.324881);
  effic2pol->SetParameter(1,0.00565381);
  effic2pol->SetParameter(2,-3.23633e-05);
  effic2pol->SetParameter(3,9.72523e-08);
  effic2pol->SetParameter(4,-1.01013e-10);
                                                                               
  effic3pol = new TF1("effic3pol4","pol4",100.,585.); 
                                                      // for -0.5 < eta < -0.3
                                                      // for  0.3 < eta <  0.5
  effic3pol->SetParameter(0,-0.306572);
  effic3pol->SetParameter(1,0.00557472);
  effic3pol->SetParameter(2,-3.33752e-05);
  effic3pol->SetParameter(3,9.83241e-08);
  effic3pol->SetParameter(4,-9.5827e-11);
                                                                                
  effic4pol = new TF1("effic4pol4","pol4",100.,600.); // for -0.3 < eta < 0.3
  effic4pol->SetParameter(0,-0.168648);
  effic4pol->SetParameter(1,0.00252021);
  effic4pol->SetParameter(2,-1.09113e-05);
  effic4pol->SetParameter(3,3.34871e-08);
  effic4pol->SetParameter(4,-3.31691e-11);

  // efficiency functions (from pt) for K+- in the TOF
  effic1polTOF = new TF1("effic1pol4TOF","pol4",0.2,2.0); 
                                                          // for -0.9 < eta < -0.7
                                                          // for  0.7 < eta <  0.9
  effic1polTOF->SetParameter(0,-0.165885);
  effic1polTOF->SetParameter(1,0.717459);
  effic1polTOF->SetParameter(2,-0.457131);
  effic1polTOF->SetParameter(3,0.284753);
  effic1polTOF->SetParameter(4,-0.105215);
                                                                               
  effic2polTOF = new TF1("effic2pol4TOF","pol4",0.2,2.4); 
                                                          // for -0.7 < eta < -0.5
                                                          // for  0.5 < eta <  0.7

  effic2polTOF->SetParameter(0,-0.165947);
  effic2polTOF->SetParameter(1,0.702475);
  effic2polTOF->SetParameter(2,-0.300313);
  effic2polTOF->SetParameter(3,0.127047);
  effic2polTOF->SetParameter(4,-0.0489395);
                                                                               

  effic3polTOF = new TF1("effic3pol4TOF","pol4",0.2,2.4); 
                                                          // for -0.5 < eta < -0.3
                                                          // for  0.3 < eta <  0.5
  effic3polTOF->SetParameter(0,-0.339516);
  effic3polTOF->SetParameter(1,1.56942);
  effic3polTOF->SetParameter(2,-1.43132);
  effic3polTOF->SetParameter(3,0.727148);
  effic3polTOF->SetParameter(4,-0.158444);
                                                                                
  effic4polTOF = new TF1("effic4pol4TOF","pol4",0.2,2.6); // for -0.3 < eta < 0.3
  effic4polTOF->SetParameter(0,-0.243435);
  effic4polTOF->SetParameter(1,1.00928);
  effic4polTOF->SetParameter(2,-0.594597);
  effic4polTOF->SetParameter(3,0.212601);
  effic4polTOF->SetParameter(4,-0.0453419);

}                                                                                             
                                            

AliHBTLLWeightsPID* AliHBTLLWeightsPID::Instance()
{                                                                                             
  if (fWeightsPID) {                                                                        
    return fWeightsPID;                                                                   
   } else {                                                                                  
   fWeightsPID = new AliHBTLLWeightsPID();                                                        
      return fWeightsPID;                                                                   
  }                                                                                         
}                                                                                             



Double_t AliHBTLLWeightsPID::GetWeightPID(const AliHBTPair* trackpair)
{

    AliHBTParticle *track1 = trackpair->Particle1();
    AliHBTParticle *track2 = trackpair->Particle2();
    
    Double_t Pt1=track1->Pt();
    Double_t eta1=track1->Eta();

    ptK->Fill(Pt1); 

          if(TMath::Abs(eta1) > 0.7 && TMath::Abs(eta1) < 0.9 && Pt1 < 0.48) {
            efficTPC1 = effic1pol->Eval(Pt1*1000.);
          }else if(TMath::Abs(eta1) > 0.5 && TMath::Abs(eta1) < 0.7 && Pt1 < 0.54) {
            efficTPC1 = effic2pol->Eval(Pt1*1000.);
          }else if(TMath::Abs(eta1) > 0.3 && TMath::Abs(eta1) < 0.5 && Pt1 < 0.585) {
            efficTPC1 = effic3pol->Eval(Pt1*1000.);
          }else if(eta1 > -0.3 && eta1 < 0.3 && Pt1 < 0.6) {
            efficTPC1 = effic4pol->Eval(Pt1*1000.);
          }

          // TOF efficiency

          if(Pt1 > 0.2) {
           if(TMath::Abs(eta1) > 0.7 && TMath::Abs(eta1) < 0.9 && Pt1 < 2.0) {
            efficTOF1 = effic1polTOF->Eval(Pt1);
           }else if(TMath::Abs(eta1) > 0.5 && TMath::Abs(eta1) < 0.7 && Pt1 < 2.4) {
            efficTOF1 = effic2polTOF->Eval(Pt1);
           }else if(TMath::Abs(eta1) > 0.3 && TMath::Abs(eta1) < 0.5 && Pt1 < 2.4) {
            efficTOF1 = effic3polTOF->Eval(Pt1);
           }else if(eta1 > -0.3 && eta1 < 0.3 && Pt1 < 2.6) {
            efficTOF1 = effic4polTOF->Eval(Pt1);
           }
          }

          Double_t rndmtpc=gRandom->Rndm();
          Double_t rndmtof=gRandom->Rndm();
          Double_t weightPID1=1.;
          if(efficTPC1 < rndmtpc && efficTOF1 < rndmtof) { weightPID1=0.;}
          ptKefftpc->Fill(Pt1,weightPID1); 


    Double_t Pt2=track2->Pt();
    Double_t eta2=track2->Eta();

// TPC efficiency

          if(TMath::Abs(eta2) > 0.7 && TMath::Abs(eta2) < 0.9 && Pt2 < 0.48) {
            efficTPC1 = effic1pol->Eval(Pt2*1000.);
          }else if(TMath::Abs(eta2) > 0.5 && TMath::Abs(eta2) < 0.7 && Pt2 < 0.54) {
            efficTPC1 = effic2pol->Eval(Pt2*1000.);
          }else if(TMath::Abs(eta2) > 0.3 && TMath::Abs(eta2) < 0.5 && Pt2 < 0.585) {
            efficTPC1 = effic3pol->Eval(Pt2*1000.);
          }else if(eta2 > -0.3 && eta2 < 0.3 && Pt2 < 0.6) {
            efficTPC1 = effic4pol->Eval(Pt2*1000.);
          }

          // TOF efficiency

          if(Pt2 > 0.2) {
           if(TMath::Abs(eta2) > 0.7 && TMath::Abs(eta2) < 0.9 && Pt2 < 2.0) {
            efficTOF1 = effic1polTOF->Eval(Pt2);
           }else if(TMath::Abs(eta2) > 0.5 && TMath::Abs(eta2) < 0.7 && Pt2 < 2.4) {
            efficTOF1 = effic2polTOF->Eval(Pt2);
           }else if(TMath::Abs(eta2) > 0.3 && TMath::Abs(eta2) < 0.5 && Pt2 < 2.4) {
            efficTOF1 = effic3polTOF->Eval(Pt2);
           }else if(eta2 > -0.3 && eta2 < 0.3 && Pt2 < 2.6) {
            efficTOF1 = effic4polTOF->Eval(Pt2);
           }
          }

          rndmtpc=gRandom->Rndm();
          rndmtof=gRandom->Rndm();

          Double_t weightPID2=1.;
          if(efficTPC1 < rndmtpc && efficTOF1 < rndmtof) { weightPID2=0.;}

          Double_t weightPID=weightPID1*weightPID2;
          ptKefftpcboth->Fill(Pt1,weightPID); 

          return weightPID;
}
