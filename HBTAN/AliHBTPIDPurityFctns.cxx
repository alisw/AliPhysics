#include "AliHBTPIDPurityFctns.h"
//_______________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTMonPhiResolutionVsPtFctn;
// class AliHBTMonThetaResolutionVsPtFctn;
//
// file: AliHBTPIDPurityFctns.cxx AliHBTPIDPurityFctns.h
//
// Caution: On 2D plots on X axis in simulated values
// That is contrary to two-particle resolutions where it is reconstructed one
//
// added by Piotr.Skowronski@cern.ch
// 
//
//////////////////////////////////////////////////////////////////////////////////


/******************************************************************/
/******************************************************************/

ClassImp(AliHBTMonPIDPurityVsPtFctn)

AliHBTMonPIDPurityVsPtFctn::AliHBTMonPIDPurityVsPtFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval),
 fGood(0x0),
 fAll(0x0)
{
  //ctor
  Rename("pidpurityvspt","PIDPurityVsPt");
  
}
/******************************************************************/

AliHBTMonPIDPurityVsPtFctn::~AliHBTMonPIDPurityVsPtFctn()
{
 //dtor
  delete fGood;
  delete fAll;
}
/******************************************************************/
void AliHBTMonPIDPurityVsPtFctn::Write()
{
 AliHBTMonitorFunction::Write();
 fGood->Write();
 fAll->Write();
}
/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Init()
{
//Initializes fuction
  if (AliHBTParticle::GetDebug()>0) Info("Init","%s",GetName());

  if (fResult == 0x0)
   {   
      Warning("Init","Function has NULL result histogram!");
      return;
   }
  
  if (fGood == 0x0)
   {
     TString numstr = fName + " Good";  //title and name of the
                                           //result histogram
     TAxis* xax = fResult->GetXaxis();
     fGood = new TH1D(numstr,numstr,xax->GetNbins(),xax->GetXmin(),xax->GetXmax());
   }

  if (fAll == 0x0)
   {
     TString numstr = fName + " All";  //title and name of the
                                           //result histogram
     TAxis* xax = fResult->GetXaxis();
     fAll = new TH1D(numstr,numstr,xax->GetNbins(),xax->GetXmin(),xax->GetXmax());
   }

  fResult->Reset();
  fResult->SetDirectory(0x0);
  fResult->Sumw2();
  fGood->Reset();
  fGood->SetDirectory(0x0);
  fGood->Sumw2();
  fAll->Reset();
  fAll->SetDirectory(0x0);
  fAll->Sumw2();

  if (AliHBTParticle::GetDebug()>0) Info("Init","%s Done.",GetName());
}

/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Rename(const Char_t * name)
{ 
  //Rename fuctions and all histograms belonging to it
  SetName(name);
  SetTitle(name);
  
  if (fResult)
   {
     TString numstr = fName + " Result";  //title and name of the result histogram
     fResult->SetName(numstr);
     fResult->SetTitle(numstr);
   }
  if (fGood)
   {
     TString numstr = fName + " Good";
     fGood->SetName(numstr);
     fGood->SetTitle(numstr);
   }
   
  if (fAll)
   {
     TString numstr = fName + " All";
     fAll->SetName(numstr);
     fAll->SetTitle(numstr);
   }
}
/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Rename(const Char_t * name, const Char_t * title)
{
 //renames and retitle the function and histograms
 
  SetName(name);
  SetTitle(title);
  
  if (fResult)
   {
     TString numstrn = fName + " Result";  //name of the result histogram
     TString numstrt = fTitle + " Result";  //title of the result histogram
     fResult->SetName(numstrn);
     fResult->SetTitle(numstrt);
   }
  if (fGood)
   {
     TString numstrn = fName + " Good";  //name of the Good histogram
     TString numstrt = fTitle + " Good";  //title of the Good histogram
     fGood->SetName(numstrn);
     fGood->SetTitle(numstrt);
   }

  if (fAll)
   {
     TString numstrn = fName + " All";  //name of the All histogram
     TString numstrt = fTitle + " All";  //title of the All histogram
     fAll->SetName(numstrn);
     fAll->SetTitle(numstrt);
   }

}
/******************************************************************/

TH1* AliHBTMonPIDPurityVsPtFctn::GetResult()
{
  //Returns the result of the fuction
  //that is histogram with effciency and contamination
  
  fResult->Divide(fGood,fAll);
  return fResult;
}
/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Process(AliHBTParticle * track,AliHBTParticle * part)
{
 //process the particle/track
 Double_t pt = part->Pt();
 fAll->Fill(pt);
 if (track->GetPdgCode() == part->GetPdgCode()) 
  {
    fGood->Fill(pt);
  }
// else
//  {
//    Info("Process","Catched pid impurity ...");
//  }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTMonPIDContaminationVsPtFctn)

AliHBTMonPIDContaminationVsPtFctn::AliHBTMonPIDContaminationVsPtFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval),
 fWrong(0x0),
 fAll(0x0)
{
  //ctor
  Rename("pidcontaminationvspt","PIDContaminationVsPt");
}
/******************************************************************/

AliHBTMonPIDContaminationVsPtFctn::~AliHBTMonPIDContaminationVsPtFctn()
{
 //dtor
  delete fWrong;
  delete fAll;
}
/******************************************************************/

void AliHBTMonPIDContaminationVsPtFctn::Write()
{
 //Writes the function results
 AliHBTMonitorFunction::Write();
 fWrong->Write();
 fAll->Write();
}
/******************************************************************/

void AliHBTMonPIDContaminationVsPtFctn::Init()
{
//Initializes fuction
  if (AliHBTParticle::GetDebug()>0) Info("Init","%s",GetName());

  if (fResult == 0x0)
   {   
      Warning("Init","Function has NULL result histogram!");
      return;
   }
  
  if (fWrong == 0x0)
   {
     TString numstr = fName + " Wrong";  //title and name of the
                                           //result histogram
     TAxis* xax = fResult->GetXaxis();
     fWrong = new TH1D(numstr,numstr,xax->GetNbins(),xax->GetXmin(),xax->GetXmax());
   }

  if (fAll == 0x0)
   {
     TString numstr = fName + " All";  //title and name of the
                                           //result histogram
     TAxis* xax = fResult->GetXaxis();
     fAll = new TH1D(numstr,numstr,xax->GetNbins(),xax->GetXmin(),xax->GetXmax());
   }
  fResult->Reset();
  fResult->SetDirectory(0x0);
  fResult->Sumw2();
  fWrong->Reset();
  fWrong->SetDirectory(0x0);
  fWrong->Sumw2();
  fAll->Reset();
  fAll->SetDirectory(0x0);
  fAll->Sumw2();
  
  if (AliHBTParticle::GetDebug()>0) Info("Init","%s Done.",GetName());
}

/******************************************************************/

void AliHBTMonPIDContaminationVsPtFctn::Rename(const Char_t * name)
{ 
  //Rename fuctions and all histograms belonging to it
  SetName(name);
  SetTitle(name);
  
  if (fResult)
   {
     TString numstr = fName + " Result";  //title and name of the result histogram
     fResult->SetName(numstr);
     fResult->SetTitle(numstr);
   }
  if (fWrong)
   {
     TString numstr = fName + " Wrong";
     fWrong->SetName(numstr);
     fWrong->SetTitle(numstr);
   }

  if (fAll)
   {
     TString numstrn = fName + " All";  //name of the All histogram
     TString numstrt = fTitle + " All";  //title of the All histogram
     fAll->SetName(numstrn);
     fAll->SetTitle(numstrt);
   }
}
/******************************************************************/

void AliHBTMonPIDContaminationVsPtFctn::Rename(const Char_t * name, const Char_t * title)
{
 //renames and retitle the function and histograms
 
  SetName(name);
  SetTitle(title);
  
  if (fResult)
   {
     TString numstrn = fName + " Result";  //name of the result histogram
     TString numstrt = fTitle + " Result";  //title of the result histogram
     fResult->SetName(numstrn);
     fResult->SetTitle(numstrt);
   }
  if (fWrong)
   {
     TString numstrn = fName + " Wrong";  //name of the Wrong histogram
     TString numstrt = fTitle + " Wrong";  //title of the Wrong histogram
     fWrong->SetName(numstrn);
     fWrong->SetTitle(numstrt);
   }
   
  if (fAll)
   {
     TString numstr = fName + " All";
     fAll->SetName(numstr);
     fAll->SetTitle(numstr);
   }

}
/******************************************************************/

TH1* AliHBTMonPIDContaminationVsPtFctn::GetResult()
{
  //Returns the result of the fuction
  //that is histogram with effciency and contamination
  
  fResult->Divide(fWrong,fAll);
  return fResult;
}
/******************************************************************/

void AliHBTMonPIDContaminationVsPtFctn::Process(AliHBTParticle * track, AliHBTParticle * part)
{
 //process the particle/track
 Double_t pt = part->Pt();
 fAll->Fill(pt);
 
 if (track->GetPdgCode() != part->GetPdgCode()) 
  {
//    Info("Process","Catched contamination");
//    track->Print();part->Print();
    fWrong->Fill(pt);
  }
}

/******************************************************************/
ClassImp(AliHBTPairPIDProbVsQInvFctn)

AliHBTPairPIDProbVsQInvFctn::AliHBTPairPIDProbVsQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qinvpidpur","Q_{inv}  Function");
}
/******************************************************************/

TH1* AliHBTPairPIDProbVsQInvFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/******************************************************************/

void AliHBTPairPIDProbVsQInvFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Fills the numerator using pair from the same event
   pair = CheckPair(pair);
   if(pair) fNumerator->Fill(pair->GetQInv(),pair->GetPIDProb());
}
/******************************************************************/

void AliHBTPairPIDProbVsQInvFctn::ProcessDiffEventParticles(AliHBTPair* pair)
 {
  //Fills the denumerator using mixed pairs
   pair = CheckPair(pair);
   if(pair) fDenominator->Fill(pair->GetQInv(),pair->GetPIDProb());
  }


ClassImp(AliHBTPairPIDProbVsQOutSQideQLongFctn)

AliHBTPairPIDProbVsQOutSQideQLongFctn::AliHBTPairPIDProbVsQOutSQideQLongFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("qoslpidpur","Pair PID Probablilty .vs. Q_{out}-Q_{side}-Q_{long} Fctn");
}
/*************************************************************/

void AliHBTPairPIDProbVsQOutSQideQLongFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
//Fills numerator
  pair  = CheckPair(pair);
  if (pair == 0x0) return;
  Double_t weight = pair->GetPIDProb();
  Double_t out = TMath::Abs(pair->GetQOutCMSLC());
  Double_t side = TMath::Abs(pair->GetQSideCMSLC());
  Double_t lon = TMath::Abs(pair->GetQLongCMSLC());
  fNumerator->Fill(out,side,lon,weight);
}
/*************************************************************/

void AliHBTPairPIDProbVsQOutSQideQLongFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Fills numerator
  pair  = CheckPair(pair);
  if (pair == 0x0) return;
  Double_t weight = pair->GetPIDProb();
  Double_t out = TMath::Abs(pair->GetQOutCMSLC());
  Double_t side = TMath::Abs(pair->GetQSideCMSLC());
  Double_t lon = TMath::Abs(pair->GetQLongCMSLC());
  fDenominator->Fill(out,side,lon,weight);
}
/*************************************************************/

TH1* AliHBTPairPIDProbVsQOutSQideQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
