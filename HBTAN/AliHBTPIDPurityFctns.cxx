#include "AliHBTPIDPurityFctns.h"
//_______________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTMonPhiResolutionVsPtFctn;
// class AliHBTMonThetaResolutionVsPtFctn;
//
// file: AliHBTPIDPurityFctns.cxx AliHBTPIDPurityFctns.h
//
// Classes for calculating PID purity, efficiency and other things connected with PID
// xxx
// xxx
//
// added by Piotr.Skowronski@cern.ch
//
//////////////////////////////////////////////////////////////////////////////////


/******************************************************************/
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

AliHBTMonPIDPurityVsPtFctn::AliHBTMonPIDPurityVsPtFctn(const AliHBTMonPIDPurityVsPtFctn& /*in*/):
 AliHBTMonTwoParticleFctn1D(),
 AliHBTCorrelFunction(),
 fGood(0x0),
 fAll(0x0)
{
  //cpy constructor
  Error("AliHBTMonPIDPurityVsPtFctn(const AliHBTMonPIDPurityVsPtFctn&",
        "Functions can not be copied because of histogram names clashes");
}

/******************************************************************/

AliHBTMonPIDPurityVsPtFctn::~AliHBTMonPIDPurityVsPtFctn()
{
 //dtor
  delete fGood;
  delete fAll;
}
/******************************************************************/

AliHBTMonPIDPurityVsPtFctn& AliHBTMonPIDPurityVsPtFctn::operator=
                                   (const AliHBTMonPIDPurityVsPtFctn& /*in*/)
{
//assigment operator
  Error("operator=","Functions can not be copied because of histogram names clashes");
  return *this;
}
/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Write()
{
//Writes a fucntion results
 AliHBTMonitorFunction::Write();
 fGood->Write();
 fAll->Write();
}
/******************************************************************/

void AliHBTMonPIDPurityVsPtFctn::Init()
{
//Initializes fuction
  if (AliVAODParticle::GetDebug()>0) Info("Init","%s",GetName());

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

  if (AliVAODParticle::GetDebug()>0) Info("Init","%s Done.",GetName());
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

void AliHBTMonPIDPurityVsPtFctn::Process(AliVAODParticle * track,AliVAODParticle * part)
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

AliHBTMonPIDContaminationVsPtFctn::AliHBTMonPIDContaminationVsPtFctn
                                          (const AliHBTMonPIDContaminationVsPtFctn& /*in*/):
 AliHBTMonTwoParticleFctn1D(),
 AliHBTCorrelFunction(),
 fWrong(0x0),
 fAll(0x0)
{
  //cpy constructor
  Error("AliHBTMonPIDContaminationVsPtFctn(const AliHBTMonPIDContaminationVsPtFctn&",
        "Functions can not be copied because of histogram names clashes");
}

/******************************************************************/

AliHBTMonPIDContaminationVsPtFctn::~AliHBTMonPIDContaminationVsPtFctn()
{
 //dtor
  delete fWrong;
  delete fAll;
}
/******************************************************************/

AliHBTMonPIDContaminationVsPtFctn& AliHBTMonPIDContaminationVsPtFctn::operator=
                                           (const AliHBTMonPIDContaminationVsPtFctn& /*in*/)
{
//assigment operator
  Error("operator=","Functions can not be copied because of histogram names clashes");
  return *this;
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
  if (AliVAODParticle::GetDebug()>0) Info("Init","%s",GetName());

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
  
  if (AliVAODParticle::GetDebug()>0) Info("Init","%s Done.",GetName());
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

void AliHBTMonPIDContaminationVsPtFctn::Process(AliVAODParticle * track, AliVAODParticle * part)
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
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTQInvCorrelFctnPerfectPID)

AliHBTQInvCorrelFctnPerfectPID::AliHBTQInvCorrelFctnPerfectPID(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qinvcfprfctpid","Q_{inv} Correlation Function Perfect PID");
}
/*************************************************************/

void AliHBTQInvCorrelFctnPerfectPID::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  fNumerator->Fill(trackpair->GetQInv());
}

/*************************************************************/
void AliHBTQInvCorrelFctnPerfectPID::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  fDenominator->Fill(trackpair->GetQInv());
}
/*************************************************************/

TH1* AliHBTQInvCorrelFctnPerfectPID::GetResult()
{
 //returns the scaled ratio
 
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTWeightQInvCorrelFctnPerfectPID)

AliHBTWeightQInvCorrelFctnPerfectPID::AliHBTWeightQInvCorrelFctnPerfectPID(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvcfprfctpid","Q_{inv} Weight Correlation Function Perfect PID");
}
/*************************************************************/

void AliHBTWeightQInvCorrelFctnPerfectPID::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  fNumerator->Fill(trackpair->GetQInv(),partpair->GetWeight());
}

/*************************************************************/
void AliHBTWeightQInvCorrelFctnPerfectPID::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  
  fDenominator->Fill(trackpair->GetQInv());
}
/*************************************************************/

TH1* AliHBTWeightQInvCorrelFctnPerfectPID::GetResult()
{
 //returns the scaled ratio
 
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}


/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTWeightQOutSQideQLongFctnPerfectPID)

AliHBTWeightQOutSQideQLongFctnPerfectPID::AliHBTWeightQOutSQideQLongFctnPerfectPID
  (Int_t nXbins, Double_t maxXval, Double_t minXval,
   Int_t nYbins, Double_t maxYval, Double_t minYval,
   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("wqoslprfctpid","Q_{out}-Q_{side}-Q_{long} Weight Fctn with perfect PID");
}
/*************************************************************/

void AliHBTWeightQOutSQideQLongFctnPerfectPID::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  Double_t weight = partpair->GetWeight();
  Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
  Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
  Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
  fNumerator->Fill(out,side,lon,weight);
}
/*************************************************************/

void AliHBTWeightQOutSQideQLongFctnPerfectPID::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
  Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
  Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
  fDenominator->Fill(out,side,lon);
}
/******************************************************************/

TH1* AliHBTWeightQOutSQideQLongFctnPerfectPID::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}


/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTQOutSQideQLongFctnPerfectPID)

AliHBTQOutSQideQLongFctnPerfectPID::AliHBTQOutSQideQLongFctnPerfectPID
   (Int_t nXbins, Double_t maxXval, Double_t minXval,
    Int_t nYbins, Double_t maxYval, Double_t minYval,
    Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("qoslprfctpid","Q_{out}-Q_{side}-Q_{long} Fctn with perfect PID");
}
/*************************************************************/

void AliHBTQOutSQideQLongFctnPerfectPID::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  
  Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
  Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
  Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
  fNumerator->Fill(out,side,lon);
}
/*************************************************************/

void AliHBTQOutSQideQLongFctnPerfectPID::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  
  Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
  Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
  Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
  fDenominator->Fill(out,side,lon);
}
/******************************************************************/

TH1* AliHBTQOutSQideQLongFctnPerfectPID::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/******************************************************************/
/******************************************************************/
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

/******************************************************************/
/******************************************************************/
/******************************************************************/

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
  Double_t out = TMath::Abs(pair->GetQOutLCMS());
  Double_t side = TMath::Abs(pair->GetQSideLCMS());
  Double_t lon = TMath::Abs(pair->GetQLongLCMS());
  fNumerator->Fill(out,side,lon,weight);
}
/*************************************************************/

void AliHBTPairPIDProbVsQOutSQideQLongFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Fills numerator
  pair  = CheckPair(pair);
  if (pair == 0x0) return;
  Double_t weight = pair->GetPIDProb();
  Double_t out = TMath::Abs(pair->GetQOutLCMS());
  Double_t side = TMath::Abs(pair->GetQSideLCMS());
  Double_t lon = TMath::Abs(pair->GetQLongLCMS());
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

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID)

AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID::AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
//ctor
//Set Axis Title
 fWriteNumAndDen = kTRUE;
 Rename("tteffptthetaphipfctPID","P_{t} \\theta \\phi Two Track Efficiency Function");
 if(fNumerator)
  {
   fNumerator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fNumerator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fNumerator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }

 if(fDenominator)
  {
   fDenominator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fDenominator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fDenominator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }
}
/******************************************************************/

void AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  Double_t x = TMath::Abs(trackpair->GetDeltaPt());
  Double_t y = TMath::Abs(trackpair->GetDeltaTheta());
  Double_t z = TMath::Abs(trackpair->GetDeltaPhi());
  fNumerator->Fill(x,y,z);
}
/******************************************************************/

void AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  if (trackpair->Particle1()->GetPdgCode() != partpair->Particle1()->GetPdgCode()) return;
  if (trackpair->Particle2()->GetPdgCode() != partpair->Particle2()->GetPdgCode()) return;

  trackpair  = CheckPair(trackpair);
  if (trackpair == 0x0) return;
//  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;

  Double_t x = TMath::Abs(trackpair->GetDeltaPt());
  Double_t y = TMath::Abs(trackpair->GetDeltaTheta());
  Double_t z = TMath::Abs(trackpair->GetDeltaPhi());
  fDenominator->Fill(x,y,z);
}
/******************************************************************/

TH1* AliHBTTwoTrackEffFctnPtThetaPhiPerfectPID::GetResult()
{
//returns ratio of numerator and denominator
 delete fRatio;
 fRatio = GetRatio(Scale());
 if(fRatio)
  {
   fRatio->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fRatio->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fRatio->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }
 return fRatio;
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTPairPIDProbVsPtThetaPhiFctn)

AliHBTPairPIDProbVsPtThetaPhiFctn::AliHBTPairPIDProbVsPtThetaPhiFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("ptthetaphipidpur","Pair PID Probablilty .vs. \\Delta P_{t}-\\Delta \\theta-\\Delta \\phi Fctn");
  if(fNumerator)
   {
    fNumerator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
    fNumerator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
    fNumerator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
   }

  if(fDenominator)
   {
    fDenominator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
    fDenominator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
    fDenominator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
   }
  
}
/*************************************************************/

void AliHBTPairPIDProbVsPtThetaPhiFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
//Fills numerator
  pair  = CheckPair(pair);
  if (pair == 0x0) return;
  Double_t weight = pair->GetPIDProb();
  Double_t pt = TMath::Abs(pair->GetDeltaPt());
  Double_t theta = TMath::Abs(pair->GetDeltaTheta());
  Double_t phi = TMath::Abs(pair->GetDeltaPhi());
  fNumerator->Fill(pt,theta,phi,weight);
}
/*************************************************************/

void AliHBTPairPIDProbVsPtThetaPhiFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Fills numerator
  pair  = CheckPair(pair);
  if (pair == 0x0) return;
  Double_t weight = pair->GetPIDProb();
  Double_t pt = TMath::Abs(pair->GetDeltaPt());
  Double_t phi = TMath::Abs(pair->GetDeltaTheta());
  Double_t theta = TMath::Abs(pair->GetDeltaPhi());
  fDenominator->Fill(pt,theta,phi,weight);
}
/*************************************************************/

TH1* AliHBTPairPIDProbVsPtThetaPhiFctn::GetResult()
{
 //returns the scaled ratio
 
 delete fRatio;
 fRatio = GetRatio(Scale());
 if(fRatio)
  {
   fRatio->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fRatio->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fRatio->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }
 return fRatio;
}

