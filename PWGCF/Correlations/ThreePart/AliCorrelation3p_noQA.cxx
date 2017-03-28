/*************************************************************************
* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

#include "AliCorrelation3p_noQA.h"
#include "AliVParticle.h"
#include "AliCFPI0.h"
#include "AliFilteredTrack.h"
#include "AliLog.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include <cerrno>
#include <memory>
#include <set>
#include "TParameter.h"
#include "TF1.h"

const double gkPii =  TMath::Pi();

ClassImp(AliCorrelation3p_noQA)

AliCorrelation3p_noQA::AliCorrelation3p_noQA(const char* name,TArrayD MBinEdges, TArrayD ZBinEdges)
  : TNamed(name?name:"AliCorrelation3p", "")
  , fHistograms(NULL)
  , fMinTriggerPt(8.0)
  , fMaxTriggerPt(15.0)
  , fMinAssociatedPt(3.)
  , fMaxAssociatedPt(fMinTriggerPt)
  , fhPhiEtaDeltaPhi12Cut1(0.5*TMath::Pi())
  , fhPhiEtaDeltaPhi12Cut2(0.25*TMath::Pi())
  , fAcceptanceCut(0.8)
  , fMixedEvent(NULL)
  , fMBinEdges(MBinEdges)
  , fZBinEdges(ZBinEdges)
  , fMultiplicity(0)
  , fVZ(0)
  , fMBin(0)
  , fVzBin(0)
  , fbinver(0)
  , fCollisionType(PbPb)
  , fTriggerType(tracks)
{
  // default constructor
}

AliCorrelation3p_noQA::AliCorrelation3p_noQA(const AliCorrelation3p_noQA& other)
  : TNamed(other)
  , fHistograms((other.fHistograms!=NULL?(new TObjArray(*other.fHistograms)):NULL))
  , fMinTriggerPt(other.fMinTriggerPt)
  , fMaxTriggerPt(other.fMaxTriggerPt)
  , fMinAssociatedPt(other.fMinAssociatedPt)
  , fMaxAssociatedPt(other.fMaxAssociatedPt)
  , fhPhiEtaDeltaPhi12Cut1(other.fhPhiEtaDeltaPhi12Cut1)
  , fhPhiEtaDeltaPhi12Cut2(other.fhPhiEtaDeltaPhi12Cut2)
  , fAcceptanceCut(other.fAcceptanceCut)
  , fMixedEvent((other.fMixedEvent!=NULL?(new AliCorrelation3p_noQA(*other.fMixedEvent)):NULL))
  , fMBinEdges(other.fMBinEdges)
  , fZBinEdges(other.fZBinEdges)
  , fMultiplicity(0)
  , fVZ(0)
  , fMBin(0)
  , fVzBin(0)
  , fbinver(other.fbinver)
  , fCollisionType(other.fCollisionType)
  , fTriggerType(other.fTriggerType)
{
  // copy constructor
}

AliCorrelation3p_noQA& AliCorrelation3p_noQA::operator=(const AliCorrelation3p_noQA& other)
{
  // assignment operator
  if (this==&other) return *this;
  this->~AliCorrelation3p_noQA();
  new (this) AliCorrelation3p_noQA(other);
  return *this;
}

AliCorrelation3p_noQA::~AliCorrelation3p_noQA()
{
  // destructor
  //
  //
  if (fHistograms) {
    delete fHistograms;
    fHistograms=NULL;
  }
  // note: mixed event is an external pointer
  fMixedEvent=NULL;
}

void AliCorrelation3p_noQA::Copy(TObject &object) const
{
  /// overloaded from TObject: copy to target object
  AliCorrelation3p_noQA* target=dynamic_cast<AliCorrelation3p_noQA*>(&object);
  if (!target) return;

  AliCorrelation3p_noQA* backupME=fMixedEvent;
  if (!target->fMixedEvent) {
    // avoid copying the mixed event object if there is no target
    const_cast<AliCorrelation3p_noQA*>(this)->fMixedEvent=NULL;
  }
  *target=*this;
  const_cast<AliCorrelation3p_noQA*>(this)->fMixedEvent=backupME;
}

int AliCorrelation3p_noQA::Init(const char* arguments)
{
  /// init class and create histograms
  const char* key=NULL;
  const TString delimiter(" ");
  TStringToken token(arguments, delimiter);
  while (token.NextToken()) {
    key="minTriggerPt=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      fMinTriggerPt=param.Atof();
      continue;
    }
    key="maxTriggerPt=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      fMaxTriggerPt=param.Atof();
      continue;
    }
    key="minAssociatedPt=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      fMinAssociatedPt=param.Atof();
      continue;
    }
    key="maxAssociatedPt=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      fMaxAssociatedPt=param.Atof();
      continue;
    }
     key="collisiontype=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      if(param.CompareTo("pp")==0)   fCollisionType=pp;
      if(param.CompareTo("PbPb")==0) fCollisionType=PbPb;
      if(param.CompareTo("pPb")==0)  fCollisionType=pPb;
      AliWarning(Form("Collision Type set to: %s",param.Data()));
      continue;
    }
     key="triggertype=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      if(param.CompareTo("tracks")==0) fTriggerType=tracks;
      if(param.CompareTo("pi0")==0)    fTriggerType=pi0;
      AliWarning(Form("Trigger Type set to: %s",param.Data()));
      continue;
    }
  }
  if (fHistograms) delete fHistograms;
  fHistograms=new TObjArray(GetNumberHist(kNofHistograms+2,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1));//One Extra for overflow
  if (!fHistograms) return -ENOMEM;
  fHistograms->SetOwner(kTRUE);
  TString infoTriggerPt;
  if (fMaxTriggerPt > fMinTriggerPt) infoTriggerPt.Form("%f < pt < %f", fMinTriggerPt, fMaxTriggerPt);
  else infoTriggerPt.Form("pt > %f", fMinTriggerPt);
  TString infoAssociatedPt;
  if (fMaxAssociatedPt > fMinAssociatedPt) infoAssociatedPt.Form("%f < pt < %f", fMinAssociatedPt, fMaxAssociatedPt);
  else infoAssociatedPt.Form("pt > %f", fMinAssociatedPt);
  AliInfo(Form("initializing %s for trigger %s and associated particle %s", GetName(), infoTriggerPt.Data(), infoAssociatedPt.Data()));
  // avoid the objects to be added to the global directory
  bool statusAddDirectory=TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  const double Pii=TMath::Pi();
  TObjArray* a=fHistograms;
  //histograms that are not binned
  a->AddAt(new TH2D("centvsvz","centrality vs vz",100,fMBinEdges.At(0),fMBinEdges.At(fMBinEdges.GetSize()-1),100,fZBinEdges.At(0),fZBinEdges.At(fZBinEdges.GetSize()-1)),kcentrvsvz);
  a->AddAt(new TH2D("centbinvsvzbin","centrality bin vs vz bin",fMBinEdges.GetSize()-1,0.0,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1,0,fZBinEdges.GetSize()-1),kcentrvsvzbin);
  //Create one set of histograms per mixing bin
  //set binning:
  int nbinseta;
  int nbinsphi;
  if(fbinver==0){nbinseta = 63;nbinsphi = 38;}
  if(fbinver==1){nbinseta = 31;nbinsphi = 18;}
  if(!(fbinver==0||fbinver==1)){nbinseta = 63;nbinsphi = 38;}//revert to original
  for(Int_t i=0;i<fMBinEdges.GetSize()-1;i++){
    for(Int_t j=0;j<fZBinEdges.GetSize()-1;j++){
    a->AddAt(new TH1D(GetNameHist("hNTriggers"                   ,i,j),"Number of triggers in this bin filled."          ,1   ,0      ,1                                                        	    ),GetNumberHist(kHistNTriggers              ,i,j));
    a->AddAt(new TH2D(GetNameHist("hDeltaPhiVsDeltaEta2p"	 ,i,j),"#Delta#Phi vs #Delta#eta"			 ,nbinseta ,-2.0*fAcceptanceCut,2.0*fAcceptanceCut,nbinsphi,-0.5*Pii,1.5*Pii			       ),GetNumberHist(khPhiEta	   ,i,j));//"classical" 2 particle correlation
    a->AddAt(new TH3F(GetNameHist("hDeltaPhiVsDeltaPhiVsDeltaEta",i,j),"#Delta#Phi_1 vs #Delta#Phi_2 vs #Delta#eta_{12}" ,nbinseta ,-2.0*fAcceptanceCut,2.0*fAcceptanceCut,nbinsphi,-0.5*Pii,1.5*Pii,nbinsphi ,-0.5*Pii,1.5*Pii),GetNumberHist(khPhiPhiDEta,i,j));//3d, DPhiDPhiDEta
    }
  }
  a->AddAt(new TH1D("overflow","overflow",3,0.0,1),GetNumberHist(kNofHistograms+1,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1));
  TH1::AddDirectory(statusAddDirectory);
  return 0;
}

int AliCorrelation3p_noQA::SetMultVZ(Double_t Mult, Double_t Vz)
{
  fMultiplicity=Mult;
  fVZ=Vz;
  fMBin = GetMultBin(fMultiplicity);
  fVzBin = GetZBin(fVZ);
  HistFill(kcentrvsvz,fMultiplicity,fVZ);
  HistFill(kcentrvsvzbin,fMBin,fVzBin);
  if (fMBin<0||fVzBin<0) return -1;
  return 1;
}

bool AliCorrelation3p_noQA::CheckTrigger( AliVParticle* ptrigger, bool doHistogram)
{
  // check trigger particle cuts
  if (!ptrigger) return false;
  if (ptrigger->Pt()<=fMinTriggerPt) return false;
  if (fMaxTriggerPt>fMinTriggerPt && ptrigger->Pt()>fMaxTriggerPt) return false;
  return true;
}

bool AliCorrelation3p_noQA::CheckAssociated( AliVParticle* p, bool doHistogram)
{
  // check associated particle cuts
  if (!p) return false;
  if (p->Pt()<=fMinAssociatedPt) return false;
  if (fMaxAssociatedPt>fMinAssociatedPt && p->Pt()>fMaxAssociatedPt) return false;
  return true;
}

int AliCorrelation3p_noQA::Fill( AliVParticle* ptrigger,  AliVParticle* p1,  AliVParticle* p2, const double weight)
{
  /// fill histograms from particles, fills each histogram exactly once.
//   Double_t fillweight = 1.0;
  if (!ptrigger || !p1 || !p2) {return 0;}
  if ((ptrigger->Pt()<=p1->Pt())||(ptrigger->Pt()<=p2->Pt())) {return 0;}
  const double Pii=TMath::Pi();
  // phi difference associated 1 to trigger particle
  Double_t DeltaPhi1 = ptrigger->Phi() - p1->Phi();
  if(DeltaPhi1<-0.5*Pii||DeltaPhi1>1.5*Pii){
    if (DeltaPhi1<-0.5*Pii) DeltaPhi1 += 2*Pii;
    if (DeltaPhi1>1.5*Pii)  DeltaPhi1 -= 2*Pii;
  }
  // phi difference associated 2 to trigger particle
  Double_t DeltaPhi2 = ptrigger->Phi() - p2->Phi();
  if(DeltaPhi2<-0.5*Pii||DeltaPhi2>1.5*Pii){
    if (DeltaPhi2<-0.5*Pii) DeltaPhi2 += 2*Pii;
    if (DeltaPhi2>1.5*Pii)  DeltaPhi2 -= 2*Pii;
  }
  // eta difference
  Double_t DeltaEta12 = p1->Eta()-p2->Eta();
  if(TMath::Abs(DeltaPhi1-DeltaPhi2)<1.0E-10&&TMath::Abs(DeltaEta12)<1.0E-10)	return 0;//Track duplicate, reject.
  if(TMath::Abs(p1->Eta()-ptrigger->Eta())<1.0E-10)			return 0;//Track duplicate, reject.
  if(TMath::Abs(p2->Eta()-ptrigger->Eta())<1.0E-10)			return 0;//Track duplicate, reject.
//   fillweight *= dynamic_cast<AliFilteredTrack*>(ptrigger)->GetEff();
//   fillweight *= dynamic_cast<AliFilteredTrack*>(p1)->GetEff();
//   fillweight *= dynamic_cast<AliFilteredTrack*>(p2)->GetEff();
//   fillweight *= weight;
  HistFill(GetNumberHist(khPhiPhiDEta,fMBin,fVzBin),DeltaEta12,DeltaPhi1,DeltaPhi2, weight);
  return 0;
}

int AliCorrelation3p_noQA::Fill(AliVParticle* ptrigger,AliVParticle* p1, const double weight)
{
  /// fill histograms from particles
//   Double_t fillweight = 1.0;
  if (!ptrigger || !p1) return -EINVAL;
  if (ptrigger->Pt()<=p1->Pt()) return 0;
  const double Pii=TMath::Pi();
  // phi difference associated  to trigger particle
  Double_t DeltaPhi = ptrigger->Phi() - p1->Phi();
  if(DeltaPhi<-0.5*Pii||DeltaPhi>1.5*Pii){
    if (DeltaPhi<-0.5*Pii) DeltaPhi += 2*Pii;
    if (DeltaPhi>1.5*Pii)  DeltaPhi -= 2*Pii;
  }
  // eta difference
  Double_t DeltaEta  = ptrigger->Eta() - p1->Eta();
//   fillweight *= dynamic_cast<AliFilteredTrack*>(ptrigger)->GetEff();
//   fillweight *= dynamic_cast<AliFilteredTrack*>(p1)->GetEff();
  HistFill(GetNumberHist(khPhiEta,fMBin,fVzBin),DeltaEta,DeltaPhi, weight);//2p correlation
  return 0;
}

int AliCorrelation3p_noQA::Filla(AliVParticle* p1,AliVParticle* p2,const double weight)
{
  /// fill histograms from particles
  if (!p2 || !p1) return -EINVAL;
  // phi difference associated  to trigger particle
  Double_t DeltaPhi = p1->Phi() - p2->Phi();
  if (DeltaPhi<-0.5*gkPii) DeltaPhi += 2*gkPii;
  if (DeltaPhi>1.5*gkPii)  DeltaPhi -= 2*gkPii;
  // eta difference
  Double_t DeltaEta  = p1->Eta() - p2->Eta();
  HistFill(GetNumberHist(khPhiEtaa,fMBin,fVzBin),DeltaEta,DeltaPhi, weight);//2p correlation
  return 0;
}

int AliCorrelation3p_noQA::FillTrigger(AliVParticle* ptrigger)
{
  Double_t fillweight = dynamic_cast<AliFilteredTrack*>(ptrigger)->GetEff();
  HistFill(GetNumberHist(kHistNTriggers,fMBin,fVzBin),0.5,fillweight);//Increments number of triggers by weight. Call before filling with any associated.
  return 1;
}

void AliCorrelation3p_noQA::Clear(Option_t * /*option*/)
{
  /// overloaded from TObject: cleanup
  return TObject::Clear();
}

void AliCorrelation3p_noQA::Print(Option_t *option) const
{
  /// overloaded from TObject: print info
  const char* key=NULL;
  bool bNoRuler=false;
  const TString delimiter(" ");
  TStringToken token(option, delimiter);
  while (token.NextToken()) {
    key="noruler";
    if (token.CompareTo( key)==0) {
      bNoRuler=true;
      continue;
    }
    AliWarning(Form("unknown option '%s'" ,token.Data()));
  }
  if (!bNoRuler)
    AliWarning("====================================================================" );
  TObject::Print();
  if (fHistograms) {
    fHistograms->Print();
  }
  if (fMixedEvent) {
    AliWarning("  ---- mixed event -------------------------------------------------" );
    fMixedEvent->Print(Form("%s noruler", option));
  }
}

TObject* AliCorrelation3p_noQA::FindObject(const char *name) const
{
  /// overloaded from TObject: find object by name
  if (fHistograms) {
    TObject* o=fHistograms->FindObject(name);
    if (o) return o;
  }
  return NULL;
}

TObject* AliCorrelation3p_noQA::FindObject(const TObject *obj) const
{
  /// overloaded from TObject: find object by pointer
  if (fHistograms) {
    TObject* o=fHistograms->FindObject(obj);
    if (o) return o;
  }
  return NULL;
}

void AliCorrelation3p_noQA::SaveAs(const char *filename,Option_t */*option*/) const
{
  /// overloaded from TObject: save to file
  std::unique_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fHistograms) fHistograms->Write();
  output->Close();
}

AliCorrelation3p_noQA& AliCorrelation3p_noQA::operator+=(const AliCorrelation3p_noQA& other)
{
  /// add histograms from another instance
  if (!fHistograms || !other.fHistograms) return *this;
  for (int i=0; i<GetNumberHist(kNofHistograms+2,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1); i++) {
    if (fHistograms->At(i)==NULL || other.fHistograms->At(i)==NULL) continue;
    TH1* target=reinterpret_cast<TH1*>(fHistograms->At(i));
    TH1* source=reinterpret_cast<TH1*>(other.fHistograms->At(i));
    if (!target || !source) continue;
    TString name(fHistograms->At(i)->GetName());
    if (name.CompareTo(target->GetName())!=0) {
      AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", i, source->GetName(), target->GetName()));
      continue;
    }
    if (source->IsA()!=target->IsA()) {
      AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", i, source->ClassName(), target->ClassName()));
      continue;
    }
    target->Add(source);
  }
  return *this;
}
int AliCorrelation3p_noQA::Merge(TCollection* li)
{
   // interface for the TFileMerger
   TIter next(li);
   while (TObject* o = next()) {
      AliCorrelation3p_noQA* src=dynamic_cast<AliCorrelation3p_noQA*>(o);
      if (!src) {
	AliWarning(Form("skipping incompatible object %s of type %s", o->ClassName(), o->GetName()));
	continue;
      }
      (*this)+=(*src);
   }
   return 0;
}

const char* AliCorrelation3p_noQA::GetNameHist(const char* name, Int_t MBin, Int_t ZBin) const
{
  if((MBin>(fMBinEdges.GetSize()-1))||(ZBin>(fZBinEdges.GetSize()-1)))return name;//Index out of bounds
  Double_t Mlow = fMBinEdges.At(MBin);
  Double_t Mhigh = fMBinEdges.At(MBin+1);
  Double_t Zlow = fZBinEdges.At(ZBin);
  Double_t Zhigh = fZBinEdges.At(ZBin+1);
  return Form("%sM(%4.2f)->(%4.2f)Z(%4.2f)->(%4.2f)",name,Mlow,Mhigh,Zlow,Zhigh);
}

Int_t AliCorrelation3p_noQA::GetNumberHist(Int_t khist, Int_t Mbin, Int_t ZBin) const
{
  if((Mbin>(fMBinEdges.GetSize()-1))||(ZBin>(fZBinEdges.GetSize()-1)))return -1;//Index out of bounds
  int offset = kNonbinnedhists;//Offset so the extra histograms come first.
  return offset+khist+(Mbin+ZBin*(fMBinEdges.GetSize()))*kNofHistograms;//increase the number by number of hists for each M+zbin.
}

Int_t AliCorrelation3p_noQA::GetMultBin(Double_t Mult)
{
  int Nbins = fMBinEdges.GetSize()-1;
  for (int bin =0;bin<Nbins;bin++)if(Mult>=fMBinEdges[bin]&&Mult<fMBinEdges[bin+1])return bin;
  //if nothing happened in the loop, return -1:
  return -1;
}

Int_t AliCorrelation3p_noQA::GetZBin(Double_t Zvert)
{
  int Nbins = fZBinEdges.GetSize()-1;
  for (int bin =0;bin<Nbins;bin++)if(Zvert>=fZBinEdges[bin]&&Zvert<fZBinEdges[bin+1])return bin;
  //if nothing happened in the loop, return -1:
  return -1;
}

void AliCorrelation3p_noQA::HistFill(Int_t Histn,Double_t Val1)
{
  if(Histn>=0)dynamic_cast<TH1D*>(fHistograms->At(Histn))->Fill(Val1);
  else   dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms+1,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1)))->Fill(0.5);//wrong histn
}
void AliCorrelation3p_noQA::HistFill(Int_t Histn,Double_t Val1, Double_t Val2)
{
  if(Histn>=0){
    if(dynamic_cast<TH2D*>(fHistograms->At(Histn)))dynamic_cast<TH2D*>(fHistograms->At(Histn))->Fill(Val1,Val2);
    else if (dynamic_cast<TH1D*>(fHistograms->At(Histn)))dynamic_cast<TH1D*>(fHistograms->At(Histn))->Fill(Val1,Val2);
  }
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}
void AliCorrelation3p_noQA::HistFill(Int_t Histn,Double_t Val1, Double_t Val2, Double_t Val3)
{
  if(Histn>=0){
    if(dynamic_cast<TH3F*>(fHistograms->At(Histn)))dynamic_cast<TH3F*>(fHistograms->At(Histn))->Fill(Val1,Val2, Val3);
    else if(dynamic_cast<TH2D*>(fHistograms->At(Histn)))dynamic_cast<TH2D*>(fHistograms->At(Histn))->Fill(Val1,Val2, Val3);
  }
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}
void AliCorrelation3p_noQA::HistFill(Int_t Histn,Double_t Val1, Double_t Val2, Double_t Val3, Double_t weight)
{
  if(Histn>=0)dynamic_cast<TH3F*>(fHistograms->At(Histn))->Fill(Val1,Val2,Val3, weight);
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}

TH2D* AliCorrelation3p_noQA::slice(TH3F* hist, const char* option, Int_t firstbin, Int_t lastbin, const char* name, Bool_t baverage) const
{//option should be xy,zy,yx,zx,xz or yz.
  TString o = TString(option);
  TString namestring = TString(name);
  TH2D* Slice = new TH2D();
  if(o.CompareTo("xy")==0||o.CompareTo("yx")==0){
    if(o.CompareTo("xy")==0){Slice = (TH2D*)hist->Project3D(Form("%s_xy",name));Slice->Reset("m");}
    if(o.CompareTo("yx")==0){Slice = (TH2D*)hist->Project3D(Form("%s_yx",name));Slice->Reset("m");}
    for(int x=0; x<=hist->GetNbinsX()+1;x++){
      for(int y=0;y<=hist->GetNbinsY()+1;y++){
	Double_t Content=0;
	Double_t locerr=0;
	for(int z=firstbin;z<=lastbin;z++){
	  if(!baverage){
	  Content += hist->GetBinContent(x,y,z);
	  locerr  += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr !=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr += 1.0/(binerr*binerr);
	    }
	    if(binerr ==0){
	     locerr +=1.0;
	    }
	  }
	}
	if(baverage&&locerr!=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(o.CompareTo("xy")==0) Slice->SetBinContent(y,x,Content);
	if(o.CompareTo("yx")==0) Slice->SetBinContent(x,y,Content);
	if(o.CompareTo("xy")==0) Slice->SetBinError(y,x,TMath::Sqrt(locerr));
	if(o.CompareTo("yx")==0) Slice->SetBinError(x,y,TMath::Sqrt(locerr));

      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
  }
  if(o.CompareTo("xz")==0||o.CompareTo("zx")==0){
    if(o.CompareTo("xz")==0){Slice = (TH2D*)hist->Project3D(Form("%s_xz",name));Slice->Reset("m");}
    if(o.CompareTo("zx")==0){Slice = (TH2D*)hist->Project3D(Form("%s_zx",name));Slice->Reset("m");}
    for(int x=0; x<=hist->GetNbinsX()+1;x++){
      for(int z=0;z<=hist->GetNbinsZ()+1;z++){
	Double_t Content=0;
	Double_t locerr = 0;
	for(int y=firstbin;y<lastbin;y++){
	  if(!baverage){
	  Content += hist->GetBinContent(x,y,z);
	  locerr   += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr!=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr += 1.0/(binerr*binerr);
	    }
	  }
	}
	if(baverage&&locerr!=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(baverage&&locerr==0){
	  Content = 0;
	  locerr = 0;
	}
	if(o.CompareTo("xz")==0) Slice->SetBinContent(z,x,Content);
	if(o.CompareTo("zx")==0) Slice->SetBinContent(x,z,Content);
	if(o.CompareTo("xz")==0) Slice->SetBinError(z,x,TMath::Sqrt(locerr));
	if(o.CompareTo("zx")==0) Slice->SetBinError(x,z,TMath::Sqrt(locerr));
      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
    Double_t stats[6];
    Slice->GetStats(stats);
    Slice->PutStats(stats);
  }
  if(o.CompareTo("yz")==0||o.CompareTo("zy")==0){
    if(o.CompareTo("yz")==0){Slice = (TH2D*)hist->Project3D(Form("%s_yz",name));Slice->Reset("m");}
    if(o.CompareTo("zy")==0){Slice = (TH2D*)hist->Project3D(Form("%s_zy",name));Slice->Reset("m");}

    for(int y=0; y<=hist->GetNbinsY()+1;y++){
      for(int z=0;z<=hist->GetNbinsZ()+1;z++){
	Double_t Content=0;
	Double_t locerr  =0;
	for(int x=firstbin;x<=lastbin;x++){
	  if(!baverage){
	    Content += hist->GetBinContent(x,y,z);
	    locerr  += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr !=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr  += 1.0/(binerr*binerr);
	    }
	  }
	}
	if(baverage&&locerr !=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(baverage&&locerr==0){
	  Content = 0;
	  locerr = 1.0;
	}
	if((o.CompareTo("yz")==0)&&!(Content==0)) Slice->SetBinContent(z,y,Content);
	if((o.CompareTo("zy")==0)&&!(Content==0)) Slice->SetBinContent(y,z,Content);
	if((o.CompareTo("yz")==0)&&!(Content==0)) Slice->SetBinError(z,y,TMath::Sqrt(locerr));
	if((o.CompareTo("zy")==0)&&!(Content==0)) Slice->SetBinError(y,z,TMath::Sqrt(locerr));
      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
    Double_t stats[6];
    Slice->GetStats(stats);
    Slice->PutStats(stats);
  }
  return Slice;
}

void AliCorrelation3p_noQA::AddSlice(TH3F* hist,TH2D* AddTo, const char* option, Int_t firstbin, Int_t lastbin, const char* name, Bool_t baverage) const
{//option should be xy,zy,yx,zx,xz or yz.
  TString o = TString(option);
  TString namestring = TString(name);
  TH2D* Slice = new TH2D();
  if(o.CompareTo("xy")==0||o.CompareTo("yx")==0){
    if(o.CompareTo("xy")==0){Slice = (TH2D*)hist->Project3D(Form("%s_xy",name));Slice->Reset("m");}
    if(o.CompareTo("yx")==0){Slice = (TH2D*)hist->Project3D(Form("%s_yx",name));Slice->Reset("m");}
    for(int x=0; x<=hist->GetNbinsX()+1;x++){
      for(int y=0;y<=hist->GetNbinsY()+1;y++){
	Double_t Content=0;
	Double_t locerr=0;
	for(int z=firstbin;z<=lastbin;z++){
	  if(!baverage){
	  Content += hist->GetBinContent(x,y,z);
	  locerr  += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr !=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr += 1.0/(binerr*binerr);
	    }
	    if(binerr ==0){
	     locerr +=1.0;
	    }
	  }
	}
	if(baverage&&locerr!=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(o.CompareTo("xy")==0) Slice->SetBinContent(y,x,Content);
	if(o.CompareTo("yx")==0) Slice->SetBinContent(x,y,Content);
	if(o.CompareTo("xy")==0) Slice->SetBinError(y,x,TMath::Sqrt(locerr));
	if(o.CompareTo("yx")==0) Slice->SetBinError(x,y,TMath::Sqrt(locerr));

      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
  }
  if(o.CompareTo("xz")==0||o.CompareTo("zx")==0){
    if(o.CompareTo("xz")==0){Slice = (TH2D*)hist->Project3D(Form("%s_xz",name));Slice->Reset("m");}
    if(o.CompareTo("zx")==0){Slice = (TH2D*)hist->Project3D(Form("%s_zx",name));Slice->Reset("m");}
    for(int x=0; x<=hist->GetNbinsX()+1;x++){
      for(int z=0;z<=hist->GetNbinsZ()+1;z++){
	Double_t Content=0;
	Double_t locerr = 0;
	for(int y=firstbin;y<lastbin;y++){
	  if(!baverage){
	  Content += hist->GetBinContent(x,y,z);
	  locerr   += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr!=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr += 1.0/(binerr*binerr);
	    }
	  }
	}
	if(baverage&&locerr!=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(baverage&&locerr==0){
	  Content = 0;
	  locerr = 0;
	}
	if(o.CompareTo("xz")==0) Slice->SetBinContent(z,x,Content);
	if(o.CompareTo("zx")==0) Slice->SetBinContent(x,z,Content);
	if(o.CompareTo("xz")==0) Slice->SetBinError(z,x,TMath::Sqrt(locerr));
	if(o.CompareTo("zx")==0) Slice->SetBinError(x,z,TMath::Sqrt(locerr));
      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
    Double_t stats[6];
    Slice->GetStats(stats);
    Slice->PutStats(stats);
  }
  if(o.CompareTo("yz")==0||o.CompareTo("zy")==0){
    if(o.CompareTo("yz")==0){Slice = (TH2D*)hist->Project3D(Form("%s_yz",name));Slice->Reset("m");}
    if(o.CompareTo("zy")==0){Slice = (TH2D*)hist->Project3D(Form("%s_zy",name));Slice->Reset("m");}

    for(int y=0; y<=hist->GetNbinsY()+1;y++){
      for(int z=0;z<=hist->GetNbinsZ()+1;z++){
	Double_t Content=0;
	Double_t locerr  =0;
	for(int x=firstbin;x<=lastbin;x++){
	  if(!baverage){
	    Content += hist->GetBinContent(x,y,z);
	    locerr  += hist->GetBinError(x,y,z)*hist->GetBinError(x,y,z);
	  }
	  if(baverage){//Get the average
	    Double_t binerr = hist->GetBinError(x,y,z);
	    if(binerr !=0){
	      Content += hist->GetBinContent(x,y,z)/(binerr*binerr);
	      locerr  += 1.0/(binerr*binerr);
	    }
	  }
	}
	if(baverage&&locerr !=0){
	  Content = Content/locerr;//normalize
	  locerr = 1.0/locerr;
	}
	if(baverage&&locerr==0){
	  Content = 0;
	  locerr = 1.0;
	}
	if((o.CompareTo("yz")==0)&&!(Content==0)) Slice->SetBinContent(z,y,Content);
	if((o.CompareTo("zy")==0)&&!(Content==0)) Slice->SetBinContent(y,z,Content);
	if((o.CompareTo("yz")==0)&&!(Content==0)) Slice->SetBinError(z,y,TMath::Sqrt(locerr));
	if((o.CompareTo("zy")==0)&&!(Content==0)) Slice->SetBinError(y,z,TMath::Sqrt(locerr));
      }
    }
    Slice->SetEntries(Slice->GetEffectiveEntries());
    Double_t stats[6];
    Slice->GetStats(stats);
    Slice->PutStats(stats);
  }
  AddTo->Add(Slice);
  delete Slice;
}

TH2D* AliCorrelation3p_noQA::DetaDphiAss(TH3F* hist, const char* name)
{
  Double_t Pii = TMath::Pi();
  TH2D * DetaDphiAss = (TH2D*)hist->Project3D(Form("%s_yx",name));DetaDphiAss->Reset("m");
  for(int x = 1;x<hist->GetNbinsX();x++){for(int y=1;y<hist->GetNbinsY();y++){for(int z=1;z<hist->GetNbinsZ();z++){
    Double_t content = hist->GetBinContent(x,y,z);
    Double_t DPhi1 = hist->GetYaxis()->GetBinCenter(y);
    Double_t DPhi2 = hist->GetZaxis()->GetBinCenter(z);
    Double_t DPhi12 = DPhi1 - DPhi2;
    while(DPhi12<-0.5*Pii||DPhi12>1.5*Pii){
      if (DPhi12<-0.5*Pii) DPhi12 += 2*Pii;
      if (DPhi12>1.5*Pii)  DPhi12 -= 2*Pii;
    }
    Int_t Dphibin = DetaDphiAss->GetYaxis()->FindBin(DPhi12);
    content += DetaDphiAss->GetBinContent(x,Dphibin);
    DetaDphiAss->SetBinContent(x,Dphibin,content);
  }}}
  for(int x = 1;x<DetaDphiAss->GetNbinsX();x++){for(int y=1;y<DetaDphiAss->GetNbinsY();y++){
    Double_t content = DetaDphiAss->GetBinContent(x,y);
    if(content<0.5){continue;}
    Double_t error = 1.0/TMath::Sqrt(content);
    DetaDphiAss->SetBinError(x,y,error);
  }}
  return DetaDphiAss;
}


TH2D* AliCorrelation3p_noQA::DeltaEtaCut(TH3F* hist, const char* option, const char* name, Bool_t baverage) const
{
  //option can be: sameside    = if deltaPhi_1<pi/2, deltaPhi_2<pi/2 and  if deltaPhi_1>pi/2, deltaPhi_2>pi/2
  //		   lesspi2     = DeltaPhi_[12]<pi/2
  //		   lesspi4     = DeltaPhi_[12]<pi/4
  TString o = TString(option);
  TAxis* dphi1axis=hist->GetYaxis();
  TAxis* dphi2axis=hist->GetZaxis();
  TH2D* Result = (TH2D*)hist->Project3D(Form("%s_yx",name));
  Result->Reset("m");
  for(int deta12 = 0;deta12<=hist->GetNbinsX()+1;deta12++){
    for(int dPhi1 = 0;dPhi1<=hist->GetNbinsY()+1;dPhi1++){
      Double_t Content = 0.0;
      Double_t errorloc = 0.0;
      Bool_t bissameside=kFALSE;
      Bool_t bislesspi2=kFALSE;
      Bool_t bislesspi4=kFALSE;
      for(int dPhi2=1;dPhi2<=hist->GetNbinsZ();dPhi2++){
	bissameside = o.CompareTo("sameside");
	bissameside=bissameside&&(((dphi1axis->GetBinCenter(dPhi1) <TMath::Pi()/2.0)&&(dphi2axis->GetBinCenter(dPhi2) <TMath::Pi()/2.0))||((dphi1axis->GetBinCenter(dPhi1) >=TMath::Pi()/2.0)&&(dphi2axis->GetBinCenter(dPhi2) >=TMath::Pi()/2.0)));
	bislesspi2 = o.CompareTo("lesspi2");
	bislesspi2 = bislesspi2&&(TMath::Abs(dphi1axis->GetBinCenter(dPhi1)-dphi2axis->GetBinCenter(dPhi2)) <TMath::Pi()/2.0);
	bislesspi4 = o.CompareTo("lesspi4");
	bislesspi4 = bislesspi4&&(TMath::Abs(dphi1axis->GetBinCenter(dPhi1)-dphi2axis->GetBinCenter(dPhi2)) <TMath::Pi()/4.0);
	if(bissameside||bislesspi2||bislesspi4){
	  if(!baverage){
	    Content += hist->GetBinContent(deta12,dPhi1,dPhi2);
	    errorloc+= hist->GetBinError(deta12,dPhi1,dPhi2)*hist->GetBinError(deta12,dPhi1,dPhi2);
	  }
	  if(baverage){
	    Double_t binerror = hist->GetBinError(deta12,dPhi1,dPhi2);
	    if(binerror>0){
	      Content += hist->GetBinContent(deta12,dPhi1,dPhi2)/(binerror*binerror);
	      errorloc += 1/(binerror*binerror);
	    }
	    else{
	      errorloc += 1;
	    }
	  }
	}
      }
      if(baverage&&errorloc!=0){
	Content = Content/errorloc;
	errorloc = 1.0/errorloc;
      }
      Result->SetBinContent(deta12,dPhi1,Content);
      Result->SetBinError(deta12,dPhi1,TMath::Sqrt(errorloc));
    }
  }
  Result->SetEntries(Result->GetEffectiveEntries());
  return Result;
}

TCanvas* AliCorrelation3p_noQA::Makecanvas(TH1D* histtopl, TH1D* histtopm, TH1D* histtopr,TH1D* histmidl,TH1D* histmidm, TH1D* histmidr,TH1D* histbotl,TH1D* histbotm, TH1D* histbotr, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(3,3);
  Canvas->cd(1);
  gPad->SetLogy(1);
  if(!Stats) histtopl->SetStats(0);
  histtopl->Draw("E");
  Canvas->cd(2);
  gPad->SetLogy(0);
  if(!Stats) histtopm->SetStats(0);
  histtopm->Draw("E");
  Canvas->cd(3);
  gPad->SetLogy(0);
  if(!Stats) histtopr->SetStats(0);
  histtopr->Draw("E");
  Canvas->cd(4);
  gPad->SetLogy(1);
  if(!Stats) histmidl->SetStats(0);
  histmidl->Draw("E");
  Canvas->cd(5);
  gPad->SetLogy(0);
  if(!Stats) histmidm->SetStats(0);
  histmidm->Draw("E");
  Canvas->cd(6);
  gPad->SetLogy(0);
  if(!Stats) histmidr->SetStats(0);
  histmidr->Draw("E");
  Canvas->cd(7);
  gPad->SetLogy(1);
  if(!Stats) histbotl->SetStats(0);
  histbotl->Draw("E");
  Canvas->cd(8);
  gPad->SetLogy(0);
  if(!Stats) histbotm->SetStats(0);
  histbotm->Draw("E");
  Canvas->cd(9);
  gPad->SetLogy(0);
  if(!Stats) histbotr->SetStats(0);
  histbotr->Draw("E");
  return Canvas;
}

TCanvas* AliCorrelation3p_noQA::Makecanvas(TH2D* histtopl, TH2D* histtopr, TH2D* histbotl, TH2D* histbotr, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  Canvas->cd(1);
  if(!Stats)histtopl->SetStats(0);
  histtopl->Draw("surf3");
  Canvas->cd(2);
  if(!Stats)histtopr->SetStats(0);
  histtopr->Draw("surf3");
  Canvas->cd(3);
  if(!Stats)histbotl->SetStats(0);
  histbotl->Draw("surf3");
  Canvas->cd(4);
  if(!Stats)histbotr->SetStats(0);
  histbotr->Draw("surf3");
  return Canvas;
}

TCanvas* AliCorrelation3p_noQA::Makecanvas(TH2D* hist, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  if(!Stats)hist->SetStats(0);
  Canvas->cd(1);
  hist->Draw("surf2");
  Canvas->cd(2);
  hist->Draw("surf3");
  Canvas->cd(3);
  Int_t binpihn = hist->GetYaxis()->FindBin(0+0.2);
  Int_t binpiln = hist->GetYaxis()->FindBin(0-0.2);
  TH1D* projY = hist->ProjectionX(Form("%s%s",hist->GetName(),"_nearside"),binpiln,binpihn);
  if(!Stats) projY->SetStats(0);
  projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
  projY->SetTitle("Integral over the near side peak with #Delta#eta_{12} = 0#pm 0.2");
  projY->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projY->SetTitleSize(0.04,"xy");
  projY->SetTitleOffset(1.05,"xy");
  projY->Draw("E");
  Canvas->cd(4);
  Int_t binpih = hist->GetYaxis()->FindBin(TMath::Pi()+0.2);
  Int_t binpil = hist->GetYaxis()->FindBin(TMath::Pi()-0.2);
  TH1D* projX = hist->ProjectionX(Form("%s%s",hist->GetName(),"_px"),binpil,binpih);
  if(!Stats) projX->SetStats(0);
  projX->SetTitle("Integral over a slice of the away side peak around with #Delta#eta_{12} = #pi#pm 0.2");
  projX->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projX->SetTitleSize(0.04,"xy");
  projX->SetTitleOffset(1.05,"xy");
  projX->Draw("E");
  return Canvas;
}

Double_t AliCorrelation3p_noQA::FindScalingfactor(const char* scalingmethod,TH2D* sighist,TH2D*mixhist)
{
  TString selector=TString(scalingmethod);
  Double_t scalingfactor = 0;
  if(selector.CompareTo("ppdata")==0){
      Double_t mixedintegral = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegral = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegral!=0) scalingfactor=sigintegral/mixedintegral;
      if(mixedintegral==0 && sigintegral!=0) scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("ppdatascaled")==0){
      Double_t mixedintegralsc = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegralsc = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegralsc!=0) scalingfactor=sigintegralsc/mixedintegralsc;
      if(mixedintegralsc==0&&sigintegralsc!=0) scalingfactor=0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("ppdata2d")){
      Double_t mixedintegral2p = mixhist->Integral(0,mixhist->GetNbinsX(),0,4);
      Double_t sigintegral2p = sighist->Integral(0,sighist->GetNbinsX(),0,4);
      if(mixedintegral2p!=0) scalingfactor = sigintegral2p/mixedintegral2p;
      if(mixedintegral2p==0&&sigintegral2p!=0)  scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("PbPbdata")==0){
      Double_t mixedintegral = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegral = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegral!=0) scalingfactor=sigintegral/mixedintegral;
      if(mixedintegral==0 && sigintegral!=0) scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("PbPbdatascaled")==0){
      Double_t mixedintegralsc = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegralsc = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegralsc!=0) scalingfactor=sigintegralsc/mixedintegralsc;
      if(mixedintegralsc==0&&sigintegralsc!=0) scalingfactor=0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("PbPbdata2d")){
      Double_t mixedintegral2p = mixhist->Integral(0,mixhist->GetNbinsX(),0,4);
      Double_t sigintegral2p = sighist->Integral(0,sighist->GetNbinsX(),0,4);
      if(mixedintegral2p!=0) scalingfactor = sigintegral2p/mixedintegral2p;
      if(mixedintegral2p==0&&sigintegral2p!=0)  scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("generator")==0){
      Double_t mixedintegral = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegral = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegral!=0) scalingfactor=sigintegral/mixedintegral;
      if(mixedintegral==0 && sigintegral!=0) scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("generatorscaled")==0){
      Double_t mixedintegralsc = mixhist->Integral(0,8,mixhist->GetYaxis()->FindBin(4),mixhist->GetNbinsY()) + mixhist->Integral(mixhist->GetXaxis()->FindBin(4),mixhist->GetNbinsX(),0,8);
      Double_t sigintegralsc = sighist->Integral(0,8,sighist->GetYaxis()->FindBin(4),sighist->GetNbinsY()) + sighist->Integral(sighist->GetXaxis()->FindBin(4),sighist->GetNbinsX(),0,8);
      if(mixedintegralsc!=0) scalingfactor=sigintegralsc/mixedintegralsc;
      if(mixedintegralsc==0&&sigintegralsc!=0) scalingfactor=0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  if(selector.CompareTo("generator2d")){
      Double_t mixedintegral2p = mixhist->Integral(0,mixhist->GetNbinsX(),0,4);
      Double_t sigintegral2p = sighist->Integral(0,sighist->GetNbinsX(),0,4);
      if(mixedintegral2p!=0) scalingfactor = sigintegral2p/mixedintegral2p;
      if(mixedintegral2p==0&&sigintegral2p!=0)  scalingfactor = 0.0;
      if(mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0))>10)scalingfactor = 1.0/mixhist->GetBinContent(mixhist->GetXaxis()->FindBin(0.0), mixhist->GetYaxis()->FindBin(0.0));
      else scalingfactor=0;
      return scalingfactor;
  }
  return scalingfactor;
}

Double_t AliCorrelation3p_noQA::GetPoint(TH1* hist, Double_t xpoint, Double_t ypoint, Double_t zpoint)
{
  //Returns the content of the bin with the corresponding points.
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  Double_t bincontent =0;
  if(hist1d){
    bincontent = hist1d->GetBinContent(hist1d->GetXaxis()->FindBin(xpoint));
  }
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  if(hist2d){
    bincontent = hist2d->GetBinContent(hist2d->GetXaxis()->FindBin(xpoint), hist2d->GetYaxis()->FindBin(ypoint));
  }
  TH2F* hist2df = dynamic_cast<TH2F*>(hist);
  if(hist2df){
    bincontent = hist2df->GetBinContent(hist2df->GetXaxis()->FindBin(xpoint), hist2df->GetYaxis()->FindBin(ypoint));
  }
  TH3D* hist3d = dynamic_cast<TH3D*>(hist);
  if(hist3d){
    bincontent = hist3d->GetBinContent(hist3d->GetXaxis()->FindBin(xpoint), hist3d->GetYaxis()->FindBin(ypoint),hist3d->GetZaxis()->FindBin(zpoint));
  }
  TH3F* hist3df = dynamic_cast<TH3F*>(hist);
  if(hist3df){
    bincontent = hist3df->GetBinContent(hist3df->GetXaxis()->FindBin(xpoint), hist3df->GetYaxis()->FindBin(ypoint),hist3df->GetZaxis()->FindBin(zpoint));
  }
  return bincontent;
}


void AliCorrelation3p_noQA::AddHists(Bool_t isAverage, TH1* histtoadd, TH1* addedhist)
{
  TH1D* hist11d = dynamic_cast<TH1D*>(histtoadd);
  TH1D* hist21d = dynamic_cast<TH1D*>(addedhist);
  TH2D* hist12d = dynamic_cast<TH2D*>(histtoadd);
  TH2D* hist22d = dynamic_cast<TH2D*>(addedhist);
  TH3F* hist13d = dynamic_cast<TH3F*>(histtoadd);
  TH3F* hist23d = dynamic_cast<TH3F*>(addedhist);
  if(hist11d&&hist21d){
    Int_t nbinsx1 = hist11d->GetNbinsX();
    if(nbinsx1!=hist21d->GetNbinsX()){AliWarning("The histograms do not match! TH1D* with different x dimensions.");
    return ;}
    for(int x=0; x<=nbinsx1+1;x++){
      Double_t content1 = histtoadd->GetBinContent(x);
      Double_t error1   = histtoadd->GetBinError(x);
      Double_t content2 = addedhist->GetBinContent(x);
      Double_t error2   = addedhist->GetBinError(x);
      Double_t result = 0;
      Double_t resulte = 0;
      if(!isAverage){
	result = content1 + content2;
	resulte = TMath::Sqrt(error1*error1 + error2*error2);
      }
      if(isAverage){
	if(content1>1.0e-10&&content2>1.0e-10){
	  result = content1/(error1*error1) + content2/(error2*error2);
	  resulte = 1/(error1*error1) + 1/(error2*error2);
	}
	if(content2>1.0e-10&&content1<1.0e-10){
	  result = content2/(error2*error2);
	  resulte =  1/(error2*error2) + 1.0;
	}
	if(content1>1.0e-10&&content2<1.0e-10){
	  result = content1/(error1*error1);
	  resulte = 1/(error1*error1);
	}
	if(resulte!=0){
	  result = result/resulte;
	  resulte = 1.0/TMath::Sqrt(resulte);
	}
	else{
	  result = 0.0;
	}
      }
      if(result!=0.0)histtoadd->SetBinContent(x,result);
      if(result!=0.0)histtoadd->SetBinError(x,resulte);
    }
  }
  else if(hist12d&&hist22d){
    Int_t nbinsx1 = hist12d->GetNbinsX();
    if(nbinsx1!=hist22d->GetNbinsX()){AliWarning("The histograms do not match! TH2D* with different x dimensions.");
    return ;}
    Int_t nbinsy1 = hist12d->GetNbinsY();
    if(nbinsy1!=hist22d->GetNbinsY()){AliWarning("The histograms do not match! TH2D* with different y dimensions.");
    return ;}
    for(int x=0; x<=nbinsx1+1;x++){
      for(int y=0; y<=nbinsy1+1;y++){
	Double_t content1 = histtoadd->GetBinContent(x,y);
	Double_t error1   = histtoadd->GetBinError(x,y);
	Double_t content2 = addedhist->GetBinContent(x,y);
	Double_t error2   = addedhist->GetBinError(x,y);
	Double_t result = 0;
	Double_t resulte = 0;
	if(!isAverage){
	result = content1 + content2;
	resulte = TMath::Sqrt(error1*error1 + error2*error2);
      }
      if(isAverage){
	if(content1>1.0e-10&&content2>1.0e-10){
	  result = content1/(error1*error1) + content2/(error2*error2);
	  resulte = 1/(error1*error1) + 1/(error2*error2);
	}
	if(content2>1.0e-10&&content1<1.0e-10){
	  result = content2/(error2*error2);
	  resulte =  1/(error2*error2) + 1.0;
	}
	if(content1>1.0e-10&&content2<1.0e-10){
	  result = content1/(error1*error1);
	  resulte = 1/(error1*error1);
	}
	if(resulte!=0){
	  result = result/resulte;
	  resulte = 1.0/TMath::Sqrt(resulte);
	}
	else{
	  result = 0.0;
	}
      }
      if(result!=0.0)histtoadd->SetBinContent(x,y,result);
      if(result!=0.0)histtoadd->SetBinError(x,y,resulte);
        }
      }
    }
  else if(hist13d&&hist23d){
    Int_t nbinsx1 = hist13d->GetNbinsX();
    if(nbinsx1!=hist23d->GetNbinsX()){AliWarning("The histograms do not match! TH3D* with different x dimensions.");
    return ;}
    Int_t nbinsy1 = hist13d->GetNbinsY();
    if(nbinsy1!=hist23d->GetNbinsY()){AliWarning("The histograms do not match! TH3D* with different y dimensions.");
    return ;}
    Int_t nbinsz1 = hist13d->GetNbinsZ();
    if(nbinsz1!=hist23d->GetNbinsZ()){AliWarning("The histograms do not match! TH2D* with different z dimensions.");
    return ;}
    for(int x=0; x<=nbinsx1+1;x++){
      for(int y=0; y<=nbinsy1+1;y++){
	for(int z=0;z<=nbinsz1+1;z++){
	  Double_t content1 = histtoadd->GetBinContent(x,y,z);
	  Double_t error1   = histtoadd->GetBinError(x,y,z);
	  Double_t content2 = addedhist->GetBinContent(x,y,z);
	  Double_t error2   = addedhist->GetBinError(x,y,z);
	  Double_t result = 0;
	  Double_t resulte = 0;
	 if(!isAverage){
	result = content1 + content2;
	resulte = TMath::Sqrt(error1*error1 + error2*error2);
      }
      if(isAverage){
	if(content1>1.0e-10&&content2>1.0e-10){
	  result = content1/(error1*error1) + content2/(error2*error2);
	  resulte = 1/(error1*error1) + 1/(error2*error2);
	}
	if(content2>1.0e-10&&content1<1.0e-10){
	  result = content2/(error2*error2);
	  resulte =  1/(error2*error2) + 1.0;
	}
	if(content1>1.0e-10&&content2<1.0e-10){
	  result = content1/(error1*error1);
	  resulte = 1/(error1*error1);
	}
	if(resulte!=0){
	  result = result/resulte;
	  resulte = 1.0/TMath::Sqrt(resulte);
	}
	else{
	  result = 0.0;
	}
      }
      if(result!=0.0)histtoadd->SetBinContent(x,y,z,result);
      if(result!=0.0)histtoadd->SetBinError(x,y,z,resulte);
        }
      }
    }
  }
  else{
    AliWarning("The histograms do not match!");
    return;
  }
}

TH1 * AliCorrelation3p_noQA::PrepareHist(int HistLocation,const char* HistName,const char* title, const char* xaxis, const char* yaxis,const char* zaxis,bool mixed){
  TH1* Hist;
  if(!mixed)Hist= dynamic_cast<TH1*>(fHistograms->At(HistLocation)->Clone(HistName));
  if(mixed) Hist= dynamic_cast<TH1*>(fMixedEvent->fHistograms->At(HistLocation)->Clone(HistName));
  Hist->SetTitle(title);
  Hist->GetXaxis()->SetTitle(xaxis);
  Hist->GetYaxis()->SetTitle(yaxis);
  if(TString(zaxis).CompareTo("")!=0) Hist->GetZaxis()->SetTitle(zaxis);
  if(dynamic_cast<TH1D*>(Hist)){
    Hist->SetTitleSize(0.045,"xy");
    Hist->SetTitleOffset(1.2,"y");
  }
  if(dynamic_cast<TH2D*>(Hist)){
    Hist->SetTitleSize(0.045,"xyz");
    Hist->SetTitleOffset(1.2,"xy");
    Hist->SetTitleOffset(0.9,"z");
  }
  if(dynamic_cast<TH3F*>(Hist)){
    Hist->SetTitleSize(0.045,"xyz");
    Hist->SetTitleOffset(1.2,"xy");
    Hist->SetTitleOffset(0.9,"z");
  }
  return Hist;
}

TH1 * AliCorrelation3p_noQA::PrepareHist(TH1* Hist,const char* title, const char* xaxis, const char* yaxis,const char* zaxis, bool scale, TParameter<double> * par){
  Double_t LeastContent = 1.0;//Least amount in mixed to scale.
  Hist->SetTitle(title);
  Hist->GetXaxis()->SetTitle(xaxis);
  Hist->GetYaxis()->SetTitle(yaxis);
  if(TString(zaxis).CompareTo("")!=0) Hist->GetZaxis()->SetTitle(zaxis);
  if(dynamic_cast<TH1D*>(Hist)){
    Hist->SetTitleSize(0.045,"xy");
    Hist->SetTitleOffset(1.2,"y");
  }
  if(dynamic_cast<TH2D*>(Hist)){
    Hist->SetTitleSize(0.045,"xyz");
    Hist->SetTitleOffset(1.2,"xy");
    Hist->SetTitleOffset(0.9,"z");
    if(scale){
      if(par) par->SetVal(GetPoint(Hist,0.0,0.0));
      if(GetPoint(Hist,0.0,0.0)>LeastContent){Hist->Scale(1.0/GetPoint(Hist,0.0,0.0));}
      else Hist->Scale(0.0);
    }
  }
  if(dynamic_cast<TH3F*>(Hist)){
    Hist->SetTitleSize(0.045,"xyz");
    Hist->SetTitleOffset(1.2,"xy");
    Hist->SetTitleOffset(0.9,"z");
  }
  return Hist;
}

int AliCorrelation3p_noQA::MakeResultsFile(const char* scalingmethod, bool recreate,bool all)
{//This function makes a new file that contains all the histograms in all versions possible.
  TFile * outfile;
  TString dir = TString(scalingmethod);
  TDirectory * mixeddir=NULL;
  if(recreate){ outfile = new TFile("results.root","RECREATE");outfile->cd();AliWarning("Recreate");}
  else{
    outfile = new TFile("results.root","UPDATE");
  }
  if(TString("").CompareTo(scalingmethod)!=0){
    if(dir.Contains("/")){
      TObjArray *dirs=dir.Tokenize("/");
      if(dirs->GetEntries()>1){
	TDirectory * processdir;
	processdir = outfile->GetDirectory((dirs->At(0))->GetName());
	if(!processdir)processdir = outfile->mkdir((dirs->At(0))->GetName());
	processdir->cd();
	mixeddir = processdir->GetDirectory((dirs->At(1))->GetName());
	if(mixeddir){processdir->rmdir((dirs->At(1))->GetName());}
	mixeddir = processdir->mkdir((dirs->At(1))->GetName());
	mixeddir->cd();
      }
    }
    else{
      mixeddir = outfile->GetDirectory(dir.Data());
      if(mixeddir) outfile->rmdir(dir.Data());
      mixeddir = outfile->mkdir(dir.Data());
      dir.Append("/");
      mixeddir->cd();
    }
  }
  else{
    mixeddir = outfile->GetDirectory("");
    dir.Append("/");
    mixeddir->cd();
  }

  TDirectory * dirmzbin = NULL;
  TDirectory * dirmzbinsig=NULL;
  TDirectory * dirmzbinmixed=NULL;
  TDirectory * dirmzbindiv=NULL;
  TCanvas    * tempcanvas=NULL;
  //Hists and directories for the total stuff:
  TDirectory * totbinstats = gDirectory->mkdir("bin_stats");
  TH1D * HistpT = new TH1D();
  TH1D * HistPhi = new TH1D();
  TH1D * HistEta = new TH1D();
  TH1D * HistTriggerpT = new TH1D();
  TH1D * HistTriggerPhi = new TH1D();
  TH1D * HistTriggerEta = new TH1D();
  TH1D * HistAssociatedpT = new TH1D();
  TH1D * HistAssociatedPhi = new TH1D();
  TH1D * HistAssociatedEta = new TH1D();  Long_t NTriggers=0;
  Bool_t setAverage = kTRUE;
  for(int mb =0;mb<fMBinEdges.GetSize()-1;mb++){
    for(int zb=0;zb<fZBinEdges.GetSize()-1;zb++){
      if(mixeddir) mixeddir->cd();
      else outfile->cd();
      dirmzbin = gDirectory->mkdir(GetNameHist("Bin",mb,zb));
//       binstats = dirmzbin->mkdir("bin_stats");
      dirmzbinsig = dirmzbin->mkdir("same_event");
      dirmzbinmixed = dirmzbin->mkdir("mixed_event");
      dirmzbindiv = dirmzbin->mkdir("divided");
/////////same event histograms
	dirmzbinsig->cd();
	TH1D* scalinghist= dynamic_cast<TH1D*>(PrepareHist(GetNumberHist(kHistNTriggers,mb,zb),"number_of_triggers","Total number of triggers filled","","# triggers"));
	scalinghist->Write();
	TH3F* DPHIDPHIDETA = dynamic_cast<TH3F*>(PrepareHist(GetNumberHist(khPhiPhiDEta,mb,zb),"DPhi_1_DPhi_2_DEta_12","#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]"));
	DPHIDPHIDETA->Write();
	TH2D* DPhi12DEta12 = dynamic_cast<TH2D*>(PrepareHist(DetaDphiAss(DPHIDPHIDETA,"DPhi_12_DEta_12"),"#Delta#Phi_{12} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{12} [rad]"));
	DPhi12DEta12->Write();
/////////DPHIDPHI histograms:
	  TH2D* DPHIDPHI3 = slice(DPHIDPHIDETA,"yz",1,DPHIDPHIDETA->GetNbinsX(),"DPhi_1_DPHI_2",kFALSE);
	  PrepareHist(DPHIDPHI3,"#Delta#Phi_{1} vs #Delta#Phi_{2}","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","# Pairs");
	  DPHIDPHI3->Write("DPhi_1_DPHI_2");
	  TH2D* DPHIDPHI3near =slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-0.4),DPHIDPHIDETA->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near",kFALSE);
	  PrepareHist(DPHIDPHI3near,"#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","# Pairs");
	  DPHIDPHI3near->Write("DPhi_1_DPHI_2_near");
	  TH2D* DPHIDPHI3mid = slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-1),DPHIDPHIDETA->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid",kFALSE);
	  AddSlice(DPHIDPHIDETA,DPHIDPHI3mid,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETA->GetXaxis()->FindBin(1),"temphist1",kFALSE);
	  PrepareHist(DPHIDPHI3mid,"#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","# Pairs");
	  DPHIDPHI3mid->Write("DPhi_1_DPHI_2_mid");
	  TH2D* DPHIDPHI3far = slice(DPHIDPHIDETA,"yz",1,DPHIDPHIDETA->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far",kFALSE);
	  AddSlice(DPHIDPHIDETA,DPHIDPHI3far,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(1)+1,DPHIDPHIDETA->GetNbinsX(),"temphist2",kFALSE);
	  PrepareHist(DPHIDPHI3far,"#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","# Pairs");
	  DPHIDPHI3far->Write("DPhi_1_DPHI_2_far");
	  tempcanvas= Makecanvas(DPHIDPHI3,DPHIDPHI3near,DPHIDPHI3mid,DPHIDPHI3far,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
///////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3 = slice(DPHIDPHIDETA,"yx",1,DPHIDPHIDETA->GetNbinsZ(),"DPhi_1_DEta_12",kFALSE);
	  PrepareHist(DPHIDETA12_3,"#Delta#Phi_{1} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","# Pairs");
	  DPHIDETA12_3->Write("DPhi_1_DEta_12");
	  tempcanvas= Makecanvas(DPHIDETA12_3,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3 =DeltaEtaCut(DPHIDPHIDETA,"sameside","DPhi_1_DEta_12_SameSide",kFALSE);
	  PrepareHist(DPHIDETA12SameSide_3,"#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","# Pairs");
	  DPHIDETA12SameSide_3->Write("DPhi_1_DEta_12_SameSide");
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3,"DPHIDEta12_SameSide_3d",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA = dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(khPhiEta,mb,zb))->Clone("DPhi_DEta"));
	  PrepareHist(DPHIDETA,"#Delta#Phi vs #Delta#eta","#Delta#eta []","#Delta#Phi [rad]","# Associated");
	  DPHIDETA->Write();
	  tempcanvas= Makecanvas(DPHIDETA,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
/////////mixed event histograms
	dirmzbinmixed->cd();
	TParameter<double> * scale = new TParameter<double>("scale",0.0);
	TH1D* scalinghistm= dynamic_cast<TH1D*>(PrepareHist(GetNumberHist(kHistNTriggers,mb,zb),"number_of_triggers","Total number of times the mixed histogram was filled with a trigger","","# triggers","",true));
	scalinghistm->Write();
	TH3F* DPHIDPHIDETAm = dynamic_cast<TH3F*>(PrepareHist(GetNumberHist(khPhiPhiDEta,mb,zb),"DPhi_1_DPhi_2_DEta_12","#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]",true));
	DPHIDPHIDETAm->Write();
	TH2D* DPhi12DEta12m = DetaDphiAss(DPHIDPHIDETAm,"DPhi_12_DEta_12");
	PrepareHist(DPhi12DEta12m,"#Delta#Phi_{12} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{12} [rad]","",true);
	DPhi12DEta12m->Write("DPhi_12_DEta_12");
//////////DPHIDPHI histograms:
	  TH2D* DPHIDPHI3m = slice(DPHIDPHIDETAm,"yz",1,DPHIDPHIDETAm->GetNbinsX(),"DPhi_1_DPHI_2",kFALSE);
	  PrepareHist(DPHIDPHI3m,"#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","",true, scale);
	  DPHIDPHI3m->Write("DPhi_1_DPHI_2");
	  scale->Write("DPhi_1_DPHI_2_scale");
	  scale->SetVal(0.0);
	  TH2D* DPHIDPHI3nearm = slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4),DPHIDPHIDETAm->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near",kFALSE);
	  PrepareHist(DPHIDPHI3nearm,"#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","",true,scale);
	  DPHIDPHI3nearm->Write("DPhi_1_DPHI_2_near");
	  scale->Write("DPhi_1_DPHI_2_near_scale");
	  scale->SetVal(0.0);
	  TH2D* DPHIDPHI3midm = slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-1),DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid",kFALSE);
	  AddSlice(DPHIDPHIDETAm,DPHIDPHI3midm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETAm->GetXaxis()->FindBin(1),"temphist1",kFALSE);
	  PrepareHist(DPHIDPHI3midm,"#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","",true,scale);
	  DPHIDPHI3midm->Write("DPhi_1_DPHI_2_mid");
	  scale->Write("DPhi_1_DPHI_2_mid_scale");
	  scale->SetVal(0.0);
	  TH2D* DPHIDPHI3farm = slice(DPHIDPHIDETAm,"yz",1,DPHIDPHIDETAm->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far",kFALSE);
	  AddSlice(DPHIDPHIDETAm,DPHIDPHI3farm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(1)+1,DPHIDPHIDETAm->GetNbinsX(),"temphist2",kFALSE);
	  PrepareHist(DPHIDPHI3farm,"#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1","#Delta#Phi_{1} [rad]","#Delta#Phi_{2} [rad]","",true,scale);
	  DPHIDPHI3farm->Write("DPhi_1_DPHI_2_far");
	  scale->Write("DPhi_1_DPHI_2_far_scale");
	  scale->SetVal(0.0);
	  tempcanvas= Makecanvas(DPHIDPHI3m,DPHIDPHI3nearm,DPHIDPHI3midm,DPHIDPHI3farm,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
///////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3m = slice(DPHIDPHIDETAm,"yx",1,DPHIDPHIDETAm->GetNbinsZ(),"DPhi_1_DEta_12",kFALSE);
	  PrepareHist(DPHIDETA12_3m,"#Delta#Phi_{1} vs #Delta#eta_{12}","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","",true,scale);
	  DPHIDETA12_3m->Write("DPhi_1_DEta_12");
	  scale->Write("DPhi_1_DEta_12_scale");
	  scale->SetVal(0.0);
	  tempcanvas= Makecanvas(DPHIDETA12_3m,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3m =DeltaEtaCut(DPHIDPHIDETAm,"sameside","DPhi_1_DEta_12_SameSide",kFALSE);
	  PrepareHist(DPHIDETA12SameSide_3m,"#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side","#Delta#eta_{12} []","#Delta#Phi_{1} [rad]","",true,scale);
	  DPHIDETA12SameSide_3m->Write("DPhi_1_DEta_12_SameSide");
	  scale->Write("DPhi_1_DEta_12_SameSide_scale");
	  scale->SetVal(0.0);
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3m,"DPHIDEta12_SameSide",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETAm = dynamic_cast<TH2D*>(fMixedEvent->fHistograms->At(GetNumberHist(khPhiEta,mb,zb))->Clone("DPhi_DEta"));
	  PrepareHist(DPHIDETAm,"#Delta#Phi vs #Delta#eta","#Delta#eta []","#Delta#Phi [rad]","",true,scale);
	  DPHIDETAm->Write();
	  scale->Write("DPhi_DEta_scale");
	  scale->SetVal(0.0);
	  tempcanvas= Makecanvas(DPHIDETAm,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  delete scale;

      Double_t resultscalingfactor = 1.0;//Scale the result with 1/ntriggers
      if(scalinghist->Integral()!=0)resultscalingfactor=1.0/scalinghist->Integral();
      else resultscalingfactor=0.0;
      bool empty = false;
      if(GetPoint(DPHIDPHI3m,0.0,0.0)==0){empty = true;}
	NTriggers+=scalinghist->Integral();
	dirmzbindiv->cd();
	//the same hists as in same and mixed.
	TH1D* scalinghistdiv = (TH1D*)scalinghist->Clone("number_of_triggers");
	if(!empty)scalinghistdiv->Divide(scalinghistm);
	else scalinghistdiv->Scale(0.0);
	scalinghistdiv->Write();
	delete scalinghistdiv;
 	TH3F* DPHIDPHIDETAdiv = (TH3F*)DPHIDPHIDETA->Clone("DPhi_1_DPhi_2_DEta_12");
	if(!empty&&GetPoint(DPHIDPHIDETAm,0.0,0.0,0.0)!=0){DPHIDPHIDETAm->Scale(1.0/GetPoint(DPHIDPHIDETAm,0.0,0.0,0.0));DPHIDPHIDETAdiv->Divide(DPHIDPHIDETAm);}
	else DPHIDPHIDETAdiv->Scale(0.0);
	if(setAverage)DPHIDPHIDETAdiv->Scale(resultscalingfactor);
	DPHIDPHIDETAdiv->Write();
	delete DPHIDPHIDETAdiv;
	TH2D* DPhi12DEta12div = (TH2D*)DPhi12DEta12->Clone("DPhi_12_DEta_12");
	if(!empty)DPhi12DEta12div->Divide(DPhi12DEta12m);
	else DPhi12DEta12div->Scale(0.0);
	if(setAverage)DPhi12DEta12div->Scale(resultscalingfactor);
	DPhi12DEta12div->Write();
	delete DPhi12DEta12div;
	//DPHIDPHI histograms:
	  TH2D* DPHIDPHI3div = (TH2D*)DPHIDPHI3->Clone("DPhi_1_DPHI_2");
	  if(!empty)DPHIDPHI3div->Divide(DPHIDPHI3m);
	  else DPHIDPHI3div->Scale(0.0);
	  if(setAverage) DPHIDPHI3div->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDPHI3div->Scale(resultscalingfactor);
	  DPHIDPHI3div->Write();
	  TH2D* DPHIDPHI3neardiv = (TH2D*)DPHIDPHI3near->Clone("DPhi_1_DPHI_2_near");
	  if(!empty)DPHIDPHI3neardiv->Divide(DPHIDPHI3nearm);
	  else DPHIDPHI3neardiv->Scale(0.0);
	  if(setAverage) DPHIDPHI3neardiv->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDPHI3neardiv->Scale(resultscalingfactor);
	  DPHIDPHI3neardiv->Write();
	  TH2D* DPHIDPHI3middiv = (TH2D*)DPHIDPHI3mid->Clone("DPhi_1_DPHI_2_mid");
	  if(!empty)DPHIDPHI3middiv->Divide(DPHIDPHI3midm);
	  else DPHIDPHI3middiv->Scale(0.0);
	  if(setAverage) DPHIDPHI3middiv->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDPHI3middiv->Scale(resultscalingfactor);
	  DPHIDPHI3middiv->Write();
	  TH2D* DPHIDPHI3fardiv = (TH2D*)DPHIDPHI3far->Clone("DPhi_1_DPHI_2_far");
	  if(!empty)DPHIDPHI3fardiv->Divide(DPHIDPHI3farm);
	  else DPHIDPHI3fardiv->Scale(0.0);
	  if(setAverage) DPHIDPHI3fardiv->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDPHI3fardiv->Scale(resultscalingfactor);
	  DPHIDPHI3fardiv->Write();
	  tempcanvas= Makecanvas(DPHIDPHI3div,DPHIDPHI3neardiv,DPHIDPHI3middiv,DPHIDPHI3fardiv,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  delete DPHIDPHI3div;delete DPHIDPHI3neardiv;delete DPHIDPHI3middiv;delete DPHIDPHI3fardiv;
////////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3div = (TH2D*) DPHIDETA12_3->Clone("DPhi_1_DEta_12");
	  if(!empty)DPHIDETA12_3div->Divide(DPHIDETA12_3m);
	  else DPHIDETA12_3div->Scale(0.0);
	  if(setAverage) DPHIDETA12_3div->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDETA12_3div->Scale(resultscalingfactor);
	  DPHIDETA12_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3div,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  delete DPHIDETA12_3div;
	  TH2D* DPHIDETA12SameSide_3div = (TH2D*)DPHIDETA12SameSide_3->Clone("DPhi_1_DEta_12_SameSide");
	  if(!empty)DPHIDETA12SameSide_3div->Divide(DPHIDETA12SameSide_3m);
	  else DPHIDETA12SameSide_3div->Scale(0.0);
	  if(setAverage) DPHIDETA12SameSide_3div->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDETA12SameSide_3div->Scale(resultscalingfactor);
	  DPHIDETA12SameSide_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3div,"DPHIDEta12_SameSide",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  delete DPHIDETA12SameSide_3div;
	  TH2D* DPHIDETAdiv = (TH2D*)DPHIDETA->Clone("DPhi_DEta");
	  if(!empty)DPHIDETAdiv->Divide(DPHIDETAm);
	  else DPHIDETAdiv->Scale(0.0);
	  if(setAverage) DPHIDETAdiv->Scale(resultscalingfactor);
	  if(!setAverage)DPHIDETAdiv->Scale(resultscalingfactor);
	  DPHIDETAdiv->Write();
	  tempcanvas= Makecanvas(DPHIDETAdiv,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  delete DPHIDETAdiv;
      delete scalinghist;delete scalinghistm;
      delete DPHIDPHIDETA;delete DPHIDPHIDETAm;
      delete DPhi12DEta12;delete DPhi12DEta12m;
      delete DPHIDPHI3;delete DPHIDPHI3near;delete DPHIDPHI3mid;delete DPHIDPHI3far;delete DPHIDPHI3m;delete DPHIDPHI3nearm;delete DPHIDPHI3midm;delete DPHIDPHI3farm;
      delete DPHIDETA12_3;
      delete DPHIDETA12SameSide_3;
      delete DPHIDETA12_3m;
      delete DPHIDETA12SameSide_3m;
      delete DPHIDETA;delete DPHIDETAm;
    }
  }
  {//Binstats
    totbinstats->cd();
    HistpT->Write();
    HistPhi->Write();
    HistEta->Write();
    HistTriggerpT->Write();
    HistTriggerPhi->Write();
    HistTriggerEta->Write();
    HistAssociatedpT->Write();
    HistAssociatedPhi->Write();
    HistAssociatedEta->Write();
    tempcanvas = Makecanvas(HistpT,HistPhi,HistEta,HistTriggerpT,HistTriggerPhi,HistTriggerEta,HistAssociatedpT,HistAssociatedPhi,HistAssociatedEta,"samestatscanvas",kFALSE);
    tempcanvas->Write();
    delete tempcanvas;
    tempcanvas=NULL;
    delete HistpT;delete HistPhi;delete HistEta;delete HistTriggerpT;delete HistTriggerPhi;delete HistTriggerEta;delete HistAssociatedpT;delete HistAssociatedPhi;delete HistAssociatedEta;
  }
  outfile->Close();
  delete outfile;
  return 1;
}
