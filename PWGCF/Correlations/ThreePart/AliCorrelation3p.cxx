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

#include "AliCorrelation3p.h"
#include "AliVParticle.h"
#include "AliCFPI0.h"
#include "AliLog.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include <iostream>
#include <cerrno>
#include <memory>
#include <set>

using namespace std;

ClassImp(AliCorrelation3p)

AliCorrelation3p::AliCorrelation3p(const char* name,TArrayD MBinEdges, TArrayD ZBinEdges)
  : TNamed(name?name:"AliCorrelation3p", "")
  , fHistograms(NULL)
  , fMinTriggerPt(8.0)
  , fMaxTriggerPt(15.0)
  , fMinAssociatedPt(3.)
  , fMaxAssociatedPt(fMinTriggerPt)
  , fhPhiEtaDeltaPhi12Cut1(0.5*TMath::Pi())
  , fhPhiEtaDeltaPhi12Cut2(0.25*TMath::Pi())
  , fMixedEvent(NULL)
  , fMBinEdges(MBinEdges)
  , fZBinEdges(ZBinEdges)
  , fMultiplicity(0)
  , fVZ(0)
  , fMBin(0)
  , fVzBin(0)
  , fCollisionType(PbPb)
  , fTriggerType(tracks)
{
  // default constructor
}

AliCorrelation3p::AliCorrelation3p(const AliCorrelation3p& other)
  : TNamed(other)
  , fHistograms((other.fHistograms!=NULL?(new TObjArray(*other.fHistograms)):NULL))
  , fMinTriggerPt(other.fMinTriggerPt)
  , fMaxTriggerPt(other.fMaxTriggerPt)
  , fMinAssociatedPt(other.fMinAssociatedPt)
  , fMaxAssociatedPt(other.fMaxAssociatedPt)
  , fhPhiEtaDeltaPhi12Cut1(other.fhPhiEtaDeltaPhi12Cut1)
  , fhPhiEtaDeltaPhi12Cut2(other.fhPhiEtaDeltaPhi12Cut2)
  , fMixedEvent((other.fMixedEvent!=NULL?(new AliCorrelation3p(*other.fMixedEvent)):NULL))
  , fMBinEdges(other.fMBinEdges)
  , fZBinEdges(other.fZBinEdges)
  , fMultiplicity(other.fMultiplicity)
  , fVZ(other.fVZ)
  , fMBin(other.fMBin)
  , fVzBin(other.fVzBin)
  , fCollisionType(other.fCollisionType)
  , fTriggerType(other.fTriggerType)
{
  // copy constructor
}

AliCorrelation3p& AliCorrelation3p::operator=(const AliCorrelation3p& other)
{
  // assignment operator
  if (this==&other) return *this;
  this->~AliCorrelation3p();
  new (this) AliCorrelation3p(other);
  return *this;
}

AliCorrelation3p::~AliCorrelation3p()
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

void AliCorrelation3p::Copy(TObject &object) const
{
  /// overloaded from TObject: copy to target object
  AliCorrelation3p* target=dynamic_cast<AliCorrelation3p*>(&object);
  if (!target) return;

  AliCorrelation3p* backupME=fMixedEvent;
  if (!target->fMixedEvent) {
    // avoid copying the mixed event object if there is no target
    const_cast<AliCorrelation3p*>(this)->fMixedEvent=NULL;
  }
  *target=*this;
  const_cast<AliCorrelation3p*>(this)->fMixedEvent=backupME;
}

int AliCorrelation3p::Init(const char* arguments)
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
      cout << "Collision Type set to: "<<param<<endl;
      continue;
    }
     key="triggertype=";
    if (token.BeginsWith(key)) {
      TString param=token;
      param.ReplaceAll(key, "");
      if(param.CompareTo("tracks")==0) fTriggerType=tracks;
      if(param.CompareTo("pi0")==0)    fTriggerType=pi0;
      cout << "Trigger Type set to: "<<param<<endl;
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
  a->AddAt(new TH1D("hTpTc","Trigger pT",100,0.0,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt)),kHistpTTriggerallbins);
  a->AddAt(new TH1D("hTphic","Trigger phi",270,-.5*Pii ,2.5*Pii),kHistPhiTriggerallbins);
  a->AddAt(new TH1D("hTetac","Trigger eta",100,-3.0,3.0),kHistEtaTriggerallbins);
  a->AddAt(new TH1D("hApTc","Associated pT",100,0.0,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt)),kHistpTAssociatedallbins);
  a->AddAt(new TH1D("hAphic","Associated phi",270,-.5*Pii ,2.5*Pii),kHistPhiAssociatedallbins);
  a->AddAt(new TH1D("hAetac","Associated eta",100,-3.0,3.0),kHistEtaAssociatedallbins);
  //Create one set of histograms per mixing bin
  for(Int_t i=0;i<fMBinEdges.GetSize()-1;i++){
    for(Int_t j=0;j<fZBinEdges.GetSize()-1;j++){
    a->AddAt(new TH1D(GetNameHist("hpT"				     ,i,j),"pT"							    			,100,0.0     ,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt)),GetNumberHist(kHistpT				,i,j));
    a->AddAt(new TH1D(GetNameHist("hphi"			     ,i,j),"phi"						    			,270,-.5*Pii ,2.5*Pii						       ),GetNumberHist(kHistPhi				,i,j));
    a->AddAt(new TH1D(GetNameHist("heta"			     ,i,j),"eta"						   			,100,-3.0    ,3.0						       ),GetNumberHist(kHistEta				,i,j));
    a->AddAt(new TH1D(GetNameHist("hTpT"			     ,i,j),"Trigger pT"						   			,100,0.0     ,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt)),GetNumberHist(kHistTriggerpT			,i,j));
    a->AddAt(new TH1D(GetNameHist("hTphi"			     ,i,j),"Trigger phi"					   			,270,-.5*Pii ,2.5*Pii						       ),GetNumberHist(kHistTriggerPhi			,i,j));
    a->AddAt(new TH1D(GetNameHist("hTeta"			     ,i,j),"Trigger eta"					      			,100,-3.0    ,3.0						       ),GetNumberHist(kHistTriggerEta			,i,j));
    a->AddAt(new TH1D(GetNameHist("hApT"			     ,i,j),"Associated pT"					      			,100,0.0     ,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt)),GetNumberHist(kHistAssociatedpT		,i,j));
    a->AddAt(new TH1D(GetNameHist("hAphi"			     ,i,j),"Associated phi"					      			,270,-.5*Pii ,2.5*Pii						       ),GetNumberHist(kHistAssociatedPhi		,i,j));
    a->AddAt(new TH1D(GetNameHist("hAeta"			     ,i,j),"Associated eta"					      			,100,-3.0    ,3.0						       ),GetNumberHist(kHistAssociatedEta		,i,j));
    a->AddAt(new TH1D(GetNameHist("hNAssoc"			     ,i,j),"Number of associated"							,1000,0	     ,2000						       ),GetNumberHist(kHistNassoc			,i,j));
    a->AddAt(new TH1D(GetNameHist("hNTriggers"                       ,i,j),"Number of triggers in this bin filled."                                     ,1   ,0      ,1                                                        ),GetNumberHist(kHistNTriggers                   ,i,j));
    a->AddAt(new TH2D(GetNameHist("hDeltaPhiVsDeltaEta2p"	     ,i,j),"#Delta#Phi vs #Delta#eta"				     			,63 ,-1.6,1.6    ,36	    ,-0.5*Pii,1.5*Pii			       ),GetNumberHist(khPhiEta				,i,j));//"classical" 2 particle correlation
    a->AddAt(new TH3D(GetNameHist("hDeltaPhiVsDeltaPhiVsDeltaEta"    ,i,j),"#Delta#Phi_1 vs #Delta#Phi_2 vs #Delta#eta_{12}"                            ,63 ,-1.6,1.6    ,36        ,-0.5*Pii,1.5*Pii  ,36 ,-0.5*Pii,1.5*Pii   ),GetNumberHist(khPhiPhiDEta                     ,i,j));//3d, DPhiDPhiDEta
    a->AddAt(new TH3D(GetNameHist("hDeltaPhiVsDeltaPhiVsDEtaScaled"  ,i,j),"#Delta#Phi_1 vs #Delta#Phi_2 vs #Delta#eta_{12} scaled with # associated"   ,63 ,-1.6,1.6    ,36        ,-0.5*Pii,1.5*Pii  ,36 ,-0.5*Pii,1.5*Pii   ),GetNumberHist(khPhiPhiDEtaScaled               ,i,j));//3d, DPhiDPhiDEta scaled
    a->AddAt(new TH1D(GetNameHist("khQAtocheckadressing"             ,i,j),"Will be filled once per event. Should match the centvzbin histogram."       ,1  ,0 ,2                                                              ),GetNumberHist(khQAtocheckadressing             ,i,j));
    }
  }
  a->AddAt(new TH1D("overflow","overflow",3,0.0,1),GetNumberHist(kNofHistograms+1,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1));
  TH1::AddDirectory(statusAddDirectory);
  return 0;
}

int AliCorrelation3p::SetMultVZ(Double_t Mult, Double_t Vz)
{
  fMultiplicity=Mult;
  fVZ=Vz;
  fMBin = GetMultBin(fMultiplicity);
  fVzBin = GetZBin(fVZ);
  HistFill(kcentrvsvz,fMultiplicity,fVZ);
  HistFill(kcentrvsvzbin,fMBin,fVzBin);
  if (fMBin<0||fVzBin<0) return -1;
  HistFill(GetNumberHist(khQAtocheckadressing,fMBin,fVzBin),1.0);
  return 1;
}

bool AliCorrelation3p::CheckTrigger( AliVParticle* ptrigger, bool doHistogram)
{
  // check trigger particle cuts
  if (!ptrigger) return false;
  AliCFPI0 *phostrigger = dynamic_cast<AliCFPI0*>(ptrigger);
  if ((phostrigger!=NULL)&&(fTriggerType==tracks)) return false;//We only want tracks as triggers
  if ((phostrigger==NULL)&&(fTriggerType==pi0)) return false;//We only want pi0s as triggers
  if (ptrigger->Pt()<=fMinTriggerPt) return false;
  if (fMaxTriggerPt>fMinTriggerPt && ptrigger->Pt()>fMaxTriggerPt) return false;
  if (doHistogram) {
    HistFill(GetNumberHist(kHistpT,fMBin,fVzBin),ptrigger->Pt());
    HistFill(GetNumberHist(kHistPhi,fMBin,fVzBin),ptrigger->Phi());
    HistFill(GetNumberHist(kHistEta,fMBin,fVzBin),ptrigger->Eta());
    HistFill(GetNumberHist(kHistTriggerpT,fMBin,fVzBin),ptrigger->Pt());
    HistFill(GetNumberHist(kHistTriggerPhi,fMBin,fVzBin),ptrigger->Phi());
    HistFill(GetNumberHist(kHistTriggerEta,fMBin,fVzBin),ptrigger->Eta());
    HistFill(kHistpTTriggerallbins,ptrigger->Pt());
    HistFill(kHistPhiTriggerallbins,ptrigger->Phi());
    HistFill(kHistEtaTriggerallbins,ptrigger->Eta());    
  }
  return true;
}

bool AliCorrelation3p::CheckAssociated( AliVParticle* p, const AliVParticle* /*ptrigger*/, bool doHistogram)
{
  // check associated particle cuts
  if (!p) return false;
  AliCFPI0 *phosp = dynamic_cast<AliCFPI0*>(p);
  if (phosp != NULL) return false;//We only want tracks as associated particles
  if (p->Pt()<=fMinAssociatedPt) return false;
  if (fMaxAssociatedPt>fMinAssociatedPt && p->Pt()>fMaxAssociatedPt) return false;
  if (doHistogram) {
    HistFill(GetNumberHist(kHistpT,fMBin,fVzBin),p->Pt());
    HistFill(GetNumberHist(kHistPhi,fMBin,fVzBin),p->Phi());
    HistFill(GetNumberHist(kHistEta,fMBin,fVzBin),p->Eta());
    HistFill(GetNumberHist(kHistAssociatedpT,fMBin,fVzBin),p->Pt());
    HistFill(GetNumberHist(kHistAssociatedPhi,fMBin,fVzBin),p->Phi());
    HistFill(GetNumberHist(kHistAssociatedEta,fMBin,fVzBin),p->Eta());
    HistFill(kHistpTAssociatedallbins,p->Pt());
    HistFill(kHistPhiAssociatedallbins,p->Phi());
    HistFill(kHistEtaAssociatedallbins,p->Eta());    
  }
  return true;
}

int AliCorrelation3p::Fill(const AliVParticle* ptrigger, const AliVParticle* p1, const AliVParticle* p2, const int weight)
{
  /// fill histograms from particles, fills each histogram exactly once.
  if (!ptrigger || !p1 || !p2) {cout << "failed fill"<<endl;return -EINVAL;}
  const double Pii=TMath::Pi();
  HistFill(GetNumberHist(kHistNassoc,fMBin,fVzBin),weight);
  // phi difference associated 1 to trigger particle
  Double_t DeltaPhi1 = ptrigger->Phi() - p1->Phi();
  while(DeltaPhi1<-0.5*Pii||DeltaPhi1>1.5*Pii){
    if (DeltaPhi1<-0.5*Pii) DeltaPhi1 += 2*Pii;
    if (DeltaPhi1>1.5*Pii)  DeltaPhi1 -= 2*Pii;
  }
  // phi difference associated 2 to trigger particle
  Double_t DeltaPhi2 = ptrigger->Phi() - p2->Phi();
  while(DeltaPhi2<-0.5*Pii||DeltaPhi2>1.5*Pii){
    if (DeltaPhi2<-0.5*Pii) DeltaPhi2 += 2*Pii;
    if (DeltaPhi2>1.5*Pii)  DeltaPhi2 -= 2*Pii;
  }
  // eta difference
  Double_t DeltaEta12 = p1      ->Eta() - p2->Eta();
  HistFill(GetNumberHist(khPhiPhiDEta,fMBin,fVzBin),DeltaEta12,DeltaPhi1,DeltaPhi2);
  if(weight>1)  HistFill(GetNumberHist(khPhiPhiDEtaScaled,fMBin,fVzBin),DeltaEta12,DeltaPhi1,DeltaPhi2,1.0/(weight-1));
  return 0;
}

int AliCorrelation3p::Fill(const AliVParticle* ptrigger, const AliVParticle* p1)
{
  /// fill histograms from particles
  if (!ptrigger || !p1) return -EINVAL;
  const double Pii=TMath::Pi();
  // phi difference associated  to trigger particle
  Double_t DeltaPhi = ptrigger->Phi() - p1->Phi();
  while(DeltaPhi<-0.5*Pii||DeltaPhi>1.5*Pii){
    if (DeltaPhi<-0.5*Pii) DeltaPhi += 2*Pii;
    if (DeltaPhi>1.5*Pii)  DeltaPhi -= 2*Pii;
  }
  // eta difference
  Double_t DeltaEta  = ptrigger->Eta() - p1->Eta();
  HistFill(GetNumberHist(khPhiEta,fMBin,fVzBin),DeltaEta,DeltaPhi);//2p correlation
  return 0;
}

int AliCorrelation3p::FillTrigger()
{
  HistFill(GetNumberHist(kHistNTriggers,fMBin,fVzBin),0.5);//Increments number of triggers by one. Call before filling with any associated.
  return 1;
}

void AliCorrelation3p::Clear(Option_t * /*option*/)
{
  /// overloaded from TObject: cleanup
  return TObject::Clear();
}

void AliCorrelation3p::Print(Option_t *option) const
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
    cout << "unknown option '" << token << "'" << endl;
  }
  if (!bNoRuler)
    cout << "====================================================================" << endl;
  TObject::Print();
  cout << "p=" << this << endl;
  if (fHistograms) {
    fHistograms->Print();
  }
  if (fMixedEvent) {
    cout << "  ---- mixed event -------------------------------------------------" << endl;
    fMixedEvent->Print(Form("%s noruler", option));
  }
}

TObject* AliCorrelation3p::FindObject(const char *name) const
{
  /// overloaded from TObject: find object by name
  if (fHistograms) {
    TObject* o=fHistograms->FindObject(name);
    if (o) return o;
  }
  return NULL;
}

TObject* AliCorrelation3p::FindObject(const TObject *obj) const
{
  /// overloaded from TObject: find object by pointer
  if (fHistograms) {
    TObject* o=fHistograms->FindObject(obj);
    if (o) return o;
  }
  return NULL;
}

void AliCorrelation3p::SaveAs(const char *filename,Option_t */*option*/) const
{
  /// overloaded from TObject: save to file
  std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fHistograms) fHistograms->Write();
  output->Close();
}

AliCorrelation3p& AliCorrelation3p::operator+=(const AliCorrelation3p& other)
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
int AliCorrelation3p::Merge(TCollection* li)
{
   // interface for the TFileMerger
   TIter next(li);
   while (TObject* o = next()) {
      AliCorrelation3p* src=dynamic_cast<AliCorrelation3p*>(o);
      if (!src) {
	AliWarning(Form("skipping incompatible object %s of type %s", o->ClassName(), o->GetName()));
	continue;
      }
      (*this)+=(*src);
   }
   return 0;
}

const char* AliCorrelation3p::GetNameHist(const char* name, Int_t MBin, Int_t ZBin) const
{
  if((MBin>(fMBinEdges.GetSize()-1))||(ZBin>(fZBinEdges.GetSize()-1)))return name;//Index out of bounds
  Double_t Mlow = fMBinEdges.At(MBin);
  Double_t Mhigh = fMBinEdges.At(MBin+1);
  Double_t Zlow = fZBinEdges.At(ZBin);
  Double_t Zhigh = fZBinEdges.At(ZBin+1);
  return Form("%sM(%4.2f)->(%4.2f)Z(%4.2f)->(%4.2f)",name,Mlow,Mhigh,Zlow,Zhigh);
}

Int_t AliCorrelation3p::GetNumberHist(Int_t khist, Int_t Mbin, Int_t ZBin) const
{
  if((Mbin>(fMBinEdges.GetSize()-1))||(ZBin>(fZBinEdges.GetSize()-1)))return -1;//Index out of bounds
  int offset = kNonbinnedhists;//Offset so the extra histograms come first.
  return offset+khist+(Mbin+ZBin*(fMBinEdges.GetSize()))*kNofHistograms;//increase the number by number of hists for each M+zbin.
}

Int_t AliCorrelation3p::GetMultBin(Double_t Mult)
{
  int Nbins = fMBinEdges.GetSize()-1;
  for (int bin =0;bin<Nbins;bin++)if(Mult>=fMBinEdges[bin]&&Mult<fMBinEdges[bin+1])return bin;
  //if nothing happened in the loop, return -1:
  return -1;
}

Int_t AliCorrelation3p::GetZBin(Double_t Zvert)
{
  int Nbins = fZBinEdges.GetSize()-1;
  for (int bin =0;bin<Nbins;bin++)if(Zvert>=fZBinEdges[bin]&&Zvert<fZBinEdges[bin+1])return bin;
  //if nothing happened in the loop, return -1:
  return -1;
}

void AliCorrelation3p::HistFill(Int_t Histn,Double_t Val1)
{
  if(Histn>=0)dynamic_cast<TH1D*>(fHistograms->At(Histn))->Fill(Val1);
  else   dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms+1,fMBinEdges.GetSize()-1,fZBinEdges.GetSize()-1)))->Fill(0.5);//wrong histn
}
void AliCorrelation3p::HistFill(Int_t Histn,Double_t Val1, Double_t Val2)
{
  if(Histn>=0)dynamic_cast<TH2D*>(fHistograms->At(Histn))->Fill(Val1,Val2);
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}
void AliCorrelation3p::HistFill(Int_t Histn,Double_t Val1, Double_t Val2, Double_t Val3)
{
  if(Histn>=0)dynamic_cast<TH3D*>(fHistograms->At(Histn))->Fill(Val1,Val2, Val3);
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}
void AliCorrelation3p::HistFill(Int_t Histn,Double_t Val1, Double_t Val2, Double_t Val3, Double_t weight)
{
  if(Histn>=0)dynamic_cast<TH3D*>(fHistograms->At(Histn))->Fill(Val1,Val2,Val3, weight);
  else dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kNofHistograms,fMBinEdges.GetSize()-1,fZBinEdges.GetSize())+1))->Fill(0.5);//wrong histn
}

TH2D* AliCorrelation3p::slice(TH3D* hist, const char* option, Int_t firstbin, Int_t lastbin, const char* name, Bool_t baverage) const
{//option should be xy,zy,yx,zx,xz or yz.
  TString o = TString(option);
  TString namestring = TString(name);
  TH2D* Slice;
  if(o.CompareTo("xy")==0||o.CompareTo("yx")==0){
    if(o.CompareTo("xy")==0){Slice = (TH2D*)hist->Project3D("xy")->Clone(name);Slice->Reset("m");}
    if(o.CompareTo("yx")==0){Slice = (TH2D*)hist->Project3D("yx")->Clone(name);Slice->Reset("m");}
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
    if(o.CompareTo("xz")==0){Slice = (TH2D*)hist->Project3D("xz")->Clone(name);Slice->Reset("m");}
    if(o.CompareTo("zx")==0){Slice = (TH2D*)hist->Project3D("zx")->Clone(name);Slice->Reset("m");}
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
    if(o.CompareTo("yz")==0){Slice = (TH2D*)hist->Project3D("yz")->Clone(name);Slice->Reset("m");}
    if(o.CompareTo("zy")==0){Slice = (TH2D*)hist->Project3D("zy")->Clone(name);Slice->Reset("m");}

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

TH2D* AliCorrelation3p::DeltaEtaCut(TH3D* hist, const char* option, const char* name, Bool_t baverage) const
{
  //option can be: sameside    = if deltaPhi_1<pi/2, deltaPhi_2<pi/2 and  if deltaPhi_1>pi/2, deltaPhi_2>pi/2
  //		   lesspi2     = DeltaPhi_[12]<pi/2
  //		   lesspi4     = DeltaPhi_[12]<pi/4
  TString o = TString(option); 
  TAxis* dphi1axis=hist->GetYaxis();
  TAxis* dphi2axis=hist->GetZaxis();
  TH2D* Result = (TH2D*)hist->Project3D("yx")->Clone(name);
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
	bislesspi2 = bislesspi2&&(abs(dphi1axis->GetBinCenter(dPhi1)-dphi2axis->GetBinCenter(dPhi2)) <TMath::Pi()/2.0);
	bislesspi4 = o.CompareTo("lesspi4");
	bislesspi4 = bislesspi4&&(abs(dphi1axis->GetBinCenter(dPhi1)-dphi2axis->GetBinCenter(dPhi2)) <TMath::Pi()/4.0);
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

TCanvas* AliCorrelation3p::Makecanvas(TH1D* histtopl, TH1D* histtopm, TH1D* histtopr,TH1D* histmidl,TH1D* histmidm, TH1D* histmidr,TH1D* histbotl,TH1D* histbotm, TH1D* histbotr, const char* name, Bool_t Stats)
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

TCanvas* AliCorrelation3p::Makecanvas(TH2D* histtopl, TH2D* histtopr, TH2D* histbotl, TH2D* histbotr, const char* name, Bool_t Stats)
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

TCanvas* AliCorrelation3p::Makecanvas(TH2D* hist, const char* name, Bool_t Stats)
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

Double_t AliCorrelation3p::FindScalingfactor(const char* scalingmethod,TH2D* sighist,TH2D*mixhist)
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

Double_t AliCorrelation3p::GetPoint(TH1* hist, Double_t xpoint, Double_t ypoint, Double_t zpoint)
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


void AliCorrelation3p::AddHists(Bool_t isAverage, TH1* histtoadd, TH1* addedhist)
{
  TH1D* hist11d = dynamic_cast<TH1D*>(histtoadd);
  TH1D* hist21d = dynamic_cast<TH1D*>(addedhist);
  TH2D* hist12d = dynamic_cast<TH2D*>(histtoadd);
  TH2D* hist22d = dynamic_cast<TH2D*>(addedhist);  
  TH3D* hist13d = dynamic_cast<TH3D*>(histtoadd);
  TH3D* hist23d = dynamic_cast<TH3D*>(addedhist);   
  if(hist11d&&hist21d){
    Int_t nbinsx1 = hist11d->GetNbinsX();
    if(nbinsx1!=hist21d->GetNbinsX()){cout <<"The histograms do not match! TH1D* with different x dimensions."<<endl;
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
    if(nbinsx1!=hist22d->GetNbinsX()){cout <<"The histograms do not match! TH2D* with different x dimensions."<<endl;
    return ;}
    Int_t nbinsy1 = hist12d->GetNbinsY();
    if(nbinsy1!=hist22d->GetNbinsY()){cout <<"The histograms do not match! TH2D* with different y dimensions."<<endl;
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
    if(nbinsx1!=hist23d->GetNbinsX()){cout <<"The histograms do not match! TH3D* with different x dimensions."<<endl;
    return ;}
    Int_t nbinsy1 = hist13d->GetNbinsY();
    if(nbinsy1!=hist23d->GetNbinsY()){cout <<"The histograms do not match! TH3D* with different y dimensions."<<endl;
    return ;}
    Int_t nbinsz1 = hist13d->GetNbinsZ();
    if(nbinsz1!=hist23d->GetNbinsZ()){cout <<"The histograms do not match! TH2D* with different z dimensions."<<endl;
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
    cout <<"The histograms do not match!"<<endl;
    return ;    
  }
}


int AliCorrelation3p::MakeResultsFile(const char* scalingmethod)
{//This function makes a new file that contains all the histograms in all versions possible.
  TFile * outfile;
  TString dir = TString(scalingmethod);
  TDirectory * mixeddir=NULL;
  if(dir.CompareTo("")==0){ outfile = new TFile("results.root","RECREATE");outfile->cd();cout << "Recreate"<<endl;}
  else{
    outfile = new TFile("results.root","UPDATE");
    mixeddir = outfile->mkdir(dir.Data());
    dir.Append("/");
    mixeddir->cd();
  }
  TDirectory * dirmzbin = NULL;
  TDirectory * binstats = NULL;
  TDirectory * dirmzbinsig=NULL;
  TDirectory * dirmzbinmixed=NULL;
  TDirectory * dirmzbindiv=NULL;
  TCanvas    * tempcanvas=NULL;
  //Hists and directories for the total stuff:
  TDirectory * totbinstats = gDirectory->mkdir("bin_stats");
//   TDirectory * divided = gDirectory->mkdir("divided");
  TH1D * HistpT, * HistPhi,* HistEta, * HistTriggerpT, * HistTriggerPhi, * HistTriggerEta, * HistAssociatedpT, * HistAssociatedPhi,* HistAssociatedEta;
  TH3D * hPhiPhiDEtadiv=NULL;TH3D * hPhiPhiDEtadivscaled=NULL;
  TH2D * hDeltaPhidiv=NULL;TH2D* hDeltaPhidivscaled=NULL;TH2D * hDeltaPhineardiv=NULL;TH2D * hDeltaPhineardivscaled=NULL;TH2D * hDeltaPhimiddiv=NULL;TH2D * hDeltaPhimiddivscaled=NULL;TH2D * hDeltaPhifardiv=NULL;TH2D * hDeltaPhifardivscaled = NULL;
  TH2D * hPhiEta12div=NULL;TH2D* hPhiEta12_divscaled=NULL;TH2D * hPhiEta12_cut1div=NULL;TH2D * hPhiEta12_cut2div=NULL;TH2D * hPhiEta12_samesidediv=NULL;TH2D * hPhiEta12_sameside_divscaled=NULL;TH2D * hPhiEtadiv=NULL;
  Double_t navm, nav;
  Long_t NTriggers=0;
  Bool_t setAverage = kTRUE;  
  for(int mb =0;mb<fMBinEdges.GetSize()-1;mb++){
    for(int zb=0;zb<fZBinEdges.GetSize()-1;zb++){
      if(mixeddir) mixeddir->cd();
      else outfile->cd();
      dirmzbin = gDirectory->mkdir(GetNameHist("Bin",mb,zb));
      binstats = dirmzbin->mkdir("bin_stats");
      dirmzbinsig = dirmzbin->mkdir("same_event");
      dirmzbinmixed = dirmzbin->mkdir("mixed_event");
      dirmzbindiv = dirmzbin->mkdir("divided");
      {//for bin statistics
	binstats->cd();
	TH1D* savehist1 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistpT,mb,zb))->Clone("total_pT"));
	savehist1->SetTitle("pT of triggers and associated combined");
	savehist1->GetXaxis()->SetTitle("pT [GeV/c]");
	savehist1->GetYaxis()->SetTitle("# particles");
	savehist1->Write();
	if(mb == 0&&zb == 0)HistpT = savehist1;
	else HistpT->Add(savehist1);
	TH1D* savehist2 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistPhi,mb,zb))->Clone("total_Phi"));
	savehist2->SetTitle("Phi of triggers and associated combined");
	savehist2->GetXaxis()->SetTitle("#Phi [rad]");
	savehist2->GetYaxis()->SetTitle("# particles");
	savehist2->Write();
	if(mb == 0&&zb == 0)HistPhi = savehist2;
	else HistPhi->Add(savehist2);
	TH1D*  savehist3 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistEta,mb,zb))->Clone("total_Eta"));
	savehist3->SetTitle("Eta of triggers and associated combined");
	savehist3->GetXaxis()->SetTitle("#eta []");
	savehist3->GetYaxis()->SetTitle("# particles");
	savehist3->Write();
	if(mb == 0&&zb == 0)HistEta = savehist3;
	else HistEta->Add(savehist3);
	TH1D* savehist4 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistTriggerpT,mb,zb))->Clone("trigger_pT"));
	savehist4->SetTitle("pT of triggers");
	savehist4->GetXaxis()->SetTitle("pT [GeV/c]");
	savehist4->GetYaxis()->SetTitle("# triggers");
	savehist4->Write();
	if(mb == 0&&zb == 0)HistTriggerpT = savehist4;
	else HistTriggerpT->Add(savehist4);
	TH1D* savehist5 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistTriggerPhi,mb,zb))->Clone("trigger_Phi"));	
	savehist5->SetTitle("Phi of triggers");
	savehist5->GetXaxis()->SetTitle("#Phi [rad]");
	savehist5->GetYaxis()->SetTitle("# triggers");
	savehist5->Write();
	if(mb == 0&&zb == 0)HistTriggerPhi = savehist5;
	else HistTriggerPhi->Add(savehist5);
	TH1D* savehist6 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistTriggerEta,mb,zb))->Clone("trigger_Eta"));      
	savehist6->SetTitle("Eta of triggers");
	savehist6->GetXaxis()->SetTitle("#eta []");
	savehist6->GetYaxis()->SetTitle("# triggers");
	savehist6->Write();
	if(mb == 0&&zb == 0)HistTriggerEta = savehist6;
	else HistTriggerEta->Add(savehist6);
	TH1D* savehist7 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistAssociatedpT,mb,zb))->Clone("associated_pT"));
	savehist7->SetTitle("pT of associated");
	savehist7->GetXaxis()->SetTitle("pT [GeV/c]");
	savehist7->GetYaxis()->SetTitle("# associated");
	savehist7->Write();
	if(mb == 0&&zb == 0)HistAssociatedpT = savehist7;
	else HistAssociatedpT->Add(savehist7);
	TH1D* savehist8 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistAssociatedPhi,mb,zb))->Clone("associated_Phi"));
	savehist8->SetTitle("Phi of associated");
	savehist8->GetXaxis()->SetTitle("#Phi [rad]");
	savehist8->GetYaxis()->SetTitle("# associated");
	savehist8->Write();
	if(mb == 0&&zb == 0)HistAssociatedPhi = savehist8;
	else HistAssociatedPhi->Add(savehist8);
	TH1D* savehist9 = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistAssociatedEta,mb,zb))->Clone("associated_Eta"));
	savehist9->SetTitle("Eta of associated");
	savehist9->GetXaxis()->SetTitle("#eta []");
	savehist9->GetYaxis()->SetTitle("# associated");	
	savehist9->Write();
	if(mb == 0&&zb == 0)HistAssociatedEta = savehist9;
	else HistAssociatedEta->Add(savehist9);
	tempcanvas = Makecanvas(savehist1,savehist2,savehist3,savehist4,savehist5,savehist6,savehist7,savehist8,savehist9,"samestatscanvas",kFALSE);
	tempcanvas->Write();
	delete tempcanvas;
	tempcanvas=NULL;
	
	TH1D* savehist = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistNassoc,mb,zb))->Clone("number_of_associated"));
	savehist->SetTitle("Number of associated per trigger");
	savehist->GetXaxis()->SetTitle("# associated in a given fill");
	savehist->GetYaxis()->SetTitle("# fills");		
	savehist->Write();
	nav = 0;
	Double_t norm = 0;
	for(int i=1;i<savehist->GetNbinsX();i++){
	  nav += savehist->GetBinCenter(i)*savehist->GetBinContent(i);
	  norm += savehist->GetBinContent(i);
	}
	if(norm !=0)nav = nav/norm;
	else nav = 0;
	savehist = dynamic_cast<TH1D*>(fMixedEvent->fHistograms->At(GetNumberHist(kHistNassoc,mb,zb))->Clone("number_of_associated_mixed"));
	savehist->SetTitle("Number of associated per trigger in mixed");
	savehist->GetXaxis()->SetTitle("# associated in a given fill");
	savehist->GetYaxis()->SetTitle("# fills");		
	savehist->Write();
	navm = 0;
	Double_t normm = 0;
	for(int i=1;i<savehist->GetNbinsX();i++){
	  navm += savehist->GetBinCenter(i)*savehist->GetBinContent(i);
	  normm += savehist->GetBinContent(i);
	}
	if(normm !=0)navm = navm/normm;
	else nav = 0;
	savehist = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb))->Clone("number_of_triggers"));
	savehist->SetTitle("Total number of triggers filled");
	savehist->GetXaxis()->SetTitle("");
	savehist->GetYaxis()->SetTitle("# triggers");
	savehist->Write();
      }
/////////same event histograms
	dirmzbinsig->cd();
	TH1D* scalinghist= dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb))->Clone("number_of_triggers"));
	scalinghist->SetTitle("Total number of triggers filled");
	scalinghist->GetXaxis()->SetTitle("");
	scalinghist->GetYaxis()->SetTitle("# triggers");
	scalinghist->Write();
	TH3D* DPHIDPHIDETA = dynamic_cast<TH3D*>(fHistograms->At(GetNumberHist(khPhiPhiDEta,mb,zb))->Clone("DPhi_1_DPhi_2_DEta_12"));
	DPHIDPHIDETA->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}");
	DPHIDPHIDETA->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	DPHIDPHIDETA->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	DPHIDPHIDETA->GetZaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	DPHIDPHIDETA->SetTitleSize(0.045,"xyz");
	DPHIDPHIDETA->SetTitleOffset(1.2,"xy");
	DPHIDPHIDETA->SetTitleOffset(0.9,"z");
	DPHIDPHIDETA->Write();
	TH3D* DPHIDPHIDETAscaled = dynamic_cast<TH3D*>(fHistograms->At(GetNumberHist(khPhiPhiDEtaScaled,mb,zb))->Clone("DPhi_1_DPhi_2_DEta_12scaled"));
	DPHIDPHIDETAscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}");
	DPHIDPHIDETAscaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	DPHIDPHIDETAscaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	DPHIDPHIDETAscaled->GetZaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	DPHIDPHIDETAscaled->SetTitleSize(0.045,"xyz");
	DPHIDPHIDETAscaled->SetTitleOffset(1.2,"xy");
	DPHIDPHIDETAscaled->SetTitleOffset(0.9,"z");
	DPHIDPHIDETAscaled->Write();
/////////DPHIDPHI histograms:
	  TH2D* DPHIDPHI3 = slice(DPHIDPHIDETA,"yz",1,DPHIDPHIDETA->GetNbinsX(),"DPhi_1_DPHI_2",kFALSE);
	  DPHIDPHI3->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2}");
	  DPHIDPHI3->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3->GetZaxis()->SetTitle("# Pairs");
	  DPHIDPHI3->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3->Write();
	  TH2D* DPHIDPHI3scaled = slice(DPHIDPHIDETAscaled,"yz",1,DPHIDPHIDETA->GetNbinsX(),"DPhi_1_DPHI_2_scaled",kFALSE);
	  DPHIDPHI3scaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2}");
	  DPHIDPHI3scaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3scaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3scaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDPHI3scaled->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3scaled->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3scaled->Write();
	  TH2D* DPHIDPHI3near = slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-0.4),DPHIDPHIDETA->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near",kFALSE);
	  DPHIDPHI3near->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4");
	  DPHIDPHI3near->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3near->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3near->GetZaxis()->SetTitle("# Pairs");
	  DPHIDPHI3near->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3near->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3near->Write();
	  TH2D* DPHIDPHI3nearscaled = slice(DPHIDPHIDETAscaled,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-0.4),DPHIDPHIDETA->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near_scaled",kFALSE);
	  DPHIDPHI3nearscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4");
	  DPHIDPHI3nearscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3nearscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3nearscaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDPHI3nearscaled->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3nearscaled->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3nearscaled->Write();
	  TH2D* DPHIDPHI3mid = slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-1),DPHIDPHIDETA->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid",kFALSE);
	  DPHIDPHI3mid->Add(slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETA->GetXaxis()->FindBin(1),"temphist1",kFALSE));
	  DPHIDPHI3mid->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1");
	  DPHIDPHI3mid->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3mid->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3mid->GetZaxis()->SetTitle("# Pairs");
	  DPHIDPHI3mid->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3mid->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3mid->Write();
	  TH2D* DPHIDPHI3midscaled = slice(DPHIDPHIDETAscaled,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(-1),DPHIDPHIDETA->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid_scaled",kFALSE);
	  DPHIDPHI3midscaled->Add(slice(DPHIDPHIDETAscaled,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETA->GetXaxis()->FindBin(1),"temphist3",kFALSE));
	  DPHIDPHI3midscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1");
	  DPHIDPHI3midscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3midscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3midscaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDPHI3midscaled->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3midscaled->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3midscaled->Write();
	  TH2D* DPHIDPHI3far = slice(DPHIDPHIDETA,"yz",1,DPHIDPHIDETA->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far",kFALSE);
	  DPHIDPHI3far->Add(slice(DPHIDPHIDETA,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(1)+1,DPHIDPHIDETA->GetNbinsX(),"temphist2",kFALSE));
	  DPHIDPHI3far->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1");
	  DPHIDPHI3far->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3far->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3far->GetZaxis()->SetTitle("# Pairs");
	  DPHIDPHI3far->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3far->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3far->Write();      
	  TH2D* DPHIDPHI3farscaled = slice(DPHIDPHIDETAscaled,"yz",1,DPHIDPHIDETA->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far_scaled",kFALSE);
	  DPHIDPHI3farscaled->Add(slice(DPHIDPHIDETAscaled,"yz",DPHIDPHIDETA->GetXaxis()->FindBin(1)+1,DPHIDPHIDETA->GetNbinsX(),"temphist4",kFALSE));
	  DPHIDPHI3farscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1");
	  DPHIDPHI3farscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3farscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3farscaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDPHI3farscaled->SetTitleSize(0.045,"xyz");
	  DPHIDPHI3farscaled->SetTitleOffset(1.2,"xyz");
	  DPHIDPHI3farscaled->Write();      
	  tempcanvas= Makecanvas(DPHIDPHI3,DPHIDPHI3near,DPHIDPHI3mid,DPHIDPHI3far,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  tempcanvas= Makecanvas(DPHIDPHI3scaled,DPHIDPHI3nearscaled,DPHIDPHI3midscaled,DPHIDPHI3farscaled,"DPHIDPHIscaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	  
/////////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3 = slice(DPHIDPHIDETA,"yx",1,DPHIDPHIDETA->GetNbinsZ(),"DPhi_1_DEta_12",kFALSE);
	  DPHIDETA12_3->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12}");
	  DPHIDETA12_3->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12_3->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12_3->GetZaxis()->SetTitle("# Pairs");
	  DPHIDETA12_3->SetTitleSize(0.045,"xyz");
	  DPHIDETA12_3->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12_3->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12_3scaled = slice(DPHIDPHIDETAscaled,"yx",1,DPHIDPHIDETA->GetNbinsZ(),"DPhi_1_DEta_12_scaled",kFALSE);
	  DPHIDETA12_3scaled->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12}");
	  DPHIDETA12_3scaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12_3scaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12_3scaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDETA12_3scaled->SetTitleSize(0.045,"xyz");
	  DPHIDETA12_3scaled->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12_3scaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3scaled,"DPHIDEta12_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12DPHI12L2PI_3 = DeltaEtaCut(DPHIDPHIDETA,"lesspi2","DPhi_1_DEta_12_DPHI12_LESS_2PI",kFALSE);
	  DPHIDETA12DPHI12L2PI_3->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for #Delta#Phi_{12}<2#pi");
	  DPHIDETA12DPHI12L2PI_3->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12DPHI12L2PI_3->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12DPHI12L2PI_3->GetZaxis()->SetTitle("# Pairs");
	  DPHIDETA12DPHI12L2PI_3->SetTitleSize(0.045,"xyz");
	  DPHIDETA12DPHI12L2PI_3->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12DPHI12L2PI_3->Write();     
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L2PI_3,"DPHIDEta12_DPHI12less2Pi",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12DPHI12L4PI_3 = DeltaEtaCut(DPHIDPHIDETA,"less4pi","DPhi_1_DEta_12_DPHI12_LESS_4PI",kFALSE);
	  DPHIDETA12DPHI12L4PI_3->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for #Delta#Phi_{12}<4#pi");
	  DPHIDETA12DPHI12L4PI_3->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12DPHI12L4PI_3->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12DPHI12L4PI_3->GetZaxis()->SetTitle("# Pairs");
	  DPHIDETA12DPHI12L4PI_3->SetTitleSize(0.045,"xyz");
	  DPHIDETA12DPHI12L4PI_3->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12DPHI12L4PI_3->Write();
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L4PI_3,"DPHIDEta12_DPHI12less4Pi_3d",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3 =DeltaEtaCut(DPHIDPHIDETA,"sameside","DPhi_1_DEta_12_SameSide",kFALSE);
	  DPHIDETA12SameSide_3->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side");
	  DPHIDETA12SameSide_3->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12SameSide_3->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12SameSide_3->GetZaxis()->SetTitle("# Pairs");
	  DPHIDETA12SameSide_3->SetTitleSize(0.045,"xyz");
	  DPHIDETA12SameSide_3->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12SameSide_3->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3,"DPHIDEta12_SameSide_3d",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3scaled =DeltaEtaCut(DPHIDPHIDETAscaled,"sameside","DPhi_1_DEta_12_SameSide_scaled",kFALSE);
	  DPHIDETA12SameSide_3scaled->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side");
	  DPHIDETA12SameSide_3scaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12SameSide_3scaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12SameSide_3scaled->GetZaxis()->SetTitle("# Associated");
	  DPHIDETA12SameSide_3scaled->SetTitleSize(0.045,"xyz");
	  DPHIDETA12SameSide_3scaled->SetTitleOffset(1.2,"xyz");
	  DPHIDETA12SameSide_3scaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3scaled,"DPHIDEta12_SameSide_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA = dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(khPhiEta,mb,zb))->Clone("DPhi_DEta"));
	  DPHIDETA->SetTitle("#Delta#Phi vs #Delta#eta");
	  DPHIDETA->GetXaxis()->SetTitle("#Delta#eta []");
	  DPHIDETA->GetYaxis()->SetTitle("#Delta#Phi [rad]");
	  DPHIDETA->GetZaxis()->SetTitle("# Associated");
	  DPHIDETA->SetTitleSize(0.045,"xyz");
	  DPHIDETA->SetTitleOffset(1.2,"xyz");
	  DPHIDETA->Write();
	  tempcanvas= Makecanvas(DPHIDETA,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
      
	Double_t LeastContent = 1.0;
/////////mixed event histograms
	dirmzbinmixed->cd();
	TH1D* scalinghistm= dynamic_cast<TH1D*>(fMixedEvent->fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb))->Clone("number_of_triggers"));
	scalinghistm->SetTitle("Total number of times the mixed histogram was filled with a trigger");
	scalinghistm->GetXaxis()->SetTitle("");
	scalinghistm->GetYaxis()->SetTitle("# triggers");
	scalinghistm->Write();
	TH3D* DPHIDPHIDETAm = dynamic_cast<TH3D*>(fMixedEvent->fHistograms->At(GetNumberHist(khPhiPhiDEta,mb,zb))->Clone("DPhi_1_DPhi_2_DEta_12"));
	DPHIDPHIDETAm->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}");
	DPHIDPHIDETAm->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	DPHIDPHIDETAm->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	DPHIDPHIDETAm->GetZaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	DPHIDPHIDETAm->SetTitleSize(0.045,"xyz");
	DPHIDPHIDETAm->SetTitleOffset(1.2,"xy");
	DPHIDPHIDETAm->SetTitleOffset(0.9,"z");
	DPHIDPHIDETAm->Write();
	TH3D* DPHIDPHIDETAmscaled = dynamic_cast<TH3D*>(fMixedEvent->fHistograms->At(GetNumberHist(khPhiPhiDEtaScaled,mb,zb))->Clone("DPhi_1_DPhi_2_DEta_12scaled"));
	DPHIDPHIDETAmscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} vs #Delta#eta_{12}");
	DPHIDPHIDETAmscaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	DPHIDPHIDETAmscaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	DPHIDPHIDETAmscaled->GetZaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	DPHIDPHIDETAmscaled->SetTitleSize(0.045,"xyz");
	DPHIDPHIDETAmscaled->SetTitleOffset(1.2,"xy");
	DPHIDPHIDETAmscaled->SetTitleOffset(0.9,"z");
	DPHIDPHIDETAmscaled->Write();
//////////DPHIDPHI histograms:
	  TH2D* DPHIDPHI3m = slice(DPHIDPHIDETAm,"yz",1,DPHIDPHIDETAm->GetNbinsX(),"DPhi_1_DPHI_2",kFALSE);
	  DPHIDPHI3m->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2}");
	  DPHIDPHI3m->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3m->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3m->SetTitleSize(0.045,"xy");
	  DPHIDPHI3m->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3m,0.0,0.0)>LeastContent){DPHIDPHI3m->Scale(1.0/GetPoint(DPHIDPHI3m,0.0,0.0));}
	  else DPHIDPHI3m->Scale(0.0);
	  DPHIDPHI3m->Write();
	  TH2D* DPHIDPHI3mscaled = slice(DPHIDPHIDETAmscaled,"yz",1,DPHIDPHIDETAmscaled->GetNbinsX(),"DPhi_1_DPHI_2_scaled",kFALSE);
	  DPHIDPHI3mscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2}");
	  DPHIDPHI3mscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3mscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3mscaled->SetTitleSize(0.045,"xy");
	  DPHIDPHI3mscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
  	  if(GetPoint(DPHIDPHI3mscaled,0.0,0.0)>LeastContent){DPHIDPHI3mscaled->Scale(1.0/GetPoint(DPHIDPHI3mscaled,0.0,0.0));}
	  else DPHIDPHI3mscaled->Scale(0.0);
	  DPHIDPHI3mscaled->Write();
	  TH2D* DPHIDPHI3nearm = slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4),DPHIDPHIDETAm->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near",kFALSE);
	  DPHIDPHI3nearm->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4");
	  DPHIDPHI3nearm->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3nearm->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3nearm->SetTitleSize(0.045,"xy");
	  DPHIDPHI3nearm->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3nearm,0.0,0.0)>LeastContent){DPHIDPHI3nearm->Scale(1.0/GetPoint(DPHIDPHI3nearm,0.0,0.0));}
	  else DPHIDPHI3nearm->Scale(0.0);	  
	  DPHIDPHI3nearm->Write();
	  TH2D* DPHIDPHI3nearmscaled = slice(DPHIDPHIDETAmscaled,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4),DPHIDPHIDETAm->GetXaxis()->FindBin(0.4),"DPhi_1_DPHI_2_near_scaled",kFALSE);
	  DPHIDPHI3nearmscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}<=0.4");
	  DPHIDPHI3nearmscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3nearmscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3nearmscaled->SetTitleSize(0.045,"xy");
	  DPHIDPHI3nearmscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3nearmscaled,0.0,0.0)>LeastContent){DPHIDPHI3nearmscaled->Scale(1.0/GetPoint(DPHIDPHI3nearmscaled,0.0,0.0));}
	  else DPHIDPHI3nearmscaled->Scale(0.0);
	  DPHIDPHI3nearmscaled->Write();
	  TH2D* DPHIDPHI3midm = slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-1),DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid",kFALSE);
	  DPHIDPHI3midm->Add(slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETAm->GetXaxis()->FindBin(1),"temphist1",kFALSE));
	  DPHIDPHI3midm->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1");
	  DPHIDPHI3midm->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3midm->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3midm->SetTitleSize(0.045,"xy");
	  DPHIDPHI3midm->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3midm,0.0,0.0)>LeastContent){DPHIDPHI3midm->Scale(1.0/GetPoint(DPHIDPHI3midm,0.0,0.0));}
	  else DPHIDPHI3midm->Scale(0.0);
	  DPHIDPHI3midm->Write();
	  TH2D* DPHIDPHI3midmscaled = slice(DPHIDPHIDETAmscaled,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(-1),DPHIDPHIDETAm->GetXaxis()->FindBin(-0.4)-1,"DPhi_1_DPHI_2_mid_scaled",kFALSE);
	  DPHIDPHI3midmscaled->Add(slice(DPHIDPHIDETAmscaled,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(0.4)+1,DPHIDPHIDETAm->GetXaxis()->FindBin(1),"temphist3",kFALSE));
	  DPHIDPHI3midmscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for 0.4<#Delta#eta_{12}<=1");
	  DPHIDPHI3midmscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3midmscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3midmscaled->SetTitleSize(0.045,"xy");
	  DPHIDPHI3midmscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3midmscaled,0.0,0.0)>LeastContent){DPHIDPHI3midmscaled->Scale(1.0/GetPoint(DPHIDPHI3midmscaled,0.0,0.0));}
	  else DPHIDPHI3midmscaled->Scale(0.0);	  
	  DPHIDPHI3midmscaled->Write();
	  TH2D* DPHIDPHI3farm = slice(DPHIDPHIDETAm,"yz",1,DPHIDPHIDETAm->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far",kFALSE);
	  DPHIDPHI3farm->Add(slice(DPHIDPHIDETAm,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(1)+1,DPHIDPHIDETAm->GetNbinsX(),"temphist2",kFALSE));
	  DPHIDPHI3farm->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1");
	  DPHIDPHI3farm->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3farm->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3farm->SetTitleSize(0.045,"xy");
	  DPHIDPHI3farm->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3farm,0.0,0.0)>LeastContent){DPHIDPHI3farm->Scale(1.0/GetPoint(DPHIDPHI3farm,0.0,0.0));}
	  else DPHIDPHI3farm->Scale(0.0);	  
	  DPHIDPHI3farm->Write();      
	  TH2D* DPHIDPHI3farmscaled = slice(DPHIDPHIDETAmscaled,"yz",1,DPHIDPHIDETAm->GetXaxis()->FindBin(-1)-1,"DPhi_1_DPHI_2_far_scaled",kFALSE);
	  DPHIDPHI3farmscaled->Add(slice(DPHIDPHIDETAmscaled,"yz",DPHIDPHIDETAm->GetXaxis()->FindBin(1)+1,DPHIDPHIDETAm->GetNbinsX(),"temphist4",kFALSE));
	  DPHIDPHI3farmscaled->SetTitle("#Delta#Phi_{1} vs #Delta#Phi_{2} for #Delta#eta_{12}>1");
	  DPHIDPHI3farmscaled->GetXaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDPHI3farmscaled->GetYaxis()->SetTitle("#Delta#Phi_{2} [rad]");
	  DPHIDPHI3farmscaled->SetTitleSize(0.045,"xy");
	  DPHIDPHI3farmscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDPHI3farmscaled,0.0,0.0)>LeastContent){DPHIDPHI3farmscaled->Scale(1.0/GetPoint(DPHIDPHI3farmscaled,0.0,0.0));}
	  else DPHIDPHI3farmscaled->Scale(0.0);	
	  DPHIDPHI3farmscaled->Write();      
	  tempcanvas= Makecanvas(DPHIDPHI3m,DPHIDPHI3nearm,DPHIDPHI3midm,DPHIDPHI3farm,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  tempcanvas= Makecanvas(DPHIDPHI3mscaled,DPHIDPHI3nearmscaled,DPHIDPHI3midmscaled,DPHIDPHI3farmscaled,"DPHIDPHIscaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	  
/////////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3m = slice(DPHIDPHIDETAm,"yx",1,DPHIDPHIDETAm->GetNbinsZ(),"DPhi_1_DEta_12",kFALSE);
// 	  DPHIDETA12_3m->Sumw2();
	  DPHIDETA12_3m->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12}");
	  DPHIDETA12_3m->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12_3m->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12_3m->SetTitleSize(0.045,"xy");
	  DPHIDETA12_3m->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12_3m,0.0,0.0)>LeastContent){DPHIDETA12_3m->Scale(1.0/GetPoint(DPHIDETA12_3m,0.0,0.0));}
	  else DPHIDETA12_3m->Scale(0.0);	  
	  DPHIDETA12_3m->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3m,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12_3mscaled = slice(DPHIDPHIDETAmscaled,"yx",1,DPHIDPHIDETAm->GetNbinsZ(),"DPhi_1_DEta_12_scaled",kFALSE);
	  DPHIDETA12_3mscaled->Sumw2();
	  DPHIDETA12_3mscaled->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12}");
	  DPHIDETA12_3mscaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12_3mscaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12_3mscaled->SetTitleSize(0.045,"xy");
	  DPHIDETA12_3mscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12_3mscaled,0.0,0.0)>LeastContent){DPHIDETA12_3mscaled->Scale(1.0/GetPoint(DPHIDETA12_3mscaled,0.0,0.0));}
	  else DPHIDETA12_3mscaled->Scale(0.0);	  
	  DPHIDETA12_3mscaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3mscaled,"DPHIDEta12_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12DPHI12L2PI_3m = DeltaEtaCut(DPHIDPHIDETAm,"lesspi2","DPhi_1_DEta_12_DPHI12_LESS_2PI",kFALSE);
	  DPHIDETA12DPHI12L2PI_3m->Sumw2();
	  DPHIDETA12DPHI12L2PI_3m->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for #Delta#Phi_{12}<2#pi");
	  DPHIDETA12DPHI12L2PI_3m->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12DPHI12L2PI_3m->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12DPHI12L2PI_3m->SetTitleSize(0.045,"xy");
	  DPHIDETA12DPHI12L2PI_3m->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12DPHI12L2PI_3m,0.0,0.0)>LeastContent){DPHIDETA12DPHI12L2PI_3m->Scale(1.0/GetPoint(DPHIDETA12DPHI12L2PI_3m,0.0,0.0));}
	  else DPHIDETA12DPHI12L2PI_3m->Scale(0.0);
	  DPHIDETA12DPHI12L2PI_3m->Write();      
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L2PI_3m,"DPHIDEta12_DPHI12less2Pi",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12DPHI12L4PI_3m = DeltaEtaCut(DPHIDPHIDETAm,"less4pi","DPhi_1_DEta_12_DPHI12_LESS_4PI",kFALSE);
	  DPHIDETA12DPHI12L4PI_3m->Sumw2();
	  DPHIDETA12DPHI12L4PI_3m->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for #Delta#Phi_{12}<4#pi");
	  DPHIDETA12DPHI12L4PI_3m->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12DPHI12L4PI_3m->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12DPHI12L4PI_3m->SetTitleSize(0.045,"xy");
	  DPHIDETA12DPHI12L4PI_3m->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12DPHI12L4PI_3m,0.0,0.0)>LeastContent){DPHIDETA12DPHI12L4PI_3m->Scale(1.0/GetPoint(DPHIDETA12DPHI12L4PI_3m,0.0,0.0));}
	  else DPHIDETA12DPHI12L4PI_3m->Scale(0.0);
	  DPHIDETA12DPHI12L4PI_3m->Write();
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L4PI_3m,"DPHIDEta12_DPHI12less4Pi",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3m =DeltaEtaCut(DPHIDPHIDETAm,"sameside","DPhi_1_DEta_12_SameSide",kFALSE);
	  DPHIDETA12SameSide_3m->Sumw2();
	  DPHIDETA12SameSide_3m->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side");
	  DPHIDETA12SameSide_3m->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12SameSide_3m->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12SameSide_3m->SetTitleSize(0.045,"xy");
	  DPHIDETA12SameSide_3m->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12SameSide_3m,0.0,0.0)>LeastContent){DPHIDETA12SameSide_3m->Scale(1.0/GetPoint(DPHIDETA12SameSide_3m,0.0,0.0));}
	  else DPHIDETA12SameSide_3m->Scale(0.0);	  
	  DPHIDETA12SameSide_3m->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3m,"DPHIDEta12_SameSide",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3mscaled =DeltaEtaCut(DPHIDPHIDETAmscaled,"sameside","DPhi_1_DEta_12_SameSide_scaled",kFALSE);
	  DPHIDETA12SameSide_3mscaled->Sumw2();
	  DPHIDETA12SameSide_3mscaled->SetTitle("#Delta#Phi_{1} vs #Delta#eta_{12} for both associated on the same side");
	  DPHIDETA12SameSide_3mscaled->GetXaxis()->SetTitle("#Delta#eta_{12} []");
	  DPHIDETA12SameSide_3mscaled->GetYaxis()->SetTitle("#Delta#Phi_{1} [rad]");
	  DPHIDETA12SameSide_3mscaled->SetTitleSize(0.045,"xy");
	  DPHIDETA12SameSide_3mscaled->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETA12SameSide_3mscaled,0.0,0.0)>LeastContent){DPHIDETA12SameSide_3mscaled->Scale(1.0/GetPoint(DPHIDETA12SameSide_3mscaled,0.0,0.0));}
	  else DPHIDETA12SameSide_3mscaled->Scale(0.0);	  
	  DPHIDETA12SameSide_3mscaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3mscaled,"DPHIDEta12_SameSide_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	  
	  TH2D* DPHIDETAm = dynamic_cast<TH2D*>(fMixedEvent->fHistograms->At(GetNumberHist(khPhiEta,mb,zb))->Clone("DPhi_DEta"));
	  DPHIDETAm->Sumw2();
	  DPHIDETAm->SetTitle("#Delta#Phi vs #Delta#eta");
	  DPHIDETAm->GetXaxis()->SetTitle("#Delta#eta []");
	  DPHIDETAm->GetYaxis()->SetTitle("#Delta#Phi [rad]");
	  DPHIDETAm->SetTitleSize(0.045,"xy");
	  DPHIDETAm->SetTitleOffset(1.2,"xy");
	  //Scaling
	  if(GetPoint(DPHIDETAm,0.0,0.0)>10){DPHIDETAm->Scale(1.0/GetPoint(DPHIDETAm,0.0,0.0));}
	  else DPHIDETAm->Scale(0.0);	 	  
	  DPHIDETAm->Write();
	  tempcanvas= Makecanvas(DPHIDETAm,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	

      Double_t resultscalingfactor = 1.0;//Scale the result with 1/ntriggers
      if(scalinghist->Integral()!=0)resultscalingfactor=1.0/scalinghist->Integral();
      else resultscalingfactor=0.0;
      
      if(GetPoint(DPHIDPHI3m,0.0,0.0)!=0){//signal divided by mixed if there was enough in mixed
	NTriggers+=scalinghist->Integral();
	dirmzbindiv->cd();
	//the same hists as in same and mixed.
	TH1D* scalinghistdiv = (TH1D*)scalinghist->Clone("number_of_triggers");
	scalinghistdiv->Divide(scalinghistm);
	scalinghistdiv->Write();
 	TH3D* DPHIDPHIDETAdiv = (TH3D*)DPHIDPHIDETA->Clone("DPhi_1_DPhi_2_DEta_12");
	DPHIDPHIDETAdiv->Divide(DPHIDPHIDETAm);
	if(setAverage)DPHIDPHIDETAdiv->Scale(resultscalingfactor);
	DPHIDPHIDETAdiv->Write();
	DPHIDPHIDETAdiv->SetBit(TH3D::kIsAverage,setAverage);
	if(!hPhiPhiDEtadiv) hPhiPhiDEtadiv=DPHIDPHIDETAdiv;
	else AddHists(setAverage,hPhiPhiDEtadiv,DPHIDPHIDETAdiv);
	DPHIDPHIDETAdiv->SetBit(TH3D::kIsAverage,kFALSE);
	TH3D* DPHIDPHIDETAdivscaled = (TH3D*)DPHIDPHIDETAscaled->Clone("DPhi_1_DPhi_2_DEta_12scaled");
	DPHIDPHIDETAdivscaled->Divide(DPHIDPHIDETAm);
	if(setAverage)DPHIDPHIDETAdivscaled->Scale(resultscalingfactor);
	DPHIDPHIDETAdivscaled->Write();
	DPHIDPHIDETAdivscaled->SetBit(TH3D::kIsAverage,setAverage);
	if(!hPhiPhiDEtadivscaled) hPhiPhiDEtadivscaled=DPHIDPHIDETAdivscaled;
	else AddHists(setAverage,hPhiPhiDEtadivscaled,DPHIDPHIDETAdivscaled);
	DPHIDPHIDETAdivscaled->SetBit(TH3D::kIsAverage,kFALSE);
	//DPHIDPHI histograms:
	  TH2D* DPHIDPHI3div = (TH2D*)DPHIDPHI3->Clone("DPhi_1_DPHI_2");
	  DPHIDPHI3div->Divide(DPHIDPHI3m);
	  if(setAverage) DPHIDPHI3div->Scale(resultscalingfactor);
	  DPHIDPHI3div->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhidiv) hDeltaPhidiv=DPHIDPHI3div;
	  else AddHists(setAverage,hDeltaPhidiv,DPHIDPHI3div);
	  DPHIDPHI3div->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3div->Scale(resultscalingfactor);
	  DPHIDPHI3div->Write();
	  TH2D* DPHIDPHI3divscaled = (TH2D*)DPHIDPHI3scaled->Clone("DPhi_1_DPHI_2_scaled");
	  DPHIDPHI3divscaled->Divide(DPHIDPHI3mscaled);
	  if(setAverage)DPHIDPHI3divscaled->Scale(resultscalingfactor);
	  DPHIDPHI3divscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhidivscaled) hDeltaPhidivscaled=DPHIDPHI3divscaled;
	  else AddHists(setAverage,hDeltaPhidivscaled,DPHIDPHI3divscaled);
	  DPHIDPHI3divscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3divscaled->Scale(resultscalingfactor);
	  DPHIDPHI3divscaled->Write();
	  TH2D* DPHIDPHI3neardiv = (TH2D*)DPHIDPHI3near->Clone("DPhi_1_DPHI_2_near");
	  DPHIDPHI3neardiv->Divide(DPHIDPHI3nearm);
	  if(setAverage)DPHIDPHI3neardiv->Scale(resultscalingfactor);
	  DPHIDPHI3neardiv->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhineardiv) hDeltaPhineardiv=DPHIDPHI3neardiv;
	  else AddHists(setAverage,hDeltaPhineardiv,DPHIDPHI3neardiv);
	  DPHIDPHI3neardiv->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3neardiv->Scale(resultscalingfactor);
	  DPHIDPHI3neardiv->Write();
	  TH2D* DPHIDPHI3neardivscaled = (TH2D*)DPHIDPHI3nearscaled->Clone("DPhi_1_DPHI_2_near_scaled");
	  DPHIDPHI3neardivscaled->Divide(DPHIDPHI3nearmscaled);
	  if(setAverage)DPHIDPHI3neardivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3neardivscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhineardivscaled) hDeltaPhineardivscaled=DPHIDPHI3neardivscaled;
	  else AddHists(setAverage,hDeltaPhineardivscaled,DPHIDPHI3neardivscaled);
	  DPHIDPHI3neardivscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3neardivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3neardivscaled->Write();
	  TH2D* DPHIDPHI3middiv = (TH2D*)DPHIDPHI3mid->Clone("DPhi_1_DPHI_2_mid");
	  DPHIDPHI3middiv->Divide(DPHIDPHI3midm);
	  if(setAverage)DPHIDPHI3middiv->Scale(resultscalingfactor);
	  DPHIDPHI3middiv->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhimiddiv) hDeltaPhimiddiv=DPHIDPHI3middiv;
	  else AddHists(setAverage,hDeltaPhimiddiv,DPHIDPHI3middiv);
	  DPHIDPHI3middiv->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3middiv->Scale(resultscalingfactor);
	  DPHIDPHI3middiv->Write();
	  TH2D* DPHIDPHI3middivscaled = (TH2D*)DPHIDPHI3mid->Clone("DPhi_1_DPHI_2_mid_scaled");
	  DPHIDPHI3middivscaled->Divide(DPHIDPHI3midmscaled);
	  if(setAverage)DPHIDPHI3middivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3middivscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhimiddivscaled) hDeltaPhimiddivscaled=DPHIDPHI3middivscaled;
	  else AddHists(setAverage,hDeltaPhimiddivscaled,DPHIDPHI3middivscaled);
	  DPHIDPHI3middivscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3middivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3middivscaled->Write();
	  TH2D* DPHIDPHI3fardiv = (TH2D*)DPHIDPHI3far->Clone("DPhi_1_DPHI_2_far");
	  DPHIDPHI3fardiv->Divide(DPHIDPHI3farm);
	  if(setAverage)DPHIDPHI3fardiv->Scale(resultscalingfactor);
	  DPHIDPHI3fardiv->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhifardiv) hDeltaPhifardiv=DPHIDPHI3fardiv;
	  else AddHists(setAverage,hDeltaPhifardiv,DPHIDPHI3fardiv);
	  DPHIDPHI3fardiv->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3fardiv->Scale(resultscalingfactor);
	  DPHIDPHI3fardiv->Write();      
	  TH2D* DPHIDPHI3fardivscaled = (TH2D*)DPHIDPHI3farscaled->Clone("DPhi_1_DPHI_2_far_scaled");
	  DPHIDPHI3fardivscaled->Divide(DPHIDPHI3farmscaled);
	  if(setAverage)DPHIDPHI3fardivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3fardivscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hDeltaPhifardivscaled) hDeltaPhifardivscaled=DPHIDPHI3fardivscaled;
	  else AddHists(setAverage,hDeltaPhifardivscaled,DPHIDPHI3fardivscaled);
	  DPHIDPHI3fardivscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDPHI3fardivscaled->Scale(resultscalingfactor);
	  DPHIDPHI3fardivscaled->Write();
	  tempcanvas= Makecanvas(DPHIDPHI3div,DPHIDPHI3neardiv,DPHIDPHI3middiv,DPHIDPHI3fardiv,"DPHIDPHI",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  tempcanvas= Makecanvas(DPHIDPHI3divscaled,DPHIDPHI3neardivscaled,DPHIDPHI3middivscaled,DPHIDPHI3fardivscaled,"DPHIDPHIscaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
//////////DPHI DETA HISTOGRAMS:
	  TH2D* DPHIDETA12_3div = (TH2D*) DPHIDETA12_3->Clone("DPhi_1_DEta_12");
	  DPHIDETA12_3div->Divide(DPHIDETA12_3m);
	  if(setAverage)DPHIDETA12_3div->Scale(resultscalingfactor);
	  if(!hPhiEta12div) hPhiEta12div=DPHIDETA12_3div;
	  else AddHists(setAverage,hPhiEta12div,DPHIDETA12_3div);
	  if(!setAverage)DPHIDETA12_3div->Scale(resultscalingfactor);
	  DPHIDETA12_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3div,"DPHIDEta12",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	  
	  TH2D* DPHIDETA12_3divscaled = (TH2D*) DPHIDETA12_3scaled->Clone("DPhi_1_DEta_12_scaled");
	  DPHIDETA12_3divscaled->Divide(DPHIDETA12_3mscaled);
	  if(setAverage)DPHIDETA12_3divscaled->Scale(resultscalingfactor);
	  DPHIDETA12_3divscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEta12_divscaled) hPhiEta12_divscaled=DPHIDETA12_3divscaled;
	  else AddHists(setAverage,hPhiEta12_divscaled,DPHIDETA12_3divscaled);
	  DPHIDETA12_3divscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETA12_3divscaled->Scale(resultscalingfactor);
	  DPHIDETA12_3divscaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12_3divscaled,"DPHIDEta12_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;	
	  TH2D* DPHIDETA12DPHI12L2PI_3div = (TH2D*)DPHIDETA12DPHI12L2PI_3->Clone("DPhi_1_DEta_12_DPHI12_LESS_2PI");
	  DPHIDETA12DPHI12L2PI_3div->Divide(DPHIDETA12DPHI12L2PI_3m);
	  if(setAverage)DPHIDETA12DPHI12L2PI_3div->Scale(resultscalingfactor);
	  DPHIDETA12DPHI12L2PI_3div->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEta12_cut1div) hPhiEta12_cut1div=DPHIDETA12DPHI12L2PI_3div;
	  else AddHists(setAverage,hPhiEta12_cut1div,DPHIDETA12DPHI12L2PI_3div);
	  DPHIDETA12DPHI12L2PI_3div->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETA12DPHI12L2PI_3div->Scale(resultscalingfactor);
	  DPHIDETA12DPHI12L2PI_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L2PI_3div,"DPHIDEta12_DPHI12less2Pi",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12DPHI12L4PI_3div = (TH2D*)DPHIDETA12DPHI12L4PI_3->Clone("DPhi_1_DEta_12_DPHI12_LESS_4PI");
	  DPHIDETA12DPHI12L4PI_3div->Divide(DPHIDETA12DPHI12L4PI_3m);
	  if(setAverage)DPHIDETA12DPHI12L4PI_3div->Scale(resultscalingfactor);
	  DPHIDETA12DPHI12L4PI_3div->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEta12_cut2div) hPhiEta12_cut2div=DPHIDETA12DPHI12L4PI_3div;
	  else AddHists(setAverage,hPhiEta12_cut2div,DPHIDETA12DPHI12L4PI_3div);
	  DPHIDETA12DPHI12L4PI_3div->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETA12DPHI12L4PI_3div->Scale(resultscalingfactor);
	  DPHIDETA12DPHI12L4PI_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12DPHI12L4PI_3div,"DPHIDEta12_DPHI12less4Pi",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3div = (TH2D*)DPHIDETA12SameSide_3->Clone("DPhi_1_DEta_12_SameSide");
	  DPHIDETA12SameSide_3div->Divide(DPHIDETA12SameSide_3m);
	  if(setAverage)DPHIDETA12SameSide_3div->Scale(resultscalingfactor);
	  DPHIDETA12SameSide_3div->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEta12_samesidediv) hPhiEta12_samesidediv=DPHIDETA12SameSide_3div;
	  else AddHists(setAverage,hPhiEta12_samesidediv,DPHIDETA12SameSide_3div);
	  DPHIDETA12SameSide_3div->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETA12SameSide_3div->Scale(resultscalingfactor);
	  DPHIDETA12SameSide_3div->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3div,"DPHIDEta12_SameSide",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETA12SameSide_3divscaled = (TH2D*)DPHIDETA12SameSide_3scaled->Clone("DPhi_1_DEta_12_SameSide_scaled");
	  DPHIDETA12SameSide_3divscaled->Divide(DPHIDETA12SameSide_3mscaled);
	  if(setAverage)DPHIDETA12SameSide_3divscaled->Scale(resultscalingfactor);
	  DPHIDETA12SameSide_3divscaled->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEta12_sameside_divscaled) hPhiEta12_sameside_divscaled=DPHIDETA12SameSide_3divscaled;
	  else AddHists(setAverage,hPhiEta12_sameside_divscaled,DPHIDETA12SameSide_3divscaled);
	  DPHIDETA12SameSide_3divscaled->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETA12SameSide_3divscaled->Scale(resultscalingfactor);
	  DPHIDETA12SameSide_3divscaled->Write();
	  tempcanvas= Makecanvas(DPHIDETA12SameSide_3divscaled,"DPHIDEta12_SameSide_scaled",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
	  TH2D* DPHIDETAdiv = (TH2D*)DPHIDETA->Clone("DPhi_DEta");
	  DPHIDETAdiv->Divide(DPHIDETAm);
	  if(setAverage)DPHIDETAdiv->Scale(resultscalingfactor);
	  DPHIDETAdiv->SetBit(TH2D::kIsAverage,setAverage);
	  if(!hPhiEtadiv) hPhiEtadiv=DPHIDETAdiv;
	  else AddHists(setAverage,hPhiEtadiv,DPHIDETAdiv);
	  DPHIDETAdiv->SetBit(TH2D::kIsAverage,kFALSE);
	  if(!setAverage)DPHIDETAdiv->Scale(resultscalingfactor);
	  DPHIDETAdiv->Write();
	  tempcanvas= Makecanvas(DPHIDETAdiv,"DPHIDEta",kFALSE);
	  tempcanvas->Write();
	  delete tempcanvas;
      }
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
  }
  outfile->Close();
  delete outfile;
  return 1;
}
