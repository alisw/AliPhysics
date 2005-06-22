// $Id$

/* Kt Jet finder class ala hep-ex/0005012 
   for CDF/D0 Run II */

#include <Riostream.h>
#include <vector>
#include <set>

#include <TClonesArray.h>
#include <TIterator.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TObjArray.h>

#include "AliJFPreCluster.h"
#include "AliJFCluster.h"
#include "AliJFClusterDifference.h"
#include "AliJFKtJet.h"
#include "AliJFKtJetFinder.h"

ClassImp(AliJFKtJetFinder)

AliJFKtJetFinder::AliJFKtJetFinder(Int_t n) : AliJFJetFinder(n)
{
}

AliJFKtJetFinder::~AliJFKtJetFinder()
{
  Clean();
}

Bool_t AliJFKtJetFinder::IsAcceptedParticle(TParticle *p)
{
  if(p->GetStatusCode()%100!=1) return kFALSE;

  Int_t pcode=p->GetPdgCode();  

  if((!fEM) && ((pcode==11)||(pcode==-11)||(pcode==22))) return kFALSE;

  TParticlePDG *pdg=p->GetPDG();
  Float_t ch=pdg->Charge(); 
  if((!fCharged)&&(ch)) return kFALSE;
  if((!fNeutral)&&(!ch)) return kFALSE;

  Float_t eta=p->Eta();
  if((eta<fEtaMin)||(eta>fEtaMax)) return kFALSE;

  Float_t phi=p->Phi();
  if((phi<fPhiMin)||(phi>fPhiMax)) return kFALSE;

  Float_t pt=p->Pt();
  if((pt<fPtMin)||(pt>fPtMax)) return kFALSE;

  return kTRUE;
}

Int_t AliJFKtJetFinder::Init(TClonesArray *particles)
{ //create precluster according to particles 
  //each precluster corresponds to a particle
  //in the accepted range
  if(particles==NULL) return -1;

  TIterator *iter=particles->MakeIterator();
  TParticle *p;
  Int_t ret=0;

  while((p=(TParticle*)iter->Next()) != NULL){
    if(IsAcceptedParticle(p)){
      ret++;

      //fPreClusterList.push_back(new AliJFPreCluster(p));
      fPreClusterList.push_back(new AliJFPreCluster(p->Px(),p->Py(),p->Pz(),-1,p));
    }
  } 

#if 0
  for(vector<AliJFPreCluster*>::iterator c=fPreClusterList.begin();c!=fPreClusterList.end();c++){
    cout << *(*c) << endl;
  }
  //exit(1);
#endif

  for(vector<AliJFPreCluster*>::iterator i=fPreClusterList.begin();i!=fPreClusterList.end();i++){
    fClusterList.push_back(new AliJFCluster(*(*i)));
  }

#if 0
  for(vector<AliJFCluster*>::iterator c=fClusterList.begin();c!=fClusterList.end();c++){
    cout << *(*c) << endl;
  }
  //exit(1);
#endif

  AliJFClusterDifference diff;
  for(vector<AliJFCluster*>::iterator i=fClusterList.begin();i!=fClusterList.end();i++){
    for(vector<AliJFCluster*>::iterator j=i;j!=fClusterList.end();j++){
      diff.SetValues(*i,*j);
      fClusterDiffSet.insert(diff);
    }
  }

#if 0
  for(multiset<AliJFClusterDifference>::iterator pos=fClusterDiffSet.begin();pos!=fClusterDiffSet.end();pos++){
    cout << *pos << endl;
  }
  //exit(1);
#endif

  return ret;
}

Int_t AliJFKtJetFinder::Run()
{
  //loop over stored diffence objects until set is empty
  while(!fClusterDiffSet.empty()){
    //get first element of sorted set
    multiset<AliJFClusterDifference>::iterator pos=fClusterDiffSet.begin();

    //cout << pos->GetDij() << endl;
    if((!pos->IsValidPointer())||(!pos->IsValidEntry())){ //delete old value left in set
      fClusterDiffSet.erase(pos);

      continue;
    }

    //get information on stored difference
    AliJFCluster *i=pos->GetI();
    AliJFCluster *j=pos->GetJ();
    bool isjetfound=pos->IsDiagonal();
    fClusterDiffSet.erase(pos); //take it out

    if(isjetfound){ //found jet
      i->MarkIsJet();

      if(fNJets==fNJetsMax){
	fNJetsMax+=fNJets;
	fJets.Expand(fNJetsMax);
      }

      fJets.AddAt(new AliJFKtJet(),fNJets);
      fJet=(AliJFKtJet*)fJets.At(fNJets);
      fJet->SetNJet(++fNJets);
      fJet->AddJet(i);
      fJet->Update();

    } else { //combine cluster and make new difference objects
      i->CombineCluster(*j);

      AliJFClusterDifference diff; //create and insert new difference objects
      for(vector<AliJFCluster*>::iterator c=fClusterList.begin();c!=fClusterList.end();c++){
	if((*c)->IsValid()) {
	  diff.SetValues(i,*c);
	  fClusterDiffSet.insert(diff);
	}
      }
    } //end if combine cluster    
  } //end while

  fJets.Sort();
  //fJets.Print();

  return fNJets;
}

void AliJFKtJetFinder::Debug()
{
  for(vector<AliJFPreCluster*>::iterator c=fPreClusterList.begin();c!=fPreClusterList.end();c++){
    cout << *(*c) << endl;
  }

  for(vector<AliJFCluster*>::iterator c=fClusterList.begin();c!=fClusterList.end();c++){
    //(*c)->Print();
  }

}

void AliJFKtJetFinder::Clean()
{
  if(!fPreClusterList.empty()){
    for(vector<AliJFPreCluster*>::iterator i=fPreClusterList.begin();i!=fPreClusterList.end();i++){
      delete (*i);
    }
    fPreClusterList.erase(fPreClusterList.begin(),fPreClusterList.end());
  }
  if(!fClusterList.empty()){
    for(vector<AliJFCluster*>::iterator i=fClusterList.begin();i!=fClusterList.end();i++){
      delete (*i);
    }
    fClusterList.erase(fClusterList.begin(),fClusterList.end());
  }

  if(!fClusterDiffSet.empty()) fClusterDiffSet.clear();

  AliJFJetFinder::Clean();
}
