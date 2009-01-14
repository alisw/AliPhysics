#include "iostream"
#include "AliFMDInput.h"
#include "AliFMDHit.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliTrackReference.h"
#include "AliFMDStripIndex.h"

//Script to compare the hits and the FMD track references for one event.
//To run:
//>gSystem->Load("libFMDutil")
//>.L CompareTrackRefsAndHits.C++
//>ReadHits t
//>t.Run()
//Note that the order of hits and trackrefs is different.

class ReadHits : public AliFMDInput{

private:
  Int_t nHits;
  
public:
  
  
  
  ReadHits(){
    AddLoad(kKinematics);
    AddLoad(kHits);
    nHits = 0;
  }

  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p) 
  {
    nHits++;
    std::cout<<hit->Px()<<"   "<<hit->Py()<<"   "<<hit->Pz()<<std::endl;
    std::cout<<hit->Detector()<<"   "<<hit->Ring()<<"    "<<hit->Sector()<<"    "<<hit->Strip()<<std::endl;
    return kTRUE;
  }
  Bool_t Finish()
  {
    Int_t nTracks = 0;
    TFile* f=TFile::Open("TrackRefs.root");
    TTree* tree = (TTree*)f->Get("Event0/TreeTR");
    
    
    TClonesArray* array=new TClonesArray("AliTrackReference");
    
    tree->SetBranchAddress("TrackReferences",&array);
    
    UShort_t det,sec,strip;
    Char_t ring;
    for (int i=0; i<tree->GetEntries(); i++) {
      tree->GetEvent(i);
    for(Int_t j = 0; j <array->GetEntriesFast();j++) {
      
      AliTrackReference* track = static_cast<AliTrackReference*>(array->At(j));
      
      if(track->DetectorId()==AliTrackReference::kFMD) {
	nTracks++;
	AliFMDStripIndex::Unpack(track->UserId(),det,ring,sec,strip);
	std::cout<<track->Px()<<"   "<<track->Py()<<"   "<<track->Pz()<<"   "<<track->UserId()<<endl;
	std::cout<<det<<"   "<<ring<<"    "<<sec<<"    "<<strip<<std::endl;
      }
      
    }
    
    }
    std::cout<<nTracks<<"     "<<nHits<<std::endl;
    return kTRUE;
  }
  
};
