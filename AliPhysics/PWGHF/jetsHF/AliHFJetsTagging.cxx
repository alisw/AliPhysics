/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* Manager class for HF jet analysis   */

/* Mailto: andrea.rossi@cern.ch,       *
 *         elena.bruna@cern.ch,        *
 *         min.jung.kweon@cern.ch,     *
 *         linus.feldkamp@cern.ch,     *
 *         svallero@to.infn.it         *
 *         s.lapointe@cern.ch          */

#include "AliHFJetsTagging.h"


#include <TNamed.h>

#include <TSystem.h>
#include <TApplication.h>
#include <TVector2.h>
#include <TMath.h>
#include "TClonesArray.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"


#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"


#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"

#include <vector>
#include <utility>
#include <algorithm>

ClassImp(AliHFJetsTagging)

AliHFJetsTagging::AliHFJetsTagging():
    TNamed(),
    fSelectedTracks(0x0){
      // method non implemented
    }

AliHFJetsTagging::AliHFJetsTagging(const char* name): 
  TNamed(name,name),
  fSelectedTracks(0x0){

    // method non implemented

}

AliHFJetsTagging::~AliHFJetsTagging(){
  //delete fSelectedTracks;
}

AliHFJetsTagging &AliHFJetsTagging::operator=(const AliHFJetsTagging &c)
{
  // assigment operator
  if (this != &c)
    ((AliHFJetsTagging &) c).Copy(*this);
  return *this;
}

Bool_t AliHFJetsTagging::IsEventSelectedTrackCounting(const AliAODEvent*event){
  //Additional event selection on vertexing quality
  if(!event) return kFALSE;
  if(event->IsPileupFromSPD(5,0.8, 3.0, 2.0, 5.0)) return kFALSE;
  if(fabs(event->GetPrimaryVertex()->GetZ())>10.)  return kFALSE;
  const AliVVertex *trkVtx =    dynamic_cast<const AliVVertex*>(event->GetPrimaryVertex()) ;
  const AliVVertex* spdVtx =     dynamic_cast<const AliVVertex*>(event->GetPrimaryVertexSPD()) ;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks"))  return kFALSE;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ((vtxTyp.Contains("vertexer: Z") && (zRes>0.25)))return kFALSE;
  if ((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5))return kFALSE;  
  if(trkVtx->GetNContributors()<2) return kFALSE;
  if(spdVtx->GetNContributors()<1) return kFALSE;
  return kTRUE;
}

AliAODMCParticle*  AliHFJetsTagging::IsMCJetParton(const TClonesArray *arrayMC,const AliEmcalJet *jet, Double_t radius){
  Int_t mcEntries=arrayMC->GetEntriesFast();
  Double_t ptpart=-1;
  Double_t dR=-99;
  const Int_t arraySize=99;

  Int_t countpart[arraySize],countpartcode[arraySize],maxInd=-1,count=0;
  Double_t maxPt=0;
  for(Int_t ii=0;ii<arraySize;ii++){
    countpart[ii]=0;
    countpartcode[ii]=0;
  }
  
  for(Int_t ii=0;ii<mcEntries;ii++){
    AliAODMCParticle* part =  (AliAODMCParticle*)  arrayMC->At(ii);
    if(!part)continue;
    Int_t pdgcode=part->GetPdgCode();
    if(abs(pdgcode)==21 || ( abs(pdgcode)>=1 && abs(pdgcode)<=5)){
      ptpart=part->Pt();
      dR = jet->DeltaR(part);

      if(dR<radius){
        if(abs(pdgcode)==5){
	  return part;
	}
        else{
	  if (count >arraySize-1) return 0x0; 
          countpartcode[count]=pdgcode;
          countpart[count]=ii;
          if(ptpart>maxPt){
            maxPt=ptpart;
            maxInd=ii;
          }
          count++;
        }
      }
    }
  }
  for(Int_t ij=0;ij<count;ij++){
      if(abs(countpartcode[ij])==4) 
	return (AliAODMCParticle*)arrayMC->At(countpart[ij]);
    }
  if(maxInd>-1){
    AliAODMCParticle* partmax = (AliAODMCParticle*)arrayMC->At(maxInd);
    return partmax;
  }
  return 0x0;
}

AliMCParticle*  AliHFJetsTagging::IsMCJetParton(const AliMCEvent *mcevent,const AliEmcalJet *jet, Double_t radius){
	//Printf("%s:%i",__FUNCTION__,__LINE__);

  Int_t mcEntries=mcevent->GetNumberOfTracks();
	//Printf("%s:%i",__FUNCTION__,__LINE__);

  Double_t ptpart=-1;
  Double_t dR=-99;
  const Int_t arraySize=99;
	//Printf("%s:%i",__FUNCTION__,__LINE__);

  Int_t countpart[arraySize],countpartcode[arraySize],maxInd=-1,count=0;
  Double_t maxPt=0;
  for(Int_t ii=0;ii<arraySize;ii++){
    countpart[ii]=0;
    countpartcode[ii]=0;
  }
	//Printf("%s:%i",__FUNCTION__,__LINE__);

  for(Int_t ii=0;ii<mcEntries;ii++){
		//Printf("%s:%i",__FUNCTION__,__LINE__);

    AliMCParticle* part =  (AliMCParticle*)  mcevent->GetTrack(ii);
    if(!part)continue;
	//Printf("%s:%i",__FUNCTION__,__LINE__);

    Int_t pdgcode=part->PdgCode();
	//Printf("%s:%i",__FUNCTION__,__LINE__);

    if(abs(pdgcode)==21 || ( abs(pdgcode)>=1 && abs(pdgcode)<=5)){
		//Printf("%s:%i",__FUNCTION__,__LINE__);

    	ptpart=part->Pt();
		//Printf("%s:%i",__FUNCTION__,__LINE__);

      dR = jet->DeltaR(part);
		//Printf("%s:%i",__FUNCTION__,__LINE__);

      if(dR<radius){
        if(abs(pdgcode)==5){
	  return part;
	}
        else{
	  if (count >arraySize-1) return 0x0;
          countpartcode[count]=pdgcode;
          countpart[count]=ii;
          if(ptpart>maxPt){
            maxPt=ptpart;
            maxInd=ii;
          }
          count++;
        }
      }
    }
  }
  for(Int_t ij=0;ij<count;ij++){
      if(abs(countpartcode[ij])==4)
	return (AliMCParticle*)mcevent->GetTrack(countpart[ij]);
    }
  if(maxInd>-1){
    AliMCParticle* partmax = (AliMCParticle*)mcevent->GetTrack(maxInd);
    return partmax;
  }
  return 0x0;
}


AliAODMCParticle* AliHFJetsTagging::IsMCJetMeson(const TClonesArray *arrayMC, const AliEmcalJet *jet, Double_t radius){

  Int_t MCentries=arrayMC->GetEntries();
  Double_t ptpart=-1;

  AliAODMCParticle* parton3 = 0;
  AliAODMCParticle* charm_parton=0; 
 
  Double_t dR=-99;
  Int_t maxInd=-1;
  Double_t maxPt=0;
  for(Int_t ii=0;ii<MCentries;ii++){
    AliAODMCParticle* part =  (AliAODMCParticle*)  arrayMC->At(ii);
    if(!part)continue;
    Int_t pdgcode=part->GetPdgCode();
       if(IsDMeson(pdgcode) || IsBMeson(pdgcode)){
	 
      ptpart=part->Pt();
      dR = jet->DeltaR(part);
      
      if(dR<radius){
        if(IsBMeson(pdgcode)){  
          return part;
        }
        else {
	  if(IsDMeson(pdgcode)) {
	    charm_parton=part;
	  }
	 if(ptpart>maxPt){
            maxPt=ptpart;
            maxInd=ii;
	 }
        }
     }
  }
 }   
    if(charm_parton) return charm_parton;

  if(maxInd>-1){
    AliAODMCParticle* partmax = (AliAODMCParticle*)arrayMC->At(maxInd);
    return partmax;
  }
  return parton3;
}



AliMCParticle* AliHFJetsTagging::IsMCJetMeson(const AliMCEvent *mcevent, const AliEmcalJet *jet, Double_t radius){

  Int_t MCentries=mcevent->GetNumberOfTracks();
  Double_t ptpart=-1;

  AliMCParticle* parton3 = 0;
  AliMCParticle* charm_parton=0;

  Double_t dR=-99;
  Int_t maxInd=-1;
  Double_t maxPt=0;
  for(Int_t ii=0;ii<MCentries;ii++){
    AliMCParticle* part =  (AliMCParticle*)  mcevent->GetTrack(ii);
    if(!part)continue;
    Int_t pdgcode=part->PdgCode();
       if(IsDMeson(pdgcode) || IsBMeson(pdgcode)){

      ptpart=part->Pt();
      dR = jet->DeltaR(part);

      if(dR<radius){
        if(IsBMeson(pdgcode)){
          return part;
        }
        else {
	  if(IsDMeson(pdgcode)) {
	    charm_parton=part;
	  }
	 if(ptpart>maxPt){
            maxPt=ptpart;
            maxInd=ii;
	 }
        }
     }
  }
 }
    if(charm_parton) return charm_parton;

  if(maxInd>-1){
    AliMCParticle* partmax =  (AliMCParticle*)  mcevent->GetTrack(maxInd);
    return partmax;
  }
  return parton3;
}


Bool_t AliHFJetsTagging::IsBMeson(Int_t pc){
	int bPdG[] = {511,521,10511,10521,513,523,10513,10523,20513,20523,20513,20523,515,525,531,
	10531,533,10533,20533,535,541,10541,543,10543,20543,545,551,10551,100551,
	110551,200551,210551,553,10553,20553,30553,100553,110553,120553,130553,200553,210553,220553,
	300553,9000533,9010553,555,10555,20555,100555,110555,120555,200555,557,100557, 5122, 5112, 5212, 5222,5114,5214,5224,5132,5232,5312,5322,5314,5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,
5534,5542,5544,5554};
	for(int i=0;i< (int)(sizeof(bPdG)/sizeof(int));++i) if(abs(pc) == bPdG[i]) return true;
        return false;
}
Bool_t AliHFJetsTagging::IsDMeson(Int_t pc){
	int bPdG[] = {411,421,10411,10421,413,423,10413,10423,20431,20423,415,
	425,431,10431,433,10433,20433,435,441,10441,100441,443,10443,20443,
	100443,30443,9000443,9010443,9020443,445,100445, 4122,4222,4212,4112,4224,4214,4114,4232,4132,4322,4312,4324,4314,4332,4334,4412,4422,4414,4424,4432,4434,4444
        };
	for(int i=0;i< (int)(sizeof(bPdG)/sizeof(int));++i) if(abs(pc) == bPdG[i]) return true;
        return false;
}

Double_t AliHFJetsTagging::RelativePhi(Double_t mphi,Double_t vphi) {
  Double_t dPhi = mphi-vphi;
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return dPhi;
}

AliAODMCParticle* AliHFJetsTagging::GetAODMCParticleFromAodTrack(const TClonesArray *arrayMC,const AliAODTrack* track){
  AliAODMCParticle* gentrack=0x0;
  if(!arrayMC || !track) return 0x0;
  Int_t label = TMath::Abs(track->GetLabel());
  if (label > arrayMC->GetEntriesFast()) return 0x0;
  gentrack = dynamic_cast<AliAODMCParticle*>(arrayMC->At(label));
  if(!gentrack) return 0x0;
  return gentrack;
}

Bool_t AliHFJetsTagging::GetSignedRPhiImpactParameter(const AliVEvent *event, const AliVTrack *track, const AliEmcalJet *jet,
                                                      Double_t &signedimpactparameter, Double_t d[2], Double_t cov[3]) {
  Int_t skipped[2];
  Float_t diamondcovxy[3];
  Double_t posAtDCA[2] = {-999,-999};
  Double_t covar[3]= {-999,-999,-999};
  Double_t pV[3]		= {0.,0.,0.};
  Double_t pTrack[3] 	= {0.,0.,0.};
  const Double_t kBeampiperadius=2.6;
  Double_t absIP = 0;
  Double_t scalarP = 0;

  if(!event|| !track ||!jet) return kFALSE;
  AliVVertex * vertex	= (AliVVertex *)event->GetPrimaryVertex();
  if(!vertex) return kFALSE;
  if(vertex->GetNContributors()< 100) {
    //Recalculate  primary vertex without "track", if number of contributers to the prim. vertex is less than 100
    
    AliVertexerTracks 	vertexer(event->GetMagneticField());
    vertexer.SetITSMode();
    vertexer.SetMinClusters(4);
    skipped[0] 		= (Int_t)track->GetID();
    vertexer.SetSkipTracks(1,skipped);	
    vertexer.SetConstraintOn();
    event->GetDiamondCovXY(diamondcovxy);
    Double_t bpos[3]=	{event->GetDiamondX(), event->GetDiamondY(),0.};
    Double_t bcov[6]=	{diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
    AliESDVertex *diamond = new AliESDVertex(bpos,bcov,1.,1);
    vertexer.SetVtxStart((AliESDVertex*)diamond);
    delete diamond;
    diamond=0x0;
    vertex = (AliVVertex*)vertexer.FindPrimaryVertex(event);
  }
  //Propagate to primary vertex
   AliExternalTrackParam betp;
   betp.CopyFromVTrack((AliVTrack*)track);
    if(!(betp.PropagateToDCA(vertex,event->GetMagneticField(),kBeampiperadius, posAtDCA, covar)))
      return kFALSE;
    vertex->GetXYZ(pV);
    betp.GetXYZ(pTrack);
    
    Double_t imparvector[3]={pTrack[0] - pV[0],pTrack[1] - pV[1],pTrack[2] - pV[2]};
    absIP = sqrt(imparvector[0]*imparvector[0]+ imparvector[1]*imparvector[1] + imparvector[2]*imparvector[2]);
    
    scalarP = (imparvector[0]*jet->Px()+imparvector[1]*jet->Py()+imparvector[2]*jet->Pz())/(absIP*jet->P());
    signedimpactparameter = scalarP>=0 ? TMath::Abs(posAtDCA[0]) : -1*TMath::Abs(posAtDCA[0]);
    d[0] = posAtDCA[0];
    d[1] = posAtDCA[1];
    cov[0] = covar[0];
    cov[1] = covar[1];
    cov[2] = covar[2];

    return kTRUE;
}


std::vector< std::pair<Double_t,Int_t> > AliHFJetsTagging::GetSortedListOfSignedRPhiImpactParameters(const AliVEvent    *event,
                                                                                                     const AliEmcalJet  *jet,
                                                                                                     const TClonesArray *tracksAccepted){

  typedef std::pair<Double_t, Int_t> doubleidx_pair;
  std::vector<doubleidx_pair> pair_list;
  
  if ( !(event) || !(jet) || !(tracksAccepted) || (tracksAccepted->GetEntriesFast() < 1) )
    return pair_list;

  for(Int_t i_entry = 0; i_entry < jet->GetNumberOfTracks(); i_entry++) {
    
    AliVTrack *track = dynamic_cast<AliVTrack *>( tracksAccepted->At(i_entry) );
    if ( ! track ) {
      continue;
    }
    
    Double_t dv[2] = {0., 0.};
    Double_t covv[3]={0., 0., 0.};
    Double_t signedImpactParameter = 0.;
    if (GetSignedRPhiImpactParameter(event, track, jet, signedImpactParameter, dv, covv) )
      pair_list.push_back(std::make_pair(signedImpactParameter, i_entry));
  }
  std::stable_sort(pair_list.begin() , pair_list.end() , sort_descend());

  return pair_list;
}
