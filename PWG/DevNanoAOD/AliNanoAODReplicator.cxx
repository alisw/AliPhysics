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

// $Id: AliNanoAODReplicator.cxx 56492 2012-05-15 18:42:47Z pcrochet $

// Implementation of a branch replicator 
// to produce nanoAODs.
//
//
// This replicator is in charge of replicating the tracks,vertices,headers
// branches of a standard AOD or ESD file into a nanoAODs 
// (AliAOD.Special.root)
// 
// The class was inspired by AliAODMuonReplicator
// 
// Author: Michele Floris, michele.floris@cern.ch


class AliESDv0;
class AliESDVertex;
class AliAODVertex;
class AliAODRecoDecay;

#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTZERO.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
#include "AliAnalysisCuts.h"
#include "TF1.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include <iostream>
#include <cassert>
#include "AliESDtrack.h"
#include "TObjArray.h"
#include "AliAnalysisFilter.h"
#include "AliNanoAODTrack.h"

#include <TFile.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TList.h>
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliVertexerTracks.h"
#include "AliKFVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAnalysisFilter.h"

#include "AliNanoAODReplicator.h"
#include "TH1.h"
#include "TCanvas.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODCustomSetter.h"

using std::cout;
using std::endl;

ClassImp(AliNanoAODReplicator)

//_____________________________________________________________________________
AliNanoAODReplicator::AliNanoAODReplicator() :
AliAODBranchReplicator(), 
  fTrackCut(0), fTracks(0x0), fHeader(0x0), fNTracksVariables(0), // FIXME: Start using cuts, and check if fNTracksVariables is needed
  fVertices(0x0), 
  fList(0x0),
  fMCParticles(0x0),
  fMCHeader(0x0),
  fMCMode(0),
  fLabelMap(),
  fParticleSelected(),
  fVarList(""),
  fVarListHeader(""),
  fCustomSetter(0){
  // Default ctor. we need it to avoid instantiating a wrong mapping when reading from file 
  }

AliNanoAODReplicator::AliNanoAODReplicator(const char* name, const char* title,
					   const char * varlist,
					   AliAnalysisCuts* trackCut,
					   Int_t mcMode
					     ) :
  AliAODBranchReplicator(name,title), 

  fTrackCut(trackCut), fTracks(0x0), fHeader(0x0), fNTracksVariables(0), // FIXME: Start using cuts, and check if fNTracksVariables is needed
  fVertices(0x0), 
  fList(0x0),
  fMCParticles(0x0),
  fMCHeader(0x0),
  fMCMode(mcMode),
  fLabelMap(),
  fParticleSelected(),
  fVarList(varlist),
  fVarListHeader(""),// FIXME: this should be set to a meaningful value: add an arg to the constructor
  fCustomSetter(0)
{
  // default ctor
  AliNanoAODTrackMapping * tm =new AliNanoAODTrackMapping(fVarList);
  fNTracksVariables = tm->GetSize();  
  //  tm->Print();
}

//_____________________________________________________________________________
AliNanoAODReplicator::~AliNanoAODReplicator()
{
  // dtor
  delete fTrackCut;
  delete fList;
}

//_____________________________________________________________________________
void AliNanoAODReplicator::SelectParticle(Int_t i)
{
  // taking the absolute values here, need to take care 
  // of negative daughter and mother
  // IDs when setting!
  
  if (!IsParticleSelected(TMath::Abs(i)))
  {
    fParticleSelected.Add(TMath::Abs(i),1);    
  }
}

//_____________________________________________________________________________
Bool_t AliNanoAODReplicator::IsParticleSelected(Int_t i)  
{
  // taking the absolute values here, need to take 
  // care with negative daughter and mother
  // IDs when setting!
  return (fParticleSelected.GetValue(TMath::Abs(i))==1);
}


//_____________________________________________________________________________
void AliNanoAODReplicator::CreateLabelMap(const AliAODEvent& source)
{  
  //
  // this should be called once all selections are done 
  // This method associates to the original index of the mc particle
  // (i) the new one (j). J runs only over particles which are
  // actually kept.
  //
  
  fLabelMap.Delete();
  
  TClonesArray* mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  Int_t j(0);
  Int_t i(0); // We need i, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
  
  TIter next(mcParticles);

  while ( next() )
  {
    if (IsParticleSelected(i))
    {
      fLabelMap.Add(i,j++);
      //      std::cout << i <<  "->" << j-1 << std::endl;
    }
    ++i;
  }  


}

//_____________________________________________________________________________
Int_t AliNanoAODReplicator::GetNewLabel(Int_t i) 
{
  // Gets the label from the new created Map
  // Call CreatLabelMap before
  // otherwise only 0 returned
  return fLabelMap.GetValue(TMath::Abs(i));
}

//_____________________________________________________________________________
void AliNanoAODReplicator::FilterMC(const AliAODEvent& source)
{
  // Filter MC information


  AliAODMCHeader* mcHeader(0x0);
  TClonesArray* mcParticles(0x0);
  
  fParticleSelected.Delete();

  //  std::cout << "MC Mode: " << fMCMode << ", Tracks " << fTracks->GetEntries() << std::endl;
  
  if ( fMCMode>=2 && !fTracks->GetEntries() ) {
    return;
  }
  // for fMCMode==1 we only copy MC information for events where there's at least one muon track
    
  mcHeader = static_cast<AliAODMCHeader*>(source.FindListObject(AliAODMCHeader::StdBranchName()));
  
  if ( mcHeader ) 
    {
      *fMCHeader = *mcHeader;
    }
  
  
  mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  if ( mcParticles && fMCMode>=2 )
    {
      // keep all primaries
      TIter nextPart(mcParticles);
      static Int_t iev = -1; // FIXME: remove this (debug)
      iev ++;
      AliAODMCParticle * prim = 0;
      Int_t iprim = 0;  // We need iprim, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
      // also select all charged primaries 
      while ((prim = (AliAODMCParticle*) nextPart())) {
	if(prim->IsPhysicalPrimary() && prim->Charge()) SelectParticle(iprim);
	// FIXME DEBUG
	if(iev == 2009) {
	  // std::cout << "IEV " << iev << std::endl;
	  // std::cout << " PART " << iprim << " " << prim->IsPhysicalPrimary() <<","<<prim->Charge() << "=" << IsParticleSelected(iprim) <<  std::endl;
	  if(iprim == 15) {
	    prim->Print();
	  }
	  
	}
	iprim++;
      } 

      // loop on (kept) tracks to find their ancestors
      TIter nextTRACK(fTracks);
      AliNanoAODTrack* track;
    
      while ( ( track = static_cast<AliNanoAODTrack*>(nextTRACK()) ) )
	{
	  Int_t label = TMath::Abs(track->GetLabel()); 
      
	  while ( label >= 0 ) 
	    {
	      SelectParticle(label);
	      AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(mcParticles->UncheckedAt(label));
	      if (!mother)
		{
		  AliError(Form("Got a null mother ! Check that ! (label %d",label)); // FIXME: I think this error is not needed
		  label = -1;
		}
	      else
		{
		  label = mother->GetMother();// do not only keep particles which created a track, but all their mothers
		}
	    }
	}
    
      CreateLabelMap(source);
    
      // Actual filtering and label remapping (shamelessly taken for the implementation of AliAODHandler::StoreMCParticles)
      TIter nextMC(mcParticles);
      AliAODMCParticle* p;
      Int_t nmc(0);  // We need nmc, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
      Int_t nmcout(0);
    
      while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
	{
	  AliAODMCParticle c(*p);

	  if ( IsParticleSelected(nmc) )
	    {
	      // 
	      Int_t d0 =  p->GetDaughter(0);
	      Int_t d1 =  p->GetDaughter(1);
	      Int_t m =   p->GetMother();
        
	      // other than for the track labels, negative values mean
	      // no daughter/mother so preserve it
        
	      if(d0<0 && d1<0)
		{
		  // no first daughter -> no second daughter
		  // nothing to be done
		  // second condition not needed just for sanity check at the end
		  c.SetDaughter(0,d0);
		  c.SetDaughter(1,d1);
		} else if(d1 < 0 && d0 >= 0) 
		{
		  // Only one daughter
		  // second condition not needed just for sanity check at the end
		  if(IsParticleSelected(d0))
		    {
		      c.SetDaughter(0,GetNewLabel(d0));
		    } else 
		    {
		      c.SetDaughter(0,-1);
		    }
		  c.SetDaughter(1,d1);
		}
	      else if (d0 > 0 && d1 > 0 )
		{
		  // we have two or more daughters loop on the stack to see if they are
		  // selected
		  Int_t d0tmp = -1;
		  Int_t d1tmp = -1;
		  for (int id = d0; id<=d1;++id)
		    {
		      if (IsParticleSelected(id))
			{
			  if(d0tmp==-1)
			    {
			      // first time
			      d0tmp = GetNewLabel(id);
			      d1tmp = d0tmp; // this is to have the same schema as on the stack i.e. with one daugther d0 and d1 are the same 
			    }
			  else d1tmp = GetNewLabel(id);
			}
		    }
		  c.SetDaughter(0,d0tmp);
		  c.SetDaughter(1,d1tmp);
		} else 
		{
		  AliFatal(Form("Unxpected indices %d %d",d0,d1));
		}
        
	      if ( m < 0 )
		{
		  c.SetMother(m);
		} else 
		{
		  if (IsParticleSelected(m)) 
		    {
		      c.SetMother(GetNewLabel(m));              
		    }
		  // else // FIXME: re-enable this checj. Sometimes it gets here. Still to be understood why
		  //   {
		  //     //		      AliError(Form("PROBLEM Mother not selected %d", m));              
		  //   }
		}
        
	      new ((*fMCParticles)[nmcout++]) AliAODMCParticle(c);
	    }
      
	  ++nmc;        
	} //closes loop over MC particles
    
      // now remap the tracks...
    
      TIter nextTrack(fTracks);
      AliNanoAODTrack* t;
      //      std::cout << "Remapping tracks" << std::endl;
    
      while ( ( t = dynamic_cast<AliNanoAODTrack*>(nextTrack()) ) )
	{
	  
	  t->SetLabel(GetNewLabel(t->GetLabel()));
	}
    
    } // closes fMCMode == 1
  else if ( mcParticles ) 
    {
      // simple copy of input MC particles to ouput MC particles
      TIter nextMC(mcParticles);
      AliAODMCParticle* p;
      Int_t nmcout(0);
    
      while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
	{
	  new ((*fMCParticles)[nmcout++]) AliAODMCParticle(*p);
	}
    }
  
  AliDebug(1,Form("input mc %d output mc %d",
                  mcParticles ? mcParticles->GetEntries() : 0,
                  fMCParticles ? fMCParticles->GetEntries() : 0));
  
  Printf("input mc %d output mc %d",
                  mcParticles ? mcParticles->GetEntries() : 0,
                  fMCParticles ? fMCParticles->GetEntries() : 0);
  

}

// //_____________________________________________________________________________
TList* AliNanoAODReplicator::GetList() const
{
  // return (and build if not already done) our internal list of managed objects
  
  if (!fList)
    {
      fList = new TList;
      fList->SetOwner(kTRUE);

      fTracks = new TClonesArray("AliNanoAODTrack");      
      fTracks->SetName("tracks"); // TODO: consider the possibility to use a different name to distinguish in AliAODEvent
      fList->Add(fTracks);    

      fHeader = new AliNanoAODHeader(3);// TODO: to be customized
      fHeader->SetName("header"); // TODO: consider the possibility to use a different name to distinguish in AliAODEvent
      fList->Add(fHeader);    


      fVertices = new TClonesArray("AliAODVertex",2);
      fVertices->SetName("vertices");    
    
        
      fList->Add(fVertices);
    
      if ( fMCMode > 0 )
	{
	  fMCHeader = new AliAODMCHeader;    
	  fMCParticles = new TClonesArray("AliAODMCParticle",1000);
	  fMCParticles->SetName(AliAODMCParticle::StdBranchName());
	  fList->Add(fMCHeader);
	  fList->Add(fMCParticles);
	}
    }
  return fList;
}

//_____________________________________________________________________________
void AliNanoAODReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  // Replicate (and filter if filters are there) the relevant parts we're interested in AODEvent
  

  // assert(fTracks!=0x0);
  
  //*fTZERO = *(source.GetTZEROData());
  
  

  fTracks->Clear("C");			
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  if (fMCMode > 0){
    if(!fMCHeader) {
      AliFatal(Form("fMCMode = %d, but MC header not found", fMCMode));
    }
    fMCHeader->Reset();
    if(!fMCParticles){
      AliFatal(Form("fMCMode = %d, but MC particles not found", fMCMode));
    }
    fMCParticles->Clear("C");
  }
  Int_t ntracks(0);
  Int_t input(0);

  AliAODVertex *vtx = source.GetPrimaryVertex();

  // TODO: implement header!
  //  *fHeader = *source.GetHeader();
  if(fCustomSetter){
    // Set custom variables in the header if the callback is set
    fCustomSetter->SetNanoAODHeader(&source, fHeader);
  }

  const Int_t entries = source.GetNumberOfTracks();
  if(entries<=0) return;

  for(Int_t j=0; j<entries; j++){
    
    AliVTrack *track = (AliVTrack*)source.GetTrack(j);
    
    AliAODTrack *aodtrack =(AliAODTrack*)track;// FIXME DYNAMIC CAST?
    if(fTrackCut && !fTrackCut->IsSelected(aodtrack)) continue;

    AliNanoAODTrack * special = new((*fTracks)[ntracks++]) AliNanoAODTrack (aodtrack, fVarList);
    
    if(fCustomSetter) fCustomSetter->SetNanoAODTrack(aodtrack, special);
  }  
  //----------------------------------------------------------
  
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) )
    {
      AliAODVertex* tmp = v->CloneWithoutRefs();
      AliAODVertex* copiedVertex = new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
      
      copiedVertex->SetNContributors(v->GetNContributors()); 
      
      delete tmp;
    }
  
  
  AliDebug(1,Form("input mu tracks=%d tracks=%d vertices=%d",
                  input,fTracks->GetEntries(),fVertices->GetEntries())); 
  
  
  // Finally, deal with MC information, if needed
  
  if ( fMCMode > 0 ) {
    FilterMC(source);      
  }
  

}



//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------
	

//-----------------------------------------------------------------------------

// AliAODVertex* AliNanoAODReplicator::PrimaryVertex(const TObjArray *trkArray,
// 						   AliAODEvent &event) const
// {
//   // Returns primary vertex to be used for this candidate
//   //AliCodeTimerAuto("",0);

//   AliESDVertex *vertexESD = 0;
//   AliAODVertex *vertexAOD = 0;


//   if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) { 
//     // primary vertex from the input event
    
//     vertexESD = new AliESDVertex(*fV1);

//   } else {
//     // primary vertex specific to this candidate

//     Int_t nTrks = trkArray->GetEntriesFast();
//     AliVertexerTracks *vertexer = new AliVertexerTracks(event.GetMagneticField());

//     if(fRecoPrimVtxSkippingTrks) { 
//       // recalculating the vertex
      
//       if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
// 	Float_t diamondcovxy[3];
// 	event.GetDiamondCovXY(diamondcovxy);
// 	Double_t pos[3]={event.GetDiamondX(),event.GetDiamondY(),0.};
// 	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
// 	AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
// 	vertexer->SetVtxStart(diamond);
// 	delete diamond; diamond=NULL;
// 	if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
// 	  vertexer->SetOnlyFitter();
//       }
//       Int_t skipped[1000];
//       Int_t nTrksToSkip=0,id;
//       AliExternalTrackParam *t = 0;
//       for(Int_t i=0; i<nTrks; i++) {
// 	t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
// 	id = (Int_t)t->GetID();
// 	if(id<0) continue;
// 	skipped[nTrksToSkip++] = id;
//       }
//       // TEMPORARY FIX
//       // For AOD, skip also tracks without covariance matrix
//       if(fInputAOD) {
// 	Double_t covtest[21];
// 	for(Int_t j=0; j<event.GetNumberOfTracks(); j++) {
// 	  AliVTrack *vtrack = (AliVTrack*)event.GetTrack(j);
// 	  if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
// 	    id = (Int_t)vtrack->GetID();
// 	    if(id<0) continue;
// 	    skipped[nTrksToSkip++] = id;
// 	  }
// 	}
//       }
//       for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
//       //
//       vertexer->SetSkipTracks(nTrksToSkip,skipped);
//       vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 
      
//     } else if(fRmTrksFromPrimVtx && nTrks>0) { 
//       // removing the prongs tracks
      
//       TObjArray rmArray(nTrks);
//       UShort_t *rmId = new UShort_t[nTrks];
//       AliESDtrack *esdTrack = 0;
//       AliESDtrack *t = 0;
//       for(Int_t i=0; i<nTrks; i++) {
// 	t = (AliESDtrack*)trkArray->UncheckedAt(i);
// 	esdTrack = new AliESDtrack(*t);
// 	rmArray.AddLast(esdTrack);
// 	if(esdTrack->GetID()>=0) {
// 	  rmId[i]=(UShort_t)esdTrack->GetID();
// 	} else {
// 	  rmId[i]=9999;
// 	}
//       }
//       Float_t diamondxy[2]={event.GetDiamondX(),event.GetDiamondY()};
//       vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
//       delete [] rmId; rmId=NULL;
//       rmArray.Delete();
      
//     }

//     if(!vertexESD) return vertexAOD;
//     if(vertexESD->GetNContributors()<=0) { 
//       //AliDebug(2,"vertexing failed"); 
//       delete vertexESD; vertexESD=NULL;
//       return vertexAOD;
//     }

//     delete vertexer; vertexer=NULL;

//   }

//   // convert to AliAODVertex
//   Double_t pos[3],cov[6],chi2perNDF;
//   vertexESD->GetXYZ(pos); // position
//   vertexESD->GetCovMatrix(cov); //covariance matrix
//   chi2perNDF = vertexESD->GetChi2toNDF();
//   delete vertexESD; vertexESD=NULL;

//   vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);

//   return vertexAOD;
// }

//_____________________________________________________________________________



// //---------------------------------------------------------------------------

void AliNanoAODReplicator::Terminate(){

}
