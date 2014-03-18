/*************************************************************************
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

//
// Implementation of a branch replicator 
// to produce aods with only few branches.
//
// This replicator is in charge of replicating the nuclei primary vertices
// tracks identified as nuclei with Z>=2, secondary vertices in form of 
// AliAODRecoDecayLF2Prong and their daughter tracks.
// These informations are stored into a reduced AODs (AliAOD.NuclEx.root) 
// 
// The vertices are filtered so that only the primary vertices make it
// to the output aods.
//
// The secondary vertices are recreated here, as a AliAODRecoDecayLF2Prong
// plus cuts that select secondary vericesoutside the primary vertex

// Authors: S. Bufalino (stefania.bufalino@cern.ch)
//          R. Lea      (ramona.lea@cern.ch)
// Based on AliAODMuonReplicator.cxx 

//NOTE : This is a test on MC : no PID response only MC truth + select only 3LH
// daughters : NOT INTENDED for any efficiency!!

class AliESDv0;
class AliESDVertex;
class AliAODVertex;
class AliAODRecoDecay;

#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTZERO.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
//#include "AliAnalysisCuts.h"
#include "TF1.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
//#include "AliPIDResponse.h"
#include <iostream>
#include <cassert>
#include "AliESDtrack.h"
#include "TObjArray.h"
#include "AliAnalysisFilter.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayLF.h"
#include "AliAODRecoDecayLF2Prong.h"

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
//#include "AliAnalysisVertexingLF.h"
#include "AliAnalysisManager.h"
#include "AliAODMCNuclExReplicator.h"
#include "TH1.h"
#include "TCanvas.h"
#include "AliInputEventHandler.h"

using std::cout;
using std::endl;

ClassImp(AliAODMCNuclExReplicator)

//_____________________________________________________________________________
AliAODMCNuclExReplicator::AliAODMCNuclExReplicator(const char* name, const char* title,
						   Int_t mcMode, 
						   Int_t nsigmaTrk1, Int_t partType1,
						   Int_t nsigmaTrk2, Int_t partType2
						   )
  :AliAODBranchReplicator(name,title),
  fBzkG(0.),
  fCosMin(),
  fDCAtracksMin(),
  fRmax(),
  fRmin(),
  fDNmin(),
  fDPmin(), 
  fHeader(0x0),
  fVertices(0x0), 
  fNuclei(0x0),
  fSecondaryVerices(0x0), 
  fDaughterTracks(0x0),
  fList(0x0),
  fMCParticles(0x0),
  fMCHeader(0x0),
  fMCMode(mcMode),
  fLabelMap(),
  fParticleSelected(),
  fReplicateHeader(kTRUE), //replicateHeader //to be fixed
  fnSigmaTrk1(nsigmaTrk1),
  fnSigmaTrk2(nsigmaTrk2),
  fpartType1(partType1),
  fpartType2(partType2),
  fSecVtxWithKF(kFALSE),
  fVertexerTracks(0x0),
  fV1(0x0),
  fAODMapSize(0),
  fAODMap(0)
  
{
  // default ctor
}

//_____________________________________________________________________________
AliAODMCNuclExReplicator::~AliAODMCNuclExReplicator()
{
  // destructor
  // delete fTrackCut;
  // delete fVertexCut;
  delete fList;
}

//_____________________________________________________________________________
void AliAODMCNuclExReplicator::SelectParticle(Int_t i)
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
Bool_t AliAODMCNuclExReplicator::IsParticleSelected(Int_t i)  
{
  // taking the absolute values here, need to take 
  // care with negative daughter and mother
  // IDs when setting!
  return (fParticleSelected.GetValue(TMath::Abs(i))==1);
}


//_____________________________________________________________________________
void AliAODMCNuclExReplicator::CreateLabelMap(const AliAODEvent& source)
{  
  //
  // this should be called once all selections are done 
  //
  
  fLabelMap.Delete();
  
  TClonesArray* mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  Int_t i(0);
  Int_t j(0);
  
  TIter next(mcParticles);
  
  while ( next() ) 
  {
    if (IsParticleSelected(i))
    {
      fLabelMap.Add(i,j++);
    }
    ++i;
  }
}

//_____________________________________________________________________________
Int_t AliAODMCNuclExReplicator::GetNewLabel(Int_t i) 
{
  // Gets the label from the new created Map
  // Call CreatLabelMap before
  // otherwise only 0 returned
  return fLabelMap.GetValue(TMath::Abs(i));
}


//_____________________________________________________________________________
TList* AliAODMCNuclExReplicator::GetList() const
{
  // return (and build if not already done) our internal list of managed objects
  if (!fList)
    {
    
      if ( fReplicateHeader )
	{
	  fHeader = new AliAODHeader;
	}
      
      
      fSecondaryVerices = new TClonesArray("AliAODRecoDecayLF2Prong",30);
      fSecondaryVerices->SetName("SecondaryVertices");    
    
      fVertices = new TClonesArray("AliAODVertex",2);
      fVertices->SetName("vertices");    
    
      fNuclei = new TClonesArray("AliAODTrack",30);
      fNuclei->SetName("Nuclei");
    
      fDaughterTracks = new TClonesArray("AliAODTrack",30);
      fDaughterTracks->SetName("DaughterTracks");

    
      fList = new TList;
      fList->SetOwner(kTRUE);
  
      fList->Add(fHeader);
      fList->Add(fVertices);
      fList->Add(fNuclei);
      fList->Add(fSecondaryVerices);
      fList->Add(fDaughterTracks);
  
      
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
void AliAODMCNuclExReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  // Replicate (and filter if filters are there) the relevant parts we're interested in AODEvent
  
  //  cout<<"-------------------->QUESTO"<<endl;

  //-----------------------------------------------
  // AliPIDResponse
  
  // AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    
  //--------------------------------------------------------

  //  printf("Execute NuclEx Replicator\n");

  //---------------------------------

  if (fReplicateHeader)
    {
      *fHeader = *(source.GetHeader());
    }
    
  fVertices->Clear("C");			
    
  fNuclei->Clear("C");

  fSecondaryVerices->Clear("C");

  fDaughterTracks->Clear("C");
  
  //----------------------------------
  
  //retrive MC infos
  
  TClonesArray *arrayMC = 0;
  AliAODMCHeader *mcHeader=0;
  Int_t mumpdg=-100;
 
  arrayMC = (TClonesArray*) source.GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!arrayMC) {
    Printf("Error: MC particles branch not found!\n");
    return;
  }
  // if(arrayMC)
  //   cout<<"Ho caricato array mc"<<endl;
  
  mcHeader =  (AliAODMCHeader*)source.GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSEDs::UserExec: MC header branch not found!\n");
    return;
  } 
  // if(mcHeader)
  //   cout<<"Ho caricato MC header"<<endl;

  //  cout<<"Centrality AOD source: "<<source.GetHeader()->GetCentrality()<<endl;

  Int_t nsv(0);
  //  Int_t nnuclei(0);
  Int_t ntracksD(0);

  Int_t input(0);
  Double_t xdummy,ydummy;

  AliAODRecoDecayLF2Prong *io2Prong  = 0;

  TObjArray *twoTrackArray    = new TObjArray(2);
  Double_t dispersion;

  // cout<<"Qui"<<endl;
  //cout<<source.GetMagneticField()<<endl;

  AliAODVertex *vtx = source.GetPrimaryVertex();
						
  //  cout<<"Source "<<source<<endl;
  //cout<<"vtx: "<<vtx<<endl;

  // A Set of very loose cut for a weak strange decay
  
  fCosMin       = 0.99;
  fDCAtracksMin = 1;
  fRmax         = 200.;
  fRmin         = 0.1;
  fDNmin        = 0.05;
  fDPmin        = 0.05;

  //----------------------------------------------------------

  //  Int_t nindices=0;
  UShort_t *indices = 0;
  const Int_t entries = source.GetNumberOfTracks();

  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  fV1 = new AliESDVertex(pos,cov,100.,100,vtx->GetName());
  //  cout<<"fV1 pricipal loop: "<<fV1<<endl;
  
  if(entries<=0) return;
  indices = new UShort_t[entries];
  memset(indices,0,sizeof(UShort_t)*entries);
  fAODMapSize = 100000;
  fAODMap = new Int_t[fAODMapSize];
  memset(fAODMap,0,sizeof(Int_t)*fAODMapSize);
  //  cent=((AliAODEvent&)source)->GetCentrality();
  
  //-------------------------------------------------------------

  if(vtx->GetNContributors()<1) {
    
    // SPD vertex cut
    vtx =source.GetPrimaryVertexSPD();
    
    if(vtx->GetNContributors()<1) {
      Info("AliAnalysisTaskHelium3Pi","No good vertex, skip event");
      return; // NO GOOD VERTEX, SKIP EVENT 
    }
  }
  
  Double_t xPrimaryVertex=0.,yPrimaryVertex=0.;
  xPrimaryVertex=vtx->GetX();
  yPrimaryVertex=vtx->GetY();
  
  fBzkG=source.GetMagneticField();
  fVertexerTracks=new AliVertexerTracks(fBzkG);

  Double_t TrackNumber = source.GetNumberOfTracks();
  Int_t label =-1;

  //Tracks arrays
  
  TArrayI Track0(TrackNumber);        //Pions                                                                          
  Int_t nTrack0=0;
  
  TArrayI Track1(TrackNumber);        //Helium3
  Int_t nTrack1=0;

  for(Int_t j=0; j<TrackNumber; j++){

    //    cout<<"Inside loop tracks"<<endl;

  
    AliVTrack *track = (AliVTrack*)source.GetTrack(j);
    
    AliAODTrack *aodtrack =(AliAODTrack*)track;

    //-----------------------------------------------------------
    //Track cuts 
    if(aodtrack->GetTPCNcls() < 70 )continue;
    if(aodtrack->Chi2perNDF() > 4 )continue;
    
    if (!aodtrack->IsOn(AliAODTrack::kTPCrefit)) continue;
    if (!aodtrack->IsOn(AliAODTrack::kTPCin)) continue;
    if (aodtrack->IsOn(AliAODTrack::kITSpureSA)) continue;

    //---------------------------------------------------------------
     
    Double_t mom = aodtrack->P();
    
    if(mom<0.150)continue;

    label = TMath::Abs(aodtrack->GetLabel());

    AliAODMCParticle *part = (AliAODMCParticle*) arrayMC->At(label);
    
    Int_t PDGCode=part->GetPdgCode();
    
    Int_t mumid = part->GetMother();
    
    if(mumid>-1){
      AliAODMCParticle *mother = (AliAODMCParticle*) arrayMC->At(mumid);
      mumpdg = mother->GetPdgCode();
    }
    
    //    if(mumpdg == 1010010030 ||mumpdg == -1010010030 ){

    if(mumpdg == 1010010030){

      //	if(PDGCode==-211 || PDGCode==+211){  
      if(PDGCode==-211){  
	Track0[nTrack0++]=j;
      }
      
      //	if(PDGCode==1000020030 ||PDGCode==-1000020030 ){
      if(PDGCode==1000020030){
	Track1[nTrack1++]=j;
	// 	new((*fNuclei)[nnuclei++]) AliAODTrack(*aodtrack);
      }
    }
  }
  
  //Set Track Daughters
  Track0.Set(nTrack0);
  Track1.Set(nTrack1);
  
  
  // cout<<"Track loops..."<<endl;
  // cout<<"npos "<<nTrack1<<endl;
  // cout<<"nneg "<<nTrack0<<endl;


  AliAODTrack *postrackAOD = 0;
  AliAODTrack *negtrackAOD = 0;
 
  AliESDtrack *postrack = 0;
  AliESDtrack *negtrack = 0;

  Bool_t isOk=kFALSE;

  for (Int_t i=0; i<nTrack1; i++){                            //! He Tracks Loop
    
    Int_t Track1idx=Track1[i];

    AliVTrack *trackpos = (AliVTrack*)source.GetTrack(Track1idx);

    postrackAOD = (AliAODTrack*)trackpos;
    postrack = new AliESDtrack(trackpos);

    //--- MC infos

    Int_t labelpos = TMath::Abs(postrack->GetLabel());
    AliAODMCParticle *partPos = (AliAODMCParticle*) arrayMC->At(labelpos);
    Int_t mumidPos = partPos->GetMother();

    //------------------------------
    
    for (Int_t k=0; k <nTrack0 ; k++) {                           //! Pion Tracks Loop
      
      Int_t Track0idx=Track0[k];
      
      AliVTrack *trackneg = (AliVTrack*)source.GetTrack(Track0idx);
      negtrackAOD =(AliAODTrack*)trackneg;
      negtrack = new AliESDtrack(trackneg); 
      
      //--- MC infos
      
      Int_t labelneg = TMath::Abs(negtrack->GetLabel());
      AliAODMCParticle *partNeg = (AliAODMCParticle*) arrayMC->At(labelneg);
      Int_t mumidNeg = partNeg->GetMother();
      
      
      //------------------------------
      //  if(mumidPos == mumidNeg && mumidNeg > 0){
      isOk=kFALSE;
      
      if(mumidPos == mumidNeg && mumidNeg > 0)
	isOk = kTRUE;
      
      twoTrackArray->AddAt(negtrack,0);
      twoTrackArray->AddAt(postrack,1);
      
      Double_t dcap1n1 = postrack->GetDCA(negtrack,fBzkG,xdummy,ydummy);
      
      Double_t dcap1toPV = TMath::Abs(postrack->GetD(xPrimaryVertex, yPrimaryVertex,fBzkG));
      Double_t dcan1toPV = TMath::Abs(negtrack->GetD(xPrimaryVertex, yPrimaryVertex,fBzkG));
      
      if(dcap1n1>fDCAtracksMin)continue;
      if((xdummy+ydummy)>fRmax )continue;
      if((xdummy+ydummy)< fRmin)continue;
	
      if ( dcan1toPV < fDNmin)               
	if ( dcap1toPV < fDPmin) continue;   
      
      //      cout<<"dcap1n1: "<<dcap1n1<<endl;
	
      AliAODVertex *vertexp1n1 = ReconstructSecondaryVertex(twoTrackArray,dispersion,kTRUE);
	
      if(!vertexp1n1) {
	  
	twoTrackArray->Clear();
	delete negtrack;
	negtrack=NULL; 
	continue; 
      }
	
      io2Prong = Make2Prong(twoTrackArray,source,vertexp1n1,dcap1n1);
	
      if(io2Prong->CosPointingAngle()<fCosMin)continue;
	
      
      // AliAODTrack *trk0 = (AliAODTrack*)io2Prong->GetDaughter(0);
      // AliAODTrack *trk1 = (AliAODTrack*)io2Prong->GetDaughter(1);
	
      // cout<<"**********************************************"<<endl;
      // cout<<trk0/*->GetID()*/<<" "<<negtrackAOD->GetID()<<endl;
      // cout<<trk1/*->GetID()*/<<" "<<postrackAOD->GetID()<<endl;
      // cout<<"d0 io2Prong: "<<io2Prong->GetProngID(1)<<endl;
      // cout<<"d1 io2Prong: "<<io2Prong->GetProngID(0)<<endl;
      // cout<<"**********************************************"<<endl;
	
      //      rd =  new((*fSecondaryVerices)[nsv++]) AliAODRecoDecayLF2Prong(*io2Prong);
      if(isOk){
	new((*fSecondaryVerices)[nsv++]) AliAODRecoDecayLF2Prong(*io2Prong);
	new((*fDaughterTracks)[ntracksD++]) AliAODTrack(*negtrackAOD);
	new((*fDaughterTracks)[ntracksD++]) AliAODTrack(*postrackAOD);
	
      }
      // rd->SetSecondaryVtx(vertexp1n1);
      // vertexp1n1->SetParent(rd);
      // AddRefs(vertexp1n1,rd,source,twoTrackArray);
      
      delete negtrack; 
      negtrack = NULL;
      
      delete vertexp1n1;
      vertexp1n1= NULL;
      continue;
      //   }
    }
    
    delete postrack; 
    postrack = NULL;
    
  }
  
  //----------------------------------------------------------
 
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) )
    {
      if ( v->GetType() == AliAODVertex::kPrimary     ||
	   v->GetType() == AliAODVertex::kMainSPD     ||
	   v->GetType() == AliAODVertex::kPileupSPD   ||
	   v->GetType() == AliAODVertex::kPileupTracks||
	   v->GetType() == AliAODVertex::kMainTPC  ) 
	{
	  	  
	  AliAODVertex* tmp = v->CloneWithoutRefs();
	  //AliAODVertex* copiedVertex = new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
	  new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
	  
	  // to insure the main vertex retains the ncontributors information
	  // (which is otherwise computed dynamically from
	  // references to tracks, which we do not keep in muon aods...)
	  // we set it here
	  
	  //copiedVertex->SetNContributors(v->GetNContributors()); 
	  
	  //  fVertices->Delete();
	  //	  delete copiedVertex;
	  delete tmp;
	  //	  printf("....Prendo il vertice primario...\n");
	}
      //       printf("....Prendo il vertice primario...\n");
    }
  
  //printf("....Done NuclEx Replicator...\n");
  
  AliDebug(1,Form("input mu tracks=%d tracks=%d vertices=%d nnuclei=%d",
                  input,fSecondaryVerices->GetEntries(),fVertices->GetEntries(),fNuclei->GetEntries())); 
  
  // cout<<"Delete..."<<endl;
  // delete foPion;
  // foPion = NULL;
  //cout<<"Delete 1"<<  endl;
  
  if(io2Prong) {delete io2Prong; io2Prong=NULL;}
  //cout<<"Delete 2"<<  endl;
  twoTrackArray->Delete();  delete twoTrackArray;
  //  cout<<"Delete 3"<<  endl;
  // vtx->Delete();  delete vtx;
  //  cout<<"Delete 4"<<  endl;
  if(fV1) { delete fV1; fV1=NULL; }
  //  cout<<"Delete 5"<<  endl;
  if(fAODMap) { delete [] fAODMap; fAODMap=NULL; }
  delete []indices;
  //  cout<<"Delete 6"<<  endl;
  delete fVertexerTracks;

  //  cout<<"Mi rompo alla fine. Perche???"<<endl;
  
} 

//-----------------------------------------------------------------------------

AliAODVertex *AliAODMCNuclExReplicator::ReconstructSecondaryVertex(TObjArray *trkArray,
								Double_t &dispersion,Bool_t useTRefArray) const

{
  // Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle
  //AliCodeTimerAuto("",0);

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  if(!fSecVtxWithKF) { // AliVertexerTracks

    fVertexerTracks->SetVtxStart(fV1);
    vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(trkArray);

    if(!vertexESD) return vertexAOD;


    Double_t vertRadius2=vertexESD->GetXv()*vertexESD->GetXv()+vertexESD->GetYv()*vertexESD->GetYv();
    if(vertRadius2>200.){
      // vertex outside beam pipe, reject candidate to avoid propagation through material
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }

  } else { // Kalman Filter vertexer (AliKFParticle)

    AliKFParticle::SetField(fBzkG);

    AliKFVertex vertexKF;

    Int_t nTrks = trkArray->GetEntriesFast();
    for(Int_t i=0; i<nTrks; i++) {
      AliESDtrack *esdTrack = (AliESDtrack*)trkArray->At(i);
      AliKFParticle daughterKF(*esdTrack,211);
      vertexKF.AddDaughter(daughterKF);
    }
    vertexESD = new AliESDVertex(vertexKF.Parameters(),
				 vertexKF.CovarianceMatrix(),
				 vertexKF.GetChi2(),
				 vertexKF.GetNContributors());

  }
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; 
  vertexESD=NULL;

  Int_t nprongs= (useTRefArray ? 0 : trkArray->GetEntriesFast());
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);

  //  cout<<"------------------> Reconstruct vertexAOD: "<<vertexAOD<<endl;

  return vertexAOD;

 
}

//-----------------------------------------------------------------------------

AliAODRecoDecayLF2Prong* AliAODMCNuclExReplicator::Make2Prong(TObjArray *twoTrackArray,const AliAODEvent &evento,
							   AliAODVertex *secVert,Double_t dca)


{ 

  //cout<<"Inside make 2 prong"<<endl;

  Int_t charge[2];
  charge[0]=1; //it was -1 //Ramona
  charge[1]=2;
      
  // From AliAnalysisVertexingLF.cxx Method:Make2Prongs
  
  fBzkG = evento.GetMagneticField();
      
  Double_t px[2],py[2],pz[2],d0[2],d0err[2];
  // Also this has been changed //Ramona
  AliESDtrack *negtrack = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack *postrack = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  //cout<<"negtrack: "<<negtrack<<" postrack: "<<postrack<<endl;
  //cout<<"kVeryBig: "<<kVeryBig<<endl;
  //cout<<"secVert: "<<secVert<<endl;

  // // propagate tracks to secondary vertex, to compute inv. mass
  
  postrack->PropagateToDCA(secVert,fBzkG,kVeryBig);
  negtrack->PropagateToDCA(secVert,fBzkG,kVeryBig);
  
  Double_t momentum[3];
  
  negtrack->GetPxPyPz(momentum);
  px[0] = charge[0]*momentum[0]; 
  py[0] = charge[0]*momentum[1]; 
  pz[0] = charge[0]*momentum[2]; 

  postrack->GetPxPyPz(momentum);
  px[1] = charge[1]*momentum[0]; 
  py[1] = charge[1]*momentum[1]; 
  pz[1] = charge[1]*momentum[2]; 
  
  //cout<< px[0] <<" "<< " "<< py[0] << " "<< pz[0]<<endl;
  //  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2]; 
  //cout<< px[1] <<" "<< " "<< py[1] << " "<< pz[1]<<endl;
  //px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2]; 
  

  // primary vertex to be used by this candidate
  AliAODVertex *primVertexAOD  = evento.GetPrimaryVertex();
  //cout<<"primVertexAOD "<<primVertexAOD<<endl;
  if(!primVertexAOD) return 0x0;
      
  Double_t d0z0[2],covd0z0[3];

  d0z0[0] = -999.;
  d0z0[1] = -999.;
  covd0z0[0] = -999.;
  covd0z0[1] = -999.;
  covd0z0[2] = -999.;

  d0[0] = d0z0[0];
  d0err[0] = TMath::Sqrt(TMath::Abs(covd0z0[0]));

  d0[1] = d0z0[0];
  d0err[1] = TMath::Sqrt(TMath::Abs(covd0z0[0]));
  
  // create the object AliAODRecoDecayLF2Prong
  //  AliAODRecoDecayLF2Prong *the2Prong = new AliAODRecoDecayLF2Prong(secVert,px,py,pz,d0,d0err,dcap1n1);
  AliAODRecoDecayLF2Prong *the2Prong = new AliAODRecoDecayLF2Prong(secVert,px,py,pz,d0,d0err,dca);
  the2Prong->SetOwnPrimaryVtx(primVertexAOD);
  
  //  the2Prong->SetProngIDs(2,{-999,999});
  // UShort_t id0[2]={99999,999999};
  // the2Prong->SetProngIDs(2,id0);

  UShort_t id[2]={(UShort_t)postrack->GetID(),(UShort_t)negtrack->GetID()};
  the2Prong->SetProngIDs(2,id);

  // cout<<"\n\n\nMake 2 Prong: id[0]"<<id[0]<<" id[1]: "<<id[1]<<endl;
  // cout<<"Get: 1 "<<the2Prong->GetProngID(0)<<" 2 "<<the2Prong->GetProngID(1)<<endl;
  // cout<<"Get: 3 "<<the2Prong->GetProngID(2)<<" 4 "<<the2Prong->GetProngID(3)<<endl<<endl<<endl;
  //delete primVertexAOD; primVertexAOD=NULL;
  
  if(postrack->Charge()!=0 && negtrack->Charge()!=0) {
      
    AddDaughterRefs(secVert,(AliAODEvent&)evento,twoTrackArray);
    //    AddDaughterRefs(secVert,(AliAODEvent*)evento,twoTrackArray);
      
  }
  
  return the2Prong;  

  delete the2Prong;
}

//----------------------------------------------------------------------------
void AliAODMCNuclExReplicator::AddDaughterRefs(AliAODVertex *v,
					    const AliAODEvent &event,
					    const TObjArray *trkArray) const

{
  // Add the AOD tracks as daughters of the vertex (TRef)
  // AliCodeTimerAuto("",0);
  // cout<<"Inside  AddDaughterRefs "<<endl;

  Int_t nDg = v->GetNDaughters();
  
  TObject *dg = 0;
  if(nDg) dg = v->GetDaughter(0);
  if(dg) return; // daughters already added
  
  Int_t nTrks = trkArray->GetEntriesFast();
  
  AliExternalTrackParam *track = 0;
  AliAODTrack *aodTrack = 0;
  Int_t id;
  
  for(Int_t i=0; i<nTrks; i++) {
    track = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
  
    id = (Int_t)track->GetID();
    //    printf("---> %d\n",id);
    if(id<0) continue; // this track is a AliAODRecoDecay
  
    aodTrack = (AliAODTrack*)event.GetTrack(fAODMap[id]);
    v->AddDaughter(aodTrack);
  }
  return;
}
//----------------------------------------------------------------------------
	
void AliAODMCNuclExReplicator::AddRefs(AliAODVertex *v,AliAODRecoDecayLF *rd,
				    const AliAODEvent &event,
				    const TObjArray *trkArray) const

{
  // Add the AOD tracks as daughters of the vertex (TRef)
  // Also add the references to the primary vertex and to the cuts
  //AliCodeTimerAuto("",0);
  
  AddDaughterRefs(v,event,trkArray);
  rd->SetPrimaryVtxRef((AliAODVertex*)event.GetPrimaryVertex());
  return;
}	

//---------------------------------------------------------------------------

void AliAODMCNuclExReplicator::Terminate(){

}
