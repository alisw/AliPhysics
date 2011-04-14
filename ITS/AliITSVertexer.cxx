#include "AliLog.h"
#include "AliMultiplicity.h"
#include "AliITSgeomTGeo.h"
#include "AliITSVertexer.h"
#include "AliITSLoader.h"
#include "AliITSMultReconstructor.h"
#include "AliITSRecPointContainer.h"
#include "AliRunLoader.h"

const Float_t AliITSVertexer::fgkPipeRadius = 3.0;

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSvertexer3D, AliITSVertexerZ.              //       
//                  AliITSVertexerCosmics                           //
//////////////////////////////////////////////////////////////////////

/* $Id$ */

//______________________________________________________________________
AliITSVertexer::AliITSVertexer():AliVertexer(),
fLadders(NULL), 
fLadOnLay2(0),
fComputeMultiplicity(kFALSE),
fDetTypeRec(NULL),
fMinTrackletsForPilup(0),
fIsPileup(0),
fNTrpuv(-2),
fZpuv(-9999999.),
fNoVertices(0),
fVertArray(NULL),
fFirstEvent(0),
fLastEvent(-1)
{
  // Default Constructor
  SetLaddersOnLayer2();
  SetMinTrackletsForPilup();
  for(Int_t i=0; i<kNSPDMod;i++) fUseModule[i]=kTRUE;
}

//______________________________________________________________________
AliITSVertexer::~AliITSVertexer() {
  // Destructor
  if(fLadders) delete [] fLadders;
  if (fNoVertices > 0){
    delete []fVertArray;
    fVertArray = NULL;
    fNoVertices = 0;
  }
}

//______________________________________________________________________
void AliITSVertexer::ResetVertex(){
  // Resets vertex related data members
  if(fNoVertices > 0){
    if(fVertArray) delete []fVertArray;
    fVertArray = NULL;
    fNoVertices = 0;
  }
  fIsPileup=kFALSE;
  fNTrpuv=-2;
  fZpuv=-99999.;

}
//______________________________________________________________________
void AliITSVertexer::FindMultiplicity(TTree *itsClusterTree){
  // Invokes AliITSMultReconstructor to determine the
  // charged multiplicity in the pixel layers
  if(fMult){delete fMult; fMult = 0;}

  Bool_t success=kTRUE;
  Bool_t cosmics=kFALSE; 
  if(!fCurrentVertex)success=kFALSE;
  if(fCurrentVertex && fCurrentVertex->GetNContributors()<1)success=kFALSE;
  if(fCurrentVertex && strstr(fCurrentVertex->GetTitle(),"cosmics")) {
    success=kFALSE; 
    cosmics=kTRUE;
  } 

  // get the FastOr bit mask
  TBits fastOrFiredMap = fDetTypeRec->GetFastOrFiredMap();
  TBits firedChipMap = fDetTypeRec->GetFiredChipMap(itsClusterTree);
 
  AliITSMultReconstructor multReco;

  if(!success){
    if(!cosmics) {     
      AliDebug(1,"Tracklets multiplicity not determined because the primary vertex was not found");
      AliDebug(1,"Just counting the number of cluster-fired chips on the SPD layers");
    }
    if (!itsClusterTree) {
      AliError(" Invalid ITS cluster tree !\n");
      return;
    }
    multReco.LoadClusterFiredChips(itsClusterTree);
    Short_t nfcL1 = multReco.GetNFiredChips(0);
    Short_t nfcL2 = multReco.GetNFiredChips(1);
    fMult = new AliMultiplicity(0,0,0,0,0,0,0,0,0,0,0,nfcL1,nfcL2,fastOrFiredMap);
    fMult->SetFiredChipMap(firedChipMap);
    AliITSRecPointContainer* rcont = AliITSRecPointContainer::Instance();
    fMult->SetITSClusters(0,rcont->GetNClustersInLayer(1,itsClusterTree));
    for(Int_t kk=2;kk<=6;kk++){
      fMult->SetITSClusters(kk-1,rcont->GetNClustersInLayerFast(kk));
    }
    return;
  }

  if (!itsClusterTree) {
    AliError(" Invalid ITS cluster tree !\n");
    return;
  }
  Double_t vtx[3];
  fCurrentVertex->GetXYZ(vtx);
  Float_t vtxf[3];
  for(Int_t i=0;i<3;i++)vtxf[i]=vtx[i];
  multReco.SetHistOn(kFALSE);
  multReco.Reconstruct(itsClusterTree,vtxf,vtxf);
  Int_t notracks=multReco.GetNTracklets();
  Float_t *tht = new Float_t [notracks];
  Float_t *phi = new Float_t [notracks];
  Float_t *dtht = new Float_t [notracks];
  Float_t *dphi = new Float_t [notracks];
  Int_t *labels = new Int_t[notracks];
  Int_t *labelsL2 = new Int_t[notracks];
  for(Int_t i=0;i<multReco.GetNTracklets();i++){
    tht[i] = multReco.GetTracklet(i)[0];
    phi[i] =  multReco.GetTracklet(i)[1];
    dtht[i] = multReco.GetTracklet(i)[3];
    dphi[i] = multReco.GetTracklet(i)[2];
    labels[i] = static_cast<Int_t>(multReco.GetTracklet(i)[4]);
    labelsL2[i] = static_cast<Int_t>(multReco.GetTracklet(i)[5]);
  }
  Int_t nosingleclus=multReco.GetNSingleClusters();
  Float_t *ths = new Float_t [nosingleclus];
  Float_t *phs = new Float_t [nosingleclus];
  Int_t *labelss = new Int_t [nosingleclus];
  for(Int_t i=0;i<nosingleclus;i++){
    ths[i] = multReco.GetCluster(i)[0];
    phs[i] = multReco.GetCluster(i)[1];
    labelss[i] = (Int_t)multReco.GetCluster(i)[2];
  }
  Short_t nfcL1 = multReco.GetNFiredChips(0);
  Short_t nfcL2 = multReco.GetNFiredChips(1);
  fMult = new AliMultiplicity(notracks,tht,phi,dtht,dphi,labels,labelsL2,nosingleclus,ths,phs,labelss,nfcL1,nfcL2,fastOrFiredMap);
  fMult->SetFiredChipMap(firedChipMap);
  AliITSRecPointContainer* rcont = AliITSRecPointContainer::Instance();
  fMult->SetITSClusters(0,rcont->GetNClustersInLayer(1,itsClusterTree));
  for(Int_t kk=2;kk<=6;kk++){
    fMult->SetITSClusters(kk-1,rcont->GetNClustersInLayerFast(kk));
  }
  delete [] tht;
  delete [] phi;
  delete [] dtht;
  delete [] dphi;
  delete [] ths;
  delete [] phs;
  delete [] labels;
  delete [] labelsL2;
  delete [] labelss;

  return;
}

//______________________________________________________________________
void AliITSVertexer::SetLaddersOnLayer2(Int_t ladwid){
  // Calculates the array of ladders on layer 2 to be used with a 
  // given ladder on layer 1
  if(ladwid == fLadOnLay2 && fLadders)return;
  fLadOnLay2=ladwid;
  Int_t ladtot1=AliITSgeomTGeo::GetNLadders(1);
  if(fLadders) delete [] fLadders;
  fLadders=new UShort_t[ladtot1];


  Double_t pos1[3],pos2[3];
  Int_t mod1=AliITSgeomTGeo::GetModuleIndex(2,1,1);
  AliITSgeomTGeo::GetTranslation(mod1,pos1);  // position of the module in the MRS 
  Double_t phi0=TMath::ATan2(pos1[1],pos1[0]);
  if(phi0<0) phi0+=2*TMath::Pi();
  Int_t mod2=AliITSgeomTGeo::GetModuleIndex(2,2,1);
  AliITSgeomTGeo::GetTranslation(mod2,pos2);
  Double_t phi2=TMath::ATan2(pos2[1],pos2[0]); 
  if(phi2<0) phi2+=2*TMath::Pi();
  Double_t deltaPhi= phi0-phi2; // phi width of a layer2 module

  for(Int_t i= 0; i<ladtot1;i++){
    Int_t modlad= AliITSgeomTGeo::GetModuleIndex(1,i+1,1);
    Double_t posmod[3];
    AliITSgeomTGeo::GetTranslation(modlad,posmod);
    Double_t phimod=TMath::ATan2(posmod[1],posmod[0]); 
    if(phimod<0) phimod+=2*TMath::Pi();
    Double_t phi1= phimod+deltaPhi*double(fLadOnLay2);
    if(phi1<0) phi1+=2*TMath::Pi();
    if(phi1>2*TMath::Pi()) phi1-=2*TMath::Pi();
    Double_t philad1=phi0-phi1;
    UShort_t lad1;
    Double_t ladder1=(philad1)/(deltaPhi) +1.; 
    if(ladder1<1){ladder1=40+ladder1;}
    lad1=int(ladder1+0.5);
    fLadders[i]=lad1;
  }
}


//______________________________________________________________________
void AliITSVertexer::Init(TString filename){
  // Initialize the vertexer in case of
  // analysis of an entire file
  AliRunLoader *rl = AliRunLoader::Instance();
  if(!rl){
    AliFatal("Run Loader not found");
    return;
  }
  if (fLastEvent < 0) SetLastEvent(rl->GetNumberOfEvents()-1);

  AliITSLoader* itsloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!filename.Contains("default"))itsloader->SetVerticesFileName(filename);
  if(!filename.Contains("null"))itsloader->LoadVertices("recreate");
}

//______________________________________________________________________
void AliITSVertexer::WriteCurrentVertex(){
  // Write the current AliVertex object to file fOutFile
  AliRunLoader *rl = AliRunLoader::Instance();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  fCurrentVertex->SetName("Vertex");
  //  const char * name = fCurrentVertex->GetName();
  //  itsLoader->SetVerticesContName(name);
  Int_t rc = itsLoader->PostVertex(fCurrentVertex);
  rc = itsLoader->WriteVertices();
}

//______________________________________________________________________
void AliITSVertexer::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent

  AliRunLoader *rl = AliRunLoader::Instance();
  AliITSLoader* itsloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsloader->LoadRecPoints("read");
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    rl->GetEvent(i);
    TTree* cltree = itsloader->TreeR();
    FindVertexForCurrentEvent(cltree);
    if(fCurrentVertex){
      WriteCurrentVertex();
    }
    else {
      AliDebug(1,Form("Vertex not found for event %d",i));
    }
  }
}
