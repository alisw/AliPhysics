#include <AliESDVertex.h>
#include "AliITSgeomTGeo.h"
#include "AliITSVertexer.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliMultiplicity.h"
#include "AliITSMultReconstructor.h"

const Float_t AliITSVertexer::fgkPipeRadius = 3.0;

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexerIons AliITSvertexer3D             //
//                  AliITSVertexerCosmics                           //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSVertexer::AliITSVertexer():AliVertexer(),
fLadders(), 
fLadOnLay2(0)	 {
  // Default Constructor
  SetLaddersOnLayer2();
}

AliITSVertexer::AliITSVertexer(TString filename):AliVertexer(),
fLadders(), 
fLadOnLay2(0)
{
  // Standard constructor
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  if(!rl){
    Fatal("AliITSVertexer","Run Loader not found");
  }
  /*
  if(rl->LoadgAlice()){
    Fatal("AliITSVertexer","The AliRun object is not available - nothing done");
  }
  */
  fCurrentVertex  = 0;   
  SetFirstEvent(0);
  SetLastEvent(0);
  //  rl->LoadHeader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!filename.Contains("default"))itsLoader->SetVerticesFileName(filename);
  if(!filename.Contains("null"))itsLoader->LoadVertices("recreate");
  itsLoader->LoadRecPoints();
  //  Int_t lst;
  SetLastEvent(rl->GetNumberOfEvents()-1);
  /*
  if(rl->TreeE()){
    lst = static_cast<Int_t>(rl->TreeE()->GetEntries());
    SetLastEvent(lst-1);
  }
  */
  SetLaddersOnLayer2();
}

//______________________________________________________________________
AliITSVertexer::AliITSVertexer(const AliITSVertexer &vtxr) : AliVertexer(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexer","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSVertexer& AliITSVertexer::operator=(const AliITSVertexer& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
AliITSVertexer::~AliITSVertexer() {
  // Destructor
 if(fLadders) delete [] fLadders;
}

//______________________________________________________________________
void AliITSVertexer::FindMultiplicity(Int_t evnumber){
  // Invokes AliITSMultReconstructor to determine the
  // charged multiplicity in the pixel layers
  if(fMult){delete fMult; fMult = 0;}
  Bool_t success=kTRUE;
  if(!fCurrentVertex)success=kFALSE;
  if(fCurrentVertex && fCurrentVertex->GetNContributors()<1)success=kFALSE;
  if(!success){
    AliWarning("Tracklets multiplicity not determined because the primary vertex was not found");
    return;
  }
  AliITSMultReconstructor* multReco = new AliITSMultReconstructor();
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  multReco->SetGeometry(itsLoader->GetITSgeom());
  itsLoader->LoadRecPoints();
  rl->GetEvent(evnumber);
  TTree* itsClusterTree = itsLoader->TreeR();
  if (!itsClusterTree) {
    AliError(" Can't get the ITS cluster tree !\n");
    return;
  }
  Double_t vtx[3];
  fCurrentVertex->GetXYZ(vtx);
  Float_t vtxf[3];
  for(Int_t i=0;i<3;i++)vtxf[i]=vtx[i];
  multReco->SetHistOn(kFALSE);
  multReco->Reconstruct(itsClusterTree,vtxf,vtxf);
  Int_t notracks=multReco->GetNTracklets();
  Float_t *tht = new Float_t [notracks];
  Float_t *phi = new Float_t [notracks];
  Float_t *dphi = new Float_t [notracks];
  Int_t *labels = new Int_t[notracks];
  for(Int_t i=0;i<multReco->GetNTracklets();i++){
    tht[i] = multReco->GetTracklet(i)[0];
    phi[i] =  multReco->GetTracklet(i)[1];
    dphi[i] = multReco->GetTracklet(i)[2];
    labels[i] = static_cast<Int_t>(multReco->GetTracklet(i)[3]);
  }
  Int_t nosingleclus=multReco->GetNSingleClusters();
  Float_t *ths = new Float_t [nosingleclus];
  Float_t *phs = new Float_t [nosingleclus];
  for(Int_t i=0;i<nosingleclus;i++){
    ths[i] = multReco->GetCluster(i)[0];
    phs[i] =  multReco->GetCluster(i)[1];
  }
  fMult = new AliMultiplicity(notracks,tht,phi,dphi,labels,nosingleclus,ths,phs);
  delete [] tht;
  delete [] phi;
  delete [] dphi;
  delete [] ths;
  delete [] phs;
  delete [] labels;
  itsLoader->UnloadRecPoints();
  delete multReco;
  return;
}

//______________________________________________________________________
void AliITSVertexer::SetLaddersOnLayer2(Int_t ladwid){
  // Calculates the array of ladders on layer 2 to be used with a 
  // given ladder on layer 1
  fLadOnLay2=ladwid;
  //  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  //  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  //  AliITSgeom* geom = itsLoader->GetITSgeom();
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
void AliITSVertexer::WriteCurrentVertex(){
  // Write the current AliVertex object to file fOutFile
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  fCurrentVertex->SetName("Vertex");
  //  const char * name = fCurrentVertex->GetName();
  //  itsLoader->SetVerticesContName(name);
  Int_t rc = itsLoader->PostVertex(fCurrentVertex);
  rc = itsLoader->WriteVertices();
}

