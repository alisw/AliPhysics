#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include "AliGeomManager.h"
#include "AliITSDetTypeRec.h"
#include "AliITSInitGeometry.h"
#include "AliITSMeanVertexer.h"
#include "AliITSRecPointContainer.h"
#include "AliITSLoader.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliRunLoader.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexer3DTapan.h"
#include "AliESDVertex.h"
#include "AliMeanVertex.h"
//#include "AliCodeTimer.h"

const Int_t AliITSMeanVertexer::fgkMaxNumOfEvents = 10000;
ClassImp(AliITSMeanVertexer)

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Class to compute vertex position using SPD local reconstruction   //
// An average vertex position using all the events                   //
// is built and saved                                                //
// Individual vertices can be optionally stored                      //
// Origin: M.Masera  (masera@to.infn.it)                             //
// Usage:                                                            //
// Class used by the ITSSPDVertexDiamondda.cxx detector algorithm    //
///////////////////////////////////////////////////////////////////////

/* $Id$ */

//______________________________________________________________________
AliITSMeanVertexer::AliITSMeanVertexer(Bool_t mode):TObject(),
fDetTypeRec(NULL),
fVertexXY(NULL),
fVertexZ(NULL),
fNoEventsContr(0),
fTotContributors(0.),
fAverContributors(0.), 
fFilterOnContributors(0),
fMode(mode),
fVertexer(NULL),
fAccEvents(fgkMaxNumOfEvents), 
fVertArray("AliESDVertex",fgkMaxNumOfEvents),
fClu0(NULL),
fIndex(0),
fErrXCut(0.),
fRCut(0.),
fLowSPD0(0),
fHighSPD0(0),
fMultH(NULL),
fErrXH(NULL),
fMultHa(NULL),
fErrXHa(NULL),
fDistH(NULL)
{
  // Default Constructor
  Reset();

  // Histograms initialization
  const Float_t xLimit = 5.0, yLimit = 5.0, zLimit = 50.0;
  const Float_t xDelta = 0.02, yDelta = 0.02, zDelta = 0.2;
  fVertexXY = new TH2F("VertexXY","Vertex Diamond (Y vs X)",
		       2*(Int_t)(xLimit/xDelta),-xLimit,xLimit,
		       2*(Int_t)(yLimit/yDelta),-yLimit,yLimit);
  fVertexXY->SetXTitle("X , cm");
  fVertexXY->SetYTitle("Y , cm");
  fVertexXY->SetOption("colz");
  fVertexZ  = new TH1F("VertexZ"," Longitudinal Vertex Profile",
		       2*(Int_t)(zLimit/zDelta),-zLimit,zLimit);
  fVertexZ->SetXTitle("Z , cm");
  fErrXH = new TH1F("errorX","Error X - before cuts",100,0.,99.);
  fMultH = new TH1F("multh","mult on layer 1 before cuts",2400,0.,7200.);
  fErrXHa = new TH1F("errorXa","Error X - after Filter",100,0.,99.);
  fErrXHa->SetLineColor(kRed);
  fMultHa = new TH1F("multha","mult on layer 1 - after Filter",2400,0.,7200.);
  fMultHa->SetLineColor(kRed);
  fDistH = new TH1F("disth","distance (mm)",100,0.,30.);
  fClu0 = new UInt_t [fgkMaxNumOfEvents];
  for(Int_t i=0;i<fgkMaxNumOfEvents;i++)fClu0[i]=0;
  fAccEvents.ResetAllBits(kTRUE);
}

//______________________________________________________________________
void AliITSMeanVertexer::Reset(){
  // Reset averages 
  for(Int_t i=0;i<3;i++){
    fWeighPosSum[i] = 0.;
    fWeighSigSum[i] = 0.;
    fAverPosSum[i] = 0.;
    fWeighPos[i] = 0.;
    fWeighSig[i] = 0.;
    fAverPos[i] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSq[i][j] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSqSum[i][j] = 0.;
    fNoEventsContr=0;
    fTotContributors = 0.;
    if(fVertexXY)fVertexXY->Reset();
    if(fVertexZ)fVertexZ->Reset();
  }
}
//______________________________________________________________________
Bool_t AliITSMeanVertexer::Init() {
  // Initialize filters
  // Initialize geometry
  // Initialize ITS classes
 
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::ApplyAlignObjsFromCDB("ITS")) return kFALSE;

  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  if (!geom) return kFALSE;
  printf("Geometry name: %s \n",(initgeom.GetGeometryName()).Data());

  fDetTypeRec = new AliITSDetTypeRec();
  fDetTypeRec->SetLoadOnlySPDCalib(kTRUE);
  fDetTypeRec->SetITSgeom(geom);
  fDetTypeRec->SetDefaults();
  fDetTypeRec->SetDefaultClusterFindersV2(kTRUE);

  // Initialize filter values to their defaults
  SetFilterOnContributors();
  SetCutOnErrX();
  SetCutOnR();
  SetCutOnCls();

  // Instatiate vertexer
  if (!fMode) {
    fVertexer = new AliITSVertexer3DTapan(1000);
  }
  else {
    fVertexer = new AliITSVertexer3D();
    fVertexer->SetDetTypeRec(fDetTypeRec);
    AliITSVertexer3D* alias = (AliITSVertexer3D*)fVertexer;
    alias->SetWideFiducialRegion(40.,0.5);
    alias->SetNarrowFiducialRegion(0.5,0.5);
    alias->SetDeltaPhiCuts(0.5,0.025);
    alias->SetDCACut(0.1);
    alias->SetPileupAlgo(3);
    fVertexer->SetComputeMultiplicity(kFALSE);
  }
  return kTRUE;
}

//______________________________________________________________________
AliITSMeanVertexer::~AliITSMeanVertexer() {
  // Destructor
  delete fDetTypeRec;
  delete fVertexXY;
  delete fVertexZ;
  delete fMultH;
  delete fErrXH;
  delete fMultHa;
  delete fErrXHa;
  delete fDistH;
  delete fVertexer;
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Reconstruct(AliRawReader *rawReader){
  // Performs SPD local reconstruction
  // and vertex finding
  // returns true in case a vertex is found

  // Run SPD cluster finder
  static Int_t evcount = -1;
  if(evcount <0){
    evcount++;
    AliInfo(Form("Low and high cuts on SPD L0 clusters %d , %d \n",fLowSPD0,fHighSPD0));
    AliInfo(Form("Reconstruct: cut on errX %f \n",fErrXCut));
  }
//  AliCodeTimerAuto("",0);
  AliITSRecPointContainer::Instance()->PrepareToRead();
  TTree* clustersTree = new TTree("TreeR", "Reconstructed Points Container"); //make a tree
  fDetTypeRec->DigitsToRecPoints(rawReader,clustersTree,"SPD");

  Bool_t vtxOK = kFALSE;
  UInt_t nl1=0;
  AliESDVertex *vtx = fVertexer->FindVertexForCurrentEvent(clustersTree);
  if (!fMode) {
    if (TMath::Abs(vtx->GetChi2()) < 0.1) vtxOK = kTRUE;
  }
  else {
    AliITSRecPointContainer* rpcont= AliITSRecPointContainer::Instance();
    nl1=rpcont->GetNClustersInLayerFast(1);
    if(Filter(vtx,nl1)) vtxOK = kTRUE;
    /*    if(vtx){
      if(vtxOK){
	printf("The vertex is OK\n");
      }
      else {
	printf("The vertex is NOT OK\n");
      }
      vtx->PrintStatus();
    }
    else {
      printf("The vertex was not reconstructed\n");
      }  */
  }
  delete clustersTree;
  if (vtxOK && fMode){
    new(fVertArray[fIndex])AliESDVertex(*vtx);
    fClu0[fIndex]=nl1;
    fAccEvents.SetBitNumber(fIndex);
    fIndex++;
    if(fIndex>=fgkMaxNumOfEvents){
      if(ComputeMean(kFALSE)){
	if(ComputeMean(kTRUE))AliInfo("Mean vertex computed");
      }
    }
  }
  if (vtx) delete vtx;

  return vtxOK;
}

//______________________________________________________________________
void AliITSMeanVertexer::WriteVertices(const char *filename){
  // Compute mean vertex and
  // store it along with the histograms
  // in a file
  
  TFile fmv(filename,"update");
  Bool_t itisOK = kFALSE;
  if(ComputeMean(kFALSE)){
    if(ComputeMean(kTRUE)){
      Double_t cov[6];
      cov[0] =  fAverPosSq[0][0];  // variance x
      cov[1] =  fAverPosSq[0][1];  // cov xy
      cov[2] =  fAverPosSq[1][1];  // variance y
      cov[3] =  fAverPosSq[0][2];  // cov xz
      cov[4] =  fAverPosSq[1][2];  // cov yz
      cov[5] =  fAverPosSq[2][2];  // variance z
      AliMeanVertex mv(fWeighPos,fWeighSig,cov,fNoEventsContr,0,0.,0.);
      mv.SetTitle("Mean Vertex");
      mv.SetName("MeanVertex");
      // we have to add chi2 here
      AliESDVertex vtx(fWeighPos,cov,0,TMath::Nint(fAverContributors),"MeanVertexPos");

      mv.Write(mv.GetName(),TObject::kOverwrite);
      vtx.Write(vtx.GetName(),TObject::kOverwrite);
      itisOK = kTRUE;
    }
  }

  if(!itisOK){
    AliError(Form("Evaluation of mean vertex not possible. Number of used events = %d",fNoEventsContr));
  }

  fVertexXY->Write(fVertexXY->GetName(),TObject::kOverwrite);
  fVertexZ->Write(fVertexZ->GetName(),TObject::kOverwrite);
  fMultH->Write(fMultH->GetName(),TObject::kOverwrite);
  fErrXH->Write(fErrXH->GetName(),TObject::kOverwrite);
  fMultHa->Write(fMultHa->GetName(),TObject::kOverwrite);
  fErrXHa->Write(fErrXHa->GetName(),TObject::kOverwrite);
  fDistH->Write(fDistH->GetName(),TObject::kOverwrite);
  fmv.Close();
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Filter(AliESDVertex *vert,UInt_t mult){
  // Apply selection criteria to events
  Bool_t status = kFALSE;
  if(!vert)return status;
  // Remove vertices reconstructed with vertexerZ
  if(strcmp(vert->GetName(),"SPDVertexZ") == 0) return status;
  Int_t ncontr = vert->GetNContributors();
  AliDebug(1,Form("Number of contributors = %d",ncontr));
  Double_t ex = vert->GetXRes();
  fMultH->Fill(mult);
  fErrXH->Fill(ex*1000.);
  if(ncontr>fFilterOnContributors && (mult>fLowSPD0 && mult<fHighSPD0) && 
     ((ex*1000.)<fErrXCut)) {
    status = kTRUE;
    fMultHa->Fill(mult);
    fErrXHa->Fill(ex*1000.);
  }
  return status;
}

//______________________________________________________________________
void AliITSMeanVertexer::AddToMean(AliESDVertex *vert){
  // update mean vertex
  Double_t currentPos[3],currentSigma[3];
  vert->GetXYZ(currentPos);
  vert->GetSigmaXYZ(currentSigma);
  Bool_t goon = kTRUE;
  for(Int_t i=0;i<3;i++)if(currentSigma[i] == 0.)goon = kFALSE;
  if(!goon)return;
  for(Int_t i=0;i<3;i++){
    fWeighPosSum[i]+=currentPos[i]/currentSigma[i]/currentSigma[i];
    fWeighSigSum[i]+=1./currentSigma[i]/currentSigma[i];
    fAverPosSum[i]+=currentPos[i];
  }
  for(Int_t i=0;i<3;i++){
    for(Int_t j=i;j<3;j++){
      fAverPosSqSum[i][j] += currentPos[i] * currentPos[j];
    }
  }

  fTotContributors+=vert->GetNContributors();
  fVertexXY->Fill(currentPos[0],currentPos[1]);
  fVertexZ->Fill(currentPos[2]);

  fNoEventsContr++;
}

//______________________________________________________________________
 Bool_t AliITSMeanVertexer::ComputeMean(Bool_t killOutliers){
   // compute mean vertex
   // called once with killOutliers = kFALSE and then again with
   // killOutliers = kTRUE

   static Bool_t firstime = kTRUE;
   if(fNoEventsContr>2 && !killOutliers)return kTRUE;  //ComputeMean is
                             // executed only once with argument kFALSE
   Double_t wpos[3];
   for(Int_t i=0;i<3;i++)wpos[i]=fWeighPos[i];

   if(killOutliers && firstime){
     firstime = kFALSE;
     Reset();
   }

   for(Int_t i=0;i<fVertArray.GetEntries();i++){
     AliESDVertex* vert = NULL;

     // second pass
     if(killOutliers && fAccEvents.TestBitNumber(i)){
       vert=(AliESDVertex*)fVertArray[i];
       Double_t dist=(vert->GetXv()-wpos[0])*(vert->GetXv()-wpos[0]);
       dist+=(vert->GetYv()-wpos[1])*(vert->GetYv()-wpos[1]);
       dist=sqrt(dist)*10.;    // distance in mm
       fDistH->Fill(dist);
       if(dist>fRCut)fAccEvents.SetBitNumber(i,kFALSE);
     }
     
     if(!fAccEvents.TestBitNumber(i))continue;
     vert=(AliESDVertex*)fVertArray[i];
     AddToMean(vert);
   }
   if(fNoEventsContr < 2) return kFALSE;
   Double_t nevents = fNoEventsContr;
   for(Int_t i=0;i<3;i++){
     fWeighPos[i] = fWeighPosSum[i]/fWeighSigSum[i]; 
     fWeighSig[i] = 1./TMath::Sqrt(fWeighSigSum[i]);
     fAverPos[i] = fAverPosSum[i]/nevents;
   } 
   for(Int_t i=0;i<3;i++){
     for(Int_t j=i;j<3;j++){
       fAverPosSq[i][j] = fAverPosSqSum[i][j]/(nevents -1.);
       fAverPosSq[i][j] -= nevents/(nevents -1.)*fAverPos[i]*fAverPos[j];
     }
   } 
   if(killOutliers)ResetArray();
   fAverContributors = fTotContributors/nevents;
   return kTRUE;
 }
