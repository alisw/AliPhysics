/*
 .L $AliPhysics_SRC/PWGPP/AliESDtoolTRDTOF.cxx++
 */

#include "AliESDtoolTRDTOF.h"
#include "AliESDtrack.h"
#include "TTree.h"
#include "TH1F.h"


AliESDtoolTRDTOF::AliESDtoolTRDTOF():
  fESDtrack(0),
  fTree(0),
  fVecSec(7),
  fVecZ(7), fVecDet(7), fVecdSec(7), fVecdEdge(7), fVecActive(7), fVecStatus(7),fVecDeadZ(7),fVecDeadR(7), fVecDeadDet0(7),fVecDeadDet(7)
{
   CacheTRDGeom();
}

/// \brief cache TRD geometry
void AliESDtoolTRDTOF::CacheTRDGeom(){
  for (Int_t iLayer=0; iLayer<6; iLayer++){
    for (Int_t iStack=0;iStack<5; iStack++){
      fZBoundary0(iLayer,iStack)=geom.GetRow0(iLayer,iStack,0);
      fZBoundary1(iLayer,iStack)=-geom.GetRow0(iLayer,4-iStack,0);
    }
  }
  for (Int_t iLayer=0; iLayer<6; iLayer++) fZBoundary1(iLayer,5)=-fZBoundary0(iLayer,0);
}

/// \brief build active detector map usinng skimmed V0 data

void AliESDtoolTRDTOF::MakeActiveMapFromV0(TTree *treeV0, Int_t nPoints, Double_t threshold) {0
  fHistoDetector=new TH1F("fHistoDetector","fHistoDetector",kNDetectors,0,kNDetectors);
  for (Int_t iLayer=0; iLayer<6; iLayer++){
    TString what=TString::Format("GetDet(%d)>>+fHistoDetector",iLayer);
    TString where=TString::Format("isOK0&&track0.GetTRDtrkltOccupancy(%d)>10",iLayer);
    treeV0->Draw(what.Data(),where.Data(),"goff",nPoints);
  }
  treeV0->Draw("GetDet(6)>>+fHistoDetector","isOK0&&isTOFIn0","goff",nPoints);
  Double_t entriesMedian=TMath::Median(10,fHistoDetector->GetArray());
  for (Int_t iDet=0;iDet<kNDetectors; iDet++){
    fActiveMap(iDet) = (fHistoDetector->GetBinContent(iDet+1)>threshold*entriesMedian);
  }
}


/// \brief cache per track information
/// status code:
///        0 - tracklet  found
///        1 - tracklet in active zone but not found
///        2 - track in dead zone
///        3 - track in not active detector
Double_t AliESDtoolTRDTOF::ChacheTrackInfo(){
  //treeV0->GetEntry(entry);
  const Double_t tanPhi=TMath::Tan(TMath::Pi()/18.);
  const Double_t edgeF0=TMath::Sqrt(1+tanPhi*tanPhi); // sqrt(1+TMath::Tan(TMath::Pi()/18.)**2)
  nFindable=0;
  nFound=0;
  for (Int_t iLayer=0; iLayer<7; iLayer++){
    Double_t radius=geom.GetTime0(iLayer%6)-3;  // to add constant
    if (iLayer==6) radius=375;
    Double_t phi=fESDtrack->GetInnerParam()->GetParameterAtRadius(radius*edgeF0,5,7);
    Double_t gZ=fESDtrack->GetInnerParam()->GetParameterAtRadius(radius*edgeF0,5,2);
    if (phi<0) phi+=TMath::TwoPi();
    Double_t sector=9.*phi/TMath::Pi();
    Int_t iStack=0;
    for (iStack=0; iStack<5; iStack++) {
      if (gZ > fZBoundary1(iLayer, iStack)) break;
    }

    if (iStack>4) iStack=4;
    Int_t iDet=0;
    if (iLayer<6)  iDet=geom.GetDetector(iLayer,iStack,sector);    // TRD numbering
    if (iLayer==6) iDet=540+sector*5+iStack;                       // TOF numbering
    fVecDet[iLayer]=iDet;
    fVecSec[iLayer]=sector;
    fVecZ[iLayer]=gZ;
    Double_t dSec=sector-int(sector);
    Double_t dEdge=(dSec<0.5) ? dSec*radius*tanPhi:(1-dSec)*radius*tanPhi;
    fVecdSec[iLayer]=dSec;
    fVecdEdge[iLayer]=dEdge;
    fVecStatus[iLayer]=0;
    fVecActive[iLayer]=0;
    //
    Bool_t isActive=0;
    Bool_t isFindable = kFALSE;
    Bool_t isActiveZ=(gZ<fZBoundary0(iLayer,iStack) && gZ>fZBoundary1(iLayer, iStack));
    Bool_t isDeadR=TMath::Abs(dEdge)<kMarginR;
    Bool_t isDeadZ=(!(gZ<fZBoundary0(iLayer,iStack)-kMarginZ && gZ>fZBoundary1(iLayer, iStack)+kMarginZ));   /// add savety margin 1 cm because of misalignemnt for TOF not sure of z bounaries
    Bool_t isDeadDet=fActiveMap(iDet)<0.5;
    Int_t status=0;
    if (iLayer<6) {
      fVecActive[iLayer]=(fESDtrack->GetTRDtrkltOccupancy(iLayer)>occuCut);
      if ( fVecActive[iLayer] > 0) {
        nFindable++;
        nFound++;
        isFindable=kTRUE;
        status=0;
        isActive=kTRUE;
      } else {
        isFindable=! (isDeadR || isDeadZ || isDeadDet);
        if (nFindable) nFindable++;
      }
    }else{
      fVecActive[iLayer]=fESDtrack->IsOn(0x2000);
      isActive=fESDtrack->IsOn(0x2000);
      isFindable=! (isDeadR || isDeadZ || isDeadDet);
    }
    fVecDeadZ[iLayer]=isDeadZ;
    fVecDeadR[iLayer]=isDeadR;
    fVecDeadDet[iLayer]=isDeadDet;
    if (fVecActive[iLayer]>0.5) continue;
    if (isFindable){
      fVecStatus[iLayer]=1;   // not found  but active
      continue;
    }
    if (isDeadDet) {fVecStatus[iLayer]=3; continue;}   // not found  but active
    fVecStatus[iLayer]=2; // dead zones
  }
  return 1;
}

/// initialize input V0 tree
/// \param tree  - input tree with V0s
void AliESDtoolTRDTOF::InitTreeV0(TTree * tree) {
  if (tree== nullptr) {
    fTree = AliXRDPROOFtoolkit::MakeChain("filtered.list", "V0s", 0, 1000);
  }else{
    fTree =tree;
  }
  
  fTree->SetAlias("isITSOn0","(track0.fFlags&0x1)>0");
  fTree->SetAlias("isITSOn1","(track1.fFlags&0x1)>0");
  fTree->SetAlias("idProton0", "type==4&&abs(tpcNsigma0.fElements[4])<2"); // stong criteria to suppress backround
  fTree->SetAlias("idProton1", "type==2&&abs(tpcNsigma1.fElements[4])<2"); // stong criteria to suppress backround
  fTree->SetAlias("isTOFIn0","abs(tofClInfo0.fElements[3])<20");
  fTree->SetAlias("isTOFGold0","abs(tofClInfo0.fElements[3])<2");
  fTree->SetAlias("isTOFProton0","abs(tofNsigma0.fElements[4])<6");
  fTree->SetAlias("isTOFIn1","abs(tofClInfo1.fElements[3])<20");
  fTree->SetAlias("isTOFProton1","abs(tofNsigma1.fElements[4])<6");
  fTree->SetAlias("dTRD0","(9*track0.fIp.GetParameterAtRadius(305,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(305,5,7)/pi+18)");
  fTree->SetAlias("dTRD1","(9*track0.fIp.GetParameterAtRadius(318,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(318,5,7)/pi+18)");
  //fTree->Draw("(9*track0.fIp.GetParameterAtRadius(300,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(300,5,7)/pi+18)","idProton0&&track0.Pt()>0.5","",100000);
  fTree->SetBranchAddress("track0.",&track);
  fTree->SetAlias("isOK0","abs(track0.fP[4])<2.5&&abs(track0.fP[3])<1&&chacheInfo(Entry$)");
  //treeV0->SetAlias("notActive")
}
