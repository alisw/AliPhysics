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

/* mailto: svallero@to.infn.it, s.lapointe@cern.ch */

/* Class defining containers for the HF b-jets analysis */

#include "AliHFJetsTaggingVertex.h"
#include "AliHFJetsContainerVertex.h"
#include "AliHFJetsContainer.h"
#include "AliEmcalJet.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TList.h"
#include "TCollection.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliLog.h"
#include "TCanvas.h"
#include "TF1.h"
#include "AliTHn.h"
#include "THn.h"
#include "TClonesArray.h"

ClassImp(AliHFJetsContainerVertex)

AliHFJetsContainerVertex::AliHFJetsContainerVertex(): 
  AliHFJetsContainer("",kTRUE),
  fType(kJets)
{
  fTagger=new AliHFJetsTaggingVertex("tagger");
  // dummy
}

AliHFJetsContainerVertex::AliHFJetsContainerVertex(const char* name,ContType contType): 
  AliHFJetsContainer(name,kTRUE),
  fType(kJets),
  fTagger(0x0)
{

  // Constructor
  fTagger=new AliHFJetsTaggingVertex("tagger");
  CreateContainerVertex(contType); 
}

//----------------------------------------------------------------
AliHFJetsContainerVertex::AliHFJetsContainerVertex(const AliHFJetsContainerVertex &c) :
  AliHFJetsContainer("",kTRUE)
{

  // AliHFJetsContainerVertex copy constructor

  ((AliHFJetsContainerVertex &) c).Copy(*this);
}

//----------------------------------------------------------------
AliHFJetsContainerVertex &AliHFJetsContainerVertex::operator=(const AliHFJetsContainerVertex &c)
{
  // assigment operator

  if (this != &c)
    ((AliHFJetsContainerVertex &) c).Copy(*this);

  return *this;
}
//----------------------------------------------------------------
AliHFJetsContainerVertex::~AliHFJetsContainerVertex()
{

  // Destructor
  if (fTagger){
    delete fTagger;
    fTagger=0;
    } 
}

//----------------------------------------------------------------
void AliHFJetsContainerVertex::Copy(TObject& c) const
{
  // copy function

  AliHFJetsContainerVertex& target = (AliHFJetsContainerVertex &) c;
  if (fType)
    target.fType = fType;
   
}
//----------------------------------------------------------------
void AliHFJetsContainerVertex::CreateContainerVertex(ContType contType)
{
  TString stype, vars;
    fType=contType;
  switch (contType) {
    case 0:
      stype.Form("Jets properties.");
      // AliInfo(Form("Creating AliHFJetsContainerVertex of type: %s ", stype.Data()));
      // Relevant variables
      vars.Form("nTrk;nEle;partDP;partBP;partBH;partContr;ptDP;ptBP;ptBH;areaCh");
      break;
    case 1:
      stype.Form("B-jets properties.");
      //AliInfo(Form("Creating AliHFJetsContainerVertex of type: %s ", stype.Data()));
      // Relevant variables
      vars.Form("nTrk;nRecoVtx;nEle;partDP;partContr;ptDP;areaCh;partBP;partBH");
      break;
    case 2:
      stype.Form("Primary vertexes properties.");
      //AliInfo(Form("Creating AliHFJetsContainerVertex of type: %s ", stype.Data()));
      // Relevant variables
      vars.Form("nTrk;ptVtx;mass;Lxy;chi2/ndf;partDP;nRealVtx;deltaX;deltaY;sigVtx;LxyJet;LxyJetPerp;cosTheta;nFromB;nFromBD;partBP;partBH");
      break;
    case 3:
      stype.Form("Jet vertexes proprties.");
      //AliInfo(Form("Creating AliHFJetsContainerVertex of type: %s ", stype.Data()));
      // Relevant variables
      // vars.Form("nRecoVtx;SumMass3MostDispl;Lxy1;Lxy2;Lxy3;mass1;mass2;mass3;nReal1;nReal2;nReal3;nFromB1;nFromB2;nFromB3;nFromPrD1;nFromPrD2;nFromPrD3;partDP;partBP;partBH;sigvtx1;sigvtx2;sigvtx3;sos1;sos2;sos3;eLxy1;sigsos1");
vars.Form("nRecoVtx;SumMass3MostDispl;Lxy1;mass1;nReal1;nFromB1;nFromPrD1;partDP;partBP;partBH;sigvtx1;sos1;eLxy1;sigsos1");
      break;
    case 4:
      stype.Form("Jet vertexes proprties data.");
      //AliInfo(Form("Creating AliHFJetsContainerVertex of type: %s ", stype.Data()));
      // Relevant variables
      // vars.Form("nRecoVtx;SumMass3MostDispl;Lxy1;Lxy2;Lxy3;mass1;mass2;mass3;nReal1;nReal2;nReal3;nFromB1;nFromB2;nFromB3;nFromPrD1;nFromPrD2;nFromPrD3;partDP;partBP;partBH;sigvtx1;sigvtx2;sigvtx3;sos1;sos2;sos3;eLxy1;sigsos1");
vars.Form("nRecoVtx;Lxy1;mass1;sigvtx1;eLxy1");
      break;
    default:
      AliError("Not a valid container type!");
      return;
  }

  AliInfo(Form(mage"Container type set to %s: %s" Bee, strContType(contType), stype.Data())); 
  // Get binning for each variable
  TObjArray *arr;
  //arr->SetOwner(kTRUE);
  TObjString *objstr;
  arr = vars.Tokenize(";");
  //static const Int_t nvars = arr->GetEntriesFast();
  Int_t nvars = arr->GetEntriesFast();
  Int_t nbins[nvars];       // number of bins for each variable
  const char* axistitle[nvars]; // axis title for each variable
  const char* axisname[nvars]; // axis name for each variable
  Double_t *binning[nvars]; // array of bins for each variable
  TIter next(arr);
  Int_t i = 0;
  while ((objstr=(TObjString*)next())){
    binning[i]=new Double_t[1000]; 
    GetBinningVertex(objstr->GetString(), nbins[i], binning[i], axistitle[i]);
    axisname[i]=objstr->GetName();
    i++;
    }
  
  AliHFJetsContainer::CreateCustomContainer(nvars,axisname,nbins, binning,axistitle);

  // Delete arrays 
  delete arr;
  for (Int_t k=0; k<nvars; k++){
    delete [] binning[k];
    //delete axisname[k];
    //delete axistitle[k];
 }
}

//----------------------------------------------------------------
void AliHFJetsContainerVertex::GetBinningVertex(TString var, Int_t& nBins, Double_t *bins,const char*& axistitle)
{

  // Assigns variable-specific bin edges 
  // (you can define array of bins "by hand" to allow for non uniform bin width) 
  //if (var.Contains("DP") || var.Contains("BP") || var.Contains("BH") ){
  //  if (var.Contains("part")) var.Form("idPart");
  //  else if (var.Contains("pt")) var.Form("ptPart");
  //}
 
  Float_t binmin=0., binmax=0.;    
  if (var.EqualTo("jetPt")){
      axistitle="p_{T,jet} (GeV/c)";
      nBins = 100; binmin= 0.; binmax= 100.;
      // Double_t *bins = {5.,10.,15., ...};
      // return bins;
  } else if (var.EqualTo("jetEta")){
      axistitle="#eta_{jet}";
      nBins = 20; binmin= -1.; binmax= 1.;
  } else if (var.EqualTo("jetPhi")){
      axistitle="#phi_{jet} (rad)";
      nBins = 20; binmin= -TMath::Pi(); binmax= TMath::Pi();
  } else if (var.EqualTo("nTrk")){
      axistitle="N_{trk}";
      nBins = 20; binmin= 0.99; binmax= 20.99;
  } else if (var.EqualTo("nEle")){
      axistitle="N_{ele}";
      nBins = 5; binmin= 0.; binmax= 4.99;
  } else if (var.EqualTo("partDP")){
      axistitle="ID_{parton} DP";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("partBP")){
      axistitle="ID_{parton} BP";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("partBH")){
      axistitle="ID_{parton} BH";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("partContr")){
      axistitle="Contribution";
      nBins = 10; binmin= 0.; binmax= 2.;
  } else if (var.EqualTo("ptDP")){
      axistitle="p_{T,part} DP (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  } else if (var.EqualTo("ptBP")){
      axistitle="p_{T,part} BP (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  } else if (var.EqualTo("ptBH")){
      axistitle="p_{T,part} BH (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  } else if (var.EqualTo("areaCh")){
      axistitle="Charged area";
      nBins = 50; binmin= 0.; binmax= 1.;
  } else if (var.EqualTo("ptVtx")){
      axistitle="p_{T,vtx} (GeV/c)";
      nBins = 20; binmin= 0.; binmax= 20.;
  } else if (var.EqualTo("mass")){
      axistitle="mass (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  } else if (var.EqualTo("mass1")){
      axistitle="mass_{1} (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  } else if (var.EqualTo("mass2")){
      axistitle="mass_{2} (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  } else if (var.EqualTo("mass3")){
      axistitle="mass_{3} (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  } else if (var.EqualTo("SumMass3MostDispl")){
      axistitle="#Sigma mass (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 20.;
  } else if (var.EqualTo("Lxy")){
      axistitle="L_{xy} (cm)";
      nBins = 200; binmin= -1.; binmax= 1.;
  } else if (var.EqualTo("Lxy1")){
      axistitle="L_{xy,1} (cm)";
      nBins = 50; binmin= 0.; binmax= 1.;
  } else if (var.EqualTo("Lxy2")){
      axistitle="L_{xy,2} (cm)";
      nBins = 50; binmin= 0.; binmax= 0.5;
  } else if (var.EqualTo("Lxy3")){
      axistitle="L_{xy,3} (cm)";
      nBins = 50; binmin= 0.; binmax= 0.5;
  } else if (var.EqualTo("chi2/ndf")){
      axistitle="#chi^{2}/NDF";
      nBins = 10; binmin= 0.; binmax= 10.;
  } else if (var.EqualTo("nRealVtx")){
      axistitle="N_{vtx,real}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nRecoVtx")){
      axistitle="N_{vtx,reco}";
      nBins = 21; binmin= -0.5; binmax= 20.5; 
  } else if (var.EqualTo("deltaX")){
      axistitle="#Delta x (cm)";
      nBins = 15; binmin= -0.03; binmax= 0.03;
  } else if (var.EqualTo("deltaY")){
      axistitle="#Delta y (cm)";
      nBins = 15; binmin= -0.03; binmax= 0.03;
  } else if (var.EqualTo("sigVtx")){
      axistitle="#sigma_{vtx}";
      nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("LxyJet")){
      axistitle="#L_{xy,jet} (cm)";
      nBins = 100; binmin= -1.; binmax= 1.;
  } else if (var.EqualTo("LxyJetPerp")){
      axistitle="#L_{xy,jet} #perp (cm)";
      nBins = 30; binmin= 0.; binmax= 0.2;
  } else if (var.EqualTo("cosTheta")){
      axistitle="cos#theta";
      nBins = 50; binmin= -1.; binmax= 1.;
  } else if (var.EqualTo("nFromB")){
      axistitle="N (from B)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromB1")){
      axistitle="N (from B1)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromB2")){
      axistitle="N (from B2)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromB3")){
      axistitle="N (from B3)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromBD")){
      axistitle="N (from D from B)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromPrD1")){
      axistitle="N (from prompt D1)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromPrD2")){
      axistitle="N (from prompt D2)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nFromPrD3")){
      axistitle="N (from prompt D3)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nReal1")){
      axistitle="N_{real,1}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nReal2")){
      axistitle="N_{real,2}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("nReal3")){
      axistitle="N_{real,3}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  } else if (var.EqualTo("sigvtx1")){
    axistitle="#sigma_{vertex,1}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("sigvtx2")){
    axistitle="#sigma_{vertex,2}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("sigvtx3")){
    axistitle="#sigma_{vertex,3}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  }  else if (var.EqualTo("sos1")){
    axistitle="sum of d_{0}^{2} (cm)";
    nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("sos2")){
    axistitle="sum of d_{0}^{2} (cm)";
    nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("sos3")){
    axistitle="sum of d_{0}^{2} (cm)";
    nBins = 20; binmin= 0.; binmax= 0.1;
  } else if (var.EqualTo("eLxy1")){
    axistitle="L_{xy}/#sigma_{L_{xy}}";
    nBins = 100; binmin= 0.; binmax= 100.;
  }else if (var.EqualTo("sigsos1")){
    axistitle="sum of sigd_{0}^{2}";
    nBins = 60; binmin= 0.; binmax= 30.;
  }else {
      AliError(Form("Variable %s not defined!", var.Data()));
      }
  

  // Define regular binning
  Double_t binwidth = (binmax-binmin)/(1.*nBins);
  for (Int_t j=0; j<nBins+1; j++){
     if (j==0) bins[j]= binmin;
     else bins[j] = bins[j-1]+binwidth;
     }
  //return bins;
}

void AliHFJetsContainerVertex::FillStepJets(AliHFJetsContainer::CFSteps step,Double_t mult, const AliEmcalJet *jet,Double_t p[3],Double_t contr,Double_t pt[3]){
  
  if (fType != kJets){
    AliError(Form(RED"This method is available only for container type kJets: you are trying to fill %s!" Bee, strContType(fType)));
  }

  // here the number of electrons is not filled (after nTrk)
  // Double_t point[14]={mult*1.,jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),jet->GetNumberOfTracks()*1.,0.,p[0],p[1],p[2],contr,pt[0],pt[1],pt[2],jet->EffectiveAreaCharged()};
  Double_t point[14]={mult*1.,jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),jet->GetNumberOfTracks()*1.,0.,p[0],p[1],p[2],contr,pt[0],pt[1],pt[2],jet->Area()}; // no eff area ch. for emcaljet, should check that area gives the same thing
  TArrayD *apoint=new TArrayD(14,point);
  AliHFJetsContainer::FillStep(step, apoint);
}

void AliHFJetsContainerVertex::FillStepQaVtx(AliHFJetsContainer::CFSteps step, Double_t mult, const AliEmcalJet *jet,const TClonesArray *vertices,Double_t* disp,Int_t nvtx,const AliAODVertex *primVtx,const TClonesArray *mcPart,Double_t p[3]){
  
  if (fType != kQaVtx){
    AliError(Form(RED"This method is available only for container type kQaVtx: you are trying to fill %s!" Bee, strContType(fType)));
  }
  
    Double_t xyz[3],vtxVect[3],jetP[3];
    Double_t xyzPrim[3];
    Double_t cosTheta;

    primVtx->GetXYZ(xyzPrim);
    jet->PxPyPz(jetP);

    Int_t *indexLxy=new Int_t[nvtx];
    Double_t pointVtxProp[21]={mult*1.,jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),jet->GetNumberOfTracks()*1.,0.,0.,0.,0.,p[0],0,0,0,0,0,0,0,0,0,p[1],p[2]};

     Double_t *decLengthXY=new Double_t[nvtx];
     Double_t *invMasses=new Double_t[nvtx];
     Double_t *nRealVtx=new Double_t[nvtx];
     Double_t *nFromBVtx=new Double_t[nvtx];
     Double_t *nFromPromptDVtx=new Double_t[nvtx];
     Double_t xMC,yMC;
     Double_t vtxP[3],vtxPt,signLxy;
     Int_t nvtxMC=0;

     for(Int_t jj=0;jj<nvtx;jj++){
       xMC=-99999.;
       yMC=-99999.;
       AliAODVertex *vtx=(AliAODVertex*)vertices->UncheckedAt(jj);
   
       Double_t chi2ndf=vtx->GetChi2perNDF();
       invMasses[jj]=fTagger->GetVertexInvariantMass(vtx);
       Double_t sigvert=disp[jj];

       nRealVtx[jj]=-1;
       nFromBVtx[jj]=-1;

       vtx->GetXYZ(xyz);
       vtxVect[0]=xyz[0]-xyzPrim[0];
       vtxVect[1]=xyz[1]-xyzPrim[1];
       vtxVect[2]=xyz[2]-xyzPrim[2];
       signLxy=vtxVect[0]*jetP[0]+vtxVect[1]*jetP[1];

       Double_t absJetPt=TMath::Sqrt(jetP[0]*jetP[0]+jetP[1]*jetP[1]);
       Double_t absVtxVect=TMath::Sqrt(vtxVect[0]*vtxVect[0]+vtxVect[1]*vtxVect[1]);
       cosTheta=signLxy/(absJetPt*absVtxVect);//angle between jet and Lxy

       decLengthXY[jj]=TMath::Sqrt((xyz[0]-xyzPrim[0])*(xyz[0]-xyzPrim[0])+(xyz[1]-xyzPrim[1])*(xyz[1]-xyzPrim[1]));
       if(signLxy<0.){
         decLengthXY[jj]*=-1.;
       }
        
       fTagger->GetVtxPxy(vtx,vtxP);
       vtxPt=TMath::Sqrt(vtxP[0]*vtxP[0]+vtxP[1]*vtxP[1]);
       pointVtxProp[5]=vtxPt;
       pointVtxProp[6]=invMasses[jj];
       pointVtxProp[7]=decLengthXY[jj];
       pointVtxProp[8]=chi2ndf;
       
       if(mcPart){
         Int_t nfromBandD=0,nfromD=0,nfromPromptD=0;
         fTagger->GetNTracksFromCommonVertex(vtx,mcPart,nvtxMC,xMC,yMC,nfromBandD,nfromD,nfromPromptD);
         pointVtxProp[10]=nvtxMC;
         pointVtxProp[11]=xyz[0]-xMC;
         pointVtxProp[12]=xyz[1]-yMC;
         pointVtxProp[17]=nfromBandD;
         pointVtxProp[18]=nfromD;

         nRealVtx[jj]=nvtxMC;
         nFromBVtx[jj]=nfromBandD;
         nFromPromptDVtx[jj]=nfromPromptD;
         //if(decLengthXY[jj]>0.2){
         //  }
         }
        
       pointVtxProp[13]=sigvert;
       pointVtxProp[14]=cosTheta*TMath::Abs(decLengthXY[jj]);
       pointVtxProp[15]=TMath::Sqrt(decLengthXY[jj]*decLengthXY[jj]-pointVtxProp[14]*pointVtxProp[14]);
       pointVtxProp[16]=cosTheta;
       TArrayD *apointVtxProp=new TArrayD(21,pointVtxProp);
       AliHFJetsContainer::FillStep(step, apointVtxProp);
       }
          
     delete indexLxy;
     delete decLengthXY;
     delete invMasses;

     return;
}
  
void AliHFJetsContainerVertex::FillStepBJets(AliHFJetsContainer::CFSteps step, Double_t mult, const AliEmcalJet *jet, Int_t nvtx,Double_t partonnat[3], Double_t contribution,Double_t ptpart){
  
  AliInfo("Filling container kBJets.");

  if (fType != kBJets){
    AliError(Form(RED"This method is available only for container type kBJets: you are trying to fill %s!" Bee, strContType(fType)));
  }

  // Double_t point[13]={mult*1., jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),jet->GetNumberOfTracks()*1.,nvtx*1.,0,partonnat[0],contribution,ptpart,jet->EffectiveAreaCharged(), partonnat[1], partonnat[2]};
  Double_t point[13]={mult*1., jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),jet->GetNumberOfTracks()*1.,nvtx*1.,0,partonnat[0],contribution,ptpart,jet->Area(), partonnat[1], partonnat[2]}; // on eff area ch., should check that area is giving the same thing ....todo
  
  //for (Int_t i=0; i<13; i++) Printf("Point %d %f", i, point[i]);
  TArrayD *apoint=new TArrayD(13,point);
  AliHFJetsContainer::FillStep(step, apoint);
}

void AliHFJetsContainerVertex::FillStepJetVtx(AliHFJetsContainer::CFSteps step, Double_t mult, const AliEmcalJet *jet,const TClonesArray *vertices,Int_t nvtx,const AliAODVertex *primVtx,const TClonesArray *mcPart,Double_t p[3],Double_t* disp){
  
  if (fType != kJetVtx){
    AliError(Form(RED"This method is available only for container type kJetVtx: you are trying to fill %s!" Bee, strContType(fType)));
  }
  
  Double_t point[18]={mult*1.,jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),nvtx*1.,0.,-1.,-1.,-1,-1,-1,p[0],p[1],p[2],-1,-1,-1,-1};

    if (vertices && primVtx){ 
       Double_t xyz[3],vtxVect[3],jetP[3];
       Double_t xyzPrim[3];
       //Double_t cosTheta;

       primVtx->GetXYZ(xyzPrim);
       jet->PxPyPz(jetP);

       Int_t *indexLxy=new Int_t[nvtx];
       Double_t *decLenXY=new Double_t[nvtx];
       Double_t *errdecLenXY=new Double_t[nvtx];
       Double_t *sigdecLenXY=new Double_t[nvtx];
       Double_t *invMasses=new Double_t[nvtx];
       Double_t *nRealVtx=new Double_t[nvtx];
       Double_t *nFromBVtx=new Double_t[nvtx];
       Double_t *nFromPromptDVtx=new Double_t[nvtx];
       Double_t *chi2perNDF=new Double_t[nvtx];
       Double_t *sigmavertex=new Double_t[nvtx];

	 Float_t *ipR1 = new Float_t[nvtx];
	 Float_t *ipR2 = new Float_t[nvtx];
	 Float_t *ipR3 = new Float_t[nvtx];

	 Double_t *sigipR1 = new Double_t[nvtx];
	 Double_t *sigipR2 = new Double_t[nvtx];
	 Double_t *sigipR3 = new Double_t[nvtx];

	 Double_t *sumOfSqs = new Double_t[nvtx];

	 Double_t *sigsumOfSqs = new Double_t[nvtx];

       Double_t xMC,yMC;
       Double_t vtxP[3],signLxy;
       Int_t nvtxMC=0;

       for(Int_t jj=0;jj<nvtx;jj++){
         xMC=-99999.;
         yMC=-99999.;
         AliAODVertex *vtx=(AliAODVertex*)vertices->UncheckedAt(jj);
         invMasses[jj]=fTagger->GetVertexInvariantMass(vtx);
         nRealVtx[jj]=-1;
         nFromBVtx[jj]=-1;

         vtx->GetXYZ(xyz);
         vtxVect[0]=xyz[0]-xyzPrim[0];
         vtxVect[1]=xyz[1]-xyzPrim[1];
         vtxVect[2]=xyz[2]-xyzPrim[2];
         signLxy=vtxVect[0]*jetP[0]+vtxVect[1]*jetP[1];
	 
	 chi2perNDF[jj] = vtx->GetChi2perNDF();
	 sigmavertex[jj] = disp[jj];

  AliAODTrack *track0 = (AliAODTrack*)vtx->GetDaughter(0);
  AliAODTrack *track1 = (AliAODTrack*)vtx->GetDaughter(1); 
  AliAODTrack *track2 = (AliAODTrack*)vtx->GetDaughter(2); 

    Double_t dca0aod[2]={-999.,-999.}, dca1aod[2]={-999.,-999.}, dca2aod[2]={-999.,-999.};
    Double_t cov0aod[3]={-999.,-999.,-999.}, cov1aod[3]={-999.,-999.,-999.}, cov2aod[3]={-999.,-999.,-999.};

    track0->PropagateToDCA(primVtx,0.,10000.,dca0aod,cov0aod);
    track1->PropagateToDCA(primVtx,0.,10000.,dca1aod,cov1aod);
    track2->PropagateToDCA(primVtx,0.,10000.,dca2aod,cov2aod);

	 ipR1[jj]= dca0aod[0];
	 ipR2[jj]= dca1aod[0];
	 ipR3[jj]= dca2aod[0];

	 if (cov0aod[0]) sigipR1[jj]= dca0aod[0]/TMath::Sqrt(cov0aod[0]);;
	 if (cov1aod[0]) sigipR2[jj]= dca1aod[0]/TMath::Sqrt(cov1aod[0]);
	 if (cov2aod[0]) sigipR3[jj]= dca2aod[0]/TMath::Sqrt(cov2aod[0]);	 

	 sumOfSqs[jj] = TMath::Sqrt(ipR1[jj]*ipR1[jj] + ipR2[jj]*ipR2[jj] + ipR3[jj]*ipR3[jj]);
	 sigsumOfSqs[jj] = TMath::Sqrt(sigipR1[jj]*sigipR1[jj] + sigipR2[jj]*sigipR2[jj] + sigipR3[jj]*sigipR3[jj]);
	 decLenXY[jj] = primVtx->DistanceXYToVertex(vtx);

         if(signLxy<0.){
           decLenXY[jj]*=-1.;
         }
         errdecLenXY[jj] = primVtx->ErrorDistanceXYToVertex(vtx);
	 sigdecLenXY[jj] = decLenXY[jj]/errdecLenXY[jj];

         fTagger->GetVtxPxy(vtx,vtxP);
       
         if(mcPart){
           Int_t nfromBandD=0,nfromD=0,nfromPromptD=0;
           fTagger->GetNTracksFromCommonVertex(vtx,mcPart,nvtxMC,xMC,yMC,nfromBandD,nfromD,nfromPromptD);

           nRealVtx[jj]=nvtxMC;
           nFromBVtx[jj]=nfromBandD;
           nFromPromptDVtx[jj]=nfromPromptD;
           
           }
        
         }
          
       TMath::Sort(nvtx,decLenXY,indexLxy);
       if(nvtx>0){
         point[6]=decLenXY[indexLxy[0]];
         point[7]=invMasses[indexLxy[0]];
         point[8]=nRealVtx[indexLxy[0]];
         point[9]=nFromBVtx[indexLxy[0]];
         point[10]=nFromPromptDVtx[indexLxy[0]];
	 // point[24]=ipR1[indexLxy[0]];
	 // point[25]=ipR2[indexLxy[0]];
	 // point[26]=ipR3[indexLxy[0]];
	 point[14]=sigmavertex[indexLxy[0]];
	 point[15]=sumOfSqs[indexLxy[0]];
	 point[16]=sigdecLenXY[indexLxy[0]];
	 point[17]=sigsumOfSqs[indexLxy[0]];
         }


       // for now we only take the most dispaced vertex


       // if(nvtx>1){
       //   point[7]=decLenXY[indexLxy[1]];
       //   point[10]=invMasses[indexLxy[1]];
       //   point[13]=nRealVtx[indexLxy[1]];
       //   point[16]=nFromBVtx[indexLxy[1]];
       //   point[19]=nFromPromptDVtx[indexLxy[1]];
       // 	 point[25]=sigmavertex[indexLxy[1]];
       // 	 point[28]=sumOfSqs[indexLxy[1]];
       // 	 // point[31]=errdecLenXY[indexLxy[0]];
       //   }

       // if(nvtx>2){

       //   point[8]=decLenXY[indexLxy[2]];
       //   point[11]=invMasses[indexLxy[2]];
       //   point[14]=nRealVtx[indexLxy[2]];
       //   point[17]=nFromBVtx[indexLxy[2]];
       //   point[20]=nFromPromptDVtx[indexLxy[2]];
       // 	 point[26]=sigmavertex[indexLxy[2]];
       // 	 point[29]=sumOfSqs[indexLxy[2]];
       // 	 // point[32]=errdecLenXY[indexLxy[0]];
       //   }

       // Calculate sum of inv masses of the 3 more displaced vertices
       for(Int_t ivtx=0;ivtx<3;ivtx++){
         if(nvtx>ivtx) point[5]+=invMasses[indexLxy[ivtx]];
         }
     delete indexLxy;
     delete decLenXY;
     delete invMasses;
     } // end if (vertices && primVtx) 

     TArrayD *apoint=new TArrayD(18,point);
     AliHFJetsContainer::FillStep(step, apoint);


     return;
}

void AliHFJetsContainerVertex::FillStepJetVtxData(AliHFJetsContainer::CFSteps step, Double_t mult, const AliEmcalJet *jet, const TClonesArray *vertices, Int_t nvtx,const AliAODVertex *primVtx, Double_t* disp){

  if (fType != kJetVtxData){
    AliError(Form(RED"This method is available only for container type kJetVtxData: you are trying to fill %s!" Bee, strContType(fType)));
  }
  
  Double_t point[9]={mult*1.,jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(),nvtx*1.,-1.,-1.,-1.,-1};

    if (vertices && primVtx){ 
       Double_t xyz[3],vtxVect[3],jetP[3];
       Double_t xyzPrim[3];
       //Double_t cosTheta;

       primVtx->GetXYZ(xyzPrim);
       jet->PxPyPz(jetP);

       Int_t *indexLxy=new Int_t[nvtx];
       Double_t *decLenXY=new Double_t[nvtx];
       Double_t *errdecLenXY=new Double_t[nvtx];
       Double_t *sigdecLenXY=new Double_t[nvtx];
       Double_t *invMasses=new Double_t[nvtx];
       Double_t *chi2perNDF=new Double_t[nvtx];
       Double_t *sigmavertex=new Double_t[nvtx];

	 Float_t *ipR1 = new Float_t[nvtx];
	 Float_t *ipR2 = new Float_t[nvtx];
	 Float_t *ipR3 = new Float_t[nvtx];

	 Double_t *sigipR1 = new Double_t[nvtx];
	 Double_t *sigipR2 = new Double_t[nvtx];
	 Double_t *sigipR3 = new Double_t[nvtx];

	 Double_t *sumOfSqs = new Double_t[nvtx];

	 Double_t *sigsumOfSqs = new Double_t[nvtx];

       Double_t vtxP[3],signLxy;
       Int_t nvtxMC=0;

       for(Int_t jj=0;jj<nvtx;jj++){
         AliAODVertex *vtx=(AliAODVertex*)vertices->UncheckedAt(jj);
         invMasses[jj]=fTagger->GetVertexInvariantMass(vtx);

         vtx->GetXYZ(xyz);
         vtxVect[0]=xyz[0]-xyzPrim[0];
         vtxVect[1]=xyz[1]-xyzPrim[1];
         vtxVect[2]=xyz[2]-xyzPrim[2];
         signLxy=vtxVect[0]*jetP[0]+vtxVect[1]*jetP[1];
	 
	 chi2perNDF[jj] = vtx->GetChi2perNDF();
	 sigmavertex[jj] = disp[jj];

	 decLenXY[jj] = primVtx->DistanceXYToVertex(vtx);


         if(signLxy<0.){
           decLenXY[jj]*=-1.;
         }
         errdecLenXY[jj] = primVtx->ErrorDistanceXYToVertex(vtx);
	 sigdecLenXY[jj] = decLenXY[jj]/errdecLenXY[jj];

         fTagger->GetVtxPxy(vtx,vtxP);

       }
       TMath::Sort(nvtx,decLenXY,indexLxy);
       if(nvtx>0){
         point[5]=decLenXY[indexLxy[0]];
         point[6]=invMasses[indexLxy[0]];
	 point[7]=sigmavertex[indexLxy[0]];
	 // point[15]=sumOfSqs[indexLxy[0]];
	 point[8]=sigdecLenXY[indexLxy[0]];
	 // point[17]=sigsumOfSqs[indexLxy[0]];
         }

    // for now we only take the most dispaced vertex

       // if(nvtx>1){
       //   point[7]=decLenXY[indexLxy[1]];
       //   point[10]=invMasses[indexLxy[1]];
       //   point[13]=nRealVtx[indexLxy[1]];
       //   point[16]=nFromBVtx[indexLxy[1]];
       //   point[19]=nFromPromptDVtx[indexLxy[1]];
       // 	 point[25]=sigmavertex[indexLxy[1]];
       // 	 point[28]=sumOfSqs[indexLxy[1]];
       // 	 // point[31]=errdecLenXY[indexLxy[0]];
       //   }

       // if(nvtx>2){

       //   point[8]=decLenXY[indexLxy[2]];
       //   point[11]=invMasses[indexLxy[2]];
       //   point[14]=nRealVtx[indexLxy[2]];
       //   point[17]=nFromBVtx[indexLxy[2]];
       //   point[20]=nFromPromptDVtx[indexLxy[2]];
       // 	 point[26]=sigmavertex[indexLxy[2]];
       // 	 point[29]=sumOfSqs[indexLxy[2]];
       // 	 // point[32]=errdecLenXY[indexLxy[0]];
       //   }

       // Calculate sum of inv masses of the 3 more displaced vertices
       // for(Int_t ivtx=0;ivtx<3;ivtx++){
       //   if(nvtx>ivtx) point[5]+=invMasses[indexLxy[ivtx]];
       //   }

     delete indexLxy;
     delete decLenXY;
     delete invMasses;
     } // end if (vertices && primVtx) 

     TArrayD *apoint=new TArrayD(9,point);
     AliHFJetsContainer::FillStep(step, apoint);


     return;
}
