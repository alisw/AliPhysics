#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalPythiaInfo.h"
#include "TRandom3.h"



#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalJetShapeExtra.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetShapeExtra)

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapeExtra::AliAnalysisTaskEmcalJetShapeExtra() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetShapeExtra", kTRUE),
  fContainer(0),
  fJetShapeType(kPythiaDef),
  fJetShapeSub(kNoSub),
 // fJetSelection(kInclusive),
  fh2ResponseUW(0x0),
  fPtJet(0x0),
  fNbOfConstvspT(0x0),
    fNumberOfJet(0x0),
fTreeObservableTagging(0)

{
   for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapeExtra::AliAnalysisTaskEmcalJetShapeExtra(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fJetShapeType(kPythiaDef),
fJetShapeSub(kNoSub),
//  fJetSelection(kInclusive),
  fh2ResponseUW(0x0),
  fPtJet(0x0),
  fNbOfConstvspT(0x0),
fNumberOfJet(0x0),
fTreeObservableTagging(0)
  
{
  // Standard constructor.
  for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
  
 DefineOutput(1, TList::Class());
 DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapeExtra::~AliAnalysisTaskEmcalJetShapeExtra()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskEmcalJetShapeExtra::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);



 
  fh2ResponseUW= new TH2F("fh2ResponseUW", "fh2ResponseUW", 100, 0, 200,  100, 0, 200); 
  fOutput->Add(fh2ResponseUW);
  fPtJet= new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);
  

  fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
  fOutput->Add(fNbOfConstvspT);
    
    fNumberOfJet= new TH1F("fNumberOfJet","fNumberOfJet",100,0.5,100.5);
  fOutput->Add(fNumberOfJet);
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

 
  TH1::AddDirectory(oldStatus);
  const Int_t nVar = 16;
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeObservableTagging = new TTree(nameoutput, nameoutput);
  
 
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode"; 
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ptDJet"; 
  fShapesVarNames[3] = "mJet";
  // fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[4] = "angularity";
  fShapesVarNames[5] = "circularity";
  fShapesVarNames[6] = "lesub";
  //fShapesVarNames[6] = "sigma2";

  fShapesVarNames[7] = "ptJetMatch"; 
  fShapesVarNames[8] = "ptDJetMatch"; 
  fShapesVarNames[9] = "mJetMatch";
  // fShapesVarNames[12] = "nbOfConstMatch";
  fShapesVarNames[10] = "angularityMatch";
  fShapesVarNames[11] = "circularityMatch";
  fShapesVarNames[12] = "lesubMatch";
  //fShapesVarNames[12] = "sigma2Match";
  fShapesVarNames[13]="weightPythia";
  fShapesVarNames[14]="ntrksEvt";
//  fShapesVarNames[15]="rhoVal";
//  fShapesVarNames[16]="rhoMassVal";
  fShapesVarNames[15]="ptUnsub";

   for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}
 
  PostData(1,fOutput);
  PostData(2,fTreeObservableTagging);

   delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapeExtra::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapeExtra::FillHistograms()
{
  // Fill histograms.
  //cout<<"base container"<<endl;
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
    Float_t kWeight=1;

  
  AliParticleContainer *partContAn = GetParticleContainer(0);
  TClonesArray *trackArrayAn = partContAn->GetArray();
  Int_t ntracksEvt = trackArrayAn->GetEntriesFast();
  
  
    if(jetCont) {

    jetCont->ResetCurrentID();

      Int_t count=0;
    while((jet1 = jetCont->GetNextAcceptJet())) {
         count++;
      if (!jet1) continue;
      AliEmcalJet* jet3 = 0x0;
      fPtJet->Fill(jet1->Pt());
  
 //       cout<<"Jet detector level Pt is"<<jet1->Pt()<<endl;


      
      fShapesVar[0] = 0.;


      if (fJetShapeType == kPythiaDef){

        if(!(fJetShapeSub==kConstSub)) jet3 = jet1->ClosestJet();
        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }
        
  //        cout<<"jet3 particle level pt is"<<jet3->Pt()<<endl;
        fh2ResponseUW->Fill(jet1->Pt(),jet3->Pt());
        
        
      }
      
      
        
    
    
      
      
      fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), jet1->Pt());
      if (jet1->GetNumberOfTracks() <= 1) continue;

  
      fShapesVar[1] = jet1->Pt();
      fShapesVar[2] = GetJetpTD(jet1,0);
      fShapesVar[3] = GetJetMass(jet1,0);
      fShapesVar[4] = GetJetAngularity(jet1,0);
      fShapesVar[5] = GetJetCircularity(jet1,0);
      fShapesVar[6] = GetJetLeSub(jet1,0);
 
 
      Float_t ptMatch=0., ptDMatch=0., massMatch=0.,angulMatch=0.,circMatch=0., lesubMatch=0.;
      Int_t kMatched = 0;

       if (fJetShapeType==kPythiaDef) {
         kMatched =1;
         if(fJetShapeSub==kConstSub) kMatched = 3;
        
         ptMatch=jet3->Pt();
         ptDMatch=GetJetpTD(jet3, kMatched);
         massMatch=GetJetMass(jet3,kMatched);
         //constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
         angulMatch=GetJetAngularity(jet3, kMatched);
        circMatch=GetJetCircularity(jet3, kMatched);
         lesubMatch=GetJetLeSub(jet3, kMatched);
         //sigma2Match = GetSigma2(jet2, kMatched);
       }
      


       
      if (fJetShapeType == kData ) {
        kMatched = 0;
        ptMatch=0.;
        ptDMatch=0.;
        massMatch=0.;
	//constMatch=0.;
        angulMatch=0.;
        circMatch=0.;
        lesubMatch=0.;
        //sigma2Match =0.;
        
      }
      
    

      fShapesVar[7] = ptMatch;
      fShapesVar[8] = ptDMatch;
      fShapesVar[9] = massMatch;
      fShapesVar[10] = angulMatch;
      fShapesVar[11] = circMatch;
      fShapesVar[12] = lesubMatch;
      fShapesVar[13] = kWeight;
      fShapesVar[14] = ntracksEvt;
//      fShapesVar[15] = rhoVal;
//      fShapesVar[16] = rhoMassVal;
      fShapesVar[15] = jet1->Pt();


      fTreeObservableTagging->Fill();
      





    }

        fNumberOfJet->Fill(count);

  }
  
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetMass(AliEmcalJet *jet,Int_t jetContNb){
   
        if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
            return jet->GetShapeProperties()->GetFirstOrderSubtracted();
        else
            return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::Angularity(AliEmcalJet *jet, Int_t jetContNb){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
        Double_t dr = TMath::Sqrt(dr2);
        num=num+vp1->Pt()*dr;
        den=den+vp1->Pt();
    }
    return num/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb ){
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
    return Angularity(jet, jetContNb);
    
}


//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::PTD(AliEmcalJet *jet, Int_t jetContNb ){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        num=num+vp1->Pt()*vp1->Pt();
        den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb ){
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
        return PTD(jet, jetContNb);
    
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::Circularity(AliEmcalJet *jet, Int_t jetContNb ){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    Double_t pxjet=jet->Px();
    Double_t pyjet=jet->Py();
    Double_t pzjet=jet->Pz();
    
    
    //2 general normalized vectors perpendicular to the jet
    TVector3  ppJ1(pxjet, pyjet, pzjet);
    TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
    ppJ3.SetMag(1.);
    TVector3  ppJ2(-pyjet, pxjet, 0);
    ppJ2.SetMag(1.);
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        TVector3 pp(vp1->Px(), vp1->Py(), vp1->Pz());
        
        //local frame
        TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
        TVector3 pPerp = pp - pLong;
        //projection onto the two perpendicular vectors defined above
        
        Float_t ppjX = pPerp.Dot(ppJ2);
        Float_t ppjY = pPerp.Dot(ppJ3);
        Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
        if(ppjT<=0) return 0;
        
        mxx += (ppjX * ppjX / ppjT);
        myy += (ppjY * ppjY / ppjT);
        mxy += (ppjX * ppjY / ppjT);
        nc++;
        sump2 += ppjT;}
    
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t circ=0;
    if(jev==1) circ=2*eval[0];
    if(jev==0) circ=2*eval[1];
    
    return circ;
    
    
    
}




//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb  ){
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
            return Circularity(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::LeSub(AliEmcalJet *jet, Int_t jetContNb  ){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    AliVParticle *vp2 = 0x0;
    std::vector<int> ordindex;
    ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
    //Printf("Nbof const = %d", jet->GetNumberOfTracks());
    //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
    
    if(ordindex.size()<2) return -1;
    
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
    if (!vp1){
        Printf("AliVParticle associated to Leading constituent not found");
        return -1;
    }
    
    vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
    if (!vp2){
        Printf("AliVParticle associated to Subleading constituent not found");
        return -1;
    }
    
    
    num=vp1->Pt();
    den=vp2->Pt();
    //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
    
    return num-den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb ) {
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
            return LeSub(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb){
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
    return jet->GetNumberOfTracks();
    
}


//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::Sigma2(AliEmcalJet *jet, Int_t jetContNb){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t ppt=vp1->Pt();
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        
        Double_t deta = vp1->Eta()-jet->Eta();
        mxx += ppt*ppt*deta*deta;
        myy += ppt*ppt*dphi*dphi;
        mxy -= ppt*ppt*deta*TMath::Abs(dphi);
        nc++;
        sump2 += ppt*ppt;
        
    }
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx , mxy , mxy , myy };
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t sig=0;
    if(jev==1) sig=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
    if(jev==0) sig=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
    
    return sig;
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapeExtra::GetSigma2(AliEmcalJet *jet, Int_t jetContNb){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
        return Sigma2(jet, jetContNb);
    
}


//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetShapeExtra::RelativePhi(Double_t mphi,Double_t vphi){
    
    if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
    else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
    if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
    else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
    double dphi = mphi-vphi;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapeExtra::RetrieveEventObjects() {
    //
    // retrieve event objects
    //
    if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
        return kFALSE;
    
    return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetShapeExtra::Terminate(Option_t *)
{
    // Called once at the end of the analysis.
    
    // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
    // if (!fTreeObservableTagging){
    //   Printf("ERROR: fTreeObservableTagging not available");
    //   return;
    // }
    
}

