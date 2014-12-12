//
// Jet mass structure analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TMath.h>
#include <TRandom3.h>

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
#include "AliParticleContainer.h"
#include "AliJetContainer.h"
#include "AliEmcalJetByJetCorrection.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalJetMassStructure.h"

ClassImp(AliAnalysisTaskEmcalJetMassStructure)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassStructure::AliAnalysisTaskEmcalJetMassStructure() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMassStructure", kTRUE),
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassType(kRaw),
  fRandom(0),
  fEfficiencyFixed(1.),
  fCorrType(kMeanPtR),
  fEJetByJetCorr(0),
  fh3PtDRMass(0),
  fh3PtDRRho(0),
  fh3PtDRMassCorr(0),
  fh3PtDRRhoCorr(0),
  fh2PtMass(0),
  fh2PtMassCorr(0),
  fhnMassResponse(0),
  fhnMassResponseCorr(0),
  fh3JetPtDRTrackPt(0)
{
  // Default constructor.

  fh3PtDRMass                       = new TH3F*[fNcentBins];
  fh3PtDRRho                        = new TH3F*[fNcentBins];
  fh3PtDRMassCorr                   = new TH3F*[fNcentBins];
  fh3PtDRRhoCorr                    = new TH3F*[fNcentBins];
  fh2PtMass                         = new TH2F*[fNcentBins];
  fh2PtMassCorr                     = new TH2F*[fNcentBins];
  fh3JetPtDRTrackPt                 = new TH3F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtDRMass[i]       = 0;
    fh3PtDRRho[i]        = 0;
    fh3PtDRMassCorr[i]   = 0;
    fh3PtDRRhoCorr[i]    = 0;
    fh2PtMass[i]         = 0;
    fh2PtMassCorr[i]     = 0;
    fh3JetPtDRTrackPt[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassStructure::AliAnalysisTaskEmcalJetMassStructure(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassType(kRaw),
  fRandom(0),
  fEfficiencyFixed(1.),
  fCorrType(kMeanPtR),
  fEJetByJetCorr(0),
  fh3PtDRMass(0),
  fh3PtDRRho(0),
  fh3PtDRMassCorr(0),
  fh3PtDRRhoCorr(0),
  fh2PtMass(0),
  fh2PtMassCorr(0),
  fhnMassResponse(0),
  fhnMassResponseCorr(0),
  fh3JetPtDRTrackPt(0)
{
  // Standard constructor.

  fh3PtDRMass                       = new TH3F*[fNcentBins];
  fh3PtDRRho                        = new TH3F*[fNcentBins];
  fh3PtDRMassCorr                   = new TH3F*[fNcentBins];
  fh3PtDRRhoCorr                    = new TH3F*[fNcentBins];
  fh2PtMass                         = new TH2F*[fNcentBins];
  fh2PtMassCorr                     = new TH2F*[fNcentBins];
  fh3JetPtDRTrackPt                 = new TH3F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtDRMass[i]       = 0;
    fh3PtDRRho[i]        = 0; 
    fh3PtDRMassCorr[i]   = 0;
    fh3PtDRRhoCorr[i]    = 0;
    fh2PtMass[i]         = 0;
    fh2PtMassCorr[i]     = 0;
    fh3JetPtDRTrackPt[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassStructure::~AliAnalysisTaskEmcalJetMassStructure()
{
  // Destructor.
  delete fRandom;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassStructure::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if(!fRandom) fRandom = new TRandom3(0);

  if(fCorrType==kMeanPtR && fEJetByJetCorr) fEJetByJetCorr->Init();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 200;
  const Double_t minPt = -50.;
  const Double_t maxPt = 150.;
  Double_t *binsPt = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binsPt[i]=(Double_t)minPt + (maxPt-minPt)/nBinsPt*(Double_t)i ;

  const Int_t nBinsM  = 120;
  const Double_t minM = -20.;
  const Double_t maxM = 40.;
  Double_t *binsM = new Double_t[nBinsM+1];
  for(Int_t i=0; i<=nBinsM; i++) binsM[i]=(Double_t)minM + (maxM-minM)/nBinsM*(Double_t)i ;

  const Int_t nBinsR  = 80;
  const Double_t minR = -0.005;
  const Double_t maxR = 0.795;
  Double_t *binsR = new Double_t[nBinsR+1];
  for(Int_t i=0; i<=nBinsR; i++) binsR[i]=(Double_t)minR + (maxR-minR)/nBinsR*(Double_t)i ;

  //track pt
  Double_t binWidth1 = 0.1;
  Double_t binWidth2 = 1.;
  Double_t binWidth3 = 1.;
  const Float_t ptmin1 =  0.;
  const Float_t ptmax1 =  5.;
  const Float_t ptmin2 =  ptmax1;
  const Float_t ptmax2 =  10.;
  const Float_t ptmin3 =  ptmax2 ;
  const Float_t ptmax3 =  100.;
  const Int_t nbin11 = (int)((ptmax1-ptmin1)/binWidth1);
  const Int_t nbin12 = (int)((ptmax2-ptmin2)/binWidth2)+nbin11;
  const Int_t nbin13 = (int)((ptmax3-ptmin3)/binWidth3)+nbin12;
  Int_t nBinsPtTr=nbin13;
  //Create array with low edges of each bin
  Double_t *binsPtTr=new Double_t[nBinsPtTr+1];
  for(Int_t i=0; i<=nBinsPtTr; i++) {
    if(i<=nbin11) binsPtTr[i]=(Double_t)ptmin1 + (ptmax1-ptmin1)/nbin11*(Double_t)i ;
    if(i<=nbin12 && i>nbin11) binsPtTr[i]=(Double_t)ptmin2 + (ptmax2-ptmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;
    if(i<=nbin13 && i>nbin12) binsPtTr[i]=(Double_t)ptmin3 + (ptmax3-ptmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;
  }

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 5;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt};

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = TString::Format("fh3PtDRMass_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};r;sum m",histName.Data());
    fh3PtDRMass[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,101,-0.005,1.005);
    fOutput->Add(fh3PtDRMass[i]);

    histName = TString::Format("fh3PtDRRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};r;sum #it{p}_{T}",histName.Data());
    fh3PtDRRho[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,101,-0.005,1.005);
    fOutput->Add(fh3PtDRRho[i]);

    histName = TString::Format("fh3PtDRMassCorr_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};r;sum m",histName.Data());
    fh3PtDRMassCorr[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,101,-0.005,1.005);
    fOutput->Add(fh3PtDRMassCorr[i]);

    histName = TString::Format("fh3PtDRRhoCorr_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};r;sum #it{p}_{T}",histName.Data());
    fh3PtDRRhoCorr[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,101,-0.005,1.005);
    fOutput->Add(fh3PtDRRhoCorr[i]);

    histName = TString::Format("fh2PtMass_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};#it{M}_{jet}",histName.Data());
    fh2PtMass[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtMass[i]);

    histName = TString::Format("fh2PtMassCorr_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};#it{M}_{jet}",histName.Data());
    fh2PtMassCorr[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtMassCorr[i]);

    histName = TString::Format("fh3JetPtDRTrackPt_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet};r;#it{p}_{T,track}",histName.Data());
    fh3JetPtDRTrackPt[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,binsPt,nBinsR,binsR,nBinsPtTr,binsPtTr);
    fOutput->Add(fh3JetPtDRTrackPt[i]);
  }

  histName = "fhnMassResponse";
  histTitle = Form("%s;#it{M}_{det};#it{M}_{part};#it{p}_{T,det};#it{p}_{T,part}",histName.Data());
  fhnMassResponse = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
  fOutput->Add(fhnMassResponse);

  histName = "fhnMassResponseCorr";
  histTitle = Form("%s;#it{M}_{det};#it{M}_{part};#it{p}_{T,det};#it{p}_{T,part}",histName.Data());
  fhnMassResponseCorr = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
  fOutput->Add(fhnMassResponseCorr);

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  if(binsPt)                delete [] binsPt;
  if(binsPtTr)              delete [] binsPtTr;
  if(binsM)                 delete [] binsM;
  if(binsR)                 delete [] binsR;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassStructure::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassStructure::FillHistograms()
{
  // study how jet mass builds up as function of radial distance to jet axis

  if(fCorrType==kMeanPtR && !fEJetByJetCorr) return kFALSE;

  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
      Double_t mJet1 = GetJetMass(jet1);
      fh2PtMass[fCentBin]->Fill(ptJet1,mJet1);
      AliEmcalJet *jPart = jet1->ClosestJet(); 
      //fill detector response
      if(jPart) {
        Double_t var[4] = {mJet1,jPart->M(),ptJet1,jPart->Pt()};
        fhnMassResponse->Fill(var);
      }

      // if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
      // 	continue;
     
      //Sort array based on distance to jet axis
      static Int_t indexes[9999] = {-1};
      static Float_t dr[9999] = {0};
      if(!GetSortedArray(indexes,dr,jet1)) continue;

      for(Int_t i=jet1->GetNumberOfTracks()-1; i>-1; i--) {
        AliVParticle *vp = static_cast<AliVParticle*>(jet1->TrackAt(indexes[i], fTracks));
        if(!vp) continue;
        fh3JetPtDRTrackPt[fCentBin]->Fill(ptJet1,dr[indexes[i]],vp->Pt());
      }
 
      if(fCorrType==kAnnulus) {
        TLorentzVector sumVec;      sumVec.SetPtEtaPhiM(0.,0.,0.,0.);
        TLorentzVector corrVec;     corrVec.SetPtEtaPhiM(0.,0.,0.,0.);
        TLorentzVector curVec;
        for(Int_t i=jet1->GetNumberOfTracks()-1; i>-1; i--) {
          AliVParticle *vp = static_cast<AliVParticle*>(jet1->TrackAt(indexes[i], fTracks));
          if(!vp) continue;
          
          curVec.SetPtEtaPhiM(vp->Pt(),vp->Eta(),vp->Phi(),vp->M());
          sumVec+=curVec;
          corrVec+=curVec;
          //	Printf("%d  %f",i,dr[indexes[i]]);
          fh3PtDRMass[fCentBin]->Fill(ptJet1,dr[indexes[i]],sumVec.M()/mJet1);
          fh3PtDRRho[fCentBin]->Fill(ptJet1,dr[indexes[i]],sumVec.Pt()/ptJet1);
          
          fh3PtDRMassCorr[fCentBin]->Fill(ptJet1,dr[indexes[i]],corrVec.M()/mJet1);
          fh3PtDRRhoCorr[fCentBin]->Fill(ptJet1,dr[indexes[i]],corrVec.Pt()/ptJet1);
          
          Double_t eff = GetEfficiency(vp->Pt());
          Double_t rnd = fRandom->Uniform(1.);
          if(rnd>eff) {//put particle back at random position in annulus; using now zero width for annulus
            Double_t t = TMath::TwoPi()*gRandom->Uniform(1.);
            Double_t rr = dr[indexes[i]];
            rr = fRandom->Uniform(0.,jetCont->GetJetRadius());
            Double_t deta = rr*TMath::Cos(t);
            Double_t dphi = rr*TMath::Sin(t);
            curVec.SetPtEtaPhiM(vp->Pt(),deta+jet1->Eta(),dphi+jet1->Phi(),vp->M());
            corrVec+=curVec;
            
            fh3PtDRMassCorr[fCentBin]->Fill(ptJet1,dr[indexes[i]],corrVec.M()/mJet1);
            fh3PtDRRhoCorr[fCentBin]->Fill(ptJet1,dr[indexes[i]],corrVec.Pt()/ptJet1);
          }
        }//track loop

        fh2PtMassCorr[fCentBin]->Fill(corrVec.Pt(), corrVec.M());
        if(jPart) {
          Double_t varCorr[4] = {corrVec.M(),jPart->M(),corrVec.Pt(),jPart->Pt()};
          fhnMassResponseCorr->Fill(varCorr);
        }
      }//kAnnulus method
      else if(fCorrType==kMeanPtR) {
        //jet-by-jet correction based on templates
        AliEmcalJet *jetCorr = fEJetByJetCorr->Eval(jet1,jetCont->GetParticleContainer()->GetArray());
        if(jPart && jetCorr) {
          Double_t varCorr[4] = {jetCorr->M(),jPart->M(),jetCorr->Pt(),jPart->Pt()};
          fhnMassResponseCorr->Fill(varCorr);
        }
      }//kMeanPtR method
    }//jet loop
  }//jet container
  
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassStructure::GetEfficiency(Double_t pt) {
  pt = 1.*pt;
  Double_t eff = 1.;
  if(fEfficiencyFixed<1.) return fEfficiencyFixed;
  
  return eff;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassStructure::GetSortedArray(Int_t indexes[], Float_t dr[], AliEmcalJet *jet) const
{
  // sort array
  const Int_t n = (Int_t)jet->GetNumberOfTracks();//(Int_t)array.size();
  if (n < 1) return kFALSE;
  
  for (Int_t i = 0; i < n; i++) {
    AliVParticle *vp = static_cast<AliVParticle*>(jet->TrackAt(i, fTracks));
    if(!vp) continue;
    dr[i] = jet->DeltaR(vp);
  }
  TMath::Sort(n, dr, indexes);
  return kTRUE;
}


//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassStructure::GetJetMass(AliEmcalJet *jet) {
  //calc subtracted jet mass
  if(fJetMassType==kRaw)
    return jet->M();
  else if(fJetMassType==kDeriv)
    return jet->GetSecondOrderSubtracted();
  
  return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassStructure::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetMassStructure::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

