#include "TList.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisUtils.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliVEventHandler.h"
#include "AliStack.h"
#include "TGraph.h"
#include "AliTracker.h"
#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "AliPicoTrack.h"
#include "TMath.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisTaskHFJetIPQA.h"
#include "AliExternalTrackParam.h"
#include "AliVertexerTracks.h"
#include "AliEmcalList.h"
#include "AliHFJetsTagging.h"
#include "AliESDUtils.h"
#include "AliKFParticle.h"
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "AliTriggerAnalysis.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "AliOADBContainer.h"
#include "TFile.h"
#include "TGrid.h"
using std::min;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
ClassImp(AliAnalysisTaskHFJetIPQA)
// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(): AliAnalysisTaskEmcalJet("AliAnalysisTaskHFJetIPQA", kFALSE),
  fGraphMean(0x0),
  fGraphSigmaData(0x0),
  fGraphSigmaMC(0x0),
  fGraphXi(0x0),
  fGraphOmega(0x0),
  fK0Star(0x0),
  fPhi(0x0),
  fGeant3FlukaProton(0x0),
  fGeant3FlukaAntiProton(0x0),
  fGeant3FlukaLambda(0x0),
  fGeant3FlukaAntiLambda(0x0),
  fGeant3FlukaKMinus(0x0),
  fOutput2(0x0),
  fMCArray(0x0),
  fJetCutsHF(new AliRDHFJetsCuts()),
  fAODBcont(0x0),
  fMCEvent(0x0),
  fESDTrackCut(0x0),
  fUtils(new AliAnalysisUtils()),
  fESD(kFALSE),
  fMcEvtSampled(kFALSE),
  fCorrrectionSamplingMode(kFALSE),
  fItsClustersInputGlobal(6),
  fBackgroundFactorLinus{0},
  fEtaSEvt(100),
  fPhiSEvt(100),
  fEtaBEvt(100),
  fPhiBEvt(100),
  fEtaCEvt(100),
  fPhiCEvt(100),
  fEtaUdsgEvt(100),
  fPhiUdsgEvt(100),
  fh2dAcceptedTracksEtaPhiPerLayer{0,0,0,0,0,0}
{
  fOutput2 =0x0;
  DefineOutput(1,  TList::Class()) ;
  for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1.;
}
//###############################################################################################################
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name): AliAnalysisTaskEmcalJet(name,kFALSE),
  fGraphMean(0x0),
  fGraphSigmaData(0x0),
  fGraphSigmaMC(0x0),
  fGraphXi(0x0),
  fGraphOmega(0x0),
  fK0Star(0x0),
  fPhi(0x0),
  fGeant3FlukaProton(0x0),
  fGeant3FlukaAntiProton(0x0),
  fGeant3FlukaLambda(0x0),
  fGeant3FlukaAntiLambda(0x0),
  fGeant3FlukaKMinus(0x0),
  fOutput2(0x0),
  fMCArray(0x0),
  fJetCutsHF(new AliRDHFJetsCuts()),
  fAODBcont(0x0),
  fMCEvent(0x0),
  fESDTrackCut(0x0),
  fUtils(new AliAnalysisUtils()),
  fESD(kFALSE),
  fMcEvtSampled(kFALSE),
  fCorrrectionSamplingMode(kFALSE),
  fItsClustersInputGlobal(6),
  fBackgroundFactorLinus{0},
  fEtaSEvt(100),
  fPhiSEvt(100),
  fEtaBEvt(100),
  fPhiBEvt(100),
  fEtaCEvt(100),
  fPhiCEvt(100),
  fEtaUdsgEvt(100),
  fPhiUdsgEvt(100),
  fh2dAcceptedTracksEtaPhiPerLayer{0,0,0,0,0,0}
{
  fOutput2 =0x0;
  DefineOutput(1,  TList::Class()) ;
  for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1; //set default to 1
}

// ########################################################################################
void AliAnalysisTaskHFJetIPQA::SetAODBContainer(AliOADBContainer* cont) {
  fAODBcont =(AliOADBContainer*)(cont->Clone("fAODBcont"));
  return;
}
// ########################################################################################
void AliAnalysisTaskHFJetIPQA::setN_ITSClusters_Input_global( Int_t value)
{
  fItsClustersInputGlobal = value;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::Notify()
{
  if(fCorrrectionSamplingMode)  return AliAnalysisTaskEmcalJet::Notify();
  if(!fAODBcont) {
      Printf("fAODBcont not set");
      return kFALSE;
    }
  Int_t runnr =-1;
  runnr = AliAnalysisManager::GetAnalysisManager()->GetRunFromPath();
  const char * pass = "pass4_mean";
  if(fIsPythia) pass = "pass4mc_mean";
  if(fGraphMean) fGraphMean->Delete();
  fGraphMean =  new TGraph(((TH1D*)(fAODBcont->GetObject(runnr,"",pass))));
  if(fIsPythia){
      if(fGraphSigmaMC) fGraphSigmaMC->Delete();
      fGraphSigmaMC =new  TGraph(((TH1D*)(fAODBcont->GetObject(runnr,"","pass4mc_sigma"))));
      if(fGraphSigmaData) fGraphSigmaData->Delete();
      fGraphSigmaData = new TGraph(((TH1D*)(fAODBcont->GetObject(runnr,"","pass4_sigma"))));
    }
  return AliAnalysisTaskEmcalJet::Notify();
}
// ########################################################################################  Main Loop
Bool_t AliAnalysisTaskHFJetIPQA::Run(){
  fEtaBEvt.clear();
  fPhiBEvt.clear();
  fEtaCEvt.clear();
  fPhiCEvt.clear();
  fEtaUdsgEvt.clear();
  fPhiUdsgEvt.clear();
  fEtaSEvt.clear();
  fPhiSEvt.clear();
  fMcEvtSampled = kFALSE;

  if(fESD) fIsEsd = kTRUE;
  fMCArray         = NULL;
  fMCEvent         = NULL;
  if(fIsPythia && !fESD) fMCArray= dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
  else if(fIsPythia)  fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());

  Double_t weight =1;
  Int_t nTracksInEvent = 0 ;
  nTracksInEvent  = InputEvent()->GetNumberOfTracks() ;
  AliVTrack * trackV =NULL;
  Double_t dca[2] = {-99999,-99999};
  Double_t cov[3] = {-99999,-99999,-99999};
  //Global Track Loop
  for(long itrack= 0; itrack<nTracksInEvent;++itrack)
    {
      trackV=(AliVTrack*)(InputEvent()->GetTrack(itrack));
      if(!trackV) continue;
      //Fill Track stats
      IncHist("fh1dTracksAccepeted",1);
      if(!fESD && !((AliAODTrack*)trackV)->TestFilterBit(1 << 4)) continue;
      if(!IsTrackAccepted(trackV,fItsClustersInputGlobal)) {
          IncHist("fh1dTracksAccepeted",3);
          continue;
        }
      IncHist("fh1dTracksAccepeted",2);
      if(!fCorrrectionSamplingMode)SmearTrackHybrid(trackV); //Hybrid  approach
      FillHist("fh2dAcceptedTracksEtaPhi",trackV->Eta(),trackV->Phi(),1);
      //Calculate impact parameters and fill histogramm
      Bool_t hasIPSuccess = kFALSE;
      dca[0]=-9999;
      dca[1]=-9999;
      cov[0]=-9999;
      cov[1]=-9999;
      cov[2]=-9999;

      if (CalculateTrackImpactParameter(trackV,dca,cov))hasIPSuccess =kTRUE;
      if(!hasIPSuccess)  continue;
      else {
          for(Int_t it = 0; it <1;++it){
              Double_t x= trackV->Pt();
              if(!fCorrrectionSamplingMode) SubtractMean (dca,trackV);
              if(!fCorrrectionSamplingMode && fIsPythia) {
                  //Add delta Sigma from smearing
                  Double_t xpt = trackV->Pt();
                  // if(xpt <1.5)  xpt = 1.5;
                  if(xpt >8)  xpt = 8;
                  Double_t valmc   = fGraphSigmaMC->Eval(xpt);
                  Double_t valdata =   fGraphSigmaData->Eval(xpt);

                  valmc *=1e-4;
                  valdata *=1e-4;
                  Double_t deltaSigma = valdata*valdata-valmc*valmc;
                  cov[0] =  cov[0] + deltaSigma;
                }

              weight =1;
              FillHist("fh2dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),trackV->Pt(),1);
              FillHist("fh2dTracksImpParZ",dca[1],trackV->Pt(),1);
              FillHist("fh2dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),trackV->Pt(),1);
              FillHist("fh2dTracksImpParZSignificance",GetValImpactParameter(kZSig,dca,cov),trackV->Pt(),1);
              FillHist("fh1dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),1);
              FillHist("fh1dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),1);
              FillHist("fh1dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),1);
              FillHist("fh1dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),1);
              if(fIsPythia){
                  Int_t corrpartidx =-1;
                  weight = GetMonteCarloCorrectionFactor(trackV,corrpartidx);
                  FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),trackV->Pt(),weight);
                  FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight);
                  FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight);
                  FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight);
                  FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight);
                }
            }
        }
    }
  if(fCorrrectionSamplingMode) return kTRUE;
  //Here Beginns the jetpart;
  AliJetContainer * jetconrec = 0x0;
  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer * jetcongen = 0x0;
  AliEmcalJet * jetgen  = 0x0;
  if(fIsPythia){
      jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
      if(!MatchJetsGeometricDefault()) cout << "Error running jet matching!" << endl;
      jetcongen->ResetCurrentID();
      // Fill gen. level jet histograms
      while ((jetgen = jetcongen->GetNextJet()))
        {
          if (!jetgen) continue;
          Int_t jetflavour =0;
          AliAODMCParticle* parton = NULL;

          //Only ESD right now
          Bool_t is_udgjet = kFALSE;
          jetflavour =IsMCJetPartonFast(jetgen,0.4,is_udgjet); //Event based association to save memory
          FillHist("fh1dJetGenPt",GetPtCorrectedMC(jetgen),1);
          if(jetflavour ==0)
            FillHist("fh1dJetGenPtUnidentified",GetPtCorrectedMC(jetgen),1);
          else if(jetflavour ==1)
            FillHist("fh1dJetGenPtudsg",GetPtCorrectedMC(jetgen),1);
          else if(jetflavour ==2)
            FillHist("fh1dJetGenPtc",GetPtCorrectedMC(jetgen),1);
          else if(jetflavour ==3)
            FillHist("fh1dJetGenPtb",GetPtCorrectedMC(jetgen),1);
        }

      jetcongen->ResetCurrentID();
      jetconrec->ResetCurrentID();
    }

  // loop rec level jets
  AliEmcalJet * jetrec  = 0x0;
  AliEmcalJet * jetmatched  = 0x0;
  jetconrec->ResetCurrentID();
  Double_t jetpt=0;

  while ((jetrec = jetconrec->GetNextJet()))
    {
      jetpt = jetrec->Pt();
      if(!(jetconrec->GetRhoParameter() == 0x0))
        {
          jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
        }

      if(fIsPythia){
          if (jetrec->MatchedJet()) {
              Double_t genpt = jetrec->MatchedJet()->Pt();
              if(!(jetcongen->GetRhoParameter() == 0x0)){
                  genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
                }
              FillHist("fh2dJetGenPtVsJetRecPt",genpt,jetpt,1);
            }
        }
      // make inclusive signed imp. parameter constituent histograms
      Int_t ntracks = (Int_t)jetrec->GetNumberOfTracks();
      Double_t dca[2] = {-99999,-99999};
      Double_t cov[3] = {-99999,-99999,-99999};
      Double_t sign=0;
      Int_t jetflavour=0;
      Bool_t is_udgjet = kFALSE;
      if(fIsPythia){
          jetmatched = 0x0;
          jetmatched =jetrec->MatchedJet();
          if(jetmatched){
              if(fESD ){
                  jetflavour =IsMCJetPartonFast(jetmatched,0.4,is_udgjet); //Event based association to save memory
                }

            }
        }

      FillHist("fh1dJetRecPt",jetrec->Pt(),1);
      if(fIsPythia){
          if(jetflavour==0) FillHist("fh1dJetRecPtUnidentified",jetpt,1);
          else if(jetflavour==1)FillHist("fh1dJetRecPtudsg",jetpt,1);
          else if(jetflavour==2)FillHist("fh1dJetRecPtc",jetpt,1);
          else if(jetflavour==3)FillHist("fh1dJetRecPtb",jetpt,1);
        }
      if(!(fJetCutsHF->IsJetSelected(jetrec))) continue;
      FillHist("fh1dJetRecEtaPhiAccepted",jetrec->Eta(),jetrec->Phi(),1);
      FillHist("fh1dJetRecPtAccepted",jetpt,1);
      if(fIsPythia){
          if(jetflavour==0)   FillHist("fh1dJetRecPtUnidentifiedAccepted",jetpt,1);
          else if(jetflavour==1)FillHist("fh1dJetRecPtudsgAccepted",jetpt,1);
          else if(jetflavour==2)FillHist("fh1dJetRecPtcAccepted",jetpt,1);
          else if(jetflavour==3)FillHist("fh1dJetRecPtbAccepted",jetpt,1);
        }
      if(fDoJetProbabilityAnalysis){
          //  Printf("Trying to calculate jet probability");
          Double_t jetprob = -999;
          jetprob = CalculateJetProb(jetrec);
          // Printf("%e",jetprob);
          if(jetprob >=0) {
              const char * flavour[5]  = {"Unidentified","udsg","c","b",""};

              if(fIsPythia)FillHist(Form("JetProbability_%s",flavour[jetflavour]),jetpt,jetprob);
              FillHist(Form("JetProbability_%s",flavour[4]),jetpt,jetprob);
            }
        }
      std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;
      for(Int_t itrack = 0; itrack < ntracks; ++itrack)
        {
          Double_t dcatrackjet =999;
          Double_t lineardecaylenth = 999;
          AliVTrack * trackV = (((AliVTrack*)((jetconrec->GetParticleContainer())->GetParticle(jetrec->TrackAt(itrack)))));
          if(!trackV) continue;
          if (!IsTrackAccepted(trackV,fItsClustersInputGlobal)) continue;
          Bool_t hasSIP =kFALSE;
          if(!fCorrrectionSamplingMode) SmearTrackHybrid(trackV);
          if(CalculateJetSignedTrackImpactParameter(trackV,jetrec,dca,cov,sign,dcatrackjet,lineardecaylenth))hasSIP =kTRUE;

          if(hasSIP)
            {
              if(!fCorrrectionSamplingMode) SubtractMean (dca,trackV);
              if(!fCorrrectionSamplingMode && fIsPythia) {
                  //Add delta Sigma from smearing
                  Double_t xpt = trackV->Pt();
                  if(xpt >8)  xpt = 8;

                  Double_t valmc   = fGraphSigmaMC->Eval(xpt);
                  Double_t valdata =   fGraphSigmaData->Eval(xpt);

                  valmc *=1e-4;
                  valdata *=1e-4;
                  Double_t deltaSigma = valdata*valdata-valmc*valmc;
                  cov[0] =  cov[0] + deltaSigma;
                }


              if(abs(dca[0])>1.) continue;
              if(abs(dca[1])>2.) continue;
              if(lineardecaylenth > 10.) continue;
              if (dcatrackjet > 0.07) continue;
              Double_t cursImParXY =TMath::Abs(GetValImpactParameter(kXY,dca,cov))*sign;
              Double_t cursImParXYZ =TMath::Abs(GetValImpactParameter(kXYZ,dca,cov))*sign;
              Double_t cursImParXYSig =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
              Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;
              Int_t corridx=-1;
              if(fIsPythia){
                  weight = GetMonteCarloCorrectionFactor(trackV,corridx);
                }
              else weight=1;
              FillHist("fh2dJetSignedImpParXY",jetpt,cursImParXY,weight);
              FillHist("fh2dJetSignedImpParXYZ",jetpt,cursImParXYZ,weight);
              FillHist("fh2dJetSignedImpParXYSignificance",jetpt,cursImParXYSig,weight);
              FillHist("fh2dJetSignedImpParXYZSignificance",jetpt,cursImParXYZSig,weight);
              const char * subtype [4] = {"Unidentified","udsg","c","b"};
              if(fIsPythia)
                {
                  if(is_udgjet){
                      FillHist("fh2dJetSignedImpParXYudg",trackV->Pt(),cursImParXY,weight);
                      FillHist("fh2dJetSignedImpParXYSignificanceudg",trackV->Pt(),cursImParXYSig,weight);
                    }
                  FillHist(Form("fh2dJetSignedImpParXY%s",subtype[jetflavour]),jetpt,cursImParXY,weight);
                  FillHist(Form("fh2dJetSignedImpParXYZ%s",subtype[jetflavour]),jetpt,cursImParXYZ,weight);
                  FillHist(Form("fh2dJetSignedImpParXYSignificance%s",subtype[jetflavour]),jetpt,cursImParXYSig,weight);
                  FillHist(Form("fh2dJetSignedImpParXYZSignificance%s",subtype[jetflavour]),jetpt,cursImParXYZSig,weight);
                }
              SJetIpPati a(cursImParXY, weight,kFALSE,kFALSE);
              sImpParXY.push_back(a);
              SJetIpPati b(cursImParXYZ, weight,kFALSE,kFALSE);
              sImpParXYZ.push_back(b);
              SJetIpPati c(cursImParXYSig, weight,kFALSE,kFALSE);
              sImpParXYSig.push_back(c);
              SJetIpPati d(cursImParXYZSig, weight,kFALSE,kFALSE);
              sImpParXYZSig.push_back(d);
            }
        }
      // end of track loop
      std::sort(sImpParXY.begin(),sImpParXY.end(), AliAnalysisTaskHFJetIPQA::mysort);
      std::sort(sImpParXYZ.begin(),sImpParXYZ.end(), AliAnalysisTaskHFJetIPQA::mysort);
      std::sort(sImpParXYSig.begin(),sImpParXYSig.end(), AliAnalysisTaskHFJetIPQA::mysort);
      std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(), AliAnalysisTaskHFJetIPQA::mysort);
      std::reverse(sImpParXY.begin(),sImpParXY.end());
      std::reverse(sImpParXYZ.begin(),sImpParXYZ.end());
      std::reverse(sImpParXYSig.begin(),sImpParXYSig.end());
      std::reverse(sImpParXYZSig.begin(),sImpParXYZSig.end());

      const char * subtype [4] = {"Unidentified","udsg","c","b"};
      const char * subord [3] = {"First","Second","Third"};
      const char * stype [4] = {"fh2dJetSignedImpParXY","fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYZSignificance"};

      for (Int_t ot = 0 ; ot <3 ;++ot){
          if ((int)sImpParXY.size()>ot){
              Double_t params[4] ={ sImpParXY.at(ot).first,sImpParXYZ.at(ot).first, sImpParXYSig.at(ot).first,  sImpParXYZSig.at(ot).first};
              Double_t weights[4] ={sImpParXY.at(ot).second,sImpParXYZ.at(ot).second,sImpParXYSig.at(ot).second,sImpParXYZSig.at(ot).second};
              //non flavour dependent histograms
              for (Int_t ost = 0 ; ost <4 ;++ost){
                  TString hname = Form("%s%s",stype[ost],subord[ot]);
                  FillHist(hname.Data(),jetpt,params[ost],weights[ost] );
                }
              if(fIsPythia){
                  //non flavour dependent histograms
                  for (Int_t ost = 0 ; ost <4 ;++ost){
                          TString hname = Form("%s%s%s",stype[ost],subtype[jetflavour],subord[ot]);
                         // Printf("hname %s",hname.Data());
                          FillHist(hname.Data(),jetpt,params[ost],weights[ost] );
                    }
                }
            }
        }
      sImpParXY.clear();
      sImpParXYZ.clear();
      sImpParXYSig.clear();
      sImpParXYZSig.clear();
    }
  return kTRUE;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits){
  WhyRejected =0;
  Bool_t accept=kTRUE;
  Double_t fMaxVtxZ = 10;
  RejectionBits=000;
  //Physics Selection Cut

  Bool_t isSelected = kFALSE;
  isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(!isSelected) {
      if(accept) WhyRejected=7;
      RejectionBits+=1<<kPhysicsSelection;
      accept=kFALSE;
    }
  // vertex requirements
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!(vertex && vertex->GetNContributors() >0)){
      if(!vertex){
          accept=kFALSE;
          RejectionBits+=1<<kNoVertex;
        }
      else {
          accept=kFALSE;
          RejectionBits+=1<<kNoContributors;
        }
    }else{
      //Test Vertex Chi2/NDF cut
      Double_t x2ndf = event->GetPrimaryVertex()->GetChi2perNDF();

      if(x2ndf > 2. ){
          accept=kFALSE;
          RejectionBits+=1<kVertexChi2NDF;
        }
      const AliVVertex* trkVtx = dynamic_cast<const AliVVertex*>(event->GetPrimaryVertex());
      if(!trkVtx) return kFALSE;

      const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(event->GetPrimaryVertexSPD());
      if(!spdVtx) return kFALSE;

      TString vtxTtl = trkVtx->GetTitle();
      if(!vtxTtl.Contains("VertexerTracks"))
        {
          accept=kFALSE;
          RejectionBits+=1<<kNoVertexTracks;
        }
      if(trkVtx->GetNContributors()<2)	{
          accept=kFALSE;
          RejectionBits+=1<<kTooFewVtxContrib;
        }
      Double_t cov[6] = { 0 };
      spdVtx->GetCovarianceMatrix(cov);
      // Remove all vertexer Z spd
      ///**************Non default

      if(spdVtx->IsFromVertexerZ())  accept=kFALSE;
      ///**************Non default
      /*
         if(spdVtx->IsFromVertexerZ() && (zRes > 0.25)) {
             accept=kFALSE;
             RejectionBits+=1<<kVertexZResolution;
           }
           */
      if((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)) {
          accept=kFALSE;
          RejectionBits+=1<<kDeltaVertexZ;
        }
      if(spdVtx->GetNContributors() <1) {
          accept=kFALSE;
          RejectionBits+=1<<kVertexZContrib;
        }
      if(TMath::Abs(trkVtx->GetZ())>=fMaxVtxZ) {
          RejectionBits+=1<<kZVtxOutFid;
          WhyRejected=6;
          accept=kFALSE;
        }
    }
  //Pileup

  // if(event->IsPileupFromSPD(5, 0.8, 3.0, 2.0, 5.0)) {
  if(event->IsPileupFromSPD(3, 0.8, 3.0, 2.0, 5.0)) {

      if(accept) WhyRejected=1;
      RejectionBits+=1<<kPileupSPD;
      accept=kFALSE;
    }
  /*if(event->IsPileupFromSPDInMultBins()) {
      if(accept) WhyRejected=1;
      RejectionBits+=1<<kPileupSPD;
      accept=kFALSE;
    }*/

  //Special out-of bunch pileup cuts
  // SPD Cluster vs Tracklet plot to estimate pileup effect
  Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);

  if(!(event->GetMultiplicity()))return kFALSE;//!!! TEST

  Int_t nTracklets = event->GetMultiplicity()->GetNumberOfTracklets();
  if(nClustersLayer0 + nClustersLayer1 > 65 + 4 * nTracklets){
      accept=kFALSE;
      RejectionBits+=1<<kSPDClusterCut;
    }
  if(fUtils->IsPileUpMV(event)){
      accept=kFALSE;
      RejectionBits+=1<<kMVPileup;
    }



  return accept;
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::SetUseMonteCarloWeighingLinus(TH1F *Pi0, TH1F *Eta, TH1F *EtaP, TH1F *Rho, TH1F *Phi, TH1F *Omega, TH1F *K0s, TH1F *Lambda, TH1F *ChargedPi, TH1F *ChargedKaon, TH1F *Proton, TH1F *D0, TH1F *DPlus, TH1F *DStarPlus, TH1F *DSPlus, TH1F *LambdaC, TH1F *BPlus, TH1F *B0, TH1F *LambdaB, TH1F *BStarPlus)
{
  for(Int_t i =1 ; i< Pi0->GetNbinsX()+1;++i){
      fBackgroundFactorLinus[bIdxPi0][i-1] =Pi0->GetBinContent(i);
      fBackgroundFactorLinus[bIdxEta][i-1] =Eta->GetBinContent(i);
      fBackgroundFactorLinus[bIdxEtaPrime][i-1] =EtaP->GetBinContent(i);
      fBackgroundFactorLinus[bIdxRho][i-1] =Rho->GetBinContent(i);
      fBackgroundFactorLinus[bIdxPhi][i-1] =Phi->GetBinContent(i);
      fBackgroundFactorLinus[bIdxOmega][i-1] =Omega->GetBinContent(i);
      fBackgroundFactorLinus[bIdxK0s][i-1] =K0s->GetBinContent(i);
      fBackgroundFactorLinus[bIdxLambda][i-1] =Lambda->GetBinContent(i);
      fBackgroundFactorLinus[bIdxPi][i-1] =ChargedPi->GetBinContent(i);
      fBackgroundFactorLinus[bIdxKaon][i-1] =ChargedKaon->GetBinContent(i);
      fBackgroundFactorLinus[bIdxProton][i-1] =Proton->GetBinContent(i);
      fBackgroundFactorLinus[bIdxD0][i-1] =D0->GetBinContent(i);
      fBackgroundFactorLinus[bIdxDPlus][i-1] =DPlus->GetBinContent(i);
      fBackgroundFactorLinus[bIdxDStarPlus][i-1] =DStarPlus->GetBinContent(i);
      fBackgroundFactorLinus[bIdxDSPlus][i-1] =DSPlus->GetBinContent(i);
      fBackgroundFactorLinus[bIdxLambdaC][i-1] =LambdaC->GetBinContent(i);
      fBackgroundFactorLinus[bIdxBPlus][i-1] =BPlus->GetBinContent(i);
      fBackgroundFactorLinus[bIdxB0][i-1] =B0->GetBinContent(i);
      fBackgroundFactorLinus[bIdxLambdaB][i-1] =LambdaB->GetBinContent(i);
      fBackgroundFactorLinus[bIdxBStarPlus][i-1] =BStarPlus->GetBinContent(i);
    }
  return;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::SetResFunction(TF1 *f, Int_t i, Int_t j){
int counter =0;
  for(double ip = 0 ; ip<100;ip = ip+1e-3 ){
      fResolutionFunction[i*5+j].SetPoint(counter,ip,f->Eval(ip));
      counter +=1;
    }
  return kTRUE;
}
// ######################################################################################## JEt Probability Function
Double_t AliAnalysisTaskHFJetIPQA::CalculatePSTrack(Double_t sign, Double_t significance ,Double_t trackPt,Int_t trclass)
{
  Double_t retval = 0;
  //switch resolution function based on track pt;
  Int_t ptbin=0;
  Double_t pt = trackPt;
  if(pt>1.5) ptbin=1;
  if(TMath::Abs(significance) >100) significance =100; //Limit to function definition range
  retval = sign * ((fResolutionFunction[trclass+ptbin*5])).Eval(TMath::Abs(significance));
  return retval;
}
// ######################################################################################## JEt Probability Function
Double_t AliAnalysisTaskHFJetIPQA::CalculateJetProb(AliEmcalJet *jet)
{
  if(!jet) return -9999;
  Double_t retval = -1;
  //Loop over all tracks calculate P(s) for all accepted later add looser cuts also
  Int_t ntracks = (Int_t)jet->GetNumberOfTracks();
  Double_t prodPS = 1;
  Double_t curps=-1;
  AliJetContainer * jetconrec = 0x0;
  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  Int_t jcounter =0;
  for(Int_t itrack = 0; itrack < ntracks; ++itrack)
    {
      AliVTrack * trackV = (((AliVTrack*)((jetconrec->GetParticleContainer())->GetParticle(jet->TrackAt(itrack)))));
      //class selection
      Bool_t isAccepted=kFALSE;
      Int_t track_class=0;
      if(IsTrackAccepted(trackV,6)) {isAccepted=kTRUE;track_class=0;}
      else if(IsTrackAccepted(trackV,5)) {isAccepted=kTRUE;track_class=1;}
      else if(IsTrackAccepted(trackV,4)) {isAccepted=kTRUE;track_class=2;}
      else if(IsTrackAccepted(trackV,-5)) {isAccepted=kTRUE;track_class=3;}
     // else if(IsTrackAccepted(trackV,-4)) {isAccepted=kTRUE;track_class=4;}
      if(isAccepted){
          Double_t dca[2] = {0,0};
          Double_t cov[3] = {0,0,0};
          Double_t dcajettrack =999;
          Double_t lineardecay =999;
          Double_t sign =1;
          CalculateJetSignedTrackImpactParameter(trackV,jet,dca,cov,sign,dcajettrack,lineardecay);
          curps =CalculatePSTrack(sign,GetValImpactParameter(kXYSig,dca,cov),trackV->Pt() ,track_class);
          prodPS*=(curps>=0 ? curps/2. : 1+curps/2.);
          jcounter++;
        }
    }
  Double_t sumPS =0;
  bool chan=false;
  for(Int_t j=0;j<jcounter;++j){
      double val = TMath::Power(-1 * TMath::Log(prodPS),j)/TMath::Factorial(j);;
      sumPS += val;
      chan=true;
    }
  if(!chan)return -1;
  retval =sumPS *prodPS;
  return retval;
}
// ######################################################################################## Event Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsEventSelected()	{
  if(fIsPythia){
      AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if(!mcH ){
          AliError("No MC Event Handler available");
          return kFALSE;
        }


      if(!mcH->InitOk()) return kFALSE;
      if(!mcH->TreeK()) return kFALSE ;
      if(!mcH->TreeTR()) return kFALSE;

    }
  AliAODEvent* aev = NULL;
  Int_t WhyRejected =0;
  ULong_t RejectionBits=0;



  if(!fESD)
    {
      aev = dynamic_cast<AliAODEvent*>(InputEvent());
      if(aev && aev->GetPrimaryVertex() && aev->GetPrimaryVertex()->GetNContributors()>0){
          FillHist("fh1dVertexZ",aev->GetPrimaryVertex()->GetZ(),1);
          Double_t vtxx = aev->GetPrimaryVertex()->GetX();
          Double_t vtxy = aev->GetPrimaryVertex()->GetY();
          FillHist("fh1dVertexR",vtxx,vtxy,1);

        }else return kFALSE;


      if(!(IsSelected(aev,WhyRejected,RejectionBits)))
        {
          IncHist("fh1dEventRejectionRDHFCuts",2);
          if(RejectionBits&(1<<kPhysicsSelection))     IncHist("fh1dEventRejectionRDHFCuts",3);
          else if(RejectionBits&(1<<kNoVertex))         IncHist("fh1dEventRejectionRDHFCuts",6);
          else if(RejectionBits&(1<<kNoVertexTracks))   IncHist("fh1dEventRejectionRDHFCuts",11);
          else if(RejectionBits&(1<<kNoContributors))   IncHist("fh1dEventRejectionRDHFCuts",12);
          else if(RejectionBits&(1<<kTooFewVtxContrib))  IncHist("fh1dEventRejectionRDHFCuts",7);
          else if(RejectionBits&(1<<kVertexZContrib))  IncHist("fh1dEventRejectionRDHFCuts",15);
          else if(RejectionBits&(1<<kVertexZResolution))  IncHist("fh1dEventRejectionRDHFCuts",14);
          else if(RejectionBits&(1<<kDeltaVertexZ))     IncHist("fh1dEventRejectionRDHFCuts",13);
          else if(RejectionBits&(1<<kZVtxOutFid))       IncHist("fh1dEventRejectionRDHFCuts",5);
          else if(RejectionBits&(1<<kOutsideCentrality))  IncHist("fh1dEventRejectionRDHFCuts",4);
          else if(RejectionBits&(1<<kNotSelTrigger))    IncHist("fh1dEventRejectionRDHFCuts",8);
          else if(RejectionBits&(1<<kSPDClusterCut))     IncHist("fh1dEventRejectionRDHFCuts",9);
          else if(RejectionBits&(1<<kMVPileup))          IncHist("fh1dEventRejectionRDHFCuts",10);
          return kFALSE;
        }else {
          IncHist("fh1dEventRejectionRDHFCuts",1);
          FillHist("fh1dVertexZAccepted",aev->GetPrimaryVertex()->GetZ(),1);
          Double_t vtxx = aev->GetPrimaryVertex()->GetX();
          Double_t vtxy = aev->GetPrimaryVertex()->GetY();
          FillHist("fh1dVertexRAccepted",vtxx,vtxy,1);
          FillHist("fh2dVertexChi2NDFNESDTracks",aev->GetPrimaryVertex()->GetChi2perNDF(),aev->GetNumberOfESDTracks(),1);
          return kTRUE;
        }
    }
  AliESDEvent* eev = NULL;
  if(fESD)
    {
      eev = dynamic_cast<AliESDEvent*>(InputEvent());
      if(eev && eev->GetPrimaryVertex() && eev->GetPrimaryVertex()->GetNContributors()>0){
          FillHist("fh1dVertexZ",eev->GetPrimaryVertex()->GetZ(),1);
          Double_t vtxx = eev->GetPrimaryVertex()->GetX();
          Double_t vtxy = eev->GetPrimaryVertex()->GetY();
          FillHist("fh1dVertexR",vtxx,vtxy,1);
        }else {
          return kFALSE;
        }

      if(!(IsSelected(eev,WhyRejected,RejectionBits)))
        {
          IncHist("fh1dEventRejectionRDHFCuts",2);
          if(RejectionBits&(1<<kPhysicsSelection))     IncHist("fh1dEventRejectionRDHFCuts",3);
          else if(RejectionBits&(1<<kNoVertex))         IncHist("fh1dEventRejectionRDHFCuts",6);
          else if(RejectionBits&(1<<kNoVertexTracks))   IncHist("fh1dEventRejectionRDHFCuts",11);
          else if(RejectionBits&(1<<kNoContributors))   IncHist("fh1dEventRejectionRDHFCuts",12);
          else if(RejectionBits&(1<<kTooFewVtxContrib))  IncHist("fh1dEventRejectionRDHFCuts",7);
          else if(RejectionBits&(1<<kVertexZContrib))  IncHist("fh1dEventRejectionRDHFCuts",15);
          else if(RejectionBits&(1<<kVertexZResolution))  IncHist("fh1dEventRejectionRDHFCuts",14);
          else if(RejectionBits&(1<<kDeltaVertexZ))     IncHist("fh1dEventRejectionRDHFCuts",13);
          else if(RejectionBits&(1<<kZVtxOutFid))       IncHist("fh1dEventRejectionRDHFCuts",5);
          else if(RejectionBits&(1<<kOutsideCentrality))  IncHist("fh1dEventRejectionRDHFCuts",4);
          else if(RejectionBits&(1<<kNotSelTrigger))    IncHist("fh1dEventRejectionRDHFCuts",8);
          else if(RejectionBits&(1<<kSPDClusterCut))     IncHist("fh1dEventRejectionRDHFCuts",9);
          else if(RejectionBits&(1<<kMVPileup))          IncHist("fh1dEventRejectionRDHFCuts",10);
          return kFALSE;
        }else {
          Double_t vtxx = eev->GetPrimaryVertex()->GetX();
          Double_t vtxy = eev->GetPrimaryVertex()->GetY();
          Double_t vtxz = eev->GetPrimaryVertex()->GetZ();
          FillHist("fh1dVertexXvsMultiplicity",vtxx,eev->GetNumberOfTracks(),1);
          FillHist("fh1dVertexYvsMultiplicity",vtxy,eev->GetNumberOfTracks(),1);
          FillHist("fh1dVertexZvsMultiplicity",vtxz,eev->GetNumberOfTracks(),1);
          IncHist("fh1dEventRejectionRDHFCuts",1);
          FillHist("fh1dVertexZAccepted",eev->GetPrimaryVertex()->GetZ(),1);
          FillHist("fh1dVertexRAccepted",vtxx,vtxy,1);
          FillHist("fh2dVertexChi2NDFNESDTracks",eev->GetPrimaryVertex()->GetChi2perNDF(),eev->GetNumberOfTracks(),1);
          return kTRUE;
        }
    }
  return kFALSE;
}
// ######################################################################################## Init histograms
void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){
  // AliAnalysisTaskEmcal::UserCreateOutputObjects();
  Double_t xomega[10] ={0,1.00993,1.40318,1.7428,2.04667,2.35055,2.77061,3.39623,4.40616,1000};
  Double_t yomega[10] ={2.74011,2.74011,3.16949,3.37288,4.02825,4.38983,4.23164,4.75141,3.62147,3.62147};
  Double_t xxi[18] ={0,0.723932,0.849057,0.947368,1.04568,1.14399,1.25124,1.35849,1.43893,1.60874,1.80536,2.04667,2.41311,
                     2.87786,3.50348,4.39722,5.45184,1000};
  Double_t yxi[18] ={1.20339,1.20339,1.45198,1.54237,1.76836,1.81356,1.85876,1.97175,2.10734,2.10734,2.15254,2.19774,2.22034,
                     2.19774,2.12994,1.83616,1.36158,1.36158};
  Double_t xK0s[24]= {0,0.0628272,0.162304,0.26178,0.356021,0.465969,0.570681,0.675393,0.743455,0.863874,0.947644,1.10471,1.29843,
                      1.50785,1.70681,1.91099,2.19895,2.60733,3.01571,3.40838,3.82199,4.51832,5.49215,1000};
  Double_t yK0s[24]= {1.31496,1.31496,1.19685,1.11024,1.16535,1.11811,1.07874,1.05512,
                      1.01575,0.976378,0.92126,0.889764,0.858268,0.811024,0.84252,0.858268,
                      0.811024,0.811024,0.779528,0.787402,0.795276,0.811024,0.88189,0.88189};
  Double_t xPhi[24] ={0,0.456294,0.54021,0.645105,0.76049,0.833916,0.938811,1.03322,1.15385,
                      1.23776,1.34266,1.46853,1.5472,1.6521,1.75175,1.86189,1.96154,2.09266,
                      2.27622,2.51748,2.71678,2.91084,3.2465,1000};
  Double_t yPhi[24] ={0.699725,0.699725,0.699725,0.721763,0.661157,0.683196,0.650138,
                      0.639118,0.639118,0.62259,0.61157,0.62259,0.61708,0.61708,0.606061,
                      0.644628,0.694215,0.699725,0.721763,0.787879,0.77686,0.809917,0.870523,0.870523};
  //This is MC /Data not Data over MC
  fGraphOmega = new TGraph(10,xomega,yomega);
  fGraphXi = new TGraph(18,xxi,yxi);
  fK0Star = new TGraph(24,xK0s,yK0s);
  fPhi = new TGraph(24,xPhi,yPhi);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
      AliVEventHandler *evhand = mgr->GetInputEventHandler();
      if (evhand) {
          if (evhand->InheritsFrom("AliESDInputHandler")) {
              fIsEsd = kTRUE;
              fESD=kTRUE;
            }
          else {
              fIsEsd = kFALSE;
              fESD=kFALSE;
            }
        }
      else  AliError("Event handler not found!");
    }
  else AliError("Analysis manager not found!");

  const Int_t nBins3dSignificance =500;
  const Int_t nBins3d =250;

  Double_t lowIPxy =-1.;
  Double_t highIPxy =1.;
  if (!fOutput2) fOutput2 = new TList ();
  fOutput2->SetOwner(kTRUE);
  //Make Graphs
  const Int_t gfProtonN = 9;
  Double_t gfProtonX [gfProtonN] ={0,0.534483,1.29741,2.21552,3.0819,3.92241,4.5819,5.39655,1000};
  Double_t gfProtonY [gfProtonN] ={0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964};

  const Int_t gfAntiProtonN = 18;
  Double_t gfAntiProtonX [gfAntiProtonN] ={0,0.806034,0.922414,1.09052,1.28448,1.5431,1.73707,1.89224,2.17672,2.43534,2.74569,3.06897,
                                           3.52155,3.88362,4.38793,5.03448,5.38362, 1000};
  Double_t gfAntiProtonY [gfAntiProtonN] ={0.922892,0.922892,	0.930723,	0.939157,0.94397,0.95241,0.956627,0.959639,0.964458,
                                           0.966867,0.971084,0.974096,0.978313,0.98012,0.983735,0.986747,0.989157,0.989157};


  const Int_t gfAntiLambdaN = 34;
  Double_t gfAntiLambdaX [gfAntiLambdaN] ={0.,0.55555,0.64646,0.75757,	0.84848,0.94949,1.06061,1.15152,1.24242,1.35354,1.44444,
                                           1.54545,1.66667,1.75758,1.84848,1.9596,2.09091,2.30303,2.50505,2.68687,2.90909,3.11111,
                                           3.31313,3.51515,3.69697,3.89899,4.20202,4.66667,5.21212,5.74747,6.50505,7.51515,9.0101,1000};
  Double_t gfAntiLambdaY [gfAntiLambdaN] = {0.864925,0.864925,0.895896,0.908209,0.915672,0.921269,0.926866,0.931343,0.935821,0.938806,0.942164,
                                            0.945149,0.947761,0.95,0.952612,0.954478,0.957836,0.960821,0.96306,0.965672,0.968657,0.970149,
                                            0.972015,0.973507,0.975,0.976493,0.978358,0.981343,0.983955,0.986194,0.988433,0.991045,0.991045,0.991045};

  const Int_t gfLambdaN = 2;
  Double_t gfLambdaX [gfLambdaN] =	{0.,1000};
  Double_t gfLambdaY [gfLambdaN] = {0.991045,0.991045};
  const Int_t gfKMinusN =13 ;
  Double_t gfKMinusX [gfKMinusN] =	{0,0.54741,0.74137,1.03879,1.36207,1.96983,2.52586,3.0819,3.67672,4.19397,5.03448,5.44828,1000};
  Double_t gfKMinusY [gfKMinusN] = {0,0.979518,0.983133,0.987349,0.989759,0.992169,0.993976,0.996386,0.995783,0.998193,0.99759,1,1000};

  fGeant3FlukaProton 	= new TGraph(gfProtonN,gfProtonX,gfProtonY);
  fGeant3FlukaAntiProton 	= new TGraph(gfAntiProtonN,gfAntiProtonX,gfAntiProtonY);
  fGeant3FlukaLambda   	= new TGraph(gfLambdaN,gfLambdaX,gfLambdaY);
  fGeant3FlukaAntiLambda  = new TGraph(gfAntiLambdaN,gfAntiLambdaX,gfAntiLambdaY);
  fGeant3FlukaKMinus 		= new TGraph(gfKMinusN,gfKMinusX,gfKMinusY);


  //ADD HISTOGRAMS
  TH1D * h = (TH1D*)AddHistogramm("fh1dEventRejectionRDHFCuts","fh1dEventRejectionRDHFCuts;reason;count",15,0,15);
  h->GetXaxis()->SetBinLabel(1,"Accepted");
  h->GetXaxis()->SetBinLabel(2,"Rejected");
  h->GetXaxis()->SetBinLabel(3,"DueToPhysicsSelection");
  h->GetXaxis()->SetBinLabel(4,"DueCentralitySelection");
  h->GetXaxis()->SetBinLabel(5,"ZVertexOutsideFiducialRegion");
  h->GetXaxis()->SetBinLabel(6,"IsEventRejectedDueToNotRecoVertex");
  h->GetXaxis()->SetBinLabel(7,"IsEventRejectedDueToVertexContributors");
  h->GetXaxis()->SetBinLabel(8,"DueToTrigger");
  h->GetXaxis()->SetBinLabel(9,"DueToSPDTrackletClusterCut");
  h->GetXaxis()->SetBinLabel(10,"DueToMVPileup");
  h->GetXaxis()->SetBinLabel(11,"NoVertexTracks");
  h->GetXaxis()->SetBinLabel(12,"NoContributorsVertexTracks");
  h->GetXaxis()->SetBinLabel(13,"DeltaVertexZSPDTracks");
  h->GetXaxis()->SetBinLabel(14,"ZVertexResolution");
  h->GetXaxis()->SetBinLabel(15,"VertexZContributors");



  AddHistogramm("fh1dVertexXvsMultiplicity",";cm;# ESD Tracks",1000,-1,1,200,0,200);
  AddHistogramm("fh1dVertexYvsMultiplicity",";cm;# ESD Tracks",1000,-1,1,200,0,200);
  AddHistogramm("fh1dVertexZvsMultiplicity",";cm;# ESD Tracks",1000,-1,1,200,0,200);

  //Vertex Z before and after
  AddHistogramm("fh1dVertexZ","Vertex Z before Event selection;primary vertex z (cm);count",500,-30,30);
  AddHistogramm("fh1dVertexZAccepted","Vertex Z after Event selection;primary vertex z (cm);count",500,-30,30);

  AddHistogramm("fh1dVertexR","Vertex R before Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);
  AddHistogramm("fh1dVertexRAccepted","Vertex R after Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);

  // Vertex Chi2/NDF vs TPC track multiplicity(ESD tracks)
  AddHistogramm("fh2dVertexChi2NDFNESDTracks","Vertex Chi2/NDF vs # tracks ESD;vertex #chi^{2}/NDF;# tracks esd",200,0,10,500,0,500);
  // AOD tracks accepted
  AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
  TH1D * h1 = GetHist1D("fh1dTracksAccepeted");
  h1->GetXaxis()->SetBinLabel(1,"total");
  h1->GetXaxis()->SetBinLabel(2,"accepted");
  h1->GetXaxis()->SetBinLabel(3,"rejected");
  // Tracks impact parameter histograms
  AddHistogramm("fh2dTracksImpParXY","radial imp. parameter ;impact parameter xy (cm);a.u.",5000,lowIPxy,highIPxy,500,0,100.);
  AddHistogramm("fh2dTracksImpParZ","z imp. parameter ;impact parameter xy (cm);a.u.",5000,lowIPxy,highIPxy,500,0,10.);
  AddHistogramm("fh2dTracksImpParXYSignificance","radial imp. parameter sig;impact parameter xy (cm);a.u.",5000,-30,30,500,0,100.);
  AddHistogramm("fh2dTracksImpParZSignificance","z imp. parameter ;impact parameter xy (cm);a.u.",5000,-30,30,500,0,100.);
  AddHistogramm("fh1dTracksImpParXY","2d imp. parameter ;impact parameter 2d (cm);a.u.",5000,-1,1.);
  AddHistogramm("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",nBins3d,0,1.);
  AddHistogramm("fh1dTracksImpParXYSignificance","radial imp. parameter ;impact parameter xy significance;a.u.",5000,-30,30);
  AddHistogramm ("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",nBins3dSignificance/2,0.,100.);
  AddHistogramm("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",200,-0.5,0.5,200,0.,TMath::TwoPi());
  AddHistogramm("fh2dAcceptedTracksEtaPhi","accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi());

  for(Int_t i = 0; i<6;++i){
      fh2dAcceptedTracksEtaPhiPerLayer[i] = new TH2D(Form("fh2dAcceptedTracksEtaPhiPerLayer_%i",i),"accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi());
    }

  AddHistogramm("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250);
  AddHistogramm("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250);

  if (fIsPythia){
      AddHistogramm("fh2dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",5000,lowIPxy,highIPxy,500,0,100);
      AddHistogramm("fh1dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",5000,-1,1.);
      AddHistogramm("fh1dTracksImpParXYZ_McCorr","3d imp. parameter (after correction);impact parameter 3d (cm);a.u.",5000,0,10.);
      AddHistogramm("fh1dTracksImpParXYSignificance_McCorr","radial imp. parameter (after correction);impact parameter xy significance;a.u.",5000,-30,30.);
      AddHistogramm("fh1dTracksImpParXYZSignificance_McCorr","3d imp. parameter (after correction);impact parameter 3d significance;a.u.",5000,0.,100.);
      AddHistogramm("fh1dJetGenPt","generator level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250);
      AddHistogramm("fh2dJetSignedImpParXYudg","fh2dJetSignedImpParXYudg;pt (GeV/c); count",1000,0,100,2000,-1,1);
      AddHistogramm("fh2dJetSignedImpParXYSignificanceudg","fh2dJetSignedImpParXYSignificanceudg;pt (GeV/c); count",1000,0,100,5000,-100,100);
      AddHistogramm("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtUnidentifiedAccepted","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtudsgAccepted","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtcAccepted","detector level jets;pt (GeV/c); count",500,0,250);
      AddHistogramm("fh1dJetRecPtbAccepted","detector level jets;pt (GeV/c); count",500,0,250);
    }

  const char * flavour[5]  = {"Unidentified","udsg","c","b",""};

  for (Int_t i = 0 ;i < 5;++i){
      if (!fIsPythia && (i<4)) continue;
      AddHistogramm(Form("JetProbability_%s",flavour[i]),Form("JetProbability_%s;jept",flavour[i]),500,0,250,500,0,1);
    }
  const char * base = "fh2dJetSignedImpPar";
  const char * dim[2]  = {"XY","XYZ"};
  const char * typ[2]  = {"","Significance"};
  const char * ordpar [4] = {"","First","Second","Third"};
  const char * special [1] = {"",/*"McCorr"*/};

  Int_t ptbins = 500;
  Double_t ptlow = 0;
  Double_t pthigh = 250;

  Int_t ipbins = 1000;
  Double_t iplow = -1;
  Double_t iphigh = 1;

  for (Int_t id = 0;id<2;++id)
    for (Int_t ifl = 0;ifl<5;++ifl)
      for (Int_t io = 0;io<4;++io)
        for (Int_t is = 0;is<1;++is)
          for (Int_t it = 0;it<2;++it){
              if(it==1) {
                  iplow=-30;
                  iphigh=30;
                  if(io==0 && ifl==4) ipbins =1000;
                  else  ipbins =1000;
                }else {
                  iplow=-1;
                  iphigh=1;
                  ipbins =1000;

                }

              if((fIsPythia||(!fIsPythia && ifl==4))/*&&(!fCorrrectionSamplingMode)*/)  AddHistogramm(Form("%s%s%s%s%s%s",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                      Form("%s%s%s%s%s%s;;",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                      ptbins,ptlow,pthigh,ipbins,iplow,iphigh);
            }
  // Notify();
  PostData(1, fOutput2); // Post data for ALL output slots > 0 here.
}
// ######################################################################################## Calculate impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliAODTrack * track,Double_t *impar, Double_t * cov)
{
  AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODVertex *vtxAOD = aev->GetPrimaryVertex();
  AliAODVertex *vtxAODNew=vtxAOD;
  if(!vtxAOD) return kFALSE;
  TString title=vtxAOD->GetTitle();
  if(!title.Contains("VertexerTracks")) return kFALSE;
  AliESDVertex *vtxESDNew =0x0;
  Bool_t recalculate = kFALSE;
  if( vtxAOD->GetNContributors() < 30){
      recalculate=kTRUE;
      AliVertexerTracks *vertexer = new AliVertexerTracks(aev->GetMagneticField());
      Int_t ndg = 1;
      vertexer->SetITSMode();
      vertexer->SetMinClusters(3);
      vertexer->SetConstraintOff();
      if(title.Contains("WithConstraint")) {
          Float_t diamondcovxy[3];
          aev->GetDiamondCovXY(diamondcovxy);
          Double_t pos[3]={aev->GetDiamondX(),aev->GetDiamondY(),0.};
          Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
          AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
          vertexer->SetVtxStart(diamond);
          delete diamond; diamond=NULL;
        }
      Int_t skipped[1] = {-1};
      Int_t id = (Int_t)track->GetID();
      if(id<0) return kFALSE;
      skipped[0] = id;
      vertexer->SetSkipTracks(1,skipped);
      vtxESDNew = vertexer->FindPrimaryVertex(aev);
      delete vertexer; vertexer=NULL;
      if(!vtxESDNew) return kFALSE;
      if(vtxESDNew->GetNContributors()<=0) {
          delete vtxESDNew; vtxESDNew=NULL;
          return kFALSE;
        }
      // convert to AliAODVertex
      Double_t pos[3],cova[6],chi2perNDF;
      vtxESDNew->GetXYZ(pos); // position
      vtxESDNew->GetCovMatrix(cova); //covariance matrix
      chi2perNDF = vtxESDNew->GetChi2toNDF();
      delete vtxESDNew; vtxESDNew=NULL;
      vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
    }
  // Calculate Impact Parameters
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  if(etp.PropagateToDCA(vtxAODNew,aev->GetMagneticField(),3.,impar,cov))
    {

      if(recalculate) delete vtxAODNew;
      return kTRUE;
    }
  else{
      if(recalculate)
        delete vtxAODNew;
      return kFALSE;

    }
}
// ######################################################################################## Calculate impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliESDtrack * track,Double_t *impar, Double_t * cov,Bool_t useTRUEvtx)
{
  //Printf("eev->GetMagneticField() %f",eev->GetMagneticField());
  const Double_t kBeampiperadius=3.;
  AliVVertex *vtxESDSkip = NULL;
  AliESDEvent * eev = (AliESDEvent*)InputEvent();
  if(!eev){
      AliDebug(1, "No Input event available\n");
      return kFALSE;
    }
  TString type = track->IsA()->GetName();
  Double_t    dcaD[2]={-999.,-999.},
      covD[3]={-999.,-999.,-999.};
  Bool_t isRecalcVertex(kFALSE);
  //case of ESD tracks
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(InputEvent());
  if(!esdevent) {
      AliDebug(1, "No esd event available\n");
      return kFALSE;
    }
  vtxESDSkip = (	 AliVVertex *)esdevent->GetPrimaryVertex();
  if(!vtxESDSkip) return kFALSE;
  //case ESD track: take copy constructor
  const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
  if(tmptrack ){
      if( vtxESDSkip->GetNContributors() < 30&& !useTRUEvtx){ // if vertex contributor is smaller than 30, recalculate the primary vertex
          AliVertexerTracks vertexer(eev->GetMagneticField());
          vertexer.SetITSMode();
          vertexer.SetMinClusters(3);
          Int_t skipped[1] = {track->GetID()};
          vertexer.SetSkipTracks(1,skipped);
          vertexer.SetConstraintOn();
          Float_t diamondcovxy[3];
          esdevent->GetDiamondCovXY(diamondcovxy);
          Double_t pos[3]={eev->GetDiamondX(),eev->GetDiamondY(),0};
          Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
          AliESDVertex diamond(pos,cov,1.,1);
          vertexer.SetVtxStart(&diamond);
          vtxESDSkip = vertexer.FindPrimaryVertex(eev);
          if(vtxESDSkip->GetNContributors()<1) {
              delete vtxESDSkip; vtxESDSkip=NULL;
            }
          isRecalcVertex = kTRUE;
        }
      if(vtxESDSkip){
          AliESDtrack esdtrack(*tmptrack);
          if(useTRUEvtx) {
              vtxESDSkip = (AliVVertex * )MCEvent()->GetPrimaryVertex();
            }
          if(esdtrack.PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, impar, cov)){
              // cov[0] = esdtrack.GetCovariance()[0];

              if(esdtrack.GetSigmaY2()<0. || esdtrack.GetSigmaZ2()<0.) {
                  if(isRecalcVertex) delete vtxESDSkip;
                  return kFALSE;
                  // this is insipired by the AliITStrackV2::Invariant() checks
                }
              if(isRecalcVertex) delete vtxESDSkip;
              if(abs(impar[1])>2.)return kFALSE;
              return kTRUE;
            }
          else  delete vtxESDSkip;
          if(isRecalcVertex) delete vtxESDSkip;
          return kFALSE;
        }
    }
  return kFALSE;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliVTrack *track, Double_t *impar, Double_t *cov){

  if (fESD){
      AliESDtrack * tr = (AliESDtrack*)track;
      return CalculateTrackImpactParameter(tr,impar, cov,kFALSE);
    }
  return CalculateTrackImpactParameter((AliAODTrack*)track,impar, cov);
}

// ######################################################################################## Calculate impact parameters based on MC event vertex and MC particle information (no special mass treatment)
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameterTruth(AliAODTrack * track,Double_t *impar, Double_t * cov)
{
  AliAODMCParticle *pMC = 0x0;
  AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODMCHeader* mcheader = dynamic_cast<AliAODMCHeader*>(aev->FindListObject(AliAODMCHeader::StdBranchName()));
  if (!mcheader) return kFALSE;

  TClonesArray * fMCparticles = dynamic_cast<TClonesArray*>(aev->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fMCparticles) return kFALSE;
  if(track->GetLabel()>-1)
    pMC = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(track->GetLabel()));
  if (!pMC) return kFALSE;

  Double_t pos[3]={0,0,0};
  mcheader->GetVertex(pos);

  Double_t cova[6]={0,0,0,0,0,0};
  Double_t chi2perNDF =0;
  AliAODVertex *vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
  Double_t xpart[3] = {0,0,0};
  pMC->XvYvZv(xpart);
  Double_t ppart[3] = {0,0,0};
  pMC->PxPyPz(ppart);
  Double_t cv[21] ;
  for (Int_t i=0;i<21;++i)cv[i] =0.;
  AliExternalTrackParam trackparam(xpart,ppart,cv,(TMath::Sign((Short_t)1,(Short_t)pMC->Charge())));

  if(trackparam.PropagateToDCA(vtxAODNew,aev->GetMagneticField(),3.,impar,cov))
    {
      delete vtxAODNew;
      return kTRUE;
    }
  else{
      delete vtxAODNew;
      return kFALSE;

    }

  return kFALSE;
}
// ######################################################################################## Calculate impact parameters based on MC event vertex and MC particle information (no special mass treatment)
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameterTruth(AliESDtrack * track,Double_t *impar, Double_t * cov)
{
  AliMCParticle *pMC = 0x0;
  AliESDEvent* eev = dynamic_cast<AliESDEvent*>(InputEvent());
  AliMCEvent * MCparticles = dynamic_cast<AliMCEvent*>(MCEvent());
  if (!MCparticles) return kFALSE;
  pMC = dynamic_cast<AliMCParticle*>(MCEvent()->GetTrack(abs(track->GetLabel())));
  if (!pMC) return kFALSE;
  Double_t pos[3]={0,0,0};
  MCparticles->GetPrimaryVertex()->GetXYZ(pos);
  Double_t cova[6]={0,0,0,0,0,0};
  Double_t chi2perNDF =0;
  const AliVVertex *vtxAODNew =  MCparticles->GetPrimaryVertex();
  Double_t xpart[3] = {0,0,0};
  pMC->XvYvZv(xpart);
  Double_t ppart[3] = {0,0,0};
  pMC->PxPyPz(ppart);
  Double_t cv[21] ;
  for (Int_t i=0;i<21;++i)cv[i] =0.;
  AliExternalTrackParam trackparam(xpart,ppart,cv,(TMath::Sign((Short_t)1,(Short_t)pMC->Charge())));
  if(trackparam.PropagateToDCA(vtxAODNew,eev->GetMagneticField(),3.,impar,cov))
    {
      return kTRUE;
    }
  else{
      return kFALSE;

    }
  return kFALSE;
}
// ######################################################################################## Calculate signed  impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliAODTrack * track,AliEmcalJet * jet ,Double_t *impar, Double_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength){
  AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODVertex *vtxAOD = aev->GetPrimaryVertex();
  if(!vtxAOD) return kFALSE;
  Double_t pos[3],cova[6],chi2perNDF;
  vtxAOD->GetXYZ(pos); // position
  vtxAOD->GetCovMatrix(cova); //covariance matrix

  // Calculate Impact Parameters
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  if(etp.PropagateToDCA(vtxAOD,aev->GetMagneticField(),3.,impar,cov))
    {
      //Calculate Sign
      Double_t posdcatrack[3]= {0.,0.,0.};
      etp.GetXYZ(posdcatrack);
      Double_t ipvector3[3] = { posdcatrack[0] - pos[0], posdcatrack[1] - pos[1], posdcatrack[2] - pos[2] };
      sign =TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz() );
      // Calculate decay legnth and track jet DCA against new vertex
      Double_t bpos[3] = { 0,0,0 };
      vtxAOD->GetXYZ(bpos);
      Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };
      Double_t bcv[21] = { 0 };
      AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
      Double_t xa = 0., xb = 0.;
      bjetparam.GetDCA(&etp, aev->GetMagneticField(), xa, xb);
      Double_t xyz[3] = { 0., 0., 0. };
      Double_t xyzb[3] = { 0., 0., 0. };
      bjetparam.GetXYZAt(xa, aev->GetMagneticField(), xyz);
      etp.GetXYZAt(xb, aev->GetMagneticField(), xyzb);
      Double_t  bdecaylength =
          TMath::Sqrt(
            (bpos[0] - xyz[0]) * (bpos[0] - xyz[0]) +
          (bpos[1] - xyz[1]) * (bpos[1] - xyz[1]) +
          (bpos[2] - xyz[2]) * (bpos[2] - xyz[2]));
      dcajetrack =
          TMath::Sqrt(
            (xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
          (xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
          (xyzb[2] - xyz[2]) * (xyzb[2] - xyz[2]));
      if(bdecaylength>0) lineardecaylength=bdecaylength;
      //delete vtxAOD;
      return kTRUE;
    }
  else{
      //delete vtxAODNew;
      return kFALSE;

    }
}

// ######################################################################################## Calculate signed  impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliESDtrack * track,AliEmcalJet * jet ,Double_t *impar, Double_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength){
  AliESDEvent* eev = dynamic_cast<AliESDEvent*>(InputEvent());
  //
  // Copied from AliHFEextraCuts::GetHFEImpactParameters  impact parameter (with recalculated primary vertex)
  //
  if(!eev){
      AliDebug(1, "No Input event available\n");
      return kFALSE;
    }
  TString type = track->IsA()->GetName();
  const Double_t kBeampiperadius=3.;
  Double_t dcaD[2]={-999.,-999.},
      covD[3]={-999.,-999.,-999.};

  Double_t pos[3] = {0.,0.,0.};
  Bool_t isRecalcVertex(kFALSE);

  //case of ESD tracks
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(eev);
  if(!esdevent) {
      AliDebug(1, "No esd event available\n");
      return kFALSE;
    }
  const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
  if(!vtxESDSkip) return kFALSE;
  //case ESD track: take copy constructor
  const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
  if(tmptrack){

      if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex

          AliVertexerTracks vertexer(eev->GetMagneticField());
          vertexer.SetITSMode();
          vertexer.SetMinClusters(4);
          Int_t skipped[2];
          skipped[0] = track->GetID();
          vertexer.SetSkipTracks(1,skipped);
          //diamond constraint
          vertexer.SetConstraintOn();
          Float_t diamondcovxy[3];
          esdevent->GetDiamondCovXY(diamondcovxy);
          Double_t pos[3]={eev->GetDiamondX(),eev->GetDiamondY(),0.};
          Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
          AliESDVertex diamond(pos,cov,1.,1);
          vertexer.SetVtxStart(&diamond);
          vtxESDSkip = vertexer.FindPrimaryVertex(eev);
          isRecalcVertex = kTRUE;
        }

      if(vtxESDSkip){
          AliESDtrack esdtrack(*tmptrack);
          if(esdtrack.PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, impar, cov)){
              if(esdtrack.GetSigmaY2()<0. || esdtrack.GetSigmaZ2()<0.) {
                  if(isRecalcVertex) delete vtxESDSkip;
                  return kFALSE;
                  // this is insipired by the AliITStrackV2::Invariant() checks
                }
              //*****//
              //Calculate Sign
              Double_t posdcatrack[3]= {0.,0.,0.};
              esdtrack.GetXYZ(posdcatrack);
              vtxESDSkip->GetXYZ(pos);
              Double_t ipvector3[3] = { posdcatrack[0] - pos[0], posdcatrack[1] - pos[1], posdcatrack[2] - pos[2] };
              sign =TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz() );
              // Calculate decay legnth and track jet DCA against new vertex
              Double_t bpos[3] = { 0,0,0 };
              vtxESDSkip->GetXYZ(bpos);
              Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };
              Double_t bcv[21] = { 0 };
              AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
              Double_t xa = 0., xb = 0.;
              bjetparam.GetDCA(&esdtrack, eev->GetMagneticField(), xa, xb);


              Double_t xyz[3] = { 0., 0., 0. };
              Double_t xyzb[3] = { 0., 0., 0. };
              bjetparam.GetXYZAt(xa, eev->GetMagneticField(), xyz);
              esdtrack.GetXYZAt(xb, eev->GetMagneticField(), xyzb);
              Double_t  bdecaylength =
                  TMath::Sqrt(
                    (bpos[0] - xyz[0]) * (bpos[0] - xyz[0]) +
                  (bpos[1] - xyz[1]) * (bpos[1] - xyz[1]) +
                  (bpos[2] - xyz[2]) * (bpos[2] - xyz[2]));
              dcajetrack =
                  TMath::Sqrt(
                    (xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
                  (xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
                  (xyzb[2] - xyz[2]) * (xyzb[2]- xyz[2]));
              if(bdecaylength>0) lineardecaylength=bdecaylength;


              // cov[0] = esdtrack.GetCovariance()[0];
              if(isRecalcVertex) delete vtxESDSkip;
              return kTRUE;
            }
          else{
              delete vtxESDSkip;
            }
          if(isRecalcVertex) delete vtxESDSkip;
          return kFALSE;
        }
    }

  return kFALSE;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliVTrack *track, AliEmcalJet *jet, Double_t *impar, Double_t *cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength){
  if (fESD){
      return CalculateJetSignedTrackImpactParameter(((AliESDtrack*)track),jet ,impar,  cov, sign, dcajetrack,lineardecaylength);
    }
  return CalculateJetSignedTrackImpactParameter(((AliAODTrack*)track),jet ,impar,  cov, sign, dcajetrack,lineardecaylength);
}
// ######################################################################################## Post-process ImpPar
Double_t AliAnalysisTaskHFJetIPQA::GetValImpactParameter(TTypeImpPar type,Double_t *impar, Double_t * cov)
{
  Double_t result =-999999;
  Double_t dFdx = 0;
  Double_t dFdy = 0;
  //This will be pT dependent


  switch(type){
    case kXY:
      result = impar[0];
      break;
    case kXYSig:
      result = impar[0]/TMath::Sqrt(cov[0]);
      break;
    case kXYZ:
      result = TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
      break;
    case kXYZSig:
      result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
      dFdx = impar[0]/result;
      dFdy = impar[1]/result;
      result /=TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
      break;
    case kZSig:

      result =  impar[1];
      result /=TMath::Sqrt(cov[2]);
      break;
    case kXYZSigmaOnly:
      result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
      dFdx = impar[0]/result;
      dFdy = impar[1]/result;
      result =TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
      break;

    default:
      break;
    }
  return result;
}// ########################################################################################Track Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(AliVTrack* track,Int_t n){
  //Negative n -> SPD any , positve n ->SPD both
  if(track->GetX() >3 )return kFALSE; //valid propagation to dca only in beam pipe
  if(fESD){
      fESDTrackCut->SetPtRange(1.);
      fESDTrackCut->SetRequireITSRefit();
      fESDTrackCut->SetRequireTPCRefit();
      // fESDTrackCut->SetCutOutDistortedRegionsTPC(kTRUE);
      fESDTrackCut->SetEtaRange(-0.9,0.9);
      fESDTrackCut->SetAcceptKinkDaughters(kFALSE);
      fESDTrackCut->SetDCAToVertex2D(kTRUE);
      fESDTrackCut->SetMaxChi2PerClusterTPC(4);
      fESDTrackCut->SetMinNClustersTPC(120);
      fESDTrackCut->SetMinNClustersITS(n);
      fESDTrackCut->SetMaxChi2PerClusterITS(4);
      fESDTrackCut->SetMaxRel1PtUncertainty(0.1);
      fESDTrackCut->SetMaxFractionSharedTPCClusters(0.4);
      if(n>0)fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
      else fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      if(!(fESDTrackCut->AcceptTrack((AliESDtrack*)track))) return kFALSE;
      if(track->GetNcls(0)!=abs(n)) return kFALSE;
      return kTRUE;
    }
  else {
      if(!(((AliAODTrack*)track)->TestFilterBit(1 << 4))) return kFALSE;
      if(track->Pt() < 1.)return kFALSE;
      if(fabs(track->Eta()) > 0.9)return kFALSE;
      ULong_t status = track->GetStatus();
      if(!(status & AliAODTrack::kTPCrefit))return kFALSE;
      if(!(status & AliAODTrack::kITSrefit))return kFALSE;
      Float_t cRatioTPC = track->GetTPCNclsF() > 0. ?
            static_cast<Float_t>(track->GetTPCNcls()) / static_cast<Float_t>(track->GetTPCNclsF()) :
            1.;
      if(cRatioTPC < 0.6) return kFALSE;
      if(track->GetNcls(0)!=n) return kFALSE;
      if(track->GetNcls(1)<120) return kFALSE;
      return kTRUE;
    }
  return kFALSE;
}
// ######################################################################################## Jet matching 1/4
Bool_t AliAnalysisTaskHFJetIPQA::MatchJetsGeometricDefault()
{
  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
  Double_t matchingpar1 =0.25;
  Double_t matchingpar2 =0.25;
  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;
  DoJetLoop();
  AliEmcalJet* jet1 = 0;
  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {
      AliEmcalJet *jet2 = jet1->ClosestJet();
      if (!jet2) continue;
      if (jet2->ClosestJet() != jet1) continue;
      if (jet1->ClosestJetDistance() > matchingpar1 || jet2->ClosestJetDistance() > matchingpar2) continue;
      // Matched jet found
      jet1->SetMatchedToClosest(1);
      jet2->SetMatchedToClosest(1);
    }
  return kTRUE;
}
// ######################################################################################## Jet matching 2/4
void AliAnalysisTaskHFJetIPQA::DoJetLoop()
{
  // Do the jet loop.
  Double_t minjetpt =1.;
  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;
  AliEmcalJet* jet1 = 0;
  AliEmcalJet* jet2 = 0;
  jets2->ResetCurrentID();
  while ((jet2 = jets2->GetNextJet())) jet2->ResetMatching();
  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {
      jet1->ResetMatching();
      if (jet1->MCPt() < minjetpt) continue;
      jets2->ResetCurrentID();
      while ((jet2 = jets2->GetNextJet())) {
          SetMatchingLevel(jet1, jet2, 1);
        } // jet2 loop
    } // jet1 loop
}
// ######################################################################################## Jet matching 3/4
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor( AliAODMCParticle * mcpart,Int_t &pCorr_indx){
  if(!mcpart) {
      return 1;
    }
  Int_t pPdgcode = abs(mcpart->GetPdgCode());
  Int_t pMotherLabel = 	mcpart->GetMother();
  Bool_t found = kFALSE;
  Double_t pTWeight =0;
  Int_t foundPdg =-1;
  Int_t correctionidx =-1;
  Int_t maPdgcode = 0;
  Int_t maPdgcodeOld = 0;
  AliAODMCParticle * motherOld =0x0;
  Bool_t isprim = kFALSE;
  AliAODMCParticle * mcpartMother = 0x0;
  if(pMotherLabel<2) isprim =kTRUE;
  else mcpartMother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(pMotherLabel));
  if (!mcpartMother )isprim =kTRUE;
  if(!isprim){
      maPdgcode = abs(mcpartMother->GetPdgCode());
      if(maPdgcode>0 && maPdgcode<6) isprim =kTRUE;
      if(mcpartMother->GetStatus()>21) isprim =kTRUE;
      if(maPdgcode==21) isprim =kTRUE;
      if(IsPromptBMeson(mcpartMother)) {
          if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
              found = kTRUE;
              foundPdg =maPdgcode;
            }
        }
    }
  if(1==1){
      if(pPdgcode == bProton && 	isprim ){
          //Is Primary proton
          found = kTRUE;
          foundPdg =bProton;
          correctionidx = bIdxProton;
          pTWeight =mcpart->Pt();
        }
      else if(pPdgcode == bPi&& 	isprim){
          //Is Primary charged pion
          found = kTRUE;
          foundPdg =bPi;
          correctionidx = bIdxPi;
          pTWeight =mcpart->Pt();
        }
      else if(pPdgcode == bKaon&& 	isprim ){
          //Is Primary charged kaon
          found = kTRUE;
          pTWeight =mcpart->Pt();
          foundPdg =bKaon;
          correctionidx = bIdxKaon;
        }
      else if(!isprim && ParticleIsPossibleSource(maPdgcode)){
          if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
              found = kTRUE;
              foundPdg =maPdgcode;
            }
        }
      else {
          while (mcpartMother->GetMother() >0){
              maPdgcodeOld = maPdgcode;
              motherOld =mcpartMother;
              isprim=kFALSE;
              pMotherLabel = 	mcpartMother->GetMother();
              mcpartMother = 0x0;
              mcpartMother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(pMotherLabel)));
              if (!mcpartMother )isprim =kTRUE;
              maPdgcode = abs(mcpartMother->GetPdgCode());

              if(IsPromptBMeson(mcpartMother)) {
                  if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
                      found = kTRUE;
                      foundPdg =maPdgcode;
                      break;
                    }}
              if(!isprim){
                  maPdgcode = abs(mcpartMother->GetPdgCode());
                  if(maPdgcode>0 && maPdgcode<6) isprim =kTRUE;
                  if(mcpartMother->GetStatus()>21) isprim =kTRUE;
                  if(maPdgcode==21) isprim =kTRUE;
                }
              if(!isprim && ParticleIsPossibleSource(maPdgcode)){
                  if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
                      found = kTRUE;
                      foundPdg =maPdgcode;
                      break;
                    }
                }
              else if (!isprim) continue;
              else break;
            }
        }

    }
  if (!found) {

      return 1;
    }
  Double_t factor = 1;
  //calculate pt of array entry
  // 0.1 -25.GeV 0.05 per bin 498 bins
  Double_t wpos = ((pTWeight - 0.15)/ 0.05);
  Double_t  fractpart, intpart;
  fractpart = modf (wpos , &intpart);
  if (fractpart > 0) intpart = intpart + 1;
  Int_t  bin = floor(intpart);
  if (bin > 497) bin = 497;			// above weight definition
  if (pTWeight < 0.1+ 1E-5) bin = 0; //below weight definition
  factor = fBackgroundFactorLinus[correctionidx][bin];
  //if(correctionidx==19) Printf ("%f",factor);
  if (factor <= 0) return 1;
  else
    return factor;
}
// ######################################################################################## Jet matching 3/4
Bool_t AliAnalysisTaskHFJetIPQA::IsTruePrimary(AliMCParticle * mcpart){
  //Still need for material production
  if(!mcpart) return kFALSE;
  AliMCParticle * mcmother = (AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother());
  if(!mcmother) return kTRUE;
  Int_t istatus = mcmother->Particle()->GetStatusCode();
  if(istatus >11)return kTRUE;
  Int_t ipdg =abs(mcmother->PdgCode())	;
  return kFALSE;
}
// ######################################################################################## Jet matching 3/4
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor( AliMCParticle * mcpart,Int_t &pCorr_indx){
  if(!mcpart)  return 1; // No corresponding Montecarlo particle was found
  Bool_t _particlesourcefound(kFALSE);
  Int_t  _particlesourcepdg(abs(mcpart->PdgCode()));
  Int_t  _particlesourceidx(-1);
  Double_t _particlesourcept(0);

  AliMCParticle * mcpartclone = mcpart;
  while(mcpart){//omega and xi test
      if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
      _particlesourcept = mcpart->Pt();
      _particlesourcepdg = abs(mcpart->PdgCode());
      if (IsSelectionParticleOmegaXiSigmaP(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
          _particlesourcefound = kTRUE;
          break;
        }
      mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
    }
  if (!_particlesourcefound) { //heavy mesons to improve templates
      mcpart = mcpartclone;

      while(mcpart){
          if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;

          _particlesourcept = mcpart->Pt();
          _particlesourcepdg = abs(mcpart->PdgCode());
          if (IsSelectionParticleStrange(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
              _particlesourcefound = kTRUE;
              break;
            }
          mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
        }
    }
  if (!_particlesourcefound) { //heavy mesons to improve templates
      mcpart = mcpartclone;
      while(mcpart){
          if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;

          _particlesourcept = mcpart->Pt();
          _particlesourcepdg = abs(mcpart->PdgCode());
          if (IsSelectionParticleMeson(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
              _particlesourcefound = kTRUE;
              break;
            }
          mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
        }
    }
  if (!_particlesourcefound) { //charged hadrons
      mcpart = mcpartclone;
      while(mcpart){
          if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;

          _particlesourcept = mcpart->Pt();
          _particlesourcepdg = abs(mcpart->PdgCode());
          if (IsSelectionParticleALICE(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
              _particlesourcefound = kTRUE;
              break;
            }
          mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
        }
    }
  if (!_particlesourcefound) {
      mcpart = mcpartclone;
      while(mcpart){
          if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
          _particlesourcept = mcpart->Pt();
          if (IsSelectionParticle(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
              _particlesourcefound = kTRUE;
              break;
            }
          mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
        }
    }
  if (!_particlesourcefound) return 1.;
  // Do the weighting
  Double_t factor = 1;
  if (_particlesourceidx <0) return 1;
  if (_particlesourceidx >19) return 1;
  //calculate pt of array entry
  // 0.1 -25.GeV 0.05 per bin 498 bins
  Double_t wpos = ((_particlesourcept - 0.15)/ 0.05);
  Double_t  fractpart, intpart;
  fractpart = modf (wpos , &intpart);
  if (fractpart > 0) intpart = intpart + 1;
  Int_t  bin = floor(intpart);
  if (bin > 497) bin = 497;			// above weight definition
  if (_particlesourcept < 0.1+ 1E-5) bin = 0; //below weight definition
  factor = fBackgroundFactorLinus[_particlesourceidx][bin];
  pCorr_indx = mcpart->GetLabel();
  Double_t flucafactor = 1;
  //Omega + xi

  switch(mcpart->PdgCode())
    {
    /*dirty hack to include xi and omega- K0(0,+) phi <https://arxiv.org/pdf/1208.5717v2.pdf*/
    case bPhi:
      factor=1;
      flucafactor =1./fPhi->Eval(_particlesourcept);
      break;
    case bK0S892:
      factor=1;
      flucafactor =1./fK0Star->Eval(_particlesourcept);
      break;
    case bK0S892plus:
      factor=1;

      flucafactor =1./fK0Star->Eval(_particlesourcept);
      //Printf("Factor %f",flucafactor);

      break;
    case bOmegaBaryon:
      flucafactor =fGraphOmega->Eval(_particlesourcept);
      break;
    case bXiBaryon:
      flucafactor =fGraphXi->Eval(_particlesourcept);
      break;
    case bLambda:
      flucafactor =fGeant3FlukaLambda->Eval(_particlesourcept);
      break;
    case -bLambda:
      flucafactor =fGeant3FlukaAntiLambda->Eval(_particlesourcept);
      break;
    case bProton:
      flucafactor =fGeant3FlukaProton->Eval(_particlesourcept);
      break;
    case -bProton:
      flucafactor =fGeant3FlukaAntiProton->Eval(_particlesourcept);
      break;
    case -bKaon:
      flucafactor =fGeant3FlukaKMinus->Eval(_particlesourcept);
      break;
    default:
      break;
    }

  factor*=flucafactor;
  if (factor <= 0)  factor = 1.;
  if (factor > 1E2) factor = 1;
  return factor;
}
// ########################################################################################
void AliAnalysisTaskHFJetIPQA::SmearTrackHybrid(AliVTrack *track)
{
  if (!track || !fIsPythia) return;
  if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))return;
  if (TESTBIT(track->GetITSClusterMap(),7)) return;

  AliExternalTrackParam et;et.CopyFromVTrack(track);
  Double_t *param=const_cast<Double_t*>(et.GetParameter());

  // Get MC info
  Int_t imc=track->GetLabel();
  if (imc>0) {

      AliAODMCParticle * mc = 0x0;
      fESD ? mc = new AliAODMCParticle((AliMCParticle*)MCEvent()->GetTrack(imc),0,0)
          :static_cast<AliAODMCParticle*>(fMCArray->At(imc));

      if (mc!=0x0) {
          Double_t mcx[3];
          Double_t mcp[3];
          Double_t mccv[36]={0.};
          Short_t  mcc;
          mc->XvYvZv(mcx);
          mc->PxPyPz(mcp);
          mcc=mc->Charge();
          AliExternalTrackParam mct(mcx,mcp,mccv,mcc);
          const Double_t *parammc=mct.GetParameter();
          AliVertex vtx(mcx,1.,1);


          et.PropagateToDCA(&vtx,track->GetBz(),3.);
          et.Rotate(mct.GetAlpha());
          Double_t ptmc=TMath::Abs(mc->Pt());
          Double_t spt1o =0.;
          // Apply the smearing
          Double_t d0zo  =param  [1];
          Double_t d0zmc =parammc[1];
          Double_t d0rpo =param  [0];
          Double_t d0rpmc=parammc[0];
          Double_t pt1o  =param  [4];
          Double_t pt1mc =parammc[4];
          Double_t xpt = track->Pt();
          //if(xpt <1.5)  xpt = 1.5;
          if(xpt >8)  xpt = 8;
          /*
          Double_t valmc   =   fCurrentSigmaFactorsMC[0]-(fCurrentSigmaFactorsMC[1]*TMath::Exp(-(1*(fCurrentSigmaFactorsMC[2]*(xpt+fCurrentSigmaFactorsMC[3])))));
          Double_t valdata =  fCurrentSigmaFactorsData[0]-(fCurrentSigmaFactorsData[1]*TMath::Exp(-(1*(fCurrentSigmaFactorsData[2]*(xpt+fCurrentSigmaFactorsData[3])))));
         */
          Double_t valmc   = fGraphSigmaMC->Eval(xpt);
          Double_t valdata =   fGraphSigmaData->Eval(xpt);

          valmc *=1e-4;
          valdata *=1e-4;

          // Printf("valmc %f valdata %f",valmc,valdata);
          /* Double_t x_gr_z[7] = {0,1.500000,2.500000,3.500000,4.500000,7.500000,1000};
             Double_t y_gr_z[7] = {1.080992,1.080992,1.086903,1.102976,1.143597,1.233094,1.233094};

             TGraph factorZ(7,x_gr_z,y_gr_z);
   */
          //Thid is a test
          //    Double_t dd0zo =d0zo-d0zmc;
          //        Double_t dd0zn =dd0zo *factorZ.Eval(ptmc);
          //      Double_t d0zn  =d0zmc+dd0zn;
          Double_t dd0rpo=d0rpo-d0rpmc;
          Double_t dd0rpn=dd0rpo *  valdata/valmc;
          Double_t  d0rpn=d0rpmc+dd0rpn;
          Double_t dpt1o =pt1o-pt1mc;
          Double_t dpt1n =dpt1o *(spt1o >0. ? 1: 1.);
          Double_t pt1n  =pt1mc+dpt1n;



          param[0]=d0rpn;
          //param[1]=d0zn ;
          param[4]=pt1n ;
          AliAODTrack tmp;
          Double_t x[3];
          Double_t p[3];

          et.GetXYZ(x);
          et.GetPxPyPz(p);

          if(fESD) ((AliESDtrack*)track)->SetParamOnly(et.GetX(),et.GetAlpha(),param);
          else {

              ((AliAODTrack*)track)->SetPosition(x,kFALSE);
              ((AliAODTrack*)track)->SetP(p,kTRUE);
            }


          UChar_t itsClusterMap = ((AliESDtrack*)track)->GetITSClusterMap();
          SETBIT(itsClusterMap,7);
          ((AliESDtrack*)track)->SetITSClusterMap(itsClusterMap);
          if (mc && fESD) {
              delete mc;
              mc=0x0;
            }
        }
      else Printf("MC Track not found");
      if (mc && fESD) delete mc;
      mc=0x0;
    }
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliAODMCParticle * particle ) {
  // If a particle is not a physical primary, check if it comes from weak decay
  Int_t mfl = 0;
  Int_t indexMoth = particle->GetMother();
  if(indexMoth < 0) return kFALSE; // if index mother < 0 and not a physical primary, is a non-stable product or one of the beams
  AliAODMCParticle* moth = dynamic_cast<AliAODMCParticle *>(fMCArray->At(indexMoth));
  Int_t pcodemoth = TMath::Abs(particle->PdgCode());
  Int_t codemoth = TMath::Abs(moth->PdgCode());
  // mass of the flavour
  mfl = Int_t (codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
  // if(mfl == 4|| mfl ==5) return kTRUE;
  if(pcodemoth==bPi0){
      if(codemoth == bK0s) return kTRUE;
      if(codemoth == bK0l) return kTRUE;
      if(TMath::Abs(codemoth) == bKaon) return kTRUE;
      if(TMath::Abs(codemoth) == bLambda) return kTRUE;
      if(codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113 || codemoth == bRhoPlus) return kTRUE;
    }
  else if (pcodemoth==bPhi){
      if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
    }
  else if (pcodemoth==bOmega){
      if(codemoth == 111 || codemoth == 221 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
    }
  else if (pcodemoth==bEtaPrime){
      if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
    }
  else if (pcodemoth==bRho){
      if(codemoth== 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331|| codemoth == bRhoPlus) return kTRUE;
    }
  else if (pcodemoth==bEta){
      if(codemoth == 111 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliMCParticle * particle ) {
  // If a particle is not a physical primary, check if it comes from weak decay
  Int_t mfl = 0;
  Int_t indexMoth = particle->GetMother();
  if(indexMoth < 0) return kFALSE; // if index mother < 0 and not a physical primary, is a non-stable product or one of the beams
  AliMCParticle* moth =(AliMCParticle*) MCEvent()->GetTrack(indexMoth);
  if(!moth) return kFALSE;
  Int_t codemoth = TMath::Abs(moth->PdgCode());
  Int_t checklist [9] ={310,130,3122,3212,3222,3112,3322,3312,3334};
  for (Int_t i = 0; i<9; ++i) {
      if(codemoth == checklist[i])return kTRUE;

    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliAODMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  )
{
  pT = mcpart->Pt();
  switch(pdg){
    case bBPlus:
      idx = bIdxBPlus;
      return kTRUE;
    case bB0:
      idx = bIdxB0;
      return kTRUE;
    case bLambdaB:
      idx = bIdxLambdaB;
      return kTRUE;
      break;
    case bBStarPlus:
      idx = bIdxBStarPlus;
      return kTRUE;
      break;
    default:
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  )
{
  pT = mcpart->Pt();
  switch(pdg){
    case bBPlus:
      idx = bIdxBPlus;
      return kTRUE;
    case bB0:
      idx = bIdxB0;
      return kTRUE;
    case bLambdaB:
      idx = bIdxLambdaB;
      return kTRUE;
      break;
    case bBStarPlus:
      idx = bIdxBStarPlus;
      return kTRUE;
      break;
    default:
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliAODMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT = mcpart->Pt();
  Bool_t isprim =kFALSE;
  Int_t mpdg = 0;
  if(mcpart->GetMother()>-1)mpdg = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(mcpart->GetMother())))->GetPdgCode();
  if(mpdg>0 && mpdg<6) isprim =kTRUE;
  if(mpdg==21) isprim =kTRUE;
  if (MCEvent()->Stack()->IsSecondaryFromMaterial(abs(mcpart->GetLabel())))return kFALSE;
  if(dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(mcpart->GetMother())))->GetStatus()>11) isprim =kTRUE;
  switch(pdg){
    case bPi0:
      if(!IsSecondaryFromWeakDecay(mcpart) ){
          idx = bIdxPi0;
          return kTRUE;
        }
      break;
    case bEta:
      if(!IsSecondaryFromWeakDecay(mcpart) ){
          idx = bIdxEta;
          return kTRUE;
        }
      break;
    case bEtaPrime:
      if(!IsSecondaryFromWeakDecay(mcpart) ){
          idx = bIdxEtaPrime;
          return kTRUE;
        }
      break;
    case bOmega:
      idx = bIdxOmega;
      if(!IsSecondaryFromWeakDecay(mcpart)){
          return kTRUE;
        }
      break;
    case bPhi:
      idx = bIdxPhi;
      if(!IsSecondaryFromWeakDecay(mcpart)){
          return kTRUE;
        }
      break;
    case bRho:
      idx = bIdxRho;
      if(!IsSecondaryFromWeakDecay(mcpart)){
          return kTRUE;
        }
      break;
    case bD0:
      idx = bIdxD0;
      if(IsPromptDMeson(mcpart)){
          return kTRUE;
        }
      break;
    case bDPlus:
      idx = bIdxDPlus;
      if(IsPromptDMeson(mcpart)){
          return kTRUE;
        }
      break;
    case bIdxDSPlus:
      idx = bDSPlus;
      if(IsPromptDMeson(mcpart)){
          return kTRUE;
        }
      break;
    case bDStarPlus:
      idx = bIdxDStarPlus;
      if(IsPromptDMeson(mcpart)){
          return kTRUE;
        }
      break;
    case bLambdaC:
      idx = bIdxLambdaC;
      if(IsPromptDMeson(mcpart)){
          return kTRUE;
        }
      break;
    case bLambda:
      idx = bIdxLambda;
      if(mcpart->IsPhysicalPrimary())
        {
          return kTRUE;
        }
      break;
    case bK0s:
      idx = bIdxK0s;

      if(mcpart->IsPhysicalPrimary())
        {
          return kTRUE;
        }
      break;
    case bPi:
      if(isprim){
          idx =bIdxPi;
          return kTRUE;
        }
      break;
    case bKaon:
      if(isprim){
          idx =bIdxKaon;
          return kTRUE;
        }
      break;
    default:
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT 	= mcpart->Pt();
  Int_t pdg2 = abs(mcpart->PdgCode());
  idx = -1;

  switch(pdg2){
    case bPi0:
      idx = bIdxPi0;
      if(IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bEta:
      idx = bIdxEta;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bEtaPrime:
      idx = bIdxEtaPrime;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bOmega:
      idx = bIdxOmega;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bPhi:
      idx = bIdxPhi;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bRho:
      idx = bIdxRho;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    case bRhoPlus: //Experimental assume same shape correction for neutral and charged rho
      idx = bIdxRho;
      if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
      break;
    default:
      break;
    }
  return kFALSE;
}// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleALICE( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT 	= mcpart->Pt();
  AliMCParticle * mother = 0x0;
  idx = -1;
  pdg = abs(mcpart->PdgCode());
  Bool_t pIsSecStrANGE= MCEvent()->Stack()->IsSecondaryFromWeakDecay(mcpart->GetLabel());

  if(IsSecondaryFromWeakDecay(mcpart) && pIsSecStrANGE ) return kFALSE;

  switch(pdg){
    case bProton:
      mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
      if(mother){
          if((abs(mother->PdgCode()) ==  3222) )
            return kFALSE;
        }
      idx=bIdxProton;
      return kTRUE;
      break;
    case bPi:
      idx=bIdxPi;
      return kTRUE;
      break;
    case bKaon:
      idx=bIdxKaon;
      return kTRUE;
      break;

    default:
      return kFALSE;
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleMeson( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT 	= mcpart->Pt();
  AliMCParticle * mother = 0x0;
  idx = -1;
  pdg = abs(mcpart->PdgCode());
  switch(pdg){
    case bD0:
      idx = bIdxD0;
      if(IsPromptDMeson(mcpart))return kTRUE;
      break;
    case bDPlus:
      idx = bIdxDPlus;
      if(IsPromptDMeson(mcpart))return kTRUE;
      break;
    case bDSPlus:
      idx = bIdxDSPlus;
      if(IsPromptDMeson(mcpart))return kTRUE;
      break;
    case bDStarPlus:
      idx = bIdxDStarPlus;
      if(IsPromptDMeson(mcpart))return kTRUE;
      break;
    case bLambdaC:
      idx = bIdxLambdaC;
      if(IsPromptDMeson(mcpart))return kTRUE;
      break;
    case bBPlus:
      idx = bIdxBPlus;
      if(IsPromptBMeson(mcpart))return kTRUE;
      break;
    case bB0:
      idx = bIdxB0;
      if(IsPromptBMeson(mcpart))return kTRUE;
      break;
    case bLambdaB:
      idx = bIdxLambdaB;
      if(IsPromptBMeson(mcpart))return kTRUE;
      break;
    case bBStarPlus:
      idx = bIdxBStarPlus;
      if(IsPromptBMeson(mcpart))return kTRUE;
      break;
    default:
      return kFALSE;
      break;
    }
  return kTRUE;
}
//##############################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleOmegaXiSigmaP( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT 	= mcpart->Pt();
  AliMCParticle * mother = 0x0;
  idx = -1;
  pdg = abs(mcpart->PdgCode());

  Bool_t pIsPhysicalPrimary= MCEvent()->IsPhysicalPrimary(mcpart->GetLabel());
  if (!pIsPhysicalPrimary) return kFALSE;
  switch(pdg){
    case bXiBaryon:
      idx = bIdxBStarPlus;//dummy!
      return kTRUE;
      break;
    case bOmegaBaryon:
      idx = bIdxBStarPlus;//dummy!
      return kTRUE;
      break;
    default:
      return kFALSE;
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleStrange( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
  pT 	= mcpart->Pt();
  AliMCParticle * mother = 0x0;
  idx = -1;
  pdg = abs(mcpart->PdgCode());

  Bool_t pIsPhysicalPrimary= MCEvent()->IsPhysicalPrimary(mcpart->GetLabel());


  switch(pdg){
    case bPhi:
      idx = bIdxPhi;//dummy will be overwritten in posrpos
      return kTRUE;
      break;
    case bK0S892:
      idx = bIdxK0s;//dummy will be overwritten in posrpos
      return kTRUE;
      break;
    case bK0S892plus:
      idx = bIdxK0s; //dummy will be overwritten in posrpos
      return kTRUE;
      break;
    case bK0s:
      idx = bIdxK0s;
      mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
      if(mother){
          if((abs(mother->PdgCode()) == bPhi))
            return kFALSE;
        }
      if(pIsPhysicalPrimary) return kTRUE;
      break;
    case bK0l:
      idx = bIdxK0s;
      mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
      if(mother){
          if((abs(mother->PdgCode()) == bPhi))
            return kFALSE;
        }
      if(pIsPhysicalPrimary)  return kTRUE;
      break;

    case bLambda:
      idx = bIdxLambda;

      mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
      if(mother){
          if((abs(mother->PdgCode()) ==  3312) || (abs(mother->PdgCode()) ==  3322) || (abs(mother->PdgCode()) ==  3334))
            return kFALSE;
        }
      return kTRUE;
      break;
    default:
      return kFALSE;
      break;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliAODMCParticle * part )
{
  if(!part) return kFALSE;
  Int_t pdg = TMath::Abs(part->GetPdgCode());
  if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
    {
      Int_t imo =  part->GetMother();
      AliAODMCParticle* pm = dynamic_cast<AliAODMCParticle *>(fMCArray->At(imo));
      Int_t mpdg = TMath::Abs(pm->GetPdgCode());
      if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
        return kTRUE;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliMCParticle * part )
{
  if(!part) return kFALSE;
  Int_t pdg = TMath::Abs(part->PdgCode());
  if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
    {
      Int_t imo =  part->GetMother();
      AliMCParticle* pm = dynamic_cast<AliMCParticle *>(MCEvent()->GetTrack(imo));
      Int_t mpdg = TMath::Abs(pm->PdgCode());
      if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
        return kTRUE;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliAODMCParticle * part )
{
  if(!part) return kFALSE;
  Int_t pdg = TMath::Abs(part->GetPdgCode());
  if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
    {
      Int_t imo =  part->GetMother();
      if(imo<0) return kTRUE;
      AliAODMCParticle* pm = dynamic_cast<AliAODMCParticle *>(fMCArray->At(imo));
      if(!pm) return kTRUE;
      Int_t mpdg = TMath::Abs(pm->GetPdgCode());
      if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
        return kTRUE;
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliMCParticle * part )
{
  if(!part) return kFALSE;
  Int_t pdg = TMath::Abs(part->PdgCode());
  if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
    {
      Int_t imo =  part->GetMother();
      if(imo<0) return kTRUE;
      AliMCParticle* pm = ((AliMCParticle*)fMCEvent->GetTrack(abs(imo)));
      if(!pm) return kTRUE;
      Int_t mpdg = TMath::Abs(pm->PdgCode());
      if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
        return kTRUE;
    }

  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::ParticleIsPossibleSource(Int_t pdg){
  Int_t pos[22] = {bPi0,bEta,bEtaPrime,bPhi,bRho,bOmega,bK0s,bLambda,bOmegaBaryon,bXiBaryon,bD0,bPi,bKaon,bProton,bDPlus,bDStarPlus,bDSPlus,bLambdaB,bLambdaC,bBPlus,bB0,bBStarPlus};
  for (Int_t i =0 ;i<22 ;++i){
      if (abs(pdg)==pos[i] ) return kTRUE;
    }
  return kFALSE;
}
// ######################################################################################## Jet matching 3/4
void AliAnalysisTaskHFJetIPQA::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Int_t matching)
{
  Double_t d1 = -1;
  Double_t d2 = -1;

  switch (matching) {
    case 1:
      GetGeometricalMatchingLevel(jet1,jet2,d1);
      d2 = d1;
      break;
    default:
      break;
    }
  if (d1 >= 0) {

      if (d1 < jet1->ClosestJetDistance()) {
          jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
          jet1->SetClosestJet(jet2, d1);
        }
      else if (d1 < jet1->SecondClosestJetDistance()) {
          jet1->SetSecondClosestJet(jet2, d1);
        }
    }
  if (d2 >= 0) {
      if (d2 < jet2->ClosestJetDistance()) {
          jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
          jet2->SetClosestJet(jet1, d2);
        }
      else if (d2 < jet2->SecondClosestJetDistance()) {
          jet2->SetSecondClosestJet(jet1, d2);
        }
    }
}

// ######################################################################################## Jet matching 4/4
void AliAnalysisTaskHFJetIPQA::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
  Double_t deta = jet2->Eta() - jet1->Eta();
  Double_t dphi = jet2->Phi() - jet1->Phi();
  dphi = TVector2::Phi_mpi_pi(dphi);
  d = sqrt(deta * deta + dphi * dphi);
}

// ######################################################################################## Monte Carlo correction factors
Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliVTrack* track,Int_t &pCorr_indx){
  AliAODMCParticle *pMCAOD = 0x0;
  AliMCParticle *pMCESD = 0x0;
  if(!fESD) pMCAOD = (AliAODMCParticle*)fMCArray->At(abs(track->GetLabel()));
  else pMCESD = ((AliMCParticle*)MCEvent()->GetTrack(abs(track->GetLabel())));
  if(!(pMCESD || pMCAOD)) return 1;
  Double_t val = 1;
  val=  fESD ? GetWeightFactor(pMCESD,pCorr_indx): GetWeightFactor(pMCAOD,pCorr_indx);
  if(val > 0 ){
      return val;
    }
  return 1.;
}
// ########################################################################################

Bool_t AliAnalysisTaskHFJetIPQA::mysort(const SJetIpPati& i, const SJetIpPati& j)
{
  if(i.first <= j.first)
    return kFALSE;
  else
    return kTRUE;
}
// ########################################################################################
AliAODMCParticle* AliAnalysisTaskHFJetIPQA::GetMCTrack( const AliAODTrack* _track)
{
  //
  // return MC track
  //
  if(!fMCArray) { AliError("No fMCArray"); return NULL;}
  Int_t nStack = fMCArray->GetEntriesFast();
  Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
  if(label > nStack) return NULL;
  AliAODMCParticle *mctrack = (AliAODMCParticle*)fMCArray->At(label);
  return mctrack;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsV0PhotonFromBeamPipeDaughter(const AliAODTrack* track)
{
  if(!track)return kFALSE;
  AliAODv0* v0aod = 0x0;
  Int_t posid = -1;
  Int_t negid = -1;
  Int_t trackid = -1;
  Double_t P[3];
  for(Int_t i = 0; i < InputEvent()->GetNumberOfV0s(); ++i) {
      P[0]=0.;
      P[1]=0.;
      P[2]=0.;
      v0aod = ((AliAODEvent*)InputEvent())->GetV0(i);
      if (!v0aod->GetOnFlyStatus()) continue;
      posid = v0aod->GetPosID();
      negid = v0aod->GetNegID();
      trackid = track->GetID();
      if(posid == trackid || negid == trackid) {
          return kTRUE;
          P[0] = v0aod->DecayVertexV0X();
          P[1] = v0aod->DecayVertexV0Y();
          P[2] = v0aod->DecayVertexV0Z();
          Double_t Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
          if(Radius < 800.) {
              //Try To construct gamma from daughters
              AliVTrack* posTrack = ((AliVTrack*)(InputEvent()->GetTrack(posid))) ;
              AliVTrack* negTrack = ((AliVTrack*)(InputEvent()->GetTrack(negid))) ;
              if(!posTrack) continue;
              if(!negTrack) continue;
              AliKFParticle::SetField((((AliAODEvent*)(InputEvent()))->GetMagneticField()));
              AliKFParticle pos(*posTrack,-11);
              AliKFParticle neg(*negTrack,11);
              AliKFParticle partGamma(pos,neg,kTRUE);
              // Caluclate Armenteros Podolanski cut
              Double_t pPlus[3]  ={	pos.GetPx() ,	pos.GetPy() ,	pos.GetPz() };
              Double_t pMinus[3] ={	neg.GetPx() ,	neg.GetPy() ,	neg.GetPz() };
              Double_t PV0[3] ={	partGamma.GetPx() ,	partGamma.GetPy() ,	partGamma.GetPz() };
              TVector3 pPlusV3(pPlus[0],pPlus[1],pPlus[2]);
              TVector3 pMinusV3(pMinus[0],pMinus[1],pMinus[2]);
              TVector3 pV0V3(PV0[0],PV0[1],PV0[2]);
              TVector3 v = pPlusV3.Cross(pV0V3);
              Double_t pPerp = v.Mag();
              Double_t qT =pPerp /pV0V3.Mag();
              // Calculate psi Pair value
              Double_t xipair = acos((pPlusV3 * pMinusV3)/(pPlusV3.Mag() * pMinusV3.Mag()));
              Double_t theataPlus = ((AliVTrack*)(InputEvent()->GetTrack(posid)))->Theta();
              Double_t theataMinus = ((AliVTrack*)(InputEvent()->GetTrack(negid)))->Theta();
              Double_t delta_theta =theataMinus-theataPlus;
              Double_t psi_pair = asin(	delta_theta/xipair);
              if((partGamma.Chi2()/partGamma.NDF() <30.) && (qT < 0.05) && (fabs(psi_pair) < 0.05) ) {
                  Printf("Rejected as Gamma");
                  return kTRUE;
                }
            }
          return kFALSE;
        }
    }
  return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsV0PhotonFromBeamPipeDaughter(const AliESDtrack* track)
{
  //Still need inclusion of Track + PID cuts
  if(!track)return kFALSE;
  AliESDv0* v0esd = 0x0;
  Int_t posid = -1;
  Int_t negid = -1;
  Int_t trackid = -1;
  Double_t P[3];
  for(Int_t i = 0; i < InputEvent()->GetNumberOfV0s(); ++i) {
      P[0]=0.;
      P[1]=0.;
      P[2]=0.;
      v0esd = ((AliESDEvent*)InputEvent())->GetV0(i);
      if (!v0esd->GetOnFlyStatus())  continue;
      v0esd->XvYvZv(P);
      Double_t Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
      if(Radius < 800.) {
          AliKFParticle::SetField((((AliESDEvent*)(InputEvent()))->GetMagneticField()));
          AliKFParticle pos(*(v0esd->GetParamP()),-11);
          AliKFParticle neg(*(v0esd->GetParamN()),11);
          AliKFParticle partGamma(pos,neg,kTRUE);
          // Caluclate Armenteros Podolanski cut
          Int_t NIndex =v0esd->GetNindex();
          Int_t PIndex =v0esd->GetPindex();
          Int_t PLabel = ((AliESDEvent*)InputEvent())->GetTrack(PIndex)->GetLabel();
          Int_t NLabel = ((AliESDEvent*)InputEvent())->GetTrack(NIndex)->GetLabel();
          if((PLabel==abs(track->GetLabel()))&&(!(NLabel==abs(track->GetLabel())))) return kTRUE;
          else return kFALSE;
          Double_t alpha =0;
          Double_t qT = GetArmenteros(v0esd,11,-11,alpha);
          Double_t psi_pair = 	GetPsiPair(v0esd);
          if((partGamma.GetChi2()/partGamma.GetNDF() <30.) && (qT < 0.05) && (fabs(psi_pair) < 0.05) ) {
              Printf("Rejected as Gamma");
              return kTRUE;
            }
          return kFALSE;
        }
    }
  return kFALSE;
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetArmenteros(AliESDv0 * v0, Int_t pidneg,Int_t pidpos ,Double_t &alpha){
  AliKFParticle pos(*(v0->GetParamP()),pidpos);
  AliKFParticle neg(*(v0->GetParamN()),pidneg);
  AliKFParticle partGamma(pos,neg,kTRUE);
  Double_t armenteros[2]={0,0};
  partGamma.GetArmenterosPodolanski(pos,neg, armenteros);
  alpha= armenteros[1];
  return armenteros[0];
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetPsiPair(AliESDv0 * v0){
  AliExternalTrackParam nt(*v0->GetParamN());
  AliExternalTrackParam pt(*v0->GetParamP());
  Float_t magField = fInputEvent->GetMagneticField();
  Double_t xyz[3] = {0.,0.,0.};
  v0->GetXYZ(xyz[0],xyz[1],xyz[2]);
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);
  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
  Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated
  Double_t momPosProp[3] = {0,0,0};
  Double_t momNegProp[3] = {0,0,0};
  Double_t psiPair = 4.;
  if(nt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
  if(pt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);
  Double_t pEle =
      TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos =
      TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
  Double_t scalarproduct =
      momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks
  psiPair = TMath::ASin(deltat/chipair);
  return psiPair;
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrected(const AliEmcalJet *jet)
{
  AliJetContainer * jetconrec = 0x0;
  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  if(jet && jetconrec)
    return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
  return -1.;
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrectedMC(const AliEmcalJet *jet)
{
  AliJetContainer * jetconrec = 0x0;
  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(1));
  if(jet && jetconrec)
    return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
  return -1.;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::GetESDITSMODULEINFO(AliESDtrack * track){
  Int_t nclsSDDSSD = 0;
  Int_t idet,status; Float_t xloc,zloc;

  const char * transl_status[7] ={
    "found (cluster is associated)",
    "dead (module is dead from OCDB)",
    "skipped (module or layer forced to be skipped)",
    "outinz (track out of z acceptance)",
    "nocls (no clusters in the road)",
    "norefit (cluster rejected during refit)",
    "deadzspd (holes in z in SPD)"};

  for(Int_t layer=0; layer<6; layer++) {
      if(layer>=2 && track->HasPointOnITSLayer(layer)) nclsSDDSSD++;
      status =  track->GetITSModuleIndexInfo(layer,idet,status,xloc,zloc);
      if(status<0) continue;
      if(layer>=2) idet+=240; // add n SPD modules
      if(layer>=4) idet+=260; // add n SDD modules
      if (status>0)Printf("Status: %s  Module %i",transl_status[status-1],idet);
    }
  return kTRUE;
}
//###############################################################################################################
Int_t  AliAnalysisTaskHFJetIPQA::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
{
  if(!jet) return 0;
  if(! fMcEvtSampled){
      //Sample MC stack upward once per event and oly if there are jets
      AliStack * Mstack = MCEvent()->Stack();
      for(Int_t iPrim = 0 ; iPrim<Mstack->GetNprimary();iPrim++){
          TParticle * part = (TParticle*)Mstack->Particle(iPrim);
          if(!part) return 0;
          if(!((part->GetStatusCode()==11) ||(part->GetStatusCode()==12))) continue;
          Double_t etap = part->Eta();
          Double_t phip = part->Phi();
          int pdg = (abs(part->GetPdgCode()));
          if(pdg == 5) {
              fEtaBEvt.push_back(etap);
              fPhiBEvt.push_back(phip);
            }
          else if(pdg== 4) {
              fEtaCEvt.push_back(etap);
              fPhiCEvt.push_back(phip);
            }
          else if(pdg == 3 ) {
              fEtaSEvt.push_back(etap);
              fPhiSEvt.push_back(phip);
            }
          else if(pdg== 1 ||pdg== 2 || pdg== 3 || pdg == 21) {
              fEtaUdsgEvt.push_back(etap);
              fPhiUdsgEvt.push_back(phip);
            }
        }
      fMcEvtSampled= kTRUE;
    }
  if(fEtaBEvt.size() ==0&& fEtaCEvt.size()==0&& fEtaSEvt.size()==0&&fEtaUdsgEvt.size()==0) return 0; //udsg
  Double_t etajet = jet->Eta();
  Double_t phijet = jet->Phi();
  //check for c jet
  for (Int_t icj = 0 ; icj <(Int_t)fEtaCEvt.size();++icj ){
      Double_t eta =fEtaCEvt.at(icj);
      Double_t phi =fPhiCEvt.at(icj);
      Double_t deta = etajet - eta;
      Double_t dphi = phijet-phi;
      dphi = TVector2::Phi_mpi_pi(dphi);
      Double_t  d = sqrt(deta * deta + dphi * dphi);
      if(d < radius) return 2;
    }
  //check for b jet
  for (Int_t icj = 0 ; icj <(Int_t)fEtaBEvt.size();++icj ){
      Double_t eta =fEtaBEvt.at(icj);
      Double_t phi =fPhiBEvt.at(icj);
      Double_t deta = etajet - eta;
      Double_t dphi = phijet - phi;
      dphi = TVector2::Phi_mpi_pi(dphi);
      Double_t  d = sqrt(deta * deta + dphi * dphi);
      if(d < radius) return 3;
    }
  //check for s jet
  for (Int_t icj = 0 ; icj <(Int_t)fEtaSEvt.size();++icj ){
      Double_t eta =fEtaSEvt.at(icj);
      Double_t phi =fPhiSEvt.at(icj);

      Double_t deta = etajet - eta;
      Double_t dphi = phijet - phi;
      dphi = TVector2::Phi_mpi_pi(dphi);
      Double_t  d = sqrt(deta * deta + dphi * dphi);
      if(d < radius) return 1;
    }
  for (Int_t icj = 0 ; icj <(Int_t)fEtaUdsgEvt.size();++icj ){
      Double_t eta =fEtaUdsgEvt.at(icj);
      Double_t phi =fPhiUdsgEvt.at(icj);

      Double_t deta = etajet - eta;
      Double_t dphi = phijet - phi;
      dphi = TVector2::Phi_mpi_pi(dphi);
      Double_t  d = sqrt(deta * deta + dphi * dphi);
      if(d < radius){
          is_udg =kTRUE;
          return 1;
        }
    }

  return 0;
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t w){
  TH1D * h1 =GetHist1D(name);
  h1->Fill(x,w);
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t y, Double_t w){
  TH2D * h2 =GetHist2D(name);
  h2->Fill(x,y,w);
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::IncHist(const char *name, Int_t bin){
  TH1D * h1 =GetHist1D(name);
  h1->SetBinContent(bin,h1->GetBinContent(bin)+1);
}
//###############################################################################################################
TH1 *AliAnalysisTaskHFJetIPQA::AddHistogramm(const char *name, const char *title, Int_t x, Double_t xlow, Double_t xhigh, Int_t y, Double_t ylow, Double_t yhigh){
  TObject * res = 0x0;
  res = fOutput2->FindObject(name);
  if((res)) return 0x0;

  TH1 * phist=0x0;
  if(y==0){ //TH1D*
      phist = new TH1D (name,title,x,xlow,xhigh);
    }
  else  {
      phist = new TH2D(name,title,x,xlow,xhigh,y,ylow,yhigh);
    }
  phist->Sumw2();

  fOutput2->Add(phist);
  Printf("Adding %s to output list",phist->GetName());
  return (TH1*)phist;
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::SubtractMean(Double_t val[], AliVTrack *track){
  Double_t  deltamean=fGraphMean->Eval(track->Pt() < 3. ? track->Pt() : 3 );//(fCurrentMeanFactors[0]-fCurrentMeanFactors[1]*TMath::Exp(-1*fCurrentMeanFactors[2]*(track->Pt()-fCurrentMeanFactors[3]))) *1e-4;

  val[0] -=deltamean*1e-4;;
}
//###############################################################################################################
void AliAnalysisTaskHFJetIPQA::RotateTracksAroundYAxis(double angle,double shiftx){

  //deg to rad
  //Î±(radians) = Î±(degrees) Ã— Ï€ / 180Âº
  angle = angle * TMath::Pi() /180;

  if(InputEvent()){
      // Rotate primary vertex //! dont do this since we try to recalculate the vertex
      AliESDVertex * vv = (AliESDVertex*) (InputEvent()->GetPrimaryVertex());
      /* if(vv){
            Double_t pos[3],cova[6],chi2perNDF;
            vv->GetXYZ(pos); // position
            vv->GetCovMatrix(cova); //covariance matrix
            chi2perNDF = vv->GetChi2toNDF();
            Double_t origvtxx [3]= {0};
            Double_t rotvtxx [3]= {0};

            TVector3 vtxx(pos);
            //Rotate vertex around y axis
            vtxx.GetXYZ(origvtxx);
            vtxx.RotateY(angle);
            vtxx.GetXYZ(rotvtxx);

            AliAODVertex tmpv(origvtxx,cova,chi2perNDF);
            cova[0] = tmpv.RotatedCovMatrixXX(0,angle);
            cova[1] = tmpv.RotatedCovMatrixXY(0,angle);
            cova[2] = tmpv.RotatedCovMatrixYY(0);
            cova[3] = tmpv.RotatedCovMatrixXZ(0,angle);
            cova[4] = tmpv.RotatedCovMatrixYZ(0,angle);
            cova[5] = tmpv.RotatedCovMatrixZZ(0,angle);

            AliESDVertex  tmpesd(rotvtxx, cova, vv->GetChi2() , vv->GetNContributors() , vv->GetName());
            *vv =tmpesd;

        }*/
      //Rotate Tracks
      for(int itrack = 0 ; itrack < InputEvent()->GetNumberOfTracks();++itrack){
          AliVTrack * vt = dynamic_cast<AliVTrack* > (InputEvent()->GetTrack(itrack));
          if (!vt) continue;
          double xyz[3] = {0};
          double pxpypz [3] = {0};
          vt->XvYvZv(xyz);
          vt->PxPyPz(pxpypz);
          TVector3 vp(pxpypz);
          TVector3 vx(xyz);
          //Rotate around y axis
          RotateVectorY(vp,angle);
          RotateVectorY(vx,angle);

          //update Track
          double rotxyz [3] ={0};
          double rotpxpypz [3] ={0};

          vp.GetXYZ(rotpxpypz);
          vx.GetXYZ(rotxyz);

          rotxyz[0] += shiftx;

          //Rotate cov matrix

          double cc[21] = {0};
          ((AliExternalTrackParam*)vt)->GetCovarianceXYZPxPyPz(cc);

          TMatrixD upperleft(3, 3);
          TMatrixD lowerright(3, 3);
          TMatrixD lowerleft(3, 3);

          upperleft(0,0) = cc[0];
          upperleft(0,1) = cc[1];
          upperleft(0,2) = cc[3];
          upperleft(1,0) = cc[1];
          upperleft(1,1) = cc[2];
          upperleft(1,2) = cc[4];
          upperleft(2,0) = cc[3];
          upperleft(2,1) = cc[4];
          upperleft(2,2) = cc[5];
          lowerright(0,0) = cc[9];
          lowerright(0,1) = cc[13];
          lowerright(0,2) = cc[18];
          lowerright(1,0) = cc[13];
          lowerright(1,1) = cc[14];
          lowerright(1,2) = cc[19];
          lowerright(2,0) = cc[18];
          lowerright(2,1) = cc[19];
          lowerright(2,2) = cc[20];
          lowerleft(0,0) = cc[6];
          lowerleft(0,1) = cc[10];
          lowerleft(0,2) = cc[15];
          lowerleft(1,0) = cc[7];
          lowerleft(1,1) = cc[11];
          lowerleft(1,2) = cc[16];
          lowerleft(2,0) = cc[8];
          lowerleft(2,1) = cc[12];
          lowerleft(2,2) = cc[17];

          RotateMatrixY(upperleft,angle);
          PartiallyRotateMatrixY(lowerleft,angle);
         // (lowerright,angle);

          //Repopulate cov matrix

          cc[0]=  upperleft(0,0) ;
          cc[1]=  upperleft(0,1) ;
          cc[2]=  upperleft(1,1) ;
          cc[3]=  upperleft(2,0) ;
          cc[4]=  upperleft(1,2) ;
          cc[5]=  upperleft(2,2) ;
          cc[9]  = lowerright(0,0) ;
          cc[13] = lowerright(0,1) ;
          cc[18] = lowerright(0,2);
          cc[14] = lowerright(1,1) ;
          cc[19] = lowerright(1,2) ;
          cc[20] = lowerright(2,2) ;
          cc[6]  = lowerleft(0,0) ;
          cc[10] = lowerleft(0,1)  ;
          cc[15] = lowerleft(0,2) ;
          cc[7]  = lowerleft(1,0);
          cc[11] = lowerleft(1,1) ;
          cc[16] = lowerleft(1,2);
          cc[8]  = lowerleft(2,0);
          cc[12] = lowerleft(2,1);
          cc[17] = lowerleft(2,2) ;

          Double_t ccempty[15] = {0}; AliExternalTrackParam trackparam(rotxyz,rotpxpypz,ccempty,vt->GetSign());

          trackparam.Set(rotxyz,rotpxpypz,cc,vt->GetSign());
          *((AliExternalTrackParam*)vt) = (trackparam);

        }
    }

}


void AliAnalysisTaskHFJetIPQA::RotateMatrixY(TMatrixD &m, double ang_rad)
{
  TMatrixD mroty(3,3);
  TMatrixD mroty_inverse(3,3);
  mroty(0,0) = cos(ang_rad);  mroty(0,1) = 0;mroty(0,2) = -1.*sin(ang_rad);
  mroty(1,0) = 0;  mroty(1,1) = 1;  mroty(1,2) = 0;
  mroty(2,0) = sin(ang_rad);  mroty(2,1) = 0;  mroty(2,2) = cos(ang_rad);
  mroty_inverse(0,0) = cos(ang_rad);  mroty_inverse(0,1) = 0;  mroty_inverse(0,2) = 1.*sin(ang_rad);
  mroty_inverse(1,0) = 0;  mroty_inverse(1,1) = 1;  mroty_inverse(1,2) = 0;
  mroty_inverse(2,0) = -1.* sin(ang_rad);  mroty_inverse(2,1) = 0;  mroty_inverse(2,2) = cos(ang_rad);
  m = mroty_inverse*m *mroty ;
}

void AliAnalysisTaskHFJetIPQA::PartiallyRotateMatrixY(TMatrixD &m, double ang_rad)
{
  TMatrixD mroty(3,3);
  TMatrixD mroty_inverse(3,3);
  mroty(0,0) = cos(ang_rad);  mroty(0,1) = 0;mroty(0,2) = -1.*sin(ang_rad);
  mroty(1,0) = 0;  mroty(1,1) = 1;  mroty(1,2) = 0;
  mroty(2,0) = sin(ang_rad);  mroty(2,1) = 0;  mroty(2,2) = cos(ang_rad);
  mroty_inverse(0,0) = cos(ang_rad);  mroty_inverse(0,1) = 0;  mroty_inverse(0,2) = 1.*sin(ang_rad);
  mroty_inverse(1,0) = 0;  mroty_inverse(1,1) = 1;  mroty_inverse(1,2) = 0;
  mroty_inverse(2,0) = -1.* sin(ang_rad);  mroty_inverse(2,1) = 0;  mroty_inverse(2,2) = cos(ang_rad);
  m = mroty_inverse*m  ;
}

void AliAnalysisTaskHFJetIPQA::RotateVectorY(TVector3 &v, double ang_rad)
{
  TMatrixD mroty(3,3);
  mroty(0,0) = cos(ang_rad);  mroty(0,1) = 0;mroty(0,2) = -1.*sin(ang_rad);
  mroty(1,0) = 0;  mroty(1,1) = 1;  mroty(1,2) = 0;
  mroty(2,0) = sin(ang_rad);  mroty(2,1) = 0;  mroty(2,2) = cos(ang_rad);
  v = mroty* v  ;
}
