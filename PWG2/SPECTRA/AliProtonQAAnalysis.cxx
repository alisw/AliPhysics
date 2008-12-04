/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

/* $Id: AliProtonQAAnalysis.cxx 29114 2008-10-03 16:49:02Z pchrist $ */

//-----------------------------------------------------------------
//                 AliProtonQAAnalysis class
//   This is the class to deal with the proton analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TArrayI.h>
#include <TParticle.h>

#include "AliProtonQAAnalysis.h"

#include <AliExternalTrackParam.h>
#include <AliESDEvent.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>

ClassImp(AliProtonQAAnalysis)

//____________________________________________________________________//
AliProtonQAAnalysis::AliProtonQAAnalysis() : 
  TObject(), 
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0), fMaxSigmaToVertexTPC(0),
  fMaxDCAXY(0), fMaxDCAXYTPC(0),
  fMaxDCAZ(0), fMaxDCAZTPC(0),
  fMaxConstrainChi2(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), 
  fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE), fMaxSigmaToVertexTPCFlag(kFALSE),
  fMaxDCAXYFlag(kFALSE), fMaxDCAXYTPCFlag(kFALSE),
  fMaxDCAZFlag(kFALSE), fMaxDCAZTPCFlag(kFALSE),
  fMaxConstrainChi2Flag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fESDpidFlag(kFALSE), fTPCpidFlag(kFALSE),
  fPointOnITSLayer1Flag(0), fPointOnITSLayer2Flag(0),
  fPointOnITSLayer3Flag(0), fPointOnITSLayer4Flag(0),
  fPointOnITSLayer5Flag(0), fPointOnITSLayer6Flag(0),
  fGlobalQAList(0), fQA2DList(0),
  fQAPrimaryProtonsAcceptedList(0),
  fQAPrimaryProtonsRejectedList(0),
  fQASecondaryProtonsAcceptedList(0),
  fQASecondaryProtonsRejectedList(0),
  fQAPrimaryAntiProtonsAcceptedList(0),
  fQAPrimaryAntiProtonsRejectedList(0),
  fQASecondaryAntiProtonsAcceptedList(0),
  fQASecondaryAntiProtonsRejectedList(0),
  fFunctionProbabilityFlag(kFALSE), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fUseTPCOnly(kFALSE), fUseHybridTPC(kFALSE),
  fPDGList(0), fMCProcessesList(0),
  fRunMCAnalysis(kFALSE),
  fMCProcessIdFlag(kFALSE), fMCProcessId(0),
  fMotherParticlePDGCodeFlag(kFALSE), fMotherParticlePDGCode(0),
  fAcceptedCutList(0), fRejectedCutList(0),
  fAcceptedDCAList(0), fRejectedDCAList(0),
  fRunEfficiencyAnalysis(kFALSE), fRunEfficiencyAnalysisEtaMode(kFALSE),
  fUseCutsInEfficiency(kFALSE),
  fEfficiencyList(0) {
  //Default constructor
  for(Int_t i = 0; i < 5; i++) fPartFrac[i] = 0.0;
}

//____________________________________________________________________//
AliProtonQAAnalysis::~AliProtonQAAnalysis() {
  //Default destructor
  if(fGlobalQAList) delete fGlobalQAList;
  if(fQA2DList) delete fQA2DList;
  if(fQAPrimaryProtonsAcceptedList) delete fQAPrimaryProtonsAcceptedList;
  if(fQAPrimaryProtonsRejectedList) delete fQAPrimaryProtonsRejectedList;
  if(fQASecondaryProtonsAcceptedList) delete fQASecondaryProtonsAcceptedList;
  if(fQASecondaryProtonsRejectedList) delete fQASecondaryProtonsRejectedList;
  if(fQAPrimaryAntiProtonsAcceptedList) 
    delete fQAPrimaryAntiProtonsAcceptedList;
  if(fQAPrimaryAntiProtonsRejectedList) 
    delete fQAPrimaryAntiProtonsRejectedList;
  if(fQASecondaryAntiProtonsAcceptedList) 
    delete fQASecondaryAntiProtonsAcceptedList;
  if(fQASecondaryAntiProtonsRejectedList) 
    delete fQASecondaryAntiProtonsRejectedList; 

  if(fPDGList) delete fPDGList;
  if(fMCProcessesList) delete fMCProcessesList;
  
  if(fAcceptedCutList) delete fAcceptedCutList;
  if(fRejectedCutList) delete fRejectedCutList;
  if(fAcceptedDCAList) delete fAcceptedDCAList;
  if(fRejectedDCAList) delete fRejectedDCAList;
 
  if(fEfficiencyList) delete fEfficiencyList;
}

//____________________________________________________________________//
Double_t AliProtonQAAnalysis::GetParticleFraction(Int_t i, Double_t p) {
  Double_t partFrac=0;
  if(fFunctionProbabilityFlag) {
    if(i == 0) partFrac = fElectronFunction->Eval(p);
    if(i == 1) partFrac = fMuonFunction->Eval(p);
    if(i == 2) partFrac = fPionFunction->Eval(p);
    if(i == 3) partFrac = fKaonFunction->Eval(p);
    if(i == 4) partFrac = fProtonFunction->Eval(p);
  }
  else partFrac = fPartFrac[i];

  return partFrac;
}

//____________________________________________________________________//
Bool_t AliProtonQAAnalysis::IsAccepted(AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  Double_t Pt = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

  if((fUseTPCOnly)&&(!fUseHybridTPC)) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
      dcaXY = -100.0, dcaZ = -100.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
      track->GetImpactParametersTPC(dcaXY,dcaZ);
    }
  }
  else if(fUseHybridTPC) {
     AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
      dcaXY = -100.0, dcaZ = -100.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
      track->GetImpactParameters(dcaXY,dcaZ);
    }
  }
  else{
    Pt = track->Pt();
    Px = track->Px();
    Py = track->Py();
    Pz = track->Pz();
    track->GetImpactParameters(dcaXY,dcaZ);
  }
     
  Int_t  fIdxInt[200];
  Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);

  Float_t chi2PerClusterITS = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);

  if(fPointOnITSLayer1Flag)
    if(!track->HasPointOnITSLayer(0)) return kFALSE;
  if(fPointOnITSLayer2Flag)
    if(!track->HasPointOnITSLayer(1)) return kFALSE;
  if(fPointOnITSLayer3Flag)
    if(!track->HasPointOnITSLayer(2)) return kFALSE;
  if(fPointOnITSLayer4Flag)
    if(!track->HasPointOnITSLayer(3)) return kFALSE;
  if(fPointOnITSLayer5Flag)
    if(!track->HasPointOnITSLayer(4)) return kFALSE;
  if(fPointOnITSLayer6Flag)
    if(!track->HasPointOnITSLayer(5)) return kFALSE;
  if(fMinITSClustersFlag)
    if(nClustersITS < fMinITSClusters) return kFALSE;
  if(fMaxChi2PerITSClusterFlag)
    if(chi2PerClusterITS > fMaxChi2PerITSCluster) return kFALSE; 
  if(fMinTPCClustersFlag)
    if(nClustersTPC < fMinTPCClusters) return kFALSE;
  if(fMaxChi2PerTPCClusterFlag)
    if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) return kFALSE; 
  if(fMaxCov11Flag)
    if(extCov[0] > fMaxCov11) return kFALSE;
  if(fMaxCov22Flag)
    if(extCov[2] > fMaxCov22) return kFALSE;
  if(fMaxCov33Flag)
    if(extCov[5] > fMaxCov33) return kFALSE;
  if(fMaxCov44Flag)
    if(extCov[9] > fMaxCov44) return kFALSE;
  if(fMaxCov55Flag)
    if(extCov[14] > fMaxCov55) return kFALSE;
  if(fMaxSigmaToVertexFlag)
    if(GetSigmaToVertex(track) > fMaxSigmaToVertex) return kFALSE;
  if(fMaxSigmaToVertexTPCFlag)
    if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) return kFALSE;
  if(fMaxDCAXYFlag) 
    if(TMath::Abs(dcaXY) > fMaxDCAXY) return kFALSE;
  if(fMaxDCAXYTPCFlag) 
    if(TMath::Abs(dcaXY) > fMaxDCAXYTPC) return kFALSE;
    if(fMaxDCAZFlag) 
    if(TMath::Abs(dcaZ) > fMaxDCAZ) return kFALSE;
  if(fMaxDCAZTPCFlag) 
    if(TMath::Abs(dcaZ) > fMaxDCAZTPC) return kFALSE;
  if(fMaxConstrainChi2Flag) {
    if(track->GetConstrainedChi2() > 0) 
      if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2) return kFALSE;
  }
  if(fITSRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
  if(fTPCRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
  if(fESDpidFlag)
    if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) return kFALSE;
  if(fTPCpidFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) return kFALSE;

  if((Pt < fMinPt) || (Pt > fMaxPt)) return kFALSE;
  if((Rapidity(Px,Py,Pz) < fMinY) || (Rapidity(Px,Py,Pz) > fMaxY)) 
    return kFALSE;

  return kTRUE;
}

//____________________________________________________________________//
void AliProtonQAAnalysis::FillQA(AliESDtrack* track, AliStack *stack) {
  // Checks if the track is excluded from the cuts
  Int_t nPrimaries = stack->GetNprimary();
  Int_t label = TMath::Abs(track->GetLabel());

  Double_t Pt = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;
  if((fUseTPCOnly)&&(!fUseHybridTPC)) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
      dcaXY = -100.0, dcaZ = -100.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
      track->GetImpactParametersTPC(dcaXY,dcaZ);
    }
  }
  else if(fUseHybridTPC) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
      dcaXY = -100.0, dcaZ = -100.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
      track->GetImpactParameters(dcaXY,dcaZ);
    }
  }
  else{
    Pt = track->Pt();
    Px = track->Px();
    Py = track->Py();
    Pz = track->Pz();
    track->GetImpactParameters(dcaXY,dcaZ);
  }
     
  Int_t  fIdxInt[200];
  Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);

  Float_t chi2PerClusterITS = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);
  
  //cout<<"Charge: "<<track->Charge()<<
  //" - Label/Primaries: "<<label<<"/"<<nPrimaries<<
  //" - TPC clusters: "<<nClustersTPC<<endl;
  //protons
  if(track->Charge() > 0) {
    //Primaries
    if(label <= nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  //status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  //status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  //cout<<"Primary proton rejected"<<endl;
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  //status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters) {
	  //cout<<"Primary proton accepted"<<endl;
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  //status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  //status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  //status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  //status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  //status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  //status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fMaxDCAXYFlag) {
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(11)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(11)))->Fill(dcaXY);
      }//DCA xy global tracking
      if(fMaxDCAXYTPCFlag) {
	if(dcaXY > fMaxDCAXYTPC) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXYTPC)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(12)))->Fill(dcaXY);
      }//DCA xy TPC tracking
      if(fMaxDCAZFlag) {
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(13)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(13)))->Fill(dcaZ);
      }//DCA z global tracking
      if(fMaxDCAZTPCFlag) {
	if(dcaZ > fMaxDCAZTPC) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZTPC)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(14)))->Fill(dcaZ);
      }//DCA z TPC tracking
      if(fMaxConstrainChi2Flag) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    //status = kFALSE;
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fMaxConstrainChi2)
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(16)))->Fill(0);
	//status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(17)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(18)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(19)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fPointOnITSLayer1Flag) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fPointOnITSLayer2Flag) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fPointOnITSLayer3Flag) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fPointOnITSLayer4Flag) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fPointOnITSLayer5Flag) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fPointOnITSLayer6Flag) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
    }//primary particle cut

    //Secondaries
    if(label > nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  //status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  //status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  //cout<<"Secondary proton rejected"<<endl;
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  //status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters) {
	  //cout<<"Secondary proton accepted"<<endl;
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  //status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  //status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  //status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  //status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  //status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  //status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fMaxDCAXYFlag) {
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(11)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(11)))->Fill(dcaXY);
      }//DCA xy global tracking
      if(fMaxDCAXYTPCFlag) {
	if(dcaXY > fMaxDCAXYTPC) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXYTPC)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(12)))->Fill(dcaXY);
      }//DCA xy TPC tracking
      if(fMaxDCAZFlag) {
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(13)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(13)))->Fill(dcaZ);
      }//DCA z global tracking
      if(fMaxDCAZTPCFlag) {
	if(dcaZ > fMaxDCAZTPC) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZTPC)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(14)))->Fill(dcaZ);
      }//DCA z TPC tracking
      if(fMaxConstrainChi2Flag) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    //status = kFALSE;
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fMaxConstrainChi2)
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(16)))->Fill(0);
	//status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(17)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(18)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(19)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fPointOnITSLayer1Flag) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fPointOnITSLayer2Flag) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fPointOnITSLayer3Flag) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fPointOnITSLayer4Flag) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fPointOnITSLayer5Flag) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fPointOnITSLayer6Flag) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
    }//secondary particle cut
  }//protons

  //antiprotons
  if(track->Charge() < 0) {
    //Primaries
    if(label <= nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  //status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  //status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  //cout<<"Primary antiproton rejected"<<endl;
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  //status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters) {
	  //cout<<"Primary antiproton accepted"<<endl;
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  //status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  //status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  //status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  //status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  //status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  //status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fMaxDCAXYFlag) {
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(11)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(11)))->Fill(dcaXY);
      }//DCA xy global tracking
      if(fMaxDCAXYTPCFlag) {
	if(dcaXY > fMaxDCAXYTPC) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXYTPC)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(12)))->Fill(dcaXY);
      }//DCA xy TPC tracking
      if(fMaxDCAZFlag) {
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(13)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(13)))->Fill(dcaZ);
      }//DCA z global tracking
      if(fMaxDCAZTPCFlag) {
	if(dcaZ > fMaxDCAZTPC) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZTPC)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(14)))->Fill(dcaZ);
      }//DCA z TPC tracking
      if(fMaxConstrainChi2Flag) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2) {
	    ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    //status = kFALSE;
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fMaxConstrainChi2)
	    ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(16)))->Fill(0);
	//status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(17)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(18)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(19)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fPointOnITSLayer1Flag) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fPointOnITSLayer2Flag) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fPointOnITSLayer3Flag) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fPointOnITSLayer4Flag) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fPointOnITSLayer5Flag) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fPointOnITSLayer6Flag) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
    }//primary particle cut

    //Secondaries
    if(label > nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  //status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  //status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  //cout<<"Secondary antiproton rejected"<<endl;
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  //status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters) {
	  //cout<<"Secondary antiproton accepted"<<endl;
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  //status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  //status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  //status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  //status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  //status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  //status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  //status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fMaxDCAXYFlag) {
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(11)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(11)))->Fill(dcaXY);
      }//DCA xy global tracking
      if(fMaxDCAXYTPCFlag) {
	if(dcaXY > fMaxDCAXYTPC) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXYTPC)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(12)))->Fill(dcaXY);
      }//DCA xy TPC tracking
      if(fMaxDCAZFlag) {
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(13)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(13)))->Fill(dcaZ);
      }//DCA z global tracking
      if(fMaxDCAZTPCFlag) {
	if(dcaZ > fMaxDCAZTPC) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZTPC)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(14)))->Fill(dcaZ);
      }//DCA z TPC tracking
      if(fMaxConstrainChi2Flag) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2) {
	    ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    //status = kFALSE;
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fMaxConstrainChi2)
	    ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(16)))->Fill(0);
	//status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(17)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(18)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(19)))->Fill(0);
	  //status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fPointOnITSLayer1Flag) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fPointOnITSLayer2Flag) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fPointOnITSLayer3Flag) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fPointOnITSLayer4Flag) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fPointOnITSLayer5Flag) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fPointOnITSLayer6Flag) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
    }//secondary particle cut
  }//antiprotons
}

//____________________________________________________________________//
Float_t AliProtonQAAnalysis::GetSigmaToVertex(AliESDtrack* esdTrack) {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  if((fUseTPCOnly)&&(!fUseHybridTPC))
    esdTrack->GetImpactParametersTPC(b,bCov);
  else
    esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    //AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
}

//____________________________________________________________________//
Double_t AliProtonQAAnalysis::Rapidity(Double_t Px, Double_t Py, Double_t Pz) {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  
  Double_t P = TMath::Sqrt(TMath::Power(Px,2) + 
                           TMath::Power(Py,2) + 
			   TMath::Power(Pz,2));
  Double_t energy = TMath::Sqrt(P*P + fMass*fMass);
  Double_t y = -999;
  if(energy != Pz) 
    y = 0.5*TMath::Log((energy + Pz)/(energy - Pz));

  return y;
}

//____________________________________________________________________//
void AliProtonQAAnalysis::SetRunQAAnalysis() {
  //initializes the QA lists
  //fQAHistograms = kTRUE;
  fGlobalQAList = new TList();
  fQA2DList = new TList();
  fQA2DList->SetName("fQA2DList");
  fGlobalQAList->Add(fQA2DList);
  
  fQAPrimaryProtonsAcceptedList = new TList();
  fQAPrimaryProtonsAcceptedList->SetName("fQAPrimaryProtonsAcceptedList");
  fGlobalQAList->Add(fQAPrimaryProtonsAcceptedList);
  
  fQAPrimaryProtonsRejectedList = new TList();
  fQAPrimaryProtonsRejectedList->SetName("fQAPrimaryProtonsRejectedList");
  fGlobalQAList->Add(fQAPrimaryProtonsRejectedList);
  
  fQASecondaryProtonsAcceptedList = new TList();
  fQASecondaryProtonsAcceptedList->SetName("fQASecondaryProtonsAcceptedList");
  fGlobalQAList->Add(fQASecondaryProtonsAcceptedList);
  
  fQASecondaryProtonsRejectedList = new TList();
  fQASecondaryProtonsRejectedList->SetName("fQASecondaryProtonsRejectedList");
  fGlobalQAList->Add(fQASecondaryProtonsRejectedList);
  
  fQAPrimaryAntiProtonsAcceptedList = new TList();
  fQAPrimaryAntiProtonsAcceptedList->SetName("fQAPrimaryAntiProtonsAcceptedList");
  fGlobalQAList->Add(fQAPrimaryAntiProtonsAcceptedList);
  
  fQAPrimaryAntiProtonsRejectedList = new TList();
  fQAPrimaryAntiProtonsRejectedList->SetName("fQAPrimaryAntiProtonsRejectedList");
  fGlobalQAList->Add(fQAPrimaryAntiProtonsRejectedList);
  
  fQASecondaryAntiProtonsAcceptedList = new TList();
  fQASecondaryAntiProtonsAcceptedList->SetName("fQASecondaryAntiProtonsAcceptedList");
  fGlobalQAList->Add(fQASecondaryAntiProtonsAcceptedList);
  
  fQASecondaryAntiProtonsRejectedList = new TList();
  fQASecondaryAntiProtonsRejectedList->SetName("fQASecondaryAntiProtonsRejectedList");
  fGlobalQAList->Add(fQASecondaryAntiProtonsRejectedList);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::SetQAYPtBins(Int_t nbinsY, Double_t minY, Double_t maxY,
				      Int_t nbinsPt, Double_t minPt, Double_t maxPt) {
  //Initializes the QA binning
  fNBinsY = nbinsY;
  fMinY = minY; fMaxY = maxY;
  fNBinsPt = nbinsPt;
  fMinPt = minPt; fMaxPt = maxPt;
  InitQA();
  InitCutLists();
  if(fRunMCAnalysis) InitMCAnalysis();
  if(fRunEfficiencyAnalysis) InitEfficiencyAnalysis();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitEfficiencyAnalysis() {
  //Initialization of the efficiency list - reconstruction & PID efficiency
  //Adding each monitored object in the list
  fEfficiencyList = new TList();

  //MC primary protons and antiprotons for the reconstruction efficiency
  TH2D *gHistMCYPtProtons = new TH2D("gHistMCYPtProtons",
				     ";;P_{T} [GeV/c]",
				     fNBinsY,fMinY,fMaxY,
				     fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtons->GetXaxis()->SetTitle("y");
  gHistMCYPtProtons->SetStats(kTRUE);
  gHistMCYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtons);
  TH2D *gHistMCYPtAntiProtons = new TH2D("gHistMCYPtAntiProtons",
					 ";y;P_{T} [GeV/c]",
					 fNBinsY,fMinY,fMaxY,
					 fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtons->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtons->SetStats(kTRUE);
  gHistMCYPtAntiProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtons);

  //MC secondary protons and antiprotons that come from weak decay for the reconstruction efficiency
  TH2D *gHistMCYPtProtonsFromWeak = new TH2D("gHistMCYPtProtonsFromWeak",
					     ";;P_{T} [GeV/c]",
					     fNBinsY,fMinY,fMaxY,
					     fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistMCYPtProtonsFromWeak->SetStats(kTRUE);
  gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtonsFromWeak);
  TH2D *gHistMCYPtAntiProtonsFromWeak = new TH2D("gHistMCYPtAntiProtonsFromWeak",
						 ";y;P_{T} [GeV/c]",
						 fNBinsY,fMinY,fMaxY,
						 fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtonsFromWeak->SetStats(kTRUE);
  gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtonsFromWeak);

  //MC secondary protons and antiprotons that come from hadronic interactions for the reconstruction efficiency
  TH2D *gHistMCYPtProtonsFromHadronic = new TH2D("gHistMCYPtProtonsFromHadronic",
						 ";;P_{T} [GeV/c]",
						 fNBinsY,fMinY,fMaxY,
						 fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistMCYPtProtonsFromHadronic->SetStats(kTRUE);
  gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtonsFromHadronic);
  TH2D *gHistMCYPtAntiProtonsFromHadronic = new TH2D("gHistMCYPtAntiProtonsFromHadronic",
						     ";y;P_{T} [GeV/c]",
						     fNBinsY,fMinY,fMaxY,
						     fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtonsFromHadronic->SetStats(kTRUE);
  gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtonsFromHadronic);
  
  //ESD primary protons and antiprotons for the reconstruction efficiency
  TH2D *gHistESDYPtProtons = new TH2D("gHistESDYPtProtons",
				      ";;P_{T} [GeV/c]",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDYPtProtons->SetStats(kTRUE);
  gHistESDYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtons);
  TH2D *gHistESDYPtAntiProtons = new TH2D("gHistESDYPtAntiProtons",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtons->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtons->SetStats(kTRUE);
  gHistESDYPtAntiProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtons);

  //ESD (anti)protons from weak decays for the reconstruction efficiency
  TH2D *gHistESDYPtProtonsFromWeak = new TH2D("gHistESDYPtProtonsFromWeak",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsFromWeak->SetStats(kTRUE);
  gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtonsFromWeak);
  TH2D *gHistESDYPtAntiProtonsFromWeak = new TH2D("gHistESDYPtAntiProtonsFromWeak",
						  ";;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsFromWeak->SetStats(kTRUE);
  gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtonsFromWeak);

  //ESD (anti)protons from hadronic interactions for the reconstruction efficiency
  TH2D *gHistESDYPtProtonsFromHadronic = new TH2D("gHistESDYPtProtonsFromHadronic",
						  ";;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsFromHadronic->SetStats(kTRUE);
  gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtonsFromHadronic);
  TH2D *gHistESDYPtAntiProtonsFromHadronic = new TH2D("gHistESDYPtAntiProtonsFromHadronic",
						      ";;P_{T} [GeV/c]",
						      fNBinsY,fMinY,fMaxY,
						      fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsFromHadronic->SetStats(kTRUE);
  gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtonsFromHadronic);
  
  
  //ESD reconstructed tracks that were initially protons for the PID efficiency
  TH2D *gHistESDInitYPtProtons = new TH2D("gHistESDInitYPtProtons",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDInitYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDInitYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDInitYPtProtons->SetStats(kTRUE);
  gHistESDInitYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDInitYPtProtons);
  
  //ESD reconstructed tracks that were initially protons and were identified as protons for the PID efficiency
  TH2D *gHistESDIdYPtProtons = new TH2D("gHistESDIdYPtProtons",
					";;P_{T} [GeV/c]",
					fNBinsY,fMinY,fMaxY,
					fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDIdYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDIdYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDIdYPtProtons->SetStats(kTRUE);
  gHistESDIdYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDIdYPtProtons);
 
  //ESD reconstructed tracks that were identified as protons for the PID contamination
  TH2D *gHistESDRecIdYPtProtons = new TH2D("gHistESDRecIdYPtProtons",
					   ";;P_{T} [GeV/c]",
					   fNBinsY,fMinY,fMaxY,
					   fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDRecIdYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDRecIdYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDRecIdYPtProtons->SetStats(kTRUE);
  gHistESDRecIdYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDRecIdYPtProtons);

  //ESD reconstructed tracks that were identified as protons but were initially not protons for the PID contamination
  TH2D *gHistESDContamYPtProtons = new TH2D("gHistESDContamYPtProtons",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fMinY,fMaxY,
					    fNBinsPt,fMinPt,fMaxPt);
  if(fRunEfficiencyAnalysisEtaMode) 
    gHistESDContamYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDContamYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDContamYPtProtons->SetStats(kTRUE);
  gHistESDContamYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDContamYPtProtons);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitCutLists() {
  //Initialization of the cut lists
  //Adding each monitored object in each list

  //Accepted cut list
  fAcceptedCutList = new TList();
  TH1F *gPrimaryProtonsClustersOnITSLayers = new TH1F("gPrimaryProtonsClustersOnITSLayers",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gPrimaryProtonsClustersOnITSLayers);
  TH1F *gPrimaryAntiProtonsClustersOnITSLayers = new TH1F("gPrimaryAntiProtonsClustersOnITSLayers",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gPrimaryAntiProtonsClustersOnITSLayers);
  TH1F *gSecondaryProtonsClustersOnITSLayers = new TH1F("gSecondaryProtonsClustersOnITSLayers",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gSecondaryProtonsClustersOnITSLayers);
  TH1F *gSecondaryAntiProtonsClustersOnITSLayers = new TH1F("gSecondaryAntiProtonsClustersOnITSLayers",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gSecondaryAntiProtonsClustersOnITSLayers);

  TH1F *gPrimaryProtonsNClustersITS = new TH1F("gPrimaryProtonsNClustersITS",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gPrimaryProtonsNClustersITS);
  TH1F *gPrimaryAntiProtonsNClustersITS = new TH1F("gPrimaryAntiProtonsNClustersITS",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gPrimaryAntiProtonsNClustersITS);
  TH1F *gSecondaryProtonsNClustersITS = new TH1F("gSecondaryProtonsNClustersITS",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gSecondaryProtonsNClustersITS);
  TH1F *gSecondaryAntiProtonsNClustersITS = new TH1F("gSecondaryAntiProtonsNClustersITS",";ITS Layer;Entries",6,0.5,6.5);
  fAcceptedCutList->Add(gSecondaryAntiProtonsNClustersITS);

  TH1F *gPrimaryProtonsChi2PerClusterITS = new TH1F("gPrimaryProtonsChi2PerClusterITS",
						    ";x^{2}/N_{clusters} (ITS);Entries",
						    100,0,20);
  fAcceptedCutList->Add(gPrimaryProtonsChi2PerClusterITS);
  TH1F *gPrimaryAntiProtonsChi2PerClusterITS = new TH1F("gPrimaryAntiProtonsChi2PerClusterITS",
							";x^{2}/N_{clusters} (ITS);Entries",
							100,0,20);
  fAcceptedCutList->Add(gPrimaryAntiProtonsChi2PerClusterITS);
  TH1F *gSecondaryProtonsChi2PerClusterITS = new TH1F("gSecondaryProtonsChi2PerClusterITS",
						      ";x^{2}/N_{clusters} (ITS);Entries",
						      100,0,20);
  fAcceptedCutList->Add(gSecondaryProtonsChi2PerClusterITS);
  TH1F *gSecondaryAntiProtonsChi2PerClusterITS = new TH1F("gSecondaryAntiProtonsChi2PerClusterITS",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,20);
  fAcceptedCutList->Add(gSecondaryAntiProtonsChi2PerClusterITS);

  TH1F *gPrimaryProtonsConstrainChi2 = new TH1F("gPrimaryProtonsConstrainChi2",
						";Log_{10}(#chi^{2});Entries",
						100,-10,10);
  fAcceptedCutList->Add(gPrimaryProtonsConstrainChi2);
  TH1F *gPrimaryAntiProtonsConstrainChi2 = new TH1F("gPrimaryAntiProtonsConstrainChi2",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fAcceptedCutList->Add(gPrimaryAntiProtonsConstrainChi2);
  TH1F *gSecondaryProtonsConstrainChi2 = new TH1F("gSecondaryProtonsConstrainChi2",
						  ";Log_{10}(#chi^{2});Entries",
						  100,-10,10);
  fAcceptedCutList->Add(gSecondaryProtonsConstrainChi2);
  TH1F *gSecondaryAntiProtonsConstrainChi2 = new TH1F("gSecondaryAntiProtonsConstrainChi2",
						      ";Log_{10}(#chi^{2});Entries",
						      100,-10,10);
  fAcceptedCutList->Add(gSecondaryAntiProtonsConstrainChi2);

  TH1F *gPrimaryProtonsTPCClusters = new TH1F("gPrimaryProtonsTPCClusters",
					      ";N_{clusters} (TPC);Entries",
					      100,0,200);
  fAcceptedCutList->Add(gPrimaryProtonsTPCClusters);
  TH1F *gPrimaryAntiProtonsTPCClusters = new TH1F("gPrimaryAntiProtonsTPCClusters",
						  ";N_{clusters} (TPC);Entries",
						  100,0,200);
  fAcceptedCutList->Add(gPrimaryAntiProtonsTPCClusters);
  TH1F *gSecondaryProtonsTPCClusters = new TH1F("gSecondaryProtonsTPCClusters",
						";N_{clusters} (TPC);Entries",
						100,0,200);
  fAcceptedCutList->Add(gSecondaryProtonsTPCClusters);
  TH1F *gSecondaryAntiProtonsTPCClusters = new TH1F("gSecondaryAntiProtonsTPCClusters",
						    ";N_{clusters} (TPC);Entries",
						    100,0,200);
  fAcceptedCutList->Add(gSecondaryAntiProtonsTPCClusters);

  TH1F *gPrimaryProtonsChi2PerClusterTPC = new TH1F("gPrimaryProtonsChi2PerClusterTPC",
						    ";x^{2}/N_{clusters} (TPC);Entries",
						    100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsChi2PerClusterTPC);
  TH1F *gPrimaryAntiProtonsChi2PerClusterTPC = new TH1F("gPrimaryAntiProtonsChi2PerClusterTPC",
							";x^{2}/N_{clusters} (TPC);Entries",
							100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsChi2PerClusterTPC);
  TH1F *gSecondaryProtonsChi2PerClusterTPC = new TH1F("gSecondaryProtonsChi2PerClusterTPC",
						      ";x^{2}/N_{clusters} (TPC);Entries",
						      100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsChi2PerClusterTPC);
  TH1F *gSecondaryAntiProtonsChi2PerClusterTPC = new TH1F("gSecondaryAntiProtonsChi2PerClusterTPC",
							  ";x^{2}/N_{clusters} (TPC);Entries",
							  100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsChi2PerClusterTPC);

  TH1F *gPrimaryProtonsExtCov11 = new TH1F("gPrimaryProtonsExtCov11",
					   ";#sigma_{y} [cm];Entries",
					   100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsExtCov11);
  TH1F *gPrimaryAntiProtonsExtCov11 = new TH1F("gPrimaryAntiProtonsExtCov11",
					       ";#sigma_{y} [cm];Entries",
					       100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsExtCov11);
  TH1F *gSecondaryProtonsExtCov11 = new TH1F("gSecondaryProtonsExtCov11",
					     ";#sigma_{y} [cm];Entries",
					     100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsExtCov11);
  TH1F *gSecondaryAntiProtonsExtCov11 = new TH1F("gSecondaryAntiProtonsExtCov11",
						 ";#sigma_{y} [cm];Entries",
						 100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsExtCov11);


  TH1F *gPrimaryProtonsExtCov22 = new TH1F("gPrimaryProtonsExtCov22",
					   ";#sigma_{z} [cm];Entries",
					   100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsExtCov22);
  TH1F *gPrimaryAntiProtonsExtCov22 = new TH1F("gPrimaryAntiProtonsExtCov22",
					       ";#sigma_{z} [cm];Entries",
					       100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsExtCov22);
  TH1F *gSecondaryProtonsExtCov22 = new TH1F("gSecondaryProtonsExtCov22",
					     ";#sigma_{z} [cm];Entries",
					     100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsExtCov22);
  TH1F *gSecondaryAntiProtonsExtCov22 = new TH1F("gSecondaryAntiProtonsExtCov22",
						 ";#sigma_{z} [cm];Entries",
						 100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsExtCov22);


  TH1F *gPrimaryProtonsExtCov33 = new TH1F("gPrimaryProtonsExtCov33",
					   ";#sigma_{sin(#phi)};Entries",
					   100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsExtCov33);
  TH1F *gPrimaryAntiProtonsExtCov33 = new TH1F("gPrimaryAntiProtonsExtCov33",
					       ";#sigma_{sin(#phi)};Entries",
					       100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsExtCov33);
  TH1F *gSecondaryProtonsExtCov33 = new TH1F("gSecondaryProtonsExtCov33",
					     ";#sigma_{sin(#phi)};Entries",
					     100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsExtCov33);
  TH1F *gSecondaryAntiProtonsExtCov33 = new TH1F("gSecondaryAntiProtonsExtCov33",
						 ";#sigma_{sin(#phi)};Entries",
						 100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsExtCov33);


  TH1F *gPrimaryProtonsExtCov44 = new TH1F("gPrimaryProtonsExtCov44",
					   ";#sigma_{tan(#lambda)};Entries",
					   100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsExtCov44);
  TH1F *gPrimaryAntiProtonsExtCov44 = new TH1F("gPrimaryAntiProtonsExtCov44",
					       ";#sigma_{tan(#lambda)};Entries",
					       100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsExtCov44);
  TH1F *gSecondaryProtonsExtCov44 = new TH1F("gSecondaryProtonsExtCov44",
					     ";#sigma_{tan(#lambda)};Entries",
					     100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsExtCov44);
  TH1F *gSecondaryAntiProtonsExtCov44 = new TH1F("gSecondaryAntiProtonsExtCov44",
						 ";#sigma_{tan(#lambda)};Entries",
						 100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsExtCov44);


  TH1F *gPrimaryProtonsExtCov55 = new TH1F("gPrimaryProtonsExtCov55",
					   ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					   100,0,4);
  fAcceptedCutList->Add(gPrimaryProtonsExtCov55);
  TH1F *gPrimaryAntiProtonsExtCov55 = new TH1F("gPrimaryAntiProtonsExtCov55",
					       ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					       100,0,4);
  fAcceptedCutList->Add(gPrimaryAntiProtonsExtCov55);
  TH1F *gSecondaryProtonsExtCov55 = new TH1F("gSecondaryProtonsExtCov55",
					     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					     100,0,4);
  fAcceptedCutList->Add(gSecondaryProtonsExtCov55);
  TH1F *gSecondaryAntiProtonsExtCov55 = new TH1F("gSecondaryAntiProtonsExtCov55",
						 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						 100,0,4);
  fAcceptedCutList->Add(gSecondaryAntiProtonsExtCov55);

  //DCA list
  fAcceptedDCAList = new TList();
  TH1F *gPrimaryProtonsDCAXY = new TH1F("gPrimaryProtonsDCAXY",
					";DCA_{xy} [cm];Entries",
					100,0,20);
  fAcceptedDCAList->Add(gPrimaryProtonsDCAXY);
  TH1F *gPrimaryAntiProtonsDCAXY = new TH1F("gPrimaryAntiProtonsDCAXY",
					    ";DCA_{xy} [cm];Entries",
					    100,0,20);
  fAcceptedDCAList->Add(gPrimaryAntiProtonsDCAXY);
  TH1F *gSecondaryProtonsDCAXY = new TH1F("gSecondaryProtonsDCAXY",
					  ";DCA_{xy} [cm];Entries",
					  100,0,20);
  fAcceptedDCAList->Add(gSecondaryProtonsDCAXY);
  TH1F *gSecondaryAntiProtonsDCAXY = new TH1F("gSecondaryAntiProtonsDCAXY",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);

  fAcceptedDCAList->Add(gSecondaryAntiProtonsDCAXY);
  TH1F *gPrimaryProtonsDCAZ = new TH1F("gPrimaryProtonsDCAZ",
				       ";DCA_{z} [cm];Entries",
				       100,0,20);
  fAcceptedDCAList->Add(gPrimaryProtonsDCAZ);
  TH1F *gPrimaryAntiProtonsDCAZ = new TH1F("gPrimaryAntiProtonsDCAZ",
					   ";DCA_{z} [cm];Entries",
					   100,0,20);
  fAcceptedDCAList->Add(gPrimaryAntiProtonsDCAZ);
  TH1F *gSecondaryProtonsDCAZ = new TH1F("gSecondaryProtonsDCAZ",
					 ";DCA_{z} [cm];Entries",
					 100,0,20);
  fAcceptedDCAList->Add(gSecondaryProtonsDCAZ);
  TH1F *gSecondaryAntiProtonsDCAZ = new TH1F("gSecondaryAntiProtonsDCAZ",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  fAcceptedDCAList->Add(gSecondaryAntiProtonsDCAZ);

  TH1F *gPrimaryProtonsSigmaToVertex = new TH1F("gPrimaryProtonsSigmaToVertex",
						";#sigma_{Vertex};Entries",
						100,0,10);
  fAcceptedDCAList->Add(gPrimaryProtonsSigmaToVertex);
  TH1F *gPrimaryAntiProtonsSigmaToVertex = new TH1F("gPrimaryAntiProtonsSigmaToVertex",
						    ";#sigma_{Vertex};Entries",
						    100,0,10);
  fAcceptedDCAList->Add(gPrimaryAntiProtonsSigmaToVertex);
  TH1F *gSecondaryProtonsSigmaToVertex = new TH1F("gSecondaryProtonsSigmaToVertex",
						  ";#sigma_{Vertex};Entries",
						  100,0,10);
  fAcceptedDCAList->Add(gSecondaryProtonsSigmaToVertex);
  TH1F *gSecondaryAntiProtonsSigmaToVertex = new TH1F("gSecondaryAntiProtonsSigmaToVertex",
						      ";#sigma_{Vertex};Entries",
						      100,0,10);
  fAcceptedDCAList->Add(gSecondaryAntiProtonsSigmaToVertex);

}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitQA() {
  //Initializes the QA histograms and builds the directory structure
  //if(!fQAHistograms) 
  SetRunQAAnalysis();

  //2D histograms
  //TDirectory *dir2D = gDirectory->mkdir("2D");
  //fGlobalQAList->Add(dir2D); dir2D->cd();
  TH2D *gHistYPtPrimaryProtonsPass = new TH2D("gHistYPtPrimaryProtonsPass",
					      ";y;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsPass);//y-pT of primary accepted ESD protons
  TH2D *gHistYPtPrimaryProtonsReject = new TH2D("gHistYPtPrimaryProtonsReject",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsReject);//y-pT of primary rejected ESD protons

  TH2D *gHistYPtSecondaryProtonsPass = new TH2D("gHistYPtSecondaryProtonsPass",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsPass);//y-pT of secondary accepted ESD protons
  TH2D *gHistYPtSecondaryProtonsReject = new TH2D("gHistYPtSecondaryProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsReject);//y-pT of secondary rejected ESD protons

  TH2D *gHistYPtPrimaryAntiProtonsPass = new TH2D("gHistYPtPrimaryAntiProtonsPass",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsPass);//y-pT of primary accepted ESD antiprotons
  TH2D *gHistYPtPrimaryAntiProtonsReject = new TH2D("gHistYPtPrimaryAntiProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsReject);//y-pT of primary rejected ESD antiprotons

  TH2D *gHistYPtSecondaryAntiProtonsPass = new TH2D("gHistYPtSecondaryAntiProtonsPass",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsPass);//y-pT of secondary accepted ESD antiprotons
  TH2D *gHistYPtSecondaryAntiProtonsReject = new TH2D("gHistYPtSecondaryAntiProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsReject);//y-pT of secondary rejected ESD antiprotons

  TH2D *gHistYPtPrimaryProtonsMC = new TH2D("gHistYPtPrimaryProtonsMC",
					    ";y;P_{T} [GeV/c]",
					    fNBinsY,fMinY,fMaxY,
					    fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsMC->SetStats(kTRUE);
  gHistYPtPrimaryProtonsMC->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsMC);//y-pT of primary MC protons
  TH2D *gHistYPtPrimaryAntiProtonsMC = new TH2D("gHistYPtPrimaryAntiProtonsMC",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsMC->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsMC->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsMC);//y-pT of primary MC antiprotons

  TH3F *gHistYPtPDGProtonsPass = new TH3F("gHistYPtPDGProtonsPass",
					  ";y;P_{T} [GeV/c];PDG",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt,
					  14,-0.5,13.5);
  fQA2DList->Add(gHistYPtPDGProtonsPass);//composition of secondary protons
  TH3F *gHistYPtPDGAntiProtonsPass = new TH3F("gHistYPtPDGAntiProtonsPass",
					      ";y;P_{T} [GeV/c];PDG",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt,
					      14,-0.5,13.5);
  fQA2DList->Add(gHistYPtPDGAntiProtonsPass);//composition of secondary antiprotons

  /*gDirectory->cd("../");
  //protons
  TDirectory *dirProtons = gDirectory->mkdir("Protons");
  fGlobalQAList->Add(dirProtons); dirProtons->cd();*/
  
  //________________________________________________________________//
  /*TDirectory *dirProtonsPrimary = gDirectory->mkdir("Primaries");
  dirProtonsPrimary->cd();
  TDirectory *dirProtonsPrimaryAccepted = gDirectory->mkdir("Accepted");
  dirProtonsPrimaryAccepted->cd();*/

  //Accepted primary protons
  TH1F *gPrimaryProtonsITSClustersPass = new TH1F("gPrimaryProtonsITSClustersPass",
					    ";N_{clusters} (ITS);Entries",
					    7,0,7);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsITSClustersPass);
  TH1F *gPrimaryProtonsChi2PerClusterITSPass = new TH1F("gPrimaryProtonsChi2PerClusterITSPass",
						  ";x^{2}/N_{clusters} (ITS);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsChi2PerClusterITSPass);
  TH1F *gPrimaryProtonsTPCClustersPass = new TH1F("gPrimaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsTPCClustersPass);
  TH1F *gPrimaryProtonsChi2PerClusterTPCPass = new TH1F("gPrimaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsChi2PerClusterTPCPass);
  TH1F *gPrimaryProtonsExtCov11Pass = new TH1F("gPrimaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsExtCov11Pass);
  TH1F *gPrimaryProtonsExtCov22Pass = new TH1F("gPrimaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsExtCov22Pass);
  TH1F *gPrimaryProtonsExtCov33Pass = new TH1F("gPrimaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsExtCov33Pass);
  TH1F *gPrimaryProtonsExtCov44Pass = new TH1F("gPrimaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsExtCov44Pass);
  TH1F *gPrimaryProtonsExtCov55Pass = new TH1F("gPrimaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsExtCov55Pass);
  TH1F *gPrimaryProtonsSigmaToVertexPass = new TH1F("gPrimaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsSigmaToVertexPass);
  TH1F *gPrimaryProtonsSigmaToVertexTPCPass = new TH1F("gPrimaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsSigmaToVertexTPCPass);
  TH1F *gPrimaryProtonsDCAXYPass = new TH1F("gPrimaryProtonsDCAXYPass",
					     ";DCA_{xy} [cm];Entries",
					     100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsDCAXYPass);
  TH1F *gPrimaryProtonsDCAXYTPCPass = new TH1F("gPrimaryProtonsDCAXYTPCPass",
					       ";DCA_{xy} [cm];Entries",
					       100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsDCAXYTPCPass);
  TH1F *gPrimaryProtonsDCAZPass = new TH1F("gPrimaryProtonsDCAZPass",
					   ";DCA_{z} [cm];Entries",
					   100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsDCAZPass);
  TH1F *gPrimaryProtonsDCAZTPCPass = new TH1F("gPrimaryProtonsDCAZTPCPass",
					      ";DCA_{z} [cm];Entries",
					      100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsDCAZTPCPass);
  TH1F *gPrimaryProtonsConstrainChi2Pass = new TH1F("gPrimaryProtonsConstrainChi2Pass",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsConstrainChi2Pass);
  TH1F *gPrimaryProtonsITSRefitPass = new TH1F("gPrimaryProtonsITSRefitPass",
					       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsITSRefitPass);
  TH1F *gPrimaryProtonsTPCRefitPass = new TH1F("gPrimaryProtonsTPCRefitPass",
					       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsTPCRefitPass);
  TH1F *gPrimaryProtonsESDpidPass = new TH1F("gPrimaryProtonsESDpidPass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsESDpidPass);
  TH1F *gPrimaryProtonsTPCpidPass = new TH1F("gPrimaryProtonsTPCpidPass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsTPCpidPass);
  TH1F *gPrimaryProtonsPointOnITSLayer1Pass = new TH1F("gPrimaryProtonsPointOnITSLayer1Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer1Pass);
  TH1F *gPrimaryProtonsPointOnITSLayer2Pass = new TH1F("gPrimaryProtonsPointOnITSLayer2Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer2Pass);
  TH1F *gPrimaryProtonsPointOnITSLayer3Pass = new TH1F("gPrimaryProtonsPointOnITSLayer3Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer3Pass);
  TH1F *gPrimaryProtonsPointOnITSLayer4Pass = new TH1F("gPrimaryProtonsPointOnITSLayer4Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer4Pass);
  TH1F *gPrimaryProtonsPointOnITSLayer5Pass = new TH1F("gPrimaryProtonsPointOnITSLayer5Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer5Pass);
  TH1F *gPrimaryProtonsPointOnITSLayer6Pass = new TH1F("gPrimaryProtonsPointOnITSLayer6Pass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsPointOnITSLayer6Pass);

  //Rejected primary protons
  /*gDirectory->cd("../");
  TDirectory *dirProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsPrimaryRejected->cd();*/

  TH1F *gPrimaryProtonsITSClustersReject = new TH1F("gPrimaryProtonsITSClustersReject",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsITSClustersReject);
  TH1F *gPrimaryProtonsChi2PerClusterITSReject = new TH1F("gPrimaryProtonsChi2PerClusterITSReject",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsChi2PerClusterITSReject);
  TH1F *gPrimaryProtonsTPCClustersReject = new TH1F("gPrimaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsTPCClustersReject);
  TH1F *gPrimaryProtonsChi2PerClusterTPCReject = new TH1F("gPrimaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsChi2PerClusterTPCReject);
  TH1F *gPrimaryProtonsExtCov11Reject = new TH1F("gPrimaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsExtCov11Reject);
  TH1F *gPrimaryProtonsExtCov22Reject = new TH1F("gPrimaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsExtCov22Reject);
  TH1F *gPrimaryProtonsExtCov33Reject = new TH1F("gPrimaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsExtCov33Reject);
  TH1F *gPrimaryProtonsExtCov44Reject = new TH1F("gPrimaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsExtCov44Reject);
  TH1F *gPrimaryProtonsExtCov55Reject = new TH1F("gPrimaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsExtCov55Reject);
  TH1F *gPrimaryProtonsSigmaToVertexReject = new TH1F("gPrimaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsSigmaToVertexReject);
  TH1F *gPrimaryProtonsSigmaToVertexTPCReject = new TH1F("gPrimaryProtonsSigmaToVertexTPCReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsSigmaToVertexTPCReject);
  TH1F *gPrimaryProtonsDCAXYReject = new TH1F("gPrimaryProtonsDCAXYReject",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsDCAXYReject);
  TH1F *gPrimaryProtonsDCAXYTPCReject = new TH1F("gPrimaryProtonsDCAXYTPCReject",
						 ";DCA_{xy} [cm];Entries",
						 100,0,20);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsDCAXYTPCReject);
  TH1F *gPrimaryProtonsDCAZReject = new TH1F("gPrimaryProtonsDCAZReject",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsDCAZReject);
  TH1F *gPrimaryProtonsDCAZTPCReject = new TH1F("gPrimaryProtonsDCAZTPCReject",
						";DCA_{z} [cm];Entries",
						100,0,20);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsDCAZTPCReject);
  TH1F *gPrimaryProtonsConstrainChi2Reject = new TH1F("gPrimaryProtonsConstrainChi2Reject",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsConstrainChi2Reject);
  TH1F *gPrimaryProtonsITSRefitReject = new TH1F("gPrimaryProtonsITSRefitReject",
						 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsITSRefitReject);
  TH1F *gPrimaryProtonsTPCRefitReject = new TH1F("gPrimaryProtonsTPCRefitReject",
						 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsTPCRefitReject);
  TH1F *gPrimaryProtonsESDpidReject = new TH1F("gPrimaryProtonsESDpidReject",
					       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsESDpidReject);
  TH1F *gPrimaryProtonsTPCpidReject = new TH1F("gPrimaryProtonsTPCpidReject",
					       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsTPCpidReject);
  TH1F *gPrimaryProtonsPointOnITSLayer1Reject = new TH1F("gPrimaryProtonsPointOnITSLayer1Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer1Reject);
  TH1F *gPrimaryProtonsPointOnITSLayer2Reject = new TH1F("gPrimaryProtonsPointOnITSLayer2Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer2Reject);
  TH1F *gPrimaryProtonsPointOnITSLayer3Reject = new TH1F("gPrimaryProtonsPointOnITSLayer3Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer3Reject);
  TH1F *gPrimaryProtonsPointOnITSLayer4Reject = new TH1F("gPrimaryProtonsPointOnITSLayer4Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer4Reject);
  TH1F *gPrimaryProtonsPointOnITSLayer5Reject = new TH1F("gPrimaryProtonsPointOnITSLayer5Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer5Reject);
  TH1F *gPrimaryProtonsPointOnITSLayer6Reject = new TH1F("gPrimaryProtonsPointOnITSLayer6Reject",
					     "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsPointOnITSLayer6Reject);

  //________________________________________________________________//
  /*gDirectory->cd("../../");

  TDirectory *dirProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirProtonsSecondary->cd();
  TDirectory *dirProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirProtonsSecondaryAccepted->cd();*/

  //Accepted secondary protons
  TH1F *gSecondaryProtonsITSClustersPass = new TH1F("gSecondaryProtonsITSClustersPass",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsITSClustersPass);
  TH1F *gSecondaryProtonsChi2PerClusterITSPass = new TH1F("gSecondaryProtonsChi2PerClusterITSPass",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsChi2PerClusterITSPass);
  TH1F *gSecondaryProtonsTPCClustersPass = new TH1F("gSecondaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsTPCClustersPass);
  TH1F *gSecondaryProtonsChi2PerClusterTPCPass = new TH1F("gSecondaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsChi2PerClusterTPCPass);
  TH1F *gSecondaryProtonsExtCov11Pass = new TH1F("gSecondaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsExtCov11Pass);
  TH1F *gSecondaryProtonsExtCov22Pass = new TH1F("gSecondaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsExtCov22Pass);
  TH1F *gSecondaryProtonsExtCov33Pass = new TH1F("gSecondaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsExtCov33Pass);
  TH1F *gSecondaryProtonsExtCov44Pass = new TH1F("gSecondaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsExtCov44Pass);
  TH1F *gSecondaryProtonsExtCov55Pass = new TH1F("gSecondaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsExtCov55Pass);
  TH1F *gSecondaryProtonsSigmaToVertexPass = new TH1F("gSecondaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsSigmaToVertexPass);
  TH1F *gSecondaryProtonsSigmaToVertexTPCPass = new TH1F("gSecondaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsSigmaToVertexTPCPass);
  TH1F *gSecondaryProtonsDCAXYPass = new TH1F("gSecondaryProtonsDCAXYPass",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsDCAXYPass);
  TH1F *gSecondaryProtonsDCAXYTPCPass = new TH1F("gSecondaryProtonsDCAXYTPCPass",
						 ";DCA_{xy} [cm];Entries",
						 100,0,20);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsDCAXYTPCPass);
  TH1F *gSecondaryProtonsDCAZPass = new TH1F("gSecondaryProtonsDCAZPass",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsDCAZPass);
  TH1F *gSecondaryProtonsDCAZTPCPass = new TH1F("gSecondaryProtonsDCAZTPCPass",
						";DCA_{z} [cm];Entries",
						100,0,20);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsDCAZTPCPass);
  TH1F *gSecondaryProtonsConstrainChi2Pass = new TH1F("gSecondaryProtonsConstrainChi2Pass",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsConstrainChi2Pass);
  TH1F *gSecondaryProtonsITSRefitPass = new TH1F("gSecondaryProtonsITSRefitPass",
						 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsITSRefitPass);
  TH1F *gSecondaryProtonsTPCRefitPass = new TH1F("gSecondaryProtonsTPCRefitPass",
						 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsTPCRefitPass);
  TH1F *gSecondaryProtonsESDpidPass = new TH1F("gSecondaryProtonsESDpidPass",
					       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsESDpidPass);
  TH1F *gSecondaryProtonsTPCpidPass = new TH1F("gSecondaryProtonsTPCpidPass",
					       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsTPCpidPass);
  TH1F *gSecondaryProtonsPointOnITSLayer1Pass = new TH1F("gSecondaryProtonsPointOnITSLayer1Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer1Pass);
  TH1F *gSecondaryProtonsPointOnITSLayer2Pass = new TH1F("gSecondaryProtonsPointOnITSLayer2Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer2Pass);
  TH1F *gSecondaryProtonsPointOnITSLayer3Pass = new TH1F("gSecondaryProtonsPointOnITSLayer3Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer3Pass);
  TH1F *gSecondaryProtonsPointOnITSLayer4Pass = new TH1F("gSecondaryProtonsPointOnITSLayer4Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer4Pass);
  TH1F *gSecondaryProtonsPointOnITSLayer5Pass = new TH1F("gSecondaryProtonsPointOnITSLayer5Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer5Pass);
  TH1F *gSecondaryProtonsPointOnITSLayer6Pass = new TH1F("gSecondaryProtonsPointOnITSLayer6Pass",
							 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsPointOnITSLayer6Pass);

  //Rejected secondary protons
  /*gDirectory->cd("../");
  TDirectory *dirProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsSecondaryRejected->cd();*/

  TH1F *gSecondaryProtonsITSClustersReject = new TH1F("gSecondaryProtonsITSClustersReject",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsITSClustersReject);
  TH1F *gSecondaryProtonsChi2PerClusterITSReject = new TH1F("gSecondaryProtonsChi2PerClusterITSReject",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsChi2PerClusterITSReject);
  TH1F *gSecondaryProtonsTPCClustersReject = new TH1F("gSecondaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsTPCClustersReject);
  TH1F *gSecondaryProtonsChi2PerClusterTPCReject = new TH1F("gSecondaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsChi2PerClusterTPCReject);
  TH1F *gSecondaryProtonsExtCov11Reject = new TH1F("gSecondaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsExtCov11Reject);
  TH1F *gSecondaryProtonsExtCov22Reject = new TH1F("gSecondaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsExtCov22Reject);
  TH1F *gSecondaryProtonsExtCov33Reject = new TH1F("gSecondaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsExtCov33Reject);
  TH1F *gSecondaryProtonsExtCov44Reject = new TH1F("gSecondaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsExtCov44Reject);
  TH1F *gSecondaryProtonsExtCov55Reject = new TH1F("gSecondaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsExtCov55Reject);
  TH1F *gSecondaryProtonsSigmaToVertexReject = new TH1F("gSecondaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsSigmaToVertexReject);
  TH1F *gSecondaryProtonsSigmaToVertexTPCReject = new TH1F("gSecondaryProtonsSigmaToVertexTPCReject",
							   ";#sigma_{Vertex};Entries",
							   100,0,10);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsSigmaToVertexTPCReject);
  TH1F *gSecondaryProtonsDCAXYReject = new TH1F("gSecondaryProtonsDCAXYReject",
						";DCA_{xy} [cm];Entries",
						100,0,20);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsDCAXYReject);
  TH1F *gSecondaryProtonsDCAXYTPCReject = new TH1F("gSecondaryProtonsDCAXYTPCReject",
						   ";DCA_{xy} [cm];Entries",
						   100,0,20);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsDCAXYTPCReject);
  TH1F *gSecondaryProtonsDCAZReject = new TH1F("gSecondaryProtonsDCAZReject",
					       ";DCA_{z} [cm];Entries",
					       100,0,20);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsDCAZReject);
  TH1F *gSecondaryProtonsDCAZTPCReject = new TH1F("gSecondaryProtonsDCAZTPCReject",
						  ";DCA_{z} [cm];Entries",
						  100,0,20);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsDCAZTPCReject);
  TH1F *gSecondaryProtonsConstrainChi2Reject = new TH1F("gSecondaryProtonsConstrainChi2Reject",
							";Log_{10}(#chi^{2});Entries",
							100,-10,10);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsConstrainChi2Reject);
  TH1F *gSecondaryProtonsITSRefitReject = new TH1F("gSecondaryProtonsITSRefitReject",
						   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsITSRefitReject);
  TH1F *gSecondaryProtonsTPCRefitReject = new TH1F("gSecondaryProtonsTPCRefitReject",
						   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsTPCRefitReject);
  TH1F *gSecondaryProtonsESDpidReject = new TH1F("gSecondaryProtonsESDpidReject",
						 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsESDpidReject);
  TH1F *gSecondaryProtonsTPCpidReject = new TH1F("gSecondaryProtonsTPCpidReject",
						 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsTPCpidReject);
  TH1F *gSecondaryProtonsPointOnITSLayer1Reject = new TH1F("gSecondaryProtonsPointOnITSLayer1Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer1Reject);
  TH1F *gSecondaryProtonsPointOnITSLayer2Reject = new TH1F("gSecondaryProtonsPointOnITSLayer2Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer2Reject);
  TH1F *gSecondaryProtonsPointOnITSLayer3Reject = new TH1F("gSecondaryProtonsPointOnITSLayer3Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer3Reject);
  TH1F *gSecondaryProtonsPointOnITSLayer4Reject = new TH1F("gSecondaryProtonsPointOnITSLayer4Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer4Reject);
  TH1F *gSecondaryProtonsPointOnITSLayer5Reject = new TH1F("gSecondaryProtonsPointOnITSLayer5Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer5Reject);
  TH1F *gSecondaryProtonsPointOnITSLayer6Reject = new TH1F("gSecondaryProtonsPointOnITSLayer6Reject",
							   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsPointOnITSLayer6Reject);
  

  /*gDirectory->cd("../../../");

  //antiprotons
  TDirectory *dirAntiProtons = gDirectory->mkdir("AntiProtons");
  fGlobalQAList->Add(dirAntiProtons); dirAntiProtons->cd();*/
  
  //________________________________________________________________//
  /*TDirectory *dirAntiProtonsPrimary = gDirectory->mkdir("Primaries");
  dirAntiProtonsPrimary->cd();
  TDirectory *dirAntiProtonsPrimaryAccepted = gDirectory->mkdir("Accepted");
  dirAntiProtonsPrimaryAccepted->cd();*/
  
  //Accepted primary antiprotons
  TH1F *gPrimaryAntiProtonsITSClustersPass = new TH1F("gPrimaryAntiProtonsITSClustersPass",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsITSClustersPass);
  TH1F *gPrimaryAntiProtonsChi2PerClusterITSPass = new TH1F("gPrimaryAntiProtonsChi2PerClusterITSPass",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsChi2PerClusterITSPass);
  TH1F *gPrimaryAntiProtonsTPCClustersPass = new TH1F("gPrimaryAntiProtonsTPCClustersPass",
						      ";N_{clusters} (TPC);Entries",
						      100,0,200);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsTPCClustersPass);
  TH1F *gPrimaryAntiProtonsChi2PerClusterTPCPass = new TH1F("gPrimaryAntiProtonsChi2PerClusterTPCPass",
							    ";x^{2}/N_{clusters} (TPC);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *gPrimaryAntiProtonsExtCov11Pass = new TH1F("gPrimaryAntiProtonsExtCov11Pass",
						   ";#sigma_{y} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsExtCov11Pass);
  TH1F *gPrimaryAntiProtonsExtCov22Pass = new TH1F("gPrimaryAntiProtonsExtCov22Pass",
						   ";#sigma_{z} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsExtCov22Pass);
  TH1F *gPrimaryAntiProtonsExtCov33Pass = new TH1F("gPrimaryAntiProtonsExtCov33Pass",
						   ";#sigma_{sin(#phi)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsExtCov33Pass);
  TH1F *gPrimaryAntiProtonsExtCov44Pass = new TH1F("gPrimaryAntiProtonsExtCov44Pass",
						   ";#sigma_{tan(#lambda)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsExtCov44Pass);
  TH1F *gPrimaryAntiProtonsExtCov55Pass = new TH1F("gPrimaryAntiProtonsExtCov55Pass",
						   ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsExtCov55Pass);
  TH1F *gPrimaryAntiProtonsSigmaToVertexPass = new TH1F("gPrimaryAntiProtonsSigmaToVertexPass",
							";#sigma_{Vertex};Entries",
							100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsSigmaToVertexPass);
  TH1F *gPrimaryAntiProtonsSigmaToVertexTPCPass = new TH1F("gPrimaryAntiProtonsSigmaToVertexTPCPass",
							   ";#sigma_{Vertex};Entries",
							   100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *gPrimaryAntiProtonsDCAXYPass = new TH1F("gPrimaryAntiProtonsDCAXYPass",
						";DCA_{xy} [cm];Entries",
						100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsDCAXYPass);
  TH1F *gPrimaryAntiProtonsDCAXYTPCPass = new TH1F("gPrimaryAntiProtonsDCAXYTPCPass",
						   ";DCA_{xy} [cm];Entries",
						   100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsDCAXYTPCPass);
  TH1F *gPrimaryAntiProtonsDCAZPass = new TH1F("gPrimaryAntiProtonsDCAZPass",
					       ";DCA_{z} [cm];Entries",
					       100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsDCAZPass);
  TH1F *gPrimaryAntiProtonsDCAZTPCPass = new TH1F("gPrimaryAntiProtonsDCAZTPCPass",
						  ";DCA_{z} [cm];Entries",
						  100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsDCAZTPCPass);
  TH1F *gPrimaryAntiProtonsConstrainChi2Pass = new TH1F("gPrimaryAntiProtonsConstrainChi2Pass",
							";Log_{10}(#chi^{2});Entries",
							100,-10,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsConstrainChi2Pass);
  TH1F *gPrimaryAntiProtonsITSRefitPass = new TH1F("gPrimaryAntiProtonsITSRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsITSRefitPass);
  TH1F *gPrimaryAntiProtonsTPCRefitPass = new TH1F("gPrimaryAntiProtonsTPCRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsTPCRefitPass);
  TH1F *gPrimaryAntiProtonsESDpidPass = new TH1F("gPrimaryAntiProtonsESDpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsESDpidPass);
  TH1F *gPrimaryAntiProtonsTPCpidPass = new TH1F("gPrimaryAntiProtonsTPCpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsTPCpidPass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer1Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer1Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer1Pass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer2Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer2Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer2Pass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer3Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer3Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer3Pass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer4Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer4Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer4Pass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer5Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer5Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer5Pass);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer6Pass = new TH1F("gPrimaryAntiProtonsPointOnITSLayer6Pass",
							   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsPointOnITSLayer6Pass);
  
  //Rejected primary antiprotons
  /*gDirectory->cd("../");
  TDirectory *dirAntiProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsPrimaryRejected->cd();*/
  
  TH1F *gPrimaryAntiProtonsITSClustersReject = new TH1F("gPrimaryAntiProtonsITSClustersReject",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsITSClustersReject);
  TH1F *gPrimaryAntiProtonsChi2PerClusterITSReject = new TH1F("gPrimaryAntiProtonsChi2PerClusterITSReject",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsChi2PerClusterITSReject);
  TH1F *gPrimaryAntiProtonsTPCClustersReject = new TH1F("gPrimaryAntiProtonsTPCClustersReject",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsTPCClustersReject);
  TH1F *gPrimaryAntiProtonsChi2PerClusterTPCReject = new TH1F("gPrimaryAntiProtonsChi2PerClusterTPCReject",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *gPrimaryAntiProtonsExtCov11Reject = new TH1F("gPrimaryAntiProtonsExtCov11Reject",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsExtCov11Reject);
  TH1F *gPrimaryAntiProtonsExtCov22Reject = new TH1F("gPrimaryAntiProtonsExtCov22Reject",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsExtCov22Reject);
  TH1F *gPrimaryAntiProtonsExtCov33Reject = new TH1F("gPrimaryAntiProtonsExtCov33Reject",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsExtCov33Reject);
  TH1F *gPrimaryAntiProtonsExtCov44Reject = new TH1F("gPrimaryAntiProtonsExtCov44Reject",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsExtCov44Reject);
  TH1F *gPrimaryAntiProtonsExtCov55Reject = new TH1F("gPrimaryAntiProtonsExtCov55Reject",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsExtCov55Reject);
  TH1F *gPrimaryAntiProtonsSigmaToVertexReject = new TH1F("gPrimaryAntiProtonsSigmaToVertexReject",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsSigmaToVertexReject);
  TH1F *gPrimaryAntiProtonsSigmaToVertexTPCReject = new TH1F("gPrimaryAntiProtonsSigmaToVertexTPCReject",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *gPrimaryAntiProtonsDCAXYReject = new TH1F("gPrimaryAntiProtonsDCAXYReject",
						  ";DCA_{xy} [cm];Entries",
						  100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsDCAXYReject);
  TH1F *gPrimaryAntiProtonsDCAXYTPCReject = new TH1F("gPrimaryAntiProtonsDCAXYTPCReject",
						     ";DCA_{xy} [cm];Entries",
						     100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsDCAXYTPCReject);
  TH1F *gPrimaryAntiProtonsDCAZReject = new TH1F("gPrimaryAntiProtonsDCAZReject",
						 ";DCA_{z} [cm];Entries",
						 100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsDCAZReject);
  TH1F *gPrimaryAntiProtonsDCAZTPCReject = new TH1F("gPrimaryAntiProtonsDCAZTPCReject",
						    ";DCA_{z} [cm];Entries",
						    100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsDCAZTPCReject);
  TH1F *gPrimaryAntiProtonsConstrainChi2Reject = new TH1F("gPrimaryAntiProtonsConstrainChi2Reject",
							  ";Log_{10}(#chi^{2});Entries",
							  100,-10,10);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsConstrainChi2Reject);
  TH1F *gPrimaryAntiProtonsITSRefitReject = new TH1F("gPrimaryAntiProtonsITSRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsITSRefitReject);
  TH1F *gPrimaryAntiProtonsTPCRefitReject = new TH1F("gPrimaryAntiProtonsTPCRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsTPCRefitReject);
  TH1F *gPrimaryAntiProtonsESDpidReject = new TH1F("gPrimaryAntiProtonsESDpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsESDpidReject);
  TH1F *gPrimaryAntiProtonsTPCpidReject = new TH1F("gPrimaryAntiProtonsTPCpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsTPCpidReject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer1Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer1Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer1Reject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer2Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer2Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer2Reject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer3Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer3Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer3Reject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer4Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer4Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer4Reject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer5Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer5Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer5Reject);
  TH1F *gPrimaryAntiProtonsPointOnITSLayer6Reject = new TH1F("gPrimaryAntiProtonsPointOnITSLayer6Reject",
							     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsPointOnITSLayer6Reject);
  
  //________________________________________________________________//
  /*gDirectory->cd("../../");

  TDirectory *dirAntiProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirAntiProtonsSecondary->cd();
  TDirectory *dirAntiProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirAntiProtonsSecondaryAccepted->cd();*/

  //Accepted secondary antiprotons
  TH1F *gSecondaryAntiProtonsITSClustersPass = new TH1F("gSecondaryAntiProtonsITSClustersPass",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsITSClustersPass);
  TH1F *gSecondaryAntiProtonsChi2PerClusterITSPass = new TH1F("gSecondaryAntiProtonsChi2PerClusterITSPass",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsChi2PerClusterITSPass);
  TH1F *gSecondaryAntiProtonsTPCClustersPass = new TH1F("gSecondaryAntiProtonsTPCClustersPass",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsTPCClustersPass);
  TH1F *gSecondaryAntiProtonsChi2PerClusterTPCPass = new TH1F("gSecondaryAntiProtonsChi2PerClusterTPCPass",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *gSecondaryAntiProtonsExtCov11Pass = new TH1F("gSecondaryAntiProtonsExtCov11Pass",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsExtCov11Pass);
  TH1F *gSecondaryAntiProtonsExtCov22Pass = new TH1F("gSecondaryAntiProtonsExtCov22Pass",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsExtCov22Pass);
  TH1F *gSecondaryAntiProtonsExtCov33Pass = new TH1F("gSecondaryAntiProtonsExtCov33Pass",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsExtCov33Pass);
  TH1F *gSecondaryAntiProtonsExtCov44Pass = new TH1F("gSecondaryAntiProtonsExtCov44Pass",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsExtCov44Pass);
  TH1F *gSecondaryAntiProtonsExtCov55Pass = new TH1F("gSecondaryAntiProtonsExtCov55Pass",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsExtCov55Pass);
  TH1F *gSecondaryAntiProtonsSigmaToVertexPass = new TH1F("gSecondaryAntiProtonsSigmaToVertexPass",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsSigmaToVertexPass);
  TH1F *gSecondaryAntiProtonsSigmaToVertexTPCPass = new TH1F("gSecondaryAntiProtonsSigmaToVertexTPCPass",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *gSecondaryAntiProtonsDCAXYPass = new TH1F("gSecondaryAntiProtonsDCAXYPass",
						  ";DCA_{xy} [cm];Entries",
						  100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsDCAXYPass);
  TH1F *gSecondaryAntiProtonsDCAXYTPCPass = new TH1F("gSecondaryAntiProtonsDCAXYTPCPass",
						     ";DCA_{xy} [cm];Entries",
						     100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsDCAXYTPCPass);
  TH1F *gSecondaryAntiProtonsDCAZPass = new TH1F("gSecondaryAntiProtonsDCAZPass",
						 ";DCA_{z} [cm];Entries",
						 100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsDCAZPass);
  TH1F *gSecondaryAntiProtonsDCAZTPCPass = new TH1F("gSecondaryAntiProtonsDCAZTPCPass",
						    ";DCA_{z} [cm];Entries",
						    100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsDCAZTPCPass);
  TH1F *gSecondaryAntiProtonsConstrainChi2Pass = new TH1F("gSecondaryAntiProtonsConstrainChi2Pass",
							  ";Log_{10}(#chi^{2});Entries",
							  100,-10,10);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsConstrainChi2Pass);
  TH1F *gSecondaryAntiProtonsITSRefitPass = new TH1F("gSecondaryAntiProtonsITSRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsITSRefitPass);
  TH1F *gSecondaryAntiProtonsTPCRefitPass = new TH1F("gSecondaryAntiProtonsTPCRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsTPCRefitPass);
  TH1F *gSecondaryAntiProtonsESDpidPass = new TH1F("gSecondaryAntiProtonsESDpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsESDpidPass);
  TH1F *gSecondaryAntiProtonsTPCpidPass = new TH1F("gSecondaryAntiProtonsTPCpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsTPCpidPass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer1Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer1Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer1Pass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer2Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer2Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer2Pass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer3Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer3Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer3Pass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer4Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer4Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer4Pass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer5Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer5Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer5Pass);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer6Pass = new TH1F("gSecondaryAntiProtonsPointOnITSLayer6Pass",
							     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsPointOnITSLayer6Pass);
  
  //Rejected secondary antiprotons
  /*gDirectory->cd("../");
  TDirectory *dirAntiProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsSecondaryRejected->cd();*/

  TH1F *gSecondaryAntiProtonsITSClustersReject = new TH1F("gSecondaryAntiProtonsITSClustersReject",
							  ";N_{clusters} (ITS);Entries",
							  7,0,7);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsITSClustersReject);
  TH1F *gSecondaryAntiProtonsChi2PerClusterITSReject = new TH1F("gSecondaryAntiProtonsChi2PerClusterITSReject",
								";x^{2}/N_{clusters} (ITS);Entries",
								100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsChi2PerClusterITSReject);
  TH1F *gSecondaryAntiProtonsTPCClustersReject = new TH1F("gSecondaryAntiProtonsTPCClustersReject",
							  ";N_{clusters} (TPC);Entries",
							  100,0,200);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsTPCClustersReject);
  TH1F *gSecondaryAntiProtonsChi2PerClusterTPCReject = new TH1F("gSecondaryAntiProtonsChi2PerClusterTPCReject",
								";x^{2}/N_{clusters} (TPC);Entries",
								100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *gSecondaryAntiProtonsExtCov11Reject = new TH1F("gSecondaryAntiProtonsExtCov11Reject",
						       ";#sigma_{y} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsExtCov11Reject);
  TH1F *gSecondaryAntiProtonsExtCov22Reject = new TH1F("gSecondaryAntiProtonsExtCov22Reject",
						       ";#sigma_{z} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsExtCov22Reject);
  TH1F *gSecondaryAntiProtonsExtCov33Reject = new TH1F("gSecondaryAntiProtonsExtCov33Reject",
						       ";#sigma_{sin(#phi)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsExtCov33Reject);
  TH1F *gSecondaryAntiProtonsExtCov44Reject = new TH1F("gSecondaryAntiProtonsExtCov44Reject",
						       ";#sigma_{tan(#lambda)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsExtCov44Reject);
  TH1F *gSecondaryAntiProtonsExtCov55Reject = new TH1F("gSecondaryAntiProtonsExtCov55Reject",
						       ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsExtCov55Reject);
  TH1F *gSecondaryAntiProtonsSigmaToVertexReject = new TH1F("gSecondaryAntiProtonsSigmaToVertexReject",
							    ";#sigma_{Vertex};Entries",
							    100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsSigmaToVertexReject);
  TH1F *gSecondaryAntiProtonsSigmaToVertexTPCReject = new TH1F("gSecondaryAntiProtonsSigmaToVertexTPCReject",
							       ";#sigma_{Vertex};Entries",
							       100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *gSecondaryAntiProtonsDCAXYReject = new TH1F("gSecondaryAntiProtonsDCAXYReject",
						    ";DCA_{xy} [cm];Entries",
						    100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsDCAXYReject);
  TH1F *gSecondaryAntiProtonsDCAXYTPCReject = new TH1F("gSecondaryAntiProtonsDCAXYTPCReject",
						       ";DCA_{xy} [cm];Entries",
						       100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsDCAXYTPCReject);
  TH1F *gSecondaryAntiProtonsDCAZReject = new TH1F("gSecondaryAntiProtonsDCAZReject",
						   ";DCA_{z} [cm];Entries",
						   100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsDCAZReject);
  TH1F *gSecondaryAntiProtonsDCAZTPCReject = new TH1F("gSecondaryAntiProtonsDCAZTPCReject",
						      ";DCA_{z} [cm];Entries",
						      100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsDCAZTPCReject);
  TH1F *gSecondaryAntiProtonsConstrainChi2Reject = new TH1F("gSecondaryAntiProtonsConstrainChi2Reject",
							    ";Log_{10}(#chi^{2});Entries",
							    100,-10,10);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsConstrainChi2Reject);
  TH1F *gSecondaryAntiProtonsITSRefitReject = new TH1F("gSecondaryAntiProtonsITSRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsITSRefitReject);
  TH1F *gSecondaryAntiProtonsTPCRefitReject = new TH1F("gSecondaryAntiProtonsTPCRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsTPCRefitReject);
  TH1F *gSecondaryAntiProtonsESDpidReject = new TH1F("gSecondaryAntiProtonsESDpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsESDpidReject);
  TH1F *gSecondaryAntiProtonsTPCpidReject = new TH1F("gSecondaryAntiProtonsTPCpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsTPCpidReject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer1Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer1Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer1Reject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer2Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer2Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer2Reject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer3Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer3Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer3Reject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer4Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer4Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer4Reject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer5Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer5Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer5Reject);
  TH1F *gSecondaryAntiProtonsPointOnITSLayer6Reject = new TH1F("gSecondaryAntiProtonsPointOnITSLayer6Reject",
							     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsPointOnITSLayer6Reject);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunEfficiencyAnalysis(AliStack *stack, 
						AliESDEvent *esd) {
  //Runs the efficiency code
  //MC loop
  Int_t nMCProtons = 0, nESDProtons = 0;
  for(Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if(fRunEfficiencyAnalysisEtaMode) {
      if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
    }
    else 
      if((Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;

    Int_t pdgcode = particle->GetPdgCode();
    if(TMath::Abs(pdgcode) != 2212) continue;

    if(iParticle <= stack->GetNprimary()) {
      if(pdgcode == 2212) {
	nMCProtons += 1;
	if(fRunEfficiencyAnalysisEtaMode) 
	  ((TH2D *)(fEfficiencyList->At(0)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(0)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
      }//protons
      if(pdgcode == -2212) {
	if(fRunEfficiencyAnalysisEtaMode) 
	  ((TH2D *)(fEfficiencyList->At(1)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(1)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
      }//antiprotons
    }//primaries
    else {
      //secondaries
      Int_t lPartMother = -1;
      Int_t motherPDGCode = -1;
      lPartMother = particle->GetFirstMother();
      TParticle *motherParticle = stack->Particle(lPartMother);
      if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();

      if(pdgcode == 2212) {
	if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	  if(fRunEfficiencyAnalysisEtaMode) 
	    ((TH2D *)(fEfficiencyList->At(2)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(2)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//weak decays
	if((particle->GetUniqueID() == 13)) {
	  if(fRunEfficiencyAnalysisEtaMode) 
	    ((TH2D *)(fEfficiencyList->At(4)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(4)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//hadronic interactions
      }//protons
      if(pdgcode == -2212) {
	if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	  if(fRunEfficiencyAnalysisEtaMode) 
	    ((TH2D *)(fEfficiencyList->At(3)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(3)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//weak decays
	if((particle->GetUniqueID() == 13)) {
	  if(fRunEfficiencyAnalysisEtaMode) 
	    ((TH2D *)(fEfficiencyList->At(5)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(5)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//hadronic interactions
      }//antiprotons
    }//secondaries
  
  }//MC loop

  //ESD loop
  Int_t nGoodTracks = esd->GetNumberOfTracks();
  TArrayI labelArray(nGoodTracks);
  Int_t labelCounter = 0;
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    if(!track) continue;

    Int_t label = TMath::Abs(track->GetLabel());
    if(IsLabelUsed(labelArray,label)) continue;
    labelArray.AddAt(label,labelCounter);
    labelCounter += 1;

    TParticle *particle = stack->Particle(label);
    if(!particle) continue;
    Int_t pdgcode = particle->GetPdgCode();
    if(TMath::Abs(pdgcode) != 2212) continue;
    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    
    Double_t Pt = 0.0, P = 0.0;
    Double_t probability[5];
    
    //TPC only
    if(fUseTPCOnly) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      Pt = tpcTrack->Pt();
      P = tpcTrack->P();

      if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
      if(fRunEfficiencyAnalysisEtaMode) {
	if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
      }
      else 
	if((Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
      
      if(fUseCutsInEfficiency) 
	if(!IsAccepted(track)) continue;

      //reconstructed primary (anti)protons
      if(pdgcode == 2212) {
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
	if(label <= stack->GetNprimary()) {
	  nESDProtons += 1;
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//primaries
	if(label > stack->GetNprimary()) {
	  Int_t lPartMother = -1;
	  Int_t motherPDGCode = -1;
	  lPartMother = particle->GetFirstMother();
	  TParticle *motherParticle = stack->Particle(lPartMother);
	  if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();

	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(Rapidity(particle->Px(),
								 particle->Py(),
								 particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial protons
      if(pdgcode == -2212) {	
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
	if(label <= stack->GetNprimary()) {
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//primaries
	if(label > stack->GetNprimary()) {
	  Int_t lPartMother = -1;
	  Int_t motherPDGCode = -1;
	  lPartMother = particle->GetFirstMother();
	  TParticle *motherParticle = stack->Particle(lPartMother);
	  if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();

	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(Rapidity(particle->Px(),
								 particle->Py(),
								 particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial antiprotons
      
      //pid
      track->GetTPCpid(probability);
      Double_t rcc = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	rcc += probability[i]*GetParticleFraction(i,P);
      if(rcc == 0.0) continue;
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
      Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
      if(fParticleType == 4) {
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(14)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(14)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
	if(TMath::Abs(pdgcode) == 2212) 
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(13)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(13)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	else
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(15)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(15)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
      }//identified as proton
    }//TPC only tracks
    else if(!fUseTPCOnly) {
      if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
      if(fRunEfficiencyAnalysisEtaMode) {
	if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
      }
      else {
	if((Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
      }
      
      if(fUseCutsInEfficiency) 
	if(!IsAccepted(track)) continue;

      //reconstructed primary (anti)protons
      if(pdgcode == 2212) {
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
	if(label <= stack->GetNprimary()) {
	  nESDProtons += 1;
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//primaries
	if(label > stack->GetNprimary()) {
	  Int_t lPartMother = -1;
	  Int_t motherPDGCode = -1;
	  lPartMother = particle->GetFirstMother();
	  TParticle *motherParticle = stack->Particle(lPartMother);
	  if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();

	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(Rapidity(particle->Px(),
								 particle->Py(),
								 particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial protons
      if(pdgcode == -2212) {	
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
	if(label <= stack->GetNprimary()) {
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//primaries
	if(label > stack->GetNprimary()) {
	  Int_t lPartMother = -1;
	  Int_t motherPDGCode = -1;
	  lPartMother = particle->GetFirstMother();
	  TParticle *motherParticle = stack->Particle(lPartMother);
	  if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();

	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fRunEfficiencyAnalysisEtaMode)
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(Rapidity(particle->Px(),
								 particle->Py(),
								 particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial antiprotons
      
      //pid
      track->GetESDpid(probability);
      Double_t rcc = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	rcc += probability[i]*GetParticleFraction(i,P);
      if(rcc == 0.0) continue;
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
      Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
      if(fParticleType == 4) {
	if(fRunEfficiencyAnalysisEtaMode)
	  ((TH2D *)(fEfficiencyList->At(14)))->Fill(particle->Eta(),
						   particle->Pt());
	else ((TH2D *)(fEfficiencyList->At(14)))->Fill(Rapidity(particle->Px(),
							       particle->Py(),
							       particle->Pz()),
						      particle->Pt());
	if(TMath::Abs(pdgcode) == 2212) 
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(13)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(13)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	else
	  if(fRunEfficiencyAnalysisEtaMode)
	    ((TH2D *)(fEfficiencyList->At(15)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(15)))->Fill(Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
      }//identified as proton
    }//global tracking
  }//track loop
  
  //Printf("MC protons: %d - ESD protons: %d",nMCProtons,nESDProtons);
}

//____________________________________________________________________//
Bool_t AliProtonQAAnalysis::IsLabelUsed(TArrayI labelArray, 
					Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunQAAnalysis(AliStack *stack, 
					AliESDEvent *esd) {
  //Runs the QA code
  //MC loop
  for(Int_t iParticle = 0; iParticle < stack->GetNprimary(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if((Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;

    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) 
      ((TH2D *)(fQA2DList->At(8)))->Fill(Rapidity(particle->Px(),
						  particle->Py(),
						  particle->Pz()),
					 particle->Pt());
    if(pdgcode == -2212) 
      ((TH2D *)(fQA2DList->At(9)))->Fill(Rapidity(particle->Px(),
						  particle->Py(),
						  particle->Pz()),
					 particle->Pt());
  }//MC loop

  //ESD loop
  Int_t nGoodTracks = esd->GetNumberOfTracks();
  TArrayI labelArray(nGoodTracks);
  Int_t labelCounter = 0;
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    if(!track) continue;
    
    Int_t label = TMath::Abs(track->GetLabel()); 
    if(IsLabelUsed(labelArray,label)) continue;
    labelArray.AddAt(label,labelCounter);
    labelCounter += 1;

    AliESDtrack trackTPC;
    
    //in case it's a TPC only track relate it to the proper vertex
    if((fUseTPCOnly)&&(!fUseHybridTPC)) {
      Float_t p[2],cov[3];
      track->GetImpactParametersTPC(p,cov);
      if (p[0]==0 && p[1]==0)  
	track->RelateToVertexTPC(((AliESDEvent*)esd)->GetPrimaryVertexTPC(),esd->GetMagneticField(),kVeryBig);
      if (!track->FillTPCOnlyTrack(trackTPC)) {
	continue;
      }
      track = &trackTPC ;
    }

    Double_t Pt = 0.0, P = 0.0;
    Double_t probability[5];
    Float_t dcaXY = 0.0, dcaZ = 0.0;
    Double_t nSigmaToVertex = GetSigmaToVertex(track);
    Int_t  fIdxInt[200];
    Int_t nClustersITS = track->GetITSclusters(fIdxInt);
    Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);

    Float_t chi2PerClusterITS = -1;
    if (nClustersITS!=0)
      chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
    Float_t chi2PerClusterTPC = -1;
    if (nClustersTPC!=0)
      chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);
    Double_t chi2ConstrainVertex = TMath::Log(track->GetConstrainedChi2());    
    Double_t extCov[15];
    track->GetExternalCovariance(extCov);

    //TPC only
    if(fUseTPCOnly) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      Pt = tpcTrack->Pt();
      P = tpcTrack->P();
      if(fUseHybridTPC) track->GetImpactParameters(dcaXY,dcaZ);
      else track->GetImpactParametersTPC(dcaXY,dcaZ);
      
      //pid
      track->GetTPCpid(probability);
      Double_t rcc = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	rcc += probability[i]*GetParticleFraction(i,P);
      if(rcc == 0.0) continue;
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
      Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
      if(fParticleType == 4) {
	FillQA(track, stack);
	if(IsAccepted(track)) {
	  if(label <= stack->GetNprimary()) {
            if(track->Charge() > 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(0)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(4)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(8)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(12)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(16)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(20)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(24)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(28)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(32)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(36)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(40)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(0)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(4)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(8)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(0)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	    }
            else if(track->Charge() < 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(1)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(5)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(9)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(13)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(17)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(21)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(25)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(29)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(33)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(37)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(41)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(1)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(5)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(9)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(4)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	    }
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    TParticle *particle = stack->Particle(label);
	    if(!particle) continue;

	    Int_t lPartMother = -1;
	    Int_t motherPDGCode = -1;
	    if(particle) {
	      lPartMother = particle->GetFirstMother();
	      TParticle *motherParticle = stack->Particle(lPartMother);
	      if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();
	    }

	    if(fMCProcessIdFlag)
	      if(particle->GetUniqueID() != fMCProcessId) continue;
	    if(fMotherParticlePDGCodeFlag)
	      if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	    if(track->Charge() > 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(2)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(6)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(10)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(14)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(18)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(22)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(26)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(30)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(34)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(38)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(42)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(2)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(6)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(10)))->Fill(nSigmaToVertex);
	      ((TH2D *)(fQA2DList->At(2)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	      ((TH3F *)(fQA2DList->At(10)))->Fill(Rapidity(tpcTrack->Px(),
							   tpcTrack->Py(),
							   tpcTrack->Pz()),
						  Pt,
						  ConvertPDGToInt(motherPDGCode));
	    }
            else if(track->Charge() < 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(3)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(7)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(11)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(15)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(19)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(23)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(27)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(31)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(35)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(39)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(43)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(3)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(7)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(11)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(6)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	      ((TH3F *)(fQA2DList->At(11)))->Fill(Rapidity(tpcTrack->Px(),
							   tpcTrack->Py(),
							   tpcTrack->Pz()),
						  Pt,
						  ConvertPDGToInt(motherPDGCode));
	    }
	  }//secondary particles
	}//accepted - track cuts
	else if(!IsAccepted(track)) {
	  if(label <= stack->GetNprimary()) {
            if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(1)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(5)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(3)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(7)))->Fill(Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz()),
						 Pt);
	  }//secondary particles
	}//rejected - track cuts
      }//proton check
    }//TPC only tracks
    //combined tracking
    else if(!fUseTPCOnly) {
      Pt = track->Pt();
      P = track->P();
      track->GetImpactParameters(dcaXY,dcaZ);

      //pid
      track->GetESDpid(probability);
      Double_t rcc = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	rcc += probability[i]*GetParticleFraction(i,P);
      if(rcc == 0.0) continue;
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++)
	w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
      Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
      if(fParticleType == 4) {
	FillQA(track, stack);
	if(IsAccepted(track)) {
	  if(label <= stack->GetNprimary()) {
            if(track->Charge() > 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(0)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(4)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(8)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(12)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(16)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(20)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(24)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(28)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(32)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(36)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(40)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(0)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(4)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(8)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(0)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	    }
            else if(track->Charge() < 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(1)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(5)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(9)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(13)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(17)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(21)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(25)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(29)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(33)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(37)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(41)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(1)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(5)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(9)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(4)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	    }
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    TParticle *particle = stack->Particle(label);
	    if(!particle) continue;

	    Int_t lPartMother = -1;
	    Int_t motherPDGCode = -1;
	    if(particle) {
	      lPartMother = particle->GetFirstMother();
	      TParticle *motherParticle = stack->Particle(lPartMother);
	      if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();
	    }

	    if(fMCProcessIdFlag)
	      if(particle->GetUniqueID() != fMCProcessId) continue;
	    if(fMotherParticlePDGCodeFlag)
	      if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	    if(track->Charge() > 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(2)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(6)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(10)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(14)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(18)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(22)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(26)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(30)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(34)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(38)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(42)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(2)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(6)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(10)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(2)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	      ((TH3F *)(fQA2DList->At(10)))->Fill(Rapidity(track->Px(),
							   track->Py(),
							   track->Pz()),
						  Pt,
						  ConvertPDGToInt(motherPDGCode));
	    }
            else if(track->Charge() < 0) {
	      for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
		if(track->HasPointOnITSLayer(iLayer))
		  ((TH1F *)(fAcceptedCutList->At(3)))->Fill(iLayer+1);
	      }
	      ((TH1F *)(fAcceptedCutList->At(7)))->Fill(nClustersITS);
	      ((TH1F *)(fAcceptedCutList->At(11)))->Fill(chi2PerClusterITS);
	      ((TH1F *)(fAcceptedCutList->At(15)))->Fill(chi2ConstrainVertex);
	      ((TH1F *)(fAcceptedCutList->At(19)))->Fill(nClustersTPC);
	      ((TH1F *)(fAcceptedCutList->At(23)))->Fill(chi2PerClusterTPC);
	      ((TH1F *)(fAcceptedCutList->At(27)))->Fill(extCov[0]);
	      ((TH1F *)(fAcceptedCutList->At(31)))->Fill(extCov[2]);
	      ((TH1F *)(fAcceptedCutList->At(35)))->Fill(extCov[5]);
	      ((TH1F *)(fAcceptedCutList->At(39)))->Fill(extCov[9]);
	      ((TH1F *)(fAcceptedCutList->At(43)))->Fill(extCov[14]);

	      ((TH1F *)(fAcceptedDCAList->At(3)))->Fill(TMath::Abs(dcaXY));
	      ((TH1F *)(fAcceptedDCAList->At(7)))->Fill(TMath::Abs(dcaZ));
	      ((TH1F *)(fAcceptedDCAList->At(11)))->Fill(nSigmaToVertex);
              ((TH2D *)(fQA2DList->At(6)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);

	      ((TH3F *)(fQA2DList->At(11)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt,
						 ConvertPDGToInt(motherPDGCode));
	    }
	  }//secondary particles
	}//accepted - track cuts
	else if(!IsAccepted(track)) {
	  if(label <= stack->GetNprimary()) {
            if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(1)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(5)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(3)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(7)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//secondary particles
	}//rejected - track cuts
      }//proton check
    }//combined tracking
  }//track loop
    
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitMCAnalysis() {
  //MC analysis - 3D histograms: y-pT-pdg
  fPDGList = new TList();
  TH3F *gHistYPtPDGProtons = new TH3F("gHistYPtPDGProtons",
				      ";y;P_{T} [GeV/c];PDG",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt,
				      14,-0.5,13.5);
  fPDGList->Add(gHistYPtPDGProtons);
  TH3F *gHistYPtPDGAntiProtons = new TH3F("gHistYPtPDGAntiProtons",
					  ";y;P_{T} [GeV/c];PDG",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt,
					  14,-0.5,13.5);
  fPDGList->Add(gHistYPtPDGAntiProtons);

  //MC processes
  fMCProcessesList = new TList();
  TH1F *gHistProtonsFromKLProcess = new TH1F("gHistProtonsFromKLProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromKLProcess);
  TH1F *gHistProtonsFromPionProcess = new TH1F("gHistProtonsFromPionProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromPionProcess);
  TH1F *gHistProtonsFromKSProcess = new TH1F("gHistProtonsFromKSProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromKSProcess);
  TH1F *gHistProtonsFromKaonProcess = new TH1F("gHistProtonsFromKaonProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromKaonProcess);
  TH1F *gHistProtonsFromNeutronProcess = new TH1F("gHistProtonsFromNeutronProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromNeutronProcess);
  TH1F *gHistProtonsFromProtonProcess = new TH1F("gHistProtonsFromProtonProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromProtonProcess);
  TH1F *gHistProtonsFromSigmaMinusProcess = new TH1F("gHistProtonsFromSigmaMinusProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromSigmaMinusProcess);
  TH1F *gHistProtonsFromLambda0Process = new TH1F("gHistProtonsFromLambda0Process","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromLambda0Process);
  TH1F *gHistProtonsFromSigmaPlusProcess = new TH1F("gHistProtonsFromSigmaPlusProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromSigmaPlusProcess);
  TH1F *gHistProtonsFromXiMinusProcess = new TH1F("gHistProtonsFromXiMinusProcess","",51,-0.5,50.5);
  fMCProcessesList->Add(gHistProtonsFromXiMinusProcess);
  TH1F *gHistProtonsFromXi0Process = new TH1F("gHistProtonsFromXi0Process","",51,-0.5,50.5);					    
  fMCProcessesList->Add(gHistProtonsFromXi0Process);
  TH1F *gHistProtonsFromOmegaProcess = new TH1F("gHistProtonsFromOmegaProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistProtonsFromOmegaProcess);

  TH1F *gHistAntiProtonsFromKLProcess = new TH1F("gHistAntiProtonsFromKLProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromKLProcess);
  TH1F *gHistAntiProtonsFromPionProcess = new TH1F("gHistAntiProtonsFromPionProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromPionProcess);
  TH1F *gHistAntiProtonsFromKSProcess = new TH1F("gHistAntiProtonsFromKSProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromKSProcess);
  TH1F *gHistAntiProtonsFromKaonProcess = new TH1F("gHistAntiProtonsFromKaonProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromKaonProcess);
  TH1F *gHistAntiProtonsFromNeutronProcess = new TH1F("gHistAntiProtonsFromNeutronProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromNeutronProcess);
  TH1F *gHistAntiProtonsFromProtonProcess = new TH1F("gHistAntiProtonsFromProtonProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromProtonProcess);
  TH1F *gHistAntiProtonsFromLambda0Process = new TH1F("gHistAntiProtonsFromLambda0Process","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromLambda0Process);
  TH1F *gHistAntiProtonsFromSigmaPlusProcess = new TH1F("gHistAntiProtonsFromSigmaPlusProcess","",51,-0.5,50.5); 
  fMCProcessesList->Add(gHistAntiProtonsFromSigmaPlusProcess);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunMCAnalysis(AliStack* stack) {
  //Main analysis part - MC 
  for(Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if((Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;

    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) {
      if(iParticle <= stack->GetNprimary()) 
	((TH3F *)(fPDGList->At(0)))->Fill(Rapidity(particle->Px(),
						   particle->Py(),
						   particle->Pz()),
					  particle->Pt(),0);
      else if(iParticle > stack->GetNprimary()) {
	Int_t lPartMother = particle->GetFirstMother();
	TParticle *motherParticle = stack->Particle(lPartMother);
	if(!motherParticle) continue;
	Int_t motherPDGCode = motherParticle->GetPdgCode();
	if(fMCProcessIdFlag)
	  if(particle->GetUniqueID() != fMCProcessId) continue;
	if(fMotherParticlePDGCodeFlag)
	  if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	((TH3F *)(fPDGList->At(0)))->Fill(Rapidity(particle->Px(),
						   particle->Py(),
						   particle->Pz()),
					  particle->Pt(),
					  ConvertPDGToInt(motherParticle->GetPdgCode()));
	//processes
	if(TMath::Abs(motherParticle->GetPdgCode()) == 130)
	  ((TH1F *)(fMCProcessesList->At(0)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 211)
	  ((TH1F *)(fMCProcessesList->At(1)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 310)
	  ((TH1F *)(fMCProcessesList->At(2)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 321)
	  ((TH1F *)(fMCProcessesList->At(3)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 2112)
	  ((TH1F *)(fMCProcessesList->At(4)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 2212)
	  ((TH1F *)(fMCProcessesList->At(5)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3112)
	  ((TH1F *)(fMCProcessesList->At(6)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3122)
	  ((TH1F *)(fMCProcessesList->At(7)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3222)
	  ((TH1F *)(fMCProcessesList->At(8)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3312)
	  ((TH1F *)(fMCProcessesList->At(9)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3322)
	  ((TH1F *)(fMCProcessesList->At(10)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3334)
	  ((TH1F *)(fMCProcessesList->At(11)))->Fill(particle->GetUniqueID());
      }//secondary proton
    }//pdgcode of proton

    if(pdgcode == -2212) {
      if(iParticle <= stack->GetNprimary()) 
	((TH3F *)(fPDGList->At(1)))->Fill(Rapidity(particle->Px(),
						   particle->Py(),
						   particle->Pz()),
					  particle->Pt(),0);
      else if(iParticle > stack->GetNprimary()) {
	Int_t lPartMother = particle->GetFirstMother();
	TParticle *motherParticle = stack->Particle(lPartMother);
	if(!motherParticle) continue;
	Int_t motherPDGCode = motherParticle->GetPdgCode();
	if(fMCProcessIdFlag)
	  if(particle->GetUniqueID() != fMCProcessId) continue;
	if(fMotherParticlePDGCodeFlag)
	  if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	((TH3F *)(fPDGList->At(1)))->Fill(Rapidity(particle->Px(),
						   particle->Py(),
						   particle->Pz()),
					  particle->Pt(),
					  ConvertPDGToInt(motherParticle->GetPdgCode()));

	//processes
	if(TMath::Abs(motherParticle->GetPdgCode()) == 130)
	  ((TH1F *)(fMCProcessesList->At(12)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 211)
	  ((TH1F *)(fMCProcessesList->At(13)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 310)
	  ((TH1F *)(fMCProcessesList->At(14)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 321)
	  ((TH1F *)(fMCProcessesList->At(15)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 2112)
	  ((TH1F *)(fMCProcessesList->At(16)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 2212)
	  ((TH1F *)(fMCProcessesList->At(17)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3122)
	  ((TH1F *)(fMCProcessesList->At(18)))->Fill(particle->GetUniqueID());
	if(TMath::Abs(motherParticle->GetPdgCode()) == 3222)
	  ((TH1F *)(fMCProcessesList->At(19)))->Fill(particle->GetUniqueID());
      }//secondary antiproton
    }//pdgcode of antiproton

  }//particle loop
}

//____________________________________________________________________//
Int_t AliProtonQAAnalysis::ConvertPDGToInt(Int_t pdgCode) {
  //Converts the pdg code to an int based on the following scheme:
  //1: PDG code: 130 - Name: K_L0
  //2: PDG code: 211 - Name: pi+
  //3: PDG code: 310 - Name: K_S0
  //4: PDG code: 321 - Name: K+
  //5: PDG code: 2112 - Name: neutron
  //6: PDG code: 2212 - Name: proton
  //7: PDG code: 3112 - Name: Sigma-
  //8: PDG code: 3122 - Name: Lambda0
  //9: PDG code: 3222 - Name: Sigma+
  //10: PDG code: 3312 - Name: Xi-
  //11: PDG code: 3322 - Name: Xi0
  //12: PDG code: 3334 - Name: Omega-
  Int_t code = -1;
  switch (TMath::Abs(pdgCode)) {
  case 130: {
    code = 1;
    break;
  }
  case 211: {
    code = 2;
    break;
  }
  case 310: {
    code = 3;
    break;
  }
  case 321: {
    code = 4;
    break;
  }
  case 2112: {
    code = 5;
    break;
  }
  case 2212: {
    code = 6;
    break;
  }
  case 3112: {
    code = 7;
    break;
  }
  case 3122: {
    code = 8;
    break;
  }
  case 3222: {
    code = 9;
    break;
  }
  case 3312: {
    code = 10;
    break;
  }
  case 3322: {
    code = 11;
    break;
  }
  case 3334: {
    code = 12;
    break;
  }
  default: {
    code = -1;
    break;
  }
  }//switch

  return code;
}








