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
  fUseTPCOnly(kFALSE),
  fPDGList(0), fMCProcessesList(0),
  fRunMCAnalysis(kFALSE) {
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
  if(fUseTPCOnly) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
    }
  }
  else{
    Pt = track->Pt();
    Px = track->Px();
    Py = track->Py();
    Pz = track->Pz();
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
  if(fUseTPCOnly) {
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
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
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
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
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
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
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
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
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
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
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
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
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
	if(dcaXY > fMaxDCAXY) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(12)))->Fill(dcaXY);
	  //status = kFALSE;
	}
	else if(dcaXY <= fMaxDCAXY)
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
	if(dcaZ > fMaxDCAZ) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(14)))->Fill(dcaZ);
	  //status = kFALSE;
	}
	else if(dcaZ <= fMaxDCAZ)
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
    }//secondary particle cut
  }//antiprotons

  //if((Pt < fMinPt) || (Pt > fMaxPt)) //status = kFALSE;
  //if((Rapidity(Px,Py,Pz) < fMinY) || (Rapidity(Px,Py,Pz) > fMaxY)) 
  //status = kFALSE;

  //return status;
}

//____________________________________________________________________//
Float_t AliProtonQAAnalysis::GetSigmaToVertex(AliESDtrack* esdTrack) {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  if(fUseTPCOnly) 
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
void AliProtonQAAnalysis::SetQAOn() {
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
  if(fRunMCAnalysis) InitMCAnalysis();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitQA() {
  //Initializes the QA histograms and builds the directory structure
  //if(!fQAHistograms) 
  SetQAOn();

  //2D histograms
  //TDirectory *dir2D = gDirectory->mkdir("2D");
  //fGlobalQAList->Add(dir2D); dir2D->cd();
  TH2D *gHistYPtPrimaryProtonsPass = new TH2D("gHistYPtPrimaryProtonsPass",
					      ";y;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsPass);
  TH2D *gHistYPtPrimaryProtonsReject = new TH2D("gHistYPtPrimaryProtonsReject",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsReject);

  TH2D *gHistYPtSecondaryProtonsPass = new TH2D("gHistYPtSecondaryProtonsPass",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsPass);
  TH2D *gHistYPtSecondaryProtonsReject = new TH2D("gHistYPtSecondaryProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsReject);

  TH2D *gHistYPtPrimaryAntiProtonsPass = new TH2D("gHistYPtPrimaryAntiProtonsPass",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsPass);
  TH2D *gHistYPtPrimaryAntiProtonsReject = new TH2D("gHistYPtPrimaryAntiProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsReject);

  TH2D *gHistYPtSecondaryAntiProtonsPass = new TH2D("gHistYPtSecondaryAntiProtonsPass",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsPass);
  TH2D *gHistYPtSecondaryAntiProtonsReject = new TH2D("gHistYPtSecondaryAntiProtonsReject",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsReject);

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
  TH1F *fPrimaryProtonsITSClustersPass = new TH1F("fPrimaryProtonsITSClustersPass",
					    ";N_{clusters} (ITS);Entries",
					    7,0,7);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsITSClustersPass);
  TH1F *fPrimaryProtonsChi2PerClusterITSPass = new TH1F("fPrimaryProtonsChi2PerClusterITSPass",
						  ";x^{2}/N_{clusters} (ITS);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsChi2PerClusterITSPass);
  TH1F *fPrimaryProtonsTPCClustersPass = new TH1F("fPrimaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCClustersPass);
  TH1F *fPrimaryProtonsChi2PerClusterTPCPass = new TH1F("fPrimaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsChi2PerClusterTPCPass);
  TH1F *fPrimaryProtonsExtCov11Pass = new TH1F("fPrimaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov11Pass);
  TH1F *fPrimaryProtonsExtCov22Pass = new TH1F("fPrimaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov22Pass);
  TH1F *fPrimaryProtonsExtCov33Pass = new TH1F("fPrimaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov33Pass);
  TH1F *fPrimaryProtonsExtCov44Pass = new TH1F("fPrimaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov44Pass);
  TH1F *fPrimaryProtonsExtCov55Pass = new TH1F("fPrimaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov55Pass);
  TH1F *fPrimaryProtonsSigmaToVertexPass = new TH1F("fPrimaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsSigmaToVertexPass);
  TH1F *fPrimaryProtonsSigmaToVertexTPCPass = new TH1F("fPrimaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsSigmaToVertexTPCPass);
  TH1F *fPrimaryProtonsDCAXYPass = new TH1F("fPrimaryProtonsDCAXYPass",
					     ";DCA_{xy} [cm];Entries",
					     100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsDCAXYPass);
  TH1F *fPrimaryProtonsDCAXYTPCPass = new TH1F("fPrimaryProtonsDCAXYTPCPass",
					       ";DCA_{xy} [cm];Entries",
					       100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsDCAXYTPCPass);
  TH1F *fPrimaryProtonsDCAZPass = new TH1F("fPrimaryProtonsDCAZPass",
					   ";DCA_{z} [cm];Entries",
					   100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsDCAZPass);
  TH1F *fPrimaryProtonsDCAZTPCPass = new TH1F("fPrimaryProtonsDCAZTPCPass",
					      ";DCA_{z} [cm];Entries",
					      100,0,20);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsDCAZTPCPass);
  TH1F *fPrimaryProtonsConstrainChi2Pass = new TH1F("fPrimaryProtonsConstrainChi2Pass",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsConstrainChi2Pass);
  TH1F *fPrimaryProtonsITSRefitPass = new TH1F("fPrimaryProtonsITSRefitPass",
					       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsITSRefitPass);
  TH1F *fPrimaryProtonsTPCRefitPass = new TH1F("fPrimaryProtonsTPCRefitPass",
					       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCRefitPass);
  TH1F *fPrimaryProtonsESDpidPass = new TH1F("fPrimaryProtonsESDpidPass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsESDpidPass);
  TH1F *fPrimaryProtonsTPCpidPass = new TH1F("fPrimaryProtonsTPCpidPass",
					     "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCpidPass);

  //Rejected primary protons
  /*gDirectory->cd("../");
  TDirectory *dirProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsPrimaryRejected->cd();*/

  TH1F *fPrimaryProtonsITSClustersReject = new TH1F("fPrimaryProtonsITSClustersReject",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsITSClustersReject);
  TH1F *fPrimaryProtonsChi2PerClusterITSReject = new TH1F("fPrimaryProtonsChi2PerClusterITSReject",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsChi2PerClusterITSReject);
  TH1F *fPrimaryProtonsTPCClustersReject = new TH1F("fPrimaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCClustersReject);
  TH1F *fPrimaryProtonsChi2PerClusterTPCReject = new TH1F("fPrimaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsChi2PerClusterTPCReject);
  TH1F *fPrimaryProtonsExtCov11Reject = new TH1F("fPrimaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov11Reject);
  TH1F *fPrimaryProtonsExtCov22Reject = new TH1F("fPrimaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov22Reject);
  TH1F *fPrimaryProtonsExtCov33Reject = new TH1F("fPrimaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov33Reject);
  TH1F *fPrimaryProtonsExtCov44Reject = new TH1F("fPrimaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov44Reject);
  TH1F *fPrimaryProtonsExtCov55Reject = new TH1F("fPrimaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov55Reject);
  TH1F *fPrimaryProtonsSigmaToVertexReject = new TH1F("fPrimaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsSigmaToVertexReject);
  TH1F *fPrimaryProtonsSigmaToVertexTPCReject = new TH1F("fPrimaryProtonsSigmaToVertexTPCReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsSigmaToVertexTPCReject);
  TH1F *fPrimaryProtonsDCAXYReject = new TH1F("fPrimaryProtonsDCAXYReject",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsDCAXYReject);
  TH1F *fPrimaryProtonsDCAXYTPCReject = new TH1F("fPrimaryProtonsDCAXYTPCReject",
						 ";DCA_{xy} [cm];Entries",
						 100,0,20);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsDCAXYTPCReject);
  TH1F *fPrimaryProtonsDCAZReject = new TH1F("fPrimaryProtonsDCAZReject",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsDCAZReject);
  TH1F *fPrimaryProtonsDCAZTPCReject = new TH1F("fPrimaryProtonsDCAZTPCReject",
						";DCA_{z} [cm];Entries",
						100,0,20);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsDCAZTPCReject);
  TH1F *fPrimaryProtonsConstrainChi2Reject = new TH1F("fPrimaryProtonsConstrainChi2Reject",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsConstrainChi2Reject);
  TH1F *fPrimaryProtonsITSRefitReject = new TH1F("fPrimaryProtonsITSRefitReject",
						 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsITSRefitReject);
  TH1F *fPrimaryProtonsTPCRefitReject = new TH1F("fPrimaryProtonsTPCRefitReject",
						 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCRefitReject);
  TH1F *fPrimaryProtonsESDpidReject = new TH1F("fPrimaryProtonsESDpidReject",
					       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsESDpidReject);
  TH1F *fPrimaryProtonsTPCpidReject = new TH1F("fPrimaryProtonsTPCpidReject",
					       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCpidReject);

  //________________________________________________________________//
  /*gDirectory->cd("../../");

  TDirectory *dirProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirProtonsSecondary->cd();
  TDirectory *dirProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirProtonsSecondaryAccepted->cd();*/

  //Accepted secondary protons
  TH1F *fSecondaryProtonsITSClustersPass = new TH1F("fSecondaryProtonsITSClustersPass",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsITSClustersPass);
  TH1F *fSecondaryProtonsChi2PerClusterITSPass = new TH1F("fSecondaryProtonsChi2PerClusterITSPass",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsChi2PerClusterITSPass);
  TH1F *fSecondaryProtonsTPCClustersPass = new TH1F("fSecondaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCClustersPass);
  TH1F *fSecondaryProtonsChi2PerClusterTPCPass = new TH1F("fSecondaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsChi2PerClusterTPCPass);
  TH1F *fSecondaryProtonsExtCov11Pass = new TH1F("fSecondaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov11Pass);
  TH1F *fSecondaryProtonsExtCov22Pass = new TH1F("fSecondaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov22Pass);
  TH1F *fSecondaryProtonsExtCov33Pass = new TH1F("fSecondaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov33Pass);
  TH1F *fSecondaryProtonsExtCov44Pass = new TH1F("fSecondaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov44Pass);
  TH1F *fSecondaryProtonsExtCov55Pass = new TH1F("fSecondaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov55Pass);
  TH1F *fSecondaryProtonsSigmaToVertexPass = new TH1F("fSecondaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsSigmaToVertexPass);
  TH1F *fSecondaryProtonsSigmaToVertexTPCPass = new TH1F("fSecondaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsSigmaToVertexTPCPass);
  TH1F *fSecondaryProtonsDCAXYPass = new TH1F("fSecondaryProtonsDCAXYPass",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsDCAXYPass);
  TH1F *fSecondaryProtonsDCAXYTPCPass = new TH1F("fSecondaryProtonsDCAXYTPCPass",
						 ";DCA_{xy} [cm];Entries",
						 100,0,20);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsDCAXYTPCPass);
  TH1F *fSecondaryProtonsDCAZPass = new TH1F("fSecondaryProtonsDCAZPass",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsDCAZPass);
  TH1F *fSecondaryProtonsDCAZTPCPass = new TH1F("fSecondaryProtonsDCAZTPCPass",
						";DCA_{z} [cm];Entries",
						100,0,20);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsDCAZTPCPass);
  TH1F *fSecondaryProtonsConstrainChi2Pass = new TH1F("fSecondaryProtonsConstrainChi2Pass",
						    ";Log_{10}(#chi^{2});Entries",
						    100,-10,10);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsConstrainChi2Pass);
  TH1F *fSecondaryProtonsITSRefitPass = new TH1F("fSecondaryProtonsITSRefitPass",
						 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsITSRefitPass);
  TH1F *fSecondaryProtonsTPCRefitPass = new TH1F("fSecondaryProtonsTPCRefitPass",
						 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCRefitPass);
  TH1F *fSecondaryProtonsESDpidPass = new TH1F("fSecondaryProtonsESDpidPass",
					       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsESDpidPass);
  TH1F *fSecondaryProtonsTPCpidPass = new TH1F("fSecondaryProtonsTPCpidPass",
					       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCpidPass);

  //Rejected secondary protons
  /*gDirectory->cd("../");
  TDirectory *dirProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsSecondaryRejected->cd();*/

  TH1F *fSecondaryProtonsITSClustersReject = new TH1F("fSecondaryProtonsITSClustersReject",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsITSClustersReject);
  TH1F *fSecondaryProtonsChi2PerClusterITSReject = new TH1F("fSecondaryProtonsChi2PerClusterITSReject",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsChi2PerClusterITSReject);
  TH1F *fSecondaryProtonsTPCClustersReject = new TH1F("fSecondaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCClustersReject);
  TH1F *fSecondaryProtonsChi2PerClusterTPCReject = new TH1F("fSecondaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsChi2PerClusterTPCReject);
  TH1F *fSecondaryProtonsExtCov11Reject = new TH1F("fSecondaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov11Reject);
  TH1F *fSecondaryProtonsExtCov22Reject = new TH1F("fSecondaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov22Reject);
  TH1F *fSecondaryProtonsExtCov33Reject = new TH1F("fSecondaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov33Reject);
  TH1F *fSecondaryProtonsExtCov44Reject = new TH1F("fSecondaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov44Reject);
  TH1F *fSecondaryProtonsExtCov55Reject = new TH1F("fSecondaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov55Reject);
  TH1F *fSecondaryProtonsSigmaToVertexReject = new TH1F("fSecondaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsSigmaToVertexReject);
  TH1F *fSecondaryProtonsSigmaToVertexTPCReject = new TH1F("fSecondaryProtonsSigmaToVertexTPCReject",
							   ";#sigma_{Vertex};Entries",
							   100,0,10);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsSigmaToVertexTPCReject);
  TH1F *fSecondaryProtonsDCAXYReject = new TH1F("fSecondaryProtonsDCAXYReject",
						";DCA_{xy} [cm];Entries",
						100,0,20);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsDCAXYReject);
  TH1F *fSecondaryProtonsDCAXYTPCReject = new TH1F("fSecondaryProtonsDCAXYTPCReject",
						   ";DCA_{xy} [cm];Entries",
						   100,0,20);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsDCAXYTPCReject);
  TH1F *fSecondaryProtonsDCAZReject = new TH1F("fSecondaryProtonsDCAZReject",
					       ";DCA_{z} [cm];Entries",
					       100,0,20);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsDCAZReject);
  TH1F *fSecondaryProtonsDCAZTPCReject = new TH1F("fSecondaryProtonsDCAZTPCReject",
						  ";DCA_{z} [cm];Entries",
						  100,0,20);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsDCAZTPCReject);
  TH1F *fSecondaryProtonsConstrainChi2Reject = new TH1F("fSecondaryProtonsConstrainChi2Reject",
							";Log_{10}(#chi^{2});Entries",
							100,-10,10);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsConstrainChi2Reject);
  TH1F *fSecondaryProtonsITSRefitReject = new TH1F("fSecondaryProtonsITSRefitReject",
						   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsITSRefitReject);
  TH1F *fSecondaryProtonsTPCRefitReject = new TH1F("fSecondaryProtonsTPCRefitReject",
						   "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCRefitReject);
  TH1F *fSecondaryProtonsESDpidReject = new TH1F("fSecondaryProtonsESDpidReject",
						 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsESDpidReject);
  TH1F *fSecondaryProtonsTPCpidReject = new TH1F("fSecondaryProtonsTPCpidReject",
						 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCpidReject);
  

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
  TH1F *fPrimaryAntiProtonsITSClustersPass = new TH1F("fPrimaryAntiProtonsITSClustersPass",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsITSClustersPass);
  TH1F *fPrimaryAntiProtonsChi2PerClusterITSPass = new TH1F("fPrimaryAntiProtonsChi2PerClusterITSPass",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsChi2PerClusterITSPass);
  TH1F *fPrimaryAntiProtonsTPCClustersPass = new TH1F("fPrimaryAntiProtonsTPCClustersPass",
						      ";N_{clusters} (TPC);Entries",
						      100,0,200);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCClustersPass);
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCPass = new TH1F("fPrimaryAntiProtonsChi2PerClusterTPCPass",
							    ";x^{2}/N_{clusters} (TPC);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *fPrimaryAntiProtonsExtCov11Pass = new TH1F("fPrimaryAntiProtonsExtCov11Pass",
						   ";#sigma_{y} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov11Pass);
  TH1F *fPrimaryAntiProtonsExtCov22Pass = new TH1F("fPrimaryAntiProtonsExtCov22Pass",
						   ";#sigma_{z} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov22Pass);
  TH1F *fPrimaryAntiProtonsExtCov33Pass = new TH1F("fPrimaryAntiProtonsExtCov33Pass",
						   ";#sigma_{sin(#phi)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov33Pass);
  TH1F *fPrimaryAntiProtonsExtCov44Pass = new TH1F("fPrimaryAntiProtonsExtCov44Pass",
						   ";#sigma_{tan(#lambda)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov44Pass);
  TH1F *fPrimaryAntiProtonsExtCov55Pass = new TH1F("fPrimaryAntiProtonsExtCov55Pass",
						   ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov55Pass);
  TH1F *fPrimaryAntiProtonsSigmaToVertexPass = new TH1F("fPrimaryAntiProtonsSigmaToVertexPass",
							";#sigma_{Vertex};Entries",
							100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsSigmaToVertexPass);
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCPass = new TH1F("fPrimaryAntiProtonsSigmaToVertexTPCPass",
							   ";#sigma_{Vertex};Entries",
							   100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *fPrimaryAntiProtonsDCAXYPass = new TH1F("fPrimaryAntiProtonsDCAXYPass",
						";DCA_{xy} [cm];Entries",
						100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsDCAXYPass);
  TH1F *fPrimaryAntiProtonsDCAXYTPCPass = new TH1F("fPrimaryAntiProtonsDCAXYTPCPass",
						   ";DCA_{xy} [cm];Entries",
						   100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsDCAXYTPCPass);
  TH1F *fPrimaryAntiProtonsDCAZPass = new TH1F("fPrimaryAntiProtonsDCAZPass",
					       ";DCA_{z} [cm];Entries",
					       100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsDCAZPass);
  TH1F *fPrimaryAntiProtonsDCAZTPCPass = new TH1F("fPrimaryAntiProtonsDCAZTPCPass",
						  ";DCA_{z} [cm];Entries",
						  100,0,20);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsDCAZTPCPass);
  TH1F *fPrimaryAntiProtonsConstrainChi2Pass = new TH1F("fPrimaryAntiProtonsConstrainChi2Pass",
							";Log_{10}(#chi^{2});Entries",
							100,-10,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsConstrainChi2Pass);
  TH1F *fPrimaryAntiProtonsITSRefitPass = new TH1F("fPrimaryAntiProtonsITSRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsITSRefitPass);
  TH1F *fPrimaryAntiProtonsTPCRefitPass = new TH1F("fPrimaryAntiProtonsTPCRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCRefitPass);
  TH1F *fPrimaryAntiProtonsESDpidPass = new TH1F("fPrimaryAntiProtonsESDpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsESDpidPass);
  TH1F *fPrimaryAntiProtonsTPCpidPass = new TH1F("fPrimaryAntiProtonsTPCpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCpidPass);
  
  //Rejected primary antiprotons
  /*gDirectory->cd("../");
  TDirectory *dirAntiProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsPrimaryRejected->cd();*/
  
  TH1F *fPrimaryAntiProtonsITSClustersReject = new TH1F("fPrimaryAntiProtonsITSClustersReject",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsITSClustersReject);
  TH1F *fPrimaryAntiProtonsChi2PerClusterITSReject = new TH1F("fPrimaryAntiProtonsChi2PerClusterITSReject",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsChi2PerClusterITSReject);
  TH1F *fPrimaryAntiProtonsTPCClustersReject = new TH1F("fPrimaryAntiProtonsTPCClustersReject",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCClustersReject);
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCReject = new TH1F("fPrimaryAntiProtonsChi2PerClusterTPCReject",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *fPrimaryAntiProtonsExtCov11Reject = new TH1F("fPrimaryAntiProtonsExtCov11Reject",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov11Reject);
  TH1F *fPrimaryAntiProtonsExtCov22Reject = new TH1F("fPrimaryAntiProtonsExtCov22Reject",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov22Reject);
  TH1F *fPrimaryAntiProtonsExtCov33Reject = new TH1F("fPrimaryAntiProtonsExtCov33Reject",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov33Reject);
  TH1F *fPrimaryAntiProtonsExtCov44Reject = new TH1F("fPrimaryAntiProtonsExtCov44Reject",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov44Reject);
  TH1F *fPrimaryAntiProtonsExtCov55Reject = new TH1F("fPrimaryAntiProtonsExtCov55Reject",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov55Reject);
  TH1F *fPrimaryAntiProtonsSigmaToVertexReject = new TH1F("fPrimaryAntiProtonsSigmaToVertexReject",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsSigmaToVertexReject);
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCReject = new TH1F("fPrimaryAntiProtonsSigmaToVertexTPCReject",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *fPrimaryAntiProtonsDCAXYReject = new TH1F("fPrimaryAntiProtonsDCAXYReject",
						  ";DCA_{xy} [cm];Entries",
						  100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsDCAXYReject);
  TH1F *fPrimaryAntiProtonsDCAXYTPCReject = new TH1F("fPrimaryAntiProtonsDCAXYTPCReject",
						     ";DCA_{xy} [cm];Entries",
						     100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsDCAXYTPCReject);
  TH1F *fPrimaryAntiProtonsDCAZReject = new TH1F("fPrimaryAntiProtonsDCAZReject",
						 ";DCA_{z} [cm];Entries",
						 100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsDCAZReject);
  TH1F *fPrimaryAntiProtonsDCAZTPCReject = new TH1F("fPrimaryAntiProtonsDCAZTPCReject",
						    ";DCA_{z} [cm];Entries",
						    100,0,20);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsDCAZTPCReject);
  TH1F *fPrimaryAntiProtonsConstrainChi2Reject = new TH1F("fPrimaryAntiProtonsConstrainChi2Reject",
							  ";Log_{10}(#chi^{2});Entries",
							  100,-10,10);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsConstrainChi2Reject);
  TH1F *fPrimaryAntiProtonsITSRefitReject = new TH1F("fPrimaryAntiProtonsITSRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsITSRefitReject);
  TH1F *fPrimaryAntiProtonsTPCRefitReject = new TH1F("fPrimaryAntiProtonsTPCRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCRefitReject);
  TH1F *fPrimaryAntiProtonsESDpidReject = new TH1F("fPrimaryAntiProtonsESDpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsESDpidReject);
  TH1F *fPrimaryAntiProtonsTPCpidReject = new TH1F("fPrimaryAntiProtonsTPCpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCpidReject);
  
  //________________________________________________________________//
  /*gDirectory->cd("../../");

  TDirectory *dirAntiProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirAntiProtonsSecondary->cd();
  TDirectory *dirAntiProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirAntiProtonsSecondaryAccepted->cd();*/

  //Accepted secondary antiprotons
  TH1F *fSecondaryAntiProtonsITSClustersPass = new TH1F("fSecondaryAntiProtonsITSClustersPass",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsITSClustersPass);
  TH1F *fSecondaryAntiProtonsChi2PerClusterITSPass = new TH1F("fSecondaryAntiProtonsChi2PerClusterITSPass",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsChi2PerClusterITSPass);
  TH1F *fSecondaryAntiProtonsTPCClustersPass = new TH1F("fSecondaryAntiProtonsTPCClustersPass",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCClustersPass);
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCPass = new TH1F("fSecondaryAntiProtonsChi2PerClusterTPCPass",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *fSecondaryAntiProtonsExtCov11Pass = new TH1F("fSecondaryAntiProtonsExtCov11Pass",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov11Pass);
  TH1F *fSecondaryAntiProtonsExtCov22Pass = new TH1F("fSecondaryAntiProtonsExtCov22Pass",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov22Pass);
  TH1F *fSecondaryAntiProtonsExtCov33Pass = new TH1F("fSecondaryAntiProtonsExtCov33Pass",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov33Pass);
  TH1F *fSecondaryAntiProtonsExtCov44Pass = new TH1F("fSecondaryAntiProtonsExtCov44Pass",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov44Pass);
  TH1F *fSecondaryAntiProtonsExtCov55Pass = new TH1F("fSecondaryAntiProtonsExtCov55Pass",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov55Pass);
  TH1F *fSecondaryAntiProtonsSigmaToVertexPass = new TH1F("fSecondaryAntiProtonsSigmaToVertexPass",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsSigmaToVertexPass);
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCPass = new TH1F("fSecondaryAntiProtonsSigmaToVertexTPCPass",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *fSecondaryAntiProtonsDCAXYPass = new TH1F("fSecondaryAntiProtonsDCAXYPass",
						  ";DCA_{xy} [cm];Entries",
						  100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsDCAXYPass);
  TH1F *fSecondaryAntiProtonsDCAXYTPCPass = new TH1F("fSecondaryAntiProtonsDCAXYTPCPass",
						     ";DCA_{xy} [cm];Entries",
						     100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsDCAXYTPCPass);
  TH1F *fSecondaryAntiProtonsDCAZPass = new TH1F("fSecondaryAntiProtonsDCAZPass",
						 ";DCA_{z} [cm];Entries",
						 100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsDCAZPass);
  TH1F *fSecondaryAntiProtonsDCAZTPCPass = new TH1F("fSecondaryAntiProtonsDCAZTPCPass",
						    ";DCA_{z} [cm];Entries",
						    100,0,20);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsDCAZTPCPass);
  TH1F *fSecondaryAntiProtonsConstrainChi2Pass = new TH1F("fSecondaryAntiProtonsConstrainChi2Pass",
							  ";Log_{10}(#chi^{2});Entries",
							  100,-10,10);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsConstrainChi2Pass);
  TH1F *fSecondaryAntiProtonsITSRefitPass = new TH1F("fSecondaryAntiProtonsITSRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsITSRefitPass);
  TH1F *fSecondaryAntiProtonsTPCRefitPass = new TH1F("fSecondaryAntiProtonsTPCRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCRefitPass);
  TH1F *fSecondaryAntiProtonsESDpidPass = new TH1F("fSecondaryAntiProtonsESDpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsESDpidPass);
  TH1F *fSecondaryAntiProtonsTPCpidPass = new TH1F("fSecondaryAntiProtonsTPCpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCpidPass);
  
  //Rejected secondary antiprotons
  /*gDirectory->cd("../");
  TDirectory *dirAntiProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsSecondaryRejected->cd();*/

  TH1F *fSecondaryAntiProtonsITSClustersReject = new TH1F("fSecondaryAntiProtonsITSClustersReject",
							  ";N_{clusters} (ITS);Entries",
							  7,0,7);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsITSClustersReject);
  TH1F *fSecondaryAntiProtonsChi2PerClusterITSReject = new TH1F("fSecondaryAntiProtonsChi2PerClusterITSReject",
								";x^{2}/N_{clusters} (ITS);Entries",
								100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsChi2PerClusterITSReject);
  TH1F *fSecondaryAntiProtonsTPCClustersReject = new TH1F("fSecondaryAntiProtonsTPCClustersReject",
							  ";N_{clusters} (TPC);Entries",
							  100,0,200);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCClustersReject);
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCReject = new TH1F("fSecondaryAntiProtonsChi2PerClusterTPCReject",
								";x^{2}/N_{clusters} (TPC);Entries",
								100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *fSecondaryAntiProtonsExtCov11Reject = new TH1F("fSecondaryAntiProtonsExtCov11Reject",
						       ";#sigma_{y} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov11Reject);
  TH1F *fSecondaryAntiProtonsExtCov22Reject = new TH1F("fSecondaryAntiProtonsExtCov22Reject",
						       ";#sigma_{z} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov22Reject);
  TH1F *fSecondaryAntiProtonsExtCov33Reject = new TH1F("fSecondaryAntiProtonsExtCov33Reject",
						       ";#sigma_{sin(#phi)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov33Reject);
  TH1F *fSecondaryAntiProtonsExtCov44Reject = new TH1F("fSecondaryAntiProtonsExtCov44Reject",
						       ";#sigma_{tan(#lambda)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov44Reject);
  TH1F *fSecondaryAntiProtonsExtCov55Reject = new TH1F("fSecondaryAntiProtonsExtCov55Reject",
						       ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov55Reject);
  TH1F *fSecondaryAntiProtonsSigmaToVertexReject = new TH1F("fSecondaryAntiProtonsSigmaToVertexReject",
							    ";#sigma_{Vertex};Entries",
							    100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsSigmaToVertexReject);
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCReject = new TH1F("fSecondaryAntiProtonsSigmaToVertexTPCReject",
							       ";#sigma_{Vertex};Entries",
							       100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *fSecondaryAntiProtonsDCAXYReject = new TH1F("fSecondaryAntiProtonsDCAXYReject",
						    ";DCA_{xy} [cm];Entries",
						    100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsDCAXYReject);
  TH1F *fSecondaryAntiProtonsDCAXYTPCReject = new TH1F("fSecondaryAntiProtonsDCAXYTPCReject",
						       ";DCA_{xy} [cm];Entries",
						       100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsDCAXYTPCReject);
  TH1F *fSecondaryAntiProtonsDCAZReject = new TH1F("fSecondaryAntiProtonsDCAZReject",
						   ";DCA_{z} [cm];Entries",
						   100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsDCAZReject);
  TH1F *fSecondaryAntiProtonsDCAZTPCReject = new TH1F("fSecondaryAntiProtonsDCAZTPCReject",
						      ";DCA_{z} [cm];Entries",
						      100,0,20);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsDCAZTPCReject);
  TH1F *fSecondaryAntiProtonsConstrainChi2Reject = new TH1F("fSecondaryAntiProtonsConstrainChi2Reject",
							    ";Log_{10}(#chi^{2});Entries",
							    100,-10,10);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsConstrainChi2Reject);
  TH1F *fSecondaryAntiProtonsITSRefitReject = new TH1F("fSecondaryAntiProtonsITSRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsITSRefitReject);
  TH1F *fSecondaryAntiProtonsTPCRefitReject = new TH1F("fSecondaryAntiProtonsTPCRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCRefitReject);
  TH1F *fSecondaryAntiProtonsESDpidReject = new TH1F("fSecondaryAntiProtonsESDpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsESDpidReject);
  TH1F *fSecondaryAntiProtonsTPCpidReject = new TH1F("fSecondaryAntiProtonsTPCpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCpidReject);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunQA(AliStack *stack, AliESDEvent *fESD) {
  //Runs the QA code
  Int_t nGoodTracks = fESD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    Int_t label = TMath::Abs(track->GetLabel()); 
    Double_t Pt = 0.0, P = 0.0;
    Double_t probability[5];

    if(fUseTPCOnly) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      Pt = tpcTrack->Pt();
      P = tpcTrack->P();
      
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
            if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(0)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(4)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(2)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(6)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//secondary particles
	}//cuts
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
	}//cuts
      }//proton check
    }//TPC only tracks
    else if(!fUseTPCOnly) {
      Pt = track->Pt();
      P = track->P();
      
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
            if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(0)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(4)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0)
              ((TH2D *)(fQA2DList->At(2)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
            else if(track->Charge() < 0)
              ((TH2D *)(fQA2DList->At(6)))->Fill(Rapidity(track->Px(),
							  track->Py(),
							  track->Pz()),
						 Pt);
	  }//secondary particles
	}//cuts
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
	}//cuts
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
    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
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








