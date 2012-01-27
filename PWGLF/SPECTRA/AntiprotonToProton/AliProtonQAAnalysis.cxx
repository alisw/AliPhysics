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
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TArrayI.h>
#include <TParticle.h>

#include "AliProtonQAAnalysis.h"
#include "AliProtonAnalysisBase.h"

#include <AliExternalTrackParam.h>
#include <AliESDEvent.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliESDVertex.h>
#include <AliGenEventHeader.h>
#include <AliMCEvent.h>

ClassImp(AliProtonQAAnalysis)

//____________________________________________________________________//
AliProtonQAAnalysis::AliProtonQAAnalysis() : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(0), fMinY(0), fMaxY(0), fY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0), fPt(0), fUseAsymmetricBinning(kFALSE),
  fGlobalQAList(0), fQAVertexList(0), fQA2DList(0),
  fQAPrimaryProtonsAcceptedList(0),
  fQAPrimaryProtonsRejectedList(0),
  fQASecondaryProtonsAcceptedList(0),
  fQASecondaryProtonsRejectedList(0),
  fQAPrimaryAntiProtonsAcceptedList(0),
  fQAPrimaryAntiProtonsRejectedList(0),
  fQASecondaryAntiProtonsAcceptedList(0),
  fQASecondaryAntiProtonsRejectedList(0),
  fPDGList(0), fMCProcessesList(0),
  fRunMCAnalysis(kFALSE),
  fMCProcessIdFlag(kFALSE), fMCProcessId(0),
  fMotherParticlePDGCodeFlag(kFALSE), fMotherParticlePDGCode(0),
  fAcceptedCutList(0), fRejectedCutList(0),
  fAcceptedDCAList(0), fRejectedDCAList(0),
  fRunEfficiencyAnalysis(kFALSE),
  fUseCutsInEfficiency(kFALSE),
  fEfficiencyList(0), fCutEfficiencyList(0) {
  //Default constructor
}

//____________________________________________________________________//
AliProtonQAAnalysis::~AliProtonQAAnalysis() {
  //Default destructor
  if(fProtonAnalysisBase) delete fProtonAnalysisBase;
  if(fGlobalQAList) delete fGlobalQAList;
  if(fQAVertexList) delete fQAVertexList;
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
  if(fCutEfficiencyList) delete fCutEfficiencyList;
}

//____________________________________________________________________//
void AliProtonQAAnalysis::FillQA(AliStack *const stack,
				 AliESDEvent *esd,
				 const AliESDVertex *vertex, 
				 AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  Int_t nPrimaries = stack->GetNprimary();
  Int_t label = TMath::Abs(track->GetLabel());

  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.

  if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
      tpcTrack->PropagateToDCA(vertex,
			       esd->GetMagneticField(),
			       100.,dca,cov);
    }
  }
  else{
    gPt = track->Pt();
    gPx = track->Px();
    gPy = track->Py();
    gPz = track->Pz();
    track->PropagateToDCA(vertex,
			  esd->GetMagneticField(),
			  100.,dca,cov);
  }

  //TParticle *particle = stack->Particle(label);
  //if(particle) {
  //Int_t pdgcode = particle->GetPdgCode();
  //Int_t gProcess = particle->GetUniqueID();
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
    
    //protons
    if(track->Charge() > 0) {
      //Primaries
      if(label <= nPrimaries) {
	if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
	  if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  }
	  else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
	}//ITS clusters
	if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
	  if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  }
	  else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
	}//chi2 per ITS cluster
	if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
	  if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  }
	  else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	  }
	}//TPC clusters
	if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
	  if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  }
	  else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
	}//chi2 per TPC cluster
	if(fProtonAnalysisBase->IsUsedMaxCov11()) {
	  if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  }
	  else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov22()) {
	  if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  }
	  else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov33()) {
	  if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  }
	  else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov44()) {
	  if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  }
	  else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov55()) {
	  if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  }
	  else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
	}//cov55
	if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
	  if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	  }
	  else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}//sigma to vertex
	if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
	  if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	  }
	  else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}//sigma to vertex TPC
	if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
	  if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
	  }
	  else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
	}//DCA xy global tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
	  if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
	  }
	  else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
	}//DCA xy TPC tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
	  if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
	  }
	  else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
	}//DCA z global tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
	  if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
	  }
	  else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
	}//DCA z TPC tracking
	if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
	  if(track->GetConstrainedChi2() > 0) {
	    if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	      ((TH1F *)(fQAPrimaryProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    }
	    else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	      ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	  }
	}//constrain chi2 - vertex
	if(fProtonAnalysisBase->IsUsedITSRefit()) {
	  if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(16)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(16)))->Fill(0);
	}//ITS refit
	if(fProtonAnalysisBase->IsUsedTPCRefit()) {
	  if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(17)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(17)))->Fill(0);
	}//TPC refit
	if(fProtonAnalysisBase->IsUsedESDpid()) {
	  if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(18)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(18)))->Fill(0);
	}//ESD pid
	if(fProtonAnalysisBase->IsUsedTPCpid()) {
	  if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(19)))->Fill(0);
	}
	  else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(19)))->Fill(0);
	}//TPC pid
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
	  if(!track->HasPointOnITSLayer(0)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(20)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(0))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(20)))->Fill(0);
	}//point on SPD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
	  if(!track->HasPointOnITSLayer(1)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(21)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(1))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(21)))->Fill(0);
	}//point on SPD2
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
	  if(!track->HasPointOnITSLayer(2)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(22)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(2))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(22)))->Fill(0);
	}//point on SDD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
	  if(!track->HasPointOnITSLayer(3)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(23)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(3))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(23)))->Fill(0);
	}//point on SDD2
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
	  if(!track->HasPointOnITSLayer(4)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(24)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(4))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(24)))->Fill(0);
	}//point on SSD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
	  if(!track->HasPointOnITSLayer(5)) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(25)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(5))
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(25)))->Fill(0);
	}//point on SSD2
	if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
	  if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	    ((TH1F *)(fQAPrimaryProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
	}
	  if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	    ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
	}//number of TPC points for the dE/dx
	
	/*if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	  ((TH2F *)(fQAPrimaryProtonsRejectedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	  ((TH2F *)(fQAPrimaryProtonsAcceptedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	  ((TH2F *)(fQAPrimaryProtonsRejectedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	((TH2F *)(fQAPrimaryProtonsAcceptedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));*/
      }//primary particle cut
      
      //Secondaries
      if(label > nPrimaries) {
	if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
	  if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  }
	  else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
	}//ITS clusters
	if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
	  if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  }
	  else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
	}//chi2 per ITS cluster
	if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
	  if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	    //cout<<"Secondary proton rejected"<<endl;
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  }
	  else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	    //cout<<"Secondary proton accepted"<<endl;
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	  }
	}//TPC clusters
	if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
	  if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  }
	  else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
	}//chi2 per TPC cluster
	if(fProtonAnalysisBase->IsUsedMaxCov11()) {
	  if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  }
	  else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov22()) {
	  if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  }
	  else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov33()) {
	  if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  }
	  else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov44()) {
	  if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  }
	  else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
	}//cov11
	if(fProtonAnalysisBase->IsUsedMaxCov55()) {
	  if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  }
	  else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
	}//cov55
	if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
	  if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	  }
	  else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}//sigma to vertex
	if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
	  if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	  }
	  else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}//sigma to vertex TPC
	if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
	  if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
	  }
	  else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
	}//DCA xy global tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
	  if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
	  }
	  else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
	}//DCA xy TPC tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
	  if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
	  }
	  else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
	}//DCA z global tracking
	if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
	  if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
	  }
	  else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
	}//DCA z TPC tracking
	if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
	  if(track->GetConstrainedChi2() > 0) {
	    if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	      ((TH1F *)(fQASecondaryProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	    }
	    else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	      ((TH1F *)(fQASecondaryProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	  }
	}//constrain chi2 - vertex
	if(fProtonAnalysisBase->IsUsedITSRefit()) {
	  if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(16)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(16)))->Fill(0);
	}//ITS refit
	if(fProtonAnalysisBase->IsUsedTPCRefit()) {
	  if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(17)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(17)))->Fill(0);
	}//TPC refit
	if(fProtonAnalysisBase->IsUsedESDpid()) {
	  if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(18)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(18)))->Fill(0);
	}//ESD pid
	if(fProtonAnalysisBase->IsUsedTPCpid()) {
	  if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(19)))->Fill(0);
	  }
	  else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(19)))->Fill(0);
	}//TPC pid
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
	  if(!track->HasPointOnITSLayer(0)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(20)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(0))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(20)))->Fill(0);
	}//point on SPD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
	  if(!track->HasPointOnITSLayer(1)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(21)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(1))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(21)))->Fill(0);
	}//point on SPD2
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
	  if(!track->HasPointOnITSLayer(2)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(22)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(2))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(22)))->Fill(0);
	}//point on SDD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
	  if(!track->HasPointOnITSLayer(3)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(23)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(3))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(23)))->Fill(0);
	}//point on SDD2
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
	  if(!track->HasPointOnITSLayer(4)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(24)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(4))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(24)))->Fill(0);
	}//point on SSD1
	if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
	  if(!track->HasPointOnITSLayer(5)) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(25)))->Fill(0);
	  }
	  else if(track->HasPointOnITSLayer(5))
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(25)))->Fill(0);
	}//point on SSD2
	if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
	  if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	    ((TH1F *)(fQASecondaryProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
	  }
	  if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	    ((TH1F *)(fQASecondaryProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
	}//number of TPC points for the dE/dx

	/*if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	  if(gProcess == 4)
	    ((TH2F *)(fQASecondaryProtonsRejectedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	  if(gProcess == 13)
	    ((TH2F *)(fQASecondaryProtonsRejectedList->At(29)))->Fill(gPt,TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY()) {
	  if(gProcess == 4)
	    ((TH2F *)(fQASecondaryProtonsAcceptedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	  if(gProcess == 13)
	    ((TH2F *)(fQASecondaryProtonsAcceptedList->At(29)))->Fill(gPt,TMath::Abs(dca[0]));
	}
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	  if(gProcess == 4)
	    ((TH2F *)(fQASecondaryProtonsRejectedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
	  if(gProcess == 13)
	    ((TH2F *)(fQASecondaryProtonsRejectedList->At(30)))->Fill(gPt,TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ()) {
	  if(gProcess == 4)
	    ((TH2F *)(fQASecondaryProtonsAcceptedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
	  if(gProcess == 13)
	    ((TH2F *)(fQASecondaryProtonsAcceptedList->At(30)))->Fill(gPt,TMath::Abs(dca[1]));
	    }*/
      }//secondary particle cut
    }//protons
    
    //antiprotons
    if(track->Charge() < 0) {
      //Primaries
      if(label <= nPrimaries) {
	if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
	if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	}
	else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
	if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	}
	else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
	if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	  //cout<<"Primary antiproton rejected"<<endl;
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	}
	else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
	if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	}
	else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fProtonAnalysisBase->IsUsedMaxCov11()) {
	if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	}
	else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov22()) {
	if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	}
	else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov33()) {
	if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	}
	else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov44()) {
	if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	}
	else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov55()) {
	if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	}
	else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
	if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}
	else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }//sigma to vertex
      if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
	if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}
	else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
	if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
      }//DCA xy global tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
	if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
      }//DCA xy TPC tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
      }//DCA z global tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
      }//DCA z TPC tracking
      if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	    ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	    ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fProtonAnalysisBase->IsUsedITSRefit()) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(16)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fProtonAnalysisBase->IsUsedTPCRefit()) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(17)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fProtonAnalysisBase->IsUsedESDpid()) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(18)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fProtonAnalysisBase->IsUsedTPCpid()) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(19)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
      if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
	if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
	}
	if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
      }//number of TPC points for the dE/dx

      /*if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	((TH2F *)(fQAPrimaryAntiProtonsRejectedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	((TH2F *)(fQAPrimaryAntiProtonsAcceptedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	((TH2F *)(fQAPrimaryAntiProtonsRejectedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
      ((TH2F *)(fQAPrimaryAntiProtonsAcceptedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));*/
    }//primary particle cut

    //Secondaries
    if(label > nPrimaries) {
      if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
	if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	}
	else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
	if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	}
	else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
	if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	}
	else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
	}
      }//TPC clusters
      if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
	if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	}
	else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fProtonAnalysisBase->IsUsedMaxCov11()) {
	if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	}
	else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov22()) {
	if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	}
	else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov33()) {
	if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	}
	else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov44()) {
	if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	}
	else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fProtonAnalysisBase->IsUsedMaxCov55()) {
	if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	}
	else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
	if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}
	else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }//sigma to vertex
      if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
	if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
	}
	else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
	if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
      }//DCA xy global tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
	if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
	}
	else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
      }//DCA xy TPC tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
      }//DCA z global tracking
      if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
	if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
	}
	else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
      }//DCA z TPC tracking
      if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
	if(track->GetConstrainedChi2() > 0) {
	  if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	    ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	  }
	  else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	    ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
      }//constrain chi2 - vertex
      if(fProtonAnalysisBase->IsUsedITSRefit()) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(16)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(16)))->Fill(0);
      }//ITS refit
      if(fProtonAnalysisBase->IsUsedTPCRefit()) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(17)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(17)))->Fill(0);
      }//TPC refit
      if(fProtonAnalysisBase->IsUsedESDpid()) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(18)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(18)))->Fill(0);
      }//ESD pid
      if(fProtonAnalysisBase->IsUsedTPCpid()) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(19)))->Fill(0);
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(19)))->Fill(0);
      }//TPC pid
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
	if(!track->HasPointOnITSLayer(0)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(20)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(0))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(20)))->Fill(0);
      }//point on SPD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
	if(!track->HasPointOnITSLayer(1)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(21)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(1))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(21)))->Fill(0);
      }//point on SPD2
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
	if(!track->HasPointOnITSLayer(2)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(22)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(2))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(22)))->Fill(0);
      }//point on SDD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
	if(!track->HasPointOnITSLayer(3)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(23)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(3))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(23)))->Fill(0);
      }//point on SDD2
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
	if(!track->HasPointOnITSLayer(4)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(24)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(4))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(24)))->Fill(0);
      }//point on SSD1
      if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
	if(!track->HasPointOnITSLayer(5)) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(25)))->Fill(0);
	}
	else if(track->HasPointOnITSLayer(5))
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(25)))->Fill(0);
      }//point on SSD2
      if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
	if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
	}
	if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
      }//number of TPC points for the dE/dx

      /*if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	if(gProcess == 4)
	  ((TH2F *)(fQASecondaryAntiProtonsRejectedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	if(gProcess == 13)
	  ((TH2F *)(fQASecondaryAntiProtonsRejectedList->At(29)))->Fill(gPt,TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY()) {
	if(gProcess == 4)
	  ((TH2F *)(fQASecondaryAntiProtonsAcceptedList->At(27)))->Fill(gPt,TMath::Abs(dca[0]));
	if(gProcess == 13)
	  ((TH2F *)(fQASecondaryAntiProtonsAcceptedList->At(29)))->Fill(gPt,TMath::Abs(dca[0]));
      }
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	if(gProcess == 4)
	  ((TH2F *)(fQASecondaryAntiProtonsRejectedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
	if(gProcess == 13)
	  ((TH2F *)(fQASecondaryAntiProtonsRejectedList->At(30)))->Fill(gPt,TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ()) {
	if(gProcess == 4)
	  ((TH2F *)(fQASecondaryAntiProtonsAcceptedList->At(28)))->Fill(gPt,TMath::Abs(dca[1]));
	if(gProcess == 13)
	  ((TH2F *)(fQASecondaryAntiProtonsAcceptedList->At(30)))->Fill(gPt,TMath::Abs(dca[1]));
	  }*/
    }//secondary particle cut
  }//antiprotons
    //}//if TParticle
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
void AliProtonQAAnalysis::SetQAYPtBins(Int_t nbinsY, 
				       Double_t minY, Double_t maxY,
				       Int_t nbinsPt, 
				       Double_t minPt, Double_t maxPt) {
  //Initializes the QA binning
  fNBinsY = nbinsY;
  fMinY = minY; fMaxY = maxY;
  fNBinsPt = nbinsPt;
  fMinPt = minPt; fMaxPt = maxPt;
  InitQA();
  InitCutLists();
  InitVertexQA();
  if(fRunMCAnalysis) InitMCAnalysis();
  if(fRunEfficiencyAnalysis) InitEfficiencyAnalysis();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::SetQAYPtBins(Int_t nbinsY, Double_t *gY,
				       Int_t nbinsPt, Double_t *gPt) {
  //Initializes the QA binning - asymmetric binning case
  fUseAsymmetricBinning = kTRUE;
  fNBinsY = nbinsY;
  for(Int_t i = 0; i < nbinsY; i++) fY[i] = gY[i];
  fMinY = gY[0]; fMaxY = gY[nbinsPt];
  fNBinsPt = nbinsPt;
  for(Int_t i = 0; i < nbinsPt; i++) fPt[i] = gPt[i];
  fMinPt = gPt[0]; fMaxPt = gPt[nbinsPt];
  InitQA();
  InitCutLists();
  InitVertexQA();
  if(fRunMCAnalysis) InitMCAnalysis();
  if(fRunEfficiencyAnalysis) InitEfficiencyAnalysis();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitEfficiencyAnalysis() {
  //Initialization of the efficiency list - reconstruction & PID & cut 
  //efficiency
  //Adding each monitored object in the list
  fEfficiencyList = new TList();

  //MC primary protons and antiprotons for the reconstruction efficiency
  TH2D *gHistMCYPtProtons = 0x0;
  if(fUseAsymmetricBinning)
    gHistMCYPtProtons = new TH2D("gHistMCYPtProtons",
				 ";;P_{T} [GeV/c]",
				 fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtProtons = new TH2D("gHistMCYPtProtons",
				 ";;P_{T} [GeV/c]",
				 fNBinsY,fMinY,fMaxY,
				 fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtons->GetXaxis()->SetTitle("y");
  gHistMCYPtProtons->SetStats(kTRUE);
  gHistMCYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtons);
  TH2D *gHistMCYPtAntiProtons = 0x0;
  if(fUseAsymmetricBinning)
    gHistMCYPtAntiProtons = new TH2D("gHistMCYPtAntiProtons",
				     ";;P_{T} [GeV/c]",
				     fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtAntiProtons = new TH2D("gHistMCYPtAntiProtons",
				     ";y;P_{T} [GeV/c]",
				     fNBinsY,fMinY,fMaxY,
				     fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtons->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtons->SetStats(kTRUE);
  gHistMCYPtAntiProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtons);

  //MC secondary protons and antiprotons that come from weak decay for the reconstruction efficiency
  TH2D *gHistMCYPtProtonsFromWeak = 0x0;
  if(fUseAsymmetricBinning)
    gHistMCYPtProtonsFromWeak = new TH2D("gHistMCYPtProtonsFromWeak",
					 ";;P_{T} [GeV/c]",
					 fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtProtonsFromWeak = new TH2D("gHistMCYPtProtonsFromWeak",
					 ";;P_{T} [GeV/c]",
					 fNBinsY,fMinY,fMaxY,
					 fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistMCYPtProtonsFromWeak->SetStats(kTRUE);
  gHistMCYPtProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtonsFromWeak);
  TH2D *gHistMCYPtAntiProtonsFromWeak = 0x0;
  if(fUseAsymmetricBinning)
    gHistMCYPtAntiProtonsFromWeak = new TH2D("gHistMCYPtAntiProtonsFromWeak",
					     ";;P_{T} [GeV/c]",
					     fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtAntiProtonsFromWeak = new TH2D("gHistMCYPtAntiProtonsFromWeak",
					     ";y;P_{T} [GeV/c]",
					     fNBinsY,fMinY,fMaxY,
					     fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtonsFromWeak->SetStats(kTRUE);
  gHistMCYPtAntiProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtonsFromWeak);

  //MC secondary protons and antiprotons that come from hadronic interactions for the reconstruction efficiency
  TH2D *gHistMCYPtProtonsFromHadronic = 0x0;
  if(fUseAsymmetricBinning)
    gHistMCYPtProtonsFromHadronic = new TH2D("gHistMCYPtProtonsFromHadronic",
					     ";;P_{T} [GeV/c]",
					     fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtProtonsFromHadronic = new TH2D("gHistMCYPtProtonsFromHadronic",
					     ";;P_{T} [GeV/c]",
					     fNBinsY,fMinY,fMaxY,
					     fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistMCYPtProtonsFromHadronic->SetStats(kTRUE);
  gHistMCYPtProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtProtonsFromHadronic);
  TH2D *gHistMCYPtAntiProtonsFromHadronic = 0x0;
  if(fUseAsymmetricBinning)
  gHistMCYPtAntiProtonsFromHadronic = new TH2D("gHistMCYPtAntiProtonsFromHadronic",
					       ";;P_{T} [GeV/c]",
					       fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistMCYPtAntiProtonsFromHadronic = new TH2D("gHistMCYPtAntiProtonsFromHadronic",
						 ";y;P_{T} [GeV/c]",
						 fNBinsY,fMinY,fMaxY,
						 fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistMCYPtAntiProtonsFromHadronic->SetStats(kTRUE);
  gHistMCYPtAntiProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistMCYPtAntiProtonsFromHadronic);
  
  //ESD primary protons and antiprotons for the reconstruction efficiency
  TH2D *gHistESDYPtProtons = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtProtons = new TH2D("gHistESDYPtProtons",
				  ";;P_{T} [GeV/c]",
				  fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtProtons = new TH2D("gHistESDYPtProtons",
				  ";;P_{T} [GeV/c]",
				  fNBinsY,fMinY,fMaxY,
				  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDYPtProtons->SetStats(kTRUE);
  gHistESDYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtons);
  TH2D *gHistESDYPtAntiProtons = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtAntiProtons = new TH2D("gHistESDYPtAntiProtons",
				      ";;P_{T} [GeV/c]",
				      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtAntiProtons = new TH2D("gHistESDYPtAntiProtons",
				      ";;P_{T} [GeV/c]",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtons->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtons->SetStats(kTRUE);
  gHistESDYPtAntiProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtons);

  //ESD (anti)protons from weak decays for the reconstruction efficiency
  TH2D *gHistESDYPtProtonsFromWeak = 0x0;
  if(fUseAsymmetricBinning)
  gHistESDYPtProtonsFromWeak = new TH2D("gHistESDYPtProtonsFromWeak",
					";;P_{T} [GeV/c]",
					fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtProtonsFromWeak = new TH2D("gHistESDYPtProtonsFromWeak",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsFromWeak->SetStats(kTRUE);
  gHistESDYPtProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtonsFromWeak);
  TH2D *gHistESDYPtAntiProtonsFromWeak = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtAntiProtonsFromWeak = new TH2D("gHistESDYPtAntiProtonsFromWeak",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtAntiProtonsFromWeak = new TH2D("gHistESDYPtAntiProtonsFromWeak",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsFromWeak->SetStats(kTRUE);
  gHistESDYPtAntiProtonsFromWeak->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtonsFromWeak);

  //ESD (anti)protons from hadronic interactions for the reconstruction efficiency
  TH2D *gHistESDYPtProtonsFromHadronic = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtProtonsFromHadronic = new TH2D("gHistESDYPtProtonsFromHadronic",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtProtonsFromHadronic = new TH2D("gHistESDYPtProtonsFromHadronic",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsFromHadronic->SetStats(kTRUE);
  gHistESDYPtProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtProtonsFromHadronic);
  TH2D *gHistESDYPtAntiProtonsFromHadronic = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtAntiProtonsFromHadronic = new TH2D("gHistESDYPtAntiProtonsFromHadronic",
						  ";;P_{T} [GeV/c]",
						  fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtAntiProtonsFromHadronic = new TH2D("gHistESDYPtAntiProtonsFromHadronic",
						  ";;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsFromHadronic->SetStats(kTRUE);
  gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDYPtAntiProtonsFromHadronic);
  
  //ESD reconstructed tracks that were initially protons for the PID efficiency
  TH3D *gHistESDInitYPtProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gNPoints[51];
    for(Int_t i = 0; i < 51; i++)
      gNPoints[i] = i*4; 
  gHistESDInitYPtProtons = new TH3D("gHistESDInitYPtProtons",
				    ";;P_{T} [GeV/c];N_{points}",
				    fNBinsY,fY,fNBinsPt,fPt,50,gNPoints);
  }
  else
    gHistESDInitYPtProtons = new TH3D("gHistESDInitYPtProtons",
				      ";;P_{T} [GeV/c];N_{points}",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt,
				      50,0,200);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDInitYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDInitYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDInitYPtProtons->SetStats(kTRUE);
  gHistESDInitYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDInitYPtProtons);
  
  //ESD reconstructed tracks that were initially protons and were identified as protons for the PID efficiency
  TH3D *gHistESDIdYPtProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gNPoints[51];
    for(Int_t i = 0; i < 51; i++)
      gNPoints[i] = i*4; 
    gHistESDIdYPtProtons = new TH3D("gHistESDIdYPtProtons",
				    ";;P_{T} [GeV/c];N_{points}",
				    fNBinsY,fY,fNBinsPt,fPt,50,gNPoints);
  }
  else
    gHistESDIdYPtProtons = new TH3D("gHistESDIdYPtProtons",
				    ";;P_{T} [GeV/c];N_{points}",
				    fNBinsY,fMinY,fMaxY,
				    fNBinsPt,fMinPt,fMaxPt,
				    50,0,200);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDIdYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDIdYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDIdYPtProtons->SetStats(kTRUE);
  gHistESDIdYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDIdYPtProtons);
 
  //ESD reconstructed tracks that were identified as protons for the PID contamination
  TH3D *gHistESDRecIdYPtProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gNPoints[51];
    for(Int_t i = 0; i < 51; i++)
      gNPoints[i] = i*4; 
    gHistESDRecIdYPtProtons = new TH3D("gHistESDRecIdYPtProtons",
				       ";;P_{T} [GeV/c];N_{points}",
				       fNBinsY,fY,fNBinsPt,fPt,50,gNPoints);
  }
  else
    gHistESDRecIdYPtProtons = new TH3D("gHistESDRecIdYPtProtons",
				       ";;P_{T} [GeV/c];N_{points}",
				       fNBinsY,fMinY,fMaxY,
				       fNBinsPt,fMinPt,fMaxPt,
				       50,0,200);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDRecIdYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDRecIdYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDRecIdYPtProtons->SetStats(kTRUE);
  gHistESDRecIdYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDRecIdYPtProtons);

  //ESD reconstructed tracks that were identified as protons but were initially not protons for the PID contamination
  TH3D *gHistESDContamYPtProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gNPoints[51];
    for(Int_t i = 0; i < 51; i++)
      gNPoints[i] = i*4; 
    gHistESDContamYPtProtons = new TH3D("gHistESDContamYPtProtons",
					";;P_{T} [GeV/c];N_{points}",
					fNBinsY,fY,fNBinsPt,fPt,50,gNPoints);
  }
  else
    gHistESDContamYPtProtons = new TH3D("gHistESDContamYPtProtons",
					";;P_{T} [GeV/c];N_{points}",
					fNBinsY,fMinY,fMaxY,
					fNBinsPt,fMinPt,fMaxPt,
					50,0,200);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDContamYPtProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDContamYPtProtons->GetXaxis()->SetTitle("y");
  gHistESDContamYPtProtons->SetStats(kTRUE);
  gHistESDContamYPtProtons->GetXaxis()->SetTitleColor(1);
  fEfficiencyList->Add(gHistESDContamYPtProtons);

  //==============//
  //Cut efficiency//
  //==============//
  fCutEfficiencyList = new TList();
  //Reconstructed primary tracks that were identified as protons
  TH2D *gHistESDYPtProtonsTotal = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtProtonsTotal = new TH2D("gHistESDYPtProtonsTotal",
				       ";;P_{T} [GeV/c]",
				       fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtProtonsTotal = new TH2D("gHistESDYPtProtonsTotal",
				       ";;P_{T} [GeV/c]",
				       fNBinsY,fMinY,fMaxY,
				       fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtProtonsTotal->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsTotal->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsTotal->SetStats(kTRUE);
  gHistESDYPtProtonsTotal->GetXaxis()->SetTitleColor(1);
  fCutEfficiencyList->Add(gHistESDYPtProtonsTotal);

  //Reconstructed primary tracks that were identified as antiprotons
  TH2D *gHistESDYPtAntiProtonsTotal = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtAntiProtonsTotal = new TH2D("gHistESDYPtAntiProtonsTotal",
					   ";;P_{T} [GeV/c]",
					   fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtAntiProtonsTotal = new TH2D("gHistESDYPtAntiProtonsTotal",
					   ";;P_{T} [GeV/c]",
					   fNBinsY,fMinY,fMaxY,
					   fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtAntiProtonsTotal->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsTotal->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsTotal->SetStats(kTRUE);
  gHistESDYPtAntiProtonsTotal->GetXaxis()->SetTitleColor(1);
  fCutEfficiencyList->Add(gHistESDYPtAntiProtonsTotal);
  
  //Reconstructed, survived primary tracks that were identified as protons
  TH2D *gHistESDYPtProtonsSurvived = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtProtonsSurvived = new TH2D("gHistESDYPtProtonsSurvived",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtProtonsSurvived = new TH2D("gHistESDYPtProtonsSurvived",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtProtonsSurvived->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtProtonsSurvived->GetXaxis()->SetTitle("y");
  gHistESDYPtProtonsSurvived->SetStats(kTRUE);
  gHistESDYPtProtonsSurvived->GetXaxis()->SetTitleColor(1);
  fCutEfficiencyList->Add(gHistESDYPtProtonsSurvived);
  
  //Reconstructed, survived primary tracks that were identified as antiprotons
  TH2D *gHistESDYPtAntiProtonsSurvived = 0x0;
  if(fUseAsymmetricBinning)
    gHistESDYPtAntiProtonsSurvived = new TH2D("gHistESDYPtAntiProtonsSurvived",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistESDYPtAntiProtonsSurvived = new TH2D("gHistESDYPtAntiProtonsSurvived",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistESDYPtAntiProtonsSurvived->GetXaxis()->SetTitle("#eta");
  else 
    gHistESDYPtAntiProtonsSurvived->GetXaxis()->SetTitle("y");
  gHistESDYPtAntiProtonsSurvived->SetStats(kTRUE);
  gHistESDYPtAntiProtonsSurvived->GetXaxis()->SetTitleColor(1);
  fCutEfficiencyList->Add(gHistESDYPtAntiProtonsSurvived);
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
  //eta-phi-Nclusters
  TH3D *gHistEtaPhiNClustersPrimaryProtonsPass = new TH3D("gHistEtaPhiNClustersPrimaryProtonsPass",
							  "Accepted primary protons;#eta;#phi;N_{clusters}(TPC)",
							  fNBinsY,fMinY,fMaxY,
							  100,0,360,
							  100,0,200);
  gHistEtaPhiNClustersPrimaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiNClustersPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiNClustersPrimaryProtonsPass);//eta-phi of primary accepted ESD protons
  TH3D *gHistEtaPhiNClustersPrimaryAntiProtonsPass = new TH3D("gHistEtaPhiNClustersPrimaryAntiProtonsPass",
							      "Accepted primary antiprotons;#eta;#phi;N_{clusters}(TPC)",
							      fNBinsY,fMinY,fMaxY,
							      100,0,360,
							      100,0,200);
  gHistEtaPhiNClustersPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiNClustersPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiNClustersPrimaryAntiProtonsPass);//eta-phi of primary accepted ESD antiprotons
  TH3D *gHistEtaPhiNClustersSecondaryProtonsPass = new TH3D("gHistEtaPhiNClustersSecondaryProtonsPass",
							    "Accepted secondary protons;#eta;#phi;N_{clusters}(TPC)",
							    fNBinsY,fMinY,fMaxY,
							    100,0,360,
							    100,0,200);
  gHistEtaPhiNClustersSecondaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiNClustersSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiNClustersSecondaryProtonsPass);//eta-phi of secondary accepted ESD protons
  TH3D *gHistEtaPhiNClustersSecondaryAntiProtonsPass = new TH3D("gHistEtaPhiNClustersSecondaryAntiProtonsPass",
								"Accepted secondary antiprotons;#eta;#phi;N_{clusters}(TPC)",
								fNBinsY,fMinY,fMaxY,
								100,0,360,
								100,0,200);
  gHistEtaPhiNClustersSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiNClustersSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiNClustersSecondaryAntiProtonsPass);//eta-phi of secondary accepted ESD antiprotons
  //eta-phi-chi^2 per TPC cluster
  TH3D *gHistEtaPhiChi2PerTPCClusterPrimaryProtonsPass = new TH3D("gHistEtaPhiChi2PerTPCClusterPrimaryProtonsPass",
								  "Accepted primary protons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								  fNBinsY,fMinY,fMaxY,
								  100,0,360,
								  100,0,4);
  gHistEtaPhiChi2PerTPCClusterPrimaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiChi2PerTPCClusterPrimaryProtonsPass);//eta-phi of primary accepted ESD protons
  TH3D *gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsPass = new TH3D("gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsPass",
								      "Accepted primary antiprotons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								      fNBinsY,fMinY,fMaxY,
								      100,0,360,
								      100,0,4);
  gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsPass);//eta-phi of primary accepted ESD antiprotons
  TH3D *gHistEtaPhiChi2PerTPCClusterSecondaryProtonsPass = new TH3D("gHistEtaPhiChi2PerTPCClusterSecondaryProtonsPass",
								    "Accepted secondary protons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								    fNBinsY,fMinY,fMaxY,
								    100,0,360,
								    100,0,4);
  gHistEtaPhiChi2PerTPCClusterSecondaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiChi2PerTPCClusterSecondaryProtonsPass);//eta-phi of secondary accepted ESD protons
  TH3D *gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsPass = new TH3D("gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsPass",
									"Accepted secondary antiprotons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
									fNBinsY,fMinY,fMaxY,
									100,0,360,
									100,0,4);
  gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsPass);//eta-phi of secondary accepted ESD antiprotons
  //eta-phi-number of TPC points for the dE/dx
  TH3D *gHistEtaPhiTPCdEdxNPointsPrimaryProtonsPass = new TH3D("gHistEtaPhiTPCdEdxNPointsPrimaryProtonsPass",
								  "Accepted primary protons;#eta;#phi;N_{points}(TPC)",
								  fNBinsY,fMinY,fMaxY,
								  100,0,360,
								  100,0,200);
  gHistEtaPhiTPCdEdxNPointsPrimaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiTPCdEdxNPointsPrimaryProtonsPass);//eta-phi of primary accepted ESD protons
  TH3D *gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsPass = new TH3D("gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsPass",
								      "Accepted primary antiprotons;#eta;#phi;N_{points}(TPC)",
								      fNBinsY,fMinY,fMaxY,
								      100,0,360,
								      100,0,200);
  gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsPass);//eta-phi of primary accepted ESD antiprotons
  TH3D *gHistEtaPhiTPCdEdxNPointsSecondaryProtonsPass = new TH3D("gHistEtaPhiTPCdEdxNPointsSecondaryProtonsPass",
								    "Accepted secondary protons;#eta;#phi;N_{points}(TPC)",
								    fNBinsY,fMinY,fMaxY,
								    100,0,360,
								    100,0,200);
  gHistEtaPhiTPCdEdxNPointsSecondaryProtonsPass->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiTPCdEdxNPointsSecondaryProtonsPass);//eta-phi of secondary accepted ESD protons
  TH3D *gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsPass = new TH3D("gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsPass",
									"Accepted secondary antiprotons;#eta;#phi;N_{points}(TPC)",
									fNBinsY,fMinY,fMaxY,
									100,0,360,
									100,0,200);
  gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fAcceptedCutList->Add(gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsPass);//eta-phi of secondary accepted ESD antiprotons

  TH1F *gPrimaryProtonsNPointsTPCdEdx = new TH1F("gPrimaryProtonsNPointsTPCdEdx",
					      ";N_{points} (TPC-dE/dx);Entries",
					      100,0,200);
  fAcceptedCutList->Add(gPrimaryProtonsNPointsTPCdEdx);
  TH1F *gPrimaryAntiProtonsNPointsTPCdEdx = new TH1F("gPrimaryAntiProtonsNPointsTPCdEdx",
						  ";N_{points} (TPC-dE/dx);Entries",
						  100,0,200);
  fAcceptedCutList->Add(gPrimaryAntiProtonsNPointsTPCdEdx);
  TH1F *gSecondaryProtonsNPointsTPCdEdx = new TH1F("gSecondaryProtonsNPointsTPCdEdx",
						";N_{points} (TPC-dE/dx);Entries",
						100,0,200);
  fAcceptedCutList->Add(gSecondaryProtonsNPointsTPCdEdx);
  TH1F *gSecondaryAntiProtonsNPointsTPCdEdx = new TH1F("gSecondaryAntiProtonsNPointsTPCdEdx",
						    ";N_{points} (TPC-dE/dx);Entries",
						    100,0,200);
  fAcceptedCutList->Add(gSecondaryAntiProtonsNPointsTPCdEdx);
  
  //Rejected cut list
  fRejectedCutList = new TList();
  //eta-phi-Nclusters
  TH3D *gHistEtaPhiNClustersPrimaryProtonsReject = new TH3D("gHistEtaPhiNClustersPrimaryProtonsReject",
							  "Rejected primary protons;#eta;#phi;N_{clusters}(TPC)",
							  fNBinsY,fMinY,fMaxY,
							  100,0,360,
							  100,0,200);
  gHistEtaPhiNClustersPrimaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiNClustersPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiNClustersPrimaryProtonsReject);//eta-phi of primary rejected ESD protons
  TH3D *gHistEtaPhiNClustersPrimaryAntiProtonsReject = new TH3D("gHistEtaPhiNClustersPrimaryAntiProtonsReject",
							      "Rejected primary antiprotons;#eta;#phi;N_{clusters}(TPC)",
							      fNBinsY,fMinY,fMaxY,
							      100,0,360,
							      100,0,200);
  gHistEtaPhiNClustersPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiNClustersPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiNClustersPrimaryAntiProtonsReject);//eta-phi of primary rejected ESD antiprotons
  TH3D *gHistEtaPhiNClustersSecondaryProtonsReject = new TH3D("gHistEtaPhiNClustersSecondaryProtonsReject",
							    "Rejected secondary protons;#eta;#phi;N_{clusters}(TPC)",
							    fNBinsY,fMinY,fMaxY,
							    100,0,360,
							    100,0,200);
  gHistEtaPhiNClustersSecondaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiNClustersSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiNClustersSecondaryProtonsReject);//eta-phi of secondary rejected ESD protons
  TH3D *gHistEtaPhiNClustersSecondaryAntiProtonsReject = new TH3D("gHistEtaPhiNClustersSecondaryAntiProtonsReject",
								"Rejected secondary antiprotons;#eta;#phi;N_{clusters}(TPC)",
								fNBinsY,fMinY,fMaxY,
								100,0,360,
								100,0,200);
  gHistEtaPhiNClustersSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiNClustersSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiNClustersSecondaryAntiProtonsReject);//eta-phi of secondary rejected ESD antiprotons
  //eta-phi-chi^2 per TPC cluster
  TH3D *gHistEtaPhiChi2PerTPCClusterPrimaryProtonsReject = new TH3D("gHistEtaPhiChi2PerTPCClusterPrimaryProtonsReject",
								  "Rejected primary protons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								  fNBinsY,fMinY,fMaxY,
								  100,0,360,
								  100,0,4);
  gHistEtaPhiChi2PerTPCClusterPrimaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiChi2PerTPCClusterPrimaryProtonsReject);//eta-phi of primary rejected ESD protons
  TH3D *gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsReject = new TH3D("gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsReject",
								      "Rejected primary antiprotons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								      fNBinsY,fMinY,fMaxY,
								      100,0,360,
								      100,0,4);
  gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiChi2PerTPCClusterPrimaryAntiProtonsReject);//eta-phi of primary rejected ESD antiprotons
  TH3D *gHistEtaPhiChi2PerTPCClusterSecondaryProtonsReject = new TH3D("gHistEtaPhiChi2PerTPCClusterSecondaryProtonsReject",
								    "Rejected secondary protons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
								    fNBinsY,fMinY,fMaxY,
								    100,0,360,
								    100,0,4);
  gHistEtaPhiChi2PerTPCClusterSecondaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiChi2PerTPCClusterSecondaryProtonsReject);//eta-phi of secondary rejected ESD protons
  TH3D *gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsReject = new TH3D("gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsReject",
									"Rejected secondary antiprotons;#eta;#phi;#chi^{2}/N_{clusters}(TPC)",
									fNBinsY,fMinY,fMaxY,
									100,0,360,
									100,0,4);
  gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiChi2PerTPCClusterSecondaryAntiProtonsReject);//eta-phi of secondary rejected ESD antiprotons
  //eta-phi-number of TPC points for the dE/dx
  TH3D *gHistEtaPhiTPCdEdxNPointsPrimaryProtonsReject = new TH3D("gHistEtaPhiTPCdEdxNPointsPrimaryProtonsReject",
								  "Rejected primary protons;#eta;#phi;N_{points}(TPC)",
								  fNBinsY,fMinY,fMaxY,
								  100,0,360,
								  100,0,200);
  gHistEtaPhiTPCdEdxNPointsPrimaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiTPCdEdxNPointsPrimaryProtonsReject);//eta-phi of primary rejected ESD protons
  TH3D *gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsReject = new TH3D("gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsReject",
								      "Rejected primary antiprotons;#eta;#phi;N_{points}(TPC)",
								      fNBinsY,fMinY,fMaxY,
								      100,0,360,
								      100,0,200);
  gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiTPCdEdxNPointsPrimaryAntiProtonsReject);//eta-phi of primary rejected ESD antiprotons
  TH3D *gHistEtaPhiTPCdEdxNPointsSecondaryProtonsReject = new TH3D("gHistEtaPhiTPCdEdxNPointsSecondaryProtonsReject",
								    "Rejected secondary protons;#eta;#phi;N_{points}(TPC)",
								    fNBinsY,fMinY,fMaxY,
								    100,0,360,
								    100,0,200);
  gHistEtaPhiTPCdEdxNPointsSecondaryProtonsReject->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiTPCdEdxNPointsSecondaryProtonsReject);//eta-phi of secondary rejected ESD protons
  TH3D *gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsReject = new TH3D("gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsReject",
									"Rejected secondary antiprotons;#eta;#phi;N_{points}(TPC)",
									fNBinsY,fMinY,fMaxY,
									100,0,360,
									100,0,200);
  gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fRejectedCutList->Add(gHistEtaPhiTPCdEdxNPointsSecondaryAntiProtonsReject);//eta-phi of secondary rejected ESD antiprotons

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

  //3D DCA vs pT plots
  TH3F *gHistPrimaryProtonsDCAxyEtaPt = new TH3F("gHistPrimaryProtonsDCAxyEtaPt",
						 ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						 9,-0.9,0.9,
						 6,0.45,1.05,
						 100,0,10);
  gHistPrimaryProtonsDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistPrimaryProtonsDCAxyEtaPt);
  TH3F *gHistPrimaryAntiProtonsDCAxyEtaPt = new TH3F("gHistPrimaryAntiProtonsDCAxyEtaPt",
						     ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						     9,-0.9,0.9,
						     6,0.45,1.05,
						     100,0,10);
  gHistPrimaryAntiProtonsDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistPrimaryAntiProtonsDCAxyEtaPt);
  TH3F *gHistSecondaryProtonsFromWeakDCAxyEtaPt = new TH3F("gHistSecondaryProtonsFromWeakDCAxyEtaPt",
							   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							   9,-0.9,0.9,
							   6,0.45,1.05,
							   100,0,10);
  gHistSecondaryProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryProtonsFromWeakDCAxyEtaPt);
  TH3F *gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt = new TH3F("gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       9,-0.9,0.9,
							       6,0.45,1.05,
							       100,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt);
  TH3F *gHistSecondaryProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryProtonsFromHadronicDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       9,-0.9,0.9,
							       6,0.45,1.05,
							       100,0,10);
  gHistSecondaryProtonsFromHadronicDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryProtonsFromHadronicDCAxyEtaPt);
  TH3F *gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
								   9,-0.9,0.9,
								   6,0.45,1.05,
								   100,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt);
  
  TH3F *gHistPrimaryProtonsDCAzEtaPt = new TH3F("gHistPrimaryProtonsDCAzEtaPt",
						";#eta;P_{T} [GeV/c];dca_{z} [cm]",
						9,-0.9,0.9,
						6,0.45,1.05,
						100,0,10);
  gHistPrimaryProtonsDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistPrimaryProtonsDCAzEtaPt);
  TH3F *gHistPrimaryAntiProtonsDCAzEtaPt = new TH3F("gHistPrimaryAntiProtonsDCAzEtaPt",
						    ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
						    9,-0.9,0.9,
						    6,0.45,1.05,
						    100,0,10);
  gHistPrimaryAntiProtonsDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistPrimaryAntiProtonsDCAzEtaPt);
  TH3F *gHistSecondaryProtonsFromWeakDCAzEtaPt = new TH3F("gHistSecondaryProtonsFromWeakDCAzEtaPt",
							  ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
							  9,-0.9,0.9,
							  6,0.45,1.05,
							  100,0,10);
  gHistSecondaryProtonsFromWeakDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryProtonsFromWeakDCAzEtaPt);
  TH3F *gHistSecondaryAntiProtonsFromWeakDCAzEtaPt = new TH3F("gHistSecondaryAntiProtonsFromWeakDCAzEtaPt",
							      ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
							      9,-0.9,0.9,
							      6,0.45,1.05,
							      100,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryAntiProtonsFromWeakDCAzEtaPt);
  TH3F *gHistSecondaryProtonsFromHadronicDCAzEtaPt = new TH3F("gHistSecondaryProtonsFromHadronicDCAzEtaPt",
							      ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
							      9,-0.9,0.9,
							      6,0.45,1.05,
							      100,0,10);
  gHistSecondaryProtonsFromHadronicDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryProtonsFromHadronicDCAzEtaPt);
  TH3F *gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt = new TH3F("gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt",
								  ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								  9,-0.9,0.9,
								  6,0.45,1.05,
								  100,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt->SetStats(kFALSE);
  fAcceptedDCAList->Add(gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitVertexQA() {
  //Initializes the Vertex QA histograms
  fQAVertexList = new TList();
  fQAVertexList->SetName("fQAVertexList");

  //Gen. multiplicity bins
  //Float_t xBins[24] = {0,1,2,4,6,8,10,15,20,30,40,50,75,100,
  //200,300,400,500,750,1000,1500,2000,2500,3000};
  //MC primary multiplicity (vertex efficiency calculation)
  TH1F *gHistMCPrimaryVz = new TH1F("gHistMCPrimaryVz",
				    ";V_{z} (gen.) [cm];Entries",
				    40,-20.,20.);
  fQAVertexList->Add(gHistMCPrimaryVz);
  //TH1F *gHistMCPrimaryMultiplicity = new TH1F("gHistMCPrimaryMultiplicity",
  //";N_{prim. gen.};Entries",
  //23,xBins);
  //fQAVertexList->Add(gHistMCPrimaryMultiplicity);
  
  //TPC
  TH1F *gHistTPCVz = new TH1F("gHistTPCVz",
			      ";V_{z} (gen.) [cm];Entries",
			      40,-20.,20.);
  fQAVertexList->Add(gHistTPCVz);
  //TH1F *gHistMCPrimaryMultiplicityTPC = new TH1F("gHistMCPrimaryMultiplicityTPC",
  //"Vertex TPC;N_{prim. gen.};Entries",
  //23,xBins);
  //fQAVertexList->Add(gHistMCPrimaryMultiplicityTPC);
  TH2F *gHistTPCESDVxN = new TH2F("gHistTPCESDVxN",
				 "Primary vertex TPC;V_{x} [cm];N_{contributors}",
				 100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistTPCESDVxN);
  TH2F *gHistTPCESDVyN = new TH2F("gHistTPCESDVyN",
				 "Primary vertex TPC;V_{y} [cm];N_{contributors}",
				 100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistTPCESDVyN);
  TH2F *gHistTPCESDVzN = new TH2F("gHistTPCESDVzN",
				 "Primary vertex TPC;V_{z} [cm];N_{contributors}",
				 100,-20.,20.,1000,0,5000);
  fQAVertexList->Add(gHistTPCESDVzN);
  TH1F *gHistTPCDiffVx = new TH1F("gHistTPCDiffVx",
				  ";V_{x}(rec.) - V_{x}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistTPCDiffVx);
  TH1F *gHistTPCDiffVy = new TH1F("gHistTPCDiffVy",
				  ";V_{y}(rec.) - V_{y}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistTPCDiffVy);
  TH1F *gHistTPCDiffVz = new TH1F("gHistTPCDiffVz",
				  ";V_{z}(rec.) - V_{z}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistTPCDiffVz);
  TH1F *gHistTPCResolutionVx = new TH1F("gHistTPCResolutionVx",
					";#sigma_{x} [#mu m];Entries",
					100,0.,1000000.);
  fQAVertexList->Add(gHistTPCResolutionVx);
  TH1F *gHistTPCResolutionVy = new TH1F("gHistTPCResolutionVy",
					";#sigma_{y} [#mu m];Entries",
					100,0.,1000000.);
  fQAVertexList->Add(gHistTPCResolutionVy);
  TH1F *gHistTPCResolutionVz = new TH1F("gHistTPCResolutionVz",
					";#sigma_{z} [#mu m];Entries",
					100,0.,6000.);
  fQAVertexList->Add(gHistTPCResolutionVz);
  
  //SPD
  TH1F *gHistSPDVz = new TH1F("gHistSPDVz",
			      ";V_{z} (gen.) [cm];Entries",
			      40,-20.,20.);
  fQAVertexList->Add(gHistSPDVz);
  //TH1F *gHistMCPrimaryMultiplicitySPD = new TH1F("gHistMCPrimaryMultiplicitySPD",
  //"Vertex SPD;N_{prim. gen.};Entries",
  //23,xBins);
  //fQAVertexList->Add(gHistMCPrimaryMultiplicitySPD);
  TH2F *gHistSPDESDVxN = new TH2F("gHistSPDESDVxN",
				 "Primary vertex SPD;V_{x} [cm];N_{contributors}",
				 100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistSPDESDVxN);
  TH2F *gHistSPDESDVyN = new TH2F("gHistSPDESDVyN",
				 "Primary vertex SPD;V_{y} [cm];N_{contributors}",
				 100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistSPDESDVyN);
  TH2F *gHistSPDESDVzN = new TH2F("gHistSPDESDVzN",
				 "Primary vertex SPD;V_{z} [cm];N_{contributors}",
				 100,-20.,20.,1000,0,5000);
  fQAVertexList->Add(gHistSPDESDVzN);
  TH1F *gHistSPDDiffVx = new TH1F("gHistSPDDiffVx",
				  ";V_{x}(rec.) - V_{x}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistSPDDiffVx);
  TH1F *gHistSPDDiffVy = new TH1F("gHistSPDDiffVy",
				  ";V_{y}(rec.) - V_{y}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistSPDDiffVy);
  TH1F *gHistSPDDiffVz = new TH1F("gHistSPDDiffVz",
				  ";V_{z}(rec.) - V_{z}(true) [#mu m];Entries",
				  100,-10000.,10000.);
  fQAVertexList->Add(gHistSPDDiffVz);
  TH1F *gHistSPDResolutionVx = new TH1F("gHistSPDResolutionVx",
					";#sigma_{x} [#mu m];Entries",
					100,0.,1000.);
  fQAVertexList->Add(gHistSPDResolutionVx);
  TH1F *gHistSPDResolutionVy = new TH1F("gHistSPDResolutionVy",
					";#sigma_{y} [#mu m];Entries",
					100,0.,1000.);
  fQAVertexList->Add(gHistSPDResolutionVy);
  TH1F *gHistSPDResolutionVz = new TH1F("gHistSPDResolutionVz",
					";#sigma_{z} [#mu m];Entries",
					100,0.,500.);
  fQAVertexList->Add(gHistSPDResolutionVz);
  
  //Tracks
  TH1F *gHistTracksVz = new TH1F("gHistTracksVz",
				 ";V_{z} (gen.) [cm];Entries",
				 40,-20.,20.);
  fQAVertexList->Add(gHistTracksVz);
  //TH1F *gHistMCPrimaryMultiplicityTracks = new TH1F("gHistMCPrimaryMultiplicityTracks",
  //"Vertex Tracks;N_{prim. gen.};Entries",
  //23,xBins);
  //fQAVertexList->Add(gHistMCPrimaryMultiplicityTracks);
  TH2F *gHistTracksESDVxN = new TH2F("gHistTracksESDVxN",
				     "Primary vertex Tracks;V_{x} [cm];N_{contributors}",
				     100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistTracksESDVxN);
  TH2F *gHistTracksESDVyN = new TH2F("gHistTracksESDVyN",
				    "Primary vertex Tracks;V_{y} [cm];N_{contributors}",
				    100,-10.,10.,1000,0,5000);
  fQAVertexList->Add(gHistTracksESDVyN);
  TH2F *gHistTracksESDVzN = new TH2F("gHistTracksESDVzN",
				    "Primary vertex Tracks;V_{z} [cm];N_{contributors}",
				    100,-20.,20.,1000,0,5000);
  fQAVertexList->Add(gHistTracksESDVzN);
  TH1F *gHistTracksDiffVx = new TH1F("gHistTracksDiffVx",
				     ";V_{x}(rec.) - V_{x}(true) [#mu m];Entries",
				     100,-10000.,10000.);
  fQAVertexList->Add(gHistTracksDiffVx);
  TH1F *gHistTracksDiffVy = new TH1F("gHistTracksDiffVy",
				     ";V_{y}(rec.) - V_{y}(true) [#mu m];Entries",
				     100,-10000.,10000.);
  fQAVertexList->Add(gHistTracksDiffVy);
  TH1F *gHistTracksDiffVz = new TH1F("gHistTracksDiffVz",
				     ";V_{z}(rec.) - V_{z}(true) [#mu m];Entries",
				     100,-10000.,10000.);
  fQAVertexList->Add(gHistTracksDiffVz);
  TH1F *gHistTracksResolutionVx = new TH1F("gHistTracksResolutionVx",
					   ";#sigma_{x} [#mu m];Entries",
					   100,0.,5000.);
  fQAVertexList->Add(gHistTracksResolutionVx);
  TH1F *gHistTracksResolutionVy = new TH1F("gHistTracksResolutionVy",
					   ";#sigma_{y} [#mu m];Entries",
					   100,0.,5000.);
  fQAVertexList->Add(gHistTracksResolutionVy);
  TH1F *gHistTracksResolutionVz = new TH1F("gHistTracksResolutionVz",
					   ";#sigma_{z} [#mu m];Entries",
					   100,0.,1000.);
  fQAVertexList->Add(gHistTracksResolutionVz);
}

//____________________________________________________________________//
void AliProtonQAAnalysis::InitQA() {
  //Initializes the QA histograms
  //if(!fQAHistograms) 
  SetRunQAAnalysis();

  //2D histograms
  //TDirectory *dir2D = gDirectory->mkdir("2D");
  //fGlobalQAList->Add(dir2D); dir2D->cd();
  TH2D *gHistYPtPrimaryProtonsPass = 0x0;
  if(fUseAsymmetricBinning)
  gHistYPtPrimaryProtonsPass = new TH2D("gHistYPtPrimaryProtonsPass",
					";;P_{T} [GeV/c]",
					fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryProtonsPass = new TH2D("gHistYPtPrimaryProtonsPass",
					  ";;P_{T} [GeV/c]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsPass);//y-pT of primary accepted ESD protons
  TH2D *gHistYPtPrimaryProtonsReject = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtPrimaryProtonsReject = new TH2D("gHistYPtPrimaryProtonsReject",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryProtonsReject = new TH2D("gHistYPtPrimaryProtonsReject",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fMinY,fMaxY,
					    fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryProtonsReject->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryProtonsReject->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsReject);//y-pT of primary rejected ESD protons

  TH2D *gHistYPtSecondaryProtonsPass = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtSecondaryProtonsPass = new TH2D("gHistYPtSecondaryProtonsPass",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtSecondaryProtonsPass = new TH2D("gHistYPtSecondaryProtonsPass",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fMinY,fMaxY,
					    fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitle("y");
  gHistYPtSecondaryProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsPass);//y-pT of secondary accepted ESD protons
  TH2D *gHistYPtSecondaryProtonsReject = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtSecondaryProtonsReject = new TH2D("gHistYPtSecondaryProtonsReject",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtSecondaryProtonsReject = new TH2D("gHistYPtSecondaryProtonsReject",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtSecondaryProtonsReject->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtSecondaryProtonsReject->GetXaxis()->SetTitle("y");
  gHistYPtSecondaryProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsReject);//y-pT of secondary rejected ESD protons

  TH2D *gHistYPtPrimaryAntiProtonsPass = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtPrimaryAntiProtonsPass = new TH2D("gHistYPtPrimaryAntiProtonsPass",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryAntiProtonsPass = new TH2D("gHistYPtPrimaryAntiProtonsPass",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsPass);//y-pT of primary accepted ESD antiprotons
  TH2D *gHistYPtPrimaryAntiProtonsReject = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtPrimaryAntiProtonsReject = new TH2D("gHistYPtPrimaryAntiProtonsReject",
						";;P_{T} [GeV/c]",
						fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryAntiProtonsReject = new TH2D("gHistYPtPrimaryAntiProtonsReject",
						";;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryAntiProtonsReject->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryAntiProtonsReject->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsReject);//y-pT of primary rejected ESD antiprotons

  TH2D *gHistYPtSecondaryAntiProtonsPass = 0x0;
  if(fUseAsymmetricBinning)
  gHistYPtSecondaryAntiProtonsPass = new TH2D("gHistYPtSecondaryAntiProtonsPass",
					      ";;P_{T} [GeV/c]",
					      fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtSecondaryAntiProtonsPass = new TH2D("gHistYPtSecondaryAntiProtonsPass",
						";;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtSecondaryAntiProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtSecondaryAntiProtonsPass->GetXaxis()->SetTitle("y");
  gHistYPtSecondaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsPass);//y-pT of secondary accepted ESD antiprotons
  TH2D *gHistYPtSecondaryAntiProtonsReject = 0x0;
  if(fUseAsymmetricBinning)
  gHistYPtSecondaryAntiProtonsReject = new TH2D("gHistYPtSecondaryAntiProtonsReject",
						";;P_{T} [GeV/c]",
						fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtSecondaryAntiProtonsReject = new TH2D("gHistYPtSecondaryAntiProtonsReject",
						  ";;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtSecondaryAntiProtonsReject->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtSecondaryAntiProtonsReject->GetXaxis()->SetTitle("y");
  gHistYPtSecondaryAntiProtonsReject->SetStats(kTRUE);
  gHistYPtSecondaryAntiProtonsReject->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiProtonsReject);//y-pT of secondary rejected ESD antiprotons

  TH2D *gHistYPtPrimaryProtonsMC = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtPrimaryProtonsMC = new TH2D("gHistYPtPrimaryProtonsMC",
					";;P_{T} [GeV/c]",
					fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryProtonsMC = new TH2D("gHistYPtPrimaryProtonsMC",
					";;P_{T} [GeV/c]",
					fNBinsY,fMinY,fMaxY,
					fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryProtonsMC->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryProtonsMC->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryProtonsMC->SetStats(kTRUE);
  gHistYPtPrimaryProtonsMC->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsMC);//y-pT of primary MC protons
  TH2D *gHistYPtPrimaryAntiProtonsMC = 0x0;
  if(fUseAsymmetricBinning)
    gHistYPtPrimaryAntiProtonsMC = new TH2D("gHistYPtPrimaryAntiProtonsMC",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fY,fNBinsPt,fPt);
  else
    gHistYPtPrimaryAntiProtonsMC = new TH2D("gHistYPtPrimaryAntiProtonsMC",
					    ";;P_{T} [GeV/c]",
					    fNBinsY,fMinY,fMaxY,
					    fNBinsPt,fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPrimaryAntiProtonsMC->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPrimaryAntiProtonsMC->GetXaxis()->SetTitle("y");
  gHistYPtPrimaryAntiProtonsMC->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsMC->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsMC);//y-pT of primary MC antiprotons

  TH3F *gHistYPtPDGProtonsPass = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gPDG[15] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5};
    gHistYPtPDGProtonsPass = new TH3F("gHistYPtPDGProtonsPass",
				      ";;P_{T} [GeV/c];PDG",
				      fNBinsY,fY,fNBinsPt,fPt,14,gPDG);
  }
  else
    gHistYPtPDGProtonsPass = new TH3F("gHistYPtPDGProtonsPass",
				      ";;P_{T} [GeV/c];PDG",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt,
				      14,-0.5,13.5);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPDGProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPDGProtonsPass->GetXaxis()->SetTitle("y");
  fQA2DList->Add(gHistYPtPDGProtonsPass);//composition of secondary protons
  TH3F *gHistYPtPDGAntiProtonsPass = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gPDG[15] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5};
    gHistYPtPDGAntiProtonsPass = new TH3F("gHistYPtPDGAntiProtonsPass",
					  ";;P_{T} [GeV/c];PDG",
					  fNBinsY,fY,fNBinsPt,fPt,14,gPDG);
  }
  else
    gHistYPtPDGAntiProtonsPass = new TH3F("gHistYPtPDGAntiProtonsPass",
					  ";;P_{T} [GeV/c];PDG",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt,
					  14,-0.5,13.5);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPDGAntiProtonsPass->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPDGAntiProtonsPass->GetXaxis()->SetTitle("y");
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
  TH1F *gPrimaryProtonsNumberOfTPCdEdxPointsPass = new TH1F("gPrimaryProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQAPrimaryProtonsAcceptedList->Add(gPrimaryProtonsNumberOfTPCdEdxPointsPass);

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
  TH1F *gPrimaryProtonsNumberOfTPCdEdxPointsReject = new TH1F("gPrimaryProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  fQAPrimaryProtonsRejectedList->Add(gPrimaryProtonsNumberOfTPCdEdxPointsReject);
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
  TH1F *gSecondaryProtonsNumberOfTPCdEdxPointsPass = new TH1F("gSecondaryProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQASecondaryProtonsAcceptedList->Add(gSecondaryProtonsNumberOfTPCdEdxPointsPass);

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
  TH1F *gSecondaryProtonsNumberOfTPCdEdxPointsReject = new TH1F("gSecondaryProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  fQASecondaryProtonsRejectedList->Add(gSecondaryProtonsNumberOfTPCdEdxPointsReject);  

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
  TH1F *gPrimaryAntiProtonsNumberOfTPCdEdxPointsPass = new TH1F("gPrimaryAntiProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQAPrimaryAntiProtonsAcceptedList->Add(gPrimaryAntiProtonsNumberOfTPCdEdxPointsPass);
  /*TH2F *gHistPrimaryAntiProtonsDCAxyPtPass = new TH2F("gHistPrimaryAntiProtonsDCAxyPtPass",
						      ";P_{T} [GeV/c];dca_{xy} [cm]",
						      16,0.3,1.1,
						      1000,0,10);
  gHistPrimaryAntiProtonsDCAxyPtPass->SetStats(kFALSE);
  fQAPrimaryAntiProtonsAcceptedList->Add(gHistPrimaryAntiProtonsDCAxyPtPass);
  TH2F *gHistPrimaryAntiProtonsDCAzPtPass = new TH2F("gHistPrimaryAntiProtonsDCAzPtPass",
						     ";P_{T} [GeV/c];dca_{z} [cm]",
						     16,0.3,1.1,
						     1000,0,10);
  gHistPrimaryAntiProtonsDCAzPtPass->SetStats(kFALSE);
  fQAPrimaryAntiProtonsAcceptedList->Add(gHistPrimaryAntiProtonsDCAzPtPass);*/

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
  TH1F *gPrimaryAntiProtonsNumberOfTPCdEdxPointsReject = new TH1F("gPrimaryAntiProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  fQAPrimaryAntiProtonsRejectedList->Add(gPrimaryAntiProtonsNumberOfTPCdEdxPointsReject);
  /*TH2F *gHistPrimaryAntiProtonsDCAxyPtReject = new TH2F("gHistPrimaryAntiProtonsDCAxyPtReject",
							";P_{T} [GeV/c];dca_{xy} [cm]",
							16,0.3,1.1,
							1000,0,10);
  gHistPrimaryAntiProtonsDCAxyPtReject->SetStats(kFALSE);
  fQAPrimaryAntiProtonsRejectedList->Add(gHistPrimaryAntiProtonsDCAxyPtReject);
  TH2F *gHistPrimaryAntiProtonsDCAzPtReject = new TH2F("gHistPrimaryAntiProtonsDCAzPtReject",
						       ";P_{T} [GeV/c];dca_{z} [cm]",
						       16,0.3,1.1,
						       1000,0,10);
  gHistPrimaryAntiProtonsDCAzPtReject->SetStats(kFALSE);
  fQAPrimaryAntiProtonsRejectedList->Add(gHistPrimaryAntiProtonsDCAzPtReject);*/
  
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
  TH1F *gSecondaryAntiProtonsNumberOfTPCdEdxPointsPass = new TH1F("gSecondaryAntiProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQASecondaryAntiProtonsAcceptedList->Add(gSecondaryAntiProtonsNumberOfTPCdEdxPointsPass);
  /*TH2F *gHistSecondaryAntiProtonsFromWeakDCAxyPtPass = new TH2F("gHistSecondaryAntiProtonsFromWeakDCAxyPtPass",
								";P_{T} [GeV/c];dca_{xy} [cm]",
								16,0.3,1.1,
								1000,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAxyPtPass->SetStats(kFALSE);
  fQASecondaryAntiProtonsAcceptedList->Add(gHistSecondaryAntiProtonsFromWeakDCAxyPtPass);
  TH2F *gHistSecondaryAntiProtonsFromWeakDCAzPtPass = new TH2F("gHistSecondaryAntiProtonsFromWeakDCAzPtPass",
							       ";P_{T} [GeV/c];dca_{z} [cm]",
							       16,0.3,1.1,
							       1000,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAzPtPass->SetStats(kFALSE);
  fQASecondaryAntiProtonsAcceptedList->Add(gHistSecondaryAntiProtonsFromWeakDCAzPtPass);
  TH2F *gHistSecondaryAntiProtonsFromHadronicDCAxyPtPass = new TH2F("gHistSecondaryAntiProtonsFromHadronicDCAxyPtPass",
								    ";P_{T} [GeV/c];dca_{xy} [cm]",
								    16,0.3,1.1,
								    1000,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAxyPtPass->SetStats(kFALSE);
  fQASecondaryAntiProtonsAcceptedList->Add(gHistSecondaryAntiProtonsFromHadronicDCAxyPtPass);
  TH2F *gHistSecondaryAntiProtonsFromHadronicDCAzPtPass = new TH2F("gHistSecondaryAntiProtonsFromHadronicDCAzPtPass",
								   ";P_{T} [GeV/c];dca_{z} [cm]",
								   16,0.3,1.1,
								   1000,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAzPtPass->SetStats(kFALSE);
  fQASecondaryAntiProtonsAcceptedList->Add(gHistSecondaryAntiProtonsFromHadronicDCAzPtPass);*/
  
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
  TH1F *gSecondaryAntiProtonsNumberOfTPCdEdxPointsReject = new TH1F("gSecondaryAntiProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  fQASecondaryAntiProtonsRejectedList->Add(gSecondaryAntiProtonsNumberOfTPCdEdxPointsReject);
  /*TH2F *gHistSecondaryAntiProtonsFromWeakDCAxyPtReject = new TH2F("gHistSecondaryAntiProtonsFromWeakDCAxyPtReject",
								  ";P_{T} [GeV/c];dca_{xy} [cm]",
								  16,0.3,1.1,
								  1000,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAxyPtReject->SetStats(kFALSE);
  fQASecondaryAntiProtonsRejectedList->Add(gHistSecondaryAntiProtonsFromWeakDCAxyPtReject);
  TH2F *gHistSecondaryAntiProtonsFromWeakDCAzPtReject = new TH2F("gHistSecondaryAntiProtonsFromWeakDCAzPtReject",
								 ";P_{T} [GeV/c];dca_{z} [cm]",
								 16,0.3,1.1,
								 1000,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAzPtReject->SetStats(kFALSE);
  fQASecondaryAntiProtonsRejectedList->Add(gHistSecondaryAntiProtonsFromWeakDCAzPtReject);
  TH2F *gHistSecondaryAntiProtonsFromHadronicDCAxyPtReject = new TH2F("gHistSecondaryAntiProtonsFromHadronicDCAxyPtReject",
								      ";P_{T} [GeV/c];dca_{xy} [cm]",
								      16,0.3,1.1,
								      1000,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAxyPtReject->SetStats(kFALSE);
  fQASecondaryAntiProtonsRejectedList->Add(gHistSecondaryAntiProtonsFromHadronicDCAxyPtReject);
  TH2F *gHistSecondaryAntiProtonsFromHadronicDCAzPtReject = new TH2F("gHistSecondaryAntiProtonsFromHadronicDCAzPtReject",
								     ";P_{T} [GeV/c];dca_{z} [cm]",
								     16,0.3,1.1,
								     1000,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAzPtReject->SetStats(kFALSE);
  fQASecondaryAntiProtonsRejectedList->Add(gHistSecondaryAntiProtonsFromHadronicDCAzPtReject);*/
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunReconstructionEfficiencyAnalysis(AliMCEvent *const mcEvent, 
							      AliESDEvent *esd,
							      const AliESDVertex *vertex) {
  //Run the reconstruction efficiency code (primaries & secondaries)
  AliStack *stack = mcEvent->Stack();

  Int_t nMCParticles = mcEvent->GetNumberOfTracks();
  Int_t nMCLabelCounter = 0;
  TArrayI labelMCArray(nMCParticles);

  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
    AliMCParticle *mcTrack = (AliMCParticle*) mcEvent->GetTrack(iTracks);
    if (!mcTrack) {
      Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
      continue;
    }

    Double_t vz = mcTrack->Zv();
    if (TMath::Abs(vz) > 50.) continue;//exclude particles generated out of the acceptance

    if(TMath::Abs(mcTrack->Eta()) > 1.0) continue;//acceptance
    if((mcTrack->Pt() > fMaxPt)||(mcTrack->Pt() < fMinPt)) continue;
    if(fProtonAnalysisBase->GetEtaMode()) {
      if((mcTrack->Eta() > fMaxY)|| (mcTrack->Eta() < fMinY)) continue;
    }
    else 
      if((fProtonAnalysisBase->Rapidity(mcTrack->Px(),mcTrack->Py(),mcTrack->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(mcTrack->Px(),mcTrack->Py(),mcTrack->Pz()) < fMinY)) continue;
    
    // Loop over Track References
    Bool_t labelTPC = kFALSE;
    AliTrackReference* trackRef = 0;
    for (Int_t iTrackRef = 0; iTrackRef  < mcTrack->GetNumberOfTrackReferences(); iTrackRef++) {
      trackRef = mcTrack->GetTrackReference(iTrackRef);
      if(trackRef) {
	Int_t detectorId = trackRef->DetectorId(); 
	if (detectorId == AliTrackReference::kTPC) {	    
	  labelTPC = kTRUE;
	  break;
	}
      }      
    }

    //findable tracks
    //if (labelTPC) {
      TParticle* particle = mcTrack->Particle();
      if(!particle) continue;
      Int_t pdgcode = particle->GetPdgCode();
      if(TMath::Abs(pdgcode) != 2212) continue;
      
      labelMCArray.AddAt(iTracks,nMCLabelCounter);
      nMCLabelCounter += 1;

      if(iTracks <= stack->GetNprimary()) {
	if(pdgcode == 2212) {
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(0)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										   particle->Py(),
										   particle->Pz()),
						     particle->Pt());
	}//protons
	if(pdgcode == -2212) {
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(1)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	AliMCParticle *mcMotherTrack = (AliMCParticle*) mcEvent->GetTrack(lPartMother);
	TParticle *motherParticle = mcMotherTrack->Particle();
	if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();
	
	if(pdgcode == 2212) {
	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fProtonAnalysisBase->GetEtaMode()) 
	      ((TH2D *)(fEfficiencyList->At(2)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(2)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode()) 
	      ((TH2D *)(fEfficiencyList->At(4)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(4)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//hadronic interactions
	}//protons
	if(pdgcode == -2212) {
	  if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	    if(fProtonAnalysisBase->GetEtaMode()) 
	      ((TH2D *)(fEfficiencyList->At(3)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(3)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode()) 
	      ((TH2D *)(fEfficiencyList->At(5)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(5)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								particle->Py(),
								particle->Pz()),
						       particle->Pt());
	  }//hadronic interactions
	}//antiprotons
      }//secondaries
      //}//findable tracks
  }//MC track loop

  //ESD track loop
  Bool_t iFound = kFALSE;
  Int_t mcGoods = nMCLabelCounter;
  for (Int_t k = 0; k < mcGoods; k++) {
    Int_t mcLabel = labelMCArray.At(k);
    iFound = kFALSE;

    Int_t nGoodTracks = esd->GetNumberOfTracks();
    TArrayI labelArray(nGoodTracks);
    Int_t labelCounter = 0;
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if(!track) continue;
            
      //TPC only
      if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
	AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
	if(!tpcTrack) continue;
	
	Int_t label = TMath::Abs(track->GetTPCLabel());
	if(IsLabelUsed(labelArray,label)) continue;
	labelArray.AddAt(label,labelCounter);
	labelCounter += 1;
	
	if (mcLabel != TMath::Abs(label)) continue;
	if(mcLabel != label) continue;
	if(label > stack->GetNtrack()) continue;

	TParticle *particle = stack->Particle(label);
	if(!particle) continue;
	Int_t pdgcode = particle->GetPdgCode();
	if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
	if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
	if(fProtonAnalysisBase->GetEtaMode()) {
	  if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
	}
	else 
	  if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
		
	if(fUseCutsInEfficiency) {
	  if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	  if(!fProtonAnalysisBase->IsAccepted(track)) continue;
	}	
	//reconstructed primary (anti)protons
	if(pdgcode == 2212) {
	  if(label <= stack->GetNprimary()) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
							 particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(8)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								  particle->Py(),
								  particle->Pz()),
							 particle->Pt());
	    }//weak decays
	    if((particle->GetUniqueID() == 13)) {
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							  particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								   particle->Py(),
								   particle->Pz()),
							  particle->Pt());
	    }//hadronic interactions
	  }//secondaries
	}//initial protons
	if(pdgcode == -2212) {	
	  if(label <= stack->GetNprimary()) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
							 particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(9)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								  particle->Py(),
								  particle->Pz()),
							 particle->Pt());
	    }//weak decays
	    if((particle->GetUniqueID() == 13)) {
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							  particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								   particle->Py(),
								   particle->Pz()),
							  particle->Pt());
	    }//hadronic interactions
	  }//secondaries
	}//initial antiprotons	
      }//TPC only tracks
      else {
	Int_t label = TMath::Abs(track->GetLabel());
	if(IsLabelUsed(labelArray,label)) continue;
	labelArray.AddAt(label,labelCounter);
	labelCounter += 1;
	
	if (mcLabel != TMath::Abs(label)) continue;
	if(mcLabel != label) continue;
	if(label > stack->GetNtrack()) continue;
		
	TParticle *particle = stack->Particle(label);
	if(!particle) continue;
	Int_t pdgcode = particle->GetPdgCode();
	if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
	if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
	if(fProtonAnalysisBase->GetEtaMode()) {
	  if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
	}
	else 
	  if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
	
	//Double_t probability[5];

	if(fUseCutsInEfficiency) {
	  if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	  if(!fProtonAnalysisBase->IsAccepted(track)) continue;
	}	
	//reconstructed primary (anti)protons
	if(pdgcode == 2212) {
	  if(label <= stack->GetNprimary()) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
							 particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(8)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								  particle->Py(),
								  particle->Pz()),
							 particle->Pt());
	    }//weak decays
	    if((particle->GetUniqueID() == 13)) {
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							  particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								   particle->Py(),
								   particle->Pz()),
							  particle->Pt());
	    }//hadronic interactions
	  }//secondaries
	}//initial protons
	if(pdgcode == -2212) {	
	  if(label <= stack->GetNprimary()) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
							 particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(9)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								  particle->Py(),
								  particle->Pz()),
							 particle->Pt());
	    }//weak decays
	    if((particle->GetUniqueID() == 13)) {
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							  particle->Pt());
	      else
		((TH2D *)(fEfficiencyList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
								   particle->Py(),
								   particle->Pz()),
							  particle->Pt());
	    }//hadronic interactions
	  }//secondaries
	}//initial antiprotons	
      }//global tracking
    }//track loop
    labelArray.Reset();
  }//loop over findable tracks

  labelMCArray.Reset();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunPIDEfficiencyAnalysis(AliStack *const stack, 
						   AliESDEvent *esd,
						   const AliESDVertex *vertex) {
  //Runs the PID efficiency analysis
  Int_t nGoodTracks = esd->GetNumberOfTracks();
  TArrayI labelArray(nGoodTracks);
  Int_t labelCounter = 0;
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    if(!track) continue;
    
    //TPC only
    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
    }
	
    Int_t label = TMath::Abs(track->GetLabel());
    if(IsLabelUsed(labelArray,label)) continue;
    labelArray.AddAt(label,labelCounter);
    labelCounter += 1;
    if(label > stack->GetNtrack()) continue;

    TParticle *particle = stack->Particle(label);
    if(!particle) continue;
    Int_t pdgcode = particle->GetPdgCode();
    
    Int_t nTPCpoints = track->GetTPCsignalN();

    if(fUseCutsInEfficiency) {
      if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
      if(!fProtonAnalysisBase->IsAccepted(track)) continue;
    }	
    if(TMath::Abs(pdgcode) == 2212) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH3D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						  particle->Pt(),nTPCpoints);
      else
	((TH3D *)(fEfficiencyList->At(12)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt(),nTPCpoints);
    }

    //pid
    if(fProtonAnalysisBase->IsProton(track)) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH3D *)(fEfficiencyList->At(14)))->Fill(particle->Eta(),
						  particle->Pt(),nTPCpoints);
      else ((TH3D *)(fEfficiencyList->At(14)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt(),nTPCpoints);
      if(TMath::Abs(pdgcode) == 2212) {
	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3D *)(fEfficiencyList->At(13)))->Fill(particle->Eta(),
						    particle->Pt(),nTPCpoints);
	else
	  ((TH3D *)(fEfficiencyList->At(13)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt(),nTPCpoints);
      }//properly identified as proton
      else {
	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3D *)(fEfficiencyList->At(15)))->Fill(particle->Eta(),
						    particle->Pt(),nTPCpoints);
	else
	  ((TH3D *)(fEfficiencyList->At(15)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt(),nTPCpoints);
      }//contamination
    }//identified as proton
  }//ESD track loop
  labelArray.Reset();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunCutEfficiencyAnalysis(AliStack *const stack, 
						   AliESDEvent *esd,
						   const AliESDVertex *vertex) {
  //Runs the cut efficiency analysis
  Int_t nGoodTracks = esd->GetNumberOfTracks();
  TArrayI labelArray(nGoodTracks);
  Int_t labelCounter = 0;
  Int_t nPrimaries = stack->GetNprimary();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    if(!track) continue;
    
    //TPC only
    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
    }
	
    Int_t label = TMath::Abs(track->GetLabel());
    if(IsLabelUsed(labelArray,label)) continue;
    labelArray.AddAt(label,labelCounter);
    labelCounter += 1;
    if(label > stack->GetNtrack()) continue;

    TParticle *particle = stack->Particle(label);
    if(!particle) continue;
		
    //select primaries
    if(label > nPrimaries) continue;
    //select identified protons
    if(!fProtonAnalysisBase->IsProton(track)) continue;
 
    if(track->Charge() > 0) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2F *)(fCutEfficiencyList->At(0)))->Fill(particle->Eta(),
						    particle->Pt());
      else
	((TH2F *)(fCutEfficiencyList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt());
    }

    if(track->Charge() < 0) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2F *)(fCutEfficiencyList->At(1)))->Fill(particle->Eta(),
						    particle->Pt());
      else
	((TH2F *)(fCutEfficiencyList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt());
    }

    //survived tracks
    if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
    if(!fProtonAnalysisBase->IsAccepted(track)) continue;
    if(track->Charge() > 0) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2F *)(fCutEfficiencyList->At(2)))->Fill(particle->Eta(),
						    particle->Pt());
      else
	((TH2F *)(fCutEfficiencyList->At(2)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt());
    }

    if(track->Charge() < 0) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2F *)(fCutEfficiencyList->At(3)))->Fill(particle->Eta(),
						    particle->Pt());
      else
	((TH2F *)(fCutEfficiencyList->At(3)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()),particle->Pt());
    }
  }//ESD track loop
  labelArray.Reset();
}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunEfficiencyAnalysis(AliStack *const stack, 
						AliESDEvent *esd,
						const AliESDVertex *vertex) {
  //Runs the efficiency code
  //MC loop
  Int_t nMCProtons = 0, nESDProtons = 0;
  for(Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if(fProtonAnalysisBase->GetEtaMode()) {
      if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
    }
    else 
      if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;

    Int_t pdgcode = particle->GetPdgCode();
    if(TMath::Abs(pdgcode) != 2212) continue;

    if(iParticle <= stack->GetNprimary()) {
      if(pdgcode == 2212) {
	nMCProtons += 1;
	if(fProtonAnalysisBase->GetEtaMode()) 
	  ((TH2D *)(fEfficiencyList->At(0)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
      }//protons
      if(pdgcode == -2212) {
	if(fProtonAnalysisBase->GetEtaMode()) 
	  ((TH2D *)(fEfficiencyList->At(1)))->Fill(particle->Eta(),
						   particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(2)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(2)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										   particle->Py(),
										   particle->Pz()),
						     particle->Pt());
	}//weak decays
	if((particle->GetUniqueID() == 13)) {
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(4)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(4)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										   particle->Py(),
										   particle->Pz()),
						     particle->Pt());
	}//hadronic interactions
      }//protons
      if(pdgcode == -2212) {
	if((particle->GetUniqueID() == 4)&&(TMath::Abs(motherPDGCode) == 3122)) {
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(3)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(3)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										   particle->Py(),
										   particle->Pz()),
						     particle->Pt());
	}//weak decays
	if((particle->GetUniqueID() == 13)) {
	  if(fProtonAnalysisBase->GetEtaMode()) 
	    ((TH2D *)(fEfficiencyList->At(5)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(5)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
							      particle->Py(),
							      particle->Pz()),
						     particle->Pt());
	}//hadronic interactions
      }//antiprotons
    }//secondaries
  
  }//MC loop

  //ESD track loop
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
    if(label > stack->GetNtrack()) continue;

    TParticle *particle = stack->Particle(label);
    if(!particle) continue;
    Int_t pdgcode = particle->GetPdgCode();
    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    
    Double_t gPt = 0.0, gP = 0.0;
    
    //TPC only
    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      gPt = tpcTrack->Pt();
      gP = tpcTrack->P();
      
      if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
      if(fProtonAnalysisBase->GetEtaMode()) {
	if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
      }
      else 
	if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
      
      if(fUseCutsInEfficiency) {
	if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	if(!fProtonAnalysisBase->IsAccepted(track)) continue;
      }      
      //reconstructed primary (anti)protons
      if(pdgcode == 2212) {
	if(label <= stack->GetNprimary()) {
	  nESDProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										     particle->Py(),
										     particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										      particle->Py(),
										      particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial protons
      if(pdgcode == -2212) {	
	if(label <= stack->GetNprimary()) {
	  if(fProtonAnalysisBase->GetEtaMode())
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										     particle->Py(),
										     particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										      particle->Py(),
										      particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial antiprotons
    }//TPC only tracks
    else {
      if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
      if(fProtonAnalysisBase->GetEtaMode()) {
	if((particle->Eta() > fMaxY)|| (particle->Eta() < fMinY)) continue;
      }
      else {
	if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
      }
      
      if(fUseCutsInEfficiency) {
	if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	if(!fProtonAnalysisBase->IsAccepted(track)) continue;
      }      
      //reconstructed primary (anti)protons
      if(pdgcode == 2212) {
	if(label <= stack->GetNprimary()) {
	  nESDProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(8)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										     particle->Py(),
										     particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										      particle->Py(),
										      particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial protons
      if(pdgcode == -2212) {	
	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(particle->Eta(),
						    particle->Pt());
	else
	  ((TH2D *)(fEfficiencyList->At(12)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										  particle->Py(),
										  particle->Pz()),
						    particle->Pt());
	if(label <= stack->GetNprimary()) {
	  if(fProtonAnalysisBase->GetEtaMode())
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(particle->Eta(),
						     particle->Pt());
	  else
	    ((TH2D *)(fEfficiencyList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(particle->Eta(),
						       particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(9)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										     particle->Py(),
										     particle->Pz()),
						       particle->Pt());
	  }//weak decays
	  if((particle->GetUniqueID() == 13)) {
	    if(fProtonAnalysisBase->GetEtaMode())
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(particle->Eta(),
							particle->Pt());
	    else
	      ((TH2D *)(fEfficiencyList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										      particle->Py(),
										      particle->Pz()),
							particle->Pt());
	  }//hadronic interactions
	}//secondaries
      }//initial antiprotons
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
void AliProtonQAAnalysis::RunVertexQA(AliGenEventHeader *header,
				      AliESDEvent *const esd) {
  //Runs the vertex QA
  //MC vertex
  TArrayF primaryVertex(3);
  header->PrimaryVertex(primaryVertex);

  //Int_t nPrimaries = stack->GetNprimary();
  ((TH1F *)(fQAVertexList->At(0)))->Fill(primaryVertex[2]);

  //TPC vertex
  const AliESDVertex *vertexTPC = esd->GetPrimaryVertexTPC();
  if(!vertexTPC) {
    Printf("ERROR: Could not retrieve the TPC vertex");
    return;
  }
  if(vertexTPC->GetNContributors() > 0) {
    ((TH1F *)(fQAVertexList->At(1)))->Fill(primaryVertex[2]);
    ((TH2F *)(fQAVertexList->At(2)))->Fill(vertexTPC->GetXv(),
					   vertexTPC->GetNContributors());
    ((TH2F *)(fQAVertexList->At(3)))->Fill(vertexTPC->GetYv(),
					   vertexTPC->GetNContributors());
    ((TH2F *)(fQAVertexList->At(4)))->Fill(vertexTPC->GetZv(),
					   vertexTPC->GetNContributors());
    ((TH1F *)(fQAVertexList->At(5)))->Fill((vertexTPC->GetXv()-primaryVertex[0])*10000.);
    ((TH1F *)(fQAVertexList->At(6)))->Fill((vertexTPC->GetYv()-primaryVertex[1])*10000.);
    ((TH1F *)(fQAVertexList->At(7)))->Fill((vertexTPC->GetZv()-primaryVertex[2])*10000.);
    ((TH1F *)(fQAVertexList->At(8)))->Fill(vertexTPC->GetXRes()*10000.);
    ((TH1F *)(fQAVertexList->At(9)))->Fill(vertexTPC->GetYRes()*10000.);
    ((TH1F *)(fQAVertexList->At(10)))->Fill(vertexTPC->GetZRes()*10000.);
  }//TPC vertex

  //SPD vertex
  const AliESDVertex *vertexSPD = esd->GetPrimaryVertexSPD();
  if(!vertexSPD) {
    Printf("ERROR: Could not retrieve the SPD vertex");
    return;
  }
  if(vertexSPD->GetNContributors() > 0) {
    ((TH1F *)(fQAVertexList->At(11)))->Fill(primaryVertex[2]);
    ((TH2F *)(fQAVertexList->At(12)))->Fill(vertexSPD->GetXv(),
					    vertexSPD->GetNContributors());
    ((TH2F *)(fQAVertexList->At(13)))->Fill(vertexSPD->GetYv(),
					    vertexSPD->GetNContributors());
    ((TH2F *)(fQAVertexList->At(14)))->Fill(vertexSPD->GetZv(),
					    vertexSPD->GetNContributors());
    ((TH1F *)(fQAVertexList->At(15)))->Fill((vertexSPD->GetXv()-primaryVertex[0])*10000.);
    ((TH1F *)(fQAVertexList->At(16)))->Fill((vertexSPD->GetYv()-primaryVertex[1])*10000.);
    ((TH1F *)(fQAVertexList->At(17)))->Fill((vertexSPD->GetZv()-primaryVertex[2])*10000.);
    ((TH1F *)(fQAVertexList->At(18)))->Fill(vertexSPD->GetXRes()*10000.);
    ((TH1F *)(fQAVertexList->At(19)))->Fill(vertexSPD->GetYRes()*10000.);
    ((TH1F *)(fQAVertexList->At(20)))->Fill(vertexSPD->GetZRes()*10000.);
  }//SPD vertex
  
  //Tracks vertex
  const AliESDVertex *vertexTracks = esd->GetPrimaryVertex();
  if(!vertexTracks) {
    Printf("ERROR: Could not retrieve the Tracks vertex");
    return;
  }
  if(vertexTracks->GetNContributors() > 0) {
    ((TH1F *)(fQAVertexList->At(21)))->Fill(primaryVertex[2]);
    ((TH2F *)(fQAVertexList->At(22)))->Fill(vertexTracks->GetXv(),
					    vertexTracks->GetNContributors());
    ((TH2F *)(fQAVertexList->At(23)))->Fill(vertexTracks->GetYv(),
					    vertexTracks->GetNContributors());
    ((TH2F *)(fQAVertexList->At(24)))->Fill(vertexTracks->GetZv(),
					    vertexTracks->GetNContributors());
    ((TH1F *)(fQAVertexList->At(25)))->Fill((vertexTracks->GetXv()-primaryVertex[0])*10000.);
    ((TH1F *)(fQAVertexList->At(26)))->Fill((vertexTracks->GetYv()-primaryVertex[1])*10000.);
    ((TH1F *)(fQAVertexList->At(27)))->Fill((vertexTracks->GetZv()-primaryVertex[2])*10000.);
    ((TH1F *)(fQAVertexList->At(28)))->Fill(vertexTracks->GetXRes()*10000.);
    ((TH1F *)(fQAVertexList->At(29)))->Fill(vertexTracks->GetYRes()*10000.);
    ((TH1F *)(fQAVertexList->At(30)))->Fill(vertexTracks->GetZRes()*10000.);
  }//Tracks vertex

}

//____________________________________________________________________//
void AliProtonQAAnalysis::RunQAAnalysis(AliStack *stack, 
					AliESDEvent *esd,
					const AliESDVertex *vertex) {
  //Runs the QA code
  //MC loop
  for(Int_t iParticle = 0; iParticle < stack->GetNprimary(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if(fProtonAnalysisBase->GetEtaMode()) {
      if((particle->Eta() > fMaxY)||(particle->Eta() < fMinY)) continue;
    }
    else {
      if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
    }
    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2D *)(fQA2DList->At(8)))->Fill(particle->Eta(),
					   particle->Pt());
      else
	((TH2D *)(fQA2DList->At(8)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
									 particle->Py(),
									 particle->Pz()),
					   particle->Pt());
    }
    if(pdgcode == -2212) {
      if(fProtonAnalysisBase->GetEtaMode())
	((TH2D *)(fQA2DList->At(9)))->Fill(particle->Eta(),
					   particle->Pt());
      else
	((TH2D *)(fQA2DList->At(9)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
									 particle->Py(),
									 particle->Pz()),
					   particle->Pt());
    }
  }//MC loop

  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
  //ESD track loop
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
    if(label > stack->GetNtrack()) continue;

    TParticle *particle = stack->Particle(label);
    if(!particle) continue;
    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    
    AliESDtrack trackTPC;
    
    //in case it's a TPC only track relate it to the proper vertex
    //if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)&&(!fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      //if((fUseTPCOnly)&&(!fUseHybridTPC)) {
      /*track->GetImpactParametersTPC(dca,cov);
      if (dca[0]==0 && dca[1]==0)  
	track->RelateToVertexTPC(((AliESDEvent*)esd)->GetPrimaryVertexTPC(),esd->GetMagneticField(),kVeryBig);
      if (!track->FillTPCOnlyTrack(trackTPC)) {
	continue;
      }
      track = &trackTPC ;*/
    //}
    
    Double_t gPt = 0.0, gP = 0.0;
    //Double_t probability[5];
    //Float_t dcaXY = 0.0, dcaZ = 0.0;
    Double_t nSigmaToVertex = fProtonAnalysisBase->GetSigmaToVertex(track);
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
    Int_t npointsTPCdEdx = track->GetTPCsignalN();

    //TPC only
    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      gPt = tpcTrack->Pt();
      gP = tpcTrack->P();
      /*if(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)
	track->GetImpactParameters(dcaXY,dcaZ);
	else track->GetImpactParametersTPC(dcaXY,dcaZ);*/
      tpcTrack->PropagateToDCA(vertex,
                               esd->GetMagneticField(),
                               100.,dca,cov);

      //pid
      if(fProtonAnalysisBase->IsProton(track)) {
	if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt

	FillQA(stack,esd,vertex,track);
	if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	  if(fProtonAnalysisBase->IsAccepted(track)) {
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
		((TH3D *)(fAcceptedCutList->At(44)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(48)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(52)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1D *)(fAcceptedCutList->At(56)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(0)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(4)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(8)))->Fill(nSigmaToVertex);
		((TH3F *)(fAcceptedDCAList->At(12)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		((TH3F *)(fAcceptedDCAList->At(18)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(0)))->Fill(tpcTrack->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										   tpcTrack->Py(),
										   tpcTrack->Pz()),
						     gPt);
	      }//accepted primary protons
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
		((TH3D *)(fAcceptedCutList->At(45)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(49)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(53)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1D *)(fAcceptedCutList->At(57)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(1)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(5)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(9)))->Fill(nSigmaToVertex);
		((TH3F *)(fAcceptedDCAList->At(13)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		((TH3F *)(fAcceptedDCAList->At(19)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(4)))->Fill(tpcTrack->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(4)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										   tpcTrack->Py(),
										   tpcTrack->Pz()),
						     gPt);
	      }//accepted primary antiprotons
	    }//accepted primary particles
	    else if(label > stack->GetNprimary()) {
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
		((TH3D *)(fAcceptedCutList->At(46)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(50)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(54)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1D *)(fAcceptedCutList->At(58)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(2)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(6)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(10)))->Fill(nSigmaToVertex);
		if(particle->GetUniqueID() == 4) {
		  ((TH3F *)(fAcceptedDCAList->At(14)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(20)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(particle->GetUniqueID() == 13) {
		  ((TH3F *)(fAcceptedDCAList->At(16)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(22)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(2)))->Fill(tpcTrack->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(2)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										   tpcTrack->Py(),
										   tpcTrack->Pz()),
						     gPt);
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH3F *)(fQA2DList->At(10)))->Fill(tpcTrack->Eta(),gPt,
						      ConvertPDGToInt(motherPDGCode));
		else
		  ((TH3F *)(fQA2DList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										    tpcTrack->Py(),
										    tpcTrack->Pz()),
						      gPt,
						      ConvertPDGToInt(motherPDGCode));
	      }//accepted secondary protons
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
		((TH3D *)(fAcceptedCutList->At(47)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(51)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(55)))->Fill(tpcTrack->Eta(),
							   tpcTrack->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1F *)(fAcceptedCutList->At(59)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(3)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(7)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(11)))->Fill(nSigmaToVertex);
		if(particle->GetUniqueID() == 4) {
		  ((TH3F *)(fAcceptedDCAList->At(15)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(21)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(particle->GetUniqueID() == 13) {
		  ((TH3F *)(fAcceptedDCAList->At(17)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(23)))->Fill(tpcTrack->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(6)))->Fill(tpcTrack->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										   tpcTrack->Py(),
										   tpcTrack->Pz()),
						     gPt);
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH3F *)(fQA2DList->At(11)))->Fill(tpcTrack->Eta(),gPt,
						      ConvertPDGToInt(motherPDGCode));
		else
		  ((TH3F *)(fQA2DList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										    tpcTrack->Py(),
										    tpcTrack->Pz()),
						      gPt,
						      ConvertPDGToInt(motherPDGCode));
	      }//accepted secondary antiprotons
	    }//accepted secondary particles
	  }//accepted - track cuts
	}//primary-like cut
	else {
	  if(label <= stack->GetNprimary()) {
	    if(track->Charge() > 0) {
	      ((TH3D *)(fRejectedCutList->At(0)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(4)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(8)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
	      
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(1)))->Fill(tpcTrack->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
										 tpcTrack->Py(),
										 tpcTrack->Pz()),
						   gPt);
	    }
	    else if(track->Charge() < 0) {
	      ((TH3D *)(fRejectedCutList->At(1)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(5)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(9)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(5)))->Fill(tpcTrack->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(5)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							    tpcTrack->Py(),
							    tpcTrack->Pz()),
						   gPt);
	    }
	  }//rejected primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0) {
	      ((TH3D *)(fRejectedCutList->At(2)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(6)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(10)))->Fill(tpcTrack->Eta(),
							 tpcTrack->Phi()*180./TMath::Pi(),
							 npointsTPCdEdx);
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(3)))->Fill(tpcTrack->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(3)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							    tpcTrack->Py(),
							    tpcTrack->Pz()),
						   gPt);
	    }
	    else if(track->Charge() < 0) {
	      ((TH3D *)(fRejectedCutList->At(3)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(7)))->Fill(tpcTrack->Eta(),
							tpcTrack->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(11)))->Fill(tpcTrack->Eta(),
							 tpcTrack->Phi()*180./TMath::Pi(),
							 npointsTPCdEdx);

	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(7)))->Fill(tpcTrack->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							    tpcTrack->Py(),
							    tpcTrack->Pz()),
						   gPt);
	    }
	  }//rejected secondary particles
	}//rejected - track cuts
      }//proton check
    }//TPC only tracks
    //combined tracking
    else {
      gPt = track->Pt();
      gP = track->P();
      track->PropagateToDCA(vertex,
                               esd->GetMagneticField(),
                               100.,dca,cov);
      
      //pid
      if(fProtonAnalysisBase->IsProton(track)) {
	if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt

	FillQA(stack,esd,vertex,track);
	if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	  if(fProtonAnalysisBase->IsAccepted(track)) {
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
		((TH3D *)(fAcceptedCutList->At(44)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(48)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(52)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1F *)(fAcceptedCutList->At(56)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(0)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(4)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(8)))->Fill(nSigmaToVertex);
		((TH3F *)(fAcceptedDCAList->At(12)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		((TH3F *)(fAcceptedDCAList->At(18)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(0)))->Fill(track->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										   track->Py(),
										   track->Pz()),
						     gPt);
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
		((TH3D *)(fAcceptedCutList->At(45)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(49)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(53)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1F *)(fAcceptedCutList->At(57)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(1)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(5)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(9)))->Fill(nSigmaToVertex);
		((TH3F *)(fAcceptedDCAList->At(13)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		((TH3F *)(fAcceptedDCAList->At(19)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(4)))->Fill(track->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(4)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										   track->Py(),
										   track->Pz()),
						     gPt);
	      }
	    }//primary particles
	    else if(label > stack->GetNprimary()) {
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
		((TH3D *)(fAcceptedCutList->At(46)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(50)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(54)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1F *)(fAcceptedCutList->At(58)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(2)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(6)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(10)))->Fill(nSigmaToVertex);
		if(particle->GetUniqueID() == 4) {
		  ((TH3F *)(fAcceptedDCAList->At(14)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(20)))->Fill(track->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(particle->GetUniqueID() == 13) {
		  ((TH3F *)(fAcceptedDCAList->At(16)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(22)))->Fill(track->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(2)))->Fill(track->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(2)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										   track->Py(),
										   track->Pz()),
						     gPt);
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH3F *)(fQA2DList->At(10)))->Fill(track->Eta(),gPt,
						      ConvertPDGToInt(motherPDGCode));
		else
		  ((TH3F *)(fQA2DList->At(10)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										    track->Py(),
										    track->Pz()),
						      gPt,
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
		((TH3D *)(fAcceptedCutList->At(47)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   nClustersTPC);
		((TH3D *)(fAcceptedCutList->At(51)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   chi2PerClusterTPC);
		((TH3D *)(fAcceptedCutList->At(55)))->Fill(track->Eta(),
							   track->Phi()*180./TMath::Pi(),
							   npointsTPCdEdx);
		((TH1F *)(fAcceptedCutList->At(59)))->Fill(npointsTPCdEdx);
		
		((TH1F *)(fAcceptedDCAList->At(3)))->Fill(TMath::Abs(dca[0]));
		((TH1F *)(fAcceptedDCAList->At(7)))->Fill(TMath::Abs(dca[1]));
		((TH1F *)(fAcceptedDCAList->At(11)))->Fill(nSigmaToVertex);
		if(particle->GetUniqueID() == 4) {
		  ((TH3F *)(fAcceptedDCAList->At(15)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(21)))->Fill(track->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(particle->GetUniqueID() == 13) {
		  ((TH3F *)(fAcceptedDCAList->At(17)))->Fill(track->Eta(),gPt,TMath::Abs(dca[0]));
		  ((TH3F *)(fAcceptedDCAList->At(23)))->Fill(track->Eta(),gPt,TMath::Abs(dca[1]));
		}
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH2D *)(fQA2DList->At(6)))->Fill(track->Eta(),gPt);
		else
		  ((TH2D *)(fQA2DList->At(6)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										   track->Py(),
										   track->Pz()),
						     gPt);
		if(fProtonAnalysisBase->GetEtaMode())
		  ((TH3F *)(fQA2DList->At(11)))->Fill(track->Eta(),gPt,
						      ConvertPDGToInt(motherPDGCode));
		else
		  ((TH3F *)(fQA2DList->At(11)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
										    track->Py(),
										    track->Pz()),
						      gPt,
						      ConvertPDGToInt(motherPDGCode));
	      }
	    }//secondary particles
	  }//accepted - track cuts
	}//primary-like cut
	else if((!fProtonAnalysisBase->IsAccepted(track)) || 
		(!fProtonAnalysisBase->IsPrimary(esd,vertex,track))) {
	  if(label <= stack->GetNprimary()) {
	    if(track->Charge() > 0) {
	      ((TH3D *)(fRejectedCutList->At(0)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(4)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(8)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(1)))->Fill(track->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
							    track->Py(),
							    track->Pz()),
						   gPt);
	    }
	    else if(track->Charge() < 0) {
	      ((TH3D *)(fRejectedCutList->At(1)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(5)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(9)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
						
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(5)))->Fill(track->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(5)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
							    track->Py(),
							    track->Pz()),
						   gPt);
	    }
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0) {
	      ((TH3D *)(fRejectedCutList->At(2)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(6)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(10)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(3)))->Fill(track->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(3)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
							    track->Py(),
							    track->Pz()),
						   gPt);
	    }
	    else if(track->Charge() < 0) {
	      ((TH3D *)(fRejectedCutList->At(3)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							nClustersTPC);
	      ((TH3D *)(fRejectedCutList->At(7)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							chi2PerClusterTPC);
	      ((TH3D *)(fRejectedCutList->At(11)))->Fill(track->Eta(),
							track->Phi()*180./TMath::Pi(),
							npointsTPCdEdx);
	      if(fProtonAnalysisBase->GetEtaMode())
		((TH2D *)(fQA2DList->At(7)))->Fill(track->Eta(),gPt);
	      else
		((TH2D *)(fQA2DList->At(7)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
							    track->Py(),
							    track->Pz()),
						   gPt);
	    }
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
  TH3F *gHistYPtPDGProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gPDG[15] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5};
  gHistYPtPDGProtons = new TH3F("gHistYPtPDGProtons",
				";;P_{T} [GeV/c];PDG",
				fNBinsY,fY,fNBinsPt,fPt,14,gPDG);
  }
  else
    gHistYPtPDGProtons = new TH3F("gHistYPtPDGProtons",
				  ";;P_{T} [GeV/c];PDG",
				  fNBinsY,fMinY,fMaxY,
				  fNBinsPt,fMinPt,fMaxPt,
				  14,-0.5,13.5);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPDGProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPDGProtons->GetXaxis()->SetTitle("y");
  fPDGList->Add(gHistYPtPDGProtons);
  TH3F *gHistYPtPDGAntiProtons = 0x0;
  if(fUseAsymmetricBinning) {
    Double_t gPDG[15] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5};
    gHistYPtPDGAntiProtons = new TH3F("gHistYPtPDGAntiProtons",
				      ";;P_{T} [GeV/c];PDG",
				      fNBinsY,fY,fNBinsPt,fPt,14,gPDG);
  }
  else
    gHistYPtPDGAntiProtons = new TH3F("gHistYPtPDGAntiProtons",
				      ";;P_{T} [GeV/c];PDG",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt,
				      14,-0.5,13.5);
  if(fProtonAnalysisBase->GetEtaMode()) 
    gHistYPtPDGAntiProtons->GetXaxis()->SetTitle("#eta");
  else 
    gHistYPtPDGAntiProtons->GetXaxis()->SetTitle("y");
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
void AliProtonQAAnalysis::RunMCAnalysis(AliStack* const stack) {
  //Main analysis part - MC 
  for(Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++) {
    TParticle *particle = stack->Particle(iParticle);
    if(!particle) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if(fProtonAnalysisBase->GetEtaMode()) {
      if((particle->Eta() > fMaxY)||(particle->Eta() < fMinY)) continue;
    }
    else {
      if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;
    }

    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) {
      if(iParticle <= stack->GetNprimary()) {
	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3F *)(fPDGList->At(0)))->Fill(particle->Eta(),particle->Pt(),0);
	else
	  ((TH3F *)(fPDGList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
						     particle->Py(),
						     particle->Pz()),
					    particle->Pt(),0);
      }
      else if(iParticle > stack->GetNprimary()) {
	Int_t lPartMother = particle->GetFirstMother();
	TParticle *motherParticle = stack->Particle(lPartMother);
	if(!motherParticle) continue;
	Int_t motherPDGCode = motherParticle->GetPdgCode();
	if(fMCProcessIdFlag)
	  if(particle->GetUniqueID() != fMCProcessId) continue;
	if(fMotherParticlePDGCodeFlag)
	  if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3F *)(fPDGList->At(0)))->Fill(particle->Eta(),
					    particle->Pt(),
					    ConvertPDGToInt(motherParticle->GetPdgCode()));
	else
	  ((TH3F *)(fPDGList->At(0)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
      if(iParticle <= stack->GetNprimary()) {
	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3F *)(fPDGList->At(1)))->Fill(particle->Eta(),particle->Pt(),0);
	else
	  ((TH3F *)(fPDGList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
						     particle->Py(),
						     particle->Pz()),
					    particle->Pt(),0);
      }
      else if(iParticle > stack->GetNprimary()) {
	Int_t lPartMother = particle->GetFirstMother();
	TParticle *motherParticle = stack->Particle(lPartMother);
	if(!motherParticle) continue;
	Int_t motherPDGCode = motherParticle->GetPdgCode();
	if(fMCProcessIdFlag)
	  if(particle->GetUniqueID() != fMCProcessId) continue;
	if(fMotherParticlePDGCodeFlag)
	  if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;

	if(fProtonAnalysisBase->GetEtaMode())
	  ((TH3F *)(fPDGList->At(1)))->Fill(particle->Eta(),
					    particle->Pt(),
					    ConvertPDGToInt(motherParticle->GetPdgCode()));
	else
	  ((TH3F *)(fPDGList->At(1)))->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
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
Int_t AliProtonQAAnalysis::ConvertPDGToInt(Int_t pdgCode) const {
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






