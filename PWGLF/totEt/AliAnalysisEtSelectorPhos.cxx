#include "AliAnalysisEtSelectorPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDCaloCluster.h"
#include <AliVEvent.h>
#include "TRefArray.h"
#include "AliPHOSGeometry.h"
#include "TH2I.h"
#include "TFile.h"
#include "TMath.h"
#include "AliLog.h"
#include <iostream>
ClassImp(AliAnalysisEtSelectorPhos)
AliAnalysisEtSelectorPhos::AliAnalysisEtSelectorPhos(AliAnalysisEtCuts* cuts): AliAnalysisEtSelector(cuts)
,fGeoUtils(0)
,fBadMapM2(0)
,fBadMapM3(0)
,fBadMapM4(0)
{
  
}

AliAnalysisEtSelectorPhos::~AliAnalysisEtSelectorPhos()
{

}

TRefArray* AliAnalysisEtSelectorPhos::GetClusters()
{  
  return 0;
}

Int_t AliAnalysisEtSelectorPhos::Init(Int_t runNumber)
{
  AliAnalysisEtSelector::Init(runNumber);
  
  int res = LoadGeometry();
  if(res) return -1;
  return LoadBadMaps();  
}

Bool_t AliAnalysisEtSelectorPhos::CutMinEnergy(const AliESDCaloCluster& cluster) const
{
  return cluster.E() > fCuts->GetReconstructedPhosClusterEnergyCut();
}

Bool_t AliAnalysisEtSelectorPhos::CutDistanceToBadChannel(const AliESDCaloCluster& cluster) const
{
    Float_t gPos[3];
    cluster.GetPosition(gPos);
    Int_t relId[4];
    TVector3 glVec(gPos);
    fGeoUtils->GlobalPos2RelId(glVec, relId);

    TVector3 locVec;
    fGeoUtils->Global2Local(locVec, glVec, relId[0]);
//    std::cout << fGeoUtils << std::endl;
    //std::cout << relId[0] << " " << cluster.IsPHOS() <<  std::endl;
    //std::cout << locVec[0] << " " << " " << locVec[1] << " " << locVec[2] << std::endl;
    for (Int_t x = 0; x < fBadMapM2->GetNbinsX(); x++)
    {
        for (Int_t z = 0; z < fBadMapM2->GetNbinsY(); z++)
        {
            if (relId[0] == 3)
            {
                if (fBadMapM2->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 3;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));
                    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
		    
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		      std::cout << "Module 2, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kFALSE;
                    }
                }
            }
            if (relId[0] == 2)
            {
                if (fBadMapM3->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 2;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
		    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		      std::cout << "Module 3, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kFALSE;
                    }
                }
            }
            if (relId[0] == 1)
            {
                if (fBadMapM4->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 1;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
		    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//			std::cout << "Module 4, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kFALSE;
                    }
                }
            }

        }
    }

    return kTRUE;

}

Bool_t AliAnalysisEtSelectorPhos::CutTrackMatching(const AliESDCaloCluster& cluster, Double_t &r) const
{

  
  // cluster->GetTrackDx(), cluster->GetTrackDz(), event->GetTrack(trackMatchedIndex)->Pt(), event->GetTrack(trackMatchedIndex)->Charge(), ev
  
  Int_t nTracksMatched = cluster.GetNTracksMatched();
  if(nTracksMatched == 0) return kFALSE;
  
  Int_t trackMatchedIndex = cluster.GetTrackMatchedIndex();
  if(trackMatchedIndex < 0) return kTRUE;
  
  AliVParticle *track = fEvent->GetTrack(trackMatchedIndex);
  Double_t dx = cluster.GetTrackDx();
  Double_t dz = cluster.GetTrackDz();
  Double_t pt = track->Pt();
  Int_t charge = track->Charge();
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  
  Double_t mf = fEvent->GetMagneticField(); //Positive for ++ and negative for --

  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ= -0.468318;
    if(charge>0)
      meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  r = TMath::Sqrt(rx*rx+rz*rz);
  if(r < fCuts->GetPhosTrackRCut()) return kFALSE;
  
  return kTRUE;
 
}

int AliAnalysisEtSelectorPhos::LoadGeometry()
{

  fGeoUtils = AliPHOSGeometry::GetInstance("IHEP");
    // ifstream f("badchannels.txt", ios::in);
  return 0;
}

int AliAnalysisEtSelectorPhos::LoadBadMaps()
{
TFile *f = TFile::Open("badchannels.root", "READ");

    if(!f)
    {
      std::cout << "Could not open badchannels.root" << std::endl;
      return -1;
    }
    
    fBadMapM2 = (TH2I*)f->Get("bad_channels_m2");
    if(fBadMapM2) 
    {
      std::cout << "Could not find bad_channels_m2 in badchannels.root" << std::endl;
    }
    fBadMapM3 = (TH2I*)f->Get("bad_channels_m3");
    if(fBadMapM2) 
    {
      std::cout << "Could not find bad_channels_m3 in badchannels.root" << std::endl;
    }
    
    fBadMapM4 = (TH2I*)f->Get("bad_channels_m4");
    if(fBadMapM4) 
    {
      std::cout << "Could not find bad_channels_m4 in badchannels.root" << std::endl;
    }

    return 0;
    
    
    
}
