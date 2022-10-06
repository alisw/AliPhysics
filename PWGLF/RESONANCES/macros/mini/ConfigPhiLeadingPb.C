Bool_t ConfigPhiLeadingPb(
 AliRsnMiniAnalysisTask *task,
 Bool_t isMC = kFALSE, 
 Bool_t isPP = kFALSE, 
 Double_t nSigmaPart1TPC = -1, 
 Double_t nSigmaPart2TPC = -1,
 Double_t nSigmaPart1TOF = -1, 
 Double_t nSigmaPart2TOF = -1,
 Int_t                   customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
 Int_t                    aodFilterBit=0
 )
{
    // -- Values ------------------------------------------------------------------------------------
    /* invariant mass   */ Int_t imID = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
    /* transv. momentum */ Int_t ptID = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* angel to leading */ Int_t alID = task->CreateValue(AliRsnMiniValue::kAngleLeading, kFALSE);
    /* pt of leading    */ Int_t ptlID = task->CreateValue(AliRsnMiniValue::kLeadingPt, kFALSE);
    /* multiplicity     */ Int_t multID = task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
    /* delta eta        */ Int_t detaID  = task->CreateValue(AliRsnMiniValue::kDeltaEta, kFALSE);
   /* phi of leading    */ Int_t philID = task->CreateValue(AliRsnMiniValue::kLeadingPhi, kFALSE);
   /* phi angle         */ Int_t phiID = task->CreateValue(AliRsnMiniValue::kPhi, kFALSE);

  // Cuts
   AliRsnCutSetDaughterParticle* cutSetK;
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s;

  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  cout<<"Value of custom quality--------------------"<<customQualityCutsID<<endl;

            if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
            else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");}
            else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(3.);}
            else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
            else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}
            else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.);}
            else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
            else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
            else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
            else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
            else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}
            else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4.);}
            else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}
            else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
            else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
            else if(customQualityCutsID==56){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
            else if(customQualityCutsID==60){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
            else if(customQualityCutsID==64){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}

    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma_%2.1fsigmaTOF",cutPiCandidate, nSigmaPart2TPC,nSigmaPart2TOF),trkQualityCut,cutPiCandidate,AliPID::kKaon,nSigmaPart2TPC,nSigmaPart2TOF);
 
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  

    // Defining output objects
    const Int_t dims = 6;
    Int_t useIM[dims] = {1, 0, 0, 0, isMC, isMC};
    TString name[dims] = {"Unlike", "Mixing", "LikePP", "LikeMM", "True", "Mother"};
    TString comp[dims] = {"PAIR", "MIX", "PAIR", "PAIR", "TRUE", "MOTHER"};
    TString output[dims] = {"SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE"};
    Char_t charge1[dims] = {'+', '+', '+', '-', '+', '+'};
    Char_t charge2[dims] = {'-', '-', '+', '-', '-', '-'};
    Int_t pdgCode[dims] = {333, 333, 333, 333, 333, 333};
    Double_t motherMass[dims] = {1.019461, 1.019461, 1.019461, 1.019461, 1.019461, 1.019461};

    for (Int_t i = 0; i < dims; i++)
    {
        if (!useIM[i])
            continue;
        AliRsnMiniOutput *out = task->CreateOutput(name[i].Data(), output[i].Data(), comp[i].Data());
        out->SetCutID(0, iCutK);
        out->SetCutID(1, iCutK);
        out->SetDaughter(0, AliRsnDaughter::kKaon);
        out->SetDaughter(1, AliRsnDaughter::kKaon);
        out->SetCharge(0, charge1[i]);
        out->SetCharge(1, charge2[i]);
        out->SetMotherPDG(pdgCode[i]);
        out->SetMotherMass(motherMass[i]);

        out->AddAxis(imID, 95, 0.985, 1.08);
        out->AddAxis(ptID, 1, 6., 10.);
       // if(!isPP ) out->AddAxis(multID,10,0.,100.);
       // else out->AddAxis(multID, 10, 0., 100.); 

       out->AddAxis(alID, 72, -0.5 * TMath::Pi(), 1.5 * TMath::Pi()); 
      //  out->AddAxis(ptlID, 15, 0., 30.); 
         out->AddAxis(detaID, 64, -1.6, 1.6);
       // out->AddAxis(philID, 18, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
       // out->AddAxis(phiID, 18, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
        
    }
    return kTRUE;
}

