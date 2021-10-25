Bool_t ConfigPhiLeadingPb(
 AliRsnMiniAnalysisTask *task,
 Bool_t isMC = kFALSE, 
 Bool_t isPP = kFALSE, 
 Double_t nSigmaPart1TPC = -1, 
 Double_t nSigmaPart2TPC = -1,
 Double_t nSigmaPart1TOF = -1, 
 Double_t nSigmaPart2TOF = -1
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
    TString scheme;
    AliRsnCutSet *cutSetKaon = new AliRsnCutSet("kaonCutSet", AliRsnTarget::kDaughter);
    AliRsnCutTrackQuality *trkQualityCut = new AliRsnCutTrackQuality("trackQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE, kTRUE);
    cutSetKaon->AddCut(trkQualityCut);

    if (!scheme.IsNull())
        scheme += "&";
    scheme += trkQualityCut->GetName();
 
    if (nSigmaPart1TPC >= 0)
    {
        AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK", AliPID::kKaon, AliRsnCutPIDNSigma::kTPC);
        cutKTPC->SinglePIDRange(nSigmaPart1TPC);
        cutSetKaon->AddCut(cutKTPC);
        if (!scheme.IsNull())
            scheme += "&";
        scheme += cutKTPC->GetName();
    }
    cutSetKaon->SetCutScheme(scheme.Data());
    Int_t iTrackCutK = task->AddTrackCuts(cutSetKaon);


    // Defining output objects
    const Int_t dims = 6;
    Int_t useIM[dims] = {1, 1, 0, 0, isMC, isMC};
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
        out->SetCutID(0, iTrackCutK);
        out->SetCutID(1, iTrackCutK);
        out->SetDaughter(0, AliRsnDaughter::kKaon);
        out->SetDaughter(1, AliRsnDaughter::kKaon);
        out->SetCharge(0, charge1[i]);
        out->SetCharge(1, charge2[i]);
        out->SetMotherPDG(pdgCode[i]);
        out->SetMotherMass(motherMass[i]);

        out->AddAxis(imID, 95, 0.985, 1.08);
        out->AddAxis(ptID, 1, 2., 10.);
        if(!isPP ) out->AddAxis(multID,5,0.,50.);
        else out->AddAxis(multID, 10, 0., 100.); 

       // out->AddAxis(alID, 36, -0.5 * TMath::Pi(), 1.5 * TMath::Pi()); 
        out->AddAxis(ptlID, 1, 4., 30.); 
       //  out->AddAxis(detaID, 16, -1.6, 1.6);
        out->AddAxis(philID, 18, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
        out->AddAxis(phiID, 18, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
        
    }
    return kTRUE;
}

