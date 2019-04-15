Bool_t ConfigPhiLeading(AliRsnMiniAnalysisTask *task, Double_t nSigmaKaon = -1)
{

    // -- Values ------------------------------------------------------------------------------------
    /* invariant mass   */ Int_t imID = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
    /* transv. momentum */ Int_t ptID = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* angel to leading */ Int_t alID = task->CreateValue(AliRsnMiniValue::kAngleLeading,kFALSE);

    // Cuts
    AliRsnCutTrackQuality *trkQualityCut = new AliRsnCutTrackQuality("trackQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    AliRsnCutSetDaughterParticle *cutSet = new AliRsnCutSetDaughterParticle("noPIDCutSet", trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kKaon, nSigmaKaon);
    Int_t iTrackCutK = task->AddTrackCuts(cutSet);

    // Defining output objects
    const Int_t dims = 4;
    TString name[dims] = {"Unlike", "Mixing", "LikePP", "LikeMM"};
    TString comp[dims] = {"PAIR", "MIX", "PAIR", "PAIR"};
    TString output[dims] = {"SPARSE", "SPARSE", "SPARSE", "SPARSE"};
    Char_t charge1[dims] = {'+', '+', '+', '-'};
    Char_t charge2[dims] = {'-', '-', '+', '-'};
    Int_t pdgCode[dims] = {333, 333, 333, 333};
    Double_t motherMass[dims] = {1.019461, 1.019461, 1.019461, 1.019461};

    for (Int_t i = 0; i < dims; i++)
    {
        AliRsnMiniOutput *out = task->CreateOutput(name[i].Data(), output[i].Data(), comp[i].Data());
        out->SetCutID(0, iTrackCutK);
        out->SetCutID(1, iTrackCutK);
        out->SetDaughter(0, AliRsnDaughter::kKaon);
        out->SetDaughter(1, AliRsnDaughter::kKaon);
        out->SetCharge(0, charge1[i]);
        out->SetCharge(1, charge2[i]);
        out->SetMotherPDG(pdgCode[i]);
        out->SetMotherMass(motherMass[i]);

        out->AddAxis(imID, 215, 0.985, 1.2);
        out->AddAxis(ptID, 200, 0., 20.); //default use mother pt
        out->AddAxis(alID, 100, -0.5*TMath::Pi(), 1.5*TMath::Pi()); //-pi/2, 3/2pi
    }
    return kTRUE;
}
