
// #define V0Tracks_cxx
#include "V0Tracks.h"

void V0Tracks::setFile(TString fileName) {
    inputFileString = fileName;
}

TString V0Tracks::returnFile() {
    return inputFileString;
}

void V0Tracks::getTree() {
    inputFile = new TFile(inputFileString.Data());
    
    //inputFile->GetObject("V0Tracks", tree);
    inputFile->GetObject("TreeTRDPID", tree);
    
    initializeTree();
}

void V0Tracks::initializeTree() {
  /*
    tree->SetBranchAddress("bTRDNtracklets", &bTRDNtracklets, &b_bTRDNtracklets);
    tree->SetBranchAddress("bTRDslices[48]", bTRDslices, &b_bTRDslices);
    tree->SetBranchAddress("bTRDMomentum[6]", bTRDMomentum, &b_bTRDMomentum);
    tree->SetBranchAddress("bTRDY[6]", bTRDY, &b_bTRDY);
    //tree->SetBranchAddress("bTRDNcls", &bTRDNcls, &b_bTRDNcls);
    tree->SetBranchAddress("bTRDtheta", &bTRDtheta, &b_bTRDtheta);
    tree->SetBranchAddress("bTRDglobalphi", &bTRDglobalphi, &b_bTRDglobalphi);
    tree->SetBranchAddress("bCentrality", &bCentrality, &b_bCentrality);
    tree->SetBranchAddress("bTRDsignal", &bTRDsignal, &b_bTRDsignal);
    tree->SetBranchAddress("bTRDnclsdEdx", &bTRDnclsdEdx, &b_bTRDnclsdEdx);
    tree->SetBranchAddress("bTRDnch", &bTRDnch, &b_bTRDnch);
    tree->SetBranchAddress("bPDG", &bPDG, &b_bPDG);
    tree->SetBranchAddress("bNSigmaTPC[3]", bNSigmaTPC, &b_bNSigmaTPC);
    tree->SetBranchAddress("bNSigmaTOF[3]", bNSigmaTOF, &b_bNSigmaTOF);
  */
 
  tree->SetBranchAddress("run", &brun, &b_brun);     // not required but useful
  //   tree->SetBranchAddress("TRDNtracklets",     &bTRDNtracklets, &b_bTRDNtracklets);
  tree->SetBranchAddress("TRDslices[48]",     bTRDslices,     &b_bTRDslices); //  not used for truncated mean
  tree->SetBranchAddress("TRDMomentum[6]",    bTRDMomentum,   &b_bTRDMomentum); // required for TM
  tree->SetBranchAddress("TRDY[6]",           bTRDY,          &b_bTRDY); // required in older versions/cross checks
  tree->SetBranchAddress("TRDNcls",           &bTRDNcls,      &b_bTRDNcls); // general TRDNcls information, usually similar to the TM one, used for cross checks
  tree->SetBranchAddress("TRDtheta",          &bTRDtheta,     &b_bTRDtheta); // only for cross checks
  // tree->SetBranchAddress("TRDglobalphi",      &bTRDglobalphi,  &b_bTRDglobalphi); // cross check
  // tree->SetBranchAddress("TRDthetalayer[6]",  bTRDThetaLocal, &b_bTRDthetalocal); // cross check
  // tree->SetBranchAddress("TRDeta[6]",         bTRDEtaLocal,   &b_bTRDetalocal); // cross check
  tree->SetBranchAddress("TRDphi[6]",         bTRDPhiLocal,   &b_bTRDphilocal); // cross check
  tree->SetBranchAddress("centrality",        &bCentrality,   &b_bCentrality); // required at least for PbPb
  tree->SetBranchAddress("TRDsignal",         &bTRDsignal,     &b_bTRDsignal); // truncated mean signal
  tree->SetBranchAddress("TRDnclsdEdx",       &bTRDnclsdEdx,   &b_bTRDnclsdEdx);  // truncated mean Ncls (similar to general one)
  tree->SetBranchAddress("TRDnch",            &bTRDnch,        &b_bTRDnch); // trd #of tracklets
  tree->SetBranchAddress("PDG",               &bPDG,           &b_bPDG); // PDG information (V0 cuts)
  tree->SetBranchAddress("NSigmaTPC[3]",      bNSigmaTPC,     &b_bNSigmaTPC); 
  tree->SetBranchAddress("NSigmaTOF[3]",      bNSigmaTOF,     &b_bNSigmaTOF);
  tree->SetBranchAddress("TRDTPCtgl",         &bTRDTPCtgl,    &b_bTRDTPCtgl); // TPC dip angle used to correct eta dependence
}


Int_t V0Tracks::getNumberOfEntries() {
    return tree->GetEntriesFast();
}

Int_t V0Tracks::getEntry(Int_t number) {
    return tree->GetEntry(number);
}

Int_t V0Tracks::run() {
    return brun;
}

Int_t V0Tracks::trdNTracklets() {
    return bTRDNtracklets;
}

Double_t *V0Tracks::trdSlices() {
    Double_t *array[48] = {0};
    
    for (Int_t i = 0; i < 48; i++) {
        array[i] = &bTRDslices[i];
    }
    
    return *array;
}

Int_t V0Tracks::trdNCh() {
    return bTRDnch;
}

Int_t V0Tracks::trdNCls() {
  return bTRDnclsdEdx;
}

Int_t V0Tracks::trdNClsGeneral() {
  return bTRDNcls;
}

Double_t V0Tracks::trdSig() {
    return bTRDsignal;
}

Double_t *V0Tracks::trdMom() {
    Double_t *array[6] = {0};
    
    for (Int_t i = 0; i < 6; i++) {
        array[i] = &bTRDMomentum[i];
    }
    
    return *array;
}

Double_t *V0Tracks::trdY() {
    Double_t *array[6] = {0};
    
    for (Int_t i = 0; i < 6; i++) {
        array[i] = &bTRDY[i];
    }
    
    return *array;
}

Double_t V0Tracks::trdTheta() {
    return bTRDtheta;
}

Double_t V0Tracks::trdGlobalPhi() {
    return bTRDglobalphi;
}

Double_t *V0Tracks::trdThetaLocal() {
    Double_t *array[6] = {0};
    
    for (Int_t i = 0; i < 6; i++) {
        array[i] = &bTRDThetaLocal[i];
    }
    
    return *array;
}

Double_t *V0Tracks::trdEtaLocal() {
    Double_t *array[6] = {0};
    
    for (Int_t i = 0; i < 6; i++) {
        array[i] = &bTRDEtaLocal[i];
    }
    
    return *array;
}

Double_t *V0Tracks::trdPhiLocal() {
    Double_t *array[6] = {0};
    
    for (Int_t i = 0; i < 6; i++) {
        array[i] = &bTRDPhiLocal[i];
    }
    
    return *array;
}

Double_t V0Tracks::centrality() {
  bCentrality=1;
    return bCentrality;
}

Int_t V0Tracks::pId() {
    return bPDG;
}

Float_t *V0Tracks::nSigmaTPC() {
    Float_t *array[3] = {0};
    
    for (Int_t i = 0; i < 3; i++) {
        array[i] = &bNSigmaTPC[i];
    }
    
    return *array;
}

Float_t *V0Tracks::nSigmaTOF() {
    Float_t *array[3] = {0};
    
    for (Int_t i = 0; i < 3; i++) {
        array[i] = &bNSigmaTOF[i];
    }
    
    return *array;
}

Double_t V0Tracks::trdTPCtgl() {
    return bTRDTPCtgl;
}
