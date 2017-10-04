//
//  AliAnalysisTaskTPCCalBeauty.cxx
//  
//
//  Created by Erin Gauger
//
//

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskTPCCalBeauty.h"
#include "AliKFParticle.h"
#include "AliAODMCParticle.h"

//#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"

class AliAnalysisTaskTPCCalBeauty;

using namespace std;

ClassImp(AliAnalysisTaskTPCCalBeauty)

AliAnalysisTaskTPCCalBeauty::AliAnalysisTaskTPCCalBeauty() :
    AliAnalysisTaskSE(),
    fAOD(0),
    fMCarray(0),
    fMCparticle(0),
    fpidResponse(0),
    fOutputList(0),
    fMultSelection(0),
    fCentrality(-1),
    fCentralityMin(0),
    fCentralityMax(100),
    fEMCEG1(kFALSE),
    fDCalDG1(kFALSE),
    fFlagClsTypeEMC(kTRUE),
    fFlagClsTypeDCAL(kTRUE),
    fUseTender(kTRUE),
    fFlagULS(kFALSE),
    fFlagLS(kFALSE),
    fNevents(0),
    fVtX(0),
    fVtY(0),
    fVtZ(0),
    fTrkPtB4TC(0),
    fTrkPt(0),
    fTrkP(0),
    fTrkClsPhi(0),
    fTrkClsEta(0),
    fClsPhi(0),
    fClsEta(0),
    fClsEamDCal(0),
    fClsEamEMCal(0),
    fTrkPhi(0),
    fTrkEta(0),
    fdEdx(0),
    fCentCheck(0),
    fTrigCheck(0),
    fEMCTrkMatch(0),
    fInvmassLS(0),
    fInvmassULS(0),
    fPhotonicElecYield(0),
    fULSdcaBelow(0),
    fLSdcaBelow(0),
    fElectronSprs(0)
{
    //Root IO constructor, don't allocate memory here
}
//_____________________________________________________________________
AliAnalysisTaskTPCCalBeauty::AliAnalysisTaskTPCCalBeauty(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0),
    fMCarray(0),
    fMCparticle(0),
    fpidResponse(0),
    fOutputList(0),
    fMultSelection(0),
    fCentrality(-1),
    fCentralityMin(0),
    fCentralityMax(100),
    fEMCEG1(kFALSE),
    fDCalDG1(kFALSE),
    fFlagClsTypeEMC(kTRUE),
    fFlagClsTypeDCAL(kTRUE),
    fUseTender(kTRUE),
    fFlagULS(kFALSE),
    fFlagLS(kFALSE),
    fNevents(0),
    fVtX(0),
    fVtY(0),
    fVtZ(0),
    fTrkPtB4TC(0),
    fTrkPt(0),
    fTrkP(0),
    fTrkClsPhi(0),
    fTrkClsEta(0),
    fClsPhi(0),
    fClsEta(0),
    fClsEamDCal(0),
    fClsEamEMCal(0),
    fTrkPhi(0),
    fTrkEta(0),
    fdEdx(0),
    fCentCheck(0),
    fTrigCheck(0),
    fEMCTrkMatch(0),
    fInvmassLS(0),
    fInvmassULS(0),
    fPhotonicElecYield(0),
    fULSdcaBelow(0),
    fLSdcaBelow(0),
    fElectronSprs(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTPCCalBeauty::~AliAnalysisTaskTPCCalBeauty()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;
    }
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::UserCreateOutputObjects()
{
    //Weights for pho reco?
    
    /////////////////
    // Output List //
    /////////////////
    
    //create a new TList that owns its objects
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    //create our histos and add them to the list
    fNevents = new TH1F("fNevents", "No. of Events; Counts", 3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,">2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,">2 Trks, Vtx_{z}<10cm");
    
    fVtX = new TH1F("fVtX","X Vertex Position;Vtx_{X};Counts",1000,-50,50);
    fOutputList->Add(fVtX);
    
    fVtY = new TH1F("fVtY","Y Vertex Position;Vtx_{Y};Counts",1000,-50,50);
    fOutputList->Add(fVtY);
    
    fVtZ = new TH1F("fVtZ","Z Vertex Position;Vtx_{Z};Counts",1000,-50,50);
    fOutputList->Add(fVtZ);
    
    fTrkPtB4TC = new TH1F("fTrkPtB4TC","Track p_{T} Distribution before track cuts;p_{T} (GeV/c);Counts",1000,0,100);
    fOutputList->Add(fTrkPtB4TC);
    
    fTrkPt = new TH1F("fTrkPt","Track p_{T} Distribution;p_{T} (GeV/c);Counts",1000,0,100);
    fOutputList->Add(fTrkPt);
    
    fTrkP = new TH1F("fTrkP","Track p Distribution;p (GeV/c);Counts",1000,0,100);
    fOutputList->Add(fTrkP);
    
    fTrkClsPhi = new TH1F("fTrkClsPhi","Track and Cluster #Delta #phi Distribution;#Delta #phi;Counts",100,0,6.3);
    fOutputList->Add(fTrkClsPhi);
    
    fTrkClsEta = new TH1F("fTrkClsEta","Track and Cluster #Delta #eta Distribution;#Delta #eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fTrkClsEta);
    
    fClsPhi = new TH1F("fClsPhi","Cluster #phi Distribution;#phi;Counts",100,0,6.3);
    fOutputList->Add(fClsPhi);
    
    fClsEta = new TH1F("fClsEta","Cluster #eta Distribution;#eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fClsEta);
    
    fClsEamDCal = new TH1F("fClsEamDCal","Cluster Energy after track matching to DCal;Cluster E;Counts",500,0.,50);
    fOutputList->Add(fClsEamDCal);

    fClsEamEMCal = new TH1F("fClsEamEMCal","Cluster Energy after track matching to EMCal;Cluster E;Counts",500,0.,50.);
    fOutputList->Add(fClsEamEMCal);
    
    fTrkPhi = new TH1F("fTrkPhi","Track #phi Distribution after matching;#phi;Counts",100,0,6.3);
    fOutputList->Add(fTrkPhi);
    
    fTrkEta = new TH1F("fTrkEta","Track #eta Distribution after matching;#eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fTrkEta);
    
    fdEdx = new TH1F("fdEdx","Track dE/dx Distribution;dE/dx;Counts",500,0,160);
    fOutputList->Add(fdEdx);
    
    fCentCheck = new TH1F("fCentCheck","Event Centrality Distribution;Centrality;Counts",100,0,100);
    fOutputList->Add(fCentCheck);
    
    fTrigCheck = new TH1F("fTrigCheck", "No. of Events; Counts",3,-0.5,2.5);
    fOutputList->Add(fTrigCheck);
    fTrigCheck->GetXaxis()->SetBinLabel(1,"INT7");
    fTrigCheck->GetXaxis()->SetBinLabel(2,"EG1");
    fTrigCheck->GetXaxis()->SetBinLabel(3,"DGl");
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","EMCal cluster distance from closest track", 100, -0.3, 0.3,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);
    
    fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassLS);
    
    fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassULS);
    
    fPhotonicElecYield = new TH1F("fPhotonicElecYield","Photonic Elec Yield; p_{T}(GeV/c); yield;", 60,0,30.);
    fOutputList->Add(fPhotonicElecYield);
    
    fULSdcaBelow = new TH2F("fULSdcaBelow","ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCA; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fULSdcaBelow);
    
    fLSdcaBelow = new TH2F("fLSdcaBelow","LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCA; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fLSdcaBelow);
    
    //Electron THnSparse
    Int_t bins1[5]=  {280,  160, 100, 100,  200}; // pT;nSigma;eop;m20;DCA
    Double_t xmin1[5]={ 2,   -8,   0,   0, -0.2};
    Double_t xmax1[5]={30,    8,   2,   1,  0.2};
    fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;m20;DCA;",5,bins1,xmin1,xmax1);
    fOutputList->Add(fElectronSprs);
    
    //add the list to our output file
    PostData(1, fOutputList);
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::UserExec(Option_t*)
{
    //get an AOD event from the analysis manager
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    
    //check if there actually is an event
    if(!fAOD) return;
    
    //Needed to get MC info
    //if(fAOD)fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    
    //PID initialized
    fpidResponse = fInputHandler->GetPIDResponse();
    
    ////////////////
    // Centrality //
    ////////////////
    Bool_t pass = kFALSE;
    if(fCentralityMin > -0.5){
        fCentrality = CheckCentrality(fAOD,pass);
        if(!pass)return;
    }
    fCentCheck->Fill(fCentrality);
    
    ///////////////////
    // Trigger Check //
    ///////////////////
    TString firedTrigger;
    TString TriggerEG1("EG1");
    TString TriggerDG1("DG1");
    if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
    if(fEMCEG1){
        if(!firedTrigger.Contains(TriggerEG1))return;
        fTrigCheck->Fill(1);
    }else if(fDCalDG1){
        if(!firedTrigger.Contains(TriggerDG1))return;
        fTrigCheck->Fill(2);
    }else{fTrigCheck->Fill(0);}
    
    ////////////////
    // Mag. field //
    ////////////////
    
    Int_t MagSign = 1;
    if(fAOD->GetMagneticField()<0) MagSign = -1;
    
    //////////////////////////////
    // Event Vertex & Selection //
    //////////////////////////////
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    
    //Making n track and vertex cut
    fNevents->Fill(0);
    if(NcontV<2)return;
    fNevents->Fill(1);
    if(TMath::Abs(pVtx->GetZ())>10.0) return;
    fNevents->Fill(2);
    
    fVtX->Fill(pVtx->GetX());
    fVtY->Fill(pVtx->GetY());
    fVtZ->Fill(pVtx->GetZ());
    
    //////////////////////
    // Find Kink Mother //
    //////////////////////
    Int_t numberofvertices = 100;
    if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
        AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
        if(!aodvertex) continue;
        if(aodvertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[numberofmotherkink] = idmother;
            numberofmotherkink++;
        }
    }

    ////////////////
    // track loop //
    ////////////////
    Int_t nTracks(fAOD->GetNumberOfTracks());
    for(Int_t i=0; i<nTracks; i++){
        
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track) continue;
        
        fTrkPtB4TC->Fill(track->Pt());
        
        //////////////////////
        // Apply track cuts //
        //////////////////////
        if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //global cuts with loose DCA cut
        
        Bool_t kinkmotherpass = kTRUE;
        for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
            if(track->GetID() == listofmotherkink[kinkmother]) {
                kinkmotherpass = kFALSE;
                continue;
            }
        }
        if(!kinkmotherpass) continue; //kink rejection
        
        Double_t d0z0[2]={-999,-999}, cov[3];
        Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
        if(track->GetTPCNcls() < 80) continue;
        if(track->GetITSNcls() < 3) continue;
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
        
        double phiMatchIts = track->Phi();
        
        if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)){
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        }
        
        //fill the track histograms
        fTrkPt->Fill(track->Pt());
        fTrkP->Fill(track->P());
        fTrkPhi->Fill(track->Phi());
        fTrkEta->Fill(track->Eta());
        fdEdx->Fill(track->GetTPCsignal());
        
        /////////////////////
        // Get MC PID Info //
        /////////////////////
        /*Int_t fpidSort = 3;
        Int_t pdg = -999;
        Int_t pidM = -99;
        Int_t pid_ele = -99;
        Int_t ilabelM = -1;
        Int_t ilabel = TMath::Abs(track->GetLabel()); //get MC label of track
        if(ilabel>0 && fMCarray)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
            pdg = fMCparticle->GetPdgCode(); //get pid of track
            
            if(TMath::Abs(pdg)==11)pid_ele = 1.0; //if electron...
            if(pid_ele==1.0)FindMother2(fMCparticle, ilabelM, pidM, fpidSort);//get its mom
            
        }*/
        //if(pidM==443)continue; // remove enhanced J/psi in MC !
        //if(pidM==-99)continue; // remove e from no mother !
        
        /////////////////////
        // Fill DCA Sparse //
        /////////////////////
        /*
        if(pid_ele==1.0){
        Double_t fvalueDCA[8] = {-999,-999,-999}; //pT,SortElectronsbyMom,DCA
        fvalueDCA[0] = track->Pt();
        fvalueDCA[1] = fpidSort;
        fvalueDCA[2] = d0z0[0]*track->Charge()*MagSign;
        
            fSparseDCA->Fill(fvalueDCA);}
        */
        
        ///////////////////////////
        // Match tracks to EMCal //
        ///////////////////////////
        Int_t EMCalIndex = -1;
        EMCalIndex = track->GetEMCALcluster();
        if(EMCalIndex < 0) continue;
        
        AliVCluster *clustMatch=0x0;
        clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex);
        Double_t emcphi = -999, emceta=-999;
        Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
        
        if(clustMatch && clustMatch->IsEMCAL())
        {
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) continue;
            
            /////////////////////////////////
            //Select EMCAL or DCAL clusters//
            /////////////////////////////////
            Float_t  emcx[3]; // cluster pos
            clustMatch->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            emcphi = clustpos.Phi();
            emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            if(emcphi > 1.39 && emcphi < 3.265) {
                fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
                fClsEamEMCal->Fill(clustMatch->E());
            }
            if(emcphi > 4.53 && emcphi < 5.708) {
                fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
                fClsEamDCal->Fill(clustMatch->E());
            }
            
            //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            if(fClsTypeEMC) fClsEamEMCal->Fill(clustMatch->E());
            if(fClsTypeDCAL) fClsEamDCal->Fill(clustMatch->E());
            
            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            fTrkClsPhi->Fill(fPhiDiff);
            fTrkClsEta->Fill(fEtaDiff);
            fClsPhi->Fill(emcphi);
            fClsEta->Fill(emceta);
            
            /////////////////
            // fill sparse //
            /////////////////
            
            Double_t nsigma = -999;
            nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            Double_t EovP = (clustMatch->E())/(track->P());
            Double_t M20 = clustMatch->GetM20();
            Double_t M02 = clustMatch->GetM02();
            
            if(nsigma>-1 && nsigma<3) {
                if(M20>0.01 && M20<0.35) {
                    if(EovP>0.9 && EovP<1.2){
                        InvMassCheck(i, track, d0z0, MagSign);
                    }
                }
            }
            
            Double_t fvalueElectron[5] = {-999,-999,-999,-999,-999};
            fvalueElectron[0] = track->Pt();
            fvalueElectron[1] = nsigma;
            fvalueElectron[2] = EovP;
            fvalueElectron[3] = M20;
            fvalueElectron[4] = d0z0[0]*track->Charge()*MagSign;
            /*
            fvalueElectron[5] = -999;
            if (fClsTypeEMC){
                fvalueElectron[5] = 0; //0=EMCal, 1=DCal
            }
            if (fClsTypeDCAL){
                fvalueElectron[5] = 1; //0=EMCal, 1=DCal
            }*/
            fElectronSprs->Fill(fvalueElectron);
            
        }
        
    }
    
    //save the data gathered in this iteration
    PostData(1,fOutputList);
}
//___________________________________________
Double_t AliAnalysisTaskTPCCalBeauty::CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass)
{
    //check centrality, Run 2
    if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
    }
    
    //AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    //if(!header) AliFatal("Not a standard AOD");
    //fMultiplicity = header->GetRefMultiplicity();

    if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax)){
        //fCentralityNoPass->Fill(fCentrality);
        //fMultiplicityNoPass->Fill(fMultiplicity);
        centralitypass = kFALSE;
    }else{
        //fCentralityPass->Fill(fCentrality);
        //fMultiplicityPass->Fill(fMultiplicity);
        centralitypass = kTRUE;
    }
    
    return fCentrality;
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
    // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface

    phidiff = 999;
    etadiff = 999;
    
    if (!t||!v) return;

    Double_t veta = t->GetTrackEtaOnEMCal();
    Double_t vphi = t->GetTrackPhiOnEMCal();

    Float_t pos[3] = {0};
    v->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi();
    etadiff=veta-ceta;
    phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::FindMother(AliAODMCParticle* part, int &ilabelM, int &pidM, Int_t &fpidSort)
{
    //gets the pid of mother track
    
    ilabelM = part->GetMother(); //get MC label for Mom
    
    if(ilabelM>0){
        
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get MC particle
        pidM = TMath::Abs(partM->GetPdgCode()); //ask for the Mom's pid
    
        if(pidM>400 && pidM<499){
            //need to check grandma since B->D pretty frequently
            Int_t ilabelGM = partM->GetMother();//get MC for Grandma
            if(ilabelGM>0){
                AliAODMCParticle *partGM = (AliAODMCParticle*)fMCarray->At(ilabelGM);
                Int_t pidGM = TMath::Abs(partGM->GetPdgCode());
                if(pidGM>500 && pidGM<599){
                    //cout<<"Mother is D, Grandma is B, pidGM = "<<pidGM<<endl;
                    fpidSort = 2; //bin for B daughters
                }else if(pidGM>400 && pidGM<499){
                    //cout<<"Mother is D, Grandma is D, pidGM = "<<pidGM<<endl;
                    Int_t ilabelGGM = partGM->GetMother();//get MC for Great Grandma
                    if (ilabelGGM>0) {
                        AliAODMCParticle *partGGM = (AliAODMCParticle*)fMCarray->At(ilabelGGM);
                        Int_t pidGGM = TMath::Abs(partGGM->GetPdgCode());
                        if(pidGGM>500 && pidGGM<599){
                            //cout<<"Mother is D, Grandma is D, GGM is B, pidGGM = "<<pidGGM<<endl;
                            fpidSort = 2;
                        }else{
                            fpidSort = 1;
                            //cout<<"Mother is D, Grandma is D, GGM !=B, pidGGM = "<<pidGGM<<endl;
                        }
                    }else{
                        fpidSort = 1;
                        //cout<<"Mother is D, Grandma = D, no GGM, pidGM = "<<pidGM<<endl;
                    }
                }else{
                    //cout<<"Mother is D, Grandma != D, pidGM = "<<pidGM<<endl;
                    fpidSort = 1;
                }
            }else{
                //cout<<"Mother is D, no Grandma, pidM = "<<pidM<<endl;
                fpidSort = 1;
            }
        }
        else if(pidM>500 && pidM<599){
            //cout<<"Mother is B, pidM = "<<pidM<<endl;
            fpidSort = 2; //bin for B daughters
        }
        else if(pidM==22 || pidM==111 || pidM==221){
            //cout<<"Mother is pi, eta, or gamma, pidM = "<<pidM<<endl;
            fpidSort = 0; //bin for gamma, eta, pi0
        }else{
            //cout<<"Mother is other, pidM = "<<pidM<<endl;
            fpidSort = 3;
            //cout << "--------------PDG code for Mother is: "<<pidM<<"--------------"<<endl;
        }
        
    }else{
        fpidSort = 3;
        //cout<<"No mother"<<endl;
    }
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::FindMother2(AliAODMCParticle* part, int &ilabelM, int &pidM, Int_t &fpidSort)
{
    //gets the pid of mother track
    
    //Get the mother, grandma, and great grandma pid
    ilabelM = part->GetMother(); //get MC label for Mom
    if(ilabelM>0){
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
        pidM = TMath::Abs(partM->GetPdgCode()); //ask for the Mom's pid
        
        //sort according to mother
        if(pidM>500 && pidM<599){
            fpidSort = 2; //Mom is B
        }
        else if(pidM>400 && pidM<499){
            fpidSort = 1; //Mom is D
        }
        else if(pidM==22 || pidM==111 || pidM==221){
            fpidSort = 0; //Mom is pi,eta,gamma
        }
        else{
            fpidSort = 6; //Mom is something else
        }
        
        if(pidM==443){
            fpidSort = 3; //Mom is J/psi
        }
        
        if(pidM>4000 && pidM<4999){
            fpidSort = 4; //Mom is c Baryon
        }
        
        if(pidM>5000 && pidM<5999){
            fpidSort = 5; //Mom is b Baryon
        }
        
        Int_t ilabelGM = partM->GetMother();//get MC for Grandma
        if(ilabelGM>0 && fpidSort!=2){
            AliAODMCParticle *partGM = (AliAODMCParticle*)fMCarray->At(ilabelGM); // get GMa particle
            Int_t pidGM = TMath::Abs(partGM->GetPdgCode()); //ask for grandma's pid
            
            //check if grandma is B
            if(pidGM>500 && pidGM<599){
                fpidSort = 2; //GMa is B
            }
            
            Int_t ilabelGGM = partGM->GetMother();//get MC for Great Grandma
            if(ilabelGGM>0 && fpidSort!=2){
                AliAODMCParticle *partGGM = (AliAODMCParticle*)fMCarray->At(ilabelGGM); // get GGMa particle
                Int_t pidGGM = TMath::Abs(partGGM->GetPdgCode()); //ask for ggma's pid
                
                //check if great grandma is B
                if(pidGGM>500 && pidGGM<599){
                    fpidSort = 2; //GGMa is B
                }
            }
        }
    }else{
        fpidSort = 6; //No mother
    }
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::InvMassCheck(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign)
{
    // Flags photonic electrons with inv mass cut
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0Asso[2]={-999,-999}, covAsso[3];
    Double_t DCAxyCut = 0.25, DCAzCut = 1;
    Int_t fPDGe1 = 11, fPDGe2 = 11;
    
    Double_t ptAsso=-999., nsigmaAsso=-999.;
    Int_t chargeAsso=0;
    Int_t charge=track->Charge();
    Double_t mass=-999., width = -999;
    Int_t MassCorrect;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Int_t Nuls=0, Nls=0;
    
    Int_t ntracks = fAOD->GetNumberOfTracks();
    for (int jtrack=0; jtrack<ntracks; jtrack++) {
        if (jtrack==itrack) {continue;} //asso track != selected track
        
        fFlagLS=kFALSE;
        fFlagULS=kFALSE;
        
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(jtrack));
        if(!trackAsso) continue;
        if(!trackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(trackAsso->GetTPCNcls() < 80) continue;
        
        nsigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        ptAsso = trackAsso->Pt();
        chargeAsso = trackAsso->Charge();
        
        //Some cuts on the associated track
        if(ptAsso < 0.3) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigmaAsso < -3 || nsigmaAsso > 3) continue;
        
        if(trackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0Asso, covAsso))
            if(TMath::Abs(d0z0Asso[0]) > DCAxyCut || TMath::Abs(d0z0Asso[1]) > DCAzCut) continue;
        
        if(charge>0) fPDGe1 = -11; //-11 in PDG is for positron, just to be confusing
        if(chargeAsso>0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fAOD->GetMagneticField());
        
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        MassCorrect = recg.GetMass(mass,width); //returns 1, not the mass
        if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
        
        if(fFlagLS && mass<0.1) Nls++;
        if(fFlagULS && mass<0.1) Nuls++;
        
        if (fFlagULS && mass<0.1 && track->Pt()>1) {
            fULSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
        }else if(fFlagLS && mass<0.1 && track->Pt()>1){
            fLSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
        }
        
    }
    
    fPhotonicElecYield->Fill(track->Pt(),Nuls-Nls);
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________
