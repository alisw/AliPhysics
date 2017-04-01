// $Id: AliJDataManager.cxx,v 1.13 2008/02/12 15:51:27 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJDataManager.cxx
  \brief
  \author J. Rak, D.J.Kim, B.S Chang (University of Jyvaskyla)
  \email: djkim@cc.jyu.fi
  \version $Revision: 1.13 $
  \date $Date: 2008/02/12 15:51:27 $
  */
////////////////////////////////////////////////////


#include  "AliJDataManager.h"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>

#include  "AliJConst.h"
#include  "AliJBaseTrack.h"
#include  "AliJTrack.h"
#include  "AliJMCTrack.h"
#include  "AliJPhoton.h"
#include  "AliJEventHeader.h"
// TODO #include  "AliESDVZERO.h"
#include "AliJHistogramInterface.h"

#include  "AliJCard.h"
#include  "AliJEventPool.h"

#include  "AliJRunHeader.h"
#include "AliJCorrelationInterface.h"
#include "AliJTrackCut.h"
#include <AliJRunTable.h>

//______________________________________________________________________________
AliJDataManager::AliJDataManager(AliJCard *inCard, AliJHistogramInterface *histin, AliJCorrelationInterface *corrin, Bool_t execLocal ):
    fChain(NULL), 
    fCard(inCard), 
    fhistos(histin), 
    fcorrelations(corrin),
    fRunHeader(NULL), 
    fEventHeader(NULL), 
    fEventHeaderList(NULL), 
    fTrackList(NULL), 
    fPhotonList(NULL), 
    fCellList(NULL), 
    fPhotonListRecalib(NULL), 
    fCellListRecalib(NULL), 
    fMCTrackList(NULL),
    fVZEROData(NULL),
    fRunInfoList(NULL),
    fhadronSelectionCut(0),
    fFilterMap(0),
    fFName(),
    fExecLocal(execLocal),
    fTriggerMask(0),
    fTrackCut(NULL)
{
  // constructor
    fChain = new TChain("JODTree");
    fhadronSelectionCut = int(fCard->Get("HadronSelectionCut"));
    

    fTriggerMask = fCard->Get("TriggerMask");

    fTrackCut = new AliJTrackCut;
}

//______________________________________________________________________________
AliJDataManager::~AliJDataManager(){
  // destructor
    if( fChain ) delete fChain;
    if( fTrackCut) delete fTrackCut;
}

//______________________________________________________________________________
AliJDataManager::AliJDataManager() :
    fChain(NULL), 
    fCard(NULL), 
    fhistos(NULL), 
    fcorrelations(NULL),
    fRunHeader(NULL), 
    fEventHeader(NULL), 
    fEventHeaderList(NULL), 
    fTrackList(NULL), 
    fPhotonList(NULL), 
    fCellList(NULL), 
    fPhotonListRecalib(NULL), 
    fCellListRecalib(NULL), 
    fMCTrackList(NULL),
    fVZEROData(NULL),
    fRunInfoList(NULL),
    fhadronSelectionCut(0),
    fFilterMap(0),
    fFName(),
    fExecLocal(true),
    fTriggerMask(0),
    fTrackCut(NULL)
{
  // default constructor
}

//______________________________________________________________________________
/* prohibitted
AliJDataManager::AliJDataManager(const AliJDataManager& obj) :
    fChain(obj.fChain), 
    fCard(obj.fCard), 
    fhistos(obj.fhistos), 
    fcorrelations(obj.fcorrelations),
    fRunHeader(obj.fRunHeader), 
    fEventHeader(obj.fEventHeader), 
    fEventHeaderList(obj.fEventHeaderList), 
    fTrackList(obj.fTrackList), 
    fPhotonList(obj.fPhotonList), 
    fCellList(obj.fCellList), 
    fPhotonListRecalib(obj.fPhotonListRecalib), 
    fCellListRecalib(obj.fCellListRecalib), 
    fMCTrackList(obj.fMCTrackList),
    fVZEROData(obj.fVZEROData),
    fRunInfoList(obj.fRunInfoList),
    fhadronSelectionCut(obj.fhadronSelectionCut),
    fFilterMap(obj.fFilterMap),
    fFName(obj.fFName),
    fExecLocal(obj.fExecLocal),
    fTriggerMask(obj.fTriggerMask),
    fTrackCut(obj.fTrackCut)
{
    // copy constructor TODO: proper handling of pointer data members
    JUNUSED(obj);
}

//______________________________________________________________________________
AliJDataManager& AliJDataManager::operator=(const AliJDataManager& obj){
  // equal sign TODO: contents
  JUNUSED(obj);
  return *this;
}*/


//______________________________________________________________________________

bool AliJDataManager::IsGoodEvent(){    
    // event checker
    if(fEventHeader==NULL) return false;
    int  nContributorVtx =  fEventHeader->GetVtxMult();
    double zVert    = fEventHeader->GetZVertex();
    UInt_t triggermaskJCorran = fEventHeader->GetTriggerMaskJCorran();
    bool goodVertex = kFALSE;
    bool triggerred = (IsSelectedTrigger((int) triggermaskJCorran)); // CUT1

    bool aMB = triggermaskJCorran & ( 1<< kMinBiasTriggerBitJCorran  );
    bool aCentral = triggermaskJCorran & ( 1<<kCentralTriggerBitJCorran );
    bool aSemiCentral = triggermaskJCorran & ( 1<<kSemiCentralTriggerBitJCorran );

    /*
    for( int i=0;i<32;i++ ){
        cout<<(triggermaskJCorran&(1<<i)?1:0)<<" ";
    }
    cout<<endl;
    */

    if( aMB ) fhistos->fhEvents->Fill( 10 );
    if( aCentral ) fhistos->fhEvents->Fill( 11 );
    if( aSemiCentral ) fhistos->fhEvents->Fill( 12 );
    if( aMB || aCentral ) fhistos->fhEvents->Fill( 13 );
    if( aMB || aCentral || aSemiCentral ) fhistos->fhEvents->Fill( 14 );

    for( int i=0;i < 31 ;i++ ){
        if( fEventHeader->GetTriggerMaskAlice() & BIT(i) )
            fhistos->fhEventTrigger->Fill(i);
    }
    fhistos->fhEventTrigger->Fill(41);

    /*
       if(!AliJRunTable::GetInstance().IsHeavyIon()){ //pp data
       if(fRunHeader->GetRunNumber() >= 146686 && fRunHeader->GetRunNumber() <= 146860){ //p+p 2.76
       if(!(fEventHeader->GetTriggerMaskAlice() & (1<<13))) triggerred= kFALSE; //noSDD    CINT1-B-NOPF-FASTNOTRD
    // Check with BS again about woSDD
    }
    }
    */
    // pPb run
    /*
       if(AliJRunTable::GetInstance().IsPA()){ //p+Pb
       if(
    //(fEventHeader->GetTriggerMaskAlice() & (1<<20)) ) // kINT7 MB
    ////(fEventHeader->GetTriggerMaskJCorran() & (1<<3)) ) // Emc1Gamma
    (fEventHeader->GetTriggerMaskJCorran() & (1<<4)) ) // kEMCEJE
    triggerred= kTRUE; //
    // Check with BS again about woSDD
    }
    */
    //6   CINT7-B-NOPF-ALLNOTRD
    //7   CINT7-ACE-NOPF-ALLNOTRD
    //8   CINT7-A-NOPF-ALLNOTRD
    //9   CINT7-C-NOPF-ALLNOTRD
    //17   CEMC7EG1-B-NOPF-ALLNOTRD
    //18   CEMC7EG2-B-NOPF-ALLNOTRD
    //19   CEMC7EJ1-B-NOPF-ALLNOTRD
    //20   CEMC7EJ2-B-NOPF-ALLNOTRD

    //Emc1GammaTriggerBitJCorran

    if(triggerred){
        fhistos->fhEvents->Fill( 1 );
        if(nContributorVtx==0){
            fhistos->fhVertexZTriggVtx->Fill(nContributorVtx,0.);
        }else{
            fhistos->fhEvents->Fill( 2 );
            fhistos->fhZVertRaw->Fill(zVert);
            fhistos->fhVertexZTriggVtx->Fill(nContributorVtx,zVert);
            goodVertex = (bool) fCard->VertInZRange(zVert);

            if( goodVertex )
                fhistos->fhEvents->Fill( 3 );
        }
    }
    //Trigger to be selected from the JCorran trigger mask is specified in the fCard
    //if(( IsSelectedTrigger((int) triggermaskJCorran )
    //     || ( fCard->MbTrigger((int) triggermaskJCorran ) && fCard->MixMBForPi0Mass() ))
    if(
            triggerred
            && goodVertex ){ //Cut2

        fhistos->fhEvents->Fill( 4 );
        return true;
    }else{
        return false; 
    }
}

//______________________________________________________________________________
void AliJDataManager::RegisterList(TClonesArray* listToFill, TClonesArray* listFromToFill, 
        int cBin, int zBin, particleType whatToFill){ 
    // corrType whatCorrType){

    // this is here just to silence the bloody rule checker
    JUNUSED( cBin );
    JUNUSED( zBin );
    JUNUSED(listFromToFill);
    JUNUSED( listToFill );

    int noIn=0, counter=0;
    //   Double_t pid[10]={0};
    //   bool isAliPion = 0;
    //   bool isAliKaon = 0;
    //   bool isAliProton = 0;

    switch (whatToFill) {

        case kJPhoton:
            break;
        case kJPizero:
            break;
        
        // Selection for hadrons
        case kJHadron:       
            noIn    = fTrackList->GetEntriesFast(); 
            counter = 0;
            {
                for(int ii=0;ii<noIn;ii++){ // loop for all tracks 
                    AliJTrack *cgl = (AliJTrack*)fTrackList->At(ii);
                    //if(fhadronSelectionCut == kTrackCutJFG || fhadronSelectionCut == kTrackCutHBT) cgl->SetUseTPCTrack();

                    // checking bit convention
                    for(int iTrackSelection=0; iTrackSelection<32; iTrackSelection++) {
                        if( cgl->IsFiltered(iTrackSelection) ) fhistos->fhTrackSelection->Fill( iTrackSelection );
                    }
                    for(int iTrackSelection=0; iTrackSelection<32; iTrackSelection++) {
                        if( cgl->TestBit(BIT(iTrackSelection)) ) fhistos->fhTrackSelection->Fill( iTrackSelection+50 );
                    }

                    if( 1 
                            //&& ( cgl->GetFilterMap() & GetFilterMap() ) 
                            && fTrackCut->IsSelected( cgl, fhadronSelectionCut )
                            && fTrackCut->SetMomentum( cgl, fhadronSelectionCut )
                            && fCard->IsInEtaRange(cgl->Eta()) 
                      ){  // 
                        cgl->SetID(ii);
                        cgl->SetParticleType(kJHadron);
                        new ((*listToFill)[counter++]) AliJTrack(*cgl);
                    }
                }
            }
            break;

        // Selection for protons
        case kJProton:       
            noIn    = fTrackList->GetEntriesFast(); 
            counter = 0;
            for(int ii=0;ii<noIn;ii++){ // loop for all tracks 
                AliJTrack *cgl = (AliJTrack*)fTrackList->At(ii);
                Double32_t prob = 0.9;//cgl->GetPID(AliJTrack::kProtonAliJ, AliJTrack::kTPCTOF);
                //cout << AliJTrack::kProtonAli <<"\t"<< AliJTrack::kTPCTOF <<"\t"<< prob<< endl;
                //Double32_t prob = cgl->GetPID(AliJTrack::AliJTrkPID(4), AliJTrack::AliJTrkPIDmethod(2));
                if( 1 
                        && fCard->IsInEtaRange(cgl->Eta()) 
                        && ( cgl->GetFilterMap() & GetFilterMap() ) 
                        && prob > 0.9
                  ){  // All cuts applied in the Train production stage
                    cgl->SetID(ii);
                    cgl->SetParticleType(kJHadron);
                    new ((*listToFill)[counter++]) AliJTrack(*cgl);
                }
            }
            break;
        
        // Selection for Monte Carlo truth particles
        case kJHadronMC:
            noIn    = fMCTrackList->GetEntriesFast();
            counter = 0;
            for(int ii=0;ii<noIn;ii++){ // loop for all tracks
                AliJMCTrack *cgl = (AliJMCTrack*)fMCTrackList->At(ii);
              
                // Note: The Monte Carlo track list is filled by AliJFilter.cxx class
                // It is already required in that class that MC tracks are physical primaries
                // This is done using IsPhysicalPrimary() method of AliAODMCParticle
              
                // Check that the track is a final state charged hadron
                // in the acceptance range of the analysis
                if(/*cgl->IsFinal() &&*/ cgl->IsCharged() && cgl->IsHadron()
                                  && fCard->IsInEtaRange(  cgl->Eta() )){
                    cgl->SetID(ii);
                    cgl->SetParticleType(kJHadronMC);
                    new ((*listToFill)[counter++]) AliJMCTrack(*cgl);
                }
            }
            break;

        default :
            cout<<"Unknown particle type in AliJDataManager.cxx/RegisterList()"<<endl;
            exit(0);
    }//switch for PID

    // make the indexing correct
    listToFill->Compress();
}


//______________________________________________________________________________
void AliJDataManager::ChainInputStream(const char* infileList){
    // chainer

    if( fExecLocal ){
      
        // read root nano data files in a list  
        char inFile[200];
        ifstream infiles(infileList);
        while ( infiles >> inFile){ 
            fChain->Add(inFile);
        }
        //fChain->Print();

        if(fChain->GetEntriesFast()<=0){
            cout<<"Empty chain from "<<infileList<<endl;
            exit(0);
        }
        cout<<Form("there are %d events.\n", (int)fChain->GetEntries())<<endl;

        // Load Branch
        if( fChain->FindBranch("HeaderList") )
            fChain->SetBranchAddress("HeaderList", &fEventHeaderList);
        if( fChain->FindBranch("TrackList") )
            fChain->SetBranchAddress("TrackList", &fTrackList);
        if( fChain->FindBranch("PhotonList") )
            fChain->SetBranchAddress("PhotonList", &fPhotonList);
        if( fChain->FindBranch("CaloCellList") )
            fChain->SetBranchAddress("CaloCellList", &fCellList);
        if( fChain->FindBranch("MCTrackList") )
            fChain->SetBranchAddress("MCTrackList", &fMCTrackList);

        if( AliJRunTable::GetInstance().IsHeavyIon() && fChain->FindBranch("AliESDVZERO")) { // if fevent-plane sources were stored
            fChain->SetBranchAddress("AliESDVZERO", &fVZEROData);
            // TODO for FMD, TZERO
        }

        // fChain->SetBranchAddress("AliJMCTrackList", &fMCTrackList);
        // Event Header
        //Read Run Header
        TFile* firstfile =  fChain->GetFile();
        cout<<"============================================DEBUG 10 "<<endl;
        TList *fUserInfo =(TList*)  firstfile->Get("RunInfo");
        fUserInfo->Print();
        fRunHeader=(AliJRunHeader*) fUserInfo->First();
    }
    else{
        fRunHeader=(AliJRunHeader*) fRunInfoList->First();
    }

    if(fRunHeader){
        fRunHeader->PrintOut();
    }

}


//______________________________________________________________________________
Double_t AliJDataManager::LoadEvent(int ievt){
    //clear clones array and counters
    //load the new fevent

    int v;

    v = 0;

    if( fExecLocal ){
        //fChain->GetEntry(ievt);
        TString fp;

        //if( fExperiment == kMC ){
        //  return fMCJDmg->GenerateEvent();
        //} else {
        v =  ((TTree*)fChain)->GetEntry(ievt);

        fp = fChain->GetCurrentFile()->GetPath();
        if( fp != fFName ){
            cout << fp.Data() << endl;
            fFName = fp;
        }

    }

    fEventHeader = dynamic_cast<AliJEventHeader*>(fEventHeaderList->At(0));

    return v;
    //}
} 





















