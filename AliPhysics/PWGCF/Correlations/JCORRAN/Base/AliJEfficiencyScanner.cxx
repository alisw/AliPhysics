//
//#include <Riostream.h>
// #include <TChain.h>
// #include <TVectorT.h> 
// #include <TVector3.h> 
// #include <TFile.h>
// #include <TH1.h> 
// #include <TObjArray.h>
// #include <TObjString.h>
// #include <TFormula.h>
// #include <TString.h>
// #include <TRefArray.h>
// #include <TNtuple.h>
// #include <TArrayF.h>

#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>

#include "AliJEfficiencyScanner.h" 
#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJPhoton.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"
#include <TFile.h>
#include <TRandom.h>


ClassImp(AliJEfficiencyScanner);

//______________________________________________________________________________
AliJEfficiencyScanner::AliJEfficiencyScanner() :   
    TNamed(),
    fMBTriggMask(0),
    fisIsolated(kFALSE),
    fisRelative(kTRUE),
    fisolParam(0.1),
    fisolCone(0.4),
    fTrackList(0),
    fMCTrackList(0x0),
    fEventHeader(0x0),
    fRunHeader(0x0),
    fhVertexZMC(0x0),
    fhVertexZTrigg(0x0),
    fhVertexZTriggVtx(0x0),
    fhVZRawMC(0x0),
    fhVZRecMC(0x0),
    fhVZRecAccMC(0x0),
    fh2VtxCent(0x0),
    fhL0Input(0x0),
    fhTriggerAlice(0x0),
    fhZVtxMCAll(0x0),
    fhZVtxMCTrigg(0x0),
    fhZVtxMCTriggVtx(0x0),
    fhZVtxRecAll(0x0),
    fhZVtxRecTrigg(0x0),
    fhZVtxRecTriggVtx(0x0),
    fVtxReFunc(0x0),
    fVtxMCFunc(0x0),
    fVtxRatioFunc(0x0),
    fVtxRatioMax(1)
{
    //Default constructor
    for( int ivtx=0;ivtx<kNVtxBin;ivtx++ ){
        for( int icent=0;icent<kNCentBin;icent++ ){
            fhChargedPtMCTriggVtx[ivtx][icent] = 0x0;
            for( int ifilter=0;ifilter<AliJTrackCut::kJNTrackCuts;ifilter++ ){
                for( int ipri=0;ipri<kNPrimaryStatus;ipri++ ){
                    for( int itt=0;itt<kNTrackType;itt++ ){
                        fhChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt] = 0x0;
                    }
                }
            }
        }
    }

    for(int i=0; i<AliJTrackCut::kJNTrackCuts; i++){
        fh2MultGenRawPrimary[i] = 0x0;
        fh2MultGenRawAll[i] = 0x0;
    }
        

}

//______________________________________________________________________________
AliJEfficiencyScanner::AliJEfficiencyScanner(const char *name):
    TNamed(name,name), 
    fMBTriggMask(0),
    fisIsolated(kFALSE),
    fisRelative(kTRUE),
    fisolParam(0.1),
    fisolCone(0.4),
    fTrackList(0),
    fMCTrackList(0x0),
    fEventHeader(0x0),
    fRunHeader(0x0),
    fhVertexZMC(0x0),
    fhVertexZTrigg(0x0),
    fhVertexZTriggVtx(0x0),
    fhVZRawMC(0x0),
    fhVZRecMC(0x0),
    fhVZRecAccMC(0x0),
    fh2VtxCent(0x0),
    fhL0Input(0x0),
    fhTriggerAlice(0x0),
    fhZVtxMCAll(0x0),
    fhZVtxMCTrigg(0x0),
    fhZVtxMCTriggVtx(0x0),
    fhZVtxRecAll(0x0),
    fhZVtxRecTrigg(0x0),
    fhZVtxRecTriggVtx(0x0),
    fVtxReFunc(0x0),
    fVtxMCFunc(0x0),
    fVtxRatioFunc(0x0),
    fVtxRatioMax(1)
{
    // Constructor
    for( int ivtx=0;ivtx<kNVtxBin;ivtx++ ){
        for( int icent=0;icent<kNCentBin;icent++ ){
            fhChargedPtMCTriggVtx[ivtx][icent] = 0x0;
            for( int ifilter=0;ifilter<AliJTrackCut::kJNTrackCuts;ifilter++ ){
                for( int ipri=0;ipri<kNPrimaryStatus;ipri++ ){
                    for( int itt=0;itt<kNTrackType;itt++ ){
                        fhChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt] = 0x0;
                    }
                }
            }
        }
    }
    for(int i=0; i<AliJTrackCut::kJNTrackCuts; i++){
        fh2MultGenRawPrimary[i] = 0x0;
        fh2MultGenRawAll[i] = 0x0;
    }
}

//____________________________________________________________________________
AliJEfficiencyScanner::AliJEfficiencyScanner(const AliJEfficiencyScanner& ap) :
    TNamed(ap.GetName(), ap.GetTitle()),
    fMBTriggMask(ap.fMBTriggMask),
    fisIsolated(ap.fisIsolated),
    fisRelative(ap.fisRelative),
    fisolParam(ap.fisolParam),
    fisolCone(ap.fisolCone),
    fTrackList(ap.fTrackList),
    fMCTrackList(ap.fMCTrackList),
    fEventHeader( ap.fEventHeader ),
    fRunHeader(ap.fRunHeader),
    fhVertexZMC(ap.fhVertexZMC),
    fhVertexZTrigg(ap.fhVertexZTrigg),
    fhVertexZTriggVtx(ap.fhVertexZTriggVtx),
    fhVZRawMC(ap.fhVZRawMC),
    fhVZRecMC(ap.fhVZRecMC),
    fhVZRecAccMC(ap.fhVZRecAccMC),
    fh2VtxCent(ap.fh2VtxCent),
    fhL0Input(ap.fhL0Input),
    fhTriggerAlice(ap.fhTriggerAlice),
    fhZVtxMCAll(ap.fhZVtxMCAll),
    fhZVtxMCTrigg(ap.fhZVtxMCTrigg),
    fhZVtxMCTriggVtx(ap.fhZVtxMCTriggVtx),
    fhZVtxRecAll(ap.fhZVtxRecAll),
    fhZVtxRecTrigg(ap.fhZVtxRecTrigg),
    fhZVtxRecTriggVtx(ap.fhZVtxRecTriggVtx),
    fVtxReFunc(ap.fVtxReFunc),
    fVtxMCFunc(ap.fVtxMCFunc),
    fVtxRatioFunc(ap.fVtxRatioFunc),
    fVtxRatioMax(ap.fVtxRatioMax)
{ 
    // cpy ctor
    for( int ivtx=0;ivtx<kNVtxBin;ivtx++ ){
        for( int icent=0;icent<kNCentBin;icent++ ){
            fhChargedPtMCTriggVtx[ivtx][icent] = ap.fhChargedPtMCTriggVtx[ivtx][icent];
            for( int ifilter=0;ifilter<AliJTrackCut::kJNTrackCuts;ifilter++ ){
                for( int ipri=0;ipri<kNPrimaryStatus;ipri++ ){
                    for( int itt=0;itt<kNTrackType;itt++ ){
                        fhChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt] = ap.fhChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt];
                    }
                }
            }
        }
    }
    for(int i=0; i<AliJTrackCut::kJNTrackCuts; i++){
        fh2MultGenRawPrimary[i] = ap.fh2MultGenRawPrimary[i];
        fh2MultGenRawAll[i] = ap.fh2MultGenRawAll[i];
    }

}

//_____________________________________________________________________________
AliJEfficiencyScanner& AliJEfficiencyScanner::operator = (const AliJEfficiencyScanner& ap)
{
    // assignment operator

    this->~AliJEfficiencyScanner();
    new(this) AliJEfficiencyScanner(ap);
    return *this;
}

//______________________________________________________________________________
AliJEfficiencyScanner::~AliJEfficiencyScanner()
{
    // destructor 
}

//________________________________________________________________________

void AliJEfficiencyScanner::UserCreateOutputObjects()
{  
    //=== create the jcorran outputs objects
    cout<<"DEBUG Start AliJEfficiencyScanner::UserCreateOutputObjects() "<<"\t"<<gDirectory<<endl;

    double ptbin[300] = {0};
    double pt = 0;
    int i = 0;
    for( i=0;i<300;i++ ){
        ptbin[i] = pt;
        if( pt > 50 ) break;
        if( pt < 3 ) pt+= 0.05;
        else if( pt < 5 ) pt+= 0.1;
        else if( pt < 10 ) pt+= 1;
        else pt+= 1;
    }
    cout<<"n Ptbin = "<<i<<endl;
    int nPtBin = i-1;

    fh2VtxCent = new TH2D("h2VtxCent","h2VtxCent",100,0,10,110,0,110 );
    fh2VtxCent->SetDirectory(gDirectory);

    int nVtxBin = kNVtxBin;
    int nCentBin = kNCentBin;
    if( fRunHeader->IsPP() ) nCentBin = 1;

    TString name ="";
    for( int ivtx=0;ivtx<nVtxBin;ivtx++ ){
        for( int icent=0;icent<nCentBin;icent++ ){
            fhChargedPtMC[ivtx][icent] = AddTH1D( Form("hChargedPtMC%02d%02d",ivtx,icent), 
                    new TH1D("","", nPtBin, ptbin) );
            fhChargedPtMCTriggVtx[ivtx][icent]=AddTH1D( 
                    Form("hChargedPtMCTriggVtx%02d%02d", ivtx, icent ),
                    new TH1D("", "", nPtBin, ptbin) );
            fhChargedPtMCTrigg[ivtx][icent]=AddTH1D( 
                    Form("hChargedPtMCTrigg%02d%02d", ivtx, icent ),
                    new TH1D("", "", nPtBin, ptbin) );
            //for any MC track filled with MC pt in triggered event
            name = Form("h2DChargedPtTrigg%02d%02d",ivtx,icent);
            fh2DChargedPtTrigg[ivtx][icent]=AddTH2D(name, new TH2D(name,name, nPtBin, ptbin, 20, -0.8, 0.8));
            //for any MC track filled with MC pt in triggered event with rec vertex
            name = Form("h2DChargedPtTriggVtx%02d%02d",ivtx,icent);
            fh2DChargedPtTriggVtx[ivtx][icent]=AddTH2D(name, new TH2D(name,name, nPtBin, ptbin, 20, -0.8, 0.8));

        }
    }
    fhTriggerAlice = AddTH1D("hTriggerAlice",new TH1D("","",32,0,32));
    fhL0Input = AddTH1D("hL0Input",new TH1D("","",32,0,32));
    // DCA bin
    double dcaBin[1000];
    double dbin=-50-1;
    int ndbin = 0;
    dcaBin[0] = dbin;
    double tol = 1e-5;
    while(dbin < 50-1){
        if( fabs(dbin) < 2-tol ) dbin+=0.01;
        else if( fabs(dbin) < 5-tol  ) dbin+= 0.05;
        else if( fabs(dbin) <= 10-tol ) dbin+= 0.1;
        else dbin += 1;
        if( fabs(dbin) < tol ) dbin=0;
        dcaBin[ndbin++] = dbin;
    }

    for( int ivtx=0;ivtx<nVtxBin;ivtx++ ){
        for( int icent=0;icent<nCentBin;icent++ ){
            for( int ifilter=0;ifilter<AliJTrackCut::kJNTrackCuts;ifilter++ ){

                name = Form("h2DChargedPtAll%02d%02d%02d",ivtx,icent,ifilter);
                fh2DChargedPtAll[ivtx][icent][ifilter]=AddTH2D(name, new TH2D(name,name, nPtBin, ptbin, 20, -0.8, 0.8));
                name = Form("h2DChargedPtRec%02d%02d%02d",ivtx,icent,ifilter);
                fh2DChargedPtRec[ivtx][icent][ifilter]=AddTH2D(name, new TH2D(name,name, nPtBin, ptbin, 20, -0.8, 0.8));
                for( int ipri=0;ipri<kNPrimaryStatus;ipri++ ){
                    for( int itt=0;itt<kNTrackType;itt++ ){
                        name = Form("hChargedPtMCRecoCentVtx%02d%02d%02d%02d%02d", ivtx, icent, ifilter, ipri, itt );
                        fhChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt]=AddTH1D(name, new TH1D( name, name, nPtBin, ptbin));
                    }
                }
                name = Form("hDCA2VertexXY%02d%02d%02d", ivtx, icent, ifilter );
                fhDCA2VertexXY[ivtx][icent][ifilter] = AddTH1D(name, new TH1D("","",ndbin-1, dcaBin ));
                name = Form("hDCA2VertexZ%02d%02d%02d", ivtx, icent, ifilter );
                fhDCA2VertexZ[ivtx][icent][ifilter] = AddTH1D( name, new TH1D("","",ndbin-1, dcaBin ));

            }
        }
    }
//  TH2D * fh2MultGenRawPrimary[AliJTrackCut::kJNTrackCuts];
//  TH2D * fh2MultGenRawAll[AliJTrackCut::kJNTrackCuts];

//  TH1D * fhDCA2VertexXY[kNVtxBin][kNCentBin][AliJTrackCut::kJNTrackCuts];
    for( int ifilter=0;ifilter<AliJTrackCut::kJNTrackCuts;ifilter++ ){
        name = Form("h2MultGenRawPrimary%02d", ifilter);
        fh2MultGenRawPrimary[ifilter] = AddTH2D( name, new TH2D(name, name, 100,0,100,100,0,100 ));
        name = Form("h2MultGenRawAll%02d", ifilter);
        fh2MultGenRawAll[ifilter] = AddTH2D( name, new TH2D(name, name, 100,0,100,100,0,100 ));
    }
    double   binsVertexMult[] = {0,1,2,3,4,5,10000};
    int   NbinsVertexMult  = sizeof(binsVertexMult)/sizeof(double)-1;
    double binsVertexZ[]    = {-10,-6,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,6,10};
    int   NbinsVertexZ   = sizeof(binsVertexZ)/sizeof(double)-1;
    fhVertexZMC = AddTH2D("hVertexZMC",new TH2D("hVertexZMC","hVertexZMC", NbinsVertexMult, binsVertexMult, NbinsVertexZ, binsVertexZ));
    fhVertexZTrigg = AddTH2D("hVertexZTrigg",new TH2D("hVertexZTrigg","hVertexZTrigg",  NbinsVertexMult, binsVertexMult, NbinsVertexZ, binsVertexZ));
    fhVertexZTriggVtx = AddTH2D("hVertexZTriggVtx",new TH2D("hVertexZTriggVtx","hVertexZTriggVtx",  NbinsVertexMult, binsVertexMult, NbinsVertexZ, binsVertexZ));
    fhVZRawMC = AddTH1D("hVtxZMC", new TH1D("hVtxZMC","VertexZ in MC ",200,-50,50));
    fhVZRecMC = AddTH1D("hVtxZRec", new TH1D("hVtxZRec","VertexZ Rec in MC",200,-50,50));
    fhVZRecAccMC = AddTH1D("hVtxZRecAcc", new TH1D("hVtxZRecAcc","VertexZ Rec in MC and acc",200,-50,50));

	int NvtxBins = 400;
	double vtxedge = 50;
    fhZVtxMCAll = AddTH1D("hZVtxMCAll", new TH1D("hZVtxMCAll","VetrexZ in MC for all event",NvtxBins, -vtxedge, vtxedge ));
    fhZVtxMCTrigg = AddTH1D("hZVtxMCTrigg", new TH1D("hZVtxMCTrigg","VetrexZ in MC for triggered event",NvtxBins, -vtxedge, vtxedge ));
    fhZVtxMCTriggVtx = AddTH1D("hZVtxMCTriggVtx", new TH1D("hZVtxMCTriggVtx","VetrexZ in MC for triggered and vtx reconstructed event",NvtxBins, -vtxedge, vtxedge ));

    fhZVtxRecAll = AddTH1D("hZVtxRecAll", new TH1D("hZVtxRecAll","VetrexZ in Rec for all event",NvtxBins, -vtxedge, vtxedge ));
    fhZVtxRecTrigg = AddTH1D("hZVtxRecTrigg", new TH1D("hZVtxRecTrigg","VetrexZ in Rec for triggered event",NvtxBins, -vtxedge, vtxedge ));
    fhZVtxRecTriggVtx = AddTH1D("hZVtxRecTriggVtx", new TH1D("hZVtxRecTriggVtx","VetrexZ in Rec for triggered and vtx reconstructed event",NvtxBins, -vtxedge, vtxedge ));

    double v0 = -20;
    double v1 = 20;
    fVtxReFunc = new TF1("VtxReFunc", "gaus",v0,v1);
    fVtxReFunc->SetParameters(1 ,  -7.14076e-01, 6.27110 );
    fVtxMCFunc = new TF1("VtxMCFunc", "gaus",v0,v1);
    fVtxMCFunc->SetParameters(1, -4.53674e-01,  5.27088e+00 );
    fVtxRatioFunc = new TF1("VtxRatioFunc", "VtxReFunc/VtxMCFunc",v0,v1);
    fVtxRatioMax = fVtxRatioFunc->GetMaximum();

    //cout << "Add(fAliJRunHeader) in UserCreateObject() ======= " << endl;
    cout<<"DEBUG END AliJEfficiencyScanner::UserCreateOutputObjects() "<<endl;

}

//______________________________________________________________________________

void AliJEfficiencyScanner::UserExec(Option_t *option)
{

    // Main loop
    // Called for each event
    JUNUSED(option);

    double zVtxCut = 10;
    double etaCut  = 0.8;
    double nCentBin = 20;

    AliJEventHeader * eventHeader = GetJEventHeader();

    //== TRIGGER( MB ) TODO is this?
    Bool_t triggeredEventMB  = eventHeader->GetTriggerMaskJCorran() & fMBTriggMask; // 1:kMB 8:kINT7 //TODO
    UInt_t trigAlice = eventHeader->GetTriggerMaskAlice();
/*
cout<<"Trigger"<<endl;
	for( int i=0;i<32 ;i++ ) cout<<Form("%02d",i)<<" ";
	cout<<endl;
	for( int i=0;i<32 ;i++ ) cout<<(eventHeader->GetTriggerMaskJCorran()&(1<<i)?1:0)<<"  ";
	cout<<endl;
	for( int i=0;i<32 ;i++ ) cout<<(eventHeader->GetTriggerMaskAlice()&(1<<i)?1:0)<<"  ";
	cout<<endl;
*/
    /*  
    bool isMB1   = trigAlice & (1UL<<0);
    bool isMBBG1 = trigAlice & (1UL<<3);
    bool isGFO = trigAlice & (1UL<<14);
    bool isMBBG3 = trigAlice & (1UL<<1); 
    bool isV0A   = trigAlice & (1UL<<10);
    bool isV0C   = trigAlice & (1UL<11);
    */

    // L0Input
    UInt_t l0Input = eventHeader->GetL0TriggerInputs();
    /* 
    bool isL0V0A  = l0Input & (1UL<<9);
    bool isL0V0C  = l0Input & (1UL<<10);
    bool isL0V0AND  = l0Input & (1UL<<4);
    bool isL0V0OR   = l0Input & (1UL<<1);
    bool isL0V0BeamGas   = l0Input & (1UL<<2);
    bool isL0SMB   = l0Input & (1UL<<3);
    bool isL0T0X   = l0Input & (1UL<<0);
    */


    
    for( UInt_t i=0;i<32;i++ ){
        if( l0Input  & 1UL<<i ) fhL0Input->Fill(i);
        if( trigAlice & 1UL<<i ) fhTriggerAlice->Fill(i);
/*
        if( isL0V0A != isV0A ) cout<<"J_WARN 0 : L0:V0A != Class:V0A "<<isL0V0A<<"\t"<<isV0A<<endl;
        if( isL0V0C != isV0C ) cout<<"J_WARN 1 : L0:V0C != Class:V0C "<<isL0V0C<<"\t"<<isV0C<<endl;
        if( isL0SMB != isGFO ) cout<<"J_WARN 2 : L0:SMB != Class:GFO "<<isL0SMB<<"\t"<<isGFO<<endl;
        if( isL0V0OR != isMBBG1 ) cout<<"J_WARN 3 : L0:V0OR != Class:MBBG1 "<<isL0V0OR<<"\t"<<isMBBG1<<endl;
        if( (isL0V0A||isL0V0C) != isMBBG1 ) cout<<"J_WARN 4 : L0:V0A|V0C != Class:MBBG1 "<<(isL0V0A||isL0V0C)<<"\t"<<isMBBG1<<endl;
        if( ( isL0T0X || isL0V0OR || isL0SMB && !isL0V0BeamGas ) != isMB1 ) cout<<"J_WARN 5 : L0:MB1 != Class:MB1 "<<isMB1<<endl;
        if( (isL0V0OR || isL0SMB && !isL0V0BeamGas ) != (( isMBBG1 || isGFO ) && isMB1 )) cout<<"J_WARN 6 : L0:MB != Class:MB "<<(isL0V0OR || isL0SMB && !isL0V0BeamGas)<<"\t"<<((isMBBG1 || isGFO ) && isMB1)<<endl; 
        if( (isL0V0OR || isL0SMB ) != (( isMBBG1 || isGFO ) )) cout<<"J_WARN 6 : L0:MB != Class:MB "<<(isL0V0OR || isL0SMB && !isL0V0BeamGas)<<"\t"<<((isMBBG1 || isGFO ) && isMB1)<<endl; 
*/

    }
    //* FOR TEST */
    //Bool_t triggeredEventMB  = ( isMBBG1 || isGFO ) && isMB1;
    //Bool_t triggeredEventMB  = isMB1;
    //Bool_t triggeredEventMB  = ( isMBBG1 || isGFO );
    //
    //== REAL MB with L0 INput
    //triggeredEventMB  = ( (isL0V0OR || isL0SMB) );
    //== CENTRALITY
    double centrality = eventHeader->GetCentrality();
    int iCent = int( centrality/(100./nCentBin));
    if( iCent< 0 || iCent>20 ) return;
    //if( centrality < 0 || centrality > 100 ) return;
    if( fRunHeader->IsPP() ) {iCent = 0; centrality=0;}

    //== Vertex
    //==== Reco
    int ncontributors = eventHeader->GetVtxMult();
    double zVtxRec = eventHeader->GetZVertex();
    Bool_t goodRecVertex = (ncontributors>0) && (fabs(zVtxRec)<=zVtxCut);
    int iVtx = 0;
    //==== MC
    double zVtxMC = eventHeader->GetZVertexMC();
    //==== vtx sampling
    /*
    if( zVtxMC > fVtxRatioFunc->GetXmin() && zVtxMC < fVtxRatioFunc->GetXmax() ){
        double vtxRatio = fVtxRatioFunc->Eval(zVtxMC)/fVtxRatioMax;
        double random = gRandom->Uniform(0,1);
        if( random > vtxRatio ) return;
    }
    */

    fhVZRawMC->Fill(zVtxMC);
    fhVertexZMC->Fill(ncontributors,zVtxMC);
    // for crosscheck MC input events
    fhZVtxMCAll->Fill(zVtxMC);
    fhZVtxRecAll->Fill(zVtxRec);

    if(triggeredEventMB){
        fhVertexZTrigg->Fill(ncontributors,zVtxMC);
       // for crosscheck passing trigger condition
        fhZVtxMCTrigg->Fill(zVtxMC);
        fhZVtxRecTrigg->Fill(zVtxRec);

        if(ncontributors>0){//reconstructed vertex
            fhVZRecMC->Fill(zVtxRec);
       		// for crosscheck passing trigger condition + reconstruted vtx
            fhZVtxMCTriggVtx->Fill(zVtxMC);
            fhZVtxRecTriggVtx->Fill(zVtxRec);
        }
        if(goodRecVertex){
            fhVertexZTriggVtx->Fill(ncontributors,zVtxMC);
            fhVZRecAccMC->Fill(zVtxRec); 
        }
    }

    fh2VtxCent->Fill( zVtxRec, centrality );

    // TODO ? MC or REC?

    int nRawMultPri = 0;
    int nGenMultPri[AliJTrackCut::kJNTrackCuts] = {0};
    int nGenMultAll[AliJTrackCut::kJNTrackCuts] = {0};

    //==============================================
    //  LOOP over MC
    //==============================================

    //if( fabs(zVtxMC) > zVtxCut ) return;
    int nMCTrks = 0;
    if( fRunHeader->IsMC() )
        nMCTrks = GetJMCTracks()->GetEntriesFast();
    for( int it=0; it< nMCTrks; it++ ){
        AliJMCTrack * track = GetJMCTrack( it );//Always IsPhysicalPrimary
        if( !track ) continue;
        if( !track->IsTrue( AliJMCTrack::kPrimary ) ) continue;
        // TODO need? if( ! track->IsFinal() ) continue;
        double eta = track->Eta();
         double pt = track->Pt();
        if( fabs(eta) > etaCut ) continue;
        if( ! track->IsCharged() ) continue;
        
        // Isolation check
	if( fisIsolated ){
	  // Isolated tracks must fullfill fiducial acceptance cut
	  if( fabs(eta) > etaCut - fisolCone ) continue;
	  double isolThreshold;
	  if( fisRelative ) isolThreshold = fisolParam * pt;
            else isolThreshold = fisolParam;
	  
          double isolSum = 0.0;
          for(int ia = 0; ia < nMCTrks; ia++){
            if ( ia == it ) continue;
            AliJMCTrack * assTrack = GetJMCTrack( ia );
            if( !assTrack ) continue;
	    if( !assTrack->IsTrue( AliJMCTrack::kPrimary ) ) continue;
	    if(!assTrack->IsCharged()) continue;
	    double assEta = assTrack->Eta();
            if( fabs(assEta) > etaCut ) continue;
            if( track->DeltaR(*assTrack) < fisolCone ) isolSum += assTrack->Pt();
          }
          if ( isolSum > isolThreshold ) continue;
	}
	
	nRawMultPri++;
        fhChargedPtMC[iVtx][iCent]->Fill(pt);     //ALL charged physical primaries
        if(triggeredEventMB){//triggered
            fhChargedPtMCTrigg[iVtx][iCent]->Fill(pt);
            fh2DChargedPtTrigg[iVtx][iCent]->Fill(pt, eta); 
            if(goodRecVertex){ //triggered+vertex  
                fhChargedPtMCTriggVtx[iVtx][iCent]->Fill(pt);
                fh2DChargedPtTriggVtx[iVtx][iCent]->Fill(pt, eta);
            }
        }
    }

    //--------------------------------------------------------------  
    //---- Reconstracted TRACKS 
    //--------------------------------------------------------------  
    const int nTrkCut = AliJTrackCut::kJNTrackCuts;
    if( ! triggeredEventMB ) return;
    if( ! goodRecVertex ) return;
    int nTrks   = GetJTracks()->GetEntriesFast();
    for(Int_t it = 0; it < nTrks; it++) {
        AliJTrack * track = GetJTrack(it);
        if( !track ) continue;
        bool isTriggered = track->GetFilterMap();
        if( !isTriggered ) continue;


        //== init Vars of Track
        double eta   = track->Eta();
        double etaTPC = eta;
        double etaGCG = eta;
        double ptRec = track->Pt();
        double ptTPC = track->PtTPC();
        if( fabs(ptTPC) < 1e-2 ) ptTPC=ptRec;
        else{
            TVector3 v(track->GetTPCTrack());
            etaTPC = v.Eta();
        }
        double ptGCG = track->PtGCG();
        if( fabs(ptGCG) < 1e-2 ) ptGCG=ptRec;
        else{
            TVector3 v(track->GetGCGTrack());
            etaGCG = v.Eta();
        }
        double ptMC  = -1;
        int iPrimary = kJFake;

        //== Skip tracks with out of Eta
        if( fabs(eta) > etaCut && fabs(etaTPC)>etaCut && fabs(etaGCG)>etaCut) continue;
	
        //== Find MC Info
        for( int imc=0;imc<nMCTrks;imc++ ){
            AliJMCTrack * mcTrack = GetJMCTrack(imc);
            if( !mcTrack || !mcTrack->IsTrue( AliJMCTrack::kPrimary ) ) continue;
            if( mcTrack && (TMath::Abs(track->GetLabel()) == TMath::Abs(mcTrack->GetLabel())) ){
                iPrimary = kJPhysicsPrimary;
                ptMC = mcTrack->Pt();
                break;
            }
        }

        //== FILL HIST
        for( int icut=0;icut<nTrkCut;icut++ ){

	  if( IsSelected( track ,icut) ){

		// ====== check isolation
		if( fisIsolated ){
		  // Isolated tracks must fullfill fiducial acceptance cut
		  if( fabs(eta) > etaCut - fisolCone ) continue;
		  double isolThreshold;
		  if( fisRelative ) isolThreshold = fisolParam * ptRec;
		    else isolThreshold = fisolParam;

		  double isolSum = 0.0;
		  for(int ia = 0; ia < nTrks; ia++){
		    if ( ia == it ) continue;
	    	    AliJTrack * assTrack = GetJTrack( ia );
	    	    if( !assTrack ) continue;
		    if( !IsSelected(assTrack,icut) ) continue;
		    double assEta = assTrack->Eta();
	    	    if( fabs(assEta) > etaCut ) continue;
	    	    if( track->DeltaR(*assTrack) < fisolCone ) isolSum += assTrack->Pt();
	  	  }
	  	  if ( isolSum > isolThreshold )continue;
                }
		  if( fabs(eta) < etaCut ){
                    fh2DChargedPtAll[iVtx][iCent][icut]->Fill(ptRec, eta); // IS THIS OK? Compare next line!
                    fh2DChargedPtRec[iVtx][iCent][icut]->Fill(ptMC, eta);  // IS THIS OK? Compare next line!
                    fhChargedPtMCRecoCentVtx[iVtx][iCent][icut][iPrimary][kJGlobal]->Fill( ptRec );
                    if( iPrimary == kJPhysicsPrimary ) nGenMultPri[icut]++;
                    nGenMultAll[icut]++;
                    //fhDCA2VertexXY[iVtx][iCent][icut]->Fill( track->GetDCAtoVertexXY() );
                    //fhDCA2VertexZ[iVtx][iCent][icut]->Fill( track->GetDCAtoVertexZ() );
                }
                if( fabs(etaTPC) < etaCut )
                    fhChargedPtMCRecoCentVtx[iVtx][iCent][icut][iPrimary][kJTPCOnly]->Fill( ptTPC );
                if( fabs(etaGCG) < etaCut )
                    fhChargedPtMCRecoCentVtx[iVtx][iCent][icut][iPrimary][kJGCG]->Fill( ptGCG );
                if( (icut==0 && fabs(etaTPC) < etaCut) || (icut!=0&&fabs(eta)<etaCut) )
                    fhChargedPtMCRecoCentVtx[iVtx][iCent][icut][iPrimary][kJMCTrack]->Fill( ptMC );
            }
        }
    }

    for( int icut=0;icut<nTrkCut;icut++ ){
          fh2MultGenRawPrimary[icut]->Fill( nRawMultPri, nGenMultPri[icut] );
          fh2MultGenRawAll[icut]->Fill( nRawMultPri, nGenMultAll[icut] );
    }
    //cout<<"DEBUG END AliJEfficiencyScanner::UserExec"<<endl;
}

//______________________________________________________________________________
void AliJEfficiencyScanner::Init()
{
    // Intialisation of parameters

}

//______________________________________________________________________________
void AliJEfficiencyScanner::Terminate(Option_t *)
{
    // termination
}


