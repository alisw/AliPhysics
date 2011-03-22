#include "TStopwatch.h"
#include "Riostream.h"
#include "TFile.h"

int runFlowOnTheFlyGlauber(Int_t nevents=100000, Int_t iseed=1)
{
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    gSystem->Load("libPWG2flowCommon");
    gSystem->Load("libPWG2flowTools");
    
    fMyTRandom3 = new TRandom3(iseed);   
    gRandom->SetSeed(iseed);
    TFile *outputFile = new TFile("flowOnTheFlyGlauber.root","RECREATE");
    
    TStopwatch timer;
    timer.Start();
    
    runOneCentrality(3,nevents);
    runOneCentrality(4,nevents);
    runOneCentrality(5,nevents);
    
    outputFile->Close();
    delete outputFile;
    
    cout<<endl; cout<<" ---- Fini ---- "<<endl; cout<<endl;
    timer.Stop(); timer.Print();
}

void runOneCentrality(Int_t centrality=3, Int_t nEvts=100)
{
    
    //centralities: 0: 0-5
    //              1: 5-10
    //              2: 10-20
    //              3: 20-30
    //              4: 30-40
    //              5: 40-50
    //              6: 50-60
    //              7: 60-70
    //              8: 70-80
    
    Double_t centPerc[] = {0,5,10,20,30,40,50,60,70,80};
    Double_t xCross[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t xCrossErr[9] = {0.};
    Double_t impact[] = {0., 3.5, 4.95, 6.98, 8.55, 9.88, 11.4, 12.09, 13.06, 13.97};
    Double_t npart[] = {1000, 356, 305, 220, 155, 105, 67, 40, 21, 10, 4, 0};
 
    Option_t *sysA="Pb"; 
    Option_t *sysB="Pb";
    Double_t signn=64; // inelastic nucleon nucleon cross section
    //const char *fname="GlauberMC_PbPb_ntuple.root"; // name output file
    //  AliGlauberMC::RunAndSaveNtuple(nevents,sysA,sysB,signn);
    // run the code to produce an ntuple:

    Double_t mind=0.4;
    Double_t r=6.62;
    Double_t a=0.546;
    
    AliGlauberMC glauber(sysA,sysB,signn);
    glauber.SetMinDistance(mind);
    glauber.Setr(r);
    glauber.Seta(a);
    glauber.SetDoPartProduction(kTRUE);
    glauber.SetBmin(impact[centrality]);
    glauber.SetBmax(impact[centrality+1]);    
    glauber.SetdNdEtaType(AliGlauberMC::kNBDSV);
    glauber.GetdNdEtaParam()[0] = 2.49;    //npp
    glauber.GetdNdEtaParam()[1] = 1.7;  //ratioSgm2Mu
    glauber.GetdNdEtaParam()[2] = 0.13; //xhard

    
    AliFlowAnalysisWithMCEventPlane* mcep2 = new AliFlowAnalysisWithMCEventPlane();
    mcep2->SetHarmonic(2); // default is v2
    mcep2->Init();
    AliFlowAnalysisWithMCEventPlane* mcep3 = new AliFlowAnalysisWithMCEventPlane();
    mcep3->SetHarmonic(3); // default is v3
    mcep3->Init();
    AliFlowAnalysisWithMCEventPlane* mcep4 = new AliFlowAnalysisWithMCEventPlane();
    mcep4->SetHarmonic(4); // default is v4
    mcep4->Init();
    AliFlowAnalysisWithMCEventPlane* mcep5 = new AliFlowAnalysisWithMCEventPlane();
    mcep5->SetHarmonic(5); // default is v5
    mcep5->Init();
    
    AliFlowAnalysisWithScalarProduct* sp2 = new AliFlowAnalysisWithScalarProduct();
    sp2->SetHarmonic(2);
    sp2->Init();
    AliFlowAnalysisWithScalarProduct* sp3 = new AliFlowAnalysisWithScalarProduct();
    sp3->SetHarmonic(3);
    sp3->Init();
    AliFlowAnalysisWithScalarProduct* sp4 = new AliFlowAnalysisWithScalarProduct();
    sp4->SetHarmonic(4);
    sp4->Init();
    AliFlowAnalysisWithScalarProduct* sp5 = new AliFlowAnalysisWithScalarProduct();
    sp5->SetHarmonic(5);
    sp5->Init();
    
    AliFlowAnalysisWithQCumulants* qc2 = new AliFlowAnalysisWithQCumulants();
    qc2->SetHarmonic(2);
    qc2->Init();
    AliFlowAnalysisWithQCumulants* qc3 = new AliFlowAnalysisWithQCumulants();
    qc3->SetHarmonic(3);
    qc3->Init();
    AliFlowAnalysisWithQCumulants* qc4 = new AliFlowAnalysisWithQCumulants();
    qc4->SetHarmonic(4);
    qc4->Init();
    AliFlowAnalysisWithQCumulants* qc5 = new AliFlowAnalysisWithQCumulants();
    qc5->SetHarmonic(5);
    qc5->Init();
    
    TH1F* histv2 = new TH1F(Form("v2 spread %.0f-%.0f",centPerc[centrality],centPerc[centrality+1]),"v_{2}",100,0,0.5);
    TH1F* histv3 = new TH1F(Form("v3 spread %.0f-%.0f",centPerc[centrality],centPerc[centrality+1]),"v_{3}",100,0,0.5);
    TH1F* histv4 = new TH1F(Form("v4 spread %.0f-%.0f",centPerc[centrality],centPerc[centrality+1]),"v_{4}",100,0,0.5);
    TH1F* histv5 = new TH1F(Form("v5 spread %.0f-%.0f",centPerc[centrality],centPerc[centrality+1]),"v_{5}",100,0,0.5);
    
    AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
    cutsRP->SetPtMin(0.2);
    cutsRP->SetPtMax(5.0);
    AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
    cutsPOI->SetPtMin(0.2);
    cutsPOI->SetPtMax(10.0);
    
    Printf("starting the main event loop..");
    for(Int_t i=0; i<nEvts; i++)
    {
        if (!glauber.NextEvent()) continue;
        Double_t mult = glauber.GetdNdEta();
        Double_t v2 = 0.2*glauber.GetEpsilon2Part(); histv2->Fill(v2);
        Double_t v3 = 0.18*glauber.GetEpsilon3Part(); histv3->Fill(v3);
        Double_t v4 = 0.09*glauber.GetEpsilon4Part(); histv4->Fill(v4);
        Double_t v5 = 0.03*glauber.GetEpsilon5Part(); histv5->Fill(v5);
        
        AliFlowEventSimple* event = new AliFlowEventSimple(mult,AliFlowEventSimple::kGenerate);
        event->SetMCReactionPlaneAngle(0.);
        event->AddFlow(0,v2,v3,v4,v5);
        event->TagRP(cutsRP);
        event->TagPOI(cutsPOI);
        
        // do flow analysis for various methods:
        mcep2->Make(event);
        mcep3->Make(event);
        mcep4->Make(event);
        mcep5->Make(event);
        sp2->Make(event);
        sp3->Make(event);
        sp4->Make(event);
        sp5->Make(event);
        qc2->Make(event);
        qc3->Make(event);
        qc4->Make(event);
        qc5->Make(event);
        cout <<"Event: " << i+1 << "\r"; cout.flush();
        delete event;
    } // end of for(Int_t i=0;i<nEvts;i++)
    
    // calculate the final results
    mcep2->Finish();
    mcep3->Finish();
    mcep4->Finish();
    mcep5->Finish();
    sp2->Finish();
    sp3->Finish();
    sp4->Finish();
    sp5->Finish();
    qc2->Finish();
    qc3->Finish();
    qc4->Finish();
    qc5->Finish();
    
    // open a new file which will hold the final results of all methods:
    const Int_t nMethods = 12;
    TString method[nMethods] = {"v2 MCEP","v3 MCEP","v4 MCEP","v5 MCEP","v2 SP","v3 SP","v4 SP","v5 SP","v2 QC","v3 QC","v4 QC","v5 QC"};
    TDirectoryFile *dirFileFinal[nMethods] = {NULL};
    TString fileName[nMethods];
    for(Int_t i=0; i<nMethods; i++)
    {
        // form a file name for each method:
        fileName[i]+=method[i].Data();
        fileName[i]+=Form(" %.0f-%.0f",centPerc[centrality],centPerc[centrality+1]);
        dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    }
    
    // store the final results
    mcep2->WriteHistograms(dirFileFinal[0]);
    mcep3->WriteHistograms(dirFileFinal[1]);
    mcep4->WriteHistograms(dirFileFinal[2]);
    mcep5->WriteHistograms(dirFileFinal[3]);
    sp2->WriteHistograms(dirFileFinal[4]);
    sp3->WriteHistograms(dirFileFinal[5]);
    sp4->WriteHistograms(dirFileFinal[6]);
    sp5->WriteHistograms(dirFileFinal[7]);
    qc2->WriteHistograms(dirFileFinal[8]);
    qc3->WriteHistograms(dirFileFinal[9]);
    qc4->WriteHistograms(dirFileFinal[10]);
    qc5->WriteHistograms(dirFileFinal[11]);
    
    histv2->Write();
    histv3->Write();
    histv4->Write();
    histv5->Write();
}

