#include <algorithm>

void printTrans(double trans[16])
{
    // calculate scaling factors:
    double sx = sqrt(trans[0]*trans[0]+trans[4]*trans[4]+trans[8]*trans[8]);
    double sy = sqrt(trans[1]*trans[1]+trans[5]*trans[5]+trans[9]*trans[9]);
    double sz = sqrt(trans[2]*trans[2]+trans[6]*trans[6]+trans[10]*trans[10]);
    
    // calculate rotation angles (way 1):
    //    double alpha = TMath::ATan2(trans[1],trans[0]);
    //    double beta = TMath::ATan2(-trans[2],sqrt(trans[6]*trans[6]+trans[10]*trans[10]));
    //    double gamma = TMath::ATan2(trans[6],trans[10]);
    
    // calculate rotation angles (way 2):
    //    double beta = TMath::ASin(trans[2]);
    //    double alpha = TMath::ASin(trans[1]/TMath::Cos(beta));
    //    double gamma = TMath::ASin(trans[6]/TMath::Cos(beta));
    
    
    // calculate rotation angles (way 3):
    double alpha,beta,gamma;
    Double_t d = trans[2]/sx;
    if(d>1) d=1; else if(d<-1) d=-1; // Fix numerical errors
    beta = TMath::ASin(d);
    Double_t cos = TMath::Cos(beta);
    if(TMath::Abs(cos) > 8.7e-6) {
        alpha = TMath::ATan2(trans[1], trans[0]);
        gamma = TMath::ATan2(trans[2]/sy, trans[10]/sz);
    } else {
        alpha = TMath::ATan2(trans[1]/sx, trans[5]/sy);
        gamma = 0;
    }
    
    // print:
    cout<<"translation( "<<trans[12]<<" , "<<trans[13]<<" , "<<trans[14]<<" )";
    cout<<"\tscaling:( "<<sx<<" , "<<sy<<" , "<<sz<<" )";
    //    cout<<"alpha:"<<alpha<<"\tbeta:"<<beta<<"\tgamma:"<<gamma<<" (rad)"<<endl;
    cout<<"\tangles ( "<<180/TMath::Pi()*alpha<<" , "<<180/TMath::Pi()*beta<<" , "<<180/TMath::Pi()*gamma<<" )"<<endl;
    
    if(trans[3]!=0 || trans[7]!=0 || trans[11]!=0 || trans[15]!=1){
        cout<<"something is wrong with trans matrix!!!"<<endl;
    }
}

//TCanvas *c1 = new TCanvas("c1","c1",1200,800);
int c1_iter=1;
TString outPath;

ofstream testFile("testFile.ma");

void extractGeom(string detector="")
{
    if(detector==""){
        cout<<"Give name of the detector as an argument:"<<endl;
        cout<<"TPC, TOF, ITS, MUON1, MUON2, MUON3, PHOS, HMPID, TRD, EMCAL"<<endl;
        return;
    }
    outPath = TString::Format("%s/Desktop/geom/",gSystem->Getenv("HOME"));
    

    testFile<<"//Maya ASCII 2016 scene"<<endl;
    testFile<<"//Name: muon_geom.ma"<<endl;
    testFile<<"//Last modified: Mon, Aug 17, 2015 08:42:10 AM"<<endl;
    testFile<<"//Codeset: UTF-8"<<endl;
    testFile<<"requires maya \"2016\";"<<endl;
    testFile<<"requires \"stereoCamera\" \"10.0\";"<<endl;
    testFile<<"requires \"stereoCamera\" \"10.0\";"<<endl;
    testFile<<"currentUnit -l centimeter -a degree -t film;"<<endl;
    testFile<<"fileInfo \"application\" \"maya\";"<<endl;
    testFile<<"fileInfo \"product\" \"Maya 2016\";"<<endl;
    testFile<<"fileInfo \"version\" \"2016\";"<<endl;
    testFile<<"fileInfo \"cutIdentifier\" \"201502261600-953408\";"<<endl;
    testFile<<"fileInfo \"osv\" \"Mac OS X 10.9.4\";"<<endl;
    testFile<<"fileInfo \"license\" \"student\";"<<endl;
    testFile<<"createNode transform -s -n \"persp\";"<<endl;
    testFile<<"rename -uid \"F15AEFCA-C948-1ADF-5E78-ECB41EAA6E8D\";"<<endl;
    testFile<<"setAttr \".v\" no;"<<endl;
    testFile<<"setAttr \".t\" -type \"double3\" -684.93632335393477 194.31614466048867 -719.2481819517842 ;"<<endl;
    testFile<<"setAttr \".r\" -type \"double3\" -8.9383504151193751 -90.799997482083043 5.0888874903416268e-14 ;"<<endl;
    testFile<<"createNode camera -s -n \"perspShape\" -p \"persp\";"<<endl;
    testFile<<"rename -uid \"3F80FA24-B74C-F837-2E7E-08859BD4F9BB\";"<<endl;
    testFile<<"setAttr -k off \".v\" no;"<<endl;
    testFile<<"setAttr \".fl\" 10.767029387320401;"<<endl;
    testFile<<"setAttr \".coi\" 1092.6344264363879;"<<endl;
    testFile<<"setAttr \".imn\" -type \"string\" \"persp\";"<<endl;
    testFile<<"setAttr \".den\" -type \"string\" \"persp_depth\";"<<endl;
    testFile<<"setAttr \".man\" -type \"string\" \"persp_mask\";"<<endl;
    testFile<<"setAttr \".hc\" -type \"string\" \"viewSet -p %camera\";"<<endl;
    testFile<<"createNode transform -n \"muon7\";"<<endl;

    
    if ( !testFile.is_open() ) {
        cout<<"\n\ncouldn't create output file"<<endl;
    }
    
    TFile *geomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo.root",gSystem->Getenv("ALICE_ROOT")));
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) geomFile->Get("Gentle");
    TList *list = gse->GetElements(); // list with all sub-detectors in file
    
    //    c1->Divide(3,3);
    
    // TPC
    if(detector=="TPC")
    {
        TEveGeoShapeExtract *tpc = list->FindObject("TPC");
        TEveGeoShapeExtract *tpcBarrel = tpc->GetElements()->FindObject("TPC_M_1");
        TGeoPcon *tpcShape = tpcBarrel->GetShape();
        
        // extract information about TPC from TGeoShape (TGeoPcon in this case)
        double tpc_bb_hl[3] = {tpcShape->GetDX(),tpcShape->GetDY(),tpcShape->GetDZ()};// BB half-length in X,Y,Z
        double *tpc_bb_origin = tpcShape->GetOrigin();// BB origin
        int tpc_n_seg = tpcShape->GetNsegments(); // number of segments (TPC sectors)
        int tpc_n_z = tpcShape->GetNz(); // number of planes in Z
        double tpc_phi_1 = tpcShape->GetPhi1();// starting angle
        double tpc_d_phi = tpcShape->GetDphi();// angle range
        double *tpc_r_max = tpcShape->GetRmax();// array of outer radii
        double *tpc_r_min = tpcShape->GetRmin();// array of inner radii
        double *tpc_z = tpcShape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nTPC params (polycone):"<<endl;
        cout<<"dx:"<<tpc_bb_hl[0]<<endl;
        cout<<"dy:"<<tpc_bb_hl[1]<<endl;
        cout<<"dz:"<<tpc_bb_hl[2]<<endl;
        cout<<"ox:"<<tpc_bb_origin[0]<<"\toy:"<<tpc_bb_origin[1]<<"\toz:"<<tpc_bb_origin[2]<<endl;
        cout<<"n seg:"<<tpc_n_seg<<endl;
        cout<<"n z:"<<tpc_n_z<<endl;
        cout<<"phi1:"<<tpc_phi_1<<endl;
        cout<<"dphi:"<<tpc_d_phi<<endl;
        for(int i=0;i<tpcShape->GetNz();i++){cout<<"r_min:"<<tpc_r_min[i]<<"\tr_max:"<<tpc_r_max[i]<<"\tz:"<<tpc_z[i]<<endl;}
        printTrans(tpcBarrel->GetTrans());
        cout<<"\n\n"<<endl;
        
        tpcShape->SaveAs(Form("%s/tpcShape.xml",outPath));
        tpcShape->SaveAs(Form("/Users/Jerus/Desktop/tpcShape.gdml",outPath));
    }
    else if(detector=="ITS")
    {
        // ITS
        TEveGeoShapeExtract *its = list->FindObject("ITS");
        TEveGeoShapeExtract *itsDets = its->GetElements()->FindObject("ITS_Dets");
        TEveGeoShapeExtract *its12 = itsDets->GetElements()->FindObject("IT12_1");
        TEveGeoShapeExtract *its34 = itsDets->GetElements()->FindObject("IT34_1");
        TEveGeoShapeExtract *its56 = itsDets->GetElements()->FindObject("IT56_1");
        
        TGeoTube *its12Shape = its12->GetShape();
        TGeoPcon *its34Shape = its34->GetShape();
        TGeoPcon *its56Shape = its56->GetShape();
        
        // extract information about ITS12 from TGeoShape (TGeoTube in this case)
        double its12_bb_hl[3] = {its12Shape->GetDX(),its12Shape->GetDY(),its12Shape->GetDZ()};// BB half-length in X,Y,Z
        double *its12_bb_origin = its12Shape->GetOrigin();// BB origin
        double its12_r_max = its12Shape->GetRmax();// outer radius
        double its12_r_min = its12Shape->GetRmin();// inner radius
        double its12_dz = its12Shape->GetDz();// half-lenght in Z
        
        cout<<"\n\nITS layers 1 & 2 params (tube):"<<endl;
        cout<<"dx:"<<its12_bb_hl[0]<<endl;
        cout<<"dy:"<<its12_bb_hl[1]<<endl;
        cout<<"dz:"<<its12_bb_hl[2]<<endl;
        cout<<"ox:"<<its12_bb_origin[0]<<"\toy:"<<its12_bb_origin[1]<<"\toz:"<<its12_bb_origin[2]<<endl;
        cout<<"r_min:"<<its12_r_min<<"\tr_max:"<<its12_r_max<<"\tz:"<<its12_dz<<endl;
        printTrans(its12->GetTrans());
        cout<<"\n\n"<<endl;
        
        // extract information about ITS34 from TGeoShape (TGeoPcon in this case)
        double its34_bb_hl[3] = {its34Shape->GetDX(),its34Shape->GetDY(),its34Shape->GetDZ()};// BB half-length in X,Y,Z
        double *its34_bb_origin = its34Shape->GetOrigin();// BB origin
        int its34_n_seg = its34Shape->GetNsegments(); // number of segments (ITS sectors)
        int its34_n_z = its34Shape->GetNz(); // number of planes in Z
        double its34_phi_1 = its34Shape->GetPhi1();// starting angle
        double its34_d_phi = its34Shape->GetDphi();// angle range
        double *its34_r_max = its34Shape->GetRmax();// array of outer radii
        double *its34_r_min = its34Shape->GetRmin();// array of inner radii
        double *its34_z = its34Shape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nITS layers 3 & 4 params (polycone):"<<endl;
        cout<<"dx:"<<its34_bb_hl[0]<<endl;
        cout<<"dy:"<<its34_bb_hl[1]<<endl;
        cout<<"dz:"<<its34_bb_hl[2]<<endl;
        cout<<"ox:"<<its34_bb_origin[0]<<"\toy:"<<its34_bb_origin[1]<<"\toz:"<<its34_bb_origin[2]<<endl;
        cout<<"n seg:"<<its34_n_seg<<endl;
        cout<<"n z:"<<its34_n_z<<endl;
        cout<<"phi1:"<<its34_phi_1<<endl;
        cout<<"dphi:"<<its34_d_phi<<endl;
        for(int i=0;i<its34_n_z;i++){cout<<"r_min:"<<its34_r_min[i]<<"\tr_max:"<<its34_r_max[i]<<"\tz:"<<its34_z[i]<<endl;}
        printTrans(its34->GetTrans());
        cout<<"\n\n"<<endl;
        
        // extract information about ITS56 from TGeoShape (TGeoPcon in this case)
        double its56_bb_hl[3] = {its56Shape->GetDX(),its56Shape->GetDY(),its56Shape->GetDZ()};// BB half-length in X,Y,Z
        double *its56_bb_origin = its56Shape->GetOrigin();// BB origin
        int its56_n_seg = its56Shape->GetNsegments(); // number of segments (ITS sectors)
        int its56_n_z = its56Shape->GetNz(); // number of planes in Z
        double its56_phi_1 = its56Shape->GetPhi1();// starting angle
        double its56_d_phi = its56Shape->GetDphi();// angle range
        double *its56_r_max = its56Shape->GetRmax();// array of outer radii
        double *its56_r_min = its56Shape->GetRmin();// array of inner radii
        double *its56_z = its56Shape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nITS layers 5 & 6 params (polycone):"<<endl;
        cout<<"dx:"<<its56_bb_hl[0]<<endl;
        cout<<"dy:"<<its56_bb_hl[1]<<endl;
        cout<<"dz:"<<its56_bb_hl[2]<<endl;
        cout<<"ox:"<<its56_bb_origin[0]<<"\toy:"<<its56_bb_origin[1]<<"\toz:"<<its56_bb_origin[2]<<endl;
        cout<<"n seg:"<<its56_n_seg<<endl;
        cout<<"n z:"<<its56_n_z<<endl;
        cout<<"phi1:"<<its56_phi_1<<endl;
        cout<<"dphi:"<<its56_d_phi<<endl;
        for(int i=0;i<its56_n_z;i++){cout<<"r_min:"<<its56_r_min[i]<<"\tr_max:"<<its56_r_max[i]<<"\tz:"<<its56_z[i]<<endl;}
        printTrans(its56->GetTrans());
        cout<<"\n\n"<<endl;
        
        
        its12Shape->SaveAs(Form("%s/its12Shape.xml",outPath));
        its34Shape->SaveAs(Form("%s/its34Shape.xml",outPath));
        its56Shape->SaveAs(Form("%s/its56Shape.xml",outPath));
    }
    else if(detector=="TOF")
    {
        // TRD & TOF
        TEveGeoShapeExtract *tof = list->FindObject("TRD+TOF");
        TEveGeoShapeExtract *tof_outer = tof->GetElements()->FindObject("B076_1");
        TEveGeoShapeExtract *tof_inner = tof->GetElements()->FindObject("BREF_1");
        
        TGeoPgon *tofOuterShape = tof_outer->GetShape();
        TGeoPgon *tofInnerShape = tof_inner->GetShape();
        
        
        // extract information about TOF outer from TGeoShape (TGeoPgon in this case)
        double tofOuter_bb_hl[3] = {tofOuterShape->GetDX(),tofOuterShape->GetDY(),tofOuterShape->GetDZ()};// BB half-length in X,Y,Z
        double *tofOuter_bb_origin = tofOuterShape->GetOrigin();// BB origin
        int tofOuter_n_seg = tofOuterShape->GetNsegments(); // number of segments (TPC sectors)
        int tofOuter_n_z = tofOuterShape->GetNz(); // number of planes in Z
        double tofOuter_phi_1 = tofOuterShape->GetPhi1();// starting angle
        double tofOuter_d_phi = tofOuterShape->GetDphi();// angle range
        double *tofOuter_r_max = tofOuterShape->GetRmax();// array of outer radii
        double *tofOuter_r_min = tofOuterShape->GetRmin();// array of inner radii
        double *tofOuter_z = tofOuterShape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nTOF outer barrel params (polygon):"<<endl;
        cout<<"dx:"<<tofOuter_bb_hl[0]<<endl;
        cout<<"dy:"<<tofOuter_bb_hl[1]<<endl;
        cout<<"dz:"<<tofOuter_bb_hl[2]<<endl;
        cout<<"ox:"<<tofOuter_bb_origin[0]<<"\toy:"<<tofOuter_bb_origin[1]<<"\toz:"<<tofOuter_bb_origin[2]<<endl;
        cout<<"n seg:"<<tofOuter_n_seg<<endl;
        cout<<"n z:"<<tofOuter_n_z<<endl;
        cout<<"phi1:"<<tofOuter_phi_1<<endl;
        cout<<"dphi:"<<tofOuter_d_phi<<endl;
        for(int i=0;i<tofOuter_n_z;i++){cout<<"r_min:"<<tofOuter_r_min[i]<<"\tr_max:"<<tofOuter_r_max[i]<<"\tz:"<<tofOuter_z[i]<<endl;}
        printTrans(tof_outer->GetTrans());
        cout<<"\n\n"<<endl;
        
        // extract information about TOF inner from TGeoShape (TGeoPgon in this case)
        double tofInner_bb_hl[3] = {tofInnerShape->GetDX(),tofInnerShape->GetDY(),tofInnerShape->GetDZ()};// BB half-length in X,Y,Z
        double *tofInner_bb_origin = tofInnerShape->GetOrigin();// BB origin
        int tofInner_n_seg = tofInnerShape->GetNsegments(); // number of segments (TPC sectors)
        int tofInner_n_z = tofInnerShape->GetNz(); // number of planes in Z
        double tofInner_phi_1 = tofInnerShape->GetPhi1();// starting angle
        double tofInner_d_phi = tofInnerShape->GetDphi();// angle range
        double *tofInner_r_max = tofInnerShape->GetRmax();// array of outer radii
        double *tofInner_r_min = tofInnerShape->GetRmin();// array of inner radii
        double *tofInner_z = tofInnerShape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nTOF inner barrel params (polygon):"<<endl;
        cout<<"dx:"<<tofInner_bb_hl[0]<<endl;
        cout<<"dy:"<<tofInner_bb_hl[1]<<endl;
        cout<<"dz:"<<tofInner_bb_hl[2]<<endl;
        cout<<"ox:"<<tofInner_bb_origin[0]<<"\toy:"<<tofInner_bb_origin[1]<<"\toz:"<<tofInner_bb_origin[2]<<endl;
        cout<<"n seg:"<<tofInner_n_seg<<endl;
        cout<<"n z:"<<tofInner_n_z<<endl;
        cout<<"phi1:"<<tofInner_phi_1<<endl;
        cout<<"dphi:"<<tofInner_d_phi<<endl;
        for(int i=0;i<tofInner_n_z;i++){cout<<"r_min:"<<tofInner_r_min[i]<<"\tr_max:"<<tofInner_r_max[i]<<"\tz:"<<tofInner_z[i]<<endl;}
        printTrans(tof_inner->GetTrans());
        cout<<"\n\n"<<endl;
        
        tofOuterShape->SaveAs(Form("%s/tofOuterShape.xml",outPath));
        tofInnerShape->SaveAs(Form("%s/tofInnerShape.xml",outPath));
    }
    else if(detector=="PHOS")
    {
        // PHOS
        TEveGeoShapeExtract *phos = list->FindObject("PHOS");
        TEveGeoShapeExtract *phosSegment[5];
        TGeoTrd1 *phosShape[5];//a trapezoid with only x length varying with z
        
        TEveGeoShapeExtract *test[5];
        
        for(int i=0;i<5;i++)
        {
            phosSegment[i] = (TEveGeoShapeExtract*)phos->GetElements()->FindObject(Form("PHOS_%i",i+1));
            phosShape[i] = (TGeoTrd1*)phosSegment[i]->GetShape();
            
            // extract information about PHOS segment from TGeoShape (TGeoTrd1 in this case)
            double phos_bb_hl[3] = {phosShape[i]->GetDX(),phosShape[i]->GetDY(),phosShape[i]->GetDZ()};// BB half-length in X,Y,Z
            double *phos_bb_origin = phosShape[i]->GetOrigin();// BB origin
            double phos_dx1 = phosShape[i]->GetDx1();
            double phos_dx2 = phosShape[i]->GetDx2();
            double phos_dy = phosShape[i]->GetDy();
            double phos_dz = phosShape[i]->GetDz();
            
            cout<<"\n\nPHOS "<<i+1<<" params (trapezoid):"<<endl;
            cout<<"dx:"<<phos_bb_hl[0]<<endl;
            cout<<"dy:"<<phos_bb_hl[1]<<endl;
            cout<<"dz:"<<phos_bb_hl[2]<<endl;
            cout<<"ox:"<<phos_bb_origin[0]<<"\toy:"<<phos_bb_origin[1]<<"\toz:"<<phos_bb_origin[2]<<endl;
            cout<<"dx1:"<<phos_dx1<<endl;
            cout<<"dx2:"<<phos_dx2<<endl;
            cout<<"dy:"<<phos_dy<<endl;
            cout<<"dz:"<<phos_dz<<endl;
            printTrans(phosSegment[i]->GetTrans());
            cout<<"\n\n"<<endl;
            
            phosShape[i]->SaveAs(Form("%s/phosSeg%iShape.xml",outPath,i+1));
        }
    }
    else if(detector=="MUON1")
    {
        TFile *muonGeomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo_muon.root",gSystem->Getenv("ALICE_ROOT")));
        TEveGeoShapeExtract* gsem = (TEveGeoShapeExtract*) muonGeomFile->Get("Gentle MUON");
        TList *muonList = gsem->GetElements(); // list with all sub-detectors in file
        
        TEveGeoShapeExtract *level1 = (TEveGeoShapeExtract*)muonList->FindObject("YOUT1_1");
        
        // for first element of level1 (YOUT1_1)
        TEveGeoShapeExtract *YOUT1_level2[4];
        TEveGeoShapeExtract *YOUT1_level3[4][4];
        TEveGeoShapeExtract *YOUT1_level4[4][4];
        TEveGeoShapeExtract *YOUT1_level5[4][4];
        TEveGeoPolyShape *YOUT1_muonShape[4][4];
        TGeoPgon *YOUT1_muonShapePgon;
        
        // extract shapes from YOUT1
        for(int l2=0;l2<4;l2++)
        {
            YOUT1_level2[l2] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject(Form("SC0%i_1",l2+1));
            
            for(int l3=0;l3<4;l3++)
            {
                if(l2<2)
                {
                    YOUT1_level3[l2][l3] = (TEveGeoShapeExtract*)YOUT1_level2[l2]->GetElements()->FindObject(Form("SE%i%i_1",l2+1,l3));
                    YOUT1_level4[l2][l3] = (TEveGeoShapeExtract*)YOUT1_level3[l2][l3]->GetElements()->FindObject(Form("SQM%i_%i",l2+1,l3+1));
                    YOUT1_level5[l2][l3] = (TEveGeoShapeExtract*)YOUT1_level4[l2][l3]->GetElements()->FindObject(Form("SA%iG_1",l2+1));
                    YOUT1_muonShape[l2][l3]= (TEveGeoPolyShape*)YOUT1_level5[l2][l3]->GetShape();
                }
                else
                {
                    YOUT1_level3[l2][l3] = (TEveGeoShapeExtract*)YOUT1_level2[l2]->GetElements()->FindObject(Form("SQM%i_%i",l2+1,l3+1));
                    YOUT1_level4[l2][l3] = (TEveGeoShapeExtract*)YOUT1_level3[l2][l3]->GetElements()->FindObject(Form("SC%iG1_1",l2+1));
                    YOUT1_muonShapePgon= (TGeoPgon*)YOUT1_level4[l2][l3]->GetShape();
                }
            }
        }
        // print info
        for(int l2=0;l2<4;l2++)
        {
            cout<<"\n\nl2:"<<l2<<endl;
            
            for(int l3=0;l3<4;l3++)
            {
                if(l2<2)
                {
                    // extract information about MOUN segment from TGeoShape (TEveGeoPolyShape in this case)
                    double muon_bb_hl[3] = {YOUT1_muonShape[l2][l3]->GetDX(),YOUT1_muonShape[l2][l3]->GetDY(),YOUT1_muonShape[l2][l3]->GetDZ()};// BB half-length in X,Y,Z
                    double *muon_bb_origin = YOUT1_muonShape[l2][l3]->GetOrigin();// BB origin
                    
                    cout<<"\n\nmuon "<<l3<<" params (box):"<<endl;
                    cout<<"BB half lenght ( "<<muon_bb_hl[0]<<" , "<<muon_bb_hl[1]<<" , "<<muon_bb_hl[2]<<" )"<<endl;
                    cout<<"BB lenght ( "<<2*muon_bb_hl[0]<<" , "<<2*muon_bb_hl[1]<<" , "<<2*muon_bb_hl[2]<<" )"<<endl;
                    cout<<"BB origin ("<<muon_bb_origin[0]<<" , "<<muon_bb_origin[1]<<" , "<<muon_bb_origin[2]<<" )"<<endl;
                    printTrans(YOUT1_level5[l2][l3]->GetTrans());
                    cout<<"\n\n"<<endl;
                    
                    YOUT1_muonShape[l2][l3]->SaveAs(Form("%s/muonSeg%i_%i_%iShape.xml",outPath,0,l2,l3));
                    printMuon(Form("%smuonSeg%i_%i_%iShape.xml",outPath,0,l2,l3),true,15,35,0.21);
                }
                else
                {
                    // extract information about MOUN segment from TGeoShape (TGeoPgon in this case)
                    double muon_bb_hl[3] = {YOUT1_muonShapePgon->GetDX(),YOUT1_muonShapePgon->GetDY(),YOUT1_muonShapePgon->GetDZ()};// BB half-length in X,Y,Z
                    double *muon_bb_origin = YOUT1_muonShapePgon->GetOrigin();// BB origin
                    int muon_n_seg = YOUT1_muonShapePgon->GetNsegments(); // number of segments (TPC sectors)
                    int muon_n_z = YOUT1_muonShapePgon->GetNz(); // number of planes in Z
                    double muon_phi_1 = YOUT1_muonShapePgon->GetPhi1();// starting angle
                    double muon_d_phi = YOUT1_muonShapePgon->GetDphi();// angle range
                    double *muon_r_max = YOUT1_muonShapePgon->GetRmax();// array of outer radii
                    double *muon_r_min = YOUT1_muonShapePgon->GetRmin();// array of inner radii
                    double *muon_z = YOUT1_muonShapePgon->GetZ();// array of Z positions of planes
                    
                    
                    cout<<"\n\nmuon "<<l3<<" params (polygon):"<<endl;
                    cout<<"BB half lenght ( "<<muon_bb_hl[0]<<" , "<<muon_bb_hl[1]<<" , "<<muon_bb_hl[2]<<" )"<<endl;
                    cout<<"BB lenght ( "<<2*muon_bb_hl[0]<<" , "<<2*muon_bb_hl[1]<<" , "<<2*muon_bb_hl[2]<<" )"<<endl;
                    cout<<"BB origin ("<<muon_bb_origin[0]<<" , "<<muon_bb_origin[1]<<" , "<<muon_bb_origin[2]<<" )"<<endl;
                    cout<<"n seg:"<<muon_n_seg<<endl;
                    cout<<"n z:"<<muon_n_z<<endl;
                    cout<<"phi1:"<<muon_phi_1<<endl;
                    cout<<"dphi:"<<muon_d_phi<<endl;
                    for(int i=0;i<muon_n_z;i++){cout<<"r_min:"<<muon_r_min[i]<<"\tr_max:"<<muon_r_max[i]<<"\tz:"<<muon_z[i]<<endl;}
                    printTrans(YOUT1_level4[l2][l3]->GetTrans());
                    cout<<"\n\n"<<endl;
                    
                    YOUT1_muonShapePgon->SaveAs(Form("%s/muonSeg%i_%i_%iShape.xml",outPath,0,l2,l3));
                }
                
            }
        }
    }
    else if(detector=="MUON2")
    {
        muon2();
    }
    else if(detector=="MUON3")
    {
        muon3();
    }
    else if(detector=="TRD")
    {
        TFile *trdGeomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo_trd.root",gSystem->Getenv("ALICE_ROOT")));
        TEveGeoShapeExtract* gset = (TEveGeoShapeExtract*) trdGeomFile->Get("Gentle TRD");
        TGeoPgon *trdShape = gset->GetShape();
        
        // extract information about TRD outer from TGeoShape (TGeoPgon in this case)
        double trd_bb_hl[3] = {trdShape->GetDX(),trdShape->GetDY(),trdShape->GetDZ()};// BB half-length in X,Y,Z
        double *trd_bb_origin = trdShape->GetOrigin();// BB origin
        int trd_n_seg = trdShape->GetNsegments(); // number of segments (TPC sectors)
        int trd_n_z = trdShape->GetNz(); // number of planes in Z
        double trd_phi_1 = trdShape->GetPhi1();// starting angle
        double trd_d_phi = trdShape->GetDphi();// angle range
        double *trd_r_max = trdShape->GetRmax();// array of outer radii
        double *trd_r_min = trdShape->GetRmin();// array of inner radii
        double *trd_z = trdShape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nTRD params (polygon):"<<endl;
        cout<<"dx:"<<trd_bb_hl[0]<<endl;
        cout<<"dy:"<<trd_bb_hl[1]<<endl;
        cout<<"dz:"<<trd_bb_hl[2]<<endl;
        cout<<"ox:"<<trd_bb_origin[0]<<"\toy:"<<trd_bb_origin[1]<<"\toz:"<<trd_bb_origin[2]<<endl;
        cout<<"n seg:"<<trd_n_seg<<endl;
        cout<<"n z:"<<trd_n_z<<endl;
        cout<<"phi1:"<<trd_phi_1<<endl;
        cout<<"dphi:"<<trd_d_phi<<endl;
        for(int i=0;i<trd_n_z;i++){cout<<"r_min:"<<trd_r_min[i]<<"\tr_max:"<<trd_r_max[i]<<"\tz:"<<trd_z[i]<<endl;}
        printTrans(gset->GetTrans());
        cout<<"\n\n"<<endl;
        
        trdShape->SaveAs(Form("%s/trdShape.xml",outPath));
    }
    else if(detector=="EMCAL")
    {
        TFile *emcalGeomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo_emcal.root",gSystem->Getenv("ALICE_ROOT")));
        TEveGeoShapeExtract* gsee = (TEveGeoShapeExtract*) emcalGeomFile->Get("Gentle EMCAL");
        TGeoPgon *emcalShape = gsee->GetShape();
        
        // extract information about EMCAL outer from TGeoShape (TGeoPgon in this case)
        double emcal_bb_hl[3] = {emcalShape->GetDX(),emcalShape->GetDY(),emcalShape->GetDZ()};// BB half-length in X,Y,Z
        double *emcal_bb_origin = emcalShape->GetOrigin();// BB origin
        int emcal_n_seg = emcalShape->GetNsegments(); // number of segments (TPC sectors)
        int emcal_n_z = emcalShape->GetNz(); // number of planes in Z
        double emcal_phi_1 = emcalShape->GetPhi1();// starting angle
        double emcal_d_phi = emcalShape->GetDphi();// angle range
        double *emcal_r_max = emcalShape->GetRmax();// array of outer radii
        double *emcal_r_min = emcalShape->GetRmin();// array of inner radii
        double *emcal_z = emcalShape->GetZ();// array of Z positions of planes
        
        cout<<"\n\nemcal params (polygon):"<<endl;
        cout<<"dx:"<<emcal_bb_hl[0]<<endl;
        cout<<"dy:"<<emcal_bb_hl[1]<<endl;
        cout<<"dz:"<<emcal_bb_hl[2]<<endl;
        cout<<"ox:"<<emcal_bb_origin[0]<<"\toy:"<<emcal_bb_origin[1]<<"\toz:"<<emcal_bb_origin[2]<<endl;
        cout<<"n seg:"<<emcal_n_seg<<endl;
        cout<<"n z:"<<emcal_n_z<<endl;
        cout<<"phi1:"<<emcal_phi_1<<endl;
        cout<<"dphi:"<<emcal_d_phi<<endl;
        for(int i=0;i<emcal_n_z;i++){cout<<"r_min:"<<emcal_r_min[i]<<"\tr_max:"<<emcal_r_max[i]<<"\tz:"<<emcal_z[i]<<endl;}
        printTrans(gsee->GetTrans());
        cout<<"\n\n"<<endl;
        
        emcalShape->SaveAs(Form("%s/emcalShape.xml",outPath));
    }
    else
    {
        cout<<"unknown detector name"<<endl;
        testFile.close();
        return;
    }
    testFile.close();
    
}

void muon2()
{
    //----------------------------------------------------------------------------------//
    //                                                                                  //
    //      WARNING!!!                                                                  //
    //                                                                                  //
    //      Looking into this function may cause you cry.                               //
    //      Do it on your own risk.                                                     //
    //                                                                                  //
    //----------------------------------------------------------------------------------//
    
    TFile *muonGeomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo_muon.root",gSystem->Getenv("ALICE_ROOT")));
    TEveGeoShapeExtract* gsem = (TEveGeoShapeExtract*) muonGeomFile->Get("Gentle MUON");
    TList *muonList = gsem->GetElements(); // list with all sub-detectors in file
    
    TEveGeoShapeExtract *level1 = (TEveGeoShapeExtract*)muonList->FindObject("DDIP_1");
    
    // for second element of level1 (DDIP)
    TEveGeoShapeExtract *DDIP_level2[4];
    TEveGeoShapeExtract *DDIP_level3[4][9];
    TEveGeoShapeExtract *DDIP_level4[4][9][4];
    TEveGeoPolyShape *DDIP_eveShapes[100];
    TGeoBBox *DDIP_bbShapes[100];
    TEveGeoShapeExtract *bbp[100];
    TEveGeoShapeExtract *esp[100];
    
    //DDPI
    DDIP_level2[0] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC05I_1");
    DDIP_level2[1] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC05O_1");
    DDIP_level2[2] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC06I_1");
    DDIP_level2[3] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC06O_1");
    
    // from SC05I_1
    DDIP_level3[0][0] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA13_1");
    DDIP_level3[0][1] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA14_1");
    DDIP_level3[0][2] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA12_1");
    DDIP_level3[0][3] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA15_1");
    DDIP_level3[0][4] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA11_1");
    DDIP_level3[0][5] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA16_1");
    DDIP_level3[0][6] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA10_1");
    DDIP_level3[0][7] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA17_1");
    DDIP_level3[0][8] = (TEveGeoShapeExtract*)DDIP_level2[0]->GetElements()->FindObject("SLA9_1");
    
    // from SC05O_1
    DDIP_level3[1][0] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA4_1");
    DDIP_level3[1][1] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA5_1");
    DDIP_level3[1][2] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA3_1");
    DDIP_level3[1][3] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA6_1");
    DDIP_level3[1][4] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA2_1");
    DDIP_level3[1][5] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA7_1");
    DDIP_level3[1][6] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA1_1");
    DDIP_level3[1][7] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA8_1");
    DDIP_level3[1][8] = (TEveGeoShapeExtract*)DDIP_level2[1]->GetElements()->FindObject("SLA0_1");
    
    // from SC06I_1
    DDIP_level3[2][0] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB13_1");
    DDIP_level3[2][1] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB14_1");
    DDIP_level3[2][2] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB12_1");
    DDIP_level3[2][3] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB15_1");
    DDIP_level3[2][4] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB11_1");
    DDIP_level3[2][5] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB16_1");
    DDIP_level3[2][6] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB10_1");
    DDIP_level3[2][7] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB17_1");
    DDIP_level3[2][8] = (TEveGeoShapeExtract*)DDIP_level2[2]->GetElements()->FindObject("SLB9_1");
    
    // from SC06O_1
    DDIP_level3[3][0] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB4_1");
    DDIP_level3[3][1] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB5_1");
    DDIP_level3[3][2] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB3_1");
    DDIP_level3[3][3] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB6_1");
    DDIP_level3[3][4] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB2_1");
    DDIP_level3[3][5] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB7_1");
    DDIP_level3[3][6] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB1_1");
    DDIP_level3[3][7] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB8_1");
    DDIP_level3[3][8] = (TEveGeoShapeExtract*)DDIP_level2[3]->GetElements()->FindObject("SLB0_1");
    
    // from SC05I_1
    DDIP_level4[0][0][0] = (TEveGeoShapeExtract*)DDIP_level3[0][0]->GetElements()->FindObject("SC5I_5");
    DDIP_level4[0][0][1] = (TEveGeoShapeExtract*)DDIP_level3[0][0]->GetElements()->FindObject("S05I_6");
    DDIP_level4[0][0][2] = (TEveGeoShapeExtract*)DDIP_level3[0][0]->GetElements()->FindObject("S05I_7");
    DDIP_level4[0][0][3] = (TEveGeoShapeExtract*)DDIP_level3[0][0]->GetElements()->FindObject("SB5I_8");
    
    DDIP_level4[0][1][0] = (TEveGeoShapeExtract*)DDIP_level3[0][1]->GetElements()->FindObject("SD5I_17");
    DDIP_level4[0][1][1] = (TEveGeoShapeExtract*)DDIP_level3[0][1]->GetElements()->FindObject("S05I_18");
    DDIP_level4[0][1][2] = (TEveGeoShapeExtract*)DDIP_level3[0][1]->GetElements()->FindObject("S05I_19");
    DDIP_level4[0][1][3] = (TEveGeoShapeExtract*)DDIP_level3[0][1]->GetElements()->FindObject("SB5I_20");
    
    DDIP_level4[0][2][0] = (TEveGeoShapeExtract*)DDIP_level3[0][2]->GetElements()->FindObject("SD5I_21");
    DDIP_level4[0][2][1] = (TEveGeoShapeExtract*)DDIP_level3[0][2]->GetElements()->FindObject("S05I_22");
    DDIP_level4[0][2][2] = (TEveGeoShapeExtract*)DDIP_level3[0][2]->GetElements()->FindObject("S05I_23");
    DDIP_level4[0][2][3] = (TEveGeoShapeExtract*)DDIP_level3[0][2]->GetElements()->FindObject("SB5I_24");
    
    DDIP_level4[0][3][0] = (TEveGeoShapeExtract*)DDIP_level3[0][3]->GetElements()->FindObject("S05I_33");
    DDIP_level4[0][3][1] = (TEveGeoShapeExtract*)DDIP_level3[0][3]->GetElements()->FindObject("S05I_34");
    DDIP_level4[0][3][2] = (TEveGeoShapeExtract*)DDIP_level3[0][3]->GetElements()->FindObject("S05I_35");
    DDIP_level4[0][3][3] = (TEveGeoShapeExtract*)DDIP_level3[0][3]->GetElements()->FindObject("SB5I_36");
    
    DDIP_level4[0][4][0] = (TEveGeoShapeExtract*)DDIP_level3[0][4]->GetElements()->FindObject("S05I_37");
    DDIP_level4[0][4][1] = (TEveGeoShapeExtract*)DDIP_level3[0][4]->GetElements()->FindObject("S05I_38");
    DDIP_level4[0][4][2] = (TEveGeoShapeExtract*)DDIP_level3[0][4]->GetElements()->FindObject("S05I_39");
    DDIP_level4[0][4][3] = (TEveGeoShapeExtract*)DDIP_level3[0][4]->GetElements()->FindObject("SB5I_40");
    
    DDIP_level4[0][5][0] = (TEveGeoShapeExtract*)DDIP_level3[0][5]->GetElements()->FindObject("S05I_47");
    DDIP_level4[0][5][1] = (TEveGeoShapeExtract*)DDIP_level3[0][5]->GetElements()->FindObject("S05I_48");
    DDIP_level4[0][5][2] = (TEveGeoShapeExtract*)DDIP_level3[0][5]->GetElements()->FindObject("S05I_49");
    DDIP_level4[0][5][3] = NULL;
    
    DDIP_level4[0][6][0] = (TEveGeoShapeExtract*)DDIP_level3[0][6]->GetElements()->FindObject("S05I_50");
    DDIP_level4[0][6][1] = (TEveGeoShapeExtract*)DDIP_level3[0][6]->GetElements()->FindObject("S05I_51");
    DDIP_level4[0][6][2] = (TEveGeoShapeExtract*)DDIP_level3[0][6]->GetElements()->FindObject("S05I_52");
    DDIP_level4[0][6][3] = NULL;
    
    DDIP_level4[0][7][0] = (TEveGeoShapeExtract*)DDIP_level3[0][7]->GetElements()->FindObject("S05I_57");
    DDIP_level4[0][7][1] = (TEveGeoShapeExtract*)DDIP_level3[0][7]->GetElements()->FindObject("S05I_58");
    DDIP_level4[0][7][2] = NULL;
    DDIP_level4[0][7][3] = NULL;
    
    DDIP_level4[0][8][0] = (TEveGeoShapeExtract*)DDIP_level3[0][8]->GetElements()->FindObject("S05I_59");
    DDIP_level4[0][8][1] = (TEveGeoShapeExtract*)DDIP_level3[0][8]->GetElements()->FindObject("S05I_60");
    DDIP_level4[0][8][2] = NULL;
    DDIP_level4[0][8][3] = NULL;
    
    // from SC05O_1
    DDIP_level4[1][0][0] = (TEveGeoShapeExtract*)DDIP_level3[1][0]->GetElements()->FindObject("SC5I_1");
    DDIP_level4[1][0][1] = (TEveGeoShapeExtract*)DDIP_level3[1][0]->GetElements()->FindObject("S05I_2");
    DDIP_level4[1][0][2] = (TEveGeoShapeExtract*)DDIP_level3[1][0]->GetElements()->FindObject("S05I_3");
    DDIP_level4[1][0][3] = (TEveGeoShapeExtract*)DDIP_level3[1][0]->GetElements()->FindObject("SB5I_4");
    
    DDIP_level4[1][1][0] = (TEveGeoShapeExtract*)DDIP_level3[1][1]->GetElements()->FindObject("SD5I_13");
    DDIP_level4[1][1][1] = (TEveGeoShapeExtract*)DDIP_level3[1][1]->GetElements()->FindObject("S05I_14");
    DDIP_level4[1][1][2] = (TEveGeoShapeExtract*)DDIP_level3[1][1]->GetElements()->FindObject("S05I_15");
    DDIP_level4[1][1][3] = (TEveGeoShapeExtract*)DDIP_level3[1][1]->GetElements()->FindObject("SB5I_16");
    
    DDIP_level4[1][2][0] = (TEveGeoShapeExtract*)DDIP_level3[1][2]->GetElements()->FindObject("SD5I_9");
    DDIP_level4[1][2][1] = (TEveGeoShapeExtract*)DDIP_level3[1][2]->GetElements()->FindObject("S05I_10");
    DDIP_level4[1][2][2] = (TEveGeoShapeExtract*)DDIP_level3[1][2]->GetElements()->FindObject("S05I_11");
    DDIP_level4[1][2][3] = (TEveGeoShapeExtract*)DDIP_level3[1][2]->GetElements()->FindObject("SB5I_12");
    
    DDIP_level4[1][3][0] = (TEveGeoShapeExtract*)DDIP_level3[1][3]->GetElements()->FindObject("S05I_29");
    DDIP_level4[1][3][1] = (TEveGeoShapeExtract*)DDIP_level3[1][3]->GetElements()->FindObject("S05I_30");
    DDIP_level4[1][3][2] = (TEveGeoShapeExtract*)DDIP_level3[1][3]->GetElements()->FindObject("S05I_31");
    DDIP_level4[1][3][3] = (TEveGeoShapeExtract*)DDIP_level3[1][3]->GetElements()->FindObject("SB5I_32");
    
    DDIP_level4[1][4][0] = (TEveGeoShapeExtract*)DDIP_level3[1][4]->GetElements()->FindObject("S05I_25");
    DDIP_level4[1][4][1] = (TEveGeoShapeExtract*)DDIP_level3[1][4]->GetElements()->FindObject("S05I_26");
    DDIP_level4[1][4][2] = (TEveGeoShapeExtract*)DDIP_level3[1][4]->GetElements()->FindObject("S05I_27");
    DDIP_level4[1][4][3] = (TEveGeoShapeExtract*)DDIP_level3[1][4]->GetElements()->FindObject("SB5I_28");
    
    DDIP_level4[1][5][0] = (TEveGeoShapeExtract*)DDIP_level3[1][5]->GetElements()->FindObject("S05I_44");
    DDIP_level4[1][5][1] = (TEveGeoShapeExtract*)DDIP_level3[1][5]->GetElements()->FindObject("S05I_45");
    DDIP_level4[1][5][2] = (TEveGeoShapeExtract*)DDIP_level3[1][5]->GetElements()->FindObject("S05I_46");
    DDIP_level4[1][5][3] = NULL;
    
    DDIP_level4[1][6][0] = (TEveGeoShapeExtract*)DDIP_level3[1][6]->GetElements()->FindObject("S05I_41");
    DDIP_level4[1][6][1] = (TEveGeoShapeExtract*)DDIP_level3[1][6]->GetElements()->FindObject("S05I_42");
    DDIP_level4[1][6][2] = (TEveGeoShapeExtract*)DDIP_level3[1][6]->GetElements()->FindObject("S05I_43");
    DDIP_level4[1][6][3] = NULL;
    
    DDIP_level4[1][7][0] = (TEveGeoShapeExtract*)DDIP_level3[1][7]->GetElements()->FindObject("S05I_55");
    DDIP_level4[1][7][1] = (TEveGeoShapeExtract*)DDIP_level3[1][7]->GetElements()->FindObject("S05I_56");
    DDIP_level4[1][7][2] = NULL;
    DDIP_level4[1][7][3] = NULL;
    
    DDIP_level4[1][8][0] = (TEveGeoShapeExtract*)DDIP_level3[1][8]->GetElements()->FindObject("S05I_53");
    DDIP_level4[1][8][1] = (TEveGeoShapeExtract*)DDIP_level3[1][8]->GetElements()->FindObject("S05I_54");
    DDIP_level4[1][8][2] = NULL;
    DDIP_level4[1][8][3] = NULL;
    
    // from SC06I_1
    DDIP_level4[2][0][0] = (TEveGeoShapeExtract*)DDIP_level3[2][0]->GetElements()->FindObject("SC6I_5");
    DDIP_level4[2][0][1] = (TEveGeoShapeExtract*)DDIP_level3[2][0]->GetElements()->FindObject("S06I_6");
    DDIP_level4[2][0][2] = (TEveGeoShapeExtract*)DDIP_level3[2][0]->GetElements()->FindObject("S06I_7");
    DDIP_level4[2][0][3] = (TEveGeoShapeExtract*)DDIP_level3[2][0]->GetElements()->FindObject("S06I_8");
    
    DDIP_level4[2][1][0] = (TEveGeoShapeExtract*)DDIP_level3[2][1]->GetElements()->FindObject("SD6I_17");
    DDIP_level4[2][1][1] = (TEveGeoShapeExtract*)DDIP_level3[2][1]->GetElements()->FindObject("S06I_18");
    DDIP_level4[2][1][2] = (TEveGeoShapeExtract*)DDIP_level3[2][1]->GetElements()->FindObject("S06I_19");
    DDIP_level4[2][1][3] = (TEveGeoShapeExtract*)DDIP_level3[2][1]->GetElements()->FindObject("S06I_20");
    
    DDIP_level4[2][2][0] = (TEveGeoShapeExtract*)DDIP_level3[2][2]->GetElements()->FindObject("SD6I_21");
    DDIP_level4[2][2][1] = (TEveGeoShapeExtract*)DDIP_level3[2][2]->GetElements()->FindObject("S06I_22");
    DDIP_level4[2][2][2] = (TEveGeoShapeExtract*)DDIP_level3[2][2]->GetElements()->FindObject("S06I_23");
    DDIP_level4[2][2][3] = (TEveGeoShapeExtract*)DDIP_level3[2][2]->GetElements()->FindObject("S06I_24");
    
    DDIP_level4[2][3][0] = (TEveGeoShapeExtract*)DDIP_level3[2][3]->GetElements()->FindObject("S06I_33");
    DDIP_level4[2][3][1] = (TEveGeoShapeExtract*)DDIP_level3[2][3]->GetElements()->FindObject("S06I_34");
    DDIP_level4[2][3][2] = (TEveGeoShapeExtract*)DDIP_level3[2][3]->GetElements()->FindObject("S06I_35");
    DDIP_level4[2][3][3] = (TEveGeoShapeExtract*)DDIP_level3[2][3]->GetElements()->FindObject("S06I_36");
    
    DDIP_level4[2][4][0] = (TEveGeoShapeExtract*)DDIP_level3[2][4]->GetElements()->FindObject("S06I_37");
    DDIP_level4[2][4][1] = (TEveGeoShapeExtract*)DDIP_level3[2][4]->GetElements()->FindObject("S06I_38");
    DDIP_level4[2][4][2] = (TEveGeoShapeExtract*)DDIP_level3[2][4]->GetElements()->FindObject("S06I_39");
    DDIP_level4[2][4][3] = (TEveGeoShapeExtract*)DDIP_level3[2][4]->GetElements()->FindObject("S06I_40");
    
    DDIP_level4[2][5][0] = (TEveGeoShapeExtract*)DDIP_level3[2][5]->GetElements()->FindObject("S06I_47");
    DDIP_level4[2][5][1] = (TEveGeoShapeExtract*)DDIP_level3[2][5]->GetElements()->FindObject("S06I_48");
    DDIP_level4[2][5][2] = (TEveGeoShapeExtract*)DDIP_level3[2][5]->GetElements()->FindObject("S06I_49");
    DDIP_level4[2][5][3] = NULL;
    
    DDIP_level4[2][6][0] = (TEveGeoShapeExtract*)DDIP_level3[2][6]->GetElements()->FindObject("S06I_50");
    DDIP_level4[2][6][1] = (TEveGeoShapeExtract*)DDIP_level3[2][6]->GetElements()->FindObject("S06I_51");
    DDIP_level4[2][6][2] = (TEveGeoShapeExtract*)DDIP_level3[2][6]->GetElements()->FindObject("S06I_52");
    DDIP_level4[2][6][3] = NULL;
    
    DDIP_level4[2][7][0] = (TEveGeoShapeExtract*)DDIP_level3[2][7]->GetElements()->FindObject("S06I_57");
    DDIP_level4[2][7][1] = (TEveGeoShapeExtract*)DDIP_level3[2][7]->GetElements()->FindObject("S06I_58");
    DDIP_level4[2][7][2] = NULL;
    DDIP_level4[2][7][3] = NULL;
    
    DDIP_level4[2][8][0] = (TEveGeoShapeExtract*)DDIP_level3[2][8]->GetElements()->FindObject("S06I_59");
    DDIP_level4[2][8][1] = (TEveGeoShapeExtract*)DDIP_level3[2][8]->GetElements()->FindObject("S06I_60");
    DDIP_level4[2][8][2] = NULL;
    DDIP_level4[2][8][3] = NULL;
    
    // from SC06O_1
    DDIP_level4[3][0][0] = (TEveGeoShapeExtract*)DDIP_level3[3][0]->GetElements()->FindObject("SC6I_1");
    DDIP_level4[3][0][1] = (TEveGeoShapeExtract*)DDIP_level3[3][0]->GetElements()->FindObject("S06I_2");
    DDIP_level4[3][0][2] = (TEveGeoShapeExtract*)DDIP_level3[3][0]->GetElements()->FindObject("S06I_3");
    DDIP_level4[3][0][3] = (TEveGeoShapeExtract*)DDIP_level3[3][0]->GetElements()->FindObject("S06I_4");
    
    DDIP_level4[3][1][0] = (TEveGeoShapeExtract*)DDIP_level3[3][1]->GetElements()->FindObject("SD6I_13");
    DDIP_level4[3][1][1] = (TEveGeoShapeExtract*)DDIP_level3[3][1]->GetElements()->FindObject("S06I_14");
    DDIP_level4[3][1][2] = (TEveGeoShapeExtract*)DDIP_level3[3][1]->GetElements()->FindObject("S06I_15");
    DDIP_level4[3][1][3] = (TEveGeoShapeExtract*)DDIP_level3[3][1]->GetElements()->FindObject("S06I_16");
    
    DDIP_level4[3][2][0] = (TEveGeoShapeExtract*)DDIP_level3[3][2]->GetElements()->FindObject("SD6I_9");
    DDIP_level4[3][2][1] = (TEveGeoShapeExtract*)DDIP_level3[3][2]->GetElements()->FindObject("S06I_10");
    DDIP_level4[3][2][2] = (TEveGeoShapeExtract*)DDIP_level3[3][2]->GetElements()->FindObject("S06I_11");
    DDIP_level4[3][2][3] = (TEveGeoShapeExtract*)DDIP_level3[3][2]->GetElements()->FindObject("S06I_12");
    
    DDIP_level4[3][3][0] = (TEveGeoShapeExtract*)DDIP_level3[3][3]->GetElements()->FindObject("S06I_29");
    DDIP_level4[3][3][1] = (TEveGeoShapeExtract*)DDIP_level3[3][3]->GetElements()->FindObject("S06I_30");
    DDIP_level4[3][3][2] = (TEveGeoShapeExtract*)DDIP_level3[3][3]->GetElements()->FindObject("S06I_31");
    DDIP_level4[3][3][3] = (TEveGeoShapeExtract*)DDIP_level3[3][3]->GetElements()->FindObject("S06I_32");
    
    DDIP_level4[3][4][0] = (TEveGeoShapeExtract*)DDIP_level3[3][4]->GetElements()->FindObject("S06I_25");
    DDIP_level4[3][4][1] = (TEveGeoShapeExtract*)DDIP_level3[3][4]->GetElements()->FindObject("S06I_26");
    DDIP_level4[3][4][2] = (TEveGeoShapeExtract*)DDIP_level3[3][4]->GetElements()->FindObject("S06I_27");
    DDIP_level4[3][4][3] = (TEveGeoShapeExtract*)DDIP_level3[3][4]->GetElements()->FindObject("S06I_28");
    
    DDIP_level4[3][5][0] = (TEveGeoShapeExtract*)DDIP_level3[3][5]->GetElements()->FindObject("S06I_44");
    DDIP_level4[3][5][1] = (TEveGeoShapeExtract*)DDIP_level3[3][5]->GetElements()->FindObject("S06I_45");
    DDIP_level4[3][5][2] = (TEveGeoShapeExtract*)DDIP_level3[3][5]->GetElements()->FindObject("S06I_46");
    DDIP_level4[3][5][3] = NULL;
    
    DDIP_level4[3][6][0] = (TEveGeoShapeExtract*)DDIP_level3[3][6]->GetElements()->FindObject("S06I_41");
    DDIP_level4[3][6][1] = (TEveGeoShapeExtract*)DDIP_level3[3][6]->GetElements()->FindObject("S06I_42");
    DDIP_level4[3][6][2] = (TEveGeoShapeExtract*)DDIP_level3[3][6]->GetElements()->FindObject("S06I_43");
    DDIP_level4[3][6][3] = NULL;
    
    DDIP_level4[3][7][0] = (TEveGeoShapeExtract*)DDIP_level3[3][7]->GetElements()->FindObject("S06I_55");
    DDIP_level4[3][7][1] = (TEveGeoShapeExtract*)DDIP_level3[3][7]->GetElements()->FindObject("S06I_56");
    DDIP_level4[3][7][2] = NULL;
    DDIP_level4[3][7][3] = NULL;
    
    DDIP_level4[3][8][0] = (TEveGeoShapeExtract*)DDIP_level3[3][8]->GetElements()->FindObject("S06I_53");
    DDIP_level4[3][8][1] = (TEveGeoShapeExtract*)DDIP_level3[3][8]->GetElements()->FindObject("S06I_54");
    DDIP_level4[3][8][2] = NULL;
    DDIP_level4[3][8][3] = NULL;
    
    
    DDIP_eveShapes[0] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][0][0]->GetElements()->FindObject("SC5P_1"))->GetElements()->FindObject("SC5H_1"))->GetElements()->FindObject("SC5G_1"))->GetShape();
    DDIP_eveShapes[1] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][0][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[2] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][0][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[3] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][0][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_eveShapes[4] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][1][0]->GetElements()->FindObject("SD5P_1"))->GetElements()->FindObject("SD5H_1"))->GetElements()->FindObject("SD5G_1"))->GetShape();
    DDIP_eveShapes[5] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][1][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[6] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][1][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[7] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][1][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_eveShapes[8] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][2][0]->GetElements()->FindObject("SD5P_1"))->GetElements()->FindObject("SD5H_1"))->GetElements()->FindObject("SD5G_1"))->GetShape();
    DDIP_bbShapes[0]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][2][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[1]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][2][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[2]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][2][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[3]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][3][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[4]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][3][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[5]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][3][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[6]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][3][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[7]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][4][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[8]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][4][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[9]  =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][4][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[10] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][4][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[11] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][5][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[12] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][5][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[13] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][5][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[14] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][6][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[15] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][6][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[16] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][6][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[17] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][7][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[18] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][7][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[19] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][8][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[20] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[0][8][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    
    esp[0] =DDIP_level4[0][0][0];
    esp[1] =DDIP_level4[0][0][1];
    esp[2] =DDIP_level4[0][0][2];
    esp[3] =DDIP_level4[0][0][3];
    
    esp[4] =DDIP_level4[0][1][0];
    esp[5] =DDIP_level4[0][1][1];
    esp[6] =DDIP_level4[0][1][2];
    esp[7] =DDIP_level4[0][1][3];
    
    esp[8] =DDIP_level4[0][2][0];
    bbp[0] =DDIP_level4[0][2][1];
    bbp[1] =DDIP_level4[0][2][2];
    bbp[2] =DDIP_level4[0][2][3];
    
    bbp[3] =DDIP_level4[0][3][0];
    bbp[4] =DDIP_level4[0][3][1];
    bbp[5] =DDIP_level4[0][3][2];
    bbp[6] =DDIP_level4[0][3][3];
    
    bbp[7] =DDIP_level4[0][4][0];
    bbp[8] =DDIP_level4[0][4][1];
    bbp[9] =DDIP_level4[0][4][2];
    bbp[10]=DDIP_level4[0][4][3];
    
    bbp[11]=DDIP_level4[0][5][0];
    bbp[12]=DDIP_level4[0][5][1];
    bbp[13]=DDIP_level4[0][5][2];
    
    bbp[14]=DDIP_level4[0][6][0];
    bbp[15]=DDIP_level4[0][6][1];
    bbp[16]=DDIP_level4[0][6][2];
    
    bbp[17]=DDIP_level4[0][7][0];
    bbp[18]=DDIP_level4[0][7][1];
    
    bbp[19]=DDIP_level4[0][8][0];
    bbp[20]=DDIP_level4[0][8][1];
    
    
    DDIP_eveShapes[9] = (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][0][0]->GetElements()->FindObject("SC5P_1"))->GetElements()->FindObject("SC5H_1"))->GetElements()->FindObject("SC5G_1"))->GetShape();
    DDIP_eveShapes[10]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][0][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[11]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][0][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[12]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][0][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_eveShapes[13]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][1][0]->GetElements()->FindObject("SD5P_1"))->GetElements()->FindObject("SD5H_1"))->GetElements()->FindObject("SD5G_1"))->GetShape();
    DDIP_eveShapes[14]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][1][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[15]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][1][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_eveShapes[16]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][1][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_eveShapes[17]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][2][0]->GetElements()->FindObject("SD5P_1"))->GetElements()->FindObject("SD5H_1"))->GetElements()->FindObject("SD5G_1"))->GetShape();
    DDIP_bbShapes[21] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][2][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[22] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][2][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[23] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][2][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[24] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][3][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[25] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][3][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[26] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][3][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[27] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][3][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[28] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][4][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[29] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][4][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[30] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][4][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[31] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][4][3]->GetElements()->FindObject("SB5P_1"))->GetElements()->FindObject("SB5H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[32] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][5][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[33] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][5][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[34] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][5][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[35] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][6][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[36] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][6][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[37] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][6][2]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[38] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][7][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[39] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][7][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    DDIP_bbShapes[40] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][8][0]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    DDIP_bbShapes[41] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[1][8][1]->GetElements()->FindObject("S05P_1"))->GetElements()->FindObject("S05H_1"))->GetElements()->FindObject("S05G_1"))->GetShape();
    
    esp[9] =DDIP_level4[1][0][0];
    esp[10]=DDIP_level4[1][0][1];
    esp[11]=DDIP_level4[1][0][2];
    esp[12]=DDIP_level4[1][0][3];
    
    esp[13]=DDIP_level4[1][1][0];
    esp[14]=DDIP_level4[1][1][1];
    esp[15]=DDIP_level4[1][1][2];
    esp[16]=DDIP_level4[1][1][3];
    
    esp[17]=DDIP_level4[1][2][0];
    
    bbp[21]=DDIP_level4[1][2][1];
    bbp[22]=DDIP_level4[1][2][2];
    bbp[23]=DDIP_level4[1][2][3];
    bbp[24]=DDIP_level4[1][3][0];
    bbp[25]=DDIP_level4[1][3][1];
    bbp[26]=DDIP_level4[1][3][2];
    bbp[27]=DDIP_level4[1][3][3];
    bbp[28]=DDIP_level4[1][4][0];
    bbp[29]=DDIP_level4[1][4][1];
    bbp[30]=DDIP_level4[1][4][2];
    bbp[31]=DDIP_level4[1][4][3];
    bbp[32]=DDIP_level4[1][5][0];
    bbp[33]=DDIP_level4[1][5][1];
    bbp[34]=DDIP_level4[1][5][2];
    
    bbp[35]=DDIP_level4[1][6][0];
    bbp[36]=DDIP_level4[1][6][1];
    bbp[37]=DDIP_level4[1][6][2];
    
    bbp[38]=DDIP_level4[1][7][0];
    bbp[39]=DDIP_level4[1][7][1];
    
    bbp[40]=DDIP_level4[1][8][0];
    bbp[41]=DDIP_level4[1][8][1];
    
    // SC06I_1
    DDIP_eveShapes[18]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][0][0]->GetElements()->FindObject("SC6P_1"))->GetElements()->FindObject("SC6H_1"))->GetElements()->FindObject("SC6G_1"))->GetShape();
    DDIP_eveShapes[19]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][0][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[20]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][0][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[21]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][0][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_eveShapes[22]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][1][0]->GetElements()->FindObject("SD6P_1"))->GetElements()->FindObject("SD6H_1"))->GetElements()->FindObject("SD6G_1"))->GetShape();
    DDIP_eveShapes[23]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][1][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[24]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][1][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[25]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][1][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_eveShapes[26]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][2][0]->GetElements()->FindObject("SD6P_1"))->GetElements()->FindObject("SD6H_1"))->GetElements()->FindObject("SD6G_1"))->GetShape();
    DDIP_bbShapes[42] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][2][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[43] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][2][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[44] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][2][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[45] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][3][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[46] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][3][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[47] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][3][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[48] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][3][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[49] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][4][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[50] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][4][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[51] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][4][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[52] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][4][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[53] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][5][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[54] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][5][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[55] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][5][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[56] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][6][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[57] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][6][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[58] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][6][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[59] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][7][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[60] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][7][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[61] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][8][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[62] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[2][8][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    
    esp[18]=DDIP_level4[2][0][0];
    esp[19]=DDIP_level4[2][0][1];
    esp[20]=DDIP_level4[2][0][2];
    esp[21]=DDIP_level4[2][0][3];
    
    esp[22]=DDIP_level4[2][1][0];
    esp[23]=DDIP_level4[2][1][1];
    esp[24]=DDIP_level4[2][1][2];
    esp[25]=DDIP_level4[2][1][3];
    
    esp[26]=DDIP_level4[2][2][0];
    bbp[42]=DDIP_level4[2][2][1];
    bbp[43]=DDIP_level4[2][2][2];
    bbp[44]=DDIP_level4[2][2][3];
    
    bbp[45]=DDIP_level4[2][3][0];
    bbp[46]=DDIP_level4[2][3][1];
    bbp[47]=DDIP_level4[2][3][2];
    bbp[48]=DDIP_level4[2][3][3];
    
    bbp[49]=DDIP_level4[2][4][0];
    bbp[50]=DDIP_level4[2][4][1];
    bbp[51]=DDIP_level4[2][4][2];
    bbp[52]=DDIP_level4[2][4][3];
    
    bbp[53]=DDIP_level4[2][5][0];
    bbp[54]=DDIP_level4[2][5][1];
    bbp[55]=DDIP_level4[2][5][2];
    
    bbp[56]=DDIP_level4[2][6][0];
    bbp[57]=DDIP_level4[2][6][1];
    bbp[58]=DDIP_level4[2][6][2];
    
    bbp[59]=DDIP_level4[2][7][0];
    bbp[60]=DDIP_level4[2][7][1];
    
    bbp[61]=DDIP_level4[2][8][0];
    bbp[62]=DDIP_level4[2][8][1];
    
    
    DDIP_eveShapes[27]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][0][0]->GetElements()->FindObject("SC6P_1"))->GetElements()->FindObject("SC6H_1"))->GetElements()->FindObject("SC6G_1"))->GetShape();
    DDIP_eveShapes[28]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][0][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[29]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][0][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[30]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][0][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_eveShapes[31]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][1][0]->GetElements()->FindObject("SD6P_1"))->GetElements()->FindObject("SD6H_1"))->GetElements()->FindObject("SD6G_1"))->GetShape();
    DDIP_eveShapes[32]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][1][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[33]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][1][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_eveShapes[34]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][1][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_eveShapes[35]= (TEveGeoPolyShape*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][2][0]->GetElements()->FindObject("SD6P_1"))->GetElements()->FindObject("SD6H_1"))->GetElements()->FindObject("SD6G_1"))->GetShape();
    DDIP_bbShapes[63] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][2][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[64] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][2][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[65] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][2][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[66] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][3][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[67] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][3][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[68] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][3][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[69] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][3][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[70] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][4][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[71] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][4][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[72] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][4][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[73] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][4][3]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[74] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][5][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[75] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][5][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[76] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][5][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[77] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][6][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[78] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][6][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[79] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][6][2]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[80] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][7][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[81] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][7][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    DDIP_bbShapes[82] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][8][0]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    DDIP_bbShapes[83] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)DDIP_level4[3][8][1]->GetElements()->FindObject("S06P_1"))->GetElements()->FindObject("S06H_1"))->GetElements()->FindObject("S06G_1"))->GetShape();
    
    esp[27]=DDIP_level4[3][0][0];
    esp[28]=DDIP_level4[3][0][1];
    esp[29]=DDIP_level4[3][0][2];
    esp[30]=DDIP_level4[3][0][3];
    
    esp[31]=DDIP_level4[3][1][0];
    esp[32]=DDIP_level4[3][1][1];
    esp[33]=DDIP_level4[3][1][2];
    esp[34]=DDIP_level4[3][1][3];
    
    esp[35]=DDIP_level4[3][2][0];
    
    bbp[63]=DDIP_level4[3][2][1];
    bbp[64]=DDIP_level4[3][2][2];
    bbp[65]=DDIP_level4[3][2][3];
    bbp[66]=DDIP_level4[3][3][0];
    bbp[67]=DDIP_level4[3][3][1];
    bbp[68]=DDIP_level4[3][3][2];
    bbp[69]=DDIP_level4[3][3][3];
    bbp[70]=DDIP_level4[3][4][0];
    bbp[71]=DDIP_level4[3][4][1];
    bbp[72]=DDIP_level4[3][4][2];
    bbp[73]=DDIP_level4[3][4][3];
    bbp[74]=DDIP_level4[3][5][0];
    bbp[75]=DDIP_level4[3][5][1];
    bbp[76]=DDIP_level4[3][5][2];
    
    bbp[77]=DDIP_level4[3][6][0];
    bbp[78]=DDIP_level4[3][6][1];
    bbp[79]=DDIP_level4[3][6][2];
    
    bbp[80]=DDIP_level4[3][7][0];
    bbp[81]=DDIP_level4[3][7][1];
    
    bbp[82]=DDIP_level4[3][8][0];
    bbp[83]=DDIP_level4[3][8][1];
    
    //    for(int i=42;i<84;i++)
    //    {
    //        double bb_hl[3] = {DDIP_bbShapes[i]->GetDX(),DDIP_bbShapes[i]->GetDY(),DDIP_bbShapes[i]->GetDZ()};// BB half-length in X,Y,Z
    //        double *bb_origin = DDIP_bbShapes[i]->GetOrigin();// BB origin
    //
    //        cout<<"\n\nmuon 2 bbox "<<i<<" params (box):"<<endl;
    //        cout<<"BB lenght ( "<<2*bb_hl[0]<<" , "<<2*bb_hl[1]<<" , "<<2*bb_hl[2]<<" )";
    //        cout<<" origin ("<<bb_origin[0]<<" , "<<bb_origin[1]<<" , "<<bb_origin[2]<<" )"<<endl;
    //        printTrans(bbp[i]->GetTrans());
    //        cout<<"\n\n"<<endl;
    //    }
    
    for(int i=17;i<36;i++)
    {
        double bb_hl[3] = {DDIP_eveShapes[i]->GetDX(),DDIP_eveShapes[i]->GetDY(),DDIP_eveShapes[i]->GetDZ()};// BB half-length in X,Y,Z
        double *bb_origin = DDIP_eveShapes[i]->GetOrigin();// BB origin
        
        cout<<"\n\nmuon 2 polygon "<<i<<" params (box):"<<endl;
        cout<<"BB lenght ( "<<2*bb_hl[0]<<" , "<<2*bb_hl[1]<<" , "<<2*bb_hl[2]<<" )";
        cout<<" origin ("<<bb_origin[0]<<" , "<<bb_origin[1]<<" , "<<bb_origin[2]<<" )"<<endl;
        printTrans(esp[i]->GetTrans());
        cout<<"\n\n"<<endl;
        
        TString path = TString::Format("%smuonEveShape_%i.xml",outPath.Data(),i);
        DDIP_eveShapes[i]->SaveAs(path.Data());
        printMuon(path);
    }
}



void muon3()
{
    //----------------------------------------------------------------------------------//
    //                                                                                  //
    //      WARNING!!!                                                                  //
    //                                                                                  //
    //      If muon2 may cause you cry, this function may actually drive you crazy.     //
    //      You better don't look inside.                                               //
    //                                                                                  //
    //----------------------------------------------------------------------------------//
    
    TFile *muonGeomFile = TFile::Open(Form("%s/../src/EVE/resources/geometry/gentle_geo_muon.root",gSystem->Getenv("ALICE_ROOT")));
    TEveGeoShapeExtract* gsem = (TEveGeoShapeExtract*) muonGeomFile->Get("Gentle MUON");
    TList *muonList = gsem->GetElements(); // list with all sub-detectors in file
    
    TEveGeoShapeExtract *level1 = (TEveGeoShapeExtract*)muonList->FindObject("YOUT2_1");
    
    // for third element of level1 (YOUT2_1)
    TEveGeoShapeExtract *YOUT2_level2[12];
    TEveGeoShapeExtract *YOUT2_level3[12][18];
    TEveGeoShapeExtract *YOUT2_level4[12][18][6];
    TEveGeoPolyShape *YOUT2_eveShapes[100];
    TGeoBBox *YOUT2_bbShapes[1000];
    TEveGeoShapeExtract *bbp[1000];
    TEveGeoShapeExtract *esp[100];
    
    for(int i=0;i<12;i++){
        YOUT2_level2[i]=NULL;
        for(int j=0;j<16;j++){
            YOUT2_level3[i][j]=NULL;
            for(int k=0;k<6;k++){
                YOUT2_level4[i][j][k]=NULL;
            }
        }
    }
    
    YOUT2_level2[0] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC07I_1");
    YOUT2_level2[1] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC07O_1");
    YOUT2_level2[2] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC08I_1");
    YOUT2_level2[3] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC08O_1");
    YOUT2_level2[4] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC09I_1");
    YOUT2_level2[5] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC09O_1");
    YOUT2_level2[6] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC10I_1");
    YOUT2_level2[7] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC10O_1");
    YOUT2_level2[8] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC11_1");
    YOUT2_level2[9] = (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC12_1");
    YOUT2_level2[10]= (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC13_1");
    YOUT2_level2[11]= (TEveGeoShapeExtract*)level1->GetElements()->FindObject("SC14_1");
    
    
    YOUT2_level3[0][0] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC19_1");
    YOUT2_level3[0][1] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC20_1");
    YOUT2_level3[0][2] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC18_1");
    YOUT2_level3[0][3] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC21_1");
    YOUT2_level3[0][4] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC17_1");
    YOUT2_level3[0][5] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC22_1");
    YOUT2_level3[0][6] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC16_1");
    YOUT2_level3[0][7] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC23_1");
    YOUT2_level3[0][8] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC15_1");
    YOUT2_level3[0][9] = (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC24_1");
    YOUT2_level3[0][10]= (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC14_1");
    YOUT2_level3[0][11]= (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC25_1");
    YOUT2_level3[0][12]= (TEveGeoShapeExtract*)YOUT2_level2[0]->GetElements()->FindObject("SLC13_1");
    
    YOUT2_level3[1][0] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC6_1");
    YOUT2_level3[1][1] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC7_1");
    YOUT2_level3[1][2] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC5_1");
    YOUT2_level3[1][3] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC8_1");
    YOUT2_level3[1][4] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC4_1");
    YOUT2_level3[1][5] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC9_1");
    YOUT2_level3[1][6] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC3_1");
    YOUT2_level3[1][7] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC10_1");
    YOUT2_level3[1][8] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC2_1");
    YOUT2_level3[1][9] = (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC11_1");
    YOUT2_level3[1][10]= (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC1_1");
    YOUT2_level3[1][11]= (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC12_1");
    YOUT2_level3[1][12]= (TEveGeoShapeExtract*)YOUT2_level2[1]->GetElements()->FindObject("SLC0_1");
    
    YOUT2_level3[2][0] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD19_1");
    YOUT2_level3[2][1] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD20_1");
    YOUT2_level3[2][2] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD18_1");
    YOUT2_level3[2][3] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD21_1");
    YOUT2_level3[2][4] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD17_1");
    YOUT2_level3[2][5] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD22_1");
    YOUT2_level3[2][6] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD16_1");
    YOUT2_level3[2][7] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD23_1");
    YOUT2_level3[2][8] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD15_1");
    YOUT2_level3[2][9] = (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD24_1");
    YOUT2_level3[2][10]= (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD14_1");
    YOUT2_level3[2][11]= (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD25_1");
    YOUT2_level3[2][12]= (TEveGeoShapeExtract*)YOUT2_level2[2]->GetElements()->FindObject("SLD13_1");
    
    YOUT2_level3[3][0] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD6_1");
    YOUT2_level3[3][1] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD7_1");
    YOUT2_level3[3][2] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD5_1");
    YOUT2_level3[3][3] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD8_1");
    YOUT2_level3[3][4] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD4_1");
    YOUT2_level3[3][5] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD9_1");
    YOUT2_level3[3][6] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD3_1");
    YOUT2_level3[3][7] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD10_1");
    YOUT2_level3[3][8] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD2_1");
    YOUT2_level3[3][9] = (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD11_1");
    YOUT2_level3[3][10]= (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD1_1");
    YOUT2_level3[3][11]= (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD12_1");
    YOUT2_level3[3][12]= (TEveGeoShapeExtract*)YOUT2_level2[3]->GetElements()->FindObject("SLD0_1");

    YOUT2_level3[4][0] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE19_1");
    YOUT2_level3[4][1] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE20_1");
    YOUT2_level3[4][2] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE18_1");
    YOUT2_level3[4][3] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE21_1");
    YOUT2_level3[4][4] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE17_1");
    YOUT2_level3[4][5] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE22_1");
    YOUT2_level3[4][6] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE16_1");
    YOUT2_level3[4][7] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE23_1");
    YOUT2_level3[4][8] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE15_1");
    YOUT2_level3[4][9] = (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE24_1");
    YOUT2_level3[4][10]= (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE14_1");
    YOUT2_level3[4][11]= (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE25_1");
    YOUT2_level3[4][12]= (TEveGeoShapeExtract*)YOUT2_level2[4]->GetElements()->FindObject("SLE13_1");
    
    YOUT2_level3[5][0] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE6_1");
    YOUT2_level3[5][1] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE7_1");
    YOUT2_level3[5][2] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE5_1");
    YOUT2_level3[5][3] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE8_1");
    YOUT2_level3[5][4] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE4_1");
    YOUT2_level3[5][5] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE9_1");
    YOUT2_level3[5][6] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE3_1");
    YOUT2_level3[5][7] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE10_1");
    YOUT2_level3[5][8] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE2_1");
    YOUT2_level3[5][9] = (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE11_1");
    YOUT2_level3[5][10]= (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE1_1");
    YOUT2_level3[5][11]= (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE12_1");
    YOUT2_level3[5][12]= (TEveGeoShapeExtract*)YOUT2_level2[5]->GetElements()->FindObject("SLE0_1");

    YOUT2_level3[6][0] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF19_1");
    YOUT2_level3[6][1] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF20_1");
    YOUT2_level3[6][2] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF18_1");
    YOUT2_level3[6][3] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF21_1");
    YOUT2_level3[6][4] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF17_1");
    YOUT2_level3[6][5] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF22_1");
    YOUT2_level3[6][6] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF16_1");
    YOUT2_level3[6][7] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF23_1");
    YOUT2_level3[6][8] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF15_1");
    YOUT2_level3[6][9] = (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF24_1");
    YOUT2_level3[6][10]= (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF14_1");
    YOUT2_level3[6][11]= (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF25_1");
    YOUT2_level3[6][12]= (TEveGeoShapeExtract*)YOUT2_level2[6]->GetElements()->FindObject("SLF13_1");
    
    YOUT2_level3[7][0] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF6_1");
    YOUT2_level3[7][1] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF7_1");
    YOUT2_level3[7][2] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF5_1");
    YOUT2_level3[7][3] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF8_1");
    YOUT2_level3[7][4] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF4_1");
    YOUT2_level3[7][5] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF9_1");
    YOUT2_level3[7][6] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF3_1");
    YOUT2_level3[7][7] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF10_1");
    YOUT2_level3[7][8] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF2_1");
    YOUT2_level3[7][9] = (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF11_1");
    YOUT2_level3[7][10]= (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF1_1");
    YOUT2_level3[7][11]= (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF12_1");
    YOUT2_level3[7][12]= (TEveGeoShapeExtract*)YOUT2_level2[7]->GetElements()->FindObject("SLF0_1");

    YOUT2_level3[8][0] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R5_1");
    YOUT2_level3[8][1] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L5_1");
    YOUT2_level3[8][2] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R4_1");
    YOUT2_level3[8][3] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R6_1");
    YOUT2_level3[8][4] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L4_1");
    YOUT2_level3[8][5] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L6_1");
    YOUT2_level3[8][6] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R3_1");
    YOUT2_level3[8][7] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R7_1");
    YOUT2_level3[8][8] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L3_1");
    YOUT2_level3[8][9] = (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L7_1");
    YOUT2_level3[8][10]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R2_1");
    YOUT2_level3[8][11]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R8_1");
    YOUT2_level3[8][12]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L2_1");
    YOUT2_level3[8][13]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L8_1");
    YOUT2_level3[8][14]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R1_1");
    YOUT2_level3[8][15]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0R9_1");
    YOUT2_level3[8][16]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L1_1");
    YOUT2_level3[8][17]= (TEveGeoShapeExtract*)YOUT2_level2[8]->GetElements()->FindObject("S0L9_1");

    YOUT2_level3[9][0] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R5_1");
    YOUT2_level3[9][1] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L5_1");
    YOUT2_level3[9][2] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R4_1");
    YOUT2_level3[9][3] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R6_1");
    YOUT2_level3[9][4] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L4_1");
    YOUT2_level3[9][5] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L6_1");
    YOUT2_level3[9][6] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R3_1");
    YOUT2_level3[9][7] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R7_1");
    YOUT2_level3[9][8] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L3_1");
    YOUT2_level3[9][9] = (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L7_1");
    YOUT2_level3[9][10]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R2_1");
    YOUT2_level3[9][11]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R8_1");
    YOUT2_level3[9][12]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L2_1");
    YOUT2_level3[9][13]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L8_1");
    YOUT2_level3[9][14]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R1_1");
    YOUT2_level3[9][15]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1R9_1");
    YOUT2_level3[9][16]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L1_1");
    YOUT2_level3[9][17]= (TEveGeoShapeExtract*)YOUT2_level2[9]->GetElements()->FindObject("S1L9_1");
    
    YOUT2_level3[10][0] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R5_1");
    YOUT2_level3[10][1] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L5_1");
    YOUT2_level3[10][2] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R4_1");
    YOUT2_level3[10][3] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R6_1");
    YOUT2_level3[10][4] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L4_1");
    YOUT2_level3[10][5] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L6_1");
    YOUT2_level3[10][6] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R3_1");
    YOUT2_level3[10][7] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R7_1");
    YOUT2_level3[10][8] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L3_1");
    YOUT2_level3[10][9] = (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L7_1");
    YOUT2_level3[10][10]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R2_1");
    YOUT2_level3[10][11]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R8_1");
    YOUT2_level3[10][12]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L2_1");
    YOUT2_level3[10][13]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L8_1");
    YOUT2_level3[10][14]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R1_1");
    YOUT2_level3[10][15]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2R9_1");
    YOUT2_level3[10][16]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L1_1");
    YOUT2_level3[10][17]= (TEveGeoShapeExtract*)YOUT2_level2[10]->GetElements()->FindObject("S2L9_1");
    
    YOUT2_level3[11][0] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R5_1");
    YOUT2_level3[11][1] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L5_1");
    YOUT2_level3[11][2] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R4_1");
    YOUT2_level3[11][3] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R6_1");
    YOUT2_level3[11][4] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L4_1");
    YOUT2_level3[11][5] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L6_1");
    YOUT2_level3[11][6] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R3_1");
    YOUT2_level3[11][7] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R7_1");
    YOUT2_level3[11][8] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L3_1");
    YOUT2_level3[11][9] = (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L7_1");
    YOUT2_level3[11][10]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R2_1");
    YOUT2_level3[11][11]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R8_1");
    YOUT2_level3[11][12]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L2_1");
    YOUT2_level3[11][13]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L8_1");
    YOUT2_level3[11][14]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R1_1");
    YOUT2_level3[11][15]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3R9_1");
    YOUT2_level3[11][16]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L1_1");
    YOUT2_level3[11][17]= (TEveGeoShapeExtract*)YOUT2_level2[11]->GetElements()->FindObject("S3L9_1");
    
    
    
    
    
    // from SC07I_1
    for(int i=0;i<2;i++)
    {
        YOUT2_level4[2*i][0][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(Form("S0%iI_6",7+i));
        YOUT2_level4[2*i][0][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(Form("S0%iI_7",7+i));
        YOUT2_level4[2*i][0][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(Form("S0%iI_8",7+i));
        YOUT2_level4[2*i][0][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(Form("S0%iI_9",7+i));
        YOUT2_level4[2*i][0][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(Form("S0%iI_10",7+i));

        YOUT2_level4[2*i][1][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("SD%iI_23",7+i));
        YOUT2_level4[2*i][1][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("S0%iI_24",7+i));
        YOUT2_level4[2*i][1][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("S0%iI_25",7+i));
        YOUT2_level4[2*i][1][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("S0%iI_26",7+i));
        YOUT2_level4[2*i][1][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("S0%iI_27",7+i));
        YOUT2_level4[2*i][1][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(Form("S0%iI_28",7+i));
        
        YOUT2_level4[2*i][2][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("SD%iI_29",7+i));
        YOUT2_level4[2*i][2][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("S0%iI_30",7+i));
        YOUT2_level4[2*i][2][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("S0%iI_31",7+i));
        YOUT2_level4[2*i][2][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("S0%iI_32",7+i));
        YOUT2_level4[2*i][2][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("S0%iI_33",7+i));
        YOUT2_level4[2*i][2][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(Form("S0%iI_34",7+i));
        
        YOUT2_level4[2*i][3][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(Form("S0%iI_45",7+i));
        YOUT2_level4[2*i][3][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(Form("S0%iI_46",7+i));
        YOUT2_level4[2*i][3][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(Form("S0%iI_47",7+i));
        YOUT2_level4[2*i][3][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(Form("S0%iI_48",7+i));
        YOUT2_level4[2*i][3][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(Form("S0%iI_49",7+i));
        
        YOUT2_level4[2*i][4][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(Form("S0%iI_50",7+i));
        YOUT2_level4[2*i][4][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(Form("S0%iI_51",7+i));
        YOUT2_level4[2*i][4][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(Form("S0%iI_52",7+i));
        YOUT2_level4[2*i][4][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(Form("S0%iI_53",7+i));
        YOUT2_level4[2*i][4][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(Form("S0%iI_54",7+i));
        
        YOUT2_level4[2*i][5][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(Form("S0%iI_65",7+i));
        YOUT2_level4[2*i][5][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(Form("S0%iI_66",7+i));
        YOUT2_level4[2*i][5][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(Form("S0%iI_67",7+i));
        YOUT2_level4[2*i][5][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(Form("S0%iI_68",7+i));
        YOUT2_level4[2*i][5][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(Form("S0%iI_69",7+i));

        YOUT2_level4[2*i][6][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(Form("S0%iI_70",7+i));
        YOUT2_level4[2*i][6][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(Form("S0%iI_71",7+i));
        YOUT2_level4[2*i][6][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(Form("S0%iI_72",7+i));
        YOUT2_level4[2*i][6][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(Form("S0%iI_73",7+i));
        YOUT2_level4[2*i][6][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(Form("S0%iI_74",7+i));

        YOUT2_level4[2*i][7][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(Form("S0%iI_83",7+i));
        YOUT2_level4[2*i][7][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(Form("S0%iI_84",7+i));
        YOUT2_level4[2*i][7][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(Form("S0%iI_85",7+i));
        YOUT2_level4[2*i][7][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(Form("S0%iI_86",7+i));

        YOUT2_level4[2*i][8][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(Form("S0%iI_87",7+i));
        YOUT2_level4[2*i][8][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(Form("S0%iI_88",7+i));
        YOUT2_level4[2*i][8][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(Form("S0%iI_89",7+i));
        YOUT2_level4[2*i][8][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(Form("S0%iI_90",7+i));

        YOUT2_level4[2*i][9][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(Form("S0%iI_97",7+i));
        YOUT2_level4[2*i][9][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(Form("S0%iI_98",7+i));
        YOUT2_level4[2*i][9][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(Form("S0%iI_99",7+i));

        YOUT2_level4[2*i][10][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(Form("S0%iI_100",7+i));
        YOUT2_level4[2*i][10][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(Form("S0%iI_101",7+i));
        YOUT2_level4[2*i][10][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(Form("S0%iI_102",7+i));

        YOUT2_level4[2*i][11][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][11]->GetElements()->FindObject(Form("S0%iI_107",7+i));
        YOUT2_level4[2*i][11][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][11]->GetElements()->FindObject(Form("S0%iI_108",7+i));

        YOUT2_level4[2*i][12][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][12]->GetElements()->FindObject(Form("S0%iI_109",7+i));
        YOUT2_level4[2*i][12][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][12]->GetElements()->FindObject(Form("S0%iI_110",7+i));
        
        //from SC07O_1
        YOUT2_level4[2*i+1][0][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(Form("S0%iI_1",7+i));
        YOUT2_level4[2*i+1][0][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(Form("S0%iI_2",7+i));
        YOUT2_level4[2*i+1][0][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(Form("S0%iI_3",7+i));
        YOUT2_level4[2*i+1][0][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(Form("S0%iI_4",7+i));
        YOUT2_level4[2*i+1][0][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(Form("S0%iI_5",7+i));
        
        YOUT2_level4[2*i+1][1][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("SD%iI_17",7+i));
        YOUT2_level4[2*i+1][1][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("S0%iI_18",7+i));
        YOUT2_level4[2*i+1][1][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("S0%iI_19",7+i));
        YOUT2_level4[2*i+1][1][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("S0%iI_20",7+i));
        YOUT2_level4[2*i+1][1][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("S0%iI_21",7+i));
        YOUT2_level4[2*i+1][1][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(Form("S0%iI_22",7+i));

        YOUT2_level4[2*i+1][2][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("SD%iI_11",7+i));
        YOUT2_level4[2*i+1][2][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("S0%iI_12",7+i));
        YOUT2_level4[2*i+1][2][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("S0%iI_13",7+i));
        YOUT2_level4[2*i+1][2][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("S0%iI_14",7+i));
        YOUT2_level4[2*i+1][2][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("S0%iI_15",7+i));
        YOUT2_level4[2*i+1][2][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(Form("S0%iI_16",7+i));
        
        YOUT2_level4[2*i+1][3][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(Form("S0%iI_40",7+i));
        YOUT2_level4[2*i+1][3][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(Form("S0%iI_41",7+i));
        YOUT2_level4[2*i+1][3][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(Form("S0%iI_42",7+i));
        YOUT2_level4[2*i+1][3][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(Form("S0%iI_43",7+i));
        YOUT2_level4[2*i+1][3][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(Form("S0%iI_44",7+i));

        YOUT2_level4[2*i+1][4][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(Form("S0%iI_35",7+i));
        YOUT2_level4[2*i+1][4][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(Form("S0%iI_36",7+i));
        YOUT2_level4[2*i+1][4][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(Form("S0%iI_37",7+i));
        YOUT2_level4[2*i+1][4][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(Form("S0%iI_38",7+i));
        YOUT2_level4[2*i+1][4][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(Form("S0%iI_39",7+i));
        
        YOUT2_level4[2*i+1][5][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(Form("S0%iI_60",7+i));
        YOUT2_level4[2*i+1][5][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(Form("S0%iI_61",7+i));
        YOUT2_level4[2*i+1][5][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(Form("S0%iI_62",7+i));
        YOUT2_level4[2*i+1][5][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(Form("S0%iI_63",7+i));
        YOUT2_level4[2*i+1][5][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(Form("S0%iI_64",7+i));
        
        YOUT2_level4[2*i+1][6][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(Form("S0%iI_55",7+i));
        YOUT2_level4[2*i+1][6][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(Form("S0%iI_56",7+i));
        YOUT2_level4[2*i+1][6][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(Form("S0%iI_57",7+i));
        YOUT2_level4[2*i+1][6][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(Form("S0%iI_58",7+i));
        YOUT2_level4[2*i+1][6][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(Form("S0%iI_59",7+i));
        
        YOUT2_level4[2*i+1][7][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(Form("S0%iI_79",7+i));
        YOUT2_level4[2*i+1][7][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(Form("S0%iI_80",7+i));
        YOUT2_level4[2*i+1][7][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(Form("S0%iI_81",7+i));
        YOUT2_level4[2*i+1][7][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(Form("S0%iI_82",7+i));
        
        YOUT2_level4[2*i+1][8][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(Form("S0%iI_75",7+i));
        YOUT2_level4[2*i+1][8][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(Form("S0%iI_76",7+i));
        YOUT2_level4[2*i+1][8][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(Form("S0%iI_77",7+i));
        YOUT2_level4[2*i+1][8][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(Form("S0%iI_78",7+i));
        
        YOUT2_level4[2*i+1][9][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(Form("S0%iI_94",7+i));
        YOUT2_level4[2*i+1][9][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(Form("S0%iI_95",7+i));
        YOUT2_level4[2*i+1][9][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(Form("S0%iI_96",7+i));
        
        YOUT2_level4[2*i+1][10][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(Form("S0%iI_91",7+i));
        YOUT2_level4[2*i+1][10][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(Form("S0%iI_92",7+i));
        YOUT2_level4[2*i+1][10][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(Form("S0%iI_93",7+i));
        
        YOUT2_level4[2*i+1][11][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][11]->GetElements()->FindObject(Form("S0%iI_105",7+i));
        YOUT2_level4[2*i+1][11][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][11]->GetElements()->FindObject(Form("S0%iI_106",7+i));
        
        YOUT2_level4[2*i+1][12][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][12]->GetElements()->FindObject(Form("S0%iI_103",7+i));
        YOUT2_level4[2*i+1][12][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][12]->GetElements()->FindObject(Form("S0%iI_104",7+i));
    }
    for(int i=2;i<4;i++)
    {
        YOUT2_level4[2*i][0][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(i==2 ? "S09I_6"  : "S10I_6"));
        YOUT2_level4[2*i][0][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(i==2 ? "S09I_7"  : "S10I_7"));
        YOUT2_level4[2*i][0][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(i==2 ? "S09I_8"  : "S10I_8"));
        YOUT2_level4[2*i][0][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(i==2 ? "S09I_9"  : "S10I_9"));
        YOUT2_level4[2*i][0][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][0]->GetElements()->FindObject(i==2 ? "S09I_10" : "S10I_10"));
        
        YOUT2_level4[2*i][1][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "SD9I_23" : "SD0I_23"));
        YOUT2_level4[2*i][1][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "S09I_24" : "S10I_24"));
        YOUT2_level4[2*i][1][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "S09I_25" : "S10I_25"));
        YOUT2_level4[2*i][1][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "S09I_26" : "S10I_26"));
        YOUT2_level4[2*i][1][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "S09I_27" : "S10I_27"));
        YOUT2_level4[2*i][1][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][1]->GetElements()->FindObject(i==2 ? "S09I_28" : "S10I_28"));
        
        YOUT2_level4[2*i][2][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "SD9I_29" : "SD0I_29"));
        YOUT2_level4[2*i][2][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "S09I_30" : "S10I_30"));
        YOUT2_level4[2*i][2][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "S09I_31" : "S10I_31"));
        YOUT2_level4[2*i][2][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "S09I_32" : "S10I_32"));
        YOUT2_level4[2*i][2][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "S09I_33" : "S10I_33"));
        YOUT2_level4[2*i][2][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][2]->GetElements()->FindObject(i==2 ? "S09I_34" : "S10I_34"));
        
        YOUT2_level4[2*i][3][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_47" : "S10I_47"));
        YOUT2_level4[2*i][3][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_48" : "S10I_48"));
        YOUT2_level4[2*i][3][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_49" : "S10I_49"));
        YOUT2_level4[2*i][3][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_50" : "S10I_50"));
        YOUT2_level4[2*i][3][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_51" : "S10I_51"));
        YOUT2_level4[2*i][3][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][3]->GetElements()->FindObject(i==2 ? "S09I_52" : "S10I_52"));
        
        YOUT2_level4[2*i][4][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_53" : "S10I_53"));
        YOUT2_level4[2*i][4][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_54" : "S10I_54"));
        YOUT2_level4[2*i][4][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_55" : "S10I_55"));
        YOUT2_level4[2*i][4][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_56" : "S10I_56"));
        YOUT2_level4[2*i][4][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_57" : "S10I_57"));
        YOUT2_level4[2*i][4][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][4]->GetElements()->FindObject(i==2 ? "S09I_58" : "S10I_58"));
        
        YOUT2_level4[2*i][5][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_71" : "S10I_71"));
        YOUT2_level4[2*i][5][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_72" : "S10I_72"));
        YOUT2_level4[2*i][5][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_73" : "S10I_73"));
        YOUT2_level4[2*i][5][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_74" : "S10I_74"));
        YOUT2_level4[2*i][5][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_75" : "S10I_75"));
        YOUT2_level4[2*i][5][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][5]->GetElements()->FindObject(i==2 ? "S09I_76" : "S10I_76"));
        
        YOUT2_level4[2*i][6][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_77" : "S10I_77"));
        YOUT2_level4[2*i][6][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_78" : "S10I_78"));
        YOUT2_level4[2*i][6][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_79" : "S10I_79"));
        YOUT2_level4[2*i][6][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_80" : "S10I_80"));
        YOUT2_level4[2*i][6][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_81" : "S10I_81"));
        YOUT2_level4[2*i][6][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][6]->GetElements()->FindObject(i==2 ? "S09I_82" : "S10I_82"));
        
        YOUT2_level4[2*i][7][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(i==2 ? "S09I_93" : "S10I_93"));
        YOUT2_level4[2*i][7][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(i==2 ? "S09I_94" : "S10I_94"));
        YOUT2_level4[2*i][7][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(i==2 ? "S09I_95" : "S10I_95"));
        YOUT2_level4[2*i][7][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(i==2 ? "S09I_96" : "S10I_96"));
        YOUT2_level4[2*i][7][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][7]->GetElements()->FindObject(i==2 ? "S09I_97" : "S10I_97"));
        
        YOUT2_level4[2*i][8][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(i==2 ? "S09I_98" : "S10I_98"));
        YOUT2_level4[2*i][8][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(i==2 ? "S09I_99" : "S10I_99"));
        YOUT2_level4[2*i][8][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(i==2 ? "S09I_100" : "S10I_100"));
        YOUT2_level4[2*i][8][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(i==2 ? "S09I_101" : "S10I_101"));
        YOUT2_level4[2*i][8][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][8]->GetElements()->FindObject(i==2 ? "S09I_102" : "S10I_102"));
        
        YOUT2_level4[2*i][9][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(i==2 ? "S09I_111" : "S10I_111"));
        YOUT2_level4[2*i][9][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(i==2 ? "S09I_112" : "S10I_112"));
        YOUT2_level4[2*i][9][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(i==2 ? "S09I_113" : "S10I_113"));
        YOUT2_level4[2*i][9][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i][9]->GetElements()->FindObject(i==2 ? "S09I_114" : "S10I_114"));
        
        YOUT2_level4[2*i][10][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(i==2 ? "S09I_115" : "S10I_115"));
        YOUT2_level4[2*i][10][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(i==2 ? "S09I_116" : "S10I_116"));
        YOUT2_level4[2*i][10][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(i==2 ? "S09I_117" : "S10I_117"));
        YOUT2_level4[2*i][10][3]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][10]->GetElements()->FindObject(i==2 ? "S09I_118" : "S10I_118"));
        
        YOUT2_level4[2*i][11][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][11]->GetElements()->FindObject(i==2 ? "S09I_125" : "S10I_125"));
        YOUT2_level4[2*i][11][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][11]->GetElements()->FindObject(i==2 ? "S09I_126" : "S10I_126"));
        YOUT2_level4[2*i][11][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][11]->GetElements()->FindObject(i==2 ? "S09I_127" : "S10I_127"));
        
        YOUT2_level4[2*i][12][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][12]->GetElements()->FindObject(i==2 ? "S09I_128" : "S10I_128"));
        YOUT2_level4[2*i][12][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][12]->GetElements()->FindObject(i==2 ? "S09I_129" : "S10I_129"));
        YOUT2_level4[2*i][12][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i][12]->GetElements()->FindObject(i==2 ? "S09I_130" : "S10I_130"));
        
        //from SC07O_1
        YOUT2_level4[2*i+1][0][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(i==2 ? "S09I_1" : "S10I_1"));
        YOUT2_level4[2*i+1][0][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(i==2 ? "S09I_2" : "S10I_2"));
        YOUT2_level4[2*i+1][0][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(i==2 ? "S09I_3" : "S10I_3"));
        YOUT2_level4[2*i+1][0][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(i==2 ? "S09I_4" : "S10I_4"));
        YOUT2_level4[2*i+1][0][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][0]->GetElements()->FindObject(i==2 ? "S09I_5" : "S10I_5"));
        
        YOUT2_level4[2*i+1][1][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "SD9I_17" : "SD0I_17"));
        YOUT2_level4[2*i+1][1][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "S09I_18" : "S10I_18"));
        YOUT2_level4[2*i+1][1][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "S09I_19" : "S10I_19"));
        YOUT2_level4[2*i+1][1][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "S09I_20" : "S10I_20"));
        YOUT2_level4[2*i+1][1][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "S09I_21" : "S10I_21"));
        YOUT2_level4[2*i+1][1][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][1]->GetElements()->FindObject(i==2 ? "S09I_22" : "S10I_22"));
        
        YOUT2_level4[2*i+1][2][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "SD9I_11" : "SD0I_11"));
        YOUT2_level4[2*i+1][2][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "S09I_12" : "S10I_12"));
        YOUT2_level4[2*i+1][2][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "S09I_13" : "S10I_13"));
        YOUT2_level4[2*i+1][2][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "S09I_14" : "S10I_14"));
        YOUT2_level4[2*i+1][2][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "S09I_15" : "S10I_15"));
        YOUT2_level4[2*i+1][2][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][2]->GetElements()->FindObject(i==2 ? "S09I_16" : "S10I_16"));
        
        YOUT2_level4[2*i+1][3][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_41" : "S10I_41"));
        YOUT2_level4[2*i+1][3][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_42" : "S10I_42"));
        YOUT2_level4[2*i+1][3][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_43" : "S10I_43"));
        YOUT2_level4[2*i+1][3][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_44" : "S10I_44"));
        YOUT2_level4[2*i+1][3][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_45" : "S10I_45"));
        YOUT2_level4[2*i+1][3][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][3]->GetElements()->FindObject(i==2 ? "S09I_46" : "S10I_46"));
        
        YOUT2_level4[2*i+1][4][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_35" : "S10I_35"));
        YOUT2_level4[2*i+1][4][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_36" : "S10I_36"));
        YOUT2_level4[2*i+1][4][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_37" : "S10I_37"));
        YOUT2_level4[2*i+1][4][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_38" : "S10I_38"));
        YOUT2_level4[2*i+1][4][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_39" : "S10I_39"));
        YOUT2_level4[2*i+1][4][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][4]->GetElements()->FindObject(i==2 ? "S09I_40" : "S10I_40"));
        
        YOUT2_level4[2*i+1][5][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_65" : "S10I_65"));
        YOUT2_level4[2*i+1][5][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_66" : "S10I_66"));
        YOUT2_level4[2*i+1][5][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_67" : "S10I_67"));
        YOUT2_level4[2*i+1][5][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_68" : "S10I_68"));
        YOUT2_level4[2*i+1][5][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_69" : "S10I_69"));
        YOUT2_level4[2*i+1][5][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][5]->GetElements()->FindObject(i==2 ? "S09I_70" : "S10I_70"));

        YOUT2_level4[2*i+1][6][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_59" : "S10I_59"));
        YOUT2_level4[2*i+1][6][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_60" : "S10I_60"));
        YOUT2_level4[2*i+1][6][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_61" : "S10I_61"));
        YOUT2_level4[2*i+1][6][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_62" : "S10I_62"));
        YOUT2_level4[2*i+1][6][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_63" : "S10I_63"));
        YOUT2_level4[2*i+1][6][5] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][6]->GetElements()->FindObject(i==2 ? "S09I_64" : "S10I_64"));
        
        YOUT2_level4[2*i+1][7][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(i==2 ? "S09I_88" : "S10I_88"));
        YOUT2_level4[2*i+1][7][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(i==2 ? "S09I_89" : "S10I_89"));
        YOUT2_level4[2*i+1][7][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(i==2 ? "S09I_90" : "S10I_90"));
        YOUT2_level4[2*i+1][7][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(i==2 ? "S09I_91" : "S10I_91"));
        YOUT2_level4[2*i+1][7][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][7]->GetElements()->FindObject(i==2 ? "S09I_92" : "S10I_92"));
        
        YOUT2_level4[2*i+1][8][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(i==2 ? "S09I_83" : "S10I_83"));
        YOUT2_level4[2*i+1][8][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(i==2 ? "S09I_84" : "S10I_84"));
        YOUT2_level4[2*i+1][8][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(i==2 ? "S09I_85" : "S10I_85"));
        YOUT2_level4[2*i+1][8][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(i==2 ? "S09I_86" : "S10I_86"));
        YOUT2_level4[2*i+1][8][4] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][8]->GetElements()->FindObject(i==2 ? "S09I_87" : "S10I_87"));
        
        YOUT2_level4[2*i+1][9][0] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(i==2 ? "S09I_107" : "S10I_107"));
        YOUT2_level4[2*i+1][9][1] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(i==2 ? "S09I_108" : "S10I_108"));
        YOUT2_level4[2*i+1][9][2] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(i==2 ? "S09I_109" : "S10I_109"));
        YOUT2_level4[2*i+1][9][3] = (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][9]->GetElements()->FindObject(i==2 ? "S09I_110" : "S10I_110"));

        
        YOUT2_level4[2*i+1][10][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(i==2 ? "S09I_103" : "S10I_103"));
        YOUT2_level4[2*i+1][10][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(i==2 ? "S09I_104" : "S10I_104"));
        YOUT2_level4[2*i+1][10][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(i==2 ? "S09I_105" : "S10I_105"));
        YOUT2_level4[2*i+1][10][3]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][10]->GetElements()->FindObject(i==2 ? "S09I_106" : "S10I_106"));
        
        YOUT2_level4[2*i+1][11][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][11]->GetElements()->FindObject(i==2 ? "S09I_122" : "S10I_122"));
        YOUT2_level4[2*i+1][11][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][11]->GetElements()->FindObject(i==2 ? "S09I_123" : "S10I_123"));
        YOUT2_level4[2*i+1][11][2]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][11]->GetElements()->FindObject(i==2 ? "S09I_124" : "S10I_124"));
        
        YOUT2_level4[2*i+1][12][0]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][12]->GetElements()->FindObject(i==2 ? "S09I_119" : "S10I_119"));
        YOUT2_level4[2*i+1][12][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][12]->GetElements()->FindObject(i==2 ? "S09I_120" : "S10I_120"));
        YOUT2_level4[2*i+1][12][1]= (TEveGeoShapeExtract*)YOUT2_level3[2*i+1][12]->GetElements()->FindObject(i==2 ? "S09I_121" : "S10I_121"));
    }
    for(int i=8;i<12;i++)
    {
        YOUT2_level4[i][0][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][0]->GetElements()->FindObject(Form("SC%iA_9",i-7));
        
        YOUT2_level4[i][1][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][1]->GetElements()->FindObject(Form("SC%iA_10",i-7));
        
        YOUT2_level4[i][2][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][2]->GetElements()->FindObject(Form("SC%iA_41",i-7));
        YOUT2_level4[i][2][1] = (TEveGeoShapeExtract*)YOUT2_level3[i][2]->GetElements()->FindObject(Form("SC%iA_45",i-7));
        
        YOUT2_level4[i][3][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][3]->GetElements()->FindObject(Form("SC%iA_42",i-7));
        YOUT2_level4[i][3][1] = (TEveGeoShapeExtract*)YOUT2_level3[i][3]->GetElements()->FindObject(Form("SC%iA_46",i-7));
        
        YOUT2_level4[i][4][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][4]->GetElements()->FindObject(Form("SC%iA_43",i-7));
        YOUT2_level4[i][4][1] = (TEveGeoShapeExtract*)YOUT2_level3[i][4]->GetElements()->FindObject(Form("SC%iA_47",i-7));
        
        YOUT2_level4[i][5][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][5]->GetElements()->FindObject(Form("SC%iA_44",i-7));
        YOUT2_level4[i][5][1] = (TEveGeoShapeExtract*)YOUT2_level3[i][5]->GetElements()->FindObject(Form("SC%iA_48",i-7));
        
        YOUT2_level4[i][6][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][6]->GetElements()->FindObject(Form("SC%iA_109",i-7));
        YOUT2_level4[i][7][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][7]->GetElements()->FindObject(Form("SC%iA_110",i-7));
        YOUT2_level4[i][8][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][8]->GetElements()->FindObject(Form("SC%iA_111",i-7));
        YOUT2_level4[i][9][0] = (TEveGeoShapeExtract*)YOUT2_level3[i][9]->GetElements()->FindObject(Form("SC%iA_112",i-7));
        YOUT2_level4[i][10][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][10]->GetElements()->FindObject(Form("SC%iA_173",i-7));
        YOUT2_level4[i][11][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][11]->GetElements()->FindObject(Form("SC%iA_174",i-7));
        YOUT2_level4[i][12][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][12]->GetElements()->FindObject(Form("SC%iA_175",i-7));
        YOUT2_level4[i][13][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][13]->GetElements()->FindObject(Form("SC%iA_176",i-7));
        YOUT2_level4[i][14][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][14]->GetElements()->FindObject(Form("SC%iA_237",i-7));
        YOUT2_level4[i][15][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][15]->GetElements()->FindObject(Form("SC%iA_238",i-7));
        YOUT2_level4[i][16][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][16]->GetElements()->FindObject(Form("SC%iA_239",i-7));
        YOUT2_level4[i][17][0]= (TEveGeoShapeExtract*)YOUT2_level3[i][17]->GetElements()->FindObject(Form("SC%iA_240",i-7));
    }
    
    
    int iter=0,iter2=0;
    for(int j=0;j<4;j++)
    {
        for(int i=0;i<2;i++)
        {
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][0][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][0][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][0][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][0][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][0][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            
            bbp[iter2++] = YOUT2_level4[2*j+i][0][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][0][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][0][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][0][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][0][4];
            

            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][0]->GetElements()->FindObject(j!=3 ? Form("SD%iP_1",7+j) : "SD0P_1"))->GetElements()->FindObject(j!=3 ? Form("SD%iH_1",7+j) : "SD0H_1"))->GetElements()->FindObject(j!=3 ? Form("SD%iG_1",7+j) : "SD0G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][1][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            
            bbp[iter2++] = YOUT2_level4[2*j+i][1][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][1][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][1][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][1][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][1][4];
            bbp[iter2++] = YOUT2_level4[2*j+i][1][5];
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][0]->GetElements()->FindObject(j!=3 ? Form("SD%iP_1",7+j) : "SD0P_1"))->GetElements()->FindObject(j!=3 ? Form("SD%iH_1",7+j) : "SD0H_1"))->GetElements()->FindObject(j!=3 ? Form("SD%iG_1",7+j) : "SD0G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][2][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            
            bbp[iter2++] = YOUT2_level4[2*j+i][2][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][2][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][2][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][2][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][2][4];
            bbp[iter2++] = YOUT2_level4[2*j+i][2][5];

            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            
            bbp[iter2++] = YOUT2_level4[2*j+i][3][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][3][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][3][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][3][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][3][4];
            
            if(YOUT2_level4[2*j+i][3][5]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][3][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][3][5];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][4][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][4][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][4][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][4][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][4][4];
            
            if(YOUT2_level4[2*j+i][4][5]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][4][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][4][5];
            }
            
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][5][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][5][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][5][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][5][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][5][4];
            
            if(YOUT2_level4[2*j+i][5][5]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][5][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][5][5];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][6][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][6][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][6][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][6][3];
            bbp[iter2++] = YOUT2_level4[2*j+i][6][4];
            
            if(YOUT2_level4[2*j+i][6][5]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][6][5]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][6][5];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][7][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][7][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][7][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][7][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][7][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][7][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][7][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][7][3];
            
            if(YOUT2_level4[2*j+i][7][4]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][7][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][7][4];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][8][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][8][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][8][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][8][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][8][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][8][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][8][2];
            bbp[iter2++] = YOUT2_level4[2*j+i][8][3];
            
            if(YOUT2_level4[2*j+i][8][4]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][8][4]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][8][4];
            }
            
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][9][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][9][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][9][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][9][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][9][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][9][2];
            
            if(YOUT2_level4[2*j+i][9][3]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][9][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][9][3];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][10][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][10][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][10][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][10][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][10][1];
            bbp[iter2++] = YOUT2_level4[2*j+i][10][2];
            
            if(YOUT2_level4[2*j+i][10][3]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][10][3]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][10][3];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][11][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][11][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][11][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][11][1];
            
            if(YOUT2_level4[2*j+i][11][2]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][11][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][11][2];
            }
            
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][12][0]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
            YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][12][1]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();

            bbp[iter2++] = YOUT2_level4[2*j+i][12][0];
            bbp[iter2++] = YOUT2_level4[2*j+i][12][1];
            
            if(YOUT2_level4[2*j+i][12][2]){
                YOUT2_bbShapes[iter++] =         (TGeoBBox*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)((TEveGeoShapeExtract*)YOUT2_level4[2*j+i][12][2]->GetElements()->FindObject(j!=3 ? Form("S0%iP_1",7+j) : "S10P_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iH_1",7+j) : "S10H_1"))->GetElements()->FindObject(j!=3 ? Form("S0%iG_1",7+j) : "S10G_1"))->GetShape();
                bbp[iter2++] = YOUT2_level4[2*j+i][12][2];
            }
        }
    }
    for(int i=8;i<12;i++)
    {
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][0][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][1][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][2][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][2][1]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][3][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][3][1]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][4][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][4][1]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][5][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][5][1]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][6][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][7][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][8][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][9][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][10][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][11][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][12][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][13][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][14][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][15][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][16][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        YOUT2_bbShapes[iter++] = (TGeoBBox*)((TEveGeoShapeExtract*)(((TEveGeoShapeExtract*)(YOUT2_level4[i][17][0]->GetElements()->FindObject(Form("SB%iA_1",i-7))))->GetElements()->FindObject(Form("S1%iG_1",i-7))))->GetShape();
        
        bbp[iter2++] = YOUT2_level4[i][0][0];
        bbp[iter2++] = YOUT2_level4[i][1][0];
        bbp[iter2++] = YOUT2_level4[i][2][0];
        bbp[iter2++] = YOUT2_level4[i][2][1];
        bbp[iter2++] = YOUT2_level4[i][3][0];
        bbp[iter2++] = YOUT2_level4[i][3][1];
        bbp[iter2++] = YOUT2_level4[i][4][0];
        bbp[iter2++] = YOUT2_level4[i][4][1];
        bbp[iter2++] = YOUT2_level4[i][5][0];
        bbp[iter2++] = YOUT2_level4[i][5][1];
        bbp[iter2++] = YOUT2_level4[i][6][0];
        bbp[iter2++] = YOUT2_level4[i][7][0];
        bbp[iter2++] = YOUT2_level4[i][8][0];
        bbp[iter2++] = YOUT2_level4[i][9][0];
        bbp[iter2++] = YOUT2_level4[i][10][0];
        bbp[iter2++] = YOUT2_level4[i][11][0];
        bbp[iter2++] = YOUT2_level4[i][12][0];
        bbp[iter2++] = YOUT2_level4[i][13][0];
        bbp[iter2++] = YOUT2_level4[i][14][0];
        bbp[iter2++] = YOUT2_level4[i][15][0];
        bbp[iter2++] = YOUT2_level4[i][16][0];
        bbp[iter2++] = YOUT2_level4[i][17][0];
    }
    

    for(int i=0;i<iter;i++)
    {
        double bb_hl[3] = {YOUT2_bbShapes[i]->GetDX(),YOUT2_bbShapes[i]->GetDY(),YOUT2_bbShapes[i]->GetDZ()};// BB half-length in X,Y,Z
        double *bb_origin = YOUT2_bbShapes[i]->GetOrigin();// BB origin
        double *trans =bbp[i]->GetTrans();
        
        if(bb_hl[0] < 0.000001)
        {
            TString path = TString::Format("%smuon3EveShape_%i.xml",outPath.Data(),i);
            (TEveGeoPolyShape*)YOUT2_bbShapes[i]->SaveAs(path.Data());
            getBBoxFromXml(path,bb_hl[0],bb_hl[1],bb_hl[2]);
        }
        
//        cout<<"\n\nmuon 3 bbox "<<i<<" params (box):"<<endl;
//        cout<<"BB lenght ( "<<2*bb_hl[0]<<" , "<<2*bb_hl[1]<<" , "<<2*bb_hl[2]<<" )";
//        cout<<" origin ("<<bb_origin[0]<<" , "<<bb_origin[1]<<" , "<<bb_origin[2]<<" )"<<endl;
//        printTrans(trans);
//        cout<<"\n\n"<<endl;
        
        boxToMayaFile(trans[12],trans[13],trans[14],bb_hl[0],bb_hl[1],bb_hl[2]);
    }
}

int muonIter=0;

void boxToMayaFile(double tx,double ty,double tz,double lx,double ly,double lz)// tx,ty,tz - translation, lx,ly,lz - distance from (0,0,0) to vertices in x,y,z
{
    testFile<<"createNode transform -n \"pCube"<<muonIter<<"\" -p \"muon7\";"<<endl;
    testFile<<"setAttr \".t\" -type \"double3\" "<<tx<<" "<<ty<<" "<<tz<<" ;"<<endl;
    testFile<<"createNode mesh -n \"pCubeShape"<<muonIter<<"\" -p \"pCube"<<muonIter<<"\";"<<endl;
    testFile<<"setAttr -s 8 \".vt[0:7]\"  "<<-lx<<" "<<-ly<<" "<<lz<<" "<<lx<<" "<<-ly<<" "<<lz<<" "<<-lx<<" "<<ly<<" "<<lz<<" "<<lx<<" "<<ly<<" "<<lz<<endl;
    testFile<<-lx<<" "<<ly<<" "<<-lz<<" "<<lx<<" "<<ly<<" "<<-lz<<" "<<-lx<<" "<<-ly<<" "<<-lz<<" "<<lx<<" "<<-ly<<" "<<-lz<<";"<<endl;
    testFile<<"setAttr -s 12 \".ed[0:11]\"  0 1 0 2 3 0 4 5 0 6 7 0 0 2 0 1 3 0 2 4 0 3 5 0 4 6 0 5 7 0 6 0 0 7 1 0;"<<endl;
    testFile<<"setAttr -s 6 -ch 24 \".fc[0:5]\" -type \"polyFaces\""<<endl;
    testFile<<"f 4 0 5 -2 -5"<<endl;
    testFile<<"mu 0 4 0 1 3 2"<<endl;
    testFile<<"f 4 1 7 -3 -7"<<endl;
    testFile<<"mu 0 4 2 3 5 4"<<endl;
    testFile<<"f 4 2 9 -4 -9"<<endl;
    testFile<<"mu 0 4 4 5 7 6"<<endl;
    testFile<<"f 4 3 11 -1 -11"<<endl;
    testFile<<"mu 0 4 6 7 9 8"<<endl;
    testFile<<"f 4 -12 -10 -8 -6"<<endl;
    testFile<<"mu 0 4 1 10 11 3"<<endl;
    testFile<<"f 4 10 4 6 8"<<endl;
    testFile<<"mu 0 4 12 0 2 13;"<<endl;
    
    muonIter++;
}

void printMuon(const char* path, bool reduce = false, double r_min=0, double r_max=10000, double z=0)
{
    cout<<path<<endl;
    ifstream inFile(path);
    vector<double> inArray;
    if (inFile.is_open())
    {
        string line;
        int from,to,from2,to2;
        while(inFile.good())
        {
            getline(inFile,line);
            if(line.find("<Double_t v=")!=string::npos)
            {
                if(line.find("cnt=")!=string::npos)
                {
                    from = line.find("<Double_t v=\"")+13;
                    to = line.find("\"\ cnt");
                    from2 = line.find("cnt=\"")+5;
                    to2 = line.find_last_of("\"");
                    for(int i=0;i<atoi(line.substr(from2,to2-from2).c_str());i++)
                    {
                        inArray.push_back(atof(line.substr(from,to-from).c_str()));
                    }
                }
                else
                {
                    from = line.find("<Double_t v=\"")+13;
                    to = line.find_last_of("\"");
                    inArray.push_back(atof(line.substr(from,to-from).c_str()));
                }
            }
        }
        inFile.close();
    }
    else{cout<<"Unable to open file"<<endl;}
    
    vector<double> x1,x2;
    vector<double> y1,y2;
    vector<double> z1,z2;
    TGraph2D *graph = new TGraph2D();
    
    for(int i=1;i<(inArray.size()+1)/3;i++)
    {
        if(inArray[i*3+2]>0){
            x1.push_back(inArray[i*3]);
            y1.push_back(inArray[i*3+1]);
            z1.push_back(inArray[i*3+2]);
        }
        else if(inArray[i*3+2]<=0)
        {
            x2.push_back(inArray[i*3]);
            y2.push_back(inArray[i*3+1]);
            z2.push_back(inArray[i*3+2]);
        }
        
        graph->SetPoint(i+1,inArray[i*3],inArray[i*3+1],inArray[i*3+2]);
    }
    
    if(reduce)
    {
        TGraph *gr = new TGraph();
        
        vector<double> x1_reduced;
        vector<double> y1_reduced;
        
        for (int i=0; i<x1.size(); i++)
        {
            if((x1[i]*x1[i]+y1[i]*y1[i]<r_min*r_min) || (x1[i]*x1[i]+y1[i]*y1[i]>r_max*r_max))
            {
                x1_reduced.push_back(x1[i]);
                y1_reduced.push_back(y1[i]);
            }
        }
        
        for(int i=0;i<x1_reduced.size();i++)
        {
            gr->SetPoint(i,x1_reduced[i],y1_reduced[i]);
            cout<<"( "<<x1_reduced[i]<<" , "<<y1_reduced[i]<<" )"<<endl;
        }
        
        //    graph->Draw("AP0");
        //    c1->cd(c1_iter++);
        //    gr->Draw("Ap*");
        cout<<"\n\nZ = "<<inArray[2]<<endl;
    }
    else
    {
        cout<<"Z>0:"<<endl;
        for (int i=0; i<x1.size(); i++) {
            cout<<"( "<<x1[i]<<" , "<<y1[i]<<" , "<<z1[i]<<" )"<<endl;
        }
        cout<<"Z<=0:"<<endl;
        for (int i=0; i<x2.size(); i++) {
            cout<<"( "<<x2[i]<<" , "<<y2[i]<<" , "<<z2[i]<<" )"<<endl;
        }
    }
}

void getBBoxFromXml(const char* path, double &hlx,double &hly,double &hlz)
{
    cout<<path<<endl;
    ifstream inFile(path);
    vector<double> inArray;
    if (inFile.is_open())
    {
        string line;
        int from,to,from2,to2;
        while(inFile.good())
        {
            getline(inFile,line);
            if(line.find("<Double_t v=")!=string::npos)
            {
                if(line.find("cnt=")!=string::npos)
                {
                    from = line.find("<Double_t v=\"")+13;
                    to = line.find("\"\ cnt");
                    from2 = line.find("cnt=\"")+5;
                    to2 = line.find_last_of("\"");
                    for(int i=0;i<atoi(line.substr(from2,to2-from2).c_str());i++)
                    {
                        inArray.push_back(atof(line.substr(from,to-from).c_str()));
                    }
                }
                else
                {
                    from = line.find("<Double_t v=\"")+13;
                    to = line.find_last_of("\"");
                    inArray.push_back(atof(line.substr(from,to-from).c_str()));
                }
            }
        }
        inFile.close();
    }
    else{cout<<"Unable to open file"<<endl;}
    
    hlx = fabs(inArray[3]);
    hly = fabs(inArray[4]);
    hlz = fabs(inArray[5]);
}



