#include <Riostream.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1F.h>
#include <TVector.h>
#include <TRefArray.h>
#include <TArrayS.h>
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TH2I.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TFile.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TBranch.h"
#include <iostream>
#include "overhead.C"
//list of all declarations
void extract(char* dc);
TFile* fin;
TFile* fout;
TList* l;
int z;
int hmtloop, particleloop, childloop, cutloop;
TH3F* s[16];
TH2D* H[8];
TH2F* M[8];
char s2[100];//character array for Lambda
char s3[100];//character array for hmt
TH1D *h[30],*b,*d,*M2;
//list of all declarations

void cinvmass4(int n,char* dc= (char*)"pikx"){
    //gROOT->LoadMacro("~/Documents/Houston/resonances/Corey/overhead.C");
    //gROOT->LoadMacro("~/Desktop/xiTest/overhead.C");
    if(!n) extract(dc);
    
    return;
}


void extract(char* dc){
    //does the signal extraction for the various decay channels using a common framework
    int j,i1,i2;
    int scheme,nm;//there are three different schemes used to do the background subtraction, which depends on the decay channel
    int nm2;
    double nmix;
    int hmtflag=0;//if hmtflag=1, look for hmt events.
    int fitflag=0;
    int cutflag=1;//if cutflag=1, begin alternate cut process.
    //int cutmax=0;
    int cutloop=0;
    int cutstart=0;
    int MCflag=0;//if MCflag=1, MC information is run
    //int nsigma=1;
    char DC2[500];//more information related to the decay channel
    char DCLambda[500];//more information related to the decay channel
    char DChmt[500];//more information related to the decay channel
    char DCcut[500];//more information related to the decay channel
    //char m0[13][100],m1[8][100];//fragments of the names of the histograms
    char m0[27][100],m1[16][100];//fragments of the names of the histograms
    
    if (cutflag == 0) {
        cout << "error cutflag must be set for now" << endl;
        return 0;
    }
    else if((MCflag == 1) && (hmtflag == 1)){
        cout << "error with MC and hmt flag" << endl;
        hmtflag=0;
    }
    
    /*if ((nsigma == 1) and (cutflag <= 1)) {
        cout << "error with nsigma cuts" << endl;
        cutstart=1;
        cutflag=2;
    }*/
    
    for (hmtloop=0; hmtloop <= hmtflag; hmtloop++) {//begin hmtloop
        j=GetDC(dc); jDC=j; TitlesDC();
        for (particleloop=0; particleloop <= 1; particleloop ++) {//begin particle loop
            if ((jDC == Klkx) || (jDC == Klk0)) {
                if (particleloop == 0) {
                    sprintf(dc,"Lambdakx");
                    cutstart=0;
                    cutflag=0;//2
                }
                else if(particleloop == 1){
                    cutstart=0;
                    cutflag=0;
                    //jDC=Klk0;
                    sprintf(dc,"Lambdak0");
                }
            }
            else if(jDC == Klpi){
                sprintf(dc,"Lambdapi");
                cutloop=0;
                particleloop=2;
            }
            
            j=GetDC(dc); jDC=j; TitlesDC();
            
            if(jDC==Kpikx || jDC==Kpkx){
                scheme=0;
                nm=6;
                sprintf(m0[0],"np");//negative-positive histogram, same event
                sprintf(m0[1],"pn");//positive-negative histogram, same event
                sprintf(m0[2],"mixnp");//negative-positive histogram, mixed events
                sprintf(m0[3],"mixpn");//positive-negative histogram, mixed events
                sprintf(m0[4],"nn");//negative-negative histogram, same event
                sprintf(m0[5],"pp");//positive-positive histogram, same event
                sprintf(m0[8],"unlike");//unlike-charge histogram, same event (sum of np and pn)
                sprintf(m0[9],"mixunlike");//unlike-charge histogram for event mixing, (sum of mixnp and mixpn)
                sprintf(m0[10],"sum");//alternate like-charge background histogram, same event (sum of nn and pp)
                sprintf(m0[11],"gm");//like-charge background histogram, same event (2*sqrt(nn*pp))
            }else if(jDC==Kpik0 || jDC==Kkxk0 || jDC==Kpk0 || jDC==Klk0){
                scheme=1;
                nm=4;
                sprintf(m0[0],"n0");//negative-neutral histogram, same event
                sprintf(m0[1],"p0");//positive-neutral histogram, same event
                sprintf(m0[2],"mixn0");//negative-neutral histogram, mixed events
                sprintf(m0[3],"mixp0");//positive-neutral histogram, mixed events
                sprintf(m0[8],"unlike");//sum of n0 and p0
                sprintf(m0[9],"mixunlike");//sum of mixn0 and mixp0
                if (MCflag == 1) {
                    sprintf(m0[14],"0p_gen");//positive/
                    sprintf(m0[15],"0a_gen");//anti
                    sprintf(m0[16],"0p_true");//positive
                    sprintf(m0[17],"0a_true");//anti
                    sprintf(m0[18],"0p_trueMM");//positive
                    sprintf(m0[19],"0a_trueMM");//anti
                    sprintf(m0[20],"0p_res");//positive
                    sprintf(m0[21],"0a_res");//anti
                    sprintf(m0[22],"sum_gen");//generated
                    sprintf(m0[23],"sum_true");//true
                    sprintf(m0[24],"sum_trueMM");//trueMM
                    sprintf(m0[25],"sum_res");//resolution
                    sprintf(m0[26],"sum_trueMM_gen");//acceptance*efficency
                }
                //NOTE: Lambda is considered "positive" and anti-Lambda is considered "negative" for this macro
            }else if(jDC==Klpi || jDC==Klkx || jDC==Klp){
                scheme=2;
                nm=8;
                sprintf(m0[0],"np");
                sprintf(m0[1],"pn");
                sprintf(m0[2],"mixnp");
                sprintf(m0[3],"mixpn");
                sprintf(m0[4],"nn");
                sprintf(m0[5],"pp");
                sprintf(m0[6],"mixnn");//negative-negative histogram, mixed events
                sprintf(m0[7],"mixpp");//positive-positive histogram, mixed events
                sprintf(m0[8],"unlike");
                sprintf(m0[9],"mixunlike");
                sprintf(m0[10],"sum");
                sprintf(m0[11],"gm");
                sprintf(m0[12],"mixlike");//like-charge histogram for event mixing (sum of mixnn and mixpp)
                sprintf(m0[13],"mixall");//(sum of mixunlike and mixlike)
                if (MCflag == 1) {
                    sprintf(m0[14],"p_gen");//positive/
                    sprintf(m0[15],"m_gen");//minus/
                    sprintf(m0[16],"p_true");//positive/
                    sprintf(m0[17],"m_true");//minus/
                    sprintf(m0[18],"p_trueMM");//positive/
                    sprintf(m0[19],"m_trueMM");//minus/
                    sprintf(m0[20],"p_res");//positive/
                    sprintf(m0[21],"m_res");//minus/
                    sprintf(m0[22],"sum_gen");//generated
                    sprintf(m0[23],"sum_true");//true
                    sprintf(m0[24],"sum_trueMM");//trueMM
                    sprintf(m0[25],"sum_res");//resolution
                    sprintf(m0[26],"sum_trueMM_gen");//acceptance*efficency
                }
            }
            
            //store the names of the histograms in the output file of the analysis task that correspond to the histograms defined above
            if(jDC==Kpikx){
                ///sprintf(DC2,"pikx_ppData_pikx");//old format
                sprintf(DC2,"pikx_TPC_pikx");
                //sprintf(m1[0],"UnlikeMP");//old format
                sprintf(m1[0],"UnlikeMP_pikx_TPC");
                sprintf(m1[1],"UnlikePM_pikx_TPC");
                sprintf(m1[2],"MixingMP_pikx_TPC");
                sprintf(m1[3],"MixingPM_pikx_TPC");
                sprintf(m1[4],"LikeMM_pikx_TPC");
                sprintf(m1[5],"LikePP_pikx_TPC");
            }else if(jDC==Kpik0){
                sprintf(DC2,"pik0_ppData_pik0");
                sprintf(m1[0],"K0Pim");
                sprintf(m1[1],"K0Pip");
                sprintf(m1[2],"K0PimMix");
                sprintf(m1[3],"K0PipMix");
                nm=4;
            }else if(jDC==Kkxk0){
                sprintf(DC2,"kxk0_ppData_kxk0");
                sprintf(m1[0],"K0Km");
                sprintf(m1[1],"K0Kp");
                sprintf(m1[2],"K0KmMix");
                sprintf(m1[3],"K0KpMix");
                nm=4;
            }else if(jDC==Kpkx){
                //sprintf(DC2,"pkx_ppData_pkx");///old format
                sprintf(DC2,"pkx_TPC_pkx");
                //sprintf(m1[0],"UnlikeMP");//old format
                sprintf(m1[0],"UnlikeMP_pkx_TPC");
                sprintf(m1[1],"UnlikePM_pkx_TPC");
                sprintf(m1[2],"MixingMP_pkx_TPC");
                sprintf(m1[3],"MixingPM_pkx_TPC");
                sprintf(m1[4],"LikeMM_pkx_TPC");
                sprintf(m1[5],"LikePP_pkx_TPC");
            }else if(jDC==Kpk0){
                sprintf(DC2,"pk0_ppData_pk0");
                sprintf(m1[0],"K0Pm");
                sprintf(m1[1],"K0Pp");
                sprintf(m1[2],"K0PmMix");
                sprintf(m1[3],"K0PpMix");
                nm=4;
            }else if(jDC==Klpi){
                sprintf(DC2,"Lambdapi_Lambdapi");
                sprintf(m1[0],"LambdaaPip_Lambdapi");
                sprintf(m1[1],"LambdapPim_Lambdapi");
                sprintf(m1[2],"LambdaaPipMix_Lambdapi");
                sprintf(m1[3],"LambdapPimMix_Lambdapi");
                sprintf(m1[4],"LambdaaPim_Lambdapi");
                sprintf(m1[5],"LambdapPip_Lambdapi");
                sprintf(m1[6],"LambdaaPimMix_Lambdapi");
                sprintf(m1[7],"LambdapPipMix_Lambdapi");
            }else if(jDC==Klkx){
                sprintf(DCLambda,"Lambdakx");
                sprintf(DChmt,"hmt");
                sprintf(DCcut,"cut1");
                ///sprintf(DC2,"Lambdakx_ppData_Lambdakx");///old format
                if (hmtloop == 1) {
                    sprintf(DC2,"Lambdakx_hmt_Lambdakx");
                    sprintf(m1[0],"LambdaaKp_Lambdakx_hmt");
                    sprintf(m1[1],"LambdapKm_Lambdakx_hmt");
                    sprintf(m1[2],"LambdaaKpMix_Lambdakx_hmt");
                    sprintf(m1[3],"LambdapKmMix_Lambdakx_hmt");
                    sprintf(m1[4],"LambdaaKm_Lambdakx_hmt");
                    sprintf(m1[5],"LambdapKp_Lambdakx_hmt");
                    sprintf(m1[6],"LambdaaKmMix_Lambdakx_hmt");
                    sprintf(m1[7],"LambdapKpMix_Lambdakx_hmt");
                }
                else
                {
                    sprintf(DC2,"Lambdakx_Lambdakx");
                    sprintf(m1[0],"LambdaaKp_Lambdakx");
                    sprintf(m1[1],"LambdapKm_Lambdakx");
                    sprintf(m1[2],"LambdaaKpMix_Lambdakx");
                    sprintf(m1[3],"LambdapKmMix_Lambdakx");
                    sprintf(m1[4],"LambdaaKm_Lambdakx");
                    sprintf(m1[5],"LambdapKp_Lambdakx");
                    sprintf(m1[6],"LambdaaKmMix_Lambdakx");
                    sprintf(m1[7],"LambdapKpMix_Lambdakx");
                }
                
            }else if(jDC==Klk0){
                sprintf(DCLambda,"Lambdak0");
                sprintf(DChmt,"hmt");
                sprintf(DCcut,"cut1");
                ///sprintf(DC2,"Lambdak0_ppData_Lambdak0");///old format
                if (hmtloop == 1) {
                    sprintf(DC2,"Lambdak0_hmt_Lambdak0");
                    //sprintf(m1[0],"LambdaaK0");///old format
                    sprintf(m1[0],"LambdaaK0_Lambdak0_hmt");
                    sprintf(m1[1],"LambdapK0_Lambdak0_hmt");
                    sprintf(m1[2],"LambdaaK0Mix_Lambdak0_hmt");
                    sprintf(m1[3],"LambdapK0Mix_Lambdak0_hmt");
                }
                else
                {
                    sprintf(DC2,"Lambdak0_Lambdak0");
                    //sprintf(m1[0],"LambdaaK0");///old format
                    sprintf(m1[0],"LambdaaK0_Lambdak0");
                    sprintf(m1[1],"LambdapK0_Lambdak0");
                    sprintf(m1[2],"LambdaaK0Mix_Lambdak0");
                    sprintf(m1[3],"LambdapK0Mix_Lambdak0");
                }
                
                nm=4;
            }else if(jDC==Klp){
                sprintf(DC2,"Lambdap_ppData_Lambdap");
                sprintf(m1[0],"LambdaaPp");
                sprintf(m1[1],"LambdapPm");
                sprintf(m1[2],"LambdaaPpMix");
                sprintf(m1[3],"LambdapPmMix");
                sprintf(m1[4],"LambdaaPm");
                sprintf(m1[5],"LambdapPp");
                sprintf(m1[6],"LambdaaPmMix");
                sprintf(m1[7],"LambdapPpMix");
            }
            
            for (childloop = 0; childloop <= 0; childloop ++) {//begin childloop//0,14 // 0,23
                /*if ((childloop >= 12) && (hmtloop == 1)) {//0,1 //12,1
                 continue;
                 }*/
                
                if (MCflag == 1) {
                    //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults994merge.root"));//output file from the analysis task
                    //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV7MC.root"));//output file from the analysis task
                    fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV9MC.root"));//output file from the analysis task
                }
                else{
                    if (hmtloop == 0) {
                            //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults992merge.root"));//output file from the analysis task
                            //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV7Normal.root"));//output file from the analysis task
                        fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV9Normal.root"));//output file from the analysis task
                        //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultspPb502.root"));//output file from the analysis task
                    }
                    else{
                        //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults993merge.root"));//output file from the analysis task
                        //fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV7hmt.root"));//output file from the analysis task
                        fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV9hmt.root"));//output file from the analysis task

                    }
                }
                /*
                 if (childloop == 0) {
                 TFile* fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV6.root"));//output file from the analysis task
                 }
                 else if((childloop >= 1) && (childloop <= 4)){//1,10
                 TFile* fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults822child%i.root",childloop));//output file from the analysis task
                 }
                 else if(childloop == 5){//11
                 TFile* fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults822merge.root"));//output file from the analysis task
                 }
                 else if((childloop >= 12) && (childloop <= 22)){
                 TFile* fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults804child%i.root",childloop-11));//output file from the analysis task
                 }
                 else if(childloop == 23){
                 TFile* fin=TFile::Open(Form("~/Desktop/xiTest/AnalysisResults804merge.root"));//output file from the analysis task
                 }*/
                
                if(!fin) return;
                
                for (cutloop=cutstart; cutloop <= cutflag; cutloop++) {//begin cutloop
                    nm2=0;
                    
                    if (MCflag == 1) {//994
                        sprintf(s0,"RsnOut_%s",DC0);//get the TList object corresponding to the desired decay channel
                        if (cutloop >= 1) {
                            sprintf(s0,"RsnOut_%s_cut%i",DC0,cutloop);//get the TList object corresponding to the desired decay channel
                        }
                    }
                    else{
                        if (hmtloop == 0) {//992
                            sprintf(s0,"RsnOut_%s",DC0);//get the TList object corresponding to the desired decay channel
                            if (cutloop >= 1) {
                                sprintf(s0,"RsnOut_%s_cut%i",DC0,cutloop);//get the TList object corresponding to the desired decay channel
                            }
                        }
                        else{//993
                            sprintf(s0,"RsnOut_%s_%s",DC0,DChmt);//get the TList object corresponding to the desired decay channel
                            if (cutloop >= 1) {
                                sprintf(s0,"RsnOut_%s_%s_cut%i",DC0,DChmt,cutloop);//get the TList object corresponding to the desired decay channel
                            }
                        }
                    }
                    
                    l=(TList*) fin->Get(s0);
                    if(!l){cerr<<"missing list "<<s0<<endl; return;}
                    cout << "list name " << s0 << endl;
                    if(jDC==Klkx){
                        if (MCflag == 1) {//994
                            sprintf(DC2,"Lambdakx_Lambdakx");
                            sprintf(m1[0],"LambdaaKp_Lambdakx");
                            sprintf(m1[1],"LambdapKm_Lambdakx");
                            sprintf(m1[2],"LambdaaKpMix_Lambdakx");
                            sprintf(m1[3],"LambdapKmMix_Lambdakx");
                            sprintf(m1[4],"LambdaaKm_Lambdakx");
                            sprintf(m1[5],"LambdapKp_Lambdakx");
                            sprintf(m1[6],"LambdaaKmMix_Lambdakx");
                            sprintf(m1[7],"LambdapKpMix_Lambdakx");
                            sprintf(m1[8],"Xi1820_p_gen_Lambdakx");
                            sprintf(m1[9],"Xi1820_m_gen_Lambdakx");
                            sprintf(m1[10],"Xi1820_p_true_Lambdakx");
                            sprintf(m1[11],"Xi1820_m_true_Lambdakx");
                            sprintf(m1[12],"Xi1820_p_trueMM_Lambdakx");
                            sprintf(m1[13],"Xi1820_m_trueMM_Lambdakx");
                            sprintf(m1[14],"Xi1820_p_res_Lambdakx");
                            sprintf(m1[15],"Xi1820_m_res_Lambdakx");
                            nm2=8;
                            if (cutloop >= 1) {
                                sprintf(DC2,"Lambdakx_cut%i_Lambdakx",cutloop);
                                sprintf(m1[0],"LambdaaKp_Lambdakx_cut%i",cutloop);
                                sprintf(m1[1],"LambdapKm_Lambdakx_cut%i",cutloop);
                                sprintf(m1[2],"LambdaaKpMix_Lambdakx_cut%i",cutloop);
                                sprintf(m1[3],"LambdapKmMix_Lambdakx_cut%i",cutloop);
                                sprintf(m1[4],"LambdaaKm_Lambdakx_cut%i",cutloop);
                                sprintf(m1[5],"LambdapKp_Lambdakx_cut%i",cutloop);
                                sprintf(m1[6],"LambdaaKmMix_Lambdakx_cut%i",cutloop);
                                sprintf(m1[7],"LambdapKpMix_Lambdakx_cut%i",cutloop);
                                sprintf(m1[8],"Xi1820_p_gen_Lambdakx_cut%i",cutloop);
                                sprintf(m1[9],"Xi1820_m_gen_Lambdakx_cut%i",cutloop);
                                sprintf(m1[10],"Xi1820_p_true_Lambdakx_cut%i",cutloop);
                                sprintf(m1[11],"Xi1820_m_true_Lambdakx_cut%i",cutloop);
                                sprintf(m1[12],"Xi1820_p_trueMM_Lambdakx_cut%i",cutloop);
                                sprintf(m1[13],"Xi1820_m_trueMM_Lambdakx_cut%i",cutloop);
                                sprintf(m1[14],"Xi1820_p_res_Lambdakx_cut%i",cutloop);
                                sprintf(m1[15],"Xi1820_m_res_Lambdakx_cut%i",cutloop);
                            }
                        }
                        else{
                            if (hmtloop == 0) {//992
                                sprintf(DC2,"Lambdakx_Lambdakx");
                                sprintf(m1[0],"LambdaaKp_Lambdakx");
                                sprintf(m1[1],"LambdapKm_Lambdakx");
                                sprintf(m1[2],"LambdaaKpMix_Lambdakx");
                                sprintf(m1[3],"LambdapKmMix_Lambdakx");
                                sprintf(m1[4],"LambdaaKm_Lambdakx");
                                sprintf(m1[5],"LambdapKp_Lambdakx");
                                sprintf(m1[6],"LambdaaKmMix_Lambdakx");
                                sprintf(m1[7],"LambdapKpMix_Lambdakx");
                                if (cutloop >= 1) {
                                    sprintf(DC2,"Lambdakx_cut%i_Lambdakx",cutloop);
                                    sprintf(m1[0],"LambdaaKp_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[1],"LambdapKm_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[2],"LambdaaKpMix_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[3],"LambdapKmMix_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[4],"LambdaaKm_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[5],"LambdapKp_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[6],"LambdaaKmMix_Lambdakx_cut%i",cutloop);
                                    sprintf(m1[7],"LambdapKpMix_Lambdakx_cut%i",cutloop);
                                }
                            }
                            else{//993
                                sprintf(DC2,"Lambdakx_hmt_Lambdakx");
                                sprintf(m1[0],"LambdaaKp_Lambdakx_hmt");
                                sprintf(m1[1],"LambdapKm_Lambdakx_hmt");
                                sprintf(m1[2],"LambdaaKpMix_Lambdakx_hmt");
                                sprintf(m1[3],"LambdapKmMix_Lambdakx_hmt");
                                sprintf(m1[4],"LambdaaKm_Lambdakx_hmt");
                                sprintf(m1[5],"LambdapKp_Lambdakx_hmt");
                                sprintf(m1[6],"LambdaaKmMix_Lambdakx_hmt");
                                sprintf(m1[7],"LambdapKpMix_Lambdakx_hmt");
                                if (cutloop >= 1) {
                                    sprintf(DC2,"Lambdakx_hmt_cut%i_Lambdakx",cutloop);
                                    sprintf(m1[0],"LambdaaKp_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[1],"LambdapKm_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[2],"LambdaaKpMix_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[3],"LambdapKmMix_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[4],"LambdaaKm_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[5],"LambdapKp_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[6],"LambdaaKmMix_Lambdakx_hmt_cut%i",cutloop);
                                    sprintf(m1[7],"LambdapKpMix_Lambdakx_hmt_cut%i",cutloop);
                                }
                            }
                        }
                        ///sprintf(DC2,"Lambdakx_ppData_Lambdakx");///old format
                    }else if(jDC==Klk0){
                        
                        if (MCflag == 1) {//994
                            sprintf(DC2,"Lambdak0_Lambdak0");
                            sprintf(m1[0],"LambdaaK0_Lambdak0");
                            sprintf(m1[1],"LambdapK0_Lambdak0");
                            sprintf(m1[2],"LambdaaK0Mix_Lambdak0");
                            sprintf(m1[3],"LambdapK0Mix_Lambdak0");
                            sprintf(m1[4],"Xi1820_0p_gen_Lambdak0");
                            sprintf(m1[5],"Xi1820_0a_gen_Lambdak0");
                            sprintf(m1[6],"Xi1820_0p_true_Lambdak0");
                            sprintf(m1[7],"Xi1820_0a_true_Lambdak0");
                            sprintf(m1[8],"Xi1820_0p_trueMM_Lambdak0");
                            sprintf(m1[9],"Xi1820_0a_trueMM_Lambdak0");
                            sprintf(m1[10],"Xi1820_0p_res_Lambdak0");
                            sprintf(m1[11],"Xi1820_0a_res_Lambdak0");
                            nm2=8;
                            if (cutloop >= 1) {
                                sprintf(DC2,"Lambdak0_cut%i_Lambdak0",cutloop);
                                sprintf(m1[0],"LambdaaK0_Lambdak0_cut%i",cutloop);
                                sprintf(m1[1],"LambdapK0_Lambdak0_cut%i",cutloop);
                                sprintf(m1[2],"LambdaaK0Mix_Lambdak0_cut%i",cutloop);
                                sprintf(m1[3],"LambdapK0Mix_Lambdak0_cut%i",cutloop);
                                sprintf(m1[4],"Xi1820_0p_gen_Lambdak0_cut%i",cutloop);
                                sprintf(m1[5],"Xi1820_0a_gen_Lambdak0_cut%i",cutloop);
                                sprintf(m1[6],"Xi1820_0p_true_Lambdak0_cut%i",cutloop);
                                sprintf(m1[7],"Xi1820_0a_true_Lambdak0_cut%i",cutloop);
                                sprintf(m1[8],"Xi1820_0p_trueMM_Lambdak0_cut%i",cutloop);
                                sprintf(m1[9],"Xi1820_0a_trueMM_Lambdak0_cut%i",cutloop);
                                sprintf(m1[10],"Xi1820_0p_res_Lambdak0_cut%i",cutloop);
                                sprintf(m1[11],"Xi1820_0a_res_Lambdak0_cut%i",cutloop);
                            }
                        }
                        else{
                            if (hmtloop == 0) {//992
                                sprintf(DC2,"Lambdak0_Lambdak0");
                                sprintf(m1[0],"LambdaaK0_Lambdak0");
                                sprintf(m1[1],"LambdapK0_Lambdak0");
                                sprintf(m1[2],"LambdaaK0Mix_Lambdak0");
                                sprintf(m1[3],"LambdapK0Mix_Lambdak0");
                                if (cutloop >= 1) {
                                    sprintf(DC2,"Lambdak0_cut%i_Lambdak0",cutloop);
                                    sprintf(m1[0],"LambdaaK0_Lambdak0_cut%i",cutloop);
                                    sprintf(m1[1],"LambdapK0_Lambdak0_cut%i",cutloop);
                                    sprintf(m1[2],"LambdaaK0Mix_Lambdak0_cut%i",cutloop);
                                    sprintf(m1[3],"LambdapK0Mix_Lambdak0_cut%i",cutloop);
                                }
                            }
                            else{//993
                                sprintf(DC2,"Lambdak0_hmt_Lambdak0");
                                sprintf(m1[0],"LambdaaK0_Lambdak0_hmt");
                                sprintf(m1[1],"LambdapK0_Lambdak0_hmt");
                                sprintf(m1[2],"LambdaaK0Mix_Lambdak0_hmt");
                                sprintf(m1[3],"LambdapK0Mix_Lambdak0_hmt");
                                if (cutloop >= 1) {
                                    sprintf(DC2,"Lambdak0_hmt_cut%i_Lambdak0",cutloop);
                                    sprintf(m1[0],"LambdaaK0_Lambdak0_hmt_cut%i",cutloop);
                                    sprintf(m1[1],"LambdapK0_Lambdak0_hmt_cut%i",cutloop);
                                    sprintf(m1[2],"LambdaaK0Mix_Lambdak0_hmt_cut%i",cutloop);
                                    sprintf(m1[3],"LambdapK0Mix_Lambdak0_hmt_cut%i",cutloop);
                                }
                            }
                        }
                        nm=4;
                    }
                    
                    /// THnSparse* s[8];///old format
                    //s[16];//TH3F*
                    
                    //H[8];//TH2D*
                    
                    for(j=0;j<(nm+nm2);j++){
                        sprintf(s0,"%s_%s",DC2,m1[j]);
                        //s[j]=(THnSparse*) l->FindObject(s0);//get the multi-dimensional histograms
                        s[j]=(TH3F*) l->FindObject(s0);//get the multi-dimensional histograms
                        if(!s[j]){cerr<<"missing histogram "<<s0<<endl; return;}
                    }
                    
                    //multiplicity centrality incomplete
                    /*
                    sprintf(s0,"MultiVsCent");
                    M[0]=(TH2F*) l->FindObject(s0);//get the multi-dimensional histograms
                    if(!M[0]){cerr<<"missing histogram "<<s0<<endl; return;}
                    M[0]->GetXaxis()->SetRange(M[0]->GetXaxis()->FindBin(cbmin[cb]+1.e-7),M[0]->GetXaxis()->FindBin(cbmax[cb]-1.e-7));//set axis ranges to project in multiplicity slices
                    M2=(TH1F*) M[0]->ProjectionX("_px",s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                     */
                    //multiplicity centrality incomplete
                    
                    //char s2[100];//character array for Lambda
                    //char s3[100];//character array for hmt
                    if (jDC == Klkx) {
                        sprintf(s2,"LambdaKX");
                    }
                    else if(jDC == Klk0){
                        sprintf(s2,"LambdaK0");
                    }
                    else if(jDC == Klpi){
                        sprintf(s2,"LambdaPi");
                    }
                    
                    if (MCflag == 1) {
                        sprintf(s3,"MC");
                    }
                    else{
                        if (hmtloop == 0) {
                            sprintf(s3,"Normal");
                        }
                        else if(hmtloop == 1){
                            sprintf(s3,"hmt");
                        }
                    }
                    
                    if (cutloop >= 1) {
                        fout=new TFile(Form("~/Desktop/xitest/plotsAllV9/%s/%s/cut%i/histogramsAll.root",s3,s2,cutloop),"RECREATE","HistoFile");//create the output file
                    }
                    else{
                        fout=new TFile(Form("~/Desktop/xitest/plotsAllV9/%s/%s/histogramsAll.root",s3,s2),"RECREATE","HistoFile");//create the output file
                    }
                    
                    //TH1D *h[30],*b,*d;
                    for(cb=0;cb<ncb && ptbins();cb++) for(pb=0;pb<npb;pb++){//loop over multiplicity and pT bins
                        if (hmtloop == 1) {
                            ncb=1;
                        }
                        
                        for(j=0;j<30;j++) h[j]=0;
                        
                        for(j=0;j<nm;j++){
                            s[j]->GetYaxis()->SetRange(s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4));//set axis ranges to project in pT slices
                            s[j]->GetZaxis()->SetRange(s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7));//set axis ranges to project in multiplicity slices
                            ///h[j]=(TH1D*) s[j]->Projection(0,"e");//project histograms
                            h[j]=(TH1D*) s[j]->ProjectionX("_px",s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                            ////h[j]=(TH1D*) s[j]->ProjectionY("_py",s[j]->GetXaxis()->FindBin(1.60),s[j]->GetXaxis()->FindBin(2.40),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                            
                            h[j]->SetName(Form("mass_%s_c%i_p%i",m0[j],cb,pb));
                            //cerr << Form("mass_%s_c%i_p%i",m0[j],cb,pb) << endl;
                        }
                        if (MCflag == 1) {//return here
                            for (j=nm; j<(nm+nm2); j++) {
                                //if using 1-d -> pt
                                if ((j-nm) < 6) {
                                    s[j]->GetXaxis()->SetRange(s[j]->GetXaxis()->FindBin(massmin[0]+1.e-4),s[j]->GetXaxis()->FindBin(massmax[0]-1.e-4));//set axis ranges to project in mass slices
                                    s[j]->GetZaxis()->SetRange(s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7));//set axis ranges to project in multiplicity slices
                                    /*h[j]=(TH1D*) s[j]->ProjectionY("_py",s[j]->GetXaxis()->FindBin(massmin[0]+1.e-4),s[j]->GetXaxis()->FindBin(massmax[0]-1.e-4),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                                    
                                    h[j]->SetName(Form("pt_%s_c%i_m%i",m0[14+(j-nm)],cb,pb));//m%i, pb is only kept to keep loop simple*/
                                    h[j]=(TH1D*) s[j]->ProjectionX("_px",s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                                    
                                    h[j]->SetName(Form("mass_%s_c%i_p%i",m0[14+(j-nm)],cb,pb));//m%i, pb is only kept to keep loop simple
                                    
                                    //h[j]->SetName(Form("mass_%s_c%i_p%i",m0[14+(j-nm)],cb,pb));
                                    //cerr << Form("mass_%s_c%i_p%i",m0[14+(j-nm)],cb,pb) << endl;
                                }
                                else{
                                    s[j]->GetYaxis()->SetRange(s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4));//set axis ranges to project in pT slices
                                    s[j]->GetZaxis()->SetRange(s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7));//set axis ranges to project in multiplicity slices
                                    h[j]=(TH1D*) s[j]->ProjectionX("_px",s[j]->GetYaxis()->FindBin(ptmin[pb]+1.e-4),s[j]->GetYaxis()->FindBin(ptmax[pb]-1.e-4),s[j]->GetZaxis()->FindBin(cbmin[cb]+1.e-7),s[j]->GetZaxis()->FindBin(cbmax[cb]-1.e-7),"e");//project histograms
                                    
                                    h[j]->SetName(Form("mass_%s_c%i_p%i",m0[14+(j-nm)],cb,pb));
                                    //cerr << Form("mass_%s_c%i_p%i",m0[14+(j-nm)],cb,pb) << endl;
                                }
                                
                            }
                        }
                        
                        h[8+nm2]=(TH1D*) h[0]->Clone(Form("mass_%s_c%i_p%i",m0[8],cb,pb));//unlike
                        h[8+nm2]->Add(h[1],1.);
                        
                        ////cout << "mean (pT) calculation for Cen " << cbmin[cb] << "-" << cbmax[cb] << " and pt " << ptmin[pb] << "-" << ptmax[pb] << " is " << h[8+nm2]->GetMean() << endl;
                        
                        h[9+nm2]=(TH1D*) h[2]->Clone(Form("mass_%s_c%i_p%i",m0[9],cb,pb));//mixunlike
                        h[9+nm2]->Add(h[3],1.);
                        
                        if(!scheme || scheme==2){
                            h[10+nm2]=(TH1D*) h[4]->Clone(Form("mass_%s_c%i_p%i",m0[10],cb,pb));//sum
                            h[10+nm2]->Add(h[5],1.);
                            
                            h[11+nm2]=(TH1D*) h[4]->Clone(Form("mass_%s_c%i_p%i",m0[11],cb,pb));//gm
                            for(j=0;j<=h[11+nm2]->GetNbinsX()+1;j++){
                                h[11+nm2]->SetBinContent(j,2.*sqrt(h[4]->GetBinContent(j)*h[5]->GetBinContent(j)));
                                h[11+nm2]->SetBinError(j,h[10+nm2]->GetBinError(j));
                            }
                        }
                        
                        if(scheme==2){
                            h[12+nm2]=(TH1D*) h[6]->Clone(Form("mass_%s_c%i_p%i",m0[12],cb,pb));//mixlike
                            h[12+nm2]->Add(h[7],1.);
                            
                            h[13+nm2]=(TH1D*) h[9]->Clone(Form("mass_%s_c%i_p%i",m0[13],cb,pb));//mixall
                            h[13+nm2]->Add(h[12],1.);
                        }
                        
                        if (MCflag == 1) {
                            //to have 1-d -> pt plots
                            h[22]=(TH1D*) h[nm]->Clone(Form("pt_%s_c%i_m%i",m0[22],cb,pb));//generated
                            h[22]->Add(h[nm+1],1.);
                            
                            h[23]=(TH1D*) h[nm+2]->Clone(Form("pt_%s_c%i_m%i",m0[23],cb,pb));//true
                            h[23]->Add(h[nm+3],1.);
                            
                            h[24]=(TH1D*) h[nm+4]->Clone(Form("pt_%s_c%i_m%i",m0[24],cb,pb));//trueMM
                            h[24]->Add(h[nm+5],1.);
                            
                            //h[25]=(TH1D*) h[nm+6]->Clone(Form("mass_%s_c%i_p%i",m0[25],cb,pb));
                            h[25]=(TH1D*) h[nm+6]->Clone(Form("pt_%s_c%i_p%i",m0[25],cb,pb));//resolution
                            h[25]->Add(h[nm+7],1.);
                            
                            //rebin according to Anders
                            h[24]->Rebin(20);
                            h[22]->Rebin(20);
                            h[26]=(TH1D*) h[24]->Clone(Form("pt_%s_c%i_m%i",m0[26],cb,pb));//acceptance*efficiency
                            h[26]->Divide(h[22]);
                            
                        }
                        
                        fout->cd();
                        for(j=0;j<30;j++) if(h[j]) h[j]->Write();//save histograms
                        
                        //-----
                        
                        if(!scheme || scheme==2){
                            for(j=0;j<7;j++){
                                nmix=1.;
                                //calculate the difference between h[i1] and h[i2], like-charge backgrounds
                                if(!j){i1=0; i2=4;}//np-nn
                                else if(j==1){i1=0; i2=5;}//np-pp
                                else if(j==2){i1=0; i2=11; nmix=0.5;}//np-gm/2
                                else if(j==3){i1=1; i2=4;}//pn-nn
                                else if(j==4){i1=1; i2=5;}//pn-pp
                                else if(j==5){i1=1; i2=11; nmix=0.5;}//pn-gm/2
                                else if(j==6){i1=8; i2=11;}//unlike-gm
                                
                                d=(TH1D*) h[i1]->Clone(Form("mass_%sQ%s_c%i_p%i",m0[i1],m0[i2],cb,pb));//"Q" identifies the histogram after background subtraction
                                d->Add(h[i2],-nmix);//error if i2 histogram has higher values then i1 histogram
                                /*
                                 if (fitflag == 1) {
                                 double par[5]=0;
                                 for (int z = 0; z < 5; z++) {
                                 par[z]=0;
                                 }
                                 if (particleloop <= 1) {
                                 TF1 *fitfunction = new TF1("fitfunction",function,1.7,2.0,5);
                                 //fitfunction->SetParLimits(0, 0.0, 4000.0);
                                 fitfunction->SetParLimits(1,0.01,0.10);
                                 fitfunction->SetParLimits(2,1.8,1.85);
                                 //fitfunction->SetParLimits(3,-2000.,8000.);
                                 //fitfunction->SetParLimits(4,-5000.,5000.);
                                 d->Fit("fitfunction","Q","SAME",1.7,2.0);
                                 for (z=0; z < 5; z++) {
                                 cout << "parmeter " << z << " " << fitfunction->GetParameter(z) << endl;
                                 }
                                 //cout << "chi2 " << fitfunction->GetChisquare() << endl;
                                 }
                                 }*/
                                fout->cd();
                                d->Write();
                            }
                        }
                        
                        //-----
                        
                        for(j=0;j<12;j++){
                            //calculate the difference between h[i1] and h[i2], mixed-event background
                            if(!j){i1=0; i2=2;}//np-mixnp
                            else if(j==1){i1=0; i2=9;}//np-mixunlike
                            else if(j==2){i1=1; i2=3;}//pn-mixpn
                            else if(j==3){i1=1; i2=9;}//pn-mixunlike
                            else if(j==4){i1=8; i2=9;}//unlike-mixunlike
                            else if(j==5){i1=4; i2=6;}//nn-mixnn
                            else if(j==6){i1=4; i2=12;}//nn-mixlike
                            else if(j==7){i1=5; i2=7;}//pp-mixpp
                            else if(j==8){i1=5; i2=12;}//pp-mixlike
                            else if(j==9){i1=10; i2=12;}//sum-mixlike
                            else if(j==10){i1=9; i2=13;}//unlike-mixall
                            else if(j==11){i1=8; i2=11;}//unlike-gm
                            if(scheme!=2 && j>=5) continue;
                            
                            for(jNR=0;jNR<nNR && NRbins();jNR++){
                                b=(TH1D*) h[i2]->Clone(Form("mass_%sN%s%i_c%i_p%i",m0[i1],m0[i2],jNR,cb,pb));//"N" identifies histogram h[i2] normalized to histogram h[i1] in the region jNR
                                d=(TH1D*) h[i1]->Clone(Form("mass_%sQ%s%i_c%i_p%i",m0[i1],m0[i2],jNR,cb,pb));//"Q" identifies the histogram after background subtraction
                                normalize_mix(h[i1],h[i2],b,d);//fill b and d
                                fout->cd();
                                /*if (fitflag == 1) {
                                 double par[5]=0;
                                 for (int z = 0; z < 5; z++) {
                                 par[z]=0;
                                 }
                                 if (particleloop <= 1) {
                                 TF1 *fitfunction = new TF1("fitfunction",function,1.7,2.0,5);
                                 TF1 *fitline = new TF1("fitline",fline,1.7,2.0,2);//new
                                 //fitfunction->SetParLimits(2, 0.0, 4000.0);
                                 fitfunction->SetParLimits(3,0.01,0.10);
                                 fitfunction->SetParLimits(4,1.8,1.85);
                                 //fitfunction->SetParLimits(3,-2000.,8000.);
                                 //fitfunction->SetParLimits(4,-5000.,5000.);
                                 //bool reject;
                                 //reject = kTRUE;
                                 d->Fit("fitline","1","SAME",1.7,2.0);
                                 //reject = kFALSE;
                                 
                                 d->Fit("fitfunction","Q","SAME",1.7,2.0);
                                 for (z=0; z < 5; z++) {
                                 cout << "parmeter " << z << " " << fitfunction->GetParameter(z) << endl;
                                 }
                                 //cout << "chi2 " << fitfunction->GetChisquare() << endl;
                                 }
                                 
                                 }*/
                                
                                b->Write();
                                d->Write();
                            }
                        }
                    }
                }//end cutloop
                
                fin->Close();
                fout->Close();
            }//end childloop
            
        }//end particle loop
        
    }//end hmtloop
    return;
}
