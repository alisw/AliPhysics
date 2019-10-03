/* 14 June 2015
Introduced the scaling for width f the ZDC distribution also. 
Earlier the means were scaled byut there was a difference of about 4% in the width.
This can cause difference in the value of asymmetry parameters + and -ve values. 
The sacling is notw performed for sigma after scaling with the mean. 
So that the two scaled distributions should be same. The motivation is same as that 
for scaling with the mean in the first place.
For this the QA plots for each run are read from the root file, the distribution is fit and the mean and sigma is obtained. Then the new ZDC scaled values are calculated and alpha ZDC per event is calculated filled in a 
Histo. The histo is also saved in the same root file.  

This files makes all the QA plots for each run. 
*/

TFile * f,*f2;
TTree * tree, *asymTree, *FBTree;
Char_t fname_dat[30];

const Int_t NMaxTrack = 8000;
const Int_t NMaxTracklet = 8000;
Int_t fRunNumber, fNumberOfTracklets,fNumberOfTracks;
//Int_t fProjSpectP,fProjSpectN,fTargSpectP,fTargSpectN;
Float_t fCentPercentile[4];
Float_t fVertexX,fVertexY,fVertexZ;
Float_t fzdcn1en,fzdcn2en,fzdcp1en,fzdcp2en,fzdcEM1en,fzdcEM2en;
Float_t fTrackletEta[NMaxTracklet],fTrackletPhi[NMaxTracklet];
Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack];
Int_t fTrackCharge[NMaxTrack];
Int_t fTrackFilterMap[NMaxTrack];
Float_t fV0Mult[64];

Int_t nentries;

TH2F * hVertexXY, *hzdcn1vsn2;
TH1F *hVertexZ, *hNTracklets, *hNTRacks, *hCent[4];
TH1F * hNTracks, *hNTrackslets;
TH1D * hTrackEta, *hTrackPhi, *hTrackletEta, *hTrackletPhi;
TH2F * hNTracks_Cent, *hNTracklets_Cent, *hZDCn_Cent;

// This holds unscaled mean and sigma
Float_t fzdcn1Mean, fzdcn1Sigma, fzdcn2Mean, fzdcn2Sigma;

// This holds scaled mean and sigma
Float_t zdcn1Mean, zdcn1Sigma, zdcn2Mean, zdcn2Sigma;

Float_t avzdcn1en=1.0,avzdcn2en=1.0, avzdcEM1en=1.0,avzdcEM2en=1.0;
Float_t avV0AMult=1.0, avV0CMult=1.0;
Float_t zdcasym  =1.0, zdcemasym = 1.0, trkasym=1.0, trkletasym = 1.0;
Float_t PartAsym = 1.0, SpectAsym = 1.0, NPartAsym = 1.0, NSpectAsym = 1.0;
Float_t v0asym = 1.0, pTasym = 1.0;

// If you change this.. remember to change the size of array inthe FBTree also 
const Int_t netagap = 10;
const Float_t etabinsz = 0.4;
Int_t nposTracks, nnegTracks, nposTracklets, nnegTracklets;
Int_t nposTrksInDEta[netagap], nnegTrksInDEta[netagap];
TProfile *hPosTrksInDEta[netagap], *hNegTrksInDEta[netagap];
Float_t  avPosNTracks=1.0, avNegNTracks=1.0, avPosNTracklets=1.0, avNegNTracklets=1.0;

TH1F *hZDCN1, *hZDCN2, *hZDCEM1, *hZDCEM2, *hPosNTracks, * hNegNTracks, * hPosNTracklets, * hNegNTracklets;
TH1F * hPospT, *hNegpT;
TH1F * hV0A, *hV0C, *hV0AScaled ,*hV0CScaled;
TH1F * hZDCN1Scaled, *hZDCN2Scaled, *hZDCEM1Scaled, *hZDCEM2Scaled;
TH1F * hZDCN1SigmaScaled, *hZDCN2SigmaScaled;
TH1F * hNTracks1Scaled, *hNTracks2Scaled;
TH1F * hNTracklets1Scaled, *hNTracklets2Scaled;
TH2F * hNTracksvsCent, *hNTrackletsvsCent, *hZDCn_Cent;
//TH2F * hNPart_ZDCN, *hNnSpect_ZDCN;

TH1F * hVertexZ4ZDCAsym[20];
TH1F * hNCh[20];
TH1F * hTrackEta4ZDCAsym[20];
TProfile * hpZDCAsym_pT;
TProfile * hpZDCAsym_NCh;
TProfile * hpZDCAsym_Vz;
TProfile * hpZDCAsym_Cent;
Int_t NEvCut_Cent=0, NEvCut_VertexZ=0, NEvCut_VertexXY=0, NEvCut_NTrack=0, NEvCut_NTracklet=0, NEvTotal=0, NEvAccepted=0
  , NEvCut_ZDC=0;

Float_t centmin=0.0,centmax=10.0;
Float_t vzmin = -5.0, vzmax = 5.0;
Float_t vxymin = -0.3, vxymax = 0.3;
Float_t trkletetamin = -1.0, trkletetamax = 1.0;
Float_t trketamin = -0.9, trketamax = 0.9;
Float_t DEta1min[netagap],DEta1max[netagap],DEta2min[netagap],DEta2max[netagap];
/*
Float_t DEta1min[netagap] ={-0.5,-0.6,-0.7,-0.8,-0.9};
Float_t DEta1max[netagap] ={-0.0,-0.1,-0.2,-0.3,-0.4}; 
Float_t DEta2max[netagap] ={0.5,0.6,0.7,0.8,0.9};
Float_t DEta2min[netagap] ={0.0,0.1,0.2,0.3,0.4}; 
*/


Int_t nfile;
//Int_t nrun = 2;
Int_t nrun = 60;
Int_t RunNumber[]={
  
  137231, 137232, 137236, 137243, 137432,
  137434, 137440, 137441, 137443, 137541,
  137544, 137686, 137691, 137692, 137704,
  137718, 137722, 137724, 137751, 137752,
  137844, 137848,
  
  138190, 138192, 138197, 138201, 138225,
  138275, 138364, 138396, 138438, 138439,
  138442, 138469, 138534, 138578, 138579,
  138582, 138583,
  
  139028,  139029,139036, 139037, 139038,
  139105, 139107, 139173, 139309, 139310,
  139314, 139328, 139329, 139360, 139437,
  139438, 139465, 139503, 139505, 139507,
  139510
};

Int_t RunMin, RunMax; 
Int_t Debug=0;

void readAOD4QA(Float_t incentmin = 10.0, Float_t incentmax = 20.0, Int_t inRunMin = 0, Int_t inRunMax = 65){
  
  centmin = incentmin;
  centmax = incentmax;
 
  RunMin = inRunMin;
  RunMax = inRunMax;
  
  //  TFile  *f = new TFile("data/AOD.137161.root");
  BookHisto();

  
  cout<<" Studying events in centrality bin "<<centmin<<" to "<<centmax<<endl;    

  for(Int_t irun = RunMin;irun<=RunMax;irun++){
    InitRun();
    BookAsymTree(RunNumber[irun]);
    nfile = 1;
    if(RunNumber[irun]==138275) nfile = 5;
    if(RunNumber[irun]==138364) nfile = 3;
    if(RunNumber[irun]==138396) nfile = 2;
    if(RunNumber[irun]==138442) nfile = 2;
    if(RunNumber[irun]==138534) nfile = 4;
    if(RunNumber[irun]==138578) nfile = 2;
    if(RunNumber[irun]==138583) nfile = 2;
    
    FillQAPlots(RunNumber[irun],nfile);
    cout<<" Back to main"<<endl;
    PrintCuts(1);
    GetAverages(RunNumber[irun]);
    
    readAverages(RunNumber[irun]);
    FillScaledQAPlots(RunNumber[irun],nfile);
    GetAllAsymmetry(RunNumber[irun],nfile);
  }

}

void readAverages(Int_t runNumber){
  
  Char_t fname_dat[120];
  sprintf(fname_dat,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/AsymRun%dCent%3.3dto%3.3d.dat",centmin,centmax,runNumber,centmin,centmax);
  cout<<" Reading data file "<<fname_dat<<" for averages "<<endl;
  ifstream fin(fname_dat);
  Int_t runnumber, nevents;
  fin>>runnumber>>nevents;
  if(runnumber!=runNumber) {
    cout<<" Probably not the correct averages. Please check!"<<endl;
    cout<<" RunNumber being processed="<<runNumber<<endl;
  }
  fin>>avzdcn1en>>avzdcn2en;
  fin>>avPosNTracks>>avNegNTracks;
  fin>>avPosNTracklets>>avNegNTracklets;
  fin>>avV0AMult>>avV0CMult;
  
  cout<<"<ZDCn1en> = "<<avzdcn1en<<" <ZDCn2en> = "<<avzdcn2en<<endl; 
  cout<<" <NposTracks>="<<avPosNTracks<<" <NNegTracks>="<<avNegNTracks<<endl; 
  cout<<" <NposTracklets>="<<avPosNTracklets<<" <NNegTracklets>="<<avNegNTracklets<<endl; 
  cout<<" avV0AMult ="<<avV0AMult<<" avV0CMult = "<<avV0CMult<<endl;
}

void PrintCuts(Int_t ntimes = 1){
  cout<<" Total number of events read = "<<NEvTotal/ntimes<<endl;
  if(ntimes==1)
  cout<<" Total number of events accepted ="<<hVertexZ->GetEntries()<<endl;
  if(ntimes==2)
      cout<<" Total number of events accepted ="<<hNTracks1Scaled->GetEntries()<<endl;
  cout<<" Number of Events Cut due to Centrality ("<<centmin<<" / "<<centmax<<")="<<NEvCut_Cent/ntimes<<endl;
  cout<<" Number of Events Cut due to VertexZ ("<<vzmin<<" / "<<vzmax<<")="<<NEvCut_VertexZ/ntimes<<endl;
  cout<<" Number of Events Cut due to VertexXY ("<<vxymin<<" / "<<vxymax<<")="<<NEvCut_VertexXY/ntimes<<endl;
  cout<<" Number of Events Cut due to NTracks or NTracklets being ("<<0<<")="<<NEvCut_NTrack/ntimes<<endl;
  cout<<" Number of Events Cut due to ZDC (Neutron / em) being ("<<0<<")="<<NEvCut_ZDC/ntimes<<endl;
}


Float_t GetAsymmetry( Double_t mult1, Double_t meanmult1, Double_t mult2, Double_t meanmult2)
{
  Double_t asym = 1.0;
  //  if(((mult1/meanmult1)+(mult2/meanmult2))>0)
  asym = ((mult1/meanmult1)-(mult2/meanmult2))/((mult1/meanmult1)+(mult2/meanmult2));
  //  cout<<" input = "<<mult1<<" "<<meanmult1<<" "<<mult2<<" "<<meanmult2<<" "<<asym<<endl;
  return asym;
} 

void BookTree(TTree* tree){

  if(!tree) return;
  
  tree->SetBranchAddress("fRunNumber",&fRunNumber);
  //  tree->SetBranchAddress("fProjSpectP",&fProjSpectP);
  //  tree->SetBranchAddress("fProjSpectN",&fProjSpectN);
  //  tree->SetBranchAddress("fTargSpectP",&fTargSpectP);
  //  tree->SetBranchAddress("fTargSpectN",&fTargSpectN);
  tree->SetBranchAddress("fNumberOfTracks",&fNumberOfTracks);
  tree->SetBranchAddress("fNumberOfTracklets",&fNumberOfTracklets);
  tree->SetBranchAddress("fVertexX",&fVertexX);
  tree->SetBranchAddress("fVertexY",&fVertexY);
  tree->SetBranchAddress("fVertexZ",&fVertexZ);
  tree->SetBranchAddress("fCentPercentile",fCentPercentile);
  tree->SetBranchAddress("fzdcn1en",&fzdcn1en);
  tree->SetBranchAddress("fzdcn2en",&fzdcn2en);
  tree->SetBranchAddress("fzdcp1en",&fzdcp1en);
  tree->SetBranchAddress("fzdcp2en",&fzdcp2en);
  tree->SetBranchAddress("fzdcEM1en",&fzdcEM1en);
  tree->SetBranchAddress("fzdcEM2en",&fzdcEM2en);
  tree->SetBranchAddress("fTrackletEta",fTrackletEta);
  tree->SetBranchAddress("fTrackletPhi",fTrackletPhi);
  tree->SetBranchAddress("fTrackEta",fTrackEta);
  tree->SetBranchAddress("fTrackPhi",fTrackPhi);
  tree->SetBranchAddress("fTrackPt",fTrackPt);
  tree->SetBranchAddress("fTrackCharge",fTrackCharge);
  tree->SetBranchAddress("fTrackFilterMap",fTrackFilterMap);
  tree->SetBranchAddress("fV0Mult",fV0Mult);
}

void BookHisto(){
  
/*
  hNPart_ZDCN = new TH2F("hNPart_ZDCN","hNPart_ZDCN",100,0.0,80000,100,0,425);
  hNPart_ZDCN->SetXTitle(" Neutron ZDC Energy A+C");
  hNPart_ZDCN->SetYTitle("Number of Participants");

  hNnSpect_ZDCN = new TH2F("hNnSpect_ZDCN","hNnSpect_ZDCN",100,0.0,80000,100,0,275);
  hNnSpect_ZDCN->SetXTitle(" Neutron ZDC Energy A+C");
  hNnSpect_ZDCN->SetYTitle("Number of Neutron Spectators");
*/
  hNTracks_Cent = new TH2F("hNTracks_Cent","hNTracks_Cent",10,centmin,centmax,140,0.0,7000);
  hNTracklets_Cent  = new TH2F("hNTracklets_Cent","hNTracklets_Cent",10,centmin,centmax,100,0.0,5000);
  hZDCn_Cent = new TH2F("hZDCn_Cent","hZDCn_Cent",10,centmin,centmax,50,0.0,100000);
  
  hVertexXY= new TH2F("hVertexXY","hVertexXY",60,-0.3,0.3,60,-0.3,0.3);
  hVertexXY->SetXTitle("Vx");
  hVertexXY->SetYTitle("Vy");
  
  hVertexZ= new TH1F("hVertexZ","hVertexZ",40,-10.0,10.0);
  hVertexZ->SetXTitle("Vz");
  hVertexZ->SetYTitle("Number of Events");
  
  hzdcn1vsn2= new TH2F("hzdcn1vsn2","hzdcn1vsn2",100,0.0,100000,100,0.0,100000);
  hzdcn1vsn2->SetXTitle("hzdcn1");
  hzdcn1vsn2->SetYTitle("hzdcn2");
  
  /*
  hpzdcn1vsn2= new TProfile("hpzdcn1vsn2","hpzdcn1vsn2",100,0.0,100000);
  hpzdcn1vsn2->SetXTitle("hzdcn1");
  hpzdcn1vsn2->SetYTitle("hzdcn2");
  */

  
  hNTracks1Scaled = new TH1F("hNTracks1Scaled","hNTracks1Scaled",400,0.0,8.0);
  hNTracks1Scaled->SetXTitle("Number of +ve #eta Tracks Scaled");
  hNTracks1Scaled->SetYTitle("Number of Events");
  
  hNTracklets1Scaled = new TH1F("hNTracklets1Scaled","hNTracklets1Scaled",400,0.0,8.0);
  hNTracklets1Scaled->SetXTitle("Number of +ve #eta Tracklets Scaled");
  hNTracklets1Scaled->SetYTitle("Number of Events");
  
  hNTracks2Scaled = new TH1F("hNTracks2Scaled","hNTracks2Scaled",400,0.0,8.0);
  hNTracks2Scaled->SetXTitle("Number of -ve #eta Tracks Scaled");
  hNTracks2Scaled->SetYTitle("Number of Events");
  
  hNTracklets2Scaled = new TH1F("hNTracklets2Scaled","hNTracklets2Scaled",400,0.0,8.0);
  hNTracklets2Scaled->SetXTitle("Number of -ve #eta Tracklets Scaled");
  hNTracklets2Scaled->SetYTitle("Number of Events");
  

  Char_t *centName[4] = {"hCent_V0M","hCent_V0MEq","hCent_TRK","hCent_ZEMvsZDC"};

  for(Int_t icent = 0;icent<4;icent++){ 
    hCent[icent] = new TH1F(centName[icent],centName[icent],10,centmin,centmax);
    Char_t Xlabel[40];
    sprintf(Xlabel,"Centrality by %s",centName[icent]);
    //    hCent[icent]->SetXTitle("Centrality Percentile");
    hCent[icent]->SetXTitle(Xlabel);
    hCent[icent]->SetYTitle("Number of Events");
  }
  
  
  for(Int_t iasym = 0;iasym<20;iasym++){
    Char_t hname[30];
    sprintf(hname,"hNCh_ZDCAsymBin%d",iasym+1);
    hNCh[iasym] = new TH1F(hname,hname,200,0.0,2000.0);
    sprintf(hname,"hVertexZ_ZDCAsymBin%d",iasym+1);
    hVertexZ4ZDCAsym[iasym] = new TH1F(hname,hname,20,-10.0,10.0);
    sprintf(hname,"hTrackEta_ZDCAsymBin%d",iasym+1);
    hTrackEta4ZDCAsym[iasym] = new TH1F(hname,hname,20,-1.0,1.0);
  }

  hTrackEta = new TH1D("hTrackEta","hTrackEta",30,-1.5,1.5);
  hTrackEta->SetXTitle(" #eta of the Tracks");
  hTrackEta->SetYTitle("Number of Tracks");
    
  hTrackPhi = new TH1D("hTrackPhi","hTrackPhi",30,0.0,TMath::Pi());
  hTrackPhi->SetXTitle(" #Phi of the Tracks");
  hTrackPhi->SetYTitle("Number of Tracks");
    
  hTrackletEta = new TH1D("hTrackletEta","hTrackletEta",30,-1.5,1.5);
  hTrackletEta->SetXTitle(" #eta of the Tracklets");
  hTrackletEta->SetYTitle("Number of Tracklets");

  hTrackletPhi = new TH1D("hTrackletPhi","hTrackletPhi",30,0.0,TMath::Pi());
  hTrackletPhi->SetXTitle(" #Phi of the Tracklets");
  hTrackletPhi->SetYTitle("Number of Trackles");
    
  hZDCN1 = new TH1F("hZDCN1","hZDCN1",100,0.0,100000);
  hZDCN1->SetXTitle("neutron energy in ZDC1");
  hZDCN1->SetYTitle("Number of Events");
  
  hZDCN2 = new TH1F("hZDCN2","hZDCN2",100,0.0,100000);
  hZDCN2->SetXTitle("neutron energy in ZDC2");
  hZDCN2->SetYTitle("Number of Events");
  
  hV0A = new TH1F("hV0A","hV0A",500,0.0,25000.0);
  hV0A->SetXTitle("V0 A Multiplicity");
  hV0A->SetYTitle("Number of Events");

  hV0C = new TH1F("hV0C","hV0C",500,0.0,25000.0);
  hV0C->SetXTitle("V0 C Multiplicity");
  hV0C->SetYTitle("Number of Events");

  hV0AScaled = new TH1F("hV0AScaled","hV0AScaled",400,0.0,5.0);
  hV0AScaled->SetXTitle("V0 A Scaled Multiplicity");
  hV0AScaled->SetYTitle("Number of Events");

  hV0CScaled = new TH1F("hV0CScaled","hV0CScaled",400,0.0,5.0);
  hV0CScaled->SetXTitle("V0 C Scaled Multiplicity");
  hV0CScaled->SetYTitle("Number of Events");

  hZDCN1Scaled = new TH1F("hZDCN1Scaled","hZDCN1Scaled",100,0.0,5.0);
  hZDCN1Scaled->SetXTitle("Scaled neutron energy in ZDC1");
  hZDCN1Scaled->SetYTitle("Number of Events");
  
  hZDCN2Scaled = new TH1F("hZDCN2Scaled","hZDCN2Scaled",100,0.0,5.0);
  hZDCN2Scaled->SetXTitle("Scaled neutron energy in ZDC2");
  hZDCN2Scaled->SetYTitle("Number of Events");

  hZDCN1SigmaScaled = new TH1F("hZDCN1SigmaScaled","hZDCN1SigmaScaled",100,0.0,5.0);
  hZDCN1SigmaScaled->SetXTitle("Scaled neutron energy in ZDC1  Sigma Scaling applied");
  hZDCN1SigmaScaled->SetYTitle("Number of Events");
  
  hZDCN2SigmaScaled = new TH1F("hZDCN2SigmaScaled","hZDCN2SigmaScaled",100,0.0,5.0);
  hZDCN2SigmaScaled->SetXTitle("Scaled neutron energy in ZDC2 Sigma Scaling Applied");
  hZDCN2SigmaScaled->SetYTitle("Number of Events");
  
  hZDCEM1 = new TH1F("hZDCEM1","hZDCEM1",100,0.0,10000);
  hZDCEM1->SetXTitle("EM energy in ZDC1");
  hZDCEM1->SetYTitle("Number of Events");
  
  hZDCEM2 = new TH1F("hZDCEM2","hZDCEM2",100,0.0,10000);
  hZDCEM2->SetXTitle("EM energy in ZDC2");
  hZDCEM2->SetYTitle("Number of Events");

  hZDCEM1Scaled = new TH1F("hZDCEM1Scaled","hZDCEM1Scaled",100,0.0,10000);
  hZDCEM1Scaled->SetXTitle("EM energy in ZDC1 Scaled");
  hZDCEM1Scaled->SetYTitle("Number of Events");
  
  hZDCEM2Scaled = new TH1F("hZDCEM2Scaled","hZDCEM2Scaled",100,0.0,10000);
  hZDCEM2Scaled->SetXTitle("EM energy in ZDC2 Scaled");
  hZDCEM2Scaled->SetYTitle("Number of Events");

  hNTracks = new TH1F("hNTracks","hNTracks",350,0.0,7000);
  hNTracks->SetXTitle("Number of Tracks");
  hNTracks->SetYTitle("Number of Events");
  
  hNTracklets = new TH1F("hNTracklets","hNTracklets",225,0.0,4500);
  hNTracklets->SetXTitle("Number of Tracklets");
  hNTracklets->SetYTitle("Number of Events");

  hPosNTracks = new TH1F("hPosNTracks","hPosNTracks",350,0.0,3500);
  hPosNTracks->SetXTitle("Number of +ve eta tracks");
  hPosNTracks->SetYTitle("Number of Events");
  
  hNegNTracks = new TH1F("hNegNTracks","hNegNTracks",350,0.0,3500);
  hNegNTracks->SetXTitle("Number of -ve eta tracks");
  hNegNTracks->SetYTitle("Number of Events");
  
  hPospT = new TH1F("hPospT","hPospT",350,0.0,3500);
  hPospT->SetXTitle("Sum of pT of +ve eta tracks");
  hPospT->SetYTitle("Number of Events");
  
  hNegpT = new TH1F("hNegpT","hNegpT",350,0.0,3500);
  hNegpT->SetXTitle("Sum of pT of -ve eta tracks");
  hNegpT->SetYTitle("Number of Events");
  

  hPosNTracklets = new TH1F("hPosNTracklets","hPosNTracklets",200,0.0,2000);
  hPosNTracklets->SetXTitle("Number of +ve eta tracklets");
  hPosNTracklets->SetYTitle("Number of Events");
  
  hNegNTracklets = new TH1F("hNegNTracklets","hNegNTracklets",200,0.0,2000);
  hNegNTracklets->SetXTitle("Number of -ve eta tracklets");
  hNegNTracklets->SetYTitle("Number of Events");
  
  hpZDCAsym_pT = new TProfile("hpZDCAsym_pT","hpZDCAsym_pT",10,-1.0,1.0);
  hpZDCAsym_pT->SetXTitle("ZDCAsymmetry");
  hpZDCAsym_pT->SetYTitle("<mean pT>");
  
  hpZDCAsym_Cent = new TProfile("hpZDCAsym_Cent","hpZDCAsym_Cent",10,-1.0,1.0);
  hpZDCAsym_Cent->SetXTitle("ZDCAsymmetry");
  hpZDCAsym_Cent->SetYTitle("<mean Centrality>");
  
  hpZDCAsym_NCh = new TProfile("hpZDCAsym_NCh","hpZDCAsym_NCh",10,-1.0,1.0);
  hpZDCAsym_NCh->SetXTitle("ZDCAsymmetry");
  hpZDCAsym_NCh->SetYTitle("<mean NCh>");

  hpZDCAsym_Vz = new TProfile("hpZDCAsym_Vz","hpZDCAsym_Vz",10,-1.0,1.0);
  hpZDCAsym_Vz->SetXTitle("ZDCAsymmetry");
  hpZDCAsym_Vz->SetYTitle("<mean Vz>");
}

void BookAsymTree(Int_t runNumber){

  Char_t fname[120];
  sprintf(fname,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/AsymRun%dCent%3.3dto%3.3d.root",centmin,centmax,runNumber,centmin,centmax);
  f2 = new TFile(fname,"RECREATE");
  
  asymTree = new TTree("asymTree","asymTree");
  //  fRunNumber,fCentPercentile,fVertexX,fVertexY, fVertexZ, fzdcn1en,fzdcn2en,fNumberOfTracks,fNumberOfTracklets,fV0Mult
  asymTree->Branch("fRunNumber",&fRunNumber,"fRunNumber/I");
  asymTree->Branch("fNumberOfTracks",&fNumberOfTracks,"fNumberOfTracks/I");
  asymTree->Branch("fNumberOfTracklets",&fNumberOfTracklets,"fNumberOfTracklets/I");
  asymTree->Branch("fVertexX",&fVertexX,"fVertexX/F");
  asymTree->Branch("fVertexY",&fVertexY,"fVertexY/F");
  asymTree->Branch("fVertexZ",&fVertexZ,"fVertexZ/F");
  asymTree->Branch("fCentPercentile",fCentPercentile,"fCentPercentile[4]/F");
  asymTree->Branch("fzdcn1en",&fzdcn1en,"fzdcn1en/F");
  asymTree->Branch("fzdcn2en",&fzdcn2en,"fzdcn2en/F");
  asymTree->Branch("fzdcEM1en",&fzdcEM1en,"fzdcEM1en/F");
  asymTree->Branch("fzdcEM2en",&fzdcEM2en,"fzdcEM2en/F");
  asymTree->Branch("fV0Mult",fV0Mult,"fV0Mult[64]/F");
  

  asymTree->Branch("SpectAsym",&SpectAsym,"SpectAsym/F");
  asymTree->Branch("PartAsym",&PartAsym,"PartAsym/F");
  asymTree->Branch("NPartAsym",&NPartAsym,"NPartAsym/F");
  asymTree->Branch("NSpectAsym",&NSpectAsym,"NSpectAsym/F");
  asymTree->Branch("zdcasym",&zdcasym,"zdcasym/F");
  asymTree->Branch("zdcemasym",&zdcemasym,"zdcemasym/F");
  asymTree->Branch("trkasym",&trkasym,"Trkasym/F");
  asymTree->Branch("trkletasym",&trkletasym,"trkletasym/F");
  asymTree->Branch("v0asym",&v0asym,"v0asym/F");
  asymTree->Branch("pTasym",&pTasym,"pTasym/F");
  FBTree = new TTree("FBTree","FBTree");
  FBTree->Branch("nposTracks",&nposTracks,"nposTracks/I");
  FBTree->Branch("nnegTracks",&nnegTracks,"nnegTracks/I");
  //  FBTree->Branch("nposTracklets",&nposTracklets,"nposTracklets/I");
  //  FBTree->Branch("nnegTracklets",&nnegTracklets,"nnegTracklets/I");
  //  FBTree->Branch("nposTrksInDEta",nposTrksInDEta,"nposTrksInDEta[10]/I");
  //  FBTree->Branch("nnegTrksInDEta",nnegTrksInDEta,"nnegTrksInDEta[10]/I");
  
  return;
  
}

void GetAverages(Int_t runNumber){
 
  avzdcn1en = hZDCN1->GetMean();
  avzdcn2en = hZDCN2->GetMean();
  avzdcEM1en = hZDCEM1->GetMean();
  avzdcEM2en = hZDCEM2->GetMean();
  avV0AMult = hV0A->GetMean();
  avV0CMult = hV0C->GetMean();
  avPosNTracks = hPosNTracks->GetMean();
  avNegNTracks = hNegNTracks->GetMean();
  avPosNTracklets = hPosNTracklets->GetMean();
  avNegNTracklets = hNegNTracklets->GetMean();
  cout<<" avzdcn1en = "<<avzdcn1en<<" avzdcn2en = "<<avzdcn2en<<" avzdcEM1en = "<<avzdcEM1en<<" avzdcEM2en = "<<avzdcEM2en<<endl;
  cout<<" avPosNTracks = "<<avPosNTracks<<" avNegNTracks = "<<avNegNTracks<<" avPosNTracklets = "<<avPosNTracklets<<" avNegNTracklets = "<<avNegNTracklets<<endl;
  cout<<" avV0A ="<<avV0AMult<<" avV0C ="<<avV0CMult<<endl;

  Char_t fname_dat[120];
  sprintf(fname_dat,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/AsymRun%dCent%3.3dto%3.3d.dat",centmin,centmax,runNumber,centmin,centmax);
  cout<<" Creating file "<<fname_dat<<" for storing averages"<<endl; 
  ofstream fout(fname_dat,ios::out);
  fout<<runNumber<<" "<<hZDCN1->GetEntries()<<endl;
  fout<<avzdcn1en<<" "<<avzdcn2en<<endl;
  fout<<avPosNTracks<<" "<<avNegNTracks<<endl;
  fout<<avPosNTracklets<<" "<<avNegNTracklets<<endl;
  fout<<avV0AMult<<" "<<avV0CMult<<endl;
  return;
}


void GetAllAsymmetry(Int_t runNumber, Int_t nFile){
  
  hZDCN1SigmaScaled->Reset();
  hZDCN2SigmaScaled->Reset();
  hpZDCAsym_Vz->Reset();
  hpZDCAsym_Cent->Reset();
  hV0AScaled->Reset();
  hV0CScaled->Reset();
  hpZDCAsym_pT->Reset();
  hpZDCAsym_NCh->Reset();
  
  
  hNTracks1Scaled->Reset();
  hNTracks2Scaled->Reset();
  
  hNTracklets1Scaled->Reset();
  hNTracklets2Scaled->Reset();

  for(Int_t iasym = 0;iasym<20;iasym++){
    hTrackEta4ZDCAsym[iasym]->Reset();
    hVertexZ4ZDCAsym[iasym]->Reset();
    hNCh[iasym]->Reset();
  }
      
  // rereading the tree to get ZDC,tracklet and track asym and filling it in another tree 
  cout<<"<ZDCn1en> = "<<avzdcn1en<<" <ZDCn2en> = "<<avzdcn2en<<endl; 
  cout<<" <NposTracks>="<<avPosNTracks<<" <NNegTracks>="<<avNegNTracks<<endl; 
  cout<<" <NposTracklets>="<<avPosNTracklets<<" <NNegTracklets>="<<avNegNTracklets<<endl; 
  cout<<" avV0AMult ="<<avV0AMult<<" avV0CMult = "<<avV0CMult<<endl;

  cout<<" Going to get Asymmetries now"<<endl;
  for(Int_t ifile = 0;ifile<nFile;ifile++){

    Char_t fname[85];
    if(nfile==1){    
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i.root",runNumber);
    }
    if(nfile>1){
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i%c.root",runNumber,97+ifile);
    }
    cout<<" Reading file "<<fname<<" for Run "<<runNumber<<endl;
    
    f = new TFile(fname);
    tree = (TTree*)f->Get("fEventTree");
    BookTree(tree);
    nentries = tree->GetEntries();
    cout<<" Entries in the tree["<<ifile+1<<"] are "<<nentries<<endl;
    //    nentries = 10000;
    for(Int_t ient = 0;ient<nentries;ient++){
      if(ient%10000==0) cout<<" reading event number:"<<ient<<" from file["<<ifile+1<<"]"<<endl;
      
      zdcasym  =1.0, zdcemasym = 1.0, trkasym=1.0, trkletasym = 1.0 , v0asym =1.0, pTasym = 1.0;
      PartAsym = 1.0, SpectAsym = 1.0, NPartAsym = 1.0, NSpectAsym = 1.0;

      tree->GetEntry(ient);
      // Event Selection
      if(SelectEvent()==kFALSE) continue;
      // Get HIJING Asymmetry
      //Int_t fProjSpectP,fProjSpectN,fTargSpectP,fTargSpectN;
      /*
      NSpectAsym = (fProjSpectN*1.0 - fTargSpectN*1.0)/(fProjSpectN*1.0 + fTargSpectN*1.0);
      SpectAsym = (fProjSpectN*1.0 +fProjSpectP*1.0 - fTargSpectP*1.0 -fTargSpectN*1.0)/(fProjSpectN*1.0 + fTargSpectN*1.0 + fProjSpectP*1.0 + fTargSpectP*1.0);
      PartAsym = (-fProjSpectN*1.0 -fProjSpectP*1.0 +fTargSpectP*1.0 +fTargSpectN*1.0)/(208.0+208.0 -fProjSpectN*1.0 - fTargSpectN*1.0 - fProjSpectP*1.0 - fTargSpectP*1.0);
      NPartAsym = (-fProjSpectN*1.0 + fTargSpectN*1.0)/(126.0+126.0 -fProjSpectN*1.0 - fTargSpectN*1.0);
      */
      // Getting zdc asymmetry
      //AliZDCReconstructor.cxx:   (ZN1,ZP1) or (ZNC, ZPC) or RIGHT refers to side C (RB26)   
      //AliZDCReconstructor.cxx:   (ZN2,ZP2) or (ZNA, ZPA) or LEFT refers to side A (RB24) 
      // Side A is +ve Z axis side
      // Side C is -ve Z axis side

      if(Debug)     cout<<" input for scaling"<<fzdcn1en<<" "<<fzdcn1Mean<<" "<<zdcn1Sigma<<" "<<fzdcn2en<<" "<<fzdcn2Mean<<" "<<zdcn2Sigma<<endl;
      //      Float_t zdcn1Scaled = GetScaled(fzdcn1en/fzdcn1Mean, zdcn1Mean, zdcn2Sigma, zdcn1Sigma);
      //      Float_t zdcn2Scaled = fzdcn2en/fzdcn2Mean;
      Float_t zdcn2Scaled = GetScaled(fzdcn2en/fzdcn2Mean, zdcn2Mean, zdcn1Sigma, zdcn2Sigma);
      Float_t zdcn1Scaled = fzdcn1en/fzdcn1Mean;
      if(Debug)cout<<"zdc1 and 2 after scaling are = "<<zdcn1Scaled<<" "<<zdcn2Scaled<<endl;
      //      zdcasym = GetAsymmetry(fzdcn2en,avzdcn2en,fzdcn1en,avzdcn1en);
      zdcasym = GetAsymmetry(zdcn2Scaled,1.0,zdcn1Scaled,1.0);
      if(Debug)cout<<" zdcasym = "<<zdcasym<<endl;
      zdcemasym = GetAsymmetry(fzdcEM2en,avzdcEM2en,fzdcEM1en,avzdcEM1en);
      if(zdcasym>1.0 || zdcasym<-1.0) {
	cout<<"z2,avz2,z1,avz1="<<fzdcn2en<<","<<avzdcn2en<<","<<fzdcn1en<<","<<avzdcn1en<<endl;
	cout<<"zdcn2Scaled, zdcn1Scaled = asym "<<zdcn2Scaled<<" "<<zdcn1Scaled<<" = "<<zdcasym<<endl;
	continue;
      }
      Int_t iasym = (Int_t)((zdcasym+1.0)/0.1);
      hZDCN1SigmaScaled->Fill(zdcn1Scaled);
      hZDCN2SigmaScaled->Fill(zdcn2Scaled);
      
      hpZDCAsym_Vz->Fill(zdcasym,fVertexZ);
      hpZDCAsym_Cent->Fill(zdcasym,fCentPercentile[0]);
      
      //      cout<<" Obtained ZDC Asym"<<endl;

      // GetV0 Scaled mult and V0 asym
      Float_t V0AMult = 0, V0CMult = 0;
      for(Int_t iv = 0;iv<64;iv++){
	if(iv<32){
	  V0CMult = V0CMult + fV0Mult[iv];
	}else{
	  V0AMult = V0AMult + fV0Mult[iv];
	}
      }
      Float_t V0AMultScaled = V0AMult/avV0AMult;
      Float_t V0CMultScaled = V0CMult/avV0CMult;
      v0asym = (V0AMultScaled - V0CMultScaled)/(V0AMultScaled + V0CMultScaled);
      hV0AScaled->Fill(V0AMultScaled);
      hV0CScaled->Fill(V0CMultScaled);

      //      cout<<" Obtained V0 asym"<<endl;

      // Track loop
      nposTracks = 0, nnegTracks=0;
      for(Int_t ietagap  = 0;ietagap<netagap;ietagap++){
	nposTrksInDEta[ietagap] = 0; nnegTrksInDEta[ietagap] = 0;
      }
      Int_t NCh=0;
      Float_t SumpTPos = 0, SumpTNeg = 0;

      /*
Please see comment elsewhere inthis program on why this array is starting from 1 and not zero as it should be.
This involves a loss of first track but thatis the best we can do. Will be corrected in future passes of data.
       */
      //      cout<<" Number of tracks  ="<<fNumberOfTracks<<endl;
      for(Int_t itr = 1;itr<fNumberOfTracks;itr++){
	if(fTrackFilterMap[itr]==0) continue;

	if(BitCheck((UInt_t)fTrackFilterMap[itr],9)==kFALSE && BitCheck((UInt_t)fTrackFilterMap[itr],8)==kFALSE) continue;
	if(Debug) cout<<" TrackEta["<<itr<<"]="<<fTrackEta[itr]<<endl;
	if(fTrackEta[itr]>0.0 && fTrackEta[itr]<trketamax) {
	  nposTracks++;
	  SumpTPos+=fTrackPt[itr];
	  if(Debug) cout<<" This is "<<nposTracks<<"th +ve track"<<endl;
	}
	if(fTrackEta[itr]<0.0 && fTrackEta[itr]>trketamin) {
	  nnegTracks++;
	  SumpTNeg+=fTrackPt[itr];
	  if(Debug) cout<<" This is "<<nnegTracks<<"th -ve track"<<endl;
	}
	//	cout<<" TrackpT eta="<<fTrackPt[itr]<<" "<<fTrackEta[itr]<<" zdcasym="<<zdcasym<<" iasym="<<iasym<<endl;

	hpZDCAsym_pT->Fill(zdcasym,fTrackPt[itr]);
	if(iasym<20 && iasym>=0){
	  hTrackEta4ZDCAsym[iasym]->Fill(fTrackEta[itr]);
	}
	if(TMath::Abs(fTrackEta[itr])<0.5) {
	  NCh++;
	  //	  cout<<" fTrackEta ["<<NCh<<"]="<<fTrackEta[itr]<<" out of "<<fNumberOfTracks<<endl;
	  
	}
      }
      //      cout<<" Obtained Tracks +/-"<<nposTracks<<" "<<nnegTracks<<endl;
      if(nposTracks==0 || nnegTracks==0) {
	cout<<" Warning empty Tracks +/-"<<nposTracks<<" "<<nnegTracks<<endl;
	continue;
      }

      if(Debug)	cout<<nnegTracks<<" "<<nposTracks<<endl;
      
      hpZDCAsym_NCh->Fill(zdcasym,NCh);
      
      hNCh[iasym]->Fill(NCh);
      if(iasym<20 && iasym>=0){
	hVertexZ4ZDCAsym[iasym]->Fill(fVertexZ);
      }

      hNTracks1Scaled->Fill(nposTracks/avPosNTracks);
      hNTracks2Scaled->Fill(nnegTracks/avNegNTracks);
      trkasym = GetAsymmetry(nposTracks,avPosNTracks,nnegTracks,avNegNTracks);


      Float_t meanpTPos = SumpTPos/nposTracks;
      Float_t meanpTNeg = SumpTNeg/nnegTracks;
      pTasym = (meanpTPos-meanpTNeg)/(meanpTPos+meanpTNeg);

      //      cout<<" Obtained Track Asym"<<endl;
      // Tracklet loop
      nposTracklets = 0, nnegTracklets=0;
      for(Int_t itr = 0;itr<fNumberOfTracklets;itr++){
	if(fTrackletEta[itr]>0.0 && fTrackletEta[itr]<trkletetamax) nposTracklets++;
	if(fTrackletEta[itr]<0.0 && fTrackletEta[itr]>trkletetamin) nnegTracklets++;
      }
      if(nposTracklets==0 || nnegTracklets==0) {
	cout<<" Warning empty Tracklets +/-"<<nposTracklets<<" "<<nnegTracklets<<endl;
	continue;
      }
      hNTracklets1Scaled->Fill(nposTracklets/avPosNTracklets);
      hNTracklets2Scaled->Fill(nnegTracklets/avNegNTracklets);
      trkletasym = GetAsymmetry(nposTracklets,avPosNTracklets,nnegTracklets,avNegNTracklets);

      //    cout<<" ZDC Asym is ="<<zdcasym<<" trkasym = "<<trkasym<<" trkletasym = "<<trkletasym<<endl;

      asymTree->Fill();
      FBTree->Fill();
    }
    f->Close();
  }
  //Saving the asymmetry parameters in a tree in a different root file
  if(hZDCN1SigmaScaled->GetEntries()==0)  return;
  hZDCN1SigmaScaled->Fit("gaus");
  hZDCN2SigmaScaled->Fit("gaus");
  
  f2->cd();
  asymTree->Write();
  FBTree->Write();
  hZDCN1SigmaScaled->Write();
  hZDCN2SigmaScaled->Write();
  hNTracks1Scaled->Write();
  hNTracks2Scaled->Write();
  hNTracklets1Scaled->Write();
  hNTracklets2Scaled->Write();
  hV0AScaled->Write();
  hV0CScaled->Write();
  
  hpZDCAsym_pT->Write();
  hpZDCAsym_NCh->Write();
  hpZDCAsym_Vz->Write();
  hpZDCAsym_Cent->Write();
  
  for(Int_t iasym = 0;iasym<20;iasym++){
    hNCh[iasym]->Write();
    hVertexZ4ZDCAsym[iasym]->Write();
    hTrackEta4ZDCAsym[iasym]->Write();
  }
  f2->Close();
  return;
}


void FillScaledQAPlots(Int_t runNumber,Int_t nFile){

  // Reset the old plots since they contain the data of older runs which has olready been saved in their root files
  // One has to change the structure tsuch that these plots are not filled in the GetAsymmetry
  // Or one has to remove it from here and fill it only in get Asymmetry
  hZDCN1Scaled->Reset();
  hZDCN2Scaled->Reset();

  cout<<" I am in FillScaled"<<endl;
  //Obtain the scaling parameters from the old QA file for that run and that centrality here.
  // Open the file for saving the new histrograms in the same file. 
  Char_t fname_QA[120];
  sprintf(fname_QA,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/QARun%dCent%3.3dto%3.3d.root",centmin,centmax,runNumber,centmin,centmax);
  TFile * f3 = new TFile(fname_QA,"UPDATE");
  TH1F * hzdcn1 = (TH1F*)f3->Get("hZDCN1");
  hzdcn1->Fit("gaus");
  fzdcn1Mean = hzdcn1->GetFunction("gaus")->GetParameter(1);
  fzdcn1Sigma = hzdcn1->GetFunction("gaus")->GetParameter(2);
  
  TH1F * hzdcn2 = (TH1F*)f3->Get("hZDCN2");
  hzdcn2->Fit("gaus");
  fzdcn2Mean = hzdcn2->GetFunction("gaus")->GetParameter(1);
  fzdcn2Sigma = hzdcn2->GetFunction("gaus")->GetParameter(2);

  cout<<" fzdcn1Mean fzdcn1Sigma = "<<fzdcn1Mean<<" "<<fzdcn1Sigma<<endl;
  cout<<" fzdcn2Mean fzdcn2Sigma = "<<fzdcn2Mean<<" "<<fzdcn2Sigma<<endl;
  // Now reading the AOD file and applying the scaling.
  Int_t nSelectEvent = 0;  
  for(Int_t ifile = 0;ifile<nFile;ifile++){
    Char_t fname[85];
    if(nfile==1){    
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i.root",runNumber);
    }
    if(nfile>1){
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i%c.root",runNumber,97+ifile);
    }
    cout<<" Reading file "<<fname<<" for Run "<<runNumber<<endl;
    f = new TFile(fname);
    tree = (TTree*)f->Get("fEventTree");
    BookTree(tree);
    nentries = tree->GetEntries();
    cout<<" Entries in the tree are "<<nentries<<endl;
    
    for(Int_t ient = 0;ient<nentries;ient++){
      tree->GetEntry(ient);
      if(ient%10000==0) cout<<" reading event number :"<<ient<<" Number of Selected events is "<<nSelectEvent<<endl;
      // Event Selection
      if(SelectEvent()==kFALSE) continue;
      nSelectEvent++;

      //      cout<<" Simple scaling fzdcn1, fzdcn2 = "<<fzdcn1en/zdcn1Mean<<" "<< fzdcn2en/zdcn2Mean<<endl;
      
      hZDCN1Scaled->Fill(fzdcn1en/fzdcn1Mean);
      hZDCN2Scaled->Fill(fzdcn2en/fzdcn2Mean);
    }
    f->Close();
  }

  hZDCN1Scaled->Fit("gaus");
  zdcn1Mean = hZDCN1Scaled->GetFunction("gaus")->GetParameter(1);
  zdcn1Sigma = hZDCN1Scaled->GetFunction("gaus")->GetParameter(2);
  hZDCN2Scaled->Fit("gaus");
  zdcn2Mean = hZDCN2Scaled->GetFunction("gaus")->GetParameter(1);
  zdcn2Sigma = hZDCN2Scaled->GetFunction("gaus")->GetParameter(2);
  cout<<" zdcn1Mean zdcn1Sigma after simple scaling = "<<zdcn1Mean<<" "<<zdcn1Sigma<<endl;
  cout<<" zdcn2Mean zdcn2Sigma after simple scaling = "<<zdcn2Mean<<" "<<zdcn2Sigma<<endl;

  cout<<" Number of Events Selected = "<<nSelectEvent<<endl;
  
  f3->cd();
  hZDCN1Scaled->Write();
  hZDCN2Scaled->Write();
  f3->Close();
  
  return;
}

Float_t GetScaled(Float_t value, Float_t mean, Float_t sigma1, Float_t sigma2){
  // Here sigma1 is the value that we want to scale sigma2 to.
  // So after scaling the scaled distribution of 2 has the same sigma as the scaled distribution of 1
  // So the scaling of sigma is only applied to 2 and not to 1
  // 1 is only scaled by the mean.
  
  Float_t scaled = (value/mean)-1;
  scaled  = (scaled*sigma1/sigma2)+1;
  return scaled;
}

void FillQAPlots(Int_t runNumber,Int_t nFile){

  // Reset the old plots since they contain the data of older runs which has olready been saved in their root files
  hZDCN1->Reset();
  hZDCN2->Reset();
  hZDCEM1->Reset();
  hZDCEM2->Reset();
  hV0A->Reset();
  hV0C->Reset();
  hPosNTracks->Reset();
  hNegNTracks->Reset();
  hPosNTracklets->Reset();
  hNegNTracklets->Reset();
  hVertexXY->Reset();
  hVertexZ->Reset();
  hzdcn1vsn2->Reset();
  hNTracks->Reset();
  hNTracklets->Reset();
  hTrackEta->Reset();
  hTrackPhi->Reset();
  hTrackletEta->Reset();
  hTrackletPhi->Reset();
  for(Int_t icent=0;icent<4;icent++){
    hCent[icent]->Reset();
  }
  hNTracks_Cent->Reset();
  hNTracklets_Cent->Reset();

  Int_t nSelectEvent = 0;
  for(Int_t ifile = 0;ifile<nFile;ifile++){
    Char_t fname[85];
    if(nfile==1){    
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i.root",runNumber);
    }
    if(nfile>1){
      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/DATA/AOD.%i%c.root",runNumber,97+ifile);
    }
    cout<<" Reading file "<<fname<<" for Run "<<runNumber<<endl;
    
    f = new TFile(fname);
    tree = (TTree*)f->Get("fEventTree");
    BookTree(tree);
    nentries = tree->GetEntries();
    cout<<" Entries in the tree are "<<nentries<<endl;
    //    nentries = 10000;
    
    
    for(Int_t ient = 0;ient<nentries;ient++){
      tree->GetEntry(ient);
      if(ient%10000==0) cout<<" reading event number :"<<ient<<" Number of Selected events is "<<nSelectEvent<<endl;
      // Event Selection
      if(SelectEvent()==kFALSE) continue;
      nSelectEvent++;
      hVertexXY->Fill(fVertexX,fVertexY);
      hVertexZ->Fill(fVertexZ);
      hzdcn1vsn2->Fill(fzdcn1en,fzdcn2en);
      //      hpzdcn1vsn2->Fill(fzdcn1en,fzdcn2en);
      hNTracklets->Fill(fNumberOfTracklets);
      hNTracks->Fill(fNumberOfTracks);
      for(Int_t icent = 0;icent<4;icent++){
	hCent[icent]->Fill(fCentPercentile[icent]);
      }
      hZDCN1->Fill(fzdcn1en);
      hZDCN2->Fill(fzdcn2en);
      hZDCEM1->Fill(fzdcEM1en);
      hZDCEM2->Fill(fzdcEM2en);
      // loop over all V0A and V0C mult and add them
      Float_t V0AMult = 0, V0CMult = 0;
      for(Int_t iv = 0;iv<64;iv++){
	if(iv<32){
	  V0CMult = V0CMult + fV0Mult[iv];
	}else{
	  V0AMult = V0AMult + fV0Mult[iv];
	}
      }
      hV0A->Fill(V0AMult);
      hV0C->Fill(V0CMult);
      
      // Track loop to get mean +ve and -ve tracks and Tracklets
      
      
      nposTracks = 0; nnegTracks=0;
      if(fNumberOfTracks>NMaxTrack) cout<<" Warning .. Increase NMaxTrack array"<<endl;
      
      /*
	RR
	Because of a bug in the AliAnalysisFBMultFluct.cxx Program the first track is carrying junk value. 
	The first track ( itr ==0) is actually fed into itr=1 and so on. So the first value of all track arrays is junk.
	This has been corrected now and should be okay in future runs
      */
      //      for(Int_t itr = 0;itr<fNumberOfTracks;itr++){
      
      
      for(Int_t itr = 1;itr<fNumberOfTracks;itr++){
	//	if(Debug)
	//	  cout<<" fTrackFilterMap["<<itr<<"]="<<fTrackFilterMap[itr]<<" "<<(UInt_t)fTrackFilterMap[itr]<<"  ";
	// printf("%12b\n",(UInt_t)fTrackFilterMap[itr]);
	if(BitCheck((UInt_t)fTrackFilterMap[itr],9)==kFALSE && BitCheck((UInt_t)fTrackFilterMap[itr],8)==kFALSE) continue;
	//	cout<<" Got Track NUmber "<<itr<<" out of "<<fNumberOfTracks<<endl;
	if(fTrackEta[itr]>trketamin && fTrackEta[itr]<trketamax){
	  if(fTrackEta[itr]>0.0) nposTracks++;
	  if(fTrackEta[itr]<0.0 ) nnegTracks++;
	  hTrackEta->Fill(fTrackEta[itr]);
	  hTrackPhi->Fill(fTrackPhi[itr]);
	}
      }
      
      if(nposTracks==0 || nnegTracks==0) continue;
      hPosNTracks->Fill(nposTracks);
      hNegNTracks->Fill(nnegTracks);
      // Tracklet loop
      nposTracklets = 0, nnegTracklets=0;
      if(fNumberOfTracklets>NMaxTracklet) cout<<" Warning .. Increase NMaxTracklet array"<<endl;
      
      for(Int_t itr = 0;itr<fNumberOfTracklets;itr++){
	if(fTrackletEta[itr]>trkletetamin && fTrackletEta[itr]<trkletetamax){
	  if(fTrackletEta[itr]>0.0) nposTracklets++;
	  if(fTrackletEta[itr]<0.0 ) nnegTracklets++;
	  hTrackletEta->Fill(fTrackletEta[itr]);
	  hTrackletPhi->Fill(fTrackletPhi[itr]);
	}
      }
      hPosNTracklets->Fill(nposTracklets);
      hNegNTracklets->Fill(nnegTracklets);
      hNTracks_Cent->Fill(fCentPercentile[0],nposTracks+nnegTracks);
      hNTracklets_Cent->Fill(fCentPercentile[0],nposTracklets+nnegTracklets);
    }
    
  }
  cout<<" Number of Events Selected = "<<nSelectEvent<<endl;
  
  Char_t fname_QA[120];
  sprintf(fname_QA,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/QARun%dCent%3.3dto%3.3d.root",centmin,centmax,runNumber,centmin,centmax);
  TFile * f3 = new TFile(fname_QA,"RECREATE");
  f3->cd();
  hZDCN1->Write();
  hZDCN2->Write();
  hZDCEM1->Write();
  hZDCEM2->Write();
  hV0A->Write();
  hV0C->Write();
  hPosNTracks->Write();
  hNegNTracks->Write();
  hPosNTracklets->Write();
  hNegNTracklets->Write();
  hVertexXY->Write();
  hVertexZ->Write();
  hzdcn1vsn2->Write();
  //  hpzdcn1vsn2->Write();
  hNTracklets->Write();
  hNTracks->Write();
  for(Int_t icent=0;icent<4;icent++){
    hCent[icent]->Write();
  }

  hNTracks_Cent->Write();
  hNTracklets_Cent->Write();
  hTrackEta->Write();
  hTrackPhi->Write();
  hTrackletEta->Write();
  hTrackletPhi->Write();
  f3->Close();
  
  return;
}


void InitRun(){
  
  
  avzdcn1en=1.0;avzdcn2en=1.0; avzdcEM1en=1.0;avzdcEM2en=1.0;
  avV0AMult=1.0; avV0CMult=1.0;
  zdcasym  =1.0; zdcemasym = 1.0; trkasym=1.0; trkletasym = 1.0;

  v0asym = 1.0; pTasym = 1.0;
  NEvCut_Cent=0; NEvCut_VertexZ=0; NEvCut_VertexXY=0; NEvCut_NTrack=0;
  NEvCut_NTracklet=0;
  NEvTotal=0; NEvAccepted=0 ; NEvCut_ZDC=0;
  avPosNTracks=1.0; avNegNTracks=1.0; avPosNTracklets=1.0; avNegNTracklets=1.0;
}

void DrawQAPlots(){
  cout<<" Drawing QA plots"<<endl;

  TCanvas * cQA = new TCanvas("cQA","cQA",1000,700);
  cQA->Divide(3,3);
  cQA->cd(1);
  hVertexXY->Draw();
  cQA->cd(2);
  hVertexZ->Draw();
  cQA->cd(3);
  cQA->GetPad(3)->SetLogy();
  hNTracklets->Draw();
  cQA->cd(4);
  hCent[0]->Draw();
  cQA->cd(5);
  hCent[1]->Draw();
  cQA->cd(6);
  hCent[2]->Draw();
  cQA->cd(7);
  hCent[3]->Draw();
  cQA->cd(8);
  hzdcn1vsn2->Draw("Box");
  cQA->cd(9);
  cQA->GetPad(9)->SetLogy();
  hNTracks->Draw();

  /*  
  TCanvas * cMult = new TCanvas("cMult","cMult",800,500);
  cMult->Divide(2,1);
  cMult->cd(1);
  cMult->GetPad(1)->SetLogy();
  hPosNTracks->Draw();
  hNegNTracks->SetLineColor(2);
  hNegNTracks->Draw("same");
  cMult->cd(2);
  cMult->GetPad(2)->SetLogy();
  hPosNTracklets->Draw();
  hNegNTracklets->SetLineColor(2);
  hNegNTracklets->Draw("same");
  */  
  
}


Bool_t SelectEvent(){
  
  NEvTotal++;
  if(fCentPercentile[0] <= centmin || fCentPercentile[0] > centmax) {
    NEvCut_Cent++;
    //    cout<<" fCentPercentile Cut = "<<fCentPercentile[0]<<endl;
    return kFALSE;
  }
  if(fVertexZ<vzmin || fVertexZ>vzmax ){
    NEvCut_VertexZ++;
    return kFALSE;
  }
  if(fVertexX<vxymin || fVertexX>vxymax || fVertexY<vxymin || fVertexY>vxymax ){
    NEvCut_VertexXY++;
    return kFALSE;
  }
  
  if(fNumberOfTracks<=0 || fNumberOfTracklets<=0){
    NEvCut_NTrack++;
    return kFALSE;
  }
  
  if(fzdcn1en<=0 || fzdcn2en<=0) {
    //    cout<<" Empty ZDC neutrons : "<<fzdcn1en<<" "<<fzdcn2en<<endl;
    NEvCut_ZDC++;
    return kFALSE;
  }
  if(fzdcEM1en<=0 || fzdcEM2en<=0) {
    //    cout<<" Empty ZDC EMCalorimeter : "<<fzdcEM1en<<" "<<fzdcEM2en<<endl;
    NEvCut_ZDC++;
    return kFALSE;
  }
  
  return kTRUE;
  
}


// This program checks the bitnary value for each power of two in a given number. 
// The first number s the given number and the second is the position ( i.e. bit number)
// in binary which needs to be checked
Bool_t  BitCheck(UInt_t filterbit,Int_t checkbit ){
  //  Int_t two[14]={0};
  Int_t powerOf2, remainder;
  for(Int_t i =13;i>=0;i--){
    powerOf2=pow(2,i);
    if(filterbit>=powerOf2) {
      //      two[i]=1;
      filterbit = filterbit-powerOf2;
      //      cout<<" I have "<<i<<" flag"<<" filterbit = "<<filterbit<<" power of 2 ="<<powerOf2<<endl;
      if(checkbit==i) {
	//	cout<<filterbit<<" I have filterbit "<<i<<endl;
	return kTRUE;
      }
      if (filterbit==0) return kFALSE;
    }
  }
}
