#include "TGeoManager.h"

/*
The user needs to provide a vector int with the list of bad channels (by tower ID).
This macro will create one rootfile for each desired list of badchannels. These root files can be used
with CreateEMCAL_OADB_BadChannels, in order to create the OADB file.
*/

// ******* Create Histograms for BadChannels according to the EMCAL SM's

void CreateEMCAL_BadChannelsHistos()
{
gSystem->Load("libOADB");  
gSystem->Load("libEMCALbase");
gSystem->Load("libEMCALUtils");
gSystem->Load("libEMCALrec");

AliEMCALGeometry   *fEMCALGeo=new AliEMCALGeometry();               //! EMCAL geometry
AliEMCALRecoUtils  *fEMCALRecoUtils=new AliEMCALRecoUtils();         //! Pointer to EMCAL utilities for clusterization
TGeoManager::Import("geometry.root");
fEMCALGeo =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

// *********** List of Bad Channels *****

// LHC11a (2.76 pp) pass1 from EMCALTenderSupply list (September 3rd 2011) AliRootTrunk revision 51405
const Int_t nTowers11a1=89;
Int_t hotChannels11a1[nTowers11a1]={74, 103, 152, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 368, 369, 370, 371, 372, 373, 374,375, 376, 377, 378, 379, 380, 381, 382, 383, 917, 1275, 1288, 1519, 1595, 1860, 1967, 2022, 2026, 2047, 2117, 2298, 2540, 2776, 3135, 3764, 6095, 6111, 6481, 6592, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7371, 7425, 7430, 7457, 7491, 7709, 8352, 8353, 8356, 8357, 8808, 8810, 8812, 8814, 9056, 9769, 9815, 9837};
Write_histo(nTowers11a1,hotChannels11a1,"11a_pass1",fEMCALGeo,fEMCALRecoUtils);

// LHC11a (2.76 pp) pass2 from EMCALTenderSupply list (September 3rd 2011) AliRootTrunk revision 51405
const Int_t nTowers11a2=24;
Int_t hotChannels11a2[nTowers11a2]= {74, 103, 152, 917, 1059, 1175, 1276, 1288, 1376, 1382, 1595, 2022, 2026, 2210, 2540, 2778, 2793, 3135, 3764, 5767, 6481, 7371, 7878, 9769};
Write_histo(nTowers11a2,hotChannels11a2,"11a_pass2",fEMCALGeo,fEMCALRecoUtils);

// LHC11c (7 TeV pp) from EMCALTenderSupply list (September 3rd 2011) AliRootTrunk revision 51405
const Int_t nTowers11c=231;
Int_t hotChannels11c[nTowers11c]={74, 103, 152, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 368, 369, 370, 371, 372, 373, 374,375, 376, 377, 378, 379, 380, 381, 382, 383, 917, 1059, 1160, 1263, 1275, 1276, 1288, 1376, 1382, 1384, 1414,1519, 1595, 1860, 1720, 1912, 1961, 1967,  2016, 2017, 2019, 2022, 2023, 2026, 2027, 2047, 2065, 2071, 2074, 2079, 2112, 2114, 2115, 2116, 2117, 2120, 2123, 2145, 2160, 2190, 2298, 2506, 2540, 2776, 2778, 3135, 3544, 3567, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3675, 3676, 3677, 3678, 3679, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725, 3726, 3727, 3764, 4026, 4121, 4157, 4208, 4209, 4568, 5560, 5767, 5969, 6065, 6076, 6084, 6087, 6095, 6111, 6128, 6129, 6130, 6131, 6133, 6136, 6137, 6138, 6139, 6140, 6141, 6142, 6143, 6340, 6369, 6425, 6481, 6561, 6592, 6678, 6907, 6925, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7089, 7371, 7375, 7425, 7457, 7430, 7491, 7709, 7874,  8273, 8352, 8353, 8354, 8356, 8357, 8362, 8808, 8769, 8810, 8811, 8812, 8814, 9056, 9217, 9302, 9365, 9389, 9815, 9769, 9833, 9837, 9850, 9895, 9997, 10082, 10086, 10087, 10091, 10095, 10112, 10113, 10114, 10115, 10116, 10117, 10118, 10119, 10120, 10121, 10122, 10123, 10576, 10718, 10723, 10918, 10919, 10922, 10925, 10926, 10927, 11276, 1136};
Write_histo(nTowers11c,hotChannels11c,"11c",fEMCALGeo,fEMCALRecoUtils);
 
}

void Write_histo(const Int_t nTowers,Int_t hotChannels[nTowers],const char *name="pass1",AliEMCALGeometry   *fEMCALGeo,AliEMCALRecoUtils  *fEMCALRecoUtils){

TH2I *hBadChannels[10];  
char fileName11[50];
char nameSM[50];  
sprintf(fileName11,"BadChannels2011_%s.root",name);

for(int i=0;i<10;i++){
	  sprintf(nameSM,"EMCALBadChannelMap_Mod%d",i);
	  hBadChannels[i]=new TH2I(nameSM,nameSM,48,0,48,24,0,24);
	  }
TFile *f=new TFile(fileName11,"RECREATE");
cout<<endl<<"********** "<<name<<" ********** "<<endl;   
  
Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;  
for(Int_t i=0; i<nTowers; i++){
	fEMCALGeo->GetCellIndex(hotChannels[i],nSupMod,nModule,nIphi,nIeta);
	fEMCALGeo->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
	fEMCALRecoUtils->SetEMCALChannelStatus(nSupMod, ieta, iphi);
	cout<<endl<<"NSupMod: "<<nSupMod<<" ieta: "<<ieta<<" iphi:"<<iphi<<endl; 
	hBadChannels[nSupMod]->SetBinContent(ieta,iphi,1);
	}
for(int i=0;i<10;i++){
  hBadChannels[i]->Write();  
  delete hBadChannels[i];
}

f->Close();
  
}