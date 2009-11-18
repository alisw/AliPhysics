//
// Macro adds ether ladder, module or channel to the SSD bad channels map.
// Input parameter: file name with original SSD bad channels map (optional),  
//                  equipment identifier (DDL),
//                  Slot number, ADC number (optional) and channel number (optional)
// Author: Oleksandr.Borysov@cern.ch
//


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <time.h>
#include <Riostream.h>
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TString.h"
#include "AliLog.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSSDv1.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSModuleDaSSD.h"
#include "AliITSHandleDaSSD.h"
#endif


#define EQUIPMENTTODDLMASK        0xFF
#define NUMBEROFSSDMODULESPERSLOT 12
#define MAXSSDDDLID               15
#define MINSSDMODULEID            500

Bool_t readBCMap(const Char_t *filename, AliITSBadChannelsSSDv2*& bcl);
Int_t addModuleBCMap(AliITSBadChannelsSSDv2 *bcl, const Int_t EqId, const Int_t SlotId, const Int_t adc);
Int_t addLadderBCMap(AliITSBadChannelsSSDv2 *bcl, const Int_t EqId, const Int_t SlotId);
Int_t  addChannelBCMap(AliITSBadChannelsSSDv2 *bcl, const Int_t EqId, const Int_t SlotId, const Int_t adc, const Int_t strp);
void  Usage (void);
void drawBadChannelMapDAQDB(const char* filename);
 
//class gives an access to the protected array with SSD_DDL_Map in case the geometry is not initialized
class AliITSRawStreamSSDtmp : public AliITSRawStreamSSDv1 {
public:
  AliITSRawStreamSSDtmp(AliRawReader* rawReader): AliITSRawStreamSSDv1(rawReader) {
    Setv11HybridDDLMapping();
    AliInfo("Using SSD DDL Map initialized by AliITSRawStreamSSDv1::Setv11HybridDDLMapping()");
  }
  virtual ~AliITSRawStreamSSDtmp() {};
  Int_t GetModId(Int_t iDDL, Int_t iModule) {return fgkDDLModuleMap[iDDL][iModule];}
  ClassDef(AliITSRawStreamSSDtmp, 0)
};

//_________________________________________________________//
Int_t patchSSDBadChannelsMap(const Char_t *fname = 0, 
			     const Int_t EqId = -1, 
			     const Int_t SlotId = -1, 
			     const Int_t adc = -1, 
			     const Int_t strp = -1) {
  //Macro to patch the bad channel list
  TString    bcfname, pathstr;
  AliITSBadChannelsSSDv2  *bc = 0;
  if (EqId < 0) { Usage(); cerr << "Equipment number (that is DDL) must be specified! Exit.\n"; return 0; } 
  else if ((EqId & EQUIPMENTTODDLMASK) > MAXSSDDDLID) { 
    cerr << "Icorrect equipment number (that is DDL)! Exit.\n"; 
    return -1;
  }  
  if (SlotId < 1 || SlotId > 9) { 
    cerr << "Slot number must be specified (in the range 1 - 9)! Exit.\n"; 
	Usage();   
    return -2; 
  }
  if (!fname || fname[0]==0) {
    cout << "File name with bad channels map is not specified, an empty map is used!\n";
    bc = new AliITSBadChannelsSSDv2();
    if (!bc) { cerr << "Error creating the AliITSBadChannelsSSDv2 object! Exit.\n"; return -1; }
    pathstr = "";
  }
  else {
	pathstr.Form(fname);
    if (!readBCMap(fname, bc)) {
      cerr << "Error reading file " << fname << " with Static Bad Channels Map! Exit.\n";
      return -3;
    }
  }  
  if (adc < 0) addLadderBCMap(bc, EqId, SlotId);
  else if ((adc < 6) || (adc > 7 && (adc - 2) < NUMBEROFSSDMODULESPERSLOT) ) {
    if (strp < 0) addModuleBCMap(bc, EqId, SlotId, adc);
    else if (strp < AliITSModuleDaSSD::GetStripsPerModuleConst()) addChannelBCMap(bc, EqId, SlotId, adc, strp);
         else {
	        cerr << "Incorrect number for Strip. Exit\n";
	        Usage(); 
            if (bc) delete bc;
            return -5;
        }
  }
  else {
    cerr << "Incorrect number for ADC. Exit\n";
    if (bc) delete bc;
    return -4;
  }
  bcfname = pathstr(0, pathstr.Last('/')+1);
  bcfname.Append(Form("ssdbcmap_%i.root", time(NULL)));
  TFile *bcfile = new TFile (bcfname.Data(),"RECREATE");
  if (bcfile->IsZombie()) {
    cerr << "Error open file " << bcfname.Data() << " for writing new bad channels map!\n";
    if (bc) delete bc;
    return -1;
  }
  bcfile->WriteTObject(bc);
  bcfile->Close();
  delete bcfile;
  cout << "New SSD bad channels map was saved in file " << bcfname.Data() << endl;
  if (bc) delete bc;
  return 0;
}

//_________________________________________________________//
Bool_t readBCMap(const Char_t *filename, AliITSBadChannelsSSDv2*& bcl) {
  // Reads Static Bad Channels Map from the file
  TFile *bcfile;
  if (!filename) {
    cout << "No file name is specified for Static Bad Channels Map!\n";
    return kFALSE;
  } 
  cout << "Reading SSD Bad Channels Map from the file " << filename << endl;
  bcfile = new TFile(filename, "READ");
  if (bcfile->IsZombie()) {
    cerr << "Error reading file " << filename << " with Static Bad Channels Map!\n";
    return kFALSE;
  }
  bcfile->GetObject("AliITSBadChannelsSSDv2;1", bcl);
  if (!bcl) {
    cout << "Error bcl == NULL!\n";
    bcfile->Close();
    delete bcfile;
    return kFALSE;
  }
  bcfile->Close();
  delete bcfile;
  return kTRUE;
}

//_________________________________________________________//
Int_t addModuleBCMap(AliITSBadChannelsSSDv2 *bcl, 
		     const Int_t EqId, const Int_t SlotId, const Int_t adc) {
  // Add module to bad channels map.
  const Char_t     isbad = 3;
  Int_t            ddl, mn, modid;
  AliRawReaderRoot       *rwr = 0;
  AliITSRawStreamSSDtmp  *rst = 0;
  rwr = new AliRawReaderRoot();
  rst = new AliITSRawStreamSSDtmp(rwr);
  ddl = EqId & EQUIPMENTTODDLMASK;
  mn = (SlotId - 1) * NUMBEROFSSDMODULESPERSLOT + (adc<8 ? adc : adc-2);
  modid = rst->GetModId(ddl, mn);
  modid -= MINSSDMODULEID;
  if (modid < 0) return 0;
  for (Int_t strind = 0; strind < AliITSModuleDaSSD::GetPNStripsPerModule(); strind++) {
      bcl->AddBadChannelP(modid, strind, isbad);
      bcl->AddBadChannelN(modid, strind, isbad);
  }
  delete rst;
  delete rwr;
  cout << "Module added:  ModId = " << modid + MINSSDMODULEID << "; Ddl/Ad/Adc: " << ddl 
       << " / " << SlotId << " / " << adc << endl;
  return 0;
}

//_________________________________________________________//
Int_t addLadderBCMap(AliITSBadChannelsSSDv2 *bcl, 
		     const Int_t EqId, const Int_t SlotId) {
  // Add ladder to bad channels map.
  for (Int_t adc = 0; adc < NUMBEROFSSDMODULESPERSLOT; adc++)
    addModuleBCMap(bcl, EqId, SlotId, (adc<6 ? adc : adc+2));
  return 0;
}

//_________________________________________________________//
Int_t  addChannelBCMap(AliITSBadChannelsSSDv2 *bcl, 
		       const Int_t EqId, const Int_t SlotId, 
                       const Int_t adc, const Int_t strp) {
  // Add strip to bad channels map.
  const Char_t     isbad = 3;
  Int_t            ddl, mn, modid;         
  AliRawReaderRoot       *rwr = 0;
  AliITSRawStreamSSDtmp  *rst = 0;
  rwr = new AliRawReaderRoot();
  rst = new AliITSRawStreamSSDtmp(rwr);
  ddl = EqId & EQUIPMENTTODDLMASK;
  mn = (SlotId - 1) * NUMBEROFSSDMODULESPERSLOT + (adc<8 ? adc : adc-2);
  modid = rst->GetModId(ddl, mn);  
  modid -= MINSSDMODULEID;
  if (modid < 0) { cout << "There is no module with given Equipment, Slot, adc.\n" ; return 0; }
  if (strp < AliITSModuleDaSSD::GetPNStripsPerModule() ) {
    bcl->AddBadChannelP(modid, strp, isbad);
  } else {
    bcl->AddBadChannelN(modid, (AliITSChannelDaSSD::GetMaxStripIdConst() - strp), isbad);
  }
  delete rst;
  delete rwr;
  cout << "Channel added (ModId/Ddl/Ad/Adc/Strip): " << modid + MINSSDMODULEID << " / " << ddl << " / " << SlotId 
       << " / " << adc << " / " << strp << endl;
  return 0;
}

//_________________________________________________________//
void  Usage (void) { 
  //Usage function
  cout << "Usage: PatchSSDBadChannelsMap(bc_fname /* can be \"\" */, EqipmentId, SlotId, adc /*optional*/, strip /*optional*/)\n"; 
}

//_______________________________________//
void drawBadChannelMapDAQDB(const char* filename) {
  //Draws the 2D plots of the bad channels maps
  const Int_t fgkSSDMODULES = 1698;
  static const Int_t fgkDefaultNStripsSSD = 768;
  gStyle->SetPalette(1,0);

  TH2F *fHistDAQDBPSideBadChannelMapLayer5 = new TH2F("fHistDAQDBPSideBadChannelMapLayer5",
						      "Layer 5;N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
  fHistDAQDBPSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistDAQDBPSideBadChannelMapLayer5->SetStats(kFALSE);
  fHistDAQDBPSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistDAQDBPSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistDAQDBPSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistDAQDBPSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistDAQDBPSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistDAQDBPSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistDAQDBPSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (p-side)[%]");

  TH2F *fHistDAQDBNSideBadChannelMapLayer5 = new TH2F("fHistDAQDBNSideBadChannelMapLayer5",
						      "Layer 5;N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
  fHistDAQDBNSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistDAQDBNSideBadChannelMapLayer5->SetStats(kFALSE);
  fHistDAQDBNSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistDAQDBNSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistDAQDBNSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistDAQDBNSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistDAQDBNSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistDAQDBNSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistDAQDBNSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
  
  TH2F *fHistDAQDBPSideBadChannelMapLayer6 = new TH2F("fHistDAQDBPSideBadChannelMapLayer6",
						      "Layer 6;N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
  fHistDAQDBPSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistDAQDBPSideBadChannelMapLayer6->SetStats(kFALSE);
  fHistDAQDBPSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistDAQDBPSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistDAQDBPSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistDAQDBPSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistDAQDBPSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistDAQDBPSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistDAQDBPSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (p-side)[%]");

  TH2F *fHistDAQDBNSideBadChannelMapLayer6 = new TH2F("fHistDAQDBNSideBadChannelMapLayer6",
						      "Layer 6;N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
  fHistDAQDBNSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistDAQDBNSideBadChannelMapLayer6->SetStats(kFALSE);
  fHistDAQDBNSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistDAQDBNSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistDAQDBNSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistDAQDBNSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistDAQDBNSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistDAQDBNSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistDAQDBNSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (n-side)[%]");

  //===============================//
  TFile *f = TFile::Open(filename);
  if(!f) {
    Printf("File poiter not valid");
    return;
  }

  if(!f->IsOpen()) {
    Printf("The file was not found");
    return;
  }
  //===============================//

  AliITSBadChannelsSSDv2 *badChannelsSSD = new AliITSBadChannelsSSDv2();
  badChannelsSSD = dynamic_cast<AliITSBadChannelsSSDv2 *>(f->Get("AliITSBadChannelsSSDv2"));

  //_____________________________________________________________________________//
  Int_t nPSideChannelsTotal = 0, nNSideChannelsTotal = 0;
  Int_t nBadPSideChannelsTotal = 0, nBadNSideChannelsTotal = 0;
  Int_t nBadPSideChannels = 0, nBadNSideChannels = 0;
  Int_t layer = 0, ladder = 0, module = 0;
  Int_t nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
  Int_t nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;
  //_____________________________________________________________________________//

  for(Int_t i = 0; i < fgkSSDMODULES; i++) {
    //for(Int_t i = 0; i < 1; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
    nBadPSideChannels = 0, nBadNSideChannels = 0;
    nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
    nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;

    Int_t badChannel = 0;
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      badChannel = (Int_t)(badChannelsSSD->GetBadChannelP(i,j));
      //cout<<"Module: "<<i+500<< " Strip: "<<j<<" - "<<badChannel<<endl;
      if(badChannel != 0) {
	if(layer == 5)
	  nPSideChannelsLayer5 += 1;
	if(layer == 6)
	  nPSideChannelsLayer6 += 1;
	nBadPSideChannels += 1;
      }
      badChannel = (Int_t)(badChannelsSSD->GetBadChannelN(i,j));
      //cout<<"Module: "<<i+500<< " Strip: "<<fgkDefaultNStripsSSD+j+1<<" - "<<badChannel<<endl;
      if(badChannel != 0) {
	if(layer == 5)                                                    
	  nNSideChannelsLayer5 += 1;
	if(layer == 6)
	  nNSideChannelsLayer6 += 1;
	nBadNSideChannels += 1;
      }
    }
    if(layer == 5) {
      if(nPSideChannelsLayer5 > 0)
	fHistDAQDBPSideBadChannelMapLayer5->Fill(module,499+ladder,
						 100.*nPSideChannelsLayer5/fgkDefaultNStripsSSD);
      else fHistDAQDBPSideBadChannelMapLayer5->Fill(module,499+ladder,0.0001);
      if(nNSideChannelsLayer5 > 0)
	fHistDAQDBNSideBadChannelMapLayer5->Fill(module,499+ladder,
						 100.*nNSideChannelsLayer5/fgkDefaultNStripsSSD);
      else fHistDAQDBNSideBadChannelMapLayer5->Fill(module,499+ladder,0.0001);
    }//layer 5
    if(layer == 6) {
      if(nPSideChannelsLayer6 > 0) 
	fHistDAQDBPSideBadChannelMapLayer6->Fill(module,599+ladder,
						 100.*nPSideChannelsLayer6/fgkDefaultNStripsSSD);
      else fHistDAQDBPSideBadChannelMapLayer6->Fill(module,599+ladder,0.0001);
      if(nNSideChannelsLayer6 > 0) 
	fHistDAQDBNSideBadChannelMapLayer6->Fill(module,599+ladder,
						 100.*nNSideChannelsLayer6/fgkDefaultNStripsSSD);
      else fHistDAQDBNSideBadChannelMapLayer6->Fill(module,599+ladder,0.0001);
    }//layer 6
      
    nBadPSideChannelsTotal += nBadPSideChannels;
    nBadNSideChannelsTotal += nBadNSideChannels;
    nPSideChannelsTotal += fgkDefaultNStripsSSD;
    nNSideChannelsTotal += fgkDefaultNStripsSSD;
  }

  cout<<"================================="<<endl;
  cout<<"Bad p-Side channels: "<<100.*nBadPSideChannelsTotal/nPSideChannelsTotal<<endl;
  cout<<"Bad n-Side channels: "<<100.*nBadNSideChannelsTotal/nNSideChannelsTotal<<endl;
  cout<<"================================="<<endl;

  TCanvas *cBadChannelDAQDB = new TCanvas("cBadChannelDAQDB",
					  "Bad channel list - DAQ DB",
					  0,0,900,900);
  cBadChannelDAQDB->SetHighLightColor(10); cBadChannelDAQDB->SetFillColor(10); 
  cBadChannelDAQDB->Divide(2,2);

  cBadChannelDAQDB->cd(1)->SetBottomMargin(.2); 
  cBadChannelDAQDB->cd(1)->SetLeftMargin(.15);
  cBadChannelDAQDB->cd(1)->SetRightMargin(.2);
  cBadChannelDAQDB->cd(1)->SetGridx(); cBadChannelDAQDB->cd(1)->SetGridy();
  cBadChannelDAQDB->cd(1); fHistDAQDBPSideBadChannelMapLayer5->Draw("colz"); 
  cBadChannelDAQDB->cd(2)->SetBottomMargin(.2); 
  cBadChannelDAQDB->cd(2)->SetLeftMargin(.15);
  cBadChannelDAQDB->cd(2)->SetRightMargin(.2);
  cBadChannelDAQDB->cd(2)->SetGridx(); cBadChannelDAQDB->cd(2)->SetGridy();
  cBadChannelDAQDB->cd(2); fHistDAQDBPSideBadChannelMapLayer6->Draw("colz");
  cBadChannelDAQDB->cd(3)->SetBottomMargin(.2); 
  cBadChannelDAQDB->cd(3)->SetLeftMargin(.15);
  cBadChannelDAQDB->cd(3)->SetRightMargin(.2);
  cBadChannelDAQDB->cd(3)->SetGridx(); cBadChannelDAQDB->cd(3)->SetGridy();
  cBadChannelDAQDB->cd(3); fHistDAQDBNSideBadChannelMapLayer5->Draw("colz"); 
  cBadChannelDAQDB->cd(4)->SetBottomMargin(.2); 
  cBadChannelDAQDB->cd(4)->SetLeftMargin(.15);
  cBadChannelDAQDB->cd(4)->SetRightMargin(.2);
  cBadChannelDAQDB->cd(4)->SetGridx(); cBadChannelDAQDB->cd(4)->SetGridy();
  cBadChannelDAQDB->cd(4); fHistDAQDBNSideBadChannelMapLayer6->Draw("colz");

  return;
}
