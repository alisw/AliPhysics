#ifndef __CINT__
#include <AliSysInfo.h>
#include <TFile.h>
#include <Riostream.h>
#include <AliDataFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <AliOADBContainer.h>
#endif

void testContainerOrig()
{
  const TString namecontainer = "AliEMCALRunDepTempCalibCorrections";
  //const TString namecontainer = "AliEMCALBadChannels";
  TFile *reader = TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data());
  //TFile *reader = TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data());
  AliOADBContainer *cont2 = static_cast<AliOADBContainer *>(reader->Get(namecontainer.Data()));
  cont2->Print();
  //cont2->GetDefaultList()->Delete();
  //cont2->WriteToFile("t.root");
  cont2->SetOwner(0);
  TObject *o1=cont2->GetObject(129983);
  o1->Print();
  AliOADBObjCache *c=cont2->GetObjectCache(129983);
  delete cont2;

  c->GetObject(129983)->Print();
  TObject *o2=c->GetObject(129983);
  if (o2) 
    o2->Print();
}

void testContainerOrig2() {
  // test implementation:
  // open 10 times the calibration entry for the temperature calibration
  // 2 ways: 
  //   1st: use AliOADBContainer::InitFromFile
  //   2nd: Open container by hand (with TFile::Open) and delete container by hand

  int N=100;
  int globalcounter = 0;
  const TString namecontainer = "AliEMCALRunDepTempCalibCorrections";
  if (1) {
    AliSysInfo::AddStamp("start", globalcounter++);
    std::cout << "Starting test InitFromFile" << std::endl;
    for(int i = 0; i < N; i++) {
      gSystem->Sleep(1000);
      std::cout << "Test instance: " << i << std::endl;
      AliOADBContainer cont(namecontainer.Data());
      cont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data(), namecontainer);
      cont.SetOwner(0);
      AliSysInfo::AddStamp("test_initfromfile", globalcounter++, i);
      // Close all files
      for(TIter fileiter = TIter(gROOT->GetListOfFiles()).Begin(); fileiter != TIter::End(); ++fileiter){
	TFile *closefile = static_cast<TFile *>(*fileiter);
	cout << "File name " << closefile->GetName() << endl;
	if(closefile){
	  closefile->Close();
	  delete closefile;
	} 
      }
      std::cout << "End test InitFromFile" << std::endl;
    }
    return;
  }
  if (1) {
    std::cout << "Starting test Manual" << std::endl;
    AliSysInfo::AddStamp("start_test_manual", globalcounter++);
    for(int i = 0; i < N; i++) {
      std::cout << "Test instance: " << i << std::endl;
      gSystem->Sleep(1000);
      AliSysInfo::AddStamp("test_manual_before", globalcounter++, i);
      TFile *reader = TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data());
      AliOADBContainer *cont2 = static_cast<AliOADBContainer *>(reader->Get(namecontainer.Data()));
      cont2->SetOwner(0);
      TList *l1=cont2->GetDefaultList();
      //for (Int_t i=0;i<l1->GetEntries();++i) printf("%i %p\n",i,l1->At(i));
      AliSysInfo::AddStamp("test_manual_open_cont", globalcounter++, i);
      delete reader;
      delete cont2;
      AliSysInfo::AddStamp("test_manual_after", globalcounter++, i);
      AliSysInfo::AddStamp("end", globalcounter++);
      std::cout << "End test Manual" << std::endl;
    }
  }
}

