#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.H>
#include "AliITSGainSSDv2.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#endif
Int_t GainCalibration(AliCDBManager * man, Float_t* gainP, Float_t* gainN) 
{
	const Int_t fgkSSDMODULES = 1698;
  	static const Int_t fgkDefaultNStripsSSD = 768;
  	AliITSGainSSDv2 *gainSSD = new AliITSGainSSDv2();
  	AliCDBEntry *entryGainSSD = man->Get("ITS/Calib/GainSSD");
  	TObject *empty = (TObject *)entryGainSSD->GetObject();
  	TString objectname = empty->GetName();
  	if(objectname=="AliITSGainSSDv2") 
	{
    		cout<<"Reading the new format of the calibration file..."<<endl;
    		gainSSD = (AliITSGainSSDv2 *)entryGainSSD->GetObject();
  	}
	else
	{
		cout<<"Now OCDB object"<<endl;
		return 0 ;
	}
  	for (Int_t i = 0; i < fgkSSDMODULES; i++) 
	{
  			//for (Int_t i = 0; i < 1; i++) {
    			//AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
   		for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) 
	 	{
     	 		gainP[i*fgkDefaultNStripsSSD+j] = gainSSD->GetGainP(i,j);
      			gainN[i*fgkDefaultNStripsSSD+j]= gainSSD->GetGainN(i,j);
     	
    		}//strip loop
  	}//module loop
	return 1;
}


void UpdateSSDGainOCDB(TString gainfile,Float_t globalMPV=84.0,const char* type = "alien", Int_t runNumber = 0, Int_t runRangemin=0 ,Int_t runRangemax= AliCDBRunRange::Infinity()) 
{
	if((runRangemax<runRangemin))
	{
		cout<<"wrong run range "<<endl;
		return;
	}	
  	const Int_t fgkPNStripsPerModule = 768;
  	const Int_t fgkNumberOfSSDModules = 1698;
	Float_t gainPold[fgkPNStripsPerModule* fgkNumberOfSSDModules];
	Float_t gainNold[fgkPNStripsPerModule* fgkNumberOfSSDModules];
	
	for(int i=0;i<fgkPNStripsPerModule* fgkNumberOfSSDModules;i++)
	{
			gainPold[i]=0.0;
			gainNold[i]=0.0;
	}
	
	TString gType = type;
  	AliCDBManager * man = AliCDBManager::Instance();
	if(gType == "alien") 
   		 man->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB/");
  	else if(gType == "local") 
    		man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  	else 
	{
    		cout<<"Allowed types: local or alien!!!"<<endl;
    		return;
  	}
	  man->SetRun(runNumber);
	if(!GainCalibration(man,gainPold,gainNold))
		return;
		
	char buffer[256];
	ifstream ingainfile(gainfile);
	if(!ingainfile.good())
		return;
	
	
	 ingainfile.getline(buffer,256);
	cout<<buffer<<endl;
	
	
	Int_t module=0;
	Int_t gainflag=0;
	Float_t gainP=0;
	Float_t gainN=0;
	Float_t fMPV=0;
	Float_t corrP[1698];
	Float_t corrN[1698];
	
	for(int jj=0;jj<1698;jj++)
    	{
		ingainfile>>module;
		ingainfile>>gainflag;
		ingainfile>>gainP;
		ingainfile>>gainN;
		ingainfile>>fMPV;
		if(gainflag==1)
		{
			corrP[module]=gainP*globalMPV/fMPV;
			corrN[module]=gainN*globalMPV/fMPV;
			
		}
		else
		{
			corrP[module]=1.0;
			corrN[module]=1.0;
		}
			
	}
	
	
  	TObjArray *array = new TObjArray();
  	array->SetName("Gain");
  	AliITSGainSSDv2  *mc;
	mc = new AliITSGainSSDv2();
  	for(Int_t mind = 0; mind < fgkNumberOfSSDModules; mind++) 
	{
    		for (Int_t i = 0; i < fgkPNStripsPerModule; i++) 
		{
			mc->AddGainP(mind,i, gainPold[mind* fgkPNStripsPerModule+i]*corrP[mind]);
     			mc->AddGainN(mind,i,gainNold[mind* fgkPNStripsPerModule+i]*corrN[mind]);			
    		}
  	}

	man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	man->SetRun(0);
	//AliCDBId id("ITS/Calib/GainSSD",0,AliCDBRunRange::Infinity());
	AliCDBId id("ITS/Calib/GainSSD",runRangemin,runRangemax);
	AliCDBMetaData *metadata= new AliCDBMetaData();
	metadata->SetResponsible("Marek.Chojnacki@cern.ch");
	metadata->SetComment("Default values for the SSD gain calibration");
	//man->Put(array,id,metadata);
	man->Put(mc,id,metadata);
	 TFile* fout=TFile::Open("AliITSGainSSD_new.root","recreate");
	mc ->Write();	
	 	fout->Close();
	
}  
