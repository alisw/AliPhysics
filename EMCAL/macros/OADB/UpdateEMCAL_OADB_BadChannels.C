/*
This a modification in a code sent by Marco Bregant, which, originally, created the OADB for misalignment matrices
In this macro, the histograms with BadChannels Factors are loaded and some TObjarrays are filled with these histograms.
At the end, a OADB container is created receiving these arrays.
The histograms with badchannels can be made by a list of towers ID ( see the macro CreateEMCAL_OADB_BadChannelsHistos.C ) or 
directly from the OCDB ( see AliEMCALOCDBTenderConverter.cxx from Jiri Kral )
*/
// ******* Create OADB Container for EMCal Bad Channels

TObjArray GetHistoObject( Int_t runNum, Bool_t LocalOCDB=0){
	Int_t i;
	char buf[100];

	TH2D *histo;

	TFile *outFile;
	
	AliCDBManager *man;
	AliCDBStorage *stor;
	AliCaloCalibPedestal *ped;

	// created the OCDB manager
	man = AliCDBManager::Instance();
	if(LocalOCDB) stor = man->GetStorage("local://$HOME/OADB/OCDB/"); //if you download the OCDB files locally.
	else {  
	  man->SetDefaultStorage("raw://");
	  man->SetRun(runNum);
	  stor = man->GetDefaultStorage();
	}
	ped = (AliCaloCalibPedestal*)(stor->Get("EMCAL/Calib/Pedestals", runNum)->GetObject());

	// get the array of histos
	TObjArray map = ped->GetDeadMap();
	for(int i = 0; i < map.GetEntries(); i++ ){
		histo = (TH2D*)(map[i]);
		printf("\n !!! EMCALBadChannelMap_Mod%d",i );
		histo->SetName( Form("EMCALBadChannelMap_Mod%d", i ));
		histo->SetTitle( Form("EMCALBadChannelMap_Mod%d", i ));
	}

	return map;
}

void UpdateEMCAL_OADB_BadChannels(const char *fileNameOADB="$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root",Bool_t LocalOCDB=1)
{

 ///   gSystem->Load("libOADB");      
    //gSystem->Load("libOADB"); //comment if compiled

    //Create OADB container for BadChannels
    AliOADBContainer *con=new AliOADBContainer("");
    con->InitFromFile(Form(fileNameOADB),"AliEMCALBadChannels"); 

//AliOADBContainer* con = new AliOADBContainer("AliEMCALBadChannels");

    //! List of brand new arrays
    //LHC13g
    TObjArray *array = new TObjArray(12);		
    array->SetName("BadChannels");
   // 'AliEMCALOCDBTenderConverter.cxx(176326,"BadChannels2012_12a1.root")
    char name[30];
    
    TObjArray array12a1 = GetHistoObject(176326,LocalOCDB);// Run176326_176700_v5_s0.root
    TObjArray array12a2 = GetHistoObject(176701,LocalOCDB);// Run176701_176714_v5_s0.root
    TObjArray array12a3 = GetHistoObject(176715,LocalOCDB);// Run176715_176729_v5_s0.root
    TObjArray array12a4 = GetHistoObject(176730,LocalOCDB);// Run176730_176858_v5_s0.root
    TObjArray array12a5 = GetHistoObject(176859,LocalOCDB);// Run176859_177295_v5_s0.root
    
    TObjArray array12b1 = GetHistoObject(177320,LocalOCDB);// Run177320_999999999_v2_s0 --> 177381-177383
    TObjArray array12b2 = GetHistoObject(177384,LocalOCDB);// Run177384_178220_v5_s0.root
    TObjArray array12b3 = GetHistoObject(177444,LocalOCDB);// Run177444_177682_v7_s0.root
    TObjArray array12b4 = GetHistoObject(177844,LocalOCDB);// Run177844_177849_v7_s0.root     
    TObjArray array12b5 = GetHistoObject(178220,LocalOCDB);// Run177384_178220_v9_s0.root

    TObjArray array12c1  = GetHistoObject(179569,LocalOCDB);// 
    TObjArray array12c2  = GetHistoObject(180000,LocalOCDB);// 
    TObjArray array12c3  = GetHistoObject(180127,LocalOCDB);// 
    TObjArray array12c4  = GetHistoObject(180134,LocalOCDB);// 
    TObjArray array12c5  = GetHistoObject(180195,LocalOCDB);// 
    TObjArray array12c6  = GetHistoObject(180201,LocalOCDB);// 
    TObjArray array12c7  = GetHistoObject(180289,LocalOCDB);// 
    TObjArray array12c8  = GetHistoObject(181617,LocalOCDB);// 
    TObjArray array12c9  = GetHistoObject(182106,LocalOCDB);// 
    TObjArray array12c10 = GetHistoObject(182624,LocalOCDB);// 
    TObjArray array12c11 = GetHistoObject(182730,LocalOCDB);// 

    
    TObjArray array12d1  = GetHistoObject(183913,LocalOCDB);// 
    TObjArray array12d2  = GetHistoObject(184482,LocalOCDB);// 
    TObjArray array12d3  = GetHistoObject(185456,LocalOCDB);// 
    TObjArray array12d4  = GetHistoObject(185909,LocalOCDB);// 
    TObjArray array12d5  = GetHistoObject(186036,LocalOCDB);// 

    TObjArray array12e  = GetHistoObject(186365,LocalOCDB);// 

    TObjArray array12f1  = GetHistoObject(186668,LocalOCDB);// 
    TObjArray array12f2  = GetHistoObject(187534,LocalOCDB);// 
    TObjArray array12f3  = GetHistoObject(188122,LocalOCDB);// 

    TObjArray array12g  = GetHistoObject(188356,LocalOCDB);// 

    TObjArray array12h1   = GetHistoObject(189122,LocalOCDB);// 
    TObjArray array12h2   = GetHistoObject(189301,LocalOCDB);//
    TObjArray array12h3   = GetHistoObject(189518,LocalOCDB);//
    TObjArray array12h4   = GetHistoObject(189699,LocalOCDB);//
    TObjArray array12h5   = GetHistoObject(190335,LocalOCDB);//
    TObjArray array12h6   = GetHistoObject(190393,LocalOCDB);//
    
    TObjArray array12h7   = GetHistoObject(190449,LocalOCDB);//
    TObjArray array12h8   = GetHistoObject(190895,LocalOCDB);//
    TObjArray array12h9   = GetHistoObject(191022,LocalOCDB);//
    TObjArray array12h10  = GetHistoObject(191129,LocalOCDB);//
    TObjArray array12h11  = GetHistoObject(191163,LocalOCDB);//
    TObjArray array12h12  = GetHistoObject(191449,LocalOCDB);//
    TObjArray array12h13  = GetHistoObject(191485,LocalOCDB);//

    TObjArray array12i1   = GetHistoObject(192745,LocalOCDB);//
    TObjArray array12i1b  = GetHistoObject(193000,LocalOCDB);//   
    TObjArray array12i2   = GetHistoObject(193301,LocalOCDB);//
    TObjArray array12i3   = GetHistoObject(193750,LocalOCDB);//


        
    con->AddDefaultObject(&array12a1);
    con->AddDefaultObject(&array12a2);
    con->AddDefaultObject(&array12a3);
    con->AddDefaultObject(&array12a4);
    con->AddDefaultObject(&array12a5);
    
    con->AddDefaultObject(&array12b1);
    con->AddDefaultObject(&array12b2);
    con->AddDefaultObject(&array12b3);
    con->AddDefaultObject(&array12b4);
    con->AddDefaultObject(&array12b5);

    con->AddDefaultObject(&array12c1);
    con->AddDefaultObject(&array12c2);
    con->AddDefaultObject(&array12c3);
    con->AddDefaultObject(&array12c4);
    con->AddDefaultObject(&array12c5);
    con->AddDefaultObject(&array12c6);
    con->AddDefaultObject(&array12c7);
    con->AddDefaultObject(&array12c8);
    con->AddDefaultObject(&array12c9);
    con->AddDefaultObject(&array12c10);
    con->AddDefaultObject(&array12c11);
    
    con->AddDefaultObject(&array12d1);
    con->AddDefaultObject(&array12d2);
    con->AddDefaultObject(&array12d3);
    con->AddDefaultObject(&array12d4);
    con->AddDefaultObject(&array12d5);

    con->AddDefaultObject(&array12e);
    
    con->AddDefaultObject(&array12d1);
    con->AddDefaultObject(&array12d2);
    con->AddDefaultObject(&array12d3);
    con->AddDefaultObject(&array12d4);
    con->AddDefaultObject(&array12d5);

    con->AddDefaultObject(&array12f1);
    con->AddDefaultObject(&array12f2);
    con->AddDefaultObject(&array12f3);

    con->AddDefaultObject(&array12g); 
    
    con->AddDefaultObject(&array12h1);
    con->AddDefaultObject(&array12h2);
    con->AddDefaultObject(&array12h3);
    con->AddDefaultObject(&array12h4);
    con->AddDefaultObject(&array12h5);    
    con->AddDefaultObject(&array12h6);    
    con->AddDefaultObject(&array12h7);    
    con->AddDefaultObject(&array12h8);    
    con->AddDefaultObject(&array12h9);    
    con->AddDefaultObject(&array12h10);    
    con->AddDefaultObject(&array12h11);    
    con->AddDefaultObject(&array12h12);    
    con->AddDefaultObject(&array12h13);    

    con->AddDefaultObject(&array12i1);
    con->AddDefaultObject(&array12i1b);
    con->AddDefaultObject(&array12i2);
    con->AddDefaultObject(&array12i3);
  
    con->UpdateObject(con->GetIndexForRun(176326),&array12a1,176326,176700); // Run176326_176700_v5_s0.root 
    con->AppendObject(&array12a2,176701,176714); //	Run176701_176714_v5_s0.root
    con->AppendObject(&array12a3,176715,176729); // Run176715_176729_v4_s0.root
    con->AppendObject(&array12a4,176730,176858); // Run176730_176858_v5_s0.root
    con->AppendObject(&array12a5,176859,177295); // Run176869_177295_v5_s0.root
 
    
   con->UpdateObject(con->GetIndexForRun(177383),&array12b1,177320,177383); // Run177320_999999999_v2_s0.root 
   con->UpdateObject(con->GetIndexForRun(177384),&array12b2,177384,177443); // Run177384_178220_v9_s0.root
   con->UpdateObject(con->GetIndexForRun(177444),&array12b3,177444,177682); // Run177444_177682_v10_s0.root
   con->UpdateObject(con->GetIndexForRun(177844),&array12b4,177844,177849); // Run177844_177849_v10_s0.root
   con->UpdateObject(con->GetIndexForRun(177850),&array12b5,177850,178220); // Run177384_178220_v9_s0.root

   con->UpdateObject(con->GetIndexForRun(179569),&array12c1 ,179569,179999); // Run176859_177295_v6_s0.root
   con->AppendObject(&array12c2 ,180000,180126); // Run180000_180126_v4_s0.root
   con->UpdateObject(con->GetIndexForRun(180127),&array12c3 ,180127,180133); // Run180127_180133_v5_s0.root
   con->AppendObject(&array12c4 ,180134,180194); // Run180134_180194_v5_s0.root
   con->UpdateObject(con->GetIndexForRun(180195),&array12c5 ,180195,180200); // Run180195_180200_v5_s0.root
   con->AppendObject(&array12c6 ,180201,180288); // Run180201_180288_v4_s0.root
   con->UpdateObject(con->GetIndexForRun(180289),&array12c7 ,180289,181616); // Run180289_181616_v5_s0.root
   con->AppendObject(&array12c8 ,181617,182105); // Run181617_182105_v5_s0.root
   con->AppendObject(&array12c9 ,182106,182623); // Run182106_182623_v5_s0.root
   con->AppendObject(&array12c10,182624,182729); // Run182624_182729_v5_s0.root
   con->UpdateObject(con->GetIndexForRun(182744),&array12c11,182730,182744); // Run182730_182744_v5_s0.root
   
   
   con->UpdateObject(con->GetIndexForRun(183913),&array12d1,183913,184481); // Run183913_184481_v11_s0.root
   con->UpdateObject(con->GetIndexForRun(184482),&array12d2,184482,185455); // def Run183913_186320_v10_s0.root 
   con->UpdateObject(con->GetIndexForRun(185456),&array12d3,185456,185784); // Run185456_185784_v11_s0.root 
   con->UpdateObject(con->GetIndexForRun(185909),&array12d4,185909,186035); // Run185909_186035_v11_s0.root
 
   con->UpdateObject(con->GetIndexForRun(186036),&array12d5,186036,186320); // def Run183913_186320_v10_s0.root  
 
   con->UpdateObject(con->GetIndexForRun(186365),&array12e,186365,186602); // Run186365_186602_v4_s0.root 
 
   con->UpdateObject(con->GetIndexForRun(186668),&array12f1,186668,187533); // def Run186668_188123_v4_s0.root
   con->AppendObject(&array12f2,187534,187562); // Run187534_187562_v5_s0.root
   con->AppendObject(&array12f3,187563,188123); // def Run186668_188123_v4_s0.root
  
   con->UpdateObject(con->GetIndexForRun(188356),&array12g,188356,188503);  // Run188356_188503_v5_s0.root
   
   con->UpdateObject(con->GetIndexForRun(189122),&array12h1,189122,189246); // Run189122_190392_v4_s0.root
   con->AppendObject(&array12h2,189301,189494); // Run189301_189494_v5_s0.root
   con->AppendObject(&array12h3,189518,189698); // Run189518_189698_v5_s0.root
   con->AppendObject(&array12h4,189699,190334); // Run189699_190334_v5_s0.root
   con->AppendObject(&array12h5,190335,190392); // Run190335_190392_v5_s0.root
   con->AppendObject(&array12h6,190393,190425); // Run190393_190425_v5_s0.root
   con->AppendObject(&array12h7,190449,190894); // def Run190393_192732_v4_s0.root
   con->AppendObject(&array12h8,190895,190984); // Run190895_190984_v5_s0.root
   con->AppendObject(&array12h9,191022,191121); // def  Run190393_192732_v4_s0.root
   con->AppendObject(&array12h10,191129,191159); // Run191129_191159_v5_s0.root
   con->AppendObject(&array12h11,191163,191448); // def  Run190393_192732_v4_s0.root
   con->AppendObject(&array12h12,191449,191451); // Run191449_191451_v5_s0.root
   con->AppendObject(&array12h13,191485,192732); // def  Run190393_192732_v4_s0.root
  
 
   con->UpdateObject(con->GetIndexForRun(192745),&array12i1,192745,192824); //Run192745_193300_v4_s0.root
   con->UpdateObject(con->GetIndexForRun(193300),&array12i1b,192825,193300); //Run192745_193300_v4_s0.root
   con->UpdateObject(con->GetIndexForRun(193004),&array12i2,193301,193749); //Run193301_193749_v4_s0.root
   con->AppendObject(&array12i3,193750,193766); //Run193750_193766_v4_s0.root

   con->WriteToFile("BetaBadChannels.root");   
}

    /*
179569-179999
180000-180126
180127-180133
180134-180194
180195-180200
180201-180288
180289-181616
181617-182105
182106-182623
182624-182729
182730-182744
*/
  /*
    Default: 183913-186320
run 183913-184481
run185456-185784
run185909-186035
186036-186320(def)
*/

/*
 Default: 183913-186320
run 183913-184481
184482-185455(def)
run185456-185784
run185909-186035
186036-186320(def)   
*/
/*

1 file runs186365 to 186602
LHC12f:
2 files:
186668 to 188123
187534-187562
LHC12g
1 file: run 188356 to 188503
*/


/*
LHC12h:
189122-189392
189393-192732 (default)

189301-189474
189518-189698
189699-190334
190335-190425
190449-190894 (def)
190895-190984
191022-191121 (def)
191129-191159
191163-191448 (def)
191449-191451
191485-192732 (def)


LHC12i
192745-193300
193301-193749
193750-193766
*/
