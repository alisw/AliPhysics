/*
This a modification in a code sent by Marco Bregant, which, originally, created the OADB for misalignment matrices
In this macro, the histograms with BadChannels Factors are loaded and some TObjarrays are filled with these histograms.
At the end, a OADB container is created receiving these arrays.
The histograms with badchannels can be made by a list of towers ID ( see the macro CreateEMCAL_OADB_BadChannelsHistos.C ) or 
directly from the OCDB ( see AliEMCALOCDBTenderConverter.cxx from Jiri Kral )
*/
// ******* Create OADB Container for EMCal Bad Channels
void CreateEMCAL_OADB_BadChannels(const char *fileName10_1	=	"BadChannels2010_1.root",
				  const char *fileName10_2	=	"BadChannels2010_2.root",
				  const char *fileName10_3	=	"BadChannels2010_3.root",
				  const char *fileName11a	=	"BadChannels2011_11a.root", 
				  const char *fileName11b	=	"BadChannels2011_11b.root",
				  const char *fileName11c1	=	"BadChannels2011_11c1.root", 
				  const char *fileName11c2	=	"BadChannels2011_11c2.root", //Subperiod with different histograms
				  const char *fileName11c3	=	"BadChannels2011_11c3.root",
				  const char *fileName11d	=	"BadChannels2011_11d.root",
				  const char *fileName11e	=	"BadChannels2011_11e.root",
				  const char *fileName11f	=	"BadChannels2011_11f.root",
				  const char *fileName11h	=	"BadChannels2011_11h.root",

				  const char *fileName12a1	=	"BadChannels2012_12a1.root",
				  const char *fileName12a2	=	"BadChannels2012_12a2.root",

				  const char *fileName12b1	=	"BadChannels2012_12b1.root",
				  const char *fileName12b2	=	"BadChannels2012_12b2.root",
				  const char *fileName12b3	=	"BadChannels2012_12b3.root",
				  const char *fileName12b4	=	"BadChannels2012_12b4.root",
				  
				  const char *fileName12c1	=	"BadChannels2012_12c1.root",
				  const char *fileName12c2	=	"BadChannels2012_12c2.root",
				  const char *fileName12c3	=	"BadChannels2012_12c3.root",
				  const char *fileName12c4	=	"BadChannels2012_12c4.root",
				  
				  const char *fileName12d1	=	"BadChannels2012_12d1.root",
				  const char *fileName12d2	=	"BadChannels2012_12d2.root",
				  const char *fileName12d3	=	"BadChannels2012_12d3.root",
				  const char *fileName12d4	=	"BadChannels2012_12d4.root",
				  
				  const char *fileName12e	=	"BadChannels2012_12e.root",
				  const char *fileName12f	=	"BadChannels2012_12f.root",
				  const char *fileName12g	=	"BadChannels2012_12g.root",
				  const char *fileName12h	=	"BadChannels2012_12h.root",

		  
				)
{

    gSystem->Load("libOADB");      
    //gSystem->Load("libOADB"); //comment if compiled
      
    //LHC10
    TObjArray array10_1(12);		array10_1.SetName("BadChannels10_1");
    TObjArray array10_2(12);		array10_2.SetName("BadChannels10_2");
    TObjArray array10_3(12);		array10_3.SetName("BadChannels10_3");
   
    //LHC11a
    TObjArray array11a(12);		array11a.SetName("BadChannels11a");

    //LHC11b
    TObjArray array11b(12);		array11b.SetName("BadChannels11b");
    
    //LHC11c
    TObjArray array11c1(12);		array11c1.SetName("BadChannels11c1");
    TObjArray array11c2(12);		array11c2.SetName("BadChannels11c2");
    TObjArray array11c3(12);		array11c3.SetName("BadChannels11c3");

    //LHC11d
    TObjArray array11d(12);		array11d.SetName("BadChannels11d");

    //LHC11e
    TObjArray array11e(12);		array11e.SetName("BadChannels11e");

    //LHC11f
    TObjArray array11f(12);		array11f.SetName("BadChannels11f");

    //LHC11h
    TObjArray array11h(12);		array11h.SetName("BadChannels11h");

    //LHC12a
    TObjArray array12a1(12);		array12a1.SetName("BadChannels12a1");
    TObjArray array12a2(12);		array12a2.SetName("BadChannels12a2");

    //LHC12b
    TObjArray array12b1(12);		array12b1.SetName("BadChannels12b1");
    TObjArray array12b2(12);		array12b2.SetName("BadChannels12b2");
    TObjArray array12b3(12);		array12b3.SetName("BadChannels12b3");
    TObjArray array12b4(12);		array12b4.SetName("BadChannels12b4");
    
    //LHC12c
    TObjArray array12c1(12);		array12c1.SetName("BadChannels12c1");
    TObjArray array12c2(12);		array12c2.SetName("BadChannels12c2");
    TObjArray array12c3(12);		array12c3.SetName("BadChannels12c3");
    TObjArray array12c4(12);		array12c4.SetName("BadChannels12c4");
    
    //LHC12d
    TObjArray array12d1(12);		array12d1.SetName("BadChannels12d1");
    TObjArray array12d2(12);		array12d2.SetName("BadChannels12d2");
    TObjArray array12d3(12);		array12d3.SetName("BadChannels12d3");
    TObjArray array12d4(12);		array12d4.SetName("BadChannels12d4");
    
    //LHC12e
    TObjArray array12e(12);		array12e.SetName("BadChannels12e");
    
    //LHC12f
    TObjArray array12f(12);		array12f.SetName("BadChannels12f");

    //LHC12g
    TObjArray array12g(12);		array12g.SetName("BadChannels12g");

    //LHC12h
    TObjArray array12h(12);		array12h.SetName("BadChannels12h");
    
    TFile	*f10_1	=	new TFile(fileName10_1,"read"); 
    TFile	*f10_2	=	new TFile(fileName10_2,"read"); 
    TFile	*f10_3	=	new TFile(fileName10_3,"read"); 
    TFile	*f11a	=	new TFile(fileName11a,"read");
    TFile	*f11b	=	new TFile(fileName11b,"read");
    TFile	*f11c1	=	new TFile(fileName11c1,"read");
    TFile	*f11c2	=	new TFile(fileName11c2,"read");
    TFile	*f11c3	=	new TFile(fileName11c3,"read");
    TFile	*f11d	=	new TFile(fileName11d,"read");
    TFile	*f11e	=	new TFile(fileName11e,"read");
    TFile	*f11f	=	new TFile(fileName11f,"read");
    TFile	*f11h	=	new TFile(fileName11h,"read");
    TFile	*f12a1	=	new TFile(fileName12a1,"read");
    TFile	*f12a2	=	new TFile(fileName12a2,"read");
    TFile	*f12b1	=	new TFile(fileName12b1,"read");
    TFile	*f12b2	=	new TFile(fileName12b2,"read");
    TFile	*f12b3	=	new TFile(fileName12b3,"read");
    TFile	*f12b4	=	new TFile(fileName12b4,"read");
    TFile	*f12c1	=	new TFile(fileName12c1,"read");
    TFile	*f12c2	=	new TFile(fileName12c2,"read");
    TFile	*f12c3	=	new TFile(fileName12c3,"read");
    TFile	*f12c4	=	new TFile(fileName12c4,"read");
    TFile	*f12d1	=	new TFile(fileName12d1,"read");
    TFile	*f12d2	=	new TFile(fileName12d2,"read");
    TFile	*f12d3	=	new TFile(fileName12d3,"read");
    TFile	*f12d4	=	new TFile(fileName12d4,"read");
    TFile	*f12e	=	new TFile(fileName12e,"read");
    TFile	*f12f	=	new TFile(fileName12f,"read");
    TFile	*f12g	=	new TFile(fileName12g,"read");
    TFile	*f12h	=	new TFile(fileName12h,"read");
    
    //Create OADB container for BadChannels
    AliOADBContainer* con = new AliOADBContainer("AliEMCALBadChannels");

    char name[30];
    for (Int_t mod=0;mod<12;mod++){
	cout<<"SM "<< mod<<endl;
	
	if (mod<4)  {
	    sprintf(name,"EMCALBadChannelMap_Mod%d",mod);
	    cout<<"BadChannels2010:"<<name<<endl;
	    array10_1.Add(f10_1->Get(name));
	    array10_2.Add(f10_2->Get(name));
	    array10_3.Add(f10_3->Get(name));

	}
	    sprintf(name,"EMCALBadChannelMap_Mod%d",mod);
	    cout<<"BadChannels 2011 and 2012:"<<name<<endl;
	    array11a.Add(f11a->Get(name));
	    array11b.Add(f11b->Get(name));
	    array11c1.Add(f11c1->Get(name));
	    array11c2.Add(f11c2->Get(name));
	    array11c3.Add(f11c3->Get(name));
	    array11d.Add(f11d->Get(name));
	    array11e.Add(f11e->Get(name));
	    array11f.Add(f11f->Get(name));
	    array11h.Add(f11h->Get(name));
	    array12a1.Add(f12a1->Get(name));
	    array12a2.Add(f12a2->Get(name));    
	    array12b1.Add(f12b1->Get(name));
	    array12b2.Add(f12b2->Get(name));
	    array12b3.Add(f12b3->Get(name));
	    array12b4.Add(f12b4->Get(name));
	    array12c1.Add(f12c1->Get(name));
	    array12c2.Add(f12c2->Get(name));
	    array12c3.Add(f12c3->Get(name));
	    array12c4.Add(f12c4->Get(name));
	    array12d1.Add(f12d1->Get(name));
	    array12d2.Add(f12d2->Get(name));
	    array12d3.Add(f12d3->Get(name));
	    array12d4.Add(f12d4->Get(name));
	    array12e.Add(f12e->Get(name));
	    array12f.Add(f12f->Get(name));
	    array12g.Add(f12g->Get(name));
	    array12h.Add(f12h->Get(name));
   } //mod
    
    con->AddDefaultObject(&array10_1);
    con->AddDefaultObject(&array10_2);
    con->AddDefaultObject(&array10_3);
    
    con->AddDefaultObject(&array11a);
    con->AddDefaultObject(&array11b);
    con->AddDefaultObject(&array11c1);
    con->AddDefaultObject(&array11c2);
    con->AddDefaultObject(&array11c3);
    con->AddDefaultObject(&array11d);
    con->AddDefaultObject(&array11e);
    con->AddDefaultObject(&array11f);
    con->AddDefaultObject(&array11h);

    con->AddDefaultObject(&array12a1);
    con->AddDefaultObject(&array12a2);   
    con->AddDefaultObject(&array12b1);
    con->AddDefaultObject(&array12b2);
    con->AddDefaultObject(&array12b3);
    con->AddDefaultObject(&array12b4);
    con->AddDefaultObject(&array12c1);
    con->AddDefaultObject(&array12c2);
    con->AddDefaultObject(&array12c3);
    con->AddDefaultObject(&array12c4);
    con->AddDefaultObject(&array12d1);
    con->AddDefaultObject(&array12d2);
    con->AddDefaultObject(&array12d3);
    con->AddDefaultObject(&array12d4);     
    con->AddDefaultObject(&array12e);   
    con->AddDefaultObject(&array12f);   
    con->AddDefaultObject(&array12g);   
    con->AddDefaultObject(&array12h);   
    
    //Establishing run number with the correct objects
    con->AppendObject(&array10_1,112000,120742);
    con->AppendObject(&array10_2,120743,121984);
    con->AppendObject(&array10_1,121985,124186);
    con->AppendObject(&array10_3,124187,125296);
    con->AppendObject(&array10_1,125297,140000);

    con->AppendObject(&array11a,144871,146860);
    con->AppendObject(&array11b,148531,150629);
    con->AppendObject(&array11c1,151636,153569);
    con->AppendObject(&array11c2,153570,154733);
    con->AppendObject(&array11c3,154734,155384);
    con->AppendObject(&array11d,156477,159635);
    con->AppendObject(&array11e,160670,162740);
    con->AppendObject(&array11f,162933,165746);
    con->AppendObject(&array11h,166529,170673);

    con->AppendObject(&array12a1,172439,176325);
    con->AppendObject(&array12a2,176326,177295);
    
    con->AppendObject(&array12b1,177381,177383);
    con->AppendObject(&array12b2,177384,177443);
    con->AppendObject(&array12b3,177444,177682);
    con->AppendObject(&array12b2,177683,177843); 
    con->AppendObject(&array12b4,177844,177849);
    con->AppendObject(&array12b2,177850,178220);
  
    con->AppendObject(&array12c1,179569,180126);
    con->AppendObject(&array12c2,180127,180194);
    con->AppendObject(&array12c3,180195,180200);
    con->AppendObject(&array12c4,180289,182740); 
    con->AppendObject(&array12c1,182741,182744);
   
    con->AppendObject(&array12d1,183913,184481);
    con->AppendObject(&array12d2,184482,185455);
    con->AppendObject(&array12d3,185456,185784);
    con->AppendObject(&array12d4,185909,186035); 
    con->AppendObject(&array12d2,186036,186320); 

    con->AppendObject(&array12e,186365,186602); 
    con->AppendObject(&array12f,186668,188123); 
    con->AppendObject(&array12g,188356,188503); 
    con->AppendObject(&array12h,189122,192732); 

    con->WriteToFile("BetaBadChannels.root");   
}

