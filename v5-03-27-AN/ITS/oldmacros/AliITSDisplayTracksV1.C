Int_t AliITSDisplayTracksV1(const char *fname = "galice.root", Int_t evNum = 0) {

	/*********************************************************************
	 *                                                                   *
	 *  Macro to plot reconstructed tracks with the Kalman Filter V1 on  *
	 *  top of the ITS detailed geometry.                                *
	 *                                                                   *
	 *  Authors: Angela Badala' and Roberto Barbera                      *
	 *                                                                   *
	 *********************************************************************/


	// First of all, here are put some variable declarations
	// that are useful in the following part:
	Int_t nparticles; // number of particles
	// ITS module coordinates [layer = 1, ladder = 2, det = 3] and absolute ID[0] of module [0]
	TArrayI ID(4);
	Int_t nmodules, dtype; // Total number of modules and module type (SSD, SPD, SDD)
	Float_t *x = 0, *y = 0, *z = 0; // Arrays where to store read coords
	Bool_t *St = 0; // Status of the track (hits only)

        // create the Canvas

        c1 = new TCanvas("c1","Track display",50,50,800,800);
	c1 -> Divide(2,2,0.001,0.001);

	// It's necessary to load the gAlice shared libs
	// if they aren't already stored in memory...
	if (gClassTable->GetID("AliRun") < 0) {
   	gROOT->LoadMacro("loadlibs.C");
		loadlibs();
	}
  // Anyway, this macro needs to read a gAlice file, so it
  // clears the gAlice object if there is already one in memory...
  else {
		if(gAlice){
			delete gAlice;
			gAlice = 0;
		}
	}

	// Now is opened the Root input file containing Geometry, Kine and Hits
  // by default its name must be "galice.root".
  // When the file is opened, its contens are shown.
	TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
	if (!file) file = new TFile(fname);
		file->ls();
	
	// Then, the macro gets the AliRun object from file.
	// If this object is not present, an error occurs
	// and the execution is stopped.
	// Anyway, this operation needs some time,
	// don't worry about an apparent absence of output and actions...
	cout << "\nSearching in '" << fname << "' for an AliRun object ... " << flush;
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice)
		cout << "FOUND!" << endl;
	else {
		cout<<"NOT FOUND! The Macro can't continue!" << endl;
		return 0;
	}
	
	// Then, the macro selects the event number specified. Default is 0.
	nparticles = gAlice->GetEvent(evNum);
	cout << "\nNumber of particles   = " << nparticles << endl;
	if (!nparticles) {
		cout << "With no particles I can't do much... Goodbye!" << endl;
		return 0;
	}
	
	// The next step is filling the ITS from the AliRun object.
	AliITS *ITS = (AliITS*)gAlice->GetModule("ITS");
  ITS->InitModules(-1, nmodules);
  cout << "Number of ITS modules = " << nmodules << endl;
	cout << "\nFilling modules (it takes a while, now)..." << flush;
	ITS->FillModules(0, 0, nmodules, " ", " ");
  cout << "DONE!" << endl;


	// Gets geometry...
	AliITSgeom *gm  = ITS->GetITSgeom();
	if (!gAlice) {
		cout << "A problem occurred when getting the geometry!!!" << endl;
		return;
	}

	
      
//inserire qui la geometria

    TNode *node, *top;
    
    const Int_t kColorITS=kRed;
    //
    //top = gAlice->GetGeometry()->GetNode("alice");


  // DETAILED GEOMETRY

  TNode *sub1node, *sub2node, *sub3node, *sub4node, *sub5node;

  // Define some variables for SPD

  Float_t dits[100];
  Float_t dits1[3], di101[3], di107[3], di10b[3];  // for layer 1 
  Float_t di103[3], di10a[3];                      // for layer 1
  Float_t dits2[3], di1d1[3], di1d7[3], di20b[3];  // for layer 2
  Float_t di1d3[3], di20a[3];                      // for layer 2  

  Float_t ddet1=200.;     // total detector thickness on layer 1 (micron)
  Float_t dchip1=200.;    // total chip thickness on layer 1 (micron)
  
  Float_t ddet2=200.;     // total detector thickness on layer 2 (micron)                         
  Float_t dchip2=200.;    // total chip thickness on layer 2 (micron)
  
  Float_t dbus=300.;      // total bus thickness on both layers (micron)
  
  ddet1  = ddet1*0.0001/2.; // conversion from tot length in um to half in cm
  ddet2  = ddet2*0.0001/2.; // conversion from tot length in um to half in cm	
  dchip1 = dchip1*0.0001/2.;// conversion from tot length in um to half in cm	
  dchip2 = dchip2*0.0001/2.;// conversion from tot length in um to half in cm	
  dbus   = dbus*0.0001/2.;  // conversion from tot length in um to half in cm       
		
  Float_t deltax, deltay; 

  Int_t thickness = 1;
  Int_t option    = 2;


  // Define some variables for SDD
  // SDD detector ladder

  Float_t ySDD;
  Float_t I302dits[3], I402dits[3], I004dits[3], I005dits[3];
  Float_t Y_SDD_sep = 0.20;
  Float_t Z_SDD_lay3[6] = {18.55, 10.95, 3.70, -3.70, -11.20, -18.35};
  Float_t Z_SDD_lay4[8] = {25.75, 18.60, 11.00, 3.70, -3.70, -11.20, -18.45, -26.05};

  // Rotation matrices
    
  // SPD - option 'a' 
  
  if (option == 1) {  
  
     new TRotMatrix("itsrot238","itsrot238",90.0,144.0,90.0,234.0,0.0,0.0);
     new TRotMatrix("itsrot236","itsrot236",90.0,180.013702,90.0,270.013702,0.0,0.0);
     new TRotMatrix("itsrot239","itsrot239",90.0,216.0,90.0,306.0,0.0,0.0);     
     new TRotMatrix("itsrot233","itsrot233",90.0,252.000504,90.0,342.000488,0.0,0.0 );     
     new TRotMatrix("itsrot240","itsrot240",90.0,288.0,90.0,18.0,0.0,0.0);
     new TRotMatrix("itsrot241","itsrot241",90.0,324.0,90.0,54.0,0.0,0.0);          
     new TRotMatrix("itsrot242","itsrot242",90.0,36.0,90.0,126.0,0.0,0.0); 
     new TRotMatrix("itsrot234","itsrot234",90.0,71.9991,90.0,161.9991,0.0,0.0);     
     new TRotMatrix("itsrot243","itsrot243",90.0,108.0,90.0,198.0,0.0,0.0);  
     new TRotMatrix("itsrot244","itsrot244",90.0,180.0,90.0,270.0,0.0,0.0);
     new TRotMatrix("itsrot245","itsrot245",90.0,162.0,90.0,252.0,0.0,0.0);
     new TRotMatrix("itsrot246","itsrot246",90.0,310.0,90.0,40.0,0.0,0.0);
     new TRotMatrix("itsrot247","itsrot247",90.0,319.0,90.0,49.0,0.0,0.0);
     new TRotMatrix("itsrot248","itsrot248",90.0,328.0,90.0,58.0,0.0,0.0);
     new TRotMatrix("itsrot249","itsrot249",90.0,337.0,90.0,67.0,0.0,0.0);     
                 
  }   

  // SPD - option 'b' (this is the default)  

  if (option == 2) {  
  
     new TRotMatrix("itsrot233","itsrot233",90.0,252.000504,90.0,342.000488,0.0,0.0);
     new TRotMatrix("itsrot244","itsrot244",90.0,216.0,90.0,306.0,0.0,0.0);
     new TRotMatrix("itsrot236","itsrot236",90.0,180.013702,90.0,270.013702,0.0,0.0);  
     new TRotMatrix("itsrot245","itsrot245",90.0,36.0,90.0,126.0,0.0,0.0);     
     new TRotMatrix("itsrot234","itsrot234",90.0,71.9991,90.0,161.9991,0.0,0.0);  
     new TRotMatrix("itsrot246","itsrot246",90.0,108.0,90.0,198.0,0.0,0.0);    
     new TRotMatrix("itsrot247","itsrot247",90.0,144.0,90.0,234.0,0.0,0.0);
     new TRotMatrix("itsrot248","itsrot248",90.0,288.0,90.0,18.0,0.0,0.0);     
     new TRotMatrix("itsrot249","itsrot249",90.0,324.0,90.0,54.0,0.0,0.0);       
     new TRotMatrix("itsrot238","itsrot238",90.0,180.0,90.0,270.0,0.0,0.0);
     new TRotMatrix("itsrot239","itsrot239",90.0,162.0,90.0,252.0,0.0,0.0);     
     new TRotMatrix("itsrot240","itsrot240",90.0,310.0,90.0,40.0,0.0,0.0);
     new TRotMatrix("itsrot241","itsrot241",90.0,319.0,90.0,49.0,0.0,0.0);
     new TRotMatrix("itsrot242","itsrot242",90.0,328.0,90.0,58.0,0.0,0.0);
     new TRotMatrix("itsrot243","itsrot243",90.0,337.0,90.0,67.0,0.0,0.0);

  }   
     
  // SDD
  
  new TRotMatrix("itsrot321","itsrot321",90.0,12.86,90.0,102.86,0.0,0.0);	 
  new TRotMatrix("itsrot333","itsrot333",90.0,38.57,90.0,128.57,0.0,0.0);
  new TRotMatrix("itsrot336","itsrot336",90.0,64.29,90.0,154.29,0.0,0.0);	
  new TRotMatrix("itsrot350","itsrot350",90.0,90.0,90.0,180.0,0.0,0.0);    
  new TRotMatrix("itsrot313","itsrot313",90.0,115.71,90.0,205.71,0.0,0.0);   
  new TRotMatrix("itsrot311","itsrot311",90.0,141.43,90.0,231.43,0.0,0.0);
  new TRotMatrix("itsrot310","itsrot310",90.0,167.14,90.0,257.14,0.0,0.0);  
  new TRotMatrix("itsrot386","itsrot386",90.0,192.86,90.0,282.86,0.0,0.0);    
  new TRotMatrix("itsrot309","itsrot309",90.0,218.57,90.0,308.57,0.0,0.0);  
  new TRotMatrix("itsrot308","itsrot308",90.0,244.29,90.0,334.29,0.0,0.0);  
  new TRotMatrix("itsrot356","itsrot356",90.0,270.0,90.0,0.0,0.0,0.0);   
  new TRotMatrix("itsrot307","itsrot307",90.0,295.71,90.0,25.71,0.0,0.0);  
  new TRotMatrix("itsrot306","itsrot306",90.0,321.43,90.0,51.43,0.0,0.0); 
  new TRotMatrix("itsrot305","itsrot305",90.0,347.14,90.0,77.14,0.0,0.0);		
  new TRotMatrix("itsrot335","itsrot335",90.0,8.18,90.0,98.18,0.0,0.0); 
  new TRotMatrix("itsrot332","itsrot332",90.0,24.55,90.0,114.55,0.0,0.0);  
  new TRotMatrix("itsrot331","itsrot331",90.0,40.91,90.0,130.91,0.0,0.0);	 
  new TRotMatrix("itsrot366","itsrot366",90.0,57.27,90.0,147.27,0.0,0.0);	
  new TRotMatrix("itsrot330","itsrot330",90.0,73.64,90.0,163.64,0.0,0.0);	   
  new TRotMatrix("itsrot350","itsrot350",90.0,90.0,90.0,180.0,0.0,0.0);    
  new TRotMatrix("itsrot329","itsrot329",90.0,106.36,90.0,196.36,0.0,0.0);  
  new TRotMatrix("itsrot328","itsrot328",90.0,122.73,90.0,212.73,0.0,0.0);  
  new TRotMatrix("itsrot327","itsrot327",90.0,139.09,90.0,229.09,0.0,0.0);  
  new TRotMatrix("itsrot326","itsrot326",90.0,155.45,90.0,245.45,0.0,0.0); 
  new TRotMatrix("itsrot325","itsrot325",90.0,171.82,90.0,261.82,0.0,0.0);  
  new TRotMatrix("itsrot324","itsrot324",90.0,188.18,90.0,278.18,0.0,0.0);   
  new TRotMatrix("itsrot323","itsrot323",90.0,204.55,90.0,294.55,0.0,0.0);   
  new TRotMatrix("itsrot322","itsrot322",90.0,220.91,90.0,310.91,0.0,0.0);  
  new TRotMatrix("itsrot320","itsrot320",90.0,237.27,90.0,327.27,0.0,0.0);  
  new TRotMatrix("itsrot319","itsrot319",90.0,253.64,90.0,343.64,0.0,0.0);  
  new TRotMatrix("itsrot318","itsrot318",90.0,270.0,90.0,360.0,0.0,0.0);  
  new TRotMatrix("itsrot317","itsrot317",90.0,286.36,90.0,16.36,0.0,0.0);  
  new TRotMatrix("itsrot316","itsrot316",90.0,302.73,90.0,32.73,0.0,0.0);	
  new TRotMatrix("itsrot315","itsrot315",90.0,319.09,90.0,49.09,0.0,0.0);	
  new TRotMatrix("itsrot314","itsrot314",90.0,335.45,90.0,65.45,0.0,0.0); 
  new TRotMatrix("itsrot334","itsrot334",90.0,351.82,90.0,81.82,0.0,0.0);	 
      
  //SSD 
  
  new TRotMatrix("itsrot504","itsrot504",90.0,127.06,90.0,217.06,0.0,0.0);  
  new TRotMatrix("itsrot505","itsrot505",90.0,116.47,90.0,206.47,0.0,0.0);  
  new TRotMatrix("itsrot506","itsrot506",90.0,105.88,90.0,195.88,0.0,0.0);  
  new TRotMatrix("itsrot507","itsrot507",90.0,95.29,90.0,185.29,0.0,0.0);  
  new TRotMatrix("itsrot508","itsrot508",90.0,84.71,90.0,174.71,0.0,0.0);
  new TRotMatrix("itsrot509","itsrot509",90.0,74.12,90.0,164.12,0.0,0.0);
  new TRotMatrix("itsrot510","itsrot510",90.0,63.53,90.0,153.53,0.0,0.0);  
  new TRotMatrix("itsrot511","itsrot511",90.0,52.94,90.0,142.94,0.0,0.0);
  new TRotMatrix("itsrot512","itsrot512",90.0,42.35,90.0,132.35,0.0,0.0);
  new TRotMatrix("itsrot513","itsrot513",90.0,31.76,90.0,121.76,0.0,0.0); 
  new TRotMatrix("itsrot653","itsrot653",90.0,21.18,90.0,111.18,0.0,0.0); 
  new TRotMatrix("itsrot514","itsrot514",90.0,10.59,90.0,100.59,0.0,0.0);  
  new TRotMatrix("itsrot515","itsrot515",90.0,349.41,90.0,79.41,0.0,0.0);  
  new TRotMatrix("itsrot516","itsrot516",90.0,338.82,90.0,68.82,0.0,0.0);  
  new TRotMatrix("itsrot517","itsrot517",90.0,328.24,90.0,58.24,0.0,0.0);  
  new TRotMatrix("itsrot518","itsrot518",90.0,317.65,90.0,47.65,0.0,0.0);
  new TRotMatrix("itsrot519","itsrot519",90.0,307.06,90.0,37.06,0.0,0.0);
  new TRotMatrix("itsrot520","itsrot520",90.0,296.47,90.0,26.47,0.0,0.0);  
  new TRotMatrix("itsrot521","itsrot521",90.0,285.88,90.0,15.88,0.0,0.0);
  new TRotMatrix("itsrot522","itsrot522",90.0,275.29,90.0,5.29,0.0,0.0);
  new TRotMatrix("itsrot523","itsrot523",90.0,264.71,90.0,354.71,0.0,0.0); 
  new TRotMatrix("itsrot524","itsrot524",90.0,254.12,90.0,344.12,0.0,0.0);  
  new TRotMatrix("itsrot525","itsrot525",90.0,243.53,90.0,333.53,0.0,0.0);  
  new TRotMatrix("itsrot526","itsrot526",90.0,232.94,90.0,322.94,0.0,0.0);  
  new TRotMatrix("itsrot527","itsrot527",90.0,222.35,90.0,312.35,0.0,0.0);  
  new TRotMatrix("itsrot528","itsrot528",90.0,211.76,90.0,301.76,0.0,0.0);
  new TRotMatrix("itsrot618","itsrot618",90.0,201.18,90.0,291.18,0.0,0.0); 
  new TRotMatrix("itsrot529","itsrot529",90.0,190.59,90.0,280.59,0.0,0.0); 
  new TRotMatrix("itsrot533","itsrot533",90.0,180.0,90.0,270.0,0.0,0.0);   
  new TRotMatrix("itsrot530","itsrot530",90.0,169.41,90.0,259.41,0.0,0.0);  
  new TRotMatrix("itsrot531","itsrot531",90.0,158.82,90.0,248.82,0.0,0.0);  
  new TRotMatrix("itsrot501","itsrot501",90.0,148.24,90.0,238.24,0.0,0.0);
  new TRotMatrix("itsrot503","itsrot503",90.0,137.65,90.0,227.65,0.0,0.0);         
  new TRotMatrix("itsrot532","itsrot532",90.0,360.0,90.0,90.0,0.0,0.0);
  new TRotMatrix("itsrot560","itsrot560",90.0,85.26,90.0,175.26,0.0,0.0);  
  new TRotMatrix("itsrot561","itsrot561",90.0,94.74,90.0,184.74,0.0,0.0);
  new TRotMatrix("itsrot562","itsrot562",90.0,104.21,90.0,194.21,0.0,0.0);
  new TRotMatrix("itsrot563","itsrot563",90.0,113.68,90.0,203.68,0.0,0.0); 
  new TRotMatrix("itsrot564","itsrot564",90.0,123.16,90.0,213.16,0.0,0.0);  
  new TRotMatrix("itsrot565","itsrot565",90.0,132.63,90.0,222.63,0.0,0.0);  
  new TRotMatrix("itsrot566","itsrot566",90.0,142.11,90.0,232.11,0.0,0.0);  
  new TRotMatrix("itsrot567","itsrot567",90.0,151.58,90.0,241.58,0.0,0.0);  
  new TRotMatrix("itsrot568","itsrot568",90.0,161.05,90.0,251.05,0.0,0.0);
  new TRotMatrix("itsrot569","itsrot569",90.0,170.53,90.0,260.53,0.0,0.0);
  new TRotMatrix("itsrot533","itsrot533",90.0,180.0,90.0,270.0,0.0,0.0); 
  new TRotMatrix("itsrot534","itsrot534",90.0,189.47,90.0,279.47,0.0,0.0);  
  new TRotMatrix("itsrot535","itsrot535",90.0,198.95,90.0,288.95,0.0,0.0);  
  new TRotMatrix("itsrot623","itsrot623",90.0,208.42,90.0,298.42,0.0,0.0);  
  new TRotMatrix("itsrot537","itsrot537",90.0,217.89,90.0,307.89,0.0,0.0);  
  new TRotMatrix("itsrot538","itsrot538",90.0,227.37,90.0,317.37,0.0,0.0);
  new TRotMatrix("itsrot539","itsrot539",90.0,236.84,90.0,326.84,0.0,0.0);
  new TRotMatrix("itsrot540","itsrot540",90.0,246.32,90.0,336.32,0.0,0.0);  
  new TRotMatrix("itsrot541","itsrot541",90.0,255.79,90.0,345.79,0.0,0.0);
  new TRotMatrix("itsrot542","itsrot542",90.0,265.26,90.0,355.26,0.0,0.0);
  new TRotMatrix("itsrot543","itsrot543",90.0,274.74,90.0,4.74,0.0,0.0); 
  new TRotMatrix("itsrot544","itsrot544",90.0,284.21,90.0,14.21,0.0,0.0);  
  new TRotMatrix("itsrot545","itsrot545",90.0,293.68,90.0,23.68,0.0,0.0);  
  new TRotMatrix("itsrot546","itsrot546",90.0,303.16,90.0,33.16,0.0,0.0);  
  new TRotMatrix("itsrot547","itsrot547",90.0,312.63,90.0,42.63,0.0,0.0);  
  new TRotMatrix("itsrot548","itsrot548",90.0,322.11,90.0,52.11,0.0,0.0);
  new TRotMatrix("itsrot549","itsrot549",90.0,331.58,90.0,61.58,0.0,0.0);
  new TRotMatrix("itsrot550","itsrot550",90.0,341.05,90.0,71.05,0.0,0.0);  
  new TRotMatrix("itsrot551","itsrot551",90.0,350.53,90.0,80.53,0.0,0.0);
  new TRotMatrix("itsrot552","itsrot552",90.0,9.47,90.0,99.47,0.0,0.0);
  new TRotMatrix("itsrot553","itsrot553",90.0,18.95,90.0,108.95,0.0,0.0);
  new TRotMatrix("itsrot620","itsrot620",90.0,28.42,90.0,118.42,0.0,0.0);  
  new TRotMatrix("itsrot555","itsrot555",90.0,37.89,90.0,127.89,0.0,0.0);  
  new TRotMatrix("itsrot556","itsrot556",90.0,47.37,90.0,137.37,0.0,0.0);  
  new TRotMatrix("itsrot557","itsrot557",90.0,56.84,90.0,146.84,0.0,0.0);  
  new TRotMatrix("itsrot558","itsrot558",90.0,66.32,90.0,156.32,0.0,0.0);
  new TRotMatrix("itsrot559","itsrot559",90.0,75.79,90.0,165.79,0.0,0.0);       
  
  
  // --- Define SPD (option 'a') volumes ----------------------------
  
  // SPD - option 'a' 
  // (this is NOT the default)
  
  if (option == 1) { 
    
    dits1[0] = 0.64;
    dits1[1] = ddet1;
    dits1[2] = 3.48;
    new TBRIK("ITS1","ITS1","void",dits1[0],dits1[1],dits1[2]);
    
    dits2[0] = 0.64;
    dits2[1] = ddet2;
    dits2[2] = 3.48;
    new TBRIK("ITS2","ITS2","void",dits2[0],dits2[1],dits2[2]);    
    
    di101[0] = 0.705;
    di101[1] = ddet1;
    di101[2] = 3.536;
    new TBRIK("I101","I101","void",di101[0],di101[1],di101[2]);
    
    di1d1[0] = 0.705;
    di1d1[1] = ddet2;
    di1d1[2] = 3.536;
    new TBRIK("I1D1","I1D1","void",di1d1[0],di1d1[1],di1d1[2]);    
    
    di103[0] = 0.793;
    di103[1] = ddet1+dchip1;
    di103[2] = 3.536;
    new TBRIK("I103","I103","void",di103[0],di103[1],di103[2]);    
    
    di1d3[0] = 0.793;
    di1d3[1] = ddet2+dchip2;
    di1d3[2] = 3.536;
    new TBRIK("I1D3","I1D3","void",di1d3[0],di1d3[1],di1d3[2]);        
        
    di10a[0] = 0.843;
    di10a[1] = ddet1+dchip1+dbus+0.0025;  
    di10a[2] = 19.344;
    new TBRIK("I10A","I10A","void",di10a[0],di10a[1],di10a[2]); 
    
    di20a[0] = 0.843;
    di20a[1] = ddet2+dchip2+dbus+0.0025;  
    di20a[2] = 19.344;
    new TBRIK("I20A","I20A","void",di20a[0],di20a[1],di20a[2]);     

    dits[0] = 3.7;
    dits[1] = 7.7;
    dits[2] = 24;
    dits[3] = 57;
    dits[4] = 100;
    new TTUBS("I12A","I12A","void",dits[0],dits[1],dits[2],dits[3],dits[4]);
     
    dits[0] = 3.7;
    dits[1] = 7.75;
    dits[2] = 26.1;
    new TTUBE("IT12","IT12","void",dits[0],dits[1],dits[2]);  
    
  }
  
  // --- Define SPD (option 'b') volumes ----------------------------
  
  // SPD - option 'b' 
  // (this is the default)

  if (option == 2) {
    
    dits1[0] = 0.64;
    dits1[1] = ddet1;
    dits1[2] = 3.48;
    new TBRIK("ITS1","ITS1","void",dits1[0],dits1[1],dits1[2]);
    
    dits2[0] = 0.64;
    dits2[1] = ddet2;
    dits2[2] = 3.48;
    new TBRIK("ITS2","ITS2","void",dits2[0],dits2[1],dits2[2]);    
    
    di101[0] = 0.705;
    di101[1] = ddet1;
    di101[2] = 3.536;
    new TBRIK("I101","I101","void",di101[0],di101[1],di101[2]);
    
    di1d1[0] = 0.705;
    di1d1[1] = ddet2;
    di1d1[2] = 3.536;
    new TBRIK("I1D1","I1D1","void",di1d1[0],di1d1[1],di1d1[2]);    
    
    di107[0] = 0.793;
    di107[1] = ddet1+dchip1;
    di107[2] = 3.536;
    new TBRIK("I107","I107","void",di107[0],di107[1],di107[2]);    
    
    di1d7[0] = 0.7975;
    di1d7[1] = ddet2+dchip2;
    di1d7[2] = 3.536;
    new TBRIK("I1D7","I1D7","void",di1d7[0],di1d7[1],di1d7[2]);        
        
    di10b[0] = 0.843;
    di10b[1] = ddet1+dchip1+dbus+0.0025;  
    di10b[2] = 19.344;
    new TBRIK("I10B","I10B","void",di10b[0],di10b[1],di10b[2]); 
    
    di20b[0] = 0.843;
    di20b[1] = ddet2+dchip2+dbus+0.0025;  
    di20b[2] = 19.344;
    new TBRIK("I20B","I20B","void",di20b[0],di20b[1],di20b[2]);     

    dits[0] = 3.7;
    dits[1] = 7.7;
    dits[2] = 24;
    dits[3] = 57;
    dits[4] = 100;
    new TTUBS("I12B","I12B","void",dits[0],dits[1],dits[2],dits[3],dits[4]);
     
    dits[0] = 3.7;
    dits[1] = 7.75;
    dits[2] = 26.1;
    new TTUBE("IT12","IT12","void",dits[0],dits[1],dits[2]);  
         

  }

  // --- Define SDD volumes ------------------------------------------
  
  TPCON *it34 = new TPCON("IT34","IT34","void",0.,360.,6); 
  it34->DefineSection(0,-34.6,23.49,28.); 
  it34->DefineSection(1,-23.65,23.49,28.); 
  it34->DefineSection(2,-23.65,14.59,28.); 
  it34->DefineSection(3,23.65,14.59,28.); 
  it34->DefineSection(4,23.65,23.49,28.); 
  it34->DefineSection(5,34.6,23.49,28.);   
  
  I302dits[0] = 3.6250;
  I302dits[1] = 0.0150;
  I302dits[2] = 4.3794; 
  new TBRIK("I302","I302","void",I302dits[0],I302dits[1],I302dits[2]);

  I004dits[0] = I302dits[0]+0.005;
  I004dits[1] = 2*I302dits[1]+Y_SDD_sep/2.;
  I004dits[2] = TMath::Abs(Z_SDD_lay3[0]);
  if (I004dits[2] < TMath::Abs(Z_SDD_lay3[5])) {
    I004dits[2] = TMath::Abs(Z_SDD_lay3[5]);
  }
  I004dits[2] = I004dits[2] + I302dits[2];  
  new TBRIK("I004","I004","void",I004dits[0],I004dits[1],I004dits[2]); 
  
  dits[0] = 3.50850;
  dits[1] = 0.01499; 
  dits[2] = 3.76320;  
  new TBRIK("ITS3","ITS3","void",dits[0],dits[1],dits[2]);    
 
  I402dits[0] = 3.6250;
  I402dits[1] = 0.0150;
  I402dits[2] = 4.3794; 
  new TBRIK("I402","I402","void",I402dits[0],I402dits[1],I402dits[2]);

  I005dits[0] = I402dits[0]+0.005;
  I005dits[1] = 2*I402dits[1]+Y_SDD_sep/2.;
  I005dits[2] = TMath::Abs(Z_SDD_lay4[0]);
  if (I005dits[2] < TMath::Abs(Z_SDD_lay4[7])) {
    I005dits[2] = TMath::Abs(Z_SDD_lay4[7]);
  }
  I005dits[2] = I005dits[2] + I402dits[2];  
  new TBRIK("I005","I005","void",I005dits[0],I005dits[1],I005dits[2]);   

  dits[0] = 3.50850;
  dits[1] = 0.01499; 
  dits[2] = 3.76320;
  new TBRIK("ITS4","ITS4","void",dits[0],dits[1],dits[2]);

  
  // --- Define SSD volumes ------------------------------------------
  

  TPCON *it56 = new TPCON("IT56","IT56","void",0.,360.,6); 
  it56->DefineSection(0,-57.45,43.6,48.); 
  it56->DefineSection(1,-49.15,43.6,48.); 
  it56->DefineSection(2,-49.15,36.9,48.); 
  it56->DefineSection(3,50.55,36.9,48.); 
  it56->DefineSection(4,50.55,43.6,48.); 
  it56->DefineSection(5,57.45,43.6,48.);    

  dits[0] = 3.75;
  dits[1] = 0.045;
  dits[2] = 43.3;
  new TBRIK("I565","I565","void",dits[0],dits[1],dits[2]);  

  dits[0] = 3.75;
  dits[1] = 0.045;
  dits[2] = 50.975;
  new TBRIK("I569","I569","void",dits[0],dits[1],dits[2]);  
  
  dits[0] = 3.75;
  dits[1] = 0.015;
  dits[2] = 2.1;
  new TBRIK("I562","I562","void",dits[0],dits[1],dits[2]);	
  
  dits[0] = 3.75;
  dits[1] = 0.015;
  dits[2] = 2.1;
  new TBRIK("I566","I566","void",dits[0],dits[1],dits[2]);	  

  dits[0] = 3.65;
  dits[1] = 0.015;
  dits[2] = 2;
  new TBRIK("ITS5","ITS5","void",dits[0],dits[1],dits[2]); 

  dits[0] = 3.65;
  dits[1] = 0.015;
  dits[2] = 2;
  new TBRIK("ITS6","ITS6","void",dits[0],dits[1],dits[2]);  

  //
  
  //top->cd();


  // --- Place SSD volumes into their mother volume    

    // Place IT56 in Alice
    node = new TNode("IT56","IT56","IT56",0.,0.,0.,"");
    node->SetLineColor(kColorITS);
    node->SetVisibility(0);
    node->cd();
       //
       // Place copy #1 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",0.,38.445,0.,"");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #2 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-7.0924,37.9412,0.,"itsrot514");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #3 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-13.8879,35.8489,0.,"itsrot653");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #4 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-20.3195,32.817,0.,"itsrot513");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #5 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-25.9002,28.4112,0.,"itsrot512");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #6 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-30.8022,23.2608,0.,"itsrot511");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #7 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-34.4146,17.1364,0.,"itsrot510");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #8 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-37.1249,10.563,0.,"itsrot509");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #9 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-38.281,3.5473,0.,"itsrot508");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #10 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-38.4338,-3.5614,0.,"itsrot507");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #11 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-36.9774,-10.521,0.,"itsrot506");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #12 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-34.5519,-17.2048,0.,"itsrot505");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #13 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-30.6798,-23.1683,0.,"itsrot504");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #14 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-26.0036,-28.5246,0.,"itsrot503");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #15 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-20.2387,-32.6866,0.,"itsrot501");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #16 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-13.9434,-35.992,0.,"itsrot531");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #17 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",-7.0643,-37.7904,0.,"itsrot530");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #18 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",0.,-38.5984,0.,"itsrot533");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #19 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",7.0642,-37.7904,0.,"itsrot529");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #20 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",13.9433,-35.992,0.,"itsrot618");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #21 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",20.2387,-32.6866,0.,"itsrot528");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #22 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",26.0036,-28.5246,0.,"itsrot527");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #23 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",30.6798,-23.1683,0.,"itsrot526");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #24 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",34.5519,-17.2048,0.,"itsrot525");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #25 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",36.9774,-10.521,0.,"itsrot524");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #26 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",38.4338,-3.5614,0.,"itsrot523");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #27 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",38.281,3.5472,0.,"itsrot522");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #28 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",37.125,10.5629,0.,"itsrot521");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #29 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",34.4146,17.1364,0.,"itsrot520");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #30 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",30.8022,23.2607,0.,"itsrot519");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #31 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",25.9002,28.4112,0.,"itsrot518");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #32 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",20.3195,32.817,0.,"itsrot517");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #33 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",13.888,35.8489,0.,"itsrot516");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #34 of I565 in IT56
       //
       sub1node = new TNode("I565","I565","I565",7.0925,37.9412,0.,"itsrot515");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,41.1546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,37.2246,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,33.3146,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,29.3846,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,25.4746,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,21.5446,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,17.6346,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,13.7046,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,9.7946,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,5.8645,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,1.9546,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-1.9754,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-5.8855,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-9.8154,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-13.7254,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-17.6555,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-21.5655,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-25.4954,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-29.4054,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-33.3354,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,0.03,-37.2454,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS5 directly in I565
          //
	  sub2node = new TNode("ITS5","ITS5","ITS5",0.,-0.03,-41.1554,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #1 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-14.139,41.1856,0.,"itsrot553");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #2 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-20.7978,38.431,0.,"itsrot620");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #3 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-26.7459,34.3631,0.,"itsrot555");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #4 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-32.1494,29.5956,0.,"itsrot556");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #5 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-36.4544,23.8169,0.,"itsrot557");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #6 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-40.0172,17.5532,0.,"itsrot558");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #7 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-42.2125,10.6897,0.,"itsrot559");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #8 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-43.5484,3.6085,0.,"itsrot560");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #9 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-43.3963,-3.5959,0.,"itsrot561");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #10 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-42.3606,-10.7271,0.,"itsrot562");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #11 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-39.8773,-17.4918,0.,"itsrot563");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #12 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-36.5823,-23.9004,0.,"itsrot564");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #13 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-32.0371,-29.4922,0.,"itsrot565");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #14 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-26.8397,-34.4836,0.,"itsrot566");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #15 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-20.7251,-38.2967,0.,"itsrot567");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #16 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-14.1886,-41.33,0.,"itsrot568");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #17 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-7.1673,-42.9511,0.,"itsrot569");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #18 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",0.,-43.6977,0.,"itsrot533");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #19 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",7.1673,-42.9511,0.,"itsrot534");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #20 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",14.1886,-41.33,0.,"itsrot535");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #21 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",20.7251,-38.2967,0.,"itsrot623");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #22 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",26.8397,-34.4836,0.,"itsrot537");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #23 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",32.0371,-29.4922,0.,"itsrot538");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #24 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",36.5822,-23.9004,0.,"itsrot539");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #25 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",39.8773,-17.4918,0.,"itsrot540");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #26 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",42.3606,-10.7272,0.,"itsrot541");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #27 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",43.3963,-3.5959,0.,"itsrot542");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #28 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",43.5484,3.6085,0.,"itsrot543");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #29 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",42.2125,10.6897,0.,"itsrot544");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #30 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",40.0172,17.5532,0.,"itsrot545");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #31 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",36.4544,23.8169,0.,"itsrot546");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #32 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",32.1494,29.5956,0.,"itsrot547");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #33 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",26.7459,34.3631,0.,"itsrot548");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #34 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",20.7978,38.431,0.,"itsrot549");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #35 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",14.139,41.1856,0.,"itsrot550");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #36 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",7.1924,43.1017,0.,"itsrot551");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #37 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",0.,43.545,0.,"");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();
       //
       // Place copy #38 of I569 in IT56
       //
       sub1node = new TNode("I569","I569","I569",-7.1924,43.1017,0.,"itsrot552");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,46.9203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();
	  //
	  // Place copy #2 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,43.0103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #3 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,39.1003,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #4 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,35.1903,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #5 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,31.2803,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #6 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,27.3703,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #7 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,23.4603,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #8 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,19.5503,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #9 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,15.6403,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #10 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,11.7303,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #11 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,7.8203,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #12 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,3.9103,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #13 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,0.0003,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);  
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #14 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-3.9097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #15 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-7.8197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #16 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-11.7297,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #17 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-15.6397,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #18 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-19.5497,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #19 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-23.4597,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);   
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #20 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-27.3697,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #21 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-31.2797,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #22 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-35.1897,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #23 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-39.0997,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #24 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,0.03,-43.0097,"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1); 
          //fNodes->Add(sub2node);
	  sub1node->cd();
	  //
	  // Place copy #25 of ITS6 in I569
          //
	  sub2node = new TNode("ITS6","ITS6","ITS6",0.,-0.03,-46.9197,"itsrot532");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
	  sub1node->cd();	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  
       //fNodes->Add(sub1node);
       node->cd();



    //fNodes->Add(node);

    c1_1->cd();
    c1_1->SetTheta(90);
    c1_1->SetPhi(0);
    node->Draw();
    
    c1_2->cd();
    c1_2->SetTheta(0);
    c1_2->SetPhi(0);   
    node->Draw();
    
    c1_3->cd();
    c1_3->SetTheta(0);
    c1_3->SetPhi(90); 
    node->Draw();    
        
    c1_4->cd();
    c1_4->SetTheta(20);
    c1_4->SetPhi(140);        
    node->Draw();

    c1->Update();

  // --- Place SDD volumes into their mother volume 

    // Place IT34 in Alice
    node = new TNode("IT34","IT34","IT34",0.,0.,0.,"");
    node->SetLineColor(kColorITS);
    node->SetVisibility(0);
    node->cd();
       //
       // Place copy #1 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-3.2777,14.3607,0.,"itsrot321");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #2 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-9.5581,11.9855,0.,"itsrot333");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #3 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-13.2713,6.3911,0.,"itsrot336");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #4 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-15.33,0.,0.,"itsrot350");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #5 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-13.2713,-6.3911,0.,"itsrot313");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #6 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-9.5581,-11.9855,0.,"itsrot311");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #7 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",-3.2777,-14.3607,0.,"itsrot310");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #8 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",3.4112,-14.9456,0.,"itsrot386");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #9 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",9.184,-11.5164,0.,"itsrot309");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #10 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",13.8119,-6.6514,0.,"itsrot308");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #11 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",14.73,0.,0.,"itsrot356");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #12 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",13.8119,6.6514,0.,"itsrot307");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #13 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",9.184,11.5164,0.,"itsrot306");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #14 of I004 in IT34
       //
       sub1node = new TNode("I004","I004","I004",3.4113,14.9456,0.,"itsrot305");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,ySDD,Z_SDD_lay3[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS3 directly in I004
          //
	  ySDD = Y_SDD_sep/2.+I302dits[1];
	  sub2node = new TNode("ITS3","ITS3","ITS3",0.,-ySDD,Z_SDD_lay3[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #1 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-3.3629,23.3895,-0.15,"itsrot335");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #2 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-10.0447,21.9949,-0.15,"itsrot332");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #3 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-15.4744,17.8584,-0.15,"itsrot331");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #4 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-20.3415,13.0727,-0.15,"itsrot366");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #5 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-22.6728,6.6573,-0.15,"itsrot330");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #6 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-24.18,0.,-0.15,"itsrot350");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #7 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-22.6728,-6.6573,-0.15,"itsrot329");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #8 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-20.3415,-13.0727,-0.15,"itsrot328");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #9 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-15.4744,-17.8584,-0.15,"itsrot327");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #10 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-10.0447,-21.9949,-0.15,"itsrot326");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #11 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",-3.3629,-23.3895,-0.15,"itsrot325");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #12 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",3.4412,-23.9339,-0.15,"itsrot324");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #13 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",9.8163,-21.4946,-0.15,"itsrot323");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #14 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",15.8345,-18.274,-0.15,"itsrot322");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #15 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",19.8788,-12.7753,-0.15,"itsrot320");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #16 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",23.2005,-6.8123,-0.15,"itsrot319");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #17 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",23.63,0.,-0.15,"itsrot318");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #18 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",23.2005,6.8123,-0.15,"itsrot317");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #19 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",19.8788,12.7753,-0.15,"itsrot316");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #20 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",15.8345,18.274,-0.15,"itsrot315");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #21 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",9.8163,21.4946,-0.15,"itsrot314");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       //
       // Place copy #22 of I005 in IT34
       //
       sub1node = new TNode("I005","I005","I005",3.4412,23.9339,-0.15,"itsrot334");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();
          //
          // Place copy #1 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[0],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #2 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[1],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #3 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[2],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #4 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[3],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #5 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[4],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #6 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[5],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #7 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,ySDD,Z_SDD_lay4[6],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
          //
          // Place copy #8 of ITS4 directly in I005
          //
	  ySDD = -(Y_SDD_sep/2.+I402dits[1]);
	  sub2node = new TNode("ITS4","ITS4","ITS4",0.,-ySDD,Z_SDD_lay4[7],"");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(1);
          //fNodes->Add(sub2node);
          sub1node->cd();        
       //fNodes->Add(sub1node);
       node->cd();       
       
              
    //fNodes->Add(node);
             
    //node->Draw("SAME");
    //c1->Update();


  // --- Place SPD (option 'a') volumes into their mother volume 
  
  // SPD - option 'a' 
  // (this is NOT the default)

  if (option == 1) {


  }
  
  

  // --- Place SPD (option 'b') volumes into their mother volume 
  
  // SPD - option 'b' 
  // (this is the default)

  if (option == 2) { 
  
    // Place IT12 in Alice
    //
    node = new TNode("IT12","IT12","IT12",0.,0.,0.,"");
    node->SetLineColor(kColorITS);
    node->SetVisibility(0);
    node->cd();    
       //
       // Place copy #1 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #2 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot245");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #3 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot234");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #4 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot246");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #5 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot247");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #6 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot236");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #7 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot244");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #8 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot233");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #9 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot248");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       //
       // Place copy #10 of I12B in IT12
       //
       sub1node = new TNode("I12B","I12B","I12B",0.,0.,0.,"itsrot249");
       sub1node->SetLineColor(kColorITS);
       sub1node->SetVisibility(0);
       sub1node->cd();    
          //
	  // Place copy #1 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(252.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(252.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",1.4531+deltax,3.8152+deltay,0.,"itsrot239");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I10B in I12B
	  //
	  deltax=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Cos(270.*TMath::Pi()/180.);  
          deltay=((ddet1-0.01/2.)+(dchip1-0.015/2.))*TMath::Sin(270.*TMath::Pi()/180.);
          sub2node = new TNode("I10B","I10B","I10B",0.203+deltax,3.8206+deltay,0.,"itsrot238");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I107 in I10B
             //
	     sub3node = new TNode("I107","I107","I107",-0.0455,-di10b[1]+di107[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I101 in I107
                //
		sub4node = new TNode("I101","I101","I101",0.,ddet1,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS1 in I101
	           //
		   sub5node = new TNode("ITS1","ITS1","ITS1",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #1 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(40.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(40.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",3.0174+deltax,6.5143+deltay,0.,"itsrot240");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #2 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(49.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(49.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",1.9612+deltax,6.9062+deltay,0.,"itsrot241");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #3 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(58.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(58.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",0.8567+deltax,7.1279+deltay,0.,"itsrot242");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
          //
	  // Place copy #4 of I20B in I12B
	  //
	  deltax=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Cos(67.*TMath::Pi()/180.);  
          deltay=((ddet2-0.01/2.)+(dchip2-0.015/2.))*TMath::Sin(67.*TMath::Pi()/180.);
          sub2node = new TNode("I20B","I20B","I20B",-0.2689+deltax,7.1742+deltay,0.,"itsrot243");
          sub2node->SetLineColor(kColorITS);
          sub2node->SetVisibility(0);
	  sub2node->cd();
             //
	     // Place copy #1 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
		// Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();
		   //		    
	           // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
	     sub2node->cd(); 
	     //
	     // Place copy #2 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #3 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-3.536,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	     //
	     // Place copy #4 of I1D7 in I20B
             //
	     sub3node = new TNode("I1D7","I1D7","I1D7",-0.0455,-di20b[1]+di1d7[1],-10.708,"");
             sub3node->SetLineColor(kColorITS);
             sub3node->SetVisibility(0);
	     sub3node->cd();
	        //
	        // Place copy #1 of I1D1 in I1D7
                //
		sub4node = new TNode("I1D1","I1D1","I1D1",0.,ddet2,0.,"");
                sub4node->SetLineColor(kColorITS);
                sub4node->SetVisibility(0);
	        sub4node->cd();		    
	           //
		   // Place copy #1 of ITS2 in I1D1
	           //
		   sub5node = new TNode("ITS2","ITS2","ITS2",0.,0.,0.,"");
                   sub5node->SetLineColor(kColorITS);                   
                   //fNodes->Add(sub5node);
	           sub4node->cd();   
	        //fNodes->Add(sub4node);  
	     sub3node->cd(); 
	     //fNodes->Add(sub3node);
             sub2node->cd(); 
	  //fNodes->Add(sub2node);   	
          sub1node->cd(); 
       //fNodes->Add(sub1node);
       node->cd(); 
       
    //fNodes->Add(node)
    
    //node->Draw("SAME");
    //c1->Update();

  } 
                  
	// Gets Tree of RecPoints...
	TTree *TR = gAlice->TreeR();
	if (!TR) {
		cout << "A problem occurred when getting the TreeR!!!" << endl;
		return;
	}
        AliITSRecPoint *recp=0;	
        TClonesArray  *recPoints = ITS->RecPoints();

	// Gets found tracks...
       const char *filename="itstracks.root";	
       TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
       if (!file) file = new TFile(filename);

       char tname[30];   

//  get tracks from ITSiotrack.root
  Int_t nev=0;	
  sprintf(tname,"TreeT%d",nev);
   TTree *tracktree=(TTree*)file->Get(tname);
   TBranch *tbranch=tracktree->GetBranch("ITStracks");
   //cout<<" nev = "<<nev<<"\n";   
   //cout<<" open the file \n"; 
	
   Int_t nentry=tracktree->GetEntries();
	
    Float_t idpoint[6],idmodule[6]; 

   TObjArray tarray(nentry);
   printf("Found tracks = %d\n",nentry);
	
   for (Int_t i=0; i<nentry; i++) {
      AliITSIOTrack *iotrack=new AliITSIOTrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
	tarray.AddLast(iotrack);
   }		 
	

	  AliITSIOTrack *iotrack;
	  Float_t global[3], local[3];	
	  Float_t idpoint[6], idmodule[6]; 
	  Float_t xp[7], yp[7], zp[7];	
	  Float_t pxtr, pytr, pztr;
	  Float_t x0tr, y0tr, z0tr;	  	   
   for (Int_t i=0; i<nentry; i++) {
     AliITSIOTrack *iotrack=new AliITSIOTrack;      	
	 iotrack=(AliITSIOTrack*)tarray.UncheckedAt(i);
	 if(!iotrack) continue;
     Int_t labITS=iotrack->GetLabel();
            pxtr=iotrack->GetPx();   
	    pytr=iotrack->GetPy();   
	    pztr=iotrack->GetPz();   
            x0tr=iotrack->GetX();   
	    y0tr=iotrack->GetY();   
	    z0tr=iotrack->GetZ();	    
     //cout << labITS << " ok" << endl;
    // Int_t labTPC=iotrack->GetTPCLabel();
		xp[0] = xp[1] = xp[2] = xp[3] = xp[4] = xp[5] = xp[6] = 0.0;
		yp[0] = yp[1] = yp[2] = yp[3] = yp[4] = yp[5] = yp[6] = 0.0;
		zp[0] = zp[1] = zp[2] = zp[3] = zp[4] = zp[5] = zp[6] = 0.0;

	  Int_t np=0;				  


          xp[0]=x0tr;
	  xp[0]=x0tr;
          zp[0]=x0tr;
	  
	  for(Int_t ipp=0; ipp<6; ipp++){
	  idpoint[ipp]= iotrack->GetIdPoint(ipp);
	  idmodule[ipp] = iotrack->GetIdModule(ipp);
	  //cout << idmodule[ipp] << " " << idpoint[ipp] << endl;
	 if(idmodule[ipp]>=0){  
	 ITS->ResetRecPoints(); 	 
    gAlice->TreeR()->GetEvent(idmodule[ipp]); 
     Int_t npoints=recPoints->GetEntries();	
	 recp =(AliITSRecPoint*)recPoints->UncheckedAt(idpoint[ipp]);
	  local[0]=recp->GetX();
	  local[1]=0.;
	  local[2]= recp->GetZ();
	  //if(i==32) cout<<" local ="<<local[0]<<" "<<local[1]<<" "<<local[2]<<"\n";
      gm->LtoG(Int_t(idmodule[ipp]),local,global);
	  //if(i==32) cout<<" global ="<<global[0]<<" "<<global[1]<<" "<<global[2]<<"\n";		
		xp[np+1]=global[0];
		yp[np+1]=global[1];
		zp[np+1]=global[2];
	   np++;			
	  } 	  	 	  
	  }
	  
	  
		
		TPolyLine3D *plxy = new TPolyLine3D(np+1, xp, yp, zp);
		plxy->SetLineColor(kBlack);

                c1_1->cd(); 
		plxy->Draw();

		
	        c1_2->cd(); 
		plxy->Draw();

		
		c1_3->cd(); 
		plxy->Draw();

		
	        c1_4->cd(); 
		plxy->Draw();

		c1->Update();
		
				
		//cout<<" i  np = "<<i<<"  "<<np<<"\n";getchar();
		/*
		TPolyLine *plxz = new TPolyLine(6, x, z);
		plxz->SetLineColor(kRed);
		plxz->Draw();
		cxz->Update();			  
      */
	  }




}

