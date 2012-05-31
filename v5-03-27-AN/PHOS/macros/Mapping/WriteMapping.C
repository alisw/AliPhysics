WriteMapping()
{
  // Create the RCU mapping files for PHOS
  // The geometrical mapping within one FEE card is read
  // from the file CSP2ALTRO.dat prepared manually beforehand.
  //
  // The hardware address of the FEE channels is a 12-bit word formed 
  // from the branch number, FEE number, ALTRO chip number and ALTRO channel
  // as follows:
  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // |0|0|0|0| |       |     |       |
  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  //          ^   FEE    chip  channel
  //       branch
  //
  // Max address: 1 1110 100 1111 = 3919
  // 
  // Author: Yuri Kharlov
  // Date  : 29 October 2009
  // $Id$


  char string[128];
  UInt_t xcell,zcell,csp,altro,chanHG,chanLG;
  UInt_t hwAddress, maxHWAddress[4], nHWaddress[4];
  const char *FEEmapping[2] = {"CSP2ALTRO_new.dat", "CSP2ALTRO_old.dat"};
  Int_t map=0;

  for (Int_t module=0; module<5; module++) {
    printf("\n\n====> Mapping file is created for the module %d\n",module);
    if      (module <4) map=0;
    else if (module==4) map=1;
    FILE *fd = fopen(FEEmapping[map],"r");
    if (fd != 0) 
      printf("Input file %s is opened successfully\n",FEEmapping[map]);
    else {
      printf("Cannot open the input file %s\n",FEEmapping[map]);
      return -1;
    }
    
    FILE *fRCU[4];

    for (Int_t iRCU=0; iRCU<4; iRCU++) {
      TString rcuFileName = Form("Mod%dRCU%d.data.unsorted",module,iRCU);
      fRCU[iRCU] = fopen(rcuFileName.Data(),"w");
      maxHWAddress[iRCU]=0;
      nHWaddress[iRCU]=0;
    }
    
    while (fgets(string,128,fd)) {
      if (string[0]=='*') {
	continue;
      }
      sscanf(string,"%d %d %d %d %d %d",
	     &xcell,&zcell,&csp,&altro,&chanHG,&chanLG);
      for (Int_t iRCU=0; iRCU<4; iRCU++) {
	for (Int_t iBranch=0; iBranch<2; iBranch++) {
	  for (Int_t iFEE=1; iFEE<=14; iFEE++) {
	    // High gain
	    hwAddress = chanHG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	    if (hwAddress > maxHWAddress[iRCU]) maxHWAddress[iRCU]=hwAddress;
	    nHWaddress[iRCU]++;
	    fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",
		    hwAddress,
		    xcell+iRCU*16,
		    zcell+27+(1-2*iBranch)+(iFEE-1)*2*(1-2*iBranch),1);
	    // Low gain
	    hwAddress = chanLG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	    if (hwAddress > maxHWAddress[iRCU]) maxHWAddress[iRCU]=hwAddress;
	    nHWaddress[iRCU]++;
	    fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",
		    hwAddress,
		    xcell+iRCU*16,
		    zcell+27+(1-2*iBranch)+(iFEE-1)*2*(1-2*iBranch),0);
	  }
	}
      }
    }
    printf("End of input file\n");
    fclose(fd);

    for (Int_t iRCU=0; iRCU<4; iRCU++) {
      fclose(fRCU[iRCU]);
    }
    
    // Post-process the RCU mapping files
    
    for (Int_t iRCU=0; iRCU<4; iRCU++) {
      
      // Add the number of channels and maximum HW address
      
      TString rcuFileName = Form("Mod%dRCU%d.data",module,iRCU);
      fRCU[iRCU] = fopen(rcuFileName.Data(),"w");
      fprintf(fRCU[iRCU],"%d\n%d\n",nHWaddress[iRCU]+256,maxHWAddress[iRCU]);
      
      // TRU mapping
      
      for (int i=0; i<128; i++) {
	fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",i,iRCU,i,2);
      }
      for (int i=0; i<128; i++) {
	fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",i+2048,iRCU,i+2048,2);
      }
      fclose(fRCU[iRCU]);
      
      // Sort HW addresses
      
      TString cmd = Form("sort -n Mod%dRCU%d.data.unsorted >> %s", 
			 module,iRCU,rcuFileName.Data());
      gSystem->Exec(cmd);
      
      cmd = Form("rm -f Mod%dRCU%d.data.unsorted", module,iRCU);
      gSystem->Exec(cmd);
      
      printf("RCU mapping file %s is created\n",rcuFileName.Data());
    }
  }
}
