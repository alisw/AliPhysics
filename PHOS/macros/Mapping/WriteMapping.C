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
  // Date  : 15 August 2009

  char string[128];
  UInt_t xcell,zcell,csp,altro,chanHG,chanLG;
  UInt_t hwAddress, maxHWAddress[4]={0,0,0,0}, nHWaddress[4]={0,0,0,0};
  FILE *fd = fopen("CSP2ALTRO.dat","r");
  if (fd != 0) 
    printf("Input file CSP2ALTRO.dat is opened successfully\n");
  else {
    printf("Cannot open the input file CSP2ALTRO.dat\n");
    return -1;
  }

  FILE *fRCU[4];

  for (Int_t iRCU=0; iRCU<4; iRCU++) {
    TString rcuFileName = Form("RCU%d.data.unsorted",iRCU);
    fRCU[iRCU] = fopen(rcuFileName.Data(),"w");
  }

  while (fgets(string,128,fd)) {
    if (string[0]=='*') {
//       printf("%s",string);
      continue;
    }
    sscanf(string,"%d %d %d %d %d %d",&xcell,&zcell,&csp,&altro,&chanHG,&chanLG);
//     printf("%2d %2d %2d %2d %2d %2d\n",
// 	     xcell,zcell,csp,altro,chanHG,chanLG);
    for (Int_t iRCU=0; iRCU<4; iRCU++) {
      for (Int_t iBranch=0; iBranch<2; iBranch++) {
	for (Int_t iFEE=1; iFEE<=14; iFEE++) {

	  hwAddress = chanHG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	  if (hwAddress > maxHWAddress[iRCU]) maxHWAddress[iRCU]=hwAddress;
	  nHWaddress[iRCU]++;
	  fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",
 		  hwAddress,xcell+iRCU*16,zcell+27+(2*iBranch-1)+(iFEE-1)*2*(2*iBranch-1),1);

	  hwAddress = chanLG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	  if (hwAddress > maxHWAddress[iRCU]) maxHWAddress[iRCU]=hwAddress;
	  nHWaddress[iRCU]++;
	  fprintf(fRCU[iRCU],"%4d %4d %4d %4d\n",
 		  hwAddress,xcell+iRCU*16,zcell+27+(2*iBranch-1)+(iFEE-1)*2*(2*iBranch-1),1);
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

    TString rcuFileName = Form("RCU%d.data",iRCU);
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

    TString cmd = Form("sort -n RCU%d.data.unsorted >> %s", iRCU,rcuFileName.Data());
    gSystem->Exec(cmd);

    cmd = Form("rm -f RCU%d.data.unsorted", iRCU);
    gSystem->Exec(cmd);

    printf("RCU mapping file %s is created\n",rcuFileName.Data());
  }

}
