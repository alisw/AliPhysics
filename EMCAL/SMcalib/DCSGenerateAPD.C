// constants
static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
static const int fgkEmCalCols = 48; // number of columns per module for EMCAL

const int NRCU = 2; // per SM
const int NBranch = 2; // per RCU
const int NFEC = 9; // per branch, labelled 1..9
const int NCSP = 32; // per FEC

// conversion between DAC (0-0x3ff) and HV values (V):
// hv = hvmin + prop*DAC; values from PHOS manual, and used in Houston/Catania
const float hvmin = 209.9;
const float prop = 0.2022; 

// some global variables
Float_t biasVoltage[NRCU][NBranch][NFEC][NCSP]; 
int towerCol[NRCU][NBranch][NFEC][NCSP]; 
int towerRow[NRCU][NBranch][NFEC][NCSP]; 

//__________________________________________________________
void Tower2FEEBiasInfo(const char *inputFileName)
{
  ifstream inputFile(inputFileName);
  int ic, ir;
  Float_t ival;
  int ircu, ibranch, card, icsp;
  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {
      inputFile >> ic >> ir >> ival;

      // could check here that ic && ir match with icol && irow, but should not be needed

      // translate to FEE type indices
      Tower2FEEMap(ic, ir, 
		   &ircu, &ibranch, &card, &icsp);

      // debug
      /*
      printf("ic %d ir %d ircu %d ibranch %d card %d icsp %d\n",
	     ic, ir, ircu, ibranch, card, icsp);
      */

      // store value
      biasVoltage[ircu][ibranch][card][icsp] = ival;
      towerCol[ircu][ibranch][card][icsp] = ic;
      towerRow[ircu][ibranch][card][icsp] = ir;
    }
  }

  inputFile.close();

  return;
}

//__________________________________________________________
void Tower2FEEMap(const int icol, const int irow,
		  int *ircu, int *ibranch, int *card, int *icsp)
{ /*
    If you are interested in where these magic numbers come from -
    See mapping info on 
    http://dsilverm.web.cern.ch/dsilverm/mapping/emcal_mapping.html
    http://dsilverm.web.cern.ch/dsilverm/mapping/ppt/Coordinates_and_Mapping.pdf
   */

  // each FEC covers a 4x8 tower area
  int C = irow/8; // Cable bundle
  int FEC = C*12 + icol/4; // FEC in 0..35 range

  *ircu = FEC / 18; // 18 FEC per RCU
  *ibranch = (FEC%18) / 9;
  *card = FEC % 9;  

  // columns and rows within an FEC area
  int tCol = icol%4;
  int tRow = irow%8;

  // The mapping to CSP is a bit complicated so I also define two more help variables here..
  // which T-card?
  int TCard = tCol/2; // 0=Top (even StripModules), 1=Bottom (odd StripModules)
  int locCol = tCol%2;  // local column inside T-card 

  *icsp = (7 - tRow) + locCol*16 + TCard*8;
}

/* Main method.. */
//__________________________________________________________
void DCSGenerateAPD(const char *inputFileName,
		    const char *outputDir,
		    const int readBack=0,
		    const int RCUFWVersion=1, // default as of May 2009 is RCU FWv2
		    const int FEEBCVersion=1) // default as of May 2009 is PCM4
{

  // set up which bias voltage should be applicable for which CSP..
  Tower2FEEBiasInfo(inputFileName);

  // general setup block: try to keep magic numbers/registers collected here
  const char *branch_str[] = { "A", "B"};
  const int nRCUFW = 2;
  const char *rcufw_str[] = { "v1", "v2"};
  const int rcu_addr_start[nRCUFW] = {0x7000, 0x0000}; 
  const int read_header[nRCUFW] = {0x520000, 0x020000}; 
  const int write_header[nRCUFW] = {0x620000, 0x220000}; 
  const int voltage_value_start[nRCUFW] = {0x700000, 0x200000}; 
  const int rcu_endmem[nRCUFW] = {0x390000, 0x3F0000}; 
  const int rcu_exeseq[nRCUFW] = {0x0, 0x5304}; 
  const int rcu_result_reg[nRCUFW] = {0x6000, 0x2000}; 

  const int nFEEBC = 2;
  const char *feebc_str[] = { "PCM3", "PCM4"};
  const int trailer_offset[nFEEBC] = {0x48, 0x68};
  const int update_voltage_register = 0x1E;
  const int max_dac_value = 0x3FF;

  printf("DCSGenerateAPD: RCU FW %s, FEE BC %s\n", 
	 rcufw_str[RCUFWVersion], feebc_str[FEEBCVersion]);

  // resulting voltage settings should be good within a few volts 
  cout << " HV-DAC prop. constant = " << prop << endl;
  char iv_dac_setting[100]; 
  
  char cfile[200];

  FILE* fout_setbias_card[NRCU][NBranch][NFEC];
  FILE* fout_readbias_card[NRCU][NBranch][NFEC];
  
  // end of setup, let's go..
  
  int rcu_addr_card = rcu_addr_start[RCUFWVersion];
  int csp_addr = trailer_offset[RCUFWVersion];
  int word = 0;
  char comment[400];
  
  int rcu_addr_read = rcu_addr_start[RCUFWVersion]; // we'll also write the readbias file in the same loop, so
  // need a separate index also

  for (int rcu=0; rcu<NRCU; rcu++) {
    for (int branch=0; branch<NBranch; branch++) {
      for (int ifec=0; ifec<NFEC; ifec++) {
	int card = ifec;
	int icard = ifec+1;

	sprintf(cfile,"%s/set_rcu_%d_bias_branch_%s_FEC_%d.scr",
		outputDir, rcu, 
		branch_str[branch], icard);
	fout_setbias_card[rcu][branch][card] = fopen(cfile, "w");

	sprintf(cfile,"%s/read_rcu_%d_bias_branch_%s_FEC_%d.scr", 
		outputDir, rcu,
		branch_str[branch], icard);
	fout_readbias_card[rcu][branch][card] = fopen(cfile, "w");

	rcu_addr_card = rcu_addr_start[RCUFWVersion];
	rcu_addr_read = rcu_addr_start[RCUFWVersion];

	for (int icsp = 0; icsp<NCSP; icsp++) {
	  
	  /* 
	     some funkiness to address the CSPs correctly follows here. 
	     DS verified this with section 16.1 "Bias voltage programming", table 8
	     of H. Muller's PHOS manual (version from Jan 2007) 
	  */ 
	  if (icsp<16) { csp_addr = trailer_offset[RCUFWVersion] + icsp; }
	  else { csp_addr = trailer_offset[RCUFWVersion] - 1 - (icsp%16); }
	  if (icsp >= 24) csp_addr += 0x20;

	  // what does the desired voltage (in V) correspond to in DAC?
	  int iv_dac = (int)( (biasVoltage[rcu][branch][card][icsp] - hvmin)/prop + 0.5); // round-off
	  if (iv_dac > max_dac_value) iv_dac = max_dac_value;
	  sprintf(iv_dac_setting,"%06X", voltage_value_start[RCUFWVersion] + iv_dac);

	  // set up instructions that should be written
	  word = write_header[RCUFWVersion] | (branch << 16) | (icard << 12) | (csp_addr);

	  // write a long comment with all info for this CSP
	  sprintf(comment, "# RCU %d, Branch %s, FEC %d, CSP %02d - Tower Col %02d, Row %02d ", 
		  rcu, branch_str[branch], icard, icsp,
		  towerCol[rcu][branch][card][icsp],
		  towerRow[rcu][branch][card][icsp]
		  );  
	
	  fprintf(fout_setbias_card[rcu][branch][card], "w 0x%04X 0x%6X   %s\n",
		  rcu_addr_card, word, comment);
	  rcu_addr_card++;

	  fprintf(fout_setbias_card[rcu][branch][card], "w 0x%04X 0x%s   # Set Voltage: %4.1f V, DAC %d (hex: %03X)\n", 
		  rcu_addr_card, iv_dac_setting, 
		  biasVoltage[rcu][branch][card][icsp], 
		  iv_dac, iv_dac
		  );
	  rcu_addr_card++;

	  // slighly modified comment for read command - include voltage info
	  sprintf(comment, "# RCU %d, Branch %s, FEC %d, CSP %02d - Tower Col %02d, Row %02d : %4.1f V, DAC %d (hex: %03X)", 
		  rcu, branch_str[branch], icard, icsp,
		  towerCol[rcu][branch][card][icsp],
		  towerRow[rcu][branch][card][icsp],
		  biasVoltage[rcu][branch][card][icsp], 
		  iv_dac, iv_dac
		  );  

	  word = read_header[RCUFWVersion] | (branch << 16) | (icard << 12) | (csp_addr);
	  fprintf(fout_readbias_card[rcu][branch][card], "w 0x%04X 0x%06X  %s\n", rcu_addr_read, word, comment);
	  rcu_addr_read++;
	} // csp loop
	
	// after CSP per card; send update command
	word = write_header[RCUFWVersion] | (branch << 16) | (icard << 12) | update_voltage_register;
	fprintf(fout_setbias_card[rcu][branch][card],"w 0x%04X 0x%06X   # Update Voltages\n", 
		rcu_addr_card, word); 
	rcu_addr_card++;
	fprintf(fout_setbias_card[rcu][branch][card],"w 0x%04X 0x%06X   \n", 
		rcu_addr_card, voltage_value_start[RCUFWVersion]);
	rcu_addr_card++;

	// also put ending for the individual card files:
	fprintf(fout_setbias_card[rcu][branch][card],"w 0x%04X 0x%06X   # End of the instruction memory\n", 
		rcu_addr_card, rcu_endmem[RCUFWVersion]);
	rcu_addr_card++;
      
	fprintf(fout_setbias_card[rcu][branch][card],"wait 100 us\n");
	fprintf(fout_setbias_card[rcu][branch][card],"w 0x%X 0x0           # execute and update registers\n",
		rcu_exeseq[RCUFWVersion]);
	if (RCUFWVersion == 0) { // specialty for old FW version
	  fprintf(fout_setbias_card[rcu][branch][card],"wait 100 us\n");
	  fprintf(fout_setbias_card[rcu][branch][card],"r 0x7800            # error checking\n");
	  fprintf(fout_setbias_card[rcu][branch][card],"w 0x6c01 0x 0       # clear registers\n");
	}

	// in case we want to check what was written
	if (readBack) {
	  fprintf(fout_setbias_card[rcu][branch][card],"wait 100 us\n");
	  fprintf(fout_setbias_card[rcu][branch][card],"b %s      # read-back the values also\n", cfile);
	  fprintf(fout_setbias_card[rcu][branch][card],"wait 100 us\n");
	}

	// close down output files (set)
	fclose(fout_setbias_card[rcu][branch][card]);
	
	// readbias ending
	fprintf(fout_readbias_card[rcu][branch][card],"w 0x%04X 0x%06X \n", 
		rcu_addr_read, rcu_endmem[RCUFWVersion]);
	rcu_addr_read++;

	fprintf(fout_readbias_card[rcu][branch][card],"wait 1 us\n");
	fprintf(fout_readbias_card[rcu][branch][card],"w 0x%X 0x0           # execute and update registers\n",
		rcu_exeseq[RCUFWVersion]);

	fprintf(fout_readbias_card[rcu][branch][card],"wait 1 us\n");
	int nRead = NCSP*2; // new FW reads back i/o 
	if (RCUFWVersion == 0) { // specialty for old FW version
	  fprintf(fout_readbias_card[rcu][branch][card],"r 0x7800            # error checking\n");
	  fprintf(fout_readbias_card[rcu][branch][card],"wait 1 us\n");
	  nRead = NCSP;
	}
	fprintf(fout_readbias_card[rcu][branch][card],"r 0x%04X %d( \n", rcu_result_reg[RCUFWVersion], nRead);
	if (RCUFWVersion == 0) { // specialty for old FW version
	  fprintf(fout_readbias_card[rcu][branch][card],"wait 1 us\n");
	  fprintf(fout_readbias_card[rcu][branch][card],"r 0x7800            # error checking\n");
	  fprintf(fout_readbias_card[rcu][branch][card],"w 0x6c01 0x 0       # clear registers\n");
	}

	// close down output files (read)
	fclose(fout_readbias_card[rcu][branch][card]);

      } // card=FEC
    } // branch
  } // rcu

}
