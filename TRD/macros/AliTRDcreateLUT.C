/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that translates an arbitrary function and sacle factor into a look-up table (LUT)         |
|                                                                                                    |
| Author: 10/2016, Benjamin Brudnyj, Jochen Klein                                                    |
\*--------------------------------------------------------------------------------------------------*/

Int_t GetPID(Int_t q)
{
  const Float_t a1 = 0.100251; //0.0889764;
  const Float_t a2 = 0.0193284; //0.0193284;
  const Float_t b1 = -108.571; //-101.892;
  const Float_t b2 = 113.847; //1.06657899999999998e+02;
  const Int_t thr  = 2750; //2988;

  TF1 elLUT("elLUT","pol1",0,thr);
  TF1 nuclLUT("nuclLUT","pol1",thr,20005);
  elLUT.SetParameters(b1, a1);
  nuclLUT.SetParameters(b2, a2);

  Int_t result;

  if (q < thr) {
    result = ((Int_t) (a1 * q + b1));
    if (TMath::Abs(result - elLUT.Eval(q)) >= 1)
      printf("error: %i <-> %i\n", result, elLUT.Eval(q));
  }
  else {
    result = ((Int_t) (a2 * q + b2));
    if (TMath::Abs(result - nuclLUT.Eval(q)) >= 1)
      printf("error: %i <-> %i\n", result, nuclLUT.Eval(q));
  }

  if (result < 0)
    result = 0;
  else if (result > 255)
    result = 255;

  return result;
}

void AliTRDcreateLUT()
{
  const Int_t Bins            = 2008;

  const char *LUT_file_name   = "lhc11dv3en";
  const char *file_type       = "tcs";

  ULong64_t scaleq0           = 429496730;
  ULong64_t scaleq1           = 0;

  Int_t addrPID[Bins]         = {0};

  for(Int_t i = 0; i < Bins; i++) {
    Int_t q = ((1ull << 32) * (i + .5)) / scaleq0;

    if (i != (((q * scaleq0)>> 16)>>16))
      printf("error: bin %i -> %i \n", i, (((q * scaleq0)>> 16)>>16));

    addrPID[i/4] |= (GetPID(q) & 0xff) << ((i % 4) * 8);
  }

  ofstream LUTfile;
  LUTfile.open(Form("%s.%s",LUT_file_name,file_type));

  LUTfile << setw(8) << "tracklet" << setw(8) << "scaleq0" << setw(14) << Form("%u;\n",scaleq0);
  LUTfile << setw(8) << "tracklet" << setw(8) << "scaleq1" << setw(14) << Form("%u;\n",scaleq1);
  LUTfile << setw(8) << "write   " << setw(8) << "127, "   << setw(14) << "0xc029,  " << setw(13) << "2008;\n";
  LUTfile << setw(8) << "write   " << setw(8) << "127, "   << setw(14) << "0xc02b,  " << setw(13) << "2008;\n";

  for(Int_t i = 0; i < Bins/4; i++)
    LUTfile << setw(8) << "write   " << setw(8) << "127, "
            << setw(14) << Form("0x%x,  ",0xc100 + i)
            << setw(13) << Form("%i;\n", addrPID[i]);

  LUTfile.close();
}
