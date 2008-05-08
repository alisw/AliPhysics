/*
Currently the same average correction (1.7% change per degree C) is assumed
for all APDs, i.e. all columns and rows.
Below is an example for how the correction can be found for a run or within 
a single runnumber.
*/

//____________________________________________________________________
void CalibTestBeam(const int runno=1838) 
{ 

  gSystem->Load("AliEMCALCalibTestBeam_cxx");

  AliEMCALCalibTestBeam *calibtb;
  calibtb = new AliEMCALCalibTestBeam(runno);

  double hours = calibtb->GetLengthOfRunInHours();

  // This commented out printf statement gave the log that is included at the bottom of this macro..
  /*
  printf("Run: %4i NEvents:%6i RunDuration(hours): %5.3f  # Run /", 
	 runno, calibtb->GetNEvents(), hours);
  printf(" Temp. data #  RangeInHours: %5.3f RangeInDegrees: %4.2f Nmeas: %4i \n",
	 calibtb->GetRangeOfTempMeasureInHours(),
	 calibtb->GetRangeOfTempMeasureInDegrees(),
	 calibtb->GetNTempVal() );
  */

  calibtb->Print(); // dump the basic info

  // example of how correction changes over time
  for (int i=0; i<hours; i++) {
    int secSinceStart = i*3600;
    double temperature = calibtb->GetTemperature(secSinceStart);
    double corr = calibtb->GetCorrection(temperature);
    cout << " hour " << i
	 << " temperature " << temperature
	 << " correction " << corr
	 << endl;
    // the ADC/amplitude for a single tower should then be modified a la
    // correctedAmplitude = corr * amplitude;
    // Since the APD yield increases as the temperature decreases, the correction
    //goes in the reverse direction. 
    //[correction value decreases as temperature decreases]     
  }

}

/* 
// Log from printf commands above; gives a summary for the testbeam runs
// - can illustrate which runs may be of interest for a calibrated analysis
Run:  147 NEvents: 10635 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  148 NEvents:  5150 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  156 NEvents:   103 RunDuration(hours): 0.007  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  176 NEvents:  1676 RunDuration(hours): 0.252  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  177 NEvents: 21337 RunDuration(hours): 0.233  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  178 NEvents: 51593 RunDuration(hours): 0.258  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  179 NEvents: 53204 RunDuration(hours): 0.237  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  180 NEvents: 51649 RunDuration(hours): 0.222  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  181 NEvents: 23073 RunDuration(hours): 0.224  # Run / Temp. data #  RangeInHours: 0.142 RangeInDegrees: 0.04 Nmeas:   34 
Run:  182 NEvents: 63759 RunDuration(hours): 0.336  # Run / Temp. data #  RangeInHours: 0.325 RangeInDegrees: 0.08 Nmeas:   87 
Run:  185 NEvents: 10113 RunDuration(hours): 0.185  # Run / Temp. data #  RangeInHours: 0.179 RangeInDegrees: 0.04 Nmeas:   65 
Run:  186 NEvents: 40287 RunDuration(hours): 0.213  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.04 Nmeas:   76 
Run:  187 NEvents: 46700 RunDuration(hours): 0.200  # Run / Temp. data #  RangeInHours: 0.195 RangeInDegrees: 0.03 Nmeas:   35 
Run:  188 NEvents: 52807 RunDuration(hours): 0.204  # Run / Temp. data #  RangeInHours: 0.202 RangeInDegrees: 0.04 Nmeas:   41 
Run:  189 NEvents: 51946 RunDuration(hours): 0.200  # Run / Temp. data #  RangeInHours: 0.199 RangeInDegrees: 0.04 Nmeas:   65 
Run:  190 NEvents: 46997 RunDuration(hours): 0.186  # Run / Temp. data #  RangeInHours: 0.166 RangeInDegrees: 0.03 Nmeas:   38 
Run:  191 NEvents: 30849 RunDuration(hours): 0.241  # Run / Temp. data #  RangeInHours: 0.228 RangeInDegrees: 0.03 Nmeas:   41 
Run:  192 NEvents: 19844 RunDuration(hours): 0.116  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.03 Nmeas:   39 
Run:  193 NEvents: 51412 RunDuration(hours): 0.219  # Run / Temp. data #  RangeInHours: 0.172 RangeInDegrees: 0.01 Nmeas:   28 
Run:  194 NEvents: 37973 RunDuration(hours): 0.215  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.03 Nmeas:   47 
Run:  195 NEvents: 51491 RunDuration(hours): 0.229  # Run / Temp. data #  RangeInHours: 0.217 RangeInDegrees: 0.04 Nmeas:   48 
Run:  196 NEvents: 50927 RunDuration(hours): 0.217  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.04 Nmeas:   73 
Run:  197 NEvents: 50242 RunDuration(hours): 0.202  # Run / Temp. data #  RangeInHours: 0.200 RangeInDegrees: 0.03 Nmeas:   72 
Run:  198 NEvents: 50036 RunDuration(hours): 0.205  # Run / Temp. data #  RangeInHours: 0.202 RangeInDegrees: 0.04 Nmeas:   74 
Run:  199 NEvents: 52872 RunDuration(hours): 0.230  # Run / Temp. data #  RangeInHours: 0.228 RangeInDegrees: 0.03 Nmeas:   86 
Run:  200 NEvents: 33405 RunDuration(hours): 0.203  # Run / Temp. data #  RangeInHours: 0.202 RangeInDegrees: 0.03 Nmeas:   57 
Run:  201 NEvents: 51321 RunDuration(hours): 0.197  # Run / Temp. data #  RangeInHours: 0.197 RangeInDegrees: 0.03 Nmeas:   62 
Run:  202 NEvents: 51287 RunDuration(hours): 0.250  # Run / Temp. data #  RangeInHours: 0.246 RangeInDegrees: 0.03 Nmeas:   71 
Run:  203 NEvents: 53594 RunDuration(hours): 0.219  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.03 Nmeas:   67 
Run:  204 NEvents:  4180 RunDuration(hours): 0.205  # Run / Temp. data #  RangeInHours: 0.204 RangeInDegrees: 0.03 Nmeas:   51 
Run:  205 NEvents: 51367 RunDuration(hours): 0.197  # Run / Temp. data #  RangeInHours: 0.195 RangeInDegrees: 0.04 Nmeas:   38 
Run:  206 NEvents: 52201 RunDuration(hours): 0.333  # Run / Temp. data #  RangeInHours: 0.316 RangeInDegrees: 0.03 Nmeas:   80 
Run:  208 NEvents: 78719 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.206 RangeInDegrees: 0.03 Nmeas:   45 
Run:  209 NEvents: 50344 RunDuration(hours): 0.146  # Run / Temp. data #  RangeInHours: 0.142 RangeInDegrees: 0.03 Nmeas:   24 
Run:  210 NEvents: 49665 RunDuration(hours): 0.139  # Run / Temp. data #  RangeInHours: 0.108 RangeInDegrees: 0.03 Nmeas:   18 
Run:  211 NEvents: 58654 RunDuration(hours): 0.169  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.03 Nmeas:   19 
Run:  212 NEvents: 52056 RunDuration(hours): 0.140  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.03 Nmeas:   13 
Run:  213 NEvents: 58463 RunDuration(hours): 0.167  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.01 Nmeas:   15 
Run:  214 NEvents: 50331 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.04 Nmeas:   23 
Run:  215 NEvents: 51810 RunDuration(hours): 0.136  # Run / Temp. data #  RangeInHours: 0.134 RangeInDegrees: 0.03 Nmeas:   23 
Run:  216 NEvents: 50795 RunDuration(hours): 0.143  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   17 
Run:  217 NEvents: 53316 RunDuration(hours): 0.138  # Run / Temp. data #  RangeInHours: 0.074 RangeInDegrees: 0.03 Nmeas:   24 
Run:  218 NEvents: 50791 RunDuration(hours): 0.130  # Run / Temp. data #  RangeInHours: 0.123 RangeInDegrees: 0.03 Nmeas:   31 
Run:  219 NEvents: 10105 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   34 
Run:  220 NEvents: 54045 RunDuration(hours): 0.139  # Run / Temp. data #  RangeInHours: 0.134 RangeInDegrees: 0.04 Nmeas:   43 
Run:  221 NEvents: 52132 RunDuration(hours): 0.246  # Run / Temp. data #  RangeInHours: 0.229 RangeInDegrees: 0.04 Nmeas:   68 
Run:  222 NEvents: 58015 RunDuration(hours): 0.158  # Run / Temp. data #  RangeInHours: 0.153 RangeInDegrees: 0.03 Nmeas:   46 
Run:  223 NEvents: 58835 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   48 
Run:  224 NEvents: 52261 RunDuration(hours): 0.136  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   38 
Run:  226 NEvents:  2546 RunDuration(hours): 0.140  # Run / Temp. data #  RangeInHours: 0.137 RangeInDegrees: 0.04 Nmeas:   46 
Run:  227 NEvents: 51614 RunDuration(hours): 0.141  # Run / Temp. data #  RangeInHours: 0.134 RangeInDegrees: 0.03 Nmeas:   38 
Run:  228 NEvents: 60129 RunDuration(hours): 0.152  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.03 Nmeas:   42 
Run:  229 NEvents: 23132 RunDuration(hours): 0.133  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.03 Nmeas:   28 
Run:  230 NEvents: 54656 RunDuration(hours): 0.162  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.03 Nmeas:   20 
Run:  231 NEvents: 53620 RunDuration(hours): 0.154  # Run / Temp. data #  RangeInHours: 0.101 RangeInDegrees: 0.03 Nmeas:   19 
Run:  232 NEvents: 45971 RunDuration(hours): 0.144  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.03 Nmeas:   21 
Run:  233 NEvents: 52110 RunDuration(hours): 0.136  # Run / Temp. data #  RangeInHours: 0.119 RangeInDegrees: 0.03 Nmeas:    8 
Run:  234 NEvents: 51988 RunDuration(hours): 0.146  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.03 Nmeas:   17 
Run:  235 NEvents: 55136 RunDuration(hours): 0.148  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.03 Nmeas:   11 
Run:  236 NEvents: 50203 RunDuration(hours): 0.160  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.03 Nmeas:   11 
Run:  237 NEvents: 50342 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.04 Nmeas:   31 
Run:  238 NEvents: 51686 RunDuration(hours): 0.130  # Run / Temp. data #  RangeInHours: 0.088 RangeInDegrees: 0.03 Nmeas:   25 
Run:  239 NEvents: 50941 RunDuration(hours): 0.128  # Run / Temp. data #  RangeInHours: 0.121 RangeInDegrees: 0.05 Nmeas:   27 
Run:  240 NEvents: 52142 RunDuration(hours): 0.142  # Run / Temp. data #  RangeInHours: 0.141 RangeInDegrees: 0.03 Nmeas:   31 
Run:  242 NEvents: 53125 RunDuration(hours): 0.131  # Run / Temp. data #  RangeInHours: 0.125 RangeInDegrees: 0.03 Nmeas:   28 
Run:  244 NEvents: 54734 RunDuration(hours): 0.137  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   41 
Run:  245 NEvents: 58234 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.143 RangeInDegrees: 0.04 Nmeas:   42 
Run:  246 NEvents: 57829 RunDuration(hours): 0.133  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   54 
Run:  252 NEvents: 53474 RunDuration(hours): 0.135  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   45 
Run:  253 NEvents: 55012 RunDuration(hours): 0.136  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.02 Nmeas:   41 
Run:  254 NEvents: 51290 RunDuration(hours): 0.657  # Run / Temp. data #  RangeInHours: 0.652 RangeInDegrees: 0.04 Nmeas:  206 
Run:  255 NEvents: 52488 RunDuration(hours): 0.619  # Run / Temp. data #  RangeInHours: 0.608 RangeInDegrees: 0.06 Nmeas:  140 
Run:  256 NEvents: 50741 RunDuration(hours): 0.676  # Run / Temp. data #  RangeInHours: 0.643 RangeInDegrees: 0.11 Nmeas:   80 
Run:  257 NEvents: 65102 RunDuration(hours): 0.585  # Run / Temp. data #  RangeInHours: 0.573 RangeInDegrees: 0.10 Nmeas:  131 
Run:  258 NEvents: 54041 RunDuration(hours): 0.502  # Run / Temp. data #  RangeInHours: 0.499 RangeInDegrees: 0.09 Nmeas:  172 
Run:  259 NEvents: 52277 RunDuration(hours): 0.225  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.05 Nmeas:   55 
Run:  260 NEvents: 55272 RunDuration(hours): 0.196  # Run / Temp. data #  RangeInHours: 0.179 RangeInDegrees: 0.04 Nmeas:   33 
Run:  261 NEvents: 52085 RunDuration(hours): 0.194  # Run / Temp. data #  RangeInHours: 0.182 RangeInDegrees: 0.04 Nmeas:   29 
Run:  262 NEvents: 53085 RunDuration(hours): 0.213  # Run / Temp. data #  RangeInHours: 0.191 RangeInDegrees: 0.04 Nmeas:   30 
Run:  263 NEvents:  8627 RunDuration(hours): 0.306  # Run / Temp. data #  RangeInHours: 0.305 RangeInDegrees: 0.06 Nmeas:   43 
Run:  264 NEvents: 52565 RunDuration(hours): 0.199  # Run / Temp. data #  RangeInHours: 0.193 RangeInDegrees: 0.04 Nmeas:   40 
Run:  265 NEvents: 52178 RunDuration(hours): 0.196  # Run / Temp. data #  RangeInHours: 0.189 RangeInDegrees: 0.06 Nmeas:   51 
Run:  266 NEvents: 38866 RunDuration(hours): 0.211  # Run / Temp. data #  RangeInHours: 0.204 RangeInDegrees: 0.05 Nmeas:   75 
Run:  267 NEvents: 19195 RunDuration(hours): 0.196  # Run / Temp. data #  RangeInHours: 0.195 RangeInDegrees: 0.05 Nmeas:   74 
Run:  268 NEvents: 41034 RunDuration(hours): 0.203  # Run / Temp. data #  RangeInHours: 0.184 RangeInDegrees: 0.05 Nmeas:   59 
Run:  269 NEvents: 34584 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.213 RangeInDegrees: 0.05 Nmeas:   40 
Run:  270 NEvents: 47310 RunDuration(hours): 0.212  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.04 Nmeas:   38 
Run:  271 NEvents: 52082 RunDuration(hours): 0.209  # Run / Temp. data #  RangeInHours: 0.197 RangeInDegrees: 0.06 Nmeas:   24 
Run:  272 NEvents: 14027 RunDuration(hours): 0.107  # Run / Temp. data #  RangeInHours: 0.036 RangeInDegrees: 0.01 Nmeas:   13 
Run:  273 NEvents:  4151 RunDuration(hours): 0.112  # Run / Temp. data #  RangeInHours: 0.047 RangeInDegrees: 0.01 Nmeas:   11 
Run:  274 NEvents: 20897 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.119 RangeInDegrees: 0.04 Nmeas:   19 
Run:  275 NEvents: 27215 RunDuration(hours): 0.121  # Run / Temp. data #  RangeInHours: 0.088 RangeInDegrees: 0.03 Nmeas:   25 
Run:  276 NEvents: 61525 RunDuration(hours): 0.226  # Run / Temp. data #  RangeInHours: 0.188 RangeInDegrees: 0.04 Nmeas:   38 
Run:  277 NEvents: 52810 RunDuration(hours): 0.247  # Run / Temp. data #  RangeInHours: 0.222 RangeInDegrees: 0.04 Nmeas:   36 
Run:  278 NEvents: 43383 RunDuration(hours): 0.183  # Run / Temp. data #  RangeInHours: 0.181 RangeInDegrees: 0.01 Nmeas:   16 
Run:  282 NEvents: 51926 RunDuration(hours): 0.149  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   20 
Run:  283 NEvents:   101 RunDuration(hours): 0.031  # Run / Temp. data #  RangeInHours: 0.020 RangeInDegrees: 0.01 Nmeas:   10 
Run:  284 NEvents:  1116 RunDuration(hours): 0.148  # Run / Temp. data #  RangeInHours: 0.137 RangeInDegrees: 0.01 Nmeas:   24 
Run:  285 NEvents: 50097 RunDuration(hours): 0.176  # Run / Temp. data #  RangeInHours: 0.164 RangeInDegrees: 0.01 Nmeas:   10 
Run:  286 NEvents: 51549 RunDuration(hours): 0.174  # Run / Temp. data #  RangeInHours: 0.152 RangeInDegrees: 0.03 Nmeas:   40 
Run:  287 NEvents: 25814 RunDuration(hours): 0.143  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   32 
Run:  288 NEvents: 18810 RunDuration(hours): 0.088  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.03 Nmeas:   30 
Run:  289 NEvents: 49670 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   47 
Run:  290 NEvents: 52972 RunDuration(hours): 0.143  # Run / Temp. data #  RangeInHours: 0.137 RangeInDegrees: 0.03 Nmeas:   40 
Run:  291 NEvents: 51061 RunDuration(hours): 0.131  # Run / Temp. data #  RangeInHours: 0.126 RangeInDegrees: 0.03 Nmeas:   45 
Run:  292 NEvents: 50850 RunDuration(hours): 0.151  # Run / Temp. data #  RangeInHours: 0.137 RangeInDegrees: 0.03 Nmeas:   45 
Run:  293 NEvents: 52556 RunDuration(hours): 0.154  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.03 Nmeas:   25 
Run:  294 NEvents: 48078 RunDuration(hours): 0.162  # Run / Temp. data #  RangeInHours: 0.155 RangeInDegrees: 0.03 Nmeas:   45 
Run:  295 NEvents: 33861 RunDuration(hours): 0.141  # Run / Temp. data #  RangeInHours: 0.138 RangeInDegrees: 0.03 Nmeas:   13 
Run:  296 NEvents:   100 RunDuration(hours): 0.007  # Run / Temp. data #  RangeInHours: 0.005 RangeInDegrees: 0.03 Nmeas:    4 
Run:  306 NEvents: 48478 RunDuration(hours): 0.101  # Run / Temp. data #  RangeInHours: 0.101 RangeInDegrees: 0.03 Nmeas:   28 
Run:  307 NEvents: 36646 RunDuration(hours): 0.092  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.03 Nmeas:   20 
Run:  308 NEvents: 51892 RunDuration(hours): 0.093  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.03 Nmeas:   32 
Run:  309 NEvents: 53709 RunDuration(hours): 0.085  # Run / Temp. data #  RangeInHours: 0.072 RangeInDegrees: 0.03 Nmeas:   27 
Run:  310 NEvents: 51979 RunDuration(hours): 0.103  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.03 Nmeas:   33 
Run:  311 NEvents: 54786 RunDuration(hours): 0.096  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.03 Nmeas:   20 
Run:  312 NEvents: 54612 RunDuration(hours): 0.100  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.03 Nmeas:   19 
Run:  313 NEvents: 30663 RunDuration(hours): 0.106  # Run / Temp. data #  RangeInHours: 0.099 RangeInDegrees: 0.03 Nmeas:   19 
Run:  314 NEvents: 42169 RunDuration(hours): 0.086  # Run / Temp. data #  RangeInHours: 0.081 RangeInDegrees: 0.04 Nmeas:   20 
Run:  315 NEvents: 50477 RunDuration(hours): 0.106  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.04 Nmeas:   19 
Run:  316 NEvents: 20385 RunDuration(hours): 0.058  # Run / Temp. data #  RangeInHours: 0.023 RangeInDegrees: 0.01 Nmeas:    6 
Run:  317 NEvents: 56260 RunDuration(hours): 0.102  # Run / Temp. data #  RangeInHours: 0.077 RangeInDegrees: 0.03 Nmeas:   22 
Run:  318 NEvents: 51883 RunDuration(hours): 0.109  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.04 Nmeas:   19 
Run:  319 NEvents: 51847 RunDuration(hours): 0.094  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.04 Nmeas:   17 
Run:  321 NEvents:  4825 RunDuration(hours): 0.073  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.03 Nmeas:   16 
Run:  322 NEvents: 51768 RunDuration(hours): 0.089  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.03 Nmeas:   10 
Run:  323 NEvents: 53528 RunDuration(hours): 0.091  # Run / Temp. data #  RangeInHours: 0.058 RangeInDegrees: 0.03 Nmeas:    9 
Run:  324 NEvents: 52830 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.03 Nmeas:   23 
Run:  325 NEvents: 59117 RunDuration(hours): 0.604  # Run / Temp. data #  RangeInHours: 0.598 RangeInDegrees: 0.08 Nmeas:  112 
Run:  326 NEvents: 50611 RunDuration(hours): 0.136  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.01 Nmeas:   28 
Run:  327 NEvents:  8846 RunDuration(hours): 0.237  # Run / Temp. data #  RangeInHours: 0.222 RangeInDegrees: 0.04 Nmeas:   33 
Run:  328 NEvents:  1000 RunDuration(hours): 0.020  # Run / Temp. data #  RangeInHours: 0.004 RangeInDegrees: 0.01 Nmeas:    3 
Run:  329 NEvents:  6267 RunDuration(hours): 0.272  # Run / Temp. data #  RangeInHours: 0.253 RangeInDegrees: 0.04 Nmeas:   38 
Run:  330 NEvents:  1000 RunDuration(hours): 0.019  # Run / Temp. data #  RangeInHours: 0.011 RangeInDegrees: 0.03 Nmeas:    6 
Run:  331 NEvents: 74222 RunDuration(hours): 0.132  # Run / Temp. data #  RangeInHours: 0.128 RangeInDegrees: 0.03 Nmeas:   50 
Run:  332 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  333 NEvents: 50440 RunDuration(hours): 0.092  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.03 Nmeas:   47 
Run:  334 NEvents: 16291 RunDuration(hours): 0.045  # Run / Temp. data #  RangeInHours: 0.041 RangeInDegrees: 0.04 Nmeas:   18 
Run:  335 NEvents: 40947 RunDuration(hours): 0.097  # Run / Temp. data #  RangeInHours: 0.094 RangeInDegrees: 0.03 Nmeas:   35 
Run:  336 NEvents: 51571 RunDuration(hours): 0.088  # Run / Temp. data #  RangeInHours: 0.078 RangeInDegrees: 0.03 Nmeas:   27 
Run:  337 NEvents:   101 RunDuration(hours): 0.004  # Run / Temp. data #  RangeInHours: 0.004 RangeInDegrees: 0.01 Nmeas:    2 
Run:  338 NEvents: 56037 RunDuration(hours): 0.095  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.01 Nmeas:   22 
Run:  339 NEvents:  2213 RunDuration(hours): 0.026  # Run / Temp. data #  RangeInHours: 0.022 RangeInDegrees: 0.03 Nmeas:   12 
Run:  340 NEvents: 53673 RunDuration(hours): 0.083  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.01 Nmeas:   12 
Run:  341 NEvents: 56073 RunDuration(hours): 0.097  # Run / Temp. data #  RangeInHours: 0.052 RangeInDegrees: 0.01 Nmeas:   15 
Run:  342 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  343 NEvents:  1932 RunDuration(hours): 0.042  # Run / Temp. data #  RangeInHours: 0.034 RangeInDegrees: 0.01 Nmeas:   15 
Run:  344 NEvents: 18875 RunDuration(hours): 0.047  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.01 Nmeas:    7 
Run:  345 NEvents:  3198 RunDuration(hours): 0.021  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  346 NEvents: 18164 RunDuration(hours): 0.048  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.01 Nmeas:    9 
Run:  347 NEvents: 54032 RunDuration(hours): 0.102  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.03 Nmeas:    7 
Run:  348 NEvents: 51956 RunDuration(hours): 0.104  # Run / Temp. data #  RangeInHours: 0.036 RangeInDegrees: 0.03 Nmeas:   11 
Run:  349 NEvents: 53885 RunDuration(hours): 0.133  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.01 Nmeas:   17 
Run:  351 NEvents:  4781 RunDuration(hours): 0.085  # Run / Temp. data #  RangeInHours: 0.051 RangeInDegrees: 0.03 Nmeas:   11 
Run:  352 NEvents: 52620 RunDuration(hours): 0.118  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.03 Nmeas:   22 
Run:  353 NEvents: 74054 RunDuration(hours): 0.186  # Run / Temp. data #  RangeInHours: 0.164 RangeInDegrees: 0.03 Nmeas:   37 
Run:  354 NEvents: 51289 RunDuration(hours): 0.129  # Run / Temp. data #  RangeInHours: 0.052 RangeInDegrees: 0.03 Nmeas:   21 
Run:  355 NEvents: 25884 RunDuration(hours): 0.120  # Run / Temp. data #  RangeInHours: 0.110 RangeInDegrees: 0.03 Nmeas:   31 
Run:  356 NEvents: 54892 RunDuration(hours): 0.219  # Run / Temp. data #  RangeInHours: 0.208 RangeInDegrees: 0.04 Nmeas:   61 
Run:  357 NEvents: 53013 RunDuration(hours): 0.201  # Run / Temp. data #  RangeInHours: 0.186 RangeInDegrees: 0.04 Nmeas:   49 
Run:  358 NEvents: 50385 RunDuration(hours): 0.239  # Run / Temp. data #  RangeInHours: 0.238 RangeInDegrees: 0.05 Nmeas:   89 
Run:  359 NEvents:   101 RunDuration(hours): 0.007  # Run / Temp. data #  RangeInHours: 0.007 RangeInDegrees: 0.01 Nmeas:    5 
Run:  360 NEvents: 54098 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.101 RangeInDegrees: 0.03 Nmeas:   44 
Run:  361 NEvents:  2501 RunDuration(hours): 0.055  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.03 Nmeas:   21 
Run:  362 NEvents: 51452 RunDuration(hours): 0.105  # Run / Temp. data #  RangeInHours: 0.094 RangeInDegrees: 0.04 Nmeas:   29 
Run:  363 NEvents: 10975 RunDuration(hours): 0.040  # Run / Temp. data #  RangeInHours: 0.034 RangeInDegrees: 0.03 Nmeas:   16 
Run:  364 NEvents: 31438 RunDuration(hours): 0.087  # Run / Temp. data #  RangeInHours: 0.071 RangeInDegrees: 0.03 Nmeas:   20 
Run:  365 NEvents: 31919 RunDuration(hours): 0.076  # Run / Temp. data #  RangeInHours: 0.054 RangeInDegrees: 0.03 Nmeas:   21 
Run:  366 NEvents: 57265 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.101 RangeInDegrees: 0.03 Nmeas:   27 
Run:  367 NEvents: 15179 RunDuration(hours): 0.082  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.03 Nmeas:   20 
Run:  368 NEvents:  3898 RunDuration(hours): 0.026  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.03 Nmeas:    9 
Run:  370 NEvents: 37343 RunDuration(hours): 0.095  # Run / Temp. data #  RangeInHours: 0.077 RangeInDegrees: 0.03 Nmeas:   11 
Run:  371 NEvents: 40086 RunDuration(hours): 0.098  # Run / Temp. data #  RangeInHours: 0.088 RangeInDegrees: 0.03 Nmeas:   12 
Run:  372 NEvents: 13832 RunDuration(hours): 0.037  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  373 NEvents:  4876 RunDuration(hours): 0.025  # Run / Temp. data #  RangeInHours: 0.022 RangeInDegrees: 0.03 Nmeas:    6 
Run:  374 NEvents: 20465 RunDuration(hours): 0.054  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.01 Nmeas:    6 
Run:  375 NEvents: 11132 RunDuration(hours): 0.031  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.03 Nmeas:    6 
Run:  376 NEvents:  7957 RunDuration(hours): 0.028  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  377 NEvents: 10649 RunDuration(hours): 0.033  # Run / Temp. data #  RangeInHours: 0.025 RangeInDegrees: 0.01 Nmeas:   12 
Run:  378 NEvents: 54986 RunDuration(hours): 0.139  # Run / Temp. data #  RangeInHours: 0.138 RangeInDegrees: 0.04 Nmeas:   27 
Run:  379 NEvents: 52148 RunDuration(hours): 0.130  # Run / Temp. data #  RangeInHours: 0.119 RangeInDegrees: 0.03 Nmeas:   21 
Run:  380 NEvents: 49235 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.128 RangeInDegrees: 0.03 Nmeas:   43 
Run:  382 NEvents: 51852 RunDuration(hours): 0.150  # Run / Temp. data #  RangeInHours: 0.148 RangeInDegrees: 0.03 Nmeas:   41 
Run:  383 NEvents:  8026 RunDuration(hours): 0.039  # Run / Temp. data #  RangeInHours: 0.036 RangeInDegrees: 0.03 Nmeas:   17 
Run:  384 NEvents: 51270 RunDuration(hours): 0.151  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.04 Nmeas:   52 
Run:  385 NEvents: 51207 RunDuration(hours): 0.150  # Run / Temp. data #  RangeInHours: 0.148 RangeInDegrees: 0.04 Nmeas:   51 
Run:  386 NEvents:   300 RunDuration(hours): 0.017  # Run / Temp. data #  RangeInHours: 0.013 RangeInDegrees: 0.01 Nmeas:    6 
Run:  387 NEvents: 51714 RunDuration(hours): 0.311  # Run / Temp. data #  RangeInHours: 0.309 RangeInDegrees: 0.05 Nmeas:   92 
Run:  389 NEvents: 50567 RunDuration(hours): 0.307  # Run / Temp. data #  RangeInHours: 0.305 RangeInDegrees: 0.04 Nmeas:   82 
Run:  390 NEvents: 47199 RunDuration(hours): 0.545  # Run / Temp. data #  RangeInHours: 0.542 RangeInDegrees: 0.05 Nmeas:  115 
Run:  391 NEvents: 51072 RunDuration(hours): 0.535  # Run / Temp. data #  RangeInHours: 0.522 RangeInDegrees: 0.03 Nmeas:  119 
Run:  392 NEvents: 53732 RunDuration(hours): 0.174  # Run / Temp. data #  RangeInHours: 0.171 RangeInDegrees: 0.03 Nmeas:   63 
Run:  394 NEvents:  3941 RunDuration(hours): 0.070  # Run / Temp. data #  RangeInHours: 0.067 RangeInDegrees: 0.01 Nmeas:   23 
Run:  395 NEvents: 27986 RunDuration(hours): 0.092  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.03 Nmeas:   33 
Run:  396 NEvents: 50391 RunDuration(hours): 0.143  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   46 
Run:  397 NEvents:  1720 RunDuration(hours): 0.056  # Run / Temp. data #  RangeInHours: 0.052 RangeInDegrees: 0.03 Nmeas:   19 
Run:  398 NEvents: 14245 RunDuration(hours): 0.047  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.03 Nmeas:   21 
Run:  399 NEvents: 20089 RunDuration(hours): 0.055  # Run / Temp. data #  RangeInHours: 0.052 RangeInDegrees: 0.03 Nmeas:    4 
Run:  400 NEvents: 51821 RunDuration(hours): 0.083  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.01 Nmeas:    8 
Run:  401 NEvents: 36029 RunDuration(hours): 0.073  # Run / Temp. data #  RangeInHours: 0.059 RangeInDegrees: 0.01 Nmeas:   12 
Run:  403 NEvents: 52036 RunDuration(hours): 0.098  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.03 Nmeas:   25 
Run:  404 NEvents: 52061 RunDuration(hours): 0.090  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.03 Nmeas:   30 
Run:  405 NEvents: 53503 RunDuration(hours): 0.108  # Run / Temp. data #  RangeInHours: 0.099 RangeInDegrees: 0.03 Nmeas:   25 
Run:  406 NEvents: 50966 RunDuration(hours): 0.108  # Run / Temp. data #  RangeInHours: 0.056 RangeInDegrees: 0.01 Nmeas:   16 
Run:  407 NEvents: 42508 RunDuration(hours): 0.147  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.03 Nmeas:   41 
Run:  408 NEvents: 58948 RunDuration(hours): 0.154  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.03 Nmeas:   44 
Run:  409 NEvents: 50475 RunDuration(hours): 0.139  # Run / Temp. data #  RangeInHours: 0.136 RangeInDegrees: 0.03 Nmeas:   35 
Run:  410 NEvents: 41365 RunDuration(hours): 0.195  # Run / Temp. data #  RangeInHours: 0.191 RangeInDegrees: 0.03 Nmeas:   54 
Run:  412 NEvents: 53020 RunDuration(hours): 0.217  # Run / Temp. data #  RangeInHours: 0.215 RangeInDegrees: 0.03 Nmeas:   73 
Run:  413 NEvents: 51091 RunDuration(hours): 0.083  # Run / Temp. data #  RangeInHours: 0.074 RangeInDegrees: 0.03 Nmeas:   29 
Run:  414 NEvents: 48794 RunDuration(hours): 0.082  # Run / Temp. data #  RangeInHours: 0.074 RangeInDegrees: 0.01 Nmeas:   26 
Run:  415 NEvents:   217 RunDuration(hours): 0.024  # Run / Temp. data #  RangeInHours: 0.020 RangeInDegrees: 0.01 Nmeas:    7 
Run:  416 NEvents: 57979 RunDuration(hours): 0.095  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.03 Nmeas:   34 
Run:  417 NEvents: 31012 RunDuration(hours): 0.094  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.01 Nmeas:   28 
Run:  418 NEvents: 51301 RunDuration(hours): 0.093  # Run / Temp. data #  RangeInHours: 0.090 RangeInDegrees: 0.03 Nmeas:   31 
Run:  419 NEvents:  2548 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.020 RangeInDegrees: 0.01 Nmeas:    9 
Run:  420 NEvents:  5668 RunDuration(hours): 0.041  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.01 Nmeas:   15 
Run:  421 NEvents: 52702 RunDuration(hours): 0.077  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.03 Nmeas:   31 
Run:  422 NEvents: 54053 RunDuration(hours): 0.090  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.03 Nmeas:   30 
Run:  423 NEvents: 49922 RunDuration(hours): 0.101  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.03 Nmeas:   30 
Run:  424 NEvents: 27601 RunDuration(hours): 0.060  # Run / Temp. data #  RangeInHours: 0.054 RangeInDegrees: 0.01 Nmeas:   11 
Run:  425 NEvents: 42854 RunDuration(hours): 0.091  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.03 Nmeas:   23 
Run:  426 NEvents: 50442 RunDuration(hours): 0.104  # Run / Temp. data #  RangeInHours: 0.098 RangeInDegrees: 0.03 Nmeas:   34 
Run:  427 NEvents: 51584 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.101 RangeInDegrees: 0.03 Nmeas:   25 
Run:  433 NEvents:  1000 RunDuration(hours): 0.012  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  434 NEvents:   322 RunDuration(hours): 0.023  # Run / Temp. data #  RangeInHours: 0.014 RangeInDegrees: 0.01 Nmeas:    4 
Run:  435 NEvents: 51561 RunDuration(hours): 0.159  # Run / Temp. data #  RangeInHours: 0.152 RangeInDegrees: 0.04 Nmeas:   46 
Run:  436 NEvents: 52887 RunDuration(hours): 0.165  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.04 Nmeas:   39 
Run:  437 NEvents: 50281 RunDuration(hours): 0.231  # Run / Temp. data #  RangeInHours: 0.224 RangeInDegrees: 0.05 Nmeas:   52 
Run:  438 NEvents: 53126 RunDuration(hours): 0.359  # Run / Temp. data #  RangeInHours: 0.357 RangeInDegrees: 0.06 Nmeas:   94 
Run:  439 NEvents: 50801 RunDuration(hours): 0.670  # Run / Temp. data #  RangeInHours: 0.664 RangeInDegrees: 0.09 Nmeas:  210 
Run:  440 NEvents: 50436 RunDuration(hours): 0.663  # Run / Temp. data #  RangeInHours: 0.654 RangeInDegrees: 0.08 Nmeas:  187 
Run:  441 NEvents:   207 RunDuration(hours): 0.016  # Run / Temp. data #  RangeInHours: 0.016 RangeInDegrees: 0.01 Nmeas:    7 
Run:  443 NEvents: 26880 RunDuration(hours): 0.060  # Run / Temp. data #  RangeInHours: 0.022 RangeInDegrees: 0.03 Nmeas:    9 
Run:  444 NEvents: 12077 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    5 
Run:  445 NEvents: 20589 RunDuration(hours): 0.037  # Run / Temp. data #  RangeInHours: 0.036 RangeInDegrees: 0.03 Nmeas:   12 
Run:  446 NEvents: 19884 RunDuration(hours): 0.046  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  447 NEvents: 10522 RunDuration(hours): 0.046  # Run / Temp. data #  RangeInHours: 0.033 RangeInDegrees: 0.01 Nmeas:   12 
Run:  448 NEvents: 30263 RunDuration(hours): 0.053  # Run / Temp. data #  RangeInHours: 0.013 RangeInDegrees: 0.01 Nmeas:    4 
Run:  449 NEvents:  1000 RunDuration(hours): 0.008  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  450 NEvents: 51450 RunDuration(hours): 0.091  # Run / Temp. data #  RangeInHours: 0.067 RangeInDegrees: 0.01 Nmeas:   17 
Run:  451 NEvents: 50419 RunDuration(hours): 0.089  # Run / Temp. data #  RangeInHours: 0.088 RangeInDegrees: 0.03 Nmeas:    7 
Run:  452 NEvents: 50235 RunDuration(hours): 0.102  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.02 Nmeas:   13 
Run:  453 NEvents: 51935 RunDuration(hours): 0.103  # Run / Temp. data #  RangeInHours: 0.079 RangeInDegrees: 0.01 Nmeas:   25 
Run:  454 NEvents: 50796 RunDuration(hours): 0.141  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   18 
Run:  455 NEvents: 55357 RunDuration(hours): 0.154  # Run / Temp. data #  RangeInHours: 0.154 RangeInDegrees: 0.03 Nmeas:   16 
Run:  456 NEvents: 50191 RunDuration(hours): 0.138  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.03 Nmeas:   23 
Run:  457 NEvents: 51977 RunDuration(hours): 0.250  # Run / Temp. data #  RangeInHours: 0.215 RangeInDegrees: 0.01 Nmeas:   22 
Run:  458 NEvents: 53359 RunDuration(hours): 0.232  # Run / Temp. data #  RangeInHours: 0.213 RangeInDegrees: 0.03 Nmeas:   39 
Run:  459 NEvents:  9808 RunDuration(hours): 0.194  # Run / Temp. data #  RangeInHours: 0.159 RangeInDegrees: 0.03 Nmeas:   24 
Run:  460 NEvents: 35928 RunDuration(hours): 0.814  # Run / Temp. data #  RangeInHours: 0.802 RangeInDegrees: 0.04 Nmeas:  175 
Run:  461 NEvents: 50098 RunDuration(hours): 1.081  # Run / Temp. data #  RangeInHours: 1.065 RangeInDegrees: 0.04 Nmeas:  159 
Run:  462 NEvents:   300 RunDuration(hours): 0.013  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  463 NEvents:  6893 RunDuration(hours): 0.038  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  464 NEvents: 14639 RunDuration(hours): 0.042  # Run / Temp. data #  RangeInHours: 0.034 RangeInDegrees: 0.01 Nmeas:   10 
Run:  465 NEvents: 14679 RunDuration(hours): 0.037  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  466 NEvents: 64133 RunDuration(hours): 0.159  # Run / Temp. data #  RangeInHours: 0.105 RangeInDegrees: 0.01 Nmeas:   18 
Run:  468 NEvents: 21303 RunDuration(hours): 0.230  # Run / Temp. data #  RangeInHours: 0.226 RangeInDegrees: 0.04 Nmeas:   38 
Run:  469 NEvents: 10023 RunDuration(hours): 0.132  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.04 Nmeas:   27 
Run:  498 NEvents:   100 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  500 NEvents:   100 RunDuration(hours): 0.013  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  501 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  504 NEvents: 27256 RunDuration(hours): 0.080  # Run / Temp. data #  RangeInHours: 0.077 RangeInDegrees: 0.03 Nmeas:   29 
Run:  509 NEvents: 27357 RunDuration(hours): 0.086  # Run / Temp. data #  RangeInHours: 0.070 RangeInDegrees: 0.01 Nmeas:   29 
Run:  510 NEvents: 27144 RunDuration(hours): 0.094  # Run / Temp. data #  RangeInHours: 0.087 RangeInDegrees: 0.01 Nmeas:   29 
Run:  511 NEvents: 18566 RunDuration(hours): 0.077  # Run / Temp. data #  RangeInHours: 0.068 RangeInDegrees: 0.01 Nmeas:   18 
Run:  512 NEvents: 28193 RunDuration(hours): 0.107  # Run / Temp. data #  RangeInHours: 0.099 RangeInDegrees: 0.03 Nmeas:   19 
Run:  514 NEvents:  9751 RunDuration(hours): 0.137  # Run / Temp. data #  RangeInHours: 0.107 RangeInDegrees: 0.01 Nmeas:    6 
Run:  517 NEvents: 27848 RunDuration(hours): 0.070  # Run / Temp. data #  RangeInHours: 0.014 RangeInDegrees: 0.01 Nmeas:    5 
Run:  518 NEvents:  9337 RunDuration(hours): 0.077  # Run / Temp. data #  RangeInHours: 0.045 RangeInDegrees: 0.01 Nmeas:   13 
Run:  519 NEvents: 28910 RunDuration(hours): 0.153  # Run / Temp. data #  RangeInHours: 0.081 RangeInDegrees: 0.03 Nmeas:   21 
Run:  520 NEvents: 28223 RunDuration(hours): 0.123  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.01 Nmeas:   19 
Run:  522 NEvents: 27777 RunDuration(hours): 0.156  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   23 
Run:  523 NEvents: 25658 RunDuration(hours): 0.116  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.02 Nmeas:   17 
Run:  524 NEvents: 25822 RunDuration(hours): 0.068  # Run / Temp. data #  RangeInHours: 0.042 RangeInDegrees: 0.03 Nmeas:    8 
Run:  525 NEvents: 25610 RunDuration(hours): 0.085  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.01 Nmeas:   15 
Run:  526 NEvents: 27187 RunDuration(hours): 0.087  # Run / Temp. data #  RangeInHours: 0.054 RangeInDegrees: 0.03 Nmeas:   16 
Run:  528 NEvents: 25753 RunDuration(hours): 0.074  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.01 Nmeas:   11 
Run:  529 NEvents: 26006 RunDuration(hours): 0.068  # Run / Temp. data #  RangeInHours: 0.065 RangeInDegrees: 0.03 Nmeas:   16 
Run:  531 NEvents: 30464 RunDuration(hours): 0.077  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.04 Nmeas:   28 
Run:  533 NEvents: 36569 RunDuration(hours): 0.096  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.03 Nmeas:   35 
Run:  537 NEvents: 16833 RunDuration(hours): 0.075  # Run / Temp. data #  RangeInHours: 0.060 RangeInDegrees: 0.03 Nmeas:   26 
Run:  538 NEvents: 23636 RunDuration(hours): 0.068  # Run / Temp. data #  RangeInHours: 0.065 RangeInDegrees: 0.03 Nmeas:   17 
Run:  539 NEvents:  6687 RunDuration(hours): 0.034  # Run / Temp. data #  RangeInHours: 0.033 RangeInDegrees: 0.03 Nmeas:   16 
Run:  540 NEvents: 25726 RunDuration(hours): 0.074  # Run / Temp. data #  RangeInHours: 0.069 RangeInDegrees: 0.03 Nmeas:   25 
Run:  541 NEvents: 32451 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.03 Nmeas:   39 
Run:  542 NEvents:  9853 RunDuration(hours): 0.041  # Run / Temp. data #  RangeInHours: 0.031 RangeInDegrees: 0.03 Nmeas:    5 
Run:  543 NEvents: 31250 RunDuration(hours): 0.084  # Run / Temp. data #  RangeInHours: 0.070 RangeInDegrees: 0.03 Nmeas:   17 
Run:  544 NEvents:  5516 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.027 RangeInDegrees: 0.03 Nmeas:    7 
Run:  545 NEvents: 32074 RunDuration(hours): 0.086  # Run / Temp. data #  RangeInHours: 0.069 RangeInDegrees: 0.01 Nmeas:   19 
Run:  546 NEvents: 31865 RunDuration(hours): 0.128  # Run / Temp. data #  RangeInHours: 0.072 RangeInDegrees: 0.03 Nmeas:   10 
Run:  547 NEvents: 28104 RunDuration(hours): 0.129  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.03 Nmeas:   19 
Run:  548 NEvents: 30692 RunDuration(hours): 0.090  # Run / Temp. data #  RangeInHours: 0.079 RangeInDegrees: 0.01 Nmeas:   28 
Run:  550 NEvents: 14984 RunDuration(hours): 0.192  # Run / Temp. data #  RangeInHours: 0.148 RangeInDegrees: 0.03 Nmeas:   22 
Run:  551 NEvents: 54727 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.01 Nmeas:   14 
Run:  552 NEvents: 52772 RunDuration(hours): 0.139  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.01 Nmeas:   31 
Run:  553 NEvents: 56816 RunDuration(hours): 0.190  # Run / Temp. data #  RangeInHours: 0.175 RangeInDegrees: 0.03 Nmeas:   10 
Run:  554 NEvents: 57409 RunDuration(hours): 0.263  # Run / Temp. data #  RangeInHours: 0.247 RangeInDegrees: 0.03 Nmeas:   34 
Run:  555 NEvents: 57952 RunDuration(hours): 0.631  # Run / Temp. data #  RangeInHours: 0.628 RangeInDegrees: 0.03 Nmeas:  171 
Run:  556 NEvents: 50058 RunDuration(hours): 0.781  # Run / Temp. data #  RangeInHours: 0.769 RangeInDegrees: 0.03 Nmeas:  128 
Run:  557 NEvents: 20579 RunDuration(hours): 1.015  # Run / Temp. data #  RangeInHours: 1.011 RangeInDegrees: 0.05 Nmeas:  162 
Run:  558 NEvents: 19662 RunDuration(hours): 0.824  # Run / Temp. data #  RangeInHours: 0.791 RangeInDegrees: 0.04 Nmeas:  161 
Run:  559 NEvents:  8231 RunDuration(hours): 0.431  # Run / Temp. data #  RangeInHours: 0.428 RangeInDegrees: 0.03 Nmeas:   90 
Run:  560 NEvents: 27072 RunDuration(hours): 0.463  # Run / Temp. data #  RangeInHours: 0.460 RangeInDegrees: 0.04 Nmeas:  133 
Run:  561 NEvents: 16976 RunDuration(hours): 0.151  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.03 Nmeas:   50 
Run:  563 NEvents: 56112 RunDuration(hours): 0.166  # Run / Temp. data #  RangeInHours: 0.163 RangeInDegrees: 0.03 Nmeas:   30 
Run:  564 NEvents: 18587 RunDuration(hours): 0.226  # Run / Temp. data #  RangeInHours: 0.219 RangeInDegrees: 0.03 Nmeas:   95 
Run:  565 NEvents: 52877 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.215 RangeInDegrees: 0.03 Nmeas:   67 
Run:  567 NEvents: 47556 RunDuration(hours): 0.231  # Run / Temp. data #  RangeInHours: 0.226 RangeInDegrees: 0.04 Nmeas:   66 
Run:  568 NEvents: 50026 RunDuration(hours): 0.214  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.04 Nmeas:   65 
Run:  569 NEvents: 17593 RunDuration(hours): 0.132  # Run / Temp. data #  RangeInHours: 0.124 RangeInDegrees: 0.03 Nmeas:   42 
Run:  570 NEvents: 35607 RunDuration(hours): 0.202  # Run / Temp. data #  RangeInHours: 0.197 RangeInDegrees: 0.03 Nmeas:   68 
Run:  572 NEvents:  8040 RunDuration(hours): 0.056  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.01 Nmeas:   13 
Run:  576 NEvents: 34003 RunDuration(hours): 0.099  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.03 Nmeas:   29 
Run:  577 NEvents: 33595 RunDuration(hours): 0.085  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.03 Nmeas:   31 
Run:  578 NEvents: 25288 RunDuration(hours): 0.055  # Run / Temp. data #  RangeInHours: 0.051 RangeInDegrees: 0.04 Nmeas:   10 
Run:  580 NEvents: 51674 RunDuration(hours): 0.155  # Run / Temp. data #  RangeInHours: 0.152 RangeInDegrees: 0.04 Nmeas:   38 
Run:  581 NEvents: 50596 RunDuration(hours): 0.218  # Run / Temp. data #  RangeInHours: 0.206 RangeInDegrees: 0.05 Nmeas:   43 
Run:  582 NEvents: 53112 RunDuration(hours): 0.334  # Run / Temp. data #  RangeInHours: 0.334 RangeInDegrees: 0.04 Nmeas:   72 
Run:  583 NEvents: 35267 RunDuration(hours): 0.441  # Run / Temp. data #  RangeInHours: 0.392 RangeInDegrees: 0.03 Nmeas:   81 
Run:  584 NEvents: 24993 RunDuration(hours): 1.060  # Run / Temp. data #  RangeInHours: 1.060 RangeInDegrees: 0.07 Nmeas:  190 
Run:  585 NEvents: 10436 RunDuration(hours): 0.027  # Run / Temp. data #  RangeInHours: 0.011 RangeInDegrees: 0.02 Nmeas:    4 
Run:  586 NEvents: 29150 RunDuration(hours): 0.073  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.03 Nmeas:   14 
Run:  587 NEvents: 17533 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.005 RangeInDegrees: 0.02 Nmeas:    4 
Run:  589 NEvents: 23292 RunDuration(hours): 0.061  # Run / Temp. data #  RangeInHours: 0.059 RangeInDegrees: 0.03 Nmeas:   12 
Run:  590 NEvents: 19941 RunDuration(hours): 0.062  # Run / Temp. data #  RangeInHours: 0.033 RangeInDegrees: 0.01 Nmeas:    8 
Run:  591 NEvents:  8077 RunDuration(hours): 0.025  # Run / Temp. data #  RangeInHours: 0.024 RangeInDegrees: 0.01 Nmeas:    5 
Run:  593 NEvents: 28203 RunDuration(hours): 0.133  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  594 NEvents: 24296 RunDuration(hours): 0.074  # Run / Temp. data #  RangeInHours: 0.071 RangeInDegrees: 0.03 Nmeas:   17 
Run:  595 NEvents: 50120 RunDuration(hours): 0.223  # Run / Temp. data #  RangeInHours: 0.134 RangeInDegrees: 0.01 Nmeas:   36 
Run:  597 NEvents: 49993 RunDuration(hours): 0.298  # Run / Temp. data #  RangeInHours: 0.265 RangeInDegrees: 0.03 Nmeas:   29 
Run:  598 NEvents: 21491 RunDuration(hours): 0.284  # Run / Temp. data #  RangeInHours: 0.281 RangeInDegrees: 0.04 Nmeas:   66 
Run:  599 NEvents:  5316 RunDuration(hours): 0.081  # Run / Temp. data #  RangeInHours: 0.074 RangeInDegrees: 0.01 Nmeas:   21 
Run:  600 NEvents: 23519 RunDuration(hours): 0.834  # Run / Temp. data #  RangeInHours: 0.782 RangeInDegrees: 0.05 Nmeas:  111 
Run:  601 NEvents: 18085 RunDuration(hours): 0.048  # Run / Temp. data #  RangeInHours: 0.045 RangeInDegrees: 0.03 Nmeas:    9 
Run:  605 NEvents:   526 RunDuration(hours): 0.075  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.01 Nmeas:    8 
Run:  606 NEvents:   529 RunDuration(hours): 0.075  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  608 NEvents:   532 RunDuration(hours): 0.065  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    1 
Run:  609 NEvents:   518 RunDuration(hours): 0.047  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.03 Nmeas:    7 
Run:  613 NEvents:   540 RunDuration(hours): 0.050  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.03 Nmeas:   16 
Run:  614 NEvents:   532 RunDuration(hours): 0.030  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.01 Nmeas:    8 
Run:  615 NEvents:   515 RunDuration(hours): 0.028  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    6 
Run:  616 NEvents:   521 RunDuration(hours): 0.027  # Run / Temp. data #  RangeInHours: 0.023 RangeInDegrees: 0.01 Nmeas:    7 
Run:  617 NEvents:   501 RunDuration(hours): 0.032  # Run / Temp. data #  RangeInHours: 0.031 RangeInDegrees: 0.01 Nmeas:    5 
Run:  618 NEvents:   500 RunDuration(hours): 0.026  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    6 
Run:  619 NEvents:   499 RunDuration(hours): 0.031  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  620 NEvents:   499 RunDuration(hours): 0.029  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  623 NEvents:   510 RunDuration(hours): 0.032  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  625 NEvents:   499 RunDuration(hours): 0.050  # Run / Temp. data #  RangeInHours: 0.045 RangeInDegrees: 0.01 Nmeas:    8 
Run:  627 NEvents:   499 RunDuration(hours): 0.384  # Run / Temp. data #  RangeInHours: 0.354 RangeInDegrees: 0.03 Nmeas:   40 
Run:  631 NEvents: 14776 RunDuration(hours): 0.140  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   42 
Run:  632 NEvents: 34720 RunDuration(hours): 0.108  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    4 
Run:  633 NEvents: 52789 RunDuration(hours): 0.200  # Run / Temp. data #  RangeInHours: 0.191 RangeInDegrees: 0.03 Nmeas:   40 
Run:  636 NEvents: 25611 RunDuration(hours): 0.521  # Run / Temp. data #  RangeInHours: 0.516 RangeInDegrees: 0.05 Nmeas:  170 
Run:  637 NEvents: 25051 RunDuration(hours): 0.437  # Run / Temp. data #  RangeInHours: 0.435 RangeInDegrees: 0.05 Nmeas:  153 
Run:  641 NEvents:  3514 RunDuration(hours): 0.063  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.03 Nmeas:   23 
Run:  642 NEvents: 25412 RunDuration(hours): 0.458  # Run / Temp. data #  RangeInHours: 0.449 RangeInDegrees: 0.05 Nmeas:  126 
Run:  643 NEvents: 25688 RunDuration(hours): 0.432  # Run / Temp. data #  RangeInHours: 0.424 RangeInDegrees: 0.05 Nmeas:   75 
Run:  644 NEvents: 17860 RunDuration(hours): 0.192  # Run / Temp. data #  RangeInHours: 0.184 RangeInDegrees: 0.03 Nmeas:   38 
Run:  646 NEvents: 56035 RunDuration(hours): 0.378  # Run / Temp. data #  RangeInHours: 0.343 RangeInDegrees: 0.04 Nmeas:   54 
Run:  647 NEvents: 21702 RunDuration(hours): 0.141  # Run / Temp. data #  RangeInHours: 0.112 RangeInDegrees: 0.03 Nmeas:   22 
Run:  648 NEvents: 55923 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.04 Nmeas:   38 
Run:  649 NEvents: 56863 RunDuration(hours): 0.288  # Run / Temp. data #  RangeInHours: 0.282 RangeInDegrees: 0.03 Nmeas:   71 
Run:  650 NEvents: 50921 RunDuration(hours): 0.181  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.03 Nmeas:   46 
Run:  652 NEvents: 50630 RunDuration(hours): 1.610  # Run / Temp. data #  RangeInHours: 1.609 RangeInDegrees: 0.06 Nmeas:  467 
Run:  653 NEvents: 51346 RunDuration(hours): 1.592  # Run / Temp. data #  RangeInHours: 1.576 RangeInDegrees: 0.04 Nmeas:  393 
Run:  654 NEvents: 26929 RunDuration(hours): 1.449  # Run / Temp. data #  RangeInHours: 1.448 RangeInDegrees: 0.04 Nmeas:  390 
Run:  658 NEvents: 25276 RunDuration(hours): 1.303  # Run / Temp. data #  RangeInHours: 1.298 RangeInDegrees: 0.06 Nmeas:  255 
Run:  659 NEvents: 25141 RunDuration(hours): 0.399  # Run / Temp. data #  RangeInHours: 0.397 RangeInDegrees: 0.03 Nmeas:   80 
Run:  660 NEvents: 50791 RunDuration(hours): 0.414  # Run / Temp. data #  RangeInHours: 0.278 RangeInDegrees: 0.03 Nmeas:   48 
Run:  661 NEvents:  9913 RunDuration(hours): 0.059  # Run / Temp. data #  RangeInHours: 0.056 RangeInDegrees: 0.03 Nmeas:   14 
Run:  662 NEvents: 43244 RunDuration(hours): 0.233  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.04 Nmeas:   55 
Run:  663 NEvents: 51203 RunDuration(hours): 0.173  # Run / Temp. data #  RangeInHours: 0.079 RangeInDegrees: 0.03 Nmeas:   11 
Run:  668 NEvents: 30225 RunDuration(hours): 0.082  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.01 Nmeas:    8 
Run:  669 NEvents: 22406 RunDuration(hours): 0.047  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  671 NEvents: 52410 RunDuration(hours): 0.121  # Run / Temp. data #  RangeInHours: 0.054 RangeInDegrees: 0.03 Nmeas:    8 
Run:  672 NEvents: 52965 RunDuration(hours): 0.189  # Run / Temp. data #  RangeInHours: 0.159 RangeInDegrees: 0.03 Nmeas:   35 
Run:  673 NEvents: 50592 RunDuration(hours): 0.237  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.01 Nmeas:   28 
Run:  674 NEvents: 50906 RunDuration(hours): 0.349  # Run / Temp. data #  RangeInHours: 0.309 RangeInDegrees: 0.03 Nmeas:   26 
Run:  675 NEvents: 25554 RunDuration(hours): 0.335  # Run / Temp. data #  RangeInHours: 0.271 RangeInDegrees: 0.03 Nmeas:   37 
Run:  677 NEvents: 25186 RunDuration(hours): 1.124  # Run / Temp. data #  RangeInHours: 1.112 RangeInDegrees: 0.06 Nmeas:  181 
Run:  683 NEvents: 22353 RunDuration(hours): 0.049  # Run / Temp. data #  RangeInHours: 0.040 RangeInDegrees: 0.01 Nmeas:   11 
Run:  684 NEvents:  9033 RunDuration(hours): 0.021  # Run / Temp. data #  RangeInHours: 0.011 RangeInDegrees: 0.01 Nmeas:    4 
Run:  685 NEvents: 23601 RunDuration(hours): 0.056  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.01 Nmeas:   10 
Run:  688 NEvents: 16756 RunDuration(hours): 0.044  # Run / Temp. data #  RangeInHours: 0.027 RangeInDegrees: 0.01 Nmeas:    6 
Run:  690 NEvents: 25058 RunDuration(hours): 0.066  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.03 Nmeas:    5 
Run:  694 NEvents: 13774 RunDuration(hours): 0.032  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  696 NEvents: 20049 RunDuration(hours): 0.069  # Run / Temp. data #  RangeInHours: 0.040 RangeInDegrees: 0.03 Nmeas:   12 
Run:  697 NEvents: 29928 RunDuration(hours): 0.094  # Run / Temp. data #  RangeInHours: 0.094 RangeInDegrees: 0.03 Nmeas:   31 
Run:  698 NEvents: 48078 RunDuration(hours): 0.301  # Run / Temp. data #  RangeInHours: 0.106 RangeInDegrees: 0.03 Nmeas:   12 
Run:  699 NEvents: 11470 RunDuration(hours): 0.090  # Run / Temp. data #  RangeInHours: 0.047 RangeInDegrees: 0.01 Nmeas:    5 
Run:  702 NEvents: 25813 RunDuration(hours): 0.751  # Run / Temp. data #  RangeInHours: 0.748 RangeInDegrees: 0.03 Nmeas:  196 
Run:  704 NEvents: 40467 RunDuration(hours): 0.277  # Run / Temp. data #  RangeInHours: 0.271 RangeInDegrees: 0.01 Nmeas:   62 
Run:  705 NEvents: 25140 RunDuration(hours): 0.273  # Run / Temp. data #  RangeInHours: 0.168 RangeInDegrees: 0.01 Nmeas:   31 
Run:  706 NEvents: 20255 RunDuration(hours): 0.806  # Run / Temp. data #  RangeInHours: 0.673 RangeInDegrees: 0.04 Nmeas:  107 
Run:  707 NEvents: 25411 RunDuration(hours): 0.362  # Run / Temp. data #  RangeInHours: 0.255 RangeInDegrees: 0.01 Nmeas:   45 
Run:  708 NEvents: 51752 RunDuration(hours): 0.331  # Run / Temp. data #  RangeInHours: 0.266 RangeInDegrees: 0.03 Nmeas:   51 
Run:  710 NEvents: 19790 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.085 RangeInDegrees: 0.02 Nmeas:   22 
Run:  711 NEvents: 50051 RunDuration(hours): 0.223  # Run / Temp. data #  RangeInHours: 0.202 RangeInDegrees: 0.04 Nmeas:   36 
Run:  712 NEvents: 50285 RunDuration(hours): 0.179  # Run / Temp. data #  RangeInHours: 0.170 RangeInDegrees: 0.04 Nmeas:   42 
Run:  713 NEvents: 50131 RunDuration(hours): 0.131  # Run / Temp. data #  RangeInHours: 0.128 RangeInDegrees: 0.04 Nmeas:   42 
Run:  714 NEvents: 53060 RunDuration(hours): 0.107  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.03 Nmeas:   42 
Run:  715 NEvents: 51052 RunDuration(hours): 0.254  # Run / Temp. data #  RangeInHours: 0.251 RangeInDegrees: 0.05 Nmeas:   91 
Run:  725 NEvents: 50284 RunDuration(hours): 0.353  # Run / Temp. data #  RangeInHours: 0.339 RangeInDegrees: 0.04 Nmeas:   81 
Run:  726 NEvents: 50544 RunDuration(hours): 0.499  # Run / Temp. data #  RangeInHours: 0.473 RangeInDegrees: 0.05 Nmeas:   99 
Run:  727 NEvents: 50122 RunDuration(hours): 0.227  # Run / Temp. data #  RangeInHours: 0.215 RangeInDegrees: 0.04 Nmeas:   28 
Run:  728 NEvents: 54338 RunDuration(hours): 0.423  # Run / Temp. data #  RangeInHours: 0.412 RangeInDegrees: 0.04 Nmeas:   61 
Run:  734 NEvents: 50699 RunDuration(hours): 0.556  # Run / Temp. data #  RangeInHours: 0.551 RangeInDegrees: 0.04 Nmeas:   57 
Run:  735 NEvents: 58675 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.01 Nmeas:    6 
Run:  736 NEvents:  8606 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.059 RangeInDegrees: 0.03 Nmeas:    7 
Run:  737 NEvents: 53085 RunDuration(hours): 0.131  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    4 
Run:  738 NEvents: 51301 RunDuration(hours): 0.156  # Run / Temp. data #  RangeInHours: 0.156 RangeInDegrees: 0.01 Nmeas:   32 
Run:  739 NEvents: 50599 RunDuration(hours): 0.207  # Run / Temp. data #  RangeInHours: 0.193 RangeInDegrees: 0.03 Nmeas:    9 
Run:  740 NEvents: 52193 RunDuration(hours): 0.316  # Run / Temp. data #  RangeInHours: 0.273 RangeInDegrees: 0.02 Nmeas:   63 
Run:  741 NEvents: 52081 RunDuration(hours): 0.692  # Run / Temp. data #  RangeInHours: 0.673 RangeInDegrees: 0.03 Nmeas:   89 
Run:  742 NEvents: 29070 RunDuration(hours): 1.781  # Run / Temp. data #  RangeInHours: 1.692 RangeInDegrees: 0.04 Nmeas:  312 
Run:  743 NEvents: 69391 RunDuration(hours): 0.227  # Run / Temp. data #  RangeInHours: 0.224 RangeInDegrees: 0.02 Nmeas:   65 
Run:  744 NEvents: 65518 RunDuration(hours): 0.214  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.01 Nmeas:   28 
Run:  745 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  746 NEvents: 32497 RunDuration(hours): 0.107  # Run / Temp. data #  RangeInHours: 0.076 RangeInDegrees: 0.03 Nmeas:    8 
Run:  747 NEvents: 40365 RunDuration(hours): 0.143  # Run / Temp. data #  RangeInHours: 0.126 RangeInDegrees: 0.03 Nmeas:   28 
Run:  748 NEvents: 12925 RunDuration(hours): 1.638  # Run / Temp. data #  RangeInHours: 1.614 RangeInDegrees: 0.13 Nmeas:  415 
Run:  749 NEvents:   513 RunDuration(hours): 0.051  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.03 Nmeas:   13 
Run:  750 NEvents:   506 RunDuration(hours): 0.049  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.01 Nmeas:    8 
Run:  752 NEvents: 50348 RunDuration(hours): 0.166  # Run / Temp. data #  RangeInHours: 0.135 RangeInDegrees: 0.03 Nmeas:   28 
Run:  753 NEvents: 52399 RunDuration(hours): 0.177  # Run / Temp. data #  RangeInHours: 0.159 RangeInDegrees: 0.03 Nmeas:   21 
Run:  754 NEvents: 47597 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.03 Nmeas:   35 
Run:  755 NEvents: 51224 RunDuration(hours): 0.201  # Run / Temp. data #  RangeInHours: 0.181 RangeInDegrees: 0.03 Nmeas:   62 
Run:  756 NEvents: 29957 RunDuration(hours): 0.163  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.03 Nmeas:   17 
Run:  757 NEvents: 20325 RunDuration(hours): 0.072  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.01 Nmeas:   14 
Run:  758 NEvents: 52361 RunDuration(hours): 0.276  # Run / Temp. data #  RangeInHours: 0.247 RangeInDegrees: 0.01 Nmeas:   46 
Run:  759 NEvents: 52559 RunDuration(hours): 0.198  # Run / Temp. data #  RangeInHours: 0.124 RangeInDegrees: 0.03 Nmeas:   20 
Run:  760 NEvents: 53425 RunDuration(hours): 0.238  # Run / Temp. data #  RangeInHours: 0.229 RangeInDegrees: 0.03 Nmeas:   44 
Run:  761 NEvents: 53233 RunDuration(hours): 0.225  # Run / Temp. data #  RangeInHours: 0.197 RangeInDegrees: 0.03 Nmeas:   45 
Run:  762 NEvents: 19260 RunDuration(hours): 0.244  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   31 
Run:  763 NEvents: 35023 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.01 Nmeas:   23 
Run:  764 NEvents: 52405 RunDuration(hours): 0.174  # Run / Temp. data #  RangeInHours: 0.163 RangeInDegrees: 0.03 Nmeas:   19 
Run:  765 NEvents: 50148 RunDuration(hours): 0.194  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.04 Nmeas:   28 
Run:  766 NEvents: 27259 RunDuration(hours): 0.117  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.03 Nmeas:   18 
Run:  767 NEvents: 25401 RunDuration(hours): 0.086  # Run / Temp. data #  RangeInHours: 0.079 RangeInDegrees: 0.03 Nmeas:   21 
Run:  768 NEvents: 50831 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.161 RangeInDegrees: 0.04 Nmeas:   49 
Run:  769 NEvents: 52636 RunDuration(hours): 0.218  # Run / Temp. data #  RangeInHours: 0.217 RangeInDegrees: 0.06 Nmeas:   79 
Run:  770 NEvents: 30463 RunDuration(hours): 0.135  # Run / Temp. data #  RangeInHours: 0.132 RangeInDegrees: 0.03 Nmeas:   38 
Run:  771 NEvents: 21762 RunDuration(hours): 0.074  # Run / Temp. data #  RangeInHours: 0.069 RangeInDegrees: 0.03 Nmeas:   21 
Run:  772 NEvents: 51941 RunDuration(hours): 0.158  # Run / Temp. data #  RangeInHours: 0.154 RangeInDegrees: 0.04 Nmeas:   41 
Run:  773 NEvents: 32301 RunDuration(hours): 0.171  # Run / Temp. data #  RangeInHours: 0.163 RangeInDegrees: 0.04 Nmeas:   37 
Run:  775 NEvents: 20720 RunDuration(hours): 0.195  # Run / Temp. data #  RangeInHours: 0.191 RangeInDegrees: 0.04 Nmeas:   34 
Run:  776 NEvents: 51204 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.03 Nmeas:   22 
Run:  778 NEvents: 36376 RunDuration(hours): 0.304  # Run / Temp. data #  RangeInHours: 0.298 RangeInDegrees: 0.04 Nmeas:   69 
Run:  780 NEvents: 15679 RunDuration(hours): 0.057  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.03 Nmeas:   12 
Run:  781 NEvents: 50642 RunDuration(hours): 0.207  # Run / Temp. data #  RangeInHours: 0.195 RangeInDegrees: 0.04 Nmeas:   45 
Run:  782 NEvents: 39197 RunDuration(hours): 0.142  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.04 Nmeas:   41 
Run:  784 NEvents: 35270 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.04 Nmeas:   33 
Run:  785 NEvents: 32336 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.088 RangeInDegrees: 0.03 Nmeas:   23 
Run:  786 NEvents: 39189 RunDuration(hours): 0.127  # Run / Temp. data #  RangeInHours: 0.106 RangeInDegrees: 0.03 Nmeas:   20 
Run:  787 NEvents: 35543 RunDuration(hours): 0.121  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.04 Nmeas:   21 
Run:  789 NEvents: 35853 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.105 RangeInDegrees: 0.03 Nmeas:   16 
Run:  790 NEvents:  2858 RunDuration(hours): 0.055  # Run / Temp. data #  RangeInHours: 0.052 RangeInDegrees: 0.04 Nmeas:    9 
Run:  791 NEvents: 30618 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.108 RangeInDegrees: 0.03 Nmeas:   18 
Run:  792 NEvents: 35922 RunDuration(hours): 0.112  # Run / Temp. data #  RangeInHours: 0.112 RangeInDegrees: 0.03 Nmeas:   26 
Run:  794 NEvents: 51098 RunDuration(hours): 0.160  # Run / Temp. data #  RangeInHours: 0.159 RangeInDegrees: 0.04 Nmeas:   43 
Run:  795 NEvents: 58855 RunDuration(hours): 0.217  # Run / Temp. data #  RangeInHours: 0.213 RangeInDegrees: 0.04 Nmeas:   40 
Run:  796 NEvents: 50530 RunDuration(hours): 0.176  # Run / Temp. data #  RangeInHours: 0.175 RangeInDegrees: 0.03 Nmeas:   47 
Run:  797 NEvents: 28315 RunDuration(hours): 0.250  # Run / Temp. data #  RangeInHours: 0.238 RangeInDegrees: 0.03 Nmeas:   45 
Run:  800 NEvents: 53653 RunDuration(hours): 0.178  # Run / Temp. data #  RangeInHours: 0.177 RangeInDegrees: 0.03 Nmeas:   53 
Run:  801 NEvents:  5781 RunDuration(hours): 0.067  # Run / Temp. data #  RangeInHours: 0.067 RangeInDegrees: 0.03 Nmeas:   11 
Run:  802 NEvents: 50113 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.01 Nmeas:   15 
Run:  804 NEvents: 48506 RunDuration(hours): 0.209  # Run / Temp. data #  RangeInHours: 0.128 RangeInDegrees: 0.01 Nmeas:   12 
Run:  805 NEvents: 57239 RunDuration(hours): 0.241  # Run / Temp. data #  RangeInHours: 0.119 RangeInDegrees: 0.01 Nmeas:    6 
Run:  806 NEvents: 64997 RunDuration(hours): 0.206  # Run / Temp. data #  RangeInHours: 0.177 RangeInDegrees: 0.03 Nmeas:   12 
Run:  807 NEvents: 52250 RunDuration(hours): 0.377  # Run / Temp. data #  RangeInHours: 0.336 RangeInDegrees: 0.04 Nmeas:   67 
Run:  808 NEvents: 49935 RunDuration(hours): 0.369  # Run / Temp. data #  RangeInHours: 0.186 RangeInDegrees: 0.03 Nmeas:   41 
Run:  809 NEvents: 51291 RunDuration(hours): 0.362  # Run / Temp. data #  RangeInHours: 0.357 RangeInDegrees: 0.03 Nmeas:   69 
Run:  810 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  813 NEvents: 50240 RunDuration(hours): 1.507  # Run / Temp. data #  RangeInHours: 1.504 RangeInDegrees: 0.12 Nmeas:  492 
Run:  815 NEvents: 52909 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.099 RangeInDegrees: 0.03 Nmeas:   17 
Run:  816 NEvents: 50772 RunDuration(hours): 0.206  # Run / Temp. data #  RangeInHours: 0.202 RangeInDegrees: 0.04 Nmeas:   32 
Run:  817 NEvents: 50817 RunDuration(hours): 0.154  # Run / Temp. data #  RangeInHours: 0.152 RangeInDegrees: 0.03 Nmeas:   21 
Run:  818 NEvents: 50848 RunDuration(hours): 0.196  # Run / Temp. data #  RangeInHours: 0.182 RangeInDegrees: 0.03 Nmeas:   30 
Run:  819 NEvents: 51048 RunDuration(hours): 0.305  # Run / Temp. data #  RangeInHours: 0.264 RangeInDegrees: 0.04 Nmeas:   39 
Run:  820 NEvents: 50783 RunDuration(hours): 0.578  # Run / Temp. data #  RangeInHours: 0.556 RangeInDegrees: 0.05 Nmeas:   73 
Run:  822 NEvents: 25797 RunDuration(hours): 1.243  # Run / Temp. data #  RangeInHours: 1.242 RangeInDegrees: 0.13 Nmeas:  287 
Run:  823 NEvents: 54005 RunDuration(hours): 0.146  # Run / Temp. data #  RangeInHours: 0.141 RangeInDegrees: 0.04 Nmeas:   40 
Run:  824 NEvents: 37619 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.04 Nmeas:   39 
Run:  825 NEvents: 16934 RunDuration(hours): 0.049  # Run / Temp. data #  RangeInHours: 0.047 RangeInDegrees: 0.01 Nmeas:   16 
Run:  826 NEvents: 51412 RunDuration(hours): 0.165  # Run / Temp. data #  RangeInHours: 0.159 RangeInDegrees: 0.04 Nmeas:   47 
Run:  827 NEvents: 51384 RunDuration(hours): 0.219  # Run / Temp. data #  RangeInHours: 0.200 RangeInDegrees: 0.04 Nmeas:   37 
Run:  828 NEvents: 50799 RunDuration(hours): 0.311  # Run / Temp. data #  RangeInHours: 0.276 RangeInDegrees: 0.04 Nmeas:   66 
Run:  829 NEvents: 57694 RunDuration(hours): 0.658  # Run / Temp. data #  RangeInHours: 0.646 RangeInDegrees: 0.09 Nmeas:   87 
Run:  830 NEvents: 25128 RunDuration(hours): 0.936  # Run / Temp. data #  RangeInHours: 0.929 RangeInDegrees: 0.12 Nmeas:  162 
Run:  831 NEvents: 17367 RunDuration(hours): 0.061  # Run / Temp. data #  RangeInHours: 0.058 RangeInDegrees: 0.03 Nmeas:   17 
Run:  832 NEvents: 10099 RunDuration(hours): 0.034  # Run / Temp. data #  RangeInHours: 0.025 RangeInDegrees: 0.03 Nmeas:    6 
Run:  833 NEvents: 26551 RunDuration(hours): 0.059  # Run / Temp. data #  RangeInHours: 0.058 RangeInDegrees: 0.03 Nmeas:   18 
Run:  834 NEvents: 51006 RunDuration(hours): 0.149  # Run / Temp. data #  RangeInHours: 0.148 RangeInDegrees: 0.04 Nmeas:   55 
Run:  835 NEvents: 50876 RunDuration(hours): 0.147  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.04 Nmeas:   56 
Run:  836 NEvents: 14162 RunDuration(hours): 0.071  # Run / Temp. data #  RangeInHours: 0.069 RangeInDegrees: 0.03 Nmeas:   20 
Run:  837 NEvents: 38686 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.04 Nmeas:   54 
Run:  838 NEvents: 50717 RunDuration(hours): 0.186  # Run / Temp. data #  RangeInHours: 0.184 RangeInDegrees: 0.03 Nmeas:   54 
Run:  839 NEvents: 51510 RunDuration(hours): 0.375  # Run / Temp. data #  RangeInHours: 0.359 RangeInDegrees: 0.06 Nmeas:   88 
Run:  840 NEvents: 50457 RunDuration(hours): 0.578  # Run / Temp. data #  RangeInHours: 0.571 RangeInDegrees: 0.08 Nmeas:  123 
Run:  841 NEvents: 17973 RunDuration(hours): 0.749  # Run / Temp. data #  RangeInHours: 0.742 RangeInDegrees: 0.20 Nmeas:  152 
Run:  842 NEvents:  6349 RunDuration(hours): 0.472  # Run / Temp. data #  RangeInHours: 0.462 RangeInDegrees: 0.41 Nmeas:  100 
Run:  843 NEvents:   100 RunDuration(hours): 0.012  # Run / Temp. data #  RangeInHours: 0.011 RangeInDegrees: 0.04 Nmeas:    5 
Run:  844 NEvents:   126 RunDuration(hours): 0.019  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.01 Nmeas:    5 
Run:  846 NEvents:   537 RunDuration(hours): 0.015  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  847 NEvents:   702 RunDuration(hours): 0.011  # Run / Temp. data #  RangeInHours: 0.004 RangeInDegrees: 0.01 Nmeas:    3 
Run:  848 NEvents:   677 RunDuration(hours): 0.010  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  850 NEvents:   714 RunDuration(hours): 0.007  # Run / Temp. data #  RangeInHours: 0.007 RangeInDegrees: 0.03 Nmeas:    5 
Run:  852 NEvents:   577 RunDuration(hours): 0.010  # Run / Temp. data #  RangeInHours: 0.004 RangeInDegrees: 0.03 Nmeas:    3 
Run:  853 NEvents:   520 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.009 RangeInDegrees: 0.03 Nmeas:    6 
Run:  854 NEvents:   592 RunDuration(hours): 0.008  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    1 
Run:  855 NEvents:   627 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  856 NEvents:   686 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    1 
Run:  857 NEvents:   541 RunDuration(hours): 0.005  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  858 NEvents:   521 RunDuration(hours): 0.008  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  859 NEvents:   520 RunDuration(hours): 0.006  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  860 NEvents:   571 RunDuration(hours): 0.005  # Run / Temp. data #  RangeInHours: 0.004 RangeInDegrees: 0.01 Nmeas:    2 
Run:  861 NEvents:   545 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  862 NEvents:   523 RunDuration(hours): 0.007  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  863 NEvents:   517 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.007 RangeInDegrees: 0.01 Nmeas:    3 
Run:  864 NEvents:   522 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  865 NEvents:   581 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.007 RangeInDegrees: 0.01 Nmeas:    3 
Run:  867 NEvents:   520 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  868 NEvents:   515 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    1 
Run:  869 NEvents:   525 RunDuration(hours): 0.010  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  870 NEvents:   513 RunDuration(hours): 0.009  # Run / Temp. data #  RangeInHours: 0.006 RangeInDegrees: 0.01 Nmeas:    4 
Run:  872 NEvents: 37857 RunDuration(hours): 0.168  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.03 Nmeas:   39 
Run:  873 NEvents: 74045 RunDuration(hours): 0.321  # Run / Temp. data #  RangeInHours: 0.319 RangeInDegrees: 0.04 Nmeas:   48 
Run:  874 NEvents: 10137 RunDuration(hours): 0.657  # Run / Temp. data #  RangeInHours: 0.639 RangeInDegrees: 0.07 Nmeas:  122 
Run:  875 NEvents: 10921 RunDuration(hours): 1.040  # Run / Temp. data #  RangeInHours: 1.036 RangeInDegrees: 0.09 Nmeas:  132 
Run:  876 NEvents: 10201 RunDuration(hours): 0.537  # Run / Temp. data #  RangeInHours: 0.486 RangeInDegrees: 0.03 Nmeas:   69 
Run:  877 NEvents: 14540 RunDuration(hours): 0.579  # Run / Temp. data #  RangeInHours: 0.549 RangeInDegrees: 0.03 Nmeas:   67 
Run:  878 NEvents: 13703 RunDuration(hours): 0.624  # Run / Temp. data #  RangeInHours: 0.581 RangeInDegrees: 0.03 Nmeas:  121 
Run:  879 NEvents: 15508 RunDuration(hours): 0.386  # Run / Temp. data #  RangeInHours: 0.173 RangeInDegrees: 0.01 Nmeas:   22 
Run:  880 NEvents: 28384 RunDuration(hours): 1.444  # Run / Temp. data #  RangeInHours: 1.426 RangeInDegrees: 0.05 Nmeas:  250 
Run:  881 NEvents: 27519 RunDuration(hours): 0.656  # Run / Temp. data #  RangeInHours: 0.654 RangeInDegrees: 0.04 Nmeas:  162 
Run:  882 NEvents: 27854 RunDuration(hours): 1.062  # Run / Temp. data #  RangeInHours: 1.058 RangeInDegrees: 0.08 Nmeas:  374 
Run:  883 NEvents: 46021 RunDuration(hours): 1.825  # Run / Temp. data #  RangeInHours: 1.774 RangeInDegrees: 0.17 Nmeas:  323 
Run:  884 NEvents: 25094 RunDuration(hours): 2.521  # Run / Temp. data #  RangeInHours: 2.511 RangeInDegrees: 0.22 Nmeas:  647 
Run:  885 NEvents: 37741 RunDuration(hours): 1.110  # Run / Temp. data #  RangeInHours: 1.109 RangeInDegrees: 0.13 Nmeas:  147 
Run:  887 NEvents: 14606 RunDuration(hours): 0.464  # Run / Temp. data #  RangeInHours: 0.457 RangeInDegrees: 0.03 Nmeas:   87 
Run:  888 NEvents:  8105 RunDuration(hours): 0.425  # Run / Temp. data #  RangeInHours: 0.424 RangeInDegrees: 0.04 Nmeas:  126 
Run:  889 NEvents: 25109 RunDuration(hours): 0.453  # Run / Temp. data #  RangeInHours: 0.451 RangeInDegrees: 0.05 Nmeas:  148 
Run:  890 NEvents: 17666 RunDuration(hours): 0.120  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.03 Nmeas:   41 
Run:  891 NEvents: 56289 RunDuration(hours): 0.936  # Run / Temp. data #  RangeInHours: 0.852 RangeInDegrees: 0.09 Nmeas:  156 
Run:  892 NEvents: 50081 RunDuration(hours): 0.777  # Run / Temp. data #  RangeInHours: 0.764 RangeInDegrees: 0.08 Nmeas:  108 
Run:  893 NEvents: 49479 RunDuration(hours): 0.834  # Run / Temp. data #  RangeInHours: 0.825 RangeInDegrees: 0.11 Nmeas:  198 
Run:  894 NEvents: 51593 RunDuration(hours): 0.856  # Run / Temp. data #  RangeInHours: 0.849 RangeInDegrees: 0.10 Nmeas:  247 
Run:  895 NEvents:  6069 RunDuration(hours): 0.501  # Run / Temp. data #  RangeInHours: 0.497 RangeInDegrees: 0.08 Nmeas:   93 
Run:  897 NEvents: 47149 RunDuration(hours): 0.855  # Run / Temp. data #  RangeInHours: 0.852 RangeInDegrees: 0.51 Nmeas:  211 
Run:  898 NEvents:   500 RunDuration(hours): 0.048  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.01 Nmeas:   10 
Run:  899 NEvents:   517 RunDuration(hours): 0.068  # Run / Temp. data #  RangeInHours: 0.065 RangeInDegrees: 0.03 Nmeas:   18 
Run:  900 NEvents: 48432 RunDuration(hours): 0.996  # Run / Temp. data #  RangeInHours: 0.982 RangeInDegrees: 0.18 Nmeas:  323 
Run:  901 NEvents: 52438 RunDuration(hours): 0.817  # Run / Temp. data #  RangeInHours: 0.811 RangeInDegrees: 0.12 Nmeas:  141 
Run:  902 NEvents: 51407 RunDuration(hours): 0.066  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.01 Nmeas:    4 
Run:  903 NEvents:  5598 RunDuration(hours): 0.052  # Run / Temp. data #  RangeInHours: 0.045 RangeInDegrees: 0.01 Nmeas:   10 
Run:  904 NEvents: 50272 RunDuration(hours): 0.072  # Run / Temp. data #  RangeInHours: 0.009 RangeInDegrees: 0.01 Nmeas:    4 
Run:  905 NEvents: 61007 RunDuration(hours): 0.130  # Run / Temp. data #  RangeInHours: 0.124 RangeInDegrees: 0.03 Nmeas:   19 
Run:  906 NEvents: 52457 RunDuration(hours): 0.259  # Run / Temp. data #  RangeInHours: 0.134 RangeInDegrees: 0.03 Nmeas:   35 
Run:  907 NEvents: 33059 RunDuration(hours): 0.961  # Run / Temp. data #  RangeInHours: 0.953 RangeInDegrees: 0.03 Nmeas:  258 
Run:  908 NEvents: 25126 RunDuration(hours): 0.935  # Run / Temp. data #  RangeInHours: 0.933 RangeInDegrees: 0.03 Nmeas:  224 
Run:  911 NEvents: 51162 RunDuration(hours): 0.249  # Run / Temp. data #  RangeInHours: 0.229 RangeInDegrees: 0.04 Nmeas:   45 
Run:  912 NEvents: 49387 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.110 RangeInDegrees: 0.03 Nmeas:   38 
Run:  914 NEvents: 13878 RunDuration(hours): 0.034  # Run / Temp. data #  RangeInHours: 0.031 RangeInDegrees: 0.03 Nmeas:   15 
Run:  916 NEvents: 24869 RunDuration(hours): 0.049  # Run / Temp. data #  RangeInHours: 0.043 RangeInDegrees: 0.01 Nmeas:   18 
Run:  917 NEvents: 16574 RunDuration(hours): 0.024  # Run / Temp. data #  RangeInHours: 0.022 RangeInDegrees: 0.03 Nmeas:    8 
Run:  920 NEvents: 24220 RunDuration(hours): 0.043  # Run / Temp. data #  RangeInHours: 0.040 RangeInDegrees: 0.03 Nmeas:   19 
Run:  921 NEvents: 24390 RunDuration(hours): 0.038  # Run / Temp. data #  RangeInHours: 0.029 RangeInDegrees: 0.03 Nmeas:   11 
Run:  923 NEvents: 24988 RunDuration(hours): 0.739  # Run / Temp. data #  RangeInHours: 0.735 RangeInDegrees: 0.04 Nmeas:  244 
Run:  924 NEvents: 51531 RunDuration(hours): 0.247  # Run / Temp. data #  RangeInHours: 0.247 RangeInDegrees: 0.03 Nmeas:  106 
Run:  927 NEvents: 54359 RunDuration(hours): 0.075  # Run / Temp. data #  RangeInHours: 0.065 RangeInDegrees: 0.03 Nmeas:   22 
Run:  930 NEvents: 50473 RunDuration(hours): 0.074  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.01 Nmeas:   18 
Run:  931 NEvents: 51423 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.03 Nmeas:   52 
Run:  932 NEvents: 50505 RunDuration(hours): 0.232  # Run / Temp. data #  RangeInHours: 0.228 RangeInDegrees: 0.03 Nmeas:   91 
Run:  933 NEvents: 25612 RunDuration(hours): 0.731  # Run / Temp. data #  RangeInHours: 0.729 RangeInDegrees: 0.05 Nmeas:  252 
Run:  934 NEvents: 13725 RunDuration(hours): 0.041  # Run / Temp. data #  RangeInHours: 0.033 RangeInDegrees: 0.03 Nmeas:   13 
Run:  935 NEvents: 36274 RunDuration(hours): 0.049  # Run / Temp. data #  RangeInHours: 0.040 RangeInDegrees: 0.03 Nmeas:   14 
Run:  936 NEvents: 49615 RunDuration(hours): 0.067  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.03 Nmeas:   19 
Run:  937 NEvents: 36786 RunDuration(hours): 0.096  # Run / Temp. data #  RangeInHours: 0.094 RangeInDegrees: 0.03 Nmeas:   17 
Run:  938 NEvents: 15117 RunDuration(hours): 0.032  # Run / Temp. data #  RangeInHours: 0.031 RangeInDegrees: 0.03 Nmeas:   12 
Run:  940 NEvents: 50776 RunDuration(hours): 0.248  # Run / Temp. data #  RangeInHours: 0.238 RangeInDegrees: 0.04 Nmeas:   66 
Run:  941 NEvents: 25011 RunDuration(hours): 0.128  # Run / Temp. data #  RangeInHours: 0.083 RangeInDegrees: 0.03 Nmeas:   22 
Run:  942 NEvents: 13575 RunDuration(hours): 0.044  # Run / Temp. data #  RangeInHours: 0.038 RangeInDegrees: 0.03 Nmeas:   13 
Run:  943 NEvents: 13894 RunDuration(hours): 0.058  # Run / Temp. data #  RangeInHours: 0.049 RangeInDegrees: 0.02 Nmeas:   11 
Run:  944 NEvents: 24691 RunDuration(hours): 0.064  # Run / Temp. data #  RangeInHours: 0.061 RangeInDegrees: 0.03 Nmeas:    8 
Run:  947 NEvents:  5681 RunDuration(hours): 0.029  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.01 Nmeas:    2 
Run:  948 NEvents: 30381 RunDuration(hours): 0.052  # Run / Temp. data #  RangeInHours: 0.051 RangeInDegrees: 0.03 Nmeas:   17 
Run:  949 NEvents: 49764 RunDuration(hours): 0.068  # Run / Temp. data #  RangeInHours: 0.041 RangeInDegrees: 0.03 Nmeas:    5 
Run:  950 NEvents: 49347 RunDuration(hours): 0.065  # Run / Temp. data #  RangeInHours: 0.060 RangeInDegrees: 0.03 Nmeas:   18 
Run:  951 NEvents: 49145 RunDuration(hours): 0.093  # Run / Temp. data #  RangeInHours: 0.067 RangeInDegrees: 0.01 Nmeas:    7 
Run:  955 NEvents: 51179 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.074 RangeInDegrees: 0.03 Nmeas:   22 
Run:  956 NEvents: 49720 RunDuration(hours): 0.346  # Run / Temp. data #  RangeInHours: 0.296 RangeInDegrees: 0.04 Nmeas:   52 
Run:  958 NEvents: 26172 RunDuration(hours): 0.113  # Run / Temp. data #  RangeInHours: 0.089 RangeInDegrees: 0.03 Nmeas:   10 
Run:  959 NEvents: 38422 RunDuration(hours): 0.771  # Run / Temp. data #  RangeInHours: 0.704 RangeInDegrees: 0.04 Nmeas:   56 
Run:  960 NEvents:   521 RunDuration(hours): 0.048  # Run / Temp. data #  RangeInHours: 0.016 RangeInDegrees: 0.01 Nmeas:    6 
Run:  961 NEvents:   552 RunDuration(hours): 0.050  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  966 NEvents:  1408 RunDuration(hours): 0.082  # Run / Temp. data #  RangeInHours: 0.079 RangeInDegrees: 0.04 Nmeas:   14 
Run:  968 NEvents:   552 RunDuration(hours): 0.030  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  969 NEvents:   551 RunDuration(hours): 0.032  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  970 NEvents:   648 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  974 NEvents:   591 RunDuration(hours): 0.037  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  975 NEvents:   503 RunDuration(hours): 0.036  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  976 NEvents:   505 RunDuration(hours): 0.031  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run:  978 NEvents:   503 RunDuration(hours): 0.030  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1000 NEvents:   529 RunDuration(hours): 0.013  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1001 NEvents: 25197 RunDuration(hours): 0.401  # Run / Temp. data #  RangeInHours: 0.397 RangeInDegrees: 0.27 Nmeas:   79 
Run: 1002 NEvents: 25107 RunDuration(hours): 0.454  # Run / Temp. data #  RangeInHours: 0.413 RangeInDegrees: 0.27 Nmeas:   89 
Run: 1003 NEvents:  1384 RunDuration(hours): 0.239  # Run / Temp. data #  RangeInHours: 0.237 RangeInDegrees: 0.27 Nmeas:   65 
Run: 1006 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1007 NEvents: 30988 RunDuration(hours): 0.274  # Run / Temp. data #  RangeInHours: 0.203 RangeInDegrees: 0.27 Nmeas:    8 
Run: 1010 NEvents:  1149 RunDuration(hours): 0.150  # Run / Temp. data #  RangeInHours: 0.141 RangeInDegrees: 0.13 Nmeas:   30 
Run: 1011 NEvents: 31325 RunDuration(hours): 0.309  # Run / Temp. data #  RangeInHours: 0.288 RangeInDegrees: 0.27 Nmeas:   22 
Run: 1015 NEvents:  3232 RunDuration(hours): 0.458  # Run / Temp. data #  RangeInHours: 0.452 RangeInDegrees: 0.27 Nmeas:   81 
Run: 1016 NEvents: 25065 RunDuration(hours): 0.223  # Run / Temp. data #  RangeInHours: 0.212 RangeInDegrees: 0.27 Nmeas:   63 
Run: 1017 NEvents:  3357 RunDuration(hours): 0.490  # Run / Temp. data #  RangeInHours: 0.456 RangeInDegrees: 0.27 Nmeas:   28 
Run: 1018 NEvents:  6212 RunDuration(hours): 0.444  # Run / Temp. data #  RangeInHours: 0.406 RangeInDegrees: 0.14 Nmeas:   39 
Run: 1019 NEvents: 25439 RunDuration(hours): 0.221  # Run / Temp. data #  RangeInHours: 0.217 RangeInDegrees: 0.14 Nmeas:   26 
Run: 1020 NEvents:  3042 RunDuration(hours): 0.429  # Run / Temp. data #  RangeInHours: 0.424 RangeInDegrees: 0.27 Nmeas:   91 
Run: 1021 NEvents: 25650 RunDuration(hours): 0.233  # Run / Temp. data #  RangeInHours: 0.219 RangeInDegrees: 0.27 Nmeas:   57 
Run: 1022 NEvents:  3026 RunDuration(hours): 0.446  # Run / Temp. data #  RangeInHours: 0.441 RangeInDegrees: 0.27 Nmeas:  123 
Run: 1023 NEvents: 25749 RunDuration(hours): 0.236  # Run / Temp. data #  RangeInHours: 0.235 RangeInDegrees: 0.27 Nmeas:   70 
Run: 1024 NEvents:  3029 RunDuration(hours): 0.452  # Run / Temp. data #  RangeInHours: 0.448 RangeInDegrees: 0.27 Nmeas:  136 
Run: 1025 NEvents: 25965 RunDuration(hours): 0.228  # Run / Temp. data #  RangeInHours: 0.221 RangeInDegrees: 0.27 Nmeas:   67 
Run: 1027 NEvents:  3059 RunDuration(hours): 0.454  # Run / Temp. data #  RangeInHours: 0.441 RangeInDegrees: 0.27 Nmeas:   73 
Run: 1028 NEvents: 25131 RunDuration(hours): 0.231  # Run / Temp. data #  RangeInHours: 0.223 RangeInDegrees: 0.27 Nmeas:   50 
Run: 1029 NEvents:  3050 RunDuration(hours): 0.443  # Run / Temp. data #  RangeInHours: 0.441 RangeInDegrees: 0.27 Nmeas:   91 
Run: 1030 NEvents:  5726 RunDuration(hours): 0.216  # Run / Temp. data #  RangeInHours: 0.215 RangeInDegrees: 0.27 Nmeas:   57 
Run: 1031 NEvents:   107 RunDuration(hours): 0.066  # Run / Temp. data #  RangeInHours: 0.055 RangeInDegrees: 0.27 Nmeas:   16 
Run: 1032 NEvents:  8705 RunDuration(hours): 0.194  # Run / Temp. data #  RangeInHours: 0.180 RangeInDegrees: 0.27 Nmeas:   59 
Run: 1033 NEvents: 29679 RunDuration(hours): 0.289  # Run / Temp. data #  RangeInHours: 0.283 RangeInDegrees: 0.27 Nmeas:   88 
Run: 1034 NEvents:  3033 RunDuration(hours): 0.441  # Run / Temp. data #  RangeInHours: 0.440 RangeInDegrees: 0.27 Nmeas:  137 
Run: 1035 NEvents: 25334 RunDuration(hours): 0.387  # Run / Temp. data #  RangeInHours: 0.384 RangeInDegrees: 0.27 Nmeas:   77 
Run: 1036 NEvents:  3464 RunDuration(hours): 0.511  # Run / Temp. data #  RangeInHours: 0.487 RangeInDegrees: 0.27 Nmeas:  152 
Run: 1037 NEvents: 29863 RunDuration(hours): 0.266  # Run / Temp. data #  RangeInHours: 0.265 RangeInDegrees: 0.27 Nmeas:   79 
Run: 1038 NEvents:  4631 RunDuration(hours): 0.698  # Run / Temp. data #  RangeInHours: 0.628 RangeInDegrees: 0.27 Nmeas:  118 
Run: 1039 NEvents: 31923 RunDuration(hours): 0.299  # Run / Temp. data #  RangeInHours: 0.123 RangeInDegrees: 0.14 Nmeas:   18 
Run: 1040 NEvents:  3923 RunDuration(hours): 0.586  # Run / Temp. data #  RangeInHours: 0.555 RangeInDegrees: 0.27 Nmeas:  107 
Run: 1041 NEvents: 54416 RunDuration(hours): 0.486  # Run / Temp. data #  RangeInHours: 0.228 RangeInDegrees: 0.14 Nmeas:   24 
Run: 1042 NEvents:  3097 RunDuration(hours): 0.439  # Run / Temp. data #  RangeInHours: 0.406 RangeInDegrees: 0.13 Nmeas:   62 
Run: 1043 NEvents: 30357 RunDuration(hours): 0.371  # Run / Temp. data #  RangeInHours: 0.315 RangeInDegrees: 0.27 Nmeas:   56 
Run: 1044 NEvents:  3122 RunDuration(hours): 0.459  # Run / Temp. data #  RangeInHours: 0.173 RangeInDegrees: 0.14 Nmeas:    6 
Run: 1045 NEvents: 26457 RunDuration(hours): 0.247  # Run / Temp. data #  RangeInHours: 0.214 RangeInDegrees: 0.13 Nmeas:   29 
Run: 1046 NEvents:  3105 RunDuration(hours): 0.481  # Run / Temp. data #  RangeInHours: 0.302 RangeInDegrees: 0.13 Nmeas:   63 
Run: 1047 NEvents: 25487 RunDuration(hours): 0.242  # Run / Temp. data #  RangeInHours: 0.199 RangeInDegrees: 0.26 Nmeas:   10 
Run: 1048 NEvents:  3024 RunDuration(hours): 0.434  # Run / Temp. data #  RangeInHours: 0.432 RangeInDegrees: 0.13 Nmeas:  122 
Run: 1049 NEvents: 25396 RunDuration(hours): 0.239  # Run / Temp. data #  RangeInHours: 0.235 RangeInDegrees: 0.13 Nmeas:   39 
Run: 1050 NEvents:  3024 RunDuration(hours): 0.435  # Run / Temp. data #  RangeInHours: 0.309 RangeInDegrees: 0.26 Nmeas:   14 
Run: 1051 NEvents: 25408 RunDuration(hours): 0.237  # Run / Temp. data #  RangeInHours: 0.007 RangeInDegrees: 0.13 Nmeas:    4 
Run: 1052 NEvents:  3027 RunDuration(hours): 0.410  # Run / Temp. data #  RangeInHours: 0.381 RangeInDegrees: 0.41 Nmeas:   69 
Run: 1053 NEvents: 25701 RunDuration(hours): 0.187  # Run / Temp. data #  RangeInHours: 0.183 RangeInDegrees: 0.14 Nmeas:   40 
Run: 1054 NEvents:  3028 RunDuration(hours): 0.376  # Run / Temp. data #  RangeInHours: 0.368 RangeInDegrees: 0.13 Nmeas:   88 
Run: 1056 NEvents: 25501 RunDuration(hours): 0.182  # Run / Temp. data #  RangeInHours: 0.171 RangeInDegrees: 0.28 Nmeas:   34 
Run: 1057 NEvents:  3128 RunDuration(hours): 0.448  # Run / Temp. data #  RangeInHours: 0.443 RangeInDegrees: 0.27 Nmeas:  129 
Run: 1059 NEvents: 25415 RunDuration(hours): 0.232  # Run / Temp. data #  RangeInHours: 0.226 RangeInDegrees: 0.13 Nmeas:   64 
Run: 1060 NEvents:  3023 RunDuration(hours): 0.467  # Run / Temp. data #  RangeInHours: 0.464 RangeInDegrees: 0.28 Nmeas:  131 
Run: 1061 NEvents: 25221 RunDuration(hours): 0.232  # Run / Temp. data #  RangeInHours: 0.223 RangeInDegrees: 0.27 Nmeas:   55 
Run: 1062 NEvents:  3039 RunDuration(hours): 0.452  # Run / Temp. data #  RangeInHours: 0.452 RangeInDegrees: 0.14 Nmeas:   90 
Run: 1086 NEvents:  3024 RunDuration(hours): 0.435  # Run / Temp. data #  RangeInHours: 0.427 RangeInDegrees: 0.27 Nmeas:  108 
Run: 1087 NEvents: 25256 RunDuration(hours): 0.236  # Run / Temp. data #  RangeInHours: 0.223 RangeInDegrees: 0.27 Nmeas:   38 
Run: 1088 NEvents:  3069 RunDuration(hours): 0.446  # Run / Temp. data #  RangeInHours: 0.413 RangeInDegrees: 0.27 Nmeas:   52 
Run: 1089 NEvents: 27219 RunDuration(hours): 0.249  # Run / Temp. data #  RangeInHours: 0.212 RangeInDegrees: 0.27 Nmeas:   25 
Run: 1090 NEvents:  3022 RunDuration(hours): 0.452  # Run / Temp. data #  RangeInHours: 0.409 RangeInDegrees: 0.27 Nmeas:   31 
Run: 1092 NEvents: 25946 RunDuration(hours): 0.309  # Run / Temp. data #  RangeInHours: 0.296 RangeInDegrees: 0.14 Nmeas:   26 
Run: 1094 NEvents:  3045 RunDuration(hours): 0.492  # Run / Temp. data #  RangeInHours: 0.381 RangeInDegrees: 0.27 Nmeas:   40 
Run: 1095 NEvents: 31347 RunDuration(hours): 0.299  # Run / Temp. data #  RangeInHours: 0.265 RangeInDegrees: 0.14 Nmeas:   29 
Run: 1096 NEvents:  3080 RunDuration(hours): 0.463  # Run / Temp. data #  RangeInHours: 0.447 RangeInDegrees: 0.14 Nmeas:   50 
Run: 1097 NEvents: 30677 RunDuration(hours): 0.288  # Run / Temp. data #  RangeInHours: 0.276 RangeInDegrees: 0.14 Nmeas:   23 
Run: 1098 NEvents:  3057 RunDuration(hours): 0.455  # Run / Temp. data #  RangeInHours: 0.359 RangeInDegrees: 0.14 Nmeas:   24 
Run: 1099 NEvents: 29857 RunDuration(hours): 0.307  # Run / Temp. data #  RangeInHours: 0.238 RangeInDegrees: 0.27 Nmeas:   12 
Run: 1100 NEvents:  3175 RunDuration(hours): 0.474  # Run / Temp. data #  RangeInHours: 0.445 RangeInDegrees: 0.27 Nmeas:  104 
Run: 1101 NEvents: 27187 RunDuration(hours): 0.275  # Run / Temp. data #  RangeInHours: 0.231 RangeInDegrees: 0.27 Nmeas:   55 
Run: 1102 NEvents:  3024 RunDuration(hours): 0.467  # Run / Temp. data #  RangeInHours: 0.432 RangeInDegrees: 0.27 Nmeas:   72 
Run: 1103 NEvents: 28033 RunDuration(hours): 0.267  # Run / Temp. data #  RangeInHours: 0.135 RangeInDegrees: 0.27 Nmeas:    4 
Run: 1105 NEvents:  4876 RunDuration(hours): 0.787  # Run / Temp. data #  RangeInHours: 0.783 RangeInDegrees: 0.40 Nmeas:  118 
Run: 1106 NEvents: 25386 RunDuration(hours): 0.263  # Run / Temp. data #  RangeInHours: 0.068 RangeInDegrees: 0.13 Nmeas:   10 
Run: 1107 NEvents:  3299 RunDuration(hours): 0.538  # Run / Temp. data #  RangeInHours: 0.530 RangeInDegrees: 0.26 Nmeas:   82 
Run: 1108 NEvents: 36196 RunDuration(hours): 0.354  # Run / Temp. data #  RangeInHours: 0.329 RangeInDegrees: 0.28 Nmeas:   47 
Run: 1109 NEvents:  3319 RunDuration(hours): 0.476  # Run / Temp. data #  RangeInHours: 0.393 RangeInDegrees: 0.27 Nmeas:   65 
Run: 1111 NEvents: 26807 RunDuration(hours): 0.282  # Run / Temp. data #  RangeInHours: 0.185 RangeInDegrees: 0.28 Nmeas:   49 
Run: 1112 NEvents:  2584 RunDuration(hours): 0.478  # Run / Temp. data #  RangeInHours: 0.472 RangeInDegrees: 0.28 Nmeas:   59 
Run: 1113 NEvents: 27998 RunDuration(hours): 0.294  # Run / Temp. data #  RangeInHours: 0.270 RangeInDegrees: 0.27 Nmeas:   18 
Run: 1114 NEvents:  3015 RunDuration(hours): 0.481  # Run / Temp. data #  RangeInHours: 0.422 RangeInDegrees: 0.28 Nmeas:   76 
Run: 1115 NEvents: 25266 RunDuration(hours): 0.244  # Run / Temp. data #  RangeInHours: 0.214 RangeInDegrees: 0.14 Nmeas:   12 
Run: 1116 NEvents:  3026 RunDuration(hours): 0.611  # Run / Temp. data #  RangeInHours: 0.593 RangeInDegrees: 0.27 Nmeas:  132 
Run: 1117 NEvents: 26259 RunDuration(hours): 0.220  # Run / Temp. data #  RangeInHours: 0.187 RangeInDegrees: 0.14 Nmeas:    8 
Run: 1118 NEvents:  3009 RunDuration(hours): 0.469  # Run / Temp. data #  RangeInHours: 0.366 RangeInDegrees: 0.14 Nmeas:   22 
Run: 1119 NEvents: 30291 RunDuration(hours): 0.306  # Run / Temp. data #  RangeInHours: 0.237 RangeInDegrees: 0.14 Nmeas:   23 
Run: 1121 NEvents:  3076 RunDuration(hours): 0.474  # Run / Temp. data #  RangeInHours: 0.466 RangeInDegrees: 0.14 Nmeas:   42 
Run: 1122 NEvents: 25838 RunDuration(hours): 0.248  # Run / Temp. data #  RangeInHours: 0.246 RangeInDegrees: 0.14 Nmeas:   13 
Run: 1123 NEvents:  3008 RunDuration(hours): 0.452  # Run / Temp. data #  RangeInHours: 0.326 RangeInDegrees: 0.27 Nmeas:    4 
Run: 1124 NEvents: 25222 RunDuration(hours): 0.255  # Run / Temp. data #  RangeInHours: 0.235 RangeInDegrees: 0.27 Nmeas:   16 
Run: 1125 NEvents:  3022 RunDuration(hours): 0.484  # Run / Temp. data #  RangeInHours: 0.441 RangeInDegrees: 0.14 Nmeas:  117 
Run: 1126 NEvents: 16993 RunDuration(hours): 0.400  # Run / Temp. data #  RangeInHours: 0.224 RangeInDegrees: 0.28 Nmeas:    6 
Run: 1128 NEvents: 10245 RunDuration(hours): 0.132  # Run / Temp. data #  RangeInHours: 0.128 RangeInDegrees: 0.13 Nmeas:   44 
Run: 1129 NEvents:  3358 RunDuration(hours): 0.591  # Run / Temp. data #  RangeInHours: 0.562 RangeInDegrees: 0.26 Nmeas:   30 
Run: 1131 NEvents: 25332 RunDuration(hours): 0.307  # Run / Temp. data #  RangeInHours: 0.269 RangeInDegrees: 0.27 Nmeas:   61 
Run: 1132 NEvents: 28460 RunDuration(hours): 0.371  # Run / Temp. data #  RangeInHours: 0.343 RangeInDegrees: 0.26 Nmeas:    8 
Run: 1133 NEvents:  3001 RunDuration(hours): 0.596  # Run / Temp. data #  RangeInHours: 0.573 RangeInDegrees: 0.27 Nmeas:  118 
Run: 1134 NEvents: 25624 RunDuration(hours): 0.310  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.14 Nmeas:    2 
Run: 1135 NEvents:  3005 RunDuration(hours): 0.598  # Run / Temp. data #  RangeInHours: 0.573 RangeInDegrees: 0.27 Nmeas:  142 
Run: 1138 NEvents: 25660 RunDuration(hours): 0.302  # Run / Temp. data #  RangeInHours: 0.271 RangeInDegrees: 0.14 Nmeas:   14 
Run: 1139 NEvents:  3004 RunDuration(hours): 0.599  # Run / Temp. data #  RangeInHours: 0.518 RangeInDegrees: 0.27 Nmeas:   65 
Run: 1141 NEvents: 25241 RunDuration(hours): 0.294  # Run / Temp. data #  RangeInHours: 0.286 RangeInDegrees: 0.27 Nmeas:   72 
Run: 1146 NEvents:  1287 RunDuration(hours): 0.539  # Run / Temp. data #  RangeInHours: 0.523 RangeInDegrees: 0.40 Nmeas:  100 
Run: 1147 NEvents:  3103 RunDuration(hours): 0.479  # Run / Temp. data #  RangeInHours: 0.450 RangeInDegrees: 0.27 Nmeas:   90 
Run: 1148 NEvents: 25587 RunDuration(hours): 0.347  # Run / Temp. data #  RangeInHours: 0.338 RangeInDegrees: 0.27 Nmeas:   74 
Run: 1149 NEvents:  3004 RunDuration(hours): 0.454  # Run / Temp. data #  RangeInHours: 0.447 RangeInDegrees: 0.14 Nmeas:   58 
Run: 1150 NEvents: 26426 RunDuration(hours): 0.365  # Run / Temp. data #  RangeInHours: 0.356 RangeInDegrees: 0.27 Nmeas:   61 
Run: 1151 NEvents:   826 RunDuration(hours): 0.282  # Run / Temp. data #  RangeInHours: 0.279 RangeInDegrees: 0.27 Nmeas:   64 
Run: 1153 NEvents:  2872 RunDuration(hours): 0.719  # Run / Temp. data #  RangeInHours: 0.625 RangeInDegrees: 0.27 Nmeas:   66 
Run: 1154 NEvents:   345 RunDuration(hours): 0.111  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1155 NEvents:  3056 RunDuration(hours): 0.217  # Run / Temp. data #  RangeInHours: 0.192 RangeInDegrees: 0.28 Nmeas:   25 
Run: 1157 NEvents: 22207 RunDuration(hours): 0.142  # Run / Temp. data #  RangeInHours: 0.032 RangeInDegrees: 0.13 Nmeas:    6 
Run: 1159 NEvents:  5103 RunDuration(hours): 0.175  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.13 Nmeas:    2 
Run: 1160 NEvents: 25159 RunDuration(hours): 0.318  # Run / Temp. data #  RangeInHours: 0.288 RangeInDegrees: 0.14 Nmeas:   51 
Run: 1161 NEvents:  5642 RunDuration(hours): 0.334  # Run / Temp. data #  RangeInHours: 0.130 RangeInDegrees: 0.26 Nmeas:    6 
Run: 1162 NEvents: 25069 RunDuration(hours): 0.295  # Run / Temp. data #  RangeInHours: 0.287 RangeInDegrees: 0.28 Nmeas:   51 
Run: 1163 NEvents: 10016 RunDuration(hours): 0.489  # Run / Temp. data #  RangeInHours: 0.477 RangeInDegrees: 0.40 Nmeas:   23 
Run: 1165 NEvents: 25512 RunDuration(hours): 0.291  # Run / Temp. data #  RangeInHours: 0.265 RangeInDegrees: 0.13 Nmeas:   45 
Run: 1166 NEvents: 10073 RunDuration(hours): 0.432  # Run / Temp. data #  RangeInHours: 0.425 RangeInDegrees: 0.27 Nmeas:  120 
Run: 1167 NEvents: 25188 RunDuration(hours): 0.283  # Run / Temp. data #  RangeInHours: 0.278 RangeInDegrees: 0.14 Nmeas:   65 
Run: 1168 NEvents: 10036 RunDuration(hours): 0.407  # Run / Temp. data #  RangeInHours: 0.393 RangeInDegrees: 0.13 Nmeas:   44 
Run: 1169 NEvents: 25408 RunDuration(hours): 0.262  # Run / Temp. data #  RangeInHours: 0.005 RangeInDegrees: 0.13 Nmeas:    4 
Run: 1170 NEvents: 10319 RunDuration(hours): 0.310  # Run / Temp. data #  RangeInHours: 0.301 RangeInDegrees: 0.14 Nmeas:   12 
Run: 1171 NEvents: 25210 RunDuration(hours): 0.264  # Run / Temp. data #  RangeInHours: 0.240 RangeInDegrees: 0.28 Nmeas:   18 
Run: 1172 NEvents: 11766 RunDuration(hours): 0.365  # Run / Temp. data #  RangeInHours: 0.347 RangeInDegrees: 0.13 Nmeas:   74 
Run: 1174 NEvents: 27407 RunDuration(hours): 0.273  # Run / Temp. data #  RangeInHours: 0.258 RangeInDegrees: 0.27 Nmeas:   30 
Run: 1175 NEvents: 29285 RunDuration(hours): 0.852  # Run / Temp. data #  RangeInHours: 0.742 RangeInDegrees: 0.28 Nmeas:   97 
Run: 1177 NEvents: 34505 RunDuration(hours): 0.357  # Run / Temp. data #  RangeInHours: 0.352 RangeInDegrees: 0.13 Nmeas:   66 
Run: 1178 NEvents: 26300 RunDuration(hours): 0.800  # Run / Temp. data #  RangeInHours: 0.594 RangeInDegrees: 0.27 Nmeas:   18 
Run: 1180 NEvents: 27202 RunDuration(hours): 0.286  # Run / Temp. data #  RangeInHours: 0.281 RangeInDegrees: 0.13 Nmeas:   80 
Run: 1181 NEvents: 11822 RunDuration(hours): 0.352  # Run / Temp. data #  RangeInHours: 0.299 RangeInDegrees: 0.13 Nmeas:   42 
Run: 1182 NEvents: 31085 RunDuration(hours): 0.309  # Run / Temp. data #  RangeInHours: 0.135 RangeInDegrees: 0.26 Nmeas:    8 
Run: 1183 NEvents: 13826 RunDuration(hours): 0.402  # Run / Temp. data #  RangeInHours: 0.395 RangeInDegrees: 0.13 Nmeas:   96 
Run: 1185 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1187 NEvents: 11527 RunDuration(hours): 0.429  # Run / Temp. data #  RangeInHours: 0.416 RangeInDegrees: 0.27 Nmeas:  108 
Run: 1188 NEvents: 34751 RunDuration(hours): 0.365  # Run / Temp. data #  RangeInHours: 0.334 RangeInDegrees: 0.27 Nmeas:   26 
Run: 1190 NEvents: 11924 RunDuration(hours): 0.437  # Run / Temp. data #  RangeInHours: 0.399 RangeInDegrees: 0.14 Nmeas:   35 
Run: 1191 NEvents: 35408 RunDuration(hours): 0.374  # Run / Temp. data #  RangeInHours: 0.361 RangeInDegrees: 0.27 Nmeas:   88 
Run: 1193 NEvents: 12591 RunDuration(hours): 0.455  # Run / Temp. data #  RangeInHours: 0.294 RangeInDegrees: 0.14 Nmeas:   67 
Run: 1194 NEvents: 34288 RunDuration(hours): 0.447  # Run / Temp. data #  RangeInHours: 0.416 RangeInDegrees: 0.27 Nmeas:   46 
Run: 1195 NEvents:   862 RunDuration(hours): 0.065  # Run / Temp. data #  RangeInHours: 0.018 RangeInDegrees: 0.14 Nmeas:    4 
Run: 1196 NEvents: 12175 RunDuration(hours): 0.653  # Run / Temp. data #  RangeInHours: 0.584 RangeInDegrees: 0.27 Nmeas:  165 
Run: 1197 NEvents: 28079 RunDuration(hours): 0.393  # Run / Temp. data #  RangeInHours: 0.388 RangeInDegrees: 0.27 Nmeas:  119 
Run: 1198 NEvents: 11119 RunDuration(hours): 0.522  # Run / Temp. data #  RangeInHours: 0.475 RangeInDegrees: 0.27 Nmeas:   84 
Run: 1199 NEvents: 39929 RunDuration(hours): 0.534  # Run / Temp. data #  RangeInHours: 0.244 RangeInDegrees: 0.27 Nmeas:    4 
Run: 1200 NEvents:  3403 RunDuration(hours): 0.353  # Run / Temp. data #  RangeInHours: 0.326 RangeInDegrees: 0.14 Nmeas:   68 
Run: 1201 NEvents: 10126 RunDuration(hours): 0.517  # Run / Temp. data #  RangeInHours: 0.308 RangeInDegrees: 0.14 Nmeas:   51 
Run: 1202 NEvents: 15848 RunDuration(hours): 0.709  # Run / Temp. data #  RangeInHours: 0.666 RangeInDegrees: 0.40 Nmeas:   97 
Run: 1203 NEvents: 37367 RunDuration(hours): 0.499  # Run / Temp. data #  RangeInHours: 0.383 RangeInDegrees: 0.40 Nmeas:   47 
Run: 1204 NEvents: 11030 RunDuration(hours): 0.475  # Run / Temp. data #  RangeInHours: 0.326 RangeInDegrees: 0.14 Nmeas:   44 
Run: 1205 NEvents: 31235 RunDuration(hours): 0.410  # Run / Temp. data #  RangeInHours: 0.402 RangeInDegrees: 0.40 Nmeas:   56 
Run: 1206 NEvents: 11660 RunDuration(hours): 0.162  # Run / Temp. data #  RangeInHours: 0.147 RangeInDegrees: 0.27 Nmeas:   19 
Run: 1207 NEvents: 10324 RunDuration(hours): 0.140  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.28 Nmeas:    8 
Run: 1208 NEvents: 10957 RunDuration(hours): 0.174  # Run / Temp. data #  RangeInHours: 0.166 RangeInDegrees: 0.27 Nmeas:   21 
Run: 1209 NEvents: 13643 RunDuration(hours): 0.173  # Run / Temp. data #  RangeInHours: 0.166 RangeInDegrees: 0.28 Nmeas:   52 
Run: 1210 NEvents: 10447 RunDuration(hours): 0.131  # Run / Temp. data #  RangeInHours: 0.078 RangeInDegrees: 0.14 Nmeas:    6 
Run: 1212 NEvents:  9583 RunDuration(hours): 0.121  # Run / Temp. data #  RangeInHours: 0.109 RangeInDegrees: 0.27 Nmeas:   19 
Run: 1219 NEvents:  3203 RunDuration(hours): 0.574  # Run / Temp. data #  RangeInHours: 0.568 RangeInDegrees: 0.28 Nmeas:  159 
Run: 1222 NEvents:  4016 RunDuration(hours): 0.341  # Run / Temp. data #  RangeInHours: 0.334 RangeInDegrees: 0.27 Nmeas:   84 
Run: 1227 NEvents:108465 RunDuration(hours): 0.954  # Run / Temp. data #  RangeInHours: 0.911 RangeInDegrees: 0.28 Nmeas:  135 
Run: 1228 NEvents:106767 RunDuration(hours): 0.923  # Run / Temp. data #  RangeInHours: 0.907 RangeInDegrees: 0.27 Nmeas:  156 
Run: 1229 NEvents:  8052 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.13 Nmeas:    2 
Run: 1230 NEvents:105105 RunDuration(hours): 0.707  # Run / Temp. data #  RangeInHours: 0.429 RangeInDegrees: 0.41 Nmeas:   52 
Run: 1232 NEvents:  1411 RunDuration(hours): 0.692  # Run / Temp. data #  RangeInHours: 0.413 RangeInDegrees: 0.14 Nmeas:   47 
Run: 1234 NEvents: 17736 RunDuration(hours): 0.640  # Run / Temp. data #  RangeInHours: 0.634 RangeInDegrees: 0.27 Nmeas:  164 
Run: 1235 NEvents:   751 RunDuration(hours): 0.613  # Run / Temp. data #  RangeInHours: 0.569 RangeInDegrees: 0.14 Nmeas:   90 
Run: 1236 NEvents: 11090 RunDuration(hours): 0.716  # Run / Temp. data #  RangeInHours: 0.667 RangeInDegrees: 0.28 Nmeas:   28 
Run: 1237 NEvents: 94737 RunDuration(hours): 2.243  # Run / Temp. data #  RangeInHours: 2.155 RangeInDegrees: 0.28 Nmeas:  329 
Run: 1239 NEvents: 51154 RunDuration(hours): 0.503  # Run / Temp. data #  RangeInHours: 0.368 RangeInDegrees: 0.27 Nmeas:   54 
Run: 1240 NEvents: 26343 RunDuration(hours): 0.368  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.13 Nmeas:    2 
Run: 1241 NEvents:   831 RunDuration(hours): 0.058  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1246 NEvents: 14936 RunDuration(hours): 0.186  # Run / Temp. data #  RangeInHours: 0.172 RangeInDegrees: 0.13 Nmeas:   37 
Run: 1252 NEvents:115552 RunDuration(hours): 1.846  # Run / Temp. data #  RangeInHours: 1.811 RangeInDegrees: 0.53 Nmeas:  244 
Run: 1257 NEvents:147505 RunDuration(hours): 2.265  # Run / Temp. data #  RangeInHours: 2.234 RangeInDegrees: 0.53 Nmeas:  406 
Run: 1258 NEvents: 50045 RunDuration(hours): 0.753  # Run / Temp. data #  RangeInHours: 0.696 RangeInDegrees: 0.28 Nmeas:   64 
Run: 1259 NEvents: 18420 RunDuration(hours): 0.409  # Run / Temp. data #  RangeInHours: 0.407 RangeInDegrees: 0.28 Nmeas:   55 
Run: 1261 NEvents:180299 RunDuration(hours): 1.816  # Run / Temp. data #  RangeInHours: 1.813 RangeInDegrees: 0.28 Nmeas:  394 
Run: 1262 NEvents: 25911 RunDuration(hours): 0.347  # Run / Temp. data #  RangeInHours: 0.322 RangeInDegrees: 0.28 Nmeas:   33 
Run: 1266 NEvents:158602 RunDuration(hours): 2.505  # Run / Temp. data #  RangeInHours: 2.488 RangeInDegrees: 0.27 Nmeas:  600 
Run: 1267 NEvents: 15589 RunDuration(hours): 0.220  # Run / Temp. data #  RangeInHours: 0.107 RangeInDegrees: 0.13 Nmeas:   33 
Run: 1269 NEvents: 46911 RunDuration(hours): 0.623  # Run / Temp. data #  RangeInHours: 0.584 RangeInDegrees: 0.53 Nmeas:   79 
Run: 1270 NEvents: 94751 RunDuration(hours): 1.829  # Run / Temp. data #  RangeInHours: 1.625 RangeInDegrees: 0.26 Nmeas:  140 
Run: 1271 NEvents:  2210 RunDuration(hours): 0.492  # Run / Temp. data #  RangeInHours: 0.484 RangeInDegrees: 0.28 Nmeas:  124 
Run: 1272 NEvents:  5469 RunDuration(hours): 0.752  # Run / Temp. data #  RangeInHours: 0.695 RangeInDegrees: 0.28 Nmeas:   60 
Run: 1274 NEvents: 14631 RunDuration(hours): 0.775  # Run / Temp. data #  RangeInHours: 0.760 RangeInDegrees: 0.40 Nmeas:   78 
Run: 1275 NEvents: 81797 RunDuration(hours): 2.662  # Run / Temp. data #  RangeInHours: 2.659 RangeInDegrees: 0.67 Nmeas:  330 
Run: 1276 NEvents: 44514 RunDuration(hours): 0.738  # Run / Temp. data #  RangeInHours: 0.547 RangeInDegrees: 0.27 Nmeas:  102 
Run: 1277 NEvents:  9638 RunDuration(hours): 0.071  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1279 NEvents: 77173 RunDuration(hours): 0.493  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1280 NEvents: 80723 RunDuration(hours): 0.630  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1282 NEvents:192866 RunDuration(hours): 2.085  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1284 NEvents: 37477 RunDuration(hours): 0.405  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1285 NEvents: 43349 RunDuration(hours): 0.291  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1286 NEvents: 45894 RunDuration(hours): 0.258  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1287 NEvents: 20042 RunDuration(hours): 0.138  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1288 NEvents: 73196 RunDuration(hours): 0.492  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1290 NEvents: 47199 RunDuration(hours): 0.295  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1291 NEvents: 60560 RunDuration(hours): 0.391  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1292 NEvents: 73891 RunDuration(hours): 0.499  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1293 NEvents: 31185 RunDuration(hours): 0.232  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1294 NEvents: 64023 RunDuration(hours): 0.437  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1295 NEvents: 39689 RunDuration(hours): 0.284  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1296 NEvents:200963 RunDuration(hours): 1.206  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1297 NEvents: 28617 RunDuration(hours): 0.206  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1299 NEvents:173623 RunDuration(hours): 1.143  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1300 NEvents:200231 RunDuration(hours): 1.422  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1301 NEvents:200574 RunDuration(hours): 1.420  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1302 NEvents: 74198 RunDuration(hours): 0.624  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1303 NEvents: 59552 RunDuration(hours): 0.684  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1305 NEvents: 67037 RunDuration(hours): 0.584  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1306 NEvents: 23299 RunDuration(hours): 0.184  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1307 NEvents:126993 RunDuration(hours): 1.204  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1308 NEvents: 27455 RunDuration(hours): 0.372  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1309 NEvents: 29006 RunDuration(hours): 0.423  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1310 NEvents: 30925 RunDuration(hours): 0.254  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1311 NEvents: 69399 RunDuration(hours): 1.910  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1313 NEvents: 26915 RunDuration(hours): 0.308  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1314 NEvents:107957 RunDuration(hours): 1.091  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1315 NEvents:  3077 RunDuration(hours): 0.043  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1317 NEvents:154191 RunDuration(hours): 1.334  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1318 NEvents: 35752 RunDuration(hours): 0.313  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1319 NEvents: 10924 RunDuration(hours): 0.318  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1320 NEvents:152805 RunDuration(hours): 1.434  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1321 NEvents: 63206 RunDuration(hours): 1.149  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1322 NEvents: 45282 RunDuration(hours): 0.636  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1326 NEvents:193502 RunDuration(hours): 1.933  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1327 NEvents:202316 RunDuration(hours): 2.102  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1328 NEvents:105343 RunDuration(hours): 1.003  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1329 NEvents:  4428 RunDuration(hours): 0.206  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1330 NEvents: 74064 RunDuration(hours): 0.729  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1331 NEvents: 34511 RunDuration(hours): 0.324  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1332 NEvents:103719 RunDuration(hours): 1.077  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1333 NEvents:101781 RunDuration(hours): 0.908  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1334 NEvents:119933 RunDuration(hours): 0.767  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1335 NEvents:100127 RunDuration(hours): 0.646  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1336 NEvents:100408 RunDuration(hours): 0.652  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1337 NEvents:103089 RunDuration(hours): 1.012  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1338 NEvents:106865 RunDuration(hours): 1.086  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1339 NEvents:170643 RunDuration(hours): 2.380  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1341 NEvents: 35461 RunDuration(hours): 0.636  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1342 NEvents: 30524 RunDuration(hours): 0.423  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1343 NEvents:  7146 RunDuration(hours): 0.084  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1344 NEvents: 28788 RunDuration(hours): 0.388  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1345 NEvents:102595 RunDuration(hours): 0.786  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1348 NEvents: 41805 RunDuration(hours): 0.316  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1349 NEvents:173261 RunDuration(hours): 1.177  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1350 NEvents:201103 RunDuration(hours): 1.900  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1351 NEvents: 29345 RunDuration(hours): 0.231  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1352 NEvents: 75251 RunDuration(hours): 0.557  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1353 NEvents:100798 RunDuration(hours): 0.718  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1354 NEvents: 62611 RunDuration(hours): 1.649  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1355 NEvents: 38380 RunDuration(hours): 1.085  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1356 NEvents: 91298 RunDuration(hours): 0.759  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1357 NEvents:104323 RunDuration(hours): 0.828  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1358 NEvents: 42980 RunDuration(hours): 0.463  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1360 NEvents:  3916 RunDuration(hours): 0.099  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1361 NEvents: 61209 RunDuration(hours): 0.347  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1362 NEvents:115612 RunDuration(hours): 0.669  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1363 NEvents:100629 RunDuration(hours): 0.530  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1364 NEvents: 43734 RunDuration(hours): 0.284  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1372 NEvents:     0 RunDuration(hours): 0.000  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1373 NEvents: 25906 RunDuration(hours): 0.275  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1376 NEvents: 25438 RunDuration(hours): 0.193  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1377 NEvents:100593 RunDuration(hours): 0.742  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1378 NEvents:100380 RunDuration(hours): 1.289  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1379 NEvents:102086 RunDuration(hours): 1.240  # Run / Temp. data #  RangeInHours: 0.265 RangeInDegrees: 0.14 Nmeas:   25 
Run: 1380 NEvents: 71985 RunDuration(hours): 0.583  # Run / Temp. data #  RangeInHours: 0.576 RangeInDegrees: 0.26 Nmeas:  126 
Run: 1381 NEvents: 29046 RunDuration(hours): 0.238  # Run / Temp. data #  RangeInHours: 0.199 RangeInDegrees: 0.27 Nmeas:   39 
Run: 1382 NEvents:109118 RunDuration(hours): 0.763  # Run / Temp. data #  RangeInHours: 0.755 RangeInDegrees: 0.27 Nmeas:  102 
Run: 1384 NEvents:137343 RunDuration(hours): 1.949  # Run / Temp. data #  RangeInHours: 1.946 RangeInDegrees: 0.26 Nmeas:  498 
Run: 1386 NEvents:200514 RunDuration(hours): 1.357  # Run / Temp. data #  RangeInHours: 1.333 RangeInDegrees: 0.27 Nmeas:  268 
Run: 1387 NEvents:203715 RunDuration(hours): 1.387  # Run / Temp. data #  RangeInHours: 1.383 RangeInDegrees: 0.40 Nmeas:  295 
Run: 1388 NEvents:201235 RunDuration(hours): 1.654  # Run / Temp. data #  RangeInHours: 1.650 RangeInDegrees: 0.40 Nmeas:  267 
Run: 1389 NEvents:200237 RunDuration(hours): 1.376  # Run / Temp. data #  RangeInHours: 1.356 RangeInDegrees: 0.27 Nmeas:  275 
Run: 1390 NEvents:  8352 RunDuration(hours): 0.071  # Run / Temp. data #  RangeInHours: 0.024 RangeInDegrees: 0.14 Nmeas:    6 
Run: 1391 NEvents: 94079 RunDuration(hours): 1.116  # Run / Temp. data #  RangeInHours: 0.975 RangeInDegrees: 0.26 Nmeas:  186 
Run: 1392 NEvents:100138 RunDuration(hours): 0.847  # Run / Temp. data #  RangeInHours: 0.832 RangeInDegrees: 0.26 Nmeas:  234 
Run: 1393 NEvents:100544 RunDuration(hours): 0.787  # Run / Temp. data #  RangeInHours: 0.641 RangeInDegrees: 0.52 Nmeas:   72 
Run: 1394 NEvents: 73909 RunDuration(hours): 0.982  # Run / Temp. data #  RangeInHours: 0.923 RangeInDegrees: 0.40 Nmeas:  123 
Run: 1396 NEvents: 25714 RunDuration(hours): 0.201  # Run / Temp. data #  RangeInHours: 0.170 RangeInDegrees: 0.27 Nmeas:    9 
Run: 1397 NEvents: 19978 RunDuration(hours): 0.201  # Run / Temp. data #  RangeInHours: 0.193 RangeInDegrees: 0.13 Nmeas:   47 
Run: 1398 NEvents: 81356 RunDuration(hours): 0.725  # Run / Temp. data #  RangeInHours: 0.722 RangeInDegrees: 0.80 Nmeas:  115 
Run: 1399 NEvents: 69802 RunDuration(hours): 0.748  # Run / Temp. data #  RangeInHours: 0.691 RangeInDegrees: 0.27 Nmeas:   91 
Run: 1400 NEvents: 20567 RunDuration(hours): 0.385  # Run / Temp. data #  RangeInHours: 0.352 RangeInDegrees: 0.27 Nmeas:   46 
Run: 1401 NEvents: 12751 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.125 RangeInDegrees: 0.27 Nmeas:   37 
Run: 1402 NEvents:   372 RunDuration(hours): 0.123  # Run / Temp. data #  RangeInHours: 0.121 RangeInDegrees: 0.14 Nmeas:   32 
Run: 1406 NEvents: 48419 RunDuration(hours): 0.678  # Run / Temp. data #  RangeInHours: 0.637 RangeInDegrees: 0.27 Nmeas:   25 
Run: 1406 NEvents: 48419 RunDuration(hours): 0.678  # Run / Temp. data #  RangeInHours: 0.637 RangeInDegrees: 0.27 Nmeas:   25 
Run: 1411 NEvents:  5619 RunDuration(hours): 0.061  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.13 Nmeas:    2 
Run: 1413 NEvents:100439 RunDuration(hours): 0.933  # Run / Temp. data #  RangeInHours: 0.911 RangeInDegrees: 0.79 Nmeas:  155 
Run: 1414 NEvents:  5908 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.27 Nmeas:   12 
Run: 1416 NEvents:  6554 RunDuration(hours): 0.170  # Run / Temp. data #  RangeInHours: 0.163 RangeInDegrees: 0.14 Nmeas:   35 
Run: 1417 NEvents: 85891 RunDuration(hours): 0.845  # Run / Temp. data #  RangeInHours: 0.787 RangeInDegrees: 0.28 Nmeas:  107 
Run: 1420 NEvents:101653 RunDuration(hours): 0.831  # Run / Temp. data #  RangeInHours: 0.700 RangeInDegrees: 0.28 Nmeas:  123 
Run: 1422 NEvents: 50277 RunDuration(hours): 1.174  # Run / Temp. data #  RangeInHours: 1.121 RangeInDegrees: 0.27 Nmeas:   79 
Run: 1423 NEvents:101652 RunDuration(hours): 0.754  # Run / Temp. data #  RangeInHours: 0.673 RangeInDegrees: 0.26 Nmeas:  121 
Run: 1424 NEvents:100375 RunDuration(hours): 1.473  # Run / Temp. data #  RangeInHours: 1.231 RangeInDegrees: 0.14 Nmeas:  223 
Run: 1428 NEvents: 86676 RunDuration(hours): 1.398  # Run / Temp. data #  RangeInHours: 1.058 RangeInDegrees: 0.28 Nmeas:  172 
Run: 1434 NEvents: 16916 RunDuration(hours): 0.241  # Run / Temp. data #  RangeInHours: 0.234 RangeInDegrees: 0.15 Nmeas:   61 
Run: 1436 NEvents:100341 RunDuration(hours): 0.748  # Run / Temp. data #  RangeInHours: 0.721 RangeInDegrees: 0.13 Nmeas:  130 
Run: 1437 NEvents:100081 RunDuration(hours): 0.757  # Run / Temp. data #  RangeInHours: 0.637 RangeInDegrees: 0.26 Nmeas:   50 
Run: 1438 NEvents:100256 RunDuration(hours): 1.279  # Run / Temp. data #  RangeInHours: 1.136 RangeInDegrees: 0.26 Nmeas:   16 
Run: 1439 NEvents: 50961 RunDuration(hours): 0.599  # Run / Temp. data #  RangeInHours: 0.455 RangeInDegrees: 0.14 Nmeas:   24 
Run: 1440 NEvents: 54984 RunDuration(hours): 0.485  # Run / Temp. data #  RangeInHours: 0.437 RangeInDegrees: 0.14 Nmeas:   56 
Run: 1441 NEvents: 27339 RunDuration(hours): 0.452  # Run / Temp. data #  RangeInHours: 0.045 RangeInDegrees: 0.14 Nmeas:   10 
Run: 1442 NEvents: 14324 RunDuration(hours): 0.132  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.27 Nmeas:   28 
Run: 1443 NEvents: 81966 RunDuration(hours): 0.619  # Run / Temp. data #  RangeInHours: 0.605 RangeInDegrees: 0.92 Nmeas:  152 
Run: 1445 NEvents: 38340 RunDuration(hours): 0.347  # Run / Temp. data #  RangeInHours: 0.339 RangeInDegrees: 0.27 Nmeas:   61 
Run: 1449 NEvents: 78695 RunDuration(hours): 0.781  # Run / Temp. data #  RangeInHours: 0.545 RangeInDegrees: 0.27 Nmeas:   96 
Run: 1451 NEvents: 59992 RunDuration(hours): 0.746  # Run / Temp. data #  RangeInHours: 0.256 RangeInDegrees: 0.27 Nmeas:   25 
Run: 1454 NEvents:106056 RunDuration(hours): 0.992  # Run / Temp. data #  RangeInHours: 0.988 RangeInDegrees: 0.27 Nmeas:   52 
Run: 1455 NEvents:113503 RunDuration(hours): 0.837  # Run / Temp. data #  RangeInHours: 0.767 RangeInDegrees: 0.13 Nmeas:  147 
Run: 1456 NEvents: 87970 RunDuration(hours): 0.706  # Run / Temp. data #  RangeInHours: 0.623 RangeInDegrees: 0.39 Nmeas:   60 
Run: 1458 NEvents: 12052 RunDuration(hours): 0.204  # Run / Temp. data #  RangeInHours: 0.181 RangeInDegrees: 0.13 Nmeas:   38 
Run: 1459 NEvents: 29715 RunDuration(hours): 0.250  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1460 NEvents: 75266 RunDuration(hours): 0.579  # Run / Temp. data #  RangeInHours: 0.429 RangeInDegrees: 0.27 Nmeas:   78 
Run: 1461 NEvents:100212 RunDuration(hours): 0.759  # Run / Temp. data #  RangeInHours: 0.603 RangeInDegrees: 0.14 Nmeas:  106 
Run: 1462 NEvents: 70607 RunDuration(hours): 1.373  # Run / Temp. data #  RangeInHours: 1.049 RangeInDegrees: 0.28 Nmeas:  198 
Run: 1463 NEvents: 44109 RunDuration(hours): 0.772  # Run / Temp. data #  RangeInHours: 0.728 RangeInDegrees: 0.28 Nmeas:   26 
Run: 1464 NEvents:100634 RunDuration(hours): 1.176  # Run / Temp. data #  RangeInHours: 1.163 RangeInDegrees: 0.14 Nmeas:   79 
Run: 1466 NEvents: 52155 RunDuration(hours): 0.450  # Run / Temp. data #  RangeInHours: 0.423 RangeInDegrees: 0.27 Nmeas:  109 
Run: 1467 NEvents: 48230 RunDuration(hours): 0.326  # Run / Temp. data #  RangeInHours: 0.319 RangeInDegrees: 0.27 Nmeas:   96 
Run: 1468 NEvents:101082 RunDuration(hours): 0.716  # Run / Temp. data #  RangeInHours: 0.691 RangeInDegrees: 0.28 Nmeas:  107 
Run: 1469 NEvents: 20765 RunDuration(hours): 0.205  # Run / Temp. data #  RangeInHours: 0.022 RangeInDegrees: 0.27 Nmeas:    4 
Run: 1470 NEvents: 99922 RunDuration(hours): 0.754  # Run / Temp. data #  RangeInHours: 0.650 RangeInDegrees: 0.26 Nmeas:   33 
Run: 1472 NEvents:  5834 RunDuration(hours): 0.172  # Run / Temp. data #  RangeInHours: 0.166 RangeInDegrees: 0.26 Nmeas:   47 
Run: 1473 NEvents:101870 RunDuration(hours): 0.771  # Run / Temp. data #  RangeInHours: 0.762 RangeInDegrees: 0.26 Nmeas:  180 
Run: 1474 NEvents:100353 RunDuration(hours): 0.734  # Run / Temp. data #  RangeInHours: 0.642 RangeInDegrees: 0.28 Nmeas:   32 
Run: 1475 NEvents:100084 RunDuration(hours): 1.053  # Run / Temp. data #  RangeInHours: 0.955 RangeInDegrees: 0.27 Nmeas:  144 
Run: 1476 NEvents:100349 RunDuration(hours): 1.023  # Run / Temp. data #  RangeInHours: 0.998 RangeInDegrees: 0.27 Nmeas:   70 
Run: 1490 NEvents:104651 RunDuration(hours): 0.881  # Run / Temp. data #  RangeInHours: 0.867 RangeInDegrees: 0.66 Nmeas:   95 
Run: 1492 NEvents:107056 RunDuration(hours): 0.855  # Run / Temp. data #  RangeInHours: 0.832 RangeInDegrees: 0.26 Nmeas:   96 
Run: 1493 NEvents:100111 RunDuration(hours): 1.879  # Run / Temp. data #  RangeInHours: 1.858 RangeInDegrees: 0.40 Nmeas:  222 
Run: 1494 NEvents:100577 RunDuration(hours): 0.811  # Run / Temp. data #  RangeInHours: 0.762 RangeInDegrees: 0.27 Nmeas:  155 
Run: 1495 NEvents:111226 RunDuration(hours): 1.411  # Run / Temp. data #  RangeInHours: 1.392 RangeInDegrees: 0.41 Nmeas:  205 
Run: 1499 NEvents:100380 RunDuration(hours): 0.888  # Run / Temp. data #  RangeInHours: 0.877 RangeInDegrees: 0.26 Nmeas:   48 
Run: 1502 NEvents: 66326 RunDuration(hours): 0.653  # Run / Temp. data #  RangeInHours: 0.568 RangeInDegrees: 0.27 Nmeas:   98 
Run: 1503 NEvents: 34813 RunDuration(hours): 0.300  # Run / Temp. data #  RangeInHours: 0.292 RangeInDegrees: 0.26 Nmeas:   75 
Run: 1504 NEvents: 44857 RunDuration(hours): 0.739  # Run / Temp. data #  RangeInHours: 0.648 RangeInDegrees: 0.26 Nmeas:   53 
Run: 1505 NEvents: 19783 RunDuration(hours): 0.298  # Run / Temp. data #  RangeInHours: 0.280 RangeInDegrees: 0.28 Nmeas:   24 
Run: 1506 NEvents: 36238 RunDuration(hours): 0.333  # Run / Temp. data #  RangeInHours: 0.323 RangeInDegrees: 0.27 Nmeas:   91 
Run: 1507 NEvents:100276 RunDuration(hours): 0.859  # Run / Temp. data #  RangeInHours: 0.840 RangeInDegrees: 0.40 Nmeas:   75 
Run: 1508 NEvents:100588 RunDuration(hours): 0.820  # Run / Temp. data #  RangeInHours: 0.809 RangeInDegrees: 0.26 Nmeas:  152 
Run: 1509 NEvents:101204 RunDuration(hours): 3.704  # Run / Temp. data #  RangeInHours: 3.672 RangeInDegrees: 0.53 Nmeas:  676 
Run: 1510 NEvents:100584 RunDuration(hours): 2.980  # Run / Temp. data #  RangeInHours: 2.968 RangeInDegrees: 0.40 Nmeas:  598 
Run: 1515 NEvents: 39480 RunDuration(hours): 1.697  # Run / Temp. data #  RangeInHours: 1.690 RangeInDegrees: 0.27 Nmeas:  400 
Run: 1516 NEvents: 20240 RunDuration(hours): 1.536  # Run / Temp. data #  RangeInHours: 1.491 RangeInDegrees: 0.40 Nmeas:  264 
Run: 1519 NEvents: 26328 RunDuration(hours): 0.968  # Run / Temp. data #  RangeInHours: 0.955 RangeInDegrees: 0.40 Nmeas:  208 
Run: 1520 NEvents: 51237 RunDuration(hours): 1.821  # Run / Temp. data #  RangeInHours: 1.820 RangeInDegrees: 0.79 Nmeas:  373 
Run: 1521 NEvents:  8208 RunDuration(hours): 0.348  # Run / Temp. data #  RangeInHours: 0.282 RangeInDegrees: 0.53 Nmeas:   32 
Run: 1522 NEvents: 10279 RunDuration(hours): 0.396  # Run / Temp. data #  RangeInHours: 0.386 RangeInDegrees: 0.27 Nmeas:   86 
Run: 1523 NEvents: 10329 RunDuration(hours): 0.885  # Run / Temp. data #  RangeInHours: 0.874 RangeInDegrees: 0.40 Nmeas:  181 
Run: 1527 NEvents: 10235 RunDuration(hours): 0.308  # Run / Temp. data #  RangeInHours: 0.298 RangeInDegrees: 0.26 Nmeas:   96 
Run: 1528 NEvents: 10016 RunDuration(hours): 0.252  # Run / Temp. data #  RangeInHours: 0.246 RangeInDegrees: 0.14 Nmeas:   30 
Run: 1529 NEvents: 10838 RunDuration(hours): 0.268  # Run / Temp. data #  RangeInHours: 0.264 RangeInDegrees: 0.27 Nmeas:   62 
Run: 1530 NEvents: 10747 RunDuration(hours): 0.123  # Run / Temp. data #  RangeInHours: 0.119 RangeInDegrees: 0.27 Nmeas:   40 
Run: 1532 NEvents: 10168 RunDuration(hours): 0.137  # Run / Temp. data #  RangeInHours: 0.136 RangeInDegrees: 0.27 Nmeas:   44 
Run: 1533 NEvents: 10620 RunDuration(hours): 0.127  # Run / Temp. data #  RangeInHours: 0.097 RangeInDegrees: 0.27 Nmeas:   24 
Run: 1536 NEvents: 10386 RunDuration(hours): 0.145  # Run / Temp. data #  RangeInHours: 0.136 RangeInDegrees: 0.14 Nmeas:   30 
Run: 1537 NEvents: 10366 RunDuration(hours): 0.129  # Run / Temp. data #  RangeInHours: 0.121 RangeInDegrees: 0.14 Nmeas:   29 
Run: 1538 NEvents: 10657 RunDuration(hours): 0.129  # Run / Temp. data #  RangeInHours: 0.117 RangeInDegrees: 0.26 Nmeas:   22 
Run: 1539 NEvents: 10222 RunDuration(hours): 0.125  # Run / Temp. data #  RangeInHours: 0.121 RangeInDegrees: 0.26 Nmeas:   36 
Run: 1540 NEvents: 10314 RunDuration(hours): 0.134  # Run / Temp. data #  RangeInHours: 0.125 RangeInDegrees: 0.26 Nmeas:   46 
Run: 1541 NEvents: 10314 RunDuration(hours): 0.126  # Run / Temp. data #  RangeInHours: 0.123 RangeInDegrees: 0.26 Nmeas:   52 
Run: 1542 NEvents: 10446 RunDuration(hours): 0.109  # Run / Temp. data #  RangeInHours: 0.107 RangeInDegrees: 0.26 Nmeas:   45 
Run: 1544 NEvents: 10030 RunDuration(hours): 0.106  # Run / Temp. data #  RangeInHours: 0.105 RangeInDegrees: 0.26 Nmeas:   43 
Run: 1545 NEvents: 10665 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.26 Nmeas:   44 
Run: 1546 NEvents: 10200 RunDuration(hours): 0.114  # Run / Temp. data #  RangeInHours: 0.092 RangeInDegrees: 0.26 Nmeas:   32 
Run: 1547 NEvents: 10514 RunDuration(hours): 0.118  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.26 Nmeas:   33 
Run: 1548 NEvents: 10339 RunDuration(hours): 0.117  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.26 Nmeas:   25 
Run: 1549 NEvents: 10371 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.108 RangeInDegrees: 0.14 Nmeas:   27 
Run: 1550 NEvents: 10233 RunDuration(hours): 0.121  # Run / Temp. data #  RangeInHours: 0.116 RangeInDegrees: 0.14 Nmeas:   31 
Run: 1551 NEvents: 10152 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.096 RangeInDegrees: 0.26 Nmeas:   19 
Run: 1552 NEvents: 10249 RunDuration(hours): 0.115  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.14 Nmeas:   30 
Run: 1553 NEvents: 10243 RunDuration(hours): 0.122  # Run / Temp. data #  RangeInHours: 0.108 RangeInDegrees: 0.14 Nmeas:   27 
Run: 1554 NEvents: 10311 RunDuration(hours): 0.118  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.14 Nmeas:   21 
Run: 1555 NEvents: 10228 RunDuration(hours): 0.117  # Run / Temp. data #  RangeInHours: 0.110 RangeInDegrees: 0.14 Nmeas:   25 
Run: 1556 NEvents: 10224 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.112 RangeInDegrees: 0.14 Nmeas:   26 
Run: 1557 NEvents: 10265 RunDuration(hours): 0.117  # Run / Temp. data #  RangeInHours: 0.112 RangeInDegrees: 0.14 Nmeas:   33 
Run: 1558 NEvents: 10031 RunDuration(hours): 0.119  # Run / Temp. data #  RangeInHours: 0.103 RangeInDegrees: 0.14 Nmeas:   20 
Run: 1559 NEvents:  6324 RunDuration(hours): 0.117  # Run / Temp. data #  RangeInHours: 0.098 RangeInDegrees: 0.27 Nmeas:   32 
Run: 1560 NEvents:  5053 RunDuration(hours): 0.067  # Run / Temp. data #  RangeInHours: 0.063 RangeInDegrees: 0.27 Nmeas:   12 
Run: 1562 NEvents: 10241 RunDuration(hours): 0.161  # Run / Temp. data #  RangeInHours: 0.154 RangeInDegrees: 0.27 Nmeas:   39 
Run: 1564 NEvents: 10043 RunDuration(hours): 0.157  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.27 Nmeas:   47 
Run: 1565 NEvents: 10274 RunDuration(hours): 0.168  # Run / Temp. data #  RangeInHours: 0.161 RangeInDegrees: 0.27 Nmeas:   37 
Run: 1566 NEvents:  1665 RunDuration(hours): 0.251  # Run / Temp. data #  RangeInHours: 0.244 RangeInDegrees: 0.27 Nmeas:   94 
Run: 1568 NEvents:  8762 RunDuration(hours): 0.140  # Run / Temp. data #  RangeInHours: 0.139 RangeInDegrees: 0.27 Nmeas:   54 
Run: 1569 NEvents: 10021 RunDuration(hours): 0.153  # Run / Temp. data #  RangeInHours: 0.143 RangeInDegrees: 0.27 Nmeas:   59 
Run: 1570 NEvents: 10138 RunDuration(hours): 0.171  # Run / Temp. data #  RangeInHours: 0.166 RangeInDegrees: 0.27 Nmeas:   61 
Run: 1571 NEvents: 10627 RunDuration(hours): 0.188  # Run / Temp. data #  RangeInHours: 0.179 RangeInDegrees: 0.27 Nmeas:   58 
Run: 1572 NEvents: 10036 RunDuration(hours): 0.152  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.27 Nmeas:   54 
Run: 1573 NEvents: 10281 RunDuration(hours): 0.153  # Run / Temp. data #  RangeInHours: 0.144 RangeInDegrees: 0.27 Nmeas:   45 
Run: 1574 NEvents: 11094 RunDuration(hours): 0.209  # Run / Temp. data #  RangeInHours: 0.204 RangeInDegrees: 0.27 Nmeas:   48 
Run: 1575 NEvents: 10208 RunDuration(hours): 0.161  # Run / Temp. data #  RangeInHours: 0.157 RangeInDegrees: 0.27 Nmeas:   41 
Run: 1576 NEvents: 10245 RunDuration(hours): 0.196  # Run / Temp. data #  RangeInHours: 0.186 RangeInDegrees: 0.27 Nmeas:   49 
Run: 1577 NEvents: 10237 RunDuration(hours): 0.165  # Run / Temp. data #  RangeInHours: 0.154 RangeInDegrees: 0.27 Nmeas:   34 
Run: 1580 NEvents: 11510 RunDuration(hours): 0.168  # Run / Temp. data #  RangeInHours: 0.148 RangeInDegrees: 0.14 Nmeas:   29 
Run: 1581 NEvents: 11004 RunDuration(hours): 0.159  # Run / Temp. data #  RangeInHours: 0.155 RangeInDegrees: 0.27 Nmeas:   42 
Run: 1582 NEvents: 10169 RunDuration(hours): 0.165  # Run / Temp. data #  RangeInHours: 0.151 RangeInDegrees: 0.14 Nmeas:   32 
Run: 1584 NEvents: 10464 RunDuration(hours): 0.149  # Run / Temp. data #  RangeInHours: 0.142 RangeInDegrees: 0.14 Nmeas:   34 
Run: 1585 NEvents: 10270 RunDuration(hours): 0.217  # Run / Temp. data #  RangeInHours: 0.211 RangeInDegrees: 0.27 Nmeas:   54 
Run: 1586 NEvents: 10084 RunDuration(hours): 0.148  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.27 Nmeas:   39 
Run: 1587 NEvents: 10571 RunDuration(hours): 0.163  # Run / Temp. data #  RangeInHours: 0.161 RangeInDegrees: 0.27 Nmeas:   42 
Run: 1588 NEvents: 10093 RunDuration(hours): 0.148  # Run / Temp. data #  RangeInHours: 0.142 RangeInDegrees: 0.27 Nmeas:   51 
Run: 1589 NEvents: 10072 RunDuration(hours): 0.153  # Run / Temp. data #  RangeInHours: 0.150 RangeInDegrees: 0.27 Nmeas:   62 
Run: 1590 NEvents: 10050 RunDuration(hours): 0.144  # Run / Temp. data #  RangeInHours: 0.143 RangeInDegrees: 0.27 Nmeas:   56 
Run: 1591 NEvents: 10219 RunDuration(hours): 0.147  # Run / Temp. data #  RangeInHours: 0.146 RangeInDegrees: 0.27 Nmeas:   67 
Run: 1562 NEvents: 10241 RunDuration(hours): 0.161  # Run / Temp. data #  RangeInHours: 0.154 RangeInDegrees: 0.27 Nmeas:   39 
Run: 1593 NEvents: 10518 RunDuration(hours): 0.158  # Run / Temp. data #  RangeInHours: 0.155 RangeInDegrees: 0.27 Nmeas:   58 
Run: 1594 NEvents: 10936 RunDuration(hours): 0.173  # Run / Temp. data #  RangeInHours: 0.171 RangeInDegrees: 0.27 Nmeas:   69 
Run: 1595 NEvents: 11166 RunDuration(hours): 0.164  # Run / Temp. data #  RangeInHours: 0.163 RangeInDegrees: 0.27 Nmeas:   67 
Run: 1596 NEvents: 12470 RunDuration(hours): 0.172  # Run / Temp. data #  RangeInHours: 0.168 RangeInDegrees: 0.27 Nmeas:   65 
Run: 1597 NEvents: 11433 RunDuration(hours): 0.156  # Run / Temp. data #  RangeInHours: 0.153 RangeInDegrees: 0.27 Nmeas:   66 
Run: 1598 NEvents: 10873 RunDuration(hours): 0.156  # Run / Temp. data #  RangeInHours: 0.152 RangeInDegrees: 0.27 Nmeas:   61 
Run: 1599 NEvents: 15678 RunDuration(hours): 0.211  # Run / Temp. data #  RangeInHours: 0.209 RangeInDegrees: 0.27 Nmeas:   80 
Run: 1612 NEvents:  5116 RunDuration(hours): 0.083  # Run / Temp. data #  RangeInHours: 0.034 RangeInDegrees: 0.26 Nmeas:    4 
Run: 1616 NEvents: 28954 RunDuration(hours): 0.220  # Run / Temp. data #  RangeInHours: 0.207 RangeInDegrees: 0.26 Nmeas:   48 
Run: 1617 NEvents: 64719 RunDuration(hours): 1.139  # Run / Temp. data #  RangeInHours: 1.109 RangeInDegrees: 0.40 Nmeas:  148 
Run: 1618 NEvents: 20970 RunDuration(hours): 0.161  # Run / Temp. data #  RangeInHours: 0.157 RangeInDegrees: 0.27 Nmeas:   33 
Run: 1619 NEvents:100298 RunDuration(hours): 1.376  # Run / Temp. data #  RangeInHours: 1.372 RangeInDegrees: 0.53 Nmeas:  304 
Run: 1621 NEvents:101082 RunDuration(hours): 1.363  # Run / Temp. data #  RangeInHours: 1.361 RangeInDegrees: 0.40 Nmeas:  213 
Run: 1624 NEvents: 97532 RunDuration(hours): 0.828  # Run / Temp. data #  RangeInHours: 0.823 RangeInDegrees: 0.40 Nmeas:  230 
Run: 1625 NEvents:102285 RunDuration(hours): 0.913  # Run / Temp. data #  RangeInHours: 0.910 RangeInDegrees: 0.40 Nmeas:  281 
Run: 1626 NEvents: 74963 RunDuration(hours): 0.586  # Run / Temp. data #  RangeInHours: 0.574 RangeInDegrees: 0.27 Nmeas:  160 
Run: 1628 NEvents: 31851 RunDuration(hours): 0.252  # Run / Temp. data #  RangeInHours: 0.251 RangeInDegrees: 0.27 Nmeas:   99 
Run: 1629 NEvents:101554 RunDuration(hours): 0.988  # Run / Temp. data #  RangeInHours: 0.981 RangeInDegrees: 0.40 Nmeas:  274 
Run: 1631 NEvents: 56688 RunDuration(hours): 0.488  # Run / Temp. data #  RangeInHours: 0.484 RangeInDegrees: 0.27 Nmeas:  186 
Run: 1632 NEvents: 50379 RunDuration(hours): 0.387  # Run / Temp. data #  RangeInHours: 0.381 RangeInDegrees: 0.14 Nmeas:  102 
Run: 1634 NEvents: 72259 RunDuration(hours): 0.508  # Run / Temp. data #  RangeInHours: 0.504 RangeInDegrees: 0.27 Nmeas:  150 
Run: 1635 NEvents: 32790 RunDuration(hours): 0.353  # Run / Temp. data #  RangeInHours: 0.350 RangeInDegrees: 0.27 Nmeas:  152 
Run: 1636 NEvents:104430 RunDuration(hours): 0.703  # Run / Temp. data #  RangeInHours: 0.700 RangeInDegrees: 0.27 Nmeas:  265 
Run: 1637 NEvents: 28746 RunDuration(hours): 0.324  # Run / Temp. data #  RangeInHours: 0.322 RangeInDegrees: 0.27 Nmeas:  118 
Run: 1638 NEvents: 26941 RunDuration(hours): 0.211  # Run / Temp. data #  RangeInHours: 0.204 RangeInDegrees: 0.27 Nmeas:   76 
Run: 1639 NEvents: 97542 RunDuration(hours): 0.660  # Run / Temp. data #  RangeInHours: 0.643 RangeInDegrees: 0.27 Nmeas:  188 
Run: 1641 NEvents:114782 RunDuration(hours): 0.794  # Run / Temp. data #  RangeInHours: 0.792 RangeInDegrees: 0.27 Nmeas:  215 
Run: 1643 NEvents:  7554 RunDuration(hours): 0.279  # Run / Temp. data #  RangeInHours: 0.276 RangeInDegrees: 0.27 Nmeas:   83 
Run: 1645 NEvents: 91988 RunDuration(hours): 1.099  # Run / Temp. data #  RangeInHours: 1.092 RangeInDegrees: 0.27 Nmeas:  398 
Run: 1646 NEvents: 56854 RunDuration(hours): 0.858  # Run / Temp. data #  RangeInHours: 0.858 RangeInDegrees: 0.27 Nmeas:  338 
Run: 1648 NEvents: 46481 RunDuration(hours): 1.078  # Run / Temp. data #  RangeInHours: 1.074 RangeInDegrees: 0.27 Nmeas:  391 
Run: 1653 NEvents: 39363 RunDuration(hours): 0.495  # Run / Temp. data #  RangeInHours: 0.476 RangeInDegrees: 0.27 Nmeas:  138 
Run: 1655 NEvents: 18714 RunDuration(hours): 0.467  # Run / Temp. data #  RangeInHours: 0.464 RangeInDegrees: 0.27 Nmeas:  136 
Run: 1658 NEvents: 42062 RunDuration(hours): 0.763  # Run / Temp. data #  RangeInHours: 0.758 RangeInDegrees: 0.27 Nmeas:  267 
Run: 1660 NEvents:100312 RunDuration(hours): 0.816  # Run / Temp. data #  RangeInHours: 0.812 RangeInDegrees: 0.40 Nmeas:  194 
Run: 1661 NEvents:101471 RunDuration(hours): 1.153  # Run / Temp. data #  RangeInHours: 1.152 RangeInDegrees: 0.40 Nmeas:  329 
Run: 1662 NEvents:100208 RunDuration(hours): 1.716  # Run / Temp. data #  RangeInHours: 1.704 RangeInDegrees: 0.52 Nmeas:  291 
Run: 1665 NEvents: 69235 RunDuration(hours): 1.426  # Run / Temp. data #  RangeInHours: 1.403 RangeInDegrees: 0.40 Nmeas:  180 
Run: 1666 NEvents: 32150 RunDuration(hours): 0.519  # Run / Temp. data #  RangeInHours: 0.511 RangeInDegrees: 0.27 Nmeas:  120 
Run: 1667 NEvents:100562 RunDuration(hours): 1.052  # Run / Temp. data #  RangeInHours: 1.049 RangeInDegrees: 0.27 Nmeas:  252 
Run: 1668 NEvents:101319 RunDuration(hours): 0.684  # Run / Temp. data #  RangeInHours: 0.679 RangeInDegrees: 0.27 Nmeas:  175 
Run: 1669 NEvents: 52340 RunDuration(hours): 0.728  # Run / Temp. data #  RangeInHours: 0.717 RangeInDegrees: 0.27 Nmeas:   58 
Run: 1670 NEvents: 50900 RunDuration(hours): 0.490  # Run / Temp. data #  RangeInHours: 0.486 RangeInDegrees: 0.27 Nmeas:   63 
Run: 1671 NEvents:100194 RunDuration(hours): 0.961  # Run / Temp. data #  RangeInHours: 0.961 RangeInDegrees: 0.27 Nmeas:  284 
Run: 1673 NEvents:100028 RunDuration(hours): 0.674  # Run / Temp. data #  RangeInHours: 0.664 RangeInDegrees: 0.27 Nmeas:  170 
Run: 1674 NEvents:101811 RunDuration(hours): 0.759  # Run / Temp. data #  RangeInHours: 0.749 RangeInDegrees: 0.27 Nmeas:   79 
Run: 1675 NEvents: 19991 RunDuration(hours): 0.304  # Run / Temp. data #  RangeInHours: 0.303 RangeInDegrees: 0.27 Nmeas:   73 
Run: 1676 NEvents: 81986 RunDuration(hours): 0.836  # Run / Temp. data #  RangeInHours: 0.830 RangeInDegrees: 0.40 Nmeas:  158 
Run: 1677 NEvents:102539 RunDuration(hours): 0.731  # Run / Temp. data #  RangeInHours: 0.726 RangeInDegrees: 0.40 Nmeas:  366 
Run: 1678 NEvents: 73316 RunDuration(hours): 0.891  # Run / Temp. data #  RangeInHours: 0.359 RangeInDegrees: 0.27 Nmeas:  138 
Run: 1680 NEvents: 16493 RunDuration(hours): 0.176  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1681 NEvents: 10572 RunDuration(hours): 0.085  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1682 NEvents:100473 RunDuration(hours): 1.213  # Run / Temp. data #  RangeInHours: 1.132 RangeInDegrees: 0.27 Nmeas:  313 
Run: 1683 NEvents:101318 RunDuration(hours): 0.651  # Run / Temp. data #  RangeInHours: 0.642 RangeInDegrees: 0.40 Nmeas:  164 
Run: 1685 NEvents: 27133 RunDuration(hours): 0.197  # Run / Temp. data #  RangeInHours: 0.186 RangeInDegrees: 0.27 Nmeas:   70 
Run: 1686 NEvents: 72988 RunDuration(hours): 0.388  # Run / Temp. data #  RangeInHours: 0.379 RangeInDegrees: 0.27 Nmeas:   95 
Run: 1687 NEvents:100395 RunDuration(hours): 0.529  # Run / Temp. data #  RangeInHours: 0.529 RangeInDegrees: 0.27 Nmeas:  171 
Run: 1688 NEvents:100303 RunDuration(hours): 0.741  # Run / Temp. data #  RangeInHours: 0.738 RangeInDegrees: 0.27 Nmeas:  232 
Run: 1689 NEvents: 53407 RunDuration(hours): 0.856  # Run / Temp. data #  RangeInHours: 0.850 RangeInDegrees: 0.27 Nmeas:  285 
Run: 1690 NEvents:   200 RunDuration(hours): 0.014  # Run / Temp. data #  RangeInHours: 0.013 RangeInDegrees: 0.14 Nmeas:    5 
Run: 1691 NEvents:   300 RunDuration(hours): 0.011  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    1 
Run: 1692 NEvents: 30119 RunDuration(hours): 0.639  # Run / Temp. data #  RangeInHours: 0.636 RangeInDegrees: 0.27 Nmeas:  238 
Run: 1693 NEvents: 20629 RunDuration(hours): 0.516  # Run / Temp. data #  RangeInHours: 0.514 RangeInDegrees: 0.27 Nmeas:  129 
Run: 1694 NEvents:100192 RunDuration(hours): 1.449  # Run / Temp. data #  RangeInHours: 1.432 RangeInDegrees: 0.53 Nmeas:  357 
Run: 1695 NEvents:  6483 RunDuration(hours): 0.423  # Run / Temp. data #  RangeInHours: 0.415 RangeInDegrees: 0.40 Nmeas:   59 
Run: 1696 NEvents: 50676 RunDuration(hours): 0.672  # Run / Temp. data #  RangeInHours: 0.661 RangeInDegrees: 0.66 Nmeas:   88 
Run: 1697 NEvents:   200 RunDuration(hours): 0.010  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1698 NEvents:   123 RunDuration(hours): 0.042  # Run / Temp. data #  RangeInHours: 0.002 RangeInDegrees: 0.14 Nmeas:    2 
Run: 1699 NEvents: 50894 RunDuration(hours): 0.718  # Run / Temp. data #  RangeInHours: 0.706 RangeInDegrees: 0.27 Nmeas:   73 
Run: 1700 NEvents: 20263 RunDuration(hours): 0.428  # Run / Temp. data #  RangeInHours: 0.378 RangeInDegrees: 0.26 Nmeas:   21 
Run: 1701 NEvents: 80294 RunDuration(hours): 1.219  # Run / Temp. data #  RangeInHours: 1.206 RangeInDegrees: 0.26 Nmeas:  200 
Run: 1702 NEvents:103318 RunDuration(hours): 0.610  # Run / Temp. data #  RangeInHours: 0.590 RangeInDegrees: 0.14 Nmeas:   33 
Run: 1704 NEvents: 48170 RunDuration(hours): 0.266  # Run / Temp. data #  RangeInHours: 0.218 RangeInDegrees: 0.14 Nmeas:   40 
Run: 1705 NEvents: 52292 RunDuration(hours): 0.276  # Run / Temp. data #  RangeInHours: 0.260 RangeInDegrees: 0.40 Nmeas:   44 
Run: 1716 NEvents:100127 RunDuration(hours): 1.027  # Run / Temp. data #  RangeInHours: 1.013 RangeInDegrees: 0.40 Nmeas:  201 
Run: 1719 NEvents:102420 RunDuration(hours): 0.739  # Run / Temp. data #  RangeInHours: 0.735 RangeInDegrees: 0.40 Nmeas:  152 
Run: 1720 NEvents: 73859 RunDuration(hours): 0.443  # Run / Temp. data #  RangeInHours: 0.439 RangeInDegrees: 0.27 Nmeas:   55 
Run: 1721 NEvents: 27735 RunDuration(hours): 0.158  # Run / Temp. data #  RangeInHours: 0.151 RangeInDegrees: 0.27 Nmeas:   35 
Run: 1722 NEvents:100140 RunDuration(hours): 0.650  # Run / Temp. data #  RangeInHours: 0.648 RangeInDegrees: 0.40 Nmeas:  156 
Run: 1723 NEvents:101267 RunDuration(hours): 0.570  # Run / Temp. data #  RangeInHours: 0.563 RangeInDegrees: 0.27 Nmeas:  122 
Run: 1724 NEvents: 37523 RunDuration(hours): 0.514  # Run / Temp. data #  RangeInHours: 0.507 RangeInDegrees: 0.27 Nmeas:  163 
Run: 1725 NEvents:  1772 RunDuration(hours): 0.116  # Run / Temp. data #  RangeInHours: 0.114 RangeInDegrees: 0.27 Nmeas:   39 
Run: 1726 NEvents: 49943 RunDuration(hours): 0.553  # Run / Temp. data #  RangeInHours: 0.552 RangeInDegrees: 0.27 Nmeas:  207 
Run: 1727 NEvents: 57963 RunDuration(hours): 1.015  # Run / Temp. data #  RangeInHours: 1.013 RangeInDegrees: 0.27 Nmeas:  308 
Run: 1728 NEvents: 23036 RunDuration(hours): 0.328  # Run / Temp. data #  RangeInHours: 0.157 RangeInDegrees: 0.27 Nmeas:   54 
Run: 1729 NEvents: 55008 RunDuration(hours): 0.490  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1734 NEvents: 38810 RunDuration(hours): 0.415  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1735 NEvents: 12936 RunDuration(hours): 0.174  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1736 NEvents:  2839 RunDuration(hours): 0.120  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1737 NEvents: 48668 RunDuration(hours): 0.458  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1738 NEvents: 50410 RunDuration(hours): 0.492  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1739 NEvents: 40178 RunDuration(hours): 0.360  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1741 NEvents: 50776 RunDuration(hours): 0.431  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1742 NEvents: 50111 RunDuration(hours): 0.412  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1743 NEvents: 50174 RunDuration(hours): 0.410  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
Run: 1744 NEvents: 50969 RunDuration(hours): 0.472  # Run / Temp. data #  RangeInHours: 0.000 RangeInDegrees: 0.00 Nmeas:    0 
*/
