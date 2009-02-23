#!/bin/bash
rm analysis_{mc,esd,esd_raw}.root correction_map.root
root.exe -b -q 'run.C(2,"'$1'",-1,0,0,2,1,"")'
aliroot correct.C
