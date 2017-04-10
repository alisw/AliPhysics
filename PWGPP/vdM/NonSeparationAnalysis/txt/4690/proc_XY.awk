#
# sort LHCdata_4690.txt | awk -F ';' -f proc_XY.awk
#
function unixsec(d1,d2) {
    dd1=d1;
    gsub("\\."," ",dd1);
    split(dd1,b," ");
    split(d2,a,"\\.");
    gsub(":"," ",a[1]);
    ts=sprintf("%s %s %s %s",b[3],b[2],b[1],a[1]);
    gsub("[.:]"," ",ts);
    return sprintf("%.3f", mktime(ts)+a[2]/1e9)
}

BEGIN {
    counter=0;
    flag=0;
}
/Nominal_Separation/ {
    sep = $6;
}

/Acquisition_Flag;;1/ {
    t0=unixsec($1,$2)+3600;
    flag=1
}
/Acquisition_Flag;;0/ {
    flag=0
    t1=unixsec($1,$2)+3600;
    printf("%.3f %.3f %12.8f # %s %s\n", t0,t1,1*sep, $1,$2);
}

/Step_Progress;;60/ {
    doPrint=1;
}

/Step_Progress;;0/ {
    if (doPrint) {
	t1=unixsec($1,$2)+3600;
#	printf("%.3f %.3f %s\n", t0,t1,$0);
    }
    doPrint=0;
}

