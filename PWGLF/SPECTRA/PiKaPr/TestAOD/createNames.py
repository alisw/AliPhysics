#!/usr/bin/python

def main():

    # print header
    output = list()
    output.append ("namespace AliSpectraNameSpace\n")
    output.append ("{\n")
    output.append ("   const char * kHistName[] =\n")
    output.append ("   {\n")

    # print histogram names
    ifile  = open("Histograms.h", "rb")
    for line in ifile:
        lineNoWS = line.strip()
        if(not lineNoWS.startswith("k")): # skip everything which is not an entry in the enum
            continue
        if("=" in lineNoWS): # skip histogram type delimeters
            continue
        col=line.split(",")
        
        output.append("     \"h"+col[0].strip()[1:]+"\",\n");

    output.append ("   };\n")

    # write file
    outfile = open("HistogramNames.h", "w")
    outfile.write("//This file was generated automatically, please do not edit!!\n\n");
    outfile.writelines(output)
    outfile.close()


## def skipLines(lineNoWS):
##     beginningsToSkip = ["//", "{", "namespace", "enum"]
##     for entry in beginningsToSkip:
##         if lineNoWS.startswith(entry):
##             return 1 
##     return 0
    
    

#######################################################################
if __name__ == "__main__":
    main()
