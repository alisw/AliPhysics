import sys
import os
import subprocess
workdirectory="/alice/cern.ch/user/b/blim/Xi1530Extractor"
gridworkdirectory="Xi1530Extractor"
BASE_PATH=os.path.dirname(os.path.abspath(__file__))
inputfiledirectory=BASE_PATH + "/"
savepath=BASE_PATH+"/data/"
commandlist = ["submit", "download", "check", "local"]
Mode = 5 #0 for submit, 1 for download, 2 for file check, 3 for local job with missing jobs
if not os.path.exists(savepath):
	os.makedirs(savepath)

if(len(sys.argv) < 2):
	print("Usage: python grid.py submit (or download or check or local)")
	exit()

if str(sys.argv[1]) in commandlist:
	Mode = commandlist.index(str(sys.argv[1]))
else:
	print("Choose from: "+", ".join(commandlist))
	exit()

# Configuration #
NormRangeStudy=1
BinCoutStudy=1
fitstudy=1
default=1
cutvar=1
fitvar=1


# Grid Input files #
inputdata = "AnalysisResults_LHC15fi_16deghijklop_17cefgijklmor_SYS_light_fixed.root"
inputRsnMCdata = "AnalysisResults_Xi1530LHC18c6b_RsnMC_SYS_fixed.root"
inputGenMCdata = "AnalysisResults_Xi1530LHC16_GenMC_final.root"

Inputfiles=[
	"extract.jdl",
	"extract.sh",
	#inputdata,
	#inputRsnMCdata,
	#inputGenMCdata,
	"DrawXi1530.C"
]

# Global variables #
totaljobs = 0
totalskip = 1
complete = 1
missingjobs = []

#Multi_bins = ["0 100", "0 10", "10 30", "30 50", "50 70", "70 100", "0 0.01", "0.01 0.1"]
Multi_bins = ["0 100", "0 10", "10 30", "30 50", "50 70", "70 100"]
#Multi_bins = ["30 100"]
#Multi_bins = ["0.01 0.1", "0 0.01"]
CutSysematic_bins = ["DefaultOption", "TPCNsigmaXi1530PionLoose", "TPCNsigmaXi1530PionTight", "TPCNsigmaXiLoose", "TPCNsigmaXiTight", "Xi1530PionZVertexLoose", "Xi1530PionZVertexTight", "DCADistLambdaDaughtersLoose", "DCADistLambdaDaughtersTight", "DCADistXiDaughtersLoose", "DCADistXiDaughtersTight", "DCADistLambdaPVLoose", "DCADistLambdaPVTight", "V0CosineOfPointingAngleLoose", "V0CosineOfPointingAngleTight", "CascadeCosineOfPointingAngleLoose", "CascadeCosineOfPointingAngleTight", "XiMassWindowLoose", "XiMassWindowTight", "XiTrackCut"]
FitSysematic_bins = ["LikeSignBkg", "BkgFit"]
Normvariations=["NormVarm", "NormVarp", "NormVarLm", "NormVarLp", "NormVarRm", "NormVarRp"]
Fitvariations=["FitVarLm", "FitVarLp", "FitVarRm", "FitVarRp", "FitVarBothm", "FitVarBothp"]
Bincountvariations =["BinCountLm", "BinCountRp", "BinCountBothp", "BinCountLp", "BinCountRm", "BinCountBothm"]
Bincountvariations_caution=["BinCountLp", "BinCountRm", "BinCountBothm"]


## Some functions ##
def if_exist(filename):
	exists = os.path.isfile(savepath+filename)
	if exists:
		return 1
	else:
		return 0

def download(args):
	global totaljobs
	global totalskip
	totaljobs += 1
	file = "AnalysisResults_Extracted_" + str(args[0]) + "_Multi_" + "%0.2f"%float(args[1]) + "-" + "%0.2f"%float(args[2]) + "_" + str(args[3]) + str(args[4]) + ".root"
	
	if(if_exist(file)):
		command = "file found, skip download: " + file
		rproc = subprocess.Popen('echo ' + command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		totalskip += 1
		return rproc
	
	folderpath = str(args[0]) + str(args[1]) + str(args[2]) + str(args[3]) + str(args[4])
	filepath = "alien://" + workdirectory + "/out/" + folderpath + "/" + file
	localpath = savepath + file
	command = """echo 'TGrid::Connect(\"alien://\");TFile::Cp(\"""" + filepath + """\", \"""" + localpath + """\");' | root -b -l | grep -v Welcome | grep -v Trying """
	#print(command)
	rproc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	outtext = rproc.communicate()[1]
	if "No more images to try" in str(outtext):
		missingjobs.append(args)
	return rproc

def submit(args):
	global totaljobs
	totaljobs += 1
	options = str(args[0]) + " " + str(args[1]) + " " + str(args[2]) + " " + str(args[3]) + " " + str(args[4])
	command = """echo 'TGrid::Connect(\"alien://\");gGrid->Cd(\"""" + gridworkdirectory + """\");TGridResult *res = gGrid->Command(\"submit extract.jdl """ + options + """\");cout << res->GetKey(0,"jobId") << endl;' | root -b -l | grep -v Welcome | grep -v Trying """
	#command = "alien_submit alien:" + workdirectory + "/extract.jdl " + options #old way
	#print(command)
	rproc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	return rproc

def filecheck(args):
	global totaljobs
	global totalskip
	totaljobs += 1
	file = "AnalysisResults_Extracted_" + str(args[0]) + "_Multi_" + "%0.2f"%float(args[1]) + "-" + "%0.2f"%float(args[2]) + "_" + str(args[3]) + str(args[4]) + ".root"
	
	if(if_exist(file)):
		command = file+ ": available" 
		rproc = subprocess.Popen('echo ' + command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		totalskip += 1
		return rproc
	else:
		missingjobs.append(args)
		command = file+ ": Missing ------------------------------------------------------" 
		rproc = subprocess.Popen('echo ' + command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		return rproc

def localjob(args):
	global totaljobs
	global totalskip
	totaljobs += 1
	file = "AnalysisResults_Extracted_" + str(args[0]) + "_Multi_" + "%0.2f"%float(args[1]) + "-" + "%0.2f"%float(args[2]) + "_" + str(args[3]) + str(args[4]) + ".root"
	
	if(if_exist(file)):
		command = file+ ": available" 
		rproc = subprocess.Popen('echo ' + command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		totalskip += 1
		return rproc
	else:
		command = """echo '-l -b -q DrawXi1530.C\(""" + str(args[0]) + ","+ str(args[1]) + ","+ str(args[2]) + """,\\\""""+ str(args[3]) + """\\\","""+ str(args[4])+ "\)' >> " + inputfiledirectory + "parallel.txt"
		rproc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
		return rproc

def doit(args):
	if(Mode == 0):
		return submit(args)
	if(Mode == 1):
		return download(args)
	if(Mode == 2):
		return filecheck(args)
	if(Mode == 3):
		return localjob(args)

##

####
if Mode is 3: #Initialize
	subprocess.call("echo '' > " + inputfiledirectory + "parallel.txt", shell=True, stderr=subprocess.PIPE)

submitjobs = []
if Mode is 0:
	for pushfile in Inputfiles:
		command = """echo 'TGrid::Connect(\"alien://\");TFile::Cp(\"""" + inputfiledirectory + pushfile + """\", \"alien://""" + workdirectory + """/""" + pushfile + """\");' | root -b -l | grep -v Welcome | grep -v Trying """
		submitjobs.append(subprocess.Popen(command, shell=True))
	for proc in submitjobs:
		output = proc.communicate()
		if output[0]:
			print(output)
	

procs = []
procs_out = []
while(complete):
	for bins in Multi_bins:
		if(default):
			fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], "Default", 1]
			proc = doit(fileinfo)
			fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], "MCcheck", 1]
			proc = doit(fileinfo)
			procs.append(proc)
		if(NormRangeStudy):
			for normvariation in Normvariations:
				for optionnumber in range(1,5):
					fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], normvariation, optionnumber]
					proc = doit(fileinfo)
					procs.append(proc)
			for proc in procs:
				proc.wait()
		if(fitstudy):
			for fitvariation in Fitvariations:
				for optionnumber in range(1,6):
					fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], fitvariation, optionnumber]
					proc = doit(fileinfo)
					procs.append(proc)
			for proc in procs:
				proc.wait()
		if(BinCoutStudy):
			for bcvariation in Bincountvariations:
				for optionnumber in range(1,10):
					if(bcvariation in Bincountvariations_caution):
						if(optionnumber > 4):
							continue
						else:
							fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], bcvariation, optionnumber]
							proc = doit(fileinfo)
							procs.append(proc)
					else:
						fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], bcvariation, optionnumber]
						proc = doit(fileinfo)
						procs.append(proc)
			for proc in procs:
				proc.wait()
		if(cutvar):
			for sys in range(1,len(CutSysematic_bins)):
				fileinfo = [str(sys+1), bins.split(" ")[0], bins.split(" ")[1], CutSysematic_bins[sys], 1]
				proc = doit(fileinfo)
				procs.append(proc)
			for proc in procs:
				proc.wait()
		if(fitvar):
			for sys in FitSysematic_bins:
				fileinfo = [1, bins.split(" ")[0], bins.split(" ")[1], sys, 1]
				proc = doit(fileinfo)
				procs.append(proc)
		for proc in procs:
				proc.wait()
	# waiting for the finishing
	complete = totaljobs - totalskip
	print("total " + str(totaljobs) + "jobs / skip: " + str(totalskip)+ "\n")
	break
print("total "+str(len(missingjobs))+" are missing, try to resubmit them")
print(missingjobs)

jobIds = []
proc = []
if Mode is 1 or 2:
	answer = raw_input("Resubmit? [yN]: ")
	if(answer == "y"):
		takearest = 0
		for argument in missingjobs:
			takearest += 1
			proc = submit(argument)
			procs.append(proc)
		for proc in procs:
			output = proc.communicate()[0]
			jobIds.append(int(s) for output in str.split() if output.isdigit())
		
		print("rebumitted jobs: ")	
		for jobid in jobIds:
			print(jobid)
	if(answer == "N"):
		print("Ok, done.")
	else:
		print("?")

if Mode is 3:
	subprocess.call("sed '1d' " + inputfiledirectory + "parallel.txt > tmpfile; mv tmpfile " + inputfiledirectory + "parallel.txt", shell=True, stderr=subprocess.PIPE)
	subprocess.call("cd " + inputfiledirectory, shell=True, stderr=subprocess.PIPE)
	answer = raw_input("Start local job? [yN]: ")
	if(answer == "y"):
		subprocess.call("cat parallel.txt | xargs -P1 -L 1 root", shell=True, stderr=subprocess.PIPE)
	if(answer == "N"):
		print("Ok, done.")
	else:
		print("?")
