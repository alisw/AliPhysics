#!/usr/bin/env python
debug = True # enable trace

def trace(x):
	global debug
	if debug: print(x)

trace("loading...")

from itertools import combinations, combinations_with_replacement
from glob import glob
from math import *
import operator
from os.path import basename
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.linear_model
import sklearn.feature_selection
import datetime

def prec_from_pathname(path):
	if   '2k' in path: return 0.002
	elif '5k' in path: return 0.005
	else: raise AssertionError('Unknown field strengh: %s' % path)

# ['x', 'y', 'z', 'xx', 'xy', 'xz', 'yy', ...]
def combinatrial_vars(vars_str='xyz', length=3):
	term_list = []
	for l in range(length):
		term_list.extend([''.join(v) for v in combinations_with_replacement(list(vars_str), 1 + l)])
	return term_list

# product :: a#* => [a] -> a
def product(xs):
	return reduce(operator.mul, xs, 1) # foldl in Haskell

# (XYZ, "xx") -> XX
def term(dataframe, vars_str):
	return product(map(lambda x: dataframe[x], list(vars_str)))

# (f(X), Y) -> (max deviation, max%, avg dev, avg%)
def deviation_stat(fX, Y, prec=0.005):
	dev = np.abs(fX - Y)
	(max_dev, avg_dev) = (dev.max(axis=0), dev.mean(axis=0))
	(max_pct, avg_pct) = (max_dev / prec * 100, avg_dev / prec * 100)
	return (max_dev, max_pct, avg_dev, avg_pct)

# IO Df
def load_samples(path, cylindrical_axis=True, absolute_axis=True, genvars=[]):
	sample_cols = ['x', 'y', 'z', 'Bx', 'By', 'Bz']
	df = pd.read_csv(path, sep=' ', names=sample_cols)

	if cylindrical_axis:
		df['r']    = np.sqrt(df.x**2 + df.y**2)
		df['p']    = np.arctan2(df.y, df.x)
		df['Bt']   = np.sqrt(df.Bx**2 + df.By**2)
		df['Bpsi'] = np.arctan2(df.By, df.Bx) - np.arctan2(df.y, df.x)
		df['Br']   = df.Bt * np.cos(df.Bpsi)
		df['Bp']   = df.Bt * np.sin(df.Bpsi)

	if absolute_axis:
		df['X']    = np.abs(df.x)
		df['Y']    = np.abs(df.y)
		df['Z']    = np.abs(df.z)

	for var in genvars:
		df[var] = term(df, var)

	return df

def choose(vars, df1, df2):
	X1 = df1.loc[:, vars].as_matrix()
	X2 = df2.loc[:, vars].as_matrix()
	return (X1, X2)

# IO ()
def run_analysis_for_all_fields():
	sample_set = glob("dat_z22/*2k*.sample.dat")
	test_set   = glob("dat_z22/*2k*.test.dat")
	#print(sample_set, test_set)
	assert(len(sample_set) == len(test_set) and len(sample_set) > 0)

	result = pd.DataFrame()
	for i, sample_file in enumerate(sample_set):
		trace("run_analysis('%s', '%s')" % (sample_file, test_set[i]))
		df = run_analysis(sample_file, test_set[i])
		result = result.append(df, ignore_index=True)
	write_header(result)

def run_analysis(sample_file = 'dat_z22/tpc2k-z0-q2.sample.dat',
                 test_file   = 'dat_z22/tpc2k-z0-q2.test.dat'):
	global precision, df, test, lr, la, xvars_full, xvars, yvars, X, Y, Xtest, Ytest, ana_result
	precision = prec_from_pathname(sample_file)
	assert(precision == prec_from_pathname(test_file))

	xvars_full = combinatrial_vars('xyz', 3)[3:] # variables except x, y, z upto 3 dims

	trace("reading training samples... " + sample_file)
	df = load_samples(sample_file, genvars=xvars_full)

	trace("reading test samples..." + test_file)
	test = load_samples(test_file, genvars=xvars_full)

	trace("linear regression fit...")
	lr = sklearn.linear_model.LinearRegression()
	#ri = sklearn.linear_model.RidgeCV()
	#la = sklearn.linear_model.LassoCV()
	fs = sklearn.feature_selection.RFE(lr, 1, verbose=0)

	#xvars = ['x','y','z','xx','yy','zz','xy','yz','xz','xzz','yzz']
	#xvars = ["xx", "yy", "zz", 'x', 'y', 'z', 'xzz', 'yzz']
	#xvars = ['xxxr', 'xrrX', 'zzrX', 'p', 'xyrr', 'xzzr', 'xrrY', 'xzrX', 'xxxz', 'xzzr']
	#xvars=['x', 'xzz', 'xyz', 'yz', 'yy', 'zz', 'xy', 'xx', 'z', 'y', 'xz', 'yzz']
	yvars = ['Bx', 'By', 'Bz']
	#yvars = ['Bz']
	(Y, Ytest) = choose(yvars, df, test)
	#(Y, Ytest) = (df['Bz'], test['Bz'])

	xvars = combinatrial_vars('xyz', 3) # use all terms upto 3rd power
	(X, Xtest) = choose(xvars, df, test)

	for y in yvars:
		fs.fit(X, df[y])
		res = pd.DataFrame({ "term": xvars, "rank": fs.ranking_ })
		trace(y)
		trace(res.sort_values(by = "rank"))
	#xvars=list(res.sort_values(by="rank")[:26]['term'])

	lr.fit(X, Y)
	trace(', '.join(yvars) + " = 1 + " + ' + '.join(xvars))
	test_dev = deviation_stat(lr.predict(Xtest), Ytest, prec=precision)
	#for i in range(len(yvars)):
	#	arr = [lr.intercept_[i]] + lr.coef_[i]
	#	arr = [ str(x) for x in arr ]
	#	print(yvars[i] + " = { " + ', '.join(arr) + " }")
	#	print("deviation stat [test]: max %.2e (%.1f%%) avg %.2e (%.1f%%)" %
	#		   ( test_dev[0][i], test_dev[1][i], test_dev[2][i], test_dev[3][i] ))

	(sample_score, test_score) = (lr.score(X, Y), lr.score(Xtest, Ytest))
	trace("linear regression R^2 [train data]: %.8f" % sample_score)
	trace("linear regression R^2 [test data] : %.8f" % test_score)

	return pd.DataFrame(
	 	{ "xvars": [xvars],
	 	  "yvars": [yvars],
	 	  "max_dev": [test_dev[0]],
	 	  "max%": [test_dev[1]],
	 	  "avg_dev": [test_dev[2]],
	 	  "avg%": [test_dev[3]],
	 	  "sample_score": [sample_score],
	 	  "score": [test_score],
	 	  "coeffs": [lr.coef_],
	 	  "intercept": [lr.intercept_],
	 	  "sample_file": [sample_file],
	 	  "test_file": [test_file],
	 	  "precision": [precision],
	 	  "volume_id": [volume_id_from_path(sample_file)]
	 	})

def volume_id_from_path(path):
	return basename(path)\
	       .replace('.sample.dat', '')\
	       .replace('-', '_')
	       
def get_location_by_volume_id(id):
	if 'its' in id: r_bin = 0
	if 'tpc' in id: r_bin = 1
	if 'tof' in id: r_bin = 2
	if 'tofext' in id: r_bin = 3
	if 'cal' in id: r_bin = 4

	z_bin = int(id.split('_')[1][1:]) # "tofext2k_z0_q4" -> 0

	if 'q1' in id: quadrant = 0
	if 'q2' in id: quadrant = 1
	if 'q3' in id: quadrant = 2
	if 'q4' in id: quadrant = 3

	return r_bin, z_bin, quadrant

def write_header(result):
	#result.to_csv("magfield_params.csv")
	#result.to_html("magfield_params.html")
	print("# This file was generated from sysid.py at " + str(datetime.datetime.today()))
	print("# " + ', '.join(result.iloc[0].yvars) + " = 1 + " + ' + '.join(result.iloc[0].xvars))
	print("# barrel r: 0 < its < 80 < tpc < 250 < tof < 400 < tofext < 423 < cal < 500")
	print("# barrel z: -550 < z < 550")
	print("# phi: 0 < q1 < 0.5pi < q2 < pi < q3 < 1.5pi < q4 < 2pi")
	print("# header: Rbin Zbin Quadrant Nval_per_compoment(=20)")
	print("# data: Nval_per_compoment x floats")
	#print("# R^2: coefficient of determination in multiple linear regression. [0,1]")
	
	print("")

	for index, row in result.iterrows():
		#print("// ** %s  -  R^2 %s" % (row.volume_id, row.score))
		print("#" + row.volume_id)
		r_bin, z_bin, quadrant = get_location_by_volume_id(row.volume_id)
		print("%s %s %s 20" % (r_bin, z_bin, quadrant))
		for i, yvar in enumerate(row.yvars):
			name = row.volume_id #+ '_' + yvar.lower()
			print("# precision: tgt %.2e max %.2e (%.1f%%) avg %.2e (%.1f%%)" %
				  (row['precision'], row['max_dev'][i], row['max%'][i], row['avg_dev'][i], row['avg%'][i]))
			coef = [row['intercept'][i]] + list(row['coeffs'][i])
			arr = [ "%.5e" % x for x in coef ]
			body = ' '.join(arr)
			#decl = "const double[] %s = { %s };\n" % (name, body)
			#print(decl)
			print(body)
		print("")

#write_header(run_analysis())
run_analysis_for_all_fields()

#for i in range(10):
#	for xvars in combinations(xvars_full, i+1):
	#(X, Xtest) = choose(xvars, df, test)
	#lr.fit(X, Y)
	#ri.fit(X, Y)
	#la.fit(X, Y)
	#fs.fit(X, Y)
	#print xvars
	#(sample_score, test_score) = (lr.score(X, Y), lr.score(Xtest, Ytest))
	#print("linear R^2[sample] %.8f" % sample_score)
	#print("linear R^2[test]   %.8f" % test_score)
	#(sample_score2, test_score2) = (la.score(X, Y), la.score(Xtest, Ytest))
	#print("lasso  R^2[sample] %.8f" % sample_score2)
	#print("lasso  R^2[test]   %.8f" % test_score2)
	#print(la.coef_)

#for i in range(len(yvars)):
#	print(yvars[i])
#	print(pd.DataFrame({"Name": xvars, "Params": lr.coef_[i]}).sort_values(by='Params'))
#	print("+ %e" % lr.intercept_[i])

#sample_dev = deviation_stat(lr.predict(X), Y, prec=precision)
#test_dev   = deviation_stat(lr.predict(Xtest), Ytest, prec=precision)
#test_dev2  = deviation_stat(la.predict(Xtest), Ytest, prec=precision)
#print("[sample] max %.2e (%.1f%%) avg %.2e (%.1f%%)" % sample_dev)
#print("[test]   max %.2e (%.1f%%) avg %.2e (%.1f%%)" % test_dev  )
#print("lasso [test]   max %.2e (%.1f%%) avg %.2e (%.1f%%)" % test_dev2  )
