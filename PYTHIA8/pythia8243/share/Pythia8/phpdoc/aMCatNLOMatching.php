<html>
<head>
<title>aMC@NLO Matching</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='aMCatNLOMatching.php'>
 
<h2>aMC@NLO Matching</h2> 
 
The aMC@NLO package [<a href="Bibliography.php#refFri02" target="page">Fri02</a>] attempts to automate the MC@NLO matching 
procedure [<a href="Bibliography.php#refFri02" target="page">Fri02</a>]. MC@NLO interprets the parton shower as NLO 
subtraction method, and removes unwanted parton-shower contributions by 
extending the subtraction scheme used to generate NLO fixed-order results. 
Upon showering, an NLO accurate prediction for inclusive observables is 
achieved. This makes MC@NLO a convenient NLO+PS matching scheme. A 
consistent extended subtraction in the NLO fixed-order result makes 
analytic knowledge of the shower emission probability necessary. Once 
this is known, interfacing the (parton-shower specific) NLO calculation 
with the shower is straightforward. 
 
<p/> 
To allow for a fast, automatic generation of shower subtractions that are 
used in the fixed-order calculation, Pythia allows to generate emissions 
with a "global" recoil scheme, in which the recoil of an emission is shared 
among all final state particles. When using aMC@NLO, this global recoil must 
be switched on. Please see 
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>Timelike Showers</a> for details and options. 
 
<p/> 
A minimal set of settings necessary for a consistent treatment of aMC@NLO 
inputs is 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>SpaceShower:pTmaxMatch = 1</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>SpaceShower:pTmaxFudge = 1.</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>TimeShower:pTmaxMatch = 1</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>TimeShower:pTmaxFudge = 1.</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; 
<code>SpaceShower:MEcorrections = off</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; 
<code>TimeShower:MEcorrections = off</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>TimeShower:globalRecoil = on</code> 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; 
<code>TimeShower:weightGluonToQuark = 1</code> 
<p/> 
 
and further (process-specific) settings related global recoils. 
 
<p/> 
Some comments are in order. The settings 
<code>SpaceShower:pTmaxMatch = 1</code>, 
<code>SpaceShower:pTmaxFudge = 1.</code>, 
<code>TimeShower:pTmaxMatch = 1</code>, 
<code>TimeShower:pTmaxFudge = 1.</code> are included to ensure that the 
correct parton shower starting scale (i.e. the scale set when generating the 
subtractions in MC@NLO) is used within Pythia. Note that the last three 
options are default in Pythia8, and that the first option differs from the 
default only if the input state does not contain final state partons. 
 
<p/> 
Matrix element corrections to the parton shower splitting kernels have to be 
switched off by <code>SpaceShower:MEcorrections = off</code> and 
 <code>TimeShower:MEcorrections = off</code> . This is necessary because 
the matrix element corrections are not suitable for showers in the 
global recoil scheme, and because it is not viable to include process-specific 
shower probabilities in an automatic framewrok like aMC@NLO. 
 
<p/> 
<code>TimeShower:globalRecoil = on</code> is necessary. Formally, it 
is allowed to switch back to a local recoil treatment beyond the first 
proposed emission of any of the hard scattering partons in Born-type events. 
Pythia offers three choices at which stage the global recoil is dropped in 
favour of a local strategy. It is necessary to supplement the setting 
<code>TimeShower:globalRecoil = on</code> by additional settings specifying 
which global recoil strategy should be used. As these choices are up to the 
user, please consult <?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>Timelike Showers</a> 
for details on these options. 
<p/> 
 
<p/> 
Finally, <code>TimeShower:weightGluonToQuark = 1</code> is not default any 
longer, but was it at the time the subtractions were first implemented, 
and so is required  for consistency until further notice. 
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
