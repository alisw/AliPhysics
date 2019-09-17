<html>
<head>
<title>ProMC files</title>
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

<form method='post' action='ProMCFiles.php'>
 
<h2>ProMC Files</h2> 
 
<a href="http://atlaswww.hep.anl.gov/asc/promc/">ProMC</a> 
[<a href="Bibliography.php#refChe14" target="page">Che14</a>] is a library 
for the storage of Monte Carlo event records, or other data, in a very 
compact binary form. It provides routines for fast input to and output 
from these compact data files. It uses "varints" as a way to store and 
compress an integer using a variable number of bytes, based on Google's 
platform- and language-neutral Protocol Buffers. Real numbers are converted 
to integers, e.g. by specifying a smallest unit of energy, momentum and 
length. Thereby a low-energy particle can be represented by a smaller 
number of bytes. 
 
<p/> 
The current PYTHIA linking and interface is based on ProMC version 1.5; 
earlier version will not do. Once you have to installed the ProMC library, 
you should configure PYTHIA with 
<pre> 
  ./configure --with-promc=/path/to/ProMC 
</pre> 
and recompile the PYTHIA library. As usual more fine-grained options 
are available to set paths to binaries, libraries and headers separately, 
if need be. 
 
<p/> 
The <code>examples/main46.cc</code> sample program illustrates how to 
write PYTHIA events onto a ProMC file. 
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
