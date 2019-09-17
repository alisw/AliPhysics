<html>
<head>
<title>Random-Number Seed</title>
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

<form method='post' action='RandomNumberSeed.php'>
 
<h2>Random-Number Seed</h2> 
 
The seed of the random number generator can be set as follows: 
 
<br/><br/><strong>Random:setSeed</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Indicates whether a user-set seed should be used every time the 
<code>Pythia::init</code> routine is called. If off, the random number 
generator is initialized with its default seed at the beginning 
of the run, and never again. If on, each new <code>Pythia::init</code> 
call (should several be made in the same run) results in the random 
number being re-initialized, thereby possibly starting over with the 
same sequence, if you do not watch out. 
   
 
<br/><br/><table><tr><td><strong>Random:seed  </td><td></td><td> <input type="text" name="2" value="-1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>maximum = 900000000</code>)</td></tr></table>
The seed to be used, if <code>setSeed</code> is on.<br/> 
A negative value gives the default seed,<br/> 
a value 0 gives a random seed based on the time, and<br/> 
a value between 1 and 900,000,000 a unique different random number 
sequence. 
   
 
<p/> 
For  on random numbers see <?php $filepath = $_GET["filepath"];
echo "<a href='RandomNumbers.php?filepath=".$filepath."' target='page'>";?>here</a>. 
This includes methods to save and restore the state of the generator, 
and some preprogrammed methods to generate non-uniform random numbers. 
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "off")
{
$data = "Random:setSeed = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "-1")
{
$data = "Random:seed = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
