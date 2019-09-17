<html>
<head>
<title>Random Numbers</title>
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

<form method='post' action='RandomNumbers.php'>
 
<h2>Random Numbers</h2> 
<ol id="toc">
  <li><a href="#section0">Internal random numbers</a></li>
  <li><a href="#section1">External random numbers</a></li>
  <li><a href="#section2">MIXMAX random numbers</a></li>
  <li><a href="#section3">The methods</a></li>
</ol>

 
This page describes the random-number generator in PYTHIA and 
how it can be replaced by an external one. 
 
<a name="section0"></a> 
<h3>Internal random numbers</h3> 
 
The <code>Rndm</code> class generates random numbers, using the 
Marsaglia-Zaman-Tsang algorithm [<a href="Bibliography.php#refMar90" target="page">Mar90</a>]. 
 
<p/> 
Random numbers <code>R</code> uniformly distributed in 
<code>0 &lt; R &lt; 1</code> are obtained with 
<pre> 
   Rndm::flat(); 
</pre> 
There are also methods to generate according to an exponential, to 
<i>x * exp(-x)</i>, to a Gaussian, or picked among a set of 
possibilities, which make use of <code>flat()</code>. 
 
<p/> 
If the random number generator is not initialized before, it will be 
so the first time it is asked to generate a random number, and then 
with the default seed, 19780503. This means that, by default, all runs 
will use identically the same random number sequence. This is 
convenient for debugging purposes, but dangerous if you intend to 
run several "identical" jobs to boost statistics. You can initialize, 
or reinitialize, with your own choice of seed with a 
<pre> 
   Rndm::init(seed); 
</pre> 
Here values <code>0 &lt; seed &lt; 900 000 000</code> gives so many 
different random number sequences, while <code>seed = 0</code> will call 
the <code>Stdlib time(0)</code> function to provide a "random" 
<code>seed</code>, and <code>seed &lt; 0</code> will revert back to 
the default <code>seed</code>. 
 
<p/> 
The <code>Pythia</code> class defines <?php $filepath = $_GET["filepath"];
echo "<a href='RandomNumberSeed.php?filepath=".$filepath."' target='page'>";?>a 
flag and a mode</a>, that allows the <code>seed</code> to be set in 
the <code>Pythia::init</code> call. That would be the standard way for a 
user to pick the random number sequence in a run. 
 
<a name="section1"></a> 
<h3>External random numbers</h3> 
 
<code>RndmEngine</code> is a base class for the external handling of 
random-number generation. The user-written derived class is called 
if a pointer to it has been handed in with the 
<code>pythia.rndmEnginePtr()</code> method. Since the default 
Marsaglia-Zaman-Tsang algorithm is quite good, chances are that any 
replacement would be a step down, but this may still be required by 
consistency with other program elements in big experimental frameworks. 
 
<p/> 
There is only one pure virtual method in <code>RndmEngine</code>, to 
generate one random number flat in the range between 0 and 1: 
<pre> 
  virtual double flat() = 0; 
</pre> 
Note that methods for initialization are not provided in the base 
class, in part since input parameters may be specific to the generator 
used, in part since initialization can as well be taken care of 
externally to the <code>Pythia</code> code. 
 
<p/> 
An example illustrating how to run with an external random number 
generator is provided in <code>main23.cc</code>. 
 
<a name="section2"></a> 
<h3>MIXMAX random numbers</h3> 
 
The MIXMAX class of random number generators utilizes 
matrix-recursion based on Anosov-Kolmogorov C-K systems, with the 
ability to create a large number of statistically independent 
sequences of random numbers based on different initial seeds. This is 
particularly advantageous in creating statistically independent 
samples when running a large number of parallel jobs, each with a 
different initial seed. In the plugin 
header <code>Pythia8Plugins/MixMax.h</code> an implementation of a 
MIXMAX random number generator is provided [<a href="Bibliography.php#refSav91" target="page">Sav91</a>,<a href="Bibliography.php#refSav15" target="page">Sav15</a>], 
courtesy of Konstantin Savvidy, as well as a PYTHIA interface through 
the <code>MixMaxRndm</code> class. 
 
In this implementation a dimensionality of 17 is used, as this has 
been found to provide faster access to large numbers of independent 
sequences. A timing comparison between the external MIXMAX random 
number generator, and the default internal PYTHIA random number 
generator is provided in the example <code>main23.cc</code>. The 
MIXMAX random number generator is found to be comparable in speed to 
the default generator. The primary methods of 
the <code>MixMaxRndm</code> class are given here. 
 
<a name="anchor1"></a>
<p/><strong> MixMaxRndm::MixMaxRndm(int seed0, int seed1, int seed2, int seed3) &nbsp;</strong> <br/>
for the given four 32-bit seed numbers. The sequence of numbers 
produced from this set of seeds is guaranteed not to collide with 
another if at least one bit of the four seeds is different, and, less 
than <i>10^100</i> random numbers are thrown. 
   
 
<a name="section3"></a> 
<h3>The methods</h3> 
 
We here collect a more complete and formal overview of 
the <code>Rndm</code> class methods. 
 
<a name="anchor2"></a>
<p/><strong> Rndm::Rndm() &nbsp;</strong> <br/>
construct a random number generator, but does not initialize it. 
   
 
<a name="anchor3"></a>
<p/><strong> Rndm::Rndm(int seed) &nbsp;</strong> <br/>
construct a random number generator, and initialize it for the 
given seed number. 
   
 
<a name="anchor4"></a>
<p/><strong> bool Rndm::rndmEnginePtr( RndmEngine* rndmPtr) &nbsp;</strong> <br/>
pass in pointer for external random number generation. 
   
 
<a name="anchor5"></a>
<p/><strong> void Rndm::init(int seed = 0) &nbsp;</strong> <br/>
initialize, or reinitialize, the random number generator for the given 
seed number. Not necessary if the seed was already set in the constructor. 
   
 
<a name="anchor6"></a>
<p/><strong> double Rndm::flat() &nbsp;</strong> <br/>
generate next random number uniformly between 0 and 1. 
   
 
<a name="anchor7"></a>
<p/><strong> double Rndm::exp() &nbsp;</strong> <br/>
generate random numbers according to <i>exp(-x)</i>. 
   
 
<a name="anchor8"></a>
<p/><strong> double Rndm::xexp() &nbsp;</strong> <br/>
generate random numbers according to <i>x exp(-x)</i>. 
   
 
<a name="anchor9"></a>
<p/><strong> double Rndm::gauss() &nbsp;</strong> <br/>
generate random numbers according to <i>exp(-x^2/2)</i>. 
   
 
<a name="anchor10"></a>
<p/><strong> pair&lt;double, double&gt; Rndm::gauss2() &nbsp;</strong> <br/>
generate a pair of random numbers according to 
<i>exp( -(x^2 + y^2) / 2)</i>. Is faster than two calls 
to <code>gauss()</code>. 
   
 
<a name="anchor11"></a>
<p/><strong> int Rndm::pick(const vector&lt;double&gt;&amp; prob) &nbsp;</strong> <br/>
pick one option among vector of (positive) probabilities. 
   
 
<a name="anchor12"></a>
<p/><strong> bool Rndm::dumpState(string fileName) &nbsp;</strong> <br/>
save the current state of the random number generator to a binary 
file. This involves two integers and 100 double-precision numbers. 
Intended for debug purposes. Note that binary files may be 
platform-dependent and thus not transportable. 
   
 
<a name="anchor13"></a>
<p/><strong> bool Rndm::readState(string fileName) &nbsp;</strong> <br/>
set the state of the random number generator by reading in a binary 
file saved by the above command. Comments as above. 
   
 
<a name="anchor14"></a>
<p/><strong> virtual double RndmEngine::flat() &nbsp;</strong> <br/>
if you want to construct an external random number generator 
(or generator interface) then you must implement this method 
in your class derived from the <code>RndmEningen</code> base class, 
to give a random number between 0 and 1. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
