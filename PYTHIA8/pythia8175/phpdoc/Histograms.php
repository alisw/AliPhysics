<html>
<head>
<title>Histograms</title>
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

<form method='post' action='Histograms.php'>
    
<h2>Histograms</h2>

The <code>Hist</code> class gives a simple implementation of 
one-dimensional histograms, useful for quick-and-dirty testing, 
without the need to link to more sophisticated packages. 
For this reason it is used in many of the
<?php $filepath = $_GET["filepath"];
echo "<a href='SampleMainPrograms.php?filepath=".$filepath."' target='page'>";?>sample main programs</a>
found in the <code>examples</code> subdirectory.

<h3>Basic principles</h3>

We here provide a simple overview of what is involved.
As a first step you need to declare a histogram, with name,
title, number of bins and <i>x</i> range (from, to).
<pre>
   Hist ZpT( "Z0 pT spectrum", 100, 0., 100.);
</pre>
Alternatively you can first declare it and later define it:
<pre>
   Hist ZpT;
   ZpT.book( "Z0 pT spectrum", 100, 0., 100.);
</pre>

Once declared, its contents can be added by repeated calls to 
<code>fill</code>,
<pre>
   ZpT.fill( 22.7, 1.); 
</pre>
where the first argument is the <i>x</i> value and the second the 
weight. Since the weight defaults to 1 the last argument could have 
been omitted in this case.   

<p/>
A set of overloaded operators have been defined, so that histograms 
can be added, subtracted, divided or multiplied by each other. Then the
contents are modified accordingly bin by bin. Thus the relative
deviation between two histograms <code>data</code> and 
<code>theory</code> can be found as
<pre>
  diff = (data - theory) / (data + theory);
</pre>
assuming that <code>diff</code>, <code>data</code> and <code>theory</code>
have been booked with the same number of bins and <i>x</i> range. That 
responsibility rests on the user; some checks are made for compatibility, 
but not enough to catch all possible mistakes. 

<p/>
Also overloaded operations with double real numbers are available. 
Again these four operations are defined bin by bin, i.e. the 
corresponding amount is added to, subtracted from, multiplied by or
divided by each bin. The double number can come before or after the
histograms, with obvious results. Thus the inverse of a histogram 
<code>result</code> is given by <code>1. / result</code>. 
The two kind of operations can be combined, e.g.
<pre>
  allpT = ZpT + 2. * WpT
</pre>
Finally, also the <code>+=, -+, *=, /=</code> are overloaded, with 
the right-hand side being either a histogram or a real number. 

<h3>Output format</h3>

<p/>
A histogram can be printed by making use of the overloaded &lt;&lt; 
operator, e.g.:
<pre>
   cout &lt;&lt; ZpT;
</pre>
The printout format is inspired by the old HBOOK one. To understand 
how to read this format, consider the simplified example 
<pre>
                                    
        3.50*10^ 2  9                     
        3.00*10^ 2  X   7               
        2.50*10^ 2  X  1X               
        2.00*10^ 2  X6 XX                
        1.50*10^ 2  XX5XX                 
        1.00*10^ 2  XXXXX                
        0.50*10^ 2  XXXXX        

          Contents 
            *10^ 2  31122
            *10^ 1  47208
            *10^ 0  79373

          Low edge  -- 
            *10^ 1  10001 
            *10^ 0  05050
</pre>
The key feature is that the <code>Contents</code> and  
<code>Low edge</code> have to be read vertically. For instance, 
the first bin has the contents 
<code>3 * 10^2 + 4 * 10^1 + 7 * 10^0 = 347</code>. Correspondingly,
the other bins have contents 179, 123, 207 and 283. The first bin 
stretches from <code>-(1 * 10^1 + 0 * 10^0) = -10</code> to the 
beginning of the second bin, at <code>-(0 * 10^1 + 5 * 10^0) = -5</code>. 

<p/>
The visual representation above the contents give a simple impression 
of the shape. An <code>X</code> means that the contents are filled up 
to this level, a digit in the topmost row the fraction to which the 
last level is filled. So the 9 of the first column indicates this bin 
is filled 9/10 of the way from <code>3.00*10^2 = 300</code> to 
<code>3.50*10^2 = 350</code>, i.e. somewhere close to 345, 
or more precisely in the range 342.5 to 347.5.

<p/>
The printout also provides some other information, such as the
number of entries, i.e. how many times the histogram has been filled,
the total weight inside the histogram, the total weight in underflow 
and overflow, and the mean value and root-mean-square width (disregarding
underflow and overflow). The mean and width assumes that all the
contents is in the middle of the respective bin. This is especially
relevant when you plot a integer quantity, such as a multiplicity.
Then it makes sense to book with limits that are half-integers, e.g.
<pre>
   Hist multMPI( "number of multiparton interactions", 20, -0.5, 19.5);
</pre>
so that the bins are centered at 0, 1, 2, ..., respectively.  This also 
avoids ambiguities which bin gets to be filled if entries are
exactly at the border between two bins. Also note that the 
<code>fill( xValue)</code> method automatically performs a cast 
to double precision where necessary, i.e. <code>xValue</code> 
can be an integer. 

<h3>The methods</h3>

We here collect a more complete and formal overview of the methods.
   
<a name="method1"></a>
<p/><strong>Hist::Hist() &nbsp;</strong> <br/>
declare a histogram, but does not define it.
  

<a name="method2"></a>
<p/><strong>Hist::Hist(string title, int numberOfBins, double xMin, double xMax) &nbsp;</strong> <br/>
declare and define a histogram, where
<br/><code>argument</code><strong> title </strong>  :  
is a string with the title of the histogram at output,
  
<br/><code>argument</code><strong> numberOfBins </strong>  :  
is the number of bin the <i>x</i> range will be subdivided into, 
limited to be at most 1000,
  
<br/><code>argument</code><strong> xMin </strong>  :  
is the lower edge of the histogram,
  
<br/><code>argument</code><strong> xMax </strong>  :  
is the upper edge of the histogram.
  
  
   
<a name="method3"></a>
<p/><strong>Hist::Hist(const Hist& h) &nbsp;</strong> <br/>
creates an identical copy of the histogram in the argument,
including bin contents.
  
   
<a name="method4"></a>
<p/><strong>Hist::Hist(string title, const Hist& h) &nbsp;</strong> <br/>
creates an identical copy of the histogram in the argument,
including bin contents, except that a new title is provided
as first argument.
  
   
<a name="method5"></a>
<p/><strong>Hist& Hist::operator=(const Hist& h) &nbsp;</strong> <br/>
copies all properties of the histogram in the argument, 
except that the original histogram title is retained. 
  

<a name="method6"></a>
<p/><strong>void Hist::book(string title, int numberOfBins, double xMin, double xMax) &nbsp;</strong> <br/>
define a histogram that previously was only declared; 
see above for the meaning of the arguments.
  

<a name="method7"></a>
<p/><strong>void Hist::name(string title) &nbsp;</strong> <br/>
change the title of a histogram, but keep other properties unchanged.
  

<a name="method8"></a>
<p/><strong>void Hist::null() &nbsp;</strong> <br/>
reset bin contents, but keep other histogram properties unchanged.
  

<a name="method9"></a>
<p/><strong>void Hist::fill(double xValue, double weight) &nbsp;</strong> <br/>
fill the histogram, where 
<br/><code>argument</code><strong> xValue </strong>  : 
is the <i>x</i> position where the filling should occur, and
  
<br/><code>argument</code><strong> weight </strong> (<code>default = <strong>1.</strong></code>) : 
is the amount of weight to be added at this <i>x</i> value.
  
  

<a name="method10"></a>
<p/><strong>friend ostream& operator<<(ostream& os, const Hist& h) &nbsp;</strong> <br/>
appends a simple histogram printout (see above for format) to the 
<code>ostream</code>, while leaving the histogram object itself
unchanged. At most 100 columns are allowed to be displayed. 
If the number of bins is larger than 100 then the contents of 
adjacent bins are added to give the value in each column. (Two by two
up to 200 bins, three by three up to 300, and so on, with the very
last column possibly summing fewer rows than the others.) 
  

<a name="method11"></a>
<p/><strong>void Hist::table(ostream& os = cout) &nbsp;</strong> <br/>
  
<strong>void Hist::table(string fileName) &nbsp;</strong> <br/>
print a two-column table, where the first column gives the center of 
each bin and the second one the corresponding bin contents. The desired
output stream or file name can be provided as argument. The former
is more flexible (e.g., it allows easy append to an existing file), 
whereas the latter is simpler for the case that each histogram should 
be a file of its own. The table may be useful for plotting e.g. with 
Gnuplot.
  

<a name="method12"></a>
<p/><strong>friend void table(const Hist& h1, const Hist& h2, ostream& os = cout) &nbsp;</strong> <br/>
  
<strong>friend void table(const Hist& h1, const Hist& h2, string fileName) &nbsp;</strong> <br/>
print a three-column table, where the first column gives the center of 
each bin and the second and third ones the corresponding bin contents
of the two histograms. Only works if the two histograms have the same
x axis (within a tiny tolerance), else nothing will be done.
  

<a name="method13"></a>
<p/><strong>double Hist::getBinContent(int iBin) &nbsp;</strong> <br/>
return the value in bin <code>iBin</code>, ranging from 1 through 
<code>numberOfBins</code>, with <code>0</code> for underflow and 
<code>numberOfBins + 1</code> for overflow.
  

<a name="method14"></a>
<p/><strong>int Hist::getEntries() &nbsp;</strong> <br/>
return the number of entries, i.e. the number of time that 
<code>fill(...)</code> has been called.
  

<a name="method15"></a>
<p/><strong>bool Hist::sameSize(const Hist& h) &nbsp;</strong> <br/>
checks that the number of bins and upper and lower limits are the 
same as in the histogram in the argument.
  

<a name="method16"></a>
<p/><strong>void Hist::takeLog(bool tenLog = true) &nbsp;</strong> <br/>
by default take 10-logarithm of current contents bin by bin. With 
optional argument <code>false</code> instead take <i>e</i>-logarithm 
of contents bin by bin. If to be used, then right before the
histogram is output. 
  

<a name="method17"></a>
<p/><strong>void Hist::takeSqrt() &nbsp;</strong> <br/>
take square root of current contents bin by bin, with negative contents 
set to zero.
  

<a name="method18"></a>
<p/><strong>Hist& Hist::operator+=(const Hist& h) &nbsp;</strong> <br/>
  
<strong>Hist& Hist::operator-=(const Hist& h) &nbsp;</strong> <br/>
adds or subtracts the current histogram by the contents of the 
histogram in the argument if <code>sameSize(...)</code> is true, 
else does nothing. 
  

<a name="method19"></a>
<p/><strong>Hist& Hist::operator*=(const Hist& h) &nbsp;</strong> <br/>
  
<strong>Hist& Hist::operator/=(const Hist& h) &nbsp;</strong> <br/>
multiplies or divides the current histogram by the contents of the 
histogram in the argument if <code>sameSize(...)</code> is true, 
else does nothing. 
  

<a name="method20"></a>
<p/><strong>Hist& Hist::operator+=(double f) &nbsp;</strong> <br/>
  
<strong>Hist& Hist::operator-=(double f) &nbsp;</strong> <br/>
adds or subtracts each bin content by the common offset <i>f</i>. 
  

<a name="method21"></a>
<p/><strong>Hist& Hist::operator*=(double f) &nbsp;</strong> <br/>
  
<strong>Hist& Hist::operator*=(double f) &nbsp;</strong> <br/>
multiplies or divides each bin content by the common factor <i>f</i>. 
  

<a name="method22"></a>
<p/><strong>friend Hist operator+(double f, const Hist& h1) &nbsp;</strong> <br/>
  
<strong>friend Hist operator+(const Hist& h1, double f) &nbsp;</strong> <br/>
  
<strong>friend Hist operator+(const Hist& h1, const Hist h2) &nbsp;</strong> <br/>
add a constant to a histogram or two histograms to each other, bin by bin.
  

<a name="method23"></a>
<p/><strong>friend Hist operator-(double f, const Hist& h1) &nbsp;</strong> <br/>
  
<strong>friend Hist operator-(const Hist& h1, double f) &nbsp;</strong> <br/>
  
<strong>friend Hist operator-(const Hist& h1, const Hist h2) &nbsp;</strong> <br/>
subtract a histogram from a constant, a constant from a histogram,
or two histograms from each other, bin by bin.
  

<a name="method24"></a>
<p/><strong>friend Hist operator*(double f, const Hist& h1) &nbsp;</strong> <br/>
  
<strong>friend Hist operator*(const Hist& h1, double f) &nbsp;</strong> <br/>
  
<strong>friend Hist operator*(const Hist& h1, const Hist h2) &nbsp;</strong> <br/>
multiply a constant by a histogram or two histograms by each other, 
bin by bin.
  

<a name="method25"></a>
<p/><strong>friend Hist operator/(double f, const Hist& h1) &nbsp;</strong> <br/>
  
<strong>friend Hist operator/(const Hist& h1, double f) &nbsp;</strong> <br/>
  
<strong>friend Hist operator/(const Hist& h1, const Hist h2) &nbsp;</strong> <br/>
divide a constant by a histogram, a histogram by a constant,
or two histograms by each other, bin by bin.
  

</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
