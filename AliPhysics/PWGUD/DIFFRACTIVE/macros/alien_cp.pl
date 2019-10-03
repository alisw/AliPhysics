#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use List::Util 'shuffle';
use threads;
use Thread::Queue;

my $SRCDIR = shift or do{  print "USAGE : $0 <alien_folder> [zipfile] [Target] \n";
    print "EXAMPLE : ./$0 /alice/cern.ch/user/b/bschang/PROD23/ root_archive.zip jcorran.root\n";
    exit 1;
};
my $ZIPNAME = shift || "root_archive.zip";
my $TARGET = shift || "jcorran.root";
my $REPEAT = shift || 1;
my $SLEEP = shift || 600;
my $output = "./";

for( 1 .. $REPEAT ){
#### Make source List
my $com = "alien_find -x t $SRCDIR $ZIPNAME 2>/dev/null";
my $xml = `$com`;
my $result = file_from_xml( $xml );
my $list = $result->{list};
my $size_sum  = 0;
for my $item ( @$list ){
	$size_sum += $item->{size};
}
my $nSrc = @$list;
print "## Number of files = ",$nSrc,"\n";
print "## Total Size is ".($size_sum)."\t".($size_sum/1024/1024)."M\n";

### Make Target List

my $com_nLocal = "find $output/$SRCDIR -name $TARGET | wc -l ";
my $nList = `$com_nLocal`; chomp $nList;
my $nTodo = $nSrc-$nList;
my $nDone :shared;
$nDone=0;
printf "## Downloaded/Src = $nList / $nSrc = %.1f%%, TODO : %d\n", $nList/$nSrc*100, $nSrc-$nList;
print "\n";
sleep 2;


my $q = new Thread::Queue;
my $threadcount = 20;
my $targetcount;

my @fail;
my @success;



# Make all Queue
my @list = shuffle( @$list );
#my @list =  @$list;
my $nJob = 0;
for ( 0 .. $#list ){
	my $item = $list[$_];
	$item->{local_dir} = $output.'/'.$item->{dirname};
    next unless ( check_sync($item) ) ;
    $item->{iJob} = $nJob++;
	$q->enqueue( $_+1 );
}

print "$nJob\n";

$threadcount = @$list if @$list <= $threadcount; 


for ( 0..($threadcount - 1)){
	$q->enqueue(undef);
	threads->new(\&do_it, $_);
}

foreach (threads->list){
	$_->join;
}


$nList = `$com_nLocal`; chomp $nList;
my $nTodo2 =  $nSrc-$nList;
printf "## Downloaded/Src = $nList / $nSrc = %.1f%%, TODO : %d\n", $nList/$nSrc*100, $nTodo2;


sub do_it{
	my $index;
	while( $index  = $q->dequeue){
        $index -= 1;
		my $item = $list[$index];
		my $check = 1; #check_sync( $item );

		unless ( $check ) {
			#	print "++++++++++ current file is same one\n";
		} else { 
			system ( 'mkdir -p '.$item->{local_dir} ); # Make DIR

            # COPY
            my $name = $item->{name};
            my $local_dir = $item->{local_dir};
            my $com = 'echo "TGrid::Connect(\"alien://\");TFile::Cp(\"alien://'.$item->{turl}.'\", \"'."$item->{local_dir}/$name".'\");" | root -b -l | grep -v Welcome | grep -v Trying';
            system $com;
            system( "cd $item->{local_dir};yes A | unzip -qq -u $item->{name}; rm $item->{name}" ) if -e $item->{local_dir}.'/'.$item->{name};
            $nDone++;
            print "\n============= $item->{iJob} Done/Todo = $nDone / $nTodo  ===============\n";
            my $check = check_sync( $item );
            if( $check ){
                push @success, $item;
            }else{
                push @fail, $item;
            }
        }

    }
}

sub file_from_xml {
    my $xml = shift;
    my @files = ( $xml =~ /<file(.*?)>/mg );
    my @list;
    for my $file ( @files ){
        push @list, { map { (/(\S*)="(.*?)"/) } split ' ', $file };
        (undef, $list[-1]{dirname},$list[-1]{suffix} ) = fileparse( $list[-1]{lfn} );
        #print $list[-1]{dirname},"\n";
    }
    @list = sort { $a->{'lfn'} cmp $b->{'lfn'} } @list;
    return { list=>\@list };
}

sub file_from_list {
    my $list= shift;
    my @files = split "\n", $list;
    my @list;
    for my $file ( @files ){
        push @list,  { lfn => $file };
        (undef, $list[-1]{dirname},$list[-1]{suffix} ) = fileparse( $list[-1]{lfn} );
        print $list[-1]{dirname},"\n";
    }
    return { list=>\@list };
}
sub check_sync {
    my ($item) = @_;
    print "### No File : $item->{local_dir}\n" and return 1 unless -e $item->{local_dir}."/".$TARGET;
    my $local_file = "$item->{local_dir}/$ZIPNAME";
    if( -e $local_file ){
      my $local_file_size = ( stat $local_file )[7];
      my $grid_file_size = $item->{size};
      print "### Different Size : $local_file = $local_file_size / $grid_file_size \n" and return undef
        unless $local_file_size ==$grid_file_size;
    }
    return undef;
}

print scalar localtime,"\n";

if ( $REPEAT > 1 ){
  for( -($SLEEP/10)..0 ){
    sleep 10;
    print $_,scalar(localtime),"\n";
  }

}


}
