#!/usr/local/bin/perl
# Modified from http://www.inference.phy.cam.ac.uk/mackay/perl/
# A 2 H . p 
#
# reads in an Afile and writes out a raw ascii 0/1 file 
# see also A2dat.p 
# and raw2A.p which is the inverse

# A2H.p < DSC.273.82.A > DSC.tmp

#  How to override min and max values with command line
eval "\$$1=\$2" while $ARGV[0]=~ /^(\w+)=(.*)/ && shift;

$_ = <>;
($L,$M) = split ; 
print STDERR $L , " " , $M , "\n" ; 
$_ = <>;
$_ = <>;
$_ = <>;
for ( $l = 1 ; $l <= $L ; $l ++ ) {
    $_ = <>;
    @y = split ;
    for ( $u = 0 ; $u <= $#y ; $u ++ ){
	$ty = $y[$u] ;
	$ones{$l,$ty} =  1 ;
    }
}				# 
my $filename = 'H';
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
foreach  $m ( 1..$M ) {
    foreach  $l ( 1..$L ) {
	if ( $ones{$l,$m} ) {
	    print $fh "1" ;
	}else { print $fh "0" ; }
	print $fh " " ;
    }
    print $fh "\n" ;
}
close $fh;