#!/usr/local/bin/perl -w
# downloaded from http://www.inference.phy.cam.ac.uk/mackay/perl/
# G2H.p 
# reads in any G matrix, does
# triangulation then deduces the H matrix.
#
# see also code/cyclic/README
# see also bin/DS.p
# example usage:
# G2H.p N=14 K=10    Gfile=G14.4

# another example:
# G2H.p N=3 K=3 Gfile=H writeH=G
#
# where H is
# 1 1 0
# 0 1 1
# 1 0 1



# special feature added Mon 18/3/02
$writeGminusH = 0 ; # To see how to do this see ~/qecc/README line 153
# or see qecc/DDS 
$Hrows=-1 ;  # specify the number of independent rows expected(?)

$writeH="H" ; # the stem of the pbm file where H is written.
$weight=0; 
$N=7 ; $K=4 ;
$N=0;$K=0; 
$Gfile="G742" ;
$verbose = 1 ; 
$cyclic = 0 ; # if 1, then it rotates the vector to make the rest of the matrix
$stopatM = 0 ;
$printG=0;

eval "\$$1=\$2" while @ARGV && $ARGV[0]=~ /^(\w+)=(.*)/ && shift;
print STDERR "This program uses terms K and M assuming that it is given a K*N G matrix; \n if you give it an H matrix, please interchange M and K.\n" ; 

print STDERR "If you want to convert G<->H for any LARGE matrices then please use\n";
print STDERR "makecode2 instead (mackay/code/MNC/) \n" ;
print STDERR "e.g.  makecode2 -Ain GHC/ram17.5b  -G GHC/Gram17.5b \n" ;

if ( $writeGminusH ) {
    if ($Hrows<0) { print "ERROR, you must specify Hrows=21 or whatever -\n To see how to do this see ~/qecc/README line 153 \n" ;
		    exit(0);
		}
}

$Kactual = 0 ; 
if ( $Gfile ) {
    print "# " if ( $verbose ) ;
    print STDOUT "$Gfile " ;
    print "\n" if ( $verbose ) ;
    open ( H , "< $Gfile" ) ; 
    $i = 0 ; $k = 0 ; 
    if ( $cyclic ) {
	$_ = <H> ; 
	@r = split ;
	$Nr = $#r + 1 ; if ( !$K ) { $K = $Nr ;  }
	for ( $k = 1 ; $k <= $K ; $k ++ ) {
	    print STDOUT "# "  if (( $verbose > 2 )   ) ;
	    $nnn = $Nr - $k  ; 
	    for ( $weight=0, $n = 1 ; $n <= $Nr ; $n ++ ) {
		$nnn ++ ; while ( $nnn >= $Nr ) { $nnn -= $Nr ; }
		$G[$i] = $r[$nnn] ;
		print STDOUT "$G[$i] "  if (( $verbose > 2 ) || $printG ) ;
		$i ++ ;
		$weight +=  $r[$nnn] ;
	    }
	    print STDOUT "\n"  if  (( $verbose > 2 ) || $printG ) ;
	} 
	$k = $K ; 
    } else {
	while ( <H> ) {
	    @r = split ;
	    print STDOUT "# "  if ( $verbose > 2 ) ; 
	    $Nr = $#r + 1 ; 
	    for ( $n = 1 ; $n <= $Nr ; $n ++ ) {
		$G[$i] = $r[$n-1] ;
		print STDOUT "$G[$i] "  if ( $verbose > 2 ) ;
		$i ++ ;
	    }
	    print STDOUT "\n"  if ( $verbose > 2 ) ; $k ++ ;
	}
	print STDOUT "# read in $k rows ($K)\n" ; 
    }
    close ( H ) ;
}
# G[0..N*K-1] vector contains the whole matrix
if ( ($K && !($k == $K)) || ($N && !($Nr == $N) ) ) {
    print "WARNING: K N mismatch $K $N <-> $k $Nr\n" ; 
    exit(0) ; 
}
    $K = $k ; 
    $N = $Nr ; 
$Mactual = $N - $K ; 
if ( $verbose ) {
    print "K = $K, N = $N\n" ; 
}
if($printG>2){exit(0);}
##############################################################################
#              Make trellis-oriented
##############################################################################
for ( $n = 1 ; $n <= $N ; $n ++ ) {
    $taken[$n] = 0 ;
}
# $i $k are reserved symbols, so is $l , $h , $e
$i = 0 ; $iorig = 0 ; 
for ( $k = 1 ; $k <= $K ; $k ++  ) { # run through rows
    if ( $verbose ) {
	if ( !($k%50) ) { print "\n$k" ; }
	print "." ; 
    }
    do { # repeatedly read G[k], modifying it until its first 1 is cool.
	$leftspan = $N+1 ;
	for ( $n = 1 ; $n <= $N ; $n ++ ) { # read in the new codeword k
	    $g[$n] = $G[$iorig] ; $iorig ++ ;
	    if ( $g[$n] ) {
		if ( $n < $leftspan ) { $leftspan = $n ; }
	    }
	}
	$n = $leftspan ; 
	if ( $leftspan > $N ) {
	    print " zero weight row $k ($Kactual)\n" if ( $verbose ) ; $Mactual ++ ;
	    $rowiszero[$k] = 1 ; 
	    $repeat = 0 ;
	}
	elsif ( $taken[$n] ) { # we have a clash and need to do a subtraction
	    $repeat = 1 ; 
	    $kkill = $taken[$n] ; $ikill = ($kkill-1) * $N ; $i = $iorig - $N ; 
# cut out the old guy from this new guy.
	    print " row $k first 1-element $n clashes with $kkill\n"  if ( $verbose > 2 ) ; 
	    for ( $n = 1 ; $n <= $N ; $n ++ ) {
		$G[$i] = ( $G[$ikill] + $G[$i] ) % 2; # modify G
		$ikill ++ ; $i++ ;
	    } 
	    $iorig -= $N ; 
	} else {
	    $rowiszero[$k] = 0 ; 
	    print " row $k first 1-element $n\n"  if ( $verbose > 2 ) ; 
	    $repeat = 0 ; 
	    $taken[$n] = $k ; $k2n[$k] = $n ;
	    $Kactual ++ ;
	}
    } while ( $repeat ) ; 
}
# G is now in trellis-oriented form on its left but not its right.

print "\n"  if ( $verbose ) ; 
&reportM ; 

if ( $stopatM ) { exit(0) ; }
print "half way point\n" ;
$i  = 0 ; 
for ( $k = 1 ; $k <= $K ; $k ++  ) {
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	print $G[$i] , " "  if ( $verbose > 1 ) ; $i ++ ;
    }
    print "\n"  if ( $verbose > 1 ) ;
}
# reverse the process , subtracting low rows from high rows.
for ( $k = $K ; $k >= 1 ; $k -- ) { 
    $i = $N*($k-1) ; 
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	$g[$n] = $G[$i] ; $i ++ ;
    } 
    if ( $k2n[$k] ) {
	$thisn = $k2n[$k] ;
	print " k = $k maps to n = $thisn\n"  if ( $verbose > 2 ) ;
	## BELIEVED that this loop could be replaced by 
##	for ( $kk = $k-1 ; $kk >= 1 ; $kk -- ) {
	    ## was
	for ( $kk = $K ; $kk >= 1 ; $kk -- ) {
	    ## but THIS is correct!!!
	    if ( $kk == $k ) { next ; }
	    $ii = ($kk-1)*$N + $thisn-1 ;
	    if ( $G[ $ii ] ) {
		print " subtracting row $k from row $kk\n" if ( $verbose > 2 ) ; 
		$ii = ($kk-1)*$N ;
		for ( $n = 1 ; $n <= $N ; $n ++ ) {
		    $G[$ii] = ( $g[$n] + $G[$ii] ) % 2 ; $ii++ ;
		} 
	    } else {
		if ( $verbose > 4 ) {
		    print "not subtracting row $k from row $kk cos G[$ii] = $G[$ii]\n" ; 
		}
	    }
	}
    }
}

# Mon 18/3/02
# for my qecc project on self-orthogonal codes,
# this is an interesting point at which to extract the matrix that generates
# the code apart from the dual.
# Assuming the top Mdual lines provided were actually the dual,
# we can print out the remaining non-zero lines

print "done to systematic form: N=$N, M=$Mactual\n" ;
if ( $verbose > 1 ) {
    $i  = 0 ; 
    for ( $k = 1 ; $k <= $K ; $k ++  ) {
	for ( $n = 1 ; $n <= $N ; $n ++ ) {
	    print $G[$i] , " " ; $i ++ ;
	}
	print "\n" ;
    }
}


if ( $writeGminusH ) {
    if ($Hrows<0) { print "ERROR, you must specify Hrows=21 or whatever\n" ;}
    else {
	$MM = 0 ; 
	print "Writing Code Minus Dual\n"  ;
	$hout = "" ;
	$minweight = $N ;  # prepare to find smallest weight in this matrix
	for ( $k = $Hrows+1 ; $k <= $K ; $k ++  ) {
	    if ( !$rowiszero[$k] ) {
		print " $k"; $MM ++ ;
		$i = ($k-1) * $N ;  

		$tweight = 0 ;
		for ( $n = 1 ; $n <= $N ; $n ++ ) {
		    $hout .=  $G[$i]." " ; $i ++ ;
		    $tweight += $G[$i] ;
		}
		$hout .= "\n" ;
		if ($tweight < $minweight ) { $minweight = $tweight ; }
	    }
	}
	print "\n $MM rows, minimum weight $minweight\n" ; 
	if ( !($writeGminusH =~ /\.pbm$/) ) {
	    open ( H , "> $writeGminusH" ) ;
	    print H $hout ;
	    close(H) ; 
	    print "written GminusH matrix to $writeGminusH\n" ;
	    $writeGminusH .= ".pbm" ; 
	}
	
	# add the pbm header
	$hout = "P1\n$N $MM\n" . $hout ;
	open ( H , "> $writeGminusH" ) ;
	print H $hout ;
	close(H) ; 
	print "xv $writeGminusH\n" ;
    }
}



$M = $N - $K ; 
$i  = 0 ; 
for ( $k = 1 ; $k <= $Mactual ; $k ++  ) {
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	$H[$i] = 0 ;  $i ++ ;
    }
}


# now allocate the untaken columns to rows of H
$m = 0 ; 
for ( $n = 1 ; $n <= $N ; $n ++ ) {
    if ( !$taken[$n] ) {
	$m ++ ;
	$n2m[$n] = $m ;
	# set H( m , n ) = 1 ; 
#	printf "m=%d n=%d %d %d\n" , $m,$n,  ($m-1)*$N + ($n-1) , $H[ ($m-1)*$N + ($n-1) ] ; 
	$H[ ($m-1)*$N + ($n-1) ] = 1 ; 
    }
    else { $n2m[$n] = 0 ; 
#	printf "m=%d n=%d\n" , $m,$n ; 
       }
}
$M = $m ;

# run through all elements of G. 
$i = 0 ; 
for ( $k = 1 ; $k <= $K ; $k ++  ) {
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	if ( $G[$i] ) {
	    if ( $taken[$n] ) { # do nothing, it's systematic bit 
	    } else {
                 # set H(m(n), k2n(k) ) = 1 
		$m = $n2m[$n] ; $nn = $k2n[$k] ; 
		$H[ ($m-1)*$N + ($nn-1) ] = 1 ; 
	    }
	}
	$i ++ ;
    }
}
for ( $n = 1 ; $n <= $N ; $n ++ ) {
    $pint[$n] = 0 ; 
}
print "H matrix in systematic form\n" ;
$hout = "" ; 
$i  = 0 ; $two = 1 ; 
for ( $k = 1 ; $k <= $M ; $k ++  ) {
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	if ( $H[$i] ) { $pint[$n] += $two ; }
	$hout .= $H[$i]." " ; $i ++ ;
    }
    $hout .= "\n" ;
    $two *= 2 ; 
}
print $hout if($verbose>1) ; 
if ( $writeH ) {
    if ( !($writeH =~ /\.pbm$/) ) {
	open ( H , "> $writeH" ) ;
	print H $hout ;
	close(H) ; 
	print "written $N $M  matrix to $writeH\n" ;
	$writeH .= ".pbm" ; 
    }
    # add the pbm header
    $hout = "P1\n$N $M\n" . $hout ;
    open ( H , "> $writeH" ) ;
    print H $hout ;
    close(H) ; 
    print "xv $writeH\n" ;
}
if ( $verbose > 1 ) {
    print "H matrix in decimal\n" ;
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	print  $pint[$n] , " " ; 
    }
}
print "\n" ;
&reportM ; 
sub reportM { 
    $Kactual = $N - $Mactual ; 
    print "N: $N\tw: $weight\tM: $Mactual\tK: $Kactual\n" ;
}


