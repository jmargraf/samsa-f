#!/usr/bin/perl
use strict;
use warnings;

our $name = $ARGV[0];
our $basis_set = $ARGV[1];
our $dim = $ARGV[2];

our $natoms;
our @atom;
our @xyz;

our $charge;
our $mult;

  readxyz();
  RunPsi();  

sub readxyz{
  my $scratch;
  my @temp;
  my $m;
  my $ang2bohr=1.889725989;

  open(IN,"<","$::name.inp");

# Skip keyword line
  $scratch=<IN> ;

# Read number of Atoms to $natoms
  $::natoms=<IN> ;
  chomp $::natoms ;

# Read Charge + Mult
#  $scratch=<IN> ;
#  chomp $scratch ;
#  @temp =split / +/, $scratch ;
#  $charge = $temp[0];
#  $mult   = $temp[1];

# Read  $natoms coordinates to xyz
  for($m=1;$m<=$natoms;$m++){
    $scratch=<IN> ;
    chomp $scratch ;
    @temp =split / +/, $scratch ;
    $atom[$m]=  $temp[0] ;
    $xyz[$m][1]=$temp[1] *$ang2bohr ;
    $xyz[$m][2]=$temp[2] *$ang2bohr ;
    $xyz[$m][3]=$temp[3] *$ang2bohr ;
  }

  close IN;                 

}



sub RunPsi{
  my $i;

#  open(LOG,">>","$name.out");
#  print LOG "  Running Psi4 Integral Package \n";
#  close LOG;

  open(PSI,">","$name\_PSI4.inp");
  print PSI "import numpy as np \n";
  print PSI " \n";
#  print PSI "memory 250 mb \n";
  print PSI " \n";
  print PSI "set globals {  \n";
  if($basis_set eq "min"){
    print PSI "  basis sto-3g \n";
  }elsif($basis_set eq "svp"){
    print PSI "  basis def2-SV(P) \n";
  }elsif($basis_set eq "po2"){
    print PSI "  basis 6-311G** \n";
  }elsif($basis_set eq "tzp"){
    print PSI "  basis def2-TZVP \n";
  }elsif($basis_set eq "tzd"){
    print PSI "  basis def2-TZVPD \n";
  }elsif($basis_set eq "qzp"){
    print PSI "  basis def2-QZVP \n";
  }
  print PSI "  units bohr  \n";
  print PSI "}\n";
  print PSI " \n";
  print PSI "molecule $name {\n";
  for($i=1;$i<=$natoms;$i++){
    print PSI "  $atom[$i] $xyz[$i][1] $xyz[$i][2] $xyz[$i][3] \n";
  }
  print PSI "}\n";
  print PSI " \n";
#  print PSI "ref_wfn = psi4.new_wavefunction($name, psi4.get_global_option('BASIS')) \n";
  print PSI "ref_wfn = psi4.core.Wavefunction.build($name, psi4.core.get_global_option('BASIS')) \n";
  print PSI "mints = MintsHelper(ref_wfn.basisset()) \n";
  print PSI " \n";
  print PSI "S = mints.ao_overlap() \n";
  print PSI "np_S = np.array(S) \n";
  print PSI "print(\"OVERLAP\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %.15f' % (i,j,np_S[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "T = mints.ao_kinetic() \n";
  print PSI "np_T = np.array(T) \n";
  print PSI "print(\"KINETIC\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %.15f' % (i,j,np_T[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "V = mints.ao_potential() \n";
  print PSI "np_V = np.array(V) \n";
  print PSI "print(\"POTENTIAL\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %.15f' % (i,j,np_V[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "ERI = mints.ao_eri() \n";
  print PSI "np_ERI = np.array(ERI) \n";
  print PSI "print(\"REPULSION\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "k=0 \n";
  print PSI "count = 0 \n";
  print PSI " \n";
  print PSI "while (i < $dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    k=0 \n";
  print PSI "    while (k < $dim): \n";
  print PSI "      l=0 \n";
  print PSI "      while( l<=k ): \n";
  print PSI "        ij = i*(i+1)/2+j \n";
  print PSI "        kl = k*(k+1)/2+l \n";
  print PSI "        if ij>=kl: \n";
  print PSI "          print '  %d  %d  %d  %d   %.15f' % (i,j,k,l,np_ERI[i][j][k][l]) \n";
  print PSI "          count = count + 1 \n";
  print PSI "        l=l+1 \n";
  print PSI "      k=k+1 \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI " \n";

  close PSI;

  `psi4 $name\_PSI4.inp $name\_PSI4.out > ints`;

#  open(LOG,">>","$name.out");
#  print LOG "    ... done  \n";
#  close LOG;
}

