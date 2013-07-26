#!/usr/bin/perl -w 
my $ra_to_track = 8.0; 
my $dec_to_track = -26.0;
my $freq = 150;

my $time_start = 6; 
my $time_stop = 10;
my $step = 0.5;
my $T = 512;

my $time = $time_start;
open (LOG,">run.log");
while ($time < $time_stop) {
  print LOG "construct_sky.pl -O oob.txt -T $T -p 1 -N 1 -n 1 -b 60 -f $freq -R $ra_to_track -D $dec_to_track -t $time \n";
  
  my $offset = ($time - $ra_to_track) * 15;
  my $dec = $dec_to_track + 26.4331 ;
  print LOG "test_rts -s -o \"$offset:$dec\" -r 2 -F 30 -S -t 1 -b sky_V_$freq -q $freq\n";

  print LOG "mv regrid_I.fits regrid_I_$time.fits\n";
  print LOG "mv regrid_Q.fits regrid_Q_$time.fits\n";
  print LOG "mv regrid_U.fits regrid_U_$time.fits\n";
  print LOG "mv regrid_V.fits regrid_V_$time.fits\n";
  print LOG "mv sky_V_$freq*.uvfits sky_$time" . "_$freq.000_00.uvfits\n";
  
  $time = $time+$step;
}

close LOG;


