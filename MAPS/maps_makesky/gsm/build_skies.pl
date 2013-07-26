#!/usr/bin/perl -w 

my $num=192;
my $freq_inc = 0.04;
my $freq = 320.0;
my $i=0;
system("date > readme.sky");
open (LOG,">readme.sky");
print LOG "\n$num channels starting at $freq attempted\n";
close LOG;

system("tar -cvf justins_sky.tar readme.sky");

for ($i=0; $i< $num ;$i = $i+1) {
  
  my $input_I = "input_I.dat";
 
  open (GSM,"| diffuse_sky") or die "Cannot fork diffuse_sky. Can you run it from the command line?";

  print GSM $freq . "\n";
  print GSM $input_I . "\n";
  close(GSM);
  
  my $freqlab = sprintf("%7.3f",$freq);
  my $outname = "sky_" . "$freqlab";

  system("convert_hpx -f $freqlab -F -i $input_I -o healpix.fits");
  system("rm hpx.fits");
  system("./HPXcvt healpix.fits hpx.fits");
  my $cleanup = $outname . "_?.fits " . $outname . "_PA.fits " . $outname . "_DEG.fits";
 
  system("rm $cleanup");

  system("wmap_pol -I hpx.fits -i ./wmap -o $outname");
  
   
  die;

  $freq = $freq + $freq_inc;

}

