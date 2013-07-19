#!/usr/bin/perl

use strict;

my $FREQ = shift;
my $WMAP_PATH = shift;
my $SIM = $ENV{SIM};

my $WMAP = "$WMAP_PATH/wmap_band_smth_.fits";

print "Processing sky for " . $FREQ . " MHz WMAP dir is " . $WMAP . "\n";

unlink("wmap_I.fits");
unlink("wmap_Q.fits");
unlink("wmap_U.fits");
unlink("wmap_NOBS.fits");

system("$SIM/bin/HPXcvt -b1 -cg $WMAP wmap_I.fits");
system("$SIM/bin/HPXcvt -b2 -cg $WMAP wmap_Q.fits");
system("$SIM/bin/HPXcvt -b3 -cg $WMAP wmap_U.fits");
system("$SIM/bin/HPXcvt -b4 -cg $WMAP wmap_NOBS.fits");

unlink("sky_I.fits");
unlink("sky_Q.fits");
unlink("sky_U.fits");
unlink("sky_V.fits");
unlink("sky_PA.fits");
unlink("sky_DEG.fits");
unlink("angelica.fits");
unlink("allsky.dat");

my $pid = open (GSM,"| diffuse_sky") or die "Cannot fork diffuse_sky. Can you run it from the command line?";

print GSM $FREQ . "\n";

my $gsm_filename = "allsky.dat";

print GSM $gsm_filename . "\n";
close(GSM);


my  $command_line = "convert_hpx -F -f $FREQ -i $gsm_filename -o angelica.fits";
system($command_line);

unlink("old_sky_I.fits");

system("$SIM/bin/HPXcvt -cg angelica.fits old_sky_I.fits");
system("wmap_pol -I old_sky_I.fits -i $WMAP_PATH/wmap -o sky");
#### Comment out these lines if you donot have trace_pol
system("trace_pol -q sky_Q.fits -u sky_U.fits -r $WMAP_PATH/faraday.fits -f $FREQ -c 4 -s 4");

system("cp test_u.fits sky_U.fits");
system("cp test_q.fits sky_Q.fits");
######################################################
unlink("old_sky_I.fits");

unlink("sky_I_$FREQ.fits");
unlink("sky_Q_$FREQ.fits");
unlink("sky_U_$FREQ.fits");
unlink("sky_V_$FREQ.fits");

my $pid1 = fork();
if ($pid1 == 0) {
  system("convert_hpx -f $FREQ -H -i sky_I.fits -o sky_I_$FREQ.fits");
  return 0;
}
else {
  sleep 1;
  my $pid2 = fork();
  if ($pid2 == 0) {
    system("convert_hpx -f $FREQ -H -i sky_Q.fits -o sky_Q_$FREQ.fits");
    return 0;
  }
  else {
    sleep 1;
    my $pid3 = fork();
    if ($pid3 == 0) {
      system("convert_hpx -f $FREQ -H -i sky_U.fits -o sky_U_$FREQ.fits");
      return 0;
    }
    else {
      sleep 1;
      my $pid4 = fork();

      if ($pid4 == 0) {
	system("convert_hpx -f $FREQ -H -i sky_V.fits -o sky_V_$FREQ.fits") ;
	return 0;
      }
      else {
    

	waitpid($pid1,0);
	waitpid($pid2,0);
	waitpid($pid3,0);
	waitpid($pid4,0);
      }
    }
  }
}





