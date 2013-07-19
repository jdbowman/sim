sub build_obs_spec {
use Astro::Time;
use Astro::Misc;
# The input hash should have all the info we need
# it has hopefully been passed by reference

  my (%in_spec) = %{$_[0]};

# The output hash
  %obs_spec = (	"FOV_center_RA","0:0:0",
      		"FOV_center_Dec","0:0:0",
		"FOV_size_RA", 0.0,
		"FOV_size_Dec", 0.0,
		"Corr_int_time", 2.0,
		"Corr_chan_bw", 0.0,
		"Scan_start","yyyy:ddd:hh:mm:ss",
		"Scan_duration", 0.0,
		"Channel", "fff:bw",
		"PNT_center_RA","0:0:0",
		"PNT_center_Dec","0:0:0");

 # what is the RA/Dec of the zenith
 # I know the LST - hence that is RA at zenith
 # the latitude of the observatory is the Dec of that point

  $turn = $time/24.0;
  
  $obs_spec{"FOV_center_RA"} = turn2str($turn,'H',0);
  
  $fov = $in_spec{"FOV.asec"};

  $deg = $in_spec{"Lat.deg"} + $in_spec{"field.off"};
  $obs_spec{"FOV_center_Dec"} = deg2str($deg,'D',0);
  $obs_spec{"FOV_size_RA"} = $in_spec{"FOV.asec"};
  $obs_spec{"FOV_size_Dec"} = $in_spec{"FOV.asec"};
  $obs_spec{"PNT_center_RA"} = $in_spec{"PNT_center_RA"};
  $obs_spec{"PNT_center_Dec"} = $in_spec{"PNT_center_Dec"};
#  $obs_spec{"Corr_int_time"} = $in_spec{"int.time.secs"};
  $obs_spec{"Corr_chan_bw"} = $in_spec{"ch.bw.MHz"};

# Scan start is interesting as it has to be in GMT and probably
# should match LST.
# there is a GHA mode which I will try first
# recall LHA = LST - RA(star)
# GHA = GST - RA(star)

# LHA = GHA - LONG(west)
#
# LST = GST - LONG (WEST)
#
#
# Lets go straight to UT

 ($sec,$min,$hour,$day,$month,$year,$weekday,$yearday,$daylight) = localtime(time);
# my $long = $in_spec{'Lon.deg'};
# $mjd = lst2mjd($turn,$yearday,$year+1900,$long);
# ($day,$year,$ut) =  mjd2dayno($mjd);
# $obs_spec{"Scan_start"} = sprintf("$year:%03d:%s",$day,turn2str($ut,'H',2));

# Using GHA mode 
 my $default = -7.8237889;
 my $gha_offset = $default;
 
 $obs_spec{"Scan_start"} = "GHA ". $gha_offset . " ";

 $obs_spec{"Scan_duration"} = $in_spec{"int.time.secs"};
# Channel is easy too
 $obs_spec{"Channel"} = sprintf ("%7.3f:%f",$in_spec{"freq.MHz"},$in_spec{"ch.bw.MHz"});


# Now we need to write it out
 my $specname = $in_spec{"NAME"};
 open (SPEC,">$specname") or die "Cannot open obs. spec. file.: $!";
 print SPEC "//Auto generated Obs. Spec. File \n";
 print SPEC "FOV_center_RA = $obs_spec{'FOV_center_RA'} \n";
 print SPEC "FOV_center_Dec = $obs_spec{'FOV_center_Dec'} \n";

 if ($obs_spec{"PNT_center_RA"} < 24) {
   
   print SPEC "PNT_center_RA = $obs_spec{'PNT_center_RA'} \n";
   print SPEC "PNT_center_Dec = $obs_spec{'PNT_center_Dec'} \n";
 
 }
 
 print SPEC "FOV_size_RA = $obs_spec{'FOV_size_RA'} \n";
 print SPEC "FOV_size_Dec = $obs_spec{'FOV_size_Dec'} \n";
 print SPEC "Corr_int_time = $obs_spec{'Corr_int_time'} \n"; 
 print SPEC "Corr_chan_bw = $obs_spec{'Corr_chan_bw'} \n";
 print SPEC "\n";
 print SPEC "Scan_start = $obs_spec{'Scan_start'} \n";
 print SPEC "Scan_duration = $obs_spec{'Scan_duration'} \n";
 print SPEC "Channel = $obs_spec{'Channel'} \n";
 print SPEC "Endscan"; 

 close (SPEC);
 return ;
}
1;


 
