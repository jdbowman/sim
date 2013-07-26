#!/usr/bin/perl

# This is a perl script to generate test skys for visgen
# It uses the GSM diffuse sky model ans a list of point sources 
# to produce skys with parameters selected on the command line
# 
# this is marshalling a number of excutables:
#
# build_pol  ---   this produces the all-sky 
#	 	   sky at any frequency
#
# maps_makesky --- this produces actual input sky from the output
#		   of build_pol
#
# test_points  --- If you want to add a source list (from the pipeline package)
#
# Steve Ord sord@cfa.harvard.edu October 2007


use Getopt::Std;
$sim_bin = $ENV{'SIM_BIN'};
$sim_dir = $ENV{'SIM'};

require $sim_bin . "/build_obs_spec.pl";

if (@ARGV == 0) {
print "-a = 400 unless defined : altitude of observatory (meters)\n";
print "-b = 60 unless defined  : Resolution (arcsec) \n";
print "-c = 0.032 unless defined : channel width (MHz) \n";
print "-d = 0 unless defined : offset of field centre from zenith (deg) \n";
print "-D = 99 unless defined : declination of pointing centre (-90 - 90 ) \n";
print "-f = 200 unless defined :  base frequency channel 0 (MHz) \n";

print "-l = -26.4331 unless defined : latitude of observatory (degrees) \n";

print "-L = 242.643 unless defined  : longitude(west) of the observatory \n";

print "-M = 0 unless defined : shall we make a movie \n";

print "-N = 300 unless defined : Number of time steps \n";

print "-n = 16 unless defined  : number of channels \n";

print "-o = channel unless defined : output root name \n";

print "-O = add OOB sources via OOB file \n";
  
print "-p = 0 unless defined : polarised sky \n";

print "-Q = 0 unless defined : DO NOT use Stokes Q in sky \n";

print "-R = 25 unless defined : RA of pointing centre (0 -> 24) \n";

print "-s = 2.0 unless defined  : time cadence (seconds) \n";

print "-S = source_list.txt unless defined  : adding sources \n";

print "-t = 17.0 unless defined : start LST (hours) \n";

print "-T = 512 unless defined : number of telescopes \n";

print "-u = 1 unless defined : run maps2uvfits? \n";

print "-U = 0 unless defined: do not use Stokes U in sky \n";
 
print "-v = 1 unless defined : run visgen? \n";

print "-w = 180 unless defined : FOV (degrees) \n";

print "-W = PATH to WMAP data \n";
die "no command line arguments \n";
}

getopts("a:b:c:d:D:f:l:N:n:M:o:O:p:Q:R:U:s:S:t:T:u:v:w:W:",\%args);
$spec{"NAME"} = "auto_obs_spec";

$args{a} = 400 unless defined $args{a}; # altitude of observatory

$args{b} = 60 unless defined $args{b}; #Resolution (arcsec)
$spec{"res.asec"}=$args{b};
$args{c} = 0.032 unless defined $args{c}; # channel width (MHz)
$spec{"ch.bw.MHz"}=$args{c};

$args{d} = 0.0 unless defined $args{d}; #field off zenith by ... (degrees)
$spec{"field.off"}=$args{d};

$args{D} = 100.0 unless defined $args{D};  #beam declination (-90 90)
$spec{"PNT_center_Dec"}=$args{D};

$args{f} = 200 unless defined $args{f}; # base frequency channel 0 (MHz)
$spec{"freq.MHz"}=$args{f};

$args{l} = -26.4331 unless defined $args{l}; # latitude of observatory (degrees)
$spec{"Lat.deg"} = $args{l};

$args{L} = 242.643 unless defined $args{L}; # longitude(west) of the observatory
$spec{"Lon.deg"} = $args{L};

$args{M} = 0 unless defined $args{M}; #shall we make a movie
$make_movie = $args{M};
$args{N} = 300 unless defined $args{N}; # Number of time steps
$spec{"N.time.steps"} = $args{N};

$args{n} = 16 unless defined $args{n}; # number of channels
$spec{"N.Ch."} = $args{n};

$args{o} = "channel" unless defined $args{o}; # output root name
$args{O} = "oob.txt" unless defined $args{O};

$args{p} = -1 unless defined $args{p}; # polarisation approx percentage

$args{R} = 25 unless defined $args{R}; # RA of Beam
$spec{"PNT_center_RA"} = $args{R};

$args{s} = 2.0 unless defined $args{s}; # time cadence (seconds)
$spec{"int.time.secs"} = $args{s};

$args{S} = "source_list.txt" unless defined $args{S}; # adding sources

$args{t} = 17.0 unless defined $args{t}; # start LST (hours)
$spec{'LST.hours'} = $args{t};
$args{T} = 512 unless defined $args{T}; # number of telescopes
$args{u} = 1 unless defined $args{u}; # run maps2uvfits?
$args{v} = 1 unless defined $args{v}; # run visgen?
$args{w} = 180 unless defined $args{w}; # FOV (degrees)
$args{W} = "$sim_dir/maps_makesky/gsm/" unless defined $args{W}; 
#lets build a hash

my $diffuse = 1;
my $poln = 1;
my @skies;
my @polnum;

print "Channel width " . $args{c} . "\n";
print "Start frequency " . $args{f} . "\n";
print "Telescope latitude " . $args{l} . "\n";
print "Number of channels  " . $args{n} . "\n";
print "LST at start  " . $args{t} . "\n";
print "Length of integration " . $args{N}*$args{s} . "\n";
print "This will result in " . $args{N}*$args{n} . " files " . "\n";

# run for each time step
# run gsm as many times as you can cover the frequency range required

$chanN = 0;
open(LOGIT,">construct.log") or die "Cannot open construct.log: $!";

while ($chanN < $args{n}) {

  $freq = $args{f} + ($chanN * $args{c});
  if ($diffuse==1 && $poln==0) {

    $pid = open (GSM,"| diffuse_sky") or die "Cannot fork diffuse_sky. Can you run it from the command line?";

    print GSM $freq . "\n";

    $gsm_filename = "allsky.dat";

    print GSM $gsm_filename . "\n";
    close(GSM);


    stat($gsm_filename) or die "Could not stat the GSM file, does it exist: $!";	
# Now we can convert this to a fits file

#$convert_filename = sprintf "%s_%07.3f_%02.0f.fits",$args{o},($args{f}+$chanN*$args{c}),($timeN*$args{s});


    $convert_filename = "!allsky.fits";    

# dont need to write to it so I can just use system to launch it

    $command_line = "convert_hpx -f $args{f} -i $gsm_filename -o $convert_filename";
  
    print LOGIT $command_line . "\n";

    $status = system("$command_line");
    die "convert_hpx failed" unless $status == 0;

    push(@skies,"allsky.fits");

  }
  elsif ($diffuse == 1 && $poln == 1) {
    if (-e "sky_I_$freq.fits") {
      break;
    }
    else {
      $command_line = "build_pol.pl $freq $args{W}";
      system($command_line);
    } 
    push(@skies,"sky_I_$freq.fits");
    push(@skies,"sky_Q_$freq.fits");
    push(@skies,"sky_U_$freq.fits");
    push(@skies,"sky_V_$freq.fits");
    
    push(@polnum,"I");
    push(@polnum,"Q");
    push(@polnum,"U");
    push(@polnum,"V");

  }
  
  $timeN = 0;

  while ($timeN < $args{N}) {

# now we have to run maps_makesky on the output
# it takes a number of arguments that may have to be calculated

# -i input_sky_fits_file 
# -o output_image_name 
# -w width(and height)_output_image
# -t LST (hours)
# -l latitude (degrees)
# -s text_source_file_name
# -a angluar_size_of_pixels (arcsec)
# -x (flip X/RA axis) -d (turns debugging on) 
    print @skies;
    unlink("metafile");
    my $skynum = 0;
      $fov = $args{w}*60.0*60.0; 
# Requested FOV in arc seconds
      $dim = $fov/$args{b}; 
# FOV / Resolution (fixed)

# it is probably best if this is a power of two

#      $powerof2 = log($dim)/log(2);

#     $power = int ($powerof2);

#      $power = $power+1;

#      $newdim = 2**$power;

# old school
       $newdim = $dim;

# change resolution for fixed FOV

      my $newb = $fov/$newdim;

#   so 

      $arg_dim = " -w " . $newdim; 
# new width in pixels

#    $fov = $newdim * $args{b}; 
# new FOV in arcsecs

      $spec{"FOV.asec"} = $fov;

  $timeinhours = $timeN*$args{s}/(60.0*60.0);

      $time = $args{t} + $timeinhours;

      if ($time >= 24.0) {
	$time = $time - 24.0;
      }


 
    foreach my $sky (@skies) {
    
      $cmd_makesky = "maps_makesky ";

# we know the input filename because we built it
      my $thesky = $sky;
      print "Processing $thesky";
      
      $arg_input = " -i $thesky"; 

      $cmd_makesky = $cmd_makesky . $arg_input;

# we get to choose the output name
      (my $skyname, my $skyextn) = split(/\./,$thesky);
      $out_filename = sprintf "%s_%07.3f_%02.0f",$skyname,($args{f}+$chanN*$args{c}),($timeN*$args{s});

      $arg_output = " -o " . $out_filename . ".fits";

      $cmd_makesky = $cmd_makesky . $arg_output;

# width and height of output image. Well this wants to cover at 
# at least the whole primary beam of the antenna. 
# probably more.

# The next stage is im2uv which maks a UVFITs file this just needs to 
# know the area of the FOV.

# finally the obs_spec file also needs the FOV and the RA and dec
# of the field centre
# this will need calcualting based upon the information supplied here
# so it may be worthwhile constructing the obs_spec files here too
# FOV
     $cmd_makesky = $cmd_makesky . $arg_dim;

# lets not forget the angular size of a pixel
#changing back to old school
      $arg_ang = " -a " . $args{b};

      $cmd_makesky = $cmd_makesky . $arg_ang;

# we know time from start LST + $timeN*cadence
# it wants the time in hours so we have to check we do not run over 24

          $arg_time = sprintf (" -t %8.5f" , $time);;

      $cmd_makesky = $cmd_makesky . $arg_time;

# we know the latitude from our own command line

      $field_cntr_kludge = $args{l} + $args{d};
      $arg_latitude = " -l " . $field_cntr_kludge;

      $cmd_makesky = $cmd_makesky . $arg_latitude . " -x ";

# Is there a source file
      $cmd_list = "test_points -t ";
      $do_list = 0;
      if (defined $args{S}) {
	stat($args{S}) or warn "Could not stat the source list file, does it exist: $!";	
	if (stat($args{S})) {

	  $arg_source = " -p " . $args{S};
	  $cmd_list = $cmd_list . $arg_source . " " ;
	  $do_list = 1;
	}
	else {
	  $do_list = 0;
	}
      }


# therefore are maps_makesky line is
      if ($diffuse == 1 ) {
	print LOGIT $cmd_makesky . "\n";

	$status = system("$cmd_makesky");
      }
    
    
# Now to add the sources as maps_makesky doesnt really know how
# LOSim does know which is how I should really do it
# but at the moment I'm just going to run the test_points executable 

# Now we need the im2uv with appropriate flux scaling
      $im2uv = "MAPS_im2uv -t 4 ";

# Need the area of the FOV in radians
# assume 180.0 degree FOV
    
      my $area = 2*3.141;
      my $flux_norm = 2.0*1.38065e-23*1e26/(299.792/$freq)**2;
      $scaling = $area*$flux_norm;

# need the input file


      if ($diffuse == 1) {

	if ($do_list == 1) {

	  $cmd_list = $cmd_list . $out_filename . ".fits";

	  print LOGIT $cmd_list . "\n";
	  system("$cmd_list");
	  $tempname = "new." . $out_filename;
	  $out_filename = $tempname;

	}

	$cmdline = $im2uv . "-i " . $out_filename . ".fits " . " -o " . $out_filename . ".dat " . " -n " . $scaling;
      
	print LOGIT $cmdline . "\n";
	system("$cmdline");

      }
      if ($make_movie == 1) {
	$jpeg_cmd = "fits2jpeg -fits " . $out_filename . ".fits -jpeg " . $timeN . ".jpg ";
	system("$jpeg_cmd");
	$convert = "convert " . $timeN . ".jpg " . $timeN . ".gif";
	system("$convert");  
      }
      open(META,">>metafile") or die "cannot open metafile";
      my $metaline = $polnum[$skynum] . " " . $out_filename . ".dat";
      print META "$metaline\n";
      close(META);
      $skynum=$skynum+1;
    }
# this should probably not be run by a human - check the args
    $spec{"NAME"} = $out_filename . ".spec";
    &build_obs_spec(\%spec); 
 
    # Now I'm going to run visgen here it may be worthwhile forking to make use of multiple CPU
    if (diffuse == 0) {
	$out_filename = "no_diffuse";
    }
    $VISGEN{'NAME'} = $out_filename;
    $VISGEN{'SITE'} = "MWA_ED";
    
    if ($args{T} == 512) {
      $array_name = $sim_dir . "/array/mwa_random_500_crossdipole_gp.txt";
    }
    else {
      $array_name = $sim_dir. "/array/mwa_32_crossdipole_gp.txt";
    }

    $VISGEN{'ARRAY'} = $array_name;

    if ($args{p} >= 0) {
      $VISGEN{'GRIDFILE'} = "none";
      $VISGEN{'METAFILE'} = "metafile";
    }
    else {
      $VISGEN{'GRIDFILE'} = $out_filename . ".dat";
      $VISGEN{'METAFILE'} = "none";
    }

    $VISGEN{'OBS_SPEC'} = $spec{'NAME'};

    $VISGEN_CMD = "visgen -N -n $VISGEN{'NAME'} -s $VISGEN{'SITE'} -A $VISGEN{'ARRAY'} -V $VISGEN{'OBS_SPEC'}";
    if ($diffuse == 1 ) { 
      if ($VISGEN{'METAFILE'} =~ /none/) {
	$VISGEN_CMD = $VISGEN_CMD . " -G $VISGEN{'GRIDFILE'}";
      }
      else {
	$VISGEN_CMD = $VISGEN_CMD . " -M $VISGEN{'METAFILE'}";
      }
    }
    else {
     $VISGEN_CMD = $VISGEN_CMD . " -Z ";
    }

    if (defined $args{O}) {
      if (stat($args{O})) {
	$VISGEN_CMD = $VISGEN_CMD . " -O " . $args{O};
      }
      else {
	warn "Cannot stat $args{O}:";
      }
    }
    
    $VISGEN_CMD = $VISGEN_CMD . "   -m 0 > /dev/null ";

    print LOGIT " $VISGEN_CMD \n";
    
    if ($args{v} == 1) {
      system("$VISGEN_CMD");
    }
    $longitude_east = 360 - $args{L};
    $MAPS2UVFITS = "maps2uvfits " . $out_filename . ".vis " . $out_filename .".uvfits" . " $args{l} $longitude_east $args{a} " . $VISGEN{'ARRAY'} ;
    
    print LOGIT $MAPS2UVFITS . "\n";;
    
    if ($args{u} == 1) {
      system("$MAPS2UVFITS");
    }
    (my $rr,my $pp,my $ff,my $fff,my $tt) = split ("_",$out_filename);
    
    my $cleanupname = $rr . "_[IQUV]_" . $ff . "_" . $fff . "_" . $tt;
    chomp($cleanupname);
    my $cleanupline = "rm " . $cleanupname . ".dat " . $cleanupname . ".fits";
   print ($cleanupline);   
   system($cleanupline);

    $timeN = $timeN + 1;
  
  } # this has produced all the time steps



  $chanN = $chanN + 1;

} # this has produced a whole bunch of channels
close(LOGIT);
# convert

# run maps_makesky - but add some point sources at specific places

# calculate the time range required

# loop over the time range to build the skys





