#!/usr/bin/perl

#this reads the 408 file into memory
#

  use Astro::Time;

  $file = $ENV{'FILE408'};

  my @ra_coord;
  my @dec_coord;
  my @pol_angle_eq;
  my @pol_angle_gal;
  my @pol_temp;
  my @l_coord;
  my @b_coord;
    
  my $max_elements = 0;
  
  open (POLN,"<$file") or die "Cannot open Survey File";
  my $line_no = 0;
  my $last_deg = 0;
  while (<POLN>) {

    $line=$_;
    ($page,$deg,$ra,$pol_eq,$temp,$pol_gal,$l,$b) = split(/,/,$line);
    $line_no=$line_no + 1;
    if ($line_no == 1) {
      $last_deg = $deg;
    }
    if ($deg =~ /A/ ) {
      $deg = $last_deg;
    }
    else {
      $last_deg = $deg;
    }
    $line_no=$line_no + 1;

#    print "PAGE=$page DEC=$deg RA=$ra EQ = $pol_eq TEMP=$temp GAL=$pol_gal L=$l B=$b\n"; 
  
# Got all the details in

# Lets parse the coodinates as angles

# RA is 4 character hhhh need to break that up
    $hours = substr($ra,0,2);
    $minutes = substr($ra,2,2);
    $ra_string = "$hours:$minutes:00";
    ($degrees,$dminutes) = split(/ /,$deg);
    
    if ($degrees == 0) {
      $max_elements++;
    }
    if ($dminutes) {
      $dec_string = "$degrees:$dminutes:00";
    }
    else {
      $dec_string = "$degrees:00:00";
    }
    
    push(@ra_coord,$ra_string);
    push(@dec_coord,$dec_string);
    push(@pol_angle_eq,$pol_eq);
    push(@pol_angle_gal,$pol_gal);
    push(@pol_temp,$temp);
    push(@l_coord,$l);
    push(@b_coord,$b);
    
  }
# first we need to know the ranges
  $num_entries = @pol_temp;
  $max_dec = $dec_coord[$num_entries-1];
  $min_dec = $dec_coord[0];

  ($maxdeg,$min) = split (/":"/,$max_dec);
  ($mindeg,$min) = split (/":"/,$min_dec);
  print $maxdeg," -> ",$mindeg, "\n";
  $degrange = $maxdeg-$mindeg;
  
  print"\n---------\n";
  print "There are ", $num_entries, " elements from ", $max_dec," to ", $min_dec, " ($degrange) degrees\n";
  print "The zero declination strip has ", $max_elements, " elements\n";
 
  print "The image plane will be ", $max_elements, " by ", $num_entries, "\n";
  @image_array;
  @angle_array;
  $ra_interval = 2;
  $dec_interval = 1.0;
  $nside = 180;

  for ($ra=0;$ra<360;$ra = $ra+$ra_interval) {
   for($dec = -89.5;$dec<90; $dec=$dec+$dec_interval) {
     push(@image_array,0.0);
     push(@angle_array,0.0);
   }
  }
  $i=0;
  foreach $entry (@pol_temp) {
# what degree element are we in   
    $dec_c = $dec_coord[$i];
    $deg=0;
    $min=0;
    ($deg,$min) = split (/":"/,$dec_c);
    if ($min != 0) {
      $deg = $deg + ($min/60.0);
    }
    $element = ($deg + 89.5) * (1.0/$dec_interval);
    $element = int ($element);
# what ra element are we in
    $turns = str2turn($ra_coord[$i],'H');
    $angle = $turns*360;
    $ra_element = int($angle*(1.0/$ra_interval));

# we are going for column ordering so the element is
# (ra element * 360) + dec element

    $position = ($element * $nside) + $ra_element;

# should change it to flux here but wont
    $image_array[$position] = $entry;
#    $image_array[$position] = $i;
    $angle_array[$position] = $pol_angle_eq[$i];
    print "$i : ra:", $ra_coord[$i], " dec:",$dec_coord[$i], "will place at ", $ra_element,",",$element, " or pos: ",$position,"\n";
    $i = $i + 1; 
  }
  open (OUTP,">polntemp.dat") or die "cannot open output file";
  open (OUTA,">polnangle.dat") or die "cannot open output file";
  
  foreach $entry (@image_array) {
    print OUTP $entry . "\n";
  }

  close (OUTP);
 
  foreach $entry (@angle_array) {
    print OUTA $entry . "\n";
  }

  close (OUTA);
  

    
