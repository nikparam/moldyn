#!/usr/bin/perl -w
#
#              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                           MEWTON-X controler
#               Mario Barbatti, M. Ruckenbauer, 2005-2013
#                     G09 added by R. Crespo-Otero
#              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#====================================================================================
#                              MAIN PROGRAM
#====================================================================================
#
# pragmas *********************************************************************
#
use lib join('/',$ENV{"NX"},"lib") ;
use colib_perl;
use hybrid_module;
use strict;
use warnings;
use Cwd;

#
# variable declarations *******************************************************
my ($mld,$BASEDIR,$TP,$RS,$RT,$DEBG,$OD,$ctd,$ctdyn,$dynot,$dynmd,$v0,$dynvlt,$anmod);
my ($JAD,$JNAD);
my ($i,$mdle,$retval);
my (@prmt);
my ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tini,$tmax,$nintc,$found_cigrad);
my ($mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,$etot_drift);
my ($getphase_def,$ms_def,$ms,$Ms,$kross,$cprog,$coptda,$current,$cascade,$never_state,$include_pair,$e_ci,$wpout,$integ);
my ($n_lines_geom,$n_lines_veloc,$nndamp,$nnat);
my ($tprevious,$nstatdynprevious,$tnp,$tn,$tlasthop,$timeafterhop);
my ($getphase,$type,$found,$found1,$found2,$koj,$adjmom);
my ($cpar,$memcol,$cirestart,$reduce_tol,$mc_conv,$ci_conv,$quad_conv,$ivmode,$mocoef,$niter,$td_st,$ld_thr,$kind_g09);
my ($n_lines_wfr,@re_coef,@im_coef,$n_st,$pop,$d_pop,$grb);
my ($progname,$methodname,$prt_mo,$nad_exec,$jobtype);
my ($ktherm,$thermostat,$thermkind,$temp,$lvp,$gamma,$tout);
my ($radius,$thinp,$kts,$lts,$nstherm,$colprob);
my ($vdoth,$ci_cons,$cio_options,$cisc_options,$idalton);
my ($hybrid_columbus, $hybrid_tinker, $hybrid_turbomole, $hybrid_analytical, $hybrid_dftbp, @hybcoljobs,$jobdir,$jobdir1,$jobdir2);
my ($jobs_r, $natoms_r, $propjob_r, $nadregions_r, $jobn, $method, $mcontr_file, $submethod);
my ($gauchk);
my ($g09chk);
my ($g09rwf);
my ($g09com); 
my ($boundary_conditions,$boundary_type,$sphere_radius,@sphere_center);
my ($nfrozen,$frozen_cart,$thisstep);
my ($GAMESS,$gamesspar,$credits,$verno,$ncpus,$run_gamess,$scr,$scrdir);
my ($kross_d,$cprog_d,$coptda_d,$cascade_d,$current_d,$never_state_d,$include_pair_d,$e_ci_d,$ci_cons_d,$cio_options_d,$cisc_options_d,$idalton_d);
my $epsilon = 1E-9;
sub interpolation($$$$);

# setting some variables ******************************************************
set_variables();

# preparing directories *******************************************************
prepare_dirs();

# copying files ***************************************************************
copy_files();

chdir("${BASEDIR}/${TP}/");

# initialize outputs **********************************************************
initialize_outputs();

# check input files ***********************************************************
check_files();
  # In this subroutine:
  #   check input files
  #   check initial parameters
  #   check programs for non-adiabatic dynamics
  #   check initial and final times
  #   check inconsistent inputs
  #   check phase calculation
  #   check number of states
  #   check initial wave function
  #   Create file for counting number of hoppings
  #   Damped dynamics
  #   Move redundant files
#   Some of these have been moved to the run-<method>.pl modules
#   #   Copy ab initio program input
#   #   check if there is either a mocoef file or scf keyword present
  #   Check thermostat
#   #   check third-part program

# %%%%%%%%%%%%%%%%%%%%%%%%%%  Step 0     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Execute third-party program in step 0
exec_prog("initial run");

# freeze gradient components if requested
freeze_gradient();

# freeze NAD-vector components if requested
freeze_nadv();

# Compute internal coordinates and perform an optional analysis in step 0
   analysis();
#fp internal_coordinate();

# Execute surface hoping routine in step 0
sh_step_0();

# Initialize molecular orbitals for Columbus
initial_mocoef();

# %%%%%%%%%%%%%%%%%%%%%%%%%%  Main loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_main_loop();

# %%%%%%%%%%%%%%%%%%%%%%%%%%  End loop   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ($t != $tini){
  type_of_dyn();
}
titlestep();

print_STDOUT("$mdle End of dynamics \n");

end_nx();

#
#====================================================================================
#                              END OF MAIN PROGRAM
#====================================================================================
#
#====================================================================================
#                              SUB ROUTINES
#====================================================================================
#
sub set_variables{
$mdle   = "moldyn.pl:";

$mld    = $ENV{"NX"};       # NX environment
$BASEDIR    = &getcwd();    # Basedir
$TP     = "TEMP";           # TEMP
$RS     = "RESULTS";        # RESULTS
$DEBG   = "DEBUG";          # DEBUG
$RT     = "INFO_RESTART";   # RESTART
$JAD    = "JOB_AD";
$JNAD   = "JOB_NAD";

$ctd    = "control.d";
$ctdyn  = "control.dyn";
$dynot  = ">>$BASEDIR/$RS/dyn.out";
$dynmd  = ">>$BASEDIR/$RS/dyn.mld";
$v0     = ">veloc0";
$dynvlt = "dynvlt";
$anmod  = "analytical.model";
$thinp  = "therm.inp";

$getphase_def = 1;   # Default of getphase. Also in sh.f90 and escalar.f90.
$ms_def       = 20;  # Default of ms. Also in sh.f90.

undef($frozen_cart);
$frozen_cart='';
$nfrozen=0;

$lvprt=getkeyword("${BASEDIR}/control.dyn","lvprt",1);

}
#
#
#====================================================================================
sub prepare_dirs{
#
# if RESULTS, DEBUG, and RESTART are soft linked they should be used
# TEMP is alway erased and created from scratch
$nxrestart = getkeyword($ctdyn,"nxrestart","0");

unless (-e "$BASEDIR/$RS"){
  system ("mkdir $BASEDIR/$RS");
}
if ($nxrestart == 0){
  system("rm -rf $BASEDIR/$RS/*");
}
#
unless (-e "$BASEDIR/$DEBG"){
  system("mkdir $BASEDIR/$DEBG");
}
if ($nxrestart == 0){
  system("rm -rf $BASEDIR/$DEBG/*");
}
#
unless (-e "$BASEDIR/$RT"){
  system ("mkdir $BASEDIR/$RT");
}
if ($nxrestart == 0){
  system("rm -rf $BASEDIR/$RT/*");
}
#
system("rm -rf $BASEDIR/$TP");
system("mkdir $BASEDIR/$TP");
}
#
#====================================================================================
#
sub copy_files{
 my ($jobadexst, $jobnadexst, $jobrtdir, $jobexdir);
# Copy input files
 system("cp -f $BASEDIR/control.dyn $BASEDIR/$TP/.");
 if ($nxrestart == 0){
   $OD=$BASEDIR;
 }else{
   $OD=$RT;
   if (!-s $RT){exit_error("$ctdyn says nxrestart = 1, but $RT is empty.",$mdle);}
   $prog=getkeyword("${BASEDIR}/control.dyn","prog",-1);
   if($prog < 0){die "no program prog = ?? in control.dyn\n";}
   $progname=prog_name($prog);
   $methodname=method_name($prog);
 }
 system("cp -rf $OD/JOB_* $OD/geom $OD/veloc $BASEDIR/$TP/.");
 if ($OD ne $RT){
 system("cp -rf $OD/JOB_* $OD/geom $OD/veloc $OD/control.dyn $BASEDIR/$RT/.");
 }
# copy the nad_vectors for hybrid jobs to TEMP/JOBEX_<jobn>.<method>/
 if(($nxrestart != 0) and ($progname eq "hybrid"))
 {
   $jobadexst = -1;
   $jobnadexst = -1;
   if( -e "$JAD/hybrid.control") {
     $jobadexst=1;
   }
   if (-e "$JNAD/hybrid.control") {
     $jobnadexst=1;
   }
   if($jobadexst > 0)
   {
     chdir($JAD);
   }
   elsif($jobnadexst > 0)
   {
     chdir($JNAD);
   }
   else {die "file hybrid.control seems not to exist!\n";}
  # read in hybrid.control -> get definitions for jobs
   ($jobs_r,$natoms_r,$propjob_r,$nadregions_r)=&read_control_file();
   for($i=1;$i<@$jobs_r;$i++)
   {
     if (defined($$jobs_r[$i]))
     {
       $jobn=$i; # this is only for comfort, I could use $i directly
       $method=lc($$jobs_r[$jobn][1]);
       if($method eq 'columbus')
       {
         if ($jobnadexst > 0)
         {
           mkdir("$BASEDIR/$TP/JOBEX_$jobn.$method");
           $jobrtdir="$BASEDIR/$RT/hybrid.$jobn.$method";
           $jobexdir="$BASEDIR/$TP/JOBEX_$jobn.$method";
           if(-e "$jobrtdir/nad_vectors")
           {
             system("cp -f $jobrtdir/nad_vectors $jobexdir/");
           }
           if(-e "$jobrtdir/pcnad_vectors")
           {
             system("cp -f $jobrtdir/pcnad_vectors $jobexdir/")
           }
         }
       }
       elsif($method eq 'tinker')
       {
         # nothing else to be done
       }
       elsif($method eq 'turbomole')
       {
         if($jobnadexst > 0)
         {
           mkdir("$BASEDIR/$TP/JOBEX_$jobn.$method");
           $jobrtdir="$BASEDIR/$RT/hybrid.$jobn.$method";
           $jobexdir="$BASEDIR/$TP/JOBEX_$jobn.$method";
           if(-e "$jobrtdir/nad_vectors")
           {
             system("cp -f $jobrtdir/nad_vectors $jobexdir/");
           }
         }
       }
       elsif($method eq 'analytical')
       {
         if($jobnadexst > 0)
         {
           mkdir("$BASEDIR/$TP/JOBEX_$jobn.$method");
           $jobrtdir="$BASEDIR/$RT/hybrid.$jobn.$method";
           $jobexdir="$BASEDIR/$TP/JOBEX_$jobn.$method";
           if(-e "$jobrtdir/nad_vectors")
           {
             system("cp -f $jobrtdir/nad_vectors $jobexdir/");
           }
           if(-e "$jobrtdir/pcnad_vectors")
           {
             system("cp -f $jobrtdir/pcnad_vectors $jobexdir/");
           }
         }
       }
       elsif($method eq 'dftb-plus')
       {
         if($jobnadexst > 0)
         {
           mkdir("$BASEDIR/$TP/JOBEX_$jobn.$method");
           $jobrtdir="$BASEDIR/$RT/hybrid.$jobn.$method";
           $jobexdir="$BASEDIR/$TP/JOBEX_$jobn.$method";
           if(-e "$jobrtdir/nad_vectors")
           {
             system("cp -f $jobrtdir/nad_vectors $jobexdir/");
           }
           if(-e "$jobrtdir/pcnad_vectors")
           {
             system("cp -f $jobrtdir/pcnad_vectors $jobexdir/");
           }
         }
       }
       else
       {
         warn "method $method unknown\n";
         die "currently implemented: columbus, tinker, turbomole, analytical, dftb-plus\n";
       }
     }
   }
   chdir(${BASEDIR});
 }

 if (-e "$OD/sh.inp")      {system("cp -f $OD/sh.inp       $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/sh.inp       $BASEDIR/$RT/.");}}
 if (-e "$OD/jiri.inp")    {system("cp -f $OD/jiri.inp     $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/jiri.inp     $BASEDIR/$RT/.");}}
 if (-e "$OD/wf.inp")      {system("cp -f $OD/wf.inp       $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/wf.inp       $BASEDIR/$RT/.");}}
 if (-e "$OD/turbomole.par")  {system("cp -f $OD/turbomole.par   $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/turbomole.par   $BASEDIR/$RT/.");}}
 if (-e "$OD/analyt.par")  {system("cp -f $OD/analyt.par   $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/analyt.par   $BASEDIR/$RT/.");}}
 if (-e "$OD/columbus.par"){system("cp -f $OD/columbus.par $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/columbus.par $BASEDIR/$RT/.");}}
 if (-e "$OD/dftb.par")    {system("cp -f $OD/dftb.par     $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/dftb.par     $BASEDIR/$RT/.");}}
 if (-e "$OD/gau.par")     {system("cp -f $OD/gau.par      $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/gau.par      $BASEDIR/$RT/.");}}
 if (-e "$OD/g09.par")     {system("cp -f $OD/g09.par      $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/g09.par      $BASEDIR/$RT/.");}}
 if (-e "$OD/+dftb+.par")   {system("cp -f $OD/dftb+.par    $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/dftb+.par     $BASEDIR/$RT/.");}}
 if (-e "$OD/gamess.par")     {system("cp -f $OD/gamess.par      $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/gamess.par      $BASEDIR/$RT/.");}}
 if (-e "$OD/rndsrst")     {system("cp -f $OD/rndsrst      $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/rndsrst      $BASEDIR/$RT/.");}}
 if (-e "$OD/$thinp")      {system("cp -f $OD/$thinp       $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/$thinp       $BASEDIR/$RT/.");}}
 if (-e "$OD/rndtm")       {system("cp -f $OD/rndtm        $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/rndtm        $BASEDIR/$RT/.");}}
 if (-e "$OD/nad_vectors")       {system("cp -f $OD/nad_vectors        $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/nad_vectors  $BASEDIR/$RT/.");}}
 if (-e "$OD/intcfl")      {system("cp -f $OD/intcfl       $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/intcfl       $BASEDIR/$RT/.");}}
 if (-e "$OD/boundaries.inp") {system("cp -f $OD/boundaries.inp  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/boundaries.inp     $BASEDIR/$RT/.");}}
 if (-e "$OD/freeze.inp") {system("cp -f $OD/freeze.inp  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/freeze.inp     $BASEDIR/$RT/.");}}
 if (-e "$OD/therm.freeze") {system("cp -f $OD/therm.freeze  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/therm.freeze     $BASEDIR/$RT/.");}}
 if (-e "$OD/print_options.inp") {system("cp -f $OD/print_options.inp  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/print_options.inp     $BASEDIR/$RT/.");}}
 if (-e "$OD/freezing_options.inp") {system("cp -f $OD/freezing_options.inp  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/freezing_options.inp     $BASEDIR/$RT/.");}}
 if (-e "$OD/stopsign.inp") {system("cp -f $OD/stopsign.inp  $BASEDIR/$TP/.");
   if ($OD ne $RT){system("cp -f $OD/stopsign.inp     $BASEDIR/$RT/.");}}
# rndsrst is the random number state. It allows to reinitialize the previous
# pseudo-random sequence from some specific point. Useful for restartings.
# rndtm is the equivalent for the thermostat.
#
# Move to TEMP directory
 chdir("$BASEDIR/$TP");
#
# to deal with old NX inputs: replace &end by / in namelists
 namelist_leg("control.dyn",$lvprt);
 if (-s "sh.inp"){namelist_leg("sh.inp",$lvprt);}
 if (-s "jiri.inp"){namelist_leg("jiri.inp",$lvprt);}
 if (-s "therm.inp"){namelist_leg("therm.inp",$lvprt);}
}
#
#====================================================================================
#
sub initialize_outputs{
  if ($nxrestart != 0){
    print_STDOUT("                           ---- Restarting NX job -----\n");
  }
  title_nx($mld);
  print_STDOUT("\nSTARTING MOLECULAR DYNAMICS   \n\n");
  # Print paths
  print_STDOUT("$mdle paths:\nBase:    $BASEDIR\nTemp:    $BASEDIR/$TP\n");
  print_STDOUT("Results: $BASEDIR/$RS\nDebug:   $BASEDIR/$DEBG\nRestart: $BASEDIR/$RT\n");
  print_STDOUT("NEWTON-X: $mld \n");
}
#
#====================================================================================
#
sub check_files{
my (@JOBS,$flag);
my $epsilon = 1E-9;

# check input files
open(F1,"geom") or die "$mdle Cannot open geom!\nsystems answer was: $!";
close(F1);

open(F2,"veloc") or die "$mdle Cannot open veloc!\nsystems answer was: $!";
close(F2);

open(F3,"control.dyn") or die"$mdle Can't open control.dyn!\nsys answer was:$!";
close(F3);

print_STDOUT("$mdle Input files: Checked (1) \n \n");

print_STDOUT("$mdle Reading intial parameters \n \n");

# read control.d (this block is also read in run-col.pl)
# $retval = callprogsp($mld,"inp",$mdle);
# This is now a routine in lib/colib_perl.pm instead of a fortran program.
# Now also comments are enabled and the namelist format is not necessary any more.
&rewrite_ctd();

&read_ctd();

# Check programs for non-adiabatic dynamics
&nad_progs();

print_STDOUT("\tInitial parameters:\n");
print_STDOUT("\tNat       = $nat\n");
print_STDOUT("\tistep     = $istep\n");
print_STDOUT("\tnstat     = $nstat\n");
print_STDOUT("\tnstatdyn  = $nstatdyn\n");
print_STDOUT("\tndamp     = $ndamp\n");
print_STDOUT("\tkt        = $kt\n");
print_STDOUT("\tdt        = $dt\n");
print_STDOUT("\tt         = $t\n");
print_STDOUT("\ttmax      = $tmax\n");
print_STDOUT("\tnint      = $nintc\n");
print_STDOUT("\tmem       = $mem\n");
print_STDOUT("\tnxrestart = $nxrestart\n");
print_STDOUT("\tthres     = $thres\n");
print_STDOUT("\tkillstat  = $killstat\n");
print_STDOUT("\ttimekill  = $timekill\n");
print_STDOUT("\tprog      = $prog\n");
print_STDOUT("\tlvprt     = $lvprt\n");
print_STDOUT("\tetot_jump = $etot_jump\n");
print_STDOUT("\tetot_drift= $etot_drift\n\n");
                                            #

if (($nad_exec eq "yes") and ($thres > 0)){
   print_STDOUT("Non-adiabatic dynamics was requested.\n\n");
   $jobtype="nonadiabatic";
}elsif(($nad_exec eq "no") or ($thres == 0)){
   print_STDOUT("Adiabatic dynamics was requested.\n\n");
   $jobtype="adiabatic";
}

# write file with initial time
open(TINI,">tini") or die "$mdle Cannot open tini\nsystems answer was: $!";
$tini = $t;
print TINI "$tini \n";
close(TINI);

# Check dt
if ($dt == 0.0){
  exit_error("Timestep dt=$dt. Fix it and run again.",$mdle);
}

# Check tmax
if ($tmax < $t){
  exit_error("tmax ($tmax) cannot be smaller than t ($t).
     Check $ctd and run again.",$mdle);
}

# write the information whether to print or not to auxp
open (AP,">auxp") or die "$mdle: cannot write auxp\n$!\n";
$thisstep=sprintf("%9.2f",${t}/${dt});
# print STDOUT "DEBUGMR: writing auxp fmodulo($thisstep,$kt) = ",&fmodulo($thisstep,$kt),"\n";
$flag=&fmodulo($thisstep,$kt);
if($flag < $epsilon) {$flag = 0;}
printf AP "%9.2f\n",$flag;
close (AP);

# Check if JOB_AD and JOB_NAD exist. Kill when necessary. If JOB_AD does not
# exist, but JOB_NAD does, copy NAD as AD and put this
# information in kind_of_JOB to be read by run-col.pl
# koj = 0 - Only JOB_AD or JOB_AD+JOB_NAD
# koj = 1 - Only JOB_NAD
if (!-e "$JAD"){
   if ($nad_exec eq "yes"){
      if (!-e "$JNAD"){
         exit_error("$JAD and $JNAD directories\n     were not found!
     Prepare the input files for $progname calculation",$mdle);
      }
      else {
        system("cp -rf $JNAD $JAD");
        print_STDOUT(
"\n     ****************************************************
     $mdle $JAD was not found.
     $JNAD will be used.
     If you want adiabatic dynamics, you must have $JAD.
     NX will continue with nonadiabatic dynamics.
     **************************************************** \n");
        if ($thres==0){
          exit_error("ERROR: When JOB_AD does not exist, thres
     must be set to some large value in control dyn
     (e.g., thres = 100 eV).",$mdle);
        }
        open(KJ,">kind_of_JOB") or die "$mdle Cannot open kind_of_JOB to write!";
        $koj = 1;
        print KJ " $koj";
        close(KJ);
     }
  }
  if (($progname eq "turbomole")
   or ($progname eq "dftb")){
     exit_error("$JAD directory was not found!
    Prepare the input files for $progname calculation",$mdle);
  }
}
else {
   open(KJ,">kind_of_JOB") or die "$mdle Cannot open kind_of_JOB to write!";
   $koj = 0;
   print KJ " $koj";
   close(KJ);
}

# check number of lines in geom and veloc
$n_lines_geom=0;
$n_lines_veloc=0;

open(F1,"geom" ) or die "$mdle Cannot open geom\nsystems answer was: $!";
while (<F1>) {$n_lines_geom++; }
close(F1) or die "$mdle Could not close geom\nsystems answer was: $!";
open(F2,"veloc") or die "$mdle Cannot open veloc\nsystems answer was: $!";
while (<F2>) {$n_lines_veloc++;}
close(F2) or die "$mdle Could not close veloc\nsystems answer was: $!";

if ($nat != $n_lines_geom ){
  exit_error("Number of atoms in geom and control.d are different!",$mdle);
}
if ($nat != $n_lines_veloc){
  exit_error("Number of atoms in veloc and control.dyn are different!",$mdle);
}

# Change defaults for TDDFT surface hopping
if (($jobtype eq "nonadiabatic") and ($progname eq "turbomole")){
  if (!-s "sh.inp"){
    open(SHI,">sh.inp") or die "Cannot open sh.inp to write\nsystems answer was: $!";
    print SHI "&shinp \n";
    print SHI " vdoth = 1\n";
    print SHI " adjmom=-1\n";
    print SHI "/ \n";
    close(SHI);
  }
  if (!-s "jiri.inp"){
    open(JINP,">jiri.inp" ) or die "$mdle Cannot open jiri.inp to write\nsystems answer was: $!";
    print JINP "&jirinp \n never_state = 1\n/ \n";
    close(JINP);
  }
}

if (($jobtype eq "nonadiabatic") and ($progname eq "g09")){
  if (!-s "sh.inp"){
    open(SHI,">sh.inp") or die "Cannot open sh.inp to write\nsystems answer was: $!";
    print SHI "&shinp \n";
    print SHI " vdoth = 1\n";
    print SHI " adjmom=-1\n";
    print SHI "/ \n";
    close(SHI);
  }
  if (!-s "jiri.inp"){
    open(JINP,">jiri.inp" ) or die "$mdle Cannot open jiri.inp to write\nsystems answer was: $!";
    print JINP "&jirinp \n never_state = 1\n/ \n";
    close(JINP);
  }
}

# jiri.inp file
if (!-e "jiri.inp"){
  open(JINP,">jiri.inp" ) or die "$mdle Cannot open jiri.inp to write\nsystems answer was: $!";
  print JINP " &jirinp \n/ \n";
  close(JINP);
}

# These defaults (below) are also defined in typeofdyn.f90 and write_pair.f90
($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d,$e_ci_d,$ci_cons_d,
$cio_options_d,$cisc_options_d,$idalton_d,$cprog_d,$coptda_d)=load_defaults("jiri.inp",$prog);
$kross        = getkeyword("jiri.inp","kross",       $kross_d);
$cascade      = getkeyword("jiri.inp","cascade",     $cascade_d);
$current      = getkeyword("jiri.inp","current",     $current_d);
$never_state  = getkeyword("jiri.inp","never_state", $never_state_d);
$include_pair = getkeyword("jiri.inp","include_pair",$include_pair_d);
$e_ci         = getkeyword("jiri.inp","e_ci",        $e_ci_d);
$ci_cons      = getkeyword("jiri.inp","ci_cons",     $ci_cons_d);
$cio_options  = getkeyword("jiri.inp","cio_options", $cio_options_d);
$cisc_options = getkeyword("jiri.inp","cisc_options",$cisc_options_d);
$idalton      = getkeyword("jiri.inp","idalton",     $idalton_d);
$cprog        = getkeyword("jiri.inp","cprog",       $cprog_d);  
$coptda       = getkeyword("jiri.inp","coptda",      $coptda_d); 

print_STDOUT("\n$mdle Parameters for nonadiabatic dynamics (jiri.inp):\n");
print_STDOUT("\tkross        = $kross\n");
print_STDOUT("\tcascade      = $cascade\n");
print_STDOUT("\tcurrent      = $current\n");
print_STDOUT("\tnever_state  = $never_state\n");
print_STDOUT("\tinclude_pair = $include_pair\n");
print_STDOUT("\te_ci         = $e_ci\n");
print_STDOUT("\tci_cons      = $ci_cons\n");
print_STDOUT("\tcio_options  = $cio_options\n");
print_STDOUT("\tcisc_options = $cisc_options\n");
print_STDOUT("\tidalton      = $idalton\n");
print_STDOUT("\tcprog        = $cprog\n");
print_STDOUT("\tcoptda       = $coptda\n");

# Write first transmomin file
if ($jobtype eq "nonadiabatic"){
  write_transmomin();
}

# Phase information
if (!-e "sh.inp"){
  open(SHI,">sh.inp") or
    die "Cannot open sh.inp to write\nsystems answer was: $!";
  print SHI "&shinp \n";
  print SHI " getphase = $getphase_def \n";
  print SHI "/ \n";
  close(SHI);
}
$found=searchkeyword("sh.inp","getphase");
if ($found == 0){
  system("cp -f sh.inp sh.inp.bk");
  open(SBK,"sh.inp.bk") or
    die"$mdle Cannot open sh.inp.bk\nsystems answer was: $!";
  open(SHI,">sh.inp") or
    die "$mdle Cannot open sh.inp to write\nsystems answer was: $!";
  while(<SBK>)  {
    chomp $_;
    if (/shinp/i)    {
        print SHI "$_ \n";
      # write default value of getphase
        print SHI " getphase = $getphase_def\n";
    }
    else    {
      print SHI "$_ \n";
    }
  }
  close(SHI);
  close(SBK);
}

# check phase calculation
$getphase = &keyinfile("sh.inp","getphase");
if ($progname eq "columbus" or $hybrid_columbus > 0){
  if ($getphase == 0)  {
    print_STDOUT("\n$mdle getphase = 0:\n");
    print_STDOUT("$mdle Phase will be computed from CI-vectors overlap\n\n");
  }
  if ($getphase == 1)  {
    print_STDOUT("\n$mdle getphase = 1:\n");
    print_STDOUT("$mdle Phase will be computed from h-vectors scalar prod.\n\n");
  }
}

# check adjustment of momentum after hopping
$found=searchkeyword("sh.inp","adjmom");
if ($found != 0)  {
  $adjmom = &keyinfile("sh.inp","adjmom");
  if ($adjmom == 1){
   print_STDOUT("$mdle *** WARNING ***: The meaning of adjmom keyword in sh.inp\n");
   print_STDOUT("         changed since NX version 0.11b. Check the new documentation.\n");
   print_STDOUT("         If you want momentum adjustment along h vector, set adjmom=0.\n\n");
  }
}

# Check adjmom / vdoth consistency
$vdoth  = getkeyword("sh.inp","vdoth","0");
$adjmom = getkeyword("sh.inp","adjmom","0");
if (($adjmom >= 0) and ($vdoth != 0)){
  print_STDOUT("$mdle *** WARNING ***: With vdoth = $vdoth adjmom cannot be >= 0. Assuming adjmom = -1.\n");
}

# Check ms / vdoth consistency
$vdoth  = getkeyword("sh.inp","vdoth","0");
$ms     = getkeyword("sh.inp","ms",$ms_def);
if (($vdoth == -1) and ($ms != 0)){
  print_STDOUT("$mdle *** WARNING ***: With vdoth = $vdoth ms must be 0. Assuming ms = 0.\n\n");
}

# Check integrator / vdoth consistency
$vdoth  = getkeyword("sh.inp","vdoth","0");
$integ  = getkeyword("sh.inp","integrator","5");
if (($vdoth == -1) and ($integ != 6)){
  print_STDOUT("$mdle *** WARNING ***: With vdoth = $vdoth integrator must be 6. Assuming integrator = 6.\n\n");
}
if (($vdoth != -1) and ($integ == 6)){
  print_STDOUT("$mdle *** WARNING ***: With integrator = $integ vdoth must be -1. Assuming vdoth = -1.\n\n");
}

# check number of states
if ($nstat < $nstatdyn ){
  exit_error("nstat must be larger or equal nstatdyn.",$mdle);
}

# check initial wave function

if (-e "wf.inp"){
  system("cp -f wf.inp wfrun");
  # Is the number of states right?
  $n_lines_wfr=0;
  open(WFR,"wfrun") or die "Cannot open wfrun to read.";
  while(<WFR>){$n_lines_wfr++;}
  close(WFR);
  if ($nstat != $n_lines_wfr){
    exit_error("Number of states in wf.inp and control.d are different!",$mdle);
  }
  # Put one blank before each line in wfrun. Just to avoid reading erros.
  open(WFR,"wfrun") or die "Cannot open wfrun to read.";
  open(WFRB,">wfrun.b") or die "Cannot open wfrun.b to write.";
  while(<WFR>){
    print WFRB " $_";
  }
  close(WFR);
  close(WFRB);
  system("mv -f wfrun.b wfrun");
  # Is the normalization right?
  open(WFR,"wfrun") or die "Cannot open wfrun to read.";
  $n_st=1;
  $pop=0.0;
  while(<WFR>){
    chomp;
    ($grb,$re_coef[$n_st],$im_coef[$n_st])=split(/\s+/,$_);
    $pop=$pop+$re_coef[$n_st]**2+$im_coef[$n_st]**2;
    $n_st++;
  }
  close(WFR);
  $d_pop=abs(1-$pop);
  if ($d_pop >= 0.05){ # This is the tolerance for the population normalization.
    exit_error("Coefficients in wf.inp are not normalized
                  Total population = $d_pop",$mdle);
  }
  print_STDOUT("$mdle Using TD-coeffcients read from wf.inp. \n\n");
}
else{
  print_STDOUT("$mdle Using default TD-coefficients. \n\n");
  open(WF,">wfrun");
  $i=1;
  while ($i <= $nstat)  {
    if ($i == $nstatdyn)    {
      print WF "1.0 0.0 \n";
    }
    else    {
      print WF "0.0 0.0 \n";
    }
    $i++;
  }
}
#
print_STDOUT("$mdle Input files: checked (2) \n \n");

# Create file for counting number of hoppings
if (!-e "nhopsold"){
  open(NHO,">nhopsold") or die "$mdle Cannot open nhopsold.\nsystems answer was: $!";
  print NHO "0 0 \n";
  close(NHO)
}

# Damped dynamics
$nndamp=$ndamp;
if ($nndamp == 1){
  print_STDOUT("$mdle Damped dynamics \n \n");
  $nnat=$nat;
  open(VELNULL,$v0) or
    die "$mdle Cannot open $v0!\nsystems answer was: $!";
  for ($i=1; $i <= $nnat; $i++)  {
    print VELNULL "0.00 0.00 0.00 \n";
  }
    system("cp -f veloc0 veloc");
}

# Move redundant files
if ((-e $JAD)  and (-e $JNAD)){
  @JOBS=($JAD,$JNAD);
}elsif((-e $JAD)  and (!-e $JNAD)){
  @JOBS=($JAD);
}elsif((!-e $JAD) and (-e $JNAD)){
  @JOBS=($JNAD);
}
foreach(@JOBS){
  if (-e "$BASEDIR/$TP/$_/geom" ){
    system("mv -f $BASEDIR/$TP/$_/geom  $BASEDIR/$TP/$_/geom.org" );
  }
  if (-e "$BASEDIR/$TP/$_/coord"){
    system("mv -f $BASEDIR/$TP/$_/coord $BASEDIR/$TP/$_/coord.org");
  }
  if ($progname eq "columbus" or $hybrid_columbus > 0) {
    if ($progname eq "columbus")
    {
      $jobdir="$BASEDIR/$TP/$_";
      if(-e "$jobdir/control.run")
      {
        system("cp -f $jobdir/control.run $jobdir/control.run.org");
      }
      if ($koj == 0)
      {
        if(-e "$jobdir/ciudgin")
        {
          system("cp -f $jobdir/ciudgin $jobdir/ciudgin.org");
        }
      }
    }
    else
    {
      for ($i=0;$i<@hybcoljobs;$i++)
      {
        $jobdir="$BASEDIR/$TP/$_/JOB_$hybcoljobs[$i].columbus";
        if (-e "$jobdir/control.run")
        {
          system("cp -f $jobdir/control.run $jobdir/control.run.org");
        }
        if ($koj == 0)
        {
          if(-e "$jobdir/ciudgin")
          {
            system("cp -f $jobdir/ciudgin $jobdir/ciudgin.org");
          }
        }
      }
    }# if($progname ... else ...
  }
  if (-e "$BASEDIR/$TP/$_/transmomin"){
    system("mv -f $BASEDIR/$TP/$_/transmomin $BASEDIR/$TP/$_/transmomin.org");
  }
  if (-e "$BASEDIR/$TP/$_/transmomin"){
    system("mv -f $BASEDIR/$TP/$_/transmomin $BASEDIR/$TP/$_/transmomin.org");
  }
}
if ($lvprt>=2) {
print_STDOUT("$mdle Redundant files in TEMP/JOB_AD & TEMP/JOB_NAD saved as *.org.\n\n");
}

# Set initial type for dynamics at timestep 0
$type = 1;

# This section is transferred to run-col.pl and executed there if $t=$tini
# -> see run-col.pl:&check_columbus()
# matruc 24apr08
# if ($progname eq "columbus") {
#   # check if there is a mocoef file present
#   unless (-e "$BASEDIR/$TP/$JAD/mocoef"){
#       exit_error("no mocoef file present",$mdle);
#   }
#   # check number of iterations
#   $niter    = getkeyword("$BASEDIR/$TP/$JAD/control.run","niter","1");
#   if ($niter != 1){
#     exit_error("NITER must be equal 1 in control.run (COLUMBUS input)",$mdle);
#   }
#   # Type of gradient
#   $found_cigrad = 0;
#   $found_cigrad=searchkeyword("$BASEDIR/$TP/$JAD/control.run","cigrad");
#   if (($found_cigrad == 0) and (-e $JNAD)){
#     $found_cigrad=searchkeyword("$BASEDIR/$TP/$JNAD/control.run","nadcoupl");
#   }
#   if (($found_cigrad == 0) and ($nstat != 1)){
#     if ($tini != $tmax){
#        exit_error("NSTAT > 1 is allowed only for CI gradient or nonadiabatic coupling \n(cigrad or nadcoupl in control.run)",$mdle);
#     }
#   }
# }

# check boundaries
($boundary_conditions,$boundary_type,$sphere_radius,@sphere_center)=check_boundary();

# check freezing of coordinates
($nfrozen,$frozen_cart)=check_freezing();

# check restrictions on writing
&check_cutback();

# check thermostat
check_thermostat();

# Check third-party programs
check_thirdparty();

}
#
#====================================================================================
#
sub read_ctd{
# read control.d
  my ($curdir,@line);
  $hybrid_columbus = 0;
  $hybrid_tinker   = 0;
  $hybrid_turbomole = 0;
  $hybrid_analytical = 0;
  @hybcoljobs = ();
  $curdir = &getcwd();
  open(CT, $ctd) or die "$mdle Cannot open $ctd!\nsystems answer was: $!";
  $_=<CT>;
  chomp;s/\s+//g;
  @prmt=split /,/;
  close(CT);
  ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc,
  $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,$etot_drift)=@prmt;

  $progname=prog_name($prog);
  $methodname=method_name($prog);

  if($progname eq "hybrid")
  {
    if( -e "$JAD/hybrid.control") {chdir("$JAD");}
    elsif (-e "$JNAD/hybrid.control") {chdir("$JNAD");}
    else {die "file hybrid.control seems not to exist!\n";}
    # read in hybrid.control -> get definitions for jobs and which job to
    # take the nad-vectors, and properties from
    ($jobs_r,$natoms_r,$propjob_r,$nadregions_r)=&read_control_file();
    for($i=1;$i<@$jobs_r;$i++)
    {
      if (defined($$jobs_r[$i]))
      {
        $jobn=$i; # this is only for comfort, i could use $i directly
        $method=lc($$jobs_r[$jobn][1]);
        print_STDOUT("job $jobn, region $$jobs_r[$jobn][0], method $$jobs_r[$jobn][1]\n");
        if($method eq 'columbus')
        {
          $hybrid_columbus = 1;
          push(@hybcoljobs,$jobn);
        }
        elsif($method eq 'tinker')
        {
          $hybrid_tinker = 1;
        }
        elsif($method eq 'turbomole')
        {
          $hybrid_turbomole = 1;
          if($jobn=$$propjob_r)
          {
            $mcontr_file="job_$jobn.$method.control";
            open( MCONTR, $mcontr_file ) or die "cannot open $mcontr_file\n$!\n";
            while (<MCONTR>)
            {
              chomp;
              s/^\s//;
              if (/^method/i)
              {
                s/\s+//g;
                @line=split/=/;
                if ( $line[1] eq "ricc-2" )
                {
                  $submethod = 2.0;
                  $hybrid_turbomole = 2;
                }
                elsif ( $line[1] eq "td-dft" )
                {
                  $submethod = 2.1;
                  $hybrid_turbomole = 2;
                }
                else { die " submethod $line[1] for Turbomole calculation unknown\n";}
              }
            }
            close MCONTR or warn "could not close $mcontr_file\n$!\n";
            if($submethod == 0)
            {
              die "You have to define a submethod for Turbomole (td-dft/ricc-2)!\n";
            }
          }
        }
        elsif($method eq 'analytical')
        {
          $hybrid_analytical = 1;
        }
        elsif($method eq 'dftb-plus')
        {
          $hybrid_dftbp = 1;
        }
        else
        {
          warn "method $method unknown\n";
          die "currently implemented: columbus, tinker, turbomole, analytical, dftbp\n";
        }
      }
    }
    chdir(${curdir});
  }
}
#
#====================================================================================
#
sub nad_progs{
my (@line);
# List of programs for which non-adiabatic dynamics routines should be called
# for Turbomole in hybrid jobs this is indicated by $hybrid_turbomole=2
  if ( ($progname eq "columbus")
    or ($hybrid_columbus > 0)
    or ($hybrid_analytical > 0)
    or ($hybrid_turbomole > 1)
    or ($progname eq "mopac")
    or ($progname eq "analytical")
    or ($progname eq "gau")
    or ($progname eq "g09")
    or ($progname eq "aces2")
    or ($progname eq "gamess")
    or ($progname eq "turbomole") )
  {
    $nad_exec = "yes";
  }
  else
  {
    $nad_exec = "no";
  }
}
#
#====================================================================================
#
sub check_boundary{
  my ($boundary_conditions,$boundary_type);
  my ($sphere_radius,@sphere_center);
  $boundary_conditions = 0;
  $boundary_type = '';

  if ( -f "boundaries.inp")
  {
    $boundary_conditions = 1;
    print_STDOUT("\n$mdle Boundary restrictions used\n");
    open (BNDINP, "boundaries.inp") or die "$mdle Cannot open boundaries.inp\n$!\n";
    while(<BNDINP>)
    {
      chomp; s/^\s+//;
      if (/^insphere/i)
      {
        if ($boundary_type ne ''){die "$mdle Multiple boundary definitions!\n";}
        $boundary_type = "sphere";
        $_ = <BNDINP>;
        chomp; s/^\s+//;
        ($sphere_radius,@sphere_center) = split /\s+/;
        print_STDOUT("\tBoundary of type: $boundary_type\n");
        print_STDOUT("\tRadius (bohr): $boundary_type\n");
        print_STDOUT("\tCenter (bohr): @sphere_center\n\n");
      }
    }
    if ($boundary_type eq '') { die "$mdle Invalid input for boundary\n";}
  }
  else
  {
    print_STDOUT "\n$mdle No boundary restriction\n";
  }
  return $boundary_conditions,$boundary_type,$sphere_radius,@sphere_center;
}
#
#====================================================================================
#
sub check_freezing{
  my ($atomlist,$thermlist,@line,$nfrozen,$nthermfreeze);
  my ($i, $foundkey, $value, @list);
  # some restrictions in useage:
  my $vdoth  = getkeyword("sh.inp","vdoth","0");

  $nfrozen=0;
  $nthermfreeze=0;
  undef($atomlist);
  $atomlist='';
  undef($thermlist);
  $thermlist='';

# Look for the keyword to freeze carthesians in freezing_options.inp
# If present read the string and make the number sequence.
# Later check if there is also freeze.inp and die if both inputs
# are present (freeze.inp still usable for backward compatiblity).
  $foundkey = 0;
  $foundkey = &searchkeyword("freezing_options.inp","freeze_cartesians");
  if($foundkey)
  {
    open(CTD,"freezing_options.inp") or die "cannot read freezing_options.inp\n$!\n";
    while(<CTD>)
    {
      chomp; s/\s+//;
      if(/^freeze_cartesians/)
      {
        s/#*$//;
        s/^freeze_cartesians=//;
        @list=&make_num_sequence($_);
      }

      foreach $i (@list)
      {
        if(vec($atomlist,$i,1) == 0)
        {
          vec($atomlist,$i,1)=1;
          $nfrozen++;
        }
      }
    }
  }

# Look for the keyword to disable the thermostat for some atoms in
# freezing_options.inp.
# If present read the string and make the number sequence.
# Later check if there is also therm.freeze and die if both inputs
# are present (therm.freeze still usable for backward compatiblity).
  $foundkey = 0;
  $foundkey = &searchkeyword("freezing_options.inp","no_thermostate");
  if($foundkey)
  {
    open(CTD,"freezing_options.inp") or die "cannot read freezing_options.inp\n$!\n";
    while(<CTD>)
    {
      chomp; s/\s+//;
      if(/^no_thermostate/)
      {
        s/#*$//;
        s/^no_thermostate=//;
        @list=&make_num_sequence($_);
      }

      foreach $i (@list)
      {
        if(vec($atomlist,$i,1) == 0)
        {
          vec($thermlist,$i,1)=1;
        }
      }
    }
  }

  if ( -f "freeze.inp" )
  {
    $foundkey = 0;
    $foundkey = &searchkeyword("freezing_options.inp","freeze_cartesians");
    if($foundkey) { die "specify frozen atoms *either* in freezing_options.inp *or* in freeze.inp\n";}
    open(FRZINP, "freeze.inp") or die "$mdle Cannot read freeze.inp\n$!\n";
    while(<FRZINP>){
      chomp; s/^\s+//;
      @line=split/\s+/;
      for($i=0;$i<@line;$i++)
      {
        if ( $line[$i] =~ /\D/ )
        {
          die "Don't know how to use input\n$_ in freeze.inp\n";
        }
        if ( $line[$i] > $nat )
        {
          die "Cannot freeze atom $_ (only $nat in input)\n";
        }
        else
        {
          if(vec($atomlist,$line[$i],1) == 0)
          {
            vec($atomlist,$line[$i],1)=1;
            $nfrozen++;
          }
        }
      }
    }
  }
  if($nfrozen > 0)
  {
    print_STDOUT("\n$mdle Freezing of cartesians used\n");
    # look if v_init is zero for all frozen atoms
    open(VEL,"veloc") or die "\$mdle cannot read veloc\n$!\n";
    $i=1;
    while(<VEL>)
    {
      if ( vec($atomlist,$i,1) == 1 )
      {
        chomp;s/^\s+//;@line=split/\s+/;
        if( ($line[0]**2+$line[1]**2+$line[2]**2) > 0)
      	{
          die "atom number $i is chosen for freezing, but has nonzero velocity\n";
        }
      }
      $i++;
    }
    print_STDOUT("$nfrozen atoms frozen\n");
  }
  if ( -f "therm.freeze" )
  {
    $foundkey = 0;
    $foundkey = &searchkeyword("freezing_options.inp","no_thermostate");
    if($foundkey) { die "specify non-thermostat atoms *either* in freezing_options.inp *or* in therm.freeze\n";}
    print_STDOUT("\n$mdle Thermostate not used for some atoms\n");
    open(FRZINP, "therm.freeze") or die "$mdle Cannot read therm.freeze\n$!\n";
    while(<FRZINP>){
      chomp; s/^\s+//;
      @line=split/\s+/;
      for($i=0;$i<@line;$i++)
      {
        if ( $line[$i] =~ /\D/ )
      	{
          die "Don't know how to use input\n$_ in therm.freeze\n";
        }
        if ( $line[$i] > $nat )
      	{
          die "Cannot freeze atom $_ (only $nat in input)\n";
        }
        else
        {
          if(vec($atomlist,$line[$i],1) == 0)
          {
            vec($thermlist,$line[$i],1)=1;
          }
        }
      }
    }
  }
  if ( ($atomlist ne '') or ($thermlist ne ''))
  {
    for($i=1;$i<=$nat;$i++)
    {
      if((vec($atomlist,$i,1) == 1) or (vec($thermlist,$i,1) == 1))
      {
        $nthermfreeze++;
      }
    }
    open(FRI,">freeze.i") or die "\n$mdle cannot write freeze.i\n$!\n";
    print FRI "$nthermfreeze\n";
    for($i=1;$i<=$nat;$i++)
    {
      if((vec($atomlist,$i,1) == 1) or (vec($thermlist,$i,1) == 1))
      {
        printf FRI "%5d ",$i;
      }
    }
    print FRI "\n";
    close FRI;
  }
  else
  {
    open(FRI,">freeze.i") or die "\n$mdle cannot write freeze.i\n$!\n";
    print FRI "0\n0\n";
    close FRI;
    return(0,'');
  }
  return($nfrozen,$atomlist);
}
#
#====================================================================================
#
sub check_cutback{
  my($foundkey,@list,$i);
  $foundkey = 0;
  $foundkey = &searchkeyword("print_options.inp","write_pos");
  if($foundkey)
  {
    open(CTD,"print_options.inp") or die "cannot read print_options.inp\n$!\n";
    while(<CTD>)
    {
      chomp; s/\s+//;
      if(/^write_pos/)
      {
        s/#*$//;
        s/^write_pos=//;
        if( /^all$/ )
        {
          # do nothing - this is the default
        }
        else
        {
          open(CP,">cutback_pos") or die "cannot write cutback_pos\n$!\n";
          @list=&make_num_sequence($_);
          printf CP ("%6d\n",$#list+1);
          foreach $i (@list)
          {
            if($i>0){ printf CP (" %6d",$i);}
          }
          print CP "\n";
          close(CP) or warn "WARNING: could not close cutback_pos properly\n$!\n";
        }
      }
    }
    close(CTD);
  }

  $foundkey = 0;
  $foundkey = &searchkeyword("print_options.inp","no_velocity");
  if($foundkey<1)
  {
    open(WVF,">write_vel") or die "cannot write to write_vel\n$!\n";
    print WVF "1\n";
    close WVF;
  }
}
#
#====================================================================================
#
sub check_thermostat{
# check thermostat
  $thermostat = "off";
  $thermkind = "none";
  if (-s $thinp){
    $ktherm    = getkeyword($thinp,"ktherm","1");
    if ($ktherm != 0){
      $thermostat = "on";
      if ($ktherm == 1){
        $thermkind = "Andersen";
      }
      if ($ktherm == 2){
        $thermkind = "Andersen-Lowe";
      }
    }
    print_STDOUT("$mdle The thermostat is $thermostat.\n\n");
    # If default is changed here, change it also in thermostat.f90 and nxinp
    $kts      = getkeyword($thinp,"kts","1");
    $lts      = getkeyword($thinp,"lts","-1");
    $nstherm  = getkeyword($thinp,"nstherm","1");
    $temp     = getkeyword($thinp,"temp","300");
    $gamma    = getkeyword($thinp,"gamma","0.2");
    $radius   = getkeyword($thinp,"radius","10");
    $lvp      = getkeyword($thinp,"lvp","1");
    print_STDOUT("\tInitial parameters for thermostat:\n");
    print_STDOUT("\tktherm    = $ktherm ($thermkind) \n");
    print_STDOUT("\tkts       = $kts    \n");
    print_STDOUT("\tlts       = $lts    \n");
    print_STDOUT("\tnstherm   = $nstherm\n");
    print_STDOUT("\ttemp      = $temp   \n");
    print_STDOUT("\tgamma     = $gamma  \n");
    print_STDOUT("\tradius    = $radius \n");
    print_STDOUT("\tlvp       = $lvp    \n\n");
    $colprob = $dt*$gamma;
    print_STDOUT("\tCollision probability dt*gamma = $colprob \n\n");
    if ( ($kts >= $lts) and ($lts != -1) ){
        exit_error("kts cannot be larger than or equal to lts. Check $thinp.",$mdle);
    }
    if ( $temp < 0 ) {
        exit_error("Temperature cannot be negative. Check $thinp.");
    }
    if ( $gamma < 0 ) {
        exit_error("Collision frequency gamma cannot be negative. Check $thinp.");
    }
    if ( $radius < 0 ) {
        exit_error("Collision radius cannot be negative. Check $thinp.");
    }
    if ( $colprob > 1 ){
        exit_error("Collision probability dt*gamma is larger than 1. Check $thinp.",$mdle);
    }
  }
}
#
#====================================================================================
#
sub check_thirdparty{
my(@line);
# check third-party program
  my ($analytpar,$path,$am1_file,$am1_alpha,$am1_beta,$am1_kx,$am1_ky);
  my ($am1_delta,$am1_x1,$am1_x2,$am1_x3,$am1_gamma);
  my ($dftbpar,$dftb_exec,$other_state,$mult,@g,$g_vers);
  my ($cbus,$dftb,$gau,$g09,$g09r,$value_found,$value_test,@control);
  if ($progname eq "analytical") {
    print_STDOUT("$mdle Dynamics with ANALYTICAL model\n\n");
# #     This is transferred to run-analyt.pl
#     $analytpar = "analyt.par";
#     $anmod = getkeyword($analytpar,"anmod","analytical.model");
#     $path  = getkeyword($analytpar,"path",$mld);
#     if (!-x "$path/$anmod")  {
#       exit_error("File $anmod does not exist or it is not an executable.",$mdle);
#     }
#     if (($path eq $mld) && ($anmod eq "analytical.model")){
#       print_STDOUT(" 2D conical intersection model of A. Ferretti, G. Granucci, A. Lami,\n");
#       print_STDOUT(" M. Persico, and G. Villani, J. Chem. Phys. 104, 5517 (1996). \n\n");
#       $am1_file  = "$JNAD/con_int.dat";
#       $am1_alpha = getkeyword($am1_file,"alpha","3.0");
#       $am1_beta  = getkeyword($am1_file,"beta", "1.5");
#       $am1_kx    = getkeyword($am1_file,"kx",   "0.02");
#       $am1_ky    = getkeyword($am1_file,"ky",   "0.10");
#       $am1_delta = getkeyword($am1_file,"delta","0.01");
#       $am1_x1    = getkeyword($am1_file,"x1",   "4.0");
#       $am1_x2    = getkeyword($am1_file,"x2",   "3.0");
#       $am1_x3    = getkeyword($am1_file,"x3",   "3.0");
#       $am1_gamma = getkeyword($am1_file,"gamma","0.04");
#       print_STDOUT("\talpha = $am1_alpha  \n");
#       print_STDOUT("\tbeta  = $am1_beta   \n");
#       print_STDOUT("\tkx    = $am1_kx     \n");
#       print_STDOUT("\tky    = $am1_ky     \n");
#       print_STDOUT("\tdelta = $am1_delta  \n");
#       print_STDOUT("\tx1    = $am1_x1     \n");
#       print_STDOUT("\tx2    = $am1_x2     \n");
#       print_STDOUT("\tx3    = $am1_x3     \n");
#       print_STDOUT("\tgamma = $am1_gamma  \n\n");
#     }
  }
  if ($progname eq "columbus") {
    $cbus  = $ENV{"COLUMBUS"};  # Columbus environment
    if (!-e "$cbus/runc"){
      exit_error("Cannot find Columbus program. Check \$COLUMBUS variable.",$mdle);
    }
    print_STDOUT("$mdle Dynamics with COLUMBUS:\n");
    $cpar = "columbus.par";
  # if default is changed here, change also in run-col and readai_col
    $memcol     = getkeyword($cpar,"mem"       ,"200");
    $cirestart  = getkeyword($cpar,"cirestart" ,  "0");
    $reduce_tol = getkeyword($cpar,"reduce_tol",  "1");
    $quad_conv  = getkeyword($cpar,"quad_conv" , "60");
    $mc_conv    = getkeyword($cpar,"mc_conv"   ,  "0");
    $ci_conv    = getkeyword($cpar,"ci_conv"   ,  "0");
    $ivmode     = getkeyword($cpar,"ivmode"    ,  "8");
    $mocoef     = getkeyword($cpar,"mocoef"    ,  "4"); #this is needed in write_restart_info()
    $prt_mo     = getkeyword($cpar,"prt_mo"    , "20");
    print_STDOUT("\tmem        = $memcol     \n");
    print_STDOUT("\tcirestart  = $cirestart  \n");
    print_STDOUT("\treduce_tol = $reduce_tol \n");
    print_STDOUT("\tquad_conv  = $quad_conv  \n");
    print_STDOUT("\tmc_conv    = $mc_conv    \n");
    print_STDOUT("\tci_conv    = $ci_conv    \n");
    print_STDOUT("\tivmode     = $ivmode     \n");
    print_STDOUT("\tmocoef     = $mocoef     \n");
    print_STDOUT("\tprt_mo     = $prt_mo   \n\n");
    if ($reduce_tol == 2){
       print_STDOUT("\t ATTENTION! reduce_tol = 2: The energies of states whose nonadiabatic couplings\n");
       print_STDOUT("\t            are not computed will be printed but they should not be used in any\n");
       print_STDOUT("\t            quantitative analysis \n\n");
    }
  }
  if ($progname eq "turbomole") {
    print_STDOUT("$mdle Dynamics with TURBOMOLE\n");
    @control=();
    @control=read_control($JAD);
    $value_found=get_control("ricc2",@control);
    print_STDOUT("$mdle Method: $methodname\n");
    $value_test="";
    if ($methodname eq "turbomole-ricc2"){
       $value_test="cc2";
    }elsif($methodname eq "turbomole-riadc2"){
       $value_test="adc";
    }
    if ($nstat > 1){ # I added this if here to allow GS dynamics. MB May 2015
      if (($methodname eq "turbomole-ricc2") or ($methodname eq "turbomole-riadc2")){
        if ($value_found !~ /$value_test/){
          exit_error(" Method is $methodname,\nbut $value_test is not in \$ricc2 group of control. Check Turbomole input.",$mdle);
        }
      }
    }
  }
  if ($progname eq "aces2") {
    print_STDOUT("$mdle Dynamics with ACES2\n");
  }
  if ($progname eq "mopac") {
    print_STDOUT("$mdle Dynamics with MOPAC\n");
  }
  if ($progname eq "dftb") {
    print_STDOUT("$mdle Dynamics with DFTB\n");
    check_dftb($dftb_exec,$other_state,$mult);
  }
  if ($progname eq "gau") {
    $g_vers     = getkeyword("gau.par","g_vers",   "09");
    $mocoef     = getkeyword("gau.par","mocoef",    "1");
    $prt_mo     = getkeyword("gau.par","prt_mo"  , "20");
    print_STDOUT("$mdle Dynamics with GAUSSIAN G$g_vers:\n\n");
    print_STDOUT("\tg_vers     = $g_vers     \n");
    print_STDOUT("\tmocoef     = $mocoef     \n");
    print_STDOUT("\tprt_mo     = $prt_mo   \n\n");
    $gauchk="gaussian.chk";
  }
  if ($progname eq "g09") {
    $g09r  = $ENV{"g09root"};  # G09 environment
    if (!defined($g09r)){
      exit_error("\$g09root variable is not defined. Check it and run again.",$mdle);
    }  
    print_STDOUT("$mdle Dynamics with GAUSSIAN 09:\n\n");
    $mocoef     = getkeyword("g09.par","mocoef",  "0");   # the default is computing the guess orbitals at every step 
    $prt_mo     = getkeyword("g09.par","prt_mo", "20");
    $td_st      = getkeyword("g09.par","td_st", "0");
    $ld_thr     = getkeyword("g09.par","ld_thr", "14"); #linear dependence threshold for use the IOP(3/59=N) in gaussian  
    $kind_g09   = getkeyword("g09.par","kind_g09", "0");
    print_STDOUT("\tmocoef     = $mocoef   \n");
    print_STDOUT("\tprt_mo     = $prt_mo   \n");
    print_STDOUT("\ttd_st      = $td_st    \n");     
    print_STDOUT("\tld_thr     = $ld_thr   \n");
    print_STDOUT("\tkind_g09   = $kind_g09 \n");
    $g09chk="gaussian.chk";
    $g09com="gaussian.com";
    $g09rwf="gaussian.rwf";
#    inputg09();      # read the gaussian input file
  }
  if ($progname eq "tinker") {
    print_STDOUT("$mdle Dynamics with TINKER\n\n");
  }
  if ($progname eq "dftb+") {
    print_STDOUT("$mdle Dynamics with DFTB+\n\n");
#     check_dftbp($dftbp_exec,$other_state,$mult);
  }
  if ($progname eq "gamess") {
    $credits=read_credits("$progname");
    print_STDOUT("$credits");
    # if default is changed here, change also in run-gamess.pl
    $GAMESS     = $ENV{"GAMESS"};
    print_STDOUT("GAMESS path: $GAMESS\n\n");
    $gamesspar  = "gamess.par";
    $verno      = getkeyword("$gamesspar","verno","00");
    $ncpus      = getkeyword("$gamesspar","ncpus",1);
    $run_gamess = getkeyword("$gamesspar","run_gamess","rungms");
    $scr        = getkeyword("$gamesspar","scr",1);
    $mocoef     = getkeyword("$gamesspar","mocoef","1");
    print_STDOUT("$mdle Dynamics with GAMESS:\n");
    print_STDOUT("\tverno      = $verno      \n");
    print_STDOUT("\tncpus      = $ncpus      \n");
    print_STDOUT("\trun_gamess = $run_gamess \n");
    print_STDOUT("\tscr        = $scr        \n");
    print_STDOUT("\tmocoef     = $mocoef     \n");
    if ($scr == 0){
      $scrdir=find_scr_gamess($GAMESS,$run_gamess);
    }elsif ($scr == 1){
      $scrdir="$BASEDIR/$TP/SCR";
      make_gamess_scr("$scrdir");
    }
    print STDOUT "\nGAMESS SCR is $scrdir\n";
    test_gamess($GAMESS,$verno,$run_gamess,$scrdir);
#   modify GAMESS rungms file before each trajectory begins.
    $retval = callprogsp($mld,"create_rungms_nx.pl $run_gamess $scrdir",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }
  }
  if ($progname eq "hybrid") {
    print_STDOUT("$mdle Dynamics with HYBRID gradients using:\n");
    open(HC,"$JAD/hybrid.control") or die "cannot open $JAD/hybrid.control\n$!\n";
    while(<HC>){
      chomp; s/^\s+//;
      if(/^job /){
        @line=split/\s+/;
        print_STDOUT("\t job:\t $line[1]\n\t\t method:\t $line[3]\n\t\t regions:\t $line[2]\n\t\t factor:\t $line[4]\n");
      }
    }
    $cpar   = "columbus.par";
    $mocoef = getkeyword($cpar,"mocoef"    ,  "4"); #this is needed in write_restart_info()
  }
}
#
#====================================================================================
#
sub check_dftb{
  # Check DFTB inputs
  my ($dftbpar,$dftb_exec,$other_state,$mult,@g);
  my ($dftb);

  $dftb         = $ENV{"DFTB"};  # DFTB environment
  $dftbpar      = "dftb.par";
  $dftb_exec    = getkeyword($dftbpar,"dftb_exec"   ,"dftb");
  $other_state  = getkeyword($dftbpar,"other_state" ,  "1");
  $mult         = getkeyword($dftbpar,"mult"        ,  "S");
  print_STDOUT("\tdftb_exec   = $dftb_exec   \n");
  print_STDOUT("\tother_state = $other_state \n");
  print_STDOUT("\tmult        = $mult        \n\n");

  if (!-e "$dftb/$dftb_exec"){
    exit_error("Cannot find DFTB program. Check \$DFTB variable \n and the executable name (dftb_exec).",$mdle);
  }
  if (($other_state != 0) and ($other_state != 1)){
    exit_error("other_state keyword must be 0 or 1",$mdle);
  }
  if ((lc($mult) ne "s") and (lc($mult) ne "t")){
    exit_error("mult keyword must be S or T",$mdle);
  }
  if (!-s "$JAD/dftb.in"){
    exit_error("DFTB inout $JAD/dftb.in does not exist or has zero size.",$mdle);
  }
  if (!-s "$JAD/in.gen"){
    exit_error("DFTB input $JAD/in.gen does not exist or has zero size.",$mdle);
  }
  open(IG,"$JAD/in.gen") or die "Cannot read $JAD/in.gen!";
  $_=<IG>;
  close(IG);
  chomp;
  $_ =~ s/^\s*//;         # remove leading blanks
  $_ =~ s/\s*$//;         # remove trailing blanks
  @g = split(/\s+/,$_);
  if ($g[0] != $nat){
    exit_error("Inconsistent input: number of atoms \n in $JAD/in.gen and control.dyn must be the same.",$mdle);
  }
  open(IG,"$JAD/dftb.in") or die "Cannot read $JAD/dftb.in!";
  $_=<IG>;
  chomp;
  $_ =~ s/^\s*//;         # remove leading blanks
  $_ =~ s/\s*$//;         # remove trailing blanks
  @g = split(/\s+/,$_);
  if ($g[0] != 4){
    exit_error("Inconsistent input: Prepare $JAD/dftb.in \n for conjugated gardient calculation (code 4).",$mdle);
  }
  $_=<IG>;
  chomp;
  $_ =~ s/^\s*//;         # remove leading blanks
  $_ =~ s/\s*$//;         # remove trailing blanks
  # clean the quotation marks
  chop($_);
  $_=substr($_,1);
  $g[0]=$_;
  if ($g[0] ne "in.gen"){
    exit_error("Inconsistent input: Input file name \n in $JAD/dftb.in must be in.gen. Now it is $g[0].",$mdle);
  }
  $_=<IG>;
  $_=<IG>;
  $_=<IG>;
  chomp;
  $_ =~ s/^\s*//;         # remove leading blanks
  $_ =~ s/\s*$//;         # remove trailing blanks
  @g = split(/\s+/,$_);
  if ($nstatdyn >= 2){
    if ($g[0] != $nstat-1){
        exit_error("Inconsistent input: Number of excited \n states in $JAD/dftb.in must be nstat-1. Now it is $g[0].",$mdle);
    }
    if ($g[1] != $nstatdyn-1){
        exit_error("Inconsistent input: State of interest \n in $JAD/dftb.in must be nstatdyn-1. Now it is $g[1].",$mdle);
    }
    # clean the quotation marks
    chop($g[2]);
    $g[2]=substr($g[2],1);
    if (lc($g[2]) ne lc($mult)){
        exit_error("Inconsistent input: Multiplicity in $JAD/dftb.in \n and in dftb.par must be the same.",$mdle);
    }
  }elsif($nstatdyn == 1){
    if ($g[0]!~/^\'/){    # If the string does not start with simple quotation mark
        exit_error("Inconsistent input: $JAD/dftb.in must be an \n input for ground state calculation.",$mdle);
    }
  }
  while(<IG>){
      chomp;
      $_ =~ s/^\s*//;         # remove leading blanks
      $_ =~ s/\s*$//;         # remove trailing blanks
      @g = split(/\s+/,$_);
  }
  if ($g[4] != 1){       # check number of optimization steps
    exit_error("Inconsistent input: $JAD/dftb.in must be an\n input for one step optimization. Now it is $g[4].",$mdle);
  }
  close(IG);
}
#
#====================================================================================
#
sub check_dftbp
{
}
#
#====================================================================================
#
sub exec_prog{
    my ($run_flag);
    ($run_flag)=@_;

    #
    # run_flag = "initial run" => run in step 0
    #          = "first run"   => first run in each step (except 0)
    #          = "second run"  => second run in each step (after hopping)
    #
    open(WV,">which_run_is_that") or die ":( which_run_is_that";
    print WV "$run_flag\n";
    close(WV);
    if ($lvprt >= 2){print_STDOUT("$mdle $run_flag\n");}

    # Copy third-party program inputs
    $vdoth = getkeyword("sh.inp","vdoth","0");
    if (( $type == 1 ) or ($vdoth != 0)) {
      system("cp -rf $BASEDIR/$TP/$JAD/* ." );
    } else {
      system("cp -rf $BASEDIR/$TP/$JNAD/* .");
    }

    # Execute analytical model
    if ($progname eq "analytical") {                        # if PROG = 0
      $retval = callprog($mld,"run-analyt.pl",$mdle);       # Call analytical model
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    # Execute third-party program
    if ($progname eq "columbus")   {                        # if PROG = 1
      $retval = callprog($mld,"run-col.pl",$mdle);                    # call COLUMBUS
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "turbomole")  {                        # if PROG = 2
      $retval = callprog($mld,"run-turbo.pl",$mdle);                  # Call TURBBOMOLE
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "aces2")      {                        # if PROG = 3
      $retval = callprog($mld,"run-a2.pl",$mdle);                     # Call ACES2
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "mopac")      {                        # if PROG = 4
      $retval = callprog($mld,"run-mopac.pl",$mdle);                  # Call MOPAC
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "dftb")      {                         # if PROG = 5
      $retval = callprog($mld,"run-dftb.pl",$mdle);                   # Call DFTB
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "gau")      {                          # if PROG = 6
      $retval = callprog($mld,"run-gau-cas.pl",$mdle);                # Call GAUSSIAN
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "g09")      {                          # if PROG = 6.5
      $retval = callprog($mld,"run-g09.pl",$mdle);                    # Call GAUSSIAN 09
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "tinker")      {                       # if PROG = 7
      $retval = callprog($mld,"run-tinker.pl",$mdle);                 # Call TINKER
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "dftb+")      {                         # if PROG = 8
      $retval = callprog($mld,"run-dftb+.pl",$mdle);                   # Call DFTB+
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "gamess")      {                       # if PROG = 10
      $retval = callprog($mld,"run-gamess.pl",$mdle);                 # Call GAMESS
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
    if ($progname eq "hybrid")   {                          # if PROG = 20
      $retval = callprog($mld,"run-hybrid.pl",$mdle);                 # use HYBRID GRADIENTS
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
    }
}
#
#====================================================================================
#
sub sh_step_0{
# Execute surface hoping routine in step 0
  if ($jobtype eq "nonadiabatic"){

    print_STDOUT("$mdle Surface hopping with SH routine.\n\n");

    # interpolation
    $found=searchkeyword("sh.inp","ms");
    if ($found == 0)  {
      $Ms = 0;
    }
    else  {
      $Ms = keyinfile("sh.inp","ms");
    }
    if ($Ms != 0)  {
      if ($istep == 0) {
          system("cp -f epot epotint1");
      }
    }
      if (!-e "dtms")  {
          open(DTMS,">dtms") or
          exit_error("Cannot open dtms to write\nsystems answer was: $!",$mdle);
          print DTMS "$dt $dt 0 $Ms $lvprt \n";
          close(DTMS);
      }

    # run sh
    $retval = callprogsp($mld,"sh >> sh.log",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }
    #system("cp -f tprob* sh.out $BASEDIR/$RS/.");

    # print initial information
    print_STDOUT("\n$mdle Nonadiabatic dynamics information: \n");
    open(SHC,"shcard") || die "Cannot open shcard to read!";
    while(<SHC>){
      print_STDOUT($_);
    }
    print_STDOUT("\n");
    close(SHC);

  }
}
#
#====================================================================================
#
sub initial_mocoef{
# Initial molecular orbitals for columbus
  if ($progname eq "columbus" or $hybrid_columbus > 0)
  {
    if ($mocoef == 0)
    {
      # The original mocoef will be used.
      if ($lvprt>=2) {print_STDOUT("$mdle using the original mocoef\n\n");}
    }
    else
    {
      # The original mocoef will not be used from now on: move it.
      if($progname eq "columbus")
      {
        if (-e "$BASEDIR/$TP/$JAD/mocoef")
        {
          if ($lvprt>=2) {print_STDOUT("$mdle saving $BASEDIR/$TP/$JAD/mocoef as .org\n\n");}
          system("mv -f $BASEDIR/$TP/$JAD/mocoef $BASEDIR/$TP/$JAD/mocoef.org");
        }
        if (-e "$BASEDIR/$TP/$JNAD/mocoef")
        {
          if ($lvprt>=2) {print_STDOUT("$mdle saving $BASEDIR/$TP/$JNAD/mocoef as .org\n\n");}
          system("mv -f $BASEDIR/$TP/$JNAD/mocoef $BASEDIR/$TP/$JNAD/mocoef.org");
        }
      }
      else
      {
        for ($i=0;$i<@hybcoljobs;$i++)
        {
          $jobdir="$BASEDIR/$TP/$JAD/JOB_$hybcoljobs[$i].columbus";
          if (-e "$jobdir/mocoef")
          {
            if ($lvprt>=2) {print_STDOUT("$mdle saving $jobdir/mocoef.org\n\n");}
            system("mv -f $jobdir/mocoef $jobdir/mocoef.org");
          }
          $jobdir="$BASEDIR/$TP/$JNAD/JOB_$hybcoljobs[$i].columbus";
          if (-e "$jobdir/mocoef")
          {
            if ($lvprt>=2) {print_STDOUT("$mdle saving $jobdir/mocoef.org\n\n");}
            system("mv -f $jobdir/mocoef $jobdir/mocoef.org");
          }
          $jobdir="$BASEDIR/$TP/JOB_$hybcoljobs[$i].columbus";
          if (-e "$jobdir/mocoef")
          {
            if ($lvprt>=2) {print_STDOUT("$mdle saving $jobdir/mocoef.org\n\n");}
            system("mv -f $jobdir/mocoef $jobdir/mocoef.org");
          }
        }
      }# if($progname ... else ...
    } # if($mocoef ...
  } # if ($progname ...
#
# Initial molecular orbitals for Gaussian
  if ($progname eq "gau"){
     if ($mocoef == 0){  # Original chk will be used
       if ($lvprt>=2) {print_STDOUT("$mdle using the original $gauchk file.\n\n");}
     }else{
       if (-e "$BASEDIR/$TP/$JAD/$gauchk"){
          system("cp -f $BASEDIR/$TP/$JAD/$gauchk $BASEDIR/$TP/.");
          system("rm -f $BASEDIR/$TP/$JAD/$gauchk");
       }
       if (-e "$BASEDIR/$TP/$JNAD/$gauchk"){
          system("cp -f $BASEDIR/$TP/$JNAD/$gauchk $BASEDIR/$TP/.");
          system("rm -f $BASEDIR/$TP/$JNAD/$gauchk");
       }
     }
     unless (-e "$BASEDIR/$TP/$gauchk"){
       exit_error("No $gauchk file present.",$mdle);
     }
  }

  if ($progname eq "g09"){
     if ($mocoef !=0){  # chk from the original step
       if ($lvprt>=2) {print_STDOUT("$mdle using the original $g09chk file.\n\n");}
         if (-e "$BASEDIR/$TP/$JAD/$g09chk") { 
           system("cp -f $BASEDIR/$TP/$JAD/$g09chk $BASEDIR/$TP/$g09chk.ini"); 
         }
         if (-e "$BASEDIR/$TP/$JNAD/$g09chk") { 
           system("cp -f $BASEDIR/$TP/$JNAD/$g09chk $BASEDIR/$TP/$g09chk.ini");           
         }
     }else{
       if (-e "$BASEDIR/$TP/$JAD/$g09chk"){
          system("cp -f $BASEDIR/$TP/$JAD/$g09chk $BASEDIR/$TP/.");
          system("rm -f $BASEDIR/$TP/$JAD/$g09chk");
       }
       if (-e "$BASEDIR/$TP/$JNAD/$g09chk"){
          system("cp -f $BASEDIR/$TP/$JNAD/$g09chk $BASEDIR/$TP/.");
          system("rm -f $BASEDIR/$TP/$JNAD/$g09chk");
       }
     }
#     unless (-e "$BASEDIR/$TP/$g09chk"){
#       exit_error("No $g09chk file present.",$mdle);
#     }
  }
} # sub inital_mocoef
#
#====================================================================================
#
sub analysis{
# this program runs an optional arbitrary analysis program "analysis_NX" if
#    this file (or link) is present
# after that, it calls "internal_coordinate"
if (-x "$BASEDIR/NX_analysis"){
#     print_STDOUT("executing NX_analysis\n");
      $retval = callprog($BASEDIR,"NX_analysis",$mdle);       # Call external analysis program
      if ($retval != 0){
         print_STDOUT("NX_analysis failed!\n");
         print_STDOUT("If you do not want to perform this optional analysis please remove the file NX_analysis\n");
         die "$mdle is dying new (NX_analysis)";
      }
}
elsif (-e "$BASEDIR/NX_analysis"){
   print_STDOUT("Warning: NX_analysis exists but is not executable!\n");
}
#
internal_coordinate();
}
#
#====================================================================================
#
sub internal_coordinate{
  my ($ic,$icc,$ind,$intc,$retval);
  #if ($progname ne "columbus"){   # For Columbus this is done in readai_col.pl
    if (-s "intcfl"){
      if ($istep == 0){
        if ($lvprt>=2) {print_STDOUT("$mdle intcfl was found. Internal coordinates will be computed.\n");}
      }
      $ic ="$BASEDIR/$RS/intec";
      $icc="intgeom";
      $retval = callprogsp($mld,"nx2int",$mdle);
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
      open(OUT, ">>$ic" ) or die "Cannot open $ic!";
      open(INP, "$icc" ) or warn "Cannot open $icc!";
      print OUT "Time $t \n";
      print OUT "     internal coordinates\n\n";
      $ind=1;
      while(<INP>){
        chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
        $intc=$_;
        printf OUT "%4d %13.7f\n",$ind,$intc;
        $ind++;
      }
      print OUT "\n";
      close(OUT);
      close(INP);
    }else{
      if ($istep == 0){
        if ($lvprt>=2) {print_STDOUT("$mdle intcfl was not found. Internal coordinates will NOT be computed.\n");}
      }
    }
  #}
}
#
#====================================================================================
#
sub do_main_loop{
 my ($oldsurf,$newsurf,$flag,$thisstep);
 if ($lvprt>=2) {print_STDOUT("$mdle entering main loop \n");}

 while ($t < $tmax)
 {
    # Get type-of-dynamics information
    type_of_dyn();

    # write parameters:
    titlestep();

    # Keep values for later comparison
    $tprevious=$t;
    $nstatdynprevious=$nstatdyn;
# MBR Changed system for decision whether to write results or not
# "auxp" is written once *here* based on the time and not on the
# step. This prevents for discontinuities in the timeline upon restart.
# All other programs simply read auxp and act accordingly.
    open (AP,">auxp") or die "$mdle: cannot write auxp\n$!\n";
    $thisstep=sprintf("%9.2f",(${t}/${dt})+1);
    $flag=&fmodulo($thisstep,$kt);
    if($flag < $epsilon) {$flag = 0;}
    printf AP "%5.2f\n",$flag;
    close (AP);

    # Execute moldyn01: r(t+dt)
    $retval = callprogsp($mld,"moldyn01",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }

    # Execute third-party program: get energies, gradients, etc.
    exec_prog("first run");

    # freeze gradient components if requested
    freeze_gradient();

    # freeze NAD vector components of requested
    freeze_nadv();

    # Execute moldyn02: v(t+dt)
    $retval = callprogsp($mld,"moldyn02",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }

    # Freeze carthesian coordinates if requested
    freeze_veloc();

    # Execute surface hopping routine
    if ($jobtype eq "nonadiabatic") {
        interpolation($dt,$ms_def,$nstat,$istep);
    }

    # If inclusion in boundaries is chosen do it now
    apply_boundary();

    # Execute thermostat
    exec_therm();

    # Execute moldyn03: outputs
    $retval = callprogsp($mld,"moldyn03",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }

    # Read control.d: update variables
    $oldsurf=$nstatdyn;
    read_ctd();
    $newsurf=$nstatdyn;

    # Do once more if molecule hops
    # Save geometry and velocity of point of hopping in RESULTS
    if ($oldsurf != $newsurf){
      print_STDOUT("Hopping $oldsurf -> $newsurf occurred: $progname will be called again.\n");
      print_STDOUT(" === Time of hopping is $t === \n"); #GREEN
      print_STDOUT("Saving geom and veloc files to DEBUG directory as hopp_geom and hopp_veloc.\n"); #GREEN
      system("cp -f geom $BASEDIR/$DEBG/hopp_geom.$oldsurf.$newsurf.$t"); #GREEN/MBR
      system("cp -f veloc $BASEDIR/$DEBG/hopp_veloc.$oldsurf.$newsurf.$t"); #GREEN/MBR

      exec_prog("second run");
      # freeze gradient components if requested
      freeze_gradient();

      # freeze NAD vector components of requested
      freeze_nadv();
    }

    # Compute internal coordinates and perform an option analysis
      analysis();
#fp    internal_coordinate();

    # Write restart information
    write_restart_info();

    # Check whether job must continue or be Killed
    kill_job();

 } # End main loop (while)

# If tini = tmax, only execute moldyn01 and stop

 if ($tini == $tmax){
    # Get type of dynamics information
    type_of_dyn();

    # write parameters:
    titlestep();

    # Keep values for later comparison
    $tprevious=$t;
    $nstatdynprevious=$nstatdyn;

    # Execute moldyn01: r(t+dt)
    $retval = callprogsp($mld,"moldyn01",$mdle);
    if ($retval != 0){
      die "$mdle is dying now\n";
    }
 }

}
#
#====================================================================================
#
sub apply_boundary{
    if ( $boundary_conditions == 1 )
    {
      if ($boundary_type eq "sphere")
      {
        if ($lvprt >= 2)
        {
          print_STDOUT("Boundary: $boundary_type centered at @sphere_center with radius $sphere_radius\n");
        }
        open (SPHINP,">sphereinp") or die "cannot write sphereinp\n$!\n";
        printf SPHINP "%10.6f %14.8f % 14.8f % 14.8f\n",$sphere_radius,@sphere_center;
        callprogsp($mld,"inclusion_sphere",$mdle);
        if ($retval !=0){
          warn "inclusion_sphere returned code $retval\n";
          die "$mdle is dying now\n";
        }
      }
    }
}
#
#====================================================================================
#
sub freeze_gradient
{
  my($i,$k,$nextfreeze,@temp,@gradient);
  if ( $nfrozen > 0 )
  {
    system("cp grad tmpgrad");
    open(TMP,"<tmpgrad") or die "\n$mdle cannot read tmpgrad\n$!\n";
    open(GRD,">grad") or die "\n$mdle cannot write grad\n$!\n";
    print STDOUT "freeze: writing grad\n";
    if($lvprt>=2)
    {
      print STDOUT "freeze: freezing gradients of  atoms:\n";
      for($i=1;$i<=$nat;$i++)
      {
        if(vec($frozen_cart,$i,1) == 1) {printf STDOUT "%5d ",$i;}
      }
      print STDOUT "\n";
    }
    for($i=1;$i<=$nat;$i++)
    {
      $_=<TMP>;
      if ( vec($frozen_cart,$i,1) == 1 )
      {
        printf GRD "% 17.10E % 17.10E % 17.10E\n",0,0,0;
        if ($lvprt>=2)
        {
          printf STDOUT "f % 17.10E % 17.10E % 17.10E\n",0,0,0;
        }
      }
      else
      {
        chomp;s/^\s+//;s/D/E/g;@gradient=split/\s+/;
        printf GRD "% 17.10E % 17.10E % 17.10E\n",$gradient[0],$gradient[1],$gradient[2];
        if ($lvprt>=2)
        {
          printf STDOUT "  % 17.10E % 17.10E % 17.10E\n",$gradient[0],$gradient[1],$gradient[2];
        }
      }
    }
    close GRD;
    close TMP;

    if($thres > 0)
    {
      system("cp grad.all tmpgrad");
      open(TMP,"<tmpgrad") or die "\n$mdle cannot read tmpgrad\n$!\n";
      open(GRD,">grad.all") or die "\n$mdle cannot write grad.all\n$!\n";
      print STDOUT "freeze: writing grad.all\n";
      for($k=1;$k<=$nstat;$k++)
      {
        if($lvprt>=2)
        {
          print STDOUT "freeze: freezing gradients of  atoms:\n";
          for($i=1;$i<=$nat;$i++)
          {
            if(vec($frozen_cart,$i,1) == 1) {printf STDOUT "%5d ",$i;}
          }
          print STDOUT "\nstate $k\n";
        }
        for($i=1;$i<=$nat;$i++)
        {
          $_=<TMP>;
          if (vec($frozen_cart,$i,1) == 1)
          {
            printf GRD "% 17.10E % 17.10E % 17.10E\n",0,0,0;
            if ($lvprt>=2)
            {
              printf STDOUT "f % 17.10E % 17.10E % 17.10E\n",0,0,0;
            }
          }
          else
          {
            chomp;s/^\s+//;s/D/E/g;@gradient=split/\s+/;
            printf GRD "% 17.10E % 17.10E % 17.10E\n",$gradient[0],$gradient[1],$gradient[2];
            if ($lvprt>=2)
            {
              printf STDOUT "  % 17.10E % 17.10E % 17.10E\n",$gradient[0],$gradient[1],$gradient[2];
            }
          }
        }
      }
      close GRD;
      close TMP;
    }
  }
}
#
#====================================================================================
#
sub freeze_nadv
{
  my($i,$j,$k,@nadvec);
  my @nadfiles=('nad_vectors','nadv','oldh','newh');
  if ( $thres > 0 )
  {
    if ( $nfrozen > 0 )
    {
      foreach(@nadfiles){
        if( -e $_ )
        {
          system("cp $_ tmpnadv");
          open(TMP,"<tmpnadv") or die "\n$mdle cannot read tmpnadv\n$!\n";
          open(NAD,">$_") or die "\n$mdle cannot write $_\n$!\n";
          print STDOUT "freeze: writing $_\n";
          for($k=2;$k<=$nstat;$k++)
          {
            for($j=1;$j<$k;$j++)
            {
              if($lvprt>=2)
              {
                print STDOUT "freeze: freezing NAD vectors of atoms:\n";
                for($i=1;$i<=$nat;$i++)
                {
                  if(vec($frozen_cart,$i,1) == 1) {printf STDOUT "%5d ",$i;}
                }
                print STDOUT "\n  states $j $k\n";
              }
              for($i=1;$i<=$nat;$i++)
              {
                $_=<TMP>;
                if ( vec($frozen_cart,$i,1) == 1 )
                {
                  printf NAD "% 20.8F % 20.8F % 20.8F\n",0,0,0;
                  if ($lvprt>=2)
                  {
                    printf STDOUT "f % 20.8F % 20.8F % 20.8F\n",0,0,0;
                  }
                }
                else
                {
                  chomp;s/^\s+//;s/D/E/g;@nadvec=split/\s+/;
                  printf NAD "% 20.8F % 20.8F % 20.8F\n",$nadvec[0],$nadvec[1],$nadvec[2];
                  if ($lvprt>=2)
                  {
                    printf STDOUT "  % 20.8F % 20.8F % 20.8F\n",$nadvec[0],$nadvec[1],$nadvec[2];
                  }
                }
              }
            }
          }
          close NAD;
          close TMP;
        }
      }
    }
  }
}
#
#====================================================================================
#
sub freeze_veloc
{
  my($i,@vel);
  if ( $nfrozen > 0 )
  {
    system("cp veloc tmpveloc");
    open(TMP,"<tmpveloc") or die "\n$mdle cannot read tmpveloc\n$!\n";
    open(VEL,">veloc") or die "\n$mdle cannot write veloc\n$!\n";
    print STDOUT "freeze: writing veloc\n";
    if($lvprt>=2)
    {
      print STDOUT "freeze: freezing velocity vectors of atoms:\n";
      for($i=1;$i<=$nat;$i++)
      {
        if(vec($frozen_cart,$i,1) == 1) {printf STDOUT "%5d ",$i;}
      }
      print STDOUT "\n";
    }
    for($i=1;$i<=$nat;$i++)
    {
      $_=<TMP>;
      if ( vec($frozen_cart,$i,1) == 1)
      {
        printf VEL "% 24.16e% 24.16e % 24.16e\n",0,0,0;
        if ($lvprt>=2){
          printf STDOUT "f % 24.16e% 24.16e % 24.16e\n",0,0,0;
        }
      }
      else
      {
        chomp;s/^\s+//;s/D/E/g;@vel=split/\s+/;
        printf VEL "% 24.16e% 24.16e % 24.16e\n",$vel[0],$vel[1],$vel[2];
        if ($lvprt>=2)
        {
          printf STDOUT "  % 24.16e% 24.16e % 24.16e\n",$vel[0],$vel[1],$vel[2];
        }
      }
    }
    close VEL;
    close TMP;
  }
}
#
#====================================================================================
#
sub type_of_dyn{
  # DEFINE TYPE OF DYNAMICS:
  # TYPE = 1: adiabadic    - ndyn
  # TYPE = 2: nonadiabatic - ndyn-1, ndyn
  # TYPE = 3: nonadiabatic - ndyn  , ndyn+1
  # TYPE = 4: nonadiabatic - ndyn-1, ndyn  , ndyn+1
  #
  # Type of dynamics is an obsolete feature. Now it is used only for distinguish between
  # classical (type = 1) and mixed QC (type != 1) steps in some parts of the code.
  # Type_of_dynamics is, however, kept in the code because of its convenient output information.
  # The control of which states will have the non-adiabatic coupling calculated is
  # done by write_pair program based on the variables in control.dyn (thres) and in jiri.inp
  # (kross, cascade, etc.)
  #
  open(TPP,">type_prev") or die "$mdle Cannot open type_prev\nsystems answer was: $!";
  print TPP "$type";
  close(TPP);
  if ($lvprt>=2) {print_STDOUT("$mdle Type = $type \n");}
  system("rm -f chg");      # if type change, this file is recreated
  $retval = callprogsp($mld,"typeofdyn",$mdle);
  if ($retval != 0){
    die "$mdle is dying now\n";
  }
  open(TOD,"typeofdyn.out") or die "$mdle Cannot open typeofdyn.out\nsystems answer was: $!";
  $type=<TOD>;
  close(TOD);
  chomp($type);
  # keep surface information also in file
  if (-s "curr_surf"){
    system("cp -f curr_surf previous_surf");
  }
  open(CS,">curr_surf") or die "$mdle Cannot open curr_surf!";
  print CS "$nstatdyn\n";
  close(CS);
  write_transmomin();
}
#
#====================================================================================
#
sub write_transmomin{
  # write transmomin
  # Write_pair program produces for every time step a new transmomin file containg
  # all pairs of states for which the Nad. coupling should be computed.
  print_STDOUT(" write_pair:",$istep,$kt);
  $wpout="wp.out";
  $retval = callprogsp($mld,"write_pair > $wpout",$mdle);
if ($retval != 0){
  die "$mdle is dying now\n";
}
  if (-s "$wpout"){
   open(WO,"$wpout") or die "Cannot open $wpout!";
   while(<WO>){print_STDOUT($_,$istep,$kt);}
   close(WO);
   system("rm -f $wpout");
  }
  print_STDOUT("\n",$istep,$kt);
}
#
#====================================================================================
#
sub titlestep{
  if ($lvprt>=2){
    print_STDOUT("$mdle Dynamics control: \n");
    print_STDOUT("$nat $istep $nstat $nstatdyn $ndamp $kt $dt $t $tmax \n");
    print_STDOUT("$nintc $mem $nxrestart $thres $killstat $timekill $prog $lvprt\n");
    print_STDOUT("$etot_jump $etot_drift\n\n");
  }
  print_STDOUT(
"       =======================================================================\n");
  print_STDOUT(
"            FINISHING STEP $istep, TIME $t fs on SURFACE $nstatdyn            \n");
  print_STDOUT(
"       =======================================================================\n\n");
}
#
#====================================================================================
#
sub exec_therm{
    if ($thermostat eq "on"){
      $tout="thermostat.out";
      $retval = callprogsp($mld,"thermostat>>$tout",$mdle);
if ($retval != 0){
  die "$mdle is dying now\n";
}
      if (-s "$tout"){
       open(TO,"$tout") or die "Cannot open $tout!";
       while(<TO>){print_STDOUT($_);}
       close(TO);
       system("rm -f $tout");
      }
      print_STDOUT("$mdle Thermostat was executed. \n");
      system("cp -f v_therm veloc");
    }
}
#
#====================================================================================
#
sub write_restart_info{
  my ($nxrestart,$JD);
  # Write Restart information
  system("cp -f geom veloc $BASEDIR/$RT/.");
  if (-e "wfrun"){
    system("cp -f wfrun $BASEDIR/$RT/wf.inp");
  }
  if (-e "rndsrst"){
    system("cp -f rndsrst $BASEDIR/$RT/.");
  }
  if (-e "rndtm"){
    system("cp -f rndtm $BASEDIR/$RT/.");
  }
  if (-e "nad_vectors"){
    system("cp -f nad_vectors $BASEDIR/$RT/.");
  }
  if (-e "therm.log"){
    system("cp -f therm.log $BASEDIR/$DEBG/.");
  }
  #update keywords
  changekeyword("$BASEDIR/$RT/control.dyn","$BASEDIR/$RT/control.dyn","nstatdyn",$nstatdyn);
  $found=searchkeyword("$BASEDIR/$RT/control.dyn","t");
  if ($found == 0){
    addkeyword("$BASEDIR/$RT/control.dyn","$BASEDIR/$RT/control.dyn","t",$t);
  }else{
    changekeyword("$BASEDIR/$RT/control.dyn","$BASEDIR/$RT/control.dyn","t",$t);
  }
  #
  $nxrestart = 1;
  $found=searchkeyword("$BASEDIR/$RT/control.dyn","nxrestart");
  if ($found == 0){
    addkeyword("$BASEDIR/$RT/control.dyn","$BASEDIR/$RT/control.dyn","nxrestart",$nxrestart);
  }else{
    changekeyword("$BASEDIR/$RT/control.dyn","$BASEDIR/$RT/control.dyn","nxrestart",$nxrestart);
  }
  #
  open(RTI,">$BASEDIR/$RT/restart.inf")
       or warn "Cannot open $BASEDIR/$RT/restart.inf to write\n\n";
  print RTI "  NEWTON-X restart information: \n";
  print RTI "  Time (fs) = $t\n  Timestep = $istep\n  Current state = $nstatdyn\n";
  if ($progname eq "columbus" or $hybrid_columbus > 0)
  {
    if ($lvprt >=2)
    {
      print_STDOUT("writing restart info for a columbus job\n");
      print_STDOUT("progname is $progname, vdoth is $vdoth, type is $type, mocoef is $mocoef\n");
    }

    if ($progname eq "columbus")
    {
      if ($mocoef == 0)
      {
        if ($lvprt>=2) {print_STDOUT("$mdle Do not copy mocoef into $RT \n");}
      }
      else
      {
        if (-s "$BASEDIR/$RT/mocoef")
        {
          $vdoth  = getkeyword("sh.inp","vdoth","0");
          if (($type != 1) and ($vdoth == 0))
          {
            $JD=$JNAD;
          }
          else
          {
            if (-s "$BASEDIR/$RT/$JAD"){
              $JD=$JAD;
            }else{
              $JD=$JNAD;
            }
          }
          system("mv -f $BASEDIR/$RT/mocoef $BASEDIR/$RT/$JD/.");
          print RTI "  New mocoef file was copied to $JD \n";
        }
      }
    }
    else
    {
      if ($mocoef == 0)
      {
        if ($lvprt>=2) {print_STDOUT("$mdle Do not copy mocoef into $RT \n");}
      }
      else
      {
        for ($i=0;$i<@hybcoljobs;$i++)
        {
          if (-s "$BASEDIR/$RT/hybrid.$hybcoljobs[$i].columbus/mocoef")
          {
            $vdoth  = getkeyword("sh.inp","vdoth","0");
            if (($type != 1) and ($vdoth == 0))
            {
              system("cp -f $BASEDIR/$RT/hybrid.$hybcoljobs[$i].columbus/mocoef $BASEDIR/$RT/$JNAD/JOB_$hybcoljobs[$i].columbus/.");
              print RTI "  New mocoef of job $hybcoljobs[$i] was copied to $JNAD/JOB_$hybcoljobs[$i].columbus \n";
            }
            else
            {
              system("cp -f $BASEDIR/$RT/hybrid.$hybcoljobs[$i].columbus/mocoef $BASEDIR/$RT/$JAD/JOB_$hybcoljobs[$i].columbus/.");
              print RTI "  New mocoef of job $hybcoljobs[$i] was copied to $JAD/JOB_$hybcoljobs[$i].columbus \n";
            }
          }
        }
      }
    }
  }

  if ($progname eq "gau"){
    if ($mocoef == 0){
      if ($lvprt>=2) {print_STDOUT("$mdle Do not copy $gauchk into $RT \n");}
    }else{
      if (-s "$BASEDIR/$RT/$gauchk"){
        if ($type != 1) {
          system("mv -f $BASEDIR/$RT/$gauchk $BASEDIR/$RT/$JNAD/.");
          print RTI "  New $gauchk file was copied to $JNAD \n";
        }else{
          system("mv -f $BASEDIR/$RT/$gauchk $BASEDIR/$RT/$JAD/.");
          print RTI "  New $gauchk file was copied to $JAD \n";
        }
      }
    }
  }

  if ($progname eq "g09"){
    if ($mocoef != 0){
 #     if ($lvprt>=2) {print_STDOUT("$mdle Do not copy $g09chk into $RT \n");}
#    }else{
      if (-s "$BASEDIR/$JNAD/$g09chk"){
 #       if ($type != 1) {
          system("cp -f $BASEDIR/$JNAD/$g09chk $BASEDIR/$RT/$g09chk.ini");
#          system("mv -f $BASEDIR/$RT/$g09rwf $BASEDIR/$RT/$JNAD/."); 
#          print RTI "  New $g09chk file was copied to $JNAD \n";
        }elsif (-s "$BASEDIR/$JAD/$g09chk"){
          system("cp -f $BASEDIR/$JAD/$g09chk $BASEDIR/$RT/$g09chk.ini");
#          system("mv -f $BASEDIR/$RT/$g09rwf $BASEDIR/$RT/$JAD/.");
#          print RTI "  New $g09rwf file was copied to $JAD \n";
  #      }
      }
    }
  }


  if ($progname eq "gamess"){
    if ($mocoef == 0){
      if ($lvprt>=2) {print_STDOUT("$mdle Do not copy gamessinput.vec into $RT \n");}
    }else{
      if (-s "$BASEDIR/$RT/gamessinput.vec"){
        if ($type != 1) {
          system("mv -f $BASEDIR/$RT/gamessinput.vec $BASEDIR/$RT/$JNAD/.");
          print RTI "  New gamessinput.vec file was copied to $JNAD \n";
        }else{
          system("mv -f $BASEDIR/$RT/gamessinput.vec $BASEDIR/$RT/$JAD/.");
          print RTI "  New gamessinput.vec file was copied to $JAD \n";
        }
      }
      if (-s "$BASEDIR/$RT/gamessinput.nacme"){
        if ($type != 1) {
          system("mv -f $BASEDIR/$RT/gamessinput.nacme $BASEDIR/$RT/$JNAD/.");
          print RTI "  New gamessinput.nacme file was copied to $JNAD \n";
        }else{
          system("mv -f $BASEDIR/$RT/gamessinput.nacme $BASEDIR/$RT/$JAD/.");
          print RTI "  New gamessinput.nacme file was copied to $JAD \n";
        }
      }
    }
  }


  print RTI "  Do not forget of setting the new value of tmax in control.dyn. \n";
  close(RTI);
}
#
#====================================================================================
#
sub  kill_job{
# Kill job if tn+1 is equal t (it should never happen)
# Kill job if system is on state killstat more than timekill femtoseconds.
 my ($ss_result);

 # Test time variation
 $tnp=$tprevious;
 $tn=$t;
 if ($tnp == $tn) {
   exit_error("Time tn equal tn+1",$mdle);
 }

 # Test state changing
 if ($nstatdyn != $nstatdynprevious){
    open(TH,">time_hop") or die "$mdle Can't open time_hop\nsys answer was: $!";
    print TH "$t \n";
    close(TH);
 }

 # Check time after hopping
 if ($timekill != 0){
    if (-e "time_hop") {
      open(TH,"time_hop") or die "$mdle Can't open time_hop\nsys answer was: $!";
      $tlasthop = <TH>;
      close(TH);
      chomp($tlasthop);
      $timeafterhop = $t-$tlasthop;
      if (($nstatdyn == $killstat) && ($timeafterhop >= $timekill)) {
        open(DOT,$dynot) or die "$mdle Cannot open $dynot";
        print DOT "\n\t-------------------------------------------------------";
        print DOT                 "-------\n";
        print DOT "\tSystem is on $killstat state more than $timekill";
        print DOT                 " femtoseconds.\n";
        print DOT "\tFinishing dynamics.\n";
        print DOT "\t---------------------------------------------------------";
        print DOT                 "-----\n";
        close DOT;
        last;
      }
    }
 }

# Run and check stop-sign
 if (-s "stopsign.inp"){
   $retval = callprogsp($mld,"stop-sign.pl",$mdle);
   if ($retval != 0){
     die "$mdle is dying now\n";
   }
   $ss_result="";
   open(SS,"ss-result") or warn "$mdle Cannot open ss-result.";
   while(<SS>){
     $ss_result=$ss_result.$_;
   }
   close(SS);
   print_STDOUT("\n$mdle STOP Conditions:\n$ss_result\n\n");
   if ($ss_result =~ /%STOP/){
     print_STDOUT("$mdle STOP Conditions were satisfied. Trajectory will be terminated.\n\n");
     last;
   }
 }
}
#
#====================================================================================
#
sub interpolation($$$$){
  my ($dt,$ms_def,$nstat,$istep)=@_;
  my ($Ms,$Delta_ts,$j,$nstatdyn,$oldsurf,$newsurf,$idum,$rec_veloc,
      $irk,$mom,$adjmom,$surf_in_first_point);
  sub interpolate_energ_1($$);
  sub vel_after_hop($$$$$$$);

 # read control.d
  open(CT, $ctd) or die "$mdle Cannot open $ctd!\nsystems answer was: $!\n";
  $_=<CT>;
  chomp;s/\s+//g;
  @prmt=split /,/;
  close(CT);
  ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc,
  $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,$etot_drift)=@prmt;

  $dt=shift(@_);

 # Read MS in sh.inp
  $found=searchkeyword("sh.inp","ms");
  if ($found == 0)  {
    $Ms = $ms_def;
  }  else  {
    $Ms = keyinfile("sh.inp","ms");
  }

 # backup files
   system("cp -f nad_vectors nad_vectors.bk");
   system("cp -f epot epot.bk");
   system("cp -f veloc veloc.bk");

  $rec_veloc = 0;
  if ($Ms == 0){                                               # do NOT interpolate

    print_STDOUT("\n$mdle Integration of the TDSE \n");
    $retval = callprogsp($mld,"sh >> sh.log",$mdle);           # surface hopping routine
    if ($retval != 0){
      die "$mdle is dying now\n";
    }

  } else {                                                     # interpolate

    if ($lvprt>=2){print_STDOUT("\n$mdle Substeps in the integration of the TDSE: \n");}
    interpolate_energ_1($nstat,$istep);                        # prepare data for interpolation

    $surf_in_first_point = $nstatdyn;                          # current surface in the 1st substep

    # Main interpolation loop
    for ($j=1; $j<=$Ms; $j++) {

      $Delta_ts = $j*$dt/$Ms;
      open(MS,">dtms") or die"Cannot open dtms to write!!\nsystems answer was: $!";
      print MS     " $Delta_ts $dt $j $Ms $lvprt \n";
      if ($lvprt>=2){print_STDOUT("$mdle MS=$Ms dt=$dt Delta_ts=$Delta_ts \n");}
      close(MS);

      $retval = callprogsp($mld,"interpol",$mdle);             # interpolate v and h
      if ($retval != 0){
        die "$mdle is dying now\n";
      }
      $retval = callprogsp($mld,"polint >> polint.log",$mdle); # interpolate E
      if ($retval != 0){
        die "$mdle is dying now\n";
      }

      open(CT, $ctd) or die "$mdle Cannot open $ctd!\nsystems answer was: $!";
      $_=<CT>;
      chomp;s/\s+//g;
      @prmt=split /,/;
      close(CT);
      $oldsurf=$prmt[3];                                       # store current state before surface hoppinp routine

      $retval = callprogsp($mld,"sh >> sh.log",$mdle);         # surface hopping routine
      if ($retval != 0){
        die "$mdle is dying now\n";
      }

      open(CT, $ctd) or die "$mdle Cannot open $ctd!\nsystems answer was: $!";
      $_=<CT>;
      chomp;s/\s+//g;
      @prmt=split /,/;
      close(CT);
      $newsurf=$prmt[3];                                       # update the information about the current surface

      # update the information about hopping event
      open(IK, "irk") or die "$mdle Cannot open irk!\nsystems answer was: $!";
      $_=<IK>;
      chomp;s/^\s+//;
      @prmt = split /\s+/;
      close(IK);
      ($irk,$mom,$adjmom) = @prmt;

      # in the last sub-timestep check what happened
      if ($j == $Ms){
        if ($surf_in_first_point != $newsurf){                 # check change of state
          $irk = 1;
          $oldsurf = $surf_in_first_point;
          $rec_veloc = 1;
          if ($lvprt>=2){print_STDOUT("$mdle Adjusting velocity.\n");}
        } else {
          $irk = 0;
        }
      }

      if ($irk > 0.5){                                         # If state change do the following
         vel_after_hop($Ms,$j,$irk,$mom,$adjmom,$newsurf,$oldsurf);

         open(IPL,">>interpol.log") or die "Cannot open interpol.log to append.";
         if ($irk == 1){
           #system("cp -f veloc oldveloc");  # after surface hopping, v is kept constant
           #system("cp -f veloc newveloc");  # until to start the new step.
           #print IPL "\n Surface Hopping: veloc won't change in the next subtimesteps.\n\n";
           #$rec_veloc = 1;
         }
         # This option is not working with interpolation!!!!
         if ($irk == 2){
           if ($mom == 1){
              print IPL "\n Frustrated Hopping. Keep the velocity. \n\n";
           }
           if ($mom == -1){
              print IPL "\n Frustrated Hopping. Invert the velocity. \n\n";
           }
           $rec_veloc = 2;
         }
         close(IPL);

      } # endif ($irk > 0.5)

      if ($lvprt >= 2) {system("cp -f interpol*.log polint* $BASEDIR/$DEBG/.");}

    } # End of main interpolation loop

  }

  # Recover backuped files after interpolation
  if ($Ms != 0){
    system("cp -f epot.bk epot");
    if ($rec_veloc == 0) {
      system("cp -f veloc.bk veloc");
      if ($lvprt>=2){print_STDOUT("$mdle Recovering veloc.\n");}
    }
    if ($rec_veloc == 1){
      system("cp -f veloc.adj veloc");
      if ($lvprt>=2){print_STDOUT("$mdle Copying adjusted veloc.\n");}
    }
    if ($rec_veloc == 2){                                      # In case of frustrated hoppoing, keep velocity
       if ($mom == 1){
          system("cp -f veloc.bk veloc");
       }
    }
  }

}
#
#====================================================================================
#
sub vel_after_hop($$$$$$$){
# After hopping (irk = 1) the velocities in the new surface (newsurf) are
# adjusted by sh program. This routine takes these velocities and projects
# them to the first and last sub-timestep of the current timestep, using
# the gradient of newsurf.
#   veloc_first = veloc - acceleration(newsurf)*j*dt/Ms
#   veloc_last  = veloc + acceleration(newsurf)*(Ms-j)*dt/Ms
# These projections are further used as oldveloc and newveloc to get the
# interpolated values. Note that if the grad(newsuf) is not available,
# it is assumed to be 0, and the velocity after hopping is kept constant up
# to the end of the timestep.
# These interpolated values are used only during the timestep.
# At the last subtimestep the original velocity is recovered and depending
# whether there is a change of surface or not (irk = 1 or 0), the original
# velocity is adjusted to conserve the energy.
# All these procedures are controlled by moldyn04.
# MB, Dec-2006.

   my ($Ms,$j,$irk,$mom,$adjmom,$newsurf,$oldsurf)=@_;

   #print_STDOUT("\n\n %%%%%% $Ms,$j,$irk,$mom,$newsurf  %%%%%%\n \n ");
   open(MD4,">md04.inp") || die "Cannot open md04.inp to write!";
   print MD4 "$Ms $j $irk $mom $adjmom $newsurf $oldsurf";
   close(MD4);

   $retval = callprogsp($mld,"moldyn04",$mdle);
if ($retval != 0){
  die "$mdle is dying now\n";
}

}
#
#====================================================================================
#
sub interpolate_energ_1($$){
   my ($nstat,$istep)=@_;
   my ($ns,@epot,$par,@newepot,@oldepot,$vx,$vy,$vz,$mute);

   if ($istep == 0) {
      system("cp -f epot epotint1");
   }
   if ($istep == 1) {
      system("cp -f oldepot epotint1");
      system("cp -f newepot epotint2");
   }
   if ($istep == 2) {
      system("cp -f newepot epotint3");
   }
   if ($istep == 3) {
      system("cp -f newepot epotint4");
   }
   if ($istep >= 4) {
      system("cp -f epotint2 epotint1");
      system("cp -f epotint3 epotint2");
      system("cp -f epotint4 epotint3");
      system("cp -f newepot  epotint4");
   }

 if ($lvprt >= 2) {
   if ($istep >= 1) {

        open(IPL,">>interpol_e.log") or die "Cannot open interpol.log to append.";

                print IPL "\n ........................................ ";
                print IPL "\n ..............  EPOT 1 ................. \n \n";

        open(NE,"epotint1") or die "Cannot open epotint1 to read!";
        $ns = 1;
        while(<NE>){
                chomp;
                $newepot[$ns]=$_;
                printf IPL ("%15.6F \n",$newepot[$ns]);
                $ns++;
        }
        close(NE);


        print IPL "\n ..............  EPOT 2 ................. \n \n";
        open(NE,"epotint2") or die "Cannot open epotint2 to read!";
        $ns = 1;
        while(<NE>){
                chomp;
                $newepot[$ns]=$_;
                printf IPL ("%15.6F \n",$newepot[$ns]);
                $ns++;
        }
        close(NE);

        if ($istep >= 2) {
        print IPL "\n ..............  EPOT 3 ................. \n \n";
        open(NE,"epotint3") or die "Cannot open epotint3 to read!";
        $ns = 1;
        while(<NE>){
                chomp;
                $newepot[$ns]=$_;
                printf IPL ("%15.6F \n",$newepot[$ns]);
                $ns++;
        }
        close(NE);

        if ($istep >= 3) {
        print IPL "\n ..............  EPOT 4 ................. \n \n";
        open(NE,"epotint4") or die "Cannot open epotint4 to read!";
        $ns = 1;
        while(<NE>){
                chomp;
                $newepot[$ns]=$_;
                printf IPL ("%15.6F \n",$newepot[$ns]);
                $ns++;
        }
        close(NE);

        }
        }
        print IPL "\n .... STARTING INTERPOLATION OF E ....... \n \n";
        close(IPL);

   }

# Write information about velocities

   open(OV,"oldveloc")       || die "Cannot open oldveloc to read.";
   open(IPL,">>interpol.log") || die "Cannot open interpol.log to append.";
   print IPL "\n ........................................ ";
   print IPL "\n .............. PREVIOUS V ............... \n \n";
           while (<OV>) {
              chomp;
              ($mute,$vx, $vy, $vz)=split(/\s+/,$_);
              printf IPL ( "%14.7F %14.7F %14.7F\n",$vx,$vy,$vz);
           }
   close(OV);

   open(NV,"newveloc")       || die "Cannot open newveloc to read.";
   print IPL "\n ........................................ ";
   print IPL "\n .............. CURRENT V ............... \n \n";
           while (<NV>) {
              chomp;
              ($mute,$vx, $vy, $vz)=split(/\s+/,$_);
              printf IPL ( "%14.7F %14.7F %14.7F\n",$vx,$vy,$vz);
           }
   close(NV);
   print IPL "\n .. STARTING INTERPOLATION OF V and H ... \n \n";
   close(IPL);

 }
}
#====================================================================================
sub fmodulo
{
  # takes two numbers and returns the modulo it does this also for negative
  # numbers and non-integers
  my($a,$b,$f,$sign,$epsilon);
  $epsilon = 1E-9;
  $a=shift(@_);
  $b=shift(@_);
#   printf STDOUT "DEBUGMR: function fmoldulo: % 30.20f, %29.20f\n", $a, $b;
  $f=0;
  if($a>=0){$sign=+1;}
  else{$sign=-1;}
  if($b<=0){die "only modulo a positive number\n $a % $b impossible\n";}
  # due to numeric problems the >= has to be done manually with an epsilon
  while((abs($a-($sign*$b))<$epsilon) or (abs($a)>$b))
  {
    $a -= ($b*$sign);
    $f++;
  }
#   printf STDOUT "DEBUGMR: result % 30.20f\n", $a;
  return $a;
}

#====================================================================================
#                              END OF SUB ROUTINES
#====================================================================================
