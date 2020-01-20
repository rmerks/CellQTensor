#!/usr/bin/perl -w  +



$shape_number   = 9;

$anisotropy     = 2;


$num_of_adhesions = 4;
$num_of_bounpoint_per_line = $anisotropy*20;

$num_of_gridpoints_on_line_in_box = 20;


$num_of_gridpoints_on_line_exp = 512; #Must be a whole multiple of num_of_gridpoints_on_line_in_box-1
$num_of_iteration = 550000;
$period = 50000;

$friction = 1.0;
$rotational_friction = 1.0;
$elasticity = 1.0; #always 1.0, only the ratio of anchoring/elasticity matters.


$phase_field = 1.0;
$defect_core = 0.15;

$rotation=0.0;

$q_cutoff = 0.10; # This is the value used in the equilibrium function to determine if the system is in equilibrium.
$lowordercutoff = 0.20;

$timestepbulkratio = 0.001;
$time_stepratio = 0.002;
$spacingupdateratio = 0.1;



$seed = -1; #Watch out! Only seed -1 gives 'random' start (different each run).



$directory = ""; #Insert the directory of the c-code (cyto.c) here.

$lambda = 1.0;
$anchoring = 10.0;

$alpha = 0.0;
$sigma = 0.5;


while ($alpha<1.5){

#Use a while/for loop to do a parametersweep.

system("mkdir /$directory/shape$shape_number");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha");
system("mkdir /$directory/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation");


$i = 1;
while ($i<3){

open (FILE, '>input.AN.2019.08.13.dat');  #simulation parameters

print FILE "$shape_number\n";
print FILE "$anisotropy\n";
print FILE "$num_of_adhesions\n";
print FILE "$num_of_bounpoint_per_line\n";
print FILE "$num_of_gridpoints_on_line_in_box\n";
print FILE "$num_of_gridpoints_on_line_exp\n";
print FILE "$num_of_iteration\n";
print FILE "$period\n";

print FILE "$alpha\n";
print FILE "$sigma\n";
print FILE "$lambda\n";
print FILE "$friction\n";
print FILE "$rotational_friction\n";
print FILE "$elasticity\n";
print FILE "$anchoring\n";
print FILE "$phase_field\n";
print FILE "$defect_core\n";
print FILE "$rotation\n";

print FILE "$q_cutoff\n";
print FILE "$lowordercutoff\n";
print FILE "$noise\n";
print FILE "$timestepbulkratio\n";
print FILE "$time_stepratio\n";
print FILE "$spacingupdateratio\n";

print FILE "$simulation_exp\n";

print FILE "$seed\n";

close (FILE);

system("make");
system("cyto >Prints.dat");
system("cat Prints.dat");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv bulk_*.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv vertex_*.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv force_*.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv global_*.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv Prints.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");
system("mv input.AN.2019.08.13.dat /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation/run$i");

$i=$i+1;
}

$alpha=$alpha+0.5;
$sigma=$sigma-0.25;
}

