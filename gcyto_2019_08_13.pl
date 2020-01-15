#!/usr/bin/perl -w  +



$shape_number   = 9;

$anisotropy     = 2;


#while ($shape_number<106){
#while ($anisotropy<2.5){

$num_of_adhesions = 4;
$num_of_bounpoint_per_line = $anisotropy*20;

$num_of_gridpoints_on_line_in_box = 20;
#while ($num_of_gridpoints_on_line_in_box<100){

$num_of_gridpoints_on_line_exp = 512; #Must be a whole multiple of num_of_gridpoints_on_line_in_box-1
$num_of_iteration = 550000;
$period = 50000;    

$simulation_exp = 0; 
#Tells us if we are comparing the simulations to experimental data. 
#$simulation_exp = 0; no comparison, using alpha & sigma values put in here.
#$simulation_exp = 1; comparison, using alpha & sigma values calculated in cyto.c from experimental data.
#$simulation_exp = 2; no comparison, using alpha & sigma values calculated in cyto.c from experimental data.


#while($alpha<2.5){

$friction = 1.0;
$rotational_friction = 1.0;
$elasticity = 1.0; #always 1.0, only the ratio of anchoring/elasticity matters.

#while($anchoring<41.0){

$phase_field = 1.0;
$defect_core = 0.15;

#while($defect_core<2.1){



$rotation=0.0;

#while ($rotation<2.1) {

$q_cutoff = 0.10; # This is the value used in the equilibrium function to determine if the sytem is in equilibrium.
$lowordercutoff = 0.20;
$noise = 0.05; # Is currently not being used.

$timestepbulkratio = 0.001;
$time_stepratio = 0.002;
$spacingupdateratio = 0.1;



$seed = -1; #Watch out! Only seed -1 gives 'random' start (different each run).



$version = "Bulk_8v1";
$directory = "data1/koen/cell_mechanics/code_continuum/jeremy"; 

$lambda = 1.0;
$anchoring = 10.0;


#while ($anchoring<11){

$alpha = 0.0; 
$sigma = 0.5;


while ($alpha<1.5){

#Use a while/for loop to do a parametersweep.

if ($simulation_exp == 1) {
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring");
}

else {
system("mkdir /$directory/$version/shape$shape_number");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma");
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha"); 
system("mkdir /$directory/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/anisotropy$anisotropy/anch$anchoring/sigma$sigma/alpha$alpha/rotation$rotation");
}

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



if ($simulation_exp == 1) {
system("make");
system("cyto >Prints.dat");
system("cat Prints.dat");
system("mkdir /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv bulk_*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv vertex_*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv force_*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv global_*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv Prints.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv input.AN.2018.12.05.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv exp_b*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
system("mv exp_comparison*.dat /data1/Jeremy/Code/$version/shape$shape_number/iterations$num_of_iteration/pointsonline$num_of_gridpoints_on_line_in_box/defectcore$defect_core/experimental/anch$anchoring/run$i");
}

else {
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
}

$i=$i+1;
}

$alpha=$alpha+0.5;
$sigma=$sigma-0.25;
}

#if ($anchoring <0.5){
#$anchoring=$anchoring+1.0;
#}
#else {
#$anchoring =$anchoring+9.0;
#}

#}

#$sigma = $sigma -0.5;
#$alpha = $alpha+1.0;

#$sigma = $sigma +0.125;
#$alpha = $alpha+1.0;
#}



#$anchoring=$anchoring+1.0;
#}

#$alpha=$alpha+1.0;
#$sigma=$sigma+0.125;
#}

#$shape_number=$shape_number+1;
#}
#$anisotropy=$anisotropy+0.5;
#}
#$num_of_gridpoints_on_line_in_box=$num_of_gridpoints_on_line_in_box+10;
#}





# The following variables were put into the C file, as they will not change per simulation.
#$num_of_pixels = 512;
#$size_of_pixel = 0.138/1000000; #Divide bij 10^6 because initial value of 0.138 was n um
#$gamma_data = 0.40; #Is dimensionless
#$sigmaratio = 14.7/1000000; #Divide bij 10^6 because initial value of 14.7 was in um

