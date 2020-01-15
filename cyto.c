/******************************************************************
**                                                               **
**                              Cell                             **
**             							 **
**                                                               ** 
** Program: cyto.c                                               **
** Version: Bulk_7v3                                      **
** Last Change: 2018/12/05                                       **
**								 **
*******************************************************************/
    /* 
     * Uses S in calculations of border. Updated jackknife with correct forefactor of (n-1)/(n)
     */

#ifdef SINGLE
#define REAL float
#else 
#define REAL double
#endif

#include "cyto.h"

#define delta(x,y) ((x == y ? 1 : 0))
#define mod(x,y) ((x % y + y) % y)
#define max(x,y) ((x > y) ? x : y)
#define sign(x) (x >= 0 ? 1 : -1 )

Vertex *vertex;
Grid *grid;
Grid_exp *grid_exp;
Grid_old *grid_old;
Lines *lines;
Forces *forces;
Border_exp *border_exp;

long period;
long shape_number;
long num_of_bounpoint;
long num_of_bounpoint_per_line;
long num_of_edgepoint;
long num_of_iteration;
long num_of_adhesions;
long num_of_adhesions_custom;
long num_of_gridpoints;
long num_of_gridpoints_on_line;
long size_of_margin = 0;
long num_of_gridpoints_on_line_in_box;
long num_of_gridpoints_on_line_x;

long num_of_lines;
long seed;
long neighbors[4];
long num_of_pixels;

long num_of_gridpoints_on_line_exp;
long num_of_gridpoints_exp;

long num_of_borderpoints;

long simulation_exp; //Tells us if we are comparing the simulations to experimental data. Has value 1 if so.


double sigma;
double alpha;
double lambda;
double friction;
double rotational_friction;
double elasticity;
double anchoring;
double noise;
double pi = 3.14159265359;
double anisotropy;
double lattice_spacing;
double adhesions[30][2];
double timestepbulk;
double timestepbulkratio;
double time_step;
double time_stepratio; //If this number is too large, this might lead to  "all surrounding pixels were not in the cell before."
double spacingupdateratio;
double order_par_global;
double wall_order_par;
double total_force;
double total_force_x;
double total_force_y;
double average_director;
double rotation;
double size_of_pixel;
double gamma_data;
double sigmaratio_initial;
double sigmaratio_scaled;
double lowordercutoff;
double q_cutoff;
int shape_difference_total=0;
int q_difference_total=0;

//Included in the Qxx and Qxy code:
double phase_field;
double defect_core;
int p = 0;
int exporting;
long iteration;


double ran2(long *);

time_t get_cpu_time();

void end();
void run();
void init();
void find_orientation();
void update_boundary();
void find_adhesion_forces();
void determine_bulkedge();
void determine_bulk();
void update_bulk();
void init_boundary();
void init_bulk(); 

void allocate_memory();
void export_conf();
//void progress_bar(long);
void get_time(CPU_Time *);

void compare_data();
void equilibrium();

void init_input_exp();
void export_conf_exp();
void comparison_exp();

long bulkpoints_sim=0, bulkpoints_exp=0, bulkpoints_shared=0;
double deltasquaredphi=0.0,deltasquaredq1=0.0,deltasquaredq2=0.0,deltasquaredq=0.0,deltasquaredangle=0.0,deltasquaredanglecom=0.0,deltasquaredanglecom2=0.0;
double deltasquaredjkq1,deltasquaredjkq2,deltasquaredjkq,errordeltasquared=0.0;
// Have to be defined here because they are needed in both comparison_exp and export_conf_exp.

int main()
{
        (void) get_cpu_time();
        
	init();
	run();
    if(simulation_exp==1){
        compare_data();
    }
	end();

  	return EXIT_SUCCESS;
}

/*******************************************************************/

void init()
{
        int i,minnumofgridpointsonline;
        double x0,y0,deltax;


        
        FILE *myFile;
        myFile =fopen("input.AN.2019.08.13.dat", "r");
        
        //read file to find the parameter values
        fscanf(myFile, "%ld\n", &shape_number);
        fscanf(myFile, "%lf\n", &anisotropy);
        fscanf(myFile, "%ld\n", &num_of_adhesions);
        fscanf(myFile, "%ld\n", &num_of_bounpoint_per_line);
        fscanf(myFile,"%ld\n",&num_of_gridpoints_on_line_in_box);
        fscanf(myFile,"%ld\n",&num_of_gridpoints_on_line_exp);
        fscanf(myFile,"%ld\n",&num_of_iteration);
        fscanf(myFile,"%ld\n",&period);
        
        fscanf(myFile,"%lf\n",&alpha);
        fscanf(myFile,"%lf\n",&sigma);
        fscanf(myFile,"%lf\n",&lambda);
        fscanf(myFile,"%lf\n",&friction);
        fscanf(myFile,"%lf\n",&rotational_friction);
        fscanf(myFile,"%lf\n",&elasticity);
        fscanf(myFile,"%lf\n",&anchoring);
        fscanf(myFile,"%lf\n",&phase_field);
        fscanf(myFile,"%lf\n",&defect_core);
        fscanf(myFile,"%lf\n",&rotation);
        
        fscanf(myFile,"%lf\n",&q_cutoff);
        fscanf(myFile,"%lf\n",&lowordercutoff);
        fscanf(myFile,"%lf\n",&noise);
        fscanf(myFile,"%lf\n",&timestepbulkratio);
        fscanf(myFile,"%lf\n",&time_stepratio);
        fscanf(myFile,"%lf\n",&spacingupdateratio); 
        
        fscanf(myFile,"%ld\n",&simulation_exp); 
        
        fscanf(myFile,"%ld\n",&seed);

        
        if(seed<0){                     //If I want the program to take a seed from the clock, then I'll make seed in the pl file negative.
            seed = time(NULL);
        }
     
        
       
        //Now depending on the shape number, I will overrule the value for the num_of_adhesions:
        if(shape_number==1){
            num_of_adhesions = 4;
            //size_of_margin = 3;
        }
        else if(shape_number==3){
            num_of_adhesions = 13;
        }
        else if(shape_number==2){                                         //If this option is chosen, it is important to define the adhesion points in clockwise manner (otherwise things might go wrong with minus sign later)            
            for(i=0;i<num_of_adhesions;i++){
                printf("Input coordinates of adhesion point ");
                printf("%i",i);
                printf(": ");
                (void) scanf("%lf\t%lf",&adhesions[i][0],&adhesions[i][1]);
            }
        }
        else if(shape_number==8){
            num_of_adhesions = 4;
        }
        else if(shape_number==9){
            num_of_adhesions = 4;
        }
        else if(shape_number==10){
            num_of_adhesions = 4;
        }
        else if(shape_number==11){
            num_of_adhesions = 3;
        }
        
        
        else if(shape_number >= 100){
            lambda=1.0;
            num_of_pixels = 512;
            size_of_pixel = 0.138; // This is in um

            
            if(shape_number == 100){
            num_of_adhesions = 21;
            gamma_data= 0.40;
            sigmaratio_initial=14.7; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
            else if(shape_number==101){
            num_of_adhesions = 25;
            gamma_data= 0.95;
            sigmaratio_initial=10.8; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
            else if(shape_number==102){
            num_of_adhesions = 23;
            gamma_data= 0.46;
            sigmaratio_initial=18.0; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
            else if(shape_number==103){
            num_of_adhesions = 15;
            gamma_data= 0.52;
            sigmaratio_initial=13.4; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
            else if(shape_number==104){
            num_of_adhesions = 19;
            gamma_data= 0.25;
            sigmaratio_initial=15.7; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
            else if(shape_number==105){
            num_of_adhesions = 8;
            gamma_data= 0.75;
            sigmaratio_initial=12.6; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
                        else if(shape_number==200){
            num_of_adhesions = 19;
            gamma_data= 0.25;
            sigmaratio_initial=15.7; //This is in um. This is not a nice notation, have to find better way to write power law 10^-6
            }
                       
        }

            
        
        num_of_bounpoint = num_of_adhesions*num_of_bounpoint_per_line;
        num_of_edgepoint = num_of_bounpoint_per_line;
             
       //I calculate the minimal amount of grid points on line you need in order for the bulk not to get messy, and print a statement if the chosen one was smaller than the minimum.
        if(anchoring!= 0){
            deltax = elasticity/(4*anchoring)*spacingupdateratio*0.5*pi; //Deltax needs to be smaller than about K/4W (elasticity/(4*anchoring))*0.1*0.5Pi. Otherwise the updates on the boundary on phi are too large.
            minnumofgridpointsonline = floor(1.0/deltax)+1;
        }
        if(anchoring !=0 && num_of_gridpoints_on_line_in_box<minnumofgridpointsonline){
            printf("%s%i%s%i%s\n","The number of gridpoints on a line (",num_of_gridpoints_on_line_in_box,") is smaller than the minimum required (",minnumofgridpointsonline,")");
        }
        
        
        //Now use the obtained parameters to calculate some extra parameters:
        num_of_gridpoints_on_line = num_of_gridpoints_on_line_in_box +2*size_of_margin;
        num_of_gridpoints_on_line_x = floor((num_of_gridpoints_on_line_in_box+0.0)*anisotropy)+2*size_of_margin; //the floor is because of troubles with integers and doubles.
        num_of_gridpoints = (num_of_gridpoints_on_line)*(num_of_gridpoints_on_line_x); //I need more gridpoints in the x-direction
	num_of_lines = (num_of_gridpoints_on_line_x)*(num_of_gridpoints_on_line-1);          //So there are also more vertical lines. this is the number of vertical (or horizontal lines), not all lines. However, for now I'll only use the vertical lines.
        lattice_spacing = 1.0/(num_of_gridpoints_on_line-2*size_of_margin-1.0);
        num_of_gridpoints_exp = num_of_gridpoints_on_line_exp*num_of_gridpoints_on_line_exp;
        
        neighbors[0] = 1;
        neighbors[1] = -1;
        neighbors[2] = num_of_gridpoints_on_line;
        neighbors[3] = -num_of_gridpoints_on_line;
        
        
        


	allocate_memory();
	init_boundary();
    find_adhesion_forces();
    init_bulk(); //Define a grid on which the cell lives. On this grid we will solve the equation for the director/Q-tensor.
    determine_bulkedge();
    determine_bulk();
    
    
  
        
}

/*******************************************************************/

void run()
{
	long i;
	
    iteration=0;
    
	for (i=0; i<=num_of_iteration; i++){

                if (i%period==0){
                    export_conf();
                    exporting=1;
                    equilibrium();
                    //export_conf_exp();

                }
                else{
                    exporting=0;
                }
                //progress_bar(i);
                //update_bulk();
                find_orientation();
                update_boundary();
                find_adhesion_forces();
                determine_bulkedge();
                determine_bulk();
                update_bulk();
                
                iteration++;
	}
	        printf("\nalpha: %.10f\n", alpha);
        printf("\nsigma: %.10f\n", sigma);
        printf("\nlambda: %.10f\n", lambda);
}

/*******************************************************************/

void allocate_memory()
{
	long i,num_of_lines_double;
        
        num_of_lines_double = 2* num_of_lines + num_of_gridpoints_on_line_x; //The last term is because there are now more horizontal lines than vertical lines.
	
        vertex = (Vertex *) malloc((num_of_bounpoint+100)*sizeof(Vertex));
        grid = (Grid *) malloc(num_of_gridpoints*sizeof(Grid));
        grid_exp = (Grid_exp *) malloc(num_of_gridpoints_exp*sizeof(Grid_exp));
        grid_old = (Grid_old *) malloc(num_of_gridpoints_exp*sizeof(Grid_old));
        border_exp = (Border_exp *) malloc(num_of_gridpoints_exp*sizeof(Border_exp));
       
        lines = (Lines *) malloc(num_of_lines_double*sizeof(Lines)); //In this version I will also store horizontal lines that cross the boundary. Although I do not need those to determine what is the bulk, these are useful for storing the local angle of the boundary.
        forces = (Forces *) malloc(num_of_adhesions*2*sizeof(Forces));

}

/*******************************************************************/

void compare_data()
{

        
        init_input_exp();
        comparison_exp();
        export_conf_exp();

           
}

void init_input_exp()
{

    
    long i,j,k,l,m,n,o;
    
    //long num_of_gridpoints_on_line_exp, num_of_gridpoints_on_line_exp_total;
    long startpoint, counts, pix_per_point;
    double sumq1, sumq2;
    long num_of_gridpoints_noborder;
    char inputfilename_orientation[64];
    char inputfilename_bulk[64];
    char inputfilename_border[64];
    char inputfilename_coherency[64];
    

    if(shape_number == 100){
        sprintf(inputfilename_orientation,"shape100cell253orientation.txt");
        sprintf(inputfilename_bulk,"shape100cell253bulk.txt");
        sprintf(inputfilename_border,"shape100cell253border.txt");
        sprintf(inputfilename_coherency,"shape100cell253coherencyv3.txt");
    }
    
    if(shape_number == 101){
        sprintf(inputfilename_orientation,"shape101cell258orientation.txt");
        sprintf(inputfilename_bulk,"shape101cell258bulk.txt");
        sprintf(inputfilename_border,"shape101cell258border.txt");
        sprintf(inputfilename_coherency,"shape101cell258coherencyv3.txt");
    }
    if(shape_number == 102){
        sprintf(inputfilename_orientation,"shape102cell267orientation.txt");
        sprintf(inputfilename_bulk,"shape102cell267bulk.txt");
        sprintf(inputfilename_border,"shape102cell267border.txt");
        sprintf(inputfilename_coherency,"shape102cell267coherencyv3.txt");
    }
    
    if(shape_number == 103){
        sprintf(inputfilename_orientation,"shape103cell303orientation.txt");
        sprintf(inputfilename_bulk,"shape103cell303bulk.txt");
        sprintf(inputfilename_border,"shape103cell303border.txt");
        sprintf(inputfilename_coherency,"shape103cell303coherencyv3.txt");
    }
    
    if(shape_number == 104){
        sprintf(inputfilename_orientation,"shape104cell312orientation.txt");
        sprintf(inputfilename_bulk,"shape104cell312bulk.txt");
        sprintf(inputfilename_border,"shape104cell312border.txt");
        sprintf(inputfilename_coherency,"shape104cell312coherencyv3.txt");
    }
    
    if(shape_number == 105){
        sprintf(inputfilename_orientation,"shape105cell792orientation.txt");
        sprintf(inputfilename_bulk,"shape105cell792bulk.txt");
        sprintf(inputfilename_border,"shape105cell792border.txt");
        sprintf(inputfilename_coherency,"shape105cell792coherencyv3.txt");
    }

    

        
    
       
    FILE *myFile1;
        myFile1 =fopen(inputfilename_orientation, "r");
    FILE *myFile2;
        myFile2 =fopen(inputfilename_bulk, "r");
    FILE *myFile3;
        myFile3 =fopen(inputfilename_border, "r");
    FILE *myFile4;
        myFile4 =fopen(inputfilename_coherency, "r");
        
        //read file to find the parameter values
    for(l=0; l<(num_of_gridpoints_exp); l++){
        fscanf(myFile1, "%ld\n", &grid_exp[l].phi_exp_degrees);
        fscanf(myFile2, "%ld\n", &grid_exp[l].bulk_exp);
        fscanf(myFile3, "%ld\n", &border_exp[l].border);
        fscanf(myFile4, "%lf\n", &grid_exp[l].coherency);
        //printf("\n%d\n", grid_exp[l].bulk_exp);
        
    }
    
    pix_per_point = num_of_gridpoints_on_line_exp/(num_of_gridpoints_on_line_in_box-1); // This is the width of each point after the rescaling, in pixels. ONLY WORKS IF THE SIZE OF THE MARGIN = 0!
    
    num_of_borderpoints=0;
    
    for(n=0; n<num_of_gridpoints_exp; n++){
        if(grid_exp[n].bulk_exp == 0){
            grid_exp[n].coherency = 0.0;
        }            
    }
    
    for(n=0; n<num_of_gridpoints_exp; n++){
        if(border_exp[n].border == 1){
            border_exp[num_of_borderpoints].x = ((n%num_of_gridpoints_on_line_exp)+0.5)*(1.0/num_of_gridpoints_on_line_exp);
            border_exp[num_of_borderpoints].y = ((n/num_of_gridpoints_on_line_exp)+0.5)*(1.0/num_of_gridpoints_on_line_exp);
            num_of_borderpoints++;
        }            
    }
    
    for(m=0; m<(num_of_gridpoints_exp); m++){
        grid_exp[m].phi_exp = (pi/180.0)*(180-grid_exp[m].phi_exp_degrees); //(180-grid_exp[m].phi_exp_degrees) mirrors the orientation, as the start point for our simulation is not the same as the origin in the raw data.
    }
        
    num_of_gridpoints_noborder=num_of_gridpoints_on_line_in_box*num_of_gridpoints_on_line_in_box;  
    
    if(num_of_gridpoints_exp<num_of_gridpoints_noborder){
        printf("\nThe amount of experimental data points is smaller than the amount in the simulation! Number of experimental points:\t%d\t Number of simulation points:\t%d\n", num_of_gridpoints_exp,num_of_gridpoints);
    }
    
    if(num_of_gridpoints_on_line_exp%((num_of_gridpoints_on_line_in_box-1)/2) != 0){
        printf("\nThe number of bounpoint per line must conform to the folllowing:num_of_gridpoints_on_line_exp must be straight divisible by ((num_of_gridpoints_on_line_in_box-1)/2) ! \n Number of experimental points:\t%d\t num_of_gridpoints_on_line_in_box:\t%d\n", num_of_gridpoints_on_line_exp,num_of_gridpoints_on_line_in_box);
    }
    
    for(i=0; i<num_of_gridpoints_noborder; i++){
        counts=0;
        sumq1=0;
        sumq2=0;
        if(i%num_of_gridpoints_on_line_in_box==0 || i%num_of_gridpoints_on_line_in_box==(num_of_gridpoints_on_line_in_box-1) || i<num_of_gridpoints_on_line_in_box || i>(num_of_gridpoints_noborder-num_of_gridpoints_on_line_in_box)){
            grid_exp[i].bulk = 0;
            grid_exp[i].q1 = 100.0;
            grid_exp[i].q2 = 100.0;
            grid_exp[i].s = 100.0;
            grid_exp[i].phi = 100.0;
            grid_exp[i].counts = 0;
            
        }
        else{
            for(j=0; j<pix_per_point; j++){
                for(k=0; k<pix_per_point; k++){
                    startpoint = num_of_gridpoints_on_line_exp*(0.5*pix_per_point) + (0.5*pix_per_point) + ((i/num_of_gridpoints_on_line_in_box)-1)*(num_of_gridpoints_on_line_exp*pix_per_point) + ((i%num_of_gridpoints_on_line_in_box)-1)*pix_per_point + j + k*num_of_gridpoints_on_line_exp; //0.5*pix_per_point is the edge we're cutting off on each side of the box, in pixels.
                    if(grid_exp[startpoint].bulk_exp==1){
                        counts++;
                        sumq1 = sumq1 + 0.5*grid_exp[startpoint].coherency*cos(2.0*grid_exp[startpoint].phi_exp); //S=1
                        sumq2 = sumq2 + 0.5*grid_exp[startpoint].coherency*sin(2.0*grid_exp[startpoint].phi_exp); //S=1
                    }
                }
           


            }
            
                
        }
        
        if(counts > 0.5*(pix_per_point*pix_per_point)){ //This is the cutoff
                grid_exp[i].bulk = 1;
                grid_exp[i].q1 = sumq1/counts;
                grid_exp[i].q2 = sumq2/counts;
                grid_exp[i].s = 2.0*sqrt(grid_exp[i].q1*grid_exp[i].q1+grid_exp[i].q2*grid_exp[i].q2);
                if(grid_exp[i].q1 > 0){
                    grid_exp[i].phi = 0.5*atan(grid_exp[i].q2/grid_exp[i].q1);
                }
                else if(grid_exp[i].q1 < 0){
                    if(grid_exp[i].q2 > 0){
                        grid_exp[i].phi = 0.5*(atan(grid_exp[i].q2/grid_exp[i].q1)+pi);
                    }
                    else if(grid_exp[i].q2 < 0){
                        grid_exp[i].phi = 0.5*(atan(grid_exp[i].q2/grid_exp[i].q1)-pi);
                    }
                    else{
                        grid_exp[i].phi = 0.5*pi;
                    }
                }
                else{
                    if(grid_exp[i].q2 > 0){
                        grid_exp[i].phi = 0.25*pi;
                    }
                    else{
                        grid_exp[i].phi = 0.75*pi;
                    }
                }    
                grid_exp[i].counts = counts;
                //printf("\nThis point is in the bulk. Point\t%ld\t With Q1 and Q2\t%.10f\t%.10f\t and counts\t%ld\n", i,grid_exp[i].q1,grid_exp[i].q2,grid_exp[i].counts);
            }
                
            else{
                    grid_exp[i].bulk = 0;
                    grid_exp[i].q1 = 100.0;
                    grid_exp[i].q2 = 100.0;
                    grid_exp[i].s = 100.0;
                    grid_exp[i].phi = 100.0;
                    grid_exp[i].counts = 0;
            }
        
    }
    /*
    printf("\nnum_of_gridpoints, num_of_gridpoints_exp, pix_per_point:\n");    
    printf("%ld\n", num_of_gridpoints);
    printf("%ld\n", num_of_gridpoints_exp);
    printf("%ld\n", pix_per_point);*/
       
}

void comparison_exp()
{
    long i,j;
    double comx=0.0,comy=0.0;
    double discomx, discomy, discom, discomtotal=0.0, discomtotal2=0.0;
    

    
    for(j=0; j<num_of_gridpoints; j++){
        if(grid_exp[j].bulk == 1){
            bulkpoints_exp++;

        }
        if(grid[j].bulk == 1){
            bulkpoints_sim++;
        }
        if(grid_exp[j].bulk == 1 && grid[j].bulk == 1){
            grid[j].sharedbulk = 1;
            bulkpoints_shared++;
            comx=comx+(j%num_of_gridpoints_on_line);
            comy=comy+(j/num_of_gridpoints_on_line);
        }
        else{
            grid[j].sharedbulk = 0;
        }
    }
    
    comx=comx/(1.0*bulkpoints_shared);
    comy=comy/(1.0*bulkpoints_shared);
    //printf("\nCOM: %.10f\t%.10f\n", comx,comy);
    
    for(i=0; i<num_of_gridpoints; i++){// comparisons
        if(grid[i].sharedbulk == 1){
            deltasquaredphi = deltasquaredphi + ((grid_exp[i].phi-grid[i].phi)*(grid_exp[i].phi-grid[i].phi));
            deltasquaredq1 = deltasquaredq1 + ((grid_exp[i].q1-grid[i].q1)*(grid_exp[i].q1-grid[i].q1));
            deltasquaredq2 = deltasquaredq2 + ((grid_exp[i].q2-grid[i].q2)*(grid_exp[i].q2-grid[i].q2));
            deltasquaredangle = deltasquaredangle + (sin(grid_exp[i].phi-grid[i].phi)*sin(grid_exp[i].phi-grid[i].phi));
            discomx=fabs((i%num_of_gridpoints_on_line)*1.0-comx);
            discomy=fabs((i/num_of_gridpoints_on_line)*1.0-comy);
            discom=sqrt(discomx*discomx+discomy*discomy);
            deltasquaredanglecom = deltasquaredanglecom + discom*(sin(grid_exp[i].phi-grid[i].phi)*sin(grid_exp[i].phi-grid[i].phi));
            deltasquaredanglecom2 = deltasquaredanglecom2 + discom*discom*(sin(grid_exp[i].phi-grid[i].phi)*sin(grid_exp[i].phi-grid[i].phi));
            discomtotal=discomtotal+discom;
            discomtotal2=discomtotal2+(discom*discom);

                //printf("\nCOM: \t%ld\t %.10f\t%.10f\t%.10f\n", i,discomx,discomy,discom);
        }
    }
    
    deltasquaredphi = deltasquaredphi/(1.0*bulkpoints_shared);
    deltasquaredq1 = deltasquaredq1/(1.0*bulkpoints_shared);
    deltasquaredq2 = deltasquaredq2/(1.0*bulkpoints_shared);
    deltasquaredq = deltasquaredq1+deltasquaredq2;
    deltasquaredangle = deltasquaredangle/(1.0*bulkpoints_shared);
    deltasquaredanglecom = deltasquaredanglecom/discomtotal;
    deltasquaredanglecom2 = deltasquaredanglecom2/discomtotal2;
    
    //Jackknife method:
    
    int bulkpointcounter;

    for(i=0; i<bulkpoints_shared; i++){
        bulkpointcounter=0;
        deltasquaredjkq1=0.0;
        deltasquaredjkq2=0.0;
        deltasquaredjkq=0.0;
        for(j=0; j<num_of_gridpoints; j++){
            if(grid[j].sharedbulk == 1){
                if(i!=bulkpointcounter){
                    deltasquaredjkq1 = deltasquaredjkq1 + ((grid_exp[j].q1-grid[j].q1)*(grid_exp[j].q1-grid[j].q1));
                    deltasquaredjkq2 = deltasquaredjkq2 + ((grid_exp[j].q2-grid[j].q2)*(grid_exp[j].q2-grid[j].q2));
                }
                bulkpointcounter++;
            }
        }
        deltasquaredjkq1 = deltasquaredjkq1/(1.0*(bulkpoints_shared-1));
        deltasquaredjkq2 = deltasquaredjkq2/(1.0*(bulkpoints_shared-1));
        deltasquaredjkq = deltasquaredjkq1+deltasquaredjkq2;
        //printf("\nJackknife deltasquared index (%ld): %.10f\n", i,deltasquaredjkq);
        errordeltasquared = errordeltasquared + (deltasquaredjkq-deltasquaredq)*(deltasquaredjkq-deltasquaredq);
    }
    
    errordeltasquared = sqrt(((bulkpoints_shared-1)/(1.0*bulkpoints_shared))*errordeltasquared);
    
    
    
    

    printf("\nNumber of bulkpoints in simulation: %ld\n", bulkpoints_sim);
    printf("\nNumber of bulkpoints in experiment: %ld\n", bulkpoints_exp);
    printf("\nNumber of bulkpoint shared: %ld\n", bulkpoints_shared);
    printf("\ndeltasquared(Q): %.10f\n", deltasquaredq);
    printf("\ndeltasquaredjk(Q): %.10f\n", deltasquaredjkq);
    printf("\nerror deltasquared(Q) using jackknife: %.10f\n", errordeltasquared);
    
}

void init_boundary()
{
    
	double x, y, dx, dy;
    long a, b, c;
    long i,j,k;
    long l, m, n, o, p, q;
    double a1;
    double rotationpi, size;
    double center_x, center_y;
    double first_adhesion[1][2];
    double first_adhesion_rotated[1][2];
    double adhesions_rotated[num_of_adhesions][2];
    double adhesions_final[num_of_adhesions][2];
    double adhesions_relative[num_of_adhesions][2]; 
    double rotation_increment, first_rotation;  
    double adhesions_exp[num_of_adhesions][2];
    double adhesions_exp_temp[num_of_adhesions][2];
    double adhesions_scaled[num_of_adhesions][2];
    double max_x, max_y, min_x, min_y, bordersize, size_x, size_y, size_max;
    long num_of_pixels, cellnr;
    double adhesions_mirrored[num_of_adhesions][2];
    double mirroring, min_total;
    double timestepbulkhelp,time_stephelp;
    double corner;

    
    size = 1.0; //add
    center_x = 0.5; // add
    center_y = 0.5; // add
    rotationpi = rotation*pi;


	// Construct a square boundary
	
        if(shape_number == 1){
            num_of_edgepoint = num_of_bounpoint/4; 
	
            dx = 0.99998*anisotropy/(num_of_edgepoint-1.0);                                //I use numbers very close to 1 and 0 so that the initial edge will not exactly be on grid points. Determinebulk works better in that case.
            dy = 0.99998/(num_of_edgepoint-1.0);
	
            for (i=0; i<num_of_edgepoint; i++){
            
		vertex[i].x = (0.00001*anisotropy)+i*dx;
		vertex[i].y = 0.99999;
		vertex[i].boundary = 1;
                vertex[i].lam = lambda;
		
		vertex[num_of_edgepoint+i].x = 0.99999*anisotropy;
		vertex[num_of_edgepoint+i].y = 0.99999-i*dy;
		vertex[num_of_edgepoint+i].boundary = 1;
                vertex[num_of_edgepoint+i].lam = lambda;
		
		vertex[2*num_of_edgepoint+i].x = 0.99999*anisotropy-i*dx;
		vertex[2*num_of_edgepoint+i].y = 0.00001;
		vertex[2*num_of_edgepoint+i].boundary = 1;
                vertex[2*num_of_edgepoint+i].lam = lambda;
                
		vertex[3*num_of_edgepoint+i].x = 0.00001*anisotropy;
		vertex[3*num_of_edgepoint+i].y = 0.00001+i*dy;
		vertex[3*num_of_edgepoint+i].boundary = 1;
                vertex[3*num_of_edgepoint+i].lam = lambda;
                
                if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                    vertex[i].adhesion =1;
                    vertex[num_of_edgepoint+i].adhesion =1;
                    vertex[2*num_of_edgepoint+i].adhesion =1;
                    vertex[3*num_of_edgepoint+i].adhesion =1;
                }
                else{
                    vertex[i].adhesion = 0;
                    vertex[num_of_edgepoint+i].adhesion =0;
                    vertex[2*num_of_edgepoint+i].adhesion =0;
                    vertex[3*num_of_edgepoint+i].adhesion =0;
                }
            }
            
        }
        
        if(shape_number ==2){
            printf("%i\n",num_of_edgepoint);
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions[k][0]-adhesions[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions[k][1]-adhesions[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
        
        if(shape_number ==3){ //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.
            num_of_adhesions = 13;
            
            adhesions[0][0] = 0.74;
            adhesions[0][1] = 0.068;
            adhesions[1][0] = 0.68;
            adhesions[1][1] = 0.125;
            adhesions[2][0] = 0.505;
            adhesions[2][1] = 0.286;
            adhesions[3][0] = 0.192;
            adhesions[3][1] =0.309;
            adhesions[4][0] = 0.255;
            adhesions[4][1] =0.551;
            adhesions[5][0] = 0.037;
            adhesions[5][1] =0.878;
            adhesions[6][0] = 0.045;
            adhesions[6][1] =0.963;
            adhesions[7][0] = 0.133;
            adhesions[7][1] = 0.995;
            adhesions[8][0] = 0.386;
            adhesions[8][1] = 0.742;
            adhesions[9][0] = 0.607;
            adhesions[9][1] =0.719;
            adhesions[10][0] = 0.919;
            adhesions[10][1] =0.742;
            adhesions[11][0] =0.979; 
            adhesions[11][1] =0.668;
            adhesions[12][0] = 0.85;
            adhesions[12][1] =0.162;
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions[k][0]-adhesions[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions[k][1]-adhesions[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
        
        if(shape_number ==4){ //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.
            num_of_adhesions = 3;
            size = 1.0; //add
            
            adhesions[0][0] = 0.95;
            adhesions[0][1] = 0.10;
            adhesions[1][0] = 0.05;
            adhesions[1][1] = 0.10;
            adhesions[2][0] = 0.50;
            adhesions[2][1] = 0.88;

            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions[k][0]-adhesions[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions[k][1]-adhesions[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }  

    
        if(shape_number == 5){ //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            rotation_increment = 2.0*pi/num_of_adhesions;
            printf("\n%i\t%.10f", num_of_adhesions, rotation_increment);
            first_rotation = 0.5*pi;

            first_adhesion[0][0] = 0; //The initial adhesion point that we will rotate later. Coordinates (0,size).
            first_adhesion[0][1] = size;
            
            first_adhesion_rotated[0][0] = first_adhesion[0][0]*cos(first_rotation) - first_adhesion[0][1]*sin(first_rotation); //Rotate the initial adhesion point over an angle of first_rotation.
            first_adhesion_rotated[0][1] = first_adhesion[0][0]*sin(first_rotation) + first_adhesion[0][1]*cos(first_rotation);
            
            a1 = 0.0;
            
            for(a=0;a<num_of_adhesions;a++){ //Generate all the remaining adhesion points around (0,0).
                adhesions_rotated[a][0] = first_adhesion_rotated[0][0]*cos(a1*rotation_increment) + first_adhesion_rotated[0][1]*sin(a1*rotation_increment); 
                adhesions_rotated[a][1] = -first_adhesion_rotated[0][0]*sin(a1*rotation_increment) + first_adhesion_rotated[0][1]*cos(a1*rotation_increment);
                
                printf("\nrotated");
                printf("\n%.10f\t",a1*rotation_increment);
                printf("\n%.10f\t%.10f", adhesions_rotated[a][0], adhesions_rotated[a][1]);
                a1++;
            }            
            
            for(b=0;b<num_of_adhesions;b++){ //Center adheion point around (center_x,center_y).
                adhesions_final[b][0] = adhesions_rotated[b][0] + center_x;
                adhesions_final[b][1] = adhesions_rotated[b][1] + center_y;
            }
            

            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions_final[k][0]-adhesions_final[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions_final[k][1]-adhesions_final[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions_final[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions_final[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
        
        if(shape_number == 6){ //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            num_of_adhesions = 14;
            
            adhesions[0][0] = 0.74;
            adhesions[0][1] = 0.068;
            adhesions[1][0] = 0.68;
            adhesions[1][1] = 0.125;
            adhesions[2][0] = 0.505;
            adhesions[2][1] = 0.286;
            adhesions[3][0] = 0.192;
            adhesions[3][1] =0.309;
            adhesions[4][0] = 0.255;
            adhesions[4][1] =0.551;
            adhesions[5][0] = 0.037;
            adhesions[5][1] =0.878;
            adhesions[6][0] = 0.045;
            adhesions[6][1] =0.963;
            adhesions[7][0] = 0.133;
            adhesions[7][1] = 0.995;
            adhesions[8][0] = 0.386;
            adhesions[8][1] = 0.742;
            adhesions[9][0] = 0.607;
            adhesions[9][1] =0.719;
            adhesions[10][0] = 0.919;
            adhesions[10][1] =0.742;
            adhesions[11][0] =0.979; 
            adhesions[11][1] =0.668;
            adhesions[12][0] = 0.815;
            adhesions[12][1] = 0.306;
            adhesions[13][0] = 0.85;
            adhesions[13][1] =0.162;

            
            for(a=0;a<num_of_adhesions;a++){
                adhesions_relative[a][0] = (adhesions[a][0] - center_x)*size; //Scale according to size and move to origin for rotation.
                adhesions_relative[a][1] = (adhesions[a][1] - center_y)*size;
            }
            
            for(b=0;b<num_of_adhesions;b++){
                adhesions_rotated[b][0] = adhesions_relative[b][0]*cos(rotationpi) + adhesions_relative[b][1]*sin(rotationpi);
                adhesions_rotated[b][1] = -adhesions_relative[b][0]*sin(rotationpi) + adhesions_relative[b][1]*cos(rotationpi);
                printf("\n%.10f\t%.10f", adhesions_rotated[b][0], adhesions_rotated[b][1]);
            }

            for(c=0;c<num_of_adhesions;c++){
                adhesions_final[c][0] = adhesions_rotated[c][0]+center_x;
                adhesions_final[c][1] = adhesions_rotated[c][1]+center_y;
            }

 
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions_final[k][0]-adhesions_final[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions_final[k][1]-adhesions_final[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions_final[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions_final[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
    
        if(shape_number == 7){ // This is the test function for importing the experimental data.
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            size = 1.0; //add
            center_x = 0.5; // add
            center_y = 0.5; // add
            num_of_adhesions = 3;        
            
            
            adhesions_exp[0][0] = 50; 
            adhesions_exp[0][1] = 200;
            adhesions_exp[1][0] = 200;
            adhesions_exp[1][1] = 250;
            adhesions_exp[2][0] = 250;
            adhesions_exp[2][1] = 200;
            
            for(k=0;k<num_of_adhesions;k++){
                adhesions[k][0] = adhesions_exp[k][0]/num_of_pixels;
                adhesions[k][1] = adhesions_exp[k][1]/num_of_pixels;
                printf("\nadhesions new\t%.10f\t%.10f", adhesions[k][0], adhesions[k][1]);
            } // Rescaling from num_of_pixels pixels to a box of size 1.
            
            max_x = 0;
            max_y = 0;
            min_x = 1;
            min_y = 1;
            
            for(l=0;l<num_of_adhesions;l++){
                if(adhesions[l][0]<min_x){
                    min_x=adhesions[l][0];                    
                }
                if(adhesions[l][0]>max_x){
                    max_x=adhesions[l][0]; 
                }
                if(adhesions[l][1]<min_y){
                    min_y=adhesions[l][1]; 
                }
                if(adhesions[l][1]>max_y){
                    max_y=adhesions[l][1]; 
                }
                printf("\nmax_x,max_y,min_x,min_y\t%.10f\t%.10f\t%.10f\t%.10f", max_x, max_y, min_x, min_y);
            } //Finding furthest values for x and y
            
            size_x = max_x-min_x;
            size_y = max_y-min_y;
            bordersize = 0.1;
            
            

            
            if(size_x>size_y){
                size_max=size_x;
                printf("\nscaled to X");
                for(m=0;m<num_of_adhesions;m++){            
                    adhesions_scaled[m][0]=bordersize+(((adhesions[m][0]-min_x)/size_max)*(1.0-2.0*bordersize));
                    adhesions_scaled[m][1]=bordersize+((((adhesions[m][1]-min_y)/size_max)+0.5-(0.5*size_y))*(1.0-2.0*bordersize));// The factor of 0.5-(0.5*size_y) is to move the center of cell to the center of the box.
                    printf("\nadhesions_scaled\t%.10f\t%.10f", adhesions_scaled[m][0], adhesions_scaled[m][1]);
                }
            }
            else{
                size_max=size_y;  
                printf("\nscaled to Y");
            } //Finding largest total size
            
            printf("\nsize_x, size_y, size_max\t%.10f\t%.10f\t%.10f\t%.10f", size_x, size_y, size_max);   
            /*
            for(m=0;m<num_of_adhesions;m++){            
                adhesions_scaled[m][0]=bordersize+(((adhesions[m][0]-min_x)/size_max)*(1.0-2.0*bordersize));
                adhesions_scaled[m][1]=bordersize+(((adhesions[m][1]-min_y)/size_max)*(1.0-2.0*bordersize));
                printf("\nadhesions_scaled\t%.10f\t%.10f", adhesions_scaled[m][0], adhesions_scaled[m][1]);
            }*/

                
                
                
                
            for(a=0;a<num_of_adhesions;a++){
                adhesions_relative[a][0] = (adhesions_scaled[a][0] - center_x)*size; //Scale according to size and move to origin for rotation.
                adhesions_relative[a][1] = (adhesions_scaled[a][1] - center_y)*size;
            }
            // adhesions_exp, max_x, min_x, max_y, max_y, size_x, size_y, size_max, border, min_xy, max_xy

            
            for(b=0;b<num_of_adhesions;b++){
                adhesions_rotated[b][0] = adhesions_relative[b][0]*cos(rotationpi) + adhesions_relative[b][1]*sin(rotationpi);
                adhesions_rotated[b][1] = -adhesions_relative[b][0]*sin(rotationpi) + adhesions_relative[b][1]*cos(rotationpi);
                printf("\n%.10f\t%.10f", adhesions_rotated[b][0], adhesions_rotated[b][1]);
            }

            for(c=0;c<num_of_adhesions;c++){
                adhesions_final[c][0] = adhesions_rotated[c][0]+center_x;
                adhesions_final[c][1] = adhesions_rotated[c][1]+center_y;
            }

 
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions_final[k][0]-adhesions_final[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions_final[k][1]-adhesions_final[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions_final[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions_final[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
        
        if(shape_number == 8){ // I'm currently using this function for the boxes.

            num_of_adhesions = 4;
            corner=0.5*(1.0/(num_of_gridpoints_on_line_in_box-1));
            
            adhesions[0][0] = corner;
            adhesions[0][1] = corner;
            adhesions[1][0] = corner;
            adhesions[1][1] = 1-corner;
            adhesions[2][0] = 1-corner;
            adhesions[2][1] = 1-corner;
            adhesions[3][0] = 1-corner;
            adhesions[3][1] = corner;
            
            
            for(a=0;a<num_of_adhesions;a++){
                adhesions_relative[a][0] = (adhesions[a][0] - center_x)*size; //Scale according to size and move to origin for rotation.
                adhesions_relative[a][1] = (adhesions[a][1] - center_y)*size;
            }
            
            for(b=0;b<num_of_adhesions;b++){
                adhesions_rotated[b][0] = adhesions_relative[b][0]*cos(rotationpi) + adhesions_relative[b][1]*sin(rotationpi);
                adhesions_rotated[b][1] = -adhesions_relative[b][0]*sin(rotationpi) + adhesions_relative[b][1]*cos(rotationpi);
                printf("\n%.10f\t%.10f", adhesions_rotated[b][0], adhesions_rotated[b][1]);
            }

            for(c=0;c<num_of_adhesions;c++){
                adhesions_final[c][0] = adhesions_rotated[c][0]+center_x;
                adhesions_final[c][1] = adhesions_rotated[c][1]+center_y;
            }

 
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions_final[k][0]-adhesions_final[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions_final[k][1]-adhesions_final[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions_final[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions_final[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
        }
        
        
        if(shape_number == 9){
            
            num_of_adhesions = 4;
            corner=0.5*(1.0/(num_of_gridpoints_on_line_in_box-1));
            num_of_edgepoint = num_of_bounpoint/4; 
	
            dx = (1.0-2.0*corner)*anisotropy/(num_of_edgepoint-1.0);                                //I use numbers very close to 1 and 0 so that the initial edge will not exactly be on grid points. Determinebulk works better in that case.
            dy = (1.0-2.0*corner)/(num_of_edgepoint-1.0);
	
            for (i=0; i<num_of_edgepoint; i++){
            
                vertex[i].x = corner+i*dx;
                vertex[i].y = (1-corner);
                vertex[i].boundary = 1;
                vertex[i].lam = lambda;
		
                vertex[num_of_edgepoint+i].x = (1-corner)+((anisotropy-1)*(1.0-2.0*corner));
                vertex[num_of_edgepoint+i].y = (1-corner)-i*dy;
                vertex[num_of_edgepoint+i].boundary = 1;
                vertex[num_of_edgepoint+i].lam = lambda;
		
                vertex[2*num_of_edgepoint+i].x = (1-corner)+((anisotropy-1)*(1.0-2.0*corner))-i*dx;
                vertex[2*num_of_edgepoint+i].y = corner;
                vertex[2*num_of_edgepoint+i].boundary = 1;
                vertex[2*num_of_edgepoint+i].lam = lambda;
                
                vertex[3*num_of_edgepoint+i].x = corner;
                vertex[3*num_of_edgepoint+i].y = corner+i*dy;
                vertex[3*num_of_edgepoint+i].boundary = 1;
                vertex[3*num_of_edgepoint+i].lam = lambda;
                
                if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                    vertex[i].adhesion =1;
                    vertex[num_of_edgepoint+i].adhesion =1;
                    vertex[2*num_of_edgepoint+i].adhesion =1;
                    vertex[3*num_of_edgepoint+i].adhesion =1;
                }
                else{
                    vertex[i].adhesion = 0;
                    vertex[num_of_edgepoint+i].adhesion =0;
                    vertex[2*num_of_edgepoint+i].adhesion =0;
                    vertex[3*num_of_edgepoint+i].adhesion =0;
                }
            }
            
        }
        
        if(shape_number == 10){
            
            num_of_adhesions = 4;
            corner=0.5*(1.0/(num_of_gridpoints_on_line_in_box-1));
            num_of_edgepoint = num_of_bounpoint/4; 
	
            dx = (1.0-2.0*corner)*anisotropy/(num_of_edgepoint-1.0);                                //I use numbers very close to 1 and 0 so that the initial edge will not exactly be on grid points. Determinebulk works better in that case.
            dy = (1.0-2.0*corner)/anisotropy/(num_of_edgepoint-1.0);
	
            for (i=0; i<num_of_edgepoint; i++){
            
                vertex[i].x = corner+i*dx;
                vertex[i].y = (1.0-2.0*corner)/anisotropy+corner;
                vertex[i].boundary = 1;
                vertex[i].lam = lambda;
		
                vertex[num_of_edgepoint+i].x = (1-corner)+((anisotropy-1)*(1.0-2.0*corner));
                vertex[num_of_edgepoint+i].y = (1.0-2.0*corner)/anisotropy+corner-i*dy;
                vertex[num_of_edgepoint+i].boundary = 1;
                vertex[num_of_edgepoint+i].lam = lambda;
		
                vertex[2*num_of_edgepoint+i].x = (1-corner)+((anisotropy-1)*(1.0-2.0*corner))-i*dx;
                vertex[2*num_of_edgepoint+i].y = corner;
                vertex[2*num_of_edgepoint+i].boundary = 1;
                vertex[2*num_of_edgepoint+i].lam = lambda;
                
                vertex[3*num_of_edgepoint+i].x = corner;
                vertex[3*num_of_edgepoint+i].y = corner+i*dy;
                vertex[3*num_of_edgepoint+i].boundary = 1;
                vertex[3*num_of_edgepoint+i].lam = lambda;
                
                if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                    vertex[i].adhesion =1;
                    vertex[num_of_edgepoint+i].adhesion =1;
                    vertex[2*num_of_edgepoint+i].adhesion =1;
                    vertex[3*num_of_edgepoint+i].adhesion =1;
                }
                else{
                    vertex[i].adhesion = 0;
                    vertex[num_of_edgepoint+i].adhesion =0;
                    vertex[2*num_of_edgepoint+i].adhesion =0;
                    vertex[3*num_of_edgepoint+i].adhesion =0;
                }
            }
            
        }
        
        if(shape_number == 11){ // I'm currently using this function for the boxes.

            num_of_adhesions = 3;
            corner=0.5*(1.0/(num_of_gridpoints_on_line_in_box-1));
            
            adhesions[0][0] = corner;
            adhesions[0][1] = corner;
            adhesions[1][0] = 0.5;
            adhesions[1][1] = 0.5*sqrt(3.0)*(1.0-2.0*corner)+corner;
            adhesions[2][0] = 1-corner;
            adhesions[2][1] = corner;

	
            for (i=0; i<num_of_edgepoint; i++){
            
                vertex[i].x = corner+i*(0.5-corner)/(num_of_edgepoint-1.0);
                vertex[i].y = corner+i*0.5*sqrt(3)*(1.0-2.0*corner)/(num_of_edgepoint-1.0);
                vertex[i].boundary = 1;
                vertex[i].lam = lambda;
		
                vertex[num_of_edgepoint+i].x = 0.5+i*(0.5-corner)/(num_of_edgepoint-1.0);
                vertex[num_of_edgepoint+i].y = 0.5*sqrt(3.0)*(1.0-2.0*corner)+corner-i*0.5*sqrt(3)*(1.0-2.0*corner)/(num_of_edgepoint-1.0);
                vertex[num_of_edgepoint+i].boundary = 1;
                vertex[num_of_edgepoint+i].lam = lambda;
		
                vertex[2*num_of_edgepoint+i].x = (1-corner)-i*(1.0-2.0*corner)/(num_of_edgepoint-1.0);
                vertex[2*num_of_edgepoint+i].y = corner;
                vertex[2*num_of_edgepoint+i].boundary = 1;
                vertex[2*num_of_edgepoint+i].lam = lambda;

                
                if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                    vertex[i].adhesion =1;
                    vertex[num_of_edgepoint+i].adhesion =1;
                    vertex[2*num_of_edgepoint+i].adhesion =1;
                }
                else{
                    vertex[i].adhesion = 0;
                    vertex[num_of_edgepoint+i].adhesion =0;
                    vertex[2*num_of_edgepoint+i].adhesion =0;
                }
            }
            
        }
    
        
        
                    /*
            // Test cell chosen simply because it was number 1. Was not used in the paper.
            num_of_adhesions = 19;
            adhesions_exp[0][0] = 79; 
            adhesions_exp[0][1] = 235;
            adhesions_exp[1][0] = 107;
            adhesions_exp[1][1] = 167;
            adhesions_exp[2][0] = 98;
            adhesions_exp[2][1] = 214;
            adhesions_exp[3][0] = 106; 
            adhesions_exp[3][1] = 240;
            adhesions_exp[4][0] = 152;
            adhesions_exp[4][1] = 150;
            adhesions_exp[5][0] = 135;
            adhesions_exp[5][1] = 244;
            adhesions_exp[6][0] = 170; 
            adhesions_exp[6][1] = 129;
            adhesions_exp[7][0] = 163;
            adhesions_exp[7][1] = 249;
            adhesions_exp[8][0] = 204;
            adhesions_exp[8][1] = 87;
            adhesions_exp[9][0] = 226; 
            adhesions_exp[9][1] = 285;
            adhesions_exp[10][0] = 261;
            adhesions_exp[10][1] = 245;
            adhesions_exp[11][0] = 337;
            adhesions_exp[11][1] = 308;
            adhesions_exp[12][0] = 328; 
            adhesions_exp[12][1] = 354;
            adhesions_exp[13][0] = 363;
            adhesions_exp[13][1] = 384;
            adhesions_exp[14][0] = 400;
            adhesions_exp[14][1] = 416;
            adhesions_exp[15][0] = 420;
            adhesions_exp[15][1] = 399;
            adhesions_exp[16][0] = 429;
            adhesions_exp[16][1] = 422;
            adhesions_exp[17][0] = 455;
            adhesions_exp[17][1] = 428;
            adhesions_exp[18][0] = 484;
            adhesions_exp[18][1] = 434;
            */
         if(shape_number >= 100){ // Applies the same rules to any experimental cell tupe (nr >= 100)
            
            num_of_pixels = 512;  
             
            if(shape_number == 100){ // This is cell 253 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.
            
            /*Cell 253 - also in the paper
            γ:      0.400612±0.006628
            λ_/σ:   14.681733±0.089017 μm
            pxsize: 0.138000 μm
            adhesions_exp[0][0] = 132; 
            adhesions_exp[0][1] = 306;
            adhesions_exp[1][0] = 142;
            adhesions_exp[1][1] = 232;
            adhesions_exp[2][0] = 147;
            adhesions_exp[2][1] = 285;
            adhesions_exp[3][0] = 156; 
            adhesions_exp[3][1] = 208;
            adhesions_exp[4][0] = 170;
            adhesions_exp[4][1] = 185;
            adhesions_exp[5][0] = 177;
            adhesions_exp[5][1] = 424;
            adhesions_exp[6][0] = 183; 
            adhesions_exp[6][1] = 160;
            adhesions_exp[7][0] = 197;
            adhesions_exp[7][1] = 184;
            adhesions_exp[8][0] = 202;
            adhesions_exp[8][1] = 420;
            adhesions_exp[9][0] = 226; 
            adhesions_exp[9][1] = 376;
            adhesions_exp[10][0] = 240;
            adhesions_exp[10][1] = 209;
            adhesions_exp[11][0] = 240;
            adhesions_exp[11][1] = 353;
            adhesions_exp[12][0] = 285; 
            adhesions_exp[12][1] = 328;
            adhesions_exp[13][0] = 297;
            adhesions_exp[13][1] = 207;
            adhesions_exp[14][0] = 296;
            adhesions_exp[14][1] = 303;
            adhesions_exp[15][0] = 312;
            adhesions_exp[15][1] = 278;
            adhesions_exp[16][0] = 335;
            adhesions_exp[16][1] = 224;
            adhesions_exp[17][0] = 369;
            adhesions_exp[17][1] = 186;
            adhesions_exp[18][0] = 375;
            adhesions_exp[18][1] = 193;
            adhesions_exp[19][0] = 405;
            adhesions_exp[19][1] = 158;
            adhesions_exp[20][0] = 418;
            adhesions_exp[20][1] = 177;
            */
            
            adhesions_exp[0][0] = 132; 
            adhesions_exp[0][1] = 306;
            adhesions_exp[1][0] = 147;
            adhesions_exp[1][1] = 285;
            adhesions_exp[2][0] = 142;
            adhesions_exp[2][1] = 232;
            adhesions_exp[3][0] = 156; 
            adhesions_exp[3][1] = 208;
            adhesions_exp[4][0] = 170;
            adhesions_exp[4][1] = 185;
            adhesions_exp[5][0] = 183; 
            adhesions_exp[5][1] = 160;
            adhesions_exp[6][0] = 197;
            adhesions_exp[6][1] = 184;
            adhesions_exp[7][0] = 240;
            adhesions_exp[7][1] = 209;          
            adhesions_exp[8][0] = 297;
            adhesions_exp[8][1] = 207;         
            adhesions_exp[9][0] = 369;
            adhesions_exp[9][1] = 186;           
            adhesions_exp[10][0] = 405;
            adhesions_exp[10][1] = 158;
            adhesions_exp[11][0] = 418;
            adhesions_exp[11][1] = 177;
            adhesions_exp[12][0] = 375;
            adhesions_exp[12][1] = 193; 
            adhesions_exp[13][0] = 335;
            adhesions_exp[13][1] = 224;
            adhesions_exp[14][0] = 312;
            adhesions_exp[14][1] = 278;
            adhesions_exp[15][0] = 296;
            adhesions_exp[15][1] = 303;            
            adhesions_exp[16][0] = 285; 
            adhesions_exp[16][1] = 328;
            adhesions_exp[17][0] = 240;
            adhesions_exp[17][1] = 353;
            adhesions_exp[18][0] = 226; 
            adhesions_exp[18][1] = 376;
            adhesions_exp[19][0] = 202;
            adhesions_exp[19][1] = 420;
            adhesions_exp[20][0] = 177;
            adhesions_exp[20][1] = 424;
            }
            
            
            
        else if(shape_number == 101){ // This is cell 258 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.
            
            /*Cell 258 - also in the paper
            γ:      0.947604±0.012780
            λ_/σ:   10.828967±0.052359 μm
            pxsize: 0.138000 μm
            adhesions_exp[0][0] = 30; 
            adhesions_exp[0][1] = 186;
            adhesions_exp[1][0] = 54;
            adhesions_exp[1][1] = 189;
            adhesions_exp[2][0] = 73;
            adhesions_exp[2][1] = 208;
            adhesions_exp[3][0] = 73; 
            adhesions_exp[3][1] = 256;
            adhesions_exp[4][0] = 142;
            adhesions_exp[4][1] = 227;
            adhesions_exp[5][0] = 156;
            adhesions_exp[5][1] = 208;
            adhesions_exp[6][0] = 210; 
            adhesions_exp[6][1] = 106;
            adhesions_exp[7][0] = 214;
            adhesions_exp[7][1] = 250;
            adhesions_exp[8][0] = 230;
            adhesions_exp[8][1] = 277;
            adhesions_exp[9][0] = 235; 
            adhesions_exp[9][1] = 59;
            adhesions_exp[10][0] = 242;
            adhesions_exp[10][1] = 300;
            adhesions_exp[11][0] = 249;
            adhesions_exp[11][1] = 82;
            adhesions_exp[12][0] = 293; 
            adhesions_exp[12][1] = 157;
            adhesions_exp[13][0] = 299;
            adhesions_exp[13][1] = 300;
            adhesions_exp[14][0] = 321;
            adhesions_exp[14][1] = 155;
            adhesions_exp[15][0] = 342;
            adhesions_exp[15][1] = 324;
            adhesions_exp[16][0] = 351;
            adhesions_exp[16][1] = 204;
            adhesions_exp[17][0] = 350;
            adhesions_exp[17][1] = 251;
            adhesions_exp[18][0] = 365;
            adhesions_exp[18][1] = 225;
            adhesions_exp[19][0] = 363;
            adhesions_exp[19][1] = 277;
            adhesions_exp[20][0] = 370;
            adhesions_exp[20][1] = 373;
            adhesions_exp[21][0] = 379;
            adhesions_exp[21][1] = 302;
            adhesions_exp[22][0] = 396;
            adhesions_exp[22][1] = 370;
            adhesions_exp[23][0] = 410;
            adhesions_exp[23][1] = 346;
            adhesions_exp[24][0] = 452;
            adhesions_exp[24][1] = 322;
            */
            
            adhesions_exp[0][0] = 30; 
            adhesions_exp[0][1] = 186;
            adhesions_exp[1][0] = 54;
            adhesions_exp[1][1] = 189;
            adhesions_exp[2][0] = 156;
            adhesions_exp[2][1] = 208;
            adhesions_exp[3][0] = 210; 
            adhesions_exp[3][1] = 106;
            adhesions_exp[4][0] = 235; 
            adhesions_exp[4][1] = 59;
            adhesions_exp[5][0] = 249;
            adhesions_exp[5][1] = 82;
            adhesions_exp[6][0] = 293; 
            adhesions_exp[6][1] = 157;            
            adhesions_exp[7][0] = 321;
            adhesions_exp[7][1] = 155;            
            adhesions_exp[8][0] = 351;
            adhesions_exp[8][1] = 204;            
            adhesions_exp[9][0] = 365;
            adhesions_exp[9][1] = 225;            
            adhesions_exp[10][0] = 350;
            adhesions_exp[10][1] = 251;            
            adhesions_exp[11][0] = 363;
            adhesions_exp[11][1] = 277;
            adhesions_exp[12][0] = 379;
            adhesions_exp[12][1] = 302;            
            adhesions_exp[13][0] = 452;
            adhesions_exp[13][1] = 322;                        
            adhesions_exp[14][0] = 410;
            adhesions_exp[14][1] = 346;
            adhesions_exp[15][0] = 396;
            adhesions_exp[15][1] = 370;            
            adhesions_exp[16][0] = 370;
            adhesions_exp[16][1] = 373;           
            adhesions_exp[17][0] = 342;
            adhesions_exp[17][1] = 324;
            adhesions_exp[18][0] = 299;
            adhesions_exp[18][1] = 300;
            adhesions_exp[19][0] = 242;
            adhesions_exp[19][1] = 300;
            adhesions_exp[20][0] = 230;
            adhesions_exp[20][1] = 277;
            adhesions_exp[21][0] = 214;
            adhesions_exp[21][1] = 250;
            adhesions_exp[22][0] = 142;
            adhesions_exp[22][1] = 227;
            adhesions_exp[23][0] = 73; 
            adhesions_exp[23][1] = 256;
            adhesions_exp[24][0] = 73;
            adhesions_exp[24][1] = 208;
            }
            
        else if(shape_number == 102){ // This is cell 267 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            //Cell 267 - also in the paper
            
            
            /*
            adhesions_exp[0][0] = 53; 
            adhesions_exp[0][1] = 161;
            adhesions_exp[1][0] = 79;
            adhesions_exp[1][1] = 115;
            adhesions_exp[2][0] = 97;
            adhesions_exp[2][1] = 230;
            adhesions_exp[3][0] = 96; 
            adhesions_exp[3][1] = 280;
            adhesions_exp[4][0] = 104;
            adhesions_exp[4][1] = 119;
            adhesions_exp[5][0] = 111;
            adhesions_exp[5][1] = 255;
            adhesions_exp[6][0] = 132; 
            adhesions_exp[6][1] = 118;
            adhesions_exp[7][0] = 143;
            adhesions_exp[7][1] = 351;
            adhesions_exp[8][0] = 205;
            adhesions_exp[8][1] = 89;
            adhesions_exp[9][0] = 209; 
            adhesions_exp[9][1] = 417;
            adhesions_exp[10][0] = 242;
            adhesions_exp[10][1] = 68;
            adhesions_exp[11][0] = 253;
            adhesions_exp[11][1] = 86;
            adhesions_exp[12][0] = 264; 
            adhesions_exp[12][1] = 466;
            adhesions_exp[13][0] = 312;
            adhesions_exp[13][1] = 132;
            adhesions_exp[14][0] = 319;
            adhesions_exp[14][1] = 369;
            adhesions_exp[15][0] = 338;
            adhesions_exp[15][1] = 141;
            adhesions_exp[16][0] = 343;
            adhesions_exp[16][1] = 323;
            adhesions_exp[17][0] = 352;
            adhesions_exp[17][1] = 159;
            adhesions_exp[18][0] = 358;
            adhesions_exp[18][1] = 298;
            adhesions_exp[19][0] = 369;
            adhesions_exp[19][1] = 271;
            adhesions_exp[20][0] = 399;
            adhesions_exp[20][1] = 180;
            adhesions_exp[21][0] = 427;
            adhesions_exp[21][1] = 222;
            adhesions_exp[22][0] = 465;
            adhesions_exp[22][1] = 195;
            */
            
            adhesions_exp[0][0] = 53; 
            adhesions_exp[0][1] = 161;
            adhesions_exp[1][0] = 79;
            adhesions_exp[1][1] = 115;
            adhesions_exp[2][0] = 104;
            adhesions_exp[2][1] = 119;
            adhesions_exp[3][0] = 132; 
            adhesions_exp[3][1] = 118; 
            adhesions_exp[4][0] = 205;
            adhesions_exp[4][1] = 89;
            adhesions_exp[5][0] = 242;
            adhesions_exp[5][1] = 68;
            adhesions_exp[6][0] = 253;
            adhesions_exp[6][1] = 86;
            adhesions_exp[7][0] = 312;
            adhesions_exp[7][1] = 132;
            adhesions_exp[8][0] = 338;
            adhesions_exp[8][1] = 141;
            adhesions_exp[9][0] = 352;
            adhesions_exp[9][1] = 159;
            adhesions_exp[10][0] = 399;
            adhesions_exp[10][1] = 180;
            adhesions_exp[11][0] = 465;
            adhesions_exp[11][1] = 195;
            adhesions_exp[12][0] = 427;
            adhesions_exp[12][1] = 222;
            adhesions_exp[13][0] = 369;
            adhesions_exp[13][1] = 271;
            adhesions_exp[14][0] = 358;
            adhesions_exp[14][1] = 298;
            adhesions_exp[15][0] = 343;
            adhesions_exp[15][1] = 323;
            adhesions_exp[16][0] = 319;
            adhesions_exp[16][1] = 369;
            adhesions_exp[17][0] = 264; 
            adhesions_exp[17][1] = 466;
            adhesions_exp[18][0] = 209; 
            adhesions_exp[18][1] = 417;
            adhesions_exp[19][0] = 143;
            adhesions_exp[19][1] = 351;
            adhesions_exp[20][0] = 96; 
            adhesions_exp[20][1] = 280;
            adhesions_exp[21][0] = 111;
            adhesions_exp[21][1] = 255;
            adhesions_exp[22][0] = 97;
            adhesions_exp[22][1] = 230;
            
            
            
            
            }
        else if(shape_number == 103){ // This is cell 303 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.
        
            //Cell 303 - also in the paper
            /*
            adhesions_exp[0][0] = 215; 
            adhesions_exp[0][1] = 145;
            adhesions_exp[1][0] = 164;
            adhesions_exp[1][1] = 333;
            adhesions_exp[2][0] = 222;
            adhesions_exp[2][1] = 171;
            adhesions_exp[3][0] = 172; 
            adhesions_exp[3][1] = 359;
            adhesions_exp[4][0] = 237;
            adhesions_exp[4][1] = 226;
            adhesions_exp[5][0] = 197;
            adhesions_exp[5][1] = 366;
            adhesions_exp[6][0] = 315; 
            adhesions_exp[6][1] = 147;
            adhesions_exp[7][0] = 278;
            adhesions_exp[7][1] = 285;
            adhesions_exp[8][0] = 375;
            adhesions_exp[8][1] = 88;
            adhesions_exp[9][0] = 394; 
            adhesions_exp[9][1] = 68;
            adhesions_exp[10][0] = 351;
            adhesions_exp[10][1] = 282;
            adhesions_exp[11][0] = 428;
            adhesions_exp[11][1] = 101;
            adhesions_exp[12][0] = 415; 
            adhesions_exp[12][1] = 148;
            adhesions_exp[13][0] = 452;
            adhesions_exp[13][1] = 283;
            adhesions_exp[14][0] = 471;
            adhesions_exp[14][1] = 264;
            
            */
            adhesions_exp[0][0] = 215; 
            adhesions_exp[0][1] = 145;
            adhesions_exp[1][0] = 315; 
            adhesions_exp[1][1] = 147;            
            adhesions_exp[2][0] = 375;
            adhesions_exp[2][1] = 88;
            adhesions_exp[3][0] = 394; 
            adhesions_exp[3][1] = 68;            
            adhesions_exp[4][0] = 428;
            adhesions_exp[4][1] = 101;            
            adhesions_exp[5][0] = 415; 
            adhesions_exp[5][1] = 148;
            adhesions_exp[6][0] = 471;
            adhesions_exp[6][1] = 264;
            adhesions_exp[7][0] = 452;
            adhesions_exp[7][1] = 283;
            adhesions_exp[8][0] = 351;
            adhesions_exp[8][1] = 282;
            adhesions_exp[9][0] = 278;
            adhesions_exp[9][1] = 285;
            adhesions_exp[10][0] = 197;
            adhesions_exp[10][1] = 366;
            adhesions_exp[11][0] = 172; 
            adhesions_exp[11][1] = 359;
            adhesions_exp[12][0] = 164;
            adhesions_exp[12][1] = 333;
            adhesions_exp[13][0] = 237;
            adhesions_exp[13][1] = 226;
            adhesions_exp[14][0] = 222;
            adhesions_exp[14][1] = 171;
            }
            
        else if(shape_number == 104){ // This is cell 312 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            //Cell 312 - also in the paper
            /*
            adhesions_exp[0][0] = 95; 
            adhesions_exp[0][1] = 364;
            adhesions_exp[1][0] = 107;
            adhesions_exp[1][1] = 418;
            adhesions_exp[2][0] = 127;
            adhesions_exp[2][1] = 397;
            adhesions_exp[3][0] = 263; 
            adhesions_exp[3][1] = 58;
            adhesions_exp[4][0] = 214;
            adhesions_exp[4][1] = 246;
            adhesions_exp[5][0] = 175;
            adhesions_exp[5][1] = 385;
            adhesions_exp[6][0] = 271; 
            adhesions_exp[6][1] = 84;
            adhesions_exp[7][0] = 258;
            adhesions_exp[7][1] = 131;
            adhesions_exp[8][0] = 247;
            adhesions_exp[8][1] = 179;
            adhesions_exp[9][0] = 193; 
            adhesions_exp[9][1] = 365;
            adhesions_exp[10][0] = 336;
            adhesions_exp[10][1] = 54;
            adhesions_exp[11][0] = 323;
            adhesions_exp[11][1] = 98;
            adhesions_exp[12][0] = 261; 
            adhesions_exp[12][1] = 333;
            adhesions_exp[13][0] = 356;
            adhesions_exp[13][1] = 323;
            adhesions_exp[14][0] = 363;
            adhesions_exp[14][1] = 58;
            adhesions_exp[15][0] = 331;
            adhesions_exp[15][1] = 228;
            adhesions_exp[16][0] = 338;
            adhesions_exp[16][1] = 255;
            adhesions_exp[17][0] = 354;
            adhesions_exp[17][1] = 308;
            adhesions_exp[18][0] = 374;
            adhesions_exp[18][1] = 290;
            */
            
//The old and incorrect order of the adhesions:
            /*
            adhesions_exp[0][0] = 95; 
            adhesions_exp[0][1] = 364;
            adhesions_exp[1][0] = 107;
            adhesions_exp[1][1] = 418;
            adhesions_exp[2][0] = 127;
            adhesions_exp[2][1] = 397;            
            adhesions_exp[3][0] = 175;
            adhesions_exp[3][1] = 385;
            adhesions_exp[4][0] = 193; 
            adhesions_exp[4][1] = 365;
            adhesions_exp[5][0] = 261; 
            adhesions_exp[5][1] = 333;
            adhesions_exp[6][0] = 356;
            adhesions_exp[6][1] = 323;
            adhesions_exp[7][0] = 354;
            adhesions_exp[7][1] = 308;
            adhesions_exp[8][0] = 374;
            adhesions_exp[8][1] = 290;
            adhesions_exp[9][0] = 338;
            adhesions_exp[9][1] = 255;
            adhesions_exp[10][0] = 331;
            adhesions_exp[10][1] = 228;
            adhesions_exp[11][0] = 323;
            adhesions_exp[11][1] = 98;
            adhesions_exp[12][0] = 363;
            adhesions_exp[12][1] = 58;
            adhesions_exp[13][0] = 336;
            adhesions_exp[13][1] = 54;
            adhesions_exp[14][0] = 271; 
            adhesions_exp[14][1] = 84;
            adhesions_exp[15][0] = 263; 
            adhesions_exp[15][1] = 58;
            adhesions_exp[16][0] = 258;
            adhesions_exp[16][1] = 131;
            adhesions_exp[17][0] = 247;
            adhesions_exp[17][1] = 179;
            adhesions_exp[18][0] = 214;
            adhesions_exp[18][1] = 246;
            }*/
            
            adhesions_exp[0][0] = 214;
            adhesions_exp[0][1] = 246;
            adhesions_exp[1][0] = 247;
            adhesions_exp[1][1] = 179;
            adhesions_exp[2][0] = 258;
            adhesions_exp[2][1] = 131;
            adhesions_exp[3][0] = 263; 
            adhesions_exp[3][1] = 58;
            adhesions_exp[4][0] = 271; 
            adhesions_exp[4][1] = 84;
            adhesions_exp[5][0] = 336;
            adhesions_exp[5][1] = 54;
            adhesions_exp[6][0] = 363;
            adhesions_exp[6][1] = 58;
            adhesions_exp[7][0] = 323;
            adhesions_exp[7][1] = 98;
            adhesions_exp[8][0] = 331;
            adhesions_exp[8][1] = 228;
            adhesions_exp[9][0] = 338;
            adhesions_exp[9][1] = 255;
            adhesions_exp[10][0] = 374;
            adhesions_exp[10][1] = 290;
            adhesions_exp[11][0] = 354;
            adhesions_exp[11][1] = 308;
            adhesions_exp[12][0] = 356;
            adhesions_exp[12][1] = 323;
            adhesions_exp[13][0] = 261; 
            adhesions_exp[13][1] = 333;
            adhesions_exp[14][0] = 193; 
            adhesions_exp[14][1] = 365;
            adhesions_exp[15][0] = 175;
            adhesions_exp[15][1] = 385;
            adhesions_exp[16][0] = 127;
            adhesions_exp[16][1] = 397; 
            adhesions_exp[17][0] = 107;
            adhesions_exp[17][1] = 418;
            adhesions_exp[18][0] = 95; 
            adhesions_exp[18][1] = 364;
        }

                        
        else if(shape_number == 105){ // This is cell 792 from the data
            //Make sure the anisotropy is 1 such that it does not go wrong with the size of the lattice.

            //Cell 792 - also in the paper
            /*
            adhesions_exp[0][0] = 125; 
            adhesions_exp[0][1] = 267;
            adhesions_exp[1][0] = 167;
            adhesions_exp[1][1] = 243;
            adhesions_exp[2][0] = 222;
            adhesions_exp[2][1] = 341;
            adhesions_exp[3][0] = 238; 
            adhesions_exp[3][1] = 219;
            adhesions_exp[4][0] = 278;
            adhesions_exp[4][1] = 342;
            adhesions_exp[5][0] = 293;
            adhesions_exp[5][1] = 269;
            adhesions_exp[6][0] = 307; 
            adhesions_exp[6][1] = 196;
            adhesions_exp[7][0] = 363;
            adhesions_exp[7][1] = 196;
            */
            adhesions_exp[0][0] = 125; 
            adhesions_exp[0][1] = 267;
            adhesions_exp[1][0] = 167;
            adhesions_exp[1][1] = 243;
            adhesions_exp[2][0] = 238; 
            adhesions_exp[2][1] = 219;
            adhesions_exp[3][0] = 307; 
            adhesions_exp[3][1] = 196;
            adhesions_exp[4][0] = 363;
            adhesions_exp[4][1] = 196;
            adhesions_exp[5][0] = 293;
            adhesions_exp[5][1] = 269;
            adhesions_exp[6][0] = 278;
            adhesions_exp[6][1] = 342;
            adhesions_exp[7][0] = 222;
            adhesions_exp[7][1] = 341;

            }
            
            else if(shape_number == 200){ // This is the test cell
                        

        }

                  
                                      
                        
                        
        else{
                printf("\nExperimental cell shape number not defined");
        }
                        
                        
        for(p=0;p<num_of_adhesions;p++){
                adhesions_exp_temp[(num_of_adhesions-1)-p][0]=adhesions_exp[p][0];
                adhesions_exp_temp[(num_of_adhesions-1)-p][1]=adhesions_exp[p][1];
        } //For inverting the order of the adhesion points.
                
        for(q=0;q<num_of_adhesions;q++){
            adhesions_exp[q][0]=adhesions_exp_temp[q][0];
            adhesions_exp[q][1]=adhesions_exp_temp[q][1];
        } //For inverting the order of the adhesion points.
                       
        for(k=0;k<num_of_adhesions;k++){
            if(adhesions_exp[k][0]>512){
                printf("\nThe following X adhesion is out of bound, with value \t%.10f\t%.10f", k, adhesions_exp[k][0]);    
            }
            if(adhesions_exp[k][1]>512){
                printf("\nThe following Y adhesion is out of bound, with value \t%.10f\t%.10f", k, adhesions_exp[k][1]);    
            }
            adhesions[k][0] = adhesions_exp[k][0]/num_of_pixels;
            adhesions[k][1] = adhesions_exp[k][1]/num_of_pixels;
                //printf("\nadhesions new\t%.10f\t%.10f", adhesions[k][0], adhesions[k][1]);
        } // Rescaling from num_of_pixels pixels to a box of size 1.

            mirroring = 0;  // 1 for horizontal, 2 for vertical mirroring
            
            for(n=0;n<num_of_adhesions;n++){
                if(mirroring==0){
                    adhesions_mirrored[n][0] = adhesions[n][0];
                    adhesions_mirrored[n][1] = adhesions[n][1];
                }
                else if(mirroring==1){
                    adhesions_mirrored[n][0] = 1.0 - adhesions[n][0];
                    adhesions_mirrored[n][1] = adhesions[n][1];
                    //printf("\nmirrored vertically");
                }
                else if(mirroring==2){
                    adhesions_mirrored[n][0] = adhesions[n][0];
                    adhesions_mirrored[n][1] = 1.0 - adhesions[n][1];
                    //printf("\nmirrored horizontally");
                }
                else if(mirroring==3){
                    adhesions_mirrored[n][0] = 1.0 - adhesions[n][0];
                    adhesions_mirrored[n][1] = 1.0 - adhesions[n][1];
                    //printf("\nmirrored horizontally+vertically");
                }
                else{
                    printf("\nError with mirroring of coordinates");
                }
            } // Mirroring if required.
            

            if(simulation_exp==0 || simulation_exp == 2){
                max_x = 0;
                max_y = 0;
                min_x = 1;
                min_y = 1;
            
                for(l=0;l<num_of_adhesions;l++){
                    if(adhesions_mirrored[l][0]<min_x){
                        min_x=adhesions_mirrored[l][0];                    
                    }
                    if(adhesions_mirrored[l][0]>max_x){
                        max_x=adhesions_mirrored[l][0]; 
                    }
                    if(adhesions_mirrored[l][1]<min_y){
                        min_y=adhesions_mirrored[l][1]; 
                    }
                    if(adhesions_mirrored[l][1]>max_y){
                        max_y=adhesions_mirrored[l][1]; 
                    }
                    //printf("\nmax_x,max_y,min_x,min_y\t%.10f\t%.10f\t%.10f\t%.10f", max_x, max_y, min_x, min_y);
                } //Finding furthest values for x and y
            
            
                size_x = max_x-min_x;
                size_y = max_y-min_y;
                bordersize = 0.18;
            
                if(min_x<min_y){
                    min_total=min_x;
                }
                else{
                    min_total=min_y;
                }

            
                if(size_x>size_y){
                    size_max=size_x;
                    //printf("\nscaled to X");
                }
                else{
                    size_max=size_y;  
                    //printf("\nscaled to Y");
                }
            
                for(m=0;m<num_of_adhesions;m++){
                    adhesions_scaled[m][0]=bordersize+((adhesions_mirrored[m][0]-min_x)*(1/size_max)*(1.0-2.0*bordersize));
                    adhesions_scaled[m][1]=bordersize+((adhesions_mirrored[m][1]-min_y)*(1/size_max)*(1.0-2.0*bordersize));

                //printf("\nadhesions_scaled\t%.10f\t%.10f", adhesions_scaled[m][0], adhesions_scaled[m][1]);
                }

            
            //printf("\nsize_x, size_y, size_max\t%.10f\t%.10f\t%.10f\t%.10f", size_x, size_y, size_max);   

                
                
                
                
                for(a=0;a<num_of_adhesions;a++){
                    adhesions_relative[a][0] = (adhesions_scaled[a][0] - center_x)*size; //Scale according to wanted size and move to origin for rotation.
                    adhesions_relative[a][1] = (adhesions_scaled[a][1] - center_y)*size;
                }
                // adhesions_exp, max_x, min_x, max_y, max_y, size_x, size_y, size_max, border, min_xy, max_xy

            
                for(b=0;b<num_of_adhesions;b++){
                    adhesions_rotated[b][0] = adhesions_relative[b][0]*cos(rotationpi) + adhesions_relative[b][1]*sin(rotationpi);
                    adhesions_rotated[b][1] = -adhesions_relative[b][0]*sin(rotationpi) + adhesions_relative[b][1]*cos(rotationpi);
                    //printf("\n%.10f\t%.10f", adhesions_rotated[b][0], adhesions_rotated[b][1]);
                }

                for(c=0;c<num_of_adhesions;c++){
                    adhesions_final[c][0] = adhesions_rotated[c][0]+center_x;
                    adhesions_final[c][1] = adhesions_rotated[c][1]+center_y;
                }
            }
            
            else if(simulation_exp==1){
                for(c=0;c<num_of_adhesions;c++){
                    adhesions_final[c][0] = adhesions[c][0];
                    adhesions_final[c][1] = adhesions[c][1];
                }
                bordersize=0.0;
                size_max=1.0;
            }

 
            
            for(j=0;j<num_of_adhesions;j++){
                for (i=0; i<num_of_edgepoint; i++){
                    if(j+1>=num_of_adhesions){                         //Make sure that the last edge loops back onto the first adhesion point.
                        k = 0;
                    }
                    else{
                        k = j+1;
                    }
                    dx = (adhesions_final[k][0]-adhesions_final[j][0])/(num_of_edgepoint-1.0);
                    dy = (adhesions_final[k][1]-adhesions_final[j][1])/(num_of_edgepoint-1.0);
                    vertex[j*num_of_edgepoint+i].x = adhesions_final[j][0]+dx*i;
                    vertex[j*num_of_edgepoint+i].y = adhesions_final[j][1]+dy*i;
                    
                    if(i==0 || i==num_of_edgepoint-1){                      //Defining the adhesions
                        vertex[j*num_of_edgepoint+i].adhesion =1;
                    }
                    else{
                        vertex[j*num_of_edgepoint+i].adhesion =0;
                    }
                }
            }
            
            if(simulation_exp==1 || simulation_exp==2){
            sigmaratio_scaled = sigmaratio_initial*((1.0-2.0*bordersize)/size_max/size);
            sigma=(num_of_pixels*size_of_pixel)/sigmaratio_scaled;
            alpha=((1.0-gamma_data)/gamma_data)*sigma;
            }
        }
        
        timestepbulkhelp = timestepbulkratio*rotational_friction*lattice_spacing*lattice_spacing/elasticity; //Calculate here the allowed timestepbulk. The ratio timestepbulk*elasticity/(rotational_friction*lattice_spacing*lattice_spacing) gives the update in phi coming from the neighboring bulk points. I want this to be smaller than or equal to some number, which I call timestepbulkratio.
        
        if(alpha+sigma>0){ // This loop was previously in init(), but had to be moved down since values of alpha and sigma are only correctly defined after init_boundary().
            time_stephelp = time_stepratio*num_of_adhesions*friction/((num_of_bounpoint-num_of_adhesions)*(sigma+alpha)); //Calculate the allowed time_step for the boundary. Here the number time_step*(of the order sigma + alpha)/friction gives the spacial update. I want that to bes smaller than the distance between two vertex points by a ratio time_stepratio. Here I approximate the distance between two vertices as 1/number of vertex per boundary, but this could be refined if necessary.
                    //printf("\n%.10f\n",time_stephelp);
            
            //Now I want to make sure that one of these time steps it the multiple of the other, such that I can make sure the times of the two processes actually match.
            if(timestepbulkhelp>time_stephelp){
                //timestepbulk = time_stephelp;
                //time_step = time_stephelp;
                printf("\nWarning: timestepbulkhelp>time_stephelp. The bulk is updated more often than necessary, this might be inefficient. timestepbulkhelp,time_stephelp = %.10f\t%.10f\n",timestepbulkhelp,time_stephelp);
            }
                //timestepbulk = timestepbulkhelp;
                //time_step = timestepbulkhelp;

        }
            timestepbulk = timestepbulkhelp;
            time_step = timestepbulkhelp;

        

        
        //printf("\n\n",timestepbulk);

}
/*******************************************************************/

void init_bulk()
{
    int i,j,k,m;
    double xvalue,yvalue,increment;
    
    increment = 1.0/(num_of_gridpoints_on_line-2*size_of_margin-1.0);
    
    for(i=0;i<num_of_gridpoints_on_line;i++){
        for(j=0;j<num_of_gridpoints_on_line_x;j++){
            grid[num_of_gridpoints_on_line_x*i+j].x = -(size_of_margin+0.0)*increment+increment*(j+0.0); //Make sure that the grid has a margin of size_of_margin points on either side of 0 and 1.
            grid[num_of_gridpoints_on_line_x*i+j].y = -(size_of_margin+0.0)*increment+increment*(i+0.0);
            
            
            
            
            
            
            
            // Included q1 and q2
                   
            
            
            
            
            
            
            grid[num_of_gridpoints_on_line_x*i+j].phi = 0.5*pi*(2*ran2(&seed)-1);
            //grid[num_of_gridpoints_on_line*i+j].phi =atan((-grid[num_of_gridpoints_on_line*i+j].x)/(grid[num_of_gridpoints_on_line*i+j].y+0.00001));
            if(grid[num_of_gridpoints_on_line_x*i+j].phi > 0.5*pi){
                grid[num_of_gridpoints_on_line_x*i+j].phi = grid[num_of_gridpoints_on_line_x*i+j].phi - pi;
            }
            else if(grid[num_of_gridpoints_on_line_x*i+j].phi <= -0.5*pi){
                grid[num_of_gridpoints_on_line_x*i+j].phi = grid[num_of_gridpoints_on_line_x*i+j].phi + pi;
            }
            grid[num_of_gridpoints_on_line_x*i+j].s = ran2(&seed);
            grid[num_of_gridpoints_on_line_x*i+j].loworder = 0;
            grid[num_of_gridpoints_on_line_x*i+j].q1 = grid[num_of_gridpoints_on_line_x*i+j].s*cos(2.0*grid[num_of_gridpoints_on_line_x*i+j].phi)*0.5;
            grid[num_of_gridpoints_on_line_x*i+j].q2 = grid[num_of_gridpoints_on_line_x*i+j].s*sin(2.0*grid[num_of_gridpoints_on_line_x*i+j].phi)*0.5;
        }
    }
    
    //the following procedure creates a data structure lines that defines an array of lines in between the array points, which are used later for determining which points are part of the bulk (used in determine_bulkedge and determine_bulk).
    //Additionally, this data structure is used in order to store the angle of the boundary between two grid points, which is required for implementing anchoring boundary conditions.
    
    //The first forloop fills the first half of the datastructure lines with vertical lines.
    for(k=0;k<num_of_lines;k++){
        lines[k].edgecross = 0;
        xvalue = -(size_of_margin+0.0)*increment+(k%num_of_gridpoints_on_line_x)*increment;
        yvalue = -(size_of_margin+0.0)*increment+increment*(k-k%num_of_gridpoints_on_line_x)/(num_of_gridpoints_on_line_x+0.0);
        lines[k].x = xvalue;
        lines[k].y = yvalue;
    }
    
    
    //The second forloop fills the second half with horizontal lines. Note that both horizontal and vertical lines are counted from left to right first, and then from bottom to top.
    for(k=num_of_lines;k<num_of_lines+(num_of_gridpoints_on_line_x-1)*(num_of_gridpoints_on_line);k++){
        lines[k].edgecross = 0;
        xvalue = -(size_of_margin+0.0)*increment+((k-num_of_lines)%num_of_gridpoints_on_line_x)*increment;
        yvalue = -(size_of_margin+0.0)*increment+increment*((k-num_of_lines)-(k-num_of_lines)%num_of_gridpoints_on_line_x)/(num_of_gridpoints_on_line_x+0.0);
        lines[k].x = xvalue;
        lines[k].y = yvalue;
    }
        
}

/*******************************************************************/
//I define a separate function that determines the local angle of the bulk theta at a given point on the cell edge. In this version I take the value of the bulk closest to the edgepoint.
void find_orientation()
{
        int i,gridpointnumber1,gridpointnumber2,gridpointnumber3,gridpointnumber4,gridpoint1x,gridpoint1y,gridpoint2x,gridpoint2y,gridpoint3x,gridpoint3y,gridpoint4x,gridpoint4y;
        int bulktrack[4];
        double xrescaled,yrescaled,increment,dhelp1,dhelp2,dmin,distance1,distance2,distance3,distance4;
        
        increment = 1.0/(num_of_gridpoints_on_line-2.0*size_of_margin-1.0);
        
        for(i=0;i<num_of_bounpoint;i++){
            //First determine which bulk points are around this vertex point. 
            xrescaled = vertex[i].x/increment + size_of_margin;
            yrescaled = vertex[i].y/increment + size_of_margin;
            gridpoint1x= floor(xrescaled);
            gridpoint1y= floor(yrescaled);
            gridpoint2x = gridpoint1x;
            gridpoint2y = gridpoint1y+1;
            gridpoint3x = gridpoint1x+1;
            gridpoint3y = gridpoint1y+1;
            gridpoint4x = gridpoint1x+1;
            gridpoint4y = gridpoint1y;
            
            gridpointnumber1 = num_of_gridpoints_on_line*gridpoint1y+gridpoint1x;
            gridpointnumber2 = num_of_gridpoints_on_line*gridpoint2y+gridpoint2x;
            gridpointnumber3 = num_of_gridpoints_on_line*gridpoint3y+gridpoint3x;
            gridpointnumber4 = num_of_gridpoints_on_line*gridpoint4y+gridpoint4x;
            
            
            //Next determine which point is both in the bulk (bulk ==1) and is the closest to the vertex point. I define the distance very large if bulk ==0 so that I will not find that point.
            if(grid[gridpointnumber1].bulk==0){
                distance1=1000;
            }
            else{
                distance1=(gridpoint1x-xrescaled)*(gridpoint1x-xrescaled)+(gridpoint1y-yrescaled)*(gridpoint1y-yrescaled);
            }
            if(grid[gridpointnumber2].bulk==0){
                distance2=1000;
            }
            else{
                distance2=(gridpoint2x-xrescaled)*(gridpoint2x-xrescaled)+(gridpoint2y-yrescaled)*(gridpoint2y-yrescaled);
            }
            if(grid[gridpointnumber3].bulk==0){
                distance3=1000;
            }
            else{
                distance3=(gridpoint3x-xrescaled)*(gridpoint3x-xrescaled)+(gridpoint3y-yrescaled)*(gridpoint3y-yrescaled);
            }
            if(grid[gridpointnumber4].bulk==0){
                distance4=1000;
            }
            else{
                distance4=(gridpoint4x-xrescaled)*(gridpoint4x-xrescaled)+(gridpoint4y-yrescaled)*(gridpoint4y-yrescaled);
            }
            
            //Now find which distance is the smallest:
            if(distance1<=distance2){
                dhelp1 = distance1;
            }
            else{
                dhelp1 = distance2;
            }
            if(dhelp1<=distance3){
                dhelp2 = dhelp1;
            }
            else{
                dhelp2 = distance3;
            }
            if(dhelp2<=distance4){
                dmin = dhelp2;
            }
            else{
                dmin = distance4;
            }
            

            
            //Which point did belong to that distance? Assign that local angle to the vertex.
            
            
            if(dmin == 1000){
                vertex[i].theta = NAN;
                vertex[i].q1 = 0.0;
                vertex[i].q2 = 0.0;
                vertex[i].s = 0.0;
            }
            else if(dmin == distance1){
                vertex[i].theta = grid[gridpointnumber1].phi;
                vertex[i].q1 = grid[gridpointnumber1].q1;
                vertex[i].q2 = grid[gridpointnumber1].q2;
                vertex[i].s = grid[gridpointnumber1].s;
            }
            else if(dmin == distance2){
                vertex[i].theta = grid[gridpointnumber2].phi;
                vertex[i].q1 = grid[gridpointnumber2].q1;
                vertex[i].q2 = grid[gridpointnumber2].q2;
                vertex[i].s = grid[gridpointnumber2].s;
            }
            else if(dmin == distance3){
                vertex[i].theta = grid[gridpointnumber3].phi;
                vertex[i].q1 = grid[gridpointnumber3].q1;
                vertex[i].q2 = grid[gridpointnumber3].q2;
                vertex[i].s = grid[gridpointnumber3].s;
            }
            else if(dmin == distance4){
                vertex[i].theta = grid[gridpointnumber4].phi;
                vertex[i].q1 = grid[gridpointnumber4].q1;
                vertex[i].q2 = grid[gridpointnumber4].q2;
                vertex[i].s = grid[gridpointnumber4].s;
            }
            else{
                printf("Error: one of the four distances should be the minimum.");
            }
        }
    
}

/*******************************************************************/
//In find_orientation() I found for each boundary vertex the local bulk value of the angle. Here I'll use that bulk value to influence the pull on the boundary.
void update_boundary()
{        
        //here we update the shape of the boundary and the positions of the mesh points.
    
        double dx, dy, x_new, y_new,forcemag,forcex,forcey,dr_i,dr_iplus,phi_j,y,drx_i,drx_inorm,drx_iplus,drx_iplusnorm,dry_i,dry_inorm,dry_iplus,dry_iplusnorm,kappatimesnx,kappatimesny;
        //Define an auxiliary array here where we can store the new positions of the edgepoints:
        double phi_i,phi_imin,delta_phi_test,av_phi,delta_phi,adott,adotnprelim,angle_i,angle_iplus,anglediffhelp,anglediff,adotni,adotnimin,adotnav,adotti,adottimin,adottav,correction,minimum,sumax,sumay,anormx,anormy,anorm;
	double auxarray[num_of_bounpoint][2];
	long i,j,k,l,m,n;
        double lambdahelp[num_of_bounpoint];
        double adotn[num_of_bounpoint];
        double normalx[num_of_bounpoint];
        double normaly[num_of_bounpoint];
        double sumaxplus,sumayplus,anormplus,anormplusx,anormplusy,adotnavplus,adottavplus;
        
        
        
        //Before starting the big for-loop, I calculate the tensions in the adhesion points 0:
        for (j=0;j<num_of_bounpoint;j=j+num_of_edgepoint){
           
           lambdahelp[j] = lambda;
        }
        
        //In this new version, I calculate the value of lambda0 later by taking the minimum value of lambda on a given arc as lambda-.
        
        
        
        
	for (i=0; i<num_of_bounpoint; i++){ 
                
                //The adhesion points should never be displaced, so we only update if the adhesion value is 0.
            
                //if(i%num_of_edgepoint!=0){
                if(vertex[i].adhesion ==0){
                
                //First we calculate the absolute distances, x and y differences between points i and i-1 and between points i+1 and i, and x and y components of unit vectors along the two lines.
                drx_i=vertex[i].x-vertex[i-1].x; 
                dry_i=vertex[i].y-vertex[i-1].y;
                
                //This shouldn't matter anymore with the doubly defined adhesion points.
                if(i==num_of_bounpoint-1){                  //This otherwise goes wrong for the last boundary point, since it's neighbour is not (0,1) but a point in the bulk.
                    drx_iplus=vertex[0].x-vertex[i].x;
                    dry_iplus=vertex[0].y-vertex[i].y;
                }
                else{
                    drx_iplus=vertex[i+1].x-vertex[i].x;
                    dry_iplus=vertex[i+1].y-vertex[i].y;
                }
                dr_i=sqrt(drx_i*drx_i+dry_i*dry_i);
                dr_iplus=sqrt(drx_iplus*drx_iplus+dry_iplus*dry_iplus);
                drx_inorm=drx_i/dr_i;
                dry_inorm=dry_i/dr_i;
                drx_iplusnorm=drx_iplus/dr_iplus;
                dry_iplusnorm=dry_iplus/dr_iplus;
                
                //Then we calculate the product kappa*n by "fitting" a circle through the three points and calculating kappa and n. Then we calculate kappa and the components of the normal vector. It is important for the sign of kappa that the vertexs are defined in the direction corresponding to the analytic theory (anti-clockwise around the circle/ellipse, so clockwise around the cell)
                
                kappatimesnx = 2*(drx_inorm*dry_iplusnorm-drx_iplusnorm*dry_inorm)*(dr_i*dry_iplusnorm+dr_iplus*dry_inorm)/((drx_i+drx_iplus)*(drx_i+drx_iplus)+(dry_i+dry_iplus)*(dry_i+dry_iplus));
                kappatimesny = 2*(drx_inorm*dry_iplusnorm-drx_iplusnorm*dry_inorm)*(-dr_i*drx_iplusnorm-dr_iplus*drx_inorm)/((drx_i+drx_iplus)*(drx_i+drx_iplus)+(dry_i+dry_iplus)*(dry_i+dry_iplus));
                
                //Now I need to find the sign of the value of kappa, before calculating kappa itself. I do this by realizing that if ny has the same sign as tx, then kappa is positive. If ny has an opposite sign to tx, then kappa is negative. So the same goes for kappatimesny: if kappatimesny has the same sign as tx, then kappa is positive. If kappatimesny has an opposite sign to tx, then kappa is negative. Furthermore, tx is equal to the above defined. 
                //Don't use this, but first try something simpler:
                
                //To determine the sign of the curvature, we determine the angle between two consecutive vectors. First the angle of the first vector:
                if(drx_i==0){
                    if(dry_i>0){
                        angle_i = 0.5*pi;
                    }
                    else{
                        angle_i = -0.5*pi;
                    }
                }
                else if(drx_i>0){
                    angle_i = atan(dry_i/drx_i);
                }
                else{
                    angle_i = atan(dry_i/drx_i)+pi;
                }
                //Then the angle of the second vector:
                if(drx_iplus==0){
                    if(dry_iplus>0){
                        angle_iplus = 0.5*pi;
                    }
                    else{
                        angle_iplus = -0.5*pi;
                    }
                }
                else if(drx_iplus>0){
                    angle_iplus = atan(dry_iplus/drx_iplus);
                }
                else{
                    angle_iplus = atan(dry_iplus/drx_iplus)+pi;
                }
                //Then calculate the difference in angles. Since this difference should be small, I make sure it is always between - pi and pi. 
                anglediffhelp = angle_iplus - angle_i;
                if(anglediffhelp>pi){
                    anglediff = anglediffhelp - 2*pi;
                }
                else if(anglediffhelp<-pi){
                    anglediff = anglediffhelp + 2*pi;
                }
                else{
                    anglediff = anglediffhelp;
                }
                //Now, if the difference in angles is positive, then the curvature is negative and vice versa:
                if(anglediff>=0){
                    vertex[i].kappa = -sqrt(kappatimesnx*kappatimesnx+kappatimesny*kappatimesny);
                }
                else{
                    vertex[i].kappa = sqrt(kappatimesnx*kappatimesnx+kappatimesny*kappatimesny);
                }
                
                if(vertex[i].kappa!=0){
                    normalx[i] = -kappatimesnx/vertex[i].kappa;             //I store normalx and normaly in an array because I need their values outside this forloop.
                    normaly[i] = -kappatimesny/vertex[i].kappa;
                }
                else{                                         //In case the curvature is 0, we calculate the normal vector by just taking that of one of the two segments.
                    normalx[i] = dry_inorm;
                    normaly[i] = -drx_inorm;
                }                
                
                // Before calculating forces, we calculate the value of the line tension for all edge points, given the tension in adhesion point 0 defined above:
                
                //Calculate here the average pull on a rod, by simpling adding the two local directions in a vectorial way. Take as the normal direction the normal direction of the rod between vertex i-1 and i (so normalx = dry_inorm and normaly = -drx_inorm). The result of the adding of two vectors should be normalized such that aj is still a unit vector.
                //First calculate average a. 
                
                double q1avg,q2avg,thetaavg,savg;
                double q1avgplus,q2avgplus,thetaavgplus,savgplus;
                
                q1avg=(vertex[i].q1+vertex[i-1].q1)/2;
                q2avg=(vertex[i].q2+vertex[i-1].q2)/2;
                savg=2.0*sqrt(q1avg*q1avg+q2avg*q2avg);
                            
                if(savg!=0.0){
                    if(q1avg > 0){
                    thetaavg = 0.5*atan(q2avg/q1avg);
                    }
                    else if(q1avg < 0){
                        if(q2avg > 0){
                        thetaavg = 0.5*(atan(q2avg/q1avg)+pi);
                        }
                        else if(q2avg < 0){
                            thetaavg = 0.5*(atan(q2avg/q1avg)-pi);
                        }
                        else{
                            thetaavg = 0.5*pi;
                        }
                    }
                    else{
                        if(q2avg > 0){
                            thetaavg = 0.25*pi;
                        }
                        else{
                            thetaavg = 0.75*pi;
                        }
                    }
                    anormx = cos(thetaavg);
                    anormy = sin(thetaavg);
                    adotnav = anormx*dry_inorm+anormy*(-drx_inorm);
                    adottav = anormx*drx_inorm+anormy*dry_inorm;
            
                
                    //Now using all this, calculate a temporary lambda file for all the lambda values.
                    lambdahelp[i] = lambdahelp[i-1]-alpha*savg*adotnav*adottav*dr_i;
                }
                else{
                    //printf("Error:S is zero for vertex calculation.\n");
                    
                    lambdahelp[i] = lambdahelp[i-1];
                }
                
                
                //Now also store lambda in the ending adhesion point:
                
                if(vertex[i+1].adhesion ==1){
                    
                    q1avgplus=(vertex[i+1].q1+vertex[i].q1)/2;
                    q2avgplus=(vertex[i+1].q2+vertex[i].q2)/2;
                    savgplus=2.0*sqrt(q1avgplus*q1avgplus+q2avgplus*q2avgplus);
                            
                    if(savgplus!=0.0){
                        if(q1avgplus > 0){
                        thetaavgplus = 0.5*atan(q2avgplus/q1avgplus);
                        }
                        else if(q1avgplus < 0){
                            if(q2avgplus > 0){
                            thetaavgplus = 0.5*(atan(q2avgplus/q1avgplus)+pi);
                            }
                            else if(q2avgplus < 0){
                                thetaavgplus = 0.5*(atan(q2avgplus/q1avgplus)-pi);
                            }
                            else{
                                thetaavgplus = 0.5*pi;
                            }
                        }
                        else{
                            if(q2avgplus > 0){
                                thetaavgplus = 0.25*pi;
                            }
                            else{
                                thetaavgplus = 0.75*pi;
                            }
                        }
                        anormplusx = cos(thetaavgplus);
                        anormplusy = sin(thetaavgplus);
                        adotnavplus = anormplusx*dry_iplusnorm+anormplusy*(-drx_iplusnorm);
                        adottavplus = anormplusx*drx_iplusnorm+anormplusy*dry_iplusnorm;
                        lambdahelp[i+1] = lambdahelp[i]-alpha*savgplus*adotnavplus*adottavplus*dr_iplus;
                    }
                    else{
                        lambdahelp[i+1] = lambdahelp[i];
                                        //printf("Error:S is zero (adhesion point) for vertex calculation.\n");
                    
                    }
                
                                

                }
                
                //Before I finish the first forloop, I use the calculationgs above to store adotn, in a temporary array. adotn is the dot product between a and n at the vertex point itself. The sign of the dot product does not matter because it is squared:
                    if(vertex[i].s==0.0){
                        adotn[i]=0.0;
                    }
                    else{
                        adotn[i] = cos(vertex[i].theta)*normalx[i]+sin(vertex[i].theta)*normaly[i];
                    }
                }
        }
        
        //Now before I proceed on calculating forces, I now first correct the values of lambda by identifying the minimum of all lambda values in an arc as lambda-.
        for(m=0;m<num_of_bounpoint;m=m+num_of_edgepoint){
                minimum = lambdahelp[m];
                for(n=1;n<num_of_edgepoint;n++){
                        if(lambdahelp[m+n]<minimum){
                            minimum = lambdahelp[m+n];
                        }
                }
                correction = lambda - minimum;
                for(n=0;n<num_of_edgepoint;n++){
                    vertex[m+n].lam = lambdahelp[m+n]+correction;
                }
        }
                        
                
        //Using the correct values of lambda, I then calculate the magnitude of the force and apply that to the different directions.
        for(l=0;l<num_of_bounpoint;l++){
            
            //if(l%num_of_edgepoint!=0){
            if(vertex[l].adhesion==0){
                forcemag = vertex[l].kappa*vertex[l].lam+sigma+alpha*vertex[l].s*adotn[l]*adotn[l];
		forcex = forcemag*normalx[l];
		forcey = forcemag*normaly[l];
                
                //Finally, I calculate the displacement using force = some constant times velocity, so 
                
                dx = time_step*forcex/friction;
                dy = time_step*forcey/friction;
                if(dx>1.0/lattice_spacing||dy>1.0/lattice_spacing){
                    printf("Error: dx or dy is larger than 1/lattice_spacing.");
                }
		
		// Construct a new position and store this is an auxiliary array.
		
		auxarray[l][0] = dx+vertex[l].x; 
		auxarray[l][1] = dy+vertex[l].y;
            }
            else{
                auxarray[l][0] = vertex[l].x;
                auxarray[l][1] = vertex[l].y;
                    /*if(i%num_of_edgepoint==num_of_edgepoint-1){
                        vertex[i].lam = vertex[i-1].lam+alpha*(vertex[i].x-vertex[i-1].x)*(vertex[i].y-vertex[i-1].y)/sqrt((vertex[i].x-vertex[i-1].x)*(vertex[i].x-vertex[i-1].x)+(vertex[i].y-vertex[i-1].y)*(vertex[i].y-vertex[i-1].y));
                    }*/
                }
        }
             
        
	for (j=0; j<num_of_bounpoint;j++){  //Now update the new coordinates by copying the auxiliary array to the edge coordinates.
		vertex[j].x = auxarray[j][0];
		vertex[j].y = auxarray[j][1];
			
	}        
}
/*******************************************************************/
//Using the shape, calculate the net forces on the adhesion points.
void find_adhesion_forces()
{
    long counter,i;
    double magnitude,deltax,deltay,deltar;
    counter =0;
    total_force=0.0;
    total_force_x=0.0;
    total_force_y=0.0;
    
    for(i = 0; i<num_of_bounpoint;i++){
        if(vertex[i].adhesion==1){
            magnitude = vertex[i].lam;
            forces[counter].x = vertex[i].x;
            forces[counter].y = vertex[i].y;
            if(vertex[i+1].adhesion==1 || i==num_of_bounpoint-1){ //All adhesion points 1
                deltax = -vertex[i].x+vertex[i-1].x;
                deltay = -vertex[i].y+vertex[i-1].y;
                deltar = sqrt(deltax*deltax+deltay*deltay);
                forces[counter].fx = magnitude*deltax/deltar;
                forces[counter].fy = magnitude*deltay/deltar; 
                forces[counter].type = 1;;
            }
            else if(vertex[i-1].adhesion==1 || i==0){           //All adhesion points 0
                deltax = vertex[i+1].x-vertex[i].x;
                deltay = vertex[i+1].y-vertex[i].y;
                deltar = sqrt(deltax*deltax+deltay*deltay);
                forces[counter].fx = magnitude*deltax/deltar;
                forces[counter].fy = magnitude*deltay/deltar;
                forces[counter].type = 0;

            }
            total_force_x=total_force_x+forces[counter].fx;
            total_force_y=total_force_y+forces[counter].fy;
            if (exporting == 1){
            //printf("\nTOTAL FORCE X Y counter %.10f\t%.10f \t %d",  total_force_x,  total_force_y,counter);
            };
            counter = counter +1;
        }
    }
    total_force=total_force_x+total_force_y;
}

/*******************************************************************/
//In this function we determine which lines in the data structure "lines" cross the boundary. Those define the edge of the bulk. This data is used in determine_bulk to determine which grid points are part of the bulk.
//To implement anchoring boundary conditions, I add the functionality here that the local angle of the boundary crossing the line is also stored in the array lines as variable psi.
//to implement anchoring boundary conditions, this function now also determines when horizontal lines are crossed, although that is not necessary to calculate which points are in the bulk.
void determine_bulkedge()
{ 
    int i,j,x,y,yline,xline;
    long xstart,xstop,ystart,ystop;
    double x_i,y_i,x_iplus,y_iplus,yvalue,xvalue,increment;
    double psi_i;
    
    increment = 1.0/(num_of_gridpoints_on_line-2*size_of_margin-1.0);
    
     for(j=0;j<(num_of_gridpoints_on_line_x)*(num_of_gridpoints_on_line-1)+(num_of_gridpoints_on_line_x-1)*(num_of_gridpoints_on_line);j++){         //Approx 2*numlines, but for anisotropic case there is a different number of horizontal and vertical lines.
         lines[j].edgecross=0;
         lines[j].psi = 200;
     }
     for(i=0;i<num_of_bounpoint-1;i++){
         
        x_i=vertex[i].x;
        y_i=vertex[i].y;
        
        if(i==num_of_bounpoint-1){     
            x_iplus = vertex[0].x;
            y_iplus = vertex[0].y;
        }
        else{
            x_iplus = vertex[i+1].x;
            y_iplus = vertex[i+1].y;
        }
        
        //For the remainder of this function, only do the following if i and i+1 are not infact the same point (which occurs at the adhesions):
        //if(x_i != x_iplus && y_i != y_iplus){
        
            //Next, determine the value of psi_i, the angle that should be stored in the linestructure.
            if(x_i==x_iplus){
                psi_i = 0.5*pi; //-0.5 pi is not interesting since I require the angle to be between -0.5 pi and 0.5pi.
            }
            else{
                psi_i = atan((y_iplus-y_i)/(x_iplus-x_i));
            }
            
            //The next if and later else statement are for determining when the boundary line crossed a vertical grid line:
            
            
            if(x_iplus>x_i){
                xstart=floor((x_i+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);          //The grid is has a margin, and is not just from 0 to 1.
                if(xstart<0 || xstart >num_of_gridpoints_on_line_x-1){
                    printf("Error: xstart should be between 0 and the number of gridpoints on a line in x-direction.");
                }
                xstop=floor((x_iplus+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);
                if(xstop<0 || xstop >num_of_gridpoints_on_line_x-1){
                    printf("Error: xstop should be between 0 and the number of gridpoints on a line in x-direction.");
                }
                if(xstart<xstop){
                    for(x=xstart;x<xstop;x++){
                        yvalue = y_i + (y_iplus-y_i)*(((x-(size_of_margin+0.0))/(num_of_gridpoints_on_line-2*(size_of_margin+0.0)-1.0))-x_i)/(x_iplus-x_i);     //In this expression I should not simply use x, because x is defined here as an integer. I should use the coordinate corresponding to x.
                        if(yvalue<-0.01 ||yvalue>1.0+0.01){
                            printf("Error: y should be between 0 and 1. Message 1.");
                            printf("\t%.10f\n",yvalue);
                        }
                        yline = floor((yvalue+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0)));          //The grid is has a margin, and is not just from -1 to 1.
                        if(yline*num_of_gridpoints_on_line_x+x>num_of_lines-1){
                            printf("Error: the num_of_lines array seems to be too small.");
                        }
                        //printf("\t%i\t%i\t%i\n",yline,yline*num_of_gridpoints_on_line+x,num_of_lines);
                        
                        
                        if(lines[yline*num_of_gridpoints_on_line_x+x].edgecross == 0){
                            lines[yline*num_of_gridpoints_on_line_x+x].edgecross = 1;
                            lines[yline*num_of_gridpoints_on_line_x+x].psi = psi_i;
                        }
                        else{
                            lines[yline*num_of_gridpoints_on_line_x+x].edgecross = 0;     //In the unlikely case that two edgerods pass through the same gridline, the value must be reset to 0.
                            lines[yline*num_of_gridpoints_on_line_x+x].psi = 200;
                        }
                    }
                    
                }
                else if (xstart>xstop){
                    printf("Error: xstart cannot be bigger than xstop.");
                }
            }
            else if (x_iplus<x_i){
                xstart=floor((x_iplus+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);       //The grid is has a margin, and is not just from -1 to 1.
                if(xstart<0 || xstart >num_of_gridpoints_on_line_x-1){
                    printf("Error: xstart should be between 0 and the number of gridpoints on a line in x-direction.");
                }
                xstop=floor((x_i+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);
                if(xstop<0 || xstop >num_of_gridpoints_on_line_x-1){
                    printf("Error: xstop should be between 0 and the number of gridpoints on a line in x-direction.");
                }
                if(xstart<xstop){
                    for(x=xstart;x<xstop;x++){
                        yvalue = y_i + (y_iplus-y_i)*(((x-(size_of_margin+0.0))/(num_of_gridpoints_on_line-2*(size_of_margin+0.0)-1.0))-x_i)/(x_iplus-x_i);
                        if(yvalue<-0.01 ||yvalue>1.0+ 0.01){
                            printf("Error: y should be between 0 and 1. Message 2.");
                            printf("\t%.10f\t%.10f\t%.10f\t%i\t%.10f\t%.10f\n",yvalue,y_i,y_iplus,x,x_i,x_iplus);
                        }
                        yline = floor((yvalue+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0)));   //The grid is has a margin, and is not just from -1 to 1.
                        if(yline*num_of_gridpoints_on_line_x+x>num_of_lines-1){
                            printf("Error: the num_of_lines array seems to be too small.");
                        }
                        //printf("%i\t%i\t%i\n",yline,yline*num_of_gridpoints_on_line+x,num_of_lines);
                                        
                        if(lines[yline*num_of_gridpoints_on_line_x+x].edgecross == 0){
                            lines[yline*num_of_gridpoints_on_line_x+x].edgecross = 1;
                            lines[yline*num_of_gridpoints_on_line_x+x].psi = psi_i;
                        }
                        else{
                            lines[yline*num_of_gridpoints_on_line_x+x].edgecross = 0;     //In the unlikely case that two edgerods pass through the same gridline, the value must be reset to 0.
                            lines[yline*num_of_gridpoints_on_line_x+x].psi = 200;
                        }
                        
                    }
                }
                else if (xstart>xstop){
                    printf("Error: xstart cannot be bigger than xstop.");
                }
            }
            
            //Next, I calculate whether the boundary line crosses a horizontal gridline. 
            
            if(y_iplus>y_i){
                ystart=floor((y_i+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);          //The grid is has a margin, and is not just from -1 to 1.
                if(ystart<0 || ystart >num_of_gridpoints_on_line-1){
                    printf("Error: ystart should be between 0 and the number of gridpoints on a line.");
                }
                ystop=floor((y_iplus+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);
                if(ystop<0 || ystop >num_of_gridpoints_on_line-1){
                    printf("Error: ystop should be between 0 and the number of gridpoints on a line.");
                }
                if(ystart<ystop){
                    for(y=ystart;y<ystop;y++){
                        xvalue = x_i + (x_iplus-x_i)*(((y-(size_of_margin+0.0))/(num_of_gridpoints_on_line-2*(size_of_margin+0.0)-1.0))-y_i)/(y_iplus-y_i);     //In this expression I should not simply use y, because y is defined here as an integer. I should use the coordinate corresponding to y.
                        if(xvalue<-0.01 ||xvalue>1.0*anisotropy+0.01){
                            printf("Error: x should be between 0 and anisotropy. Message 1.");
                            printf("\t%.10f\n",yvalue);
                        }
                        xline = floor((xvalue+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0)));          //The grid is has a margin, and is not just from -1 to 1.
                        if(y*num_of_gridpoints_on_line_x+xline>num_of_lines-1){
                            printf("Error: the num_of_lines array seems to be too small.");
                        }
                        //printf("\t%i\t%i\t%i\n",xline,xline*num_of_gridpoints_on_line+y,num_of_lines);
                        if(lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross == 0){
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross = 1;
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].psi = psi_i;
                        }
                        else{
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross = 0;     //In the unlikely case that two edgerods pass through the same gridline, the value must be reset to 0.
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].psi = 200;
                        }
                    }
                }
                else if (ystart>ystop){
                    printf("Error: ystart cannot be bigger than ystop.");
                }
            }
            else if (y_iplus<y_i){
                ystart=floor((y_iplus+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);       //The grid is has a margin, and is not just from -1 to 1.
                if(ystart<0 || ystart >num_of_gridpoints_on_line-1){
                    printf("Error: ystart should be between 0 and 100.");
                }
                ystop=floor((y_i+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0))+1.0);
                if(ystop<0 || ystop >num_of_gridpoints_on_line-1){
                    printf("Error: ystop should be between 0 and 100.");
                }
                if(ystart<ystop){
                    for(y=ystart;y<ystop;y++){
                        xvalue = x_i + (x_iplus-x_i)*(((y-(size_of_margin+0.0))/(num_of_gridpoints_on_line-2*(size_of_margin+0.0)-1.0))-y_i)/(y_iplus-y_i);
                        if(xvalue<-0.01 ||xvalue>1.0*anisotropy + 0.01){
                            printf("Error: x should be between 0 and anisotropy. Message 2.");
                            printf("\t%.10f\t%.10f\t%.10f\t%i\t%.10f\t%.10f\n",yvalue,y_i,y_iplus,x,x_i,x_iplus);
                        }
                        xline = floor((xvalue+(size_of_margin+0.0)*increment)*(num_of_gridpoints_on_line-1.0-2*(size_of_margin+0.0)));   //The grid is has a margin, and is not just from -1 to 1.
                        if(y*num_of_gridpoints_on_line_x+xline>num_of_lines-1){
                            printf("Error: the num_of_lines array seems to be too small.");
                        }
                        //printf("%i\t%i\t%i\n",yline,yline*num_of_gridpoints_on_line+x,num_of_lines);
                        if(lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross == 0){
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross = 1;
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].psi = psi_i;
                        }
                        else{
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].edgecross = 0;     //In the unlikely case that two edgerods pass through the same gridline, the value must be reset to 0.
                            lines[num_of_lines+y*num_of_gridpoints_on_line_x+xline].psi = 200;
                        }
                    }
                }
                else if (ystart>ystop){
                    printf("Error: ystart cannot be bigger than xstop.");
                }
            }
            
            
            
        
        }
        //Finally, check if all values of the angle psi are within the correct range:
        for(j=0;j<(num_of_gridpoints_on_line_x)*(num_of_gridpoints_on_line-1)+(num_of_gridpoints_on_line_x-1)*(num_of_gridpoints_on_line);j++){
            if(lines[j].edgecross==1){
                if(lines[j].psi > 0.5*pi || lines[j].psi <= -0.5*pi){
                    printf("Error: the angle psi is out of bound.");
                }
            }
        }
     
            

}



/*******************************************************************/
//In this function I used the function determine_bulkedge to find which points are part of the bulk. The idea: I'll walk from bottom to top for every column and change the value of bulk when I pass the boundary.

void determine_bulk()
{
    int i,j,bulkswap;
    
    for(i=0;i<num_of_gridpoints_on_line_x;i++){
        grid[i].bulk=0; //Inital value = 0 (not in the bulk)
        for(j=0;j<num_of_gridpoints_on_line-1;j++){
            bulkswap = lines[j*num_of_gridpoints_on_line_x+i].edgecross;
            if(bulkswap ==0){
                grid[i+num_of_gridpoints_on_line_x*(j+1)].bulk = grid[i+num_of_gridpoints_on_line_x*j].bulk;
            }
            else if(bulkswap ==1){
                if(grid[i+num_of_gridpoints_on_line_x*j].bulk==0){
                    grid[i+num_of_gridpoints_on_line_x*(j+1)].bulk =1;
                }
                else if(grid[i+num_of_gridpoints_on_line_x*j].bulk==1){
                    grid[i+num_of_gridpoints_on_line_x*(j+1)].bulk=0;
                }
                else{
                    printf("Error: the grid.bulk value should be either 0 or 1.");
                }
            }
            else{
                printf("Error: bulkswap should be either 0 or 1");
            }
            
        }
        
        
    }
    
    
}
/*******************************************************************/
void update_bulk()
{
    double q1array[num_of_gridpoints],q2array[num_of_gridpoints],newq1array[num_of_gridpoints],newq2array[num_of_gridpoints];
    double cosavx,sinavx,cosavy,sinavy,phiavx,phiavy,deltaphix,deltaphiy,deltaphixhelp,deltaphiyhelp,newphi;
    double q1avx,q2avx,q1avy,q2avy;
    double phiavxdouble,phiavxdoublehelp,phiavydouble,phiavydoublehelp;
    double sumangle,sumsin,sumcos,phidouble,phidoublehelp;
    double sumq1,sumq2,newq1,newq2;
    double normalderivative,phiavxhelp,phiavyhelp,targetangle,correctionanglehelp,correctionangle,normalderivative1,normalderivative2,targetangle1,targetangle2,correctionangle1,correctionangle2,correctionanglehelp1,correctionanglehelp2;
    double normalderivative_q1,normalderivative_q2,normalderivative_q1_left,normalderivative_q1_right,normalderivative_q2_left,normalderivative_q2_right,normalderivative_q1_top,normalderivative_q1_bottom,normalderivative_q2_top,normalderivative_q2_bottom; //new variables with Qxx and Qxy
    double correction_q1,correction_q2,correction_q1_left,correction_q1_right,correction_q2_left,correction_q2_right,correction_q1_top,correction_q1_bottom,correction_q2_top,correction_q2_bottom;
    int i,j,k,l,m,n,o,anglecount;
    double sum_cos,sum_sin,gridcount,doubcosav,doubsinav,help_director;
    double sumcossquared,edgecounter,cossquared;
    double testorder;


    
    // Might not need some of these ^
    
    p++;

    
    for(j=0;j<num_of_gridpoints;j++){
        q1array[j]=0;
        q2array[j]=0;
    }
    
    for(i=0;i<num_of_gridpoints;i++){
        if(grid[i].bulk == 1){             // I only want to update the value of q1 an q2 if the gridpoint is actually part of the bulk.  
            
            //In the full program, the edge of the cell might expand outwards. When it does so, the angle at the expanding edge is not defined, but according to the code takes the value "100" that I give to points that are not inside the cell. Therefore I first write a function here that replaces this 100 with a value that is commensurate with the local angles around it: 

            if(grid[i].q1 ==100 && grid[i].q2 ==100){
                //printf("The cell is expanding.");
                anglecount = 0;
                sumq1 = 0.0;
                sumq2 = 0.0;
                sumangle = 0.0;
                
                for(k=0;k<4;k++){
                    if(grid[i+neighbors[k]].q1 != 100 && grid[i+neighbors[k]].q2 != 100){
                        anglecount = anglecount + 1; // Hoeveel buren er eerst wel in e cell zaten (phi != 100)
                        
                        sumq1 = sumq1 + grid[i+neighbors[k]].q1;
                        sumq2 = sumq2 + grid[i+neighbors[k]].q2;
                      
                    }
                }
                if(anglecount == 0){
                    printf("Error: all surrounding points were not in the cell before. Was the update too fast?");
                    break;
                }
                else if(anglecount == 4){
                    printf("Error: all surrounding points are part of the cell. Cell edge is probably not smooth");
                }
                else if(anglecount == 1){
                    newq1array[i] = sumq1;
                    newq2array[i] = sumq2;
                }
                else if(anglecount == 2){
                    newq1array[i] = sumq1*0.5;
                    newq2array[i] = sumq2*0.5;
                }
                else if(anglecount ==3){ //In this case, just use the value of the bulk where the edge is coming from:      
                    newq1array[i] = sumq1/3.0;
                    newq2array[i] = sumq2/3.0;
                   //printf("\n%s\n","Anglecount = 3");
                //An alternative would be to average over the Qxx and Qxy values of the other 3 points.
                    
                    /*if(grid[i+1].q1 ==100){
                        grid[i].q1 = grid[i-1].q1;
                        grid[i].q2 = grid[i-1].q2;
                    }
                    if(grid[i-1].q1 ==100){
                        grid[i].q1 = grid[i+1].q1;
                        grid[i].q2 = grid[i+1].q2;
                    }
                    if(grid[i+num_of_gridpoints_on_line_x].q1 ==100){
                        grid[i].q1 = grid[i-num_of_gridpoints_on_line_x].q1;
                        grid[i].q2 = grid[i-num_of_gridpoints_on_line_x].q2;
                    }
                    if(grid[i-num_of_gridpoints_on_line_x].q1 ==100){
                        grid[i].q1 = grid[i+num_of_gridpoints_on_line_x].q1;
                        grid[i].q2 = grid[i+num_of_gridpoints_on_line_x].q2;
                    }*/
                }
                
                else{
                    printf("Anglecount is out of bound");
                }
                                

                
            }
        }
    }
    
    for(i=0;i<num_of_gridpoints;i++){
        if(grid[i].bulk == 1){             // I only want to update the value of q1 an q2 if the gridpoint is actually part of the bulk.  
            
            //In the full program, the edge of the cell might expand outwards. When it does so, the angle at the expanding edge is not defined, but according to the code takes the value "100" that I give to points that are not inside the cell. Therefore I first write a function here that replaces this 100 with a value that is commensurate with the local angles around it: 

            if(grid[i].q1 ==100 && grid[i].q2 ==100){
                    grid[i].q1 = newq1array[i];
                    grid[i].q2 = newq2array[i];
                    if(grid[i].q1*grid[i].q1+grid[i].q2*grid[i].q2 > 0.275625){
                        printf("The newly included point has S > 1.05");
                        for(k=0;k<4;k++){
                        printf("\n%.10f\t%.10f\t%i\n",grid[i+neighbors[k]].q1,grid[i+neighbors[k]].q2,anglecount);
                }
            }
                
            }
            

        }
    }
    i=0;
    
    for(i=0;i<num_of_gridpoints;i++){
        if(grid[i].bulk == 1){             // I only want to update the value of q1 an q2 if the gridpoint is actually part of the bulk.  
            testorder = 2.0*sqrt(grid[i].q1*grid[i].q1+grid[i].q2*grid[i].q2);
            
            if(testorder > 1.05){
                printf("Error: S > 1.05 for point already in bulk");
                printf("\n%i\t%.10f\t%.10f\t%.10f\n",i,testorder, grid[i].q1, grid[i].q2);
            }                

            //Next, I'm going to update the bulk
            
   
            if(grid[i+1].bulk==1 && grid[i-1].bulk ==1){
                
                q1avx = grid[i+1].q1 + grid[i-1].q1;
                q2avx = grid[i+1].q2 + grid[i-1].q2;
            }
            
            else if(grid[i+1].bulk==1){        
                //Here I now implement a ghost point at i-1, and use an anchoring boundary_condition to determine the value at i-1
                //First I'll define the targetangle, which is stored in the datastructure lines.                
                if(grid[i-1].bulk!=0){
                    printf("error: grid i-1 !=0");
                }
                
                if(lines[num_of_lines+i-1].edgecross ==1){
                    targetangle = lines[num_of_lines+i-1].psi; //We add num_of_lines here because these are horizontal lines.
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 1.\n");
                }
                
                // Need to produce values for q1avx and q2avx.
                
                correction_q1 = cos(2.0*targetangle)*0.5;
                correction_q2 = sin(2.0*targetangle)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                                
                normalderivative_q1 = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1);
                normalderivative_q2 = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2);
                
                q1avx = 2.0*grid[i+1].q1+2.0*lattice_spacing*normalderivative_q1;
                q2avx = 2.0*grid[i+1].q2+2.0*lattice_spacing*normalderivative_q2;
            }
            
            else if(grid[i-1].bulk==1){
                
                if(lines[num_of_lines+i].edgecross ==1){
                    targetangle = lines[num_of_lines+i].psi; //We add num_of_lines here because these are horizontal lines.
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 2.\n");
                }
                
                // Need to produce values for q1avx and q2avx.
                
                correction_q1 = cos(2.0*targetangle)*0.5;
                correction_q2 = sin(2.0*targetangle)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                                
                normalderivative_q1 = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1);
                normalderivative_q2 = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2);
                
                q1avx = 2.0*grid[i-1].q1+2.0*lattice_spacing*normalderivative_q1;
                q2avx = 2.0*grid[i-1].q2+2.0*lattice_spacing*normalderivative_q2;
                
            }
            else{ //Otherwise, the pixel has no neighbors either to the left or to the right. This means a boundary condition should be implied on both:
                
                if(lines[num_of_lines+i-1].edgecross ==1){
                    targetangle1 = lines[num_of_lines+i-1].psi;
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 3.\n");
                }
                
                if(lines[num_of_lines+i].edgecross ==1){
                    targetangle2 = lines[num_of_lines+i].psi;
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 4.\n");                    
                }

                correction_q1_left = cos(2.0*targetangle1)*0.5;
                correction_q2_left = sin(2.0*targetangle1)*0.5;
                
                correction_q1_right = cos(2.0*targetangle2)*0.5;
                correction_q2_right = sin(2.0*targetangle2)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                
                normalderivative_q1_left = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1_left);
                normalderivative_q2_left = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2_left);
                
                normalderivative_q1_right = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1_right);
                normalderivative_q2_right = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2_right);
                
                q1avx = 2.0*grid[i].q1+lattice_spacing*(normalderivative_q1_left+normalderivative_q1_right);
                q2avx = 2.0*grid[i].q2+lattice_spacing*(normalderivative_q2_left+normalderivative_q2_right);
                
            }
            //Now do the same for the y-direction. 
            
            if(grid[i+num_of_gridpoints_on_line_x].bulk==1 && grid[i-num_of_gridpoints_on_line_x].bulk ==1){
                q1avy = grid[i+num_of_gridpoints_on_line_x].q1 + grid[i-num_of_gridpoints_on_line_x].q1;
                q2avy = grid[i+num_of_gridpoints_on_line_x].q2 + grid[i-num_of_gridpoints_on_line_x].q2;
            }
            
            //Note that in the current code if both x and y directions have a ghost point then the boundary condition is implemented twice. I need to correct for that later.
                      
            else if(grid[i+num_of_gridpoints_on_line_x].bulk==1){      //Here I now implement a ghost point at j-1, and use an anchoring boundary_condition to determine the value at j-1.
                //First I'll define the targetangle, which is stored in the datastructure lines.                
                if(lines[i-num_of_gridpoints_on_line_x].edgecross ==1){
                    targetangle = lines[i-num_of_gridpoints_on_line_x].psi; //We do not add num_of_lines here because these are vertical lines.
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 5.\n");
                }
                
                // Need to produce values for q1avx and q2avx.
                
                correction_q1 = cos(2.0*targetangle)*0.5;
                correction_q2 = sin(2.0*targetangle)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                                
                normalderivative_q1 = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1);
                normalderivative_q2 = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2);
                
                q1avy = 2.0*grid[i+num_of_gridpoints_on_line_x].q1+2.0*lattice_spacing*normalderivative_q1;
                q2avy = 2.0*grid[i+num_of_gridpoints_on_line_x].q2+2.0*lattice_spacing*normalderivative_q2;

            }
            else if(grid[i-num_of_gridpoints_on_line_x].bulk==1){
                if(lines[i].edgecross ==1){
                    targetangle = lines[i].psi; //We do not add num_of_lines here because these are vertical lines.
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 6.\n");
                }
                
                // Need to produce values for q1avx and q2avx.
                
                correction_q1 = cos(2.0*targetangle)*0.5;
                correction_q2 = sin(2.0*targetangle)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                                
                normalderivative_q1 = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1);
                normalderivative_q2 = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2);
                
                q1avy = 2.0*grid[i-num_of_gridpoints_on_line_x].q1+2.0*lattice_spacing*normalderivative_q1;
                q2avy = 2.0*grid[i-num_of_gridpoints_on_line_x].q2+2.0*lattice_spacing*normalderivative_q2;
                
            }
            else{ //Otherwise, the pixel has no neighbors either to the top or to the bottom. This means a boundary condition should be implied on both:
                
                if(lines[i-num_of_gridpoints_on_line_x].edgecross ==1){
                    targetangle1 = lines[i-num_of_gridpoints_on_line_x].psi;
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 7.\n");
                }
                
                if(lines[i].edgecross ==1){
                    targetangle2 = lines[i].psi;
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Message 8.\n");
                    
                }

                correction_q1_top = cos(2.0*targetangle1)*0.5;
                correction_q2_top = sin(2.0*targetangle1)*0.5;
                
                correction_q1_bottom = cos(2.0*targetangle2)*0.5;
                correction_q2_bottom = sin(2.0*targetangle2)*0.5;
                
                // correction_q1 = Qxx0, correction_q2 = Qxy0
                // Qxx0 = (S/2)*cos(2*targetangle) with S=1
                
                normalderivative_q1_top = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1_top);
                normalderivative_q2_bottom = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2_top);
                
                normalderivative_q1_top = -2.0*anchoring/elasticity*(grid[i].q1-correction_q1_bottom);
                normalderivative_q2_bottom = -2.0*anchoring/elasticity*(grid[i].q2-correction_q2_bottom);
                
                q1avy = 2.0*grid[i].q1+lattice_spacing*(normalderivative_q1_top+normalderivative_q1_bottom);
                q2avy = 2.0*grid[i].q2+lattice_spacing*(normalderivative_q2_top+normalderivative_q2_bottom);
                
                
                //phiavx = 300;
            
            }      
            // No noise added to Qxx and Qxy calcs
            

                //printf("q1avx newq1calc:");
                //printf("\n%i\t%.10f\t%0.10f\n",i,grid[i].q2,q1avx);
            /*if(p == 20041 && i == 334){
                printf("q1avx = \t%.10f",q1avx);
                printf("q1avy = \t%.10f",q1avy);
                printf("q2avx = \t%.10f",q2avx);
                printf("q2avy = \t%.10f",q2avy);
            }*/
                
            newq1 = grid[i].q1+timestepbulk*elasticity/rotational_friction*((q1avx+q1avy-4.0*grid[i].q1)/(lattice_spacing*lattice_spacing)-2.0*grid[i].q1/(defect_core*defect_core)*(4.0*(grid[i].q1*grid[i].q1+grid[i].q2*grid[i].q2)-phase_field));
                //printf("newq1:");
                //printf("\n%i\t%.10f\t%0.10f\n",i,grid[i].q2,newq1);
                //franck_constant = elasticity.
            newq2 = grid[i].q2+timestepbulk*elasticity/rotational_friction*((q2avx+q2avy-4.0*grid[i].q2)/(lattice_spacing*lattice_spacing)-2.0*grid[i].q2/(defect_core*defect_core)*(4.0*(grid[i].q1*grid[i].q1+grid[i].q2*grid[i].q2)-phase_field));

            q1array[i] = newq1;
            q2array[i] = newq2;
            
            // NEEDS UPDATE 
            /*
            else if(phiavx != 300 && phiavy == 300){
                newphi = grid[i].phi + timestepbulk*elasticity/(rotational_friction*lattice_spacing*lattice_spacing)*(2*deltaphix)+noise*(2*ran2(&seed)-1);
            }
            else if(phiavx == 300 && phiavy != 300){
                newphi = grid[i].phi + timestepbulk*elasticity/(rotational_friction*lattice_spacing*lattice_spacing)*(2*deltaphiy)+noise*(2*ran2(&seed)-1);
            }
            else{
                printf("Error: this pixel has no neighbors in any direction.");
            }
            if (newphi >pi/2){
                phiarray[i] = newphi-pi;
            }
            else if (newphi <=-pi/2){
                phiarray[i] = newphi + pi;
            }
            else{
                phiarray[i] = newphi;
            }*/

        }
        else{
            q1array[i] = 100;   //To make sure I don't get confused by the gridvalues outside the bulk, I give these gridpoints a value outside of the domain -pi,pi.
            q2array[i] = 100;
        }
    }
    
    for(l=0;l<num_of_gridpoints;l++){
        grid[l].q1 = q1array[l];
        grid[l].q2 = q2array[l];
        if(grid[l].bulk ==1 && (grid[l].q1 > 0.6 || grid[l].q1 < -0.6 || grid[l].q2 > 0.6 || grid[l].q2 < -0.6)){ //change
            printf("\n Error: the angle is out of bound. Qxx or Qxy > 0.6\n");
            printf("\n Caclulation nr %i\n",p);
        }
    }
    
    for(o=0;o<num_of_gridpoints;o++){ // Calculates the values for phi and S from q1 and q2.
        grid[o].loworder = 0;
        if(grid[o].q1 != 100 && grid[o].q2 != 100){
            grid[o].s = 2.0*sqrt(grid[o].q1*grid[o].q1+grid[o].q2*grid[o].q2);
            if(grid[o].s < lowordercutoff){
                grid[o].loworder = 1;
            }
            
            if(grid[o].s == 0){
                printf("\nS = 0\n");
            }
            else if(grid[o].s > 1.05){
                printf("\nS > 1.05");
                printf("\n%i\t%.10f\t%.10f\t%.10f\t%i\n",o,grid[o].s,grid[o].q1,grid[o].q2,p);
            }
                
            if(grid[o].q1 > 0){
                grid[o].phi = 0.5*atan(grid[o].q2/grid[o].q1);
            }
            else if(grid[o].q1 < 0){
                if(grid[o].q2 > 0){
                    grid[o].phi = 0.5*(atan(grid[o].q2/grid[o].q1)+pi);
                }
                else if(grid[o].q2 < 0){
                    grid[o].phi = 0.5*(atan(grid[o].q2/grid[o].q1)-pi);
                }
                else{
                    grid[o].phi = 0.5*pi;
                }
            }
            else{
                if(grid[o].q2 > 0){
                    grid[o].phi = 0.25*pi;
                }
                else{
                    grid[o].phi = 0.75*pi;
                }
            }    
        }
        else{
            grid[o].s = 100;
            grid[o].phi = 100;
        }
    }
            
        
        
    //Now the new configuration is entirely defined. Now I want to calculate the (global) order parameter corresponding to this configuration, and the average director.
    
    

    sum_cos = 0;
    sum_sin = 0;
    gridcount = 0;
    for(m = 0; m<num_of_gridpoints; m++){
        if(grid[m].bulk ==1){
            sum_cos = sum_cos + cos(2*grid[m].phi);
            sum_sin = sum_sin + sin(2*grid[m].phi);
            gridcount = gridcount + 1.0;
        }
    }
    doubcosav = sum_cos/gridcount;
    doubsinav = sum_sin/gridcount;
    order_par_global = sqrt(doubcosav*doubcosav+doubsinav*doubsinav);
    
    //Now the director:
    if(sum_cos ==0){
        if(sum_sin>0){
            average_director = 0.25*pi;
        }
        else{
            average_director = -0.25*pi;
        }
    }
    else{
        help_director = atan(sum_sin/sum_cos);
        if(sum_cos<0){
            if(sum_sin>=0){
                help_director = help_director + pi;
            }
            else{
                help_director = help_director - pi;
            }
        }
        average_director = 0.5*help_director;
    }  
    
    //I also want to calculate the wall order parameter, which measures to what degree the edge of the bulk alligns to the targetangle.
    
    sumcossquared = 0.0;
    edgecounter = 0.0;
    for(n= 0; n<num_of_gridpoints; n++){
        if(grid[n].bulk ==1){
            //Note, in the following I check in four directions whether there is a neighbour that is not inside the cell. I ignore the possibility of having two neighbours with different targetangles. In practice this should only be a problem in the corners of the cell.
            if(grid[n-1].bulk ==0){
                if(lines[num_of_lines+n-1].edgecross ==1){
                    targetangle = lines[num_of_lines+n-1].psi; //We add num_of_lines here because these are horizontal lines.
                    if (targetangle > 0.5*pi || targetangle < -0.5*pi){
                        printf("Error: something is wrong with calculating targetangle.");
                    }
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Calculating wall parameter, message 1.\n");
                }
                cossquared = cos(grid[n].phi-targetangle)*cos(grid[n].phi-targetangle);
                sumcossquared = sumcossquared + cossquared;
                edgecounter = edgecounter + 1.0;
            }
            else if(grid[n+1].bulk ==0){
                if(lines[num_of_lines+n].edgecross ==1){
                    targetangle = lines[num_of_lines+n].psi; //We add num_of_lines here because these are horizontal lines.
                    if (targetangle > 0.5*pi || targetangle < -0.5*pi){
                        printf("Error: something is wrong with calculating targetangle.");
                    }
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Calculating wall parameter, message 2.\n");
                }
                cossquared = cos(grid[n].phi-targetangle)*cos(grid[n].phi-targetangle);
                sumcossquared = sumcossquared + cossquared;
                edgecounter = edgecounter + 1.0;
            }
            else if(grid[n-num_of_gridpoints_on_line_x].bulk ==0){
                if(lines[n-num_of_gridpoints_on_line_x].edgecross ==1){
                    targetangle = lines[n-num_of_gridpoints_on_line_x].psi; //We add num_of_lines here because these are horizontal lines.
                    if (targetangle > 0.5*pi || targetangle < -0.5*pi){
                        printf("Error: something is wrong with calculating targetangle.");
                    }
                }
                else{
                    printf("Error: there is a mistake in the linesarray. Calculating wall parameter, message 3.\n");
                }
                cossquared = cos(grid[n].phi-targetangle)*cos(grid[n].phi-targetangle);
                sumcossquared = sumcossquared + cossquared;
                edgecounter = edgecounter + 1.0;
            }
            else if(grid[n+num_of_gridpoints_on_line_x].bulk ==0){
                if(lines[n].edgecross ==1){
                    targetangle = lines[n].psi; //We add num_of_lines here because these are horizontal lines.
                    if (targetangle > 0.5*pi || targetangle < -0.5*pi){
                        printf("\nError: something is wrong with calculating targetangle.");
                    }
                }
                else{
                    printf("\nError: there is a mistake in the linesarray. Calculating wall parameter, message 4.\n");
                }
                cossquared = cos(grid[n].phi-targetangle)*cos(grid[n].phi-targetangle);
                sumcossquared = sumcossquared + cossquared;
                edgecounter = edgecounter + 1.0;
                //printf("%.10f\t%.10f\t%.10f\n",grid[n].phi,targetangle,cossquared);
            }
        }
    }
    wall_order_par = 2.0*sumcossquared/edgecounter - 1.0;
    //printf("%.10f\t%.10f\n",sumcossquared,edgecounter);

}

/*******************************************************************/

void export_conf()
{
	static long local_time=0;
        char f_na[32];
        char f_na2[32];
        char f_na3[32];
        char f_na4[32];
        long i;
        
    
    
    
    
    // To do: Split boundary and bulk points in new file
    
    
    
    
	FILE *f_ou;	
        sprintf(f_na,"bulk_%ld.dat",local_time);
        f_ou = fopen(f_na,"w");

        /*for (i = 0; i<num_of_bounpoint; i++){
                fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\n",vertex[i].x,vertex[i].y,vertex[i].lam,vertex[i].kappa);
        } */  
        
        /*for (i = 0; i<num_of_gridpoints; i++){
                fprintf(f_ou,"%.10f\t%.10f\t%.10f\n",grid[i].x,grid[i].y,grid[i].phi);
        }*/
        
        /*for (i=0;i<2*num_of_lines;i++){
            fprintf(f_ou,"%.10f\t%.10f\t%i\t%.10f\n",lines[i].x,lines[i].y,lines[i].edgecross,lines[i].psi);
        }*/
        
        //This I normally want to store

        for (i=0;i<num_of_gridpoints;i++){
            fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%i\t%.10f\t%.10f\t%.10f\t%ld\n",grid[i].x,grid[i].y,grid[i].phi,grid[i].bulk,grid[i].q1,grid[i].q2,grid[i].s,grid[i].loworder);
        }
        
	fclose(f_ou);
    
    FILE *f_ou2;	
        sprintf(f_na2,"vertex_%ld.dat",local_time);
        f_ou2 = fopen(f_na2,"w");

        /*for (i = 0; i<num_of_bounpoint; i++){
                fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\n",vertex[i].x,vertex[i].y,vertex[i].lam,vertex[i].kappa);
        } */  
        
        /*for (i = 0; i<num_of_gridpoints; i++){
                fprintf(f_ou,"%.10f\t%.10f\t%.10f\n",grid[i].x,grid[i].y,grid[i].phi);
        }*/
        
        /*for (i=0;i<2*num_of_lines;i++){
            fprintf(f_ou,"%.10f\t%.10f\t%i\t%.10f\n",lines[i].x,lines[i].y,lines[i].edgecross,lines[i].psi);
        }*/
        
        //This I normally want to store
        for (i=0;i<num_of_bounpoint;i++){
            fprintf(f_ou2,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",vertex[i].x,vertex[i].y,vertex[i].lam,vertex[i].kappa,vertex[i].theta);
        }
        
	fclose(f_ou2);
        
    FILE *f_ou3;
        
        sprintf(f_na3,"force_%ld.dat",local_time);
        f_ou3 = fopen(f_na3,"w");
        
        for(i=0;i<2*num_of_adhesions;i++){
            fprintf(f_ou3,"%.10f\t%.10f\t%i\t%.10f\t%.10f\n",forces[i].x,forces[i].y,forces[i].type,forces[i].fx,forces[i].fy);
        }
        
    fclose(f_ou3);
    
    FILE *f_ou4;	
        sprintf(f_na4,"global_%ld.dat",local_time);
        f_ou4 = fopen(f_na4,"w");


        fprintf(f_ou4,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%ld\t%ld\n",alpha, sigma,order_par_global,wall_order_par,average_director,total_force_x,total_force_y,total_force,q_difference_total,shape_difference_total);
        
	fclose(f_ou4);
    
    
        
        local_time++;
}

void export_conf_exp()
{
	static long local_time=0;
        char f_na5[32];
        char f_na6[32];
        char f_na7[32];
        char f_na8[32];

        long i,j,n;

	FILE *f_ou5;	
        sprintf(f_na5,"exp_bulk_shape_%ld_linepoints_%ld.dat",shape_number, num_of_gridpoints_on_line_in_box);
        f_ou5 = fopen(f_na5,"w");


        for (i=0;i<num_of_gridpoints;i++){
            fprintf(f_ou5,"%.10f\t%.10f\t%.10f\t%i\t%.10f\t%.10f\t%.10f\n",grid[i].x,grid[i].y,grid_exp[i].phi,grid_exp[i].bulk,grid_exp[i].q1,grid_exp[i].q2,grid_exp[i].s); // Still need to change all values to grid_exp and not grid.
        }
        
	fclose(f_ou5);
    
    FILE *f_ou6;	
        sprintf(f_na6,"exp_border_shape_%ld_linepoints_%ld.dat",shape_number, num_of_gridpoints_on_line_in_box);
        f_ou6 = fopen(f_na6,"w");


        for (j=0;j<num_of_borderpoints;j++){
            fprintf(f_ou6,"%.10f\t%.10f\n",border_exp[j].x,border_exp[j].y); // Still need to change all values to grid_exp and not grid.
        }
        
	fclose(f_ou6);
    
    FILE *f_ou7;	
        sprintf(f_na7,"exp_comparison_%ld_linepoints_%ld.dat",shape_number, num_of_gridpoints_on_line_in_box);
        f_ou7 = fopen(f_na7,"w");

        fprintf(f_ou7,"%.10f\t%ld\t%ld\t%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",anchoring,bulkpoints_sim,bulkpoints_exp,bulkpoints_shared,deltasquaredq,deltasquaredangle,deltasquaredanglecom,deltasquaredanglecom2, errordeltasquared);
        
	fclose(f_ou7);
    
    
    FILE *f_ou8;	
        sprintf(f_na8,"coherence_shape_%ld.dat",shape_number);
        f_ou8 = fopen(f_na8,"w");
        for(n=0; n<num_of_gridpoints_exp; n++){
            if(n%num_of_gridpoints_on_line_exp==0 && n!=0){
                fprintf(f_ou8,"\n");
            }
        fprintf(f_ou8,"%.3f\t",grid_exp[n].coherency); // Still need to change all values to grid_exp and not grid.
        

        }
                fclose(f_ou8);
        //local_time++;
}

/*******************************************************************/

void equilibrium()
{
    
    int i,j;
    int shape_difference, q_difference;
    double q1dif,q2dif,q_cutoff;

    q_cutoff = 0.1;

    if(iteration==0){
        for(i=0;i<num_of_gridpoints;i++){
            grid_old[i].bulk=grid[i].bulk;
            grid_old[i].q1=grid[i].q1;
            grid_old[i].q2=grid[i].q2;
        }
    }

    else{
        shape_difference=0;
        q_difference=0;
        for(j=0;j<num_of_gridpoints;j++){
            if(grid_old[j].bulk!=grid[j].bulk){
                shape_difference++;                
            }
            if((grid_old[j].bulk==1)&&(grid[j].bulk==1)){
                q1dif=sqrt((grid_old[j].q1-grid[j].q1)*(grid_old[j].q1-grid[j].q1));
                q2dif=sqrt((grid_old[j].q2-grid[j].q2)*(grid_old[j].q2-grid[j].q2));
                if((q1dif>q_cutoff)||(q2dif>q_cutoff)){
                    q_difference++;
                }
            }
            grid_old[j].bulk=grid[j].bulk;
            grid_old[j].q1=grid[j].q1;
            grid_old[j].q2=grid[j].q2;
        }
        
        if(shape_difference>0){
            printf("\nDifference in the shape: %ld gridpoints have changed in being part of bulk after %ld iterations.\n", shape_difference, iteration);
            shape_difference_total=iteration/period; // This is the highest period where this is the case.
        }
        if(q_difference>0){
            printf("\nDifference in the Q: %ld gridpoints have changed in Q more than %.10f after %ld iterations.\n", q_difference, q_cutoff, iteration);
            q_difference_total=iteration/period; // This is the highest period where this is the case.
        }
    }
    
    
}


/*******************************************************************/

void end()
{
	long i;
	
        CPU_Time cpu_time;
        
	FILE *f_ou;
	
	f_ou = fopen("ou.dat","w");
	for (i=0; i<num_of_bounpoint; i++){
		fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\n",vertex[i].x,vertex[i].y,vertex[i].lam,vertex[i].kappa);
	}
	/*for (i = 0; i<num_of_triangles; i++){
                fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",vertex[triangle[i].vid[0]].x,vertex[triangle[i].vid[0]].y,vertex[triangle[i].vid[1]].x,vertex[triangle[i].vid[1]].y,vertex[triangle[i].vid[2]].x,vertex[triangle[i].vid[2]].y,counter[i],neighbourarray[i][0],neighbourarray[i][1],neighbourarray[i][2],neighbourarray[i][3],neighbourarray[i][4],neighbourarray[i][5]);
        }*/
	fclose(f_ou);
        
	get_time(&cpu_time);
        printf("\n");
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
		cpu_time.d,cpu_time.h,
	 	cpu_time.m,cpu_time.s);
  	printf("\n");
}

/*******************************************************************/
/*void progress_bar(long i)
{
	long j, ticks, percent, num_of_ticks=20;
	char progress[32];

	percent = 100*(i+1)/num_of_iteration;
	ticks = num_of_ticks*percent/100;
	
	for (j=0; j<20; j++){
		progress[j] = ((j<=ticks)?'#':' ');
	}
	progress[20] = '\0';
	
	printf("Progress: [%s] %ld%%\r",progress,percent);
	fflush(stdout);
}*/

/*******************************************************************/

void get_time(CPU_Time *cpu_time)
{
  	long sec_in_day, sec_in_ora, sec_in_min;

	time_t t = get_cpu_time();

  	sec_in_day = 86400;
  	sec_in_ora = 3600;
  	sec_in_min = 60;

  	cpu_time->d = t/sec_in_day; t = t%sec_in_day;
  	cpu_time->h = t/sec_in_ora; t = t%sec_in_ora;
  	cpu_time->m = t/sec_in_min; t = t%sec_in_min;
  	cpu_time->s = t;
}

/*******************************************************************/

time_t get_cpu_time()
{
	static int first_call=1;
	static time_t t_ini;

	time_t elap, t_now;
		
	if (first_call){
		time(&t_ini);
		first_call=0;
		elap=0;
	}
	else {
		time(&t_now);
		elap=t_now-t_ini;
	}	

	return elap;
}
/*******************************************************************/

double ran2(long *idum)
{
   	int j;
  	long k;
  	static long idum2=123456789;
  	static long iy=0;
  	static long iv[NTAB];
  	double temp;
 
  	if (*idum <= 0) { 
    	if (-(*idum) < 1) *idum=1; 
    	else *idum = -(*idum);
    	idum2 = (*idum);
    	for (j=NTAB+7; j>=0; j--) {
      		k = (*idum)/IQ1;
      		*idum = IA1*(*idum-k*IQ1)-k*IR1;
      		if (*idum < 0) *idum += IM1;
      		if (j < NTAB) iv[j] = *idum;
    	}
    	iy = iv[0];
  	}
  	k = (*idum)/IQ1; 
  	*idum = IA1*(*idum-k*IQ1)-k*IR1; 
  	if (*idum < 0) *idum += IM1; 
  	k = idum2/IQ2;
  	idum2 = IA2*(idum2-k*IQ2)-k*IR2; 
  	if (idum2 < 0) idum2 += IM2;
  	j = iy/NDIV; 
  	iy = iv[j]-idum2;
  	iv[j] = *idum; 
  	if (iy < 1) iy += IMM1;
  	if ((temp = AM*iy) > RNMX) return RNMX; 
  	else return temp;
}

