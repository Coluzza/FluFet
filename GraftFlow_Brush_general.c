/*This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "./lib/param_parser.h"

#define DBL_EPSILON 2.2204460492503131E-16
#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE 1e-10



#define N1 (1.0/Protein->N)
#define S 21
#define ATOM_N 0
#define ATOM_CA 1
#define ATOM_C 2
#define ATOM_O 3
#define ATOM_H 4

#define PROTEIN 			 Param->Nres
#define COLLOID 			 PROTEIN+1   
#define LINKER  			 PROTEIN+2 
#define ANCHOR_LINKER  PROTEIN+3 
#define POLYMER 			 PROTEIN+4 
#define ANCHOR_POLYMER PROTEIN+5 
#define SRD 					 PROTEIN+6 
#define END_POLYMER 	 PROTEIN+7
 
#define END_NEUTRAL 0
#define END_NEUTRAL_ACCEPTOR 1
#define END_NEUTRAL_DONOR 2
#define END_NEUTRAL_ACCEPTOR_DONOR 3
#define END_POSITIVE_ACCEPTOR 4
#define END_POSITIVE_DONOR 5
#define END_POSITIVE_ACCEPTOR_DONOR 6
#define END_NEGATIVE_ACCEPTOR 7
#define END_NEGATIVE_DONOR 8
#define END_NEGATIVE_ACCEPTOR_DONOR 9
#define POL_END_TYPES 10


#define ALL 1
#define PROTEIN_ONLY 0
#define BRUSH_ONLY 2
#define NO_POL_ANCHORS 3

#define LAMMPS_SCALE 1
#define XTC_SCALE 10


#define ATOM_IGNORE 3
#define ANCHOR_TYPE 5
#define MOL_ATOM_TYPES 7
#define NATOM 5
#define MOL_DIST 10

#define YES 1
#define NO 0

#define CYLINDER 0
#define SPHERE 1
#define SLIT 2
#define CHANNEL 3

#define TRANS 0
#define ANCHOR 1



#define TRUE 1
#define FALSE 0


double dround(double x) { return floor(x+0.5); }
double SQR(double x) { return x*x; }
double sqrlen(double v[3]) {
	double d2 = 0.0;
	int i;
	for(i=0;i<3;i++)
	d2 += SQR(v[i]);
	return d2;
}
void get_mi_vector(double res[3], double a[3], double b[3])
	{
	int i;
	
	for(i=0;i<3;i++) {
		res[i] = a[i] - b[i];
		#ifdef PARTIAL_PERIODIC
		if (PERIODIC(i))
		#endif
		// res[i] -= dround(res[i]*box_l_i[i])*box_l[i];
	}
}

void vector_product(double a[3], double b[3], double c[3]) {
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
	return ;
}
double scalar(double a[3], double b[3]) {
	double d2 = 0.0;
	int i;
	for(i=0;i<3;i++)
	d2 += a[i]*b[i];
	return d2;
}




struct part{
	double x,y,z;
	int id;
	int flag;
	int res;
	int lammps_idx;
	int prot;
};

/*******SHEA Parameters*****/




struct GlobalSistem{
	double Bbond[6],DotProd[3];
	double ECA_Range;
	double ECA_Range2;
	double Half_Box_x;
	double Half_Box_y;
	double Half_Box_z;
	double box_x;
	double box_y;
	double box_z;
	int N_Atoms, N_Angles, N_Dihedrals,N_Bonds;
	int N_Mol;
	int Tot_Atoms, Tot_Angles, Tot_Dihedrals,Tot_Bonds;
	int  Tot_Angles_Types, Tot_Dihedrals_Types;
	int tot_atom_type;
	int tot_bonds_type;
	int polymers_bonds_type;
	double radius,inradius;
	double chi_p[21];	
	int bondtype;
	double MOL_SHIFT;
};

struct Proteins{
	int N;
	int Nres;
	int NMasked;
	int NSkip;
	int Nread;
	struct part *Atom;
	struct part *Atom_read;
	struct part *Linker;
	struct part *Linker_Anchors;
	char Atoms_name[NATOM][10];
	char Amminoacids[21][4];
	double max_dist;
	double *Average_r;
	double *Count_r2;
	struct part *gmass;
	int *Ignore;
	int *Mask;
	int *Mask_wall;
	double FIRST_AtomHx,FIRST_AtomHy,FIRST_AtomHz;
	int bondtype_head;
	char *group_id;
};

struct Brushes{
	struct part *Anchors;
	struct part *Atom;
	double maxx;
	double minx;
	double N_Chain_Frac[POL_END_TYPES];
	int N_Anchor;
	int *indeces;
};

struct Interaction_Param{
	char *seq;
	char *mask;
	char *mask_wall_file;
	char *mask_prot;
	char *ignore;
	char *pdb_filename;
	char *prefix;	
	double Linker_Size;
	int Linker_Length;
	double Polymer_Size;
	int Polymer_Length;
	int N_Polymer;
	int N_Proteins;
	double Polymer_Fraction;
	double Box_Scale_x;
	double Box_Scale_yz;
	double chi_s;
	long seed;
	unsigned int useed;
	double Wall_Attr;
	double Wall_Attr_Gen;
	double E_Scale;
	double init_temp;
	double final_temp;
	double bin_flow_temp;
	int steps_flow_temp;
	double init_flow_temp;
	double final_flow_temp;
	int random_rotation_flag;
	double pressure;
	int mask_flag;
	int protein_flag;
	int mask_wall_flag;
	int brush_flag;
	int Geometry;
	char* Geometry_type;
	int Simul_type;
	char *Simul_string;
	int Flow_Flag;
	double FRACTIONS[POL_END_TYPES];
	int N_Colloids;
	int Nres;
	double Coll_R;
	double F_Neutral;
	double F_Neutral_Acceptor;
	double F_Neutral_Donor;
	double F_Neutral_Acceptor_Donor;
	double F_Positive_Acceptor;
	double F_Positive_Donor;
	double F_Positive_Acceptor_Donor;
	double F_Negative_Acceptor;
	double F_Negative_Donor;
	double F_Negative_Acceptor_Donor;
	int GPU;
	char GPUFF[5];
};





double PI,PI_2,PI2,PI_180;
unsigned int Nparam;



void sort2(unsigned long n, double arr[], double brr[]);
int partint 		(double e);

void initialize_parameters (struct Interaction_Param *Param,struct param **ptr_parameters, unsigned int *n_param);

double ran3	(long *idum);



/*Functions to write the atomparam file for lammps*/

void Init_Topology (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins *Protein, struct Brushes *Brush, int type);
void write_topology (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush,FILE *f_atoms,double scale, int type);
void write_masses (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins *Protein);

/*Functions to write the protein section of  atomparam file for lammps*/
void create_protein_copies (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers);
void write_protein_atoms (struct GlobalSistem *sist, struct Interaction_Param *Param, FILE *f_atoms, struct Proteins **Protein_pointers, double scale);
void write_linker_atoms (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, double scale);
void write_protein_dihedrals (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers);
void write_protein_angles (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers);
void write_linker_bonds (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers);
void write_protein_bonds (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param,  struct Proteins **Protein_pointers);
void proteins_equidistant(struct Interaction_Param *Param, struct part *Sphere_Points);
		
/*Functions to write the brush section of  atomparam file for lammps*/
void write_brush_atoms 		(struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type);
void write_brush_bonds 		(struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Brushes *Brush,int type);
void grow_cylinder_brush	(struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type);
void grow_plane_brush 		(struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type);
void grow_sphere_brush 		(struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type);

void write_binary (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein,FILE *fp);
void write_pdb (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein,FILE *fp);
void write_sequence (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp);




/*Functions to read the protein files*/
void allocatePDB (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void deep_copy_proteins(struct Proteins *src,struct Proteins *dest);
void copy_part (struct part *src,struct part *dest, int N);
void create_Protein (struct GlobalSistem *sist, struct Interaction_Param *Param,struct Proteins *Protein);
void readProtCoordinates (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void readIgnoreBonds (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void readMasks (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void test_protein (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void reorderProtCoordinates (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);

/*Functions to manipulate the protein coordiantes*/
void protein_properties (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void protein_gmass (struct GlobalSistem *sist, struct Proteins *Protein);
void protein_rotate_shift (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void Angle_generator (struct GlobalSistem *sist, struct Proteins *Protein);
void Rescaler(struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void Ref (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein);
void Sheadecoder (char *Dec, int N);
void fastadecoder (int *Dec, char *Enc, int N);
void Rotation (double *,double *,double *,double ,double ,double ,double ,int );
double calc_dihedralf_angle_LAMMPS(double *p1, double *p2, double *p3, double *p4);

/*Functions to create the brush coordinates*/
void create_brush_anchor_cylinder 	(struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush);
void create_brush_anchor_plane 		(struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush);
void create_brush_anchor_sphere   	(struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush);
void RandomAnchorsType(struct Interaction_Param *Param, int N, int Nselect, struct Brushes *Brush);
/*Functions to write the in.file for lammps*/
void write_infile (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush, FILE *fp, int prev_step, int cur_step);


void write_variables (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp, struct Proteins *Protein);
void write_groups (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush,FILE *fp);
void write_pair_coeff (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush, FILE *fp);
void write_protein_coeff(struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp);
void write_angles (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp);
void write_fix (struct Interaction_Param *Param, FILE *fp);
void write_fix_prot (struct Interaction_Param *Param, FILE *fp);
void write_harmonicbonds (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp);

void write_1_minimize (FILE *fp, struct Interaction_Param *Param);
void write_2_minimize (FILE *fp, struct Interaction_Param *Param);
void write_3_minimize (FILE *fp, struct Interaction_Param *Param);
void write_4_docking (struct GlobalSistem *sist, FILE *fp, struct Interaction_Param *Param,struct Proteins *Protein);
void write_5_anchored (FILE *fp,struct Interaction_Param *Param);
void write_fixbrooks (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers, FILE *fp);
void write_6_bound_init_T (FILE *fp,struct Interaction_Param *Param);
void write_6_bound_final_T (FILE *fp,struct Interaction_Param *Param);
void write_7_srd (FILE *fp,struct Interaction_Param *Param, int restart_flag, double flow_temp,int srd_step);

/*Functions to write the tcl scripts for VMD*/

void write_tcl_header (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp,double scale);
void write_tcl_protein (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp,double scale);
void write_tcl_all (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp,double scale);
void write_tcl_flow (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp,double scale);
int main(int argc, char** argv){
	
	int i,j,k;
	char line[1024];
	FILE *fp=NULL;
	
	struct Interaction_Param *Param;
	struct GlobalSistem *sist=NULL;
	struct Proteins *Protein=NULL;
	struct Proteins **Protein_pointers=NULL;
	struct Brushes *Brush=NULL;
	struct param *parameters=NULL;
	
	
	sist= (struct GlobalSistem *)calloc(1,sizeof(struct GlobalSistem));
	Param= (struct Interaction_Param *)calloc(1,sizeof(struct Interaction_Param));
	
	
	Brush=(struct Brushes*)calloc(1,sizeof(struct Brushes));
	
	PI=atan(1)*4;
	PI_2=atan(1)*2;
	PI2=atan(1)*8;
	PI_180=atan(1)*4/180;
	
	sist->ECA_Range=12.0;
	sist->ECA_Range2=12*12;
	
	
	sist->chi_p[1]=1.8;
	sist->chi_p[2]=2.5;
	sist->chi_p[3]=-3.5;
	sist->chi_p[4]=-3.5;
	sist->chi_p[5]=2.8;
	sist->chi_p[6]=-0.4;
	sist->chi_p[7]=-3.2;
	sist->chi_p[8]=4.5;
	sist->chi_p[9]=-3.9;
	sist->chi_p[10]=3.8;
	sist->chi_p[11]=1.9;
	sist->chi_p[12]=-3.5;
	sist->chi_p[13]=1.6;
	sist->chi_p[14]=-3.5;
	sist->chi_p[15]=-4.5;
	sist->chi_p[16]=-0.8;
	sist->chi_p[17]=-0.7;
	sist->chi_p[18]=4.2;
	sist->chi_p[19]=-0.9;
	sist->chi_p[20]=-1.3;
	
	
	sist->Bbond[0]=1.4500000;//sist->RbondCaN= 1.45000;
	sist->Bbond[1]=1.5200000;//sist->RbondCCa=1.52000;
	sist->Bbond[2]=1.2300000;//sist->RbondCO= 1.23000;
	sist->Bbond[3]=1.3300000;//sist->RbondCN= 1.33000;
	sist->Bbond[4]=1.0000000;//sist->RbondNH= 1.00000;
	
	
	
	
	printf("Constants Defined\n");fflush(NULL);
	
	Param->mask=NULL;
	Param->pdb_filename=NULL;
	Param->prefix=NULL;
	
	initialize_parameters(Param,&parameters,&Nparam);
	printf("Param Init Nparam=%d\n",Nparam);fflush(NULL);
	get_parameters("param.dat",parameters,Nparam);
	printf("Masks \n");
	printf("Protein->Mask  -%s-\n",Param->mask);fflush(NULL);
	printf("Protein_Mask_Wall  -%s-\n",Param->mask_wall_file);fflush(NULL);
	printf("Protein_Mask_Protonation  -%s-\n",Param->mask_prot);fflush(NULL);
	printf("Ignore_Bond_List  -%s-\n",Param->ignore);fflush(NULL);
	
	
	printf("Simul type \n");
	printf("Protein->Geometry_type  -%s-\n",Param->Geometry_type);fflush(NULL);
	printf("Simul_string  -%s-\n",Param->Simul_string);fflush(NULL);
	printf("Protein_Mask_Protonation  -%s-\n",Param->mask_prot);fflush(NULL);

	
	if (strcmp(Param->Geometry_type, "CYLINDER") == 0) Param->Geometry=CYLINDER;
	if (strcmp(Param->Geometry_type, "SPHERE") == 0) Param->Geometry=SPHERE;
	if (strcmp(Param->Geometry_type, "SLIT") == 0) Param->Geometry=SLIT;
	if (strcmp(Param->Geometry_type, "CHANNEL") == 0) Param->Geometry=CHANNEL;
	
	if (strcmp(Param->Simul_string, "TRANS") == 0) Param->Simul_type=TRANS;
	if (strcmp(Param->Simul_string, "ANCHOR") == 0) Param->Simul_type=ANCHOR;
	
	if (Param->GPU==YES){
		strcpy(Param->GPUFF, "/gpu"); 
	}else{
		strcpy(Param->GPUFF, ""); 
	} 
	
	switch(Param->Geometry){
		case CYLINDER:
		case SLIT:
		case CHANNEL:
			if(Param->N_Colloids != 0){
				printf("WARNING!!! %s Simul is not compatible with Colloids yet. Setting number of colloids to 0\n", Param->Geometry_type); fflush(NULL);
				Param->N_Colloids = 0;
			}
			if(Param->Simul_type == ANCHOR && Param->N_Proteins > 1){
				printf("WARNING!!! %s-ANCHOR Simul is not compatible with more than 1 protein. Setting number of proteins to 1\n", Param->Geometry_type); fflush(NULL);
				Param->N_Proteins = 1;
			}
			break;
		case SPHERE:
			if(Param->N_Colloids == 0){
				printf("WARNING!!! There should be at least one Colloid to run SPHERE type simulation\n"); fflush(NULL);
				exit(1);
			}
			if(Param->Simul_type == ANCHOR && Param->N_Colloids > 1){
				printf("WARNING!!! ANCHOR Simul is not compatible with more than one Colloid yet. Setting number of colloids to 1\n"); fflush(NULL);
				Param->N_Colloids = 1;
			}
			break;
	}
	
	
		Protein=(struct Proteins*)calloc(1,sizeof(struct Proteins));
		Protein_pointers=(struct Proteins**)calloc(Param->N_Proteins,sizeof(struct Proteins*));
		Protein_pointers[0]=Protein;
	
	Param->Linker_Size=2*Param->Linker_Size;
	Param->Polymer_Size=Param->Linker_Size;
	
	Param->seed=-Param->useed;
	
	
	Param->FRACTIONS[END_NEUTRAL]=Param->F_Neutral;
	Param->FRACTIONS[END_NEUTRAL_ACCEPTOR]=Param->F_Neutral_Acceptor;
	Param->FRACTIONS[END_NEUTRAL_DONOR]=Param->F_Neutral_Donor;
	Param->FRACTIONS[END_NEUTRAL_ACCEPTOR_DONOR]=Param->F_Neutral_Acceptor_Donor;
	Param->FRACTIONS[END_POSITIVE_ACCEPTOR]=Param->F_Positive_Acceptor;
	Param->FRACTIONS[END_POSITIVE_DONOR]=Param->F_Positive_Donor;
	Param->FRACTIONS[END_POSITIVE_ACCEPTOR_DONOR]=Param->F_Positive_Acceptor_Donor;
	Param->FRACTIONS[END_NEGATIVE_ACCEPTOR]=Param->F_Negative_Acceptor;
	Param->FRACTIONS[END_NEGATIVE_DONOR]=Param->F_Negative_Donor;
	Param->FRACTIONS[END_NEGATIVE_ACCEPTOR_DONOR]=Param->F_Negative_Acceptor_Donor;
	
	
	
	printf("Param Read\n");fflush(NULL);
	if(Param->N_Proteins>0){
		create_Protein(sist,Param,Protein);
		Param->Nres=Protein->Nres;
	}else{
		Param->protein_flag=NO;
	}
	
	
	
	printf("Protein Created\n");fflush(NULL);
	
		// Definining LInker length
	if(Param->Linker_Length>3){
		sist->MOL_SHIFT=(Param->Linker_Length)*Param->Linker_Size;
	}else{
		sist->MOL_SHIFT=3*Param->Linker_Size;
	}
	if(sist->MOL_SHIFT<4.5) sist->MOL_SHIFT=4.5;
	
	double polymer_scale=2*Param->Polymer_Length*Param->Polymer_Size;
	double ref_scale=polymer_scale;
	double prot_scale=0;

	
	
	
	if(Param->Box_Scale_x<0){
		sist->box_x =-Param->Box_Scale_x;		
	}else{
		if(Param->protein_flag==YES) {
			prot_scale=Param->Box_Scale_x*Protein->max_dist+2*(sist->MOL_SHIFT*1.5);
			if(prot_scale>polymer_scale){
				ref_scale=prot_scale;
			}else{
				ref_scale=polymer_scale;
			}
		}
		sist->box_x=ref_scale;
	}
	
		
	if(Param->Box_Scale_yz<0){
		switch (Param->Geometry) {
			case SPHERE:
			sist->box_y=sist->box_z=-Param->Box_Scale_x;
			break;
			case CYLINDER:
			case SLIT:
			case CHANNEL:
				sist->box_z=sist->box_y=-Param->Box_Scale_yz;
			break;
			default:
				sist->box_z=sist->box_y=-Param->Box_Scale_yz;
			break;
		}	
	}else{
		if(Param->protein_flag==YES) {
			prot_scale=Param->Box_Scale_yz*Protein->max_dist+2*(sist->MOL_SHIFT*1.5);
			if(prot_scale>polymer_scale){
				ref_scale=prot_scale;
			}else{
				ref_scale=polymer_scale;
			}	
		}
		sist->box_y=sist->box_z=ref_scale;
	}
			
	



		
	sist->Half_Box_x=sist->box_x/2.0;
	sist->Half_Box_y=sist->box_y/2.0;
	sist->Half_Box_z=sist->box_z/2.0;
	if(Param->Geometry==CYLINDER){
		sist->radius=sist->box_y/2.;
		sist->inradius=sist->radius-(sist->MOL_SHIFT);
	}else{
		sist->radius=sist->box_y*2.;
		sist->inradius=sist->radius;
	}
	double area_factor = 4.0 / (Param->Polymer_Size * Param->Polymer_Size);
	switch (Param->Geometry) {
		case SLIT:
		case CHANNEL:
			Param->N_Polymer = (int)(Param->Polymer_Fraction * area_factor * sist->box_x * sist->box_y);
			break;
		case CYLINDER:
			Param->N_Polymer = (int)(Param->Polymer_Fraction * area_factor * sist->box_x * 2 * sist->radius);
			break;
		case SPHERE:
			Param->N_Polymer = (int)(Param->Polymer_Fraction * area_factor * Param->Coll_R * Param->Coll_R);
			break;
	}
	
	
	if(Param->N_Polymer>0){
		Param->brush_flag=YES;	
		
		for(k=0;k<POL_END_TYPES;k++){
			if(Param->FRACTIONS[k]>0){
				Brush->N_Chain_Frac[k]=lround(Param->N_Polymer*Param->FRACTIONS[k]);
			}else{
				Brush->N_Chain_Frac[k]=0;
			}
		}
		Param->N_Polymer=0;
		for(k=0;k<POL_END_TYPES;k++){
			Param->N_Polymer+=Brush->N_Chain_Frac[k]; //Correctly summing all the polymer types
		}
		
			
	}else{
		Param->brush_flag=NO;
	}
		printf("CALCUALTED NUMBER OF POLYMERS Param->N_Polymer=%d; Param->Polymer_Fraction=%lf, sist->box_x=%lf, radius=%lf, Param->Polymer_Size=%lf\n",Param->N_Polymer,Param->Polymer_Fraction,sist->box_x,sist->radius,Param->Polymer_Size);fflush(NULL);
	
	printf("Parameter Read\n");fflush(NULL);
	
	
	if(Param->protein_flag==YES) protein_rotate_shift(sist,Param,Protein);
	
	if(Param->protein_flag==YES) protein_properties(sist,Param,Protein);
	
	if(Param->Simul_type==TRANS){
		/*if(Param->protein_flag==YES){
			Brush->minx=Protein->max_dist; //This leaves space at the beginning for the protein to entr the pore
			Brush->maxx=sist->box_x-Param->Polymer_Length*Param->Polymer_Size; //This leaves space at the end for the polymer brush to extend
		}else{
			Brush->minx=0; //This  doesn't leaves space at the beginning for the protein to entr the pore
			Brush->maxx=sist->box_x; //This doesn't leaves space at the end for the polymer brush to extend
		}*/
		Brush->minx=0; //This doesn't leaves space at the beginning for the protein to entr the pore
		Brush->maxx=sist->box_x; //This  doesn't leaves space at the end for the polymer brush to extend
	}
	printf("CALCUALTED BRUSH minx %lf and max %lf\n",Brush->minx,Brush->maxx);fflush(NULL);
	
	
	if(Param->protein_flag==YES){
		create_protein_copies(sist,Param,Protein_pointers);
	}
	
	switch (Param->Geometry) {
		case CYLINDER:
			create_brush_anchor_cylinder(sist, Param, Protein_pointers, Brush);
			break;
		case SLIT:
		case CHANNEL:
			create_brush_anchor_plane(sist, Param, Protein_pointers, Brush);
			break;
		case SPHERE:
			create_brush_anchor_sphere(sist, Param, Protein_pointers, Brush);
			break;
	}
	
	
	sprintf(line,"%s-out.atomparam",Param->prefix);
	fp=fopen(line,"w");
	
	write_topology(sist,Param,Protein_pointers,Brush,fp,LAMMPS_SCALE,ALL);
	
	Param->steps_flow_temp=fabs(Param->final_flow_temp-Param->init_flow_temp)/Param->bin_flow_temp;
	Param->bin_flow_temp=(Param->final_flow_temp-Param->init_flow_temp)/Param->steps_flow_temp; // In this way we have the directionality of the steps
	printf("Param->steps_flow_temp=%d\n",Param->steps_flow_temp);
	
	
	fclose(fp);
	j=0;
	int prev=0;
	for(i=1;i<=8+Param->steps_flow_temp+1;i++){	
		printf("Write step i=%d/%d\n",i,8+Param->steps_flow_temp);
		if(i>3){			
			if(i>=8){
				 if(Param->Flow_Flag==YES) {
					 sprintf(line,"in.%s-%d",Param->prefix,i);
						fp=fopen(line,"w");
					 write_infile (sist, Param,Protein_pointers,Brush,fp,prev,i);
					 fclose(fp);
					 prev=i;
				 }
			 }else{
				if(Param->protein_flag==YES){
					if(i==4){
						if((Param->Simul_type==ANCHOR)&&(Param->protein_flag==YES)){
							sprintf(line,"in.%s-%d",Param->prefix,i);
							fp=fopen(line,"w");
							write_infile (sist, Param,Protein_pointers,Brush,fp,prev,i);
							fclose(fp);
							prev=i;
						}
					}else{
						sprintf(line,"in.%s-%d",Param->prefix,i);
						fp=fopen(line,"w");
						write_infile (sist, Param,Protein_pointers,Brush,fp,prev,i);
						fclose(fp);
						prev=i;
					}
				}
			 }
		}else{
			sprintf(line,"in.%s-%d",Param->prefix,i);
			fp=fopen(line,"w");
			write_infile (sist, Param,Protein_pointers,Brush,fp,prev,i);
			fclose(fp);
			prev=i;
		}
		
	}
	
	sprintf(line,"%s-out-xtc-all.topology",Param->prefix);
	fp=fopen(line,"w");
	
	write_topology(sist,Param,Protein_pointers,Brush,fp,XTC_SCALE,ALL);
	
	fclose(fp);
	
	
	sprintf(line,"%s-out-xtc-flow.topology",Param->prefix);
	fp=fopen(line,"w");
	
	write_topology(sist,Param,Protein_pointers,Brush,fp,XTC_SCALE,NO_POL_ANCHORS);
	
	fclose(fp);
	
	sprintf(line,"%s-out-protein.topology",Param->prefix);
	fp=fopen(line,"w");
	
	write_topology(sist,Param,Protein_pointers,Brush,fp,LAMMPS_SCALE,PROTEIN);
	
	fclose(fp);
	
	
	sprintf(line,"%s-out_simul.tcl",Param->prefix);
	fp=fopen(line,"w");
	
	write_tcl_header(sist,Param,fp,LAMMPS_SCALE);
	fclose(fp);
	
	sprintf(line,"%s-out_protein.tcl",Param->prefix);
	fp=fopen(line,"w");
	
	write_tcl_protein(sist,Param,fp,LAMMPS_SCALE);
	fclose(fp);
	
	sprintf(line,"%s-out_all.tcl",Param->prefix);
	fp=fopen(line,"w");
	
	write_tcl_all(sist,Param,fp,XTC_SCALE);
	fclose(fp);
	
	sprintf(line,"%s-out_flow.tcl",Param->prefix);
	fp=fopen(line,"w");
	write_tcl_flow(sist,Param,fp,XTC_SCALE);
	
	fclose(fp);
	return(0);
}


void initialize_parameters (struct Interaction_Param *Param,struct param **ptr_parameters, unsigned int *n_param){
	int err5=1;
	
	*n_param = 42; 
	struct param * parameters;
	parameters = (struct param*) calloc (*n_param,sizeof(struct param));
	*ptr_parameters = parameters;
	if (Param->pdb_filename == NULL )
		{
		Param->pdb_filename = (char*)calloc(2048,sizeof(char));
	}
	if (Param->prefix == NULL )
		{
		Param->prefix = (char*)calloc(2048,sizeof(char));
	}
	if (Param->mask == NULL )
		{
		Param->mask = (char*)calloc(2048,sizeof(char));
	}
	if (Param->mask_wall_file == NULL )
		{
		Param->mask_wall_file = (char*)calloc(2048,sizeof(char));
	}
	if (Param->mask_prot == NULL )
		{
		Param->mask_prot = (char*)calloc(2048,sizeof(char));
	}
	if (Param->ignore == NULL )
		{
		Param->ignore = (char*)calloc(2048,sizeof(char));
	}
	if (Param->Geometry_type == NULL )
		{
		Param->Geometry_type = (char*)calloc(2048,sizeof(char));
	}
	if (Param->Simul_string == NULL )
		{
		Param->Simul_string = (char*)calloc(2048,sizeof(char));
	}
	
	unsigned  int i=0;
	parameters [i].name    = "E_Scale";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->E_Scale;
	*parameters[i].u.dval = 4;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Random_Rotation";
	parameters [i].flag    = 1;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->random_rotation_flag;
	*parameters[i].u.ival = YES;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Linker_Length";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->Linker_Length;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Flow";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->Flow_Flag;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Linker_Size";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Linker_Size;
	*parameters[i].u.dval = 3;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Polymer_Fraction";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Polymer_Fraction;
	*parameters[i].u.dval = 3;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "N_Proteins";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->N_Proteins;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Polymer_Length";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->Polymer_Length;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}	
	parameters [i].name    = "Polymer_Size";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Polymer_Size;
	*parameters[i].u.dval = 3;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Pressure";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->pressure;
	*parameters[i].u.dval = 0.0001;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "N_Colloids";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->N_Colloids;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Colloid_Radius";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Coll_R;
	*parameters[i].u.dval = 20;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Init_Temp";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->init_temp;
	*parameters[i].u.dval = 1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Final_Temp";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->final_temp;
	*parameters[i].u.dval = 1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Init_Flow_Temp";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->init_flow_temp;
	*parameters[i].u.dval = 1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Final_Flow_Temp";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->final_flow_temp;
	*parameters[i].u.dval = 1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Bin_Flow_Temp";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->bin_flow_temp;
	*parameters[i].u.dval = 1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "chi_s";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->chi_s;
	*parameters[i].u.dval = -1;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Protein_PDB";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->pdb_filename;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Geometry";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->Geometry_type;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Simul_Type";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->Simul_string;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	
	parameters [i].name    = "Protein_Mask";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->mask;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Protein_Mask_Wall";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->mask_wall_file;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Protein_Mask_Protonation";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->mask_prot;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Ignore_Bond_List";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->ignore;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Prefix";
	parameters [i].flag    = 0;
	parameters [i].type    = P_STRING;
	parameters [i].u.sval  = &Param->prefix;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Box_Scale_x";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Box_Scale_x;
	*parameters[i].u.dval = -1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Box_Scale_yz";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Box_Scale_yz;
	*parameters[i].u.dval = -1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Wall_Attr";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Wall_Attr;
	*parameters[i].u.dval = 10.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Wall_Attr_Gen";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->Wall_Attr_Gen;
	*parameters[i].u.dval = 4.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "Seed";
	parameters [i].flag    = 1;
	parameters [i].type    = P_UINT;
	parameters [i].u.uval  = &Param->useed;
	*parameters[i].u.uval = 1; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Neutral";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Neutral;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	
	parameters [i].name    = "F_Neutral_Acceptor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Neutral_Acceptor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Neutral_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Neutral_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Neutral_Acceptor_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Neutral_Acceptor_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Positive_Acceptor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Positive_Acceptor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Positive_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Positive_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Positive_Acceptor_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Positive_Acceptor_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Negative_Acceptor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Negative_Acceptor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Negative_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Negative_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "F_Negative_Acceptor_Donor";
	parameters [i].flag    = 1;
	parameters [i].type    = P_DOUBLE;
	parameters [i].u.dval  = &Param->F_Negative_Acceptor_Donor;
	*parameters[i].u.dval = 1.0; //Default value used in optmization
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
	parameters [i].name    = "GPU";
	parameters [i].flag    = 0;
	parameters [i].type    = P_INT;
	parameters [i].u.ival  = &Param->GPU;
	i++;
	if(*n_param<i){ printf("Too many paramters n_param= %u i= %u at code line %d\n",*n_param,i,__LINE__);
		fflush(NULL);
		exit(1);
	}
}


void Angle_generator (struct GlobalSistem *sist, struct Proteins* Protein){
	int pr,i,prt,j,k;
	struct part *Temp_Chain=NULL;
	double theta,alpha;
	double dx,dy,dz,r2;
	double dx1,dy1,dz1,r;
	int lenght=20;
	
	Temp_Chain=(struct part *)calloc(lenght,sizeof(struct part));
	
	
	
	Temp_Chain[0].x=0.0;
	Temp_Chain[0].y=0.0;
	Temp_Chain[0].z=0.0;
	Temp_Chain[0].id=ATOM_N;
	
	alpha=45.0;
	for(i=1;i<lenght;i++){
		
		
		if(i<lenght){
			theta=((90.0-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[0]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[0]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_CA;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+111.0-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[1]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[1]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_C;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+111.0+121.1-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[2]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[2]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_O;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+121.9+119.5-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-4].x+sist->Bbond[4]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-4].y+sist->Bbond[4]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-4].z;
			Temp_Chain[i].id=ATOM_H;
			
			Protein->FIRST_AtomHx=Temp_Chain[i].x;
			Protein->FIRST_AtomHy=Temp_Chain[i].y;
			Protein->FIRST_AtomHz=Temp_Chain[i].z;
			
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0-115.6+111.0-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-3].x+sist->Bbond[3]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-3].y+sist->Bbond[3]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-3].z;
			Temp_Chain[i].id=ATOM_N;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+121.9-115.6+111.0-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[0]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[0]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_CA;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+(121.9-115.6)-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[1]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[1]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_C;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((123.2-90.0+121.9-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-1].x+sist->Bbond[2]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-1].y+sist->Bbond[2]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-1].z;
			Temp_Chain[i].id=ATOM_O;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((-90.0+118.2+(121.9-115.6+111.0)-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-4].x+sist->Bbond[4]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-4].y+sist->Bbond[4]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-4].z;
			Temp_Chain[i].id=ATOM_H;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+121.9-alpha)/180.0)*PI;
			Temp_Chain[i].x=Temp_Chain[i-3].x+sist->Bbond[3]*cos(theta);
			Temp_Chain[i].y=Temp_Chain[i-3].y+sist->Bbond[3]*sin(theta);
			Temp_Chain[i].z=Temp_Chain[i-3].z;
			Temp_Chain[i].id=ATOM_N;
			//printf("%d %d %lf\n",i,Temp_Chain[i].id,theta);
			
		}
		//printf("##############\n");
	}
	
	
	
	
	dx=(Temp_Chain[ATOM_O].x-Temp_Chain[ATOM_C].x);
	dy=(Temp_Chain[ATOM_O].y-Temp_Chain[ATOM_C].y);
	dz=(Temp_Chain[ATOM_O].z-Temp_Chain[ATOM_C].z);
	
	dx1=(Temp_Chain[ATOM_H+NATOM].x-Temp_Chain[ATOM_C].x);
	dy1=(Temp_Chain[ATOM_H+NATOM].y-Temp_Chain[ATOM_C].y);
	dz1=(Temp_Chain[ATOM_H+NATOM].z-Temp_Chain[ATOM_C].z);
	
	
	sist->DotProd[1]=dx*dx1+dy*dy1+dz*dz1;//Angle OCH
	
	
	dx=(Temp_Chain[ATOM_N+NATOM].x-Temp_Chain[ATOM_C].x);
	dy=(Temp_Chain[ATOM_N+NATOM].y-Temp_Chain[ATOM_C].y);
	dz=(Temp_Chain[ATOM_N+NATOM].z-Temp_Chain[ATOM_C].z);
	
	
	
	sist->DotProd[0]=dx*dx1+dy*dy1+dz*dz1;// Angle HCN
	
	dx=(Temp_Chain[ATOM_CA].x-Temp_Chain[ATOM_N].x);
	dy=(Temp_Chain[ATOM_CA].y-Temp_Chain[ATOM_N].y);
	dz=(Temp_Chain[ATOM_CA].z-Temp_Chain[ATOM_N].z);
	
	dx1=(Temp_Chain[ATOM_C].x-Temp_Chain[ATOM_N].x);
	dy1=(Temp_Chain[ATOM_C].y-Temp_Chain[ATOM_N].y);
	dz1=(Temp_Chain[ATOM_C].z-Temp_Chain[ATOM_N].z);
	
	sist->DotProd[2]=dx*dx1+dy*dy1+dz*dz1;// Angle CNCa 
	
	
	
	return;
}

void Ref (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	int pr,i,prt,j,k;
	double theta,alpha;
	double dx,dy,dz,r2;
	double dx1,dy1,dz1,r;
	int lenght=20;
	//Atom=(struct part *)calloc(lenght,sizeof(struct part));
	
	
	lenght=Protein->N;
	Protein->Atom[0].x=0.0;
	Protein->Atom[0].y=0.0;
	Protein->Atom[0].z=0.0;
	Protein->Atom[0].id=ATOM_N;
	
	alpha=45.0;
	for(i=1;i<lenght;i++){
		
		
		if(i<lenght){
			theta=((90.0-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[0]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[0]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_CA;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+111.0-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[1]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[1]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_C;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+111.0+121.1-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[2]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[2]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_O;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+121.9+119.5-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-4].x+sist->Bbond[4]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-4].y+sist->Bbond[4]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-4].z;
			Protein->Atom[i].id=ATOM_H;
			
			Protein->FIRST_AtomHx=Protein->Atom[i].x;
			Protein->FIRST_AtomHy=Protein->Atom[i].y;
			Protein->FIRST_AtomHz=Protein->Atom[i].z;
			
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0-115.6+111.0-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-3].x+sist->Bbond[3]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-3].y+sist->Bbond[3]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-3].z;
			Protein->Atom[i].id=ATOM_N;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+121.9-115.6+111.0-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[0]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[0]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_CA;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((90.0+(121.9-115.6)-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[1]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[1]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_C;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((123.2-90.0+121.9-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-1].x+sist->Bbond[2]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-1].y+sist->Bbond[2]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-1].z;
			Protein->Atom[i].id=ATOM_O;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((-90.0+118.2+(121.9-115.6+111.0)-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-4].x+sist->Bbond[4]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-4].y+sist->Bbond[4]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-4].z;
			Protein->Atom[i].id=ATOM_H;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			i++;
		}
		if(i<lenght){
			theta=((270.0+121.9-alpha)/180.0)*PI;
			Protein->Atom[i].x=Protein->Atom[i-3].x+sist->Bbond[3]*cos(theta);
			Protein->Atom[i].y=Protein->Atom[i-3].y+sist->Bbond[3]*sin(theta);
			Protein->Atom[i].z=Protein->Atom[i-3].z;
			Protein->Atom[i].id=ATOM_N;
			//printf("%d %d %lf\n",i,Protein->Atom[i].id,theta);
			
		}
		//printf("##############\n");
	}
	
	
	
	
	dx=(Protein->Atom[ATOM_O].x-Protein->Atom[ATOM_C].x);
	dy=(Protein->Atom[ATOM_O].y-Protein->Atom[ATOM_C].y);
	dz=(Protein->Atom[ATOM_O].z-Protein->Atom[ATOM_C].z);
	
	dx1=(Protein->Atom[ATOM_H+NATOM].x-Protein->Atom[ATOM_C].x);
	dy1=(Protein->Atom[ATOM_H+NATOM].y-Protein->Atom[ATOM_C].y);
	dz1=(Protein->Atom[ATOM_H+NATOM].z-Protein->Atom[ATOM_C].z);
	
	
	sist->DotProd[1]=dx*dx1+dy*dy1+dz*dz1;//Angle OCH
	
	
	dx=(Protein->Atom[ATOM_N+NATOM].x-Protein->Atom[ATOM_C].x);
	dy=(Protein->Atom[ATOM_N+NATOM].y-Protein->Atom[ATOM_C].y);
	dz=(Protein->Atom[ATOM_N+NATOM].z-Protein->Atom[ATOM_C].z);
	
	
	
	sist->DotProd[0]=dx*dx1+dy*dy1+dz*dz1;// Angle HCN
	
	dx=(Protein->Atom[ATOM_CA].x-Protein->Atom[ATOM_N].x);
	dy=(Protein->Atom[ATOM_CA].y-Protein->Atom[ATOM_N].y);
	dz=(Protein->Atom[ATOM_CA].z-Protein->Atom[ATOM_N].z);
	
	dx1=(Protein->Atom[ATOM_C].x-Protein->Atom[ATOM_N].x);
	dy1=(Protein->Atom[ATOM_C].y-Protein->Atom[ATOM_N].y);
	dz1=(Protein->Atom[ATOM_C].z-Protein->Atom[ATOM_N].z);
	
	sist->DotProd[2]=dx*dx1+dy*dy1+dz*dz1;// Angle CNCa 
	
	
	
	return;
}

void Rescaler (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	int pr,i,prt,j,k;
	double dx,dy,dz,r2,r;
	double Etouch,EtouchO;
	
	
	
	//printf("################## Rescaling\n");
	
	
	for(k=0;k<Protein->N;k++){
		if(Protein->Atom[k].id==ATOM_H){
			
			dx=(Protein->Atom[k-4].x-Protein->Atom[k].x);
			dy=(Protein->Atom[k-4].y-Protein->Atom[k].y);
			dz=(Protein->Atom[k-4].z-Protein->Atom[k].z);
			//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
			r2=(dx*dx+ dy*dy + dz*dz);
			r=sqrt(r2);
			if(fabs(r-sist->Bbond[4])>1e-10){
				//printf("H-N rescale  %d %30.20lf \n",k,r-sist->Bbond[4]);fflush(NULL);
				Protein->Atom[k].x-=Protein->Atom[k-4].x;
				Protein->Atom[k].y-=Protein->Atom[k-4].y;
				Protein->Atom[k].z-=Protein->Atom[k-4].z;
				
				Protein->Atom[k].x*=sist->Bbond[4]/r;
				Protein->Atom[k].y*=sist->Bbond[4]/r;
				Protein->Atom[k].z*=sist->Bbond[4]/r;
				
				Protein->Atom[k].x+=Protein->Atom[k-4].x;
				Protein->Atom[k].y+=Protein->Atom[k-4].y;
				Protein->Atom[k].z+=Protein->Atom[k-4].z;
				
				
				
				dx=(Protein->Atom[k-4].x-Protein->Atom[k].x);
				dy=(Protein->Atom[k-4].y-Protein->Atom[k].y);
				dz=(Protein->Atom[k-4].z-Protein->Atom[k].z);
				//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
				r2=(dx*dx+ dy*dy + dz*dz);
				r=sqrt(r2);
				if(fabs(r-sist->Bbond[4])>1e-10){
					printf("H-N wrong  %d %30.20lf \n",k,r-sist->Bbond[4]);fflush(NULL);
					exit(EXIT_FAILURE);
				}
			}
			
		}
		
		if(Protein->Atom[k].id==ATOM_O){
			dx=(Protein->Atom[k-1].x-Protein->Atom[k].x);
			dy=(Protein->Atom[k-1].y-Protein->Atom[k].y);
			dz=(Protein->Atom[k-1].z-Protein->Atom[k].z);
			//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
			r2=(dx*dx+ dy*dy + dz*dz);
			r=sqrt(r2);
			if(fabs(r-sist->Bbond[2])>1e-10){
				//printf("O-C rescale  %d %30.20lf \n",k,r-sist->Bbond[2]);fflush(NULL);
				Protein->Atom[k].x-=Protein->Atom[k-1].x;
				Protein->Atom[k].y-=Protein->Atom[k-1].y;
				Protein->Atom[k].z-=Protein->Atom[k-1].z;
				
				Protein->Atom[k].x*=sist->Bbond[2]/r;
				Protein->Atom[k].y*=sist->Bbond[2]/r;
				Protein->Atom[k].z*=sist->Bbond[2]/r;
				
				Protein->Atom[k].x+=Protein->Atom[k-1].x;
				Protein->Atom[k].y+=Protein->Atom[k-1].y;
				Protein->Atom[k].z+=Protein->Atom[k-1].z;
				
				
				
				dx=(Protein->Atom[k-1].x-Protein->Atom[k].x);
				dy=(Protein->Atom[k-1].y-Protein->Atom[k].y);
				dz=(Protein->Atom[k-1].z-Protein->Atom[k].z);
				//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
				r2=(dx*dx+ dy*dy + dz*dz);
				r=sqrt(r2);
				if(fabs(r-sist->Bbond[2])>1e-10){
					printf("O-C wrong  %d %30.20lf \n",k,r-sist->Bbond[2]);fflush(NULL);
					exit(EXIT_FAILURE);
				}
			}
		}
		if(Protein->Atom[k].id==ATOM_N){
			
			dx=(Protein->Atom[k+1].x-Protein->Atom[k].x);
			dy=(Protein->Atom[k+1].y-Protein->Atom[k].y);
			dz=(Protein->Atom[k+1].z-Protein->Atom[k].z);
			//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
			r2=(dx*dx+ dy*dy + dz*dz);
			r=sqrt(r2);
			if(fabs(r-sist->Bbond[0])>1e-10){
				//printf("CA-N rescale  %d %30.20lf \n",k,r-sist->Bbond[0]);fflush(NULL);
				dx*=1.0-sist->Bbond[0]/r;
				dy*=1.0-sist->Bbond[0]/r;
				dz*=1.0-sist->Bbond[0]/r;
				
				for(j=k+1;j<Protein->N;j++){
					Protein->Atom[j].x-=dx;
					Protein->Atom[j].y-=dy;
					Protein->Atom[j].z-=dz;
				}
				
				dx=(Protein->Atom[k+1].x-Protein->Atom[k].x);
				dy=(Protein->Atom[k+1].y-Protein->Atom[k].y);
				dz=(Protein->Atom[k+1].z-Protein->Atom[k].z);
				//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
				r2=(dx*dx+ dy*dy + dz*dz);
				r=sqrt(r2);
				if(fabs(r-sist->Bbond[0])>1e-10){
					printf("CA-N wrong  %d %30.20lf \n",k,r-sist->Bbond[0]);fflush(NULL);
					exit(EXIT_FAILURE);
				}
			}
		}
		if(Protein->Atom[k].id==ATOM_CA){
			
			dx=(Protein->Atom[k+1].x-Protein->Atom[k].x);
			dy=(Protein->Atom[k+1].y-Protein->Atom[k].y);
			dz=(Protein->Atom[k+1].z-Protein->Atom[k].z);
			//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
			r2=(dx*dx+ dy*dy + dz*dz);
			r=sqrt(r2);
			if(fabs(r-sist->Bbond[1])>1e-10){
				//printf("C-CA rescale  %d %30.20lf \n",k,r-sist->Bbond[1]);fflush(NULL);
				
				dx*=1.0-sist->Bbond[1]/r;
				dy*=1.0-sist->Bbond[1]/r;
				dz*=1.0-sist->Bbond[1]/r;
				
				for(j=k+1;j<Protein->N;j++){
					Protein->Atom[j].x-=dx;
					Protein->Atom[j].y-=dy;
					Protein->Atom[j].z-=dz;
				}
				
				
				dx=(Protein->Atom[k+1].x-Protein->Atom[k].x);
				dy=(Protein->Atom[k+1].y-Protein->Atom[k].y);
				dz=(Protein->Atom[k+1].z-Protein->Atom[k].z);
				//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
				r2=(dx*dx+ dy*dy + dz*dz);
				r=sqrt(r2);
				if(fabs(r-sist->Bbond[1])>1e-10){
					printf("C-CA wrong  %d %30.20lf \n",k,r-sist->Bbond[1]);fflush(NULL);
					exit(EXIT_FAILURE);
				}
			}
			
		}
		if(k>0){//otherwise has problems with the last N where k+3 >Plength
			if(Protein->Atom[k].id==ATOM_N){
				
				dx=(Protein->Atom[k-3].x-Protein->Atom[k].x);
				dy=(Protein->Atom[k-3].y-Protein->Atom[k].y);
				dz=(Protein->Atom[k-3].z-Protein->Atom[k].z);
				//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
				r2=(dx*dx+ dy*dy + dz*dz);
				r=sqrt(r2);
				if(fabs(r-sist->Bbond[3])>1e-10){
					//printf("C-N rescale  %d %30.20lf \n",k,r-sist->Bbond[3]);fflush(NULL);
					
					dx*=1.0-sist->Bbond[3]/r;
					dy*=1.0-sist->Bbond[3]/r;
					dz*=1.0-sist->Bbond[3]/r;
					for(j=k;j<Protein->N;j++){
						Protein->Atom[j].x+=dx;
						Protein->Atom[j].y+=dy;
						Protein->Atom[j].z+=dz;
					}
					
					
					dx=(Protein->Atom[k-3].x-Protein->Atom[k].x);
					dy=(Protein->Atom[k-3].y-Protein->Atom[k].y);
					dz=(Protein->Atom[k-3].z-Protein->Atom[k].z);
					//dx=P_Dist(dx);dy=P_Dist(dy);dz=P_Dist(dz);
					r2=(dx*dx+ dy*dy + dz*dz);
					r=sqrt(r2);
					if(fabs(r-sist->Bbond[3])>1e-10){
						printf("C-N wrong  %d %30.20lf \n",k,r-sist->Bbond[3]);fflush(NULL);
						exit(EXIT_FAILURE);
					}
				}
			}
		}
		
		
		
	}
	//printf("################## Rescaling ok\n");
	
	
	
	
	
	return;
	
	
}


double calc_dihedralf_angle_LAMMPS(double *p1, double *p2, double *p3, double *p4){
	int i;
	
	double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
	double cosphi,sinphi;
	// modified new variable to store angle
	double anglephi;
	// 1st bond
	vb1x = p1[0] - p2[0];
	vb1y = p1[1] - p2[1];
	vb1z = p1[2] - p2[2];
	
	// 2nd bond
	
	vb2x = p3[0] - p2[0];
	vb2y = p3[1] - p2[1];
	vb2z = p3[2] - p2[2];
	
	vb2xm = -vb2x;
	vb2ym = -vb2y;
	vb2zm = -vb2z;
	
	// 3rd bond
	
	vb3x = p4[0] - p3[0];
	vb3y = p4[1] - p3[1];
	vb3z = p4[2] - p3[2];
	
	// c,s calculation
	
	ax = vb1y*vb2zm - vb1z*vb2ym;
	ay = vb1z*vb2xm - vb1x*vb2zm;
	az = vb1x*vb2ym - vb1y*vb2xm;
	bx = vb3y*vb2zm - vb3z*vb2ym;
	by = vb3z*vb2xm - vb3x*vb2zm;
	bz = vb3x*vb2ym - vb3y*vb2xm;
	
	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
	rg = sqrt(rgsq);
	
	rginv = ra2inv = rb2inv = 0.0;
	if (rg > 0) rginv = 1.0/rg;
	if (rasq > 0) ra2inv = 1.0/rasq;
	if (rbsq > 0) rb2inv = 1.0/rbsq;
	rabinv = sqrt(ra2inv*rb2inv);
	
	cosphi = (ax*bx + ay*by + az*bz)*rabinv;
	sinphi = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);
	
	//added to see if it can be used as criterium
	// printf("The sin is %lf\n",sinphi);
	/* Calculate dihedral angle */
	
	//modified to have right 2pi range
	anglephi=acos(cosphi);
	if (sinphi < 0){anglephi=-anglephi;}
	return (anglephi);
	//if( scalar(aXb, c) < 0.0 ) *phi = (2.0*PI) - *phi;
	
}






void Sheadecoder (char *Dec, int N){
	int i;
	
	for(i=0;i<N;i++){
		switch (Dec[i]) {
			case 'B':
			Dec[i]=1;
			break;
			case 'L':
			Dec[i]=2;
			break;
			case 'N':
			Dec[i]=3;
			break;
			default:
			printf("errore nel codice della Senquenza Residuo[%d]=%c unkwon\n",i,Dec[i]);
			fflush(NULL);
			exit(0);
		} 
	}
}




void fastadecoder (int *Dec, char *Enc, int N){
	int group=0,i=0;
	int Seq_Length=0;
	int err5=1;
	
	
	for(i=0;i<N;i++){
		switch (Enc[i]) {
			case 'A':
			Dec[i]=1;
			break;
			case 'C':
			Dec[i]=2;
			break;
			case 'D':
			Dec[i]=3;
			break;
			case 'E':
			Dec[i]=4;
			break;
			case 'F':
			Dec[i]=5;
			break;
			case 'G':
			Dec[i]=6;
			break;
			case 'H':
			Dec[i]=7;
			break;
			case 'I':
			Dec[i]=8;
			break;
			case 'K':
			Dec[i]=9;
			break;
			case 'L':
			Dec[i]=10;
			break;
			case 'M':
			Dec[i]=11;
			break;
			case 'N':
			Dec[i]=12;
			break;
			case 'P':
			Dec[i]=13;
			break;
			case 'Q':
			Dec[i]=14;
			break;
			case 'R':
			Dec[i]=15;
			break;
			case 'S':
			Dec[i]=16;
			break;
			case 'T':
			Dec[i]=17;
			break;
			case 'V':
			Dec[i]=18;
			break;
			case 'W':
			Dec[i]=19;
			break;
			case 'Y':
			Dec[i]=20;
			break;
			case 'X':
			Dec[i]=1;
			break;
			default:
			printf("%s:%d Errore nel codice della Senquenza Residuo Enc[%d]=%c unkwon\n",__FILE__,__LINE__,i,Enc[i]);
			fflush(NULL);
			exit(0);
		}
	}
	
}


/**
* @brief Rotate a particle around a give axis
* @param[in] sist: structure containing the simulation parameters and status
* @param[in,out] atom: particle to rotate
* @param[in] alpha: rotation angle
* @param[in] u1,u2,u3: x,y and z components of rotation axis u
* @param[in] flag: tells the function whether to compute the rotation matrix (1)
* or used the stored values (0).
*
* @warning: there is no check on whether the matrix has been stored already or not.
* @warning: performs a rotation of amplitude alpha (previously was 2*alpha)
*/
/*
*/
void Rotation (double *x,double *y,double *z,double alpha,double u1,double u2,double u3,int flag)
	{
	double v1,v2,v3;
	double o1,o2,o3;
	static double a11,a12,a13,a21,a22,a23,a31,a32,a33;
	double nx,ny,nz;
	int err5=1;
	
	double theta=alpha;
	double CsTh, SnTh,OnemCsTh;
	static double modulo2,modulo;
	
	//let's try to reduce numerical inaccuracies..(LT)
	v1=*x;
	v2=*y;
	v3=*z;
	//printf("BEFORE INSIDE|| %lf %lf %lf\n",v1,v2,v3);fflush(NULL);
	if(flag==1){
		CsTh = cos(theta);
		SnTh = sin(theta);
		OnemCsTh = 1 - CsTh;
		
		
		modulo2=u1*u1+u2*u2+u3*u3;
		modulo=sqrt(modulo2);
		u1=u1/modulo;
		u2=u2/modulo;
		u3=u3/modulo;
		
		//This matrix is multiplied by modulo, one must then divide by it.
		/*
		a11= modulo*CsTh+u1*u1*OnemCsTh/modulo;
		a12= (u1*u2*OnemCsTh/modulo-u3*SnTh);
		a13= (u1*u3*OnemCsTh/modulo+u2*SnTh);
		
		a21= (u1*u2*OnemCsTh/modulo+u3*SnTh);
		a22= modulo*CsTh+u2*u2*OnemCsTh/modulo;
		a23= (u2*u3*OnemCsTh/modulo-u1*SnTh);
		
		a31= (u1*u3*OnemCsTh/modulo-u2*SnTh);
		a32= (u2*u3*OnemCsTh/modulo+u1*SnTh);
		a33= modulo*CsTh+u3*u3*OnemCsTh/modulo;
		*/
		a11= CsTh+u1*u1*OnemCsTh;
		a12= (u1*u2*OnemCsTh-u3*SnTh);
		a13= (u1*u3*OnemCsTh+u2*SnTh);
		
		a21= (u1*u2*OnemCsTh+u3*SnTh);
		a22= CsTh+u2*u2*OnemCsTh;
		a23= (u2*u3*OnemCsTh-u1*SnTh);
		
		a31= (u1*u3*OnemCsTh-u2*SnTh);
		a32= (u2*u3*OnemCsTh+u1*SnTh);
		a33= CsTh+u3*u3*OnemCsTh;
	}
	
	o1=a11*v1+a12*v2+a13*v3;
	o2=a21*v1+a22*v2+a23*v3;
	o3=a31*v1+a32*v2+a33*v3;
	//printf("AFTER INSIDE|| %lf %lf %lf\n",o1,o2,o3);fflush(NULL);
	*x=o1;
	*y=o2;
	*z=o3;
	//printf("TRANSFER INSIDE|| %lf %lf %lf\n",*x,*y,*z);fflush(NULL);
	//atom->x=nx/modulo;
	//	atom->y=ny/modulo;
	//	atom->z=nz/modulo;
	
	return;
}
void Init_Topology (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins *Protein, struct Brushes *Brush, int type){
	int idx,i,k;
	sist->Tot_Atoms=0;
	/******************* ATOMS *****************/
	if(Param->protein_flag==YES) sist->Tot_Atoms=Param->N_Proteins*PROTEIN+Param->N_Proteins*(Protein->NMasked*(Param->Linker_Length));
	
	int polymer_atoms = (type == ALL) ? Param->N_Polymer * Param->Polymer_Length : Param->N_Polymer * (Param->Polymer_Length - 1);
	switch (Param->Geometry) {
		case SPHERE:
			polymer_atoms += Param->N_Colloids; // Includes Anchors
			sist->Tot_Atoms += polymer_atoms * Param->N_Colloids;
			break;
		case CYLINDER:
		case SLIT:
		case CHANNEL:
			sist->Tot_Atoms += polymer_atoms;
			break;
	}
	printf("Calculated number of Atoms\n");
	printf("Protein Atoms %d\n",PROTEIN);
	printf("Linker Atoms %d*%d=%d\n",Protein->NMasked,(Param->Linker_Length),Protein->NMasked*(Param->Linker_Length));
	printf("Brush Atoms %d*%d=%d\n",Param->N_Polymer,Param->Polymer_Length,Param->N_Polymer*Param->Polymer_Length);
	fflush(NULL);
	sist->Tot_Bonds=0;
	//PROTEIN BONDS
	if(Param->protein_flag==YES) {
		sist->Tot_Bonds=Protein->Nres-1;
		sist->Tot_Bonds-=Protein->NSkip;
		sist->Tot_Bonds*=Param->N_Proteins;
	}
	
	//LINKER BONDS
	sist->Tot_Bonds+=Protein->NMasked*Param->Linker_Length*Param->N_Proteins;
	//BRUSH BONDS
	switch(Param->Geometry){
		case SPHERE:
			if(type==ALL) sist->Tot_Bonds+=Param->N_Polymer*(Param->Polymer_Length-1)*Param->N_Colloids; //Includes Anchors
			if(type==NO_POL_ANCHORS) sist->Tot_Bonds+=Param->N_Polymer*(Param->Polymer_Length-2)*Param->N_Colloids;
		break;
		case CHANNEL:
		case SLIT:
		case CYLINDER:
			if(type==ALL) sist->Tot_Bonds+=Param->N_Polymer*(Param->Polymer_Length-1); //Includes Anchors
			if(type==NO_POL_ANCHORS) sist->Tot_Bonds+=Param->N_Polymer*(Param->Polymer_Length-2);
		break;
	}
	sist->Tot_Dihedrals=0;
	sist->Tot_Angles=0;
	if(Param->protein_flag==YES){
		for(i=ATOM_CA;i<Protein->N;i+=NATOM){
			if(i+3*NATOM<Protein->N){
				if((Protein->Ignore[i]==0)&&(Protein->Ignore[i+NATOM]==0)&&(Protein->Ignore[i+2*NATOM]==0)&&(Protein->Ignore[i+3*NATOM]==0)){
					sist->Tot_Dihedrals++;
				}
			}
			if(i+2*NATOM<Protein->N){
				if((Protein->Ignore[i]==0)&&(Protein->Ignore[i+NATOM]==0)&&(Protein->Ignore[i+2*NATOM]==0)){
					sist->Tot_Angles++;
				}
			}
		}
	}
	
	sist->Tot_Dihedrals_Types=sist->Tot_Dihedrals;
	sist->Tot_Angles_Types=sist->Tot_Angles;
	
	
	sist->Tot_Dihedrals*=Param->N_Proteins;
	sist->Tot_Angles*=Param->N_Proteins;
	
	sist->tot_atom_type=SRD+POL_END_TYPES;
	
	
	if(Param->protein_flag==YES) sist->tot_bonds_type=Protein->Nres-1;
	if(Param->protein_flag==YES) sist->tot_bonds_type-=Protein->NSkip;
	sist->tot_bonds_type++;
	sist->polymers_bonds_type=sist->tot_bonds_type;
	

	
	return;
}

void create_brush_anchor_cylinder (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush){
	int idx,i,k,m;
	double theta;
	double dist=0;
	double dx,dy,dz;
	int flag=0;

	idx=0;
	Brush->Anchors=(struct part*)calloc(Param->N_Polymer,sizeof(struct part));
	for(k=0;k<Param->N_Polymer;k++){
		//Place anchor points on the surface of the cylinder
		theta=ran3(&Param->seed)*2*PI;
		Brush->Anchors[idx].x=ran3(&Param->seed)*(sist->box_x-4*Param->Polymer_Size)+2*Param->Polymer_Size;		
		Brush->Anchors[idx].y=(sist->inradius+Param->Polymer_Size*1.1/2)*cos(theta);
		Brush->Anchors[idx].z=(sist->inradius+Param->Polymer_Size*1.1/2)*sin(theta);
		Brush->Anchors[idx].flag=-1;
		//Test overalp with Protein
		flag=1;
		switch(Param->Simul_type){
			case TRANS:
				if((Brush->Anchors[idx].x<Brush->minx)||(Brush->Anchors[idx].x>Brush->maxx)) flag=0; // This removes the polymers at the entrnace and exit of the pore to allow the entrance of the protein
			if((Param->mask_flag==YES)&&(Param->protein_flag==YES)){
					
					for(m=0;m<Param->N_Proteins;m++){
						for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
							dx=(Brush->Anchors[idx].x-Protein_pointers[m]->Atom[i].x); //Remove the brush before the proteins initial position
							


							if(dx<0) flag=0;

						}
					}
				}
				
				
				
			break;
			case ANCHOR:
				if((Param->mask_flag==YES)&&(Param->protein_flag==YES)){
					
					for(m=0;m<Param->N_Proteins;m++){
						for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
							dx=(Protein_pointers[m]->Atom[i].x-Brush->Anchors[idx].x);
							dy=(Protein_pointers[m]->Atom[i].y-Brush->Anchors[idx].y-sist->Half_Box_y);


							dist=sqrt(dx*dx+dy*dy);


							if((dist<Param->Polymer_Size*2)&&(Brush->Anchors[idx].z<0)) flag=0;

						}
					}
				}
			break;		
		}
		if(flag==1){
			Brush->Anchors[idx].lammps_idx=idx*Param->Polymer_Length;	
			Brush->Anchors[idx].flag=-1;
			idx++;
		}
		
	}
	Brush->N_Anchor=Param->N_Polymer=idx;
	
	printf("Before Fraction Param->N_Polymer=%d\n",Param->N_Polymer);fflush(NULL);
	
	for(k=0;k<POL_END_TYPES;k++){
		if(Param->FRACTIONS[k]>0){
			Brush->N_Chain_Frac[k]=floor(Param->N_Polymer*Param->FRACTIONS[k]);
		}else{
			Brush->N_Chain_Frac[k]=0;
		}
	}
	Param->N_Polymer=0;
	for(k=0;k<POL_END_TYPES;k++){
		Param->N_Polymer+=Brush->N_Chain_Frac[k]; //Correctly summing all the polymer types
	}
	
	printf("After Fraction Param->N_Polymer=%d\n",Param->N_Polymer);fflush(NULL);
	RandomAnchorsType(Param,Brush->N_Anchor,Param->N_Polymer,Brush);


	
	return;
}


void create_brush_anchor_plane (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush){
	int idx,i,k,m;
	double theta;
	double dist=0;
	double dx,dy,dz;
	int flag=0;

	idx=0;
	Brush->Anchors=(struct part*)calloc(Param->N_Polymer,sizeof(struct part));
	for(k=0;k<Param->N_Polymer;k++){
		//Place anchor points on the surface of the cylinder
		theta=ran3(&Param->seed)*2*PI;
		Brush->Anchors[idx].x=ran3(&Param->seed)*(sist->box_x-4*Param->Polymer_Size)+2*Param->Polymer_Size;		
		Brush->Anchors[idx].y=ran3(&Param->seed)*(sist->box_y-4*Param->Polymer_Size)+2*Param->Polymer_Size;
		Brush->Anchors[idx].z=sist->MOL_SHIFT-Param->Polymer_Size*1.1/2;
		Brush->Anchors[idx].flag=-1;
		//Test overalp with Protein
		flag=1;
		switch(Param->Simul_type){
			case TRANS:
				if((Brush->Anchors[idx].x<Brush->minx)||(Brush->Anchors[idx].x>Brush->maxx)) flag=0; // This removes the polymers at the entrnace and exit of the pore to allow the entrance of the protein
			/*if((Param->mask_flag==YES)&&(Param->protein_flag==YES)){
					
					for(m=0;m<Param->N_Proteins;m++){
						for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
							dx=(Brush->Anchors[idx].x-Protein_pointers[m]->Atom[i].x); //Remove the brush before the proteins initial position
							


							if(dx<0) flag=0;

						}
					}
				}*/
				
				
				
			break;
			case ANCHOR:
				if((Param->mask_flag==YES)&&(Param->protein_flag==YES)){
					
					for(m=0;m<Param->N_Proteins;m++){
						for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
							dx=(Protein_pointers[m]->Atom[i].x-Brush->Anchors[idx].x);
							dy=(Protein_pointers[m]->Atom[i].y-Brush->Anchors[idx].y);
							//printf("Protein_pointers[%d]->Atom[%d].x=%lf Protein_pointers[%d]->Atom[%d].y=%lf\n",m,i,Protein_pointers[m]->Atom[i].x,m,i,Protein_pointers[m]->Atom[i].y);fflush(NULL);
							//printf("Brush->Anchors[%d].x=%lf Brush->Anchors[%d].y=%lf\n",idx,Brush->Anchors[idx].x,idx,Brush->Anchors[idx].y);fflush(NULL);

							dist=sqrt(dx*dx+dy*dy);

							
							if(dist<Param->Polymer_Size*2) flag=0;
							//printf("dist=%lf flag=%d\n",dist,flag);fflush(NULL);
						}
					}
				}
			break;		
		}
		if(flag==1){
			Brush->Anchors[idx].lammps_idx=idx*Param->Polymer_Length;	
			Brush->Anchors[idx].flag=-1;
			idx++;
		}
		
	}
	Brush->N_Anchor=Param->N_Polymer=idx;
	
	printf("Before Fraction Param->N_Polymer=%d\n",Param->N_Polymer);fflush(NULL);
	
	for(k=0;k<POL_END_TYPES;k++){
		if(Param->FRACTIONS[k]>0){
			Brush->N_Chain_Frac[k]=floor(Param->N_Polymer*Param->FRACTIONS[k]);
		}else{
			Brush->N_Chain_Frac[k]=0;
		}
	}
	Param->N_Polymer=0;
	for(k=0;k<POL_END_TYPES;k++){
		Param->N_Polymer+=Brush->N_Chain_Frac[k]; //Correctly summing all the polymer types
	}
	
	printf("After Fraction Param->N_Polymer=%d\n",Param->N_Polymer);fflush(NULL);
	RandomAnchorsType(Param,Brush->N_Anchor,Param->N_Polymer,Brush);


	
	return;
}

void selectRandomAnchors(struct Interaction_Param *Param, int N, int Nselect, struct Brushes *Brush) {
	
	if(Nselect>=N) return;
	
    int numbers[N];
    int selected[Nselect];
    int i, j;

    // Initialize the array with numbers from 0 to 9
    for (i = 0; i < N; i++) {
        numbers[i] = i;
    }

  
    // Randomly select 5 distinct numbers
    for (i = 0; i < Nselect; i++) {
        int randIndex = ran3(&Param->seed)*(Nselect - i);

        // Swap the selected random number to the last available position in the array
        int temp = numbers[Nselect - i - 1];
        numbers[Nselect - i - 1] = numbers[randIndex];
        numbers[randIndex] = temp;

        // Store the selected number in the output array
        selected[i] = numbers[Nselect - i - 1];
    }
		for (i = 0; i < Nselect; i++) {
			Brush->Anchors[selected[i]].flag=1;
		}
  	return;
}


void RandomAnchorsType(struct Interaction_Param *Param, int N, int Nselect, struct Brushes *Brush) {
	
	
	
    int numbers[N];
    int selected[N];
    int i, j,k;
		int runningTotal;
		
		
    // Initialize the array with numbers from 0 to 9
    for (i = 0; i < N; i++) {
        numbers[i] = i;
    }
		Brush->indeces=(int *)calloc(Nselect,sizeof(int));
  
    // Randomly select 5 distinct numbers
    for (i = 0; i < Nselect; i++) {
        int randIndex = ran3(&Param->seed)*(Nselect - i);

        // Swap the selected random number to the last available position in the array
        int temp = numbers[Nselect - i - 1];
        numbers[Nselect - i - 1] = numbers[randIndex];
        numbers[randIndex] = temp;

        // Store the selected number in the output array
        selected[i] = numbers[Nselect - i - 1];
    }
		runningTotal = 0;
    for (k = 0; k < POL_END_TYPES; k++) {
        for (i = runningTotal; i < runningTotal + Brush->N_Chain_Frac[k]; i++) {
            Brush->Anchors[selected[i]].flag = k;  // Assuming types are 1-indexed
						Brush->indeces[i]=selected[i];
        }
        runningTotal += Brush->N_Chain_Frac[k];
    }
  	return;
}


void create_brush_anchor_sphere (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush){
	

	int N_points=Param->N_Polymer+1;
	int N_count=0;
	double a=4*PI/N_points; //Average area per point
	double d=sqrt(a); //average distances between points
	int M_theta=lround(PI/d),M_phi;
	double d_theta=PI/M_theta;
	double d_phi=a/d_theta;
	double theta,phi;
	int i=0,j=0;
	int l,k,m;
	double temp_x,temp_y,temp_z;
	double dir_x,dir_y,dir_z,norm;
	char type_chain[4];
	char type_end1[3];
	char type_end2[3];
	struct part *Anchors=NULL;
	double proj_x,proj_y,proj_z;
	int idx=0;
	int flag=0;
	double dx,dy,dz,dist;
	
	
	//Counting number of polymers
	for(i=0;i<M_theta;i++){
		theta=PI*(i+0.5)/M_theta;
		N_count+=lround(2*PI*sin(theta)/d_phi);	
	}
	Brush->N_Anchor=N_count;
	N_count=0;
	
	Anchors=(struct part*)calloc(Brush->N_Anchor,sizeof(struct part));
	//Fille the sphere surface with hard points
	for(i=0;i<M_theta;i++){
		theta=PI*(i+0.5)/M_theta;
		M_phi=lround(2*PI*sin(theta)/d_phi);
		for(j=0;j<M_phi;j++){
			
			phi=2*PI*j/M_phi;
			Anchors[N_count].x=Param->Coll_R*sin(theta)*cos(phi);
			Anchors[N_count].y=Param->Coll_R*sin(theta)*sin(phi);
			Anchors[N_count].z=Param->Coll_R*cos(theta);
			Anchors[N_count].flag=-1;
			N_count++;
		}
		
	}
	
	Brush->Anchors=(struct part*)calloc(Brush->N_Anchor,sizeof(struct part));
	//Test overalp with Protein
	idx=0;
	for(k=0;k<Brush->N_Anchor;k++){
		flag=1;
		switch(Param->Simul_type){

			case ANCHOR:
				if((Param->mask_flag==YES)&&(Param->protein_flag==YES)){
					
					for(m=0;m<Param->N_Proteins;m++){
						for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){

							dir_x=(Protein_pointers[m]->Atom[i].x-sist->Half_Box_x); //Vector Residue-Center of Colloid
							dir_y=(Protein_pointers[m]->Atom[i].y-sist->Half_Box_y);
							dir_z=(Protein_pointers[m]->Atom[i].z-sist->Half_Box_z);

							norm=sqrt(dir_x*dir_x+dir_y*dir_y+dir_z*dir_z);
							dir_x/=norm;
							dir_y/=norm;
							dir_z/=norm;

							proj_x=sist->Half_Box_x+dir_x*(Param->Coll_R); // Projected position of the residue over the colloid surface
							proj_y=sist->Half_Box_y+dir_y*(Param->Coll_R);
							proj_z=sist->Half_Box_z+dir_z*(Param->Coll_R);

							dx=(proj_x-Anchors[k].x-sist->Half_Box_x);
							dy=(proj_y-Anchors[k].y-sist->Half_Box_y);
							dz=(proj_z-Anchors[k].z-sist->Half_Box_z);

							dist=sqrt(dx*dx+dy*dy+dz*dz);


								if((dist<Param->Polymer_Size*2)) flag=0;

							}
					}
				}
			break;		
		}
		if(flag==1){
			Brush->Anchors[idx].x=Anchors[k].x;
			Brush->Anchors[idx].y=Anchors[k].y;
			Brush->Anchors[idx].z=Anchors[k].z;
			Brush->Anchors[idx].lammps_idx=idx*Param->Polymer_Length;	
			Brush->Anchors[idx].flag=-1;
			idx++;
		}

	}

	
	Brush->N_Anchor=idx;
	//randomly select the brsuh grafting points
	if(Param->N_Polymer>Brush->N_Anchor) Param->N_Polymer=Brush->N_Anchor;

	
	for(k=0;k<POL_END_TYPES;k++){
		if(Param->FRACTIONS[k]>0){
			Brush->N_Chain_Frac[k]=floor(Param->N_Polymer*Param->FRACTIONS[k]);
		}else{
			Brush->N_Chain_Frac[k]=0;
		}
	}
	Param->N_Polymer=0;
	for(k=0;k<POL_END_TYPES;k++){
		Param->N_Polymer+=Brush->N_Chain_Frac[k]; //Correctly summing all the polymer types
	}
	
	RandomAnchorsType(Param,Brush->N_Anchor,Param->N_Polymer,Brush);
	
	
	
	printf("create_brush_anchor_sphere Param->N_Polymer=%d\n",Param->N_Polymer);fflush(NULL);
	
	
	
	
	return;
}



void write_topology (struct GlobalSistem *sist,struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush,FILE *f_atoms, double scale, int type){
	int idx=0;
	
	
	
	
	Init_Topology(sist,Param,Protein_pointers[0],Brush,type);
	
	fprintf(f_atoms,"LAMMPS Description\n\n");
	fprintf(f_atoms,"%d atoms\n",sist->Tot_Atoms); 
	fprintf(f_atoms,"%d bonds\n",sist->Tot_Bonds);
	fprintf(f_atoms,"%d angles\n",sist->Tot_Angles);
	fprintf(f_atoms,"%d dihedrals\n",sist->Tot_Dihedrals);
	fprintf(f_atoms,"00 impropers\n\n");	
	fprintf(f_atoms,"%d atom types\n",sist->tot_atom_type); // Protein LInker Polymer Brush and the srd particles
	fprintf(f_atoms,"%d bond types\n",sist->tot_bonds_type); // Inclides the bond type that connects the linker moelcules to the protein 
	fprintf(f_atoms,"%d angle types\n",sist->Tot_Angles_Types);
	fprintf(f_atoms,"%d dihedral types\n",sist->Tot_Dihedrals_Types);
	fprintf(f_atoms,"00 improper types\n\n");
	
	
	
	
	fprintf(f_atoms,"0 %lf xlo xhi\n0 %lf ylo yhi\n0 %lf zlo zhi\n\n",sist->box_x*scale,sist->box_y*scale,sist->box_z*scale);
	//modified to make sure masses are set for all atom types in GO model
	fprintf(f_atoms,"\n\n");
	
	fprintf(f_atoms,"Masses\n\n");
	write_masses(sist,f_atoms,Param,Protein_pointers[0]);
	fprintf(f_atoms,"\n\n");
	
	
	fprintf(f_atoms,"Atoms\n\n");
	sist->N_Atoms=1;
	sist->N_Mol=1;
	if(Param->protein_flag==YES){
		write_protein_atoms(sist,Param,f_atoms,Protein_pointers,scale);
		if(Param->mask_flag==YES){
			write_linker_atoms(sist,f_atoms,Param,Protein_pointers,scale);				
		}
	}
	
	
	
	if((type==ALL)||(type==NO_POL_ANCHORS)) write_brush_atoms(sist,f_atoms,Param,Protein_pointers,Brush,scale,type); 
	
	fprintf(f_atoms,"\n\n");
	
	fprintf(f_atoms,"Bonds\n\n");
	sist->N_Bonds=1;
	sist->bondtype=1;
	if(Param->protein_flag==YES) write_protein_bonds(sist,f_atoms, Param, Protein_pointers);
	sist->polymers_bonds_type=sist->bondtype;
	if((Param->mask_flag==YES)&&(Param->protein_flag==YES)) write_linker_bonds(sist,f_atoms,Param,Protein_pointers);				
	
	
	if((type==ALL)||(type==NO_POL_ANCHORS)) write_brush_bonds(sist,f_atoms,Param,Brush,type);
	fprintf(f_atoms,"\n\n");
	if(Param->protein_flag==YES){
		fprintf(f_atoms,"Angles\n\n");
		write_protein_angles(sist,f_atoms,Param,Protein_pointers);
		fprintf(f_atoms,"\n\n");


		fprintf(f_atoms,"Dihedrals\n\n");
		write_protein_dihedrals(sist,f_atoms,Param,Protein_pointers);
		fprintf(f_atoms,"\n\n");
	}
	
	return;
}


void write_masses (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins *Protein){
	int idx=0;
	int i;
	
		
	
		for(i=1;i<=PROTEIN;i++){
			fprintf(f_atoms,"%d 1.0\n",i);

		}
	fprintf(f_atoms,"%d 100.0\n",COLLOID        );
	fprintf(f_atoms,"%d 1.0\n",LINKER        );
	fprintf(f_atoms,"%d 1.0\n",ANCHOR_LINKER );
	fprintf(f_atoms,"%d 1.0\n",POLYMER       );
	fprintf(f_atoms,"%d 1.0\n",ANCHOR_POLYMER);
	fprintf(f_atoms,"%d 1.0\n",SRD           );
	for(i=END_POLYMER;i<END_POLYMER+POL_END_TYPES;i++){
		fprintf(f_atoms,"%d 1.0\n", i);
	}
	return;
	
}	

void deep_copy_proteins(struct Proteins *src,struct Proteins *dest) {
    
		int i;
		
		
    // Copy simple fields
    dest->N = src->N;
    dest->Nres = src->Nres;
    dest->NMasked = src->NMasked;
    dest->NSkip = src->NSkip;
    dest->Nread = src->Nread;
    dest->max_dist = src->max_dist;
    dest->FIRST_AtomHx = src->FIRST_AtomHx;
    dest->FIRST_AtomHy = src->FIRST_AtomHy;
    dest->FIRST_AtomHz = src->FIRST_AtomHz;
    dest->bondtype_head = src->bondtype_head;

    // Copy arrays of basic types
   	strcpy(dest->Atoms_name[ATOM_N],"N  ");
		strcpy(dest->Atoms_name[ATOM_CA],"CA ");
		strcpy(dest->Atoms_name[ATOM_C],"C  ");
		strcpy(dest->Atoms_name[ATOM_O],"O  ");
		strcpy(dest->Atoms_name[ATOM_H],"H  ");


		strcpy(dest->Amminoacids[1],"ALA");
		strcpy(dest->Amminoacids[2],"CYS");
		strcpy(dest->Amminoacids[3],"ASP");
		strcpy(dest->Amminoacids[4],"GLU");
		strcpy(dest->Amminoacids[5],"PHE");
		strcpy(dest->Amminoacids[6],"GLY");
		strcpy(dest->Amminoacids[7],"HIS");
		strcpy(dest->Amminoacids[8],"ILE");
		strcpy(dest->Amminoacids[9],"LYS");
		strcpy(dest->Amminoacids[10],"LEU");
		strcpy(dest->Amminoacids[11],"MET");
		strcpy(dest->Amminoacids[12],"ASN");
		strcpy(dest->Amminoacids[13],"PRO");
		strcpy(dest->Amminoacids[14],"GLN");
		strcpy(dest->Amminoacids[15],"ARG");
		strcpy(dest->Amminoacids[16],"SER");
		strcpy(dest->Amminoacids[17],"THR");
		strcpy(dest->Amminoacids[18],"VAL");
		strcpy(dest->Amminoacids[19],"TRP");
		strcpy(dest->Amminoacids[20],"TYR");
		

    // Deep copy the arrays of struct part and double
		dest->Linker = (struct part *) calloc(src->N ,sizeof(struct part));
    dest->Atom = (struct part *) calloc(src->N ,sizeof(struct part));
    dest->Atom_read = (struct part *) calloc(src->Nread ,sizeof(struct part));
    
    dest->Linker_Anchors = (struct part *) calloc(src->N ,sizeof(struct part));
    dest->gmass = (struct part *) calloc(1 ,sizeof(struct part));
    dest->Average_r = (double *) calloc(src->Nres ,sizeof(double));
    dest->Count_r2 = (double *) calloc(src->Nres ,sizeof(double));
    dest->Ignore = (int *) calloc(src->N ,sizeof(int));
    dest->Mask = (int *) calloc(src->N ,sizeof(int));
    dest->Mask_wall = (int *) calloc(src->N ,sizeof(int));
		dest->group_id = (char *) calloc(src->N*10 ,sizeof(char));
		
    if (!dest->Atom || !dest->Atom_read || !dest->Linker || !dest->Linker_Anchors ||
        !dest->gmass || !dest->Average_r || !dest->Count_r2 ||
        !dest->Ignore || !dest->Mask || !dest->Mask_wall) {
        // Handle memory allocation failure, e.g., free allocated memory and exit
    }
		copy_part (src->Atom,dest->Atom, src->N);
		copy_part (src->Atom_read,dest->Atom_read, src->Nread);
		copy_part (src->Linker,dest->Linker, src->N);
		copy_part (src->Linker_Anchors,dest->Linker_Anchors, src->N);
		copy_part (src->gmass,dest->gmass, 1);
    
		for(i=0;i<src->Nres;i++){
			dest->Average_r[i]=src->Average_r[i];
			dest->Count_r2[i]=src->Count_r2[i];
		}
		for(i=0;i<src->N;i++){
			dest->Ignore[i]=src->Ignore[i];
			dest->Mask[i]=src->Mask[i];
			dest->Mask_wall[i]=src->Mask_wall[i];		
		}
   strcpy(dest->group_id,src->group_id);
    
		
    return;
}
void copy_part (struct part *src,struct part *dest, int N){
	int i;
	
	for(i=0;i<N;i++){
		dest[i].x= 					src[i].x;
		dest[i].y= 					src[i].y;
		dest[i].z= 					src[i].z;
		dest[i].id=					src[i].id;
		dest[i].flag=				src[i].flag;
		dest[i].res= 				src[i].res;
		dest[i].lammps_idx=	src[i].lammps_idx;
		dest[i].prot=				src[i].prot;		
	}
	
	return;
}

	
void write_protein_atoms (struct GlobalSistem *sist, struct Interaction_Param *Param, FILE *f_atoms, struct Proteins **Protein_pointers, double scale){	
	int i=1;
	int ii=1;
	int k;
	static struct part *Sphere_Points=NULL;
	double x,y,z;
	double MOL_SHIFT=0;
	char temp_string[10];
	double min_z_prot=10000,max_z_prot=-1000;

	
		switch(Param->Simul_type){
				case TRANS: // This includes the interaction of a free protein with a NP					
					for(k=0;k<Param->N_Proteins;k++){						
						ii=1;
							for(i=ATOM_CA;i<Protein_pointers[k]->N;i+=NATOM){
								x=Protein_pointers[k]->Atom[i].x; //Shifts the protein at the entrance of the pore
								y=Protein_pointers[k]->Atom[i].y;
								z=Protein_pointers[k]->Atom[i].z;	
								Protein_pointers[k]->Atom[i].lammps_idx=sist->N_Atoms;
								/***********Format Atom_Id Molecule_Id Atom_Type Charge Coordinates******************/
								fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",Protein_pointers[k]->Atom[i].lammps_idx,sist->N_Mol,ii,x*scale,y*scale,z*scale);
								sprintf(temp_string, "%d ", sist->N_Atoms);
       		 			strcat(Protein_pointers[k]->group_id, temp_string);
								sist->N_Atoms++;
								ii++;
							}														
					}		
				break;
				case ANCHOR:	 // This only makes sense for a single protein fixated on a surface. At the momnet only a cylinder						
					switch(Param->Geometry){
						case SPHERE:
							for(k=0;k<Param->N_Proteins;k++){
								ii=1;
								for(i=ATOM_CA;i<Protein_pointers[k]->N;i+=NATOM){	
									x=Protein_pointers[k]->Atom[i].x; // Shifts the protein around the colloid in the middle of the box
									y=Protein_pointers[k]->Atom[i].y;
									z=Protein_pointers[k]->Atom[i].z;
									Protein_pointers[k]->Atom[i].lammps_idx=sist->N_Atoms;
									/***********Format Atom_Id Molecule_Id Atom_Type Charge Coordinates******************/
									fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",Protein_pointers[k]->Atom[i].lammps_idx,sist->N_Mol,ii,x*scale,y*scale,z*scale);
									sprintf(temp_string, "%d ", sist->N_Atoms);
									strcat(Protein_pointers[k]->group_id, temp_string);
									sist->N_Atoms++;
									ii++;																																			
								}	
							}
						break;	
						case CYLINDER:
							for(i=ATOM_CA;i<Protein_pointers[0]->N;i+=NATOM){
								if(min_z_prot>Protein_pointers[0]->Atom[i].z) min_z_prot=Protein_pointers[0]->Atom[i].z;
								if(max_z_prot<Protein_pointers[0]->Atom[i].z) max_z_prot=Protein_pointers[0]->Atom[i].z;
							}

							z = min_z_prot - (sist->radius - sist->inradius) - sist->MOL_SHIFT - Param->Linker_Size * 1.1; // We first shift the lowest CA to 0 then we add the radius-inradius gap, then the linker length and then Param->Linker_Size*1.1
							for(i=0;i<Protein_pointers[0]->N;i++){
								Protein_pointers[0]->Atom[i].z -= z;
							}

							for(i=ATOM_CA;i<Protein_pointers[0]->N;i+=NATOM){	
								x = Protein_pointers[0]->Atom[i].x; // Shifts the protein around the colloid in the middle of the box
								y = Protein_pointers[0]->Atom[i].y;
								z = Protein_pointers[0]->Atom[i].z;
								Protein_pointers[0]->Atom[i].lammps_idx = sist->N_Atoms;
								/***********Format Atom_Id Molecule_Id Atom_Type Charge Coordinates******************/
								fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",Protein_pointers[0]->Atom[i].lammps_idx,sist->N_Mol,ii,x*scale,y*scale,z*scale);
								sprintf(temp_string, "%d ", sist->N_Atoms);
								strcat(Protein_pointers[0]->group_id, temp_string);
								sist->N_Atoms++;
								ii++;
							}	
						break;
						case SLIT:
						case CHANNEL:
							for(i=ATOM_CA;i<Protein_pointers[0]->N;i+=NATOM){
								if(min_z_prot>Protein_pointers[0]->Atom[i].z) min_z_prot=Protein_pointers[0]->Atom[i].z;
								if(max_z_prot<Protein_pointers[0]->Atom[i].z) max_z_prot=Protein_pointers[0]->Atom[i].z;
							}

							z = min_z_prot -2*sist->MOL_SHIFT - Param->Linker_Size * 1.1; // We first shift the lowest CA to 0 then we add the radius-inradius gap, then the linker length and then Param->Linker_Size*1.1
							for(i=0;i<Protein_pointers[0]->N;i++){
								Protein_pointers[0]->Atom[i].z -= z;
							}

							for(i=ATOM_CA;i<Protein_pointers[0]->N;i+=NATOM){	
								x = Protein_pointers[0]->Atom[i].x; // Shifts the protein around the colloid in the middle of the box
								y = Protein_pointers[0]->Atom[i].y;
								z = Protein_pointers[0]->Atom[i].z;
								Protein_pointers[0]->Atom[i].lammps_idx = sist->N_Atoms;
								/***********Format Atom_Id Molecule_Id Atom_Type Charge Coordinates******************/
								fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",Protein_pointers[0]->Atom[i].lammps_idx,sist->N_Mol,ii,x*scale,y*scale,z*scale);
								sprintf(temp_string, "%d ", sist->N_Atoms);
								strcat(Protein_pointers[0]->group_id, temp_string);
								sist->N_Atoms++;
								ii++;
							}	
						break;
					}												
				break;				
		}
		sist->N_Mol++;
	
	return;	
}



void create_protein_copies (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers){	
	int i=1;
	int ii=1;
	int k;
	static struct part *Sphere_Points=NULL;
	double x,y,z;
	double MOL_SHIFT=0;
	char temp_string[10];
	double new_x,new_y,new_z;
	int flag=0;

	MOL_SHIFT=(Param->Coll_R+Protein_pointers[0]->max_dist/2.+1.5);
	//MOL_SHIFT=Param->Coll_R;
	if(Sphere_Points==NULL){
		Sphere_Points=(struct part *)calloc(1000,sizeof(struct part));
		Sphere_Points[0].x=0;
 		Sphere_Points[0].y=0;
 		Sphere_Points[0].z=1;
		if(Param->N_Proteins>1) proteins_equidistant(Param,Sphere_Points);
	}

	
		switch(Param->Simul_type){
			case TRANS: // This includes the interaction of a free protein with a N
					for(k=0;k<Param->N_Proteins;k++){
							if(k>0) {
								Protein_pointers[k]=(struct Proteins*)calloc(1,sizeof(struct Proteins));
								deep_copy_proteins(Protein_pointers[0],Protein_pointers[k]);
							}
							protein_rotate_shift(sist,Param,Protein_pointers[k]);
					}
					for(k=0;k<Param->N_Proteins;k++){
						for(i=ATOM_CA;i<Protein_pointers[k]->N;i+=NATOM){	
							switch(Param->Geometry){
							case SPHERE:
								if(Param->N_Colloids==2){
									Protein_pointers[k]->Atom[i].x=(Protein_pointers[k]->Atom[i].x); //Shifts the protein at the entrance of the pore created by the colloids
									Protein_pointers[k]->Atom[i].y=(Protein_pointers[k]->Atom[i].y);
									Protein_pointers[k]->Atom[i].z=(Protein_pointers[k]->Atom[i].z+sist->Half_Box_z);
								}
								if(Param->N_Colloids==1){
									Protein_pointers[k]->Atom[i].x=(Protein_pointers[k]->Atom[i].x); //Shifts the protein to a corner
									Protein_pointers[k]->Atom[i].y=(Protein_pointers[k]->Atom[i].y);
									Protein_pointers[k]->Atom[i].z=(Protein_pointers[k]->Atom[i].z);
								}
							break;
							case CYLINDER:
								Protein_pointers[k]->Atom[i].x=(Protein_pointers[k]->Atom[i].x); //Shifts the protein at the entrance of the pore or in the middle of the box
								Protein_pointers[k]->Atom[i].y=(Protein_pointers[k]->Atom[i].y+sist->Half_Box_y);
								Protein_pointers[k]->Atom[i].z=(Protein_pointers[k]->Atom[i].z+sist->Half_Box_z);
								break;							
							case SLIT:
								Protein_pointers[k]->Atom[i].x=(Protein_pointers[k]->Atom[i].x); //Shifts the protein at the entrance of the pore or in the middle of the box
								Protein_pointers[k]->Atom[i].y=(Protein_pointers[k]->Atom[i].y+sist->Half_Box_y);
								Protein_pointers[k]->Atom[i].z=(Protein_pointers[k]->Atom[i].z+sist->Half_Box_z);
								break;
							case CHANNEL:
								Protein_pointers[k]->Atom[i].x=(Protein_pointers[k]->Atom[i].x); //Shifts the protein at the entrance of the pore or in the middle of the box
								Protein_pointers[k]->Atom[i].y=(Protein_pointers[k]->Atom[i].y+sist->Half_Box_y);
								Protein_pointers[k]->Atom[i].z=(Protein_pointers[k]->Atom[i].z+sist->Half_Box_z);
								break;
							}								
						}
							
					}		
			break;
			case ANCHOR:
				switch (Param->Geometry)
				{
				case SPHERE:
					for (k = 0; k < Param->N_Proteins; k++)
					{
						if (k > 0)
						{
							Protein_pointers[k] = (struct Proteins *)calloc(1, sizeof(struct Proteins));
							deep_copy_proteins(Protein_pointers[0], Protein_pointers[k]);
						}
						protein_rotate_shift(sist, Param, Protein_pointers[k]);
					}
					for (k = 0; k < Param->N_Proteins; k++)
					{
						for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
						{

							Protein_pointers[k]->Atom[i].x = (Protein_pointers[k]->Atom[i].x + Sphere_Points[k].x * MOL_SHIFT + sist->Half_Box_x); // Shifts the protein around the colloid in the middle of the box
							Protein_pointers[k]->Atom[i].y = (Protein_pointers[k]->Atom[i].y + Sphere_Points[k].y * MOL_SHIFT + sist->Half_Box_y);
							Protein_pointers[k]->Atom[i].z = (Protein_pointers[k]->Atom[i].z + Sphere_Points[k].z * MOL_SHIFT + sist->Half_Box_z);
						}
					}
					break;
				case CYLINDER:
					k = 0;
					while (k < Param->N_Proteins)
					{
						double offset_x = (ran3(&Param->seed))*sist->box_x + sist->Half_Box_x;
						double offset_y = (ran3(&Param->seed))*sist->box_y + sist->Half_Box_y;
						flag = 0;
						for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
						{
							new_x = (Protein_pointers[k]->Atom[i].x + offset_x);
							new_y = (Protein_pointers[k]->Atom[i].y + offset_y);
							new_z = (Protein_pointers[k]->Atom[i].z);
							if (sqrt(new_y * new_y + new_z * new_z) > sist->inradius)
							{
								flag = 1;
							}
						}
						if (flag == 0)
						{
							for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
							{
								new_x = (Protein_pointers[k]->Atom[i].x + offset_x);
								new_y = (Protein_pointers[k]->Atom[i].y + offset_y);
								new_z = (Protein_pointers[k]->Atom[i].z);
								Protein_pointers[k]->Atom[i].x = new_x; // Shifts the protein in the middle of the box
								Protein_pointers[k]->Atom[i].y = new_y;
								Protein_pointers[k]->Atom[i].z = new_z;
							}
							k++;
						}
					}
					break;
				case SLIT:
					for (k = 0; k < Param->N_Proteins; k++)
					{
						double offset_x = (ran3(&Param->seed))*sist->box_x + sist->Half_Box_x;
						double offset_y = (ran3(&Param->seed))*sist->box_y + sist->Half_Box_y;
						for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
						{
							Protein_pointers[k]->Atom[i].x = (Protein_pointers[k]->Atom[i].x + offset_x); // Shifts the protein in the middle of the box
							Protein_pointers[k]->Atom[i].y = (Protein_pointers[k]->Atom[i].y + offset_y);
							Protein_pointers[k]->Atom[i].z = (Protein_pointers[k]->Atom[i].z);
						}
					}
					break;
				case CHANNEL:
					k = 0;
					while (k < Param->N_Proteins)
					{
						double offset_x = (ran3(&Param->seed))*sist->box_x + sist->Half_Box_x;
						double offset_y = (ran3(&Param->seed))*sist->box_y + sist->Half_Box_y;
						flag = 0;
						for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
						{
							new_x = (Protein_pointers[k]->Atom[i].x + offset_x);
							new_y = (Protein_pointers[k]->Atom[i].y + offset_y);
							new_z = (Protein_pointers[k]->Atom[i].z);
							if (new_y > sist->box_y-sist->MOL_SHIFT || new_y < sist->MOL_SHIFT )
							{
								flag = 1;
							}
						}
						if (flag == 0)
						{
							for (i = ATOM_CA; i < Protein_pointers[k]->N; i += NATOM)
							{
								new_x = (Protein_pointers[k]->Atom[i].x + offset_x);
								new_y = (Protein_pointers[k]->Atom[i].y + offset_y);
								new_z = (Protein_pointers[k]->Atom[i].z);
								Protein_pointers[k]->Atom[i].x = new_x; // Shifts the protein in the middle of the box
								Protein_pointers[k]->Atom[i].y = new_y;
								Protein_pointers[k]->Atom[i].z = new_z;
							}
							k++;
						}
					}
					break;
				}
				break;
		}
		sist->N_Mol++;
	
	return;	
}

void proteins_equidistant(struct Interaction_Param *Param, struct part *Sphere_Points){
	//distribute the proteins uniformly around the sphere. Taken from Deserno2004
	int N_points=Param->N_Proteins;
	double r=1.;
	double r2=r*r;
	int N_count=0;
	double a=4*PI*r2/N_points;
	double d=sqrt(a);
	int M_theta=lround(PI/d),M_phi;
	double d_theta=PI/M_theta;
	double d_phi=a/d_theta;
	double theta,phi;
	int i=0,j=0;
	int l,k;
	double temp_x,temp_y,temp_z;
	
	
	
	
	for(i=0;i<M_theta;i++){
		theta=PI*(i+0.5)/M_theta;
		M_phi=lround(2*PI*sin(theta)/d_phi);
		for(j=0;j<M_phi;j++){
			
			phi=2*PI*j/M_phi;
			Sphere_Points[N_count].x=r*sin(theta)*cos(phi);
			Sphere_Points[N_count].y=r*sin(theta)*sin(phi);
			Sphere_Points[N_count].z=r*cos(theta);
			N_count++;
		}
		
	}
	
	if(Param->N_Proteins>N_count) Param->N_Proteins=N_count;
	
	
	return;
}


void write_linker_atoms (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, double scale){
	double module;
	double translatex,translatey,translatez;
	int i,k,m;
	int Linker;
	char temp_string[10];
	
	
	for(m=0;m<Param->N_Proteins;m++){
		/******************* LINKER CONTRUCTION *****************/
		for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
			if(Protein_pointers[m]->Mask[i]==1) {	
				Protein_pointers[m]->Linker_Anchors[i].lammps_idx=Protein_pointers[m]->Atom[i].lammps_idx;
				Protein_pointers[m]->Linker[i].lammps_idx=sist->N_Atoms;

				//Growing the linker from the protein anchor
				for(k=1;k<=Param->Linker_Length-1;k++){

						switch(Param->Geometry){
							case SPHERE:
								translatex=sist->Half_Box_x-Protein_pointers[m]->Atom[i].x;
								translatey=sist->Half_Box_y-Protein_pointers[m]->Atom[i].y;
								translatez=sist->Half_Box_z-Protein_pointers[m]->Atom[i].z;
							break;
							case CYLINDER:										
								translatex=0;
								translatey=0;
								translatez=-Protein_pointers[m]->Atom[i].z;
							break;
							case SLIT:										
								translatex=0;
								translatey=0;
								translatez=-Protein_pointers[m]->Atom[i].z;
							break;
							case CHANNEL:										
								translatex=0;
								translatey=0;
								translatez=-Protein_pointers[m]->Atom[i].z;
							break;
						}
					module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

					translatex=(translatex*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].x;
					translatey=(translatey*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].y;
					translatez=(translatez*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].z;

					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,LINKER,translatex*scale,translatey*scale,translatez*scale);
					
					
					sprintf(temp_string, "%d ", sist->N_Atoms);
       		strcat(Protein_pointers[m]->group_id, temp_string);
					sist->N_Atoms++;
				}
				//Adding the dangling anchor that will stick to the surface covalently
				k=Param->Linker_Length;
				switch(Param->Geometry){
					case SPHERE:
						translatex=sist->Half_Box_x-Protein_pointers[m]->Atom[i].x;
						translatey=sist->Half_Box_y-Protein_pointers[m]->Atom[i].y;
						translatez=sist->Half_Box_z-Protein_pointers[m]->Atom[i].z;
					break;
					case CYLINDER:										
						translatex=0;
						translatey=0;
						translatez=-Protein_pointers[m]->Atom[i].z;
					break;
					case SLIT:										
						translatex=0;
						translatey=0;
						translatez=-Protein_pointers[m]->Atom[i].z;
					break;
					case CHANNEL:										
						translatex=0;
						translatey=0;
						translatez=-Protein_pointers[m]->Atom[i].z;
					break;
				}

				module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

				translatex=(translatex*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].x;
				translatey=(translatey*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].y;
				translatez=(translatez*Param->Linker_Size*k)/module+Protein_pointers[m]->Atom[i].z;

				fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,ANCHOR_LINKER,translatex*scale,translatey*scale,translatez*scale);
				sprintf(temp_string, "%d ", sist->N_Atoms);
       	strcat(Protein_pointers[m]->group_id, temp_string);
				sist->N_Atoms++;

				sist->N_Mol++;

			}
		}		
	}
	return;
	
}
void write_brush_atoms (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type){
	int k;
	
	
	
	switch(Param->Geometry){
			case CYLINDER:
				grow_cylinder_brush(sist,f_atoms,Param,Protein_pointers,Brush,scale,type);					
			break;
			case SLIT:
				grow_plane_brush(sist,f_atoms,Param,Protein_pointers,Brush,scale,type);					
			break;
			case CHANNEL:
				grow_plane_brush(sist,f_atoms,Param,Protein_pointers,Brush,scale,type);					
			break;
			case SPHERE:						
				grow_sphere_brush(sist,f_atoms,Param,Protein_pointers,Brush,scale,type);													
			break;				
	}
	
	
	
	return;	
}


void grow_cylinder_brush (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type){
	
	int i,k;
	int end_type;
	double module;
	double new_x,new_y,new_z;
	double direction_x,direction_y,direction_z;
	struct part *anchor=NULL;
	
	
	for(k=0;k<Param->N_Polymer;k++){
		if(Brush->Anchors[k].flag>-1){
			anchor=&Brush->Anchors[k];
			Brush->Anchors[k].lammps_idx=sist->N_Atoms;
			end_type=Brush->Anchors[k].flag;
			
			if(type==ALL) {  //Includes Anchors
				fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,ANCHOR_POLYMER,anchor->x*scale,(anchor->y+sist->Half_Box_y)*scale,(anchor->z+sist->Half_Box_z)*scale);
				sist->N_Atoms++;
			}
		
			for(i=1;i<Param->Polymer_Length;i++){
				direction_x = 0.0;
				direction_y = -anchor->y;
				direction_z = -anchor->z;

				module=sqrt(direction_x*direction_x+direction_y*direction_y+direction_z*direction_z);

				direction_x /= module;
				direction_y /= module;
				direction_z /= module;
				new_x = anchor->x ;
				if(i==1){
					new_y = anchor->y + direction_y*Param->Polymer_Size*i*2.1;
					new_z = anchor->z + direction_z*Param->Polymer_Size*i*2.1;
				}else{
					new_y = anchor->y + direction_y*Param->Polymer_Size*2;
					new_z = anchor->z + direction_z*Param->Polymer_Size*2;
				}
				if(i==Param->Polymer_Length-1){
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,END_POLYMER+end_type,new_x*scale,(new_y+sist->Half_Box_y)*scale,(new_z+sist->Half_Box_z)*scale);
				}else{
					
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,POLYMER,new_x*scale,(new_y+sist->Half_Box_y)*scale,(new_z+sist->Half_Box_z)*scale);
				
				}
				sist->N_Atoms++;

			}
			sist->N_Mol++;
		
		}	
	}
	return;
	
}

void grow_plane_brush (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type){
	
	int i,k;
	int end_type;
	double module;
	double new_x,new_y,new_z;
	double direction_x,direction_y,direction_z;
	struct part *anchor=NULL;
	
	
	for(k=0;k<Param->N_Polymer;k++){
		if(Brush->Anchors[k].flag>-1){
			anchor=&Brush->Anchors[k];
			Brush->Anchors[k].lammps_idx=sist->N_Atoms;
			end_type=Brush->Anchors[k].flag;
			
			if(type==ALL) {  //Includes Anchors
				fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,ANCHOR_POLYMER,anchor->x*scale,anchor->y*scale,anchor->z*scale);
				sist->N_Atoms++;
			}
		
			for(i=1;i<Param->Polymer_Length;i++){
				direction_x = 0.0;
				direction_y = 0.0;
				direction_z = anchor->z;

				module=sqrt(direction_x*direction_x+direction_y*direction_y+direction_z*direction_z);

				direction_x /= module;
				direction_y /= module;
				direction_z /= module;
				new_x = anchor->x ;
				if(i==1){
					new_y = anchor->y + direction_y*Param->Polymer_Size*i*2.1;
					new_z = anchor->z + direction_z*Param->Polymer_Size*i*2.1;
				}else{
					new_y = anchor->y + direction_y*Param->Polymer_Size*2;
					new_z = anchor->z + direction_z*Param->Polymer_Size*2;
				}
				if(i==Param->Polymer_Length-1){
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,END_POLYMER+end_type,new_x*scale,new_y*scale,new_z*scale);
				}else{
					
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,POLYMER,new_x*scale,new_y*scale,new_z*scale);
				
				}
				sist->N_Atoms++;

			}
			sist->N_Mol++;
		
		}	
	}
	return;
	
}

void grow_sphere_brush (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush, double scale, int type){
	
	int i,k,l;
	int end_type;
	double module;
	double new_x,new_y,new_z;
	double direction_x,direction_y,direction_z;
	struct part *anchor=NULL;
	struct part Colloid;
	
	
	for(l=0;l<Param->N_Colloids;l++){
		if(l==0){
			Colloid.x=sist->Half_Box_x;
			Colloid.y=sist->Half_Box_y;
			Colloid.z=sist->Half_Box_z;
		}else{
			Colloid.x=0;
			Colloid.y=0;
			Colloid.z=0;
		}
		
		
		fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,COLLOID,Colloid.x*scale,Colloid.y*scale,Colloid.z*scale);
		
		sist->N_Atoms++;
		sist->N_Mol++;	
		for(k=0;k<Param->N_Polymer;k++){
			if(Brush->Anchors[k].flag>-1){
				anchor=&Brush->Anchors[k];
				if(l==0) Brush->Anchors[k].lammps_idx=sist->N_Atoms; //counting only for the first colloid
				end_type=Brush->Anchors[k].flag;

				if(type==ALL) {  //Includes Anchors
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,ANCHOR_POLYMER,(anchor->x+Colloid.x)*scale,(anchor->y+Colloid.y)*scale,(anchor->z+Colloid.z)*scale);
					sist->N_Atoms++;
				}

				for(i=1;i<Param->Polymer_Length;i++){
					direction_x = anchor->x;
					direction_y = anchor->y;
					direction_z = anchor->z;

					module=sqrt(direction_x*direction_x+direction_y*direction_y+direction_z*direction_z);

					direction_x /= module;
					direction_y /= module;
					direction_z /= module;
					/*new_x = anchor->x + direction_x*Param->Polymer_Size*i+Colloid.x;
					new_y = anchor->y + direction_y*Param->Polymer_Size*i+Colloid.y;
					new_z = anchor->z + direction_z*Param->Polymer_Size*i+Colloid.z;*/

					new_x = anchor->x + direction_x*Param->Polymer_Size+Colloid.x; //We create a compressed brush so that it naturally expands instead of waiting until it contracts
					new_y = anchor->y + direction_y*Param->Polymer_Size+Colloid.y;
					new_z = anchor->z + direction_z*Param->Polymer_Size+Colloid.z;
					


					if(i==Param->Polymer_Length-1){
						fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,END_POLYMER+end_type,new_x*scale,new_y*scale,new_z*scale);
					}else{
						fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,POLYMER,new_x*scale,new_y*scale,new_z*scale);
					}
					sist->N_Atoms++;

				}
				sist->N_Mol++;

			}	
		}
	}
	
	return;
	
}





/* void write_brush_atoms_martini (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins *Protein, struct Brushes *Brush, double scale, int type){
	int k;
	
	double anchor_x,anchor_y,anchor_z;
	
	
	
	
	for(k=0;k<Param->N_Polymer;k++){
		if(Brush->Anchors[k].flag==1){
			//Recall anchor points on the surface of the cylinder		
			anchor_x=Brush->Anchors[k].x;
			anchor_y=Brush->Anchors[k].y;
			anchor_z=Brush->Anchors[k].z;

			Brush->Anchors[k].lammps_idx=sist->N_Atoms;
			//Shifting the anchors of half box becase are generated on the surface of a cylinder centered in 0
			if(type==ALL) {  //Includes Anchors
					fprintf(f_atoms,"%d %d %d 0 %lf %lf %lf\n",sist->N_Atoms,sist->N_Mol,ANCHOR_POLYMER,anchor_x*scale,(anchor_y+sist->Half_Box_y)*scale,(anchor_z+sist->Half_Box_z)*scale);
				sist->N_Atoms++;
			}
			switch(Param->Geometry){
					case CYLINDER:
						grow_cylinder_brush(sist,f_atoms,Param,Protein,&Brush->Anchors[k],scale,type);					
					break;
					case SPHERE:						
						grow_sphere_brush(sist,f_atoms,Param,Protein,&Brush->Anchors[k],scale, type);													
					break;				
			}
			sist->N_Mol++;
		
		}
		
	}
	return;	
}

void grow_sphere_brush_Martini (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins *Protein, struct part *anchor, double scale){
	
sprintf(line,"%s.gro1",Param->prefix);
	fp_martini=fopen(line,"w");
	sprintf(line,"%s.gro3",Param->prefix);
	fp_martini2=fopen(line,"w");
	fprintf(fp_martini,"GROningen MAchine for Chemical Simulation\n");
	int tot_chain_gromacs=Param->N_Polymer*Param->N_Chains+1;
	for(k=4;k<MOL_ATOM_TYPES;k++){//Adding charged groups ends 
			tot_chain_gromacs+=2*sist->N_Chain_Frac[k];
	}
	int tot_gromacs=tot_chain_gromacs;
	for(k=7;k<MOL_ATOM_TYPES;k++){ //adding extra  hydrogen to acid groups ends from the count
		tot_gromacs+=sist->N_Chain_Frac[k];
	}
	fprintf(fp_martini,"%d\n",tot_gromacs-1+Param->Taken_Fillers);
	
	//fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n","1COR","AU",1,sist->Colloid.x/10,sist->Colloid.y/10,sist->Colloid.z/10);	
	int index2=1;
	chain_index=0;
	printf("LAMMPS  Growth of the brush chains\n");fflush(NULL);
	for(k=0;k<MOL_ATOM_TYPES;k++){ //Loop over the end types
		
		if(k==0||k==1||k==2||k==3) { //Neutral group 
			sprintf(type_chain,"PEN");
			sprintf(type_end1,"XX");
			sprintf(type_end2,"XX");
		} 
		if(k==4||k==5||k==6) {//Basic group 
			sprintf(type_chain,"PEB");			
			sprintf(type_end1,"D");
			sprintf(type_end2,"DP");
		}
		if(k==7||k==8||k==9) { //Acid group 
			sprintf(type_chain,"PEA");
			sprintf(type_end1,"D");
			sprintf(type_end2,"DP");
		}
		printf("LAMMPS  Growth type= %d\n",k);fflush(NULL);
		if(FRACTIONS[k]>0){
			for(i=0;i<sist->N_Chain_Frac[k];i++){ //Growth of the brush chains

				sist->Graft[chain_index].x=Param->Colloid_Radius*sist->Sphere_Anchor[chain_index].x+sist->Colloid.x;
				sist->Graft[chain_index].y=Param->Colloid_Radius*sist->Sphere_Anchor[chain_index].y+sist->Colloid.y;
				sist->Graft[chain_index].z=Param->Colloid_Radius*sist->Sphere_Anchor[chain_index].z+sist->Colloid.z;
				sist->Graft[chain_index].lammpstype=index;
				//fprintf(fp,"%d %d %d 0 %lf %lf %lf\n",index,3,N/NATOM+1,sist->Graft[chain_index].x,sist->Graft[chain_index].y,sist->Graft[chain_index].z);	
				//index++;
				sprintf(line,"%d%s",chain_index+2,type_chain);
				if(Param->N_Polymer>1){				
					for(j=1;j<Param->N_Polymer;j++){
						translatex=sist->Graft[chain_index].x-sist->Colloid.x;
						translatey=sist->Graft[chain_index].y-sist->Colloid.y;
						translatez=sist->Graft[chain_index].z-sist->Colloid.z;

						module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

						translatex=(translatex*Param->Polymer_Size*j)/module+sist->Graft[chain_index].x;
						translatey=(translatey*Param->Polymer_Size*j)/module+sist->Graft[chain_index].y;
						translatez=(translatez*Param->Polymer_Size*j)/module+sist->Graft[chain_index].z;
						fprintf(fp,"%d %d %d 0 %lf %lf %lf\n",index,3+k,N/NATOM+2,translatex,translatey,translatez);


						if(j==1) 	fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,"OH",index_martini-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);
						if(j>1) 	fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,"C2O",index_martini-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);

						index++;
						index_martini++;
					}
				}
				translatex=sist->Graft[chain_index].x-sist->Colloid.x;
				translatey=sist->Graft[chain_index].y-sist->Colloid.y;
				translatez=sist->Graft[chain_index].z-sist->Colloid.z;

				module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

				translatex=(translatex*Param->Polymer_Size*(Param->N_Polymer))/module+sist->Graft[chain_index].x;
				translatey=(translatey*Param->Polymer_Size*(Param->N_Polymer))/module+sist->Graft[chain_index].y;
				translatez=(translatez*Param->Polymer_Size*(Param->N_Polymer))/module+sist->Graft[chain_index].z;
				fprintf(fp,"%d %d %d 0 %lf %lf %lf\n",index,3+k,N/NATOM+3+k,translatex,translatey,translatez); // End functionalized monomer
				fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,"P2",index_martini-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);
				index++;
				index_martini++;
				if(k>3){ //Charged group PEA and PEB 

					translatex=sist->Graft[chain_index].x-sist->Colloid.x;
					translatey=sist->Graft[chain_index].y-sist->Colloid.y;
					translatez=sist->Graft[chain_index].z-sist->Colloid.z;

					module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

					translatex=(translatex*Param->Polymer_Size*(Param->N_Polymer+1))/module+sist->Graft[chain_index].x;
					translatey=(translatey*Param->Polymer_Size*(Param->N_Polymer+1))/module+sist->Graft[chain_index].y;
					translatez=(translatez*Param->Polymer_Size*(Param->N_Polymer+1))/module+sist->Graft[chain_index].z;

					fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,type_end1,index_martini-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);
					index_martini++;

					translatex=sist->Graft[chain_index].x-sist->Colloid.x;
					translatey=sist->Graft[chain_index].y-sist->Colloid.y;
					translatez=sist->Graft[chain_index].z-sist->Colloid.z;

					module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

					translatex=(translatex*Param->Polymer_Size*(Param->N_Polymer+2))/module+sist->Graft[chain_index].x;
					translatey=(translatey*Param->Polymer_Size*(Param->N_Polymer+2))/module+sist->Graft[chain_index].y;
					translatez=(translatez*Param->Polymer_Size*(Param->N_Polymer+2))/module+sist->Graft[chain_index].z;

					fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,type_end2,index_martini-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);
					index_martini++;
				 }
				if(k>6){ //Acid group PEA ends with an  hydrogen
					translatex=sist->Graft[chain_index].x-sist->Colloid.x;
					translatey=sist->Graft[chain_index].y-sist->Colloid.y;
					translatez=sist->Graft[chain_index].z-sist->Colloid.z;

					module=sqrt(translatex*translatex+translatey*translatey+translatez*translatez);

					translatex=(translatex*Param->Polymer_Size*(Param->N_Polymer+3))/module+sist->Graft[chain_index].x;
					translatey=(translatey*Param->Polymer_Size*(Param->N_Polymer+3))/module+sist->Graft[chain_index].y;
					translatez=(translatez*Param->Polymer_Size*(Param->N_Polymer+3))/module+sist->Graft[chain_index].z;
					sprintf(line,"%d%s",Param->N_Chains+index2+1,"H  ");
					fprintf(fp_martini2,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,"POS",tot_chain_gromacs+index2+Param->Taken_Fillers,translatex/10,translatey/10,translatez/10);
					index2++;
				}
				sist->Graft[chain_index].end=index;	


				sist->Graft[chain_index].prot=k;
				chain_index++;
			}
		}
	}
	fprintf(fp,"\n\n");
	fprintf(fp_martini2," %lf %lf %lf\n",sist->box_x/10,sist->box_y/10,sist->box_z/10);
	fclose(fp_martini);
	fclose(fp_martini2);
	
	sprintf(line,"%s.gro2",Param->prefix);
	fp_martini=fopen(line,"w");
	for(i=0;i<Param->N_Fillers;i++){
		if(sist->Colloid_Fillers[i].flag==1){
			translatex=Param->Colloid_Radius*sist->Colloid_Fillers[i].x+sist->Colloid.x;
			translatey=Param->Colloid_Radius*sist->Colloid_Fillers[i].y+sist->Colloid.y;
			translatez=Param->Colloid_Radius*sist->Colloid_Fillers[i].z+sist->Colloid.z;
			
			sprintf(line,"%d%s",chain_index+2,"SUR");
			fprintf(fp_martini,"%8s%7s%5d%8.3lf%8.3lf%8.3lf\n",line,"OH",index-sist->Colloid.lammpstype+1,translatex/10,translatey/10,translatez/10);
			index++;
			chain_index++;
		}
	}
	
	
	fclose(fp_martini);
	
	return;
}	 */
void write_protein_bonds (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers){
	
	int i,k;
	int Protein_id;
	int init_bondtype=sist->bondtype;
	for(k=0;k<Param->N_Proteins;k++){
		sist->bondtype=init_bondtype;
		Protein_pointers[k]->bondtype_head=sist->bondtype;
		for(i=ATOM_CA;i<Protein_pointers[k]->N-NATOM;i+=NATOM){
			Protein_id=Protein_pointers[k]->Atom[i].lammps_idx;
			if(Protein_pointers[k]->Ignore[i]==0){ //Check if the protein bonds should be ignored or not
				fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->bondtype,Protein_id,Protein_id+1);
				sist->N_Bonds++;
				sist->bondtype++;
			}
		}	
		
	}
	return;
}


void write_linker_bonds (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers){

	int i,k,m;
	int idx=0;
	for(m=0;m<Param->N_Proteins;m++){
		for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
			if(Protein_pointers[m]->Mask[i]==1) {	
				fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,Protein_pointers[m]->Linker_Anchors[i].lammps_idx,Protein_pointers[m]->Linker[i].lammps_idx);
				idx=Protein_pointers[m]->Linker[i].lammps_idx;
				sist->N_Bonds++;
				for(k=0;k<Param->Linker_Length-1;k++){
					fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,idx+k,idx+k+1);
					sist->N_Bonds++;
				}
			}	
		}
	}
	return;
}


	
void write_brush_bonds (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Brushes *Brush,int type){
	int k,i,l;
	int idx;
	int length_corr=1; 
	int last=0;
	
if(type==ALL) length_corr=0; //Includes Anchors
	
	switch(Param->Geometry){
	case SPHERE:
		if(type==ALL){
		last=Param->Polymer_Length;
	}else{
		last=Param->Polymer_Length-1;
	}
	
	for(l=0;l<Param->N_Colloids;l++){
		for(k=0;k<Param->N_Polymer;k++){
			if(Brush->Anchors[k].flag>-1){	
				idx=Brush->Anchors[k].lammps_idx+l*(Param->N_Polymer*last+1);
				
				for(i=0;i<Param->Polymer_Length-1-length_corr;i++){
					fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,idx+i,idx+i+1);
					sist->N_Bonds++;
				}
			}
		}
	}
	break;
	case CYLINDER:
 
   	
  
     for(k=0;k<Param->N_Polymer;k++){
       if(Brush->Anchors[k].flag>-1){  
         idx=Brush->Anchors[k].lammps_idx;
         
         for(i=0;i<Param->Polymer_Length-1-length_corr;i++){
           fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,idx+i,idx+i+1);
           sist->N_Bonds++;
         }
       }
     }
	 break;
	 case SLIT:
     for(k=0;k<Param->N_Polymer;k++){
       if(Brush->Anchors[k].flag>-1){  
         idx=Brush->Anchors[k].lammps_idx;
         
         for(i=0;i<Param->Polymer_Length-1-length_corr;i++){
           fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,idx+i,idx+i+1);
           sist->N_Bonds++;
         }
       }
     }
   
	break;
	case CHANNEL:
     for(k=0;k<Param->N_Polymer;k++){
       if(Brush->Anchors[k].flag>-1){  
         idx=Brush->Anchors[k].lammps_idx;
         
         for(i=0;i<Param->Polymer_Length-1-length_corr;i++){
           fprintf(f_atoms,"%d %d %d %d\n",sist->N_Bonds,sist->polymers_bonds_type,idx+i,idx+i+1);
           sist->N_Bonds++;
         }
       }
     }
   
	break;
}
	
	
 return;
}

void write_protein_angles (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers){
	int i,k;
	int Protein_id;
	int N_Angles=1;
	int Angles_types=1;
	for(k=0;k<Param->N_Proteins;k++){
		Angles_types=1;
		for(i=ATOM_CA;i<Protein_pointers[k]->N-2*NATOM;i+=NATOM){
				
			if((Protein_pointers[k]->Ignore[i]==0)&&(Protein_pointers[k]->Ignore[i+NATOM]==0)&&(Protein_pointers[k]->Ignore[i+2*NATOM]==0)){	
				Protein_id=Protein_pointers[k]->Atom[i].lammps_idx;
				fprintf(f_atoms,"%d %d %d %d %d\n",N_Angles,Angles_types,Protein_id,Protein_id+1,Protein_id+2);			
				N_Angles++;
				Angles_types++;
			}
		}		
	}
	return;
}
	
void write_protein_dihedrals (struct GlobalSistem *sist, FILE *f_atoms, struct Interaction_Param *Param, struct Proteins **Protein_pointers){	
	int i,ii,k;
	int N_Dihedrals=1;
	int Protein_id;
	int Angles_types=1;
	ii=1;
	for(k=0;k<Param->N_Proteins;k++){
		Angles_types=1;
		for(i=ATOM_CA;i<Protein_pointers[k]->N;i+=NATOM){
			if(i+3*NATOM<Protein_pointers[k]->N){			
				if((Protein_pointers[k]->Ignore[i]==0)&&(Protein_pointers[k]->Ignore[i+NATOM]==0)&&(Protein_pointers[k]->Ignore[i+2*NATOM]==0)&&(Protein_pointers[k]->Ignore[i+3*NATOM]==0)){		
					Protein_id=Protein_pointers[k]->Atom[i].lammps_idx;
					fprintf(f_atoms,"%d %d %d %d %d %d\n",N_Dihedrals,Angles_types,Protein_id,Protein_id+1,Protein_id+2,Protein_id+3);
					N_Dihedrals++;
					Angles_types++;
				}
			}
			ii++;
		}
	}
	return;
}


void allocatePDB (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	char line[1024],lett[10],old_lett_c,ammino[10],variation[10];
	FILE *fp;
	int i,p,j;
	double x,y,z;
	Protein->Nres=0;
	
	fp=fopen(Param->pdb_filename,"r");
	if ( fp == NULL) {
		printf ("File %s not found!!!!\n",Param->pdb_filename);fflush(NULL);
		Param->protein_flag=NO;
		Protein->Nres=1;
	}else{
		Param->protein_flag=YES;
		//N=-1;
		while(fgets(line,sizeof(line),fp) !=NULL){
			old_lett_c=variation[0];
			for(i=0;i<10;i++){
				variation[i]=ammino[i]=lett[i]=0;
			}
			
			p=sscanf (line,"%s %s %lf %lf %lf %s\n",lett,ammino,&x,&y,&z,variation);
			
			if((variation[0]=='B')&&(old_lett_c==0)) variation[0]=0;
			if(p>=5){
				if((lett[0]=='C')&&(lett[1]=='A')&&((variation[0]=='A')||(variation[0]==0))) Protein->Nres++;
			}
		}
		fclose(fp);
	}
	printf("PDB Read\n");fflush(NULL);
	
	Protein->N=Protein->Nres*NATOM;
	Protein->Atom=(struct part *)calloc(Protein->N,sizeof(struct part));
	Protein->Atom_read=(struct part *)calloc(Protein->N,sizeof(struct part));
	Protein->gmass=(struct part*)calloc(1,sizeof(struct part));
	Protein->group_id=(char *)calloc(Protein->N*10,sizeof(char));
	
	
	strcpy(Protein->Atoms_name[ATOM_N],"N  ");
	strcpy(Protein->Atoms_name[ATOM_CA],"CA ");
	strcpy(Protein->Atoms_name[ATOM_C],"C  ");
	strcpy(Protein->Atoms_name[ATOM_O],"O  ");
	strcpy(Protein->Atoms_name[ATOM_H],"H  ");
	
	
	strcpy(Protein->Amminoacids[1],"ALA");
	strcpy(Protein->Amminoacids[2],"CYS");
	strcpy(Protein->Amminoacids[3],"ASP");
	strcpy(Protein->Amminoacids[4],"GLU");
	strcpy(Protein->Amminoacids[5],"PHE");
	strcpy(Protein->Amminoacids[6],"GLY");
	strcpy(Protein->Amminoacids[7],"HIS");
	strcpy(Protein->Amminoacids[8],"ILE");
	strcpy(Protein->Amminoacids[9],"LYS");
	strcpy(Protein->Amminoacids[10],"LEU");
	strcpy(Protein->Amminoacids[11],"MET");
	strcpy(Protein->Amminoacids[12],"ASN");
	strcpy(Protein->Amminoacids[13],"PRO");
	strcpy(Protein->Amminoacids[14],"GLN");
	strcpy(Protein->Amminoacids[15],"ARG");
	strcpy(Protein->Amminoacids[16],"SER");
	strcpy(Protein->Amminoacids[17],"THR");
	strcpy(Protein->Amminoacids[18],"VAL");
	strcpy(Protein->Amminoacids[19],"TRP");
	strcpy(Protein->Amminoacids[20],"TYR");
	
	printf("Atoms allocated\n");fflush(NULL);
	
	
	for(i=0;i<Protein->N;i+=NATOM){
		for(j=0;j<NATOM;j++){
			Protein->Atom[i+j].x=0;
			Protein->Atom[i+j].y=0;
			Protein->Atom[i+j].z=0;
			Protein->Atom[i+j].id=j;
		}
	}
	printf("Atoms Init\n");fflush(NULL);
	
	
	
	
	return;
}

void readMasks(struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein) {
    char line[1024];
    int ID_Masked = 0;
    int NInt = 0;
    int i, p;
    FILE *fp = NULL;

    Protein->Mask = (int *)calloc(Protein->N, sizeof(int));
    printf("Protein->Mask Read -%s-\n", Param->mask);
    fflush(NULL);
    fp = fopen(Param->mask, "r");

    if ((fp == NULL) || (Param->Simul_type == TRANS)) {
        if (fp == NULL) printf("Protein->Mask File %s not found. No cross linked residues to the wall\n", Param->mask);
        if (Param->Simul_type == TRANS) printf("WARNING! Ignoring Protein->Mask File %s in TRANS simulations\n", Param->mask);
        Param->mask_flag = NO;
    } else {
        Param->mask_flag = YES;
        if (fgets(line, sizeof(line), fp) != NULL) {
            p = sscanf(line, "%d", &Protein->NMasked);
            printf("Protein->Mask Size %d\n", Protein->NMasked);
            fflush(NULL);
            for (i = 0; i < Protein->NMasked; i++) {
                if (fgets(line, sizeof(line), fp) != NULL) {
                    p = sscanf(line, "%d", &ID_Masked);
                    printf("Protein->Mask idx %d %d\n", i, (ID_Masked - 1) * NATOM + ATOM_CA);
                    fflush(NULL);
                    Protein->Mask[(ID_Masked - 1) * NATOM + ATOM_CA] = 1;
                }
            }
        }
        fclose(fp);
    }

    printf("Mask_wall_file Read -%s-\n\n", Param->mask_wall_file);
    fflush(NULL);

    Protein->Mask_wall = (int *)calloc(Protein->N, sizeof(int));
    fp = fopen(Param->mask_wall_file, "r");
    if (fp == NULL) {
        printf("Protein->Mask File %s not found. No specific interactions with the wall will be added (e.g. charges)\n", Param->mask_wall_file);
        Param->mask_wall_flag = NO;
    } else {
        Param->mask_wall_flag = YES;
        if (fgets(line, sizeof(line), fp) != NULL) {
            p = sscanf(line, "%d", &NInt);
            for (i = 0; i < NInt; i++) {
                if (fgets(line, sizeof(line), fp) != NULL) {
                    p = sscanf(line, "%d", &ID_Masked);
                    Protein->Mask_wall[(ID_Masked - 1) * NATOM + ATOM_CA] = 1;
                }
            }
        }
        fclose(fp);
    }
    printf("Mask Wall Done\n\n");
    fflush(NULL);

    Protein->Linker_Anchors = (struct part*)calloc(Protein->N, sizeof(struct part));
    Protein->Linker = (struct part*)calloc(Protein->N, sizeof(struct part));

    fp = fopen(Param->mask_prot, "r");
    if (fp == NULL) {
        printf("WARNING! Mask Protonation File %s not found. Ignore protonations\n", Param->mask);
        fflush(NULL);
        for (i = ATOM_CA; i < Protein->N; i += NATOM) {
            Protein->Atom[i].prot = END_NEUTRAL;
        }
    } else {
        for (i = ATOM_CA; i < Protein->N; i += NATOM) {
            if (fgets(line, sizeof(line), fp) != NULL) {
                p = sscanf(line, "%d", &Protein->Atom[i].prot);
            }
        }
        fclose(fp);
    }

    printf("Protein->Mask Protonation Done\n");
    fflush(NULL);

    return;
}


void readIgnoreBonds(struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein) {
    int i, p;
    int ID_Skip = 0;
    char line[1024];
    FILE *fp = NULL;

    Protein->Ignore = (int *)calloc(Protein->N, sizeof(int));
    fp = fopen(Param->ignore, "r");
    if (fp == NULL) {
        printf("Ignore File %s not found. No bonds will be ignored\n", Param->ignore);
        for (i = ATOM_CA; i < Protein->N - NATOM; i += NATOM) {
            if (Protein->Ignore[i] != 0) {
                printf("There are ignored bonds and that is not possible\n");
                exit(0);
            }
        }
    } else {
        if (fgets(line, sizeof(line), fp) != NULL) {
            p = sscanf(line, "%d", &Protein->NSkip);
            for (i = 0; i < Protein->NSkip; i++) {
                if (fgets(line, sizeof(line), fp) != NULL) {
                    p = sscanf(line, "%d", &ID_Skip);
                    Protein->Ignore[(ID_Skip - 1) * NATOM + ATOM_CA] = 1;
                }
            }
        }
        fclose(fp);
    }
    printf("Ignore Done\n");
    fflush(NULL);
}







void readProtCoordinates (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	char line[1024],lett[10],old_lett_c,ammino[10],variation[10];
	int *Index_Seq=NULL;
	char **Seq=NULL;
	int i,p,j;
	int ii;
	int NN;
	double x,y,z;
	FILE *fp=NULL;
	
	Index_Seq=(int *)calloc(Protein->N,sizeof(int));
	Seq=(char **)calloc(Protein->N,sizeof(char*));
	for(i=0;i<Protein->N;i++){
		Seq[i]=(char *)calloc(4,sizeof(char));
	}
	
	fp=fopen(Param->pdb_filename,"r");
	if ( fp == NULL) {
		printf ("File %s not found!!!!\n",Param->pdb_filename);
		fflush(NULL);
		exit(EXIT_FAILURE);
	}else{
		NN=0;
		i=0;j=0;
		
		while(fgets(line,sizeof(line),fp) !=NULL){
			old_lett_c=variation[0];
			for(i=0;i<10;i++){
				variation[i]=ammino[i]=lett[i]=0;
			}
			p=sscanf (line,"%s %s %lf %lf %lf %s\n",lett,ammino,&x,&y,&z,variation);
			if((variation[0]=='B')&&(old_lett_c==0)) variation[0]=0;
			#ifdef PROGRESS
			printf ("-> %d %c%c %c %c %s\n",NN,lett[0],lett[1],variation[0],old_lett_c,ammino);fflush(NULL);
			#endif
			if((p>=5)&&((variation[0]=='A')||(variation[0]==0))){
				switch (lett[0]){
					case 'N':
					if(lett[1]==0){						
						Protein->Atom[NN].x=x;
						Protein->Atom[NN].y=y;
						Protein->Atom[NN].z=z;
						Protein->Atom[NN].id=ATOM_N;
						//printf ("%d/%d %d %s\n",NN,N,NN*5/4,Atoms_name[Protein->Atom[NN].id]);
						for(j=0;j<NATOM;j++){
							Seq[NN*NATOM/4+j][0]=ammino[0];
							Seq[NN*NATOM/4+j][1]=ammino[1];
							Seq[NN*NATOM/4+j][2]=ammino[2];
							for(ii=1;ii<21;ii++){
								
								if(strncmp(Protein->Amminoacids[ii],Seq[NN*NATOM/4+j],3)==0){
									//printf("NN=%d ii=%d %s  %s\n",NN,ii,Amminoacids[ii],Seq[NN*NATOM/4+j]);
									Index_Seq[NN*NATOM/4+j]=ii;
									ii=21;
								}
							}
						}
						
						
						NN++;
						
						
					}
					break;
					case 'C':
					if(lett[1]=='A'){
						
						
						Protein->Atom[NN].x=x;
						Protein->Atom[NN].y=y;
						Protein->Atom[NN].z=z;
						Protein->Atom[NN].id=ATOM_CA;
						//printf ("%d %s\n",NN,Atoms_name[Protein->Atom[NN].id]);
						
						NN++;
						
					}else{
						if(lett[1]==0){
							
							
							
							
							Protein->Atom[NN].x=x;
							Protein->Atom[NN].y=y;
							Protein->Atom[NN].z=z;
							Protein->Atom[NN].id=ATOM_C;
							//printf ("%d %s\n",NN,Atoms_name[Protein->Atom[NN].id]);
							
							NN++;
							
						}
					}
					break;
					case 'O':
					if(lett[1]==0){
						
						
						
						Protein->Atom[NN].x=x;
						Protein->Atom[NN].y=y;
						Protein->Atom[NN].z=z;
						Protein->Atom[NN].id=ATOM_O;
						//printf ("%d %s\n",NN,Atoms_name[Protein->Atom[NN].id]);
						
						NN++;
						
						
					}
					break;
					/*case 'H':
					if(lett[1]==' '){
						Protein->Atom[NN].x=x;
						Protein->Atom[NN].y=y;
						Protein->Atom[NN].z=z;
						Protein->Atom[NN].id=ATOM_H;
						printf ("%d %s\n",NN,Atoms_name[Protein->Atom[NN].id]);
						NN++;
						
						
						
					}
					break;*/
				}
			}	
			
			
		}
		 fclose(fp);
	}
	Protein->Nread=NN;
	
	for(i=0;i<Protein->N;i++){
		
		Protein->Atom_read[i].x=Protein->Atom[i].x;
		Protein->Atom_read[i].y=Protein->Atom[i].y;
		Protein->Atom_read[i].z=Protein->Atom[i].z;
		Protein->Atom_read[i].id=Protein->Atom[i].id;
	}
	

	for(i=0;i<Protein->N;i++){
		Protein->Atom[i].res=Index_Seq[i];
	}
	Angle_generator(sist,Protein);
	
	printf("Angles Generated\n");fflush(NULL);
	reorderProtCoordinates (sist,Param,Protein);
}


void reorderProtCoordinates (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	double Order_x[6],Order_y[6],Order_z[6],Order_id[6];
	struct part *Atom2=NULL;
	int index;
	int i,j,k;
	
	Atom2=(struct part *)calloc(Protein->N,sizeof(struct part));
	
	for(i=0;i<partint((double)(Protein->Nread)/4.0);i++){
		for(j=1;j<NATOM;j++){
			Order_id[j]=Protein->Atom[i*4+j-1].id;
			Order_x[j]=Protein->Atom[i*4+j-1].x;
			Order_y[j]=Protein->Atom[i*4+j-1].y;
			Order_z[j]=Protein->Atom[i*4+j-1].z;
			
		}
		sort2(4,Order_id,Order_x);
		for(j=1;j<NATOM;j++){
			Order_id[j]=Protein->Atom[i*4+j-1].id;
			
		}
		sort2(4,Order_id,Order_y);
		for(j=1;j<NATOM;j++){
			Order_id[j]=Protein->Atom[i*4+j-1].id;
			
		}
		sort2(4,Order_id,Order_z);
		for(j=1;j<NATOM;j++){
			Order_id[j]=Protein->Atom[i*4+j-1].id;
			
		}
		sort2(4,Order_id,Order_id);
		
		
		for(j=1;j<NATOM;j++){
			Protein->Atom[i*4+j-1].id=(int)(Order_id[j]);
			Protein->Atom[i*4+j-1].x=Order_x[j];
			Protein->Atom[i*4+j-1].y=Order_y[j];
			Protein->Atom[i*4+j-1].z=Order_z[j];
			
		}
		
	}
	k=0;
	for(i=0;i<partint((double)(Protein->Nread)/4.0);i++){
		for(j=0;j<4;j++){
			index=i*4+j;
			
			Atom2[k].id=Protein->Atom[index].id;
			Atom2[k].x=Protein->Atom[index].x;
			Atom2[k].y=Protein->Atom[index].y;
			Atom2[k].z=Protein->Atom[index].z;
			k++;
		}
		Atom2[k].x=0;
		Atom2[k].y=0;
		Atom2[k].z=0;
		Atom2[k].id=ATOM_H;
		k++;
	}
	for(i=0;i<Protein->N;i++){
		Protein->Atom[i].id=Atom2[i].id;
		Protein->Atom[i].x=Atom2[i].x;
		Protein->Atom[i].y=Atom2[i].y;
		Protein->Atom[i].z=Atom2[i].z;
	}
	
	
	for(i=0;i<partint((double)(Protein->N)/NATOM);i++){
		for(j=1;j<6;j++){
			Order_id[j]=Protein->Atom[i*NATOM+j-1].id;
			Order_x[j]=Protein->Atom[i*NATOM+j-1].x;
			Order_y[j]=Protein->Atom[i*NATOM+j-1].y;
			Order_z[j]=Protein->Atom[i*NATOM+j-1].z;
			
			
		}
		sort2(NATOM,Order_id,Order_x);
		for(j=1;j<6;j++){
			Order_id[j]=Protein->Atom[i*NATOM+j-1].id;
			
		}
		sort2(NATOM,Order_id,Order_y);
		for(j=1;j<6;j++){
			Order_id[j]=Protein->Atom[i*NATOM+j-1].id;
			
		}
		sort2(NATOM,Order_id,Order_z);
		for(j=1;j<6;j++){
			Order_id[j]=Protein->Atom[i*NATOM+j-1].id;
			
		}
		sort2(NATOM,Order_id,Order_id);
		
		
		for(j=1;j<6;j++){
			Protein->Atom[i*NATOM+j-1].id=(int)(Order_id[j]);
			Protein->Atom[i*NATOM+j-1].x=Order_x[j];
			Protein->Atom[i*NATOM+j-1].y=Order_y[j];
			Protein->Atom[i*NATOM+j-1].z=Order_z[j];
			
		}

	}
	
	return;
}





void test_protein (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	int Error=0;
	int i,j;
	
	for(i=0;i<partint((double)(Protein->N)/NATOM);i++){		
		j=0;
		if(Protein->Atom[i*NATOM+j].id!=ATOM_N){
			printf("Wrong N id %d %d for protein %s\n",j,i*NATOM+j,Param->pdb_filename);
			fflush(NULL);
			Error=1;
		}
		//printf("i= %d Seq= %s Index_Seq= %d\n",i*NATOM+j,Seq[i*NATOM+j],Index_Seq[i*NATOM+j]);
		j++;
		if(Protein->Atom[i*NATOM+j].id!=ATOM_CA){
			printf("Wrong CA id %d %d for protein %s\n",j,i*NATOM+j,Param->pdb_filename);
			fflush(NULL);
			Error=1;
		}
		//Index_Seq[i*NATOM+j]=Index_Seq[i*NATOM];
		//strcpy(Seq[i*NATOM+j],Seq[i*NATOM]);
		//printf("i= %d Seq= %s Index_Seq= %d\n",i*NATOM+j,Seq[i*NATOM+j],Index_Seq[i*NATOM+j]);
		j++;
		if(Protein->Atom[i*NATOM+j].id!=ATOM_C){
			printf("Wrong C id %d %d for protein %s\n",j,i*NATOM+j,Param->pdb_filename);
			fflush(NULL);
			Error=1;
		}
		//Index_Seq[i*NATOM+j]=Index_Seq[i*NATOM];
		//strcpy(Seq[i*NATOM+j],Seq[i*NATOM]);
		//printf("i= %d Seq= %s Index_Seq= %d\n",i*NATOM+j,Seq[i*NATOM+j],Index_Seq[i*NATOM+j]);
		j++;
		if(Protein->Atom[i*NATOM+j].id!=ATOM_O){
			printf("Wrong O id %d %d for protein %s\n",j,i*NATOM+j,Param->pdb_filename);
			fflush(NULL);
			Error=1;
		}
		//Index_Seq[i*NATOM+j]=Index_Seq[i*NATOM];
		//strcpy(Seq[i*NATOM+j],Seq[i*NATOM]);
		//printf("i= %d Seq= %s Index_Seq= %d\n",i*NATOM+j,Seq[i*NATOM+j],Index_Seq[i*NATOM+j]);
		j++;
		if(Protein->Atom[i*NATOM+j].id!=ATOM_H){
			printf("Wrong H id %d %d for protein %s\n",j,i*NATOM+j,Param->pdb_filename);
			fflush(NULL);
			Error=1;
		}
		//Index_Seq[i*NATOM+j]=Index_Seq[i*NATOM];
		//strcpy(Seq[i*NATOM+j],Seq[i*NATOM]);
		//printf("i= %d Seq= %s Index_Seq= %d\n",i*NATOM+j,Seq[i*NATOM+j],Index_Seq[i*NATOM+j]);	
	}
	if(Error==1){
		for(i=0;i<Protein->N;i++){	
			printf ("%d %s %lf %lf %lf\n",i,Protein->Atoms_name[Protein->Atom[i].id],Protein->Atom[i].x,Protein->Atom[i].y,Protein->Atom[i].z);
		}
		for(i=0;i<Protein->N;i++){	
			printf ("---READ---- %d %s %lf %lf %lf\n",i,Protein->Atoms_name[Protein->Atom_read[i].id],Protein->Atom_read[i].x,Protein->Atom_read[i].y,Protein->Atom_read[i].z);
		}
		fflush(NULL);
	}
	
	
	return;
}

void write_sequence (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp){
	
	int i;
	char sequenza[]="OACDEFGHIKLMNPQRSTVWY";
	
	if (fp == NULL) {
		printf ("Could not create %s-out.seq\n",Param->prefix);
		fflush(NULL);
		exit(EXIT_FAILURE);
	}else{
		/*for(i=0;i<Protein->N;i++){
			printf("i=%d %s\n",i,Seq[i]);
		}*/
		for(i=ATOM_CA;i<Protein->N;i+=NATOM){
			fprintf(fp,"%c",sequenza[Protein->Atom[i].res]);
		}
		fprintf(fp,"\n");
	}
	
	
	
	return;
	
}


void write_pdb (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein,FILE *fp){
	
	int i,kk;
	
	if ( fp == NULL) {
		printf ("Could not create %s-out.pdb file\n",Param->prefix);
		fflush(NULL);
		exit(EXIT_FAILURE);
	}else{
		
		
		kk=1;
		for(i=0;i<Protein->N;i++){
			fprintf(fp,"ATOM  %5d  %s %s  %4d    %8.3lf%8.3lf%8.3lf\n",i+1,Protein->Atoms_name[Protein->Atom[i].id],Protein->Amminoacids[Protein->Atom[i].res],kk,
			Protein->Atom[i].x,
			Protein->Atom[i].y,
			Protein->Atom[i].z);
			if(Protein->Atom[i].id==ATOM_H) kk++;
		}
	
	}
	
	
	
	return;
	
}
void write_binary (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein,FILE *fp){
	
	int kk,i,id;
	double x,y,z;
	
	if ( fp == NULL) {
		printf ("Could not create %s-out.bin file\n",Param->prefix);
		fflush(NULL);
		exit(EXIT_FAILURE);
	}else{
	
		
		kk=1;
		for(i=0;i<Protein->N;i++){
			
			x=Protein->Atom[i].x;
			y=Protein->Atom[i].y;
			z=Protein->Atom[i].z;
			id=Protein->Atom[i].id;
			fwrite(&x,sizeof(double),1,fp);
			fwrite(&y,sizeof(double),1,fp);
			fwrite(&z,sizeof(double),1,fp);
			fwrite(&id,sizeof(int),1,fp);
			
			if(Protein->Atom[i].id==ATOM_H) kk++;
		}
		
	}
	
	
	
	return;
	
}

void protein_gmass (struct GlobalSistem *sist, struct Proteins *Protein){
	int i;
	
	Protein->gmass->x=0;
	Protein->gmass->y=0;
	Protein->gmass->z=0;
	
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){	
	
		Protein->gmass->x+=Protein->Atom[i].x;
		Protein->gmass->y+=Protein->Atom[i].y;
		Protein->gmass->z+=Protein->Atom[i].z;
		
	}	
	Protein->gmass->x/=Protein->Nres;
	Protein->gmass->y/=Protein->Nres;
	Protein->gmass->z/=Protein->Nres;
	
	return;
	
}



void protein_rotate_shift (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	
		
	int rot_flag=1;
	double Rot_axis_X;
	double Rot_axis_Y;
	double Rot_axis_Z;
	double Rot_angle;
	double s=0;
	double V1=0;
	double V2=0;
	double min_z_prot=10000;
	double max_z_prot=-10000;
	double min_y_prot=10000;
	double max_y_prot=-10000;
	double min_x_prot=10000;
	double max_x_prot=-10000;
	double translatez;
	int i;
	
	protein_gmass(sist,Protein);
	
	/**************** RANDOM ORIENTATION OF THE PROTEIN *****************/
	do{
		V1=2*ran3(&Param->seed)-1;
		V2=2*ran3(&Param->seed)-1;
		s=V1*V1+V2*V2;
	}while(s>=1);
	
	Rot_axis_X=2*V1*sqrt(1-s);
	Rot_axis_Y=2*V2*sqrt(1-s);
	Rot_axis_Z= (1-2*s);
	Rot_angle=((2*PI*ran3(&Param->seed))-PI);
	
	//printf("%lf %lf %lf || %lf\n",Rot_axis_X,Rot_axis_Y,Rot_axis_Z,Rot_angle);fflush(NULL);
	for(i=0;i<Protein->N;i++){
		Protein->Atom[i].x-=Protein->gmass->x;
		Protein->Atom[i].y-=Protein->gmass->y;
		Protein->Atom[i].z-=Protein->gmass->z;
		
		//printf("BEFORE || %lf %lf %lf\n",Protein->Atom[i].x,Protein->Atom[i].y,Protein->Atom[i].z);fflush(NULL);
		if(Param->random_rotation_flag == YES) {
			Rotation(&Protein->Atom[i].x,&Protein->Atom[i].y,&Protein->Atom[i].z,Rot_angle,Rot_axis_X,Rot_axis_Y,Rot_axis_Z,rot_flag);
		}
		//printf("AFTER || %lf %lf %lf\n",Protein->Atom[i].x,Protein->Atom[i].y,Protein->Atom[i].z);fflush(NULL);
		//Rotation_Unwrap(sist,atom,Rot->alpha,Rot->RotX,Rot->RotY,Rot->RotZ,rot_flag);
		rot_flag=0;
		
		
		Protein->Atom[i].x+=Protein->gmass->x;
		Protein->Atom[i].y+=Protein->gmass->y;
		Protein->Atom[i].z+=Protein->gmass->z;
	}
	/******************* NEW CMASS *************/
	
	protein_gmass(sist,Protein);
	
	/**************** SHIFT PROTEIN TO BOX CENTER *****************/
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		if(min_z_prot>Protein->Atom[i].z) min_z_prot=Protein->Atom[i].z;
		if(max_z_prot<Protein->Atom[i].z) max_z_prot=Protein->Atom[i].z;
		
		if(min_y_prot>Protein->Atom[i].y) min_y_prot=Protein->Atom[i].y;
		if(max_y_prot<Protein->Atom[i].y) max_y_prot=Protein->Atom[i].y;
		
		if(min_x_prot>Protein->Atom[i].x) min_x_prot=Protein->Atom[i].x;
		if(max_x_prot<Protein->Atom[i].x) max_x_prot=Protein->Atom[i].x;
	}
	
	for(i=0;i<Protein->N;i++){
		
		switch(Param->Simul_type){
					case TRANS:
						//Protein->Atom[i].x-=min_x_prot-3; // Place the protein at the origin
						Protein->Atom[i].x-=Protein->gmass->x; // Place the protein at the origin
						Protein->Atom[i].y-=Protein->gmass->y;
						Protein->Atom[i].z-=Protein->gmass->z;
					break;
					case ANCHOR:
						if(Param->Geometry==CYLINDER){
							Protein->Atom[i].x-=Protein->gmass->x-sist->Half_Box_x;
							Protein->Atom[i].y-=Protein->gmass->y-sist->Half_Box_y;
							Protein->Atom[i].z-=Protein->gmass->z-sist->Half_Box_z;
						}else{
							Protein->Atom[i].x-=Protein->gmass->x; // Place the protein at the origin
							Protein->Atom[i].y-=Protein->gmass->y;
							Protein->Atom[i].z-=Protein->gmass->z;
						}
					break;				
			}
		
		
		
		
	}
	
	
	return;
}


void protein_properties (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein){
	
	int i,j;
	int ii,jj;
	double dx_CA_CA,dy_CA_CA,dz_CA_CA;
	
	double r2;
	
	/**************** TOTAL CONTACTS*****************/
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		
		for(j=ATOM_CA;j<Protein->N;j+=NATOM){
			if((j!=i)&&(j!=i+NATOM)&&(j!=i-NATOM)&&(j!=i-2*NATOM)&&(j!=i+2*NATOM)){
				dx_CA_CA=Protein->Atom[j].x-Protein->Atom[i].x;
				dy_CA_CA=Protein->Atom[j].y-Protein->Atom[i].y;
				dz_CA_CA=Protein->Atom[j].z-Protein->Atom[i].z;
				
				r2=dx_CA_CA*dx_CA_CA+dy_CA_CA*dy_CA_CA+dz_CA_CA*dz_CA_CA;
				if(r2>Protein->max_dist){
					Protein->max_dist=r2;
				}
				
			}
		}
		
	}	
	Protein->max_dist=sqrt(Protein->max_dist);
	
	Protein->Average_r=(double *)calloc(Protein->Nres,sizeof(double));
	Protein->Count_r2=(double *)calloc(Protein->Nres,sizeof(double));
	ii=1;
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		Protein->Average_r[ii-1]=0;
		Protein->Count_r2[ii-1]=DBL_EPSILON;
		jj=1;
		for(j=ATOM_CA;j<Protein->N;j+=NATOM){
			if(j>=i){
				dx_CA_CA=Protein->Atom[j].x-Protein->Atom[i].x;
				dy_CA_CA=Protein->Atom[j].y-Protein->Atom[i].y;
				dz_CA_CA=Protein->Atom[j].z-Protein->Atom[i].z;
				
				r2=dx_CA_CA*dx_CA_CA+dy_CA_CA*dy_CA_CA+dz_CA_CA*dz_CA_CA;
				if((r2<sist->ECA_Range2)&&(((ii>jj+ATOM_IGNORE)||(ii<jj-ATOM_IGNORE))&&((jj>ii+ATOM_IGNORE)||(jj<ii-ATOM_IGNORE)))&&(jj!=ii)){
					Protein->Average_r[ii-1]+=sqrt(r2);
					Protein->Count_r2[ii-1]++;
				}
			}
			
			jj++;
		}
		Protein->Average_r[ii-1]/=Protein->Count_r2[ii-1];
		//printf("Protein->Average_r[%d]= %e %e\n",ii-1,Protein->Average_r[ii-1],Protein->Count_r2[ii-1]);
		ii++;
	}		
	
	
	return;
	
}


void create_Protein (struct GlobalSistem *sist, struct Interaction_Param *Param,struct Proteins *Protein){
	
	FILE *fp=NULL;
	char line[1024];
	allocatePDB  (sist,Param,Protein);
	if(Param->protein_flag==YES){
		readProtCoordinates (sist,Param,Protein);  


		//Save Sequence
		sprintf(line,"%s-out.seq",Param->prefix);
		fp=fopen(line,"w");
		write_sequence(sist,Param,Protein,fp);
		fclose(fp);

		sprintf(line,"%s-out.bin",Param->prefix);
		fp=fopen(line,"w");
		write_binary(sist,Param,Protein,fp);
		fclose(fp);

		sprintf(line,"%s-out.pdb",Param->prefix);
		fp=fopen(line,"w");
		write_pdb(sist,Param,Protein,fp);
		fclose(fp);

		protein_properties(sist,Param,Protein);

		/* Test Protein*/


		test_protein(sist,Param,Protein);






		readIgnoreBonds (sist,Param,Protein);
		readMasks (sist,Param,Protein);
	}
	return;
}

void write_infile (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush, FILE *fp, int prev_step, int cur_step){
	
	int k;
	
	fprintf(fp, "units lj\n");
  	fprintf(fp, "atom_style full\n");
  	fprintf(fp, "atom_modify first forcegroup #needed for SRDmod\n");
  	fprintf(fp, "dimension 3\n");
  	fprintf(fp, "boundary p p p\n");
	
	if(cur_step==1)  fprintf(fp, "read_data %s-out.atomparam\n",Param->prefix);
	if(prev_step==1) fprintf(fp,"read_restart restart_1_minimize.*\n");
	if(prev_step==2){
		if(Param->N_Polymer>0){
			fprintf(fp,"read_restart restart_2_unfolding.*\n");
		}else{
			fprintf(fp, "read_data %s-out.atomparam\n",Param->prefix);
		}
	}
	if(prev_step==3) fprintf(fp,"read_restart restart_3_folding.*\n");
	if(prev_step==4) fprintf(fp,"read_restart restart_4a_docking.*\n");
	if(prev_step==5) fprintf(fp,"read_restart restart_5_anchored.*\n");
	if(prev_step==6) fprintf(fp,"read_restart restart_6_bound_lowT.*\n");
	if(prev_step==7) fprintf(fp,"read_restart restart_6_bound_highT.*\n");
	if(prev_step>=8) fprintf(fp,"read_restart restart_7_flow.*\n");
	
	write_groups(sist,Param,Protein_pointers,Brush,fp);
	write_variables(sist,Param,fp,Protein_pointers[0]);
	
	
	
	// harmonicbonds
	write_harmonicbonds(sist,Param,Protein_pointers[0],fp);
	
	
			
	if(Param->protein_flag==YES) write_angles(sist,Param,Protein_pointers[0],fp);
	
	fprintf(fp, "neighbor 0.3 bin\n");
  	fprintf(fp, "neigh_modify delay 0 every 10 check yes\n");
  	fprintf(fp, "neigh_modify include forcegroup\n");
  	fprintf(fp, "comm_modify mode single cutoff 80 group forcegroup vel yes\n\n");
	
	
	write_fix(Param,fp);
	if(cur_step==1){
		 write_1_minimize(fp,Param);
		 return;
	 }
	
	
	if(cur_step==2) {
		write_2_minimize(fp,Param);
		 return;
	 }
	

	if(Param->Simul_type==ANCHOR) write_fix_prot(Param,fp);
	write_harmonicbonds(sist,Param,Protein_pointers[0],fp);
	write_pair_coeff(sist,Param,Protein_pointers,Brush,fp);
	
	if(cur_step==3) {
		write_3_minimize(fp,Param);
		 return;
	 }
	
	
	if(cur_step==4) {
		for(k=0;k<Param->N_Proteins;k++){
		 write_4_docking(sist,fp,Param,Protein_pointers[k]);
		}	
		return;
	 }
	
	if(cur_step==5) {
		if(Param->protein_flag==YES) write_5_anchored(fp,Param);
		 return;
	 }
	
	if(Param->protein_flag==YES) write_fixbrooks(sist,Param,Protein_pointers,fp);
	
	if(cur_step==6){
		 if(Param->protein_flag==YES) write_6_bound_init_T(fp,Param);
		 return;
	 }
	
	if(cur_step==7){  
		if(Param->protein_flag==YES) write_6_bound_final_T(fp,Param);
		 return;
	 }
	
	
	fprintf(fp, "variable prot_k equal ${eps_h}*20\n");
	fprintf(fp, "variable linker_k equal ${eps_h}*20\n");
	write_harmonicbonds(sist,Param,Protein_pointers[0],fp);


	
	if(cur_step==8){
		 write_7_srd(fp,Param,NO,Param->init_flow_temp,cur_step-7);
		 return;
	 }
	if(cur_step>8){
		if(cur_step-7<=Param->steps_flow_temp+1){
		 	write_7_srd(fp,Param,YES,Param->init_flow_temp+Param->bin_flow_temp*(cur_step-8),cur_step-7); //linear change in the Temperature
		}else{
			write_7_srd(fp,Param,YES,1,cur_step-7); //linear change in the Temperature
		}
		 return;
	 }
	
	return;
	
}




void write_variables (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp, struct Proteins *Protein){
	
	int i;
	char sequenza[]="OACDEFGHIKLMNPQRSTVWY";
	
	fprintf(fp,"	# Hard Core potential variables\n");
	fprintf(fp,"variable seed 	equal %u\n",Param->useed);
	fprintf(fp,"variable pressure 	equal %lf\n",Param->pressure);
	fprintf(fp,"variable init_temp 	equal %lf\n",Param->init_temp);
	fprintf(fp,"variable final_temp 	equal %lf\n",Param->final_temp);
	fprintf(fp,"variable init_flow_temp 	equal %lf\n",Param->init_flow_temp);
	fprintf(fp,"variable final_flow_temp 	equal %lf\n",Param->final_flow_temp);
	
	fprintf(fp,"variable chi_s equal %lf #4.5 Hydrophobic ; 1.5 Moderate Phobic; -1 Philic\n\n",Param->chi_s);
			
	fprintf(fp,"variable sigma_C equal 3.000\n");
	fprintf(fp,"variable sigma_linker equal %lf\n",Param->Linker_Size);
	fprintf(fp,"variable sigma_anchor equal ${sigma_linker}\n");
	
	fprintf(fp,"variable sigma_C_anchor equal ${sigma_C}*0.5+${sigma_anchor}*0.5 # sigma_C+sigma_C_anchor/2\n");
	fprintf(fp,"variable sigma_C_linker equal ${sigma_C}*0.5+${sigma_linker}*0.5 # sigma_C+sigma_C_linker/2\n");
	fprintf(fp,"variable sigma_anchor_linker equal ${sigma_anchor}*0.5+${sigma_linker}*0.5 # sigma_anchor+sigma_anchor_linker/2\n\n");
	
	fprintf(fp,"variable cutoff_C equal ${sigma_C}*1.12246204830938\n");

	
	fprintf(fp,"variable range_mol equal ${sigma_C}*5\n");
	fprintf(fp,"# GO potential variables\n");
	fprintf(fp,"variable alpha equal 2.0*3.8\n");
	fprintf(fp,"variable eps_h equal %lf\n",Param->E_Scale);
	fprintf(fp,"variable GO equal -${eps_h}\n");
	fprintf(fp,"variable GO2 equal ${eps_h}\n");
	fprintf(fp,"variable prot_k equal ${eps_h}*5\n");
	fprintf(fp,"variable angle_k equal ${eps_h}*6\n");
	fprintf(fp,"variable angle_k2 equal ${eps_h}*6\n");
	fprintf(fp,"variable torsional_1 equal ${eps_h}*6\n");
	fprintf(fp,"variable torsional_2 equal ${eps_h}*5.4\n");
	
	fprintf(fp,"variable range_yukawa equal ${sigma_C}*10\n");
	fprintf(fp,"variable debye equal 1/7.0\n");
	fprintf(fp,"variable yukawa_strength_attr equal -0.025*${eps_h}\n"); // The original value was -0.1
	fprintf(fp,"variable yukawa_strength_repul equal 0\n");
	fprintf(fp,"variable hb_strength equal 25*${eps_h}\n"); // The original value was 100
	
	fprintf(fp,"variable cutoff_go equal ${sigma_C}*4.0\n");
	
	fprintf(fp,"# Linker Variables\n");
	fprintf(fp,"variable linker_k equal ${eps_h}\n"); // Stronger will intially detach the protein form the surface we make stronger later
	if(Param->protein_flag==YES){		
		fprintf(fp,"# SEQ = ");
		for(i=ATOM_CA;i<Protein->N;i+=NATOM){
			fprintf(fp,"%c",sequenza[Protein->Atom[i].res]);
		}
		fprintf(fp,"\n");
	}
	
	
	fprintf(fp,"	# SRD PORE\n");
	
	fprintf(fp,"variable        xlo equal 0.0\n");
	fprintf(fp,"variable        xhi equal %lf\n",sist->box_x);
	fprintf(fp,"variable        ylo equal 0.0\n");
	fprintf(fp,"variable        yhi equal %lf\n",sist->box_y);
	fprintf(fp,"variable        zlo equal 0.0\n");
	fprintf(fp,"variable        zhi equal %lf\n",sist->box_z);
	fprintf(fp,"variable        radius equal %lf\n",sist->radius);
	fprintf(fp,"variable        inradius equal  %lf\n",sist->inradius);
	fprintf(fp,"variable 	xlo_tube equal ${xlo}-400.0\n");
	fprintf(fp,"variable 	xhi_tube equal ${xhi}+400.0\n");
	fprintf(fp,"variable 	ylo_tube equal ${ylo}-400\n");
	fprintf(fp,"variable 	yhi_tube equal ${yhi}+400\n");
	fprintf(fp,"variable 	ylo_intube equal ${ylo}+%lf\n",sist->MOL_SHIFT);
	fprintf(fp,"variable 	yhi_intube equal ${yhi}-%lf\n",sist->MOL_SHIFT);
	fprintf(fp,"variable 	zlo_intube equal ${zlo}+%lf\n",sist->MOL_SHIFT);
	fprintf(fp,"variable 	zhi_intube equal ${zhi}-%lf\n",sist->MOL_SHIFT);
	
	fprintf(fp,"#Create regions of the cylinder\n");
	if(Param->Geometry==CYLINDER){ 
		fprintf(fp,"region poresrd cylinder x ${radius} ${radius} ${inradius} ${xlo} ${xhi}\n");
		fprintf(fp,"region poreinner cylinder x ${radius} ${radius} ${inradius} ${xlo_tube} ${xhi_tube}\n");
	 	fprintf(fp,"region poreouter cylinder x ${radius} ${radius} ${inradius} ${xlo_tube} ${xhi_tube} side out #This is the outsied of the inner pore and is used to define the intersect in which the anchor will be trapped; in order to make a tybe in LAMMPS the cylinder mus tbe longer than the box to erase with the PBC the base walls\n");
		fprintf(fp,"region porewall cylinder x	${radius} ${radius} ${radius} ${xlo_tube} ${xhi_tube} #The actual impentrable wall\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"region anchorwall intersect 2 porewall poreouter #This is th region where the anchor points once they enter are trapped\n");
			fprintf(fp,"group docked dynamic anchor region anchorwall #Dynamic group of the anchor inside the region anchorwall\n"); 
		}
		fprintf(fp,"group  free dynamic forcegroup region poreinner #Dynamic group of the force group in the simulation pore  region poreinner\n"); 
	}
	if(Param->Geometry==SPHERE){
		fprintf(fp,"region colloidinner sphere %lf %lf %lf %lf\n",sist->Half_Box_x,sist->Half_Box_y,sist->Half_Box_z,Param->Coll_R-3);
		fprintf(fp,"region colloidouter sphere %lf %lf %lf %lf side out\n",sist->Half_Box_x,sist->Half_Box_y,sist->Half_Box_z,Param->Coll_R-3);	
		//fprintf(fp,"region colloid sphere %lf %lf %lf %lf side out\n",sist->Half_Box_x,sist->Half_Box_y,sist->Half_Box_z,Param->Coll_R);		
		if(Param->mask_flag==YES) {			
			fprintf(fp,"group docked dynamic anchor region colloidinner #Dynamic group of the anchor inside the region anchorwall\n"); 
		}
		fprintf(fp,"group  free dynamic forcegroup region colloidouter #Dynamic group of the force group in the simulation pore  region poreinner\n"); 
	}
	
	if(Param->Geometry==SLIT){ 
		fprintf(fp,"region blocksrd block ${xlo} ${xhi} ${ylo} ${yhi} ${zlo_intube} ${zhi_intube}\n");
		fprintf(fp,"region blockinner block ${xlo_tube} ${xhi_tube} ${ylo_tube} ${yhi_tube} ${zlo_intube} ${zhi_intube}\n");
	 	fprintf(fp,"region blockouter block ${xlo_tube} ${xhi_tube} ${ylo_tube} ${yhi_tube} ${zlo_intube} ${zhi_intube} side out #This is the outsied of the inner pore and is used to define the intersect in which the anchor will be trapped; in order to make a tybe in LAMMPS the cylinder mus tbe longer than the box to erase with the PBC the base walls\n");
		fprintf(fp,"region blockwall block ${xlo_tube} ${xhi_tube} ${ylo_tube} ${yhi_tube} ${zlo} ${zhi} #The actual impentrable wall\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"region anchorwall intersect 2 blockwall blockouter #This is th region where the anchor points once they enter are trapped\n");
			fprintf(fp,"group docked dynamic anchor region anchorwall #Dynamic group of the anchor inside the region anchorwall\n"); 
		}
		fprintf(fp,"group  free dynamic forcegroup region blockinner #Dynamic group of the force group in the simulation block  region blockinner\n"); 
	}
	if(Param->Geometry==CHANNEL){ 
		fprintf(fp,"region blocksrd block ${xlo} ${xhi} ${ylo_intube} ${yhi_intube} ${zlo_intube} ${zhi_intube}\n");
		fprintf(fp,"region blockinner block ${xlo_tube} ${xhi_tube} ${ylo_intube} ${yhi_intube} ${zlo_intube} ${zhi_intube}\n");
	 	fprintf(fp,"region blockouter block ${xlo_tube} ${xhi_tube} ${ylo_intube} ${yhi_intube} ${zlo_intube} ${zhi_intube} side out #This is the outsied of the inner pore and is used to define the intersect in which the anchor will be trapped; in order to make a tybe in LAMMPS the cylinder mus tbe longer than the box to erase with the PBC the base walls\n");
		fprintf(fp,"region blockwall block ${xlo_tube} ${xhi_tube} ${ylo} ${yhi} ${zlo} ${zhi} #The actual impentrable wall\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"region anchorwall intersect 2 blockwall blockouter #This is th region where the anchor points once they enter are trapped\n");
			fprintf(fp,"group docked dynamic anchor region anchorwall #Dynamic group of the anchor inside the region anchorwall\n"); 
		}
		fprintf(fp,"group  free dynamic forcegroup region blockinner #Dynamic group of the force group in the simulation block  region blockinner\n"); 
	}
	

	
		

	
	fprintf(fp,"# Region grafting wall\n");
	
	
	return;
	
}
void write_groups (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers, struct Brushes *Brush,FILE *fp){
	
	int res_flag=0;	
	int i,j,ii,k,m;
	
	
	if(Param->protein_flag==YES){
	for(i=1;i<=21;i++){

		res_flag=0;
		for(m=0;m<Param->N_Proteins;m++){
			for(j=ATOM_CA;j<Protein_pointers[m]->N;j+=NATOM){
				if(Protein_pointers[m]->Atom[j].res==i) res_flag=1;		
			}
		}
		if(res_flag==1){
			fprintf(fp,"group %s id ",Protein_pointers[0]->Amminoacids[i]);
			for(m=0;m<Param->N_Proteins;m++){
				for(j=ATOM_CA;j<Protein_pointers[m]->N;j+=NATOM){
					if(Protein_pointers[m]->Atom[j].res==i) fprintf(fp,"%d ",Protein_pointers[m]->Atom[j].lammps_idx);

				}
			}
			fprintf(fp,"\n");
		}
	}

	fprintf(fp,"\n");
	
		if(Param->mask_wall_flag==YES){
			fprintf(fp,"group wall_int id ");
			for(m=0;m<Param->N_Proteins;m++){
				for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
					if(Protein_pointers[m]->Mask_wall[i]==YES) {
						fprintf(fp,"%d ",Protein_pointers[m]->Atom[i].lammps_idx);	

					}
				}
			}
			fprintf(fp,"\n");
		}
	fprintf(fp,"group temp id 1\n"); // This is a dummy group to be able to alway erase it and redefin it
	fprintf(fp,"group temp delete\n");
	fprintf(fp,"group temp id ");
	ii=1;
	for(m=0;m<Param->N_Proteins;m++){
		for(i=ATOM_CA;i<Protein_pointers[m]->N;i+=NATOM){
			fprintf(fp,"%d ",Protein_pointers[m]->Atom[i].lammps_idx);
		}
	}
	fprintf(fp,"\n");
	fprintf(fp,"group linker type %d\n",LINKER);
	fprintf(fp,"group protein union temp linker\n");
	fprintf(fp,"group linker delete\n");
	fprintf(fp,"group temp delete\n");
	
}else{
	fprintf(fp,"group protein type %d\n",LINKER);
}
	
	
	
	fprintf(fp,"group anchor type %d\n",ANCHOR_LINKER);
	fprintf(fp,"group polymer type %d ",POLYMER);
	for(k=END_POLYMER;k<END_POLYMER+POL_END_TYPES;k++){
		fprintf(fp,"%d ",k);
	}
	fprintf(fp,"\n");
	fprintf(fp,"group anchor_polymer type %d %d\n",ANCHOR_POLYMER,COLLOID);	
	fprintf(fp,"group forcegroup union protein polymer anchor\n");
	fprintf(fp,"group temp type %d\n",COLLOID);
	fprintf(fp,"group dock_group union protein anchor temp\n");	
	fprintf(fp,"group temp delete\n");
	

	return;
}

void write_pair_coeff (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers,struct Brushes *Brush, FILE *fp){
	
	
	int k,k1;
	

	if(Param->protein_flag==YES){
	//modified: fixed wrong sigma and shifting
		
			fprintf(fp,"pair_style \t	 hybrid/overlay  sigmoidal ${cutoff_go} soft%s ${sigma_C}  yukawa%s ${debye}  ${range_yukawa}\n\n",Param->GPUFF,Param->GPUFF);
		
		
	}else{
		
			fprintf(fp,"pair_style \t	 hybrid/overlay soft%s ${sigma_C}  yukawa%s ${debye}  ${range_yukawa}\n\n",Param->GPUFF,Param->GPUFF);
		
		
	}
	fprintf(fp,"pair_coeff * * soft%s 1000 \n",Param->GPUFF);
	
	
	if(Param->protein_flag==YES) write_protein_coeff(sist,Param,Protein_pointers[0],fp);
	
	
	// All vs 

	fprintf(fp,"pair_coeff * %d soft%s 1000 %lf\n",COLLOID,Param->GPUFF,Param->Coll_R+1.5);  // COLLOID LInker 
	fprintf(fp,"pair_coeff * %d soft%s 1000 ${sigma_C_linker}\n",LINKER,Param->GPUFF);  // LINKER LInker 
	fprintf(fp,"pair_coeff * %d soft%s 1000 ${sigma_C_linker}\n",ANCHOR_LINKER,Param->GPUFF);  // ANCHOR_LINKER Anchor 
	//fprintf(fp,"pair_coeff %d %d soft%s 0 ${sigma_C_linker}\n",COLLOID,ANCHOR_LINKER,Param->GPUFF);  // ANCHOR_LINKER Anchor 
	fprintf(fp,"pair_coeff * %d soft%s 100 ${sigma_C_linker}\n",POLYMER,Param->GPUFF);  // POLYMER Polymer
	fprintf(fp,"pair_coeff * %d soft%s 100 ${sigma_C_linker}\n",ANCHOR_POLYMER,Param->GPUFF);  // ANCHOR_POLYMER Anchor Polymer
	fprintf(fp,"pair_coeff * %d soft%s 100 ${sigma_C_linker}\n",SRD,Param->GPUFF);  // SRD Srd

	
		for (k = END_POLYMER; k < END_POLYMER + POL_END_TYPES; k++)
		{
			fprintf(fp, "pair_coeff * %d soft%s 0 ${sigma_C_linker}\n", k, Param->GPUFF); // Ends of the polymers
		}

		// Charge interactions between polymer ends
		for (k = END_POLYMER + END_POSITIVE_ACCEPTOR; k < END_POLYMER + END_POSITIVE_ACCEPTOR + 3; k++)
		{
			for (k1 = k; k1 < END_POLYMER + END_POSITIVE_ACCEPTOR + 3; k1++)
			{
				// fprintf(fp,"pair_coeff %d %d yukawa%s ${yukawa_strength_repul} \n",k,k1); //Ends Ends
				fprintf(fp, "pair_coeff %d %d soft%s 100 ${sigma_C_linker} \n", k, k1, Param->GPUFF); // Ends Ends
			}
			for (k1 = END_POLYMER + END_NEGATIVE_ACCEPTOR; k1 < END_POLYMER + END_NEGATIVE_ACCEPTOR + 3; k1++)
			{
				fprintf(fp, "pair_coeff %d %d yukawa%s ${yukawa_strength_attr} \n", k, k1, Param->GPUFF); // Ends Ends
			}
		}

		for (k = END_POLYMER + END_NEGATIVE_ACCEPTOR; k < END_POLYMER + END_NEGATIVE_ACCEPTOR + 3; k++)
		{
			for (k1 = k; k1 < END_POLYMER + END_NEGATIVE_ACCEPTOR + 3; k1++)
			{
				// fprintf(fp,"pair_coeff %d %d yukawa%s ${yukawa_strength_repul} \n",k,k1); //Ends Ends
				fprintf(fp, "pair_coeff %d %d soft%s 100 ${sigma_C_linker} \n", k, k1, Param->GPUFF); // Ends Ends
			}
		}
	
	
	fprintf(fp,"pair_modify shift yes\n");
	fprintf(fp,"\n\n");
	
	return;
}



void write_protein_coeff(struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp){
	int ii,jj;
	int i,j;
	double dx_CA_CA,dy_CA_CA,dz_CA_CA;
	double r2,r;
	
	ii=1;
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		jj=1;
		for(j=ATOM_CA;j<Protein->N;j+=NATOM){
			if(jj>=ii){
				if((((ii>jj+ATOM_IGNORE)||(ii<jj-ATOM_IGNORE))&&((jj>ii+ATOM_IGNORE)||(jj<ii-ATOM_IGNORE)))&&(jj!=ii)){
					dx_CA_CA=Protein->Atom[j].x-Protein->Atom[i].x;
					dy_CA_CA=Protein->Atom[j].y-Protein->Atom[i].y;
					dz_CA_CA=Protein->Atom[j].z-Protein->Atom[i].z;
					
					r2=dx_CA_CA*dx_CA_CA+dy_CA_CA*dy_CA_CA+dz_CA_CA*dz_CA_CA;
					if((r2<sist->ECA_Range2)){
						r=sqrt(r2);
						//fprintf(fp,"pair_coeff %d %d soft 1 ${bonds_NN}\n",ii,jj);
						//fprintf(fp,"pair_coeff %d %d sigmoidal ${GO} ${alpha} %lf %lf\n",ii,jj,sqrt(Protein->Average_r[ii-1]/Count_r2[ii-1]),sqrt(Protein->Average_r[ii-1]/Count_r2[ii-1])*2.0);
						if(r2<Protein->Average_r[ii-1]*Protein->Average_r[ii-1]){
							fprintf(fp,"pair_coeff %d %d sigmoidal ${GO} ${alpha} %lf %lf\n",ii,jj,r,r*2.0);
						}else{
							fprintf(fp,"pair_coeff %d %d sigmoidal ${GO2} ${alpha} %lf %lf\n",ii,jj,Protein->Average_r[ii-1],Protein->Average_r[ii-1]*2.0);
						}
					}else{
						//fprintf(fp,"pair_coeff %d %d soft%s 100 \n",ii,jj,Param->GPUFF);
						//fprintf(fp,"pair_coeff %d %d soft 1 ${bonds_NN}\n",ii,jj);
					}
				}else{
					//fprintf(fp,"pair_coeff %d %d soft%s 0 \n",ii,jj,Param->GPUFF);
					//fprintf(fp,"pair_coeff %d %d none\n",ii,jj);
				}
			}
			
			jj++;
		}
		fprintf(fp,"print \"Pair %d Ok\"\n",ii);
		ii++;
	}
	if(Param->N_Polymer>0){
		ii = 1;
		for (i = ATOM_CA; i < Protein->N; i += NATOM)
		{
			// HBonds with Polymer ENDS
			if ((Protein->Atom[i].prot == END_NEUTRAL_DONOR) || (Protein->Atom[i].prot == END_POSITIVE_DONOR) || (Protein->Atom[i].prot == END_NEGATIVE_DONOR))
			{
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEUTRAL_ACCEPTOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_POSITIVE_ACCEPTOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEGATIVE_ACCEPTOR + END_POLYMER, 3., 3 * 2.0);
			}
			if ((Protein->Atom[i].prot == END_NEUTRAL_ACCEPTOR) || (Protein->Atom[i].prot == END_POSITIVE_ACCEPTOR) || (Protein->Atom[i].prot == END_NEGATIVE_ACCEPTOR))
			{
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEUTRAL_DONOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_POSITIVE_DONOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEGATIVE_DONOR + END_POLYMER, 3., 3 * 2.0);
			}

			// Charges with Polymer ENDS
			if ((Protein->Atom[i].prot == END_POSITIVE_ACCEPTOR) || (Protein->Atom[i].prot == END_POSITIVE_DONOR) || (Protein->Atom[i].prot == END_POSITIVE_ACCEPTOR_DONOR))
			{
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEGATIVE_ACCEPTOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEGATIVE_DONOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_NEGATIVE_ACCEPTOR_DONOR + END_POLYMER, 3., 3 * 2.0);
			}
			if ((Protein->Atom[i].prot == END_NEGATIVE_ACCEPTOR) || (Protein->Atom[i].prot == END_NEGATIVE_DONOR) || (Protein->Atom[i].prot == END_NEGATIVE_ACCEPTOR_DONOR))
			{
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_POSITIVE_ACCEPTOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_POSITIVE_DONOR + END_POLYMER, 3., 3 * 2.0);
				fprintf(fp, "pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n", ii, END_POSITIVE_ACCEPTOR_DONOR + END_POLYMER, 3., 3 * 2.0);
			}
			/*if((Param->Geometry==SPHERE)&&(Param->mask_wall_flag==YES)&&(Protein->Mask_wall[i]==YES)) {
				fprintf(fp,"pair_coeff %d %d sigmoidal ${hb_strength} ${alpha} %lf %lf\n",ii,COLLOID,3.,3*2.0);
			}*/
			ii++;
		}
	}
	
	return;
	
}




void write_angles (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp){
	
	double p1[3],p2[3],p3[3],p4[3];
	double dx_CA_CA_bw,dy_CA_CA_bw,dz_CA_CA_bw;
	double dx_CA_CA_fw,dy_CA_CA_fw,dz_CA_CA_fw;
	double r_fw,r_bw;
	double angolo,phi;
	int Bond_index=1;
	int ii,i;
	
	fprintf(fp,"dihedral_style\t go\n\n");
	ii=1;
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		Protein->Atom[i].flag=0;
	}
	for(i=ATOM_CA;i<Protein->N;i+=NATOM){
		if(i+3*NATOM<Protein->N){
			if((Protein->Ignore[i]==0)&&(Protein->Ignore[i+NATOM]==0)&&(Protein->Ignore[i+2*NATOM]==0)&&(Protein->Ignore[i+3*NATOM]==0)){
				p1[0]=Protein->Atom[i].x;
				p1[1]=Protein->Atom[i].y;
				p1[2]=Protein->Atom[i].z;
				
				p2[0]=Protein->Atom[i+NATOM].x;
				p2[1]=Protein->Atom[i+NATOM].y;
				p2[2]=Protein->Atom[i+NATOM].z;
				
				p3[0]=Protein->Atom[i+2*NATOM].x;
				p3[1]=Protein->Atom[i+2*NATOM].y;
				p3[2]=Protein->Atom[i+2*NATOM].z;
				
				p4[0]=Protein->Atom[i+3*NATOM].x;
				p4[1]=Protein->Atom[i+3*NATOM].y;
				p4[2]=Protein->Atom[i+3*NATOM].z;
				
				phi=calc_dihedralf_angle_LAMMPS(p1,p2,p3,p4)*180./PI;
				
				
				if((fabs(fabs(phi)-45)<10)||(fabs(fabs(phi)-150)<30)) {
					
					fprintf(fp," dihedral_coeff %d ${torsional_1} %lf\n",Bond_index,phi);
					
					
				}else{ // Random
					
					fprintf(fp," dihedral_coeff %d ${torsional_2} %lf\n",Bond_index,phi);
					Protein->Atom[i].flag=1;
					Protein->Atom[i+NATOM].flag=1;
					Protein->Atom[i+2*NATOM].flag=1;
					Protein->Atom[i+3*NATOM].flag=1;
				}
				Bond_index++;
			}
		}else{
			// To avoid double counting the last three beads are designed to have zero energy contribution
		}
	}	
	
	fprintf(fp,"\n\n");
	
	fprintf(fp,"angle_style\t harmonic\n\n");
	ii=1;
	
	for(i=ATOM_CA;i<Protein->N-2*NATOM;i+=NATOM){
		
		if((Protein->Ignore[i]==0)&&(Protein->Ignore[i+NATOM]==0)&&(Protein->Ignore[i+2*NATOM]==0)){
			dx_CA_CA_bw=Protein->Atom[i].x-Protein->Atom[i+NATOM].x;
			dy_CA_CA_bw=Protein->Atom[i].y-Protein->Atom[i+NATOM].y;
			dz_CA_CA_bw=Protein->Atom[i].z-Protein->Atom[i+NATOM].z;
			
			dx_CA_CA_fw=Protein->Atom[i+2*NATOM].x-Protein->Atom[i+NATOM].x;
			dy_CA_CA_fw=Protein->Atom[i+2*NATOM].y-Protein->Atom[i+NATOM].y;
			dz_CA_CA_fw=Protein->Atom[i+2*NATOM].z-Protein->Atom[i+NATOM].z;
			
			r_fw=sqrt(dx_CA_CA_fw*dx_CA_CA_fw+dy_CA_CA_fw*dy_CA_CA_fw+dz_CA_CA_fw*dz_CA_CA_fw);
			r_bw=sqrt(dx_CA_CA_bw*dx_CA_CA_bw+dy_CA_CA_bw*dy_CA_CA_bw+dz_CA_CA_bw*dz_CA_CA_bw);
			angolo=acos((dx_CA_CA_fw*dx_CA_CA_bw+dy_CA_CA_fw*dy_CA_CA_bw+dz_CA_CA_fw*dz_CA_CA_bw)/(r_fw*r_bw))*180./PI;
			if(Protein->Atom[i+NATOM].flag==0){
				fprintf(fp,"angle_coeff %d ${angle_k} %15.10lf\n",ii,angolo);
			}else{
				fprintf(fp,"angle_coeff %d ${angle_k2} %15.10lf\n",ii,angolo);
			}
			
			
			ii++;
		}
	}
	return;
}
void write_fix (struct Interaction_Param *Param, FILE *fp){
	
	if(Param->Geometry==CYLINDER) {
		fprintf(fp,"fix wall_anchor anchor wall/region porewall lj126 10 ${sigma_anchor} ${cutoff_C}\n"); // Repulsive  wall to trap anchors
		if(Param->brush_flag==YES) fprintf(fp,"fix wall_polymer polymer wall/region poreinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
		if(Param->protein_flag==YES) fprintf(fp,"fix wall_prot_rep protein wall/region poreinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
	
		if(Param->mask_wall_flag==YES) fprintf(fp,"fix wall_prot_attr wall_int wall/region poreinner lj126 %lf ${sigma_C} ${range_mol}\n",Param->Wall_Attr); // Attractive special interactions to the wall
	
	}

	if(Param->Geometry==SLIT) {
		fprintf(fp,"fix wall_anchor anchor wall/region blockwall lj126 10 ${sigma_anchor} ${cutoff_C}\n"); // Repulsive  wall to trap anchors
		if(Param->brush_flag==YES) fprintf(fp,"fix wall_polymer polymer wall/region blockinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
		if(Param->protein_flag==YES) fprintf(fp,"fix wall_prot_rep protein wall/region blockinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
	
		if(Param->mask_wall_flag==YES) fprintf(fp,"fix wall_prot_attr wall_int wall/region blockinner lj126 %lf ${sigma_C} ${range_mol}\n",Param->Wall_Attr); // Attractive special interactions to the wall
	
	}
	if(Param->Geometry==CHANNEL) {
		fprintf(fp,"fix wall_anchor anchor wall/region blockwall lj126 10 ${sigma_anchor} ${cutoff_C}\n"); // Repulsive  wall to trap anchors
		if(Param->brush_flag==YES) fprintf(fp,"fix wall_polymer polymer wall/region blockinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
		if(Param->protein_flag==YES) fprintf(fp,"fix wall_prot_rep protein wall/region blockinner lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
	
		if(Param->mask_wall_flag==YES) fprintf(fp,"fix wall_prot_attr wall_int wall/region blockinner lj126 %lf ${sigma_C} ${range_mol}\n",Param->Wall_Attr); // Attractive special interactions to the wall
	
	}
	
	fprintf(fp, "velocity      forcegroup create 1  ${seed}\n");
  if(Param->brush_flag==YES) fprintf(fp, "velocity anchor_polymer zero linear\n");
  if(Param->brush_flag==YES) fprintf(fp, "fix freeze_anchor anchor_polymer setforce 0.0 0.0 0.0\n");
	
	return;
}
void write_fix_prot (struct Interaction_Param *Param, FILE *fp){
	
	//We place the colloidal wall here becasue we need to remove the clashes with the colloid
	if(Param->Geometry==SPHERE) { // We write it sperately for the spheres becasue I notices that sometimes there are clashes with the colloid and is better to minimize first		
		if(Param->protein_flag==YES) fprintf(fp,"fix wall_prot_rep protein wall/region colloidouter lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall	
		if(Param->mask_wall_flag==YES) fprintf(fp,"fix wall_prot_attr wall_int wall/region colloidouter lj126 %lf ${sigma_C} ${range_mol}\n",Param->Wall_Attr); // Attractive special interactions to the wall
		//fprintf(fp,"fix wall_anchor anchor wall/region colloidouter lj126 10 ${sigma_anchor} ${cutoff_C}\n"); // Repulsive  wall to trap anchors
		if(Param->brush_flag==YES) fprintf(fp,"fix wall_polymer polymer wall/region colloidouter lj126 10.0 ${sigma_C} ${cutoff_C}\n"); // Repulsive to the wall
	}
	
	return;
}
void write_harmonicbonds (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins *Protein, FILE *fp){
	
	int ii,jj;
	int i,j;
	double dx_CA_CA,dy_CA_CA,dz_CA_CA;
	double r2;
	if(Param->protein_flag==YES){
	fprintf(fp,"bond_style hybrid gromos harmonic\n\n");
	fprintf(fp,"special_bonds lj 0 0 0\n");
	//fprintf(fp,"bond_coeff 1 ${prot_k} ${sigma}\n");	
	//fprintf(fp,"\n\n");
	
	ii=0;
	for(i=ATOM_CA;i<Protein->N-NATOM;i+=NATOM){
		
		if(Protein->Ignore[i]==0){
			dx_CA_CA=Protein->Atom[i+NATOM].x-Protein->Atom[i].x;
			dy_CA_CA=Protein->Atom[i+NATOM].y-Protein->Atom[i].y;
			dz_CA_CA=Protein->Atom[i+NATOM].z-Protein->Atom[i].z;
			
			r2=dx_CA_CA*dx_CA_CA+dy_CA_CA*dy_CA_CA+dz_CA_CA*dz_CA_CA;
			fprintf(fp,"bond_coeff %d gromos ${prot_k} %15.10lf # FENE Extra params 2.0  %15.10lf \n",Protein->bondtype_head+ii,sqrt(r2),sqrt(r2));
			
			
			ii++;
		}
	}
	fprintf(fp,"bond_coeff %d harmonic ${linker_k} ${sigma_linker}\n",sist->polymers_bonds_type);		
}else{
	fprintf(fp,"bond_style harmonic\n\n");
	fprintf(fp,"bond_coeff %d ${linker_k} ${sigma_linker}\n",sist->polymers_bonds_type);		
}
	
	fprintf(fp,"\n\n");
	return;
}


void write_1_minimize (FILE *fp, struct Interaction_Param *Param){
	
	fprintf(fp,"#Group size write_1_minimize\n");
	if(Param->Simul_type==ANCHOR){
		fprintf(fp,"variable natoms equal \"count(free)\"\n");
		fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"variable natoms equal \"count(docked)\"\n");
			fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
		}
	}
	
	if(Param->N_Polymer>0){
    fprintf(fp, "pair_style soft%s ${sigma_C}\n",Param->GPUFF);
		fprintf(fp, "velocity      polymer create 10.0  ${seed}\n");
    fprintf(fp, "pair_coeff * * 1.0\n");
		if(Param->Geometry==SPHERE)  fprintf(fp,"pair_coeff * %d 1000 %lf\n",COLLOID,Param->Coll_R+1.5);  // COLLOID LInker 
    fprintf(fp, "variable prefactor equal ramp(1,100)\n");
    fprintf(fp, "pair_modify shift yes\n");
    fprintf(fp, "fix 4 forcegroup adapt 1 pair soft a * * v_prefactor\n");
    fprintf(fp, "fix 2 polymer nve\n");
    fprintf(fp, "fix 3 polymer langevin 10 10 1 ${seed} zero yes\n");
    fprintf(fp, "fix 5 polymer momentum 10 linear 1 1 1 angular rescale\n");
    fprintf(fp, "dump 1 all xtc 10 1_minimize.lammpstrj.xtc\n");
	 fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
    if(Param->protein_flag==YES) fprintf(fp, "dump 2 dock_group custom 1000 1_minimize.lammpstrj id  type  x y z  ix iy iz\n");
    fprintf(fp, "thermo          1000\n");
    fprintf(fp, "timestep        0.01\n");
		fprintf(fp, "restart 1000 restart_1_minimize.1 restart_1_minimize.2\n");
    fprintf(fp, "run 20000\n");
    fprintf(fp, "write_restart restart_1_minimize.3\n");
    fprintf(fp, "undump 1\n");
   	if(Param->protein_flag==YES)  fprintf(fp, "undump 2\n");
		fprintf(fp, "unfix 2\n");
		fprintf(fp, "unfix 3\n");
		fprintf(fp, "unfix 5\n");
    fprintf(fp, "unfix 4\n");

	}
		
		

	
	return;
}
void write_2_minimize (FILE *fp, struct Interaction_Param *Param){	
	
	
	fprintf(fp,"#Group size write_2_minimize\n");
	if(Param->Simul_type==ANCHOR){
			fprintf(fp,"variable natoms equal \"count(free)\"\n");
			fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"variable natoms equal \"count(docked)\"\n");
			fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
		}
	}



	fprintf(fp, "dump 1 all xtc 1000 2_unfolding.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
  if(Param->protein_flag==YES) fprintf(fp, "dump 2 dock_group custom 1000 2_unfolding.lammpstrj id  type  x y z  ix iy iz\n");
  fprintf(fp, "velocity      forcegroup create 10.0  ${seed}\n");

  fprintf(fp, "pair_style soft%s ${sigma_C}\n",Param->GPUFF);

  fprintf(fp, "pair_coeff * * 100\n");
	if(Param->Geometry==SPHERE)  fprintf(fp,"pair_coeff * %d 1000 %lf\n",COLLOID,Param->Coll_R+1.5);  // COLLOID LInker 
  fprintf(fp, "pair_modify shift yes\n");
  fprintf(fp, "fix 2 forcegroup nve\n");
  fprintf(fp, "fix 3 forcegroup langevin 10.0 10.0 1 ${seed} zero yes\n");
  fprintf(fp, "fix 5 forcegroup momentum 10 linear 1 1 1 angular rescale\n");
  fprintf(fp, "thermo          1000\n");
  fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_2_unfolding.1 restart_2_unfolding.2\n");
  fprintf(fp, "run 3000\n");
  fprintf(fp, "write_restart restart_2_unfolding.3\n");
  fprintf(fp, "undump 1\n");
  if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
	fprintf(fp, "unfix 5\n");
  fprintf(fp, "velocity      forcegroup create 1  ${seed}\n");
	fprintf(fp, "variable linker_k equal ${eps_h}*5\n");
	
	return;
}
void write_3_minimize (FILE *fp, struct Interaction_Param *Param){	
	
	
	fprintf(fp,"#Group size write_3_minimize\n");
			if(Param->Simul_type==ANCHOR){
	fprintf(fp,"variable natoms equal \"count(free)\"\n");
	fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
	if(Param->mask_flag==YES) {
	fprintf(fp,"variable natoms equal \"count(docked)\"\n");
	fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
}
}
	if(Param->protein_flag==YES){
		fprintf(fp, "compute cc_protein all pair sigmoidal evdwl\n");
	}
	fprintf(fp, "\n");
	if(Param->protein_flag==YES){
		fprintf(fp, "thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cc_protein\n");
	}else{
		fprintf(fp, "thermo_style custom step temp ebond eangle edihed epair pe ke etotal\n");
	}
	fprintf(fp, "thermo 1000\n");
	fprintf(fp, "\n");
	fprintf(fp, "dump 1 all xtc   100  3_folding.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp, "dump 2 dock_group custom   100  3_folding.lammpstrj  id  type  x y z  ix iy iz\n");
	fprintf(fp, "fix 2 forcegroup nve\n");
	fprintf(fp, "fix 3 forcegroup langevin 1 1 1 ${seed} zero yes\n");
	fprintf(fp, "fix 5 forcegroup momentum 1 linear 1 1 1 angular rescale\n");
	fprintf(fp, "thermo          5000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_3_folding.1 restart_3_folding.2\n");
	fprintf(fp, "run 500000\n");
	fprintf(fp, "write_restart restart_3_folding.3\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
	fprintf(fp, "unfix 5\n");
	
	return;
}
void write_4_docking (struct GlobalSistem *sist, FILE *fp, struct Interaction_Param *Param, struct Proteins *Protein){
	
	double force_x,force_y,force_z,module;
	
	switch(Param->Geometry){
		case SPHERE:
		protein_gmass(sist,Protein);
		fprintf(fp,"group temp id %s",Protein->group_id);
		force_x=sist->Half_Box_x-Protein->gmass->x; 
		force_y=sist->Half_Box_y-Protein->gmass->y; 
		force_z=sist->Half_Box_z-Protein->gmass->z; 
		module=sqrt(force_x*force_x+force_y*force_y+force_z*force_z);
		fprintf(fp,"\n");	
		fprintf(fp,"#Protein gmass %lf %lf %lf\n",Protein->gmass->x,Protein->gmass->y,Protein->gmass->z);	
		fprintf(fp,"#Sphere position %lf %lf %lf\n",sist->Half_Box_x,sist->Half_Box_y,sist->Half_Box_z);	
				
		fprintf(fp, "fix dock temp addforce %lf %lf %lf\n",force_x*10/module,force_y*10/module,force_z*10/module);
		fprintf(fp, "fix dock anchor addforce %lf %lf %lf\n",force_x*15/module,force_y*15/module,force_z*15/module);
		break;
		case CYLINDER:	
			fprintf(fp,"group temp union dock_group\n");
			fprintf(fp,"fix dock temp addforce 0.0 0.0 -0.5\n");
			fprintf(fp,"fix dock anchor addforce 0.0 0.0 -1.0\n");
			
		break;
		case SLIT:	
			fprintf(fp,"group temp union dock_group\n");
			fprintf(fp,"fix dock temp addforce 0.0 0.0 -0.5\n");
			fprintf(fp,"fix dock anchor addforce 0.0 0.0 -1.0\n");
			
		break;
		case CHANNEL:	
			fprintf(fp,"group temp union dock_group\n");
			fprintf(fp,"fix dock temp addforce 0.0 0.0 -0.5\n");
			fprintf(fp,"fix dock anchor addforce 0.0 0.0 -1.0\n");
			
		break;
	}
	fprintf(fp, "\n");
	
	fprintf(fp,"#Group size write_4a_minimize\n");
	if(Param->Simul_type==ANCHOR){
		fprintf(fp,"variable natoms equal \"count(free)\"\n");
		fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
		if(Param->mask_flag==YES) {
			fprintf(fp,"variable natoms equal \"count(docked)\"\n");
			fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
		}
	}
	fprintf(fp, "dump 1 all xtc   1000  4a_docking.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp, "dump 2 temp custom   1000  4a_docking.lammpstrj  id  type  x y z  ix iy iz\n");
	fprintf(fp, "dump_modify 2 append yes\n");  // Ensure appending for the lammpstrj file
	fprintf(fp, "fix 2 free nve\n");
	fprintf(fp, "fix 3 free langevin 1 1 1 ${seed} zero yes\n");
	//fprintf(fp, "fix 5 free momentum 1 linear 1 1 1 angular rescale\n");
	fprintf(fp, "thermo          1000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_4a_docking.1 restart_4a_docking.2\n");
	fprintf(fp, "run 200000\n");
	fprintf(fp, "write_restart restart_4a_docking.3\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
	/*switch(Param->Geometry){
		case SPHERE:
		break;
		case CYLINDER:
		break;
		case SLIT:
		case CHANNEL:
		fprintf(fp, "set group docked z %lf\n",sist->MOL_SHIFT-Param->Polymer_Size*1.1/2);
		break;
	}*/

	//fprintf(fp, "unfix 5\n");
	
	/*switch(Param->Geometry){
		case SPHERE:		
		fprintf(fp, "fix dock temp addforce %lf %lf %lf\n",force_x/module,force_y/module,force_z/module);
		break;
		case CYLINDER:	
			fprintf(fp, "fix dock temp addforce 0.0 0.0 -15.0\n");
		break;
	}
	
	
	fprintf(fp, "\n");
	
	fprintf(fp,"#Group size write_4b_minimize\n");
	fprintf(fp,"variable natoms equal \"count(free)\"\n");
	fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
	fprintf(fp,"variable natoms equal \"count(docked)\"\n");
	fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
	
	fprintf(fp, "dump 1 all xtc   1000  4b_docking.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp, "dump 2 temp custom   1000  4b_docking.lammpstrj  id  type  x y z  ix iy iz\n");
	fprintf(fp, "dump_modify 2 append yes\n");  // Ensure appending for the lammpstrj file
	fprintf(fp, "fix 2 free nve\n");
	fprintf(fp, "fix 3 free langevin 1 1 1 ${seed} zero yes\n");
	//fprintf(fp, "fix 5 free momentum 1 linear 1 1 1 angular rescale\n");
	fprintf(fp, "thermo          1000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "run 200000\n");
	fprintf(fp, "write_restart restart_4b_docking\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
	//fprintf(fp, "unfix 5\n");*/
	
	fprintf(fp, "unfix dock\n");

	fprintf(fp, "group temp delete\n");
		
	return;
	
}

void write_5_anchored (FILE *fp,struct Interaction_Param *Param){
	
	
	fprintf(fp,"#Group size write_5_anchored\n");
			if(Param->Simul_type==ANCHOR){
	fprintf(fp,"variable natoms equal \"count(free)\"\n");
	fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
	if(Param->mask_flag==YES) {
	fprintf(fp,"variable natoms equal \"count(docked)\"\n");
	fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
}
}
	fprintf(fp,"dump 1 all xtc  1000  5_anchored.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp,"dump 2 dock_group custom  1000  5_anchored.lammpstrj  id  type  x y z  ix iy iz\n");
	if(Param->mask_flag==YES) fprintf(fp,"velocity anchor zero linear\n");
	if(Param->Geometry==CYLINDER) if(Param->mask_flag==YES) fprintf(fp,"fix 6 docked setforce 0.0 0.0 0.0 region anchorwall\n");
	if(Param->Geometry==SLIT) if(Param->mask_flag==YES) fprintf(fp,"fix 6 docked setforce 0.0 0.0 0.0 region anchorwall\n");
	if(Param->Geometry==CHANNEL) if(Param->mask_flag==YES) fprintf(fp,"fix 6 docked setforce 0.0 0.0 0.0 region anchorwall\n");

	fprintf(fp, "fix 2 free nve\n");
	fprintf(fp, "fix 3 free langevin 1 1 1 ${seed} zero yes\n");
	
	fprintf(fp, "thermo          1000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_5_anchored.1 restart_5_anchored.2\n");
	fprintf(fp, "run 200000\n");
	fprintf(fp, "write_restart restart_5_anchored.3\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
		
	return;
	
}
void write_fixbrooks (struct GlobalSistem *sist, struct Interaction_Param *Param, struct Proteins **Protein_pointers, FILE *fp){
	int i,j,m;
	int res_flag;
	
	for(i=1;i<21;i++){
		res_flag=0;
		for(m=0;m<Param->N_Proteins;m++){
			for(j=ATOM_CA;j<Protein_pointers[m]->N;j+=NATOM){
				if(Protein_pointers[m]->Atom[j].res==i) res_flag=1;		
			}
		}
		if(res_flag==1){
			if(Param->Geometry==CYLINDER)  fprintf(fp,"fix wall_%s %s wall/region poreinner brooks %lf ${chi_s} %lf ${sigma_C} ${cutoff_C}\n",Protein_pointers[0]->Amminoacids[i],Protein_pointers[0]->Amminoacids[i],Param->Wall_Attr_Gen,sist->chi_p[i]); // Hydrophobic to the wall
			if(Param->Geometry==SLIT)  fprintf(fp,"fix wall_%s %s wall/region blockinner brooks %lf ${chi_s} %lf ${sigma_C} ${cutoff_C}\n",Protein_pointers[0]->Amminoacids[i],Protein_pointers[0]->Amminoacids[i],Param->Wall_Attr_Gen,sist->chi_p[i]); // Hydrophobic to the wall
			if(Param->Geometry==CHANNEL)  fprintf(fp,"fix wall_%s %s wall/region blockinner brooks %lf ${chi_s} %lf ${sigma_C} ${cutoff_C}\n",Protein_pointers[0]->Amminoacids[i],Protein_pointers[0]->Amminoacids[i],Param->Wall_Attr_Gen,sist->chi_p[i]); // Hydrophobic to the wall
			if(Param->Geometry==SPHERE)    fprintf(fp,"fix wall_%s %s wall/region colloidouter brooks %lf ${chi_s} %lf ${sigma_C} ${cutoff_C}\n",Protein_pointers[0]->Amminoacids[i],Protein_pointers[0]->Amminoacids[i],Param->Wall_Attr_Gen,sist->chi_p[i]);
		}
	}
	
	return;
	
}
void write_6_bound_init_T (FILE *fp, struct Interaction_Param *Param){
	
	
	
	fprintf(fp,"#Group size write_6_bound_init_T\n");
			if(Param->Simul_type==ANCHOR){
	fprintf(fp,"variable natoms equal \"count(free)\"\n");
	fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
	if(Param->mask_flag==YES) {
	fprintf(fp,"variable natoms equal \"count(docked)\"\n");
	fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
}
}
	
	fprintf(fp, "dump 1 all xtc   1000  6_bound_lowT.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp, "dump 2 dock_group custom   1000  6_bound_lowT.lammpstrj  id  type  x y z  ix iy iz\n");
	
	
	
	fprintf(fp, "fix 2 free nve\n");
	fprintf(fp, "fix 3 free langevin 1 ${init_temp} ${final_temp} ${seed} zero yes\n");
	
	
	fprintf(fp, "thermo          1000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_6_bound_lowT.1 restart_6_bound_lowT.2\n");
	fprintf(fp, "run 100000\n");
	fprintf(fp, "write_restart restart_6_bound_lowT.3\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");

	
	
	
	return;
}


void write_6_bound_final_T (FILE *fp, struct Interaction_Param *Param){

	
	fprintf(fp,"#Group size write_6_bound_final_T\n");
			if(Param->Simul_type==ANCHOR){
	fprintf(fp,"variable natoms equal \"count(free)\"\n");
	fprintf(fp,"print \"Number of atoms in the group free : ${natoms}\"\n");
	if(Param->mask_flag==YES) {
	fprintf(fp,"variable natoms equal \"count(docked)\"\n");
	fprintf(fp,"print \"Number of atoms in the group docked : ${natoms}\"\n");
}
}
	
	fprintf(fp, "dump 1 all xtc   1000  6_bound_highT.lammpstrj.xtc\n");
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	if(Param->protein_flag==YES) fprintf(fp, "dump 2 dock_group custom   1000  6_bound_highT.lammpstrj  id  type  x y z  ix iy iz\n");
	
	
	fprintf(fp, "fix 2 free nve\n");
	fprintf(fp, "fix 3 free langevin ${final_temp} ${init_flow_temp} 1 ${seed} zero yes\n");
	
	
	fprintf(fp, "thermo          1000\n");
	fprintf(fp, "timestep        0.002\n");
	fprintf(fp, "restart 1000 restart_6_bound_highT.1 restart_6_bound_highT.2\n");
	fprintf(fp, "run 100000\n");
	fprintf(fp, "write_restart restart_6_bound_highT.3\n");
	fprintf(fp, "undump 1\n");
	if(Param->protein_flag==YES) fprintf(fp, "undump 2\n");
	fprintf(fp, "unfix 2\n");
	fprintf(fp, "unfix 3\n");
	
	return;
}

void write_7_srd (FILE *fp,struct Interaction_Param *Param, int restart_flag, double flow_temp, int srd_step){
	printf("write_7_srd step=%d/%d\n",srd_step,Param->steps_flow_temp);
	
	fprintf(fp,"atom_modify first forcegroup\n");
	fprintf(fp,"neigh_modify include forcegroup\n");
	fprintf(fp,"variable flow_temp 	equal %lf\n",flow_temp);
	
	if(restart_flag==NO) fprintf(fp,"lattice        sc 0.5\n");
	
	if(Param->Simul_type==ANCHOR){
		fprintf(fp,"variable natoms equal \"count(free)\"\n");
		fprintf(fp,"print \"Number of atoms in the groupfree : ${natoms}\"\n");
	}
	if(Param->Geometry==CYLINDER){		
		if(restart_flag==NO) fprintf(fp,"create_atoms   %d region poresrd\n",SRD);
		fprintf(fp,"group solvent type %d\n",SRD);
		//fprintf(fp,"group moving union solvent polymer protein\n");
		fprintf(fp,"group  free delete\n");
		fprintf(fp,"group  free region poreinner\n");
		fprintf(fp,"group dock_group delete\n");
		fprintf(fp,"group temp type %d\n",COLLOID);
		fprintf(fp,"group dock_group union protein anchor polymer temp\n");	
		fprintf(fp,"group temp delete\n");
	}
	if(Param->Geometry==SLIT){		
		if(restart_flag==NO) fprintf(fp,"create_atoms   %d region blocksrd\n",SRD);
		fprintf(fp,"group solvent type %d\n",SRD);
		//fprintf(fp,"group moving union solvent polymer protein\n");
		fprintf(fp,"group  free delete\n");
		fprintf(fp,"group  free region blockinner\n");
		fprintf(fp,"group dock_group delete\n");
		fprintf(fp,"group temp type %d\n",COLLOID);
		fprintf(fp,"group dock_group union protein anchor polymer temp\n");	
		fprintf(fp,"group temp delete\n");
	}
	if(Param->Geometry==CHANNEL){		
		if(restart_flag==NO) fprintf(fp,"create_atoms   %d region blocksrd\n",SRD);
		fprintf(fp,"group solvent type %d\n",SRD);
		//fprintf(fp,"group moving union solvent polymer protein\n");
		fprintf(fp,"group  free delete\n");
		fprintf(fp,"group  free region blockinner\n");
		fprintf(fp,"group dock_group delete\n");
		fprintf(fp,"group temp type %d\n",COLLOID);
		fprintf(fp,"group dock_group union protein anchor polymer temp\n");	
		fprintf(fp,"group temp delete\n");
	}
	if(Param->Geometry==SPHERE){		
		if(restart_flag==NO) fprintf(fp,"create_atoms   %d box\n",SRD);
		fprintf(fp,"group solvent type %d\n",SRD);
		//fprintf(fp,"group moving union solvent polymer protein\n");
		fprintf(fp,"group  free delete\n");
		fprintf(fp,"group  free union forcegroup solvent\n");
		fprintf(fp,"group dock_group delete\n");
		fprintf(fp,"group temp type %d\n",COLLOID);
		fprintf(fp,"group dock_group union protein anchor polymer temp\n");	
		fprintf(fp,"group temp delete\n");
		
	}
	if(Param->Simul_type==ANCHOR){
		fprintf(fp,"variable natoms equal \"count(free)\"\n");
		fprintf(fp,"print \"Number of atoms in the groupfree : ${natoms}\"\n");
	}
	fprintf(fp, "#initial velocity of protein can be changed\n");
	fprintf(fp, "velocity    solvent create ${flow_temp} ${seed}\n");
	fprintf(fp, "#velocity    free set NULL NULL NULL units box\n\n");

	fprintf(fp, "timestep    0.002\n");
	fprintf(fp, "fix             s1 free srd 90 NULL ${flow_temp} 2.0 49894 rescale no tstat yes shift yes 1234 beadgroup forcegroup  fflow ${pressure} ${pressure} rcyl ${inradius}    virtual 5 alpha 120\n");
	
	fprintf(fp, "#just for orientation\n");
	if(Param->protein_flag==YES) fprintf(fp, "compute cm1 protein com\n");
	if(Param->protein_flag==YES) fprintf(fp, "compute rg protein gyration\n");
	fprintf(fp, "#control variables\n");
	if(Param->protein_flag==YES){
		fprintf(fp, "thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cm1[1] c_rg\n");
	}else{
		fprintf(fp, "thermo_style custom step temp ebond eangle edihed epair pe ke etotal\n");
	}
	fprintf(fp, "thermo        1000\n");
	fprintf(fp, "#write out gyration tensor\n");
	if(Param->protein_flag==YES) fprintf(fp, "fix            r1 protein ave/time 1000 1 1000 c_rg[1] c_rg[2] c_rg[3] c_rg[4] c_rg[5] c_rg[6] file tensor.data\n");
	fprintf(fp, "#write out flow profile if needed\n");
	fprintf(fp, "# Set a variable for x spacing\n");
  	fprintf(fp, "variable xspacing equal 10.0\n");
  	fprintf(fp, "# Define a compute that chunks atoms into 3D bins along x, y, z axes\n");
  	fprintf(fp, "compute cc_solvent solvent chunk/atom bin/3d x 0.0 ${xspacing} y 0.0 0.2 z 0.0 0.2 units box\n");
  	fprintf(fp, "# Set a fix that averages chunk properties over a series of timesteps and outputs to file\n");
  	fprintf(fp, "fix v1 solvent ave/chunk 100 1000 100000 cc_solvent vx vy vz temp norm all file vel_srd_f1t\n");
	fprintf(fp, "dump 1 dock_group xtc    5000  7_flow.lammpstrj_%d.xtc\n",srd_step);
	fprintf(fp, "dump_modify 1 append yes unwrap yes pbc yes\n");  // Ensure appending for the xtc file
	fprintf(fp, "restart 1000 restart_7_flow.1 restart_7_flow.2\n");
	if((srd_step==Param->steps_flow_temp+1)||(srd_step==Param->steps_flow_temp+2)){
		fprintf(fp, "run            10000000\n");
	}else{
		fprintf(fp, "run            10000\n");
	}
	fprintf(fp, "write_restart restart_7_flow.3\n");
	
	return;
	
}



void write_tcl_header (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp, double scale){
	
	
	
	fprintf(fp,"package require text\n");
	fprintf(fp,"set molname [topo readlammpsdata %s-out.atomparam full]\n",Param->prefix);
	fprintf(fp,"pbc box\n");
	switch(Param->Geometry){
		case CYLINDER:
		fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->radius*scale);	
		fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->inradius*scale);
		break;
	}
	
	return;
	
}
void write_tcl_protein (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp, double scale){
	
	fprintf(fp,"package require text\n");
	fprintf(fp,"set molname [topo readlammpsdata %s-out-protein.topology full]\n",Param->prefix);
	fprintf(fp,"pbc box\n");
	switch(Param->Geometry){
		case CYLINDER:
	fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->radius*scale);	
	fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->inradius*scale);
	break;
	}
	fprintf(fp, "mol addfile \"1_minimize.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"2_unfolding.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"3_folding.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"4a_docking.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"4b_docking.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"5_anchored.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"6_bound_lowT.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"6_bound_highT.lammpstrj\" type lammpstrj waitfor all molid $molname\n");
return;
}
void write_tcl_all (struct GlobalSistem *sist,struct Interaction_Param *Param, FILE *fp, double scale){
	
	fprintf(fp,"package require text\n");
	fprintf(fp,"set molname [topo readlammpsdata %s-out-xtc-all.topology full]\n",Param->prefix);
	fprintf(fp,"pbc box\n");
		switch(Param->Geometry){
		case CYLINDER:
	fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->radius*scale);	
	fprintf(fp,"draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30 \n",sist->radius*scale,sist->radius*scale,sist->box_x*scale,sist->radius*scale,sist->radius*scale,sist->inradius*scale);
	break;
	}
	fprintf(fp, "mol addfile \"1_minimize.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"2_unfolding.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"3_folding.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"4a_docking.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"4b_docking.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"5_anchored.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"6_bound_lowT.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");
	fprintf(fp, "mol addfile \"6_bound_highT.lammpstrj.xtc\" type xtc waitfor all molid $molname\n");

	return;
}
void write_tcl_flow(struct GlobalSistem *sist, struct Interaction_Param *Param, FILE *fp, double scale) {
    // Load necessary package
    fprintf(fp, "package require text\n\n");

    // Load the topology file
    fprintf(fp, "# Load the topology file\n");
    fprintf(fp, "set molname [topo readlammpsdata \"%s-out-xtc-flow.topology\" full]\n", Param->prefix);
    fprintf(fp, "pbc box\n\n");

    // Handle different geometries (only CYLINDER for now)
    switch (Param->Geometry) {
        case CYLINDER:
            fprintf(fp, "# Draw cylinders for visualization\n");
            fprintf(fp, "draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30\n", 
                    sist->radius * scale, sist->radius * scale, sist->box_x * scale, 
                    sist->radius * scale, sist->radius * scale, sist->radius * scale);
            fprintf(fp, "draw cylinder {0 %lf %lf} {%lf %lf %lf} radius %lf resolution 30\n", 
                    sist->radius * scale, sist->radius * scale, sist->box_x * scale, 
                    sist->radius * scale, sist->radius * scale, sist->inradius * scale);
            break;
    }

    // Retrieve and sort XTC files numerically
    fprintf(fp, "\n# Load and sort XTC files numerically\n");
    fprintf(fp, "set xtc_files [glob -nocomplain 7_flow.lammpstrj_*.xtc]\n\n");

    // Define a sorting function for numerical order
    fprintf(fp, "proc num_sort {a b} {\n");
    fprintf(fp, "    regexp {lammpstrj_(\\d+)\\.xtc$} $a -> numA\n");
    fprintf(fp, "    regexp {lammpstrj_(\\d+)\\.xtc$} $b -> numB\n");
    fprintf(fp, "    if {![string is integer -strict $numA] || ![string is integer -strict $numB]} {\n");
    fprintf(fp, "        return [string compare $a $b]  ;# Fallback to lexicographic if extraction fails\n");
    fprintf(fp, "    }\n");
    fprintf(fp, "    return [expr {$numA - $numB}]\n");
    fprintf(fp, "}\n\n");

    // Sort the files numerically
    fprintf(fp, "set sorted_xtc_files [lsort -command num_sort $xtc_files]\n\n");

    // Ensure files exist before processing
    fprintf(fp, "if {[llength $sorted_xtc_files] == 0} {\n");
    fprintf(fp, "    puts \"No matching XTC files found.\"\n");
    fprintf(fp, "} else {\n");

    // Load files in correct order
    fprintf(fp, "    foreach xtc_file $sorted_xtc_files {\n");
    fprintf(fp, "        puts \"Loading file: $xtc_file\"\n");
    fprintf(fp, "        mol addfile $xtc_file type xtc waitfor all molid $molname\n");
    fprintf(fp, "    }\n");
    fprintf(fp, "    puts \"All matching XTC files loaded successfully.\"\n");
    fprintf(fp, "}\n");

    return;
}
