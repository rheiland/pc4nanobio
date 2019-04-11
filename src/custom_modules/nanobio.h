/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#ifndef nanobio__h_
#define nanobio__h_
 

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

double Hill_function( double input, double Hill_power , double EC_50 );
// output = (1-param)*input1 + param*input2 ( 0 <= param <= 1 )
void interpolate( double* pOutput , double input1, double input2, double param ); 
void interpolate( std::vector<double>* pOutput , std::vector<double>& input1, std::vector<double> input2, double param );

class Nanoparticle_Definition
{
 private:
 public:
	std::string name = "NP-preform"; 
	std::string units = "dimensionless";
	int ID = 0; 
	int substrate_ID = 4; 
	
	// tissue pharmacokinetics 
	double diffusion_coefficient = 1e3;
	double survival_lifetime = 90000; // 1500 hours 
	double ECM_binding_rate = 0.0; 
	double ECM_unbinding_rate = 0.0; 
	double ECM_saturation_concentration = 0.1; 
	
	double degradation_rate = 1e-5; // ~ 1500 hour survival time ; set to 1/survival_lifetime
	double blood_concentration = 0.0; // 
	double blood_tissue_transfer_rate = 0.0; // 
	
	void read_from_pugixml( pugi::xml_node parent_node ); // done 
		// missing: blood stuff 
};

class Basic_Phenotype
{
 public:
	std::string name = "live"; 
 
	// cycle 
	double max_birth_rate = 9.26e-4; // 1/min  // used 

	double o2_proliferation_saturation = 38; // value at which extra o2 does not increase proliferation // used 
	double o2_proliferation_threshold = 5.0; // value at which o2 is sufficient for proliferation // used 
	double o2_reference = 38; // physioxic reference value, in the linked reference Phenotype // used 

	double glucose_proliferation_reference = 1.0; // used
	double glucose_proliferation_saturation = 1.0; // used  
	double glucose_proliferation_threshold = 0.0; // used 
	
	// necrosis 
	double o2_necrosis_threshold = 5.0; // used 
	double o2_necrosis_max = 2.5;  // used 
	double max_necrosis_rate = 0.00278; // 1.0 / (6.0 * 60.0); //used 
	
	// metabolism effects 
	// these two must sum to 1.0 
	// PhysiCell we re-weight these as needed
	double relative_aerobic_effects = 1.0;  // used 
	double relative_glycolytic_effects = 0.0; // used 
	
	// apoptosis 
	double apoptosis_rate = 9.26e-5; // 1/min // used 
	
	// motility 
	bool is_motile = false; // used 
	double motility_bias = 0.5; // used 
	int motility_gradient_substrate_index = 0; // used in custom 
	bool negative_taxis = false; // used in custom 
	double motility_speed = 1.0; // used 
	double motility_persistence_time = 0.0;  // used 

	// mechanics 
	double max_relative_adhesion_distance = 1.5; // used 
	double adhesion_strength = 0.0; // used 
	double repulsion_strength = 0.0; // used 
	
	// hypoxic response (for later)
	
	double o2_hypoxic_threshold = 15.0; // value at which hypoxic signaling starts // used
	double o2_hypoxic_response = 8.0; // value at which omics changes are observed  // usd
	double o2_hypoxic_saturation = 4.0; // value at which hypoxic signalign saturates // used  
	
	// secretion and uptake 
	
	std::vector<double> uptake_rates = {};
	std::vector<double> secretion_rates = {}; 
	std::vector<double> saturation_densities = {}; 
	
	// Basic_Phenotype();
	
	void read_from_pugixml( pugi::xml_node parent_node ); // done 
};

class Pharmacodynamics
{
 private:
 public: 
	int effect_model_type = 0; 
		// 0: simple, 1: intermediate, 2: detailed 
	double EC_50 = 0.5; 
	double Hill_power = 1.0; 
	
	double mechanistic_response_rate = 1.0; 
	double mechanistic_deactivation_rate = 0.0; 
	
/*	
	double enable_influx_efflux = false; 
	double influx = 1.0; 
	double efflux = 0.0; 
	double cell_saturation_concentration = 1.0; 
*/
	
	bool enable_active_influx = true; 
	double relative_max_internal_concentration = 2.0; 
	double internalization_rate = 0.0058;
	double reference_external_concentration = 1.0; 

	
	Basic_Phenotype max_response; // done 
	
	// cycle effects 
	bool cycle_enabled = false; 
	// apoptosis effects 
	bool apoptosis_enabled = false; 
	// metabolism effects 
	bool metabolism_enabled = false; 
	// motility effects 
	bool motility_enabled = false; 
	// mechanics effects
	bool mechanics_enabled = false; 
	// secretion effects
	bool secretion_enabled = false; 
	
	// Pharmacodynamics();
	
	void read_from_pugixml( pugi::xml_node parent_node ); // done 

	// where do I find the corresponding drug in the cell's custom data? 
	
	int cell_Ei = 6.6e-8;  // adding a custom data of Ei source for making pH stable
	
	int cell_NP_i = 0;
	int cell_NP_AUC_i = 0; 
	int cell_response_i = 0;
	int cell_effect_i = 0;
	
	// where do I find the corresponding drug in the BioFVM substrates? 
	int BioFVM_substrate_ID = 0; 

	void advance_internalization( Cell* pCell, Phenotype& phenotype , double dt ); // done
	void advance_mechanistic_response_model( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	
	void (*effect_model)( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt ) = NULL;
	
	void advance_PD_model( Cell* pCell, Phenotype& phenotype, double dt ); 
};

void simple_effect_model( Pharmacodynamics& PD , Cell* pCell, Phenotype& phenotype, double dt ); // done 
void AUC_effect_model( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt ); // done 
void mechanistic_effect_model( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt ); // done 

class NP_transformation
{
 public:
	// NP transition rates
	int NP_condition_substrate_index = 0; 
	int NP_start_index = 0; 
	int NP_end_index = 1; 
	int NP_substrate_start_index = 4;
	int NP_substrate_end_index = 5; 
	
	double condition_1 = 0.0;
	double rate_1 = 0.0; 
	
	double condition_2 = 1.0;
	double rate_2 = 1.0; 	
	
	double get_transition_rate( int voxel_index ); 
	
	void read_from_pugixml(pugi::xml_node parent_node );   
};

class Nanobio_Parameters
{
 private: 
 public:
	// tumor properties 
	double tumor_radius = 350.0; 
	std::string tumor_cell_line_name = "untreated cancer"; 
 
	// baseline phenotype
	Basic_Phenotype untreated_phenotype; // read 
 
	// nanoparticle properties 
	Nanoparticle_Definition NP_preform; // read 
	Nanoparticle_Definition NP_reconfigured; // read 
	
	// pharmacodynamics
	Pharmacodynamics NP_preform_PD;  // read 
	Pharmacodynamics NP_reconfigured_PD; // not read 
	
	// NP transition rates
	NP_transformation NP1_to_NP2; 
	NP_transformation NP2_to_NP1; 
	
	// therapy 
	std::vector<int> biofvm_indices = {4}; 
	std::vector<double> activation_times = {9e9}; 
	std::vector<double> doses = {0}; 
	std::vector<double> clearance_rates = {0};  
	
	
	void read_from_pugixml( void ); // done 
};

extern Nanobio_Parameters nanobio_parameters; 

void parse_nanobio_parameters(void); 





// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_nanoparticle_response( Cell* pCell, Phenotype& phenotype, double dt ); // done 
void tumor_cell_nanoparticle_response( Cell* pCell , Phenotype& phenotype, double dt ); // not done 

// set the tumor cell properties, then call the function 
// to set up the tumor cells 
void create_cell_types( void ); // done, I think  

void setup_tissue(); // done 

// set up the microenvironment to include: 
// oxygen, glucose, H+ ions, and nanoparticle 
void setup_microenvironment( void );  // done! 

std::vector<std::string> nanobio_coloring_function( Cell* );

void NP_ECM_binding_model( double dt ); 

void chemotaxis_function( Cell* pCell , Phenotype& phenotype , double dt ); // done 

void apply_therapies( void ); 

void NP_transformations( double dt );

void add_Ei_source( double dt ); // adding Ei source function !!!!

void update_pH( double dt ); // adding update_pH function !!!

#endif
