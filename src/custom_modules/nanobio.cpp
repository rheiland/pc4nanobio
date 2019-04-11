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

#include "./nanobio.h"

Nanobio_Parameters nanobio_parameters; 

void Nanoparticle_Definition::read_from_pugixml( pugi::xml_node parent_node )
{
	name = parent_node.attribute( "name" ).value();
	units = parent_node.attribute( "units" ).value();
	ID = parent_node.attribute( "ID" ).as_int(); 
	substrate_ID = parent_node.attribute( "substrate_ID" ).as_int(); 
	
	// read pharmacokinetics 
	
	pugi::xml_node node = xml_find_node( parent_node , "pharmacokinetics" ); 
	
	diffusion_coefficient = xml_get_double_value( node , "diffusion_coefficient" );
	survival_lifetime = xml_get_double_value( node , "survival_lifetime" );
	ECM_binding_rate = xml_get_double_value( node , "ECM_binding_rate" );
	ECM_unbinding_rate = xml_get_double_value( node , "ECM_unbinding_rate" );
	ECM_saturation_concentration = xml_get_double_value( node , "ECM_saturation_concentration" );
	
	degradation_rate = 1.0 / (1e-17 + survival_lifetime); 
	
	node = node.parent(); 
	
	return; 
}

// pass the XML node belonging to a basic_phenotype 
void Basic_Phenotype::read_from_pugixml( pugi::xml_node parent_node )
{
	pugi::xml_node node = parent_node; 
	
	// phenotype
	name = node.attribute("name").value(); 
	
	// cycle 
	node = xml_find_node( parent_node, "cycle" ); 
	max_birth_rate = xml_get_double_value( node , "max_birth_rate" ); 
	// cycle parameters 
	node = xml_find_node( node, "parameters" );
	
	o2_proliferation_saturation = xml_get_double_value( node , "o2_proliferation_saturation" );
	o2_proliferation_threshold = xml_get_double_value( node , "o2_proliferation_threshold" );
	o2_reference = xml_get_double_value( node , "o2_reference" );
	
	glucose_proliferation_reference = xml_get_double_value( node , "glucose_proliferation_reference" );
	glucose_proliferation_saturation = xml_get_double_value( node , "glucose_proliferation_saturation" );
	glucose_proliferation_threshold = xml_get_double_value( node, "glucose_proliferation_threshold" ); 
	
	// necrosis 
	node = xml_find_node( parent_node, "necrosis" ); 
	max_necrosis_rate = xml_get_double_value( node , "max_necrosis_rate" ); 
	// necrosis parameters 
	node = xml_find_node( node, "parameters" );
 	o2_necrosis_threshold = xml_get_double_value( node , "o2_necrosis_threshold" ); 
	o2_necrosis_max = xml_get_double_value( node, "o2_necrosis_max" );
	
	// apoptosis 
	node = xml_find_node( parent_node, "apoptosis" ); 
	apoptosis_rate = xml_get_double_value( node , "apoptosis_rate" ); 	
	
	// metabolism 
	// metabolism parameters 
	node = xml_find_node( parent_node, "metabolism" ); 
	node = xml_find_node( node, "parameters" );
	relative_aerobic_effects = xml_get_double_value( node , "relative_aerobic_effects" ); 	
	relative_glycolytic_effects = xml_get_double_value( node , "relative_glycolytic_effects" ); 	
	double sum = relative_aerobic_effects + relative_glycolytic_effects; 
	relative_aerobic_effects /= sum; 
	relative_glycolytic_effects /= sum; 
	
	// motility  
	node = xml_find_node( parent_node, "motility" ); 
	is_motile = xml_get_bool_value( node, "is_motile" ); 
	motility_bias = xml_get_double_value( node , "bias" ); 
	motility_gradient_substrate_index = xml_get_int_value( node , "motility_gradient_substrate_index" ); 
	negative_taxis = xml_get_bool_value( node , "negative_taxis" ); 
	motility_speed = xml_get_double_value( node , "speed" ); 
	motility_persistence_time = xml_get_double_value( node , "persistence_time" ); 
	
	// mechanics 
	node = xml_find_node( parent_node, "mechanics" ); 
	max_relative_adhesion_distance = xml_get_double_value( node , "max_relative_adhesion_distance" ); 
	adhesion_strength = xml_get_double_value( node , "adhesion_strength" ); 
	repulsion_strength = xml_get_double_value( node , "repulsion_strength" ); 
	
	// hypoxic response (for future)
	node = xml_find_node( parent_node, "hypoxia" ); 
	o2_hypoxic_threshold = xml_get_double_value( node , "o2_hypoxic_threshold" ); 
	o2_hypoxic_response = xml_get_double_value( node , "o2_hypoxic_response" ); 
	o2_hypoxic_saturation = xml_get_double_value( node , "o2_hypoxic_saturation" ); 
	
	// secretion and uptake 
	node = xml_find_node( parent_node, "secretion" ); 
	node = xml_find_node( node , "substrate" ); 
	while( node )
	{
		uptake_rates.push_back( xml_get_double_value( node , "uptake_rate" ) ); 
		secretion_rates.push_back( xml_get_double_value( node , "secretion_rate" ) ); 
		saturation_densities.push_back( xml_get_double_value( node , "saturation_density" ) ); 
		
		node = node.next_sibling( "substrate" ); 
	}
	
	return; 
}

void Nanobio_Parameters::read_from_pugixml(void)
{
	pugi::xml_node node; 
	pugi::xml_node node1;
	
	node = xml_find_node( physicell_config_root , "tumor" );
	
	// tumor 
	tumor_radius = xml_get_double_value( node , "radius" );
	node = node.parent(); 

	// get cancer cell line properties 
	
	node = xml_find_node( physicell_config_root , "cell_definitions" ); 
	node = xml_find_node( node , "cell_definition" ); // first (and only definition is untreated tumor) 
	
	tumor_cell_line_name = node.attribute("name" ).value(); 
	
	node = xml_find_node( node , "basic_phenotype" ); 
	untreated_phenotype.read_from_pugixml( node );
	
	// now, get nanoparticle properties, including pharmacodynamics
		
	node = xml_find_node( physicell_config_root , "nanoparticle_definitions" ); 
	node = xml_find_node( node, "nanoparticle" ); 
	NP_preform.read_from_pugixml( node ); 
	node = xml_find_node( node, "pharmacodynamics" ); 
	NP_preform_PD.read_from_pugixml( node ); 
	node = node.parent(); 
	
	// add proper custom variables to cell default for this NP1
	NP_preform_PD.cell_NP_i = cell_defaults.custom_data.add_variable( "NP1" , "dimensionless" , 0.0 );
	NP_preform_PD.cell_NP_AUC_i = cell_defaults.custom_data.add_variable( "NP1_AUC" , "dimensionless" , 0.0 );
	NP_preform_PD.cell_response_i = cell_defaults.custom_data.add_variable( "NP1 response" , "dimensionless" , 0.0 );
	NP_preform_PD.cell_effect_i = cell_defaults.custom_data.add_variable( "NP1 effect" , "dimensionless" , 0.0 );
	
	NP_preform_PD.BioFVM_substrate_ID =  NP_preform.substrate_ID;

	// move on to NP2
	
	node = node.next_sibling( "nanoparticle" ); 
	NP_reconfigured.read_from_pugixml( node ); 
	node = xml_find_node( node, "pharmacodynamics" ); 
	NP_reconfigured_PD.read_from_pugixml( node ); 
	node = node.parent(); 

	// add proper custom variables to cell default for this NP2
	NP_reconfigured_PD.cell_NP_i = cell_defaults.custom_data.add_variable( "NP2" , "dimensionless" , 0.0 );
	NP_reconfigured_PD.cell_NP_AUC_i = cell_defaults.custom_data.add_variable( "NP2_AUC" , "dimensionless" , 0.0 );
	NP_reconfigured_PD.cell_response_i = cell_defaults.custom_data.add_variable( "NP2 response" , "dimensionless" , 0.0 );
	NP_reconfigured_PD.cell_effect_i = cell_defaults.custom_data.add_variable( "NP2 effect" , "dimensionless" , 0.0 );
	
	NP_reconfigured_PD.BioFVM_substrate_ID =  NP_reconfigured.substrate_ID;
	
	// now, get transformation(s)
	 
	 
	node = xml_find_node( physicell_config_root , "nanoparticle_definitions" ); 
	node = xml_find_node( node , "transformations" );
	node = xml_find_node( node , "transformation" ); 
	NP1_to_NP2.read_from_pugixml( node ); 
	
	node = node.next_sibling( "transformation" ); 
	NP2_to_NP1.read_from_pugixml( node ); 
	
	// now, set up therapy 
	
	node = xml_find_node( physicell_config_root , "treatment" ); 
	node = xml_find_node( node, "substrates" ); 
	node1 = xml_find_node( node, "substrate" ); 
	
	int k=0; 
	
	biofvm_indices.clear();
	activation_times.clear();
	doses.clear(); 
	clearance_rates.clear(); 
	
	std::cout << std::endl << "Parsing treatment information ... " << std::endl; 
	
	while( node1 )
	{
		biofvm_indices.push_back( node1.attribute("ID").as_int() ); 
		activation_times.push_back( xml_get_double_value( node1, "activation_time" ) ); 
		doses.push_back( xml_get_double_value( node1, "dose" ) ); 
		clearance_rates.push_back( xml_get_double_value( node1, "clearance_rate" ) ); 
		
		std::cout << "\tApplying treatment to substrate " << biofvm_indices[k] << " with dose = " << doses[k] << " at time " << activation_times[k] << " cleared at rate " << clearance_rates[k] << std::endl; 
		
		k++; 
		node1 = node1.next_sibling( "substrate" ); 
	}
	
	return; 
}

void parse_nanobio_parameters(void)
{
	nanobio_parameters.read_from_pugixml();
	
	return; 
}


void create_cell_types( void ) // done 
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// set functions 
	cell_defaults.functions.update_migration_bias = chemotaxis_function; 
	
	// use default proliferation and death 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation = nanobio_parameters.untreated_phenotype.o2_proliferation_saturation;  
	cell_defaults.parameters.o2_proliferation_threshold = nanobio_parameters.untreated_phenotype.o2_proliferation_threshold;  
	cell_defaults.parameters.o2_reference = nanobio_parameters.untreated_phenotype.o2_reference; 
	
	
	cell_defaults.custom_data.add_variable( "glucose_proliferation_reference" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.glucose_proliferation_reference ); 
	cell_defaults.custom_data.add_variable( "glucose_proliferation_saturation" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.glucose_proliferation_saturation ); 
	cell_defaults.custom_data.add_variable( "glucose_proliferation_threshold" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.glucose_proliferation_threshold ); 

	// set default uptake and secretion 
	
	cell_defaults.phenotype.secretion.secretion_rates = nanobio_parameters.untreated_phenotype.secretion_rates; 
	cell_defaults.phenotype.secretion.uptake_rates = nanobio_parameters.untreated_phenotype.uptake_rates; 
	cell_defaults.phenotype.secretion.saturation_densities = nanobio_parameters.untreated_phenotype.saturation_densities; 
	
	// modify the uptake rates for the NPs 
	cell_defaults.phenotype.secretion.uptake_rates[ nanobio_parameters.NP_preform_PD.BioFVM_substrate_ID] 
		= nanobio_parameters.NP_preform_PD.internalization_rate; // r
	cell_defaults.phenotype.secretion.uptake_rates[ nanobio_parameters.NP_preform_PD.BioFVM_substrate_ID] 
		*= nanobio_parameters.NP_preform_PD.relative_max_internal_concentration; // r*R 
		
	cell_defaults.phenotype.secretion.uptake_rates[ nanobio_parameters.NP_reconfigured_PD.BioFVM_substrate_ID] 
		= nanobio_parameters.NP_reconfigured_PD.internalization_rate; // r
	cell_defaults.phenotype.secretion.uptake_rates[ nanobio_parameters.NP_reconfigured_PD.BioFVM_substrate_ID] 
		*= nanobio_parameters.NP_reconfigured_PD.relative_max_internal_concentration; // r*R 
	
	// // set reference phenotype properties 
	// cycling
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) = nanobio_parameters.untreated_phenotype.max_birth_rate; 

	// mechanics 
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = nanobio_parameters.untreated_phenotype.adhesion_strength; 
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = nanobio_parameters.untreated_phenotype.repulsion_strength; 
	cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance = nanobio_parameters.untreated_phenotype.max_relative_adhesion_distance; 
	
	// apoptosis 
	cell_defaults.phenotype.death.rates[apoptosis_index] = nanobio_parameters.untreated_phenotype.apoptosis_rate; 
	
	// necrosis 
	cell_defaults.parameters.max_necrosis_rate = nanobio_parameters.untreated_phenotype.max_necrosis_rate; 	
	cell_defaults.parameters.necrosis_type = PhysiCell_constants::stochastic_necrosis;;
	cell_defaults.parameters.o2_necrosis_threshold = nanobio_parameters.untreated_phenotype.o2_necrosis_threshold;  
	cell_defaults.parameters.o2_necrosis_max = nanobio_parameters.untreated_phenotype.o2_necrosis_max; 
	
	// metabolism
	cell_defaults.custom_data.add_variable( "relative_aerobic_effects" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.relative_aerobic_effects ); 
	cell_defaults.custom_data.add_variable( "relative_glycolytic_effects" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.relative_glycolytic_effects ); 
	
	// motility 
	cell_defaults.phenotype.motility.migration_bias = nanobio_parameters.untreated_phenotype.motility_bias; 
	cell_defaults.phenotype.motility.migration_speed = nanobio_parameters.untreated_phenotype.motility_speed; 
	cell_defaults.phenotype.motility.persistence_time = nanobio_parameters.untreated_phenotype.motility_persistence_time; 
	cell_defaults.phenotype.motility.is_motile = nanobio_parameters.untreated_phenotype.is_motile; 

	cell_defaults.custom_data.add_variable( "motility_gradient_substrate_index" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.motility_gradient_substrate_index  ); 
	cell_defaults.custom_data.add_variable( "negative_taxis" , "dimensionless", 
		nanobio_parameters.untreated_phenotype.negative_taxis ); 
	
	// hypoxic response 
	cell_defaults.parameters.o2_hypoxic_threshold = nanobio_parameters.untreated_phenotype.o2_hypoxic_threshold;  
	cell_defaults.parameters.o2_hypoxic_response = nanobio_parameters.untreated_phenotype.o2_hypoxic_response;  
	cell_defaults.parameters.o2_hypoxic_saturation = nanobio_parameters.untreated_phenotype.o2_hypoxic_saturation; 
	
	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_nanoparticle_response; 
	
	cell_defaults.name = nanobio_parameters.tumor_cell_line_name; 
	cell_defaults.type = 0; 
	
	// now, add custom variables for the PD 
	// this is done in the NP setup, so they "know" where their damage is stored 
	
	return; 
}


void setup_microenvironment( void )
{
	pugi::xml_node node; 

	// set domain parameters

	// extract from XML 
	
	node = xml_find_node( physicell_config_root , "domain" );

	double xmin = xml_get_double_value( node , "x_min" );
	double xmax = xml_get_double_value( node , "x_max" );
	double ymin = xml_get_double_value( node , "y_min" );
	double ymax = xml_get_double_value( node , "y_max" );
	double zmin = xml_get_double_value( node , "z_min" );
	double zmax = xml_get_double_value( node , "z_max" );
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 
	
	node = node.parent(); 
	
	default_microenvironment_options.X_range = {xmin, xmax}; 
	default_microenvironment_options.Y_range = {ymin, ymax}; 
	default_microenvironment_options.simulate_2D = true; 
	
	default_microenvironment_options.dx = dx; 
	default_microenvironment_options.dy = dy; 
	default_microenvironment_options.dz = dz; 
	
	// configure the non-nanoparticle substrates 
	
	// override oxygen as first. we'll set it ourselves. 
	// Don't let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = false; 

	
	node = xml_find_node( physicell_config_root , "substrates" ); 
	
	pugi::xml_node node1 = node.child( "variable" ); // xml_find_node( node , "variable" ); 
	node = node1; 
	int i = 0; 
	while( node )
	{
		// name and units 
		std::string name = node.attribute( "name" ).value(); 
		std::string units = node.attribute( "units" ).value(); 
		
		if( i == 0 )
		{ microenvironment.set_density( 0, name, units ); }
		else
		{ microenvironment.add_density( name, units ); }
		
		node1 = xml_find_node( node, "physical_parameter_set" ); 
		
		microenvironment.diffusion_coefficients[i] = 
			xml_get_double_value( node1, "diffusion_coefficient" ); 
		microenvironment.decay_rates[i] = 
			xml_get_double_value( node1, "decay_rate" ); 
			
		// if this is one of our nanoparticles, use those properties.
		if( i == nanobio_parameters.NP_preform.substrate_ID )
		{
			microenvironment.set_density( i, nanobio_parameters.NP_preform.name, nanobio_parameters.NP_preform.units ); 
			microenvironment.diffusion_coefficients[i] = nanobio_parameters.NP_preform.diffusion_coefficient; 
			microenvironment.decay_rates[i] = nanobio_parameters.NP_preform.degradation_rate;
		}
		if( i == nanobio_parameters.NP_reconfigured.substrate_ID )
		{
			microenvironment.set_density( i, nanobio_parameters.NP_reconfigured.name, nanobio_parameters.NP_reconfigured.units ); 
			microenvironment.diffusion_coefficients[i] = nanobio_parameters.NP_reconfigured.diffusion_coefficient; 
			microenvironment.decay_rates[i] = nanobio_parameters.NP_reconfigured.degradation_rate;
		}
		
//		node1 = node1.parent(); 
		
		node = node.next_sibling( "variable" ); 
		i++; 
	}
	
	// calculate gradients? 
	
	default_microenvironment_options.calculate_gradients = xml_get_bool_value( node, "calculate_gradients" );  

	// now add nanoparticles 
	
	/* the following code is no longer executed */

	node = xml_find_node( physicell_config_root , "nanoparticles" ); 
	node1 = xml_find_node( node , "nanoparticle" ); 
	while( node1 )
	{
		// name and units 
		std::string name = node1.attribute( "name" ).value(); 
		std::string units = node1.attribute( "units" ).value(); 

		// am I already in BioFVM? 
		int biofvm_index = node1.attribute( "substrate_ID" ).as_int(); 
		if( microenvironment.density_names[biofvm_index] == "NP1" )
		{ std::cout << "hiya!" << std::endl;}
		

		microenvironment.add_density( name, units );
		
		
		node1 = xml_find_node( node1 ,  "pharmacokinetics" ); 
		microenvironment.diffusion_coefficients[i] = 
			xml_get_double_value( node1, "diffusion_coefficient" ); 
		microenvironment.decay_rates[i] = 1.0 / 
			xml_get_double_value( node1, "survival_lifetime" ); 
			
		node1 = node1.parent(); 
		
		node1 = node1.next_sibling( "nanoparticle" ); 
		i++; 
	}
	
	// set Dirichlet conditions 

	i=0; 
	default_microenvironment_options.outer_Dirichlet_conditions = true;

	// substrates first 
	node = xml_find_node( physicell_config_root , "tissue_conditions" ); 
	node = xml_find_node( node , "substrates" ); 
	node1 = xml_find_node( node , "substrate" ); 
	while( node1 )
	{
		default_microenvironment_options.Dirichlet_condition_vector[i] = 
			node1.text().as_double(); 
			
		default_microenvironment_options.Dirichlet_activation_vector[i] = 
			node1.attribute( "enabled" ).as_bool(); 
		
		node1 = node1.next_sibling( "substrate" ); 
		i++; 
	}
	node = node.parent(); 
	
	// now nanoparticles 
	node = xml_find_node( node , "nanoparticles" ); 
	node1 = xml_find_node( node , "nanoparticle" ); 
	while( node1 )
	{
		default_microenvironment_options.Dirichlet_condition_vector[i] = 
			node1.text().as_double(); 
			
		default_microenvironment_options.Dirichlet_activation_vector[i] = 
			node1.attribute( "enabled" ).as_bool(); 
		
		node1 = node1.next_sibling( "nanoparticle" ); 
		i++; 
	}
	
	initialize_microenvironment(); 	

	return; 
}	

void setup_tissue( void )
{
	pugi::xml_node node = xml_find_node( physicell_config_root , "domain" );
	
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = nanobio_parameters.tumor_radius; 
	
	Cell* pCell = NULL; 
	
	double center_x = 0.5*(  xml_get_double_value( node , "x_min" ) + xml_get_double_value( node , "x_max" ) );
	double center_y = 0.5*(  xml_get_double_value( node , "y_min" ) + xml_get_double_value( node , "y_max" ) );
		
	double x = 0; 
	double x_outer =  tumor_radius; 
	double y = 0; 
	
	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(); // tumor cell 
			pCell->assign_position( center_x+x , center_y+y , 0.0 );
			/*
			pCell->custom_data[0] = NormalRandom( 1.0, 0.33 );
			if( pCell->custom_data[0] < 0.0 )
			{ pCell->custom_data[0] = 0.0; }
			if( pCell->custom_data[0] > 2.0 )
			{ pCell->custom_data[0] = .0; }
			*/
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( center_x+x , center_y-y , 0.0 );
				/*
				pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
				if( pCell->custom_data[0] < 0.0 )
				{ pCell->custom_data[0] = 0.0; }
				if( pCell->custom_data[0] > 2.0 )
				{ pCell->custom_data[0] = .0; }
				*/
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( center_x-x , center_y+y , 0.0 );
				/*
				pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
				if( pCell->custom_data[0] < 0.0 )
				{ pCell->custom_data[0] = 0.0; }
				if( pCell->custom_data[0] > 2.0 )
				{ pCell->custom_data[0] = .0; }
				*/
				
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( center_x-x , center_y-y , 0.0 );
					/*
					pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
					if( pCell->custom_data[0] < 0.0 )
					{ pCell->custom_data[0] = 0.0; }
					if( pCell->custom_data[0] > 2.0 )
					{ pCell->custom_data[0] = .0; }
					*/
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	/*
	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[0]; 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		sum +=  ( (*all_cells)[i]->custom_data[0] - mean )*( (*all_cells)[i]->custom_data[0] - mean ); 
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
	
	std::cout << std::endl << "Oncoprotein summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 
	*/
	
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_nanoparticle_response( Cell* pCell, Phenotype& phenotype, double dt )
{
	// supported cycle models:
		// advanced_Ki67_cycle_model= 0;
		// basic_Ki67_cycle_model=1
		// live_cells_cycle_model = 5; 
	
	if( phenotype.death.dead == true )
	{ return; }
	
	// set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
	static bool indices_initiated = false; 
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	static int necrosis_index; 
	static int apoptosis_index; 
	
	static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	static int glucose_i = pCell->get_microenvironment()->find_density_index( "glucose" ); 
	
	static int glucose_prolif_ref_i = pCell->custom_data.find_variable_index( "glucose_proliferation_reference" );
	static int glucose_prolif_sat_i = pCell->custom_data.find_variable_index( "glucose_proliferation_saturation" );
	static int glucose_prolif_threshold_i = pCell->custom_data.find_variable_index( "glucose_proliferation_threshold" );
	
	static int r1_i = pCell->custom_data.find_variable_index( "relative_aerobic_effects" );
	static int r2_i = pCell->custom_data.find_variable_index( "relative_glycolytic_effects" );
		
	if( indices_initiated == false )
	{
		// Ki67 models
		
		if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_negative );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );  
			
			if( phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive );
				indices_initiated = true; 
			}
			if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
				indices_initiated = true; 
			}
		}
		
		// live model 
			
		if( phenotype.cycle.model().code ==  PhysiCell_constants::live_cells_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );  
			
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			indices_initiated = true; 
		}
		
		// cytometry models 
		
		if( phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_separated_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::G0G1_phase );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );  

			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::S_phase );
			indices_initiated = true; 
		}		
		
	}
	
	// don't continue if we never "figured out" the current cycle model. 
	if( indices_initiated == false )
	{
		return; 
	}

	// sample the microenvironment to get the pO2 value 
	
	double pO2 = (pCell->nearest_density_vector())[oxygen_i]; // PhysiCell_constants::oxygen_index]; 
	double g = (pCell->nearest_density_vector())[glucose_i]; // PhysiCell_constants::oxygen_index]; 
	
	int n = pCell->phenotype.cycle.current_phase_index(); 
	
	// function:
	//  birth_rate_max*(  r1*(pO2-pO2_min)/(pO2_ref-pO2_min)*g + r2*(g-g_min)/(g_ref-g_min) ); 
	// mult_o2 = max( 0 , min( (pO2-pO2_min)/(pO2_ref-pO2_min) , 1 ))
	// mult_g = max( 0, mind( (g-g_min)/(g_ref-g_min) ,1 ))
	// mult = ( r1*mult_o2*mult_g + r2*mult_g ); 
	
	// this multiplier is for linear interpolation of the oxygen value 
	double O2 = 1.0;
	if( pO2 < pCell->parameters.o2_proliferation_saturation )
	{
		O2 = ( pO2 - pCell->parameters.o2_proliferation_threshold ) 
			/ ( pCell->parameters.o2_proliferation_saturation - pCell->parameters.o2_proliferation_threshold );
	}
	if( pO2 < pCell->parameters.o2_proliferation_threshold )
	{ 
		O2 = 0.0; 
	}
	double GLUCOSE = 1.0;
	if( g < pCell->custom_data[glucose_prolif_sat_i] )
	{
		GLUCOSE = ( g - pCell->custom_data[glucose_prolif_threshold_i] ) 
			/ ( pCell->custom_data[glucose_prolif_sat_i] - pCell->custom_data[glucose_prolif_threshold_i]);
	}
	if( g < pCell->custom_data[glucose_prolif_threshold_i] )
	{ 
		GLUCOSE = 0.0; 
	}
	
	double r1 = pCell->custom_data[r1_i] ; 
	double r2 = pCell->custom_data[r2_i] ; 
	
	double multiplier = O2; // O2 
	multiplier *= r1; // r1*O2
	multiplier += r2; // r1*O2 + r2 
	multiplier *= GLUCOSE; // r1*O2*GLUCOSE + r2*GLUCOSE 
	
	// now, update the appropriate cycle transition rate 
	
	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier * 
		pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
	
	// Update necrosis rate
	multiplier = 0.0;
	if( pO2 < pCell->parameters.o2_necrosis_threshold )
	{
		multiplier = ( pCell->parameters.o2_necrosis_threshold - pO2 ) 
			/ ( pCell->parameters.o2_necrosis_threshold - pCell->parameters.o2_necrosis_max );
	}
	if( pO2 < pCell->parameters.o2_necrosis_max )
	{ 
		multiplier = 1.0; 
	}	
	
	// now, update the necrosis rate 
	
	pCell->phenotype.death.rates[necrosis_index] = multiplier * pCell->parameters.max_necrosis_rate; 
	
	// therapy response here! 
	

	// influx / efflux, effect
	
	nanobio_parameters.NP_preform_PD.advance_internalization( pCell, phenotype , dt ); 
	nanobio_parameters.NP_preform_PD.effect_model( (nanobio_parameters.NP_preform_PD), pCell, phenotype , dt ); 
	nanobio_parameters.NP_reconfigured_PD.advance_internalization( pCell, phenotype , dt ); 
	nanobio_parameters.NP_reconfigured_PD.effect_model( (nanobio_parameters.NP_reconfigured_PD), pCell, phenotype , dt ); 
 	
	if( pCell->custom_data[ nanobio_parameters.NP_preform_PD.cell_effect_i ] > 
		pCell->custom_data[ nanobio_parameters.NP_reconfigured_PD.cell_effect_i ] )
	{
		nanobio_parameters.NP_preform_PD.advance_PD_model( pCell, phenotype, dt );
	}
	else
	{
		nanobio_parameters.NP_reconfigured_PD.advance_PD_model( pCell, phenotype, dt ); 
	}

	return; 
}

std::vector<std::string> nanobio_coloring_function( Cell* pCell )
{
	static bool color_initialized = false; 
	static int custom_data_index = 0;
	static double blue_value = 0; 
	static double yellow_value = 1; 
	static double diff = 1.0; 
	
	if( color_initialized == false )
	{
		pugi::xml_node node = xml_find_node( physicell_config_root , "visualization" ); 
		custom_data_index = xml_get_int_value( node , "custom_data_index" ); 
		blue_value = xml_get_double_value( node , "blue_value" ); 
		yellow_value = xml_get_double_value( node , "yellow_value" ); 
		
		diff = yellow_value - blue_value; 
		
		
		std::cout << "Shading cells according to custom varaiable " 
			<< custom_data_index << " : " 
			<< pCell->custom_data.variables[ custom_data_index ].name << std::endl; 
		std::cout << "\tIf this is not what you expect, change: " << std::endl << std::endl 
			<< "\t<visualization>\n\t\t<custom_data_index />\n\t</visualization>" 
			<< std::endl << std::endl;  
		
		color_initialized = true; 
	}
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ return output; } 
	
	// live cells are green, but shaded by internal NP value 
	if( pCell->phenotype.death.dead == false )
	{
		// nanobio_parameters.NP_preform_PD.cell_NP_i; 
		
		double scaled_value = pCell->custom_data[custom_data_index];
		
		scaled_value -= blue_value; 
		scaled_value /= diff;

		if( scaled_value > 1.0 )
		{ scaled_value = 1.0; } 
		if( scaled_value < 0.0 )
		{ scaled_value = 0.0; } 
		
		int SV = (int) round( scaled_value * 255.0 ); 
		
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", SV, SV, 255-SV );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

void NP_ECM_binding_model( double dt )
{
	
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels(); i++ )
	{
		// calculate binding 
		
		// calculate unbinding 
		
		// get it done! 
		
		
	}
	
	return; 
}


void chemotaxis_function( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int motility_gradient_substrate_index_i = pCell->custom_data.find_variable_index( "motility_gradient_substrate_index" ); 
	static int negative_taxis_i = pCell->custom_data.find_variable_index( "negative_taxis" ); 

	if( phenotype.motility.is_motile == false || phenotype.motility.migration_bias < 1e-5 )
	{ return; }		

	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( (int) pCell->custom_data[motility_gradient_substrate_index_i] );
	normalize( &( phenotype.motility.migration_bias_direction ) );	
	if( (bool) pCell->custom_data[negative_taxis_i] )
	{ phenotype.motility.migration_bias_direction *= -1.0; }
	
	return; 
}


void Pharmacodynamics::read_from_pugixml( pugi::xml_node parent_node )
{
	// basic form of effects model, and import/export parameters 
	pugi::xml_node node = parent_node; 
	effect_model_type = xml_get_int_value( node, "effect_model" ); 
	EC_50 = xml_get_double_value( node, "EC_50" ); 
	Hill_power = xml_get_double_value( node, "Hill_power" ); 
	mechanistic_response_rate = xml_get_double_value( node, "mechanistic_response_rate" ); 
	mechanistic_deactivation_rate = xml_get_double_value( node, "mechanistic_deactivation_rate" ); 
		
	if( effect_model_type == 0 )
	{ effect_model = simple_effect_model; }
	if( effect_model_type == 1 )
	{ effect_model = AUC_effect_model; }
	if( effect_model_type == 2 )
	{ effect_model = mechanistic_effect_model; }

	enable_active_influx = xml_get_bool_value( node, "enable_active_influx" ); 
	relative_max_internal_concentration = xml_get_double_value( node , "relative_max_internal_concentration" ); 
	reference_external_concentration = xml_get_double_value( node, "reference_external_concentration"); 
	internalization_rate = xml_get_double_value( node , "internalization_rate" );
	
	// determine which parts of the PD model are enabled
	node = xml_find_node( parent_node , "enabled_responses" );
	cycle_enabled = xml_get_bool_value( node , "cycle" ); 
	apoptosis_enabled = xml_get_bool_value( node , "apoptosis" ); 
	metabolism_enabled = xml_get_bool_value( node , "metabolism" ); 
	motility_enabled = xml_get_bool_value( node , "motility" ); 
	mechanics_enabled = xml_get_bool_value( node , "mechanics" ); 
	secretion_enabled = xml_get_bool_value( node , "secretions" ); 
	node = node.parent(); 
	
	// now read the max response phenotype 
	node = xml_find_node( parent_node , "max_response" ); 
	node = xml_find_node( node , "basic_phenotype" ); 
	max_response.read_from_pugixml( node ); 
	node = node.parent();
	node = node.parent(); 
	
	return; 
}

void Pharmacodynamics::advance_internalization( Cell* pCell, Phenotype& phenotype , double dt )
{
	if( enable_active_influx )
	{ 
		// backwards Euler
		double temp = internalization_rate; // r 
		temp *= (pCell->nearest_density_vector())[BioFVM_substrate_ID]; // r*rho 
		temp *= dt; // dt*r*rho 
		
		pCell->custom_data[ cell_NP_i ] += temp; // n(i) = n(i-1) + dt*r*rho 
		temp += 1.0; // 1 + dt*r*rho 
		pCell->custom_data[ cell_NP_i ] /= temp; // n(i) = ( n(i-1) + dt*r*rho )/( 1 + dt*r*rho );
			// ==> (  n(i) - n(i-1) )/dt = r*( 1 - n(i) )*rho(i-1) 
			
		// if rho > 1 by numerical error, then rho grows without bound for backwards euler! 
//		if( pCell->custom_data[ cell_NP_i ] > 1.0 )
//		{ pCell->custom_data[ cell_NP_i ] = 1.0; }

		// set the tissue-level uptake parameter 
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] = 1.0; // 1.0 
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] -= pCell->custom_data[ cell_NP_i ]; // 1 - rho_i 
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] *= internalization_rate; // r*(1-rho_i)
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] *= relative_max_internal_concentration; // r*R*(1-rho_i) 
	}
	else
	{
		// set the internal drug level equal to the extracellular value 
		 pCell->custom_data[ cell_NP_i ] = 	(pCell->nearest_density_vector())[BioFVM_substrate_ID]; 
		 
		// set the tissue-level uptake parameter 
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] = internalization_rate; // r
		phenotype.secretion.uptake_rates[BioFVM_substrate_ID] *= relative_max_internal_concentration; // r*R
	}
	
	// update AUC 
	
	pCell->custom_data[ cell_NP_AUC_i ] += dt* pCell->custom_data[ cell_NP_i ];
	
	return; 
}


//void Pharmacodynamics::simple_effect_model( Cell* pCell, Phenotype& phenotype, double dt )
void simple_effect_model( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt )
{
	pCell->custom_data[ PD.cell_effect_i ] = Hill_function( pCell->custom_data[ PD.cell_NP_i ] , PD.Hill_power , PD.EC_50 ); 
	
	return; 
}

// void Pharmacodynamics::AUC_effect_model( Cell* pCell, Phenotype& phenotype, double dt )
void AUC_effect_model( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt )
{
	pCell->custom_data[ PD.cell_effect_i ] = Hill_function( pCell->custom_data[ PD.cell_NP_AUC_i ] , PD.Hill_power , PD.EC_50 ); 
	
	return; 
}

// void Pharmacodynamics::mechanistic_effect_model( Cell* pCell, Phenotype& phenotype, double dt )
void mechanistic_effect_model( Pharmacodynamics& PD, Cell* pCell, Phenotype& phenotype, double dt )
{
	pCell->custom_data[ PD.cell_effect_i ] = Hill_function( pCell->custom_data[ PD.cell_NP_AUC_i ] , PD.Hill_power , PD.EC_50 ); 
	
	return; 
}

void Pharmacodynamics::advance_mechanistic_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// dResponse/dt = response*rate*NP - deactivation_rate*Response; 
	
	double temp = pCell->custom_data[ cell_NP_i ] * mechanistic_response_rate;  //  response_rate*NP
	temp -= pCell->custom_data[ cell_response_i ] * mechanistic_deactivation_rate; // response_rate*NP - deactivation_rate*Response
	temp *= dt; // dt*( response_rate*NP - deactivation_rate*Response )
	
	pCell->custom_data[ cell_response_i ] += temp; 
	
	pCell->custom_data[ cell_effect_i ] = Hill_function( pCell->custom_data[ cell_response_i ] , Hill_power , EC_50 ); 
	
	return; 
}


double Hill_function( double input, double Hill_power , double EC_50 )
{
	double temp = input; // c
	temp /= EC_50; // c/c0
	temp = std::pow( temp, Hill_power ); // (c/c0)^n 
	double output = temp; // (c/c0)^n
	temp += 1.0; // 1 + (c/c0)^n 
	output /= temp; // // (c/c0)^n / ( 1 + (c/c0)^n ) 
	
	return output; 
}

void Pharmacodynamics::advance_PD_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// influx / efflux 
	
//	advance_internalization( pCell, phenotype , dt ); 
	
	// effect 
	
//	effect_model( *this, pCell, phenotype , dt ); 
	
	// alterations to phenotype 
	
	static int start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );	
	static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	
	// cell cycle 
	
	if( cycle_enabled )
	{
		interpolate( 
			&(phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index)) , // overwrite this
			phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) , // unperturbed
			max_response.max_birth_rate ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
	}
	
	// apoptosis 
	
	if( apoptosis_enabled )
	{
		interpolate( 
			&(phenotype.death.rates[apoptosis_index]) , // overwrite this
			nanobio_parameters.untreated_phenotype.apoptosis_rate , // unperturbed
			max_response.apoptosis_rate ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
	}
	
	// motility 

	if( motility_enabled )
	{
		interpolate( 
			&(phenotype.motility.migration_bias) , // overwrite this
			nanobio_parameters.untreated_phenotype.motility_bias , // unperturbed
			max_response.motility_bias ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.motility.migration_speed) , // overwrite this
			nanobio_parameters.untreated_phenotype.motility_speed , // unperturbed
			max_response.motility_speed ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.motility.persistence_time) , // overwrite this
			nanobio_parameters.untreated_phenotype.motility_persistence_time , // unperturbed
			max_response.motility_persistence_time ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 

		if( pCell->custom_data[ cell_effect_i ] >= 0.5 )
		{
			static int motility_gradient_substrate_index_i = pCell->custom_data.find_variable_index( "motility_gradient_substrate_index" ); 
			static int negative_taxis_i = pCell->custom_data.find_variable_index( "negative_taxis" ); 
			
			phenotype.motility.is_motile = max_response.is_motile; 
			pCell->custom_data[motility_gradient_substrate_index_i] = max_response.motility_gradient_substrate_index; 
			pCell->custom_data[negative_taxis_i] = max_response.negative_taxis; 
			
		}
		else
		{
			static int motility_gradient_substrate_index_i = pCell->custom_data.find_variable_index( "motility_gradient_substrate_index" ); 
			static int negative_taxis_i = pCell->custom_data.find_variable_index( "negative_taxis" ); 

			phenotype.motility.is_motile = nanobio_parameters.untreated_phenotype.is_motile; 
			pCell->custom_data[motility_gradient_substrate_index_i] = nanobio_parameters.untreated_phenotype.motility_gradient_substrate_index; 
			pCell->custom_data[negative_taxis_i] = nanobio_parameters.untreated_phenotype.negative_taxis; 
		}
	}
	
	if( mechanics_enabled )
	{
		interpolate( 
			&(phenotype.mechanics.relative_maximum_adhesion_distance) , // overwrite this
			nanobio_parameters.untreated_phenotype.max_relative_adhesion_distance , // unperturbed
			max_response.max_relative_adhesion_distance ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.mechanics.cell_cell_adhesion_strength) , // overwrite this
			nanobio_parameters.untreated_phenotype.adhesion_strength , // unperturbed
			max_response.adhesion_strength ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.mechanics.cell_cell_repulsion_strength) , // overwrite this
			nanobio_parameters.untreated_phenotype.repulsion_strength , // unperturbed
			max_response.repulsion_strength ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
	}
	
	if( secretion_enabled )
	{
		interpolate( 
			&(phenotype.secretion.uptake_rates) , // overwrite this
			nanobio_parameters.untreated_phenotype.uptake_rates , // unperturbed
			max_response.uptake_rates ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.secretion.secretion_rates) , // overwrite this
			nanobio_parameters.untreated_phenotype.secretion_rates , // unperturbed
			max_response.secretion_rates ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
		
		interpolate( 
			&(phenotype.secretion.saturation_densities) , // overwrite this
			nanobio_parameters.untreated_phenotype.saturation_densities , // unperturbed
			max_response.saturation_densities ,  // max perturb 
			pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
	}
	
	if( metabolism_enabled )
	{
		static int r1_i = pCell->custom_data.find_variable_index( "relative_aerobic_effects" );
		static int r2_i = pCell->custom_data.find_variable_index( "relative_glycolytic_effects" );

		interpolate( 
		&(pCell->custom_data[r1_i]) , // overwrite this
		nanobio_parameters.untreated_phenotype.relative_aerobic_effects , // unperturbed
		max_response.relative_aerobic_effects ,  // max perturb 
		pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 		

		interpolate( 
		&(pCell->custom_data[r2_i]) , // overwrite this
		nanobio_parameters.untreated_phenotype.relative_glycolytic_effects , // unperturbed
		max_response.relative_glycolytic_effects ,  // max perturb 
		pCell->custom_data[ cell_effect_i ]  // interpolation parameter  
		); 
	}
	
	return; 
}

// output = (1-param)*input1 + param*input2 ( 0 <= param <= 1 )
void interpolate( double* pOutput , double input1, double input2, double param )
{
	*pOutput = input2; // input2
	*pOutput -= input1; // input2-input1
	*pOutput *= param; // param*(input2-input1)
	*pOutput += input1; // input1 + param*(input2-input1) 
		// == (1-param)*input1 + param*input2 
		
	return; 
}

void interpolate( std::vector<double>* pOutput , std::vector<double>& input1, std::vector<double> input2, double param )
{
	*pOutput = input2; // input2
	*pOutput -= input1; // input2-input1
	*pOutput *= param; // param*(input2-input1)
	*pOutput += input1; // input1 + param*(input2-input1) 
		// == (1-param)*input1 + param*input2 
		
	return; 
}



void apply_therapies( void ) 
{
	static double therapy_dt = 5; // update therapy every 5 minutes
	static double next_therapy_time = 0.0; 
	
	static double tolerance = 0.01 * diffusion_dt; 
	
	static int therapy_update_time_multiplier = 0; 

	// DON'T CHECK FOR THERAPY EVERY 0.01 MINUTES! 
	// Let's check it every 5 minutes instead
	
	if( PhysiCell_globals.current_time <= next_therapy_time - tolerance)
	{ return; }
	else
	{ 
		therapy_update_time_multiplier++;
		next_therapy_time = therapy_update_time_multiplier*therapy_dt; 
	}
	
	for( int k=0; k < nanobio_parameters.biofvm_indices.size() ; k++ )
	{
		double elapsed_therapy_time = PhysiCell_globals.current_time - nanobio_parameters.activation_times[k]; 
		
		int index = nanobio_parameters.biofvm_indices[k]; 
		// If therapy is started, this time is non-negative. 
		if( elapsed_therapy_time >= -tolerance )
		{
			// compute the current boundary / blood value, 
			// based on the model that the drug is globally cleared (e.g., in the renal system)
			double dose = nanobio_parameters.doses[k] * exp( -nanobio_parameters.clearance_rates[k] * elapsed_therapy_time ); 
			
			// make sure the kth condition is set to be Dirichlet 
			microenvironment.set_substrate_dirichlet_activation(index,true); 
			
			#pragma omp parallel for 
			for( int i=0 ; i < microenvironment.mesh.voxels.size() ;i++ )
			{
				if( microenvironment.mesh.voxels[i].is_Dirichlet == true )
				{ microenvironment.update_dirichlet_node( i, index, dose ); }
			}
			
		}
		
	}
	
	return; 
}

void NP_transformation::read_from_pugixml(pugi::xml_node parent_node )   
{
	NP_condition_substrate_index = xml_get_int_value( parent_node, "condition_substrate_ID" );
	
	NP_start_index = xml_get_int_value( parent_node, "start_NP_ID" );
	NP_end_index = xml_get_int_value( parent_node, "end_NP_ID" );
	
	NP_substrate_start_index = xml_get_int_value( parent_node, "start_substrate_ID" );
	NP_substrate_end_index = xml_get_int_value( parent_node, "end_substrate_ID" );
	
	condition_1 = xml_get_double_value( parent_node , "condition1" ); 
	rate_1 = xml_get_double_value( parent_node, "rate1" );
	condition_2 = xml_get_double_value( parent_node , "condition2" ); 
	rate_2 = xml_get_double_value( parent_node, "rate2" );
	
	std::cout << "substrate " << NP_substrate_start_index << " to " << NP_substrate_end_index  
		<< " varies with substrate " << NP_condition_substrate_index << std::endl 
	<< "\t" << condition_1 << " : " << rate_1 << std::endl 
	<< "\t" << condition_2 << " : " << rate_2 << std::endl << std::endl; 

	return;
};

double NP_transformation::get_transition_rate( int voxel_index )
{
	static double diff = condition_2 - condition_1; 
	
	double output = 0.0; 
	double condition = microenvironment(voxel_index)[ NP_condition_substrate_index ]; 
	double param = condition;
	param -= condition_2; 
	param /= diff; 
	if( param > 1 )
	{ param = 1.0; }
	if( param < 0.0 )
	{ param = 0.0; } 

	
	interpolate( &output , rate_1 , rate_2 , param ); 
	
	return output; 
}

void NP_transformations( double dt )
{
	#pragma omp parallel for 
	for( int k=0 ; k < microenvironment.number_of_voxels(); k++ )
	{
		if( microenvironment.is_dirichlet_node( k ) == false )
		{
			// NP1 --> NP2 
			double rate = nanobio_parameters.NP1_to_NP2.get_transition_rate( k ); 
			double temp = rate; 
			temp *= dt; 
			temp += 1.0; 
			microenvironment(k)[ nanobio_parameters.NP1_to_NP2.NP_substrate_start_index ] /= temp; 
			
			temp = rate; 
			temp *= dt; 
			temp *= microenvironment(k)[ nanobio_parameters.NP1_to_NP2.NP_substrate_start_index ]; 
			
			microenvironment(k)[ nanobio_parameters.NP1_to_NP2.NP_substrate_end_index ] += temp; 
			
			// NP2 --> NP1 
			
			rate = nanobio_parameters.NP2_to_NP1.get_transition_rate( k ); 
			temp = rate; 
			temp *= dt; 
			temp += 1.0; 
			microenvironment(k)[ nanobio_parameters.NP2_to_NP1.NP_substrate_start_index ] /= temp; 
			
			temp = rate; 
			temp *= dt; 
			temp *= microenvironment(k)[ nanobio_parameters.NP2_to_NP1.NP_substrate_start_index ]; 
			
			microenvironment(k)[ nanobio_parameters.NP2_to_NP1.NP_substrate_end_index ] += temp; 	
			
		}
	}
	
	return; 
}
