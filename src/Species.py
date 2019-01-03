# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:58:38 2018

@author: Lützenkirchen, Heberling, Jara
"""
import activity_function_module as afm

class Species:
    '''
        properties:
            name            e.g. String that labels the species like CO2, H+, Ca2+ or Q, A, Whatsoever, and etc. 
        methods:
    '''
    # Constructor
    def __init__(self, name):
        self.name = name  
# Aquous species (Aq_species) is a daughter class of Species
class Aq_Species (Species):
    '''
        properties:
            f_act_coef                  e.g. Type of function for calculating the activity coefficient
            charge                      e.g. Value of charge for instance H+ is 1, Mg+2 is 2, Cl- is -1, ...
            gfw                         e.g. The gram formula weight for the species 1 gram of X equal to 1 mol of X
            ionsizepar                  e.g. Parameter that can be needed for calculating the ionic strenght
            devionsizepar               e.g. Parameter that can be needed for calculating the ionic strenght
            charge_element              e.g. Special properties that indicates that the species in consideration is charge, so what in PhreeQc would be E.
        methods:
            set_charge(charge)
            set_gfw(molar_weight)
            set_ionsizeparameter(a)
            set_deviationionsizeparameter(b)
            Set_f_activity_coefficient (String_name_func)
            it_is_charge(boolean)
            log_coefficient_activity( ionic_strength, c=0, A=0, B = 0, a = 0, b = 0, z = 0)
            dgamma_dionicstrength(actcoeff, ionic_strength, z=0, A=0)
    '''
    # Constructor
    def __init__(self, name):
        '''Sets the inputted parameter name (e.g 'NaCl', 'Mg+2', 'Cl-') into the property name.'''
        Species.__init__(self, name)
        self.f_act_coef = 'Davis'
    # Setters
    def set_charge(self, charge):
        '''Sets the inputted parameter charge (e.g 0, 2, -1) into the property charge.'''
        self.charge = charge
    
    def set_gfw(self, molar_weight):
        '''Sets the inputted parameter molar_weight (e.g 58.44, 48.61, 35.4530) into the property gfw.
           It is assumed that the unit is g/mol.
        '''
        self.gfw = molar_weight
    
    
    def set_ionsizeparameter(self, a):
        '''
            The ion size paramete is set. The units of the input must be angstroms --> 1e-8 agnstrom = 1 cm
            and it will be converted to cm
        '''
        self.ionsizepar = a*1e-8
    
    def set_deviationionsizeparameter(self, b):
        '''
            The deviation ion size paramete is set.
            I am unsure about the units of this parameter. So far I leave it like that, it might change in the future
        '''
        self.devionsizepar = b
        
    def Set_f_activity_coefficient (self, String_name_func):
        '''
        Davis, Deby-Huckel
        '''
        
        self.f_act_coef = String_name_func
        
    #Getters
    
    # Special things
    def it_is_charge(self, boolean):
        self.charge_element = boolean
        
        
        ################### Calculations ############################
        
    def log_coefficient_activity(self, ionic_strength, c=0, A=0, B = 0, a = 0, b = 0, z = 0):
        if hasattr(self, 'charge'):
            z = self.charge
        if hasattr(self, 'ionsizepar'):
            a = self.ionsizepar
        if hasattr(self, 'devionsizepar'):
            b = self.devionsizepar
        log_act_coef = afm.f_log_list(self.f_act_coef, ionic_strength, c, A, B, a, b, z)
        return log_act_coef
    
    def dgamma_dionicstrength(self, actcoeff, ionic_strength, z=0, A=0):
        '''
            returns the derivative value of the activity regarding the inoic strength (dgamma/dionic_strength)
        '''
        if hasattr(self, 'charge'):
            z = self.charge
        return afm.dgamma_dionicstrength(self.f_act_coef, actcoeff, ionic_strength, z, A)
    
# Surface Species (Surf_species) is a daughter class of Species
        
class Surf_species (Species):
    '''
        properties:
            type_sorption                               e.g. Sorption can be of different types therefore a name must be given so far the types considered in the code are: 
                                                                                                                            CCM             (Which stands for constant capacitance model - Westall 1980)
            T_solid                                     e.g. This attribute/property should only be defined to primary species. It set the value of the total sum of primary and secondary species of a solid. Component value for solid.  T = stoichiometric*primary+stoichiometric*secondary                                                                  
            sp_surf_area                                e.g. The specific surface area in m^2/g  needed to calculate T_sigma
            solid_concentration_or_grams                e.g. The specific surface area in g or g/L or g/kgwater  needed to calculate T_sigma
            C1                                          i.e. For CCM model or for the TLM
            C2                                          i.e. For the TLM
            names_Boltz_psi                             e.g. name of the electrostatic potential, necessary for assigning names in stoichiometric and component matrices
            type_relation                               i.e. it is a flag that indicates that the following primary species share the electrostatic potential with other species, e.g. Hfo_w shared it with Hfo_s
            zq                                          i.e. Charge of the solid, for instance Hfo_sOH2+ has a charge of +1.
        methods:
            set_type_sorption (type_sorpt)
            moles_component_solid(mols_or_mols_liter_Tsolid)
            specific_surface_area(Specific_Area)
            set_type_relation (species_related)
            create_name_electrostatic_potential_type (type_sorpt)
            solid_concentration(Solid_or_SolidConcentration)
            set_capacitance_1 (Capacitance_value)
            set_capacitance_2 (Capacitance_value)
            set_surface_charge (z_q)
    '''
    # Constructor
    def __init__(self, name):
        '''Sets the inputted parameter name (e.g 'SurfOH', ...) into the property name.'''
        Species.__init__(self, name)
    

    # METHODS ONLY FOR PRIMARY SPECIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    # I am unsure if that is the proper way to define the following 
    # Type
    def set_type_sorption (self, type_sorpt):
        '''
            Sets the type of sorption reaction to take into account:
                - Non-electrostatic []
                - Constant capacity [Schindler et al. (1976)]
                - diffuse_layer [Dzombak and Morel (1990)]
                - three layer exchange [Davies and Leckie (1978)]
                - dependent         --> these would be use for weak, strong, and other type of assembles. Namely, one surface potential and different phases. A species can be primary but electrostatic dependent
        '''
        self.type_sorption = type_sorpt
        self.names_Boltz_psi = self.create_name_electrostatic_potential_type (type_sorpt)
    
    def set_type_relation (self, species_related):
        '''
            This is a flag that points to another species. 
            e.g.
            Imaging just one surface with two types of sorption sites: weak and strong, both will share the same electrostatic potential but when solving the problem we need only one electrostatic potential
            Hence, one species will point to the name of the ohter. The decision of which one points to each other, is given by the user.
        '''
        self.type_relation = species_related
    
    def create_name_electrostatic_potential_type (self, type_sorpt):
        '''
        Sets the name of the boltzmanelectrostatic potential according to the type of sorption:
                - Non-electrostatic []
                - Constant capacity [Schindler et al. (1976)]
                - diffuse_layer [Dzombak and Morel (1990)]
                - three layer exchange [Davies and Leckie (1978)]
                - dependent         --> these would be use for weak, strong, and other type of assembles. Namely, one surface potential and different phases. A species can be primary but electrostatic dependent
        '''
        if type_sorpt == 'CCM':
            return self.name + '_Boltzf_psi0'               # Boltzf_psi0 stands for Boltzman factor of the psi which is equal to
        elif type_sorpt == 'DLM':
            return self.name + '_Boltzf_psi_diffuselayer'
        elif type_sorpt == 'TLM':
            return [self.name + '_Boltzf_psi_0', self.name + '_Boltzf_psi_beta', self.name + '_Boltzf_psi_diffuselayer']
        
    # In order to solve the speciation problem with sorption, the number of total moles of the component must be known.
    # These value is not the number of moles of the solid species but the number of moles of the sum of all solid species related to the primary
    # These is what in CHEPROO would be the unknown u_sorptspecies or what Westall call the T_sigma
    def moles_component_solid(self, mols_or_mols_liter_Tsolid):
        '''
            These value is not the number of moles of the solid species but the number of moles of the sum of all solid species related to the primary
            These is what in CHEPROO would be the unknown u_sorptspecies or what Westall call the T_SOH in CCM example
            It is need to solve the Newton-Raphson algorithm
            It has to do with the solid unknown or chemical part, not with the electrostatic part
            It shoudl be mols or moles/liter (mol/kgw), etc. it is somehow determined by the Newton-Raphson used, and assumpition of the user. These class does not check units.
        '''
        
        self.T_solid = float(mols_or_mols_liter_Tsolid)
    
    def specific_surface_area(self, Specific_Area):
        '''
            Parameter needed to calculate the T_sigma (or the electrostatic surface component value )for the solution of the Jacobian
            It is assume here that the units are square meter/gram
            It should be defined only for electrostatic independent primary species
        '''
        self.sp_surf_area = float(Specific_Area)
        
    def solid_concentration(self, Solid_or_SolidConcentration):
        '''
            Parameter needed to calculate the T_sigma (or the electrostatic surface component value )for the solution of the Jacobian
            It is assume here that the units are square gram or grams/l or grams/kgwater. User must be aware of the proper units
            It should be defined only for electrostatic independent primary species
        '''
        self.solid_concentration_or_grams = float(Solid_or_SolidConcentration)
        
    def set_capacitance_1 (self, Capacitance_value):
        '''
            These parameter is defined to the CCM model or to the Triple layer model.
            In each model has a different meaning. Although it is obviously a Capacitance, therefore has units.
            The units are C/m^2
        '''
        
        self.C1 = float(Capacitance_value)
        
    def set_capacitance_2 (self, Capacitance_value):
        '''
            These parameter is defined to the TLM model 
            It is obviously a Capacitance, therefore has units.
            The units are C/m^2
        '''
        
        self.C2 = float(Capacitance_value)
    
    def set_surface_charge (self, z_q):
        '''
            I am not 100% of the use of these variable, but it is reasonable to add the information if we assume that 
            Hfo_wOH2+ will have a charge of +1, or Hfo_wSO4- will have a charge of -1, also other species can have other values
            It is questionable the use but Bethe seems to use this approach, in its DLM.
        '''
        self.zq = float(z_q)