# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:00:16 2018

@author: DaniJ

This module is supossed to contain the algorithms and information of Chemical speciation plus sorption.
It is a daughter of Database_SC but it can be used without a database. 
[If feasible (question of time), I will keep it apart]
"""

from Database_SC import Database_SC
import numpy as np

class ChemSys_Surf (Database_SC):
    '''
        ChemSys is a daughter class from Database_SC which is a daughter class of Database. Hence, they depend on these parameters.        
            #Note for myself and other contributors, if you add or delete properties or methods of the class, documeted it here. Otherwise, it is a little caos (regarding my own experience)
            properties:
                Faraday_constant
                temperature
                dielectric_constant
                A_activitypar
                B_activitypar
            methods:
                set_Faraday_constant
                set_temperature
                set_dielectric_constant
                calculate_dielectric_constant
                calculate_A_activitypar
                calculate_B_activitypar
                
        NOTE: Remark that ChemSys_Surf is a daughter class from Database_SC. Therefore, in order to create the pseudo S matrix (The stoichiometric matrix that does not contain the surface potential as unknown). Methods like ...
                ... set_names_aq_primary_species (names_aq_pri_sp), set_names_aq_secondary_species (names_aq_sec_sp), set_names_sorpt_primary_species (names_sorpt_pri_sp), set_names_sorpt_secondary_species (names_sorpt_sec_sp), set_aq_list_pri_class (list_aq_pri_sp), ...
                ... set_aq_list_sec_class (list_aq_sec_sp) can be used and must be used. However, it has to be check that the input given is in accordance with the own system, that can be done by ???????????
    '''
    # Constructor
    def __init__(self):
        self.Faraday_constant = 96485.3328959 # C/mol
        self.temperature = (273.15+25)   # It assumed that initially we are at T=25°C and we assume atmospheric pressure for dielectric and other constants
        self.universal_gas_constant = 8.314472  # J/(K*mol)
        self.calculate_dielectric_constant()
        pass
    # Creating first pseudoS
    
    #Setters
    
    # Faraday constant    
    def set_Faraday_constant (self, new_value):
        '''
            The Faraday constant is instantiated with the class. The Faraday constant has the value 96485.33289(59) C mol−1 [Obtained from WIKI: https://en.wikipedia.org/wiki/Faraday_constant]
            The constant is the relationship between the elementary charge or the magnitude of the charge of an electron ['e'] and the Avogrado constant (The number of particles in a mol) [NA]
            F = e * NA
            
            e ≈ 1.60217662×10−19 C
            NA ≈ 6.02214086×1023 mol−1
            
            Note of one of the authors: I do not think that it should be modified but maybe someone what to play with the value
        '''
        
        self.Faraday_constant = new_value
    
    # Temperature
    def set_temperature(self, new_T):
        '''
            Temperature is supposed to be given in kelvins.
        '''
        
        self.temperature = new_T
     
    # Universal gas constant
    def set_universal_gas_constant (self, r_value):
        '''
            Set the universal gas constant
        '''
        self.universal_gas_constant = r_value
        
     # dielectric constant
    def set_dielectric_constant (self, e_c):
         '''
             Set the dielectric constant of water
         '''
        self.dielectric_constant = e_c
        
        
        
        
        
        
        
    # Calculations

    # Dielectric constant of water    
    def calculate_dielectric_constant(self):
        '''
            Calculates the dielectric constant
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
        '''
        self.dielectric_constant = 2727.586 + 0.6224107*self.temperature - 466.9151*np.log(self.temperature) - (52000.87/self.temperature)
        
    def calculate_A_activitypar (self):
        '''
            Calculates the parameter A of the Debye Hueckel equation
            The units are supossed to be kg^(1/2)/mol^(1/2)
            Actually if you want the L/mol is possible to divide by the square of the density to obtain such value
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
        '''
        A = 1.82483e6*np.sqrt(self.waterdensity)
        B = (self.temperature*self.dielectric_constant)**(3/2)
        self.A_activitypar = A/B
        
    def calculate_B_activitypar (self):
        '''
            Calculates the parameter A of the Debye Hueckel equation
            The units are supossed to be kg^(1/2)/mol^(1/2)*cm
            Actually if you want the L/mol is possible to divide by the square of the density to obtain such value
            The extra-calculations are baased on the book section 1.1.2.6 Calculation of activity coefficient -- Groundwater Geochemistry --- Broder J. Merkel, Britta Planer-Friedrich
            Here the equation is a bit different than that given in the book. The book takes the equation from 
            Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures; II, Debye Huckel parameters for activity coefficients and relative partial molal properties
            The differences is 10-8 and is related to the fact that they uses angstroms instead of cm
        '''
        A = 50.29158649e8*np.sqrt(self.waterdensity)
        B = np.sqrt(self.temperature*self.dielectric_constant)
        self.B_activitypar = A/B