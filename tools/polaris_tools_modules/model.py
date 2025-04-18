#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.math import Math
from polaris_tools_modules.base import Model
from polaris_tools_custom.model import *


class ModelChooser:
    """The ModelChooser class provides the chosen model.
    """

    def __init__(self, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own model, add its name to the dictionary
            and write a class with its options as a derived class of class Model.

        Args:
            parse_args (ArgumentParser) : Provides all parameters chosen
            by user when executing PolarisTools.
        """
        self.parse_args = parse_args

        # Get math module
        self.math = Math()

        # dict: Dictionary with all usable models
        self.model_dict = {
            'default': Model,
            'disk': Disk,
            'sphere': Sphere,
        }
        update_model_dict(self.model_dict)

    def get_module(self):
        """Chooses model class from user input

            Note:
                Parameters set by PolarisTools overwrite preset values in the
                separate model classes.

            Returns:
                Instance of chosen model.
        """
        if self.parse_args.model_name in self.model_dict.keys():
            model = self.model_dict[self.parse_args.model_name]()
        elif self.parse_args.model_name is not None:
            raise ValueError('Model name not known! You can add a new model in model.py.')
        else:
            model = self.model_dict['default']()

        # Set user input variables
        if 'grid_type' in vars(self.parse_args).keys():
            model.update_parameter(self.parse_args.extra_parameter)
            if self.parse_args.grid_type is not None:
                model.parameter['grid_type'] = self.parse_args.grid_type
            if self.parse_args.gas_mass is not None:
                model.parameter['gas_mass'] = self.math.parse(
                    self.parse_args.gas_mass, 'mass')
            if self.parse_args.inner_radius is not None:
                model.parameter['inner_radius'] = self.math.parse(
                    self.parse_args.inner_radius, 'length')
            if self.parse_args.outer_radius is not None:
                model.parameter['outer_radius'] = self.math.parse(
                    self.parse_args.outer_radius, 'length')
            if self.parse_args.z_max is not None:
                model.cylindrical_parameter['z_max'] = self.math.parse(
                    self.parse_args.z_max, 'length')
            if self.parse_args.n_r is not None:
                model.spherical_parameter['n_r'] = self.parse_args.n_r
                model.cylindrical_parameter['n_r'] = self.parse_args.n_r
            if self.parse_args.n_ph is not None:
                model.spherical_parameter['n_ph'] = self.parse_args.n_ph
                model.cylindrical_parameter['n_ph'] = self.parse_args.n_ph
            if self.parse_args.n_th is not None:
                model.spherical_parameter['n_th'] = self.parse_args.n_th
            if self.parse_args.n_z is not None:
                model.cylindrical_parameter['n_z'] = self.parse_args.n_z
            if self.parse_args.sf_r is not None:
                model.spherical_parameter['sf_r'] = self.parse_args.sf_r
                model.cylindrical_parameter['sf_r'] = self.parse_args.sf_r
            if self.parse_args.sf_ph is not None:
                model.spherical_parameter['sf_ph'] = self.parse_args.sf_ph
                model.cylindrical_parameter['sf_ph'] = self.parse_args.sf_ph
            if self.parse_args.sf_th is not None:
                model.spherical_parameter['sf_th'] = self.parse_args.sf_th
            if self.parse_args.sf_z is not None:
                model.cylindrical_parameter['sf_z'] = self.parse_args.sf_z
        elif 'distance' in vars(self.parse_args).keys():
            if self.parse_args.distance is not None:
                model.parameter['distance'] = self.math.parse(
                    self.parse_args.distance, 'length')
    
        # Set the grid extent if global extent is set
        if model.parameter['grid_type'] == 'octree' and model.parameter['outer_radius'] is not None:
            model.octree_parameter['sidelength'] = 2. * \
                model.parameter['outer_radius']
        elif model.parameter['grid_type'] == 'spherical':
            if model.parameter['inner_radius'] is not None:
                model.spherical_parameter['inner_radius'] = model.parameter['inner_radius']
            if model.parameter['outer_radius'] is not None:
                model.spherical_parameter['outer_radius'] = model.parameter['outer_radius']
        elif model.parameter['grid_type'] == 'cylindrical':
            if model.parameter['inner_radius'] is not None:
                model.cylindrical_parameter['inner_radius'] = model.parameter['inner_radius']
            if model.parameter['outer_radius'] is not None:
                model.cylindrical_parameter['outer_radius'] = model.parameter['outer_radius']
                if model.cylindrical_parameter['z_max'] is None:
                    model.cylindrical_parameter['z_max'] = model.parameter['outer_radius']

        return model


class Disk(Model):
    """Power-law or viscous accretion disk model.
    Andrews et al. (2009), ApJ 700, 1502
    Kwon et al. (2015), ApJ 808, 102
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['gas_mass'] = 1e-3 * self.math.const['M_sun']
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100. * self.math.const['au']
        
        # In the case of a spherical grid
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 181
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03
        # sf_th = 1 is sinus; sf_th > 1 is exp with step width sf_th; rest is linear
        self.spherical_parameter['sf_th'] = 1.0
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 100
        self.cylindrical_parameter['n_z'] = 181
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.03
        # sf_z = -1 is using scale height; sf_z = 1 is sinus;
        # sf_z > 1 is exp with step width sf_z and rest is linear
        self.cylindrical_parameter['sf_z'] = -1
        
        # Default disk parameter
        self.parameter['ref_radius'] = 100. * self.math.const['au']
        self.parameter['ref_scale_height'] = 10. * self.math.const['au']
        self.parameter['beta'] = 1.1
        self.parameter['alpha'] = 3.0 * (self.parameter['beta'] - 0.5)
        self.parameter['tapered_gamma'] = None

        self.custom_parameter_list = [
            'ref_radius',
            'ref_scale_height',
            'beta',
            'alpha',
            'tapered_gamma',
        ]

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.custom_parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f'  - Info: {param} updated to {self.parameter[param]}')

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.default_disk_density(self.position,
                                                     inner_radius=self.parameter['inner_radius'],
                                                     outer_radius=self.parameter['outer_radius'],
                                                     ref_radius=self.parameter['ref_radius'],
                                                     ref_scale_height=self.parameter['ref_scale_height'],
                                                     alpha=self.parameter['alpha'], beta=self.parameter['beta'],
                                                     tapered_gamma=self.parameter['tapered_gamma'])
        return gas_density

    def get_scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        scale_height = self.math.default_disk_scale_height(radius,
                                                           ref_radius=self.parameter['ref_radius'],
                                                           ref_scale_height=self.parameter['ref_scale_height'],
                                                           beta=self.parameter['beta'])
        return scale_height


class Sphere(Model):
    """A sphere model with constant density
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the sphere model
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100. * self.math.const['au']
        
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 91
        self.spherical_parameter['n_ph'] = 91
        self.spherical_parameter['sf_r'] = 1.03
        self.parameter['gas_mass'] = 1e-4 * self.math.const['M_sun']
        
        self.parameter['mag_field_geometry'] = 'toroidal'
        self.parameter['mag_field_strength'] = 1e-10

        self.custom_parameter_list = [
            'mag_field_geometry',
            'mag_field_strength'
        ]

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.const_sphere_density(self.position,
                                                     outer_radius=self.parameter['outer_radius'],
                                                     inner_radius=self.parameter['inner_radius'])
        return gas_density

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        if self.parameter['mag_field_geometry'] == 'toroidal':
            magnetic_field = self.math.toroidal_mag_field(
                position=self.position, mag_field_strength=self.parameter['mag_field_strength'])
        elif self.parameter['mag_field_geometry'] == 'vertical':
            magnetic_field = self.math.simple_mag_field(
                mag_field_strength=self.parameter['mag_field_strength'], axis='z')
        elif self.parameter['mag_field_geometry'] == 'radial':
            magnetic_field = self.math.radial_mag_field(
                mag_field_strength=self.parameter['mag_field_strength'], position=self.position)
        else:
            magnetic_field = [0, 0, 0]
            print(f'  - Warning: invalid magnetic field geometry: {self.parameter["mag_field_geometry"]}')
        return magnetic_field

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        model parameter that depend on other parameter.
        """
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.custom_parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f'  - Info: {param} updated to {self.parameter[param]}')
