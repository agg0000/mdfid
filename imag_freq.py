#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from abc import ABC, abstractmethod

class freq_file( ABC ):

    def __init__( self ):
        self = self

    def build( self, file_name ):
        self.file_name   = file_name
        self.natom       = 0
        self.freqs       = []
        self.elements    = []
        self.freq_coords = []
        self.init_coords = []

    @abstractmethod
    def get_freq_norm( self ):
        pass

    def remove_imag( self, scale, mix ):
        min_freq = min( self.freqs )
        if min_freq > 0:
            print( "Warning: There is no imag freq!" )
            return

        if mix:
            out_array = np.zeros( ( 2, self.natom, 3 ) )

        for i, freq in enumerate( self.freqs ):
            if freq < 0:
                for n, pm in enumerate([-1, 1]):
                    if mix:
                        out_array[n] += pm * scale * self.freq_coords[i]
                    else:
                        new_coord = self.init_coords + pm * scale * self.freq_coords[i]

                        out_file = "{}{}{}".format( self.file_name.split( "." )[0], i, n )
                        self.out_imag_xyz( self.natom, self.elements, new_coord, out_file )

        if mix:
            for j in range( 2 ):
                out_file = "{}{}".format( self.file_name.split( "." )[0], j )
                
                new_coord = self.init_coords + out_array[j]
                self.out_imag_xyz( self.natom, self.elements, new_coord, out_file )

        return

    def out_xyz( self ):
        out_name = self.file_name.split( "." )[0] + "_opt"

        self.out_imag_xyz( self.natom, self.elements, self.init_coords, out_name )

    def show_vib( self ):
        vib_name = self.file_name.split( "." )[0]
        self.out_vib_xyz( self.natom, self.elements, self.init_coords, self.freqs, self.freq_coords, vib_name )

    @staticmethod
    def out_imag_xyz( natom, elements, coords, out_name ):
        if len( coords ) != natom or len( elements ) != natom:
            print( "Error: coords or elements len is not eq to natom!" )
            exit(21)

        with open("{}.xyz".format( out_name ), "w+") as xyz_file:
            xyz_file.write( "{}\n\n".format( natom ) )

            for i in range(natom):
                xyz_file.write( " {0:<6s} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f}\n".format( elements[i], coords[i] ) )

        return

    @staticmethod
    def out_vib_xyz( natom, elements, coords, frequences, freq_coords, out_name ):
        if len( coords ) != natom or len( elements ) != natom:
            print( "Error: coords or elements len is not eq to natom!" )
            exit(21)

        with open("{}_vib.xyz".format( out_name ), "w+") as xyz_file:
            for n, freq in enumerate( frequences ):
                xyz_file.write( "{}\n".format( natom ) )
                xyz_file.write( "  Frequencies: {:10.4f} cm^-1\n".format( frequences[n] ) )

                for i in range( natom ):
                    xyz_file.write( " {0:<6s} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f} {2[0]:9.4f} {2[1]:9.4f} {2[2]:9.4f}\n".format( elements[i], coords[i], freq_coords[n][i] ) )

        return

class gau_freq_file( freq_file ):

    def get_freq_norm( self ):
        with open( self.file_name ) as input_file:
            file_line = input_file.readlines()
        
        coord_line = 0
        for i, line in enumerate( file_line ):
            if "Standard orientation" in line or "Input orientation" in line:
                coord_line = i
        
            if "NAtoms" in line:
                self.natom = int( line.split()[1] )
        
        if not self.natom:
            print( "Error in read number of atoms" )
            exit(1)
        
        freq_read = False
        ele_read  = True
        freq_coords = []
        elements = []
        for i, line in enumerate( file_line ):
            if "Charge" in line and "Multiplicity" in line and ele_read:
                ele_read  = False
                for j in range( self.natom ):
                    ele = file_line[i + j + 1].split()
                    elements.append(ele[0])

                self.elements = elements

            if "Frequencies" in line:
                freq_read = True
        
                freq_line = line.split()[2:]
                line_n_freq = len( freq_line )
        
                self.freqs += [ float( freq ) for freq in freq_line ]
                for j in range( line_n_freq ):
                    coord = []
                    for k in range( self.natom ):
                        norm_coord = file_line[i + 5 + k].split()[2:]
                        coord.append( norm_coord[3 * j: 3 * j + 3] )
        
                    freq_coords.append( coord )
        
        init_coords = []
        if coord_line:
            coords = []
            for j in range( self.natom ):
                coord = file_line[coord_line + j + 5].split()
                coords.append( coord[3: 6] )

            init_coords = coords
        else:
            print( "Error in read orientation" )
            exit(3)

        if not freq_read:
            print( "Error in read frequencies" )
            exit(2)
        
        self.freq_coords = np.array( freq_coords, dtype = float )
        self.init_coords = np.array( init_coords, dtype = float )

class bdf_freq_file( freq_file ):

    def get_freq_norm( self ):
        with open( self.file_name ) as input_file:
            file_line = input_file.readlines()
        
        for line in file_line:
            if "Number of atoms" in line:
                self.natom = int( line.split()[-1] )
                break
        
        if not self.natom:
            print( "Error in read number of atoms" )
            exit(1)
        
        init_coords = []
        freq_coords = []
        freq_read   = False
        norm_mode   = True
        for i, line in enumerate( file_line ):
            if "Cartesian coordinates (Angstrom)" in line:
                coord    = []
                elements = []
                for j in range( self.natom ):
                    ele = file_line[i + j + 4].split()
                    coord.append( ele[3: 6] )

                    elements.append( ele[1] )

                init_coords   = coord
                self.elements = elements
        
            if "Frequencies" in line and norm_mode:
                freq_read = True
        
                freq_line = line.split()[1:]
                line_n_freq = len( freq_line )
        
                self.freqs += [ float( freq ) for freq in freq_line ]
                for j in range( line_n_freq ):
                    coord = []
                    for k in range( self.natom ):
                        norm_coord = file_line[i + 4 + k].split()[2:]
                        coord.append( norm_coord[3 * j: 3 * j + 3] )
        
                    freq_coords.append( coord )

            if "Results of translations and rotations" in line:
                norm_mode = False
        
        if not freq_read:
            print( "Error in read frequencies" )
            exit(2)
        
        self.freq_coords = np.array( freq_coords, dtype = float )
        self.init_coords = np.array( init_coords, dtype = float )

class orca_freq_file( freq_file ):

    def get_freq_norm( self ):
        with open( self.file_name ) as input_file:
            file_line = input_file.readlines()
        
        degrees_of_freedom = 0
        for line in file_line:
            if "The number of degrees" in line:
                degrees_of_freedom = int( line.split()[-1] )

            if "Number of atoms" in line:
                self.natom = int( line.split()[-1] )
        
        if degrees_of_freedom:
            remove_tran_rota = 3 * self.natom - degrees_of_freedom
        else:
            remove_tran_rota = 6
            if "The molecule is recognized as being linear" in line:
                remove_tran_rota = 5

            degrees_of_freedom = 3 * self.natom - remove_tran_rota

        if not self.natom:
            print( "Error in read number of atoms" )
            exit(1)
        
        init_coords = []
        freq_coords = []
        freq_read   = False
        for i, line in enumerate( file_line ):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coord    = []
                elements = []
                for j in range( self.natom ):
                    ele = file_line[i + j + 2].split()
                    coord.append( ele[1: 4] )

                    elements.append( ele[0] )

                init_coords   = coord
                self.elements = elements
        
            if "VIBRATIONAL FREQUENCIES" in line:
                freq_read = True

                ii = i + 5 + remove_tran_rota
                for j in range(degrees_of_freedom):
                    freq_line = file_line[ii + j].split()[1]
                    freq = float( freq_line )
                    self.freqs.append( freq )
        
            if "NORMAL MODES" in line:
                if remove_tran_rota == 5:
                    coord = []
                    for j in range( 3 * self.natom ):
                        norm_coord = file_line[i + 8 + j].split()[-1]
                        coord.append( norm_coord )

                    coord_array = np.array( coord, dtype = float ).reshape(self.natom, 3)
                    freq_coords.append( coord_array )

                j = i + 8 + 3 * self.natom
                while file_line[j].strip():
                    line_n_freq = len( file_line[j].split() )

                    for k in range( line_n_freq ):
                        coord = []
                        for l in range( 3 * self.natom ):
                            norm = file_line[j + 1 + l].split()
                            norm_coord = norm[k + 1]
                            coord.append( norm_coord )

                        coord_array = np.array( coord, dtype = float ).reshape(self.natom, 3)
                        freq_coords.append( coord_array )

                    j += 3 * self.natom + 1
        
        if not freq_read:
            print( "Error in read frequencies" )
            exit(2)
        
        self.freq_coords = np.array( freq_coords, dtype = float )
        self.init_coords = np.array( init_coords, dtype = float )

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option( '-s', dest = 's', type = float, default = 0.1,   help = 'scale for norm to coord' )
    parser.add_option( '-f', dest = 'f', type = str,   default = "gau", help = 'what soft output file'   )
    
    parser.add_option( '-m', dest = 'm', action = "store_true", default = False, help = "if mix imag freq" )
    parser.add_option( '-c', dest = 'c', action = "store_true", default = False, help = "only output coordinate" )
    parser.add_option( '-w', dest = 'w', action = "store_true", default = False, help = "show the vibration with xyz file in Jmol" )
    
    (options, args) = parser.parse_args()
    mix   = options.m
    soft  = options.f
    coord = options.c
    scale = options.s

    vib_show = options.w

    key_class = { "gau":  gau_freq_file(), 
                  "bdf":  bdf_freq_file(), 
                  "orca": orca_freq_file() }
    freq_get = key_class.get( soft )

    freq_get.build( args[0] )
    freq_get.get_freq_norm()

    if coord:
        freq_get.out_xyz()
    else:
        if not vib_show:
            freq_get.remove_imag( scale, mix )
        else:
            freq_get.show_vib()

