#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from abc import ABC, abstractmethod

class freq_file( ABC ):

    def __init__( self ):
        self = self

    def build( self, file_name, scale, mix ):
        self.mix         = mix
        self.scale       = scale
        self.file_name   = file_name
        self.natom       = 0
        self.freqs       = []
        self.elements    = []
        self.freq_coords = []
        self.init_coords = []

    @abstractmethod
    def get_freq_norm( self ):
        pass

    def remove_imag( self ):
        min_freq = min( self.freqs )
        if min_freq > 0:
            print( "Warning: There is no imag freq!" )
            return

        for i, freq in enumerate( self.freqs ):
            if freq < 0:
                for n, pm in enumerate([-1, 1]):
                    new_coord = self.init_coords + pm * self.scale * self.freq_coords[i]
                    out_file = "{}{}{}".format( self.file_name.split( "." )[0], i, n )

                    self.out_xyz( self.natom, self.elements, new_coord, out_file )

        return

    @staticmethod
    def out_xyz( natom, elements, coords, out_name ):
        if len( coords ) != natom:
            print( "Error: coords len is not eq to natom!" )
            exit(21)

        with open("{}.xyz".format( out_name ), "w+") as xyz_file:
            xyz_file.write( "{}\n\n".format( natom ) )

            for i in range(natom):
                xyz_file.write( " {0:<7s} {1[0]:10.5f} {1[1]:10.5f} {1[2]:10.5f}\n".format( elements[i], coords[i] ) )

        return

class gau_freq_file( freq_file ):

    def get_freq_norm( self ):
        with open( self.file_name ) as input_file:
            file_line = input_file.readlines()
        
        for line in file_line:
            if "NAtom" in line:
                self.natom = int( line.split()[1] )
                break
        
        if not self.natom:
            print( "error in read number of atoms" )
            exit(1)
        
        first = True
        freq_read   = False
        init_coords = []
        freq_coords = []
        for i, line in enumerate( file_line ):
            if first and "Charge " in line and "Multiplicity " in line:
                first = False
                for j in range( self.natom ):
                    ele = file_line[i + j + 1].split()
                    init_coords.append( ele[1: 4] )

                    self.elements.append( ele[0] )
        
            if "Frequencies" in line:
                freq_read = True
        
                freq_line = line.split()[2:]
                line_n_freq = len( freq_line )
        
                self.freqs += [ float( freq ) for freq in freq_line ]
                for j in range( line_n_freq ):
                    coord = []
                    for k in range( self.natom ):
                        norm_coord = file_line[i + 8 + k].split()[2:]
                        coord.append( norm_coord[3 * j: 3 * j + 3] )
        
                    freq_coords.append( coord )
        
        if not freq_read:
            print( "error in read frequencies" )
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
            print( "error in read number of atoms" )
            exit(1)
        
        first = True
        freq_read   = False
        init_coords = []
        freq_coords = []
        for i, line in enumerate( file_line ):
            if first and "GEOMETRY" == line.strip():
                first = False
                k = 2 if "Read" in file_line[i + 1] else 1
                for j in range( self.natom ):
                    ele = file_line[i + j + k].split()
                    init_coords.append( ele[1: 4] )

                    self.elements.append( ele[0] )
        
            if "Frequencies" in line:
                freq_read = True
        
                freq_line = line.split()[2:]
                line_n_freq = len( freq_line )
        
                self.freqs += [ float( freq ) for freq in freq_line ]
                for j in range( line_n_freq ):
                    coord = []
                    for k in range( self.natom ):
                        norm_coord = file_line[i + 4 + k].split()[2:]
                        coord.append( norm_coord[3 * j: 3 * j + 3] )
        
                    freq_coords.append( coord )
        
        if not freq_read:
            print( "error in read frequencies" )
            exit(2)
        
        self.freq_coords = np.array( freq_coords, dtype = float )
        self.init_coords = np.array( init_coords, dtype = float )

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option( '-s', dest = 's', type = float, default = 0.1,   help = 'scale for norm to coord' )
    parser.add_option( '-f', dest = 'f', type = str,   default = "gau", help = 'what soft output file'   )
    
    parser.add_option( '-m', dest = 'm', action = "store_true", default = False, help = "if mix imag freq" )
    
    (options, args) = parser.parse_args()
    mix   = options.m
    soft  = options.f
    scale = options.s

    key_class = { "gau": gau_freq_file(), "bdf": bdf_freq_file() }
    freq_get = key_class.get( soft )

    freq_get.build( args[0], scale, mix )
    freq_get.get_freq_norm()

    freq_get.remove_imag()

