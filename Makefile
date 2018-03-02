#*********************************************************************
#    Makefile (version 0.2)                                          
#                                                                    
#    A simple Makefile to make the LAkernel package and its programs.
#
#    Reference:
#    H. Saigo, J.-P. Vert, T. Akutsu and N. Ueda, "Protein homology 
#    detection using string alignment kernels", Bioinformatics, 
#    vol.20, p.1682-1689, 2004.
#                                                                    
#                                                                    
#    Copyright 2003 Jean-Philippe Vert                                 
#
#    This file is part of LAkernel.
#
#    LAkernel is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    Foobar is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Foobar; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# *********************************************************************/

#
# The settings to use for the gcc compiler.  If you don't have gcc,
# just set the CC value to your compiler and make sure optimization
# is turned on (it significantly affects the running time).
#
CC=gcc
CFLAGS= -lm -g -Ofast -Wall -Wshadow -lstdc++

all: LZWkernel_all_vs_all LZWdist_all_vs_all

LZWkernel_all_vs_all: lzw-kernel_all_vs_all.o seqio.o
	$(CC)  -o lzw-kernel_all_vs_all  lzw-kernel_all_vs_all.o  seqio.o $(CFLAGS)

LZWdist_all_vs_all:  lzw-dist_all_vs_all.o seqio.o
	$(CC)  -o lzw-dist_all_vs_all  lzw-dist_all_vs_all.o  seqio.o $(CFLAGS)


LZWkernel_all_vs_all.o:  seqio.h
LZWdist_all_vs_all.o:  seqio.h
