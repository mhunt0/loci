#!/bin/bash
###############################################################################
#
# Copyright 2008, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################

INSTALL_PATH=$INSTALL_DIR/$LOCI_INSTALL_DIR

echo INSTALL_PATH = $INSTALL_PATH

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib
mkdir -p $INSTALL_PATH/bin

echo Installing Library Files
cp Tools/libTools.so $INSTALL_PATH/lib
cp System/libLoci.so $INSTALL_PATH/lib
cp FVMMod/fvm_m.so $INSTALL_PATH/lib
cp ParMetis-3.1/METISLib/libmetis.so $INSTALL_PATH/lib
cp ParMetis-3.1/ParMETISLib/libparmetis.so $INSTALL_PATH/lib
cp sprng/libsprng.so $INSTALL_PATH/lib

echo Installing Loci Tools
cp lpp/lpp $INSTALL_PATH/bin
cp FVMtools/cobalt2xdr $INSTALL_PATH/bin
cp FVMtools/extract $INSTALL_PATH/bin
cp FVMtools/pb $INSTALL_PATH/bin
cp FVMtools/plot3d2xdr $INSTALL_PATH/bin
cp FVMtools/solidMesh2xdr $INSTALL_PATH/bin
cp FVMtools/xdr2cobalt $INSTALL_PATH/bin
cp FVMtools/xdr2vog $INSTALL_PATH/bin
cp FVMtools/vog2xdr $INSTALL_PATH/bin
cp FVMtools/vog2surf $INSTALL_PATH/bin
cp FVMtools/ugrid2vog $INSTALL_PATH/bin
cp FVMtools/cfd++2xdr $INSTALL_PATH/bin
cp FVMtools/vogmerge $INSTALL_PATH/bin
cp FVMtools/vogcheck $INSTALL_PATH/bin
cp FVMadapt/refmesh $INSTALL_PATH/bin
cp FVMadapt/marker $INSTALL_PATH/bin

echo cp Loci.conf comp.conf sys.conf $INSTALL_PATH
cp Loci.conf comp.conf sys.conf $INSTALL_PATH

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include
cp include/*.lh $INSTALL_PATH/include
cp include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb ; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done

mkdir -p $INSTALL_PATH/docs
mkdir -p $INSTALL_PATH/docs/1D-Diffusion
mkdir -p $INSTALL_PATH/docs/heat
mkdir -p $INSTALL_PATH/docs/Datatypes

cp Tutorial/docs/tutorial.pdf $INSTALL_PATH/docs
cp Tutorial/1D-Diffusion/Makefile $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.loci $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.lh $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.cc $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.vars $INSTALL_PATH/docs/1D-Diffusion

cp Tutorial/heat/Makefile $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.loci $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.lh $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.vog $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.vars $INSTALL_PATH/docs/heat

cp Tutorial/Datatypes/Makefile $INSTALL_PATH/docs/Datatypes
cp Tutorial/Datatypes/*.cc $INSTALL_PATH/docs/Datatypes

chmod -R a+rX $INSTALL_PATH
