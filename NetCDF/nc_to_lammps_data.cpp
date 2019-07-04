/*
Copyright (c) 2017 Lars Pastewka

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Compile with:
Linux: g++ -std=c++11 nc_to_lammps_data.cpp -o nc_to_lammps_data `nc-config --libs --cxx4libs`
Linux, static linking: g++ -std=c++11 nc_to_lammps_data.cpp -o nc_to_lammps_data -Wl,-Bstatic `nc-config --libs --cxx4libs` -Wl,-Bdynamic -lhdf5 -lhdf5_hl -lcurl
MaxOS: c++ -std=c++11 -stdlib=libc++ nc_to_lammps_data.cpp -o nc_to_lammps_data `nc-config --libs --cxx4libs`
*/

#include <stdio.h>

#include <cmath>

#include <netcdf>

using namespace netCDF;
using namespace exceptions;

const size_t chunksize = 1000*1000;
const double tol = 1e-12;

void progressbar(int x, int maxx, int len=60)
{
    double xx = double(x)/(maxx-1);
    int xxi = xx*len;
    if (xxi > len)  xxi = len;
    fprintf(stderr, "\r|");
    for (int i = 0; i < xxi; ++i)  fprintf(stderr, "#");
    for (int i = 0; i < len-xxi; ++i)  fprintf(stderr, "-");
    fprintf(stderr, "| %i%% (%i/%i)", int(xx*100), x, maxx);
}

void nc2data(char *ncfn, char *datafn, int frame_index)
{
    NcFile nc(ncfn, NcFile::read);

    const NcDim &cell_spatial = nc.getDim("cell_spatial"); 
    const NcDim &cell_angular = nc.getDim("cell_angular"); 
    const NcDim &frame = nc.getDim("frame");
    const NcDim &atom = nc.getDim("atom");
    const NcDim &spatial = nc.getDim("spatial"); 

    const NcVar &cell_lengths = nc.getVar("cell_lengths");
    const NcVar &cell_angles = nc.getVar("cell_angles");
    const NcVar &cell_origin = nc.getVar("cell_origin");
    const NcVar &ids = nc.getVar("id");
    NcVar types = nc.getVar("type");
    if (types.isNull()) types = nc.getVar("atom_types");
    const NcVar &coordinates = nc.getVar("coordinates");
    const NcVar &velocities = nc.getVar("velocities");

    if (cell_spatial.getSize() != 3) {
        printf("'cell_spatial' dimension differs from 3.\n");
        return;
    }
    if (cell_angular.getSize() != 3) {
        printf("'cell_angular' dimension differs from 3.\n");
        return;
    }
    if (spatial.getSize() != 3) {
        printf("'spatial' dimension differs from 3.\n");
        return;
    }

    const auto &cell_lengths_dims = cell_lengths.getDims();
    if (cell_lengths_dims.size() != 2 || cell_lengths_dims[0] != frame ||
        cell_lengths_dims[1] != cell_spatial) {
        printf("'cell_lengths' variable must have dimensions "
               "'(frame, cell_spatial)'.\n");
        return;
    }
    const auto &cell_angles_dims = cell_angles.getDims();
    if (cell_angles_dims.size() != 2 || cell_angles_dims[0] != frame ||
        cell_angles_dims[1] != cell_angular) {
        printf("'cell_angles' variable must have dimensions "
               "'(frame, cell_angular)'.\n");
        return;
    }
    const auto &ids_dims = ids.getDims();
    if (ids_dims.size() != 2 || ids_dims[0] != frame || ids_dims[1] != atom) {
        printf("'id' variable must have dimensions '(frame, atom)'.\n");
        return;
    }
    const auto &types_dims = types.getDims();
    if (types_dims.size() != 2 || types_dims[0] != frame ||
        types_dims[1] != atom) {
        printf("'type' variable must have dimensions '(frame, atom)'.\n");
        return;
    }
    const auto &coordinates_dims = coordinates.getDims();
    if (coordinates_dims.size() != 3 || coordinates_dims[0] != frame ||
        coordinates_dims[1] != atom || coordinates_dims[2] != spatial) {
        printf("'coordinates' variable must have dimensions "
               "'(frame, atom, spatial)'.");
        return;
    }
    if (!velocities.isNull()) {
        const auto &velocities_dims = velocities.getDims();
        if (velocities_dims.size() != 3 || velocities_dims[0] != frame ||
            velocities_dims[1] != atom || velocities_dims[2] != spatial) {
            printf("'velocities' variable must have dimensions "
                   "'(frame, atom, spatial)'.");
            return;
        }
    }

    size_t nframes = frame.getSize();
    size_t natoms = atom.getSize();

    if (frame_index < 0) frame_index = nframes-frame_index;

    std::vector<int> intScalarChunk1(chunksize), intScalarChunk2(chunksize);
    std::vector<double> doubleVectorChunk(chunksize*3);

    FILE *f = fopen(datafn, "w");
    fprintf(f, "Converted from %s by nc2data\n\n", ncfn);
    fprintf(f, "%li  atoms\n", natoms);

    std::vector<size_t> start {size_t(frame_index), 0, 0};
    std::vector<size_t> count {1, 3, 3};

    int natomtypes = 0;
    printf("Determining number of atom types...\n");
    size_t chunkstart = 0;
    while (chunkstart < natoms) {
        progressbar(chunkstart, natoms);
        start[1] = chunkstart;
        count[1] = std::min(chunksize, natoms-chunkstart);
        types.getVar(start, count, intScalarChunk1.data());

        for (int i = 0; i < count[1]; ++i) {
            natomtypes = std::max(natomtypes, intScalarChunk1[i]);
        }

        chunkstart += chunksize;
    }
    progressbar(natoms, natoms);
    printf("\n");
    fprintf(f, "%i  atom types\n", natomtypes);

    double l[3], a[3];
    start[1] = 0;
    count[1] = 3;
    cell_lengths.getVar(start, count, l);
    cell_angles.getVar(start, count, a);

    double xlo = 0.0, ylo = 0.0, zlo = 0.0, xhi, yhi, zhi, xy, xz, yz;

    xhi = l[0];
    xy = std::cos(a[2]*M_PI/180) * l[1];
    yhi = std::sin(a[2]*M_PI/180) * l[1];
    xz = std::cos(a[1]*M_PI/180) * l[2];
    yz = (l[1] * l[2] * std::cos(a[0]*M_PI/180) - xy * xz) / yhi;
    zhi = std::sqrt(l[2]*l[2] - xz*xz - yz*yz);

    if (!cell_origin.isNull()) {
        cell_origin.getVar(start, count, a);
        xlo = a[0];
        ylo = a[1];
        zlo = a[2];
    }
    xhi += xlo;
    yhi += ylo;
    zhi += zlo;

    if (std::abs(l[0]) < tol || std::abs(l[1]) < tol || std::abs(l[2]) < tol) {
        double sxlo, sylo, szlo, sxhi, syhi, szhi;
        printf("Determining shrink-wrapped box...\n");
        start[1] = 0;
        count[1] = 1;
        double first_coordinate[3];
        types.getVar(start, count, first_coordinate);
        sxlo = sxhi = first_coordinate[0];
        sylo = syhi = first_coordinate[1];
        szlo = szhi = first_coordinate[2];
        size_t chunkstart = 0;
        while (chunkstart < natoms) {
            progressbar(chunkstart, natoms);
            start[1] = chunkstart;
            count[1] = std::min(chunksize, natoms-chunkstart);
            coordinates.getVar(start, count, doubleVectorChunk.data());
    
            for (int i = 0; i < count[1]; ++i) {
                sxlo = std::min(sxlo, doubleVectorChunk[3*i]);
                sylo = std::min(sylo, doubleVectorChunk[3*i+1]);
                szlo = std::min(szlo, doubleVectorChunk[3*i+2]);
                sxhi = std::max(sxhi, doubleVectorChunk[3*i]);
                syhi = std::max(syhi, doubleVectorChunk[3*i+1]);
                szhi = std::max(szhi, doubleVectorChunk[3*i+2]);
            }
    
            chunkstart += chunksize;
        }
        progressbar(natoms, natoms);
        printf("\n");
        if (std::abs(l[0]) < tol) {
            xlo = sxlo;
            xhi = sxhi;
        }
        if (std::abs(l[1]) < tol) {
            ylo = sylo;
            yhi = syhi;
        }
        if (std::abs(l[2]) < tol) {
            zlo = szlo;
            zhi = szhi;
        }
    }

    fprintf(f, "%f %f  xlo xhi\n", xlo, xhi);
    fprintf(f, "%f %f  ylo yhi\n", ylo, yhi);
    fprintf(f, "%f %f  zlo zhi\n", zlo, zhi);
    if (std::abs(xy) > tol || std::abs(xz) > tol || std::abs(yz) > tol) {
        fprintf(f, "%f %f %f  xy xz yz\n", xy, xz, yz);
    }
    fprintf(f, "\n\n");

    fprintf(f, "Atoms\n\n");
    printf("Writing coordinates...\n");
    chunkstart = 0;
    while (chunkstart < natoms) {
        progressbar(chunkstart, natoms);
        start[1] = chunkstart;
        count[1] = std::min(chunksize, natoms-chunkstart);
        ids.getVar(start, count, intScalarChunk1.data());
        types.getVar(start, count, intScalarChunk2.data());
        coordinates.getVar(start, count, doubleVectorChunk.data());

        for (int i = 0; i < count[1]; ++i) {
            fprintf(f, "%i %i %f %f %f\n", intScalarChunk1[i],
                    intScalarChunk2[i], doubleVectorChunk[3*i],
                    doubleVectorChunk[3*i+1], doubleVectorChunk[3*i+2]);
        }

        chunkstart += chunksize;
    }
    progressbar(natoms, natoms);
    printf("\n");

    if (!velocities.isNull()) {
        fprintf(f, "\n\nVelocities\n\n");
        printf("Writing velocities...\n");
        chunkstart = 0;
        while (chunkstart < natoms) {
            progressbar(chunkstart, natoms);
            start[1] = chunkstart;
            count[1] = std::min(chunksize, natoms-chunkstart);
            ids.getVar(start, count, intScalarChunk1.data());
            velocities.getVar(start, count, doubleVectorChunk.data());
    
            for (int i = 0; i < count[1]; ++i) {
                fprintf(f, "%i %f %f %f\n", intScalarChunk1[i],
                        doubleVectorChunk[3*i], doubleVectorChunk[3*i+1],
                        doubleVectorChunk[3*i+2]);
            }
    
            chunkstart += chunksize;
        }
        progressbar(natoms, natoms);
        printf("\n");   
    }

    fclose(f);
}

void syntax()
{
    printf("Syntax: nc2data <NetCDF-file> <frame-number> <data-file>\n"
           "Convert AMBER NetCDF trajectory to LAMMPS data file.\n");
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        syntax();
        return 2;
    }
    try {
        nc2data(argv[1], argv[3], atoi(argv[2]));
    } catch(NcException& e) {
        std::cout << "CONVERSION FAILED" << std::endl
                  << "-----------------" << std::endl
                  << e.what() << std::endl;
        return 2;
    }
    return 0;
}
