/*
Copyright (c) 2017 Lars Pastewka
              2022 Hannes Holey

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

Linux: g++ -std=c++11 ncsplit.cpp -o ncsplit `nc-config --libs --cxx4libs`
Linux, static linking: g++ -std=c++11 ncsplit.cpp -o ncsplit -Wl,-Bstatic `nc-config --libs --cxx4libs` -Wl,-Bdynamic -lhdf5 -lhdf5_hl -lcurl
MaxOS: c++ -std=c++11 -stdlib=libc++ ncsplit.cpp -o ncsplit `nc-config --libs --cxx4libs`
*/

#include <stdio.h>

#include <cmath>

#include <netcdf>

using namespace netCDF;
using namespace exceptions;

const size_t chunksize = 1000*1000;
const double tol = 1e-12;

void progressbar(int x, int maxx, int len=60, const char *prefix=NULL)
{
    double xx = double(x)/(maxx-1);
    int xxi = xx*len;
    if (xxi > len)  xxi = len;
    if (prefix)
        fprintf(stderr, "\r%s |", prefix);
    else
        fprintf(stderr, "\r|");
    for (int i = 0; i < xxi; ++i)  fprintf(stderr, "#");
    for (int i = 0; i < len-xxi; ++i)  fprintf(stderr, "-");
    fprintf(stderr, "| %i%% (%i/%i)", int(xx*100), x, maxx);
}

template <typename T>
void copyvar(const char *name, const NcVar &in_var,
             std::vector<size_t> &in_start, const NcVar &out_var,
             std::vector<size_t> &out_start, std::vector<size_t> &count,
             size_t natoms, T *buffer, int atom_dim=1)
{
    size_t chunkstart = 0;
    while (chunkstart < natoms) {
        progressbar(chunkstart, natoms, 60, name);
        in_start[atom_dim] = chunkstart;
        out_start[atom_dim] = chunkstart;
        count[atom_dim] = std::min(chunksize, natoms-chunkstart);

        in_var.getVar(in_start, count, buffer);
        out_var.putVar(out_start, count, buffer);

        chunkstart += chunksize;
    }
}

void split(char *infn, char *outpref, int nchunks)
{
    NcFile innc(infn, NcFile::read);

    const NcDim &cell_spatial = innc.getDim("cell_spatial");
    const NcDim &cell_angular = innc.getDim("cell_angular");
    const NcDim &frame = innc.getDim("frame");
    const NcDim &atom = innc.getDim("atom");
    const NcDim &spatial = innc.getDim("spatial");

    const NcVar &cell_lengths = innc.getVar("cell_lengths");
    const NcVar &cell_angles = innc.getVar("cell_angles");
    const NcVar &cell_origin = innc.getVar("cell_origin");
    const NcVar &ids = innc.getVar("id");
    NcVar types = innc.getVar("type");
    if (types.isNull()) types = innc.getVar("atom_types");
    const NcVar &coordinates = innc.getVar("coordinates");
    const NcVar &velocities = innc.getVar("velocities");

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
    if (!ids.isNull()) {
        const auto &ids_dims = ids.getDims();
        if (ids_dims.size() != 2 || ids_dims[0] != frame || ids_dims[1] != atom) {
            printf("'id' variable must have dimensions '(frame, atom)'.\n");
            return;
        }
    }
    if (!types.isNull()) {
        const auto &types_dims = types.getDims();
        if (types_dims.size() != 2 || types_dims[0] != frame ||
            types_dims[1] != atom) {
            printf("'type' variable must have dimensions '(frame, atom)'.\n");
            return;
        }
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

    for (int i = 0; i < nchunks; i++){

        std::string suffix = ".nc";
        auto zero_padded_num = std::string(3 - std::to_string(i).length(), '0') + std::to_string(i);
        std::string outfn{outpref};
        outfn += zero_padded_num;
        outfn += suffix;

        std::cout << outfn << std::endl;

        NcFile outnc(outfn, NcFile::newFile, NcFile::classic64);
        outnc.putAtt("Conventions", "AMBER");
        outnc.putAtt("ConventionVersion", "1.0");
        outnc.putAtt("program", "ncsplit");
        outnc.putAtt("programVersion", "0.1");

        const NcDim &out_cell_spatial = outnc.addDim("cell_spatial", cell_spatial.getSize());
        const NcDim &out_cell_angular = outnc.addDim("cell_angular", cell_angular.getSize());
        const NcDim &out_frame = outnc.addDim("frame");
        const NcDim &out_atom = outnc.addDim("atom", atom.getSize());
        const NcDim &out_spatial = outnc.addDim("spatial", spatial.getSize());

        const NcVar &out_cell_lengths = outnc.addVar("cell_lengths", cell_lengths.getType(), {out_frame, out_cell_spatial});
        const NcVar &out_cell_angles = outnc.addVar("cell_angles", cell_angles.getType(), {out_frame, out_cell_angular});
        const NcVar &out_cell_origin = outnc.addVar("cell_origin", cell_origin.getType(), {out_frame, out_cell_spatial});
        NcVar out_ids;
        if (!ids.isNull()) {
            out_ids = outnc.addVar("id", ids.getType(), {out_frame, out_atom});
        }
        NcVar out_types;
        if (!types.isNull()) {
            out_types = outnc.addVar("atom_types", types.getType(), {out_frame, out_atom});
        }
        const NcVar &out_coordinates = outnc.addVar("coordinates", coordinates.getType(), {out_frame, out_atom, out_spatial});
        NcVar out_velocities;
        if (!velocities.isNull()) {
          out_velocities = outnc.addVar("velocities", velocities.getType(), {out_frame, out_atom, out_spatial});
        }

        double vec[3];
        std::vector<int> intScalarChunk(chunksize);
        std::vector<double> doubleVectorChunk(chunksize*3);

        size_t out_frame_index = 0;

        size_t start {i * nframes / nchunks};
        size_t end {(i + 1) * nframes / nchunks};

        for (size_t frame_index = start; frame_index < end; frame_index++) {
            printf("\33[2K\r=== FRAME %zu ===\n", frame_index);

            std::vector<size_t> start {frame_index, 0, 0};
            std::vector<size_t> out_start {out_frame_index, 0, 0};
            std::vector<size_t> count {1, 3, 3};

            cell_lengths.getVar(start, count, vec);
            out_cell_lengths.putVar(out_start, count, vec);
            cell_angles.getVar(start, count, vec);
            out_cell_angles.putVar(out_start, count, vec);
            cell_origin.getVar(start, count, vec);
            out_cell_origin.putVar(out_start, count, vec);

            if (!ids.isNull())
                copyvar("id", ids, start, out_ids, out_start, count, natoms,
                        intScalarChunk.data());
            if (!types.isNull())
                copyvar("atom_types", types, start, out_types, out_start, count,
                        natoms, intScalarChunk.data());
            copyvar("coordinates", coordinates, start, out_coordinates,
                    out_start, count, natoms, doubleVectorChunk.data());
            if (!velocities.isNull())
                copyvar("velocities", velocities, start, out_velocities,
                        out_start, count, natoms, doubleVectorChunk.data());

            out_frame_index++;

        }
    }
}

void syntax()
{
    printf("Syntax: ncsplit <NetCDF-input-file> <every> <NetCDF-output-file>\n"
           "Keep only every n-th frame in an AMBER NetCDF trajectory file.");
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        syntax();
        return 2;
    }
    try {
        split(argv[1], argv[3], atoi(argv[2]));
    } catch(NcException& e) {
        std::cout << "SPLIT FAILED" << std::endl
                  << "---------------" << std::endl
                  << e.what() << std::endl;
        return 2;
    }
    return 0;
}
