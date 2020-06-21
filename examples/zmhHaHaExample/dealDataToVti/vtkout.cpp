#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <cstdlib>

using namespace plb;

typedef double T;

void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<float>& geometry, std::string valueTag)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);
    std::auto_ptr<MultiScalarField3D<float> > slice = generateMultiScalarField<float>(geometry, sliceBox);
    plb_ifstream geometryFile(fNameIn.c_str());
    for (plint iX=0; iX<nx; ++iX) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << fNameIn << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryFile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
    }

    {
        VtkImageOutput3D<T> vtkOut(valueTag, 1.0);
        vtkOut.writeData<float>(*copyConvert<float,T>(geometry, geometry.getBoundingBox()), valueTag, 1.0);
    }

}

int main(int argc, char **argv)
{
    plbInit(&argc, &argv);

    if (argc!=7) {
        pcout << "Error missing some input parameter\n";
        pcout << "The structure is :\n";
        pcout << "1. Input file name.\n";
        pcout << "2. Output directory name.\n";
        pcout << "3. number of cells in X direction.\n";
        pcout << "4. number of cells in Y direction.\n";
        pcout << "5. number of cells in Z direction.\n";
        pcout << "6. value tag\n";
        pcout << "Example: " << argv[0] << " palabos.dat tmp/ 201 201 201 abc\n";
        exit (EXIT_FAILURE);
    }
    std::string fNameIn  = argv[1];
    std::string fNameOut = argv[2];

    const plint nx = atoi(argv[3]);
    const plint ny = atoi(argv[4]);
    const plint nz = atoi(argv[5]);

    std::string valueTag = argv[6];

    global::directories().setOutputDir(fNameOut+"/");

    pcout << "Reading the geometry file." << std::endl;
    MultiScalarField3D<float> geometry(nx,ny,nz);
    readGeometry(fNameIn, fNameOut, geometry, valueTag);

    return 0;
}
