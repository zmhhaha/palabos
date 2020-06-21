#include "palabos2D.h"
#include "palabos2D.hh"

#include <vector>
#include <cmath>
#include <cstdlib>

using namespace plb;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

// This function object returns a zero velocity, and a pressure which decreases
//   linearly in z-direction. It is used to initialize the particle populations.
class PressureGradient {
public:
    PressureGradient(T deltaP_, plint ny_) : deltaP(deltaP_), ny(ny_)
    { }
    void operator() (plint iX, plint iY, T& density, Array<T,2>& velocity) const
    {
        velocity.resetToZero();
        density = (T)1 - deltaP*DESCRIPTOR<T>::invCs2 / (T)(ny-1) * (T)(ny-1-iY);
    }
private:
    T deltaP;
    plint ny;
};

void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField2D<int>& geometry)
{
    plb_ifstream geometryFile(fNameIn.c_str());
    if (!geometryFile.is_open()) {
        pcout << "Error: could not open geometry file " << fNameIn << std::endl;
        exit(EXIT_FAILURE);
    }
    geometryFile >> geometry;

    {
        VtkImageOutput2D<T> vtkOut("porousMedium", 1.0);
        vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
    }
}

void porousMediaSetup(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* boundaryCondition,
        MultiScalarField2D<int>& geometry, T deltaP)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    pcout << "Definition of inlet/outlet." << std::endl;
    Box2D inlet (1,nx-2, ny-1,ny-1);
    boundaryCondition->addPressureBoundary1P(inlet, lattice);
    setBoundaryDensity(lattice, inlet, (T) 1.);

    Box2D outlet(1,nx-2, 0,0);
    boundaryCondition->addPressureBoundary1N(outlet, lattice);
    setBoundaryDensity(lattice, outlet, (T) 1. - deltaP*DESCRIPTOR<T>::invCs2);

    pcout << "Definition of the geometry." << std::endl;
    // Where "geometry" evaluates to 1, use bounce-back.
    defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
    // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
    defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

    pcout << "Initilization of rho and u." << std::endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), PressureGradient(deltaP, ny) );

    lattice.initialize();
    delete boundaryCondition;
}

void writeGifs(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif(createFileName("u", iter, 6),
            *computeVelocityNorm(lattice, Box2D(0,nx-1, 0,ny-1)),
            imSize, imSize );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", 1.);
}

T computePermeability(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, T nu, T deltaP, Box2D domain )
{
    pcout << "Computing the permeability." << std::endl;

    // Compute only the y-direction of the velocity (direction of the flow).
    plint yComponent = 1;
    plint ny = lattice.getNy();

    T meanU = -computeAverage(*computeVelocityComponent(lattice, domain, yComponent));

    pcout << "Average velocity     = " << meanU                         << std::endl;
    pcout << "Lattice viscosity nu = " << nu                            << std::endl;
    pcout << "Grad P               = " << deltaP/(T)(ny-1)              << std::endl;
    pcout << "Permeability         = " << nu*meanU / (deltaP/(T)(ny-1)) << std::endl;

    return meanU;
}

int main(int argc, char **argv)
{
    plbInit(&argc, &argv);

    if (argc!=6) {
        pcout << "Error missing some input parameter\n";
        pcout << "The structure is :\n";
        pcout << "1. Input file name.\n";
        pcout << "2. Output directory name.\n";
        pcout << "3. number of cells in X direction.\n";
        pcout << "4. number of cells in Y direction.\n";
        pcout << "5. Delta P .\n";
        pcout << "Example: " << argv[0] << " palabos.dat tmp/ 201 201 0.00005\n";
        exit (EXIT_FAILURE);
    }
    std::string fNameIn  = argv[1];
    std::string fNameOut = argv[2];

    const plint nx = atoi(argv[3]);
    const plint ny = atoi(argv[4]);
    const T deltaP = atof(argv[5]);

    global::directories().setOutputDir(fNameOut+"/");

    const T omega = 1.0;
    const T nu    = ((T)1/omega- (T)0.5)/DESCRIPTOR<T>::invCs2;

    pcout << "Creation of the lattice." << std::endl;
    MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx,ny, new BGKdynamics<T,DESCRIPTOR>(omega));
    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    pcout << "Reading the geometry file." << std::endl;
    MultiScalarField2D<int> geometry(nx,ny);
    readGeometry(fNameIn, fNameOut, geometry);

    pcout << "nu = " << nu << std::endl;
    pcout << "deltaP = " << deltaP << std::endl;
    pcout << "omega = " << omega << std::endl;
    pcout << "nx = " << lattice.getNx() << std::endl;
    pcout << "ny = " << lattice.getNy() << std::endl;

    porousMediaSetup(lattice, createLocalBoundaryCondition2D<T,DESCRIPTOR>(), geometry, deltaP);

    // The value-tracer is used to stop the simulation once is has converged.
    // 1st parameter:velocity
    // 2nd parameter:size
    // 3rd parameters:threshold
    // 1st and second parameters ae used for the length of the time average (size/velocity)
    util::ValueTracer<T> converge(1.0,1000.0,1.0e-4);

    pcout << "Simulation begins" << std::endl;
    plint iT=0;

    const plint maxT = 100000;
    for (;iT<maxT; ++iT) {
        if (iT % 20 == 0) {
            pcout << "Iteration " << iT << std::endl;
        }
        if (iT % 50 == 0 && iT>0) {
            writeGifs(lattice,iT);
        }

        lattice.collideAndStream();
        converge.takeValue(getStoredAverageEnergy(lattice),true);

        if (converge.hasConverged()) {
            break;
        }
    }

    pcout << "End of simulation at iteration " << iT << std::endl;

    pcout << "Permeability:" << std::endl << std::endl;
    computePermeability(lattice, nu, deltaP, lattice.getBoundingBox());
    pcout << std::endl;

    pcout << "Writing VTK file ..." << std::endl << std::endl;
    writeVTK(lattice,iT);
    pcout << "Finished!" << std::endl << std::endl;

    return 0;
}
————————————————
版权声明：本文为CSDN博主「胖虎一号」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/qq_28632981/article/details/106643221
