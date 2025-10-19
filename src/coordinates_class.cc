// coordinates_class.cc
////
#include<iostream>
#include<omp.h>
#include<fftw3.h>
#include<boost/multi_array.hpp>
#include"coordinates_class.h"
Coordinates::Coordinates(const int howmany, const int Lx, const int Ly, const int Lz, const int padding_y)
  : _padding_y(padding_y)
  , _howmany(howmany)
  , _L{Lx, Ly,              Lz,      Lx*Ly*Lz}
  , _d{Lx, Ly+padding_y, 2*(Lz/2+1), 1}
  , _z{Lx, Ly+padding_y,    Lz/2+1,  1}
  , _invN(1.0/_L[3])
  , _dipole(boost::extents[howmany][Lx][Ly+padding_y][ 2*(Lz/2+1) ]) // SoA!
{
  for (int i=0; i<3; ++i) _d[3] *= _d[i];
  for (int i=0; i<3; ++i) _z[3] *= _z[i];

  if (fftw_init_threads()) {
    fftw_plan_with_nthreads(omp_get_max_threads());
  } else {
    throw std::runtime_error("Could not fftw_init_threads().");
  }
  _plan_r2c_3_in = fftw_plan_many_dft_r2c(3, _L, _howmany, // rank, n[], howmany
                                          _dipole.data(),  _d, 1 , _d[3],
         reinterpret_cast<fftw_complex *>(_dipole.data()), _z, 1 , _z[3], FFTW_MEASURE);
  _plan_c2r_3_in = fftw_plan_many_dft_c2r(3, _L, _howmany,
         reinterpret_cast<fftw_complex *>(_dipole.data()), _z, 1 , _z[3],
                                          _dipole.data(),  _d, 1 , _d[3], FFTW_MEASURE);
}


Coordinates::~Coordinates()
{
  fftw_destroy_plan(_plan_r2c_3_in);
  fftw_destroy_plan(_plan_c2r_3_in);
  std::cout << "An instance of Coordinates class is destructed." << std::endl << std::flush;
}

void Coordinates::forward_backward(const int n_times)
{
  auto start = std::chrono::system_clock::now();
  for (int it=0; it<n_times; ++it) {
    fftw_execute_dft_r2c(_plan_r2c_3_in,
                         _dipole.data(),
                         reinterpret_cast<fftw_complex *>(_dipole.data()));
    fftw_execute_dft_c2r(_plan_c2r_3_in,
                         reinterpret_cast<fftw_complex *>(_dipole.data()),
                         _dipole.data());
#   pragma omp parallel for
    for (int i=0; i<3*_d[3]; ++i) {
      _dipole.data()[i] *= _invN;
    }
  }
  auto end = std::chrono::system_clock::now();
  report_timing(end-start, "_plan_r2c_3_in and _plan_c2r_3_in");
}

void Coordinates::set_values()
{
  for (int id=0; id<3; ++id) {
# pragma omp parallel for
    for (int ix=0; ix<_L[0]; ++ix) {
      for (int iy=0; iy<_L[1]; ++iy) {
        for (int iz=0; iz<_L[2]; ++iz) {
          _dipole[id][ix][iy][iz] = id + 0.0001*ix +  0.00000001*iy +  0.000000000001*iz;
        }
      }
    }
  }
}

void Coordinates::show_some()
{
  for (int iz=0; iz<5; ++iz) {
    printf("%14.12f\n", _dipole[1][2][3][iz]);
  }
}
