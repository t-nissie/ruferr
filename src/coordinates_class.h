// coordinates_class.h
////
#pragma once
#include <chrono>
#include<fftw3.h>
#include<boost/multi_array.hpp>
typedef boost::multi_array<double,4> double_array4_t;

class Coordinates
{
private:
  int _padding_y, _howmany, _L[4], _d[4], _z[4];
  double _invN;
  double_array4_t _dipole;
  fftw_plan _plan_c2r_3_in = NULL;
  fftw_plan _plan_r2c_3_in = NULL;

  template<typename _Rep, typename _Period>
  static void report_timing(const std::chrono::duration<_Rep, _Period>& dur,
                            const char *s)
  {
    std::cout << std::setw(16) << std::left << s << ": "
              << std::setw( 5) << std::right
              << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()
              << " ms" << std::endl;
  }

public:
  Coordinates(const int howmany, const int Lx, const int Ly, const int Lz, const int padding_y);
  ~Coordinates();
  void forward_backward(const int n_times);
  void set_values();
  void show_some();
};
