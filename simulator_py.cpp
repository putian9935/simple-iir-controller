#include "simulator_pc.cpp"
#include <cstdlib>
#include <cstring>

constexpr size_t tot_size = sizeof readings / sizeof(double);
static double errors[tot_size];

template <typename Simulator> void simulate(Simulator sim) {
  reading_index = 0;
  for (size_t i = 0; i < tot_size; ++i) {
    sim.update();
    errors[i] = sim.last_error;
  }
};

double output2input(double const); // forward declaration of plant's DC transfer

template <int len_zeroes, int len_poles>
void simulate_c_interf(double *zeroes, double *poles, double overall_gain) {
  static double local_zeroes[len_zeroes], local_poles[len_poles];
  memcpy(local_zeroes, zeroes, len_zeroes * sizeof(double));
  memcpy(local_poles, poles, len_poles * sizeof(double));
  simulate(make_iir_cascade_simulator(local_zeroes, local_poles, overall_gain,
                                      output2input));
}

/******************************/

constexpr int lz = 1; // size of zeroes
constexpr int lp = 2; // size of poles

// the DC transfer function of the plant
double output2input(double const x) { return 4.1e-4 * x * 6 / 65536; }

/******************************/

extern "C" {

__declspec(dllexport) int get_len_zeroes() { return lz; }
__declspec(dllexport) int get_len_poles() { return lp; }

__declspec(dllexport) double *simulate(double *zeroes, double *poles,
                                       double overall_gain) {
  simulate_c_interf<lz, lp>(zeroes, poles, overall_gain);
  return errors;
}
}
