#include <iostream>
#include <vector>
extern "C" {
__declspec(dllexport) int
#ifdef WIN32
    __stdcall
#else
    __cdecl
#endif
    calc(double *temp, const size_t length) {
  double timeStep = 0.1;
  double currentTime = 0.;
  double totalTime = 50.;

  double totalTemp = 0.;
  std::vector<double> tempChange(length, 0.);

  temp[20] = 100.;

  for (currentTime = 0.; currentTime < totalTime; currentTime += timeStep) {
    for (std::size_t i = 1; i + 1 < length; ++i) {
      tempChange[i] = 0.1 * (temp[i + 1] - 2. * temp[i] + temp[i - 1]);
    }

    totalTemp = 0.;

    for (std::size_t i = 0; i < length; ++i) {
      temp[i] += tempChange[i] * timeStep;

      totalTemp += temp[i];
    }
  }

  return 0;
}
}