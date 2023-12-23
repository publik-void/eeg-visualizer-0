#include <string>
#include <sstream>

template <typename T>
std::string fancy(T const x) {
  std::stringstream s;
  s << x;
  return s.str();
}

std::string fancy(float const x) { return fancy<float>(x); }
std::string fancy(double const x) { return fancy<double>(x); }
