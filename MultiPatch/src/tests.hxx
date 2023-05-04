#ifndef MULTIPATCH_TESTS_HXX
#define MULTIPATCH_TESTS_HXX

#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <sstream>

namespace MultiPatchTests {

/**
 * Seed value of the random number generation engines.
 */
constexpr const int random_seed = 100;

/**
 * The "grid spacing" used in finite difference operators
 */
constexpr const CCTK_REAL fd_delta = 1.0e-3;

/**
 * The floating point comparison tolerance when testing the equality of exact
 * and FD computed derivatives
 */
constexpr const CCTK_REAL fd_comp_tol = 1.0e-7; // sqrt(machine eps)

/**
 * Test two floating point values for approximate equality.
 *
 * This function is based on
 * https://docs.julialang.org/en/v1/base/math/#Base.isapprox
 *
 * @param x The first number to compare.
 * @param y The second number to compare.
 * @param atol The absolute tolerance for the comparison.
 * @param rtol The relative tolerance of the comparison.
 * @return True if x ~ y, false otherwise.
 */
template <typename fp_type>
CCTK_DEVICE CCTK_HOST inline bool isapprox(fp_type x, fp_type y,
                                           fp_type atol = 0.0) {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol =
      atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon());
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

/**
 * Test if -boundary < variable < boundary. Fails if variable ~= +/- boundary.
 *
 * @param variable The variable to test.
 * @param boundary The absolute value of the region boundary
 * @return A boolean indicating if the variable is within the region.
 */
template <typename T>
CCTK_DEVICE CCTK_HOST inline bool within(T variable, T boundary) {
  return (variable > -boundary && !isapprox(variable, -boundary)) &&
         (variable < boundary && !isapprox(variable, boundary));
}

/**
 * Test if -boundary ~= variable or variable ~= boundary.
 *
 * @param variable The variable to test.
 * @param boundary The absolute value of the region boundary
 * @return A boolean indicating if the variable is within the region.
 */
template <typename T>
CCTK_DEVICE CCTK_HOST inline bool at_boundary(T variable, T boundary) {
  return (isapprox(variable, boundary) || isapprox(variable, -boundary));
}

/**
 * Tag representing possible colors to apply to strings.
 */
enum class string_color { green, red };

/**
 * Formats a string to be colored in ANSI compatible terminals.
 *
 * @tparam color The color to use.
 * @param string The string to color.
 * @return The colored string.
 */
template <string_color color> std::string colored(const std ::string &str) {
  std::ostringstream output;

  if constexpr (color == string_color::red) {
    output << "\033[31;1m";
    output << str;
    output << "\033[0m";
  } else if constexpr (color == string_color::green) {
    output << "\033[32;1m";
    output << str;
    output << "\033[0m";
  }

  return output.str();
}

const auto PASSED = colored<string_color::green>("PASSED");
const auto FAILED = colored<string_color::red>("FAILED");

/**
 * Determines the direction that a finite difference derivative will be
 * performed
 */
enum class fd_direction { x = 0, y = 1, z = 2 };

/**
 * Computes the second order accurate finite difference derivative of a
 * function that takes a vector as input and produces another vector as output
 * in a specified direction.
 *
 * @param function The function to derivate
 * @param point The point where the derivative is to be computed.
 * @tparam dir The direction of the derivative.
 * @return The derivative of function in the specified direction and point.
 */
template <fd_direction dir, typename vector_t, typename function_t>
inline vector_t fd_4(const function_t &function, vector_t point) {

  vector_t point_p_1d = point;
  vector_t point_p_2d = point;
  vector_t point_m_1d = point;
  vector_t point_m_2d = point;

  point_p_1d(static_cast<int>(dir)) += fd_delta;
  point_p_2d(static_cast<int>(dir)) += 2 * fd_delta;
  point_m_1d(static_cast<int>(dir)) -= fd_delta;
  point_m_2d(static_cast<int>(dir)) -= 2 * fd_delta;

  const auto f_p_1d = function(point_p_1d);
  const auto f_p_2d = function(point_p_2d);
  const auto f_m_1d = function(point_m_1d);
  const auto f_m_2d = function(point_m_2d);
  return (f_m_2d - 8 * f_m_1d + 8 * f_p_1d - f_p_2d) / (12 * fd_delta);
}

template <fd_direction dir_inner, fd_direction dir_outer, typename vector_t,
          typename function_t>
inline vector_t fd2_4(const function_t &function, vector_t point) {

  vector_t point_p_1d = point;
  vector_t point_p_2d = point;
  vector_t point_m_1d = point;
  vector_t point_m_2d = point;

  point_p_1d(static_cast<int>(dir_outer)) += fd_delta;
  point_p_2d(static_cast<int>(dir_outer)) += 2 * fd_delta;
  point_m_1d(static_cast<int>(dir_outer)) -= fd_delta;
  point_m_2d(static_cast<int>(dir_outer)) -= 2 * fd_delta;

  const auto f_p_1d = fd_4<dir_inner>(function, point_p_1d);
  const auto f_p_2d = fd_4<dir_inner>(function, point_p_2d);
  const auto f_m_1d = fd_4<dir_inner>(function, point_m_1d);
  const auto f_m_2d = fd_4<dir_inner>(function, point_m_2d);

  return (f_m_2d - 8 * f_m_1d + 8 * f_p_1d - f_p_2d) / (12 * fd_delta);
}

} // namespace MultiPatchTests

#endif // MULTIPATCH_TESTS_HXX
