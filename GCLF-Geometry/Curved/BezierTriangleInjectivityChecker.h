#pragma once
#include <algorithm>
#include <queue>
#include <vector>
#include "BezierTriangle.h"

namespace GCLF
{
namespace Geometry
{
/**
 * @brief Utility class to check validity of bezier triangles.
 *
 */

class BezierTriangleInjectivityChecker
{
private:
  double min_area = 0.0;
#ifdef ALL_EXACT
  int max_subdiv = 200; // number of subdivisions for conservative injectivity test
#else
  int max_subdiv = 15;  // more does not make much sense in double precision
#endif

  int current_level = 0;
  std::vector<double> v1, v2;
  int degreeB;
  int degreeC;

  //(degree+1)*i - ((i*(i-1))/2) + j;
  std::vector<int> indB;
  std::vector<int> indC;

  double b1(int i, int j) { return v1[indB[i] + j]; }
  double b2(int i, int j) { return v2[indB[i] + j]; }

  // For below four functions, constant factor does not matter: *3.0
  double xu(int i, int j) { return (b1(i + 1, j) - b1(i, j)); }
  double xv(int i, int j) { return (b1(i, j + 1) - b1(i, j)); }
  double yu(int i, int j) { return (b2(i + 1, j) - b2(i, j)); }
  double yv(int i, int j) { return (b2(i, j + 1) - b2(i, j)); }

  double P(
    const std::vector<double>& C,
    int i, int j, int n1, int n2, int n3
  );

  void bisect(
    const std::vector<double>& C,
    int degree,
    const std::vector<int>& ind,
    std::vector<double>& C0,
    std::vector<double>& C1
  );

  // TODO: dynamically generate fac_ depending on degree.
  static const double fac_[];

  double facT(int n)
  {
    // if (n == 0) return 1.0;
    // if (n == 1) return 1.0;
    // if (n == 2) return 2.0;
    // if (n == 3) return 6.0;
    // if (n == 4) return 24.0;
    // if (n == 5) return 120.0;
    // return fac(n - 1) * n;
    // NOTE: in most cases, subscript n won't overflow. 
    return fac_[n];
  }

public:
  BezierTriangleInjectivityChecker() = default;

  BezierTriangleInjectivityChecker(int _degree) { set_degree(_degree); }

  void set_degree(int degree);
  void set_max_subdiv(int subdiv) { max_subdiv = subdiv; }

  int is_definitely_injective_by_determinant(
    const std::vector<double>& _b1,
    const std::vector<double>& _b2,
    bool via_bisect = true
  );

  int is_definitely_injective_by_determinant(
    const Points ctrlpnts,
    bool via_bisect = true
  );

  int level() { return current_level; }
};

}// namespace Geometry
}// namespace GCLF