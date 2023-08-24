#include "BezierTriangleInjectivityChecker.h"

namespace GCLF
{
namespace Geometry
{
const double BezierTriangleInjectivityChecker::fac_[] =
//0    1    2    3    4     5      6      7       8       9        10
{ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320., 362880., 3628800. };

void BezierTriangleInjectivityChecker::set_degree(int degree)
{
  degreeB = degree;
  degreeC = 2 * (degree - 1);

  indB = std::vector<int>(degreeB + 1);
  for (int i = 0, val = 0, add = degreeB + 1; i <= degreeB; i++)
  {
    indB[i] = val;
    val += add;
    add--;
  }
  indC = std::vector<int>(degreeC + 1);
  for (int i = 0, val = 0, add = degreeC + 1; i <= degreeC; i++)
  {
    indC[i] = val;
    val += add;
    add--;
  }
}

double BezierTriangleInjectivityChecker::P(
  const std::vector<double>& C, int i, int j, int n1, int n2, int n3)
{
  if (n1 > 0)
    return 0.5 * (P(C, i, j, n1 - 1, n2, n3) + P(C, i - 1, j + 1, n1 - 1, n2, n3));
  if (n2 > 0)
    return 0.5 * (P(C, i, j, n1, n2 - 1, n3) + P(C, i, j - 1, n1, n2 - 1, n3));
  if (n3 > 0)
    return 0.5 * (P(C, i, j, n1, n2, n3 - 1) + P(C, i + 1, j, n1, n2, n3 - 1));
  return C[indC[i] + j];
}


void BezierTriangleInjectivityChecker::bisect(
  const std::vector<double>& C,
  int degree,
  const std::vector<int>& ind,
  std::vector<double>& C0,
  std::vector<double>& C1)
{
  std::vector<double> dCrow(degree + 1);
  for (int i = 0; i <= degree; i++)
  {
    for (int j = 0; j <= degree - i; j++)
      dCrow[j] = C[ind[i] + j];
    C0[ind[0] + i] = dCrow[degree - i];
    C1[ind[0] + degree - i] = dCrow[0];
    for (int l = 1; l <= degree - i; l++)
    {
      for (int j = 0; j <= degree - i - l; j++)
        dCrow[j] = 0.5 * (dCrow[j] + dCrow[j + 1]);
      C0[ind[l] + i] = dCrow[degree - i - l];
      C1[ind[l] + degree - i - l] = dCrow[0];
    }
  }
}


/// @brief Check if a plannar bezier triangle is definitely injective by determinant.
/// @param _b1 the x component of control points.
/// @param _b2 the y component of control points.
/// @param via_bisect use 1-2 subdivision if true or 1-4 subdivision if false.
/// @return integer indicating injectivity.
/// @retval 0 : determinant non-positive at some (corner) point.
/// @retval 1 : definitely injective.
/// @retval 2 : injectivity could not be certified.
/// @note control points are expected in lexicographic index order,
/// i.e b003, b012, b021, b030, b102, b111, b120, b201, b210, b300.
int BezierTriangleInjectivityChecker::is_definitely_injective_by_determinant(
  const std::vector<double>& _b1, const std::vector<double>& _b2, bool via_bisect)
{
  v1 = _b1;
  v2 = _b2;

  // compute quartic control values
  int derivDegree = degreeB - 1;
  std::vector<double> C(((degreeC + 1) * (degreeC + 2)) / 2, 0.0);
  for (int i = 0; i <= degreeC; i++)
    for (int j = 0; j <= degreeC - i; j++)
    {
      for (int ri = 0; ri <= derivDegree; ri++)
        for (int rj = 0; rj <= derivDegree - ri; rj++)
        {
          int si = i - ri;
          int sj = j - rj;
          int rk = derivDegree - ri - rj;
          int sk = derivDegree - si - sj;
          if (si < 0)
            continue;
          if (sj < 0)
            continue;
          if (si + sj > derivDegree)
            continue;
          if (degreeC - ri - rj - si - sj != degreeC - i - j)
            continue;
          double lijk = (facT(i) * facT(j) * facT(degreeC - i - j))
            / (facT(ri) * facT(rj) * facT(rk) * facT(si) * facT(sj)
              * facT(sk)); // constant factor does not matter:  * 1.0/6.0;
          // (general: 1.0/facT(degreeB))
          C[indC[i] + j] += lijk * (xu(ri, rj) * yv(si, sj) - xv(ri, rj) * yu(si, sj));
        }
    }

  std::queue<std::pair<std::vector<double>, int>> q;
  q.push(std::make_pair(C, 0));
  while (!q.empty())
  {
    std::vector<double> C = q.front().first;
    current_level = q.front().second;
    q.pop();

    double mi = *std::min_element(C.begin(), C.end());
    if (mi > min_area)
      continue; // determinant positive everywhere in this triangle

    double mic
      = std::min(C[indC[degreeC] + 0],
        std::min(C[indC[0] + degreeC], C[indC[0] + 0])); // minimum of corner values
    if (mic <= min_area)
    {
      return 0; // determinant non-positive at some (corner) point
    }

    if (current_level >= max_subdiv)
    {
      return 2; // injectivity could not be certified 
    }

    if (via_bisect)
    {
      // two subtriangles
      std::vector<double> C0(C.size());
      std::vector<double> C1(C.size());
      bisect(C, degreeC, indC, C0, C1);
      q.push(std::make_pair(C0, current_level + 1));
      q.push(std::make_pair(C1, current_level + 1));
    }
    else
    {
      // four subtriangles
      std::vector<double> C0, C1, C2, C3;
      for (int i = 0; i <= degreeC; i++)
        for (int m = 0; m <= i; m++)
          C0.push_back(P(C, i, 0, m, 0, degreeC - i));
      for (int i = 0; i <= degreeC; i++)
        for (int m = 0; m <= degreeC - i; m++)
          C1.push_back(P(C, i, degreeC - i, i, m, 0));
      for (int i = 0; i <= degreeC; i++)
        for (int m = 0; m <= degreeC - i; m++)
          C2.push_back(P(C, 0, i, 0, i, m));
      for (int i = 0; i <= degreeC; i++)
        for (int j = 0; j <= degreeC - i; j++)
          C3.push_back(P(C, i, j, i, j, degreeC - i - j));

      q.push(std::make_pair(C0, current_level + 1));
      q.push(std::make_pair(C1, current_level + 1));
      q.push(std::make_pair(C2, current_level + 1));
      q.push(std::make_pair(C3, current_level + 1));
    }
  }
  return 1;
}

/// @brief A simple wrap. Check if a plannar bezier triangle is definitely injective by determinant.
/// @param ctrlpnts Control points of a plannar Bezier triangle.
/// They come from gclf::BezierTriangle.
/// @param via_bisect use 1-2 subdivision if true or 1-4 subdivision if false.
/// @return integer indicating injectivity.
/// @retval 0 : determinant non-positive at some (corner) point.
/// @retval 1 : definitely injective.
/// @retval 2 : injectivity could not be certified.
int BezierTriangleInjectivityChecker::is_definitely_injective_by_determinant(
  const Points ctrlpnts,
  bool via_bisect)
{
  // Adjust the order of control points.
  // b1 and b2 are the mirror reversed order of ctrlpnts.
  // An example of ctrlpnts:
  // b300(0),
  // b201(1), b210(2),
  // b102(3), b111(4), b120(5),
  // b003(6), b012(7), b021(8), b030(9)
  // An example of b1 and b2:
  // b300(9),
  // b201(7), b210(8),
  // b102(4), b111(5), b120(6),
  // b003(0), b012(1), b021(2), b030(3)
  // We need to correct the reversed order.
  // The mirror transformation is not important.
  std::vector<double> b1, b2;
  b1.reserve(ctrlpnts.size()); b2.reserve(ctrlpnts.size());
  for (auto iter = ctrlpnts.rbegin(); iter != ctrlpnts.rend(); iter++)
  {
    b1.push_back(iter->x());
    b2.push_back(iter->y());
  }

  return is_definitely_injective_by_determinant(b1, b2, via_bisect);
}

}// namespace Geometry
}// namespace GCLF