#include "implicit_point.h"

#pragma intrinsic(fabs)

// Uncomment the following to activate overflow/underflow checks
//#define CHECK_FOR_XYZERFLOWS

int same_side_filtered(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double px, double py, double pz)
{
   double bx_ax = bx - ax;
   double by_ay = by - ay;
   double bz_az = bz - az;
   double cx_ax = cx - ax;
   double cy_ay = cy - ay;
   double cz_az = cz - az;
   double tmi_0 = bx_ax * cy_ay;
   double tmi_1 = by_ay * cx_ax;
   double i01 = tmi_0 - tmi_1;
   double tmj_0 = bx_ax * cz_az;
   double tmj_1 = bz_az * cx_ax;
   double j02 = tmj_0 - tmj_1;
   double tmk_0 = by_ay * cz_az;
   double tmk_1 = bz_az * cy_ay;
   double k12 = tmk_0 - tmk_1;
   double px_ax = px - ax;
   double py_ay = py - ay;
   double pz_az = pz - az;
   double tmii_0 = bx_ax * py_ay;
   double tmii_1 = by_ay * px_ax;
   double ii01 = tmii_0 - tmii_1;
   double tmjj_0 = bx_ax * pz_az;
   double tmjj_1 = bz_az * px_ax;
   double jj02 = tmjj_0 - tmjj_1;
   double tmkk_0 = by_ay * pz_az;
   double tmkk_1 = bz_az * py_ay;
   double kk12 = tmkk_0 - tmkk_1;
   double di = i01 * ii01;
   double dj = j02 * jj02;
   double dk = k12 * kk12;
   double dij = di + dj;
   double dijk = dij + dk;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(bx_ax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(by_ay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bz_az)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cx_ax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cy_ay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cz_az)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(px_ax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(py_ay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pz_az)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.376676550535195e-14;
   if (dijk > epsilon) return IP_Sign::POSITIVE;
   if (-dijk > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int same_side_interval(interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz, interval_number px, interval_number py, interval_number pz)
{
   setFPUModeToRoundUP();
   interval_number bx_ax(bx - ax);
   interval_number by_ay(by - ay);
   interval_number bz_az(bz - az);
   interval_number cx_ax(cx - ax);
   interval_number cy_ay(cy - ay);
   interval_number cz_az(cz - az);
   interval_number tmi_0(bx_ax * cy_ay);
   interval_number tmi_1(by_ay * cx_ax);
   interval_number i01(tmi_0 - tmi_1);
   interval_number tmj_0(bx_ax * cz_az);
   interval_number tmj_1(bz_az * cx_ax);
   interval_number j02(tmj_0 - tmj_1);
   interval_number tmk_0(by_ay * cz_az);
   interval_number tmk_1(bz_az * cy_ay);
   interval_number k12(tmk_0 - tmk_1);
   interval_number px_ax(px - ax);
   interval_number py_ay(py - ay);
   interval_number pz_az(pz - az);
   interval_number tmii_0(bx_ax * py_ay);
   interval_number tmii_1(by_ay * px_ax);
   interval_number ii01(tmii_0 - tmii_1);
   interval_number tmjj_0(bx_ax * pz_az);
   interval_number tmjj_1(bz_az * px_ax);
   interval_number jj02(tmjj_0 - tmjj_1);
   interval_number tmkk_0(by_ay * pz_az);
   interval_number tmkk_1(bz_az * py_ay);
   interval_number kk12(tmkk_0 - tmkk_1);
   interval_number di(i01 * ii01);
   interval_number dj(j02 * jj02);
   interval_number dk(k12 * kk12);
   interval_number dij(di + dj);
   interval_number dijk(dij + dk);
   setFPUModeToRoundNEAR();

   if (!dijk.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return dijk.sign();
}

int same_side_exact(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double px, double py, double pz)
{
   expansionObject o;
   double bx_ax[2];
   o.two_Diff(bx, ax, bx_ax);
   double by_ay[2];
   o.two_Diff(by, ay, by_ay);
   double bz_az[2];
   o.two_Diff(bz, az, bz_az);
   double cx_ax[2];
   o.two_Diff(cx, ax, cx_ax);
   double cy_ay[2];
   o.two_Diff(cy, ay, cy_ay);
   double cz_az[2];
   o.two_Diff(cz, az, cz_az);
   double tmi_0[8];
   int tmi_0_len = o.Gen_Product(2, bx_ax, 2, cy_ay, tmi_0);
   double tmi_1[8];
   int tmi_1_len = o.Gen_Product(2, by_ay, 2, cx_ax, tmi_1);
   double i01[16];
   int i01_len = o.Gen_Diff(tmi_0_len, tmi_0, tmi_1_len, tmi_1, i01);
   double tmj_0[8];
   int tmj_0_len = o.Gen_Product(2, bx_ax, 2, cz_az, tmj_0);
   double tmj_1[8];
   int tmj_1_len = o.Gen_Product(2, bz_az, 2, cx_ax, tmj_1);
   double j02[16];
   int j02_len = o.Gen_Diff(tmj_0_len, tmj_0, tmj_1_len, tmj_1, j02);
   double tmk_0[8];
   int tmk_0_len = o.Gen_Product(2, by_ay, 2, cz_az, tmk_0);
   double tmk_1[8];
   int tmk_1_len = o.Gen_Product(2, bz_az, 2, cy_ay, tmk_1);
   double k12[16];
   int k12_len = o.Gen_Diff(tmk_0_len, tmk_0, tmk_1_len, tmk_1, k12);
   double px_ax[2];
   o.two_Diff(px, ax, px_ax);
   double py_ay[2];
   o.two_Diff(py, ay, py_ay);
   double pz_az[2];
   o.two_Diff(pz, az, pz_az);
   double tmii_0[8];
   int tmii_0_len = o.Gen_Product(2, bx_ax, 2, py_ay, tmii_0);
   double tmii_1[8];
   int tmii_1_len = o.Gen_Product(2, by_ay, 2, px_ax, tmii_1);
   double ii01[16];
   int ii01_len = o.Gen_Diff(tmii_0_len, tmii_0, tmii_1_len, tmii_1, ii01);
   double tmjj_0[8];
   int tmjj_0_len = o.Gen_Product(2, bx_ax, 2, pz_az, tmjj_0);
   double tmjj_1[8];
   int tmjj_1_len = o.Gen_Product(2, bz_az, 2, px_ax, tmjj_1);
   double jj02[16];
   int jj02_len = o.Gen_Diff(tmjj_0_len, tmjj_0, tmjj_1_len, tmjj_1, jj02);
   double tmkk_0[8];
   int tmkk_0_len = o.Gen_Product(2, by_ay, 2, pz_az, tmkk_0);
   double tmkk_1[8];
   int tmkk_1_len = o.Gen_Product(2, bz_az, 2, py_ay, tmkk_1);
   double kk12[16];
   int kk12_len = o.Gen_Diff(tmkk_0_len, tmkk_0, tmkk_1_len, tmkk_1, kk12);
   double di_p[128], *di = di_p;
   int di_len = o.Gen_Product_With_PreAlloc(i01_len, i01, ii01_len, ii01, &di, 128);
   double dj_p[128], *dj = dj_p;
   int dj_len = o.Gen_Product_With_PreAlloc(j02_len, j02, jj02_len, jj02, &dj, 128);
   double dk_p[128], *dk = dk_p;
   int dk_len = o.Gen_Product_With_PreAlloc(k12_len, k12, kk12_len, kk12, &dk, 128);
   double dij_p[128], *dij = dij_p;
   int dij_len = o.Gen_Sum_With_PreAlloc(di_len, di, dj_len, dj, &dij, 128);
   double dijk_p[128], *dijk = dijk_p;
   int dijk_len = o.Gen_Sum_With_PreAlloc(dij_len, dij, dk_len, dk, &dijk, 128);

   double return_value = dijk[dijk_len - 1];
   if (dijk_p != dijk) free(dijk);
   if (dij_p != dij) free(dij);
   if (dk_p != dk) free(dk);
   if (dj_p != dj) free(dj);
   if (di_p != di) free(di);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

int same_side(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double px, double py, double pz)
{
   int ret;
   ret = same_side_filtered(ax, ay, az, bx, by, bz, cx, cy, cz, px, py, pz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = same_side_interval(ax, ay, az, bx, by, bz, cx, cy, cz, px, py, pz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return same_side_exact(ax, ay, az, bx, by, bz, cx, cy, cz, px, py, pz);
}

