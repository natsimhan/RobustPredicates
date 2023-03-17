<?php declare(strict_types=1);

namespace Natsimhan\RobustPredicates;

/**
 * Note: unlike J. Shewchuk's original code, all the functions in this library
 * assume y axis is oriented downwards ↓, so the semantics are different.
 *
 * @author        : original : Volodymyr Agafonkin @mourner
 * @author        : PHP's adaptation : Jonathan BURON aka Natsimhan <jonathan@ouebsson.fr>
 * @source        : https://github.com/mourner/robust-predicates/blob/master/src/orient2d.js
 */
final class Orient2d {

  private const CCW_ERR_BOUND_A = (3 + 16 * Utils::EPSILON) * Utils::EPSILON;
  private const CCW_ERR_BOUND_B = (2 + 12 * Utils::EPSILON) * Utils::EPSILON;
  private const CCW_ERR_BOUND_C = (9 + 64 * Utils::EPSILON) * Utils::EPSILON * Utils::EPSILON;

  /**
   * Returns a positive value if the points a, b, and c occur in
   * counterclockwise order (c lies to the left of the directed line defined by
   * points a and b).
   *
   * Returns a negative value if they occur in clockwise order (c lies to the
   * right of the directed line ab).
   *
   * Returns zero if they are collinear.
   *
   * The result is also an approximation of twice the signed area of the
   * triangle defined by the three points.
   *
   * @param float $ax
   * @param float $ay
   * @param float $bx
   * @param float $by
   * @param float $cx
   * @param float $cy
   * @return float
   */
  public static function orient2d(
    float $ax, float $ay,
    float $bx, float $by,
    float $cx, float $cy,
  ): float {
    $detleft = ($ay - $cy) * ($bx - $cx);
    $detright = ($ax - $cx) * ($by - $cy);
    $det = $detleft - $detright;

    if($detleft === 0.0 || $detright === 0.0 || ($detleft > 0.0) !== ($detright > 0.0)) {
      return $det;
    }

    $detsum = abs($detleft + $detright);
    if(abs($det) >= self::CCW_ERR_BOUND_A * $detsum) {
      return $det;
    }

    return -self::orient2dadapt($ax, $ay, $bx, $by, $cx, $cy, $detsum);
  }

  /**
   * Simple, approximate, non-robust versions of above predicates orient2d().
   * Use when robustness isn't needed.
   *
   * @param float $ax
   * @param float $ay
   * @param float $bx
   * @param float $by
   * @param float $cx
   * @param float $cy
   * @return float
   */
  public static function orient2dfast(
    float $ax, float $ay,
    float $bx, float $by,
    float $cx, float $cy,
  ): float {
    return ($ay - $cy) * ($bx - $cx) - ($ax - $cx) * ($by - $cy);
  }

  private static function orient2dadapt(float $ax, float $ay, float $bx, float $by, float $cx, float $cy, float $detsum): float {
    $B = Utils::vec(4);
    $C1 = Utils::vec(8);
    $C2 = Utils::vec(12);
    $D = Utils::vec(16);
    $u = Utils::vec(4);

    $acx = $ax - $cx;
    $bcx = $bx - $cx;
    $acy = $ay - $cy;
    $bcy = $by - $cy;

    $s1 = $acx * $bcy;
    $c = Utils::SPLITTER * $acx;
    $ahi = $c - ($c - $acx);
    $alo = $acx - $ahi;
    $c = Utils::SPLITTER * $bcy;
    $bhi = $c - ($c - $bcy);
    $blo = $bcy - $bhi;
    $s0 = $alo * $blo - ($s1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $t1 = $acy * $bcx;
    $c = Utils::SPLITTER * $acy;
    $ahi = $c - ($c - $acy);
    $alo = $acy - $ahi;
    $c = Utils::SPLITTER * $bcx;
    $bhi = $c - ($c - $bcx);
    $blo = $bcx - $bhi;
    $t0 = $alo * $blo - ($t1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $_i = $s0 - $t0;
    $bvirt = $s0 - $_i;
    $B[0] = $s0 - ($_i + $bvirt) + ($bvirt - $t0);
    $_j = $s1 + $_i;
    $bvirt = $_j - $s1;
    $_0 = $s1 - ($_j - $bvirt) + ($_i - $bvirt);
    $_i = $_0 - $t1;
    $bvirt = $_0 - $_i;
    $B[1] = $_0 - ($_i + $bvirt) + ($bvirt - $t1);
    $u3 = $_j + $_i;
    $bvirt = $u3 - $_j;
    $B[2] = $_j - ($u3 - $bvirt) + ($_i - $bvirt);
    $B[3] = $u3;

    $det = Utils::estimate(4, $B);
    $errbound = self::CCW_ERR_BOUND_B * $detsum;
    if($det >= $errbound || -$det >= $errbound) {
      return $det;
    }

    $bvirt = $ax - $acx;
    $acxtail = $ax - ($acx + $bvirt) + ($bvirt - $cx);
    $bvirt = $bx - $bcx;
    $bcxtail = $bx - ($bcx + $bvirt) + ($bvirt - $cx);
    $bvirt = $ay - $acy;
    $acytail = $ay - ($acy + $bvirt) + ($bvirt - $cy);
    $bvirt = $by - $bcy;
    $bcytail = $by - ($bcy + $bvirt) + ($bvirt - $cy);

    if($acxtail === 0.0 && $acytail === 0.0 && $bcxtail === 0.0 && $bcytail === 0.0) {
      return $det;
    }

    $errbound = self::CCW_ERR_BOUND_C * $detsum + Utils::RESULT_ERR_BOUND * abs($det);
    $det += ($acx * $bcytail + $bcy * $acxtail) - ($acy * $bcxtail + $bcx * $acytail);
    if($det >= $errbound || -$det >= $errbound) {
      return $det;
    }

    $s1 = $acxtail * $bcy;
    $c = Utils::SPLITTER * $acxtail;
    $ahi = $c - ($c - $acxtail);
    $alo = $acxtail - $ahi;
    $c = Utils::SPLITTER * $bcy;
    $bhi = $c - ($c - $bcy);
    $blo = $bcy - $bhi;
    $s0 = $alo * $blo - ($s1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $t1 = $acytail * $bcx;
    $c = Utils::SPLITTER * $acytail;
    $ahi = $c - ($c - $acytail);
    $alo = $acytail - $ahi;
    $c = Utils::SPLITTER * $bcx;
    $bhi = $c - ($c - $bcx);
    $blo = $bcx - $bhi;
    $t0 = $alo * $blo - ($t1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $_i = $s0 - $t0;
    $bvirt = $s0 - $_i;
    $u[0] = $s0 - ($_i + $bvirt) + ($bvirt - $t0);
    $_j = $s1 + $_i;
    $bvirt = $_j - $s1;
    $_0 = $s1 - ($_j - $bvirt) + ($_i - $bvirt);
    $_i = $_0 - $t1;
    $bvirt = $_0 - $_i;
    $u[1] = $_0 - ($_i + $bvirt) + ($bvirt - $t1);
    $u3 = $_j + $_i;
    $bvirt = $u3 - $_j;
    $u[2] = $_j - ($u3 - $bvirt) + ($_i - $bvirt);
    $u[3] = $u3;
    $C1len = Utils::sum(4, $B, 4, $u, $C1);

    $s1 = $acx * $bcytail;
    $c = Utils::SPLITTER * $acx;
    $ahi = $c - ($c - $acx);
    $alo = $acx - $ahi;
    $c = Utils::SPLITTER * $bcytail;
    $bhi = $c - ($c - $bcytail);
    $blo = $bcytail - $bhi;
    $s0 = $alo * $blo - ($s1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $t1 = $acy * $bcxtail;
    $c = Utils::SPLITTER * $acy;
    $ahi = $c - ($c - $acy);
    $alo = $acy - $ahi;
    $c = Utils::SPLITTER * $bcxtail;
    $bhi = $c - ($c - $bcxtail);
    $blo = $bcxtail - $bhi;
    $t0 = $alo * $blo - ($t1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $_i = $s0 - $t0;
    $bvirt = $s0 - $_i;
    $u[0] = $s0 - ($_i + $bvirt) + ($bvirt - $t0);
    $_j = $s1 + $_i;
    $bvirt = $_j - $s1;
    $_0 = $s1 - ($_j - $bvirt) + ($_i - $bvirt);
    $_i = $_0 - $t1;
    $bvirt = $_0 - $_i;
    $u[1] = $_0 - ($_i + $bvirt) + ($bvirt - $t1);
    $u3 = $_j + $_i;
    $bvirt = $u3 - $_j;
    $u[2] = $_j - ($u3 - $bvirt) + ($_i - $bvirt);
    $u[3] = $u3;
    $C2len = Utils::sum($C1len, $C1, 4, $u, $C2);

    $s1 = $acxtail * $bcytail;
    $c = Utils::SPLITTER * $acxtail;
    $ahi = $c - ($c - $acxtail);
    $alo = $acxtail - $ahi;
    $c = Utils::SPLITTER * $bcytail;
    $bhi = $c - ($c - $bcytail);
    $blo = $bcytail - $bhi;
    $s0 = $alo * $blo - ($s1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $t1 = $acytail * $bcxtail;
    $c = Utils::SPLITTER * $acytail;
    $ahi = $c - ($c - $acytail);
    $alo = $acytail - $ahi;
    $c = Utils::SPLITTER * $bcxtail;
    $bhi = $c - ($c - $bcxtail);
    $blo = $bcxtail - $bhi;
    $t0 = $alo * $blo - ($t1 - $ahi * $bhi - $alo * $bhi - $ahi * $blo);
    $_i = $s0 - $t0;
    $bvirt = $s0 - $_i;
    $u[0] = $s0 - ($_i + $bvirt) + ($bvirt - $t0);
    $_j = $s1 + $_i;
    $bvirt = $_j - $s1;
    $_0 = $s1 - ($_j - $bvirt) + ($_i - $bvirt);
    $_i = $_0 - $t1;
    $bvirt = $_0 - $_i;
    $u[1] = $_0 - ($_i + $bvirt) + ($bvirt - $t1);
    $u3 = $_j + $_i;
    $bvirt = $u3 - $_j;
    $u[2] = $_j - ($u3 - $bvirt) + ($_i - $bvirt);
    $u[3] = $u3;
    $Dlen = Utils::sum($C2len, $C2, 4, $u, $D);

    return $D[$Dlen - 1];
  }

}
