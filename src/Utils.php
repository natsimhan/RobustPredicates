<?php declare(strict_types=1);

namespace Natsimhan\RobustPredicates;

/**
 * @author        : original : Volodymyr Agafonkin @mourner
 * @author        : PHP's adaptation : Jonathan BURON aka Natsimhan <jonathan@ouebsson.fr>
 * @source https://github.com/mourner/robust-predicates/blob/master/src/util.js
 */
final class Utils {

  public const EPSILON = 1.1102230246251565e-16;
  public const SPLITTER = 134217729;
  public const RESULT_ERR_BOUND = (3 + 8 * self::EPSILON) * self::EPSILON;

  //

  /**
   * fast_expansion_sum_zeroelim routine from original code
   *
   * @param int $elen
   * @param array<float> $e
   * @param int $flen
   * @param array<float> $f
   * @param array<float> $h
   * @return int
   */
  public static function sum(int $elen, array $e, int $flen, array $f, array &$h): int {
    $enow = $e[0];
    $fnow = $f[0];
    $eindex = 0;
    $findex = 0;
    if(($fnow > $enow) === ($fnow > -$enow)) {
      $Q = $enow;
      $enow = $e[++$eindex];
    } else {
      $Q = $fnow;
      $fnow = $f[++$findex];
    }
    $hindex = 0;
    if($eindex < $elen && $findex < $flen) {
      if(($fnow > $enow) === ($fnow > -$enow)) {
        $Qnew = $enow + $Q;
        $hh = $Q - ($Qnew - $enow);
        $enow = $e[++$eindex];
      } else {
        $Qnew = $fnow + $Q;
        $hh = $Q - ($Qnew - $fnow);
        $fnow = $f[++$findex];
      }
      $Q = $Qnew;
      if($hh !== 0.0) {
        $h[$hindex++] = $hh;
      }
      while($eindex < $elen && $findex < $flen) {
        if(($fnow > $enow) === ($fnow > -$enow)) {
          $Qnew = $Q + $enow;
          $bvirt = $Qnew - $Q;
          $hh = $Q - ($Qnew - $bvirt) + ($enow - $bvirt);
          $enow = $e[++$eindex];
        } else {
          $Qnew = $Q + $fnow;
          $bvirt = $Qnew - $Q;
          $hh = $Q - ($Qnew - $bvirt) + ($fnow - $bvirt);
          $fnow = $f[++$findex];
        }
        $Q = $Qnew;
        if($hh !== 0.0) {
          $h[$hindex++] = $hh;
        }
      }
    }
    while($eindex < $elen) {
      $Qnew = $Q + $enow;
      $bvirt = $Qnew - $Q;
      $hh = $Q - ($Qnew - $bvirt) + ($enow - $bvirt);
      $enow = $e[++$eindex];
      $Q = $Qnew;
      if($hh !== 0.0) {
        $h[$hindex++] = $hh;
      }
    }
    while($findex < $flen) {
      $Qnew = $Q + $fnow;
      $bvirt = $Qnew - $Q;
      $hh = $Q - ($Qnew - $bvirt) + ($fnow - $bvirt);
      $fnow = $f[++$findex];
      $Q = $Qnew;
      if($hh !== 0.0) {
        $h[$hindex++] = $hh;
      }
    }
    if($Q !== 0.0 || $hindex === 0) {
      $h[$hindex++] = $Q;
    }
    return $hindex;
  }

  /**
   * @param int $elen
   * @param array<float> $e
   * @return float
   */
  public static function estimate(int $elen, array $e): float {
    $Q = $e[0];
    for($i = 1; $i < $elen; $i++) {
      $Q += $e[$i];
    }
    return $Q;
  }

  /**
   * Return `new Float64Array(n)` in original JS code.
   *
   * @param int $n
   * @return array<float>
   */
  public static function vec(int $n): array {
    return array_fill(0, $n, 0.0);
  }
}
