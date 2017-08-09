package de.sb.toolbox.math;

import static java.lang.Math.E;
import static java.lang.Math.PI;
import java.util.function.UnaryOperator;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.math.Complex.MutableDoublePrecision;
import de.sb.toolbox.math.Complex.MutableSinglePrecision;
import javafx.scene.image.WritableImage;
import javafx.scene.paint.Color;


/**
 * This facade provides additional non-polymorphic operations for mutable complex arguments, both for single and double
 * precision.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public final class ComplexMath {

	/**
	 * Prevents external instantiation.
	 */
	private ComplexMath () {}


	//************************//
	// exponential operations //
	//************************//

	/**
	 * Returns <tt>z<sup>1/n</sup></tt>, i.e. the principal branch solution for the finitely multi-valued {@code integer} root
	 * of a {@code complex} operand {@code z}, using de Moivre's formula. The multi-valued nature of this operation implies a
	 * discontinuity, which is located at negative real values of the operand {@code z}, except for <tt>n=&plusmn;1</tt> where
	 * this operation is single-valued, and {@code n=0} where no valid branches exist (results in {@code NaN}).<br />
	 * <br />
	 * The principal branch solution is the one whose argument (positive or negative) is closest to zero; it's real part is
	 * guaranteed to be positive for {@code |n|>1}. The solutions for the other branches can be obtained by multiplying the
	 * principal branch solution with <tt>e<sup>i&middot;&tau;&middot;k/n</sup></i></tt>, using branch index
	 * <tt>k&isin;&#x2124;</tt> within range <tt>[0,|n|[</tt>. When combining this operation with it's inverse, as in
	 * <tt>(z<sup>1/n</sup>)<sup>n</sup> = (z<sup>n</sup>)<sup>1/n</sup> = z</tt>, use<br />
	 * <tt>k = (floor(|z|&middot;n/&tau;+&half;)%|n|</tt>, adding {@code |n|} if {@code k<0}.<br />
	 * <br />
	 * Note that fractional exponentiation may be performed by first using this operation to calculate the fraction's divisor
	 * root of {@code z} (including choice of the proper branch), and subsequently raising the result to the power of the
	 * fraction's dividend using {@link #pow(T, long)}.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @param n the integer root
	 * @return the principal value <tt>z<sup>1/n</sup> = </tt> <tt>|z|<sup>1/n</sup>&middot;e<sup>i&middot;&phi;(z)/n</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T root (final T z, long n) throws NullPointerException {
		if (n == (byte) n) {
			switch ((byte) n) {
				case -3:
					return z.clone().cbrt().inv();
				case -2:
					return z.clone().sqrt().inv();
				case +2:
					return z.clone().sqrt();
				case +3:
					return z.clone().cbrt();
				default:
					break;
			}
		}
		return pow(z, 1f / n);
	}


	/**
	 * Returns <tt>z<sup>1/n</sup></tt>, i.e. the principal branch solution for the finitely multi-valued {@code integer} root
	 * of a {@code complex} operand {@code z}, using de Moivre's formula. The multi-valued nature of this operation implies a
	 * discontinuity, which is located at negative real values of the operand {@code z}, except for <tt>n=&plusmn;1</tt> where
	 * this operation is single-valued, and {@code n=0} where no valid branches exist (results in {@code NaN}).<br />
	 * <br />
	 * The principal branch solution is the one whose argument (positive or negative) is closest to zero; it's real part is
	 * guaranteed to be positive for {@code |n|>1}. The solutions for the other branches can be obtained by multiplying the
	 * principal branch solution with <tt>e<sup>i&middot;&tau;&middot;k/n</sup></i></tt>, using branch index
	 * <tt>k&isin;&#x2124;</tt> within range <tt>[0,|n|[</tt>. When combining this operation with it's inverse, as in
	 * <tt>(z<sup>1/n</sup>)<sup>n</sup> = (z<sup>n</sup>)<sup>1/n</sup> = z</tt>, use<br />
	 * <tt>k = (floor(|z|&middot;n/&tau;+&half;)%|n|</tt>, adding {@code |n|} if {@code k<0}.<br />
	 * <br />
	 * Note that fractional exponentiation may be performed by first using this operation to calculate the fraction's divisor
	 * root of {@code z} (including choice of the proper branch), and subsequently raising the result to the power of the
	 * fraction's dividend using {@link #pow(T, long)}.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @param n the integer root
	 * @return the principal value <tt>z<sup>1/n</sup> = </tt> <tt>|z|<sup>1/n</sup>&middot;e<sup>i&middot;&phi;(z)/n</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T root (final T z, long n) throws NullPointerException {
		if (n == (byte) n) {
			switch ((byte) n) {
				case -3:
					return z.clone().cbrt().inv();
				case -2:
					return z.clone().sqrt().inv();
				case +2:
					return z.clone().sqrt();
				case +3:
					return z.clone().cbrt();
				default:
					break;
			}
		}
		return pow(z, 1d / n);
	}


	/**
	 * Returns <tt>z<sup>n</sup></tt>, i.e. a {@code complex} base raised to the power of an {@code integer} exponent, using de
	 * Moivre's formula. The special case <tt>z<sup>0</sup></tt> calculates to {@code 1} in all cases, even for zero, infinite
	 * and {@code NaN} values of {@code z}, as defined in {@link Math#pow}.<br />
	 * <br />
	 * Note that this operation is finitely cyclic, with
	 * <tt>(z&middot;e<sup>i&middot;&tau;&middot;k/n</sup>)<sup>n</sup> = z<sup>n</sup></tt> for every integer {@code k} within
	 * range {@code [0,|n|[}. This cyclic behavior causes the root operation (as the inverse of this operation) to become
	 * multi-valued in <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to
	 * become multi-valued in <tt>&#x211D;</tt>.<br />
	 * <br />
	 * Also note that fractional exponentiation may be performed by first using {@link #root(T, long)} to calculate the
	 * fraction's divisor root of {@code z} (including choice of the proper branch), and subsequently using this operation to
	 * raise the result to the power of the fraction's dividend.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex base
	 * @param n the integer exponent
	 * @return the value <tt>z<sup>n</sup> = </tt> <tt>|z|<sup>n</sup>&middot;e<sup>i&middot;&phi;(z)&middot;n</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T pow (final T z, final long n) throws NullPointerException {
		if (n == (byte) n) {
			switch ((byte) n) {
				case -3:
					return z.clone().cb().inv();
				case -2:
					return z.clone().sq().inv();
				case +2:
					return z.clone().sq();
				case +3:
					return z.clone().cb();
				default:
					break;
			}
		}
		return pow(z, (float) n);
	}


	/**
	 * Returns <tt>z<sup>n</sup></tt>, i.e. a {@code complex} base raised to the power of an {@code integer} exponent, using de
	 * Moivre's formula. The special case <tt>z<sup>0</sup></tt> calculates to {@code 1} in all cases, even for zero, infinite
	 * and {@code NaN} values of {@code z}, as defined in {@link Math#pow}.<br />
	 * <br />
	 * Note that this operation is finitely cyclic, with
	 * <tt>(z&middot;e<sup>i&middot;&tau;&middot;k/n</sup>)<sup>n</sup> = z<sup>n</sup></tt> for every integer {@code k} within
	 * range {@code [0,|n|[}. This cyclic behavior causes the root operation (as the inverse of this operation) to become
	 * multi-valued in <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to
	 * become multi-valued in <tt>&#x211D;</tt>.<br />
	 * <br />
	 * Also note that fractional exponentiation may be performed by first using {@link #root(T, long)} to calculate the
	 * fraction's divisor root of {@code z} (including choice of the proper branch), and subsequently using this operation to
	 * raise the result to the power of the fraction's dividend.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex base
	 * @param n the integer exponent
	 * @return the value <tt>z<sup>n</sup> = </tt> <tt>|z|<sup>n</sup>&middot;e<sup>i&middot;&phi;(z)&middot;n</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T pow (final T z, final long n) throws NullPointerException {
		if (n == (byte) n) {
			switch ((byte) n) {
				case -3:
					return z.clone().cb().inv();
				case -2:
					return z.clone().sq().inv();
				case +2:
					return z.clone().sq();
				case +3:
					return z.clone().cb();
				default:
					break;
			}
		}
		return pow(z, (double) n);
	}


	/**
	 * Returns <tt>z<sup>r</sup></tt>, i.e. a {@code complex} base raised to the power of a {@code real} exponent, using de
	 * Moivre's formula. The special case <tt>z<sup>0</sup></tt> calculates to {@code 1} in all cases, even for zero, infinite
	 * and {@code NaN} values of {@code z}, as defined in {@link Math#pow}.<br />
	 * <br />
	 * This operation can be both infinitely cyclic and infinitely multi-valued, depending on the value of the given exponent.
	 * Therefore, it is practically unsuitable for branch analysis, for which both {@link #root(T, long)} and
	 * {@link #pow(T, long)} should be preferred. Note that the multi-valued nature of this operation implies a discontinuity,
	 * which is located at negative real values of the operand {@code z}, except for <tt>r=&plusmn;1</tt> where this operation
	 * is single-valued, and <tt>r&rarr;&plusmn;&infin;</tt> where no valid branches exist (results in {@code NaN}).
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex base
	 * @param r the real exponent
	 * @return the value <tt>z<sup>n</sup> = </tt> <tt>|z|<sup>r</sup>&middot;e<sup>i&middot;&phi;(z)&middot;r</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T pow (final T z, final float r) throws NullPointerException {
		if (r == Double.NEGATIVE_INFINITY | r == Double.POSITIVE_INFINITY) return z.isZero() ? z.clone().setCartesian(Float.NaN, Float.NaN) : z.clone().setCartesian(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
		if (r == -1) return z.clone().inv();
		if (r == 0) return z.clone().setCartesian(1, 0);
		if (r == -1) return z;
		return z.clone().setPolar((float) Math.pow(z.abs(), r), z.arg() * r);
	}


	/**
	 * Returns <tt>z<sup>r</sup></tt>, i.e. a {@code complex} base raised to the power of a {@code real} exponent, using de
	 * Moivre's formula. The special case <tt>z<sup>0</sup></tt> calculates to {@code 1} in all cases, even for zero, infinite
	 * and {@code NaN} values of {@code z}, as defined in {@link Math#pow}.<br />
	 * <br />
	 * This operation can be both infinitely cyclic and infinitely multi-valued, depending on the value of the given exponent.
	 * Therefore, it is practically unsuitable for branch analysis, for which both {@link #root(T, long)} and
	 * {@link #pow(T, long)} should be preferred. Note that the multi-valued nature of this operation implies a discontinuity,
	 * which is located at negative real values of the operand {@code z}, except for <tt>r=&plusmn;1</tt> where this operation
	 * is single-valued, and <tt>r&rarr;&plusmn;&infin;</tt> where no valid branches exist (results in {@code NaN}).
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex base
	 * @param r the real exponent
	 * @return the value <tt>z<sup>n</sup> = </tt> <tt>|z|<sup>r</sup>&middot;e<sup>i&middot;&phi;(z)&middot;r</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T pow (final T z, final double r) throws NullPointerException {
		if (r == Double.NEGATIVE_INFINITY | r == Double.POSITIVE_INFINITY) return z.isZero() ? z.clone().setCartesian(Double.NaN, Double.NaN) : z.clone().setCartesian(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		if (r == -1) return z.clone().inv();
		if (r == 0) return z.clone().setCartesian(1, 0);
		if (r == +1) return z;
		return z.clone().setPolar(Math.pow(z.abs(), r), z.arg() * r);
	}


	/**
	 * Returns {@code log(z)}, i.e. the principal branch solution for the infinitely multi-valued natural logarithm (base
	 * {@code e}) of a {@code complex} operand. The multi-valued nature implies a discontinuity, which is located at negative
	 * real values of the operand {@code z} .</br />
	 * <br />
	 * The principal branch solution is the one whose imaginary part is guaranteed to be within range <tt>[-&pi;,+&pi;[</tt>.
	 * The solutions for the other branches may be obtained by adding <tt>i&middot;&tau;&middot;k</tt> (branch Index {@code k}
	 * &isin; &#x2124;) to the principal branch solution. Regardless of the branch chosen, the real part of the result is
	 * guaranteed to be positive.<br />
	 * <br />
	 * When combining this operation with it's inverse, as in <tt>log(exp(z)) = exp(log(z)) =
	 * z</tt>, choose branch <tt>k = floor(&image;(z)/&tau; + &half;)</tt> for correct branch alignment.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the principal value <tt>log(z) = log|z|+i&middot;&phi;(z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T log (final T z) throws NullPointerException {
		return z.clone().setCartesian((float) Math.log(z.abs()), z.arg());
	}


	/**
	 * Returns {@code log(z)}, i.e. the principal branch solution for the infinitely multi-valued natural logarithm (base
	 * {@code e}) of a {@code complex} operand. The multi-valued nature implies a discontinuity, which is located at negative
	 * real values of the operand {@code z} .</br />
	 * <br />
	 * The principal branch solution is the one whose imaginary part is guaranteed to be within range <tt>[-&pi;,+&pi;[</tt>.
	 * The solutions for the other branches may be obtained by adding <tt>i&middot;&tau;&middot;k</tt> (branch Index {@code k}
	 * &isin; &#x2124;) to the principal branch solution. Regardless of the branch chosen, the real part of the result is
	 * guaranteed to be positive.<br />
	 * <br />
	 * When combining this operation with it's inverse, as in <tt>log(exp(z)) = exp(log(z)) =
	 * z</tt>, choose branch <tt>k = floor(&image;(z)/&tau; + &half;)</tt> for correct branch alignment.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the principal value <tt>log(z) = log|z|+i&middot;&phi;(z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T log (final T z) throws NullPointerException {
		return z.clone().setCartesian(Math.log(z.abs()), z.arg());
	}


	/**
	 * Returns <tt>log<sub>b</sub>(z)</tt>, i.e. the principal branch solution for the infinitely multi-valued {@code real}
	 * based logarithm of a {@code complex} operand. The multi-valued nature implies a discontinuity, which is located at
	 * negative real values of the operand {@code z}. The exact behavior of this operation depends on the nature of the base
	 * {@code b}:
	 * <ul>
	 * <li>if a positive base {@code b} equals the operand {@code z}, then the result is always one.</li>
	 * <li>if the base {@code b} is {@code 0}, then the result is zero.</li>
	 * <li>if the base {@code b} is {@code +1}, then the result is infinite and/or NaN.</li>
	 * <li>if the base {@code b} is Euler's number {@code e}, then see {@link #log(T)} for details.</li>
	 * <li>if the base {@code b} is any other positive number, then the principal branch solution is the one whose imaginary
	 * part is guaranteed to be within range <i>[-&pi;/log(b),+&pi;/log(b)[</i>. The solutions for the other branches may be
	 * obtained by adding <i>i&middot;&tau;&middot;k/log(b)</i> (branch index <tt>k &isin; &#x2124;</tt>) to the principal
	 * branch solution. Regardless of the branch chosen, the real part of the result is guaranteed to be positive.</li>
	 * <li>if the base {@code b} is {@code -1}, then the solutions for the other branches may be obtained by adding
	 * <tt>2&middot;k</tt> (with branch index <tt>k &isin; &#x2124;</tt>) to the principal branch solution.</li>
	 * <li>if the base {@code b} is any other negative number, then the solutions for the other branches are also obtainable,
	 * but the branch correction becomes more complicated.</li>
	 * </ul>
	 * When combining this operation with it's inverse, as in <i>log(b, exp(b,z)) = exp(b, log(b,z)) = z</i>, choose
	 * <ul>
	 * <li>if the base {@code b} is positive: <tt>k = floor(&half;+&image;(z)/&tau;&middot;log(b))</tt></li>
	 * <li>if the base {@code b = -1}:
	 * <ul>
	 * <li><tt>k = floor(&half;(&real;(z)+1))</tt> if &real;(z) < -1</li>
	 * <li><tt>k = ceil(&half;(&real;(z)-1))</tt> if &real;(z) >= -1</li>
	 * </ul>
	 * </ul>
	 * @param <T> the declaration type of the complex argument and result
	 * @param b the real base
	 * @param z the complex operand
	 * @return the principal value <tt>log<sub>b</sub>(z) = log(z)/log(b) = </tt> <tt>(log|z|+i&middot;&phi;(z))/log(b)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T log (final float b, final T z) throws NullPointerException {
		if (b == (float) E) return log(z);
		if (z.isReal() && z.re() == b) return z.clone().setCartesian(1, 0);

		final float logB = (float) Math.log(Math.abs(b));
		final T logZ = log(z);
		return b >= 0 ? logZ.div(logB) : logZ.div(z.clone().setCartesian(logB, (float) -PI));
	}


	/**
	 * Returns <tt>log<sub>b</sub>(z)</tt>, i.e. the principal branch solution for the infinitely multi-valued {@code real}
	 * based logarithm of a {@code complex} operand. The multi-valued nature implies a discontinuity, which is located at
	 * negative real values of the operand {@code z}. The exact behavior of this operation depends on the nature of the base
	 * {@code b}:
	 * <ul>
	 * <li>if a positive base {@code b} equals the operand {@code z}, then the result is always one.</li>
	 * <li>if the base {@code b} is {@code 0}, then the result is zero.</li>
	 * <li>if the base {@code b} is {@code +1}, then the result is infinite and/or NaN.</li>
	 * <li>if the base {@code b} is Euler's number {@code e}, then see {@link #log(T)} for details.</li>
	 * <li>if the base {@code b} is any other positive number, then the principal branch solution is the one whose imaginary
	 * part is guaranteed to be within range <i>[-&pi;/log(b),+&pi;/log(b)[</i>. The solutions for the other branches may be
	 * obtained by adding <i>i&middot;&tau;&middot;k/log(b)</i> (branch index <tt>k &isin; &#x2124;</tt>) to the principal
	 * branch solution. Regardless of the branch chosen, the real part of the result is guaranteed to be positive.</li>
	 * <li>if the base {@code b} is {@code -1}, then the solutions for the other branches may be obtained by adding
	 * <tt>2&middot;k</tt> (with branch index <tt>k &isin; &#x2124;</tt>) to the principal branch solution.</li>
	 * <li>if the base {@code b} is any other negative number, then the solutions for the other branches are also obtainable,
	 * but the branch correction becomes more complicated.</li>
	 * </ul>
	 * When combining this operation with it's inverse, as in <i>log(b, exp(b,z)) = exp(b, log(b,z)) = z</i>, choose
	 * <ul>
	 * <li>if the base {@code b} is positive: <tt>k = floor(&half;+&image;(z)/&tau;&middot;log(b))</tt></li>
	 * <li>if the base {@code b = -1}:
	 * <ul>
	 * <li><tt>k = floor(&half;(&real;(z)+1))</tt> if &real;(z) < -1</li>
	 * <li><tt>k = ceil(&half;(&real;(z)-1))</tt> if &real;(z) >= -1</li>
	 * </ul>
	 * </ul>
	 * @param <T> the declaration type of the complex argument and result
	 * @param b the real base
	 * @param z the complex operand
	 * @return the principal value <tt>log<sub>b</sub>(z) = log(z)/log(b) = </tt> <tt>(log|z|+i&middot;&phi;(z))/log(b)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T log (final double b, final T z) throws NullPointerException {
		if (b == E) return log(z);
		if (z.isReal() && z.re() == b) return z.clone().setCartesian(1, 0);

		final double logB = Math.log(Math.abs(b));
		final T logZ = log(z);
		return b >= 0 ? logZ.div(logB) : logZ.div(z.clone().setCartesian(logB, -PI));
	}


	/**
	 * Returns <tt>log<sub>z<sub>1</sub></sub>(z<sub>2</sub>)</tt>, i.e. the principal branch solution for the infinitely
	 * multi-valued {@code complex} based logarithm of a {@code complex} operand. The multi-valued nature implies a
	 * discontinuity, which is located at negative real values of the operand <tt>z<sub>2</sub></tt>. The exact behavior of this
	 * operation depends on the nature of the base <tt>z<sub>1</sub></tt>:
	 * <ul>
	 * <li>if the base z<sub>1</sub> equals the operand <tt>z<sub>2</sub></tt>, then the result is always one.</li>
	 * <li>if the base z<sub>1</sub> is Euler's number {@code e}, then see {@link #log(T)} for details.</li>
	 * <li>if the base z<sub>1</sub> is a {@code real} number, then see {@link #log(float, T)} for details.</li>
	 * <li>if the base z<sub>1</sub> is a {@code complex} number, the branching behavior of this operation is similar to
	 * {@link #log(float, T)} when passing negative bases. However, branch determination and correction become much more
	 * complicated.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex arguments and result
	 * @param z1 the complex base
	 * @param z2 the complex operand
	 * @return the principal value <tt>log<sub>z<sub>1</sub></sub>(z<sub>2</sub>) =
	 *         log(z<sub>2</sub>)/log(z<sub>1</sub>) = 
	 *         (log|z<sub>2</sub>|+i&middot;&phi;(z<sub>2</sub>)) /
	 *         (log|z<sub>1</sub>|+i&middot;&phi;(z<sub>1</sub>))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T log (final T z1, final T z2) throws NullPointerException {
		if (z1.isReal()) return log(z1.re(), z2);
		if (z1 == z2 || z1.equals(z2)) return z2.clone().setCartesian(1, 0);
		return log(z2).div(log(z1));
	}


	/**
	 * Returns <tt>log<sub>z<sub>1</sub></sub>(z<sub>2</sub>)</tt>, i.e. the principal branch solution for the infinitely
	 * multi-valued {@code complex} based logarithm of a {@code complex} operand. The multi-valued nature implies a
	 * discontinuity, which is located at negative real values of the operand <tt>z<sub>2</sub></tt>. The exact behavior of this
	 * operation depends on the nature of the base <tt>z<sub>1</sub></tt>:
	 * <ul>
	 * <li>if the base z<sub>1</sub> equals the operand <tt>z<sub>2</sub></tt>, then the result is always one.</li>
	 * <li>if the base z<sub>1</sub> is Euler's number {@code e}, then see {@link #log(T)} for details.</li>
	 * <li>if the base z<sub>1</sub> is a {@code real} number, then see {@link #log(double, T)} for details.</li>
	 * <li>if the base z<sub>1</sub> is a {@code complex} number, the branching behavior of this operation is similar to
	 * {@link #log(double, T)} when passing negative bases. However, branch determination and correction become much more
	 * complicated.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex arguments and result
	 * @param z1 the complex base
	 * @param z2 the complex operand
	 * @return the principal value <tt>log<sub>z<sub>1</sub></sub>(z<sub>2</sub>) =
	 *         log(z<sub>2</sub>)/log(z<sub>1</sub>) = 
	 *         (log|z<sub>2</sub>|+i&middot;&phi;(z<sub>2</sub>)) /
	 *         (log|z<sub>1</sub>|+i&middot;&phi;(z<sub>1</sub>))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T log (final T z1, final T z2) throws NullPointerException {
		if (z1.isReal()) return log(z1.re(), z2);
		if (z1 == z2 || z1.equals(z2)) return z2.clone().setCartesian(1, 0);
		return log(z2).div(log(z1));
	}


	/**
	 * Returns <tt>e<sup>z</sup></tt>, i.e. Euler's number raised to the power of a {@code complex} operand.<br />
	 * <br />
	 * Note that this operation is infinitely cyclic, with <tt>e<sup>z+i&middot;&tau;&middot;k</sup> = e<sup>z</sup></tt> for
	 * any integer k. This behavior causes the logarithm (as the inverse of this operation) to become multi-valued in
	 * <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to become multi-valued
	 * in <tt>&#x211D;</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex exponent
	 * @return the value <tt>e<sup>z</sup> = </tt> <tt>e<sup>&real;(z)</sup>&middot;e<sup>i&middot;&image;(z)</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T exp (final T z) throws NullPointerException {
		return z.clone().setPolar((float) Math.exp(z.re()), z.im());
	}


	/**
	 * Returns <tt>e<sup>z</sup></tt>, i.e. Euler's number raised to the power of a {@code complex} operand.<br />
	 * <br />
	 * Note that this operation is infinitely cyclic, with <tt>e<sup>z+i&middot;&tau;&middot;k</sup> = e<sup>z</sup></tt> for
	 * any integer k. This behavior causes the logarithm (as the inverse of this operation) to become multi-valued in
	 * <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to become multi-valued
	 * in <tt>&#x211D;</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex exponent
	 * @return the value <tt>e<sup>z</sup> = </tt> <tt>e<sup>&real;(z)</sup>&middot;e<sup>i&middot;&image;(z)</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T exp (final T z) throws NullPointerException {
		return z.clone().setPolar(Math.exp(z.re()), z.im());
	}


	/**
	 * Returns <tt>b<sup>z</sup></tt>, i.e. a {@code real} base raised to the power of a {@code complex} exponent. The result of
	 * <tt>1<sup>z</sup></tt> and <tt>b<sup>0</sup></tt> is always one, by definition including the special case
	 * <tt>0<sup>0</sup></tt> . The result of <tt>0<sup>z</sup></tt> is zero unless {@code z} is zero.<br />
	 * <br />
	 * Note that this operation is infinitely cyclic except for special bases like zero, one or infinity, and special exponents
	 * like 0 or infinity. This behavior causes the logarithm (as the inverse of this operation) to become multi-valued in
	 * <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to become multi-valued
	 * in <tt>&#x211D;</tt>. The exact cycle behavior depends on the nature of the base {@code b}:
	 * <ul>
	 * <li>if the given base {@code b} is strictly positive except {@code +1}, then the cycling is simply vertical, with<br />
	 * <i>exp(b, z + i&middot;&tau;&middot;k / log(b)) = exp(b, z)</i> for any integer {@code k}.</li>
	 * <li>if the given base is {@code b=-1}, then the cycling is simply horizontal, with <i>exp(b, z + 2&middot;k) = exp(b,
	 * z)</i> for any integer {@code k}.</li>
	 * <li>if the given base {@code b} is strictly negative except {@code -1}, then the cycling is obliquely rotated from the
	 * vertical by an angle between ]0,&pi;[, depending on the absolute value of {@code b}.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex argument and result
	 * @param b the (positive) real base
	 * @param z the complex exponent
	 * @return the value <tt>b<sup>z</sup> = </tt><tt>e<sup>log(b)&middot;z</sup> = </tt>
	 *         <tt>e<sup>log(b)&middot;&real;(z)</sup>&middot;e<sup>i&middot;log(b)&middot;&image;(z)</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T exp (final float b, final T z) throws NullPointerException {
		final T clone = z.clone();

		if (b < 0) return exp(clone.setCartesian(b, 0), z);
		if (b == 0) return clone.setCartesian(z.isZero() ? 1 : 0, 0);
		if (b == 1) return clone.setCartesian(1, 0);
		if (b == (float) E) return exp(z);

		final float logB = (float) Math.log(b);
		return clone.setPolar((float) Math.exp(logB * z.re()), logB * z.im());
	}


	/**
	 * Returns <tt>b<sup>z</sup></tt>, i.e. a {@code real} base raised to the power of a {@code complex} exponent. The result of
	 * <tt>1<sup>z</sup></tt> and <tt>b<sup>0</sup></tt> is always one, by definition including the special case
	 * <tt>0<sup>0</sup></tt> . The result of <tt>0<sup>z</sup></tt> is zero unless {@code z} is zero.<br />
	 * <br />
	 * Note that this operation is infinitely cyclic except for special bases like zero, one or infinity, and special exponents
	 * like 0 or infinity. This behavior causes the logarithm (as the inverse of this operation) to become multi-valued in
	 * <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the arc sine function to become multi-valued
	 * in <tt>&#x211D;</tt>. The exact cycle behavior depends on the nature of the base {@code b}:
	 * <ul>
	 * <li>if the given base {@code b} is strictly positive except {@code +1}, then the cycling is simply vertical, with<br />
	 * <i>exp(b, z + i&middot;&tau;&middot;k / log(b)) = exp(b, z)</i> for any integer {@code k}.</li>
	 * <li>if the given base is {@code b=-1}, then the cycling is simply horizontal, with <i>exp(b, z + 2&middot;k) = exp(b,
	 * z)</i> for any integer {@code k}.</li>
	 * <li>if the given base {@code b} is strictly negative except {@code -1}, then the cycling is obliquely rotated from the
	 * vertical by an angle between ]0,&pi;[, depending on the absolute value of {@code b}.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex argument and result
	 * @param b the (positive) real base
	 * @param z the complex exponent
	 * @return the value <tt>b<sup>z</sup> = </tt><tt>e<sup>log(b)&middot;z</sup> = </tt>
	 *         <tt>e<sup>log(b)&middot;&real;(z)</sup>&middot;e<sup>i&middot;log(b)&middot;&image;(z)</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T exp (final double b, final T z) throws NullPointerException {
		final T clone = z.clone();

		if (b < 0) return exp(clone.setCartesian(b, 0), z);
		if (b == 0) return clone.setCartesian(z.isZero() ? 1 : 0, 0);
		if (b == 1) return clone.setCartesian(1, 0);
		if (b == E) return exp(z);

		final double logB = Math.log(b);
		return clone.setPolar(Math.exp(logB * z.re()), logB * z.im());
	}


	/**
	 * Returns <tt>z<sub>1</sub><sup>z<sub>2</sub></sup></tt>, i.e. the principal branch solution for a {@code complex} base
	 * raised to the power of a {@code complex} exponent, which is infinitely multi-valued in <tt>&#x2102;</tt>. The exact
	 * behavior of this operation depends on the nature of the operands:
	 * <ul>
	 * <li>If the base <tt>z<sub>1</sub></tt> is positive {@code real} number, then see {@link #exp(float, T)} for details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is an {@code integer} number, then see {@link #pow(T, int)} for details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is or approximates a unit fraction, then see {@link #root(T, int)} for
	 * details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is a {@code real} number, but neither an integer nor a unit fraction, then it
	 * is technically nevertheless a dyadic fraction in lowest terms because any IEEE floating point number is. Therefore, this
	 * operation behaves as if {@link #root(T, int)} and {@link #pow(T, int)} had been combined, but without the guarantee that
	 * the principal branch solution returned is the one whose argument (positive or negative) is closest to zero.</li>
	 * <li>Otherwise this operation displays a more delicate behavior: If the base is set constant and the exponent is varied,
	 * then this operation shows a cyclic pattern similarly to {@link #pow(T, int)}, but rotated by varying degrees. If the base
	 * is varied and the exponent set constant, then this operation displays a pattern similarly to {@link #root(T, int)}, but
	 * rotated by varying degrees resulting in cycle or spiral patterns. If z<sub>1</sub> = z <sub>2</sub> or z<sub>1</sub> = -z
	 * <sub>2</sub>, then this operation displays a cyclic pattern similarly to {@link #pow(T, int)}, but bended around the
	 * origin. If z<sub>1</sub> = z <sub>2</sub><sup>-1</sup>, then this operation displays a cyclic shamrock pattern.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex arguments and result
	 * @param z1 the complex base
	 * @param z2 the complex exponent
	 * @return the value <tt>z<sub>1</sub><sup>z<sub>2</sub></sup></tt> = </tt>
	 *         <tt>e<sup>log(z<sub>1</sub>)&middot;z<sub>2</sub></sup> = </tt>
	 *         <tt>e<sup>(log|z<sub>1</sub>|+i&middot;&phi;(z<sub>1</sub>)) &middot;
	 *         (&real;(z<sub>2</sub>)+i&middot;&image;(z<sub>2</sub>))</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T exp (final T z1, final T z2) throws NullPointerException {
		if (z1.isReal()) {
			final float re = z1.re();
			if (re >= 0) return exp(re, z2);
		}

		if (z2.isReal()) {
			final float re = z2.re();
			if (re == (long) re) return pow(z1, (long) re);
			final float inv = 1 / re;
			if (inv == (long) inv) return root(z1, (long) inv);
		}

		return exp(log(z1).mul(z2));
	}


	/**
	 * Returns <tt>z<sub>1</sub><sup>z<sub>2</sub></sup></tt>, i.e. the principal branch solution for a {@code complex} base
	 * raised to the power of a {@code complex} exponent, which is infinitely multi-valued in <tt>&#x2102;</tt>. The exact
	 * behavior of this operation depends on the nature of the operands:
	 * <ul>
	 * <li>If the base <tt>z<sub>1</sub></tt> is positive {@code real} number, then see {@link #exp(double, T)} for details.
	 * </li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is an {@code integer} number, then see {@link #pow(T, int)} for details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is or approximates a unit fraction, then see {@link #root(T, int)} for
	 * details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is a {@code real} number, but neither an integer nor a unit fraction, then it
	 * is technically nevertheless a dyadic fraction in lowest terms because any IEEE floating point number is. Therefore, this
	 * operation behaves as if {@link #root(T, int)} and {@link #pow(T, int)} had been combined, but without the guarantee that
	 * the principal branch solution returned is the one whose argument (positive or negative) is closest to zero.</li>
	 * <li>Otherwise this operation displays a more delicate behavior: If the base is set constant and the exponent is varied,
	 * then this operation shows a cyclic pattern similarly to {@link #pow(T, int)}, but rotated by varying degrees. If the base
	 * is varied and the exponent set constant, then this operation displays a pattern similarly to {@link #root(T, int)}, but
	 * rotated by varying degrees resulting in cycle or spiral patterns. If z<sub>1</sub> = z <sub>2</sub> or z<sub>1</sub> = -z
	 * <sub>2</sub>, then this operation displays a cyclic pattern similarly to {@link #pow(T, int)}, but bended around the
	 * origin. If z<sub>1</sub> = z <sub>2</sub><sup>-1</sup>, then this operation displays a cyclic shamrock pattern.</li>
	 * </ul>
	 * @param <T> the declaration type of the complex arguments and result
	 * @param z1 the complex base
	 * @param z2 the complex exponent
	 * @return the value <tt>z<sub>1</sub><sup>z<sub>2</sub></sup></tt> = </tt>
	 *         <tt>e<sup>log(z<sub>1</sub>)&middot;z<sub>2</sub></sup> = </tt>
	 *         <tt>e<sup>(log|z<sub>1</sub>|+i&middot;&phi;(z<sub>1</sub>)) &middot;
	 *         (&real;(z<sub>2</sub>)+i&middot;&image;(z<sub>2</sub>))</sup></tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T exp (final T z1, final T z2) throws NullPointerException {
		if (z1.isReal()) {
			final double re = z1.re();
			if (re >= 0) return exp(re, z2);
		}

		if (z2.isReal()) {
			final double re = z2.re();
			if (re == (long) re) return pow(z1, (long) re);
			final double inv = 1 / re;
			if (inv == (long) inv) return root(z1, (long) inv);
		}

		return exp(log(z1).mul(z2));
	}


	//**************************//
	// trigonometric operations //
	//**************************//

	/**
	 * Returns <tt>sin(z) = -i&middot;&half;(e<sup>+iz</sup>-e<sup>-iz</sup>) = </tt>
	 * <tt>sin(&real;(z))&middot;cosh(&image;(z)) + i&middot;cos(&real;(z))&middot;sinh(&image;(z))</tt>
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>sin(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) + </tt>
	 *         <tt>i&middot;cos(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T sin (final T z) throws NullPointerException {
		final float sin = (float) Math.sin(z.re()), cos = (float) Math.cos(z.re());
		final float pex = (float) Math.exp(z.im()), nex = 1 / pex;
		return z.clone().setCartesian(+.5f * (pex + nex) * sin, +.5f * (pex - nex) * cos);
	}


	/**
	 * Returns <tt>sin(z) = -i&middot;&half;(e<sup>+iz</sup>-e<sup>-iz</sup>) = </tt>
	 * <tt>sin(&real;(z))&middot;cosh(&image;(z)) + i&middot;cos(&real;(z))&middot;sinh(&image;(z))</tt>
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>sin(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) + </tt>
	 *         <tt>i&middot;cos(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T sin (final T z) throws NullPointerException {
		//		return z.clone().setCartesian(+Math.sin(z.re()) * Math.cosh(z.im()), +Math.cos(z.re()) * Math.sinh(z.im()));
		final double sin = Math.sin(z.re()), cos = Math.cos(z.re());
		final double pex = Math.exp(z.im()), nex = 1 / pex;
		return z.clone().setCartesian(+.5d * (pex + nex) * sin, +.5d * (pex - nex) * cos);
	}


	/**
	 * Returns <tt>cos(z) = &half;(e<sup>+iz</sup>+e<sup>-iz</sup>) = </tt>
	 * <tt>cos(&real;(z))&middot;cosh(&image;(z)) - i&middot;sin(&real;(z))&middot;sinh(&image;(z))</tt>
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>cos(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) + </tt>
	 *         <tt>i&middot;sin(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T cos (final T z) throws NullPointerException {
		final float sin = (float) Math.sin(z.re()), cos = (float) Math.cos(z.re());
		final float pex = (float) Math.exp(z.im()), nex = 1 / pex;
		return z.clone().setCartesian(+.5f * (pex + nex) * cos, -.5f * (pex - nex) * sin);
	}


	/**
	 * Returns <tt>cos(z) = &half;(e<sup>+iz</sup>+e<sup>-iz</sup>) = </tt>
	 * <tt>cos(&real;(z))&middot;cosh(&image;(z)) - i&middot;sin(&real;(z))&middot;sinh(&image;(z))</tt>
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>cos(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) - </tt>
	 *         <tt>i&middot;sin(&real;(z))&middot;&half;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T cos (final T z) throws NullPointerException {
		// return z.clone().setCartesian(+Math.cos(z.re()) * Math.cosh(z.im()),-Math.sin(z.re()) * Math.sinh(z.im()));
		final double sin = Math.sin(z.re()), cos = Math.cos(z.re());
		final double pex = Math.exp(z.im()), nex = 1 / pex;
		return z.clone().setCartesian(+.5d * (pex + nex) * cos, -.5d * (pex - nex) * sin);
	}


	/**
	 * Returns <tt>tan(z) = -i&middot;(e<sup>+iz</sup>-e<sup>-iz</sup>)/(e<sup>+iz</sup>+e<sup>-iz</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <br />
	 *         <tt>(sin(&real;(z))&middot;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) + </tt>
	 *         <tt>cos(&real;(z))&middot;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)) / </tt><br />
	 *         <tt>(cos(&real;(z))&middot;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) - </tt>
	 *         <tt>sin(&real;(z))&middot;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T tan (final T z) throws NullPointerException {
		final float sin = (float) Math.sin(z.re()), cos = (float) Math.cos(z.re());
		final float pex = (float) Math.exp(z.im()), nex = 1 / pex;
		final T dividend = z.clone().setCartesian(+sin * (pex + nex), +cos * (pex - nex));
		final T divisor = z.clone().setCartesian(+cos * (pex + nex), -sin * (pex - nex));
		return dividend.div(divisor);
	}


	/**
	 * Returns <tt>tan(z) = -i&middot;(e<sup>+iz</sup>-e<sup>-iz</sup>)/(e<sup>+iz</sup>+e<sup>-iz</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <br />
	 *         <tt>(sin(&real;(z))&middot;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) + </tt>
	 *         <tt>cos(&real;(z))&middot;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>)) / </tt><br />
	 *         <tt>(cos(&real;(z))&middot;(e<sup>+&image;(z)</sup>+e<sup>-&image;(z)</sup>) - </tt>
	 *         <tt>sin(&real;(z))&middot;(e<sup>+&image;(z)</sup>-e<sup>-&image;(z)</sup>))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T tan (final T z) throws NullPointerException {
		final double sin = Math.sin(z.re()), cos = Math.cos(z.re());
		final double pex = Math.exp(z.im()), nex = 1 / pex;
		final T dividend = z.clone().setCartesian(+sin * (pex + nex), +cos * (pex - nex));
		final T divisor = z.clone().setCartesian(+cos * (pex + nex), -sin * (pex - nex));
		return dividend.div(divisor);
	}


	/**
	 * Returns the {@code arc sine} of {@code z}, i.e. <tt>asin(z) = </tt>
	 * <tt>-i&middot;log(i&middot;z &plusmn; (1-z<sup>2</sup>)<sup>&half;</sup>) = </tt><br />
	 * <tt>&half;&pi;-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>)) = &half;&pi;-acos(z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>&half;&pi;-acos(z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T asin (final T z) throws NullPointerException {
		return acos(z).neg().add((float) PI / 2);
	}


	/**
	 * Returns the {@code arc sine} of {@code z}, i.e. <tt>asin(z) = </tt>
	 * <tt>-i&middot;log(i&middot;z &plusmn; (1-z<sup>2</sup>)<sup>&half;</sup>)) = </tt><br />
	 * <tt>&half;&pi;-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>)) = &half;&pi;-acos(z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>&half;&pi;-acos(z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T asin (final T z) throws NullPointerException {
		return acos(z).neg().add(PI / 2);
	}


	/**
	 * Returns the {@code arc cosine} of {@code z}, i.e. <tt>acos(z) = </tt>
	 * <tt>-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>) = &half;&pi;-asin(z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T acos (final T z) throws NullPointerException {
		final T sqrt = z.clone().sq().sub(1).sqrt();
		if (z.re() > 0 & z.im() < 0 | z.re() < 0 & z.im() > 0) sqrt.neg();
		return log(sqrt.add(z)).idiv();
	}


	/**
	 * Returns the {@code arc cosine} of {@code z}, i.e. <tt>acos(z) = </tt>
	 * <tt>-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>) = = &half;&pi;-asin(z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T acos (final T z) throws NullPointerException {
		final T sqrt = z.clone().sq().sub(1).sqrt();
		if (z.re() > 0 & z.im() < 0 | z.re() < 0 & z.im() > 0) sqrt.neg();
		return log(sqrt.add(z)).idiv();
	}


	/**
	 * Returns the {@code arc tangent} of {@code z}, i.e. <tt>atan(z) = &half;i&middot;log((1-iz) / 1+iz)).
	 * &#64;param <T> the declaration type of the complex argument and result
	 * &#64;param z the complex operand
	 * @return the value <tt>&half;i&middot;log((1-iz) / 1+iz))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T atan (final T z) throws NullPointerException {
		return log(z.clone().idiv().add(1).div(z.clone().imul().add(1))).imul().mul(.5f);
	}


	/**
	 * Returns the {@code arc tangent} of {@code z}, i.e. <tt>atan(z) = &half;i&middot;log((1-iz) / 1+iz)).
	 * &#64;param <T> the declaration type of the complex argument and result
	 * &#64;param z the complex operand
	 * @return the value <tt>&half;i&middot;log((1-iz) / 1+iz))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T atan (final T z) throws NullPointerException {
		return log(z.clone().idiv().add(1).div(z.clone().imul().add(1))).imul().mul(.5d);
	}


	//***********************//
	// hyperbolic operations //
	//***********************//

	/**
	 * Returns the hyperbolic {@code sine} of {@code z}, i.e. <tt>sinh(z) = &half;(e<sup>z</sup>-e<sup>-z</sup>) = </tt>
	 * <tt>i&middot;sin(-i&middot;z) = -i&middot;sin(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;sin(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #sin(T)
	 */
	static public <T extends MutableSinglePrecision<T>> T sinh (final T z) throws NullPointerException {
		return sin(z.clone().imul()).idiv();
	}


	/**
	 * Returns the hyperbolic {@code sine} of {@code z}, i.e. <tt>sinh(z) = &half;(e<sup>z</sup>-e<sup>-z</sup>) = </tt>
	 * <tt>i&middot;sin(-i&middot;z) = -i&middot;sin(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;sin(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #sin(T)
	 */
	static public <T extends MutableDoublePrecision<T>> T sinh (final T z) throws NullPointerException {
		return sin(z.clone().imul()).idiv();
	}


	/**
	 * Returns the hyperbolic {@code cosine} of {@code z}, i.e.
	 * <tt>cosh(z) = &half;(e<sup>z</sup>+e<sup>-z</sup>) = cos(-i&middot;z) = cos(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>cos(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #cos(T)
	 */
	static public <T extends MutableSinglePrecision<T>> T cosh (final T z) throws NullPointerException {
		return cos(z.clone().imul());
	}


	/**
	 * Returns the hyperbolic {@code cosine} of {@code z}, i.e.
	 * <tt>cosh(z) = &half;(e<sup>z</sup>+e<sup>-z</sup>) = cos(-i&middot;z) = cos(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>cos(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #cos(T)
	 */
	static public <T extends MutableDoublePrecision<T>> T cosh (final T z) throws NullPointerException {
		return cos(z.clone().imul());
	}


	/**
	 * Returns the hyperbolic {@code tangent} of {@code z}, i.e. <tt>tanh(z) = -i&middot;sin(i&middot;z)/cos(i&middot;z) = </tt>
	 * <br />
	 * <tt>-i&middot;tan(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;tan(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #tan(T)
	 */
	static public <T extends MutableSinglePrecision<T>> T tanh (final T z) throws NullPointerException {
		return tan(z.clone().imul()).idiv();
	}


	/**
	 * Returns the hyperbolic {@code tangent} of {@code z}, i.e. <tt>tanh(z) = -i&middot;sin(i&middot;z)/cos(i&middot;z) = </tt>
	 * <br />
	 * <tt>-i&middot;tan(i&middot;z)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>-i&middot;tan(i&middot;z)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @see #tan(T)
	 */
	static public <T extends MutableDoublePrecision<T>> T tanh (final T z) throws NullPointerException {
		return tan(z.clone().imul()).idiv();
	}


	/**
	 * Returns the hyperbolic {@code area sine} of {@code z}, i.e.
	 * <tt>asinh(z) = log(z &plusmn; (z<sup>2</sup>+1)<sup>&half;</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>log(z + (z<sup>2</sup>+1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T asinh (final T z) throws NullPointerException {
		return log(z.clone().sq().add(1).sqrt().add(z));
	}


	/**
	 * Returns the hyperbolic {@code area sine} of {@code z}, i.e.
	 * <tt>asinh(z) = log(z &plusmn; (z<sup>2</sup>+1)<sup>&half;</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>log(z + (z<sup>2</sup>+1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T asinh (final T z) throws NullPointerException {
		return log(z.clone().sq().add(1).sqrt().add(z));
	}


	/**
	 * Returns the hyperbolic {@code area cosine} of {@code z}, i.e.
	 * <tt>acosh(z) = log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>) = </tt><br />
	 * <tt>log(z + (z+1)<sup>&half;</sup>&middot;(z-1)<sup>&half;</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>log(z + (z+1)<sup>&half;</sup>&middot;(z-1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T acosh (final T z) throws NullPointerException {
		return log(z.clone().add(1).sqrt().mul(z.clone().sub(1).sqrt()).add(z));
	}


	/**
	 * Returns the hyperbolic {@code area cosine} of {@code z}, i.e.
	 * <tt>acosh(z) = log(z &plusmn; (z<sup>2</sup>-1)<sup>&half;</sup>) = </tt><br />
	 * <tt>log(z + (z+1)<sup>&half;</sup>&middot;(z-1)<sup>&half;</sup>)</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>log(z + (z+1)<sup>&half;</sup>&middot;(z-1)<sup>&half;</sup>)</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T acosh (final T z) throws NullPointerException {
		return log(z.clone().add(1).sqrt().mul(z.clone().sub(1).sqrt()).add(z));
	}


	/**
	 * Returns the hyperbolic {@code area tangent} of {@code z}, i.e.<br />
	 * <tt>atanh(z) = &half;log((1+z) / (1-z))</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>&half;log((1+z) / (1-z))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableSinglePrecision<T>> T atanh (final T z) throws NullPointerException {
		return log(z.clone().add(1).div(z.clone().neg().add(1))).mul(.5f);
	}


	/**
	 * Returns the hyperbolic {@code area tangent} of {@code z}, i.e.<br />
	 * <tt>atanh(z) = &half;log((1+z) / (1-z))</tt>.
	 * @param <T> the declaration type of the complex argument and result
	 * @param z the complex operand
	 * @return the value <tt>&half;log((1+z) / (1-z))</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T extends MutableDoublePrecision<T>> T atanh (final T z) throws NullPointerException {
		return log(z.clone().add(1).div(z.clone().neg().add(1))).mul(.5d);
	}


	//*******************//
	// vector operations //
	//*******************//

	/**
	 * Returns a clone (deep copy) of the given complex vector.
	 * @param <T> the component type of the argument and result arrays
	 * @param vector the complex vector
	 * @return the clone
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws ArrayStoreException if the cloning of any vector element results in an incompatibe vector component type
	 */
	static public <T extends Complex<T>> T[] clone (final T[] vector) throws NullPointerException, ArrayStoreException {
		final T[] clone = vector.clone();
		for (int index = 0; index < vector.length; ++index) {
			final T value = clone[index];
			if (value != null) clone[index] = value.clone();
		}
		return clone;
	}


	/**
	 * Returns an image created from plotting the given function as a color wheel over the specified complex value range.
	 * @param function the complex function used to map complex values
	 * @param left the minimum real part of any complex operand value within the plot
	 * @param low the minimum imaginary part of any complex operand value within the plot
	 * @param right the maximum real part of any complex operand value within the plot
	 * @param high the maximum imaginary part of any complex operand value within the plot
	 * @param magnification the number of pixels to be created per distance unit.
	 * @param saturation the saturation of each pixel within range [0.0, 1.0]
	 * @return the plotted image
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public WritableImage plotColorWheel (final UnaryOperator<MutableSinglePrecision<?>> function, final float left, final float low, final float right, final float high, float magnification, final float saturation) throws NullPointerException {
		if (low > high | left > right | magnification <= 0 | saturation < 0 | saturation > 1) throw new IllegalArgumentException();

		final float resolution = 1 / magnification;
		final int width = (int) Math.round((right - left) * magnification);
		final int height = (int) Math.round((high - low) * magnification);
		final WritableImage image = new WritableImage(width, height);

		final FloatCartesianComplex z1 = new FloatCartesianComplex();
		for (int row = 0; row < height; ++row) {
			for (int column = 0; column < width; ++column) {
				final float re = left + column * resolution;
				final float im = high - row * resolution;

				z1.setCartesian(re, im);
				final MutableSinglePrecision<?> z2 = function.apply(z1);
				if (z2.isNaN()) z2.setCartesian(0, 0);

				// mapping: the origin is black, −1 is red, 1 is cyan, and a point at infinity is white:
				final double hue = 180 + (180 / PI) * z2.arg();
				final double brightness = 1 - Math.pow(.5d, z2.abs());

				final Color color = Color.hsb(hue, saturation, brightness);
				image.getPixelWriter().setColor(column, row, color);
			}
		}

		return image;
	}


	/**
	 * Returns an image created from plotting the given function as a color wheel over the specified complex value range.
	 * @param function the complex function used to map complex values
	 * @param left the minimum real part of any complex operand value within the plot
	 * @param low the minimum imaginary part of any complex operand value within the plot
	 * @param right the maximum real part of any complex operand value within the plot
	 * @param high the maximum imaginary part of any complex operand value within the plot
	 * @param magnification the number of pixels to be created per distance unit.
	 * @param saturation the saturation of each pixel within range [0.0, 1.0]
	 * @return the plotted image
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public WritableImage plotColorWheel (final UnaryOperator<MutableDoublePrecision<?>> function, final double left, final double low, final double right, final double high, double magnification, final float saturation) throws NullPointerException {
		if (low > high | left > right | magnification <= 0 | saturation < 0 | saturation > 1) throw new IllegalArgumentException();

		final double resolution = 1 / magnification;
		final int width = (int) Math.round((right - left) * magnification);
		final int height = (int) Math.round((high - low) * magnification);
		final WritableImage image = new WritableImage(width, height);

		final DoubleCartesianComplex z1 = new DoubleCartesianComplex();
		for (int row = 0; row < height; ++row) {
			for (int column = 0; column < width; ++column) {
				final double re = left + column * resolution;
				final double im = high - row * resolution;

				z1.setCartesian(re, im);
				final MutableDoublePrecision<?> z2 = function.apply(z1);
				if (z2.isNaN()) z2.setCartesian(0, 0);

				// mapping: the origin is black, +1 is red, -1 is cyan, and a point at infinity is white:
				final double hue = (180 / PI) * z2.arg();
				final double brightness = 1 - Math.pow(.5d, z2.abs());
				assert hue >= 0 & hue <= 1 & brightness >= 0 & brightness <= 1;

				final Color color = Color.hsb(hue, saturation, brightness);
				image.getPixelWriter().setColor(column, row, color);
			}
		}

		return image;
	}
}