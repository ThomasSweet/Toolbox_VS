package de.sb.toolbox.math;

import static java.lang.Double.doubleToRawLongBits;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.log;
import static java.lang.Math.scalb;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.System.arraycopy;
import static de.sb.toolbox.math.IntMath.floorLog2;
import static de.sb.toolbox.util.ArraySupport.swap;
import de.sb.toolbox.Copyright;


/**
 * This facade adds additional mathematical operations for {@code 64-bit} floating-point arithmetics.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public final class DoubleMath {

	/**
	 * An empty (and therefore both immutable and cacheable) vector.
	 */
	static public final double[] EMPTY = new double[0];

	/**
	 * <tt>&tau; = 2*&pi;</tt>. See the <a href="http://tauday.com/tau-manifesto">Tau-Manifesto</a> why TAU is the superior
	 * choice for the circle constant.
	 */
	static public final double TAU = 2 * PI;

	/**
	 * The binary square root, i.e. <tt>2<sup>&half;</sup></tt>.
	 */
	static public final double SQRT2 = 1.4142135623730951;

	/**
	 * The binary logarithm, i.e. <tt>log(2)</tt>.
	 */

	static public final double LOG2 = 0.6931471805599453094;
	/**
	 * The binary mask for double-precision floating-point sign access.
	 */

	static public final long MASK_DOUBLE_SIGN = 0x8000000000000000L;

	/**
	 * The binary mask for double-precision exponent access.
	 */
	static public final long MASK_DOUBLE_EXPONENT = 0x7ff0000000000000L;

	/**
	 * The binary mask for double-precision floating-point mantissa access.
	 */
	static public final long MASK_DOUBLE_MANTISSA = 0x000fffffffffffffL;


	/**
	 * Prevents external instantiation.
	 */
	private DoubleMath () {}


	//*******************//
	// scalar operations //
	//*******************//

	/**
	 * Clips the given value to the norm range <tt>[-1,+1]</tt>.
	 * @param value the value
	 * @return the value within norm range <tt>[-1,+1]</tt>
	 */
	static public double clip (final double value) {
		return Math.max(-1, Math.min(+1, value));
	}


	/**
	 * Returns the linear interpolation between two operands.
	 * @param ratio the ratio {@code r}, usually within range <tt>[0,1]</tt>
	 * @param left the left operand <tt>x<sub>L</sub></tt>
	 * @param right the right operand <tt>x<sub>R</sub></tt>
	 * @return the value <tt>(1-r)&middot;x<sub>L</sub> + r&middot;x<sub>R</sub></tt>
	 */
	static public double interpolate (final double ratio, final double left, final double right) throws IllegalArgumentException {
		return left + ratio * (right - left);
	}


	/**
	 * Returns the square of the given operand {@code x}.
	 * @param x the operand value
	 * @return the value <tt>x<sup>2</sup></tt>
	 */
	static public double sq (final double x) {
		return x * x;
	}


	/**
	 * Returns the cube of the given operand {@code x}.
	 * @param x the operand value
	 * @return the value <tt>x<sup>3</sup></tt>
	 */
	static public double cb (final double x) {
		return x * x * x;
	}


	/**
	 * Returns the {@code Euclidean modulo} of the given dividend and divisor, based on the {@code truncated modulo} returned by
	 * the {@code %} operator. Note that both the truncated modulo and the divisor absolute are precalculated in order to profit
	 * from both supersymmetric computation and branch-free conditional assignment (<i>CMOV</i>).
	 * @param x the dividend value
	 * @param y the divisor value
	 * @return the value <tt>x mod<sub>e</sub> y</tt> within range {@code [0,|y|[}, or {@code NaN} if the given divisor is zero
	 */
	static public double mod (final double x, final double y) {
		final double mod = x % y, abs = Math.abs(y);
		return x >= 0 ? mod : abs + mod;
	}


	/**
	 * Returns the {@code binary logarithm} of the given operand.
	 * @param operand the operand
	 * @return the value <tt>ceil(log<sub>2</sub>(x))</tt>, with {@code -1} representing {@code negative infinity}
	 * @throws ArithmeticException if the given operand is strictly negative
	 */
	static public double log2 (final double x) throws ArithmeticException {
		return log(x) * (1 / LOG2);
	}



	//********************//
	// complex operations //
	//********************//

	/**
	 * Returns the {@code real} part <tt>&real;(z)</tt> of the given complex number <tt>z = x*e<sup>i&middot;y</sup></tt>.
	 * @param x the absolute value <tt>|z|</tt>
	 * @param y the argument angle <tt>&phi;(z)</tt>
	 * @return the real value <tt>&real;(z) = x&middot;cos(y)</tt>
	 */
	static public double re (final double x, final double y) {
		if (y == +0) return +x;
		if (y == -PI) return -x;
		if (y == +PI / 2 | y == -PI / 2) return 0;
		return x * cos(y);
	}


	/**
	 * Returns the {@code imaginary} part <tt>&image;(z)</tt> of the given complex number <tt>z = x*e<sup>i&middot;y</sup></tt>.
	 * @param x the absolute value <tt>|z|</tt>
	 * @param y the argument angle <tt>&phi;(z)</tt>
	 * @return the imaginary value <tt>&image;(z) = x&middot;sin(y)</tt>
	 */
	static public double im (final double x, final double y) {
		if (y == +PI / 2) return +x;
		if (y == -PI / 2) return -x;
		if (y == +0 | y == -PI) return 0;
		return x * sin(y);
	}


	/**
	 * Returns the <tt>absolute value |z|</tt> of the given complex number <tt>z = x+i&middot;y</tt>.
	 * @param x the real part <tt>&real;(z)</tt>
	 * @param y the imaginary part <tt>&image;(z)</tt>
	 * @return the absolute value <tt>|z| = (&real;(z)<sup>2</sup> +
	 *         &image;(z)<sup>2</sup>)<sup>&half;</sup></tt>
	 */
	static public double abs (final double x, final double y) {
		if (x == 0) return Math.abs(y);
		if (y == 0) return Math.abs(x);
		return sqrt(x * x + y * y);
	}


	/**
	 * Returns the <tt>argument angle &phi;(z)</tt> of the given complex number {@code z = x+iy}.
	 * @param x the real part <tt>&real;(z)</tt>
	 * @param y the imaginary part <tt>&image;(z)</tt>
	 * @return the argument angle <tt> &phi; = atan2(&image;(z), &real;(z))</tt>
	 */
	static public double arg (final double x, final double y) {
		if (y == 0) return (doubleToRawLongBits(x) & MASK_DOUBLE_SIGN) == 0 ? TAU * +0 / 8 : TAU * -4 / 8;
		if (x == 0) return y > 0 ? TAU * +2 / 8 : TAU * -2 / 8;
		if (x == +y) return x > 0 ? TAU * +1 / 8 : TAU * -3 / 8;
		if (x == -y) return x > 0 ? TAU * -1 / 8 : TAU * +3 / 8;
		return atan2(y, x);
	}


	//*******************//
	// vector operations //
	//*******************//

	/**
	 * Returns the linear interpolation between two vector operands. Note that the left operand carries the result, and that the
	 * right operands is also modified unless {@code r=1}.
	 * @param ratio the ratio {@code r}, usually within range <tt>[0,1]</tt>
	 * @param left the left vector operand <tt>X<sub>L</sub></tt>
	 * @param right the right vector operand <tt>X<sub>R</sub></tt>
	 * @return the value <tt>(1-r)&middot;X<sub>L</sub> + r&middot;X<sub>R</sub></tt>
	 */
	static public double[] interpolate (final double ratio, final double[] left, final double[] right) throws IllegalArgumentException {
		return add(mul(left, 1 - ratio), mul(right, ratio));
	}


	/**
	 * Returns the {@code smallest} coordinate of the given vector, i.e. the one with a value closest to negative infinity. If
	 * the vector has zero length, the result is zero.
	 * @param vector the operand vector
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public double min (final double... vector) throws NullPointerException {
		if (vector.length == 0) return 0;

		double min = vector[0];
		for (int index = 1; index < vector.length; ++index) {
			min = Math.min(min, vector[index]);
		}
		return min;
	}


	/**
	 * Returns the {@code largest} coordinate of the given vector, i.e. the one with a value closest to positive infinity. If
	 * the vector has zero length, the result is zero.
	 * @param vector the operand vector
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public double max (final double... vector) throws NullPointerException {
		if (vector.length == 0) return 0;

		double max = vector[0];
		for (int index = 1; index < vector.length; ++index) {
			max = Math.max(max, vector[index]);
		}
		return max;
	}


	/**
	 * Returns the {@code absolute value} (the Euclidean norm) of the given vector.
	 * @param vector the operand vector
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public double abs (final double... vector) throws NullPointerException {
		return sqrt(dot(vector, vector));
	}


	/**
	 * Returns the {@code argument angle} between the given vectors.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the angle in radians within range <tt>[-&pi;,&pi;[</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given arguments do not share the same length
	 */
	static public double arg (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		final double dotProduct = dot(left, right);
		final double absProduct = abs(left) * abs(right);
		final double arg = acos(dotProduct / absProduct);
		return arg < (PI / 2) ? arg : arg - PI;
	}


	/**
	 * {@code Assigns} the state of the right vector to the left vector, and returns the latter.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the modified left operand
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given arguments do not share the same length
	 */
	static public double[] assign (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length) throw new IllegalArgumentException();
		arraycopy(right, 0, left, 0, left.length);
		return left;
	}


	/**
	 * {@code Adds} the right vector to the left vector, assigns the result to the latter, and returns it.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the modified left operand
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given arguments do not share the same length
	 */
	static public double[] add (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length) throw new IllegalArgumentException();
		for (int index = 0; index < left.length; ++index)
			left[index] += right[index];
		return left;
	}


	/**
	 * {@code Subtracts} the right vector from the left vector, assigns the result to the latter, and returns it.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the modified left operand
	 * @throws NullPointerException if any of the given vectors is {@code null}
	 * @throws IllegalArgumentException if the two vectors do not share the same length
	 */
	static public double[] sub (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length) throw new IllegalArgumentException();
		for (int index = 0; index < left.length; ++index)
			left[index] -= right[index];
		return left;
	}


	/**
	 * {@code Scalar-multiplies} the given given value with the given vector, assigns the result to the latter, and returns it.
	 * @param left the left vector operand
	 * @param right the right scalar operand
	 * @return the modified left operand
	 * @throws NullPointerException if the given vector is {@code null}
	 */
	static public double[] mul (final double[] left, final double right) throws NullPointerException {
		for (int index = 0; index < left.length; ++index)
			left[index] *= right;
		return left;
	}


	/**
	 * Returns the result of {@code dot-multiplying} the left vector with the right vector.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the scalar {@code dot-product} of both vectors
	 * @throws NullPointerException if any of the given vectors is {@code null}
	 * @throws IllegalArgumentException if the two vectors do not share the same length
	 */
	static public double dot (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length) throw new IllegalArgumentException();

		double result = 0;
		for (int index = 0; index < left.length; ++index) {
			result += left[index] * right[index];
		}
		return result;
	}


	/**
	 * Assigns the result of {@code cross-multiplying} the left vector with the right vector to the former, and returns it.
	 * @param left the left vector operand
	 * @param right the right vector operand
	 * @return the modified left operand
	 * @throws NullPointerException if any of the given vectors is {@code null}
	 * @throws IllegalArgumentException if the two vectors do not share the same length, or the length is neither {@code 3} nor
	 *         {@code 7}
	 */
	static public double[] cross (final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length) throw new IllegalArgumentException();
		switch (left.length) {
			case 0: {
				break;
			}
			case 1: {// v[i] = 0
				left[0] = 0;
				break;
			}
			case 3: {// v[i] = [(i+1) % 3] x [(i+2) % 3]
				final double l0 = left[0], l1 = left[1], l2 = left[2];
				final double r0 = right[0], r1 = right[1], r2 = right[2];
				left[0] = l1 * r2 - l2 * r1;
				left[1] = l2 * r0 - l0 * r2;
				left[2] = l0 * r1 - l1 * r0;
				break;
			}
			case 7: {// v[i] = [(i+1) % 7] x [(i+3) % 7] + [(i+2) % 7] x [(i+6) % 7] + [(i+4) % 7] x [(i+5) % 7]
				final double l0 = left[0], l1 = left[1], l2 = left[2], l3 = left[3], l4 = left[4], l5 = left[5], l6 = left[6];
				final double r0 = right[0], r1 = right[1], r2 = right[2], r3 = right[3], r4 = right[4], r5 = right[5], r6 = right[6];
				left[0] = l1 * r3 - l3 * r1 + l2 * r6 - l6 * r2 + l4 * r5 - l5 * r4;
				left[1] = l2 * r4 - l4 * r2 + l3 * r0 - l0 * r3 + l5 * r6 - l6 * r5;
				left[2] = l3 * r5 - l5 * r3 + l4 * r1 - l1 * r4 + l6 * r0 - l0 * r6;
				left[3] = l4 * r6 - l6 * r4 + l5 * r2 - l2 * r5 + l0 * r1 - l1 * r0;
				left[4] = l5 * r0 - l0 * r5 + l6 * r3 - l3 * r6 + l1 * r2 - l2 * r1;
				left[5] = l6 * r1 - l1 * r6 + l0 * r4 - l4 * r0 + l2 * r3 - l3 * r2;
				left[6] = l0 * r2 - l2 * r0 + l1 * r5 - l5 * r1 + l3 * r4 - l4 * r3;
			}
			default: {
				throw new IllegalArgumentException();
			}
		}
		return left;
	}


	/**
	 * Assigns the respective values <tt>v[j]'=v[j]+w[k]</tt> and <tt>w[k]'=v[j]-w[k]</tt>.
	 * @param left the left vector operand <tt>v</tt>
	 * @param leftIndex the vector index <tt>j</tt>
	 * @param right the right vector operand <tt>w</tt>
	 * @param rightIndex the vector index <tt>k</tt>
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void mux (final double[] left, final int leftIndex, final double[] right, final int rightIndex) throws NullPointerException, IndexOutOfBoundsException {
		final double leftValue = left[leftIndex], rightValue = right[rightIndex];
		left[leftIndex] = leftValue + rightValue;
		right[rightIndex] = leftValue - rightValue;
	}


	//**********************//
	// transform operations //
	//**********************//

	/**
	 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of <tt>&half;N</tt> <i>braided</i> complex
	 * numbers.
	 * @param inverse whether or not an {@code inverse} fourier transform shall be performed
	 * @param separate whether or not <i>channel separation</i> shall be performed
	 * @param vector an array of <tt>&half;N</tt> braided complex numbers in Cartesian form, alternating even indexed real parts
	 *        with odd indexed imaginary ones
	 * @throws NullPointerException if the given vector is {@code null}
	 * @throws IllegalArgumentException if the given vector's length is not an even power of two
	 */
	static public void fft (final boolean inverse, final boolean separate, final double[] vector) throws NullPointerException, IllegalArgumentException {
		if (inverse) {
			for (int index = 1; index < vector.length; index += 2)
				vector[index] *= -1;
		} else if (separate) {
			entangle(inverse, vector);
		}

		fft(vector);

		if (inverse) {
			final double norm = scalb(2d, -floorLog2(vector.length));
			for (int index = 0; index < vector.length; index += 2) {
				vector[index] *= +norm;
				vector[index + 1] *= -norm;
			}
			if (separate) entangle(inverse, vector);
		}
	}


	/**
	 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of <tt>&half;N</tt> <i>braided</i> complex
	 * numbers. Note that an <i>inverse</i> transform can be performed in two ways:
	 * <ul>
	 * <li>by conjugating both the argument and result of this operation, and additionally norming the result by {@code 2/N}.
	 * </li>
	 * <li>by swapping the real and imaginary parts of both the argument and the result.</i>
	 * </ul>
	 * <i>Implementation notice</i>: This is an extremely optimized variant of the {@code Cooley-Tukey} algorithm, based on
	 * <i>perfect shuffle</i> and <i>sine/cosine</i> function tables for caching. Note that the JIT-compiler will inline the
	 * lookups within these tables (and any other static methods called). Logarithm-based iteration additionally eliminates
	 * expensive divisions, leaving only cheap addition, subtraction, multiplication and shift operations to be performed. The
	 * number of index variables is reduced to increase the chance that the JIT-compiler can keep them in registers (to avoid
	 * expensive cache or worse memory lookups), while derived indices are cheaply recalculated on demand taking advantage of
	 * modern super-scalar processor designs.
	 * @param vector an array of <tt>&half;N</tt> braided complex numbers in Cartesian form, alternating even indexed real parts
	 *        with odd indexed imaginary ones
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument's length is not an even power of two
	 */
	static public void fft (final double[] vector) throws NullPointerException, IllegalArgumentException {
		final int magnitude = floorLog2(vector.length) - 1;
		if (vector.length != 2 << magnitude) throw new IllegalArgumentException();
		if (magnitude == 0) return;

		for (final FunctionTables.SwapEntry entry : FunctionTables.getPerfectShuffleTable(magnitude)) {
			final int left = 2 * entry.getLeft(), right = 2 * entry.getRight();
			swap(vector, left, vector, right);
			swap(vector, left + 1, vector, right + 1);
		}

		final FunctionTables.Trigonometric trigonometricTable = FunctionTables.getTrigonometricTable(magnitude);
		for (int depth = 0; depth < magnitude; ++depth) {
			for (int offset = 0; offset < (1 << depth); ++offset) {
				final int angleIndex = offset << (magnitude - depth - 1);
				final double cos = trigonometricTable.cos(angleIndex);
				final double sin = trigonometricTable.sin(angleIndex);

				for (int left = 2 * offset, right = left + (2 << depth); left < (2 << magnitude); left += 4 << depth, right += 4 << depth) {
					double re = vector[right], im = vector[right + 1];
					final double taoRe = cos * re - sin * im, taoIm = cos * im + sin * re;
					re = vector[left];
					im = vector[left + 1];
					vector[right] = re - taoRe;
					vector[right + 1] = im - taoIm;
					vector[left] = re + taoRe;
					vector[left + 1] = im + taoIm;
				}
			}
		}
	}


	/**
	 * If inverse is {@code true}, the given <i>natural spectrum</i> is detangled into a <i>two-channel spectrum</i> using the
	 * following operations for each corresponding complex spectrum entry (r,l,t &isin; &#x2102;): <br />
	 * <tt>forEach(l,r): r' = i&middot;(r<sup>*</sup>-l); l' = r<sup>*</sup>+l</tt>
	 * <p>
	 * If inverse is {@code false}, the given two-channel spectrum is entangled into a <i>natural spectrum</i> using the
	 * following operations (again r,l,t &isin; &#x2102;):<br />
	 * <tt>r,l,t &isin; &#x2102;: r' = &half;(l-i&middot;r)<sup>*</sup>; l' = &half;(l+i&middot;r)</tt>
	 * </p>
	 * <p>
	 * Additionally, the array elements at index <tt>1</tt> and <tt>&frac12;N</tt> are swapped, which implies that after
	 * detangling the
	 * <ul>
	 * <li>spectrum entries <tt>0</tt> and <tt>1</tt> both belong to the left channel, representing the signed amplitudes of the
	 * frequencies <tt>f<sub>0</sub></tt> and <tt>f<sub>Nyquist</sub></tt></li>
	 * <li>spectrum entries <tt>&frac12;N</tt> and <tt>&frac12;N+1</tt> both belong to the right channel, again representing the
	 * signed amplitudes of the frequencies <tt>f<sub>0</sub></tt> and <tt>f<sub>Nyquist</sub></tt></li>
	 * </ul>
	 * Detangling of a natural spectrum created by an iFFT operation is required in order to cleanly separate it's left and
	 * right halves; in it's detangled state, a spectrum's left half is solely influenced by the real parts of the values that
	 * went into the iFFT, and the right half is solely influenced by the corresponding imaginary parts.
	 * </p>
	 * @param spectrum an array of <tt>&half;N</tt> braided complex numbers in Cartesian form, alternating even indexed real
	 *        parts with odd indexed imaginary ones
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument's length is odd
	 */
	static private void entangle (final boolean inverse, final double[] spectrum) {
		if ((spectrum.length & 1) == 1) throw new IllegalArgumentException();

		final int half = spectrum.length / 2;
		final double swap = spectrum[half];
		spectrum[half] = spectrum[1];
		spectrum[1] = swap;

		if (inverse) {
			for (int leftIndex = 2, rightIndex = spectrum.length - 2; leftIndex < half; leftIndex += 2, rightIndex -= 2) {
				final double leftRe = spectrum[leftIndex], leftIm = spectrum[leftIndex + 1];
				final double rightRe = spectrum[rightIndex], rightIm = spectrum[rightIndex + 1];
				spectrum[leftIndex + 0] = +leftRe + rightRe;
				spectrum[leftIndex + 1] = +leftIm - rightIm;
				spectrum[rightIndex + 0] = +leftIm + rightIm;
				spectrum[rightIndex + 1] = -leftRe + rightRe;
			}
		} else {
			for (int leftIndex = 2, rightIndex = spectrum.length - 2; leftIndex < half; leftIndex += 2, rightIndex -= 2) {
				final double leftRe = spectrum[leftIndex], leftIm = spectrum[leftIndex + 1];
				final double rightRe = spectrum[rightIndex], rightIm = spectrum[rightIndex + 1];
				spectrum[leftIndex + 0] = (+leftRe - rightIm) / 2;
				spectrum[leftIndex + 1] = (+leftIm + rightRe) / 2;
				spectrum[rightIndex + 0] = (+leftRe + rightIm) / 2;
				spectrum[rightIndex + 1] = (-leftIm + rightRe) / 2;
			}
		}
	}
}