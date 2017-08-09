package de.sb.toolbox.math;

import static de.sb.toolbox.math.FloatMath.MASK_FLOAT_EXPONENT;
import static de.sb.toolbox.math.FloatMath.MASK_FLOAT_MANTISSA;
import static de.sb.toolbox.math.FloatMath.MASK_FLOAT_SIGN;
import static de.sb.toolbox.math.IntMath.floorLog2;
import static de.sb.toolbox.util.ArraySupport.swap;
import static java.lang.Float.floatToRawIntBits;
import static java.lang.Float.intBitsToFloat;
import static java.lang.Math.PI;
import static java.lang.Math.scalb;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.math.Complex.AbstractSinglePrecision;
import de.sb.toolbox.math.Complex.MutableSinglePrecision;
import de.sb.toolbox.math.FunctionTables.SwapEntry;


/**
 * This adapter/facade adds mathematical operations for complex values packed into 64bit {@code long} types. The resulting
 * values are bit-fields of type {@code long} in {@code Cartesian} form, i.e. {@code z = x + i&middot;y}, and contain the
 * {@code real part} as a 32bit floating-point number ( {@code float}) within it's lower halfword, and the
 * {@code imaginary part} as a 32bit floating-point number in it's upper halfword.<br />
 * <br />
 * The main advantage of this design is that complex values are represented as primitive types, not objects. Therefore, they
 * will be initialized automatically to zero instead of requiring explicit initialization. Standard cloning works both nicely
 * and fast, and there are no issues with intermediate objects. Additionally, there are no indexing issues as each complex value
 * will only claim one index within an array (instead of two if real and imaginary parts are handled separately).<br />
 * <br />
 * The downside of this design is that the real and imaginary parts need to be unpacked before most operations, and complex
 * results need to be packed; this is slower than array access as it involves bit shifting and masking, and additionall limits
 * accuracy to 32bits. Note that the Java "long" operators are NOT suitable for anything sensible with regard to packed values!
 * Therefore, it is usually best to handle them as "black-box" bit-fields, and use the operations provided herein to manipulate
 * them.<br />
 * <br />
 * Note that accuracy of all basic operators (add, sub, mul, div) is within one unit in the last place (1ulp), instead of the
 * loss of up to half the precision that is usually associated with floating-point multiplication and division operations! This
 * is achieved by performing multiplication and division operations in double precision internally, for the price of slightly
 * reduced speed. Finally, note that this class supports the Tau manifest in replacing any occurrences of the half-circle
 * constant <tt>Pi (&pi;)</tt> with the full-circle constant <tt>Tau (&tau;)</tt>, defined as <tt>&tau; = 2&middot;π</tt>.
 * Finally, note that this documentation uses operator symbols for all four complex coordinates:
 * <ul>
 * <li><tt>&real;(z)</tt> is the real part of z in Cartesian coordinates</li>
 * <li><tt>&image;(z)</tt> is the imaginary part of z in Cartesian coordinates</li>
 * <li><tt>|z| = abs(z)</tt> is the absolute value of z in Polar coordinates.</li>
 * <li><tt>&phi;(z) = arg(z)</tt> is the argument angle of z in Polar coordinates.</li>
 * </ul>
 * .
 */
@Copyright(year = 2010, holders = "Sascha Baumeister")
public final class PrimitiveCartesian extends AbstractSinglePrecision<PrimitiveCartesian>implements MutableSinglePrecision<PrimitiveCartesian> {
	static private final long serialVersionUID = 1;
	static private final long MASK_HALFWORD = 0xffffffffL;
	static public final long ZERO = fromCartesian(0, 0);
	static public final long ONE = fromCartesian(1, 0);
	static public final long E = fromCartesian((float) Math.E, 0);
	static public final long TAU = fromCartesian((float) (2 * Math.PI), 0);
	static public final long I = fromCartesian(0, 1);

	private final long[] vector;
	private final int index;


	//****************************//
	// complex adapter operations //
	//****************************//

	/**
	 * Creates a new instance based on the value zero.
	 */
	public PrimitiveCartesian () {
		this(new long[] { 0 }, 0);
	}


	/**
	 * Creates a new instance based on the given value.
	 * @param value the complex value in {@code primitive-cartesian} representation
	 */
	public PrimitiveCartesian (final long value) {
		this(new long[] { value }, 0);
	}


	/**
	 * Creates a new instance based on the given vector element.
	 * @param vector the complex vector in {@code primitive-cartesian} representation
	 * @param index the element index
	 * @throws NullPointerException if the given vector is {@code null}
	 * @throws IllegalArgumentException if the given index is out of range
	 */
	public PrimitiveCartesian (final long[] vector, final int index) throws NullPointerException, IllegalArgumentException {
		if (index < 0 | index >= vector.length) throw new IllegalArgumentException();
		this.vector = vector;
		this.index = index;
	}


	/**
	 * {@inheritDoc} Note that this implementation creates a copy based on a single value array in order to decouple the clone
	 * from the receiver's array reference.
	 */
	public PrimitiveCartesian clone () {
		return new PrimitiveCartesian(this.get());
	}


	/**
	 * Returns the {@code primitive-cartesian} representation of this complex number {@code z}
	 */
	public long get () {
		return this.vector[this.index];
	}


	/**
	 * Performs the state assignment <tt>z<sub>1</sub> = z<sub>2</sub></tt>.
	 * @param z the complex operand <tt>z<sub>2</sub></tt> in {@code primitive-cartesian} representation
	 * @return the modified receiver
	 */
	public PrimitiveCartesian set (final long z) {
		this.vector[this.index] = z;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isReal () {
		return isReal(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isImaginary () {
		return isImaginary(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isZero () {
		return isZero(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isInfinite () {
		return isInfinite(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isNaN () {
		return isNaN(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public float re () {
		return re(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public float im () {
		return im(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public float abs () {
		return abs(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public float arg () {
		return arg(this.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian conj () {
		return this.set(conj(this.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian neg () {
		return this.set(neg(this.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian imul () {
		return this.set(imul(this.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian idiv () {
		return this.set(idiv(this.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian inv () {
		return this.set(inv(this.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian setCartesian (final float re, final float im) {
		return this.set(fromCartesian(re, im));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian setPolar (final float abs, final float arg) {
		return this.set(fromPolar(abs, arg));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian add (final float value) {
		return this.set(add(this.get(), value));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian sub (final float value) {
		return this.set(sub(this.get(), value));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian mul (final float value) {
		return this.set(mul(this.get(), value));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian div (final float value) {
		return this.set(div(this.get(), value));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian set (final PrimitiveCartesian value) throws NullPointerException {
		return this.set(value.get());
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian add (final PrimitiveCartesian value) throws NullPointerException {
		return this.set(add(this.get(), value.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian sub (final PrimitiveCartesian value) throws NullPointerException {
		return this.set(sub(this.get(), value.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian mul (final PrimitiveCartesian value) throws NullPointerException {
		return this.set(mul(this.get(), value.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian div (final PrimitiveCartesian value) throws NullPointerException {
		return this.set(div(this.get(), value.get()));
	}


	/**
	 * {@inheritDoc}
	 */
	public void mux (final PrimitiveCartesian value) throws NullPointerException {
		final long z1 = this.get(), z2 = value.get();
		final float re1 = re(z1), im1 = im(z1), re2 = re(z2), im2 = im(z2);
		this.setCartesian(re1 + re2, im1 + im2);
		value.setCartesian(re1 - re2, im1 - im2);
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian sq () {
		this.vector[this.index] = sq(this.get());
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian sqrt () {
		this.vector[this.index] = sqrt(this.get());
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian cb () {
		this.vector[this.index] = cb(this.get());
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public PrimitiveCartesian cbrt () {
		this.vector[this.index] = cbrt(this.get());
		return this;
	}


	//****************************//
	// primitive value operations //
	//****************************//

	/**
	 * Returns the {@code primitive-cartesian} representation of the given complex number {@code z}. The resulting value is a
	 * bit-field of primitive type {@code long}, and contains the real part <tt>&real;(z)</tt> as a 32bit floating-point number
	 * within it's lower halfword, and the imaginary part <tt>&image;(z)</tt> as a 32bit floating-point number in it's upper
	 * halfword.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the {@code primitive-cartesian} representation of <tt>z=&real;(z)+i&middot;&image;(z)</tt>
	 */
	static public long fromComplex (final SinglePrecision<?> z) {
		if (z instanceof PrimitiveCartesian) {
			final PrimitiveCartesian other = (PrimitiveCartesian) z;
			return other.vector[other.index];
		}
		return fromCartesian(z.re(), z.im());
	}


	/**
	 * Returns the {@code primitive-cartesian} representation of the complex number <tt>z=&real;(z)+i&middot;&image;</tt>. The
	 * resulting value is a bit-field of primitive type {@code long}, and contains the real part <tt>&real;(z)</tt> as a 32bit
	 * floating-point number within it's lower halfword, and the imaginary part <tt>&image;(z)</tt> as a 32bit floating-point
	 * number in it's upper halfword.
	 * @param re the real part <tt>&real;(z)</tt>
	 * @param im the imaginary part <tt>&image;(z)</tt>
	 * @return the {@code primitive-cartesian} representation of <tt>z=&real;(z)+i&middot;&image;(z)</tt>
	 */
	static public long fromCartesian (final float re, final float im) {
		return ((long) floatToRawIntBits(im) << Float.SIZE) | (floatToRawIntBits(re) & MASK_HALFWORD);
	}


	/**
	 * Returns the {@code primitive-cartesian} representation of the complex number
	 * <tt>z = |z|&middot;e<sup>i&middot;&phi;(z)</sup></tt> , using Euler's formula for the conversion. The resulting value is
	 * a bit-field of primitive type {@code long}, and contains the real part <tt>&real;(z) = |z|&middot;cos(&phi;(z))</tt> as a
	 * 32bit floating-point number within it's lower halfword, and the imaginary part
	 * <tt>&image;(z) = |z|&middot;sin(&phi;(z))</tt> as a 32bit floating-point number in it's upper halfword.
	 * @param abs the absolute value <tt>|z|</tt>
	 * @param arg the argument angle <tt>&phi;(z)</tt>
	 * @return the {@code primitive-cartesian} representation of <tt>z = |z|&middot;e<sup>i&middot;&phi;(z)</sup></tt>
	 */
	static public long fromPolar (final float abs, final float arg) {
		return fromCartesian(FloatMath.re(abs, arg), FloatMath.im(abs, arg));
	}


	/**
	 * Returns the {@code Cartesian} text representation of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the text representation of <tt>z = &real;(z)+i&middot;&image;(z)</tt>
	 */
	static public String toCartesianString (final long z) {
		if (isNaN(z)) return "NaN";
		if (isInfinite(z)) return "INFINITY";
		return String.format(CARTESIAN_FORMAT, re(z), im(z));
	}


	/**
	 * Returns the {@code Polar} text representation of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the text representation of <tt>z = |z|&middot;e<sup>i&middot;&phi;(z)</sup></tt>, with an argument angle
	 *         normalized to multiples of <tt>&pi;</tt>
	 */
	static public String toPolarString (final long z) {
		if (isNaN(z)) return "NaN";
		if (isInfinite(z)) return "INFINITY";
		return String.format(POLAR_FORMAT, abs(z), arg(z) * (1 / PI));
	}


	/**
	 * Returns {@code true} if the given operand {@code z} represents a {@code real} number including zero, {@code false}
	 * otherwise. The imaginary part &image;(z) of real numbers is {@code 0}, and their argument angle
	 * <tt>((&phi;(z)+&pi;) mod<sub>e</sub> &tau;)-&pi;</tt> is {@code 0} or <tt>&plusmn;&pi;</tt> (except for zero).
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return whether or not the given operand represents a real number.
	 */
	static public boolean isReal (final long z) {
		return (z & ((MASK_HALFWORD ^ MASK_FLOAT_SIGN) << Float.SIZE)) == 0;
	}


	/**
	 * Returns {@code true} if the given operand {@code z} represents an {@code imaginary} number including zero, {@code false}
	 * otherwise. The real part &real;(z) of imaginary numbers is {@code 0}, and their argument angle
	 * <tt>((&phi;(z)+&pi;) mod<sub>e</sub> &tau;)-&pi;</tt> is <tt>&plusmn;&half;&pi;</tt> (except for zero).
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return whether or not the given operand represents an imaginary number.
	 */
	static public boolean isImaginary (final long z) {
		return (z & (MASK_HALFWORD ^ MASK_FLOAT_SIGN)) == 0;
	}


	/**
	 * Returns {@code true} if the given operand represents {@code zero}, {@code false} otherwise.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return whether or not the given operand represents zero.
	 */
	static public boolean isZero (final long z) {
		return (z & ~(MASK_FLOAT_SIGN | (MASK_FLOAT_SIGN << Float.SIZE))) == 0;
	}


	/**
	 * Returns {@code true} if the given operand represents {@code infinity}, {@code false} otherwise.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return whether or not the given operand represents zero.
	 */
	static public boolean isInfinite (final long z) {
		final long y = z >>> Float.SIZE;
		return (z & (MASK_FLOAT_EXPONENT | MASK_FLOAT_MANTISSA)) == MASK_FLOAT_EXPONENT | (y & (MASK_FLOAT_EXPONENT | MASK_FLOAT_MANTISSA)) == MASK_FLOAT_EXPONENT;
	}


	/**
	 * Returns {@code true} if the given operand represents {@code NaN}, {@code false} otherwise.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return whether or not the given operand does not represent a number.
	 */
	static public boolean isNaN (final long z) {
		final long y = z >>> Float.SIZE;
		return ((z & MASK_FLOAT_EXPONENT) == MASK_FLOAT_EXPONENT & (z & MASK_FLOAT_MANTISSA) != 0) | ((y & MASK_FLOAT_EXPONENT) == MASK_FLOAT_EXPONENT & (y & MASK_FLOAT_MANTISSA) != 0);
	}


	/**
	 * Returns the {@code real} part <tt>&real;(z)</tt> of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the real value <tt>&real;(z)</tt>
	 */
	static public float re (final long z) {
		return intBitsToFloat((int) z);
	}


	/**
	 * Returns the {@code imaginary} part <tt>&image;(z)</tt> of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the imaginary value <tt>&image;(z)</tt>
	 */
	static public float im (final long z) {
		return intBitsToFloat((int) (z >>> Float.SIZE));
	}


	/**
	 * Returns the absolute value <tt>|z|</tt> of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the absolute value <tt>|z| = (&real;(z)<sup>2</sup>+&image;(z)<sup>2</sup>)<sup>&half;</sup></tt>
	 */
	static public float abs (final long z) {
		return FloatMath.abs(re(z), im(z));
	}


	/**
	 * Returns the argument angle <tt>&phi;(z)</tt> of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the argument angle <tt>&phi;(z) = atan2(&image;(z),&real;(z))</tt>
	 */
	static public float arg (final long z) {
		return FloatMath.arg(re(z), im(z));
	}


	/**
	 * Returns the {@code conjugate} value <tt>z' = <span style="text-decoration: overline;">z</span></tt> of the given operand
	 * {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value <tt><span style="text-decoration: overline;">z</span> = </tt> <tt>+&real;(z)-i&middot;&image;(z)</tt>
	 *         in {@code primitive-cartesian} representation
	 */
	static public long conj (final long z) {
		return z ^ (MASK_FLOAT_SIGN << Float.SIZE);
	}


	/**
	 * Returns the {@code negative} value <tt>z'= -z</tt> of the given operand {@code z}. This effectively reflects the number
	 * on the origin.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value <tt>-&real;(z) - i&middot;&image;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long neg (final long z) {
		return z ^ (MASK_FLOAT_SIGN | (MASK_FLOAT_SIGN << Float.SIZE));
	}


	/**
	 * Returns the {@code imaginary product} <tt>z'= +i&middot;z</tt> of the given operand {@code z}. This effectively turns the
	 * operand by 90&deg; counter-clockwise around the origin.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value of <tt>+i&middot;z = -&image;(z)+i&middot;&real;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long imul (final long z) {
		return (z << Float.SIZE) | ((z >>> Float.SIZE) ^ MASK_FLOAT_SIGN);
	}


	/**
	 * Returns the {@code imaginary quotient} <tt>z'= z/i = -i&middot;z</tt> of the given operand {@code z}. This effectively
	 * turns the operand by 90&deg; clockwise around the origin.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value of <tt>-i&middot;z = +&image;(z)-i&middot;&real;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long idiv (final long z) {
		return ((z ^ MASK_FLOAT_SIGN) << Float.SIZE) | (z >>> Float.SIZE);
	}


	/**
	 * Returns the {@code reciprocal} <tt>z'= 1/z = z<sup>-1</sup></tt> of the given operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value of <tt>z<sup>-1</sup> = </tt>
	 *         <tt><span style="text-decoration: overline;">z</span>/|z|<sup>2</sup></tt> in {@code primitive-cartesian}
	 *         representation
	 */
	static public long inv (final long z) {
		final float x = re(z), y = im(z), norm = 1 / (x * x + y * y);
		return fromCartesian(x * norm, -y * norm);
	}


	/**
	 * Returns the scalar {@code sum} <tt>z' = z+v</tt>.
	 * @param z the complex operand to the left in {@code primitive-cartesian} representation
	 * @param v the scalar operand to the right
	 * @return the value <tt>z+v = &real;(z)+v+i&middot;&image;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long add (final long z, final float v) {
		return (z & ~MASK_HALFWORD) | (floatToRawIntBits(re(z) + v) & MASK_HALFWORD);
	}


	/**
	 * Returns the scalar {@code difference} <tt>z' = z-v</tt>.
	 * @param z the complex operand to the left in {@code primitive-cartesian} representation
	 * @param v the scalar operand to the right
	 * @return the value <tt>z-v = &real;(z)-v+i&middot;&image;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long sub (final long z, final float v) {
		return (z & ~MASK_HALFWORD) | (floatToRawIntBits(re(z) - v) & MASK_HALFWORD);
	}


	/**
	 * Returns the scalar {@code product} <tt>z' = z&middot;v</tt>.
	 * @param z the complex operand to the left in {@code primitive-cartesian} representation
	 * @param v the scalar operand to the right
	 * @return the value <tt>z&middot;v = &real;(z)&middot;v+i&middot;&image;(z)&middot;v</tt> in {@code primitive-cartesian}
	 *         representation
	 */
	static public long mul (final long z, final float v) {
		return fromCartesian(re(z) * v, im(z) * v);
	}


	/**
	 * Returns the scalar {@code quotient} <tt>z' = z/v</tt>.
	 * @param z the complex operand to the left in {@code primitive-cartesian} representation
	 * @param v the scalar operand to the right
	 * @return the value <tt>z/v = </tt> <tt>&real;(z)&middot;v<sup>-1</sup>+i&middot;&image;(z)&middot;v<sup>-1</sup></tt> in
	 *         {@code primitive-cartesian} representation
	 */
	static public long div (final long z, final float v) {
		return mul(z, 1 / v);
	}


	/**
	 * Returns the {@code sum} <tt>z' = z<sub>1</sub>+z<sub>2</sub></tt>.
	 * @param z1 the complex operand to the left in {@code primitive-cartesian} representation
	 * @param z2 the complex operand to the right in {@code primitive-cartesian} representation
	 * @return the value <tt>z<sub>1</sub>+z<sub>2</sub> = </tt>
	 *         <tt>&real;(z<sub>1</sub>)+&real;(z<sub>2</sub>)+i&middot;(&image;(z<sub>1</sub>)+&image;(z<sub>2</sub>))</tt> in
	 *         {@code primitive-cartesian} representation
	 */
	static public long add (final long z1, final long z2) {
		return fromCartesian(re(z1) + re(z2), im(z1) + im(z2));
	}


	/**
	 * Returns the {@code difference} <tt>z' = z<sub>1</sub>-z<sub>2</sub></tt>.
	 * @param z1 the complex operand to the left in {@code primitive-cartesian} representation
	 * @param z2 the complex operand to the right in {@code primitive-cartesian} representation
	 * @return the value <tt>z<sub>1</sub>-z<sub>2</sub> = </tt>
	 *         <tt>&real;(z<sub>1</sub>)-&real;(z<sub>2</sub>)+i&middot;(&image;(z<sub>1</sub>)-&image;(z<sub>2</sub>))</tt> in
	 *         {@code primitive-cartesian} representation
	 */
	static public long sub (final long z1, final long z2) {
		return fromCartesian(re(z1) - re(z2), im(z1) - im(z2));
	}


	/**
	 * Returns the {@code product} <tt>z' = z<sub>1</sub>&middot;z<sub>2</sub></tt>.
	 * @param z1 the complex operand to the left in {@code primitive-cartesian} representation
	 * @param z2 the complex operand to the right in {@code primitive-cartesian} representation
	 * @return the value <tt>z<sub>1</sub>&middot;z<sub>2</sub> = </tt>
	 *         <tt>&real;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)-&image;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>) + </tt>
	 *         <tt>i&middot;(&image;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)+&real;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>))</tt>
	 *         in {@code primitive-cartesian} representation
	 */
	static public long mul (final long z1, final long z2) {
		final float re1 = re(z1), im1 = im(z1), re2 = re(z2), im2 = im(z2);
		return fromCartesian(re1 * re2 - im1 * im2, im1 * re2 + re1 * im2);
	}


	/**
	 * Returns the {@code quotient} <tt>z' = z<sub>1</sub>/z<sub>2</sub> = </tt>
	 * <tt>z<sub>1</sub>&middot;&#x305z<sub>2</sub>&middot;|&#x305z<sub>2</sub>|<sup>-2</sup></tt>.
	 * @param z1 the complex operand to the left in {@code primitive-cartesian} representation to the left
	 * @param z2 the complex operand to the right in {@code primitive-cartesian} representation to the right
	 * @return the value <tt>z<sub>1</sub>/z<sub>2</sub> = </tt>
	 *         <tt>(&real;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)+&image;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>) + </tt>
	 *         <tt>i&middot;(&image;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)-&real;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>))) / </tt>
	 *         <tt>(&real;(z<sub>2</sub>)<sup>2</sup>+&image;(z<sub>2</sub>)<sup>2</sup>)</tt> in {@code primitive-cartesian}
	 *         representation
	 */
	static public long div (final long z1, final long z2) {
		final float re1 = re(z1), im1 = im(z1), re2 = re(z2), im2 = im(z2);
		final float norm = 1 / (re2 * re2 + im2 * im2);
		return fromCartesian((re1 * re2 + im1 * im2) * norm, (im1 * re2 - re1 * im2) * norm);
	}


	/**
	 * Returns <tt>z<sup>2</sup> = z&middot;z</tt>, i.e. the square of a {@code complex} operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value <tt>z&middot;z = &real;(z)<sup>2</sup>-&image;(z)<sup>2</sup> +
	 *         i&middot;2&middot;&real;(z)&middot;&image;(z)</tt> in {@code primitive-cartesian} representation
	 */
	static public long sq (final long z) {
		final float re = re(z), im = im(z);
		return fromCartesian(re * re - im * im, 2 * re * im);
	}


	/**
	 * Returns <tt>z<sup>½</sup> = |z|<sup>&half;</sup>&middot;e<sup>&half;&phi;(z)&middot;i</sup></tt>
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the principal value <tt>|z|<sup>&half;</sup>&middot;e<sup>&half;&phi;(z)&middot;i</sup> = </tt>
	 *         <tt>(&half;&middot;(|z|+&real;(z)))<sup>&half;</sup> +
	 *         i&middot;sign(&image;(z))&middot;(&half;&middot;(|z|-&real;(z)))<sup>&half;</sup></tt> in
	 *         {@code primitive-cartesian} representation
	 */
	static public long sqrt (final long z) {
		final float re = re(z), im = im(z);
		if (re >= 0 & im == 0) return fromCartesian((float) Math.sqrt(re), 0);

		final float abs = FloatMath.abs(re, im);
		final float sqrtRe = (float) Math.sqrt((abs + re) / 2), sqrtIm = (float) Math.sqrt((abs - re) / 2);
		return fromCartesian(sqrtRe, im < 0 ? -sqrtIm : sqrtIm);
	}


	/**
	 * Returns <tt>z<sup>3</sup> = z&middot;z&middot;z</tt>, i.e. the cube of a {@code complex} operand {@code z}.
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the value <tt>z&middot;z&middot;z = </tt> <tt>&real;(z)<sup>3</sup>-3&real;(z)&image;(z)<sup>2</sup> + </tt>
	 *         <tt>i&middot;(3&real;(z)<sup>2</sup>&image;(z)-&image;(z)<sup>3</sup>)</tt> in {@code primitive-cartesian}
	 *         representation
	 */
	static public long cb (final long z) {
		final float re = re(z), im = im(z);
		return fromCartesian(re * re * re - 3 * re * im * im, 3 * re * re * im - im * im * im);
	}


	/**
	 * Returns <tt>z<sup>&#x2153;</sup> = </tt> <tt>|z|<sup>&#x2153;</sup>&middot;e<sup>&#x2153;&phi;(z)&middot;i</sup></tt>
	 * @param z the complex operand in {@code primitive-cartesian} representation
	 * @return the principal value <tt>|z|<sup>&#x2153;</sup>&middot;e<sup>&#x2153;&phi;(z)&middot;i</sup></tt> in
	 *         {@code primitive-cartesian} representation
	 */
	static public long cbrt (final long z) {
		final float re = re(z);
		if (isReal(z)) return fromCartesian((float) Math.cbrt(re), 0);

		final float abs = abs(z), arg = arg(z);
		return fromPolar((float) Math.cbrt(abs), arg / 3);
	}


	/**
	 * {@inheritDoc}<br />
	 * <br />
	 * Note that this implementation always throws an exception.
	 */
	public void fft (final boolean inverse, final boolean separate, final PrimitiveCartesian[] vector) throws UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}


	/**
	 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of {@code N} complex numbers.
	 * @param inverse whether or not an {@code inverse} fourier transform shall be performed
	 * @param separate whether or not <i>channel separation</i> shall be performed
	 * @param vector an array of complex numbers
	 * @throws NullPointerException if the given vector is {@code null}
	 * @throws IllegalArgumentException if the given vector's length is not a power of two
	 * @see FunctionTables#getPerfectShuffleTable(int)
	 * @see FunctionTables#getTrigonometricTable(int)
	 */
	static public void fft (final boolean inverse, final boolean separate, final long[] vector) throws NullPointerException, IllegalArgumentException {
		if (inverse) {
			for (int index = 0; index < vector.length; ++index)
				vector[index] = conj(vector[index]);
		} else if (separate) {
			entangle(inverse, vector);
		}

		fft(vector);

		if (inverse) {
			final float norm = scalb(1f, -floorLog2(vector.length));
			for (int index = 0; index < vector.length; ++index)
				vector[index] = conj(mul(vector[index], norm));
			if (separate) entangle(inverse, vector);
		}
	}


	/**
	 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of {@code N} complex numbers. Note that an
	 * <i>inverse</i> transform can be performed in two ways:
	 * <ul>
	 * <li>by conjugating both the argument and result of this operation, and additionally norming the result by {@code 1/N}.
	 * </li>
	 * <li>by swapping the real and imaginary parts of both the argument and the result, and additionally norming the result by
	 * {@code 1/N}.</i>
	 * </ul>
	 * <i>Implementation notice</i>: This is an extremely optimized variant of the {@code Cooley-Tukey} algorithm, based on
	 * <i>perfect shuffle</i> and <i>sine/cosine</i> function tables for caching. Note that the JIT-compiler will inline the
	 * lookups within these tables (and any other static methods called). Logarithm-based iteration additionally eliminates
	 * expensive divisions, leaving only cheap addition, subtraction, multiplication and shift operations to be performed. The
	 * number of index variables is reduced to increase the chance that the JIT-compiler can keep them in registers (to avoid
	 * expensive cache or worse memory lookups), while derived indices are cheaply recalculated on demand taking advantage of
	 * modern super-scalar processor designs.
	 * @param vector an array of <tt>N = 2<sup>magnitude</sup></tt> complex numbers in {@code primitive-cartesian}
	 *        representation
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument's length is not a power of two
	 */
	static public void fft (final long[] vector) throws NullPointerException, IllegalArgumentException {
		final int magnitude = floorLog2(vector.length);
		if (vector.length != 1 << magnitude) throw new IllegalArgumentException();
		if (magnitude == 0) return;

		for (final SwapEntry entry : FunctionTables.getPerfectShuffleTable(magnitude)) {
			swap(vector, entry.getLeft(), vector, entry.getRight());
		}

		final FunctionTables.Trigonometric trigonometricTable = FunctionTables.getTrigonometricTable(magnitude);
		for (int depth = 0; depth < magnitude; ++depth) {
			for (int offset = 0; offset < 1 << depth; ++offset) {
				final int angle = offset << magnitude - depth - 1;
				final long unit = fromCartesian((float) trigonometricTable.cos(angle), (float) trigonometricTable.sin(angle));

				for (int left = offset, right = left + (1 << depth); right < 1 << magnitude; left += 2 << depth, right += 2 << depth) {
					final long product = mul(vector[right], unit);
					vector[right] = sub(vector[left], product);
					vector[left] = add(vector[left], product);
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
	 * @param inverse {@code true} for detangling, {@code false} for entangling
	 * @param spectrum the spectrum
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given spectrum's length is odd
	 */
	static private void entangle (final boolean inverse, final long[] spectrum) {
		if ((spectrum.length & 1) == 1) throw new IllegalArgumentException();

		final int halfLength = spectrum.length / 2;
		final long swap = spectrum[halfLength];
		spectrum[halfLength] = spectrum[1];
		spectrum[1] = swap;

		if (inverse) {
			for (int asc = 1, desc = spectrum.length - 1; asc < desc; ++asc, --desc) {
				final long left = spectrum[asc], right = conj(spectrum[desc]);
				spectrum[desc] = imul(sub(right, left));
				spectrum[asc] = add(right, left);
			}
		} else {
			for (int asc = 1, desc = spectrum.length - 1; asc < desc; ++asc, --desc) {
				final long left = spectrum[asc], right = imul(spectrum[desc]);
				spectrum[desc] = mul(conj(sub(left, right)), -.5f);
				spectrum[asc] = mul(add(left, right), +.5f);
			}
		}
	}
}