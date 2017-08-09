package de.sb.toolbox.math;

import de.sb.toolbox.Copyright;


/**
 * Instances of this interface model <tt>complex numbers z&isin;&#x2102;</tt>, either represented in {@code Cartesian} form
 * <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>.
 * @param <T> the declaration type of the method receiver
 */
@Copyright(year = 2008, holders = "Sascha Baumeister")
public interface Complex<T extends Complex<T>> extends Cloneable {

	/**
	 * A format string for printing {@Cartesian} coordinates.
	 */
	static final String CARTESIAN_FORMAT = "%.10G%+.10G\u00B7i";

	/**
	 * A format string for printing {@Polar} coordinates.
	 */
	static final String POLAR_FORMAT = "%.10G\u00B7e^(%+.10G\u00B7\u03C0i)";


	/**
	 * Returns a clone of this complex number {@code z}.
	 * @return the clone {@code z' = z}
	 */
	T clone ();


	/**
	 * Returns {@code true} if this complex number represents a {@code real} number including {@code zero}, {@code false}
	 * otherwise. The imaginary part &image;(z) of real numbers is {@code 0}, and their argument angle
	 * <tt>((&phi;(z)+&pi;) mod<sub>e</sub> &tau;)-&pi;</tt> is {@code 0} or <tt>&plusmn;&pi;</tt> (except for zero).
	 * @return whether or not this complex number represents a real number.
	 */
	boolean isReal ();


	/**
	 * Returns {@code true} if this complex number represents an {@code imaginary} number including {@code zero}, {@code false}
	 * otherwise. The real part &real;(z) of imaginary numbers is {@code 0}, and their argument angle
	 * <tt>((&phi;(z)+&pi;) mod<sub>e</sub> &tau;)-&pi;</tt> is <tt>&plusmn;&half;&pi;</tt> (except for zero).
	 * @return whether or not this complex number represents an imaginary number.
	 */
	boolean isImaginary ();


	/**
	 * Returns {@code true} if this complex number represents {@code zero}, {@code false} otherwise.
	 * @return whether or not this complex number represents zero.
	 */
	boolean isZero ();


	/**
	 * Returns {@code true} if this complex number represents {@code infinity}, {@code false} otherwise.
	 * @return whether or not this complex number represents {@code infinity}.
	 */
	boolean isInfinite ();


	/**
	 * Returns {@code true} if this complex number represents {@code NaN}, {@code false} otherwise.
	 * @return whether or not this complex number represents {@code NaN}.
	 */
	boolean isNaN ();


	/**
	 * Returns the {@code Cartesian} text representation of this complex number {@code z}.
	 * @return the text representation of <tt>z = &real;(z)+&image;(z)&middot;i</tt>
	 */
	String toCartesianString ();


	/**
	 * Returns the {@code Polar} text representation of this complex number {@code z}.
	 * @return the text representation of <tt>z = |z|&middot;e<sup>i&middot;&phi;(z)</sup></tt>, with an argument angle
	 *         normalized to multiples of <tt>&pi;</tt>
	 */
	String toPolarString ();


	/**
	 * Returns a string representation of the given complex value, consisting of both {@code Cartesian} and {@code Polar}
	 * coordinates.
	 * @param z the complex value
	 * @return the string representation
	 * @see Object#toString()
	 */
	static public String toString (final Complex<?> z) {
		if (z.isNaN()) return "NaN";
		if (z.isInfinite()) return "INFINITY";
		return String.format("(%s == %s)", z.toCartesianString(), z.toPolarString());
	}



	/**
	 * Instances of this interface model {@code single-precision} complex numbers <tt>z &isin; &#x2102;</tt>, either represented
	 * in {@code Cartesian} form <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form
	 * <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>. Note that instances may be {@code immutable}.
	 * @param <T> the declaration type of the method receiver
	 */
	static interface SinglePrecision<T extends SinglePrecision<T>> extends Complex<T>, Comparable<SinglePrecision<?>> {

		/**
		 * Returns the {@code real} part <tt>&real;(z)</tt> of this complex number {@code z}.
		 * @return the real value <tt>&real;(z)</tt>
		 */
		float re ();


		/**
		 * Returns the {@code imaginary} part <tt>&image;(z)</tt> of this complex number {@code z}.
		 * @return the imaginary value <tt>&image;(z)</tt>
		 */
		float im ();


		/**
		 * Returns the absolute value <tt>|z|</tt> of this complex number {@code z}.
		 * @return the absolute value <tt>|z| = (&real;(z)<sup>2</sup>+&image;(z)<sup>2</sup>)<sup>&half;</sup></tt>
		 */
		float abs ();


		/**
		 * Returns the argument angle <tt>&phi;(z)</tt> of this complex number {@code z}.
		 * @return the argument angle <tt>&phi;(z) = atan2(&image;(z),&real;(z))</tt>
		 */
		float arg ();


		/**
		 * Compares this object with the specified operand for order, with {@code null} values considered smaller than any other
		 * value except {@code null}, and {@code NaN} values considered greater than any other value except {@code NaN}. The
		 * operation returns a strictly negative integer if <tt>&real;(this) < &real;(z)</tt>, or a stictly positive integer if
		 * <tt>&real;(this) > &real;(z)</tt>. If both real parts are equal, then the operation returns a strictly negative
		 * integer if <tt>&image;(this) < &image;(z)</tt>, or a strictly positive integer if <tt>&image;(this) > &image;(z)</tt>
		 * . The operation returns {@code zero} if and only if both the receiver and the given operand share the same real and
		 * imaginary parts.<br />
		 * <br />
		 * Note that the {@code natural order} thus imposed on complex values is not a {@code total order} in strict
		 * mathematical terms, as complex numbers are proven to produce contradictions under any kind of order. Also note that
		 * the definition above implies that {@code INFINITE} values are distinguishable by their real parts, imaginary parts
		 * and their respective signs, while they are indistinguishable for {@code equality} operations based on a
		 * {@code Riemann Sphere} model.
		 * @param z the {@code single precision} operand to be compared.
		 * @return a strictly negative integer, zero, or a strictly positive integer as this value is less than, equal to, or
		 *         greater than the specified value
		 * @see #compare(SinglePrecision, SinglePrecision)
		 */
		default int compareTo (final SinglePrecision<?> z) {
			return compare(this, z);
		}


		/**
		 * Returns a hash code for a {@code single precision} complex operand. {@code null} values always return {@code zero},
		 * {@code INFINITY} values always return {@code Float.hashCode(POSITIVE_INFINITY)}, and {@code NaN} values always return
		 * {@code Float.hashCode(NaN)}.
		 * @param z the complex operand
		 * @return the hash code
		 * @see Object#hashCode()
		 * @see #equals(SinglePrecision, SinglePrecision)
		 */
		static int hashCode (final SinglePrecision<?> z) {
			if (z == null) return 0;
			if (z.isNaN()) return Float.hashCode(Float.NaN);
			if (z.isInfinite()) return Float.hashCode(Float.POSITIVE_INFINITY);
			return Float.hashCode(z.re()) ^ Float.hashCode(z.im());
		}


		/**
		 * Indicates whether or not two {@code double precision} complex values are equal. Complex values are compared using an
		 * extended {@code Riemann sphere} model: two values are considered equal if they are both {@code null}, both
		 * {@code INFINITY}, or share the same real and imaginary parts. {@code NaN} values are never consider equal.
		 * @param left the left complex operand
		 * @param right the right complex operand
		 * @return whether or not the given complex values are equal
		 * @see Object#equals(Object)
		 * @see #hashCode()
		 */
		static boolean equals (final SinglePrecision<?> left, final SinglePrecision<?> right) {
			if (left == null | right == null) return left == null & right == null;
			if (left.isNaN() | right.isNaN()) return false;
			if (left.isInfinite() | right.isInfinite()) return left.isInfinite() & right.isInfinite();
			return left.re() == right.re() & left.im() == right.im();
		}


		/**
		 * Compares this object with the specified operand for order, with {@code null} values considered smaller than any other
		 * value except {@code null}, and {@code NaN} values considered greater than any other value except {@code NaN}. The
		 * operation returns a strictly negative integer if <tt>&real;(z<sub>1</sub>) < &real;(z<sub>2</sub>)</tt>, or a stictly
		 * positive integer if <tt>&real;(z<sub>1</sub>) > &real;(z<sub>2</sub>)</tt>. If both real parts are equal, then the
		 * operation returns a strictly negative integer if <tt>&image;(z<sub>1</sub>) < &image;(z<sub>2</sub>)</tt>, or a
		 * strictly positive integer if <tt>&image;(z<sub>1</sub>) > &image;(z<sub>2</sub>)</tt>. The operation returns
		 * {@code zero} if and only if both the receiver and the given operand share the same real and imaginary parts.<br />
		 * <br />
		 * Note that the {@code natural order} thus imposed on complex values is not a {@code total order} in strict
		 * mathematical terms, as complex numbers are proven to produce contradictions under any kind of order. Also note that
		 * the definition above implies that {@code INFINITE} values are distinguishable by their real parts, imaginary parts
		 * and their respective signs, while they are indistinguishable for {@code equality} checks based on a
		 * {@code Riemann Sphere} model.
		 * @param z1 the left {@code single precision} operand to be compared.
		 * @param z2 the right {@code single precision} operand to be compared.
		 * @return a strictly negative integer, zero, or a strictly positive integer as <tt>z<sub>1</sub></tt> is less than,
		 *         equal to, or greater than <tt>z<sub>2</sub></tt>
		 */
		static int compare (final SinglePrecision<?> z1, final SinglePrecision<?> z2) {
			if (z1 == null | z2 == null) return z1 == z2 ? 0 : (z1 == null ? -1 : +1);
			if (z1.isNaN() | z2.isNaN()) return z1.isNaN() == z2.isNaN() ? 0 : (z1.isNaN() ? +1 : -1);
			final int compare = Float.compare(z1.re(), z2.re());
			return compare == 0 ? Float.compare(z1.im(), z2.im()) : compare;
		}
	}



	/**
	 * Instances of this interface model {@code double-precision} complex numbers <tt>z &isin; &#x2102;</tt>, either represented
	 * in {@code Cartesian} form <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form
	 * <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>. Note that instances may be {@code immutable}.
	 * @param <T> the declaration type of the method receiver
	 */
	static interface DoublePrecision<T extends DoublePrecision<T>> extends Complex<T>, Comparable<DoublePrecision<?>> {

		/**
		 * Returns the {@code real} part <tt>&real;(z)</tt> of this complex number {@code z}.
		 * @return the real value <tt>&real;(z)</tt>
		 */
		double re ();


		/**
		 * Returns the {@code imaginary} part <tt>&image;(z)</tt> of this complex number {@code z}.
		 * @return the imaginary value <tt>&image;(z)</tt>
		 */
		double im ();


		/**
		 * Returns the {@code absolute} value <tt>|z|</tt> of this complex number {@code z}.
		 * @return the absolute value <tt>|z| = (&real;(z)<sup>2</sup>+&image;(z)<sup>2</sup>)<sup>&half;</sup></tt>
		 */
		double abs ();


		/**
		 * Returns the {@code argument} angle <tt>&phi;(z)</tt> of this complex number {@code z}.
		 * @return the argument angle <tt>&phi;(z) = atan2(&image;(z),&real;(z))</tt>
		 */
		double arg ();


		/**
		 * Compares this object with the specified operand for order, with {@code null} values considered smaller than any other
		 * value except {@code null}, and {@code NaN} values considered greater than any other value except {@code NaN}. The
		 * operation returns a strictly negative integer if <tt>&real;(this) < &real;(z)</tt>, or a stictly positive integer if
		 * <tt>&real;(this) > &real;(z)</tt>. If both real parts are equal, then the operation returns a strictly negative
		 * integer if <tt>&image;(this) < &image;(z)</tt>, or a strictly positive integer if <tt>&image;(this) > &image;(z)</tt>
		 * . The operation returns {@code zero} if and only if both the receiver and the given operand share the same real and
		 * imaginary parts.<br />
		 * <br />
		 * Note that the {@code natural order} thus imposed on complex values is not a {@code total order} in strict
		 * mathematical terms, as complex numbers are proven to produce contradictions under any kind of order. Also note that
		 * the definition above implies that {@code INFINITE} values are distinguishable by their real parts, imaginary parts
		 * and their respective signs, while they are indistinguishable for {@code equality} operations based on a
		 * {@code Riemann Sphere} model.
		 * @param z the {@code double precision} operand to be compared.
		 * @return a strictly negative integer, zero, or a strictly positive integer as this value is less than, equal to, or
		 *         greater than the specified value
		 * @see #compare(DoublePrecision, DoublePrecision)
		 */
		default int compareTo (final DoublePrecision<?> z) {
			return compare(this, z);
		}


		/**
		 * Returns a hash code for a {@code double precision} complex operand. {@code null} values always return {@code zero},
		 * {@code INFINITY} values always return {@code Double.hashCode(POSITIVE_INFINITY)}, and {@code NaN} values always
		 * return {@code Double.hashCode(NaN)}.
		 * @param z the complex operand
		 * @return the hash code
		 * @see Object#hashCode()
		 * @see #equals(DoublePrecision, DoublePrecision)
		 */
		static int hashCode (final DoublePrecision<?> z) {
			if (z == null) return 0;
			if (z.isNaN()) return Double.hashCode(Double.NaN);
			if (z.isInfinite()) return Double.hashCode(Double.POSITIVE_INFINITY);
			return Double.hashCode(z.re()) ^ Double.hashCode(z.im());
		}


		/**
		 * Indicates whether or not two {@code double precision} complex values are equal. Complex values are compared using an
		 * extended {@code Riemann sphere} model: two values are considered equal if they are both {@code null}, both
		 * {@code INFINITY}, or share the same real and imaginary parts. {@code NaN} values are never consider equal.
		 * @param left the left complex operand
		 * @param right the right complex operand
		 * @return whether or not the given complex values are equal
		 * @see Object#equals(Object)
		 * @see #hashCode()
		 */
		static boolean equals (final DoublePrecision<?> left, final DoublePrecision<?> right) {
			if (left == null | right == null) return left == null & right == null;
			if (left.isNaN() | right.isNaN()) return false;
			if (left.isInfinite() | right.isInfinite()) return left.isInfinite() & right.isInfinite();
			return left.re() == right.re() & left.im() == right.im();
		}


		/**
		 * Compares this object with the specified operand for order, with {@code null} values considered smaller than any other
		 * value except {@code null}, and {@code NaN} values considered greater than any other value except {@code NaN}. The
		 * operation returns a strictly negative integer if <tt>&real;(z<sub>1</sub>) < &real;(z<sub>2</sub>)</tt>, or a stictly
		 * positive integer if <tt>&real;(z<sub>1</sub>) > &real;(z<sub>2</sub>)</tt>. If both real parts are equal, then the
		 * operation returns a strictly negative integer if <tt>&image;(z<sub>1</sub>) < &image;(z<sub>2</sub>)</tt>, or a
		 * strictly positive integer if <tt>&image;(z<sub>1</sub>) > &image;(z<sub>2</sub>)</tt>. The operation returns
		 * {@code zero} if and only if both the receiver and the given operand share the same real and imaginary parts.<br />
		 * <br />
		 * Note that the {@code natural order} thus imposed on complex values is not a {@code total order} in strict
		 * mathematical terms, as complex numbers are proven to produce contradictions under any kind of order. Also note that
		 * the definition above implies that {@code INFINITE} values are distinguishable by their real parts, imaginary parts
		 * and their respective signs, while they are indistinguishable for {@code equality} checks based on a
		 * {@code Riemann Sphere} model.
		 * @param z1 the left {@code double precision} operand to be compared.
		 * @param z2 the right {@code double precision} operand to be compared.
		 * @return a strictly negative integer, zero, or a strictly positive integer as <tt>z<sub>1</sub></tt> is less than,
		 *         equal to, or greater than <tt>z<sub>2</sub></tt>
		 */
		static int compare (final DoublePrecision<?> z1, final DoublePrecision<?> z2) {
			if (z1 == null | z2 == null) return z1 == z2 ? 0 : (z1 == null ? -1 : +1);
			if (z1.isNaN() | z2.isNaN()) return z1.isNaN() == z2.isNaN() ? 0 : (z1.isNaN() ? +1 : -1);
			final int compare = Double.compare(z1.re(), z2.re());
			return compare == 0 ? Double.compare(z1.im(), z2.im()) : compare;
		}
	}



	/**
	 * Instances of this interface model {@code mutable} complex numbers <tt>z &isin; &#x2102;</tt>, either represented in
	 * {@code Cartesian} form <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form
	 * <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>.
	 * @param <T> the declaration type of the method receiver
	 */
	static interface Mutable<T extends Mutable<T>> extends Complex<T> {

		/**
		 * Assigns the {@code conjugated} value <tt>z'= <span style="text-decoration: overline;">z</span>
		 * = +&real;(z)-i&middot;&image;(z)</tt> to this complex number {@code z}. This effectively reflects the number on the
		 * x-axis.
		 * @return the modified receiver
		 */
		T conj ();


		/**
		 * Assigns the {@code negative} value <tt>z'= -z = -&real;(z)-i&middot;&image;(z)</tt> to this complex number {@code z}.
		 * This effectively reflects the number on the origin.
		 * @return the modified receiver
		 */
		T neg ();


		/**
		 * Assigns the {@code imaginary product} value <tt>z'= +i&middot;z = -&image;(z)+i&middot;&real;(z)</tt> to this complex
		 * number {@code z}. This effectively turns the number by 90&deg; counter-clockwise around the origin.
		 * @return the modified receiver
		 */
		T imul ();


		/**
		 * Assigns the {@code imaginary quotient} value <tt>z'= z/i = -i&middot;z = +&image;(z)-i&middot;&real;(z)</tt> to this
		 * complex number {@code z}. This effectively turns the number by 90&deg; clockwise around the origin.
		 * @return the modified receiver
		 */
		T idiv ();


		/**
		 * Assigns the {@code reciprocal} value <tt>z'= z<sup>-1</sup> = 1/z =
		 * <span style="text-decoration: overline;">z</span>/|z|<sup>2</sup></tt> to this complex number {@code z}.
		 * @return the modified receiver
		 */
		T inv ();


		/**
		 * Performs the state assignment <tt>z -> z<sup>2</sup> = z&middot;z</tt>. Note that this operation is dual-cyclic, with
		 * (-z)<sup>2</sup> = z<sup>2</sup>. This cyclic behavior causes the square root operation (as the inverse of this
		 * operation) to become dual-valued in <tt>&#x2102;</tt>, similarly to the sine function's cyclic behavior causing the
		 * arc sine function to become multi-valued in <tt>&#x211D;</tt>.
		 * @return the modified receiver
		 */
		T sq ();


		/**
		 * Performs the state assignment <tt>z -> z<sup>&half;</sup> = </tt>
		 * <tt>|z|<sup>&half;</sup>&middot;e<sup>&half;&phi;(z)&middot;i</sup></tt>. The dual-valued nature of this operation
		 * implies a discontinuity, which is located at negative real values of the operand {@code z}. For receivers with
		 * argument angles within range <tt>[-&pi;,+&pi;[</tt>, this operation always results in the the principal branch
		 * solution, i.e. the one whose real part is positive. The solution for the other branch can be obtained by multiplying
		 * the principal branch solution with <tt>e<sup>i&middot;&pi;</sup></tt>, i.e. negating it.
		 * @return the modified receiver
		 */
		T sqrt ();


		/**
		 * Performs the state assignment <tt>z -> z<sup>3</sup> = z&middot;z&middot;z</tt>. Note that this operation is
		 * trial-cyclic, with (-z)<sup>2</sup> = z<sup>2</sup>. This cyclic behavior causes the cube root operation (as the
		 * inverse of this operation) to become trial-valued in <tt>&#x2102;</tt>, similarly to the sine function's cyclic
		 * behavior causing the arc sine function to become multi-valued in <tt>&#x211D;</tt>.
		 * @return the modified receiver
		 */
		T cb ();


		/**
		 * Performs the state assignment <tt>z -> z<sup>&#x2153;</sup> = </tt>
		 * <tt>|z|<sup>&#x2153;</sup>&middot;e<sup>&#x2153;&phi;(z)&middot;i</sup></tt>. The trial-valued nature of this
		 * operation implies a discontinuity, which is located at negative real values of the operand {@code z}. For receivers
		 * that represent real numbers, this operation always results in a real number. For any others within argument angle
		 * range <tt>[-&pi;,+&pi;[</tt>, it results in the principal branch solution, i.e. the one whose argument (positive or
		 * negative) is closest to zero; it's real part is guaranteed to be positive. The solutions for the other branches can
		 * be obtained by multiplying the principal branch solution with <i>e <sup>&#x2153;k&middot;&tau;&middot;i</sup></i>.
		 * @return the modified receiver
		 */
		T cbrt ();


		/**
		 * Performs the state assignment <tt>z<sub>1</sub> = z<sub>2</sub></tt>.
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @return the modified receiver
		 * @throws NullPointerException if the given argument is {@code null}
		 */
		T set (T value) throws NullPointerException;


		/**
		 * Performs the state assignment <tt>z<sub>1</sub> += z<sub>2</sub> = </tt>
		 * <tt>&real;(z<sub>1</sub>)+&real;(z<sub>2</sub>)+i&middot;(&image;(z<sub>1</sub>)+&image;(z<sub>2</sub>))</tt>.
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @return the modified receiver
		 * @throws NullPointerException if the given argument is {@code null}
		 */
		T add (T value) throws NullPointerException;


		/**
		 * Performs the state assignment <tt>z<sub>1</sub> -= z<sub>2</sub> = </tt>
		 * <tt>&real;(z<sub>1</sub>)-&real;(z<sub>2</sub>)+i&middot;(&image;(z<sub>1</sub>)-&image;(z<sub>2</sub>))</tt>.
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @return the modified receiver
		 * @throws NullPointerException if the given argument is {@code null}
		 */
		T sub (T value) throws NullPointerException;


		/**
		 * Performs the state assignment <tt>z<sub>1</sub> *= z<sub>2</sub> = </tt>
		 * <tt>&real;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)-&image;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>) + </tt>
		 * <tt>i&middot;(&image;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)+&real;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>))</tt>
		 * .
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @return the modified receiver
		 * @throws NullPointerException if the given argument is {@code null}
		 */
		T mul (T value) throws NullPointerException;


		/**
		 * Performs the state assignment <tt>z<sub>1</sub> /= z<sub>2</sub> = </tt>
		 * <tt>(&real;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)+&image;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>) + </tt>
		 * <tt>i&middot;(&image;(z<sub>1</sub>)&middot;&real;(z<sub>2</sub>)-&real;(z<sub>1</sub>)&middot;&image;(z<sub>2</sub>))) / </tt>
		 * <tt>(&real;(z<sub>2</sub>)<sup>2</sup>+&image;(z<sub>2</sub>)<sup>2</sup>)</tt>.
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @return the modified receiver
		 * @throws NullPointerException if the given argument is {@code null}
		 */
		T div (T value) throws NullPointerException;


		/**
		 * Assigns the values <tt>z<sub>1</sub>' = z<sub>1</sub>+z<sub>2</sub></tt> and
		 * <tt>z<sub>2</sub>' = z<sub>1</sub>-z<sub>2</sub></tt> to the receiver <tt>z<sub>1</sub></tt> and the given value
		 * <tt>z<sub>2</sub></tt> respectively.
		 * @param <T> the declaration type of the method receiver
		 * @param value the complex operand <tt>z<sub>2</sub></tt>
		 * @throws NullPointerException if any of the given arguments is {@code null}
		 */
		void mux (T value) throws NullPointerException;


		/**
		 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of {@code N} complex numbers.<br />
		 * <br />
		 * <i>Implementation notice</i>: Note that an <i>inverse</i> transform can be performed in two ways:
		 * <ul>
		 * <li>by conjugating both the argument and result of this operation, and additionally norming the result by {@code 1/N}
		 * .</li>
		 * <li>by swapping the real and imaginary parts of both the argument and the result, and additionally norming the result
		 * by {@code 1/N}.</i>
		 * </ul>
		 * The method receiver can be used to store temporary data like complex units, but implementations are not required to
		 * do so. Apart from that it's sole purpose is to resolve polymorphic implementations of this algorithm. Note that there
		 * is no default implementation because runtime tests have shown a drop in performance of around 40% if virtual method
		 * calls need to be resolved during the operation; more restrictive type declatations for the complex numbers involved
		 * (like the method signature allowing overriding with the implementing type as array component type) allow the JIT
		 * compiler to inline code more aggressively. Reuse of the referenced function tables is highly recommended as well.
		 * @param <T> the declaration type of the method receiver
		 * @param inverse whether or not an {@code inverse} fourier transform shall be performed
		 * @param separate whether or not {@code channel separation} shall be performed
		 * @param vector an array of <tt>N = 2<sup>magnitude</sup></tt> complex numbers
		 * @throws NullPointerException if the given vector is {@code null}
		 * @throws IllegalArgumentException if the given vector has no power of two length
		 * @throws UnsupportedOperationException if the type does not support fourier transform or channel separation
		 * @see FunctionTables#getPerfectShuffleTable(int)
		 * @see FunctionTables#getTrigonometricTable(int)
		 */
		void fft (final boolean inverse, final boolean separate, final T[] vector) throws NullPointerException, IllegalArgumentException, UnsupportedOperationException;
	}



	/**
	 * Instances of this interface model {@code mutable single precision} complex numbers <tt>z &isin; &#x2102;</tt>, either
	 * represented in {@code Cartesian} form <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form
	 * <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>.
	 * @param <T> the declaration type of the method receiver
	 */
	static interface MutableSinglePrecision<T extends MutableSinglePrecision<T>> extends Mutable<T>, SinglePrecision<T> {

		/**
		 * Performs the state assignment <tt>z = &real;(z)+i&middot;&image;(z)</tt>.
		 * @param re the real part <tt>&real;(z)</tt>
		 * @param im the imaginary part <tt>&image;(z)</tt>
		 * @return the modified receiver
		 */
		T setCartesian (final float re, final float im);


		/**
		 * Performs the state assignment <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup> = </tt>
		 * <tt>|z|&middot;cos(&phi;(z)) + i&middot;|z|&middot;sin(&phi;(z))</tt>.
		 * @param abs the absolute value <tt>|z|</tt>
		 * @param arg the argument angle <tt>&phi;(z)</tt>
		 * @return the modified receiver
		 */
		T setPolar (final float abs, final float arg);


		/**
		 * Performs the state assignment <tt>z += v</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T add (float value);


		/**
		 * Performs the state assignment <tt>z -= v</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T sub (float value);


		/**
		 * Performs the state assignment <tt>z *= v</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T mul (float value);


		/**
		 * Performs the state assignment <tt>z /= v</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T div (float value);
	}



	/**
	 * Instances of this interface model {@code mutable double precision} complex numbers <tt>z &isin; &#x2102;</tt>, either
	 * represented in {@code Cartesian} form <tt>z = &real;(z)+&image;(z)&middot;i</tt>, or in {@code Polar} form
	 * <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup></tt>.
	 * @param <T> the declaration type of the method receiver
	 */
	static interface MutableDoublePrecision<T extends MutableDoublePrecision<T>> extends Mutable<T>, DoublePrecision<T> {

		/**
		 * Performs the state assignment <tt>z = &real;(z)+i&middot;&image;(z)</tt>.
		 * @param re the real part <tt>&real;(z)</tt>
		 * @param im the imaginary part <tt>&image;(z)</tt>
		 * @return the modified receiver
		 */
		T setCartesian (final double re, final double im);


		/**
		 * Performs the state assignment <tt>z = |z|&middot;e<sup>&phi;(z)&middot;i</sup> = </tt>
		 * <tt>|z|&middot;cos(&phi;(z)) + i&middot;|z|&middot;sin(&phi;(z))</tt>.
		 * @param abs the absolute value <tt>|z|</tt>
		 * @param arg the argument angle <tt>&phi;(z)</tt>
		 * @return the modified receiver
		 */
		T setPolar (final double abs, final double arg);


		/**
		 * Performs the state assignment <tt>z += v = &real;(z)+v+i&middot;&image;(z)</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T add (double value);


		/**
		 * Performs the state assignment <tt>z -= v = &real;(z)-v+i&middot;&image;(z)</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T sub (double value);


		/**
		 * Performs the state assignment <tt>z *= v = &real;(z)&middot;v+i&middot;&image;(z)&middot;v</tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T mul (double value);


		/**
		 * Performs the state assignment <tt>z /= v = </tt>
		 * <tt>&real;(z)&middot;v<sup>-1</sup>+i&middot;&image;(z)&middot;v<sup>-1</sup></tt>.
		 * @param value the scalar value {@code v}
		 * @return the modified receiver
		 */
		T div (double value);
	}



	/**
	 * This abstract class models {@code complex numbers}.
	 * @param <T> the declaration type of the method receiver
	 */
	static public abstract class AbstractComplex<T extends AbstractComplex<T>> extends Number implements Complex<T> {
		static private final long serialVersionUID = 1L;


		/**
		 * {@inheritDoc}
		 */
		@Override
		public T clone () {
			try {
				@SuppressWarnings("unchecked")
				final T clone = (T) super.clone();
				return clone;
			} catch (final CloneNotSupportedException exception) {
				throw new AssertionError(exception);
			}
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public String toString () {
			return Complex.toString(this);
		}


		/**
		 * Returns the value of the specified number as a {@code double}, which may involve rounding.
		 * @return the numeric value represented by this object after conversion to type {@code double}.
		 */
		public abstract double doubleValue ();


		/**
		 * Returns the value of the specified number as a {@code float}, which may involve rounding.
		 * @return the numeric value represented by this object after conversion to type {@code float}.
		 */
		public abstract float floatValue ();


		/**
		 * Returns the real value <tt>&real;(z)</tt> of the specified number {@code z} as a {@code long}, which may involve
		 * rounding or truncation.
		 * @return the numeric value represented by this object after conversion to type {@code long}.
		 */
		public abstract long longValue ();


		/**
		 * Returns the real value <tt>&real;(z)</tt> of the specified number {@code z} as a {@code int}, which may involve
		 * rounding or truncation.
		 * @return the numeric value represented by this object after conversion to type {@code int}.
		 */
		public int intValue () {
			return (int) this.longValue();
		}


		/**
		 * Returns the real value <tt>&real;(z)</tt> of the specified number {@code z} as a {@code short}, which may involve
		 * rounding or truncation.
		 * @return the numeric value represented by this object after conversion to type {@code short}.
		 */
		public short shortValue () {
			return (short) this.longValue();
		}


		/**
		 * Returns the real value <tt>&real;(z)</tt> of the specified number {@code z} as a {@code byte}, which may involve
		 * rounding or truncation.
		 * @return the numeric value represented by this object after conversion to type {@code byte}.
		 */
		public byte byteValue () {
			return (byte) this.longValue();
		}
	}



	/**
	 * This abstract class models {@code complex numbers} that store single precision coordinates, without making any
	 * assumptions on whether subclasses store said coordinates in {@code Cartesian} or {@code Polar} form.
	 * @param <T> the declaration type of the method receiver
	 */
	static public abstract class AbstractSinglePrecision<T extends AbstractSinglePrecision<T>> extends AbstractComplex<T>implements SinglePrecision<T> {
		static private final long serialVersionUID = 1L;


		/**
		 * Returns a hash code for this complex operand. Note that {@code INFINITY} values always return
		 * {@code Float.hashCode(POSITIVE_INFINITY)}, and {@code NaN} values always return {@code Float.hashCode(NaN)}.
		 * @return the hash code
		 * @see Object#hashCode()
		 * @see #equals(Object)
		 */
		@Override
		public int hashCode () {
			return SinglePrecision.hashCode(this);
		}


		/**
		 * Indicates whether or not two {@code single precision} complex values are equal. Complex values are compared using an
		 * extended {@code Riemann sphere} model: two values are considered equal if they are both complex single precision
		 * values, and either both {@code INFINITY}, or sharing the same real and imaginary parts. {@code null} and {@code NaN}
		 * values are never considered equal.
		 * @param object the object to compare with
		 * @return whether or not the receiver and the given object are equal
		 * @see Object#equals(Object)
		 * @see #hashCode()
		 */
		@Override
		public boolean equals (final Object object) {
			if (!(object instanceof SinglePrecision)) return false;
			return SinglePrecision.equals(this, (SinglePrecision<?>) object);
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public double doubleValue () {
			return this.floatValue();
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public float floatValue () {
			return this.re();
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public long longValue () {
			return (long) this.floatValue();
		}


		/**
		 * {@inheritDoc}
		 */
		public String toCartesianString () {
			if (this.isNaN()) return "NaN";
			if (this.isInfinite()) return "INFINITY";
			return String.format(CARTESIAN_FORMAT, this.re(), this.im());
		};


		/**
		 * {@inheritDoc}
		 */
		public String toPolarString () {
			if (this.isNaN()) return "NaN";
			if (this.isInfinite()) return "INFINITY";
			return String.format(POLAR_FORMAT, this.abs(), this.arg() / Math.PI);
		}
	}



	/**
	 * This abstract class models {@code complex numbers} that store double precision coordinates, without making any
	 * assumptions on whether subclasses store said coordinates in {@code Cartesian} or {@code Polar} form.
	 * @param <T> the declaration type of the method receiver
	 */
	static public abstract class AbstractDoublePrecision<R extends AbstractDoublePrecision<R>> extends AbstractComplex<R>implements DoublePrecision<R> {
		static private final long serialVersionUID = 1L;


		/**
		 * Returns a hash code for this complex operand. Note that {@code INFINITY} values always return
		 * {@code Double.hashCode(POSITIVE_INFINITY)}, and {@code NaN} values always return {@code Double.hashCode(NaN)}.
		 * @return the hash code
		 * @see Object#hashCode()
		 * @see #equals(Object)
		 */
		@Override
		public int hashCode () {
			return DoublePrecision.hashCode(this);
		}


		/**
		 * Indicates whether or not two {@code double precision} complex values are equal. Complex values are compared using an
		 * extended {@code Riemann sphere} model: two values are considered equal if they are both complex single precision
		 * values, and either both {@code INFINITY}, or sharing the same real and imaginary parts. {@code null} and {@code NaN}
		 * values are never considered equal.
		 * @param object the object to compare with
		 * @return whether or not the receiver and the given object are equal
		 * @see Object#equals(Object)
		 * @see #hashCode()
		 */
		@Override
		public boolean equals (final Object object) {
			if (!(object instanceof DoublePrecision)) return false;
			return DoublePrecision.equals(this, (DoublePrecision<?>) object);
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public double doubleValue () {
			return this.re();
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public float floatValue () {
			return (float) this.doubleValue();
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public long longValue () {
			return (long) this.doubleValue();
		}


		/**
		 * {@inheritDoc}
		 */
		public String toCartesianString () {
			if (this.isNaN()) return "NaN";
			if (this.isInfinite()) return "INFINITY";
			return String.format(CARTESIAN_FORMAT, this.re(), this.im());
		};


		/**
		 * {@inheritDoc}
		 */
		public String toPolarString () {
			if (this.isNaN()) return "NaN";
			if (this.isInfinite()) return "INFINITY";
			return String.format(POLAR_FORMAT, this.abs(), this.arg() / Math.PI);
		}
	}
}