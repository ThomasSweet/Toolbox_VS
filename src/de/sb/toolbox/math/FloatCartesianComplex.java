package de.sb.toolbox.math;

import static de.sb.toolbox.math.IntMath.floorLog2;
import static java.lang.Math.cos;
import static java.lang.Math.scalb;
import static java.lang.Math.sin;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.math.FunctionTables.SwapEntry;


/**
 * This class implements mutable {@code complex numbers} that store single precision {@Cartesian} coordinates.
 */
@Copyright(year = 2008, holders = "Sascha Baumeister")
public final class FloatCartesianComplex extends Complex.AbstractSinglePrecision<FloatCartesianComplex>implements Complex.MutableSinglePrecision<FloatCartesianComplex> {
	static private final long serialVersionUID = 1;
	static private final ForkJoinPool FORK_JOIN_POOL = new ForkJoinPool(Runtime.getRuntime().availableProcessors());
	static private final int PARALLEL_MAGNITUDE_THRESHOLD = 14;

	private float re;
	private float im;


	/**
	 * Construct a new instance with a real and an imaginary part of zero.
	 */
	public FloatCartesianComplex () {
		this(0, 0);
	}


	/**
	 * Construct a new instance with the given real part and an imaginary part of zero.
	 * @param real the real part
	 */
	public FloatCartesianComplex (final float re) {
		this(re, 0);
	}


	/**
	 * Construct a new instance with the given real and imaginary parts.
	 * @param re the real part
	 * @param im the imaginary part
	 */
	public FloatCartesianComplex (final float re, final float im) {
		this.re = re;
		this.im = im;
	}


	/**
	 * Construct a new instance with the given complex value.
	 * @param value the complex value
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	public FloatCartesianComplex (final Complex.SinglePrecision<?> value) throws NullPointerException {
		this(value.re(), value.im());
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex setCartesian (final float re, final float im) {
		this.re = re;
		this.im = im;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex setPolar (final float abs, final float arg) {
		this.re = abs * (float) cos(arg);
		this.im = abs * (float) sin(arg);
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isReal () {
		return this.im == 0;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isImaginary () {
		return this.re == 0;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isZero () {
		return this.re == 0 & this.im == 0;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isInfinite () {
		return Float.isInfinite(this.re) | Float.isInfinite(this.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isNaN () {
		return Float.isNaN(this.re) | Float.isNaN(this.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public float re () {
		return this.re;
	}


	/**
	 * {@inheritDoc}
	 */
	public float im () {
		return this.im;
	}


	/**
	 * {@inheritDoc}
	 */
	public float abs () {
		return FloatMath.abs(this.re, this.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public float arg () {
		return FloatMath.arg(this.re, this.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex conj () {
		this.im = -this.im;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex neg () {
		return this.setCartesian(-this.re, -this.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex imul () {
		return this.setCartesian(-this.im, +this.re);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex idiv () {
		return this.setCartesian(+this.im, -this.re);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex inv () {
		final float re = this.re, im = this.im, norm = 1 / (re * re + im * im);
		return this.setCartesian(+re * norm, -im * norm);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex add (final float value) {
		this.re += value;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex sub (final float value) {
		this.re -= value;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex mul (final float value) {
		this.re *= value;
		this.im *= value;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex div (final float value) {
		return this.mul(1 / value);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex set (final FloatCartesianComplex value) throws NullPointerException {
		return this.setCartesian(value.re, value.im);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex add (final FloatCartesianComplex value) throws NullPointerException {
		this.re += value.re;
		this.im += value.im;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex sub (final FloatCartesianComplex value) throws NullPointerException {
		this.re -= value.re;
		this.im -= value.im;
		return this;
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex mul (final FloatCartesianComplex value) throws NullPointerException {
		final float re1 = this.re, im1 = this.im, re2 = value.re, im2 = value.im;
		return this.setCartesian(re1 * re2 - im1 * im2, re1 * im2 + im1 * re2);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex div (final FloatCartesianComplex value) throws NullPointerException {
		final float re1 = this.re, im1 = this.im, re2 = value.re, im2 = value.im, norm = 1 / (re2 * re2 + im2 * im2);
		return this.setCartesian((re1 * re2 + im1 * im2) * norm, (im1 * re2 - re1 * im2) * norm);
	}


	/**
	 * {@inheritDoc}
	 */
	public void mux (final FloatCartesianComplex value) throws NullPointerException {
		final float re1 = this.re, im1 = this.im, re2 = value.re, im2 = value.im;
		this.setCartesian(re1 + re2, im1 + im2);
		value.setCartesian(re1 - re2, im1 - im2);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex sq () {
		final float re = this.re, im = this.im;
		return this.setCartesian(re * re - im * im, 2 * re * im);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex sqrt () {
		final float re = this.re, im = this.im;
		if (re >= 0 & im == 0) {
			this.re = (float) Math.sqrt(re);
			return this;
		}

		final float abs = this.abs(), signum = im >= 0 ? +1 : -1;
		return this.setCartesian((float) Math.sqrt((abs + re) / 2), (float) Math.sqrt((abs - re) / 2) * signum);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex cb () {
		final float re = this.re, im = this.im;
		return this.setCartesian(re * re * re - 3 * re * im * im, 3 * re * re * im - im * im * im);
	}


	/**
	 * {@inheritDoc}
	 */
	public FloatCartesianComplex cbrt () {
		if (this.isReal()) {
			this.re = (float) Math.cbrt(this.re);
			return this;
		}

		final float abs = (float) Math.cbrt(this.abs());
		final float arg = this.arg() / 3;
		return this.setPolar(abs, arg);
	}


	/**
	 * {@inheritDoc}
	 * @see #fft(FloatCartesianComplex[], int, int)
	 */
	@Override
	public void fft (final boolean inverse, final boolean separate, final FloatCartesianComplex[] vector) throws NullPointerException, IllegalArgumentException {
		final int magnitude = floorLog2(vector.length);
		if (1 << magnitude != vector.length) throw new IllegalArgumentException();
		if (magnitude == 0) return;

		if (inverse) {
			for (int index = 0; index < vector.length; ++index)
				vector[index].conj();
		} else if (separate) {
			entangle(inverse, vector);
		}

		if (vector.length < PARALLEL_MAGNITUDE_THRESHOLD) {
			this.fft(magnitude, vector);
		} else {
			FORK_JOIN_POOL.invoke(new RecursiveFourierTransformer(magnitude, vector));
		}

		if (inverse) {
			final float norm = (float) scalb(1d, -magnitude);
			for (int index = 0; index < vector.length; ++index)
				vector[index].mul(norm).conj();
			if (separate) entangle(inverse, vector);
		}
	}


	/**
	 * Performs an <i>in-place Fast Fourier Transform</i> of the given vector of {@code N} complex numbers. Note that an
	 * <i>inverse</i> transform can be performed in two ways:
	 * <ul>
	 * <li>by conjugating both the argument and result of this operation, and additionally norming the result by {@code 1/N}.
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
	 * @param magnitude the value <tt>log<sub>2</sub>(N)</tt>
	 * @param vector an array of <tt>N = 2<sup>magnitude</sup></tt> complex numbers
	 * @throws NullPointerException if the given vector is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if the given vector's length is not a power of two
	 */
	private void fft (final int magnitude, final FloatCartesianComplex[] vector) throws NullPointerException, ArrayIndexOutOfBoundsException {
		assert 1 << magnitude == vector.length;

		final FloatCartesianComplex tao = new FloatCartesianComplex();
		for (final SwapEntry entry : FunctionTables.getPerfectShuffleTable(magnitude)) {
			tao.set(vector[entry.getLeft()]);
			vector[entry.getLeft()].set(vector[entry.getRight()]);
			vector[entry.getRight()].set(tao);
		}

		final FunctionTables.Trigonometric trigonometricTable = FunctionTables.getTrigonometricTable(magnitude);
		for (int depth = 0; depth < magnitude; ++depth) {
			for (int offset = 0; offset < 1 << depth; ++offset) {
				final int angleIndex = offset << magnitude - depth - 1;
				this.setCartesian((float) trigonometricTable.cos(angleIndex), (float) trigonometricTable.sin(angleIndex));

				for (int index = offset; index < 1 << magnitude; index += 2 << depth) {
					vector[index].mux(vector[index + (1 << depth)].mul(this));
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
	static private void entangle (final boolean inverse, final FloatCartesianComplex[] spectrum) throws NullPointerException, IllegalArgumentException {
		if ((spectrum.length & 1) == 1) throw new IllegalArgumentException();

		final int halfLength = spectrum.length / 2;
		final FloatCartesianComplex left = new FloatCartesianComplex();
		spectrum[halfLength].set(spectrum[1]);
		spectrum[1].set(left);

		final FloatCartesianComplex right = new FloatCartesianComplex();
		if (inverse) {
			for (int asc = 1, desc = spectrum.length - 1; asc < desc; ++asc, --desc) {
				left.set(spectrum[asc]);
				right.set(spectrum[desc].conj());
				spectrum[desc].sub(left).imul();
				spectrum[asc].add(right);
			}
		} else {
			for (int asc = 1, desc = spectrum.length - 1; asc < desc; ++asc, --desc) {
				left.set(spectrum[asc]);
				right.set(spectrum[desc].imul());
				spectrum[desc].sub(left).conj().mul(-.5f);
				spectrum[asc].add(right).mul(+.5f);
			}
		}
	}



	/**
	 * Resultless {@link ForkJoinTask} modelling a parallel-recursive FFT task, based on an optimized algorithm of (Danielson &
	 * Lanczos}, 1942.
	 */
	static private class RecursiveFourierTransformer extends RecursiveAction {
		static private final long serialVersionUID = 1L;

		private final int magnitude;
		private final FloatCartesianComplex[] vector;


		/**
		 * Creates a new instance.
		 * @param magnitude the value <tt>log<sub>2</sub>(N)</tt>
		 * @param vector an array of <tt>N = 2<sup>magnitude</sup></tt> complex numbers
		 */
		public RecursiveFourierTransformer (final int magnitude, final FloatCartesianComplex[] vector) {
			super();
			this.magnitude = magnitude;
			this.vector = vector;
		}


		/**
		 * {@inheritDoc}
		 */
		protected void compute () {
			final FloatCartesianComplex unit = new FloatCartesianComplex();
			if (this.magnitude < PARALLEL_MAGNITUDE_THRESHOLD) {
				unit.fft(this.magnitude, this.vector);
				return;
			}

			// prepare stage: divide vector into even and odd indexed parts
			final int half = 1 << this.magnitude - 1;
			final FloatCartesianComplex[] even = new FloatCartesianComplex[half], odd = new FloatCartesianComplex[half];
			for (int index = 0; index < half; ++index) {
				even[index] = this.vector[2 * index + 0];
				odd[index] = this.vector[2 * index + 1];
			}

			// divide stage: transform partial terms
			final RecursiveFourierTransformer evenTask = new RecursiveFourierTransformer(this.magnitude - 1, even);
			final RecursiveFourierTransformer oddTask = new RecursiveFourierTransformer(this.magnitude - 1, odd);
			oddTask.fork();
			evenTask.compute();
			oddTask.join();

			// conquer stage: recombine partial results
			final FunctionTables.Trigonometric trigonometricTable = FunctionTables.getTrigonometricTable(this.magnitude);
			for (int index = 0; index < half; ++index) {
				this.vector[index] = even[index];
				this.vector[index + half] = odd[index];
				unit.setCartesian((float) trigonometricTable.cos(index), (float) trigonometricTable.sin(index));
				even[index].mux(odd[index].mul(unit));
			}
		}
	}
}