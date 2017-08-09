package de.sb.toolbox.util;

import static java.lang.Integer.highestOneBit;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import de.sb.toolbox.Copyright;


/**
 * Instances model buffers capable of collecting primitive values and returning them as a primitive array, without the need for
 * an intermediate boxed collection. Note that implementations are not <i>thread safe</i>, therefore all operations require
 * external synchronization if an instance is shared between threads.
 * @param <A> the array type
 */
@Copyright(year = 2017, holders = "Sascha Baumeister")
public interface ArrayBuffer<A> {
	static final int INITIAL_CAPACITY = 0xffff;


	/**
	 * Returns the content of this buffer as an array.
	 * @return the array content
	 */
	A get ();


	/**
	 * Appends the given array's values to this buffer.
	 * @param values the values array
	 * @throws NullPointerException if the given array is {@code null}
	 * @throws NegativeArraySizeException if appending the given values requires a buffer larger than {@link Integer#MAX_VALUE
	 *         MAX_VALUE}
	 */
	void put (A values) throws NullPointerException, NegativeArraySizeException;


	/**
	 * Appends a portion of the given array's values to this buffer.
	 * @param values the values array
	 * @param offset the array offset at which to start retrieving samples
	 * @param amount the number of values to append
	 * @throws NullPointerException if the given array is {@code null}
	 * @throws IllegalArgumentException if the given offset or amount is negative, or if their sum exceeds the given array's
	 *         length
	 * @throws NegativeArraySizeException if appending the given values requires a buffer larger than {@link Integer#MAX_VALUE
	 *         MAX_VALUE}
	 */
	void put (A values, int offset, int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException;



	/**
	 * Array buffer implementation for <i>64-bit floating-point</i> values.
	 */
	final class OfDouble implements ArrayBuffer<double[]> {
		private double[] buffer = new double[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public double[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final double value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final double[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final double[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>32-bit floating-point</i> values.
	 */
	final class OfFloat implements ArrayBuffer<float[]> {
		private float[] buffer = new float[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public float[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final float value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final float[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final float[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>64-bit integer</i> values.
	 */
	final class OfLong implements ArrayBuffer<long[]> {
		private long[] buffer = new long[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public long[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final long value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final long[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final long[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>32-bit integer</i> values.
	 */
	final class OfInt implements ArrayBuffer<int[]> {
		private int[] buffer = new int[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public int[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final int value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final int[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final int[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>16-bit integer</i> values.
	 */
	final class OfShort implements ArrayBuffer<short[]> {
		private short[] buffer = new short[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public short[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final short value) {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final short[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final short[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>8-bit integer</i> values.
	 */
	final class OfByte implements ArrayBuffer<byte[]> {
		private byte[] buffer = new byte[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public byte[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final byte value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final byte[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final byte[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>character</i> values.
	 */
	final class OfChar implements ArrayBuffer<char[]> {
		private char[] buffer = new char[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public char[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final char value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final char[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final char[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}



	/**
	 * Array buffer implementation for <i>boolean</i> values.
	 */
	final class OfBoolean implements ArrayBuffer<boolean[]> {
		private boolean[] buffer = new boolean[INITIAL_CAPACITY];
		private int position = 0;


		/**
		 * {@inheritDoc}
		 */
		public boolean[] get () {
			return copyOf(this.buffer, this.position);
		}


		/**
		 * Appends the given value to this buffer.
		 * @param value the value
		 * @throws NegativeArraySizeException if appending the given value requires a buffer larger than
		 *         {@link Integer#MAX_VALUE MAX_VALUE}
		 */
		public void put (final boolean value) throws NegativeArraySizeException {
			if (this.position == this.buffer.length) this.buffer = copyOf(this.buffer, (this.buffer.length << 1) | 1);
			this.buffer[this.position++] = value;
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final boolean[] values) throws NullPointerException, NegativeArraySizeException {
			this.put(values, 0, values.length);
		}


		/**
		 * {@inheritDoc}
		 */
		public void put (final boolean[] values, final int offset, final int amount) throws NullPointerException, IllegalArgumentException, NegativeArraySizeException {
			if (offset < 0 | amount < 0 | values.length < offset + amount) throw new IllegalArgumentException();

			if (this.position + amount > this.buffer.length) this.buffer = copyOf(this.buffer, 2 * highestOneBit(this.position + amount) - 1);
			arraycopy(values, offset, this.buffer, this.position, amount);
			this.position += amount;
		}
	}
}