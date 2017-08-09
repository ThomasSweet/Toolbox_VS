package de.sb.toolbox.math;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.scalb;
import static de.sb.toolbox.math.TypeMath.HALF_BYTES;
import static de.sb.toolbox.math.TypeMath.INT_24_BYTES;
import static de.sb.toolbox.math.TypeMath.INT_24_MAX_VALUE;
import static de.sb.toolbox.math.TypeMath.INT_24_MIN_VALUE;
import static de.sb.toolbox.math.TypeMath.INT_24_SIZE;
import static de.sb.toolbox.math.TypeMath.getHalf;
import static de.sb.toolbox.math.TypeMath.getInt24;
import static de.sb.toolbox.math.TypeMath.putHalf;
import static de.sb.toolbox.math.TypeMath.putInt24;
import java.nio.BufferOverflowException;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.ReadOnlyBufferException;
import de.sb.toolbox.Copyright;


/**
 * Instances of this enum marshal/unmarshal <i>single-precision</i> floating-point values normed to range <tt>[-1,+1]</tt>
 * to/from a variety of binary representations, both integer (signed/unsigned) and floating-point based, and both for
 * <i>little-endian</i> and <i>big-endian</i> byte order. The marshaled values are handled by a <i>byte buffer</i> instance
 * which already provides many of the conversions required.
 */
@Copyright(year = 2017, holders = "Sascha Baumeister")
public enum FloatMarshaler {

	/**
	 * 8-bit signed integer marshaler. The marshaled range is <tt>[-2<sup>7</sup>,+2<sup>7</sup>[</tt>, while the unmarshaled
	 * (norm) range is clipped to <tt>[-1,+1]</tt>.
	 */
	INT_8 (Byte.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(buffer.get(), 1 - Byte.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.put((byte) round(max(Byte.MIN_VALUE, min(Byte.MAX_VALUE, scalb(value, Byte.SIZE - 1)))));
		}
	},

	/**
	 * 8-bit unsigned integer marshaler. The marshaled range is <tt>[0,+2<sup>8</sup>[</tt>, while the unmarshaled (norm) range
	 * is clipped to <tt>[-1,+1]</tt>.
	 */
	UINT_8 (Byte.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(buffer.get() - Byte.MIN_VALUE, 1 - Byte.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.put((byte) round((Byte.MAX_VALUE + .5f) * max(0, min(2, value + 1))));
		}
	},

	/**
	 * 16-bit signed integer marshaler. The marshaled range is <tt>[-2<sup>15</sup>,+2<sup>15</sup>[</tt>, while the unmarshaled
	 * (norm) range is clipped to <tt>[-1,+1]</tt>.
	 */
	INT_16 (Short.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(buffer.getShort(), 1 - Short.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putShort((short) round(max(Short.MIN_VALUE, min(Short.MAX_VALUE, scalb(value, Short.SIZE - 1)))));
		}
	},

	/**
	 * 16-bit unsigned integer marshaler. The marshaled range is <tt>[0,+2<sup>16</sup>[</tt>, while the unmarshaled (norm)
	 * range is clipped to <tt>[-1,+1]</tt>.
	 */
	UINT_16 (Short.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(buffer.getShort() - Short.MIN_VALUE, 1 - Short.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putShort((short) round((Short.MAX_VALUE + .5f) * max(0, min(2, value + 1))));
		}
	},

	/**
	 * 24-bit signed integer marshaler. The marshaled range is <tt>[-2<sup>23</sup>,+2<sup>23</sup>[</tt>, while the unmarshaled
	 * (norm) range is clipped to <tt>[-1,+1]</tt>.
	 */
	INT_24 (INT_24_BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(getInt24(buffer), 1 - INT_24_SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			putInt24(buffer, round(max(INT_24_MIN_VALUE, min(INT_24_MAX_VALUE, scalb(value, INT_24_SIZE - 1)))));
		}
	},

	/**
	 * 24-bit unsigned integer marshaler. The marshaled range is <tt>[0,+2<sup>24</sup>[</tt>, while the unmarshaled (norm)
	 * range is clipped to <tt>[-1,+1]</tt>.
	 */
	UINT_24 (INT_24_BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(getInt24(buffer) - INT_24_MIN_VALUE, 1 - INT_24_SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putInt(round((INT_24_MAX_VALUE + .5f) * max(0, min(2, value + 1))));
		}
	},

	/**
	 * 32-bit signed integer marshaler. The marshaled range is <tt>[-2<sup>31</sup>,+2<sup>31</sup>[</tt>, while the unmarshaled
	 * (norm) range is clipped to <tt>[-1,+1]</tt>.
	 */
	INT_32 (Integer.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb((float) buffer.getInt(), 1 - Integer.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putInt(round(max(Integer.MIN_VALUE, min(Integer.MAX_VALUE, scalb(value, Integer.SIZE - 1)))));
		}
	},

	/**
	 * 32-bit unsigned integer marshaler. The marshaled range is <tt>[0,+2<sup>32</sup>[</tt>, while the unmarshaled (norm)
	 * range is clipped to <tt>[-1,+1]</tt>.
	 */
	UINT_32 (Integer.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return scalb(buffer.getInt() - (long) Integer.MIN_VALUE, 1 - Integer.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putInt(round((Integer.MAX_VALUE + .5f) * max(0, min(2, value + 1))));
		}
	},

	/**
	 * 64-bit signed integer marshaler. The marshaled range is <tt>[-2<sup>63</sup>,+2<sup>63</sup>[</tt>, while the unmarshaled
	 * (norm) range is clipped to <tt>[-1,+1]</tt>.
	 */
	INT_64 (Long.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return (float) scalb((double) buffer.getLong(), 1 - Long.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putLong(round(max(Long.MIN_VALUE, min(Long.MAX_VALUE, scalb((double) value, Long.SIZE - 1)))));
		}
	},

	/**
	 * 64-bit unsigned integer marshaler. The marshaled range is <tt>[0,+2<sup>64</sup>[</tt>, while the unmarshaled (norm)
	 * range is clipped to <tt>[-1,+1]</tt>.
	 */
	UINT_64 (Long.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return (float) scalb(buffer.getLong() - (double) Long.MIN_VALUE, 1 - Long.SIZE);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putLong(round((Long.MAX_VALUE + .5) * max(0, min(2, value + 1))));
		}
	},

	/**
	 * 16-bit floating-point (IEEE 754-2008 "half precision") marshaler. Neither the unmarshaled nor the marshaled range is
	 * clipped, but values with absolutes equaling or exceeding <tt>2<sup>16</sup></tt> will become <tt>&plusmn;&infin;</tt>.
	 */
	HALF (HALF_BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return getHalf(buffer);
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			putHalf(buffer, value);
		}
	},

	/**
	 * 32-bit floating-point (IEEE 754 "single precision") marshaler. Neither the unmarshaled nor the marshaled range is
	 * clipped, but values with absolutes equaling or exceeding <tt>2<sup>128</sup></tt> will become <tt>&plusmn;&infin;</tt>.
	 */
	SINGLE (Float.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return buffer.getFloat();
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putFloat(value);
		}
	},

	/**
	 * 64-bit floating-point (IEEE 754 "double precision") marshaler. Neither the unmarshaled nor the marshaled range is
	 * clipped, and values with absolutes equaling or exceeding <tt>2<sup>1024</sup></tt> are already <tt>&plusmn;&infin;</tt>.
	 */
	DOUBLE (Double.BYTES) {
		public float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
			return (float) buffer.getDouble();
		}

		public void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
			buffer.putDouble(value);
		}
	};


	private final int valueSize;


	/**
	 * Creates a new instance with the given value size.
	 * @param valueSize the size of marshaled values in bytes
	 */
	private FloatMarshaler (final int valueSize) {
		this.valueSize = valueSize;
	}


	/**
	 * Returns the value size.
	 * @return the size of marshaled values in bytes
	 */
	public int valueSize () {
		return valueSize;
	}


	/**
	 * Reads a number of bytes from the given buffer, in it's current byte order and at it's current position, incrementing the
	 * latter; <i>unmarshals</i> the resulting bytes and returns the resulting value.
	 * @param buffer the byte buffer
	 * @return the unmarshaled value read
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws BufferUnderflowException if there are too few bytes remaining in the given buffer
	 */
	public abstract float get (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException;


	/**
	 * <i>Marshals</i> the given value and writes the resulting bytes into the given buffer, in it's current byte order and at
	 * it's current position, incrementing the latter. Note that the given value is clipped to norm range <tt>[-1,+1]</tt> if
	 * this marshaler creates integer values.
	 * @param buffer the byte buffer
	 * @param value the value to be marshaled and written
	 * @throws NullPointerException if the given buffer is {@code null}
	 * @throws BufferOverflowException if there are too few bytes remaining in the given buffer
	 * @throws ReadOnlyBufferException if the given buffer is read-only
	 */
	public abstract void put (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException;


	/**
	 * Returns the marshaler matching the given value properties.
	 * @param floatingPoint {@code true} if the marshaled values are <i>floating-point</i> numbers, {@code false} if they are
	 *        <i>integer</i> numbers
	 * @param unsigned {@code true} if the marshaled values are <i>unsigned</i>, {@code false} if they are <i>signed</i> (only
	 *        relevant for integer values)
	 * @param valueSize the size of marshaled values in bytes
	 * @return the matching marshaler
	 * @throws IllegalArgumentException if the given property combination is not supported
	 */
	static public FloatMarshaler valueOf (final boolean floatingPoint, final boolean unsigned, final int valueSize) throws IllegalArgumentException {
		if (floatingPoint) {
			switch (valueSize) {
				case HALF_BYTES:
					return HALF;
				case Float.BYTES:
					return SINGLE;
				case Double.BYTES:
					return DOUBLE;
				default:
					break;
			}
		} else {
			if (unsigned) {
				switch (valueSize) {
					case Byte.BYTES:
						return UINT_8;
					case Short.BYTES:
						return UINT_16;
					case INT_24_BYTES:
						return UINT_24;
					case Integer.BYTES:
						return UINT_32;
					case Long.BYTES:
						return UINT_64;
					default:
						break;
				}
			} else {
				switch (valueSize) {
					case Byte.BYTES:
						return INT_8;
					case Short.BYTES:
						return INT_16;
					case INT_24_BYTES:
						return INT_24;
					case Integer.BYTES:
						return INT_32;
					case Long.BYTES:
						return INT_64;
					default:
						break;
				}
			}
		}

		throw new IllegalArgumentException();
	}
}