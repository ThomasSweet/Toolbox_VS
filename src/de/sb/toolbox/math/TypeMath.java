package de.sb.toolbox.math;

import static java.lang.Float.NEGATIVE_INFINITY;
import static java.lang.Float.NaN;
import static java.lang.Float.POSITIVE_INFINITY;
import static java.lang.Float.floatToIntBits;
import static java.lang.Float.intBitsToFloat;
import static java.lang.Integer.numberOfLeadingZeros;
import static java.nio.ByteOrder.BIG_ENDIAN;
import java.nio.BufferOverflowException;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.ReadOnlyBufferException;


/**
 * This fassade provides constants and methods useful when dealing with the <i>virtual</i> types {@code int24} and {@code half}.
 * The <i>int24</i> type represents 24-bit two-complement integer values, and is physically represented by 32-bit <i>int</i>
 * values. The <i>half</i> type represents 16-bit floating-point values organized according to the <i>IEEE 754-2008</i>
 * "half-precision" bit layout, and is physically represented by {@code short} values (bit layout), and {@code float} values
 * (content). Note that <i>half-precision</i> values can represent the integer range
 * <tt>[-2<sup>11</sup>,+2<sup>11</sup>] = [-2048,+2048]</tt> without any gaps, i.e. all values within this range can be
 * represented exactly.
 */
public final class TypeMath {

	/**
	 * The number of bits used to represent an {@code int24} value.
	 */
	static public final int INT_24_BYTES = 3;

	/**
	 * The number of bytes used to represent an {@code int24} value.
	 */
	static public final int INT_24_SIZE = INT_24_BYTES * Byte.SIZE;

	/**
	 * A constant holding the {@code int} conversion of the minimum value an {@code int24} can have, -2<sup>23</sup>.
	 */
	static public final int INT_24_MIN_VALUE = -(1 << INT_24_SIZE - 1);

	/**
	 * A constant holding the {@code int} conversion of the maximum value an {@code int24} can have, +2<sup>31</sup>-1.
	 */
	static public final int INT_24_MAX_VALUE = -INT_24_MIN_VALUE - 1;

	/**
	 * The number of bits used to represent a {@code half} value.
	 */
	static public final int HALF_SIZE = Short.SIZE;

	/**
	 * The number of bytes used to represent a {@code half} value.
	 */
	static public final int HALF_BYTES = Short.BYTES;

	/**
	 * A constant holding the {@code float} conversion of the smallest positive normal value of the <i>virtual type</i>
	 * {@code half}, 2<sup>-14</sup>. It is equal to {@code Half.shortBitsToFloat(0x0400)}.
	 */
	static public final float HALF_MIN_NORMAL = 0x1p-14f;

	/**
	 * A constant holding the {@code float} conversion of the smallest positive nonzero value of <i>virtual type</i>
	 * {@code half}, 2<sup>-24</sup>. It is equal to {@code Half.shortBitsToFloat(0x0001)}.
	 */
	static public final float HALF_MIN_VALUE = 0x1p-24f;

	/**
	 * A constant holding the {@code float} conversion of the largest positive finite value of <i>virtual type</i> {@code half},
	 * <tt>(2-2<sup>-10</sup>)&middot;2<sup>15</sup>=65504</tt>. It is equal to {@code Half.shortBitsToHalf(0x7bff)}. Note that
	 * any single-precision value within range <tt>[65504,65536[</tt> that is converted to and from {@code half} will become
	 * this number, while larger values become infinity.
	 */
	static public final float HALF_MAX_VALUE = 0x1.ffcp+15f;

	/**
	 * Maximum exponent a finite {@code half} variable may have.
	 */
	static public final int HALF_MAX_EXPONENT = HALF_SIZE - 1;

	/**
	 * Minimum exponent a normalized {@code half} variable may have.
	 */
	static public final int HALF_MIN_EXPONENT = 2 - HALF_SIZE;

	/**
	 * Minimum exponent a non-zero {@code half} variable may have.
	 */
	static public final int HALF_MIN_SUB_EXPONENT = -HALF_SIZE - 8;

	/**
	 * The bit representation of the special value <tt>+&infin;</tt>.
	 */
	static private final short HALF_POSITIVE_INFINITY_BITS = (short) 0x7c00;

	/**
	 * The <i>canonical</i> bit representation of the special value <tt>NaN</tt>.
	 */
	static private final short HALF_NAN_BITS = (short) 0x7e00;

	/**
	 * The number of bits used to represent a {@code half} value's exponent.
	 */
	static private final int HALF_EXPONENT_SIZE = 5;

	/**
	 * The number of bits used to represent a {@code half} value's significand.
	 */
	static private final int HALF_SIGNIFICAND_SIZE = HALF_SIZE - HALF_EXPONENT_SIZE - 1;

	/**
	 * The bit mask used to extract a {@code half} value's exponent.
	 */
	static private final int HALF_EXPONENT_MASK = (1 << HALF_EXPONENT_SIZE) - 1;

	/**
	 * The bit mask used to extract a {@code half} value's significand.
	 */
	static private final int HALF_SIGNIFICAND_MASK = (1 << HALF_SIGNIFICAND_SIZE) - 1;

	/**
	 * The number of bits used to represent a {@code float} value's exponent.
	 */
	static private final int FLOAT_EXPONENT_SIZE = 8;

	/**
	 * The number of bits used to represent a {@code float} value's significand.
	 */
	static private final int FLOAT_SIGNIFICAND_SIZE = Float.SIZE - FLOAT_EXPONENT_SIZE - 1;

	/**
	 * The bit mask used to extract a {@code float} value's exponent.
	 */
	static private final int FLOAT_EXPONENT_MASK = (1 << FLOAT_EXPONENT_SIZE) - 1;


	/**
	 * Prevents external instantiation.
	 */
	private TypeMath () {}


	//***************************//
	// 24-bit integer operations //
	//***************************//

	/**
	 * Returns a representation of the specified integer value according to a 24-bit two-complement integer bit layout.
	 * @param value a <i>32-bit</i> integer number
	 * @return the bits that represent the corresponding <i>24-bit</i> integer number
	 */
	static public int int24ToIntBits (final int value) {
		return value & 0xffffff;
	}


	/**
	 * Returns the {@code int} value corresponding to a given bit representation. The argument is considered to be a
	 * representation of a 24-bit integer in two-complement bit layout; the leading 8 bits are ignored.
	 * @param int24Bits the bits that represent a <i>24-bit</i> integer number
	 * @return the corresponding <i>32-bit</i> integer number
	 */
	static public int intBitsToInt24 (final int int24Bits) {
		return int24Bits << Byte.SIZE >> Byte.SIZE;
	}


	/**
	 * Reads the next three bytes at the given buffer's current position, composing them into an int value according to the
	 * current byte order, and then increments the position by three.
	 * @param buffer the buffer
	 * @return the <i>24-bit integer</i> value at the buffer's current position
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws BufferUnderflowException if there are fewer than three bytes remaining in the given buffer
	 */
	static public int getInt24 (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
		return buffer.order() == BIG_ENDIAN
			? (buffer.get() << 16) | ((buffer.get() & 0xff) << 8) | (buffer.get() & 0xff)
			: (buffer.get() & 0xff) | ((buffer.get() & 0xff) << 8) | (buffer.get() << 16);
	}


	/**
	 * Writes three bytes containing the 24 low-bits of the given int value, in the current byte order, into the given buffer at
	 * it's current position, and then increments the position by three.
	 * @param buffer the buffer
	 * @param value the <i>24-bit integer</i> value (the 8 high-bits are ignored)
	 * @throws NullPointerException if the given buffer is {@code null}
	 * @throws BufferOverflowException if there are fewer than three bytes remaining in the given buffer
	 * @throws ReadOnlyBufferException if the given buffer is read-only
	 */
	static public void putInt24 (final ByteBuffer buffer, final int value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
		if (buffer.order() == BIG_ENDIAN) {
			buffer.put((byte) (value >>> 16));
			buffer.put((byte) (value >>> 8));
			buffer.put((byte) value);
		} else {
			buffer.put((byte) value);
			buffer.put((byte) (value >>> 8));
			buffer.put((byte) (value >>> 16));
		}
	}


	//**********************************//
	// 16-bit half-precision operations //
	//**********************************//

	/**
	 * Returns {@code true} if the argument represents a finite floating-point value; returns {@code false} otherwise (for NaN
	 * and infinity arguments).
	 * @param halfBits the bits that represent a <i>half-precision</i> floating-point number
	 * @return {@code true} if the argument is a finite floating-point value, {@code false} otherwise.
	 */
	static public boolean isFiniteHalf (final short halfBits) {
		return halfExponent(halfBits) != HALF_MAX_EXPONENT + 1;
	}


	/**
	 * Returns {@code true} if the argument represents an infinite floating-point value, {@code false} otherwise.
	 * @param halfBits the bits that represent a <i>half-precision</i> floating-point number
	 * @return {@code true} if the argument is <tt>&plusmn;&infin;</tt>; {@code false} otherwise.
	 */
	static public boolean isInfiniteHalf (final short halfBits) {
		return halfBits == HALF_POSITIVE_INFINITY_BITS | halfBits == (HALF_POSITIVE_INFINITY_BITS | Short.MIN_VALUE);
	}


	/**
	 * Returns {@code true} if the argument represents a Not-a-Number (NaN) value, {@code false} otherwise.
	 * @param halfBits the bits that represent a <i>half-precision</i> floating-point number
	 * @return {@code true} if the argument is {@code NaN}; {@code false} otherwise.
	 */
	static public boolean isNanHalf (final short halfBits) {
		return !isFiniteHalf(halfBits) & !isInfiniteHalf(halfBits);
	}


	/**
	 * Returns the absolute value of a {@code half} value. If the argument is not negative, the argument is returned. If the
	 * argument is negative, the negation of the argument is returned. Special cases:
	 * <ul>
	 * <li>If the argument is positive zero or negative zero, the result is positive zero.
	 * <li>If the argument is infinite, the result is positive infinity.
	 * <li>If the argument is NaN, the result is NaN.
	 * </ul>
	 * @param halfBits the bits that represent the <i>half-precision</i> floating-point number <tt>x</tt>
	 * @return the bits that represent the <i>half-precision</i> floating-point number <tt>|x|</tt>
	 */
	static public short halfAbs (final short halfBits) {
		return (short) (halfBits & Short.MAX_VALUE);
	}


	/**
	 * Returns the signum function of the argument; <tt>&plusmn;0h</tt> if the argument represents zero, {@code +1h} if the
	 * argument represents a value greater than zero, {@code -1h} if the argument represents a value less than zero. Special
	 * Cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is <i>canonical</i> NaN.
	 * <li>If the argument is positive zero or negative zero, then the result is the same as the argument.
	 * </ul>
	 * @param halfBits the bits that represent the <i>half-precision</i> floating-point number <tt>x</tt>
	 * @return the bits that represent the <i>half-precision</i> floating-point number <tt>sgn(x)</tt>
	 */
	static public short halfSignum (final short halfBits) {
		if (halfBits == 0 | halfBits == Short.MIN_VALUE) return halfBits;
		if (isNanHalf(halfBits)) return HALF_NAN_BITS;
		return (halfBits >>> Integer.SIZE - 1) == 1 ? (short) 0xbc00 : 0x3c00;
	}


	/**
	 * Returns the unbiased exponent used in the representation of a {@code half}. Special cases:
	 * <ul>
	 * <li>If the argument is NaN or infinite, then the result is {@link #HALF_MAX_EXPONENT} + 1.
	 * <li>If the argument is zero or subnormal, then the result is {@link #HALF_MIN_EXPONENT} -1.
	 * </ul>
	 * @param halfBits the bits that represent a <i>half-precision</i> floating-point number
	 * @return the unbiased exponent of the argument
	 */
	static public int halfExponent (final short halfBits) {
		return ((halfBits >>> HALF_SIGNIFICAND_SIZE) & HALF_EXPONENT_MASK) - HALF_MAX_EXPONENT;
	}


	/**
	 * Returns a representation of the specified floating-point value according to the <i>IEEE 754-2008</i> floating-point
	 * "half format" bit layout
	 * <p>
	 * Bit 15 (the bit that is selected by the mask {@code 0x8000}) represents the sign of the floating-point number. Bits 14-10
	 * (the bits that are selected by the mask {@code 0x7c00}) represent the exponent. Bits 9-0 (the bits that are selected by
	 * the mask {@code 0x03ff}) represent the significand (sometimes called the mantissa) of the floating-point number.
	 * <p>
	 * If the argument is or exceeds <tt>+2<sup>16</sup></tt> (including <tt>+&infin;</tt>), the result is {@code 0x7c00}.
	 * <p>
	 * If the argument is or deceeds <tt>-2<sup>16</sup></tt> (including <tt>-&infin;</tt>), the result is {@code 0xfc00}.
	 * <p>
	 * If the argument is {@code NaN}, the result is {@code 0x7e00}.
	 * <p>
	 * In all cases, the result is a short that, when given to the {@link #shortBitsToHalf(short)} method, will produce a
	 * floating-point value that is the same as the argument to {@code floatToShortBits} (except all NaN values are collapsed to
	 * a single <i>canonical</i> NaN value).
	 * @param value a <i>single-precision</i> floating-point number
	 * @return the bits that represent the corresponding <i>half-precision</i> floating-point number
	 */
	static public short halfToShortBits (final float value) {
		if (Float.isNaN(value)) return HALF_NAN_BITS;

		final int intBits = floatToIntBits(value);
		final int sign = (intBits & Integer.MIN_VALUE) >>> HALF_SIZE;
		final int significand = (intBits >>> FLOAT_SIGNIFICAND_SIZE - HALF_SIGNIFICAND_SIZE) & HALF_SIGNIFICAND_MASK;
		final int exponent = ((intBits >>> FLOAT_SIGNIFICAND_SIZE) & FLOAT_EXPONENT_MASK) - Float.MAX_EXPONENT;
		if (exponent < HALF_MIN_SUB_EXPONENT) return (short) sign;
		if (exponent < HALF_MIN_EXPONENT) return (short) (sign | (1 << exponent - HALF_MIN_SUB_EXPONENT) | (significand >>> HALF_MIN_EXPONENT - exponent));
		if (exponent <= HALF_MAX_EXPONENT) return (short) (sign | (exponent + HALF_MAX_EXPONENT << HALF_SIGNIFICAND_SIZE) | significand);
		return (short) (sign | HALF_POSITIVE_INFINITY_BITS);
	}


	/**
	 * Returns the {@code float} value corresponding to a given bit representation. The argument is considered to be a
	 * representation of a floating-point value according to the <i>IEEE 754-2008</i> floating-point "half format" bit layout.
	 * <p />
	 * If the argument is {@code 0x7c00}, the result is positive infinity.
	 * <p />
	 * If the argument is {@code 0xfc00}, the result is negative infinity.
	 * <p />
	 * If the argument is any value in the range {@code 0x7c01} through {@code 0x7fff} or in the range {@code 0xfc01} through
	 * {@code 0xffff}, the result is a <i>canonical</i> NaN.
	 * <p />
	 * In all other cases, let <i>s</i>, <i>e</i>, and <i>m</i> be three values that can be computed from the argument:<pre>
	 * int s = (((bits & 0xffff) >> 15) == 0) ? +1 : -1;
	 * int e = (bits >> 10) & 0x1f;
	 * long m = (e == 0) ? (bits & 0x3ff) << 1 : (bits & 0x3ff) | 0x400;</pre>Then the floating-point result equals the value of
	 * the mathematical expression <i>s</i>&middot; <i>m</i>&middot;2<sup><i>e</i>-25</sup>.
	 * @param halfBits the bits that represent a <i>half-precision</i> floating-point number
	 * @return the corresponding <i>single-precision</i> floating-point number
	 */
	static public float shortBitsToHalf (final short halfBits) {
		switch (halfBits) {
			case 0:
				return +0f;
			case Short.MIN_VALUE:
				return -0f;
			case HALF_POSITIVE_INFINITY_BITS:
				return POSITIVE_INFINITY;
			case HALF_POSITIVE_INFINITY_BITS | Short.MIN_VALUE:
				return NEGATIVE_INFINITY;
			default: {
				final int intBits = halfBits;
				int exponent = ((intBits >>> HALF_SIGNIFICAND_SIZE) & HALF_EXPONENT_MASK) + (Float.MAX_EXPONENT - HALF_MAX_EXPONENT);
				int significand = intBits & HALF_SIGNIFICAND_MASK;
				switch (exponent) {
					case Float.MAX_EXPONENT + HALF_MAX_EXPONENT + 1:
						return NaN;
					case Float.MAX_EXPONENT + HALF_MIN_EXPONENT - 1:
						final int shift = (HALF_SIGNIFICAND_SIZE - Float.SIZE + 1) + numberOfLeadingZeros(significand);
						exponent += 1 - shift;
						significand = (significand << shift) ^ (1 << HALF_SIGNIFICAND_SIZE);
					default:
						return intBitsToFloat((intBits & Integer.MIN_VALUE) | (exponent << FLOAT_SIGNIFICAND_SIZE) | (significand << FLOAT_SIGNIFICAND_SIZE - HALF_SIGNIFICAND_SIZE));
				}
			}
		}
	}


	/**
	 * Reads the next two bytes at the given buffer's current position, composing them into a <i>half-precision
	 * floating-point</i> value according to the current byte order, and then increments the position by two.
	 * @param buffer the buffer
	 * @return the <i>half-precision floating-point</i> value at the buffer's current position
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws BufferUnderflowException if there are fewer than two bytes remaining in the given buffer
	 */
	static public float getHalf (final ByteBuffer buffer) throws NullPointerException, BufferUnderflowException {
		return shortBitsToHalf(buffer.getShort());
	}


	/**
	 * Writes two bytes containing the <i>half-precision</i> bit representation of the given float value, in the current byte
	 * order, into the given buffer at it's current position, and then increments the position by two.
	 * @param buffer the buffer
	 * @param value the <i>half-precision floating-point</i> value
	 * @throws NullPointerException if the given buffer is {@code null}
	 * @throws BufferOverflowException if there are fewer than two bytes remaining in the given buffer
	 * @throws ReadOnlyBufferException if the given buffer is read-only
	 */
	static public void putHalf (final ByteBuffer buffer, final float value) throws NullPointerException, BufferOverflowException, ReadOnlyBufferException {
		buffer.putShort(halfToShortBits(value));
	}
}