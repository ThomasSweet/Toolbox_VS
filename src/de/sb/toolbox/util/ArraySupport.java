package de.sb.toolbox.util;

import static java.lang.Math.max;
import static java.lang.reflect.Array.newInstance;
import static de.sb.toolbox.util.Iterables.newList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.PrimitiveIterator;
import de.sb.toolbox.Copyright;


/**
 * This facade provides array operations for element iterating, swapping, braiding, and unbraiding.
 */
@Copyright(year = 2017, holders = "Sascha Baumeister")
public final class ArraySupport {

	/**
	 * Prevents external instantiation.
	 */
	private ArraySupport () {}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final boolean[] leftVector, final int leftIndex, final boolean[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final boolean value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final char[] leftVector, final int leftIndex, final char[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final char value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final byte[] leftVector, final int leftIndex, final byte[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final byte value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final short[] leftVector, final int leftIndex, final short[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final short value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final int[] leftVector, final int leftIndex, final int[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final int value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final long[] leftVector, final int leftIndex, final long[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final long value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final float[] leftVector, final int leftIndex, final float[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final float value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public void swap (final double[] leftVector, final int leftIndex, final double[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final double value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * Swaps the position of the given elements within the given vectors. Note that both the given vectors and the given indices
	 * may be the same.
	 * @param <T> the vector component type
	 * @param leftVector the left vector operand
	 * @param leftIndex the left index
	 * @param rightVector the right vector operand
	 * @param rightIndex the right index
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayIndexOutOfBoundsException if any of the given indices is out of bounds
	 */
	static public <T> void swap (final T[] leftVector, final int leftIndex, final T[] rightVector, final int rightIndex) throws NullPointerException, ArrayIndexOutOfBoundsException {
		final T value = leftVector[leftIndex];
		leftVector[leftIndex] = rightVector[rightIndex];
		rightVector[rightIndex] = value;
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final boolean[] left, final boolean[] right, final boolean[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final char[] left, final char[] right, final char[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final byte[] left, final byte[] right, final byte[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final short[] left, final short[] right, final short[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final int[] left, final int[] right, final int[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final long[] left, final long[] right, final long[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final float[] left, final float[] right, final float[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void braid (final double[] left, final double[] right, final double[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Braids} the elements of the given left and right operands by alternately copying their elements into the given
	 * braid vector in sequence.
	 * @param <T> the vector component type
	 * @param left the left vector operand, containing the {@code even} indexed elements for the braided vector
	 * @param right the right vector operand, containing the {@code odd} indexed elements for the braided vector
	 * @param braid the braided result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public <T> void braid (final T[] left, final T[] right, final T[] braid) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			braid[index] = left[index >> 1];
			braid[index + 1] = right[index >> 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final boolean[] braid, final boolean[] left, final boolean[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final char[] braid, final char[] left, final char[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final byte[] braid, final byte[] left, final byte[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final short[] braid, final short[] left, final short[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final int[] braid, final int[] left, final int[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final long[] braid, final long[] left, final long[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final float[] braid, final float[] left, final float[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public void unbraid (final double[] braid, final double[] left, final double[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * {@code Unbraides} the elements of the given braided vector operand by copying it's even indexed elements into the given
	 * left result vector, and it's odd indexed elements into the given right result vector.
	 * @param <T> the vector component type
	 * @param braid the braided vector operand, with even indexed elements for the left vector, and odd indexed elements for the
	 *        right vector
	 * @param left the unbraided left result vector
	 * @param right the unbraided right result vector
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given left and right vectors don't share the same length, or if the given braid
	 *         vector isn't double their size
	 */
	static public <T> void unbraid (final T[] braid, final T[] left, final T[] right) throws NullPointerException, IllegalArgumentException {
		if (left.length != right.length | left.length + right.length != braid.length) throw new IllegalArgumentException();

		for (int index = 0; index < braid.length; index += 2) {
			left[index >> 1] = braid[index];
			right[index >> 1] = braid[index + 1];
		}
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public boolean[][] wrap (final boolean... vector) throws NullPointerException {
		final boolean[][] result = new boolean[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public char[][] wrap (final char... vector) throws NullPointerException {
		final char[][] result = new char[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public byte[][] wrap (final byte... vector) throws NullPointerException {
		final byte[][] result = new byte[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public short[][] wrap (final short... vector) throws NullPointerException {
		final short[][] result = new short[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public int[][] wrap (final int... vector) throws NullPointerException {
		final int[][] result = new int[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public long[][] wrap (final long... vector) throws NullPointerException {
		final long[][] result = new long[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public float[][] wrap (final float... vector) throws NullPointerException {
		final float[][] result = new float[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	static public double[][] wrap (final double... vector) throws NullPointerException {
		final double[][] result = new double[vector.length][1];
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Wraps each coordinate of the given vector into it's own single-coordinate vector, and returns the resulting matrix.
	 * @param <T> the vector component type
	 * @param vector the operand vector
	 * @return the matrix created
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	@SafeVarargs
	static public <T> T[][] wrap (final T... vector) throws NullPointerException {
		@SuppressWarnings("unchecked")
		final T[][] result = (T[][]) newInstance(vector.getClass().getComponentType(), vector.length, 1);
		for (int index = 0; index < vector.length; ++index)
			result[index][0] = vector[index];
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public boolean[] unwrap (final boolean[]... matrix) throws NullPointerException, IllegalArgumentException {
		final boolean[] result = new boolean[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final boolean[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public char[] unwrap (final char[]... matrix) throws NullPointerException, IllegalArgumentException {
		final char[] result = new char[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final char[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public byte[] unwrap (final byte[]... matrix) throws NullPointerException, IllegalArgumentException {
		final byte[] result = new byte[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final byte[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public short[] unwrap (final short[]... matrix) throws NullPointerException, IllegalArgumentException {
		final short[] result = new short[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final short[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public int[] unwrap (final int[]... matrix) throws NullPointerException, IllegalArgumentException {
		final int[] result = new int[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final int[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public long[] unwrap (final long[]... matrix) throws NullPointerException, IllegalArgumentException {
		final long[] result = new long[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final long[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public float[] unwrap (final float[]... matrix) throws NullPointerException, IllegalArgumentException {
		final float[] result = new float[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final float[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	static public double[] unwrap (final double[]... matrix) throws NullPointerException, IllegalArgumentException {
		final double[] result = new double[matrix.length];
		for (int index = 0; index < matrix.length; ++index) {
			final double[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Unwraps each single-coordinate element of the given matrix, and returns the resulting vector.
	 * @param <T> the vector component type
	 * @param matrix the operand matrix
	 * @return the vector created
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if any of the given matrix elements has a non-unary length
	 */
	@SafeVarargs
	static public <T> T[] unwrap (final T[]... matrix) throws NullPointerException, IllegalArgumentException {
		@SuppressWarnings("unchecked")
		final Class<T> componentType = (Class<T>) matrix.getClass().getComponentType().getComponentType();

		@SuppressWarnings("unchecked")
		final T[] result = (T[]) newInstance(componentType, matrix.length);
		for (int index = 0; index < matrix.length; ++index) {
			final T[] element = matrix[index];
			if (element.length != 1) throw new IllegalArgumentException();
			result[index] = element[0];
		}
		return result;
	}


	/**
	 * Returns the given iterator's elements as an array of the given component type.
	 * @param iterator the iterator
	 * @param componentType the resulting array's component type
	 * @return the array created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws ArrayStoreException if any iterator element is incompatible with the given component type
	 */
	static public <T> T[] newArray (final Class<T> componentType, final Iterator<? extends T> iterator) throws NullPointerException, ArrayStoreException {
		final List<? extends T> list = newList(iterator);

		@SuppressWarnings("unchecked")
		final T[] array = list.toArray((T[]) newInstance(componentType, list.size()));
		return array;
	}


	/**
	 * Returns the given array wrapped into a modifiable list iterator.
	 * @param array the array
	 * @return the list iterator adapter created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T> ListIterator<T> newListIterator (final T[] array) {
		return new ListIterator<T>() {
			private int position = -1;


			public boolean hasNext () {
				return this.position + 1 < array.length;
			}


			public int nextIndex () {
				return this.position + 1;
			}


			public T next () throws NoSuchElementException {
				if (this.hasNext()) return array[++this.position];
				throw new NoSuchElementException();
			}


			public int previousIndex () {
				return max(this.position - 1, -1);
			}


			public boolean hasPrevious () {
				return this.position > 0;
			}


			public T previous () throws NoSuchElementException {
				if (this.hasPrevious()) return array[--this.position];
				throw new NoSuchElementException();
			}


			public void remove () throws UnsupportedOperationException {
				throw new UnsupportedOperationException();
			}


			public void set (T element) throws IllegalStateException {
				try {
					array[this.position] = element;
				} catch (final ArrayIndexOutOfBoundsException exception) {
					throw new IllegalStateException(exception);
				}
			}


			public void add (final T element) throws UnsupportedOperationException {
				throw new UnsupportedOperationException();
			}
		};
	}


	/**
	 * Returns the given array wrapped into an unmodifiable iterator.
	 * @param array the array
	 * @return the iterator adapter created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public <T> Iterator<T> newIterator (final T[] array) {
		return new Iterator<T>() {
			private int position;


			public boolean hasNext () {
				return this.position != array.length;
			}


			public T next () {
				if (this.hasNext()) return array[this.position++];
				throw new NoSuchElementException();
			}
		};
	}


	/**
	 * Returns the given array wrapped into an unmodifiable int iterator.
	 * @param array the array
	 * @return the iterator adapter created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public PrimitiveIterator.OfInt newIterator (final int[] array) {
		return new PrimitiveIterator.OfInt() {
			private int position;


			public boolean hasNext () {
				return this.position != array.length;
			}


			public int nextInt () {
				if (this.hasNext()) return array[this.position++];
				throw new NoSuchElementException();
			}
		};
	}


	/**
	 * Returns the given array wrapped into an unmodifiable long iterator.
	 * @param array the array
	 * @return the iterator adapter created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public PrimitiveIterator.OfLong newIterator (final long[] array) {
		return new PrimitiveIterator.OfLong() {
			private int position;


			public boolean hasNext () {
				return this.position != array.length;
			}


			public long nextLong () {
				if (this.hasNext()) return array[this.position++];
				throw new NoSuchElementException();
			}
		};
	}


	/**
	 * Returns the given array wrapped into an unmodifiable double iterator.
	 * @param array the array
	 * @return the iterator adapter created
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	static public PrimitiveIterator.OfDouble newIterator (final double[] array) {
		return new PrimitiveIterator.OfDouble() {
			private int position;


			public boolean hasNext () {
				return this.position != array.length;
			}


			public double nextDouble () {
				if (this.hasNext()) return array[this.position++];
				throw new NoSuchElementException();
			}
		};
	}
}