package de.sb.toolbox.validation;

import java.util.function.Predicate;
import de.sb.toolbox.Copyright;


/**
 * Simple predicate validator for texts representing floating-point numbers within a defined range.
 */
@Copyright(year = 2016, holders = "Sascha Baumeister")
public class DoubleValidator implements Predicate<String> {

	private final double min, max;


	/**
	 * Creates an unbounded validator.
	 */
	public DoubleValidator () {
		this(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
	}


	/**
	 * Creates a validator that accepts values between min and max.
	 * @param min the min value (included)
	 * @param max the max value (included)
	 * @throws IllegalArgumentException if min or max is NaN, or if min exceeds max
	 */
	public DoubleValidator (final double min, final double max) {
		if (Double.isNaN(min) | Double.isNaN(max) | min > max) throw new IllegalArgumentException();

		this.min = min;
		this.max = max;
	}


	/**
	 * Returns the min value.
	 * @return the allowed minimum for parsed values
	 */
	public double getMin () {
		return this.min;
	}


	/**
	 * Returns the max value.
	 * @return the allowed maximum for parsed values
	 */
	public double getMax () {
		return this.max;
	}


	/**
	 * {@inheritDoc}
	 * @param text the (optional) text to be validated
	 * @return whether or not the given text validates to a valid double value within the given range
	 */
	@Override
	public boolean test (final String text) {
		try {
			final double value = Double.parseDouble(text);
			return value >= this.min & value <= this.max;
		} catch (final Exception exception) {
			return false;
		}
	}
}
