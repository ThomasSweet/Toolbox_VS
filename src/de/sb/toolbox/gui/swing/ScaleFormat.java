package de.sb.toolbox.gui.swing;

import java.math.BigDecimal;
import java.text.FieldPosition;
import java.text.Format;
import java.text.NumberFormat;
import java.text.ParsePosition;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.Reflections;


/**
 * Convenience wrapper class for {@link NumberFormat} instances the shifts the decimal scale of values by a given amount. For
 * example, this is useful for business applications that handle price information as integral (denominator) values to avoid
 * floating point rounding problems, but need to display or edit the information in nominal currency.
 */
@Copyright(year = 2012, holders = "Sascha Baumeister")
public class ScaleFormat extends Format {
	static private final long serialVersionUID = 1L;

	private final NumberFormat delegate;
	private final int scale;


	/**
	 * Creates a new instance based on the given delegate. Note that the given scale may be negative.
	 * @param delegate the underlying {@link NumberFormat}
	 * @param scale the decimal left shift for formatting, and the decimal right shift during parsing
	 * @throws NullPointerException if the given delegate is {@code null}
	 */
	public ScaleFormat (final NumberFormat delegate, final int scale) {
		if (delegate == null) throw new NullPointerException();

		this.delegate = delegate;
		this.scale = scale;
	}


	/**
	 * {@inheritDoc}
	 * @throws NullPointerException if buffer or position is {@code null}
	 * @throws IllegalArgumentException if the Format cannot format the given object
	 */
	@Override
	public StringBuffer format (final Object object, final StringBuffer buffer, final FieldPosition position) {
		if (!(object instanceof Number)) throw new IllegalArgumentException("Cannot format given Object as a Number");

		final BigDecimal value = Reflections.cast((Number) object, BigDecimal.class);
		final BigDecimal scaledValue = value.scaleByPowerOfTen(this.scale);
		return this.delegate.format(scaledValue, buffer, position);
	}


	/**
	 * {@inheritDoc}
	 * @throws NullPointerException if position is {@code null}
	 */
	@Override
	public Object parseObject (final String source, final ParsePosition position) {
		final Number number = this.delegate.parse(source, position);
		if (number == null) return null;

		final BigDecimal decimal = Reflections.cast(number, BigDecimal.class);
		final BigDecimal scaledDecimal = decimal.scaleByPowerOfTen(-this.scale);
		try {
			return scaledDecimal.longValueExact();
		} catch (final ArithmeticException exception) {
			return scaledDecimal.doubleValue();
		}
	}
}