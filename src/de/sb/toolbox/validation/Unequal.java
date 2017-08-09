package de.sb.toolbox.validation;

import java.lang.annotation.ElementType;
import java.lang.annotation.Repeatable;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.util.Objects;
import javax.validation.Constraint;
import javax.validation.ConstraintValidator;
import javax.validation.ConstraintValidatorContext;
import javax.validation.Payload;
import de.sb.toolbox.Copyright;


/**
 * The annotated element must have a text representation that differs from the annotation's value, which is always the case with
 * {@code null}. Accepts any type.
 * @see Object#toString()
 */
@Target({ ElementType.FIELD, ElementType.METHOD, ElementType.CONSTRUCTOR, ElementType.PARAMETER, ElementType.ANNOTATION_TYPE })
@Retention(RetentionPolicy.RUNTIME)
@Repeatable(Unequal.List.class)
@Constraint(validatedBy = Unequal.Validator.class)
@Copyright(year = 2015, holders = "Sascha Baumeister")
public @interface Unequal {

	/**
	 * Allows several {@link Unequal} annotations on the same element.
	 */
	@Target({ ElementType.FIELD, ElementType.METHOD, ElementType.CONSTRUCTOR, ElementType.PARAMETER, ElementType.ANNOTATION_TYPE })
	@Retention(RetentionPolicy.RUNTIME)
	static @interface List {
		Unequal[]value();
	}


	/**
	 * @return the value
	 */
	String value();


	String message() default "{de.sb.toolbox.validation.Unequal.message}";


	Class<?>[]groups() default {};


	Class<? extends Payload>[]payload() default {};



	/**
	 * Validator for the {@link Unequal} annotation.
	 */
	static class Validator implements ConstraintValidator<Unequal,Object> {
		private volatile String value;


		/**
		 * {@inheritDoc}
		 */
		public void initialize (final Unequal annotation) {
			this.value = annotation.value();
		}


		/**
		 * {@inheritDoc}
		 */
		public boolean isValid (final Object object, final ConstraintValidatorContext context) {
			return !Objects.equals(this.value, object.toString());
		}
	}
}