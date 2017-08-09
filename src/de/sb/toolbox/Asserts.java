package de.sb.toolbox;

import java.util.Objects;


/**
 * Additional helper methods for JUnit test cases.
 */
@Copyright(year = 2015, holders = "Sascha Baumeister")
public class Asserts {

	/**
	 * Interface describing code blocks, which are related to {@link java.lang.Runnable} and
	 * {@link java.util.concurrent.Callable}. Similarly to {@code Callable}, but in opposition to {@code Runnable}, code block
	 * excution may declare any exception type. Similarly to {@code Runnable}, but in opposition to {@code Callable}, code block
	 * execution cannot declare a result.
	 */
	@FunctionalInterface
	static public interface Block {

		/**
		 * Executes this code block.
		 * @throws Throwable if there is a problem.
		 */
		void execute () throws Throwable;
	}


	/**
	 * Simplifies checks for the occurrence of exceptions in code blocks.
	 * @param block the code block, usually defined using a lambda expression
	 * @param expectedExceptionType the expected exception type
	 */
	public static void assertThrows (final Block block, final Class<? extends Throwable> expectedExceptionType) {
		Objects.requireNonNull(expectedExceptionType);
		try {
			block.execute();
		} catch (final Throwable exception) {
			if (expectedExceptionType.isInstance(exception)) return;
			final String message = String.format("Unexpected exception, expected <%s> but was <%s>", expectedExceptionType.getName(), exception.getClass().getName());
			throw new AssertionError(message, exception);
		}
		throw new AssertionError(String.format("Expected exception: <%s>", expectedExceptionType.getName()));
	}


	/**
	 * Simplifies checks for the non-occurrence of exceptions in code blocks.
	 * @param block the code block, usually defined using a lambda expression
	 * @param expectedExceptionType the expected exception type
	 */
	public static void assertThrowsNot (final Block block, final Class<? extends Throwable> unexpectedExceptionType) {
		Objects.requireNonNull(unexpectedExceptionType);
		try {
			block.execute();
		} catch (final Throwable exception) {
			if (unexpectedExceptionType.isInstance(exception)) {
				throw new AssertionError(String.format("Unexpected exception: <%s>", exception.getClass().getName()), exception);
			}
		}
	}
}