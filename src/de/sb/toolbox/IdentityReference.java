package de.sb.toolbox;

/**
 * Reusable class referencing a single object, or {@code null}, checking for object identity when evaluating reference equality.
 * Any collection containing identity references as elements is effectively an identity collection of the referenced objects.
 * @param <T> the object type
 */
@Copyright(year = 2015, holders = "Sascha Baumeister")
public class IdentityReference<T> {
	private volatile T object = null;


	/**
	 * Returns the referenced object.
	 * @return the object
	 */
	public T get () {
		return this.object;
	}


	/**
	 * Sets the referenced object.
	 * @param object the object
	 */
	public void put (final T object) {
		this.object = object;
	}


	/**
	 * {@inheritDoc}
	 * @return the same hash code for the referenced object as would be returned by the default method {@link Object#hashCode()}
	 *         , whether it is overridden or not.
	 */
	@Override
	public int hashCode () {
		return System.identityHashCode(this.object);
	}


	/**
	 * {@inheritDoc}
	 * @return whether or not the given object is a reference, and both it and the receiver reference the identical object.
	 */
	@Override
	public boolean equals (final Object object) {
		if (object == null || !(object instanceof IdentityReference)) return false;

		final IdentityReference<?> reference = (IdentityReference<?>) object;
		return this.object == reference.object;
	}
}